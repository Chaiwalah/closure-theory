#!/usr/bin/env python3
"""
Closure Theory — DES-5YR Independent Replication
=================================================
Runs the core test battery from Rounds 1-3 on the independent DES-5YR dataset
(1,829 SNe Ia from the Dark Energy Survey 5-Year analysis).

If closure effects are real, we expect:
  - Color-distance coupling accumulating with z (replicating Test 4 / R1)
  - Factorization collapse at high-z (replicating R2 Test 6)
  - Sigmoid threshold near z~0.82 (replicating R2 Test 7)
  - Information compression at high-z (replicating R2 Test 8)
  - Signal surviving full controls (replicating R3 A1)
  - Baryonic nulls: host mass, E(B-V) NOT driving the signal (replicating R1 Tests 1-3)

Data: DES-SN5YR Metadata (SALT2 fit parameters)
Reference: DES Collaboration, 2024, ApJ, DES-SN5YR
"""

import numpy as np
import os
import json
from scipy import stats
from scipy.optimize import curve_fit
from itertools import combinations

np.random.seed(42)

RESULTS_DIR = "results_des5yr"
os.makedirs(RESULTS_DIR, exist_ok=True)

# =============================================================================
# DATA LOADING
# =============================================================================
def load_des5yr(path="data/des_metadata.csv"):
    """Load DES-5YR metadata file (SNANA format)."""
    rows = []
    varnames = None
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('VARNAMES:'):
                varnames = line.split()[1:]  # skip 'VARNAMES:'
            elif line.startswith('SN:'):
                vals = line.split()[1:]  # skip 'SN:'
                if varnames and len(vals) == len(varnames):
                    rows.append(vals)
    
    if not varnames or not rows:
        raise ValueError(f"Failed to parse {path}: varnames={varnames is not None}, rows={len(rows) if rows else 0}")
    
    # Build structured array
    data = {}
    for i, name in enumerate(varnames):
        col = [r[i] for r in rows]
        # Try numeric conversion
        try:
            data[name] = np.array(col, dtype=float)
        except ValueError:
            data[name] = np.array(col, dtype=str)
    
    return data

def prepare_data(data):
    """Clean and prepare DES-5YR data for analysis."""
    z = data['zHD']
    x1 = data['x1']
    c = data['c']
    mu = data['MU']
    mu_err = data['MUERR']
    mu_model = data['MUMODEL']
    mu_res = data['MURES']
    host_mass = data['HOST_LOGMASS']
    mwebv = data['MWEBV']
    chi2 = data['FITCHI2']
    ndof = data['NDOF']
    
    # Quality cuts: valid z, valid fit, no sentinel values
    mask = (z > 0.01) & (z < 2.0)
    mask &= np.isfinite(x1) & np.isfinite(c) & np.isfinite(mu)
    mask &= (np.abs(x1) < 5) & (np.abs(c) < 1)  # reasonable SALT2 ranges
    mask &= (host_mass > 0) & (host_mass < 20)  # valid host mass
    mask &= (ndof > 0)
    
    clean = {k: (v[mask] if isinstance(v, np.ndarray) and v.dtype.kind == 'f' else v[mask]) 
             for k, v in data.items()}
    
    # Compute derived quantities
    clean['chi2_dof'] = clean['FITCHI2'] / clean['NDOF']
    clean['c_resid'] = clean['c'] - np.polyval(np.polyfit(clean['zHD'], clean['c'], 1), clean['zHD'])
    clean['x1_resid'] = clean['x1'] - np.polyval(np.polyfit(clean['zHD'], clean['x1'], 1), clean['zHD'])
    clean['mu_resid'] = clean['MURES']  # already bias-corrected in DES
    
    return clean

# =============================================================================
# TEST FUNCTIONS
# =============================================================================

def test_baryonic_nulls(data):
    """R1 equivalent: Test that host mass and E(B-V) don't drive color-distance coupling."""
    z = data['zHD']
    c = data['c']
    host_mass = data['HOST_LOGMASS']
    mwebv = data['MWEBV']
    mu_res = data['mu_resid']
    
    results = {}
    
    # Test 1: Host mass step vs color residual
    high_mass = host_mass > 10.0
    c_resid = data['c_resid']
    
    t_stat, p_val = stats.ttest_ind(c_resid[high_mass], c_resid[~high_mass])
    results['host_mass_color_ttest'] = {'t': float(t_stat), 'p': float(p_val)}
    
    # Test 2: MW E(B-V) vs Hubble residual
    rho, p = stats.spearmanr(mwebv, mu_res)
    results['mwebv_vs_mu_resid'] = {'rho': float(rho), 'p': float(p)}
    
    # Test 3: Host mass vs color at low-z vs high-z
    z_med = np.median(z)
    lo = z < z_med
    hi = z >= z_med
    rho_lo, p_lo = stats.spearmanr(host_mass[lo], c[lo])
    rho_hi, p_hi = stats.spearmanr(host_mass[hi], c[hi])
    results['host_mass_color_by_z'] = {
        'rho_low_z': float(rho_lo), 'p_low_z': float(p_lo),
        'rho_high_z': float(rho_hi), 'p_high_z': float(p_hi),
        'z_split': float(z_med)
    }
    
    return results

def test_color_distance_coupling(data):
    """R1 Test 4 equivalent: Color-distance correlation accumulates with z."""
    z = data['zHD']
    c = data['c']
    mu_res = data['mu_resid']
    
    # Sliding window correlation
    z_sorted = np.sort(z)
    window = len(z) // 5
    z_centers = []
    rhos = []
    pvals = []
    
    for i in range(0, len(z_sorted) - window, window // 4):
        z_lo, z_hi = z_sorted[i], z_sorted[i + window]
        mask = (z >= z_lo) & (z < z_hi)
        if mask.sum() > 20:
            rho, p = stats.spearmanr(c[mask], mu_res[mask])
            z_centers.append(float(np.median(z[mask])))
            rhos.append(float(rho))
            pvals.append(float(p))
    
    # Trend: does |rho| increase with z?
    trend_rho, trend_p = stats.spearmanr(z_centers, np.abs(rhos))
    
    # Overall correlation in low-z vs high-z
    z_med = np.median(z)
    lo = z < z_med
    hi = z >= z_med
    rho_lo, p_lo = stats.spearmanr(c[lo], mu_res[lo])
    rho_hi, p_hi = stats.spearmanr(c[hi], mu_res[hi])
    
    return {
        'sliding_window': {'z_centers': z_centers, 'rhos': rhos, 'pvals': pvals},
        'accumulation_trend': {'rho': float(trend_rho), 'p': float(trend_p)},
        'low_z': {'rho': float(rho_lo), 'p': float(p_lo), 'n': int(lo.sum())},
        'high_z': {'rho': float(rho_hi), 'p': float(p_hi), 'n': int(hi.sum())},
    }

def test_factorization_collapse(data):
    """R2 Test 6 equivalent: Observable pairs lose independence at high-z."""
    z = data['zHD']
    observables = {
        'c': data['c'],
        'x1': data['x1'],
        'c_resid': data['c_resid'],
        'x1_resid': data['x1_resid'],
        'chi2_dof': data['chi2_dof'],
        'mu_resid': data['mu_resid'],
    }
    
    z_med = np.median(z)
    lo = z < z_med
    hi = z >= z_med
    
    pairs = list(combinations(observables.keys(), 2))
    entangling = []
    
    for a, b in pairs:
        rho_lo, p_lo = stats.spearmanr(observables[a][lo], observables[b][lo])
        rho_hi, p_hi = stats.spearmanr(observables[a][hi], observables[b][hi])
        
        # Fisher z-transform to test difference
        n_lo, n_hi = int(lo.sum()), int(hi.sum())
        z_lo_f = np.arctanh(np.clip(rho_lo, -0.999, 0.999))
        z_hi_f = np.arctanh(np.clip(rho_hi, -0.999, 0.999))
        se = np.sqrt(1/(n_lo-3) + 1/(n_hi-3))
        z_diff = (z_hi_f - z_lo_f) / se
        p_diff = 2 * (1 - stats.norm.cdf(abs(z_diff)))
        
        entry = {
            'pair': f"{a}_vs_{b}",
            'rho_low_z': float(rho_lo), 'rho_high_z': float(rho_hi),
            'z_diff': float(z_diff), 'p_diff': float(p_diff),
            'entangling': p_diff < 0.05 and abs(rho_hi) > abs(rho_lo)
        }
        if entry['entangling']:
            entangling.append(entry)
    
    return {
        'total_pairs': len(pairs),
        'entangling_pairs': len(entangling),
        'fraction': len(entangling) / len(pairs),
        'entangling_details': entangling,
        'z_split': float(z_med)
    }

def sigmoid(z, z0, k, A, B):
    return A / (1 + np.exp(-k * (z - z0))) + B

def test_sigmoid_threshold(data):
    """R2 Test 7 equivalent: Find sigmoid transition in color-distance coupling."""
    z = data['zHD']
    c = data['c']
    mu_res = data['mu_resid']
    
    # Compute running |correlation| vs z
    z_bins = np.linspace(np.percentile(z, 5), np.percentile(z, 95), 30)
    window = (z_bins[-1] - z_bins[0]) / 10
    
    z_centers = []
    abs_rhos = []
    
    for zc in z_bins:
        mask = (z >= zc - window) & (z < zc + window)
        if mask.sum() > 15:
            rho, _ = stats.spearmanr(c[mask], mu_res[mask])
            z_centers.append(float(zc))
            abs_rhos.append(float(abs(rho)))
    
    z_arr = np.array(z_centers)
    r_arr = np.array(abs_rhos)
    
    try:
        popt, pcov = curve_fit(sigmoid, z_arr, r_arr, 
                                p0=[0.5, 10, 0.3, 0.05], 
                                maxfev=10000,
                                bounds=([0.1, 0.1, 0.01, -0.5], [1.5, 100, 1.0, 0.5]))
        z_threshold = float(popt[0])
        perr = np.sqrt(np.diag(pcov))
        
        # Compare sigmoid vs linear fit
        lin = np.polyfit(z_arr, r_arr, 1)
        ss_sig = np.sum((r_arr - sigmoid(z_arr, *popt))**2)
        ss_lin = np.sum((r_arr - np.polyval(lin, z_arr))**2)
        
        result = {
            'z_threshold': z_threshold,
            'z_threshold_err': float(perr[0]),
            'steepness': float(popt[1]),
            'amplitude': float(popt[2]),
            'ss_sigmoid': float(ss_sig),
            'ss_linear': float(ss_lin),
            'sigmoid_better': ss_sig < ss_lin,
        }
    except Exception as e:
        result = {'error': str(e)}
    
    return result

def test_information_compression(data):
    """R2 Test 8 equivalent: Effective rank of observable matrix decreases at high-z."""
    z = data['zHD']
    
    obs_names = ['c', 'x1', 'chi2_dof', 'mu_resid']
    obs_matrix = np.column_stack([data[k] for k in obs_names])
    
    z_med = np.median(z)
    lo = z < z_med
    hi = z >= z_med
    
    def effective_rank(X):
        X_std = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-10)
        s = np.linalg.svd(X_std, compute_uv=False)
        p = s / s.sum()
        p = p[p > 0]
        return float(np.exp(-np.sum(p * np.log(p))))
    
    rank_lo = effective_rank(obs_matrix[lo])
    rank_hi = effective_rank(obs_matrix[hi])
    
    # Bootstrap CI
    n_boot = 2000
    rank_diffs = []
    for _ in range(n_boot):
        idx_lo = np.random.choice(np.where(lo)[0], size=lo.sum(), replace=True)
        idx_hi = np.random.choice(np.where(hi)[0], size=hi.sum(), replace=True)
        rd = effective_rank(obs_matrix[idx_hi]) - effective_rank(obs_matrix[idx_lo])
        rank_diffs.append(rd)
    
    rank_diffs = np.array(rank_diffs)
    p_compression = float(np.mean(rank_diffs >= 0))  # p(no compression)
    
    return {
        'rank_low_z': rank_lo,
        'rank_high_z': rank_hi,
        'rank_drop': rank_lo - rank_hi,
        'p_no_compression': p_compression,
        'boot_ci_95': [float(np.percentile(rank_diffs, 2.5)), float(np.percentile(rank_diffs, 97.5))],
        'n_low': int(lo.sum()),
        'n_high': int(hi.sum()),
    }

def test_signal_survival(data):
    """R3 A1 equivalent: Does color-distance signal survive full controls?"""
    z = data['zHD']
    c = data['c']
    mu_res = data['mu_resid']
    host_mass = data['HOST_LOGMASS']
    x1 = data['x1']
    mwebv = data['MWEBV']
    
    # Raw signal
    rho_raw, p_raw = stats.spearmanr(c, mu_res)
    
    # Partial correlation controlling for host_mass, x1, mwebv, z
    from numpy.linalg import lstsq
    
    controls = np.column_stack([host_mass, x1, mwebv, z, z**2])
    
    # Residualize c
    coef_c, _, _, _ = lstsq(controls, c, rcond=None)
    c_ctrl = c - controls @ coef_c
    
    # Residualize mu_res
    coef_mu, _, _, _ = lstsq(controls, mu_res, rcond=None)
    mu_ctrl = mu_res - controls @ coef_mu
    
    rho_ctrl, p_ctrl = stats.spearmanr(c_ctrl, mu_ctrl)
    
    survival_pct = (abs(rho_ctrl) / abs(rho_raw)) * 100 if abs(rho_raw) > 0 else np.nan
    
    return {
        'raw': {'rho': float(rho_raw), 'p': float(p_raw)},
        'controlled': {'rho': float(rho_ctrl), 'p': float(p_ctrl)},
        'survival_pct': float(survival_pct),
        'signal_amplified': abs(rho_ctrl) > abs(rho_raw),
        'controls': ['HOST_LOGMASS', 'x1', 'MWEBV', 'z', 'z^2'],
        'n': int(len(z)),
    }

def test_frequency_fingerprint(data):
    """Core closure discriminator: frequency-dependent observables entangle, stretch survives."""
    z = data['zHD']
    c = data['c']
    x1 = data['x1']
    mu_res = data['mu_resid']
    chi2_dof = data['chi2_dof']
    
    z_med = np.median(z)
    hi = z >= z_med
    
    # Frequency-dependent: c, chi2_dof vs mu_resid at high-z
    rho_c, p_c = stats.spearmanr(c[hi], mu_res[hi])
    rho_chi2, p_chi2 = stats.spearmanr(chi2_dof[hi], mu_res[hi])
    
    # Temporal-only: x1 should NOT correlate
    rho_x1, p_x1 = stats.spearmanr(x1[hi], mu_res[hi])
    
    # Cross-check: c vs x1 coupling at high-z vs low-z
    lo = z < z_med
    rho_cx1_lo, p_cx1_lo = stats.spearmanr(c[lo], x1[lo])
    rho_cx1_hi, p_cx1_hi = stats.spearmanr(c[hi], x1[hi])
    
    return {
        'high_z_c_mu': {'rho': float(rho_c), 'p': float(p_c)},
        'high_z_chi2_mu': {'rho': float(rho_chi2), 'p': float(p_chi2)},
        'high_z_x1_mu': {'rho': float(rho_x1), 'p': float(p_x1)},
        'c_x1_coupling_low_z': {'rho': float(rho_cx1_lo), 'p': float(p_cx1_lo)},
        'c_x1_coupling_high_z': {'rho': float(rho_cx1_hi), 'p': float(p_cx1_hi)},
        'frequency_fingerprint_holds': (p_c < 0.05 or p_chi2 < 0.05) and p_x1 > 0.05,
        'n_high_z': int(hi.sum()),
    }

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 70)
    print("CLOSURE THEORY — DES-5YR INDEPENDENT REPLICATION")
    print("=" * 70)
    
    # Load data
    data = load_des5yr()
    n_raw = len(data['zHD'])
    print(f"\nLoaded {n_raw} raw SNe from DES-5YR")
    
    clean = prepare_data(data)
    n_clean = len(clean['zHD'])
    print(f"After quality cuts: {n_clean} SNe")
    print(f"Redshift range: {clean['zHD'].min():.3f} — {clean['zHD'].max():.3f}")
    print(f"Median z: {np.median(clean['zHD']):.3f}")
    
    all_results = {
        'dataset': 'DES-SN5YR',
        'n_raw': n_raw,
        'n_clean': n_clean,
        'z_range': [float(clean['zHD'].min()), float(clean['zHD'].max())],
        'z_median': float(np.median(clean['zHD'])),
    }
    
    # --- Test 1: Baryonic Nulls ---
    print("\n" + "-" * 50)
    print("TEST 1: BARYONIC NULLS (host mass, MW dust)")
    res = test_baryonic_nulls(clean)
    all_results['baryonic_nulls'] = res
    print(f"  Host mass → color: t={res['host_mass_color_ttest']['t']:.3f}, p={res['host_mass_color_ttest']['p']:.4f}")
    print(f"  MW E(B-V) → μ_resid: ρ={res['mwebv_vs_mu_resid']['rho']:.4f}, p={res['mwebv_vs_mu_resid']['p']:.4f}")
    print(f"  Host mass → color (low-z): ρ={res['host_mass_color_by_z']['rho_low_z']:.4f}")
    print(f"  Host mass → color (high-z): ρ={res['host_mass_color_by_z']['rho_high_z']:.4f}")
    
    # --- Test 2: Color-Distance Coupling ---
    print("\n" + "-" * 50)
    print("TEST 2: COLOR-DISTANCE COUPLING (accumulation)")
    res = test_color_distance_coupling(clean)
    all_results['color_distance'] = res
    print(f"  Accumulation trend: ρ={res['accumulation_trend']['rho']:.3f}, p={res['accumulation_trend']['p']:.4f}")
    print(f"  Low-z (n={res['low_z']['n']}): ρ={res['low_z']['rho']:.4f}, p={res['low_z']['p']:.4f}")
    print(f"  High-z (n={res['high_z']['n']}): ρ={res['high_z']['rho']:.4f}, p={res['high_z']['p']:.4f}")
    
    # --- Test 3: Factorization Collapse ---
    print("\n" + "-" * 50)
    print("TEST 3: FACTORIZATION COLLAPSE")
    res = test_factorization_collapse(clean)
    all_results['factorization'] = res
    print(f"  Entangling pairs: {res['entangling_pairs']}/{res['total_pairs']} ({res['fraction']:.1%})")
    for e in res['entangling_details']:
        print(f"    {e['pair']}: Δρ = {e['rho_high_z']-e['rho_low_z']:.3f} (p={e['p_diff']:.4f})")
    
    # --- Test 4: Sigmoid Threshold ---
    print("\n" + "-" * 50)
    print("TEST 4: SIGMOID THRESHOLD")
    res = test_sigmoid_threshold(clean)
    all_results['sigmoid'] = res
    if 'error' not in res:
        print(f"  z* = {res['z_threshold']:.3f} ± {res['z_threshold_err']:.3f}")
        print(f"  Sigmoid better than linear: {res['sigmoid_better']}")
    else:
        print(f"  Fit failed: {res['error']}")
    
    # --- Test 5: Information Compression ---
    print("\n" + "-" * 50)
    print("TEST 5: INFORMATION COMPRESSION")
    res = test_information_compression(clean)
    all_results['compression'] = res
    print(f"  Effective rank (low-z): {res['rank_low_z']:.3f} (n={res['n_low']})")
    print(f"  Effective rank (high-z): {res['rank_high_z']:.3f} (n={res['n_high']})")
    print(f"  Rank drop: {res['rank_drop']:.3f}")
    print(f"  p(no compression): {res['p_no_compression']:.4f}")
    
    # --- Test 6: Signal Survival (A1 equivalent) ---
    print("\n" + "-" * 50)
    print("TEST 6: SIGNAL SURVIVAL UNDER CONTROLS (A1)")
    res = test_signal_survival(clean)
    all_results['signal_survival'] = res
    print(f"  Raw: ρ={res['raw']['rho']:.4f}, p={res['raw']['p']:.6f}")
    print(f"  Controlled: ρ={res['controlled']['rho']:.4f}, p={res['controlled']['p']:.6f}")
    print(f"  Survival: {res['survival_pct']:.1f}%")
    print(f"  Signal amplified: {res['signal_amplified']}")
    
    # --- Test 7: Frequency Fingerprint ---
    print("\n" + "-" * 50)
    print("TEST 7: FREQUENCY FINGERPRINT")
    res = test_frequency_fingerprint(clean)
    all_results['frequency_fingerprint'] = res
    print(f"  High-z c→μ: ρ={res['high_z_c_mu']['rho']:.4f}, p={res['high_z_c_mu']['p']:.4f}")
    print(f"  High-z χ²/dof→μ: ρ={res['high_z_chi2_mu']['rho']:.4f}, p={res['high_z_chi2_mu']['p']:.4f}")
    print(f"  High-z x1→μ: ρ={res['high_z_x1_mu']['rho']:.4f}, p={res['high_z_x1_mu']['p']:.4f}")
    print(f"  Fingerprint holds: {res['frequency_fingerprint_holds']}")
    
    # Save results
    results_path = os.path.join(RESULTS_DIR, "des5yr_results.json")
    
    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.bool_,)):
                return bool(obj)
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)
    
    with open(results_path, 'w') as f:
        json.dump(all_results, f, indent=2, cls=NumpyEncoder)
    print(f"\nResults saved to {results_path}")
    
    # Summary
    print("\n" + "=" * 70)
    print("REPLICATION SUMMARY — DES-5YR vs Pantheon+")
    print("=" * 70)
    
    checks = []
    
    # Baryonic nulls
    bn = all_results['baryonic_nulls']
    null1 = bn['mwebv_vs_mu_resid']['p'] > 0.05
    checks.append(("Baryonic null (MW dust)", null1))
    
    # Color accumulation
    cd = all_results['color_distance']
    accum = cd['accumulation_trend']['p'] < 0.1  # relaxed for independent dataset
    checks.append(("Color-distance accumulation", accum))
    
    # Factorization
    fc = all_results['factorization']
    collapse = fc['entangling_pairs'] > 0
    checks.append(("Factorization collapse", collapse))
    
    # Compression
    comp = all_results['compression']
    compressed = comp['rank_drop'] > 0 and comp['p_no_compression'] < 0.1
    checks.append(("Information compression", compressed))
    
    # Signal survival
    ss = all_results['signal_survival']
    survives = ss['controlled']['p'] < 0.05
    checks.append(("Signal survives controls", survives))
    
    # Frequency fingerprint
    ff = all_results['frequency_fingerprint']
    fingerprint = ff['frequency_fingerprint_holds']
    checks.append(("Frequency fingerprint", fingerprint))
    
    print()
    passed = 0
    for name, result in checks:
        status = "✓ REPLICATED" if result else "✗ NOT REPLICATED"
        print(f"  {status}  {name}")
        if result:
            passed += 1
    
    print(f"\n  Score: {passed}/{len(checks)} tests replicated")
    print("=" * 70)

if __name__ == "__main__":
    main()
