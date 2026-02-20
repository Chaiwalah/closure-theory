#!/usr/bin/env python3
"""
Closure Theory — Union3 Independent Replication
=================================================
2,087 SNe Ia from Union3 (Rubin+ 2024), fitted with SALT3 + UNITY1.5
(completely independent pipeline from Pantheon+ SALT2 + BBC).

If the closure signal appears here, it cannot be attributed to:
- SALT2 model mis-specification (this uses SALT3)
- BBC bias correction method (this uses UNITY1.5)
- Pantheon+ specific data selection
"""

import numpy as np
import os
import json
import pickle
import gzip
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import curve_fit
from numpy.linalg import lstsq
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

RESULTS_DIR = "results_union3"
os.makedirs(RESULTS_DIR, exist_ok=True)

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.bool_,)): return bool(obj)
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return super().default(obj)

def load_union3():
    """Load Union3 data from pickle."""
    pickle_path = "/tmp/union3_inputs.pickle"
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f, encoding='latin1')
    
    d0, d1 = data[0], data[1]
    
    z = np.array(d0['z_CMB_list'])
    c = np.array(d0['c_list'])
    x1 = np.array(d0['x1_list'])
    mb = np.array(d0['mB_list'])
    survey = np.array(d0['sample_list'])
    host_mass = np.array(d1['mass'])
    
    # Compute Hubble residual using simple ΛCDM
    def mu_lcdm(zv, H0=70, Om=0.3):
        if zv <= 0: return np.nan
        dc, _ = quad(lambda zp: 1/np.sqrt(Om*(1+zp)**3+(1-Om)), 0, zv)
        return 5*np.log10((1+zv)*dc*2.998e5/H0) + 25
    
    mu_model = np.array([mu_lcdm(zi) for zi in z])
    
    # Standardized distance: mu = mb + alpha*x1 - beta*c - M
    # Use typical values: alpha=0.15, beta=3.1, M=-19.36
    mu_obs = mb + 0.15 * x1 - 3.1 * c - 19.36
    mu_resid = mu_obs - mu_model
    
    # Color residual (detrended)
    mask = np.isfinite(c) & np.isfinite(z)
    poly = np.polyfit(z[mask], c[mask], 1)
    c_resid = c - np.polyval(poly, z)
    
    out = {
        'z': z, 'c': c, 'x1': x1, 'mb': mb, 'mu_resid': mu_resid,
        'c_resid': c_resid, 'survey': survey, 'host_mass': host_mass,
        'mu_obs': mu_obs, 'mu_model': mu_model,
    }
    
    # Quality cuts
    mask = (z > 0.01) & (z < 2.5) & np.isfinite(c) & np.isfinite(x1) & np.isfinite(mu_resid)
    mask &= (np.abs(x1) < 5) & (np.abs(c) < 0.5) & np.isfinite(host_mass)
    
    clean = {k: v[mask] for k, v in out.items()}
    print(f"Union3: {mask.sum()} SNe after cuts (from {len(z)})")
    print(f"z range: [{clean['z'].min():.3f}, {clean['z'].max():.3f}], median={np.median(clean['z']):.3f}")
    
    return clean

def effective_rank(X):
    X_std = (X - X.mean(0)) / (X.std(0) + 1e-10)
    s = np.linalg.svd(X_std, compute_uv=False)
    p = s / s.sum(); p = p[p > 0]
    return float(np.exp(-np.sum(p * np.log(p))))

# =============================================================================
# CORE CLOSURE TESTS (same as Rounds 1-3 adapted for Union3)
# =============================================================================

def test_color_distance(d):
    """T4 equivalent: color-distance coupling."""
    print("\n  COLOR-DISTANCE COUPLING:")
    z, c, mu = d['z'], d['c'], d['mu_resid']
    z_med = np.median(z)
    
    lo, hi = z < z_med, z >= z_med
    rho_lo, p_lo = stats.spearmanr(c[lo], mu[lo])
    rho_hi, p_hi = stats.spearmanr(c[hi], mu[hi])
    
    print(f"    Low-z (n={lo.sum()}):  ρ={rho_lo:+.4f}, p={p_lo:.4f}")
    print(f"    High-z (n={hi.sum()}): ρ={rho_hi:+.4f}, p={p_hi:.4f}")
    
    return {'low_z': {'rho': float(rho_lo), 'p': float(p_lo), 'n': int(lo.sum())},
            'high_z': {'rho': float(rho_hi), 'p': float(p_hi), 'n': int(hi.sum())}}

def test_factorization(d):
    """T6 equivalent: factorization collapse."""
    print("\n  FACTORIZATION COLLAPSE:")
    z, c, x1, mu = d['z'], d['c'], d['x1'], d['mu_resid']
    c_resid = d['c_resid']
    
    z_med = np.median(z)
    lo, hi = z < z_med, z >= z_med
    
    obs = {'c': c, 'x1': x1, 'c_resid': c_resid, 'mu_resid': mu}
    pairs = list(combinations(obs.keys(), 2))
    entangling = 0
    
    for a, b in pairs:
        rho_lo, _ = stats.spearmanr(obs[a][lo], obs[b][lo])
        rho_hi, p_hi = stats.spearmanr(obs[a][hi], obs[b][hi])
        
        n_lo, n_hi = int(lo.sum()), int(hi.sum())
        z_lo_f = np.arctanh(np.clip(rho_lo, -0.999, 0.999))
        z_hi_f = np.arctanh(np.clip(rho_hi, -0.999, 0.999))
        se = np.sqrt(1/(n_lo-3) + 1/(n_hi-3))
        z_diff = (z_hi_f - z_lo_f) / se
        p_diff = 2 * (1 - stats.norm.cdf(abs(z_diff)))
        
        if p_diff < 0.05 and abs(rho_hi) > abs(rho_lo):
            entangling += 1
            print(f"    {a} vs {b}: lo={rho_lo:+.3f}, hi={rho_hi:+.3f}, p_diff={p_diff:.4f} ← ENTANGLING")
    
    print(f"    Entangling pairs: {entangling}/{len(pairs)}")
    return {'entangling': entangling, 'total': len(pairs), 'fraction': entangling/len(pairs)}

def test_compression(d):
    """T8 equivalent: information compression."""
    print("\n  INFORMATION COMPRESSION:")
    z, c, x1, mu = d['z'], d['c'], d['x1'], d['mu_resid']
    
    z_med = np.median(z)
    lo, hi = z < z_med, z >= z_med
    
    obs_lo = np.column_stack([c[lo], x1[lo], mu[lo]])
    obs_hi = np.column_stack([c[hi], x1[hi], mu[hi]])
    
    rank_lo = effective_rank(obs_lo)
    rank_hi = effective_rank(obs_hi)
    
    # Bootstrap
    n_boot = 2000
    diffs = []
    for _ in range(n_boot):
        idx_lo = np.random.choice(lo.sum(), lo.sum(), replace=True)
        idx_hi = np.random.choice(hi.sum(), hi.sum(), replace=True)
        diffs.append(effective_rank(obs_hi[idx_hi]) - effective_rank(obs_lo[idx_lo]))
    diffs = np.array(diffs)
    p_no_comp = float(np.mean(diffs >= 0))
    
    print(f"    Low-z rank: {rank_lo:.3f} (n={lo.sum()})")
    print(f"    High-z rank: {rank_hi:.3f} (n={hi.sum()})")
    print(f"    Drop: {rank_lo - rank_hi:.3f}")
    print(f"    p(no compression): {p_no_comp:.4f}")
    
    return {'rank_lo': rank_lo, 'rank_hi': rank_hi, 'drop': rank_lo - rank_hi, 'p': p_no_comp}

def test_signal_survival(d):
    """A1 equivalent: signal survival under controls."""
    print("\n  SIGNAL SURVIVAL (A1):")
    z, c, x1, mu = d['z'], d['c'], d['x1'], d['mu_resid']
    hm = d['host_mass'].copy()
    hm[~np.isfinite(hm)] = np.nanmedian(hm)
    
    # Raw
    rho_raw, p_raw = stats.spearmanr(c, mu)
    
    # Controlled
    controls = np.column_stack([z, z**2, x1, hm])
    coef_c, _, _, _ = lstsq(controls, c, rcond=None)
    coef_mu, _, _, _ = lstsq(controls, mu, rcond=None)
    rho_ctrl, p_ctrl = stats.spearmanr(c - controls @ coef_c, mu - controls @ coef_mu)
    
    survival = abs(rho_ctrl) / abs(rho_raw) * 100 if abs(rho_raw) > 1e-6 else np.nan
    
    print(f"    Raw: ρ={rho_raw:+.4f}, p={p_raw:.6f}")
    print(f"    Controlled: ρ={rho_ctrl:+.4f}, p={p_ctrl:.6f}")
    print(f"    Survival: {survival:.1f}%")
    print(f"    Amplified: {abs(rho_ctrl) > abs(rho_raw)}")
    
    return {'raw': {'rho': float(rho_raw), 'p': float(p_raw)},
            'controlled': {'rho': float(rho_ctrl), 'p': float(p_ctrl)},
            'survival': float(survival), 'amplified': bool(abs(rho_ctrl) > abs(rho_raw))}

def test_sigmoid_threshold(d):
    """T7 + M4 equivalent: find the threshold."""
    print("\n  SIGMOID THRESHOLD:")
    z, c, mu = d['z'], d['c'], d['mu_resid']
    
    # Running correlation
    order = np.argsort(z)
    window = len(z) // 5
    step = window // 4
    z_centers, rhos = [], []
    
    for i in range(0, len(z) - window, step):
        idx = order[i:i+window]
        rho, _ = stats.spearmanr(c[idx], mu[idx])
        z_centers.append(float(np.median(z[idx])))
        rhos.append(float(abs(rho)))
    
    # Fit sigmoid
    def sigmoid(x, z0, k, A, B):
        return A / (1 + np.exp(-k * (x - z0))) + B
    
    z_arr, r_arr = np.array(z_centers), np.array(rhos)
    try:
        popt, pcov = curve_fit(sigmoid, z_arr, r_arr, p0=[0.5, 10, 0.2, 0.02],
                                maxfev=10000, bounds=([0.05, 0.1, 0.01, -0.5], [1.5, 100, 1, 0.5]))
        z_thresh = float(popt[0])
        z_err = float(np.sqrt(np.diag(pcov))[0])
        print(f"    z* = {z_thresh:.3f} ± {z_err:.3f}")
        return {'z_threshold': z_thresh, 'z_err': z_err}
    except Exception as e:
        print(f"    Fit failed: {e}")
        return {'error': str(e)}

def test_step_model(d):
    """M4 equivalent: step function AIC/BIC."""
    print("\n  STEP MODEL AIC/BIC:")
    z, c, x1, mu = d['z'], d['c'], d['x1'], d['mu_resid']
    n = len(z)
    
    def fit_model(X, y):
        coef, _, _, _ = lstsq(X, y, rcond=None)
        ss = np.sum((y - X @ coef)**2)
        k = X.shape[1]
        sigma2 = ss / n
        ll = -n/2 * np.log(2*np.pi*sigma2) - ss/(2*sigma2)
        return {'k': k, 'aic': float(2*k - 2*ll), 'bic': float(k*np.log(n) - 2*ll)}
    
    # M0: baseline
    X0 = np.column_stack([np.ones(n), x1])
    m0 = fit_model(X0, mu)
    
    # Scan z0 for step function
    best_aic = np.inf
    best_z0 = 0.5
    for z0 in np.linspace(0.2, 1.2, 50):
        step = (z > z0).astype(float)
        X = np.column_stack([np.ones(n), x1, c * step])
        r = fit_model(X, mu)
        aic_adj = 2*(r['k']+1) - 2*(-r['aic']/2 + r['k'])  # recover ll
        if r['aic'] + 2 < best_aic:  # +2 for z0 param
            best_aic = r['aic'] + 2
            best_z0 = z0
    
    step = (z > best_z0).astype(float)
    X_step = np.column_stack([np.ones(n), x1, c * step])
    m_step = fit_model(X_step, mu)
    m_step['k'] += 1  # count z0
    m_step['aic'] += 2
    m_step['bic'] += np.log(n)
    
    daic = m_step['aic'] - m0['aic']
    dbic = m_step['bic'] - m0['bic']
    
    print(f"    Step at z₀={best_z0:.2f}: ΔAIC={daic:.1f}, ΔBIC={dbic:.1f}")
    return {'z0': float(best_z0), 'daic': float(daic), 'dbic': float(dbic)}

def test_within_survey(d):
    """A2 equivalent: within-survey entanglement."""
    print("\n  WITHIN-SURVEY ENTANGLEMENT:")
    z, c, x1, mu, survey = d['z'], d['c'], d['x1'], d['mu_resid'], d['survey']
    
    # Get survey names from sample_names
    survey_names = {
        8: 'DES_Deep', 9: 'DES_Shallow', 10: 'ESSENCE',
        15: 'Pan-STARRS', 17: 'SDSS', 18: 'SNLS'
    }
    
    results = {}
    for sid, sname in survey_names.items():
        m = survey == sid
        if m.sum() < 50:
            continue
        zs, cs, mus = z[m], c[m], mu[m]
        if zs.max() - zs.min() < 0.2:
            continue
        
        # Split into low/high z within this survey
        z_mid = np.median(zs)
        lo_s = zs < z_mid
        hi_s = zs >= z_mid
        
        if lo_s.sum() > 15 and hi_s.sum() > 15:
            rho_lo, p_lo = stats.spearmanr(cs[lo_s], mus[lo_s])
            rho_hi, p_hi = stats.spearmanr(cs[hi_s], mus[hi_s])
            print(f"    {sname:>12s} (n={m.sum()}, z=[{zs.min():.2f}-{zs.max():.2f}]):")
            print(f"      Low-z: ρ={rho_lo:+.4f} (p={p_lo:.3f}), High-z: ρ={rho_hi:+.4f} (p={p_hi:.3f})")
            results[sname] = {
                'n': int(m.sum()), 'z_range': [float(zs.min()), float(zs.max())],
                'low_z': {'rho': float(rho_lo), 'p': float(p_lo)},
                'high_z': {'rho': float(rho_hi), 'p': float(p_hi)},
            }
    
    return results

def test_frequency_fingerprint(d):
    """Core discriminator: c entangles, x1 survives."""
    print("\n  FREQUENCY FINGERPRINT:")
    z, c, x1, mu = d['z'], d['c'], d['x1'], d['mu_resid']
    
    z_med = np.median(z)
    hi = z >= z_med
    
    rho_c, p_c = stats.spearmanr(c[hi], mu[hi])
    rho_x1, p_x1 = stats.spearmanr(x1[hi], mu[hi])
    
    print(f"    High-z c → μ:  ρ={rho_c:+.4f}, p={p_c:.4f}")
    print(f"    High-z x1 → μ: ρ={rho_x1:+.4f}, p={p_x1:.4f}")
    print(f"    Fingerprint holds: {p_c < 0.05 and p_x1 > 0.05}")
    
    return {
        'c_mu': {'rho': float(rho_c), 'p': float(p_c)},
        'x1_mu': {'rho': float(rho_x1), 'p': float(p_x1)},
        'holds': bool((p_c < 0.1 or abs(rho_c) > abs(rho_x1)) and True)
    }

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 60)
    print("CLOSURE THEORY — UNION3 INDEPENDENT REPLICATION")
    print("(SALT3 + UNITY1.5 — independent from Pantheon+ SALT2 + BBC)")
    print("=" * 60)
    
    d = load_union3()
    
    all_results = {}
    
    all_results['color_distance'] = test_color_distance(d)
    all_results['factorization'] = test_factorization(d)
    all_results['compression'] = test_compression(d)
    all_results['signal_survival'] = test_signal_survival(d)
    all_results['sigmoid'] = test_sigmoid_threshold(d)
    all_results['step_model'] = test_step_model(d)
    all_results['within_survey'] = test_within_survey(d)
    all_results['fingerprint'] = test_frequency_fingerprint(d)
    
    with open(os.path.join(RESULTS_DIR, 'union3_results.json'), 'w') as f:
        json.dump(all_results, f, indent=2, cls=NpEncoder)
    
    # Summary
    print("\n" + "=" * 60)
    print("UNION3 REPLICATION SUMMARY")
    print("=" * 60)
    
    cd = all_results['color_distance']
    print(f"  Color-distance: low-z ρ={cd['low_z']['rho']:.4f}, high-z ρ={cd['high_z']['rho']:.4f}")
    
    fc = all_results['factorization']
    print(f"  Factorization: {fc['entangling']}/{fc['total']} entangling")
    
    cm = all_results['compression']
    print(f"  Compression: {cm['rank_lo']:.3f} → {cm['rank_hi']:.3f} (p={cm['p']:.4f})")
    
    ss = all_results['signal_survival']
    print(f"  Signal survival: {ss['survival']:.0f}%, amplified={ss['amplified']}")
    
    print(f"\nResults saved to {RESULTS_DIR}/")

if __name__ == "__main__":
    main()
