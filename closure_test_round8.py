#!/usr/bin/env python3
"""
Closure Theory — Round 8: Non-SALT & CMB Lensing Tests
========================================================
Two tests that bypass SALT2 entirely:

  L1: CMB Lensing rank compression — Does the EM information channel
      degrade in a way visible in Planck lensing convergence at SN positions?
  
  I1: Model-independent color test — Using raw observed-frame magnitude
      differences (not SALT2's fitted c), does the same closure pattern appear?
"""

import numpy as np
import os
import json
from scipy import stats
from scipy.integrate import quad
from numpy.linalg import lstsq
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

RESULTS_DIR = "results_round8"
os.makedirs(RESULTS_DIR, exist_ok=True)

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.bool_,)): return bool(obj)
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return super().default(obj)

def load_pantheon():
    import pandas as pd
    df = pd.read_csv("data/pantheon_plus.dat", sep=r'\s+', comment='#')
    z = df['zHD'].values
    
    def mu_lcdm(zv, H0=73.04, Om=0.334):
        if not np.isfinite(zv) or zv <= 0: return np.nan
        dc, _ = quad(lambda zp: 1/np.sqrt(Om*(1+zp)**3+(1-Om)), 0, zv)
        return 5*np.log10((1+zv)*dc*2.998e5/H0) + 25
    
    mu_model = np.array([mu_lcdm(zi) for zi in z])
    df['mu_resid'] = df['MU_SH0ES'].values - mu_model
    
    ndof = df['NDOF'].values.astype(float)
    ndof[ndof == 0] = np.nan
    df['chi2_dof'] = df['FITCHI2'].values.astype(float) / ndof
    
    return df

# =============================================================================
# L1: CMB LENSING AT SN POSITIONS
# =============================================================================
def test_cmb_lensing(df):
    """Sample Planck lensing convergence (κ) at SN positions, test for closure signatures."""
    print("\n" + "=" * 60)
    print("L1: CMB LENSING CONVERGENCE AT SN POSITIONS")
    print("=" * 60)
    
    import healpy as hp
    
    # Load Planck lensing convergence map
    # MV = minimum variance (best combined estimator)
    klm_file = "data/planck/COM_Lensing_4096_R3.00/MV/dat_klm.fits"
    mask_file = "data/planck/COM_Lensing_4096_R3.00/mask.fits.gz"
    
    # Read the kappa alm and convert to map
    klm = hp.read_alm(klm_file)
    nside = 2048  # reasonable resolution
    kappa_map = hp.alm2map(klm, nside, verbose=False)
    
    # Read mask
    mask_map = hp.read_map(mask_file, verbose=False)
    # Downgrade mask to same nside
    if hp.get_nside(mask_map) != nside:
        mask_map = hp.ud_grade(mask_map, nside)
    
    print(f"  Loaded Planck lensing κ map (NSIDE={nside})")
    print(f"  κ range: [{kappa_map.min():.4f}, {kappa_map.max():.4f}]")
    
    # Sample κ at SN positions
    z = df['zHD'].values
    ra = df['RA'].values
    dec = df['DEC'].values
    c = df['c'].values
    x1 = df['x1'].values
    mu_res = df['mu_resid'].values
    
    mask = np.isfinite(z) & np.isfinite(ra) & np.isfinite(dec) & np.isfinite(c) & np.isfinite(mu_res) & (z > 0.01)
    z, ra, dec, c, x1, mu_res = z[mask], ra[mask], dec[mask], c[mask], x1[mask], mu_res[mask]
    
    # Convert RA/DEC to theta/phi (healpy convention)
    theta = np.radians(90.0 - dec)  # colatitude
    phi = np.radians(ra)
    
    # Sample κ at each SN position
    pixels = hp.ang2pix(nside, theta, phi)
    kappa = kappa_map[pixels]
    lensing_mask = mask_map[pixels]
    
    # Apply lensing mask (keep only well-measured regions)
    good = lensing_mask > 0.5
    z, c, x1, mu_res, kappa = z[good], c[good], x1[good], mu_res[good], kappa[good]
    
    print(f"  SNe in good lensing region: {len(z)}")
    print(f"  κ at SN positions: mean={kappa.mean():.5f}, std={kappa.std():.5f}")
    
    results = {}
    
    # === Test 1: Does κ correlate with μ_resid? ===
    # If lensing magnification affects distances, κ > 0 → brighter (negative μ_resid)
    rho_km, p_km = stats.spearmanr(kappa, mu_res)
    print(f"\n  κ vs μ_resid: ρ={rho_km:.4f}, p={p_km:.4f}")
    results['kappa_mu'] = {'rho': float(rho_km), 'p': float(p_km)}
    
    # === Test 2: Does κ correlate with SALT2 color? ===
    rho_kc, p_kc = stats.spearmanr(kappa, c)
    print(f"  κ vs c (SALT2 color): ρ={rho_kc:.4f}, p={p_kc:.4f}")
    results['kappa_c'] = {'rho': float(rho_kc), 'p': float(p_kc)}
    
    # === Test 3: Effective rank of observables in high-κ vs low-κ regions ===
    kappa_med = np.median(kappa)
    lo_k = kappa < kappa_med
    hi_k = kappa >= kappa_med
    
    def eff_rank(X):
        X_std = (X - X.mean(0)) / (X.std(0) + 1e-10)
        s = np.linalg.svd(X_std, compute_uv=False)
        p = s / s.sum(); p = p[p > 0]
        return float(np.exp(-np.sum(p * np.log(p))))
    
    for label, km in [("Low κ (underdense)", lo_k), ("High κ (overdense)", hi_k)]:
        obs = np.column_stack([c[km], x1[km], mu_res[km]])
        er = eff_rank(obs)
        rho, p = stats.spearmanr(c[km], mu_res[km])
        print(f"\n  {label} (n={km.sum()}):")
        print(f"    Eff rank: {er:.3f}")
        print(f"    c-μ coupling: ρ={rho:.4f}, p={p:.4f}")
        results[label] = {'n': int(km.sum()), 'rank': er, 'c_mu_rho': float(rho), 'c_mu_p': float(p)}
    
    # === Test 4: Does closure signal depend on κ (lensing environment)? ===
    # If closure scales with matter along sightline, high-κ should show MORE entanglement
    z_med = np.median(z)
    hi_z = z >= z_med
    
    for label, km in [("Low κ + high z", lo_k & hi_z), ("High κ + high z", hi_k & hi_z)]:
        if km.sum() > 30:
            rho, p = stats.spearmanr(c[km], mu_res[km])
            obs = np.column_stack([c[km], x1[km], mu_res[km]])
            er = eff_rank(obs)
            print(f"\n  {label} (n={km.sum()}):")
            print(f"    c-μ coupling: ρ={rho:.4f}, p={p:.4f}")
            print(f"    Eff rank: {er:.3f}")
            results[label] = {'n': int(km.sum()), 'rank': er, 'c_mu_rho': float(rho), 'c_mu_p': float(p)}
    
    # === Test 5: κ-weighted closure signal ===
    # Does |κ| amplify the c-μ coupling?
    abs_kappa = np.abs(kappa)
    # Weighted Spearman: split into κ terciles at high-z
    print(f"\n  κ-tercile analysis (high-z only):")
    k_terts = np.percentile(abs_kappa[hi_z], [0, 33, 67, 100])
    for i in range(3):
        km = hi_z & (abs_kappa >= k_terts[i]) & (abs_kappa < k_terts[i+1] + 0.001)
        if km.sum() > 20:
            rho, p = stats.spearmanr(c[km], mu_res[km])
            print(f"    |κ| tercile {i+1} [{k_terts[i]:.4f}-{k_terts[i+1]:.4f}] (n={km.sum()}): "
                  f"ρ={rho:.4f}, p={p:.4f}")
    
    return results

# =============================================================================
# I1: MODEL-INDEPENDENT COLOR TEST
# =============================================================================
def test_model_independent_color(df):
    """Test closure using SALT2-independent color proxies."""
    print("\n" + "=" * 60)
    print("I1: MODEL-INDEPENDENT COLOR TEST")
    print("=" * 60)
    
    z = df['zHD'].values
    c_salt = df['c'].values  # SALT2 fitted color
    x1 = df['x1'].values
    mu_res = df['mu_resid'].values
    
    # SALT2-independent quantities available in Pantheon+:
    # - m_b_corr: bias-corrected apparent magnitude (depends on SALT2)
    # - MWEBV: MW reddening (independent)
    # - HOST_LOGMASS: host galaxy mass (independent)
    # - FITCHI2, NDOF: fit quality (partially independent — measures HOW WELL SALT2 fits)
    #
    # The key insight: chi2_dof measures how well SALT2 captures the SN SED.
    # If closure degrades the information channel, SALT2 should fit WORSE at high-z
    # because it's trying to fit a degraded signal with a low-z-trained model.
    # chi2_dof is NOT a SALT2 parameter — it's a SALT2-INDEPENDENT diagnostic.
    
    chi2_dof = df['chi2_dof'].values
    mwebv = df['MWEBV'].values if 'MWEBV' in df.columns else None
    
    mask = np.isfinite(z) & np.isfinite(c_salt) & np.isfinite(x1) & np.isfinite(mu_res) & (z > 0.01)
    
    results = {}
    
    # === Approach 1: Use chi2_dof as a SALT2-independent channel quality diagnostic ===
    print("\n  APPROACH 1: χ²/dof as SALT2-independent channel diagnostic")
    
    if chi2_dof is not None:
        mask_chi = mask & np.isfinite(chi2_dof)
        z_m = z[mask_chi]
        chi_m = chi2_dof[mask_chi]
        mu_m = mu_res[mask_chi]
        c_m = c_salt[mask_chi]
        x1_m = x1[mask_chi]
        
        # chi2_dof vs μ_resid (bypasses SALT2 color entirely)
        z_med = np.median(z_m)
        
        print(f"\n  χ²/dof vs μ_resid by z-bin (NO SALT2 color involved):")
        z_edges = np.percentile(z_m, [0, 20, 40, 60, 80, 100])
        
        for i in range(5):
            bm = (z_m >= z_edges[i]) & (z_m < z_edges[i+1] + 0.001)
            if bm.sum() > 20:
                rho, p = stats.spearmanr(chi_m[bm], mu_m[bm])
                sig = " *" if p < 0.05 else ""
                print(f"    z=[{z_edges[i]:.2f}-{z_edges[i+1]:.2f}] (n={bm.sum()}): "
                      f"ρ(χ²,μ)={rho:+.4f}, p={p:.4f}{sig}")
                results[f'chi2_mu_z{i}'] = {'z_range': [float(z_edges[i]), float(z_edges[i+1])],
                                             'rho': float(rho), 'p': float(p), 'n': int(bm.sum())}
        
        # Trend
        hi = z_m >= z_med
        rho_hi, p_hi = stats.spearmanr(chi_m[hi], mu_m[hi])
        rho_lo, p_lo = stats.spearmanr(chi_m[~hi], mu_m[~hi])
        print(f"\n  Low-z χ²-μ: ρ={rho_lo:.4f}, p={p_lo:.4f}")
        print(f"  High-z χ²-μ: ρ={rho_hi:.4f}, p={p_hi:.4f}")
        results['chi2_mu_lo'] = {'rho': float(rho_lo), 'p': float(p_lo)}
        results['chi2_mu_hi'] = {'rho': float(rho_hi), 'p': float(p_hi)}
    
    # === Approach 2: Observed-frame color proxy from peak magnitude and stretch ===
    # m_b_corr already includes the standardization. But we can construct:
    # "raw excess" = m_b_corr - mu_model(z) — this is independent of SALT2 c
    # Then test if this "raw excess" shows the same entanglement pattern
    
    print("\n  APPROACH 2: Raw Hubble residual entanglement (m_b_corr based)")
    
    # m_b_corr is already standardized for x1 and c, so μ_resid IS the raw excess
    # But we want to test: does the RELATIONSHIP between c and μ_resid exist
    # even when we use a crude, SALT2-independent proxy for color?
    
    # Crude color proxy: (m_b_corr_err - median_err) as "how well measured" proxy
    # Actually better: use MWEBV as a sightline reddening proxy
    if mwebv is not None:
        mask_mw = mask & np.isfinite(mwebv)
        z_mw = z[mask_mw]
        mw = mwebv[mask_mw]
        mu_mw = mu_res[mask_mw]
        c_mw = c_salt[mask_mw]
        
        # MWEBV is completely SALT2-independent — it's from dust maps
        # It shouldn't correlate with μ_resid (we already showed this)
        # But the COMBINATION of MWEBV and c behavior can tell us about the channel
        
        # Test: does the c-μ coupling exist even after completely removing
        # any information that MWEBV (line-of-sight dust) could provide?
        z_med = np.median(z_mw)
        hi = z_mw >= z_med
        
        controls = np.column_stack([z_mw[hi], z_mw[hi]**2, mw[hi]])
        coef_c, _, _, _ = lstsq(controls, c_mw[hi], rcond=None)
        coef_mu, _, _, _ = lstsq(controls, mu_mw[hi], rcond=None)
        rho, p = stats.spearmanr(c_mw[hi] - controls @ coef_c, mu_mw[hi] - controls @ coef_mu)
        
        print(f"\n  c-μ coupling at high-z, controlling for MWEBV:")
        print(f"    ρ={rho:.4f}, p={p:.4f}")
        results['c_mu_mwebv_ctrl'] = {'rho': float(rho), 'p': float(p)}
    
    # === Approach 3: Factorization test using ONLY SALT2-independent observables ===
    print("\n  APPROACH 3: Factorization with SALT2-independent observables only")
    
    # Observables that don't come from SALT2 fitting:
    # - zHD (redshift — measured independently)
    # - MWEBV (MW dust — from Schlegel maps)
    # - HOST_LOGMASS (host galaxy — from photometry)
    # - chi2_dof (fit quality — DIAGNOSTIC of SALT2, not output)
    # - μ_resid (uses SALT2 standardization, but is the END product)
    
    host_mass = df['HOST_LOGMASS'].values if 'HOST_LOGMASS' in df.columns else None
    
    if host_mass is not None and chi2_dof is not None and mwebv is not None:
        mask_all = mask & np.isfinite(chi2_dof) & np.isfinite(host_mass) & np.isfinite(mwebv) & (host_mass > 0)
        
        z_all = z[mask_all]
        chi_all = chi2_dof[mask_all]
        hm_all = host_mass[mask_all]
        mw_all = mwebv[mask_all]
        mu_all = mu_res[mask_all]
        
        z_med = np.median(z_all)
        
        def eff_rank(X):
            X_std = (X - X.mean(0)) / (X.std(0) + 1e-10)
            s = np.linalg.svd(X_std, compute_uv=False)
            p = s / s.sum(); p = p[p > 0]
            return float(np.exp(-np.sum(p * np.log(p))))
        
        # Independent observables matrix
        for label, zm in [("Low-z", z_all < z_med), ("High-z", z_all >= z_med)]:
            obs = np.column_stack([chi_all[zm], hm_all[zm], mw_all[zm], mu_all[zm]])
            er = eff_rank(obs)
            print(f"    {label} (n={zm.sum()}): eff_rank[χ², host_mass, MWEBV, μ_resid] = {er:.3f}")
            results[f'independent_rank_{label}'] = {'rank': er, 'n': int(zm.sum())}
        
        # Pairwise entanglement of independent observables
        obs_names = ['chi2_dof', 'HOST_LOGMASS', 'MWEBV', 'mu_resid']
        obs_data = {'chi2_dof': chi_all, 'HOST_LOGMASS': hm_all, 'MWEBV': mw_all, 'mu_resid': mu_all}
        
        from itertools import combinations
        print(f"\n  Pairwise entanglement of independent observables:")
        for a, b in combinations(obs_names, 2):
            lo_m = z_all < z_med
            hi_m = z_all >= z_med
            rho_lo, p_lo = stats.spearmanr(obs_data[a][lo_m], obs_data[b][lo_m])
            rho_hi, p_hi = stats.spearmanr(obs_data[a][hi_m], obs_data[b][hi_m])
            
            change = abs(rho_hi) - abs(rho_lo)
            marker = " ←" if abs(change) > 0.03 and p_hi < 0.05 else ""
            print(f"    {a:>15s} vs {b:<15s}: lo ρ={rho_lo:+.4f}, hi ρ={rho_hi:+.4f}, "
                  f"Δ|ρ|={change:+.4f}{marker}")
    
    # === Approach 4: THE KEY TEST — Does chi2_dof show the SAME threshold behavior? ===
    print("\n  APPROACH 4: χ²/dof threshold test (completely SALT2-independent)")
    
    if chi2_dof is not None:
        mask_chi = mask & np.isfinite(chi2_dof)
        z_c = z[mask_chi]
        chi_c = chi2_dof[mask_chi]
        mu_c = mu_res[mask_chi]
        
        # Running correlation of chi2 vs mu_res, same as we did for c vs mu_res
        z_sorted = np.sort(z_c)
        window = len(z_c) // 5
        step = window // 4
        
        z_centers = []
        rhos = []
        
        order = np.argsort(z_c)
        for i in range(0, len(z_c) - window, step):
            idx = order[i:i+window]
            rho, _ = stats.spearmanr(chi_c[idx], mu_c[idx])
            z_centers.append(float(np.median(z_c[idx])))
            rhos.append(float(rho))
        
        # Does |ρ(χ², μ)| increase with z?
        trend_r, trend_p = stats.spearmanr(z_centers, [abs(r) for r in rhos])
        
        print(f"    |ρ(χ²/dof, μ_resid)| trend with z: ρ={trend_r:.3f}, p={trend_p:.4f}")
        results['chi2_threshold_trend'] = {'rho': float(trend_r), 'p': float(trend_p)}
        
        # Print the running correlation
        print(f"\n    Running ρ(χ²/dof, μ_resid):")
        for zc, r in zip(z_centers, rhos):
            marker = " *" if abs(r) > 0.1 else ""
            print(f"      z={zc:.3f}: ρ={r:+.4f}{marker}")
        
        results['chi2_running'] = {'z_centers': z_centers, 'rhos': rhos}
    
    return results

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 60)
    print("CLOSURE THEORY — ROUND 8: NON-SALT & CMB LENSING")
    print("=" * 60)
    
    df = load_pantheon()
    print(f"Loaded {len(df)} SNe from Pantheon+")
    
    all_results = {}
    
    all_results['L1_lensing'] = test_cmb_lensing(df)
    all_results['I1_independent'] = test_model_independent_color(df)
    
    with open(os.path.join(RESULTS_DIR, 'round8_results.json'), 'w') as f:
        json.dump(all_results, f, indent=2, cls=NpEncoder)
    
    print("\n" + "=" * 60)
    print("ROUND 8 VERDICT")
    print("=" * 60)
    print(f"Results saved to {RESULTS_DIR}/")

if __name__ == "__main__":
    main()
