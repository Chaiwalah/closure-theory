#!/usr/bin/env python3
"""
FRB Cross-Match — Declination-Controlled Reanalysis
====================================================
Tests whether the x1-DM signal from closure_test_frb.py survives
after controlling for CHIME's declination bias.

CHIME only sees dec > -10° with sensitivity peaking at +70-80°.
This creates spatially non-uniform DM sampling that could confound
SN property correlations.
"""

import numpy as np
import os
import json
from scipy import stats
from numpy.linalg import lstsq

np.random.seed(42)

RESULTS_DIR = "results_frb_dec"
os.makedirs(RESULTS_DIR, exist_ok=True)

def load_frb_crossmatch():
    """Load the cross-matched FRB-SN data from the previous run."""
    import pandas as pd
    
    # Load Pantheon+
    pdata = pd.read_csv('data/pantheon_plus.dat', sep=r'\s+', comment='#')
    
    # Load CHIME FRB
    frb_data = pd.read_csv('data/chimefrbcat1.csv', sep='\t')
    frb_data = frb_data.replace(-9999, np.nan).replace('-9999', np.nan)
    
    return pdata, frb_data

def angular_sep(ra1, dec1, ra2, dec2):
    """Angular separation in degrees."""
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    cos_sep = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)
    return np.degrees(np.arccos(np.clip(cos_sep, -1, 1)))

def crossmatch(pdata, frb_data, max_sep=10.0):
    """Cross-match SNe with FRBs within max_sep degrees, assign nearest FRB DM."""
    frb_ra = frb_data['ra'].values.astype(float)
    frb_dec = frb_data['dec'].values.astype(float)
    frb_dm = frb_data['bonsai_dm'].values.astype(float)
    
    valid = np.isfinite(frb_ra) & np.isfinite(frb_dec) & np.isfinite(frb_dm) & (frb_dm > 0)
    frb_ra = frb_ra[valid]
    frb_dec = frb_dec[valid]
    frb_dm = frb_dm[valid]
    
    print(f"Valid FRBs: {len(frb_ra)}")
    
    sn_ra = pdata['RA'].values
    sn_dec = pdata['DEC'].values if 'DEC' in pdata else pdata['DECL'].values
    
    matched = {'sn_idx': [], 'frb_idx': [], 'sep': [], 'dm': [],
               'frb_dec': [], 'sn_dec': []}
    
    for i in range(len(sn_ra)):
        if not (np.isfinite(sn_ra[i]) and np.isfinite(sn_dec[i])):
            continue
        seps = angular_sep(sn_ra[i], sn_dec[i], frb_ra, frb_dec)
        j = np.argmin(seps)
        if seps[j] <= max_sep:
            matched['sn_idx'].append(i)
            matched['frb_idx'].append(j)
            matched['sep'].append(seps[j])
            matched['dm'].append(frb_dm[j])
            matched['frb_dec'].append(frb_dec[j])
            matched['sn_dec'].append(sn_dec[i])
    
    for k in matched:
        matched[k] = np.array(matched[k])
    
    return matched

def main():
    print("=" * 70)
    print("FRB CROSS-MATCH — DECLINATION-CONTROLLED REANALYSIS")
    print("=" * 70)
    
    pdata, frb_data = load_frb_crossmatch()
    matched = crossmatch(pdata, frb_data)
    
    idx = matched['sn_idx'].astype(int)
    dm = matched['dm']
    frb_dec = matched['frb_dec']
    sn_dec = matched['sn_dec']
    sep = matched['sep']
    
    z = pdata['zHD'].values[idx]
    x1 = pdata['x1'].values[idx]
    c = pdata['c'].values[idx]
    
    # Compute Hubble residual
    from scipy.integrate import quad
    
    def mu_lcdm(z_val, H0=70, Om=0.3):
        def integrand(zp):
            return 1.0 / np.sqrt(Om*(1+zp)**3 + (1-Om))
        dc, _ = quad(integrand, 0, z_val)
        dl = (1+z_val) * dc * 2.998e5 / H0
        return 5*np.log10(dl) + 25
    
    mu_model = np.array([mu_lcdm(zi) for zi in z])
    if 'MU_SH0ES' in pdata.columns:
        mu_obs = pdata['MU_SH0ES'].values[idx]
    elif 'm_b_corr' in pdata.columns:
        mu_obs = pdata['m_b_corr'].values[idx]
    else:
        mu_obs = pdata['mB'].values[idx] + 0.15*x1 - 3.1*c - 19.36
    
    mu_res = mu_obs - mu_model
    
    n = len(z)
    print(f"\nMatched sightlines: {n}")
    print(f"Dec range (FRB): {frb_dec.min():.1f}° to {frb_dec.max():.1f}°")
    print(f"Dec range (SN): {sn_dec.min():.1f}° to {sn_dec.max():.1f}°")
    
    results = {}
    
    # --- Original signal (no dec control) ---
    print("\n" + "-" * 50)
    print("ORIGINAL (no declination control)")
    rho_x1, p_x1 = stats.spearmanr(x1, dm)
    rho_c, p_c = stats.spearmanr(c, dm)
    print(f"  x1 vs DM: ρ={rho_x1:.4f}, p={p_x1:.4f}")
    print(f"  c vs DM:  ρ={rho_c:.4f}, p={p_c:.4f}")
    results['original'] = {
        'x1_dm': {'rho': float(rho_x1), 'p': float(p_x1)},
        'c_dm': {'rho': float(rho_c), 'p': float(p_c)},
        'n': n
    }
    
    # --- Declination-controlled partial correlation ---
    print("\n" + "-" * 50)
    print("DECLINATION-CONTROLLED (partial correlation)")
    
    # Control variables: frb_dec, sn_dec, angular sep, z
    controls = np.column_stack([frb_dec, sn_dec, sep, z])
    
    def partial_corr(x, y, C):
        """Partial Spearman correlation controlling for C."""
        # Residualize both x and y against C using OLS
        coef_x, _, _, _ = lstsq(C, x, rcond=None)
        coef_y, _, _, _ = lstsq(C, y, rcond=None)
        x_resid = x - C @ coef_x
        y_resid = y - C @ coef_y
        return stats.spearmanr(x_resid, y_resid)
    
    rho_x1_ctrl, p_x1_ctrl = partial_corr(x1, dm, controls)
    rho_c_ctrl, p_c_ctrl = partial_corr(c, dm, controls)
    rho_mu_ctrl, p_mu_ctrl = partial_corr(mu_res, dm, controls)
    
    print(f"  x1 vs DM (ctrl dec,sep,z): ρ={rho_x1_ctrl:.4f}, p={p_x1_ctrl:.4f}")
    print(f"  c vs DM (ctrl dec,sep,z):  ρ={rho_c_ctrl:.4f}, p={p_c_ctrl:.4f}")
    print(f"  μ_res vs DM (ctrl):        ρ={rho_mu_ctrl:.4f}, p={p_mu_ctrl:.4f}")
    
    results['dec_controlled'] = {
        'x1_dm': {'rho': float(rho_x1_ctrl), 'p': float(p_x1_ctrl)},
        'c_dm': {'rho': float(rho_c_ctrl), 'p': float(p_c_ctrl)},
        'mu_dm': {'rho': float(rho_mu_ctrl), 'p': float(p_mu_ctrl)},
    }
    
    # --- Declination-binned analysis ---
    print("\n" + "-" * 50)
    print("DECLINATION-BINNED ANALYSIS")
    
    dec_bins = np.percentile(sn_dec, [0, 33, 67, 100])
    for i in range(len(dec_bins)-1):
        mask = (sn_dec >= dec_bins[i]) & (sn_dec < dec_bins[i+1] + 0.01)
        if mask.sum() < 20:
            continue
        rho, p = stats.spearmanr(x1[mask], dm[mask])
        print(f"  Dec [{dec_bins[i]:.0f}°, {dec_bins[i+1]:.0f}°] (n={mask.sum()}): x1-DM ρ={rho:.4f}, p={p:.4f}")
    
    # --- Redshift-binned with dec control ---
    print("\n" + "-" * 50)
    print("REDSHIFT BINS WITH DEC CONTROL")
    
    z_threshold = 0.82
    lo = z < z_threshold
    hi = z >= z_threshold
    
    for label, mask in [("z<0.82", lo), ("z≥0.82", hi)]:
        if mask.sum() < 10:
            print(f"  {label}: n={mask.sum()} (too few)")
            continue
        C = np.column_stack([frb_dec[mask], sn_dec[mask], sep[mask]])
        rho, p = partial_corr(x1[mask], dm[mask], C)
        print(f"  {label} (n={mask.sum()}): x1-DM ρ={rho:.4f}, p={p:.4f} [dec-controlled]")
    
    results['z_binned_dec_ctrl'] = {}
    for label, mask in [("z_lt_082", lo), ("z_ge_082", hi)]:
        if mask.sum() >= 10:
            C = np.column_stack([frb_dec[mask], sn_dec[mask], sep[mask]])
            rho, p = partial_corr(x1[mask], dm[mask], C)
            results['z_binned_dec_ctrl'][label] = {'rho': float(rho), 'p': float(p), 'n': int(mask.sum())}
    
    # --- Permutation test with dec control ---
    print("\n" + "-" * 50)
    print("PERMUTATION TEST (dec-controlled)")
    
    # Residualize x1 and dm against dec controls
    coef_x1, _, _, _ = lstsq(controls, x1, rcond=None)
    coef_dm, _, _, _ = lstsq(controls, dm, rcond=None)
    x1_resid = x1 - controls @ coef_x1
    dm_resid = dm - controls @ coef_dm
    
    observed_rho = stats.spearmanr(x1_resid, dm_resid)[0]
    n_perm = 5000
    count = 0
    for _ in range(n_perm):
        perm_rho = stats.spearmanr(np.random.permutation(x1_resid), dm_resid)[0]
        if abs(perm_rho) >= abs(observed_rho):
            count += 1
    p_perm = count / n_perm
    
    print(f"  Observed ρ (residualized): {observed_rho:.4f}")
    print(f"  Permutation p: {p_perm:.4f}")
    
    results['permutation_dec_ctrl'] = {
        'observed_rho': float(observed_rho),
        'p_perm': float(p_perm),
        'n_perm': n_perm
    }
    
    # --- Summary ---
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    x1_killed = p_x1_ctrl > 0.05
    c_still_null = p_c_ctrl > 0.05
    
    if x1_killed:
        print("  ✓ x1-DM signal DISAPPEARS after declination control")
        print("    → Confirms spatial confound hypothesis")
    else:
        print("  ✗ x1-DM signal SURVIVES declination control")
        print("    → May indicate real baryonic coupling")
    
    if c_still_null:
        print("  ✓ c-DM remains NULL (frequency-dependent observable immune)")
    
    print(f"\n  Verdict: {'SPATIAL CONFOUND CONFIRMED' if x1_killed else 'REAL SIGNAL — NEEDS INVESTIGATION'}")
    
    results['verdict'] = 'spatial_confound' if x1_killed else 'real_signal'
    
    with open(os.path.join(RESULTS_DIR, 'frb_dec_results.json'), 'w') as f:
        json.dump(results, f, indent=2, default=float)
    
    print(f"\nResults saved to {RESULTS_DIR}/frb_dec_results.json")

if __name__ == "__main__":
    main()
