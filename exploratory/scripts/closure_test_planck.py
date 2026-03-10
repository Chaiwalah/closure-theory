#!/usr/bin/env python3
"""
Closure Theory — Round 8: Planck CMB Lensing Tests
====================================================
Cross-correlate Planck PR3 lensing convergence (κ) with SN Ia sightlines.

Tests:
  L1: κ at SN positions vs μ_resid — does integrated matter explain residuals?
  L2: κ at SN positions vs color residual — does lensing correlate with color anomaly?
  L3: κ high-z vs low-z split — does κ-coupling turn on at z*≈0.82?
  L4: κ as control variable — does closure signal survive after controlling for κ?
  L5: κ vs factorization — does lensing convergence predict observable entanglement?

Closure prediction:
  - L1: Weak or null (lensing magnification is ~0.01 mag, much smaller than closure signal)
  - L2: NULL (color anomaly is not caused by intervening matter)
  - L3: No z-dependent κ-coupling (lensing is a smooth function of z)
  - L4: Signal SURVIVES (closure is not lensing)
  - L5: NULL (lensing doesn't cause factorization collapse)

Data:
  - Planck PR3 lensing: data/planck/COM_Lensing_4096_R3.00/MV/dat_klm.fits
  - Mean field: MV/mf_klm.fits  
  - Mask: mask.fits.gz
  - Pantheon+: data/pantheon_plus.dat
"""

import numpy as np
import pandas as pd
import healpy as hp
from scipy import stats
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ── Paths ──
BASE = Path(__file__).parent
PLANCK_DIR = BASE / "data" / "planck" / "COM_Lensing_4096_R3.00"
RESULTS_DIR = BASE / "results_planck"
RESULTS_DIR.mkdir(exist_ok=True)

def load_pantheon():
    """Load Pantheon+ with standard cuts."""
    df = pd.read_csv(BASE / "data" / "pantheon_plus.dat", sep=r'\s+', comment='#')
    for col in ['zHD', 'MU_SH0ES', 'MU_SH0ES_ERR_DIAG', 'c', 'x1', 'RA', 'DEC']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna(subset=['zHD', 'MU_SH0ES', 'c', 'x1', 'RA', 'DEC'])
    df = df[df['zHD'] > 0.01].copy()
    
    # Compute μ_resid (subtract simple Hubble law fit)
    from numpy.polynomial import polynomial as P
    log_z = np.log10(df['zHD'].values)
    mu = df['MU_SH0ES'].values
    mask_fit = df['zHD'].values < 0.5  # fit in linear regime
    coeffs = P.polyfit(log_z[mask_fit], mu[mask_fit], 2)
    mu_model = P.polyval(log_z, coeffs)
    df['mu_resid'] = mu - mu_model
    
    # Color residual: remove linear z-trend from c
    slope, intercept, _, _, _ = stats.linregress(df['zHD'], df['c'])
    df['c_resid'] = df['c'] - (slope * df['zHD'] + intercept)
    
    return df


def load_kappa_map():
    """Load Planck PR3 MV lensing convergence κ map.
    
    The data is stored as alm (spherical harmonic coefficients).
    We subtract the mean field, then convert to a real-space HEALPix map.
    """
    # Read lensing alm and mean field
    dat_alm = hp.read_alm(str(PLANCK_DIR / "MV" / "dat_klm.fits"))
    mf_alm = hp.read_alm(str(PLANCK_DIR / "MV" / "mf_klm.fits"))
    
    # Subtract mean field
    kappa_alm = dat_alm - mf_alm
    
    # Convert phi (lensing potential) alm to kappa (convergence) alm
    # κ_lm = l(l+1)/2 * φ_lm
    lmax = hp.Alm.getlmax(len(kappa_alm))
    l_arr = np.arange(lmax + 1)
    fl = l_arr * (l_arr + 1) / 2.0
    fl[0] = 0  # monopole
    kappa_alm = hp.almxfl(kappa_alm, fl)
    
    # Convert to real-space map at NSIDE=2048 (good enough, original is 4096)
    nside = 2048
    kappa_map = hp.alm2map(kappa_alm, nside, verbose=False)
    
    # Load mask
    mask = hp.read_map(str(PLANCK_DIR / "mask.fits.gz"), verbose=False)
    # Downgrade mask to match
    if hp.get_nside(mask) != nside:
        mask = hp.ud_grade(mask, nside)
    
    return kappa_map, mask, nside


def get_kappa_at_positions(kappa_map, mask, nside, ra_deg, dec_deg):
    """Sample κ map at SN positions. Return κ values and validity mask."""
    # Convert RA/DEC to theta/phi (healpy convention)
    theta = np.radians(90.0 - dec_deg)  # colatitude
    phi = np.radians(ra_deg)
    
    # Get pixel indices
    pix = hp.ang2pix(nside, theta, phi)
    
    kappa_vals = kappa_map[pix]
    mask_vals = mask[pix] > 0.5  # valid pixels
    
    return kappa_vals, mask_vals


def test_L1_kappa_vs_mu_resid(df, kappa, valid):
    """L1: Does integrated matter (κ) explain distance residuals?"""
    m = valid & np.isfinite(kappa) & np.isfinite(df['mu_resid'].values)
    k = kappa[m]
    mu_r = df['mu_resid'].values[m]
    
    rho, p = stats.spearmanr(k, mu_r)
    
    # Also partial: control for z
    from numpy.linalg import lstsq
    z = df['zHD'].values[m]
    # Residualize both against z
    A = np.column_stack([z, np.ones(len(z))])
    k_resid = k - A @ lstsq(A, k, rcond=None)[0]
    mu_resid2 = mu_r - A @ lstsq(A, mu_r, rcond=None)[0]
    rho_partial, p_partial = stats.spearmanr(k_resid, mu_resid2)
    
    return {
        'test': 'L1: κ vs μ_resid',
        'N': int(m.sum()),
        'rho': round(rho, 4),
        'p': round(p, 6),
        'rho_partial_z': round(rho_partial, 4),
        'p_partial_z': round(p_partial, 6),
        'prediction': 'Weak or null (lensing magnification ~0.01 mag)',
    }


def test_L2_kappa_vs_c_resid(df, kappa, valid):
    """L2: Does lensing correlate with color anomaly?"""
    m = valid & np.isfinite(kappa) & np.isfinite(df['c_resid'].values)
    k = kappa[m]
    c_r = df['c_resid'].values[m]
    
    rho, p = stats.spearmanr(k, c_r)
    
    # Partial: control for z
    z = df['zHD'].values[m]
    A = np.column_stack([z, np.ones(len(z))])
    from numpy.linalg import lstsq
    k_resid = k - A @ lstsq(A, k, rcond=None)[0]
    c_resid2 = c_r - A @ lstsq(A, c_r, rcond=None)[0]
    rho_partial, p_partial = stats.spearmanr(k_resid, c_resid2)
    
    return {
        'test': 'L2: κ vs c_resid',
        'N': int(m.sum()),
        'rho': round(rho, 4),
        'p': round(p, 6),
        'rho_partial_z': round(rho_partial, 4),
        'p_partial_z': round(p_partial, 6),
        'prediction': 'NULL (color anomaly not caused by intervening matter)',
    }


def test_L3_kappa_z_split(df, kappa, valid):
    """L3: Does κ-coupling turn on at z*≈0.82?"""
    m = valid & np.isfinite(kappa) & np.isfinite(df['c_resid'].values)
    k = kappa[m]
    c_r = df['c_resid'].values[m]
    z = df['zHD'].values[m]
    
    z_thresh = 0.82
    lo = z < z_thresh
    hi = z >= z_thresh
    
    results = {}
    for label, mask_z in [('low_z', lo), ('high_z', hi)]:
        if mask_z.sum() > 10:
            rho, p = stats.spearmanr(k[mask_z], c_r[mask_z])
            results[f'rho_{label}'] = round(rho, 4)
            results[f'p_{label}'] = round(p, 6)
            results[f'N_{label}'] = int(mask_z.sum())
        else:
            results[f'rho_{label}'] = None
            results[f'N_{label}'] = int(mask_z.sum())
    
    # Also test κ vs μ_resid split
    mu_r = df['mu_resid'].values[m]
    for label, mask_z in [('low_z', lo), ('high_z', hi)]:
        if mask_z.sum() > 10:
            rho, p = stats.spearmanr(k[mask_z], mu_r[mask_z])
            results[f'rho_mu_{label}'] = round(rho, 4)
            results[f'p_mu_{label}'] = round(p, 6)
    
    results['test'] = 'L3: κ coupling z-split at z*=0.82'
    results['prediction'] = 'No z-dependent κ-coupling (lensing is smooth function of z)'
    return results


def test_L4_closure_survives_kappa(df, kappa, valid):
    """L4: Does closure signal survive after controlling for κ?
    
    Test: c_resid vs μ_resid correlation at high-z, before and after κ control.
    If closure is real and not lensing, signal should survive.
    """
    z_thresh = 0.82
    m = (valid & np.isfinite(kappa) & 
         np.isfinite(df['c_resid'].values) & 
         np.isfinite(df['mu_resid'].values) &
         (df['zHD'].values >= z_thresh))
    
    if m.sum() < 20:
        return {'test': 'L4: Closure survives κ control', 'N': int(m.sum()), 'status': 'insufficient data'}
    
    k = kappa[m]
    c_r = df['c_resid'].values[m]
    mu_r = df['mu_resid'].values[m]
    z = df['zHD'].values[m]
    
    # Before κ control
    rho_before, p_before = stats.spearmanr(c_r, mu_r)
    
    # After κ control: residualize c_resid and μ_resid against κ and z
    from numpy.linalg import lstsq
    A = np.column_stack([k, z, np.ones(len(k))])
    c_controlled = c_r - A @ lstsq(A, c_r, rcond=None)[0]
    mu_controlled = mu_r - A @ lstsq(A, mu_r, rcond=None)[0]
    
    rho_after, p_after = stats.spearmanr(c_controlled, mu_controlled)
    
    survival = abs(rho_after) / abs(rho_before) * 100 if abs(rho_before) > 0 else None
    
    return {
        'test': 'L4: Closure signal survives κ control (z≥0.82)',
        'N': int(m.sum()),
        'rho_before': round(rho_before, 4),
        'p_before': round(p_before, 6),
        'rho_after_kappa_ctrl': round(rho_after, 4),
        'p_after_kappa_ctrl': round(p_after, 6),
        'survival_pct': round(survival, 1) if survival else None,
        'prediction': 'Signal SURVIVES (closure ≠ lensing)',
    }


def test_L5_kappa_vs_x1(df, kappa, valid):
    """L5: κ vs x1 — lensing shouldn't correlate with stretch either.
    Serves as additional null check."""
    m = valid & np.isfinite(kappa) & np.isfinite(df['x1'].values)
    k = kappa[m]
    x1 = df['x1'].values[m]
    
    rho, p = stats.spearmanr(k, x1)
    
    return {
        'test': 'L5: κ vs x1 (null check)',
        'N': int(m.sum()),
        'rho': round(rho, 4),
        'p': round(p, 6),
        'prediction': 'NULL (stretch immune to both lensing and closure)',
    }


def main():
    print("=" * 70)
    print("CLOSURE THEORY — ROUND 8: PLANCK CMB LENSING TESTS")
    print("=" * 70)
    
    print("\nLoading Pantheon+ data...")
    df = load_pantheon()
    print(f"  → {len(df)} SNe after cuts")
    
    print("\nLoading Planck PR3 lensing convergence map...")
    kappa_map, mask, nside = load_kappa_map()
    print(f"  → NSIDE={nside}, {np.sum(mask > 0.5)} valid pixels")
    
    print("\nSampling κ at SN positions...")
    kappa_vals, valid = get_kappa_at_positions(
        kappa_map, mask, nside, df['RA'].values, df['DEC'].values
    )
    print(f"  → {valid.sum()}/{len(df)} SNe in valid lensing region")
    print(f"  → κ range: [{kappa_vals[valid].min():.4f}, {kappa_vals[valid].max():.4f}]")
    print(f"  → κ mean: {kappa_vals[valid].mean():.6f}, std: {kappa_vals[valid].std():.4f}")
    
    # Reset index for alignment
    df = df.reset_index(drop=True)
    
    all_results = []
    
    # Run tests
    tests = [
        ("L1", test_L1_kappa_vs_mu_resid),
        ("L2", test_L2_kappa_vs_c_resid),
        ("L3", test_L3_kappa_z_split),
        ("L4", test_L4_closure_survives_kappa),
        ("L5", test_L5_kappa_vs_x1),
    ]
    
    for tid, func in tests:
        print(f"\n{'─' * 50}")
        print(f"Running {tid}...")
        result = func(df, kappa_vals, valid)
        all_results.append(result)
        for k, v in result.items():
            print(f"  {k}: {v}")
    
    # Save results
    with open(RESULTS_DIR / "planck_results.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY — ROUND 8: PLANCK CMB LENSING")
    print(f"{'=' * 70}")
    for r in all_results:
        name = r['test']
        if 'p' in r:
            print(f"  {name}: ρ={r.get('rho','?')}, p={r.get('p','?')}")
        elif 'survival_pct' in r:
            print(f"  {name}: survival={r.get('survival_pct','?')}%")
        else:
            print(f"  {name}: see details above")
    
    print(f"\nResults saved to {RESULTS_DIR}/")
    print("Done.")


if __name__ == "__main__":
    main()
