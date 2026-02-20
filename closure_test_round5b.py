#!/usr/bin/env python3
"""
Closure Theory — Round 5B: Anisotropy Stress Tests
====================================================
Rigorous null tests for the S1 anisotropy detection.

Tests:
  A) Stratified permutation null (survey × z-bin × SNR-bin × host-mass-bin)
  B) Spherical harmonic model selection (ℓ=0,1,2 with AIC/BIC + cross-validation)
  C) Anisotropy by redshift regime (low-z weak, high-z strong?)
  D) CMB dipole alignment — the make-or-break figure
"""

import numpy as np
import os
import json
from scipy import stats
from scipy.optimize import curve_fit, minimize
from itertools import combinations
from numpy.linalg import lstsq
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

RESULTS_DIR = "results_round5b"
os.makedirs(RESULTS_DIR, exist_ok=True)

# =============================================================================
# DATA LOADING
# =============================================================================
def load_and_prepare():
    import pandas as pd
    df = pd.read_csv("data/pantheon_plus.dat", sep=r'\s+', comment='#')
    
    z = df['zHD'].values
    c = df['c'].values
    x1 = df['x1'].values
    
    # Hubble residual
    from scipy.integrate import quad
    def mu_lcdm(zv, H0=73.04, Om=0.334):
        if not np.isfinite(zv) or zv <= 0: return np.nan
        dc, _ = quad(lambda zp: 1/np.sqrt(Om*(1+zp)**3+(1-Om)), 0, zv)
        return 5*np.log10((1+zv)*dc*2.998e5/H0) + 25
    
    mu_model = np.array([mu_lcdm(zi) for zi in z])
    
    if 'MU_SH0ES' in df.columns:
        mu_res = df['MU_SH0ES'].values - mu_model
    else:
        mu_res = df['m_b_corr'].values - mu_model
    
    # Color residual
    mask_fit = np.isfinite(c) & np.isfinite(z)
    poly = np.polyfit(z[mask_fit], c[mask_fit], 1)
    c_resid = c - np.polyval(poly, z)
    
    # Survey ID
    survey = df['IDSURVEY'].values if 'IDSURVEY' in df.columns else np.ones(len(z))
    
    # Host mass
    host_mass = df['HOST_LOGMASS'].values if 'HOST_LOGMASS' in df.columns else np.full(len(z), 10.0)
    
    # SNR proxy (use MUERR inverse as SNR proxy)
    if 'MU_SH0ES_ERR_DIAG' in df.columns:
        snr_proxy = 1.0 / (df['MU_SH0ES_ERR_DIAG'].values + 0.001)
    elif 'm_b_corr_err_DIAG' in df.columns:
        snr_proxy = 1.0 / (df['m_b_corr_err_DIAG'].values + 0.001)
    else:
        snr_proxy = np.ones(len(z)) * 10
    
    # RA/DEC
    ra = df['RA'].values
    dec = df['DEC'].values
    
    # chi2/dof
    chi2_dof = np.full(len(z), np.nan)
    if 'NDOF' in df.columns and 'FITCHI2' in df.columns:
        ndof = df['NDOF'].values
        valid = ndof > 0
        chi2_dof[valid] = df['FITCHI2'].values[valid] / ndof[valid]
    
    # Quality mask
    mask = (np.isfinite(z) & np.isfinite(c) & np.isfinite(x1) & 
            np.isfinite(mu_res) & np.isfinite(ra) & np.isfinite(dec) &
            (z > 0.01) & (z < 2.5))
    
    out = {
        'z': z[mask], 'c': c[mask], 'x1': x1[mask], 'mu_res': mu_res[mask],
        'c_resid': c_resid[mask], 'ra': ra[mask], 'dec': dec[mask],
        'survey': survey[mask], 'host_mass': host_mass[mask],
        'snr_proxy': snr_proxy[mask], 'chi2_dof': chi2_dof[mask],
    }
    print(f"Loaded {mask.sum()} SNe after quality cuts")
    print(f"z range: {out['z'].min():.3f} — {out['z'].max():.3f}")
    return out

def partial_r_controlled(c, mu, z, x1, host_mass):
    """Compute partial Spearman r(c, mu | z, z², x1, host_mass)."""
    hm = host_mass.copy()
    hm[~np.isfinite(hm)] = np.nanmedian(hm)
    controls = np.column_stack([z, z**2, x1, hm])
    
    coef_c, _, _, _ = lstsq(controls, c, rcond=None)
    coef_mu, _, _, _ = lstsq(controls, mu, rcond=None)
    c_ctrl = c - controls @ coef_c
    mu_ctrl = mu - controls @ coef_mu
    
    rho, p = stats.spearmanr(c_ctrl, mu_ctrl)
    return rho, p

# =============================================================================
# CMB dipole direction (from Planck)
# =============================================================================
CMB_DIPOLE_RA = 167.99   # degrees (l=264.02°, b=48.25° → RA~168°, Dec~-7°)
CMB_DIPOLE_DEC = -6.98

def cos_angle_from_cmb(ra, dec):
    """Cosine of angle between each SN and the CMB dipole direction."""
    ra1 = np.radians(CMB_DIPOLE_RA)
    dec1 = np.radians(CMB_DIPOLE_DEC)
    ra2 = np.radians(ra)
    dec2 = np.radians(dec)
    return np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)

# =============================================================================
# A) STRATIFIED PERMUTATION NULL
# =============================================================================
def test_stratified_permutation(d):
    """Shuffle c within (survey × z-bin × SNR-bin × host-mass-bin), recompute anisotropy."""
    print("\n" + "=" * 60)
    print("A) STRATIFIED PERMUTATION NULL")
    print("=" * 60)
    
    z, c, mu_res = d['z'], d['c'], d['mu_res']
    ra, dec = d['ra'], d['dec']
    survey, host_mass, snr = d['survey'], d['host_mass'], d['snr_proxy']
    x1 = d['x1']
    
    # Create bins
    z_bins = np.digitize(z, np.percentile(z, [0, 25, 50, 75, 100]))
    snr_bins = np.digitize(snr, np.percentile(snr, [0, 33, 67, 100]))
    hm = host_mass.copy()
    hm[~np.isfinite(hm)] = np.nanmedian(hm)
    hm_bins = np.digitize(hm, [0, 10, 11, 20])
    
    # Strata = survey × z_bin × snr_bin × hm_bin
    strata = survey * 10000 + z_bins * 100 + snr_bins * 10 + hm_bins
    unique_strata = np.unique(strata)
    
    n_strata = len(unique_strata)
    strata_sizes = [np.sum(strata == s) for s in unique_strata]
    print(f"  Strata: {n_strata}, sizes: min={min(strata_sizes)}, median={int(np.median(strata_sizes))}, max={max(strata_sizes)}")
    
    # Compute real anisotropy metric: variance of partial-r across sky patches
    cos_theta = cos_angle_from_cmb(ra, dec)
    
    def compute_anisotropy_metric(c_vals):
        """Compute partial-r in CMB-aligned hemispheres."""
        toward = cos_theta > 0
        away = cos_theta <= 0
        
        r_toward, _ = partial_r_controlled(c_vals[toward], mu_res[toward], 
                                            z[toward], x1[toward], host_mass[toward])
        r_away, _ = partial_r_controlled(c_vals[away], mu_res[away],
                                          z[away], x1[away], host_mass[away])
        return abs(r_toward - r_away)  # anisotropy = difference between hemispheres
    
    # Also: variance of partial-r across RA patches
    def compute_patch_variance(c_vals):
        """Variance of partial-r across 6 RA patches."""
        rs = []
        ra_edges = np.linspace(0, 360, 7)
        for i in range(6):
            m = (ra >= ra_edges[i]) & (ra < ra_edges[i+1])
            if m.sum() > 50:
                r, _ = partial_r_controlled(c_vals[m], mu_res[m], z[m], x1[m], host_mass[m])
                rs.append(r)
        return np.var(rs) if len(rs) > 2 else 0
    
    real_hemi = compute_anisotropy_metric(c)
    real_var = compute_patch_variance(c)
    print(f"\n  Real hemisphere asymmetry: {real_hemi:.4f}")
    print(f"  Real patch variance: {real_var:.6f}")
    
    # Permutation test
    n_perm = 2000
    perm_hemi = []
    perm_var = []
    
    for p_idx in range(n_perm):
        c_shuffled = c.copy()
        for s in unique_strata:
            idx = np.where(strata == s)[0]
            if len(idx) > 1:
                c_shuffled[idx] = np.random.permutation(c_shuffled[idx])
        
        perm_hemi.append(compute_anisotropy_metric(c_shuffled))
        perm_var.append(compute_patch_variance(c_shuffled))
    
    perm_hemi = np.array(perm_hemi)
    perm_var = np.array(perm_var)
    
    p_hemi = float(np.mean(perm_hemi >= real_hemi))
    p_var = float(np.mean(perm_var >= real_var))
    
    print(f"\n  Permutation results ({n_perm} perms, stratified by survey×z×SNR×host_mass):")
    print(f"    Hemisphere asymmetry: p = {p_hemi:.4f}")
    print(f"    Patch variance: p = {p_var:.4f}")
    
    if p_hemi < 0.05:
        print("    → Hemisphere asymmetry SURVIVES stratified null ✓")
    else:
        print("    → Hemisphere asymmetry consistent with null")
    
    if p_var < 0.05:
        print("    → Patch variance SURVIVES stratified null ✓")
    else:
        print("    → Patch variance consistent with null")
    
    return {
        'real_hemisphere_asymmetry': float(real_hemi),
        'real_patch_variance': float(real_var),
        'p_hemisphere': p_hemi,
        'p_patch_variance': p_var,
        'n_perm': n_perm,
        'n_strata': n_strata,
        'perm_hemi_95': float(np.percentile(perm_hemi, 95)),
        'perm_var_95': float(np.percentile(perm_var, 95)),
    }

# =============================================================================
# B) SPHERICAL HARMONIC MODEL SELECTION
# =============================================================================
def test_spherical_harmonics(d):
    """Fit ℓ=0,1,2 spherical harmonics to closure strength map, compare AIC/BIC."""
    print("\n" + "=" * 60)
    print("B) SPHERICAL HARMONIC MODEL SELECTION")
    print("=" * 60)
    
    z, c, mu_res = d['z'], d['c'], d['mu_res']
    ra, dec = d['ra'], d['dec']
    x1, host_mass = d['x1'], d['host_mass']
    
    # Only high-z (where closure lives)
    hi = z > np.median(z)
    
    # Compute closure strength in sky pixels
    # Use HEALPix-like equal-area tiling: 12 base pixels
    # Approximate with RA/DEC grid
    n_ra_bins, n_dec_bins = 6, 3
    ra_edges = np.linspace(0, 360, n_ra_bins + 1)
    dec_edges = np.linspace(-90, 90, n_dec_bins + 1)
    
    pixel_ra = []
    pixel_dec = []
    pixel_r = []
    pixel_n = []
    pixel_r_err = []
    
    for i in range(n_ra_bins):
        for j in range(n_dec_bins):
            m = hi & (ra >= ra_edges[i]) & (ra < ra_edges[i+1]) & (dec >= dec_edges[j]) & (dec < dec_edges[j+1])
            if m.sum() > 30:
                r, p = partial_r_controlled(c[m], mu_res[m], z[m], x1[m], host_mass[m])
                pixel_ra.append((ra_edges[i] + ra_edges[i+1]) / 2)
                pixel_dec.append((dec_edges[j] + dec_edges[j+1]) / 2)
                pixel_r.append(r)
                pixel_n.append(m.sum())
                # Approximate error on r
                pixel_r_err.append(1.0 / np.sqrt(m.sum() - 3) if m.sum() > 3 else 1.0)
    
    pixel_ra = np.array(pixel_ra)
    pixel_dec = np.array(pixel_dec)
    pixel_r = np.array(pixel_r)
    pixel_n = np.array(pixel_n)
    pixel_r_err = np.array(pixel_r_err)
    
    n_pixels = len(pixel_r)
    print(f"  Sky pixels with >30 high-z SNe: {n_pixels}")
    print(f"  Pixel r values: {pixel_r}")
    
    # Convert to unit vectors
    ra_rad = np.radians(pixel_ra)
    dec_rad = np.radians(pixel_dec)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    zz = np.sin(dec_rad)
    
    # Real spherical harmonics basis functions
    # ℓ=0: Y00 = 1
    # ℓ=1: Y1m = x, y, z  (3 params)
    # ℓ=2: Y2m = xy, xz, yz, x²-y², 2z²-x²-y²  (5 params)
    
    Y0 = np.ones((n_pixels, 1))
    Y1 = np.column_stack([x, y, zz])
    Y2 = np.column_stack([x*y, x*zz, y*zz, x**2 - y**2, 2*zz**2 - x**2 - y**2])
    
    results = {}
    
    for label, basis in [("L=0 (monopole)", Y0), 
                          ("L=0+1 (dipole)", np.column_stack([Y0, Y1])),
                          ("L=0+1+2 (quadrupole)", np.column_stack([Y0, Y1, Y2]))]:
        k = basis.shape[1]  # number of parameters
        
        if n_pixels <= k:
            print(f"  {label}: Too few pixels ({n_pixels}) for {k} params")
            continue
        
        # Weighted least squares
        W = np.diag(1.0 / pixel_r_err**2)
        coef = np.linalg.solve(basis.T @ W @ basis, basis.T @ W @ pixel_r)
        residuals = pixel_r - basis @ coef
        ss_res = np.sum(residuals**2 / pixel_r_err**2)
        
        # Log-likelihood (Gaussian)
        ll = -0.5 * ss_res - 0.5 * n_pixels * np.log(2 * np.pi) - np.sum(np.log(pixel_r_err))
        
        aic = 2 * k - 2 * ll
        bic = k * np.log(n_pixels) - 2 * ll
        
        results[label] = {
            'k': k, 'ss_res': float(ss_res), 'll': float(ll),
            'aic': float(aic), 'bic': float(bic),
            'r_squared': float(1 - np.sum(residuals**2) / np.sum((pixel_r - pixel_r.mean())**2))
        }
        
        print(f"  {label}: k={k}, R²={results[label]['r_squared']:.3f}, AIC={aic:.1f}, BIC={bic:.1f}")
    
    # Model comparison
    if len(results) >= 2:
        models = list(results.keys())
        best_aic = min(models, key=lambda m: results[m]['aic'])
        best_bic = min(models, key=lambda m: results[m]['bic'])
        results['best_aic'] = best_aic
        results['best_bic'] = best_bic
        print(f"\n  Best by AIC: {best_aic}")
        print(f"  Best by BIC: {best_bic}")
    
    # Leave-one-out cross-validation
    print("\n  Leave-one-out cross-validation:")
    for label, basis in [("L=0 (monopole)", Y0),
                          ("L=0+1 (dipole)", np.column_stack([Y0, Y1])),
                          ("L=0+1+2 (quadrupole)", np.column_stack([Y0, Y1, Y2]))]:
        k = basis.shape[1]
        if n_pixels <= k + 1:
            continue
        
        cv_errors = []
        for leave_out in range(n_pixels):
            train_mask = np.ones(n_pixels, dtype=bool)
            train_mask[leave_out] = False
            
            B_train = basis[train_mask]
            r_train = pixel_r[train_mask]
            
            try:
                coef_cv = np.linalg.lstsq(B_train, r_train, rcond=None)[0]
                pred = basis[leave_out] @ coef_cv
                cv_errors.append((pixel_r[leave_out] - pred)**2)
            except:
                pass
        
        cv_mse = np.mean(cv_errors)
        results[f'cv_mse_{label}'] = float(cv_mse)
        print(f"    {label}: CV-MSE = {cv_mse:.6f}")
    
    return results

# =============================================================================
# C) ANISOTROPY BY REDSHIFT REGIME
# =============================================================================
def test_anisotropy_by_z(d):
    """Does anisotropy appear only where closure is strong?"""
    print("\n" + "=" * 60)
    print("C) ANISOTROPY BY REDSHIFT REGIME")
    print("=" * 60)
    
    z, c, mu_res = d['z'], d['c'], d['mu_res']
    ra, dec = d['ra'], d['dec']
    x1, host_mass = d['x1'], d['host_mass']
    
    cos_theta = cos_angle_from_cmb(ra, dec)
    
    results = {}
    
    # Split at multiple z boundaries
    z_splits = [np.percentile(z, 25), np.median(z), np.percentile(z, 75)]
    
    for z_split in z_splits:
        lo = z < z_split
        hi = z >= z_split
        
        for label, mask in [(f"z<{z_split:.2f}", lo), (f"z≥{z_split:.2f}", hi)]:
            if mask.sum() < 100:
                continue
            
            # CMB hemisphere split
            toward = mask & (cos_theta > 0)
            away = mask & (cos_theta <= 0)
            
            if toward.sum() > 30 and away.sum() > 30:
                r_toward, p_toward = partial_r_controlled(c[toward], mu_res[toward], 
                                                          z[toward], x1[toward], host_mass[toward])
                r_away, p_away = partial_r_controlled(c[away], mu_res[away],
                                                      z[away], x1[away], host_mass[away])
                asym = abs(r_toward - r_away)
                
                entry = {
                    'r_toward': float(r_toward), 'p_toward': float(p_toward),
                    'r_away': float(r_away), 'p_away': float(p_away),
                    'asymmetry': float(asym),
                    'n_toward': int(toward.sum()), 'n_away': int(away.sum()),
                }
                results[label] = entry
                
                flag = " ←" if asym > 0.05 else ""
                print(f"  {label:>12s}: toward ρ={r_toward:+.4f} (p={p_toward:.3f}), "
                      f"away ρ={r_away:+.4f} (p={p_away:.3f}), asymmetry={asym:.4f}{flag}")
    
    # Key check: is asymmetry stronger at high-z?
    z_med = np.median(z)
    lo_key = f"z<{z_med:.2f}"
    hi_key = f"z≥{z_med:.2f}"
    
    if lo_key in results and hi_key in results:
        asym_lo = results[lo_key]['asymmetry']
        asym_hi = results[hi_key]['asymmetry']
        results['closure_regime_check'] = {
            'low_z_asymmetry': asym_lo,
            'high_z_asymmetry': asym_hi,
            'ratio': asym_hi / (asym_lo + 1e-10),
            'stronger_at_high_z': asym_hi > asym_lo,
        }
        print(f"\n  Low-z asymmetry: {asym_lo:.4f}")
        print(f"  High-z asymmetry: {asym_hi:.4f}")
        
        if asym_hi > 2 * asym_lo:
            print("  → Anisotropy TURNS ON in closure regime ✓")
        elif asym_hi > asym_lo:
            print("  → Anisotropy stronger at high-z (moderate)")
        else:
            print("  → ⚠️ Anisotropy present at low-z too (systematics concern)")
    
    return results

# =============================================================================
# D) CMB DIPOLE ALIGNMENT — THE MAKE-OR-BREAK FIGURE
# =============================================================================
def test_cmb_dipole_alignment(d):
    """The money plot: closure strength vs angle from CMB dipole."""
    print("\n" + "=" * 60)
    print("D) CMB DIPOLE ALIGNMENT — THE MAKE-OR-BREAK")
    print("=" * 60)
    
    z, c, mu_res = d['z'], d['c'], d['mu_res']
    ra, dec = d['ra'], d['dec']
    x1, host_mass = d['x1'], d['host_mass']
    survey, snr = d['survey'], d['snr_proxy']
    
    cos_theta = cos_angle_from_cmb(ra, dec)
    z_med = np.median(z)
    
    results = {'high_z': {}, 'low_z': {}, 'permutation': {}}
    
    # Bin by cos(theta) from CMB dipole
    n_bins = 8
    cos_edges = np.linspace(-1, 1, n_bins + 1)
    
    for regime, z_mask, regime_label in [("high_z", z >= z_med, "HIGH-Z (closure regime)"),
                                          ("low_z", z < z_med, "LOW-Z (no closure)")]:
        print(f"\n  {regime_label}:")
        bin_centers = []
        bin_r = []
        bin_n = []
        
        for i in range(n_bins):
            m = z_mask & (cos_theta >= cos_edges[i]) & (cos_theta < cos_edges[i+1])
            if m.sum() > 30:
                r, p = partial_r_controlled(c[m], mu_res[m], z[m], x1[m], host_mass[m])
                cc = (cos_edges[i] + cos_edges[i+1]) / 2
                bin_centers.append(float(cc))
                bin_r.append(float(r))
                bin_n.append(int(m.sum()))
                print(f"    cos(θ)=[{cos_edges[i]:+.2f},{cos_edges[i+1]:+.2f}] "
                      f"(n={m.sum():4d}): partial r = {r:+.4f}")
        
        results[regime]['bin_centers'] = bin_centers
        results[regime]['bin_r'] = bin_r
        results[regime]['bin_n'] = bin_n
        
        # Trend: does closure strength correlate with CMB dipole direction?
        if len(bin_centers) > 3:
            trend_r, trend_p = stats.spearmanr(bin_centers, bin_r)
            results[regime]['trend_rho'] = float(trend_r)
            results[regime]['trend_p'] = float(trend_p)
            print(f"    Trend with CMB dipole: ρ={trend_r:.4f}, p={trend_p:.4f}")
    
    # PERMUTATION BANDS for high-z
    print(f"\n  Computing permutation bands (stratified)...")
    hi_mask = z >= z_med
    
    # Create strata
    z_bins_perm = np.digitize(z, np.percentile(z[hi_mask], [0, 25, 50, 75, 100]))
    snr_bins_perm = np.digitize(snr, np.percentile(snr, [0, 33, 67, 100]))
    hm = host_mass.copy()
    hm[~np.isfinite(hm)] = np.nanmedian(hm)
    hm_bins_perm = np.digitize(hm, [0, 10, 11, 20])
    strata = survey * 10000 + z_bins_perm * 100 + snr_bins_perm * 10 + hm_bins_perm
    unique_strata = np.unique(strata)
    
    n_perm = 1000
    perm_r_by_bin = {i: [] for i in range(n_bins)}
    
    for _ in range(n_perm):
        c_shuf = c.copy()
        for s in unique_strata:
            idx = np.where(strata == s)[0]
            if len(idx) > 1:
                c_shuf[idx] = np.random.permutation(c_shuf[idx])
        
        for i in range(n_bins):
            m = hi_mask & (cos_theta >= cos_edges[i]) & (cos_theta < cos_edges[i+1])
            if m.sum() > 30:
                r, _ = partial_r_controlled(c_shuf[m], mu_res[m], z[m], x1[m], host_mass[m])
                perm_r_by_bin[i].append(r)
    
    # Compute bands
    perm_bands = []
    exceedances = []
    print(f"\n  Permutation bands vs real (high-z):")
    real_bins = results['high_z']['bin_r']
    real_centers = results['high_z']['bin_centers']
    
    for i, (center, real_r) in enumerate(zip(real_centers, real_bins)):
        if len(perm_r_by_bin[i]) > 0:
            perm_arr = np.array(perm_r_by_bin[i])
            p2_5 = np.percentile(perm_arr, 2.5)
            p97_5 = np.percentile(perm_arr, 97.5)
            outside = real_r < p2_5 or real_r > p97_5
            exceedances.append(outside)
            perm_bands.append({'center': center, 'lo': float(p2_5), 'hi': float(p97_5), 
                              'real': real_r, 'outside': bool(outside)})
            marker = " ***" if outside else ""
            print(f"    cos(θ)={center:+.2f}: real={real_r:+.4f}, "
                  f"95% band=[{p2_5:+.4f}, {p97_5:+.4f}]{marker}")
    
    results['permutation']['bands'] = perm_bands
    results['permutation']['bins_outside_95'] = int(sum(exceedances))
    results['permutation']['total_bins'] = len(exceedances)
    
    print(f"\n  Bins outside 95% permutation band: {sum(exceedances)}/{len(exceedances)}")
    
    return results

# =============================================================================
# E) ROBUST STATISTICS CHECK
# =============================================================================
def test_robust_metrics(d):
    """Verify anisotropy with multiple coupling metrics."""
    print("\n" + "=" * 60)
    print("E) ROBUSTNESS: MULTIPLE COUPLING METRICS")
    print("=" * 60)
    
    z, c, mu_res = d['z'], d['c'], d['mu_res']
    ra, dec = d['ra'], d['dec']
    x1, host_mass = d['x1'], d['host_mass']
    
    cos_theta = cos_angle_from_cmb(ra, dec)
    z_med = np.median(z)
    hi = z >= z_med
    
    toward = hi & (cos_theta > 0)
    away = hi & (cos_theta <= 0)
    
    results = {}
    
    # Metric 1: Partial Spearman r
    r_t, _ = partial_r_controlled(c[toward], mu_res[toward], z[toward], x1[toward], host_mass[toward])
    r_a, _ = partial_r_controlled(c[away], mu_res[away], z[away], x1[away], host_mass[away])
    results['partial_r'] = {'toward': float(r_t), 'away': float(r_a), 'diff': float(abs(r_t - r_a))}
    print(f"  Partial Spearman r: toward={r_t:+.4f}, away={r_a:+.4f}, Δ={abs(r_t-r_a):.4f}")
    
    # Metric 2: Explained variance (R²) of c→μ_res regression
    from numpy.linalg import lstsq
    for label, mask in [("toward", toward), ("away", away)]:
        hm = host_mass[mask].copy()
        hm[~np.isfinite(hm)] = np.nanmedian(hm)
        X = np.column_stack([c[mask], z[mask], z[mask]**2, x1[mask], hm])
        coef, _, _, _ = lstsq(X, mu_res[mask], rcond=None)
        pred = X @ coef
        ss_res = np.sum((mu_res[mask] - pred)**2)
        ss_tot = np.sum((mu_res[mask] - mu_res[mask].mean())**2)
        r2 = 1 - ss_res / ss_tot
        
        X_noc = np.column_stack([z[mask], z[mask]**2, x1[mask], hm])
        coef_noc, _, _, _ = lstsq(X_noc, mu_res[mask], rcond=None)
        pred_noc = X_noc @ coef_noc
        ss_res_noc = np.sum((mu_res[mask] - pred_noc)**2)
        r2_noc = 1 - ss_res_noc / ss_tot
        
        delta_r2 = r2 - r2_noc  # how much variance color explains beyond controls
        results[f'r2_{label}'] = float(delta_r2)
    
    print(f"  ΔR² (color contribution): toward={results['r2_toward']:.6f}, away={results['r2_away']:.6f}")
    
    # Metric 3: Mutual information (bootstrap)
    def bootstrap_mi(x, y, n_boot=500, n_bins=10):
        """Bootstrapped discretized mutual information."""
        mi_vals = []
        for _ in range(n_boot):
            idx = np.random.choice(len(x), len(x), replace=True)
            xb, yb = x[idx], y[idx]
            # Discretize
            xd = np.digitize(xb, np.linspace(xb.min(), xb.max(), n_bins))
            yd = np.digitize(yb, np.linspace(yb.min(), yb.max(), n_bins))
            # Joint and marginal counts
            joint = np.zeros((n_bins+1, n_bins+1))
            for xi, yi in zip(xd, yd):
                joint[xi, yi] += 1
            joint /= joint.sum()
            px = joint.sum(axis=1)
            py = joint.sum(axis=0)
            
            mi = 0
            for a in range(n_bins+1):
                for b in range(n_bins+1):
                    if joint[a,b] > 0 and px[a] > 0 and py[b] > 0:
                        mi += joint[a,b] * np.log2(joint[a,b] / (px[a] * py[b]))
            mi_vals.append(mi)
        return np.mean(mi_vals), np.std(mi_vals)
    
    mi_t, mi_t_err = bootstrap_mi(c[toward], mu_res[toward])
    mi_a, mi_a_err = bootstrap_mi(c[away], mu_res[away])
    results['mi'] = {
        'toward': float(mi_t), 'toward_err': float(mi_t_err),
        'away': float(mi_a), 'away_err': float(mi_a_err)
    }
    print(f"  Mutual info: toward={mi_t:.4f}±{mi_t_err:.4f}, away={mi_a:.4f}±{mi_a_err:.4f}")
    
    # Summary
    metrics_agree = (
        (results['partial_r']['toward'] > results['partial_r']['away']) ==
        (results['r2_toward'] > results['r2_away']) ==
        (mi_t > mi_a)
    )
    results['metrics_agree'] = bool(metrics_agree)
    
    stronger_dir = "toward CMB dipole" if abs(r_t) > abs(r_a) else "away from CMB dipole"
    print(f"\n  Stronger coupling direction: {stronger_dir}")
    print(f"  All metrics agree on direction: {metrics_agree}")
    
    return results

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 60)
    print("CLOSURE THEORY — ROUND 5B: ANISOTROPY STRESS TESTS")
    print("=" * 60)
    
    d = load_and_prepare()
    all_results = {}
    
    # A) Stratified permutation
    all_results['A_permutation'] = test_stratified_permutation(d)
    
    # B) Spherical harmonics
    all_results['B_harmonics'] = test_spherical_harmonics(d)
    
    # C) Anisotropy by z regime
    all_results['C_z_regime'] = test_anisotropy_by_z(d)
    
    # D) CMB dipole alignment
    all_results['D_cmb_dipole'] = test_cmb_dipole_alignment(d)
    
    # E) Robust metrics
    all_results['E_robust'] = test_robust_metrics(d)
    
    # Save
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.bool_,)): return bool(obj)
            if isinstance(obj, (np.integer_,)): return int(obj)
            if isinstance(obj, (np.floating,)): return float(obj)
            if isinstance(obj, np.ndarray): return obj.tolist()
            return super().default(obj)
    
    with open(os.path.join(RESULTS_DIR, 'round5b_results.json'), 'w') as f:
        json.dump(all_results, f, indent=2, cls=NpEncoder)
    
    # VERDICT
    print("\n" + "=" * 60)
    print("ROUND 5B VERDICT")
    print("=" * 60)
    
    a = all_results['A_permutation']
    print(f"  A) Stratified perm: hemisphere p={a['p_hemisphere']:.4f}, patch p={a['p_patch_variance']:.4f}")
    
    b = all_results['B_harmonics']
    print(f"  B) Best model: AIC→{b.get('best_aic','?')}, BIC→{b.get('best_bic','?')}")
    
    c_res = all_results['C_z_regime']
    if 'closure_regime_check' in c_res:
        cr = c_res['closure_regime_check']
        print(f"  C) High-z asym={cr['high_z_asymmetry']:.4f} vs low-z={cr['low_z_asymmetry']:.4f} "
              f"(ratio={cr['ratio']:.1f}×)")
    
    d_res = all_results['D_cmb_dipole']
    dp = d_res['permutation']
    print(f"  D) CMB dipole: {dp['bins_outside_95']}/{dp['total_bins']} bins outside 95% band")
    
    e = all_results['E_robust']
    print(f"  E) Metrics agree: {e['metrics_agree']}")
    
    print(f"\nResults saved to {RESULTS_DIR}/")

if __name__ == "__main__":
    main()
