#!/usr/bin/env python3
"""
closure_orientation_test.py — Orientation / Continuum Anisotropy Kill Test
==========================================================================

GPT's challenge: the "indestructible" Hβ divergence might be continuum
anisotropy + orientation bias in a flux-limited survey, not cosmology.

FWHM split ≈ inclination split (FWHM ~ v_vir * sin(θ))
EW = F_line / F_cont, if F_cont is anisotropic (thin disk cos(θ)),
then EW becomes orientation-dependent.

TWO KILL TESTS:

Test 1: EW RATIO CANCELLATION
  Take ratios like EW(MgII)/EW(Hβ) or EW(CIV)/EW(CIII).
  Ratios share the same continuum → anisotropy cancels.
  - If B(z) COLLAPSES in ratios → continuum/orientation is the driver
  - If B(z) SURVIVES in ratios → real line physics / epoch effect

Test 2: LINE DISPERSION vs FWHM SPLIT
  Line dispersion (sigma) is less inclination-dependent than FWHM.
  If B(z) disappears when splitting by sigma instead of FWHM,
  then the signal was orientation-coupled.

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr
from astropy.io import fits
from numpy.polynomial import polynomial as P
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

print("Loading DR16Q catalog...")
f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z_all = d['Z_DR16Q']
lbol = d['LOGLBOL']

# ============================================================
# HELPER
# ============================================================
def compute_B_trend_generic(values, z_vals, split_vals, split_median, n_bins=8):
    """Compute B(z) for arbitrary values split by split_vals at split_median."""
    if len(values) < 2000:
        return None
    
    n_bins = min(n_bins, len(values) // 3000)
    if n_bins < 3:
        n_bins = 3
    z_edges = np.percentile(z_vals, np.linspace(0, 100, n_bins + 1))
    z_edges = np.unique(np.round(z_edges, 3))
    
    z_list, B_list = [], []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_mask = (z_vals >= z_lo) & (z_vals < z_hi)
        nl = bin_mask & (split_vals < split_median)
        bl = bin_mask & (split_vals >= split_median)
        
        if np.sum(nl) < 50 or np.sum(bl) < 50:
            continue
        
        all_v = values[bin_mask]
        nl_v = values[nl]
        bl_v = values[bl]
        
        var_total = np.var(all_v)
        if var_total < 1e-10:
            continue
        
        n_nl, n_bl = len(nl_v), len(bl_v)
        n_tot = n_nl + n_bl
        mu_all = np.mean(all_v)
        mu_nl, mu_bl = np.mean(nl_v), np.mean(bl_v)
        
        var_between = (n_nl/n_tot)*(mu_nl-mu_all)**2 + (n_bl/n_tot)*(mu_bl-mu_all)**2
        B = var_between / var_total
        
        z_list.append(np.mean(z_vals[bin_mask]))
        B_list.append(B)
    
    if len(z_list) < 3:
        return None
    
    rho, p = spearmanr(z_list, B_list)
    return {'rho_B': rho, 'p_B': p, 'n_bins': len(z_list), 'B_first': B_list[0], 'B_last': B_list[-1]}


# ============================================================
# TEST 1: EW RATIO CANCELLATION
# ============================================================
print(f"\n{'=' * 100}")
print("TEST 1: EW RATIO CANCELLATION — Does B(z) survive in line ratios?")
print("If continuum anisotropy drives divergence, ratios should KILL it")
print("=" * 100)

# We need lines visible at the same redshifts
# MgII (z=0.4-2.5) overlaps with CIII (z=0.8-3.0) at z=0.8-2.5
# Hβ (z=0-1) overlaps with MgII at z=0.4-1.0 (limited)

mgii = d['MGII_BR']
ciii = d['CIII_BR']
civ = d['CIV']
hbeta = d['HBETA_BR']
siiv = d['SIIV_OIV']

# Ratio 1: MgII / CIII (z = 0.8-2.5, both visible)
print(f"\n--- Ratio: log(EW_MgII / EW_CIII) at z=0.8-2.5 ---")

mask_mc = ((z_all >= 0.8) & (z_all < 2.5) &
           (mgii[:, 2] > 0) & (mgii[:, 4] > 100) & (mgii[:, 4] < 30000) &
           (ciii[:, 2] > 0) & (ciii[:, 4] > 100) & (ciii[:, 4] < 30000) &
           np.isfinite(mgii[:, 2]) & np.isfinite(ciii[:, 2]) &
           (lbol > 40) & (lbol < 50))

log_ratio_mc = np.log10(mgii[mask_mc, 2] / ciii[mask_mc, 2])
z_mc = z_all[mask_mc]
fwhm_mgii_mc = mgii[mask_mc, 4]
fwhm_med_mc = np.median(fwhm_mgii_mc)

# Residualize ratio against L_bol
lbol_mc = lbol[mask_mc]
coef = P.polyfit(lbol_mc, log_ratio_mc, 1)
log_ratio_mc_resid = log_ratio_mc - P.polyval(lbol_mc, coef)

# Raw EW for comparison
log_ew_mgii = np.log10(mgii[mask_mc, 2])
coef_ew = P.polyfit(lbol_mc, log_ew_mgii, 1)
log_ew_mgii_resid = log_ew_mgii - P.polyval(lbol_mc, coef_ew)

print(f"  N = {np.sum(mask_mc)}")

r_raw = compute_B_trend_generic(log_ew_mgii_resid, z_mc, fwhm_mgii_mc, fwhm_med_mc)
r_ratio = compute_B_trend_generic(log_ratio_mc_resid, z_mc, fwhm_mgii_mc, fwhm_med_mc)

if r_raw and r_ratio:
    print(f"  Raw MgII EW:       ρ_B = {r_raw['rho_B']:+.3f}, p = {r_raw['p_B']:.4f}")
    print(f"  Ratio MgII/CIII:   ρ_B = {r_ratio['rho_B']:+.3f}, p = {r_ratio['p_B']:.4f}")
    if abs(r_ratio['rho_B']) < abs(r_raw['rho_B']) * 0.5:
        print(f"  🔥 RATIO COLLAPSES B(z) — continuum anisotropy IS the driver")
    elif abs(r_ratio['rho_B']) > abs(r_raw['rho_B']) * 0.75:
        print(f"  ✓ RATIO PRESERVES B(z) — real line physics, NOT orientation")
    else:
        print(f"  ⚠️  Partial reduction")

# Ratio 2: CIV / CIII (z = 1.5-3.0, both visible)
print(f"\n--- Ratio: log(EW_CIV / EW_CIII) at z=1.5-3.0 ---")

mask_cc = ((z_all >= 1.5) & (z_all < 3.0) &
           (civ[:, 2] > 0) & (civ[:, 4] > 100) & (civ[:, 4] < 30000) &
           (ciii[:, 2] > 0) & (ciii[:, 4] > 100) & (ciii[:, 4] < 30000) &
           np.isfinite(civ[:, 2]) & np.isfinite(ciii[:, 2]) &
           (lbol > 40) & (lbol < 50))

log_ratio_cc = np.log10(civ[mask_cc, 2] / ciii[mask_cc, 2])
z_cc = z_all[mask_cc]
fwhm_civ_cc = civ[mask_cc, 4]
fwhm_med_cc = np.median(fwhm_civ_cc)

lbol_cc = lbol[mask_cc]
coef = P.polyfit(lbol_cc, log_ratio_cc, 1)
log_ratio_cc_resid = log_ratio_cc - P.polyval(lbol_cc, coef)

log_ew_civ = np.log10(civ[mask_cc, 2])
coef_ew = P.polyfit(lbol_cc, log_ew_civ, 1)
log_ew_civ_resid = log_ew_civ - P.polyval(lbol_cc, coef_ew)

print(f"  N = {np.sum(mask_cc)}")

r_raw_c = compute_B_trend_generic(log_ew_civ_resid, z_cc, fwhm_civ_cc, fwhm_med_cc)
r_ratio_c = compute_B_trend_generic(log_ratio_cc_resid, z_cc, fwhm_civ_cc, fwhm_med_cc)

if r_raw_c and r_ratio_c:
    print(f"  Raw CIV EW:        ρ_B = {r_raw_c['rho_B']:+.3f}, p = {r_raw_c['p_B']:.4f}")
    print(f"  Ratio CIV/CIII:    ρ_B = {r_ratio_c['rho_B']:+.3f}, p = {r_ratio_c['p_B']:.4f}")
    if abs(r_ratio_c['rho_B']) < abs(r_raw_c['rho_B']) * 0.5:
        print(f"  🔥 RATIO COLLAPSES B(z) — continuum anisotropy IS the driver")
    elif abs(r_ratio_c['rho_B']) > abs(r_raw_c['rho_B']) * 0.75:
        print(f"  ✓ RATIO PRESERVES B(z) — real line physics, NOT orientation")
    else:
        print(f"  ⚠️  Partial reduction")

# Ratio 3: CIV / SiIV (z = 1.5-4.0, both high-IP)
print(f"\n--- Ratio: log(EW_CIV / EW_SiIV) at z=1.5-4.0 ---")

mask_cs = ((z_all >= 1.5) & (z_all < 4.0) &
           (civ[:, 2] > 0) & (civ[:, 4] > 100) & (civ[:, 4] < 30000) &
           (siiv[:, 2] > 0) & (siiv[:, 4] > 100) & (siiv[:, 4] < 30000) &
           np.isfinite(civ[:, 2]) & np.isfinite(siiv[:, 2]) &
           (lbol > 40) & (lbol < 50))

log_ratio_cs = np.log10(civ[mask_cs, 2] / siiv[mask_cs, 2])
z_cs = z_all[mask_cs]
fwhm_civ_cs = civ[mask_cs, 4]
fwhm_med_cs = np.median(fwhm_civ_cs)

lbol_cs = lbol[mask_cs]
coef = P.polyfit(lbol_cs, log_ratio_cs, 1)
log_ratio_cs_resid = log_ratio_cs - P.polyval(lbol_cs, coef)

print(f"  N = {np.sum(mask_cs)}")

r_ratio_cs = compute_B_trend_generic(log_ratio_cs_resid, z_cs, fwhm_civ_cs, fwhm_med_cs)
if r_ratio_cs:
    print(f"  Ratio CIV/SiIV:    ρ_B = {r_ratio_cs['rho_B']:+.3f}, p = {r_ratio_cs['p_B']:.4f}")

# Ratio 4: MgII / SiIV (cross-IP ratio, z=1.5-2.5)
print(f"\n--- Ratio: log(EW_MgII / EW_SiIV) at z=1.5-2.5 (cross-IP) ---")

mask_ms = ((z_all >= 1.5) & (z_all < 2.5) &
           (mgii[:, 2] > 0) & (mgii[:, 4] > 100) & (mgii[:, 4] < 30000) &
           (siiv[:, 2] > 0) & (siiv[:, 4] > 100) & (siiv[:, 4] < 30000) &
           np.isfinite(mgii[:, 2]) & np.isfinite(siiv[:, 2]) &
           (lbol > 40) & (lbol < 50))

log_ratio_ms = np.log10(mgii[mask_ms, 2] / siiv[mask_ms, 2])
z_ms = z_all[mask_ms]
fwhm_mgii_ms = mgii[mask_ms, 4]
fwhm_med_ms = np.median(fwhm_mgii_ms)

lbol_ms = lbol[mask_ms]
coef = P.polyfit(lbol_ms, log_ratio_ms, 1)
log_ratio_ms_resid = log_ratio_ms - P.polyval(lbol_ms, coef)

print(f"  N = {np.sum(mask_ms)}")

r_ratio_ms = compute_B_trend_generic(log_ratio_ms_resid, z_ms, fwhm_mgii_ms, fwhm_med_ms)
if r_ratio_ms:
    print(f"  Ratio MgII/SiIV:   ρ_B = {r_ratio_ms['rho_B']:+.3f}, p = {r_ratio_ms['p_B']:.4f}")


# ============================================================
# TEST 2: SIGMA SPLIT vs FWHM SPLIT
# ============================================================
print(f"\n\n{'=' * 100}")
print("TEST 2: SIGMA SPLIT vs FWHM SPLIT")
print("Line dispersion (sigma, col 3) is less inclination-dependent than FWHM (col 4)")
print("If B(z) disappears with sigma split → orientation was driving it")
print("=" * 100)

for line_name, col_name, z_min, z_max in [
    ('MgII_BR', 'MGII_BR', 0.4, 2.5),
    ('CIV', 'CIV', 1.5, 4.0),
    ('Hβ_BR', 'HBETA_BR', 0.3, 1.0),
]:
    ld = d[col_name]
    
    mask = ((z_all >= z_min) & (z_all < z_max) &
            (ld[:, 4] > 100) & (ld[:, 4] < 30000) &
            (ld[:, 2] > 0) & (ld[:, 2] < 5000) &
            (ld[:, 3] > 10) & (ld[:, 3] < 100) &  # sigma (log scale)
            np.isfinite(ld[:, 2]) & np.isfinite(ld[:, 3]) & np.isfinite(ld[:, 4]) &
            (lbol > 40) & (lbol < 50))
    
    N = np.sum(mask)
    if N < 2000:
        print(f"\n  {line_name}: insufficient data (N={N})")
        continue
    
    z_v = z_all[mask]
    log_ew = np.log10(ld[mask, 2])
    fwhm_v = ld[mask, 4]
    sigma_v = ld[mask, 3]  # log(sigma_km/s)
    
    # Residualize against L_bol
    lbol_v = lbol[mask]
    coef = P.polyfit(lbol_v, log_ew, 1)
    log_ew_resid = log_ew - P.polyval(lbol_v, coef)
    
    fwhm_med = np.median(fwhm_v)
    sigma_med = np.median(sigma_v)
    
    r_fwhm = compute_B_trend_generic(log_ew_resid, z_v, fwhm_v, fwhm_med)
    r_sigma = compute_B_trend_generic(log_ew_resid, z_v, sigma_v, sigma_med)
    
    print(f"\n  {line_name} (N={N}):")
    if r_fwhm:
        print(f"    FWHM split:  ρ_B = {r_fwhm['rho_B']:+.3f}, p = {r_fwhm['p_B']:.4f}")
    if r_sigma:
        print(f"    Sigma split:  ρ_B = {r_sigma['rho_B']:+.3f}, p = {r_sigma['p_B']:.4f}")
    
    if r_fwhm and r_sigma:
        if abs(r_sigma['rho_B']) < abs(r_fwhm['rho_B']) * 0.5:
            print(f"    🔥 Sigma split KILLS B(z) — orientation WAS driving it")
        elif abs(r_sigma['rho_B']) > abs(r_fwhm['rho_B']) * 0.75:
            print(f"    ✓ Sigma split PRESERVES B(z) — NOT orientation")
        else:
            print(f"    ⚠️  Partial reduction")


# ============================================================
# VERDICT
# ============================================================
print(f"\n\n{'=' * 100}")
print("FINAL VERDICT: ORIENTATION vs COSMOLOGY")
print("=" * 100)
print("""
If ratios collapse B(z) AND sigma kills it → Orientation + continuum anisotropy
If ratios preserve B(z) AND sigma preserves it → Real epoch physics / cosmological
If mixed → Both contribute (orientation modulates, epoch underlies)
""")
