#!/usr/bin/env python3
"""
ROGUE-WAVE vs SMOOTH OPERATOR TEST
====================================
Tests whether correlation degradation is:
  A) Cumulative/smooth (additive with total void path length)
  B) Rogue-hit (dominated by single deepest void segment)

Uses Planck CMB lensing κ as matter density proxy along each quasar sightline.
κ at different smoothing scales gives us void structure information.

Key idea: 
- κ_large_scale = total integrated matter (proxy for L_void when inverted)
- κ_small_scale - κ_large_scale = local structure ("roughness")
- min(κ across scales) = deepest void proxy (D_max)
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import stats
from scipy.optimize import curve_fit
import healpy as hp
import json
import os
import warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_rogue_wave'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("ROGUE-WAVE vs SMOOTH OPERATOR TEST")
print("=" * 70)

# ============================================================
# 1. Load DR16Q quasar catalog
# ============================================================
print("\n[1] Loading DR16Q catalog...")
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

z = data['Z_DR16Q']
ra = data['RA']
dec = data['DEC']
snr = data['SN_MEDIAN_ALL']

# Line measurements — DR16Q_prop has 6 columns per line:
# [wave_peak, wave_center, EW, sigma, flux_peak, flux_total]
# Use col 4 (flux_peak) as primary, col 2 (EW) as secondary
# Use col 5 for error proxy (flux_total scatter)

def get_line(name):
    """Extract flux and crude error from DR16Q_prop 6-column line format"""
    raw = data[name]
    flux = raw[:, 4]   # peak flux
    ew = raw[:, 2]     # equivalent width
    # Error: use the _ERR column col 4 if available
    err_name = name + '_ERR'
    if err_name in [c.name for c in hdu[1].columns]:
        err_raw = data[err_name]
        flux_err = err_raw[:, 4] if err_raw.ndim > 1 else err_raw
    else:
        flux_err = np.abs(flux) * 0.1  # fallback 10%
    return flux, ew, flux_err

oiii, oiii_ew, oiii_err = get_line('OIII5007')
hbeta, hbeta_ew, hbeta_err = get_line('HBETA')
oii, oii_ew, oii_err = get_line('OII3728')
sii, sii_ew, sii_err = get_line('SII6718')
nii, nii_ew, nii_err = get_line('NII6585')
halpha, halpha_ew, halpha_err = get_line('HALPHA')

print(f"  Total objects: {len(z)}")

# ============================================================
# 2. Load Planck lensing convergence map
# ============================================================
print("\n[2] Loading Planck lensing κ map...")

# Read alm and reconstruct at multiple smoothing scales
alm_file = '/root/clawd/COM_Lensing_4096_R3.00/MV/dat_klm.fits'
mask_file = '/root/clawd/COM_Lensing_4096_R3.00/mask.fits'

# Read the convergence alm
klm = hp.read_alm(alm_file)
mask = hp.read_map(mask_file, dtype=float)

# Reconstruct κ at multiple smoothing scales
NSIDE = 1024  # ~3.4 arcmin resolution
smoothing_scales = [10, 30, 60, 120, 240]  # arcmin

print(f"  Reconstructing κ maps at {len(smoothing_scales)} smoothing scales...")
kappa_maps = {}
for fwhm_arcmin in smoothing_scales:
    fwhm_rad = np.radians(fwhm_arcmin / 60.0)
    # Apply smoothing in alm space and reconstruct
    klm_smooth = hp.almxfl(klm, hp.gauss_beam(fwhm_rad, lmax=hp.Alm.getlmax(len(klm))))
    kmap = hp.alm2map(klm_smooth, NSIDE, verbose=False)
    kappa_maps[fwhm_arcmin] = kmap
    print(f"    Scale {fwhm_arcmin}': σ(κ) = {np.std(kmap):.6f}")

# Downgrade mask to NSIDE
mask_ds = hp.ud_grade(mask, NSIDE)

# ============================================================
# 3. Quality cuts + compute diagnostic ratios
# ============================================================
print("\n[3] Applying quality cuts...")

# For the rogue-wave test, we need objects with:
# - Good redshift (z > 0.3 for enough path)
# - Good SNR
# - Valid line measurements for at least one diagnostic ratio
# - Not in masked Planck region

# Convert RA/DEC to healpix pixels
theta = np.radians(90.0 - dec)
phi = np.radians(ra)
pixels = hp.ang2pix(NSIDE, theta, phi)

# Mask cut
in_mask = mask_ds[pixels] > 0.5

# Basic quality
good_z = (z > 0.3) & (z < 2.5)
good_snr = snr > 3

# For diagnostic damage: use OIII/Hbeta ratio (well-measured, mid-ladder)
# and SII (if available, high-sensitivity end)
valid_oiii_hb = (oiii > 0) & (hbeta > 0) & (oiii_err > 0) & (hbeta_err > 0)
# Flatten err arrays if multi-dim
if oiii_err.ndim > 1:
    oiii_err_1d = oiii_err[:, 4] if oiii_err.shape[1] > 4 else oiii_err[:, 0]
else:
    oiii_err_1d = oiii_err
if hbeta_err.ndim > 1:
    hbeta_err_1d = hbeta_err[:, 4] if hbeta_err.shape[1] > 4 else hbeta_err[:, 0]
else:
    hbeta_err_1d = hbeta_err
valid_oiii_hb &= (oiii_err_1d > 0) & (hbeta_err_1d > 0)
valid_oiii_hb &= (oiii / np.maximum(oiii_err_1d, 1e-30) > 3) & (hbeta / np.maximum(hbeta_err_1d, 1e-30) > 3)

# Also try OII if available
valid_oii = (oii > 0) & (oii_err.flatten()[:len(oii)] > 0 if oii_err.ndim > 1 else oii_err > 0)

# NII/Halpha for locked control
if nii_err.ndim > 1:
    nii_err_1d = nii_err[:, 4] if nii_err.shape[1] > 4 else nii_err[:, 0]
else:
    nii_err_1d = nii_err
if halpha_err.ndim > 1:
    halpha_err_1d = halpha_err[:, 4] if halpha_err.shape[1] > 4 else halpha_err[:, 0]
else:
    halpha_err_1d = halpha_err
valid_nii_ha = (nii > 0) & (halpha > 0) & (nii_err_1d > 0) & (halpha_err_1d > 0)
valid_nii_ha &= (nii / np.maximum(nii_err_1d, 1e-30) > 3) & (halpha / np.maximum(halpha_err_1d, 1e-30) > 3)

# Primary sample: OIII/Hbeta diagnostic
sel_primary = good_z & good_snr & in_mask & valid_oiii_hb
print(f"  Primary sample (OIII/Hβ): {np.sum(sel_primary)}")

# Locked control: NII/Halpha
sel_locked = good_z & good_snr & in_mask & valid_nii_ha
print(f"  Locked control (NII/Hα): {np.sum(sel_locked)}")

# SII sample (fewer objects, low-z only)
valid_sii = (sii > 0) & (sii_err > 0) & (sii / sii_err > 3)
sel_sii = good_z & good_snr & in_mask & valid_sii & valid_oiii_hb
print(f"  SII sample (SII + OIII/Hβ): {np.sum(sel_sii)}")

# ============================================================
# 4. Compute per-object κ features along each sightline
# ============================================================
print("\n[4] Computing sightline κ features...")

def compute_sightline_features(pix_indices):
    """For each sightline, compute:
    - kappa_mean: average κ across scales (total integrated matter)
    - kappa_min: minimum κ across scales (deepest void proxy = D_max)
    - kappa_range: max - min κ across scales (roughness)
    - kappa_fine: κ at finest scale (10')
    - kappa_coarse: κ at coarsest scale (240')
    - L_void_proxy: sum of (1 - κ/κ_max) across scales where κ < median
    """
    n = len(pix_indices)
    features = np.zeros((n, len(smoothing_scales)))
    
    for i, scale in enumerate(smoothing_scales):
        features[:, i] = kappa_maps[scale][pix_indices]
    
    kappa_mean = np.mean(features, axis=1)
    kappa_min = np.min(features, axis=1)
    kappa_max = np.max(features, axis=1)
    kappa_range = kappa_max - kappa_min  # roughness
    kappa_fine = features[:, 0]   # 10' = most local
    kappa_coarse = features[:, -1]  # 240' = most integrated
    
    # L_void proxy: how much of the multiscale profile is underdense
    median_per_scale = np.median(features, axis=0)
    void_mask = features < median_per_scale[np.newaxis, :]
    L_void = np.sum(void_mask.astype(float) * np.abs(features - median_per_scale[np.newaxis, :]), axis=1)
    
    # D_max proxy: the single deepest underdensity relative to median
    D_max = np.min(features - median_per_scale[np.newaxis, :], axis=1)
    
    return {
        'kappa_mean': kappa_mean,
        'kappa_min': kappa_min,
        'kappa_range': kappa_range,
        'kappa_fine': kappa_fine,
        'kappa_coarse': kappa_coarse,
        'L_void': L_void,
        'D_max': D_max,  # most negative = deepest void
    }

pix_primary = pixels[sel_primary]
feat_primary = compute_sightline_features(pix_primary)

pix_locked = pixels[sel_locked]
feat_locked = compute_sightline_features(pix_locked)

print(f"  κ features computed for {len(pix_primary)} primary + {len(pix_locked)} locked sightlines")

# ============================================================
# 5. Compute "damage score" per object
# ============================================================
print("\n[5] Computing damage scores...")

# Damage = deviation of diagnostic ratio from the local z-bin trend
# In a healthy universe, OIII/Hbeta should have a predictable z-dependence
# "Damage" = how much an object's ratio deviates after removing z-trend

z_primary = z[sel_primary]
ratio_primary = np.log10(oiii[sel_primary] / hbeta[sel_primary])

# Remove z-trend using running median
z_bins_edge = np.percentile(z_primary, np.linspace(0, 100, 21))
z_bins_center = 0.5 * (z_bins_edge[:-1] + z_bins_edge[1:])
z_bin_idx = np.digitize(z_primary, z_bins_edge) - 1
z_bin_idx = np.clip(z_bin_idx, 0, 19)

# Per-bin median and MAD
ratio_residual = np.zeros_like(ratio_primary)
ratio_scatter = np.zeros_like(ratio_primary)
for b in range(20):
    in_bin = z_bin_idx == b
    if np.sum(in_bin) < 10:
        ratio_residual[in_bin] = np.nan
        continue
    med = np.median(ratio_primary[in_bin])
    mad = 1.4826 * np.median(np.abs(ratio_primary[in_bin] - med))
    if mad == 0:
        mad = np.std(ratio_primary[in_bin])
    ratio_residual[in_bin] = (ratio_primary[in_bin] - med) / max(mad, 1e-6)
    ratio_scatter[in_bin] = mad

# Damage = |residual| (how far from local trend)
damage = np.abs(ratio_residual)

# Also compute for locked ratio (NII/Halpha) as control
z_locked = z[sel_locked]
ratio_locked = np.log10(nii[sel_locked] / halpha[sel_locked])

z_bins_locked = np.digitize(z_locked, np.percentile(z_locked, np.linspace(0, 100, 21))) - 1
z_bins_locked = np.clip(z_bins_locked, 0, 19)
ratio_resid_locked = np.zeros_like(ratio_locked)
for b in range(20):
    in_bin = z_bins_locked == b
    if np.sum(in_bin) < 10:
        ratio_resid_locked[in_bin] = np.nan
        continue
    med = np.median(ratio_locked[in_bin])
    mad = 1.4826 * np.median(np.abs(ratio_locked[in_bin] - med))
    if mad == 0:
        mad = np.std(ratio_locked[in_bin])
    ratio_resid_locked[in_bin] = (ratio_locked[in_bin] - med) / max(mad, 1e-6)

damage_locked = np.abs(ratio_resid_locked)

print(f"  Primary damage: mean={np.nanmean(damage):.3f}, std={np.nanstd(damage):.3f}")
print(f"  Locked damage:  mean={np.nanmean(damage_locked):.3f}, std={np.nanstd(damage_locked):.3f}")

# ============================================================
# 6. THE CORE TEST: Smooth vs Rogue-Hit
# ============================================================
print("\n" + "=" * 70)
print("CORE TEST: SMOOTH vs ROGUE-HIT")
print("=" * 70)

# Remove NaN
valid = np.isfinite(damage) & np.isfinite(feat_primary['L_void']) & np.isfinite(feat_primary['D_max'])
dam = damage[valid]
L_void = feat_primary['L_void'][valid]
D_max = feat_primary['D_max'][valid]  # most negative = deepest void
kappa_mean = feat_primary['kappa_mean'][valid]
kappa_range = feat_primary['kappa_range'][valid]
z_valid = z_primary[valid]

print(f"\n  Valid objects: {np.sum(valid)}")

# --- Model 1: Smooth (damage ~ a * L_void + b * z) ---
print("\n--- Model 1: SMOOTH (damage ~ L_void + z) ---")
from sklearn.linear_model import LinearRegression
X_smooth = np.column_stack([L_void, z_valid])
reg_smooth = LinearRegression().fit(X_smooth, dam)
pred_smooth = reg_smooth.predict(X_smooth)
r2_smooth = 1 - np.sum((dam - pred_smooth)**2) / np.sum((dam - np.mean(dam))**2)
rho_smooth, p_smooth = stats.spearmanr(L_void, dam)
print(f"  L_void coefficient: {reg_smooth.coef_[0]:.6f}")
print(f"  z coefficient: {reg_smooth.coef_[1]:.6f}")
print(f"  R² (smooth model): {r2_smooth:.6f}")
print(f"  Spearman(L_void, damage): ρ = {rho_smooth:.4f}, p = {p_smooth:.2e}")

# --- Model 2: Rogue-Hit (damage ~ b * |D_max| + c * z) ---
print("\n--- Model 2: ROGUE-HIT (damage ~ |D_max| + z) ---")
D_max_abs = np.abs(D_max)  # deeper void = larger value
X_rogue = np.column_stack([D_max_abs, z_valid])
reg_rogue = LinearRegression().fit(X_rogue, dam)
pred_rogue = reg_rogue.predict(X_rogue)
r2_rogue = 1 - np.sum((dam - pred_rogue)**2) / np.sum((dam - np.mean(dam))**2)
rho_rogue, p_rogue = stats.spearmanr(D_max_abs, dam)
print(f"  |D_max| coefficient: {reg_rogue.coef_[0]:.6f}")
print(f"  z coefficient: {reg_rogue.coef_[1]:.6f}")
print(f"  R² (rogue model): {r2_rogue:.6f}")
print(f"  Spearman(|D_max|, damage): ρ = {rho_rogue:.4f}, p = {p_rogue:.2e}")

# --- Model 3: Combined (both features) ---
print("\n--- Model 3: COMBINED (L_void + |D_max| + z) ---")
X_combined = np.column_stack([L_void, D_max_abs, z_valid])
reg_combined = LinearRegression().fit(X_combined, dam)
pred_combined = reg_combined.predict(X_combined)
r2_combined = 1 - np.sum((dam - pred_combined)**2) / np.sum((dam - np.mean(dam))**2)
print(f"  L_void coefficient: {reg_combined.coef_[0]:.6f}")
print(f"  |D_max| coefficient: {reg_combined.coef_[1]:.6f}")
print(f"  z coefficient: {reg_combined.coef_[2]:.6f}")
print(f"  R² (combined): {r2_combined:.6f}")

# --- Model comparison ---
print("\n--- MODEL COMPARISON ---")
n = len(dam)
k_smooth = 2
k_rogue = 2
k_combined = 3
aic_smooth = n * np.log(np.sum((dam - pred_smooth)**2)/n) + 2*k_smooth
aic_rogue = n * np.log(np.sum((dam - pred_rogue)**2)/n) + 2*k_rogue
aic_combined = n * np.log(np.sum((dam - pred_combined)**2)/n) + 2*k_combined
print(f"  AIC Smooth:   {aic_smooth:.1f}")
print(f"  AIC Rogue:    {aic_rogue:.1f}")
print(f"  AIC Combined: {aic_combined:.1f}")
print(f"  ΔAIC (Rogue - Smooth): {aic_rogue - aic_smooth:.1f}")

if aic_rogue < aic_smooth - 10:
    winner = "ROGUE-HIT"
elif aic_smooth < aic_rogue - 10:
    winner = "SMOOTH"
else:
    winner = "INDISTINGUISHABLE"
print(f"  >>> Winner: {winner}")

# ============================================================
# 7. DISTRIBUTION SHAPE TEST: Unimodal vs Mixture
# ============================================================
print("\n" + "=" * 70)
print("DISTRIBUTION SHAPE TEST")
print("=" * 70)

# Test whether damage distributions become bimodal/heavy-tailed at high-z
# Split into z-bins and test distribution shape

z_terciles = np.percentile(z_valid, [0, 33, 67, 100])
labels = ['Low-z', 'Mid-z', 'High-z']

from scipy.stats import kurtosis, skew, normaltest, anderson

print(f"\n  {'Bin':<10} {'N':>6} {'Mean':>8} {'Std':>8} {'Skew':>8} {'Kurt':>8} {'Normal_p':>10} {'Anderson':>10}")
print("  " + "-" * 80)

shape_results = []
for i in range(3):
    in_bin = (z_valid >= z_terciles[i]) & (z_valid < z_terciles[i+1])
    d = dam[in_bin]
    if len(d) < 20:
        continue
    sk = skew(d)
    ku = kurtosis(d)  # excess kurtosis
    _, norm_p = normaltest(d)
    ad = anderson(d, 'norm')
    
    shape_results.append({
        'bin': labels[i],
        'n': len(d),
        'mean': float(np.mean(d)),
        'std': float(np.std(d)),
        'skew': float(sk),
        'kurtosis': float(ku),
        'normal_p': float(norm_p),
        'anderson_stat': float(ad.statistic),
    })
    
    print(f"  {labels[i]:<10} {len(d):>6} {np.mean(d):>8.3f} {np.std(d):>8.3f} {sk:>8.3f} {ku:>8.3f} {norm_p:>10.2e} {ad.statistic:>10.2f}")

# Kurtosis trend
if len(shape_results) == 3:
    kurt_trend = [s['kurtosis'] for s in shape_results]
    print(f"\n  Kurtosis trend (low→high z): {kurt_trend[0]:.3f} → {kurt_trend[1]:.3f} → {kurt_trend[2]:.3f}")
    if kurt_trend[2] > kurt_trend[0] + 0.5:
        print("  >>> HEAVIER TAILS at high-z — consistent with ROGUE-HIT / intermittent events")
    elif kurt_trend[2] < kurt_trend[0] - 0.5:
        print("  >>> LIGHTER TAILS at high-z — consistent with SMOOTH averaging")
    else:
        print("  >>> Kurtosis roughly stable — inconclusive on rogue vs smooth")

# ============================================================
# 8. QUINTILE DEEP-DIVE: Damage vs void depth
# ============================================================
print("\n" + "=" * 70)
print("QUINTILE ANALYSIS: DAMAGE vs VOID DEPTH")
print("=" * 70)

# D_max quintiles
dmax_quintiles = np.percentile(D_max, [0, 20, 40, 60, 80, 100])
print(f"\n  D_max quintile edges: {[f'{v:.4f}' for v in dmax_quintiles]}")
print(f"  (More negative = deeper void)")

print(f"\n  {'Quintile':<12} {'D_max range':>20} {'N':>6} {'Mean damage':>12} {'Std':>8}")
print("  " + "-" * 65)

quintile_means = []
for q in range(5):
    in_q = (D_max >= dmax_quintiles[q]) & (D_max < dmax_quintiles[q+1])
    if q == 4:  # include upper edge
        in_q = (D_max >= dmax_quintiles[q]) & (D_max <= dmax_quintiles[q+1])
    d = dam[in_q]
    quintile_means.append(float(np.mean(d)))
    print(f"  Q{q+1} ({'deepest' if q==0 else 'densest' if q==4 else 'mid'}){'':<5} [{dmax_quintiles[q]:>8.4f}, {dmax_quintiles[q+1]:>8.4f}] {len(d):>6} {np.mean(d):>12.4f} {np.std(d):>8.4f}")

# Is the deepest-void quintile the most damaged?
rho_quint, p_quint = stats.spearmanr(range(5), quintile_means)
print(f"\n  Quintile trend (Spearman): ρ = {rho_quint:.3f}, p = {p_quint:.3f}")
if quintile_means[0] > quintile_means[4]:
    print("  >>> DEEPEST VOIDS = MOST DAMAGE — consistent with operator")
else:
    print("  >>> DENSEST regions = MOST DAMAGE — opposite to operator prediction")

# ============================================================
# 9. LOCKED CONTROL: Same test on NII/Halpha
# ============================================================
print("\n" + "=" * 70)
print("CONTROL: LOCKED RATIO (NII/Hα) — SHOULD SHOW NO VOID DEPENDENCE")
print("=" * 70)

valid_lock = np.isfinite(damage_locked) & np.isfinite(feat_locked['D_max'])
dam_lock = damage_locked[valid_lock]
D_max_lock = feat_locked['D_max'][valid_lock]
L_void_lock = feat_locked['L_void'][valid_lock]

rho_lock_dmax, p_lock_dmax = stats.spearmanr(np.abs(D_max_lock), dam_lock)
rho_lock_lvoid, p_lock_lvoid = stats.spearmanr(L_void_lock, dam_lock)
print(f"  Locked sample: {np.sum(valid_lock)} objects")
print(f"  Spearman(|D_max|, locked damage): ρ = {rho_lock_dmax:.4f}, p = {p_lock_dmax:.2e}")
print(f"  Spearman(L_void, locked damage):  ρ = {rho_lock_lvoid:.4f}, p = {p_lock_lvoid:.2e}")
if abs(rho_lock_dmax) < 0.01 and p_lock_dmax > 0.05:
    print("  >>> CONTROL PASSED — locked ratios unaffected by void structure")
else:
    print(f"  >>> CONTROL {'MARGINAL' if abs(rho_lock_dmax) < 0.03 else 'FAILED'}")

# ============================================================
# 10. Z-STRATIFIED ROGUE TEST
# ============================================================
print("\n" + "=" * 70)
print("Z-STRATIFIED: DOES ROGUE BEHAVIOR GROW WITH DISTANCE?")
print("=" * 70)

z_edges = np.percentile(z_valid, [0, 25, 50, 75, 100])
print(f"\n  {'z-bin':<15} {'N':>6} {'ρ(Dmax,dam)':>12} {'p':>10} {'ρ(Lvoid,dam)':>13} {'p':>10}")
print("  " + "-" * 70)

stratified_results = []
for i in range(4):
    in_bin = (z_valid >= z_edges[i]) & (z_valid < z_edges[i+1])
    if i == 3:
        in_bin = (z_valid >= z_edges[i]) & (z_valid <= z_edges[i+1])
    if np.sum(in_bin) < 30:
        continue
    rho_d, p_d = stats.spearmanr(np.abs(D_max[in_bin]), dam[in_bin])
    rho_l, p_l = stats.spearmanr(L_void[in_bin], dam[in_bin])
    label = f"[{z_edges[i]:.2f}, {z_edges[i+1]:.2f}]"
    print(f"  {label:<15} {np.sum(in_bin):>6} {rho_d:>12.4f} {p_d:>10.2e} {rho_l:>13.4f} {p_l:>10.2e}")
    stratified_results.append({
        'z_lo': float(z_edges[i]),
        'z_hi': float(z_edges[i+1]),
        'n': int(np.sum(in_bin)),
        'rho_dmax': float(rho_d),
        'p_dmax': float(p_d),
        'rho_lvoid': float(rho_l),
        'p_lvoid': float(p_l),
    })

# Does the effect grow?
if len(stratified_results) >= 3:
    rho_trend = [s['rho_dmax'] for s in stratified_results]
    print(f"\n  ρ(D_max) trend across z: {' → '.join([f'{r:.4f}' for r in rho_trend])}")

# ============================================================
# 11. BIMODALITY TEST (Hartigans' dip test approximation)
# ============================================================
print("\n" + "=" * 70)
print("BIMODALITY TEST (damage distribution)")
print("=" * 70)

# Use kernel density estimation to check for multiple modes
from scipy.signal import argrelextrema

for i, label in enumerate(labels):
    in_bin = (z_valid >= z_terciles[i]) & (z_valid < z_terciles[i+1])
    d = dam[in_bin]
    if len(d) < 50:
        continue
    
    # KDE
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(d, bw_method=0.2)
    x_grid = np.linspace(np.min(d), np.percentile(d, 99), 500)
    density = kde(x_grid)
    
    # Count peaks
    peaks = argrelextrema(density, np.greater, order=20)[0]
    n_peaks = len(peaks)
    
    print(f"  {label}: {n_peaks} mode(s) detected in damage distribution (N={len(d)})")
    if n_peaks > 1:
        peak_locs = x_grid[peaks]
        print(f"    Peak locations: {[f'{p:.2f}' for p in peak_locs]}")
        print(f"    >>> MULTIMODAL — supports MIXTURE / ROGUE-HIT model")

# ============================================================
# SAVE RESULTS
# ============================================================
print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

results = {
    'test': 'Rogue-Wave vs Smooth Operator',
    'n_primary': int(np.sum(sel_primary)),
    'n_locked': int(np.sum(sel_locked)),
    'n_valid': int(np.sum(valid)),
    'models': {
        'smooth': {
            'R2': float(r2_smooth),
            'AIC': float(aic_smooth),
            'rho_Lvoid_damage': float(rho_smooth),
            'p_Lvoid_damage': float(p_smooth),
        },
        'rogue': {
            'R2': float(r2_rogue),
            'AIC': float(aic_rogue),
            'rho_Dmax_damage': float(rho_rogue),
            'p_Dmax_damage': float(p_rogue),
        },
        'combined': {
            'R2': float(r2_combined),
            'AIC': float(aic_combined),
        },
        'winner': winner,
        'delta_AIC': float(aic_rogue - aic_smooth),
    },
    'shape': shape_results,
    'quintiles': {
        'means': quintile_means,
        'rho': float(rho_quint),
        'p': float(p_quint),
        'deepest_most_damaged': bool(quintile_means[0] > quintile_means[4]),
    },
    'locked_control': {
        'rho_dmax': float(rho_lock_dmax),
        'p_dmax': float(p_lock_dmax),
        'rho_lvoid': float(rho_lock_lvoid),
        'p_lvoid': float(p_lock_lvoid),
    },
    'stratified': stratified_results,
}

with open(f'{OUTDIR}/rogue_wave_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n  Results saved to {OUTDIR}/rogue_wave_results.json")

# ============================================================
# VERDICT
# ============================================================
print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)

print(f"""
  Model comparison:
    Smooth (L_void):     R² = {r2_smooth:.6f}, AIC = {aic_smooth:.1f}
    Rogue-Hit (|D_max|): R² = {r2_rogue:.6f}, AIC = {aic_rogue:.1f}
    Combined:            R² = {r2_combined:.6f}, AIC = {aic_combined:.1f}
    
  Winner: {winner} (ΔAIC = {aic_rogue - aic_smooth:.1f})
  
  Distribution shape:
    Kurtosis trend: {' → '.join([f"{s['kurtosis']:.2f}" for s in shape_results])}
    
  Locked control:
    ρ(|D_max|, locked damage) = {rho_lock_dmax:.4f} (p = {p_lock_dmax:.2e})
    
  Quintile monotonicity:
    ρ = {rho_quint:.3f} (p = {p_quint:.3f})
    Deepest voids most damaged: {quintile_means[0] > quintile_means[4]}
""")

print("Done.")
