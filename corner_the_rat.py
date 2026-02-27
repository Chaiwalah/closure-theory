#!/usr/bin/env python3
"""
CORNER THE RAT — Three decisive tests
========================================
1. D1 diagnostic: Is the phase transition real or contamination?
2. Metric swap: Does HB_FWHM predict damage in ratios that DON'T use Hβ?
3. 2D state map: (FWHM, Ledd) heatmap of damage-vs-z slopes
"""

import numpy as np
from astropy.io import fits
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_corner_rat'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("CORNER THE RAT")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

def fix(arr): return arr.astype(arr.dtype.newbyteorder('='))

z = fix(data['Z_DR16Q'])
snr = fix(data['SN_MEDIAN_ALL'])
ebv = fix(data['EBV'])
oiii_flux = fix(data['OIII5007'][:, 4])
oiii_ew = fix(data['OIII5007'][:, 2])
oiii_fwhm = fix(data['OIII5007'][:, 3])
hbeta_flux = fix(data['HBETA'][:, 4])
hbeta_ew = fix(data['HBETA'][:, 2])
hbeta_fwhm = fix(data['HBETA'][:, 3])
hbeta_br_flux = fix(data['HBETA_BR'][:, 4])
hbeta_br_fwhm = fix(data['HBETA_BR'][:, 3])
nii_flux = fix(data['NII6585'][:, 4])
halpha_flux = fix(data['HALPHA'][:, 4])
halpha_fwhm = fix(data['HALPHA'][:, 3])
sii_flux = fix(data['SII6718'][:, 4])
oii_flux = fix(data['OII3728'][:, 4])
logLbol = fix(data['LOGLBOL'])
logMbh = fix(data['LOGMBH'])
logLedd = fix(data['LOGLEDD_RATIO'])

# Fit quality stats
hb_stat = fix(data['HB_COMP_STAT'])
ha_stat = fix(data['HA_COMP_STAT'])

# FeII
feii_ew = fix(data['FEII_OPT_EW'])

N_total = len(z)

# Standard quality cuts
good = ((z > 0.3) & (z < 2.0) & (snr > 3) &
        (oiii_flux > 0) & (hbeta_flux > 0) &
        np.isfinite(oiii_flux) & np.isfinite(hbeta_flux) &
        np.isfinite(logLbol) & np.isfinite(logMbh) & np.isfinite(logLedd) &
        (logLbol > 0) & (logMbh > 0) &
        np.isfinite(hbeta_fwhm) & (hbeta_fwhm > 0))

idx = np.where(good)[0]
N = len(idx)
Z = z[idx]
print(f"  Base sample: {N}")

# Log ratio + simple z-detrended damage (no GBR — keep it fast and transparent)
log_ratio = np.log10(oiii_flux[idx] / hbeta_flux[idx])

def compute_damage(ratio, redshift, n_bins=20):
    """Z-bin normalized absolute residual"""
    edges = np.percentile(redshift, np.linspace(0, 100, n_bins + 1))
    zb = np.clip(np.digitize(redshift, edges) - 1, 0, n_bins - 1)
    damage = np.full(len(ratio), np.nan)
    for b in range(n_bins):
        m = zb == b
        if np.sum(m) < 10: continue
        med = np.median(ratio[m])
        mad = max(1.4826 * np.median(np.abs(ratio[m] - med)), 1e-6)
        damage[m] = np.abs((ratio[m] - med) / mad)
    return damage

damage = compute_damage(log_ratio, Z)
HB_FWHM = np.log10(np.maximum(hbeta_fwhm[idx], 1))

# Decile edges for HB_FWHM
fwhm_edges = np.percentile(HB_FWHM, np.linspace(0, 100, 11))
fwhm_bin = np.clip(np.digitize(HB_FWHM, fwhm_edges) - 1, 0, 9)

# ============================================================
# TEST 1: D1 DIAGNOSTIC — What is D1?
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: D1 DIAGNOSTIC — Is the phase transition real?")
print("=" * 70)

d1 = fwhm_bin == 0  # narrowest 10%
d2_10 = fwhm_bin > 0  # rest

print(f"\n  D1 (narrowest 10%): {np.sum(d1)} objects")
print(f"  D2-D10: {np.sum(d2_10)} objects")

# A) What does D1 look like?
print(f"\n  --- D1 Properties ---")
print(f"  HB_FWHM range: [{np.min(hbeta_fwhm[idx][d1]):.1f}, {np.max(hbeta_fwhm[idx][d1]):.1f}] km/s")
print(f"  Median HB_FWHM: {np.median(hbeta_fwhm[idx][d1]):.1f} km/s")
print(f"  log(OIII/Hβ) median: {np.median(log_ratio[d1]):.3f}")
print(f"  log(OIII/Hβ) std: {np.std(log_ratio[d1]):.3f}")
print(f"  logLedd median: {np.median(logLedd[idx][d1]):.3f}")
print(f"  logMbh median: {np.median(logMbh[idx][d1]):.3f}")
print(f"  SNR median: {np.median(snr[idx][d1]):.1f}")
print(f"  E(B-V) median: {np.median(ebv[idx][d1]):.4f}")

# B) Broad Hβ detection in D1
hb_br_detected = (hbeta_br_flux[idx] > 0) & np.isfinite(hbeta_br_flux[idx])
d1_br_frac = np.sum(d1 & hb_br_detected) / np.sum(d1)
rest_br_frac = np.sum(d2_10 & hb_br_detected) / np.sum(d2_10)
print(f"\n  Broad Hβ detection rate:")
print(f"    D1: {100*d1_br_frac:.1f}%")
print(f"    D2-D10: {100*rest_br_frac:.1f}%")

# C) OIII dominance (NLR-dominated?)
oiii_dominance = log_ratio  # high = OIII >> Hβ = NLR dominated
print(f"\n  OIII/Hβ > 3 (NLR dominated):")
print(f"    D1: {100*np.mean(log_ratio[d1] > np.log10(3)):.1f}%")
print(f"    D2-D10: {100*np.mean(log_ratio[d2_10] > np.log10(3)):.1f}%")

# D) FeII strength (EV1 marker)
feii_valid = np.isfinite(feii_ew[idx]) & (feii_ew[idx] > 0)
if np.sum(d1 & feii_valid) > 10:
    print(f"\n  FeII EW (EV1 marker):")
    print(f"    D1 median: {np.median(feii_ew[idx][d1 & feii_valid]):.1f}")
    print(f"    D2-D10 median: {np.median(feii_ew[idx][d2_10 & feii_valid]):.1f}")

# E) Fit quality
print(f"\n  Hβ fit stat:")
print(f"    D1 median: {np.median(hb_stat[idx][d1]):.2f}")
print(f"    D2-D10 median: {np.median(hb_stat[idx][d2_10]):.2f}")

# F) NLS1 check: narrow Hβ + high Ledd + strong FeII = NLS1
nls1_like = d1 & (logLedd[idx] > -0.5)
print(f"\n  NLS1-like (D1 + logLedd > -0.5): {np.sum(nls1_like)} ({100*np.mean(nls1_like[d1]):.1f}% of D1)")

# G) REPEAT THE SLOPE TEST WITHOUT D1
print(f"\n  --- Slopes WITHOUT D1 (D2-D10 only) ---")
for d in range(1, 10):
    in_d = fwhm_bin == d
    sl, _, _, p, _ = stats.linregress(Z[in_d], damage[in_d])
    label = "NARROW" if d < 3 else "BROAD" if d > 6 else "mid"
    print(f"  D{d+1}: slope = {sl:+.4f}, p = {p:.2e}, N = {np.sum(in_d)} {label}")

# H) Does D1 slope survive with ONLY broad Hβ detections?
d1_br = d1 & hb_br_detected
if np.sum(d1_br) > 50:
    sl, _, _, p, _ = stats.linregress(Z[d1_br], damage[d1_br])
    print(f"\n  D1 with broad Hβ detected ONLY: slope = {sl:+.4f}, p = {p:.2e}, N = {np.sum(d1_br)}")
else:
    print(f"\n  D1 with broad Hβ: only {np.sum(d1_br)} objects — too few")

# ============================================================
# TEST 2: METRIC SWAP — Does FWHM predict damage in non-Hβ ratios?
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: METRIC SWAP — Cross-ratio validation")
print("=" * 70)

# For each alternative ratio, compute damage and test HB_FWHM interaction
alt_ratios = {}

# NII/Hα (locked-ish)
nii_ha_good = (nii_flux[idx] > 0) & (halpha_flux[idx] > 0) & np.isfinite(nii_flux[idx]) & np.isfinite(halpha_flux[idx])
if np.sum(nii_ha_good) > 500:
    alt_ratios['NII/Hα'] = (np.log10(nii_flux[idx][nii_ha_good] / halpha_flux[idx][nii_ha_good]),
                             Z[nii_ha_good], HB_FWHM[nii_ha_good], logLedd[idx][nii_ha_good])

# SII/OIII (high sensitivity / locked)
sii_oiii_good = (sii_flux[idx] > 0) & (oiii_flux[idx] > 0) & np.isfinite(sii_flux[idx])
if np.sum(sii_oiii_good) > 500:
    alt_ratios['SII/OIII'] = (np.log10(sii_flux[idx][sii_oiii_good] / oiii_flux[idx][sii_oiii_good]),
                               Z[sii_oiii_good], HB_FWHM[sii_oiii_good], logLedd[idx][sii_oiii_good])

# OII/OIII (both forbidden, different sensitivity)
oii_oiii_good = (oii_flux[idx] > 0) & (oiii_flux[idx] > 0) & np.isfinite(oii_flux[idx])
if np.sum(oii_oiii_good) > 500:
    alt_ratios['OII/OIII'] = (np.log10(oii_flux[idx][oii_oiii_good] / oiii_flux[idx][oii_oiii_good]),
                               Z[oii_oiii_good], HB_FWHM[oii_oiii_good], logLedd[idx][oii_oiii_good])

# OIII_EW / HB_EW ratio (EW-based, different from flux)
ew_good = (oiii_ew[idx] > 0) & (hbeta_ew[idx] > 0) & np.isfinite(oiii_ew[idx]) & np.isfinite(hbeta_ew[idx])
if np.sum(ew_good) > 500:
    alt_ratios['OIII_EW/HB_EW'] = (np.log10(oiii_ew[idx][ew_good] / hbeta_ew[idx][ew_good]),
                                     Z[ew_good], HB_FWHM[ew_good], logLedd[idx][ew_good])

print(f"\n  {'Ratio':<20} {'N':>8} {'ρ(FWHM×z, dam)':>16} {'p':>12} {'ρ(Ledd×z, dam)':>16} {'p':>12} {'Interpretation'}")
print("  " + "-" * 95)

for name, (ratio, z_alt, fwhm_alt, ledd_alt) in alt_ratios.items():
    dam_alt = compute_damage(ratio, z_alt)
    valid = np.isfinite(dam_alt)
    
    z_std = (z_alt[valid] - np.mean(z_alt[valid])) / np.std(z_alt[valid])
    fwhm_std = (fwhm_alt[valid] - np.mean(fwhm_alt[valid])) / np.std(fwhm_alt[valid])
    ledd_std = (ledd_alt[valid] - np.mean(ledd_alt[valid])) / np.std(ledd_alt[valid])
    
    rho_f, p_f = stats.spearmanr(fwhm_std * z_std, dam_alt[valid])
    rho_l, p_l = stats.spearmanr(ledd_std * z_std, dam_alt[valid])
    
    if abs(rho_f) > 0.02 and p_f < 0.001:
        interp = "<<< HANDSHAKE TRANSFERS"
    elif p_f < 0.05:
        interp = "marginal"
    else:
        interp = "no transfer"
    
    print(f"  {name:<20} {np.sum(valid):>8} {rho_f:>+16.4f} {p_f:>12.2e} {rho_l:>+16.4f} {p_l:>12.2e} {interp}")

# Also test OIII/Hβ for reference
z_std_ref = (Z - np.mean(Z)) / np.std(Z)
fwhm_std_ref = (HB_FWHM - np.mean(HB_FWHM)) / np.std(HB_FWHM)
ledd_std_ref = (logLedd[idx] - np.mean(logLedd[idx])) / np.std(logLedd[idx])
valid_ref = np.isfinite(damage)
rho_f_ref, p_f_ref = stats.spearmanr(fwhm_std_ref[valid_ref] * z_std_ref[valid_ref], damage[valid_ref])
rho_l_ref, p_l_ref = stats.spearmanr(ledd_std_ref[valid_ref] * z_std_ref[valid_ref], damage[valid_ref])
print(f"  {'OIII/Hβ (ref)':<20} {np.sum(valid_ref):>8} {rho_f_ref:>+16.4f} {p_f_ref:>12.2e} {rho_l_ref:>+16.4f} {p_l_ref:>12.2e} <<< REFERENCE")

# ============================================================
# TEST 3: 2D STATE MAP — (FWHM, Ledd) heatmap
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: 2D STATE MAP — (HB_FWHM, logLedd) slope heatmap")
print("=" * 70)

valid_2d = np.isfinite(damage) & np.isfinite(HB_FWHM) & np.isfinite(logLedd[idx])
Z_2d = Z[valid_2d]
D_2d = damage[valid_2d]
F_2d = HB_FWHM[valid_2d]
L_2d = logLedd[idx][valid_2d]

# 7×7 grid (adaptive bins)
n_grid = 7
f_edges = np.percentile(F_2d, np.linspace(0, 100, n_grid + 1))
l_edges = np.percentile(L_2d, np.linspace(0, 100, n_grid + 1))

f_bin = np.clip(np.digitize(F_2d, f_edges) - 1, 0, n_grid - 1)
l_bin = np.clip(np.digitize(L_2d, l_edges) - 1, 0, n_grid - 1)

print(f"\n  Grid: {n_grid}×{n_grid} = {n_grid**2} cells")
print(f"  FWHM edges: {[f'{e:.2f}' for e in f_edges]}")
print(f"  Ledd edges: {[f'{e:.2f}' for e in l_edges]}")

slope_map = np.full((n_grid, n_grid), np.nan)
p_map = np.full((n_grid, n_grid), np.nan)
n_map = np.zeros((n_grid, n_grid), dtype=int)

print(f"\n  {'':>6}", end='')
for j in range(n_grid):
    l_mid = 0.5*(l_edges[j] + l_edges[j+1])
    print(f"  L={l_mid:+.1f}", end='')
print()
print("  " + "-" * (8 + n_grid * 10))

for i in range(n_grid):
    f_mid = 0.5*(f_edges[i] + f_edges[i+1])
    print(f"  F={f_mid:.1f}", end='')
    for j in range(n_grid):
        cell = (f_bin == i) & (l_bin == j)
        n_cell = np.sum(cell)
        n_map[i, j] = n_cell
        if n_cell < 30:
            print(f"  {'---':>7}", end='')
            continue
        sl, _, _, p, _ = stats.linregress(Z_2d[cell], D_2d[cell])
        slope_map[i, j] = sl
        p_map[i, j] = p
        marker = '*' if p < 0.01 else '.' if p < 0.05 else ' '
        sign = '+' if sl > 0 else '-'
        print(f"  {sign}{abs(sl):.3f}{marker}", end='')
    print()

# Find the sign flip ridge
print(f"\n  --- Sign Flip Analysis ---")
pos_cells = np.sum(slope_map > 0)
neg_cells = np.sum(slope_map < 0)
total_cells = np.sum(np.isfinite(slope_map))
print(f"  Positive slopes: {pos_cells}/{total_cells}")
print(f"  Negative slopes: {neg_cells}/{total_cells}")

# Where is the boundary?
for i in range(n_grid):
    for j in range(n_grid):
        if np.isfinite(slope_map[i, j]) and slope_map[i, j] < -0.3:
            f_mid = 0.5*(f_edges[i] + f_edges[i+1])
            l_mid = 0.5*(l_edges[j] + l_edges[j+1])
            print(f"  Strong negative: FWHM={f_mid:.2f}, Ledd={l_mid:.2f}, slope={slope_map[i,j]:+.3f}, N={n_map[i,j]}")

# ============================================================
# TEST 4: MATCHED STATE COMPARISON (low-z vs high-z at fixed state)
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: MATCHED STATE — Same accretion state, different z")
print("=" * 70)

# Split into low-z and high-z
z_med = np.median(Z)
lo_z = Z < z_med
hi_z = Z >= z_med

# For matching: use HB_FWHM + logLedd + logLbol + SNR
from sklearn.neighbors import NearestNeighbors

match_features = ['HB_FWHM', 'logLedd', 'logLbol', 'SNR']
X_lo = np.column_stack([HB_FWHM[lo_z], logLedd[idx][lo_z], logLbol[idx][lo_z], snr[idx][lo_z]])
X_hi = np.column_stack([HB_FWHM[hi_z], logLedd[idx][hi_z], logLbol[idx][hi_z], snr[idx][hi_z]])

# Drop NaN
valid_lo = np.all(np.isfinite(X_lo), axis=1) & np.isfinite(damage[lo_z])
valid_hi = np.all(np.isfinite(X_hi), axis=1) & np.isfinite(damage[hi_z])

X_lo_c = X_lo[valid_lo]
X_hi_c = X_hi[valid_hi]
D_lo = damage[lo_z][valid_lo]
D_hi = damage[hi_z][valid_hi]

# Standardize
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler().fit(np.vstack([X_lo_c, X_hi_c]))
X_lo_s = scaler.transform(X_lo_c)
X_hi_s = scaler.transform(X_hi_c)

# Find nearest neighbor in high-z for each low-z object
nn = NearestNeighbors(n_neighbors=1, metric='euclidean')
nn.fit(X_hi_s)
distances, indices = nn.kneighbors(X_lo_s)

# Quality threshold: only keep good matches
dist_thresh = np.percentile(distances, 90)  # top 90% closest
good_match = distances.ravel() < dist_thresh

D_lo_matched = D_lo[good_match]
D_hi_matched = D_hi[indices.ravel()[good_match]]

print(f"  Low-z: {np.sum(valid_lo)}, High-z: {np.sum(valid_hi)}")
print(f"  Matched pairs: {np.sum(good_match)}")
print(f"  Match distance (median): {np.median(distances[good_match]):.3f}")

# Compare damage distributions
print(f"\n  Matched damage comparison:")
print(f"    Low-z mean damage:  {np.mean(D_lo_matched):.4f} ± {np.std(D_lo_matched)/np.sqrt(len(D_lo_matched)):.4f}")
print(f"    High-z mean damage: {np.mean(D_hi_matched):.4f} ± {np.std(D_hi_matched)/np.sqrt(len(D_hi_matched)):.4f}")

t_stat, t_p = stats.ttest_rel(D_lo_matched, D_hi_matched)
ks_stat, ks_p = stats.ks_2samp(D_lo_matched, D_hi_matched)
wilcox_stat, wilcox_p = stats.wilcoxon(D_lo_matched - D_hi_matched)

print(f"    Paired t-test: t = {t_stat:.3f}, p = {t_p:.2e}")
print(f"    KS test: D = {ks_stat:.4f}, p = {ks_p:.2e}")
print(f"    Wilcoxon: p = {wilcox_p:.2e}")

diff = np.mean(D_hi_matched) - np.mean(D_lo_matched)
if t_p < 0.01 and diff > 0:
    print(f"    >>> HIGH-Z MORE DAMAGED at matched state — exposure effect SURVIVES!")
elif t_p < 0.01 and diff < 0:
    print(f"    >>> LOW-Z MORE DAMAGED — reversed! Selection/evolution?")
else:
    print(f"    >>> NO SIGNIFICANT DIFFERENCE at matched state")

# By susceptibility tercile
print(f"\n  --- Matched comparison by susceptibility tercile ---")
S_lo = HB_FWHM[lo_z][valid_lo][good_match] + logLedd[idx][lo_z][valid_lo][good_match]
s_terc = np.percentile(S_lo, [0, 33, 67, 100])

for t, label in enumerate(['Low S', 'Mid S', 'High S']):
    in_t = (S_lo >= s_terc[t]) & (S_lo < s_terc[t+1] if t < 2 else True)
    if np.sum(in_t) < 30: continue
    d_lo_t = D_lo_matched[in_t]
    d_hi_t = D_hi_matched[in_t]
    diff_t = np.mean(d_hi_t) - np.mean(d_lo_t)
    _, p_t = stats.wilcoxon(d_lo_t - d_hi_t)
    print(f"    {label}: Δdamage = {diff_t:+.4f}, p = {p_t:.2e}, N = {np.sum(in_t)}")

print("\n\nDone.")
