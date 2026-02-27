#!/usr/bin/env python3
"""
SUSCEPTIBILITY COMPARISON: Circular vs Clean vs Middle Ground
==============================================================
S1: log(OIII_EW) + log(HB_FWHM) + logLedd  [CIRCULAR — includes damage channel]
S2: log(HB_FWHM) + logLedd                   [MIDDLE — kinematic + AGN state only]
S3: logMbh + logLedd + logLbol + BLR_FWHM    [CLEAN — no Hβ/OIII at all]
S4: log(HB_FWHM) only                        [PUREST — just kinematics]

Run identical operator map on all four. Compare.
"""

import numpy as np
from astropy.io import fits
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_susceptibility_comparison'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("SUSCEPTIBILITY COMPARISON")
print("=" * 70)

# Load
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

def fix(arr): return arr.astype(arr.dtype.newbyteorder('='))

z = fix(data['Z_DR16Q'])
snr = fix(data['SN_MEDIAN_ALL'])
ebv = fix(data['EBV'])
oiii_flux = fix(data['OIII5007'][:, 4])
hbeta_flux = fix(data['HBETA'][:, 4])
oiii_ew = fix(data['OIII5007'][:, 2])
hbeta_ew = fix(data['HBETA'][:, 2])
oiii_fwhm = fix(data['OIII5007'][:, 3])
hbeta_fwhm = fix(data['HBETA'][:, 3])
civ_fwhm = fix(data['CIV'][:, 3])
mgii_fwhm = fix(data['MGII'][:, 3])
logLbol = fix(data['LOGLBOL'])
logMbh = fix(data['LOGMBH'])
logLedd = fix(data['LOGLEDD_RATIO'])

# Quality cuts
good = ((z > 0.3) & (z < 2.0) & (snr > 3) &
        (oiii_flux > 0) & (hbeta_flux > 0) &
        np.isfinite(oiii_flux) & np.isfinite(hbeta_flux) &
        np.isfinite(logLbol) & np.isfinite(logMbh) & np.isfinite(logLedd) &
        (logLbol > 0) & (logMbh > 0) &
        np.isfinite(oiii_ew) & np.isfinite(hbeta_ew) &
        np.isfinite(oiii_fwhm) & np.isfinite(hbeta_fwhm) &
        (oiii_fwhm > 0) & (hbeta_fwhm > 0) &
        (oiii_ew > 0) & (hbeta_ew > 0))

idx = np.where(good)[0]
N = len(idx)
Z = z[idx]

# Source model
print(f"\n[1] Source model (N={N})...")
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_predict

log_ratio = np.log10(oiii_flux[idx] / hbeta_flux[idx])
source_X = np.column_stack([
    snr[idx], np.log10(np.maximum(oiii_ew[idx], 1e-3)),
    np.log10(np.maximum(hbeta_ew[idx], 1e-3)),
    np.log10(np.maximum(oiii_fwhm[idx], 1)),
    np.log10(np.maximum(hbeta_fwhm[idx], 1)),
    logLbol[idx], logMbh[idx], logLedd[idx], ebv[idx],
])

gbr = GradientBoostingRegressor(n_estimators=100, max_depth=3, learning_rate=0.1,
                                  subsample=0.8, random_state=42)
ratio_pred = cross_val_predict(gbr, source_X, log_ratio, cv=5)
residual = log_ratio - ratio_pred
r2 = 1 - np.sum(residual**2) / np.sum((log_ratio - np.mean(log_ratio))**2)
print(f"  R² = {r2:.4f}")

# Z-bin damage
n_bins = 20
z_edges = np.percentile(Z, np.linspace(0, 100, n_bins + 1))
z_bin = np.clip(np.digitize(Z, z_edges) - 1, 0, n_bins - 1)
damage = np.full(N, np.nan)
for b in range(n_bins):
    m = z_bin == b
    if np.sum(m) < 10: continue
    med = np.median(residual[m])
    mad = max(1.4826 * np.median(np.abs(residual[m] - med)), 1e-6)
    damage[m] = np.abs((residual[m] - med) / mad)

# ============================================================
# Build susceptibilities
# ============================================================
print("\n[2] Building susceptibilities...")

def standardize(x):
    v = np.isfinite(x)
    out = np.full_like(x, np.nan)
    out[v] = (x[v] - np.mean(x[v])) / np.std(x[v])
    return out

# S1: CIRCULAR (original that gave ρ=1.000)
S1 = (np.log10(np.maximum(oiii_ew[idx], 1e-3)) +
      np.log10(np.maximum(hbeta_fwhm[idx], 1)) +
      logLedd[idx])
S1 = standardize(S1)

# S2: MIDDLE GROUND (HB_FWHM + logLedd — kinematic + AGN state)
S2 = (np.log10(np.maximum(hbeta_fwhm[idx], 1)) + logLedd[idx])
S2 = standardize(S2)

# S3: CLEAN (no Hβ/OIII at all)
blr_fwhm = np.full(N, np.nan)
cv = (civ_fwhm[idx] > 0) & np.isfinite(civ_fwhm[idx])
mv = (~cv) & (mgii_fwhm[idx] > 0) & np.isfinite(mgii_fwhm[idx])
blr_fwhm[cv] = np.log10(civ_fwhm[idx][cv])
blr_fwhm[mv] = np.log10(mgii_fwhm[idx][mv])
S3_parts = [standardize(logMbh[idx]), standardize(logLedd[idx]), standardize(logLbol[idx])]
S3 = np.nanmean(S3_parts, axis=0)
blr_v = np.isfinite(blr_fwhm)
blr_std = standardize(blr_fwhm)
S3[blr_v] = (S3[blr_v] * 3 + blr_std[blr_v]) / 4

# S4: PUREST (just HB_FWHM — kinematic only)
S4 = standardize(np.log10(np.maximum(hbeta_fwhm[idx], 1)))

# S5: OIII_EW only (the suspected circular driver — isolate it)
S5 = standardize(np.log10(np.maximum(oiii_ew[idx], 1e-3)))

susceptibilities = {
    'S1_circular (OIII_EW+HB_FWHM+Ledd)': S1,
    'S2_middle (HB_FWHM+Ledd)': S2,
    'S3_clean (Mbh+Ledd+Lbol+BLR_FWHM)': S3,
    'S4_purest (HB_FWHM only)': S4,
    'S5_isolated (OIII_EW only)': S5,
}

# ============================================================
# Run operator map on each
# ============================================================
print("\n" + "=" * 70)
all_results = {}

for s_name, S in susceptibilities.items():
    print(f"\n{'='*60}")
    print(f"  {s_name}")
    print(f"{'='*60}")
    
    valid = np.isfinite(S) & np.isfinite(damage)
    Z_v = Z[valid]
    D_v = damage[valid]
    S_v = S[valid]
    
    # Decile bins
    s_edges = np.percentile(S_v, np.linspace(0, 100, 11))
    s_bin = np.clip(np.digitize(S_v, s_edges) - 1, 0, 9)
    
    slopes = []
    for d in range(10):
        in_d = s_bin == d
        if np.sum(in_d) < 50:
            slopes.append(np.nan)
            continue
        sl, _, _, p, _ = stats.linregress(Z_v[in_d], D_v[in_d])
        slopes.append(sl)
    
    valid_sl = [(i, slopes[i]) for i in range(10) if np.isfinite(slopes[i])]
    rho, p_rho = stats.spearmanr([v[0] for v in valid_sl], [v[1] for v in valid_sl])
    
    # Interaction
    exposure = (Z_v - np.mean(Z_v)) / np.std(Z_v)
    interact = S_v * exposure
    rho_int, p_int = stats.spearmanr(interact, D_v)
    
    print(f"  Slopes: {' → '.join([f'{s:+.3f}' for s in slopes if np.isfinite(s)])}")
    print(f"  Slope vs decile: ρ = {rho:+.3f}, p = {p_rho:.4f}")
    print(f"  Interaction S×z: ρ = {rho_int:+.4f}, p = {p_int:.2e}")
    print(f"  N = {np.sum(valid)}")
    
    # Sign inversion check
    neg_slopes = sum(1 for s in slopes if np.isfinite(s) and s < 0)
    pos_slopes = sum(1 for s in slopes if np.isfinite(s) and s > 0)
    first_slope = next((s for s in slopes if np.isfinite(s)), 0)
    last_slope = next((s for s in reversed(slopes) if np.isfinite(s)), 0)
    sign_flip = (first_slope < 0 and last_slope > 0) or (first_slope > 0 and last_slope < 0)
    
    print(f"  Sign inversion: {'YES' if sign_flip else 'NO'} ({neg_slopes} neg, {pos_slopes} pos)")
    print(f"  First decile slope: {first_slope:+.4f}, Last: {last_slope:+.4f}")
    
    # Tercile crossing curves
    s_terc = np.percentile(S_v, [0, 33, 67, 100])
    terc_slopes = []
    for t in range(3):
        in_t = (S_v >= s_terc[t]) & (S_v < s_terc[t+1] if t < 2 else True)
        sl, _, _, p, se = stats.linregress(Z_v[in_t], D_v[in_t])
        terc_slopes.append(sl)
        label = ['Low', 'Mid', 'High'][t]
        print(f"    {label} S tercile: slope = {sl:+.4f}, p = {p:.2e}, N = {np.sum(in_t)}")
    
    all_results[s_name] = {
        'slopes': [float(s) for s in slopes if np.isfinite(s)],
        'rho_decile': float(rho),
        'p_decile': float(p_rho),
        'rho_interaction': float(rho_int),
        'p_interaction': float(p_int),
        'sign_inversion': sign_flip,
        'tercile_slopes': [float(s) for s in terc_slopes],
    }

# ============================================================
# COMPARISON TABLE
# ============================================================
print("\n" + "=" * 70)
print("COMPARISON TABLE")
print("=" * 70)
print(f"\n  {'Susceptibility':<45} {'ρ(decile)':>10} {'p':>10} {'ρ(S×z)':>10} {'p':>12} {'Sign flip':>10}")
print("  " + "-" * 100)

for s_name, r in all_results.items():
    print(f"  {s_name:<45} {r['rho_decile']:>+10.3f} {r['p_decile']:>10.4f} {r['rho_interaction']:>+10.4f} {r['p_interaction']:>12.2e} {'YES' if r['sign_inversion'] else 'no':>10}")

# ============================================================
# DECOMPOSITION: What fraction of S1's power comes from OIII_EW?
# ============================================================
print("\n" + "=" * 70)
print("DECOMPOSITION: Which ingredient drives S1?")
print("=" * 70)

# Correlation between each ingredient and the interaction signal
valid_all = np.isfinite(S1) & np.isfinite(damage)
Z_va = Z[valid_all]
D_va = damage[valid_all]
exposure_va = (Z_va - np.mean(Z_va)) / np.std(Z_va)

ingredients = {
    'log(OIII_EW)': standardize(np.log10(np.maximum(oiii_ew[idx], 1e-3)))[valid_all],
    'log(HB_FWHM)': standardize(np.log10(np.maximum(hbeta_fwhm[idx], 1)))[valid_all],
    'logLedd': standardize(logLedd[idx])[valid_all],
    'log(HB_EW)': standardize(np.log10(np.maximum(hbeta_ew[idx], 1e-3)))[valid_all],
    'log(OIII_FWHM)': standardize(np.log10(np.maximum(oiii_fwhm[idx], 1)))[valid_all],
    'logMbh': standardize(logMbh[idx])[valid_all],
    'logLbol': standardize(logLbol[idx])[valid_all],
    'SNR': standardize(snr[idx])[valid_all],
}

print(f"\n  {'Ingredient':<20} {'ρ(ingred×z, damage)':>20} {'p':>12} {'ρ(ingred, damage)':>18} {'p':>12}")
print("  " + "-" * 85)

for name, vals in ingredients.items():
    v = np.isfinite(vals)
    interact_i = vals[v] * exposure_va[v]
    rho_i, p_i = stats.spearmanr(interact_i, D_va[v])
    rho_m, p_m = stats.spearmanr(vals[v], D_va[v])
    marker = " <<<" if abs(rho_i) > 0.02 else ""
    print(f"  {name:<20} {rho_i:>+20.4f} {p_i:>12.2e} {rho_m:>+18.4f} {p_m:>12.2e}{marker}")

# ============================================================
# THE KEY QUESTION: Does HB_FWHM alone show the sign inversion?
# ============================================================
print("\n" + "=" * 70)
print("KEY QUESTION: Does HB_FWHM alone produce sign inversion?")
print("=" * 70)

S_hbf = S4  # HB_FWHM only
valid_hbf = np.isfinite(S_hbf) & np.isfinite(damage)
Z_hbf = Z[valid_hbf]
D_hbf = damage[valid_hbf]
S_hbf_v = S_hbf[valid_hbf]

s_edges_hbf = np.percentile(S_hbf_v, np.linspace(0, 100, 11))
s_bin_hbf = np.clip(np.digitize(S_hbf_v, s_edges_hbf) - 1, 0, 9)

print(f"\n  {'Decile':<8} {'HB_FWHM range':>25} {'N':>6} {'slope':>10} {'p':>12} {'mean_dam':>10}")
print("  " + "-" * 80)

for d in range(10):
    in_d = s_bin_hbf == d
    if np.sum(in_d) < 50: continue
    sl, _, _, p, _ = stats.linregress(Z_hbf[in_d], D_hbf[in_d])
    s_lo, s_hi = s_edges_hbf[d], s_edges_hbf[d+1]
    label = "NARROW" if d < 3 else "BROAD" if d > 6 else "mid"
    print(f"  D{d+1:<7} [{s_lo:>+10.3f}, {s_hi:>+10.3f}] {np.sum(in_d):>6} {sl:>+10.4f} {p:>12.2e} {np.mean(D_hbf[in_d]):>10.3f} {label}")

# Save
with open(f'{OUTDIR}/susceptibility_comparison.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print(f"\n  Saved to {OUTDIR}/")
print("\nDone.")
