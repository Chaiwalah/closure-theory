#!/usr/bin/env python3
"""
OPERATOR MAP: Susceptibility-Gated Exposure Law
=================================================
The one-plot version of the operator:
  "Complexity turns distance into damage"

Bin by susceptibility deciles, fit damage vs z within each,
show that the exposure slope steepens with susceptibility.
"""

import numpy as np
from astropy.io import fits
from scipy import stats
from scipy.optimize import curve_fit
import json, os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_operator_map'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("OPERATOR MAP: SUSCEPTIBILITY-GATED EXPOSURE LAW")
print("=" * 70)

# ============================================================
# 1. Load & prepare (same pipeline as tail_residual_test)
# ============================================================
print("\n[1] Loading data...")
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

def fix(arr):
    return arr.astype(arr.dtype.newbyteorder('='))

z = fix(data['Z_DR16Q'])
snr = fix(data['SN_MEDIAN_ALL'])
ebv = fix(data['EBV'])

oiii_flux = fix(data['OIII5007'][:, 4])
hbeta_flux = fix(data['HBETA'][:, 4])
nii_flux = fix(data['NII6585'][:, 4])
halpha_flux = fix(data['HALPHA'][:, 4])

oiii_ew = fix(data['OIII5007'][:, 2])
hbeta_ew = fix(data['HBETA'][:, 2])
oiii_fwhm = fix(data['OIII5007'][:, 3])
hbeta_fwhm = fix(data['HBETA'][:, 3])

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
print(f"  Clean sample: {N}")

Z = z[idx]
SNR = snr[idx]
OIII = oiii_flux[idx]
HB = hbeta_flux[idx]
OIII_EW = oiii_ew[idx]
HB_EW = hbeta_ew[idx]
OIII_FWHM = oiii_fwhm[idx]
HB_FWHM = hbeta_fwhm[idx]
LLBOL = logLbol[idx]
LMBH = logMbh[idx]
LLEDD = logLedd[idx]

log_ratio = np.log10(OIII / HB)

# ============================================================
# 2. Source model (GBR, cross-validated)
# ============================================================
print("\n[2] Fitting source model...")
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_predict

source_X = np.column_stack([
    SNR,
    np.log10(np.maximum(OIII_EW, 1e-3)),
    np.log10(np.maximum(HB_EW, 1e-3)),
    np.log10(np.maximum(OIII_FWHM, 1)),
    np.log10(np.maximum(HB_FWHM, 1)),
    LLBOL, LMBH, LLEDD, fix(data['EBV'][idx]),
])

gbr = GradientBoostingRegressor(n_estimators=200, max_depth=4, learning_rate=0.1,
                                  subsample=0.8, random_state=42)
ratio_pred = cross_val_predict(gbr, source_X, log_ratio, cv=5)
residual = log_ratio - ratio_pred

r2 = 1 - np.sum(residual**2) / np.sum((log_ratio - np.mean(log_ratio))**2)
print(f"  Source model R² = {r2:.4f}")
print(f"  Residual std = {np.std(residual):.4f}")

# ============================================================
# 3. Susceptibility score
# ============================================================
print("\n[3] Computing susceptibility...")

# Multiple susceptibility definitions to test robustness
# S1: OIII EW complexity + Hβ FWHM + Eddington ratio (from interaction test)
S1 = (np.log10(np.maximum(OIII_EW, 1e-3)) +
      np.log10(np.maximum(HB_FWHM, 1)) +
      LLEDD)
S1 = (S1 - np.mean(S1)) / np.std(S1)

# S2: Information richness = log(OIII_EW/HB_EW) + log(OIII_FWHM/HB_FWHM) 
# (how different the two lines are = more diagnostic info)
S2 = (np.log10(np.maximum(OIII_EW, 1e-3) / np.maximum(HB_EW, 1e-3)) +
      np.log10(np.maximum(OIII_FWHM, 1) / np.maximum(HB_FWHM, 1)))
S2 = (S2 - np.mean(S2)) / np.std(S2)

# S3: Pure Eddington ratio (simplest physical proxy for "how active/complex")
S3 = (LLEDD - np.mean(LLEDD)) / np.std(LLEDD)

susceptibilities = {'S1_composite': S1, 'S2_info_richness': S2, 'S3_eddington': S3}

# ============================================================
# 4. RESIDUAL DAMAGE (z-bin normalized)
# ============================================================
print("\n[4] Computing residual damage...")

n_bins = 20
z_edges = np.percentile(Z, np.linspace(0, 100, n_bins + 1))
z_bin = np.digitize(Z, z_edges) - 1
z_bin = np.clip(z_bin, 0, n_bins - 1)

damage = np.full(N, np.nan)
for b in range(n_bins):
    m = z_bin == b
    if np.sum(m) < 10:
        continue
    med = np.median(residual[m])
    mad = 1.4826 * np.median(np.abs(residual[m] - med))
    if mad < 1e-6:
        mad = np.std(residual[m])
    damage[m] = np.abs((residual[m] - med) / mad)

valid = np.isfinite(damage)
print(f"  Valid: {np.sum(valid)}")

# ============================================================
# 5. THE OPERATOR MAP: Susceptibility deciles × exposure slope
# ============================================================
print("\n" + "=" * 70)
print("THE OPERATOR MAP")
print("=" * 70)

for s_name, S in susceptibilities.items():
    print(f"\n{'='*50}")
    print(f"  Susceptibility: {s_name}")
    print(f"{'='*50}")
    
    # Decile bins
    s_edges = np.percentile(S[valid], np.linspace(0, 100, 11))
    s_bin = np.digitize(S[valid], s_edges) - 1
    s_bin = np.clip(s_bin, 0, 9)
    
    Z_v = Z[valid]
    D_v = damage[valid]
    S_v = S[valid]
    
    print(f"\n  {'Decile':<8} {'S range':>20} {'N':>6} {'slope(dam vs z)':>16} {'p':>12} {'mean_dam':>10} {'Interpretation'}")
    print("  " + "-" * 90)
    
    slopes = []
    slope_errs = []
    s_centers = []
    mean_damages = []
    
    for d in range(10):
        in_dec = s_bin == d
        if np.sum(in_dec) < 30:
            slopes.append(np.nan)
            slope_errs.append(np.nan)
            s_centers.append(np.nan)
            mean_damages.append(np.nan)
            continue
        
        z_dec = Z_v[in_dec]
        dam_dec = D_v[in_dec]
        
        # Linear fit: damage = a*z + b
        slope, intercept, r, p, se = stats.linregress(z_dec, dam_dec)
        rho, rho_p = stats.spearmanr(z_dec, dam_dec)
        
        s_lo, s_hi = s_edges[d], s_edges[d+1]
        s_center = 0.5 * (s_lo + s_hi)
        
        slopes.append(slope)
        slope_errs.append(se)
        s_centers.append(s_center)
        mean_damages.append(np.mean(dam_dec))
        
        label = "LOW suscept" if d < 3 else "HIGH suscept" if d > 6 else "mid"
        print(f"  D{d+1:<7} [{s_lo:>+8.3f}, {s_hi:>+8.3f}] {np.sum(in_dec):>6} {slope:>+16.4f} {p:>12.2e} {np.mean(dam_dec):>10.3f} {label}")
    
    # The money test: does slope correlate with susceptibility?
    valid_slopes = np.array([(i, slopes[i]) for i in range(10) if np.isfinite(slopes[i])])
    if len(valid_slopes) >= 5:
        dec_idx = valid_slopes[:, 0]
        dec_slopes = valid_slopes[:, 1]
        rho_slope, p_slope = stats.spearmanr(dec_idx, dec_slopes)
        
        print(f"\n  SLOPE vs SUSCEPTIBILITY DECILE:")
        print(f"    Spearman ρ = {rho_slope:.3f}, p = {p_slope:.4f}")
        print(f"    Slopes: {' → '.join([f'{s:+.4f}' for s in slopes if np.isfinite(s)])}")
        
        if rho_slope > 0.3 and p_slope < 0.05:
            print(f"    >>> ✅ OPERATOR CONFIRMED: exposure slope STEEPENS with susceptibility!")
            print(f"    >>> 'Complexity turns distance into damage'")
        elif rho_slope > 0 and p_slope < 0.1:
            print(f"    >>> ⚠️ MARGINAL: trend in right direction but not definitive")
        elif rho_slope < -0.3 and p_slope < 0.05:
            print(f"    >>> ❌ REVERSED: low-susceptibility sources show MORE distance effect")
        else:
            print(f"    >>> ⏸️ INCONCLUSIVE: no clear trend")

# ============================================================
# 6. SIGMOID FITS per susceptibility tercile
# ============================================================
print("\n" + "=" * 70)
print("SIGMOID FITS BY SUSCEPTIBILITY TERCILE")
print("=" * 70)

def sigmoid(z, a, z0, k, b):
    return a / (1 + np.exp(-k * (z - z0))) + b

S_use = S1  # primary susceptibility
s_tercile_edges = np.percentile(S_use[valid], [0, 33, 67, 100])
tercile_labels = ['Low susceptibility', 'Mid susceptibility', 'High susceptibility']

sigmoid_results = []
for t in range(3):
    in_t = (S_use[valid] >= s_tercile_edges[t]) & (S_use[valid] < s_tercile_edges[t+1])
    if t == 2:
        in_t = (S_use[valid] >= s_tercile_edges[t])
    
    z_t = Z[valid][in_t]
    d_t = damage[valid][in_t]
    
    # Bin into 10 z-bins for cleaner fit
    z_bin_edges = np.percentile(z_t, np.linspace(0, 100, 11))
    z_bin_centers = []
    d_bin_means = []
    d_bin_stds = []
    for b in range(10):
        bm = (z_t >= z_bin_edges[b]) & (z_t < z_bin_edges[b+1])
        if b == 9:
            bm = (z_t >= z_bin_edges[b])
        if np.sum(bm) < 10:
            continue
        z_bin_centers.append(np.mean(z_t[bm]))
        d_bin_means.append(np.mean(d_t[bm]))
        d_bin_stds.append(np.std(d_t[bm]) / np.sqrt(np.sum(bm)))
    
    z_bc = np.array(z_bin_centers)
    d_bm = np.array(d_bin_means)
    
    # Linear fit
    sl, inter, r, p, se = stats.linregress(z_bc, d_bm)
    
    # Try sigmoid
    try:
        popt, pcov = curve_fit(sigmoid, z_bc, d_bm, p0=[0.2, 1.0, 5.0, 0.7], maxfev=5000)
        sig_pred = sigmoid(z_bc, *popt)
        ss_sig = np.sum((d_bm - sig_pred)**2)
        ss_lin = np.sum((d_bm - (sl*z_bc + inter))**2)
        f_stat = ((ss_lin - ss_sig) / 2) / (ss_sig / max(len(z_bc) - 4, 1))
        sig_fit = f"z₀={popt[1]:.2f}, k={popt[2]:.2f}, F={f_stat:.1f}"
    except:
        sig_fit = "failed"
        f_stat = 0
    
    print(f"\n  {tercile_labels[t]} (N={np.sum(in_t)}):")
    print(f"    Linear slope: {sl:+.4f} ± {se:.4f} (p={p:.2e})")
    print(f"    Sigmoid: {sig_fit}")
    print(f"    z-binned means: {' → '.join([f'{d:.3f}' for d in d_bm])}")
    
    sigmoid_results.append({
        'tercile': tercile_labels[t],
        'n': int(np.sum(in_t)),
        'linear_slope': float(sl),
        'linear_p': float(p),
        'z_centers': [float(x) for x in z_bc],
        'damage_means': [float(x) for x in d_bm],
    })

# ============================================================
# 7. LOCKED CONTROL: Same analysis on NII/Hα
# ============================================================
print("\n" + "=" * 70)
print("LOCKED CONTROL: NII/Hα OPERATOR MAP")
print("=" * 70)

nii_good = (good & (nii_flux > 0) & (halpha_flux > 0) & 
            np.isfinite(nii_flux) & np.isfinite(halpha_flux))
nii_idx = np.where(nii_good)[0]

if len(nii_idx) > 500:
    Z_nii = z[nii_idx]
    nii_ratio = np.log10(nii_flux[nii_idx] / halpha_flux[nii_idx])
    S_nii = S1[np.isin(idx, nii_idx)]  # susceptibility for matching objects
    
    # Wait — S1 is indexed differently. Recompute for NII sample
    S_nii_raw = (np.log10(np.maximum(oiii_ew[nii_idx], 1e-3)) +
                 np.log10(np.maximum(hbeta_fwhm[nii_idx], 1)) +
                 logLedd[nii_idx])
    # Only use objects where these are valid
    nii_valid = np.isfinite(S_nii_raw) & np.isfinite(Z_nii) & np.isfinite(nii_ratio)
    
    if np.sum(nii_valid) > 200:
        Z_nv = Z_nii[nii_valid]
        R_nv = nii_ratio[nii_valid]
        S_nv = S_nii_raw[nii_valid]
        S_nv = (S_nv - np.mean(S_nv)) / np.std(S_nv)
        
        # Z-detrend
        nii_zbin = np.digitize(Z_nv, np.percentile(Z_nv, np.linspace(0, 100, 11))) - 1
        nii_zbin = np.clip(nii_zbin, 0, 9)
        nii_damage = np.full(len(Z_nv), np.nan)
        for b in range(10):
            m = nii_zbin == b
            if np.sum(m) < 5:
                continue
            med = np.median(R_nv[m])
            mad = 1.4826 * np.median(np.abs(R_nv[m] - med))
            if mad < 1e-6:
                mad = np.std(R_nv[m])
            nii_damage[m] = np.abs((R_nv[m] - med) / mad)
        
        nii_v = np.isfinite(nii_damage)
        
        # Susceptibility terciles
        s_terc = np.percentile(S_nv[nii_v], [0, 33, 67, 100])
        print(f"  Locked sample: {np.sum(nii_v)} objects")
        for t, label in enumerate(['Low S', 'Mid S', 'High S']):
            in_t = (S_nv[nii_v] >= s_terc[t]) & (S_nv[nii_v] < s_terc[t+1])
            if t == 2:
                in_t = S_nv[nii_v] >= s_terc[t]
            if np.sum(in_t) < 30:
                continue
            sl, inter, r, p, se = stats.linregress(Z_nv[nii_v][in_t], nii_damage[nii_v][in_t])
            print(f"  {label}: slope = {sl:+.4f} ± {se:.4f}, p = {p:.2e}, N = {np.sum(in_t)}")
    else:
        print(f"  Too few valid NII objects ({np.sum(nii_valid)})")
else:
    print(f"  Too few NII objects ({len(nii_idx)})")

# ============================================================
# 8. CROSS-SUSCEPTIBILITY: Does the pattern hold for SII?
# ============================================================
print("\n" + "=" * 70)
print("CROSS-CHECK: SII (highest sensitivity doublet)")
print("=" * 70)

sii_flux_arr = fix(data['SII6718'][:, 4])
sii_good = (good & (sii_flux_arr > 0) & np.isfinite(sii_flux_arr))
sii_idx = np.where(sii_good)[0]

if len(sii_idx) > 500:
    Z_sii = z[sii_idx]
    # SII/OIII ratio (high-sensitivity / locked)
    sii_ratio = np.log10(sii_flux_arr[sii_idx] / oiii_flux[sii_idx])
    
    S_sii = (np.log10(np.maximum(oiii_ew[sii_idx], 1e-3)) +
             np.log10(np.maximum(hbeta_fwhm[sii_idx], 1)) +
             logLedd[sii_idx])
    sii_valid = np.isfinite(S_sii) & np.isfinite(sii_ratio) & np.isfinite(Z_sii)
    
    if np.sum(sii_valid) > 200:
        Z_sv = Z_sii[sii_valid]
        R_sv = sii_ratio[sii_valid]
        S_sv = S_sii[sii_valid]
        S_sv = (S_sv - np.mean(S_sv)) / np.std(S_sv)
        
        # Z-detrend
        sii_zbin = np.digitize(Z_sv, np.percentile(Z_sv, np.linspace(0, 100, 11))) - 1
        sii_zbin = np.clip(sii_zbin, 0, 9)
        sii_damage = np.full(len(Z_sv), np.nan)
        for b in range(10):
            m = sii_zbin == b
            if np.sum(m) < 5:
                continue
            med = np.median(R_sv[m])
            mad = 1.4826 * np.median(np.abs(R_sv[m] - med))
            if mad < 1e-6:
                mad = np.std(R_sv[m])
            sii_damage[m] = np.abs((R_sv[m] - med) / mad)
        
        sii_v = np.isfinite(sii_damage)
        
        # Interaction test
        exposure_sii = (Z_sv[sii_v] - np.mean(Z_sv[sii_v])) / np.std(Z_sv[sii_v])
        interact_sii = S_sv[sii_v] * exposure_sii
        rho_sii, p_sii = stats.spearmanr(interact_sii, sii_damage[sii_v])
        
        print(f"  SII sample: {np.sum(sii_v)} objects")
        print(f"  Interaction (S×z) vs SII damage: ρ = {rho_sii:.4f}, p = {p_sii:.2e}")
        
        # Tercile slopes
        s_terc = np.percentile(S_sv[sii_v], [0, 33, 67, 100])
        for t, label in enumerate(['Low S', 'Mid S', 'High S']):
            in_t = (S_sv[sii_v] >= s_terc[t]) & (S_sv[sii_v] < s_terc[t+1])
            if t == 2:
                in_t = S_sv[sii_v] >= s_terc[t]
            if np.sum(in_t) < 30:
                continue
            sl, inter, r, p, se = stats.linregress(Z_sv[sii_v][in_t], sii_damage[sii_v][in_t])
            print(f"  {label}: slope = {sl:+.4f} ± {se:.4f}, p = {p:.2e}, N = {np.sum(in_t)}")
    else:
        print(f"  Too few valid SII objects ({np.sum(sii_valid)})")
else:
    print(f"  Too few SII objects ({len(sii_idx)})")

# ============================================================
# SAVE
# ============================================================
print("\n" + "=" * 70)
print("SAVING")
print("=" * 70)

# Save the slopes for plotting
with open(f'{OUTDIR}/operator_map_results.json', 'w') as f:
    json.dump({
        'test': 'Operator Map — Susceptibility-Gated Exposure',
        'source_model_R2': float(r2),
        'sigmoid_terciles': sigmoid_results,
    }, f, indent=2)

print(f"  Saved to {OUTDIR}/operator_map_results.json")
print("\nDone.")
