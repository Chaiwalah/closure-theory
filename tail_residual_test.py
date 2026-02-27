#!/usr/bin/env python3
"""
TAIL-CONDITIONAL TEST v2 — SOURCE-RESIDUALIZED
================================================
Stage A: Build source model to predict OIII/Hβ from source properties
Stage B: Define damage as RESIDUAL after source subtraction
Stage C: Rerun tail-conditional on residuals — does Path AUC rise?
Stage D: Interaction test — does residual scale with susceptibility × exposure?
"""

import numpy as np
from astropy.io import fits
from scipy import stats
import healpy as hp
import json, os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_tail_residual'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("TAIL-CONDITIONAL v2 — SOURCE-RESIDUALIZED")
print("=" * 70)

# ============================================================
# 1. Load everything
# ============================================================
print("\n[1] Loading data...")
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

# Fix big-endian for pandas compatibility
def fix_endian(arr):
    return arr.astype(arr.dtype.newbyteorder('='))

z = fix_endian(data['Z_DR16Q'])
ra = fix_endian(data['RA'])
dec = fix_endian(data['DEC'])
snr = fix_endian(data['SN_MEDIAN_ALL'])
ebv = fix_endian(data['EBV'])

# Lines (col 4 = peak flux, col 2 = EW, col 3 = sigma/FWHM)
oiii_flux = fix_endian(data['OIII5007'][:, 4])
hbeta_flux = fix_endian(data['HBETA'][:, 4])
nii_flux = fix_endian(data['NII6585'][:, 4])
halpha_flux = fix_endian(data['HALPHA'][:, 4])
sii_flux = fix_endian(data['SII6718'][:, 4])
oii_flux = fix_endian(data['OII3728'][:, 4])

oiii_ew = fix_endian(data['OIII5007'][:, 2])
hbeta_ew = fix_endian(data['HBETA'][:, 2])
oiii_fwhm = fix_endian(data['OIII5007'][:, 3])
hbeta_fwhm = fix_endian(data['HBETA'][:, 3])

logLbol = fix_endian(data['LOGLBOL'])
logMbh = fix_endian(data['LOGMBH'])
logLedd = fix_endian(data['LOGLEDD_RATIO'])

# Planck κ
print("  Loading Planck κ...")
klm = hp.read_alm('/root/clawd/COM_Lensing_4096_R3.00/MV/dat_klm.fits')
mask = hp.read_map('/root/clawd/COM_Lensing_4096_R3.00/mask.fits', dtype=float)
NSIDE = 1024

scales = [10, 30, 60, 120, 240]
kappa_maps = {}
for fwhm in scales:
    fwhm_rad = np.radians(fwhm / 60.0)
    klm_s = hp.almxfl(klm, hp.gauss_beam(fwhm_rad, lmax=hp.Alm.getlmax(len(klm))))
    kappa_maps[fwhm] = hp.alm2map(klm_s, NSIDE, verbose=False)

mask_ds = hp.ud_grade(mask, NSIDE)
theta = np.radians(90.0 - dec)
phi = np.radians(ra)
pixels = hp.ang2pix(NSIDE, theta, phi)
in_mask = mask_ds[pixels] > 0.5

# ============================================================
# 2. Quality cuts
# ============================================================
print("\n[2] Quality cuts...")
good = ((z > 0.3) & (z < 2.0) & (snr > 3) & in_mask & 
        (oiii_flux > 0) & (hbeta_flux > 0) & 
        np.isfinite(oiii_flux) & np.isfinite(hbeta_flux) &
        np.isfinite(logLbol) & np.isfinite(logMbh) & np.isfinite(logLedd) &
        (logLbol > 0) & (logMbh > 0) &
        np.isfinite(oiii_ew) & np.isfinite(hbeta_ew) &
        np.isfinite(oiii_fwhm) & np.isfinite(hbeta_fwhm) &
        (oiii_fwhm > 0) & (hbeta_fwhm > 0))

idx = np.where(good)[0]
N = len(idx)
print(f"  Clean sample: {N}")

# Extract arrays for clean sample
Z = z[idx]
SNR = snr[idx]
EBV = ebv[idx]
OIII = oiii_flux[idx]
HB = hbeta_flux[idx]
OIII_EW = oiii_ew[idx]
HB_EW = hbeta_ew[idx]
OIII_FWHM = oiii_fwhm[idx]
HB_FWHM = hbeta_fwhm[idx]
LLBOL = logLbol[idx]
LMBH = logMbh[idx]
LLEDD = logLedd[idx]
PIX = pixels[idx]

# κ features
kappa = {}
for s in scales:
    kappa[s] = kappa_maps[s][PIX]
kappa_mean = np.mean([kappa[s] for s in scales], axis=0)
kappa_min = np.min([kappa[s] for s in scales], axis=0)
kappa_range = np.max([kappa[s] for s in scales], axis=0) - kappa_min
kappa_grad = kappa[10] - kappa[240]

# Target: log(OIII/Hβ)
log_ratio = np.log10(OIII / HB)

# ============================================================
# 3. STAGE A — Source model (predict ratio from source properties)
# ============================================================
print("\n" + "=" * 70)
print("STAGE A: SOURCE MODEL")
print("=" * 70)

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_predict
from sklearn.preprocessing import StandardScaler

# Source features ONLY (no path, no z beyond what's needed)
# We deliberately EXCLUDE z from the source model to preserve path signal
# But we include z-dependent source properties (Lbol, Mbh, Ledd which already encode z)
source_X = np.column_stack([
    SNR,
    np.log10(np.maximum(OIII_EW, 1e-3)),
    np.log10(np.maximum(HB_EW, 1e-3)),
    np.log10(np.maximum(OIII_FWHM, 1)),
    np.log10(np.maximum(HB_FWHM, 1)),
    LLBOL,
    LMBH,
    LLEDD,
    EBV,
])
source_names = ['SNR', 'log_OIII_EW', 'log_HB_EW', 'log_OIII_FWHM', 'log_HB_FWHM',
                'logLbol', 'logMbh', 'logLedd', 'EBV']

# Use cross-validated predictions to avoid overfitting
print("  Fitting GBR source model (cross-validated)...")
gbr = GradientBoostingRegressor(n_estimators=200, max_depth=4, learning_rate=0.1,
                                  subsample=0.8, random_state=42)
ratio_pred_gbr = cross_val_predict(gbr, source_X, log_ratio, cv=5)

# Also fit a simple linear model for comparison
scaler = StandardScaler()
source_X_scaled = scaler.fit_transform(source_X)
ratio_pred_lin = cross_val_predict(LinearRegression(), source_X_scaled, log_ratio, cv=5)

# Residuals
resid_gbr = log_ratio - ratio_pred_gbr
resid_lin = log_ratio - ratio_pred_lin

r2_gbr = 1 - np.sum(resid_gbr**2) / np.sum((log_ratio - np.mean(log_ratio))**2)
r2_lin = 1 - np.sum(resid_lin**2) / np.sum((log_ratio - np.mean(log_ratio))**2)

print(f"  Linear source model R² = {r2_lin:.4f}")
print(f"  GBR source model R² = {r2_gbr:.4f}")
print(f"  Residual std (GBR): {np.std(resid_gbr):.4f}")
print(f"  Residual std (Lin): {np.std(resid_lin):.4f}")

# Use GBR residuals (better source subtraction = purer path signal)
residual = resid_gbr

# ============================================================
# 4. STAGE B — Residual damage score
# ============================================================
print("\n" + "=" * 70)
print("STAGE B: RESIDUAL DAMAGE")
print("=" * 70)

# Z-bin normalize the RESIDUAL (remove any remaining z-trend)
n_bins = 20
z_edges = np.percentile(Z, np.linspace(0, 100, n_bins + 1))
z_bin = np.digitize(Z, z_edges) - 1
z_bin = np.clip(z_bin, 0, n_bins - 1)

damage_resid = np.full(N, np.nan)
for b in range(n_bins):
    m = z_bin == b
    if np.sum(m) < 10:
        continue
    med = np.median(residual[m])
    mad = 1.4826 * np.median(np.abs(residual[m] - med))
    if mad < 1e-6:
        mad = np.std(residual[m])
    damage_resid[m] = np.abs((residual[m] - med) / mad)

valid = np.isfinite(damage_resid)
print(f"  Valid residual damage scores: {np.sum(valid)}")
print(f"  Mean: {np.nanmean(damage_resid):.3f}, Std: {np.nanstd(damage_resid):.3f}")

# ============================================================
# 5. STAGE C — Tail-conditional on residuals
# ============================================================
print("\n" + "=" * 70)
print("STAGE C: TAIL PREDICTION ON RESIDUALS")
print("=" * 70)

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_auc_score

# Label tails within z-bins
tail_5 = np.zeros(N, dtype=bool)
tail_1 = np.zeros(N, dtype=bool)
for b in range(n_bins):
    m = (z_bin == b) & valid
    if np.sum(m) < 20:
        continue
    t5 = np.percentile(damage_resid[m], 95)
    t1 = np.percentile(damage_resid[m], 99)
    tail_5[m & (damage_resid >= t5)] = True
    tail_1[m & (damage_resid >= t1)] = True

print(f"  Top 5% tail: {np.sum(tail_5)} objects")
print(f"  Top 1% tail: {np.sum(tail_1)} objects")

# Feature matrices
source_feats = source_X[valid]
path_feats = np.column_stack([kappa_mean[valid], kappa_min[valid], 
                               kappa_range[valid], kappa_grad[valid],
                               kappa[10][valid], kappa[30][valid], 
                               kappa[60][valid], kappa[120][valid], kappa[240][valid]])
path_names = ['κ_mean', 'κ_min', 'κ_range', 'κ_grad', 'κ_10', 'κ_30', 'κ_60', 'κ_120', 'κ_240']

all_feats = np.column_stack([source_feats, path_feats])
all_names = source_names + path_names

results = {}

for tail_arr, tail_name in [(tail_5[valid], 'Top 5%'), (tail_1[valid], 'Top 1%')]:
    y = tail_arr.astype(int)
    if y.sum() < 10:
        print(f"\n  {tail_name}: too few tail objects ({y.sum()}), skipping")
        continue
    
    print(f"\n--- {tail_name} RESIDUAL tail ---")
    
    for feat_name, X_raw in [
        ('Source', source_feats.copy()),
        ('Path', path_feats.copy()),
        ('All', all_feats.copy()),
    ]:
        # Impute NaN
        for j in range(X_raw.shape[1]):
            nan_m = np.isnan(X_raw[:, j]) | np.isinf(X_raw[:, j])
            if nan_m.any():
                X_raw[nan_m, j] = np.nanmedian(X_raw[:, j])
        
        sc = StandardScaler()
        X_sc = sc.fit_transform(X_raw)
        
        clf = LogisticRegression(max_iter=1000, class_weight='balanced', C=0.1)
        try:
            cv_auc = cross_val_score(clf, X_sc, y, cv=5, scoring='roc_auc')
            auc = np.mean(cv_auc)
            auc_std = np.std(cv_auc)
        except:
            auc, auc_std = 0.5, 0.0
        
        print(f"  {feat_name:<10} CV-AUC: {auc:.3f} ± {auc_std:.3f}")
        
        results[f'{tail_name}_{feat_name}'] = {
            'cv_auc': float(auc),
            'cv_auc_std': float(auc_std),
        }
        
        # Top features for All model
        if feat_name == 'All':
            clf.fit(X_sc, y)
            coef_abs = np.abs(clf.coef_[0])
            top_idx = np.argsort(coef_abs)[::-1][:8]
            print(f"  Top predictors:")
            for rank, i in enumerate(top_idx):
                name = all_names[i] if i < len(all_names) else f'feat_{i}'
                print(f"    {rank+1}. {name:<20} coef={clf.coef_[0][i]:+.4f}")
    
    # Verdict
    src_auc = results[f'{tail_name}_Source']['cv_auc']
    path_auc = results[f'{tail_name}_Path']['cv_auc']
    print(f"\n  Source AUC: {src_auc:.3f} vs Path AUC: {path_auc:.3f}")
    if path_auc > src_auc + 0.02:
        print(f"  >>> PATH NOW DOMINATES after source subtraction!")
    elif src_auc > path_auc + 0.02:
        print(f"  >>> SOURCE STILL DOMINATES even after subtraction")
    else:
        print(f"  >>> ROUGHLY EQUAL / both near chance")

# ============================================================
# 6. STAGE D — Interaction test: susceptibility × exposure
# ============================================================
print("\n" + "=" * 70)
print("STAGE D: SUSCEPTIBILITY × EXPOSURE INTERACTION")
print("=" * 70)

# Define susceptibility: how "complex" is the source?
# Use OIII EW breadth + Hβ FWHM + Eddington ratio as proxy
susceptibility = (np.log10(np.maximum(OIII_EW, 1e-3)) + 
                  np.log10(np.maximum(HB_FWHM, 1)) + 
                  LLEDD)
susceptibility = (susceptibility - np.mean(susceptibility)) / np.std(susceptibility)

# Exposure: redshift (path length proxy)
exposure = (Z - np.mean(Z)) / np.std(Z)

# Interaction term
interaction = susceptibility * exposure

# Predict residual damage from: susceptibility, exposure, interaction
X_interact = np.column_stack([susceptibility[valid], exposure[valid], interaction[valid]])
y_damage = damage_resid[valid]
y_valid = np.isfinite(y_damage)

X_i = X_interact[y_valid]
y_i = y_damage[y_valid]

from sklearn.linear_model import LinearRegression as LR
reg = LR().fit(X_i, y_i)

print(f"  Susceptibility coef: {reg.coef_[0]:.4f}")
print(f"  Exposure coef:       {reg.coef_[1]:.4f}")
print(f"  Interaction coef:    {reg.coef_[2]:.4f}")
print(f"  R²:                  {reg.score(X_i, y_i):.6f}")

# Spearman of interaction term with residual damage
rho_int, p_int = stats.spearmanr(interaction[valid][y_valid], y_i)
print(f"\n  Spearman(susceptibility×exposure, residual damage): ρ = {rho_int:.4f}, p = {p_int:.2e}")

if p_int < 0.01 and rho_int > 0:
    print("  >>> INTERACTION SIGNIFICANT — operator couples susceptibility to path length!")
elif p_int < 0.05:
    print("  >>> MARGINAL interaction signal")
else:
    print("  >>> No significant interaction detected")

# ============================================================
# 7. LOCKED CONTROL on residuals
# ============================================================
print("\n" + "=" * 70)
print("LOCKED CONTROL: NII/Hα residual vs path")
print("=" * 70)

nii_good = (nii_flux[idx] > 0) & (halpha_flux[idx] > 0)
nii_good &= np.isfinite(nii_flux[idx]) & np.isfinite(halpha_flux[idx])
nii_idx = nii_good & valid

if np.sum(nii_idx) > 100:
    nii_ratio = np.log10(nii_flux[idx][nii_idx] / halpha_flux[idx][nii_idx])
    nii_kappa = kappa_mean[nii_idx]
    nii_z = Z[nii_idx]
    
    # Z-detrend
    nii_resid = np.zeros_like(nii_ratio)
    nii_zbin = np.digitize(nii_z, np.percentile(nii_z, np.linspace(0, 100, 11))) - 1
    nii_zbin = np.clip(nii_zbin, 0, 9)
    for b in range(10):
        m = nii_zbin == b
        if np.sum(m) < 5:
            continue
        med = np.median(nii_ratio[m])
        nii_resid[m] = nii_ratio[m] - med
    
    rho_nii_k, p_nii_k = stats.spearmanr(nii_kappa, np.abs(nii_resid))
    rho_nii_z, p_nii_z = stats.spearmanr(nii_z, np.abs(nii_resid))
    
    print(f"  Locked sample: {np.sum(nii_idx)} objects")
    print(f"  Spearman(κ_mean, |NII/Hα resid|): ρ = {rho_nii_k:.4f}, p = {p_nii_k:.2e}")
    print(f"  Spearman(z, |NII/Hα resid|):      ρ = {rho_nii_z:.4f}, p = {p_nii_z:.2e}")
    if abs(rho_nii_k) < 0.02 and p_nii_k > 0.05:
        print("  >>> LOCKED CONTROL PASSED")
    else:
        print(f"  >>> LOCKED CONTROL {'MARGINAL' if abs(rho_nii_k) < 0.05 else 'NEEDS ATTENTION'}")
else:
    print(f"  Only {np.sum(nii_idx)} NII/Hα objects — insufficient")

# ============================================================
# 8. KURTOSIS TREND on residuals
# ============================================================
print("\n" + "=" * 70)
print("KURTOSIS TREND (residual damage)")
print("=" * 70)

from scipy.stats import kurtosis as kurt_func

z_terciles = np.percentile(Z[valid], [0, 33, 67, 100])
for i, label in enumerate(['Low-z', 'Mid-z', 'High-z']):
    m = (Z[valid] >= z_terciles[i]) & (Z[valid] < z_terciles[i+1])
    d = damage_resid[valid][m]
    d = d[np.isfinite(d)]
    if len(d) < 50:
        continue
    k = kurt_func(d)
    print(f"  {label}: N={len(d)}, kurtosis={k:.3f}, mean={np.mean(d):.3f}, std={np.std(d):.3f}")

# ============================================================
# SAVE
# ============================================================
final = {
    'test': 'Tail-Conditional v2 (Source-Residualized)',
    'source_model_R2_gbr': float(r2_gbr),
    'source_model_R2_lin': float(r2_lin),
    'n_objects': int(np.sum(valid)),
    'tail_classification': results,
    'interaction': {
        'susceptibility_coef': float(reg.coef_[0]),
        'exposure_coef': float(reg.coef_[1]),
        'interaction_coef': float(reg.coef_[2]),
        'R2': float(reg.score(X_i, y_i)),
        'spearman_rho': float(rho_int),
        'spearman_p': float(p_int),
    },
}

with open(f'{OUTDIR}/tail_residual_results.json', 'w') as f:
    json.dump(final, f, indent=2)

print(f"\n  Saved to {OUTDIR}/tail_residual_results.json")
print("\nDone.")
