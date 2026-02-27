#!/usr/bin/env python3
"""
TAIL-CONDITIONAL MAPPING TEST
================================
Instead of predicting mean damage, predict WHO lands in the tails.

If we can predict tail membership from physical/path variables,
we've found the operator's coordinate.

Also investigates the low-z bimodality to check for mundane regime splits.
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import stats
from scipy.stats import gaussian_kde, kurtosis, skew
import healpy as hp
import json
import os
import warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_tail_conditional'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("TAIL-CONDITIONAL MAPPING TEST")
print("=" * 70)

# ============================================================
# 1. Load data
# ============================================================
print("\n[1] Loading DR16Q + Planck κ...")

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

z = data['Z_DR16Q']
ra = data['RA']
dec = data['DEC']
snr = data['SN_MEDIAN_ALL']
ebv = data['EBV']

# Line measurements (col indices: 0=wave_peak, 1=wave_center, 2=EW, 3=sigma, 4=flux_peak, 5=flux_total)
def get_flux(name):
    return data[name][:, 4]

def get_ew(name):
    return data[name][:, 2]

def get_fwhm(name):
    return data[name][:, 3]

oiii_flux = get_flux('OIII5007')
hbeta_flux = get_flux('HBETA')
nii_flux = get_flux('NII6585')
halpha_flux = get_flux('HALPHA')
sii_flux = get_flux('SII6718')
oii_flux = get_flux('OII3728')

oiii_ew = get_ew('OIII5007')
hbeta_ew = get_ew('HBETA')
oiii_fwhm = get_fwhm('OIII5007')
hbeta_fwhm = get_fwhm('HBETA')

# Additional properties
logLbol = data['LOGLBOL']
logMbh = data['LOGMBH']
logLedd = data['LOGLEDD_RATIO']

# Planck κ at multiple scales
print("  Loading Planck κ maps...")
klm = hp.read_alm('/root/clawd/COM_Lensing_4096_R3.00/MV/dat_klm.fits')
mask = hp.read_map('/root/clawd/COM_Lensing_4096_R3.00/mask.fits', dtype=float)
NSIDE = 1024

scales = [10, 30, 60, 120, 240]
kappa_maps = {}
for fwhm_arcmin in scales:
    fwhm_rad = np.radians(fwhm_arcmin / 60.0)
    klm_s = hp.almxfl(klm, hp.gauss_beam(fwhm_rad, lmax=hp.Alm.getlmax(len(klm))))
    kappa_maps[fwhm_arcmin] = hp.alm2map(klm_s, NSIDE, verbose=False)

mask_ds = hp.ud_grade(mask, NSIDE)

theta = np.radians(90.0 - dec)
phi = np.radians(ra)
pixels = hp.ang2pix(NSIDE, theta, phi)
in_mask = mask_ds[pixels] > 0.5

print(f"  Objects in Planck mask: {np.sum(in_mask)}")

# ============================================================
# 2. Build feature matrix
# ============================================================
print("\n[2] Building feature matrix...")

# Quality cuts
good = (z > 0.3) & (z < 2.0) & (snr > 3) & in_mask
good &= (oiii_flux > 0) & (hbeta_flux > 0)
good &= np.isfinite(oiii_flux) & np.isfinite(hbeta_flux)

# Compute diagnostic ratio + damage score
log_ratio = np.full(len(z), np.nan)
log_ratio[good] = np.log10(oiii_flux[good] / hbeta_flux[good])

# Z-bin residual (damage)
z_good = z[good]
ratio_good = log_ratio[good]

n_bins = 20
z_edges = np.percentile(z_good, np.linspace(0, 100, n_bins + 1))
z_bin_idx = np.digitize(z_good, z_edges) - 1
z_bin_idx = np.clip(z_bin_idx, 0, n_bins - 1)

damage = np.full(len(z_good), np.nan)
for b in range(n_bins):
    in_bin = z_bin_idx == b
    if np.sum(in_bin) < 10:
        continue
    med = np.median(ratio_good[in_bin])
    mad = 1.4826 * np.median(np.abs(ratio_good[in_bin] - med))
    if mad < 1e-6:
        mad = np.std(ratio_good[in_bin])
    damage[in_bin] = np.abs((ratio_good[in_bin] - med) / mad)

# Build feature table for good objects
idx_good = np.where(good)[0]
pix_good = pixels[idx_good]

features = pd.DataFrame({
    'z': z_good,
    'snr': snr[idx_good],
    'ebv': ebv[idx_good],
    'oiii_flux': oiii_flux[idx_good],
    'hbeta_flux': hbeta_flux[idx_good],
    'oiii_ew': oiii_ew[idx_good],
    'hbeta_ew': hbeta_ew[idx_good],
    'oiii_fwhm': oiii_fwhm[idx_good],
    'hbeta_fwhm': hbeta_fwhm[idx_good],
    'log_ratio': ratio_good,
    'logLbol': logLbol[idx_good],
    'logMbh': logMbh[idx_good],
    'logLedd': logLedd[idx_good],
    'damage': damage,
})

# Add κ features
for scale in scales:
    features[f'kappa_{scale}'] = kappa_maps[scale][pix_good]

# Derived κ features
features['kappa_mean'] = features[[f'kappa_{s}' for s in scales]].mean(axis=1)
features['kappa_min'] = features[[f'kappa_{s}' for s in scales]].min(axis=1)
features['kappa_max'] = features[[f'kappa_{s}' for s in scales]].max(axis=1)
features['kappa_range'] = features['kappa_max'] - features['kappa_min']
features['kappa_gradient'] = features['kappa_10'] - features['kappa_240']  # local - global

# Galactic latitude (proxy for MW foreground)
features['galactic_b'] = np.abs(90.0 - np.degrees(theta[idx_good]))

# Flux ratio (not log)
features['flux_ratio'] = oiii_flux[idx_good] / hbeta_flux[idx_good]

# EW ratio
valid_ew = (oiii_ew[idx_good] > 0) & (hbeta_ew[idx_good] > 0)
features['ew_ratio'] = np.where(valid_ew, oiii_ew[idx_good] / hbeta_ew[idx_good], np.nan)

# Drop NaN damage
features = features.dropna(subset=['damage'])
print(f"  Feature matrix: {len(features)} objects × {len(features.columns)} features")

# ============================================================
# 3. TAIL LABELING
# ============================================================
print("\n[3] Labeling tails...")

# Within each z-bin, label top 5% and top 1% as "tail"
features['z_bin'] = pd.cut(features['z'], bins=10, labels=False)

for pct, label in [(5, 'tail_5pct'), (1, 'tail_1pct')]:
    features[label] = False
    for b in features['z_bin'].unique():
        mask_b = features['z_bin'] == b
        if mask_b.sum() < 20:
            continue
        threshold = np.percentile(features.loc[mask_b, 'damage'], 100 - pct)
        features.loc[mask_b & (features['damage'] >= threshold), label] = True

print(f"  Top 5% tail: {features['tail_5pct'].sum()} objects ({100*features['tail_5pct'].mean():.1f}%)")
print(f"  Top 1% tail: {features['tail_1pct'].sum()} objects ({100*features['tail_1pct'].mean():.1f}%)")

# ============================================================
# 4. TAIL PREDICTION (Logistic Regression)
# ============================================================
print("\n" + "=" * 70)
print("TAIL-CONDITIONAL CLASSIFICATION")
print("=" * 70)

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_auc_score

# Feature groups
source_features = ['snr', 'logLbol', 'logMbh', 'logLedd', 'oiii_flux', 'hbeta_flux', 
                   'oiii_ew', 'hbeta_ew', 'oiii_fwhm', 'hbeta_fwhm', 'flux_ratio']
path_features = ['kappa_mean', 'kappa_min', 'kappa_max', 'kappa_range', 'kappa_gradient',
                 'kappa_10', 'kappa_30', 'kappa_60', 'kappa_120', 'kappa_240']
nuisance_features = ['ebv', 'galactic_b']
all_features = source_features + path_features + nuisance_features

# Clean: drop rows with NaN/inf in any feature
feat_clean = features.dropna(subset=all_features + ['tail_5pct'])
feat_clean = feat_clean.replace([np.inf, -np.inf], np.nan).dropna(subset=all_features)
print(f"\n  Clean sample: {len(feat_clean)} objects")

results = {}

for tail_label, tail_name in [('tail_5pct', 'Top 5%'), ('tail_1pct', 'Top 1%')]:
    print(f"\n--- Predicting {tail_name} tail ---")
    y = feat_clean[tail_label].astype(int).values
    
    for feat_group_name, feat_cols in [
        ('Source only', source_features),
        ('Path only', path_features),
        ('Nuisance only', nuisance_features),
        ('Source + Path', source_features + path_features),
        ('All features', all_features),
    ]:
        X = feat_clean[feat_cols].values.copy()
        
        # Handle any remaining NaN
        col_means = np.nanmean(X, axis=0)
        for j in range(X.shape[1]):
            nan_mask = np.isnan(X[:, j])
            X[nan_mask, j] = col_means[j]
        
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        clf = LogisticRegression(max_iter=1000, class_weight='balanced', C=0.1)
        
        # Cross-validated AUC
        cv_auc = cross_val_score(clf, X_scaled, y, cv=5, scoring='roc_auc')
        
        # Fit full model for coefficients
        clf.fit(X_scaled, y)
        auc_full = roc_auc_score(y, clf.predict_proba(X_scaled)[:, 1])
        
        print(f"  {feat_group_name:<20} CV-AUC: {np.mean(cv_auc):.3f} ± {np.std(cv_auc):.3f}  (full: {auc_full:.3f})")
        
        results[f'{tail_name}_{feat_group_name}'] = {
            'cv_auc_mean': float(np.mean(cv_auc)),
            'cv_auc_std': float(np.std(cv_auc)),
            'full_auc': float(auc_full),
            'n_features': len(feat_cols),
        }
        
        # Print top features for "All features" model
        if feat_group_name == 'All features':
            coef_abs = np.abs(clf.coef_[0])
            top_idx = np.argsort(coef_abs)[::-1][:10]
            print(f"  Top predictors of {tail_name} tail:")
            for rank, idx in enumerate(top_idx):
                feat_name = all_features[idx]
                coef = clf.coef_[0][idx]
                print(f"    {rank+1}. {feat_name:<20} coef={coef:+.4f}")

# ============================================================
# 5. SOURCE vs PATH: Who predicts tails better?
# ============================================================
print("\n" + "=" * 70)
print("SOURCE vs PATH: WHO PREDICTS THE TAILS?")
print("=" * 70)

for tail_name in ['Top 5%', 'Top 1%']:
    source_auc = results[f'{tail_name}_Source only']['cv_auc_mean']
    path_auc = results[f'{tail_name}_Path only']['cv_auc_mean']
    combined_auc = results[f'{tail_name}_Source + Path']['cv_auc_mean']
    
    print(f"\n  {tail_name}:")
    print(f"    Source alone:  AUC = {source_auc:.3f}")
    print(f"    Path alone:    AUC = {path_auc:.3f}")
    print(f"    Combined:      AUC = {combined_auc:.3f}")
    
    if source_auc > path_auc + 0.02:
        print(f"    >>> SOURCE DOMINATES — tail objects are intrinsically different")
    elif path_auc > source_auc + 0.02:
        print(f"    >>> PATH DOMINATES — propagation/environment drives the tails")
    else:
        print(f"    >>> ROUGHLY EQUAL — both contribute or neither explains much")
    
    if combined_auc > max(source_auc, path_auc) + 0.01:
        print(f"    >>> ADDITIVE — source and path carry independent info about tails")

# ============================================================
# 6. LOW-Z BIMODALITY INVESTIGATION
# ============================================================
print("\n" + "=" * 70)
print("LOW-Z BIMODALITY INVESTIGATION")
print("=" * 70)

# Split low-z sample and see what distinguishes the two modes
low_z = features[(features['z'] >= 0.3) & (features['z'] < 0.55)].copy()
print(f"\n  Low-z sample (0.3 < z < 0.55): {len(low_z)} objects")

# Fit 2-component Gaussian mixture
from sklearn.mixture import GaussianMixture
gm = GaussianMixture(n_components=2, random_state=42)
gm.fit(low_z[['damage']].values)
labels = gm.predict(low_z[['damage']].values)

low_z['mode'] = labels
mode_counts = low_z['mode'].value_counts()
print(f"  Mode 0: {mode_counts.get(0, 0)} objects (mean damage = {low_z[low_z['mode']==0]['damage'].mean():.3f})")
print(f"  Mode 1: {mode_counts.get(1, 0)} objects (mean damage = {low_z[low_z['mode']==1]['damage'].mean():.3f})")

# What distinguishes the modes?
print(f"\n  Feature comparison between modes:")
compare_cols = ['z', 'snr', 'ebv', 'logLbol', 'logMbh', 'logLedd', 
                'oiii_flux', 'hbeta_flux', 'oiii_fwhm', 'hbeta_fwhm',
                'kappa_mean', 'kappa_min', 'kappa_range', 'galactic_b']

print(f"  {'Feature':<20} {'Mode 0 med':>12} {'Mode 1 med':>12} {'KS stat':>10} {'KS p':>12} {'Interpretation'}")
print("  " + "-" * 80)

bimodal_results = []
for col in compare_cols:
    v0 = low_z[low_z['mode'] == 0][col].dropna()
    v1 = low_z[low_z['mode'] == 1][col].dropna()
    if len(v0) < 10 or len(v1) < 10:
        continue
    ks_stat, ks_p = stats.ks_2samp(v0, v1)
    med0, med1 = np.median(v0), np.median(v1)
    
    if ks_p < 0.001:
        interp = "*** HIGHLY DIFFERENT"
    elif ks_p < 0.05:
        interp = "* different"
    else:
        interp = "similar"
    
    print(f"  {col:<20} {med0:>12.3f} {med1:>12.3f} {ks_stat:>10.3f} {ks_p:>12.2e} {interp}")
    bimodal_results.append({
        'feature': col, 'med_0': float(med0), 'med_1': float(med1),
        'ks_stat': float(ks_stat), 'ks_p': float(ks_p),
    })

# Check if it's a PLATE/MJD/FIBER split (pipeline regime)
print(f"\n  Pipeline regime check:")
plate = data['PLATE'][idx_good]
mjd = data['MJD'][idx_good]
fiber = data['FIBERID'][idx_good]

low_z_idx = low_z.index
low_z['plate'] = plate[low_z_idx - idx_good[0]].values if len(idx_good) > 0 else 0
low_z['mjd'] = mjd[low_z_idx - idx_good[0]].values if len(idx_good) > 0 else 0

# Check if modes correlate with BOSS vs SDSS-I/II
# BOSS fibers are 2" vs 3" for legacy
boss_flag = data['IF_BOSS_SDSS'][idx_good]

low_z_boss = boss_flag[low_z.index.values - features.index.values[0]]
n_boss_0 = np.sum((labels == 0) & (low_z_boss == 1))
n_boss_1 = np.sum((labels == 1) & (low_z_boss == 1))
n_sdss_0 = np.sum((labels == 0) & (low_z_boss == 0))
n_sdss_1 = np.sum((labels == 1) & (low_z_boss == 0))
print(f"  Mode 0: BOSS={n_boss_0}, SDSS-legacy={n_sdss_0}")
print(f"  Mode 1: BOSS={n_boss_1}, SDSS-legacy={n_sdss_1}")

# Chi-squared test
contingency = np.array([[n_boss_0, n_sdss_0], [n_boss_1, n_sdss_1]])
if contingency.min() > 0:
    chi2, chi_p, _, _ = stats.chi2_contingency(contingency)
    print(f"  Chi² for BOSS/SDSS vs mode: χ² = {chi2:.1f}, p = {chi_p:.2e}")
    if chi_p < 0.01:
        print("  >>> PIPELINE REGIME SPLIT — bimodality is likely BOSS vs SDSS aperture difference!")
    else:
        print("  >>> Not a simple BOSS/SDSS split")

# ============================================================
# 7. Z-STRATIFIED TAIL ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("Z-STRATIFIED: DO TAIL PREDICTORS CHANGE WITH DISTANCE?")
print("=" * 70)

z_quartiles = np.percentile(feat_clean['z'], [0, 25, 50, 75, 100])
for i in range(4):
    z_lo, z_hi = z_quartiles[i], z_quartiles[i+1]
    in_bin = (feat_clean['z'] >= z_lo) & (feat_clean['z'] < z_hi)
    if i == 3:
        in_bin = (feat_clean['z'] >= z_lo) & (feat_clean['z'] <= z_hi)
    
    sub = feat_clean[in_bin]
    if len(sub) < 100:
        continue
    
    y = sub['tail_5pct'].astype(int).values
    if y.sum() < 5:
        continue
    
    # Source vs Path
    for feat_group, feat_cols in [('Source', source_features), ('Path', path_features)]:
        X = sub[feat_cols].values.copy()
        col_means = np.nanmean(X, axis=0)
        for j in range(X.shape[1]):
            nan_mask = np.isnan(X[:, j])
            X[nan_mask, j] = col_means[j]
        
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        clf = LogisticRegression(max_iter=1000, class_weight='balanced', C=0.1)
        
        try:
            cv_auc = cross_val_score(clf, X_scaled, y, cv=3, scoring='roc_auc')
            auc = np.mean(cv_auc)
        except:
            auc = 0.5
        
        if feat_group == 'Source':
            source_auc_z = auc
        else:
            path_auc_z = auc
    
    print(f"  z=[{z_lo:.2f}, {z_hi:.2f}] N={len(sub):>5}  Source AUC={source_auc_z:.3f}  Path AUC={path_auc_z:.3f}  {'PATH↑' if path_auc_z > source_auc_z + 0.02 else 'SOURCE↑' if source_auc_z > path_auc_z + 0.02 else 'TIED'}")

# ============================================================
# SAVE
# ============================================================
print("\n" + "=" * 70)
print("SAVING")
print("=" * 70)

all_results = {
    'test': 'Tail-Conditional Mapping',
    'n_objects': len(feat_clean),
    'classification': results,
    'bimodal_investigation': bimodal_results,
}

with open(f'{OUTDIR}/tail_conditional_results.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print(f"  Saved to {OUTDIR}/tail_conditional_results.json")
print("\nDone.")
