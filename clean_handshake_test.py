#!/usr/bin/env python3
"""
CLEAN HANDSHAKE TEST (GPT Tests A2 + E1 + D1)
================================================
A2: Susceptibility WITHOUT any OIII/Hβ ingredients (no circularity)
E1: Susceptibility-stratified exposure curves (the one-plot test)
D1: Copula test (joint structure without marginals)
B3: NLR vs BLR spatial extent test
"""

import numpy as np
from astropy.io import fits
from scipy import stats
from scipy.optimize import curve_fit
import json, os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_clean_handshake'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("CLEAN HANDSHAKE TEST")
print("=" * 70)

# ============================================================
# 1. Load
# ============================================================
print("\n[1] Loading...")
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

def fix(arr):
    return arr.astype(arr.dtype.newbyteorder('='))

z = fix(data['Z_DR16Q'])
snr = fix(data['SN_MEDIAN_ALL'])
ebv = fix(data['EBV'])

# Lines
oiii_flux = fix(data['OIII5007'][:, 4])
hbeta_flux = fix(data['HBETA'][:, 4])
oiii_ew = fix(data['OIII5007'][:, 2])
hbeta_ew = fix(data['HBETA'][:, 2])
oiii_fwhm = fix(data['OIII5007'][:, 3])
hbeta_fwhm = fix(data['HBETA'][:, 3])

# CIV, MgII for clean susceptibility
civ_flux = fix(data['CIV'][:, 4])
civ_fwhm = fix(data['CIV'][:, 3])
civ_ew = fix(data['CIV'][:, 2])
mgii_flux = fix(data['MGII'][:, 4])
mgii_fwhm = fix(data['MGII'][:, 3])
mgii_ew = fix(data['MGII'][:, 2])

# Hα for BLR
halpha_flux = fix(data['HALPHA'][:, 4])
halpha_fwhm = fix(data['HALPHA'][:, 3])

# NII, SII for NLR
nii_flux = fix(data['NII6585'][:, 4])
sii_flux = fix(data['SII6718'][:, 4])

# Properties
logLbol = fix(data['LOGLBOL'])
logMbh = fix(data['LOGMBH'])
logLedd = fix(data['LOGLEDD_RATIO'])

print(f"  Total: {len(z)}")

# ============================================================
# 2. Source model (same as before, GBR on source features)
# ============================================================
print("\n[2] Source model for OIII/Hβ residuals...")
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_predict

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
log_ratio = np.log10(oiii_flux[idx] / hbeta_flux[idx])

source_X = np.column_stack([
    snr[idx],
    np.log10(np.maximum(oiii_ew[idx], 1e-3)),
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
print(f"  R² = {r2:.4f}, N = {N}")

# Z-bin normalized damage
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
# 3. TEST A2: CLEAN SUSCEPTIBILITY (no OIII/Hβ ingredients)
# ============================================================
print("\n" + "=" * 70)
print("TEST A2: CLEAN SUSCEPTIBILITY (excludes OIII/Hβ)")
print("=" * 70)

# Build susceptibility from ONLY non-OIII/Hβ properties:
# MBH, Ledd, continuum luminosity, EBV
# If CIV/MgII available, add their FWHM (BLR kinematics)

# S_clean: pure AGN state, no line ratio ingredients
S_clean = np.full(N, np.nan)
valid_clean = (np.isfinite(logMbh[idx]) & np.isfinite(logLedd[idx]) & 
               np.isfinite(logLbol[idx]))

# Add CIV FWHM if available
civ_fwhm_valid = (civ_fwhm[idx] > 0) & np.isfinite(civ_fwhm[idx])
mgii_fwhm_valid = (mgii_fwhm[idx] > 0) & np.isfinite(mgii_fwhm[idx])

# For objects with CIV: use CIV FWHM as BLR complexity
# For objects without: use MgII FWHM
# For objects with neither: use only MBH + Ledd

blr_fwhm = np.full(N, np.nan)
blr_fwhm[civ_fwhm_valid] = np.log10(civ_fwhm[idx][civ_fwhm_valid])
no_civ = ~civ_fwhm_valid & mgii_fwhm_valid
blr_fwhm[no_civ] = np.log10(mgii_fwhm[idx][no_civ])

# Clean susceptibility = MBH + Ledd + BLR FWHM (if available) + Lbol
# Standardized and combined
components = []
names = []

# Always available
components.append(logMbh[idx])
names.append('logMbh')
components.append(logLedd[idx])
names.append('logLedd')
components.append(logLbol[idx])
names.append('logLbol')

# Standardize each and sum (equal weight)
S_clean = np.zeros(N)
n_components = 0
for comp in components:
    valid_c = np.isfinite(comp)
    if np.sum(valid_c) < N * 0.5:
        continue
    comp_std = (comp - np.nanmean(comp)) / np.nanstd(comp)
    comp_std[~valid_c] = 0
    S_clean += comp_std
    n_components += 1

# Add BLR FWHM where available
blr_valid = np.isfinite(blr_fwhm)
if np.sum(blr_valid) > N * 0.3:
    blr_std = np.zeros(N)
    blr_std[blr_valid] = (blr_fwhm[blr_valid] - np.nanmean(blr_fwhm[blr_valid])) / np.nanstd(blr_fwhm[blr_valid])
    S_clean += blr_std
    n_components += 1
    names.append('BLR_FWHM')

S_clean /= n_components
valid_s = np.isfinite(S_clean) & np.isfinite(damage)

print(f"  Clean susceptibility from: {', '.join(names)}")
print(f"  Valid objects: {np.sum(valid_s)}")
print(f"  S_clean range: [{np.min(S_clean[valid_s]):.2f}, {np.max(S_clean[valid_s]):.2f}]")

# ============================================================
# 4. TEST E1: Susceptibility-stratified exposure curves
# ============================================================
print("\n" + "=" * 70)
print("TEST E1: SUSCEPTIBILITY-STRATIFIED EXPOSURE CURVES")
print("=" * 70)

Z_v = Z[valid_s]
D_v = damage[valid_s]
S_v = S_clean[valid_s]

# Decile bins on susceptibility
s_edges = np.percentile(S_v, np.linspace(0, 100, 11))
s_bin = np.clip(np.digitize(S_v, s_edges) - 1, 0, 9)

print(f"\n  {'Decile':<8} {'S range':>20} {'N':>6} {'slope':>10} {'p':>12} {'mean_dam':>10}")
print("  " + "-" * 75)

slopes = []
s_centers = []

for d in range(10):
    in_d = s_bin == d
    if np.sum(in_d) < 50:
        slopes.append(np.nan)
        s_centers.append(np.nan)
        continue
    
    sl, inter, r, p, se = stats.linregress(Z_v[in_d], D_v[in_d])
    s_lo, s_hi = s_edges[d], s_edges[d+1]
    
    slopes.append(sl)
    s_centers.append(0.5*(s_lo+s_hi))
    
    label = "LOW" if d < 3 else "HIGH" if d > 6 else "mid"
    print(f"  D{d+1:<7} [{s_lo:>+8.3f}, {s_hi:>+8.3f}] {np.sum(in_d):>6} {sl:>+10.4f} {p:>12.2e} {np.mean(D_v[in_d]):>10.3f} {label}")

valid_sl = [(i, slopes[i]) for i in range(10) if np.isfinite(slopes[i])]
if len(valid_sl) >= 5:
    rho, p_rho = stats.spearmanr([v[0] for v in valid_sl], [v[1] for v in valid_sl])
    print(f"\n  SLOPE vs DECILE: ρ = {rho:.3f}, p = {p_rho:.4f}")
    print(f"  Slopes: {' → '.join([f'{slopes[i]:+.4f}' for i in range(10) if np.isfinite(slopes[i])])}")
    
    if rho > 0.5 and p_rho < 0.05:
        print("  >>> ✅ CLEAN HANDSHAKE CONFIRMED — curves fan out with susceptibility")
    elif rho < -0.5 and p_rho < 0.05:
        print("  >>> ❌ REVERSED — same as before, low susceptibility shows more z-effect")
    else:
        print(f"  >>> {'⚠️ MARGINAL' if abs(rho) > 0.3 else '⏸️ INCONCLUSIVE'}")

# ============================================================
# 5. THE CROSSING CURVES (binned exposure plot per tercile)
# ============================================================
print("\n" + "=" * 70)
print("CROSSING CURVES: Damage vs z by susceptibility tercile")
print("=" * 70)

s_terc = np.percentile(S_v, [0, 33, 67, 100])
tercile_labels = ['Low S (clean)', 'Mid S (clean)', 'High S (clean)']

for t in range(3):
    in_t = (S_v >= s_terc[t]) & (S_v < s_terc[t+1] if t < 2 else True)
    z_t = Z_v[in_t]
    d_t = D_v[in_t]
    
    # 8 z-bins
    zb_edges = np.percentile(z_t, np.linspace(0, 100, 9))
    z_centers = []
    d_means = []
    d_errs = []
    for b in range(8):
        bm = (z_t >= zb_edges[b]) & (z_t < zb_edges[b+1] if b < 7 else True)
        if np.sum(bm) < 10: continue
        z_centers.append(np.mean(z_t[bm]))
        d_means.append(np.mean(d_t[bm]))
        d_errs.append(np.std(d_t[bm]) / np.sqrt(np.sum(bm)))
    
    sl, inter, r, p, se = stats.linregress(z_centers, d_means)
    print(f"\n  {tercile_labels[t]} (N={np.sum(in_t)}):")
    print(f"    slope = {sl:+.4f}, p = {p:.2e}")
    print(f"    z-binned: {' → '.join([f'{d:.3f}' for d in d_means])}")

# ============================================================
# 6. INTERACTION TEST with clean susceptibility
# ============================================================
print("\n" + "=" * 70)
print("INTERACTION: S_clean × z → residual damage")
print("=" * 70)

exposure = (Z_v - np.mean(Z_v)) / np.std(Z_v)
interact = S_v * exposure

from sklearn.linear_model import LinearRegression
X_int = np.column_stack([S_v, exposure, interact])
reg = LinearRegression().fit(X_int, D_v)

rho_int, p_int = stats.spearmanr(interact, D_v)

print(f"  Susceptibility coef: {reg.coef_[0]:+.4f}")
print(f"  Exposure coef:       {reg.coef_[1]:+.4f}")
print(f"  Interaction coef:    {reg.coef_[2]:+.4f}")
print(f"  R²:                  {reg.score(X_int, D_v):.6f}")
print(f"  Spearman(S×z, damage): ρ = {rho_int:.4f}, p = {p_int:.2e}")

if p_int < 1e-10:
    print("  >>> INTERACTION SURVIVES with clean susceptibility!")
elif p_int < 0.01:
    print("  >>> Significant but weaker")
else:
    print("  >>> NOT SIGNIFICANT with clean susceptibility")

# ============================================================
# 7. TEST D1: COPULA TEST (rank-space dependence)
# ============================================================
print("\n" + "=" * 70)
print("TEST D1: COPULA — Joint structure without marginals")
print("=" * 70)

# Convert OIII and Hβ to ranks within z-bins, then measure rank correlation
# If rank correlation collapses while marginals (rank distributions) stay uniform,
# that's a pure joint-structure operator

z_terciles_cop = np.percentile(Z, [0, 33, 67, 100])
cop_labels = ['Low-z', 'Mid-z', 'High-z']

print(f"\n  {'z-bin':<10} {'N':>6} {'Pearson(flux)':>14} {'Spearman(flux)':>15} {'Kendall(flux)':>14} {'Copula MI':>10}")
print("  " + "-" * 75)

for t in range(3):
    in_t = (Z >= z_terciles_cop[t]) & (Z < z_terciles_cop[t+1] if t < 2 else True)
    o = oiii_flux[idx][in_t]
    h = hbeta_flux[idx][in_t]
    
    valid_oh = (o > 0) & (h > 0) & np.isfinite(o) & np.isfinite(h)
    o_v = o[valid_oh]
    h_v = h[valid_oh]
    
    if len(o_v) < 50: continue
    
    # Raw correlations
    pearson_r, _ = stats.pearsonr(np.log10(o_v), np.log10(h_v))
    spearman_r, _ = stats.spearmanr(o_v, h_v)
    kendall_r, _ = stats.kendalltau(o_v, h_v)
    
    # Rank transform (copula)
    o_rank = stats.rankdata(o_v) / (len(o_v) + 1)
    h_rank = stats.rankdata(h_v) / (len(h_v) + 1)
    
    # Copula MI approximation: use rank correlation as proxy
    # True MI would need KDE but rank Spearman captures monotone dependence
    cop_spearman, _ = stats.spearmanr(o_rank, h_rank)
    
    print(f"  {cop_labels[t]:<10} {len(o_v):>6} {pearson_r:>14.4f} {spearman_r:>15.4f} {kendall_r:>14.4f} {cop_spearman:>10.4f}")

# Now do the RESIDUAL copula (after source model subtraction)
print(f"\n  Residual copula (source-subtracted):")
print(f"  {'z-bin':<10} {'N':>6} {'Spearman(resid,z)':>18} {'p':>12}")
print("  " + "-" * 50)

for t in range(3):
    in_t = (Z >= z_terciles_cop[t]) & (Z < z_terciles_cop[t+1] if t < 2 else True)
    r_t = residual[in_t]
    z_t = Z[in_t]
    valid_r = np.isfinite(r_t)
    
    # Rank-transform residual
    r_rank = stats.rankdata(r_t[valid_r]) / (np.sum(valid_r) + 1)
    z_rank = stats.rankdata(z_t[valid_r]) / (np.sum(valid_r) + 1)
    
    rho_cop, p_cop = stats.spearmanr(r_rank, z_rank)
    print(f"  {cop_labels[t]:<10} {np.sum(valid_r):>6} {rho_cop:>18.4f} {p_cop:>12.2e}")

# ============================================================
# 8. TEST B3: NLR vs BLR — Spatial extent test
# ============================================================
print("\n" + "=" * 70)
print("TEST B3: NLR vs BLR SPATIAL EXTENT")
print("=" * 70)

# NLR: [OIII] (narrow, spatially extended ~kpc)
# BLR: Hβ broad, MgII, CIV (compact, ~light-days)
# If handshake is about spatial extent, BLR should be more susceptible

# Compare damage in OIII-only vs Hβ-broad-only residuals
# We already have OIII/Hβ. Let's also check:
# - Hβ_broad vs OIII (BLR/NLR ratio)
# - CIV vs MgII (both BLR but different ionization)

hbeta_br_flux = fix(data['HBETA_BR'][:, 4])
hbeta_br_fwhm = fix(data['HBETA_BR'][:, 3])

# Objects with both NLR and BLR measurements
nlr_blr_good = (good & 
                (hbeta_br_flux > 0) & np.isfinite(hbeta_br_flux) &
                (oiii_flux > 0) & np.isfinite(oiii_flux))
nlr_blr_idx = np.where(nlr_blr_good)[0]

if len(nlr_blr_idx) > 500:
    Z_nb = z[nlr_blr_idx]
    
    # NLR damage: OIII narrow component scatter
    oiii_nb = oiii_flux[nlr_blr_idx]
    hb_br_nb = hbeta_br_flux[nlr_blr_idx]
    
    # BLR/NLR ratio
    blr_nlr_ratio = np.log10(hb_br_nb / oiii_nb)
    valid_bn = np.isfinite(blr_nlr_ratio)
    
    if np.sum(valid_bn) > 200:
        # Z-detrend
        zbn = np.digitize(Z_nb[valid_bn], np.percentile(Z_nb[valid_bn], np.linspace(0, 100, 11))) - 1
        zbn = np.clip(zbn, 0, 9)
        bn_damage = np.full(np.sum(valid_bn), np.nan)
        for b in range(10):
            m = zbn == b
            if np.sum(m) < 5: continue
            med = np.median(blr_nlr_ratio[valid_bn][m])
            mad = max(1.4826 * np.median(np.abs(blr_nlr_ratio[valid_bn][m] - med)), 1e-6)
            bn_damage[m] = np.abs((blr_nlr_ratio[valid_bn][m] - med) / mad)
        
        bn_v = np.isfinite(bn_damage)
        
        # Does BLR/NLR ratio damage correlate with z?
        rho_bn, p_bn = stats.spearmanr(Z_nb[valid_bn][bn_v], bn_damage[bn_v])
        print(f"  BLR/NLR sample: {np.sum(bn_v)} objects")
        print(f"  Spearman(z, BLR/NLR damage): ρ = {rho_bn:.4f}, p = {p_bn:.2e}")
        
        # Interaction with clean susceptibility
        S_nb = S_clean[np.isin(idx, nlr_blr_idx[valid_bn][bn_v])]
        if len(S_nb) == np.sum(bn_v):
            z_nb_std = (Z_nb[valid_bn][bn_v] - np.mean(Z_nb[valid_bn][bn_v])) / np.std(Z_nb[valid_bn][bn_v])
            int_nb = S_nb * z_nb_std
            rho_int_nb, p_int_nb = stats.spearmanr(int_nb, bn_damage[bn_v])
            print(f"  Interaction(S×z) vs BLR/NLR damage: ρ = {rho_int_nb:.4f}, p = {p_int_nb:.2e}")
        
        # Compare: NLR-only damage vs BLR-only damage
        # OIII scatter (NLR)
        oiii_resid = np.full(np.sum(valid_bn), np.nan)
        oiii_log = np.log10(oiii_nb[valid_bn])
        for b in range(10):
            m = zbn == b
            if np.sum(m) < 5: continue
            med = np.median(oiii_log[m])
            mad = max(1.4826 * np.median(np.abs(oiii_log[m] - med)), 1e-6)
            oiii_resid[m] = np.abs((oiii_log[m] - med) / mad)
        
        # Hβ broad scatter (BLR)
        hb_resid = np.full(np.sum(valid_bn), np.nan)
        hb_log = np.log10(hb_br_nb[valid_bn])
        for b in range(10):
            m = zbn == b
            if np.sum(m) < 5: continue
            med = np.median(hb_log[m])
            mad = max(1.4826 * np.median(np.abs(hb_log[m] - med)), 1e-6)
            hb_resid[m] = np.abs((hb_log[m] - med) / mad)
        
        ov = np.isfinite(oiii_resid) & np.isfinite(hb_resid)
        if np.sum(ov) > 100:
            rho_nlr, p_nlr = stats.spearmanr(Z_nb[valid_bn][ov], oiii_resid[ov])
            rho_blr, p_blr = stats.spearmanr(Z_nb[valid_bn][ov], hb_resid[ov])
            print(f"\n  NLR ([OIII]) scatter vs z: ρ = {rho_nlr:.4f}, p = {p_nlr:.2e}")
            print(f"  BLR (Hβ broad) scatter vs z: ρ = {rho_blr:.4f}, p = {p_blr:.2e}")
            if abs(rho_blr) > abs(rho_nlr) + 0.02:
                print("  >>> BLR more affected — spatial extent/geometry matters")
            elif abs(rho_nlr) > abs(rho_blr) + 0.02:
                print("  >>> NLR more affected — NOT microlensing/size")
            else:
                print("  >>> Similar — no spatial extent preference")
else:
    print(f"  Too few BLR objects ({len(nlr_blr_idx)})")

# ============================================================
# 9. CIV vs MgII (both BLR, different ionization)
# ============================================================
print("\n" + "=" * 70)
print("CIV vs MgII: Same region, different ionization")
print("=" * 70)

civ_mgii_good = (good & 
                 (civ_flux > 0) & np.isfinite(civ_flux) &
                 (mgii_flux > 0) & np.isfinite(mgii_flux))
cm_idx = np.where(civ_mgii_good)[0]

if len(cm_idx) > 500:
    Z_cm = z[cm_idx]
    cm_ratio = np.log10(civ_flux[cm_idx] / mgii_flux[cm_idx])
    valid_cm = np.isfinite(cm_ratio)
    
    if np.sum(valid_cm) > 200:
        # Z-detrend
        zcm = np.digitize(Z_cm[valid_cm], np.percentile(Z_cm[valid_cm], np.linspace(0, 100, 11))) - 1
        zcm = np.clip(zcm, 0, 9)
        cm_damage = np.full(np.sum(valid_cm), np.nan)
        for b in range(10):
            m = zcm == b
            if np.sum(m) < 5: continue
            med = np.median(cm_ratio[valid_cm][m])
            mad = max(1.4826 * np.median(np.abs(cm_ratio[valid_cm][m] - med)), 1e-6)
            cm_damage[m] = np.abs((cm_ratio[valid_cm][m] - med) / mad)
        
        cm_v = np.isfinite(cm_damage)
        rho_cm, p_cm = stats.spearmanr(Z_cm[valid_cm][cm_v], cm_damage[cm_v])
        print(f"  CIV/MgII sample: {np.sum(cm_v)} objects")
        print(f"  Spearman(z, CIV/MgII damage): ρ = {rho_cm:.4f}, p = {p_cm:.2e}")
else:
    print(f"  Too few CIV+MgII objects ({len(cm_idx)})")

# ============================================================
# SAVE
# ============================================================
print("\n\nDone.")
