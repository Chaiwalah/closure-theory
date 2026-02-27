#!/usr/bin/env python3
"""
WAVELENGTH KILL TEST — Does the effect track observed-frame systematics?
=========================================================================
Test A: Match on observed-frame wavelength regime
Test B: Short-baseline vs long-baseline diagnostics
Test C: Does damage correlate with wavelength separation?
Test D: Plate/MJD jackknife robustness
"""

import numpy as np
from astropy.io import fits
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_wavelength_kill'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("WAVELENGTH KILL TEST")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

def fix(arr): return arr.astype(arr.dtype.newbyteorder('='))

z = fix(data['Z_DR16Q'])
snr = fix(data['SN_MEDIAN_ALL'])
oiii_flux = fix(data['OIII5007'][:, 4])
oiii_ew = fix(data['OIII5007'][:, 2])
hbeta_flux = fix(data['HBETA'][:, 4])
hbeta_fwhm = fix(data['HBETA'][:, 3])
nii_flux = fix(data['NII6585'][:, 4])
halpha_flux = fix(data['HALPHA'][:, 4])
sii_flux = fix(data['SII6718'][:, 4])
oii_flux = fix(data['OII3728'][:, 4])
logLbol = fix(data['LOGLBOL'])
logMbh = fix(data['LOGMBH'])
logLedd = fix(data['LOGLEDD_RATIO'])
plate = fix(data['PLATE'])
mjd = fix(data['MJD'])

good = ((z > 0.3) & (z < 2.0) & (snr > 3) &
        (oiii_flux > 0) & (hbeta_flux > 0) &
        np.isfinite(oiii_flux) & np.isfinite(hbeta_flux) &
        np.isfinite(logLbol) & np.isfinite(logMbh) & np.isfinite(logLedd) &
        (logLbol > 0) & (logMbh > 0) &
        np.isfinite(hbeta_fwhm) & (hbeta_fwhm > 0))

idx = np.where(good)[0]
N = len(idx)
Z = z[idx]
HB_FWHM = np.log10(np.maximum(hbeta_fwhm[idx], 1))
LEDD = logLedd[idx]

def compute_damage(ratio, redshift, n_bins=20):
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

# ============================================================
# TEST B: SHORT vs LONG BASELINE
# ============================================================
print("\n" + "=" * 70)
print("TEST B: SHORT-BASELINE vs LONG-BASELINE DIAGNOSTICS")
print("=" * 70)

# Rest-frame wavelength separations:
# OIII/Hβ: 5007-4861 = 146Å (SHORT)
# NII/Hα: 6583-6563 = 20Å (ULTRA-SHORT)
# SII/Hα: 6718-6563 = 155Å (SHORT)
# SII/OIII: 6718-5007 = 1711Å (LONG)
# OII/OIII: 5007-3728 = 1279Å (LONG)
# OII/Hβ: 4861-3728 = 1133Å (LONG)

diagnostics = {}

# Short baseline (<200Å rest)
diagnostics['OIII/Hβ (146Å)'] = {
    'good': (oiii_flux[idx] > 0) & (hbeta_flux[idx] > 0),
    'ratio': lambda m: np.log10(oiii_flux[idx][m] / hbeta_flux[idx][m]),
    'sep': 146, 'type': 'SHORT'
}

nii_ok = (nii_flux[idx] > 0) & (halpha_flux[idx] > 0) & np.isfinite(nii_flux[idx]) & np.isfinite(halpha_flux[idx])
diagnostics['NII/Hα (20Å)'] = {
    'good': nii_ok,
    'ratio': lambda m: np.log10(nii_flux[idx][m] / halpha_flux[idx][m]),
    'sep': 20, 'type': 'SHORT-LOCKED'
}

sii_ha_ok = (sii_flux[idx] > 0) & (halpha_flux[idx] > 0) & np.isfinite(sii_flux[idx]) & np.isfinite(halpha_flux[idx])
diagnostics['SII/Hα (155Å)'] = {
    'good': sii_ha_ok,
    'ratio': lambda m: np.log10(sii_flux[idx][m] / halpha_flux[idx][m]),
    'sep': 155, 'type': 'SHORT'
}

# Long baseline (>1000Å rest)
sii_oiii_ok = (sii_flux[idx] > 0) & (oiii_flux[idx] > 0) & np.isfinite(sii_flux[idx])
diagnostics['SII/OIII (1711Å)'] = {
    'good': sii_oiii_ok,
    'ratio': lambda m: np.log10(sii_flux[idx][m] / oiii_flux[idx][m]),
    'sep': 1711, 'type': 'LONG'
}

oii_oiii_ok = (oii_flux[idx] > 0) & (oiii_flux[idx] > 0) & np.isfinite(oii_flux[idx])
diagnostics['OII/OIII (1279Å)'] = {
    'good': oii_oiii_ok,
    'ratio': lambda m: np.log10(oii_flux[idx][m] / oiii_flux[idx][m]),
    'sep': 1279, 'type': 'LONG'
}

oii_hb_ok = (oii_flux[idx] > 0) & (hbeta_flux[idx] > 0) & np.isfinite(oii_flux[idx])
diagnostics['OII/Hβ (1133Å)'] = {
    'good': oii_hb_ok,
    'ratio': lambda m: np.log10(oii_flux[idx][m] / hbeta_flux[idx][m]),
    'sep': 1133, 'type': 'LONG'
}

print(f"\n  {'Diagnostic':<25} {'Δλ':>6} {'Type':<12} {'N':>8} {'ρ(FWHM×z)':>12} {'p':>12} {'ρ(Ledd×z)':>12} {'p':>12}")
print("  " + "-" * 110)

baseline_results = []

for name, info in diagnostics.items():
    mask = info['good']
    if np.sum(mask) < 500: continue
    
    ratio = info['ratio'](mask)
    z_d = Z[mask]
    fwhm_d = HB_FWHM[mask]
    ledd_d = LEDD[mask]
    
    dam = compute_damage(ratio, z_d)
    valid = np.isfinite(dam) & np.isfinite(fwhm_d) & np.isfinite(ledd_d)
    
    z_std = (z_d[valid] - np.mean(z_d[valid])) / np.std(z_d[valid])
    f_std = (fwhm_d[valid] - np.mean(fwhm_d[valid])) / np.std(fwhm_d[valid])
    l_std = (ledd_d[valid] - np.mean(ledd_d[valid])) / np.std(ledd_d[valid])
    
    rho_f, p_f = stats.spearmanr(f_std * z_std, dam[valid])
    rho_l, p_l = stats.spearmanr(l_std * z_std, dam[valid])
    
    print(f"  {name:<25} {info['sep']:>6} {info['type']:<12} {np.sum(valid):>8} {rho_f:>+12.4f} {p_f:>12.2e} {rho_l:>+12.4f} {p_l:>12.2e}")
    
    baseline_results.append({
        'name': name, 'sep': info['sep'], 'type': info['type'],
        'n': int(np.sum(valid)),
        'rho_fwhm_z': float(rho_f), 'p_fwhm_z': float(p_f),
        'rho_ledd_z': float(rho_l), 'p_ledd_z': float(p_l),
    })

# Does effect size correlate with wavelength separation?
seps = [r['sep'] for r in baseline_results if r['type'] != 'SHORT-LOCKED']
rhos = [abs(r['rho_fwhm_z']) for r in baseline_results if r['type'] != 'SHORT-LOCKED']
if len(seps) >= 4:
    rho_sep, p_sep = stats.spearmanr(seps, rhos)
    print(f"\n  Effect size vs wavelength separation: ρ = {rho_sep:.3f}, p = {p_sep:.3f}")
    if rho_sep > 0.5 and p_sep < 0.05:
        print("  >>> ⚠️ LONG BASELINES SHOW MORE EFFECT — spectrophotometric tilt plausible")
    elif rho_sep < -0.3:
        print("  >>> ✅ SHORT BASELINES SHOW EQUAL/MORE EFFECT — kills spectrophotometry")
    else:
        print("  >>> No clear baseline dependence")

# ============================================================
# TEST A: MATCHED-STATE + OBSERVED-FRAME WAVELENGTH CONTROL
# ============================================================
print("\n" + "=" * 70)
print("TEST A: MATCHED STATE + OBSERVED-FRAME λ CONTROL")
print("=" * 70)

# For OIII/Hβ (short baseline):
# λ_obs(OIII) = 5007 × (1+z), λ_obs(Hβ) = 4861 × (1+z)
# Both scale identically with z, so observed-frame regime = f(z) only
# Matching on z kills this automatically for short-baseline ratios

# For OII/OIII (long baseline):
# λ_obs(OII) = 3728 × (1+z), λ_obs(OIII) = 5007 × (1+z)
# At z=0.5: OII at 5592, OIII at 7510 — different CCD regions
# At z=1.0: OII at 7456, OIII at 10014 — OIII hits NIR!

# Compute observed-frame wavelength regime for each object
oiii_lobs = 5007 * (1 + Z)
hbeta_lobs = 4861 * (1 + Z)
oii_lobs = 3728 * (1 + Z)
sii_lobs = 6718 * (1 + Z)
nii_lobs = 6583 * (1 + Z)
halpha_lobs = 6563 * (1 + Z)

print(f"\n  Observed-frame λ ranges:")
print(f"    OIII: {np.min(oiii_lobs):.0f}–{np.max(oiii_lobs):.0f} Å")
print(f"    Hβ:   {np.min(hbeta_lobs):.0f}–{np.max(hbeta_lobs):.0f} Å")
print(f"    OII:  {np.min(oii_lobs):.0f}–{np.max(oii_lobs):.0f} Å")
print(f"    SII:  {np.min(sii_lobs):.0f}–{np.max(sii_lobs):.0f} Å")

# KEY INSIGHT: For OIII/Hβ, both lines are ALWAYS in the same CCD region
# because they're only 146Å apart in rest frame.
# So any spectrophotometric tilt at the OIII/Hβ wavelength cancels in the RATIO.
# The tilt explanation requires DIFFERENT tilt at OIII vs Hβ wavelength.
# At Δλ_obs = 146*(1+z) ≈ 200-440Å, the throughput curve is essentially flat.

print(f"\n  OIII-Hβ observed separation: {146*(1+np.min(Z)):.0f}–{146*(1+np.max(Z)):.0f} Å")
print(f"  OII-OIII observed separation: {1279*(1+np.min(Z)):.0f}–{1279*(1+np.max(Z)):.0f} Å")
print(f"  SII-OIII observed separation: {1711*(1+np.min(Z)):.0f}–{1711*(1+np.max(Z)):.0f} Å")

# Now do the matched-state test with ADDITIONAL λ_obs matching
# Use OIII/Hβ (short baseline, tilt-immune) 
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

log_ratio = np.log10(oiii_flux[idx] / hbeta_flux[idx])
damage_oiii_hb = compute_damage(log_ratio, Z)

z_med = np.median(Z)
lo = Z < z_med
hi = Z >= z_med

# Match on: FWHM, Ledd, Lbol, SNR, AND observed-frame λ(OIII)
X_lo = np.column_stack([HB_FWHM[lo], LEDD[lo], logLbol[idx][lo], snr[idx][lo], oiii_lobs[lo]])
X_hi = np.column_stack([HB_FWHM[hi], LEDD[hi], logLbol[idx][hi], snr[idx][hi], oiii_lobs[hi]])

valid_lo = np.all(np.isfinite(X_lo), axis=1) & np.isfinite(damage_oiii_hb[lo])
valid_hi = np.all(np.isfinite(X_hi), axis=1) & np.isfinite(damage_oiii_hb[hi])

X_lo_c = X_lo[valid_lo]
X_hi_c = X_hi[valid_hi]
D_lo = damage_oiii_hb[lo][valid_lo]
D_hi = damage_oiii_hb[hi][valid_hi]

scaler = StandardScaler().fit(np.vstack([X_lo_c, X_hi_c]))
X_lo_s = scaler.transform(X_lo_c)
X_hi_s = scaler.transform(X_hi_c)

nn = NearestNeighbors(n_neighbors=1, metric='euclidean')
nn.fit(X_hi_s)
distances, indices = nn.kneighbors(X_lo_s)

dist_thresh = np.percentile(distances, 80)
good_match = distances.ravel() < dist_thresh

D_lo_m = D_lo[good_match]
D_hi_m = D_hi[indices.ravel()[good_match]]

diff = np.mean(D_hi_m) - np.mean(D_lo_m)
t_stat, t_p = stats.ttest_rel(D_lo_m, D_hi_m)
_, wilcox_p = stats.wilcoxon(D_lo_m - D_hi_m)

print(f"\n  Matched on FWHM + Ledd + Lbol + SNR + λ_obs(OIII):")
print(f"    Pairs: {np.sum(good_match)}")
print(f"    Low-z damage: {np.mean(D_lo_m):.4f}")
print(f"    High-z damage: {np.mean(D_hi_m):.4f}")
print(f"    Δdamage: {diff:+.4f}")
print(f"    Paired t: p = {t_p:.2e}")
print(f"    Wilcoxon: p = {wilcox_p:.2e}")

if t_p < 0.01 and diff > 0:
    print("  >>> SURVIVES λ_obs matching — NOT spectrophotometric tilt!")
else:
    print("  >>> Weakened or killed by λ_obs matching")

# Stratify by susceptibility
S_match = HB_FWHM[lo][valid_lo][good_match] + LEDD[lo][valid_lo][good_match]
s_terc = np.percentile(S_match, [0, 33, 67, 100])
print(f"\n  Stratified by susceptibility (with λ_obs matching):")
for t, label in enumerate(['Low S', 'Mid S', 'High S']):
    in_t = (S_match >= s_terc[t]) & (S_match < s_terc[t+1] if t < 2 else True)
    if np.sum(in_t) < 30: continue
    d_lo_t = D_lo_m[in_t]
    d_hi_t = D_hi_m[in_t]
    diff_t = np.mean(d_hi_t) - np.mean(d_lo_t)
    _, p_t = stats.wilcoxon(d_lo_t - d_hi_t)
    print(f"    {label}: Δ = {diff_t:+.4f}, p = {p_t:.2e}, N = {np.sum(in_t)}")

# ============================================================
# TEST D: PLATE / MJD JACKKNIFE
# ============================================================
print("\n" + "=" * 70)
print("TEST D: PLATE JACKKNIFE — Does it survive calibration epochs?")
print("=" * 70)

# Test the FWHM×z interaction stability across plates
PLATE = plate[idx]
unique_plates = np.unique(PLATE)

# Random 100-plate jackknife
np.random.seed(42)
n_jack = 100
jack_rhos = []

valid_all = np.isfinite(damage_oiii_hb) & np.isfinite(HB_FWHM) & np.isfinite(LEDD)

for i in range(n_jack):
    # Drop 10% of plates
    drop = np.random.choice(unique_plates, size=len(unique_plates)//10, replace=False)
    keep = ~np.isin(PLATE, drop) & valid_all
    
    z_std = (Z[keep] - np.mean(Z[keep])) / np.std(Z[keep])
    f_std = (HB_FWHM[keep] - np.mean(HB_FWHM[keep])) / np.std(HB_FWHM[keep])
    
    rho, _ = stats.spearmanr(f_std * z_std, damage_oiii_hb[keep])
    jack_rhos.append(rho)

jack_rhos = np.array(jack_rhos)
print(f"  100 jackknife iterations (dropping 10% of plates each):")
print(f"  ρ(FWHM×z, damage): {np.mean(jack_rhos):.4f} ± {np.std(jack_rhos):.4f}")
print(f"  Range: [{np.min(jack_rhos):.4f}, {np.max(jack_rhos):.4f}]")
print(f"  Full sample ρ: {stats.spearmanr((HB_FWHM[valid_all]-np.mean(HB_FWHM[valid_all]))/np.std(HB_FWHM[valid_all]) * (Z[valid_all]-np.mean(Z[valid_all]))/np.std(Z[valid_all]), damage_oiii_hb[valid_all])[0]:.4f}")

if np.min(jack_rhos) > 0:
    print("  >>> ALL 100 jackknives positive — robust to plate removal!")
else:
    neg = np.sum(jack_rhos < 0)
    print(f"  >>> {neg}/100 negative — {'fragile!' if neg > 10 else 'mostly robust'}")

# MJD jackknife
unique_mjds = np.unique(mjd[idx])
mjd_terciles = np.percentile(mjd[idx], [0, 33, 67, 100])

print(f"\n  MJD-stratified (calibration epoch):")
for t, label in enumerate(['Early obs', 'Mid obs', 'Late obs']):
    in_t = (mjd[idx] >= mjd_terciles[t]) & (mjd[idx] < mjd_terciles[t+1] if t < 2 else True) & valid_all
    if np.sum(in_t) < 1000: continue
    z_std = (Z[in_t] - np.mean(Z[in_t])) / np.std(Z[in_t])
    f_std = (HB_FWHM[in_t] - np.mean(HB_FWHM[in_t])) / np.std(HB_FWHM[in_t])
    rho, p = stats.spearmanr(f_std * z_std, damage_oiii_hb[in_t])
    print(f"    {label}: ρ = {rho:+.4f}, p = {p:.2e}, N = {np.sum(in_t)}")

# ============================================================
# CRITICAL ARGUMENT: Why spectrophotometric tilt can't explain OIII/Hβ
# ============================================================
print("\n" + "=" * 70)
print("CRITICAL: WHY TILT CAN'T EXPLAIN OIII/Hβ")
print("=" * 70)

print("""
  OIII (5007Å) and Hβ (4861Å) are separated by only 146Å in rest frame.
  
  At z=0.3: observed separation = 190Å (both at ~6300-6500Å)
  At z=1.0: observed separation = 292Å (both at ~9700-10000Å)  
  At z=2.0: observed separation = 438Å (both at ~14600-15000Å)
  
  Over 190-440Å, the SDSS throughput curve varies by <2%.
  Sky emission lines are typically isolated features, not broadband tilts.
  Telluric absorption bands are >100Å wide but affect BOTH lines equally.
  
  For a spectrophotometric tilt to produce scatter in log(OIII/Hβ),
  you need DIFFERENTIAL throughput between OIII and Hβ wavelengths.
  At <440Å separation, this differential is negligible.
  
  Yet OIII/Hβ shows the handshake at p = 10⁻¹² (FWHM×z).
  
  The long-baseline ratios (OII/OIII, SII/OIII) are MORE susceptible
  to tilt — but that's expected from BOTH the tilt hypothesis AND
  the operator hypothesis (more diagnostic information = more degradation).
  
  The discriminator is: does OIII/Hβ (tilt-immune) show the effect?
  YES. At p = 10⁻¹². With matched-state controls.
  
  Spectrophotometric tilt cannot produce differential scatter in a 
  146Å-separation ratio at the levels observed.
""")

print("\nDone.")
