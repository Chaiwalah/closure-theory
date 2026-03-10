#!/usr/bin/env python3
"""
KILL C1 — Pipeline Artifact Test
==================================
C1 says: SDSS pipeline measures EW worse at high-z (SNR, continuum placement,
template drift) while FWHM is more robust. The asymmetry is measurement, not physics.

TO KILL C1 we need independent evidence. What we CAN do right now:

1. INTERNAL CONSISTENCY: If it's pipeline, the effect should correlate with
   pipeline quality indicators (SNR, chi2, plate quality) AFTER controlling for z.
   If it doesn't correlate with pipeline quality, C1 is wounded.

2. ERROR BAR TEST: SDSS provides EW_err and FWHM_err. If C1 is right,
   the EW degradation should vanish when we restrict to low-error measurements.

3. MULTI-EPOCH TEST: Some quasars were observed multiple times by SDSS.
   If C1 is right, the "degradation" should be INCONSISTENT between epochs.
   If A1 is right, it should be CONSISTENT (same sightline = same effect).

4. BOSS vs eBOSS: DR16Q contains quasars from both surveys with different
   spectrographs. If the same quasars show the same pattern, C1 is weakened.

5. WAVELENGTH COVERAGE TEST: At different z, different lines fall in different
   parts of the CCD with different sensitivity. If C1 is right, the "degradation"
   should depend on observed wavelength. If A1, it should depend on REST wavelength
   (through P).
"""

import numpy as np
from scipy import stats
from astropy.io import fits
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = Path('results_kill_C1')
RESULTS_DIR.mkdir(exist_ok=True)

f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z = d['Z_DR16Q']
snr = d['SN_MEDIAN_ALL']

# Check what quality columns exist
all_cols = [c.name for c in f[1].columns]
print("Available quality columns:")
quality_cols = [c for c in all_cols if any(k in c.upper() for k in ['CHI', 'QUAL', 'SNR', 'SN_', 'PLATE', 'MJD', 'FIBER', 'INSTRUMENT', 'SURVEY', 'BOSS'])]
for c in sorted(quality_cols)[:30]:
    print(f"  {c}")

lines = {}
for name in ['OIII5007', 'HBETA', 'CIV', 'CIII_ALL', 'MGII', 'LYA',
             'NII6585', 'SII6718', 'OII3728', 'HALPHA']:
    col = d[name]
    if col.ndim == 2 and col.shape[1] >= 6:
        lines[name] = {'ew': col[:, 2], 'ew_err': col[:, 3], 
                       'fwhm': col[:, 4], 'fwhm_err': col[:, 5]}

# Try to get survey/instrument info
survey = None
for cname in ['SURVEY', 'INSTRUMENT', 'SOURCE_Z', 'SPECTRO']:
    if cname in all_cols:
        survey = d[cname]
        print(f"\nFound survey column: {cname}")
        if hasattr(survey, 'dtype') and survey.dtype.kind in ('U', 'S', 'O'):
            unique_vals = np.unique(survey[:100])
            print(f"  Sample values: {unique_vals[:10]}")
        break

plate = d['PLATE'] if 'PLATE' in all_cols else None
mjd = d['MJD'] if 'MJD' in all_cols else None

f.close()

print(f"\n{'='*80}")
print("🔬 KILL C1 — PIPELINE ARTIFACT TESTS")
print(f"{'='*80}")

# =====================================================================
# TEST 1: ERROR BAR CONTROL
# =====================================================================
print(f"\n{'='*80}")
print("TEST 1: ERROR BAR CONTROL")
print("If C1 is right, restricting to low EW_err should eliminate the effect.")
print("If A1 is right, low-error measurements should show the same pattern.")
print("=" * 80)

test_lines = [
    ('CIV', 'CIII_ALL', 'CIV-CIII'),
    ('OIII5007', 'HBETA', 'OIII-Hβ'),
]

for name_a, name_b, label in test_lines:
    if name_a not in lines or name_b not in lines:
        continue
    
    ew_a, ew_b = lines[name_a]['ew'], lines[name_b]['ew']
    ew_err_a, ew_err_b = lines[name_a]['ew_err'], lines[name_b]['ew_err']
    fwhm_a, fwhm_b = lines[name_a]['fwhm'], lines[name_b]['fwhm']
    
    print(f"\n  {label}:")
    
    # Relative error
    rel_err_a = np.abs(ew_err_a / (ew_a + 1e-10))
    rel_err_b = np.abs(ew_err_b / (ew_b + 1e-10))
    
    for err_label, err_cut in [("ALL data", 999), ("EW err < 20%", 0.20), ("EW err < 10%", 0.10), ("EW err < 5%", 0.05)]:
        err_mask = (rel_err_a < err_cut) & (rel_err_b < err_cut)
        
        z_bins = [(0.3, 0.7), (0.7, 1.0), (1.0, 1.5), (1.5, 2.5)]
        print(f"\n    {err_label}:")
        print(f"    {'z-bin':<10} {'r(EW)':<10} {'r(FWHM)':<10} {'N':<8}")
        
        for z_lo, z_hi in z_bins:
            mask = err_mask & (z > z_lo) & (z < z_hi)
            valid_ew = mask & np.isfinite(ew_a) & np.isfinite(ew_b) & (ew_a != 0) & (ew_b != 0)
            valid_fw = mask & np.isfinite(fwhm_a) & np.isfinite(fwhm_b) & (fwhm_a > 0) & (fwhm_b > 0)
            
            if valid_ew.sum() < 30:
                print(f"    {z_lo:.1f}-{z_hi:.1f}    (too few: {valid_ew.sum()})")
                continue
            
            r_ew, _ = stats.pearsonr(ew_a[valid_ew], ew_b[valid_ew])
            r_fw = np.nan
            if valid_fw.sum() >= 30:
                r_fw, _ = stats.pearsonr(fwhm_a[valid_fw], fwhm_b[valid_fw])
            
            fw_str = f"{r_fw:+.3f}" if np.isfinite(r_fw) else "N/A"
            print(f"    {z_lo:.1f}-{z_hi:.1f}    {r_ew:+.4f}   {fw_str}   {valid_ew.sum()}")


# =====================================================================
# TEST 2: SNR QUARTILE INDEPENDENCE
# =====================================================================
print(f"\n\n{'='*80}")
print("TEST 2: SNR QUARTILE INDEPENDENCE")
print("If C1 is right, the degradation should be STRONGEST in low-SNR data.")
print("If A1 is right, it should be present at ALL SNR levels.")
print("=" * 80)

snr_quartiles = np.nanpercentile(snr[np.isfinite(snr)], [25, 50, 75])
snr_labels = [
    (0, snr_quartiles[0], "Q1 (lowest SNR)"),
    (snr_quartiles[0], snr_quartiles[1], "Q2"),
    (snr_quartiles[1], snr_quartiles[2], "Q3"),
    (snr_quartiles[2], 9999, "Q4 (highest SNR)"),
]

for name_a, name_b, pair_label in test_lines:
    if name_a not in lines or name_b not in lines:
        continue
    
    ew_a, ew_b = lines[name_a]['ew'], lines[name_b]['ew']
    fwhm_a, fwhm_b = lines[name_a]['fwhm'], lines[name_b]['fwhm']
    
    print(f"\n  {pair_label}:")
    print(f"  {'SNR quartile':<20} {'z<0.8 r(EW)':<14} {'z>1.0 r(EW)':<14} {'Δr(EW)':<10} {'z>1.0 r(FW)':<14}")
    
    for snr_lo, snr_hi, snr_label in snr_labels:
        snr_mask = (snr > snr_lo) & (snr < snr_hi)
        
        # Low z
        m_low = snr_mask & (z > 0.3) & (z < 0.8) & np.isfinite(ew_a) & np.isfinite(ew_b) & (ew_a != 0) & (ew_b != 0)
        # High z
        m_high = snr_mask & (z > 1.0) & (z < 2.5) & np.isfinite(ew_a) & np.isfinite(ew_b) & (ew_a != 0) & (ew_b != 0)
        m_high_f = snr_mask & (z > 1.0) & (z < 2.5) & np.isfinite(fwhm_a) & np.isfinite(fwhm_b) & (fwhm_a > 0) & (fwhm_b > 0)
        
        r_low = r_high = r_fw_high = np.nan
        if m_low.sum() >= 30:
            r_low, _ = stats.pearsonr(ew_a[m_low], ew_b[m_low])
        if m_high.sum() >= 30:
            r_high, _ = stats.pearsonr(ew_a[m_high], ew_b[m_high])
        if m_high_f.sum() >= 30:
            r_fw_high, _ = stats.pearsonr(fwhm_a[m_high_f], fwhm_b[m_high_f])
        
        delta = r_high - r_low if np.isfinite(r_high) and np.isfinite(r_low) else np.nan
        
        lo_str = f"{r_low:+.3f}" if np.isfinite(r_low) else "N/A"
        hi_str = f"{r_high:+.3f}" if np.isfinite(r_high) else "N/A"
        d_str = f"{delta:+.4f}" if np.isfinite(delta) else "N/A"
        fw_str = f"{r_fw_high:+.3f}" if np.isfinite(r_fw_high) else "N/A"
        
        print(f"  {snr_label:<20} {lo_str:<14} {hi_str:<14} {d_str:<10} {fw_str:<14}")


# =====================================================================
# TEST 3: OBSERVED vs REST WAVELENGTH
# =====================================================================
print(f"\n\n{'='*80}")
print("TEST 3: OBSERVED vs REST WAVELENGTH DEPENDENCE")
print("C1 predicts: degradation depends on where line falls on CCD (observed λ)")
print("A1 predicts: degradation depends on line identity (rest λ / susceptibility P)")
print("=" * 80)
print("""
Key test: CIV at z=1.5 and MgII at z=0.0 both appear at ~3900Å on CCD.
If C1 is right, they should show SIMILAR degradation (same detector region).
If A1 is right, CIV should be MORE degraded (higher P) regardless of where it falls.
""")

# CIV rest=1549, at z=1.5 → observed 3872Å
# MgII rest=2798, at z=0.4 → observed 3917Å  
# Very similar observed wavelength!

if 'CIV' in lines and 'MGII' in lines:
    ew_civ = lines['CIV']['ew']
    fwhm_civ = lines['CIV']['fwhm']
    ew_mgii = lines['MGII']['ew']
    fwhm_mgii = lines['MGII']['fwhm']
    
    # CIV at z~1.5 (observed ~3900Å)
    civ_mask = (z > 1.3) & (z < 1.7) & np.isfinite(ew_civ) & (ew_civ != 0) & np.isfinite(fwhm_civ) & (fwhm_civ > 0)
    # MgII at z~0.4 (observed ~3900Å)  
    mgii_mask = (z > 0.3) & (z < 0.5) & np.isfinite(ew_mgii) & (ew_mgii != 0) & np.isfinite(fwhm_mgii) & (fwhm_mgii > 0)
    
    print(f"  CIV at z=1.3-1.7 (observed ~3870-4180Å): N={civ_mask.sum()}")
    print(f"  MgII at z=0.3-0.5 (observed ~3640-4200Å): N={mgii_mask.sum()}")
    
    # Measure EW scatter (CV) as degradation proxy
    civ_ew_cv = np.std(ew_civ[civ_mask]) / np.abs(np.mean(ew_civ[civ_mask]))
    mgii_ew_cv = np.std(ew_mgii[mgii_mask]) / np.abs(np.mean(ew_mgii[mgii_mask]))
    
    civ_fw_cv = np.std(fwhm_civ[civ_mask]) / np.mean(fwhm_civ[civ_mask])
    mgii_fw_cv = np.std(fwhm_mgii[mgii_mask]) / np.mean(fwhm_mgii[mgii_mask])
    
    print(f"\n  Same observed wavelength region (~3900Å):")
    print(f"  CIV (z~1.5, rest 1549Å): EW CV={civ_ew_cv:.3f}, FWHM CV={civ_fw_cv:.3f}, ratio={civ_ew_cv/civ_fw_cv:.2f}")
    print(f"  MgII (z~0.4, rest 2798Å): EW CV={mgii_ew_cv:.3f}, FWHM CV={mgii_fw_cv:.3f}, ratio={mgii_ew_cv/mgii_fw_cv:.2f}")
    
    if civ_ew_cv / civ_fw_cv > mgii_ew_cv / mgii_fw_cv:
        print(f"\n  CIV is MORE degraded (EW/FWHM ratio higher) despite same CCD region")
        print(f"  → Degradation follows LINE IDENTITY (P), not detector position")
        print(f"  → C1 (pipeline artifact) CANNOT explain this")
    else:
        print(f"\n  MgII equally or more degraded — C1 still viable")

    # Now compare CIV at DIFFERENT observed wavelengths but same line
    print(f"\n  CIV at different observed λ (same line, different z):")
    for z_lo, z_hi in [(1.0, 1.3), (1.3, 1.7), (1.7, 2.2), (2.2, 3.0)]:
        m = (z > z_lo) & (z < z_hi) & np.isfinite(ew_civ) & (ew_civ != 0) & np.isfinite(fwhm_civ) & (fwhm_civ > 0)
        if m.sum() < 100:
            continue
        obs_lam = 1549 * (1 + (z_lo + z_hi) / 2)
        cv_ew = np.std(ew_civ[m]) / np.abs(np.mean(ew_civ[m]))
        cv_fw = np.std(fwhm_civ[m]) / np.mean(fwhm_civ[m])
        print(f"    z={z_lo:.1f}-{z_hi:.1f} (obs ~{obs_lam:.0f}Å): EW CV={cv_ew:.3f}, FWHM CV={cv_fw:.3f}, ratio={cv_ew/cv_fw:.2f}, N={m.sum()}")


# =====================================================================
# TEST 4: PLATE-LEVEL CONSISTENCY  
# =====================================================================
print(f"\n\n{'='*80}")
print("TEST 4: PLATE-LEVEL CONSISTENCY")
print("If C1 is right, the effect should vary between plates (instrument/conditions).")
print("If A1 is right, it should be consistent (same sky = same effect).")
print("=" * 80)

if plate is not None and 'CIV' in lines:
    ew_civ = lines['CIV']['ew']
    fwhm_civ = lines['CIV']['fwhm']
    
    # Get plates with enough high-z quasars
    above_z = (z > 1.5) & (z < 2.5) & np.isfinite(ew_civ) & (ew_civ != 0) & np.isfinite(fwhm_civ) & (fwhm_civ > 0)
    
    unique_plates = np.unique(plate[above_z])
    plate_stats = []
    
    for p in unique_plates:
        m = above_z & (plate == p)
        if m.sum() < 50:
            continue
        
        cv_ew = np.std(ew_civ[m]) / np.abs(np.mean(ew_civ[m]))
        cv_fw = np.std(fwhm_civ[m]) / np.mean(fwhm_civ[m])
        plate_stats.append({'plate': int(p), 'n': int(m.sum()), 'cv_ew': cv_ew, 'cv_fw': cv_fw, 'ratio': cv_ew / cv_fw})
    
    if len(plate_stats) > 10:
        ratios = [p['ratio'] for p in plate_stats]
        print(f"  Plates with 50+ high-z CIV quasars: {len(plate_stats)}")
        print(f"  EW/FWHM CV ratio across plates:")
        print(f"    Mean: {np.mean(ratios):.3f}")
        print(f"    Std:  {np.std(ratios):.3f}")
        print(f"    CV of ratio: {np.std(ratios)/np.mean(ratios):.3f}")
        print(f"    Min: {np.min(ratios):.3f}, Max: {np.max(ratios):.3f}")
        
        if np.std(ratios) / np.mean(ratios) < 0.3:
            print(f"\n  → EW/FWHM degradation ratio is CONSISTENT across plates")
            print(f"  → C1 weakened: pipeline conditions vary by plate but effect doesn't")
        else:
            print(f"\n  → Significant plate-to-plate variation — C1 still viable")


print(f"\n\n{'='*80}")
print("SUMMARY — C1 KILL STATUS")
print("=" * 80)
print("""
C1 (pipeline artifact) is killed if:
1. Effect persists in low-error subsamples ✓/✗
2. Effect is SNR-independent ✓/✗  
3. Lines at same CCD position but different identity show different degradation ✓/✗
4. Effect is plate-consistent ✓/✗

The DEFINITIVE C1 kill requires an independent pipeline (DESI line catalog).
These tests can wound C1 but not fully eliminate it.
""")

with open(RESULTS_DIR / 'column_list.txt', 'w') as fout:
    for c in sorted(all_cols):
        fout.write(c + '\n')
print(f"Column list saved to {RESULTS_DIR}/column_list.txt")
