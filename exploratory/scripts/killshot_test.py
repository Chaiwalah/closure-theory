#!/usr/bin/env python3
"""
KILLSHOT TEST — Gemini's adversarial recommendation
====================================================
Test: Does the variance of branch-locked line ratios increase with z?

We don't have the OIII] 1660/1666 doublet, but we CAN test:
1. OIII 5007 EW/FWHM ratio variance vs z (locked line — should be FLAT)
2. NII 6583 EW/FWHM ratio variance vs z (locked line — should be FLAT)  
3. SII 6716 EW/FWHM ratio variance vs z (unlocked, density-diagnostic — should RISE)
4. CIV EW/FWHM ratio variance vs z (unlocked, resonance — should RISE)

Also test: EW coefficient of variation (CV = std/mean) vs z
for locked vs unlocked lines.

If locked lines have FLAT variance and unlocked lines have RISING variance,
that's the branching-ratio invariance test adapted to available data.
"""

import numpy as np
from astropy.io import fits
from scipy.stats import spearmanr, kurtosis
import warnings
warnings.filterwarnings('ignore')

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data
z = data['Z_DR16Q']

lines = {
    # LOCKED (branch-ratio immune, P≈0)
    'OIII_5007': ('OIII5007', 'LOCKED'),
    'NII_6583':  ('NII6585',  'LOCKED'),
    # UNLOCKED (susceptible)
    'SII_6716':  ('SII6718',  'UNLOCKED'),
    'CIV_1549':  ('CIV',      'UNLOCKED'),
    'MGII_2798': ('MGII',     'UNLOCKED'),
    'LYA_1216':  ('LYA',      'UNLOCKED'),
    'CIII_1909': ('CIII_ALL',  'UNLOCKED'),
    'HBETA':     ('HBETA',     'UNLOCKED'),
}

print("=" * 80)
print("KILLSHOT TEST: BRANCH-LOCKED vs UNLOCKED LINE VARIANCE EVOLUTION")
print("Locked lines: variance should be FLAT with z (immune to medium)")
print("Unlocked lines: variance should RISE with z (medium scrambles them)")
print("=" * 80)

results = {}

for name, (col, lock_status) in lines.items():
    try:
        c = data[col]
    except:
        continue
    if len(c.shape) < 2 or c.shape[1] < 5:
        continue
    
    ew = c[:, 2]
    fwhm = c[:, 4]
    
    valid = (np.isfinite(ew) & np.isfinite(fwhm) & 
             (ew > 0) & (fwhm > 0) & (ew < 500) & (fwhm < 30000) & (z > 0.05))
    
    if valid.sum() < 3000:
        continue
    
    zv = z[valid]
    ewv = ew[valid]
    fwv = fwhm[valid]
    
    # Equal-count bins
    idx = np.argsort(zv)
    n = len(idx)
    nbins = min(15, n // 500)
    bs = n // nbins
    
    zmeds = []
    ew_cv = []      # coefficient of variation
    ew_kurt = []    # kurtosis
    fwhm_cv = []
    fwhm_kurt = []
    ew_iqr = []     # interquartile range / median (robust CV)
    
    for i in range(nbins):
        sl = idx[i*bs:(i+1)*bs]
        zmeds.append(np.median(zv[sl]))
        
        ew_bin = ewv[sl]
        fw_bin = fwv[sl]
        
        ew_cv.append(np.std(ew_bin) / np.mean(ew_bin))
        fwhm_cv.append(np.std(fw_bin) / np.mean(fw_bin))
        ew_kurt.append(kurtosis(ew_bin))
        fwhm_kurt.append(kurtosis(fw_bin))
        ew_iqr.append((np.percentile(ew_bin, 75) - np.percentile(ew_bin, 25)) / np.median(ew_bin))
    
    zmeds = np.array(zmeds)
    ew_cv = np.array(ew_cv)
    ew_kurt = np.array(ew_kurt)
    fwhm_cv = np.array(fwhm_cv)
    fwhm_kurt = np.array(fwhm_kurt)
    ew_iqr = np.array(ew_iqr)
    
    # Trends
    r_cv, p_cv = spearmanr(zmeds, ew_cv)
    r_kurt, p_kurt = spearmanr(zmeds, ew_kurt)
    r_iqr, p_iqr = spearmanr(zmeds, ew_iqr)
    r_fwhm_cv, p_fwhm_cv = spearmanr(zmeds, fwhm_cv)
    
    results[name] = {
        'lock': lock_status,
        'N': valid.sum(),
        'z_range': (float(zmeds[0]), float(zmeds[-1])),
        'ew_cv_trend_r': float(r_cv),
        'ew_cv_trend_p': float(p_cv),
        'ew_kurt_trend_r': float(r_kurt),
        'ew_kurt_trend_p': float(p_kurt),
        'ew_iqr_trend_r': float(r_iqr),
        'ew_iqr_trend_p': float(p_iqr),
        'fwhm_cv_trend_r': float(r_fwhm_cv),
    }
    
    rises = "RISES" if r_cv > 0.3 and p_cv < 0.05 else "FLAT" if abs(r_cv) < 0.3 else "FALLS"
    kurt_rises = "RISES" if r_kurt > 0.3 and p_kurt < 0.05 else "FLAT" if abs(r_kurt) < 0.3 else "FALLS"
    
    expected = "should be FLAT" if lock_status == "LOCKED" else "should RISE"
    match = "✓" if (lock_status == "LOCKED" and rises != "RISES") or \
                    (lock_status == "UNLOCKED" and rises == "RISES") else \
             "?" if rises == "FLAT" else "✗"
    
    print(f"\n{name} [{lock_status}] — N={valid.sum():,} z=[{zmeds[0]:.2f}–{zmeds[-1]:.2f}]")
    print(f"  EW CV trend:   r={r_cv:+.3f} p={p_cv:.2e} → {rises} ({expected}) {match}")
    print(f"  EW kurtosis:   r={r_kurt:+.3f} p={p_kurt:.2e} → {kurt_rises}")
    print(f"  EW IQR trend:  r={r_iqr:+.3f} p={p_iqr:.2e}")
    print(f"  FWHM CV trend: r={r_fwhm_cv:+.3f} p={p_fwhm_cv:.2e}")

# ============================================================
# SUMMARY: LOCKED vs UNLOCKED
# ============================================================
print("\n" + "=" * 80)
print("SUMMARY: LOCKED vs UNLOCKED COMPARISON")
print("=" * 80)

locked_cvs = [r['ew_cv_trend_r'] for n, r in results.items() if r['lock'] == 'LOCKED']
unlocked_cvs = [r['ew_cv_trend_r'] for n, r in results.items() if r['lock'] == 'UNLOCKED']
locked_kurts = [r['ew_kurt_trend_r'] for n, r in results.items() if r['lock'] == 'LOCKED']
unlocked_kurts = [r['ew_kurt_trend_r'] for n, r in results.items() if r['lock'] == 'UNLOCKED']

print(f"\n  LOCKED lines (should be flat):")
print(f"    Mean EW CV trend:      r = {np.mean(locked_cvs):+.3f}")
print(f"    Mean EW kurtosis trend: r = {np.mean(locked_kurts):+.3f}")

print(f"\n  UNLOCKED lines (should rise):")
print(f"    Mean EW CV trend:      r = {np.mean(unlocked_cvs):+.3f}")
print(f"    Mean EW kurtosis trend: r = {np.mean(unlocked_kurts):+.3f}")

gap_cv = np.mean(unlocked_cvs) - np.mean(locked_cvs)
gap_kurt = np.mean(unlocked_kurts) - np.mean(locked_kurts)

print(f"\n  GAP (unlocked - locked):")
print(f"    CV trend gap:      {gap_cv:+.3f}")
print(f"    Kurtosis trend gap: {gap_kurt:+.3f}")

if gap_kurt > 0.3:
    print(f"\n  🎯 UNLOCKED lines show MORE kurtosis growth than LOCKED lines")
    print(f"  This is consistent with refractive decoherence prediction.")
    print(f"  Branching-locked lines are immune. Diagnostic lines are not.")
elif gap_kurt > 0:
    print(f"\n  ⚠️  Trend in predicted direction but gap is small.")
else:
    print(f"\n  ✗ No clear separation. Needs investigation.")

# ============================================================
# THE TABLE Gemini asked for
# ============================================================
print(f"\n{'='*80}")
print(f"GEMINI'S VERDICT TABLE — Adapted")
print(f"{'='*80}")
print(f"\n{'Line':>12} {'Lock':>8} {'EW_CV':>8} {'EW_k':>8} {'FWHM_CV':>8} {'Prediction':>12} {'Result':>8}")
print("-" * 70)
for name in ['OIII_5007','NII_6583','SII_6716','CIV_1549','MGII_2798','LYA_1216','CIII_1909','HBETA']:
    if name not in results: continue
    r = results[name]
    pred = "FLAT" if r['lock'] == 'LOCKED' else "RISE"
    cv_dir = "↑" if r['ew_cv_trend_r'] > 0.3 else "→" if abs(r['ew_cv_trend_r']) < 0.3 else "↓"
    k_dir = "↑" if r['ew_kurt_trend_r'] > 0.3 else "→" if abs(r['ew_kurt_trend_r']) < 0.3 else "↓"
    
    print(f"{name:>12} {r['lock']:>8} {cv_dir:>8} {k_dir:>8} "
          f"{'→' if abs(r['fwhm_cv_trend_r'])<0.3 else '↑':>8} {pred:>12}")

print(f"\n{'='*80}")
print("DONE — Killshot test complete.")
print(f"{'='*80}")

hdu.close()
