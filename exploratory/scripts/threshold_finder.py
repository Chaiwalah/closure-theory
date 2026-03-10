#!/usr/bin/env python3
"""
THRESHOLD FINDER — Let the data speak.
No forced sigmoids. No preconceptions.
Just: how does EVERYTHING change with z?
Look for breakpoints, regime changes, gradients.
"""

import numpy as np
from astropy.io import fits
from scipy.stats import spearmanr, pearsonr, ks_2samp
import json
import os
import warnings
warnings.filterwarnings('ignore')

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data
z = data['Z_DR16Q']

lines = [
    ('LYA', 'LYA', 1216),
    ('CIV', 'CIV', 1549),
    ('CIII', 'CIII_ALL', 1909),
    ('MGII', 'MGII', 2798),
    ('OII', 'OII3728', 3728),
    ('HBETA', 'HBETA', 4861),
    ('OIII', 'OIII5007', 5007),
    ('HALPHA', 'HALPHA', 6563),
    ('NII', 'NII6585', 6585),
    ('SII', 'SII6718', 6718),
]

def get_ew_fwhm(col_name):
    try:
        col = data[col_name]
        if len(col.shape) < 2 or col.shape[1] < 5:
            return None, None, None
        ew = col[:, 2]
        fwhm = col[:, 4]
        valid = (np.isfinite(ew) & np.isfinite(fwhm) & 
                 (ew > 0) & (fwhm > 0) & (ew < 500) & (fwhm < 30000) & (z > 0.05))
        return ew, fwhm, valid
    except:
        return None, None, None

# ============================================================
# METHOD 1: SLIDING WINDOW BREAKPOINT DETECTION
# For each line, slide a window through z and compute 
# EW/FWHM ratio in each window. Where does the ratio change?
# ============================================================
print("=" * 70)
print("METHOD 1: EW/FWHM RATIO — SLIDING WINDOW")
print("Simple ratio. If medium eats EW but not FWHM, ratio drops with z.")
print("=" * 70)

for name, col_name, rest_wav in lines:
    ew, fwhm, valid = get_ew_fwhm(col_name)
    if ew is None or valid.sum() < 3000:
        continue
    
    z_v = z[valid]
    ratio = ew[valid] / fwhm[valid]
    
    # Sort by z
    order = np.argsort(z_v)
    z_s = z_v[order]
    r_s = ratio[order]
    
    # Sliding window: 5000 objects wide, step 1000
    window = min(5000, len(z_s) // 5)
    step = window // 5
    
    z_med = []
    ratio_med = []
    ratio_iqr = []
    
    for start in range(0, len(z_s) - window, step):
        chunk = r_s[start:start+window]
        z_med.append(np.median(z_s[start:start+window]))
        ratio_med.append(np.median(chunk))
        ratio_iqr.append(np.percentile(chunk, 75) - np.percentile(chunk, 25))
    
    z_med = np.array(z_med)
    ratio_med = np.array(ratio_med)
    ratio_iqr = np.array(ratio_iqr)
    
    # Find maximum gradient (biggest change)
    if len(ratio_med) > 3:
        gradient = np.gradient(ratio_med, z_med)
        max_drop_idx = np.argmin(gradient)
        max_rise_idx = np.argmax(gradient)
        
        # Overall trend
        r_trend, p_trend = spearmanr(z_med, ratio_med)
        
        print(f"\n{name} ({rest_wav}Å) N={valid.sum():,} | z=[{z_s[0]:.2f}–{z_s[-1]:.2f}]")
        print(f"  Median EW/FWHM: {ratio_med[0]:.4f} → {ratio_med[-1]:.4f} (Δ={ratio_med[-1]-ratio_med[0]:+.4f})")
        print(f"  Trend: r={r_trend:+.3f} p={p_trend:.2e}")
        print(f"  Steepest drop at z={z_med[max_drop_idx]:.3f} (gradient={gradient[max_drop_idx]:.4f})")
        print(f"  Steepest rise at z={z_med[max_rise_idx]:.3f} (gradient={gradient[max_rise_idx]:.4f})")

# ============================================================
# METHOD 2: SPLIT-HALF KS TEST
# For every possible z-split, compare low-z vs high-z 
# distributions. Where is the split most different?
# ============================================================
print("\n" + "=" * 70)
print("METHOD 2: OPTIMAL SPLIT — Where do two halves differ most?")
print("KS test on EW distribution, scanning z splits.")
print("=" * 70)

for name, col_name, rest_wav in lines:
    ew, fwhm, valid = get_ew_fwhm(col_name)
    if ew is None or valid.sum() < 5000:
        continue
    
    z_v = z[valid]
    ew_v = ew[valid]
    fwhm_v = fwhm[valid]
    
    # Scan z-splits
    z_splits = np.linspace(np.percentile(z_v, 15), np.percentile(z_v, 85), 30)
    
    ks_ew = []
    ks_fwhm = []
    
    for z_cut in z_splits:
        lo = z_v < z_cut
        hi = z_v >= z_cut
        if lo.sum() < 500 or hi.sum() < 500:
            ks_ew.append(0)
            ks_fwhm.append(0)
            continue
        
        stat_ew, _ = ks_2samp(ew_v[lo], ew_v[hi])
        stat_fwhm, _ = ks_2samp(fwhm_v[lo], fwhm_v[hi])
        ks_ew.append(stat_ew)
        ks_fwhm.append(stat_fwhm)
    
    ks_ew = np.array(ks_ew)
    ks_fwhm = np.array(ks_fwhm)
    
    # The DIFFERENCE between EW divergence and FWHM divergence
    # If medium eats EW not FWHM: EW KS should peak but FWHM KS shouldn't
    ks_diff = ks_ew - ks_fwhm
    
    best_ew = z_splits[np.argmax(ks_ew)]
    best_fwhm = z_splits[np.argmax(ks_fwhm)]
    best_diff = z_splits[np.argmax(ks_diff)]
    
    print(f"\n{name} ({rest_wav}Å)")
    print(f"  Best EW split:       z = {best_ew:.3f} (KS = {ks_ew.max():.4f})")
    print(f"  Best FWHM split:     z = {best_fwhm:.3f} (KS = {ks_fwhm.max():.4f})")
    print(f"  Best EW-FWHM gap:    z = {best_diff:.3f} (gap = {ks_diff.max():.4f})")
    print(f"  EW diverges more?    {'YES' if ks_ew.max() > ks_fwhm.max() else 'NO'} ({ks_ew.max()/max(ks_fwhm.max(),0.001):.2f}×)")

# ============================================================
# METHOD 3: RUNNING KURTOSIS
# Kurtosis measures tail heaviness. If medium adds 
# stochastic damage, kurtosis should increase.
# ============================================================
print("\n" + "=" * 70)
print("METHOD 3: RUNNING KURTOSIS OF EW")
print("Stochastic damage → heavier tails → kurtosis rises")
print("=" * 70)

from scipy.stats import kurtosis

for name, col_name, rest_wav in lines:
    ew, fwhm, valid = get_ew_fwhm(col_name)
    if ew is None or valid.sum() < 5000:
        continue
    
    z_v = z[valid]
    ew_v = ew[valid]
    fwhm_v = fwhm[valid]
    
    order = np.argsort(z_v)
    z_s = z_v[order]
    ew_s = ew_v[order]
    fwhm_s = fwhm_v[order]
    
    window = min(5000, len(z_s) // 5)
    step = window // 5
    
    z_med = []
    kurt_ew = []
    kurt_fwhm = []
    
    for start in range(0, len(z_s) - window, step):
        z_med.append(np.median(z_s[start:start+window]))
        kurt_ew.append(kurtosis(ew_s[start:start+window]))
        kurt_fwhm.append(kurtosis(fwhm_s[start:start+window]))
    
    z_med = np.array(z_med)
    kurt_ew = np.array(kurt_ew)
    kurt_fwhm = np.array(kurt_fwhm)
    
    r_ew, p_ew = spearmanr(z_med, kurt_ew)
    r_fwhm, p_fwhm = spearmanr(z_med, kurt_fwhm)
    
    # Where does kurtosis peak?
    peak_ew = z_med[np.argmax(kurt_ew)]
    peak_fwhm = z_med[np.argmax(kurt_fwhm)]
    
    print(f"\n{name} ({rest_wav}Å)")
    print(f"  EW kurtosis trend:   r={r_ew:+.3f} p={p_ew:.2e} | peak at z={peak_ew:.3f}")
    print(f"  FWHM kurtosis trend: r={r_fwhm:+.3f} p={p_fwhm:.2e} | peak at z={peak_fwhm:.3f}")
    
    if r_ew > 0 and p_ew < 0.05 and (r_fwhm <= 0 or p_fwhm > 0.05):
        print(f"  → EW tails fatten with z, FWHM doesn't. SELECTIVE DAMAGE ✓")

# ============================================================
# METHOD 4: CROSS-LINE COHERENCE
# Do different lines' EWs stay correlated with each other?
# If medium damages them, inter-line coherence should decay.
# ============================================================
print("\n" + "=" * 70)
print("METHOD 4: INTER-LINE EW COHERENCE vs z")
print("Do different lines' EWs track each other? Does that break?")
print("=" * 70)

# Pick pairs of lines that coexist in z-range
pairs = [
    ('MGII', 'CIII_ALL', 'MGII-CIII'),
    ('CIV', 'CIII_ALL', 'CIV-CIII'),
    ('HBETA', 'OIII5007', 'HBETA-OIII'),
    ('MGII', 'CIV', 'MGII-CIV'),
    ('CIV', 'LYA', 'CIV-LYA'),
]

for col1, col2, pair_name in pairs:
    try:
        c1 = data[col1]
        c2 = data[col2]
    except:
        continue
    
    if len(c1.shape) < 2 or len(c2.shape) < 2:
        continue
    
    ew1 = c1[:, 2]
    ew2 = c2[:, 2]
    
    valid = (np.isfinite(ew1) & np.isfinite(ew2) & 
             (ew1 > 0) & (ew2 > 0) & (ew1 < 500) & (ew2 < 500) & (z > 0.05))
    
    if valid.sum() < 3000:
        continue
    
    z_v = z[valid]
    ew1_v = ew1[valid]
    ew2_v = ew2[valid]
    
    # Bin by z
    z_edges = np.linspace(np.percentile(z_v, 5), np.percentile(z_v, 95), 15)
    
    corrs = []
    z_mids = []
    
    for i in range(len(z_edges)-1):
        mask = (z_v >= z_edges[i]) & (z_v < z_edges[i+1])
        if mask.sum() < 100:
            continue
        r, p = spearmanr(ew1_v[mask], ew2_v[mask])
        corrs.append(r)
        z_mids.append((z_edges[i] + z_edges[i+1]) / 2)
    
    if len(z_mids) < 5:
        continue
    
    corrs = np.array(corrs)
    z_mids = np.array(z_mids)
    
    r_trend, p_trend = spearmanr(z_mids, corrs)
    
    print(f"\n{pair_name}: N={valid.sum():,}")
    print(f"  Coherence: {corrs[0]:.3f} → {corrs[-1]:.3f} (Δ={corrs[-1]-corrs[0]:+.3f})")
    print(f"  Trend: r={r_trend:+.3f} p={p_trend:.2e}")
    
    if r_trend < -0.3 and p_trend < 0.05:
        # Find where coherence drops fastest
        grad = np.gradient(corrs, z_mids)
        steepest = z_mids[np.argmin(grad)]
        print(f"  → COHERENCE DECAYING. Steepest drop at z={steepest:.3f}")

# ============================================================
# METHOD 5: RAW NUMBERS — What does EACH z-bin look like?
# Just print the data. Look with human eyes.
# ============================================================
print("\n" + "=" * 70)
print("METHOD 5: RAW DATA TABLE — MGII (biggest sample)")
print("Just look at the numbers.")
print("=" * 70)

ew, fwhm, valid = get_ew_fwhm('MGII')
z_v = z[valid]
ew_v = ew[valid]
fwhm_v = fwhm[valid]

z_edges = np.linspace(0.5, 2.5, 21)
print(f"\n{'z_bin':>8} {'N':>7} {'med_EW':>8} {'med_FWHM':>9} {'EW/FWHM':>8} {'std_EW':>8} {'std_FWHM':>9} {'EW_kurt':>8} {'FWHM_kurt':>10}")
print("-" * 88)

for i in range(len(z_edges)-1):
    mask = (z_v >= z_edges[i]) & (z_v < z_edges[i+1])
    if mask.sum() < 100:
        continue
    
    ew_bin = ew_v[mask]
    fwhm_bin = fwhm_v[mask]
    
    z_mid = (z_edges[i] + z_edges[i+1]) / 2
    print(f"{z_mid:8.2f} {mask.sum():7d} {np.median(ew_bin):8.2f} {np.median(fwhm_bin):9.1f} "
          f"{np.median(ew_bin)/np.median(fwhm_bin):8.4f} {np.std(ew_bin):8.2f} {np.std(fwhm_bin):9.1f} "
          f"{kurtosis(ew_bin):8.2f} {kurtosis(fwhm_bin):10.2f}")

print("\n" + "=" * 70)
print("DONE. Now look at the data with human eyes.")
print("=" * 70)

hdu.close()
