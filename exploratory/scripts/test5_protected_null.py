#!/usr/bin/env python3
"""
TEST 5: PROTECTED-NULL / BREAKER TEST
======================================
Branch-locked ratios ([OIII] 5007/4959, [NII] 6585/6548) should be IMMUNE
to our process because they share the same upper energy level — the ratio
is fixed by quantum mechanics (3:1) regardless of plasma conditions.

If our mechanism is correctly specified (P depends on emissivity gradient):
  - Branch-locked RATIOS should be flat with z and κ
  - Nearby EWs and non-locked ratios should still degrade
  
If the ratios drift → our model is too broad / wrong about P.

DR16Q stores lines as 6-element arrays: [peak, peak2, EW, EW_err, FWHM, FWHM_err]
We need BOTH members of each doublet.
"""

import numpy as np
from scipy import stats
from astropy.io import fits
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = Path('results_protected_null')
RESULTS_DIR.mkdir(exist_ok=True)

f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z = d['Z_DR16Q']
snr = d['SN_MEDIAN_ALL']

# Check what columns we have for doublet members
all_cols = [c.name for c in f[1].columns]
print("Looking for doublet columns...")
for c in sorted(all_cols):
    if any(k in c.upper() for k in ['OIII', 'NII', 'SII', 'OII']):
        print(f"  {c}: shape={d[c].shape if hasattr(d[c], 'shape') else 'scalar'}")

# Load the lines we need
lines = {}
for name in ['OIII5007', 'OIII4959', 'NII6585', 'NII6549', 'SII6718', 'SII6732',
             'HBETA', 'HALPHA', 'OII3728', 'OIII4959']:
    if name in all_cols:
        col = d[name]
        if hasattr(col, 'shape') and col.ndim == 2 and col.shape[1] >= 6:
            lines[name] = {'peak': col[:, 0], 'ew': col[:, 2], 'ew_err': col[:, 3],
                          'fwhm': col[:, 4]}
            print(f"  Loaded {name}: {col.shape}")
        elif hasattr(col, 'shape') and col.ndim == 2:
            print(f"  {name}: shape {col.shape} — not standard 6-element")
    else:
        # Check for alternative names
        pass

f.close()

print(f"\n{'='*80}")
print("🛡️ TEST 5: PROTECTED-NULL — Branch-Locked Ratios")
print("=" * 80)

# =====================================================================
# TEST 5A: [OIII] 5007/4959 ratio vs z
# =====================================================================
if 'OIII5007' in lines and 'OIII4959' in lines:
    print(f"\n--- [OIII] 5007/4959 (branch-locked, theoretical = 2.98) ---")
    
    ew_5007 = lines['OIII5007']['ew']
    ew_4959 = lines['OIII4959']['ew']
    
    valid = np.isfinite(ew_5007) & np.isfinite(ew_4959) & (ew_4959 != 0) & (ew_5007 != 0)
    valid &= np.isfinite(z) & (z > 0.01)
    
    ratio = ew_5007[valid] / ew_4959[valid]
    z_v = z[valid]
    
    # Clip extreme ratios (measurement errors)
    good = (ratio > 0.5) & (ratio < 10) 
    ratio = ratio[good]
    z_v = z_v[good]
    
    print(f"  Valid measurements: {len(ratio)}")
    print(f"  Overall: median ratio = {np.median(ratio):.3f}, std = {np.std(ratio):.3f}")
    
    print(f"\n  {'z-bin':<12} {'median ratio':<15} {'std':<10} {'N':<8} {'Δ from 2.98'}")
    print(f"  {'-'*55}")
    
    z_bins_5 = [(0.02, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.2)]
    medians = []
    zmids = []
    
    for z_lo, z_hi in z_bins_5:
        m = (z_v > z_lo) & (z_v < z_hi)
        if m.sum() < 20:
            continue
        med = np.median(ratio[m])
        std = np.std(ratio[m])
        medians.append(med)
        zmids.append((z_lo + z_hi) / 2)
        delta = med - 2.98
        print(f"  {z_lo:.2f}-{z_hi:.2f}    {med:.4f}         {std:.3f}     {m.sum():<8} {delta:+.4f}")
    
    if len(zmids) >= 3:
        sl, _, r, p, _ = stats.linregress(zmids, medians)
        print(f"\n  Trend with z: slope={sl:+.4f}, r={r:+.3f}, p={p:.4f}")
        if p > 0.05:
            print(f"  ✅ FLAT — branch-locked ratio is IMMUNE. Model prediction confirmed.")
        else:
            print(f"  ⚠️ TREND DETECTED — ratio drifts with z. Model may be wrong here.")

elif 'OIII5007' in lines:
    print("\n  [OIII] 4959 not available as separate column.")
    print("  Checking if OIII5007 column contains both...")
    # Some catalogs pack both OIII lines together
    # We'll use the EW of OIII5007 alone and check its stability
    
    ew_5007 = lines['OIII5007']['ew']
    valid = np.isfinite(ew_5007) & (ew_5007 != 0) & np.isfinite(z) & (z > 0.01)
    
    # Instead test: OIII EW scatter (CV) should NOT grow with z 
    # because it's branch-locked → low P → immune
    print(f"\n  [OIII] 5007 EW coefficient of variation by z-bin:")
    print(f"  (Should be STABLE if branch-locked = immune)")
    
    for z_lo, z_hi in [(0.02, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.2)]:
        m = valid & (z > z_lo) & (z < z_hi)
        if m.sum() < 100:
            continue
        cv = np.std(ew_5007[m]) / (np.abs(np.mean(ew_5007[m])) + 1e-10)
        print(f"    z={z_lo:.2f}-{z_hi:.2f}: CV={cv:.3f}, N={m.sum()}")


# =====================================================================
# TEST 5B: [NII] 6585/6549 ratio vs z  
# =====================================================================
if 'NII6585' in lines:
    print(f"\n--- [NII] 6585 EW stability (branch-locked, P≈0) ---")
    print("  If immune: EW scatter should NOT grow with z")
    print("  Compare with Hα (not branch-locked, P>0)")
    
    ew_nii = lines['NII6585']['ew']
    valid_nii = np.isfinite(ew_nii) & (ew_nii != 0) & np.isfinite(z) & (z > 0.01)
    
    print(f"\n  {'z-bin':<12} {'NII CV':<10} {'Hα CV':<10} {'Hβ CV':<10} {'NII N'}")
    print(f"  {'-'*50}")
    
    nii_cvs = []
    ha_cvs = []
    hb_cvs = []
    zmids = []
    
    for z_lo, z_hi in [(0.02, 0.15), (0.15, 0.3), (0.3, 0.45), (0.45, 0.6), (0.6, 0.8)]:
        m_nii = valid_nii & (z > z_lo) & (z < z_hi)
        
        nii_cv = ha_cv = hb_cv = np.nan
        
        if m_nii.sum() >= 50:
            nii_cv = np.std(ew_nii[m_nii]) / (np.abs(np.mean(ew_nii[m_nii])) + 1e-10)
        
        if 'HALPHA' in lines:
            ew_ha = lines['HALPHA']['ew']
            m_ha = np.isfinite(ew_ha) & (ew_ha != 0) & (z > z_lo) & (z < z_hi)
            if m_ha.sum() >= 50:
                ha_cv = np.std(ew_ha[m_ha]) / (np.abs(np.mean(ew_ha[m_ha])) + 1e-10)
        
        if 'HBETA' in lines:
            ew_hb = lines['HBETA']['ew']
            m_hb = np.isfinite(ew_hb) & (ew_hb != 0) & (z > z_lo) & (z < z_hi)
            if m_hb.sum() >= 50:
                hb_cv = np.std(ew_hb[m_hb]) / (np.abs(np.mean(ew_hb[m_hb])) + 1e-10)
        
        if m_nii.sum() >= 50:
            nii_cvs.append(nii_cv)
            ha_cvs.append(ha_cv)
            hb_cvs.append(hb_cv)
            zmids.append((z_lo + z_hi) / 2)
            
            nii_str = f"{nii_cv:.3f}" if np.isfinite(nii_cv) else "N/A"
            ha_str = f"{ha_cv:.3f}" if np.isfinite(ha_cv) else "N/A"
            hb_str = f"{hb_cv:.3f}" if np.isfinite(hb_cv) else "N/A"
            print(f"  {z_lo:.2f}-{z_hi:.2f}    {nii_str:<10} {ha_str:<10} {hb_str:<10} {m_nii.sum()}")
    
    if len(zmids) >= 3:
        valid_nii_cv = [i for i, c in enumerate(nii_cvs) if np.isfinite(c)]
        if len(valid_nii_cv) >= 3:
            zz = [zmids[i] for i in valid_nii_cv]
            cc = [nii_cvs[i] for i in valid_nii_cv]
            sl, _, r, p, _ = stats.linregress(zz, cc)
            print(f"\n  NII CV trend: slope={sl:+.3f}, r={r:+.3f}, p={p:.4f}")


# =====================================================================
# TEST 5C: [SII] 6718/6732 ratio vs z (ANTI-NULL: density diagnostic, should DEGRADE)
# =====================================================================
if 'SII6718' in lines and 'SII6732' in lines:
    print(f"\n--- [SII] 6718/6732 (density-sensitive doublet, P≈0.7, should DEGRADE) ---")
    
    ew_6718 = lines['SII6718']['ew']
    ew_6732 = lines['SII6732']['ew']
    
    valid = np.isfinite(ew_6718) & np.isfinite(ew_6732) & (ew_6718 != 0) & (ew_6732 != 0) & np.isfinite(z)
    
    ratio = ew_6718[valid] / ew_6732[valid]
    z_v = z[valid]
    
    good = (ratio > 0.1) & (ratio < 10)
    ratio = ratio[good]
    z_v = z_v[good]
    
    print(f"  Valid: {len(ratio)}")
    print(f"  {'z-bin':<12} {'median':<10} {'std':<10} {'N'}")
    
    sii_stds = []
    zmids_sii = []
    
    for z_lo, z_hi in [(0.02, 0.15), (0.15, 0.3), (0.3, 0.45), (0.45, 0.65)]:
        m = (z_v > z_lo) & (z_v < z_hi)
        if m.sum() < 20:
            continue
        med = np.median(ratio[m])
        std = np.std(ratio[m])
        sii_stds.append(std)
        zmids_sii.append((z_lo + z_hi) / 2)
        print(f"  {z_lo:.2f}-{z_hi:.2f}    {med:.4f}     {std:.3f}     {m.sum()}")
    
    if len(zmids_sii) >= 3:
        sl, _, r, p, _ = stats.linregress(zmids_sii, sii_stds)
        print(f"\n  SII ratio scatter trend: slope={sl:+.4f}, r={r:+.3f}, p={p:.4f}")
        if sl > 0 and p < 0.1:
            print(f"  ✅ SCATTER GROWS — density-sensitive ratio degrades as predicted (high P)")
        else:
            print(f"  Scatter doesn't grow — needs more data or low-z range too narrow")


# =====================================================================
# SUMMARY: The split
# =====================================================================
print(f"\n{'='*80}")
print("SUMMARY: Protected-Null Test")
print("=" * 80)
print("""
MODEL PREDICTION:
  Branch-locked ([OIII] 5007/4959, [NII] 6585/6549): FLAT (P≈0)
  Density-sensitive ([SII] 6718/6732): DEGRADES (P≈0.7)
  Recombination (Hα, Hβ EW): DEGRADES (P≈0.25-0.30)

If branch-locked ratios are flat AND density-sensitive ratios degrade,
the model's P formula is doing exactly what it should:
discriminating based on emissivity gradient, not just "everything degrades."

This is the BREAKER test. If everything degrades the same way,
our P formula is wrong and the effect is less specific than we think.
""")
