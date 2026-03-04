#!/usr/bin/env python3
"""
TEST 4: REPEAT-SPECTRA INVARIANCE TEST
========================================
Same source, same path, different instrument state.
If the effect is real (tied to line-of-sight), the "damage score" should be
REPEATABLE between epochs. If it's pipeline noise, it should vary with
plate/fiber/MJD.

DR16Q has PLATE, MJD, FIBERID. Objects observed on multiple plates
are repeat observations. Find them and test repeatability.

Actually: DR16Q is one row per quasar (best observation). We need the
full spectroscopic catalog to find repeats. BUT we can use PLATE+MJD+FIBER
to identify objects observed on different plates (some quasars are in
plate overlaps). 

Alternative approach: Use the DR16Q SDSS_NAME to find duplicates
(same object, different observation), or check if there's a NSPEC column.
"""

import numpy as np
from scipy import stats
from astropy.io import fits
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = Path('results_repeat_spectra')
RESULTS_DIR.mkdir(exist_ok=True)

f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z = d['Z_DR16Q']
ra = d['RA']
dec = d['DEC']
plate = d['PLATE']
mjd = d['MJD']
snr = d['SN_MEDIAN_ALL']

all_cols = [c.name for c in f[1].columns]
print("Checking for repeat/duplicate indicators...")
for c in sorted(all_cols):
    if any(k in c.upper() for k in ['NSPEC', 'NOBS', 'DUP', 'REPEAT', 'PRIM', 'SDSS_NAME', 'THING_ID']):
        print(f"  {c}")

# Check for SDSS_NAME or similar identifier
sdss_name = None
for cname in ['SDSS_NAME', 'SDSSJ', 'NAME']:
    if cname in all_cols:
        sdss_name = d[cname]
        print(f"  Found name column: {cname}")
        break

lines = {}
for name in ['MGII', 'CIV', 'CIII_ALL', 'HBETA', 'OIII5007', 'LYA']:
    col = d[name]
    if col.ndim == 2 and col.shape[1] >= 6:
        lines[name] = {'ew': col[:, 2], 'fwhm': col[:, 4]}

f.close()

print(f"\n{'='*80}")
print("🔁 TEST 4: REPEAT-SPECTRA INVARIANCE")
print("=" * 80)

# =====================================================================
# APPROACH 1: Find close positional duplicates (plate overlap regions)
# =====================================================================
print("\nFinding positional duplicates (same object, different plate)...")
print("Objects within 1 arcsec of each other on different plates = repeats")

# Bin by rounded position
bin_size = 0.005  # ~18 arcsec bins
pos_bins = defaultdict(list)
for i in range(len(ra)):
    if np.isfinite(ra[i]) and np.isfinite(dec[i]):
        key = (round(ra[i] / bin_size), round(dec[i] / bin_size))
        pos_bins[key].append(i)

# Find groups with multiple entries on different plates
match_tol = 1.0 / 3600  # 1 arcsec
repeat_groups = []

for key, indices in pos_bins.items():
    if len(indices) < 2:
        continue
    # Check all pairs
    for a in range(len(indices)):
        for b in range(a+1, len(indices)):
            i, j = indices[a], indices[b]
            dist = np.sqrt((ra[i] - ra[j])**2 + ((dec[i] - dec[j]) * np.cos(np.radians(dec[i])))**2)
            if dist < match_tol and plate[i] != plate[j]:
                repeat_groups.append((i, j))

print(f"  Found {len(repeat_groups)} repeat pairs (same position, different plate)")

if len(repeat_groups) < 50:
    print("  Too few plate-overlap repeats. Trying wider tolerance (2 arcsec)...")
    match_tol = 2.0 / 3600
    repeat_groups = []
    for key, indices in pos_bins.items():
        if len(indices) < 2:
            continue
        for nk_ra in (-1, 0, 1):
            for nk_dec in (-1, 0, 1):
                nkey = (key[0] + nk_ra, key[1] + nk_dec)
                for j_idx in pos_bins.get(nkey, []):
                    for i_idx in indices:
                        if i_idx >= j_idx:
                            continue
                        if plate[i_idx] == plate[j_idx]:
                            continue
                        dist = np.sqrt((ra[i_idx] - ra[j_idx])**2 + 
                                      ((dec[i_idx] - dec[j_idx]) * np.cos(np.radians(dec[i_idx])))**2)
                        if dist < match_tol:
                            repeat_groups.append((i_idx, j_idx))
    
    # Deduplicate
    repeat_groups = list(set(repeat_groups))
    print(f"  Found {len(repeat_groups)} repeat pairs at 2 arcsec")

# =====================================================================
# APPROACH 2: Use the catalog as-is but test plate-level consistency
# =====================================================================
# Even without true repeats, we can test: for quasars on the same plate
# at the same z, is the damage score more consistent within-plate than
# between-plate? If pipeline dominates, within-plate should be tighter.

print(f"\n{'='*80}")
print("APPROACH 2: PLATE-LEVEL CONSISTENCY TEST")
print("If damage is LOS-locked: between-plate scatter at same z ≈ within-plate scatter")
print("If damage is pipeline: within-plate < between-plate (plates have consistent biases)")
print("=" * 80)

# Focus on MgII (pipeline-confirmed)
if 'MGII' in lines:
    ew_mg = lines['MGII']['ew']
    fw_mg = lines['MGII']['fwhm']
    
    # Damage score = log(|EW|/FWHM) — how much EW "sticks out" relative to kinematics
    valid = np.isfinite(ew_mg) & (ew_mg != 0) & np.isfinite(fw_mg) & (fw_mg > 0)
    valid &= np.isfinite(z) & (z > 0.5) & (z < 2.0) & np.isfinite(snr)
    
    damage = np.full(len(z), np.nan)
    damage[valid] = np.log10(np.abs(ew_mg[valid]) / fw_mg[valid])
    
    # For narrow z-bins, compare within-plate vs between-plate variance
    print(f"\n  MgII damage score: within-plate vs between-plate variance")
    print(f"  {'z-bin':<10} {'σ_within':<12} {'σ_between':<12} {'ratio':<10} {'Plates':<8} {'N'}")
    print(f"  {'-'*60}")
    
    for z_lo, z_hi in [(0.5, 0.7), (0.7, 0.9), (0.9, 1.1), (1.1, 1.4), (1.4, 2.0)]:
        z_mask = valid & (z > z_lo) & (z < z_hi)
        
        if z_mask.sum() < 200:
            continue
        
        # Group by plate
        plate_groups = defaultdict(list)
        for idx in np.where(z_mask)[0]:
            plate_groups[plate[idx]].append(damage[idx])
        
        # Only use plates with 10+ quasars
        plate_means = []
        within_vars = []
        for p, vals in plate_groups.items():
            if len(vals) >= 10:
                plate_means.append(np.mean(vals))
                within_vars.append(np.var(vals))
        
        if len(plate_means) < 5:
            continue
        
        sigma_between = np.std(plate_means)
        sigma_within = np.sqrt(np.mean(within_vars))
        ratio = sigma_between / sigma_within if sigma_within > 0 else np.nan
        
        tag = ""
        if ratio < 0.1:
            tag = "→ Pipeline-dominated (plates have own bias)"
        elif ratio > 0.5:
            tag = "→ LOS-dominated (damage varies within plate)"
        else:
            tag = "→ Mixed"
        
        print(f"  {z_lo:.1f}-{z_hi:.1f}    {sigma_within:.4f}      {sigma_between:.4f}      {ratio:.3f}     {len(plate_means):<8} {z_mask.sum()}")
        print(f"  {'':<10} {tag}")


# =====================================================================
# If we found repeats, test them directly
# =====================================================================
if len(repeat_groups) >= 20:
    print(f"\n{'='*80}")
    print(f"DIRECT REPEAT TEST: {len(repeat_groups)} pairs")
    print("=" * 80)
    
    for lname in ['MGII', 'CIV']:
        if lname not in lines:
            continue
        
        ew = lines[lname]['ew']
        fw = lines[lname]['fwhm']
        
        ew_diffs = []
        fw_diffs = []
        snr_diffs = []
        
        for i, j in repeat_groups:
            if np.isfinite(ew[i]) and np.isfinite(ew[j]) and ew[i] != 0 and ew[j] != 0:
                if np.isfinite(fw[i]) and np.isfinite(fw[j]) and fw[i] > 0 and fw[j] > 0:
                    d_i = np.log10(np.abs(ew[i]) / fw[i])
                    d_j = np.log10(np.abs(ew[j]) / fw[j])
                    ew_diffs.append(d_i - d_j)
                    
                    fw_diffs.append(np.log10(fw[i]) - np.log10(fw[j]))
                    
                    if np.isfinite(snr[i]) and np.isfinite(snr[j]):
                        snr_diffs.append(snr[i] - snr[j])
        
        if len(ew_diffs) >= 10:
            ew_diffs = np.array(ew_diffs)
            fw_diffs = np.array(fw_diffs)
            
            print(f"\n  {lname} repeat pairs (N={len(ew_diffs)}):")
            print(f"    Damage score Δ: mean={np.mean(ew_diffs):+.4f}, std={np.std(ew_diffs):.4f}")
            print(f"    FWHM Δ:         mean={np.mean(fw_diffs):+.4f}, std={np.std(fw_diffs):.4f}")
            
            # ICC-like: is damage more repeatable than expected from noise?
            icc_damage = 1 - np.var(ew_diffs) / (2 * np.var([d for pair in repeat_groups for d in [np.log10(np.abs(ew[pair[0]]) / fw[pair[0]]) if np.isfinite(ew[pair[0]]) and ew[pair[0]] != 0 and fw[pair[0]] > 0 else np.nan]]))
            print(f"    Damage ICC (approx): {icc_damage:.3f}")
            
            if np.std(ew_diffs) < 0.1:
                print(f"    ✅ Damage is HIGHLY REPEATABLE — LOS-locked, not pipeline noise")
            else:
                print(f"    ⚠️ Damage has significant epoch-to-epoch variation")


print(f"\n{'='*80}")
print("SUMMARY")
print("=" * 80)
print("""
A real propagation effect should be:
  - Object-locked (same sightline = same damage, regardless of observation epoch)
  - Within-plate variance >> between-plate variance (damage varies with sky, not instrument)

A pipeline artifact should be:
  - Epoch-locked (varies with SNR, plate, fiber)
  - Within-plate variance << between-plate variance (each plate has its own bias)
""")
