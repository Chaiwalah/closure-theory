#!/usr/bin/env python3
"""
TEST 5 REDO: Protected-Null using Rakshit+2020
================================================
Rakshit has SEPARATE measurements for doublet members.
Test: branch-locked ratios should be FLAT with z.

From Rakshit ReadMe:
  EWNII6549: bytes 903-913
  EWNII6585: bytes 953-963
  EWSII6718: bytes 1003-1013
  EWSII6732: bytes 1053-1064  ← THIS IS THE KEY
  EWOIII4959C: bytes 1341-1352
  EWOIII5007C: bytes 1367-1378
  
Branch-locked pairs (same upper level):
  [OIII] 5007/4959 — ratio = 2.98 (A-coefficients)
  [NII] 6585/6549 — ratio = 2.95

Density-sensitive pair (different upper levels):  
  [SII] 6718/6732 — ratio varies 0.4-1.4 with n_e
"""

import numpy as np
from scipy import stats
import gzip
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = Path('results_protected_null_v2')
RESULTS_DIR.mkdir(exist_ok=True)

print("=" * 80)
print("🛡️ TEST 5v2: PROTECTED-NULL with Rakshit+2020 Doublets")
print("=" * 80)

# Column positions (0-indexed bytes)
col_defs = {
    'z':              (94, 102),
    'snr':            (103, 111),
    'ew_nii6549':     (902, 913),
    'ew_nii6585':     (952, 963),
    'ew_sii6718':     (1002, 1013),
    'ew_sii6732':     (1052, 1064),
    'ew_oiii4959c':   (1340, 1352),
    'ew_oiii5007c':   (1366, 1378),
    'ew_hbeta_br':    (1392, 1404),  # broad Hβ EW (anti-null)
    'ew_halpha_br':   (1130, 1142),  # broad Hα EW
    'fwhm_hbeta_br':  (1452, 1464),  # broad Hβ FWHM
}

print("Loading Rakshit+2020 (extracting doublet columns)...")
rak = {k: [] for k in col_defs}

with gzip.open('/root/clawd/data/sdss/rakshit2020.dat.gz', 'rt', encoding='ascii', errors='replace') as fh:
    for i, line in enumerate(fh):
        if i % 100000 == 0:
            print(f"  Read {i} lines...")
        for col_name, (start, end) in col_defs.items():
            try:
                val = float(line[start:end].strip())
            except (ValueError, IndexError):
                val = np.nan
            rak[col_name].append(val)

for k in rak:
    rak[k] = np.array(rak[k])
    if k != 'z' and k != 'snr':
        rak[k][rak[k] <= -998] = np.nan

z = rak['z']
n = len(z)
print(f"  Loaded {n} quasars")

# =====================================================================
# TEST 5A: [OIII] 5007/4959 ratio vs z (BRANCH-LOCKED)
# =====================================================================
print(f"\n{'='*80}")
print("[OIII] 5007C / 4959C — Branch-locked (theoretical ratio ≈ 2.98)")
print("Should be FLAT with z if our model is correct (P ≈ 0)")
print("=" * 80)

ew_5007 = rak['ew_oiii5007c']
ew_4959 = rak['ew_oiii4959c']

valid = np.isfinite(ew_5007) & np.isfinite(ew_4959) & (ew_4959 != 0) & (ew_5007 != 0)
valid &= np.isfinite(z) & (z > 0.01)

ratio = ew_5007[valid] / ew_4959[valid]
z_v = z[valid]

# Clip to reasonable range
good = (ratio > 0.5) & (ratio < 10)
ratio_g = ratio[good]
z_g = z_v[good]

print(f"  Valid measurements: {len(ratio_g)}")
print(f"  Overall: median = {np.median(ratio_g):.3f}, std = {np.std(ratio_g):.3f}")

z_bins = [(0.02, 0.2), (0.2, 0.35), (0.35, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.2)]
print(f"\n  {'z-bin':<12} {'median':<10} {'std':<10} {'N':<8}")
print(f"  {'-'*40}")

meds = []
stds = []
zmids = []
for z_lo, z_hi in z_bins:
    m = (z_g > z_lo) & (z_g < z_hi)
    if m.sum() < 20:
        continue
    med = np.median(ratio_g[m])
    std = np.std(ratio_g[m])
    meds.append(med)
    stds.append(std)
    zmids.append((z_lo + z_hi) / 2)
    print(f"  {z_lo:.2f}-{z_hi:.2f}    {med:.4f}     {std:.3f}     {m.sum()}")

if len(zmids) >= 3:
    sl_med, _, r_med, p_med, _ = stats.linregress(zmids, meds)
    sl_std, _, r_std, p_std, _ = stats.linregress(zmids, stds)
    print(f"\n  Median trend: slope={sl_med:+.4f}, r={r_med:+.3f}, p={p_med:.4f}")
    print(f"  Scatter trend: slope={sl_std:+.4f}, r={r_std:+.3f}, p={p_std:.4f}")
    if p_med > 0.05 and p_std > 0.05:
        print(f"  ✅ FLAT — [OIII] branch ratio is IMMUNE to z. P≈0 confirmed.")
    elif p_std < 0.05 and sl_std > 0:
        print(f"  ⚠️ Scatter grows but median stable — measurement noise?")
    else:
        print(f"  ⚠️ Trend detected — needs investigation")


# =====================================================================
# TEST 5B: [NII] 6585/6549 ratio vs z (BRANCH-LOCKED)
# =====================================================================
print(f"\n{'='*80}")
print("[NII] 6585 / 6549 — Branch-locked (theoretical ratio ≈ 2.95)")
print("=" * 80)

ew_6585 = rak['ew_nii6585']
ew_6549 = rak['ew_nii6549']

valid = np.isfinite(ew_6585) & np.isfinite(ew_6549) & (ew_6549 != 0) & (ew_6585 != 0)
valid &= np.isfinite(z) & (z > 0.01)

ratio = ew_6585[valid] / ew_6549[valid]
z_v = z[valid]

good = (ratio > 0.5) & (ratio < 10)
ratio_g = ratio[good]
z_g = z_v[good]

print(f"  Valid: {len(ratio_g)}, median = {np.median(ratio_g):.3f}")

z_bins_nii = [(0.02, 0.15), (0.15, 0.3), (0.3, 0.45), (0.45, 0.65)]
meds_n = []
zmids_n = []
for z_lo, z_hi in z_bins_nii:
    m = (z_g > z_lo) & (z_g < z_hi)
    if m.sum() < 20:
        continue
    med = np.median(ratio_g[m])
    meds_n.append(med)
    zmids_n.append((z_lo + z_hi) / 2)
    print(f"  z={z_lo:.2f}-{z_hi:.2f}: median={med:.4f}, std={np.std(ratio_g[m]):.3f}, N={m.sum()}")

if len(zmids_n) >= 3:
    sl, _, r, p, _ = stats.linregress(zmids_n, meds_n)
    print(f"  Trend: slope={sl:+.4f}, r={r:+.3f}, p={p:.4f}")


# =====================================================================
# TEST 5C: [SII] 6718/6732 ratio vs z (DENSITY-SENSITIVE = should DEGRADE)
# =====================================================================
print(f"\n{'='*80}")
print("[SII] 6718 / 6732 — Density-sensitive (P ≈ 0.7, should DEGRADE)")
print("Different upper levels, n_crit = 1400 vs 3600 cm⁻³")
print("=" * 80)

ew_6718 = rak['ew_sii6718']
ew_6732 = rak['ew_sii6732']

valid = np.isfinite(ew_6718) & np.isfinite(ew_6732) & (ew_6718 != 0) & (ew_6732 != 0)
valid &= np.isfinite(z) & (z > 0.01)

ratio = ew_6718[valid] / ew_6732[valid]
z_v = z[valid]

good = (ratio > 0.1) & (ratio < 10)
ratio_g = ratio[good]
z_g = z_v[good]

print(f"  Valid: {len(ratio_g)}, median = {np.median(ratio_g):.3f}")

z_bins_sii = [(0.02, 0.1), (0.1, 0.2), (0.2, 0.35), (0.35, 0.55)]
stds_s = []
zmids_s = []
for z_lo, z_hi in z_bins_sii:
    m = (z_g > z_lo) & (z_g < z_hi)
    if m.sum() < 20:
        continue
    med = np.median(ratio_g[m])
    std = np.std(ratio_g[m])
    stds_s.append(std)
    zmids_s.append((z_lo + z_hi) / 2)
    print(f"  z={z_lo:.2f}-{z_hi:.2f}: median={med:.4f}, std={std:.3f}, N={m.sum()}")

if len(zmids_s) >= 3:
    sl, _, r, p, _ = stats.linregress(zmids_s, stds_s)
    print(f"  Scatter trend: slope={sl:+.4f}, r={r:+.3f}, p={p:.4f}")
    if sl > 0:
        print(f"  ✅ Scatter GROWS — density-sensitive ratio degrades. High P confirmed.")
    else:
        print(f"  Scatter doesn't grow — may need higher z (these are all sub-threshold)")


# =====================================================================
# COMPARISON TABLE
# =====================================================================
print(f"\n{'='*80}")
print("COMPARISON: Protected vs Vulnerable Ratios")
print("=" * 80)
print("""
  Observable              Type                P    Prediction     Result
  ─────────────────────────────────────────────────────────────────────
  [OIII] 5007/4959        Branch-locked       0    FLAT           ?
  [NII] 6585/6549         Branch-locked       0    FLAT           ?
  [SII] 6718/6732         Density-sensitive   0.7  DEGRADES       ?
  Hβ broad EW             Recombination       0.3  DEGRADES       (done earlier)
  
If the split is clean (locked=flat, sensitive=degrades), the P formula
is correctly discriminating based on emissivity gradients.
""")
