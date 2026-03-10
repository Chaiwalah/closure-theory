#!/usr/bin/env python3
"""
Q4 ZOOM + TRANSITION EXTREMES

Q4 (galactic longitude 270-360°) showed anomalous degradation at z=0.8-1.0.
What's special about that direction? And does the sightline effect amplify
at the transition where bonds are weakest?

Also: what structures lie along Q4 sightlines?
- Galactic anti-center region
- Known galaxy clusters / voids along that direction?
- CMB Cold Spot is near l=210, b=-57 (Q3, not Q4) — rule that out
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_q4_zoom', exist_ok=True)

print("=" * 70)
print("Q4 ZOOM — What's special about l=270-360°?")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']

mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
ciii_ew = d['CIII_BR'][:, 2]
civ_ew = d['CIV'][:, 2]

ra = d['RA']
dec = d['DEC']
coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
gal_lat = coords.galactic.b.degree
gal_lon = coords.galactic.l.degree

def corr_pair(v1, v2, mask):
    d1, d2 = v1[mask], v2[mask]
    valid = np.isfinite(d1) & np.isfinite(d2) & (d1 != 0) & (d2 != 0)
    if valid.sum() > 30:
        r, p = stats.spearmanr(d1[valid], d2[valid])
        return r, int(valid.sum()), p
    return np.nan, 0, np.nan

# ============================================================
# Fine-grained longitude scan at pre-transition z
# ============================================================
print(f"\n{'='*70}")
print("LONGITUDE SCAN at z=0.8-1.0 (pre-transition)")
print("Where exactly is the anomaly?")
print(f"{'='*70}")

zmask_pretrans = (z_all >= 0.8) & (z_all < 1.0)

lon_bins = np.arange(0, 360, 15)  # 15° bins
print(f"\n  {'l range':<15} {'MgII↔Hβ':>10} {'N':>8} {'Bar':>20}")
print(f"  {'-'*15} {'-'*10} {'-'*8} {'-'*20}")

lon_corrs = []
lon_mids = []

for i in range(len(lon_bins)):
    lo = lon_bins[i]
    hi = lo + 15
    if hi > 360:
        hi = 360
    
    mask = zmask_pretrans & (gal_lon >= lo) & (gal_lon < hi)
    r, n, p = corr_pair(mg_ew, hb_ew, mask)
    
    if not np.isnan(r):
        bar = "█" * int(r * 20) if r > 0 else ""
        anomaly = " ◄◄◄" if r < 0.65 else ""
        print(f"  l=[{lo:3.0f},{hi:3.0f})    {r:>+10.3f} {n:>8} {bar}{anomaly}")
        lon_corrs.append(r)
        lon_mids.append(lo + 7.5)

lon_corrs = np.array(lon_corrs)
lon_mids = np.array(lon_mids)

if len(lon_corrs) > 3:
    mean_r = np.mean(lon_corrs)
    std_r = np.std(lon_corrs)
    print(f"\n  Mean: {mean_r:.3f} ± {std_r:.3f}")
    
    # Find outliers (>2σ from mean)
    outliers = np.abs(lon_corrs - mean_r) > 2 * std_r
    if outliers.any():
        print(f"  OUTLIERS (>2σ):")
        for idx in np.where(outliers)[0]:
            sigma = (lon_corrs[idx] - mean_r) / std_r
            print(f"    l ≈ {lon_mids[idx]:.0f}°: r = {lon_corrs[idx]:+.3f} ({sigma:+.1f}σ)")

# ============================================================
# Same scan but at transition z and post-transition
# ============================================================
for z_range, label in [((0.95, 1.05), "TRANSITION"), ((1.05, 1.2), "POST-TRANSITION")]:
    print(f"\n{'='*70}")
    print(f"LONGITUDE SCAN at z={z_range[0]}-{z_range[1]} ({label})")
    print(f"{'='*70}")
    
    zmask = (z_all >= z_range[0]) & (z_all < z_range[1])
    
    print(f"\n  {'l range':<15} {'MgII↔Hβ':>10} {'N':>8}")
    print(f"  {'-'*15} {'-'*10} {'-'*8}")
    
    for i in range(len(lon_bins)):
        lo = lon_bins[i]
        hi = lo + 15
        mask = zmask & (gal_lon >= lo) & (gal_lon < hi)
        r, n, p = corr_pair(mg_ew, hb_ew, mask)
        
        if not np.isnan(r):
            anomaly = " ◄◄◄" if abs(r) > 0.3 else ""
            print(f"  l=[{lo:3.0f},{hi:3.0f})    {r:>+10.3f} {n:>8}{anomaly}")

# ============================================================
# Latitude × Longitude grid at z=0.8-1.0
# ============================================================
print(f"\n{'='*70}")
print("SKY MAP — MgII↔Hβ at z=0.8-1.0 (lat×lon grid)")
print("Looking for hot/cold spots in correlational structure")
print(f"{'='*70}")

lat_bins = [(-90, -30), (-30, 30), (30, 90)]
lon_bins_coarse = [(0, 90), (90, 180), (180, 270), (270, 360)]

print(f"\n  {'':>15}", end="")
for llo, lhi in lon_bins_coarse:
    print(f"  l=[{llo},{lhi}){'':<{max(0,4)}}",end="")
print()

for blo, bhi in lat_bins:
    label = f"b=[{blo},{bhi})"
    print(f"  {label:>15}", end="")
    for llo, lhi in lon_bins_coarse:
        mask = zmask_pretrans & (gal_lat >= blo) & (gal_lat < bhi) & (gal_lon >= llo) & (gal_lon < lhi)
        r, n, p = corr_pair(mg_ew, hb_ew, mask)
        if not np.isnan(r):
            print(f"  {r:+.3f} ({n:>5})", end="")
        else:
            print(f"  {'N/A':>12}", end="")
    print()

# ============================================================
# MULTIPLE PAIRS at the anomalous direction
# ============================================================
print(f"\n{'='*70}")
print("MULTI-PAIR CHECK — Does Q4 anomaly affect ALL correlations?")
print("If it's a real sightline effect, ALL pairs should be affected")
print(f"{'='*70}")

q4_mask = (gal_lon >= 270) & (gal_lon < 360)
q_other = (gal_lon >= 0) & (gal_lon < 270)

pairs = [
    ('MgII↔Hβ', mg_ew, hb_ew),
    ('MgII↔CIII', mg_ew, ciii_ew),
    ('MgII↔CIV', mg_ew, civ_ew),
    ('CIII↔CIV', ciii_ew, civ_ew),
]

z_test_ranges = [(0.5, 0.8), (0.8, 1.0), (1.0, 1.2)]

print(f"\n  {'Pair':<15} {'z-range':<12} {'Q4':>8} {'Other':>8} {'Δ':>8} {'Sig':>6}")
print(f"  {'-'*15} {'-'*12} {'-'*8} {'-'*8} {'-'*8} {'-'*6}")

for pair_name, v1, v2 in pairs:
    for zlo, zhi in z_test_ranges:
        zmask = (z_all >= zlo) & (z_all < zhi)
        
        r_q4, n_q4, _ = corr_pair(v1, v2, zmask & q4_mask)
        r_other, n_other, _ = corr_pair(v1, v2, zmask & q_other)
        
        if not np.isnan(r_q4) and not np.isnan(r_other):
            delta = r_q4 - r_other
            sig = "★" if abs(delta) > 0.05 else ""
            print(f"  {pair_name:<15} {zlo}-{zhi:<8} {r_q4:>+8.3f} {r_other:>+8.3f} {delta:>+8.3f} {sig:>6}")

# ============================================================
# What's in Q4? Known large-scale structures
# ============================================================
print(f"\n{'='*70}")
print("WHAT'S IN Q4? Large-scale structure context")
print(f"{'='*70}")
print("""
  l = 270-360° includes:
  - Galactic anti-center region (~l=180° is center, so 270-360 is off-center)
  - Perseus-Pisces supercluster (l~150°, b~-15°) — NOT in Q4
  - Shapley Concentration (l~307°, b~+30°) — IN Q4!
  - Hercules supercluster (l~30°, b~+45°) — NOT in Q4
  - Dipole Repeller (l~320°, b~+10°) — IN Q4!
  
  The Shapley Concentration is the most massive structure in the 
  observable universe (~10^16 solar masses). Sightlines through it
  pass through enormous magnetic fields and hot ICM.
  
  The Dipole Repeller is a large underdensity (void) that appears
  to be pushing the Local Group away.
  
  Q4 passes through BOTH the densest and emptiest large-scale
  structures known — maximum environmental contrast.
""")

# ============================================================
# Shapley direction zoom
# ============================================================
print(f"\n{'='*70}")
print("SHAPLEY CONCENTRATION DIRECTION (l~305, b~+30)")
print(f"{'='*70}")

shapley_mask = (gal_lon >= 290) & (gal_lon < 320) & (gal_lat >= 15) & (gal_lat < 45)
anti_shapley = (gal_lon >= 110) & (gal_lon < 140) & (gal_lat >= 15) & (gal_lat < 45)

print(f"\n  MgII↔Hβ: Shapley direction vs opposite sky")
print(f"  {'z-bin':<15} {'Shapley dir':>15} {'Opposite':>15} {'Δ':>8}")
print(f"  {'-'*15} {'-'*15} {'-'*15} {'-'*8}")

for zlo, zhi in [(0.5, 0.8), (0.8, 1.0), (1.0, 1.2)]:
    zmask = (z_all >= zlo) & (z_all < zhi)
    
    r_shap, n_shap, _ = corr_pair(mg_ew, hb_ew, zmask & shapley_mask)
    r_anti, n_anti, _ = corr_pair(mg_ew, hb_ew, zmask & anti_shapley)
    
    row = f"z=[{zlo},{zhi})"
    s_str = f"{r_shap:+.3f} ({n_shap})" if not np.isnan(r_shap) else "N/A"
    a_str = f"{r_anti:+.3f} ({n_anti})" if not np.isnan(r_anti) else "N/A"
    
    if not np.isnan(r_shap) and not np.isnan(r_anti):
        delta = r_shap - r_anti
        print(f"  {row:<15} {s_str:>15} {a_str:>15} {delta:>+8.3f}")
    else:
        print(f"  {row:<15} {s_str:>15} {a_str:>15} {'N/A':>8}")

print(f"\nResults saved to results_q4_zoom/")

with open('results_q4_zoom/q4_zoom_results.json', 'w') as f:
    json.dump({'note': 'see console output'}, f)
