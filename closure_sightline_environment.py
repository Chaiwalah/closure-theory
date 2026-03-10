#!/usr/bin/env python3
"""
SIGHTLINE ENVIRONMENT TESTS — Does the PATH matter, not just the distance?

If the wave equilibrates with its environment, then:
- Sightlines through magnetic-rich regions → MORE degradation
- Sightlines through voids → LESS degradation
- At the SAME z.

Tests:
1. GALACTIC LATITUDE: High |b| (through galactic pole, less MW magnetic field)
   vs Low |b| (through galactic plane, more MW magnetic field)
   → Different degradation at same z?

2. RADIO-LOUD vs RADIO-QUIET: Radio-loud quasars live in denser environments
   (more magnetic fields, more IGM) → Different correlation structure?

3. ECLIPTIC LATITUDE: Controls for solar system magnetic environment
   → If this matters, something very local is affecting "cosmological" correlations

4. SKY POSITION DEPENDENCE: Split sky into quadrants
   → Correlational structure should be isotropic if z is the only variable
   → Anisotropy = path matters

5. DUST EXTINCTION (E(B-V)): More dust = more ISM along sightline
   → If dust column predicts correlational degradation beyond z, the medium matters
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_sightline', exist_ok=True)

print("=" * 70)
print("SIGHTLINE ENVIRONMENT — Does the PATH change correlations?")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']

# Key observables
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
ciii_ew = d['CIII_BR'][:, 2]
civ_ew = d['CIV'][:, 2]

# Sky coordinates
ra = d['RA']
dec = d['DEC']

# Convert to galactic coordinates
print("\nConverting coordinates...")
coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
gal_lat = coords.galactic.b.degree  # galactic latitude
gal_lon = coords.galactic.l.degree
ecl_lat = coords.barycentrictrueecliptic.lat.degree  # ecliptic latitude

print(f"Galactic latitude range: {gal_lat.min():.1f} to {gal_lat.max():.1f}")
print(f"Ecliptic latitude range: {ecl_lat.min():.1f} to {ecl_lat.max():.1f}")

# ============================================================
# Helper: compute MgII↔Hβ correlation in a subsample
# ============================================================
def corr_mghb(mask):
    v1, v2 = mg_ew[mask], hb_ew[mask]
    valid = np.isfinite(v1) & np.isfinite(v2) & (v1 != 0) & (v2 != 0)
    if valid.sum() > 30:
        r, p = stats.spearmanr(v1[valid], v2[valid])
        return r, valid.sum()
    return np.nan, 0

def corr_mgciii(mask):
    v1, v2 = mg_ew[mask], ciii_ew[mask]
    valid = np.isfinite(v1) & np.isfinite(v2) & (v1 != 0) & (v2 != 0)
    if valid.sum() > 30:
        r, p = stats.spearmanr(v1[valid], v2[valid])
        return r, valid.sum()
    return np.nan, 0

def corr_mgciv(mask):
    v1, v2 = mg_ew[mask], civ_ew[mask]
    valid = np.isfinite(v1) & np.isfinite(v2) & (v1 != 0) & (v2 != 0)
    if valid.sum() > 30:
        r, p = stats.spearmanr(v1[valid], v2[valid])
        return r, valid.sum()
    return np.nan, 0

# ============================================================
# TEST 1: GALACTIC LATITUDE
# ============================================================
print(f"\n{'='*70}")
print("TEST 1: GALACTIC LATITUDE — Magnetic field of Milky Way")
print("High |b| = through galactic pole (less MW B-field)")
print("Low |b| = through galactic plane (more MW B-field)")
print(f"{'='*70}")

abs_b = np.abs(gal_lat)
b_cuts = [(20, 40, "Low |b| (20-40°)"), (40, 60, "Mid |b| (40-60°)"), (60, 90, "High |b| (60-90°)")]

z_bins = [(0.5, 0.8), (0.8, 1.0), (1.0, 1.2), (1.2, 1.5), (1.5, 2.0)]

print(f"\n  MgII↔Hβ:")
print(f"  {'z-bin':<15}", end="")
for _, _, label in b_cuts:
    print(f"  {label:>20}", end="")
print(f"  {'Δ(hi-lo)':>10}")
print(f"  {'-'*15}", end="")
for _ in b_cuts:
    print(f"  {'-'*20}", end="")
print(f"  {'-'*10}")

sightline_results = {}

for zlo, zhi in z_bins:
    zmask = (z_all >= zlo) & (z_all < zhi)
    row = f"z=[{zlo},{zhi})"
    print(f"  {row:<15}", end="")
    
    corrs_by_b = []
    for blo, bhi, label in b_cuts:
        mask = zmask & (abs_b >= blo) & (abs_b < bhi)
        r, n = corr_mghb(mask)
        corrs_by_b.append(r)
        r_str = f"{r:+.3f} (N={n})" if not np.isnan(r) else "N/A"
        print(f"  {r_str:>20}", end="")
    
    if len(corrs_by_b) >= 2 and not np.isnan(corrs_by_b[0]) and not np.isnan(corrs_by_b[-1]):
        delta = corrs_by_b[-1] - corrs_by_b[0]
        print(f"  {delta:>+10.3f}")
        sightline_results[row] = {
            'low_b': float(corrs_by_b[0]),
            'high_b': float(corrs_by_b[-1]),
            'delta': float(delta),
        }
    else:
        print(f"  {'N/A':>10}")

# ============================================================
# TEST 2: ECLIPTIC LATITUDE (solar system environment)
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: ECLIPTIC LATITUDE — Solar system magnetic environment")
print("If THIS matters, something very local affects 'cosmological' correlations")
print(f"{'='*70}")

abs_ecl = np.abs(ecl_lat)
ecl_cuts = [(0, 30, "Low |β| (0-30°)"), (30, 60, "Mid |β| (30-60°)"), (60, 90, "High |β| (60-90°)")]

print(f"\n  MgII↔Hβ:")
print(f"  {'z-bin':<15}", end="")
for _, _, label in ecl_cuts:
    print(f"  {label:>20}", end="")
print(f"  {'Δ(hi-lo)':>10}")
print(f"  {'-'*15}", end="")
for _ in ecl_cuts:
    print(f"  {'-'*20}", end="")
print(f"  {'-'*10}")

for zlo, zhi in z_bins:
    zmask = (z_all >= zlo) & (z_all < zhi)
    row = f"z=[{zlo},{zhi})"
    print(f"  {row:<15}", end="")
    
    corrs = []
    for elo, ehi, label in ecl_cuts:
        mask = zmask & (abs_ecl >= elo) & (abs_ecl < ehi)
        r, n = corr_mghb(mask)
        corrs.append(r)
        r_str = f"{r:+.3f} (N={n})" if not np.isnan(r) else "N/A"
        print(f"  {r_str:>20}", end="")
    
    if len(corrs) >= 2 and not np.isnan(corrs[0]) and not np.isnan(corrs[-1]):
        delta = corrs[-1] - corrs[0]
        print(f"  {delta:>+10.3f}")
    else:
        print(f"  {'N/A':>10}")

# ============================================================
# TEST 3: SKY QUADRANT ANISOTROPY
# ============================================================
print(f"\n{'='*70}")
print("TEST 3: SKY QUADRANT — Is correlational degradation isotropic?")
print("If z is the only variable, all quadrants should agree")
print(f"{'='*70}")

# Split by galactic longitude quadrants
quad_cuts = [
    (0, 90, "Q1 (l=0-90)"),
    (90, 180, "Q2 (l=90-180)"),
    (180, 270, "Q3 (l=180-270)"),
    (270, 360, "Q4 (l=270-360)"),
]

print(f"\n  MgII↔Hβ:")
print(f"  {'z-bin':<15}", end="")
for _, _, label in quad_cuts:
    print(f"  {label:>15}", end="")
print(f"  {'Range':>8}")
print(f"  {'-'*15}", end="")
for _ in quad_cuts:
    print(f"  {'-'*15}", end="")
print(f"  {'-'*8}")

for zlo, zhi in z_bins:
    zmask = (z_all >= zlo) & (z_all < zhi)
    row = f"z=[{zlo},{zhi})"
    print(f"  {row:<15}", end="")
    
    corrs = []
    for qlo, qhi, label in quad_cuts:
        mask = zmask & (gal_lon >= qlo) & (gal_lon < qhi)
        r, n = corr_mghb(mask)
        corrs.append(r)
        r_str = f"{r:+.3f}" if not np.isnan(r) else "N/A"
        print(f"  {r_str:>15}", end="")
    
    valid_corrs = [c for c in corrs if not np.isnan(c)]
    if len(valid_corrs) >= 2:
        rng = max(valid_corrs) - min(valid_corrs)
        print(f"  {rng:>8.3f}")
    else:
        print(f"  {'N/A':>8}")

# ============================================================
# TEST 4: DUST EXTINCTION (GAL_EXT from DR16Q)
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: DUST EXTINCTION — Does ISM column predict degradation?")
print(f"{'='*70}")

# Check if extinction column exists
ext_cols = [c for c in d.dtype.names if 'EXT' in c.upper() or 'EBV' in c.upper() or 'DUST' in c.upper()]
print(f"  Available extinction columns: {ext_cols}")

if 'EXTINCTION' in d.dtype.names:
    ext = d['EXTINCTION']
    # SDSS extinction is multi-band array
    if ext.ndim == 2:
        ext_r = ext[:, 2]  # r-band extinction
    else:
        ext_r = ext
    
    ext_valid = np.isfinite(ext_r) & (ext_r > 0)
    ext_lo = np.percentile(ext_r[ext_valid], 25)
    ext_hi = np.percentile(ext_r[ext_valid], 75)
    
    print(f"  r-band extinction: median={np.median(ext_r[ext_valid]):.3f}, Q1={ext_lo:.3f}, Q3={ext_hi:.3f}")
    
    print(f"\n  MgII↔Hβ at fixed z, split by dust column:")
    print(f"  {'z-bin':<15} {'Low dust':>15} {'High dust':>15} {'Δ':>8}")
    print(f"  {'-'*15} {'-'*15} {'-'*15} {'-'*8}")
    
    for zlo, zhi in z_bins:
        zmask = (z_all >= zlo) & (z_all < zhi) & ext_valid
        
        mask_lo = zmask & (ext_r <= ext_lo)
        mask_hi = zmask & (ext_r >= ext_hi)
        
        r_lo, n_lo = corr_mghb(mask_lo)
        r_hi, n_hi = corr_mghb(mask_hi)
        
        row = f"z=[{zlo},{zhi})"
        lo_str = f"{r_lo:+.3f} ({n_lo})" if not np.isnan(r_lo) else "N/A"
        hi_str = f"{r_hi:+.3f} ({n_hi})" if not np.isnan(r_hi) else "N/A"
        
        if not np.isnan(r_lo) and not np.isnan(r_hi):
            delta = r_hi - r_lo
            print(f"  {row:<15} {lo_str:>15} {hi_str:>15} {delta:>+8.3f}")
        else:
            print(f"  {row:<15} {lo_str:>15} {hi_str:>15} {'N/A':>8}")

# ============================================================
# TEST 5: TRANSITION ZONE — Does environment matter MORE at z≈1?
# ============================================================
print(f"\n{'='*70}")
print("TEST 5: ENVIRONMENT × TRANSITION — Does the path matter more near z₀?")
print("If wave equilibration is the mechanism, sightline effects should")
print("be STRONGEST near the transition (z≈1) where bonds are weakest")
print(f"{'='*70}")

# Use galactic latitude as environment proxy
# Compare |b| effect at different z
z_fine = [(0.5, 0.7), (0.7, 0.85), (0.85, 0.95), (0.95, 1.05), 
          (1.05, 1.15), (1.15, 1.3), (1.3, 1.6), (1.6, 2.0)]

print(f"\n  MgII↔Hβ: High |b| (>50°) vs Low |b| (<35°)")
print(f"  {'z-bin':<15} {'Low |b|':>12} {'High |b|':>12} {'Δ':>8} {'|Δ| bar':>20}")
print(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*8} {'-'*20}")

deltas = []
z_mids = []

for zlo, zhi in z_fine:
    zmask = (z_all >= zlo) & (z_all < zhi)
    
    mask_lo_b = zmask & (abs_b < 35)
    mask_hi_b = zmask & (abs_b > 50)
    
    r_lo, n_lo = corr_mghb(mask_lo_b)
    r_hi, n_hi = corr_mghb(mask_hi_b)
    
    row = f"z=[{zlo},{zhi})"
    
    if not np.isnan(r_lo) and not np.isnan(r_hi):
        delta = r_hi - r_lo
        bar = "█" * int(abs(delta) * 100)
        sign = "+" if delta > 0 else "-"
        print(f"  {row:<15} {r_lo:>+12.3f} {r_hi:>+12.3f} {delta:>+8.3f} {bar}")
        deltas.append(delta)
        z_mids.append((zlo + zhi) / 2)
    else:
        print(f"  {row:<15} {'N/A':>12} {'N/A':>12}")

if len(deltas) > 3:
    # Does the environment effect peak at the transition?
    deltas = np.array(deltas)
    z_mids = np.array(z_mids)
    peak_idx = np.argmax(np.abs(deltas))
    print(f"\n  ★ LARGEST environment effect at z ≈ {z_mids[peak_idx]:.2f} (Δ = {deltas[peak_idx]:+.3f})")
    
    # Does it correlate with transition proximity?
    dist_from_transition = np.abs(z_mids - 1.045)
    r_env, p_env = stats.spearmanr(dist_from_transition, np.abs(deltas))
    print(f"  Environment effect vs distance from z₀: ρ = {r_env:+.3f} (p = {p_env:.3e})")
    if r_env < -0.3 and p_env < 0.05:
        print(f"  ⚠️ Environment matters MORE near the transition — wave equilibration supported!")

# Save
with open('results_sightline/sightline_results.json', 'w') as f:
    json.dump(sightline_results, f, indent=2, default=str)

print(f"\n\nResults saved to results_sightline/")
