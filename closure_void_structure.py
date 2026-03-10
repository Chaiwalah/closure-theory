#!/usr/bin/env python3
"""
VOID vs STRUCTURE — The test everyone agreed on.

Does foreground structure PROTECT and void BREAK?

Approach: Use DR16Q quasars at z=0.1-0.5 as foreground structure tracers.
For each transition-zone quasar (z=0.8-1.15), count how many foreground
quasars lie within some angular radius. More foreground = denser sightline.

Also: use Gemini's microlensing idea — check for intervening absorbers
in DR16Q metadata if available.

Tests:
1. Foreground quasar density vs early breaking
2. Absorber count vs early breaking (if data available)
3. Foreground density vs late holding
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
from scipy.spatial import cKDTree
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_void_structure', exist_ok=True)

print("=" * 70)
print("VOID vs STRUCTURE — Does dense foreground protect?")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
logmbh = d['LOGMBH']
loglbol = d['LOGLBOL']
ra_q = d['RA']
dec_q = d['DEC']

coords_all = SkyCoord(ra=ra_q*u.degree, dec=dec_q*u.degree, frame='icrs')
gal_lat = coords_all.galactic.b.degree
gal_lon = coords_all.galactic.l.degree

# ============================================================
# Build foreground density map using low-z quasars
# ============================================================
print("\nBuilding foreground density from z=0.1-0.5 quasars...")

fg_mask = (z_all >= 0.1) & (z_all < 0.5)
fg_ra = np.radians(ra_q[fg_mask])
fg_dec = np.radians(dec_q[fg_mask])

# Convert to 3D cartesian for KD-tree (on unit sphere)
fg_x = np.cos(fg_dec) * np.cos(fg_ra)
fg_y = np.cos(fg_dec) * np.sin(fg_ra)
fg_z = np.sin(fg_dec)
fg_xyz = np.column_stack([fg_x, fg_y, fg_z])

tree_fg = cKDTree(fg_xyz)
print(f"Foreground tracers: {fg_mask.sum()}")

# ============================================================
# Recreate groups (same as before)
# ============================================================
print("Recomputing rank agreements...")

trans_mask = (z_all >= 0.75) & (z_all < 1.15)
valid = trans_mask & np.isfinite(mg_ew) & np.isfinite(hb_ew) & (mg_ew != 0) & (hb_ew != 0)
z_v = z_all[valid]
mg_v = mg_ew[valid]
hb_v = hb_ew[valid]

window = 0.05
rank_agreement = np.full(valid.sum(), np.nan)

for i in range(valid.sum()):
    z_i = z_v[i]
    local = (z_v >= z_i - window) & (z_v < z_i + window)
    if local.sum() < 30:
        continue
    mg_pct = stats.percentileofscore(mg_v[local], mg_v[i]) / 100
    hb_pct = stats.percentileofscore(hb_v[local], hb_v[i]) / 100
    rank_agreement[i] = 1 - abs(mg_pct - hb_pct)

valid_ra = np.isfinite(rank_agreement)
valid_indices = np.where(valid)[0]

early_zone = (z_v >= 0.80) & (z_v < 0.95)
late_zone = (z_v >= 1.00) & (z_v < 1.15)

ra_early = rank_agreement[early_zone & valid_ra]
ra_late = rank_agreement[late_zone & valid_ra]

early_thresh = np.percentile(ra_early, 20)
late_thresh = np.percentile(ra_late, 80)

early_breakers = early_zone & valid_ra & (rank_agreement <= early_thresh)
early_normal = early_zone & valid_ra & (rank_agreement > np.percentile(ra_early, 40)) & (rank_agreement < np.percentile(ra_early, 60))
late_holders = late_zone & valid_ra & (rank_agreement >= late_thresh)
late_normal = late_zone & valid_ra & (rank_agreement > np.percentile(ra_late, 40)) & (rank_agreement < np.percentile(ra_late, 60))

eb_idx = valid_indices[early_breakers]
en_idx = valid_indices[early_normal]
lh_idx = valid_indices[late_holders]
ln_idx = valid_indices[late_normal]

print(f"Groups: EB={len(eb_idx)}, EN={len(en_idx)}, LH={len(lh_idx)}, LN={len(ln_idx)}")

# ============================================================
# Compute foreground density for each group
# ============================================================
print(f"\n{'='*70}")
print("Computing foreground quasar density (2° radius)...")
print(f"{'='*70}")

# Angular radius in 3D chord length: chord = 2*sin(θ/2)
radius_deg = 2.0
chord = 2 * np.sin(np.radians(radius_deg) / 2)

def fg_density(indices):
    """Count foreground quasars within radius for each object"""
    ra_rad = np.radians(ra_q[indices])
    dec_rad = np.radians(dec_q[indices])
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    xyz = np.column_stack([x, y, z])
    
    counts = tree_fg.query_ball_point(xyz, r=chord)
    return np.array([len(c) for c in counts])

fg_eb = fg_density(eb_idx)
fg_en = fg_density(en_idx)
fg_lh = fg_density(lh_idx)
fg_ln = fg_density(ln_idx)

print(f"\n  EARLY ZONE (z=0.8-0.95):")
print(f"  Breakers fg density: median = {np.median(fg_eb):.0f} (mean = {np.mean(fg_eb):.1f})")
print(f"  Normal fg density:   median = {np.median(fg_en):.0f} (mean = {np.mean(fg_en):.1f})")
_, p_early = stats.mannwhitneyu(fg_eb, fg_en, alternative='two-sided')
_, p_early_less = stats.mannwhitneyu(fg_eb, fg_en, alternative='less')
print(f"  MWU p (two-sided): {p_early:.4f}")
print(f"  MWU p (breakers < normal): {p_early_less:.4f}")

if np.median(fg_eb) < np.median(fg_en):
    print(f"  → Breakers have LESS foreground structure → VOID HYPOTHESIS SUPPORTED")
elif np.median(fg_eb) > np.median(fg_en):
    print(f"  → Breakers have MORE foreground structure → void hypothesis challenged")
else:
    print(f"  → Same foreground density")

print(f"\n  LATE ZONE (z=1.0-1.15):")
print(f"  Holders fg density: median = {np.median(fg_lh):.0f} (mean = {np.mean(fg_lh):.1f})")
print(f"  Normal fg density:  median = {np.median(fg_ln):.0f} (mean = {np.mean(fg_ln):.1f})")
_, p_late = stats.mannwhitneyu(fg_lh, fg_ln, alternative='greater')
print(f"  MWU p (holders > normal): {p_late:.4f}")

# ============================================================
# Quintile analysis — monotonic relationship?
# ============================================================
print(f"\n{'='*70}")
print("QUINTILE ANALYSIS — Is foreground density monotonically related to breaking?")
print(f"{'='*70}")

# For all early-zone objects, split by fg density quintiles
all_early = early_zone & valid_ra
ae_idx = valid_indices[all_early]
fg_all_early = fg_density(ae_idx)
ra_all_early = rank_agreement[all_early]

quintile_edges = np.percentile(fg_all_early, [0, 20, 40, 60, 80, 100])
print(f"\n  Fg density quintile edges: {[f'{e:.0f}' for e in quintile_edges]}")

print(f"\n  {'Fg density range':<20} {'Mean rank agr':>15} {'Breaker frac':>15} {'N':>8}")
print(f"  {'-'*20} {'-'*15} {'-'*15} {'-'*8}")

mean_ras = []
breaker_fracs = []

for i in range(5):
    lo, hi = quintile_edges[i], quintile_edges[i+1]
    mask = (fg_all_early >= lo) & (fg_all_early <= hi + 0.1)
    
    ra_bin = ra_all_early[mask]
    mean_ra = np.mean(ra_bin)
    breaker_frac = (ra_bin <= early_thresh).mean()
    
    mean_ras.append(mean_ra)
    breaker_fracs.append(breaker_frac)
    
    label = f"fg=[{lo:.0f},{hi:.0f}]"
    print(f"  {label:<20} {mean_ra:>15.3f} {100*breaker_frac:>14.1f}% {mask.sum():>8}")

# Monotonicity test
r_mono, p_mono = stats.spearmanr(range(5), breaker_fracs)
print(f"\n  Monotonicity (fg density vs breaker fraction): ρ = {r_mono:+.3f} (p = {p_mono:.3f})")

if r_mono < -0.5 and p_mono < 0.1:
    print(f"  ★ MORE foreground → FEWER breakers → STRUCTURE PROTECTS")
elif r_mono > 0.5 and p_mono < 0.1:
    print(f"  ★ MORE foreground → MORE breakers → structure DAMAGES")
else:
    print(f"  ~ No clear monotonic relationship")

# ============================================================
# Check intervening absorber data in DR16Q
# ============================================================
print(f"\n{'='*70}")
print("ABSORBER CHECK — Does DR16Q have intervening absorber data?")
print(f"{'='*70}")

# Check for absorber-related columns
abs_cols = [c for c in d.dtype.names if any(x in c.upper() for x in ['ABS', 'NABS', 'BAL', 'DLA', 'INTER', 'NQSO'])]
print(f"  Potentially relevant columns: {abs_cols}")

for col in abs_cols[:10]:
    try:
        vals = d[col]
        if vals.ndim == 0 or len(vals.shape) > 1:
            continue
        valid_vals = vals[np.isfinite(vals)] if np.issubdtype(vals.dtype, np.floating) else vals
        print(f"    {col}: dtype={vals.dtype}, sample={vals[:3]}")
    except:
        pass

# BAL = Broad Absorption Line quasars (internal, not intervening)
# Check for N_abs or similar
if 'BAL_FLAG' in d.dtype.names or 'BAL_PROB' in d.dtype.names:
    bal_col = 'BAL_FLAG' if 'BAL_FLAG' in d.dtype.names else 'BAL_PROB'
    print(f"\n  Found {bal_col} — checking if BAL quasars are overrepresented in breakers")
    
    bal_eb = d[bal_col][eb_idx]
    bal_en = d[bal_col][en_idx]
    
    # Handle different types
    if np.issubdtype(bal_eb.dtype, np.floating):
        bal_eb_v = bal_eb[np.isfinite(bal_eb)]
        bal_en_v = bal_en[np.isfinite(bal_en)]
        print(f"    Breakers mean {bal_col}: {np.mean(bal_eb_v):.3f}")
        print(f"    Normal mean {bal_col}:   {np.mean(bal_en_v):.3f}")

# ============================================================
# Try different angular radii for foreground density
# ============================================================
print(f"\n{'='*70}")
print("SCALE TEST — Does the foreground effect depend on angular scale?")
print(f"{'='*70}")

for radius in [0.5, 1.0, 2.0, 5.0, 10.0]:
    chord_r = 2 * np.sin(np.radians(radius) / 2)
    
    # Quick count for breakers and normals
    ra_rad_eb = np.radians(ra_q[eb_idx[:500]])  # subsample for speed
    dec_rad_eb = np.radians(dec_q[eb_idx[:500]])
    xyz_eb = np.column_stack([np.cos(dec_rad_eb)*np.cos(ra_rad_eb),
                               np.cos(dec_rad_eb)*np.sin(ra_rad_eb),
                               np.sin(dec_rad_eb)])
    
    ra_rad_en = np.radians(ra_q[en_idx[:500]])
    dec_rad_en = np.radians(dec_q[en_idx[:500]])
    xyz_en = np.column_stack([np.cos(dec_rad_en)*np.cos(ra_rad_en),
                               np.cos(dec_rad_en)*np.sin(ra_rad_en),
                               np.sin(dec_rad_en)])
    
    counts_eb = [len(c) for c in tree_fg.query_ball_point(xyz_eb, r=chord_r)]
    counts_en = [len(c) for c in tree_fg.query_ball_point(xyz_en, r=chord_r)]
    
    med_eb = np.median(counts_eb)
    med_en = np.median(counts_en)
    _, p = stats.mannwhitneyu(counts_eb, counts_en, alternative='two-sided')
    
    direction = "◄ LESS" if med_eb < med_en else "► MORE" if med_eb > med_en else "= SAME"
    sig = "★" if p < 0.05 else ""
    print(f"  r={radius:>4.1f}°: Breakers={med_eb:>5.0f}, Normal={med_en:>5.0f}, p={p:.3f} {direction} {sig}")

# Save
results = {
    'fg_density_breakers_median': float(np.median(fg_eb)),
    'fg_density_normals_median': float(np.median(fg_en)),
    'fg_density_holders_median': float(np.median(fg_lh)),
    'p_early': float(p_early),
    'monotonicity_rho': float(r_mono),
    'monotonicity_p': float(p_mono),
}

with open('results_void_structure/void_structure_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to results_void_structure/")
