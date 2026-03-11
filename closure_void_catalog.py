#!/usr/bin/env python3
"""
COSMIC VOID CATALOG CROSS-MATCH

Uses Mao+2017 SDSS DR12 BOSS void catalog (VizieR J/ApJ/835/161)
to test whether early breakers preferentially align with void sightlines.

Also tests against known major cosmic structures (voids + superclusters)
as a sanity check.
"""

import numpy as np
from astropy.io import fits, votable
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
from scipy.spatial import cKDTree
import json, os, warnings, subprocess, io
warnings.filterwarnings('ignore')

os.makedirs('results_void_catalog', exist_ok=True)

print("=" * 70)
print("COSMIC VOID CATALOG CROSS-MATCH")
print("=" * 70)

# ============================================================
# Download Mao+2017 void catalog from VizieR
# ============================================================
print("\nDownloading Mao+2017 SDSS DR12 BOSS void catalog...")

result = subprocess.run([
    'curl', '-sL', '--max-time', '120',
    'https://vizier.cds.unistra.fr/viz-bin/votable?-source=J/ApJ/835/161/table1&-out.max=50000&-out.form=votable'
], capture_output=True)

vot = votable.parse(io.BytesIO(result.stdout))
table = vot.get_first_table().to_table()

void_ra = np.array(table['RAJ2000'], dtype=float)
void_dec = np.array(table['DEJ2000'], dtype=float)
void_radius = np.array(table['Reff'], dtype=float)  # Mpc/h effective radius
void_z = np.array(table['z'], dtype=float)

# Filter valid
valid_v = np.isfinite(void_ra) & np.isfinite(void_dec) & np.isfinite(void_radius) & np.isfinite(void_z)
void_ra = void_ra[valid_v]
void_dec = void_dec[valid_v]
void_radius = void_radius[valid_v]
void_z = void_z[valid_v]

print(f"  Loaded {len(void_ra)} voids")
print(f"  Redshift range: {void_z.min():.3f} - {void_z.max():.3f}")
print(f"  Radius range: {void_radius.min():.1f} - {void_radius.max():.1f} Mpc/h")
print(f"  Median radius: {np.median(void_radius):.1f} Mpc/h")

# Convert void angular sizes (approximate)
# At z=0.4, 1 Mpc/h ~ 0.3 degrees (very rough)
# void_angular_radius ~ void_radius * 0.3 / (void_z * 3000) degrees (rough)

# ============================================================
# Load quasar data
# ============================================================
print("\nLoading DR16Q quasars...")
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
ra_q = d['RA']
dec_q = d['DEC']

# ============================================================
# Recreate groups
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
# Build void sky density map
# ============================================================
print(f"\n{'='*70}")
print("Computing void proximity for each quasar sightline...")
print(f"{'='*70}")

# Use only foreground voids (z < 0.7, safely in front of our quasars)
fg_voids = void_z < 0.7
print(f"  Foreground voids (z < 0.7): {fg_voids.sum()}")

# Build KD-tree of void centers
void_ra_rad = np.radians(void_ra[fg_voids])
void_dec_rad = np.radians(void_dec[fg_voids])
void_x = np.cos(void_dec_rad) * np.cos(void_ra_rad)
void_y = np.cos(void_dec_rad) * np.sin(void_ra_rad)
void_z_cart = np.sin(void_dec_rad)
void_xyz = np.column_stack([void_x, void_y, void_z_cart])

tree_void = cKDTree(void_xyz)

# For each quasar, count nearby voids within different radii
def void_proximity(indices, radius_deg=5.0):
    """Count voids within angular radius of each quasar."""
    chord = 2 * np.sin(np.radians(radius_deg) / 2)
    
    ra_rad = np.radians(ra_q[indices])
    dec_rad = np.radians(dec_q[indices])
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    xyz = np.column_stack([x, y, z])
    
    counts = tree_void.query_ball_point(xyz, r=chord)
    return np.array([len(c) for c in counts])

# Also compute distance to nearest void
def nearest_void_dist(indices):
    """Angular distance to nearest void center (degrees)."""
    ra_rad = np.radians(ra_q[indices])
    dec_rad = np.radians(dec_q[indices])
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    xyz = np.column_stack([x, y, z])
    
    dists, _ = tree_void.query(xyz)
    # Convert chord distance to angle
    return np.degrees(2 * np.arcsin(dists / 2))

# ============================================================
# TEST 1: Void count within various radii
# ============================================================
print(f"\n{'='*70}")
print("TEST 1: Void count near quasar sightlines")
print(f"{'='*70}")

for radius in [2.0, 5.0, 10.0]:
    v_eb = void_proximity(eb_idx, radius)
    v_en = void_proximity(en_idx, radius)
    v_lh = void_proximity(lh_idx, radius)
    v_ln = void_proximity(ln_idx, radius)
    
    _, p = stats.mannwhitneyu(v_eb, v_en, alternative='two-sided')
    _, p_greater = stats.mannwhitneyu(v_eb, v_en, alternative='greater')
    
    direction = "MORE" if np.median(v_eb) > np.median(v_en) else "LESS" if np.median(v_eb) < np.median(v_en) else "SAME"
    sig = "★" if p < 0.05 else ""
    
    print(f"\n  r={radius}°:")
    print(f"    Breakers: median={np.median(v_eb):.0f}, mean={np.mean(v_eb):.1f}")
    print(f"    Normal:   median={np.median(v_en):.0f}, mean={np.mean(v_en):.1f}")
    print(f"    MWU p (two-sided): {p:.4f} → Breakers have {direction} void neighbors {sig}")
    print(f"    MWU p (breakers > normal): {p_greater:.4f}")

# ============================================================
# TEST 2: Distance to nearest void
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: Distance to nearest void center")
print(f"{'='*70}")

dist_eb = nearest_void_dist(eb_idx)
dist_en = nearest_void_dist(en_idx)
dist_lh = nearest_void_dist(lh_idx)
dist_ln = nearest_void_dist(ln_idx)

print(f"\n  EARLY ZONE:")
print(f"  Breakers nearest void: median = {np.median(dist_eb):.2f}°")
print(f"  Normal nearest void:   median = {np.median(dist_en):.2f}°")
_, p_dist = stats.mannwhitneyu(dist_eb, dist_en, alternative='two-sided')
_, p_dist_less = stats.mannwhitneyu(dist_eb, dist_en, alternative='less')
print(f"  MWU p (two-sided): {p_dist:.4f}")
print(f"  MWU p (breakers closer to voids): {p_dist_less:.4f}")

print(f"\n  LATE ZONE:")
print(f"  Holders nearest void: median = {np.median(dist_lh):.2f}°")
print(f"  Normal nearest void:  median = {np.median(dist_ln):.2f}°")
_, p_dist_lh = stats.mannwhitneyu(dist_lh, dist_ln, alternative='two-sided')
print(f"  MWU p (two-sided): {p_dist_lh:.4f}")

# ============================================================
# TEST 3: Quintile analysis — void count vs breaker fraction
# ============================================================
print(f"\n{'='*70}")
print("QUINTILE ANALYSIS — Void proximity vs breaker fraction")
print(f"{'='*70}")

# Use 5° radius for quintile analysis
all_early = early_zone & valid_ra
ae_idx = valid_indices[all_early]
v_all_early = void_proximity(ae_idx, 5.0)
ra_all_early = rank_agreement[all_early]

quintile_edges = np.percentile(v_all_early, [0, 20, 40, 60, 80, 100])
print(f"\n  Void count (5°) quintile edges: {[f'{e:.0f}' for e in quintile_edges]}")

print(f"\n  {'Void count range':<20} {'Mean rank agr':>15} {'Breaker frac':>15} {'N':>8}")
print(f"  {'-'*20} {'-'*15} {'-'*15} {'-'*8}")

breaker_fracs = []
for i in range(5):
    lo, hi = quintile_edges[i], quintile_edges[i+1]
    if i < 4:
        m = (v_all_early >= lo) & (v_all_early < hi)
    else:
        m = (v_all_early >= lo) & (v_all_early <= hi)
    
    if m.sum() == 0:
        breaker_fracs.append(np.nan)
        continue
    
    ra_bin = ra_all_early[m]
    mean_ra = np.mean(ra_bin)
    breaker_frac = (ra_bin <= early_thresh).mean()
    breaker_fracs.append(breaker_frac)
    
    label = f"voids=[{lo:.0f},{hi:.0f}]"
    print(f"  {label:<20} {mean_ra:>15.3f} {100*breaker_frac:>14.1f}% {m.sum():>8}")

valid_bf = [x for x in breaker_fracs if not np.isnan(x)]
r_mono, p_mono = stats.spearmanr(range(len(valid_bf)), valid_bf)
print(f"\n  Monotonicity: ρ = {r_mono:+.3f} (p = {p_mono:.3f})")

# ============================================================
# TEST 4: Known major structures check
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: Known cosmic structures — breaker excess/deficit")
print(f"{'='*70}")

known_structures = [
    # Name, RA, DEC, radius_deg, type
    ("Boötes Void", 218.0, 46.0, 15.0, "void"),
    ("Eridanus Supervoid", 51.0, -20.0, 10.0, "void"),
    ("KBC Void (Local)", 180.0, 0.0, 30.0, "void"),
    ("Shapley Concentration", 202.0, -31.0, 10.0, "supercluster"),
    ("Sloan Great Wall", 195.0, 7.5, 15.0, "supercluster"),
    ("Coma Cluster", 194.9, 27.98, 5.0, "supercluster"),
    ("Corona Borealis SC", 229.0, 27.0, 5.0, "supercluster"),
    ("Hercules SC", 242.0, 17.0, 5.0, "supercluster"),
]

for name, ra_s, dec_s, r_s, stype in known_structures:
    # Count breakers and normals in this direction
    sc = SkyCoord(ra=ra_s*u.degree, dec=dec_s*u.degree)
    
    qc_eb = SkyCoord(ra=ra_q[eb_idx]*u.degree, dec=dec_q[eb_idx]*u.degree)
    qc_en = SkyCoord(ra=ra_q[en_idx]*u.degree, dec=dec_q[en_idx]*u.degree)
    
    sep_eb = sc.separation(qc_eb).degree
    sep_en = sc.separation(qc_en).degree
    
    n_eb = (sep_eb < r_s).sum()
    n_en = (sep_en < r_s).sum()
    
    # Normalize by total
    frac_eb = n_eb / len(eb_idx)
    frac_en = n_en / len(en_idx)
    
    if frac_en > 0:
        excess = (frac_eb - frac_en) / frac_en * 100
    else:
        excess = 0
    
    marker = "🕳️" if stype == "void" else "🌌"
    sig = ""
    if abs(excess) > 20:
        sig = " ★"
    
    print(f"  {marker} {name:25s}: EB={n_eb:>4}, EN={n_en:>4}, excess={excess:+.1f}%{sig}")

# Save results
results = {
    'n_voids_total': int(len(void_ra)),
    'n_voids_foreground': int(fg_voids.sum()),
    'dist_breakers_median': float(np.median(dist_eb)),
    'dist_normals_median': float(np.median(dist_en)),
    'p_distance': float(p_dist),
    'monotonicity_rho': float(r_mono),
    'monotonicity_p': float(p_mono),
}

with open('results_void_catalog/void_catalog_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to results_void_catalog/")
