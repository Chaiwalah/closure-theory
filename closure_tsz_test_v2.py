#!/usr/bin/env python3
"""
tSZ STRUCTURAL STATE TEST v2

Uses Planck PSZ2 cluster catalog (1,653 clusters) to construct a
sightline-based "ordered state" metric: proximity to virialized
massive clusters (high tSZ = ordered, bound space).

Two metrics per sightline:
1. N_clusters: number of PSZ2 clusters within angular radius
2. Y_integrated: sum of Compton-y of nearby clusters (weighted by distance)

Tests whether proximity to ordered structures predicts preservation
of observable correlations.
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
from scipy.spatial import cKDTree
import json, os, sys, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_tsz', exist_ok=True)

print("=" * 70)
print("tSZ STRUCTURAL STATE TEST v2")
print("PSZ2 cluster proximity vs observable degradation")
print("=" * 70)
sys.stdout.flush()

# ============================================================
# 1. LOAD PSZ2 CATALOG
# ============================================================
print("\n[1] Loading Planck PSZ2 catalog...")
sys.stdout.flush()

psz2 = fits.open('/root/clawd/projects/closure-theory/data/planck_psz2.fits')[1].data
ra_cl = psz2['RAJ2000'].astype(float)
dec_cl = psz2['DEJ2000'].astype(float)
y_cl = psz2['Y5R500'].astype(float)  # Compton-y within 5*R500
snr_cl = psz2['SNR'].astype(float)

# Filter valid clusters
valid_cl = np.isfinite(ra_cl) & np.isfinite(dec_cl) & np.isfinite(y_cl) & (y_cl > 0)
ra_cl = ra_cl[valid_cl]
dec_cl = dec_cl[valid_cl]
y_cl = y_cl[valid_cl]
snr_cl = snr_cl[valid_cl]

print(f"  Valid clusters: {len(ra_cl)}")
print(f"  Y5R500 range: [{y_cl.min():.2e}, {y_cl.max():.2e}]")
print(f"  SNR range: [{snr_cl.min():.1f}, {snr_cl.max():.1f}]")

# Build 3D tree on unit sphere for fast angular matching
coords_cl = SkyCoord(ra=ra_cl*u.degree, dec=dec_cl*u.degree, frame='icrs')
x_cl = np.cos(np.radians(dec_cl)) * np.cos(np.radians(ra_cl))
y_cl_3d = np.cos(np.radians(dec_cl)) * np.sin(np.radians(ra_cl))
z_cl_3d = np.sin(np.radians(dec_cl))
tree_cl = cKDTree(np.column_stack([x_cl, y_cl_3d, z_cl_3d]))

def compute_cluster_metrics(ra_arr, dec_arr, search_radius_deg=5.0):
    """Compute cluster proximity metrics for array of positions."""
    # Convert search radius to Cartesian distance on unit sphere
    search_rad_cart = 2 * np.sin(np.radians(search_radius_deg) / 2)
    
    x = np.cos(np.radians(dec_arr)) * np.cos(np.radians(ra_arr))
    y = np.cos(np.radians(dec_arr)) * np.sin(np.radians(ra_arr))
    z = np.sin(np.radians(dec_arr))
    
    points = np.column_stack([x, y, z])
    
    n_clusters = np.zeros(len(ra_arr), dtype=int)
    y_sum = np.zeros(len(ra_arr))
    y_nearest = np.full(len(ra_arr), np.nan)
    dist_nearest = np.full(len(ra_arr), np.inf)
    
    # Batch query
    neighbors = tree_cl.query_ball_point(points, search_rad_cart)
    
    for i, nbrs in enumerate(neighbors):
        n_clusters[i] = len(nbrs)
        if len(nbrs) > 0:
            y_sum[i] = y_cl[nbrs].sum()
            dists = np.sqrt(np.sum((points[i] - np.column_stack([x_cl[nbrs], y_cl_3d[nbrs], z_cl_3d[nbrs]]))**2, axis=1))
            nearest_idx = np.argmin(dists)
            y_nearest[i] = y_cl[nbrs[nearest_idx]]
            dist_nearest[i] = np.degrees(2 * np.arcsin(dists[nearest_idx] / 2))
    
    return n_clusters, y_sum, y_nearest, dist_nearest

sys.stdout.flush()

# ============================================================
# 2. SN Ia ANALYSIS
# ============================================================
print("\n[2] SN Ia cluster proximity analysis...")
sys.stdout.flush()

pan_data = np.genfromtxt('/root/clawd/projects/closure-theory/data/pantheon_plus.dat',
                         names=True, dtype=None, encoding='utf-8')

ra_sn = pan_data['RA'].astype(float)
dec_sn = pan_data['DEC'].astype(float)
z_sn = pan_data['zHD'].astype(float)
c_sn = pan_data['c'].astype(float)
x1_sn = pan_data['x1'].astype(float)
mu_corr = pan_data['m_b_corr'].astype(float)

valid_sn = np.isfinite(ra_sn) & np.isfinite(dec_sn) & np.isfinite(z_sn) & \
           np.isfinite(c_sn) & np.isfinite(x1_sn) & (z_sn > 0.01)

ra_v = ra_sn[valid_sn]
dec_v = dec_sn[valid_sn]
z_v = z_sn[valid_sn]
c_v = c_sn[valid_sn]
x1_v = x1_sn[valid_sn]
mu_v = mu_corr[valid_sn]

# Test multiple search radii
for radius in [2.0, 5.0, 10.0]:
    n_cl, y_sum, y_near, d_near = compute_cluster_metrics(ra_v, dec_v, radius)
    
    rho_nc, p_nc = stats.spearmanr(n_cl, c_v)
    rho_yc, p_yc = stats.spearmanr(y_sum, c_v)
    rho_nm, p_nm = stats.spearmanr(n_cl, mu_v)
    
    has_cluster = n_cl > 0
    print(f"\n  Radius = {radius}° (SNe near cluster: {has_cluster.sum()}/{len(ra_v)})")
    print(f"    N_cluster vs Color:  ρ = {rho_nc:+.4f}, p = {p_nc:.2e}")
    print(f"    Y_sum vs Color:      ρ = {rho_yc:+.4f}, p = {p_yc:.2e}")
    print(f"    N_cluster vs μ_corr: ρ = {rho_nm:+.4f}, p = {p_nm:.2e}")

# Use 5° for main analysis
n_cl_sn, y_sum_sn, _, d_near_sn = compute_cluster_metrics(ra_v, dec_v, 5.0)
sys.stdout.flush()

# Quintile by Y_sum
print(f"\n  Quintile analysis (5° radius, Y_sum):")
# Many will have y_sum=0, so split into "near cluster" vs "far from cluster"
near = y_sum_sn > 0
far = y_sum_sn == 0
print(f"    Near cluster: {near.sum()} SNe, mean(c)={c_v[near].mean():.4f}, std(c)={c_v[near].std():.4f}")
print(f"    Far from cl:  {far.sum()} SNe, mean(c)={c_v[far].mean():.4f}, std(c)={c_v[far].std():.4f}")

if near.sum() > 20 and far.sum() > 20:
    t_stat_c, t_p_c = stats.ttest_ind(c_v[near], c_v[far])
    t_stat_x, t_p_x = stats.ttest_ind(x1_v[near], x1_v[far])
    t_stat_m, t_p_m = stats.ttest_ind(mu_v[near], mu_v[far])
    print(f"    Color diff:   t={t_stat_c:.3f}, p={t_p_c:.3e}")
    print(f"    Stretch diff: t={t_stat_x:.3f}, p={t_p_x:.3e}")
    print(f"    μ_corr diff:  t={t_stat_m:.3f}, p={t_p_m:.3e}")
    
    # Variance comparison
    f_stat_c = c_v[far].var() / c_v[near].var()
    f_stat_x = x1_v[far].var() / x1_v[near].var()
    print(f"    Var ratio color (far/near):   {f_stat_c:.4f}")
    print(f"    Var ratio stretch (far/near): {f_stat_x:.4f}")
    print(f"    Theory: far > near (more scatter in disordered sightlines)")

sys.stdout.flush()

# ============================================================
# 3. QUASAR ANALYSIS — Breakers vs Holders
# ============================================================
print("\n[3] Quasar breaker/holder cluster proximity...")
sys.stdout.flush()

hdu_q = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
dq = hdu_q[1].data
z_q = dq['Z_DR16Q']
mg_ew = dq['MGII_BR'][:, 2]
hb_ew = dq['HBETA_BR'][:, 2]
ra_q = dq['RA']
dec_q = dq['DEC']

trans_mask = (z_q >= 0.75) & (z_q < 1.15) & np.isfinite(mg_ew) & np.isfinite(hb_ew) & \
             (mg_ew != 0) & (hb_ew != 0)

z_t = z_q[trans_mask]
mg_t = mg_ew[trans_mask]
hb_t = hb_ew[trans_mask]
ra_t = ra_q[trans_mask]
dec_t = dec_q[trans_mask]

print(f"  Transition zone quasars: {trans_mask.sum()}")

# Compute rank agreement
print("  Computing rank agreements...")
sys.stdout.flush()
window = 0.05
n_trans = len(z_t)
rank_agr = np.full(n_trans, np.nan)

# Sample for speed
np.random.seed(42)
sample_idx = np.random.choice(n_trans, min(50000, n_trans), replace=False)

for i in sample_idx:
    z_i = z_t[i]
    local = (z_t >= z_i - window) & (z_t < z_i + window)
    if local.sum() < 30:
        continue
    mg_pct = stats.percentileofscore(mg_t[local], mg_t[i]) / 100
    hb_pct = stats.percentileofscore(hb_t[local], hb_t[i]) / 100
    rank_agr[i] = 1 - abs(mg_pct - hb_pct)

valid_ra = np.isfinite(rank_agr)
early_zone = (z_t >= 0.80) & (z_t < 0.95)
late_zone = (z_t >= 1.00) & (z_t < 1.15)

ra_early_vals = rank_agr[early_zone & valid_ra]
ra_late_vals = rank_agr[late_zone & valid_ra]

if len(ra_early_vals) > 100 and len(ra_late_vals) > 100:
    early_thresh = np.percentile(ra_early_vals, 20)
    late_thresh = np.percentile(ra_late_vals, 80)
    
    breakers = early_zone & valid_ra & (rank_agr < early_thresh)
    holders = late_zone & valid_ra & (rank_agr > late_thresh)
    
    print(f"  Breakers: {breakers.sum()}, Holders: {holders.sum()}")
    
    # Cluster proximity for breakers vs holders
    n_cl_b, y_sum_b, _, d_near_b = compute_cluster_metrics(ra_t[breakers], dec_t[breakers], 5.0)
    n_cl_h, y_sum_h, _, d_near_h = compute_cluster_metrics(ra_t[holders], dec_t[holders], 5.0)
    
    print(f"\n  Breakers: mean N_cl={n_cl_b.mean():.3f}, mean Y_sum={y_sum_b.mean():.4e}, near_cl={( n_cl_b>0).sum()}/{breakers.sum()}")
    print(f"  Holders:  mean N_cl={n_cl_h.mean():.3f}, mean Y_sum={y_sum_h.mean():.4e}, near_cl={(n_cl_h>0).sum()}/{holders.sum()}")
    
    stat_u_n, p_u_n = stats.mannwhitneyu(n_cl_b, n_cl_h, alternative='two-sided')
    stat_u_y, p_u_y = stats.mannwhitneyu(y_sum_b, y_sum_h, alternative='two-sided')
    print(f"  N_cluster Mann-Whitney p = {p_u_n:.4e}")
    print(f"  Y_sum Mann-Whitney p = {p_u_y:.4e}")
    print(f"  Theory: holders should be nearer to clusters (ordered sightlines)")
    
    if n_cl_h.mean() > n_cl_b.mean():
        print(f"  ✅ Holders see more clusters")
    else:
        print(f"  ❌ Breakers see more clusters (wrong direction)")
        
sys.stdout.flush()

# ============================================================
# 4. FULL QUASAR COUPLING vs CLUSTER PROXIMITY
# ============================================================
print("\n[4] Quasar coupling vs cluster proximity (all transition zone)...")
sys.stdout.flush()

# Random sample of 50K quasars
sample_q = np.random.choice(n_trans, min(50000, n_trans), replace=False)
ra_qs = ra_t[sample_q]
dec_qs = dec_t[sample_q]
mg_qs = mg_t[sample_q]
hb_qs = hb_t[sample_q]

n_cl_q, y_sum_q, _, _ = compute_cluster_metrics(ra_qs, dec_qs, 5.0)

log_mg = np.log10(np.abs(mg_qs) + 1)
log_hb = np.log10(np.abs(hb_qs) + 1)

# Split by N_cluster
near_q = n_cl_q > 0
far_q = n_cl_q == 0

if near_q.sum() > 100 and far_q.sum() > 100:
    r_near, p_near = stats.spearmanr(log_mg[near_q], log_hb[near_q])
    r_far, p_far = stats.spearmanr(log_mg[far_q], log_hb[far_q])
    
    print(f"  Near cluster ({near_q.sum()}): ρ(MgII,Hβ) = {r_near:.4f}")
    print(f"  Far from cl  ({far_q.sum()}): ρ(MgII,Hβ) = {r_far:.4f}")
    print(f"  Δρ = {r_near - r_far:+.4f}")
    print(f"  Theory: near > far (ordered preserves coupling)")
    
    if r_near > r_far:
        print(f"  ✅ Correct direction")
    else:
        print(f"  ❌ Wrong direction")

# Quintile by N_clusters (for those with any)
if near_q.sum() > 200:
    q_bins_n = np.percentile(n_cl_q[near_q], [25, 50, 75])
    print(f"\n  Among cluster-proximate quasars:")
    print(f"  {'N_cl bin':<15s} {'N':>6s} {'ρ(MgII,Hβ)':>12s}")
    print("  " + "-" * 40)
    
    for lo, hi, label in [(0, q_bins_n[0], f'1-{int(q_bins_n[0])}'),
                           (q_bins_n[0], q_bins_n[1], f'{int(q_bins_n[0])}-{int(q_bins_n[1])}'),
                           (q_bins_n[1], q_bins_n[2], f'{int(q_bins_n[1])}-{int(q_bins_n[2])}'),
                           (q_bins_n[2], 999, f'>{int(q_bins_n[2])}')]:
        mask = near_q & (n_cl_q > lo) & (n_cl_q <= hi)
        if mask.sum() > 30:
            r, _ = stats.spearmanr(log_mg[mask], log_hb[mask])
            print(f"  {label:<15s} {mask.sum():6d} {r:12.4f}")

sys.stdout.flush()

# ============================================================
# 5. SHAPLEY CONCENTRATION TEST
# ============================================================
print("\n[5] Shapley Concentration test...")
sys.stdout.flush()

# Shapley center: RA~202°, Dec~-31° (l~312°, b~+30°), z~0.048
shapley_ra, shapley_dec = 202.0, -31.0
# Eridanus supervoid: RA~50°, Dec~-20°
eridanus_ra, eridanus_dec = 50.0, -20.0

# Find quasars near Shapley vs Eridanus sightlines (within 15°)
coords_qt = SkyCoord(ra=ra_qs*u.degree, dec=dec_qs*u.degree, frame='icrs')
shapley_coord = SkyCoord(ra=shapley_ra*u.degree, dec=shapley_dec*u.degree, frame='icrs')
eridanus_coord = SkyCoord(ra=eridanus_ra*u.degree, dec=eridanus_dec*u.degree, frame='icrs')

sep_shapley = coords_qt.separation(shapley_coord).degree
sep_eridanus = coords_qt.separation(eridanus_coord).degree

shapley_mask = sep_shapley < 15
eridanus_mask = sep_eridanus < 15

if shapley_mask.sum() > 50 and eridanus_mask.sum() > 50:
    r_shapley, _ = stats.spearmanr(log_mg[shapley_mask], log_hb[shapley_mask])
    r_eridanus, _ = stats.spearmanr(log_mg[eridanus_mask], log_hb[eridanus_mask])
    
    print(f"  Shapley sightlines  ({shapley_mask.sum()}): ρ(MgII,Hβ) = {r_shapley:.4f}")
    print(f"  Eridanus sightlines ({eridanus_mask.sum()}): ρ(MgII,Hβ) = {r_eridanus:.4f}")
    print(f"  Δρ = {r_shapley - r_eridanus:+.4f}")
    print(f"  Theory: Shapley (ordered) > Eridanus (disordered)")
    
    if r_shapley > r_eridanus:
        print(f"  ✅ Ordered sightlines preserve coupling better")
    else:
        print(f"  ❌ Wrong direction")
else:
    print(f"  Shapley: {shapley_mask.sum()}, Eridanus: {eridanus_mask.sum()} — insufficient")

sys.stdout.flush()

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY — tSZ STRUCTURAL STATE TEST v2")
print("=" * 70)

print(f"\nPSZ2 clusters used: {len(ra_cl)}")

print(f"\nSN Ia (N={len(ra_v)}):")
for radius in [2.0, 5.0, 10.0]:
    n_cl_r, y_sum_r, _, _ = compute_cluster_metrics(ra_v, dec_v, radius)
    rho_r, p_r = stats.spearmanr(n_cl_r, c_v)
    print(f"  {radius}°: N_cl vs color ρ = {rho_r:+.4f} (p={p_r:.2e})")

if near_q.sum() > 100 and far_q.sum() > 100:
    print(f"\nQuasar coupling (transition zone):")
    print(f"  Near cluster: ρ = {r_near:.4f}")
    print(f"  Far from cl:  ρ = {r_far:.4f}")
    print(f"  Δρ = {r_near - r_far:+.4f}")

if shapley_mask.sum() > 50 and eridanus_mask.sum() > 50:
    print(f"\nShapley vs Eridanus:")
    print(f"  Shapley:  ρ = {r_shapley:.4f}")
    print(f"  Eridanus: ρ = {r_eridanus:.4f}")
    print(f"  Δρ = {r_shapley - r_eridanus:+.4f}")

print("\n--- COMPARISON WITH KILL GRID ---")
print("  RM variance:    p = 0.97  (DEAD)")
print("  Planck κ:       p = 0.71  (DEAD)")
print("  Dust:           wrong direction (DEAD)")
print("  Foreground density: ρ = -0.20 (weak, not monotonic)")
if near_q.sum() > 100 and far_q.sum() > 100:
    delta_rho = r_near - r_far
    if abs(delta_rho) > 0.01:
        print(f"  Cluster proximity: Δρ = {delta_rho:+.4f} ← NEW")
    else:
        print(f"  Cluster proximity: Δρ = {delta_rho:+.4f} (weak)")

# Save results
results = {
    'psz2_clusters': int(len(ra_cl)),
    'sn_near_far': {
        'near_n': int(near.sum()),
        'far_n': int(far.sum()),
    },
}

with open('results_tsz/tsz_v2_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nResults saved to results_tsz/tsz_v2_results.json")
print("=" * 70)
