#!/usr/bin/env python3
"""
P2 + P3 + SPATIAL AUTOCORRELATION BATTERY
Three tests in one run to conserve resources.

P2: Spatial autocorrelation of breaker map (are patches clustered or random?)
P3: Known structure alignment (do breaker patches coincide with known LSS features?)
P4: Density gradient test (does shear/gradient matter more than density itself?)
"""
import numpy as np
from scipy import stats
from scipy.spatial import cKDTree
import healpy as hp
import json, os, sys, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_battery', exist_ok=True)

print("="*60)
print("P2/P3/P4 TEST BATTERY")
print("="*60)
sys.stdout.flush()

# ============================================================
# LOAD DATA
# ============================================================
print("\n[1] Loading data...")
sys.stdout.flush()
d = np.genfromtxt('data/pantheon_plus.dat', names=True, dtype=None, encoding='utf-8')
ra = d['RA'].astype(float); dec = d['DEC'].astype(float)
z = d['zHD'].astype(float); c = d['c'].astype(float)
x1 = d['x1'].astype(float); mu = d['m_b_corr'].astype(float)
valid = np.isfinite(ra)&np.isfinite(dec)&np.isfinite(z)&np.isfinite(c)&(z>0.01)
ra,dec,z,c,x1,mu = ra[valid],dec[valid],z[valid],c[valid],x1[valid],mu[valid]
N = len(ra)

# Breaker classification
is_breaker = np.zeros(N, dtype=bool)
for i in range(N):
    local = np.abs(z - z[i]) < 0.05
    if local.sum() < 10: continue
    mu_c, sig_c = c[local].mean(), c[local].std()
    if sig_c > 0 and np.abs(c[i] - mu_c) > sig_c:
        is_breaker[i] = True
print(f"  {N} SNe, {is_breaker.sum()} breakers")

# ============================================================
# P2: SPATIAL AUTOCORRELATION
# ============================================================
print("\n[2] P2: Spatial autocorrelation of breakers...")
sys.stdout.flush()

# Convert to 3D unit sphere for angular distances
x3 = np.cos(np.radians(dec))*np.cos(np.radians(ra))
y3 = np.cos(np.radians(dec))*np.sin(np.radians(ra))
z3 = np.sin(np.radians(dec))
coords = np.column_stack([x3, y3, z3])
tree = cKDTree(coords)

# For angular bins, compute breaker-breaker pair fraction vs random
# Moran's I statistic for spatial autocorrelation
angle_bins = [2, 5, 10, 15, 20, 30]  # degrees
print(f"  Computing pair statistics at angular scales: {angle_bins}°")

results_p2 = {}
b_frac_global = is_breaker.mean()

for ang in angle_bins:
    r_cart = 2 * np.sin(np.radians(ang) / 2)
    pairs = tree.query_pairs(r_cart)
    if len(pairs) == 0:
        continue
    
    # For each pair, count breaker-breaker, breaker-normal, normal-normal
    bb, bn, nn = 0, 0, 0
    for i, j in pairs:
        bi, bj = is_breaker[i], is_breaker[j]
        if bi and bj: bb += 1
        elif bi or bj: bn += 1
        else: nn += 1
    
    total = bb + bn + nn
    # Expected breaker-breaker if random
    exp_bb = b_frac_global**2 * total
    # Excess
    excess = (bb - exp_bb) / exp_bb if exp_bb > 0 else 0
    
    results_p2[f'{ang}deg'] = {
        'pairs': total, 'bb': bb, 'bn': bn, 'nn': nn,
        'expected_bb': round(exp_bb, 1),
        'excess': round(excess, 4)
    }
    
    direction = "CLUSTERED" if excess > 0.05 else "ANTI-CLUSTERED" if excess < -0.05 else "RANDOM"
    print(f"  {ang}°: {total} pairs, BB={bb} (exp={exp_bb:.0f}), excess={excess:+.3f} → {direction}")

sys.stdout.flush()

# Moran's I on HEALPix breaker fraction map
NSIDE = 8
NPIX = hp.nside2npix(NSIDE)
theta = np.radians(90 - dec); phi = np.radians(ra)
pix = hp.ang2pix(NSIDE, theta, phi)

count_map = np.zeros(NPIX)
break_map = np.zeros(NPIX)
for i in range(N):
    count_map[pix[i]] += 1
    break_map[pix[i]] += is_breaker[i]

frac_map = np.full(NPIX, np.nan)
good = count_map >= 3
frac_map[good] = break_map[good] / count_map[good]

# Moran's I for HEALPix neighbors
good_idx = np.where(good)[0]
n_good = len(good_idx)
frac_good = frac_map[good_idx]
mean_f = frac_good.mean()
dev = frac_good - mean_f

# Weight matrix: HEALPix neighbors
W_sum = 0
numerator = 0
for ii, p in enumerate(good_idx):
    neighbors = hp.get_all_neighbours(NSIDE, p)
    neighbors = [n for n in neighbors if n >= 0 and n in good_idx]
    for nb in neighbors:
        jj = np.where(good_idx == nb)[0]
        if len(jj) > 0:
            jj = jj[0]
            numerator += dev[ii] * dev[jj]
            W_sum += 1

denom = np.sum(dev**2)
if W_sum > 0 and denom > 0:
    morans_I = (n_good / W_sum) * (numerator / denom)
else:
    morans_I = 0

# Expected under null
E_I = -1 / (n_good - 1)
print(f"\n  Moran's I = {morans_I:.4f} (expected under null: {E_I:.4f})")
print(f"  {'POSITIVE autocorrelation (breakers cluster)' if morans_I > 0.05 else 'NO significant autocorrelation' if abs(morans_I) < 0.05 else 'NEGATIVE (breakers anti-cluster)'}")

sys.stdout.flush()

# ============================================================
# P3: KNOWN STRUCTURE ALIGNMENT
# ============================================================
print("\n[3] P3: Alignment with known large-scale structures...")
sys.stdout.flush()

# Test breaker rate near known LSS features vs away
structures = {
    'Sloan_Great_Wall': {'ra': 195, 'dec': 5, 'z_range': (0.04, 0.12)},
    'CfA_Great_Wall': {'ra': 195, 'dec': 30, 'z_range': (0.02, 0.04)},
    'Shapley_SCL': {'ra': 202, 'dec': -31, 'z_range': (0.03, 0.06)},
    'Hercules_SCL': {'ra': 241, 'dec': 18, 'z_range': (0.03, 0.04)},
    'Bootes_Void': {'ra': 218, 'dec': 46, 'z_range': (0.04, 0.07)},
    'Dipole_Repeller': {'ra': 230, 'dec': -15, 'z_range': (0.04, 0.06)},
    'Cold_Spot': {'ra': 49.5, 'dec': -19.5, 'z_range': (0, 1)},  # CMB cold spot direction
}

from astropy.coordinates import SkyCoord
import astropy.units as u

sn_coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)

print(f"  {'Structure':<20s} {'Near(15°)':<10s} {'Break%':<8s} {'Far Break%':<10s} {'Δ':>8s}")
print("  " + "-"*60)

for name, info in structures.items():
    center = SkyCoord(ra=info['ra']*u.deg, dec=info['dec']*u.deg)
    sep = sn_coords.separation(center).degree
    near = sep < 15
    far = sep >= 15
    
    if near.sum() < 5: continue
    
    br_near = is_breaker[near].mean()
    br_far = is_breaker[far].mean()
    delta = br_near - br_far
    
    print(f"  {name:<20s} {near.sum():<10d} {100*br_near:<8.1f} {100*br_far:<10.1f} {100*delta:+8.1f}%")

# Test: are breakers preferentially behind/near walls vs voids?
# Use supergalactic plane as a proxy for the cosmic web spine
from astropy.coordinates import SkyCoord, Galactocentric
# Supergalactic coords
sg_coords = sn_coords.supergalactic
sg_b = sg_coords.sgb.degree  # supergalactic latitude

near_plane = np.abs(sg_b) < 15  # near supergalactic plane (web spine)
far_plane = np.abs(sg_b) >= 30  # away from plane (voids)

br_plane = is_breaker[near_plane].mean()
br_void = is_breaker[far_plane].mean()

print(f"\n  Supergalactic plane test:")
print(f"    Near plane (|SGB|<15°): {near_plane.sum()} SNe, {100*br_plane:.1f}% breakers")
print(f"    Far from plane (|SGB|>30°): {far_plane.sum()} SNe, {100*br_void:.1f}% breakers")
print(f"    Δ = {100*(br_plane - br_void):+.1f}%")

if near_plane.sum() > 20 and far_plane.sum() > 20:
    chi2_sg, p_sg = stats.chi2_contingency(
        [[is_breaker[near_plane].sum(), (~is_breaker[near_plane]).sum()],
         [is_breaker[far_plane].sum(), (~is_breaker[far_plane]).sum()]]
    )[:2]
    print(f"    χ² = {chi2_sg:.2f}, p = {p_sg:.4f}")

sys.stdout.flush()

# ============================================================
# P4: DENSITY GRADIENT TEST
# ============================================================
print("\n[4] P4: Local density gradient (angular clustering) test...")
sys.stdout.flush()

# For each SN, compute local SN angular density and density VARIANCE
# Idea: it's not how many neighbors, but how UNEVEN the distribution is
density_5deg = np.zeros(N)
density_var_5deg = np.zeros(N)

r5 = 2 * np.sin(np.radians(5) / 2)
neighbors_5 = tree.query_ball_point(coords, r5)

for i in range(N):
    nbrs = [j for j in neighbors_5[i] if j != i]
    density_5deg[i] = len(nbrs)
    
    if len(nbrs) >= 4:
        # Angular variance: how unevenly distributed are neighbors?
        nbr_coords = coords[nbrs]
        # Compute angular separations from SN i to all neighbors
        dots = np.clip(np.dot(nbr_coords, coords[i]), -1, 1)
        angles = np.arccos(dots)
        density_var_5deg[i] = np.std(angles)

# Correlations
has_nbrs = density_5deg >= 4
rho_dens, p_dens = stats.spearmanr(density_5deg[has_nbrs], is_breaker[has_nbrs].astype(float))
rho_grad, p_grad = stats.spearmanr(density_var_5deg[has_nbrs], is_breaker[has_nbrs].astype(float))

print(f"  SN angular density (5°) vs breaker: ρ = {rho_dens:+.4f}, p = {p_dens:.3e}")
print(f"  SN angular gradient (5°) vs breaker: ρ = {rho_grad:+.4f}, p = {p_grad:.3e}")
print(f"  Theory: gradient should matter MORE than density")

if abs(rho_grad) > abs(rho_dens):
    print(f"  ✅ Gradient > density ({abs(rho_grad):.4f} > {abs(rho_dens):.4f})")
else:
    print(f"  ❌ Density > gradient ({abs(rho_dens):.4f} > {abs(rho_grad):.4f})")

# Quintile analysis of density gradient
quintiles = np.percentile(density_var_5deg[has_nbrs], [20, 40, 60, 80])
print(f"\n  Breaker rate by angular gradient quintile:")
dv = density_var_5deg[has_nbrs]
ib = is_breaker[has_nbrs]
rates = []
for lo, hi, label in [(0, quintiles[0], 'Q1(smooth)'),
                       (quintiles[0], quintiles[1], 'Q2'),
                       (quintiles[1], quintiles[2], 'Q3'),
                       (quintiles[2], quintiles[3], 'Q4'),
                       (quintiles[3], 999, 'Q5(rough)')]:
    mask = (dv >= lo) & (dv < hi)
    if mask.sum() > 10:
        rate = ib[mask].mean()
        rates.append(rate)
        print(f"    {label}: {mask.sum()} SNe, {100*rate:.1f}% breakers")

if len(rates) >= 4:
    rho_q, _ = stats.spearmanr(range(len(rates)), rates)
    print(f"  Quintile trend: ρ = {rho_q:+.3f}")

sys.stdout.flush()

# ============================================================
# BONUS: GALACTIC vs EXTRAGALACTIC AXIS
# ============================================================
print("\n[5] Bonus: Galactic latitude dependence...")
sys.stdout.flush()

gal_coords = sn_coords.galactic
gal_b = gal_coords.b.degree

for bcut in [15, 20, 30, 45]:
    low_b = np.abs(gal_b) < bcut
    high_b = np.abs(gal_b) >= bcut
    if low_b.sum() > 20 and high_b.sum() > 20:
        br_low = is_breaker[low_b].mean()
        br_high = is_breaker[high_b].mean()
        print(f"  |b|<{bcut}°: {100*br_low:.1f}% breakers ({low_b.sum()} SNe) | |b|≥{bcut}°: {100*br_high:.1f}% ({high_b.sum()}) | Δ={100*(br_low-br_high):+.1f}%")

sys.stdout.flush()

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"\nP2 Spatial Autocorrelation:")
print(f"  Moran's I = {morans_I:.4f} (expect {E_I:.4f} under null)")
print(f"  Breaker pairs excess at small angles: see above")

print(f"\nP3 Structure Alignment:")
print(f"  Supergalactic plane: Δ = {100*(br_plane - br_void):+.1f}%")

print(f"\nP4 Density Gradient:")
print(f"  Density vs breaker: ρ = {rho_dens:+.4f}")
print(f"  Gradient vs breaker: ρ = {rho_grad:+.4f}")

# Save
results = {
    'p2_morans_I': round(morans_I, 4),
    'p2_morans_expected': round(E_I, 4),
    'p2_pair_excess': results_p2,
    'p3_supergalactic_delta': round(br_plane - br_void, 4),
    'p3_supergalactic_p': round(p_sg, 6) if 'p_sg' in dir() else None,
    'p4_density_rho': round(rho_dens, 4),
    'p4_gradient_rho': round(rho_grad, 4),
    'verdict_p2': 'CLUSTERED' if morans_I > 0.05 else 'RANDOM',
    'verdict_p4': 'GRADIENT>DENSITY' if abs(rho_grad) > abs(rho_dens) else 'DENSITY>GRADIENT'
}

with open('results_battery/battery_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nSaved to results_battery/battery_results.json")
print("="*60)
