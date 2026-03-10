#!/usr/bin/env python3
"""
closure_mixture_decomp.py — Law of Total Variance Decomposition
================================================================

Var(Z|t) = E[Var(Z|C,t)] + Var(E[Z|C,t])
         = within-channel    + between-channel

Try every available split variable C to find the latent channel:
1. Host mass (high/low)
2. Survey (IDSURVEY)
3. Color sign (red/blue)
4. Stretch sign (fast/slow)
5. MW extinction (dusty/clean)
6. Galactic latitude (high/low |b|)

For each split: compute within-group and between-group variance.
If a split dramatically reduces within-group scatter → found C.

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr
import json
from pathlib import Path

np.random.seed(42)

# Load data
data = []
with open('data/pantheon_plus.dat', 'r') as f:
    header = f.readline().strip().split()
    for line in f:
        parts = line.strip().split()
        if len(parts) >= len(header):
            row = {}
            for i, h in enumerate(header):
                try: row[h] = float(parts[i])
                except: row[h] = parts[i]
            data.append(row)

seen = set()
sne = []
for row in data:
    z = row.get('zHD', 0)
    if z > 0.01 and z < 2.5 and abs(row.get('c',-999)) < 0.3 and abs(row.get('x1',-999)) < 3:
        key = (row['CID'], round(z, 5))
        if key not in seen:
            seen.add(key)
            sne.append(row)

# Add galactic latitude
def gal_b(ra, dec):
    ra_r, dec_r = np.radians(ra), np.radians(dec)
    ra_ngp, dec_ngp = np.radians(192.8595), np.radians(27.1284)
    sin_b = np.sin(dec_r)*np.sin(dec_ngp) + np.cos(dec_r)*np.cos(dec_ngp)*np.cos(ra_r-ra_ngp)
    return np.degrees(np.arcsin(np.clip(sin_b, -1, 1)))

for s in sne:
    s['abs_gal_b'] = abs(gal_b(s['RA'], s['DEC']))

print(f"Loaded {len(sne)} unique SNe Ia\n")

z_edges = [0.01, 0.1, 0.2, 0.3, 0.5, 0.82, 2.5]

# ============================================================
# For each observable Z and each split C, compute decomposition
# ============================================================

observables = {
    'mB': lambda s: s['mB'],
    'x1': lambda s: s['x1'],
    'c': lambda s: s['c'],
    'm_b_corr': lambda s: s['m_b_corr'],
}

# Split functions — return True for "group A"
def make_median_split(key):
    vals = [s.get(key, 0) for s in sne if s.get(key, -999) > -990]
    med = np.median(vals)
    return lambda s: s.get(key, -999) <= med, f"{key}≤{med:.2f}"

splits = []

# Host mass
mass_vals = [s['HOST_LOGMASS'] for s in sne if s.get('HOST_LOGMASS', 0) > 5]
mass_med = np.median(mass_vals)
splits.append(("Host mass", lambda s: s.get('HOST_LOGMASS', 0) > 5 and s['HOST_LOGMASS'] <= mass_med))

# Color sign
splits.append(("Color sign (c<0)", lambda s: s['c'] < 0))

# Stretch sign  
splits.append(("Stretch sign (x1<0)", lambda s: s['x1'] < 0))

# MW extinction
mwebv_med = np.median([s['MWEBV'] for s in sne])
splits.append(("MW dust (low)", lambda s: s['MWEBV'] <= mwebv_med))

# Galactic latitude
b_med = np.median([s['abs_gal_b'] for s in sne])
splits.append(("|b| (high lat)", lambda s: s['abs_gal_b'] > b_med))

# Survey — use most common survey as split
surveys = [s['IDSURVEY'] for s in sne]
from collections import Counter
survey_counts = Counter(surveys)
top_survey = survey_counts.most_common(1)[0][0]
splits.append((f"Survey={top_survey}", lambda s: s['IDSURVEY'] == top_survey))

print("=" * 90)
print("LAW OF TOTAL VARIANCE DECOMPOSITION")
print("Var(Z) = E[Var(Z|C)] + Var(E[Z|C]) = within + between")
print("=" * 90)

# For each z-bin and each split, decompose variance of x1
# (x1 showed strongest contraction)

print(f"\n  Observable: x1 (stretch — strongest contraction signal)")
print(f"\n{'z_range':<12} {'Split C':<22} {'Var_tot':>8} {'Within':>8} {'Between':>8} {'%Between':>9} {'p(A)':>6}")
print("-" * 80)

best_splits = {}  # track which split explains most between-variance

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 20:
        continue
    
    z_label = f"[{z_lo:.2f},{z_hi:.2f})"
    x1_all = np.array([s['x1'] for s in bin_sne])
    var_total = np.var(x1_all)
    
    best_between_frac = 0
    best_split_name = ""
    
    for split_name, split_fn in splits:
        grp_a = [s for s in bin_sne if split_fn(s)]
        grp_b = [s for s in bin_sne if not split_fn(s)]
        
        if len(grp_a) < 5 or len(grp_b) < 5:
            continue
        
        x1_a = np.array([s['x1'] for s in grp_a])
        x1_b = np.array([s['x1'] for s in grp_b])
        
        p_a = len(grp_a) / len(bin_sne)
        p_b = 1 - p_a
        
        # Within-channel: E[Var(Z|C)]
        within = p_a * np.var(x1_a) + p_b * np.var(x1_b)
        
        # Between-channel: Var(E[Z|C])
        between = p_a * p_b * (np.mean(x1_a) - np.mean(x1_b))**2
        
        # Check: within + between should ≈ var_total
        pct_between = between / var_total * 100 if var_total > 0 else 0
        
        print(f"{z_label:<12} {split_name:<22} {var_total:>8.4f} {within:>8.4f} {between:>8.4f} {pct_between:>8.1f}% {p_a:>6.2f}")
        
        if pct_between > best_between_frac:
            best_between_frac = pct_between
            best_split_name = split_name
    
    best_splits[z_label] = (best_split_name, best_between_frac)
    print(f"{'':<12} {'→ BEST: '+best_split_name:<22} {'':<8} {'':<8} {'':<8} {best_between_frac:>8.1f}%")
    print()

# ============================================================
# SAME FOR COLOR
# ============================================================
print(f"\n  Observable: c (color — stays flat)")
print(f"\n{'z_range':<12} {'Split C':<22} {'Var_tot':>8} {'Within':>8} {'Between':>8} {'%Between':>9}")
print("-" * 75)

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 20:
        continue
    
    z_label = f"[{z_lo:.2f},{z_hi:.2f})"
    c_all = np.array([s['c'] for s in bin_sne])
    var_total = np.var(c_all)
    
    for split_name, split_fn in splits[:3]:  # just top 3 splits
        grp_a = [s for s in bin_sne if split_fn(s)]
        grp_b = [s for s in bin_sne if not split_fn(s)]
        
        if len(grp_a) < 5 or len(grp_b) < 5:
            continue
        
        c_a = np.array([s['c'] for s in grp_a])
        c_b = np.array([s['c'] for s in grp_b])
        
        p_a = len(grp_a) / len(bin_sne)
        p_b = 1 - p_a
        
        within = p_a * np.var(c_a) + p_b * np.var(c_b)
        between = p_a * p_b * (np.mean(c_a) - np.mean(c_b))**2
        pct = between / var_total * 100 if var_total > 0 else 0
        
        print(f"{z_label:<12} {split_name:<22} {var_total:>8.6f} {within:>8.6f} {between:>8.6f} {pct:>8.1f}%")

# ============================================================
# DOES THE BETWEEN-CHANNEL FRACTION CHANGE WITH Z?
# ============================================================
print(f"\n\n" + "=" * 90)
print("BETWEEN-CHANNEL FRACTION vs z — Does the mixture shift?")
print("=" * 90)

# For the best split (host mass), track between-fraction vs z
print(f"\n  Using Host Mass split:")
print(f"  {'z':>6} {'Var_tot':>8} {'Within':>8} {'Between':>8} {'%Between':>9} {'N_A':>5} {'N_B':>5}")
print(f"  " + "-" * 55)

z_fine_edges = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.82, 1.0, 2.5]
z_btwn = []
pct_btwn = []

for i in range(len(z_fine_edges) - 1):
    z_lo, z_hi = z_fine_edges[i], z_fine_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 15:
        continue
    
    grp_a = [s for s in bin_sne if s.get('HOST_LOGMASS', 0) > 5 and s['HOST_LOGMASS'] <= mass_med]
    grp_b = [s for s in bin_sne if s.get('HOST_LOGMASS', 0) > 5 and s['HOST_LOGMASS'] > mass_med]
    
    if len(grp_a) < 3 or len(grp_b) < 3:
        continue
    
    x1_all = np.array([s['x1'] for s in bin_sne])
    x1_a = np.array([s['x1'] for s in grp_a])
    x1_b = np.array([s['x1'] for s in grp_b])
    
    var_tot = np.var(x1_all)
    p_a = len(grp_a) / (len(grp_a) + len(grp_b))
    p_b = 1 - p_a
    
    within = p_a * np.var(x1_a) + p_b * np.var(x1_b)
    between = p_a * p_b * (np.mean(x1_a) - np.mean(x1_b))**2
    pct = between / var_tot * 100 if var_tot > 0 else 0
    
    zc = np.mean([s['zHD'] for s in bin_sne])
    z_btwn.append(zc)
    pct_btwn.append(pct)
    
    print(f"  {zc:>6.2f} {var_tot:>8.4f} {within:>8.4f} {between:>8.4f} {pct:>8.1f}% {len(grp_a):>5} {len(grp_b):>5}")

if len(z_btwn) >= 4:
    rho_btwn, p_btwn = spearmanr(z_btwn, pct_btwn)
    print(f"\n  %Between vs z: ρ = {rho_btwn:+.3f}, p = {p_btwn:.4f}")
    
    if abs(rho_btwn) > 0.4:
        print(f"  → Between-channel fraction CHANGES with z → mixture IS shifting")
    else:
        print(f"  → Between-channel fraction STABLE → not a shifting mixture")

# ============================================================
# GAUSSIAN MIXTURE MODEL (2 components)
# ============================================================
print(f"\n\n" + "=" * 90)
print("GAUSSIAN MIXTURE MODEL: Find latent C from data alone")
print("=" * 90)

try:
    from sklearn.mixture import GaussianMixture
    
    print(f"\n  Fitting 2-component GMM to (x1, c) per z-bin:")
    print(f"  {'z':>6} {'N':>4} {'μ1_x1':>7} {'μ2_x1':>7} {'Δμ_x1':>7} {'σ1_x1':>7} {'σ2_x1':>7} {'p1':>5} {'BIC_2-BIC_1':>12}")
    print(f"  " + "-" * 75)
    
    z_gmm = []
    bic_diff = []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
        if len(bin_sne) < 30:
            continue
        
        X = np.column_stack([[s['x1'] for s in bin_sne], [s['c'] for s in bin_sne]])
        
        gmm1 = GaussianMixture(n_components=1, random_state=42).fit(X)
        gmm2 = GaussianMixture(n_components=2, random_state=42).fit(X)
        
        bic1 = gmm1.bic(X)
        bic2 = gmm2.bic(X)
        delta_bic = bic2 - bic1  # negative = 2 components better
        
        means = gmm2.means_
        covs = gmm2.covariances_
        weights = gmm2.weights_
        
        # Sort by x1 mean
        idx = np.argsort(means[:, 0])
        
        zc = np.mean([s['zHD'] for s in bin_sne])
        z_gmm.append(zc)
        bic_diff.append(delta_bic)
        
        sig1 = np.sqrt(covs[idx[0]][0,0])
        sig2 = np.sqrt(covs[idx[1]][0,0])
        delta_mu = means[idx[1],0] - means[idx[0],0]
        
        better = "★2>1" if delta_bic < -10 else ""
        
        print(f"  {zc:>6.2f} {len(bin_sne):>4} {means[idx[0],0]:>7.2f} {means[idx[1],0]:>7.2f} {delta_mu:>7.2f} {sig1:>7.2f} {sig2:>7.2f} {weights[idx[0]]:>5.2f} {delta_bic:>12.1f} {better}")
    
    if len(z_gmm) >= 3:
        rho_bic, p_bic_val = spearmanr(z_gmm, bic_diff)
        print(f"\n  ΔBIC vs z: ρ = {rho_bic:+.3f}, p = {p_bic_val:.4f}")
        print(f"  (negative ΔBIC = 2 components preferred)")
        
        n_prefer_2 = sum(1 for b in bic_diff if b < -10)
        print(f"  Bins preferring 2 components: {n_prefer_2}/{len(bic_diff)}")

except ImportError:
    print("  sklearn not available — skipping GMM")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_mixture_decomp")
results_dir.mkdir(exist_ok=True)

results = {'best_splits': {k: {'split': v[0], 'pct_between': v[1]} for k, v in best_splits.items()}}
with open(results_dir / "mixture_decomp.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {results_dir}/")
