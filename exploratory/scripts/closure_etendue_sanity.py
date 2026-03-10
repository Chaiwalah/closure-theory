#!/usr/bin/env python3
"""
closure_etendue_sanity.py — Étendue Artifact Kill Test
========================================================

GPT's challenge: Var × Peak might be conserved by construction if 
Peak_Density uses adaptive bandwidth KDE (h ∝ σ → peak ∝ 1/σ → Var×Peak ∝ σ).

Test with:
1. Fixed-bin histogram (same bins across all z)
2. Fixed-bandwidth KDE (same h across all z)  
3. Fitted Gaussian mixture peak
4. Collision probability ∫f² (Rényi-2, estimator-independent)
5. Shannon entropy H = -∫f ln f

If étendue survives all methods → real conservation law
If only adaptive KDE → artifact

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr, gaussian_kde
from scipy.optimize import minimize_scalar
from astropy.io import fits
from numpy.polynomial import polynomial as P
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

print("Loading DR16Q catalog...")
f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z_all = d['Z_DR16Q']
lbol = d['LOGLBOL']

def compute_peak_methods(x, fixed_bins, fixed_bw):
    """Compute peak density using multiple methods."""
    results = {}
    
    # Method 1: Fixed-bin histogram
    counts, edges = np.histogram(x, bins=fixed_bins, density=True)
    results['peak_fixedbin'] = np.max(counts)
    
    # Method 2: Fixed-bandwidth KDE
    try:
        kde_fixed = gaussian_kde(x, bw_method=fixed_bw / np.std(x))
        grid = np.linspace(np.percentile(x, 1), np.percentile(x, 99), 200)
        results['peak_fixedbw'] = np.max(kde_fixed.evaluate(grid))
    except:
        results['peak_fixedbw'] = np.nan
    
    # Method 3: Adaptive KDE (Silverman — the one we used before)
    try:
        kde_adapt = gaussian_kde(x)
        grid = np.linspace(np.percentile(x, 1), np.percentile(x, 99), 200)
        results['peak_adaptive'] = np.max(kde_adapt.evaluate(grid))
    except:
        results['peak_adaptive'] = np.nan
    
    # Method 4: Collision probability ∫f² (Rényi-2)
    try:
        kde = gaussian_kde(x)
        grid = np.linspace(np.percentile(x, 0.5), np.percentile(x, 99.5), 500)
        f_vals = kde.evaluate(grid)
        dx = grid[1] - grid[0]
        results['collision_prob'] = np.sum(f_vals**2) * dx
    except:
        results['collision_prob'] = np.nan
    
    # Method 5: Shannon entropy
    try:
        kde = gaussian_kde(x)
        grid = np.linspace(np.percentile(x, 0.5), np.percentile(x, 99.5), 500)
        f_vals = kde.evaluate(grid)
        dx = grid[1] - grid[0]
        f_pos = f_vals[f_vals > 1e-10]
        results['shannon_H'] = -np.sum(f_pos * np.log(f_pos)) * dx
    except:
        results['shannon_H'] = np.nan
    
    # Method 6: Fitted single Gaussian peak (1/sqrt(2πσ²))
    results['peak_gaussian'] = 1.0 / (np.std(x) * np.sqrt(2 * np.pi))
    
    return results


# ============================================================
# TEST ON MgII_BR (the one with p=1.000)
# ============================================================
print(f"\n{'=' * 130}")
print("MgII_BR: ÉTENDUE SANITY CHECK — Is Var × Peak conserved by construction?")
print("=" * 130)

mgii = d['MGII_BR']
mask = ((z_all >= 0.4) & (z_all < 2.5) &
        (mgii[:, 4] > 100) & (mgii[:, 4] < 30000) &
        (mgii[:, 2] > 0) & (mgii[:, 2] < 5000) &
        np.isfinite(mgii[:, 2]) & np.isfinite(mgii[:, 4]) &
        (lbol > 40) & (lbol < 50))

z_v = z_all[mask]
log_ew = np.log10(mgii[mask, 2])
lbol_v = lbol[mask]
coef = P.polyfit(lbol_v, log_ew, 1)
log_ew_r = log_ew - P.polyval(lbol_v, coef)

# Fixed bins: use global range, 50 bins
global_min, global_max = np.percentile(log_ew_r, 0.5), np.percentile(log_ew_r, 99.5)
fixed_bins = np.linspace(global_min, global_max, 51)

# Fixed bandwidth: use global std * Silverman factor
global_std = np.std(log_ew_r)
N_total = len(log_ew_r)
fixed_bw = global_std * (4 / (3 * N_total)) ** 0.2  # Silverman on full sample

n_zbins = 8
z_edges = np.percentile(z_v, np.linspace(0, 100, n_zbins + 1))
z_edges = np.unique(np.round(z_edges, 3))

# Collect results
methods = ['peak_fixedbin', 'peak_fixedbw', 'peak_adaptive', 'peak_gaussian', 'collision_prob']
method_labels = ['Fixed-bin hist', 'Fixed-bw KDE', 'Adaptive KDE', 'Gaussian peak', 'Collision ∫f²']

print(f"\n  {'z':>6} {'N':>6} {'Var':>8} | ", end='')
for ml in method_labels:
    print(f"{'Peak('+ml+')':>16} {'Var×P':>8} | ", end='')
print()
print("  " + "-" * 160)

z_list = []
etendue_results = {m: [] for m in methods}
var_list = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bm = (z_v >= z_lo) & (z_v < z_hi)
    n = np.sum(bm)
    if n < 500:
        continue
    
    zc = np.mean(z_v[bm])
    x = log_ew_r[bm]
    var = np.var(x)
    
    peaks = compute_peak_methods(x, fixed_bins, fixed_bw)
    
    z_list.append(zc)
    var_list.append(var)
    
    print(f"  {zc:>6.3f} {n:>6} {var:>8.5f} | ", end='')
    for m in methods:
        p = peaks[m]
        et = var * p if not np.isnan(p) else np.nan
        etendue_results[m].append(et)
        print(f"{p:>16.4f} {et:>8.4f} | ", end='')
    print()

# Correlations
print(f"\n  {'Method':>20} | {'ρ(étendue,z)':>14} {'p':>8} | {'CV(étendue)':>12} | Verdict")
print("  " + "-" * 80)

for m, ml in zip(methods, method_labels):
    vals = etendue_results[m]
    if len(vals) >= 3 and not any(np.isnan(vals)):
        rho, p = spearmanr(z_list, vals)
        cv = np.std(vals) / np.mean(vals) * 100
        verdict = "🔥 CONSERVED" if abs(rho) < 0.3 else "ARTIFACT" if abs(rho) > 0.5 else "WEAK"
        print(f"  {ml:>20} | {rho:>+14.3f} {p:>8.4f} | {cv:>11.1f}% | {verdict}")

# Also test Var × exp(-2H) (GPT's suggestion)
print(f"\n  GPT's suggested invariants:")
shannon_list = []
renyi2_list = []
fisher_list = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bm = (z_v >= z_lo) & (z_v < z_hi)
    if np.sum(bm) < 500:
        continue
    x = log_ew_r[bm]
    
    # Shannon
    kde = gaussian_kde(x)
    grid = np.linspace(np.percentile(x, 0.5), np.percentile(x, 99.5), 500)
    f_vals = kde.evaluate(grid)
    dx = grid[1] - grid[0]
    f_pos = f_vals[f_vals > 1e-10]
    H = -np.sum(f_pos * np.log(f_pos)) * dx
    shannon_list.append(H)
    
    # Collision prob (Rényi-2)
    cp = np.sum(f_vals**2) * dx
    renyi2_list.append(cp)

# Var × e^{-2H}
var_arr = np.array(var_list)
inv1 = var_arr * np.exp(-2 * np.array(shannon_list))
inv2 = var_arr * np.array(renyi2_list)

rho1, p1 = spearmanr(z_list, inv1)
rho2, p2 = spearmanr(z_list, inv2)
cv1 = np.std(inv1) / np.mean(inv1) * 100
cv2 = np.std(inv2) / np.mean(inv2) * 100

print(f"  {'Var × e^(-2H)':>20} | {rho1:>+14.3f} {p1:>8.4f} | {cv1:>11.1f}% | {'🔥 CONSERVED' if abs(rho1) < 0.3 else 'NOT CONSERVED'}")
print(f"  {'Var × ∫f²':>20} | {rho2:>+14.3f} {p2:>8.4f} | {cv2:>11.1f}% | {'🔥 CONSERVED' if abs(rho2) < 0.3 else 'NOT CONSERVED'}")

# ============================================================
# ALSO TEST ON CIII (the other conserved one)
# ============================================================
print(f"\n\n{'=' * 130}")
print("CIII_BR: SANITY CHECK")
print("=" * 130)

ciii = d['CIII_BR']
mask_c = ((z_all >= 0.8) & (z_all < 3.0) &
          (ciii[:, 4] > 100) & (ciii[:, 4] < 30000) &
          (ciii[:, 2] > 0) & (ciii[:, 2] < 5000) &
          np.isfinite(ciii[:, 2]) & np.isfinite(ciii[:, 4]) &
          (lbol > 40) & (lbol < 50))

z_c = z_all[mask_c]
log_ew_c = np.log10(ciii[mask_c, 2])
lbol_c = lbol[mask_c]
coef_c = P.polyfit(lbol_c, log_ew_c, 1)
log_ew_cr = log_ew_c - P.polyval(lbol_c, coef_c)

global_min_c, global_max_c = np.percentile(log_ew_cr, 0.5), np.percentile(log_ew_cr, 99.5)
fixed_bins_c = np.linspace(global_min_c, global_max_c, 51)
global_std_c = np.std(log_ew_cr)
fixed_bw_c = global_std_c * (4 / (3 * len(log_ew_cr))) ** 0.2

z_edges_c = np.percentile(z_c, np.linspace(0, 100, n_zbins + 1))
z_edges_c = np.unique(np.round(z_edges_c, 3))

z_list_c, var_list_c = [], []
et_fixedbin_c, et_fixedbw_c, et_gaussian_c = [], [], []

for i in range(len(z_edges_c) - 1):
    z_lo, z_hi = z_edges_c[i], z_edges_c[i+1]
    bm = (z_c >= z_lo) & (z_c < z_hi)
    if np.sum(bm) < 500:
        continue
    zc = np.mean(z_c[bm])
    x = log_ew_cr[bm]
    var = np.var(x)
    peaks = compute_peak_methods(x, fixed_bins_c, fixed_bw_c)
    
    z_list_c.append(zc)
    var_list_c.append(var)
    et_fixedbin_c.append(var * peaks['peak_fixedbin'])
    et_fixedbw_c.append(var * peaks['peak_fixedbw'])
    et_gaussian_c.append(var * peaks['peak_gaussian'])

print(f"\n  {'Method':>20} | {'ρ(étendue,z)':>14} {'p':>8} | Verdict")
print("  " + "-" * 55)
for label, vals in [('Fixed-bin hist', et_fixedbin_c), ('Fixed-bw KDE', et_fixedbw_c), ('Gaussian peak', et_gaussian_c)]:
    rho, p = spearmanr(z_list_c, vals)
    print(f"  {label:>20} | {rho:>+14.3f} {p:>8.4f} | {'🔥 CONSERVED' if abs(rho) < 0.3 else 'ARTIFACT' if abs(rho) > 0.5 else 'WEAK'}")

# ============================================================
# THE GAUSSIAN IDENTITY CHECK
# ============================================================
print(f"\n\n{'=' * 130}")
print("GAUSSIAN IDENTITY CHECK: Var × (1/√(2πVar)) = √(Var/(2π)) — this ALWAYS grows with Var")
print("For Gaussian, Var×Peak = σ/√(2π). NOT constant unless σ is constant.")
print("So if Var×Peak IS constant, the distribution is NOT Gaussian — it's doing something special.")
print("=" * 130)

print(f"\n  MgII Var×GaussianPeak values: {['%.4f' % v for v in [v * (1/(np.sqrt(2*np.pi)*np.sqrt(v))) for v in var_list]]}")
print(f"  = √(Var/2π) = {['%.4f' % np.sqrt(v/(2*np.pi)) for v in var_list]}")
rho_g, p_g = spearmanr(z_list, [np.sqrt(v/(2*np.pi)) for v in var_list])
print(f"  ρ(√(Var/2π), z) = {rho_g:+.3f}, p = {p_g:.4f}")
print(f"  If this is NOT conserved but adaptive KDE étendue IS → the non-Gaussianity is doing the work")
