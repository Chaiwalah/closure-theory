#!/usr/bin/env python3
"""
TEMPORAL EATING TEST — Does the operator eat light curve shape?
================================================================
If the operator eats relational information during propagation,
temporal coherence (light curve shape) should degrade with distance
beyond what stretch/color correction accounts for.

Tests:
A) Light curve fit quality (chi2/dof) vs redshift — does fitting get worse?
B) Stretch (x1) residual scatter vs redshift — does shape diversity inflate?
C) Color (c) residual scatter vs redshift — does color scatter inflate?
D) Hubble residual structure vs redshift — development constraint signature
E) x1-c correlation vs redshift — does temporal-chromatic coupling degrade?
"""

import numpy as np
import pandas as pd
from scipy import stats
import os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_temporal'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("TEMPORAL EATING TEST — Light Curve Shape Degradation")
print("=" * 70)

df = pd.read_csv('/root/clawd/projects/closure-theory/results/pantheon_with_y.csv')
print(f"  Pantheon+ sample: {len(df)} SNe Ia")

# Clean
df = df[df['zHD'] > 0.01].copy()  # remove very local
df = df[np.isfinite(df['x1']) & np.isfinite(df['c']) & np.isfinite(df['zHD'])].copy()
df = df[np.isfinite(df['FITCHI2']) & np.isfinite(df['NDOF']) & (df['NDOF'] > 0)].copy()
df['chi2_dof_val'] = df['FITCHI2'] / df['NDOF']

print(f"  After cleaning: {len(df)} SNe Ia")
print(f"  z range: {df['zHD'].min():.4f} - {df['zHD'].max():.4f}")

# ============================================================
# TEST A: FIT QUALITY vs REDSHIFT
# ============================================================
print("\n" + "=" * 70)
print("TEST A: Light curve fit quality (χ²/dof) vs redshift")
print("=" * 70)
print("  If temporal coherence degrades, fits should get WORSE at high-z")

z = df['zHD'].values
chi2dof = df['chi2_dof_val'].values

# Bin by redshift
n_bins = 8
edges = np.percentile(z, np.linspace(0, 100, n_bins + 1))
results_A = []
for i in range(n_bins):
    mask = (z >= edges[i]) & (z < edges[i+1])
    if i == n_bins - 1: mask = (z >= edges[i]) & (z <= edges[i+1])
    zbin = z[mask]
    chi2bin = chi2dof[mask]
    results_A.append({
        'z_mid': np.median(zbin),
        'chi2_median': np.median(chi2bin),
        'chi2_mean': np.mean(chi2bin),
        'chi2_std': np.std(chi2bin),
        'N': np.sum(mask)
    })
    print(f"  z={np.median(zbin):.3f}: χ²/dof median={np.median(chi2bin):.3f}, "
          f"mean={np.mean(chi2bin):.3f}, std={np.std(chi2bin):.3f}, N={np.sum(mask)}")

sl, _, _, p, _ = stats.linregress(z, chi2dof)
print(f"\n  Overall trend: slope = {sl:+.4f}, p = {p:.2e}")

# Spearman (robust)
rho, p_s = stats.spearmanr(z, chi2dof)
print(f"  Spearman: ρ = {rho:+.4f}, p = {p_s:.2e}")

# ============================================================
# TEST B: STRETCH (x1) SCATTER vs REDSHIFT
# ============================================================
print("\n" + "=" * 70)
print("TEST B: Stretch (x1) scatter vs redshift")
print("=" * 70)
print("  If temporal shape degrades, x1 dispersion should increase with z")

x1 = df['x1'].values
x1_err = df['x1ERR'].values

results_B = []
for i in range(n_bins):
    mask = (z >= edges[i]) & (z < edges[i+1])
    if i == n_bins - 1: mask = (z >= edges[i]) & (z <= edges[i+1])
    x1bin = x1[mask]
    x1_err_bin = x1_err[mask]
    # Intrinsic scatter: observed - measurement noise
    obs_var = np.var(x1bin)
    noise_var = np.median(x1_err_bin**2)
    intr_scatter = np.sqrt(max(obs_var - noise_var, 0))
    results_B.append({
        'z_mid': np.median(z[mask]),
        'x1_std': np.std(x1bin),
        'x1_intrinsic': intr_scatter,
        'x1_err_med': np.median(x1_err_bin),
        'N': np.sum(mask)
    })
    print(f"  z={np.median(z[mask]):.3f}: x1 std={np.std(x1bin):.3f}, "
          f"intrinsic={intr_scatter:.3f}, err={np.median(x1_err_bin):.3f}, N={np.sum(mask)}")

# ============================================================
# TEST C: COLOR (c) SCATTER vs REDSHIFT
# ============================================================
print("\n" + "=" * 70)
print("TEST C: Color (c) scatter vs redshift")
print("=" * 70)

c = df['c'].values
c_err = df['cERR'].values

results_C = []
for i in range(n_bins):
    mask = (z >= edges[i]) & (z < edges[i+1])
    if i == n_bins - 1: mask = (z >= edges[i]) & (z <= edges[i+1])
    cbin = c[mask]
    c_err_bin = c_err[mask]
    obs_var = np.var(cbin)
    noise_var = np.median(c_err_bin**2)
    intr_scatter = np.sqrt(max(obs_var - noise_var, 0))
    results_C.append({
        'z_mid': np.median(z[mask]),
        'c_std': np.std(cbin),
        'c_intrinsic': intr_scatter,
        'c_err_med': np.median(c_err_bin),
        'N': np.sum(mask)
    })
    print(f"  z={np.median(z[mask]):.3f}: c std={np.std(cbin):.3f}, "
          f"intrinsic={intr_scatter:.3f}, err={np.median(c_err_bin):.3f}, N={np.sum(mask)}")

# ============================================================
# TEST D: HUBBLE RESIDUAL SCATTER vs REDSHIFT
# ============================================================
print("\n" + "=" * 70)
print("TEST D: Hubble residual scatter vs redshift")
print("=" * 70)
print("  After all corrections, is there EXTRA scatter at high-z?")

if 'mu_resid' in df.columns:
    mu_res = df['mu_resid'].values
    
    results_D = []
    for i in range(n_bins):
        mask = (z >= edges[i]) & (z < edges[i+1])
        if i == n_bins - 1: mask = (z >= edges[i]) & (z <= edges[i+1])
        resid = mu_res[mask]
        results_D.append({
            'z_mid': np.median(z[mask]),
            'resid_std': np.std(resid),
            'resid_mad': 1.4826 * np.median(np.abs(resid - np.median(resid))),
            'N': np.sum(mask)
        })
        print(f"  z={np.median(z[mask]):.3f}: σ(μ_resid)={np.std(resid):.4f}, "
              f"MAD={1.4826*np.median(np.abs(resid-np.median(resid))):.4f}, N={np.sum(mask)}")

# ============================================================
# TEST E: x1-c CORRELATION vs REDSHIFT (THE CLOSURE TEST)
# ============================================================
print("\n" + "=" * 70)
print("TEST E: Stretch-Color coupling vs redshift (temporal-chromatic)")
print("=" * 70)
print("  x1 (temporal shape) × c (chromatic) = relational information")
print("  If operator eats correlations, this coupling should degrade with z")

results_E = []
for i in range(n_bins):
    mask = (z >= edges[i]) & (z < edges[i+1])
    if i == n_bins - 1: mask = (z >= edges[i]) & (z <= edges[i+1])
    if np.sum(mask) < 20: continue
    r, p_val = stats.spearmanr(x1[mask], c[mask])
    results_E.append({
        'z_mid': np.median(z[mask]),
        'rho': r,
        'p': p_val,
        'N': np.sum(mask)
    })
    print(f"  z={np.median(z[mask]):.3f}: ρ(x1,c) = {r:+.4f}, p = {p_val:.2e}, N={np.sum(mask)}")

if len(results_E) >= 3:
    z_mids = [r['z_mid'] for r in results_E]
    rhos = [r['rho'] for r in results_E]
    trend_rho, _, _, trend_p, _ = stats.linregress(z_mids, rhos)
    print(f"\n  Coupling trend: slope = {trend_rho:+.4f}, p = {trend_p:.2e}")
    print(f"  Low-z ρ = {rhos[0]:+.4f}, High-z ρ = {rhos[-1]:+.4f}")
    delta_pct = (rhos[-1] - rhos[0]) / abs(rhos[0]) * 100 if rhos[0] != 0 else 0
    print(f"  Change: {delta_pct:+.1f}%")

# ============================================================
# TEST F: DEVELOPMENT CONSTRAINT — Selection at high-z
# ============================================================
print("\n" + "=" * 70)
print("TEST F: Development constraint — population narrowing")
print("=" * 70)
print("  If operator constrains what CAN be observed, high-z population")
print("  should be NARROWER (only 'simple' SNe survive selection)")

for param, label in [(x1, 'x1 (stretch)'), (c, 'c (color)')]:
    print(f"\n  {label}:")
    lo_z = param[(z > 0.01) & (z < 0.1)]
    hi_z = param[z > 0.5]
    if len(lo_z) > 10 and len(hi_z) > 10:
        print(f"    Low-z (z<0.1):  mean={np.mean(lo_z):+.3f}, std={np.std(lo_z):.3f}, N={len(lo_z)}")
        print(f"    High-z (z>0.5): mean={np.mean(hi_z):+.3f}, std={np.std(hi_z):.3f}, N={len(hi_z)}")
        # Levene test for equal variances
        stat_l, p_l = stats.levene(lo_z, hi_z)
        print(f"    Levene test (equal variance): F={stat_l:.2f}, p={p_l:.2e}")
        if np.std(hi_z) < np.std(lo_z):
            print(f"    → High-z is NARROWER ({100*(1-np.std(hi_z)/np.std(lo_z)):.1f}% less dispersed)")
        else:
            print(f"    → High-z is WIDER ({100*(np.std(hi_z)/np.std(lo_z)-1):.1f}% more dispersed)")

# ============================================================
# TEST G: x1-c COVARIANCE MATRIX vs Z
# ============================================================
print("\n" + "=" * 70)
print("TEST G: Joint x1-c covariance structure vs redshift")
print("=" * 70)

if 'COV_x1_c' in df.columns:
    cov_x1c = df['COV_x1_c'].values
    for i in range(n_bins):
        mask = (z >= edges[i]) & (z < edges[i+1])
        if i == n_bins - 1: mask = (z >= edges[i]) & (z <= edges[i+1])
        print(f"  z={np.median(z[mask]):.3f}: median COV(x1,c) = {np.median(cov_x1c[mask]):.6f}, N={np.sum(mask)}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"  A) Fit quality trend:      slope = {sl:+.4f}, ρ = {rho:+.4f}")
if results_E:
    print(f"  E) x1-c coupling:         low-z ρ = {results_E[0]['rho']:+.3f} → high-z ρ = {results_E[-1]['rho']:+.3f}")
print(f"\n  Key question: Does temporal-chromatic coupling degrade like spectral coupling?")
print(f"  If yes → operator eats time-encoded information too")

print("\n✅ Done. Results in", OUTDIR)
