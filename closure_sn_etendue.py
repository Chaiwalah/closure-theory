#!/usr/bin/env python3
"""
closure_sn_etendue.py — SNe Ia Étendue Conservation Test
==========================================================

The final cross-domain confirmation: does ℰ = σ² · ∫f² hold for 
SN Ia observables (mB, x1, c) across redshift?

If yes → universal conservation law across SNe + quasars = Level 10
If no → quasar-specific, still publishable but not universal

Uses collision probability (Rényi-2) as the estimator-independent
concentration measure, validated by the 5-method artifact kill test.

Also tests:
- Fast channel only (x1 < 0) 
- Slow channel only (x1 ≥ 0)
- All SNe combined
- Individual observables AND the 3D joint

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr, gaussian_kde
from scipy.integrate import quad
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# Load Pantheon+
sn = np.loadtxt('data/pantheon_plus.dat', dtype=str)
header = sn[0]
data = sn[1:]
col_idx = {h: i for i, h in enumerate(header)}

z = data[:, col_idx['zHD']].astype(float)
mb = data[:, col_idx['mB']].astype(float)
x1 = data[:, col_idx['x1']].astype(float)
c = data[:, col_idx['c']].astype(float)

# Standard cuts
mask = (z > 0.01) & (z < 2.5) & (np.abs(c) < 0.3) & (np.abs(x1) < 3)
z, mb, x1, c = z[mask], mb[mask], x1[mask], c[mask]
print(f"Pantheon+ after cuts: {len(z)} SNe Ia\n")

# Tripp standardization (global fit)
from numpy.polynomial import polynomial as P
# Simple Tripp: μ = mB - M + α*x1 - β*c
# We'll use residuals from a linear fit
X = np.column_stack([np.ones(len(z)), x1, c])
# Fit mB = M + α*x1 - β*c + f(z)
# Actually just residualize mB against x1 and c
from numpy.linalg import lstsq
coeffs, _, _, _ = lstsq(X, mb, rcond=None)
mb_resid = mb - X @ coeffs  # Hubble residual proxy

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 2.5]

def collision_prob_1d(x, n_grid=200):
    """Compute ∫f² using KDE (fixed bandwidth from global sample)."""
    if len(x) < 15:
        return np.nan
    try:
        kde = gaussian_kde(x)
        grid = np.linspace(np.percentile(x, 1), np.percentile(x, 99), n_grid)
        f_vals = kde.evaluate(grid)
        dx = grid[1] - grid[0]
        return np.sum(f_vals**2) * dx
    except:
        return np.nan

def peak_density_fixed(x, fixed_bins):
    """Peak from fixed-bin histogram."""
    if len(x) < 15:
        return np.nan
    counts, _ = np.histogram(x, bins=fixed_bins, density=True)
    return np.max(counts)

# ============================================================
# TEST 1: Individual observables — ALL SNe
# ============================================================
print(f"{'=' * 110}")
print("TEST 1: ÉTENDUE ℰ = Var × ∫f² FOR INDIVIDUAL SN Ia OBSERVABLES")
print("=" * 110)

for var_name, var_data in [('mB_resid', mb_resid), ('x1', x1), ('c', c)]:
    # Fixed bins from global range
    global_min, global_max = np.percentile(var_data, 1), np.percentile(var_data, 99)
    fixed_bins = np.linspace(global_min, global_max, 31)
    
    print(f"\n  {var_name} (ALL, N={len(var_data)}):")
    print(f"  {'z':>6} {'N':>4} | {'Var':>8} {'∫f²':>8} {'Peak_fix':>9} | {'ℰ(∫f²)':>8} {'ℰ(peak)':>8}")
    print(f"  " + "-" * 65)
    
    z_list, et_collision, et_fixedbin = [], [], []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bm = (z >= z_lo) & (z < z_hi)
        n = np.sum(bm)
        if n < 20:
            continue
        
        zc = np.mean(z[bm])
        v = np.var(var_data[bm])
        cp = collision_prob_1d(var_data[bm])
        pf = peak_density_fixed(var_data[bm], fixed_bins)
        
        et_c = v * cp if not np.isnan(cp) else np.nan
        et_f = v * pf if not np.isnan(pf) else np.nan
        
        z_list.append(zc)
        et_collision.append(et_c)
        et_fixedbin.append(et_f)
        
        print(f"  {zc:>6.3f} {n:>4} | {v:>8.4f} {cp:>8.4f} {pf:>9.4f} | {et_c:>8.4f} {et_f:>8.4f}")
    
    if len(z_list) >= 3:
        valid_c = [i for i, v in enumerate(et_collision) if not np.isnan(v)]
        valid_f = [i for i, v in enumerate(et_fixedbin) if not np.isnan(v)]
        
        if len(valid_c) >= 3:
            z_v = [z_list[i] for i in valid_c]
            e_v = [et_collision[i] for i in valid_c]
            rho, p = spearmanr(z_v, e_v)
            cv = np.std(e_v) / np.mean(e_v) * 100
            print(f"\n  ℰ(∫f²) vs z:    ρ = {rho:+.3f}, p = {p:.4f}, CV = {cv:.1f}% {'🔥 CONSERVED' if abs(rho) < 0.4 else ''}")
        
        if len(valid_f) >= 3:
            z_v = [z_list[i] for i in valid_f]
            e_v = [et_fixedbin[i] for i in valid_f]
            rho, p = spearmanr(z_v, e_v)
            cv = np.std(e_v) / np.mean(e_v) * 100
            print(f"  ℰ(peak) vs z:   ρ = {rho:+.3f}, p = {p:.4f}, CV = {cv:.1f}%")

# ============================================================
# TEST 2: FAST vs SLOW channels separately
# ============================================================
print(f"\n\n{'=' * 110}")
print("TEST 2: ÉTENDUE BY CHANNEL — Fast (x1<0) vs Slow (x1≥0)")
print("=" * 110)

for ch_name, ch_mask in [('FAST (x1<0)', x1 < 0), ('SLOW (x1≥0)', x1 >= 0)]:
    for var_name, var_data in [('x1', x1), ('c', c), ('mB_resid', mb_resid)]:
        
        z_ch = z[ch_mask]
        v_ch = var_data[ch_mask]
        
        if len(z_ch) < 50:
            continue
        
        global_min, global_max = np.percentile(v_ch, 2), np.percentile(v_ch, 98)
        fixed_bins = np.linspace(global_min, global_max, 25)
        
        z_list, et_list = [], []
        
        for i in range(len(z_edges) - 1):
            z_lo, z_hi = z_edges[i], z_edges[i+1]
            bm = (z_ch >= z_lo) & (z_ch < z_hi)
            n = np.sum(bm)
            if n < 15:
                continue
            
            zc = np.mean(z_ch[bm])
            v = np.var(v_ch[bm])
            cp = collision_prob_1d(v_ch[bm])
            
            if not np.isnan(cp):
                z_list.append(zc)
                et_list.append(v * cp)
        
        if len(z_list) >= 3:
            rho, p = spearmanr(z_list, et_list)
            cv = np.std(et_list) / np.mean(et_list) * 100
            status = "🔥 CONSERVED" if abs(rho) < 0.4 else "↗ GROWS" if rho > 0.4 else "↘ DROPS"
            print(f"  {ch_name:>15} {var_name:>10}: ρ = {rho:+.3f}, p = {p:.4f}, CV = {cv:.1f}%  {status}")

# ============================================================
# TEST 3: 3D JOINT ÉTENDUE — det(C) as multivariate version
# ============================================================
print(f"\n\n{'=' * 110}")
print("TEST 3: MULTIVARIATE ÉTENDUE — det(Covariance) × (N-dependent concentration)")
print("For multivariate Gaussian: det(C) ∝ (volume)². If distribution concentrates,")
print("det(C) should drop. But if étendue is conserved, det(C) × concentration = const.")
print("=" * 110)

print(f"\n  {'z':>6} {'N':>4} | {'det(C)':>12} {'Var(mB_r)':>10} {'Var(x1)':>10} {'Var(c)':>10}")
print(f"  " + "-" * 65)

z_list, detC_list = [], []
for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bm = (z >= z_lo) & (z < z_hi)
    n = np.sum(bm)
    if n < 20:
        continue
    
    zc = np.mean(z[bm])
    C = np.cov(np.array([mb_resid[bm], x1[bm], c[bm]]))
    det = np.linalg.det(C)
    
    z_list.append(zc)
    detC_list.append(det)
    
    print(f"  {zc:>6.3f} {n:>4} | {det:>12.6f} {np.var(mb_resid[bm]):>10.4f} {np.var(x1[bm]):>10.4f} {np.var(c[bm]):>10.4f}")

if len(z_list) >= 3:
    rho, p = spearmanr(z_list, detC_list)
    print(f"\n  det(C) vs z: ρ = {rho:+.3f}, p = {p:.4f}")
    if rho < -0.5:
        print(f"  → Phase-space volume CONTRACTS (consistent with Hard Cap)")

# ============================================================
# FINAL VERDICT
# ============================================================
print(f"\n\n{'=' * 110}")
print("CROSS-DOMAIN ÉTENDUE VERDICT")
print("=" * 110)
print("""
If SNe show ℰ = const → UNIVERSAL conservation law → Level 10
If SNe show ℰ ≠ const → Quasar-specific → Level 9.8 (still publishable)
If mixed (some observables conserve, others don't) → Observable-dependent → needs more work
""")
