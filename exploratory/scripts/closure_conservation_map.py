#!/usr/bin/env python3
"""
closure_conservation_map.py — Untangling the Conservation Threads
==================================================================

Three threads Humza identified:
1. The "amplified" locked signal — is it the TRUE energy footprint?
2. Conservation: locked + diagnostic combined should be constant per epoch
3. Any discrepancy = what we call "dark energy"

Tests:
A. Does α(z) (stretch coefficient, locked) RISE as β(z) (color, diagnostic) FALLS?
B. Is total information content (locked + diagnostic) conserved with z?
C. Does the conservation residual match the dark energy parameterization?

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import spearmanr, pearsonr
import json
from pathlib import Path

np.random.seed(42)

# ============================================================
# LOAD PANTHEON+ DATA
# ============================================================

data = []
with open('data/pantheon_plus.dat', 'r') as f:
    header = f.readline().strip().split()
    for line in f:
        parts = line.strip().split()
        if len(parts) >= len(header):
            row = {}
            for i, h in enumerate(header):
                try:
                    row[h] = float(parts[i])
                except:
                    row[h] = parts[i]
            data.append(row)

# Filter and deduplicate
seen = set()
sne = []
for row in data:
    z = row.get('zHD', 0)
    c = row.get('c', -999)
    x1 = row.get('x1', -999)
    mB = row.get('mB', -999)
    cErr = row.get('cERR', -999)
    x1Err = row.get('x1ERR', -999)
    
    if z > 0.01 and z < 2.5 and abs(c) < 0.3 and abs(x1) < 3 and cErr > 0 and x1Err > 0:
        key = (row['CID'], round(z, 5))
        if key not in seen:
            seen.add(key)
            sne.append(row)

print(f"Loaded {len(sne)} unique SNe Ia")

# ============================================================
# TEST A: α(z) vs β(z) — Does locked rise as diagnostic falls?
# ============================================================
print("\n" + "=" * 70)
print("TEST A: α(z) vs β(z) — Locked/Diagnostic Channel Transfer")
print("=" * 70)

# Compute α and β in redshift bins
z_edges = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.55, 0.7, 0.85, 1.0, 1.5, 2.5]

alphas = []
betas = []
alpha_errs = []
beta_errs = []
z_centers = []
n_per_bin = []

print(f"\n{'z_range':<14} {'N':>5} {'α':>8} {'α_err':>8} {'β':>8} {'β_err':>8} {'α+β':>8}")
print("-" * 70)

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_data = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    
    if len(bin_data) < 10:
        continue
    
    colors = np.array([s['c'] for s in bin_data])
    stretches = np.array([s['x1'] for s in bin_data])
    mBs = np.array([s['mB'] for s in bin_data])
    
    # Fit: mB = M + β*c - α*x1 (note: α has opposite sign in Tripp)
    A = np.column_stack([np.ones(len(bin_data)), colors, stretches])
    try:
        params, res, rank, sv = np.linalg.lstsq(A, mBs, rcond=None)
        M_fit, beta_fit, neg_alpha_fit = params
        alpha_fit = -neg_alpha_fit  # Tripp convention: brighter for larger x1
        
        # Bootstrap errors
        n_boot = 300
        a_boots, b_boots = [], []
        for _ in range(n_boot):
            idx = np.random.choice(len(bin_data), len(bin_data), replace=True)
            try:
                p = np.linalg.lstsq(A[idx], mBs[idx], rcond=None)[0]
                b_boots.append(p[1])
                a_boots.append(-p[2])
            except:
                pass
        
        a_err = np.std(a_boots) if a_boots else 999
        b_err = np.std(b_boots) if b_boots else 999
        
        z_center = np.mean([s['zHD'] for s in bin_data])
        
        alphas.append(alpha_fit)
        betas.append(beta_fit)
        alpha_errs.append(a_err)
        beta_errs.append(b_err)
        z_centers.append(z_center)
        n_per_bin.append(len(bin_data))
        
        print(f"[{z_lo:.2f},{z_hi:.2f}){'':<5} {len(bin_data):>5} {alpha_fit:>8.3f} {a_err:>8.3f} {beta_fit:>8.3f} {b_err:>8.3f} {alpha_fit+beta_fit:>8.3f}")
    except:
        pass

z_arr = np.array(z_centers)
alpha_arr = np.array(alphas)
beta_arr = np.array(betas)

# Test: does α rise as β falls?
if len(z_arr) >= 4:
    rho_alpha_z, p_alpha_z = spearmanr(z_arr, alpha_arr)
    rho_beta_z, p_beta_z = spearmanr(z_arr, beta_arr)
    rho_alpha_beta, p_alpha_beta = spearmanr(alpha_arr, beta_arr)
    
    print(f"\n  α vs z: ρ = {rho_alpha_z:+.3f}, p = {p_alpha_z:.4f}")
    print(f"  β vs z: ρ = {rho_beta_z:+.3f}, p = {p_beta_z:.4f}")
    print(f"  α vs β: ρ = {rho_alpha_beta:+.3f}, p = {p_alpha_beta:.4f}")
    
    if rho_alpha_z > 0 and rho_beta_z < 0:
        print(f"\n  🔥 α RISES while β FALLS → locked channel AMPLIFIED as diagnostic compressed!")
        print(f"     This is the conservation transfer: diagnostic → locked")
    elif rho_beta_z < 0 and abs(rho_alpha_z) < 0.3:
        print(f"\n  β falls but α is flat → diagnostic degrades, locked stays constant")
        print(f"     No amplification of locked channel")
    else:
        print(f"\n  Unexpected pattern — needs investigation")

# ============================================================
# TEST B: Total "Information Content" Conservation
# ============================================================
print("\n" + "=" * 70)
print("TEST B: Total Information Content — Is α + β Conserved?")
print("=" * 70)

# α + β is a proxy for total standardization power
# If conservation: α + β = constant across z
total = alpha_arr + beta_arr

rho_total_z, p_total_z = spearmanr(z_arr, total)
print(f"\n  α + β vs z: ρ = {rho_total_z:+.3f}, p = {p_total_z:.4f}")
print(f"  Mean α + β = {np.mean(total):.3f} ± {np.std(total):.3f}")
print(f"  CV = {np.std(total)/np.mean(total)*100:.1f}%")

if abs(rho_total_z) < 0.3 and p_total_z > 0.1:
    print(f"\n  ★ α + β is CONSISTENT with constant → total info CONSERVED")
elif rho_total_z < -0.3:
    print(f"\n  ⚠️ α + β DROPS with z → total info NOT conserved")
    print(f"     The residual may be what we call 'dark energy'")
elif rho_total_z > 0.3:
    print(f"\n  ⚠️ α + β RISES with z → unexpected")

# ============================================================
# TEST C: Weighted Information Content (q-weighted)
# ============================================================
print("\n" + "=" * 70)
print("TEST C: q-Weighted Information — α/q_stretch + β/q_color")  
print("=" * 70)

# q values from our derivation
q_color = 1.0
q_stretch = 0.039

# If conservation holds in q-space:
# The total information carried, weighted by susceptibility, should be conserved
# Stretch carries info at q=0.039 (nearly locked)
# Color carries info at q=1.0 (fully diagnostic)

# Normalized contributions
I_locked = alpha_arr * q_stretch  # locked channel information proxy
I_diag = beta_arr * q_color       # diagnostic channel information proxy
I_total = I_locked + I_diag

rho_locked_z, p_locked = spearmanr(z_arr, I_locked)
rho_diag_z, p_diag = spearmanr(z_arr, I_diag)
rho_total_z2, p_total2 = spearmanr(z_arr, I_total)

print(f"\n  I_locked (α×q_stretch) vs z: ρ = {rho_locked_z:+.3f}, p = {p_locked:.4f}")
print(f"  I_diag (β×q_color) vs z:     ρ = {rho_diag_z:+.3f}, p = {p_diag:.4f}")
print(f"  I_total vs z:                ρ = {rho_total_z2:+.3f}, p = {p_total2:.4f}")

print(f"\n  {'z':>6} {'I_lock':>8} {'I_diag':>8} {'I_total':>8} {'ratio':>8}")
print(f"  " + "-" * 45)
for i in range(len(z_arr)):
    ratio = I_diag[i] / I_locked[i] if I_locked[i] != 0 else 999
    print(f"  {z_arr[i]:>6.2f} {I_locked[i]:>8.4f} {I_diag[i]:>8.3f} {I_total[i]:>8.3f} {ratio:>8.1f}")

# ============================================================
# TEST D: Scatter as Information — Variance Transfer
# ============================================================
print("\n" + "=" * 70)
print("TEST D: Variance Transfer — Color scatter vs Stretch scatter")
print("=" * 70)

# If diagnostic info is being lost, color SCATTER should increase with z
# (noisier measurements as signal degrades)
# If locked is amplified, stretch scatter should... stay constant or decrease

print(f"\n{'z_range':<14} {'N':>5} {'σ(c)':>8} {'σ(x1)':>8} {'σ(c)/σ(x1)':>12}")
print("-" * 55)

color_scatters = []
stretch_scatters = []
z_scat_centers = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_data = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    
    if len(bin_data) < 15:
        continue
    
    colors = np.array([s['c'] for s in bin_data])
    stretches = np.array([s['x1'] for s in bin_data])
    
    # Intrinsic scatter (subtract measurement error in quadrature)
    c_errs = np.array([s['cERR'] for s in bin_data])
    x1_errs = np.array([s['x1ERR'] for s in bin_data])
    
    raw_c_var = np.var(colors)
    raw_x1_var = np.var(stretches)
    mean_c_err_var = np.mean(c_errs**2)
    mean_x1_err_var = np.mean(x1_errs**2)
    
    intrinsic_c = np.sqrt(max(0, raw_c_var - mean_c_err_var))
    intrinsic_x1 = np.sqrt(max(0, raw_x1_var - mean_x1_err_var))
    
    z_center = np.mean([s['zHD'] for s in bin_data])
    ratio = intrinsic_c / intrinsic_x1 if intrinsic_x1 > 0 else 999
    
    color_scatters.append(intrinsic_c)
    stretch_scatters.append(intrinsic_x1)
    z_scat_centers.append(z_center)
    
    print(f"[{z_lo:.2f},{z_hi:.2f}){'':<5} {len(bin_data):>5} {intrinsic_c:>8.4f} {intrinsic_x1:>8.3f} {ratio:>12.3f}")

if len(z_scat_centers) >= 4:
    rho_c, p_c = spearmanr(z_scat_centers, color_scatters)
    rho_x1, p_x1 = spearmanr(z_scat_centers, stretch_scatters)
    
    print(f"\n  σ(color) vs z: ρ = {rho_c:+.3f}, p = {p_c:.4f}")
    print(f"  σ(stretch) vs z: ρ = {rho_x1:+.3f}, p = {p_x1:.4f}")
    
    if rho_c > 0:
        print(f"  → Color scatter INCREASES with z (diagnostic channel getting noisier)")
    if rho_x1 < 0:
        print(f"  → Stretch scatter DECREASES with z (locked channel getting tighter)")
    elif abs(rho_x1) < 0.2:
        print(f"  → Stretch scatter FLAT with z (locked channel immune)")

# ============================================================
# TEST E: Energy Footprint — Does mB residual match dark energy?
# ============================================================
print("\n" + "=" * 70)
print("TEST E: Conservation Residual vs Dark Energy")
print("=" * 70)

# If total info is NOT perfectly conserved, the residual should match
# the dark energy signal (w ≈ -0.77 from our earlier Tripp analysis)

# Use α + β trend to compute "missing information" at each z
if len(z_arr) >= 4:
    # Fit linear trend to α + β
    from numpy.polynomial import polynomial as P
    coeffs = np.polyfit(z_arr, total, 1)
    total_fit = np.polyval(coeffs, z_arr)
    residual = total - total_fit
    
    # Compare trend to our predicted dark energy compression
    # From closure_tripp_bias.py: Δw ≈ -0.23 from compression
    # This should show up as a specific slope in (α + β) vs z
    
    slope = coeffs[0]
    intercept = coeffs[1]
    
    print(f"  Linear fit: (α + β) = {intercept:.3f} + {slope:+.3f} × z")
    print(f"  At z=0: α + β = {intercept:.3f}")
    print(f"  At z=1: α + β = {intercept + slope:.3f}")
    print(f"  Change: {slope:.3f} ({slope/intercept*100:+.1f}%)")
    
    if slope < 0:
        print(f"\n  Total standardization power DECREASES with z")
        print(f"  Lost power per unit z: {abs(slope):.3f}")
        print(f"  This 'missing' standardization ≈ information cost of looking backward")
    elif abs(slope) < 0.5:
        print(f"\n  Total standardization power roughly CONSTANT")
        print(f"  Conservation approximately holds")

# ============================================================
# SAVE
# ============================================================

results_dir = Path("results_conservation_map")
results_dir.mkdir(exist_ok=True)

results = {
    'z_centers': [float(z) for z in z_arr],
    'alphas': [float(a) for a in alpha_arr],
    'betas': [float(b) for b in beta_arr],
    'alpha_plus_beta': [float(t) for t in total],
    'test_A': {
        'alpha_vs_z_rho': float(rho_alpha_z),
        'beta_vs_z_rho': float(rho_beta_z),
        'alpha_vs_beta_rho': float(rho_alpha_beta),
    },
    'test_B': {
        'total_vs_z_rho': float(rho_total_z),
        'total_vs_z_p': float(p_total_z),
    }
}

with open(results_dir / "conservation_map.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {results_dir}/")
print("\n" + "=" * 70)
print("END — Conservation Map")
print("=" * 70)
