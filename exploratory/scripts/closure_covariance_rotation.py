#!/usr/bin/env python3
"""
closure_covariance_rotation.py — GPT's Four-Point Diagnostic
==============================================================

GPT asks:
1. Does conditional variance of distance stay flat while color leverage drops?
2. Do the principal axes of the (c, x1, mB) covariance rotate with z?
3. Does the surviving low-noise axis align with stretch (geometric) vs color (diagnostic)?
4. Does foreground environment change the rotation rate?

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
                try:
                    row[h] = float(parts[i])
                except:
                    row[h] = parts[i]
            data.append(row)

seen = set()
sne = []
for row in data:
    z = row.get('zHD', 0)
    c = row.get('c', -999)
    x1 = row.get('x1', -999)
    if z > 0.01 and z < 2.5 and abs(c) < 0.3 and abs(x1) < 3:
        key = (row['CID'], round(z, 5))
        if key not in seen:
            seen.add(key)
            sne.append(row)

# Galactic latitude
def gal_b(ra_deg, dec_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    ra_ngp, dec_ngp = np.radians(192.8595), np.radians(27.1284)
    sin_b = np.sin(dec)*np.sin(dec_ngp) + np.cos(dec)*np.cos(dec_ngp)*np.cos(ra-ra_ngp)
    return np.degrees(np.arcsin(np.clip(sin_b, -1, 1)))

for s in sne:
    s['abs_gal_b'] = abs(gal_b(s['RA'], s['DEC']))

print(f"Loaded {len(sne)} unique SNe Ia\n")

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.82, 1.0, 2.5]

# ============================================================
# Q1: Conditional variance flat while color leverage drops?
# ============================================================
print("=" * 80)
print("Q1: Conditional Variance vs Color Leverage")
print("=" * 80)

print(f"\n{'z_range':<12} {'N':>4} {'R²_c':>7} {'R²_x1':>7} {'Var(m|c,x1)':>12} {'leverage_c':>10} {'leverage_x1':>11}")
print("-" * 75)

z_centers = []
R2c_list = []
R2x_list = []
condvar_list = []
leverage_c_list = []
leverage_x1_list = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 15:
        continue
    
    mB = np.array([s['mB'] for s in bin_sne])
    x1 = np.array([s['x1'] for s in bin_sne])
    c = np.array([s['c'] for s in bin_sne])
    
    ss_tot = np.sum((mB - np.mean(mB))**2)
    if ss_tot == 0: continue
    
    # R² from color alone
    A_c = np.column_stack([np.ones(len(bin_sne)), c])
    p_c = np.linalg.lstsq(A_c, mB, rcond=None)[0]
    R2_c = 1 - np.sum((mB - A_c @ p_c)**2) / ss_tot
    
    # R² from stretch alone
    A_x = np.column_stack([np.ones(len(bin_sne)), x1])
    p_x = np.linalg.lstsq(A_x, mB, rcond=None)[0]
    R2_x = 1 - np.sum((mB - A_x @ p_x)**2) / ss_tot
    
    # Conditional variance from both
    A_both = np.column_stack([np.ones(len(bin_sne)), c, x1])
    p_both = np.linalg.lstsq(A_both, mB, rcond=None)[0]
    resid = mB - A_both @ p_both
    cond_var = np.var(resid)
    
    # Leverage: fraction of total standardization from each
    R2_both = 1 - np.sum(resid**2) / ss_tot
    lev_c = R2_c / R2_both if R2_both > 0 else 0
    lev_x = R2_x / R2_both if R2_both > 0 else 0
    
    z_center = np.mean([s['zHD'] for s in bin_sne])
    z_centers.append(z_center)
    R2c_list.append(R2_c)
    R2x_list.append(R2_x)
    condvar_list.append(cond_var)
    leverage_c_list.append(lev_c)
    leverage_x1_list.append(lev_x)
    
    print(f"[{z_lo:.2f},{z_hi:.2f}){'':<3} {len(bin_sne):>4} {R2_c:>7.3f} {R2_x:>7.3f} {cond_var:>12.5f} {lev_c:>10.3f} {lev_x:>11.3f}")

z_arr = np.array(z_centers)
cv_arr = np.array(condvar_list)
lc_arr = np.array(leverage_c_list)
lx_arr = np.array(leverage_x1_list)

rho_cv, p_cv = spearmanr(z_arr, cv_arr)
rho_lc, p_lc = spearmanr(z_arr, lc_arr)
rho_lx, p_lx = spearmanr(z_arr, lx_arr)

print(f"\n  Var(mB|c,x1) vs z:  ρ = {rho_cv:+.3f}, p = {p_cv:.4f}")
print(f"  Leverage_c vs z:    ρ = {rho_lc:+.3f}, p = {p_lc:.4f}")
print(f"  Leverage_x1 vs z:   ρ = {rho_lx:+.3f}, p = {p_lx:.4f}")

if rho_cv < -0.3 and rho_lc < -0.3:
    print(f"\n  🔥 Conditional variance DROPS while color leverage DROPS")
    print(f"     → Population tightens faster than color standardization weakens")
elif abs(rho_cv) < 0.3 and rho_lc < -0.3:
    print(f"\n  🔥 Conditional variance FLAT while color leverage drops")
    print(f"     → Other channels compensate exactly")

# ============================================================
# Q2: Do principal axes of (mB, c, x1) covariance ROTATE with z?
# ============================================================
print(f"\n\n" + "=" * 80)
print("Q2: Principal Axis Rotation of (mB, x1, c) Covariance")
print("=" * 80)

# For each z bin, compute eigenvectors of 3D covariance
# Track the angle of the first principal component relative to the mB axis
# And the projection of each eigenvector onto the c and x1 directions

print(f"\n{'z_range':<12} {'N':>4} | {'PC1→mB':>7} {'PC1→c':>7} {'PC1→x1':>7} | {'PC2→mB':>7} {'PC2→c':>7} {'PC2→x1':>7} | {'λ1/λ3':>7}")
print("-" * 90)

pc1_c_list = []
pc1_x1_list = []
pc1_mB_list = []
anisotropy_list = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 15:
        continue
    
    mB = np.array([s['mB'] for s in bin_sne])
    x1 = np.array([s['x1'] for s in bin_sne])
    c = np.array([s['c'] for s in bin_sne])
    
    # Standardize before PCA
    mB_s = (mB - np.mean(mB)) / (np.std(mB) + 1e-10)
    x1_s = (x1 - np.mean(x1)) / (np.std(x1) + 1e-10)
    c_s = (c - np.mean(c)) / (np.std(c) + 1e-10)
    
    features = np.column_stack([mB_s, x1_s, c_s])
    C = np.cov(features, rowvar=False)
    
    eigenvalues, eigenvectors = np.linalg.eigh(C)
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    # PC1 and PC2 projections onto (mB, x1, c) axes
    pc1 = eigenvectors[:, 0]  # [mB, x1, c] components
    pc2 = eigenvectors[:, 1]
    
    # Take absolute values (sign ambiguity in eigenvectors)
    pc1 = np.abs(pc1)
    pc2 = np.abs(pc2)
    
    anisotropy = eigenvalues[0] / eigenvalues[2] if eigenvalues[2] > 0 else 999
    
    z_center = np.mean([s['zHD'] for s in bin_sne])
    pc1_mB_list.append(pc1[0])
    pc1_x1_list.append(pc1[1])
    pc1_c_list.append(pc1[2])
    anisotropy_list.append(anisotropy)
    
    print(f"[{z_lo:.2f},{z_hi:.2f}){'':<3} {len(bin_sne):>4} | {pc1[0]:>7.3f} {pc1[2]:>7.3f} {pc1[1]:>7.3f} | {pc2[0]:>7.3f} {pc2[2]:>7.3f} {pc2[1]:>7.3f} | {anisotropy:>7.1f}")

pc1_c = np.array(pc1_c_list)
pc1_x1 = np.array(pc1_x1_list)
pc1_mB = np.array(pc1_mB_list)
aniso = np.array(anisotropy_list)

rho_pc1c, p_pc1c = spearmanr(z_arr[:len(pc1_c)], pc1_c)
rho_pc1x, p_pc1x = spearmanr(z_arr[:len(pc1_x1)], pc1_x1)
rho_pc1m, p_pc1m = spearmanr(z_arr[:len(pc1_mB)], pc1_mB)
rho_aniso, p_aniso = spearmanr(z_arr[:len(aniso)], aniso)

print(f"\n  PC1→color vs z:     ρ = {rho_pc1c:+.3f}, p = {p_pc1c:.4f}")
print(f"  PC1→stretch vs z:   ρ = {rho_pc1x:+.3f}, p = {p_pc1x:.4f}")
print(f"  PC1→mB vs z:        ρ = {rho_pc1m:+.3f}, p = {p_pc1m:.4f}")
print(f"  Anisotropy (λ1/λ3): ρ = {rho_aniso:+.3f}, p = {p_aniso:.4f}")

# ============================================================
# Q3: Does the surviving low-noise axis align with stretch?
# ============================================================
print(f"\n\n" + "=" * 80)
print("Q3: Which axis is the surviving low-noise direction?")
print("=" * 80)

# The LAST eigenvector (smallest eigenvalue) is the low-noise axis
# = the direction the population is TIGHTEST
# If this aligns with mB → distances are geometrically constrained
# If it aligns with color → color has collapsed (diagnostic dead)

print(f"\n{'z_range':<12} {'N':>4} | {'PC3→mB':>7} {'PC3→c':>7} {'PC3→x1':>7} | {'Interpretation'}")
print("-" * 80)

pc3_c_list = []
pc3_x1_list = []
pc3_mB_list = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 15:
        continue
    
    mB = np.array([s['mB'] for s in bin_sne])
    x1 = np.array([s['x1'] for s in bin_sne])
    c = np.array([s['c'] for s in bin_sne])
    
    mB_s = (mB - np.mean(mB)) / (np.std(mB) + 1e-10)
    x1_s = (x1 - np.mean(x1)) / (np.std(x1) + 1e-10)
    c_s = (c - np.mean(c)) / (np.std(c) + 1e-10)
    
    features = np.column_stack([mB_s, x1_s, c_s])
    C = np.cov(features, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(C)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvectors = eigenvectors[:, idx]
    
    pc3 = np.abs(eigenvectors[:, 2])  # smallest eigenvalue direction
    
    # Interpret which axis dominates the tight direction
    dominant = ["mB (geometry)", "x1 (stretch/timing)", "c (color/thermo)"][np.argmax(pc3)]
    
    pc3_mB_list.append(pc3[0])
    pc3_x1_list.append(pc3[1])
    pc3_c_list.append(pc3[2])
    
    print(f"[{z_lo:.2f},{z_hi:.2f}){'':<3} {len(bin_sne):>4} | {pc3[0]:>7.3f} {pc3[2]:>7.3f} {pc3[1]:>7.3f} | tight axis → {dominant}")

pc3_c = np.array(pc3_c_list)
pc3_x1 = np.array(pc3_x1_list)

rho_tight_c, p_tight_c = spearmanr(z_arr[:len(pc3_c)], pc3_c)
rho_tight_x, p_tight_x = spearmanr(z_arr[:len(pc3_x1)], pc3_x1)

print(f"\n  Tight axis → color vs z:   ρ = {rho_tight_c:+.3f}, p = {p_tight_c:.4f}")
print(f"  Tight axis → stretch vs z: ρ = {rho_tight_x:+.3f}, p = {p_tight_x:.4f}")

if rho_tight_c > 0.3:
    print(f"\n  🔥 At high z, the TIGHTEST direction is COLOR")
    print(f"     → Color variation collapses → diagnostic channel DEAD")
    print(f"     → Population converges in thermodynamic properties")
elif rho_tight_x > 0.3:
    print(f"\n  At high z, the tightest direction is STRETCH")
    print(f"  → Geometric timing becomes the constrained direction")

# ============================================================
# Q4: Does foreground change rotation rate?
# ============================================================
print(f"\n\n" + "=" * 80)
print("Q4: Does Foreground Environment Change the Rotation?")
print("=" * 80)

# Split by MWEBV and host mass, check if PC structure differs
for z_lo, z_hi in [(0.01, 0.3), (0.3, 0.82), (0.01, 0.82)]:
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 40:
        continue
    
    mwebv_med = np.median([s['MWEBV'] for s in bin_sne])
    
    for label, subset in [("Clean (low MW)", [s for s in bin_sne if s['MWEBV'] <= mwebv_med]),
                           ("Dusty (high MW)", [s for s in bin_sne if s['MWEBV'] > mwebv_med])]:
        if len(subset) < 15:
            continue
        
        mB = np.array([s['mB'] for s in subset])
        x1 = np.array([s['x1'] for s in subset])
        c = np.array([s['c'] for s in subset])
        
        mB_s = (mB - np.mean(mB)) / (np.std(mB) + 1e-10)
        x1_s = (x1 - np.mean(x1)) / (np.std(x1) + 1e-10)
        c_s = (c - np.mean(c)) / (np.std(c) + 1e-10)
        
        features = np.column_stack([mB_s, x1_s, c_s])
        C = np.cov(features, rowvar=False)
        eigenvalues, eigenvectors = np.linalg.eigh(C)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        
        pc1 = np.abs(eigenvectors[:, 0])
        pc3 = np.abs(eigenvectors[:, 2])
        aniso = eigenvalues[0] / eigenvalues[2]
        
        print(f"  [{z_lo:.2f},{z_hi:.2f}) {label:<18} N={len(subset):>4} | PC1: mB={pc1[0]:.2f} c={pc1[2]:.2f} x1={pc1[1]:.2f} | PC3: mB={pc3[0]:.2f} c={pc3[2]:.2f} x1={pc3[1]:.2f} | λ1/λ3={aniso:.1f}")
    
    # Also split by host mass
    masses = [s['HOST_LOGMASS'] for s in bin_sne if s.get('HOST_LOGMASS', 0) > 0]
    if masses:
        mass_med = np.median(masses)
        for label, subset in [("Low mass host", [s for s in bin_sne if s.get('HOST_LOGMASS', 0) > 0 and s['HOST_LOGMASS'] <= mass_med]),
                               ("High mass host", [s for s in bin_sne if s.get('HOST_LOGMASS', 0) > 0 and s['HOST_LOGMASS'] > mass_med])]:
            if len(subset) < 15:
                continue
            
            mB = np.array([s['mB'] for s in subset])
            x1 = np.array([s['x1'] for s in subset])
            c = np.array([s['c'] for s in subset])
            
            mB_s = (mB - np.mean(mB)) / (np.std(mB) + 1e-10)
            x1_s = (x1 - np.mean(x1)) / (np.std(x1) + 1e-10)
            c_s = (c - np.mean(c)) / (np.std(c) + 1e-10)
            
            features = np.column_stack([mB_s, x1_s, c_s])
            C = np.cov(features, rowvar=False)
            eigenvalues, eigenvectors = np.linalg.eigh(C)
            idx = np.argsort(eigenvalues)[::-1]
            eigenvalues = eigenvalues[idx]
            eigenvectors = eigenvectors[:, idx]
            
            pc1 = np.abs(eigenvectors[:, 0])
            pc3 = np.abs(eigenvectors[:, 2])
            aniso = eigenvalues[0] / eigenvalues[2]
            
            print(f"  [{z_lo:.2f},{z_hi:.2f}) {label:<18} N={len(subset):>4} | PC1: mB={pc1[0]:.2f} c={pc1[2]:.2f} x1={pc1[1]:.2f} | PC3: mB={pc3[0]:.2f} c={pc3[2]:.2f} x1={pc3[1]:.2f} | λ1/λ3={aniso:.1f}")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_covariance_rotation")
results_dir.mkdir(exist_ok=True)

results = {
    'z_centers': [float(z) for z in z_arr],
    'condvar': [float(c) for c in condvar_list],
    'leverage_c': [float(l) for l in leverage_c_list],
    'leverage_x1': [float(l) for l in leverage_x1_list],
    'pc1_color': [float(p) for p in pc1_c_list],
    'pc1_stretch': [float(p) for p in pc1_x1_list],
    'pc3_color': [float(p) for p in pc3_c_list],
    'pc3_stretch': [float(p) for p in pc3_x1_list],
}

with open(results_dir / "covariance_rotation.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {results_dir}/")
print("\n" + "=" * 80)
print("SUMMARY — GPT's Four Questions Answered")
print("=" * 80)
print(f"  Q1: Var(mB|c,x1) vs z:    ρ = {rho_cv:+.3f} (p={p_cv:.3f})")
print(f"      Color leverage vs z:   ρ = {rho_lc:+.3f} (p={p_lc:.3f})")
print(f"  Q2: PC1→color vs z:       ρ = {rho_pc1c:+.3f} (p={p_pc1c:.3f})")
print(f"      Anisotropy vs z:       ρ = {rho_aniso:+.3f} (p={p_aniso:.3f})")
print(f"  Q3: Tight axis→color vs z: ρ = {rho_tight_c:+.3f} (p={p_tight_c:.3f})")
print(f"      Tight axis→stretch:    ρ = {rho_tight_x:+.3f} (p={p_tight_x:.3f})")
print(f"  Q4: See environment splits above")
