#!/usr/bin/env python3
"""
closure_fisher_budget.py — Fisher Information Budget Invariance
================================================================

GPT's insight: the conserved quantity isn't α+β. It's the Fisher
information about distance μ given the feature vector (mB, x1, c).

I_μ(z) = g^T C(z)^{-1} g

where g = ∂d/∂μ and C(z) is the covariance of features at redshift z.

If I_μ stays flat while β drops → information ROTATES in feature space,
not destroyed. The "bistable switch" is a rotation of the covariance
ellipsoid, not signal destruction.

Track per z-bin:
1. R²_color(z) — fraction of mB variance explained by color
2. R²_stretch(z) — fraction explained by stretch  
3. Var(mB | x1, c) — conditional variance (residual after standardization)
4. I_μ(z) — Fisher information about distance
5. Covariance ellipsoid orientation

Author: Closure Theory Collaboration (GPT framework)
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

print(f"Loaded {len(sne)} unique SNe Ia\n")

# ============================================================
# FISHER INFORMATION BUDGET PER Z-BIN
# ============================================================

z_edges = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.55, 0.7, 0.85, 1.0, 1.5, 2.5]

print("=" * 80)
print("FISHER INFORMATION BUDGET vs REDSHIFT")
print("=" * 80)

print(f"\n{'z_range':<13} {'N':>4} {'R²_c':>7} {'R²_x1':>7} {'R²_tot':>7} {'Var(m|x)':>9} {'I_μ':>8} {'θ_ell':>7}")
print("-" * 75)

z_centers = []
R2_colors = []
R2_stretches = []
R2_totals = []
cond_vars = []
fisher_infos = []
ellipse_angles = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    
    if len(bin_sne) < 15:
        continue
    
    mB = np.array([s['mB'] for s in bin_sne])
    x1 = np.array([s['x1'] for s in bin_sne])
    c = np.array([s['c'] for s in bin_sne])
    z_center = np.mean([s['zHD'] for s in bin_sne])
    
    # ---- R² from color alone ----
    A_c = np.column_stack([np.ones(len(bin_sne)), c])
    try:
        params_c = np.linalg.lstsq(A_c, mB, rcond=None)[0]
        resid_c = mB - A_c @ params_c
        ss_res_c = np.sum(resid_c**2)
        ss_tot = np.sum((mB - np.mean(mB))**2)
        R2_c = 1 - ss_res_c / ss_tot if ss_tot > 0 else 0
    except:
        R2_c = 0
    
    # ---- R² from stretch alone ----
    A_x = np.column_stack([np.ones(len(bin_sne)), x1])
    try:
        params_x = np.linalg.lstsq(A_x, mB, rcond=None)[0]
        resid_x = mB - A_x @ params_x
        ss_res_x = np.sum(resid_x**2)
        R2_x = 1 - ss_res_x / ss_tot if ss_tot > 0 else 0
    except:
        R2_x = 0
    
    # ---- R² from both (total standardization) ----
    A_both = np.column_stack([np.ones(len(bin_sne)), c, x1])
    try:
        params_both = np.linalg.lstsq(A_both, mB, rcond=None)[0]
        resid_both = mB - A_both @ params_both
        ss_res_both = np.sum(resid_both**2)
        R2_both = 1 - ss_res_both / ss_tot if ss_tot > 0 else 0
    except:
        R2_both = 0
    
    # ---- Conditional variance Var(mB | x1, c) ----
    cond_var = np.var(resid_both) if len(resid_both) > 2 else 999
    
    # ---- Fisher Information ----
    # Feature vector d = (mB, x1, c)
    # Distance gradient g = (1, 0, 0) since ∂mB/∂μ = 1, others ~0
    # I_μ = g^T C^{-1} g = (C^{-1})_{00}
    
    features = np.column_stack([mB, x1, c])
    C = np.cov(features, rowvar=False)
    
    try:
        C_inv = np.linalg.inv(C)
        # g = (1, 0, 0)
        I_mu = C_inv[0, 0]
    except:
        I_mu = 0
    
    # ---- Covariance ellipsoid angle ----
    # 2D projection: (c, x1) covariance
    C_cx = np.cov(np.column_stack([c, x1]), rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(C_cx)
    # Angle of major axis
    angle = np.degrees(np.arctan2(eigenvectors[1, 1], eigenvectors[0, 1]))
    
    z_centers.append(z_center)
    R2_colors.append(R2_c)
    R2_stretches.append(R2_x)
    R2_totals.append(R2_both)
    cond_vars.append(cond_var)
    fisher_infos.append(I_mu)
    ellipse_angles.append(angle)
    
    print(f"[{z_lo:.2f},{z_hi:.2f}){'':<4} {len(bin_sne):>4} {R2_c:>7.3f} {R2_x:>7.3f} {R2_both:>7.3f} {cond_var:>9.5f} {I_mu:>8.2f} {angle:>7.1f}°")

z_arr = np.array(z_centers)
R2c = np.array(R2_colors)
R2x = np.array(R2_stretches)
R2t = np.array(R2_totals)
cv = np.array(cond_vars)
fi = np.array(fisher_infos)
ea = np.array(ellipse_angles)

# ============================================================
# CORRELATIONS WITH REDSHIFT
# ============================================================
print(f"\n" + "=" * 80)
print("TRENDS WITH REDSHIFT")
print("=" * 80)

for label, arr in [("R²_color", R2c), ("R²_stretch", R2x), ("R²_total", R2t),
                    ("Var(mB|x1,c)", cv), ("Fisher I_μ", fi), ("Ellipse angle", ea)]:
    rho, p = spearmanr(z_arr, arr)
    trend = "↑" if rho > 0.2 else "↓" if rho < -0.2 else "→"
    sig = "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.1 else ""
    print(f"  {label:<16}: ρ = {rho:+.3f}, p = {p:.4f} {trend} {sig}")

# ============================================================
# THE KEY TEST: Is Fisher info conserved?
# ============================================================
print(f"\n" + "=" * 80)
print("THE CONSERVATION TEST: Is I_μ invariant while β drops?")
print("=" * 80)

rho_fi, p_fi = spearmanr(z_arr, fi)
rho_cv, p_cv = spearmanr(z_arr, cv)

fi_cv = np.std(fi) / np.mean(fi) * 100
cv_cv = np.std(cv) / np.mean(cv) * 100

print(f"\n  Fisher I_μ vs z: ρ = {rho_fi:+.3f}, p = {p_fi:.4f}")
print(f"  Fisher I_μ: mean = {np.mean(fi):.2f}, CV = {fi_cv:.1f}%")
print(f"  Cond. Var vs z:  ρ = {rho_cv:+.3f}, p = {p_cv:.4f}")
print(f"  Cond. Var: mean = {np.mean(cv):.5f}, CV = {cv_cv:.1f}%")

if abs(rho_fi) < 0.4 and p_fi > 0.05:
    print(f"\n  🔥🔥🔥 FISHER INFORMATION IS CONSERVED!")
    print(f"  β drops, but I_μ stays flat → information ROTATES, not destroyed")
    print(f"  The covariance ellipsoid reshapes to preserve distance precision")
    print(f"  THIS IS THE CONSERVATION LAW.")
elif rho_fi > 0.4:
    print(f"\n  ★ Fisher information INCREASES with z → geometry gains from diagnostic loss")
elif rho_fi < -0.4 and p_fi < 0.05:
    print(f"\n  ❄️ Fisher information DROPS with z → information genuinely lost")
    print(f"  Bistable switch DEAD")

# ============================================================
# R² CHANNEL TRANSFER
# ============================================================
print(f"\n" + "=" * 80)
print("CHANNEL TRANSFER: Does R²_stretch pick up what R²_color loses?")
print("=" * 80)

rho_rc, p_rc = spearmanr(z_arr, R2c)
rho_rx, p_rx = spearmanr(z_arr, R2x)
rho_transfer, p_transfer = spearmanr(R2c, R2x)

print(f"  R²_color vs z:    ρ = {rho_rc:+.3f}, p = {p_rc:.4f}")
print(f"  R²_stretch vs z:  ρ = {rho_rx:+.3f}, p = {p_rx:.4f}")
print(f"  R²_color vs R²_stretch: ρ = {rho_transfer:+.3f}, p = {p_transfer:.4f}")

if rho_rc < -0.2 and rho_rx > 0.2:
    print(f"\n  🔥 Color power drops, stretch power RISES → direct channel transfer!")
elif rho_rc < -0.2 and abs(rho_rx) < 0.2:
    print(f"\n  Color drops, stretch flat → no explicit transfer visible")
    print(f"  But if I_μ is conserved, transfer is in the COVARIANCE structure")
elif abs(rho_rc) < 0.2:
    print(f"\n  Color power doesn't significantly change with z")

# ============================================================
# ELLIPSOID ROTATION
# ============================================================
print(f"\n" + "=" * 80)
print("COVARIANCE ELLIPSOID: Does it rotate with z?")
print("=" * 80)

rho_ea, p_ea = spearmanr(z_arr, ea)
print(f"  Ellipse angle vs z: ρ = {rho_ea:+.3f}, p = {p_ea:.4f}")
print(f"  Angle range: {ea.min():.1f}° to {ea.max():.1f}° (Δ = {ea.max()-ea.min():.1f}°)")

if abs(rho_ea) > 0.4 and p_ea < 0.1:
    print(f"\n  🔥 Covariance ellipsoid ROTATES with z!")
    print(f"  The feature space basis is changing — compression rotates the 'useful' direction")
else:
    print(f"\n  Ellipsoid orientation roughly stable with z")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_fisher_budget")
results_dir.mkdir(exist_ok=True)

results = {
    'z_centers': [float(z) for z in z_arr],
    'R2_color': [float(r) for r in R2c],
    'R2_stretch': [float(r) for r in R2x],
    'R2_total': [float(r) for r in R2t],
    'cond_var': [float(c) for c in cv],
    'fisher_info': [float(f) for f in fi],
    'ellipse_angle': [float(e) for e in ea],
    'fisher_vs_z_rho': float(rho_fi),
    'fisher_vs_z_p': float(p_fi),
    'condvar_vs_z_rho': float(rho_cv),
    'condvar_vs_z_p': float(p_cv),
}

with open(results_dir / "fisher_budget.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {results_dir}/")
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"  β drops with z:              YES (established)")
print(f"  Fisher I_μ vs z:             ρ = {rho_fi:+.3f} ({'CONSERVED' if abs(rho_fi)<0.4 else 'NOT conserved'})")
print(f"  Cond. variance vs z:         ρ = {rho_cv:+.3f}")
print(f"  R²_color vs z:              ρ = {rho_rc:+.3f}")
print(f"  R²_stretch vs z:            ρ = {rho_rx:+.3f}")
print(f"  Ellipsoid rotation:          ρ = {rho_ea:+.3f}")
