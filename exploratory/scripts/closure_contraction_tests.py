#!/usr/bin/env python3
"""
closure_contraction_tests.py — The Fork in the Road
=====================================================

GPT's three tests to separate:
A. Population contraction (less diversity at high z)
B. Propagation compression (medium acts on signal)

TEST 1: Eigenvalue contraction law — fit λ_i(z) to exponential decay
TEST 2: Orthogonal subspace — does shrinkage happen in (c,x1) not mB?
TEST 3: Within-class persistence — control for host properties

The decisive question: does contraction persist at FIXED host type?
If yes → epoch/propagation. If no → demographics.

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr
from scipy.optimize import curve_fit
from scipy.integrate import quad
import json
from pathlib import Path

np.random.seed(42)

# Cosmology for Sigma(z)
H0_cgs = 67.4e5 / 3.086e24
Omega_m, Omega_L, Omega_b = 0.315, 0.685, 0.0493
m_p = 1.6726e-24
rho_crit = 1.878e-29
n_H0 = Omega_b * rho_crit * (67.4/100)**2 / m_p * 0.76

def E(z): return np.sqrt(Omega_m*(1+z)**3 + Omega_L)
def Sigma(z):
    if z <= 0: return 0
    r, _ = quad(lambda zz: n_H0*(1+zz)**2*2.998e10/(H0_cgs*E(zz)), 0, z)
    return r

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

print(f"Loaded {len(sne)} unique SNe Ia\n")

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.82, 1.0, 2.5]

# ============================================================
# TEST 1: EIGENVALUE CONTRACTION LAW
# ============================================================
print("=" * 80)
print("TEST 1: Eigenvalue Contraction — λ_i(z) vs Σ(z)")
print("=" * 80)

# Compute standardized covariance eigenvalues per z-bin
z_centers = []
eigenval_sets = []  # list of [λ1, λ2, λ3] per bin
sigma_vals = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 20:
        continue
    
    mB = np.array([s['mB'] for s in bin_sne])
    x1 = np.array([s['x1'] for s in bin_sne])
    c = np.array([s['c'] for s in bin_sne])
    
    # RAW covariance (not standardized — we want to see actual scatter)
    C = np.cov(np.column_stack([mB, x1, c]), rowvar=False)
    eigvals = np.sort(np.linalg.eigvalsh(C))[::-1]
    
    z_center = np.mean([s['zHD'] for s in bin_sne])
    sig = Sigma(z_center)
    
    z_centers.append(z_center)
    eigenval_sets.append(eigvals)
    sigma_vals.append(sig)

z_arr = np.array(z_centers)
sig_arr = np.array(sigma_vals)
ev = np.array(eigenval_sets)

print(f"\n{'z':>6} {'Σ(z)':>12} {'λ₁ (mB)':>10} {'λ₂':>10} {'λ₃ (tight)':>12} {'λ₁/λ₃':>8}")
print("-" * 65)
for i in range(len(z_arr)):
    print(f"{z_arr[i]:>6.2f} {sig_arr[i]:>12.2e} {ev[i,0]:>10.4f} {ev[i,1]:>10.4f} {ev[i,2]:>12.6f} {ev[i,0]/ev[i,2]:>8.1f}")

# Fit exponential contraction: λ(Σ) = λ_∞ + (λ_0 - λ_∞) exp(-2k Σ)
def exp_contraction(sig, lam0, lam_inf, k):
    return lam_inf + (lam0 - lam_inf) * np.exp(-2 * k * sig)

print(f"\n  Exponential contraction fits: λ_i(Σ) = λ_∞ + (λ_0 − λ_∞)·exp(−2k·Σ)")
for j, label in enumerate(["λ₁ (largest)", "λ₂ (middle)", "λ₃ (smallest)"]):
    lam = ev[:, j]
    rho, p = spearmanr(sig_arr, lam)
    
    try:
        popt, _ = curve_fit(exp_contraction, sig_arr, lam, 
                           p0=[lam[0], lam[-1], 1e-22],
                           bounds=([0, 0, 0], [100, 100, 1e-18]),
                           maxfev=5000)
        lam0, lam_inf, k_fit = popt
        pred = exp_contraction(sig_arr, *popt)
        ss_res = np.sum((lam - pred)**2)
        ss_tot = np.sum((lam - np.mean(lam))**2)
        r_sq = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        
        print(f"  {label}: λ₀={lam0:.4f}, λ_∞={lam_inf:.6f}, k={k_fit:.2e}, R²={r_sq:.3f}")
        print(f"           ρ(λ,Σ) = {rho:+.3f}, p = {p:.4f}")
    except Exception as e:
        print(f"  {label}: fit failed ({e}), ρ(λ,Σ) = {rho:+.3f}, p = {p:.4f}")

# ============================================================
# TEST 2: ORTHOGONAL SUBSPACE — Where does the shrinkage happen?
# ============================================================
print(f"\n\n" + "=" * 80)
print("TEST 2: Where Does Shrinkage Happen? mB vs (c, x1) Subspace")
print("=" * 80)

print(f"\n{'z':>6} {'Var(mB)':>10} {'Var(c)':>10} {'Var(x1)':>10} {'Var_⊥':>10} {'Var_mB/Var_⊥':>13}")
print("-" * 65)

var_mB_list = []
var_c_list = []
var_x1_list = []
var_perp_list = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 20:
        continue
    
    mB = np.array([s['mB'] for s in bin_sne])
    x1 = np.array([s['x1'] for s in bin_sne])
    c = np.array([s['c'] for s in bin_sne])
    
    v_mB = np.var(mB)
    v_c = np.var(c)
    v_x1 = np.var(x1)
    v_perp = v_c + v_x1  # nuisance subspace variance
    
    z_center = np.mean([s['zHD'] for s in bin_sne])
    var_mB_list.append(v_mB)
    var_c_list.append(v_c)
    var_x1_list.append(v_x1)
    var_perp_list.append(v_perp)
    
    ratio = v_mB / v_perp if v_perp > 0 else 999
    print(f"{z_center:>6.2f} {v_mB:>10.4f} {v_c:>10.6f} {v_x1:>10.4f} {v_perp:>10.4f} {ratio:>13.3f}")

vm = np.array(var_mB_list)
vc = np.array(var_c_list)
vx = np.array(var_x1_list)
vp = np.array(var_perp_list)

rho_vm, p_vm = spearmanr(z_arr[:len(vm)], vm)
rho_vc, p_vc = spearmanr(z_arr[:len(vc)], vc)
rho_vx, p_vx = spearmanr(z_arr[:len(vx)], vx)
rho_vp, p_vp = spearmanr(z_arr[:len(vp)], vp)

print(f"\n  Var(mB) vs z:  ρ = {rho_vm:+.3f}, p = {p_vm:.4f}")
print(f"  Var(c) vs z:   ρ = {rho_vc:+.3f}, p = {p_vc:.4f}")
print(f"  Var(x1) vs z:  ρ = {rho_vx:+.3f}, p = {p_vx:.4f}")
print(f"  Var_⊥ vs z:    ρ = {rho_vp:+.3f}, p = {p_vp:.4f}")

if abs(rho_vp) > abs(rho_vm) and rho_vp < 0:
    print(f"\n  🔥 Nuisance subspace (c, x1) shrinks FASTER than mB")
    print(f"     → Population homogenizes in thermodynamic properties")
    print(f"     → Distances stay diverse while source states converge")
elif rho_vm < -0.3 and rho_vp < -0.3:
    print(f"\n  Both mB and nuisance shrink — global contraction")

# ============================================================
# TEST 3: WITHIN-CLASS PERSISTENCE — The decisive fork
# ============================================================
print(f"\n\n" + "=" * 80)
print("TEST 3: Within-Class Persistence — Source Evolution vs Propagation")
print("=" * 80)

# Split by host mass (low vs high), then check if contraction persists WITHIN each class
# If it persists → epoch/propagation effect
# If it disappears → demographics (different host populations at different z)

mass_med = np.median([s['HOST_LOGMASS'] for s in sne if s.get('HOST_LOGMASS', 0) > 5])

# Also split by stretch (x1 > 0 vs x1 < 0) as a light-curve subclass
# And by color (c > 0 vs c < 0)

splits = [
    ("Low mass (M<{:.1f})".format(mass_med), 
     lambda s: s.get('HOST_LOGMASS', 0) > 5 and s['HOST_LOGMASS'] <= mass_med),
    ("High mass (M>{:.1f})".format(mass_med),
     lambda s: s.get('HOST_LOGMASS', 0) > 5 and s['HOST_LOGMASS'] > mass_med),
    ("Fast decliners (x1<0)",
     lambda s: s['x1'] < 0),
    ("Slow decliners (x1>0)",
     lambda s: s['x1'] > 0),
    ("Blue (c<0)",
     lambda s: s['c'] < 0),
    ("Red (c>0)",
     lambda s: s['c'] > 0),
]

print(f"\n  Testing: does Var(c) + Var(x1) still shrink with z within each subclass?")
print(f"\n  {'Subclass':<28} {'N':>5} {'ρ(Var_⊥,z)':>12} {'p':>8} {'Persists?':<12}")
print(f"  " + "-" * 70)

persistence_results = []

for label, cut_fn in splits:
    subset = [s for s in sne if cut_fn(s)]
    
    if len(subset) < 50:
        continue
    
    # Compute Var_perp in z bins within this subset
    sub_z = []
    sub_vp = []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_sub = [s for s in subset if z_lo <= s['zHD'] < z_hi]
        if len(bin_sub) < 10:
            continue
        
        c_vals = np.array([s['c'] for s in bin_sub])
        x1_vals = np.array([s['x1'] for s in bin_sub])
        vp = np.var(c_vals) + np.var(x1_vals)
        
        sub_z.append(np.mean([s['zHD'] for s in bin_sub]))
        sub_vp.append(vp)
    
    if len(sub_z) >= 4:
        rho_sub, p_sub = spearmanr(sub_z, sub_vp)
        persists = "YES ✓" if rho_sub < -0.3 else "NO ✗" if rho_sub > 0.3 else "UNCLEAR"
        print(f"  {label:<28} {len(subset):>5} {rho_sub:>+12.3f} {p_sub:>8.4f} {persists:<12}")
        
        persistence_results.append({
            'label': label,
            'N': len(subset),
            'rho': float(rho_sub),
            'p': float(p_sub),
            'persists': rho_sub < -0.3,
        })

# Count
n_persist = sum(1 for r in persistence_results if r['persists'])
n_total = len(persistence_results)

print(f"\n  VERDICT: {n_persist}/{n_total} subclasses show persistent contraction")

if n_persist >= n_total * 0.7:
    print(f"\n  🔥🔥🔥 CONTRACTION PERSISTS WITHIN CONTROLLED SUBCLASSES")
    print(f"  → NOT demographics or population mixing")
    print(f"  → Epoch effect or propagation compression")
    print(f"  → Humza's 'hard cap' framework SUPPORTED")
elif n_persist >= n_total * 0.4:
    print(f"\n  ⚠️ Mixed results — some persist, some don't")
    print(f"  → Likely a COMBINATION of source evolution + epoch effect")
elif n_persist < n_total * 0.4:
    print(f"\n  ❄️ Contraction mostly disappears within subclasses")
    print(f"  → Primarily demographic/population mixing")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_contraction")
results_dir.mkdir(exist_ok=True)

results = {
    'eigenvalue_trends': {
        'z_centers': [float(z) for z in z_arr],
        'sigma_vals': [float(s) for s in sig_arr],
        'eigenvalues': ev.tolist(),
    },
    'variance_decomposition': {
        'var_mB_rho': float(rho_vm),
        'var_c_rho': float(rho_vc),
        'var_x1_rho': float(rho_vx),
        'var_perp_rho': float(rho_vp),
    },
    'persistence': persistence_results,
    'n_persist': n_persist,
    'n_total': n_total,
}

with open(results_dir / "contraction_tests.json", 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults saved to {results_dir}/")
print("\n" + "=" * 80)
print("THE FORK: Source Evolution vs Propagation")
print("=" * 80)
print(f"  Eigenvalue contraction: ρ(λ₁,z) = {spearmanr(z_arr, ev[:,0])[0]:+.3f}")
print(f"  Nuisance shrinks faster: Var_⊥ ρ = {rho_vp:+.3f} vs Var(mB) ρ = {rho_vm:+.3f}")
print(f"  Within-class persistence: {n_persist}/{n_total}")
print(f"  → {'EPOCH/PROPAGATION' if n_persist >= n_total*0.7 else 'MIXED' if n_persist >= n_total*0.4 else 'DEMOGRAPHICS'}")
