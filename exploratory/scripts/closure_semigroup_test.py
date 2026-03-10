#!/usr/bin/env python3
"""
closure_semigroup_test.py — Contraction Semigroup Test
=======================================================

GPT's model: C(τ) = exp(-Kτ) C_0 exp(-K^T τ)

If K is diagonal: λ_i(τ) = λ_i(0) exp(-2k_i τ)

Test: do the INTRINSIC variances (after error subtraction)
follow exponential decay in z (or Σ, or lookback time)?

If yes: contraction semigroup is the correct mathematical description.
The contraction rates k_mB, k_x1, k_c tell us which channels
contract fastest.

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

# Cosmology
H0 = 67.4
Omega_m, Omega_L = 0.315, 0.685
Omega_b = 0.0493
m_p = 1.6726e-24
rho_crit = 1.878e-29
n_H0 = Omega_b * rho_crit * (H0/100)**2 / m_p * 0.76
H0_cgs = H0 * 1e5 / 3.086e24
c_cgs = 2.998e10

def E(z): return np.sqrt(Omega_m*(1+z)**3 + Omega_L)

def Sigma(z):
    if z <= 0: return 0
    r, _ = quad(lambda zz: n_H0*(1+zz)**2*c_cgs/(H0_cgs*E(zz)), 0, z)
    return r

def lookback_time(z):
    """Lookback time in Gyr"""
    r, _ = quad(lambda zz: 1/((1+zz)*E(zz)), 0, z)
    return r / (H0 * 1e5 / 3.086e24 / 3.156e16)  # convert to Gyr

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

# Compute intrinsic variances per z-bin
z_centers = []
var_mB_intr = []
var_x1_intr = []
var_c_intr = []
sig_vals = []
t_lookback = []

print(f"{'z':>6} {'t_lb(Gyr)':>10} {'Σ':>12} | {'Var_i(mB)':>10} {'Var_i(x1)':>10} {'Var_i(c)':>10}")
print("-" * 70)

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 15:
        continue
    
    mB = np.array([s['mB'] for s in bin_sne])
    x1 = np.array([s['x1'] for s in bin_sne])
    c = np.array([s['c'] for s in bin_sne])
    mBerr = np.array([s['mBERR'] for s in bin_sne])
    x1err = np.array([s['x1ERR'] for s in bin_sne])
    cerr = np.array([s['cERR'] for s in bin_sne])
    
    # Intrinsic = observed - mean measurement error
    vi_mB = max(0.001, np.var(mB) - np.mean(mBerr**2))
    vi_x1 = max(0.001, np.var(x1) - np.mean(x1err**2))
    vi_c = max(0.001, np.var(c) - np.mean(cerr**2))
    
    zc = np.mean([s['zHD'] for s in bin_sne])
    sig = Sigma(zc)
    tlb = lookback_time(zc)
    
    z_centers.append(zc)
    var_mB_intr.append(vi_mB)
    var_x1_intr.append(vi_x1)
    var_c_intr.append(vi_c)
    sig_vals.append(sig)
    t_lookback.append(tlb)
    
    print(f"{zc:>6.2f} {tlb:>10.2f} {sig:>12.2e} | {vi_mB:>10.4f} {vi_x1:>10.4f} {vi_c:>10.6f}")

z_arr = np.array(z_centers)
sig_arr = np.array(sig_vals)
t_arr = np.array(t_lookback)
vm = np.array(var_mB_intr)
vx = np.array(var_x1_intr)
vc = np.array(var_c_intr)

# ============================================================
# FIT: λ(τ) = λ_∞ + (λ_0 - λ_∞) exp(-2k τ)
# ============================================================
print(f"\n" + "=" * 80)
print("CONTRACTION SEMIGROUP FITS")
print("=" * 80)

def exp_decay(tau, lam0, lam_inf, k):
    return lam_inf + (lam0 - lam_inf) * np.exp(-2 * k * tau)

def simple_exp(tau, lam0, k):
    return lam0 * np.exp(-2 * k * tau)

# Try three τ variables: z, Σ(z), lookback time
for tau_name, tau_arr in [("z", z_arr), ("Σ(z)", sig_arr), ("t_lookback", t_arr)]:
    print(f"\n  τ = {tau_name}:")
    print(f"  {'Variable':<12} {'λ₀':>8} {'λ_∞':>8} {'k':>12} {'R²':>8} {'R²_simple':>10}")
    print(f"  " + "-" * 60)
    
    for label, var_arr in [("Var_i(mB)", vm), ("Var_i(x1)", vx), ("Var_i(c)", vc)]:
        # Full model: λ_∞ + (λ_0 - λ_∞) exp(-2k τ)
        try:
            if tau_name == "Σ(z)":
                p0 = [var_arr[0], var_arr[-1], 1e-22]
                bounds = ([0, 0, 0], [10, 10, 1e-18])
            elif tau_name == "t_lookback":
                p0 = [var_arr[0], var_arr[-1], 0.1]
                bounds = ([0, 0, 0], [10, 10, 10])
            else:
                p0 = [var_arr[0], var_arr[-1], 1.0]
                bounds = ([0, 0, 0], [10, 10, 50])
            
            popt, _ = curve_fit(exp_decay, tau_arr, var_arr, p0=p0, bounds=bounds, maxfev=10000)
            pred = exp_decay(tau_arr, *popt)
            ss_res = np.sum((var_arr - pred)**2)
            ss_tot = np.sum((var_arr - np.mean(var_arr))**2)
            r_sq = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        except:
            popt = [0, 0, 0]
            r_sq = -999
        
        # Simple model: λ_0 exp(-2k τ)
        try:
            if tau_name == "Σ(z)":
                p0s = [var_arr[0], 1e-22]
                boundss = ([0, 0], [10, 1e-18])
            elif tau_name == "t_lookback":
                p0s = [var_arr[0], 0.1]
                boundss = ([0, 0], [10, 10])
            else:
                p0s = [var_arr[0], 1.0]
                boundss = ([0, 0], [10, 50])
            
            popt_s, _ = curve_fit(simple_exp, tau_arr, var_arr, p0=p0s, bounds=boundss, maxfev=10000)
            pred_s = simple_exp(tau_arr, *popt_s)
            ss_res_s = np.sum((var_arr - pred_s)**2)
            r_sq_s = 1 - ss_res_s/ss_tot if ss_tot > 0 else 0
        except:
            r_sq_s = -999
        
        print(f"  {label:<12} {popt[0]:>8.4f} {popt[1]:>8.4f} {popt[2]:>12.4e} {r_sq:>8.3f} {r_sq_s:>10.3f}")

# ============================================================
# CONTRACTION RATE COMPARISON
# ============================================================
print(f"\n" + "=" * 80)
print("CONTRACTION RATES: Which channel contracts fastest?")
print("=" * 80)

# Use z as τ (simplest, most intuitive)
rates = {}
for label, var_arr in [("mB", vm), ("x1", vx), ("c", vc)]:
    try:
        popt, _ = curve_fit(exp_decay, z_arr, var_arr, 
                           p0=[var_arr[0], var_arr[-1], 1.0],
                           bounds=([0, 0, 0], [10, 10, 50]),
                           maxfev=10000)
        rates[label] = {'lam0': popt[0], 'lam_inf': popt[1], 'k': popt[2]}
        
        pred = exp_decay(z_arr, *popt)
        ss_res = np.sum((var_arr - pred)**2)
        ss_tot = np.sum((var_arr - np.mean(var_arr))**2)
        r_sq = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        rates[label]['r_sq'] = r_sq
    except:
        rates[label] = {'k': 0, 'r_sq': -1}

print(f"\n  Channel  k (contraction rate)   R²      Half-life (Δz)")
print(f"  " + "-" * 55)
for label in ["mB", "x1", "c"]:
    k = rates[label]['k']
    r_sq = rates[label]['r_sq']
    half_life = np.log(2) / (2 * k) if k > 0 else float('inf')
    print(f"  {label:<8} {k:>20.3f}   {r_sq:>6.3f}   {half_life:>10.2f}")

# Check anisotropy of K
if all(rates[l]['k'] > 0 for l in ['mB', 'x1', 'c']):
    k_vals = [rates[l]['k'] for l in ['mB', 'x1', 'c']]
    print(f"\n  k_mB : k_x1 : k_c = {k_vals[0]:.3f} : {k_vals[1]:.3f} : {k_vals[2]:.3f}")
    print(f"  Ratio k_x1/k_mB = {k_vals[1]/k_vals[0]:.2f}")
    print(f"  Ratio k_x1/k_c  = {k_vals[1]/k_vals[2]:.2f}")
    
    if k_vals[1] > k_vals[0] and k_vals[1] > k_vals[2]:
        print(f"\n  🔥 x1 (stretch/timing) contracts FASTEST")
        print(f"     → The 'geometric timing' channel has the highest contraction rate")
        print(f"     → This is the OPPOSITE of what 'locked channels immune' predicts")
        print(f"     → Or: stretch is NOT purely locked — it has diagnostic content too")

# ============================================================
# ATTRACTOR MANIFOLD: What is x*?
# ============================================================
print(f"\n" + "=" * 80)
print("ATTRACTOR MANIFOLD: Where does the population converge TO?")
print("=" * 80)

print(f"\n{'z':>6} {'⟨mB⟩':>8} {'⟨x1⟩':>8} {'⟨c⟩':>8}")
print("-" * 35)

mean_mB = []
mean_x1 = []
mean_c = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 15:
        continue
    
    mm = np.mean([s['mB'] for s in bin_sne])
    mx = np.mean([s['x1'] for s in bin_sne])
    mc = np.mean([s['c'] for s in bin_sne])
    
    mean_mB.append(mm)
    mean_x1.append(mx)
    mean_c.append(mc)
    
    print(f"{np.mean([s['zHD'] for s in bin_sne]):>6.2f} {mm:>8.3f} {mx:>8.3f} {mc:>8.4f}")

rho_mx, p_mx = spearmanr(z_arr[:len(mean_x1)], mean_x1)
rho_mc, p_mc = spearmanr(z_arr[:len(mean_c)], mean_c)

print(f"\n  ⟨x1⟩ vs z: ρ = {rho_mx:+.3f}, p = {p_mx:.4f}")
print(f"  ⟨c⟩ vs z:  ρ = {rho_mc:+.3f}, p = {p_mc:.4f}")
print(f"\n  Attractor x* ≈ (mB→distance-dependent, x1→{np.mean(mean_x1[-2:]):.2f}, c→{np.mean(mean_c[-2:]):.4f})")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_semigroup")
results_dir.mkdir(exist_ok=True)

results = {
    'rates': {k: {kk: float(vv) for kk, vv in v.items()} for k, v in rates.items()},
    'z_centers': [float(z) for z in z_arr],
    'intrinsic_vars': {
        'mB': [float(v) for v in vm],
        'x1': [float(v) for v in vx],
        'c': [float(v) for v in vc],
    }
}

with open(results_dir / "semigroup.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {results_dir}/")
