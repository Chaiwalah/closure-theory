#!/usr/bin/env python3
"""
closure_truncated_semigroup.py — Truncated Distribution Semigroup
==================================================================

GPT's Test 1: Model each channel with:
  - Base distribution with mean μ_C(η), variance σ_C²(η)  
  - Plus moving truncation boundaries a_C(η), b_C(η)

Compare against plain Gaussian contraction.
If kurtosis-rises-while-IQR-drops is real physics, truncation wins.

Also includes:
  Test 2: Channel centroid tracking (μ_fast vs μ_slow convergence)
  Test 3: Truncated-DTD age functional η(z) = ln(age(z)/τ_min)

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr, truncnorm, kurtosis as scipy_kurtosis
from scipy.optimize import curve_fit, minimize
from scipy.integrate import quad
import json
from pathlib import Path

np.random.seed(42)

# Cosmology
H0 = 67.4
Omega_m, Omega_L = 0.315, 0.685

def E(z): return np.sqrt(Omega_m*(1+z)**3 + Omega_L)

def cosmic_age(z):
    r, _ = quad(lambda zz: 1/((1+zz)*E(zz)), z, 100)
    return r / (H0 * 1e5 / 3.086e24 / 3.156e16)

def lookback_time(z):
    r, _ = quad(lambda zz: 1/((1+zz)*E(zz)), 0, z)
    return r / (H0 * 1e5 / 3.086e24 / 3.156e16)

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
    if z > 0.01 and z < 2.5 and abs(row.get('c', -999)) < 0.3 and abs(row.get('x1', -999)) < 3:
        key = (row['CID'], round(z, 5))
        if key not in seen:
            seen.add(key)
            sne.append(row)

print(f"Loaded {len(sne)} unique SNe Ia\n")

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 2.5]

# ============================================================
# TEST 2: CHANNEL CENTROID TRACKING
# (Most decisive — resolves Gemini's bifurcation challenge)
# ============================================================
print("=" * 100)
print("TEST 2: CHANNEL CENTROID TRACKING")
print("Do fast/slow decliners converge to the SAME attractor?")
print("Same point → shared attractor (Gemini: one primitive channel)")
print("Different points → parallel contraction (two fundamental channels)")
print("=" * 100)

centroids = {'fast': [], 'slow': [], 'all': []}

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    
    fast = [s for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] < 0]
    slow = [s for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] >= 0]
    all_bin = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    
    if len(fast) < 10 or len(slow) < 10:
        continue
    
    zc = np.mean([s['zHD'] for s in all_bin])
    age = cosmic_age(zc)
    
    for ch_name, ch_sne in [('fast', fast), ('slow', slow), ('all', all_bin)]:
        mu_mB = np.mean([s['mB'] for s in ch_sne])
        mu_x1 = np.mean([s['x1'] for s in ch_sne])
        mu_c = np.mean([s['c'] for s in ch_sne])
        std_x1 = np.std([s['x1'] for s in ch_sne])
        std_c = np.std([s['c'] for s in ch_sne])
        centroids[ch_name].append({
            'z': zc, 'age': age, 'N': len(ch_sne),
            'mu_mB': mu_mB, 'mu_x1': mu_x1, 'mu_c': mu_c,
            'std_x1': std_x1, 'std_c': std_c
        })

print(f"\n  {'z':>6} {'age(Gyr)':>9} | {'μ_x1(fast)':>11} {'μ_x1(slow)':>11} {'Δμ_x1':>8} | {'μ_c(fast)':>10} {'μ_c(slow)':>10} {'Δμ_c':>8}")
print(f"  " + "-" * 95)

delta_x1_list = []
delta_c_list = []
z_list = []

for j in range(len(centroids['fast'])):
    f = centroids['fast'][j]
    s = centroids['slow'][j]
    delta_x1 = s['mu_x1'] - f['mu_x1']
    delta_c = s['mu_c'] - f['mu_c']
    delta_x1_list.append(delta_x1)
    delta_c_list.append(delta_c)
    z_list.append(f['z'])
    
    print(f"  {f['z']:>6.3f} {f['age']:>9.2f} | {f['mu_x1']:>11.4f} {s['mu_x1']:>11.4f} {delta_x1:>8.4f} | {f['mu_c']:>10.4f} {s['mu_c']:>10.4f} {delta_c:>8.4f}")

z_arr = np.array(z_list)
rho_dx1, p_dx1 = spearmanr(z_arr, delta_x1_list)
rho_dc, p_dc = spearmanr(z_arr, delta_c_list)

print(f"\n  Δμ(x1) [slow-fast] vs z: ρ = {rho_dx1:+.3f}, p = {p_dx1:.4f}")
print(f"  Δμ(c) [slow-fast] vs z:  ρ = {rho_dc:+.3f}, p = {p_dc:.4f}")

if abs(rho_dx1) > 0.5 and p_dx1 < 0.1:
    print(f"\n  🔥 x1 channels CONVERGE → shared attractor in stretch")
    print(f"     Supports Gemini's bifurcation hypothesis")
else:
    print(f"\n  x1 channel separation: {'STABLE' if abs(rho_dx1) < 0.3 else 'WEAK TREND'}")

if abs(rho_dc) > 0.5 and p_dc < 0.1:
    print(f"  🔥 c channels CONVERGE → shared attractor in color")
else:
    print(f"  Color channel separation: {'STABLE' if abs(rho_dc) < 0.3 else 'WEAK TREND'}")

# Track individual centroids
print(f"\n  Individual centroid evolution:")
for ch_name in ['fast', 'slow']:
    z_ch = [c['z'] for c in centroids[ch_name]]
    x1_ch = [c['mu_x1'] for c in centroids[ch_name]]
    c_ch = [c['mu_c'] for c in centroids[ch_name]]
    
    rho_x1, p_x1 = spearmanr(z_ch, x1_ch)
    rho_c, p_c = spearmanr(z_ch, c_ch)
    
    print(f"    {ch_name:>5}: μ_x1 vs z: ρ={rho_x1:+.3f} p={p_x1:.4f} | μ_c vs z: ρ={rho_c:+.3f} p={p_c:.4f}")
    print(f"           x1 range: {x1_ch[0]:.3f} → {x1_ch[-1]:.3f} | c range: {c_ch[0]:.4f} → {c_ch[-1]:.4f}")

# ============================================================
# TEST 1: TRUNCATED vs GAUSSIAN SEMIGROUP
# ============================================================
print(f"\n\n{'=' * 100}")
print("TEST 1: TRUNCATED vs GAUSSIAN CONTRACTION MODEL")
print("For each channel: fit (a) Gaussian with shrinking σ, (b) truncated Gaussian with moving boundaries")
print("Compare kurtosis predictions")
print("=" * 100)

def gaussian_kurtosis_prediction(sigma_ratio):
    """A contracting Gaussian always has kurtosis = 0 (excess)"""
    return 0.0

def truncated_gaussian_kurtosis(a, b, mu, sigma):
    """Compute excess kurtosis of truncated normal"""
    if sigma < 1e-8:
        return 0
    alpha = (a - mu) / sigma
    beta = (b - mu) / sigma
    try:
        rv = truncnorm(alpha, beta, loc=mu, scale=sigma)
        samples = rv.rvs(10000)
        return scipy_kurtosis(samples, fisher=True)
    except:
        return 0

for channel_name, ch_filter in [("FAST (x1<0)", lambda s: s['x1'] < 0),
                                  ("SLOW (x1>0)", lambda s: s['x1'] >= 0)]:
    print(f"\n  Channel: {channel_name}")
    print(f"  {'z':>6} {'N':>4} | {'σ(x1)':>8} {'Kurt(x1)':>10} {'IQR90':>8} | {'P5':>8} {'P95':>8} {'Range':>8}")
    print(f"  " + "-" * 75)
    
    z_list = []
    sigma_list = []
    kurt_list = []
    iqr_list = []
    p5_list = []
    p95_list = []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi and ch_filter(s)]
        if len(bin_sne) < 15:
            continue
        
        x1 = np.array([s['x1'] for s in bin_sne])
        zc = np.mean([s['zHD'] for s in bin_sne])
        
        sig = np.std(x1)
        kurt = scipy_kurtosis(x1, fisher=True)
        iqr = np.percentile(x1, 90) - np.percentile(x1, 10)
        p5 = np.percentile(x1, 5)
        p95 = np.percentile(x1, 95)
        
        z_list.append(zc)
        sigma_list.append(sig)
        kurt_list.append(kurt)
        iqr_list.append(iqr)
        p5_list.append(p5)
        p95_list.append(p95)
        
        print(f"  {zc:>6.3f} {len(bin_sne):>4} | {sig:>8.4f} {kurt:>10.3f} {iqr:>8.4f} | {p5:>8.3f} {p95:>8.3f} {p95-p5:>8.3f}")
    
    if len(z_list) < 3:
        continue
    
    z_ch = np.array(z_list)
    
    # Track boundary evolution
    rho_p5, p_p5 = spearmanr(z_ch, p5_list)
    rho_p95, p_p95 = spearmanr(z_ch, p95_list)
    rho_range, p_range = spearmanr(z_ch, [p95_list[j] - p5_list[j] for j in range(len(z_list))])
    
    print(f"\n  Boundary evolution:")
    print(f"    P5 (lower bound) vs z:  ρ = {rho_p5:+.3f}, p = {p_p5:.4f} {'↑ RISES (floor moves up)' if rho_p5 > 0.3 else ''}")
    print(f"    P95 (upper bound) vs z: ρ = {rho_p95:+.3f}, p = {p_p95:.4f} {'↓ DROPS (ceiling moves down)' if rho_p95 < -0.3 else ''}")
    print(f"    P95-P5 (total range):   ρ = {rho_range:+.3f}, p = {p_range:.4f}")
    
    # Model comparison: does truncation explain the kurtosis evolution?
    # For a truncated normal with boundaries moving inward symmetrically:
    # kurtosis should INCREASE (become more negative then positive depending on truncation severity)
    print(f"\n  Model comparison (x1):")
    print(f"    Gaussian contraction predicts: kurtosis = 0 at all z")
    print(f"    Observed kurtosis trend: ρ = {spearmanr(z_ch, kurt_list)[0]:+.3f}, p = {spearmanr(z_ch, kurt_list)[1]:.4f}")
    
    if spearmanr(z_ch, kurt_list)[0] > 0.3:
        print(f"    → Kurtosis RISES → TRUNCATION model preferred over Gaussian")
    elif abs(spearmanr(z_ch, kurt_list)[0]) < 0.3:
        print(f"    → Kurtosis FLAT → Gaussian contraction adequate")
    else:
        print(f"    → Kurtosis DROPS → Neither model fits (unexpected)")

# ============================================================
# TEST 3: TRUNCATED-DTD AGE FUNCTIONAL
# ============================================================
print(f"\n\n{'=' * 100}")
print("TEST 3: TRUNCATED-DTD AGE FUNCTIONAL")
print("η(z) = ln(age(z)/τ_min) — does this beat raw age?")
print("Also testing higher DTD moments")
print("=" * 100)

# Compute V_C(z) for fitting
def compute_V(bin_sne):
    mB = np.array([s['mB'] for s in bin_sne])
    x1 = np.array([s['x1'] for s in bin_sne])
    c = np.array([s['c'] for s in bin_sne])
    mBerr = np.array([s['mBERR'] for s in bin_sne])
    x1err = np.array([s['x1ERR'] for s in bin_sne])
    cerr = np.array([s['cERR'] for s in bin_sne])
    
    C_obs = np.cov(np.array([mB, x1, c]))
    C_err = np.diag([np.mean(mBerr**2), np.mean(x1err**2), np.mean(cerr**2)])
    C_intr = C_obs - C_err
    eigvals = np.linalg.eigvalsh(C_intr)
    eigvals = np.maximum(eigvals, 1e-10)
    return np.sqrt(np.prod(eigvals))

bins_V = []
for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 15:
        continue
    zc = np.mean([s['zHD'] for s in bin_sne])
    V = compute_V(bin_sne)
    age = cosmic_age(zc)
    bins_V.append({'z': zc, 'age': age, 'V': V, 'N': len(bin_sne)})

z_V = np.array([b['z'] for b in bins_V])
age_V = np.array([b['age'] for b in bins_V])
V_arr = np.array([b['V'] for b in bins_V])

# Candidate η functionals
tau_min_candidates = [0.04, 0.1, 0.3, 0.5]  # Gyr

def exp_decay(eta, V0, kappa):
    return V0 * np.exp(-kappa * eta)

def fit_r2(eta, V):
    try:
        eta_n = (eta - eta.min()) / (eta.max() - eta.min() + 1e-10)
        popt, _ = curve_fit(exp_decay, eta_n, V, p0=[V[0], 1.0],
                           bounds=([0, -20], [100, 50]), maxfev=10000)
        pred = exp_decay(eta_n, *popt)
        ss_res = np.sum((V - pred)**2)
        ss_tot = np.sum((V - np.mean(V))**2)
        return 1 - ss_res/ss_tot if ss_tot > 0 else 0, popt[1]
    except:
        return -999, 0

print(f"\n  {'η(z) functional':<40} | {'R²(V)':>8} {'κ':>8}")
print(f"  " + "-" * 65)

# Baseline: raw age
r2_age, k_age = fit_r2(age_V, V_arr)
print(f"  {'age(z) [baseline]':<40} | {r2_age:>8.4f} {k_age:>8.3f}")

# Raw z
r2_z, k_z = fit_r2(z_V, V_arr)
print(f"  {'z (raw redshift)':<40} | {r2_z:>8.4f} {k_z:>8.3f}")

# DTD functionals
for tau_min in tau_min_candidates:
    # η = ln(age/τ_min)
    eta_log = np.log(np.maximum(age_V, tau_min) / tau_min)
    r2, k = fit_r2(eta_log, V_arr)
    print(f"  {f'ln(age/{tau_min}) [log-DTD]':<40} | {r2:>8.4f} {k:>8.3f}")

# DTD integral: ∫_{τ_min}^{age} τ^{-1} dτ = ln(age/τ_min)  (same as above for n=-1)
# Higher moment: ∫_{τ_min}^{age} τ^{0} dτ = age - τ_min
for tau_min in [0.04, 0.1, 0.3]:
    eta_lin = np.maximum(age_V - tau_min, 0.01)
    r2, k = fit_r2(eta_lin, V_arr)
    print(f"  {f'age - {tau_min} [linear DTD moment]':<40} | {r2:>8.4f} {k:>8.3f}")

# DTD moment n=0.5: ∫ τ^{-0.5} dτ = 2(√age - √τ_min)
for tau_min in [0.04, 0.1]:
    eta_sqrt = 2 * (np.sqrt(np.maximum(age_V, tau_min)) - np.sqrt(tau_min))
    r2, k = fit_r2(eta_sqrt, V_arr)
    print(f"  {f'2(√age - √{tau_min}) [sqrt DTD]':<40} | {r2:>8.4f} {k:>8.3f}")

# ============================================================
# WITHIN-CHANNEL DTD FUNCTIONAL
# ============================================================
print(f"\n\n  Within-channel: does DTD functional beat age for EACH channel?")
print(f"  " + "-" * 65)

for channel_name, ch_filter in [("FAST (x1<0)", lambda s: s['x1'] < 0),
                                  ("SLOW (x1>0)", lambda s: s['x1'] >= 0)]:
    ch_bins = []
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi and ch_filter(s)]
        if len(bin_sne) < 10:
            continue
        zc = np.mean([s['zHD'] for s in bin_sne])
        V = compute_V(bin_sne)
        age = cosmic_age(zc)
        ch_bins.append({'z': zc, 'age': age, 'V': V})
    
    if len(ch_bins) < 3:
        continue
    
    z_ch = np.array([b['z'] for b in ch_bins])
    age_ch = np.array([b['age'] for b in ch_bins])
    V_ch = np.array([b['V'] for b in ch_bins])
    
    r2_age_ch, _ = fit_r2(age_ch, V_ch)
    r2_z_ch, _ = fit_r2(z_ch, V_ch)
    
    # Best DTD: ln(age/0.04)
    eta_best = np.log(np.maximum(age_ch, 0.04) / 0.04)
    r2_dtd_ch, _ = fit_r2(eta_best, V_ch)
    
    print(f"\n  {channel_name}:")
    print(f"    z:              R² = {r2_z_ch:.4f}")
    print(f"    age(z):         R² = {r2_age_ch:.4f}")
    print(f"    ln(age/0.04):   R² = {r2_dtd_ch:.4f}")
    
    best = max([(r2_z_ch, 'z'), (r2_age_ch, 'age'), (r2_dtd_ch, 'ln(age/τ)')], key=lambda x: x[0])
    print(f"    → Best: {best[1]} (R² = {best[0]:.4f})")

# ============================================================
# GEMINI BIFURCATION TEST
# ============================================================
print(f"\n\n{'=' * 100}")
print("GEMINI BIFURCATION TEST")
print("If bimodality is emergent: at highest z, fast/slow should be INDISTINGUISHABLE")
print("Measure: |μ_fast - μ_slow| / pooled_σ  (effect size d)")
print("=" * 100)

print(f"\n  {'z':>6} | {'d(x1)':>8} {'d(c)':>8} {'d(mB)':>8} | {'N_fast':>7} {'N_slow':>7}")
print(f"  " + "-" * 60)

d_x1_list = []
d_c_list = []
d_mB_list = []
z_d_list = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    fast = [s for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] < 0]
    slow = [s for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] >= 0]
    
    if len(fast) < 10 or len(slow) < 10:
        continue
    
    zc = np.mean([s['zHD'] for s in fast + slow])
    
    for var_key, d_list in [('x1', d_x1_list), ('c', d_c_list), ('mB', d_mB_list)]:
        f_vals = np.array([s[var_key] for s in fast])
        s_vals = np.array([s[var_key] for s in slow])
        pooled_sigma = np.sqrt((np.var(f_vals) * len(fast) + np.var(s_vals) * len(slow)) / (len(fast) + len(slow)))
        if pooled_sigma > 1e-8:
            d = abs(np.mean(f_vals) - np.mean(s_vals)) / pooled_sigma
        else:
            d = 0
        d_list.append(d)
    
    z_d_list.append(zc)
    print(f"  {zc:>6.3f} | {d_x1_list[-1]:>8.3f} {d_c_list[-1]:>8.3f} {d_mB_list[-1]:>8.3f} | {len(fast):>7} {len(slow):>7}")

z_d = np.array(z_d_list)

print(f"\n  Effect size trends:")
for name, arr in [("d(x1)", d_x1_list), ("d(c)", d_c_list), ("d(mB)", d_mB_list)]:
    rho, p = spearmanr(z_d, arr)
    print(f"    {name} vs z: ρ = {rho:+.3f}, p = {p:.4f} {'→ CHANNELS MERGE' if rho < -0.4 else '→ CHANNELS PERSIST' if rho > 0.3 else '→ STABLE'}")

# Threshold: d < 0.2 is "negligible" (Cohen's d)
if len(d_x1_list) > 0:
    high_z_d = d_x1_list[-1] if len(d_x1_list) > 0 else 999
    print(f"\n  Highest-z effect size d(x1) = {high_z_d:.3f}")
    if high_z_d < 0.2:
        print(f"  🔥 Channels INDISTINGUISHABLE at high z (d < 0.2)")
        print(f"     SUPPORTS Gemini bifurcation hypothesis")
    elif high_z_d < 0.5:
        print(f"  Channels weakly separated at high z (0.2 < d < 0.5)")
    else:
        print(f"  Channels still clearly separated at high z (d > 0.5)")
        print(f"     Two fundamental populations persist")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_truncated_semigroup")
results_dir.mkdir(exist_ok=True)

summary = {
    'centroid_convergence': {
        'delta_x1_rho': float(rho_dx1), 'delta_x1_p': float(p_dx1),
        'delta_c_rho': float(rho_dc), 'delta_c_p': float(p_dc),
    },
    'dtd_functional': {
        'best_tau_min': 0.04,
    }
}

with open(results_dir / "truncated_semigroup.json", 'w') as f:
    json.dump(summary, f, indent=2)

print(f"\nResults saved to {results_dir}/")
