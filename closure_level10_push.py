#!/usr/bin/env python3
"""
closure_level10_push.py — Three Tests to Close Level 10
=========================================================

Test A: Knee test — P95 boundary vs age. Smooth exponential or hinge?
        Hinge = pathway conditioning time confirmed.

Test B: Survey-split test — Does P95 migration toward zero appear
        identically in every survey subsample?
        If yes, it's physics not cadence.

Test C: H₀ bias quantification — Refit Hubble diagram using only
        slow-channel SNe. Does H₀ shift?

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr, linregress
from scipy.optimize import curve_fit, minimize
from scipy.integrate import quad
import json
from pathlib import Path

np.random.seed(42)

# Cosmology
H0 = 67.4
Omega_m, Omega_L = 0.315, 0.685
c_km = 2.998e5  # km/s

def E(z): return np.sqrt(Omega_m*(1+z)**3 + Omega_L)

def cosmic_age(z):
    r, _ = quad(lambda zz: 1/((1+zz)*E(zz)), z, 100)
    return r / (H0 * 1e5 / 3.086e24 / 3.156e16)

def luminosity_distance(z):
    """Luminosity distance in Mpc for flat LCDM"""
    dc, _ = quad(lambda zz: c_km / (H0 * E(zz)), 0, z)
    return dc * (1 + z)

def distance_modulus(z):
    dL = luminosity_distance(z)
    if dL <= 0:
        return 0
    return 5 * np.log10(dL) + 25

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

# Survey ID mapping (IDSURVEY in Pantheon+)
survey_names = {
    1: 'SDSS', 4: 'SNLS', 5: 'CSP', 10: 'DES', 15: 'PS1',
    35: 'LOWZ', 36: 'LOWZ', 37: 'LOWZ', 38: 'LOWZ',
    40: 'LOWZ', 41: 'LOWZ', 42: 'LOWZ', 43: 'LOWZ',
    44: 'HST', 45: 'HST', 46: 'HST', 48: 'HST',
    50: 'SDSS', 51: 'SNLS', 56: 'PS1',
    61: 'DES', 62: 'DES', 63: 'DES', 64: 'DES', 65: 'DES',
    100: 'LOWZ', 101: 'CFA1', 104: 'CFA2', 106: 'CFA3',
    107: 'CFA3', 108: 'CFA4', 109: 'CFA4', 110: 'CSP',
    111: 'CSP', 112: 'CSP', 113: 'CSP',
    150: 'FOUND', 151: 'FOUND',
}

def get_survey(sn):
    sid = int(sn.get('IDSURVEY', 0))
    return survey_names.get(sid, f'S{sid}')

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 2.5]

# ============================================================
# TEST A: KNEE TEST — P95 boundary vs age
# ============================================================
print("=" * 100)
print("TEST A: KNEE TEST — P95(x1) boundary vs cosmic age")
print("Smooth exponential vs hinge model")
print("Hinge = pathway conditioning time (WD cooling / binary evolution gate)")
print("=" * 100)

# Get fast-channel P95 and P5 per z-bin
fast_boundaries = []
for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    fast = [s for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] < 0]
    if len(fast) < 15:
        continue
    
    x1 = np.array([s['x1'] for s in fast])
    zc = np.mean([s['zHD'] for s in fast])
    age = cosmic_age(zc)
    
    fast_boundaries.append({
        'z': zc, 'age': age, 'N': len(fast),
        'P5': np.percentile(x1, 5),
        'P10': np.percentile(x1, 10),
        'P90': np.percentile(x1, 90),
        'P95': np.percentile(x1, 95),
        'mean': np.mean(x1),
        'std': np.std(x1),
    })

age_arr = np.array([b['age'] for b in fast_boundaries])
z_arr_f = np.array([b['z'] for b in fast_boundaries])
p95_arr = np.array([b['P95'] for b in fast_boundaries])
p5_arr = np.array([b['P5'] for b in fast_boundaries])

# Normalize age for fitting (huge numbers from cosmic_age function)
age_norm = age_arr / age_arr[0]  # relative to lowest-z bin

print(f"\n  {'z':>6} {'age_rel':>8} {'N':>4} | {'P5':>8} {'P95':>8} {'mean':>8} {'σ':>8}")
print(f"  " + "-" * 60)
for b in fast_boundaries:
    print(f"  {b['z']:>6.3f} {b['age']/age_arr[0]:>8.4f} {b['N']:>4} | {b['P5']:>8.3f} {b['P95']:>8.3f} {b['mean']:>8.3f} {b['std']:>8.3f}")

# Model 1: Smooth exponential P95(age) = a + b*exp(-k*age)
def smooth_exp(age, a, b, k):
    return a + b * np.exp(-k * age)

# Model 2: Hinge/knee: P95 = a + b*max(0, age - age_knee)
def hinge(age, a, b, age_knee):
    return a + b * np.maximum(0, age - age_knee)

# Model 3: Linear
def linear(age, a, b):
    return a + b * age

print(f"\n  Fitting P95(x1) vs age_relative:")

# Fit smooth exponential
try:
    popt_exp, _ = curve_fit(smooth_exp, age_norm, p95_arr,
                             p0=[p95_arr[-1], p95_arr[0]-p95_arr[-1], 1.0],
                             bounds=([-3, -3, 0], [0, 3, 50]), maxfev=10000)
    pred_exp = smooth_exp(age_norm, *popt_exp)
    ss_res_exp = np.sum((p95_arr - pred_exp)**2)
    ss_tot = np.sum((p95_arr - np.mean(p95_arr))**2)
    r2_exp = 1 - ss_res_exp/ss_tot if ss_tot > 0 else 0
    print(f"    Exponential: R² = {r2_exp:.4f}, a={popt_exp[0]:.3f}, b={popt_exp[1]:.3f}, k={popt_exp[2]:.3f}")
except Exception as e:
    r2_exp = -1
    print(f"    Exponential: fit failed ({e})")

# Fit hinge
best_hinge_r2 = -999
best_hinge_params = None
for age_knee_trial in np.linspace(age_norm.min() + 0.05, age_norm.max() - 0.05, 20):
    try:
        popt_h, _ = curve_fit(lambda a, aa, bb: hinge(a, aa, bb, age_knee_trial),
                               age_norm, p95_arr,
                               p0=[p95_arr[-1], 0.1], maxfev=5000)
        pred_h = hinge(age_norm, popt_h[0], popt_h[1], age_knee_trial)
        ss_h = np.sum((p95_arr - pred_h)**2)
        r2_h = 1 - ss_h/ss_tot if ss_tot > 0 else 0
        if r2_h > best_hinge_r2:
            best_hinge_r2 = r2_h
            best_hinge_params = (popt_h[0], popt_h[1], age_knee_trial)
    except:
        pass

if best_hinge_params:
    # Convert knee age back to z
    knee_age_rel = best_hinge_params[2]
    knee_age_abs = knee_age_rel * age_arr[0]
    print(f"    Hinge: R² = {best_hinge_r2:.4f}, knee at age_rel={knee_age_rel:.3f}")
else:
    best_hinge_r2 = -1
    print(f"    Hinge: fit failed")

# Fit linear
try:
    popt_lin, _ = curve_fit(linear, age_norm, p95_arr, p0=[p95_arr[0], -0.1])
    pred_lin = linear(age_norm, *popt_lin)
    ss_lin = np.sum((p95_arr - pred_lin)**2)
    r2_lin = 1 - ss_lin/ss_tot if ss_tot > 0 else 0
    print(f"    Linear: R² = {r2_lin:.4f}, slope={popt_lin[1]:.4f}")
except:
    r2_lin = -1
    print(f"    Linear: fit failed")

print(f"\n  Model comparison for P95(age):")
models_a = [("Exponential", r2_exp), ("Hinge", best_hinge_r2), ("Linear", r2_lin)]
for name, r2 in sorted(models_a, key=lambda x: -x[1]):
    marker = " ← BEST" if r2 == max(m[1] for m in models_a) else ""
    print(f"    {name:<15} R² = {r2:.4f}{marker}")

if best_hinge_r2 > r2_exp + 0.02:
    print(f"\n  🔥 HINGE wins → pathway conditioning time CONFIRMED")
    print(f"     Knee at age_relative ≈ {best_hinge_params[2]:.3f}")
elif r2_exp > best_hinge_r2 + 0.02:
    print(f"\n  Smooth exponential wins → continuous truncation (no sharp gate)")
else:
    print(f"\n  Models indistinguishable with current binning")

# Same for P5
print(f"\n  P5 boundary comparison:")
try:
    popt_p5_lin, _ = curve_fit(linear, age_norm, p5_arr, p0=[p5_arr[0], 0.1])
    pred_p5_lin = linear(age_norm, *popt_p5_lin)
    ss_tot_p5 = np.sum((p5_arr - np.mean(p5_arr))**2)
    r2_p5_lin = 1 - np.sum((p5_arr - pred_p5_lin)**2)/ss_tot_p5 if ss_tot_p5 > 0 else 0
    print(f"    P5 linear: R² = {r2_p5_lin:.4f}, slope = {popt_p5_lin[1]:.4f}")
    print(f"    P95 linear slope / P5 linear slope = {popt_lin[1]/popt_p5_lin[1]:.2f}" if abs(popt_p5_lin[1]) > 1e-6 else "    P5 slope ~ 0")
except:
    pass

# ============================================================
# TEST B: SURVEY-SPLIT TEST
# ============================================================
print(f"\n\n{'=' * 100}")
print("TEST B: SURVEY-SPLIT — Does P95 migration appear in every survey?")
print("If yes → physics. If varies by survey → cadence/selection artifact.")
print("=" * 100)

# Group SNe by survey
survey_counts = {}
for s in sne:
    sv = get_survey(s)
    survey_counts[sv] = survey_counts.get(sv, 0) + 1

# Only test surveys with enough fast decliners
print(f"\n  Survey sizes (fast decliners only):")
fast_by_survey = {}
for s in sne:
    if s['x1'] < 0:
        sv = get_survey(s)
        if sv not in fast_by_survey:
            fast_by_survey[sv] = []
        fast_by_survey[sv].append(s)

for sv in sorted(fast_by_survey.keys(), key=lambda x: -len(fast_by_survey[x])):
    print(f"    {sv:>8}: {len(fast_by_survey[sv]):>5} fast decliners")

# For surveys with enough data, compute P95 trend
print(f"\n  P95(x1) trend by survey (fast decliners, z-binned):")
print(f"  {'Survey':<10} {'N_fast':>7} {'z_range':>12} | {'ρ(P95,z)':>10} {'p':>8} | {'Verdict':>15}")
print(f"  " + "-" * 75)

survey_z_edges = [0.01, 0.1, 0.2, 0.4, 1.0, 2.5]  # coarser bins for per-survey

consistent_count = 0
total_tested = 0

for sv in sorted(fast_by_survey.keys(), key=lambda x: -len(fast_by_survey[x])):
    sv_fast = fast_by_survey[sv]
    if len(sv_fast) < 30:
        continue
    
    z_sv = np.array([s['zHD'] for s in sv_fast])
    z_range = f"{z_sv.min():.2f}-{z_sv.max():.2f}"
    
    # Bin within survey
    z_list, p95_list = [], []
    for i in range(len(survey_z_edges) - 1):
        z_lo, z_hi = survey_z_edges[i], survey_z_edges[i+1]
        bin_sne = [s for s in sv_fast if z_lo <= s['zHD'] < z_hi]
        if len(bin_sne) < 8:
            continue
        x1 = np.array([s['x1'] for s in bin_sne])
        z_list.append(np.mean([s['zHD'] for s in bin_sne]))
        p95_list.append(np.percentile(x1, 95))
    
    if len(z_list) < 3:
        print(f"  {sv:<10} {len(sv_fast):>7} {z_range:>12} | {'<3 bins':>10} {'':>8} | {'SKIP':>15}")
        continue
    
    rho, p = spearmanr(z_list, p95_list)
    total_tested += 1
    
    if rho > 0.3:
        verdict = "RISES ✓"
        consistent_count += 1
    elif rho < -0.3:
        verdict = "DROPS ✗"
    else:
        verdict = "FLAT ?"
    
    print(f"  {sv:<10} {len(sv_fast):>7} {z_range:>12} | {rho:>10.3f} {p:>8.4f} | {verdict:>15}")

print(f"\n  Consistency: {consistent_count}/{total_tested} surveys show P95 rising with z")
if total_tested > 0:
    if consistent_count / total_tested >= 0.5:
        print(f"  ✓ Majority consistent → likely PHYSICS")
    else:
        print(f"  ✗ Mixed results → possible survey-dependent effects")

# ============================================================
# TEST C: H₀ BIAS QUANTIFICATION
# ============================================================
print(f"\n\n{'=' * 100}")
print("TEST C: H₀ BIAS — Slow-only vs All SNe Hubble diagram")
print("If fast channel introduces bias at high z, slow-only should give different H₀")
print("=" * 100)

# Standard Tripp estimator: μ_SN = mB - M + α·x1 - β·c
# We fit for M (absolute magnitude offset) which absorbs H₀

def tripp_residual(sn, M, alpha, beta):
    mu_sn = sn['mB'] - M + alpha * sn['x1'] - beta * sn['c']
    mu_cosmo = distance_modulus(sn['zHD'])
    return mu_sn - mu_cosmo

def fit_tripp(sample, label):
    """Fit Tripp parameters and return M, alpha, beta, scatter"""
    from scipy.optimize import minimize
    
    def chi2(params):
        M, alpha, beta = params
        resids = []
        for s in sample:
            mu_sn = s['mB'] - M + alpha * s['x1'] - beta * s['c']
            mu_cosmo = distance_modulus(s['zHD'])
            err = s.get('mBERR', 0.1)
            resids.append(((mu_sn - mu_cosmo) / max(err, 0.01))**2)
        return np.sum(resids)
    
    result = minimize(chi2, x0=[-19.3, 0.15, 3.0], method='Nelder-Mead',
                      options={'maxiter': 10000})
    M, alpha, beta = result.x
    
    # Compute residuals
    resids = []
    for s in sample:
        mu_sn = s['mB'] - M + alpha * s['x1'] - beta * s['c']
        mu_cosmo = distance_modulus(s['zHD'])
        resids.append(mu_sn - mu_cosmo)
    
    resids = np.array(resids)
    scatter = np.std(resids)
    
    return M, alpha, beta, scatter, resids

# Fit all SNe
M_all, alpha_all, beta_all, scatter_all, resids_all = fit_tripp(sne, "ALL")
print(f"\n  ALL SNe (N={len(sne)}):")
print(f"    M = {M_all:.4f}, α = {alpha_all:.4f}, β = {beta_all:.4f}")
print(f"    Scatter = {scatter_all:.4f} mag")

# Fit slow only
slow_sne = [s for s in sne if s['x1'] >= 0]
M_slow, alpha_slow, beta_slow, scatter_slow, resids_slow = fit_tripp(slow_sne, "SLOW")
print(f"\n  SLOW ONLY (N={len(slow_sne)}):")
print(f"    M = {M_slow:.4f}, α = {alpha_slow:.4f}, β = {beta_slow:.4f}")
print(f"    Scatter = {scatter_slow:.4f} mag")

# Fit fast only
fast_sne = [s for s in sne if s['x1'] < 0]
M_fast, alpha_fast, beta_fast, scatter_fast, resids_fast = fit_tripp(fast_sne, "FAST")
print(f"\n  FAST ONLY (N={len(fast_sne)}):")
print(f"    M = {M_fast:.4f}, α = {alpha_fast:.4f}, β = {beta_fast:.4f}")
print(f"    Scatter = {scatter_fast:.4f} mag")

# The key: M absorbs H₀ → ΔM implies ΔH₀
# μ = mB - M + corrections = 5log(dL) + 25
# dL ∝ c/H₀ → M = M_true + 5log(H₀/70) + const
# So ΔM = 5 * Δlog(H₀) → ΔH₀/H₀ = (10^(ΔM/5) - 1)
delta_M = M_slow - M_all
delta_M_fast = M_fast - M_all
H0_ratio = 10**(delta_M / 5)
H0_implied = H0 * H0_ratio

print(f"\n  ΔM (slow - all) = {delta_M:+.4f} mag")
print(f"  ΔM (fast - all) = {delta_M_fast:+.4f} mag")
print(f"  Implied H₀ shift: {H0:.1f} → {H0_implied:.1f} km/s/Mpc ({(H0_ratio-1)*100:+.2f}%)")

if abs(delta_M) > 0.02:
    print(f"\n  🔥 SIGNIFICANT SHIFT — fast channel introduces systematic bias")
    print(f"     Slow-only H₀ differs by {abs(delta_M):.3f} mag ({abs((H0_ratio-1)*100):.1f}%)")
else:
    print(f"\n  Shift is small ({abs(delta_M):.4f} mag) — bias is minimal")

# Z-dependent bias check
print(f"\n  Z-dependent residual comparison:")
print(f"  {'z':>6} | {'⟨resid⟩_all':>12} {'⟨resid⟩_slow':>13} {'⟨resid⟩_fast':>13} | {'Δ(slow-all)':>12}")
print(f"  " + "-" * 75)

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    
    # All
    idx_all = [j for j, s in enumerate(sne) if z_lo <= s['zHD'] < z_hi]
    if len(idx_all) < 10:
        continue
    zc = np.mean([sne[j]['zHD'] for j in idx_all])
    mean_all = np.mean([resids_all[j] for j in idx_all])
    
    # Slow
    slow_bin = [s for s in slow_sne if z_lo <= s['zHD'] < z_hi]
    if len(slow_bin) < 5:
        mean_slow_r = np.nan
    else:
        slow_resids = []
        for s in slow_bin:
            mu_sn = s['mB'] - M_slow + alpha_slow * s['x1'] - beta_slow * s['c']
            mu_cosmo = distance_modulus(s['zHD'])
            slow_resids.append(mu_sn - mu_cosmo)
        mean_slow_r = np.mean(slow_resids)
    
    # Fast  
    fast_bin = [s for s in fast_sne if z_lo <= s['zHD'] < z_hi]
    if len(fast_bin) < 5:
        mean_fast_r = np.nan
    else:
        fast_resids = []
        for s in fast_bin:
            mu_sn = s['mB'] - M_fast + alpha_fast * s['x1'] - beta_fast * s['c']
            mu_cosmo = distance_modulus(s['zHD'])
            fast_resids.append(mu_sn - mu_cosmo)
        mean_fast_r = np.mean(fast_resids)
    
    delta = mean_slow_r - mean_all if not np.isnan(mean_slow_r) else np.nan
    
    print(f"  {zc:>6.3f} | {mean_all:>12.4f} {mean_slow_r:>13.4f} {mean_fast_r:>13.4f} | {delta:>12.4f}")

# Check if residual difference grows with z
all_z_resids = []
for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    
    fast_bin = [s for s in fast_sne if z_lo <= s['zHD'] < z_hi]
    slow_bin = [s for s in slow_sne if z_lo <= s['zHD'] < z_hi]
    
    if len(fast_bin) < 5 or len(slow_bin) < 5:
        continue
    
    zc = np.mean([s['zHD'] for s in fast_bin + slow_bin])
    
    fast_mu = [s['mB'] - M_all + alpha_all * s['x1'] - beta_all * s['c'] for s in fast_bin]
    slow_mu = [s['mB'] - M_all + alpha_all * s['x1'] - beta_all * s['c'] for s in slow_bin]
    fast_cosmo = [distance_modulus(s['zHD']) for s in fast_bin]
    slow_cosmo = [distance_modulus(s['zHD']) for s in slow_bin]
    
    fast_resid = np.mean(np.array(fast_mu) - np.array(fast_cosmo))
    slow_resid = np.mean(np.array(slow_mu) - np.array(slow_cosmo))
    
    all_z_resids.append({'z': zc, 'fast_resid': fast_resid, 'slow_resid': slow_resid,
                         'gap': slow_resid - fast_resid})

if len(all_z_resids) >= 3:
    z_r = np.array([r['z'] for r in all_z_resids])
    gap_r = np.array([r['gap'] for r in all_z_resids])
    rho_gap, p_gap = spearmanr(z_r, gap_r)
    
    print(f"\n  Residual gap (slow-fast) vs z: ρ = {rho_gap:+.3f}, p = {p_gap:.4f}")
    if abs(rho_gap) > 0.5 and p_gap < 0.1:
        print(f"  🔥 GAP EVOLVES WITH Z — fast channel creates z-dependent bias!")
        print(f"     This COULD contribute to H₀ tension if not accounted for")
    else:
        print(f"  Gap is stable with z — no z-dependent bias from channel mixing")

# ============================================================
# HIGH-Z vs LOW-Z H₀ SPLIT
# ============================================================
print(f"\n\n  Split-z H₀ test (low z < 0.15 vs high z > 0.15):")

for label, sample in [("ALL", sne), ("SLOW", slow_sne), ("FAST", fast_sne)]:
    low_z = [s for s in sample if s['zHD'] < 0.15]
    high_z = [s for s in sample if s['zHD'] >= 0.15]
    
    if len(low_z) < 20 or len(high_z) < 20:
        continue
    
    M_lo, _, _, scat_lo, _ = fit_tripp(low_z, f"{label}_low")
    M_hi, _, _, scat_hi, _ = fit_tripp(high_z, f"{label}_high")
    
    delta = M_hi - M_lo
    h0_shift = (10**(delta/5) - 1) * 100
    
    print(f"    {label:>5}: M_low={M_lo:.3f} M_high={M_hi:.3f} ΔM={delta:+.4f} → H₀ shift: {h0_shift:+.2f}%")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_level10")
results_dir.mkdir(exist_ok=True)

summary = {
    'test_a_knee': {
        'r2_exponential': float(r2_exp),
        'r2_hinge': float(best_hinge_r2),
        'r2_linear': float(r2_lin),
    },
    'test_c_h0': {
        'M_all': float(M_all), 'M_slow': float(M_slow), 'M_fast': float(M_fast),
        'delta_M': float(delta_M),
        'H0_implied': float(H0_implied),
    }
}

with open(results_dir / "level10.json", 'w') as f:
    json.dump(summary, f, indent=2)

print(f"\n\nResults saved to {results_dir}/")
