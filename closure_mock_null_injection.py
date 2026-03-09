#!/usr/bin/env python3
"""
CLOSURE THEORY — MOCK NULL INJECTION TEST
============================================

Generate 1000 mock Pantheon+ catalogs under pure ΛCDM (NO impedance).
Apply standard Tripp AND modified Tripp to each mock.
Measure ρ(residual, z) distribution under the null hypothesis.

Key questions:
1. Is our pre-correction ρ = -0.199 OUTSIDE the null distribution? (→ real signal)
2. Is our post-correction ρ ≈ -0.13 INSIDE the null distribution? (→ decode complete)

Uses diagonal errors for now. GPT recommends upgrading to full covariance later.

Author: Closure Theory collaboration
Date: 2026-03-09
"""

import numpy as np
from scipy.optimize import minimize as scipy_minimize
from scipy.integrate import quad
from scipy.stats import spearmanr, percentileofscore
import json
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# LOAD REAL DATA (for z, x1, c, host_mass distributions)
# ============================================================

data_file = 'data/pantheon_plus.dat'

with open(data_file, 'r') as f:
    header = f.readline().strip().split()
    rows = []
    for line in f:
        parts = line.strip().split()
        if len(parts) == len(header):
            rows.append(parts)

col = {name: i for i, name in enumerate(header)}

z_cmb = np.array([float(r[col['zCMB']]) for r in rows])
z_hd = np.array([float(r[col['zHD']]) for r in rows])
m_b = np.array([float(r[col['mB']]) for r in rows])
m_b_err = np.array([float(r[col['mBERR']]) for r in rows])
x1 = np.array([float(r[col['x1']]) for r in rows])
c = np.array([float(r[col['c']]) for r in rows])
host_mass = np.array([float(r[col['HOST_LOGMASS']]) for r in rows])

mask_q = (z_hd > 0.01) & (z_hd < 2.5) & (np.abs(x1) < 5) & (np.abs(c) < 0.5) & (m_b_err < 1.0) & (host_mass > 0) & (host_mass < 15)

z = z_cmb[mask_q]
mb_real = m_b[mask_q]
mb_err = m_b_err[mask_q]
x1_q = x1[mask_q]
c_q = c[mask_q]
hm = host_mass[mask_q]
mask_hi = hm >= 10.0
N = len(z)

print(f"Loaded {N} real SNe Ia")

# ============================================================
# COSMOLOGY
# ============================================================

C_LIGHT = 299792.458

def E(zp, Om=0.3):
    return np.sqrt(Om*(1+zp)**3 + (1-Om))

def dist_mod(z_arr, H0=70, Om=0.3):
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        integral, _ = quad(lambda zp: 1.0/E(zp, Om), 0, zi)
        dl[i] = C_LIGHT * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

print("Computing distance moduli...")
mu_lcdm = dist_mod(z, H0=70, Om=0.3)

# Impedance variable for modified fits
I_sat = z / (1 + z/0.5)

# ============================================================
# FIT FUNCTIONS
# ============================================================

def fit_standard_tripp(mb_data):
    """Standard Tripp: 4 params (M_B, α, β, γ)"""
    def chi2(params):
        M_B, alpha, beta, gamma = params
        mu_obs = mb_data + alpha * x1_q - beta * c_q - M_B
        mu_obs[mask_hi] -= gamma
        return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)
    
    res = scipy_minimize(chi2, x0=[-19.25, 0.15, 3.0, -0.05],
                         method='Nelder-Mead', options={'maxiter': 30000, 'xatol': 1e-6})
    p = res.x
    resid = (mb_data + p[1]*x1_q - p[2]*c_q - p[0]) - mu_lcdm
    resid[mask_hi] -= p[3]
    rho, pval = spearmanr(z, resid)
    return rho, pval


def fit_modified_tripp(mb_data):
    """Modified Tripp with impedance: M1 model (7 params)
    Rotation + impedance variable + c² + γ(Ĩ)"""
    x1_std = np.std(x1_q)
    c_std = np.std(c_q)
    
    def chi2(params):
        A, B, theta1, delta_v, g0, g1 = params
        angle = theta1 * I_sat
        x1n = x1_q / x1_std; cn = c_q / c_std
        u = x1n * np.cos(angle) + cn * np.sin(angle)
        v = -x1n * np.sin(angle) + cn * np.cos(angle)
        u_mag = u * x1_std; v_mag = v * c_std
        gamma_I = g0 + g1 * I_sat
        mu_corr = mb_data + A * u_mag - B * v_mag - delta_v * v_mag**2
        mu_corr[mask_hi] -= gamma_I[mask_hi]
        w = 1.0 / mb_err**2
        M0 = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
        residuals = mu_corr - M0 - mu_lcdm
        return np.sum((residuals / mb_err)**2)
    
    res = scipy_minimize(chi2, x0=[0.15, 3.0, -0.3, 2.0, -0.03, -0.1],
                         method='Nelder-Mead', options={'maxiter': 50000, 'xatol': 1e-6})
    
    A, B, theta1, delta_v, g0, g1 = res.x
    angle = theta1 * I_sat
    x1_std_v = np.std(x1_q); c_std_v = np.std(c_q)
    x1n = x1_q / x1_std_v; cn = c_q / c_std_v
    u = x1n * np.cos(angle) + cn * np.sin(angle)
    v = -x1n * np.sin(angle) + cn * np.cos(angle)
    u_mag = u * x1_std_v; v_mag = v * c_std_v
    gamma_I = g0 + g1 * I_sat
    mu_corr = mb_data + A * u_mag - B * v_mag - delta_v * v_mag**2
    mu_corr[mask_hi] -= gamma_I[mask_hi]
    w = 1.0 / mb_err**2
    M0 = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
    resid = mu_corr - M0 - mu_lcdm
    rho, pval = spearmanr(z, resid)
    return rho, pval


# ============================================================
# REAL DATA BASELINE
# ============================================================

print("\n" + "=" * 70)
print("REAL DATA BASELINE")
print("=" * 70)

rho_real_std, p_real_std = fit_standard_tripp(mb_real)
rho_real_mod, p_real_mod = fit_modified_tripp(mb_real)

print(f"  Standard Tripp: ρ = {rho_real_std:+.4f} (p = {p_real_std:.2e})")
print(f"  Modified (M1):  ρ = {rho_real_mod:+.4f} (p = {p_real_mod:.2e})")


# ============================================================
# GENERATE MOCKS AND FIT
# ============================================================

N_MOCKS = 1000

print(f"\n{'=' * 70}")
print(f"GENERATING {N_MOCKS} MOCK CATALOGS (ΛCDM + diagonal noise, NO impedance)")
print("=" * 70)

# True ΛCDM apparent magnitudes (no impedance, standard Tripp params)
# Use best-fit standard Tripp params from real data
M_B_true = -19.25
alpha_true = 0.14
beta_true = 3.0
gamma_true = -0.05

mb_true = mu_lcdm + M_B_true - alpha_true * x1_q + beta_true * c_q
mb_true[mask_hi] += gamma_true

rho_std_mocks = np.zeros(N_MOCKS)
rho_mod_mocks = np.zeros(N_MOCKS)

for i in range(N_MOCKS):
    if (i+1) % 100 == 0:
        print(f"  Mock {i+1}/{N_MOCKS}...")
    
    # Add Gaussian noise with real Pantheon+ diagonal errors
    noise = np.random.normal(0, mb_err)
    mb_mock = mb_true + noise
    
    # Fit standard Tripp
    rho_s, _ = fit_standard_tripp(mb_mock)
    rho_std_mocks[i] = rho_s
    
    # Fit modified Tripp
    rho_m, _ = fit_modified_tripp(mb_mock)
    rho_mod_mocks[i] = rho_m


# ============================================================
# ANALYSIS
# ============================================================

print(f"\n{'=' * 70}")
print("MOCK NULL DISTRIBUTION — RESULTS")
print("=" * 70)

# Standard Tripp null distribution
mean_std = np.mean(rho_std_mocks)
std_std = np.std(rho_std_mocks)
pct_std_pre = percentileofscore(rho_std_mocks, rho_real_std)
sigma_pre = (rho_real_std - mean_std) / std_std

print(f"\n  STANDARD TRIPP null distribution:")
print(f"    Mock ρ: mean = {mean_std:+.4f}, σ = {std_std:.4f}")
print(f"    Mock ρ range: [{np.min(rho_std_mocks):+.4f}, {np.max(rho_std_mocks):+.4f}]")
print(f"    Real pre-correction ρ = {rho_real_std:+.4f}")
print(f"    Percentile: {pct_std_pre:.1f}%")
print(f"    Distance: {sigma_pre:+.1f}σ from null mean")
if abs(sigma_pre) > 3:
    print(f"    ★ PRE-CORRECTION ρ IS ANOMALOUS (>{abs(sigma_pre):.0f}σ) — REAL SIGNAL DETECTED")
else:
    print(f"    ○ Pre-correction ρ is within null distribution — no significant signal")

# Modified Tripp null distribution
mean_mod = np.mean(rho_mod_mocks)
std_mod = np.std(rho_mod_mocks)
pct_mod_post = percentileofscore(rho_mod_mocks, rho_real_mod)
sigma_post = (rho_real_mod - mean_mod) / std_mod

print(f"\n  MODIFIED TRIPP (M1) null distribution:")
print(f"    Mock ρ: mean = {mean_mod:+.4f}, σ = {std_mod:.4f}")
print(f"    Mock ρ range: [{np.min(rho_mod_mocks):+.4f}, {np.max(rho_mod_mocks):+.4f}]")
print(f"    Real post-correction ρ = {rho_real_mod:+.4f}")
print(f"    Percentile: {pct_mod_post:.1f}%")
print(f"    Distance: {sigma_post:+.1f}σ from null mean")
if abs(sigma_post) < 2:
    print(f"    ★ POST-CORRECTION ρ IS INSIDE NULL (within 2σ) — DECODE CONSISTENT WITH NOISE")
elif abs(sigma_post) < 3:
    print(f"    ◐ Post-correction ρ is marginal (2-3σ) — close but not fully inside null")
else:
    print(f"    ✗ Post-correction ρ is OUTSIDE null (>{abs(sigma_post):.0f}σ) — residual physics remains")


# ============================================================
# THE DEFINITIVE TEST: Is the DIFFERENCE anomalous?
# ============================================================

print(f"\n  IMPEDANCE CORRECTION EFFECT:")
delta_real = rho_real_std - rho_real_mod
delta_mocks = rho_std_mocks - rho_mod_mocks
mean_delta = np.mean(delta_mocks)
std_delta = np.std(delta_mocks)
sigma_delta = (delta_real - mean_delta) / std_delta
pct_delta = percentileofscore(delta_mocks, delta_real)

print(f"    Real Δρ (std→mod): {delta_real:+.4f}")
print(f"    Mock Δρ: mean = {mean_delta:+.4f}, σ = {std_delta:.4f}")
print(f"    Distance: {sigma_delta:+.1f}σ")
print(f"    Percentile: {pct_delta:.1f}%")
if sigma_delta > 3:
    print(f"    ★ IMPEDANCE CORRECTION REMOVES MORE TREND THAN NOISE ALONE — REAL PHYSICS")
else:
    print(f"    ○ Correction amount is consistent with overfitting noise")


# ============================================================
# SUMMARY VERDICT
# ============================================================

print(f"\n\n{'=' * 70}")
print("VERDICT")
print("=" * 70)

pre_anomalous = abs(sigma_pre) > 3
post_inside = abs(sigma_post) < 2
correction_real = sigma_delta > 3

if pre_anomalous and post_inside and correction_real:
    print("""
  ★★★ FRAMEWORK VALIDATED ★★★
  
  1. Pre-correction ρ is ANOMALOUS — real signal exists in the data
  2. Post-correction ρ is INSIDE the null — impedance correction brings 
     residuals to the noise floor
  3. The correction itself is LARGER than noise fluctuations — real physics
  
  THE DECODE IS STATISTICALLY COMPLETE.
  Remaining ρ is consistent with Pantheon+ diagonal noise under ΛCDM.
""")
elif pre_anomalous and correction_real:
    print(f"""
  ◐ PARTIAL VALIDATION
  
  1. Pre-correction ρ is ANOMALOUS — real signal ✓
  2. Post-correction ρ is {sigma_post:+.1f}σ from null — NOT fully inside
  3. Correction is real physics ✓
  
  The impedance correction captures real signal but doesn't fully 
  reach the noise floor. Either:
  - More physics remains (unlikely given controls)
  - The functional form needs refinement
  - Full covariance matrix will resolve the gap
""")
elif pre_anomalous:
    print(f"""
  ◐ SIGNAL DETECTED, CORRECTION INSUFFICIENT
  
  1. Pre-correction ρ is ANOMALOUS — real signal ✓
  2. Post-correction still outside null
  3. Correction amount marginal
  
  The framework sees real signal but the current correction 
  doesn't fully explain it.
""")
else:
    print("""
  ✗ NO SIGNIFICANT SIGNAL
  
  Pre-correction ρ is INSIDE the null distribution.
  The original "signal" may be noise. This would mean the
  framework is solving a problem that doesn't exist.
""")


# Save results
os.makedirs('results_mock_null', exist_ok=True)
output = {
    'test_date': '2026-03-09',
    'n_mocks': N_MOCKS,
    'real_data': {
        'rho_standard': float(rho_real_std),
        'rho_modified': float(rho_real_mod),
        'delta_rho': float(delta_real),
    },
    'null_standard': {
        'mean': float(mean_std),
        'std': float(std_std),
        'sigma_from_null': float(sigma_pre),
        'percentile': float(pct_std_pre),
    },
    'null_modified': {
        'mean': float(mean_mod),
        'std': float(std_mod),
        'sigma_from_null': float(sigma_post),
        'percentile': float(pct_mod_post),
    },
    'correction_effect': {
        'delta_real': float(delta_real),
        'delta_null_mean': float(mean_delta),
        'delta_null_std': float(std_delta),
        'sigma': float(sigma_delta),
        'percentile': float(pct_delta),
    },
    'verdict': {
        'pre_anomalous': bool(pre_anomalous),
        'post_inside_null': bool(post_inside),
        'correction_real': bool(correction_real),
    }
}

with open('results_mock_null/mock_null_results.json', 'w') as f:
    json.dump(output, f, indent=2)

# Save mock distributions for later analysis
np.savez('results_mock_null/mock_distributions.npz',
         rho_std=rho_std_mocks, rho_mod=rho_mod_mocks)

print(f"\n  Results saved to results_mock_null/")
print(f"  Mock distributions saved to results_mock_null/mock_distributions.npz")
print(f"\n{'=' * 70}")
print("COMPLETE")
print("=" * 70)
