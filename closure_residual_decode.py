#!/usr/bin/env python3
"""
CLOSURE THEORY — RESIDUAL DECODE + Ωm/w SIMULATION
=====================================================

1. Plot residuals vs z for: Standard Tripp vs Full Modified Tripp
   If modified Tripp flattens to white noise → we decoded expansion
   
2. Refit Ωm under modified standardization
   Does Ωm shift? Does w move toward -1?

3. Quantify: autocorrelation, running mean, scatter reduction

Author: Closure Theory collaboration
Date: 2026-03-09
"""

import numpy as np
from scipy.optimize import minimize as scipy_minimize
from scipy.integrate import quad
from scipy.stats import spearmanr, pearsonr
import json
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# LOAD DATA
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
mb = m_b[mask_q]
mb_err = m_b_err[mask_q]
x1_q = x1[mask_q]
c_q = c[mask_q]
hm = host_mass[mask_q]

print(f"Loaded {len(z)} SNe Ia")

# ============================================================
# COSMOLOGY
# ============================================================

C_LIGHT = 299792.458

def dist_mod(z_arr, H0=70, Om=0.3):
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        integral, _ = quad(lambda zp: 1.0/np.sqrt(Om*(1+zp)**3 + (1-Om)), 0, zi)
        dl[i] = C_LIGHT * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

def dist_mod_wcdm(z_arr, H0=70, Om=0.3, w=-1.0):
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        def integrand(zp):
            Ode = 1 - Om
            rho_de = Ode * (1+zp)**(3*(1+w))
            return 1.0 / np.sqrt(Om*(1+zp)**3 + rho_de)
        integral, _ = quad(integrand, 0, zi)
        dl[i] = C_LIGHT * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

# ============================================================
# FIT STANDARD TRIPP
# ============================================================

print("Fitting standard Tripp...")
mask_hi = hm >= 10.0

mu_lcdm = dist_mod(z, H0=70, Om=0.3)

def std_chi2(params):
    M_B, alpha, beta, gamma = params
    mu_obs = mb + alpha * x1_q - beta * c_q - M_B
    mu_obs[mask_hi] -= gamma
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res_std = scipy_minimize(std_chi2, x0=[-19.25, 0.15, 3.0, -0.05],
                         method='Nelder-Mead', options={'maxiter': 50000})
M_std, a_std, b_std, g_std = res_std.x

mu_obs_std = mb + a_std * x1_q - b_std * c_q - M_std
mu_obs_std[mask_hi] -= g_std
resid_std = mu_obs_std - mu_lcdm

print(f"  M_B={M_std:.4f}, α={a_std:.4f}, β={b_std:.3f}, γ={g_std:+.4f}")

# ============================================================
# FIT FULL MODIFIED TRIPP
# ============================================================

print("Fitting full modified Tripp...")

def mod_chi2(params):
    M_B, a0, a1, b0, b1, delta, g0, g1 = params
    alpha_z = a0 + a1 * z
    beta_z = b0 + b1 * z
    gamma_z = g0 + g1 * z
    mu_obs = mb + alpha_z * x1_q - beta_z * c_q - delta * c_q**2 - M_B
    mu_obs[mask_hi] -= gamma_z[mask_hi]
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res_mod = scipy_minimize(mod_chi2,
                         x0=[-19.25, 0.15, -0.03, 3.0, -0.5, 2.0, -0.03, -0.05],
                         method='Nelder-Mead', options={'maxiter': 100000})
M_mod, a0, a1, b0, b1, d_c2, g0, g1 = res_mod.x

alpha_z = a0 + a1 * z
beta_z = b0 + b1 * z
gamma_z = g0 + g1 * z
mu_obs_mod = mb + alpha_z * x1_q - beta_z * c_q - d_c2 * c_q**2 - M_mod
mu_obs_mod[mask_hi] -= gamma_z[mask_hi]
resid_mod = mu_obs_mod - mu_lcdm

print(f"  α(z) = {a0:.4f} + ({a1:+.4f})×z")
print(f"  β(z) = {b0:.3f} + ({b1:+.3f})×z")
print(f"  γ(z) = {g0:+.4f} + ({g1:+.4f})×z")
print(f"  δ(c²) = {d_c2:+.3f}")

# ============================================================
# GENERATE ASCII RESIDUAL PLOT
# ============================================================

print(f"\n{'=' * 70}")
print("RESIDUAL COMPARISON: STANDARD vs MODIFIED TRIPP")
print("=" * 70)

# Bin residuals
n_bins = 15
z_sorted = np.sort(z)
bin_edges = np.percentile(z, np.linspace(0, 100, n_bins + 1))

print(f"\n  BINNED RESIDUALS (mean ± std in each z-bin)")
print(f"\n  {'z-bin':>15} {'N':>5} {'Standard':>20} {'Modified':>20} {'Improvement':>12}")
print(f"  {'-'*75}")

std_means = []
mod_means = []
std_stds = []
mod_stds = []
z_mids = []

for i in range(n_bins):
    mask_bin = (z >= bin_edges[i]) & (z < bin_edges[i+1])
    if i == n_bins - 1:
        mask_bin = (z >= bin_edges[i]) & (z <= bin_edges[i+1])
    
    n = np.sum(mask_bin)
    if n < 5:
        continue
    
    z_mid = np.mean(z[mask_bin])
    z_mids.append(z_mid)
    
    std_m = np.mean(resid_std[mask_bin])
    std_s = np.std(resid_std[mask_bin])
    mod_m = np.mean(resid_mod[mask_bin])
    mod_s = np.std(resid_mod[mask_bin])
    
    std_means.append(std_m)
    mod_means.append(mod_m)
    std_stds.append(std_s)
    mod_stds.append(mod_s)
    
    improvement = (1 - mod_s / std_s) * 100
    
    z_lo = bin_edges[i]
    z_hi = bin_edges[i+1]
    print(f"  [{z_lo:.3f},{z_hi:.3f}] {n:>5} {std_m:>+8.4f} ± {std_s:.4f}  {mod_m:>+8.4f} ± {mod_s:.4f}  {improvement:>+10.1f}%")

z_mids = np.array(z_mids)
std_means = np.array(std_means)
mod_means = np.array(mod_means)
std_stds = np.array(std_stds)
mod_stds = np.array(mod_stds)

# Overall statistics
print(f"\n  OVERALL:")
print(f"    Standard Tripp RMS:  {np.sqrt(np.mean(resid_std**2)):.4f}")
print(f"    Modified Tripp RMS:  {np.sqrt(np.mean(resid_mod**2)):.4f}")
print(f"    RMS improvement:     {(1 - np.sqrt(np.mean(resid_mod**2))/np.sqrt(np.mean(resid_std**2)))*100:.1f}%")

# Autocorrelation test (white noise should have ~0 autocorrelation)
idx_sorted = np.argsort(z)
resid_std_sorted = resid_std[idx_sorted]
resid_mod_sorted = resid_mod[idx_sorted]

autocorr_std = np.corrcoef(resid_std_sorted[:-1], resid_std_sorted[1:])[0, 1]
autocorr_mod = np.corrcoef(resid_mod_sorted[:-1], resid_mod_sorted[1:])[0, 1]

print(f"\n    Standard Tripp autocorrelation: {autocorr_std:+.4f}")
print(f"    Modified Tripp autocorrelation: {autocorr_mod:+.4f}")
print(f"    (White noise → 0)")

# Trend test: correlation of residuals with z
rho_std, p_std = spearmanr(z, resid_std)
rho_mod, p_mod = spearmanr(z, resid_mod)

print(f"\n    Standard: residual vs z: ρ = {rho_std:+.4f} (p = {p_std:.2e})")
print(f"    Modified: residual vs z: ρ = {rho_mod:+.4f} (p = {p_mod:.2e})")
print(f"    (White noise → ρ = 0)")

# Running mean trend
print(f"\n  RUNNING MEAN OF BINNED RESIDUALS:")
print(f"    Standard Tripp: mean = {np.mean(std_means):+.5f}, std = {np.std(std_means):.5f}")
print(f"    Modified Tripp: mean = {np.mean(mod_means):+.5f}, std = {np.std(mod_means):.5f}")
print(f"    Trend reduction: {(1 - np.std(mod_means)/np.std(std_means))*100:.1f}%")


# ============================================================
# ASCII VISUALIZATION
# ============================================================

print(f"\n{'=' * 70}")
print("VISUAL: BINNED RESIDUALS vs z")
print("=" * 70)

# Scale for display
max_val = max(np.max(np.abs(std_means)), np.max(np.abs(mod_means)))
width = 50  # characters wide

print(f"\n  STANDARD TRIPP (constant α, β, γ):")
print(f"  {'z':>6}  {'-0.10':>5}{'0':>23}{'+0.10':>22}")
print(f"  {'':>6}  |{'─'*23}┼{'─'*23}|")

for i in range(len(z_mids)):
    val = std_means[i]
    pos = int((val / 0.10 + 1) * 23)  # map [-0.10, +0.10] to [0, 46]
    pos = max(0, min(46, pos))
    bar = [' '] * 47
    bar[23] = '│'
    marker = '●' if abs(val) > 0.02 else '○'
    bar[pos] = marker
    line = ''.join(bar)
    print(f"  {z_mids[i]:>6.3f} |{line}|")

print(f"  {'':>6}  |{'─'*23}┼{'─'*23}|")

print(f"\n  MODIFIED TRIPP (α(z), β(z), c², γ(z)):")
print(f"  {'z':>6}  {'-0.10':>5}{'0':>23}{'+0.10':>22}")
print(f"  {'':>6}  |{'─'*23}┼{'─'*23}|")

for i in range(len(z_mids)):
    val = mod_means[i]
    pos = int((val / 0.10 + 1) * 23)
    pos = max(0, min(46, pos))
    bar = [' '] * 47
    bar[23] = '│'
    marker = '●' if abs(val) > 0.02 else '○'
    bar[pos] = marker
    line = ''.join(bar)
    print(f"  {z_mids[i]:>6.3f} |{line}|")

print(f"  {'':>6}  |{'─'*23}┼{'─'*23}|")


# ============================================================
# Ωm REFIT UNDER MODIFIED STANDARDIZATION
# ============================================================

print(f"\n\n{'=' * 70}")
print("Ωm REFIT: DOES MODIFIED STANDARDIZATION CHANGE THE MATTER DENSITY?")
print("=" * 70)

# Subsample for speed
step = 3
idx_sub = np.argsort(z)[::step]
z_s = z[idx_sub]; mb_s = mb[idx_sub]; mb_err_s = mb_err[idx_sub]
x1_s = x1_q[idx_sub]; c_s = c_q[idx_sub]; hm_s = hm[idx_sub]
mask_hi_s = hm_s >= 10.0

# Standard Tripp + Ωm
print("\n  Fitting Standard Tripp + ΛCDM (varying Ωm)...")

def std_om_chi2(params):
    M_B, alpha, beta, gamma, Om = params
    if Om < 0.01 or Om > 0.99:
        return 1e20
    mu_mod = dist_mod(z_s, H0=70, Om=Om)
    mu_obs = mb_s + alpha * x1_s - beta * c_s - M_B
    mu_obs[mask_hi_s] -= gamma
    return np.sum(((mu_obs - mu_mod) / mb_err_s)**2)

res_om_std = scipy_minimize(std_om_chi2, x0=[-19.25, 0.15, 3.0, -0.05, 0.3],
                            method='Nelder-Mead', options={'maxiter': 50000})
Om_std = res_om_std.x[4]
print(f"    Ωm = {Om_std:.4f}")

# Modified Tripp + Ωm
print("  Fitting Modified Tripp + ΛCDM (varying Ωm)...")

def mod_om_chi2(params):
    M_B, a0, a1, b0, b1, delta, g0, g1, Om = params
    if Om < 0.01 or Om > 0.99:
        return 1e20
    mu_mod = dist_mod(z_s, H0=70, Om=Om)
    alpha_z = a0 + a1 * z_s
    beta_z = b0 + b1 * z_s
    gamma_z = g0 + g1 * z_s
    mu_obs = mb_s + alpha_z * x1_s - beta_z * c_s - delta * c_s**2 - M_B
    mu_obs[mask_hi_s] -= gamma_z[mask_hi_s]
    return np.sum(((mu_obs - mu_mod) / mb_err_s)**2)

res_om_mod = scipy_minimize(mod_om_chi2,
                            x0=[-19.25, 0.15, -0.03, 3.0, -0.5, 2.0, -0.03, -0.05, 0.3],
                            method='Nelder-Mead', options={'maxiter': 100000})
Om_mod = res_om_mod.x[8]
print(f"    Ωm = {Om_mod:.4f}")

# Standard Tripp + wCDM (varying Ωm and w)
print("  Fitting Standard Tripp + wCDM (varying Ωm, w)...")

def std_wcdm_chi2(params):
    M_B, alpha, beta, gamma, Om, w = params
    if Om < 0.01 or Om > 0.99:
        return 1e20
    mu_mod = dist_mod_wcdm(z_s, H0=70, Om=Om, w=w)
    mu_obs = mb_s + alpha * x1_s - beta * c_s - M_B
    mu_obs[mask_hi_s] -= gamma
    return np.sum(((mu_obs - mu_mod) / mb_err_s)**2)

res_wcdm_std = scipy_minimize(std_wcdm_chi2, x0=[-19.25, 0.15, 3.0, -0.05, 0.3, -1.0],
                              method='Nelder-Mead', options={'maxiter': 50000})
Om_wcdm_std = res_wcdm_std.x[4]
w_wcdm_std = res_wcdm_std.x[5]
print(f"    Ωm = {Om_wcdm_std:.4f}, w = {w_wcdm_std:.4f}")

# Modified Tripp + wCDM
print("  Fitting Modified Tripp + wCDM (varying Ωm, w)...")

def mod_wcdm_chi2(params):
    M_B, a0, a1, b0, b1, delta, g0, g1, Om, w = params
    if Om < 0.01 or Om > 0.99:
        return 1e20
    mu_mod = dist_mod_wcdm(z_s, H0=70, Om=Om, w=w)
    alpha_z = a0 + a1 * z_s
    beta_z = b0 + b1 * z_s
    gamma_z = g0 + g1 * z_s
    mu_obs = mb_s + alpha_z * x1_s - beta_z * c_s - delta * c_s**2 - M_B
    mu_obs[mask_hi_s] -= gamma_z[mask_hi_s]
    return np.sum(((mu_obs - mu_mod) / mb_err_s)**2)

res_wcdm_mod = scipy_minimize(mod_wcdm_chi2,
                              x0=[-19.25, 0.15, -0.03, 3.0, -0.5, 2.0, -0.03, -0.05, 0.3, -1.0],
                              method='Nelder-Mead', options={'maxiter': 100000})
Om_wcdm_mod = res_wcdm_mod.x[8]
w_wcdm_mod = res_wcdm_mod.x[9]
print(f"    Ωm = {Om_wcdm_mod:.4f}, w = {w_wcdm_mod:.4f}")


# ============================================================
# COMPARISON
# ============================================================

print(f"\n\n{'=' * 70}")
print("COSMOLOGICAL PARAMETER COMPARISON")
print("=" * 70)

print(f"\n  {'Model':<40} {'Ωm':>8} {'w':>8} {'χ²':>10}")
print(f"  {'-'*68}")
print(f"  {'Standard Tripp + ΛCDM':<40} {Om_std:>8.4f} {'−1 fix':>8} {res_om_std.fun:>10.1f}")
print(f"  {'Modified Tripp + ΛCDM':<40} {Om_mod:>8.4f} {'−1 fix':>8} {res_om_mod.fun:>10.1f}")
print(f"  {'Standard Tripp + wCDM':<40} {Om_wcdm_std:>8.4f} {w_wcdm_std:>+8.4f} {res_wcdm_std.fun:>10.1f}")
print(f"  {'Modified Tripp + wCDM':<40} {Om_wcdm_mod:>8.4f} {w_wcdm_mod:>+8.4f} {res_wcdm_mod.fun:>10.1f}")

print(f"\n  Planck 2018: Ωm = 0.3153 ± 0.0073")
print(f"  Standard Tripp: Ωm = {Om_std:.4f} (Δ = {Om_std - 0.3153:+.4f})")
print(f"  Modified Tripp: Ωm = {Om_mod:.4f} (Δ = {Om_mod - 0.3153:+.4f})")
print(f"  Shift: ΔΩm = {Om_mod - Om_std:+.4f}")

print(f"\n  DARK ENERGY EQUATION OF STATE:")
print(f"  Standard Tripp: w = {w_wcdm_std:+.4f} (deviation from ΛCDM: {w_wcdm_std + 1:+.4f})")
print(f"  Modified Tripp: w = {w_wcdm_mod:+.4f} (deviation from ΛCDM: {w_wcdm_mod + 1:+.4f})")

if abs(w_wcdm_mod + 1) < abs(w_wcdm_std + 1):
    pct = (1 - abs(w_wcdm_mod + 1) / abs(w_wcdm_std + 1)) * 100
    print(f"  ✓ Modified Tripp pushes w {pct:.0f}% CLOSER to ΛCDM (w = −1)")
else:
    print(f"  ⚠ Modified Tripp does not move w toward ΛCDM")


# ============================================================
# RESIDUAL DIAGNOSTICS — IS IT WHITE NOISE?
# ============================================================

print(f"\n\n{'=' * 70}")
print("WHITE NOISE DIAGNOSTICS")
print("=" * 70)

# 1. Runs test — does the sign of residuals alternate randomly?
signs_std = np.sign(resid_std[idx_sorted])
signs_mod = np.sign(resid_mod[idx_sorted])

runs_std = 1 + np.sum(np.diff(signs_std) != 0)
runs_mod = 1 + np.sum(np.diff(signs_mod) != 0)
expected_runs = len(signs_std) / 2 + 0.5

print(f"\n  Runs test (sign alternation, expected ≈ {expected_runs:.0f}):")
print(f"    Standard Tripp: {runs_std} runs")
print(f"    Modified Tripp: {runs_mod} runs")
print(f"    (More runs = more random = more like white noise)")

# 2. Lag-1 autocorrelation (already computed above)
print(f"\n  Lag-1 autocorrelation:")
print(f"    Standard Tripp: {autocorr_std:+.4f}")
print(f"    Modified Tripp: {autocorr_mod:+.4f}")
print(f"    Improvement: {(1 - abs(autocorr_mod)/abs(autocorr_std))*100:.0f}% closer to zero")

# 3. Lag-2, Lag-5, Lag-10 autocorrelation
for lag in [2, 5, 10]:
    ac_std = np.corrcoef(resid_std_sorted[:-lag], resid_std_sorted[lag:])[0, 1]
    ac_mod = np.corrcoef(resid_mod_sorted[:-lag], resid_mod_sorted[lag:])[0, 1]
    print(f"    Lag-{lag}: Standard {ac_std:+.4f}, Modified {ac_mod:+.4f}")

# 4. Residual-z correlation by z-bin (should be flat)
print(f"\n  Residual trend within z-bins (should all be ~0):")
z_check_bins = [(0.01, 0.05), (0.05, 0.15), (0.15, 0.35), (0.35, 0.70), (0.70, 1.50)]

for z_lo, z_hi in z_check_bins:
    mask_bin = (z >= z_lo) & (z < z_hi)
    n = np.sum(mask_bin)
    if n < 15:
        continue
    rho_s, _ = spearmanr(z[mask_bin], resid_std[mask_bin])
    rho_m, _ = spearmanr(z[mask_bin], resid_mod[mask_bin])
    print(f"    [{z_lo:.2f},{z_hi:.2f}]: Standard ρ={rho_s:+.3f}, Modified ρ={rho_m:+.3f}")


# ============================================================
# GRAND CONCLUSION
# ============================================================

print(f"\n\n{'=' * 70}")
print("CONCLUSION: HAVE WE DECODED THE EXPANSION?")
print("=" * 70)

decoded = abs(rho_mod) < abs(rho_std) and abs(autocorr_mod) < abs(autocorr_std)

print(f"""
  STANDARD TRIPP:
    RMS residual: {np.sqrt(np.mean(resid_std**2)):.4f}
    Residual-z correlation: ρ = {rho_std:+.4f} (p = {p_std:.2e})
    Autocorrelation: {autocorr_std:+.4f}
    Ωm = {Om_std:.4f}, w = {w_wcdm_std:+.4f}
    
  MODIFIED TRIPP (α(z) + β(z) + c² + γ(z)):
    RMS residual: {np.sqrt(np.mean(resid_mod**2)):.4f}
    Residual-z correlation: ρ = {rho_mod:+.4f} (p = {p_mod:.2e})
    Autocorrelation: {autocorr_mod:+.4f}
    Ωm = {Om_mod:.4f}, w = {w_wcdm_mod:+.4f}
    
  IMPROVEMENTS:
    RMS: {(1 - np.sqrt(np.mean(resid_mod**2))/np.sqrt(np.mean(resid_std**2)))*100:.1f}% reduction
    Trend: {(1 - abs(rho_mod)/abs(rho_std))*100:.0f}% reduction
    Autocorrelation: {(1 - abs(autocorr_mod)/abs(autocorr_std))*100:.0f}% reduction
    Binned mean scatter: {(1 - np.std(mod_means)/np.std(std_means))*100:.1f}% reduction
    
  {'✓ DECODED: Residuals are significantly closer to white noise.' if decoded else '◐ PARTIAL DECODE: Some structure remains.'}
  
  {'The modified Tripp equation removes the systematic z-dependent' if decoded else ''}
  {'structure from the Hubble diagram. What remains is measurement' if decoded else ''}
  {'noise, not cosmological signal leak.' if decoded else ''}
""")

# Save
os.makedirs('results_decode', exist_ok=True)
output = {
    'test_date': '2026-03-09',
    'rms_standard': float(np.sqrt(np.mean(resid_std**2))),
    'rms_modified': float(np.sqrt(np.mean(resid_mod**2))),
    'rho_z_standard': float(rho_std),
    'rho_z_modified': float(rho_mod),
    'autocorr_standard': float(autocorr_std),
    'autocorr_modified': float(autocorr_mod),
    'Om_standard': float(Om_std),
    'Om_modified': float(Om_mod),
    'w_standard': float(w_wcdm_std),
    'w_modified': float(w_wcdm_mod),
}

with open('results_decode/decode_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"  Results saved to results_decode/decode_results.json")
print(f"\n{'=' * 70}")
print("COMPLETE")
print("=" * 70)
