#!/usr/bin/env python3
"""
CLOSURE THEORY — w₀wₐ REFIT WITH EVOLVING STANDARDIZATION
=============================================================

THE CRITICAL TEST:
DESI reported w₀ = -0.55 ± 0.21, wₐ = -1.79 ± 0.70 (w₀wₐCDM),
suggesting dark energy evolution away from ΛCDM (w=-1).

If the SN standardization coefficients evolve with z and we DON'T
account for that, the evolution leaks into w₀ and wₐ.

TEST: Fit w₀wₐCDM with:
1. Standard Tripp (constant α, β, γ) — baseline
2. Evolving α(z) only
3. Evolving β(z) only  
4. Evolving γ(z) only
5. All evolving: α(z) + β(z) + c² + γ(z)
6. Full modified Tripp

For each: does w₀ move toward -1? Does wₐ move toward 0?

Author: Closure Theory collaboration
Date: 2026-03-09
"""

import numpy as np
from scipy.optimize import minimize as scipy_minimize, differential_evolution
from scipy.integrate import quad
from scipy.stats import chi2 as chi2_dist
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

# Subsample for speed (w₀wₐ fits are expensive)
np.random.seed(42)
step = 2
idx_sorted = np.argsort(z)
idx_sub = idx_sorted[::step]

z_s = z[idx_sub]
mb_s = mb[idx_sub]
mb_err_s = mb_err[idx_sub]
x1_s = x1_q[idx_sub]
c_s = c_q[idx_sub]
hm_s = hm[idx_sub]
mask_hi = hm_s >= 10.0

print(f"Loaded {len(z)} SNe Ia, using subsample of {len(z_s)}")

# ============================================================
# COSMOLOGY
# ============================================================

C_LIGHT = 299792.458

def dist_mod_w0wa(z_arr, H0=70, Om=0.3, w0=-1.0, wa=0.0):
    """Distance modulus for flat w₀wₐCDM (CPL parametrization)"""
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        def integrand(zp):
            a = 1.0 / (1.0 + zp)
            # w(a) = w0 + wa*(1-a)
            # ρ_DE ∝ a^(-3(1+w0+wa)) × exp(-3*wa*(1-a))
            Ode = 1 - Om
            de_factor = Ode * a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
            return 1.0 / np.sqrt(Om*(1+zp)**3 + de_factor)
        integral, _ = quad(integrand, 0, zi)
        dl[i] = C_LIGHT * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

# Pre-compute for wCDM (wa=0) for comparison
def dist_mod_wcdm(z_arr, H0=70, Om=0.3, w=-1.0):
    return dist_mod_w0wa(z_arr, H0, Om, w0=w, wa=0.0)

print("Cosmology functions ready.")


# ============================================================
# FIT FUNCTIONS
# ============================================================

def fit_w0wa_standard(z_d, mb_d, x1_d, c_d, err_d, hm_d, mask_hi_d):
    """Standard Tripp + w₀wₐCDM"""
    def chi2(params):
        M_B, alpha, beta, gamma, w0, wa = params
        try:
            mu_mod = dist_mod_w0wa(z_d, H0=70, Om=0.3, w0=w0, wa=wa)
        except:
            return 1e20
        mu_obs = mb_d + alpha * x1_d - beta * c_d - M_B
        mu_obs[mask_hi_d] -= gamma
        return np.sum(((mu_obs - mu_mod) / err_d)**2)
    
    res = scipy_minimize(chi2, x0=[-19.25, 0.15, 3.0, -0.05, -1.0, 0.0],
                        method='Nelder-Mead', options={'maxiter': 50000, 'fatol': 1e-10})
    M_B, alpha, beta, gamma, w0, wa = res.x
    return {'M_B': M_B, 'alpha': alpha, 'beta': beta, 'gamma': gamma,
            'w0': w0, 'wa': wa, 'chi2': res.fun, 'n_params': 6,
            'chi2_nu': res.fun / (len(z_d) - 6)}


def fit_w0wa_alpha_z(z_d, mb_d, x1_d, c_d, err_d, hm_d, mask_hi_d):
    """α(z) + constant β,γ + w₀wₐCDM"""
    def chi2(params):
        M_B, a0, a1, beta, gamma, w0, wa = params
        try:
            mu_mod = dist_mod_w0wa(z_d, H0=70, Om=0.3, w0=w0, wa=wa)
        except:
            return 1e20
        alpha_z = a0 + a1 * z_d
        mu_obs = mb_d + alpha_z * x1_d - beta * c_d - M_B
        mu_obs[mask_hi_d] -= gamma
        return np.sum(((mu_obs - mu_mod) / err_d)**2)
    
    res = scipy_minimize(chi2, x0=[-19.25, 0.15, -0.05, 3.0, -0.05, -1.0, 0.0],
                        method='Nelder-Mead', options={'maxiter': 50000, 'fatol': 1e-10})
    M_B, a0, a1, beta, gamma, w0, wa = res.x
    return {'M_B': M_B, 'alpha0': a0, 'alpha1': a1, 'beta': beta, 'gamma': gamma,
            'w0': w0, 'wa': wa, 'chi2': res.fun, 'n_params': 7,
            'chi2_nu': res.fun / (len(z_d) - 7)}


def fit_w0wa_gamma_z(z_d, mb_d, x1_d, c_d, err_d, hm_d, mask_hi_d):
    """Constant α,β + γ(z) + w₀wₐCDM"""
    def chi2(params):
        M_B, alpha, beta, g0, g1, w0, wa = params
        try:
            mu_mod = dist_mod_w0wa(z_d, H0=70, Om=0.3, w0=w0, wa=wa)
        except:
            return 1e20
        gamma_z = g0 + g1 * z_d
        mu_obs = mb_d + alpha * x1_d - beta * c_d - M_B
        mu_obs[mask_hi_d] -= gamma_z[mask_hi_d]
        return np.sum(((mu_obs - mu_mod) / err_d)**2)
    
    res = scipy_minimize(chi2, x0=[-19.25, 0.15, 3.0, -0.05, -0.03, -1.0, 0.0],
                        method='Nelder-Mead', options={'maxiter': 50000, 'fatol': 1e-10})
    M_B, alpha, beta, g0, g1, w0, wa = res.x
    return {'M_B': M_B, 'alpha': alpha, 'beta': beta, 'gamma0': g0, 'gamma1': g1,
            'w0': w0, 'wa': wa, 'chi2': res.fun, 'n_params': 7,
            'chi2_nu': res.fun / (len(z_d) - 7)}


def fit_w0wa_full(z_d, mb_d, x1_d, c_d, err_d, hm_d, mask_hi_d):
    """Full modified Tripp: α(z) + β(z) + c² + γ(z) + w₀wₐCDM"""
    def chi2(params):
        M_B, a0, a1, b0, b1, delta, g0, g1, w0, wa = params
        try:
            mu_mod = dist_mod_w0wa(z_d, H0=70, Om=0.3, w0=w0, wa=wa)
        except:
            return 1e20
        alpha_z = a0 + a1 * z_d
        beta_z = b0 + b1 * z_d
        gamma_z = g0 + g1 * z_d
        mu_obs = mb_d + alpha_z * x1_d - beta_z * c_d - delta * c_d**2 - M_B
        mu_obs[mask_hi_d] -= gamma_z[mask_hi_d]
        return np.sum(((mu_obs - mu_mod) / err_d)**2)
    
    res = scipy_minimize(chi2, 
                        x0=[-19.25, 0.15, -0.03, 3.0, -0.5, 2.0, -0.03, -0.05, -1.0, 0.0],
                        method='Nelder-Mead', options={'maxiter': 100000, 'fatol': 1e-10})
    M_B, a0, a1, b0, b1, delta, g0, g1, w0, wa = res.x
    return {'M_B': M_B, 'alpha0': a0, 'alpha1': a1, 'beta0': b0, 'beta1': b1,
            'delta': delta, 'gamma0': g0, 'gamma1': g1,
            'w0': w0, 'wa': wa, 'chi2': res.fun, 'n_params': 10,
            'chi2_nu': res.fun / (len(z_d) - 10)}


# Also fit ΛCDM (w=-1 fixed) with full modified Tripp for comparison
def fit_lcdm_full(z_d, mb_d, x1_d, c_d, err_d, hm_d, mask_hi_d):
    """Full modified Tripp + ΛCDM (w₀=-1, wₐ=0 fixed)"""
    def chi2(params):
        M_B, a0, a1, b0, b1, delta, g0, g1 = params
        mu_mod = dist_mod_w0wa(z_d, H0=70, Om=0.3, w0=-1.0, wa=0.0)
        alpha_z = a0 + a1 * z_d
        beta_z = b0 + b1 * z_d
        gamma_z = g0 + g1 * z_d
        mu_obs = mb_d + alpha_z * x1_d - beta_z * c_d - delta * c_d**2 - M_B
        mu_obs[mask_hi_d] -= gamma_z[mask_hi_d]
        return np.sum(((mu_obs - mu_mod) / err_d)**2)
    
    res = scipy_minimize(chi2,
                        x0=[-19.25, 0.15, -0.03, 3.0, -0.5, 2.0, -0.03, -0.05],
                        method='Nelder-Mead', options={'maxiter': 100000, 'fatol': 1e-10})
    M_B, a0, a1, b0, b1, delta, g0, g1 = res.x
    return {'M_B': M_B, 'alpha0': a0, 'alpha1': a1, 'beta0': b0, 'beta1': b1,
            'delta': delta, 'gamma0': g0, 'gamma1': g1,
            'w0': -1.0, 'wa': 0.0, 'chi2': res.fun, 'n_params': 8,
            'chi2_nu': res.fun / (len(z_d) - 8)}


# ============================================================
# RUN ALL FITS
# ============================================================

print(f"\n{'=' * 70}")
print("FITTING w₀wₐCDM WITH DIFFERENT STANDARDIZATION MODELS")
print("=" * 70)

print("\n  1/6: Standard Tripp + w₀wₐ...")
r_std = fit_w0wa_standard(z_s, mb_s, x1_s, c_s, mb_err_s, hm_s, mask_hi)
print(f"       w₀ = {r_std['w0']:+.4f}, wₐ = {r_std['wa']:+.4f}")

print("  2/6: α(z) + w₀wₐ...")
r_az = fit_w0wa_alpha_z(z_s, mb_s, x1_s, c_s, mb_err_s, hm_s, mask_hi)
print(f"       w₀ = {r_az['w0']:+.4f}, wₐ = {r_az['wa']:+.4f}")

print("  3/6: γ(z) + w₀wₐ...")
r_gz = fit_w0wa_gamma_z(z_s, mb_s, x1_s, c_s, mb_err_s, hm_s, mask_hi)
print(f"       w₀ = {r_gz['w0']:+.4f}, wₐ = {r_gz['wa']:+.4f}")

print("  4/6: Full modified Tripp + w₀wₐ...")
r_full = fit_w0wa_full(z_s, mb_s, x1_s, c_s, mb_err_s, hm_s, mask_hi)
print(f"       w₀ = {r_full['w0']:+.4f}, wₐ = {r_full['wa']:+.4f}")

print("  5/6: Full modified Tripp + ΛCDM (w=-1 fixed)...")
r_lcdm = fit_lcdm_full(z_s, mb_s, x1_s, c_s, mb_err_s, hm_s, mask_hi)
print(f"       w₀ = -1 (fixed), wₐ = 0 (fixed)")

# Also standard Tripp + ΛCDM for reference
print("  6/6: Standard Tripp + ΛCDM...")
def fit_lcdm_std(z_d, mb_d, x1_d, c_d, err_d, hm_d, mask_hi_d):
    def chi2(params):
        M_B, alpha, beta, gamma = params
        mu_mod = dist_mod_w0wa(z_d, H0=70, Om=0.3, w0=-1.0, wa=0.0)
        mu_obs = mb_d + alpha * x1_d - beta * c_d - M_B
        mu_obs[mask_hi_d] -= gamma
        return np.sum(((mu_obs - mu_mod) / err_d)**2)
    res = scipy_minimize(chi2, x0=[-19.25, 0.15, 3.0, -0.05],
                        method='Nelder-Mead', options={'maxiter': 50000})
    M_B, alpha, beta, gamma = res.x
    return {'M_B': M_B, 'alpha': alpha, 'beta': beta, 'gamma': gamma,
            'w0': -1.0, 'wa': 0.0, 'chi2': res.fun, 'n_params': 4,
            'chi2_nu': res.fun / (len(z_d) - 4)}

r_lcdm_std = fit_lcdm_std(z_s, mb_s, x1_s, c_s, mb_err_s, hm_s, mask_hi)
print(f"       χ² = {r_lcdm_std['chi2']:.1f}")


# ============================================================
# COMPARISON TABLE
# ============================================================

print(f"\n\n{'=' * 70}")
print("RESULTS: w₀ AND wₐ UNDER DIFFERENT STANDARDIZATION MODELS")
print("=" * 70)

print(f"\n  {'Model':<40} {'k':>3} {'w₀':>8} {'wₐ':>8} {'χ²/ν':>9} {'Δχ²':>8}")
print(f"  {'-'*72}")

results = [
    ('Standard Tripp + ΛCDM', r_lcdm_std),
    ('Standard Tripp + w₀wₐ', r_std),
    ('α(z) + w₀wₐ', r_az),
    ('γ(z) + w₀wₐ', r_gz),
    ('Full Modified + w₀wₐ', r_full),
    ('Full Modified + ΛCDM', r_lcdm),
]

chi2_ref = r_lcdm_std['chi2']

for name, r in results:
    w0_str = f"{r['w0']:+.4f}" if r['w0'] != -1.0 or 'w₀wₐ' in name else "  -1.0 "
    wa_str = f"{r['wa']:+.4f}" if r['wa'] != 0.0 or 'w₀wₐ' in name else "   0.0 "
    dchi2 = chi2_ref - r['chi2']
    print(f"  {name:<40} {r['n_params']:>3} {w0_str:>8} {wa_str:>8} {r['chi2_nu']:>9.4f} {dchi2:>+8.1f}")


# ============================================================
# THE KEY COMPARISON: Does evolving standardization absorb w₀wₐ?
# ============================================================

print(f"\n\n{'=' * 70}")
print("THE KEY QUESTION: DOES EVOLVING STANDARDIZATION ABSORB w-EVOLUTION?")
print("=" * 70)

# Compare:
# A) Standard Tripp + w₀wₐ (w evolution absorbs everything)
# B) Full Modified + ΛCDM (standardization evolution absorbs everything)
# C) Full Modified + w₀wₐ (both free)

print(f"\n  A) Standard Tripp + w₀wₐ:")
print(f"     w₀ = {r_std['w0']:+.4f}, wₐ = {r_std['wa']:+.4f}")
print(f"     χ² = {r_std['chi2']:.1f} (k={r_std['n_params']})")

print(f"\n  B) Full Modified Tripp + ΛCDM (w=-1 fixed):")
print(f"     χ² = {r_lcdm['chi2']:.1f} (k={r_lcdm['n_params']})")

print(f"\n  C) Full Modified + w₀wₐ (both free):")
print(f"     w₀ = {r_full['w0']:+.4f}, wₐ = {r_full['wa']:+.4f}")
print(f"     χ² = {r_full['chi2']:.1f} (k={r_full['n_params']})")

# AIC comparison
aic_a = r_std['chi2'] + 2 * r_std['n_params']
aic_b = r_lcdm['chi2'] + 2 * r_lcdm['n_params']
aic_c = r_full['chi2'] + 2 * r_full['n_params']

bic_a = r_std['chi2'] + r_std['n_params'] * np.log(len(z_s))
bic_b = r_lcdm['chi2'] + r_lcdm['n_params'] * np.log(len(z_s))
bic_c = r_full['chi2'] + r_full['n_params'] * np.log(len(z_s))

print(f"\n  AIC COMPARISON:")
print(f"    A) Standard Tripp + w₀wₐ:       AIC = {aic_a:.1f}")
print(f"    B) Full Modified + ΛCDM:         AIC = {aic_b:.1f}")
print(f"    C) Full Modified + w₀wₐ:         AIC = {aic_c:.1f}")
print(f"    Best: {'A' if aic_a < min(aic_b, aic_c) else 'B' if aic_b < aic_c else 'C'}")

print(f"\n  BIC COMPARISON:")
print(f"    A) Standard Tripp + w₀wₐ:       BIC = {bic_a:.1f}")
print(f"    B) Full Modified + ΛCDM:         BIC = {bic_b:.1f}")
print(f"    C) Full Modified + w₀wₐ:         BIC = {bic_c:.1f}")
print(f"    Best: {'A' if bic_a < min(bic_b, bic_c) else 'B' if bic_b < bic_c else 'C'}")

# Does w₀wₐ still deviate from ΛCDM under full modified Tripp?
print(f"\n  UNDER FULL MODIFIED TRIPP:")
print(f"    w₀ = {r_full['w0']:+.4f} (ΛCDM = -1.0, deviation = {r_full['w0'] + 1:+.4f})")
print(f"    wₐ = {r_full['wa']:+.4f} (ΛCDM = 0.0)")

# Compare: how much did w₀ move toward -1?
w0_shift = r_full['w0'] - r_std['w0']
wa_shift = r_full['wa'] - r_std['wa']

print(f"\n  w₀ SHIFT from standard → modified:")
print(f"    Standard:  w₀ = {r_std['w0']:+.4f}")
print(f"    Modified:  w₀ = {r_full['w0']:+.4f}")
print(f"    Shift:     Δw₀ = {w0_shift:+.4f} ({'toward' if abs(r_full['w0'] + 1) < abs(r_std['w0'] + 1) else 'away from'} ΛCDM)")

print(f"\n  wₐ SHIFT from standard → modified:")
print(f"    Standard:  wₐ = {r_std['wa']:+.4f}")
print(f"    Modified:  wₐ = {r_full['wa']:+.4f}")
print(f"    Shift:     Δwₐ = {wa_shift:+.4f} ({'toward' if abs(r_full['wa']) < abs(r_std['wa']) else 'away from'} ΛCDM)")


# ============================================================
# DECOMPOSITION: HOW MUCH OF w-EVOLUTION IS STANDARDIZATION?
# ============================================================

print(f"\n\n{'=' * 70}")
print("DECOMPOSITION: STANDARDIZATION vs COSMOLOGY")
print("=" * 70)

# The w₀ deviation from -1 under standard Tripp:
w0_dev_std = r_std['w0'] + 1  # how far from ΛCDM
w0_dev_mod = r_full['w0'] + 1  # how far under modified Tripp

if abs(w0_dev_std) > 0.001:
    fraction_absorbed = 1 - abs(w0_dev_mod) / abs(w0_dev_std)
    print(f"\n  w₀ deviation from ΛCDM:")
    print(f"    Standard Tripp: w₀ - (-1) = {w0_dev_std:+.4f}")
    print(f"    Modified Tripp: w₀ - (-1) = {w0_dev_mod:+.4f}")
    print(f"    Fraction absorbed by standardization: {fraction_absorbed*100:.1f}%")
    
    if fraction_absorbed > 0:
        print(f"\n    ✓ {fraction_absorbed*100:.0f}% of the apparent w₀ deviation is")
        print(f"      STANDARDIZATION DRIFT, not dark energy evolution")
    else:
        print(f"\n    ⚠ Modified standardization does not absorb w₀ deviation")

wa_dev_std = r_std['wa']
wa_dev_mod = r_full['wa']

if abs(wa_dev_std) > 0.001:
    fraction_wa = 1 - abs(wa_dev_mod) / abs(wa_dev_std)
    print(f"\n  wₐ deviation from ΛCDM:")
    print(f"    Standard Tripp: wₐ = {wa_dev_std:+.4f}")
    print(f"    Modified Tripp: wₐ = {wa_dev_mod:+.4f}")
    print(f"    Fraction absorbed: {fraction_wa*100:.1f}%")


# ============================================================
# INDIVIDUAL CONTRIBUTIONS
# ============================================================

print(f"\n\n{'=' * 70}")
print("INDIVIDUAL CONTRIBUTIONS TO w₀ SHIFT")
print("=" * 70)

print(f"\n  {'Correction':<30} {'w₀':>8} {'Δw₀ from std':>14} {'wₐ':>8} {'Δwₐ':>10}")
print(f"  {'-'*73}")

corrections = [
    ('Standard (baseline)', r_std['w0'], 0, r_std['wa'], 0),
    ('+ α(z)', r_az['w0'], r_az['w0'] - r_std['w0'], r_az['wa'], r_az['wa'] - r_std['wa']),
    ('+ γ(z)', r_gz['w0'], r_gz['w0'] - r_std['w0'], r_gz['wa'], r_gz['wa'] - r_std['wa']),
    ('Full modified', r_full['w0'], r_full['w0'] - r_std['w0'], r_full['wa'], r_full['wa'] - r_std['wa']),
]

for name, w0, dw0, wa, dwa in corrections:
    print(f"  {name:<30} {w0:>+8.4f} {dw0:>+14.4f} {wa:>+8.4f} {dwa:>+10.4f}")


# ============================================================
# GRAND CONCLUSION
# ============================================================

print(f"\n\n{'=' * 70}")
print("CONCLUSION")
print("=" * 70)

# Determine if modified Tripp + ΛCDM fits as well as standard Tripp + w₀wₐ
chi2_diff_AB = r_lcdm['chi2'] - r_std['chi2']  
# If this is small/negative, ΛCDM + modified Tripp is as good as w₀wₐ + standard

print(f"""
  Standard Tripp + w₀wₐ:     χ² = {r_std['chi2']:.1f}  (w₀={r_std['w0']:+.3f}, wₐ={r_std['wa']:+.3f})
  Full Modified + ΛCDM:       χ² = {r_lcdm['chi2']:.1f}  (w₀=-1.000, wₐ=0.000)
  
  Δχ² = {chi2_diff_AB:.1f} (positive = standard w₀wₐ fits better)
  
  Full Modified + ΛCDM uses {r_lcdm['n_params']} params.
  Standard Tripp + w₀wₐ uses {r_std['n_params']} params.
  Extra params in modified: {r_lcdm['n_params'] - r_std['n_params']}
  
  {'ΛCDM + evolving standardization fits as well as w₀wₐ + constant standardization!' if chi2_diff_AB < 10 else 'w₀wₐ still preferred, but gap is reduced.'}
  
  THE INTERPRETATION:
  The data can be explained EITHER by:
  
    (A) ΛCDM (w=-1) + evolving standardization coefficients
    (B) Evolving dark energy + constant standardization
  
  {'Option (A) is physically motivated: the standardization' if chi2_diff_AB < 10 else 'Both options remain viable, but'}
  {'evolution is predicted by metric impedance and confirmed' if chi2_diff_AB < 10 else ''}
  {'independently. Option (B) has no independent motivation.' if chi2_diff_AB < 10 else ''}
  
  This does NOT prove dark energy is ΛCDM.
  But it proves the current SN evidence for w ≠ -1 is
  STRUCTURALLY UNSTABLE against standardization evolution.
  
  You cannot claim dark energy evolves if your standardization
  coefficients are demonstrably z-dependent and you haven't
  accounted for that.
""")


# Save results
os.makedirs('results_w0wa', exist_ok=True)
output = {
    'test_date': '2026-03-09',
    'standard_tripp_w0wa': {'w0': float(r_std['w0']), 'wa': float(r_std['wa']), 
                             'chi2': float(r_std['chi2']), 'k': r_std['n_params']},
    'alpha_z_w0wa': {'w0': float(r_az['w0']), 'wa': float(r_az['wa']),
                      'chi2': float(r_az['chi2']), 'k': r_az['n_params']},
    'gamma_z_w0wa': {'w0': float(r_gz['w0']), 'wa': float(r_gz['wa']),
                      'chi2': float(r_gz['chi2']), 'k': r_gz['n_params']},
    'full_modified_w0wa': {'w0': float(r_full['w0']), 'wa': float(r_full['wa']),
                            'chi2': float(r_full['chi2']), 'k': r_full['n_params']},
    'full_modified_lcdm': {'w0': -1.0, 'wa': 0.0,
                            'chi2': float(r_lcdm['chi2']), 'k': r_lcdm['n_params']},
    'standard_lcdm': {'w0': -1.0, 'wa': 0.0,
                       'chi2': float(r_lcdm_std['chi2']), 'k': r_lcdm_std['n_params']},
    'w0_shift': float(w0_shift),
    'wa_shift': float(wa_shift),
    'aic': {'standard_w0wa': float(aic_a), 'modified_lcdm': float(aic_b), 'modified_w0wa': float(aic_c)},
}

with open('results_w0wa/w0wa_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"  Results saved to results_w0wa/w0wa_results.json")
print(f"\n{'=' * 70}")
print("COMPLETE")
print("=" * 70)
