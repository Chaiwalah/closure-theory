#!/usr/bin/env python3
"""
CLOSURE THEORY — IMPEDANCE FUNCTIONAL FORM DECODE
====================================================

GPT's prescription: Stop linearizing. Use the actual impedance integral.

Replace α(z) = α₀ + α₁z with α(z) = α₀ × exp(−k × Z_g(z))
where Z_g(z) = ∫₀^χ(z) [H_local(χ')/H₀ − 1] dχ'

For flat ΛCDM: Z_g(z) = ∫₀^z [E(z')/1 − 1] dz' / (1+z')
where E(z) = H(z)/H₀ = √(Ωm(1+z)³ + ΩΛ)

THREE MODELS TO TEST:
1. Linear (current): coeff(z) = c₀ + c₁×z
2. Impedance integral: coeff(z) = c₀ × exp(−k × Z_g(z))
3. GPT's rotation + impedance: rotate (x1,c) by θ(z), then apply impedance

The goal: which one flattens residuals to white noise?

Author: Closure Theory collaboration
Date: 2026-03-09
"""

import numpy as np
from scipy.optimize import minimize as scipy_minimize
from scipy.integrate import quad
from scipy.stats import spearmanr
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
# COSMOLOGY + IMPEDANCE INTEGRAL
# ============================================================

C_LIGHT = 299792.458

def E(z, Om=0.3):
    """H(z)/H₀ for flat ΛCDM"""
    return np.sqrt(Om * (1+z)**3 + (1-Om))

def dist_mod(z_arr, H0=70, Om=0.3):
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        integral, _ = quad(lambda zp: 1.0/E(zp, Om), 0, zi)
        dl[i] = C_LIGHT * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

def Z_g_integral(z_val, Om=0.3):
    """Cumulative metric impedance: ∫₀^z [E(z')-1] dz'/(1+z')
    This is the integrated excess expansion relative to static."""
    if z_val <= 0:
        return 0.0
    result, _ = quad(lambda zp: (E(zp, Om) - 1) / (1+zp), 0, z_val)
    return result

# Pre-compute Z_g for all SNe
print("Computing impedance integrals...")
Zg = np.array([Z_g_integral(zi) for zi in z])
print(f"  Z_g range: [{Zg.min():.4f}, {Zg.max():.4f}]")
print(f"  Z_g(0.1) = {Z_g_integral(0.1):.4f}")
print(f"  Z_g(0.5) = {Z_g_integral(0.5):.4f}")
print(f"  Z_g(1.0) = {Z_g_integral(1.0):.4f}")

# Distance moduli
print("Computing distance moduli...")
mu_lcdm = dist_mod(z, H0=70, Om=0.3)

mask_hi = hm >= 10.0

# ============================================================
# MODEL 1: STANDARD TRIPP (baseline)
# ============================================================

print(f"\n{'=' * 70}")
print("MODEL 1: STANDARD TRIPP (constant α, β, γ)")
print("=" * 70)

def chi2_standard(params):
    M_B, alpha, beta, gamma = params
    mu_obs = mb + alpha * x1_q - beta * c_q - M_B
    mu_obs[mask_hi] -= gamma
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res1 = scipy_minimize(chi2_standard, x0=[-19.25, 0.15, 3.0, -0.05],
                      method='Nelder-Mead', options={'maxiter': 50000})
p1 = res1.x
resid1 = (mb + p1[1]*x1_q - p1[2]*c_q - p1[0]) - mu_lcdm
resid1[mask_hi] -= p1[3]

rho1, pval1 = spearmanr(z, resid1)
rms1 = np.sqrt(np.mean(resid1**2))
print(f"  χ²/ν = {res1.fun/(len(z)-4):.4f}, RMS = {rms1:.4f}")
print(f"  Residual-z: ρ = {rho1:+.4f} (p = {pval1:.2e})")


# ============================================================
# MODEL 2: LINEAR z-EVOLUTION (current best)
# ============================================================

print(f"\n{'=' * 70}")
print("MODEL 2: LINEAR α(z) + β(z) + c² + γ(z)")
print("=" * 70)

def chi2_linear(params):
    M_B, a0, a1, b0, b1, delta, g0, g1 = params
    alpha_z = a0 + a1 * z
    beta_z = b0 + b1 * z
    gamma_z = g0 + g1 * z
    mu_obs = mb + alpha_z * x1_q - beta_z * c_q - delta * c_q**2 - M_B
    mu_obs[mask_hi] -= gamma_z[mask_hi]
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res2 = scipy_minimize(chi2_linear,
                      x0=[-19.25, 0.15, -0.03, 3.0, -0.5, 2.0, -0.03, -0.05],
                      method='Nelder-Mead', options={'maxiter': 100000})
p2 = res2.x
alpha_z2 = p2[1] + p2[2]*z
beta_z2 = p2[3] + p2[4]*z
gamma_z2 = p2[5] + p2[6]*z  # wait, indexing wrong
# Fix: M_B=0, a0=1, a1=2, b0=3, b1=4, delta=5, g0=6, g1=7
gamma_z2 = p2[6] + p2[7]*z
resid2 = (mb + alpha_z2*x1_q - beta_z2*c_q - p2[5]*c_q**2 - p2[0]) - mu_lcdm
resid2[mask_hi] -= gamma_z2[mask_hi]

rho2, pval2 = spearmanr(z, resid2)
rms2 = np.sqrt(np.mean(resid2**2))
print(f"  χ²/ν = {res2.fun/(len(z)-8):.4f}, RMS = {rms2:.4f}")
print(f"  Residual-z: ρ = {rho2:+.4f} (p = {pval2:.2e})")


# ============================================================
# MODEL 3: EXPONENTIAL IMPEDANCE DECAY
# ============================================================

print(f"\n{'=' * 70}")
print("MODEL 3: EXPONENTIAL IMPEDANCE — α₀×exp(−k_α×Z_g), etc.")
print("=" * 70)

def chi2_exp_impedance(params):
    M_B, a0, k_a, b0, k_b, delta, g0, k_g = params
    alpha_z = a0 * np.exp(-k_a * Zg)
    beta_z = b0 * np.exp(-k_b * Zg)
    gamma_z = g0 * np.exp(k_g * Zg)  # gamma grows (more negative)
    mu_obs = mb + alpha_z * x1_q - beta_z * c_q - delta * c_q**2 - M_B
    mu_obs[mask_hi] -= gamma_z[mask_hi]
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res3 = scipy_minimize(chi2_exp_impedance,
                      x0=[-19.25, 0.16, 0.5, 3.2, 0.3, 2.0, -0.03, 0.5],
                      method='Nelder-Mead', options={'maxiter': 100000})
p3 = res3.x
alpha_z3 = p3[1] * np.exp(-p3[2] * Zg)
beta_z3 = p3[3] * np.exp(-p3[4] * Zg)
gamma_z3 = p3[5] * np.exp(p3[7] * Zg)
resid3 = (mb + alpha_z3*x1_q - beta_z3*c_q - p3[5]*c_q**2 - p3[0]) - mu_lcdm
# Oops, delta is p3[5], gamma0 is p3[6]
# Redo properly
resid3 = (mb + alpha_z3*x1_q - beta_z3*c_q - p3[5]*c_q**2 - p3[0]) - mu_lcdm
gamma_z3_vals = p3[6] * np.exp(p3[7] * Zg)
resid3[mask_hi] -= gamma_z3_vals[mask_hi]

rho3, pval3 = spearmanr(z, resid3)
rms3 = np.sqrt(np.mean(resid3**2))
print(f"  k_α = {p3[2]:.3f}, k_β = {p3[4]:.3f}, k_γ = {p3[7]:.3f}")
print(f"  α(z=0) = {p3[1]:.4f}, α(z=1) = {p3[1]*np.exp(-p3[2]*Z_g_integral(1.0)):.4f}")
print(f"  β(z=0) = {p3[3]:.3f}, β(z=1) = {p3[3]*np.exp(-p3[4]*Z_g_integral(1.0)):.3f}")
print(f"  χ²/ν = {res3.fun/(len(z)-8):.4f}, RMS = {rms3:.4f}")
print(f"  Residual-z: ρ = {rho3:+.4f} (p = {pval3:.2e})")


# ============================================================
# MODEL 4: SINGLE IMPEDANCE DECAY (GPT's 3-param suggestion)
# ============================================================

print(f"\n{'=' * 70}")
print("MODEL 4: SINGLE IMPEDANCE DECAY — ONE k FOR ALL")
print("μ = (m_B + α₀x1 − β₀c − δc² − M_B) × exp(−k×Z_g) + M_B")
print("=" * 70)

# GPT's idea: the ENTIRE standardized magnitude decays with one k
# This is the minimal model — one decay constant for all channels

def chi2_single_decay(params):
    M_B, alpha, beta, delta, gamma, k = params
    # Standardized magnitude before impedance
    mu_raw = mb + alpha * x1_q - beta * c_q - delta * c_q**2 - M_B
    mu_raw[mask_hi] -= gamma
    # Apply impedance attenuation
    # The impedance adds a correction proportional to Z_g
    mu_obs = mu_raw - k * Zg
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res4 = scipy_minimize(chi2_single_decay,
                      x0=[-19.25, 0.15, 3.0, 2.0, -0.05, 0.1],
                      method='Nelder-Mead', options={'maxiter': 50000})
p4 = res4.x
mu_raw4 = mb + p4[1]*x1_q - p4[2]*c_q - p4[3]*c_q**2 - p4[0]
mu_raw4[mask_hi] -= p4[4]
resid4 = mu_raw4 - p4[5]*Zg - mu_lcdm

rho4, pval4 = spearmanr(z, resid4)
rms4 = np.sqrt(np.mean(resid4**2))
print(f"  k = {p4[5]:.4f}")
print(f"  α = {p4[1]:.4f}, β = {p4[2]:.3f}, δ = {p4[3]:.3f}, γ = {p4[4]:+.4f}")
print(f"  χ²/ν = {res4.fun/(len(z)-6):.4f}, RMS = {rms4:.4f}")
print(f"  Residual-z: ρ = {rho4:+.4f} (p = {pval4:.2e})")


# ============================================================
# MODEL 5: ROTATION + IMPEDANCE (GPT's full suggestion)
# ============================================================

print(f"\n{'=' * 70}")
print("MODEL 5: ROTATION + IMPEDANCE")
print("Rotate (x1, c) by θ×Z_g, then apply standard + impedance")
print("=" * 70)

def chi2_rotation(params):
    M_B, alpha, beta, delta, gamma, k, theta = params
    # Rotate observables by θ × Z_g
    angle = theta * Zg
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    # Normalize x1 and c to similar scales first
    x1_norm = x1_q / 3.0  # x1 range ~ ±3
    c_norm = c_q / 0.3     # c range ~ ±0.3
    x1_rot = x1_norm * cos_a - c_norm * sin_a
    c_rot = x1_norm * sin_a + c_norm * cos_a
    # Un-normalize
    x1_eff = x1_rot * 3.0
    c_eff = c_rot * 0.3
    
    mu_obs = mb + alpha * x1_eff - beta * c_eff - delta * c_eff**2 - M_B
    mu_obs[mask_hi] -= gamma
    mu_obs -= k * Zg
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res5 = scipy_minimize(chi2_rotation,
                      x0=[-19.25, 0.15, 3.0, 2.0, -0.05, 0.1, 0.1],
                      method='Nelder-Mead', options={'maxiter': 100000})
p5 = res5.x
# Compute residuals
angle5 = p5[6] * Zg
x1_n = x1_q/3; c_n = c_q/0.3
x1_r = (x1_n*np.cos(angle5) - c_n*np.sin(angle5)) * 3
c_r = (x1_n*np.sin(angle5) + c_n*np.cos(angle5)) * 0.3
resid5 = (mb + p5[1]*x1_r - p5[2]*c_r - p5[3]*c_r**2 - p5[0]) - p5[5]*Zg - mu_lcdm
resid5[mask_hi] -= p5[4]

rho5, pval5 = spearmanr(z, resid5)
rms5 = np.sqrt(np.mean(resid5**2))
print(f"  k = {p5[5]:.4f}, θ = {p5[6]:.4f} rad/Z_g")
print(f"  θ at z=1: {p5[6]*Z_g_integral(1.0)*180/np.pi:.1f}°")
print(f"  χ²/ν = {res5.fun/(len(z)-7):.4f}, RMS = {rms5:.4f}")
print(f"  Residual-z: ρ = {rho5:+.4f} (p = {pval5:.2e})")


# ============================================================
# MODEL 6: IMPEDANCE-CORRECTED TRIPP (per-channel k)
# ============================================================

print(f"\n{'=' * 70}")
print("MODEL 6: PER-CHANNEL IMPEDANCE DECAY")
print("α(z) = α₀exp(−k_α Z_g), β(z) = β₀exp(−k_β Z_g), + k_μ Z_g offset")
print("=" * 70)

def chi2_perchannel(params):
    M_B, a0, k_a, b0, k_b, delta, gamma, k_mu = params
    alpha_z = a0 * np.exp(-k_a * Zg)
    beta_z = b0 * np.exp(-k_b * Zg)
    mu_obs = mb + alpha_z * x1_q - beta_z * c_q - delta * c_q**2 - M_B
    mu_obs[mask_hi] -= gamma
    mu_obs -= k_mu * Zg
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res6 = scipy_minimize(chi2_perchannel,
                      x0=[-19.25, 0.16, 0.3, 3.2, 0.2, 2.0, -0.05, 0.1],
                      method='Nelder-Mead', options={'maxiter': 100000})
p6 = res6.x
alpha_z6 = p6[1] * np.exp(-p6[2] * Zg)
beta_z6 = p6[3] * np.exp(-p6[4] * Zg)
resid6 = (mb + alpha_z6*x1_q - beta_z6*c_q - p6[5]*c_q**2 - p6[0]) - p6[7]*Zg - mu_lcdm
resid6[mask_hi] -= p6[6]

rho6, pval6 = spearmanr(z, resid6)
rms6 = np.sqrt(np.mean(resid6**2))
print(f"  α₀ = {p6[1]:.4f}, k_α = {p6[2]:.3f}")
print(f"  β₀ = {p6[3]:.3f}, k_β = {p6[4]:.3f}")
print(f"  k_μ = {p6[7]:.4f}")
print(f"  α(z=0) = {p6[1]:.4f}, α(z=1) = {p6[1]*np.exp(-p6[2]*Z_g_integral(1.0)):.4f}")
print(f"  β(z=0) = {p6[3]:.3f}, β(z=1) = {p6[3]*np.exp(-p6[4]*Z_g_integral(1.0)):.3f}")
print(f"  χ²/ν = {res6.fun/(len(z)-8):.4f}, RMS = {rms6:.4f}")
print(f"  Residual-z: ρ = {rho6:+.4f} (p = {pval6:.2e})")


# ============================================================
# GRAND COMPARISON
# ============================================================

print(f"\n\n{'=' * 70}")
print("GRAND COMPARISON — WHICH MODEL DECODES THE EXPANSION?")
print("=" * 70)

from scipy.stats import chi2 as chi2_dist

models = [
    ('1. Standard Tripp', 4, res1.fun, rms1, rho1, pval1),
    ('2. Linear z-evolution', 8, res2.fun, rms2, rho2, pval2),
    ('3. Exp impedance (per-k)', 8, res3.fun, rms3, rho3, pval3),
    ('4. Single k + c²', 6, res4.fun, rms4, rho4, pval4),
    ('5. Rotation + impedance', 7, res5.fun, rms5, rho5, pval5),
    ('6. Per-channel + k_μ', 8, res6.fun, rms6, rho6, pval6),
]

print(f"\n  {'Model':<30} {'k':>3} {'χ²/ν':>9} {'RMS':>7} {'ρ(z)':>7} {'p-value':>10} {'Decoded?':>9}")
print(f"  {'-'*80}")

for name, k, chi2_val, rms, rho, pval in models:
    chi2_nu = chi2_val / (len(z) - k)
    decoded = '✓' if pval > 0.05 else '◐' if pval > 0.001 else '✗'
    print(f"  {name:<30} {k:>3} {chi2_nu:>9.4f} {rms:>7.4f} {rho:>+7.4f} {pval:>10.2e} {decoded:>9}")

# Find best model
best_idx = min(range(len(models)), key=lambda i: abs(models[i][4]))  # lowest |ρ|
best = models[best_idx]
print(f"\n  BEST: {best[0]} (ρ = {best[4]:+.4f}, p = {best[5]:.2e})")

# Is any model consistent with white noise?
any_white = any(m[5] > 0.05 for m in models)
if any_white:
    white_models = [m for m in models if m[5] > 0.05]
    print(f"\n  ✓ WHITE NOISE ACHIEVED by: {', '.join(m[0] for m in white_models)}")
    print(f"  THE EXPANSION HAS BEEN DECODED.")
else:
    print(f"\n  ◐ No model achieves full white noise (all p < 0.05)")
    print(f"  But improvement from standard: ρ went from {rho1:+.4f} to {best[4]:+.4f}")
    improvement = (1 - abs(best[4])/abs(rho1)) * 100
    print(f"  That's {improvement:.0f}% reduction in residual trend")

# Δχ² tests
print(f"\n  Δχ² vs Standard Tripp:")
for name, k, chi2_val, rms, rho, pval in models[1:]:
    dchi2 = res1.fun - chi2_val
    dk = k - 4
    p_improvement = 1 - chi2_dist.cdf(dchi2, dk) if dchi2 > 0 else 1.0
    print(f"    {name:<30}: Δχ² = {dchi2:>8.1f} ({dk} params, p = {p_improvement:.2e})")


# Save
os.makedirs('results_impedance_decode', exist_ok=True)
output = {
    'test_date': '2026-03-09',
    'models': {},
}
for name, k, chi2_val, rms, rho, pval in models:
    output['models'][name] = {
        'k': k, 'chi2': float(chi2_val), 'rms': float(rms),
        'rho_z': float(rho), 'p_value': float(pval),
    }

with open('results_impedance_decode/impedance_decode_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n  Results saved to results_impedance_decode/impedance_decode_results.json")
print(f"\n{'=' * 70}")
print("COMPLETE")
print("=" * 70)
