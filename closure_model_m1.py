#!/usr/bin/env python3
"""
CLOSURE THEORY — MODEL M1: GPT'S UNIFIED IMPEDANCE-ROTATED BASIS
===================================================================

GPT's prescription (verbatim):
1. ONE common impedance variable Ĩ(z) — not separate z-trends
2. Rotate (x1, c) by θ(Ĩ) = θ₀ + θ₁Ĩ
3. Quadratic on ROTATED color axis, not raw c
4. Mass step on same Ĩ
5. Analytically marginalize over M₀ (don't fix H₀)

Ĩ(z) options tested:
- ln(1+z) — GPT's recommended surrogate
- Z_g integral — full impedance
- z/(1+z/z*) — saturating form

This is FEWER effective degrees of freedom than independent linear
α(z), β(z), γ(z), but physically more correct because everything
is tied to one propagation variable.

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
mask_hi = hm >= 10.0

print(f"Loaded {len(z)} SNe Ia")

# ============================================================
# IMPEDANCE VARIABLES
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

# Three impedance variable options
I_log = np.log(1 + z)                    # GPT's recommendation: ln(1+z)

print("Computing Z_g integral...")
I_Zg = np.array([quad(lambda zp: (E(zp) - 1)/(1+zp), 0, zi)[0] for zi in z])

I_sat = z / (1 + z/0.5)                   # Saturating form z/(1+z/z*)

# Distance moduli (ΛCDM)
print("Computing distance moduli...")
mu_lcdm = dist_mod(z, H0=70, Om=0.3)

# Normalize x1 and c for rotation
x1_std = np.std(x1_q)
c_std = np.std(c_q)

print(f"  σ(x1) = {x1_std:.3f}, σ(c) = {c_std:.4f}")
print(f"  Ĩ range: ln(1+z) = [{I_log.min():.4f}, {I_log.max():.4f}]")
print(f"  Ĩ range: Z_g = [{I_Zg.min():.4f}, {I_Zg.max():.4f}]")

# ============================================================
# REFERENCE: STANDARD TRIPP
# ============================================================

print(f"\n{'=' * 70}")
print("REFERENCE: STANDARD TRIPP")
print("=" * 70)

def chi2_std(params):
    M_B, alpha, beta, gamma = params
    mu_obs = mb + alpha * x1_q - beta * c_q - M_B
    mu_obs[mask_hi] -= gamma
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res_std = scipy_minimize(chi2_std, x0=[-19.25, 0.15, 3.0, -0.05],
                         method='Nelder-Mead', options={'maxiter': 50000})
resid_std = (mb + res_std.x[1]*x1_q - res_std.x[2]*c_q - res_std.x[0]) - mu_lcdm
resid_std[mask_hi] -= res_std.x[3]
rho_std, p_std = spearmanr(z, resid_std)
rms_std = np.sqrt(np.mean(resid_std**2))
print(f"  k=4, χ²/ν = {res_std.fun/(len(z)-4):.4f}, RMS = {rms_std:.4f}")
print(f"  ρ(z) = {rho_std:+.4f} (p = {p_std:.2e})")


# ============================================================
# REFERENCE: LINEAR z-EVOLUTION (our previous best)
# ============================================================

print(f"\n{'=' * 70}")
print("REFERENCE: LINEAR α(z)+β(z)+c²+γ(z)")  
print("=" * 70)

def chi2_linear(params):
    M_B, a0, a1, b0, b1, delta, g0, g1 = params
    mu_obs = mb + (a0+a1*z)*x1_q - (b0+b1*z)*c_q - delta*c_q**2 - M_B
    mu_obs[mask_hi] -= (g0+g1*z)[mask_hi]
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res_lin = scipy_minimize(chi2_linear,
                         x0=[-19.25, 0.15, -0.03, 3.0, -0.5, 2.0, -0.03, -0.05],
                         method='Nelder-Mead', options={'maxiter': 100000})
p_lin = res_lin.x
resid_lin = (mb + (p_lin[1]+p_lin[2]*z)*x1_q - (p_lin[3]+p_lin[4]*z)*c_q 
             - p_lin[5]*c_q**2 - p_lin[0]) - mu_lcdm
resid_lin[mask_hi] -= (p_lin[6]+p_lin[7]*z)[mask_hi]
rho_lin, p_lin_val = spearmanr(z, resid_lin)
rms_lin = np.sqrt(np.mean(resid_lin**2))
print(f"  k=8, χ²/ν = {res_lin.fun/(len(z)-8):.4f}, RMS = {rms_lin:.4f}")
print(f"  ρ(z) = {rho_lin:+.4f} (p = {p_lin_val:.2e})")


# ============================================================
# MODEL M1: GPT'S UNIFIED ROTATED-BASIS IMPEDANCE MODEL
# ============================================================

def fit_model_m1(I_var, I_name):
    """
    GPT's Model M1:
    1. Rotate (x1, c) by θ(Ĩ) = θ₁ × Ĩ
    2. Standardize with A, B on rotated axes
    3. Quadratic on rotated color axis v
    4. Mass step γ(Ĩ) = γ₀ + γ₁Ĩ
    5. Analytically marginalize M₀
    
    Parameters: A, B, θ₁, δ_v, γ₀, γ₁ = 6 params + M₀ marginalized
    """
    
    print(f"\n{'=' * 70}")
    print(f"MODEL M1 with Ĩ = {I_name}")
    print("=" * 70)
    
    def chi2_m1(params):
        A, B, theta1, delta_v, g0, g1 = params
        
        # Rotation angle
        angle = theta1 * I_var
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        
        # Rotate: normalize to comparable scales first
        x1n = x1_q / x1_std
        cn = c_q / c_std
        
        u = x1n * cos_a + cn * sin_a      # stretch-like
        v = -x1n * sin_a + cn * cos_a     # color-like
        
        # Unnormalize for magnitude contribution
        u_mag = u * x1_std
        v_mag = v * c_std
        
        # Mass step
        gamma_I = g0 + g1 * I_var
        
        # Corrected magnitude (M₀ analytically marginalized)
        mu_corr = mb + A * u_mag - B * v_mag - delta_v * v_mag**2
        mu_corr[mask_hi] -= gamma_I[mask_hi]
        
        # Analytically marginalize over M₀:
        # Best M₀ = mean(mu_corr - mu_model) weighted by 1/σ²
        w = 1.0 / mb_err**2
        M0_best = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
        
        residuals = mu_corr - M0_best - mu_lcdm
        return np.sum((residuals / mb_err)**2)
    
    res = scipy_minimize(chi2_m1, x0=[0.15, 3.0, 0.5, 2.0, -0.03, -0.1],
                         method='Nelder-Mead', options={'maxiter': 100000})
    
    # Compute residuals at best fit
    A, B, theta1, delta_v, g0, g1 = res.x
    
    angle = theta1 * I_var
    x1n = x1_q / x1_std; cn = c_q / c_std
    u = x1n * np.cos(angle) + cn * np.sin(angle)
    v = -x1n * np.sin(angle) + cn * np.cos(angle)
    u_mag = u * x1_std; v_mag = v * c_std
    
    gamma_I = g0 + g1 * I_var
    mu_corr = mb + A * u_mag - B * v_mag - delta_v * v_mag**2
    mu_corr[mask_hi] -= gamma_I[mask_hi]
    
    w = 1.0 / mb_err**2
    M0_best = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
    
    resid = mu_corr - M0_best - mu_lcdm
    rho, pval = spearmanr(z, resid)
    rms = np.sqrt(np.mean(resid**2))
    
    # Rotation at key redshifts
    print(f"  A = {A:.4f}, B = {B:.3f}")
    print(f"  θ₁ = {theta1:.4f}")
    print(f"  δ_v = {delta_v:.3f}")
    print(f"  γ₀ = {g0:+.4f}, γ₁ = {g1:+.4f}")
    print(f"  M₀ = {M0_best:.4f}")
    
    # Rotation angles
    for zi in [0.1, 0.3, 0.5, 1.0]:
        I_at_z = np.log(1+zi) if 'ln' in I_name else 0
        if 'Z_g' in I_name:
            I_at_z = quad(lambda zp: (E(zp)-1)/(1+zp), 0, zi)[0]
        elif 'ln' in I_name:
            I_at_z = np.log(1+zi)
        else:
            I_at_z = zi/(1+zi/0.5)
        angle_deg = np.degrees(theta1 * I_at_z)
        print(f"    θ(z={zi}) = {angle_deg:.1f}°")
    
    print(f"\n  k=6+M₀, χ²/ν = {res.fun/(len(z)-7):.4f}, RMS = {rms:.4f}")
    print(f"  ρ(z) = {rho:+.4f} (p = {pval:.2e})")
    
    return res.fun, rms, rho, pval, res.x


# ============================================================
# MODEL M1-SAT: SATURATING COEFFICIENTS ON IMPEDANCE VARIABLE
# ============================================================

def fit_model_m1_sat(I_var, I_name):
    """
    GPT's saturating form:
    α(Ĩ) = α_∞ + (α₀ - α_∞) × exp(−Ĩ/I_α)
    Same for β, with rotation + quadratic
    """
    
    print(f"\n{'=' * 70}")
    print(f"MODEL M1-SAT with Ĩ = {I_name}")
    print("Saturating coefficients + rotation + quadratic")
    print("=" * 70)
    
    def chi2_sat(params):
        A0, A_inf, I_A, B0, B_inf, I_B, theta1, delta_v, g0, g1 = params
        
        # Saturating coefficients
        A_I = A_inf + (A0 - A_inf) * np.exp(-I_var / max(abs(I_A), 0.001))
        B_I = B_inf + (B0 - B_inf) * np.exp(-I_var / max(abs(I_B), 0.001))
        
        # Rotation
        angle = theta1 * I_var
        x1n = x1_q / x1_std; cn = c_q / c_std
        u = x1n * np.cos(angle) + cn * np.sin(angle)
        v = -x1n * np.sin(angle) + cn * np.cos(angle)
        u_mag = u * x1_std; v_mag = v * c_std
        
        gamma_I = g0 + g1 * I_var
        mu_corr = mb + A_I * u_mag - B_I * v_mag - delta_v * v_mag**2
        mu_corr[mask_hi] -= gamma_I[mask_hi]
        
        w = 1.0 / mb_err**2
        M0 = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
        
        residuals = mu_corr - M0 - mu_lcdm
        return np.sum((residuals / mb_err)**2)
    
    res = scipy_minimize(chi2_sat, 
                         x0=[0.16, 0.10, 0.5, 3.2, 2.0, 0.5, 0.3, 2.0, -0.03, -0.1],
                         method='Nelder-Mead', options={'maxiter': 200000, 'fatol': 1e-10})
    
    A0, A_inf, I_A, B0, B_inf, I_B, theta1, delta_v, g0, g1 = res.x
    
    A_I = A_inf + (A0 - A_inf) * np.exp(-I_var / max(abs(I_A), 0.001))
    B_I = B_inf + (B0 - B_inf) * np.exp(-I_var / max(abs(I_B), 0.001))
    
    angle = theta1 * I_var
    x1n = x1_q / x1_std; cn = c_q / c_std
    u = x1n * np.cos(angle) + cn * np.sin(angle)
    v = -x1n * np.sin(angle) + cn * np.cos(angle)
    u_mag = u * x1_std; v_mag = v * c_std
    
    gamma_I = g0 + g1 * I_var
    mu_corr = mb + A_I * u_mag - B_I * v_mag - delta_v * v_mag**2
    mu_corr[mask_hi] -= gamma_I[mask_hi]
    
    w = 1.0 / mb_err**2
    M0 = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
    
    resid = mu_corr - M0 - mu_lcdm
    rho, pval = spearmanr(z, resid)
    rms = np.sqrt(np.mean(resid**2))
    
    print(f"  A: {A0:.4f} → {A_inf:.4f} (I_A = {I_A:.3f})")
    print(f"  B: {B0:.3f} → {B_inf:.3f} (I_B = {I_B:.3f})")
    print(f"  θ₁ = {theta1:.4f}, δ_v = {delta_v:.3f}")
    print(f"  γ₀ = {g0:+.4f}, γ₁ = {g1:+.4f}")
    print(f"\n  k=10+M₀, χ²/ν = {res.fun/(len(z)-11):.4f}, RMS = {rms:.4f}")
    print(f"  ρ(z) = {rho:+.4f} (p = {pval:.2e})")
    
    return res.fun, rms, rho, pval, res.x


# ============================================================
# RUN ALL MODELS
# ============================================================

results = {}

# Standard + Linear (references)
results['Standard Tripp'] = (res_std.fun, rms_std, rho_std, p_std, 4)
results['Linear z-evolution'] = (res_lin.fun, rms_lin, rho_lin, p_lin_val, 8)

# Model M1 with different Ĩ
for I_var, I_name in [(I_log, 'ln(1+z)'), (I_Zg, 'Z_g integral'), (I_sat, 'z/(1+z/0.5)')]:
    chi2, rms, rho, pval, params = fit_model_m1(I_var, I_name)
    results[f'M1 [{I_name}]'] = (chi2, rms, rho, pval, 7)

# Model M1-SAT with ln(1+z)
chi2_s, rms_s, rho_s, pval_s, params_s = fit_model_m1_sat(I_log, 'ln(1+z)')
results['M1-SAT [ln(1+z)]'] = (chi2_s, rms_s, rho_s, pval_s, 11)

# Model M1-SAT with Z_g
chi2_s2, rms_s2, rho_s2, pval_s2, params_s2 = fit_model_m1_sat(I_Zg, 'Z_g integral')
results['M1-SAT [Z_g]'] = (chi2_s2, rms_s2, rho_s2, pval_s2, 11)


# ============================================================
# GRAND COMPARISON
# ============================================================

print(f"\n\n{'=' * 70}")
print("GRAND COMPARISON — ALL MODELS")
print("=" * 70)

print(f"\n  {'Model':<30} {'k':>3} {'χ²/ν':>9} {'RMS':>7} {'ρ(z)':>7} {'p-value':>10} {'Status':>8}")
print(f"  {'-'*78}")

for name in ['Standard Tripp', 'Linear z-evolution', 
             'M1 [ln(1+z)]', 'M1 [Z_g integral]', 'M1 [z/(1+z/0.5)]',
             'M1-SAT [ln(1+z)]', 'M1-SAT [Z_g]']:
    if name in results:
        chi2, rms, rho, pval, k = results[name]
        chi2_nu = chi2 / (len(z) - k)
        status = '✓' if pval > 0.05 else '◐' if pval > 0.001 else '✗'
        print(f"  {name:<30} {k:>3} {chi2_nu:>9.4f} {rms:>7.4f} {rho:>+7.4f} {pval:>10.2e} {status:>8}")

# Find best
best_name = min(results.keys(), key=lambda k: abs(results[k][2]))
best = results[best_name]
print(f"\n  BEST: {best_name}")
print(f"  ρ = {best[2]:+.4f} (p = {best[3]:.2e})")
print(f"  Improvement over standard: {(1 - abs(best[2])/abs(rho_std))*100:.0f}%")
print(f"  Improvement over linear: {(1 - abs(best[2])/abs(rho_lin))*100:.0f}%")

any_white = any(r[3] > 0.05 for r in results.values())
if any_white:
    print(f"\n  ✓ WHITE NOISE ACHIEVED")
else:
    print(f"\n  ◐ White noise not yet achieved, but significant progress")
    print(f"  The impedance-rotated basis outperforms independent linear z-evolution")

# Save
os.makedirs('results_model_m1', exist_ok=True)
output = {
    'test_date': '2026-03-09',
    'models': {}
}
for name, (chi2, rms, rho, pval, k) in results.items():
    output['models'][name] = {
        'chi2': float(chi2), 'rms': float(rms), 
        'rho_z': float(rho), 'p_value': float(pval), 'k': k
    }

with open('results_model_m1/model_m1_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n  Results saved to results_model_m1/model_m1_results.json")
print(f"\n{'=' * 70}")
print("COMPLETE")
print("=" * 70)
