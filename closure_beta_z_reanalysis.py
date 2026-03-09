#!/usr/bin/env python3
"""
CLOSURE THEORY — β(z) REANALYSIS
==================================

THE TEST:
Replace the constant β in Tripp standardization with a z-dependent β(z).
If the Hubble tension decreases, metric impedance is real.

This is framework-INDEPENDENT. Pure empirical observation:
"Does allowing β to vary with z change the inferred H₀?"

If yes: the color coefficient absorbs a z-dependent systematic.
If no: color is well-behaved and impedance doesn't affect SN cosmology.

We test multiple β(z) models:
1. Constant β (standard Tripp — baseline)
2. Linear β(z) = β₀ + β₁×z  
3. Step function β(z) = β_lo (z<z_cut) or β_hi (z≥z_cut)
4. Power law β(z) = β₀ × (1+z)^γ
5. Impedance-motivated: β(z) = β_dust + β_imp/(1 + z/z₀)

For each model, we fit for H₀ (via M_B) and check:
- Does H₀ change?
- Does χ² improve?
- Does the mass step change?
- Do Hubble residuals still correlate with color?

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

z_hd = np.array([float(r[col['zHD']]) for r in rows])
z_cmb = np.array([float(r[col['zCMB']]) for r in rows])
m_b = np.array([float(r[col['mB']]) for r in rows])
m_b_err = np.array([float(r[col['mBERR']]) for r in rows])
x1 = np.array([float(r[col['x1']]) for r in rows])
c = np.array([float(r[col['c']]) for r in rows])
host_mass = np.array([float(r[col['HOST_LOGMASS']]) for r in rows])

# Quality cuts
mask_q = (z_hd > 0.01) & (z_hd < 2.5) & (np.abs(x1) < 5) & (np.abs(c) < 0.5) & (m_b_err < 1.0) & (host_mass > 0)

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

def luminosity_distance(z_arr, H0, Om=0.3):
    """Compute luminosity distance for flat LCDM"""
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        integral, _ = quad(lambda zp: 1.0/np.sqrt(Om*(1+zp)**3 + (1-Om)), 0, zi)
        dl[i] = C_LIGHT * (1+zi) * integral / H0
    return dl

def distance_modulus(z_arr, H0, Om=0.3):
    dl = luminosity_distance(z_arr, H0, Om)
    return 5 * np.log10(dl) + 25

# Pre-compute distance moduli on a grid for speed
print("Building distance modulus grid...")
z_grid = np.linspace(0.01, 2.5, 500)
mu_grid_70 = distance_modulus(z_grid, H0=70.0)

def mu_interp(z_arr, H0):
    """Fast distance modulus using interpolation + H0 scaling"""
    # μ(H0) = μ(70) - 5*log10(H0/70)
    mu_at_70 = np.interp(z_arr, z_grid, mu_grid_70)
    return mu_at_70 - 5*np.log10(H0/70.0)

print("Grid ready.")


# ============================================================
# β(z) MODELS
# ============================================================

def beta_constant(z, params):
    """Standard Tripp: β = const"""
    return np.full_like(z, params['beta0'])

def beta_linear(z, params):
    """β(z) = β₀ + β₁ × z"""
    return params['beta0'] + params['beta1'] * z

def beta_step(z, params):
    """β = β_lo for z < z_cut, β_hi for z ≥ z_cut"""
    beta = np.where(z < params['z_cut'], params['beta_lo'], params['beta_hi'])
    return beta

def beta_power(z, params):
    """β(z) = β₀ × (1+z)^γ"""
    return params['beta0'] * (1 + z)**params['gamma']

def beta_impedance(z, params):
    """β(z) = β_dust + β_imp / (1 + z/z₀)
    At low z: β ≈ β_dust + β_imp (full)
    At high z: β ≈ β_dust (impedance saturates, β compensates less)"""
    return params['beta_dust'] + params['beta_imp'] / (1 + z/params['z0'])


# ============================================================
# FITTING FRAMEWORK
# ============================================================

def fit_model(z_d, mb_d, x1_d, c_d, err_d, hm_d, beta_func, param_names, x0, 
              fit_H0=True, fit_mass_step=True):
    """
    Fit SN standardization with arbitrary β(z) model.
    Always fits: M_B, α
    Optionally fits: H₀, mass_step (γ)
    Plus whatever β(z) parameters the model needs.
    """
    mass_thresh = 10.0
    mask_hi = hm_d >= mass_thresh
    
    def neg_loglik(params_vec):
        idx = 0
        M_B = params_vec[idx]; idx += 1
        alpha = params_vec[idx]; idx += 1
        
        if fit_H0:
            H0 = params_vec[idx]; idx += 1
        else:
            H0 = 70.0
            
        if fit_mass_step:
            gamma = params_vec[idx]; idx += 1
        else:
            gamma = 0.0
        
        # β model parameters
        beta_params = {}
        for name in param_names:
            beta_params[name] = params_vec[idx]; idx += 1
        
        # Compute β(z) for each SN
        beta_z = beta_func(z_d, beta_params)
        
        # Standard candle equation
        mu_obs = mb_d + alpha * x1_d - beta_z * c_d - M_B
        mu_obs[mask_hi] -= gamma
        
        # Model distance modulus
        mu_mod = mu_interp(z_d, H0)
        
        # χ²
        chi2 = np.sum(((mu_obs - mu_mod) / err_d)**2)
        return chi2
    
    # Count degrees of freedom
    n_params = 2 + int(fit_H0) + int(fit_mass_step) + len(param_names)
    
    res = scipy_minimize(neg_loglik, x0=x0, method='Nelder-Mead',
                        options={'maxiter': 50000, 'fatol': 1e-10, 'xatol': 1e-10})
    
    # Extract results
    idx = 0
    results = {}
    results['M_B'] = res.x[idx]; idx += 1
    results['alpha'] = res.x[idx]; idx += 1
    if fit_H0:
        results['H0'] = res.x[idx]; idx += 1
    else:
        results['H0'] = 70.0
    if fit_mass_step:
        results['gamma'] = res.x[idx]; idx += 1
    else:
        results['gamma'] = 0.0
    
    results['beta_params'] = {}
    for name in param_names:
        results['beta_params'][name] = res.x[idx]; idx += 1
    
    results['chi2'] = res.fun
    results['chi2_nu'] = res.fun / (len(z_d) - n_params)
    results['n_params'] = n_params
    results['AIC'] = res.fun + 2 * n_params
    results['BIC'] = res.fun + n_params * np.log(len(z_d))
    
    # Compute H₀ from M_B (using absolute magnitude relation)
    # H₀ = 10^((M_B + 19.23 + 25) / 5) × c / d_L_ref
    # More directly: ΔM_B = -5 × log10(H₀/H₀_ref)
    # We calibrate against the standard Tripp M_B = -19.253 → H₀ = 73.04
    if not fit_H0:
        delta_MB = results['M_B'] - (-19.253)
        results['H0_from_MB'] = 73.04 * 10**(delta_MB / 5)
    
    # Residuals
    beta_z = beta_func(z_d, results['beta_params'])
    mu_obs = mb_d + results['alpha'] * x1_d - beta_z * c_d - results['M_B']
    mu_obs[mask_hi] -= results['gamma']
    mu_mod = mu_interp(z_d, results['H0'])
    results['residuals'] = mu_obs - mu_mod
    results['beta_z_values'] = beta_z
    
    return results


# ============================================================
# RUN ALL MODELS
# ============================================================

print(f"\n{'=' * 70}")
print("FITTING ALL β(z) MODELS")
print("=" * 70)

models = {}

# Model 1: Constant β (baseline)
print("\n  Fitting Model 1: Constant β...")
r1 = fit_model(z, mb, x1_q, c_q, mb_err, hm,
               beta_constant, ['beta0'],
               x0=[-19.25, 0.15, 73.0, -0.05, 3.0],
               fit_H0=True, fit_mass_step=True)
models['constant'] = r1
print(f"    β = {r1['beta_params']['beta0']:.3f}, H₀ = {r1['H0']:.2f}, γ = {r1['gamma']:+.4f}, χ²/ν = {r1['chi2_nu']:.4f}")

# Model 2: Linear β(z)
print("  Fitting Model 2: Linear β(z)...")
r2 = fit_model(z, mb, x1_q, c_q, mb_err, hm,
               beta_linear, ['beta0', 'beta1'],
               x0=[-19.25, 0.15, 73.0, -0.05, 3.0, -1.0],
               fit_H0=True, fit_mass_step=True)
models['linear'] = r2
print(f"    β₀ = {r2['beta_params']['beta0']:.3f}, β₁ = {r2['beta_params']['beta1']:+.3f}")
print(f"    H₀ = {r2['H0']:.2f}, γ = {r2['gamma']:+.4f}, χ²/ν = {r2['chi2_nu']:.4f}")

# Model 3: Step function β(z)
print("  Fitting Model 3: Step function β(z)...")
r3 = fit_model(z, mb, x1_q, c_q, mb_err, hm,
               beta_step, ['beta_lo', 'beta_hi', 'z_cut'],
               x0=[-19.25, 0.15, 73.0, -0.05, 3.2, 2.5, 0.10],
               fit_H0=True, fit_mass_step=True)
models['step'] = r3
print(f"    β_lo = {r3['beta_params']['beta_lo']:.3f}, β_hi = {r3['beta_params']['beta_hi']:.3f}, z_cut = {r3['beta_params']['z_cut']:.3f}")
print(f"    H₀ = {r3['H0']:.2f}, γ = {r3['gamma']:+.4f}, χ²/ν = {r3['chi2_nu']:.4f}")

# Model 4: Power law β(z)
print("  Fitting Model 4: Power law β(z)...")
r4 = fit_model(z, mb, x1_q, c_q, mb_err, hm,
               beta_power, ['beta0', 'gamma'],
               x0=[-19.25, 0.15, 73.0, -0.05, 3.5, -0.5],
               fit_H0=True, fit_mass_step=True)
models['power'] = r4
print(f"    β₀ = {r4['beta_params']['beta0']:.3f}, γ = {r4['beta_params']['gamma']:+.3f}")
print(f"    β(z=0) = {r4['beta_params']['beta0']:.3f}, β(z=1) = {r4['beta_params']['beta0'] * 2**r4['beta_params']['gamma']:.3f}")
print(f"    H₀ = {r4['H0']:.2f}, mass_step = {r4['gamma']:+.4f}, χ²/ν = {r4['chi2_nu']:.4f}")

# Model 5: Impedance-motivated β(z)
print("  Fitting Model 5: Impedance-motivated β(z)...")
r5 = fit_model(z, mb, x1_q, c_q, mb_err, hm,
               beta_impedance, ['beta_dust', 'beta_imp', 'z0'],
               x0=[-19.25, 0.15, 73.0, -0.05, 2.0, 1.5, 0.5],
               fit_H0=True, fit_mass_step=True)
models['impedance'] = r5
print(f"    β_dust = {r5['beta_params']['beta_dust']:.3f}, β_imp = {r5['beta_params']['beta_imp']:.3f}, z₀ = {r5['beta_params']['z0']:.3f}")
print(f"    β(z=0) = {r5['beta_params']['beta_dust'] + r5['beta_params']['beta_imp']:.3f}")
print(f"    β(z=1) = {r5['beta_params']['beta_dust'] + r5['beta_params']['beta_imp']/(1+1/r5['beta_params']['z0']):.3f}")
print(f"    H₀ = {r5['H0']:.2f}, mass_step = {r5['gamma']:+.4f}, χ²/ν = {r5['chi2_nu']:.4f}")


# ============================================================
# COMPARISON TABLE
# ============================================================

print(f"\n\n{'=' * 70}")
print("MODEL COMPARISON")
print("=" * 70)

print(f"\n  {'Model':<25} {'k':>3} {'χ²/ν':>10} {'AIC':>10} {'BIC':>10} {'H₀':>7} {'γ':>8} {'ΔH₀':>7}")
print(f"  {'-'*85}")

H0_base = models['constant']['H0']

for name, label in [('constant', 'Constant β'),
                     ('linear', 'Linear β(z)'),
                     ('step', 'Step β(z)'),
                     ('power', 'Power law β(z)'),
                     ('impedance', 'Impedance β(z)')]:
    r = models[name]
    delta_h0 = r['H0'] - H0_base
    print(f"  {label:<25} {r['n_params']:>3} {r['chi2_nu']:>10.4f} {r['AIC']:>10.1f} {r['BIC']:>10.1f} "
          f"{r['H0']:>7.2f} {r['gamma']:>+8.4f} {delta_h0:>+7.2f}")

# Δχ² test for nested models (linear vs constant)
dchi2_linear = models['constant']['chi2'] - models['linear']['chi2']
dof_linear = 1  # one extra parameter
from scipy.stats import chi2 as chi2_dist
p_linear = 1 - chi2_dist.cdf(dchi2_linear, dof_linear)

print(f"\n  STATISTICAL TESTS:")
print(f"    Δχ² (constant → linear): {dchi2_linear:.2f} (p = {p_linear:.4f})")
print(f"    {'SIGNIFICANT' if p_linear < 0.05 else 'not significant'} improvement with β(z)")

dchi2_power = models['constant']['chi2'] - models['power']['chi2']
p_power = 1 - chi2_dist.cdf(dchi2_power, 1)
print(f"    Δχ² (constant → power):  {dchi2_power:.2f} (p = {p_power:.4f})")
print(f"    {'SIGNIFICANT' if p_power < 0.05 else 'not significant'} improvement")

dchi2_imp = models['constant']['chi2'] - models['impedance']['chi2']
p_imp = 1 - chi2_dist.cdf(dchi2_imp, 2)
print(f"    Δχ² (constant → impedance): {dchi2_imp:.2f} (p = {p_imp:.4f})")
print(f"    {'SIGNIFICANT' if p_imp < 0.05 else 'not significant'} improvement")


# ============================================================
# H₀ SHIFT ANALYSIS
# ============================================================

print(f"\n\n{'=' * 70}")
print("H₀ SHIFT FROM β(z) CORRECTION")
print("=" * 70)

print(f"\n  Standard (constant β): H₀ = {models['constant']['H0']:.2f}")
print(f"  Linear β(z):           H₀ = {models['linear']['H0']:.2f}  (Δ = {models['linear']['H0'] - models['constant']['H0']:+.2f})")
print(f"  Power law β(z):        H₀ = {models['power']['H0']:.2f}  (Δ = {models['power']['H0'] - models['constant']['H0']:+.2f})")
print(f"  Impedance β(z):        H₀ = {models['impedance']['H0']:.2f}  (Δ = {models['impedance']['H0'] - models['constant']['H0']:+.2f})")

print(f"\n  Planck CMB:            H₀ = 67.36")
print(f"  Tension with constant β: {models['constant']['H0'] - 67.36:.2f} km/s/Mpc")

best_beta_z = min(['linear', 'power', 'impedance'], key=lambda k: models[k]['chi2_nu'])
print(f"  Tension with best β(z) [{best_beta_z}]: {models[best_beta_z]['H0'] - 67.36:.2f} km/s/Mpc")
print(f"  Tension reduction: {(models['constant']['H0'] - 67.36) - (models[best_beta_z]['H0'] - 67.36):.2f} km/s/Mpc")
print(f"  Fractional reduction: {((models['constant']['H0'] - models[best_beta_z]['H0']) / (models['constant']['H0'] - 67.36))*100:.1f}%")


# ============================================================
# MASS STEP UNDER β(z)
# ============================================================

print(f"\n\n{'=' * 70}")
print("MASS STEP UNDER β(z) MODELS")
print("=" * 70)

print(f"\n  {'Model':<25} {'γ (mass step)':>15} {'Δ from constant':>18}")
print(f"  {'-'*60}")

for name, label in [('constant', 'Constant β'),
                     ('linear', 'Linear β(z)'),
                     ('power', 'Power law β(z)'),
                     ('impedance', 'Impedance β(z)')]:
    r = models[name]
    delta_gamma = abs(r['gamma']) - abs(models['constant']['gamma'])
    print(f"  {label:<25} {r['gamma']:>+15.4f} {delta_gamma:>+18.4f}")

print(f"\n  If β(z) absorbs impedance properly:")
print(f"  - The mass step should DECREASE (less environment bias leaks into γ)")
print(f"  - H₀ should DECREASE (toward Planck)")
print(f"  - Both happening simultaneously = impedance is the mechanism")


# ============================================================
# RESIDUAL-COLOR CORRELATION UNDER β(z)
# ============================================================

print(f"\n\n{'=' * 70}")
print("RESIDUAL-COLOR CORRELATION AFTER β(z) CORRECTION")
print("=" * 70)

z_bins = [(0.01, 0.10), (0.10, 0.30), (0.30, 0.60), (0.60, 1.50)]

print(f"\n  {'z-bin':>15} {'const β':>10} {'linear':>10} {'power':>10} {'impedance':>10}")
print(f"  {'-'*58}")

for z_lo, z_hi in z_bins:
    mask_bin = (z >= z_lo) & (z < z_hi)
    n = np.sum(mask_bin)
    if n < 15:
        continue
    
    rhos = []
    for name in ['constant', 'linear', 'power', 'impedance']:
        rho, _ = spearmanr(models[name]['residuals'][mask_bin], c_q[mask_bin])
        rhos.append(rho)
    
    print(f"  [{z_lo:.2f}, {z_hi:.2f}]  {rhos[0]:>+10.3f} {rhos[1]:>+10.3f} {rhos[2]:>+10.3f} {rhos[3]:>+10.3f}")

print(f"\n  EXPECTATION: β(z) should REDUCE the residual-color correlation")
print(f"  especially at high z where constant β is most wrong")


# ============================================================
# β(z) PROFILE
# ============================================================

print(f"\n\n{'=' * 70}")
print("β(z) PROFILE — WHAT THE DATA SAYS")
print("=" * 70)

z_profile = np.array([0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 1.00, 1.50])

print(f"\n  {'z':>6} {'constant':>10} {'linear':>10} {'power':>10} {'impedance':>10}")
print(f"  {'-'*48}")

for zi in z_profile:
    b_const = beta_constant(np.array([zi]), r1['beta_params'])[0]
    b_lin = beta_linear(np.array([zi]), r2['beta_params'])[0]
    b_pow = beta_power(np.array([zi]), r4['beta_params'])[0]
    b_imp = beta_impedance(np.array([zi]), r5['beta_params'])[0]
    print(f"  {zi:>6.2f} {b_const:>10.3f} {b_lin:>10.3f} {b_pow:>10.3f} {b_imp:>10.3f}")


# ============================================================
# THE THIRD LAW? — CHECKING FOR MISSING PHYSICS
# ============================================================

print(f"\n\n{'=' * 70}")
print("HUNTING THE THIRD LAW — WHAT'S MISSING?")
print("=" * 70)

# Humza's intuition: there's something small we're not accounting for.
# Let's look at what the residuals tell us.

# After best β(z) correction, what patterns remain?
best = models[best_beta_z]
resid = best['residuals']

# 1. Residuals vs redshift — any remaining trend?
rho_rz, p_rz = spearmanr(z, resid)
print(f"\n  Residuals vs z (after {best_beta_z} β(z)):")
print(f"    Spearman ρ = {rho_rz:+.4f} (p = {p_rz:.4f})")

# 2. Residuals vs stretch — should be zero
rho_rx1, p_rx1 = spearmanr(x1_q, resid)
print(f"  Residuals vs stretch (x1):")
print(f"    Spearman ρ = {rho_rx1:+.4f} (p = {p_rx1:.4f})")

# 3. Residuals vs host mass — any remaining structure?
rho_rm, p_rm = spearmanr(hm, resid)
print(f"  Residuals vs host mass:")
print(f"    Spearman ρ = {rho_rm:+.4f} (p = {p_rm:.4f})")

# 4. Residuals vs color SQUARED — nonlinear color effects?
rho_rc2, p_rc2 = spearmanr(c_q**2, resid)
print(f"  Residuals vs c²:")
print(f"    Spearman ρ = {rho_rc2:+.4f} (p = {p_rc2:.4f})")

# 5. Residuals vs stretch × color interaction
rho_rxc, p_rxc = spearmanr(x1_q * c_q, resid)
print(f"  Residuals vs x1×c (interaction):")
print(f"    Spearman ρ = {rho_rxc:+.4f} (p = {p_rxc:.4f})")

# 6. Residuals vs stretch × z — does stretch evolve too?
rho_rx1z, p_rx1z = spearmanr(x1_q * z, resid)
print(f"  Residuals vs x1×z (stretch evolution):")
print(f"    Spearman ρ = {rho_rx1z:+.4f} (p = {p_rx1z:.4f})")

# 7. Residuals vs host mass × z — mass step evolves?
rho_rmz, p_rmz = spearmanr(hm * z, resid)
print(f"  Residuals vs mass×z (mass step evolution):")
print(f"    Spearman ρ = {rho_rmz:+.4f} (p = {p_rmz:.4f})")

# 8. Color × host mass interaction — the missing link?
rho_rcm, p_rcm = spearmanr(c_q * hm, resid)
print(f"  Residuals vs c×mass (color-mass coupling):")
print(f"    Spearman ρ = {rho_rcm:+.4f} (p = {p_rcm:.4f})")

# Find the most significant remaining pattern
remaining = [
    ('z', rho_rz, p_rz),
    ('x1', rho_rx1, p_rx1),
    ('host_mass', rho_rm, p_rm),
    ('c²', rho_rc2, p_rc2),
    ('x1×c', rho_rxc, p_rxc),
    ('x1×z', rho_rx1z, p_rx1z),
    ('mass×z', rho_rmz, p_rmz),
    ('c×mass', rho_rcm, p_rcm),
]

remaining.sort(key=lambda x: x[2])  # sort by p-value

print(f"\n  MOST SIGNIFICANT REMAINING PATTERNS (after {best_beta_z} β(z)):")
print(f"  {'Variable':<15} {'ρ':>8} {'p-value':>12} {'Significant?':>14}")
print(f"  {'-'*52}")
for name, rho, p in remaining[:5]:
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
    print(f"  {name:<15} {rho:>+8.4f} {p:>12.2e} {sig:>14}")

print(f"\n  IF significant patterns remain after β(z):")
print(f"  → There may be a THIRD correction term beyond α(z) and β(z)")
print(f"  → Humza's 'third law' intuition")
print(f"  → Could be: α(z), γ(z), or an interaction term")


# ============================================================
# ALPHA(z) CHECK — DOES STRETCH ALSO EVOLVE?
# ============================================================

print(f"\n\n{'=' * 70}")
print("BONUS: DOES α ALSO NEED z-DEPENDENCE?")
print("=" * 70)

# Fit α(z) = α₀ + α₁ × z alongside constant β
def fit_alpha_z(z_d, mb_d, x1_d, c_d, err_d, hm_d):
    mass_thresh = 10.0
    mask_hi = hm_d >= mass_thresh
    
    def chi2(params):
        M_B, alpha0, alpha1, beta, gamma, H0 = params
        alpha_z = alpha0 + alpha1 * z_d
        mu_obs = mb_d + alpha_z * x1_d - beta * c_d - M_B
        mu_obs[mask_hi] -= gamma
        mu_mod = mu_interp(z_d, H0)
        return np.sum(((mu_obs - mu_mod) / err_d)**2)
    
    res = scipy_minimize(chi2, x0=[-19.25, 0.15, -0.05, 3.0, -0.05, 73.0],
                        method='Nelder-Mead', options={'maxiter': 50000})
    
    M_B, alpha0, alpha1, beta, gamma, H0 = res.x
    return {'M_B': M_B, 'alpha0': alpha0, 'alpha1': alpha1, 'beta': beta,
            'gamma': gamma, 'H0': H0, 'chi2': res.fun, 
            'chi2_nu': res.fun / (len(z_d) - 6)}

r_az = fit_alpha_z(z, mb, x1_q, c_q, mb_err, hm)
print(f"\n  α(z) = {r_az['alpha0']:.4f} + {r_az['alpha1']:+.4f} × z")
print(f"  α(z=0) = {r_az['alpha0']:.4f}")
print(f"  α(z=1) = {r_az['alpha0'] + r_az['alpha1']:.4f}")
print(f"  H₀ = {r_az['H0']:.2f}")
print(f"  χ²/ν = {r_az['chi2_nu']:.4f}")
print(f"  Δχ² vs constant α: {models['constant']['chi2'] - r_az['chi2']:.2f}")

# Both α(z) and β(z)?
def fit_both_z(z_d, mb_d, x1_d, c_d, err_d, hm_d):
    mass_thresh = 10.0
    mask_hi = hm_d >= mass_thresh
    
    def chi2(params):
        M_B, alpha0, alpha1, beta0, beta1, gamma, H0 = params
        alpha_z = alpha0 + alpha1 * z_d
        beta_z = beta0 + beta1 * z_d
        mu_obs = mb_d + alpha_z * x1_d - beta_z * c_d - M_B
        mu_obs[mask_hi] -= gamma
        mu_mod = mu_interp(z_d, H0)
        return np.sum(((mu_obs - mu_mod) / err_d)**2)
    
    res = scipy_minimize(chi2, x0=[-19.25, 0.15, -0.05, 3.0, -1.0, -0.05, 73.0],
                        method='Nelder-Mead', options={'maxiter': 50000})
    
    M_B, alpha0, alpha1, beta0, beta1, gamma, H0 = res.x
    return {'M_B': M_B, 'alpha0': alpha0, 'alpha1': alpha1, 
            'beta0': beta0, 'beta1': beta1,
            'gamma': gamma, 'H0': H0, 'chi2': res.fun,
            'chi2_nu': res.fun / (len(z_d) - 7)}

r_both = fit_both_z(z, mb, x1_q, c_q, mb_err, hm)
print(f"\n  COMBINED: α(z) + β(z)")
print(f"  α(z) = {r_both['alpha0']:.4f} + {r_both['alpha1']:+.4f} × z")
print(f"  β(z) = {r_both['beta0']:.4f} + {r_both['beta1']:+.4f} × z")
print(f"  H₀ = {r_both['H0']:.2f}")
print(f"  γ (mass step) = {r_both['gamma']:+.4f}")
print(f"  χ²/ν = {r_both['chi2_nu']:.4f}")
print(f"  Δχ² vs constant: {models['constant']['chi2'] - r_both['chi2']:.2f}")

# Final H₀ comparison
print(f"\n\n{'=' * 70}")
print("FINAL: H₀ UNDER ALL MODELS")
print("=" * 70)

print(f"\n  {'Model':<35} {'H₀':>7} {'Tension':>9} {'χ²/ν':>10}")
print(f"  {'-'*65}")
all_results = [
    ('Constant β (standard)', models['constant']['H0'], models['constant']['chi2_nu']),
    ('Linear β(z)', models['linear']['H0'], models['linear']['chi2_nu']),
    ('Power law β(z)', models['power']['H0'], models['power']['chi2_nu']),
    ('Impedance β(z)', models['impedance']['H0'], models['impedance']['chi2_nu']),
    ('α(z) only', r_az['H0'], r_az['chi2_nu']),
    ('α(z) + β(z)', r_both['H0'], r_both['chi2_nu']),
    ('Stretch-only (no color)', 70.50, 42.99),
    ('Planck CMB', 67.36, None),
]

for name, h0, chi2_nu in all_results:
    tension = h0 - 67.36
    chi_str = f"{chi2_nu:.4f}" if chi2_nu else "—"
    print(f"  {name:<35} {h0:>7.2f} {tension:>+9.2f} {chi_str:>10}")


# Save
output = {
    'test_date': '2026-03-09',
    'models': {},
    'alpha_z': {'alpha0': r_az['alpha0'], 'alpha1': r_az['alpha1'], 'H0': r_az['H0']},
    'both_z': {'alpha0': r_both['alpha0'], 'alpha1': r_both['alpha1'],
               'beta0': r_both['beta0'], 'beta1': r_both['beta1'], 'H0': r_both['H0']},
}

for name in models:
    r = models[name]
    output['models'][name] = {
        'H0': float(r['H0']),
        'gamma': float(r['gamma']),
        'chi2_nu': float(r['chi2_nu']),
        'AIC': float(r['AIC']),
        'BIC': float(r['BIC']),
        'beta_params': {k: float(v) for k, v in r['beta_params'].items()},
    }

os.makedirs('results_beta_z', exist_ok=True)
with open('results_beta_z/beta_z_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n  Results saved to results_beta_z/beta_z_results.json")
print(f"\n{'=' * 70}")
print("COMPLETE")
print("=" * 70)
