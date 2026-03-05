#!/usr/bin/env python3
"""
closure_screen_model.py — Screen-Dominated Compression Model
=============================================================

GPT's insight: the compression isn't continuous diffuse Beer-Lambert.
It's a first-hit percolation process through discrete baryonic screens
(CGM halos, filament interfaces, galaxy group environments).

The transmission is a MIXTURE:
T(z) = E[exp(-τ_diffuse)] × [(1-P_hit(z)) + P_hit(z)×exp(-τ_screen)]

where P_hit(z) = 1 - exp(-N_screen(z)) is the probability of encountering
at least one high-impact screen along the sightline.

KEY PREDICTION: This produces a SHARPER sigmoid than pure Beer-Lambert.
If τ_screen ≥ 1, the transition from "unhit" to "hit" is very sharp.

TEST: Does the screen model reproduce k = 8.0 from physical parameters?

Author: Closure Theory Collaboration  
Date: 2026-03-05
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit, minimize_scalar
from scipy.stats import spearmanr
import json
from pathlib import Path

sigma_T = 6.6524e-25

# Cosmology
H0 = 67.4
Omega_m = 0.315
Omega_L = 0.685
Omega_b = 0.0493
m_p = 1.6726e-24
rho_crit = 1.878e-29
n_H0 = Omega_b * rho_crit * (H0/100)**2 / m_p * 0.76
H0_cgs = 67.4e5 / 3.086e24

def E(z):
    return np.sqrt(Omega_m * (1+z)**3 + Omega_L)

def dSigma_dz(z):
    return n_H0 * (1+z)**2 * 2.998e10 / (H0_cgs * E(z))

def Sigma(z):
    if z <= 0: return 0.0
    result, _ = quad(dSigma_dz, 0, z)
    return result

# ============================================================
# SCREEN MODEL PARAMETERS
# ============================================================

# Number density of "screens" (CGM halos, filament crossings)
# A typical sightline crosses ~1 massive halo CGM per unit redshift
# (from Mg II absorber statistics: dN/dz ≈ 0.5-1.5 for W_r > 0.3 Å)
# This increases with z as structure grows: dN/dz ∝ (1+z)^γ with γ ≈ 1-2

# Screen optical depth for diagnostic content
# A single CGM/filament crossing has N_H ~ 10^19 - 10^21 cm⁻²
# With Γ_Σ ~ 300 × σ_T and q=1:
# τ_screen = Γ × N_H_screen × q²
# For N_H = 10^20: τ = 300 × 6.65e-25 × 1e20 × 1 = 0.02 (too low)
# For N_H = 10^21: τ = 300 × 6.65e-25 × 1e21 × 1 = 0.20 (moderate)
# For N_H = 3×10^21 (CGM): τ = 0.60 (significant)
# For N_H = 10^22 (galaxy disk): τ = 2.0 (saturating)

# The KEY: screens don't need huge τ individually if MULTIPLE screens hit.
# But GPT's model is about the FIRST hit being decisive.

# Let's parameterize and scan:
# dN_screen/dz = n_0 × (1+z)^γ_screen
# P_hit(z) = 1 - exp(-∫₀ᶻ n_0(1+z')^γ dz')
# = 1 - exp(-n_0 × [(1+z)^(γ+1) - 1] / (γ+1))

def N_screens(z, n0, gamma):
    """Expected number of screen crossings from 0 to z"""
    return n0 * ((1+z)**(gamma+1) - 1) / (gamma+1)

def P_hit(z, n0, gamma):
    """Probability of hitting at least one screen"""
    return 1 - np.exp(-N_screens(z, n0, gamma))

def T_screen_model(z, n0, gamma, tau_screen, Gamma_diffuse):
    """Transmission through screen + diffuse model
    
    T = T_diffuse × T_screen_mixture
    T_diffuse = exp(-Gamma_diffuse × Sigma(z) × q²)
    T_screen_mixture = (1-P_hit)×1 + P_hit×exp(-tau_screen)
    """
    # Diffuse component (background Beer-Lambert)
    Sig = Sigma(z) if z > 0 else 0
    T_diff = np.exp(-Gamma_diffuse * Sig * 1.0**2)
    
    # Screen mixture component
    p = P_hit(z, n0, gamma)
    T_mix = (1 - p) * 1.0 + p * np.exp(-tau_screen)
    
    return T_diff * T_mix

def C_screen_model(z, n0, gamma, tau_screen, Gamma_diffuse):
    """Compression = 1 - Transmission"""
    return 1 - T_screen_model(z, n0, gamma, tau_screen, Gamma_diffuse)

# ============================================================
# SCAN SCREEN PARAMETERS FOR k ≈ 8.0
# ============================================================

print("=" * 70)
print("SCREEN MODEL: Can it produce k = 8.0?")
print("=" * 70)

z_grid = np.linspace(0.01, 2.5, 300)

def sigmoid_fit(z, z0, k):
    return 1 / (1 + np.exp(-k * (z - z0)))

print(f"\n{'n₀':>6} {'γ':>5} {'τ_scr':>6} {'Γ_diff':>8} | {'z₀':>6} {'k':>6} {'R²':>8} | Notes")
print("-" * 85)

best_k_diff = 100
best_params = None

for n0 in [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
    for gamma in [0.5, 1.0, 1.5, 2.0]:
        for tau_screen in [0.5, 1.0, 2.0, 3.0, 5.0]:
            for Gamma_ratio in [0, 30, 100]:
                Gamma_diff = Gamma_ratio * sigma_T
                
                # Compute compression curve
                C = np.array([C_screen_model(z, n0, gamma, tau_screen, Gamma_diff) for z in z_grid])
                
                # Skip if doesn't reach 50% by z=2
                if C[-1] < 0.4:
                    continue
                
                # Fit sigmoid
                try:
                    popt, _ = curve_fit(sigmoid_fit, z_grid, C, p0=[0.82, 5.0], maxfev=3000)
                    z0_fit, k_fit = popt
                    
                    C_fit = sigmoid_fit(z_grid, *popt)
                    ss_res = np.sum((C - C_fit)**2)
                    ss_tot = np.sum((C - np.mean(C))**2)
                    r_sq = 1 - ss_res/ss_tot if ss_tot > 0 else 0
                    
                    # Check if z₀ ≈ 0.82 and k ≈ 8
                    if 0.5 < z0_fit < 1.2 and k_fit > 3.0 and r_sq > 0.95:
                        k_diff = abs(k_fit - 8.0)
                        z0_diff = abs(z0_fit - 0.82)
                        score = k_diff + 5 * z0_diff  # weighted distance
                        
                        notes = ""
                        if k_fit > 6:
                            notes += "k>6! "
                        if abs(z0_fit - 0.82) < 0.1:
                            notes += "z₀✓ "
                        
                        if score < best_k_diff:
                            best_k_diff = score
                            best_params = (n0, gamma, tau_screen, Gamma_ratio, z0_fit, k_fit, r_sq)
                        
                        if notes:
                            print(f"{n0:>6.1f} {gamma:>5.1f} {tau_screen:>6.1f} {Gamma_ratio:>8} | {z0_fit:>6.2f} {k_fit:>6.1f} {r_sq:>8.4f} | {notes}")
                except:
                    pass

if best_params:
    n0, gamma, tau_screen, Gamma_ratio, z0_fit, k_fit, r_sq = best_params
    print(f"\n★ BEST MATCH:")
    print(f"  n₀ = {n0} screens/dz, γ = {gamma}, τ_screen = {tau_screen}, Γ_diffuse = {Gamma_ratio}×σ_T")
    print(f"  → z₀ = {z0_fit:.3f}, k = {k_fit:.2f}, R² = {r_sq:.4f}")
    print(f"  Target: z₀ = 0.82, k = 8.0")
    
    if k_fit > 6.0 and abs(z0_fit - 0.82) < 0.15:
        print(f"\n  ★★★ HOT: Screen model reproduces sharp sigmoid!")
        print(f"        k = {k_fit:.1f} vs target 8.0 (Beer-Lambert gave 2.3)")
    elif k_fit > 4.0:
        print(f"\n  ★★ WARM: Significantly sharper than Beer-Lambert (k={k_fit:.1f} vs 2.3)")
    else:
        print(f"\n  ~ Marginal improvement over Beer-Lambert")
else:
    print(f"\n  No solutions with z₀ ≈ 0.82 and k > 3 found")

# ============================================================
# DETAILED ANALYSIS OF BEST MODEL
# ============================================================

if best_params:
    n0, gamma, tau_screen, Gamma_ratio, z0_fit, k_fit, r_sq = best_params
    Gamma_diff = Gamma_ratio * sigma_T
    
    print(f"\n" + "=" * 70)
    print(f"DETAILED ANALYSIS: Best screen model")
    print(f"=" * 70)
    
    # Show P_hit(z) profile
    print(f"\n{'z':>6} {'P_hit':>8} {'N_screens':>10} {'C_screen':>10} {'C_total':>10}")
    print("-" * 50)
    for z_val in [0.1, 0.3, 0.5, 0.7, 0.82, 1.0, 1.2, 1.5, 2.0]:
        p = P_hit(z_val, n0, gamma)
        n = N_screens(z_val, n0, gamma)
        c_scr = p * (1 - np.exp(-tau_screen))  # screen contribution
        c_tot = C_screen_model(z_val, n0, gamma, tau_screen, Gamma_diff)
        print(f"{z_val:>6.2f} {p:>8.3f} {n:>10.2f} {c_scr:>10.3f} {c_tot:>10.3f}")
    
    # Where does the screen dominate over diffuse?
    print(f"\n  Physical meaning:")
    print(f"  n₀ = {n0}: ~{n0:.1f} high-impact structures per unit redshift at z=0")
    print(f"  γ = {gamma}: screen density grows as (1+z)^{gamma:.1f}")
    print(f"  τ_screen = {tau_screen}: each screen destroys {(1-np.exp(-tau_screen))*100:.0f}% of diagnostic info")
    
    z_50_hit = None
    for z_val in np.linspace(0.01, 3.0, 1000):
        if P_hit(z_val, n0, gamma) > 0.5:
            z_50_hit = z_val
            break
    
    if z_50_hit:
        print(f"  z(50% hit probability) = {z_50_hit:.2f}")
        print(f"  At z = {z_50_hit:.2f}: half of sightlines have crossed ≥1 screen")

# ============================================================
# COMPARE: Pure BL vs Screen vs Empirical
# ============================================================

print(f"\n" + "=" * 70)
print(f"COMPARISON: Beer-Lambert vs Screen vs Empirical Sigmoid")
print(f"=" * 70)

# Pure Beer-Lambert (k = 2.3)
Gamma_BL = 311 * sigma_T
C_BL = np.array([1 - np.exp(-Gamma_BL * 1.0 * Sigma(z)) for z in z_grid])

# Screen model (best)
if best_params:
    C_screen = np.array([C_screen_model(z, n0, gamma, tau_screen, Gamma_diff) for z in z_grid])

# Empirical sigmoid (k = 8.0)
C_emp = np.array([sigmoid_fit(z, 0.82, 8.0) for z in z_grid])

# Fit all three
for label, C_data in [("Beer-Lambert", C_BL), ("Screen model", C_screen if best_params else None), ("Empirical", C_emp)]:
    if C_data is None:
        continue
    try:
        popt, _ = curve_fit(sigmoid_fit, z_grid, C_data, p0=[0.8, 5.0], maxfev=3000)
        print(f"  {label:<20}: z₀ = {popt[0]:.3f}, k = {popt[1]:.2f}")
    except:
        print(f"  {label:<20}: fit failed")

# ============================================================
# THE MULTIPLICATIVITY TEST (screen model prediction)
# ============================================================

print(f"\n" + "=" * 70)
print(f"MULTIPLICATIVITY: Screen model prediction")
print(f"=" * 70)

# Screen model generically breaks multiplicativity because the
# "hit/not-hit" mixture doesn't factorize in population averages.
# 
# For a population of sightlines:
# <T(0→z)> ≠ <T(0→z_mid)> × <T(z_mid→z)>
# because the hit/not-hit status is correlated across the path.

if best_params:
    z_test = 1.5
    z_mid = 0.75
    
    # Population average transmission
    T_full = T_screen_model(z_test, n0, gamma, tau_screen, Gamma_diff)
    T_half1 = T_screen_model(z_mid, n0, gamma, tau_screen, Gamma_diff)
    
    # For the second half, we need the CONDITIONAL transmission
    # given that we've already traveled from 0 to z_mid.
    # This is where the mixture breaks: if you haven't been hit yet,
    # the second half has a different P_hit than unconditional.
    
    # Unconditional T(z_mid→z_test):
    # Would be T_screen_model(z_test-z_mid, ...) if independent
    # But the screen crossings are NOT independent segments
    
    # Simple product test:
    T_product = T_half1 * T_screen_model(z_test, n0, gamma, tau_screen, Gamma_diff) / T_screen_model(z_test, n0, gamma, tau_screen, Gamma_diff)
    # Actually, let me compute it properly:
    # The "second segment" transmission (from z_mid to z_test):
    # For screens: N_screens(z_mid→z_test) = N_screens(z_test) - N_screens(z_mid)
    N_scr_full = N_screens(z_test, n0, gamma)
    N_scr_half1 = N_screens(z_mid, n0, gamma)
    N_scr_half2 = N_scr_full - N_scr_half1
    
    P_hit_half2 = 1 - np.exp(-N_scr_half2)
    T_screen_half2 = (1-P_hit_half2) + P_hit_half2 * np.exp(-tau_screen)
    
    # Diffuse component IS multiplicative
    Sig_full = Sigma(z_test)
    Sig_mid = Sigma(z_mid)
    T_diff_half2 = np.exp(-Gamma_diff * (Sig_full - Sig_mid))
    
    T_half2_independent = T_diff_half2 * T_screen_half2
    T_product_independent = T_half1 * T_half2_independent
    
    # But the TRUE population-averaged transmission is:
    # <T(0→z)> = <T_diff(0→z) × T_screen(0→z)>
    # And <T(0→z_mid) × T(z_mid→z)> ≠ <T(0→z)> when there are correlations
    
    print(f"  z_test = {z_test}, z_mid = {z_mid}")
    print(f"  T(0→{z_test}):                    {T_full:.6f}")
    print(f"  T(0→{z_mid}) × T({z_mid}→{z_test}):  {T_product_independent:.6f}")
    print(f"  Ratio:                           {T_product_independent/T_full:.4f}")
    print(f"  Deviation from multiplicativity:  {abs(T_product_independent/T_full - 1)*100:.1f}%")
    
    # Compare to pure Beer-Lambert (always multiplicative)
    T_BL_full = np.exp(-Gamma_BL * Sig_full)
    T_BL_half1 = np.exp(-Gamma_BL * Sig_mid)
    T_BL_half2 = np.exp(-Gamma_BL * (Sig_full - Sig_mid))
    print(f"\n  Beer-Lambert ratio:              {T_BL_half1 * T_BL_half2 / T_BL_full:.4f}")
    print(f"  (always 1.0000 by construction)")
    
    if abs(T_product_independent/T_full - 1) > 0.05:
        print(f"\n  ★ Screen model BREAKS multiplicativity by {abs(T_product_independent/T_full - 1)*100:.0f}%")
        print(f"    This is a testable prediction: split sightlines should show non-multiplicative compression")

# ============================================================
# CMB COVARIANCE: GPT's A_s + ω_m direction
# ============================================================

print(f"\n" + "=" * 70)
print(f"CMB COVARIANCE: GPT's ℓ-dependent compression direction")
print(f"=" * 70)

print(f"""
  GPT's mechanism for anti-correlated H₀↓ and σ₈↑:
  
  1. ℓ-dependent compression suppresses high-ℓ more than low-ℓ
  2. Fitter response:
     a) Increases A_s × exp(-2τ) to restore lost high-ℓ power
        → A_s up → σ₈ UP (σ₈ ∝ A_s^0.5 × growth)
     b) Increases ω_m = Ω_m h² to fix peak structure at fixed θ*
        → h down → H₀ DOWN
  3. Both happen simultaneously → anti-correlated tensions
  
  Key insight: A_s and ω_m both shift UP, but:
     - A_s up → σ₈ up (amplitude effect dominates)
     - ω_m up at fixed θ* → h down → H₀ down (geometry)
  
  TT vs EE discriminant:
  - Temperature and polarization are sourced differently
  - Post-recombination scattering should produce DIFFERENT
    suppression patterns in TT vs EE
  - Planck TT gives A_L = 1.22 ± 0.08
  - Planck EE gives A_L = 0.97 ± 0.26
  - The DIFFERENCE (TT > EE) is consistent with compression
    affecting temperature more than polarization!
""")

# Planck A_L split
AL_TT = 1.22
AL_TT_err = 0.08
AL_EE = 0.97
AL_EE_err = 0.26

delta_AL = AL_TT - AL_EE
delta_AL_err = np.sqrt(AL_TT_err**2 + AL_EE_err**2)

print(f"  A_L(TT) - A_L(EE) = {delta_AL:.2f} ± {delta_AL_err:.2f}")
print(f"  Significance: {delta_AL/delta_AL_err:.1f}σ")
print(f"  Direction: TT > EE → temperature MORE compressed than polarization")

if delta_AL > 0:
    print(f"\n  ★ CONSISTENT with compression affecting T more than E")
    print(f"    (Resonant scattering affects intensity more than polarization)")
else:
    print(f"\n  ✗ Wrong direction for compression model")

# ============================================================
# SAVE RESULTS
# ============================================================

results_dir = Path("results_screen_model")
results_dir.mkdir(exist_ok=True)

results = {
    'best_params': {
        'n0': best_params[0] if best_params else None,
        'gamma': best_params[1] if best_params else None,
        'tau_screen': best_params[2] if best_params else None,
        'Gamma_diffuse': best_params[3] if best_params else None,
        'z0_fit': best_params[4] if best_params else None,
        'k_fit': best_params[5] if best_params else None,
        'r_sq': best_params[6] if best_params else None,
    },
    'AL_TT_EE_split': {
        'AL_TT': AL_TT,
        'AL_EE': AL_EE,
        'delta': float(delta_AL),
        'sigma': float(delta_AL/delta_AL_err),
    }
}

with open(results_dir / "screen_model_results.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {results_dir}/")
print("\n" + "=" * 70)
print("END")
print("=" * 70)
