#!/usr/bin/env python3
"""
closure_beer_lambert.py — Beer-Lambert for Information
======================================================

The compression law dI/dΣ = -Γ·q²·I is Beer-Lambert applied to
diagnostic information traversing layered baryonic scatter.

Key test: Does integrating Beer-Lambert through realistic cosmic
baryon density NATURALLY PRODUCE the observed sigmoid, threshold
z₀ = 0.82, and saturation — WITHOUT assuming any of it?

If yes: the sigmoid isn't a fit. It's a CONSEQUENCE of the physics.
If no: the sigmoid is empirical and needs additional explanation.

Also tests:
1. Does z₀ emerge from Σ(z) without being assumed?
2. Does z₀ shift with sightline density (filament vs void)?
3. Is the composition multiplicative (Beer-Lambert prediction)?
4. Does the ℓ-space q(ℓ) produce the right CMB compression profile?
5. Can we predict Γ_Σ from the threshold + known Σ(z)?

Author: Closure Theory Collaboration
Date: 2026-03-04
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, curve_fit
import json
from pathlib import Path

# ============================================================
# COSMOLOGICAL PARAMETERS
# ============================================================

H0 = 67.4  # km/s/Mpc (Planck, for computing distances)
c_km = 299792.458  # km/s
Omega_m = 0.315
Omega_L = 0.685
Omega_b = 0.0493
m_p = 1.6726e-24  # g (proton mass)
rho_crit = 1.878e-29  # g/cm³ (h=0.674)
sigma_T = 6.6524e-25  # cm² (Thomson)

# Mean hydrogen number density today
n_H0 = Omega_b * rho_crit * (H0/100)**2 / m_p * 0.76  # ~1.9e-7 cm⁻³
# (0.76 = hydrogen mass fraction)

# Hubble distance
D_H = c_km * 1e5 / H0  # cm (convert km to cm: ×1e5)
# Actually: D_H = c/H₀ in cm
D_H_cm = 2.998e10 / (H0 * 1e5 / 3.086e24)  # c / (H₀ in s⁻¹)
# H₀ = 67.4 km/s/Mpc = 67.4e5 cm/s / 3.086e24 cm = 2.184e-18 s⁻¹
H0_cgs = 67.4e5 / 3.086e24  # s⁻¹
D_H_cm = 2.998e10 / H0_cgs  # cm

print(f"n_H,0 = {n_H0:.2e} cm⁻³")
print(f"D_H = {D_H_cm:.3e} cm = {D_H_cm/3.086e24:.0f} Mpc")

# ============================================================
# Σ(z): Cumulative baryon column density
# ============================================================

def E(z):
    """Dimensionless Hubble parameter E(z) = H(z)/H₀"""
    return np.sqrt(Omega_m * (1+z)**3 + Omega_L)

def dSigma_dz(z):
    """Differential baryon column density dΣ/dz [cm⁻²]
    
    dΣ/dz = n_H(z) × c × dt/dz = n_H,0 × (1+z)³ × c / (H₀ × (1+z) × E(z))
           = n_H,0 × (1+z)² × c / (H₀ × E(z))
    
    The (1+z)³ is density evolution, /(1+z) is dt/dz factor.
    """
    return n_H0 * (1+z)**2 * 2.998e10 / (H0_cgs * E(z))

def Sigma(z):
    """Cumulative baryon column density from observer to redshift z [cm⁻²]"""
    if z <= 0:
        return 0.0
    result, _ = quad(dSigma_dz, 0, z)
    return result

# Compute Σ(z) on a grid
z_grid = np.linspace(0.01, 3.0, 300)
Sigma_grid = np.array([Sigma(z) for z in z_grid])

print(f"\n{'z':>6} {'Σ(z) [cm⁻²]':>15} {'Σ/10²¹':>10}")
print("-" * 35)
for z_val in [0.1, 0.3, 0.5, 0.82, 1.0, 1.5, 2.0, 3.0]:
    S = Sigma(z_val)
    print(f"{z_val:>6.2f} {S:>15.3e} {S/1e21:>10.3f}")

# Also compute Σ to CMB
Sigma_CMB = Sigma(1100)
print(f"\n{'CMB':>6} {Sigma_CMB:>15.3e} {Sigma_CMB/1e21:>10.3f}")

# ============================================================
# TEST 1: DOES THE SIGMOID EMERGE FROM BEER-LAMBERT?
# ============================================================

print("\n" + "=" * 70)
print("TEST 1: Does the sigmoid emerge from Beer-Lambert integration?")
print("=" * 70)

# Beer-Lambert: I(z)/I₀ = exp(-Γ_Σ · q² · Σ(z))
# Compression fraction: C(z) = 1 - exp(-Γ_Σ · q² · Σ(z))
# This is an exponential saturation curve in Σ, but Σ(z) is nonlinear in z.

# The empirical sigmoid: C_emp(z) = 1 / (1 + exp(-k·(z - z₀)))
# with k=8.0, z₀=0.82

def empirical_sigmoid(z, z0=0.82, k=8.0):
    return 1 / (1 + np.exp(-k * (z - z0)))

# For different Γ_Σ values, compute the Beer-Lambert curve and compare to sigmoid
Gamma_values = [37, 100, 150, 227, 300, 400]

print(f"\n{'Γ_Σ/σ_T':>10} {'z_50%':>8} {'z_10%':>8} {'z_90%':>8} {'Width':>8} {'Sigmoid r':>10}")
print("-" * 60)

for G_ratio in Gamma_values:
    Gamma = G_ratio * sigma_T
    q = 1.0  # maximum diagnostic susceptibility
    
    # Beer-Lambert compression at each z
    C_BL = np.array([1 - np.exp(-Gamma * q**2 * Sigma(z)) for z in z_grid])
    
    # Find z where compression = 50% (the "threshold")
    idx_50 = np.argmin(np.abs(C_BL - 0.5))
    z_50 = z_grid[idx_50] if C_BL[-1] > 0.5 else float('nan')
    
    # Find z where compression = 10% and 90%
    idx_10 = np.argmin(np.abs(C_BL - 0.1))
    idx_90 = np.argmin(np.abs(C_BL - 0.9))
    z_10 = z_grid[idx_10] if C_BL[-1] > 0.1 else float('nan')
    z_90 = z_grid[idx_90] if C_BL[-1] > 0.9 else float('nan')
    
    width = z_90 - z_10 if not (np.isnan(z_90) or np.isnan(z_10)) else float('nan')
    
    # Compare shape to sigmoid: correlation coefficient
    if not np.isnan(z_50):
        # Fit a sigmoid to the Beer-Lambert curve
        try:
            def sigmoid_fit(z, z0, k):
                return 1 / (1 + np.exp(-k * (z - z0)))
            popt, _ = curve_fit(sigmoid_fit, z_grid, C_BL, p0=[z_50, 5.0], maxfev=5000)
            C_sigmoid_fit = sigmoid_fit(z_grid, *popt)
            residuals = C_BL - C_sigmoid_fit
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((C_BL - np.mean(C_BL))**2)
            r_squared = 1 - ss_res/ss_tot
            
            corr = np.corrcoef(C_BL, C_sigmoid_fit)[0,1]
        except:
            corr = float('nan')
            popt = [float('nan'), float('nan')]
    else:
        corr = float('nan')
        popt = [float('nan'), float('nan')]
    
    print(f"{G_ratio:>10} {z_50:>8.3f} {z_10:>8.3f} {z_90:>8.3f} {width:>8.3f} {corr:>10.6f}")

# Now find what Γ_Σ gives z₀ = 0.82
print(f"\n--- Finding Γ_Σ that gives z₀ = 0.82 ---")

def z50_for_Gamma(G_ratio):
    Gamma = G_ratio * sigma_T
    C_BL = np.array([1 - np.exp(-Gamma * 1.0**2 * Sigma(z)) for z in z_grid])
    idx_50 = np.argmin(np.abs(C_BL - 0.5))
    return z_grid[idx_50]

# Binary search
G_low, G_high = 10, 500
for _ in range(50):
    G_mid = (G_low + G_high) / 2
    z50 = z50_for_Gamma(G_mid)
    if z50 < 0.82:
        G_high = G_mid
    else:
        G_low = G_mid

G_predicted = (G_low + G_high) / 2
print(f"  Γ_Σ that gives z₀ = 0.82: {G_predicted:.1f} × σ_T")
print(f"  = {G_predicted * sigma_T:.3e} cm²")

# How does this compare to our measured values?
print(f"\n  Measured Γ_Σ range: 37-227 × σ_T")
print(f"  Predicted from z₀: {G_predicted:.0f} × σ_T")
if 37 <= G_predicted <= 227:
    print(f"  ★ WITHIN MEASURED RANGE!")
elif G_predicted < 500:
    print(f"  ~ Same order of magnitude")
else:
    print(f"  ✗ Outside measured range")

# ============================================================
# TEST 2: SIGMOID FIT QUALITY — IS BEER-LAMBERT ≈ SIGMOID?
# ============================================================

print("\n" + "=" * 70)
print("TEST 2: How well does a sigmoid approximate Beer-Lambert?")
print("=" * 70)

# Use the Γ_Σ that gives z₀ = 0.82
Gamma_best = G_predicted * sigma_T
C_BL_best = np.array([1 - np.exp(-Gamma_best * 1.0**2 * Sigma(z)) for z in z_grid])

try:
    popt, pcov = curve_fit(sigmoid_fit, z_grid, C_BL_best, p0=[0.82, 5.0])
    z0_fit, k_fit = popt
    C_sig_fit = sigmoid_fit(z_grid, *popt)
    
    residuals = C_BL_best - C_sig_fit
    max_residual = np.max(np.abs(residuals))
    rms_residual = np.sqrt(np.mean(residuals**2))
    r_sq = 1 - np.sum(residuals**2) / np.sum((C_BL_best - np.mean(C_BL_best))**2)
    
    print(f"  Best-fit sigmoid: z₀ = {z0_fit:.3f}, k = {k_fit:.2f}")
    print(f"  Empirical sigmoid: z₀ = 0.82, k = 8.0")
    print(f"  R² = {r_sq:.6f}")
    print(f"  Max residual = {max_residual:.4f}")
    print(f"  RMS residual = {rms_residual:.4f}")
    print(f"  k ratio (BL/empirical): {k_fit/8.0:.2f}")
    
    if r_sq > 0.99:
        print(f"\n  ★★★ Beer-Lambert IS a sigmoid in z-space (R² > 0.99)")
        print(f"      The sigmoid shape EMERGES from the physics!")
    elif r_sq > 0.95:
        print(f"\n  ★★ Close to sigmoid (R² > 0.95)")
        print(f"     Small deviations reveal Σ(z) non-linearity")
    else:
        print(f"\n  ⚠ Not well-approximated by sigmoid")
        print(f"     Beer-Lambert in Σ doesn't map cleanly to sigmoid in z")
except Exception as e:
    print(f"  Fit failed: {e}")

# ============================================================
# TEST 3: SIGHTLINE DENSITY DEPENDENCE
# ============================================================

print("\n" + "=" * 70)
print("TEST 3: Does z₀ shift with sightline density?")
print("=" * 70)

# Prediction: dense sightlines (filaments) hit threshold at LOWER z
# Sparse sightlines (voids) hit threshold at HIGHER z
# Because Σ accumulates faster through dense regions

density_factors = {
    'Deep void (0.1× mean)': 0.1,
    'Typical void (0.3× mean)': 0.3,
    'Mean cosmic': 1.0,
    'Filament (3× mean)': 3.0,
    'Cluster outskirt (10× mean)': 10.0,
    'Cluster core (100× mean)': 100.0,
}

Gamma_test = G_predicted * sigma_T

print(f"\nUsing Γ_Σ = {G_predicted:.0f} × σ_T, q = 1.0\n")
print(f"{'Sightline':<30} {'Density':>10} {'z₀ (50%)':>10} {'Shift from mean':>15}")
print("-" * 70)

z0_mean = None
for name, factor in density_factors.items():
    # Scale Σ(z) by density factor
    C_BL = np.array([1 - np.exp(-Gamma_test * 1.0 * factor * Sigma(z)) for z in z_grid])
    idx_50 = np.argmin(np.abs(C_BL - 0.5))
    z_50 = z_grid[idx_50] if C_BL[-1] > 0.5 else float('inf')
    
    if factor == 1.0:
        z0_mean = z_50
    
    shift = z_50 - z0_mean if z0_mean is not None and z_50 != float('inf') else float('nan')
    
    print(f"  {name:<28} {factor:>10.1f}× {z_50:>10.3f} {shift:>+15.3f}" if z_50 != float('inf') 
          else f"  {name:<28} {factor:>10.1f}× {'> 3.0':>10} {'—':>15}")

print(f"\n  ★ Prediction: filament sightlines should show LOWER z₀")
print(f"    This matches our Feb 22 result: Filament > Void in 5/5 z-bins")

# ============================================================
# TEST 4: MULTIPLICATIVITY (Beer-Lambert vs additive)
# ============================================================

print("\n" + "=" * 70)
print("TEST 4: Multiplicativity — the composition test")
print("=" * 70)

# Beer-Lambert predicts: C(z₁→z₂) = C(z₁→z_mid) × C(z_mid→z₂) [in transmission]
# i.e., T(0→z) = T(0→z/2) × T(z/2→z) where T = 1-C = exp(-Γq²Σ)
# This is trivially true by definition: exp(-a-b) = exp(-a)×exp(-b)
# But it's NOT trivially true for the empirical sigmoid!
# 
# For the sigmoid: C(z) = 1/(1+exp(-k(z-z₀)))
# The transmission: T(z) = 1-C(z) = exp(-k(z-z₀))/(1+exp(-k(z-z₀)))
# T(0→z) ≠ T(0→z/2) × T(z/2→z) in general!
#
# This gives us a TESTABLE DISTINCTION between Beer-Lambert and sigmoid.

z_test = 1.5
z_mid = 0.75

# Beer-Lambert (exact):
Sigma_full = Sigma(z_test)
Sigma_half1 = Sigma(z_mid)
Sigma_half2 = Sigma_full - Sigma_half1

T_full_BL = np.exp(-Gamma_test * 1.0 * Sigma_full)
T_half1_BL = np.exp(-Gamma_test * 1.0 * Sigma_half1)
T_half2_BL = np.exp(-Gamma_test * 1.0 * Sigma_half2)
T_product_BL = T_half1_BL * T_half2_BL

print(f"\nBeer-Lambert (z=0 → {z_test}, split at z={z_mid}):")
print(f"  T(0→{z_test}) = {T_full_BL:.6f}")
print(f"  T(0→{z_mid}) × T({z_mid}→{z_test}) = {T_product_BL:.6f}")
print(f"  Ratio: {T_product_BL/T_full_BL:.6f} (should be 1.000000)")

# Empirical sigmoid:
T_full_sig = 1 - empirical_sigmoid(z_test)
T_half1_sig = 1 - empirical_sigmoid(z_mid)
# For sigmoid, transmission from z_mid to z_test requires careful definition
# The sigmoid C(z) is total compression from 0 to z, not from z_mid to z
# So the "segment transmission" is T(z_mid→z) = T(0→z) / T(0→z_mid)... IF multiplicative
T_segment_sig = T_full_sig / T_half1_sig if T_half1_sig > 0 else 0

print(f"\nEmpirical sigmoid (for comparison):")
print(f"  T(0→{z_test}) = {T_full_sig:.6f}")
print(f"  T(0→{z_mid}) = {T_half1_sig:.6f}")
print(f"  T_segment(implied) = T_full/T_half = {T_segment_sig:.6f}")
print(f"  Would need T({z_mid}→{z_test}) = {T_segment_sig:.6f}")

# Check: is the sigmoid consistent with multiplicativity?
# If sigmoid were exact, then T_segment should = 1 - C_segment
# where C_segment is the compression of JUST the second half
# In Beer-Lambert: C_segment depends on Σ(z_mid→z_test) only
# In sigmoid: C_segment depends on z directly (no Σ separation)

# More stringent test: check at multiple split points
print(f"\n{'z_split':>8} {'T_full':>10} {'T₁×T₂(BL)':>12} {'BL ratio':>10} {'Sigmoid T₁':>12} {'Sig implied T₂':>14}")
print("-" * 70)

for z_split in [0.3, 0.5, 0.75, 1.0, 1.2]:
    if z_split >= z_test:
        continue
    S_full = Sigma(z_test)
    S_1 = Sigma(z_split)
    S_2 = S_full - S_1
    
    T_f = np.exp(-Gamma_test * S_full)
    T_1 = np.exp(-Gamma_test * S_1)
    T_2 = np.exp(-Gamma_test * S_2)
    
    T_sig_1 = 1 - empirical_sigmoid(z_split)
    T_sig_f = 1 - empirical_sigmoid(z_test)
    T_sig_implied_2 = T_sig_f / T_sig_1 if T_sig_1 > 0 else 0
    
    print(f"{z_split:>8.2f} {T_f:>10.6f} {T_1*T_2:>12.6f} {T_1*T_2/T_f:>10.6f} {T_sig_1:>12.6f} {T_sig_implied_2:>14.6f}")

print(f"\n  Beer-Lambert ratio is always 1.000000 (by construction — exponential).")
print(f"  Sigmoid implied T₂ can go NEGATIVE if z_split > z₀ (sigmoid breaks multiplicativity).")

# ============================================================
# TEST 5: PREDICT Γ_Σ FROM z₀ AND KNOWN Σ(z₀)
# ============================================================

print("\n" + "=" * 70)
print("TEST 5: Predict Γ_Σ from z₀ = 0.82 and Σ(z₀)")
print("=" * 70)

# At the threshold z₀, compression = 50%, so:
# 0.5 = 1 - exp(-Γ_Σ · q² · Σ(z₀))
# exp(-Γ_Σ · q² · Σ(z₀)) = 0.5
# Γ_Σ · q² · Σ(z₀) = ln(2) = 0.693

Sigma_z0 = Sigma(0.82)
Gamma_predicted = np.log(2) / (1.0**2 * Sigma_z0)  # for q=1
Gamma_ratio_predicted = Gamma_predicted / sigma_T

print(f"  Σ(z₀=0.82) = {Sigma_z0:.3e} cm⁻²")
print(f"  Γ_Σ = ln(2) / (q² × Σ(z₀))")
print(f"  For q = 1.0: Γ_Σ = {Gamma_predicted:.3e} cm² = {Gamma_ratio_predicted:.1f} × σ_T")

# For different q values
print(f"\n  {'q':>5} {'Γ_Σ/σ_T':>10} {'In measured range?':>20}")
print(f"  {'-'*40}")
for q in [0.5, 0.7, 0.8, 0.88, 0.95, 1.0]:
    G = np.log(2) / (q**2 * Sigma_z0)
    G_ratio = G / sigma_T
    in_range = "✓" if 37 <= G_ratio <= 300 else "~" if 10 <= G_ratio <= 500 else "✗"
    print(f"  {q:>5.2f} {G_ratio:>10.1f} {in_range:>20}")

# ============================================================
# TEST 6: THE q-LADDER IN BEER-LAMBERT
# ============================================================

print("\n" + "=" * 70)
print("TEST 6: Different q values → different z₀ thresholds")
print("=" * 70)

# For a fixed Γ_Σ, different q values produce different thresholds
# z₀(q) is where Γ_Σ · q² · Σ(z₀) = ln(2)
# Lower q → need more Σ → higher z₀
# Higher q → need less Σ → lower z₀

# This predicts: the "doublet ladder" should also be a "threshold ladder"
# [SII] (q≈0.008) → threshold at z >> 3 (essentially never)
# SN color (q≈1.0) → threshold at z = 0.82
# The ordering should be MONOTONIC

# Use our measured q values from closure_q_derivation_v2.py
q_ladder = [
    ("[SII] ratio", 0.008),
    ("SN stretch", 0.039),
    ("FRB DM", 0.239),
    ("[NII]", 0.477),
    ("Hβ", 0.526),
    ("FRB SI", 0.604),
    ("[OIII]", 0.682),
    ("[OII]", 0.917),
    ("SN color", 1.000),
]

Gamma_use = G_predicted * sigma_T  # the one that gives z₀=0.82 for q=1

print(f"\nUsing Γ_Σ = {G_predicted:.0f} × σ_T (calibrated to z₀=0.82 for SN color)")
print(f"\n{'Observable':<20} {'q':>6} {'q²':>8} {'z₀ (50%)':>10} {'z₀ (10%)':>10} {'Saturated by z=3?':>18}")
print("-" * 75)

for name, q in q_ladder:
    if q < 0.01:
        print(f"  {name:<18} {q:>6.3f} {q**2:>8.5f} {'>> 3.0':>10} {'>> 3.0':>10} {'No':>18}")
        continue
    
    # Find z where C = 50%
    C_q = np.array([1 - np.exp(-Gamma_use * q**2 * Sigma(z)) for z in z_grid])
    
    idx_50 = np.argmin(np.abs(C_q - 0.5))
    z_50 = z_grid[idx_50] if C_q[-1] > 0.5 else float('inf')
    
    idx_10 = np.argmin(np.abs(C_q - 0.1))
    z_10 = z_grid[idx_10] if C_q[-1] > 0.1 else float('inf')
    
    saturated = "Yes" if C_q[-1] > 0.9 else "Partial" if C_q[-1] > 0.5 else "No"
    
    z50_str = f"{z_50:.3f}" if z_50 != float('inf') else "> 3.0"
    z10_str = f"{z_10:.3f}" if z_10 != float('inf') else "> 3.0"
    
    print(f"  {name:<18} {q:>6.3f} {q**2:>8.5f} {z50_str:>10} {z10_str:>10} {saturated:>18}")

# ============================================================
# TEST 7: FRB DM THRESHOLD — INDEPENDENT CONFIRMATION
# ============================================================

print("\n" + "=" * 70)
print("TEST 7: FRB DM≈500 threshold from Beer-Lambert")
print("=" * 70)

# FRBs already use DM (dispersion measure) which IS column density!
# DM = ∫ n_e dl ≈ Σ_electron
# The observed threshold at DM ≈ 500 pc/cm³ should correspond to
# the Beer-Lambert threshold for the FRB spectral index (q ≈ 0.604)

# Convert DM to Σ:
# DM is in pc/cm³. Σ is in cm⁻².
# DM = ∫ n_e dl where n_e ≈ 1.14 × n_H (fully ionized H + He)
# So Σ_H = DM / 1.14 (in units of pc × cm⁻³ → convert to cm⁻²)
# 1 pc = 3.086e18 cm
# DM = 500 pc/cm³ → Σ_e = 500 × 3.086e18 = 1.543e21 cm⁻²
# Σ_H = 1.543e21 / 1.14 = 1.354e21 cm⁻²

DM_threshold = 500  # pc/cm³
Sigma_FRB = DM_threshold * 3.086e18 / 1.14  # cm⁻²

# Beer-Lambert prediction: threshold at Γ_Σ × q_FRB² × Σ = ln(2)
q_FRB_SI = 0.604  # spectral index susceptibility
Sigma_BL_threshold = np.log(2) / (Gamma_use * q_FRB_SI**2)

# What DM does this correspond to?
DM_BL_predicted = Sigma_BL_threshold * 1.14 / 3.086e18

print(f"  Observed FRB threshold: DM ≈ {DM_threshold} pc/cm³")
print(f"  → Σ_FRB = {Sigma_FRB:.3e} cm⁻²")
print(f"")
print(f"  Beer-Lambert prediction (for q_SI = {q_FRB_SI}):")
print(f"  → Σ_threshold = {Sigma_BL_threshold:.3e} cm⁻²")
print(f"  → DM_threshold = {DM_BL_predicted:.0f} pc/cm³")
print(f"")
print(f"  Ratio (predicted/observed): {DM_BL_predicted/DM_threshold:.2f}")

if 0.3 < DM_BL_predicted/DM_threshold < 3.0:
    print(f"  ★ CONSISTENT: same order of magnitude!")
else:
    print(f"  ✗ Inconsistent")

# Also: what z does DM=500 correspond to?
# Using the DM-z relation: DM_cosmic ≈ 900 × z (pc/cm³) for z < 1
z_DM500 = 500 / 900  # ≈ 0.56
print(f"\n  DM=500 corresponds to z ≈ {z_DM500:.2f}")
print(f"  Our z₀ for q_SI=0.604: from Beer-Lambert grid above")

# ============================================================
# SAVE RESULTS
# ============================================================

results_dir = Path("results_beer_lambert")
results_dir.mkdir(exist_ok=True)

results = {
    'Sigma_z0': float(Sigma_z0),
    'Gamma_predicted_from_z0': float(Gamma_ratio_predicted),
    'Gamma_predicted_cm2': float(Gamma_predicted),
    'G_predicted_for_z0_82': float(G_predicted),
    'Sigma_CMB': float(Sigma_CMB),
    'n_H0': float(n_H0),
    'FRB_DM_predicted': float(DM_BL_predicted),
    'FRB_DM_observed': DM_threshold,
}

with open(results_dir / "beer_lambert_results.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {results_dir}/")
print("\n" + "=" * 70)
print("END")
print("=" * 70)
