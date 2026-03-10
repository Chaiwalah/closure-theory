#!/usr/bin/env python3
"""
closure_cmb_propagation.py — The Correlation Propagation Test
=============================================================

The five-observable test was wrong because it treated CMB parameters
as independent. They're not. They're extracted by the same fitter from
the same power spectrum. Compression doesn't hit them independently —
it propagates through the parameter covariance structure.

This is the CMB version of the doublet ladder.

The approach:
1. Take Planck's published parameter covariance matrix
2. Model compression as a perturbation to peak heights (diagnostic content)
3. Propagate through the fitter's degeneracy structure
4. Check if ONE compression amplitude produces the JOINT observed shifts
   in H₀, σ₈, A_L, ω_b, n_s simultaneously

The key insight (from Feb 22): the operator acts on CORRELATIONS, not signals.
The CMB parameters are correlated. One perturbation → correlated shifts.

Author: Closure Theory Collaboration
Date: 2026-03-04
"""

import numpy as np
import json
from pathlib import Path

# ============================================================
# PLANCK 2018 PARAMETER VALUES AND COVARIANCE
# ============================================================

# Planck 2018 best-fit values (TT,TE,EE+lowE+lensing, Table 2 of 1807.06209)
planck_params = {
    'H0': 67.36,           # km/s/Mpc
    'omega_b': 0.02237,    # Ω_b h²
    'omega_c': 0.1200,     # Ω_c h²
    'sigma8': 0.8111,
    'n_s': 0.9649,
    'tau': 0.0544,         # reionization optical depth
    'A_L': 1.180,          # lensing anomaly (when floated)
}

# Planck 2018 ΛCDM parameter uncertainties (1σ, from Table 2)
planck_sigma = {
    'H0': 0.54,
    'omega_b': 0.00015,
    'omega_c': 0.0012,
    'sigma8': 0.0060,
    'n_s': 0.0042,
    'tau': 0.0073,
    'A_L': 0.065,
}

# Planck published parameter correlation coefficients
# (from Planck 2018 MCMC chains, approximate from published degeneracy plots)
# The KEY correlations:
#   H₀-ω_b: +0.60 (more baryons → higher H₀ through r_d)
#   H₀-ω_c: -0.55 (more CDM → lower H₀)  
#   H₀-σ₈: +0.40 (higher H₀ → higher σ₈ through growth)
#   H₀-n_s: +0.60 (higher n_s → higher H₀ through peak positions)
#   σ₈-ω_c: +0.75 (more CDM → higher σ₈)
#   A_L-ω_b: +0.30 (A_L partially absorbs ω_b shifts)
#   H₀-A_L: -0.20 (when A_L floated, anti-correlated with H₀)
#   σ₈-A_L: -0.45 (A_L smoothing anti-correlated with amplitude)

# Parameter order: H0, omega_b, omega_c, sigma8, n_s, tau, A_L
param_names = ['H0', 'omega_b', 'omega_c', 'sigma8', 'n_s', 'tau', 'A_L']
n_params = len(param_names)

# Correlation matrix (approximate from published Planck chains)
corr = np.array([
    # H0    ω_b    ω_c    σ₈     n_s    τ      A_L
    [1.00,  0.60, -0.55,  0.40,  0.60,  0.05, -0.20],  # H0
    [0.60,  1.00, -0.30,  0.10,  0.65,  0.08,  0.30],  # ω_b
    [-0.55,-0.30,  1.00,  0.75, -0.40, -0.05,  0.10],  # ω_c
    [0.40,  0.10,  0.75,  1.00,  0.10,  0.60, -0.45],  # σ₈
    [0.60,  0.65, -0.40,  0.10,  1.00,  0.15, -0.10],  # n_s
    [0.05,  0.08, -0.05,  0.60,  0.15,  1.00, -0.15],  # τ
    [-0.20, 0.30,  0.10, -0.45, -0.10, -0.15,  1.00],  # A_L
])

# Construct covariance matrix
sigmas = np.array([planck_sigma[p] for p in param_names])
cov = np.outer(sigmas, sigmas) * corr

print("=" * 70)
print("CMB CORRELATION PROPAGATION TEST")
print("The operator acts on CORRELATIONS, not signals.")
print("=" * 70)

# ============================================================
# THE COMPRESSION MODEL
# ============================================================
# 
# Diagnostic compression suppresses CMB peak heights.
# In parameter space, this is a SPECIFIC DIRECTION of perturbation:
#
# Peak heights are primarily controlled by:
#   - ω_b (baryon loading → odd/even peak ratio) — HIGH diagnostic content
#   - τ (reionization → overall suppression) — moderate diagnostic content
#   - A_s (primordial amplitude → overall normalization) — mixed
#
# Peak positions are controlled by:
#   - θ_* = r_s/D_A (angular scale → purely geometric) — LOCKED
#   - ω_c (CDM → equality scale → peak spacing) — moderate
#
# So compression primarily perturbs ω_b (and τ, A_s).
# The fitter then propagates this through the covariance to ALL parameters.
#
# The COMPRESSION DIRECTION in parameter space:
# We perturb ω_b downward (compression makes peaks look less baryon-loaded)
# and let the covariance propagate it.

print("\n--- The Compression Direction ---")
print("Peak height suppression → fitter interprets as lower ω_b")
print("Then covariance propagates to all parameters.\n")

# ============================================================
# WHAT SHIFTS ARE OBSERVED?
# ============================================================

# "True" values (from geometric/local measurements):
true_values = {
    'H0': 73.04,      # SH0ES (local, geometric anchors)
    'omega_b': None,   # We don't have an independent measurement
    'omega_c': None,   # We don't have an independent measurement  
    'sigma8': 0.766,   # Weak lensing (geometric)
    'n_s': None,       # No independent measurement
    'tau': None,       # No independent measurement
    'A_L': 1.000,      # Should be 1.0 by definition (no anomaly)
}

# Observed shifts (Planck - True):
# These are the shifts that compression must explain
observed_shifts = {
    'H0': 67.36 - 73.04,      # = -5.68
    'sigma8': 0.8111 - 0.766,  # = +0.045
    'A_L': 1.180 - 1.000,      # = +0.180
}

print("Observed shifts (Planck_inferred - True_geometric):")
for p, shift in observed_shifts.items():
    print(f"  Δ{p} = {shift:+.3f} ({shift/planck_sigma[p]:+.1f}σ)")

# ============================================================
# THE PROPAGATION TEST
# ============================================================

# If we perturb ω_b by δ (in units of σ), the covariance propagates it:
# Δp_i = Σ_j (C_ij / C_jj) × δ_j
# Where δ is the perturbation in the ω_b direction.
#
# In terms of correlation: if ω_b shifts by N sigmas,
# parameter i shifts by r(i, ω_b) × (σ_i/σ_ωb) × Δω_b... 
# Actually simpler: Δp_i/σ_i = r(i, ω_b) × Δω_b/σ_ωb

# But it's not JUST ω_b. The compression also directly affects:
# - Peak heights → ω_b direction (primary)
# - Damping tail → τ and A_L direction (secondary)
# - Power spectrum tilt → n_s direction (tertiary)

# Let's model the compression as a VECTOR in parameter space
# that points in the "peak suppression" direction.
# This direction is a weighted combination:
#   Primary: ω_b (weight 0.7) — most diagnostic
#   Secondary: τ (weight 0.2) — affects normalization
#   Tertiary: A_L (weight -0.1) — compensates smoothing

# The compression amplitude α (in sigma units) is what we solve for.

print("\n" + "=" * 70)
print("PROPAGATION: Single compression → all parameter shifts")
print("=" * 70)

# Method: directly use the correlation structure.
# If compression acts primarily through ω_b:
# Δω_b = -α × σ_ωb  (negative: compression looks like less baryons)
# Then each parameter shifts by:
# Δp_i = r(i, ω_b) × σ_i × (-α)

# What α is needed to get the observed H₀ shift?
# ΔH₀ = r(H0, ω_b) × σ_H0 × (-α) = -5.68
# α = 5.68 / (0.60 × 0.54) = 17.5

# That's 17.5σ in ω_b — way too much. 
# This means ω_b isn't the only compression channel.

# BETTER MODEL: Compression acts as a perturbation to the DAMPING TAIL
# (high-ℓ CMB power), which the fitter decomposes into multiple parameters.
# The fitter's response to "less small-scale power" is a COMBINATION of:
#   - Lower ω_b (to reduce baryon loading)
#   - Lower H₀ (degenerate with ω_b through r_d)
#   - Higher A_L (to smooth peaks via lensing)
#   - Lower σ₈ (less power = less clustering)... wait, σ₈ goes UP in Planck

# Actually, the fitter response to damping tail suppression is subtle.
# Let me think about this more carefully.

# The CMB power spectrum at high ℓ: C_ℓ ∝ A_s × exp(-2τ) × D(ℓ) × T²(ℓ)
# where D(ℓ) is diffusion damping and T(ℓ) is the transfer function.
# If compression suppresses high-ℓ power, the fitter can absorb it via:
#   1. Increasing silk damping (decrease ω_b → decrease c_s → more diffusion)
#   2. Increasing τ (more reionization → more suppression)
#   3. Increasing A_L (more lensing smoothing)
#   4. Tilting n_s red (less small-scale power)

# Each absorption channel has a different "cost" in χ².
# The fitter finds the minimum-χ² combination.

# Published Planck results give us the ACTUAL shifts when A_L is floated:
# Without A_L: H₀ = 67.36, ω_b = 0.02237
# With A_L: H₀ = 67.9 ± 0.7, ω_b = 0.02249 ± 0.00016, A_L = 1.18 ± 0.065
# The shifts when A_L is floated:
# ΔH₀ = +0.54, Δω_b = +0.00012, ΔA_L = +0.18
# This tells us how the fitter REDISTRIBUTES when given the A_L degree of freedom.

print("\n--- Published Planck shifts when A_L is floated ---")
print("  ΔH₀ = +0.54 km/s/Mpc (H₀ goes UP slightly)")
print("  Δω_b = +0.00012 (baryons go UP slightly)")
print("  ΔA_L = +0.18 (lensing anomaly appears)")
print("  This is the fitter redistributing the compression signal.\n")

# ============================================================
# THE DOUBLET LADDER ANALOG
# ============================================================

# In the emission line doublet ladder:
#   q = 0 (locked) → zero degradation
#   q = 1 (diagnostic) → maximum degradation
#   The ORDERING is monotonic, r = -0.975
#
# For CMB parameters, each parameter has an effective q:
#   Peak POSITIONS (θ*) → q ≈ 0 (geometric, locked)
#   Peak HEIGHTS (ω_b, A_s) → q ≈ 0.9 (thermodynamic, diagnostic)
#   Damping tail (Silk damping) → q ≈ 0.8 (thermal diffusion)
#   Lensing (C_ℓ^φφ) → q ≈ 0.1 (gravitational, geometric)
#
# When compression hits the high-q channels, the fitter shifts
# the low-q channels to compensate. This creates the ILLUSION
# of different Γ_Σ for different parameters.

print("=" * 70)
print("THE CMB DOUBLET LADDER")
print("(Same law, applied to CMB parameter susceptibilities)")
print("=" * 70)

# CMB observables and their diagnostic susceptibility
cmb_ladder = [
    ("θ* (angular scale)",        0.00, "geometric ruler → LOCKED",           "unchanged"),
    ("BAO scale (D_V/r_d)",       0.05, "geometric × r_d anchor",             "r_d biased"),
    ("Lensing power (C_ℓ^φφ)",    0.10, "gravitational → mostly geometric",   "A_L = 1.18"),
    ("Peak spacing (ℓ_n+1-ℓ_n)",  0.15, "geometric + baryon loading",         "slight shift"),
    ("Matter-radiation eq (k_eq)", 0.25, "density ratio → mixed",             "H₀ ≈ 70-71"),
    ("n_s (spectral tilt)",        0.40, "inflation + transfer function",      "slight red tilt"),
    ("τ (reionization)",           0.50, "ionization state → diagnostic",      "degeneracy"),
    ("σ₈ (power amplitude)",       0.60, "normalization × growth → mixed",     "S₈ tension"),
    ("Peak height ratios (R_n)",   0.85, "baryon loading → thermodynamic",     "ω_b biased"),
    ("Silk damping (ℓ_D)",         0.90, "thermal diffusion → diagnostic",     "damping biased"),
    ("r_d (sound horizon)",        0.88, "thermodynamic integral → diagnostic","H₀ = 67"),
]

print(f"\n{'Observable':<35} {'q':>5} {'Nature':<40} {'Observed Effect'}")
print("-" * 120)
for name, q, nature, effect in cmb_ladder:
    marker = "◆" if q > 0.5 else "◇" if q > 0.2 else "○"
    print(f"{marker} {name:<33} {q:>5.2f} {nature:<40} {effect}")

print(f"\n○ = locked/geometric (q < 0.2)")
print(f"◇ = mixed (0.2 < q < 0.5)")
print(f"◆ = diagnostic (q > 0.5)")

# ============================================================
# PREDICTION: THE CMB COMPRESSION SIGNATURE
# ============================================================

print("\n" + "=" * 70)
print("PREDICTION: CMB COMPRESSION SIGNATURE")
print("=" * 70)

# If compression acts on the CMB with our measured Γ_Σ,
# the signature in ℓ-space should be:
#   - Low ℓ (large scales): UNAFFECTED (superhorizon, geometric)
#   - ℓ ≈ 200-800 (first 3 peaks): PARTIALLY affected (peak heights, not positions)
#   - ℓ > 1000 (damping tail): MAXIMALLY affected (Silk damping is diagnostic)
#
# This ℓ-dependent compression is EXACTLY what the A_L anomaly is:
# the fitter sees "too smooth" at high ℓ relative to low ℓ.

# Quantitative prediction:
# The compression at multipole ℓ scales as:
# ε(ℓ) = Γ_Σ × Σ_* × q²(ℓ)
# where q(ℓ) transitions from ~0 at low ℓ to ~0.9 at high ℓ
# with a transition around ℓ ~ 500-800 (where diffusion damping kicks in)

sigma_T = 6.6524e-25
Gamma = 227 * sigma_T  # Tripp-measured value
Sigma = 3e21  # cm⁻²

ell_values = [2, 50, 200, 500, 800, 1000, 1500, 2000, 2500]
print(f"\nUsing Γ_Σ = 227 × σ_T, Σ_* = 3×10²¹ cm⁻²")
print(f"\n{'ℓ':>6} {'q(ℓ)':>8} {'ε(ℓ)':>10} {'Peak suppression':>18} {'Primary content'}")
print("-" * 80)

for ell in ell_values:
    # q(ℓ) model: sigmoid transition from geometric to diagnostic
    # Low ℓ: geometric (ISW, SW plateau)
    # High ℓ: diagnostic (Silk damping, baryon loading)
    q_ell = 0.9 / (1 + np.exp(-(ell - 600) / 200))  # sigmoid, midpoint ℓ=600
    
    epsilon = Gamma * Sigma * q_ell**2
    suppression = 1 - np.exp(-epsilon)
    
    if ell <= 50:
        content = "SW plateau (geometric)"
    elif ell <= 300:
        content = "1st peak (mixed)"
    elif ell <= 600:
        content = "2nd-3rd peaks (diagnostic↑)"
    elif ell <= 1200:
        content = "Damping onset (diagnostic)"
    else:
        content = "Silk tail (max diagnostic)"
    
    print(f"{ell:>6} {q_ell:>8.3f} {epsilon:>10.5f} {suppression*100:>15.2f}% {content}")

# ============================================================
# THE KEY TEST: DOES ONE ε EXPLAIN THE FITTER'S RESPONSE?
# ============================================================

print("\n" + "=" * 70)
print("THE KEY TEST: One compression → fitter response")
print("=" * 70)

# The Planck fitter sees the compressed power spectrum and fits ΛCDM.
# The response pattern should be:
#
# 1. θ* unchanged (peaks don't move) → constrains ω_m h² × D_A normally ✓
# 2. Peak height ratios biased → ω_b shifts to compensate → r_d shifts → H₀ shifts
# 3. Damping tail suppressed → A_L > 1 to explain "extra smoothing"
# 4. Overall normalization affected → σ₈ shifts
#
# The RATIO of shifts should follow the covariance structure.
# Specifically:
#   ΔH₀/Δσ₈ should equal the published H₀-σ₈ degeneracy slope
#   ΔH₀/ΔA_L should equal the published H₀-A_L degeneracy slope

# Published Planck degeneracy slopes (from MCMC):
# dH₀/dσ₈ ≈ +36 (km/s/Mpc per unit σ₈, from correlation + sigmas)
dH0_dsigma8_planck = corr[0, 3] * planck_sigma['H0'] / planck_sigma['sigma8']
# dH₀/dA_L ≈ -1.7 (km/s/Mpc per unit A_L)
dH0_dAL_planck = corr[0, 6] * planck_sigma['H0'] / planck_sigma['A_L']

print(f"\nPublished degeneracy slopes (from Planck covariance):")
print(f"  dH₀/dσ₈ = {dH0_dsigma8_planck:+.1f} km/s/Mpc per unit σ₈")
print(f"  dH₀/dA_L = {dH0_dAL_planck:+.1f} km/s/Mpc per unit A_L")

# Now check: do the OBSERVED shifts follow these slopes?
ratio_H0_sigma8_obs = observed_shifts['H0'] / observed_shifts['sigma8']
ratio_H0_AL_obs = observed_shifts['H0'] / observed_shifts['A_L']

print(f"\nObserved shift ratios:")
print(f"  ΔH₀/Δσ₈ = {ratio_H0_sigma8_obs:+.1f}")
print(f"  ΔH₀/ΔA_L = {ratio_H0_AL_obs:+.1f}")

print(f"\nPredicted from covariance:")
print(f"  ΔH₀/Δσ₈ = {dH0_dsigma8_planck:+.1f}")
print(f"  ΔH₀/ΔA_L = {dH0_dAL_planck:+.1f}")

# If these match, the shifts ARE correlated through the fitter, not independent.
match_sigma8 = abs(ratio_H0_sigma8_obs - dH0_dsigma8_planck) / abs(dH0_dsigma8_planck)
match_AL = abs(ratio_H0_AL_obs - dH0_dAL_planck) / abs(dH0_dAL_planck)

print(f"\n  σ₈ slope agreement: {(1-match_sigma8)*100:.0f}%")
print(f"  A_L slope agreement: {(1-match_AL)*100:.0f}%")

# ============================================================
# THE SINGLE-PARAMETER COMPRESSION MODEL
# ============================================================

print("\n" + "=" * 70)
print("SINGLE-PARAMETER MODEL: Solve for compression amplitude")
print("=" * 70)

# Instead of solving for Γ_Σ per observable, solve for ONE compression
# amplitude that, propagated through the covariance, gives the best
# JOINT fit to all observed shifts.

# Model: the compression is a perturbation in the "peak height" direction
# in parameter space. The fitter absorbs it into the nearest ΛCDM parameters.
# 
# The peak-height direction primarily involves:
#   ω_b (weight ~0.5): peak heights = baryon loading
#   τ (weight ~0.2): overall normalization
#   n_s (weight ~0.2): tilt of power spectrum
#   A_L (weight ~-0.1): lensing smoothing absorbs residual

# Compression perturbation vector (in σ units)
# This is the direction in parameter space that "peak suppression" corresponds to
# Negative ω_b (less baryons), positive τ (more suppression), positive A_L (smooth)
compression_direction = np.array([
    0.0,    # H₀ (derived, not directly perturbed)
    -0.5,   # ω_b (primary: less baryon loading)
    0.0,    # ω_c (not directly affected)
    0.0,    # σ₈ (derived)
    -0.2,   # n_s (slight red tilt)
    0.2,    # τ (absorbs some suppression)
    0.1,    # A_L (smoothing compensator)
])

# Propagate through covariance: Δp = C × (C_direct)^{-1} × compression_direction
# More simply: the response is the compression direction PLUS its covariance-mediated
# projection onto the constrained parameters (H₀, σ₈).

# For derived parameters (H₀, σ₈), their shift is determined by their correlations
# with the directly perturbed parameters (ω_b, n_s, τ, A_L).

# Δp_derived = Σ_j corr(derived, j) × (σ_derived/σ_j) × Δp_j
# where Δp_j is the direct perturbation to parameter j.

# Let α be the amplitude (in units where the perturbation vector is unit)
# We solve for α by matching ΔH₀ = -5.68

# ΔH₀(α) = α × Σ_j [compression_direction_j × corr(H₀, j) × σ_H₀ / σ_j × σ_j]
# where the sum is over directly-perturbed parameters

# Wait, need to be more careful. The compression direction is in SIGMA units.
# Δp_j = α × compression_direction_j × σ_j
# ΔH₀ = Σ_j [corr(H₀, j) × (σ_H₀) × (Δp_j / σ_j)]
#       = α × σ_H₀ × Σ_j [corr(H₀, j) × compression_direction_j]

idx = {p: i for i, p in enumerate(param_names)}

# Compute the H₀ response coefficient
H0_response = 0
for j, pname in enumerate(param_names):
    if compression_direction[j] != 0:
        H0_response += corr[idx['H0'], j] * compression_direction[j]

sigma8_response = 0
for j, pname in enumerate(param_names):
    if compression_direction[j] != 0:
        sigma8_response += corr[idx['sigma8'], j] * compression_direction[j]

AL_response = compression_direction[idx['A_L']]  # A_L is directly perturbed

print(f"Response coefficients (per unit α):")
print(f"  ΔH₀:  {H0_response:+.3f} × σ_H₀ × α = {H0_response * planck_sigma['H0']:+.3f} × α km/s/Mpc")
print(f"  Δσ₈:  {sigma8_response:+.3f} × σ_σ₈ × α = {sigma8_response * planck_sigma['sigma8']:+.4f} × α")
print(f"  ΔA_L: {AL_response:+.3f} × σ_A_L × α = {AL_response * planck_sigma['A_L']:+.4f} × α")

# Solve for α from H₀:
if H0_response * planck_sigma['H0'] != 0:
    alpha_H0 = observed_shifts['H0'] / (H0_response * planck_sigma['H0'])
else:
    alpha_H0 = float('inf')

print(f"\nSolving for α from ΔH₀ = {observed_shifts['H0']:.2f}:")
print(f"  α = {alpha_H0:.1f}")

# Predict the OTHER shifts from this α:
pred_sigma8 = alpha_H0 * sigma8_response * planck_sigma['sigma8']
pred_AL = alpha_H0 * AL_response * planck_sigma['A_L']

# Also compute all parameter shifts
print(f"\nPredicted shifts from α = {alpha_H0:.1f}:")
print(f"  {'Parameter':<15} {'Predicted':>12} {'Observed':>12} {'Match':>8}")
print(f"  {'-'*50}")
print(f"  {'ΔH₀':<15} {alpha_H0 * H0_response * planck_sigma['H0']:>+12.2f} {observed_shifts['H0']:>+12.2f} {'(input)':>8}")
print(f"  {'Δσ₈':<15} {pred_sigma8:>+12.4f} {observed_shifts['sigma8']:>+12.4f} {'✓' if abs(pred_sigma8 - observed_shifts['sigma8']) < 0.03 else '✗':>8}")
print(f"  {'ΔA_L':<15} {pred_AL:>+12.4f} {observed_shifts['A_L']:>+12.4f} {'✓' if abs(pred_AL - observed_shifts['A_L']) < 0.1 else '✗':>8}")

# Also show what ω_b shift this implies
pred_omega_b = alpha_H0 * compression_direction[idx['omega_b']] * planck_sigma['omega_b']
print(f"  {'Δω_b':<15} {pred_omega_b:>+12.6f} {'(prediction)':>12}")
print(f"  {'True ω_b':<15} {planck_params['omega_b'] - pred_omega_b:>12.5f} {'(corrected)':>12}")

# ============================================================
# SENSITIVITY TO COMPRESSION DIRECTION
# ============================================================

print("\n" + "=" * 70)
print("SENSITIVITY: Different compression direction models")
print("=" * 70)

directions = {
    'ω_b only': np.array([0, -1, 0, 0, 0, 0, 0]),
    'ω_b + τ': np.array([0, -0.7, 0, 0, 0, 0.3, 0]),
    'ω_b + n_s': np.array([0, -0.6, 0, 0, -0.4, 0, 0]),
    'ω_b + τ + A_L': np.array([0, -0.5, 0, 0, -0.2, 0.2, 0.1]),
    'Damping tail (ω_b + τ + n_s)': np.array([0, -0.5, 0, 0, -0.3, 0.2, 0]),
    'Pure A_L': np.array([0, 0, 0, 0, 0, 0, 1]),
}

print(f"\n{'Model':<35} {'α_H₀':>8} {'Δσ₈_pred':>10} {'ΔA_L_pred':>10} {'σ₈ match':>10} {'A_L match':>10}")
print("-" * 90)

for name, direction in directions.items():
    # H₀ response
    h0_resp = sum(corr[idx['H0'], j] * direction[j] for j in range(n_params))
    s8_resp = sum(corr[idx['sigma8'], j] * direction[j] for j in range(n_params))
    al_resp = direction[idx['A_L']]
    
    if abs(h0_resp * planck_sigma['H0']) > 1e-10:
        alpha = observed_shifts['H0'] / (h0_resp * planck_sigma['H0'])
        ds8 = alpha * s8_resp * planck_sigma['sigma8']
        dal = alpha * al_resp * planck_sigma['A_L']
        
        s8_match = "✓" if abs(ds8 - observed_shifts['sigma8']) < 0.03 else "~" if abs(ds8 - observed_shifts['sigma8']) < 0.06 else "✗"
        al_match = "✓" if abs(dal - observed_shifts['A_L']) < 0.1 else "~" if abs(dal - observed_shifts['A_L']) < 0.2 else "✗"
        
        print(f"  {name:<33} {alpha:>+8.1f} {ds8:>+10.4f} {dal:>+10.4f} {s8_match:>10} {al_match:>10}")
    else:
        print(f"  {name:<33} {'∞':>8} {'—':>10} {'—':>10} {'—':>10} {'—':>10}")

print(f"\n  Targets: Δσ₈ = {observed_shifts['sigma8']:+.4f}, ΔA_L = {observed_shifts['A_L']:+.4f}")

# ============================================================
# THE ACID TEST: JOINT χ² FOR EACH MODEL
# ============================================================

print("\n" + "=" * 70)
print("JOINT χ²: Which compression direction best explains ALL shifts?")
print("=" * 70)

obs_vec = np.array([observed_shifts['H0'], observed_shifts['sigma8'], observed_shifts['A_L']])
obs_sigma = np.array([planck_sigma['H0'], planck_sigma['sigma8'], planck_sigma['A_L']])

best_chi2 = 1e10
best_model = None

for name, direction in directions.items():
    h0_resp = sum(corr[idx['H0'], j] * direction[j] for j in range(n_params))
    s8_resp = sum(corr[idx['sigma8'], j] * direction[j] for j in range(n_params))
    al_resp = direction[idx['A_L']]
    
    if abs(h0_resp * planck_sigma['H0']) < 1e-10:
        continue
    
    # Scan α for minimum joint χ²
    best_alpha = None
    best_chi2_model = 1e10
    
    for alpha in np.linspace(-50, 50, 10000):
        pred = np.array([
            alpha * h0_resp * planck_sigma['H0'],
            alpha * s8_resp * planck_sigma['sigma8'],
            alpha * al_resp * planck_sigma['A_L']
        ])
        chi2 = np.sum(((pred - obs_vec) / obs_sigma) ** 2)
        if chi2 < best_chi2_model:
            best_chi2_model = chi2
            best_alpha = alpha
    
    pred_best = np.array([
        best_alpha * h0_resp * planck_sigma['H0'],
        best_alpha * s8_resp * planck_sigma['sigma8'],
        best_alpha * al_resp * planck_sigma['A_L']
    ])
    
    print(f"\n  {name}:")
    print(f"    Best α = {best_alpha:.1f}")
    print(f"    ΔH₀ = {pred_best[0]:+.2f} (obs: {obs_vec[0]:+.2f})")
    print(f"    Δσ₈ = {pred_best[1]:+.4f} (obs: {obs_vec[1]:+.4f})")
    print(f"    ΔA_L = {pred_best[2]:+.4f} (obs: {obs_vec[2]:+.4f})")
    print(f"    Joint χ² = {best_chi2_model:.2f} (3 dof)")
    
    if best_chi2_model < best_chi2:
        best_chi2 = best_chi2_model
        best_model = name

print(f"\n★ BEST MODEL: {best_model} (χ² = {best_chi2:.2f})")

# ============================================================
# SAVE RESULTS
# ============================================================

results_dir = Path("results_cmb_propagation")
results_dir.mkdir(exist_ok=True)

results = {
    'cmb_ladder': [(name, q, nature, effect) for name, q, nature, effect in cmb_ladder],
    'observed_shifts': observed_shifts,
    'compression_directions_tested': list(directions.keys()),
    'best_model': best_model,
    'best_chi2': float(best_chi2),
    'degeneracy_slopes': {
        'dH0_dsigma8_planck': float(dH0_dsigma8_planck),
        'dH0_dAL_planck': float(dH0_dAL_planck),
        'dH0_dsigma8_observed': float(ratio_H0_sigma8_obs),
        'dH0_dAL_observed': float(ratio_H0_AL_obs),
    }
}

with open(results_dir / "cmb_propagation_results.json", 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults saved to {results_dir}/")
print("\n" + "=" * 70)
print("END")
print("=" * 70)
