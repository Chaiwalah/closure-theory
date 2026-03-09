#!/usr/bin/env python3
"""
TRANSITION SUSCEPTIBILITY INDEX — The Unified Selector
=======================================================
Compute the variance of the transition magnetic moment σ_μ²
for each ladder line using Clebsch-Gordan coefficients.

This measures the SPREAD of Zeeman frequency shifts across
all allowed m_u → m_l transitions, weighted by transition
probabilities.

Hypothesis: σ_μ² is the unified quantum property that
reproduces diagnostic sensitivity q (r = -0.975).

Author: Clawd
Date: 7 March 2026
"""

import numpy as np
from scipy import stats
from itertools import product
import json
import os

RESULTS_DIR = '/root/clawd/projects/closure-theory/results_gfactor'
os.makedirs(RESULTS_DIR, exist_ok=True)

def lande_g(J, L, S):
    """Landé g-factor from quantum numbers."""
    if J == 0:
        return 0.0
    return 1.0 + (J*(J+1) + S*(S+1) - L*(L+1)) / (2*J*(J+1))

def clebsch_gordan_sq(J_upper, m_upper, J_lower, m_lower, delta_m):
    """
    Approximate relative transition probability for 
    electric dipole (E1) transitions: |<J_l m_l|T^1_q|J_u m_u>|²
    
    Uses the 3j-symbol squared, proportional to the Wigner 3j:
    (J_u  1  J_l)²
    (-m_u q  m_l)
    
    For our purposes, the key selection rules are:
    delta_m = m_lower - m_upper = q (q = -1, 0, +1)
    |J_upper - J_lower| <= 1
    """
    q = m_lower - m_upper
    if abs(q) > 1:
        return 0.0
    
    # Use analytical 3j symbols for common cases
    # For J → J transitions (ΔJ = 0)
    if J_upper == J_lower:
        J = J_upper
        if J == 0:
            return 0.0
        if q == 0:  # π component
            return m_upper**2 / (J * (J + 1))
        elif abs(q) == 1:  # σ components
            return (J*(J+1) - m_upper*(m_upper + q)) / (2 * J * (J + 1))
        return 0.0
    
    # For J → J-1 transitions (ΔJ = -1)
    elif J_lower == J_upper - 1:
        J = J_upper
        if q == 0:  # π
            return (J**2 - m_upper**2) / (J * (2*J - 1) * (2*J + 1) / (2*J))
        elif q == +1:  # σ+
            return ((J - m_upper) * (J - m_upper - 1)) / (2 * J * (2*J + 1))
        elif q == -1:  # σ-
            return ((J + m_upper) * (J + m_upper - 1)) / (2 * J * (2*J + 1))
        return 0.0
    
    # For J → J+1 transitions (ΔJ = +1)
    elif J_lower == J_upper + 1:
        J = J_upper
        if q == 0:  # π
            return ((J + 1)**2 - m_upper**2) / ((J + 1) * (2*J + 1) * (2*J + 3) / (2*(J+1)))
        elif q == +1:  # σ+
            return ((J + m_upper + 1) * (J + m_upper + 2)) / (2 * (J+1) * (2*J + 1))
        elif q == -1:  # σ-
            return ((J - m_upper + 1) * (J - m_upper + 2)) / (2 * (J+1) * (2*J + 1))
        return 0.0
    
    return 0.0


def compute_zeeman_variance(J_upper, L_upper, S_upper, J_lower, L_lower, S_lower,
                             transition_type='E1'):
    """
    Compute the variance of the transition magnetic moment σ_μ².
    
    For each allowed m_u → m_l transition:
    - Compute the frequency shift: Δν ∝ (m_u × g_u - m_l × g_l)
    - Weight by transition probability (Clebsch-Gordan squared)
    - Compute variance of the weighted distribution
    
    This is the "spread" of the Zeeman pattern.
    """
    g_u = lande_g(J_upper, L_upper, S_upper)
    g_l = lande_g(J_lower, L_lower, S_lower)
    
    # All allowed m values
    m_upper_vals = np.arange(-J_upper, J_upper + 1, 1)
    m_lower_vals = np.arange(-J_lower, J_lower + 1, 1)
    
    shifts = []
    weights = []
    components = []
    
    for m_u in m_upper_vals:
        for m_l in m_lower_vals:
            delta_m = m_l - m_u
            
            # Selection rules for E1: Δm = 0, ±1
            if transition_type == 'E1' and abs(delta_m) > 1:
                continue
            # For M1: Δm = 0, ±1 (same as E1 for angular part)
            if transition_type == 'M1' and abs(delta_m) > 1:
                continue
            # For E2: Δm = 0, ±1, ±2
            if transition_type == 'E2' and abs(delta_m) > 2:
                continue
            
            # Frequency shift in units of μ_B × B
            shift = m_u * g_u - m_l * g_l
            
            # Weight (simplified — use uniform for E2/M1, CG for E1)
            if transition_type == 'E1':
                w = clebsch_gordan_sq(J_upper, m_u, J_lower, m_l, delta_m)
            else:
                # For forbidden lines, use simpler weighting
                # (proper calculation would need higher-order multipole matrix elements)
                w = 1.0  # uniform weight as approximation
            
            if w > 0:
                shifts.append(shift)
                weights.append(w)
                components.append({
                    'm_u': float(m_u), 'm_l': float(m_l),
                    'delta_m': int(delta_m),
                    'shift': float(shift), 'weight': float(w)
                })
    
    if not shifts:
        return 0.0, 0, g_u, g_l, []
    
    shifts = np.array(shifts)
    weights = np.array(weights)
    weights = weights / weights.sum()  # normalize
    
    # Weighted mean and variance
    mean_shift = np.sum(weights * shifts)
    variance = np.sum(weights * (shifts - mean_shift)**2)
    
    n_components = len(shifts)
    
    return variance, n_components, g_u, g_l, components


# ============================================================
# COMPUTE FOR ALL LADDER LINES
# ============================================================

print("=" * 80)
print("TRANSITION SUSCEPTIBILITY INDEX — σ_μ² (Zeeman Variance)")
print("=" * 80)

# Define transitions with full quantum numbers
# Format: (name, J_u, L_u, S_u, J_l, L_l, S_l, transition_type, degradation)

transitions = [
    # [NII] 6584: ¹D₂ → ³P₂ (forbidden, M1/E2)
    # N⁺: ground config [He] 2s² 2p²
    # ¹D₂: L=2, S=0, J=2
    # ³P₂: L=1, S=1, J=2
    ('[NII] 6584',  2, 2, 0,  2, 1, 1, 'M1', 0.000),
    
    # [OIII] 5007: ¹D₂ → ³P₂ (forbidden, M1/E2)
    # O²⁺: same term structure as N⁺
    ('[OIII] 5007', 2, 2, 0,  2, 1, 1, 'M1', 0.000),
    
    # Hβ 4861: principal series, multiple transitions
    # Dominant: 4d ²D₅/₂ → 2p ²P₃/₂ and similar
    # Hydrogen has no fine structure splitting in g (all g ≈ 1 for given J)
    # Use the dominant component: J=5/2→3/2
    ('Hβ 4861',     2.5, 2, 0.5,  1.5, 1, 0.5, 'E1', -0.038),
    # Also test J=3/2→1/2 component
    
    # [OII] 3726: ²D₃/₂ → ⁴S₃/₂ (forbidden, E2)
    ('[OII] 3726',  1.5, 2, 0.5,  1.5, 0, 1.5, 'E2', -0.179),
    
    # [OII] 3729: ²D₅/₂ → ⁴S₃/₂ (forbidden, E2)
    ('[OII] 3729',  2.5, 2, 0.5,  1.5, 0, 1.5, 'E2', -0.179),
    
    # CIV 1548: ²P₃/₂ → ²S₁/₂ (permitted, E1)
    ('CIV 1548',    1.5, 1, 0.5,  0.5, 0, 0.5, 'E1', -0.289),
    
    # CIV 1550: ²P₁/₂ → ²S₁/₂ (permitted, E1)
    ('CIV 1550',    0.5, 1, 0.5,  0.5, 0, 0.5, 'E1', -0.289),
    
    # [SII] 6716: ²D₃/₂ → ⁴S₃/₂ (forbidden, E2)
    ('[SII] 6716',  1.5, 2, 0.5,  1.5, 0, 1.5, 'E2', -0.396),
    
    # [SII] 6731: ²D₅/₂ → ⁴S₃/₂ (forbidden, E2)
    ('[SII] 6731',  2.5, 2, 0.5,  1.5, 0, 1.5, 'E2', -0.396),
]

print(f"\n{'Line':<16} {'Type':<5} {'J_u→J_l':<10} {'g_u':<7} {'g_l':<7} "
      f"{'|Δg|':<7} {'N_comp':<8} {'σ_μ²':<10} {'Degradation':<12}")
print("-" * 95)

all_results = []
for name, Ju, Lu, Su, Jl, Ll, Sl, ttype, deg in transitions:
    var, n_comp, gu, gl, components = compute_zeeman_variance(
        Ju, Lu, Su, Jl, Ll, Sl, ttype)
    
    delta_g = abs(gu - gl)
    
    print(f"{name:<16} {ttype:<5} {Ju}→{Jl:<6} {gu:<7.4f} {gl:<7.4f} "
          f"{delta_g:<7.4f} {n_comp:<8} {var:<10.4f} {deg:<12.3f}")
    
    all_results.append({
        'name': name, 'type': ttype,
        'J_upper': float(Ju), 'J_lower': float(Jl),
        'g_upper': float(gu), 'g_lower': float(gl),
        'delta_g': float(delta_g),
        'n_components': int(n_comp),
        'sigma_mu_sq': float(var),
        'degradation': float(deg),
    })

# ============================================================
# COMPUTE AVERAGED VALUES FOR THE 6-LINE LADDER
# ============================================================

print("\n" + "=" * 80)
print("6-LINE LADDER — AVERAGED DOUBLET VALUES")
print("=" * 80)

# Average doublet components for lines that are doublets
ladder_averaged = [
    {'name': '[NII] 6584', 'sigma_mu_sq': all_results[0]['sigma_mu_sq'],
     'degradation': 0.000, 'delta_g': all_results[0]['delta_g'],
     'n_comp': all_results[0]['n_components'], 'type': 'M1'},
    {'name': '[OIII] 5007', 'sigma_mu_sq': all_results[1]['sigma_mu_sq'],
     'degradation': 0.000, 'delta_g': all_results[1]['delta_g'],
     'n_comp': all_results[1]['n_components'], 'type': 'M1'},
    {'name': 'Hβ 4861', 'sigma_mu_sq': all_results[2]['sigma_mu_sq'],
     'degradation': -0.038, 'delta_g': all_results[2]['delta_g'],
     'n_comp': all_results[2]['n_components'], 'type': 'E1'},
    {'name': '[OII] 3727',
     'sigma_mu_sq': (all_results[3]['sigma_mu_sq'] + all_results[4]['sigma_mu_sq'])/2,
     'degradation': -0.179,
     'delta_g': (all_results[3]['delta_g'] + all_results[4]['delta_g'])/2,
     'n_comp': (all_results[3]['n_components'] + all_results[4]['n_components'])//2,
     'type': 'E2'},
    {'name': 'CIV 1549',
     'sigma_mu_sq': (all_results[5]['sigma_mu_sq'] + all_results[6]['sigma_mu_sq'])/2,
     'degradation': -0.289,
     'delta_g': (all_results[5]['delta_g'] + all_results[6]['delta_g'])/2,
     'n_comp': (all_results[5]['n_components'] + all_results[6]['n_components'])//2,
     'type': 'E1'},
    {'name': '[SII] 6718',
     'sigma_mu_sq': (all_results[7]['sigma_mu_sq'] + all_results[8]['sigma_mu_sq'])/2,
     'degradation': -0.396,
     'delta_g': (all_results[7]['delta_g'] + all_results[8]['delta_g'])/2,
     'n_comp': (all_results[7]['n_components'] + all_results[8]['n_components'])//2,
     'type': 'E2'},
]

print(f"\n{'Line':<16} {'Type':<5} {'|Δg|':<7} {'N_comp':<8} {'σ_μ²':<10} {'Degradation':<12}")
print("-" * 58)
for line in ladder_averaged:
    print(f"{line['name']:<16} {line['type']:<5} {line['delta_g']:<7.4f} "
          f"{line['n_comp']:<8} {line['sigma_mu_sq']:<10.4f} {line['degradation']:<12.3f}")

# ============================================================
# CORRELATION ANALYSIS
# ============================================================

print("\n" + "=" * 80)
print("CORRELATION ANALYSIS — Finding the Unified Selector")
print("=" * 80)

deg = np.array([l['degradation'] for l in ladder_averaged])
sigma_sq = np.array([l['sigma_mu_sq'] for l in ladder_averaged])
delta_g = np.array([l['delta_g'] for l in ladder_averaged])
n_comp = np.array([l['n_comp'] for l in ladder_averaged], dtype=float)

# Diagnostic sensitivity q for reference
q_values = np.array([0.0, 0.0, 0.3, 0.4, 0.6, 0.7])

tests = [
    ('Diagnostic sensitivity q', q_values),
    ('σ_μ² (Zeeman variance)', sigma_sq),
    ('σ_μ⁴ (variance squared)', sigma_sq**2),
    ('|Δg|', delta_g),
    ('|Δg|²', delta_g**2),
    ('N_components', n_comp),
    ('σ_μ² × N_comp', sigma_sq * n_comp),
    ('√(σ_μ²)', np.sqrt(sigma_sq)),
]

print(f"\n{'Variable':<30} {'Pearson r':<12} {'p-value':<12} {'Spearman r':<12} {'p-value':<12}")
print("-" * 78)

for name, vals in tests:
    if np.std(vals) == 0:
        print(f"{name:<30} {'N/A (constant)':<12}")
        continue
    r, p = stats.pearsonr(vals, deg)
    rs, ps = stats.spearmanr(vals, deg)
    marker = " ◄◄◄" if abs(r) > 0.95 else (" ◄◄" if abs(r) > 0.9 else "")
    print(f"{name:<30} {r:+.4f}      {p:.6f}    {rs:+.4f}      {ps:.6f}{marker}")

# ============================================================
# MULTIPOLE ORDER ANALYSIS
# ============================================================

print("\n" + "=" * 80)
print("MULTIPOLE ORDER ANALYSIS")
print("=" * 80)

# Check if transition type alone predicts degradation
print("\nTransition type vs degradation:")
for line in ladder_averaged:
    print(f"  {line['name']:<16} {line['type']:<5} → {line['degradation']:+.3f}")

# Assign numerical multipole order
# M1 = 1 (magnetic dipole — forbidden)
# E1 = 2 (electric dipole — permitted)
# E2 = 3 (electric quadrupole — forbidden)
multipole_map = {'M1': 1, 'E1': 2, 'E2': 3}
multipole = np.array([multipole_map[l['type']] for l in ladder_averaged], dtype=float)

r_mp, p_mp = stats.pearsonr(multipole, deg)
print(f"\nMultipole order vs degradation: r = {r_mp:+.4f}, p = {p_mp:.4f}")

# Alternative encoding: E2 > E1 > M1 (by degradation)
multipole_v2 = np.array([0, 0, 1, 2, 1, 2], dtype=float)  # M1=0, E1=1, E2=2
r_mp2, p_mp2 = stats.pearsonr(multipole_v2, deg)
print(f"Multipole (M1=0,E1=1,E2=2) vs degradation: r = {r_mp2:+.4f}, p = {p_mp2:.4f}")

# ============================================================
# COMPOSITE INDEX: |Δg|² × multipole_factor
# ============================================================

print("\n" + "=" * 80)
print("COMPOSITE INDICES — Searching for the Holy Grail")
print("=" * 80)

# Try various combinations
# The key: what separates [OII], CIV, [SII] with |Δg|=1?
# Their critical densities differ:
# [OII]: n_crit ~ 3.4 × 10³ cm⁻³
# CIV: permitted line, n_crit effectively very high
# [SII]: n_crit ~ 1.5 × 10³ cm⁻³ (6716) / 3.9 × 10³ (6731)

# Log critical density (approximate)
log_ncrit = np.array([
    5.0,   # [NII]: ~10⁵
    5.8,   # [OIII]: ~7×10⁵
    14.0,  # Hβ: permitted, ~10¹⁴ (effectively infinite)
    3.5,   # [OII]: ~3×10³
    14.0,  # CIV: permitted, effectively infinite
    3.2,   # [SII]: ~1.5×10³
])

# Inverse of log n_crit (lower crit density = more fragile)
inv_logncrit = 1.0 / log_ncrit

# Einstein A coefficient (s⁻¹) — radiative decay rate
A_coeff = np.array([
    1.01e-3,   # [NII] 6584
    6.21e-3,   # [OIII] 5007
    1.98e8,    # Hβ (permitted)
    1.30e-4,   # [OII] 3727 (avg)
    2.65e8,    # CIV 1549 (avg)
    2.60e-4,   # [SII] 6718 (avg)
])

# Radiative lifetime
tau_rad = 1.0 / A_coeff

composite_tests = [
    ('|Δg|² × 1/log(n_crit)', delta_g**2 * inv_logncrit),
    ('σ_μ² × 1/log(n_crit)', sigma_sq * inv_logncrit),
    ('|Δg|² × log(τ_rad)', delta_g**2 * np.log10(tau_rad)),
    ('σ_μ² × log(τ_rad)', sigma_sq * np.log10(tau_rad)),
    ('|Δg|² × τ_rad^0.1', delta_g**2 * tau_rad**0.1),
    ('1/log(n_crit) alone', inv_logncrit),
    ('log(τ_rad) alone', np.log10(tau_rad)),
]

print(f"\n{'Composite Index':<30} {'Pearson r':<12} {'p-value':<12}")
print("-" * 54)
for name, vals in composite_tests:
    if np.std(vals) == 0:
        continue
    r, p = stats.pearsonr(vals, deg)
    marker = " ◄◄◄ HOLY GRAIL?" if abs(r) > 0.97 else (" ◄◄◄" if abs(r) > 0.95 else "")
    print(f"{name:<30} {r:+.4f}      {p:.6f}{marker}")

# ============================================================
# SAVE RESULTS
# ============================================================

output = {
    'individual_lines': all_results,
    'ladder_averaged': ladder_averaged,
    'correlations': {
        'q_vs_deg': float(stats.pearsonr(q_values, deg)[0]),
        'sigma_sq_vs_deg': float(stats.pearsonr(sigma_sq, deg)[0]),
        'delta_g_sq_vs_deg': float(stats.pearsonr(delta_g**2, deg)[0]),
    }
}

with open(os.path.join(RESULTS_DIR, 'susceptibility_results.json'), 'w') as f:
    json.dump(output, f, indent=2, default=str)

print(f"\n\nResults saved to {RESULTS_DIR}/susceptibility_results.json")
