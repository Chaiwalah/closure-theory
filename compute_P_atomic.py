#!/usr/bin/env python3
"""
COMPUTE P FROM FIRST PRINCIPLES
================================
GPT's prediction: P = emissivity Jacobian × trapping leverage × (1 - branch lock)

Using PyNeb for atomic data + published values for escape/trapping.
NO fitting to our data. Pure atomic physics. Then compare to empirical ranking.
"""

import numpy as np

# ============================================================
# ATOMIC DATA — from PyNeb, CHIANTI, and standard references
# ============================================================

# Critical densities (cm⁻³) — where collisional de-excitation = radiative
# Lines near their critical density have steep ∂lnj/∂ln(n_e)
critical_densities = {
    'LYA':    None,      # Resonance line, not collisionally de-excited in same way
    'CIV':    None,      # Resonance, permitted
    'CIII':   None,      # Semi-forbidden (intercombination)
    'MGII':   None,      # Resonance
    'OII':    3.4e3,     # [OII] 3728 (actually doublet 3726/3729)
    'HBETA':  None,      # Recombination
    'OIII':   7.0e5,     # [OIII] 5007
    'HALPHA': None,      # Recombination
    'NII':    8.6e4,     # [NII] 6583
    'SII':    1.5e3,     # [SII] 6716 (and 3.6e3 for 6731)
}

# Typical AGN NLR/BLR densities: n_e ~ 10²-10⁴ cm⁻³ (NLR), 10⁹-10¹¹ (BLR)

# ============================================================
# COMPONENT 1: Emissivity gradient ‖∇ ln j‖
# ============================================================
# For forbidden lines: j ∝ n_e × Ω(T)/g × exp(-E/kT) / (1 + n_e/n_crit)
# ∂ln j/∂ln n_e = 1 - n_e/(n_e + n_crit) = n_crit/(n_e + n_crit)
# At n_e << n_crit: gradient → 1 (linear, not sensitive)
# At n_e ~ n_crit: gradient → 0.5 (STEEPEST response region)
# At n_e >> n_crit: gradient → 0 (suppressed)

# For [SII]: n_crit = 1500, typical NLR n_e ~ 300-3000 → RIGHT at critical density
# For [OIII]: n_crit = 700,000, typical NLR n_e ~ 300-3000 → FAR below critical density
# For [NII]: n_crit = 86,000, typical NLR n_e ~ 300-3000 → FAR below critical density

# Temperature sensitivity: ∂ln j/∂ln T
# For collisionally excited lines: j ∝ exp(-E_exc/kT)
# ∂ln j/∂ln T = E_exc/(kT) 
# Higher excitation energy → MORE temperature sensitive

# Excitation energies (eV)
excitation_energy_eV = {
    'LYA':    10.2,    # Lyα (1s→2p)
    'CIV':    8.0,     # CIV 1549
    'CIII':   6.5,     # CIII] 1909
    'MGII':   4.4,     # MgII 2798
    'OII':    3.3,     # [OII] 3728
    'HBETA':  2.55,    # Hβ (recombination, different physics)
    'OIII':   2.48,    # [OIII] 5007
    'HALPHA': 1.89,    # Hα
    'NII':    1.89,    # [NII] 6583
    'SII':    1.84,    # [SII] 6716
}

# At T = 10,000 K (typical NLR): kT = 0.86 eV
kT = 0.86  # eV at 10,000 K

# ∂ln j/∂ln T for each line
def temp_sensitivity(E_eV, kT=0.86):
    return E_eV / kT

# Density sensitivity for forbidden lines
# = n_crit / (n_e + n_crit) evaluated at typical NLR density
def density_sensitivity(n_crit, n_e_typical=1000):
    if n_crit is None:
        return 0.5  # Permitted/resonance lines: different physics
    # At typical NLR density
    grad = n_crit / (n_e_typical + n_crit)
    # But what matters is HOW MUCH j changes across the density range
    # Lines near critical density have steepest gradient
    # Compute over range n_e = 100 to 10000
    n_range = np.logspace(2, 4, 100)
    j_relative = n_range / (1 + n_range / n_crit)
    # Normalize
    j_norm = j_relative / j_relative.max()
    # Dynamic range of j across the density range
    dynamic_range = j_norm.max() / j_norm.min()
    return np.log10(dynamic_range)

# ============================================================
# COMPONENT 2: Escape/trapping leverage (E)
# ============================================================
# Resonance lines: large optical depth, many scatterings
# Forbidden lines: optically thin (by definition), E ~ 1
# Recombination lines: moderate

trapping_leverage = {
    'LYA':    10.0,    # Archetypal trapped line. τ ~ 10⁴-10⁷ in AGN
    'CIV':    5.0,     # Resonance, τ ~ 10²-10⁴
    'CIII':   2.0,     # Semi-forbidden, moderate τ
    'MGII':   3.0,     # Resonance doublet, τ ~ 10-10³
    'OII':    1.0,     # Forbidden, optically thin
    'HBETA':  1.5,     # Recombination, moderate τ in BLR
    'OIII':   1.0,     # Forbidden, optically thin
    'HALPHA': 2.0,     # Can be optically thick in dense BLR
    'NII':    1.0,     # Forbidden, optically thin
    'SII':    1.0,     # Forbidden, optically thin
}

# ============================================================
# COMPONENT 3: Branching lock (B)
# ============================================================
# Lines from same upper level: ratio is FIXED by A-coefficients
# [OIII] 5007/4959 = 2.98 (fixed)
# [NII] 6583/6548 = 3.05 (fixed)
# [SII] 6716/6731 — DIFFERENT upper levels! NOT locked!
# [OII] 3726/3729 — DIFFERENT upper levels (like SII)

branch_lock = {
    'LYA':    0.0,    # No branching lock
    'CIV':    0.0,    # Doublet but different physics
    'CIII':   0.0,    # No lock
    'MGII':   0.0,    # Doublet, some locking but weak
    'OII':    0.0,    # DIFFERENT upper levels (density diagnostic!)
    'HBETA':  0.0,    # No lock (but Balmer decrement has some)
    'OIII':   0.95,   # STRONGLY locked (5007/4959 = 2.98, invariant)
    'HALPHA': 0.3,    # Partial — Balmer series has branching constraints
    'NII':    0.95,   # STRONGLY locked (6583/6548 = 3.05, invariant)
    'SII':    0.0,    # NOT locked — different upper levels!
}

# ============================================================
# COMPUTE P_theoretical
# ============================================================
print("=" * 80)
print("COMPUTING P FROM FIRST PRINCIPLES")
print("P = (density_sensitivity + temp_sensitivity) × trapping × (1 - branch_lock)")
print("=" * 80)

lines = ['LYA', 'CIV', 'CIII', 'MGII', 'OII', 'HBETA', 'OIII', 'HALPHA', 'NII', 'SII']

results = {}
for name in lines:
    n_crit = critical_densities[name]
    E = excitation_energy_eV[name]
    
    # Density sensitivity
    dens_sens = density_sensitivity(n_crit)
    
    # Temperature sensitivity
    temp_sens = temp_sensitivity(E)
    
    # Total microstate sensitivity (sum of gradients)
    micro_sens = dens_sens + temp_sens
    
    # Trapping
    trap = trapping_leverage[name]
    
    # Branch lock
    lock = branch_lock[name]
    
    # P_theoretical
    P_theo = micro_sens * trap * (1 - lock)
    
    results[name] = {
        'dens_sens': dens_sens,
        'temp_sens': temp_sens,
        'micro_sens': micro_sens,
        'trap': trap,
        'lock': lock,
        'P_raw': P_theo,
    }

# Normalize P to [0, 1]
P_values = np.array([results[n]['P_raw'] for n in lines])
P_min, P_max = P_values.min(), P_values.max()
for name in lines:
    results[name]['P_norm'] = (results[name]['P_raw'] - P_min) / (P_max - P_min)

# Empirical P from our data (from universal_v2.py)
P_empirical = {
    'LYA':    1.000,
    'SII':    0.983,
    'CIV':    0.930,
    'MGII':   0.875,
    'OIII':   0.864,
    'OII':    0.845,
    'CIII':   0.791,
    'HBETA':  0.517,
    'HALPHA': 0.124,
    'NII':    0.000,
}

# ============================================================
# COMPARISON
# ============================================================
print(f"\n{'Line':>8} {'dens':>6} {'temp':>6} {'micro':>7} {'trap':>5} {'lock':>5} "
      f"{'P_raw':>7} {'P_theo':>7} {'P_emp':>7} {'match':>7}")
print("-" * 80)

theo_arr = []
emp_arr = []

for name in lines:
    r = results[name]
    p_e = P_empirical[name]
    p_t = r['P_norm']
    
    theo_arr.append(p_t)
    emp_arr.append(p_e)
    
    # Check if ranking matches
    match = "—"
    
    print(f"{name:>8} {r['dens_sens']:>6.2f} {r['temp_sens']:>6.2f} {r['micro_sens']:>7.2f} "
          f"{r['trap']:>5.1f} {r['lock']:>5.2f} {r['P_raw']:>7.2f} {p_t:>7.3f} {p_e:>7.3f}")

# Correlation
from scipy.stats import spearmanr, pearsonr

theo_arr = np.array(theo_arr)
emp_arr = np.array(emp_arr)

r_spearman, p_spearman = spearmanr(theo_arr, emp_arr)
r_pearson, p_pearson = pearsonr(theo_arr, emp_arr)

print(f"\n{'='*80}")
print(f"CORRELATION: THEORETICAL P vs EMPIRICAL P")
print(f"{'='*80}")
print(f"  Spearman: r = {r_spearman:.4f}  (p = {p_spearman:.4e})")
print(f"  Pearson:  r = {r_pearson:.4f}  (p = {p_pearson:.4e})")

# Rank comparison
theo_rank = np.argsort(np.argsort(-theo_arr))  # ranks (0 = highest)
emp_rank = np.argsort(np.argsort(-emp_arr))

print(f"\n  RANKING COMPARISON:")
print(f"  {'Line':>8} {'Theo rank':>10} {'Emp rank':>10} {'Δ':>5}")
print(f"  {'-'*35}")
for i, name in enumerate(lines):
    delta = abs(int(theo_rank[i]) - int(emp_rank[i]))
    marker = " ✓" if delta <= 1 else f" ✗ (off by {delta})"
    print(f"  {name:>8} {theo_rank[i]+1:>10} {emp_rank[i]+1:>10} {delta:>5}{marker}")

# Mean rank displacement
mean_disp = np.mean(np.abs(theo_rank.astype(float) - emp_rank.astype(float)))
print(f"\n  Mean rank displacement: {mean_disp:.2f}")
print(f"  (0 = perfect match, 4.5 = random)")

if r_spearman > 0.7 and p_spearman < 0.05:
    print(f"\n  🎯 THEORETICAL P PREDICTS EMPIRICAL P")
    print(f"  The mechanism is identified: emissivity Jacobian × trapping × (1 - branch lock)")
elif r_spearman > 0.5:
    print(f"\n  ⚠️  Moderate correlation. Close but needs refinement.")
else:
    print(f"\n  ✗ Weak correlation. Model needs rethinking.")

# ============================================================
# WHAT EACH COMPONENT CONTRIBUTES
# ============================================================
print(f"\n{'='*80}")
print(f"COMPONENT ANALYSIS")
print(f"{'='*80}")

for component, getter in [
    ("Density sensitivity alone", lambda r: r['dens_sens']),
    ("Temperature sensitivity alone", lambda r: r['temp_sens']),
    ("Microstate (dens+temp)", lambda r: r['micro_sens']),
    ("Trapping alone", lambda r: r['trap']),
    ("Micro × trap", lambda r: r['micro_sens'] * r['trap']),
    ("Full model (micro × trap × (1-lock))", lambda r: r['P_raw']),
]:
    vals = np.array([getter(results[n]) for n in lines])
    # Normalize
    vals_norm = (vals - vals.min()) / (vals.max() - vals.min() + 1e-10)
    r_s, p_s = spearmanr(vals_norm, emp_arr)
    print(f"  {component:>40}: r = {r_s:+.3f} (p = {p_s:.3e})")

print(f"\n{'='*80}")
print("DONE")
print(f"{'='*80}")
