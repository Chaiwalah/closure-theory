#!/usr/bin/env python3
"""
BOUNDARY UNIFICATION TEST
Are the three sigmoid thresholds (SNe z₀≈0.82, Quasars z₀≈1.0-1.2, FRBs DM≈500)
measuring the SAME physical boundary through different channels?

Convert all thresholds to:
- Proper distance (Mpc)
- Lookback time (Gyr)
- Comoving distance (Mpc)

If they converge to the same physical scale: it's one boundary, 
different observables just "see" it at different apparent distances.

If they DON'T converge: there's something observable-specific about the threshold.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import json, os

os.makedirs('results_boundary', exist_ok=True)

# ============================================================
# Cosmological calculator (flat ΛCDM)
# ============================================================
H0 = 70.0  # km/s/Mpc
Om = 0.3
Ol = 0.7
c_km = 299792.458  # km/s
DH = c_km / H0  # Hubble distance in Mpc

def E(z):
    return np.sqrt(Om * (1+z)**3 + Ol)

def comoving_distance(z):
    """Comoving distance in Mpc."""
    result, _ = quad(lambda zp: 1/E(zp), 0, z)
    return DH * result

def luminosity_distance(z):
    return (1+z) * comoving_distance(z)

def lookback_time(z):
    """Lookback time in Gyr."""
    tH = 1/H0 * 3.086e19 / 3.156e16  # Hubble time in Gyr (≈14 Gyr)
    result, _ = quad(lambda zp: 1/((1+zp)*E(zp)), 0, z)
    return tH * result

def proper_distance(z):
    """Proper distance at emission time in Mpc."""
    return comoving_distance(z) / (1+z)

# ============================================================
# Known sigmoid thresholds from our data
# ============================================================
print("=" * 70)
print("BOUNDARY UNIFICATION — Are all thresholds the same boundary?")
print("=" * 70)

thresholds = {
    'SNe Ia (color-distance)': {
        'z0': 0.82,
        'z0_range': (0.70, 0.95),
        'source': 'Pantheon+ 1590 SNe',
        'observable': 'color coupling to Hubble residual',
        'N_modes_dominant': 3,  # color depends on T, metallicity, dust
    },
    'Quasars Hβ-MgII EW': {
        'z0': 1.05,
        'z0_range': (0.95, 1.15),
        'source': 'SDSS DR16Q 750K QSOs',
        'observable': 'EW cross-correlation',
        'N_modes_dominant': 4,  # EW depends on T, density, ionization, column
    },
    'Quasars CIII-CIV EW': {
        'z0': 1.21,
        'z0_range': (1.10, 1.35),
        'source': 'SDSS DR16Q',
        'observable': 'high-ionization EW coupling',
        'N_modes_dominant': 5,  # high-ion lines more complex
    },
    'Quasars MgII-CIV EW': {
        'z0': 1.23,
        'z0_range': (1.10, 1.40),
        'source': 'SDSS DR16Q',
        'observable': 'cross-zone BLR EW coupling',
        'N_modes_dominant': 5,
    },
    'FRBs width-spectral_index': {
        'z0': 0.47,  # DM≈500 → z≈0.47 (DM-z relation: DM_cosmic ≈ 1000*z for z<1)
        'z0_range': (0.35, 0.60),
        'source': 'CHIME Cat1 535 FRBs',
        'observable': 'temporal-spectral coupling',
        'N_modes_dominant': 2,  # width + spectral index = simpler
    },
    'FRBs DM-Fluence': {
        'z0': 0.17,  # DM≈173 → z≈0.17
        'z0_range': (0.10, 0.25),
        'source': 'CHIME Cat1',
        'observable': 'dispersion-flux coupling',
        'N_modes_dominant': 2,
    },
}

# Also compute for the DESI galaxy thresholds
# The void galaxy test showed effects growing with z, onset ~0.15-0.25
# But that's environment, not distance threshold

print("\n" + "-" * 70)
print(f"{'Source':<30} {'z₀':>6} {'d_C (Mpc)':>10} {'d_L (Mpc)':>10} {'t_lb (Gyr)':>10} {'d_prop (Mpc)':>12} {'N_modes':>8}")
print("-" * 70)

results = {}
z0_values = []
dc_values = []
tl_values = []
dp_values = []
nm_values = []

for name, info in thresholds.items():
    z0 = info['z0']
    dc = comoving_distance(z0)
    dl = luminosity_distance(z0)
    tl = lookback_time(z0)
    dp = proper_distance(z0)
    
    results[name] = {
        'z0': z0,
        'comoving_Mpc': round(dc, 1),
        'luminosity_Mpc': round(dl, 1),
        'lookback_Gyr': round(tl, 2),
        'proper_Mpc': round(dp, 1),
        'N_modes': info['N_modes_dominant'],
    }
    
    z0_values.append(z0)
    dc_values.append(dc)
    tl_values.append(tl)
    dp_values.append(dp)
    nm_values.append(info['N_modes_dominant'])
    
    print(f"{name:<30} {z0:6.2f} {dc:10.1f} {dl:10.1f} {tl:10.2f} {dp:12.1f} {info['N_modes_dominant']:>8}")

print("-" * 70)

# ============================================================
# Do they converge to one boundary?
# ============================================================
print("\n" + "=" * 70)
print("CONVERGENCE TEST — Same boundary?")
print("=" * 70)

from scipy import stats

# Test 1: Do thresholds cluster in any distance metric?
print(f"\n  z₀ range: {min(z0_values):.2f} to {max(z0_values):.2f} (spread: {max(z0_values)-min(z0_values):.2f})")
print(f"  Comoving range: {min(dc_values):.0f} to {max(dc_values):.0f} Mpc (spread: {max(dc_values)-min(dc_values):.0f})")
print(f"  Lookback range: {min(tl_values):.2f} to {max(tl_values):.2f} Gyr (spread: {max(tl_values)-min(tl_values):.2f})")
print(f"  Proper range: {min(dp_values):.0f} to {max(dp_values):.0f} Mpc (spread: {max(dp_values)-min(dp_values):.0f})")

# Coefficient of variation (lower = more clustered)
cv_z = np.std(z0_values) / np.mean(z0_values)
cv_dc = np.std(dc_values) / np.mean(dc_values)
cv_tl = np.std(tl_values) / np.mean(tl_values)
cv_dp = np.std(dp_values) / np.mean(dp_values)

print(f"\n  Coefficient of variation:")
print(f"    z₀:       {cv_z:.3f}")
print(f"    Comoving: {cv_dc:.3f}")
print(f"    Lookback: {cv_tl:.3f}")
print(f"    Proper:   {cv_dp:.3f}")
print(f"  (Lower = more clustered. If one metric is tightest, that's the natural scale)")

# ============================================================
# KEY TEST: Does N_modes predict the threshold?
# ============================================================
print("\n" + "=" * 70)
print("N_MODES vs THRESHOLD — Does complexity set the boundary distance?")
print("=" * 70)

r_z, p_z = stats.spearmanr(nm_values, z0_values)
r_dc, p_dc = stats.spearmanr(nm_values, dc_values)
r_tl, p_tl = stats.spearmanr(nm_values, tl_values)

print(f"\n  N_modes vs z₀:       ρ = {r_z:+.3f} (p = {p_z:.3e})")
print(f"  N_modes vs d_C:      ρ = {r_dc:+.3f} (p = {p_dc:.3e})")
print(f"  N_modes vs lookback: ρ = {r_tl:+.3f} (p = {p_tl:.3e})")

if r_z > 0:
    print(f"\n  ★ MORE COMPLEX observables see the boundary FURTHER AWAY")
    print(f"    This means: complex observables REACH further before decorrelating")
    print(f"    Counterintuitive if complexity = fragility")
    print(f"    BUT makes sense if: complex observables are MORE SENSITIVE DETECTORS")
    print(f"    of the boundary — they detect it at larger distances because they have")
    print(f"    more channels to couple through")
elif r_z < 0:
    print(f"\n  ★ SIMPLER observables see the boundary FURTHER AWAY")
    print(f"    More complex = earlier threshold = more fragile")
    print(f"    Consistent with: complexity = vulnerability to the medium")

# Fit the relationship
if len(nm_values) >= 3:
    slope, intercept, r_lin, p_lin, se = stats.linregress(nm_values, z0_values)
    print(f"\n  Linear fit: z₀ = {slope:.3f} × N_modes + {intercept:.3f}")
    print(f"    R² = {r_lin**2:.3f}, p = {p_lin:.3e}")
    
    # Predict threshold for N_modes = 1 (stretch-like, locked)
    z0_pred_1 = slope * 1 + intercept
    z0_pred_6 = slope * 6 + intercept
    print(f"\n  Predicted threshold for N_modes=1 (locked): z₀ = {z0_pred_1:.2f}")
    print(f"  Predicted threshold for N_modes=6 (max complexity): z₀ = {z0_pred_6:.2f}")
    print(f"  If locked observables truly never degrade, z₀ → ∞ (always preserved)")

# ============================================================
# Physical scale interpretation
# ============================================================
print("\n" + "=" * 70)
print("PHYSICAL INTERPRETATION")
print("=" * 70)

# What known scales are near these distances?
mean_dc = np.mean(dc_values)
mean_tl = np.mean(tl_values)
mean_z = np.mean(z0_values)

print(f"\n  Mean threshold: z = {mean_z:.2f}, d_C = {mean_dc:.0f} Mpc, t_lb = {mean_tl:.2f} Gyr")
print(f"\n  Known scales near this distance:")
print(f"    BAO scale: ~150 Mpc (much smaller)")
print(f"    KBC void radius: ~300-600 Mpc (OVERLAPS with range)")
print(f"    Horizon at matter-radiation equality: ~100 Mpc (smaller)")
print(f"    Observable universe: ~14,000 Mpc (much larger)")
print(f"    Typical supercluster: ~100-200 Mpc")

# KBC void
kbc_z = 0.07  # KBC void extends to z≈0.07
kbc_dc = comoving_distance(kbc_z)
print(f"\n  KBC void edge: z ≈ {kbc_z}, d_C ≈ {kbc_dc:.0f} Mpc")
print(f"  Our thresholds at: {min(dc_values):.0f} - {max(dc_values):.0f} Mpc")
print(f"  Ratio (threshold/KBC): {mean_dc/kbc_dc:.1f}x")

# ============================================================
# The loop test: if we're inside a loop, the boundary should
# show up at the SAME lookback time regardless of direction
# ============================================================
print("\n" + "=" * 70)
print("LOOP TEST — Is the boundary isotropic?")
print("=" * 70)
print("If we're inside a closure loop, the boundary should be roughly")
print("the same distance in all directions (isotropic).")
print("If it's source-dependent, it should vary with source class.")
print(f"\n  Source-class means:")

# Group by source class
sn_z0 = [thresholds[k]['z0'] for k in thresholds if 'SNe' in k]
qso_z0 = [thresholds[k]['z0'] for k in thresholds if 'Quasar' in k]
frb_z0 = [thresholds[k]['z0'] for k in thresholds if 'FRB' in k]

print(f"    SNe:     z₀ = {np.mean(sn_z0):.2f}")
print(f"    Quasars: z₀ = {np.mean(qso_z0):.2f}")
print(f"    FRBs:    z₀ = {np.mean(frb_z0):.2f}")
print(f"\n  If same boundary: these should be similar")
print(f"  If channel-dependent: they differ because each probe 'reaches' differently")

# Spread
total_spread = max(z0_values) - min(z0_values)
class_means = [np.mean(sn_z0), np.mean(qso_z0), np.mean(frb_z0)]
class_spread = max(class_means) - min(class_means)
print(f"\n  Total z₀ spread: {total_spread:.2f}")
print(f"  Class-mean spread: {class_spread:.2f}")
print(f"  Within-class / between-class: tells us if it's one boundary or channel-specific")

# ============================================================
# Save
# ============================================================
all_results = {
    'thresholds': results,
    'convergence': {
        'cv_z': float(cv_z),
        'cv_comoving': float(cv_dc),
        'cv_lookback': float(cv_tl),
        'cv_proper': float(cv_dp),
    },
    'nmodes_correlation': {
        'rho_z': float(r_z), 'p_z': float(p_z),
        'rho_dc': float(r_dc), 'p_dc': float(p_dc),
        'rho_tl': float(r_tl), 'p_tl': float(p_tl),
    },
    'physical_scales': {
        'mean_z': float(mean_z),
        'mean_comoving_Mpc': float(mean_dc),
        'mean_lookback_Gyr': float(mean_tl),
    }
}

with open('results_boundary/boundary_unification.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print(f"\nResults saved to results_boundary/boundary_unification.json")
