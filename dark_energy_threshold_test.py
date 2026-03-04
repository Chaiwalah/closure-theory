#!/usr/bin/env python3
"""
DARK ENERGY THRESHOLD TEST
==========================
Does the cosmic deceleration→acceleration transition predict our sigmoid threshold?

Physics:
- ΛCDM has a specific redshift where expansion switches from decelerating to accelerating
- The cosmic web structure responds to this transition with a LAG
- We test whether this lag matches our observed sigmoid z₀ ≈ 0.82 (SNe) / 1.05 (quasars)

Key computation:
1. Deceleration parameter q(z) → find z_transition where q=0
2. Cosmic web filling factor f(z) → fraction of volume in filaments vs voids
3. DM(z) properly integrated through expanding universe
4. Phase coherence length as function of web connectivity
5. Predicted sigmoid z₀ from web coherence breakdown
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import erfc
import json
import os

# ============================================================
# ΛCDM PARAMETERS (Planck 2018)
# ============================================================
H0 = 67.4        # km/s/Mpc
Omega_m = 0.315   # matter density
Omega_L = 0.685   # dark energy density
Omega_b = 0.0493  # baryon density
c = 299792.458    # km/s
f_IGM = 0.83      # fraction of baryons in IGM (Fukugita+04)
f_e = 0.875       # electron fraction (H + He ionized)
n_e0 = 2.17e-7    # mean electron density today, cm⁻³

# Derived
DH = c / H0  # Hubble distance in Mpc

print("=" * 70)
print("DARK ENERGY THRESHOLD TEST")
print("Does cosmic web tearing predict our sigmoid z₀?")
print("=" * 70)

# ============================================================
# TEST 1: Deceleration→Acceleration Transition
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: DECELERATION → ACCELERATION TRANSITION")
print("=" * 70)

def H_over_H0(z):
    """H(z)/H₀ for flat ΛCDM"""
    return np.sqrt(Omega_m * (1+z)**3 + Omega_L)

def q(z):
    """Deceleration parameter q(z)"""
    Ez2 = Omega_m * (1+z)**3 + Omega_L
    return Omega_m * (1+z)**3 / (2 * Ez2) - Omega_L / Ez2

# Find z where q(z) = 0
z_decel = brentq(q, 0.01, 2.0)
print(f"\nDeceleration→acceleration transition: z_t = {z_decel:.4f}")
print(f"  q(z<{z_decel:.2f}) < 0: accelerating (dark energy wins)")
print(f"  q(z>{z_decel:.2f}) > 0: decelerating (gravity wins)")

# Compute q at key redshifts
z_check = [0.0, 0.3, 0.5, z_decel, 0.82, 1.05, 1.5, 2.0]
print(f"\n{'z':>6} {'q(z)':>10} {'Phase':>20}")
print("-" * 40)
for z in z_check:
    phase = "ACCELERATING" if q(z) < 0 else "DECELERATING"
    marker = " ← TRANSITION" if abs(z - z_decel) < 0.01 else ""
    marker = " ← SNe z₀" if abs(z - 0.82) < 0.01 else marker
    marker = " ← QSO z₀" if abs(z - 1.05) < 0.01 else marker
    print(f"{z:6.3f} {q(z):10.4f} {phase:>20}{marker}")

# ============================================================
# TEST 2: PROPER DM(z) INTEGRATION
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: DISPERSION MEASURE vs REDSHIFT")
print("=" * 70)

def dDM_dz(z):
    """DM contribution per unit redshift (pc/cm³ per Δz)"""
    # DM_IGM = ∫ n_e(z) × c × (1+z) / H(z) dz
    # n_e(z) = n_e0 × (1+z)³ × f_IGM × f_e
    n_e_z = n_e0 * (1+z)**3 * f_IGM * f_e
    Hz = H0 * H_over_H0(z)  # km/s/Mpc
    # Convert: c/H in Mpc, n_e in cm⁻³, need pc/cm³
    # c/H(z) in Mpc, × 1e6 pc/Mpc × n_e in cm⁻³ = pc/cm³
    dl_dz = c * (1+z) / Hz  # Mpc per unit z (proper distance element)
    return n_e_z * dl_dz * 1e6  # pc/cm³

def DM_IGM(z_max):
    """Total IGM DM to redshift z"""
    result, _ = quad(dDM_dz, 0, z_max)
    return result

# Compute DM at key redshifts
print(f"\n{'z':>6} {'DM (pc/cm³)':>14} {'Note':>30}")
print("-" * 55)
z_range = np.linspace(0.1, 2.5, 25)
dm_values = []
for z in z_range:
    dm = DM_IGM(z)
    dm_values.append(dm)

# Find z at our observed thresholds
z_at_774 = brentq(lambda z: DM_IGM(z) - 774, 0.1, 3.0)
z_at_500 = brentq(lambda z: DM_IGM(z) - 500, 0.1, 3.0)
dm_at_decel = DM_IGM(z_decel)
dm_at_sne = DM_IGM(0.82)
dm_at_qso = DM_IGM(1.05)

key_points = [
    (z_decel, dm_at_decel, "DECEL→ACCEL TRANSITION"),
    (0.82, dm_at_sne, "SNe SIGMOID z₀"),
    (z_at_774, 774, "OUR THRESHOLD (774 pc/cm³)"),
    (1.05, dm_at_qso, "QUASAR SIGMOID z₀"),
    (z_at_500, 500, "FRB THRESHOLD (500 pc/cm³)"),
]

for z, dm, note in sorted(key_points):
    print(f"{z:6.3f} {dm:14.1f} {note:>30}")

print(f"\n*** DM at deceleration transition: {dm_at_decel:.1f} pc/cm³ ***")
print(f"*** Our SNe threshold DM: 774 pc/cm³ (at z={z_at_774:.3f}) ***")
print(f"*** GAP: {774 - dm_at_decel:.1f} pc/cm³ = {(z_at_774-z_decel)/z_decel*100:.1f}% redshift lag ***")

# ============================================================
# TEST 3: COSMIC WEB FILLING FACTOR
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: COSMIC WEB CONNECTIVITY MODEL")
print("=" * 70)

def web_filling_factor(z):
    """
    Fraction of IGM volume in filaments (vs voids).
    
    Based on Press-Schechter / Sheth-Tormen collapse fraction.
    At high z, most matter is diffuse → high connectivity.
    At low z, matter has collapsed into halos/filaments → lower connectivity
    but denser structures.
    
    The CONNECTIVITY (not density) is what matters for phase coherence.
    A connected web = coherent lens network.
    A fragmented web = intermittent, lossy.
    """
    # Linear growth factor D(z) ∝ H(z) × ∫₀ᵃ da/(aH)³
    # Approximate: D(z) ≈ (1+z)⁻¹ × [Ω_m(z)]^0.55
    Omega_m_z = Omega_m * (1+z)**3 / (Omega_m * (1+z)**3 + Omega_L)
    D_z = 1/(1+z) * Omega_m_z**0.55 / (Omega_m**0.55)  # normalized to D(0)=1
    
    # Collapse fraction: fraction of matter in structures > M*
    # Using simplified Press-Schechter
    sigma_8 = 0.811  # Planck 2018
    delta_c = 1.686  # critical overdensity for collapse
    
    # sigma(M*, z) = sigma_8 × D(z) for M* scale
    sigma_z = sigma_8 * D_z
    
    # Collapsed fraction
    f_coll = erfc(delta_c / (np.sqrt(2) * sigma_z))
    
    return f_coll

def web_connectivity(z):
    """
    Phase coherence of the cosmic web lens network.
    
    High connectivity → coherent refractive transfer (information preserved)
    Low connectivity → intermittent → lossy (information destroyed)
    
    Connectivity depends on:
    1. Filling factor (how much volume is in filaments)
    2. Mean free path between filaments
    3. Whether filaments form a percolating network
    
    Percolation threshold: ~16% volume filling for random networks.
    Below this → isolated islands. Above → connected web.
    """
    f = web_filling_factor(z)
    
    # Percolation threshold
    f_perc = 0.16
    
    # Connectivity is a sigmoid around the percolation threshold
    # Above threshold: connected web, coherent transfer
    # Below threshold: isolated islands, intermittent
    sharpness = 20  # transition width
    connectivity = 1 / (1 + np.exp(-sharpness * (f - f_perc)))
    
    return connectivity, f

print(f"\n{'z':>6} {'f_coll':>8} {'f_fill':>8} {'Connected':>12} {'Phase':>20}")
print("-" * 60)

z_test = np.linspace(0.0, 2.5, 50)
connectivities = []
z_vals = []
for z in z_test:
    conn, f_fill = web_connectivity(z)
    connectivities.append(conn)
    z_vals.append(z)

# Find where connectivity drops to 50% (transition)
conn_array = np.array(connectivities)
z_array = np.array(z_vals)

# Print key redshifts
for z in [0.0, 0.3, z_decel, 0.82, 1.05, 1.5, 2.0, 2.5]:
    conn, f_fill = web_connectivity(z)
    phase = "CONNECTED" if conn > 0.5 else "FRAGMENTED"
    f_c = web_filling_factor(z)
    print(f"{z:6.3f} {f_c:8.4f} {f_fill:8.4f} {conn:12.4f} {phase:>20}")

# ============================================================
# TEST 4: PREDICTED SIGMOID FROM WEB COHERENCE
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: PREDICTED vs OBSERVED SIGMOID")
print("=" * 70)

def coherence_degradation(z):
    """
    Information degradation as function of redshift.
    
    Model: cumulative loss of phase coherence through
    the cosmic web lens network.
    
    D(z) = 1 - exp(-∫₀ᶻ κ(z') × (1-C(z')) dz')
    
    where κ(z) = opacity to coherent transfer
          C(z) = web connectivity (coherent → lossless)
          (1-C) = fragmented fraction (lossy)
    """
    def integrand(zp):
        conn, _ = web_connectivity(zp)
        # Base opacity from electron density
        n_e = n_e0 * (1+zp)**3
        # Normalize so integral gives reasonable values
        kappa = n_e / n_e0  # relative to today
        # Loss rate proportional to (1 - connectivity)
        return kappa * (1 - conn) / H_over_H0(zp)
    
    result, _ = quad(integrand, 0, z)
    return 1 - np.exp(-0.5 * result)

# Compute degradation curve
z_fine = np.linspace(0.01, 2.5, 200)
degradation = [coherence_degradation(z) for z in z_fine]
degradation = np.array(degradation)

# Fit sigmoid to find predicted z₀
# D(z) ≈ 1/(1 + exp(-k(z - z₀)))
from scipy.optimize import curve_fit

def sigmoid(z, z0, k, A):
    return A / (1 + np.exp(-k * (z - z0)))

try:
    popt, pcov = curve_fit(sigmoid, z_fine, degradation, p0=[0.8, 5, 1.0], maxfev=10000)
    z0_pred, k_pred, A_pred = popt
    
    print(f"\nPREDICTED sigmoid z₀ = {z0_pred:.4f} (steepness k = {k_pred:.2f})")
    print(f"OBSERVED  sigmoid z₀ = 0.82 (SNe Ia)")
    print(f"OBSERVED  sigmoid z₀ = 1.05 (Quasars)")
    print(f"\nΔz (predicted - SNe):    {z0_pred - 0.82:+.4f}")
    print(f"Δz (predicted - QSO):    {z0_pred - 1.05:+.4f}")
    print(f"Δz (predicted - decel):  {z0_pred - z_decel:+.4f}")
    
    # How close is the prediction?
    if abs(z0_pred - 0.82) < 0.15:
        verdict_sne = "UNCOMFORTABLY CLOSE"
    elif abs(z0_pred - 0.82) < 0.3:
        verdict_sne = "IN THE NEIGHBORHOOD"
    else:
        verdict_sne = "DOESN'T MATCH"
    
    print(f"\nSNe verdict: {verdict_sne}")
    
except Exception as e:
    print(f"Sigmoid fit failed: {e}")
    z0_pred = None

# ============================================================
# TEST 5: SENSITIVITY TO COSMOLOGICAL PARAMETERS
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: COSMOLOGICAL PARAMETER SENSITIVITY")
print("=" * 70)
print("If z₀ tracks the dark energy transition, different cosmologies")
print("should predict different z₀. This makes z₀ a NEW PROBE.")

cosmologies = [
    ("Planck 2018 (standard)", 0.315, 0.685),
    ("High Ω_m (0.35)", 0.35, 0.65),
    ("Low Ω_m (0.28)", 0.28, 0.72),
    ("Einstein-de Sitter (no DE)", 1.0, 0.0),
    ("De Sitter (pure DE)", 0.01, 0.99),
    ("DESI 2024 hint (w≠-1)", 0.295, 0.705),
]

print(f"\n{'Cosmology':>35} {'Ω_m':>6} {'Ω_Λ':>6} {'z_decel':>8} {'DM_decel':>10}")
print("-" * 70)

for name, om, ol in cosmologies:
    if om >= 1.0:
        # No acceleration in EdS
        print(f"{name:>35} {om:6.3f} {ol:6.3f} {'N/A':>8} {'N/A':>10}")
        continue
    if om < 0.02:
        z_d = 0.0
    else:
        # q(z) = 0 → Ω_m(1+z)³ = 2Ω_Λ → z = (2Ω_Λ/Ω_m)^(1/3) - 1
        z_d = (2*ol/om)**(1/3) - 1
    
    # DM at that z (approximate: DM ∝ z for low z, more complex for high z)
    # Use proper integral with these parameters
    def dDM_custom(z, om=om, ol=ol):
        Hz = H0 * np.sqrt(om * (1+z)**3 + ol)
        n_e_z = n_e0 * (1+z)**3 * f_IGM * f_e
        return n_e_z * c * (1+z) / Hz * 1e6
    
    dm_d, _ = quad(dDM_custom, 0, z_d)
    print(f"{name:>35} {om:6.3f} {ol:6.3f} {z_d:8.3f} {dm_d:10.1f}")

# ============================================================
# TEST 6: THE MONEY PLOT — TIMELINE
# ============================================================
print("\n" + "=" * 70)
print("TEST 6: COSMIC TIMELINE WITH THRESHOLDS")
print("=" * 70)

from scipy.integrate import quad as squad

def lookback_time(z):
    """Lookback time in Gyr"""
    def integrand(zp):
        return 1 / ((1+zp) * H_over_H0(zp))
    result, _ = squad(integrand, 0, z)
    # Convert: 1/H₀ in Gyr
    H0_per_Gyr = H0 * 3.24e-20 * 3.156e16  # H₀ in Gyr⁻¹
    return result / H0_per_Gyr

events = [
    (0.0, "NOW"),
    (z_decel, f"DECEL→ACCEL TRANSITION"),
    (0.82, "SNe SIGMOID z₀ (OBSERVED)"),
    (1.05, "QUASAR SIGMOID z₀ (OBSERVED)"),
    (1.65, "ANGULAR DIAMETER TURNOVER"),
    (2.0, "PEAK STAR FORMATION"),
]

if z0_pred is not None:
    events.append((z0_pred, "PREDICTED THRESHOLD (this test)"))

events.sort(key=lambda x: x[0])

print(f"\n{'z':>6} {'Lookback':>10} {'DM':>10} {'Event':>40}")
print("-" * 70)
for z, name in events:
    t = lookback_time(z)
    dm = DM_IGM(z)
    print(f"{z:6.3f} {t:8.2f} Gyr {dm:8.1f} {name:>40}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
DECELERATION→ACCELERATION TRANSITION:
  z_transition = {z_decel:.4f}
  DM at transition = {dm_at_decel:.1f} pc/cm³
  Lookback time = {lookback_time(z_decel):.2f} Gyr

OUR OBSERVED THRESHOLDS:
  SNe sigmoid z₀ = 0.82 → DM = {dm_at_sne:.1f} pc/cm³
  QSO sigmoid z₀ = 1.05 → DM = {dm_at_qso:.1f} pc/cm³
  FRB threshold  ~ 500 pc/cm³ → z ≈ {z_at_500:.3f}

PREDICTED THRESHOLD (cosmic web coherence model):
  z₀ = {z0_pred:.4f} (DM = {DM_IGM(z0_pred):.1f} pc/cm³)

KEY RESULT:
  The deceleration transition at z={z_decel:.3f} PRECEDES our SNe threshold 
  at z=0.82 by Δz={0.82-z_decel:.3f} ({(0.82-z_decel)/z_decel*100:.1f}% lag).
  
  This lag is consistent with the cosmic web RESPONDING to the expansion 
  change — structure doesn't tear instantly when dark energy starts winning.
  It takes ~{lookback_time(0.82)-lookback_time(z_decel):.2f} Gyr for the web to fragment 
  enough that spectral coherence breaks down.
  
  The different z₀ for SNe (0.82) vs quasars (1.05) could reflect different 
  SENSITIVITY to web connectivity — SNe probe broadband color (sensitive to 
  large-scale coherence) while quasars probe emission line structure 
  (sensitive to smaller-scale coherence that persists longer).
""")

if z0_pred and abs(z0_pred - 0.82) < 0.2:
    print("⚠️  THE MODEL THRESHOLD LANDS WITHIN 0.2 OF THE SNe OBSERVATION.")
    print("⚠️  This is either a coincidence or the dark energy transition")
    print("⚠️  is physically connected to spectral information degradation.")
    print("⚠️  UNCOMFORTABLY COMFORTABLE.")
elif z0_pred:
    print(f"Model prediction z₀={z0_pred:.3f} vs observed 0.82.")
    print("Further refinement needed.")

# Save results
results = {
    'z_deceleration': float(z_decel),
    'dm_at_deceleration': float(dm_at_decel),
    'z0_predicted': float(z0_pred) if z0_pred else None,
    'z0_sne_observed': 0.82,
    'z0_qso_observed': 1.05,
    'z_at_774_dm': float(z_at_774),
    'z_at_500_dm': float(z_at_500),
    'lag_z': float(0.82 - z_decel),
    'lag_percent': float((0.82 - z_decel) / z_decel * 100),
    'lag_gyr': float(lookback_time(0.82) - lookback_time(z_decel)),
}

os.makedirs('results_dark_energy', exist_ok=True)
with open('results_dark_energy/threshold_test.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nResults saved to results_dark_energy/threshold_test.json")
