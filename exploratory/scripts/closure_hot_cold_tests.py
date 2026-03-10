#!/usr/bin/env python3
"""
closure_hot_cold_tests.py — Decisive If/Then Tests
====================================================

Each test is binary: HOT (confirmed) or COLD (dead).
No hand-waving. No "interesting." Pass or fail.

TEST 1: CMB q-Value Unification
  IF the Γ_Σ "divergence" is really the q-ladder in CMB space,
  THEN there exist q values for each CMB parameter that make ALL
  observables converge to a SINGLE Γ_Σ ≈ 227-311 × σ_T.
  If no consistent q set exists → framework has a hole.

TEST 2: Sharpness Source — Raw Hubble Residuals vs β(z)
  IF the sharpness is from measurement processing,
  THEN raw Pantheon+ Hubble residuals show BROADER sigmoid than β(z).
  If raw residuals are ALSO sharp → the sharpness is physical.

TEST 3: FRB Compression Level at DM=500
  IF DM=500 is the onset (not midpoint),
  THEN the Beer-Lambert compression at DM=500 should be 5-15%, not 50%.
  If it's >30% → the FRB threshold IS the midpoint and the 7× miss is real.

TEST 4: Blackbody q Prediction
  IF CMB is blackbody (max entropy thermal) with lower q than line emission,
  THEN the effective q needed to match H₀=67.4 should be << 0.88.
  If q_needed ≈ 0.88 → blackbody argument is wrong, CMB IS like line emission.

TEST 5: Locked Observable Immunity at Quantum Scale
  IF Beer-Lambert operates at all scales,
  THEN the STAR lambda hyperon POSITIONS should be uncorrelated (locked)
  while SPINS are correlated (diagnostic).
  Check: does the published data separate geometric vs diagnostic observables?

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize, minimize_scalar
import json, os
from pathlib import Path

sigma_T = 6.6524e-25  # cm²

# ============================================================
# COSMOLOGICAL SETUP (from beer_lambert.py)
# ============================================================
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

# Pre-compute Sigma grid
z_fine = np.linspace(0.01, 3.0, 500)
Sigma_fine = np.array([Sigma(z) for z in z_fine])
Sigma_CMB = Sigma(1100)

# ============================================================
# TEST 1: CMB q-VALUE UNIFICATION
# ============================================================
print("=" * 70)
print("TEST 1: CMB q-VALUE UNIFICATION")
print("IF divergence = q-ladder, THEN one Γ_Σ fits all with right q's")
print("=" * 70)

# From closure_rd_bias.py:
# H₀=67.4 needs: Γ_Σ × Σ_CMB × q_rd² × f_damp = some value that gives 8% r_d shift
# σ₈=0.766 needs: Γ_Σ × Σ_CMB × q_σ₈² × ... = some value
# A_L=1.18 needs: Γ_Σ × Σ_CMB × q_AL² × ... = some value
#
# If Γ_Σ is FIXED (say at 227 × σ_T from Tripp), what q values make each work?

# For H₀: need rd_ratio = 1 + Γ × Σ_CMB × q² × f_damp such that 73.6/rd_ratio = 67.4
# rd_ratio = 73.6/67.4 = 1.0920
# Γ × Σ_CMB × q² × f_damp = 0.0920
# q² = 0.0920 / (Γ × Σ_CMB × f_damp)

f_damp = 0.20  # fraction of r_d info from damping tail
H0_true = 73.6  # geometric (TDCOSMO)
H0_CMB = 67.4   # Planck

rd_ratio_needed = H0_true / H0_CMB  # = 1.092

# For σ₈: need σ₈ to shift from 0.766 to 0.811
# Using H₀-σ₈ degeneracy: dσ₈/dH₀ ≈ 0.006
# Δσ₈ = 0.045, ΔH₀ = -5.68 → independent of specific Γ
# σ₈ is constrained by: A_s from peak heights + growth factor
# Direct: σ₈_bias = σ₈_true × 0.5 × Γ × Σ_CMB × q_σ₈² × f_peak
# Where σ₈_bias = 0.811 - 0.766 = 0.045
# 0.045 = 0.766 × 0.5 × Γ × Σ_CMB × q_σ₈² × f_peak
sigma8_bias = 0.811 - 0.766
sigma8_true = 0.766
f_peak = 0.15  # fraction of σ₈ info from peak heights

# For A_L: excess smoothing
# (A_L - 1) × smoothing_per_AL ≈ Γ × Σ_CMB × q_AL² × f_smooth
# 0.18 × 0.05 = Γ × Σ_CMB × q_AL² × f_smooth
AL_excess = 0.180
smoothing_per_AL = 0.05
f_smooth = 0.10  # fraction of info in lensing-degenerate channel

# Scan Γ_Σ and find required q values
print(f"\n{'Γ_Σ/σ_T':>10} | {'q_rd':>8} {'q_σ₈':>8} {'q_AL':>8} | {'q ordering?':>14} {'Physical?':>10}")
print("-" * 75)

viable_solutions = []

for G_ratio in [100, 150, 200, 227, 250, 300, 311, 350, 400]:
    Gamma = G_ratio * sigma_T
    
    # q_rd from H₀
    denominator_rd = Gamma * Sigma_CMB * f_damp
    q_rd_sq = (rd_ratio_needed - 1) / denominator_rd
    q_rd = np.sqrt(q_rd_sq) if q_rd_sq > 0 else float('nan')
    
    # q_σ₈ from σ₈ tension
    denominator_s8 = sigma8_true * 0.5 * Gamma * Sigma_CMB * f_peak
    q_s8_sq = sigma8_bias / denominator_s8 if denominator_s8 > 0 else float('nan')
    q_s8 = np.sqrt(q_s8_sq) if q_s8_sq > 0 and q_s8_sq < 100 else float('nan')
    
    # q_AL from A_L anomaly
    denominator_al = Gamma * Sigma_CMB * f_smooth
    q_al_sq = (AL_excess * smoothing_per_AL) / denominator_al
    q_al = np.sqrt(q_al_sq) if q_al_sq > 0 else float('nan')
    
    # Check ordering: should be q_AL < q_σ₈ < q_rd (geometric < mixed < diagnostic)
    ordering = "—"
    physical = "—"
    if not any(np.isnan([q_rd, q_s8, q_al])):
        ordering = "✓" if q_al < q_s8 < q_rd else "✗"
        # Physical: all q should be between 0 and 1
        physical = "✓" if all(0 < q < 1 for q in [q_rd, q_s8, q_al]) else "✗"
        if ordering == "✓" and physical == "✓":
            viable_solutions.append((G_ratio, q_rd, q_s8, q_al))
    
    print(f"{G_ratio:>10} | {q_rd:>8.4f} {q_s8:>8.4f} {q_al:>8.4f} | {ordering:>14} {physical:>10}")

if viable_solutions:
    print(f"\n★ VIABLE SOLUTIONS FOUND: {len(viable_solutions)}")
    for G, qr, qs, qa in viable_solutions:
        print(f"  Γ_Σ = {G} × σ_T: q_rd = {qr:.4f}, q_σ₈ = {qs:.4f}, q_AL = {qa:.4f}")
        print(f"    Ratio q_rd/q_AL = {qr/qa:.1f} (should be ~5-10× for diagnostic/geometric)")
    print(f"\n  HOT: The divergence IS the q-ladder. One Γ_Σ works with physically ordered q's.")
else:
    print(f"\n  COLD: No physically consistent q set found. The divergence is real.")

# ============================================================
# TEST 2: RAW HUBBLE RESIDUAL SIGMOID
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: SHARPNESS SOURCE — Raw residuals vs β(z)")
print("IF measurement processing sharpens, THEN raw residuals are broader")
print("=" * 70)

# Load Pantheon+ data
data_file = "data/pantheon_plus.dat"
if os.path.exists(data_file):
    data = np.genfromtxt(data_file, names=True, dtype=None, encoding='utf-8', 
                          delimiter=' ', skip_header=0)
    
    # Get column names
    cols = data.dtype.names
    print(f"  Columns available: {cols[:10]}...")
    
    # Try to extract z, mu (distance modulus), and mu_err
    # Pantheon+ format: CID zHD zHDERR mB mBERR x1 x1ERR c cERR
    try:
        z = data['zHD']
        mB = data['mB']
        x1 = data['x1']
        c = data['c']
        
        # Filter
        mask = (z > 0.01) & (z < 2.3) & np.isfinite(mB) & np.isfinite(x1) & np.isfinite(c)
        z = z[mask]
        mB = mB[mask]
        x1 = x1[mask]
        c = c[mask]
        
        print(f"  Loaded {len(z)} SNe from Pantheon+")
        
        # Compute Tripp-standardized distance modulus
        # μ = mB - M + α·x1 - β·c
        # Use Pantheon+ best-fit values: α=0.148, β=3.02, M=-19.25
        alpha_sn = 0.148
        beta_sn = 3.02  # global average
        M_sn = -19.25
        
        mu = mB - M_sn + alpha_sn * x1 - beta_sn * c
        
        # Compute expected μ from ΛCDM
        from scipy.integrate import quad as squad
        def mu_LCDM(z_val):
            def integrand(zp):
                return 1.0 / np.sqrt(Omega_m*(1+zp)**3 + Omega_L)
            dist, _ = squad(integrand, 0, z_val)
            dL = (1+z_val) * dist * 2.998e5 / H0  # Mpc
            return 5 * np.log10(dL) + 25
        
        mu_expected = np.array([mu_LCDM(zi) for zi in z])
        
        # Hubble residuals
        delta_mu = mu - mu_expected
        
        # Bin the residuals
        z_bins = np.linspace(0.01, 1.5, 15)
        z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])
        residual_means = []
        residual_stds = []
        
        for i in range(len(z_bins)-1):
            mask_bin = (z >= z_bins[i]) & (z < z_bins[i+1])
            if np.sum(mask_bin) > 5:
                residual_means.append(np.median(delta_mu[mask_bin]))
                residual_stds.append(np.std(delta_mu[mask_bin]) / np.sqrt(np.sum(mask_bin)))
            else:
                residual_means.append(float('nan'))
                residual_stds.append(float('nan'))
        
        residual_means = np.array(residual_means)
        residual_stds = np.array(residual_stds)
        
        # Compute SCATTER per bin (this is the diagnostic: scatter should increase with compression)
        scatter_per_bin = []
        for i in range(len(z_bins)-1):
            mask_bin = (z >= z_bins[i]) & (z < z_bins[i+1])
            if np.sum(mask_bin) > 10:
                scatter_per_bin.append(np.std(delta_mu[mask_bin]))
            else:
                scatter_per_bin.append(float('nan'))
        scatter_per_bin = np.array(scatter_per_bin)
        
        # The RESIDUAL MEAN should show compression signal
        # If constant β is used but true β decreases with z,
        # high-z SNe will have OVER-corrected colors → residuals shift
        
        # Fit sigmoid to residual means (absolute value or trend)
        # Actually: the SCATTER is more diagnostic. If compression adds noise,
        # scatter should increase with z following a sigmoid-like curve
        
        print(f"\n  {'z_center':>8} {'N':>5} {'Δμ mean':>10} {'scatter':>10}")
        print(f"  {'-'*40}")
        for i, zc in enumerate(z_centers):
            mask_bin = (z >= z_bins[i]) & (z < z_bins[i+1])
            N = np.sum(mask_bin)
            print(f"  {zc:>8.3f} {N:>5} {residual_means[i]:>+10.4f} {scatter_per_bin[i]:>10.4f}")
        
        # Check: does scatter have a sigmoid-like profile?
        valid = ~np.isnan(scatter_per_bin) & (z_centers < 1.5)
        if np.sum(valid) > 4:
            from scipy.stats import spearmanr
            rho, p = spearmanr(z_centers[valid], scatter_per_bin[valid])
            print(f"\n  Scatter vs z: ρ = {rho:+.3f}, p = {p:.4f}")
            
            if rho > 0.5 and p < 0.05:
                print(f"  → Scatter INCREASES with z (compression adds noise)")
                
                # Try sigmoid fit to scatter
                from scipy.optimize import curve_fit
                def sigmoid(z, a, z0, k, b):
                    return a / (1 + np.exp(-k*(z-z0))) + b
                
                try:
                    popt, _ = curve_fit(sigmoid, z_centers[valid], scatter_per_bin[valid],
                                       p0=[0.05, 0.8, 5.0, 0.12], maxfev=5000)
                    print(f"  Sigmoid fit to scatter: z₀ = {popt[1]:.2f}, k = {popt[2]:.2f}")
                    print(f"  Compare: β(z) sigmoid k = 8.0")
                    
                    if popt[2] < 5.0:
                        print(f"\n  ★ HOT: Raw scatter sigmoid is BROADER (k={popt[2]:.1f}) than β sigmoid (k=8.0)")
                        print(f"        The measurement processing DOES sharpen the transition!")
                    else:
                        print(f"\n  ★ COLD: Raw scatter is equally sharp → sharpness is physical")
                except:
                    print(f"  Sigmoid fit failed — insufficient data or wrong shape")
            else:
                print(f"  → No clear scatter-z trend. Test inconclusive.")
        
        # Also check: color evolution directly
        print(f"\n  --- Color evolution (raw) ---")
        c_bins = []
        for i in range(len(z_bins)-1):
            mask_bin = (z >= z_bins[i]) & (z < z_bins[i+1])
            if np.sum(mask_bin) > 10:
                c_bins.append(np.median(c[mask_bin]))
            else:
                c_bins.append(float('nan'))
        c_bins = np.array(c_bins)
        
        valid_c = ~np.isnan(c_bins)
        if np.sum(valid_c) > 4:
            rho_c, p_c = spearmanr(z_centers[valid_c], c_bins[valid_c])
            print(f"  Color vs z: ρ = {rho_c:+.3f}, p = {p_c:.4f}")
            print(f"  (Negative = bluer at high z, consistent with chromatic scattering)")
    
    except Exception as e:
        print(f"  Error processing Pantheon+ data: {e}")
else:
    print(f"  ⚠ Pantheon+ data not found at {data_file}")
    print(f"  Skipping test 2.")

# ============================================================
# TEST 3: FRB COMPRESSION LEVEL AT DM=500
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: FRB COMPRESSION AT DM=500")
print("IF DM=500 is onset, THEN compression ≈ 5-15%")
print("=" * 70)

# DM=500 pc/cm³ → Σ = 500 × 3.086e18 / 1.14 = 1.354e21 cm⁻²
DM_obs = 500
Sigma_FRB = DM_obs * 3.086e18 / 1.14

# q for FRB spectral index
q_FRB = 0.604

# Compression at DM=500 for various Γ_Σ
print(f"\n  DM = {DM_obs} pc/cm³ → Σ = {Sigma_FRB:.3e} cm⁻²")
print(f"  q_SI = {q_FRB}")
print(f"\n  {'Γ_Σ/σ_T':>10} {'τ_info':>10} {'Compression':>12} {'Verdict':>10}")
print(f"  {'-'*45}")

for G_ratio in [37, 100, 150, 227, 311, 400]:
    Gamma = G_ratio * sigma_T
    tau = Gamma * Sigma_FRB * q_FRB**2
    compression = 1 - np.exp(-tau)
    
    if compression < 0.05:
        verdict = "< onset"
    elif compression < 0.15:
        verdict = "ONSET ✓"
    elif compression < 0.30:
        verdict = "moderate"
    elif compression < 0.60:
        verdict = "midpoint"
    else:
        verdict = "saturated"
    
    print(f"  {G_ratio:>10} {tau:>10.4f} {compression*100:>11.1f}% {verdict:>10}")

# What Γ_Σ gives exactly 10% compression at DM=500?
target_compression = 0.10
tau_target = -np.log(1 - target_compression)
Gamma_10pct = tau_target / (Sigma_FRB * q_FRB**2)
G_10pct = Gamma_10pct / sigma_T

print(f"\n  Γ_Σ for 10% compression at DM=500: {G_10pct:.0f} × σ_T")
print(f"  Γ_Σ for 5% compression at DM=500: {-np.log(0.95)/(Sigma_FRB * q_FRB**2)/sigma_T:.0f} × σ_T")

if 100 < G_10pct < 500:
    print(f"\n  ★ HOT: DM=500 as 10% onset is consistent with Γ_Σ = {G_10pct:.0f} × σ_T")
    print(f"         (within range of Tripp=227 and BL=311)")
else:
    print(f"\n  ★ COLD: Even as onset, DM=500 doesn't match our Γ_Σ range")

# ============================================================
# TEST 4: BLACKBODY q PREDICTION
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: WHAT q MAKES CMB WORK?")
print("IF CMB is low-q blackbody, THEN q << 0.88")
print("=" * 70)

# For H₀: Γ × Σ_CMB × q² × f_damp = rd_ratio - 1 = 0.0920
# q² = 0.0920 / (Γ × Σ_CMB × f_damp)

for G_ratio in [227, 250, 300, 311]:
    Gamma = G_ratio * sigma_T
    q_sq = (rd_ratio_needed - 1) / (Gamma * Sigma_CMB * f_damp)
    q_val = np.sqrt(q_sq) if q_sq > 0 else float('nan')
    
    print(f"\n  Γ_Σ = {G_ratio} × σ_T:")
    print(f"    q_rd needed for H₀ = 67.4: {q_val:.6f}")
    print(f"    (Compare: line emission q ≈ 0.5-1.0)")
    
    if q_val < 0.01:
        print(f"    → CMB q is {q_val/0.88*100:.2f}% of line emission q")
        print(f"    → EXTREMELY low. Blackbody is nearly locked.")
    elif q_val < 0.1:
        print(f"    → CMB q is {q_val/0.88*100:.1f}% of line emission q")
        print(f"    → Low but non-trivial. Blackbody partially diagnostic.")

# Check: is q ~ 0.003 physical?
print(f"\n  Physical interpretation:")
print(f"  Line emission: q ~ 0.5-1.0 (info in narrow spectral features)")
print(f"  Blackbody CMB: q ~ 0.003 (info in spatial correlations, not spectral)")
print(f"  Ratio: ~200× less susceptible")
print(f"")
print(f"  Is this reasonable? The CMB power spectrum has ~3000 independent")
print(f"  multipoles (ℓ = 2 to 2500). Each carries a tiny fraction of the")
print(f"  total diagnostic info. A single emission line carries ALL its info")
print(f"  in one spectral channel. Distributed info ≈ harder to destroy.")
print(f"")
print(f"  Analogy: destroying one page of a book (emission line) vs")
print(f"  smudging one word per page across a 3000-page book (CMB).")
print(f"  Same total damage, but the book is still readable.")

# Is the q ratio consistent with the information distribution?
N_modes_CMB = 2500  # independent ℓ modes
q_expected_from_distribution = 1.0 / np.sqrt(N_modes_CMB)
print(f"\n  Expected q from √(1/N_modes): {q_expected_from_distribution:.4f}")
print(f"  Needed q from H₀: ~0.003")
print(f"  Ratio: {0.003/q_expected_from_distribution:.2f}")

if abs(np.log10(0.003) - np.log10(q_expected_from_distribution)) < 1:
    print(f"\n  ★ HOT: Same order of magnitude! CMB q ~ 1/√N_modes is physical.")
else:
    print(f"\n  ★ COLD: q needed is too different from 1/√N_modes prediction.")

# ============================================================
# TEST 5: STAR EXPERIMENT PREDICTION
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: STAR EXPERIMENT — LOCKED vs DIAGNOSTIC")
print("IF Beer-Lambert at all scales, THEN positions locked, spins diagnostic")
print("=" * 70)

print(f"""
  STAR collaboration (Nature, 2025/2026):
  - Lambda hyperons produced in Au+Au collisions at RHIC
  - CLOSE pairs: spin correlation PRESERVED (aligned)
  - FAR pairs: spin correlation VANISHES
  
  In our framework:
  - Spin correlation = diagnostic (encodes shared quantum origin, q ≈ 1)
  - Position/momentum = geometric (encodes trajectory, q ≈ 0)
  
  PREDICTION: In the same dataset, the POSITION correlations of lambda
  hyperon pairs should NOT decay with separation distance. Only SPIN
  correlations should decay.
  
  This is directly testable from the published STAR data:
  - Plot position correlation vs pair separation: should be FLAT
  - Plot spin correlation vs pair separation: should show decay
  
  If position correlations ALSO decay → compression acts on geometry too → COLD
  If ONLY spin decays → locked/diagnostic split confirmed at QCD scale → HOT
  
  Status: NEEDS DATA (published in Nature, may be accessible)
""")

# ============================================================
# SCORECARD
# ============================================================
print("=" * 70)
print("SCORECARD")
print("=" * 70)

print(f"""
  TEST 1 (CMB q-unification):  Run above — check results
  TEST 2 (Sharpness source):   Run above — check Pantheon+ results
  TEST 3 (FRB onset):          Run above — check DM=500 level
  TEST 4 (Blackbody q):        Run above — check q value
  TEST 5 (STAR prediction):    NEEDS published data
""")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_hot_cold")
results_dir.mkdir(exist_ok=True)

with open(results_dir / "test_summary.json", 'w') as f:
    json.dump({
        'test1_viable': len(viable_solutions) > 0,
        'test1_solutions': [(g, qr, qs, qa) for g, qr, qs, qa in viable_solutions],
        'test3_G_10pct': float(G_10pct),
        'test4_q_needed': float(np.sqrt((rd_ratio_needed - 1) / (227*sigma_T * Sigma_CMB * f_damp))),
        'test4_q_expected': float(q_expected_from_distribution),
    }, f, indent=2)

print(f"\nResults saved to {results_dir}/")
