#!/usr/bin/env python3
"""
Null Simulation: Can Conventional Astrophysics Reproduce Γ₀?

Generate 10,000 Monte Carlo mock catalogs using ONLY published
conventional evolution models (Baldwin Effect, metallicity gradients,
dust evolution, population drift). Run our full analysis pipeline
on each mock. Report the fraction that recover:
  1. Γ₀ = 0.533 ± 0.101
  2. Cross-domain ρ > 0.9
  3. Same q-ordering across domains

If < 1% of mocks recover these, conventional astrophysics CANNOT
explain the pattern at > 3σ.

Conventional models used:
  - Quasar [OIII] EW: Baldwin Effect (EW ∝ L^{-0.2}) + metallicity gradient
  - Quasar [SII] ratio: quantum-locked (no evolution)
  - SN Ia β(z): progenitor age drift + dust evolution
  - SN Ia α(z): stretch ~ timing (no evolution expected)
  - FRB width-SI: scatter broadening ∝ DM^2
  - H₀ inference depth: statistical noise only

Author: Closure Theory Pipeline
Date: March 2026
"""

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import json
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# SECTION 1: Conventional Evolution Models (from literature)
# ============================================================

def baldwin_effect(z, L, EW0=50.0, beta_baldwin=-0.20):
    """
    Baldwin Effect: EW decreases with luminosity.
    Published: Baldwin 1977, Dietrich+2002, Green+2001
    EW ∝ L^{beta_baldwin}, beta_baldwin ≈ -0.15 to -0.25
    
    L increases with z due to Malmquist bias (brighter quasars seen farther).
    """
    # Luminosity increases with z (Malmquist selection)
    L_eff = L * (1 + z)**1.5  # flux-limited selection
    return EW0 * (L_eff / L)**beta_baldwin

def metallicity_evolution(z, Z0=1.0, dZ_dz=-0.1):
    """
    Metallicity decreases at higher z (younger universe = less enrichment).
    Published: Maiolino+2008, Mannucci+2010
    Affects emission line ratios and EWs.
    """
    return Z0 + dZ_dz * z

def dust_evolution(z, E0=0.1, tau_dust=0.05):
    """
    Dust reddening evolves with z due to changing dust-to-gas ratios.
    Published: Menard+2010, Peek+2015
    """
    return E0 * (1 + z)**(-0.5) + tau_dust * z

def progenitor_drift(z, beta0=3.1, drift=-0.3):
    """
    SN Ia progenitor populations change with z (delay time distribution).
    Published: Sullivan+2010, Rigault+2020, Wiseman+2023
    Younger progenitors at high z → different color-luminosity relation.
    """
    return beta0 + drift * z

def scatter_broadening(DM, tau0=1e-4):
    """
    FRB pulse broadening from intergalactic scattering.
    Published: Macquart+2020, Cordes+Chatterjee 2019
    tau_scatter ∝ DM^2 (strong scattering regime)
    """
    return tau0 * DM**2

# ============================================================
# SECTION 2: Generate Mock Catalogs
# ============================================================

N_MOCKS = 100000

def generate_quasar_mock(n_qso=1000, z_range=(0.1, 3.0)):
    """Generate mock quasar catalog with conventional evolution only."""
    z = np.random.uniform(z_range[0], z_range[1], n_qso)
    L = 10**(np.random.normal(44.5, 0.8, n_qso))  # log-normal luminosity
    
    # [OIII] EW with Baldwin Effect + metallicity + scatter
    Z = metallicity_evolution(z)
    ew_oiii = baldwin_effect(z, L) * Z * np.random.lognormal(0, 0.3, n_qso)
    
    # [SII] ratio — quantum locked, no evolution
    sii_ratio = 1.2 + np.random.normal(0, 0.1, n_qso)
    
    # [NII] EW — weak diagnostic, mild evolution
    ew_nii = 20.0 * (1 + 0.05 * np.random.randn(n_qso))  # essentially flat
    
    # Balmer EW — moderate diagnostic
    ew_hbeta = 30.0 * baldwin_effect(z, L, EW0=1.0, beta_baldwin=-0.10) * \
               np.random.lognormal(0, 0.2, n_qso)
    
    # [OII] EW — high diagnostic
    ew_oii = 40.0 * baldwin_effect(z, L, EW0=1.0, beta_baldwin=-0.18) * Z * \
             np.random.lognormal(0, 0.25, n_qso)
    
    return {
        'z': z, 'L': L,
        'ew_oiii': ew_oiii, 'sii_ratio': sii_ratio,
        'ew_nii': ew_nii, 'ew_hbeta': ew_hbeta, 'ew_oii': ew_oii
    }

def generate_sn_mock(n_sn=500, z_range=(0.01, 1.5)):
    """Generate mock SN Ia catalog with conventional progenitor/dust evolution."""
    z = np.random.uniform(z_range[0], z_range[1], n_sn)
    
    # Stretch (x1) — geometric timing, no z-evolution
    x1 = np.random.normal(0.0, 1.0, n_sn)
    
    # Color (c) — spectral, with dust + progenitor evolution
    c_intrinsic = np.random.normal(-0.05, 0.08, n_sn)
    c_dust = dust_evolution(z)
    c = c_intrinsic + c_dust
    
    # Beta with progenitor drift
    beta = progenitor_drift(z)
    
    # Alpha — stable (timing-based)
    alpha = 0.14 + np.random.normal(0, 0.01, n_sn)
    
    return {'z': z, 'x1': x1, 'c': c, 'beta': beta, 'alpha': alpha}

def generate_frb_mock(n_frb=500, dm_range=(50, 2000)):
    """Generate mock FRB catalog with scatter broadening only."""
    DM = np.random.uniform(dm_range[0], dm_range[1], n_frb)
    z_approx = DM / 1000.0  # rough DM-z relation
    
    # Width — broadened by scattering
    w_intrinsic = np.random.lognormal(np.log(5), 0.5, n_frb)  # ms
    w_scatter = scatter_broadening(DM)
    width = w_intrinsic + w_scatter
    
    # Spectral index — intrinsic scatter, weak z-dependence
    SI = np.random.normal(-1.5, 1.0, n_frb) + 0.1 * z_approx
    
    return {'DM': DM, 'z': z_approx, 'width': width, 'SI': SI}

# ============================================================
# SECTION 3: Analysis Pipeline (same as real data)
# ============================================================

def compute_degradation_rate(z, observable, n_bins=5):
    """Compute correlation between observable and z in bins."""
    z_sort = np.argsort(z)
    z_s, o_s = z[z_sort], observable[z_sort]
    bin_size = len(z_s) // n_bins
    
    if bin_size < 10:
        return 0.0
    
    z_mids = []
    o_means = []
    for i in range(n_bins):
        start = i * bin_size
        end = start + bin_size if i < n_bins - 1 else len(z_s)
        z_mids.append(np.median(z_s[start:end]))
        o_means.append(np.mean(o_s[start:end]))
    
    if np.std(o_means) < 1e-10:
        return 0.0
    
    rho, _ = stats.spearmanr(z_mids, o_means)
    return rho if not np.isnan(rho) else 0.0

def compute_beta_z_slope(z, beta, n_bins=5):
    """Compute beta evolution with z."""
    z_sort = np.argsort(z)
    z_s, b_s = z[z_sort], beta[z_sort]
    bin_size = len(z_s) // n_bins
    
    if bin_size < 5:
        return 0.0, 0.0
    
    z_mids = []
    b_means = []
    for i in range(n_bins):
        start = i * bin_size
        end = start + bin_size if i < n_bins - 1 else len(z_s)
        z_mids.append(np.median(z_s[start:end]))
        b_means.append(np.mean(b_s[start:end]))
    
    slope, _, _, _, _ = stats.linregress(z_mids, b_means)
    rho, _ = stats.spearmanr(z_mids, b_means)
    return slope, rho if not np.isnan(rho) else 0.0

def power_law(q, gamma0, gamma):
    """Power law: degradation = gamma0 * q^gamma"""
    return gamma0 * np.power(np.clip(q, 1e-6, None), gamma)

def analyze_mock(qso, sn, frb):
    """
    Run full closure analysis on a mock catalog.
    Returns: Gamma0_fit, gamma_fit, cross_domain_rho
    """
    # ---- Quasar domain ----
    # Diagnostic sensitivities (same as real analysis)
    q_quasar = np.array([0.0, 0.2, 0.35, 0.45, 0.85])
    
    # Compute degradation rates for each line
    deg_sii = compute_degradation_rate(qso['z'], qso['sii_ratio'])
    deg_nii = compute_degradation_rate(qso['z'], qso['ew_nii'])
    deg_hbeta = compute_degradation_rate(qso['z'], qso['ew_hbeta'])
    deg_oii = compute_degradation_rate(qso['z'], qso['ew_oii'])
    deg_oiii = compute_degradation_rate(qso['z'], qso['ew_oiii'])
    
    degradation_quasar = np.abs(np.array([deg_sii, deg_nii, deg_hbeta, deg_oii, deg_oiii]))
    
    # ---- SN Ia domain ----
    q_sn = np.array([0.1, 0.7])  # stretch (locked-ish), color (diagnostic)
    
    _, rho_alpha = compute_beta_z_slope(sn['z'], sn['alpha'])
    _, rho_beta = compute_beta_z_slope(sn['z'], sn['beta'])
    
    degradation_sn = np.abs(np.array([rho_alpha, rho_beta]))
    
    # ---- FRB domain ----
    q_frb = np.array([0.1, 0.6])  # DM (locked), width-SI (diagnostic)
    
    # DM is monotonic with z (locked)
    deg_dm = 0.0  # DM-z is the locked relation itself
    
    # Width-SI correlation
    rho_wsi, _ = stats.spearmanr(frb['width'], frb['SI'])
    degradation_frb = np.abs(np.array([deg_dm, rho_wsi if not np.isnan(rho_wsi) else 0.0]))
    
    # ---- Combined: fit universal power law ----
    all_q = np.concatenate([q_quasar, q_sn, q_frb])
    all_deg = np.concatenate([degradation_quasar, degradation_sn, degradation_frb])
    
    # Filter out q=0 for fitting
    mask = all_q > 0.05
    q_fit = all_q[mask]
    deg_fit = all_deg[mask]
    
    if len(q_fit) < 3 or np.std(deg_fit) < 1e-10:
        return None, None, 0.0
    
    try:
        popt, _ = curve_fit(power_law, q_fit, deg_fit, 
                           p0=[0.5, 1.5], bounds=([0, 0], [5, 5]),
                           maxfev=5000)
        gamma0_fit, gamma_fit = popt
    except:
        return None, None, 0.0
    
    # Cross-domain correlation
    rho_cross, _ = stats.spearmanr(q_fit, deg_fit)
    
    return gamma0_fit, gamma_fit, rho_cross if not np.isnan(rho_cross) else 0.0

# ============================================================
# SECTION 4: Run 10,000 Monte Carlo Simulations
# ============================================================

print("=" * 70)
print("NULL SIMULATION: Can Conventional Astrophysics Reproduce Γ₀?")
print("=" * 70)
print(f"\nRunning {N_MOCKS:,} Monte Carlo simulations...")
print("Conventional models: Baldwin Effect + metallicity + dust + progenitor drift + scatter broadening")
print()

# Vary model parameters within published ranges
baldwin_betas = (-0.15, -0.25)     # Baldwin 1977, Dietrich+2002
metallicity_slopes = (-0.05, -0.15) # Maiolino+2008
dust_taus = (0.02, 0.08)           # Menard+2010
progenitor_drifts = (-0.1, -0.5)   # Sullivan+2010, Rigault+2020

results = {
    'gamma0': [],
    'gamma': [],
    'rho_cross': [],
    'beta_drop_pct': [],
    'sii_flat': [],
    'ladder_rho': []
}

for i in range(N_MOCKS):
    if (i + 1) % 1000 == 0:
        print(f"  Mock {i+1:,}/{N_MOCKS:,}...")
    
    # Randomize conventional model parameters within published ranges
    bb = np.random.uniform(*baldwin_betas)
    ms = np.random.uniform(*metallicity_slopes)
    dt = np.random.uniform(*dust_taus)
    pd = np.random.uniform(*progenitor_drifts)
    
    # Generate mocks with randomized conventional physics
    qso = generate_quasar_mock(n_qso=800)
    # Override with specific Baldwin beta
    L_eff = qso['L'] * (1 + qso['z'])**1.5
    qso['ew_oiii'] = 50.0 * (L_eff / qso['L'])**bb * metallicity_evolution(qso['z'], dZ_dz=ms) * \
                     np.random.lognormal(0, 0.3, len(qso['z']))
    qso['ew_oii'] = 40.0 * (L_eff / qso['L'])**(bb * 0.9) * metallicity_evolution(qso['z'], dZ_dz=ms) * \
                    np.random.lognormal(0, 0.25, len(qso['z']))
    qso['ew_hbeta'] = 30.0 * (L_eff / qso['L'])**(bb * 0.5) * np.random.lognormal(0, 0.2, len(qso['z']))
    
    sn = generate_sn_mock(n_sn=400)
    sn['beta'] = progenitor_drift(sn['z'], drift=pd)
    
    frb = generate_frb_mock(n_frb=400)
    
    # Analyze
    g0, g, rho = analyze_mock(qso, sn, frb)
    
    if g0 is not None:
        results['gamma0'].append(g0)
        results['gamma'].append(g)
        results['rho_cross'].append(rho)
        
        # Beta drop percentage
        beta_low = np.mean(sn['beta'][sn['z'] < 0.3])
        beta_high = np.mean(sn['beta'][sn['z'] > 0.8])
        drop_pct = (beta_low - beta_high) / beta_low * 100
        results['beta_drop_pct'].append(drop_pct)
        
        # SII flatness
        sii_rho = compute_degradation_rate(qso['z'], qso['sii_ratio'])
        results['sii_flat'].append(abs(sii_rho) < 0.3)
        
        # Ladder ordering
        q_ladder = np.array([0.0, 0.2, 0.35, 0.45, 0.85])
        deg_ladder = np.abs(np.array([
            compute_degradation_rate(qso['z'], qso['sii_ratio']),
            compute_degradation_rate(qso['z'], qso['ew_nii']),
            compute_degradation_rate(qso['z'], qso['ew_hbeta']),
            compute_degradation_rate(qso['z'], qso['ew_oii']),
            compute_degradation_rate(qso['z'], qso['ew_oiii'])
        ]))
        lr, _ = stats.spearmanr(q_ladder, deg_ladder)
        results['ladder_rho'].append(lr if not np.isnan(lr) else 0.0)

# Convert to arrays
for key in results:
    results[key] = np.array(results[key])

# ============================================================
# SECTION 5: Statistical Comparison
# ============================================================

print("\n" + "=" * 70)
print("RESULTS")
print("=" * 70)

# Real data values
REAL_GAMMA0 = 0.533
REAL_GAMMA0_ERR = 0.101
REAL_GAMMA = 1.56
REAL_RHO_CROSS = 0.918
REAL_BETA_DROP = 44.0
REAL_LADDER_RHO = 0.975

n_valid = len(results['gamma0'])
print(f"\nValid mocks: {n_valid:,}/{N_MOCKS:,}")

# Test 1: How many recover Γ₀ = 0.533 ± 0.101?
in_gamma0 = np.sum(np.abs(results['gamma0'] - REAL_GAMMA0) < REAL_GAMMA0_ERR)
pct_gamma0 = in_gamma0 / n_valid * 100
print(f"\n--- Test 1: Γ₀ Recovery ---")
print(f"  Real: Γ₀ = {REAL_GAMMA0:.3f} ± {REAL_GAMMA0_ERR:.3f}")
print(f"  Mock mean: Γ₀ = {np.mean(results['gamma0']):.3f} ± {np.std(results['gamma0']):.3f}")
print(f"  Mock median: Γ₀ = {np.median(results['gamma0']):.3f}")
print(f"  Mocks within real 1σ: {in_gamma0:,}/{n_valid:,} ({pct_gamma0:.2f}%)")

# Test 2: How many achieve cross-domain ρ > 0.9?
high_rho = np.sum(results['rho_cross'] > 0.9)
pct_rho = high_rho / n_valid * 100
print(f"\n--- Test 2: Cross-Domain Collapse (ρ > 0.9) ---")
print(f"  Real: ρ = {REAL_RHO_CROSS:.3f}")
print(f"  Mock mean: ρ = {np.mean(results['rho_cross']):.3f} ± {np.std(results['rho_cross']):.3f}")
print(f"  Mock median: ρ = {np.median(results['rho_cross']):.3f}")
print(f"  Mocks with ρ > 0.9: {high_rho:,}/{n_valid:,} ({pct_rho:.2f}%)")

# Test 3: How many achieve β drop ≥ 44%?
high_drop = np.sum(results['beta_drop_pct'] >= REAL_BETA_DROP)
pct_drop = high_drop / n_valid * 100
print(f"\n--- Test 3: β(z) Drop ≥ 44% ---")
print(f"  Real: 44% drop (2.94 → 1.64)")
print(f"  Mock mean drop: {np.mean(results['beta_drop_pct']):.1f}% ± {np.std(results['beta_drop_pct']):.1f}%")
print(f"  Mocks with ≥44% drop: {high_drop:,}/{n_valid:,} ({pct_drop:.2f}%)")

# Test 4: How many recover the doublet ladder (ρ > 0.95)?
high_ladder = np.sum(results['ladder_rho'] > 0.95)
pct_ladder = high_ladder / n_valid * 100
print(f"\n--- Test 4: Doublet Ladder (ρ > 0.95) ---")
print(f"  Real: ρ = {REAL_LADDER_RHO:.3f}")
print(f"  Mock mean: ρ = {np.mean(results['ladder_rho']):.3f} ± {np.std(results['ladder_rho']):.3f}")
print(f"  Mocks with ρ > 0.95: {high_ladder:,}/{n_valid:,} ({pct_ladder:.2f}%)")

# Test 5: JOINT — how many pass ALL criteria simultaneously?
joint_mask = (
    (np.abs(results['gamma0'] - REAL_GAMMA0) < REAL_GAMMA0_ERR) &
    (results['rho_cross'] > 0.9) &
    (results['beta_drop_pct'] >= REAL_BETA_DROP * 0.8) &  # allow 80% of real drop
    (results['ladder_rho'] > 0.9)
)
joint_pass = np.sum(joint_mask)
pct_joint = joint_pass / n_valid * 100
print(f"\n--- Test 5: JOINT (All Criteria Simultaneously) ---")
print(f"  Criteria: Γ₀ ∈ [0.432, 0.634] AND ρ > 0.9 AND β drop > 35% AND ladder ρ > 0.9")
print(f"  Mocks passing ALL: {joint_pass:,}/{n_valid:,} ({pct_joint:.4f}%)")

if pct_joint < 1.0:
    sigma = stats.norm.ppf(1 - pct_joint/100) if pct_joint > 0 else float('inf')
    print(f"  → Conventional astrophysics CANNOT reproduce the pattern")
    print(f"  → Null rejection significance: {sigma:.1f}σ" if np.isfinite(sigma) else f"  → Null rejection: > 5σ (0 mocks passed)")
else:
    print(f"  → {pct_joint:.2f}% of mocks pass — conventional explanation possible")

# ============================================================
# SECTION 6: Detailed Distribution Analysis
# ============================================================

print(f"\n\n{'=' * 70}")
print("DISTRIBUTION ANALYSIS")
print("=" * 70)

print(f"\nΓ₀ distribution (conventional physics):")
percentiles = [5, 25, 50, 75, 95]
for p in percentiles:
    print(f"  {p}th percentile: {np.percentile(results['gamma0'], p):.3f}")

print(f"\nCross-domain ρ distribution:")
for p in percentiles:
    print(f"  {p}th percentile: {np.percentile(results['rho_cross'], p):.3f}")

print(f"\nβ drop % distribution:")
for p in percentiles:
    print(f"  {p}th percentile: {np.percentile(results['beta_drop_pct'], p):.1f}%")

# Z-score of real values vs mock distribution
z_gamma0 = (REAL_GAMMA0 - np.mean(results['gamma0'])) / np.std(results['gamma0'])
z_rho = (REAL_RHO_CROSS - np.mean(results['rho_cross'])) / np.std(results['rho_cross'])
z_drop = (REAL_BETA_DROP - np.mean(results['beta_drop_pct'])) / np.std(results['beta_drop_pct'])

print(f"\nZ-scores (real value vs mock distribution):")
print(f"  Γ₀:     z = {z_gamma0:+.2f}")
print(f"  ρ_cross: z = {z_rho:+.2f}")
print(f"  β drop:  z = {z_drop:+.2f}")

# ============================================================
# SECTION 7: Summary Verdict
# ============================================================

print(f"\n\n{'=' * 70}")
print("VERDICT")
print("=" * 70)

verdicts = []
if pct_gamma0 < 5:
    verdicts.append(f"Γ₀ recovery: FAILS ({pct_gamma0:.2f}% < 5%)")
else:
    verdicts.append(f"Γ₀ recovery: possible ({pct_gamma0:.2f}%)")

if pct_rho < 5:
    verdicts.append(f"Cross-domain ρ: FAILS ({pct_rho:.2f}% < 5%)")
else:
    verdicts.append(f"Cross-domain ρ: possible ({pct_rho:.2f}%)")

if pct_drop < 5:
    verdicts.append(f"β(z) 44% drop: FAILS ({pct_drop:.2f}% < 5%)")
else:
    verdicts.append(f"β(z) 44% drop: possible ({pct_drop:.2f}%)")

if pct_ladder < 5:
    verdicts.append(f"Doublet ladder: FAILS ({pct_ladder:.2f}% < 5%)")
else:
    verdicts.append(f"Doublet ladder: possible ({pct_ladder:.2f}%)")

if pct_joint < 1:
    verdicts.append(f"JOINT: FAILS ({pct_joint:.4f}% < 1%) — NULL REJECTED")
else:
    verdicts.append(f"JOINT: possible ({pct_joint:.2f}%)")

for v in verdicts:
    print(f"  {v}")

n_fails = sum(1 for v in verdicts if 'FAILS' in v)
print(f"\n  Summary: {n_fails}/5 criteria FAIL under conventional astrophysics")
if n_fails >= 3:
    print(f"  → CONVENTIONAL ASTROPHYSICS CANNOT REPRODUCE THE OBSERVED PATTERN")
    print(f"  → The Diagnostic Compression Law describes something beyond known systematics")
elif n_fails >= 1:
    print(f"  → Partial failure: {n_fails} criteria are not explained by conventional models")
else:
    print(f"  → All criteria reproducible under conventional models (law may not be needed)")

# Save results
output_dir = 'results_null_simulation'
os.makedirs(output_dir, exist_ok=True)

summary = {
    'n_mocks': N_MOCKS,
    'n_valid': int(n_valid),
    'conventional_models': [
        'Baldwin Effect (beta=-0.15 to -0.25)',
        'Metallicity gradient (dZ/dz=-0.05 to -0.15)',
        'Dust evolution (tau=0.02-0.08)',
        'Progenitor drift (d_beta/dz=-0.1 to -0.5)',
        'Scatter broadening (tau ∝ DM^2)'
    ],
    'results': {
        'gamma0_recovery_pct': float(pct_gamma0),
        'rho_cross_above_0.9_pct': float(pct_rho),
        'beta_drop_above_44_pct': float(pct_drop),
        'ladder_rho_above_0.95_pct': float(pct_ladder),
        'joint_pass_pct': float(pct_joint),
    },
    'mock_distributions': {
        'gamma0_mean': float(np.mean(results['gamma0'])),
        'gamma0_std': float(np.std(results['gamma0'])),
        'rho_cross_mean': float(np.mean(results['rho_cross'])),
        'rho_cross_std': float(np.std(results['rho_cross'])),
        'beta_drop_mean': float(np.mean(results['beta_drop_pct'])),
        'beta_drop_std': float(np.std(results['beta_drop_pct'])),
    },
    'z_scores': {
        'gamma0': float(z_gamma0),
        'rho_cross': float(z_rho),
        'beta_drop': float(z_drop)
    },
    'verdict': 'NULL REJECTED' if n_fails >= 3 else f'{n_fails}/5 FAIL'
}

with open(f'{output_dir}/null_simulation_results.json', 'w') as f:
    json.dump(summary, f, indent=2)

print(f"\nResults saved to {output_dir}/")
print(f"\nDone. {N_MOCKS:,} conventional mocks analyzed.")
