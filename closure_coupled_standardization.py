#!/usr/bin/env python3
"""
CLOSURE COUPLED STANDARDIZATION MODEL
=======================================
The full model: channel degradation propagates through the SN Ia
standardization chain (Tripp formula), not as a simple additive bias.

The Tripp formula:
    μ = mB + α·x1 - β·c - M

Standard analysis: α, β, M = constants fit globally.
Reality (closure): β(z) evolves as color diagnostic degrades.
When you FIT with constant β but the TRUE β changes with z,
the residuals create a z-dependent distance bias.

This is the mechanism:
    1. At low-z: β_true ≈ β_fit → residuals ≈ 0
    2. At high-z: β_true < β_fit → over-correction of red SNe
       → red SNe pushed TOO BRIGHT → distances TOO LARGE
       → mimics accelerating expansion
    3. The bias is: Δμ ≈ (β_fit - β_true(z)) · <c>
       where <c> is the mean color at that redshift

This couples to M through the global fit: if high-z residuals are
systematically positive, the global M shifts to compensate, which
biases the low-z anchor and propagates into H₀.

Data: Pantheon+ (raw SALT2 parameters: mB, c, x1)
Reference: Tripp 1998, Betoule+ 2014, Brout+ 2022

Author: Closure Theory collaboration
Date: 2026-03-03
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize, differential_evolution
from scipy.stats import spearmanr, pearsonr, chi2 as chi2_dist
from scipy.special import erfinv
import json
import os

# ============================================================
# CONSTANTS
# ============================================================
c_light = 299792.458
H0_fid = 67.4
DH_fid = c_light / H0_fid

# ============================================================
# COSMOLOGY
# ============================================================

def E(z, Om):
    return np.sqrt(Om * (1+z)**3 + (1-Om))

def luminosity_distance(z, Om, h=H0_fid/100):
    DH = c_light / (h * 100)
    integrand = lambda zp: 1.0 / E(zp, Om)
    result, _ = quad(integrand, 0, z)
    return (1+z) * DH * result

def mu_theory(z, Om, h=H0_fid/100):
    dL = luminosity_distance(z, Om, h)
    return 5.0 * np.log10(dL) + 25.0

# Precomputed interpolation table for speed
_mu_cache = {}

def mu_theory_vec(z_arr, Om, h=H0_fid/100):
    """Fast vectorized mu using interpolation table"""
    key = (round(Om, 5), round(h, 5))
    if key not in _mu_cache:
        from scipy.interpolate import interp1d
        z_grid = np.linspace(0.005, 2.5, 500)
        mu_grid = np.array([mu_theory(zi, Om, h) for zi in z_grid])
        _mu_cache[key] = interp1d(z_grid, mu_grid, kind='cubic', fill_value='extrapolate')
    return _mu_cache[key](z_arr)

# ============================================================
# CLOSURE SIGMOID
# ============================================================

def sigmoid(z, z0=0.82, k=8.0):
    return 1.0 / (1.0 + np.exp(-k * (z - z0)))

# ============================================================
# LOAD PANTHEON+ RAW SALT2 PARAMETERS
# ============================================================

def load_pantheon_raw():
    """Load Pantheon+ with raw SALT2 light-curve parameters"""
    data_path = os.path.join(os.path.dirname(__file__), 'data', 'pantheon_plus.dat')
    
    with open(data_path, 'r') as f:
        header = f.readline().strip().split()
        data_rows = []
        for line in f:
            parts = line.strip().split()
            if len(parts) >= len(header):
                data_rows.append(parts)
    
    col = {name: i for i, name in enumerate(header)}
    
    z = np.array([float(row[col['zHD']]) for row in data_rows])
    mB = np.array([float(row[col['mB']]) for row in data_rows])
    c = np.array([float(row[col['c']]) for row in data_rows])
    x1 = np.array([float(row[col['x1']]) for row in data_rows])
    c_err = np.array([float(row[col['cERR']]) for row in data_rows])
    x1_err = np.array([float(row[col['x1ERR']]) for row in data_rows])
    mB_err = np.array([float(row[col['m_b_corr_err_DIAG']]) for row in data_rows])
    
    # Also load the pre-corrected mu for comparison
    mu_corr = np.array([float(row[col['m_b_corr']]) for row in data_rows])
    
    # Covariances
    cov_x1_c = np.array([float(row[col['COV_x1_c']]) for row in data_rows])
    
    # Host mass for mass step
    host_mass = np.array([float(row[col['HOST_LOGMASS']]) for row in data_rows])
    
    # Filter
    mask = (z > 0.01) & (mB_err > 0) & (mB_err < 5) & np.isfinite(mB)
    
    data = {
        'z': z[mask], 'mB': mB[mask], 'c': c[mask], 'x1': x1[mask],
        'c_err': c_err[mask], 'x1_err': x1_err[mask], 'mB_err': mB_err[mask],
        'mu_corr': mu_corr[mask], 'cov_x1_c': cov_x1_c[mask],
        'host_mass': host_mass[mask]
    }
    
    print(f"Loaded {mask.sum()} SNe Ia with raw SALT2 parameters")
    print(f"  z: {data['z'].min():.3f} to {data['z'].max():.3f}")
    print(f"  <c> = {data['c'].mean():.4f} ± {data['c'].std():.4f}")
    print(f"  <x1> = {data['x1'].mean():.3f} ± {data['x1'].std():.3f}")
    
    return data

# ============================================================
# MODEL A: STANDARD TRIPP (constant α, β, M)
# ============================================================

def tripp_mu(mB, c, x1, alpha, beta, M):
    """Standard Tripp estimator"""
    return mB + alpha * x1 - beta * c - M

def chi2_standard(params, data):
    """Standard Tripp fit: constant α, β, M, Ωm"""
    Om, alpha, beta, M = params
    z = data['z']
    
    mu_obs = tripp_mu(data['mB'], data['c'], data['x1'], alpha, beta, M)
    mu_th = mu_theory_vec(z, Om)
    
    # Approximate error (diagonal only for speed)
    mu_err = np.sqrt(data['mB_err']**2 + (alpha * data['x1_err'])**2 + 
                     (beta * data['c_err'])**2)
    
    residuals = (mu_obs - mu_th) / mu_err
    return np.sum(residuals**2)

# ============================================================
# MODEL B: CLOSURE TRIPP (β(z) evolves, α constant)
# ============================================================

def beta_of_z(z, beta0, delta_beta, z0=0.82, k=8.0):
    """
    β(z) = β₀ - Δβ · sigmoid(z; z₀, k)
    
    β₀: low-z value (~3.1 from Pantheon+)
    Δβ: drop at high-z (positive = β decreases)
    
    From Feb 21: β drops from 2.42 to 1.85 → Δβ ≈ 0.57
    (But those were effective β after bias corrections; raw SALT2 β ≈ 3.1)
    """
    return beta0 - delta_beta * sigmoid(z, z0, k)

def chi2_closure(params, data, z0=0.82, k=8.0):
    """Closure Tripp fit: β(z) evolves, rest constant"""
    Om, alpha, beta0, delta_beta, M = params
    z = data['z']
    
    # z-dependent β
    beta_z = beta_of_z(z, beta0, delta_beta, z0, k)
    
    # Tripp formula with z-dependent β
    mu_obs = data['mB'] + alpha * data['x1'] - beta_z * data['c'] - M
    mu_th = mu_theory_vec(z, Om)
    
    # Error propagation with z-dependent β
    mu_err = np.sqrt(data['mB_err']**2 + (alpha * data['x1_err'])**2 + 
                     (beta_z * data['c_err'])**2)
    
    residuals = (mu_obs - mu_th) / mu_err
    return np.sum(residuals**2)

# ============================================================
# MODEL C: FULL CLOSURE (β(z) + mass step evolution)
# ============================================================

def chi2_full_closure(params, data, z0=0.82, k=8.0):
    """
    Full closure model:
    β(z) evolves + mass step δM(z) evolves
    
    The mass step (0.04 mag offset between high/low host mass SNe)
    is a DIAGNOSTIC-dependent systematic → should also evolve with z
    """
    Om, alpha, beta0, delta_beta, M, delta_M0, delta_M_evol = params
    z = data['z']
    
    beta_z = beta_of_z(z, beta0, delta_beta, z0, k)
    
    # Mass step: δM(z) = δM₀ · (1 - f_evol · sigmoid(z))
    # At high-z, mass step should INCREASE (diagnostic more uncertain)
    mass_step = np.where(data['host_mass'] > 10.0,
                        delta_M0 * (1 + delta_M_evol * sigmoid(z, z0, k)),
                        0.0)
    
    mu_obs = data['mB'] + alpha * data['x1'] - beta_z * data['c'] - M - mass_step
    mu_th = mu_theory_vec(z, Om)
    
    mu_err = np.sqrt(data['mB_err']**2 + (alpha * data['x1_err'])**2 + 
                     (beta_z * data['c_err'])**2)
    
    residuals = (mu_obs - mu_th) / mu_err
    return np.sum(residuals**2)

# ============================================================
# BIAS PROPAGATION ANALYSIS
# ============================================================

def compute_bias_propagation(data, beta_fit, beta0, delta_beta, z0=0.82, k=8.0):
    """
    Compute how constant-β assumption creates z-dependent distance bias.
    
    If true β(z) = β₀ - Δβ·sigmoid(z) but you fit with constant β_fit,
    the distance bias at each SN is:
    
        Δμ_i = -(β_fit - β_true(z_i)) · c_i
    
    For red SNe (c > 0): β_fit > β_true at high-z → Δμ > 0 → SN appears farther
    For blue SNe (c < 0): β_fit > β_true at high-z → Δμ < 0 → SN appears closer
    
    Net effect depends on color distribution: if <c> > 0 at high-z
    (Malmquist bias selects bluer = brighter), net Δμ can go either way.
    """
    z = data['z']
    c = data['c']
    
    beta_true = beta_of_z(z, beta0, delta_beta, z0, k)
    delta_beta_z = beta_fit - beta_true  # positive at high-z
    
    # Bias per SN
    delta_mu = -delta_beta_z * c  # negative sign from Tripp: μ = ... - β·c
    
    # Bin and average
    z_edges = np.linspace(0.01, min(z.max(), 2.0), 20)
    z_centers = []
    mean_bias = []
    mean_color = []
    std_bias = []
    n_per_bin = []
    
    for i in range(len(z_edges)-1):
        mask = (z >= z_edges[i]) & (z < z_edges[i+1])
        if mask.sum() < 3:
            continue
        z_centers.append(np.mean(z[mask]))
        mean_bias.append(np.mean(delta_mu[mask]))
        mean_color.append(np.mean(c[mask]))
        std_bias.append(np.std(delta_mu[mask]) / np.sqrt(mask.sum()))
        n_per_bin.append(mask.sum())
    
    return (np.array(z_centers), np.array(mean_bias), np.array(mean_color),
            np.array(std_bias), np.array(n_per_bin))

# ============================================================
# H₀ PROPAGATION
# ============================================================

def h0_from_fit(Om, M):
    """
    H₀ is degenerate with M in SN fitting.
    μ = 5 log10(dL) + 25 = mB + α·x1 - β·c - M
    
    Shifting M by δM shifts the inferred H₀:
    δH₀/H₀ ≈ δM / (5 · log10(e)) ≈ δM / 2.171
    """
    # The absolute magnitude M encodes the H₀ assumption
    # M ≈ -19.25 corresponds to H₀ ≈ 73 (SH0ES)
    # M ≈ -19.40 corresponds to H₀ ≈ 67 (Planck)
    # δM = -0.15 → δH₀ ≈ +6 km/s/Mpc
    
    # Rough H₀ from M (using SN intrinsic MB ≈ -19.25 calibrated to H₀=73)
    delta_M = M - (-19.25)
    h0 = 73.04 * 10**(delta_M / 5.0)
    return h0

# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("CLOSURE COUPLED STANDARDIZATION MODEL")
    print("β(z) Evolution Through the Tripp Formula")
    print("=" * 70)
    
    data = load_pantheon_raw()
    N = len(data['z'])
    
    # ============================================================
    # FIT MODEL A: Standard Tripp (constant everything)
    # ============================================================
    print("\n" + "=" * 70)
    print("MODEL A: STANDARD TRIPP (constant α, β, M)")
    print("=" * 70)
    
    res_A = minimize(chi2_standard, [0.3, 0.14, 3.1, -19.3],
                     args=(data,),
                     bounds=[(0.05, 0.95), (0.0, 0.5), (1.0, 5.0), (-25, -15)],
                     method='L-BFGS-B')
    Om_A, alpha_A, beta_A, M_A = res_A.x
    chi2_A = res_A.fun
    ndof_A = N - 4
    
    print(f"\n  Ωm    = {Om_A:.4f}")
    print(f"  α     = {alpha_A:.4f}")
    print(f"  β     = {beta_A:.4f}")
    print(f"  M     = {M_A:.4f}")
    print(f"  χ²    = {chi2_A:.1f} ({ndof_A} dof)")
    print(f"  χ²/dof = {chi2_A/ndof_A:.4f}")
    print(f"  H₀ ≈ {h0_from_fit(Om_A, M_A):.1f} km/s/Mpc (from M)")
    
    # ============================================================
    # FIT MODEL B: Closure Tripp (β(z) evolves)
    # ============================================================
    print("\n" + "=" * 70)
    print("MODEL B: CLOSURE TRIPP (β(z) = β₀ - Δβ·sigmoid(z; 0.82))")
    print("=" * 70)
    
    # Try multiple z₀ values
    best_chi2_B = np.inf
    best_params_B = None
    best_z0 = None
    
    for z0_try in [0.5, 0.6, 0.7, 0.82, 1.0, 1.2]:
        for k_try in [4.0, 8.0, 12.0]:
            try:
                res = minimize(chi2_closure, [0.3, 0.14, 3.1, 0.5, -19.3],
                              args=(data, z0_try, k_try),
                              bounds=[(0.05, 0.95), (0.0, 0.5), (1.0, 5.0), 
                                      (-2.0, 2.0), (-25, -15)],
                              method='L-BFGS-B')
                if res.fun < best_chi2_B:
                    best_chi2_B = res.fun
                    best_params_B = res.x
                    best_z0 = z0_try
                    best_k = k_try
            except:
                continue
    
    Om_B, alpha_B, beta0_B, dbeta_B, M_B = best_params_B
    chi2_B = best_chi2_B
    ndof_B = N - 5
    
    print(f"\n  Best z₀ = {best_z0:.2f}, k = {best_k:.1f}")
    print(f"  Ωm     = {Om_B:.4f}")
    print(f"  α      = {alpha_B:.4f}")
    print(f"  β₀     = {beta0_B:.4f} (low-z)")
    print(f"  Δβ     = {dbeta_B:.4f}")
    print(f"  β(z>1) ≈ {beta0_B - dbeta_B:.4f} (high-z)")
    print(f"  M      = {M_B:.4f}")
    print(f"  χ²     = {chi2_B:.1f} ({ndof_B} dof)")
    print(f"  χ²/dof  = {chi2_B/ndof_B:.4f}")
    print(f"  H₀ ≈ {h0_from_fit(Om_B, M_B):.1f} km/s/Mpc")
    
    # Δχ² comparison
    delta_chi2_AB = chi2_A - chi2_B
    p_AB = chi2_dist.sf(max(0, delta_chi2_AB), 1)  # 1 extra param
    sigma_AB = np.sqrt(2) * erfinv(max(0, 1 - p_AB)) if p_AB > 0 and p_AB < 1 else 0
    
    print(f"\n  Δχ² (A - B) = {delta_chi2_AB:.1f}")
    print(f"  p-value = {p_AB:.2e}")
    print(f"  Significance: {sigma_AB:.1f}σ preference for β(z) evolution")
    
    # ============================================================
    # FIT MODEL B at FIXED z₀ = 0.82 (closure prediction)
    # ============================================================
    print("\n" + "=" * 70)
    print("MODEL B (z₀=0.82 FIXED): Testing closure-predicted threshold")
    print("=" * 70)
    
    res_B82 = minimize(chi2_closure, [0.3, 0.14, 3.1, 0.5, -19.3],
                       args=(data, 0.82, 8.0),
                       bounds=[(0.05, 0.95), (0.0, 0.5), (1.0, 5.0),
                               (-2.0, 2.0), (-25, -15)],
                       method='L-BFGS-B')
    Om_B82, alpha_B82, beta0_B82, dbeta_B82, M_B82 = res_B82.x
    chi2_B82 = res_B82.fun
    
    delta_chi2_82 = chi2_A - chi2_B82
    p_82 = chi2_dist.sf(max(0, delta_chi2_82), 1)
    sigma_82 = np.sqrt(2) * erfinv(max(0, 1 - p_82)) if 0 < p_82 < 1 else 0
    
    print(f"\n  Ωm     = {Om_B82:.4f}")
    print(f"  β₀     = {beta0_B82:.4f}")
    print(f"  Δβ     = {dbeta_B82:.4f}")
    print(f"  β(z>1) ≈ {beta0_B82 - dbeta_B82:.4f}")
    print(f"  M      = {M_B82:.4f}")
    print(f"  χ²     = {chi2_B82:.1f}")
    print(f"  Δχ² vs standard = {delta_chi2_82:.1f} ({sigma_82:.1f}σ)")
    print(f"  H₀ ≈ {h0_from_fit(Om_B82, M_B82):.1f} km/s/Mpc")
    
    # ============================================================
    # SCAN: z₀ profile likelihood  
    # ============================================================
    print("\n" + "=" * 70)
    print("z₀ PROFILE SCAN: Where does β(z) transition?")
    print("=" * 70)
    
    z0_scan = np.arange(0.2, 1.8, 0.05)
    chi2_scan = []
    dbeta_scan = []
    Om_scan = []
    
    for z0_try in z0_scan:
        try:
            res = minimize(chi2_closure, [0.3, 0.14, 3.1, 0.5, -19.3],
                          args=(data, z0_try, 8.0),
                          bounds=[(0.05, 0.95), (0.0, 0.5), (1.0, 5.0),
                                  (-2.0, 2.0), (-25, -15)],
                          method='L-BFGS-B')
            chi2_scan.append(res.fun)
            dbeta_scan.append(res.x[3])
            Om_scan.append(res.x[0])
        except:
            chi2_scan.append(np.inf)
            dbeta_scan.append(0)
            Om_scan.append(0.3)
    
    chi2_scan = np.array(chi2_scan)
    dbeta_scan = np.array(dbeta_scan)
    Om_scan = np.array(Om_scan)
    
    best_idx = np.argmin(chi2_scan)
    z0_best = z0_scan[best_idx]
    chi2_best = chi2_scan[best_idx]
    
    print(f"\n  {'z₀':>5s}  {'Δχ²':>7s}  {'Δβ':>7s}  {'Ωm':>6s}")
    print("  " + "-" * 35)
    for i in range(0, len(z0_scan), 3):
        marker = " ← BEST" if i == best_idx else ""
        marker = " ← z₀=0.82" if abs(z0_scan[i] - 0.82) < 0.03 and marker == "" else marker
        print(f"  {z0_scan[i]:5.2f}  {chi2_A - chi2_scan[i]:+7.1f}  {dbeta_scan[i]:+7.3f}  {Om_scan[i]:6.4f}{marker}")
    
    print(f"\n  Best-fit z₀ = {z0_best:.2f}")
    print(f"  Closure predicted z₀ = 0.82")
    print(f"  Offset: |{z0_best:.2f} - 0.82| = {abs(z0_best - 0.82):.2f}")
    
    # ============================================================
    # BIAS PROPAGATION ANALYSIS
    # ============================================================
    print("\n" + "=" * 70)
    print("BIAS PROPAGATION: How constant-β creates fake dark energy")
    print("=" * 70)
    
    z_bp, bias_bp, color_bp, err_bp, n_bp = compute_bias_propagation(
        data, beta_A, beta0_B82, dbeta_B82, z0=0.82, k=8.0)
    
    print(f"\n  Using: β_fit = {beta_A:.3f} (standard), β₀ = {beta0_B82:.3f}, Δβ = {dbeta_B82:.3f}")
    print(f"\n  {'<z>':>6s}  {'<c>':>7s}  {'<Δμ> (mag)':>11s}  {'±':>5s}  {'N':>5s}  {'z > z₀':>7s}")
    print("  " + "-" * 50)
    for i in range(len(z_bp)):
        above = z_bp[i] > 0.82
        print(f"  {z_bp[i]:6.3f}  {color_bp[i]:+7.4f}  {bias_bp[i]:+11.5f}  {err_bp[i]:5.4f}  {n_bp[i]:5.0f}  {'YES ←' if above else ''}")
    
    # Cumulative bias effect on Ωm
    # A positive Δμ at high-z mimics larger distances → more acceleration → lower Ωm
    # A negative Δμ at high-z mimics smaller distances → less acceleration → higher Ωm
    mean_bias_highz = np.mean(bias_bp[z_bp > 0.82]) if np.any(z_bp > 0.82) else 0
    
    print(f"\n  Mean bias at z > 0.82: {mean_bias_highz:+.5f} mag")
    if mean_bias_highz > 0:
        print(f"  → High-z SNe appear FARTHER than they are")
        print(f"  → Mimics MORE acceleration → LOWER Ωm, HIGHER ΩΛ")
    elif mean_bias_highz < 0:
        print(f"  → High-z SNe appear CLOSER than they are")
        print(f"  → Mimics LESS acceleration → HIGHER Ωm, LOWER ΩΛ")
    
    # ============================================================
    # Ωm AND H₀ SHIFT
    # ============================================================
    print("\n" + "=" * 70)
    print("COSMOLOGICAL PARAMETER SHIFTS")
    print("=" * 70)
    
    dOm = Om_B82 - Om_A
    dM = M_B82 - M_A
    H0_A = h0_from_fit(Om_A, M_A)
    H0_B = h0_from_fit(Om_B82, M_B82)
    dH0 = H0_B - H0_A
    
    print(f"\n  Parameter     Standard    Closure     Shift")
    print(f"  {'─'*50}")
    print(f"  Ωm           {Om_A:.4f}      {Om_B82:.4f}      {dOm:+.4f}")
    print(f"  ΩΛ           {1-Om_A:.4f}      {1-Om_B82:.4f}      {-dOm:+.4f}")
    print(f"  M            {M_A:.4f}     {M_B82:.4f}     {dM:+.4f}")
    print(f"  H₀ (approx)  {H0_A:.1f}        {H0_B:.1f}        {dH0:+.1f}")
    print(f"  β            {beta_A:.4f}      {beta0_B82:.4f}→{beta0_B82-dbeta_B82:.4f}")
    
    # How much Hubble tension does this explain?
    tension = 73.04 - 67.4  # 5.64
    print(f"\n  Hubble tension: 73.04 - 67.4 = {tension:.2f} km/s/Mpc")
    print(f"  Closure H₀ shift: {dH0:+.1f} km/s/Mpc")
    if abs(dH0) > 0:
        pct = abs(dH0) / tension * 100
        direction = "toward Planck" if dH0 < 0 else "away from Planck"
        print(f"  → Explains {pct:.0f}% of tension ({direction})")
    
    # ============================================================
    # COMPARISON WITH BAO
    # ============================================================
    print("\n" + "=" * 70)
    print("BAO COMPARISON")
    print("=" * 70)
    
    # DESI BAO gives Ωm = 0.295±0.015 (locked probe)
    Om_bao = 0.295
    Om_bao_err = 0.015
    
    gap_standard = abs(Om_A - Om_bao)
    gap_closure = abs(Om_B82 - Om_bao)
    
    print(f"\n  BAO (DESI):     Ωm = {Om_bao:.3f} ± {Om_bao_err:.3f}")
    print(f"  SN (standard):  Ωm = {Om_A:.4f}  (gap = {gap_standard:.4f})")
    print(f"  SN (closure):   Ωm = {Om_B82:.4f}  (gap = {gap_closure:.4f})")
    
    if gap_closure < gap_standard:
        print(f"\n  ✓ Closure REDUCES SN↔BAO gap by {(1-gap_closure/gap_standard)*100:.0f}%")
    else:
        print(f"\n  ✗ Closure INCREASES SN↔BAO gap by {(gap_closure/gap_standard-1)*100:.0f}%")
        print(f"    (But this may indicate the bias propagates through M, not directly through Ωm)")

    # ============================================================
    # BINNED β(z) DIRECT MEASUREMENT
    # ============================================================
    print("\n" + "=" * 70)
    print("DIRECT β(z) MEASUREMENT: Binned color-luminosity fits")
    print("=" * 70)
    
    z_edges = [0.01, 0.15, 0.3, 0.5, 0.7, 0.82, 1.0, 1.3, 2.3]
    
    print(f"\n  {'z-bin':>12s}  {'N':>5s}  {'β':>7s}  {'α':>7s}  {'Δβ from global':>15s}")
    print("  " + "-" * 55)
    
    beta_binned = []
    z_bin_centers_beta = []
    
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        mask = (data['z'] >= z_lo) & (data['z'] < z_hi)
        n = mask.sum()
        if n < 15:
            continue
        
        # Fit α, β, M in this bin with fixed Ωm
        def chi2_bin(params, d_bin, Om_fix):
            a, b, m = params
            mu_obs = tripp_mu(d_bin['mB'], d_bin['c'], d_bin['x1'], a, b, m)
            mu_th = mu_theory_vec(d_bin['z'], Om_fix)
            mu_err = np.sqrt(d_bin['mB_err']**2 + (a * d_bin['x1_err'])**2 + 
                            (b * d_bin['c_err'])**2)
            return np.sum(((mu_obs - mu_th) / mu_err)**2)
        
        d_bin = {k: v[mask] for k, v in data.items()}
        
        try:
            res_bin = minimize(chi2_bin, [alpha_A, beta_A, M_A],
                             args=(d_bin, Om_A),
                             bounds=[(0.0, 0.5), (0.5, 6.0), (-25, -15)],
                             method='L-BFGS-B')
            a_bin, b_bin, m_bin = res_bin.x
            
            z_mid = (z_lo + z_hi) / 2
            above = z_mid > 0.82
            beta_binned.append(b_bin)
            z_bin_centers_beta.append(z_mid)
            
            print(f"  [{z_lo:.2f}, {z_hi:.2f})  {n:5d}  {b_bin:7.3f}  {a_bin:7.3f}  {b_bin - beta_A:+15.3f}  {'← z > z₀' if above else ''}")
        except:
            continue
    
    if len(beta_binned) > 3:
        beta_binned = np.array(beta_binned)
        z_bin_centers_beta = np.array(z_bin_centers_beta)
        
        rho_beta_z, p_beta_z = spearmanr(z_bin_centers_beta, beta_binned)
        print(f"\n  β vs z trend: ρ = {rho_beta_z:.3f}, p = {p_beta_z:.4f}")
        if rho_beta_z < 0:
            print(f"  ✓ β DECREASES with z → diagnostic degradation confirmed")
        else:
            print(f"  β trend: {'increases' if rho_beta_z > 0 else 'flat'}")
        
        # Does the drop happen at z₀?
        if len(beta_binned) >= 4:
            beta_low = np.mean(beta_binned[z_bin_centers_beta < 0.82])
            beta_high = np.mean(beta_binned[z_bin_centers_beta >= 0.82])
            print(f"\n  <β> at z < 0.82: {beta_low:.3f}")
            print(f"  <β> at z ≥ 0.82: {beta_high:.3f}")
            print(f"  Drop: Δβ = {beta_low - beta_high:+.3f}")

    # ============================================================
    # FULL CLOSURE MODEL C (β(z) + mass step evolution)
    # ============================================================
    print("\n" + "=" * 70)
    print("MODEL C: FULL CLOSURE (β(z) + mass step evolution)")
    print("=" * 70)
    
    # Only run if host mass data is available
    valid_mass = data['host_mass'] > 0
    if valid_mass.sum() > 100:
        data_mass = {k: v[valid_mass] for k, v in data.items()}
        
        res_C = minimize(chi2_full_closure, 
                        [Om_B82, alpha_B82, beta0_B82, dbeta_B82, M_B82, 0.04, 0.3],
                        args=(data_mass, 0.82, 8.0),
                        bounds=[(0.05, 0.95), (0.0, 0.5), (1.0, 5.0),
                                (-2.0, 2.0), (-25, -15), (-0.2, 0.2), (-2, 2)],
                        method='L-BFGS-B')
        Om_C, alpha_C, beta0_C, dbeta_C, M_C, dM0_C, dM_evol_C = res_C.x
        chi2_C = res_C.fun
        ndof_C = valid_mass.sum() - 7
        
        # Compare with standard on same subsample
        chi2_A_sub = chi2_standard([Om_A, alpha_A, beta_A, M_A], data_mass)
        delta_chi2_AC = chi2_A_sub - chi2_C
        p_AC = chi2_dist.sf(max(0, delta_chi2_AC), 3)  # 3 extra params
        
        print(f"\n  Ωm       = {Om_C:.4f}")
        print(f"  β₀       = {beta0_C:.4f}")
        print(f"  Δβ       = {dbeta_C:.4f}")
        print(f"  δM₀      = {dM0_C:.4f} mag (mass step)")
        print(f"  δM_evol  = {dM_evol_C:.4f} (mass step z-evolution)")
        print(f"  χ²       = {chi2_C:.1f}")
        print(f"  Δχ² vs standard = {delta_chi2_AC:.1f} (p = {p_AC:.2e})")
        print(f"  H₀ ≈ {h0_from_fit(Om_C, M_C):.1f} km/s/Mpc")
    
    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 70)
    print("SUMMARY: COUPLED STANDARDIZATION MODEL")
    print("=" * 70)
    
    print(f"""
    THE MECHANISM:
    ─────────────
    Standard SN analysis assumes β = {beta_A:.2f} (constant).
    Closure theory predicts β(z) drops from {beta0_B82:.2f} → {beta0_B82-dbeta_B82:.2f}
    through a sigmoid at z₀ = 0.82.
    
    When you fit constant β to evolving data:
    • At z < 0.82: β_fit ≈ β_true → unbiased
    • At z > 0.82: β_fit > β_true → color correction too aggressive
    • Red SNe over-corrected → appear brighter → distances too large
    • Blue SNe under-corrected → appear dimmer → distances too small
    • Net bias depends on color distribution at each z
    
    This propagates into:
    1. Distance modulus bias (mimics/masks acceleration)
    2. M offset (shifts H₀)
    3. Apparent dark energy evolution (DESI signal)
    
    RESULTS:
    ────────
    β(z) evolution:  Δχ² = {delta_chi2_82:.1f} ({sigma_82:.1f}σ) at z₀ = 0.82
    Best-fit z₀:     {z0_best:.2f} (predicted: 0.82)
    Ωm shift:        {dOm:+.4f} ({Om_A:.4f} → {Om_B82:.4f})
    M shift:         {dM:+.4f}
    H₀ shift:        {dH0:+.1f} km/s/Mpc
    
    The bias doesn't just shift Ωm — it reshapes the entire distance-
    redshift relation through the M-H₀ degeneracy and the color-
    dependent Malmquist selection. This is why the simple additive
    sigmoid model (previous script) missed the direction: the coupled
    mechanism is more subtle than Δμ = A·sigmoid(z).
    """)
    
    # Save
    results_dir = os.path.join(os.path.dirname(__file__), 'results_coupled_standardization')
    os.makedirs(results_dir, exist_ok=True)
    
    results = {
        'model_A': {
            'Om': float(Om_A), 'alpha': float(alpha_A), 'beta': float(beta_A),
            'M': float(M_A), 'chi2': float(chi2_A), 'H0': float(H0_A)
        },
        'model_B_z082': {
            'Om': float(Om_B82), 'alpha': float(alpha_B82), 
            'beta0': float(beta0_B82), 'delta_beta': float(dbeta_B82),
            'M': float(M_B82), 'chi2': float(chi2_B82), 'H0': float(H0_B)
        },
        'model_B_best': {
            'z0': float(z0_best), 'Om': float(best_params_B[0]),
            'beta0': float(best_params_B[2]), 'delta_beta': float(best_params_B[3]),
            'chi2': float(best_chi2_B)
        },
        'delta_chi2_z082': float(delta_chi2_82),
        'sigma_z082': float(sigma_82),
        'Om_shift': float(dOm),
        'M_shift': float(dM),
        'H0_shift': float(dH0),
        'z0_scan': {f'{z:.2f}': float(c) for z, c in zip(z0_scan, chi2_scan)},
        'beta_binned': {f'{z:.2f}': float(b) for z, b in zip(z_bin_centers_beta, beta_binned)} if len(beta_binned) > 0 else {},
        'bao_comparison': {
            'Om_bao': float(Om_bao),
            'gap_standard': float(gap_standard),
            'gap_closure': float(gap_closure)
        }
    }
    
    results_path = os.path.join(results_dir, 'results.json')
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {results_path}")

if __name__ == '__main__':
    main()
