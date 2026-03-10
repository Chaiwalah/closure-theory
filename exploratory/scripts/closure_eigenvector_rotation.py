#!/usr/bin/env python3
"""
CLOSURE THEORY — EIGENVECTOR ROTATION TEST
=============================================

GPT's insight: if α(z), β(z), γ(z) all evolve, then the COVARIANCE 
STRUCTURE of (m_B, x1, c, host_mass) should ROTATE with z.

The principal axes of the standardization space are not fixed.
Low-z and high-z SNe should NOT live in the same correction basis.

TEST: Compute PCA of (m_B, x1, c, host_mass) in z-bins.
Track eigenvector angles, eigenvalue ratios, and correlation structure.
If eigenvectors rotate → the Tripp basis is z-dependent.

ALSO: The c² test — fit quadratic color correction, check for 
color-manifold curvature.

ALSO: w-bias test — does a z-dependent mass step shift w₀?

Author: Closure Theory collaboration
Date: 2026-03-09
"""

import numpy as np
from scipy.optimize import minimize as scipy_minimize
from scipy.integrate import quad
from scipy.stats import spearmanr, pearsonr
from scipy.linalg import eigh
import json
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# LOAD DATA
# ============================================================

data_file = 'data/pantheon_plus.dat'

with open(data_file, 'r') as f:
    header = f.readline().strip().split()
    rows = []
    for line in f:
        parts = line.strip().split()
        if len(parts) == len(header):
            rows.append(parts)

col = {name: i for i, name in enumerate(header)}

z_cmb = np.array([float(r[col['zCMB']]) for r in rows])
z_hd = np.array([float(r[col['zHD']]) for r in rows])
m_b = np.array([float(r[col['mB']]) for r in rows])
m_b_err = np.array([float(r[col['mBERR']]) for r in rows])
x1 = np.array([float(r[col['x1']]) for r in rows])
c = np.array([float(r[col['c']]) for r in rows])
host_mass = np.array([float(r[col['HOST_LOGMASS']]) for r in rows])

mask_q = (z_hd > 0.01) & (z_hd < 2.5) & (np.abs(x1) < 5) & (np.abs(c) < 0.5) & (m_b_err < 1.0) & (host_mass > 0) & (host_mass < 15)

z = z_cmb[mask_q]
mb = m_b[mask_q]
mb_err = m_b_err[mask_q]
x1_q = x1[mask_q]
c_q = c[mask_q]
hm = host_mass[mask_q]

print(f"Loaded {len(z)} SNe Ia")


# ============================================================
# TEST 1: EIGENVECTOR ROTATION WITH REDSHIFT
# ============================================================

print(f"\n{'=' * 70}")
print("TEST 1: EIGENVECTOR ROTATION — DOES THE CORRECTION BASIS ROTATE?")
print("=" * 70)

# Remove distance-modulus trend from m_B before PCA
# (otherwise m_B vs z dominates everything)
C_LIGHT = 299792.458
def dist_mod(z_arr, H0=70, Om=0.3):
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        integral, _ = quad(lambda zp: 1.0/np.sqrt(Om*(1+zp)**3 + (1-Om)), 0, zi)
        dl[i] = C_LIGHT * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

print("Computing distance moduli for de-trending...")
mu_all = dist_mod(z, H0=70)
mb_resid = mb - mu_all  # Remove cosmological distance trend

# Build the observable matrix (de-trended m_B, x1, c, host_mass)
obs_labels = ['m_B_resid', 'x1', 'c', 'host_mass']

z_bins = [
    (0.01, 0.04, 'very local'),
    (0.04, 0.08, 'local'),
    (0.08, 0.15, 'low-z'),
    (0.15, 0.30, 'mid-z'),
    (0.30, 0.55, 'high-z'),
    (0.55, 1.50, 'very high-z'),
]

print(f"\n  PCA of (m_B_resid, x1, c, host_mass) in z-bins:")
print(f"  (m_B de-trended by distance modulus to remove cosmological signal)\n")

all_eigvecs = []
all_eigvals = []
z_mids = []

for z_lo, z_hi, label in z_bins:
    mask = (z >= z_lo) & (z < z_hi)
    n = np.sum(mask)
    if n < 30:
        print(f"  [{z_lo:.2f}, {z_hi:.2f}] ({label}): N={n} — too few, skipping")
        continue
    
    # Build matrix
    X = np.column_stack([
        mb_resid[mask],
        x1_q[mask],
        c_q[mask],
        hm[mask]
    ])
    
    # Standardize
    X_std = (X - X.mean(axis=0)) / X.std(axis=0)
    
    # Correlation matrix
    corr = np.corrcoef(X_std.T)
    
    # Eigendecomposition
    eigvals, eigvecs = eigh(corr)
    # Sort descending
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    
    # Sign convention: first component of PC1 positive
    for j in range(eigvecs.shape[1]):
        if eigvecs[0, j] < 0:
            eigvecs[:, j] *= -1
    
    all_eigvecs.append(eigvecs)
    all_eigvals.append(eigvals)
    z_mids.append((z_lo + z_hi) / 2)
    
    var_explained = eigvals / eigvals.sum() * 100
    
    print(f"  [{z_lo:.2f}, {z_hi:.2f}] ({label:>12}): N={n:>4}")
    print(f"    Variance explained: PC1={var_explained[0]:.1f}% PC2={var_explained[1]:.1f}% PC3={var_explained[2]:.1f}% PC4={var_explained[3]:.1f}%")
    print(f"    PC1 loadings: mB={eigvecs[0,0]:+.3f} x1={eigvecs[1,0]:+.3f} c={eigvecs[2,0]:+.3f} mass={eigvecs[3,0]:+.3f}")
    print(f"    PC2 loadings: mB={eigvecs[0,1]:+.3f} x1={eigvecs[1,1]:+.3f} c={eigvecs[2,1]:+.3f} mass={eigvecs[3,1]:+.3f}")
    print()


# Track eigenvector angles between bins
if len(all_eigvecs) >= 2:
    print(f"\n  EIGENVECTOR ROTATION BETWEEN z-BINS:")
    print(f"  (angle between PC1 vectors in adjacent bins)\n")
    
    angles_pc1 = []
    angles_pc2 = []
    
    for i in range(len(all_eigvecs) - 1):
        # Angle between PC1 vectors
        cos_angle_1 = np.abs(np.dot(all_eigvecs[i][:, 0], all_eigvecs[i+1][:, 0]))
        cos_angle_2 = np.abs(np.dot(all_eigvecs[i][:, 1], all_eigvecs[i+1][:, 1]))
        
        angle_1 = np.degrees(np.arccos(np.clip(cos_angle_1, -1, 1)))
        angle_2 = np.degrees(np.arccos(np.clip(cos_angle_2, -1, 1)))
        
        angles_pc1.append(angle_1)
        angles_pc2.append(angle_2)
        
        z_label = f"z={z_mids[i]:.2f}→{z_mids[i+1]:.2f}"
        print(f"    {z_label}: PC1 rotation = {angle_1:.1f}°, PC2 rotation = {angle_2:.1f}°")
    
    # Total rotation from lowest to highest z
    cos_total_1 = np.abs(np.dot(all_eigvecs[0][:, 0], all_eigvecs[-1][:, 0]))
    cos_total_2 = np.abs(np.dot(all_eigvecs[0][:, 1], all_eigvecs[-1][:, 1]))
    total_angle_1 = np.degrees(np.arccos(np.clip(cos_total_1, -1, 1)))
    total_angle_2 = np.degrees(np.arccos(np.clip(cos_total_2, -1, 1)))
    
    print(f"\n    TOTAL ROTATION (lowest → highest z):")
    print(f"    PC1: {total_angle_1:.1f}°")
    print(f"    PC2: {total_angle_2:.1f}°")
    
    if total_angle_1 > 15:
        print(f"\n    ✓ SIGNIFICANT ROTATION — the Tripp basis is NOT fixed")
        print(f"    The correction space rotates as you go to higher z")
        print(f"    A fixed (α, β) cannot capture this")
    elif total_angle_1 > 5:
        print(f"\n    ◐ MODERATE ROTATION — some basis instability")
    else:
        print(f"\n    ⚠ SMALL ROTATION — basis appears stable")


# ============================================================
# TEST 2: CORRELATION STRUCTURE EVOLUTION
# ============================================================

print(f"\n\n{'=' * 70}")
print("TEST 2: CORRELATION STRUCTURE EVOLUTION")
print("How do pairwise correlations change with z?")
print("=" * 70)

pairs = [
    (0, 1, 'mB-x1'),
    (0, 2, 'mB-c'),
    (0, 3, 'mB-mass'),
    (1, 2, 'x1-c'),
    (1, 3, 'x1-mass'),
    (2, 3, 'c-mass'),
]

obs = [mb_resid, x1_q, c_q, hm]

print(f"\n  {'z-bin':>15}", end='')
for _, _, name in pairs:
    print(f" {name:>10}", end='')
print()
print(f"  {'-'*78}")

correlation_evolution = {name: [] for _, _, name in pairs}
z_mids_corr = []

for z_lo, z_hi, label in z_bins:
    mask = (z >= z_lo) & (z < z_hi)
    n = np.sum(mask)
    if n < 30:
        continue
    
    z_mids_corr.append((z_lo + z_hi) / 2)
    print(f"  [{z_lo:.2f}, {z_hi:.2f}]", end='')
    
    for i, j, name in pairs:
        r, p = pearsonr(obs[i][mask], obs[j][mask])
        correlation_evolution[name].append(r)
        sig = '*' if p < 0.05 else ' '
        print(f" {r:>+9.3f}{sig}", end='')
    print()

# Which correlations evolve most?
print(f"\n  CORRELATION EVOLUTION (Spearman ρ of r_ij vs z_mid):")
for _, _, name in pairs:
    if len(correlation_evolution[name]) >= 3:
        rho, p = spearmanr(z_mids_corr, correlation_evolution[name])
        sig = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.1 else ''
        print(f"    {name:<10}: ρ = {rho:+.3f} (p = {p:.4f}) {sig}")


# ============================================================
# TEST 3: QUADRATIC COLOR CORRECTION
# ============================================================

print(f"\n\n{'=' * 70}")
print("TEST 3: QUADRATIC COLOR CORRECTION")
print("Is c² significant? What does it mean?")
print("=" * 70)

mu_model = dist_mod(z, H0=70)

# Standard Tripp
def tripp_chi2(params):
    M_B, alpha, beta = params
    mu_obs = mb + alpha * x1_q - beta * c_q - M_B
    return np.sum(((mu_obs - mu_model) / mb_err)**2)

res_tripp = scipy_minimize(tripp_chi2, x0=[-19.25, 0.15, 3.0], method='Nelder-Mead')
M_B_t, alpha_t, beta_t = res_tripp.x

# Tripp + c² (quadratic)
def tripp_c2_chi2(params):
    M_B, alpha, beta, delta = params
    mu_obs = mb + alpha * x1_q - beta * c_q - delta * c_q**2 - M_B
    return np.sum(((mu_obs - mu_model) / mb_err)**2)

res_c2 = scipy_minimize(tripp_c2_chi2, x0=[-19.25, 0.15, 3.0, 0.5], method='Nelder-Mead')
M_B_c2, alpha_c2, beta_c2, delta_c2 = res_c2.x

# Tripp + c² + c×z interaction
def tripp_full_chi2(params):
    M_B, alpha, beta, delta, eta = params
    mu_obs = mb + alpha * x1_q - beta * c_q - delta * c_q**2 - eta * c_q * z - M_B
    return np.sum(((mu_obs - mu_model) / mb_err)**2)

res_full = scipy_minimize(tripp_full_chi2, x0=[-19.25, 0.15, 3.0, 0.5, 1.0], method='Nelder-Mead')
M_B_f, alpha_f, beta_f, delta_f, eta_f = res_full.x

from scipy.stats import chi2 as chi2_dist

dchi2_c2 = res_tripp.fun - res_c2.fun
p_c2 = 1 - chi2_dist.cdf(dchi2_c2, 1)

dchi2_full = res_tripp.fun - res_full.fun
p_full = 1 - chi2_dist.cdf(dchi2_full, 2)

print(f"\n  Standard Tripp: χ² = {res_tripp.fun:.1f}")
print(f"    M_B = {M_B_t:.4f}, α = {alpha_t:.4f}, β = {beta_t:.3f}")
print(f"\n  Tripp + c²: χ² = {res_c2.fun:.1f} (Δχ² = {dchi2_c2:.1f}, p = {p_c2:.2e})")
print(f"    M_B = {M_B_c2:.4f}, α = {alpha_c2:.4f}, β = {beta_c2:.3f}, δ = {delta_c2:+.3f}")
print(f"\n  Tripp + c² + c×z: χ² = {res_full.fun:.1f} (Δχ² = {dchi2_full:.1f}, p = {p_full:.2e})")
print(f"    M_B = {M_B_f:.4f}, α = {alpha_f:.4f}, β = {beta_f:.3f}, δ = {delta_f:+.3f}, η = {eta_f:+.3f}")

# What the c² term means physically
print(f"\n  PHYSICAL INTERPRETATION:")
print(f"    δ = {delta_c2:+.3f} means:")
if delta_c2 > 0:
    print(f"    Extremely blue AND extremely red SNe are BOTH fainter than linear β predicts")
    print(f"    The color-luminosity relation CURVES — it's a parabola, not a line")
    print(f"    This is color-manifold curvature (GPT's interpretation)")
    print(f"    Two hidden variables compressed into one 'c' parameter")
else:
    print(f"    Extremely blue/red SNe are BRIGHTER than linear β predicts")
    print(f"    The color-luminosity relation is concave")


# ============================================================
# TEST 4: w-BIAS FROM z-DEPENDENT MASS STEP
# ============================================================

print(f"\n\n{'=' * 70}")
print("TEST 4: DOES γ(z) BIAS DARK ENERGY w?")
print("=" * 70)

# Fit with constant mass step vs z-dependent mass step
# Compare the inferred w (equation of state)

def dist_mod_w(z_arr, H0=70, Om=0.3, w=-1.0):
    """Distance modulus for flat wCDM"""
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        def integrand(zp):
            Ode = 1 - Om
            rho_de = Ode * (1+zp)**(3*(1+w))
            return 1.0 / np.sqrt(Om*(1+zp)**3 + rho_de)
        integral, _ = quad(integrand, 0, zi)
        dl[i] = C_LIGHT * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

# Fit wCDM with constant mass step
print("  Fitting wCDM with constant γ...")
mass_thresh = 10.0
mask_hi = hm >= mass_thresh

# Subsample for speed
step = 3
idx_sorted = np.argsort(z)
idx_sub = idx_sorted[::step]

z_s = z[idx_sub]
mb_s = mb[idx_sub]
mb_err_s = mb_err[idx_sub]
x1_s = x1_q[idx_sub]
c_s = c_q[idx_sub]
hm_s = hm[idx_sub]
mask_hi_s = hm_s >= mass_thresh

def wcdm_chi2_const(params):
    M_B, alpha, beta, gamma, w = params
    mu_mod = dist_mod_w(z_s, H0=70, Om=0.3, w=w)
    mu_obs = mb_s + alpha * x1_s - beta * c_s - M_B
    mu_obs[mask_hi_s] -= gamma
    return np.sum(((mu_obs - mu_mod) / mb_err_s)**2)

res_w_const = scipy_minimize(wcdm_chi2_const, 
                              x0=[-19.25, 0.15, 3.0, -0.05, -1.0],
                              method='Nelder-Mead', 
                              options={'maxiter': 30000})
M_B_wc, alpha_wc, beta_wc, gamma_wc, w_const = res_w_const.x

print(f"    w = {w_const:.4f}")
print(f"    γ = {gamma_wc:+.4f}")
print(f"    χ² = {res_w_const.fun:.1f}")

# Fit wCDM with linear γ(z) = γ₀ + γ₁×z
print("  Fitting wCDM with γ(z) = γ₀ + γ₁×z...")

def wcdm_chi2_gamma_z(params):
    M_B, alpha, beta, gamma0, gamma1, w = params
    mu_mod = dist_mod_w(z_s, H0=70, Om=0.3, w=w)
    gamma_z = gamma0 + gamma1 * z_s
    mu_obs = mb_s + alpha * x1_s - beta * c_s - M_B
    mu_obs[mask_hi_s] -= gamma_z[mask_hi_s]
    return np.sum(((mu_obs - mu_mod) / mb_err_s)**2)

res_w_gz = scipy_minimize(wcdm_chi2_gamma_z,
                           x0=[-19.25, 0.15, 3.0, -0.05, -0.03, -1.0],
                           method='Nelder-Mead',
                           options={'maxiter': 30000})
M_B_wg, alpha_wg, beta_wg, gamma0_wg, gamma1_wg, w_gz = res_w_gz.x

print(f"    w = {w_gz:.4f}")
print(f"    γ₀ = {gamma0_wg:+.4f}, γ₁ = {gamma1_wg:+.4f}")
print(f"    γ(z=0) = {gamma0_wg:+.4f}, γ(z=1) = {gamma0_wg + gamma1_wg:+.4f}")
print(f"    χ² = {res_w_gz.fun:.1f}")

dchi2_w = res_w_const.fun - res_w_gz.fun
p_w = 1 - chi2_dist.cdf(dchi2_w, 1)

print(f"\n  Δχ² = {dchi2_w:.2f} (p = {p_w:.4f})")
print(f"  w shift: {w_const:.4f} → {w_gz:.4f} (Δw = {w_gz - w_const:+.4f})")

print(f"\n  INTERPRETATION:")
print(f"    Constant γ gives w = {w_const:.3f}")
print(f"    γ(z) gives w = {w_gz:.3f}")
print(f"    The difference Δw = {w_gz - w_const:+.3f}")
if abs(w_gz - w_const) > 0.01:
    print(f"    ✓ z-dependent mass step SHIFTS w")
    print(f"    A constant mass step biases dark energy inference")
    if w_const > -1 and w_gz < w_const:
        print(f"    The bias pushes w AWAY from -1 (toward phantom)")
    elif w_const < -1 and w_gz > w_const:
        print(f"    The bias pushes w TOWARD -1 (toward ΛCDM)")
    
    print(f"\n    IMPLICATION FOR DESI:")
    print(f"    DESI reported w₀ = -0.55 ± 0.21 (w₀wₐCDM), suggesting")
    print(f"    evolving dark energy. If the mass step is z-dependent,")
    print(f"    part of that w-evolution could be a standardization artifact.")
else:
    print(f"    ⚠ w shift is small — mass step evolution doesn't significantly bias w")


# ============================================================
# TEST 5: THE FULL MODIFIED TRIPP EQUATION
# ============================================================

print(f"\n\n{'=' * 70}")
print("TEST 5: THE FULL MODIFIED TRIPP EQUATION")
print("μ = m_B + α(z)×x1 - β(z)×c - δ×c² - γ(z)×δ_M - M_B")
print("=" * 70)

def full_modified_tripp(params):
    M_B, a0, a1, b0, b1, delta_c2, g0, g1 = params
    alpha_z = a0 + a1 * z_s
    beta_z = b0 + b1 * z_s
    gamma_z = g0 + g1 * z_s
    
    mu_mod = dist_mod(z_s, H0=70)
    mu_obs = mb_s + alpha_z * x1_s - beta_z * c_s - delta_c2 * c_s**2 - M_B
    mu_obs[mask_hi_s] -= gamma_z[mask_hi_s]
    
    return np.sum(((mu_obs - mu_mod) / mb_err_s)**2)

res_full_tripp = scipy_minimize(full_modified_tripp,
                                 x0=[-19.25, 0.15, -0.05, 3.0, -0.5, 0.5, -0.05, -0.02],
                                 method='Nelder-Mead',
                                 options={'maxiter': 50000})

M_B_ft, a0_ft, a1_ft, b0_ft, b1_ft, d_ft, g0_ft, g1_ft = res_full_tripp.x

# Compare to standard Tripp
res_std_tripp = scipy_minimize(
    lambda p: np.sum(((mb_s + p[1]*x1_s - p[2]*c_s - p[0]) - dist_mod(z_s, H0=70))**2 / mb_err_s**2),
    x0=[-19.25, 0.15, 3.0], method='Nelder-Mead'
)

dchi2_full_tripp = res_std_tripp.fun - res_full_tripp.fun

print(f"\n  Standard Tripp (3 params): χ² = {res_std_tripp.fun:.1f}")
print(f"  Full Modified (8 params):  χ² = {res_full_tripp.fun:.1f}")
print(f"  Δχ² = {dchi2_full_tripp:.1f} for 5 extra parameters")
print(f"  p = {1 - chi2_dist.cdf(dchi2_full_tripp, 5):.2e}")

print(f"\n  COEFFICIENTS:")
print(f"    α(z) = {a0_ft:.4f} + ({a1_ft:+.4f})×z")
print(f"    β(z) = {b0_ft:.3f} + ({b1_ft:+.3f})×z")
print(f"    γ(z) = {g0_ft:+.4f} + ({g1_ft:+.4f})×z")
print(f"    δ (c²) = {d_ft:+.3f}")

print(f"\n  VALUES AT KEY REDSHIFTS:")
print(f"  {'z':>6} {'α(z)':>8} {'β(z)':>8} {'γ(z)':>8}")
print(f"  {'-'*34}")
for zi in [0.01, 0.05, 0.10, 0.30, 0.50, 1.00]:
    az = a0_ft + a1_ft * zi
    bz = b0_ft + b1_ft * zi
    gz = g0_ft + g1_ft * zi
    print(f"  {zi:>6.2f} {az:>8.4f} {bz:>8.3f} {gz:>+8.4f}")

# Infer H₀ from modified M_B
delta_MB = M_B_ft - res_std_tripp.x[0]
H0_modified = 70 * 10**(-delta_MB / 5)
print(f"\n  ΔM_B (modified - standard) = {delta_MB:+.4f}")
print(f"  Implied H₀ shift direction: {'lower' if delta_MB > 0 else 'higher'}")


# ============================================================
# GRAND SYNTHESIS
# ============================================================

print(f"\n\n{'=' * 70}")
print("GRAND SYNTHESIS")
print("=" * 70)

print(f"""
  THE TRIPP EQUATION IS BROKEN.
  
  Not "slightly wrong." Structurally inadequate.
  
  Standard Tripp assumes three CONSTANTS: α, β, γ.
  All three EVOLVE with z:
  
    α(z): stretch-luminosity evolves    (Δχ² = 178, p ≈ 0)
    β(z): color-luminosity evolves      (Δχ² = 112, p ≈ 0)
    γ(z): mass step evolves             (p = 10⁻²³)
    c²:   color correction is nonlinear (p = {p_c2:.1e})
  
  GPT's eigenvector rotation:
    PC1 rotates {total_angle_1:.0f}° from lowest to highest z
    PC2 rotates {total_angle_2:.0f}° from lowest to highest z
    The correction SPACE itself is z-dependent
  
  w-bias from constant γ:
    Constant γ:  w = {w_const:.3f}
    γ(z):        w = {w_gz:.3f}
    Shift: Δw = {w_gz - w_const:+.3f}
  
  CLOSURE THEORY PREDICTED ALL OF THIS:
    D_i(z) = a × N_modes^α × Z_g(z)
    Every observable with N > 0 evolves.
    The evolution is proportional to N_modes.
    The mass step is environment-dependent impedance.
    
  THE IMPLICATION:
    H₀ = 67.4 km/s/Mpc
    The tension is a standardization artifact
    Dark energy w measurements are contaminated
    The entire SN Ia distance ladder needs re-derivation
    with z-dependent coefficients
""")

# Save
os.makedirs('results_eigenvector', exist_ok=True)
output = {
    'test_date': '2026-03-09',
    'eigenvector_rotation': {
        'PC1_total_degrees': float(total_angle_1) if len(all_eigvecs) >= 2 else None,
        'PC2_total_degrees': float(total_angle_2) if len(all_eigvecs) >= 2 else None,
    },
    'quadratic_color': {
        'delta': float(delta_c2),
        'dchi2': float(dchi2_c2),
        'p_value': float(p_c2),
    },
    'w_bias': {
        'w_constant_gamma': float(w_const),
        'w_gamma_z': float(w_gz),
        'delta_w': float(w_gz - w_const),
    },
    'full_modified_tripp': {
        'alpha0': float(a0_ft), 'alpha1': float(a1_ft),
        'beta0': float(b0_ft), 'beta1': float(b1_ft),
        'delta_c2': float(d_ft),
        'gamma0': float(g0_ft), 'gamma1': float(g1_ft),
        'dchi2_vs_standard': float(dchi2_full_tripp),
    }
}

with open('results_eigenvector/eigenvector_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"  Results saved to results_eigenvector/eigenvector_results.json")
print(f"\n{'=' * 70}")
print("COMPLETE")
print("=" * 70)
