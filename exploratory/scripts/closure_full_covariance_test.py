#!/usr/bin/env python3
"""
CLOSURE THEORY — TEST 1: FULL COVARIANCE MATRIX FIT
======================================================

Use the official Pantheon+ covariance matrix (stat + sys) instead of
diagonal-only errors. Refit standard Tripp and modified Tripp.
Then run mock null injection with the FULL covariance to measure
the TRUE systematic floor.

Key question: Does the off-diagonal structure create ρ ~ -0.10 to -0.13
in the null distribution? If yes, the decode is complete.

Author: Closure Theory collaboration
Date: 2026-03-09
"""

import numpy as np
from scipy.optimize import minimize as scipy_minimize
from scipy.integrate import quad
from scipy.stats import spearmanr, percentileofscore
from scipy.linalg import cho_factor, cho_solve
import json
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# LOAD PANTHEON+ DATA (pre-standardized distance moduli)
# ============================================================

# The DES/Pantheon+ public release has pre-standardized μ values
# The covariance is for these μ values, not raw mB
# So we need to use the Pantheon+ mu file that matches the covariance

# The covariance.cov file is 1701x1701 — this matches the FULL Pantheon+ sample
# First line is N, then N*N values

print("Loading covariance matrix (1701x1701)...")
with open('data/covariance.cov', 'r') as f:
    N_cov = int(f.readline().strip())
    cov_flat = np.array([float(f.readline().strip()) for _ in range(N_cov * N_cov)])

C_full = cov_flat.reshape(N_cov, N_cov)
print(f"  Covariance: {N_cov} x {N_cov}")
print(f"  Diagonal range: [{np.min(np.diag(C_full)):.6f}, {np.max(np.diag(C_full)):.6f}]")
print(f"  Diagonal mean: {np.mean(np.diag(C_full)):.6f}")
print(f"  Off-diagonal fraction: {np.mean(np.abs(C_full[np.triu_indices(N_cov, 1)])) / np.mean(np.diag(C_full)):.4f}")

# ============================================================
# LOAD RAW PANTHEON+ DATA (with x1, c, host mass)
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

z_all = np.array([float(r[col['zCMB']]) for r in rows])
z_hd_all = np.array([float(r[col['zHD']]) for r in rows])
m_b_all = np.array([float(r[col['mB']]) for r in rows])
m_b_err_all = np.array([float(r[col['mBERR']]) for r in rows])
x1_all = np.array([float(r[col['x1']]) for r in rows])
c_all = np.array([float(r[col['c']]) for r in rows])
hm_all = np.array([float(r[col['HOST_LOGMASS']]) for r in rows])

print(f"  Raw data: {len(z_all)} SNe")

# Quality cuts (same as before)
mask_q = (z_hd_all > 0.01) & (z_hd_all < 2.5) & (np.abs(x1_all) < 5) & (np.abs(c_all) < 0.5) & (m_b_err_all < 1.0) & (hm_all > 0) & (hm_all < 15)

# We need to track which indices survive the cut for covariance slicing
idx_keep = np.where(mask_q)[0]
N_keep = len(idx_keep)

z = z_all[mask_q]
mb = m_b_all[mask_q]
mb_err = m_b_err_all[mask_q]
x1_q = x1_all[mask_q]
c_q = c_all[mask_q]
hm = hm_all[mask_q]
mask_hi = hm >= 10.0

print(f"  After cuts: {N_keep} SNe")

# Slice covariance to match our quality cuts
# The covariance is for ALL 1701 SNe; we need the submatrix
if N_cov == len(z_all):
    C_sub = C_full[np.ix_(idx_keep, idx_keep)]
    print(f"  Covariance submatrix: {C_sub.shape}")
else:
    print(f"  WARNING: Covariance size {N_cov} != data size {len(z_all)}")
    print(f"  Using diagonal only as fallback")
    C_sub = np.diag(mb_err**2)

# Check if covariance is for mB or mu
# If diagonal ≈ mBERR², it's for mB
diag_ratio = np.mean(np.sqrt(np.diag(C_sub))) / np.mean(mb_err)
print(f"  √(diag) / mBERR ratio: {diag_ratio:.3f}")
if abs(diag_ratio - 1.0) < 0.5:
    print(f"  → Covariance appears to be for mB (matches mBERR)")
else:
    print(f"  → Covariance may be for μ or different units")

# ============================================================
# COSMOLOGY
# ============================================================

C_LIGHT = 299792.458

def E(zp, Om=0.3):
    return np.sqrt(Om*(1+zp)**3 + (1-Om))

def dist_mod(z_arr, H0=70, Om=0.3):
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        integral, _ = quad(lambda zp: 1.0/E(zp, Om), 0, zi)
        dl[i] = C_LIGHT * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

print("Computing distance moduli...")
mu_lcdm = dist_mod(z, H0=70, Om=0.3)

I_sat = z / (1 + z/0.5)

# Cholesky decomposition for efficient likelihood evaluation
print("Cholesky factorization...")
try:
    cho_C = cho_factor(C_sub)
    use_full_cov = True
    print("  ✓ Full covariance ready")
except Exception as e:
    print(f"  ✗ Cholesky failed: {e}")
    print(f"  Falling back to diagonal")
    use_full_cov = False


# ============================================================
# FIT WITH FULL COVARIANCE
# ============================================================

def chi2_full_cov(resid):
    """χ² using full covariance: r^T C^-1 r"""
    if use_full_cov:
        return float(resid @ cho_solve(cho_C, resid))
    else:
        return float(np.sum((resid / mb_err)**2))

def chi2_diag(resid):
    """χ² using diagonal only"""
    return float(np.sum((resid / mb_err)**2))


# --- Standard Tripp ---
print(f"\n{'=' * 70}")
print("STANDARD TRIPP — FULL COVARIANCE vs DIAGONAL")
print("=" * 70)

def fit_std(chi2_func):
    def objective(params):
        M_B, alpha, beta, gamma = params
        mu_obs = mb + alpha * x1_q - beta * c_q - M_B
        mu_obs[mask_hi] -= gamma
        resid = mu_obs - mu_lcdm
        return chi2_func(resid)
    
    res = scipy_minimize(objective, x0=[-19.25, 0.15, 3.0, -0.05],
                         method='Nelder-Mead', options={'maxiter': 50000})
    p = res.x
    resid = (mb + p[1]*x1_q - p[2]*c_q - p[0]) - mu_lcdm
    resid[mask_hi] -= p[3]
    rho, pval = spearmanr(z, resid)
    return res.fun, rho, pval, p

chi2_std_full, rho_std_full, p_std_full, params_std_full = fit_std(chi2_full_cov)
chi2_std_diag, rho_std_diag, p_std_diag, params_std_diag = fit_std(chi2_diag)

print(f"  FULL COV: χ²={chi2_std_full:.1f}, χ²/ν={chi2_std_full/(N_keep-4):.4f}, ρ={rho_std_full:+.4f} (p={p_std_full:.2e})")
print(f"  DIAGONAL: χ²={chi2_std_diag:.1f}, χ²/ν={chi2_std_diag/(N_keep-4):.4f}, ρ={rho_std_diag:+.4f} (p={p_std_diag:.2e})")


# --- Modified Tripp (M1) ---
print(f"\n{'=' * 70}")
print("MODIFIED TRIPP (M1) — FULL COVARIANCE vs DIAGONAL")
print("=" * 70)

x1_std = np.std(x1_q)
c_std = np.std(c_q)

def fit_m1(chi2_func):
    def objective(params):
        A, B, theta1, delta_v, g0, g1 = params
        angle = theta1 * I_sat
        x1n = x1_q / x1_std; cn = c_q / c_std
        u = x1n * np.cos(angle) + cn * np.sin(angle)
        v = -x1n * np.sin(angle) + cn * np.cos(angle)
        u_mag = u * x1_std; v_mag = v * c_std
        gamma_I = g0 + g1 * I_sat
        mu_corr = mb + A * u_mag - B * v_mag - delta_v * v_mag**2
        mu_corr[mask_hi] -= gamma_I[mask_hi]
        # Analytically marginalize M₀
        if use_full_cov and chi2_func == chi2_full_cov:
            ones = np.ones(N_keep)
            Cinv_ones = cho_solve(cho_C, ones)
            Cinv_r = cho_solve(cho_C, mu_corr - mu_lcdm)
            M0 = np.dot(ones, Cinv_r) / np.dot(ones, Cinv_ones)
        else:
            w = 1.0 / mb_err**2
            M0 = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
        resid = mu_corr - M0 - mu_lcdm
        return chi2_func(resid)
    
    res = scipy_minimize(objective, x0=[0.15, 3.0, -0.3, 2.0, -0.03, -0.1],
                         method='Nelder-Mead', options={'maxiter': 100000})
    
    # Compute residuals at best fit
    A, B, theta1, delta_v, g0, g1 = res.x
    angle = theta1 * I_sat
    x1n = x1_q / x1_std; cn = c_q / c_std
    u = x1n * np.cos(angle) + cn * np.sin(angle)
    v = -x1n * np.sin(angle) + cn * np.cos(angle)
    u_mag = u * x1_std; v_mag = v * c_std
    gamma_I = g0 + g1 * I_sat
    mu_corr = mb + A * u_mag - B * v_mag - delta_v * v_mag**2
    mu_corr[mask_hi] -= gamma_I[mask_hi]
    
    if use_full_cov and chi2_func == chi2_full_cov:
        ones = np.ones(N_keep)
        Cinv_ones = cho_solve(cho_C, ones)
        Cinv_r = cho_solve(cho_C, mu_corr - mu_lcdm)
        M0 = np.dot(ones, Cinv_r) / np.dot(ones, Cinv_ones)
    else:
        w = 1.0 / mb_err**2
        M0 = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
    
    resid = mu_corr - M0 - mu_lcdm
    rho, pval = spearmanr(z, resid)
    return res.fun, rho, pval, res.x

chi2_m1_full, rho_m1_full, p_m1_full, params_m1_full = fit_m1(chi2_full_cov)
chi2_m1_diag, rho_m1_diag, p_m1_diag, params_m1_diag = fit_m1(chi2_diag)

print(f"  FULL COV: χ²={chi2_m1_full:.1f}, χ²/ν={chi2_m1_full/(N_keep-7):.4f}, ρ={rho_m1_full:+.4f} (p={p_m1_full:.2e})")
print(f"  DIAGONAL: χ²={chi2_m1_diag:.1f}, χ²/ν={chi2_m1_diag/(N_keep-7):.4f}, ρ={rho_m1_diag:+.4f} (p={p_m1_diag:.2e})")


# ============================================================
# MOCK NULL WITH FULL COVARIANCE (GPT's upgrade)
# ============================================================

N_MOCKS = 200  # Fewer mocks since full-cov fitting is slower

print(f"\n{'=' * 70}")
print(f"MOCK NULL INJECTION — FULL COVARIANCE ({N_MOCKS} mocks)")
print("=" * 70)

# Generate mocks from FULL covariance (multivariate normal)
print("Generating multivariate normal draws from full covariance...")

# True mB under ΛCDM + standard Tripp
M_B_true = params_std_full[0]
alpha_true = params_std_full[1]
beta_true = params_std_full[2]
gamma_true = params_std_full[3]

mb_true = mu_lcdm + M_B_true - alpha_true * x1_q + beta_true * c_q
mb_true[mask_hi] += gamma_true

# Cholesky of covariance for generating correlated noise
try:
    L = np.linalg.cholesky(C_sub)
    can_generate = True
    print("  ✓ Cholesky decomposition for mock generation ready")
except:
    print("  ✗ Cholesky failed for generation — adding small regularization")
    C_reg = C_sub + 1e-6 * np.eye(N_keep)
    try:
        L = np.linalg.cholesky(C_reg)
        can_generate = True
        print("  ✓ Regularized Cholesky ready")
    except:
        can_generate = False
        print("  ✗ Cannot generate correlated mocks — using diagonal")

rho_std_mocks = np.zeros(N_MOCKS)
rho_m1_mocks = np.zeros(N_MOCKS)

for i in range(N_MOCKS):
    if (i+1) % 50 == 0:
        print(f"  Mock {i+1}/{N_MOCKS}...")
    
    # Generate correlated noise
    if can_generate:
        noise = L @ np.random.normal(0, 1, N_keep)
    else:
        noise = np.random.normal(0, mb_err)
    
    mb_mock = mb_true + noise
    
    # Fit standard Tripp with DIAGONAL (same as real analysis)
    def chi2_s(params):
        M_B, alpha, beta, gamma = params
        mu_obs = mb_mock + alpha * x1_q - beta * c_q - M_B
        mu_obs[mask_hi] -= gamma
        return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)
    
    res_s = scipy_minimize(chi2_s, x0=[-19.25, 0.15, 3.0, -0.05],
                           method='Nelder-Mead', options={'maxiter': 20000})
    p_s = res_s.x
    resid_s = (mb_mock + p_s[1]*x1_q - p_s[2]*c_q - p_s[0]) - mu_lcdm
    resid_s[mask_hi] -= p_s[3]
    rho_std_mocks[i], _ = spearmanr(z, resid_s)
    
    # Fit modified M1 with DIAGONAL
    def chi2_m(params):
        A, B, theta1, delta_v, g0, g1 = params
        angle = theta1 * I_sat
        x1n = x1_q / x1_std; cn = c_q / c_std
        u = x1n * np.cos(angle) + cn * np.sin(angle)
        v = -x1n * np.sin(angle) + cn * np.cos(angle)
        u_mag = u * x1_std; v_mag = v * c_std
        gamma_I = g0 + g1 * I_sat
        mu_corr = mb_mock + A * u_mag - B * v_mag - delta_v * v_mag**2
        mu_corr[mask_hi] -= gamma_I[mask_hi]
        w = 1.0 / mb_err**2
        M0 = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
        residuals = mu_corr - M0 - mu_lcdm
        return np.sum((residuals / mb_err)**2)
    
    res_m = scipy_minimize(chi2_m, x0=[0.15, 3.0, -0.3, 2.0, -0.03, -0.1],
                           method='Nelder-Mead', options={'maxiter': 30000})
    A, B, theta1, delta_v, g0, g1 = res_m.x
    angle = theta1 * I_sat
    x1n = x1_q / x1_std; cn = c_q / c_std
    u = x1n * np.cos(angle) + cn * np.sin(angle)
    v = -x1n * np.sin(angle) + cn * np.cos(angle)
    u_mag = u * x1_std; v_mag = v * c_std
    gamma_I = g0 + g1 * I_sat
    mu_corr = mb_mock + A * u_mag - B * v_mag - delta_v * v_mag**2
    mu_corr[mask_hi] -= gamma_I[mask_hi]
    w = 1.0 / mb_err**2
    M0 = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
    resid_m = mu_corr - M0 - mu_lcdm
    rho_m1_mocks[i], _ = spearmanr(z, resid_m)


# ============================================================
# ANALYSIS — FULL COVARIANCE MOCK NULL
# ============================================================

print(f"\n{'=' * 70}")
print("FULL COVARIANCE MOCK NULL — RESULTS")
print("=" * 70)

# Standard Tripp
mean_s = np.mean(rho_std_mocks)
std_s = np.std(rho_std_mocks)
sigma_pre = (rho_std_full - mean_s) / std_s
pct_pre = percentileofscore(rho_std_mocks, rho_std_full)

print(f"\n  STANDARD TRIPP (correlated mocks, diagonal fit):")
print(f"    Null ρ: mean = {mean_s:+.4f}, σ = {std_s:.4f}")
print(f"    Null ρ range: [{np.min(rho_std_mocks):+.4f}, {np.max(rho_std_mocks):+.4f}]")
print(f"    Real ρ = {rho_std_full:+.4f} ({sigma_pre:+.1f}σ from null)")

# Modified Tripp
mean_m = np.mean(rho_m1_mocks)
std_m = np.std(rho_m1_mocks)
sigma_post = (rho_m1_full - mean_m) / std_m
pct_post = percentileofscore(rho_m1_mocks, rho_m1_full)

print(f"\n  MODIFIED TRIPP M1 (correlated mocks, diagonal fit):")
print(f"    Null ρ: mean = {mean_m:+.4f}, σ = {std_m:.4f}")
print(f"    Null ρ range: [{np.min(rho_m1_mocks):+.4f}, {np.max(rho_m1_mocks):+.4f}]")
print(f"    Real ρ = {rho_m1_full:+.4f} ({sigma_post:+.1f}σ from null)")

# The KEY question: does correlated noise shift the null?
print(f"\n  KEY COMPARISON — DIAGONAL vs FULL-COV NULL:")
print(f"    Diagonal null (std):  mean = -0.0006, σ = 0.0243")
print(f"    Full-cov null (std):  mean = {mean_s:+.4f}, σ = {std_s:.4f}")
print(f"    Shift: {mean_s - (-0.0006):+.4f}")
if abs(mean_s) > 0.05:
    print(f"    ★ OFF-DIAGONAL STRUCTURE CREATES SYSTEMATIC ρ BIAS")
else:
    print(f"    ○ Off-diagonal structure does NOT create large systematic bias")


# ============================================================
# VERDICT
# ============================================================

print(f"\n{'=' * 70}")
print("VERDICT — FULL COVARIANCE TEST")
print("=" * 70)

pre_anomalous = abs(sigma_pre) > 3
post_inside = abs(sigma_post) < 2

if pre_anomalous and post_inside:
    print(f"""
  ★★★ FRAMEWORK VALIDATED WITH FULL COVARIANCE ★★★
  Pre-correction: {sigma_pre:+.1f}σ from null (ANOMALOUS)
  Post-correction: {sigma_post:+.1f}σ from null (INSIDE)
  THE DECODE IS COMPLETE.
""")
elif pre_anomalous and abs(sigma_post) < 3:
    print(f"""
  ◐ STRONG VALIDATION — MARGINAL DECODE
  Pre-correction: {sigma_pre:+.1f}σ (ANOMALOUS)
  Post-correction: {sigma_post:+.1f}σ (marginal, between 2-3σ)
  Framework captures dominant signal; residual is at boundary.
""")
elif pre_anomalous:
    print(f"""
  ◐ SIGNAL CONFIRMED, DECODE INCOMPLETE
  Pre-correction: {sigma_pre:+.1f}σ (ANOMALOUS)
  Post-correction: {sigma_post:+.1f}σ (still outside)
  Real signal exists but correction doesn't reach the floor.
""")
else:
    print(f"""
  ⚠ UNEXPECTED: Pre-correction inside full-cov null
  This would mean the "signal" is an artifact of correlated noise.
""")


# Save
os.makedirs('results_full_covariance', exist_ok=True)
output = {
    'test_date': '2026-03-09',
    'n_sne': int(N_keep),
    'covariance_size': int(N_cov),
    'full_cov_fits': {
        'std_chi2': float(chi2_std_full), 'std_rho': float(rho_std_full),
        'm1_chi2': float(chi2_m1_full), 'm1_rho': float(rho_m1_full),
    },
    'diag_fits': {
        'std_chi2': float(chi2_std_diag), 'std_rho': float(rho_std_diag),
        'm1_chi2': float(chi2_m1_diag), 'm1_rho': float(rho_m1_diag),
    },
    'full_cov_null': {
        'n_mocks': N_MOCKS,
        'std_mean': float(mean_s), 'std_sigma': float(std_s),
        'std_sigma_from_null': float(sigma_pre),
        'm1_mean': float(mean_m), 'm1_sigma': float(std_m),
        'm1_sigma_from_null': float(sigma_post),
    }
}

with open('results_full_covariance/full_cov_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n  Results saved to results_full_covariance/full_cov_results.json")
print(f"\n{'=' * 70}")
print("COMPLETE")
print("=" * 70)
