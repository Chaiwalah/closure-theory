#!/usr/bin/env python3
"""
CLOSURE THEORY — STRETCH-ONLY H₀ TEST
========================================

THE CRITICAL TEST (proposed by GPT, Gemini, Grok — all agree):

If the Hubble tension is partly a metric impedance artifact,
and color (c) is the impedance-bearing channel while stretch (x1)
is relatively protected, then:

1. H₀ from stretch-only standardization should be LOWER than stretch+color
2. The difference should be larger for local SNe (inside KBC void)
3. The difference should be smaller for distant SNe (outside KBC void)

We use the Pantheon+ SN Ia sample with SALT2 light-curve parameters.

TRIPP STANDARDIZATION:
    μ = m_B - M_B + α×x1 - β×c

STRETCH-ONLY:
    μ = m_B - M_B + α×x1

The difference in inferred H₀ between these two = the color channel's
contribution to the Hubble tension.

Author: Closure Theory collaboration
Date: 2026-03-09
"""

import numpy as np
from scipy.optimize import minimize
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# LOAD PANTHEON+ DATA
# ============================================================

print("=" * 70)
print("STRETCH-ONLY H₀ TEST — THE CRITICAL EXPERIMENT")
print("=" * 70)

data_file = 'data/pantheon_plus.dat'

# Parse the space-separated file
with open(data_file, 'r') as f:
    header = f.readline().strip().split()
    rows = []
    for line in f:
        parts = line.strip().split()
        if len(parts) == len(header):
            rows.append(parts)

print(f"\n  Loaded {len(rows)} SN Ia entries from Pantheon+")
print(f"  Columns: {len(header)}")

# Extract needed columns
col = {name: i for i, name in enumerate(header)}

# Build arrays
cid = [r[col['CID']] for r in rows]
z_hd = np.array([float(r[col['zHD']]) for r in rows])
z_cmb = np.array([float(r[col['zCMB']]) for r in rows])
m_b = np.array([float(r[col['mB']]) for r in rows])
m_b_err = np.array([float(r[col['mBERR']]) for r in rows])
x1 = np.array([float(r[col['x1']]) for r in rows])
x1_err = np.array([float(r[col['x1ERR']]) for r in rows])
c = np.array([float(r[col['c']]) for r in rows])
c_err = np.array([float(r[col['cERR']]) for r in rows])
m_b_corr = np.array([float(r[col['m_b_corr']]) for r in rows])

# Quality cuts: reasonable values
mask_quality = (z_hd > 0.01) & (z_hd < 2.5) & (np.abs(x1) < 5) & (np.abs(c) < 0.5) & (m_b_err < 1.0)
print(f"  After quality cuts: {np.sum(mask_quality)} SNe")

z = z_cmb[mask_quality]
mb = m_b[mask_quality]
mb_err = m_b_err[mask_quality]
x1_q = x1[mask_quality]
x1_err_q = x1_err[mask_quality]
c_q = c[mask_quality]
c_err_q = c_err[mask_quality]
mb_corr = m_b_corr[mask_quality]

# ============================================================
# COSMOLOGY FUNCTIONS
# ============================================================

def luminosity_distance(z, H0=70, Om=0.3):
    """Simplified flat ΛCDM luminosity distance in Mpc"""
    from scipy.integrate import quad
    c_light = 299792.458  # km/s
    
    def integrand(zp):
        return 1.0 / np.sqrt(Om * (1+zp)**3 + (1-Om))
    
    dl = np.zeros_like(z)
    for i, zi in enumerate(z):
        integral, _ = quad(integrand, 0, zi)
        dl[i] = c_light * (1+zi) * integral / H0
    return dl

def distance_modulus(z, H0=70, Om=0.3):
    """Distance modulus μ = 5 log10(dL/10pc)"""
    dl = luminosity_distance(z, H0, Om)
    return 5 * np.log10(dl) + 25

def fit_H0(z_data, mu_data, mu_err, Om=0.3):
    """Fit H0 by minimizing chi-squared"""
    from scipy.optimize import minimize_scalar
    
    def chi2(H0):
        mu_model = distance_modulus(z_data, H0, Om)
        return np.sum(((mu_data - mu_model) / mu_err)**2)
    
    result = minimize_scalar(chi2, bounds=(50, 90), method='bounded')
    
    # Error from Δχ² = 1
    chi2_min = result.fun
    H0_best = result.x
    
    # Find 1σ bounds
    from scipy.optimize import brentq
    try:
        H0_lo = brentq(lambda h: chi2(h) - chi2_min - 1, 50, H0_best)
        H0_hi = brentq(lambda h: chi2(h) - chi2_min - 1, H0_best, 90)
        sigma = (H0_hi - H0_lo) / 2
    except:
        sigma = 1.0
    
    return H0_best, sigma, chi2_min / (len(z_data) - 1)


# ============================================================
# TEST 1: FULL TRIPP vs STRETCH-ONLY (FULL SAMPLE)
# ============================================================

print(f"\n\n{'=' * 70}")
print("TEST 1: FULL TRIPP vs STRETCH-ONLY (ALL SNe)")
print("=" * 70)

# Fit Tripp parameters first
# μ = m_B + α×x1 - β×c - M_B
# We fit α, β, M_B simultaneously

def tripp_chi2(params, z_data, mb_data, x1_data, c_data, mb_err_data, Om=0.3):
    """Full Tripp standardization chi-squared"""
    M_B, alpha, beta = params
    mu_obs = mb_data + alpha * x1_data - beta * c_data - M_B
    mu_model = distance_modulus(z_data, H0=70, Om=Om)  # fix H0, float M_B
    # M_B and H0 are degenerate; we absorb H0 into M_B
    residuals = mu_obs - mu_model
    return np.sum((residuals / mb_err_data)**2)

# Use a subsample for speed (every 3rd SN, sorted by z)
idx_sorted = np.argsort(z)
step = 3
idx_sub = idx_sorted[::step]
z_sub = z[idx_sub]
mb_sub = mb[idx_sub]
x1_sub = x1_q[idx_sub]
c_sub = c_q[idx_sub]
mb_err_sub = mb_err[idx_sub]

n_sub = len(z_sub)
print(f"\n  Using subsample of {n_sub} SNe (every {step}th, sorted by z)")

# Fit full Tripp
from scipy.optimize import minimize as scipy_minimize

res_tripp = scipy_minimize(tripp_chi2, x0=[-19.3, 0.15, 3.1],
                           args=(z_sub, mb_sub, x1_sub, c_sub, mb_err_sub),
                           method='Nelder-Mead')
M_B_tripp, alpha_tripp, beta_tripp = res_tripp.x
chi2_tripp = res_tripp.fun / (n_sub - 3)

print(f"\n  FULL TRIPP FIT:")
print(f"    M_B   = {M_B_tripp:.4f}")
print(f"    α     = {alpha_tripp:.4f} (stretch coefficient)")
print(f"    β     = {beta_tripp:.4f} (color coefficient)")
print(f"    χ²/ν  = {chi2_tripp:.2f}")

# Distance moduli from Tripp
mu_tripp = mb_sub + alpha_tripp * x1_sub - beta_tripp * c_sub - M_B_tripp

# Fit stretch-only
def stretch_chi2(params, z_data, mb_data, x1_data, mb_err_data, Om=0.3):
    """Stretch-only standardization"""
    M_B, alpha = params
    mu_obs = mb_data + alpha * x1_data - M_B
    mu_model = distance_modulus(z_data, H0=70, Om=Om)
    residuals = mu_obs - mu_model
    return np.sum((residuals / mb_err_data)**2)

res_stretch = scipy_minimize(stretch_chi2, x0=[-19.3, 0.15],
                             args=(z_sub, mb_sub, x1_sub, mb_err_sub),
                             method='Nelder-Mead')
M_B_stretch, alpha_stretch = res_stretch.x
chi2_stretch = res_stretch.fun / (n_sub - 2)

print(f"\n  STRETCH-ONLY FIT:")
print(f"    M_B   = {M_B_stretch:.4f}")
print(f"    α     = {alpha_stretch:.4f}")
print(f"    χ²/ν  = {chi2_stretch:.2f}")

mu_stretch = mb_sub + alpha_stretch * x1_sub - M_B_stretch

# The M_B difference encodes the H₀ shift
# μ = m_B - M_B + corrections = 5 log10(dL/10pc)
# ΔM_B = -5 log10(H0_tripp/H0_stretch)
# H0_stretch/H0_tripp = 10^(ΔM_B/5)

delta_MB = M_B_stretch - M_B_tripp
H0_ratio = 10**(delta_MB / 5)

print(f"\n  COMPARISON:")
print(f"    ΔM_B = M_B(stretch) - M_B(tripp) = {delta_MB:+.4f}")
print(f"    H₀ ratio = 10^(ΔM_B/5) = {H0_ratio:.4f}")
print(f"    If Tripp gives H₀ = 73.0: stretch-only gives H₀ ≈ {73.0 * H0_ratio:.1f}")
print(f"    If Tripp gives H₀ = 73.6: stretch-only gives H₀ ≈ {73.6 * H0_ratio:.1f}")

direction = "LOWER" if H0_ratio < 1 else "HIGHER"
print(f"\n    Stretch-only H₀ is {direction} than Tripp H₀")
print(f"    Impedance prediction: should be LOWER ({'✓ CONFIRMED' if H0_ratio < 1 else '✗ VIOLATED'})")


# ============================================================
# TEST 2: INSIDE vs OUTSIDE KBC VOID
# ============================================================

print(f"\n\n{'=' * 70}")
print("TEST 2: KBC VOID SPLIT — INSIDE (z<0.07) vs OUTSIDE (z>0.07)")
print("=" * 70)

z_kbc = 0.07  # KBC void boundary

mask_inside = z_sub < z_kbc
mask_outside = z_sub >= z_kbc

n_inside = np.sum(mask_inside)
n_outside = np.sum(mask_outside)
print(f"\n  Inside KBC (z < {z_kbc}): {n_inside} SNe")
print(f"  Outside KBC (z ≥ {z_kbc}): {n_outside} SNe")

for label, mask in [("INSIDE KBC VOID (z < 0.07)", mask_inside), 
                     ("OUTSIDE KBC VOID (z ≥ 0.07)", mask_outside)]:
    if np.sum(mask) < 10:
        print(f"\n  {label}: Too few SNe ({np.sum(mask)}), skipping")
        continue
    
    z_m = z_sub[mask]
    mb_m = mb_sub[mask]
    x1_m = x1_sub[mask]
    c_m = c_sub[mask]
    err_m = mb_err_sub[mask]
    
    # Tripp fit
    res_t = scipy_minimize(tripp_chi2, x0=[-19.3, 0.15, 3.1],
                           args=(z_m, mb_m, x1_m, c_m, err_m),
                           method='Nelder-Mead')
    MB_t, a_t, b_t = res_t.x
    
    # Stretch-only fit
    res_s = scipy_minimize(stretch_chi2, x0=[-19.3, 0.15],
                           args=(z_m, mb_m, x1_m, err_m),
                           method='Nelder-Mead')
    MB_s, a_s = res_s.x
    
    dMB = MB_s - MB_t
    ratio = 10**(dMB / 5)
    
    print(f"\n  {label} ({np.sum(mask)} SNe):")
    print(f"    Tripp:   M_B={MB_t:.3f}, α={a_t:.3f}, β={b_t:.3f}")
    print(f"    Stretch: M_B={MB_s:.3f}, α={a_s:.3f}")
    print(f"    ΔM_B = {dMB:+.4f}, H₀ ratio = {ratio:.4f}")
    print(f"    Color contribution to H₀: {(1-ratio)*100:+.1f}%")


# ============================================================
# TEST 3: z-BINNED COLOR CONTRIBUTION
# ============================================================

print(f"\n\n{'=' * 70}")
print("TEST 3: COLOR CONTRIBUTION TO H₀ IN z-BINS")
print("=" * 70)

# Split into z-bins and measure the color channel's H₀ contribution in each
z_bins = [(0.01, 0.05), (0.05, 0.10), (0.10, 0.20), (0.20, 0.40), (0.40, 0.80)]

print(f"\n  {'z-bin':>15} {'N':>5} {'ΔM_B':>8} {'H₀ ratio':>10} {'Color ΔH₀':>10} {'In void?':>10}")
print(f"  {'-'*65}")

dMB_bins = []
z_mid_bins = []

for z_lo, z_hi in z_bins:
    mask_bin = (z_sub >= z_lo) & (z_sub < z_hi)
    n_bin = np.sum(mask_bin)
    
    if n_bin < 10:
        print(f"  [{z_lo:.2f}, {z_hi:.2f}]  {n_bin:>5}  (insufficient data)")
        continue
    
    z_b = z_sub[mask_bin]
    mb_b = mb_sub[mask_bin]
    x1_b = x1_sub[mask_bin]
    c_b = c_sub[mask_bin]
    err_b = mb_err_sub[mask_bin]
    
    # Tripp
    try:
        res_t = scipy_minimize(tripp_chi2, x0=[-19.3, 0.15, 3.1],
                               args=(z_b, mb_b, x1_b, c_b, err_b),
                               method='Nelder-Mead')
        MB_t = res_t.x[0]
    except:
        continue
    
    # Stretch-only
    try:
        res_s = scipy_minimize(stretch_chi2, x0=[-19.3, 0.15],
                               args=(z_b, mb_b, x1_b, err_b),
                               method='Nelder-Mead')
        MB_s = res_s.x[0]
    except:
        continue
    
    dMB = MB_s - MB_t
    ratio = 10**(dMB / 5)
    z_mid = (z_lo + z_hi) / 2
    in_void = "YES" if z_mid < 0.07 else "partial" if z_mid < 0.15 else "NO"
    
    # Color contribution to H₀ (in km/s/Mpc assuming baseline 73)
    color_dH0 = 73.0 * (1 - ratio)
    
    print(f"  [{z_lo:.2f}, {z_hi:.2f}]  {n_bin:>5}  {dMB:>+8.4f}  {ratio:>10.4f}  {color_dH0:>+10.2f}  {in_void:>10}")
    
    dMB_bins.append(dMB)
    z_mid_bins.append(z_mid)

if len(dMB_bins) >= 3:
    dMB_arr = np.array(dMB_bins)
    z_mid_arr = np.array(z_mid_bins)
    
    rho_z, p_z = spearmanr(z_mid_arr, dMB_arr)
    print(f"\n  Color ΔM_B vs z: Spearman ρ = {rho_z:+.3f} (p = {p_z:.4f})")
    
    print(f"\n  IMPEDANCE PREDICTION:")
    print(f"  Color contribution should be LARGER at low z (inside void)")
    print(f"  and SMALLER at high z (outside void)")
    print(f"  → ΔM_B should decrease (become less negative) with z")
    

# ============================================================
# TEST 4: HUBBLE RESIDUALS — COLOR vs STRETCH
# ============================================================

print(f"\n\n{'=' * 70}")
print("TEST 4: HUBBLE RESIDUALS — COLOR TERM vs STRETCH TERM")
print("=" * 70)

# Using best-fit Tripp parameters, decompose the correction:
# μ_obs = m_B + α×x1 - β×c - M_B
# color_correction = -β×c
# stretch_correction = α×x1
# 
# If color carries impedance, color_correction should correlate
# with z differently than stretch_correction

color_corr = -beta_tripp * c_sub
stretch_corr = alpha_tripp * x1_sub

# Bin by z and measure mean correction in each bin
print(f"\n  Using global Tripp fit: α={alpha_tripp:.3f}, β={beta_tripp:.3f}")
print()
print(f"  {'z-bin':>15} {'N':>5} {'⟨color_corr⟩':>14} {'⟨stretch_corr⟩':>16} {'⟨color⟩/⟨stretch⟩':>18}")
print(f"  {'-'*72}")

for z_lo, z_hi in z_bins:
    mask_bin = (z_sub >= z_lo) & (z_sub < z_hi)
    n_bin = np.sum(mask_bin)
    if n_bin < 5:
        continue
    
    cc = np.mean(color_corr[mask_bin])
    sc = np.mean(stretch_corr[mask_bin])
    ratio_cs = cc / sc if abs(sc) > 0.001 else float('inf')
    
    print(f"  [{z_lo:.2f}, {z_hi:.2f}]  {n_bin:>5}  {cc:>+14.4f}  {sc:>+16.4f}  {ratio_cs:>18.2f}")

# Overall correlation of corrections with z
rho_color_z, p_color_z = spearmanr(z_sub, color_corr)
rho_stretch_z, p_stretch_z = spearmanr(z_sub, stretch_corr)

print(f"\n  Color correction vs z:   ρ = {rho_color_z:+.3f} (p = {p_color_z:.4f})")
print(f"  Stretch correction vs z: ρ = {rho_stretch_z:+.3f} (p = {p_stretch_z:.4f})")

print(f"\n  PREDICTION: Color correction should evolve with z (impedance channel)")
print(f"  PREDICTION: Stretch correction should be flat (locked channel)")
if abs(rho_color_z) > abs(rho_stretch_z):
    print(f"  RESULT: Color evolves MORE with z than stretch ✓")
else:
    print(f"  RESULT: Stretch evolves MORE — unexpected")


# ============================================================
# TEST 5: SCATTER COMPARISON
# ============================================================

print(f"\n\n{'=' * 70}")
print("TEST 5: SCATTER — TRIPP vs STRETCH-ONLY")
print("=" * 70)

# Compute Hubble residuals for both standardizations
mu_model = distance_modulus(z_sub, H0=70, Om=0.3)

resid_tripp = mu_tripp - mu_model
resid_stretch = mu_stretch - mu_model

# Remove mean (absorbed into M_B)
resid_tripp -= np.mean(resid_tripp)
resid_stretch -= np.mean(resid_stretch)

scatter_tripp = np.std(resid_tripp)
scatter_stretch = np.std(resid_stretch)

print(f"\n  Hubble residual scatter:")
print(f"    Tripp (stretch+color): σ = {scatter_tripp:.4f} mag")
print(f"    Stretch-only:          σ = {scatter_stretch:.4f} mag")
print(f"    Increase: {(scatter_stretch/scatter_tripp - 1)*100:.1f}%")
print()
print(f"  Color correction reduces scatter (useful for precision)")
print(f"  But may introduce BIAS (impedance channel)")
print(f"  Trade-off: less scatter but more systematic bias")

# z-dependent scatter
print(f"\n  Scatter by z-bin:")
print(f"  {'z-bin':>15} {'σ_Tripp':>10} {'σ_Stretch':>12} {'Δσ%':>8}")
print(f"  {'-'*48}")

for z_lo, z_hi in z_bins:
    mask_bin = (z_sub >= z_lo) & (z_sub < z_hi)
    if np.sum(mask_bin) < 10:
        continue
    
    s_t = np.std(resid_tripp[mask_bin])
    s_s = np.std(resid_stretch[mask_bin])
    ds = (s_s / s_t - 1) * 100
    
    print(f"  [{z_lo:.2f}, {z_hi:.2f}]  {s_t:>10.4f}  {s_s:>12.4f}  {ds:>+8.1f}%")


# ============================================================
# SYNTHESIS
# ============================================================

print(f"\n\n{'=' * 70}")
print("SYNTHESIS: STRETCH-ONLY TEST RESULTS")
print("=" * 70)

print(f"""
  FULL SAMPLE:
    Dropping color changes M_B by {delta_MB:+.4f}
    This corresponds to H₀ shifting by factor {H0_ratio:.4f}
    If Tripp gives H₀ = 73.0 → Stretch-only gives H₀ ≈ {73.0*H0_ratio:.1f}
    Direction: {'CORRECT (lower)' if H0_ratio < 1 else 'WRONG (higher)'}
    
  INTERPRETATION:
    The color channel (β×c) in SN Ia standardization carries
    {'an impedance bias that inflates H₀' if H0_ratio < 1 else 'information that is NOT impedance-biased'}.
    
  {'This is consistent with the metric impedance prediction.' if H0_ratio < 1 else 'This challenges the simple metric impedance model.'}
  {'Removing color moves H₀ toward the CMB value.' if H0_ratio < 1 else 'Color may not be the primary impedance channel.'}
""")

# Save results
import json
output_dir = 'results_stretch_only'
os.makedirs(output_dir, exist_ok=True) if not os.path.exists('results_stretch_only') else None

output = {
    'test_date': '2026-03-09',
    'n_sne_total': int(np.sum(mask_quality)),
    'n_sne_subsample': n_sub,
    'tripp_params': {
        'M_B': float(M_B_tripp),
        'alpha': float(alpha_tripp),
        'beta': float(beta_tripp),
    },
    'stretch_only_params': {
        'M_B': float(M_B_stretch),
        'alpha': float(alpha_stretch),
    },
    'delta_MB': float(delta_MB),
    'H0_ratio': float(H0_ratio),
    'H0_tripp_73': float(73.0 * H0_ratio),
    'scatter_tripp': float(scatter_tripp),
    'scatter_stretch': float(scatter_stretch),
    'color_correction_vs_z_rho': float(rho_color_z),
    'stretch_correction_vs_z_rho': float(rho_stretch_z),
}

import os
os.makedirs('results_stretch_only', exist_ok=True)
with open('results_stretch_only/stretch_only_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"  Results saved to results_stretch_only/stretch_only_results.json")
print(f"\n{'=' * 70}")
print("TEST COMPLETE")
print("=" * 70)
