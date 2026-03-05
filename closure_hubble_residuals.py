#!/usr/bin/env python3
"""
closure_hubble_residuals.py — Hubble Residual Scatter vs z
============================================================

DISCRIMINANT TEST: Is the lost standardization power (16.2%/z) absorbed
into the geometric distance (coupled/conserved) or genuinely lost (dissipated)?

If COUPLED: Hubble residual scatter stays constant with z
  (lost diagnostic info → geometric channel → distances still accurate)
If DISSIPATED: Hubble residual scatter INCREASES with z
  (lost info → nowhere → distances get noisier)

Also: Does the residual trend correlate with the β(z) degradation?

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr
from scipy.optimize import minimize
import json
from pathlib import Path

np.random.seed(42)

# Cosmology
H0 = 70.0  # SH0ES-like for Pantheon+ corrected magnitudes
Omega_m = 0.315
c_km = 299792.458

def luminosity_distance(z, Om=Omega_m):
    """Simple numerical integration for d_L"""
    from scipy.integrate import quad
    def integrand(zz):
        return 1.0 / np.sqrt(Om * (1+zz)**3 + (1-Om))
    result, _ = quad(integrand, 0, z)
    return (1+z) * c_km / H0 * result  # in Mpc

def mu_theory(z, Om=Omega_m):
    """Theoretical distance modulus"""
    dL = luminosity_distance(z, Om)
    return 5 * np.log10(dL) + 25

# Load data
data = []
with open('data/pantheon_plus.dat', 'r') as f:
    header = f.readline().strip().split()
    for line in f:
        parts = line.strip().split()
        if len(parts) >= len(header):
            row = {}
            for i, h in enumerate(header):
                try:
                    row[h] = float(parts[i])
                except:
                    row[h] = parts[i]
            data.append(row)

seen = set()
sne = []
for row in data:
    z = row.get('zHD', 0)
    c = row.get('c', -999)
    x1 = row.get('x1', -999)
    if z > 0.01 and z < 2.5 and abs(c) < 0.3 and abs(x1) < 3:
        key = (row['CID'], round(z, 5))
        if key not in seen:
            seen.add(key)
            sne.append(row)

print(f"Loaded {len(sne)} unique SNe Ia")

# ============================================================
# Compute Hubble residuals using corrected magnitudes
# ============================================================
# m_b_corr = mB + α*x1 - β*c - M (already corrected in Pantheon+)
# Use m_b_corr directly as the standardized apparent magnitude

# Find best-fit M offset (nuisance parameter)
zs = np.array([s['zHD'] for s in sne])
m_corr = np.array([s['m_b_corr'] for s in sne])
mu_th = np.array([mu_theory(z) for z in zs])

# M = median(m_corr - mu_theory)
M_offset = np.median(m_corr - mu_th)
residuals = m_corr - mu_th - M_offset

print(f"M offset: {M_offset:.4f}")
print(f"Residual RMS: {np.std(residuals):.4f} mag")

# ============================================================
# TEST: Residual scatter vs z
# ============================================================
print("\n" + "=" * 70)
print("HUBBLE RESIDUAL SCATTER VS REDSHIFT")
print("=" * 70)

z_edges = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.55, 0.7, 0.85, 1.0, 1.5, 2.5]

z_centers = []
scatters = []
mean_resids = []
n_bins = []
mad_scatters = []

print(f"\n{'z_range':<14} {'N':>5} {'⟨Δμ⟩':>8} {'σ(Δμ)':>8} {'MAD':>8}")
print("-" * 50)

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    mask = (zs >= z_lo) & (zs < z_hi)
    bin_resid = residuals[mask]
    
    if len(bin_resid) < 5:
        continue
    
    z_center = np.mean(zs[mask])
    scatter = np.std(bin_resid)
    mean_r = np.mean(bin_resid)
    mad = np.median(np.abs(bin_resid - np.median(bin_resid))) * 1.4826  # normalized MAD
    
    z_centers.append(z_center)
    scatters.append(scatter)
    mean_resids.append(mean_r)
    n_bins.append(np.sum(mask))
    mad_scatters.append(mad)
    
    print(f"[{z_lo:.2f},{z_hi:.2f}){'':<5} {np.sum(mask):>5} {mean_r:>8.4f} {scatter:>8.4f} {mad:>8.4f}")

z_arr = np.array(z_centers)
scatter_arr = np.array(scatters)
mad_arr = np.array(mad_scatters)

rho_scatter, p_scatter = spearmanr(z_arr, scatter_arr)
rho_mad, p_mad = spearmanr(z_arr, mad_arr)

print(f"\n  σ(Δμ) vs z: ρ = {rho_scatter:+.3f}, p = {p_scatter:.4f}")
print(f"  MAD vs z:   ρ = {rho_mad:+.3f}, p = {p_mad:.4f}")

if rho_scatter > 0.3 and p_scatter < 0.1:
    print(f"\n  ❄️ Hubble residual scatter INCREASES with z")
    print(f"     → Information is DISSIPATED (not conserved)")
    print(f"     → Lost standardization power → noisier distances")
    print(f"     → BISTABLE SWITCH DISFAVORED")
elif abs(rho_scatter) < 0.3:
    print(f"\n  🔥 Hubble residual scatter is FLAT with z")
    print(f"     → Lost diagnostic info is ABSORBED into geometry")
    print(f"     → Consistent with coupled/conserved system")
    print(f"     → BISTABLE SWITCH STILL VIABLE")
elif rho_scatter < -0.3:
    print(f"\n  ★ Hubble residual scatter DECREASES with z")
    print(f"     → Distances get BETTER despite losing diagnostic power")
    print(f"     → Strong evidence for geometric absorption of lost info")

# ============================================================
# TEST: Mean residual trend (systematic bias)
# ============================================================
print("\n" + "=" * 70)
print("MEAN RESIDUAL TREND (systematic distance bias)")
print("=" * 70)

rho_mean, p_mean = spearmanr(z_arr, mean_resids)
print(f"  ⟨Δμ⟩ vs z: ρ = {rho_mean:+.3f}, p = {p_mean:.4f}")

if abs(rho_mean) > 0.3:
    # Fit linear trend
    coeffs = np.polyfit(z_arr, mean_resids, 1)
    print(f"  Linear fit: ⟨Δμ⟩ = {coeffs[1]:.4f} + {coeffs[0]:+.4f} × z")
    print(f"  At z=0: {coeffs[1]:.4f} mag")
    print(f"  At z=1: {coeffs[1]+coeffs[0]:.4f} mag")
    
    if coeffs[0] > 0:
        print(f"\n  Objects appear FAINTER than expected at high z")
        print(f"  This IS the dark energy signal — but under compression:")
        print(f"  It could be β degradation making high-z SNe look fainter")
    elif coeffs[0] < 0:
        print(f"\n  Objects appear BRIGHTER than expected at high z")

# ============================================================
# TEST: Residual scatter vs β(z) — direct link
# ============================================================
print("\n" + "=" * 70)
print("RESIDUAL SCATTER vs LOCAL β — Direct correlation")
print("=" * 70)

# For each z bin, compute both β and residual scatter
# If they're correlated, it means β degradation directly causes distance scatter

# Recompute β per bin
betas = []
z_beta_centers = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    
    if len(bin_sne) < 10:
        continue
    
    colors = np.array([s['c'] for s in bin_sne])
    stretches = np.array([s['x1'] for s in bin_sne])
    mBs = np.array([s['mB'] for s in bin_sne])
    
    A = np.column_stack([np.ones(len(bin_sne)), colors, stretches])
    try:
        params = np.linalg.lstsq(A, mBs, rcond=None)[0]
        betas.append(params[1])
        z_beta_centers.append(np.mean([s['zHD'] for s in bin_sne]))
    except:
        pass

# Match z bins between β and scatter measurements
if len(betas) >= 4 and len(scatters) >= 4:
    # Interpolate to common z grid
    common_z = sorted(set(z_centers) & set(z_beta_centers))
    if len(common_z) < 4:
        # Use nearest matching
        common_betas = []
        common_scatters = []
        for zc in z_centers:
            # Find closest beta
            dists = [abs(zc - zb) for zb in z_beta_centers]
            if min(dists) < 0.05:
                idx = np.argmin(dists)
                common_betas.append(betas[idx])
                common_scatters.append(scatters[z_centers.index(zc)])
        
        if len(common_betas) >= 4:
            rho_bs, p_bs = spearmanr(common_betas, common_scatters)
            print(f"  β vs σ(Δμ): ρ = {rho_bs:+.3f}, p = {p_bs:.4f}")
    else:
        beta_arr = np.array(betas)
        rho_bs, p_bs = spearmanr(beta_arr[:len(scatter_arr)], scatter_arr[:len(beta_arr)])
        print(f"  β vs σ(Δμ): ρ = {rho_bs:+.3f}, p = {p_bs:.4f}")
    
        if rho_bs > 0.3:
            print(f"  → As β drops, scatter drops too (counter-intuitive)")
            print(f"     Losing diagnostic power doesn't hurt distances")
        elif rho_bs < -0.3:
            print(f"  → As β drops, scatter rises (expected if dissipated)")

# ============================================================
# TEST: Raw vs Corrected magnitude scatter
# ============================================================
print("\n" + "=" * 70)
print("RAW vs CORRECTED MAGNITUDE SCATTER")
print("=" * 70)

print(f"\n{'z_range':<14} {'N':>5} {'σ(mB)_raw':>10} {'σ(mB)_corr':>11} {'ratio':>8}")
print("-" * 55)

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    
    if len(bin_sne) < 10:
        continue
    
    raw_mB = np.array([s['mB'] for s in bin_sne])
    corr_mB = np.array([s['m_b_corr'] for s in bin_sne])
    raw_mu = np.array([mu_theory(s['zHD']) for s in bin_sne])
    
    raw_resid = raw_mB - raw_mu
    corr_resid = corr_mB - raw_mu
    
    sigma_raw = np.std(raw_resid)
    sigma_corr = np.std(corr_resid)
    ratio = sigma_corr / sigma_raw if sigma_raw > 0 else 999
    
    print(f"[{z_lo:.2f},{z_hi:.2f}){'':<5} {len(bin_sne):>5} {sigma_raw:>10.4f} {sigma_corr:>11.4f} {ratio:>8.3f}")

print(f"\n  ratio = σ(corrected)/σ(raw) = standardization improvement factor")
print(f"  If ratio INCREASES with z → standardization becomes LESS effective")
print(f"  If ratio stays constant → standardization absorbs the lost info")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_hubble_residuals")
results_dir.mkdir(exist_ok=True)

results = {
    'scatter_vs_z_rho': float(rho_scatter),
    'scatter_vs_z_p': float(p_scatter),
    'mad_vs_z_rho': float(rho_mad),
    'mad_vs_z_p': float(p_mad),
    'mean_resid_vs_z_rho': float(rho_mean),
    'z_centers': [float(z) for z in z_arr],
    'scatters': [float(s) for s in scatter_arr],
    'mean_resids': [float(m) for m in mean_resids],
}

with open(results_dir / "hubble_residuals.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {results_dir}/")
print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)
if rho_scatter > 0.3 and p_scatter < 0.1:
    print("  DISSIPATED: Lost info → nowhere → distances noisier → bistable DISFAVORED")
elif abs(rho_scatter) < 0.3:
    print("  INCONCLUSIVE or CONSERVED: Scatter flat → can't rule out coupled system")
else:
    print("  CONSERVED: Scatter improves → geometric absorption → bistable SUPPORTED")
