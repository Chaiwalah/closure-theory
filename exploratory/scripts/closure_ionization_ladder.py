#!/usr/bin/env python3
"""
closure_ionization_ladder.py — Ionization Ladder Test
======================================================

Grok's prediction: the sign of B(z) for EW is determined by ionization
potential / BLR formation radius.

  HIGH ionization (inner BLR, short conditioning time) → CONVERGE (like CIV)
  LOW ionization (outer BLR, long conditioning time) → DIVERGE (like MgII)

Ionization ladder (ionization potential to create the ion):
  Line          Ion     IP (eV)    BLR zone    Prediction
  -------------------------------------------------------
  Hα_BR         H I     13.6       outer       DIVERGE
  Hβ_BR         H I     13.6       outer       DIVERGE  
  MgII_BR       Mg II   15.0       outer       DIVERGE  (confirmed)
  CIII_BR       C III   47.9       intermediate MIXED?
  SiIV+OIV      Si IV   45.1       inner       CONVERGE
  CIV           C IV    64.5       inner       CONVERGE (confirmed)
  NIII1750      N III   47.4       intermediate MIXED?
  Lyα           H I     13.6       inner(!)    SPECIAL (resonance)

Test: compute B_EW(z) slope for each line and plot against IP.
If ρ(IP, slope_B) is negative → higher IP = more convergent → prediction confirmed.

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr
from astropy.io import fits
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# Load
print("Loading DR16Q catalog...")
f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z_all = d['Z_DR16Q']
lbol = d['LOGLBOL']

# Line definitions: (column_name, z_min, z_max, ionization_potential_eV, label)
LINE_DEFS = [
    ('HALPHA_BR',  0.15, 0.55,  13.6,  'Hα_BR'),
    ('HBETA_BR',   0.30, 1.00,  13.6,  'Hβ_BR'),
    ('MGII_BR',    0.40, 2.50,  15.0,  'MgII_BR'),
    ('CIII_BR',    0.80, 3.00,  47.9,  'CIII_BR'),
    ('NIII1750',   1.50, 3.50,  47.4,  'NIII1750'),
    ('SIIV_OIV',   1.50, 4.00,  45.1,  'SiIV+OIV'),
    ('CIV',        1.50, 4.00,  64.5,  'CIV'),
]

def compute_B_EW_trend(col_name, z_min, z_max, label, ip):
    """Compute B_EW(z) trend for a single line."""
    ld = d[col_name]
    
    # Quality cuts
    mask = ((z_all >= z_min) & (z_all < z_max) & 
            (ld[:, 4] > 100) & (ld[:, 4] < 30000) &
            (ld[:, 2] > 0) & (ld[:, 2] < 5000) &
            np.isfinite(ld[:, 2]) & np.isfinite(ld[:, 4]) &
            (lbol > 40) & (lbol < 50) & np.isfinite(lbol))
    
    N_total = np.sum(mask)
    if N_total < 2000:
        return None
    
    z_valid = z_all[mask]
    ew_valid = ld[mask, 2]
    fwhm_valid = ld[mask, 4]
    lbol_valid = lbol[mask]
    
    # Residualize EW against luminosity
    log_ew = np.log10(np.maximum(ew_valid, 0.1))
    from numpy.polynomial import polynomial as P
    coef = P.polyfit(lbol_valid, log_ew, 1)
    log_ew_resid = log_ew - P.polyval(lbol_valid, coef)
    
    # Median FWHM split
    fwhm_med = np.median(fwhm_valid)
    
    # Adaptive z-bins
    n_bins = min(8, N_total // 3000)
    if n_bins < 3:
        n_bins = 3
    z_edges = np.percentile(z_valid, np.linspace(0, 100, n_bins + 1))
    z_edges = np.unique(np.round(z_edges, 3))
    
    z_list, B_list, delta_ew_list = [], [], []
    iqr_ew_list = []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_mask = (z_valid >= z_lo) & (z_valid < z_hi)
        
        nl = bin_mask & (fwhm_valid < fwhm_med)
        bl = bin_mask & (fwhm_valid >= fwhm_med)
        
        if np.sum(nl) < 50 or np.sum(bl) < 50:
            continue
        
        zc = np.mean(z_valid[bin_mask])
        
        # B(z) for residualized EW
        all_vals = log_ew_resid[bin_mask]
        nl_vals = log_ew_resid[nl]
        bl_vals = log_ew_resid[bl]
        
        var_total = np.var(all_vals)
        if var_total < 1e-10:
            continue
        
        n_nl, n_bl = len(nl_vals), len(bl_vals)
        n_tot = n_nl + n_bl
        mu_all = np.mean(all_vals)
        mu_nl, mu_bl = np.mean(nl_vals), np.mean(bl_vals)
        
        var_between = (n_nl/n_tot)*(mu_nl-mu_all)**2 + (n_bl/n_tot)*(mu_bl-mu_all)**2
        B = var_between / var_total
        
        z_list.append(zc)
        B_list.append(B)
        delta_ew_list.append(mu_bl - mu_nl)
        iqr_ew_list.append(np.percentile(all_vals, 90) - np.percentile(all_vals, 10))
    
    if len(z_list) < 3:
        return None
    
    z_arr = np.array(z_list)
    rho_B, p_B = spearmanr(z_arr, B_list)
    rho_delta, p_delta = spearmanr(z_arr, delta_ew_list)
    rho_iqr, p_iqr = spearmanr(z_arr, iqr_ew_list)
    
    return {
        'label': label,
        'ip': ip,
        'N': N_total,
        'n_bins': len(z_list),
        'z_range': (z_list[0], z_list[-1]),
        'rho_B': rho_B,
        'p_B': p_B,
        'rho_delta': rho_delta,
        'p_delta': p_delta,
        'rho_iqr': rho_iqr,
        'p_iqr': p_iqr,
        'B_first': B_list[0],
        'B_last': B_list[-1],
    }

# ============================================================
# RUN THE LADDER
# ============================================================
print(f"\n{'=' * 120}")
print("IONIZATION LADDER TEST")
print("Prediction: B_EW(z) slope correlates with ionization potential")
print("  High IP (inner BLR) → B drops (converge) → negative ρ_B")
print("  Low IP (outer BLR) → B rises (diverge) → positive ρ_B")
print(f"{'=' * 120}\n")

results = []
for col, z_min, z_max, ip, label in LINE_DEFS:
    r = compute_B_EW_trend(col, z_min, z_max, label, ip)
    if r:
        results.append(r)

# Print results table
print(f"{'Line':>12} {'IP(eV)':>8} {'N':>8} {'bins':>5} | {'ρ(B,z)':>8} {'p':>8} {'B_lo-z':>8} {'B_hi-z':>8} | {'ρ(ΔEW,z)':>10} {'p':>8} | {'ρ(IQR,z)':>10} {'p':>8} | {'Verdict':>12}")
print("-" * 130)

for r in sorted(results, key=lambda x: x['ip']):
    verdict = "DIVERGE" if r['rho_B'] > 0.3 else "CONVERGE" if r['rho_B'] < -0.3 else "FLAT"
    sig = "***" if r['p_B'] < 0.001 else "**" if r['p_B'] < 0.01 else "*" if r['p_B'] < 0.05 else ""
    
    print(f"{r['label']:>12} {r['ip']:>8.1f} {r['N']:>8} {r['n_bins']:>5} | "
          f"{r['rho_B']:>+8.3f} {r['p_B']:>8.4f} {r['B_first']:>8.4f} {r['B_last']:>8.4f} | "
          f"{r['rho_delta']:>+10.3f} {r['p_delta']:>8.4f} | "
          f"{r['rho_iqr']:>+10.3f} {r['p_iqr']:>8.4f} | "
          f"{verdict:>10}{sig}")

# ============================================================
# THE KEY TEST: ρ(IP, slope_B) across lines
# ============================================================
print(f"\n{'=' * 120}")
print("KEY TEST: Does ionization potential predict B(z) direction?")
print(f"{'=' * 120}")

ips = [r['ip'] for r in results]
rho_Bs = [r['rho_B'] for r in results]
delta_rhos = [r['rho_delta'] for r in results]
iqr_rhos = [r['rho_iqr'] for r in results]

if len(ips) >= 3:
    rho_meta, p_meta = spearmanr(ips, rho_Bs)
    rho_delta_meta, p_delta_meta = spearmanr(ips, delta_rhos)
    rho_iqr_meta, p_iqr_meta = spearmanr(ips, iqr_rhos)
    
    print(f"\n  ρ(IP, ρ_B):     {rho_meta:+.3f}, p = {p_meta:.4f}")
    print(f"  ρ(IP, ρ_ΔEW):   {rho_delta_meta:+.3f}, p = {p_delta_meta:.4f}")
    print(f"  ρ(IP, ρ_IQR):   {rho_iqr_meta:+.3f}, p = {p_iqr_meta:.4f}")
    
    if rho_meta < -0.5 and p_meta < 0.1:
        print(f"\n  🔥 IONIZATION LADDER CONFIRMED")
        print(f"     Higher IP → more convergent B(z)")
        print(f"     The BLR formation radius sets the 'delay time'")
        print(f"     One phenomenon, ionization-stratified windows")
    elif rho_meta > 0.5:
        print(f"\n  ✗ IONIZATION LADDER REVERSED — prediction fails")
    else:
        print(f"\n  ⚠️  IONIZATION LADDER INCONCLUSIVE (|ρ| < 0.5)")
    
    # Print the ladder
    print(f"\n  Ionization Ladder (sorted by IP):")
    for r in sorted(results, key=lambda x: x['ip']):
        arrow = "↗ DIVERGE" if r['rho_B'] > 0.3 else "↘ CONVERGE" if r['rho_B'] < -0.3 else "→ FLAT"
        print(f"    {r['ip']:>5.1f} eV  {r['label']:>12}  ρ_B = {r['rho_B']:+.3f}  {arrow}")

# ============================================================
# BONUS: IQR universality check
# ============================================================
print(f"\n{'=' * 120}")
print("BONUS: Does IQR shrink universally across ALL lines?")
print(f"{'=' * 120}")

n_shrink = sum(1 for r in results if r['rho_iqr'] < -0.3)
n_total = len(results)
print(f"\n  Lines with IQR shrinkage (ρ < −0.3): {n_shrink}/{n_total}")
for r in sorted(results, key=lambda x: x['rho_iqr']):
    sig = "🔥" if r['p_iqr'] < 0.05 else ""
    print(f"    {r['label']:>12}: ρ(IQR,z) = {r['rho_iqr']:+.3f}, p = {r['p_iqr']:.4f} {sig}")
