#!/usr/bin/env python3
"""
closure_eddington_control.py — GPT's Eddington Ratio Test
==========================================================

GPT's prediction: MgII EW divergence is driven by Eddington-ratio 
separation between NL and BL bins. If we residualize against L/L_Edd
(not just L_bol), the MgII divergence B_EW(z) should COLLAPSE.

CIV convergence should be UNAFFECTED (wind-driven, not Eddington).

This is the mechanism discriminator:
  - If MgII collapses → accretion-state separation, source-side
  - If MgII persists → something deeper (propagation/transfer)

Uses LOGMBH_MGII and LOGMBH_CIV for virial mass estimates,
then L/L_Edd = L_bol / L_Edd where L_Edd = 1.26e38 * M_BH/M_sun.

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr
from astropy.io import fits
from numpy.polynomial import polynomial as P
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

print("Loading DR16Q catalog...")
f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z_all = d['Z_DR16Q']
lbol = d['LOGLBOL']

# ============================================================
# HELPER: Compute B_EW(z) with configurable residualization
# ============================================================
def compute_B_trend(col_name, z_min, z_max, label, resid_vars, resid_labels):
    """
    Compute B_EW(z) after residualizing log(EW) against specified variables.
    Returns dict with ρ_B, p_B, and per-bin B values.
    """
    ld = d[col_name]
    
    mask = ((z_all >= z_min) & (z_all < z_max) & 
            (ld[:, 4] > 100) & (ld[:, 4] < 30000) &
            (ld[:, 2] > 0) & (ld[:, 2] < 5000) &
            np.isfinite(ld[:, 2]) & np.isfinite(ld[:, 4]) &
            (lbol > 40) & (lbol < 50) & np.isfinite(lbol))
    
    # Additional masking for residualization variables
    for rv in resid_vars:
        mask &= np.isfinite(rv) & (rv > -90) & (rv < 90)
    
    N = np.sum(mask)
    if N < 2000:
        return None
    
    z_v = z_all[mask]
    ew_v = ld[mask, 2]
    fwhm_v = ld[mask, 4]
    
    # Residualize log(EW)
    log_ew = np.log10(np.maximum(ew_v, 0.1))
    
    for rv in resid_vars:
        rv_masked = rv[mask]
        coef = P.polyfit(rv_masked, log_ew, 1)
        log_ew = log_ew - P.polyval(rv_masked, coef)
    
    # Median FWHM split
    fwhm_med = np.median(fwhm_v)
    
    # Z-bins
    n_bins = min(8, N // 3000)
    if n_bins < 3:
        n_bins = 3
    z_edges = np.percentile(z_v, np.linspace(0, 100, n_bins + 1))
    z_edges = np.unique(np.round(z_edges, 3))
    
    z_list, B_list, delta_list = [], [], []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_mask = (z_v >= z_lo) & (z_v < z_hi)
        nl = bin_mask & (fwhm_v < fwhm_med)
        bl = bin_mask & (fwhm_v >= fwhm_med)
        
        if np.sum(nl) < 50 or np.sum(bl) < 50:
            continue
        
        zc = np.mean(z_v[bin_mask])
        all_vals = log_ew[bin_mask]
        nl_vals = log_ew[nl]
        bl_vals = log_ew[bl]
        
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
        delta_list.append(mu_bl - mu_nl)
    
    if len(z_list) < 3:
        return None
    
    z_arr = np.array(z_list)
    rho_B, p_B = spearmanr(z_arr, B_list)
    rho_d, p_d = spearmanr(z_arr, delta_list)
    
    return {
        'label': label,
        'resid': '+'.join(resid_labels),
        'N': N,
        'rho_B': rho_B,
        'p_B': p_B,
        'rho_delta': rho_d,
        'p_delta': p_d,
        'B_first': B_list[0],
        'B_last': B_list[-1],
        'z_bins': z_list,
        'B_bins': B_list,
    }


# ============================================================
# MgII TEST
# ============================================================
print(f"\n{'=' * 100}")
print("MgII_BR: Does Eddington ratio control collapse the divergence?")
print("=" * 100)

# Eddington ratio from MgII-based BH mass
logmbh_mgii = d['LOGMBH_MGII']
logledd_mgii = lbol - logmbh_mgii - np.log10(1.26e38) + np.log10(3.846e33)  # L_Edd in erg/s
# Actually: LOGLEDD_RATIO is already in the catalog
logledd = d['LOGLEDD_RATIO']

print("\n  Control 1: Residualize against L_bol only (baseline)")
r1 = compute_B_trend('MGII_BR', 0.4, 2.5, 'MgII', [lbol], ['L_bol'])
if r1:
    print(f"    N={r1['N']}, ρ_B={r1['rho_B']:+.3f}, p={r1['p_B']:.4f}, B: {r1['B_first']:.4f}→{r1['B_last']:.4f}")

print("\n  Control 2: Residualize against L_bol + L/L_Edd (GPT's test)")
r2 = compute_B_trend('MGII_BR', 0.4, 2.5, 'MgII', [lbol, logledd], ['L_bol', 'L/L_Edd'])
if r2:
    print(f"    N={r2['N']}, ρ_B={r2['rho_B']:+.3f}, p={r2['p_B']:.4f}, B: {r2['B_first']:.4f}→{r2['B_last']:.4f}")

print("\n  Control 3: Residualize against L_bol + LOGMBH (mass directly)")
r3 = compute_B_trend('MGII_BR', 0.4, 2.5, 'MgII', [lbol, logmbh_mgii], ['L_bol', 'M_BH'])
if r3:
    print(f"    N={r3['N']}, ρ_B={r3['rho_B']:+.3f}, p={r3['p_B']:.4f}, B: {r3['B_first']:.4f}→{r3['B_last']:.4f}")

# ============================================================
# CIV TEST 
# ============================================================
print(f"\n{'=' * 100}")
print("CIV: Does Eddington ratio control affect the convergence?")
print("(GPT predicts: should be UNAFFECTED — wind-driven)")
print("=" * 100)

logmbh_civ = d['LOGMBH_CIV']

print("\n  Control 1: Residualize against L_bol only (baseline)")
c1 = compute_B_trend('CIV', 1.5, 4.0, 'CIV', [lbol], ['L_bol'])
if c1:
    print(f"    N={c1['N']}, ρ_B={c1['rho_B']:+.3f}, p={c1['p_B']:.4f}, B: {c1['B_first']:.4f}→{c1['B_last']:.4f}")

print("\n  Control 2: Residualize against L_bol + L/L_Edd")
c2 = compute_B_trend('CIV', 1.5, 4.0, 'CIV', [lbol, logledd], ['L_bol', 'L/L_Edd'])
if c2:
    print(f"    N={c2['N']}, ρ_B={c2['rho_B']:+.3f}, p={c2['p_B']:.4f}, B: {c2['B_first']:.4f}→{c2['B_last']:.4f}")

print("\n  Control 3: Residualize against L_bol + LOGMBH")
c3 = compute_B_trend('CIV', 1.5, 4.0, 'CIV', [lbol, logmbh_civ], ['L_bol', 'M_BH'])
if c3:
    print(f"    N={c3['N']}, ρ_B={c3['rho_B']:+.3f}, p={c3['p_B']:.4f}, B: {c3['B_first']:.4f}→{c3['B_last']:.4f}")

# ============================================================
# VERDICT
# ============================================================
print(f"\n{'=' * 100}")
print("VERDICT")
print("=" * 100)

if r1 and r2:
    mgii_collapse = abs(r2['rho_B']) < abs(r1['rho_B']) * 0.5
    print(f"\n  MgII divergence with L_bol only:       ρ_B = {r1['rho_B']:+.3f}")
    print(f"  MgII divergence with L_bol + L/L_Edd:  ρ_B = {r2['rho_B']:+.3f}")
    if r3:
        print(f"  MgII divergence with L_bol + M_BH:     ρ_B = {r3['rho_B']:+.3f}")
    
    if mgii_collapse:
        print(f"\n  ✓ GPT CONFIRMED: MgII divergence COLLAPSES with Eddington control")
        print(f"    → Source-side accretion-state separation is the driver")
    else:
        print(f"\n  ✗ GPT prediction FAILS: MgII divergence PERSISTS after Eddington control")
        print(f"    → Something beyond accretion state drives the divergence")

if c1 and c2:
    civ_stable = abs(c2['rho_B'] - c1['rho_B']) < 0.3
    print(f"\n  CIV convergence with L_bol only:       ρ_B = {c1['rho_B']:+.3f}")
    print(f"  CIV convergence with L_bol + L/L_Edd:  ρ_B = {c2['rho_B']:+.3f}")
    if c3:
        print(f"  CIV convergence with L_bol + M_BH:     ρ_B = {c3['rho_B']:+.3f}")
    
    if civ_stable:
        print(f"\n  ✓ GPT CONFIRMED: CIV convergence UNAFFECTED by Eddington control")
        print(f"    → Wind/transfer physics drives CIV, not accretion state")
    else:
        print(f"\n  ⚠️  CIV convergence CHANGED by Eddington control")
        print(f"    → More complex than pure wind physics")

# ============================================================
# ALSO: Run the full ladder with Eddington control
# ============================================================
print(f"\n{'=' * 100}")
print("FULL LADDER WITH EDDINGTON CONTROL")
print("Does the ionization gradient survive after L/L_Edd residualization?")
print("=" * 100)

LINE_DEFS = [
    ('HALPHA_BR',  0.15, 0.55,  13.6,  'Hα_BR'),
    ('HBETA_BR',   0.30, 1.00,  13.6,  'Hβ_BR'),
    ('MGII_BR',    0.40, 2.50,  15.0,  'MgII_BR'),
    ('CIII_BR',    0.80, 3.00,  47.9,  'CIII_BR'),
    ('NIII1750',   1.50, 3.50,  47.4,  'NIII1750'),
    ('SIIV_OIV',   1.50, 4.00,  45.1,  'SiIV+OIV'),
    ('CIV',        1.50, 4.00,  64.5,  'CIV'),
]

results_lbol = []
results_edd = []

for col, z_min, z_max, ip, label in LINE_DEFS:
    r_lbol = compute_B_trend(col, z_min, z_max, label, [lbol], ['L_bol'])
    r_edd = compute_B_trend(col, z_min, z_max, label, [lbol, logledd], ['L_bol', 'L/L_Edd'])
    if r_lbol:
        r_lbol['ip'] = ip
        results_lbol.append(r_lbol)
    if r_edd:
        r_edd['ip'] = ip
        results_edd.append(r_edd)

print(f"\n  {'Line':>12} {'IP':>6} | {'ρ_B (L_bol)':>12} {'ρ_B (+L/LEdd)':>14} | {'Change':>8}")
print(f"  " + "-" * 65)
for r1, r2 in zip(sorted(results_lbol, key=lambda x: x['ip']), 
                   sorted(results_edd, key=lambda x: x['ip'])):
    change = r2['rho_B'] - r1['rho_B']
    print(f"  {r1['label']:>12} {r1['ip']:>6.1f} | {r1['rho_B']:>+12.3f} {r2['rho_B']:>+14.3f} | {change:>+8.3f}")

if len(results_edd) >= 3:
    ips_edd = [r['ip'] for r in sorted(results_edd, key=lambda x: x['ip'])]
    rhos_edd = [r['rho_B'] for r in sorted(results_edd, key=lambda x: x['ip'])]
    rho_meta, p_meta = spearmanr(ips_edd, rhos_edd)
    print(f"\n  Ladder ρ(IP, ρ_B) after Eddington control: {rho_meta:+.3f}, p = {p_meta:.4f}")
    
    ips_lbol = [r['ip'] for r in sorted(results_lbol, key=lambda x: x['ip'])]
    rhos_lbol = [r['rho_B'] for r in sorted(results_lbol, key=lambda x: x['ip'])]
    rho_meta_orig, _ = spearmanr(ips_lbol, rhos_lbol)
    print(f"  Ladder ρ(IP, ρ_B) L_bol only:              {rho_meta_orig:+.3f}")
    
    if abs(rho_meta) < abs(rho_meta_orig) * 0.5:
        print(f"\n  ⚠️  Eddington control WEAKENS the ionization ladder")
        print(f"    → Accretion state drives part of the IP gradient")
    else:
        print(f"\n  ✓ Ionization ladder SURVIVES Eddington control")
        print(f"    → The IP gradient is NOT just an Eddington-ratio artifact")
