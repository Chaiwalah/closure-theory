#!/usr/bin/env python3
"""
closure_feii_control.py — FeII / Eigenvector 1 Control Test
=============================================================

GPT's prediction: MgII divergence is driven by BLR covering/structure,
not accretion state. FeII strength is the best proxy for BLR covering
factor and Eigenvector 1 (the primary axis of quasar diversity).

Test: Add FeII_UV_EW to the residualization of MgII EW.
  - If MgII B_EW(z) collapses → covering/structure IS the driver
  - If MgII B_EW(z) persists → something even deeper

Also test:
  - CIV with FeII control (should be less affected)
  - Full ladder with FeII control

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
feii_uv = d['FEII_UV_EW']
logledd = d['LOGLEDD_RATIO']

def compute_B_trend(col_name, z_min, z_max, label, resid_vars, resid_labels, extra_mask=None):
    ld = d[col_name]
    
    mask = ((z_all >= z_min) & (z_all < z_max) & 
            (ld[:, 4] > 100) & (ld[:, 4] < 30000) &
            (ld[:, 2] > 0) & (ld[:, 2] < 5000) &
            np.isfinite(ld[:, 2]) & np.isfinite(ld[:, 4]) &
            (lbol > 40) & (lbol < 50) & np.isfinite(lbol))
    
    for rv in resid_vars:
        mask &= np.isfinite(rv) & (rv > -999) & (rv < 999)
    
    if extra_mask is not None:
        mask &= extra_mask
    
    N = np.sum(mask)
    if N < 2000:
        return None
    
    z_v = z_all[mask]
    ew_v = ld[mask, 2]
    fwhm_v = ld[mask, 4]
    
    log_ew = np.log10(np.maximum(ew_v, 0.1))
    
    for rv in resid_vars:
        rv_m = rv[mask]
        coef = P.polyfit(rv_m, log_ew, 1)
        log_ew = log_ew - P.polyval(rv_m, coef)
    
    fwhm_med = np.median(fwhm_v)
    
    n_bins = min(8, N // 3000)
    if n_bins < 3:
        n_bins = 3
    z_edges = np.percentile(z_v, np.linspace(0, 100, n_bins + 1))
    z_edges = np.unique(np.round(z_edges, 3))
    
    z_list, B_list, delta_list, iqr_list = [], [], [], []
    
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
        iqr_list.append(np.percentile(all_vals, 90) - np.percentile(all_vals, 10))
    
    if len(z_list) < 3:
        return None
    
    rho_B, p_B = spearmanr(z_list, B_list)
    rho_d, p_d = spearmanr(z_list, delta_list)
    rho_iqr, p_iqr = spearmanr(z_list, iqr_list)
    
    return {
        'label': label, 'resid': '+'.join(resid_labels), 'N': N,
        'rho_B': rho_B, 'p_B': p_B,
        'rho_delta': rho_d, 'p_delta': p_d,
        'rho_iqr': rho_iqr, 'p_iqr': p_iqr,
        'B_first': B_list[0], 'B_last': B_list[-1],
    }

# FeII mask (UV)
feii_valid = (feii_uv > 0) & np.isfinite(feii_uv)
log_feii = np.where(feii_valid, np.log10(np.maximum(feii_uv, 0.01)), 0)

# ============================================================
# MgII CONTROLS
# ============================================================
print(f"\n{'=' * 100}")
print("MgII_BR: Progressive Residualization")
print("=" * 100)

controls = [
    ([lbol], ['L_bol'], None),
    ([lbol, logledd], ['L_bol', 'L/L_Edd'], None),
    ([lbol, log_feii], ['L_bol', 'FeII_UV'], feii_valid),
    ([lbol, logledd, log_feii], ['L_bol', 'L/L_Edd', 'FeII_UV'], feii_valid),
]

print(f"\n  {'Control':>35} {'N':>8} | {'ρ_B':>7} {'p':>8} | {'ρ_ΔEW':>7} {'p':>8} | {'ρ_IQR':>7} {'p':>8}")
print(f"  " + "-" * 100)

mgii_results = []
for resid_vars, resid_labels, extra in controls:
    r = compute_B_trend('MGII_BR', 0.4, 2.5, 'MgII', resid_vars, resid_labels, extra)
    if r:
        mgii_results.append(r)
        print(f"  {r['resid']:>35} {r['N']:>8} | {r['rho_B']:>+7.3f} {r['p_B']:>8.4f} | "
              f"{r['rho_delta']:>+7.3f} {r['p_delta']:>8.4f} | {r['rho_iqr']:>+7.3f} {r['p_iqr']:>8.4f}")

# ============================================================
# CIV CONTROLS
# ============================================================
print(f"\n{'=' * 100}")
print("CIV: Progressive Residualization")
print("=" * 100)

print(f"\n  {'Control':>35} {'N':>8} | {'ρ_B':>7} {'p':>8} | {'ρ_ΔEW':>7} {'p':>8} | {'ρ_IQR':>7} {'p':>8}")
print(f"  " + "-" * 100)

civ_results = []
for resid_vars, resid_labels, extra in controls:
    r = compute_B_trend('CIV', 1.5, 4.0, 'CIV', resid_vars, resid_labels, extra)
    if r:
        civ_results.append(r)
        print(f"  {r['resid']:>35} {r['N']:>8} | {r['rho_B']:>+7.3f} {r['p_B']:>8.4f} | "
              f"{r['rho_delta']:>+7.3f} {r['p_delta']:>8.4f} | {r['rho_iqr']:>+7.3f} {r['p_iqr']:>8.4f}")

# ============================================================
# Hβ CONTROLS (should be immune to everything)
# ============================================================
print(f"\n{'=' * 100}")
print("Hβ_BR: Progressive Residualization (control line)")
print("=" * 100)

print(f"\n  {'Control':>35} {'N':>8} | {'ρ_B':>7} {'p':>8} | {'ρ_ΔEW':>7} {'p':>8}")
print(f"  " + "-" * 100)

for resid_vars, resid_labels, extra in controls:
    r = compute_B_trend('HBETA_BR', 0.3, 1.0, 'Hβ', resid_vars, resid_labels, extra)
    if r:
        print(f"  {r['resid']:>35} {r['N']:>8} | {r['rho_B']:>+7.3f} {r['p_B']:>8.4f} | "
              f"{r['rho_delta']:>+7.3f} {r['p_delta']:>8.4f}")

# ============================================================
# FULL LADDER WITH FeII CONTROL
# ============================================================
print(f"\n{'=' * 100}")
print("FULL IONIZATION LADDER: L_bol only vs L_bol + FeII_UV")
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

results_base = []
results_feii = []
results_full = []

for col, z_min, z_max, ip, label in LINE_DEFS:
    r1 = compute_B_trend(col, z_min, z_max, label, [lbol], ['L_bol'])
    r2 = compute_B_trend(col, z_min, z_max, label, [lbol, log_feii], ['L_bol', 'FeII'], feii_valid)
    r3 = compute_B_trend(col, z_min, z_max, label, [lbol, logledd, log_feii], ['L_bol', 'L/LEdd', 'FeII'], feii_valid)
    if r1:
        r1['ip'] = ip
        results_base.append(r1)
    if r2:
        r2['ip'] = ip
        results_feii.append(r2)
    if r3:
        r3['ip'] = ip
        results_full.append(r3)

print(f"\n  {'Line':>12} {'IP':>6} | {'ρ_B (base)':>11} {'ρ_B (+FeII)':>12} {'ρ_B (full)':>11} | {'Δ(base→FeII)':>14} {'Δ(base→full)':>14}")
print(f"  " + "-" * 95)

for r1, r2, r3 in zip(sorted(results_base, key=lambda x: x['ip']),
                        sorted(results_feii, key=lambda x: x['ip']),
                        sorted(results_full, key=lambda x: x['ip'])):
    d1 = r2['rho_B'] - r1['rho_B']
    d2 = r3['rho_B'] - r1['rho_B']
    print(f"  {r1['label']:>12} {r1['ip']:>6.1f} | {r1['rho_B']:>+11.3f} {r2['rho_B']:>+12.3f} {r3['rho_B']:>+11.3f} | {d1:>+14.3f} {d2:>+14.3f}")

# Meta-correlations
if len(results_base) >= 3:
    for label, res_list in [("Base (L_bol)", results_base), 
                              ("+FeII", results_feii),
                              ("+L/LEdd+FeII", results_full)]:
        ips = [r['ip'] for r in sorted(res_list, key=lambda x: x['ip'])]
        rhos = [r['rho_B'] for r in sorted(res_list, key=lambda x: x['ip'])]
        rho_meta, p_meta = spearmanr(ips, rhos)
        print(f"\n  Ladder ρ(IP, ρ_B) [{label}]: {rho_meta:+.3f}, p = {p_meta:.4f}")

# ============================================================
# VERDICT
# ============================================================
print(f"\n{'=' * 100}")
print("VERDICT")
print("=" * 100)

if len(mgii_results) >= 2:
    base = mgii_results[0]
    feii_ctrl = [r for r in mgii_results if 'FeII' in r['resid']]
    if feii_ctrl:
        best = feii_ctrl[-1]  # fullest control
        collapse = abs(best['rho_B']) < abs(base['rho_B']) * 0.5
        print(f"\n  MgII base:  ρ_B = {base['rho_B']:+.3f} (p = {base['p_B']:.4f})")
        print(f"  MgII +FeII: ρ_B = {best['rho_B']:+.3f} (p = {best['p_B']:.4f})")
        if collapse:
            print(f"  🔥 MgII divergence COLLAPSES with FeII control")
            print(f"     → BLR covering/structure (Eigenvector 1) IS the driver")
        elif abs(best['rho_B']) < abs(base['rho_B']) * 0.75:
            print(f"  ⚠️  MgII divergence WEAKENED but not collapsed")
        else:
            print(f"  ✗ MgII divergence PERSISTS even with FeII control")
