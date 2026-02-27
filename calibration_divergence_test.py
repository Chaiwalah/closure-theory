#!/usr/bin/env python3
"""
CALIBRATION DIVERGENCE TEST
============================
If local physics is universal, redundant diagnostic ratios measuring the same
physical quantity should agree at all redshifts. If local physics is emergent
and epoch-dependent, they should systematically diverge at high-z.

Prediction:
- Locked pairs: no divergence (Δρ ≈ 0)
- Diagnostic pairs: divergence proportional to sensitivity (Δρ < 0)
- If locked pairs ALSO diverge: it's SNR/selection, not the process
"""

import numpy as np
from astropy.io import fits
from scipy.stats import spearmanr
import json, os

FITS_PATH = '/root/clawd/data/sdss/dr16q_prop.fits'
OUT_DIR = '/root/clawd/projects/closure-theory/results_calibration_divergence'
os.makedirs(OUT_DIR, exist_ok=True)

def fix(arr):
    return arr.astype(arr.dtype.newbyteorder('='))

print("Loading DR16Q...")
data = fits.open(FITS_PATH)[1].data
z = fix(data['Z_DR16Q'])

# Line arrays: [wave_peak, wave_center, EW, FWHM, flux_peak, flux_total]
# idx: 2=EW, 3=FWHM, 5=flux_total
def get_line(name, idx):
    arr = fix(data[name])
    return arr[:, idx]

# Build line measurements
L = {}
for name in ['OIII5007','OIII5007C','HBETA','HBETA_BR','SII6718','NII6585',
             'OII3728','CIV','CIII_ALL','MGII','MGII_BR','LYA','HALPHA','HALPHA_BR']:
    try:
        L[name] = {
            'ew': get_line(name, 2),
            'fwhm': get_line(name, 3),
            'flux': get_line(name, 5),
        }
    except:
        pass

print(f"Loaded {len(L)} lines")

# Define test pairs: (name, type, series_a, series_b)
pairs = []

# LOCKED: OIII 5007 vs OIII 4959 equivalent (5007C is the core component)
# Actually OIII5007/OIII5007C isn't the doublet — 5007C is narrow core vs total
# Use OIII5007 flux vs OIII5007 EW as internal consistency instead

# LOCKED: Halpha vs Halpha_BR — broad component should track total
pairs.append(('Halpha_total vs Halpha_BR flux', 'semi-locked', L['HALPHA']['flux'], L['HALPHA_BR']['flux']))

# LOCKED: MGII vs MGII_BR
pairs.append(('MgII_total vs MgII_BR flux', 'semi-locked', L['MGII']['flux'], L['MGII_BR']['flux']))

# DIAGNOSTIC: OIII vs Hbeta (ionization parameter)
pairs.append(('OIII vs Hbeta EW (ionization)', 'diagnostic', L['OIII5007']['ew'], L['HBETA']['ew']))
pairs.append(('OIII vs Hbeta flux', 'diagnostic', L['OIII5007']['flux'], L['HBETA']['flux']))

# DIAGNOSTIC: NII vs SII (low-ionization gas)  
pairs.append(('NII vs SII EW (low-ion)', 'diagnostic', L['NII6585']['ew'], L['SII6718']['ew']))
pairs.append(('NII vs SII flux', 'diagnostic', L['NII6585']['flux'], L['SII6718']['flux']))

# DIAGNOSTIC: CIV vs CIII (high-ionization BLR)
pairs.append(('CIV vs CIII EW (high-ion)', 'diagnostic', L['CIV']['ew'], L['CIII_ALL']['ew']))
pairs.append(('CIV vs CIII flux', 'diagnostic', L['CIV']['flux'], L['CIII_ALL']['flux']))

# DIAGNOSTIC: Lya vs CIV (UV)
pairs.append(('Lya vs CIV EW (UV)', 'diagnostic', L['LYA']['ew'], L['CIV']['ew']))

# DIAGNOSTIC: OIII_FWHM vs Hbeta_FWHM (kinematic)
pairs.append(('OIII vs Hbeta FWHM (kinematic)', 'diagnostic', L['OIII5007']['fwhm'], L['HBETA']['fwhm']))

# DIAGNOSTIC: Halpha vs Hbeta (Balmer decrement — should be ~2.86 in Case B)
pairs.append(('Halpha vs Hbeta flux (Balmer dec)', 'diagnostic', L['HALPHA']['flux'], L['HBETA']['flux']))

# INTERNAL: OIII EW vs OIII flux (same line, different measure)
pairs.append(('OIII EW vs flux (internal)', 'internal', L['OIII5007']['ew'], L['OIII5007']['flux']))

# INTERNAL: CIV EW vs flux
pairs.append(('CIV EW vs flux (internal)', 'internal', L['CIV']['ew'], L['CIV']['flux']))

# CROSS-REGION: NLR (OIII) vs BLR (Hbeta_BR) — different spatial zones
pairs.append(('OIII vs Hbeta_BR flux (NLR vs BLR)', 'diagnostic', L['OIII5007']['flux'], L['HBETA_BR']['flux']))

# OII vs OIII — different ionization zones
pairs.append(('OII vs OIII EW', 'diagnostic', L['OII3728']['ew'], L['OIII5007']['ew']))

print(f"\nTesting {len(pairs)} pairs...")
print("="*90)

all_results = []

for name, ptype, a, b in pairs:
    mask = np.isfinite(a) & np.isfinite(b) & (a != 0) & (b != 0) & np.isfinite(z) & (z > 0.05)
    if mask.sum() < 200:
        print(f"  {name}: SKIP (n={mask.sum()})")
        continue
    
    zm, am, bm = z[mask], a[mask], b[mask]
    
    # Adaptive z-bins
    zlo, zhi = np.percentile(zm, [2, 98])
    zbins = np.linspace(zlo, zhi, 11)
    
    bin_results = []
    for i in range(len(zbins)-1):
        sel = (zm >= zbins[i]) & (zm < zbins[i+1])
        if sel.sum() > 30:
            rho, p = spearmanr(am[sel], bm[sel])
            bin_results.append({
                'z_mid': round((zbins[i]+zbins[i+1])/2, 3),
                'n': int(sel.sum()),
                'rho': round(rho, 4),
            })
    
    if len(bin_results) < 4:
        continue
    
    zz = [r['z_mid'] for r in bin_results]
    rr = [r['rho'] for r in bin_results]
    trend, tp = spearmanr(zz, rr)
    
    rho_low = bin_results[0]['rho']
    rho_high = bin_results[-1]['rho']
    delta = rho_high - rho_low
    
    result = {
        'pair': name, 'type': ptype,
        'n_total': int(mask.sum()),
        'rho_low_z': rho_low, 'rho_high_z': rho_high,
        'delta_rho': round(delta, 4),
        'trend_rho': round(trend, 4), 'trend_p': float(f'{tp:.2e}'),
        'bins': bin_results,
    }
    all_results.append(result)
    
    flag = "🔴" if delta < -0.05 and tp < 0.05 else ("🟢" if abs(delta) < 0.02 else "🟡")
    print(f"  {flag} [{ptype:11s}] {name:40s} ρ: {rho_low:+.3f} → {rho_high:+.3f}  Δ={delta:+.4f}  trend={trend:+.3f} p={tp:.1e}")

# === SUMMARY ===
print("\n" + "="*90)
print("SUMMARY: Divergence by type")
print("="*90)

for t in ['locked', 'semi-locked', 'internal', 'diagnostic']:
    deltas = [r['delta_rho'] for r in all_results if r['type'] == t]
    if deltas:
        print(f"  {t:12s}: mean Δρ = {np.mean(deltas):+.4f}, median = {np.median(deltas):+.4f}, n={len(deltas)}")

# === THE LADDER ===
print("\n" + "="*90)
print("DIVERGENCE LADDER (sorted by Δρ)")
print("="*90)
for r in sorted(all_results, key=lambda x: x['delta_rho']):
    print(f"  Δρ={r['delta_rho']:+.4f}  [{r['type']:11s}]  {r['pair']}")

# === KEY TEST: Is locked divergence < diagnostic divergence? ===
locked = [r['delta_rho'] for r in all_results if r['type'] in ('locked', 'semi-locked')]
diag = [r['delta_rho'] for r in all_results if r['type'] == 'diagnostic']
if locked and diag:
    print(f"\n  LOCKED/SEMI mean Δρ:     {np.mean(locked):+.4f}")
    print(f"  DIAGNOSTIC mean Δρ:      {np.mean(diag):+.4f}")
    separation = np.mean(diag) - np.mean(locked)
    print(f"  Separation:              {separation:+.4f}")
    if np.mean(diag) < np.mean(locked) - 0.02:
        print("  >>> DIAGNOSTIC DIVERGES MORE THAN LOCKED — CONSISTENT WITH EMERGENT PHYSICS")
    elif abs(np.mean(diag) - np.mean(locked)) < 0.02:
        print("  >>> NO DIFFERENCE — UNIVERSAL PHYSICS OR SNR EFFECT")
    else:
        print("  >>> LOCKED DIVERGES MORE — UNEXPECTED, CHECK DATA QUALITY")

# Save
with open(os.path.join(OUT_DIR, 'results.json'), 'w') as f:
    json.dump({'pairs': all_results}, f, indent=2, default=str)

print(f"\nSaved to {OUT_DIR}/results.json")
