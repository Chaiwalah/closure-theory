#!/usr/bin/env python3
"""
DECORRELATION THRESHOLD TEST
=============================
THE key prediction: if the medium degrades EW but not FWHM,
their mutual information should DECAY with redshift.

The z at which decorrelation accelerates = the threshold.
This is INDEPENDENT of the dispersion approach.
"""

import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.stats import spearmanr, pearsonr
import json
import os
import warnings
warnings.filterwarnings('ignore')

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data
z = data['Z_DR16Q']

lines = [
    ('HALPHA', 'HALPHA', 6563),
    ('HBETA', 'HBETA', 4861),
    ('MGII', 'MGII', 2798),
    ('CIV', 'CIV', 1549),
    ('CIII', 'CIII_ALL', 1909),
    ('OII', 'OII3728', 3728),
    ('OIII', 'OIII5007', 5007),
    ('NII', 'NII6585', 6585),
    ('SII', 'SII6718', 6718),
    ('LYA', 'LYA', 1216),
]

def sigmoid_decay(x, z0, k, r_high, r_low):
    """Correlation decays from r_high to r_low with sigmoid at z0"""
    return r_low + (r_high - r_low) / (1 + np.exp(k * (x - z0)))

print("=" * 70)
print("DECORRELATION THRESHOLD TEST")
print("EW-FWHM correlation decay → maps the medium's effect")
print("=" * 70)

all_results = {}

for name, col_name, rest_wav in lines:
    try:
        col = data[col_name]
    except:
        continue
    
    if len(col.shape) < 2 or col.shape[1] < 5:
        continue
    
    ew = col[:, 2]
    fwhm = col[:, 4]
    
    valid = (np.isfinite(ew) & np.isfinite(fwhm) & 
             (ew > 0) & (fwhm > 0) & 
             (ew < 500) & (fwhm < 30000) &
             (z > 0.05))
    
    if valid.sum() < 3000:
        continue
    
    z_v = z[valid]
    ew_v = ew[valid]
    fwhm_v = fwhm[valid]
    
    # Fine binning — 30 bins
    z_min, z_max = np.percentile(z_v, [2, 98])
    n_bins = 30
    z_edges = np.linspace(z_min, z_max, n_bins + 1)
    
    corr_spearman = []
    corr_pearson = []
    z_centers = []
    n_per_bin = []
    
    for i in range(n_bins):
        mask = (z_v >= z_edges[i]) & (z_v < z_edges[i+1])
        if mask.sum() < 100:
            continue
        
        rs, ps = spearmanr(ew_v[mask], fwhm_v[mask])
        rp, pp = pearsonr(ew_v[mask], fwhm_v[mask])
        
        corr_spearman.append(rs)
        corr_pearson.append(rp)
        z_centers.append((z_edges[i] + z_edges[i+1]) / 2)
        n_per_bin.append(mask.sum())
    
    corr_s = np.array(corr_spearman)
    corr_p = np.array(corr_pearson)
    z_c = np.array(z_centers)
    
    if len(z_c) < 6:
        continue
    
    # Overall trend
    r_trend, p_trend = spearmanr(z_c, corr_s)
    
    # Try sigmoid fit
    try:
        p0 = [np.median(z_c), 3, corr_s[0], corr_s[-1]]
        popt, pcov = curve_fit(sigmoid_decay, z_c, corr_s, p0=p0, maxfev=10000)
        z0_fit = popt[0]
        k_fit = popt[1]
        r_high = popt[2]
        r_low = popt[3]
        
        pred = sigmoid_decay(z_c, *popt)
        ss_res = np.sum((corr_s - pred)**2)
        ss_tot = np.sum((corr_s - np.mean(corr_s))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        
        # Amplitude of decay
        decay_amp = r_high - r_low
        
    except:
        z0_fit = None
        r2 = 0
        decay_amp = corr_s[0] - corr_s[-1]
    
    # Compute the "first half vs second half" test
    mid = len(z_c) // 2
    r_first = np.mean(corr_s[:mid])
    r_second = np.mean(corr_s[mid:])
    
    print(f"\n{'='*50}")
    print(f"{name} ({rest_wav}Å) — N = {valid.sum():,}")
    print(f"  z range: {z_min:.2f} – {z_max:.2f}")
    print(f"  Correlation (first half):  r = {r_first:.4f}")
    print(f"  Correlation (second half): r = {r_second:.4f}")
    print(f"  Decay: Δr = {r_first - r_second:+.4f}")
    print(f"  Overall trend: r = {r_trend:.3f} (p = {p_trend:.2e})")
    
    if z0_fit is not None and 0 < z0_fit < 5:
        print(f"  Sigmoid z₀ = {z0_fit:.3f} (R² = {r2:.3f})")
        print(f"  Decay amplitude: {decay_amp:.3f}")
    
    result = {
        'rest_wav': rest_wav,
        'n_objects': int(valid.sum()),
        'z_range': [float(z_min), float(z_max)],
        'r_first_half': float(r_first),
        'r_second_half': float(r_second),
        'decay_delta_r': float(r_first - r_second),
        'trend_r': float(r_trend),
        'trend_p': float(p_trend),
        'z0_sigmoid': float(z0_fit) if z0_fit and 0 < z0_fit < 5 else None,
        'r2_sigmoid': float(r2) if z0_fit else None,
        'bin_z': [float(x) for x in z_c],
        'bin_corr': [float(x) for x in corr_s],
    }
    all_results[name] = result

# ============================================================
# CROSS-LINE COMPARISON
# ============================================================
print("\n" + "=" * 70)
print("CROSS-LINE COMPARISON")
print("=" * 70)

# Which lines show significant decorrelation?
sig_lines = [(name, d) for name, d in all_results.items() 
             if d['trend_p'] < 0.05 and d['trend_r'] < 0]

print(f"\nLines with significant EW-FWHM decorrelation (p<0.05):")
print(f"{'Line':>10} {'λ':>6} {'Δr':>8} {'trend_r':>8} {'p':>12} {'z₀':>8}")
print("-" * 58)
for name, d in sorted(sig_lines, key=lambda x: x[1]['rest_wav']):
    z0_str = f"{d['z0_sigmoid']:.3f}" if d['z0_sigmoid'] else "N/A"
    print(f"{name:>10} {d['rest_wav']:>5.0f}Å {d['decay_delta_r']:>+7.4f} {d['trend_r']:>8.3f} {d['trend_p']:>12.2e} {z0_str:>8}")

# Lines WITHOUT decorrelation (controls)
nosig = [(name, d) for name, d in all_results.items()
         if d['trend_p'] >= 0.05 or d['trend_r'] >= 0]

if nosig:
    print(f"\nLines WITHOUT decorrelation (controls):")
    for name, d in nosig:
        print(f"  {name} ({d['rest_wav']}Å): trend r = {d['trend_r']:.3f}, p = {d['trend_p']:.2e}")

# Does z₀ correlate with rest wavelength?
sig_with_z0 = [(d['rest_wav'], d['z0_sigmoid']) for _, d in sig_lines if d['z0_sigmoid']]
if len(sig_with_z0) >= 3:
    wavs = np.array([x[0] for x in sig_with_z0])
    z0s = np.array([x[1] for x in sig_with_z0])
    r_wz, p_wz = spearmanr(wavs, z0s)
    print(f"\nz₀ vs rest wavelength: r = {r_wz:.3f} (p = {p_wz:.3e})")
    if r_wz > 0:
        print("→ Bluer lines threshold EARLIER → chromatic medium ✓")
    else:
        print("→ No clear wavelength ordering")

# ============================================================
# PREDICTION TABLE
# ============================================================
print("\n" + "=" * 70)
print("PREDICTION TABLE")
print("=" * 70)
print(f"""
Dark energy transition:     z = 0.632
Web coherence prediction:   z = 0.730
SNe Ia sigmoid:             z = 0.820
Quasar EW coupling:         z = 1.050

Decorrelation thresholds found:""")

for name, d in sorted(all_results.items(), key=lambda x: x[1].get('z0_sigmoid') or 99):
    if d.get('z0_sigmoid') and d['trend_p'] < 0.1:
        status = "✓ SIGNIFICANT" if d['trend_p'] < 0.05 else "marginal"
        print(f"  {name:>10} ({d['rest_wav']}Å): z₀ = {d['z0_sigmoid']:.3f}  [{status}]")

# Save
os.makedirs('results_threshold_map', exist_ok=True)
with open('results_threshold_map/decorrelation.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print("\nSaved to results_threshold_map/decorrelation.json")
hdu.close()
