#!/usr/bin/env python3
"""
WAVE vs SIGMOID TEST — Does the handoff oscillate?

Humza's intuition: the optical→UV handoff might follow a trigonometric
pattern (cosine, sine, cotangent) rather than a monotonic sigmoid.

If periodic: standing wave / resonance in correlation space
If sigmoid: one-time phase transition
If neither: something else entirely

Test: Fit MgII's allegiance (coupling to Hβ vs CIII vs CIV) as function of z
with sigmoid, cosine, and polynomial models. Compare fits.
"""

import numpy as np
from astropy.io import fits
from scipy import stats, optimize
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_wave_test', exist_ok=True)

print("=" * 70)
print("WAVE vs SIGMOID — Does the correlation structure oscillate?")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']

# Extract key observables
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
ciii_ew = d['CIII_BR'][:, 2]
civ_ew = d['CIV'][:, 2]
ha_ew = d['HALPHA_BR'][:, 2]
lya_ew = d['LYA'][:, 2]

# Fine z-bins for high resolution
z_edges = np.arange(0.3, 2.5, 0.05)
z_mids = (z_edges[:-1] + z_edges[1:]) / 2

# ============================================================
# Compute correlations at fine z resolution
# ============================================================
pairs = [
    ('MgII↔Hβ', mg_ew, hb_ew),
    ('MgII↔CIII', mg_ew, ciii_ew),
    ('MgII↔CIV', mg_ew, civ_ew),
    ('Hα↔Hβ', ha_ew, hb_ew),
    ('CIII↔CIV', ciii_ew, civ_ew),
    ('MgII↔LyA', mg_ew, lya_ew),
]

all_series = {}

for pair_name, v1, v2 in pairs:
    corrs = []
    z_valid = []
    ns = []
    
    for i in range(len(z_edges) - 1):
        zlo, zhi = z_edges[i], z_edges[i+1]
        zmask = (z_all >= zlo) & (z_all < zhi)
        
        d1 = v1[zmask]
        d2 = v2[zmask]
        valid = np.isfinite(d1) & np.isfinite(d2) & (d1 != 0) & (d2 != 0)
        
        if valid.sum() > 50:
            r, p = stats.spearmanr(d1[valid], d2[valid])
            corrs.append(r)
            z_valid.append(z_mids[i])
            ns.append(int(valid.sum()))
    
    if len(corrs) > 5:
        all_series[pair_name] = {
            'z': np.array(z_valid),
            'r': np.array(corrs),
            'n': ns,
        }

# ============================================================
# Model fitting: sigmoid vs cosine vs polynomial
# ============================================================
print(f"\n{'='*70}")
print("MODEL COMPARISON — What shape is the handoff?")
print(f"{'='*70}")

def sigmoid(z, a, b, z0, k):
    """Sigmoid: r = a + b / (1 + exp(-k*(z-z0)))"""
    return a + b / (1 + np.exp(-k * (z - z0)))

def cosine_model(z, a, b, omega, phi):
    """Cosine: r = a + b * cos(omega*z + phi)"""
    return a + b * np.cos(omega * z + phi)

def damped_cosine(z, a, b, omega, phi, tau):
    """Damped cosine: r = a + b * cos(omega*z + phi) * exp(-z/tau)"""
    return a + b * np.cos(omega * z + phi) * np.exp(-z / tau)

for pair_name, data in all_series.items():
    z = data['z']
    r = data['r']
    
    if len(z) < 8:
        continue
    
    print(f"\n  {pair_name} ({len(z)} z-bins):")
    
    # Fit sigmoid
    try:
        popt_sig, _ = optimize.curve_fit(sigmoid, z, r, 
                                          p0=[0.8, -0.8, 1.0, -5],
                                          maxfev=10000)
        r_pred_sig = sigmoid(z, *popt_sig)
        ss_res_sig = np.sum((r - r_pred_sig)**2)
        ss_tot = np.sum((r - np.mean(r))**2)
        r2_sig = 1 - ss_res_sig / ss_tot if ss_tot > 0 else 0
        aic_sig = len(z) * np.log(ss_res_sig / len(z)) + 2 * 4
        print(f"    Sigmoid:  R²={r2_sig:.4f}  AIC={aic_sig:.1f}  z₀={popt_sig[2]:.2f}  k={popt_sig[3]:.1f}")
    except:
        r2_sig = -999
        aic_sig = 999
        print(f"    Sigmoid:  FAILED TO FIT")
    
    # Fit cosine
    try:
        popt_cos, _ = optimize.curve_fit(cosine_model, z, r, 
                                          p0=[0.4, 0.4, 3.0, 0],
                                          maxfev=10000)
        r_pred_cos = cosine_model(z, *popt_cos)
        ss_res_cos = np.sum((r - r_pred_cos)**2)
        r2_cos = 1 - ss_res_cos / ss_tot if ss_tot > 0 else 0
        aic_cos = len(z) * np.log(ss_res_cos / len(z)) + 2 * 4
        period = 2 * np.pi / abs(popt_cos[2])
        print(f"    Cosine:   R²={r2_cos:.4f}  AIC={aic_cos:.1f}  ω={popt_cos[2]:.2f}  period={period:.2f}")
    except:
        r2_cos = -999
        aic_cos = 999
        print(f"    Cosine:   FAILED TO FIT")
    
    # Fit damped cosine
    try:
        popt_dc, _ = optimize.curve_fit(damped_cosine, z, r, 
                                         p0=[0.4, 0.4, 3.0, 0, 2.0],
                                         maxfev=10000)
        r_pred_dc = damped_cosine(z, *popt_dc)
        ss_res_dc = np.sum((r - r_pred_dc)**2)
        r2_dc = 1 - ss_res_dc / ss_tot if ss_tot > 0 else 0
        aic_dc = len(z) * np.log(ss_res_dc / len(z)) + 2 * 5
        period_dc = 2 * np.pi / abs(popt_dc[2])
        print(f"    Damped:   R²={r2_dc:.4f}  AIC={aic_dc:.1f}  ω={popt_dc[2]:.2f}  period={period_dc:.2f}  τ={popt_dc[4]:.2f}")
    except:
        r2_dc = -999
        aic_dc = 999
        print(f"    Damped:   FAILED TO FIT")
    
    # Polynomial (cubic)
    try:
        coeffs = np.polyfit(z, r, 3)
        r_pred_poly = np.polyval(coeffs, z)
        ss_res_poly = np.sum((r - r_pred_poly)**2)
        r2_poly = 1 - ss_res_poly / ss_tot if ss_tot > 0 else 0
        aic_poly = len(z) * np.log(ss_res_poly / len(z)) + 2 * 4
        print(f"    Cubic:    R²={r2_poly:.4f}  AIC={aic_poly:.1f}")
    except:
        r2_poly = -999
        aic_poly = 999
        print(f"    Cubic:    FAILED TO FIT")
    
    # Winner
    models = {'sigmoid': aic_sig, 'cosine': aic_cos, 'damped_cos': aic_dc, 'cubic': aic_poly}
    winner = min(models, key=models.get)
    print(f"    ★ WINNER: {winner} (lowest AIC)")

# ============================================================
# Fine-grained correlation evolution (the actual shape)
# ============================================================
print(f"\n{'='*70}")
print("FINE-GRAINED CORRELATION EVOLUTION")
print(f"{'='*70}")

for pair_name in ['MgII↔Hβ', 'MgII↔CIII', 'MgII↔CIV']:
    if pair_name not in all_series:
        continue
    data = all_series[pair_name]
    print(f"\n  {pair_name}:")
    for z_val, r_val, n_val in zip(data['z'], data['r'], data['n']):
        bar_len = int(abs(r_val) * 30)
        bar = "█" * bar_len if r_val > 0 else "▒" * bar_len
        sign = "+" if r_val > 0 else "-"
        print(f"    z={z_val:.2f}: {sign}{abs(r_val):.3f} {bar} (N={n_val})")

# ============================================================
# Check for oscillation in residuals after sigmoid fit
# ============================================================
print(f"\n{'='*70}")
print("RESIDUAL OSCILLATION — Does the sigmoid miss periodic structure?")
print(f"{'='*70}")

for pair_name in ['MgII↔Hβ', 'MgII↔CIII']:
    if pair_name not in all_series:
        continue
    data = all_series[pair_name]
    z = data['z']
    r = data['r']
    
    if len(z) < 8:
        continue
    
    try:
        popt, _ = optimize.curve_fit(sigmoid, z, r, p0=[0.8, -0.8, 1.0, -5], maxfev=10000)
        residuals = r - sigmoid(z, *popt)
        
        # Autocorrelation of residuals (are they periodic?)
        n = len(residuals)
        autocorr = np.correlate(residuals - np.mean(residuals), 
                                residuals - np.mean(residuals), mode='full')
        autocorr = autocorr[n-1:] / autocorr[n-1]  # normalize
        
        # FFT of residuals
        fft = np.fft.rfft(residuals)
        power = np.abs(fft)**2
        freqs = np.fft.rfftfreq(n, d=np.mean(np.diff(z)))
        
        peak_freq_idx = np.argmax(power[1:]) + 1  # skip DC
        peak_freq = freqs[peak_freq_idx]
        peak_period = 1.0 / peak_freq if peak_freq > 0 else np.inf
        
        print(f"\n  {pair_name}:")
        print(f"    Residual std: {np.std(residuals):.4f}")
        print(f"    Peak frequency: {peak_freq:.2f} per z-unit")
        print(f"    Peak period: Δz = {peak_period:.2f}")
        print(f"    Autocorr lag-1: {autocorr[1]:.3f}")
        print(f"    Autocorr lag-2: {autocorr[2]:.3f}")
        
        # Is there significant oscillation?
        if abs(autocorr[1]) > 0.3:
            print(f"    ⚠️ SIGNIFICANT autocorrelation — sigmoid may miss periodic structure")
        else:
            print(f"    ✓ No significant periodicity in residuals")
    except:
        print(f"\n  {pair_name}: fit failed")

# Save
results_out = {}
for pair_name, data in all_series.items():
    results_out[pair_name] = {
        'z': data['z'].tolist(),
        'r': data['r'].tolist(),
        'n': data['n'],
    }

with open('results_wave_test/wave_test_results.json', 'w') as f:
    json.dump(results_out, f, indent=2)

print(f"\n\nResults saved to results_wave_test/")
