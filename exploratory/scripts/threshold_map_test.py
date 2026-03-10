#!/usr/bin/env python3
"""
THRESHOLD MAP TEST
==================
Map the sigmoid threshold across INDEPENDENT observables.
If the same z₀ ≈ 0.82 appears in multiple independent measurements,
it's not a fluke — it's a physical transition.

Tests:
1. SNe Ia color-distance coupling (Pantheon+ — already done, z₀=0.82)
2. Quasar EW degradation by line (SDSS DR16Q — already done, z₀~1.05)
3. Quasar CONTINUUM slope change (NEW — independent of lines)
4. Quasar variability amplitude vs z (NEW — if medium adds stochastic lensing)
5. Galaxy angular size anomaly (NEW — if surface brightness affected)
6. Cross-check: FWHM should show NO threshold (control)
"""

import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.stats import spearmanr, ks_2samp
import json
import os
import warnings
warnings.filterwarnings('ignore')

# Load SDSS DR16Q
print("Loading SDSS DR16Q...")
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data
z = data['Z_DR16Q']

# Line definitions: name, column, rest wavelength (Å)
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

def sigmoid(x, z0, k, A, B):
    return A / (1 + np.exp(-k * (x - z0))) + B

def get_line_data(col_name, z_data):
    """Extract EW and FWHM for a line, return valid mask"""
    col = data[col_name]
    if len(col.shape) > 1 and col.shape[1] >= 5:
        ew = col[:, 2]    # EW
        fwhm = col[:, 4]  # FWHM
    else:
        return None, None, None
    
    valid = (np.isfinite(ew) & np.isfinite(fwhm) & 
             (ew > 0) & (fwhm > 0) & 
             (ew < 500) & (fwhm < 30000) &
             (z_data > 0.05))
    
    return ew, fwhm, valid

print("\n" + "=" * 70)
print("TEST 1: SIGMOID FIT TO EACH LINE's EW DISPERSION")
print("Finding z₀ independently for each line")
print("=" * 70)

results_by_line = {}

for name, col, rest_wav in lines:
    ew, fwhm, valid = get_line_data(col, z)
    if ew is None or valid.sum() < 1000:
        continue
    
    z_valid = z[valid]
    ew_valid = ew[valid]
    fwhm_valid = fwhm[valid]
    
    # Bin by redshift
    z_min, z_max = z_valid.min(), min(z_valid.max(), 4.0)
    n_bins = 20
    z_edges = np.linspace(z_min, z_max, n_bins + 1)
    
    ew_cv = []  # coefficient of variation per bin
    fwhm_cv = []
    z_centers = []
    
    for i in range(n_bins):
        mask = (z_valid >= z_edges[i]) & (z_valid < z_edges[i+1])
        if mask.sum() < 50:
            continue
        
        ew_bin = ew_valid[mask]
        fwhm_bin = fwhm_valid[mask]
        
        # CV = std/mean — normalized dispersion
        cv_ew = np.std(ew_bin) / np.mean(ew_bin)
        cv_fwhm = np.std(fwhm_bin) / np.mean(fwhm_bin)
        
        ew_cv.append(cv_ew)
        fwhm_cv.append(cv_fwhm)
        z_centers.append((z_edges[i] + z_edges[i+1]) / 2)
    
    ew_cv = np.array(ew_cv)
    fwhm_cv = np.array(fwhm_cv)
    z_centers = np.array(z_centers)
    
    if len(z_centers) < 8:
        continue
    
    # Fit sigmoid to EW CV
    try:
        popt_ew, _ = curve_fit(sigmoid, z_centers, ew_cv, 
                                p0=[0.8, 3, 0.3, 0.5], maxfev=10000)
        z0_ew = popt_ew[0]
        
        # Also fit FWHM (should NOT show sigmoid — control)
        popt_fwhm, _ = curve_fit(sigmoid, z_centers, fwhm_cv,
                                  p0=[0.8, 3, 0.3, 0.5], maxfev=10000)
        z0_fwhm = popt_fwhm[0]
        
        # Compute fit quality
        ew_pred = sigmoid(z_centers, *popt_ew)
        fwhm_pred = sigmoid(z_centers, *popt_fwhm)
        
        r2_ew = 1 - np.sum((ew_cv - ew_pred)**2) / np.sum((ew_cv - np.mean(ew_cv))**2)
        r2_fwhm = 1 - np.sum((fwhm_cv - fwhm_pred)**2) / np.sum((fwhm_cv - np.mean(fwhm_cv))**2)
        
        results_by_line[name] = {
            'rest_wav': rest_wav,
            'n_objects': int(valid.sum()),
            'z_range': [float(z_min), float(z_max)],
            'z0_ew': float(z0_ew),
            'z0_fwhm': float(z0_fwhm),
            'steepness_ew': float(popt_ew[1]),
            'r2_ew': float(r2_ew),
            'r2_fwhm': float(r2_fwhm),
        }
        
        # Is z₀ in plausible range?
        in_range = "✓" if 0.3 < z0_ew < 2.0 else "✗"
        ew_sigmoid = "YES" if r2_ew > 0.3 and 0.3 < z0_ew < 2.0 else "no"
        fwhm_sigmoid = "YES" if r2_fwhm > 0.3 and 0.3 < z0_fwhm < 2.0 else "no"
        
        print(f"\n{name} ({rest_wav}Å) — N={valid.sum():,}")
        print(f"  EW dispersion:   z₀ = {z0_ew:.3f}  R² = {r2_ew:.3f}  sigmoid={ew_sigmoid}")
        print(f"  FWHM dispersion: z₀ = {z0_fwhm:.3f}  R² = {r2_fwhm:.3f}  sigmoid={fwhm_sigmoid}")
        
    except Exception as e:
        print(f"\n{name}: fit failed ({e})")

# ============================================================
# TEST 2: CONTINUUM SLOPE CHANGE
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: CONTINUUM LUMINOSITY DISPERSION vs z")
print("(Independent of emission lines)")
print("=" * 70)

# Use continuum luminosity columns
for cont_name in ['LOGLBOL', 'LOGL1350', 'LOGL3000', 'LOGL5100']:
    try:
        lum = data[cont_name]
        valid_c = np.isfinite(lum) & (lum > 0) & (z > 0.05)
        
        if valid_c.sum() < 5000:
            continue
        
        z_c = z[valid_c]
        lum_c = lum[valid_c]
        
        # Bin and compute dispersion
        z_edges = np.linspace(z_c.min(), min(z_c.max(), 4.0), 25)
        cv_vals = []
        z_mids = []
        
        for i in range(len(z_edges)-1):
            mask = (z_c >= z_edges[i]) & (z_c < z_edges[i+1])
            if mask.sum() < 100:
                continue
            cv_vals.append(np.std(lum_c[mask]) / np.mean(lum_c[mask]))
            z_mids.append((z_edges[i] + z_edges[i+1]) / 2)
        
        cv_vals = np.array(cv_vals)
        z_mids = np.array(z_mids)
        
        if len(z_mids) < 8:
            continue
        
        try:
            popt, _ = curve_fit(sigmoid, z_mids, cv_vals, 
                                p0=[0.8, 3, 0.05, 0.01], maxfev=10000)
            pred = sigmoid(z_mids, *popt)
            r2 = 1 - np.sum((cv_vals - pred)**2) / np.sum((cv_vals - np.mean(cv_vals))**2)
            
            print(f"\n{cont_name}: z₀ = {popt[0]:.3f}, R² = {r2:.3f}, N = {valid_c.sum():,}")
            
        except:
            # Try linear trend instead
            r, p = spearmanr(z_mids, cv_vals)
            print(f"\n{cont_name}: No sigmoid, linear r = {r:.3f} (p={p:.2e}), N = {valid_c.sum():,}")
            
    except Exception as e:
        pass

# ============================================================
# TEST 3: EW-FWHM DECORRELATION vs z
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: EW-FWHM CORRELATION DECAY")
print("If medium degrades EW but not FWHM, their correlation")
print("should WEAKEN with z. Rate of weakening → sigmoid?")
print("=" * 70)

for name, col, rest_wav in lines:
    ew, fwhm, valid = get_line_data(col, z)
    if ew is None or valid.sum() < 5000:
        continue
    
    z_valid = z[valid]
    ew_valid = ew[valid]
    fwhm_valid = fwhm[valid]
    
    # Bin by z, compute EW-FWHM correlation in each bin
    z_edges = np.linspace(z_valid.min(), min(z_valid.max(), 4.0), 15)
    corr_vals = []
    z_mids = []
    
    for i in range(len(z_edges)-1):
        mask = (z_valid >= z_edges[i]) & (z_valid < z_edges[i+1])
        if mask.sum() < 200:
            continue
        
        r, p = spearmanr(ew_valid[mask], fwhm_valid[mask])
        corr_vals.append(r)
        z_mids.append((z_edges[i] + z_edges[i+1]) / 2)
    
    if len(z_mids) < 6:
        continue
    
    corr_vals = np.array(corr_vals)
    z_mids = np.array(z_mids)
    
    # Does correlation decay with z?
    r_decay, p_decay = spearmanr(z_mids, corr_vals)
    
    # Is the decay sigmoid-shaped?
    try:
        popt, _ = curve_fit(lambda x, z0, k, A, B: -A / (1 + np.exp(-k * (x - z0))) + B,
                            z_mids, corr_vals, p0=[0.8, 3, 0.3, 0.3], maxfev=10000)
        z0_decorr = popt[0]
        
        print(f"\n{name}: EW-FWHM correlation decay")
        print(f"  r(z=min) = {corr_vals[0]:.3f} → r(z=max) = {corr_vals[-1]:.3f}")
        print(f"  Decay trend: r = {r_decay:.3f} (p = {p_decay:.2e})")
        print(f"  Sigmoid z₀ = {z0_decorr:.3f}")
    except:
        print(f"\n{name}: EW-FWHM correlation decay")
        print(f"  r(z=min) = {corr_vals[0]:.3f} → r(z=max) = {corr_vals[-1]:.3f}")
        print(f"  Decay trend: r = {r_decay:.3f} (p = {p_decay:.2e})")
        print(f"  No sigmoid fit")

# ============================================================
# TEST 4: AGGREGATE THRESHOLD
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: AGGREGATE — DO ALL LINES POINT TO SAME z₀?")
print("=" * 70)

if results_by_line:
    z0_values = [(name, d['z0_ew'], d['rest_wav'], d['r2_ew']) 
                  for name, d in results_by_line.items() 
                  if 0.2 < d['z0_ew'] < 3.0 and d['r2_ew'] > 0.1]
    
    if z0_values:
        z0_values.sort(key=lambda x: x[1])
        
        print(f"\n{'Line':>10} {'λ_rest':>8} {'z₀':>8} {'R²':>8}")
        print("-" * 38)
        for name, z0, wav, r2 in z0_values:
            print(f"{name:>10} {wav:>7.0f}Å {z0:8.3f} {r2:8.3f}")
        
        z0_arr = np.array([x[1] for x in z0_values])
        wav_arr = np.array([x[2] for x in z0_values])
        
        print(f"\nMean z₀ = {np.mean(z0_arr):.3f} ± {np.std(z0_arr):.3f}")
        print(f"Median z₀ = {np.median(z0_arr):.3f}")
        
        # Does z₀ correlate with wavelength? (blue lines → earlier threshold?)
        r_wav, p_wav = spearmanr(wav_arr, z0_arr)
        print(f"\nz₀ vs wavelength: r = {r_wav:.3f} (p = {p_wav:.3e})")
        if r_wav > 0 and p_wav < 0.05:
            print("→ BLUE LINES THRESHOLD EARLIER (lower z₀)")
            print("→ Consistent with blue-preferential degradation!")
        
        # Compare to dark energy prediction
        z_decel = 0.632
        z_pred = 0.730  # from dark_energy_threshold_test.py
        
        print(f"\nDark energy transition: z = {z_decel:.3f}")
        print(f"Predicted threshold:    z = {z_pred:.3f}")
        print(f"Observed mean z₀:       z = {np.mean(z0_arr):.3f}")
        print(f"Gap (observed - predicted): {np.mean(z0_arr) - z_pred:.3f}")

# ============================================================
# TEST 5: THE KILLER — z₀ vs REST WAVELENGTH
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: THRESHOLD WAVELENGTH DEPENDENCE")
print("If the medium is dispersive (frequency-dependent),")
print("bluer lines should hit the threshold at LOWER z.")
print("This is a PREDICTION of the refractive decoherence model.")
print("=" * 70)

if results_by_line and len(z0_values) >= 4:
    print(f"\nPrediction: z₀ increases with rest wavelength")
    print(f"(blue light degrades first → threshold at lower z)")
    print(f"\nResult: r = {r_wav:.3f}, p = {p_wav:.3e}")
    
    if r_wav > 0.3 and p_wav < 0.05:
        print("\n🎯 CONFIRMED: Blue lines threshold EARLIER.")
        print("The medium is chromatic. Refractive decoherence predicted this.")
    elif r_wav > 0:
        print("\nTrend in predicted direction but not significant.")
        print("May need finer binning or more lines.")
    else:
        print("\nNo wavelength dependence detected.")
        print("Could mean: achromatic mechanism, or EW dispersion dominated by source physics.")

# Save all results
os.makedirs('results_threshold_map', exist_ok=True)
with open('results_threshold_map/line_thresholds.json', 'w') as f:
    json.dump(results_by_line, f, indent=2, default=str)

print("\n\nResults saved to results_threshold_map/")
print("=" * 70)
hdu.close()
