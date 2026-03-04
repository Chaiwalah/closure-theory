#!/usr/bin/env python3
"""
UNIVERSAL EQUATION v2
=====================
Key insight from v1: wavelength alone doesn't collapse the data.
Lines differ by TYPE — diagnostic sensitivity matters.

Try: r(z) = r₀(line) - P(line) × f(z)
where P(line) = susceptibility from doublet ladder
and f(z) = universal degradation function

Also try: separate BLR vs NLR lines.
"""

import numpy as np
from astropy.io import fits
from scipy.stats import spearmanr
from scipy.optimize import curve_fit, minimize
import warnings
warnings.filterwarnings('ignore')

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data
z = data['Z_DR16Q']

# Line susceptibility from doublet ladder (Feb 23 result)
# Higher = more diagnostic = more vulnerable
line_info = {
    'LYA':    {'col': 'LYA',      'wav': 1216, 'P': 0.8, 'type': 'BLR'},
    'CIV':    {'col': 'CIV',      'wav': 1549, 'P': 0.6, 'type': 'BLR'},
    'CIII':   {'col': 'CIII_ALL', 'wav': 1909, 'P': 0.5, 'type': 'BLR'},
    'MGII':   {'col': 'MGII',     'wav': 2798, 'P': 0.3, 'type': 'BLR'},
    'OII':    {'col': 'OII3728',  'wav': 3728, 'P': 0.4, 'type': 'NLR'},
    'HBETA':  {'col': 'HBETA',    'wav': 4861, 'P': 0.3, 'type': 'BLR'},
    'OIII':   {'col': 'OIII5007', 'wav': 5007, 'P': 0.0, 'type': 'NLR'},
    'HALPHA': {'col': 'HALPHA',   'wav': 6563, 'P': 0.3, 'type': 'BLR'},
    'NII':    {'col': 'NII6585',  'wav': 6585, 'P': 0.0, 'type': 'NLR'},
    'SII':    {'col': 'SII6718',  'wav': 6718, 'P': 0.7, 'type': 'NLR'},
}

# Extract data
all_line_data = {}
for name, info in line_info.items():
    try:
        col = data[info['col']]
    except:
        continue
    if len(col.shape) < 2 or col.shape[1] < 5:
        continue
    
    ew, fw = col[:,2], col[:,4]
    v = (np.isfinite(ew) & np.isfinite(fw) & (ew>0) & (fw>0) & 
         (ew<500) & (fw<30000) & (z>0.05))
    if v.sum() < 3000:
        continue
    
    zv, ewv, fwv = z[v], ew[v], fw[v]
    idx = np.argsort(zv)
    n = len(idx)
    nbins = min(20, n // 500)
    bs = n // nbins
    
    zmeds, corrs = [], []
    for i in range(nbins):
        sl = idx[i*bs:(i+1)*bs]
        r, _ = spearmanr(ewv[sl], fwv[sl])
        zmeds.append(np.median(zv[sl]))
        corrs.append(r)
    
    all_line_data[name] = {
        'z': np.array(zmeds),
        'r': np.array(corrs),
        'wav': info['wav'],
        'P': info['P'],
        'type': info['type'],
        'N': v.sum(),
    }

# ============================================================
# APPROACH 1: Measure dr/dz for each line, see if it scales
# with any line property
# ============================================================
print("=" * 70)
print("APPROACH 1: DECORRELATION RATE vs LINE PROPERTIES")
print("Does the RATE of decorrelation scale with something?")
print("=" * 70)

rates = []
for name, d in all_line_data.items():
    # Linear fit to get slope
    if len(d['z']) < 5: continue
    slope = np.polyfit(d['z'], d['r'], 1)[0]
    rates.append({
        'name': name, 'slope': slope, 'wav': d['wav'], 
        'P': d['P'], 'type': d['type'],
        'z_range': d['z'][-1] - d['z'][0],
    })

print(f"\n{'Line':>8} {'type':>4} {'λ':>6} {'P':>5} {'dr/dz':>8} {'Direction':>12}")
print("-" * 50)
for r in sorted(rates, key=lambda x: x['slope']):
    direction = "DECAYS" if r['slope'] < -0.01 else "GROWS" if r['slope'] > 0.01 else "FLAT"
    print(f"{r['name']:>8} {r['type']:>4} {r['wav']:>5}Å {r['P']:>5.1f} {r['slope']:>+8.4f} {direction:>12}")

# Correlate slope with wavelength and susceptibility
slopes = np.array([r['slope'] for r in rates])
wavs = np.array([r['wav'] for r in rates])
Ps = np.array([r['P'] for r in rates])

r_wav, p_wav = spearmanr(wavs, slopes)
r_P, p_P = spearmanr(Ps, slopes)

print(f"\nSlope vs wavelength: r = {r_wav:+.3f} (p = {p_wav:.3e})")
print(f"Slope vs susceptibility P: r = {r_P:+.3f} (p = {p_P:.3e})")

# ============================================================
# APPROACH 2: Split BLR vs NLR
# ============================================================
print(f"\n{'='*70}")
print("APPROACH 2: BLR vs NLR LINES")
print("='*70")

blr_slopes = [r for r in rates if r['type'] == 'BLR']
nlr_slopes = [r for r in rates if r['type'] == 'NLR']

print(f"\nBLR lines (broad, complex emissivity):")
for r in blr_slopes:
    print(f"  {r['name']:>8} ({r['wav']}Å): dr/dz = {r['slope']:+.4f}")
print(f"  Mean dr/dz = {np.mean([r['slope'] for r in blr_slopes]):+.4f}")

print(f"\nNLR lines (narrow, forbidden):")
for r in nlr_slopes:
    print(f"  {r['name']:>8} ({r['wav']}Å): dr/dz = {r['slope']:+.4f}")
print(f"  Mean dr/dz = {np.mean([r['slope'] for r in nlr_slopes]):+.4f}")

# ============================================================
# APPROACH 3: Let data DEFINE susceptibility
# Instead of using our preset P values, derive P from the data
# P_empirical = -slope (lines that decay faster are more susceptible)
# ============================================================
print(f"\n{'='*70}")
print("APPROACH 3: DATA-DERIVED SUSCEPTIBILITY")
print("Let the data define which lines are most vulnerable")
print("='*70")

# Normalize slopes to [0, 1] range as empirical susceptibility
slope_arr = np.array([r['slope'] for r in rates])
# Most negative slope = most susceptible = P=1
# Most positive slope = least susceptible = P=0
P_emp = (slope_arr.max() - slope_arr) / (slope_arr.max() - slope_arr.min())

print(f"\n{'Line':>8} {'type':>4} {'λ':>6} {'dr/dz':>8} {'P_empirical':>12} {'P_preset':>10}")
print("-" * 56)
for r, p_e in zip(sorted(rates, key=lambda x: x['slope']), 
                   P_emp[np.argsort(slope_arr)]):
    print(f"{r['name']:>8} {r['type']:>4} {r['wav']:>5}Å {r['slope']:>+8.4f} {p_e:>12.3f} {r['P']:>10.1f}")

# Does empirical P correlate with wavelength?
r_emp_wav, p_emp_wav = spearmanr(wavs, P_emp)
print(f"\nEmpirical P vs wavelength: r = {r_emp_wav:+.3f} (p = {p_emp_wav:.3e})")

# ============================================================
# APPROACH 4: TWO-PARAMETER FIT
# For lines that DO decorrelate, fit r = r₀ - β × z
# Then see if β depends on wavelength
# ============================================================
print(f"\n{'='*70}")
print("APPROACH 4: DECORRELATION COEFFICIENT β(λ)")
print("For decaying lines only: r(z) = r₀ - β·z")
print("Does β scale with wavelength?")
print("='*70")

decaying = [r for r in rates if r['slope'] < -0.01]
if len(decaying) >= 3:
    betas = np.array([-r['slope'] for r in decaying])  # positive = decay rate
    wavs_d = np.array([r['wav'] for r in decaying])
    
    print(f"\n{'Line':>8} {'λ':>6} {'β':>8}")
    print("-" * 26)
    for r in sorted(decaying, key=lambda x: x['wav']):
        print(f"{r['name']:>8} {r['wav']:>5}Å {-r['slope']:>8.4f}")
    
    r_bw, p_bw = spearmanr(wavs_d, betas)
    print(f"\nβ vs wavelength: r = {r_bw:+.3f} (p = {p_bw:.3e})")
    
    if len(decaying) >= 3:
        # Fit β = c × λ^α
        log_wav = np.log(wavs_d)
        log_beta = np.log(betas)
        try:
            fit = np.polyfit(log_wav, log_beta, 1)
            alpha_fit = fit[0]
            c_fit = np.exp(fit[1])
            print(f"\nPower law fit: β = {c_fit:.6f} × λ^{alpha_fit:.3f}")
            if alpha_fit < 0:
                print(f"  → Bluer lines decay FASTER (β ∝ λ^{alpha_fit:.1f})")
                print(f"  → For plasma: expect α ≈ -2 (refractive ∝ λ²)")
                print(f"  → Observed: α = {alpha_fit:.2f}")
        except:
            pass

# ============================================================
# APPROACH 5: THE BIG PICTURE TABLE
# Everything we measured, one view
# ============================================================
print(f"\n{'='*70}")
print("THE BIG PICTURE")
print("='*70")

print(f"\n{'Line':>8} {'λ':>6} {'type':>4} {'z_range':>12} {'r_start':>8} {'r_end':>8} {'Δr':>8} {'dr/dz':>8} {'Decays?':>8}")
print("-" * 85)

for name in ['LYA','CIV','CIII','MGII','OII','HBETA','OIII','HALPHA','NII','SII']:
    if name not in all_line_data: continue
    d = all_line_data[name]
    info = line_info[name]
    slope = np.polyfit(d['z'], d['r'], 1)[0]
    decays = "YES" if slope < -0.02 else "no" if slope > 0.02 else "~flat"
    
    print(f"{name:>8} {d['wav']:>5}Å {info['type']:>4} "
          f"[{d['z'][0]:.2f}–{d['z'][-1]:.2f}] "
          f"{d['r'][0]:>8.3f} {d['r'][-1]:>8.3f} {d['r'][-1]-d['r'][0]:>+8.3f} "
          f"{slope:>+8.4f} {decays:>8}")

print(f"""
PATTERN:
  Lines that DECAY (EW-FWHM decorrelate with z):
    → LYA, CIV, HBETA, OIII, SII
  Lines that DON'T (or increase):
    → MGII, OII, CIII, HALPHA, NII
    
  This is NOT simply BLR vs NLR.
  This is NOT simply blue vs red.
  What separates them?
""")

hdu.close()
