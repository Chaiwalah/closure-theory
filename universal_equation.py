#!/usr/bin/env python3
"""
UNIVERSAL EQUATION SEARCH
=========================
Can ONE equation describe how ALL lines decorrelate with z?
Can we collapse all lines onto a single curve?

Try:
1. r(z) = a - b*z  (linear)
2. r(z) = a * exp(-z/τ)  (exponential decay)
3. r(z) = a / (1 + (z/z0)^n)  (generalized sigmoid)
4. r(z) = a * (1 - z/z_max)^p  (power law)
5. Universal curve: r(z, λ) = f(z * g(λ))  (wavelength-scaled)
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

lines = [
    ('LYA', 'LYA', 1216), ('CIV', 'CIV', 1549), ('CIII', 'CIII_ALL', 1909),
    ('MGII', 'MGII', 2798), ('OII', 'OII3728', 3728), ('HBETA', 'HBETA', 4861),
    ('OIII', 'OIII5007', 5007), ('HALPHA', 'HALPHA', 6563),
    ('NII', 'NII6585', 6585), ('SII', 'SII6718', 6718),
]

# ============================================================
# STEP 1: Extract r(EW,FWHM) vs z for every line
# ============================================================
print("Extracting EW-FWHM correlation vs z for all lines...\n")

all_line_data = {}  # name -> (z_centers, correlations, rest_wav)

for name, col_name, rest_wav in lines:
    try:
        col = data[col_name]
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
    
    # Equal-count bins (25 bins)
    idx = np.argsort(zv)
    n = len(idx)
    nbins = min(25, n // 500)
    bs = n // nbins
    
    zmeds, corrs, ns = [], [], []
    for i in range(nbins):
        sl = idx[i*bs:(i+1)*bs]
        if len(sl) < 100:
            continue
        r, _ = spearmanr(ewv[sl], fwv[sl])
        zmeds.append(np.median(zv[sl]))
        corrs.append(r)
        ns.append(len(sl))
    
    zmeds = np.array(zmeds)
    corrs = np.array(corrs)
    all_line_data[name] = (zmeds, corrs, rest_wav, np.array(ns))
    
    print(f"  {name:>8} ({rest_wav}Å): {nbins} bins, z=[{zmeds[0]:.2f}–{zmeds[-1]:.2f}], "
          f"r=[{corrs[0]:.3f}→{corrs[-1]:.3f}]")

# ============================================================
# STEP 2: Fit candidate equations to EACH line independently
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: FIT CANDIDATE EQUATIONS PER LINE")
print("=" * 70)

def f_linear(z, a, b):
    return a + b * z

def f_exp(z, a, tau, c):
    return a * np.exp(-z/tau) + c

def f_power(z, a, z0, n):
    return a / (1 + (z/z0)**n)

def f_log(z, a, b):
    return a - b * np.log(1 + z)

models = [
    ('linear',  f_linear, [0.3, -0.05], 2),
    ('exp',     f_exp,    [0.3, 2.0, 0.1], 3),
    ('power',   f_power,  [0.4, 1.0, 2.0], 3),
    ('log',     f_log,    [0.4, 0.1], 2),
]

print(f"\n{'Line':>8} {'λ':>6} | {'linear':>12} {'exp':>12} {'power':>12} {'log':>12} | {'best':>8}")
print("-" * 80)

best_fits = {}

for name, (zmeds, corrs, wav, ns) in all_line_data.items():
    results = {}
    for mname, func, p0, nparams in models:
        try:
            popt, _ = curve_fit(func, zmeds, corrs, p0=p0, maxfev=10000)
            pred = func(zmeds, *popt)
            ss_res = np.sum((corrs - pred)**2)
            ss_tot = np.sum((corrs - np.mean(corrs))**2)
            r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
            # BIC for model comparison
            n_pts = len(zmeds)
            bic = n_pts * np.log(ss_res/n_pts + 1e-15) + nparams * np.log(n_pts)
            results[mname] = {'r2': r2, 'bic': bic, 'params': popt}
        except:
            results[mname] = {'r2': -99, 'bic': 9999, 'params': None}
    
    # Print R² for each model
    r2s = {m: results[m]['r2'] for m in ['linear','exp','power','log']}
    best = max(r2s, key=r2s.get)
    
    print(f"{name:>8} {wav:>5}Å | "
          f"{r2s['linear']:>11.3f}  {r2s['exp']:>11.3f}  {r2s['power']:>11.3f}  {r2s['log']:>11.3f}  | "
          f"{best:>8}")
    
    best_fits[name] = {'wav': wav, 'results': results, 'best': best}

# ============================================================
# STEP 3: UNIVERSAL CURVE — Can we collapse all lines?
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: UNIVERSAL CURVE SEARCH")
print("Can r(z, λ) = F(z × λ^α) collapse all lines onto ONE curve?")
print("=" * 70)

# Collect all (z, r, λ) triples
all_z = []
all_r = []
all_wav = []

for name, (zmeds, corrs, wav, ns) in all_line_data.items():
    for i in range(len(zmeds)):
        all_z.append(zmeds[i])
        all_r.append(corrs[i])
        all_wav.append(wav)

all_z = np.array(all_z)
all_r = np.array(all_r)
all_wav = np.array(all_wav)

# Try: r = a - b * z * (λ/λ_ref)^α
# This says decorrelation rate depends on wavelength
lambda_ref = 2798.0  # MgII as reference

def universal_1(params, z, wav, r):
    """r = a - b * z * (λ/λ_ref)^α"""
    a, b, alpha = params
    x = z * (wav / lambda_ref)**alpha
    pred = a - b * x
    return np.sum((r - pred)**2)

def universal_2(params, z, wav, r):
    """r = a * exp(-z * (λ/λ_ref)^α / τ) + c"""
    a, tau, alpha, c = params
    x = z * (wav / lambda_ref)**alpha
    pred = a * np.exp(-x / tau) + c
    return np.sum((r - pred)**2)

def universal_3(params, z, wav, r):
    """r = a / (1 + (z*(λ/λ_ref)^α / z0)^n) + c"""
    a, z0, n, alpha, c = params
    x = z * (wav / lambda_ref)**alpha
    pred = a / (1 + (x/z0)**n) + c
    return np.sum((r - pred)**2)

# Also try: r = a - b * log(1 + z * (λ/λ_ref)^α)
def universal_4(params, z, wav, r):
    """r = a - b * log(1 + z * (λ/λ_ref)^α)"""
    a, b, alpha = params
    x = z * (wav / lambda_ref)**alpha
    pred = a - b * np.log(1 + x)
    return np.sum((r - pred)**2)

# Also try with DM instead of z
# DM ∝ z for small z, DM ∝ z^1.2 roughly
def universal_5(params, z, wav, r):
    """r = a - b * z^γ * (λ/λ_ref)^α  (power in z too)"""
    a, b, gamma, alpha = params
    x = z**gamma * (wav / lambda_ref)**alpha
    pred = a - b * x
    return np.sum((r - pred)**2)

candidates = [
    ('r = a - b·z·(λ/λ₀)^α',          universal_1, [0.35, 0.03, -0.5], 3),
    ('r = a·exp(-z·(λ/λ₀)^α/τ) + c',  universal_2, [0.2, 3.0, -0.5, 0.15], 4),
    ('r = a/(1+(z·(λ/λ₀)^α/z₀)^n)+c', universal_3, [0.3, 1.0, 2.0, -0.5, 0.1], 5),
    ('r = a - b·log(1+z·(λ/λ₀)^α)',    universal_4, [0.4, 0.1, -0.5], 3),
    ('r = a - b·z^γ·(λ/λ₀)^α',         universal_5, [0.35, 0.03, 1.0, -0.5], 4),
]

print(f"\nFitting {len(all_z)} data points across {len(all_line_data)} lines...\n")

best_universal = None
best_ss = 1e30

for desc, func, p0, nparams in candidates:
    try:
        result = minimize(func, p0, args=(all_z, all_wav, all_r), method='Nelder-Mead',
                         options={'maxiter': 50000, 'xatol': 1e-8})
        ss = result.fun
        ss_tot = np.sum((all_r - np.mean(all_r))**2)
        r2 = 1 - ss/ss_tot
        bic = len(all_z) * np.log(ss/len(all_z)) + nparams * np.log(len(all_z))
        
        print(f"  {desc}")
        print(f"    R² = {r2:.4f}  BIC = {bic:.1f}  params = {[f'{p:.4f}' for p in result.x]}")
        
        if ss < best_ss:
            best_ss = ss
            best_universal = (desc, result.x, r2, bic)
    except Exception as e:
        print(f"  {desc}: FAILED ({e})")

# ============================================================
# STEP 4: SHOW THE BEST UNIVERSAL FIT
# ============================================================
if best_universal:
    desc, params, r2, bic = best_universal
    print(f"\n{'='*70}")
    print(f"BEST UNIVERSAL FIT")
    print(f"{'='*70}")
    print(f"\n  Equation: {desc}")
    print(f"  R² = {r2:.4f}")
    print(f"  Parameters: {[f'{p:.4f}' for p in params]}")
    print(f"  λ_ref = {lambda_ref} Å (MgII)")
    
    # Show predicted vs observed for each line
    print(f"\n  Per-line residuals:")
    print(f"  {'Line':>8} {'λ':>6} {'r_obs_first':>12} {'r_pred_first':>12} {'r_obs_last':>12} {'r_pred_last':>12}")
    print(f"  {'-'*60}")
    
    for name, (zmeds, corrs, wav, ns) in all_line_data.items():
        if 'log' in desc:
            x = zmeds * (wav/lambda_ref)**params[2]
            pred = params[0] - params[1] * np.log(1 + x)
        elif 'exp' in desc:
            x = zmeds * (wav/lambda_ref)**params[2]
            pred = params[0] * np.exp(-x/params[1]) + params[3]
        elif 'γ' in desc or 'z^' in desc:
            x = zmeds**params[2] * (wav/lambda_ref)**params[3]
            pred = params[0] - params[1] * x
        elif '/(1+' in desc:
            x = zmeds * (wav/lambda_ref)**params[3]
            pred = params[0] / (1 + (x/params[1])**params[2]) + params[4]
        else:
            x = zmeds * (wav/lambda_ref)**params[2]
            pred = params[0] - params[1] * x
        
        print(f"  {name:>8} {wav:>5}Å {corrs[0]:>12.4f} {pred[0]:>12.4f} "
              f"{corrs[-1]:>12.4f} {pred[-1]:>12.4f}")

# ============================================================
# STEP 5: EQUATION WITH PHYSICAL VARIABLES
# ============================================================
print(f"\n{'='*70}")
print("STEP 5: PHYSICAL INTERPRETATION")
print("='*70")

# The scaling variable z * (λ/λ_ref)^α
# If α < 0: bluer lines (smaller λ) → LARGER effective z → decorrelate FASTER
# If α > 0: redder lines decorrelate faster (opposite to our expectation)

if best_universal:
    desc, params, r2, bic = best_universal
    
    # Find α
    if 'γ' in desc or 'z^' in desc:
        alpha = params[3]
        gamma = params[2]
        print(f"\n  Wavelength exponent α = {alpha:.4f}")
        print(f"  Redshift exponent γ = {gamma:.4f}")
    elif len(params) >= 3:
        alpha = params[2] if len(params) <= 4 else params[3]
        print(f"\n  Wavelength exponent α = {alpha:.4f}")
    
    if alpha < 0:
        print(f"  → α < 0: BLUER lines decorrelate FASTER")
        print(f"  → Consistent with blue-preferential degradation ✓")
        print(f"  → Medium is dispersive (frequency-dependent)")
    elif alpha > 0:
        print(f"  → α > 0: REDDER lines decorrelate faster")
        print(f"  → Unexpected direction — needs investigation")
    
    # Physical meaning of the scaling
    print(f"\n  The universal variable is: ξ = z × (λ/{lambda_ref:.0f})^{alpha:.2f}")
    print(f"  All lines collapse onto r(ξ) — one curve, one medium, one process.")
    
    if alpha < 0:
        # For plasma: refractive effects scale as λ² (∝ ν⁻²)
        # So we'd expect α = 2 for pure plasma refraction
        print(f"\n  For pure plasma refraction: expect α = -2 (∝ λ²)")
        print(f"  Observed α = {alpha:.2f}")
        if -2.5 < alpha < -1.5:
            print(f"  → CONSISTENT with plasma dispersion! ✓")
        elif -1.5 < alpha < 0:
            print(f"  → Weaker than pure plasma — may include geometric dilution")
        else:
            print(f"  → Stronger than pure plasma — additional mechanism?")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
hdu.close()
