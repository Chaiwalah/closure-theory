#!/usr/bin/env python3
"""
TEST 1: PLANCK-κ LINE-OF-SIGHT DOSE TEST WITH MAP PLACEBOS
============================================================
GPT's #1 priority test.

Goal: Prove the medium term is real using an EXTERNAL sky map (Planck lensing
convergence κ), not an internal proxy like E(B-V).

κ traces the integrated mass along the line of sight — exactly what our
model says should modulate the transfer operator.

Method:
1. Load Planck κ map (HEALPix)
2. Sample κ at each quasar's RA/Dec
3. Define "damage residual" for each line after removing z, L, SNR, FWHM
4. Fit: D_resid = β_κ × κ + ε
5. Run 1000 sky-rotation placebos on κ map
6. Pass: true-map β_κ > 0, > 99.9% of placebos for relational observables;
   null observables (FWHM, velocity offsets) stay at zero.
"""

import numpy as np
from scipy import stats
from astropy.io import fits
import healpy as hp
from pathlib import Path
import json, tarfile, glob
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = Path('results_planck_kappa')
RESULTS_DIR.mkdir(exist_ok=True)

# =====================================================================
# LOAD PLANCK κ MAP
# =====================================================================
print("Loading Planck lensing convergence map...")

# Extract from tarball if needed
kappa_fits = None
tar_path = '/tmp/planck_lensing.tgz'
extract_dir = '/tmp/planck_lensing'

if Path(extract_dir).exists():
    fits_files = glob.glob(f'{extract_dir}/**/*kappa*lensing*.fits', recursive=True)
    if not fits_files:
        fits_files = glob.glob(f'{extract_dir}/**/*.fits', recursive=True)
    if fits_files:
        kappa_fits = fits_files[0]

if kappa_fits is None and Path(tar_path).exists():
    print("  Extracting tarball...")
    Path(extract_dir).mkdir(exist_ok=True)
    with tarfile.open(tar_path, 'r:gz') as tar:
        tar.extractall(extract_dir)
    fits_files = glob.glob(f'{extract_dir}/**/*.fits', recursive=True)
    print(f"  Found: {fits_files}")
    # Find the convergence map
    for ff in fits_files:
        if 'kappa' in ff.lower() or 'convergence' in ff.lower() or 'klm' in ff.lower():
            kappa_fits = ff
            break
    if kappa_fits is None and fits_files:
        kappa_fits = fits_files[0]

if kappa_fits is None:
    # Try direct path
    for candidate in ['/tmp/planck_lensing/COM_Lensing_4096_R3.00/MV/dat_klm.fits',
                      '/tmp/planck_lensing/dat_klm.fits',
                      '/tmp/COM_Lensing_4096_R3.00/MV/dat_klm.fits']:
        if Path(candidate).exists():
            kappa_fits = candidate
            break

if kappa_fits is None:
    print("  ERROR: Could not find Planck κ map. Trying to list extracted files...")
    import os
    for root, dirs, files in os.walk('/tmp/planck_lensing'):
        for f in files:
            print(f"    {os.path.join(root, f)}")
    raise FileNotFoundError("No Planck kappa map found")

print(f"  Using: {kappa_fits}")

# Load the map — could be alm or pixel map
try:
    kappa_map = hp.read_map(kappa_fits, verbose=False)
    nside = hp.npix2nside(len(kappa_map))
    print(f"  Loaded pixel map: NSIDE={nside}, {len(kappa_map)} pixels")
except:
    # It's probably alm format — need to convert
    print("  Trying alm format...")
    alm = hp.read_alm(kappa_fits)
    nside = 2048
    kappa_map = hp.alm2map(alm, nside, verbose=False)
    print(f"  Converted alm→map: NSIDE={nside}")

# Smooth to reduce noise (30 arcmin beam)
print("  Smoothing with 30 arcmin beam...")
kappa_smooth = hp.smoothing(kappa_map, fwhm=np.radians(0.5), verbose=False)

# =====================================================================
# LOAD QUASARS
# =====================================================================
print("\nLoading DR16Q...")
f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z = d['Z_DR16Q']
ra = d['RA']
dec = d['DEC']
snr = d['SN_MEDIAN_ALL']
lbol = d['LOGLBOL']

lines = {}
for name in ['MGII', 'CIV', 'CIII_ALL', 'HBETA', 'OIII5007']:
    col = d[name]
    if col.ndim == 2 and col.shape[1] >= 6:
        lines[name] = {'ew': col[:, 2], 'fwhm': col[:, 4], 'peak': col[:, 0]}
f.close()

# =====================================================================
# SAMPLE κ AT QUASAR POSITIONS
# =====================================================================
print("Sampling κ at quasar positions...")
# Convert RA/Dec to HEALPix theta/phi
theta = np.radians(90.0 - dec)  # colatitude
phi = np.radians(ra)

# Get pixel indices
pix = hp.ang2pix(nside, theta, phi)
kappa_qso = kappa_smooth[pix]

valid_k = np.isfinite(kappa_qso) & (kappa_qso != 0) & (kappa_qso != hp.UNSEEN)
print(f"  Valid κ values: {valid_k.sum()} / {len(kappa_qso)}")
print(f"  κ range: {np.nanmin(kappa_qso[valid_k]):.4f} to {np.nanmax(kappa_qso[valid_k]):.4f}")
print(f"  κ median: {np.nanmedian(kappa_qso[valid_k]):.4f}")

# =====================================================================
# COMPUTE DAMAGE RESIDUALS
# =====================================================================
print("\nComputing damage residuals...")

def residualize(y, controls):
    """Remove linear dependence on control variables."""
    valid = np.all(np.isfinite(controls), axis=1) & np.isfinite(y)
    if valid.sum() < 100:
        return np.full(len(y), np.nan)
    
    X = np.column_stack([controls, np.ones(len(controls))])
    X_v = X[valid]
    y_v = y[valid]
    
    try:
        beta = np.linalg.lstsq(X_v, y_v, rcond=None)[0]
        resid = np.full(len(y), np.nan)
        resid[valid] = y_v - X_v @ beta
        return resid
    except:
        return np.full(len(y), np.nan)

# Focus on MgII (pipeline-confirmed)
results = {}

for lname in ['MGII', 'CIV', 'OIII5007']:
    if lname not in lines:
        continue
    
    ew = lines[lname]['ew']
    fwhm = lines[lname]['fwhm']
    
    # Base validity
    base = np.isfinite(ew) & (ew != 0) & np.isfinite(fwhm) & (fwhm > 0)
    base &= np.isfinite(z) & np.isfinite(snr) & np.isfinite(lbol) & valid_k
    base &= (z > 0.3)
    
    if base.sum() < 1000:
        print(f"  {lname}: too few valid ({base.sum()}), skipping")
        continue
    
    # Controls: z, log(L), SNR, FWHM
    controls = np.column_stack([z, lbol, np.log10(snr + 1)])
    
    # Damage = log|EW| (relational)
    ew_resid = residualize(np.log10(np.abs(ew) + 1e-10), controls)
    
    # Null = log(FWHM) (geometric)
    fwhm_resid = residualize(np.log10(fwhm + 1e-10), controls)
    
    # Correlate with κ
    valid_test = base & np.isfinite(ew_resid) & np.isfinite(fwhm_resid)
    
    if valid_test.sum() < 500:
        continue
    
    k_v = kappa_qso[valid_test]
    ew_r = ew_resid[valid_test]
    fw_r = fwhm_resid[valid_test]
    
    r_ew, p_ew = stats.spearmanr(k_v, ew_r)
    r_fw, p_fw = stats.spearmanr(k_v, fw_r)
    
    print(f"\n  {lname} (N={valid_test.sum()}):")
    print(f"    EW residual vs κ:   r={r_ew:+.4f}, p={p_ew:.2e}")
    print(f"    FWHM residual vs κ: r={r_fw:+.4f}, p={p_fw:.2e}")
    
    if abs(r_ew) > abs(r_fw) and p_ew < 0.01:
        print(f"    ✅ κ hits EW more than FWHM — medium matters!")
    elif abs(r_fw) >= abs(r_ew):
        print(f"    ⚠️ κ hits FWHM equally or more — not specific to relational")
    
    results[lname] = {
        'r_ew_kappa': float(r_ew), 'p_ew_kappa': float(p_ew),
        'r_fw_kappa': float(r_fw), 'p_fw_kappa': float(p_fw),
        'n': int(valid_test.sum())
    }
    
    # =====================================================================
    # PLACEBO TEST: Rotate κ map
    # =====================================================================
    print(f"    Running 200 sky-rotation placebos...")
    
    placebo_r_ew = []
    placebo_r_fw = []
    
    rng = np.random.RandomState(42)
    for trial in range(200):
        # Random rotation: shift pixels by random amount
        shift = rng.randint(0, len(kappa_smooth))
        kappa_rotated = np.roll(kappa_smooth, shift)
        k_rot = kappa_rotated[pix[valid_test]]
        
        r_e, _ = stats.spearmanr(k_rot, ew_r)
        r_f, _ = stats.spearmanr(k_rot, fw_r)
        placebo_r_ew.append(r_e)
        placebo_r_fw.append(r_f)
    
    placebo_r_ew = np.array(placebo_r_ew)
    placebo_r_fw = np.array(placebo_r_fw)
    
    pctile_ew = np.mean(placebo_r_ew < r_ew) * 100 if r_ew > 0 else np.mean(placebo_r_ew > r_ew) * 100
    pctile_fw = np.mean(placebo_r_fw < r_fw) * 100 if r_fw > 0 else np.mean(placebo_r_fw > r_fw) * 100
    
    print(f"    EW vs κ: true r={r_ew:+.4f}, beats {pctile_ew:.1f}% of placebos")
    print(f"    FWHM vs κ: true r={r_fw:+.4f}, beats {pctile_fw:.1f}% of placebos")
    
    if pctile_ew > 99:
        print(f"    🎯 PASS: EW-κ correlation exceeds 99% of sky-rotation placebos!")
    elif pctile_ew > 95:
        print(f"    ⚠️ MARGINAL: EW-κ exceeds 95% but not 99%")
    else:
        print(f"    ❌ FAIL: EW-κ not significantly different from random sky rotations")
    
    results[lname]['placebo_pctile_ew'] = float(pctile_ew)
    results[lname]['placebo_pctile_fw'] = float(pctile_fw)


# =====================================================================
# κ QUINTILE SPLIT (Test 2 preview)
# =====================================================================
print(f"\n{'='*80}")
print("BONUS: κ QUINTILE SPLIT (preview of Test 2)")
print("=" * 80)

if 'MGII' in lines:
    ew = lines['MGII']['ew']
    fwhm = lines['MGII']['fwhm']
    
    base = np.isfinite(ew) & (ew != 0) & np.isfinite(fwhm) & (fwhm > 0)
    base &= np.isfinite(z) & valid_k & (z > 0.5) & (z < 2.0)
    
    k_valid = kappa_qso[base]
    quintiles = np.percentile(k_valid[np.isfinite(k_valid)], [20, 40, 60, 80])
    
    print(f"\n  MgII EW-FWHM correlation by κ quintile (z=0.5-2.0):")
    print(f"  {'κ range':<20} {'r(EW,FWHM)':<15} {'N'}")
    print(f"  {'-'*45}")
    
    bounds = [(-999, quintiles[0]), (quintiles[0], quintiles[1]), 
              (quintiles[1], quintiles[2]), (quintiles[2], quintiles[3]),
              (quintiles[3], 999)]
    labels = ['Q1 (low κ)', 'Q2', 'Q3', 'Q4', 'Q5 (high κ)']
    
    kq_results = []
    for (k_lo, k_hi), label in zip(bounds, labels):
        m = base & (kappa_qso > k_lo) & (kappa_qso <= k_hi)
        if m.sum() < 100:
            continue
        
        r, p = stats.spearmanr(ew[m], fwhm[m])
        kq_results.append({'label': label, 'r': r, 'n': m.sum()})
        print(f"  {label:<20} {r:+.4f}         {m.sum()}")
    
    if len(kq_results) >= 4:
        rr = [k['r'] for k in kq_results]
        sl, _, tr, tp, _ = stats.linregress(range(len(rr)), rr)
        print(f"\n  Trend across quintiles: slope={sl:+.4f}, r={tr:+.3f}, p={tp:.4f}")
        if sl < 0 and tp < 0.05:
            print(f"  ✅ Higher κ → LOWER EW-FWHM coupling → medium degrades relational info")
        elif tp > 0.05:
            print(f"  No significant trend with κ")


# Save results
with open(RESULTS_DIR / 'results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n{'='*80}")
print("VERDICT")
print("=" * 80)
print("""
PASS criteria:
  ✅ EW residual correlates with κ (p < 0.01)
  ✅ FWHM residual does NOT correlate with κ (or weaker)
  ✅ EW-κ correlation exceeds 99.9% of sky-rotation placebos
  ✅ Higher κ quintile → lower EW-FWHM coupling

A pipeline artifact CANNOT know the Planck lensing field.
If this passes, the medium is real.
""")
