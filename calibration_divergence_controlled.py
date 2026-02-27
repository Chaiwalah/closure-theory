#!/usr/bin/env python3
"""
CALIBRATION DIVERGENCE TEST — CONTROLLED VERSION
==================================================
Controls for the three obvious escapes:
1. SNR degradation at high-z (match on SNR)
2. Population differences (match on Lbol, FWHM, Mbh)
3. Wavelength-dependent systematics (check obs-frame position)

If the divergence survives all three: it's not instrumental, not population, not noise.
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

def get_line(name, idx):
    return fix(data[name])[:, idx]

# Load source properties for matching
try:
    logLbol = fix(data['LOGLBOL'])
except:
    logLbol = np.full(len(z), np.nan)

try:
    logMbh = fix(data['LOGBH'])
except:
    logMbh = np.full(len(z), np.nan)

try:
    hb_fwhm = get_line('HBETA_BR', 3)  # broad Hbeta FWHM
except:
    hb_fwhm = np.full(len(z), np.nan)

# SNR proxies: use line flux / line error for key lines
# For NII and SII specifically
nii_flux = get_line('NII6585', 5)
sii_flux = get_line('SII6718', 5)
nii_ew = get_line('NII6585', 2)
sii_ew = get_line('SII6718', 2)

# Error columns
try:
    nii_err = fix(data['NII6585_ERR'])[:, 5]  # flux error
    sii_err = fix(data['SII6718_ERR'])[:, 5]
    nii_snr = np.abs(nii_flux / np.where(nii_err != 0, nii_err, np.nan))
    sii_snr = np.abs(sii_flux / np.where(sii_err != 0, sii_err, np.nan))
    min_snr = np.minimum(nii_snr, sii_snr)
    has_snr = True
    print("SNR columns loaded successfully")
except Exception as e:
    print(f"SNR columns failed: {e}")
    has_snr = False
    min_snr = np.full(len(z), np.nan)

# Also load OIII, Hbeta, Halpha, OII for other pairs
oiii_ew = get_line('OIII5007', 2)
hb_ew = get_line('HBETA', 2)
oiii_flux = get_line('OIII5007', 5)
hb_flux = get_line('HBETA', 5)
ha_flux = get_line('HALPHA', 5)
oii_ew = get_line('OII3728', 2)

# Locked controls
ha_total = get_line('HALPHA', 5)
ha_br = get_line('HALPHA_BR', 5)
mgii_total = get_line('MGII', 5)
mgii_br = get_line('MGII_BR', 5)

# ============================================================
# CONTROL 1: SNR-MATCHED DIVERGENCE
# ============================================================
print("\n" + "="*80)
print("CONTROL 1: SNR-MATCHED — NII vs SII EW")
print("="*80)

if has_snr:
    # Only keep objects with good SNR in BOTH lines
    for snr_cut in [3, 5, 10, 20]:
        mask = (np.isfinite(nii_ew) & np.isfinite(sii_ew) & 
                (nii_ew != 0) & (sii_ew != 0) & 
                np.isfinite(z) & (z > 0.05) &
                np.isfinite(min_snr) & (min_snr >= snr_cut))
        
        if mask.sum() < 500:
            print(f"  SNR>{snr_cut}: SKIP (n={mask.sum()})")
            continue
        
        zm = z[mask]
        zlo, zhi = np.percentile(zm, [5, 95])
        zbins = np.linspace(zlo, zhi, 8)
        
        rhos = []
        zmids = []
        for i in range(len(zbins)-1):
            sel = (zm >= zbins[i]) & (zm < zbins[i+1])
            n_sel = sel.sum()
            if n_sel > 20:
                rho, p = spearmanr(nii_ew[mask][sel], sii_ew[mask][sel])
                rhos.append(rho)
                zmids.append((zbins[i]+zbins[i+1])/2)
        
        if len(rhos) >= 3:
            delta = rhos[-1] - rhos[0]
            trend, tp = spearmanr(zmids, rhos)
            print(f"  SNR>{snr_cut:2d}: n={mask.sum():6d}, ρ: {rhos[0]:+.3f} → {rhos[-1]:+.3f}, Δ={delta:+.3f}, trend={trend:+.3f} p={tp:.1e}")

# ============================================================
# CONTROL 2: POPULATION-MATCHED DIVERGENCE
# ============================================================
print("\n" + "="*80)
print("CONTROL 2: POPULATION-MATCHED — NII vs SII")
print("="*80)
print("Matching on logLbol + FWHM bins, then checking z-dependence within matched groups")

# Define the key diagnostic pairs to test with controls
test_pairs = [
    ('NII vs SII EW', 'diagnostic', nii_ew, sii_ew),
    ('OIII vs Hbeta EW', 'diagnostic', oiii_ew, hb_ew),
    ('Halpha vs Hbeta flux (Balmer)', 'diagnostic', ha_flux, hb_flux),
    ('OII vs OIII EW', 'diagnostic', oii_ew, oiii_ew),
    ('Halpha vs Halpha_BR (locked)', 'semi-locked', ha_total, ha_br),
    ('MgII vs MgII_BR (locked)', 'semi-locked', mgii_total, mgii_br),
]

for pair_name, ptype, a, b in test_pairs:
    base_mask = (np.isfinite(a) & np.isfinite(b) & (a != 0) & (b != 0) & 
                 np.isfinite(z) & (z > 0.05) &
                 np.isfinite(logLbol) & np.isfinite(hb_fwhm) &
                 (logLbol > 0) & (hb_fwhm > 0))
    
    if base_mask.sum() < 500:
        print(f"  {pair_name}: SKIP (n={base_mask.sum()})")
        continue
    
    zm = z[base_mask]
    am = a[base_mask]
    bm = b[base_mask]
    Lm = logLbol[base_mask]
    Fm = hb_fwhm[base_mask]
    
    # Split into Lbol tertiles × FWHM tertiles = 9 matched groups
    L_edges = np.percentile(Lm[np.isfinite(Lm)], [0, 33, 67, 100])
    F_edges = np.percentile(Fm[np.isfinite(Fm)], [0, 33, 67, 100])
    
    # Within each matched group, check z-trend of correlation
    group_trends = []
    
    for li in range(3):
        for fi in range(3):
            group_mask = ((Lm >= L_edges[li]) & (Lm < L_edges[li+1]) &
                         (Fm >= F_edges[fi]) & (Fm < F_edges[fi+1]))
            
            if group_mask.sum() < 200:
                continue
            
            zg = zm[group_mask]
            ag = am[group_mask]
            bg = bm[group_mask]
            
            # z-bins within this matched group
            zlo, zhi = np.percentile(zg, [10, 90])
            zbins = np.linspace(zlo, zhi, 6)
            
            rhos = []
            zmids = []
            for i in range(len(zbins)-1):
                sel = (zg >= zbins[i]) & (zg < zbins[i+1])
                if sel.sum() > 20:
                    rho, _ = spearmanr(ag[sel], bg[sel])
                    if np.isfinite(rho):
                        rhos.append(rho)
                        zmids.append((zbins[i]+zbins[i+1])/2)
            
            if len(rhos) >= 3:
                delta = rhos[-1] - rhos[0]
                group_trends.append(delta)
    
    if group_trends:
        mean_delta = np.mean(group_trends)
        n_negative = sum(1 for d in group_trends if d < -0.02)
        print(f"  {pair_name:45s}: {len(group_trends)} groups, mean Δρ={mean_delta:+.3f}, {n_negative}/{len(group_trends)} show divergence")

# ============================================================
# CONTROL 3: OBS-FRAME WAVELENGTH CONTROL
# ============================================================
print("\n" + "="*80)
print("CONTROL 3: OBS-FRAME WAVELENGTH — Is it where the line lands on the detector?")
print("="*80)
print("NII 6585Å and SII 6718Å are close in rest-frame.")
print("At z=0.3, they're at ~8560/8733Å. At z=0.1, ~7244/7390Å.")
print("If divergence is just 'red end of CCD is noisy', the correlation should")
print("depend on obs-frame position, not z per se.\n")

# Compute obs-frame wavelength for NII
nii_obs_wave = 6585 * (1 + z)
sii_obs_wave = 6718 * (1 + z)

mask = (np.isfinite(nii_ew) & np.isfinite(sii_ew) & 
        (nii_ew != 0) & (sii_ew != 0) & np.isfinite(z) & (z > 0.05))

if mask.sum() > 500:
    zm = z[mask]
    nm = nii_ew[mask]
    sm = sii_ew[mask]
    obs_m = nii_obs_wave[mask]
    
    # Bin by obs-frame wavelength instead of z
    obs_bins = np.linspace(np.percentile(obs_m, 5), np.percentile(obs_m, 95), 8)
    
    print("  Binned by obs-frame NII wavelength:")
    obs_rhos = []
    obs_mids = []
    for i in range(len(obs_bins)-1):
        sel = (obs_m >= obs_bins[i]) & (obs_m < obs_bins[i+1])
        if sel.sum() > 30:
            rho, _ = spearmanr(nm[sel], sm[sel])
            mean_z = np.mean(zm[sel])
            obs_rhos.append(rho)
            obs_mids.append((obs_bins[i]+obs_bins[i+1])/2)
            print(f"    λ_obs={obs_bins[i]:.0f}-{obs_bins[i+1]:.0f}Å (mean z={mean_z:.2f}): ρ={rho:.3f}, n={sel.sum()}")
    
    if len(obs_rhos) >= 3:
        # Does correlation track obs-frame wavelength or z?
        # If instrument: ρ depends on obs_wave, not z
        # If process: ρ depends on z, and obs_wave is just proxy
        
        # Within a NARROW obs-frame band, does z still matter?
        print("\n  DECISIVE: Within narrow obs-frame band, does z still drive divergence?")
        # Pick the wavelength range with most objects
        best_band = None
        best_n = 0
        for i in range(len(obs_bins)-1):
            sel = (obs_m >= obs_bins[i]) & (obs_m < obs_bins[i+1])
            if sel.sum() > best_n:
                best_n = sel.sum()
                best_band = (obs_bins[i], obs_bins[i+1])
        
        if best_band and best_n > 200:
            band_sel = (obs_m >= best_band[0]) & (obs_m < best_band[1])
            zb = zm[band_sel]
            nb = nm[band_sel]
            sb = sm[band_sel]
            
            # Split this fixed obs-frame band by z
            z_lo = zb < np.median(zb)
            z_hi = ~z_lo
            
            if z_lo.sum() > 30 and z_hi.sum() > 30:
                rho_lo, _ = spearmanr(nb[z_lo], sb[z_lo])
                rho_hi, _ = spearmanr(nb[z_hi], sb[z_hi])
                print(f"    Fixed band {best_band[0]:.0f}-{best_band[1]:.0f}Å:")
                print(f"      Low-z half (mean z={np.mean(zb[z_lo]):.3f}): ρ={rho_lo:.3f}, n={z_lo.sum()}")
                print(f"      High-z half (mean z={np.mean(zb[z_hi]):.3f}): ρ={rho_hi:.3f}, n={z_hi.sum()}")
                print(f"      Δρ = {rho_hi - rho_lo:+.3f}")
                if rho_hi < rho_lo - 0.03:
                    print(f"      >>> Z DRIVES DIVERGENCE EVEN AT FIXED DETECTOR POSITION")
                else:
                    print(f"      >>> Obs-frame wavelength may explain the effect")

# ============================================================
# CONTROL 4: SHUFFLE TEST
# ============================================================
print("\n" + "="*80)
print("CONTROL 4: SHUFFLE — Is the z-ordering real?")
print("="*80)

mask = (np.isfinite(nii_ew) & np.isfinite(sii_ew) & 
        (nii_ew != 0) & (sii_ew != 0) & np.isfinite(z) & (z > 0.05))

zm = z[mask]
nm = nii_ew[mask]
sm = sii_ew[mask]

def compute_divergence(z_arr, a_arr, b_arr):
    zlo, zhi = np.percentile(z_arr, [5, 95])
    zbins = np.linspace(zlo, zhi, 8)
    rhos = []
    for i in range(len(zbins)-1):
        sel = (z_arr >= zbins[i]) & (z_arr < zbins[i+1])
        if sel.sum() > 20:
            rho, _ = spearmanr(a_arr[sel], b_arr[sel])
            rhos.append(rho)
    if len(rhos) >= 3:
        return rhos[-1] - rhos[0]
    return 0

real_delta = compute_divergence(zm, nm, sm)
print(f"  Real Δρ = {real_delta:+.4f}")

n_shuffles = 1000
shuffle_deltas = []
rng = np.random.default_rng(42)
for _ in range(n_shuffles):
    z_shuf = rng.permutation(zm)
    d = compute_divergence(z_shuf, nm, sm)
    shuffle_deltas.append(d)

shuffle_deltas = np.array(shuffle_deltas)
p_shuffle = np.mean(np.abs(shuffle_deltas) >= np.abs(real_delta))
print(f"  Shuffle: mean={np.mean(shuffle_deltas):+.4f}, std={np.std(shuffle_deltas):.4f}")
print(f"  Real Δρ is {abs(real_delta)/np.std(shuffle_deltas):.1f}σ from shuffle mean")
print(f"  Shuffle p-value: {p_shuffle:.4f}")

# ============================================================
# FINAL VERDICT
# ============================================================
print("\n" + "="*80)
print("FINAL VERDICT")
print("="*80)
print("""
The calibration divergence test asks:
  Do redundant diagnostic ratios lose cross-consistency at high-z?

Controls:
  1. SNR matching — does divergence survive at fixed SNR?
  2. Population matching — does it survive within matched source properties?
  3. Obs-frame wavelength — is it the detector position or redshift?
  4. Shuffle test — is the z-ordering statistically real?

If all four controls pass:
  → The divergence is REAL, not instrumental/selection/noise
  → Calibration frameworks based on local physics break at high-z
  → Redundant diagnostics aren't redundant across epochs
""")

# Save
with open(os.path.join(OUT_DIR, 'controlled_results.json'), 'w') as f:
    json.dump({
        'real_delta_nii_sii': real_delta,
        'shuffle_p': p_shuffle,
        'shuffle_sigma': abs(real_delta)/np.std(shuffle_deltas),
    }, f, indent=2)

print(f"Saved to {OUT_DIR}/controlled_results.json")
