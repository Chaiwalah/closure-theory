#!/usr/bin/env python3
"""
SPATIAL DOMAIN TEST — Third Domain: The Trampoline Has Texture
===============================================================
Test A: Sightline-dependent degradation (κ correlation)
  - At fixed z, do high-κ sightlines show more diagnostic damage?
Test C: Angular clustering of degradation
  - Do nearby quasars on the sky share similar damage patterns?
  - Diagnostic damage should cluster; locked should not.
"""

import numpy as np
from astropy.io import fits
from astropy import units as u
from scipy import stats
import healpy as hp
import os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_spatial_domain'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("SPATIAL DOMAIN TEST — The Trampoline Has Texture")
print("=" * 70)

# Load quasars
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

def fix(arr): return arr.astype(arr.dtype.newbyteorder('='))

z = fix(data['Z_DR16Q'])
ra = fix(data['RA'])
dec = fix(data['DEC'])
snr = fix(data['SN_MEDIAN_ALL'])
oiii_flux = fix(data['OIII5007'][:, 4])
oiii_ew = fix(data['OIII5007'][:, 2])
hbeta_flux = fix(data['HBETA'][:, 4])
hbeta_fwhm = fix(data['HBETA'][:, 3])
nii_flux = fix(data['NII6585'][:, 4])
halpha_flux = fix(data['HALPHA'][:, 4])
sii_flux = fix(data['SII6718'][:, 4])
oii_flux = fix(data['OII3728'][:, 4])
# Locked proxy: NII/Hα (both from same region, relatively stable ratio)
# No OIII 4959 in DR16Q prop file
logLbol = fix(data['LOGLBOL'])
logMbh = fix(data['LOGMBH'])
logLedd = fix(data['LOGLEDD_RATIO'])

# Quality cuts
good = ((z > 0.3) & (z < 2.0) & (snr > 3) &
        (oiii_flux > 0) & (hbeta_flux > 0) &
        np.isfinite(oiii_flux) & np.isfinite(hbeta_flux) &
        np.isfinite(logLbol) & np.isfinite(logMbh) & np.isfinite(logLedd) &
        (logLbol > 0) & (logMbh > 0) &
        np.isfinite(hbeta_fwhm) & (hbeta_fwhm > 0) &
        np.isfinite(ra) & np.isfinite(dec))

idx = np.where(good)[0]
N = len(idx)
Z = z[idx]
RA = ra[idx]
DEC = dec[idx]
SNR = snr[idx]

log_diag = np.log10(oiii_flux[idx] / hbeta_flux[idx])  # diagnostic
HB_FWHM = np.log10(np.maximum(hbeta_fwhm[idx], 1))

# Locked-ish ratio: NII/Hα (same region, low diagnostic sensitivity)
nii_f = nii_flux[idx]
ha_f = halpha_flux[idx]
locked_good = (nii_f > 0) & (ha_f > 0) & np.isfinite(nii_f) & np.isfinite(ha_f)
log_locked = np.full(N, np.nan)
log_locked[locked_good] = np.log10(nii_f[locked_good] / ha_f[locked_good])

print(f"  Base sample: {N}")
print(f"  With locked ratio: {np.sum(locked_good)}")

# Z-bin normalized damage
def compute_damage(ratio, redshift, n_bins=20):
    edges = np.percentile(redshift, np.linspace(0, 100, n_bins + 1))
    zb = np.clip(np.digitize(redshift, edges) - 1, 0, n_bins - 1)
    damage = np.full(len(ratio), np.nan)
    for b in range(n_bins):
        m = zb == b
        if np.sum(m) < 10: continue
        med = np.median(ratio[m])
        mad = max(1.4826 * np.median(np.abs(ratio[m] - med)), 1e-6)
        damage[m] = np.abs((ratio[m] - med) / mad)
    return damage

diag_damage = compute_damage(log_diag, Z)
locked_damage = compute_damage(log_locked[locked_good], Z[locked_good])

# ============================================================
# TEST A: SIGHTLINE κ CORRELATION
# ============================================================
print("\n" + "=" * 70)
print("TEST A: Sightline-dependent degradation (Planck κ)")
print("=" * 70)

# Load Planck lensing convergence
try:
    klm_file = '/root/clawd/COM_Lensing_4096_R3.00/MV/dat_klm.fits'
    klm = hp.read_alm(klm_file)
    NSIDE = 512
    kappa_map = hp.alm2map(klm, NSIDE, verbose=False)
    print(f"  Planck κ map loaded: NSIDE={NSIDE}, {len(kappa_map)} pixels")
    
    # Get κ for each quasar
    theta = np.radians(90 - DEC)  # co-latitude
    phi = np.radians(RA)
    pix = hp.ang2pix(NSIDE, theta, phi)
    kappa = kappa_map[pix]
    
    print(f"  κ range: [{np.min(kappa):.4f}, {np.max(kappa):.4f}]")
    print(f"  κ median: {np.median(kappa):.4f}")
    
    # Test A1: κ vs diagnostic damage at fixed z
    print("\n  --- A1: κ-damage correlation in z-bins ---")
    z_edges = np.percentile(Z, [0, 25, 50, 75, 100])
    for i in range(4):
        zbin = (Z >= z_edges[i]) & (Z < z_edges[i+1])
        if i == 3: zbin = (Z >= z_edges[i]) & (Z <= z_edges[i+1])
        valid = zbin & np.isfinite(diag_damage)
        if np.sum(valid) < 100: continue
        rho, p = stats.spearmanr(kappa[valid], diag_damage[valid])
        print(f"  z=[{z_edges[i]:.2f},{z_edges[i+1]:.2f}]: ρ(κ,damage)={rho:+.4f}, p={p:.2e}, N={np.sum(valid)}")
    
    # Overall
    valid_all = np.isfinite(diag_damage)
    rho_all, p_all = stats.spearmanr(kappa[valid_all], diag_damage[valid_all])
    print(f"  Overall: ρ(κ,diag_damage) = {rho_all:+.4f}, p = {p_all:.2e}")
    
    # Test A2: Same but for LOCKED ratio
    print("\n  --- A2: κ vs LOCKED damage (control) ---")
    locked_kappa = kappa[locked_good]
    valid_locked = np.isfinite(locked_damage)
    if np.sum(valid_locked) > 100:
        rho_l, p_l = stats.spearmanr(locked_kappa[valid_locked], locked_damage[valid_locked])
        print(f"  ρ(κ,locked_damage) = {rho_l:+.4f}, p = {p_l:.2e}, N={np.sum(valid_locked)}")
    
    # Test A3: κ quintiles — damage profile
    print("\n  --- A3: Damage by κ quintile (fixed z) ---")
    # Use middle z range for cleanest test
    zmid = (Z > 0.5) & (Z < 1.5) & np.isfinite(diag_damage)
    kappa_mid = kappa[zmid]
    damage_mid = diag_damage[zmid]
    k_edges = np.percentile(kappa_mid, np.linspace(0, 100, 6))
    print(f"  (z=0.5-1.5, N={np.sum(zmid)})")
    q_damages = []
    for q in range(5):
        qmask = (kappa_mid >= k_edges[q]) & (kappa_mid < k_edges[q+1])
        if q == 4: qmask = (kappa_mid >= k_edges[q]) & (kappa_mid <= k_edges[q+1])
        md = np.median(damage_mid[qmask])
        q_damages.append(md)
        print(f"  Q{q+1} (κ={np.median(kappa_mid[qmask]):+.4f}): median damage={md:.4f}, N={np.sum(qmask)}")
    
    rho_q, p_q = stats.spearmanr(range(5), q_damages)
    print(f"  Quintile trend: ρ = {rho_q:+.3f}, p = {p_q:.3f}")
    
    # Test A4: Susceptibility × κ interaction
    print("\n  --- A4: Does κ effect depend on susceptibility? ---")
    s_med = np.median(HB_FWHM)
    for label, smask in [("Low S (narrow FWHM)", HB_FWHM < s_med), 
                          ("High S (broad FWHM)", HB_FWHM >= s_med)]:
        vm = smask & np.isfinite(diag_damage)
        rho_s, p_s = stats.spearmanr(kappa[vm], diag_damage[vm])
        print(f"  {label}: ρ(κ,damage) = {rho_s:+.4f}, p = {p_s:.2e}, N={np.sum(vm)}")

    HAS_KAPPA = True
except Exception as e:
    print(f"  Failed to load Planck κ: {e}")
    HAS_KAPPA = False

# ============================================================
# TEST C: ANGULAR CLUSTERING OF DEGRADATION
# ============================================================
print("\n" + "=" * 70)
print("TEST C: Angular clustering of diagnostic damage")
print("=" * 70)
print("  Do nearby quasars share similar damage? Diagnostic vs locked.")

# Use HEALPix to bin quasars into sky pixels
NSIDE_SKY = 32  # ~1.8 deg pixels
theta_q = np.radians(90 - DEC)
phi_q = np.radians(RA)
pix_q = hp.ang2pix(NSIDE_SKY, theta_q, phi_q)

# For each pixel with enough quasars, compute mean damage
unique_pix = np.unique(pix_q)
print(f"  Sky pixels (NSIDE={NSIDE_SKY}): {len(unique_pix)} occupied")

# Build pixel-level damage maps
pix_diag_damage = {}
pix_locked_damage = {}
pix_z = {}
pix_count = {}

for p in unique_pix:
    mask = (pix_q == p) & np.isfinite(diag_damage)
    if np.sum(mask) < 5: continue
    pix_diag_damage[p] = np.median(diag_damage[mask])
    pix_z[p] = np.median(Z[mask])
    pix_count[p] = np.sum(mask)
    
    # Locked
    lmask = locked_good.copy()
    lmask_full = np.zeros(N, dtype=bool)
    lmask_full[locked_good] = np.isfinite(locked_damage)
    lpix = (pix_q == p) & lmask_full
    if np.sum(lpix) >= 3:
        # Get locked damage for these objects
        # Map back to locked_damage indices
        locked_idx_map = np.where(locked_good)[0]
        in_pix = np.isin(locked_idx_map, np.where(lpix)[0])
        if np.sum(in_pix) >= 3:
            pix_locked_damage[p] = np.median(locked_damage[in_pix])

valid_pix = list(pix_diag_damage.keys())
print(f"  Pixels with ≥5 quasars: {len(valid_pix)}")
print(f"  Pixels with locked data: {len(pix_locked_damage)}")

# Compute angular separation between all pixel pairs and correlation
# Use pixel centers
if len(valid_pix) > 50:
    pix_centers = np.array([hp.pix2ang(NSIDE_SKY, p) for p in valid_pix])
    pix_dd = np.array([pix_diag_damage[p] for p in valid_pix])
    pix_zz = np.array([pix_z[p] for p in valid_pix])
    
    # Z-detrend the pixel damage
    sl, ic, _, _, _ = stats.linregress(pix_zz, pix_dd)
    pix_dd_detrend = pix_dd - (sl * pix_zz + ic)
    
    # Angular separation bins
    n_pairs_per_bin = []
    corr_per_bin = []
    sep_bins = [0, 2, 5, 10, 20, 40, 90, 180]  # degrees
    
    print(f"\n  --- Angular correlation function (diagnostic) ---")
    # Sample pairs (full NxN is too expensive for >1000 pixels)
    n_pix = len(valid_pix)
    
    # Compute all pairwise separations efficiently
    # For manageable computation, subsample if needed
    if n_pix > 2000:
        np.random.seed(42)
        sub = np.random.choice(n_pix, 2000, replace=False)
    else:
        sub = np.arange(n_pix)
    
    n_sub = len(sub)
    print(f"  Using {n_sub} pixels for pair computation")
    
    # Vectorized angular separation
    theta_s = pix_centers[sub, 0]
    phi_s = pix_centers[sub, 1]
    dd_s = pix_dd_detrend[sub]
    
    for b in range(len(sep_bins) - 1):
        lo_sep = np.radians(sep_bins[b])
        hi_sep = np.radians(sep_bins[b+1])
        
        pairs_d1 = []
        pairs_d2 = []
        
        # Sample pairs within this angular bin
        for i in range(min(n_sub, 500)):
            for j in range(i+1, min(n_sub, 500)):
                # Angular separation
                cos_sep = (np.sin(np.pi/2 - theta_s[i]) * np.sin(np.pi/2 - theta_s[j]) +
                          np.cos(np.pi/2 - theta_s[i]) * np.cos(np.pi/2 - theta_s[j]) *
                          np.cos(phi_s[i] - phi_s[j]))
                cos_sep = np.clip(cos_sep, -1, 1)
                sep = np.arccos(cos_sep)
                
                if lo_sep <= sep < hi_sep:
                    pairs_d1.append(dd_s[i])
                    pairs_d2.append(dd_s[j])
        
        if len(pairs_d1) > 10:
            r, p = stats.spearmanr(pairs_d1, pairs_d2)
            print(f"  [{sep_bins[b]:3d}°-{sep_bins[b+1]:3d}°]: ρ = {r:+.4f}, p = {p:.2e}, N_pairs = {len(pairs_d1)}")
            corr_per_bin.append(r)
        else:
            print(f"  [{sep_bins[b]:3d}°-{sep_bins[b+1]:3d}°]: too few pairs ({len(pairs_d1)})")
            corr_per_bin.append(np.nan)
    
    # Now do the same for locked damage
    locked_valid_pix = [p for p in valid_pix if p in pix_locked_damage]
    if len(locked_valid_pix) > 50:
        print(f"\n  --- Angular correlation function (LOCKED — control) ---")
        lpix_centers = np.array([hp.pix2ang(NSIDE_SKY, p) for p in locked_valid_pix])
        lpix_dd = np.array([pix_locked_damage[p] for p in locked_valid_pix])
        lpix_zz = np.array([pix_z[p] for p in locked_valid_pix])
        
        sl_l, ic_l, _, _, _ = stats.linregress(lpix_zz, lpix_dd)
        lpix_dd_detrend = lpix_dd - (sl_l * lpix_zz + ic_l)
        
        n_lpix = len(locked_valid_pix)
        if n_lpix > 500: 
            lsub = np.random.choice(n_lpix, 500, replace=False)
        else:
            lsub = np.arange(n_lpix)
        
        ltheta = lpix_centers[lsub, 0]
        lphi = lpix_centers[lsub, 1]
        ldd = lpix_dd_detrend[lsub]
        
        for b in range(len(sep_bins) - 1):
            lo_sep = np.radians(sep_bins[b])
            hi_sep = np.radians(sep_bins[b+1])
            
            lp1 = []
            lp2 = []
            n_lsub = len(lsub)
            
            for i in range(min(n_lsub, 300)):
                for j in range(i+1, min(n_lsub, 300)):
                    cos_sep = (np.sin(np.pi/2 - ltheta[i]) * np.sin(np.pi/2 - ltheta[j]) +
                              np.cos(np.pi/2 - ltheta[i]) * np.cos(np.pi/2 - ltheta[j]) *
                              np.cos(lphi[i] - lphi[j]))
                    cos_sep = np.clip(cos_sep, -1, 1)
                    sep = np.arccos(cos_sep)
                    
                    if lo_sep <= sep < hi_sep:
                        lp1.append(ldd[i])
                        lp2.append(ldd[j])
            
            if len(lp1) > 10:
                r, p = stats.spearmanr(lp1, lp2)
                print(f"  [{sep_bins[b]:3d}°-{sep_bins[b+1]:3d}°]: ρ = {r:+.4f}, p = {p:.2e}, N_pairs = {len(lp1)}")
            else:
                print(f"  [{sep_bins[b]:3d}°-{sep_bins[b+1]:3d}°]: too few pairs ({len(lp1)})")

# ============================================================
# TEST C2: DIFFERENTIAL — Diagnostic minus Locked angular clustering
# ============================================================
print("\n" + "=" * 70)
print("TEST C2: Pixel-level diagnostic vs locked damage correlation")
print("=" * 70)

both_pix = [p for p in valid_pix if p in pix_locked_damage]
if len(both_pix) > 30:
    dd_both = np.array([pix_diag_damage[p] for p in both_pix])
    ld_both = np.array([pix_locked_damage[p] for p in both_pix])
    zz_both = np.array([pix_z[p] for p in both_pix])
    
    # Are diagnostic and locked spatially correlated?
    r_dl, p_dl = stats.spearmanr(dd_both, ld_both)
    print(f"  ρ(diag_damage, locked_damage) by pixel = {r_dl:+.4f}, p = {p_dl:.2e}, N = {len(both_pix)}")
    
    # Variance comparison
    print(f"  Diagnostic damage pixel variance: {np.var(dd_both):.6f}")
    print(f"  Locked damage pixel variance: {np.var(ld_both):.6f}")
    ratio_var = np.var(dd_both) / max(np.var(ld_both), 1e-10)
    print(f"  Ratio (diag/locked): {ratio_var:.2f}x")
    print(f"  If operator has spatial structure: diagnostic variance >> locked variance")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
if HAS_KAPPA:
    print(f"  A) κ-damage correlation: ρ = {rho_all:+.4f}, p = {p_all:.2e}")
    print(f"     → More matter along sightline = {'MORE' if rho_all > 0 else 'LESS'} diagnostic damage")
print(f"  C) Angular clustering: see correlation functions above")
print(f"     → If diagnostic clusters but locked doesn't = operator has spatial structure")

print("\n✅ Done. Results in", OUTDIR)
