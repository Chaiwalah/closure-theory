#!/usr/bin/env python3
"""
OPERATOR SPATIAL PROFILE — No narrative. Just profiling.
=========================================================
Strip source properties + redshift from damage.
Whatever's left = the spatial fingerprint (or lack thereof).

1. Build damage model from source + z
2. Compute residuals
3. Map residuals on the sky
4. Look for ANY structure: angular, κ, survey-dependent, hemisphere
5. Profile: what fraction of damage is source vs distance vs spatial?
"""

import numpy as np
from astropy.io import fits
from scipy import stats
import healpy as hp
import os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_operator_profile'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("OPERATOR SPATIAL PROFILE — Strip and search")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

def fix(arr): return arr.astype(arr.dtype.newbyteorder('='))

z = fix(data['Z_DR16Q'])
ra = fix(data['RA'])
dec = fix(data['DEC'])
snr = fix(data['SN_MEDIAN_ALL'])
oiii_flux = fix(data['OIII5007'][:, 4])
hbeta_flux = fix(data['HBETA'][:, 4])
hbeta_fwhm = fix(data['HBETA'][:, 3])
logLbol = fix(data['LOGLBOL'])
logMbh = fix(data['LOGMBH'])
logLedd = fix(data['LOGLEDD_RATIO'])

good = ((z > 0.3) & (z < 2.0) & (snr > 3) &
        (oiii_flux > 0) & (hbeta_flux > 0) &
        np.isfinite(oiii_flux) & np.isfinite(hbeta_flux) &
        np.isfinite(logLbol) & np.isfinite(logMbh) & np.isfinite(logLedd) &
        (logLbol > 0) & (logMbh > 0) &
        np.isfinite(hbeta_fwhm) & (hbeta_fwhm > 0) &
        np.isfinite(ra) & np.isfinite(dec))

idx = np.where(good)[0]
N = len(idx)
Z = z[idx]; RA = ra[idx]; DEC = dec[idx]; SNR_arr = snr[idx]

log_ratio = np.log10(oiii_flux[idx] / hbeta_flux[idx])
HB_FWHM = np.log10(np.maximum(hbeta_fwhm[idx], 1))
LEDD = logLedd[idx]; LBOL = logLbol[idx]; MBH = logMbh[idx]

print(f"  Sample: {N}")

# ============================================================
# STEP 1: Build source+distance damage model
# ============================================================
print("\n--- Step 1: Source+distance model ---")

# Z-bin normalized scatter (our standard damage metric)
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

damage = compute_damage(log_ratio, Z)
valid = np.isfinite(damage)

# Variance decomposition: how much damage is explained by what?
from sklearn.ensemble import GradientBoostingRegressor

# Model 1: z only
X_z = Z[valid].reshape(-1, 1)
y = damage[valid]

gbr_z = GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=42)
gbr_z.fit(X_z, y)
r2_z = gbr_z.score(X_z, y)
print(f"  R²(z only) = {r2_z:.4f}")

# Model 2: source only (FWHM, Ledd, Lbol, Mbh, SNR)
X_src = np.column_stack([HB_FWHM[valid], LEDD[valid], LBOL[valid], MBH[valid], SNR_arr[valid]])
gbr_src = GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=42)
gbr_src.fit(X_src, y)
r2_src = gbr_src.score(X_src, y)
print(f"  R²(source only) = {r2_src:.4f}")

# Model 3: source + z
X_all = np.column_stack([Z[valid], HB_FWHM[valid], LEDD[valid], LBOL[valid], MBH[valid], SNR_arr[valid]])
gbr_all = GradientBoostingRegressor(n_estimators=100, max_depth=4, random_state=42)
gbr_all.fit(X_all, y)
r2_all = gbr_all.score(X_all, y)
print(f"  R²(source + z) = {r2_all:.4f}")

# Model 4: source + z + sky position
X_sky = np.column_stack([Z[valid], HB_FWHM[valid], LEDD[valid], LBOL[valid], MBH[valid], 
                          SNR_arr[valid], RA[valid], DEC[valid]])
gbr_sky = GradientBoostingRegressor(n_estimators=100, max_depth=4, random_state=42)
gbr_sky.fit(X_sky, y)
r2_sky = gbr_sky.score(X_sky, y)
print(f"  R²(source + z + sky) = {r2_sky:.4f}")

print(f"\n  VARIANCE DECOMPOSITION:")
print(f"  Distance (z):        {r2_z*100:.2f}%")
print(f"  Source properties:    {r2_src*100:.2f}%")
print(f"  Source + distance:    {r2_all*100:.2f}%")
print(f"  Source + dist + sky:  {r2_sky*100:.2f}%")
print(f"  Sky contribution:     {(r2_sky - r2_all)*100:.2f}% (marginal)")
print(f"  Unexplained:          {(1 - r2_sky)*100:.2f}%")

# ============================================================
# STEP 2: Residuals after source+z removal
# ============================================================
print("\n--- Step 2: Residuals (damage - model) ---")

residuals = np.full(N, np.nan)
residuals[valid] = y - gbr_all.predict(X_all)

print(f"  Residual std: {np.nanstd(residuals):.4f}")
print(f"  Residual mean: {np.nanmean(residuals):.6f}")

# ============================================================
# STEP 3: Map residuals on sky
# ============================================================
print("\n--- Step 3: Sky map of residuals ---")

NSIDE = 32
theta = np.radians(90 - DEC)
phi = np.radians(RA)
pix = hp.ang2pix(NSIDE, theta, phi)

# Pixel-level residuals
unique_pix = np.unique(pix)
pix_resid = {}
pix_n = {}
for p in unique_pix:
    mask = (pix == p) & np.isfinite(residuals)
    if np.sum(mask) >= 10:
        pix_resid[p] = np.median(residuals[mask])
        pix_n[p] = np.sum(mask)

print(f"  Pixels with ≥10 quasars: {len(pix_resid)}")

resid_values = np.array(list(pix_resid.values()))
print(f"  Pixel residual range: [{np.min(resid_values):.4f}, {np.max(resid_values):.4f}]")
print(f"  Pixel residual std: {np.std(resid_values):.4f}")

# ============================================================
# STEP 4: Look for structure in residuals
# ============================================================
print("\n--- Step 4: Structure search ---")

# 4A: Hemisphere asymmetry
print("\n  4A: Hemisphere asymmetry")
north = np.isfinite(residuals) & (DEC > 0)
south = np.isfinite(residuals) & (DEC < 0)
print(f"  North: mean resid = {np.mean(residuals[north]):+.4f}, N = {np.sum(north)}")
print(f"  South: mean resid = {np.mean(residuals[south]):+.4f}, N = {np.sum(south)}")
t_ns, p_ns = stats.ttest_ind(residuals[north], residuals[south])
print(f"  t = {t_ns:.2f}, p = {p_ns:.2e}")

# East-West
east = np.isfinite(residuals) & (RA > 180)
west = np.isfinite(residuals) & (RA <= 180)
print(f"  East (RA>180): mean resid = {np.mean(residuals[east]):+.4f}, N = {np.sum(east)}")
print(f"  West (RA≤180): mean resid = {np.mean(residuals[west]):+.4f}, N = {np.sum(west)}")
t_ew, p_ew = stats.ttest_ind(residuals[east], residuals[west])
print(f"  t = {t_ew:.2f}, p = {p_ew:.2e}")

# 4B: Galactic plane proximity
print("\n  4B: Galactic latitude dependence")
# Convert to galactic coords
from astropy.coordinates import SkyCoord
from astropy import units as u
coords = SkyCoord(ra=RA*u.deg, dec=DEC*u.deg, frame='icrs')
gal_b = coords.galactic.b.deg

for label, bmin, bmax in [("Low |b| (<30°)", 0, 30), 
                            ("Mid |b| (30-60°)", 30, 60),
                            ("High |b| (>60°)", 60, 90)]:
    mask = np.isfinite(residuals) & (np.abs(gal_b) >= bmin) & (np.abs(gal_b) < bmax)
    if np.sum(mask) > 100:
        print(f"  {label}: mean = {np.mean(residuals[mask]):+.4f}, std = {np.std(residuals[mask]):.4f}, N = {np.sum(mask)}")

from astropy import units as u
rho_gb, p_gb = stats.spearmanr(np.abs(gal_b[np.isfinite(residuals)]), residuals[np.isfinite(residuals)])
print(f"  ρ(|b|, residual) = {rho_gb:+.4f}, p = {p_gb:.2e}")

# 4C: κ correlation with RESIDUALS (not raw damage)
print("\n  4C: κ vs residuals (after source+z removal)")
try:
    klm = hp.read_alm('/root/clawd/COM_Lensing_4096_R3.00/MV/dat_klm.fits')
    kappa_map = hp.alm2map(klm, 512, verbose=False)
    kappa_pix = hp.ang2pix(512, theta, phi)
    kappa = kappa_map[kappa_pix]
    
    v = np.isfinite(residuals)
    rho_k, p_k = stats.spearmanr(kappa[v], residuals[v])
    print(f"  ρ(κ, residual) = {rho_k:+.4f}, p = {p_k:.2e}")
    
    # Quintiles
    k_edges = np.percentile(kappa[v], np.linspace(0, 100, 6))
    for q in range(5):
        qm = v & (kappa >= k_edges[q]) & (kappa < k_edges[q+1])
        if q == 4: qm = v & (kappa >= k_edges[q])
        print(f"  κ Q{q+1}: mean resid = {np.mean(residuals[qm]):+.4f}, N = {np.sum(qm)}")
except Exception as e:
    print(f"  κ failed: {e}")

# 4D: Survey/plate systematic check
print("\n  4D: Is residual structure just SDSS plate systematics?")
# Check if nearby pixels have correlated residuals AFTER z-detrending
# Use Moran's I as spatial autocorrelation measure
pix_list = list(pix_resid.keys())
pix_vals = np.array([pix_resid[p] for p in pix_list])
pix_coords = np.array([hp.pix2vec(NSIDE, p) for p in pix_list])

# Simple neighbor correlation
n_neighbors = 0
neighbor_products = []
for i in range(len(pix_list)):
    neighbors = hp.get_all_neighbours(NSIDE, pix_list[i])
    for nb in neighbors:
        if nb in pix_resid:
            neighbor_products.append(pix_vals[i] * pix_resid[nb])
            n_neighbors += 1

if n_neighbors > 0:
    mean_prod = np.mean(neighbor_products)
    # Compare to random shuffle
    np.random.seed(42)
    n_shuffles = 1000
    shuffle_prods = []
    for _ in range(n_shuffles):
        shuf = np.random.permutation(pix_vals)
        shuf_dict = {p: s for p, s in zip(pix_list, shuf)}
        sp = []
        for i in range(len(pix_list)):
            neighbors = hp.get_all_neighbours(NSIDE, pix_list[i])
            for nb in neighbors:
                if nb in shuf_dict:
                    sp.append(shuf[i] * shuf_dict[nb])
        shuffle_prods.append(np.mean(sp))
    
    shuffle_mean = np.mean(shuffle_prods)
    shuffle_std = np.std(shuffle_prods)
    z_score = (mean_prod - shuffle_mean) / shuffle_std if shuffle_std > 0 else 0
    print(f"  Neighbor correlation: {mean_prod:.6f}")
    print(f"  Random expectation: {shuffle_mean:.6f} ± {shuffle_std:.6f}")
    print(f"  Z-score: {z_score:.2f}")
    print(f"  → {'SPATIAL STRUCTURE DETECTED' if abs(z_score) > 3 else 'No significant spatial structure'}")

# 4E: Residual vs pure distance (comoving) after z removal
print("\n  4E: Distance anomaly — residual vs comoving distance at fixed source")
# If operator is purely z-dependent (homogeneous), residuals should be flat vs z
# If there's a nonlinear distance effect the GBR missed...
rho_zr, p_zr = stats.spearmanr(Z[valid], residuals[valid])
print(f"  ρ(z, residual) = {rho_zr:+.4f}, p = {p_zr:.2e}")
print(f"  (Should be ~0 if model captured z-dependence properly)")

# 4F: Quadrant analysis — most detailed spatial breakdown
print("\n  4F: RA quadrant analysis")
for ra_lo, ra_hi in [(0, 90), (90, 180), (180, 270), (270, 360)]:
    mask = np.isfinite(residuals) & (RA >= ra_lo) & (RA < ra_hi)
    if np.sum(mask) > 100:
        print(f"  RA [{ra_lo:3d}-{ra_hi:3d}°]: mean = {np.mean(residuals[mask]):+.5f}, "
              f"std = {np.std(residuals[mask]):.4f}, N = {np.sum(mask)}")

# ============================================================
# STEP 5: ANOMALY HUNT — What's weird?
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Anomaly hunt")
print("=" * 70)

# 5A: Find the most anomalous sky regions
print("\n  5A: Most anomalous sky pixels")
sorted_pix = sorted(pix_resid.items(), key=lambda x: abs(x[1]), reverse=True)
print("  Top 10 most anomalous pixels:")
for p, val in sorted_pix[:10]:
    t, ph = hp.pix2ang(NSIDE, p)
    ra_p, dec_p = np.degrees(ph), 90 - np.degrees(t)
    print(f"  Pixel {p}: resid = {val:+.4f}, RA = {ra_p:.1f}°, Dec = {dec_p:.1f}°, N = {pix_n[p]}")

# 5B: Is there a "cold spot" / "hot spot" of damage?
print("\n  5B: Contiguous anomaly search")
# Find clusters of same-sign residuals
pos_pix = [p for p, v in pix_resid.items() if v > np.std(resid_values)]
neg_pix = [p for p, v in pix_resid.items() if v < -np.std(resid_values)]
print(f"  Pixels >1σ positive (excess damage): {len(pos_pix)}")
print(f"  Pixels >1σ negative (less damage): {len(neg_pix)}")

# 5C: Dipole test
print("\n  5C: Dipole in residuals")
# Fit residual = A*cos(θ) + B*sin(θ)*cos(φ) + C*sin(θ)*sin(φ) + D
pix_theta = np.array([hp.pix2ang(NSIDE, p)[0] for p in pix_list])
pix_phi = np.array([hp.pix2ang(NSIDE, p)[1] for p in pix_list])

X_dipole = np.column_stack([
    np.cos(pix_theta),
    np.sin(pix_theta) * np.cos(pix_phi),
    np.sin(pix_theta) * np.sin(pix_phi),
    np.ones(len(pix_list))
])

from numpy.linalg import lstsq
coeffs, resid_sum, _, _ = lstsq(X_dipole, pix_vals, rcond=None)
dipole_pred = X_dipole @ coeffs
r2_dipole = 1 - np.sum((pix_vals - dipole_pred)**2) / np.sum((pix_vals - np.mean(pix_vals))**2)

dipole_amp = np.sqrt(coeffs[0]**2 + coeffs[1]**2 + coeffs[2]**2)
# Direction
dipole_theta = np.arccos(coeffs[0] / dipole_amp)
dipole_phi = np.arctan2(coeffs[2], coeffs[1])
print(f"  Dipole amplitude: {dipole_amp:.5f}")
print(f"  Dipole direction: RA = {np.degrees(dipole_phi):.1f}°, Dec = {90 - np.degrees(dipole_theta):.1f}°")
print(f"  Dipole R²: {r2_dipole:.4f}")
print(f"  Monopole (offset): {coeffs[3]:.5f}")

# F-test for dipole significance
n_pix = len(pix_vals)
F_dipole = (r2_dipole / 3) / ((1 - r2_dipole) / (n_pix - 4))
from scipy.stats import f as f_dist
p_dipole = 1 - f_dist.cdf(F_dipole, 3, n_pix - 4)
print(f"  F-test: F = {F_dipole:.2f}, p = {p_dipole:.2e}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("OPERATOR PROFILE — SPATIAL COMPONENT")
print("=" * 70)
print(f"\n  VARIANCE BUDGET:")
print(f"  Source properties:  {r2_src*100:.1f}%")
print(f"  Distance (z):       {r2_z*100:.1f}%")
print(f"  Combined:           {r2_all*100:.1f}%")
print(f"  + Sky position:     {r2_sky*100:.1f}% (+{(r2_sky-r2_all)*100:.2f}%)")
print(f"  Unexplained:        {(1-r2_sky)*100:.1f}%")
print(f"\n  The operator is {r2_all*100:.0f}% source+distance, {(r2_sky-r2_all)*100:.1f}% spatial, {(1-r2_sky)*100:.0f}% unknown")

print("\n✅ Done.")
