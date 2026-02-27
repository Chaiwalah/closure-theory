#!/usr/bin/env python3
"""
ROTATION/MOTION DIAGNOSTIC — Looking for spin & motion signatures
==================================================================
The universe has no up/down. But it might have spin.
Test every reference frame we can construct:
1. CMB rest frame (our motion vector)
2. Galactic spin axis
3. Local Group velocity vector
4. Quasar inclination (face-on vs edge-on)
5. Headwind vs tailwind (are we plowing into it?)
6. Full spherical harmonic decomposition of damage
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy import stats
import healpy as hp
import os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_rotation'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("ROTATION/MOTION DIAGNOSTIC")
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

# Damage metric
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

# Source model residuals (quick linear to save time)
from numpy.linalg import lstsq
X_model = np.column_stack([Z[valid], HB_FWHM[valid], LEDD[valid], LBOL[valid], 
                            MBH[valid], SNR_arr[valid], np.ones(np.sum(valid))])
coeffs, _, _, _ = lstsq(X_model, damage[valid], rcond=None)
residuals = np.full(N, np.nan)
residuals[valid] = damage[valid] - X_model @ coeffs
print(f"  Residual std: {np.nanstd(residuals):.4f}")

# Coordinates
coords = SkyCoord(ra=RA*u.deg, dec=DEC*u.deg, frame='icrs')
gal = coords.galactic
GAL_L = gal.l.deg
GAL_B = gal.b.deg

# ============================================================
# KEY REFERENCE DIRECTIONS
# ============================================================
# CMB dipole: direction we're moving (l=264, b=48)
cmb_dir = SkyCoord(l=264*u.deg, b=48*u.deg, frame='galactic')
# Great Attractor (l=307, b=9)
ga_dir = SkyCoord(l=307*u.deg, b=9*u.deg, frame='galactic')
# Milky Way spin axis (North Galactic Pole: l=0, b=90)
mw_spin = SkyCoord(l=0*u.deg, b=90*u.deg, frame='galactic')
# Virgo cluster direction (l=284, b=74)
virgo_dir = SkyCoord(l=284*u.deg, b=74*u.deg, frame='galactic')
# Shapley supercluster
shapley_dir = SkyCoord(l=312*u.deg, b=31*u.deg, frame='galactic')
# CMB Cold Spot (l=210, b=-57)
cold_spot = SkyCoord(l=210*u.deg, b=-57*u.deg, frame='galactic')

# Angular separation from each reference direction
angle_cmb = coords.separation(cmb_dir).deg
angle_ga = coords.separation(ga_dir).deg
angle_spin = coords.separation(mw_spin).deg  # = 90 - |b| roughly
angle_virgo = coords.separation(virgo_dir).deg
angle_shapley = coords.separation(shapley_dir).deg
angle_cold = coords.separation(cold_spot).deg

# ============================================================
# TEST 1: HEADWIND vs TAILWIND (CMB dipole)
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: Headwind vs Tailwind (CMB motion)")
print("=" * 70)
print("  Are photons coming FROM our direction of motion more/less damaged?")

bins = [(0, 30, "Head-on (0-30°)"), (30, 60, "Quarter (30-60°)"),
        (60, 90, "Broadside (60-90°)"), (90, 120, "Quarter-tail (90-120°)"),
        (120, 150, "Near-tail (120-150°)"), (150, 180, "Tail (150-180°)")]

cmb_results = []
for lo, hi, label in bins:
    mask = valid & (angle_cmb >= lo) & (angle_cmb < hi)
    if np.sum(mask) < 100:
        print(f"  {label}: too few ({np.sum(mask)})")
        continue
    mean_r = np.mean(residuals[mask])
    mean_d = np.mean(damage[mask])
    cmb_results.append(((lo+hi)/2, mean_r))
    print(f"  {label}: N={np.sum(mask)}, mean_resid={mean_r:+.4f}, mean_damage={mean_d:.4f}, mean_z={np.mean(Z[mask]):.3f}")

rho_cmb, p_cmb = stats.spearmanr(angle_cmb[valid], residuals[valid])
print(f"\n  ρ(angle_from_CMB, residual) = {rho_cmb:+.4f}, p = {p_cmb:.2e}")

# ============================================================
# TEST 2: GALACTIC SPIN AXIS (are we in a centrifuge?)
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: Galactic Spin Axis")
print("=" * 70)
print("  Does damage depend on angle from MW rotation axis?")

spin_bins = [(0, 30), (30, 60), (60, 90), (90, 120), (120, 150), (150, 180)]
for lo, hi in spin_bins:
    mask = valid & (angle_spin >= lo) & (angle_spin < hi)
    if np.sum(mask) < 100: continue
    print(f"  [{lo:3d}-{hi:3d}°]: N={np.sum(mask)}, resid={np.mean(residuals[mask]):+.4f}, damage={np.mean(damage[mask]):.4f}")

rho_spin, p_spin = stats.spearmanr(angle_spin[valid], residuals[valid])
print(f"\n  ρ(angle_from_spin, residual) = {rho_spin:+.4f}, p = {p_spin:.2e}")

# ============================================================
# TEST 3: GREAT ATTRACTOR ALIGNMENT
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: Great Attractor direction")
print("=" * 70)

for lo, hi in spin_bins:
    mask = valid & (angle_ga >= lo) & (angle_ga < hi)
    if np.sum(mask) < 100: continue
    print(f"  [{lo:3d}-{hi:3d}°]: N={np.sum(mask)}, resid={np.mean(residuals[mask]):+.4f}, damage={np.mean(damage[mask]):.4f}")

rho_ga, p_ga = stats.spearmanr(angle_ga[valid], residuals[valid])
print(f"\n  ρ(angle_from_GA, residual) = {rho_ga:+.4f}, p = {p_ga:.2e}")

# ============================================================
# TEST 4: ALL REFERENCE DIRECTIONS — COMPARISON
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: Reference direction comparison")
print("=" * 70)

for name, angle in [("CMB dipole", angle_cmb), ("Great Attractor", angle_ga),
                     ("MW spin axis", angle_spin), ("Virgo cluster", angle_virgo),
                     ("Shapley SC", angle_shapley), ("CMB Cold Spot", angle_cold)]:
    rho, p = stats.spearmanr(angle[valid], residuals[valid])
    # Also test with raw damage
    rho_raw, p_raw = stats.spearmanr(angle[valid], damage[valid])
    print(f"  {name:20s}: ρ(resid) = {rho:+.4f} p={p:.2e} | ρ(raw) = {rho_raw:+.4f} p={p_raw:.2e}")

# ============================================================
# TEST 5: QUASAR INCLINATION (face-on vs edge-on)
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: Quasar inclination proxy (FWHM as viewing angle)")
print("=" * 70)
print("  Broad FWHM = edge-on (seeing the disk). Narrow = face-on (looking down barrel)")
print("  If operator cares about spin/rotation of SOURCE:")

# Does damage pattern change with FWHM in a way that depends on sky position?
# Cross-term: FWHM × angle_from_spin axis
cross_spin_fwhm = HB_FWHM * np.cos(np.radians(angle_spin))
rho_xsf, p_xsf = stats.spearmanr(cross_spin_fwhm[valid], residuals[valid])
print(f"  ρ(FWHM × cos(angle_spin), residual) = {rho_xsf:+.4f}, p = {p_xsf:.2e}")

cross_cmb_fwhm = HB_FWHM * np.cos(np.radians(angle_cmb))
rho_xcf, p_xcf = stats.spearmanr(cross_cmb_fwhm[valid], residuals[valid])
print(f"  ρ(FWHM × cos(angle_CMB), residual) = {rho_xcf:+.4f}, p = {p_xcf:.2e}")

# ============================================================
# TEST 6: SPHERICAL HARMONIC DECOMPOSITION
# ============================================================
print("\n" + "=" * 70)
print("TEST 6: Spherical harmonic power spectrum of damage")
print("=" * 70)
print("  If operator has rotational structure, specific l-modes will have excess power")

# Build HEALPix map of residuals
NSIDE = 32
theta_hp = np.radians(90 - DEC)
phi_hp = np.radians(RA)
pix_hp = hp.ang2pix(NSIDE, theta_hp, phi_hp)

# Weighted mean per pixel
resid_map = np.full(hp.nside2npix(NSIDE), hp.UNSEEN)
count_map = np.zeros(hp.nside2npix(NSIDE))

for i in range(N):
    if np.isfinite(residuals[i]):
        p = pix_hp[i]
        if resid_map[p] == hp.UNSEEN:
            resid_map[p] = 0
        resid_map[p] += residuals[i]
        count_map[p] += 1

filled = count_map > 0
resid_map[filled] /= count_map[filled]

# Power spectrum
cl = hp.anafast(resid_map, lmax=30)
ell = np.arange(len(cl))

print(f"  Filled pixels: {np.sum(filled)} / {hp.nside2npix(NSIDE)}")
print(f"\n  Power by multipole (low l = large-scale structure):")
for l in range(0, min(16, len(cl))):
    print(f"  l={l:2d}: C_l = {cl[l]:.6f}{'  ← MONOPOLE' if l==0 else '  ← DIPOLE' if l==1 else '  ← QUADRUPOLE' if l==2 else ''}")

# Ratio of low-l to high-l power
if len(cl) > 10:
    low_l = np.mean(cl[1:4])  # dipole + quadrupole + octupole
    high_l = np.mean(cl[10:20])
    print(f"\n  Low-l power (1-3): {low_l:.6f}")
    print(f"  High-l power (10-19): {high_l:.6f}")
    print(f"  Ratio (low/high): {low_l/high_l:.2f}")
    print(f"  → {'EXCESS large-scale power' if low_l/high_l > 5 else 'Normal power distribution'}")

# ============================================================
# TEST 7: THE SPIN TEST — Does damage care about angular momentum?
# ============================================================
print("\n" + "=" * 70)
print("TEST 7: Angular momentum signature")
print("=" * 70)
print("  If there's a rotational component, damage should vary with")
print("  AZIMUTHAL angle around some axis, not just polar angle FROM axis")

# For each reference axis, compute azimuthal angle and test
for name, ref_dir in [("MW spin", mw_spin), ("CMB motion", cmb_dir)]:
    # Project each quasar onto plane perpendicular to reference axis
    ref_vec = np.array(ref_dir.galactic.cartesian.xyz.value)
    q_vecs = np.array(gal.cartesian.xyz.value).T  # (N, 3)
    
    # Azimuthal angle in plane perpendicular to ref
    # Find two orthogonal vectors in the plane
    if abs(ref_vec[2]) < 0.9:
        perp1 = np.cross(ref_vec, [0, 0, 1])
    else:
        perp1 = np.cross(ref_vec, [1, 0, 0])
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(ref_vec, perp1)
    perp2 /= np.linalg.norm(perp2)
    
    # Project quasars
    proj1 = q_vecs @ perp1
    proj2 = q_vecs @ perp2
    azimuth = np.degrees(np.arctan2(proj2, proj1))  # -180 to 180
    
    # Bin by azimuth
    print(f"\n  Reference: {name}")
    az_bins = np.linspace(-180, 180, 9)  # 8 sectors of 45°
    az_resids = []
    az_mids = []
    for i in range(8):
        mask = valid & (azimuth >= az_bins[i]) & (azimuth < az_bins[i+1])
        if np.sum(mask) < 100: continue
        mr = np.mean(residuals[mask])
        az_resids.append(mr)
        az_mids.append((az_bins[i] + az_bins[i+1]) / 2)
        print(f"    [{az_bins[i]:+4.0f}° to {az_bins[i+1]:+4.0f}°]: resid = {mr:+.4f}, N = {np.sum(mask)}")
    
    if len(az_resids) >= 4:
        # Test for sinusoidal pattern (rotation signature)
        az_arr = np.array(az_mids)
        res_arr = np.array(az_resids)
        # Fit: resid = A*sin(az + phi)
        sin_comp = np.sin(np.radians(az_arr))
        cos_comp = np.cos(np.radians(az_arr))
        X_fit = np.column_stack([sin_comp, cos_comp, np.ones(len(az_arr))])
        c, _, _, _ = lstsq(X_fit, res_arr, rcond=None)
        amp = np.sqrt(c[0]**2 + c[1]**2)
        phase = np.degrees(np.arctan2(c[1], c[0]))
        pred = X_fit @ c
        r2 = 1 - np.sum((res_arr - pred)**2) / np.sum((res_arr - np.mean(res_arr))**2)
        print(f"    Sinusoidal fit: amplitude = {amp:.5f}, phase = {phase:.1f}°, R² = {r2:.3f}")
        print(f"    → {'ROTATION SIGNAL' if r2 > 0.5 and amp > 0.005 else 'No rotation signal'}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY — OPERATOR ROTATIONAL/MOTION PROFILE")
print("=" * 70)
print(f"  CMB motion:     ρ = {rho_cmb:+.4f}, p = {p_cmb:.2e}")
print(f"  MW spin:        ρ = {rho_spin:+.4f}, p = {p_spin:.2e}")
print(f"  Great Attractor: ρ = {rho_ga:+.4f}, p = {p_ga:.2e}")
print(f"  Cross FWHM×spin: ρ = {rho_xsf:+.4f}, p = {p_xsf:.2e}")
print(f"  Cross FWHM×CMB:  ρ = {rho_xcf:+.4f}, p = {p_xcf:.2e}")

print("\n✅ Done.")
