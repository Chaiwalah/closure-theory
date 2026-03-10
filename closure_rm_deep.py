#!/usr/bin/env python3
"""
RM DEEP ANALYSIS — Squeeze everything from the 2,161 direct matches

These are radio-loud DR16Q quasars with their OWN measured RM from NVSS.
Their RM includes: source intrinsic + IGM path + MW foreground.

Key result from first pass: Δ = -0.369 at z > 1 (high RM vs low RM).
Now: break it down. What's driving this? Can we subtract MW contribution?
Does the RESIDUAL (extragalactic) RM predict degradation?

Also: look at ALL our previous findings through the RM lens.
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
from scipy.optimize import curve_fit
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_rm_deep', exist_ok=True)

print("=" * 70)
print("RM DEEP — Mining the 2,161 direct matches")
print("=" * 70)

# Load DR16Q
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
ciii_ew = d['CIII_BR'][:, 2]
civ_ew = d['CIV'][:, 2]
ra_q = d['RA']
dec_q = d['DEC']

# Load NVSS RM catalog
rm_cat = Table.read('data/nvss_rm_taylor2009.fits')
coords_rm = SkyCoord(
    ra=[str(r) for r in rm_cat['RAJ2000']], 
    dec=[str(d) for d in rm_cat['DEJ2000']], 
    unit=(u.hourangle, u.degree), frame='icrs'
)
rm_val = np.array(rm_cat['RM'], dtype=float)
rm_err = np.array(rm_cat['e_RM'], dtype=float)

# Build quasar coordinates
coords_q = SkyCoord(ra=ra_q*u.degree, dec=dec_q*u.degree, frame='icrs')

# Cross match
idx_rm, sep, _ = coords_q.match_to_catalog_sky(coords_rm)
matched = sep < 30 * u.arcsec

print(f"Direct matches: {matched.sum()}")

# Extract matched sample
z_m = z_all[matched]
rm_m = rm_val[idx_rm[matched]]
rm_e = rm_err[idx_rm[matched]]
mg_m = mg_ew[matched]
hb_m = hb_ew[matched]
ciii_m = ciii_ew[matched]
civ_m = civ_ew[matched]

# Galactic coordinates for MW RM subtraction
gal_lat_m = coords_q[matched].galactic.b.degree
gal_lon_m = coords_q[matched].galactic.l.degree

abs_rm = np.abs(rm_m)

print(f"z range: {z_m.min():.2f} - {z_m.max():.2f}")
print(f"|RM| range: {abs_rm.min():.1f} - {abs_rm.max():.1f} rad/m²")
print(f"|RM| median: {np.median(abs_rm):.1f} rad/m²")

# ============================================================
# MW RM SUBTRACTION (approximate)
# ============================================================
print(f"\n{'='*70}")
print("MW RM SUBTRACTION")
print("MW RM ∝ 1/sin(|b|) at low latitudes. Subtract trend to get residual RM.")
print(f"{'='*70}")

# The MW contribution dominates at low |b|. At high |b|, RM is mostly extragalactic.
# Simple model: RM = RM_extragalactic + RM_MW * csc(|b|)
# We'll use the median RM at high |b| as baseline

high_b = np.abs(gal_lat_m) > 60
low_b = np.abs(gal_lat_m) < 30

median_rm_high_b = np.median(abs_rm[high_b])
median_rm_low_b = np.median(abs_rm[low_b])
print(f"Median |RM| at |b| > 60°: {median_rm_high_b:.1f} rad/m² (mostly extragalactic)")
print(f"Median |RM| at |b| < 30°: {median_rm_low_b:.1f} rad/m² (MW-dominated)")

# Simple residual: subtract latitude-dependent MW model
# Use running median of |RM| vs |b| as MW model
abs_b = np.abs(gal_lat_m)
b_bins = np.arange(20, 85, 5)
mw_model = np.zeros_like(abs_rm)

for i in range(len(b_bins)-1):
    mask = (abs_b >= b_bins[i]) & (abs_b < b_bins[i+1])
    if mask.sum() > 10:
        med = np.median(abs_rm[mask])
        mw_model[mask] = med

# For objects outside the binned range
mw_model[abs_b < 20] = median_rm_low_b
mw_model[abs_b >= 80] = median_rm_high_b

# Residual |RM| (excess above MW median for that latitude)
rm_residual = abs_rm - mw_model
print(f"Residual |RM| range: {rm_residual.min():.1f} to {rm_residual.max():.1f}")
print(f"Residual |RM| median: {np.median(rm_residual):.1f}")

# ============================================================
# TEST with RESIDUAL RM (MW-subtracted)
# ============================================================
print(f"\n{'='*70}")
print("TEST: Does RESIDUAL |RM| (extragalactic excess) predict degradation?")
print(f"{'='*70}")

rm_res_med = np.median(rm_residual)

def corr_pair(v1, v2, mask):
    d1, d2 = v1[mask], v2[mask]
    valid = np.isfinite(d1) & np.isfinite(d2) & (d1 != 0) & (d2 != 0)
    if valid.sum() > 15:
        r, p = stats.spearmanr(d1[valid], d2[valid])
        return r, int(valid.sum()), p
    return np.nan, 0, np.nan

z_bins = [(0.5, 0.8), (0.8, 1.0), (1.0, 1.3), (1.3, 2.0)]

print(f"\n  MgII↔Hβ: Low residual RM vs High residual RM")
print(f"  {'z-bin':<15} {'Low res|RM|':>15} {'High res|RM|':>15} {'Δ':>8}")
print(f"  {'-'*15} {'-'*15} {'-'*15} {'-'*8}")

for zlo, zhi in z_bins:
    zmask = (z_m >= zlo) & (z_m < zhi)
    lo = zmask & (rm_residual <= rm_res_med)
    hi = zmask & (rm_residual > rm_res_med)
    
    r_lo, n_lo, _ = corr_pair(mg_m, hb_m, lo)
    r_hi, n_hi, _ = corr_pair(mg_m, hb_m, hi)
    
    row = f"z=[{zlo},{zhi})"
    lo_str = f"{r_lo:+.3f} ({n_lo})" if not np.isnan(r_lo) else "N/A"
    hi_str = f"{r_hi:+.3f} ({n_hi})" if not np.isnan(r_hi) else "N/A"
    
    if not np.isnan(r_lo) and not np.isnan(r_hi):
        delta = r_hi - r_lo
        print(f"  {row:<15} {lo_str:>15} {hi_str:>15} {delta:>+8.3f}")
    else:
        print(f"  {row:<15} {lo_str:>15} {hi_str:>15} {'N/A':>8}")

# ============================================================
# RM × Z interaction — the key test
# ============================================================
print(f"\n{'='*70}")
print("RM × Z INTERACTION — Does RM matter MORE at higher z?")
print("If wave equilibration is cumulative, RM effect should GROW with z")
print(f"{'='*70}")

z_fine = [(0.4, 0.6), (0.6, 0.8), (0.8, 0.95), (0.95, 1.1), (1.1, 1.3), (1.3, 1.8)]

deltas = []
z_mids = []

print(f"\n  {'z-range':<12} {'Low |RM|':>10} {'High |RM|':>10} {'Δ':>8} {'N_lo':>6} {'N_hi':>6} {'Visual':>20}")
print(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*8} {'-'*6} {'-'*6} {'-'*20}")

for zlo, zhi in z_fine:
    zmask = (z_m >= zlo) & (z_m < zhi)
    lo = zmask & (abs_rm <= np.median(abs_rm))
    hi = zmask & (abs_rm > np.median(abs_rm))
    
    r_lo, n_lo, _ = corr_pair(mg_m, hb_m, lo)
    r_hi, n_hi, _ = corr_pair(mg_m, hb_m, hi)
    
    row = f"z=[{zlo},{zhi})"
    
    if not np.isnan(r_lo) and not np.isnan(r_hi):
        delta = r_hi - r_lo
        deltas.append(delta)
        z_mids.append((zlo + zhi) / 2)
        
        if delta < 0:
            bar = "◄" * min(int(abs(delta) * 30), 20)
        else:
            bar = "►" * min(int(abs(delta) * 30), 20)
        
        print(f"  {row:<12} {r_lo:>+10.3f} {r_hi:>+10.3f} {delta:>+8.3f} {n_lo:>6} {n_hi:>6} {bar}")
    else:
        print(f"  {row:<12} {'N/A':>10} {'N/A':>10} {'N/A':>8}")

if len(deltas) > 2:
    deltas = np.array(deltas)
    z_mids = np.array(z_mids)
    
    # Does the RM effect grow with z?
    r_growth, p_growth = stats.spearmanr(z_mids, deltas)
    print(f"\n  RM effect growth with z: ρ = {r_growth:+.3f} (p = {p_growth:.3f})")
    if r_growth < -0.5 and p_growth < 0.1:
        print(f"  ★ RM EFFECT GROWS NEGATIVELY with z — cumulative wave-medium equilibration!")

# ============================================================
# MULTIPLE LINE PAIRS — Does RM affect optical MORE than UV?
# ============================================================
print(f"\n{'='*70}")
print("WAVELENGTH TEST — Does RM affect optical (Hβ) more than UV (CIII)?")
print("Faraday ∝ λ² → optical should show MORE RM-dependent degradation")
print(f"{'='*70}")

print(f"\n  At z = 0.85-1.1 (transition zone, most constraining):")
zmask_trans = (z_m >= 0.85) & (z_m < 1.1)
lo_rm_trans = zmask_trans & (abs_rm <= np.median(abs_rm))
hi_rm_trans = zmask_trans & (abs_rm > np.median(abs_rm))

pairs = [
    ('MgII↔Hβ (opt)', mg_m, hb_m, 4861),
    ('MgII↔CIII (UV)', mg_m, ciii_m, 1909),
    ('MgII↔CIV (UV)', mg_m, civ_m, 1549),
]

print(f"  {'Pair':<20} {'λ(Å)':>6} {'Low |RM|':>10} {'High |RM|':>10} {'Δ':>8}")
print(f"  {'-'*20} {'-'*6} {'-'*10} {'-'*10} {'-'*8}")

for name, v1, v2, lam in pairs:
    r_lo, n_lo, _ = corr_pair(v1, v2, lo_rm_trans)
    r_hi, n_hi, _ = corr_pair(v1, v2, hi_rm_trans)
    
    if not np.isnan(r_lo) and not np.isnan(r_hi):
        delta = r_hi - r_lo
        print(f"  {name:<20} {lam:>6} {r_lo:>+10.3f} {r_hi:>+10.3f} {delta:>+8.3f}")
    else:
        print(f"  {name:<20} {lam:>6} {'N/A':>10} {'N/A':>10} {'N/A':>8}")

# ============================================================
# DOES RM SHIFT THE SIGMOID?
# ============================================================
print(f"\n{'='*70}")
print("SIGMOID SHIFT — Does high |RM| move the transition to LOWER z?")
print("If medium accelerates equilibration, threshold should shift earlier")
print(f"{'='*70}")

def sigmoid(z, rho0, k, z0):
    return rho0 / (1 + np.exp(k * (z - z0)))

# Compute sliding correlation for low and high RM
z_windows = np.arange(0.5, 1.4, 0.05)
window = 0.15

for label, rm_mask, color in [
    ("Low |RM|", abs_rm <= np.percentile(abs_rm, 33), "blue"),
    ("High |RM|", abs_rm >= np.percentile(abs_rm, 67), "red"),
]:
    z_pts = []
    r_pts = []
    
    for z_c in z_windows:
        mask = rm_mask & (z_m >= z_c - window) & (z_m < z_c + window)
        r, n, _ = corr_pair(mg_m, hb_m, mask)
        if not np.isnan(r) and n >= 15:
            z_pts.append(z_c)
            r_pts.append(r)
    
    z_pts = np.array(z_pts)
    r_pts = np.array(r_pts)
    
    if len(z_pts) > 5:
        try:
            popt, _ = curve_fit(sigmoid, z_pts, r_pts, p0=[0.8, 10, 1.0], maxfev=5000)
            print(f"\n  {label}: z₀ = {popt[2]:.3f}, k = {popt[1]:.1f}, ρ₀ = {popt[0]:.3f}")
            print(f"    Data points: {len(z_pts)}, z range: [{z_pts.min():.2f}, {z_pts.max():.2f}]")
            
            # Print the curve
            print(f"    {'z':>6} {'ρ_data':>8} {'ρ_fit':>8}")
            for zp, rp in zip(z_pts, r_pts):
                rf = sigmoid(zp, *popt)
                print(f"    {zp:>6.2f} {rp:>+8.3f} {rf:>+8.3f}")
                
        except Exception as e:
            print(f"\n  {label}: Sigmoid fit failed ({e})")
            print(f"    Data: {list(zip(z_pts.round(2), np.round(r_pts, 3)))}")

# ============================================================
# RM vs Z — is RM itself z-dependent? (confound check)
# ============================================================
print(f"\n{'='*70}")
print("CONFOUND CHECK: Does |RM| correlate with z?")
print("If radio-loud quasars at different z have different RM distributions,")
print("the z>1 result could be a selection effect")
print(f"{'='*70}")

r_rmz, p_rmz = stats.spearmanr(z_m, abs_rm)
print(f"\n  Spearman(|RM|, z) = {r_rmz:+.3f} (p = {p_rmz:.3e})")

for zlo, zhi in [(0.5, 0.8), (0.8, 1.0), (1.0, 1.3), (1.3, 2.0)]:
    zmask = (z_m >= zlo) & (z_m < zhi)
    if zmask.sum() > 10:
        med_rm = np.median(abs_rm[zmask])
        print(f"  z=[{zlo},{zhi}): N={zmask.sum()}, median |RM| = {med_rm:.1f}")

# Save
with open('results_rm_deep/rm_deep_results.json', 'w') as f:
    json.dump({
        'n_matched': int(matched.sum()),
        'mw_subtraction': {
            'median_rm_high_b': float(median_rm_high_b),
            'median_rm_low_b': float(median_rm_low_b)
        }
    }, f, indent=2)

print(f"\nResults saved to results_rm_deep/")
