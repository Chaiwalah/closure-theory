#!/usr/bin/env python3
"""
MATCHED-PAIR CONFOUND KILLER

The 2,161 RM-matched quasars are radio-loud — a special AGN population.
Could the RM effect be source-tagged (radio-loud AGN have different 
internal line physics) rather than path-tagged?

Kill this by matching high-RM and low-RM quasars on:
- Redshift (within Δz = 0.05)
- BH mass (within 0.3 dex)
- Luminosity (within 0.3 dex)
- Eddington ratio (within 0.3 dex)

If matched pairs still show the RM effect → path, not source.
If matched pairs kill the effect → source physics confound.
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
from scipy.spatial import cKDTree
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_rm_matched', exist_ok=True)

print("=" * 70)
print("MATCHED-PAIR CONFOUND KILLER")
print("Do radio-loud AGN with high RM differ from low RM")
print("even after matching on mass, luminosity, Eddington, z?")
print("=" * 70)

# Load DR16Q
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
ciii_ew = d['CIII_BR'][:, 2]

# BH mass and luminosity
logmbh = d['LOGMBH']
loglbol = d['LOGLBOL']
logedd = d['LOGLEDD_RATIO']

ra_q = d['RA']
dec_q = d['DEC']

# Load NVSS RM catalog
rm_cat = Table.read('data/nvss_rm_taylor2009.fits')
coords_rm = SkyCoord(
    ra=[str(r) for r in rm_cat['RAJ2000']], 
    dec=[str(dd) for dd in rm_cat['DEJ2000']], 
    unit=(u.hourangle, u.degree), frame='icrs'
)
rm_val = np.array(rm_cat['RM'], dtype=float)

# Cross match
coords_q = SkyCoord(ra=ra_q*u.degree, dec=dec_q*u.degree, frame='icrs')
idx_rm, sep, _ = coords_q.match_to_catalog_sky(coords_rm)
matched = sep < 30 * u.arcsec

# Extract matched sample with all properties
z_m = z_all[matched]
rm_m = np.abs(rm_val[idx_rm[matched]])
mg_m = mg_ew[matched]
hb_m = hb_ew[matched]
ciii_m = ciii_ew[matched]
mass_m = logmbh[matched]
lum_m = loglbol[matched]
edd_m = logedd[matched]

print(f"Total RM-matched quasars: {matched.sum()}")

# Filter for objects with valid mass + lum + eddington
has_props = (np.isfinite(mass_m) & (mass_m > 0) & 
             np.isfinite(lum_m) & (lum_m > 0) &
             np.isfinite(edd_m))
print(f"With valid BH mass + luminosity + Eddington: {has_props.sum()}")

z_p = z_m[has_props]
rm_p = rm_m[has_props]
mg_p = mg_m[has_props]
hb_p = hb_m[has_props]
ciii_p = ciii_m[has_props]
mass_p = mass_m[has_props]
lum_p = lum_m[has_props]
edd_p = edd_m[has_props]

# ============================================================
# First: Are high-RM and low-RM quasars different in source properties?
# ============================================================
print(f"\n{'='*70}")
print("PROPERTY CHECK: Do high-RM quasars differ in source properties?")
print(f"{'='*70}")

rm_median = np.median(rm_p)
lo_rm = rm_p <= rm_median
hi_rm = rm_p > rm_median

print(f"\n  {'Property':<20} {'Low |RM|':>12} {'High |RM|':>12} {'Δ':>8} {'p':>10}")
print(f"  {'-'*20} {'-'*12} {'-'*12} {'-'*8} {'-'*10}")

for name, vals in [('z', z_p), ('log M_BH', mass_p), ('log L_bol', lum_p), ('log L/L_Edd', edd_p)]:
    med_lo = np.median(vals[lo_rm])
    med_hi = np.median(vals[hi_rm])
    _, p_ks = stats.ks_2samp(vals[lo_rm], vals[hi_rm])
    print(f"  {name:<20} {med_lo:>12.2f} {med_hi:>12.2f} {med_hi-med_lo:>+8.2f} {p_ks:>10.3f}")

# ============================================================
# MATCHED PAIRS using nearest-neighbor in property space
# ============================================================
print(f"\n{'='*70}")
print("BUILDING MATCHED PAIRS")
print("For each high-RM quasar, find closest low-RM match in")
print("(z, log_M_BH, log_L_bol, log_Edd) space")
print(f"{'='*70}")

# Normalize properties to unit variance for matching
def normalize(arr):
    return (arr - np.nanmean(arr)) / np.nanstd(arr)

# Build feature vectors
features_lo = np.column_stack([
    normalize(z_p[lo_rm]),
    normalize(mass_p[lo_rm]),
    normalize(lum_p[lo_rm]),
    normalize(edd_p[lo_rm])
])

features_hi = np.column_stack([
    normalize(z_p[hi_rm]),
    normalize(mass_p[hi_rm]),
    normalize(lum_p[hi_rm]),
    normalize(edd_p[hi_rm])
])

# Remove any rows with NaN
valid_lo = np.all(np.isfinite(features_lo), axis=1)
valid_hi = np.all(np.isfinite(features_hi), axis=1)

features_lo_clean = features_lo[valid_lo]
features_hi_clean = features_hi[valid_hi]

# Indices back to the has_props subset
idx_lo_all = np.where(lo_rm)[0][valid_lo]
idx_hi_all = np.where(hi_rm)[0][valid_hi]

print(f"  Valid low-RM: {len(features_lo_clean)}")
print(f"  Valid high-RM: {len(features_hi_clean)}")

# Build KD-tree on low-RM quasars
tree = cKDTree(features_lo_clean)

# For each high-RM quasar, find nearest low-RM match
distances, nn_indices = tree.query(features_hi_clean, k=1)

# Quality filter: only keep matches with distance < 0.5 in normalized space
good_match = distances < 0.5
print(f"  Good matches (distance < 0.5σ): {good_match.sum()} / {len(distances)}")

# Also try stricter matching
very_good = distances < 0.3
print(f"  Very good matches (< 0.3σ): {very_good.sum()}")

# Get paired indices
hi_paired = idx_hi_all[good_match]
lo_paired = idx_lo_all[nn_indices[good_match]]

print(f"\n  Match quality check:")
print(f"  {'Property':<15} {'|Δ| median':>12} {'|Δ| 90th':>12}")
print(f"  {'-'*15} {'-'*12} {'-'*12}")
for name, vals in [('z', z_p), ('log M_BH', mass_p), ('log L_bol', lum_p), ('log L/L_Edd', edd_p)]:
    deltas = np.abs(vals[hi_paired] - vals[lo_paired])
    print(f"  {name:<15} {np.median(deltas):>12.3f} {np.percentile(deltas, 90):>12.3f}")

# ============================================================
# THE TEST: Do matched pairs still show the RM effect?
# ============================================================
print(f"\n{'='*70}")
print("THE TEST: MgII↔Hβ in matched pairs")
print("If RM effect survives matching → PATH, not source")
print("If RM effect dies → SOURCE PHYSICS confound")
print(f"{'='*70}")

def corr_sub(v1, v2, idx):
    d1, d2 = v1[idx], v2[idx]
    valid = np.isfinite(d1) & np.isfinite(d2) & (d1 != 0) & (d2 != 0)
    if valid.sum() > 15:
        r, p = stats.spearmanr(d1[valid], d2[valid])
        return r, int(valid.sum()), p
    return np.nan, 0, np.nan

# Overall
r_hi, n_hi, p_hi = corr_sub(mg_p, hb_p, hi_paired)
r_lo, n_lo, p_lo = corr_sub(mg_p, hb_p, lo_paired)

print(f"\n  OVERALL (all z):")
print(f"    High |RM| (matched): ρ = {r_hi:+.3f} (N={n_hi})")
print(f"    Low |RM| (matched):  ρ = {r_lo:+.3f} (N={n_lo})")
print(f"    Δ = {r_hi - r_lo:+.3f}")

# By z bin
print(f"\n  BY REDSHIFT:")
print(f"  {'z-bin':<15} {'High |RM|':>12} {'Low |RM|':>12} {'Δ':>8} {'N_hi':>6} {'N_lo':>6}")
print(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*8} {'-'*6} {'-'*6}")

z_bins = [(0.4, 0.7), (0.7, 0.9), (0.9, 1.05), (1.05, 1.3), (1.3, 2.0)]

for zlo, zhi in z_bins:
    z_mask_hi = (z_p[hi_paired] >= zlo) & (z_p[hi_paired] < zhi)
    z_mask_lo = (z_p[lo_paired] >= zlo) & (z_p[lo_paired] < zhi)
    
    hi_idx = hi_paired[z_mask_hi]
    lo_idx = lo_paired[z_mask_lo]
    
    r_h, n_h, _ = corr_sub(mg_p, hb_p, hi_idx)
    r_l, n_l, _ = corr_sub(mg_p, hb_p, lo_idx)
    
    row = f"z=[{zlo},{zhi})"
    if not np.isnan(r_h) and not np.isnan(r_l):
        delta = r_h - r_l
        print(f"  {row:<15} {r_h:>+12.3f} {r_l:>+12.3f} {delta:>+8.3f} {n_h:>6} {n_l:>6}")
    else:
        h_str = f"{r_h:+.3f}" if not np.isnan(r_h) else "N/A"
        l_str = f"{r_l:+.3f}" if not np.isnan(r_l) else "N/A"
        print(f"  {row:<15} {h_str:>12} {l_str:>12} {'N/A':>8} {n_h:>6} {n_l:>6}")

# ============================================================
# WAVELENGTH CHECK in matched pairs
# ============================================================
print(f"\n{'='*70}")
print("WAVELENGTH CHECK: Does λ² dependence survive matching?")
print(f"{'='*70}")

# At z=0.85-1.15 (transition zone)
z_trans_hi = (z_p[hi_paired] >= 0.85) & (z_p[hi_paired] < 1.15)
z_trans_lo = (z_p[lo_paired] >= 0.85) & (z_p[lo_paired] < 1.15)

pairs = [
    ('MgII↔Hβ (opt)', mg_p, hb_p, 4861),
    ('MgII↔CIII (UV)', mg_p, ciii_p, 1909),
]

print(f"\n  At z = 0.85-1.15 (matched pairs only):")
print(f"  {'Pair':<20} {'High |RM|':>10} {'Low |RM|':>10} {'Δ':>8}")
print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*8}")

for name, v1, v2, lam in pairs:
    r_h, n_h, _ = corr_sub(v1, v2, hi_paired[z_trans_hi])
    r_l, n_l, _ = corr_sub(v1, v2, lo_paired[z_trans_lo])
    
    if not np.isnan(r_h) and not np.isnan(r_l):
        delta = r_h - r_l
        print(f"  {name:<20} {r_h:>+10.3f} {r_l:>+10.3f} {delta:>+8.3f}")
    else:
        print(f"  {name:<20} {'N/A':>10} {'N/A':>10} {'N/A':>8}")

# ============================================================
# BOOTSTRAP: Is the effect significant?
# ============================================================
print(f"\n{'='*70}")
print("BOOTSTRAP SIGNIFICANCE (1000 iterations)")
print(f"{'='*70}")

# Focus on z=0.8-1.2 where the effect should be strongest
z_focus_hi = (z_p[hi_paired] >= 0.8) & (z_p[hi_paired] < 1.2)
z_focus_lo = (z_p[lo_paired] >= 0.8) & (z_p[lo_paired] < 1.2)

hi_focus = hi_paired[z_focus_hi]
lo_focus = lo_paired[z_focus_lo]

if len(hi_focus) > 20 and len(lo_focus) > 20:
    # Bootstrap: shuffle RM labels and recompute Δ
    n_boot = 1000
    observed_delta = (corr_sub(mg_p, hb_p, hi_focus)[0] - 
                      corr_sub(mg_p, hb_p, lo_focus)[0])
    
    all_focus = np.concatenate([hi_focus, lo_focus])
    null_deltas = []
    
    for _ in range(n_boot):
        np.random.shuffle(all_focus)
        mid = len(hi_focus)
        r1, _, _ = corr_sub(mg_p, hb_p, all_focus[:mid])
        r2, _, _ = corr_sub(mg_p, hb_p, all_focus[mid:mid+len(lo_focus)])
        if not np.isnan(r1) and not np.isnan(r2):
            null_deltas.append(r1 - r2)
    
    null_deltas = np.array(null_deltas)
    p_boot = np.mean(np.abs(null_deltas) >= np.abs(observed_delta))
    
    print(f"\n  Observed Δ at z=0.8-1.2: {observed_delta:+.3f}")
    print(f"  Null distribution: mean = {np.mean(null_deltas):+.3f}, std = {np.std(null_deltas):.3f}")
    print(f"  Bootstrap p-value: {p_boot:.3f}")
    print(f"  Z-score: {(observed_delta - np.mean(null_deltas)) / np.std(null_deltas):+.1f}σ")
    
    if p_boot < 0.05:
        print(f"  ★ SIGNIFICANT: RM effect SURVIVES matching at p < 0.05")
    elif p_boot < 0.1:
        print(f"  ~ MARGINAL: p < 0.10, suggestive but not conclusive")
    else:
        print(f"  ✗ NOT SIGNIFICANT: RM effect may be a source-physics confound")

# Save
with open('results_rm_matched/matched_pair_results.json', 'w') as f:
    json.dump({
        'n_hi_paired': int(len(hi_paired)),
        'n_lo_paired': int(len(lo_paired)),
        'match_quality_median_dz': float(np.median(np.abs(z_p[hi_paired] - z_p[lo_paired]))),
    }, f, indent=2)

print(f"\nResults saved to results_rm_matched/")
