#!/usr/bin/env python3
"""
RM CROSS-MATCH — The Test Everyone Agreed On

Cross-match DR16Q quasars against Taylor+2009 NVSS RM catalog (37,543 sources).
Then: do quasars behind high-|RM| sightlines show more correlational degradation
at fixed z than quasars behind low-|RM| sightlines?

Two approaches:
1. DIRECT MATCH: Find DR16Q quasars that ARE in the RM catalog (radio-loud subset)
2. NEAREST-NEIGHBOR RM: For ALL DR16Q quasars, assign the RM of the nearest
   RM-catalog source as a sightline proxy (the RM sky is smooth on ~1° scales)
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_rm_crossmatch', exist_ok=True)

print("=" * 70)
print("RM CROSS-MATCH — Does the magnetic path predict correlational degradation?")
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
# RM catalog has sexagesimal coords — use SkyCoord to parse
coords_rm = SkyCoord(
    ra=[str(r) for r in rm_cat['RAJ2000']], 
    dec=[str(d) for d in rm_cat['DEJ2000']], 
    unit=(u.hourangle, u.degree), frame='icrs'
)
ra_rm = coords_rm.ra.degree
dec_rm = coords_rm.dec.degree
rm_val = np.array(rm_cat['RM'], dtype=float)
rm_err = np.array(rm_cat['e_RM'], dtype=float)

print(f"\nDR16Q: {len(z_all)} quasars")
print(f"NVSS RM: {len(rm_val)} sources")
print(f"|RM| range: {np.abs(rm_val).min():.1f} - {np.abs(rm_val).max():.1f} rad/m²")

# Build sky coordinates
print("\nBuilding sky coordinates...")
coords_q = SkyCoord(ra=ra_q*u.degree, dec=dec_q*u.degree, frame='icrs')
# coords_rm already built above during coordinate parsing

# ============================================================
# APPROACH 1: Direct match (within 30 arcsec)
# ============================================================
print(f"\n{'='*70}")
print("APPROACH 1: Direct cross-match (30 arcsec radius)")
print("DR16Q quasars that are ALSO in the NVSS RM catalog")
print(f"{'='*70}")

idx_rm, sep, _ = coords_q.match_to_catalog_sky(coords_rm)
matched = sep < 30 * u.arcsec

print(f"Direct matches within 30\": {matched.sum()}")

if matched.sum() > 100:
    z_matched = z_all[matched]
    rm_matched = rm_val[idx_rm[matched]]
    mg_matched = mg_ew[matched]
    hb_matched = hb_ew[matched]
    ciii_matched = ciii_ew[matched]
    
    abs_rm_matched = np.abs(rm_matched)
    rm_median = np.median(abs_rm_matched)
    
    print(f"|RM| for matched quasars: median = {rm_median:.1f} rad/m²")
    
    # Split by |RM| at fixed z
    print(f"\n  MgII↔Hβ: LOW |RM| (< median) vs HIGH |RM| (> median)")
    print(f"  {'z-bin':<15} {'Low |RM|':>12} {'High |RM|':>12} {'Δ':>8}")
    print(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*8}")
    
    z_bins_direct = [(0.5, 0.8), (0.8, 1.0), (1.0, 1.3)]
    
    for zlo, zhi in z_bins_direct:
        zmask = (z_matched >= zlo) & (z_matched < zhi)
        lo_rm = zmask & (abs_rm_matched <= rm_median)
        hi_rm = zmask & (abs_rm_matched > rm_median)
        
        v1_lo, v2_lo = mg_matched[lo_rm], hb_matched[lo_rm]
        v1_hi, v2_hi = mg_matched[hi_rm], hb_matched[hi_rm]
        
        val_lo = np.isfinite(v1_lo) & np.isfinite(v2_lo) & (v1_lo != 0) & (v2_lo != 0)
        val_hi = np.isfinite(v1_hi) & np.isfinite(v2_hi) & (v1_hi != 0) & (v2_hi != 0)
        
        row = f"z=[{zlo},{zhi})"
        if val_lo.sum() > 20 and val_hi.sum() > 20:
            r_lo, _ = stats.spearmanr(v1_lo[val_lo], v2_lo[val_lo])
            r_hi, _ = stats.spearmanr(v1_hi[val_hi], v2_hi[val_hi])
            delta = r_hi - r_lo
            print(f"  {row:<15} {r_lo:>+12.3f} {r_hi:>+12.3f} {delta:>+8.3f}")
        else:
            print(f"  {row:<15} {'N/A':>12} {'N/A':>12} {'N/A':>8}")

# ============================================================
# APPROACH 2: Nearest-neighbor RM assignment (all DR16Q)
# ============================================================
print(f"\n{'='*70}")
print("APPROACH 2: Nearest-neighbor RM (sightline proxy)")
print("Assign each quasar the |RM| of nearest NVSS RM source")
print("RM sky is smooth on ~1° scales → nearby RM is a sightline proxy")
print(f"{'='*70}")

# Already have idx_rm and sep from the match
nn_rm = np.abs(rm_val[idx_rm])
nn_sep = sep.to(u.degree).value

print(f"Nearest RM source: median separation = {np.median(nn_sep):.2f}°")
print(f"90th percentile separation: {np.percentile(nn_sep, 90):.2f}°")

# Use only quasars with a reasonably close RM neighbor (< 1°)
close_enough = nn_sep < 1.0
print(f"Quasars with RM neighbor < 1°: {close_enough.sum()} ({100*close_enough.mean():.1f}%)")

# Split into |RM| terciles
nn_rm_close = nn_rm[close_enough]
rm_t1 = np.percentile(nn_rm_close, 33)
rm_t2 = np.percentile(nn_rm_close, 67)
print(f"|RM| tercile boundaries: {rm_t1:.1f}, {rm_t2:.1f} rad/m²")

def corr_masked(v1_all, v2_all, mask):
    v1, v2 = v1_all[mask], v2_all[mask]
    valid = np.isfinite(v1) & np.isfinite(v2) & (v1 != 0) & (v2 != 0)
    if valid.sum() > 30:
        r, p = stats.spearmanr(v1[valid], v2[valid])
        return r, int(valid.sum()), p
    return np.nan, 0, np.nan

# MgII↔Hβ by |RM| tercile at fixed z
z_bins = [(0.5, 0.7), (0.7, 0.85), (0.85, 0.95), (0.95, 1.05), (1.05, 1.2)]
rm_bins = [
    (0, rm_t1, f"Low |RM|(<{rm_t1:.0f})"),
    (rm_t1, rm_t2, f"Mid |RM|({rm_t1:.0f}-{rm_t2:.0f})"),
    (rm_t2, 9999, f"High |RM|(>{rm_t2:.0f})"),
]

print(f"\n  MgII↔Hβ by sightline |RM| tercile:")
print(f"  {'z-bin':<15}", end="")
for _, _, label in rm_bins:
    print(f"  {label:>20}", end="")
print(f"  {'Δ(hi-lo)':>10}")
print(f"  {'-'*15}", end="")
for _ in rm_bins:
    print(f"  {'-'*20}", end="")
print(f"  {'-'*10}")

results = {}

for zlo, zhi in z_bins:
    zmask_base = close_enough & (z_all >= zlo) & (z_all < zhi)
    
    row = f"z=[{zlo},{zhi})"
    print(f"  {row:<15}", end="")
    
    corrs = []
    for rlo, rhi, label in rm_bins:
        mask = zmask_base & (nn_rm >= rlo) & (nn_rm < rhi)
        r, n, p = corr_masked(mg_ew, hb_ew, mask)
        corrs.append(r)
        r_str = f"{r:+.3f} (N={n})" if not np.isnan(r) else "N/A"
        print(f"  {r_str:>20}", end="")
    
    if not np.isnan(corrs[0]) and not np.isnan(corrs[-1]):
        delta = corrs[-1] - corrs[0]
        print(f"  {delta:>+10.3f}")
        results[row] = {'low_rm': float(corrs[0]), 'high_rm': float(corrs[-1]), 'delta': float(delta)}
    else:
        print(f"  {'N/A':>10}")

# ============================================================
# MgII↔CIII (UV control) — should NOT be affected
# ============================================================
print(f"\n  MgII↔CIII (UV control — should be FLAT across |RM|):")
print(f"  {'z-bin':<15}", end="")
for _, _, label in rm_bins:
    print(f"  {label:>20}", end="")
print(f"  {'Δ(hi-lo)':>10}")
print(f"  {'-'*15}", end="")
for _ in rm_bins:
    print(f"  {'-'*20}", end="")
print(f"  {'-'*10}")

for zlo, zhi in z_bins:
    zmask_base = close_enough & (z_all >= zlo) & (z_all < zhi)
    
    row = f"z=[{zlo},{zhi})"
    print(f"  {row:<15}", end="")
    
    corrs = []
    for rlo, rhi, label in rm_bins:
        mask = zmask_base & (nn_rm >= rlo) & (nn_rm < rhi)
        r, n, p = corr_masked(mg_ew, ciii_ew, mask)
        corrs.append(r)
        r_str = f"{r:+.3f} (N={n})" if not np.isnan(r) else "N/A"
        print(f"  {r_str:>20}", end="")
    
    if not np.isnan(corrs[0]) and not np.isnan(corrs[-1]):
        delta = corrs[-1] - corrs[0]
        print(f"  {delta:>+10.3f}")
    else:
        print(f"  {'N/A':>10}")

# ============================================================
# CONTINUOUS: Spearman of residual ρ vs |RM|
# ============================================================
print(f"\n{'='*70}")
print("CONTINUOUS TEST: Does |RM| predict MgII↔Hβ correlation strength?")
print("Compute per-object residual from z-trend, regress against |RM|")
print(f"{'='*70}")

# For each narrow z-bin, compute rank-based "contribution to correlation"
# Simple proxy: at each z, does high |RM| → lower rank alignment?

# Split into fine z-bins, compute sliding correlation for high vs low |RM| 
z_fine = np.arange(0.5, 1.3, 0.05)
deltas_fine = []
z_mids_fine = []

for i in range(len(z_fine) - 1):
    zlo, zhi = z_fine[i], z_fine[i+1]
    zmask = close_enough & (z_all >= zlo) & (z_all < zhi)
    
    lo = zmask & (nn_rm <= rm_t1)
    hi = zmask & (nn_rm >= rm_t2)
    
    r_lo, n_lo, _ = corr_masked(mg_ew, hb_ew, lo)
    r_hi, n_hi, _ = corr_masked(mg_ew, hb_ew, hi)
    
    if not np.isnan(r_lo) and not np.isnan(r_hi):
        deltas_fine.append(r_hi - r_lo)
        z_mids_fine.append((zlo + zhi) / 2)

if deltas_fine:
    deltas_fine = np.array(deltas_fine)
    z_mids_fine = np.array(z_mids_fine)
    
    mean_delta = np.mean(deltas_fine)
    print(f"\n  Mean Δ(high_RM - low_RM) across z-bins: {mean_delta:+.3f}")
    print(f"  If negative: high |RM| sightlines → LOWER optical correlations ✓")
    print(f"  If zero: sightline |RM| doesn't matter ✗")
    
    # Does the effect increase near z₀?
    dist_from_z0 = np.abs(z_mids_fine - 1.045)
    r_trend, p_trend = stats.spearmanr(dist_from_z0, np.abs(deltas_fine))
    print(f"\n  |Δ| vs distance from z₀: ρ = {r_trend:+.3f} (p = {p_trend:.3e})")
    
    # Print the fine deltas
    print(f"\n  {'z-mid':>6} {'Δ(hi-lo)':>10} {'Bar':>30}")
    print(f"  {'-'*6} {'-'*10} {'-'*30}")
    for z_m, delta in zip(z_mids_fine, deltas_fine):
        if delta < 0:
            bar = "◄" * int(abs(delta) * 100)
        else:
            bar = "►" * int(abs(delta) * 100)
        print(f"  {z_m:>6.2f} {delta:>+10.3f} {bar}")

# ============================================================
# ABSOLUTE |RM| GRADIENT
# ============================================================
print(f"\n{'='*70}")
print("GRADIENT: Does HIGHER |RM| → monotonically LOWER optical correlation?")
print(f"{'='*70}")

# At z=0.8-1.0 (pre-transition, most data), split into |RM| quintiles
zmask_sweet = close_enough & (z_all >= 0.8) & (z_all < 1.0)
rm_in_range = nn_rm[zmask_sweet]

quintile_edges = np.percentile(rm_in_range[np.isfinite(rm_in_range)], [0, 20, 40, 60, 80, 100])
print(f"\n  |RM| quintile edges: {[f'{e:.1f}' for e in quintile_edges]}")

print(f"\n  {'|RM| range':<20} {'MgII↔Hβ':>10} {'MgII↔CIII':>10} {'N':>8}")
print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*8}")

optical_corrs = []
uv_corrs = []

for i in range(5):
    lo, hi = quintile_edges[i], quintile_edges[i+1]
    mask = zmask_sweet & (nn_rm >= lo) & (nn_rm < hi + 0.01)
    
    r_opt, n_opt, _ = corr_masked(mg_ew, hb_ew, mask)
    r_uv, n_uv, _ = corr_masked(mg_ew, ciii_ew, mask)
    
    optical_corrs.append(r_opt)
    uv_corrs.append(r_uv)
    
    label = f"|RM|=[{lo:.0f},{hi:.0f})"
    opt_str = f"{r_opt:+.3f}" if not np.isnan(r_opt) else "N/A"
    uv_str = f"{r_uv:+.3f}" if not np.isnan(r_uv) else "N/A"
    print(f"  {label:<20} {opt_str:>10} {uv_str:>10} {n_opt:>8}")

# Monotonicity test
valid_opt = [r for r in optical_corrs if not np.isnan(r)]
if len(valid_opt) >= 3:
    r_mono, p_mono = stats.spearmanr(range(len(valid_opt)), valid_opt)
    print(f"\n  Monotonicity test (optical): ρ = {r_mono:+.3f} (p = {p_mono:.3f})")
    if r_mono < -0.5 and p_mono < 0.05:
        print(f"  ★ CONFIRMED: Higher |RM| → lower optical correlation!")
    elif r_mono > 0.5:
        print(f"  ✗ INVERTED: Higher |RM| → HIGHER optical correlation (unexpected)")
    else:
        print(f"  ~ No clear monotonic trend")

valid_uv = [r for r in uv_corrs if not np.isnan(r)]
if len(valid_uv) >= 3:
    r_mono_uv, p_mono_uv = stats.spearmanr(range(len(valid_uv)), valid_uv)
    print(f"  Monotonicity test (UV): ρ = {r_mono_uv:+.3f} (p = {p_mono_uv:.3f})")

# Save
with open('results_rm_crossmatch/rm_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\n\nResults saved to results_rm_crossmatch/")
