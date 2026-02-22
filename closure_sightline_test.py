#!/usr/bin/env python3
"""
Sightline Test: Void-Crossing vs Filament-Crossing
====================================================
The PATH test: does the light's journey matter, not just its origin?

1. Get high-z quasars (z=1.5-2.5) with CIV+MgII
2. Get foreground galaxies (z=0.1-1.0) as structure tracers
3. For each quasar, count foreground galaxies within angular cone
4. "Void-crossing" = few foreground galaxies, "filament-crossing" = many
5. Compare CIV/MgII correlation for each class

If void-crossing sightlines show worse correlations → the channel medium 
degrades information during transit. Not source, not instrument. The path.
"""
import numpy as np
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
import requests, json, os, sys, warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = 'results_sightline'
os.makedirs(RESULTS_DIR, exist_ok=True)

def query_desi(sql, timeout=600):
    print(f"  Querying...", flush=True)
    r = requests.get('https://datalab.noirlab.edu/tap/sync',
                     params={'REQUEST':'doQuery','LANG':'ADQL','QUERY':sql,'FORMAT':'csv'},
                     timeout=timeout)
    if r.status_code != 200 or 'ERROR' in r.text[:500]:
        print(f"  Error: {r.text[:300]}", flush=True)
        return None
    lines = r.text.strip().split('\n')
    print(f"  Rows: {len(lines)-1}", flush=True)
    data = []
    for line in lines[1:]:
        try:
            data.append([float(p) for p in line.split(',')])
        except:
            continue
    return np.array(data) if len(data) > 100 else None

print("=" * 70, flush=True)
print("SIGHTLINE TEST: VOID-CROSSING vs FILAMENT-CROSSING", flush=True)
print("=" * 70, flush=True)

# Step 1: Get background quasars (z=1.5-2.5)
print("\n[1] Background quasars (z=1.5-2.5, CIV+MgII)...", flush=True)
sql_qso = """SELECT TOP 100000
    q.target_ra, q.target_dec, q.z,
    q.civ_1549_flux, q.civ_1549_flux_ivar,
    q.mgii_2796_flux, q.mgii_2796_flux_ivar
FROM desi_dr1.agnqso q
WHERE q.z BETWEEN 1.5 AND 2.5
  AND q.civ_1549_flux > 0 AND q.civ_1549_flux_ivar > 0
  AND q.mgii_2796_flux > 0 AND q.mgii_2796_flux_ivar > 0
  AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
  AND q.mgii_2796_flux * SQRT(q.mgii_2796_flux_ivar) > 3
ORDER BY q.random_id"""

qso = query_desi(sql_qso)
if qso is None:
    print("FATAL: QSO query failed", flush=True)
    sys.exit(1)

qso_ra = qso[:, 0]
qso_dec = qso[:, 1]
qso_z = qso[:, 2]
civ_flux = qso[:, 3]
civ_ivar = qso[:, 4]
mgii_flux = qso[:, 5]
mgii_ivar = qso[:, 6]
civ_snr = civ_flux * np.sqrt(np.maximum(civ_ivar, 1e-30))
mgii_snr = mgii_flux * np.sqrt(np.maximum(mgii_ivar, 1e-30))

print(f"  Got {len(qso)} background quasars", flush=True)

# Step 2: Get foreground galaxies as structure tracers
# Use zpix table which has all DESI targets with redshifts
print("\n[2] Foreground structure tracers (z=0.1-1.0)...", flush=True)

# Query foreground galaxies — use zpix for general galaxy positions
# We need a LOT for density estimation, but limit query size
sql_fg = """SELECT TOP 500000
    z.mean_fiber_ra, z.mean_fiber_dec, z.z
FROM desi_dr1.zpix z
WHERE z.z BETWEEN 0.1 AND 1.0
  AND z.zwarn = 0
  AND z.spectype = 'GALAXY'
  AND z.zcat_primary = 't'
ORDER BY z.random_id"""

fg = query_desi(sql_fg)

if fg is None:
    print("FATAL: Foreground query failed", flush=True)
    sys.exit(1)

fg_ra = fg[:, 0]
fg_dec = fg[:, 1]
fg_z = fg[:, 2]
print(f"  Got {len(fg)} foreground galaxies", flush=True)
print(f"  z range: {fg_z.min():.2f} - {fg_z.max():.2f}", flush=True)

# Step 3: For each quasar, count foreground galaxies within angular cone
# Use a cone of ~10 arcmin (~0.17 deg) — typical void angular scale at z~0.5
print("\n[3] Computing foreground density along each quasar sightline...", flush=True)

# Build 2D tree of foreground galaxies (RA, Dec)
fg_coords = np.column_stack([fg_ra, fg_dec])
fg_tree = cKDTree(fg_coords)

# For each quasar, count foreground galaxies within search radius
search_radius_deg = 0.25  # ~15 arcmin cone
print(f"  Search radius: {search_radius_deg:.2f} deg ({search_radius_deg*60:.0f} arcmin)", flush=True)

qso_coords = np.column_stack([qso_ra, qso_dec])

# Count neighbors (this is a 2D angular count — foreground density proxy)
print(f"  Counting foreground neighbors for {len(qso)} quasars...", flush=True)
fg_counts = fg_tree.query_ball_point(qso_coords, r=search_radius_deg)
n_fg = np.array([len(c) for c in fg_counts])

print(f"  Foreground count range: {n_fg.min()} - {n_fg.max()}", flush=True)
print(f"  Median: {np.median(n_fg):.0f}, Mean: {np.mean(n_fg):.1f}", flush=True)

# Step 4: Split by foreground density and compare correlations
print("\n" + "=" * 70, flush=True)
print("SIGHTLINE SPLIT: VOID-CROSSING vs FILAMENT-CROSSING", flush=True)
print("=" * 70, flush=True)

z_edges = np.array([1.5, 1.7, 1.9, 2.1, 2.3, 2.5])
combined_snr = np.sqrt(civ_snr**2 + mgii_snr**2)
lum = np.log10(civ_flux + mgii_flux + 1e-20)

results = []

for i in range(len(z_edges)-1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    zmask = (qso_z >= z_lo) & (qso_z < z_hi)
    if zmask.sum() < 500:
        print(f"  z=[{z_lo:.1f},{z_hi:.1f}]: {zmask.sum()} — skip", flush=True)
        continue
    
    # Median split by foreground count
    med_fg = np.median(n_fg[zmask])
    filament_m = zmask & (n_fg >= med_fg)  # many foreground = filament-crossing
    void_m = zmask & (n_fg < med_fg)  # few foreground = void-crossing
    
    # RAW
    rho_fil_raw = spearmanr(civ_flux[filament_m], mgii_flux[filament_m])[0]
    rho_void_raw = spearmanr(civ_flux[void_m], mgii_flux[void_m])[0]
    
    # MATCHED on SNR + luminosity
    f_idx, v_idx = np.where(filament_m)[0], np.where(void_m)[0]
    matched_f, matched_v = [], []
    nb = 10
    snr_bins = np.percentile(combined_snr[zmask], np.linspace(0,100,nb+1))
    lum_bins = np.percentile(lum[zmask], np.linspace(0,100,nb+1))
    for si in range(nb):
        for li in range(nb):
            f_in = f_idx[(combined_snr[f_idx]>=snr_bins[si])&(combined_snr[f_idx]<snr_bins[si+1])&
                         (lum[f_idx]>=lum_bins[li])&(lum[f_idx]<lum_bins[li+1])]
            v_in = v_idx[(combined_snr[v_idx]>=snr_bins[si])&(combined_snr[v_idx]<snr_bins[si+1])&
                         (lum[v_idx]>=lum_bins[li])&(lum[v_idx]<lum_bins[li+1])]
            nm = min(len(f_in), len(v_in))
            if nm > 0:
                matched_f.extend(np.random.choice(f_in, nm, replace=False))
                matched_v.extend(np.random.choice(v_in, nm, replace=False))
    
    if len(matched_f) < 100:
        print(f"  z=[{z_lo:.1f},{z_hi:.1f}]: matched {len(matched_f)} — skip", flush=True)
        continue
    
    matched_f, matched_v = np.array(matched_f), np.array(matched_v)
    rho_fil = spearmanr(civ_flux[matched_f], mgii_flux[matched_f])[0]
    rho_void = spearmanr(civ_flux[matched_v], mgii_flux[matched_v])[0]
    delta = rho_fil - rho_void  # positive = filament preserves better
    
    # Check matching
    fg_gap = abs(np.median(n_fg[matched_f]) - np.median(n_fg[matched_v]))
    
    print(f"  z=[{z_lo:.1f},{z_hi:.1f}] N={zmask.sum()} matched={len(matched_f)}", flush=True)
    print(f"    Raw:     filament={rho_fil_raw:.3f} void={rho_void_raw:.3f} Δ={rho_fil_raw-rho_void_raw:+.4f}", flush=True)
    print(f"    Matched: filament={rho_fil:.3f} void={rho_void:.3f} Δ={delta:+.4f}", flush=True)
    print(f"    Median fg count: filament={np.median(n_fg[matched_f]):.0f} void={np.median(n_fg[matched_v]):.0f}", flush=True)
    
    results.append({
        'z_range': f"{z_lo:.1f}-{z_hi:.1f}",
        'n_total': int(zmask.sum()),
        'n_matched': int(len(matched_f)),
        'raw_filament': round(rho_fil_raw, 4),
        'raw_void': round(rho_void_raw, 4),
        'raw_delta': round(rho_fil_raw - rho_void_raw, 4),
        'matched_filament': round(rho_fil, 4),
        'matched_void': round(rho_void, 4),
        'matched_delta': round(delta, 4),
        'median_fg_filament': int(np.median(n_fg[matched_f])),
        'median_fg_void': int(np.median(n_fg[matched_v])),
    })

# SUMMARY
print("\n" + "=" * 70, flush=True)
print("SUMMARY", flush=True)
print("=" * 70, flush=True)

if results:
    n_pos = sum(1 for r in results if r['matched_delta'] > 0)
    mean_d = np.mean([r['matched_delta'] for r in results])
    print(f"  Filament-crossing > Void-crossing: {n_pos}/{len(results)} bins", flush=True)
    print(f"  Mean Δρ (filament - void): {mean_d:+.4f}", flush=True)
    
    if n_pos >= len(results) * 0.6 and mean_d > 0.003:
        print(f"\n  🏆 SIGHTLINE TEST: PASS — the PATH matters, not just the source", flush=True)
    elif mean_d > 0:
        print(f"\n  ⚠️  SIGHTLINE TEST: SUGGESTIVE but weak", flush=True)
    else:
        print(f"\n  ❌ SIGHTLINE TEST: NO EFFECT — path doesn't matter", flush=True)
else:
    print("  No results computed", flush=True)

output = {
    'test': 'Sightline: Void-Crossing vs Filament-Crossing',
    'method': 'Foreground galaxy count within 0.25 deg cone as structure proxy',
    'background': 'DESI agnqso z=1.5-2.5 with CIV+MgII',
    'foreground': 'DESI zpix galaxies z=0.1-1.0',
    'results': results,
}
with open(f'{RESULTS_DIR}/sightline_results.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {RESULTS_DIR}/sightline_results.json", flush=True)
