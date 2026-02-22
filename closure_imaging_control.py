#!/usr/bin/env python3
"""
Imaging Depth / E(B-V) Control Test
=====================================
THE LAST DEFENSIVE TEST: Does our 3D density proxy correlate with survey quality?

If "dense" regions are just regions with deeper imaging or less extinction,
the environment effect is a survey artifact, not physics.

Tests:
1. Query DESI photometry columns (FLUX_G/R/Z ivar as depth proxy, EBV)
2. Correlate 3D density with imaging quality metrics
3. If correlated: rerun env split matching on imaging quality too
4. If uncorrelated: last mundane explanation dead

Also:
5. Galactic latitude control — does density proxy track |b|?
"""
import numpy as np
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid
import requests, json, os, sys, warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = 'results_imaging_control'
os.makedirs(RESULTS_DIR, exist_ok=True)

def query_desi(sql, timeout=600):
    print(f"  Querying DESI...", flush=True)
    r = requests.get('https://datalab.noirlab.edu/tap/sync',
                     params={'REQUEST':'doQuery','LANG':'ADQL','QUERY':sql,'FORMAT':'csv'},
                     timeout=timeout)
    if r.status_code != 200 or 'ERROR' in r.text[:500]:
        print(f"  Error: {r.text[:300]}", flush=True)
        return None, None
    lines = r.text.strip().split('\n')
    header = lines[0].split(',')
    print(f"  Columns: {header}", flush=True)
    print(f"  Rows: {len(lines)-1}", flush=True)
    data = []
    for line in lines[1:]:
        try:
            data.append([float(p) for p in line.split(',')])
        except:
            continue
    return (np.array(data), header) if len(data) > 100 else (None, header)

def compute_3d_density(ra, dec, z, n_neighbors=10):
    z_fine = np.linspace(0, max(z)*1.1, 5000)
    Hz = np.sqrt(0.3*(1+z_fine)**3 + 0.7)
    chi_fine = cumulative_trapezoid(1.0/Hz, z_fine, initial=0) * 2997.9
    chi = np.interp(z, z_fine, chi_fine)
    ra_r, dec_r = np.radians(ra), np.radians(dec)
    x = chi * np.cos(dec_r) * np.cos(ra_r)
    y = chi * np.cos(dec_r) * np.sin(ra_r)
    zc = chi * np.sin(dec_r)
    tree = cKDTree(np.column_stack([x, y, zc]))
    dists, _ = tree.query(np.column_stack([x, y, zc]), k=n_neighbors+1)
    return n_neighbors / (4/3 * np.pi * dists[:, -1]**3 + 1e-10)

print("=" * 70, flush=True)
print("IMAGING DEPTH / E(B-V) CONTROL TEST", flush=True)
print("=" * 70, flush=True)

# First: check what photometry columns exist in agnqso
print("\n[0] Checking available imaging columns...", flush=True)
sql_cols = """SELECT column_name FROM tap_schema.columns 
WHERE table_name='desi_dr1.agnqso' 
AND (column_name LIKE '%ebv%' OR column_name LIKE '%EBV%'
     OR column_name LIKE '%nobs%' OR column_name LIKE '%depth%'
     OR column_name LIKE '%psf%' OR column_name LIKE '%galdepth%'
     OR column_name LIKE '%flux_g%' OR column_name LIKE '%flux_r%' 
     OR column_name LIKE '%flux_z%' OR column_name LIKE '%flux_w%'
     OR column_name LIKE '%ivar%' OR column_name LIKE '%fibertot%'
     OR column_name LIKE '%snr%' OR column_name LIKE '%tsnr%')
ORDER BY column_name"""

cols_data, _ = query_desi(sql_cols, timeout=30)
if cols_data is None:
    # Try fetching column names differently
    print("  Trying text-based column fetch...", flush=True)
    r = requests.get('https://datalab.noirlab.edu/tap/sync',
                     params={'REQUEST':'doQuery','LANG':'ADQL',
                             'QUERY':sql_cols,'FORMAT':'csv'}, timeout=30)
    print(f"  Raw: {r.text[:500]}", flush=True)

# Let's also check the zpix table which has photometry
print("\n[0b] Checking zpix table for imaging columns...", flush=True)
sql_cols2 = """SELECT column_name FROM tap_schema.columns 
WHERE table_name='desi_dr1.zpix' 
AND (column_name LIKE '%ebv%' OR column_name LIKE '%tsnr%' 
     OR column_name LIKE '%fiber%' OR column_name LIKE '%flux%')
ORDER BY column_name"""
r2 = requests.get('https://datalab.noirlab.edu/tap/sync',
                   params={'REQUEST':'doQuery','LANG':'ADQL',
                           'QUERY':sql_cols2,'FORMAT':'csv'}, timeout=30)
print(f"  zpix columns: {r2.text[:500]}", flush=True)

# Check photometry table
print("\n[0c] Checking photometry table...", flush=True)
sql_cols3 = """SELECT column_name FROM tap_schema.columns 
WHERE table_name='desi_dr1.photometry' 
AND (column_name LIKE '%ebv%' OR column_name LIKE '%depth%' 
     OR column_name LIKE '%nobs%' OR column_name LIKE '%psf%'
     OR column_name LIKE '%gal%')
ORDER BY column_name"""
r3 = requests.get('https://datalab.noirlab.edu/tap/sync',
                   params={'REQUEST':'doQuery','LANG':'ADQL',
                           'QUERY':sql_cols3,'FORMAT':'csv'}, timeout=30)
print(f"  photometry columns: {r3.text[:500]}", flush=True)

# Also check what tables exist
print("\n[0d] Checking available tables...", flush=True)
sql_tables = """SELECT table_name FROM tap_schema.tables 
WHERE schema_name='desi_dr1' ORDER BY table_name"""
r4 = requests.get('https://datalab.noirlab.edu/tap/sync',
                   params={'REQUEST':'doQuery','LANG':'ADQL',
                           'QUERY':sql_tables,'FORMAT':'csv'}, timeout=30)
print(f"  Tables: {r4.text[:800]}", flush=True)

# Now query with what we can get — at minimum we have target_ra/dec 
# which gives us galactic coordinates for E(B-V) proxy
print("\n[1] Querying main sample with galactic latitude proxy...", flush=True)

# We can compute galactic |b| from RA/Dec as an extinction proxy
# Also try to join with photometry if it exists
sql_main = """SELECT TOP 150000
    q.target_ra, q.target_dec, q.z,
    q.civ_1549_flux, q.civ_1549_flux_ivar,
    q.mgii_2796_flux, q.mgii_2796_flux_ivar
FROM desi_dr1.agnqso q
WHERE q.z BETWEEN 1.3 AND 2.6
  AND q.civ_1549_flux > 0 AND q.civ_1549_flux_ivar > 0
  AND q.mgii_2796_flux > 0 AND q.mgii_2796_flux_ivar > 0
  AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
  AND q.mgii_2796_flux * SQRT(q.mgii_2796_flux_ivar) > 3
ORDER BY q.random_id"""

data, header = query_desi(sql_main)
if data is None:
    print("FATAL: Main query failed", flush=True)
    sys.exit(1)

ra = data[:, 0]
dec = data[:, 1]
z = data[:, 2]
civ_flux = data[:, 3]
civ_ivar = data[:, 4]
mgii_flux = data[:, 5]
mgii_ivar = data[:, 6]
civ_snr = civ_flux * np.sqrt(np.maximum(civ_ivar, 1e-30))
mgii_snr = mgii_flux * np.sqrt(np.maximum(mgii_ivar, 1e-30))
lum = np.log10(civ_flux + mgii_flux + 1e-20)

print(f"  Got {len(data)} quasars, z={z.min():.2f}-{z.max():.2f}", flush=True)

# Compute galactic coordinates from RA/Dec (for extinction proxy)
def equatorial_to_galactic(ra_deg, dec_deg):
    """Convert equatorial to galactic coordinates"""
    ra_r = np.radians(ra_deg)
    dec_r = np.radians(dec_deg)
    # J2000 galactic pole: RA=192.85948, Dec=27.12825, l_NCP=122.93192
    ra_gp = np.radians(192.85948)
    dec_gp = np.radians(27.12825)
    l_ncp = np.radians(122.93192)
    
    sin_b = np.sin(dec_gp)*np.sin(dec_r) + np.cos(dec_gp)*np.cos(dec_r)*np.cos(ra_r - ra_gp)
    b = np.arcsin(np.clip(sin_b, -1, 1))
    
    y = np.cos(dec_r)*np.sin(ra_r - ra_gp)
    x = np.cos(dec_gp)*np.sin(dec_r) - np.sin(dec_gp)*np.cos(dec_r)*np.cos(ra_r - ra_gp)
    l = l_ncp - np.arctan2(y, x)
    
    return np.degrees(l) % 360, np.degrees(b)

gal_l, gal_b = equatorial_to_galactic(ra, dec)
abs_b = np.abs(gal_b)

print(f"  Galactic |b| range: {abs_b.min():.1f} - {abs_b.max():.1f} deg", flush=True)
print(f"  Median |b|: {np.median(abs_b):.1f} deg", flush=True)

# Compute 3D density
print("\n[2] Computing 3D density...", flush=True)
density = compute_3d_density(ra, dec, z)
log_dens = np.log10(density + 1e-30)

# ============================================================
# TEST 1: Does density correlate with galactic latitude?
# ============================================================
print("\n" + "=" * 70, flush=True)
print("TEST 1: DENSITY vs GALACTIC LATITUDE (extinction proxy)", flush=True)
print("=" * 70, flush=True)

r_dens_b, p_dens_b = spearmanr(log_dens, abs_b)
print(f"  Spearman(log_density, |b|) = {r_dens_b:.4f}, p = {p_dens_b:.2e}", flush=True)

# Per z-bin
z_edges = np.array([1.3, 1.6, 1.9, 2.2, 2.6])
for i in range(len(z_edges)-1):
    zmask = (z >= z_edges[i]) & (z < z_edges[i+1])
    r_zb, p_zb = spearmanr(log_dens[zmask], abs_b[zmask])
    print(f"  z=[{z_edges[i]:.1f},{z_edges[i+1]:.1f}]: r={r_zb:.4f}, p={p_zb:.2e}", flush=True)

# ============================================================
# TEST 2: Does density correlate with SNR?
# ============================================================
print("\n" + "=" * 70, flush=True)
print("TEST 2: DENSITY vs SNR (survey depth proxy)", flush=True)
print("=" * 70, flush=True)

combined_snr = np.sqrt(civ_snr**2 + mgii_snr**2)
r_dens_snr, p_dens_snr = spearmanr(log_dens, combined_snr)
print(f"  Spearman(log_density, combined_SNR) = {r_dens_snr:.4f}, p = {p_dens_snr:.2e}", flush=True)

for i in range(len(z_edges)-1):
    zmask = (z >= z_edges[i]) & (z < z_edges[i+1])
    r_zs, p_zs = spearmanr(log_dens[zmask], combined_snr[zmask])
    print(f"  z=[{z_edges[i]:.1f},{z_edges[i+1]:.1f}]: r={r_zs:.4f}, p={p_zs:.2e}", flush=True)

# ============================================================
# TEST 3: Does density correlate with sky position (RA)?
# ============================================================
print("\n" + "=" * 70, flush=True)
print("TEST 3: DENSITY vs RA (survey footprint proxy)", flush=True)
print("=" * 70, flush=True)

r_dens_ra, p_dens_ra = spearmanr(log_dens, ra)
print(f"  Spearman(log_density, RA) = {r_dens_ra:.4f}, p = {p_dens_ra:.2e}", flush=True)

# ============================================================
# TEST 4: Dense vs Sparse AFTER matching on |b| too
# ============================================================
print("\n" + "=" * 70, flush=True)
print("TEST 4: ENV SPLIT WITH |b| + SNR + LUM MATCHING", flush=True)
print("=" * 70, flush=True)

env_results = []
for i in range(len(z_edges)-1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    zmask = (z >= z_lo) & (z < z_hi)
    if zmask.sum() < 500:
        continue
    
    med = np.median(density[zmask])
    dense_m = zmask & (density >= med)
    sparse_m = zmask & (density < med)
    d_idx, s_idx = np.where(dense_m)[0], np.where(sparse_m)[0]
    
    # Match on: SNR + luminosity + |b| (4D matching)
    matched_d, matched_s = [], []
    n_bins = 6
    snr_bins = np.percentile(combined_snr[zmask], np.linspace(0,100,n_bins+1))
    lum_bins = np.percentile(lum[zmask], np.linspace(0,100,n_bins+1))
    b_bins = np.percentile(abs_b[zmask], np.linspace(0,100,n_bins+1))
    
    for si in range(n_bins):
        for li in range(n_bins):
            for bi in range(n_bins):
                d_in = d_idx[
                    (combined_snr[d_idx]>=snr_bins[si])&(combined_snr[d_idx]<snr_bins[si+1])&
                    (lum[d_idx]>=lum_bins[li])&(lum[d_idx]<lum_bins[li+1])&
                    (abs_b[d_idx]>=b_bins[bi])&(abs_b[d_idx]<b_bins[bi+1])
                ]
                s_in = s_idx[
                    (combined_snr[s_idx]>=snr_bins[si])&(combined_snr[s_idx]<snr_bins[si+1])&
                    (lum[s_idx]>=lum_bins[li])&(lum[s_idx]<lum_bins[li+1])&
                    (abs_b[s_idx]>=b_bins[bi])&(abs_b[s_idx]<b_bins[bi+1])
                ]
                nm = min(len(d_in), len(s_in))
                if nm > 0:
                    matched_d.extend(np.random.choice(d_in, nm, replace=False))
                    matched_s.extend(np.random.choice(s_in, nm, replace=False))
    
    matched_d, matched_s = np.array(matched_d), np.array(matched_s)
    if len(matched_d) < 100:
        print(f"  z=[{z_lo:.1f},{z_hi:.1f}]: matched {len(matched_d)} — skip", flush=True)
        continue
    
    rho_d = spearmanr(civ_flux[matched_d], mgii_flux[matched_d])[0]
    rho_s = spearmanr(civ_flux[matched_s], mgii_flux[matched_s])[0]
    delta = rho_d - rho_s
    
    # Also check that |b| is actually matched
    b_diff = abs(np.median(abs_b[matched_d]) - np.median(abs_b[matched_s]))
    snr_diff = abs(np.median(combined_snr[matched_d]) - np.median(combined_snr[matched_s]))
    
    print(f"  z=[{z_lo:.1f},{z_hi:.1f}] matched={len(matched_d)} "
          f"dense={rho_d:.3f} sparse={rho_s:.3f} Δ={delta:+.4f} "
          f"|b| gap={b_diff:.1f}° SNR gap={snr_diff:.1f}", flush=True)
    
    env_results.append({
        'z_range': f"{z_lo:.1f}-{z_hi:.1f}",
        'n_matched': int(len(matched_d)),
        'dense_rho': round(rho_d, 4),
        'sparse_rho': round(rho_s, 4),
        'delta': round(delta, 4),
        'b_gap': round(b_diff, 2),
        'snr_gap': round(snr_diff, 2),
    })

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70, flush=True)
print("SUMMARY", flush=True)
print("=" * 70, flush=True)

print(f"\nDensity-|b| correlation:  r={r_dens_b:.4f} (p={p_dens_b:.2e})", flush=True)
print(f"Density-SNR correlation: r={r_dens_snr:.4f} (p={p_dens_snr:.2e})", flush=True)
print(f"Density-RA correlation:  r={r_dens_ra:.4f} (p={p_dens_ra:.2e})", flush=True)

contaminated = abs(r_dens_b) > 0.1 or abs(r_dens_snr) > 0.1
if contaminated:
    print("\n⚠️  Density proxy shows moderate correlation with imaging quality", flush=True)
    print("   → Need to verify signal survives after matching", flush=True)
else:
    print("\n✅ Density proxy is CLEAN — weak/no correlation with imaging metrics", flush=True)

if env_results:
    n_pos = sum(1 for r in env_results if r['delta'] > 0)
    mean_d = np.mean([r['delta'] for r in env_results])
    print(f"\nAfter |b|+SNR+lum matching: dense>sparse in {n_pos}/{len(env_results)} bins", flush=True)
    print(f"Mean Δρ = {mean_d:+.4f}", flush=True)
    
    if n_pos >= len(env_results) * 0.6 and mean_d > 0.005:
        print("\n🏆 VERDICT: SIGNAL SURVIVES IMAGING CONTROL — PASS", flush=True)
    elif mean_d > 0:
        print("\n⚠️  VERDICT: SIGNAL WEAKENED BUT POSITIVE — MARGINAL", flush=True)
    else:
        print("\n❌ VERDICT: SIGNAL KILLED BY IMAGING CONTROL — FAIL", flush=True)

output = {
    'test': 'Imaging Depth / E(B-V) Control',
    'density_correlations': {
        'density_vs_abs_b': {'r': round(r_dens_b, 4), 'p': float(f"{p_dens_b:.2e}")},
        'density_vs_snr': {'r': round(r_dens_snr, 4), 'p': float(f"{p_dens_snr:.2e}")},
        'density_vs_ra': {'r': round(r_dens_ra, 4), 'p': float(f"{p_dens_ra:.2e}")},
    },
    'matched_results': env_results,
}
with open(f'{RESULTS_DIR}/imaging_control_results.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {RESULTS_DIR}/imaging_control_results.json", flush=True)
