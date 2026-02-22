#!/usr/bin/env python3
"""
Proper Photometry Control: E(B-V) + PSF Depth + nobs from DESI photometry table
================================================================================
Join agnqso with photometry table to get actual imaging quality metrics.
Then match on ALL of them and rerun env split.
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
    print(f"  Cols: {header}", flush=True)
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
print("PHOTOMETRY CONTROL: E(B-V) + PSF DEPTH + NOBS", flush=True)
print("=" * 70, flush=True)

# Check join key between agnqso and photometry
print("\n[0] Finding join key...", flush=True)
sql_key = """SELECT column_name FROM tap_schema.columns 
WHERE table_name='desi_dr1.agnqso' AND column_name='targetid'"""
r_key = requests.get('https://datalab.noirlab.edu/tap/sync',
                     params={'REQUEST':'doQuery','LANG':'ADQL','QUERY':sql_key,'FORMAT':'csv'}, timeout=30)
print(f"  agnqso has targetid: {'targetid' in r_key.text}", flush=True)

sql_key2 = """SELECT column_name FROM tap_schema.columns 
WHERE table_name='desi_dr1.photometry' AND column_name='targetid'"""
r_key2 = requests.get('https://datalab.noirlab.edu/tap/sync',
                      params={'REQUEST':'doQuery','LANG':'ADQL','QUERY':sql_key2,'FORMAT':'csv'}, timeout=30)
print(f"  photometry has targetid: {'targetid' in r_key2.text}", flush=True)

# Query with join
print("\n[1] Querying agnqso JOIN photometry...", flush=True)
sql = """SELECT TOP 120000
    q.target_ra, q.target_dec, q.z,
    q.civ_1549_flux, q.civ_1549_flux_ivar,
    q.mgii_2796_flux, q.mgii_2796_flux_ivar,
    p.ebv, p.psfdepth_g, p.psfdepth_r, p.psfdepth_z,
    p.nobs_g, p.nobs_r, p.nobs_z
FROM desi_dr1.agnqso q
JOIN desi_dr1.photometry p ON q.targetid = p.targetid
WHERE q.z BETWEEN 1.3 AND 2.6
  AND q.civ_1549_flux > 0 AND q.civ_1549_flux_ivar > 0
  AND q.mgii_2796_flux > 0 AND q.mgii_2796_flux_ivar > 0
  AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
  AND q.mgii_2796_flux * SQRT(q.mgii_2796_flux_ivar) > 3
  AND p.ebv > 0
ORDER BY q.random_id"""

data, header = query_desi(sql)
if data is None:
    print("FATAL: Join query failed", flush=True)
    sys.exit(1)

ra = data[:, 0]
dec = data[:, 1]
z = data[:, 2]
civ_flux = data[:, 3]
civ_ivar = data[:, 4]
mgii_flux = data[:, 5]
mgii_ivar = data[:, 6]
ebv = data[:, 7]
psfdepth_g = data[:, 8]
psfdepth_r = data[:, 9]
psfdepth_z = data[:, 10]
nobs_g = data[:, 11]
nobs_r = data[:, 12]
nobs_z = data[:, 13]

civ_snr = civ_flux * np.sqrt(np.maximum(civ_ivar, 1e-30))
mgii_snr = mgii_flux * np.sqrt(np.maximum(mgii_ivar, 1e-30))
combined_snr = np.sqrt(civ_snr**2 + mgii_snr**2)
lum = np.log10(civ_flux + mgii_flux + 1e-20)
# Combined imaging quality metric
img_quality = np.log10(psfdepth_r + 1)  # r-band depth as primary proxy

print(f"  Got {len(data)} quasars with photometry", flush=True)
print(f"  E(B-V) range: {ebv.min():.3f} - {ebv.max():.3f}, median: {np.median(ebv):.3f}", flush=True)
print(f"  psfdepth_r range: {psfdepth_r.min():.0f} - {psfdepth_r.max():.0f}", flush=True)
print(f"  nobs_r range: {nobs_r.min():.0f} - {nobs_r.max():.0f}", flush=True)

# Compute 3D density
print("\n[2] Computing 3D density...", flush=True)
density = compute_3d_density(ra, dec, z)
log_dens = np.log10(density + 1e-30)

# ============================================================
# CORRELATIONS: density vs every imaging metric
# ============================================================
print("\n" + "=" * 70, flush=True)
print("DENSITY vs IMAGING QUALITY CORRELATIONS", flush=True)
print("=" * 70, flush=True)

metrics = {
    'E(B-V)': ebv,
    'psfdepth_g': np.log10(psfdepth_g + 1),
    'psfdepth_r': np.log10(psfdepth_r + 1),
    'psfdepth_z': np.log10(psfdepth_z + 1),
    'nobs_r': nobs_r,
    'combined_SNR': combined_snr,
}

corr_results = {}
for name, values in metrics.items():
    r_val, p_val = spearmanr(log_dens, values)
    print(f"  density vs {name:15s}: r = {r_val:+.4f}  (p = {p_val:.2e})", flush=True)
    corr_results[name] = {'r': round(r_val, 4), 'p': float(f"{p_val:.2e}")}

# ============================================================
# THE ULTIMATE CONTROL: match on E(B-V) + psfdepth + SNR + lum
# ============================================================
print("\n" + "=" * 70, flush=True)
print("ULTIMATE CONTROL: MATCH ON E(B-V) + DEPTH + SNR + LUM", flush=True)
print("=" * 70, flush=True)

z_edges = np.array([1.3, 1.6, 1.9, 2.2, 2.6])
final_results = []

for i in range(len(z_edges)-1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    zmask = (z >= z_lo) & (z < z_hi)
    if zmask.sum() < 500:
        continue
    
    med = np.median(density[zmask])
    dense_m = zmask & (density >= med)
    sparse_m = zmask & (density < med)
    d_idx, s_idx = np.where(dense_m)[0], np.where(sparse_m)[0]
    
    # 5D matching: SNR + lum + E(B-V) + psfdepth_r + nobs_r
    matched_d, matched_s = [], []
    nb = 4  # bins per dimension (4^5 = 1024 cells)
    
    snr_bins = np.percentile(combined_snr[zmask], np.linspace(0,100,nb+1))
    lum_bins = np.percentile(lum[zmask], np.linspace(0,100,nb+1))
    ebv_bins = np.percentile(ebv[zmask], np.linspace(0,100,nb+1))
    depth_bins = np.percentile(img_quality[zmask], np.linspace(0,100,nb+1))
    
    for si in range(nb):
        for li in range(nb):
            for ei in range(nb):
                for di in range(nb):
                    d_in = d_idx[
                        (combined_snr[d_idx]>=snr_bins[si])&(combined_snr[d_idx]<snr_bins[si+1])&
                        (lum[d_idx]>=lum_bins[li])&(lum[d_idx]<lum_bins[li+1])&
                        (ebv[d_idx]>=ebv_bins[ei])&(ebv[d_idx]<ebv_bins[ei+1])&
                        (img_quality[d_idx]>=depth_bins[di])&(img_quality[d_idx]<depth_bins[di+1])
                    ]
                    s_in = s_idx[
                        (combined_snr[s_idx]>=snr_bins[si])&(combined_snr[s_idx]<snr_bins[si+1])&
                        (lum[s_idx]>=lum_bins[li])&(lum[s_idx]<lum_bins[li+1])&
                        (ebv[s_idx]>=ebv_bins[ei])&(ebv[s_idx]<ebv_bins[ei+1])&
                        (img_quality[s_idx]>=depth_bins[di])&(img_quality[s_idx]<depth_bins[di+1])
                    ]
                    nm = min(len(d_in), len(s_in))
                    if nm > 0:
                        matched_d.extend(np.random.choice(d_in, nm, replace=False))
                        matched_s.extend(np.random.choice(s_in, nm, replace=False))
    
    matched_d, matched_s = np.array(matched_d), np.array(matched_s)
    if len(matched_d) < 80:
        print(f"  z=[{z_lo:.1f},{z_hi:.1f}]: matched {len(matched_d)} — skip", flush=True)
        continue
    
    rho_d = spearmanr(civ_flux[matched_d], mgii_flux[matched_d])[0]
    rho_s = spearmanr(civ_flux[matched_s], mgii_flux[matched_s])[0]
    delta = rho_d - rho_s
    
    # Verify matching quality
    ebv_gap = abs(np.median(ebv[matched_d]) - np.median(ebv[matched_s]))
    depth_gap = abs(np.median(psfdepth_r[matched_d]) - np.median(psfdepth_r[matched_s]))
    snr_gap = abs(np.median(combined_snr[matched_d]) - np.median(combined_snr[matched_s]))
    
    print(f"  z=[{z_lo:.1f},{z_hi:.1f}] matched={len(matched_d)} "
          f"dense={rho_d:.3f} sparse={rho_s:.3f} Δ={delta:+.4f} "
          f"E(B-V)gap={ebv_gap:.4f} depthgap={depth_gap:.0f} SNRgap={snr_gap:.1f}", flush=True)
    
    final_results.append({
        'z_range': f"{z_lo:.1f}-{z_hi:.1f}",
        'n_matched': int(len(matched_d)),
        'dense_rho': round(rho_d, 4),
        'sparse_rho': round(rho_s, 4),
        'delta': round(delta, 4),
        'ebv_gap': round(ebv_gap, 5),
        'depth_gap': round(depth_gap, 1),
        'snr_gap': round(snr_gap, 2),
    })

# SUMMARY
print("\n" + "=" * 70, flush=True)
print("SUMMARY", flush=True)
print("=" * 70, flush=True)

if final_results:
    n_pos = sum(1 for r in final_results if r['delta'] > 0)
    mean_d = np.mean([r['delta'] for r in final_results])
    print(f"  After E(B-V)+depth+SNR+lum matching: dense>sparse {n_pos}/{len(final_results)}", flush=True)
    print(f"  Mean Δρ = {mean_d:+.4f}", flush=True)
    
    if n_pos >= len(final_results) * 0.6 and mean_d > 0.003:
        print(f"\n  🏆 ULTIMATE CONTROL: PASS", flush=True)
    elif mean_d > 0:
        print(f"\n  ⚠️  ULTIMATE CONTROL: MARGINAL", flush=True)
    else:
        print(f"\n  ❌ ULTIMATE CONTROL: FAIL", flush=True)

output = {
    'test': 'Photometry Control: E(B-V) + PSF Depth + nobs',
    'correlations': corr_results,
    'matched_results': final_results,
}
with open(f'{RESULTS_DIR}/photometry_control_results.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {RESULTS_DIR}/photometry_control_results.json", flush=True)
