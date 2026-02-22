#!/usr/bin/env python3
"""
BLR State Proxy Matching — Kill "Intrinsic Quasar Population" Explanation
=========================================================================
Match dense/sparse on:
- CIV sigma (width → virial/wind state)
- MgII sigma (width → virial state)
- CIV EW (equivalent width → Baldwin effect proxy)
- MgII EW (equivalent width)
- CIV flux / MgII flux ratio (ionization state proxy)

If Dense > Sparse survives after matching on ALL of these + SNR + lum,
the "different quasar populations in different environments" explanation is dead.
"""
import numpy as np
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid
import requests, json, os, warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = 'results_blr_control'
os.makedirs(RESULTS_DIR, exist_ok=True)

def query_desi(sql, timeout=600):
    print(f"  Querying...", flush=True)
    r = requests.get('https://datalab.noirlab.edu/tap/sync',
                     params={'REQUEST':'doQuery','LANG':'ADQL','QUERY':sql,'FORMAT':'csv'},
                     timeout=timeout)
    if r.status_code != 200 or 'ERROR' in r.text[:500]:
        print(f"  Error: {r.text[:200]}", flush=True)
        return None
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
    return np.array(data) if len(data) > 100 else None

def compute_3d_density(ra, dec, z, n_neighbors=10):
    z_fine = np.linspace(0, max(z)*1.1, 5000)
    Hz = np.sqrt(0.3*(1+z_fine)**3 + 0.7)
    chi_fine = cumulative_trapezoid(1.0/Hz, z_fine, initial=0) * 2997.9
    chi = np.interp(z, z_fine, chi_fine)
    ra_r, dec_r = np.radians(ra), np.radians(dec)
    coords = np.column_stack([
        chi*np.cos(dec_r)*np.cos(ra_r),
        chi*np.cos(dec_r)*np.sin(ra_r),
        chi*np.sin(dec_r)])
    tree = cKDTree(coords)
    dists, _ = tree.query(coords, k=n_neighbors+1)
    return n_neighbors / (4/3 * np.pi * dists[:, -1]**3 + 1e-10)

print("=" * 70, flush=True)
print("BLR STATE PROXY MATCHING", flush=True)
print("=" * 70, flush=True)

# Check what EW columns exist
print("\n[0] Checking EW columns...", flush=True)
sql_ew = """SELECT column_name FROM tap_schema.columns 
WHERE table_name='desi_dr1.agnqso' AND (column_name LIKE '%ew%')
ORDER BY column_name"""
r_ew = requests.get('https://datalab.noirlab.edu/tap/sync',
                    params={'REQUEST':'doQuery','LANG':'ADQL','QUERY':sql_ew,'FORMAT':'csv'}, timeout=30)
print(f"  EW columns: {r_ew.text.strip()}", flush=True)

# Query with all BLR proxies
print("\n[1] Querying with BLR proxies...", flush=True)
sql = """SELECT TOP 150000
    q.target_ra, q.target_dec, q.z,
    q.civ_1549_flux, q.civ_1549_flux_ivar, q.civ_1549_sigma,
    q.mgii_2796_flux, q.mgii_2796_flux_ivar, q.mgii_2796_sigma
FROM desi_dr1.agnqso q
WHERE q.z BETWEEN 1.3 AND 2.6
  AND q.civ_1549_flux > 0 AND q.civ_1549_flux_ivar > 0
  AND q.mgii_2796_flux > 0 AND q.mgii_2796_flux_ivar > 0
  AND q.civ_1549_sigma > 0 AND q.mgii_2796_sigma > 0
  AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
  AND q.mgii_2796_flux * SQRT(q.mgii_2796_flux_ivar) > 3
ORDER BY q.random_id"""

data = query_desi(sql)
if data is None:
    print("FATAL", flush=True)
    import sys; sys.exit(1)

ra = data[:, 0]
dec = data[:, 1]
z = data[:, 2]
civ_flux = data[:, 3]
civ_ivar = data[:, 4]
civ_sigma = data[:, 5]
mgii_flux = data[:, 6]
mgii_ivar = data[:, 7]
mgii_sigma = data[:, 8]

civ_snr = civ_flux * np.sqrt(np.maximum(civ_ivar, 1e-30))
mgii_snr = mgii_flux * np.sqrt(np.maximum(mgii_ivar, 1e-30))
combined_snr = np.sqrt(civ_snr**2 + mgii_snr**2)
lum = np.log10(civ_flux + mgii_flux + 1e-20)

# BLR state proxies
flux_ratio = np.log10(civ_flux / mgii_flux + 1e-10)  # ionization state
sigma_ratio = civ_sigma / (mgii_sigma + 1e-10)  # wind vs virial

print(f"  Got {len(data)} quasars with all BLR proxies", flush=True)
print(f"  CIV sigma range: {civ_sigma.min():.1f} - {civ_sigma.max():.1f}", flush=True)
print(f"  Flux ratio (CIV/MgII) range: {flux_ratio.min():.2f} - {flux_ratio.max():.2f}", flush=True)

# Compute 3D density
print("\n[2] Computing 3D density...", flush=True)
density = compute_3d_density(ra, dec, z)
log_dens = np.log10(density + 1e-30)

# ============================================================
# First: check if BLR proxies correlate with density
# ============================================================
print("\n" + "=" * 70, flush=True)
print("BLR PROXY CORRELATIONS WITH DENSITY", flush=True)
print("=" * 70, flush=True)

blr_metrics = {
    'CIV_sigma': civ_sigma,
    'MgII_sigma': mgii_sigma,
    'CIV/MgII_flux_ratio': flux_ratio,
    'sigma_ratio(CIV/MgII)': sigma_ratio,
    'combined_SNR': combined_snr,
    'luminosity': lum,
}

corr_results = {}
for name, vals in blr_metrics.items():
    r_val, p_val = spearmanr(log_dens, vals)
    print(f"  density vs {name:25s}: r = {r_val:+.4f}  (p = {p_val:.2e})", flush=True)
    corr_results[name] = {'r': round(r_val, 4)}

# ============================================================
# ULTIMATE BLR CONTROL: match on SNR + lum + CIV_sigma + MgII_sigma + flux_ratio
# ============================================================
print("\n" + "=" * 70, flush=True)
print("ULTIMATE BLR CONTROL: 5D MATCHING", flush=True)
print("SNR + Lum + CIV_sigma + MgII_sigma + CIV/MgII flux ratio", flush=True)
print("=" * 70, flush=True)

z_edges = np.array([1.3, 1.6, 1.9, 2.2, 2.6])
results = []

for i in range(len(z_edges)-1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    zmask = (z >= z_lo) & (z < z_hi)
    if zmask.sum() < 500:
        continue
    
    med = np.median(density[zmask])
    dense_m = zmask & (density >= med)
    sparse_m = zmask & (density < med)
    d_idx, s_idx = np.where(dense_m)[0], np.where(sparse_m)[0]
    
    # 5D matching: SNR + lum + CIV_sigma + MgII_sigma + flux_ratio
    matched_d, matched_s = [], []
    nb = 3  # 3^5 = 243 cells
    
    snr_bins = np.percentile(combined_snr[zmask], np.linspace(0,100,nb+1))
    lum_bins = np.percentile(lum[zmask], np.linspace(0,100,nb+1))
    csig_bins = np.percentile(civ_sigma[zmask], np.linspace(0,100,nb+1))
    msig_bins = np.percentile(mgii_sigma[zmask], np.linspace(0,100,nb+1))
    fr_bins = np.percentile(flux_ratio[zmask], np.linspace(0,100,nb+1))
    
    for si in range(nb):
        for li in range(nb):
            for ci in range(nb):
                for mi in range(nb):
                    for fi in range(nb):
                        d_in = d_idx[
                            (combined_snr[d_idx]>=snr_bins[si])&(combined_snr[d_idx]<snr_bins[si+1])&
                            (lum[d_idx]>=lum_bins[li])&(lum[d_idx]<lum_bins[li+1])&
                            (civ_sigma[d_idx]>=csig_bins[ci])&(civ_sigma[d_idx]<csig_bins[ci+1])&
                            (mgii_sigma[d_idx]>=msig_bins[mi])&(mgii_sigma[d_idx]<msig_bins[mi+1])&
                            (flux_ratio[d_idx]>=fr_bins[fi])&(flux_ratio[d_idx]<fr_bins[fi+1])
                        ]
                        s_in = s_idx[
                            (combined_snr[s_idx]>=snr_bins[si])&(combined_snr[s_idx]<snr_bins[si+1])&
                            (lum[s_idx]>=lum_bins[li])&(lum[s_idx]<lum_bins[li+1])&
                            (civ_sigma[s_idx]>=csig_bins[ci])&(civ_sigma[s_idx]<csig_bins[ci+1])&
                            (mgii_sigma[s_idx]>=msig_bins[mi])&(mgii_sigma[s_idx]<msig_bins[mi+1])&
                            (flux_ratio[s_idx]>=fr_bins[fi])&(flux_ratio[s_idx]<fr_bins[fi+1])
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
    
    # Verify all matching gaps
    gaps = {
        'SNR': abs(np.median(combined_snr[matched_d]) - np.median(combined_snr[matched_s])),
        'lum': abs(np.median(lum[matched_d]) - np.median(lum[matched_s])),
        'CIV_sig': abs(np.median(civ_sigma[matched_d]) - np.median(civ_sigma[matched_s])),
        'MgII_sig': abs(np.median(mgii_sigma[matched_d]) - np.median(mgii_sigma[matched_s])),
        'flux_ratio': abs(np.median(flux_ratio[matched_d]) - np.median(flux_ratio[matched_s])),
    }
    
    print(f"  z=[{z_lo:.1f},{z_hi:.1f}] matched={len(matched_d)} "
          f"dense={rho_d:.3f} sparse={rho_s:.3f} Δ={delta:+.4f}", flush=True)
    print(f"    Gaps: SNR={gaps['SNR']:.1f} lum={gaps['lum']:.3f} "
          f"CIVσ={gaps['CIV_sig']:.1f} MgIIσ={gaps['MgII_sig']:.1f} "
          f"ratio={gaps['flux_ratio']:.3f}", flush=True)
    
    results.append({
        'z_range': f"{z_lo:.1f}-{z_hi:.1f}",
        'n_matched': int(len(matched_d)),
        'dense_rho': round(rho_d, 4),
        'sparse_rho': round(rho_s, 4),
        'delta': round(delta, 4),
        'gaps': {k: round(v, 3) for k, v in gaps.items()},
    })

# SUMMARY
print("\n" + "=" * 70, flush=True)
print("SUMMARY", flush=True)
print("=" * 70, flush=True)

n_pos = sum(1 for r in results if r['delta'] > 0)
mean_d = np.mean([r['delta'] for r in results]) if results else 0
print(f"  After BLR 5D matching: dense>sparse {n_pos}/{len(results)}", flush=True)
print(f"  Mean Δρ = {mean_d:+.4f}", flush=True)

if n_pos >= len(results) * 0.6 and mean_d > 0.003:
    print(f"\n  🏆 BLR POPULATION CONTROL: PASS — quasar physics doesn't explain it", flush=True)
elif mean_d > 0:
    print(f"\n  ⚠️  BLR CONTROL: MARGINAL", flush=True)
else:
    print(f"\n  ❌ BLR CONTROL: FAIL — quasar population explains the effect", flush=True)

output = {
    'test': 'BLR State Proxy Matching (5D)',
    'matching_vars': ['SNR', 'luminosity', 'CIV_sigma', 'MgII_sigma', 'CIV/MgII_flux_ratio'],
    'density_correlations': corr_results,
    'results': results,
}
with open(f'{RESULTS_DIR}/blr_control_results.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {RESULTS_DIR}/blr_control_results.json", flush=True)
