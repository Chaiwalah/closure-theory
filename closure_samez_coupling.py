#!/usr/bin/env python3
"""
Same-z Coupling Comparison
============================
Fix the confound: MgII/Hβ was at z=0.4-1.0 while CIV/MgII at z=1.3-2.6.
To isolate coupling strength from z-range, compare pairs at MATCHED z-ranges.

At z=1.3-1.9 we can measure:
- CIV/MgII (moderate coupling)
- CIV/NeV (weak coupling) 
- CIV/Hβ_broad (cross-species, if available)

If ordering holds at same z → coupling drives it, not distance.
"""
import numpy as np
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid
import requests, json, os, warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = 'results_samez_coupling'
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
print("SAME-Z COUPLING COMPARISON", flush=True)
print("=" * 70, flush=True)

# Query all objects at z=1.3-1.9 with CIV + MgII + NeV (all three)
print("\n[1] Querying objects with ALL THREE lines at z=1.3-1.9...", flush=True)
sql = """SELECT TOP 200000
    q.target_ra, q.target_dec, q.z,
    q.civ_1549_flux, q.civ_1549_flux_ivar,
    q.mgii_2796_flux, q.mgii_2796_flux_ivar,
    q.nev_3426_flux, q.nev_3426_flux_ivar
FROM desi_dr1.agnqso q
WHERE q.z BETWEEN 1.3 AND 1.9
  AND q.civ_1549_flux > 0 AND q.civ_1549_flux_ivar > 0
  AND q.mgii_2796_flux > 0 AND q.mgii_2796_flux_ivar > 0
  AND q.nev_3426_flux > 0 AND q.nev_3426_flux_ivar > 0
  AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
  AND q.mgii_2796_flux * SQRT(q.mgii_2796_flux_ivar) > 3
  AND q.nev_3426_flux * SQRT(q.nev_3426_flux_ivar) > 2
ORDER BY q.random_id"""

data = query_desi(sql)
if data is None or len(data) < 500:
    print("Not enough objects with all 3 lines. Trying pairs separately.", flush=True)
    # Fall back: query CIV/MgII and CIV/NeV separately at same z, use same density
    
    print("\n[1b] CIV/MgII at z=1.3-1.9...", flush=True)
    sql_cm = """SELECT TOP 150000
        q.target_ra, q.target_dec, q.z,
        q.civ_1549_flux, q.civ_1549_flux_ivar,
        q.mgii_2796_flux, q.mgii_2796_flux_ivar
    FROM desi_dr1.agnqso q
    WHERE q.z BETWEEN 1.3 AND 1.9
      AND q.civ_1549_flux > 0 AND q.civ_1549_flux_ivar > 0
      AND q.mgii_2796_flux > 0 AND q.mgii_2796_flux_ivar > 0
      AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
      AND q.mgii_2796_flux * SQRT(q.mgii_2796_flux_ivar) > 3
    ORDER BY q.random_id"""
    d_cm = query_desi(sql_cm)
    
    print("\n[1c] CIV/NeV at z=1.3-1.9...", flush=True)
    sql_cn = """SELECT TOP 150000
        q.target_ra, q.target_dec, q.z,
        q.civ_1549_flux, q.civ_1549_flux_ivar,
        q.nev_3426_flux, q.nev_3426_flux_ivar
    FROM desi_dr1.agnqso q
    WHERE q.z BETWEEN 1.3 AND 1.9
      AND q.civ_1549_flux > 0 AND q.civ_1549_flux_ivar > 0
      AND q.nev_3426_flux > 0 AND q.nev_3426_flux_ivar > 0
      AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
      AND q.nev_3426_flux * SQRT(q.nev_3426_flux_ivar) > 2
    ORDER BY q.random_id"""
    d_cn = query_desi(sql_cn)
    
    use_triple = False
else:
    print(f"  Got {len(data)} with all 3 lines", flush=True)
    use_triple = True

z_edges = np.array([1.3, 1.5, 1.7, 1.9])
results = []

def run_env_test(ra, dec, z, f1, f2, snr1, snr2, label, z_edges):
    density = compute_3d_density(ra, dec, z)
    lum = np.log10(f1 + f2 + 1e-20)
    combined_snr = np.sqrt(snr1**2 + snr2**2)
    
    pair_results = []
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        zmask = (z >= z_lo) & (z < z_hi)
        if zmask.sum() < 300:
            continue
        med = np.median(density[zmask])
        dense_m = zmask & (density >= med)
        sparse_m = zmask & (density < med)
        d_idx, s_idx = np.where(dense_m)[0], np.where(sparse_m)[0]
        
        # Match SNR + lum
        matched_d, matched_s = [], []
        nb = 8
        snr_bins = np.percentile(combined_snr[zmask], np.linspace(0,100,nb+1))
        lum_bins = np.percentile(lum[zmask], np.linspace(0,100,nb+1))
        for si in range(nb):
            for li in range(nb):
                d_in = d_idx[(combined_snr[d_idx]>=snr_bins[si])&(combined_snr[d_idx]<snr_bins[si+1])&
                             (lum[d_idx]>=lum_bins[li])&(lum[d_idx]<lum_bins[li+1])]
                s_in = s_idx[(combined_snr[s_idx]>=snr_bins[si])&(combined_snr[s_idx]<snr_bins[si+1])&
                             (lum[s_idx]>=lum_bins[li])&(lum[s_idx]<lum_bins[li+1])]
                nm = min(len(d_in), len(s_in))
                if nm > 0:
                    matched_d.extend(np.random.choice(d_in, nm, replace=False))
                    matched_s.extend(np.random.choice(s_in, nm, replace=False))
        
        if len(matched_d) < 80:
            continue
        matched_d, matched_s = np.array(matched_d), np.array(matched_s)
        rho_d = spearmanr(f1[matched_d], f2[matched_d])[0]
        rho_s = spearmanr(f1[matched_s], f2[matched_s])[0]
        delta = rho_d - rho_s
        print(f"    z=[{z_lo:.1f},{z_hi:.1f}] N={zmask.sum()} matched={len(matched_d)} Δ={delta:+.4f}", flush=True)
        pair_results.append({'z': (z_lo+z_hi)/2, 'delta': delta, 'n': len(matched_d)})
    
    mean_d = np.mean([r['delta'] for r in pair_results]) if pair_results else 0
    n_pos = sum(1 for r in pair_results if r['delta'] > 0)
    return {'pair': label, 'mean_delta': round(mean_d, 4), 'n_pos': n_pos, 
            'n_bins': len(pair_results), 'bins': pair_results}

if use_triple:
    ra, dec, z = data[:,0], data[:,1], data[:,2]
    civ_f, civ_iv = data[:,3], data[:,4]
    mgii_f, mgii_iv = data[:,5], data[:,6]
    nev_f, nev_iv = data[:,7], data[:,8]
    
    civ_snr = civ_f * np.sqrt(np.maximum(civ_iv, 1e-30))
    mgii_snr = mgii_f * np.sqrt(np.maximum(mgii_iv, 1e-30))
    nev_snr = nev_f * np.sqrt(np.maximum(nev_iv, 1e-30))
    
    print("\n  CIV/MgII (same objects, same z):", flush=True)
    r_cm = run_env_test(ra, dec, z, civ_f, mgii_f, civ_snr, mgii_snr, "CIV/MgII", z_edges)
    
    print("\n  CIV/NeV (same objects, same z):", flush=True)
    r_cn = run_env_test(ra, dec, z, civ_f, nev_f, civ_snr, nev_snr, "CIV/NeV", z_edges)
    
    results = [r_cm, r_cn]
else:
    if d_cm is not None:
        print("\n  CIV/MgII at z=1.3-1.9:", flush=True)
        snr1 = d_cm[:,3]*np.sqrt(np.maximum(d_cm[:,4],1e-30))
        snr2 = d_cm[:,5]*np.sqrt(np.maximum(d_cm[:,6],1e-30))
        r_cm = run_env_test(d_cm[:,0], d_cm[:,1], d_cm[:,2], d_cm[:,3], d_cm[:,5],
                           snr1, snr2, "CIV/MgII", z_edges)
        results.append(r_cm)
    
    if d_cn is not None:
        print("\n  CIV/NeV at z=1.3-1.9:", flush=True)
        snr1 = d_cn[:,3]*np.sqrt(np.maximum(d_cn[:,4],1e-30))
        snr2 = d_cn[:,5]*np.sqrt(np.maximum(d_cn[:,6],1e-30))
        r_cn = run_env_test(d_cn[:,0], d_cn[:,1], d_cn[:,2], d_cn[:,3], d_cn[:,5],
                           snr1, snr2, "CIV/NeV", z_edges)
        results.append(r_cn)

print("\n" + "=" * 70, flush=True)
print("SAME-Z COMPARISON", flush=True)
print("=" * 70, flush=True)

for r in results:
    print(f"  {r['pair']}: mean Δρ = {r['mean_delta']:+.4f} ({r['n_pos']}/{r['n_bins']} positive)", flush=True)

if len(results) >= 2:
    cm_delta = abs(results[0]['mean_delta'])
    cn_delta = abs(results[1]['mean_delta'])
    if cn_delta > cm_delta:
        print(f"\n  ✅ ORDERING CONFIRMED AT SAME Z: |CIV/NeV| ({cn_delta:.4f}) > |CIV/MgII| ({cm_delta:.4f})", flush=True)
        print(f"  Weaker coupling → stronger environment modulation, independent of distance", flush=True)
    else:
        print(f"\n  ❌ ORDERING NOT CONFIRMED: |CIV/NeV| ({cn_delta:.4f}) vs |CIV/MgII| ({cm_delta:.4f})", flush=True)

output = {'test': 'Same-z Coupling Comparison', 'z_range': '1.3-1.9', 'results': results}
with open(f'{RESULTS_DIR}/samez_results.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {RESULTS_DIR}/samez_results.json", flush=True)
