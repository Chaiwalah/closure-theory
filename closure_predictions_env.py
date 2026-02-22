#!/usr/bin/env python3
"""
Pre-Registered Prediction Test: Environment Δρ for 3 new line pairs
===================================================================
Tests whether environment modulation scales with coupling looseness:
|Δρ(HeII/CIV)| < |Δρ(CIV/MgII)| < |Δρ(MgII/Hβ)| < |Δρ(CIV/NeV)|
"""
import numpy as np
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid
import requests, json, os, sys, warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = 'results_predictions_env'
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
    if len(lines) < 50:
        return None
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
    x = chi * np.cos(dec_r) * np.cos(ra_r)
    y = chi * np.cos(dec_r) * np.sin(ra_r)
    zc = chi * np.sin(dec_r)
    tree = cKDTree(np.column_stack([x, y, zc]))
    dists, _ = tree.query(np.column_stack([x, y, zc]), k=n_neighbors+1)
    return n_neighbors / (4/3 * np.pi * dists[:, -1]**3 + 1e-10)

def env_delta_test(ra, dec, z, flux1, flux2, snr1, snr2, lum, label, z_edges):
    """Run environment split with SNR+lum matching, return mean Δρ"""
    print(f"\n  Testing {label}...", flush=True)
    density = compute_3d_density(ra, dec, z)
    
    results = []
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        zmask = (z >= z_lo) & (z < z_hi)
        n = zmask.sum()
        if n < 300:
            continue
        
        med = np.median(density[zmask])
        dense_m = zmask & (density >= med)
        sparse_m = zmask & (density < med)
        
        # SNR + lum matching
        snr_c = np.sqrt(snr1**2 + snr2**2)
        d_idx, s_idx = np.where(dense_m)[0], np.where(sparse_m)[0]
        matched_d, matched_s = [], []
        n_bins = 8
        snr_bins = np.percentile(snr_c[zmask], np.linspace(0,100,n_bins+1))
        lum_bins = np.percentile(lum[zmask], np.linspace(0,100,n_bins+1))
        for si in range(n_bins):
            for li in range(n_bins):
                d_in = d_idx[(snr_c[d_idx]>=snr_bins[si])&(snr_c[d_idx]<snr_bins[si+1])&
                             (lum[d_idx]>=lum_bins[li])&(lum[d_idx]<lum_bins[li+1])]
                s_in = s_idx[(snr_c[s_idx]>=snr_bins[si])&(snr_c[s_idx]<snr_bins[si+1])&
                             (lum[s_idx]>=lum_bins[li])&(lum[s_idx]<lum_bins[li+1])]
                nm = min(len(d_in), len(s_in))
                if nm > 0:
                    matched_d.extend(np.random.choice(d_in, nm, replace=False))
                    matched_s.extend(np.random.choice(s_in, nm, replace=False))
        
        if len(matched_d) < 80:
            continue
        matched_d, matched_s = np.array(matched_d), np.array(matched_s)
        
        rho_d = spearmanr(flux1[matched_d], flux2[matched_d])[0]
        rho_s = spearmanr(flux1[matched_s], flux2[matched_s])[0]
        delta = rho_d - rho_s
        
        print(f"    z=[{z_lo:.1f},{z_hi:.1f}] N={n} matched={len(matched_d)} "
              f"Δ={delta:+.4f}", flush=True)
        results.append({'z': (z_lo+z_hi)/2, 'delta': delta, 'n': len(matched_d)})
    
    mean_delta = np.mean([r['delta'] for r in results]) if results else 0
    n_positive = sum(1 for r in results if r['delta'] > 0)
    print(f"  {label}: mean Δρ = {mean_delta:+.4f}, positive {n_positive}/{len(results)}", flush=True)
    return {'pair': label, 'mean_delta': round(mean_delta, 4), 
            'n_positive': n_positive, 'n_bins': len(results), 'bins': results}

# ============================================================
print("=" * 70, flush=True)
print("PRE-REGISTERED PREDICTION TEST", flush=True)
print("=" * 70, flush=True)

# PAIR 0 (baseline): CIV/MgII — already tested, re-run for consistency
print("\n[0] CIV/MgII (baseline)...", flush=True)
sql0 = """SELECT TOP 150000 q.target_ra, q.target_dec, q.z,
    q.civ_1549_flux, q.civ_1549_flux_ivar,
    q.mgii_2796_flux, q.mgii_2796_flux_ivar
FROM desi_dr1.agnqso q
WHERE q.z BETWEEN 1.3 AND 2.6
  AND q.civ_1549_flux > 0 AND q.civ_1549_flux_ivar > 0
  AND q.mgii_2796_flux > 0 AND q.mgii_2796_flux_ivar > 0
  AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
  AND q.mgii_2796_flux * SQRT(q.mgii_2796_flux_ivar) > 3
ORDER BY q.random_id"""

d0 = query_desi(sql0)
if d0 is not None:
    z_edges_civ = np.array([1.3, 1.6, 1.9, 2.2, 2.6])
    snr0_1 = d0[:,3]*np.sqrt(np.maximum(d0[:,4],1e-30))
    snr0_2 = d0[:,5]*np.sqrt(np.maximum(d0[:,6],1e-30))
    lum0 = np.log10(d0[:,3]+d0[:,5]+1e-20)
    r0 = env_delta_test(d0[:,0], d0[:,1], d0[:,2], d0[:,3], d0[:,5],
                        snr0_1, snr0_2, lum0, "CIV/MgII", z_edges_civ)

# PAIR 1: HeII(4686) / CIV(1549) — predict SMALLER Δρ
print("\n[1] HeII/CIV (predict small Δρ)...", flush=True)
sql1 = """SELECT TOP 150000 q.target_ra, q.target_dec, q.z,
    q.heii_4686_flux, q.heii_4686_flux_ivar,
    q.civ_1549_flux, q.civ_1549_flux_ivar
FROM desi_dr1.agnqso q
WHERE q.z BETWEEN 0.5 AND 2.0
  AND q.heii_4686_flux > 0 AND q.heii_4686_flux_ivar > 0
  AND q.civ_1549_flux > 0 AND q.civ_1549_flux_ivar > 0
  AND q.heii_4686_flux * SQRT(q.heii_4686_flux_ivar) > 2
  AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
ORDER BY q.random_id"""

d1 = query_desi(sql1)
if d1 is not None:
    z_edges_heii = np.array([0.5, 0.8, 1.1, 1.4, 1.7, 2.0])
    snr1_1 = d1[:,3]*np.sqrt(np.maximum(d1[:,4],1e-30))
    snr1_2 = d1[:,5]*np.sqrt(np.maximum(d1[:,6],1e-30))
    lum1 = np.log10(d1[:,3]+d1[:,5]+1e-20)
    r1 = env_delta_test(d1[:,0], d1[:,1], d1[:,2], d1[:,3], d1[:,5],
                        snr1_1, snr1_2, lum1, "HeII/CIV", z_edges_heii)
else:
    print("  HeII/CIV: query failed or too few objects", flush=True)
    r1 = None

# PAIR 2: MgII(2796) / Hβ_broad — predict LARGER Δρ
print("\n[2] MgII/Hβ_broad (predict large Δρ)...", flush=True)
sql2 = """SELECT TOP 150000 q.target_ra, q.target_dec, q.z,
    q.mgii_2796_flux, q.mgii_2796_flux_ivar,
    q.hbeta_broad_flux, q.hbeta_broad_flux_ivar
FROM desi_dr1.agnqso q
WHERE q.z BETWEEN 0.4 AND 1.0
  AND q.mgii_2796_flux > 0 AND q.mgii_2796_flux_ivar > 0
  AND q.hbeta_broad_flux > 0 AND q.hbeta_broad_flux_ivar > 0
  AND q.mgii_2796_flux * SQRT(q.mgii_2796_flux_ivar) > 3
  AND q.hbeta_broad_flux * SQRT(q.hbeta_broad_flux_ivar) > 2
ORDER BY q.random_id"""

d2 = query_desi(sql2)
if d2 is not None:
    z_edges_mgii_hb = np.array([0.4, 0.55, 0.7, 0.85, 1.0])
    snr2_1 = d2[:,3]*np.sqrt(np.maximum(d2[:,4],1e-30))
    snr2_2 = d2[:,5]*np.sqrt(np.maximum(d2[:,6],1e-30))
    lum2 = np.log10(d2[:,3]+d2[:,5]+1e-20)
    r2 = env_delta_test(d2[:,0], d2[:,1], d2[:,2], d2[:,3], d2[:,5],
                        snr2_1, snr2_2, lum2, "MgII/Hbeta_broad", z_edges_mgii_hb)
else:
    print("  MgII/Hβ: query failed or too few objects", flush=True)
    r2 = None

# PAIR 3: CIV(1549) / NeV(3426) — predict LARGEST Δρ
print("\n[3] CIV/NeV (predict largest Δρ)...", flush=True)
sql3 = """SELECT TOP 150000 q.target_ra, q.target_dec, q.z,
    q.civ_1549_flux, q.civ_1549_flux_ivar,
    q.nev_3426_flux, q.nev_3426_flux_ivar
FROM desi_dr1.agnqso q
WHERE q.z BETWEEN 1.0 AND 2.5
  AND q.civ_1549_flux > 0 AND q.civ_1549_flux_ivar > 0
  AND q.nev_3426_flux > 0 AND q.nev_3426_flux_ivar > 0
  AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
  AND q.nev_3426_flux * SQRT(q.nev_3426_flux_ivar) > 2
ORDER BY q.random_id"""

d3 = query_desi(sql3)
if d3 is not None:
    z_edges_nev = np.array([1.0, 1.3, 1.6, 1.9, 2.2, 2.5])
    snr3_1 = d3[:,3]*np.sqrt(np.maximum(d3[:,4],1e-30))
    snr3_2 = d3[:,5]*np.sqrt(np.maximum(d3[:,6],1e-30))
    lum3 = np.log10(d3[:,3]+d3[:,5]+1e-20)
    r3 = env_delta_test(d3[:,0], d3[:,1], d3[:,2], d3[:,3], d3[:,5],
                        snr3_1, snr3_2, lum3, "CIV/NeV", z_edges_nev)
else:
    print("  CIV/NeV: query failed or too few objects", flush=True)
    r3 = None

# ============================================================
print("\n" + "=" * 70, flush=True)
print("PREDICTION SCORECARD", flush=True)
print("=" * 70, flush=True)

pairs = [('CIV/MgII (baseline)', r0),
         ('HeII/CIV (predict small)', r1),
         ('MgII/Hβ (predict large)', r2),
         ('CIV/NeV (predict largest)', r3)]

for name, r in pairs:
    if r:
        print(f"  {name}: Δρ = {r['mean_delta']:+.4f}  ({r['n_positive']}/{r['n_bins']} positive)", flush=True)
    else:
        print(f"  {name}: NO DATA", flush=True)

# Check ordering prediction
deltas = {}
for name, r in pairs:
    if r:
        deltas[name.split('(')[0].strip()] = abs(r['mean_delta'])

print(f"\nOrdering prediction: |HeII/CIV| < |CIV/MgII| < |MgII/Hβ| < |CIV/NeV|", flush=True)
if len(deltas) >= 3:
    sorted_pairs = sorted(deltas.items(), key=lambda x: x[1])
    print(f"Actual ordering:    {' < '.join(f'{k} ({v:.4f})' for k,v in sorted_pairs)}", flush=True)

output = {'predictions': [{'pair': name, 'result': r} for name, r in pairs]}
with open(f'{RESULTS_DIR}/prediction_results.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {RESULTS_DIR}/prediction_results.json", flush=True)
