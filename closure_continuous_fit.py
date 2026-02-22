#!/usr/bin/env python3
"""
Continuous Overdensity Model: ρ(z, δ) = ρ₀(z) + A(z)·f(δ)
============================================================
Replace binary dense/sparse with continuous density field.
Fit a transport-law-like model showing A(z) grows with z.
"""
import numpy as np
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import curve_fit
import requests, json, os, warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = 'results_continuous_fit'
os.makedirs(RESULTS_DIR, exist_ok=True)

def query_desi(sql, timeout=600):
    r = requests.get('https://datalab.noirlab.edu/tap/sync',
                     params={'REQUEST':'doQuery','LANG':'ADQL','QUERY':sql,'FORMAT':'csv'},
                     timeout=timeout)
    lines = r.text.strip().split('\n')
    data = []
    for line in lines[1:]:
        try: data.append([float(p) for p in line.split(',')])
        except: continue
    return np.array(data) if len(data) > 100 else None

def compute_3d_density(ra, dec, z, n_neighbors=10):
    z_fine = np.linspace(0, max(z)*1.1, 5000)
    Hz = np.sqrt(0.3*(1+z_fine)**3 + 0.7)
    chi_fine = cumulative_trapezoid(1.0/Hz, z_fine, initial=0) * 2997.9
    chi = np.interp(z, z_fine, chi_fine)
    ra_r, dec_r = np.radians(ra), np.radians(dec)
    coords = np.column_stack([chi*np.cos(dec_r)*np.cos(ra_r),
                               chi*np.cos(dec_r)*np.sin(ra_r),
                               chi*np.sin(dec_r)])
    tree = cKDTree(coords)
    dists, _ = tree.query(coords, k=n_neighbors+1)
    return n_neighbors / (4/3 * np.pi * dists[:, -1]**3 + 1e-10)

print("=" * 70, flush=True)
print("CONTINUOUS OVERDENSITY MODEL", flush=True)
print("=" * 70, flush=True)

# Reuse our standard query
print("\n[1] Querying...", flush=True)
sql = """SELECT TOP 150000
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

data = query_desi(sql)
print(f"  Got {len(data)} objects", flush=True)

ra, dec, z = data[:,0], data[:,1], data[:,2]
civ_f, mgii_f = data[:,3], data[:,5]

print("\n[2] Computing density...", flush=True)
density = compute_3d_density(ra, dec, z)
log_dens = np.log10(density + 1e-30)

# Normalize log_density to [0,1] within each z-bin for comparison
z_edges = np.arange(1.3, 2.7, 0.1)  # fine z-bins

print("\n[3] Computing ρ(CIV,MgII) in (z, δ) cells...", flush=True)

# For each z-bin, split into 5 density quintiles and compute Spearman
z_coarse = np.array([1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.6])
n_quint = 5

model_data = []  # (z_mid, delta_quantile, rho, n)

for i in range(len(z_coarse)-1):
    z_lo, z_hi = z_coarse[i], z_coarse[i+1]
    zmask = (z >= z_lo) & (z < z_hi)
    if zmask.sum() < 500:
        continue
    
    z_mid = (z_lo + z_hi) / 2
    dens_bin = log_dens[zmask]
    pcts = np.percentile(dens_bin, np.linspace(0, 100, n_quint+1))
    
    for q in range(n_quint):
        qmask_local = (dens_bin >= pcts[q]) & (dens_bin < pcts[q+1] + (1e-10 if q == n_quint-1 else 0))
        global_mask = zmask.copy()
        idx_in_bin = np.where(zmask)[0]
        qidx = idx_in_bin[qmask_local]
        
        if len(qidx) < 50:
            continue
        
        rho = spearmanr(civ_f[qidx], mgii_f[qidx])[0]
        delta_mid = np.median(dens_bin[qmask_local])
        
        model_data.append({
            'z_mid': z_mid,
            'log_density': round(delta_mid, 3),
            'quintile': q + 1,
            'spearman': round(rho, 4),
            'n': len(qidx),
        })

print(f"  Computed {len(model_data)} (z, δ) cells", flush=True)

# ============================================================
# FIT: ρ = ρ₀(z) + A(z) · log_density
# Within each z-bin, fit linear relation ρ vs log_density
# ============================================================
print("\n" + "=" * 70, flush=True)
print("FIT: ρ = ρ₀(z) + A(z) · log(δ)", flush=True)
print("=" * 70, flush=True)

fit_results = []
for i in range(len(z_coarse)-1):
    z_lo, z_hi = z_coarse[i], z_coarse[i+1]
    z_mid = (z_lo + z_hi) / 2
    
    cells = [c for c in model_data if abs(c['z_mid'] - z_mid) < 0.01]
    if len(cells) < 4:
        continue
    
    deltas = [c['log_density'] for c in cells]
    rhos = [c['spearman'] for c in cells]
    
    # Linear fit
    coeffs = np.polyfit(deltas, rhos, 1)
    A_z = coeffs[0]  # slope = sensitivity to density
    rho_0 = coeffs[1]  # intercept = baseline correlation
    
    # Spearman trend
    trend_r, trend_p = spearmanr(deltas, rhos)
    
    print(f"  z={z_mid:.2f}: A(z) = {A_z:+.4f}, ρ₀ = {rho_0:.3f}, "
          f"trend r={trend_r:.3f} (p={trend_p:.3f})", flush=True)
    
    fit_results.append({
        'z_mid': z_mid,
        'A_z': round(A_z, 5),
        'rho_0': round(rho_0, 4),
        'trend_r': round(trend_r, 3),
        'trend_p': round(trend_p, 4),
        'n_cells': len(cells),
    })

# Does A(z) grow with z?
if len(fit_results) >= 3:
    z_vals = [f['z_mid'] for f in fit_results]
    A_vals = [f['A_z'] for f in fit_results]
    Az_trend_r, Az_trend_p = spearmanr(z_vals, A_vals)
    
    print(f"\n  A(z) trend with z: r={Az_trend_r:.3f}, p={Az_trend_p:.4f}", flush=True)
    if Az_trend_r > 0:
        print(f"  ✅ A(z) GROWS with z — environment sensitivity increases with distance", flush=True)
    else:
        print(f"  ❌ A(z) does not grow with z", flush=True)

# SUMMARY
print("\n" + "=" * 70, flush=True)
print("SUMMARY", flush=True)
print("=" * 70, flush=True)

n_positive_A = sum(1 for f in fit_results if f['A_z'] > 0)
print(f"  A(z) > 0 in {n_positive_A}/{len(fit_results)} z-bins (positive = denser → higher ρ)", flush=True)
print(f"  A(z) trend with z: r={Az_trend_r:.3f} (p={Az_trend_p:.4f})", flush=True)

mean_A = np.mean(A_vals)
print(f"  Mean A(z): {mean_A:+.5f}", flush=True)
print(f"  Interpretation: each unit increase in log(density) raises Spearman by {mean_A:.4f}", flush=True)

output = {
    'test': 'Continuous overdensity model ρ(z,δ)',
    'model': 'ρ = ρ₀(z) + A(z)·log(δ)',
    'cells': model_data,
    'fits': fit_results,
    'A_z_trend': {'r': round(Az_trend_r, 3), 'p': round(Az_trend_p, 4)},
}
with open(f'{RESULTS_DIR}/continuous_fit_results.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {RESULTS_DIR}/continuous_fit_results.json", flush=True)
