#!/usr/bin/env python3
"""
Hardened Void Test V2 — Three additional tests:
1. Continuous overdensity gradient (quintiles, not just median split)
2. CIV sigma matching (wind/blueshift control — GPT kill #1 survivor test)
3. Kinematic channel independence (sigma correlations × environment)

Uses same 200k DESI agnqso sample as v1.
"""
import numpy as np
from scipy.stats import spearmanr, kendalltau
from scipy.spatial import cKDTree
import requests, json, os, sys, warnings, time
warnings.filterwarnings('ignore')

RESULTS_DIR = 'results_void_hardened'
os.makedirs(RESULTS_DIR, exist_ok=True)

def query_desi(sql, timeout=600):
    print(f"  Querying DESI...", flush=True)
    r = requests.get('https://datalab.noirlab.edu/tap/sync',
                     params={'REQUEST':'doQuery','LANG':'ADQL','QUERY':sql,'FORMAT':'csv'},
                     timeout=timeout)
    if r.status_code != 200 or 'ERROR' in r.text[:500]:
        print(f"  Query error: {r.text[:300]}", flush=True)
        return None, None
    lines = r.text.strip().split('\n')
    header = lines[0].split(',')
    print(f"  Columns: {header}", flush=True)
    print(f"  Rows: {len(lines)-1}", flush=True)
    if len(lines) < 50:
        return None, header
    data = []
    for line in lines[1:]:
        parts = line.split(',')
        try:
            data.append([float(p) for p in parts])
        except:
            continue
    return np.array(data) if len(data) > 100 else None, header

def compute_3d_density(ra, dec, z, n_neighbors=10):
    from scipy.integrate import cumulative_trapezoid
    z_max = max(z) * 1.1
    z_fine = np.linspace(0, z_max, 5000)
    Hz = np.sqrt(0.3*(1+z_fine)**3 + 0.7)
    chi_fine = cumulative_trapezoid(1.0/Hz, z_fine, initial=0) * 2997.9
    chi = np.interp(z, z_fine, chi_fine)
    
    ra_rad, dec_rad = np.radians(ra), np.radians(dec)
    x = chi * np.cos(dec_rad) * np.cos(ra_rad)
    y = chi * np.cos(dec_rad) * np.sin(ra_rad)
    zc = chi * np.sin(dec_rad)
    coords = np.column_stack([x, y, zc])
    
    print(f"  Building 3D tree for {len(coords)} objects...", flush=True)
    tree = cKDTree(coords)
    dists, _ = tree.query(coords, k=n_neighbors+1)
    density = n_neighbors / (4/3 * np.pi * dists[:, -1]**3 + 1e-10)
    return density

def bootstrap_corr(x, y, n_boot=500, metric='spearman'):
    vals = []
    for _ in range(n_boot):
        idx = np.random.choice(len(x), len(x), replace=True)
        if metric == 'spearman':
            vals.append(spearmanr(x[idx], y[idx])[0])
        else:
            vals.append(kendalltau(x[idx], y[idx])[0])
    return np.percentile(vals, [2.5, 50, 97.5])

# ============================================================
print("=" * 70, flush=True)
print("HARDENED VOID TEST V2 — GRADIENT + WIND CONTROL + KINEMATIC", flush=True)
print("=" * 70, flush=True)

# Query: need positions, fluxes, errors, AND sigmas (line widths)
print("\n[1] Querying DESI agnqso with sigmas...", flush=True)
sql = """
SELECT TOP 200000
    q.target_ra, q.target_dec, q.z as redshift,
    q.civ_1549_flux, q.civ_1549_flux_ivar, q.civ_1549_sigma,
    q.mgii_2796_flux, q.mgii_2796_flux_ivar, q.mgii_2796_sigma,
    q.mgii_2803_flux, q.mgii_2803_flux_ivar, q.mgii_2803_sigma,
    q.hbeta_broad_flux, q.hbeta_broad_flux_ivar, q.hbeta_broad_sigma
FROM desi_dr1.agnqso q
WHERE q.z BETWEEN 1.3 AND 2.6
  AND q.civ_1549_flux > 0
  AND q.civ_1549_flux_ivar > 0
  AND q.mgii_2796_flux > 0
  AND q.mgii_2796_flux_ivar > 0
  AND q.civ_1549_flux * SQRT(q.civ_1549_flux_ivar) > 3
  AND q.mgii_2796_flux * SQRT(q.mgii_2796_flux_ivar) > 3
ORDER BY q.random_id
"""

data, header = query_desi(sql)
if data is None:
    print("FATAL: Query failed", flush=True)
    sys.exit(1)

print(f"  Got {len(data)} quasars", flush=True)

ra = data[:, 0]
dec = data[:, 1]
z = data[:, 2]
civ_flux = data[:, 3]
civ_ivar = data[:, 4]
civ_sigma = data[:, 5]
mgii_flux = data[:, 6]
mgii_ivar = data[:, 7]
mgii_sigma = data[:, 8]
mgii2803_flux = data[:, 9]
mgii2803_ivar = data[:, 10]
mgii2803_sigma = data[:, 11]
hb_flux = data[:, 12]
hb_ivar = data[:, 13]
hb_sigma = data[:, 14]

civ_snr = civ_flux * np.sqrt(np.maximum(civ_ivar, 1e-30))
mgii_snr = mgii_flux * np.sqrt(np.maximum(mgii_ivar, 1e-30))
lum_proxy = np.log10(civ_flux + mgii_flux + 1e-20)

print(f"  z range: {z.min():.2f} - {z.max():.2f}", flush=True)

# Compute 3D density
print("\n[2] Computing 3D density...", flush=True)
density = compute_3d_density(ra, dec, z)
log_density = np.log10(density + 1e-30)

z_edges = np.array([1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.6])

# ============================================================
# TEST A: Continuous overdensity gradient (quintiles)
# ============================================================
print("\n" + "=" * 70, flush=True)
print("TEST A: CONTINUOUS OVERDENSITY GRADIENT (QUINTILES)", flush=True)
print("=" * 70, flush=True)

gradient_results = []

for i in range(len(z_edges)-1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    zmask = (z >= z_lo) & (z < z_hi)
    n_bin = zmask.sum()
    if n_bin < 500:
        continue
    
    # Split into quintiles by density
    dens_bin = density[zmask]
    pcts = np.percentile(dens_bin, [20, 40, 60, 80])
    quintile_labels = np.digitize(dens_bin, pcts)  # 0-4
    
    civ_bin = civ_flux[zmask]
    mgii_bin = mgii_flux[zmask]
    
    q_results = []
    print(f"\n  z=[{z_lo:.1f},{z_hi:.1f}] N={n_bin}", flush=True)
    for q in range(5):
        qmask = quintile_labels == q
        if qmask.sum() < 50:
            continue
        rho = spearmanr(civ_bin[qmask], mgii_bin[qmask])[0]
        ci = bootstrap_corr(civ_bin[qmask], mgii_bin[qmask], n_boot=500)
        q_results.append({
            'quintile': q+1,
            'label': ['Q1(sparsest)','Q2','Q3','Q4','Q5(densest)'][q],
            'n': int(qmask.sum()),
            'spearman': round(rho, 4),
            'ci_lo': round(ci[0], 4),
            'ci_hi': round(ci[2], 4),
        })
        print(f"    Q{q+1}: ρ={rho:.3f} [{ci[0]:.3f}, {ci[2]:.3f}] N={qmask.sum()}", flush=True)
    
    # Linear trend across quintiles
    if len(q_results) >= 4:
        rhos = [r['spearman'] for r in q_results]
        qs = list(range(1, len(rhos)+1))
        trend_r, trend_p = spearmanr(qs, rhos)
        print(f"    Gradient: r={trend_r:.3f}, p={trend_p:.4f}", flush=True)
    else:
        trend_r, trend_p = np.nan, np.nan
    
    gradient_results.append({
        'z_range': f"{z_lo:.1f}-{z_hi:.1f}",
        'n_total': int(n_bin),
        'quintiles': q_results,
        'gradient_r': round(trend_r, 4) if not np.isnan(trend_r) else None,
        'gradient_p': round(trend_p, 4) if not np.isnan(trend_p) else None,
    })

# ============================================================
# TEST B: CIV sigma matching (wind control)
# ============================================================
print("\n" + "=" * 70, flush=True)
print("TEST B: CIV SIGMA (WIND/BLUESHIFT) MATCHING", flush=True)
print("=" * 70, flush=True)
print("  If CIV wind properties differ by environment, matching on", flush=True)
print("  CIV sigma should kill the signal if it's population-driven.", flush=True)

wind_results = []

for i in range(len(z_edges)-1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    zmask = (z >= z_lo) & (z < z_hi)
    # Also require valid sigmas
    valid = zmask & (civ_sigma > 0) & (mgii_sigma > 0)
    n_valid = valid.sum()
    if n_valid < 500:
        print(f"  z=[{z_lo:.1f},{z_hi:.1f}]: {n_valid} with valid sigma — SKIP", flush=True)
        continue
    
    # Median density split
    med = np.median(density[valid])
    dense_m = valid & (density >= med)
    sparse_m = valid & (density < med)
    
    # Match on: CIV_SNR + MgII_SNR + luminosity + CIV_sigma (the wind proxy)
    d_idx = np.where(dense_m)[0]
    s_idx = np.where(sparse_m)[0]
    
    # 4D matching via binning
    n_bins = 5
    snr_combined = np.sqrt(civ_snr**2 + mgii_snr**2)
    
    vars_to_match = [
        (snr_combined, 'SNR'),
        (lum_proxy, 'Lum'),
        (civ_sigma, 'CIV_sigma'),
    ]
    
    # Simple: bin on each variable, intersect
    matched_d, matched_s = [], []
    snr_bins = np.percentile(snr_combined[valid], np.linspace(0, 100, n_bins+1))
    lum_bins = np.percentile(lum_proxy[valid], np.linspace(0, 100, n_bins+1))
    sig_bins = np.percentile(civ_sigma[valid], np.linspace(0, 100, n_bins+1))
    
    for si in range(n_bins):
        for li in range(n_bins):
            for wi in range(n_bins):
                d_in = d_idx[
                    (snr_combined[d_idx] >= snr_bins[si]) & (snr_combined[d_idx] < snr_bins[si+1]) &
                    (lum_proxy[d_idx] >= lum_bins[li]) & (lum_proxy[d_idx] < lum_bins[li+1]) &
                    (civ_sigma[d_idx] >= sig_bins[wi]) & (civ_sigma[d_idx] < sig_bins[wi+1])
                ]
                s_in = s_idx[
                    (snr_combined[s_idx] >= snr_bins[si]) & (snr_combined[s_idx] < snr_bins[si+1]) &
                    (lum_proxy[s_idx] >= lum_bins[li]) & (lum_proxy[s_idx] < lum_bins[li+1]) &
                    (civ_sigma[s_idx] >= sig_bins[wi]) & (civ_sigma[s_idx] < sig_bins[wi+1])
                ]
                nm = min(len(d_in), len(s_in))
                if nm > 0:
                    matched_d.extend(np.random.choice(d_in, nm, replace=False))
                    matched_s.extend(np.random.choice(s_in, nm, replace=False))
    
    matched_d = np.array(matched_d)
    matched_s = np.array(matched_s)
    
    if len(matched_d) < 100:
        print(f"  z=[{z_lo:.1f},{z_hi:.1f}]: matched {len(matched_d)} — too few", flush=True)
        continue
    
    rho_d = spearmanr(civ_flux[matched_d], mgii_flux[matched_d])[0]
    rho_s = spearmanr(civ_flux[matched_s], mgii_flux[matched_s])[0]
    delta = rho_d - rho_s
    
    # Bootstrap CI
    deltas_boot = []
    max_n = min(3000, len(matched_d))
    d_sub = np.random.choice(matched_d, max_n, replace=False) if len(matched_d) > max_n else matched_d
    s_sub = np.random.choice(matched_s, max_n, replace=False) if len(matched_s) > max_n else matched_s
    for _ in range(500):
        bi_d = np.random.choice(d_sub, len(d_sub), replace=True)
        bi_s = np.random.choice(s_sub, len(s_sub), replace=True)
        deltas_boot.append(spearmanr(civ_flux[bi_d], mgii_flux[bi_d])[0] - 
                          spearmanr(civ_flux[bi_s], mgii_flux[bi_s])[0])
    ci_lo, ci_hi = np.percentile(deltas_boot, [2.5, 97.5])
    
    print(f"  z=[{z_lo:.1f},{z_hi:.1f}] matched={len(matched_d)} "
          f"dense={rho_d:.3f} sparse={rho_s:.3f} Δ={delta:+.3f} CI=[{ci_lo:+.3f},{ci_hi:+.3f}]", flush=True)
    
    wind_results.append({
        'z_range': f"{z_lo:.1f}-{z_hi:.1f}",
        'n_matched': int(len(matched_d)),
        'dense_spearman': round(rho_d, 4),
        'sparse_spearman': round(rho_s, 4),
        'delta': round(delta, 4),
        'ci_lo': round(ci_lo, 4),
        'ci_hi': round(ci_hi, 4),
        'matching': 'SNR + Lum + CIV_sigma (5x5x5 grid)',
    })

# ============================================================
# TEST C: Kinematic channel × environment
# ============================================================
print("\n" + "=" * 70, flush=True)
print("TEST C: KINEMATIC CHANNEL (SIGMA CORRELATIONS × ENVIRONMENT)", flush=True)
print("=" * 70, flush=True)
print("  If spectral degrades but kinematic doesn't → different bandwidth.", flush=True)
print("  CIV_sigma vs MgII_sigma correlation, dense vs sparse.", flush=True)

kinematic_results = []

for i in range(len(z_edges)-1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    valid = (z >= z_lo) & (z < z_hi) & (civ_sigma > 0) & (mgii_sigma > 0)
    n_valid = valid.sum()
    if n_valid < 500:
        print(f"  z=[{z_lo:.1f},{z_hi:.1f}]: {n_valid} valid — SKIP", flush=True)
        continue
    
    med = np.median(density[valid])
    dense_m = valid & (density >= med)
    sparse_m = valid & (density < med)
    
    # SPECTRAL channel: flux correlations
    spec_d = spearmanr(civ_flux[dense_m], mgii_flux[dense_m])[0]
    spec_s = spearmanr(civ_flux[sparse_m], mgii_flux[sparse_m])[0]
    spec_delta = spec_d - spec_s
    
    # KINEMATIC channel: sigma correlations  
    kin_d = spearmanr(civ_sigma[dense_m], mgii_sigma[dense_m])[0]
    kin_s = spearmanr(civ_sigma[sparse_m], mgii_sigma[sparse_m])[0]
    kin_delta = kin_d - kin_s
    
    print(f"  z=[{z_lo:.1f},{z_hi:.1f}] N={n_valid}", flush=True)
    print(f"    SPECTRAL (flux):  dense={spec_d:.3f} sparse={spec_s:.3f} Δ={spec_delta:+.3f}", flush=True)
    print(f"    KINEMATIC (sigma): dense={kin_d:.3f} sparse={kin_s:.3f} Δ={kin_delta:+.3f}", flush=True)
    
    kinematic_results.append({
        'z_range': f"{z_lo:.1f}-{z_hi:.1f}",
        'n': int(n_valid),
        'spectral_dense': round(spec_d, 4),
        'spectral_sparse': round(spec_s, 4),
        'spectral_delta': round(spec_delta, 4),
        'kinematic_dense': round(kin_d, 4),
        'kinematic_sparse': round(kin_s, 4),
        'kinematic_delta': round(kin_delta, 4),
        'channel_divergence': round(spec_delta - kin_delta, 4),
    })

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70, flush=True)
print("SUMMARY", flush=True)
print("=" * 70, flush=True)

# Gradient
print("\nTEST A — Gradient:", flush=True)
n_positive_gradient = sum(1 for g in gradient_results if g['gradient_r'] and g['gradient_r'] > 0)
print(f"  Positive gradient (denser=higher ρ) in {n_positive_gradient}/{len(gradient_results)} bins", flush=True)

# Wind control
print("\nTEST B — Wind control:", flush=True)
n_wind_positive = sum(1 for w in wind_results if w['delta'] > 0)
print(f"  Dense > Sparse after CIV_sigma matching: {n_wind_positive}/{len(wind_results)} bins", flush=True)
wind_mean_delta = np.mean([w['delta'] for w in wind_results]) if wind_results else 0
print(f"  Mean Δ: {wind_mean_delta:+.4f}", flush=True)

# Kinematic
print("\nTEST C — Kinematic independence:", flush=True)
for kr in kinematic_results:
    direction = "DIVERGENT" if (kr['spectral_delta'] > 0 and kr['kinematic_delta'] <= 0) else \
                "CONVERGENT" if (kr['spectral_delta'] > 0 and kr['kinematic_delta'] > 0) else "OTHER"
    print(f"  z={kr['z_range']}: spectral Δ={kr['spectral_delta']:+.3f}, "
          f"kinematic Δ={kr['kinematic_delta']:+.3f} → {direction}", flush=True)

# Save
output = {
    'test': 'Hardened Void V2 — Gradient + Wind + Kinematic',
    'gradient': gradient_results,
    'wind_control': wind_results,
    'kinematic_independence': kinematic_results,
}

with open(f'{RESULTS_DIR}/v2_results.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to {RESULTS_DIR}/v2_results.json", flush=True)
