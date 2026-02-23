#!/usr/bin/env python3
"""
FULL-SAMPLE COPULA + ASSOCIATION TEST
======================================
Designed for PC (64GB RAM). 500K quasars × 500K fg galaxies.
Tests:
1. Spearman/Kendall Δρ with bootstrap significance (dual-SNR matched)
2. Quintile gradient across 5 density bins
3. Lower/upper tail dependence at q=0.80,0.85,0.90,0.95
4. Energy distance (copula) for completeness
5. All 3 z-bins with full statistical power
"""

import numpy as np
import json, os, time
import requests
from scipy.stats import spearmanr, kendalltau, rankdata
from scipy.spatial.distance import cdist

RESULTS_DIR = "results_copula_full"
os.makedirs(RESULTS_DIR, exist_ok=True)
TAP_URL = "https://datalab.noirlab.edu/tap/sync"

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def tap_query(query, timeout=600):
    params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    resp = requests.post(TAP_URL, data=params, timeout=timeout)
    resp.raise_for_status()
    text = resp.text.strip()
    if text.startswith('<?xml'):
        log(f"TAP ERROR: {text[:500]}")
        raise RuntimeError("TAP error")
    lines = text.split('\n')
    return lines[0].split(','), lines[1:]

def parse_csv(header, rows):
    data = {col: [] for col in header}
    for line in rows:
        vals = line.split(',')
        if len(vals) != len(header):
            continue
        try:
            for i, col in enumerate(header):
                data[col].append(float(vals[i]))
        except ValueError:
            continue
    return {col: np.array(data[col]) for col in data}

def energy_distance(X, Y, max_n=5000):
    if len(X) > max_n:
        X = X[np.random.choice(len(X), max_n, replace=False)]
    if len(Y) > max_n:
        Y = Y[np.random.choice(len(Y), max_n, replace=False)]
    xy = cdist(X, Y, 'euclidean').mean()
    xx = cdist(X, X, 'euclidean').mean()
    yy = cdist(Y, Y, 'euclidean').mean()
    return 2 * xy - xx - yy

def tail_dep(u, v, q):
    mask = v > q
    return np.mean(u[mask] > q) if mask.sum() >= 10 else np.nan

def lower_tail_dep(u, v, q):
    t = 1 - q
    mask = v < t
    return np.mean(u[mask] < t) if mask.sum() >= 10 else np.nan

# ── Data acquisition ──────────────────────────────────────────────────

def get_quasars():
    cache = f"{RESULTS_DIR}/quasar_cache_full.npz"
    if os.path.exists(cache):
        log("Loading cached quasars...")
        d = np.load(cache, allow_pickle=True)
        return {k: d[k] for k in d.files}
    
    log("Querying DESI quasars (full sample, z=1.3-2.6)...")
    # Query in batches to avoid timeout
    all_data = None
    for z_lo, z_hi in [(1.3, 1.8), (1.8, 2.2), (2.2, 2.6)]:
        query = f"""
        SELECT TOP 250000
            z.targetid, z.mean_fiber_ra as ra, z.mean_fiber_dec as dec, z.z as redshift,
            a.civ_1549_flux, a.civ_1549_flux_ivar, a.civ_1549_sigma,
            a.mgii_2796_flux, a.mgii_2796_flux_ivar, a.mgii_2796_sigma,
            a.mgii_2803_flux, a.mgii_2803_flux_ivar
        FROM desi_dr1.zpix AS z
        JOIN desi_dr1.agnqso AS a ON z.targetid = a.targetid
        WHERE z.z BETWEEN {z_lo} AND {z_hi}
          AND z.zcat_primary = 'True'
          AND z.zwarn = 0
          AND a.civ_1549_flux_ivar > 0 AND a.mgii_2796_flux_ivar > 0
          AND a.civ_1549_flux > 0 AND a.mgii_2796_flux > 0
          AND a.civ_1549_flux * SQRT(a.civ_1549_flux_ivar) > 5
          AND a.mgii_2796_flux * SQRT(a.mgii_2796_flux_ivar) > 5
          AND a.civ_1549_sigma > 0 AND a.mgii_2796_sigma > 0
        """
        log(f"  Querying z={z_lo}-{z_hi}...")
        header, rows = tap_query(query)
        batch = parse_csv(header, rows)
        log(f"    Got {len(batch['ra'])} rows")
        
        if all_data is None:
            all_data = batch
        else:
            for col in all_data:
                all_data[col] = np.concatenate([all_data[col], batch[col]])
    
    log(f"  Total: {len(all_data['ra'])} quasars")
    np.savez(cache, **all_data)
    return all_data

def get_fg_galaxies():
    cache = f"{RESULTS_DIR}/fg_galaxy_cache_full.npz"
    if os.path.exists(cache):
        log("Loading cached fg galaxies...")
        d = np.load(cache, allow_pickle=True)
        return {k: d[k] for k in d.files}
    
    log("Querying DESI fg galaxies (z<0.8)...")
    # Query in batches
    all_data = None
    for z_lo, z_hi in [(0.01, 0.3), (0.3, 0.55), (0.55, 0.8)]:
        query = f"""
        SELECT TOP 500000
            z.mean_fiber_ra as ra, z.mean_fiber_dec as dec, z.z as redshift
        FROM desi_dr1.zpix AS z
        WHERE z.z BETWEEN {z_lo} AND {z_hi}
          AND z.zcat_primary = 'True'
          AND z.zwarn = 0
          AND z.spectype = 'GALAXY'
        """
        log(f"  Querying fg z={z_lo}-{z_hi}...")
        header, rows = tap_query(query)
        batch = parse_csv(header, rows)
        log(f"    Got {len(batch['ra'])} rows")
        
        if all_data is None:
            all_data = batch
        else:
            for col in all_data:
                all_data[col] = np.concatenate([all_data[col], batch[col]])
    
    log(f"  Total: {len(all_data['ra'])} fg galaxies")
    np.savez(cache, **all_data)
    return all_data

def compute_fg_density(q_ra, q_dec, fg_ra, fg_dec, cone_deg=0.5):
    cache = f"{RESULTS_DIR}/fg_density_full.npy"
    if os.path.exists(cache):
        log("Loading cached fg density...")
        return np.load(cache)
    
    log(f"Computing fg density ({len(q_ra)} × {len(fg_ra)}, cone={cone_deg}°)...")
    cos_limit = np.float32(np.cos(np.radians(cone_deg)))
    
    fg_t = np.radians(90.0 - fg_dec).astype(np.float32)
    fg_p = np.radians(fg_ra).astype(np.float32)
    fg_xyz = np.stack([np.sin(fg_t)*np.cos(fg_p), np.sin(fg_t)*np.sin(fg_p), np.cos(fg_t)], axis=1)
    
    q_t = np.radians(90.0 - q_dec).astype(np.float32)
    q_p = np.radians(q_ra).astype(np.float32)
    q_x = np.sin(q_t)*np.cos(q_p)
    q_y = np.sin(q_t)*np.sin(q_p)
    q_z = np.cos(q_t)
    
    n_q = len(q_ra)
    counts = np.zeros(n_q, dtype=np.int32)
    
    # On 64GB, can do bigger chunks: 5000 × N_fg
    chunk = 5000
    n_chunks = (n_q + chunk - 1) // chunk
    
    for i in range(n_chunks):
        if i % 20 == 0:
            log(f"  Chunk {i+1}/{n_chunks}...")
        i0 = i * chunk
        i1 = min(i0 + chunk, n_q)
        q_xyz = np.stack([q_x[i0:i1], q_y[i0:i1], q_z[i0:i1]], axis=1)
        dots = q_xyz @ fg_xyz.T
        counts[i0:i1] = np.sum(dots >= cos_limit, axis=1)
    
    np.save(cache, counts)
    log(f"  Done. Range: {counts.min()}-{counts.max()}, median={np.median(counts):.0f}")
    return counts

# ── Analysis ──────────────────────────────────────────────────────────

def analyze_bin(civ_f, mgii_f, civ_snr, mgii_snr, density, z_label):
    """Full analysis for one z-bin."""
    n = len(civ_f)
    q_edges = np.percentile(density, [0, 20, 40, 60, 80, 100])
    
    sparse = density <= q_edges[1]
    dense = density >= q_edges[4]
    
    # Dual-SNR matching
    s_cs, d_cs = civ_snr[sparse], civ_snr[dense]
    s_ms, d_ms = mgii_snr[sparse], mgii_snr[dense]
    cl = max(np.percentile(s_cs,10), np.percentile(d_cs,10))
    ch = min(np.percentile(s_cs,90), np.percentile(d_cs,90))
    ml = max(np.percentile(s_ms,10), np.percentile(d_ms,10))
    mh = min(np.percentile(s_ms,90), np.percentile(d_ms,90))
    
    snr_ok = (civ_snr>=cl) & (civ_snr<=ch) & (mgii_snr>=ml) & (mgii_snr<=mh)
    sm = sparse & snr_ok
    dm = dense & snr_ok
    ns, nd = sm.sum(), dm.sum()
    
    log(f"\n  {z_label}: N={n}, matched sparse={ns}, dense={nd}")
    
    if ns < 100 or nd < 100:
        log(f"  Too few, skipping")
        return None
    
    results = {'n_total': int(n), 'n_sparse': int(ns), 'n_dense': int(nd)}
    
    # ── Association ──
    rho_s, _ = spearmanr(civ_f[sm], mgii_f[sm])
    rho_d, _ = spearmanr(civ_f[dm], mgii_f[dm])
    tau_s, _ = kendalltau(civ_f[sm], mgii_f[sm])
    tau_d, _ = kendalltau(civ_f[dm], mgii_f[dm])
    drho = rho_d - rho_s
    
    log(f"  Spearman: sparse={rho_s:.4f}, dense={rho_d:.4f}, Δ={drho:+.4f}")
    log(f"  Kendall:  sparse={tau_s:.4f}, dense={tau_d:.4f}, Δ={tau_d-tau_s:+.4f}")
    
    # Bootstrap Δρ
    all_c = np.concatenate([civ_f[sm], civ_f[dm]])
    all_m = np.concatenate([mgii_f[sm], mgii_f[dm]])
    boot = []
    for _ in range(1000):
        p = np.random.permutation(ns+nd)
        r1, _ = spearmanr(all_c[p[:ns]], all_m[p[:ns]])
        r2, _ = spearmanr(all_c[p[ns:]], all_m[p[ns:]])
        boot.append(r2 - r1)
    bp = np.mean(np.array(boot) >= drho)
    log(f"  Bootstrap p(Δρ): {bp:.4f} {'✅' if bp < 0.05 else '❌'}")
    
    results['association'] = {
        'rho_sparse': float(rho_s), 'rho_dense': float(rho_d),
        'tau_sparse': float(tau_s), 'tau_dense': float(tau_d),
        'delta_rho': float(drho), 'delta_tau': float(tau_d-tau_s),
        'bootstrap_p': float(bp),
    }
    
    # ── Tail dependence ──
    s_u = rankdata(civ_f[sm])/(ns+1)
    s_v = rankdata(mgii_f[sm])/(ns+1)
    d_u = rankdata(civ_f[dm])/(nd+1)
    d_v = rankdata(mgii_f[dm])/(nd+1)
    
    td = {}
    for q in [0.80, 0.85, 0.90, 0.95]:
        lu_s, lu_d = tail_dep(s_u,s_v,q), tail_dep(d_u,d_v,q)
        ll_s, ll_d = lower_tail_dep(s_u,s_v,q), lower_tail_dep(d_u,d_v,q)
        log(f"  q={q}: upper Δ={lu_d-lu_s:+.3f}, lower Δ={ll_d-ll_s:+.3f}")
        td[f"q={q}"] = {
            'upper_sparse': float(lu_s), 'upper_dense': float(lu_d),
            'lower_sparse': float(ll_s), 'lower_dense': float(ll_d),
            'delta_upper': float(lu_d-lu_s), 'delta_lower': float(ll_d-ll_s),
        }
    results['tail_dependence'] = td
    
    # ── Quintile gradient ──
    q_rhos, q_taus, q_dens = [], [], []
    for qi in range(5):
        if qi < 4:
            qm = snr_ok & (density >= q_edges[qi]) & (density < q_edges[qi+1])
        else:
            qm = snr_ok & (density >= q_edges[qi])
        nq = qm.sum()
        if nq < 50:
            continue
        r, _ = spearmanr(civ_f[qm], mgii_f[qm])
        t, _ = kendalltau(civ_f[qm], mgii_f[qm])
        md = np.mean(density[qm])
        q_rhos.append(r); q_taus.append(t); q_dens.append(md)
        log(f"  Q{qi+1} (N={nq}, d={md:.1f}): ρ={r:.4f}")
    
    if len(q_rhos) >= 4:
        rg, pg = spearmanr(q_dens, q_rhos)
        log(f"  Gradient: r={rg:+.3f} (p={pg:.3f}) {'✅' if rg > 0.8 else '⚠️'}")
        results['gradient'] = {
            'densities': [float(x) for x in q_dens],
            'spearman': [float(x) for x in q_rhos],
            'kendall': [float(x) for x in q_taus],
            'r': float(rg), 'p': float(pg),
        }
    
    # ── Copula energy distance ──
    sc = np.column_stack([s_u, s_v])
    dc = np.column_stack([d_u, d_v])
    ed = energy_distance(sc, dc, max_n=5000)
    
    combined = np.vstack([sc, dc])
    boot_ed = []
    for _ in range(200):
        p = np.random.permutation(len(combined))
        boot_ed.append(energy_distance(combined[p[:ns]], combined[p[ns:ns+nd]], max_n=3000))
    ep = np.mean(np.array(boot_ed) >= ed)
    log(f"  Copula ED: {ed:.6f}, p={ep:.4f}")
    results['copula'] = {'energy_distance': float(ed), 'bootstrap_p': float(ep)}
    
    return results

def main():
    np.random.seed(42)
    
    quasars = get_quasars()
    fg = get_fg_galaxies()
    density = compute_fg_density(quasars['ra'], quasars['dec'], fg['ra'], fg['dec'])
    
    civ_f = quasars['civ_1549_flux']
    mgii_f = quasars['mgii_2796_flux']
    civ_snr = civ_f * np.sqrt(quasars['civ_1549_flux_ivar'])
    mgii_snr = mgii_f * np.sqrt(quasars['mgii_2796_flux_ivar'])
    z = quasars['redshift']
    
    all_results = {}
    for z_lo, z_hi in [(1.3, 1.7), (1.7, 2.1), (2.1, 2.6)]:
        zl = f"z={z_lo}-{z_hi}"
        mask = (z >= z_lo) & (z < z_hi) & (density > 0)
        r = analyze_bin(civ_f[mask], mgii_f[mask], civ_snr[mask], mgii_snr[mask],
                       density[mask], zl)
        if r:
            all_results[zl] = r
    
    # ── Summary ──
    log(f"\n{'='*70}")
    log("FINAL SUMMARY")
    log(f"{'='*70}")
    for zl, r in all_results.items():
        a = r['association']
        log(f"\n{zl} (sparse={r['n_sparse']}, dense={r['n_dense']}):")
        log(f"  Δρ={a['delta_rho']:+.4f}, boot p={a['bootstrap_p']:.4f}")
        if 'gradient' in r:
            log(f"  Gradient: r={r['gradient']['r']:+.3f}")
        td = r['tail_dependence']
        log(f"  Lower tail Δλ(q=0.95): {td['q=0.95']['delta_lower']:+.3f}")
        log(f"  Copula ED p={r['copula']['bootstrap_p']:.3f}")
    
    with open(f"{RESULTS_DIR}/copula_full_results.json", 'w') as f:
        json.dump(all_results, f, indent=2)
    log(f"\nSaved to {RESULTS_DIR}/copula_full_results.json")

if __name__ == '__main__':
    main()
