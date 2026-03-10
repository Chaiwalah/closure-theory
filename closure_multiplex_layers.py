#!/usr/bin/env python3
"""
MULTIPLEX LAYER ANALYSIS — GPT's Test A

Treat the observable space as a MULTILAYER network:
  Layer 1 (Optical/Galaxy): Hα, Hβ, stellar continuum, host properties
  Layer 2 (UV/BLR): CIV, CIII, HeII, LyA
  Layer 3 (Bridge): MgII (straddles both regimes)
  Layer 4 (State): LOGMBH, LOGLEDD, LOGL3000 (global properties)

Track per z-bin:
  - Giant component per layer
  - Interlayer correlations (bridge strength)
  - Mutually connected component (nodes connected in BOTH layers)
  - Layer-specific k-core
  - Which layer is "load-bearing" at each z?
"""

import numpy as np
from astropy.io import fits
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_multiplex', exist_ok=True)

print("=" * 70)
print("MULTIPLEX LAYER ANALYSIS — Two competing lattices")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']

# Define layers
layers = {
    'optical': {
        'Halpha_EW': ('HALPHA_BR', 2),
        'Halpha_FWHM': ('HALPHA_BR', 4),
        'Halpha_logL': ('HALPHA_BR', 3),
        'Hbeta_EW': ('HBETA_BR', 2),
        'Hbeta_FWHM': ('HBETA_BR', 4),
        'Hbeta_logL': ('HBETA_BR', 3),
        'LOGL5100': ('LOGL5100', None),
    },
    'uv_blr': {
        'CIV_EW': ('CIV', 2),
        'CIV_FWHM': ('CIV', 4),
        'CIV_logL': ('CIV', 3),
        'CIII_EW': ('CIII_BR', 2),
        'CIII_FWHM': ('CIII_BR', 4),
        'CIII_logL': ('CIII_BR', 3),
        'HeII_EW': ('HEII1640_BR', 2),
        'HeII_FWHM': ('HEII1640_BR', 4),
        'LyA_EW': ('LYA', 2),
        'LyA_FWHM': ('LYA', 4),
        'LyA_logL': ('LYA', 3),
        'LOGL1350': ('LOGL1350', None),
    },
    'bridge': {
        'MgII_EW': ('MGII_BR', 2),
        'MgII_FWHM': ('MGII_BR', 4),
        'MgII_logL': ('MGII_BR', 3),
        'LOGL3000': ('LOGL3000', None),
    },
    'state': {
        'LOGMBH': ('LOGMBH', None),
        'LOGLEDD': ('LOGLEDD_RATIO', None),
    },
}

# Extract all observables
all_obs = {}
for layer_name, obs_dict in layers.items():
    for obs_name, (col, idx) in obs_dict.items():
        if col in d.dtype.names:
            if idx is not None:
                arr = d[col]
                if arr.ndim == 2 and arr.shape[1] > idx:
                    all_obs[obs_name] = arr[:, idx]
            else:
                all_obs[obs_name] = d[col]

print(f"Total observables extracted: {len(all_obs)}")
for layer_name, obs_dict in layers.items():
    obs_in_layer = [k for k in obs_dict if k in all_obs]
    print(f"  {layer_name}: {len(obs_in_layer)} — {obs_in_layer}")

# ============================================================
# Helper functions
# ============================================================
def compute_layer_graph(obs_subset, data_mask, threshold=0.3):
    """Build adjacency for a subset of observables."""
    names = list(obs_subset.keys())
    n = len(names)
    adj = np.zeros((n, n), dtype=int)
    corr = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            vi = obs_subset[names[i]][data_mask]
            vj = obs_subset[names[j]][data_mask]
            valid = np.isfinite(vi) & np.isfinite(vj) & (vi != 0) & (vj != 0)
            if valid.sum() > 30:
                r, p = stats.spearmanr(vi[valid], vj[valid])
                corr[i, j] = r
                corr[j, i] = r
                if abs(r) > threshold:
                    adj[i, j] = 1
                    adj[j, i] = 1
    
    return adj, corr, names

def giant_component_frac(adj):
    """Return fraction of nodes in giant component."""
    n = adj.shape[0]
    if n == 0:
        return 0.0
    visited = [False] * n
    max_comp = 0
    for start in range(n):
        if visited[start]:
            continue
        comp = 0
        queue = [start]
        while queue:
            node = queue.pop(0)
            if visited[node]:
                continue
            visited[node] = True
            comp += 1
            for nb in range(n):
                if adj[node, nb] and not visited[nb]:
                    queue.append(nb)
        max_comp = max(max_comp, comp)
    return max_comp / n

def mean_abs_correlation(corr):
    """Mean |correlation| of non-zero upper triangle."""
    n = corr.shape[0]
    vals = []
    for i in range(n):
        for j in range(i+1, n):
            if corr[i, j] != 0:
                vals.append(abs(corr[i, j]))
    return np.mean(vals) if vals else 0.0

def edge_density(adj):
    """Fraction of possible edges that exist."""
    n = adj.shape[0]
    max_e = n * (n - 1) // 2
    if max_e == 0:
        return 0.0
    return (adj.sum() // 2) / max_e

# ============================================================
# Interlayer coupling: correlations BETWEEN layers
# ============================================================
def interlayer_coupling(layer1_obs, layer2_obs, data_mask, threshold=0.3):
    """Compute mean |correlation| between observables in different layers."""
    names1 = list(layer1_obs.keys())
    names2 = list(layer2_obs.keys())
    
    corrs = []
    edges = 0
    total = 0
    
    for n1 in names1:
        for n2 in names2:
            v1 = layer1_obs[n1][data_mask]
            v2 = layer2_obs[n2][data_mask]
            valid = np.isfinite(v1) & np.isfinite(v2) & (v1 != 0) & (v2 != 0)
            if valid.sum() > 30:
                r, p = stats.spearmanr(v1[valid], v2[valid])
                corrs.append(abs(r))
                total += 1
                if abs(r) > threshold:
                    edges += 1
    
    return np.mean(corrs) if corrs else 0.0, edges, total

# ============================================================
# Run across z-bins
# ============================================================
z_bins = [(0.3, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 0.95), 
          (0.95, 1.05), (1.05, 1.15), (1.15, 1.3), (1.3, 1.6), (1.6, 2.0)]

# Get layer-specific observables
opt_obs = {k: all_obs[k] for k in layers['optical'] if k in all_obs}
uv_obs = {k: all_obs[k] for k in layers['uv_blr'] if k in all_obs}
bridge_obs = {k: all_obs[k] for k in layers['bridge'] if k in all_obs}
state_obs = {k: all_obs[k] for k in layers['state'] if k in all_obs}

results = {}

print(f"\n{'='*70}")
print("LAYER-BY-LAYER GIANT COMPONENT")
print(f"{'='*70}")
print(f"\n  {'z':<12} {'Optical':>10} {'UV/BLR':>10} {'Bridge':>10} {'State':>10} {'O↔UV':>10}")
print(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")

for zlo, zhi in z_bins:
    zmask = (z_all >= zlo) & (z_all < zhi)
    n = zmask.sum()
    if n < 100:
        continue
    
    # Per-layer graphs
    adj_opt, corr_opt, _ = compute_layer_graph(opt_obs, zmask)
    adj_uv, corr_uv, _ = compute_layer_graph(uv_obs, zmask)
    adj_br, corr_br, _ = compute_layer_graph(bridge_obs, zmask)
    adj_st, corr_st, _ = compute_layer_graph(state_obs, zmask)
    
    gf_opt = giant_component_frac(adj_opt)
    gf_uv = giant_component_frac(adj_uv)
    gf_br = giant_component_frac(adj_br)
    gf_st = giant_component_frac(adj_st)
    
    # Interlayer coupling
    coupling_opt_uv, edges_ou, total_ou = interlayer_coupling(opt_obs, uv_obs, zmask)
    coupling_opt_br, _, _ = interlayer_coupling(opt_obs, bridge_obs, zmask)
    coupling_uv_br, _, _ = interlayer_coupling(uv_obs, bridge_obs, zmask)
    
    zkey = f"z=[{zlo},{zhi})"
    results[zkey] = {
        'N': int(n),
        'optical_giant': float(gf_opt),
        'uv_giant': float(gf_uv),
        'bridge_giant': float(gf_br),
        'state_giant': float(gf_st),
        'coupling_opt_uv': float(coupling_opt_uv),
        'coupling_opt_bridge': float(coupling_opt_br),
        'coupling_uv_bridge': float(coupling_uv_br),
        'interlayer_edges': int(edges_ou),
        'interlayer_total': int(total_ou),
    }
    
    # Visual
    opt_bar = "█" * int(gf_opt * 10)
    uv_bar = "█" * int(gf_uv * 10)
    
    print(f"  {zkey:<12} {gf_opt:>9.0%} {gf_uv:>9.0%} {gf_br:>9.0%} {gf_st:>9.0%} {coupling_opt_uv:>9.3f}")

# ============================================================
# The handoff: where does UV overtake Optical?
# ============================================================
print(f"\n{'='*70}")
print("THE HANDOFF — Which layer is load-bearing?")
print(f"{'='*70}")

print(f"\n  {'z':<12} {'Optical':>10} {'UV/BLR':>10} {'Leader':>12} {'Gap':>8}")
print(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*12} {'-'*8}")

for zkey, data in results.items():
    opt = data['optical_giant']
    uv = data['uv_giant']
    leader = "OPTICAL" if opt > uv else "UV/BLR" if uv > opt else "TIE"
    gap = abs(opt - uv)
    print(f"  {zkey:<12} {opt:>9.0%} {uv:>9.0%} {leader:>12} {gap:>7.0%}")

# ============================================================
# Interlayer coupling evolution
# ============================================================
print(f"\n{'='*70}")
print("INTERLAYER COUPLING — Bridge strength between layers")
print(f"{'='*70}")
print(f"\n  {'z':<12} {'Opt↔UV':>10} {'Opt↔Bridge':>12} {'UV↔Bridge':>12} {'Interlayer edges':>18}")
print(f"  {'-'*12} {'-'*10} {'-'*12} {'-'*12} {'-'*18}")

for zkey, data in results.items():
    print(f"  {zkey:<12} {data['coupling_opt_uv']:>10.3f} {data['coupling_opt_bridge']:>12.3f} {data['coupling_uv_bridge']:>12.3f} {data['interlayer_edges']:>8}/{data['interlayer_total']:>3}")

# ============================================================
# MgII as the bridge: coupling to each layer
# ============================================================
print(f"\n{'='*70}")
print("MgII BRIDGE — Coupling to Optical vs UV layer")
print(f"{'='*70}")

mg_ew = all_obs.get('MgII_EW')
if mg_ew is not None:
    print(f"\n  {'z':<12} {'MgII↔Hβ':>10} {'MgII↔CIV':>10} {'MgII↔CIII':>12} {'Allegiance':>12}")
    print(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*12} {'-'*12}")
    
    for zlo, zhi in z_bins:
        zmask = (z_all >= zlo) & (z_all < zhi)
        if zmask.sum() < 100:
            continue
        
        # MgII vs optical representative (Hβ)
        hb = all_obs.get('Hbeta_EW')
        civ = all_obs.get('CIV_EW')
        ciii = all_obs.get('CIII_EW')
        
        r_hb = r_civ = r_ciii = np.nan
        
        if hb is not None:
            v1, v2 = mg_ew[zmask], hb[zmask]
            valid = np.isfinite(v1) & np.isfinite(v2) & (v1 != 0) & (v2 != 0)
            if valid.sum() > 30:
                r_hb, _ = stats.spearmanr(v1[valid], v2[valid])
        
        if civ is not None:
            v1, v2 = mg_ew[zmask], civ[zmask]
            valid = np.isfinite(v1) & np.isfinite(v2) & (v1 != 0) & (v2 != 0)
            if valid.sum() > 30:
                r_civ, _ = stats.spearmanr(v1[valid], v2[valid])
        
        if ciii is not None:
            v1, v2 = mg_ew[zmask], ciii[zmask]
            valid = np.isfinite(v1) & np.isfinite(v2) & (v1 != 0) & (v2 != 0)
            if valid.sum() > 30:
                r_ciii, _ = stats.spearmanr(v1[valid], v2[valid])
        
        # Which layer is MgII more coupled to?
        opt_coupling = abs(r_hb) if not np.isnan(r_hb) else 0
        uv_coupling = max(abs(r_civ) if not np.isnan(r_civ) else 0,
                         abs(r_ciii) if not np.isnan(r_ciii) else 0)
        allegiance = "OPTICAL" if opt_coupling > uv_coupling else "UV/BLR"
        
        zkey = f"z=[{zlo},{zhi})"
        r_hb_str = f"{r_hb:+.3f}" if not np.isnan(r_hb) else "   N/A"
        r_civ_str = f"{r_civ:+.3f}" if not np.isnan(r_civ) else "   N/A"
        r_ciii_str = f"{r_ciii:+.3f}" if not np.isnan(r_ciii) else "   N/A"
        
        print(f"  {zkey:<12} {r_hb_str:>10} {r_civ_str:>10} {r_ciii_str:>12} {allegiance:>12}")

# ============================================================
# Mutually connected component (nodes connected in BOTH layers)
# ============================================================
print(f"\n{'='*70}")
print("MUTUALLY CONNECTED — Nodes connected in BOTH optical AND UV")
print(f"{'='*70}")
print("(A node is 'mutually connected' if it has edges in both layers)")

# For this we need to combine optical + UV into a full graph and track
# which nodes have connections in EACH layer
all_layer_obs = {}
all_layer_obs.update(opt_obs)
all_layer_obs.update(uv_obs)
all_layer_obs.update(bridge_obs)
all_layer_obs.update(state_obs)

all_names = list(all_layer_obs.keys())
opt_indices = set(i for i, n in enumerate(all_names) if n in opt_obs)
uv_indices = set(i for i, n in enumerate(all_names) if n in uv_obs)
bridge_indices = set(i for i, n in enumerate(all_names) if n in bridge_obs)

print(f"\n  {'z':<12} {'Opt-connected':>15} {'UV-connected':>15} {'Both':>8} {'Mutual%':>10}")
print(f"  {'-'*12} {'-'*15} {'-'*15} {'-'*8} {'-'*10}")

for zlo, zhi in z_bins:
    zmask = (z_all >= zlo) & (z_all < zhi)
    if zmask.sum() < 100:
        continue
    
    # Full graph
    n_all = len(all_names)
    adj_full = np.zeros((n_all, n_all), dtype=int)
    
    for i in range(n_all):
        for j in range(i+1, n_all):
            vi = all_layer_obs[all_names[i]][zmask]
            vj = all_layer_obs[all_names[j]][zmask]
            valid = np.isfinite(vi) & np.isfinite(vj) & (vi != 0) & (vj != 0)
            if valid.sum() > 30:
                r, _ = stats.spearmanr(vi[valid], vj[valid])
                if abs(r) > 0.3:
                    adj_full[i, j] = 1
                    adj_full[j, i] = 1
    
    # For each node, check if it has edges to BOTH optical and UV layer nodes
    opt_connected = set()
    uv_connected = set()
    
    for i in range(n_all):
        has_opt_edge = any(adj_full[i, j] for j in opt_indices if j != i)
        has_uv_edge = any(adj_full[i, j] for j in uv_indices if j != i)
        if has_opt_edge:
            opt_connected.add(i)
        if has_uv_edge:
            uv_connected.add(i)
    
    mutual = opt_connected & uv_connected
    mutual_names = [all_names[i] for i in mutual]
    
    zkey = f"z=[{zlo},{zhi})"
    mutual_pct = len(mutual) / n_all if n_all > 0 else 0
    print(f"  {zkey:<12} {len(opt_connected):>15} {len(uv_connected):>15} {len(mutual):>8} {mutual_pct:>9.0%}")
    if mutual_names:
        print(f"  {'':>12} Mutual nodes: {mutual_names[:8]}")

# Save
with open('results_multiplex/multiplex_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\n\nResults saved to results_multiplex/")
