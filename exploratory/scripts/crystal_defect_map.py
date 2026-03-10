#!/usr/bin/env python3
"""
CRYSTAL DEFECT MAP — Where Does the Membrane Fail?
=====================================================
If the vacuum is a crystal with transport properties, defects should 
appear as sightlines where correlation degradation is anomalously 
LOW (defect = hole in membrane) or HIGH (stressed lattice/grain boundary).

Method:
1. Divide sky into HEALPix-like patches (RA/DEC grid)
2. In each patch, measure EW coupling degradation vs z
3. Compare to global sigmoid → compute residual
4. Map residuals onto sky
5. Look for spatial clustering (defects are localized)

Additional tests:
- Do outlier sightlines show ALL lines preserved (not just one)?
- Do outlier positions correlate with known cosmic structures?
- Is the scatter spatially structured or random?

Data: DR16Q (750K quasars with RA, DEC, z, emission line properties)
"""

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from astropy.io import fits
from collections import defaultdict
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

RESULTS_DIR = 'results_crystal_defect'

def load_data():
    f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
    d = f[1].data
    z = d['Z_DR16Q']
    ra = d['RA']
    dec = d['DEC']
    snr = d['SN_MEDIAN_ALL']
    
    # Extract emission line EWs
    lines = {}
    for name in ['OIII5007', 'HBETA', 'NII6585', 'SII6718', 'OII3728', 'MGII', 'CIV']:
        col = d[name]
        if col.ndim == 2 and col.shape[1] >= 6:
            lines[name] = {
                'ew': col[:, 2],
                'fwhm': col[:, 4],
                'flux': col[:, 0]
            }
    
    f.close()
    return z, ra, dec, snr, lines

def measure_ew_coupling(ew_a, ew_b, mask):
    """Pearson correlation between two EW arrays."""
    if mask.sum() < 30:
        return np.nan, np.nan, mask.sum()
    a = ew_a[mask]
    b = ew_b[mask]
    valid = np.isfinite(a) & np.isfinite(b) & (a != 0) & (b != 0)
    if valid.sum() < 30:
        return np.nan, np.nan, valid.sum()
    r, p = stats.pearsonr(a[valid], b[valid])
    return r, p, valid.sum()

def sigmoid(z, A, z0, w, C):
    return A / (1.0 + np.exp(-(z - z0) / w)) + C

print("=" * 70)
print("CRYSTAL DEFECT MAP — DR16Q 750K Quasars")
print("=" * 70)

z, ra, dec, snr, lines = load_data()
N = len(z)
print(f"Loaded {N} quasars")

# Quality cut
good = (snr > 2) & (z > 0.3) & (z < 2.5)
print(f"After quality cuts (SNR>2, 0.3<z<2.5): {good.sum()}")

# ============================================================
# STEP 1: Establish global sigmoid baseline
# ============================================================
print("\n" + "=" * 70)
print("STEP 1: Global EW Coupling Sigmoid")
print("=" * 70)

# Use OIII-HBETA as the primary diagnostic pair (best coverage at low-mid z)
# And CIV-MGII for high z
pairs = [
    ("OIII-HBETA", "OIII5007", "HBETA", 0.3, 1.0),
    ("CIV-MGII", "CIV", "MGII", 1.0, 2.5),
]

global_curves = {}

for pair_name, line_a, line_b, z_lo, z_hi in pairs:
    if line_a not in lines or line_b not in lines:
        print(f"  {pair_name}: missing data, skipping")
        continue
    
    ew_a = lines[line_a]['ew']
    ew_b = lines[line_b]['ew']
    
    mask_range = good & (z >= z_lo) & (z < z_hi)
    
    n_bins = 8
    z_edges = np.percentile(z[mask_range], np.linspace(0, 100, n_bins + 1))
    
    z_mids, rhos, ns = [], [], []
    
    print(f"\n  {pair_name} (z={z_lo}-{z_hi}):")
    for i in range(n_bins):
        mask = mask_range & (z >= z_edges[i]) & (z < z_edges[i+1])
        r, p, n = measure_ew_coupling(ew_a, ew_b, mask)
        if np.isfinite(r):
            z_mids.append(np.median(z[mask]))
            rhos.append(r)
            ns.append(n)
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"    z={z_mids[-1]:.3f}: ρ={r:+.4f}  N={n:6d}  {sig}")
    
    global_curves[pair_name] = {
        'z': np.array(z_mids),
        'rho': np.array(rhos),
        'n': np.array(ns),
        'line_a': line_a,
        'line_b': line_b,
        'z_range': (z_lo, z_hi)
    }

# ============================================================
# STEP 2: Sky Patch Grid
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Sky Patch Residuals")
print("=" * 70)

# Divide sky into patches (RA bins × DEC bins)
N_RA = 12   # 30° strips
N_DEC = 6   # ~30° strips
ra_edges = np.linspace(0, 360, N_RA + 1)
dec_edges = np.linspace(-10, 85, N_DEC + 1)  # SDSS footprint roughly

print(f"  Grid: {N_RA} RA × {N_DEC} DEC = {N_RA * N_DEC} patches")

patch_results = []
anomalies = []

for pair_name, curve in global_curves.items():
    line_a = curve['line_a']
    line_b = curve['line_b']
    z_lo, z_hi = curve['z_range']
    ew_a = lines[line_a]['ew']
    ew_b = lines[line_b]['ew']
    
    # Global mean coupling for this pair
    mask_all = good & (z >= z_lo) & (z < z_hi)
    r_global, _, n_global = measure_ew_coupling(ew_a, ew_b, mask_all)
    
    print(f"\n  {pair_name}: Global ρ = {r_global:+.4f} (N={n_global})")
    print(f"  {'Patch':>20s}  {'ρ_local':>8s}  {'ρ_global':>8s}  {'Δρ':>8s}  {'N':>6s}  {'σ':>6s}")
    
    for i_ra in range(N_RA):
        for i_dec in range(N_DEC):
            ra_lo, ra_hi = ra_edges[i_ra], ra_edges[i_ra + 1]
            dec_lo, dec_hi = dec_edges[i_dec], dec_edges[i_dec + 1]
            
            mask_patch = (mask_all & 
                         (ra >= ra_lo) & (ra < ra_hi) &
                         (dec >= dec_lo) & (dec < dec_hi))
            
            r_local, p_local, n_local = measure_ew_coupling(ew_a, ew_b, mask_patch)
            
            if np.isfinite(r_local) and n_local >= 50:
                delta_r = r_local - r_global
                # Bootstrap error estimate
                se = 1.0 / np.sqrt(n_local - 3)  # Fisher z approx
                sigma = delta_r / se if se > 0 else 0
                
                patch_results.append({
                    'pair': pair_name,
                    'ra_center': (ra_lo + ra_hi) / 2,
                    'dec_center': (dec_lo + dec_hi) / 2,
                    'ra_range': f"{ra_lo:.0f}-{ra_hi:.0f}",
                    'dec_range': f"{dec_lo:.0f}-{dec_hi:.0f}",
                    'rho_local': float(r_local),
                    'rho_global': float(r_global),
                    'delta_rho': float(delta_r),
                    'sigma': float(sigma),
                    'n': int(n_local)
                })
                
                if abs(sigma) > 2.5:
                    flag = "🔴" if delta_r > 0 else "🔵"
                    label = f"RA[{ra_lo:.0f}-{ra_hi:.0f}] DEC[{dec_lo:.0f}-{dec_hi:.0f}]"
                    print(f"  {label:>20s}  {r_local:+8.4f}  {r_global:+8.4f}  "
                          f"{delta_r:+8.4f}  {n_local:6d}  {sigma:+6.1f}  {flag}")
                    
                    anomalies.append({
                        'patch': label,
                        'pair': pair_name,
                        'delta_rho': float(delta_r),
                        'sigma': float(sigma),
                        'n': int(n_local),
                        'type': 'DEFECT (preserved)' if delta_r > 0 else 'STRESS (excess degradation)',
                        'ra_center': (ra_lo + ra_hi) / 2,
                        'dec_center': (dec_lo + dec_hi) / 2
                    })

# ============================================================
# STEP 3: Spatial Clustering Test
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Spatial Structure in Residuals")
print("=" * 70)

if len(patch_results) > 10:
    deltas = np.array([p['delta_rho'] for p in patch_results])
    ras = np.array([p['ra_center'] for p in patch_results])
    decs = np.array([p['dec_center'] for p in patch_results])
    
    # Is delta_rho spatially correlated?
    # Check: do nearby patches have similar residuals?
    n_patches = len(patch_results)
    near_pairs = []
    far_pairs = []
    
    for i in range(n_patches):
        for j in range(i+1, n_patches):
            if patch_results[i]['pair'] != patch_results[j]['pair']:
                continue
            # Angular separation (rough)
            dra = min(abs(ras[i] - ras[j]), 360 - abs(ras[i] - ras[j]))
            ddec = abs(decs[i] - decs[j])
            sep = np.sqrt(dra**2 + ddec**2)
            
            if sep < 45:  # nearby
                near_pairs.append(deltas[i] * deltas[j])
            elif sep > 90:  # far
                far_pairs.append(deltas[i] * deltas[j])
    
    if near_pairs and far_pairs:
        mean_near = np.mean(near_pairs)
        mean_far = np.mean(far_pairs)
        
        print(f"\n  Spatial correlation test:")
        print(f"    Mean Δρ product (nearby patches, <45°): {mean_near:+.6f}")
        print(f"    Mean Δρ product (far patches, >90°):    {mean_far:+.6f}")
        print(f"    Ratio: {mean_near/mean_far:.2f}" if abs(mean_far) > 1e-8 else "    Ratio: undefined")
        
        if mean_near > mean_far and mean_near > 0:
            print(f"    ✓ SPATIALLY CORRELATED — nearby patches have similar residuals")
            print(f"      → Crystal defects are LOCALIZED, not random")
        else:
            print(f"    ~ No strong spatial correlation detected")
    
    # Overall distribution
    print(f"\n  Residual distribution:")
    print(f"    Mean Δρ: {np.mean(deltas):+.4f}")
    print(f"    Std Δρ:  {np.std(deltas):.4f}")
    print(f"    Skew:    {stats.skew(deltas):.3f}")
    print(f"    Kurt:    {stats.kurtosis(deltas):.3f}")
    
    # Normality test
    _, p_norm = stats.normaltest(deltas)
    print(f"    Normal test p = {p_norm:.4f}")
    if p_norm < 0.05:
        print(f"    ⚠️ NON-NORMAL — excess structure beyond random scatter")

# ============================================================
# STEP 4: Multi-line consistency check for anomalous patches
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Multi-Line Consistency in Anomalous Patches")
print("=" * 70)

if anomalies:
    print(f"\n  {len(anomalies)} anomalous patches detected (|σ| > 2.5)")
    
    # For each anomalous patch, check if OTHER line pairs show the same direction
    for anom in anomalies:
        ra_c = anom['ra_center']
        dec_c = anom['dec_center']
        print(f"\n  {anom['patch']} ({anom['type']}):")
        print(f"    Primary: {anom['pair']} Δρ={anom['delta_rho']:+.4f} ({anom['sigma']:+.1f}σ) N={anom['n']}")
        
        # Check all available line pairs in this patch
        for pair_name, curve in global_curves.items():
            if pair_name == anom['pair']:
                continue
            
            line_a = curve['line_a']
            line_b = curve['line_b']
            z_lo, z_hi = curve['z_range']
            ew_a = lines[line_a]['ew']
            ew_b = lines[line_b]['ew']
            
            mask_all = good & (z >= z_lo) & (z < z_hi)
            r_global, _, _ = measure_ew_coupling(ew_a, ew_b, mask_all)
            
            # Same spatial patch
            mask_patch = (mask_all & 
                         (ra >= ra_c - 15) & (ra < ra_c + 15) &
                         (dec >= dec_c - 15) & (dec < dec_c + 15))
            
            r_local, p_local, n_local = measure_ew_coupling(ew_a, ew_b, mask_patch)
            
            if np.isfinite(r_local) and n_local >= 30:
                delta = r_local - r_global
                consistent = (delta > 0) == (anom['delta_rho'] > 0)
                print(f"    Cross-check: {pair_name} Δρ={delta:+.4f} N={n_local} "
                      f"{'✓ SAME DIRECTION' if consistent else '✗ OPPOSITE'}")

# ============================================================
# STEP 5: Redshift-resolved defect signature
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Anomalous Sightlines — z-Resolved")
print("=" * 70)

# For the strongest anomaly, trace coupling vs z in detail
if anomalies:
    best = max(anomalies, key=lambda x: abs(x['sigma']))
    ra_c = best['ra_center']
    dec_c = best['dec_center']
    pair_name = best['pair']
    curve = global_curves[pair_name]
    
    line_a = curve['line_a']
    line_b = curve['line_b']
    z_lo, z_hi = curve['z_range']
    ew_a = lines[line_a]['ew']
    ew_b = lines[line_b]['ew']
    
    mask_patch = (good & (z >= z_lo) & (z < z_hi) &
                  (ra >= ra_c - 15) & (ra < ra_c + 15) &
                  (dec >= dec_c - 15) & (dec < dec_c + 15))
    
    mask_global = good & (z >= z_lo) & (z < z_hi)
    
    print(f"\n  Strongest anomaly: {best['patch']} ({best['pair']}, {best['sigma']:+.1f}σ)")
    print(f"\n  {'z_bin':>12s}  {'ρ_patch':>8s}  {'ρ_global':>8s}  {'Δρ':>8s}  {'N_patch':>7s}")
    
    n_zbins = 4
    z_edges_5 = np.percentile(z[mask_patch], np.linspace(0, 100, n_zbins + 1))
    z_edges_g = np.percentile(z[mask_global], np.linspace(0, 100, n_zbins + 1))
    
    for i in range(n_zbins):
        m_p = mask_patch & (z >= z_edges_5[i]) & (z < z_edges_5[i+1])
        m_g = mask_global & (z >= z_edges_g[i]) & (z < z_edges_g[i+1])
        
        r_p, _, n_p = measure_ew_coupling(ew_a, ew_b, m_p)
        r_g, _, n_g = measure_ew_coupling(ew_a, ew_b, m_g)
        
        if np.isfinite(r_p) and np.isfinite(r_g):
            delta = r_p - r_g
            print(f"  [{z_edges_5[i]:.2f},{z_edges_5[i+1]:.2f}]  {r_p:+8.4f}  {r_g:+8.4f}  "
                  f"{delta:+8.4f}  {n_p:7d}")

# ============================================================
# STEP 6: Known structure cross-reference
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: Cross-Reference with Known Cosmic Structures")
print("=" * 70)

known_structures = [
    {"name": "CMB Cold Spot (Eridanus)", "ra": 50, "dec": -20, "type": "void"},
    {"name": "Sloan Great Wall", "ra": 195, "dec": 6, "type": "filament"},
    {"name": "CfA2 Great Wall", "ra": 195, "dec": 30, "type": "filament"},
    {"name": "Boötes Void", "ra": 218, "dec": 46, "type": "void"},
    {"name": "Coma Cluster", "ra": 194.9, "dec": 27.9, "type": "cluster"},
    {"name": "Virgo Cluster", "ra": 187.7, "dec": 12.3, "type": "cluster"},
    {"name": "Hercules Supercluster", "ra": 241, "dec": 17, "type": "cluster"},
    {"name": "Shapley Supercluster", "ra": 202, "dec": -31, "type": "cluster"},
    {"name": "Giant GRB Ring", "ra": 269, "dec": 62, "type": "ring"},
    {"name": "Huge-LQG (quasar group)", "ra": 162, "dec": 13, "type": "overdensity"},
]

if anomalies:
    print(f"\n  Checking {len(anomalies)} anomalous patches against known structures:")
    
    for anom in sorted(anomalies, key=lambda x: -abs(x['sigma'])):
        nearby = []
        for struct in known_structures:
            dra = min(abs(anom['ra_center'] - struct['ra']), 
                      360 - abs(anom['ra_center'] - struct['ra']))
            ddec = abs(anom['dec_center'] - struct['dec'])
            sep = np.sqrt(dra**2 + ddec**2)
            if sep < 30:
                nearby.append(f"{struct['name']} ({struct['type']}, {sep:.0f}°)")
        
        near_str = ", ".join(nearby) if nearby else "no known structures"
        print(f"  {anom['patch']:>30s} ({anom['sigma']:+.1f}σ {anom['type'][:7]}): {near_str}")

# ============================================================
# SAVE RESULTS
# ============================================================
import os
os.makedirs(RESULTS_DIR, exist_ok=True)

output = {
    'n_quasars': int(good.sum()),
    'n_patches': len(patch_results),
    'n_anomalous': len(anomalies),
    'patches': patch_results,
    'anomalies': anomalies,
}

with open(f'{RESULTS_DIR}/defect_map.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n  Results saved to {RESULTS_DIR}/defect_map.json")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("CRYSTAL DEFECT MAP — SUMMARY")
print("=" * 70)

print(f"""
  Total quasars analyzed: {good.sum()}
  Sky patches: {len(patch_results)}
  Anomalous patches (|σ| > 2.5): {len(anomalies)}
  
  DEFECT patches (correlation preserved): {sum(1 for a in anomalies if a['delta_rho'] > 0)}
  STRESS patches (excess degradation):    {sum(1 for a in anomalies if a['delta_rho'] < 0)}
  
  If the crystal model is correct:
  - DEFECT patches should align with known VOIDS (thinner crystal)
  - STRESS patches should align with known OVERDENSITIES (strained crystal)
  
  ...OR the opposite, depending on whether density PROTECTS or STRESSES.
  The data decides. We don't force it.
""")
