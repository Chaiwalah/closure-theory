#!/usr/bin/env python3
"""
MEMBRANE MODEL — PANTHEON+ SIGMOID THRESHOLD TEST (Test 6)
============================================================
Split Pantheon+ by host environment proxy and test whether 
the sigmoid threshold z₀ shifts with density.

Membrane prediction:
  Dense environments: z₀ ≈ 0.87-0.95 (gravity protects)
  Void environments:  z₀ ≈ 0.72-0.77 (free expansion decoheres)

Environment proxies available in Pantheon+:
  - HOST_LOGMASS: stellar mass of host galaxy (high mass ≈ denser environment)
  - MWEBV: Milky Way E(B-V) (higher = more foreground structure)
  - Galactic latitude from RA/DEC (high |b| = cleaner void-like sightline)
  
Also testing: color-distance coupling vs redshift in each subsample
to see if the sigmoid shape shifts.

ANOMALY TRACKING: Log anything weird.
"""

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# LOAD PANTHEON+
# ============================================================
def load_pantheon():
    data = {'z': [], 'c': [], 'c_err': [], 'x1': [], 'x1_err': [], 
            'mu': [], 'mu_err': [], 'mB': [], 'ra': [], 'dec': [],
            'logmass': [], 'ebv': []}
    
    with open('data/pantheon_plus.dat', 'r') as f:
        header = f.readline().split()
        for line in f:
            parts = line.split()
            if len(parts) < len(header):
                continue
            row = {header[i]: parts[i] for i in range(len(header))}
            try:
                z = float(row['zHD'])
                c = float(row['c'])
                c_err = float(row['cERR'])
                x1 = float(row['x1'])
                x1_err = float(row['x1ERR'])
                mu = float(row['m_b_corr'])
                mu_err = float(row['m_b_corr_err_DIAG'])
                mB = float(row['mB'])
                ra = float(row['RA'])
                dec = float(row['DEC'])
                logmass = float(row.get('HOST_LOGMASS', '-999'))
                ebv = float(row.get('MWEBV', '0'))
                
                if z > 0.01 and abs(c) < 0.5 and abs(x1) < 5 and c_err > 0 and c_err < 1:
                    data['z'].append(z)
                    data['c'].append(c)
                    data['c_err'].append(c_err)
                    data['x1'].append(x1)
                    data['x1_err'].append(x1_err)
                    data['mu'].append(mu)
                    data['mu_err'].append(mu_err)
                    data['mB'].append(mB)
                    data['ra'].append(ra)
                    data['dec'].append(dec)
                    data['logmass'].append(logmass)
                    data['ebv'].append(ebv)
            except (ValueError, KeyError):
                continue
    
    for k in data:
        data[k] = np.array(data[k])
    
    return data

def galactic_lat(ra_deg, dec_deg):
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    ra_ngp = np.radians(192.85948)
    dec_ngp = np.radians(27.12825)
    sin_b = (np.sin(dec_ngp) * np.sin(dec) + 
             np.cos(dec_ngp) * np.cos(dec) * np.cos(ra - ra_ngp))
    return np.degrees(np.arcsin(np.clip(sin_b, -1, 1)))

def sigmoid(z, A, z0, w, C):
    """Sigmoid: A / (1 + exp(-(z - z0) / w)) + C"""
    return A / (1.0 + np.exp(-(z - z0) / w)) + C

def measure_coupling_binned(z, c, mu, n_bins=6, label=""):
    """Measure color-distance coupling in redshift bins."""
    z_edges = np.percentile(z, np.linspace(0, 100, n_bins + 1))
    z_mids = []
    rhos = []
    pvals = []
    ns = []
    
    for i in range(n_bins):
        mask = (z >= z_edges[i]) & (z < z_edges[i+1])
        if mask.sum() < 20:
            continue
        r, p = stats.pearsonr(c[mask], mu[mask])
        z_mids.append(np.median(z[mask]))
        rhos.append(r)
        pvals.append(p)
        ns.append(mask.sum())
    
    return np.array(z_mids), np.array(rhos), np.array(pvals), np.array(ns)

def fit_sigmoid(z_mids, rhos, label=""):
    """Fit sigmoid to coupling vs z. Return z₀ and fit quality."""
    try:
        # Initial guess: coupling goes from ~0 to ~-0.3
        popt, pcov = curve_fit(sigmoid, z_mids, rhos, 
                               p0=[-0.3, 0.8, 0.2, 0.0],
                               bounds=([-1, 0.1, 0.01, -0.5], [0.5, 2.0, 1.0, 0.5]),
                               maxfev=10000)
        residuals = rhos - sigmoid(z_mids, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((rhos - np.mean(rhos))**2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return popt, r_squared
    except Exception as e:
        return None, 0

# ============================================================
# MAIN
# ============================================================
print("=" * 70)
print("MEMBRANE MODEL — PANTHEON+ ENVIRONMENT SPLIT TEST")
print("=" * 70)

data = load_pantheon()
N = len(data['z'])
print(f"\nLoaded {N} SNe Ia (z > 0.01, quality cuts applied)")

# Compute galactic latitude
gal_b = galactic_lat(data['ra'], data['dec'])

anomalies = []

# ============================================================
# SPLIT 1: HOST STELLAR MASS
# ============================================================
print("\n" + "=" * 70)
print("SPLIT 1: Host Stellar Mass (logM > 10 vs logM < 10)")
print("=" * 70)

valid_mass = data['logmass'] > 0  # -999 = missing
n_valid = valid_mass.sum()
print(f"  SNe with valid host mass: {n_valid}/{N}")

if n_valid > 100:
    mass_med = np.median(data['logmass'][valid_mass])
    # Standard mass step split at 10
    hi_mass = valid_mass & (data['logmass'] >= 10.0)
    lo_mass = valid_mass & (data['logmass'] < 10.0)
    
    print(f"  High-mass hosts (logM ≥ 10): {hi_mass.sum()}")
    print(f"  Low-mass hosts (logM < 10):  {lo_mass.sum()}")
    print(f"  Median logM: {mass_med:.2f}")
    
    for label, mask in [("HIGH-MASS", hi_mass), ("LOW-MASS", lo_mass)]:
        z_m, rho_m, p_m, n_m = measure_coupling_binned(
            data['z'][mask], data['c'][mask], data['mu'][mask], n_bins=5, label=label)
        
        print(f"\n  {label} — Color-Distance coupling by z:")
        for i in range(len(z_m)):
            sig = "***" if p_m[i] < 0.001 else "**" if p_m[i] < 0.01 else "*" if p_m[i] < 0.05 else ""
            print(f"    z={z_m[i]:.3f}: ρ={rho_m[i]:+.4f}  p={p_m[i]:.4f}  N={n_m[i]}  {sig}")
        
        popt, r2 = fit_sigmoid(z_m, rho_m, label)
        if popt is not None:
            print(f"  Sigmoid fit: z₀={popt[1]:.3f}, A={popt[0]:.3f}, w={popt[2]:.3f}, R²={r2:.3f}")
        else:
            print(f"  Sigmoid fit: FAILED")
            anomalies.append(f"SPLIT1-{label}: Sigmoid fit failed")
    
    # Compare z₀ between splits
    popt_hi, r2_hi = fit_sigmoid(*measure_coupling_binned(
        data['z'][hi_mass], data['c'][hi_mass], data['mu'][hi_mass], n_bins=5)[:2])
    popt_lo, r2_lo = fit_sigmoid(*measure_coupling_binned(
        data['z'][lo_mass], data['c'][lo_mass], data['mu'][lo_mass], n_bins=5)[:2])
    
    if popt_hi is not None and popt_lo is not None:
        delta_z0 = popt_hi[1] - popt_lo[1]
        print(f"\n  *** z₀ COMPARISON ***")
        print(f"  High-mass z₀ = {popt_hi[1]:.3f} (R²={r2_hi:.3f})")
        print(f"  Low-mass  z₀ = {popt_lo[1]:.3f} (R²={r2_lo:.3f})")
        print(f"  Δz₀ = {delta_z0:+.3f}")
        print(f"  Membrane prediction: Δz₀ > 0 (dense environments protect)")
        
        if delta_z0 > 0:
            print(f"  ✓ CONSISTENT — high-mass hosts push threshold higher")
        else:
            print(f"  ✗ OPPOSITE — high-mass hosts degrade FASTER")
            anomalies.append(f"SPLIT1: Δz₀={delta_z0:.3f} OPPOSITE to prediction")
    
    # Also check: coupling STRENGTH at matched z
    print(f"\n  Coupling strength at matched redshift:")
    z_edges_matched = [0.01, 0.3, 0.6, 0.9, 1.2]
    for i in range(len(z_edges_matched)-1):
        zlo, zhi = z_edges_matched[i], z_edges_matched[i+1]
        m_hi = hi_mass & (data['z'] >= zlo) & (data['z'] < zhi)
        m_lo = lo_mass & (data['z'] >= zlo) & (data['z'] < zhi)
        
        if m_hi.sum() > 15 and m_lo.sum() > 15:
            r_hi, p_hi = stats.pearsonr(data['c'][m_hi], data['mu'][m_hi])
            r_lo, p_lo = stats.pearsonr(data['c'][m_lo], data['mu'][m_lo])
            delta = r_hi - r_lo
            print(f"    z=[{zlo:.1f},{zhi:.1f}]: ρ_hi={r_hi:+.4f} ρ_lo={r_lo:+.4f} Δρ={delta:+.4f}  "
                  f"(N_hi={m_hi.sum()}, N_lo={m_lo.sum()})")
            
            if abs(delta) > 0.1:
                anomalies.append(f"SPLIT1 z=[{zlo},{zhi}]: Large Δρ={delta:.3f}")

# ============================================================
# SPLIT 2: GALACTIC LATITUDE (sightline proxy)
# ============================================================
print("\n" + "=" * 70)
print("SPLIT 2: Galactic Latitude (|b| > 45° vs |b| < 30°)")
print("=" * 70)

hi_lat = np.abs(gal_b) > 45  # clean sightlines (less foreground)
lo_lat = np.abs(gal_b) < 30  # dirty sightlines (more foreground)

print(f"  High |b| > 45° (clean sightlines): {hi_lat.sum()}")
print(f"  Low |b| < 30° (dirty sightlines):  {lo_lat.sum()}")

for label, mask in [("HIGH-LAT (clean)", hi_lat), ("LOW-LAT (dirty)", lo_lat)]:
    z_m, rho_m, p_m, n_m = measure_coupling_binned(
        data['z'][mask], data['c'][mask], data['mu'][mask], n_bins=5, label=label)
    
    print(f"\n  {label}:")
    for i in range(len(z_m)):
        sig = "***" if p_m[i] < 0.001 else "**" if p_m[i] < 0.01 else "*" if p_m[i] < 0.05 else ""
        print(f"    z={z_m[i]:.3f}: ρ={rho_m[i]:+.4f}  p={p_m[i]:.4f}  N={n_m[i]}  {sig}")

popt_hiL, r2_hiL = fit_sigmoid(*measure_coupling_binned(
    data['z'][hi_lat], data['c'][hi_lat], data['mu'][hi_lat], n_bins=5)[:2])
popt_loL, r2_loL = fit_sigmoid(*measure_coupling_binned(
    data['z'][lo_lat], data['c'][lo_lat], data['mu'][lo_lat], n_bins=5)[:2])

if popt_hiL is not None and popt_loL is not None:
    delta_z0_lat = popt_loL[1] - popt_hiL[1]  # Note: LOW lat = MORE foreground = denser
    print(f"\n  *** z₀ COMPARISON ***")
    print(f"  Low-lat (dirty/dense) z₀ = {popt_loL[1]:.3f} (R²={r2_loL:.3f})")
    print(f"  High-lat (clean/void) z₀ = {popt_hiL[1]:.3f} (R²={r2_hiL:.3f})")
    print(f"  Δz₀ = {delta_z0_lat:+.3f}")
    print(f"  Membrane prediction: dirty sightlines (more structure) → higher z₀")
    
    if delta_z0_lat > 0:
        print(f"  ✓ CONSISTENT")
    else:
        print(f"  ✗ OPPOSITE")
        anomalies.append(f"SPLIT2: Δz₀_lat={delta_z0_lat:.3f} OPPOSITE")

# ============================================================
# SPLIT 3: MW E(B-V) (foreground dust = structure proxy)
# ============================================================
print("\n" + "=" * 70)
print("SPLIT 3: MW E(B-V) (high dust = more foreground structure)")
print("=" * 70)

ebv_med = np.median(data['ebv'])
hi_ebv = data['ebv'] > ebv_med
lo_ebv = data['ebv'] <= ebv_med

print(f"  Median E(B-V) = {ebv_med:.4f}")
print(f"  High E(B-V) (more structure): {hi_ebv.sum()}")
print(f"  Low E(B-V) (less structure):  {lo_ebv.sum()}")

for label, mask in [("HIGH-EBV (dense)", hi_ebv), ("LOW-EBV (clean)", lo_ebv)]:
    z_m, rho_m, p_m, n_m = measure_coupling_binned(
        data['z'][mask], data['c'][mask], data['mu'][mask], n_bins=5, label=label)
    
    print(f"\n  {label}:")
    for i in range(len(z_m)):
        sig = "***" if p_m[i] < 0.001 else "**" if p_m[i] < 0.01 else "*" if p_m[i] < 0.05 else ""
        print(f"    z={z_m[i]:.3f}: ρ={rho_m[i]:+.4f}  p={p_m[i]:.4f}  N={n_m[i]}  {sig}")

popt_hiE, r2_hiE = fit_sigmoid(*measure_coupling_binned(
    data['z'][hi_ebv], data['c'][hi_ebv], data['mu'][hi_ebv], n_bins=5)[:2])
popt_loE, r2_loE = fit_sigmoid(*measure_coupling_binned(
    data['z'][lo_ebv], data['c'][lo_ebv], data['mu'][lo_ebv], n_bins=5)[:2])

if popt_hiE is not None and popt_loE is not None:
    delta_z0_ebv = popt_hiE[1] - popt_loE[1]
    print(f"\n  *** z₀ COMPARISON ***")
    print(f"  High E(B-V) (dense) z₀ = {popt_hiE[1]:.3f} (R²={r2_hiE:.3f})")
    print(f"  Low E(B-V) (clean)  z₀ = {popt_loE[1]:.3f} (R²={r2_loE:.3f})")
    print(f"  Δz₀ = {delta_z0_ebv:+.3f}")
    print(f"  Membrane prediction: high E(B-V) (denser foreground) → higher z₀")
    
    if delta_z0_ebv > 0:
        print(f"  ✓ CONSISTENT")
    else:
        print(f"  ✗ OPPOSITE")
        anomalies.append(f"SPLIT3: Δz₀_ebv={delta_z0_ebv:.3f} OPPOSITE")

# ============================================================
# FULL SAMPLE — BASELINE SIGMOID
# ============================================================
print("\n" + "=" * 70)
print("BASELINE: Full Sample Sigmoid")
print("=" * 70)

z_all, rho_all, p_all, n_all = measure_coupling_binned(
    data['z'], data['c'], data['mu'], n_bins=8)

print(f"\n  Full sample color-distance coupling:")
for i in range(len(z_all)):
    sig = "***" if p_all[i] < 0.001 else "**" if p_all[i] < 0.01 else "*" if p_all[i] < 0.05 else ""
    print(f"    z={z_all[i]:.3f}: ρ={rho_all[i]:+.4f}  p={p_all[i]:.4f}  N={n_all[i]}  {sig}")

popt_all, r2_all = fit_sigmoid(z_all, rho_all)
if popt_all is not None:
    print(f"\n  Sigmoid fit: z₀={popt_all[1]:.3f}, A={popt_all[0]:.3f}, w={popt_all[2]:.3f}, R²={r2_all:.3f}")

# ============================================================
# STRETCH (x1) — CONTROL (should NOT show environment split)
# ============================================================
print("\n" + "=" * 70)
print("CONTROL: Stretch-Distance Coupling (should be environment-independent)")
print("=" * 70)

if n_valid > 100:
    for label, mask in [("HIGH-MASS", hi_mass), ("LOW-MASS", lo_mass)]:
        z_m, rho_m, p_m, n_m = measure_coupling_binned(
            data['z'][mask], data['x1'][mask], data['mu'][mask], n_bins=5, label=label)
        
        print(f"\n  {label} — Stretch-Distance coupling:")
        for i in range(len(z_m)):
            print(f"    z={z_m[i]:.3f}: ρ={rho_m[i]:+.4f}  p={p_m[i]:.4f}  N={n_m[i]}")
        
        mean_rho = np.mean(np.abs(rho_m))
        print(f"  Mean |ρ| = {mean_rho:.4f}")
        if mean_rho > 0.1:
            anomalies.append(f"CONTROL-{label}: Stretch shows environment coupling |ρ|={mean_rho:.3f}")

# ============================================================
# ADDITIONAL: Color-color coupling and scatter analysis
# ============================================================
print("\n" + "=" * 70)
print("BONUS: Color Scatter by Environment × Redshift")
print("=" * 70)

if n_valid > 100:
    z_edges_bonus = [0.01, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5]
    print(f"\n  Color scatter (σ_c) by host mass and z:")
    print(f"  {'z_bin':>10s}  {'σ_hi':>8s}  {'σ_lo':>8s}  {'ratio':>8s}  {'N_hi':>6s}  {'N_lo':>6s}")
    
    for i in range(len(z_edges_bonus)-1):
        zlo, zhi = z_edges_bonus[i], z_edges_bonus[i+1]
        m_hi = hi_mass & (data['z'] >= zlo) & (data['z'] < zhi)
        m_lo = lo_mass & (data['z'] >= zlo) & (data['z'] < zhi)
        
        if m_hi.sum() > 10 and m_lo.sum() > 10:
            sig_hi = np.std(data['c'][m_hi])
            sig_lo = np.std(data['c'][m_lo])
            ratio = sig_hi / sig_lo if sig_lo > 0 else float('inf')
            print(f"  [{zlo:.1f},{zhi:.1f}]  {sig_hi:8.4f}  {sig_lo:8.4f}  {ratio:8.3f}  {m_hi.sum():6d}  {m_lo.sum():6d}")
            
            if ratio > 1.3 or ratio < 0.7:
                anomalies.append(f"SCATTER z=[{zlo},{zhi}]: σ_hi/σ_lo={ratio:.2f}")

# ============================================================
# ANOMALY LOG
# ============================================================
print("\n" + "=" * 70)
print("ANOMALY LOG")
print("=" * 70)

if anomalies:
    for i, a in enumerate(anomalies):
        print(f"  [{i+1}] {a}")
else:
    print("  No anomalies detected.")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

results = {
    "n_total": int(N),
    "baseline_z0": float(popt_all[1]) if popt_all is not None else None,
    "baseline_R2": float(r2_all),
    "splits": {},
    "anomalies": anomalies
}

for name, popt_d, popt_c, r2_d, r2_c, desc in [
    ("host_mass", popt_hi, popt_lo, r2_hi, r2_lo, "High-mass vs Low-mass"),
    ("gal_lat", popt_loL, popt_hiL, r2_loL, r2_hiL, "Low-lat(dense) vs High-lat(clean)"),
    ("mw_ebv", popt_hiE, popt_loE, r2_hiE, r2_loE, "High-EBV(dense) vs Low-EBV(clean)"),
]:
    if popt_d is not None and popt_c is not None:
        dz = popt_d[1] - popt_c[1]
        results["splits"][name] = {
            "z0_dense": float(popt_d[1]),
            "z0_clean": float(popt_c[1]),
            "delta_z0": float(dz),
            "R2_dense": float(r2_d),
            "R2_clean": float(r2_c),
            "consistent": bool(dz > 0)
        }
        direction = "✓" if dz > 0 else "✗"
        print(f"  {desc:40s}: Δz₀={dz:+.3f} {direction}")

with open('results_membrane_pantheon/results.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\n  Results saved to results_membrane_pantheon/results.json")
