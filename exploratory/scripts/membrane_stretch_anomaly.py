#!/usr/bin/env python3
"""
STRETCH ANOMALY INVESTIGATION
================================
The membrane model predicts stretch (x1) should be IMMUNE to 
environment/path effects because it's kinematic (velocity structure).

But we found: low-mass hosts show stretch-distance coupling 
ρ = 0.264, p = 0.001 at low z.

Questions:
1. Is this real or a selection artifact?
2. Does it depend on redshift (like color does)?
3. Does it depend on environment?
4. Does stretch couple to color? (If so, it's not independent)
5. Is stretch truly kinematic, or does it leak diagnostic information?
6. Does the coupling pattern match color or differ?

LEAVE NO STONE UNTURNED.
"""

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

def load_pantheon():
    data = {'z': [], 'c': [], 'c_err': [], 'x1': [], 'x1_err': [], 
            'mu': [], 'mu_err': [], 'mB': [], 'ra': [], 'dec': [],
            'logmass': [], 'ebv': [], 'x0': [], 'survey': []}
    
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
                x0 = float(row.get('x0', '0'))
                survey = int(row.get('IDSURVEY', '0'))
                
                if z > 0.01 and abs(c) < 0.5 and abs(x1) < 5 and c_err > 0 and c_err < 1:
                    for k, v in [('z',z),('c',c),('c_err',c_err),('x1',x1),('x1_err',x1_err),
                                 ('mu',mu),('mu_err',mu_err),('mB',mB),('ra',ra),('dec',dec),
                                 ('logmass',logmass),('ebv',ebv),('x0',x0),('survey',survey)]:
                        data[k].append(v)
            except (ValueError, KeyError):
                continue
    
    for k in data:
        data[k] = np.array(data[k])
    return data

anomalies = []
data = load_pantheon()
N = len(data['z'])
print(f"Loaded {N} SNe Ia\n")

hi_mass = data['logmass'] >= 10.0
lo_mass = data['logmass'] < 10.0

# ============================================================
# TEST 1: Stretch-Distance coupling by redshift (fine bins)
# ============================================================
print("=" * 70)
print("TEST 1: Stretch-Distance Coupling by Redshift (Fine Bins)")
print("=" * 70)

z_edges = [0.01, 0.03, 0.05, 0.08, 0.12, 0.2, 0.3, 0.5, 0.8, 1.5]
print(f"\n  {'z_bin':>12s}  {'ρ(x1,μ)':>8s}  {'p':>10s}  {'N':>5s}  {'sig':>4s}")

for i in range(len(z_edges)-1):
    mask = (data['z'] >= z_edges[i]) & (data['z'] < z_edges[i+1])
    if mask.sum() < 15:
        continue
    r, p = stats.pearsonr(data['x1'][mask], data['mu'][mask])
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
    print(f"  [{z_edges[i]:.2f},{z_edges[i+1]:.2f}]  {r:+8.4f}  {p:10.6f}  {mask.sum():5d}  {sig}")

# ============================================================
# TEST 2: Stretch-Distance by host mass (fine bins)
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: Stretch-Distance by Host Mass (Fine Bins)")
print("=" * 70)

for label, mask_env in [("HIGH-MASS (≥10)", hi_mass), ("LOW-MASS (<10)", lo_mass)]:
    print(f"\n  {label}:")
    print(f"  {'z_bin':>12s}  {'ρ(x1,μ)':>8s}  {'p':>10s}  {'N':>5s}  {'sig':>4s}")
    
    for i in range(len(z_edges)-1):
        mask = mask_env & (data['z'] >= z_edges[i]) & (data['z'] < z_edges[i+1])
        if mask.sum() < 15:
            continue
        r, p = stats.pearsonr(data['x1'][mask], data['mu'][mask])
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"  [{z_edges[i]:.2f},{z_edges[i+1]:.2f}]  {r:+8.4f}  {p:10.6f}  {mask.sum():5d}  {sig}")
        
        if p < 0.01:
            anomalies.append(f"TEST2-{label[:4]} z=[{z_edges[i]:.2f},{z_edges[i+1]:.2f}]: ρ(x1,μ)={r:.3f}, p={p:.4f}")

# ============================================================
# TEST 3: Is stretch coupled to color? (Contamination check)
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: Stretch-Color Coupling (x1 vs c)")
print("=" * 70)
print("  If x1 correlates with c, stretch isn't purely kinematic")

print(f"\n  {'z_bin':>12s}  {'ρ(x1,c)':>8s}  {'p':>10s}  {'N':>5s}  {'sig':>4s}")

for i in range(len(z_edges)-1):
    mask = (data['z'] >= z_edges[i]) & (data['z'] < z_edges[i+1])
    if mask.sum() < 15:
        continue
    r, p = stats.pearsonr(data['x1'][mask], data['c'][mask])
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
    print(f"  [{z_edges[i]:.2f},{z_edges[i+1]:.2f}]  {r:+8.4f}  {p:10.6f}  {mask.sum():5d}  {sig}")

# Full sample
r_xc, p_xc = stats.pearsonr(data['x1'], data['c'])
print(f"\n  Full sample: ρ(x1,c) = {r_xc:+.4f}, p = {p_xc:.2e}")

if abs(r_xc) > 0.1:
    anomalies.append(f"TEST3: x1-c coupling ρ={r_xc:.3f} — stretch contaminated by color")

# ============================================================
# TEST 4: Partial correlation — x1-μ controlling for c
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: Partial Correlation x1-μ | c (removing color contamination)")
print("=" * 70)

def partial_corr(x, y, z_control):
    """Partial correlation of x,y controlling for z."""
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(z_control)
    x, y, z_c = x[valid], y[valid], z_control[valid]
    
    # Regress x on z, get residuals
    slope_xz, int_xz, _, _, _ = stats.linregress(z_c, x)
    res_x = x - (slope_xz * z_c + int_xz)
    
    # Regress y on z, get residuals
    slope_yz, int_yz, _, _, _ = stats.linregress(z_c, y)
    res_y = y - (slope_yz * z_c + int_yz)
    
    return stats.pearsonr(res_x, res_y)

print(f"\n  {'z_bin':>12s}  {'ρ(x1,μ)':>8s}  {'ρ(x1,μ|c)':>10s}  {'p_partial':>10s}  {'N':>5s}")

for i in range(len(z_edges)-1):
    mask = (data['z'] >= z_edges[i]) & (data['z'] < z_edges[i+1])
    if mask.sum() < 20:
        continue
    r_raw, p_raw = stats.pearsonr(data['x1'][mask], data['mu'][mask])
    r_part, p_part = partial_corr(data['x1'][mask], data['mu'][mask], data['c'][mask])
    print(f"  [{z_edges[i]:.2f},{z_edges[i+1]:.2f}]  {r_raw:+8.4f}  {r_part:+10.4f}  {p_part:10.6f}  {mask.sum():5d}")

# By mass
print(f"\n  Partial correlation by host mass:")
for label, mask_env in [("HIGH-MASS", hi_mass), ("LOW-MASS", lo_mass)]:
    r_part, p_part = partial_corr(data['x1'][mask_env], data['mu'][mask_env], data['c'][mask_env])
    r_raw, p_raw = stats.pearsonr(data['x1'][mask_env], data['mu'][mask_env])
    print(f"  {label:12s}: ρ_raw={r_raw:+.4f}  ρ_partial={r_part:+.4f}  p={p_part:.4f}")

# ============================================================
# TEST 5: Stretch distribution by environment
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: Stretch Distribution by Environment")
print("=" * 70)

for label, mask_env in [("HIGH-MASS", hi_mass), ("LOW-MASS", lo_mass)]:
    x1_vals = data['x1'][mask_env]
    print(f"\n  {label}: mean={np.mean(x1_vals):.3f} std={np.std(x1_vals):.3f} "
          f"skew={stats.skew(x1_vals):.3f} kurt={stats.kurtosis(x1_vals):.3f} N={mask_env.sum()}")

ks_stat, ks_p = stats.ks_2samp(data['x1'][hi_mass], data['x1'][lo_mass])
print(f"\n  KS test (hi vs lo mass x1 distribution): D={ks_stat:.4f}, p={ks_p:.4f}")

if ks_p < 0.01:
    anomalies.append(f"TEST5: x1 distribution differs by host mass (KS p={ks_p:.4f})")

# By z bins
print(f"\n  Mean x1 by mass × redshift:")
print(f"  {'z_bin':>12s}  {'x1_hi':>8s}  {'x1_lo':>8s}  {'Δx1':>8s}  {'p(MW)':>10s}")

for i in range(len(z_edges)-1):
    m_hi = hi_mass & (data['z'] >= z_edges[i]) & (data['z'] < z_edges[i+1])
    m_lo = lo_mass & (data['z'] >= z_edges[i]) & (data['z'] < z_edges[i+1])
    if m_hi.sum() > 10 and m_lo.sum() > 10:
        u_stat, u_p = stats.mannwhitneyu(data['x1'][m_hi], data['x1'][m_lo], alternative='two-sided')
        dx1 = np.mean(data['x1'][m_hi]) - np.mean(data['x1'][m_lo])
        print(f"  [{z_edges[i]:.2f},{z_edges[i+1]:.2f}]  {np.mean(data['x1'][m_hi]):+8.3f}  "
              f"{np.mean(data['x1'][m_lo]):+8.3f}  {dx1:+8.3f}  {u_p:10.4f}")

# ============================================================
# TEST 6: Stretch-distance coupling AFTER removing Malmquist bias
# ============================================================
print("\n" + "=" * 70)
print("TEST 6: Malmquist Bias Check")
print("=" * 70)
print("  At high z, surveys preferentially detect brighter (higher x1) SNe")
print("  This creates artificial x1-μ correlation from selection")

# Check: does x1 correlate with z? (Malmquist signature)
r_x1z, p_x1z = stats.pearsonr(data['x1'], data['z'])
print(f"\n  ρ(x1, z) = {r_x1z:+.4f}, p = {p_x1z:.2e}")

if abs(r_x1z) > 0.05 and p_x1z < 0.01:
    anomalies.append(f"TEST6: x1-z correlation ρ={r_x1z:.3f} — Malmquist bias present")

# By mass
for label, mask_env in [("HIGH-MASS", hi_mass), ("LOW-MASS", lo_mass)]:
    r, p = stats.pearsonr(data['x1'][mask_env], data['z'][mask_env])
    print(f"  {label}: ρ(x1,z) = {r:+.4f}, p = {p:.2e}")

# x1 vs z residual: remove the Malmquist trend, then re-test
print(f"\n  After detrending x1 vs z:")
slope_x1z, int_x1z, _, _, _ = stats.linregress(data['z'], data['x1'])
x1_detrended = data['x1'] - (slope_x1z * data['z'] + int_x1z)

for label, mask_env in [("FULL", np.ones(N, bool)), ("HIGH-MASS", hi_mass), ("LOW-MASS", lo_mass)]:
    r, p = stats.pearsonr(x1_detrended[mask_env], data['mu'][mask_env])
    print(f"  {label}: ρ(x1_detrended, μ) = {r:+.4f}, p = {p:.4f}")

# ============================================================
# TEST 7: Color-distance coupling comparison (sanity check)
# ============================================================
print("\n" + "=" * 70)
print("TEST 7: Color vs Stretch Coupling Patterns Compared")
print("=" * 70)
print("  Do they behave the SAME or DIFFERENTLY with z?")
print("  Same = stretch is contaminated. Different = stretch has its own physics.")

print(f"\n  {'z_bin':>12s}  {'ρ(c,μ)':>8s}  {'ρ(x1,μ)':>8s}  {'ratio':>8s}  {'N':>5s}")

ratios = []
z_mids_7 = []

for i in range(len(z_edges)-1):
    mask = (data['z'] >= z_edges[i]) & (data['z'] < z_edges[i+1])
    if mask.sum() < 20:
        continue
    r_c, _ = stats.pearsonr(data['c'][mask], data['mu'][mask])
    r_x1, _ = stats.pearsonr(data['x1'][mask], data['mu'][mask])
    ratio = r_x1 / r_c if abs(r_c) > 0.01 else float('nan')
    print(f"  [{z_edges[i]:.2f},{z_edges[i+1]:.2f}]  {r_c:+8.4f}  {r_x1:+8.4f}  {ratio:+8.3f}  {mask.sum():5d}")
    
    if np.isfinite(ratio):
        ratios.append(ratio)
        z_mids_7.append((z_edges[i] + z_edges[i+1])/2)

if len(ratios) > 3:
    r_ratio_z, p_ratio_z = stats.pearsonr(z_mids_7, ratios)
    print(f"\n  x1/c coupling ratio trend with z: r={r_ratio_z:+.4f}, p={p_ratio_z:.4f}")
    if abs(r_ratio_z) > 0.5 and p_ratio_z < 0.05:
        anomalies.append(f"TEST7: x1/c ratio evolves with z (r={r_ratio_z:.3f}) — different physics")
        print(f"  ✓ DIFFERENT PATTERNS — stretch and color evolve differently")
    else:
        print(f"  ~ Similar patterns — stretch may be tracking color leakage")

# ============================================================
# TEST 8: x1 error vs z (measurement artifact check)
# ============================================================
print("\n" + "=" * 70)
print("TEST 8: Measurement Quality Check")
print("=" * 70)

r_err_z, p_err_z = stats.pearsonr(data['x1_err'], data['z'])
print(f"  ρ(x1_err, z) = {r_err_z:+.4f}, p = {p_err_z:.2e}")
print(f"  ρ(c_err, z)  = {stats.pearsonr(data['c_err'], data['z'])[0]:+.4f}")

# Does x1_err predict x1-mu coupling? (artifact signature)
print(f"\n  x1 error terciles:")
err_t = np.percentile(data['x1_err'], [33, 66])
for label, lo, hi in [("Low err", 0, err_t[0]), ("Mid err", err_t[0], err_t[1]), ("High err", err_t[1], 10)]:
    mask = (data['x1_err'] >= lo) & (data['x1_err'] < hi)
    r, p = stats.pearsonr(data['x1'][mask], data['mu'][mask])
    print(f"  {label:10s}: ρ(x1,μ)={r:+.4f}  p={p:.4f}  N={mask.sum()}")

# ============================================================
# TEST 9: Survey-specific stretch behavior
# ============================================================
print("\n" + "=" * 70)
print("TEST 9: Survey-Specific Stretch Coupling")
print("=" * 70)

surveys = np.unique(data['survey'])
print(f"\n  {'Survey':>8s}  {'N':>6s}  {'<z>':>6s}  {'ρ(x1,μ)':>8s}  {'p':>10s}  {'ρ(c,μ)':>8s}")

for s in sorted(surveys):
    mask = data['survey'] == s
    if mask.sum() < 30:
        continue
    r_x1, p_x1 = stats.pearsonr(data['x1'][mask], data['mu'][mask])
    r_c, p_c = stats.pearsonr(data['c'][mask], data['mu'][mask])
    print(f"  {int(s):8d}  {mask.sum():6d}  {np.mean(data['z'][mask]):6.3f}  "
          f"{r_x1:+8.4f}  {p_x1:10.4f}  {r_c:+8.4f}")

# ============================================================
# TEST 10: The killer — does x1 carry diagnostic information?
# ============================================================
print("\n" + "=" * 70)
print("TEST 10: Does Stretch Carry Diagnostic Information?")
print("=" * 70)
print("  If x1 responds to host metallicity/age, it's not purely kinematic")

# x1 vs host mass (proxy for metallicity)
r_x1m, p_x1m = stats.pearsonr(data['x1'], data['logmass'])
print(f"\n  ρ(x1, logM_host) = {r_x1m:+.4f}, p = {p_x1m:.2e}")

# Does this coupling evolve with z?
print(f"\n  x1-logM coupling by redshift:")
print(f"  {'z_bin':>12s}  {'ρ(x1,M)':>8s}  {'p':>10s}  {'N':>5s}")

for i in range(len(z_edges)-1):
    mask = (data['z'] >= z_edges[i]) & (data['z'] < z_edges[i+1])
    if mask.sum() < 20:
        continue
    r, p = stats.pearsonr(data['x1'][mask], data['logmass'][mask])
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
    print(f"  [{z_edges[i]:.2f},{z_edges[i+1]:.2f}]  {r:+8.4f}  {p:10.4f}  {mask.sum():5d}  {sig}")

if abs(r_x1m) > 0.1 and p_x1m < 0.01:
    anomalies.append(f"TEST10: x1 couples to host mass ρ={r_x1m:.3f} — NOT purely kinematic")
    print(f"\n  ⚠️ STRETCH IS NOT PURELY KINEMATIC — it carries host information")
else:
    print(f"\n  ✓ Stretch appears kinematically clean (weak host coupling)")

# ============================================================
# ANOMALY LOG
# ============================================================
print("\n" + "=" * 70)
print("ANOMALY LOG — STRETCH INVESTIGATION")
print("=" * 70)

for i, a in enumerate(anomalies):
    print(f"  [{i+1}] {a}")

if not anomalies:
    print("  No anomalies detected.")

# ============================================================
# VERDICT
# ============================================================
print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)

print("""
  The stretch anomaly investigation reveals whether x1 is truly 
  kinematic (membrane-immune) or carries diagnostic information 
  (membrane-susceptible).
  
  Key questions answered:
  - Is x1-μ coupling real or selection artifact?
  - Does x1 carry host/environment information?
  - Does x1 behave like color (same physics) or differently?
  - Does removing color contamination kill the coupling?
""")
