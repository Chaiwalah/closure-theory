#!/usr/bin/env python3
"""
CLOSURE THEORY — FRB CROSS-DOMAIN TEST v3
535 CHIME FRBs + 186 localized FRBs from FRBs/FRB repo.

Uses DM as distance proxy for CHIME (no host z).
Uses actual host z for localized sample.
Tests paired observable correlations for closure signatures.
"""

import json, os, glob
import numpy as np
from scipy import stats, optimize
import warnings
warnings.filterwarnings('ignore')

print("=" * 60)
print("CLOSURE THEORY — FRB CROSS-DOMAIN TEST v3")
print("=" * 60)

# ============================================================
# LOAD CHIME CATALOG (535 FRBs)
# ============================================================
chime_path = '/root/closure-bayesn/data/FRB_repo/frb/data/FRBs/CHIME_catalog-2021-1-27.json'
with open(chime_path) as f:
    chime = json.load(f)

print(f"\n  CHIME Catalog: {len(chime)} FRBs")

def sf(val, default=np.nan):
    if val is None: return default
    try:
        v = float(val)
        return v if v > -900 else default
    except: return default

c_dm = np.array([sf(r.get('fitburst_dm')) for r in chime])
c_dm_exc = np.array([sf(r.get('dm_excess_ne2001')) for r in chime])
c_scat = np.array([sf(r.get('scattering_time_ms')) for r in chime])
c_width = np.array([sf(r.get('pulse_width_ms')) for r in chime])
c_fluence = np.array([sf(r.get('fluence')) for r in chime])
c_flux = np.array([sf(r.get('flux')) for r in chime])
c_sp_idx = np.array([sf(r.get('spectral_index')) for r in chime])
c_sp_run = np.array([sf(r.get('spectral_running')) for r in chime])
c_peak_freq = np.array([sf(r.get('peak_frequency')) for r in chime])
c_snr = np.array([sf(r.get('fitburst_snr')) for r in chime])
c_repeater = np.array([r.get('repeater_of','') not in [None, '', 'None'] for r in chime])

v = lambda x: ~np.isnan(x) & (x != 0)
print(f"  Valid DM: {v(c_dm).sum()}")
print(f"  Valid DM_exc: {v(c_dm_exc).sum()}")
print(f"  Valid scattering: {(~np.isnan(c_scat) & (c_scat > 0)).sum()}")
print(f"  Valid width: {(~np.isnan(c_width) & (c_width > 0)).sum()}")
print(f"  Valid fluence: {(~np.isnan(c_fluence) & (c_fluence > 0)).sum()}")
print(f"  Valid spectral_idx: {v(c_sp_idx).sum()}")
print(f"  Repeaters: {c_repeater.sum()}")
print(f"  DM range: {np.nanmin(c_dm):.0f} — {np.nanmax(c_dm):.0f}")
print(f"  DM_exc range: {np.nanmin(c_dm_exc):.0f} — {np.nanmax(c_dm_exc):.0f}")

# ============================================================
# LOAD LOCALIZED FRBs (with host z + RM)
# ============================================================
frb_dir = '/root/closure-bayesn/data/FRB_repo/frb/data/FRBs'
loc_frbs = []
for fname in sorted(os.listdir(frb_dir)):
    if not fname.startswith('FRB') or not fname.endswith('.json'): continue
    with open(os.path.join(frb_dir, fname)) as f:
        d = json.load(f)
    
    def gv(obj, key):
        val = obj.get(key)
        if val is None: return np.nan
        if isinstance(val, dict): return sf(val.get('value'))
        return sf(val)
    
    dm = gv(d, 'DM')
    if np.isnan(dm): continue
    
    pulse = d.get('pulse', {}) or {}
    loc_frbs.append({
        'name': d.get('FRB',''), 'dm': dm, 'dm_ism': gv(d, 'DMISM'),
        'z': gv(d, 'z'), 'rm': gv(d, 'RM'),
        'width': gv(pulse, 'Wi'), 'tscatt': gv(pulse, 'tscatt'),
        'repeater': d.get('repeater', False)
    })

print(f"\n  Localized FRBs: {len(loc_frbs)}")
l_z = np.array([f['z'] for f in loc_frbs])
l_rm = np.array([f['rm'] for f in loc_frbs])
l_dm = np.array([f['dm'] for f in loc_frbs])
l_dm_ism = np.array([f['dm_ism'] for f in loc_frbs])
l_dm_exc = l_dm - l_dm_ism - 65
vz = ~np.isnan(l_z) & (l_z > 0)
vrm = ~np.isnan(l_rm)
print(f"  With host z: {vz.sum()}")
print(f"  With RM: {vrm.sum()}")

# ============================================================
# TEST C1: DM vs Scattering — Correlation in DM bins (CHIME)
# τ ∝ DM^α is the standard model. If closure is real, the
# power-law index or scatter around it should change at high DM.
# ============================================================
print("\n" + "=" * 60)
print("[1/6] Test C1: DM vs Scattering Time — 535 CHIME FRBs")
print("=" * 60)

mask_c1 = ~np.isnan(c_dm) & ~np.isnan(c_scat) & (c_scat > 0) & (c_dm > 0)
print(f"\n  Pairs available: {mask_c1.sum()}")

if mask_c1.sum() > 30:
    ld = np.log10(c_dm[mask_c1])
    ls = np.log10(c_scat[mask_c1])
    
    # Overall power-law
    slope, intercept, r_all, p_all, se = stats.linregress(ld, ls)
    print(f"  Overall: τ ∝ DM^{slope:.2f}, r={r_all:.3f}, p={p_all:.2e}")
    
    # By DM bins
    dm_edges = [100, 250, 400, 600, 900, 1500, 3000]
    print(f"\n  Closure predicts: power-law slope or r CHANGES at high DM")
    corr_vals = []
    slope_vals = []
    dm_mids = []
    for i in range(len(dm_edges)-1):
        m = mask_c1 & (c_dm >= dm_edges[i]) & (c_dm < dm_edges[i+1])
        if m.sum() >= 8:
            s, inter, r, p, _ = stats.linregress(np.log10(c_dm[m]), np.log10(c_scat[m]))
            corr_vals.append(r)
            slope_vals.append(s)
            dm_mids.append((dm_edges[i]+dm_edges[i+1])/2)
            print(f"  DM=[{dm_edges[i]},{dm_edges[i+1]}): "
                  f"slope={s:.2f}, r={r:.3f}, p={p:.4f}, N={m.sum()}")
    
    if len(dm_mids) >= 3:
        r_t, p_t = stats.spearmanr(dm_mids, corr_vals)
        print(f"\n  Correlation evolution with DM: r={r_t:.3f}, p={p_t:.4f}")
        r_s, p_s = stats.spearmanr(dm_mids, slope_vals)
        print(f"  Slope evolution with DM: r={r_s:.3f}, p={p_s:.4f}")

# ============================================================
# TEST C2: DM vs Width — CHIME
# ============================================================
print("\n" + "=" * 60)
print("[2/6] Test C2: DM vs Width — 535 CHIME FRBs")
print("=" * 60)

mask_c2 = ~np.isnan(c_dm) & ~np.isnan(c_width) & (c_width > 0) & (c_dm > 0)
print(f"\n  Pairs available: {mask_c2.sum()}")

if mask_c2.sum() > 30:
    dm_edges = [100, 250, 400, 600, 900, 1500, 3000]
    corr_vals2 = []
    dm_mids2 = []
    for i in range(len(dm_edges)-1):
        m = mask_c2 & (c_dm >= dm_edges[i]) & (c_dm < dm_edges[i+1])
        if m.sum() >= 8:
            r, p = stats.spearmanr(c_dm[m], c_width[m])
            corr_vals2.append(r)
            dm_mids2.append((dm_edges[i]+dm_edges[i+1])/2)
            print(f"  DM=[{dm_edges[i]},{dm_edges[i+1]}): r={r:.3f}, p={p:.4f}, N={m.sum()}")
    
    if len(dm_mids2) >= 3:
        r_t, p_t = stats.spearmanr(dm_mids2, corr_vals2)
        print(f"\n  Correlation evolution: r={r_t:.3f}, p={p_t:.4f}")

# ============================================================
# TEST C3: DM vs Spectral Index — CHIME
# Spectral index is intrinsic. Should be INDEPENDENT of DM.
# Closure: coupling emerges at high DM
# ============================================================
print("\n" + "=" * 60)
print("[3/6] Test C3: DM vs Spectral Index — CHIME")
print("=" * 60)

mask_c3 = ~np.isnan(c_dm) & ~np.isnan(c_sp_idx) & (c_dm > 0)
print(f"\n  Pairs available: {mask_c3.sum()}")

if mask_c3.sum() > 30:
    dm_edges = [100, 250, 400, 600, 900, 1500, 3000]
    for i in range(len(dm_edges)-1):
        m = mask_c3 & (c_dm >= dm_edges[i]) & (c_dm < dm_edges[i+1])
        if m.sum() >= 8:
            r, p = stats.spearmanr(c_dm[m], c_sp_idx[m])
            print(f"  DM=[{dm_edges[i]},{dm_edges[i+1]}): r={r:.3f}, p={p:.4f}, N={m.sum()}")

# ============================================================
# TEST C4: Observable Rank Compression — CHIME
# Multi-variate: [DM, τ, width, fluence, sp_idx]
# Closure: effective rank decreases at high DM
# ============================================================
print("\n" + "=" * 60)
print("[4/6] Test C4: Observable Rank Compression — CHIME")
print("=" * 60)

mask_c4 = (~np.isnan(c_dm) & ~np.isnan(c_scat) & ~np.isnan(c_width) 
           & ~np.isnan(c_fluence) & ~np.isnan(c_sp_idx)
           & (c_dm > 0) & (c_scat > 0) & (c_width > 0) & (c_fluence > 0))
print(f"\n  FRBs with all 5 observables: {mask_c4.sum()}")

if mask_c4.sum() > 30:
    dm_edges = [100, 350, 600, 1000, 3000]
    for i in range(len(dm_edges)-1):
        m = mask_c4 & (c_dm >= dm_edges[i]) & (c_dm < dm_edges[i+1])
        if m.sum() >= 10:
            data = np.column_stack([
                np.log10(c_dm[m]),
                np.log10(c_scat[m]),
                np.log10(c_width[m]),
                np.log10(c_fluence[m]),
                c_sp_idx[m]
            ])
            data = (data - data.mean(0)) / (data.std(0) + 1e-10)
            eigvals = np.linalg.eigvalsh(np.corrcoef(data.T))[::-1]
            eff_rank = np.sum(eigvals > 0.2)  # eigenvalue threshold
            var_top2 = eigvals[:2].sum() / eigvals.sum()
            mean_offdiag = (np.abs(np.corrcoef(data.T)).sum() - 5) / 20
            print(f"  DM=[{dm_edges[i]},{dm_edges[i+1]}): N={m.sum()}, "
                  f"eff_rank={eff_rank}/5, top2_var={var_top2:.3f}, "
                  f"mean|corr|={mean_offdiag:.3f}")

# ============================================================
# TEST C5: DM_exc vs |RM| by redshift — Localized FRBs
# THE MONEY PAIR: same electron column, different weighting
# ============================================================
print("\n" + "=" * 60)
print("[5/6] Test C5: DM_exc vs |RM| by Redshift — Localized FRBs")
print("=" * 60)

mask_c5 = vz & vrm & ~np.isnan(l_dm_exc)
print(f"\n  FRBs with z + RM: {mask_c5.sum()}")

if mask_c5.sum() > 10:
    r_all, p_all = stats.spearmanr(l_dm_exc[mask_c5], np.abs(l_rm[mask_c5]))
    print(f"  Overall: r={r_all:.3f}, p={p_all:.4f}")
    
    z_edges = [0, 0.1, 0.2, 0.35, 0.55, 1.5]
    corr_z = []
    z_mids = []
    for i in range(len(z_edges)-1):
        m = mask_c5 & (l_z >= z_edges[i]) & (l_z < z_edges[i+1])
        if m.sum() >= 3:
            r, p = stats.spearmanr(l_dm_exc[m], np.abs(l_rm[m]))
            corr_z.append(r)
            z_mids.append((z_edges[i]+z_edges[i+1])/2)
            print(f"  z=[{z_edges[i]},{z_edges[i+1]}): r={r:.3f}, p={p:.4f}, N={m.sum()}")
    
    if len(z_mids) >= 3:
        r_t, p_t = stats.spearmanr(z_mids, corr_z)
        print(f"\n  DM-RM correlation evolution with z: r={r_t:.3f}, p={p_t:.4f}")

# ============================================================
# TEST C6: Sigmoid Transition — Rolling Correlations (CHIME)
# ============================================================
print("\n" + "=" * 60)
print("[6/6] Test C6: Sigmoid Transition Detection — CHIME")
print("=" * 60)

def sigmoid(x, A, x0, k, B):
    return A / (1 + np.exp(-k * (x - x0))) + B

# DM vs Scattering rolling correlation
mask_roll = ~np.isnan(c_dm) & ~np.isnan(c_scat) & (c_scat > 0) & (c_dm > 0)
if mask_roll.sum() > 60:
    idx = np.argsort(c_dm[mask_roll])
    dm_s = c_dm[mask_roll][idx]
    sc_s = c_scat[mask_roll][idx]
    
    W = max(25, mask_roll.sum() // 12)
    roll_dm, roll_r = [], []
    for i in range(len(dm_s) - W):
        r, _ = stats.spearmanr(np.log10(dm_s[i:i+W]), np.log10(sc_s[i:i+W]))
        roll_dm.append(np.median(dm_s[i:i+W]))
        roll_r.append(r)
    
    roll_dm = np.array(roll_dm)
    roll_r = np.array(roll_r)
    log_roll_dm = np.log10(roll_dm)
    
    # Sigmoid fit
    try:
        popt, _ = optimize.curve_fit(sigmoid, log_roll_dm, roll_r,
                                      p0=[0.5, 2.7, 5, 0.2],
                                      bounds=([-2, 1.5, 0.1, -2], [2, 3.5, 50, 2]),
                                      maxfev=10000)
        dm_thresh = 10**popt[1]
        z_approx = (dm_thresh - 65) / 935
        
        slope, inter, r_lin, p_lin, _ = stats.linregress(log_roll_dm, roll_r)
        ss_sig = np.sum((roll_r - sigmoid(log_roll_dm, *popt))**2)
        ss_lin = np.sum((roll_r - (slope*log_roll_dm + inter))**2)
        n = len(roll_dm)
        F = ((ss_lin - ss_sig)/2) / (ss_sig/(n-4)) if ss_sig > 0 else 0
        
        print(f"\n  DM-Scattering rolling correlation (window={W}):")
        print(f"  Sigmoid DM₀ = {dm_thresh:.0f} pc/cm3 (≈ z={z_approx:.2f})")
        print(f"  k = {popt[2]:.1f}, A = {popt[0]:.3f}, B = {popt[3]:.3f}")
        print(f"  SS_sigmoid={ss_sig:.4f}, SS_linear={ss_lin:.4f}")
        print(f"  Sigmoid better: {ss_sig < ss_lin}, F={F:.2f}")
    except Exception as e:
        print(f"  DM-Scat sigmoid failed: {e}")

# DM vs Width rolling correlation
mask_roll2 = ~np.isnan(c_dm) & ~np.isnan(c_width) & (c_width > 0) & (c_dm > 0)
if mask_roll2.sum() > 60:
    idx2 = np.argsort(c_dm[mask_roll2])
    dm_s2 = c_dm[mask_roll2][idx2]
    w_s2 = c_width[mask_roll2][idx2]
    
    W2 = max(25, mask_roll2.sum() // 12)
    roll_dm2, roll_r2 = [], []
    for i in range(len(dm_s2) - W2):
        r, _ = stats.spearmanr(dm_s2[i:i+W2], w_s2[i:i+W2])
        roll_dm2.append(np.median(dm_s2[i:i+W2]))
        roll_r2.append(r)
    
    roll_dm2 = np.array(roll_dm2)
    roll_r2 = np.array(roll_r2)
    
    try:
        popt2, _ = optimize.curve_fit(sigmoid, np.log10(roll_dm2), roll_r2,
                                       p0=[0.5, 2.7, 5, 0.1],
                                       bounds=([-2, 1.5, 0.1, -2], [2, 3.5, 50, 2]),
                                       maxfev=10000)
        dm_thresh2 = 10**popt2[1]
        z_approx2 = (dm_thresh2 - 65) / 935
        
        slope2, inter2, _, _, _ = stats.linregress(np.log10(roll_dm2), roll_r2)
        ss_sig2 = np.sum((roll_r2 - sigmoid(np.log10(roll_dm2), *popt2))**2)
        ss_lin2 = np.sum((roll_r2 - (slope2*np.log10(roll_dm2) + inter2))**2)
        F2 = ((ss_lin2 - ss_sig2)/2) / (ss_sig2/(len(roll_dm2)-4)) if ss_sig2 > 0 else 0
        
        print(f"\n  DM-Width rolling correlation (window={W2}):")
        print(f"  Sigmoid DM₀ = {dm_thresh2:.0f} pc/cm3 (≈ z={z_approx2:.2f})")
        print(f"  k = {popt2[2]:.1f}")
        print(f"  SS_sigmoid={ss_sig2:.4f}, SS_linear={ss_lin2:.4f}")
        print(f"  Sigmoid better: {ss_sig2 < ss_lin2}, F={F2:.2f}")
    except Exception as e:
        print(f"  DM-Width sigmoid failed: {e}")

# DM vs Fluence rolling correlation
mask_roll3 = ~np.isnan(c_dm) & ~np.isnan(c_fluence) & (c_fluence > 0) & (c_dm > 0)
if mask_roll3.sum() > 60:
    idx3 = np.argsort(c_dm[mask_roll3])
    dm_s3 = c_dm[mask_roll3][idx3]
    fl_s3 = c_fluence[mask_roll3][idx3]
    
    W3 = max(25, mask_roll3.sum() // 12)
    roll_dm3, roll_r3 = [], []
    for i in range(len(dm_s3) - W3):
        r, _ = stats.spearmanr(dm_s3[i:i+W3], np.log10(fl_s3[i:i+W3]))
        roll_dm3.append(np.median(dm_s3[i:i+W3]))
        roll_r3.append(r)
    
    roll_dm3 = np.array(roll_dm3)
    roll_r3 = np.array(roll_r3)
    
    try:
        popt3, _ = optimize.curve_fit(sigmoid, np.log10(roll_dm3), roll_r3,
                                       p0=[-0.3, 2.7, 5, 0],
                                       bounds=([-2, 1.5, 0.1, -2], [2, 3.5, 50, 2]),
                                       maxfev=10000)
        dm_thresh3 = 10**popt3[1]
        z_approx3 = (dm_thresh3 - 65) / 935
        
        slope3, inter3, _, _, _ = stats.linregress(np.log10(roll_dm3), roll_r3)
        ss_sig3 = np.sum((roll_r3 - sigmoid(np.log10(roll_dm3), *popt3))**2)
        ss_lin3 = np.sum((roll_r3 - (slope3*np.log10(roll_dm3) + inter3))**2)
        F3 = ((ss_lin3 - ss_sig3)/2) / (ss_sig3/(len(roll_dm3)-4)) if ss_sig3 > 0 else 0
        
        print(f"\n  DM-Fluence rolling correlation (window={W3}):")
        print(f"  Sigmoid DM₀ = {dm_thresh3:.0f} pc/cm3 (≈ z={z_approx3:.2f})")
        print(f"  k = {popt3[2]:.1f}")
        print(f"  SS_sigmoid={ss_sig3:.4f}, SS_linear={ss_lin3:.4f}")
        print(f"  Sigmoid better: {ss_sig3 < ss_lin3}, F={F3:.2f}")
    except Exception as e:
        print(f"  DM-Fluence sigmoid failed: {e}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("SUMMARY")  
print("=" * 60)
print(f"""
DATA: {len(chime)} CHIME FRBs + {len(loc_frbs)} localized FRBs
ALL REAL CATALOG DATA — no simulation.

CLOSURE PREDICTIONS:
  C1: DM-Scattering power-law changes at high DM
  C2: DM-Width coupling evolves with DM
  C3: DM-SpectralIndex correlation EMERGES at high DM
  C4: Observable rank compression at high DM
  C5: DM-RM decoupling with redshift
  C6: Sigmoid transition in rolling correlations

NULL PREDICTIONS:
  × Constant τ ∝ DM^α everywhere
  × Width, fluence, sp_idx independent of DM at all ranges
  × DM-RM coupling constant with z
  × No threshold behavior
""")

results = {
    'n_chime': len(chime),
    'n_localized': len(loc_frbs),
    'source': 'CHIME Cat1 + FRBs/FRB repo',
    'real_data': True
}
with open('/root/closure-bayesn/results_frb_v3.json', 'w') as f:
    json.dump(results, f, indent=2)
print("  Results saved to results_frb_v3.json")
