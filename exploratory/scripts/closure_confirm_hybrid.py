#!/usr/bin/env python3
"""
closure_confirm_hybrid.py — Confirm or Deny Grok/Gemini Claims
================================================================

Grok claims: φ(1,Z) = αZ — color trace-negativity is LINEAR in metallicity
Gemini claims: enriched medium acts as resonator/modulator, not just source-side

Three decisive tests:

TEST 1: SLOW CHANNEL HOST-MASS SENSITIVITY
  - If slow channel is ALSO affected by host mass → medium contribution (Gemini/Grok)
  - If slow channel is IMMUNE → pure source-side (GPT)
  This is the clean separator. The slow channel has no absorbing boundary.
  If the medium is doing work, it should show up there too.

TEST 2: COLOR ÉTENDUE vs HOST MASS — IS IT LINEAR IN Z?
  - Grok says φ(1,Z) = αZ (linear). Test: bin extended-fast by host mass
  - Compute ℰ(c) drop rate as function of host mass
  - If linear → Grok confirmed. If threshold → different mechanism.

TEST 3: CROSS-CHANNEL COHERENCE
  - If medium modulates, BOTH corridors should show correlated z-trends
    in SOME observable (even if magnitudes differ)
  - Compute rank correlation of z-trends BETWEEN core and extended
  - High coherence → shared driver (medium). Low → independent (source).

Author: Closure Theory Collaboration
Date: 2026-03-06
"""

import numpy as np
from scipy.stats import spearmanr, gaussian_kde
from scipy.optimize import minimize
from scipy.stats import norm
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

sn = np.loadtxt('data/pantheon_plus.dat', dtype=str)
header = list(sn[0])
data = sn[1:]
col = {h: i for i, h in enumerate(header)}

z = data[:, col['zHD']].astype(float)
mb = data[:, col['mB']].astype(float)
x1 = data[:, col['x1']].astype(float)
c = data[:, col['c']].astype(float)
host_mass = data[:, col['HOST_LOGMASS']].astype(float)

mask = (z > 0.01) & (z < 2.5) & (np.abs(c) < 0.3) & (np.abs(x1) < 3)
z, mb, x1, c, host_mass = z[mask], mb[mask], x1[mask], c[mask], host_mass[mask]
print(f"N = {len(z)} SNe\n")

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 2.5]

def collision_prob(x):
    if len(x) < 12: return np.nan
    try:
        kde = gaussian_kde(x)
        grid = np.linspace(np.percentile(x, 1), np.percentile(x, 99), 200)
        f_vals = kde.evaluate(grid)
        dx = grid[1] - grid[0]
        return np.sum(f_vals**2) * dx
    except: return np.nan

# ============================================================
# TEST 1: SLOW CHANNEL — HOST MASS SENSITIVITY
# ============================================================
print("=" * 110)
print("TEST 1: SLOW CHANNEL (x1≥0) — Is it affected by host mass?")
print("  If YES → medium/environment contributes (Gemini/Grok win)")
print("  If NO  → pure source-side (GPT wins)")
print("=" * 110)

slow_mask = x1 >= 0
mass_valid = host_mass > 0
mass_med_slow = np.median(host_mass[mass_valid & slow_mask])

for obs_name, obs_data in [('x1', x1), ('c', c), ('mB', mb)]:
    print(f"\n  {obs_name}:")
    for mass_label, mass_cut in [('ALL', np.ones(len(z), bool)),
                                  ('LOW-mass', (host_mass < mass_med_slow) & mass_valid),
                                  ('HIGH-mass', (host_mass >= mass_med_slow) & mass_valid)]:
        z_l, et_l, var_l = [], [], []
        for i in range(len(z_edges) - 1):
            z_lo, z_hi = z_edges[i], z_edges[i+1]
            bm = (z >= z_lo) & (z < z_hi) & slow_mask & mass_cut
            n = np.sum(bm)
            if n < 15: continue
            v = np.var(obs_data[bm])
            cp = collision_prob(obs_data[bm])
            if not np.isnan(cp):
                z_l.append(np.mean(z[bm]))
                et_l.append(v * cp)
                var_l.append(v)
        
        if len(z_l) >= 3:
            rho_e, p_e = spearmanr(z_l, et_l)
            rho_v, p_v = spearmanr(z_l, var_l)
            print(f"    {mass_label:>12}: ℰ vs z: ρ = {rho_e:+.3f} p={p_e:.3f} | Var vs z: ρ = {rho_v:+.3f} p={p_v:.3f}")

# ============================================================
# TEST 2: COLOR ÉTENDUE DROP RATE vs HOST MASS (EXTENDED-FAST)
# ============================================================
print(f"\n\n{'=' * 110}")
print("TEST 2: Is color étendue loss LINEAR in host mass? (Grok: φ(1,Z)=αZ)")
print("=" * 110)

fast_mask = x1 < 0

def fit_2gauss(x):
    best = None
    best_ll = 1e10
    for w0 in [0.3, 0.5, 0.7]:
        for m1_0 in [-0.3, -0.5]:
            for m2_0 in [-1.0, -1.5]:
                try:
                    def neg_loglik(params):
                        w, m1, s1, m2, s2 = params
                        if w < 0.01 or w > 0.99 or s1 < 0.01 or s2 < 0.01: return 1e10
                        p = w * norm.pdf(x, m1, s1) + (1-w) * norm.pdf(x, m2, s2)
                        return -np.sum(np.log(np.maximum(p, 1e-30)))
                    res = minimize(neg_loglik, [w0, m1_0, 0.3, m2_0, 0.5],
                                   method='Nelder-Mead', options={'maxiter': 5000})
                    if res.fun < best_ll:
                        best_ll = res.fun
                        best = res.x
                except: pass
    if best is None: return None
    w, m1, s1, m2, s2 = best
    s1, s2 = abs(s1), abs(s2)
    if m1 < m2: w, m1, s1, m2, s2 = 1-w, m2, s2, m1, s1
    return w, m1, s1, m2, s2

# Split into mass terciles for finer resolution
mass_terciles = np.percentile(host_mass[mass_valid & fast_mask], [33.3, 66.7])
mass_bins = [
    ('LOW (T1)', (host_mass < mass_terciles[0]) & mass_valid),
    ('MID (T2)', (host_mass >= mass_terciles[0]) & (host_mass < mass_terciles[1]) & mass_valid),
    ('HIGH (T3)', (host_mass >= mass_terciles[1]) & mass_valid)
]

print(f"\n  Mass terciles: T1 < {mass_terciles[0]:.2f}, T2 = [{mass_terciles[0]:.2f}, {mass_terciles[1]:.2f}), T3 ≥ {mass_terciles[1]:.2f}")
print(f"\n  {'Mass bin':>12} | {'ℰ(c)_ext slope':>16} | {'ℰ(x1)_ext slope':>17} | {'σ_ext slope':>13}")
print(f"  " + "-" * 70)

mass_centers = []
color_slopes = []

for label, mcut in mass_bins:
    z_l = []
    et_c_ext, et_x1_ext, sig_ext_l = [], [], []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bm = (z >= z_lo) & (z < z_hi) & fast_mask & mcut
        n = np.sum(bm)
        if n < 20: continue
        
        fit = fit_2gauss(x1[bm])
        if fit is None: continue
        w_core, mu_core, sig_core, mu_ext, sig_ext = fit
        
        # Assign to extended
        prob_core = w_core * norm.pdf(x1[bm], mu_core, sig_core)
        prob_ext = (1-w_core) * norm.pdf(x1[bm], mu_ext, sig_ext)
        is_ext = (prob_ext / (prob_core + prob_ext)) > 0.5
        
        if np.sum(is_ext) < 10: continue
        
        zc = np.mean(z[bm])
        z_l.append(zc)
        
        # Color étendue of extended
        c_ext = c[bm][is_ext]
        v_c = np.var(c_ext)
        cp_c = collision_prob(c_ext)
        et_c_ext.append(v_c * cp_c if not np.isnan(cp_c) else np.nan)
        
        # x1 étendue of extended
        x1_ext = x1[bm][is_ext]
        v_x = np.var(x1_ext)
        cp_x = collision_prob(x1_ext)
        et_x1_ext.append(v_x * cp_x if not np.isnan(cp_x) else np.nan)
        
        sig_ext_l.append(sig_ext)
    
    if len(z_l) >= 3:
        valid_c = [i for i in range(len(z_l)) if not np.isnan(et_c_ext[i])]
        valid_x = [i for i in range(len(z_l)) if not np.isnan(et_x1_ext[i])]
        
        rho_c = spearmanr([z_l[i] for i in valid_c], [et_c_ext[i] for i in valid_c])[0] if len(valid_c) >= 3 else np.nan
        rho_x = spearmanr([z_l[i] for i in valid_x], [et_x1_ext[i] for i in valid_x])[0] if len(valid_x) >= 3 else np.nan
        rho_s = spearmanr(z_l, sig_ext_l)[0] if len(z_l) >= 3 else np.nan
        
        print(f"  {label:>12} | {rho_c:>+16.3f} | {rho_x:>+17.3f} | {rho_s:>+13.3f}")
        
        if not np.isnan(rho_c):
            mc = np.mean(host_mass[fast_mask & mcut & mass_valid])
            mass_centers.append(mc)
            color_slopes.append(rho_c)

if len(mass_centers) >= 3:
    rho_linear, p_linear = spearmanr(mass_centers, color_slopes)
    print(f"\n  Color étendue drop rate vs host mass: ρ = {rho_linear:+.3f}, p = {p_linear:.3f}")
    if abs(rho_linear) > 0.8:
        print(f"  → {'LINEAR (Grok confirmed)' if rho_linear < -0.5 else 'ANTI-LINEAR'}")
    else:
        print(f"  → Not clearly linear — threshold or nonlinear relationship")


# ============================================================
# TEST 3: CROSS-CORRIDOR COHERENCE
# ============================================================
print(f"\n\n{'=' * 110}")
print("TEST 3: Cross-Corridor Coherence — Shared Driver or Independent?")
print("  High coherence → shared external driver (medium)")
print("  Low coherence → independent internal processes (source)")
print("=" * 110)

z_list3 = []
var_core_x1, var_ext_x1 = [], []
var_core_c, var_ext_c = [], []
var_core_mb, var_ext_mb = [], []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bm = (z >= z_lo) & (z < z_hi) & fast_mask
    n = np.sum(bm)
    if n < 30: continue
    
    fit = fit_2gauss(x1[bm])
    if fit is None: continue
    w_core, mu_core, sig_core, mu_ext, sig_ext = fit
    
    prob_core = w_core * norm.pdf(x1[bm], mu_core, sig_core)
    prob_ext = (1-w_core) * norm.pdf(x1[bm], mu_ext, sig_ext)
    is_core = (prob_core / (prob_core + prob_ext)) > 0.5
    is_ext = ~is_core
    
    if np.sum(is_core) < 10 or np.sum(is_ext) < 10: continue
    
    zc = np.mean(z[bm])
    z_list3.append(zc)
    
    var_core_x1.append(np.var(x1[bm][is_core]))
    var_ext_x1.append(np.var(x1[bm][is_ext]))
    var_core_c.append(np.var(c[bm][is_core]))
    var_ext_c.append(np.var(c[bm][is_ext]))
    var_core_mb.append(np.var(mb[bm][is_core]))
    var_ext_mb.append(np.var(mb[bm][is_ext]))

if len(z_list3) >= 3:
    print(f"\n  Cross-corridor Var(z) correlation (do they move together?):")
    for name, vc, ve in [('x1', var_core_x1, var_ext_x1), 
                          ('c', var_core_c, var_ext_c),
                          ('mB', var_core_mb, var_ext_mb)]:
        rho, p = spearmanr(vc, ve)
        print(f"    Var({name})_core vs Var({name})_ext: ρ = {rho:+.3f}, p = {p:.3f} {'← COHERENT' if rho > 0.6 else '← INDEPENDENT' if rho < 0.3 else ''}")
    
    # Also: do their z-TRENDS correlate?
    # Detrend each against z, then correlate residuals
    from numpy.polynomial.polynomial import polyfit, polyval
    print(f"\n  Z-detrended residual correlation (shared hidden driver?):")
    for name, vc, ve in [('x1', var_core_x1, var_ext_x1),
                          ('c', var_core_c, var_ext_c)]:
        coef_c = polyfit(z_list3, vc, 1)
        coef_e = polyfit(z_list3, ve, 1)
        res_c = np.array(vc) - polyval(z_list3, coef_c)
        res_e = np.array(ve) - polyval(z_list3, coef_e)
        rho, p = spearmanr(res_c, res_e)
        print(f"    Var({name}) residuals: ρ = {rho:+.3f}, p = {p:.3f} {'← SHARED DRIVER' if rho > 0.5 else '← INDEPENDENT' if abs(rho) < 0.3 else ''}")


# ============================================================
# VERDICT
# ============================================================
print(f"\n\n{'=' * 110}")
print("VERDICT: Source-Side (GPT) vs Hybrid (Grok) vs Medium (Gemini)")
print("=" * 110)
print("""
Decision tree:
  Slow channel affected by host mass? 
    YES → Medium contributes (Gemini/Grok)
    NO  → Pure source-side (GPT)
  
  Color drop linear in Z?
    YES → Grok's φ(1,Z)=αZ confirmed
    NO  → Threshold or nonlinear
  
  Cross-corridor coherence?
    HIGH → Shared external driver (medium)
    LOW  → Independent source processes
""")
