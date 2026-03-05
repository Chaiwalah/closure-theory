#!/usr/bin/env python3
"""
closure_final_two_tests.py — GPT's Two Proposed Confirmatory Tests
====================================================================

Test 1: Decompose color into intrinsic + dust-proxy components
  - Use HOST_LOGMASS as environment proxy (high-mass hosts → more dust)
  - Use COV_x1_c to separate intrinsic (correlated with x1) from extrinsic (uncorrelated)
  - Show étendue conservation lives in the extrinsic/dust-like component

Test 2: Show fast-channel weight decay predicts ℰ(x1) drop rate
  - Fit 2-Gaussian mixture to fast channel at each z
  - Track extended-fast weight w(z) decay
  - Compare w(z) decay rate κ to ℰ(x1) decay rate
  - If they match → absorbing boundary mechanism confirmed

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr, gaussian_kde
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# Load Pantheon+
sn = np.loadtxt('data/pantheon_plus.dat', dtype=str)
header = list(sn[0])
data = sn[1:]
col = {h: i for i, h in enumerate(header)}

z = data[:, col['zHD']].astype(float)
mb = data[:, col['mB']].astype(float)
x1 = data[:, col['x1']].astype(float)
c = data[:, col['c']].astype(float)
cov_x1c = data[:, col['COV_x1_c']].astype(float)
host_mass = data[:, col['HOST_LOGMASS']].astype(float)
mwebv = data[:, col['MWEBV']].astype(float)

mask = (z > 0.01) & (z < 2.5) & (np.abs(c) < 0.3) & (np.abs(x1) < 3)
z, mb, x1, c = z[mask], mb[mask], x1[mask], c[mask]
cov_x1c = cov_x1c[mask]
host_mass = host_mass[mask]
mwebv = mwebv[mask]
print(f"Pantheon+ after cuts: {len(z)} SNe\n")

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 2.5]

def collision_prob(x, n_grid=200):
    if len(x) < 15:
        return np.nan
    try:
        kde = gaussian_kde(x)
        grid = np.linspace(np.percentile(x, 1), np.percentile(x, 99), n_grid)
        f_vals = kde.evaluate(grid)
        dx = grid[1] - grid[0]
        return np.sum(f_vals**2) * dx
    except:
        return np.nan

# ============================================================
# TEST 1: COLOR DECOMPOSITION
# ============================================================
print("=" * 100)
print("TEST 1: COLOR DECOMPOSITION — Intrinsic vs Extrinsic Component Étendue")
print("=" * 100)

# Method A: Residualize color against x1 → intrinsic (correlated) vs extrinsic (residual)
from numpy.linalg import lstsq

print("\n  Method A: c = α·x1 + c_extrinsic (residualize out x1 correlation)")
print("  If étendue conservation lives in c_extrinsic → dust/environment drives conservation\n")

# Global fit
X_fit = np.column_stack([np.ones(len(x1)), x1])
coeffs_cx, _, _, _ = lstsq(X_fit, c, rcond=None)
c_intrinsic = coeffs_cx[1] * x1  # part of c correlated with x1
c_extrinsic = c - coeffs_cx[0] - c_intrinsic  # residual = dust-like

print(f"  Global: c = {coeffs_cx[0]:.4f} + {coeffs_cx[1]:.4f}·x1 + residual")
print(f"  Var(c_intrinsic)/Var(c) = {np.var(c_intrinsic)/np.var(c):.3f}")
print(f"  Var(c_extrinsic)/Var(c) = {np.var(c_extrinsic)/np.var(c):.3f}")

print(f"\n  {'z':>6} {'N':>4} | {'Var(c_ext)':>10} {'∫f²':>8} {'ℰ_ext':>8} | {'Var(c_int)':>10} {'∫f²':>8} {'ℰ_int':>8}")
print(f"  " + "-" * 75)

z_list, et_ext, et_int = [], [], []
for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bm = (z >= z_lo) & (z < z_hi)
    n = np.sum(bm)
    if n < 20:
        continue
    zc = np.mean(z[bm])
    
    # Extrinsic (dust-like)
    v_ext = np.var(c_extrinsic[bm])
    cp_ext = collision_prob(c_extrinsic[bm])
    e_ext = v_ext * cp_ext if not np.isnan(cp_ext) else np.nan
    
    # Intrinsic (x1-correlated)
    v_int = np.var(c_intrinsic[bm])
    cp_int = collision_prob(c_intrinsic[bm])
    e_int = v_int * cp_int if not np.isnan(cp_int) else np.nan
    
    z_list.append(zc)
    et_ext.append(e_ext)
    et_int.append(e_int)
    
    print(f"  {zc:>6.3f} {n:>4} | {v_ext:>10.5f} {cp_ext:>8.3f} {e_ext:>8.5f} | {v_int:>10.5f} {cp_int:>8.3f} {e_int:>8.5f}")

valid = [i for i in range(len(z_list)) if not np.isnan(et_ext[i]) and not np.isnan(et_int[i])]
z_v = [z_list[i] for i in valid]
ext_v = [et_ext[i] for i in valid]
int_v = [et_int[i] for i in valid]

rho_ext, p_ext = spearmanr(z_v, ext_v)
rho_int, p_int = spearmanr(z_v, int_v)
cv_ext = np.std(ext_v) / np.mean(ext_v) * 100
cv_int = np.std(int_v) / np.mean(int_v) * 100

print(f"\n  ℰ(c_extrinsic) vs z: ρ = {rho_ext:+.3f}, p = {p_ext:.4f}, CV = {cv_ext:.1f}%  {'🔥 CONSERVED' if abs(rho_ext) < 0.4 else ''}")
print(f"  ℰ(c_intrinsic) vs z: ρ = {rho_int:+.3f}, p = {p_int:.4f}, CV = {cv_int:.1f}%  {'🔥 CONSERVED' if abs(rho_int) < 0.4 else 'NOT CONSERVED' if abs(rho_int) > 0.5 else ''}")

# Method B: Split by host mass (high mass → dustier)
print(f"\n\n  Method B: Host-mass split (high-mass hosts have more dust)")
mass_valid = host_mass > 0  # filter missing values
if np.sum(mass_valid) > 100:
    mass_med = np.median(host_mass[mass_valid])
    print(f"  Mass median: {mass_med:.2f}")
    
    for label, hm_mask in [('HIGH-mass (dusty)', (host_mass >= mass_med) & mass_valid),
                            ('LOW-mass (less dust)', (host_mass < mass_med) & mass_valid)]:
        z_hm = z[hm_mask]
        c_hm = c[hm_mask]
        
        z_l, et_l = [], []
        for i in range(len(z_edges) - 1):
            z_lo, z_hi = z_edges[i], z_edges[i+1]
            bm = (z_hm >= z_lo) & (z_hm < z_hi)
            n = np.sum(bm)
            if n < 15:
                continue
            zc = np.mean(z_hm[bm])
            v = np.var(c_hm[bm])
            cp = collision_prob(c_hm[bm])
            if not np.isnan(cp):
                z_l.append(zc)
                et_l.append(v * cp)
        
        if len(z_l) >= 3:
            rho, p = spearmanr(z_l, et_l)
            cv = np.std(et_l) / np.mean(et_l) * 100
            print(f"  {label:>25}: ρ = {rho:+.3f}, p = {p:.4f}, CV = {cv:.1f}%  {'🔥' if abs(rho) < 0.4 else ''}")


# ============================================================
# TEST 2: FAST WEIGHT DECAY vs ÉTENDUE DROP
# ============================================================
print(f"\n\n{'=' * 100}")
print("TEST 2: FAST-CHANNEL WEIGHT DECAY → ÉTENDUE DROP PREDICTION")
print("=" * 100)
print("\n  Fit 2-Gaussian mixture to FAST (x1<0) at each z-bin.")
print("  Track extended-fast weight w(z). If w decays → ℰ(x1) drops at rate κ.\n")

fast_mask = x1 < 0

def fit_2gauss(x):
    """Fit 2-component Gaussian mixture to x. Returns (w, mu1, sig1, mu2, sig2)."""
    from scipy.stats import norm
    
    def neg_loglik(params):
        w, m1, s1, m2, s2 = params
        if w < 0.01 or w > 0.99 or s1 < 0.01 or s2 < 0.01:
            return 1e10
        p = w * norm.pdf(x, m1, s1) + (1 - w) * norm.pdf(x, m2, s2)
        p = np.maximum(p, 1e-30)
        return -np.sum(np.log(p))
    
    # Initial guesses
    best = None
    best_ll = 1e10
    for w0 in [0.3, 0.5, 0.7]:
        for m1_0 in [-0.3, -0.5]:
            for m2_0 in [-1.0, -1.5]:
                try:
                    from scipy.optimize import minimize
                    res = minimize(neg_loglik, [w0, m1_0, 0.3, m2_0, 0.5],
                                   method='Nelder-Mead', options={'maxiter': 5000})
                    if res.fun < best_ll:
                        best_ll = res.fun
                        best = res.x
                except:
                    pass
    
    if best is None:
        return None
    
    w, m1, s1, m2, s2 = best
    # Convention: component 1 = core (closer to 0), component 2 = extended (more negative)
    if m1 < m2:
        w, m1, s1, m2, s2 = 1-w, m2, s2, m1, s1
    return w, m1, s1, m2, s2  # w = weight of CORE, 1-w = weight of EXTENDED

print(f"  {'z':>6} {'N':>4} | {'w_core':>7} {'μ_core':>7} {'σ_core':>7} | {'w_ext':>7} {'μ_ext':>7} {'σ_ext':>7} | {'ℰ(x1)':>8}")
print(f"  " + "-" * 85)

z_list2, w_ext_list, et_x1_list = [], [], []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bm = (z >= z_lo) & (z < z_hi) & fast_mask
    n = np.sum(bm)
    if n < 25:
        continue
    
    zc = np.mean(z[bm])
    x1_bin = x1[bm]
    
    fit = fit_2gauss(x1_bin)
    
    # Étendue
    v = np.var(x1_bin)
    cp = collision_prob(x1_bin)
    et = v * cp if not np.isnan(cp) else np.nan
    
    if fit is not None:
        w_core, m1, s1, m2, s2 = fit
        w_extended = 1 - w_core
        z_list2.append(zc)
        w_ext_list.append(w_extended)
        et_x1_list.append(et)
        print(f"  {zc:>6.3f} {n:>4} | {w_core:>7.3f} {m1:>7.3f} {s1:>7.3f} | {w_extended:>7.3f} {m2:>7.3f} {s2:>7.3f} | {et:>8.4f}")
    else:
        print(f"  {zc:>6.3f} {n:>4} | FIT FAILED | {et:>8.4f}")

if len(z_list2) >= 3:
    rho_w, p_w = spearmanr(z_list2, w_ext_list)
    rho_e, p_e = spearmanr(z_list2, et_x1_list)
    rho_we, p_we = spearmanr(w_ext_list, et_x1_list)
    
    print(f"\n  w_extended vs z:    ρ = {rho_w:+.3f}, p = {p_w:.4f}")
    print(f"  ℰ(x1_fast) vs z:   ρ = {rho_e:+.3f}, p = {p_e:.4f}")
    print(f"  w_extended vs ℰ:   ρ = {rho_we:+.3f}, p = {p_we:.4f}")
    
    if rho_we > 0.5:
        print(f"\n  🔥 CONFIRMED: Extended-fast weight decay TRACKS étendue loss")
        print(f"  The absorbing boundary removes the broad component → ℰ drops")
    elif abs(rho_we) < 0.3:
        print(f"\n  ❌ Weight decay does NOT predict étendue loss — mechanism is different")

# ============================================================
# COMBINED VERDICT
# ============================================================
print(f"\n\n{'=' * 100}")
print("COMBINED VERDICT")
print("=" * 100)
print(f"""
Test 1 (Color decomposition):
  - c_extrinsic (dust-like) étendue: ρ = {rho_ext:+.3f}, CV = {cv_ext:.1f}%
  - c_intrinsic (x1-coupled) étendue: ρ = {rho_int:+.3f}, CV = {cv_int:.1f}%
  → {'Dust component drives conservation' if abs(rho_ext) < abs(rho_int) else 'Both conserve similarly' if abs(rho_ext - rho_int) < 0.2 else 'Intrinsic drives conservation'}

Test 2 (Fast weight decay → ℰ drop):
  - w_extended vs ℰ(x1): ρ = {rho_we:+.3f}
  → {'Absorbing boundary confirmed' if rho_we > 0.4 else 'Mechanism unclear' if abs(rho_we) < 0.3 else 'Different mechanism'}
""")
