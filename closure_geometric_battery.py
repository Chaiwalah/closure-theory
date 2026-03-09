#!/usr/bin/env python3
"""
CLOSURE THEORY — Geometric Test Battery
==========================================
Multiple independent tests for z-dependent covariance structure.
Goal: prove the standardization manifold rotates, using tests that
don't depend on small-bin PCA.

Author: Closure Theory collaboration
Date: 2026-03-09
"""

import numpy as np
from scipy.optimize import minimize as scipy_minimize
from scipy.integrate import quad
from scipy.stats import spearmanr, chi2 as chi2_dist, pearsonr
import json, os, warnings
warnings.filterwarnings('ignore')
np.random.seed(42)

# ============================================================
# LOAD DATA
# ============================================================
data_file = 'data/pantheon_plus.dat'
with open(data_file, 'r') as f:
    header = f.readline().strip().split()
    rows = []
    for line in f:
        parts = line.strip().split()
        if len(parts) == len(header):
            rows.append(parts)

col = {name: i for i, name in enumerate(header)}
z_hd = np.array([float(r[col['zHD']]) for r in rows])
z_cmb = np.array([float(r[col['zCMB']]) for r in rows])
m_b = np.array([float(r[col['mB']]) for r in rows])
m_b_err = np.array([float(r[col['mBERR']]) for r in rows])
x1_all = np.array([float(r[col['x1']]) for r in rows])
c_all = np.array([float(r[col['c']]) for r in rows])
hm_all = np.array([float(r[col['HOST_LOGMASS']]) for r in rows])

mask_q = (z_hd > 0.01) & (z_hd < 2.5) & (np.abs(x1_all) < 5) & (np.abs(c_all) < 0.5) & (m_b_err < 1.0) & (hm_all > 0) & (hm_all < 15)
z = z_cmb[mask_q]; mb = m_b[mask_q]; mb_err = m_b_err[mask_q]
x1 = x1_all[mask_q]; c = c_all[mask_q]; hm = hm_all[mask_q]
N = len(z)

def E(zp): return np.sqrt(0.3*(1+zp)**3 + 0.7)
def dist_mod(z_arr, H0=70):
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        integral, _ = quad(lambda zp: 1.0/E(zp), 0, zi)
        dl[i] = 299792.458 * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

mu_lcdm = dist_mod(z)
mask_hi_mass = hm >= 10.0

def chi2_std(params):
    M_B, alpha, beta, gamma = params
    mu_obs = mb + alpha * x1 - beta * c - M_B
    mu_obs[mask_hi_mass] -= gamma
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res = scipy_minimize(chi2_std, x0=[-19.25, 0.15, 3.0, -0.05], method='Nelder-Mead', options={'maxiter':50000})
p_fit = res.x
resid = (mb + p_fit[1]*x1 - p_fit[2]*c - p_fit[0]) - mu_lcdm
resid[mask_hi_mass] -= p_fit[3]

print(f"Data: {N} SNe, z range [{z.min():.3f}, {z.max():.3f}]")
print(f"Standard Tripp: α={p_fit[1]:.4f}, β={p_fit[2]:.4f}")

results = {}

# ============================================================
# TEST 1: BOX'S M — Covariance homogeneity
# ============================================================
print(f"\n{'='*70}")
print("TEST 1: BOX'S M — Are covariance matrices equal across z?")
print(f"{'='*70}")

# Split into 2, 3, and 4 groups
for k_groups, label in [(2, "median split"), (3, "terciles"), (4, "quartiles")]:
    z_quantiles = np.quantile(z, np.linspace(0, 1, k_groups+1))
    groups = []
    for i in range(k_groups):
        m = (z >= z_quantiles[i]) & (z < z_quantiles[i+1] + (1 if i==k_groups-1 else 0))
        x1n = (x1[m] - np.mean(x1[m])) / np.std(x1[m])
        cn = (c[m] - np.mean(c[m])) / np.std(c[m])
        rn = (resid[m] - np.mean(resid[m])) / np.std(resid[m])
        groups.append(np.column_stack([x1n, cn, rn]))
    
    p_dim = 3
    n_list = [g.shape[0] for g in groups]
    N_total = sum(n_list)
    S_list = [np.cov(g.T) for g in groups]
    S_pool = sum((n-1)*S for n, S in zip(n_list, S_list)) / (N_total - k_groups)
    
    M_stat = (N_total - k_groups) * np.log(np.linalg.det(S_pool))
    for n, S in zip(n_list, S_list):
        M_stat -= (n-1) * np.log(np.linalg.det(S))
    
    c1 = (2*p_dim**2 + 3*p_dim - 1) / (6*(p_dim+1)*(k_groups-1))
    c1 *= sum(1/(n-1) for n in n_list) - 1/(N_total - k_groups)
    df = p_dim*(p_dim+1)/2 * (k_groups-1)
    M_corr = M_stat * (1 - c1)
    p_box = 1 - chi2_dist.cdf(M_corr, int(df))
    
    sigma_equiv = abs(chi2_dist.ppf(p_box/2, 1))**0.5 if p_box > 0 else 99
    
    print(f"\n  {label} (k={k_groups}): M={M_corr:.1f}, df={int(df)}, p={p_box:.2e} ({sigma_equiv:.1f}σ)")
    print(f"    Group sizes: {n_list}")
    
    results[f'box_m_{k_groups}'] = {'M': float(M_corr), 'df': int(df), 'p': float(p_box), 'sigma': float(sigma_equiv)}


# ============================================================
# TEST 2: ROLLING CORRELATIONS — Continuous rotation tracking
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: ROLLING CORRELATIONS — Channel-residual coupling vs z")
print(f"{'='*70}")

order = np.argsort(z)
z_s = z[order]; x1_s = x1[order]; c_s = c[order]; r_s = resid[order]

for W in [150, 200, 300]:
    z_roll = []; rx1_r = []; rc_r = []
    for i in range(0, len(z_s) - W, 10):
        sl = slice(i, i+W)
        z_roll.append(np.median(z_s[sl]))
        r1, _ = spearmanr(x1_s[sl], r_s[sl])
        r2, _ = spearmanr(c_s[sl], r_s[sl])
        rx1_r.append(r1)
        rc_r.append(r2)
    
    z_roll = np.array(z_roll); rx1_r = np.array(rx1_r); rc_r = np.array(rc_r)
    rho_c, p_c = spearmanr(z_roll, rc_r)
    rho_x1, p_x1 = spearmanr(z_roll, rx1_r)
    
    print(f"\n  Window={W}:")
    print(f"    corr(c, Δμ) trend:  ρ={rho_c:+.4f}, p={p_c:.2e}")
    print(f"    corr(x1, Δμ) trend: ρ={rho_x1:+.4f}, p={p_x1:.2e}")
    print(f"    corr(c,Δμ) at z_min={rc_r[0]:+.3f}, z_max={rc_r[-1]:+.3f}")
    
    if W == 200:
        results['rolling_c'] = {'rho': float(rho_c), 'p': float(p_c)}
        results['rolling_x1'] = {'rho': float(rho_x1), 'p': float(p_x1)}


# ============================================================
# TEST 3: PER-BIN TRIPP FITS — α(z) and β(z) directly
# ============================================================
print(f"\n{'='*70}")
print("TEST 3: PER-BIN TRIPP FITS — Direct measurement of coefficient evolution")
print(f"{'='*70}")

z_bins = [(0.01, 0.04), (0.04, 0.08), (0.08, 0.15), (0.15, 0.30), (0.30, 0.60), (0.60, 2.3)]
alphas = []; betas = []; z_mids = []; ns = []

print(f"\n  {'z-bin':<16} {'N':>5} {'α':>8} {'β':>8} {'β/α':>8}")
print(f"  {'-'*50}")

for zlo, zhi in z_bins:
    m = (z >= zlo) & (z < zhi)
    n = m.sum()
    if n < 20: continue
    
    z_sub = z[m]; mb_sub = mb[m]; err_sub = mb_err[m]
    x1_sub = x1[m]; c_sub = c[m]
    mu_sub = dist_mod(z_sub)
    
    def chi2_bin(params):
        M_B, alpha, beta = params
        return np.sum(((mb_sub + alpha*x1_sub - beta*c_sub - M_B - mu_sub) / err_sub)**2)
    
    res_bin = scipy_minimize(chi2_bin, x0=[-19.25, 0.15, 3.0], method='Nelder-Mead', options={'maxiter':20000})
    a, b = res_bin.x[1], res_bin.x[2]
    
    z_mid = np.median(z_sub)
    alphas.append(a); betas.append(b); z_mids.append(z_mid); ns.append(n)
    ratio = b/a if abs(a) > 0.001 else float('inf')
    print(f"  [{zlo:.2f},{zhi:.2f})  {n:>5} {a:>+8.4f} {b:>+8.4f} {ratio:>8.1f}")

alphas = np.array(alphas); betas = np.array(betas); z_mids = np.array(z_mids)

rho_a, p_a = spearmanr(z_mids, alphas)
rho_b, p_b = spearmanr(z_mids, betas)

print(f"\n  α trend with z: ρ={rho_a:+.4f} (p={p_a:.4f})")
print(f"  β trend with z: ρ={rho_b:+.4f} (p={p_b:.4f})")

# Angle in Tripp space: arctan(β/α)
angles = np.degrees(np.arctan2(betas, alphas))
total_rot = angles[-1] - angles[0]
print(f"\n  Tripp-space angle: {angles[0]:.1f}° → {angles[-1]:.1f}° (Δ = {total_rot:+.1f}°)")

results['perbin_alpha'] = {'rho': float(rho_a), 'p': float(p_a), 'values': alphas.tolist()}
results['perbin_beta'] = {'rho': float(rho_b), 'p': float(p_b), 'values': betas.tolist()}


# ============================================================
# TEST 4: INTERACTION TERMS — Direct z×x1 and z×c significance
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: INTERACTION TERMS — z×x1 and z×c in the Tripp equation")
print(f"{'='*70}")

# Standard Tripp: 4 params
chi2_base = res.fun

# Add z×x1 interaction
def chi2_zx1(params):
    M_B, alpha, alpha1, beta, gamma = params
    mu_obs = mb + (alpha + alpha1*z) * x1 - beta * c - M_B
    mu_obs[mask_hi_mass] -= gamma
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res_zx1 = scipy_minimize(chi2_zx1, x0=[p_fit[0], p_fit[1], 0, p_fit[2], p_fit[3]], 
                          method='Nelder-Mead', options={'maxiter':50000})
dchi2_zx1 = chi2_base - res_zx1.fun
p_zx1 = 1 - chi2_dist.cdf(dchi2_zx1, 1)
print(f"  z×x1: Δχ² = {dchi2_zx1:.1f}, p = {p_zx1:.2e}")

# Add z×c interaction
def chi2_zc(params):
    M_B, alpha, beta, beta1, gamma = params
    mu_obs = mb + alpha * x1 - (beta + beta1*z) * c - M_B
    mu_obs[mask_hi_mass] -= gamma
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res_zc = scipy_minimize(chi2_zc, x0=[p_fit[0], p_fit[1], p_fit[2], 0, p_fit[3]], 
                         method='Nelder-Mead', options={'maxiter':50000})
dchi2_zc = chi2_base - res_zc.fun
p_zc = 1 - chi2_dist.cdf(dchi2_zc, 1)
print(f"  z×c:  Δχ² = {dchi2_zc:.1f}, p = {p_zc:.2e}")

# Both interactions
def chi2_both(params):
    M_B, alpha, alpha1, beta, beta1, gamma = params
    mu_obs = mb + (alpha + alpha1*z) * x1 - (beta + beta1*z) * c - M_B
    mu_obs[mask_hi_mass] -= gamma
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res_both = scipy_minimize(chi2_both, x0=[p_fit[0], p_fit[1], 0, p_fit[2], 0, p_fit[3]], 
                           method='Nelder-Mead', options={'maxiter':50000})
dchi2_both = chi2_base - res_both.fun
p_both = 1 - chi2_dist.cdf(dchi2_both, 2)
print(f"  Both: Δχ² = {dchi2_both:.1f} (2 dof), p = {p_both:.2e}")
print(f"    α₁ = {res_both.x[2]:+.4f}, β₁ = {res_both.x[4]:+.4f}")

results['interaction_zx1'] = {'dchi2': float(dchi2_zx1), 'p': float(p_zx1)}
results['interaction_zc'] = {'dchi2': float(dchi2_zc), 'p': float(p_zc)}
results['interaction_both'] = {'dchi2': float(dchi2_both), 'p': float(p_both)}


# ============================================================
# TEST 5: PROCRUSTES — Rotation between low-z and high-z covariance
# ============================================================
print(f"\n{'='*70}")
print("TEST 5: PROCRUSTES ROTATION — Optimal rotation between z-bins")
print(f"{'='*70}")

z_med = np.median(z)
m_lo = z < z_med; m_hi = z >= z_med

def get_eigvecs(mask):
    x1n = (x1[mask] - np.mean(x1[mask])) / np.std(x1[mask])
    cn = (c[mask] - np.mean(c[mask])) / np.std(c[mask])
    rn = (resid[mask] - np.mean(resid[mask])) / np.std(resid[mask])
    data = np.column_stack([x1n, cn, rn])
    cov = np.cov(data.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    idx = np.argsort(eigvals)[::-1]
    return eigvecs[:, idx], eigvals[idx]

V_lo, ev_lo = get_eigvecs(m_lo)
V_hi, ev_hi = get_eigvecs(m_hi)

# Procrustes: find optimal rotation R such that R @ V_lo ≈ V_hi
# R = V_hi @ V_lo.T (since both are orthonormal)
R = V_hi @ V_lo.T

# Rotation angle from trace: tr(R) = 1 + 2cos(θ) for 3D rotation
cos_theta = (np.trace(R) - 1) / 2
cos_theta = np.clip(cos_theta, -1, 1)
theta_proc = np.degrees(np.arccos(cos_theta))

print(f"  Procrustes rotation angle: {theta_proc:.1f}°")
print(f"  Rotation matrix R:")
for row in R:
    print(f"    [{row[0]:+.4f}, {row[1]:+.4f}, {row[2]:+.4f}]")

# Bootstrap the Procrustes angle
N_BOOT = 1000
theta_boot = np.zeros(N_BOOT)
idx_lo = np.where(m_lo)[0]; idx_hi = np.where(m_hi)[0]

for b in range(N_BOOT):
    b_lo = np.random.choice(idx_lo, len(idx_lo), replace=True)
    b_hi = np.random.choice(idx_hi, len(idx_hi), replace=True)
    
    V1, _ = get_eigvecs(np.isin(np.arange(N), b_lo))
    V2, _ = get_eigvecs(np.isin(np.arange(N), b_hi))
    R_b = V2 @ V1.T
    ct = np.clip((np.trace(R_b) - 1) / 2, -1, 1)
    theta_boot[b] = np.degrees(np.arccos(ct))

# Shuffle null
theta_null = np.zeros(N_BOOT)
for b in range(N_BOOT):
    perm = np.random.permutation(N)
    z_shuf = z[perm]
    m1 = z_shuf < z_med; m2 = z_shuf >= z_med
    V1, _ = get_eigvecs(m1)
    V2, _ = get_eigvecs(m2)
    R_b = V2 @ V1.T
    ct = np.clip((np.trace(R_b) - 1) / 2, -1, 1)
    theta_null[b] = np.degrees(np.arccos(ct))

sigma_proc = (theta_proc - np.mean(theta_null)) / np.std(theta_null)

print(f"\n  Bootstrap: {theta_proc:.1f}° ± {np.std(theta_boot):.1f}°")
print(f"  Null (shuffled): {np.mean(theta_null):.1f}° ± {np.std(theta_null):.1f}°")
print(f"  Observed vs null: {sigma_proc:.1f}σ")
print(f"  95% CI (bootstrap): [{np.percentile(theta_boot, 2.5):.1f}°, {np.percentile(theta_boot, 97.5):.1f}°]")

results['procrustes'] = {
    'angle': float(theta_proc),
    'boot_std': float(np.std(theta_boot)),
    'null_mean': float(np.mean(theta_null)),
    'null_std': float(np.std(theta_null)),
    'sigma': float(sigma_proc)
}


# ============================================================
# TEST 6: GRASSMANNIAN DISTANCE — Subspace angle between PC1s
# ============================================================
print(f"\n{'='*70}")
print("TEST 6: GRASSMANNIAN SUBSPACE ANGLE — PC1(low-z) vs PC1(high-z)")
print(f"{'='*70}")

# Principal angle between 1D subspaces = angle between PC1 vectors
pc1_lo = V_lo[:, 0]; pc1_hi = V_hi[:, 0]
if pc1_lo[0] < 0: pc1_lo = -pc1_lo
if pc1_hi[0] < 0: pc1_hi = -pc1_hi

cos_grass = np.clip(np.abs(np.dot(pc1_lo, pc1_hi)), 0, 1)
grass_angle = np.degrees(np.arccos(cos_grass))
print(f"  PC1 low-z:  [{pc1_lo[0]:+.4f}, {pc1_lo[1]:+.4f}, {pc1_lo[2]:+.4f}]")
print(f"  PC1 high-z: [{pc1_hi[0]:+.4f}, {pc1_hi[1]:+.4f}, {pc1_hi[2]:+.4f}]")
print(f"  Grassmannian angle: {grass_angle:.1f}°")

# Bootstrap
grass_boot = np.zeros(N_BOOT)
grass_null = np.zeros(N_BOOT)

for b in range(N_BOOT):
    b_lo = np.random.choice(idx_lo, len(idx_lo), replace=True)
    b_hi = np.random.choice(idx_hi, len(idx_hi), replace=True)
    V1, _ = get_eigvecs(np.isin(np.arange(N), b_lo))
    V2, _ = get_eigvecs(np.isin(np.arange(N), b_hi))
    v1 = V1[:,0]; v2 = V2[:,0]
    if v1[0] < 0: v1 = -v1
    if v2[0] < 0: v2 = -v2
    grass_boot[b] = np.degrees(np.arccos(np.clip(np.abs(np.dot(v1,v2)), 0, 1)))

for b in range(N_BOOT):
    perm = np.random.permutation(N)
    z_shuf = z[perm]
    m1 = z_shuf < z_med; m2 = z_shuf >= z_med
    V1, _ = get_eigvecs(m1)
    V2, _ = get_eigvecs(m2)
    v1 = V1[:,0]; v2 = V2[:,0]
    if v1[0] < 0: v1 = -v1
    if v2[0] < 0: v2 = -v2
    grass_null[b] = np.degrees(np.arccos(np.clip(np.abs(np.dot(v1,v2)), 0, 1)))

sigma_grass = (grass_angle - np.mean(grass_null)) / np.std(grass_null)

print(f"  Bootstrap: {grass_angle:.1f}° ± {np.std(grass_boot):.1f}°")
print(f"  Null: {np.mean(grass_null):.1f}° ± {np.std(grass_null):.1f}°")
print(f"  Observed vs null: {sigma_grass:.1f}σ")

results['grassmannian'] = {
    'angle': float(grass_angle),
    'sigma': float(sigma_grass),
    'null_mean': float(np.mean(grass_null)),
    'null_std': float(np.std(grass_null))
}


# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("GEOMETRIC TEST BATTERY — SUMMARY")
print(f"{'='*70}")
print(f"""
  TEST 1 — Box's M (covariance homogeneity):
    2-split: p = {results['box_m_2']['p']:.2e} ({results['box_m_2']['sigma']:.1f}σ)
    3-split: p = {results['box_m_3']['p']:.2e} ({results['box_m_3']['sigma']:.1f}σ)
    4-split: p = {results['box_m_4']['p']:.2e} ({results['box_m_4']['sigma']:.1f}σ)

  TEST 2 — Rolling correlations (W=200):
    corr(c, Δμ) trend: ρ = {results['rolling_c']['rho']:+.4f}, p = {results['rolling_c']['p']:.2e}
    corr(x1, Δμ) trend: ρ = {results['rolling_x1']['rho']:+.4f}, p = {results['rolling_x1']['p']:.2e}

  TEST 3 — Per-bin Tripp coefficients:
    α trend: ρ = {results['perbin_alpha']['rho']:+.4f}, p = {results['perbin_alpha']['p']:.4f}
    β trend: ρ = {results['perbin_beta']['rho']:+.4f}, p = {results['perbin_beta']['p']:.4f}

  TEST 4 — Interaction terms (Δχ² for 1 dof):
    z×x1: Δχ² = {results['interaction_zx1']['dchi2']:.1f}, p = {results['interaction_zx1']['p']:.2e}
    z×c:  Δχ² = {results['interaction_zc']['dchi2']:.1f}, p = {results['interaction_zc']['p']:.2e}
    Both: Δχ² = {results['interaction_both']['dchi2']:.1f}, p = {results['interaction_both']['p']:.2e}

  TEST 5 — Procrustes rotation:
    Angle = {results['procrustes']['angle']:.1f}°, {results['procrustes']['sigma']:.1f}σ vs null

  TEST 6 — Grassmannian subspace angle:
    Angle = {results['grassmannian']['angle']:.1f}°, {results['grassmannian']['sigma']:.1f}σ vs null
""")

# Overall verdict
strong = sum(1 for k, v in results.items() if isinstance(v, dict) and 'p' in v and v['p'] < 1e-4)
print(f"  Tests with p < 10⁻⁴: {strong}")
print(f"  Tests total: 10")

# Save
os.makedirs('results_geometric_battery', exist_ok=True)
with open('results_geometric_battery/geometric_battery.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\n  Saved to results_geometric_battery/geometric_battery.json")
