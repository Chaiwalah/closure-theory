#!/usr/bin/env python3
"""
closure_pca_predictions.py — PCA Prediction Test + Deep Metallicity Probes
===========================================================================

PART A: PCA on fast-channel covariance evolution
  - Extract eigenmodes from (mB, x1, c) covariance at each z-bin
  - Separately for core and extended corridors
  - Test all agent predictions: N modes, dominant loading, rotation angle,
    which contracts, cross-corridor eigenvalue coherence

PART B: Deep metallicity probes
  - Host mass × corridor × z interaction effects
  - Does metallicity rotate the projection or drive the contraction?
  - Eigenmode loading shifts with host mass
  - Metallicity as a hidden control variable

Author: Closure Theory Collaboration
Date: 2026-03-06
"""

import numpy as np
from scipy.stats import spearmanr
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

def assign_components(x1_arr, fit):
    w_core, mu_core, sig_core, mu_ext, sig_ext = fit
    p_core = w_core * norm.pdf(x1_arr, mu_core, sig_core)
    p_ext = (1-w_core) * norm.pdf(x1_arr, mu_ext, sig_ext)
    return p_core / (p_core + p_ext) > 0.5

def angle_between_vectors(v1, v2):
    """Angle in degrees between two vectors."""
    cos_angle = np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1)
    return np.degrees(np.arccos(abs(cos_angle)))


# ============================================================
# PART A: PCA ON FAST-CHANNEL COVARIANCE EVOLUTION
# ============================================================
print("=" * 120)
print("PART A: PCA ON FAST-CHANNEL (mB, x1, c) — Core vs Extended Corridors")
print("=" * 120)

# Standardize globally
mb_std = (mb - np.mean(mb[fast_mask])) / np.std(mb[fast_mask])
x1_std = (x1 - np.mean(x1[fast_mask])) / np.std(x1[fast_mask])
c_std = (c - np.mean(c[fast_mask])) / np.std(c[fast_mask])

z_list = []
evals_core, evals_ext = [], []
evecs_core_1, evecs_ext_1 = [], []  # First eigenvector
evecs_core_2, evecs_ext_2 = [], []  # Second eigenvector
angles_1 = []  # Angle between first eigenvectors

print(f"\n  {'z':>6} {'N_c':>4} {'N_e':>4} | {'λ1_core':>8} {'λ2_core':>8} {'λ3_core':>8} | {'λ1_ext':>8} {'λ2_ext':>8} {'λ3_ext':>8} | {'θ_PC1':>6} {'θ_PC2':>6}")
print(f"  " + "-" * 100)

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bm = (z >= z_lo) & (z < z_hi) & fast_mask
    n = np.sum(bm)
    if n < 30: continue
    
    zc = np.mean(z[bm])
    fit = fit_2gauss(x1[bm])
    if fit is None: continue
    
    is_core = assign_components(x1[bm], fit)
    is_ext = ~is_core
    n_core, n_ext = np.sum(is_core), np.sum(is_ext)
    if n_core < 10 or n_ext < 10: continue
    
    # Build data matrices
    X_core = np.column_stack([mb_std[bm][is_core], x1_std[bm][is_core], c_std[bm][is_core]])
    X_ext = np.column_stack([mb_std[bm][is_ext], x1_std[bm][is_ext], c_std[bm][is_ext]])
    
    # Covariance and eigen decomposition
    C_core = np.cov(X_core.T)
    C_ext = np.cov(X_ext.T)
    
    vals_c, vecs_c = np.linalg.eigh(C_core)
    vals_e, vecs_e = np.linalg.eigh(C_ext)
    
    # Sort descending
    idx_c = np.argsort(vals_c)[::-1]
    idx_e = np.argsort(vals_e)[::-1]
    vals_c, vecs_c = vals_c[idx_c], vecs_c[:, idx_c]
    vals_e, vecs_e = vals_e[idx_e], vecs_e[:, idx_e]
    
    z_list.append(zc)
    evals_core.append(vals_c)
    evals_ext.append(vals_e)
    evecs_core_1.append(vecs_c[:, 0])
    evecs_ext_1.append(vecs_e[:, 0])
    evecs_core_2.append(vecs_c[:, 1])
    evecs_ext_2.append(vecs_e[:, 1])
    
    # Angle between first eigenvectors
    theta1 = angle_between_vectors(vecs_c[:, 0], vecs_e[:, 0])
    theta2 = angle_between_vectors(vecs_c[:, 1], vecs_e[:, 1])
    angles_1.append(theta1)
    
    print(f"  {zc:>6.3f} {n_core:>4} {n_ext:>4} | {vals_c[0]:>8.3f} {vals_c[1]:>8.3f} {vals_c[2]:>8.3f} | {vals_e[0]:>8.3f} {vals_e[1]:>8.3f} {vals_e[2]:>8.3f} | {theta1:>6.1f}° {theta2:>6.1f}°")

# ============================================================
# PREDICTIONS SCORECARD
# ============================================================
print(f"\n\n{'=' * 120}")
print("PREDICTIONS SCORECARD")
print("=" * 120)

# 1. N significant eigenmodes
evals_all_core = np.mean(evals_core, axis=0)
evals_all_ext = np.mean(evals_ext, axis=0)
frac_core = evals_all_core / np.sum(evals_all_core)
frac_ext = evals_all_ext / np.sum(evals_all_ext)
print(f"\n  1. N SIGNIFICANT EIGENMODES")
print(f"     Core: λ fractions = [{frac_core[0]:.3f}, {frac_core[1]:.3f}, {frac_core[2]:.3f}]")
print(f"     Ext:  λ fractions = [{frac_ext[0]:.3f}, {frac_ext[1]:.3f}, {frac_ext[2]:.3f}]")
n_sig = np.sum(frac_core > 0.15)
print(f"     → {n_sig} significant modes (>15% variance)")
print(f"     ALL AGENTS PREDICTED: 2 → {'✅ CORRECT' if n_sig == 2 else '❌ WRONG'}")

# 2. Dominant eigenmode loading
print(f"\n  2. DOMINANT EIGENMODE LOADING")
labels = ['mB', 'x1', 'c']
mean_pc1_core = np.mean([np.abs(v) for v in evecs_core_1], axis=0)
mean_pc1_ext = np.mean([np.abs(v) for v in evecs_ext_1], axis=0)
dominant_core = labels[np.argmax(mean_pc1_core)]
dominant_ext = labels[np.argmax(mean_pc1_ext)]
print(f"     Core PC1 mean |loadings|: mB={mean_pc1_core[0]:.3f}, x1={mean_pc1_core[1]:.3f}, c={mean_pc1_core[2]:.3f} → {dominant_core}")
print(f"     Ext  PC1 mean |loadings|: mB={mean_pc1_ext[0]:.3f}, x1={mean_pc1_ext[1]:.3f}, c={mean_pc1_ext[2]:.3f} → {dominant_ext}")
print(f"     Grok predicted: x1 → {'✅' if dominant_core == 'x1' or dominant_ext == 'x1' else '❌'}")
print(f"     GPT predicted: mB → {'✅' if dominant_core == 'mB' or dominant_ext == 'mB' else '❌'}")
print(f"     Gemini predicted: mB → {'✅' if dominant_core == 'mB' or dominant_ext == 'mB' else '❌'}")
print(f"     Clawd predicted: mB → {'✅' if dominant_core == 'mB' or dominant_ext == 'mB' else '❌'}")

# 3. Rotation angle
mean_angle = np.mean(angles_1)
std_angle = np.std(angles_1)
print(f"\n  3. ROTATION ANGLE BETWEEN CORRIDORS (PC1)")
print(f"     Mean angle: {mean_angle:.1f}° ± {std_angle:.1f}°")
print(f"     Per-bin angles: {['%.1f°' % a for a in angles_1]}")
print(f"     Grok predicted: ~90° → {'✅' if abs(mean_angle - 90) < 20 else '❌'}")
print(f"     GPT predicted: ~30° → {'✅' if abs(mean_angle - 30) < 15 else '❌'}")
print(f"     Gemini predicted: ~90° → {'✅' if abs(mean_angle - 90) < 20 else '❌'}")
print(f"     Clawd predicted: ~20° → {'✅' if abs(mean_angle - 20) < 15 else '❌'}")

# 4. Which eigenmode contracts
print(f"\n  4. WHICH EIGENMODE CONTRACTS WITH z")
if len(z_list) >= 3:
    lam1_core = [e[0] for e in evals_core]
    lam2_core = [e[1] for e in evals_core]
    rho1, p1 = spearmanr(z_list, lam1_core)
    rho2, p2 = spearmanr(z_list, lam2_core)
    print(f"     Core: λ1 vs z: ρ={rho1:+.3f} p={p1:.3f} | λ2 vs z: ρ={rho2:+.3f} p={p2:.3f}")
    
    lam1_ext = [e[0] for e in evals_ext]
    lam2_ext = [e[1] for e in evals_ext]
    rho1e, p1e = spearmanr(z_list, lam1_ext)
    rho2e, p2e = spearmanr(z_list, lam2_ext)
    print(f"     Ext:  λ1 vs z: ρ={rho1e:+.3f} p={p1e:.3f} | λ2 vs z: ρ={rho2e:+.3f} p={p2e:.3f}")
    
    contracts_first = (rho1 < -0.3 or rho1e < -0.3)
    contracts_second = (rho2 < -0.3 or rho2e < -0.3)
    print(f"     ALL AGENTS PREDICTED: First contracts → {'✅' if contracts_first else '❌'}")

# 5. Cross-corridor λ1 coherence
print(f"\n  5. CROSS-CORRIDOR λ1 COHERENCE")
if len(z_list) >= 3:
    rho_cross, p_cross = spearmanr(lam1_core, lam1_ext)
    print(f"     λ1_core vs λ1_ext: ρ = {rho_cross:+.3f}, p = {p_cross:.3f}")
    print(f"     Grok predicted: Negative → {'✅' if rho_cross < -0.3 else '❌'}")
    print(f"     GPT predicted: Positive → {'✅' if rho_cross > 0.3 else '❌'}")
    print(f"     Gemini predicted: Positive → {'✅' if rho_cross > 0.3 else '❌'}")
    print(f"     Clawd predicted: Weak positive → {'✅' if rho_cross > 0 else '❌'}")


# ============================================================
# PART B: DEEP METALLICITY PROBES
# ============================================================
print(f"\n\n{'=' * 120}")
print("PART B: DEEP METALLICITY PROBES — How Does Host Mass Affect the Eigenmodes?")
print("=" * 120)

mass_valid = host_mass > 0
mass_med = np.median(host_mass[mass_valid & fast_mask])

for mass_label, mass_cut in [('LOW-mass', (host_mass < mass_med) & mass_valid),
                              ('HIGH-mass', (host_mass >= mass_med) & mass_valid)]:
    print(f"\n  {mass_label} hosts:")
    z_m, angles_m = [], []
    lam1_c_m, lam1_e_m = [], []
    loadings_c_m, loadings_e_m = [], []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bm = (z >= z_lo) & (z < z_hi) & fast_mask & mass_cut
        n = np.sum(bm)
        if n < 25: continue
        
        zc = np.mean(z[bm])
        fit = fit_2gauss(x1[bm])
        if fit is None: continue
        
        is_core = assign_components(x1[bm], fit)
        is_ext = ~is_core
        if np.sum(is_core) < 8 or np.sum(is_ext) < 8: continue
        
        X_core = np.column_stack([mb_std[bm][is_core], x1_std[bm][is_core], c_std[bm][is_core]])
        X_ext = np.column_stack([mb_std[bm][is_ext], x1_std[bm][is_ext], c_std[bm][is_ext]])
        
        C_c = np.cov(X_core.T)
        C_e = np.cov(X_ext.T)
        
        vals_c, vecs_c = np.linalg.eigh(C_c)
        vals_e, vecs_e = np.linalg.eigh(C_e)
        
        idx_c = np.argsort(vals_c)[::-1]
        idx_e = np.argsort(vals_e)[::-1]
        vals_c, vecs_c = vals_c[idx_c], vecs_c[:, idx_c]
        vals_e, vecs_e = vals_e[idx_e], vecs_e[:, idx_e]
        
        theta = angle_between_vectors(vecs_c[:, 0], vecs_e[:, 0])
        
        z_m.append(zc)
        angles_m.append(theta)
        lam1_c_m.append(vals_c[0])
        lam1_e_m.append(vals_e[0])
        loadings_c_m.append(np.abs(vecs_c[:, 0]))
        loadings_e_m.append(np.abs(vecs_e[:, 0]))
    
    if len(z_m) >= 2:
        mean_theta = np.mean(angles_m)
        mean_load_c = np.mean(loadings_c_m, axis=0)
        mean_load_e = np.mean(loadings_e_m, axis=0)
        print(f"    Mean rotation angle: {mean_theta:.1f}°")
        print(f"    Core PC1 loadings: mB={mean_load_c[0]:.3f} x1={mean_load_c[1]:.3f} c={mean_load_c[2]:.3f}")
        print(f"    Ext  PC1 loadings: mB={mean_load_e[0]:.3f} x1={mean_load_e[1]:.3f} c={mean_load_e[2]:.3f}")
        if len(z_m) >= 3:
            rho_ang, p_ang = spearmanr(z_m, angles_m)
            print(f"    Rotation angle vs z: ρ = {rho_ang:+.3f}, p = {p_ang:.3f}")
            rho_lam, p_lam = spearmanr(z_m, lam1_e_m)
            print(f"    Extended λ1 vs z: ρ = {rho_lam:+.3f}, p = {p_lam:.3f}")

# Does the rotation angle differ by host mass?
print(f"\n\n  METALLICITY × ROTATION INTERACTION:")
print(f"  Does host mass change the rotation angle between corridors?")
# Collect all-mass angles and compare
if len(angles_1) >= 3:
    print(f"  ALL hosts mean angle: {np.mean(angles_1):.1f}°")


# ============================================================
# PART C: METALLICITY AS HIDDEN VARIABLE
# ============================================================
print(f"\n\n{'=' * 120}")
print("PART C: Is Metallicity a Rotation Parameter or a Contraction Parameter?")
print("=" * 120)

# If metallicity ROTATES: the angle between corridors should change with host mass
# If metallicity CONTRACTS: the eigenvalues should change but angle stays constant

print(f"\n  Test: Does the ANGLE change with host mass (rotation) or just the EIGENVALUE (contraction)?")

# Compute angle and λ1 for different host mass bins
mass_bins_fine = np.percentile(host_mass[mass_valid & fast_mask], [0, 25, 50, 75, 100])

print(f"\n  {'Mass range':>20} | {'Mean angle':>11} | {'Mean λ1_ext':>12} | {'N bins':>7}")
print(f"  " + "-" * 60)

mass_angles = []
mass_lam1s = []
mass_centers = []

for q in range(len(mass_bins_fine) - 1):
    m_lo, m_hi = mass_bins_fine[q], mass_bins_fine[q+1]
    mcut = (host_mass >= m_lo) & (host_mass < m_hi) & mass_valid
    
    angles_q, lam1s_q = [], []
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bm = (z >= z_lo) & (z < z_hi) & fast_mask & mcut
        n = np.sum(bm)
        if n < 20: continue
        
        fit = fit_2gauss(x1[bm])
        if fit is None: continue
        
        is_core = assign_components(x1[bm], fit)
        is_ext = ~is_core
        if np.sum(is_core) < 6 or np.sum(is_ext) < 6: continue
        
        X_c = np.column_stack([mb_std[bm][is_core], x1_std[bm][is_core], c_std[bm][is_core]])
        X_e = np.column_stack([mb_std[bm][is_ext], x1_std[bm][is_ext], c_std[bm][is_ext]])
        
        vals_c, vecs_c = np.linalg.eigh(np.cov(X_c.T))
        vals_e, vecs_e = np.linalg.eigh(np.cov(X_e.T))
        
        idx_c = np.argsort(vals_c)[::-1]
        idx_e = np.argsort(vals_e)[::-1]
        
        theta = angle_between_vectors(vecs_c[:, idx_c[0]], vecs_e[:, idx_e[0]])
        angles_q.append(theta)
        lam1s_q.append(vals_e[idx_e[0]])
    
    if len(angles_q) >= 1:
        ma = np.mean(angles_q)
        ml = np.mean(lam1s_q)
        mc = (m_lo + m_hi) / 2
        mass_angles.append(ma)
        mass_lam1s.append(ml)
        mass_centers.append(mc)
        print(f"  {f'[{m_lo:.1f}, {m_hi:.1f})':>20} | {ma:>11.1f}° | {ml:>12.3f} | {len(angles_q):>7}")

if len(mass_centers) >= 3:
    rho_rot, p_rot = spearmanr(mass_centers, mass_angles)
    rho_con, p_con = spearmanr(mass_centers, mass_lam1s)
    print(f"\n  Host mass vs rotation angle: ρ = {rho_rot:+.3f}, p = {p_rot:.3f}")
    print(f"  Host mass vs λ1_extended:    ρ = {rho_con:+.3f}, p = {p_con:.3f}")
    
    if abs(rho_rot) > abs(rho_con):
        print(f"\n  → METALLICITY IS A ROTATION PARAMETER (changes the projection angle)")
        print(f"    Host mass doesn't just narrow — it TILTS how the latent process maps onto observables")
    elif abs(rho_con) > abs(rho_rot):
        print(f"\n  → METALLICITY IS A CONTRACTION PARAMETER (changes the eigenvalue magnitude)")
        print(f"    Host mass narrows the state space without changing the projection geometry")
    else:
        print(f"\n  → METALLICITY DOES BOTH (rotates AND contracts)")
