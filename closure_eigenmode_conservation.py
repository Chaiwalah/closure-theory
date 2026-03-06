#!/usr/bin/env python3
"""
closure_eigenmode_conservation.py — Tests for Agents to Chew On
================================================================

TEST 1: Is Tr(Λ) conserved? (GPT's formulation)
  - Compute total eigenvalue sum at each z-bin for core and extended
  - If conserved → conservation is in eigenspace, not observable space

TEST 2: Does PC1 eigenvalue stay constant while observed c-variance changes?
  - If yes → rotation is redistributing variance, not creating/destroying it

TEST 3: PC1 loading sign structure
  - Is it pure c, or mixed c+mB+x1? 
  - Sign pattern tells you if it's a pure thermodynamic axis or mixed

TEST 4: Rotation angle vs z — does it evolve? (by host mass)
  - Low-Z: predict angle grows with z (Gemini)
  - High-Z: predict angle stable (Grok)

TEST 5: Gaia bridge — does color-dependent PC1 loading match
  the frequency fingerprint from Feb 20?
  - Split ALL SNe by color (blue vs red) and compute PC1 separately
  - If loading structure rotates with SN color → same as Gaia chromatic PM

TEST 6: Slow channel PCA — is it the same story?
  - Run identical PCA on slow channel (x1≥0)
  - If PC1 = color there too → universal, not fast-channel specific
  - If different → corridor-specific

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

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 2.5]
fast_mask = x1 < 0
slow_mask = x1 >= 0
mass_valid = host_mass > 0

# Standardize on ALL SNe
mb_mean, mb_std_v = np.mean(mb), np.std(mb)
x1_mean, x1_std_v = np.mean(x1), np.std(x1)
c_mean, c_std_v = np.mean(c), np.std(c)
mb_s = (mb - mb_mean) / mb_std_v
x1_s = (x1 - x1_mean) / x1_std_v
c_s = (c - c_mean) / c_std_v

def fit_2gauss(x):
    best, best_ll = None, 1e10
    for w0 in [0.3, 0.5, 0.7]:
        for m1 in [-0.3, -0.5]:
            for m2 in [-1.0, -1.5]:
                try:
                    def nll(p):
                        w,a,s1,b,s2 = p
                        if w<0.01 or w>0.99 or s1<0.01 or s2<0.01: return 1e10
                        return -np.sum(np.log(np.maximum(w*norm.pdf(x,a,s1)+(1-w)*norm.pdf(x,b,s2),1e-30)))
                    r = minimize(nll,[w0,m1,0.3,m2,0.5],method='Nelder-Mead',options={'maxiter':5000})
                    if r.fun < best_ll: best_ll, best = r.fun, r.x
                except: pass
    if best is None: return None
    w,a,s1,b,s2 = best; s1,s2=abs(s1),abs(s2)
    if a < b: w,a,s1,b,s2 = 1-w,b,s2,a,s1
    return w,a,s1,b,s2

def pca_bin(mb_b, x1_b, c_b):
    X = np.column_stack([mb_b, x1_b, c_b])
    C = np.cov(X.T)
    vals, vecs = np.linalg.eigh(C)
    idx = np.argsort(vals)[::-1]
    return vals[idx], vecs[:, idx]

def angle_between(v1, v2):
    cos_a = np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1)
    return np.degrees(np.arccos(abs(cos_a)))


# ============================================================
# TEST 1: Is Tr(Λ) conserved across z?
# ============================================================
print("=" * 110)
print("TEST 1: Is Tr(Λ) conserved? (Total eigenvalue sum across z)")
print("=" * 110)

for label, ch_mask in [('FAST (x1<0)', fast_mask), ('SLOW (x1≥0)', slow_mask), ('ALL', np.ones(len(z), bool))]:
    z_l, tr_l = [], []
    for i in range(len(z_edges)-1):
        bm = (z >= z_edges[i]) & (z < z_edges[i+1]) & ch_mask
        if np.sum(bm) < 20: continue
        vals, _ = pca_bin(mb_s[bm], x1_s[bm], c_s[bm])
        z_l.append(np.mean(z[bm]))
        tr_l.append(np.sum(vals))
    if len(z_l) >= 3:
        rho, p = spearmanr(z_l, tr_l)
        cv = np.std(tr_l)/np.mean(tr_l)*100
        print(f"  {label:>15}: Tr(Λ) vs z: ρ={rho:+.3f} p={p:.3f} CV={cv:.1f}% {'🔥 CONSERVED' if abs(rho)<0.4 else ''}")


# ============================================================
# TEST 2: PC1 eigenvalue vs observed c-variance
# ============================================================
print(f"\n{'=' * 110}")
print("TEST 2: Does PC1 eigenvalue stay constant while Var(c) changes?")
print("=" * 110)

z_l2, lam1_l, varc_l = [], [], []
for i in range(len(z_edges)-1):
    bm = (z >= z_edges[i]) & (z < z_edges[i+1]) & fast_mask
    if np.sum(bm) < 20: continue
    vals, _ = pca_bin(mb_s[bm], x1_s[bm], c_s[bm])
    z_l2.append(np.mean(z[bm]))
    lam1_l.append(vals[0])
    varc_l.append(np.var(c[bm]))

if len(z_l2) >= 3:
    rho_l, p_l = spearmanr(z_l2, lam1_l)
    rho_v, p_v = spearmanr(z_l2, varc_l)
    rho_lv, p_lv = spearmanr(lam1_l, varc_l)
    print(f"  λ1 vs z:     ρ={rho_l:+.3f} p={p_l:.3f}")
    print(f"  Var(c) vs z: ρ={rho_v:+.3f} p={p_v:.3f}")
    print(f"  λ1 vs Var(c): ρ={rho_lv:+.3f} p={p_lv:.3f}")
    if abs(rho_l) < 0.4 and abs(rho_v) > 0.4:
        print(f"  → λ1 STABLE while Var(c) changes → rotation redistributes, doesn't create/destroy")


# ============================================================
# TEST 3: PC1 loading sign structure
# ============================================================
print(f"\n{'=' * 110}")
print("TEST 3: PC1 Loading Sign Structure — Pure c or mixed?")
print("=" * 110)

print(f"\n  {'z':>6} {'N':>4} | {'mB':>8} {'x1':>8} {'c':>8} | dominant | {'λ1/Σλ':>6}")
print(f"  " + "-" * 60)
for i in range(len(z_edges)-1):
    bm = (z >= z_edges[i]) & (z < z_edges[i+1]) & fast_mask
    n = np.sum(bm)
    if n < 20: continue
    vals, vecs = pca_bin(mb_s[bm], x1_s[bm], c_s[bm])
    v1 = vecs[:, 0]
    labels = ['mB', 'x1', 'c']
    dom = labels[np.argmax(np.abs(v1))]
    frac = vals[0] / np.sum(vals)
    print(f"  {np.mean(z[bm]):>6.3f} {n:>4} | {v1[0]:>+8.3f} {v1[1]:>+8.3f} {v1[2]:>+8.3f} | {dom:>8} | {frac:>6.1%}")


# ============================================================
# TEST 4: Rotation angle vs z, split by host mass
# ============================================================
print(f"\n{'=' * 110}")
print("TEST 4: Rotation Angle vs z — by Host Mass")
print("=" * 110)

mass_med = np.median(host_mass[mass_valid & fast_mask])
for mlabel, mcut in [('LOW-mass', (host_mass < mass_med) & mass_valid),
                      ('HIGH-mass', (host_mass >= mass_med) & mass_valid)]:
    z_l4, ang_l4 = [], []
    for i in range(len(z_edges)-1):
        bm = (z >= z_edges[i]) & (z < z_edges[i+1]) & fast_mask & mcut
        if np.sum(bm) < 25: continue
        fit = fit_2gauss(x1[bm])
        if fit is None: continue
        w,a,s1,b,s2 = fit
        pc = w*norm.pdf(x1[bm],a,s1); pe = (1-w)*norm.pdf(x1[bm],b,s2)
        is_core = (pc/(pc+pe)) > 0.5
        if np.sum(is_core)<8 or np.sum(~is_core)<8: continue
        _, vecs_c = pca_bin(mb_s[bm][is_core], x1_s[bm][is_core], c_s[bm][is_core])
        _, vecs_e = pca_bin(mb_s[bm][~is_core], x1_s[bm][~is_core], c_s[bm][~is_core])
        ang = angle_between(vecs_c[:,0], vecs_e[:,0])
        z_l4.append(np.mean(z[bm]))
        ang_l4.append(ang)
    if len(z_l4) >= 3:
        rho, p = spearmanr(z_l4, ang_l4)
        print(f"  {mlabel:>12}: angle vs z: ρ={rho:+.3f} p={p:.3f} | angles: {['%.1f°'%a for a in ang_l4]}")


# ============================================================
# TEST 5: Color-split PCA (Gaia bridge)
# ============================================================
print(f"\n{'=' * 110}")
print("TEST 5: Gaia Bridge — Does PC1 rotate with SN color?")
print("  Split fast channel by color (blue vs red), compare PC1 loadings")
print("=" * 110)

c_med = np.median(c[fast_mask])
print(f"  Color median (fast): {c_med:.4f}")

for clabel, ccut in [('BLUE (c<med)', c < c_med), ('RED (c≥med)', c >= c_med)]:
    all_vecs = []
    for i in range(len(z_edges)-1):
        bm = (z >= z_edges[i]) & (z < z_edges[i+1]) & fast_mask & ccut
        if np.sum(bm) < 20: continue
        _, vecs = pca_bin(mb_s[bm], x1_s[bm], c_s[bm])
        all_vecs.append(vecs[:, 0])
    if all_vecs:
        mean_load = np.mean([np.abs(v) for v in all_vecs], axis=0)
        print(f"  {clabel:>20}: PC1 |loadings| mB={mean_load[0]:.3f} x1={mean_load[1]:.3f} c={mean_load[2]:.3f}")

# Angle between blue and red PC1
blue_vecs, red_vecs = [], []
z_l5, ang_l5 = [], []
for i in range(len(z_edges)-1):
    bm_b = (z >= z_edges[i]) & (z < z_edges[i+1]) & fast_mask & (c < c_med)
    bm_r = (z >= z_edges[i]) & (z < z_edges[i+1]) & fast_mask & (c >= c_med)
    if np.sum(bm_b) < 15 or np.sum(bm_r) < 15: continue
    _, vb = pca_bin(mb_s[bm_b], x1_s[bm_b], c_s[bm_b])
    _, vr = pca_bin(mb_s[bm_r], x1_s[bm_r], c_s[bm_r])
    ang = angle_between(vb[:,0], vr[:,0])
    z_l5.append(np.mean(z[(z >= z_edges[i]) & (z < z_edges[i+1]) & fast_mask]))
    ang_l5.append(ang)

if ang_l5:
    print(f"\n  Blue vs Red PC1 rotation: {['%.1f°'%a for a in ang_l5]}")
    print(f"  Mean: {np.mean(ang_l5):.1f}° (cf. Gaia chromatic PM: 8-13°)")


# ============================================================
# TEST 6: Slow Channel PCA — Same or Different?
# ============================================================
print(f"\n{'=' * 110}")
print("TEST 6: Slow Channel PCA — Is PC1 also color?")
print("=" * 110)

print(f"\n  {'z':>6} {'N':>4} | {'mB':>8} {'x1':>8} {'c':>8} | dominant | {'λ1/Σλ':>6}")
print(f"  " + "-" * 60)
for i in range(len(z_edges)-1):
    bm = (z >= z_edges[i]) & (z < z_edges[i+1]) & slow_mask
    n = np.sum(bm)
    if n < 20: continue
    vals, vecs = pca_bin(mb_s[bm], x1_s[bm], c_s[bm])
    v1 = vecs[:, 0]
    labels = ['mB', 'x1', 'c']
    dom = labels[np.argmax(np.abs(v1))]
    frac = vals[0] / np.sum(vals)
    print(f"  {np.mean(z[bm]):>6.3f} {n:>4} | {v1[0]:>+8.3f} {v1[1]:>+8.3f} {v1[2]:>+8.3f} | {dom:>8} | {frac:>6.1%}")

# Compare fast vs slow PC1
print(f"\n  Mean PC1 loadings:")
fast_loads, slow_loads = [], []
for i in range(len(z_edges)-1):
    bm_f = (z >= z_edges[i]) & (z < z_edges[i+1]) & fast_mask
    bm_s = (z >= z_edges[i]) & (z < z_edges[i+1]) & slow_mask
    if np.sum(bm_f) < 20 or np.sum(bm_s) < 20: continue
    _, vf = pca_bin(mb_s[bm_f], x1_s[bm_f], c_s[bm_f])
    _, vs = pca_bin(mb_s[bm_s], x1_s[bm_s], c_s[bm_s])
    fast_loads.append(np.abs(vf[:,0]))
    slow_loads.append(np.abs(vs[:,0]))

if fast_loads and slow_loads:
    mf = np.mean(fast_loads, axis=0)
    ms = np.mean(slow_loads, axis=0)
    print(f"  FAST: mB={mf[0]:.3f} x1={mf[1]:.3f} c={mf[2]:.3f}")
    print(f"  SLOW: mB={ms[0]:.3f} x1={ms[1]:.3f} c={ms[2]:.3f}")
    ang_fs = angle_between(mf, ms)
    print(f"  Angle between fast/slow PC1: {ang_fs:.1f}°")
