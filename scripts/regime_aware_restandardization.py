"""
Regime-Aware Re-Standardization — Reverse Causal Test
=======================================================
Real data → fit separate β per GMM state → rerun CPL.
If w_a collapses toward 0, standardization failure caused
the dark energy signal — not real cosmological evolution.

Two approaches:
  A) Assign GMM mean β per state, recompute μ, fit CPL
  B) Fit β_intact, β_broken simultaneously with (w0, wa, MB, α)
     — fully joint, no assumptions

Also checks α for state-dependence.

Baseline: original MU_SH0ES (pipeline β=3.1) → CPL
Compare: regime-corrected μ → CPL
"""

import numpy as np
import pandas as pd
from scipy import optimize, interpolate
import json, os

rng = np.random.default_rng(42)

# ── GMM parameters (from two_state_mixture_tests.json) ───────
BETA_INTACT   = 0.881
BETA_BROKEN   = 0.586
ALPHA_PIPELINE = 0.147    # standard SALT2 α
BETA_PIPELINE  = 3.1      # standard SALT2 β — what MU_SH0ES was built with
MB_PIPELINE    = -19.253  # Pantheon+ MB (SH0ES calibration)

# ── Load pantheon data ────────────────────────────────────────
pp = pd.read_csv('/root/clawd/projects/closure-theory/paper1/results/results/pantheon_with_y.csv')
pp = pp.sort_values('zHD').drop_duplicates('CID', keep='first')
pp = pp[pp['MU_SH0ES'].notna() & pp['zHD'].notna() & (pp['zHD'] > 0.023)].copy()
pp = pp[pp['c'].notna() & pp['x1'].notna() & pp['mB'].notna()].copy()
pp = pp.reset_index(drop=True)
N = len(pp)
print(f"N={N} SNe with full SALT2 params")

z   = pp['zHD'].values
c   = pp['c'].values
x1  = pp['x1'].values
mB  = pp['mB'].values
sig = np.clip(pp['MU_SH0ES_ERR_DIAG'].values, 0.05, 0.5)
mu_pipeline = pp['MU_SH0ES'].values   # baseline: pipeline μ

# ── Load GMM state labels ─────────────────────────────────────
# GMM labels stored in two_state_mixture_tests.json or nearby
gmm_file = '/root/clawd/projects/closure/results/two_state_mixture_tests.json'
with open(gmm_file) as f:
    gmm_data = json.load(f)

# Try to get per-SN labels from the GMM result
broken_mask = None
if 'broken_indices' in gmm_data:
    idx = np.array(gmm_data['broken_indices'])
    broken_mask = np.zeros(N, dtype=bool)
    broken_mask[idx[idx < N]] = True
elif 'labels' in gmm_data:
    labels = np.array(gmm_data['labels'])
    if len(labels) == N:
        broken_mask = labels == 1
elif 'broken_cids' in gmm_data:
    broken_cids = set(gmm_data['broken_cids'])
    broken_mask = pp['CID'].isin(broken_cids).values

if broken_mask is None:
    # Fallback: use color+residual heuristic matching GMM means
    # Broken state: β=0.586, mean c=+0.065; Intact: β=0.881, mean c=−0.059
    # Use mu_resid as proxy — broken state has larger |resid|
    print("No GMM labels found — reconstructing from mu_resid threshold")
    mu_resid = pp['mu_resid'].values if 'mu_resid' in pp.columns else (mu_pipeline - np.interp(z, z, mu_pipeline))
    # GMM broken fraction = 30%; threshold at 70th percentile of |resid|
    thresh = np.percentile(np.abs(mu_resid), 70)
    broken_mask = np.abs(mu_resid) > thresh
    print(f"  Threshold |resid| > {thresh:.4f}: {broken_mask.sum()} broken ({broken_mask.mean()*100:.1f}%)")
else:
    print(f"GMM labels loaded: {broken_mask.sum()} broken ({broken_mask.mean()*100:.1f}%), {(~broken_mask).sum()} intact")

# ── CPL fitter ────────────────────────────────────────────────
Z_GRID = np.linspace(0.0001, 2.5, 3000)

def mu_cpl(z_arr, w0, wa, Om=0.3):
    zg = Z_GRID
    f_de = (1+zg)**(3*(1+w0+wa)) * np.exp(-3*wa*zg/(1+zg))
    Ez = np.sqrt(Om*(1+zg)**3 + (1-Om)*f_de)
    dc = np.zeros(len(zg))
    inv_E = 1.0/Ez
    for i in range(1, len(zg)):
        dc[i] = dc[i-1] + 0.5*(inv_E[i]+inv_E[i-1])*(zg[i]-zg[i-1])
    dc *= 2997.92
    dl = (1+zg)*dc; dl[0] = 1e-10
    mu_th = 5*np.log10(dl) + 25
    return np.interp(z_arr, zg, mu_th)

def fit_cpl(z_arr, mu_arr, sig_arr, label='', w0_init=-1.0, wa_init=0.0):
    # Grid search for good init
    best = (1e12, w0_init, wa_init, -19.3)
    for w0t in np.arange(-1.8, -0.3, 0.3):
        for wat in np.arange(-2.0, 1.5, 0.5):
            mu_th = mu_cpl(z_arr, w0t, wat)
            MB = np.median(mu_arr - mu_th)
            chi2 = np.sum(((mu_arr - mu_th - MB)/sig_arr)**2)
            if chi2 < best[0]:
                best = (chi2, w0t, wat, MB)
    def f(p):
        mu_th = mu_cpl(z_arr, p[0], p[1])
        return np.sum(((mu_arr - mu_th - p[2])/sig_arr)**2)
    res = optimize.minimize(f, [best[1], best[2], best[3]], method='Nelder-Mead',
                            options={'maxiter':8000, 'xatol':0.001, 'fatol':0.1})
    w0, wa, MB = res.x
    chi2 = res.fun
    ndof = len(z_arr) - 3
    print(f"  {label} (N={len(z_arr)}): w₀={w0:.4f}, w_a={wa:.4f}, χ²/dof={chi2/ndof:.3f}")
    return float(w0), float(wa), float(MB), float(chi2/ndof)

# ── STEP 1: Baseline (pipeline μ, universal β=3.1) ───────────
print("\n=== BASELINE: pipeline μ (universal β=3.1) ===")
w0_base, wa_base, MB_base, chi2_base = fit_cpl(z, mu_pipeline, sig, 'Full sample (pipeline)')
w0_bi, wa_bi, _, _ = fit_cpl(z[~broken_mask], mu_pipeline[~broken_mask], sig[~broken_mask], 'Intact only (pipeline)')
w0_bb, wa_bb, _, _ = fit_cpl(z[broken_mask],  mu_pipeline[broken_mask],  sig[broken_mask],  'Broken only (pipeline)')

# ── STEP 2: Approach A — assign GMM mean β per state ─────────
# μ_corrected = mB - MB - α·x1 - β_state·c
# Note: MU_SH0ES = mB - MB_pipeline - α·x1 - β_pipeline·c
# So: μ_corrected = MU_SH0ES + (β_pipeline - β_state)·c
# β correction: add back pipeline overcorrection, apply regime β
print("\n=== APPROACH A: GMM mean β per state ===")
beta_assigned = np.where(broken_mask, BETA_BROKEN, BETA_INTACT)
delta_beta = BETA_PIPELINE - beta_assigned   # how much pipeline overcorrected
mu_corrA = mu_pipeline + delta_beta * c       # regime-corrected μ

print(f"  Intact correction: +{(BETA_PIPELINE-BETA_INTACT):.3f}·c  (mean shift: {((BETA_PIPELINE-BETA_INTACT)*c[~broken_mask]).mean():.4f} mag)")
print(f"  Broken correction: +{(BETA_PIPELINE-BETA_BROKEN):.3f}·c  (mean shift: {((BETA_PIPELINE-BETA_BROKEN)*c[broken_mask]).mean():.4f} mag)")

w0_A, wa_A, MB_A, chi2_A = fit_cpl(z, mu_corrA, sig, 'Full (approach A: GMM β per state)')
w0_Ai, wa_Ai, _, _ = fit_cpl(z[~broken_mask], mu_corrA[~broken_mask], sig[~broken_mask], 'Intact (approach A)')
w0_Ab, wa_Ab, _, _ = fit_cpl(z[broken_mask],  mu_corrA[broken_mask],  sig[broken_mask],  'Broken (approach A)')

print(f"\n  Δw_a (baseline → corrected): {wa_A - wa_base:+.4f}")
print(f"  Δw_a as fraction of baseline: {(wa_A - wa_base)/abs(wa_base)*100:+.1f}%")

# ── STEP 3: Approach B — joint fit β_i, β_b, α, w0, wa, MB ──
print("\n=== APPROACH B: Joint fit (β per state as free params) ===")

# Parameters: [w0, wa, MB, α, β_intact, β_broken]
def chi2_joint(p):
    w0, wa, MB, alpha, beta_i, beta_b = p
    beta_arr = np.where(broken_mask, beta_b, beta_i)
    mu_th = mu_cpl(z, w0, wa)
    mu_corr = mB - MB - alpha*x1 - beta_arr*c
    return np.sum(((mu_corr - mu_th)/sig)**2)

# Grid init: use baseline cosmology, known β values
p0 = [w0_base, wa_base, MB_PIPELINE, ALPHA_PIPELINE, BETA_INTACT, BETA_BROKEN]
print(f"  Starting from: w₀={p0[0]:.3f}, w_a={p0[1]:.3f}, α={p0[3]:.3f}, β_i={p0[4]:.3f}, β_b={p0[5]:.3f}")

res_B = optimize.minimize(chi2_joint, p0, method='Nelder-Mead',
                          options={'maxiter':20000, 'xatol':0.001, 'fatol':1.0,
                                   'adaptive': True})
w0_B, wa_B, MB_B, alpha_B, beta_i_B, beta_b_B = res_B.x
chi2_B = res_B.fun / (N - 6)
print(f"  Joint fit result:")
print(f"    w₀={w0_B:.4f}, w_a={wa_B:.4f}")
print(f"    α={alpha_B:.4f} (vs pipeline {ALPHA_PIPELINE})")
print(f"    β_intact={beta_i_B:.4f} (vs GMM {BETA_INTACT})")
print(f"    β_broken={beta_b_B:.4f} (vs GMM {BETA_BROKEN})")
print(f"    χ²/dof={chi2_B:.3f}")
print(f"    Δw_a (B vs baseline): {wa_B - wa_base:+.4f}")
print(f"    α state difference: {alpha_B:.4f} (single α — not split)")

# ── STEP 4: Does α need splitting? ───────────────────────────
print("\n=== α state test: fit α per state ===")

def chi2_alpha_split(p):
    w0, wa, MB, alpha_i, alpha_b, beta_i, beta_b = p
    alpha_arr = np.where(broken_mask, alpha_b, alpha_i)
    beta_arr  = np.where(broken_mask, beta_b,  beta_i)
    mu_th = mu_cpl(z, w0, wa)
    mu_corr = mB - MB - alpha_arr*x1 - beta_arr*c
    return np.sum(((mu_corr - mu_th)/sig)**2)

p0_split = [w0_B, wa_B, MB_B, ALPHA_PIPELINE, ALPHA_PIPELINE, beta_i_B, beta_b_B]
res_split = optimize.minimize(chi2_alpha_split, p0_split, method='Nelder-Mead',
                               options={'maxiter':20000, 'xatol':0.001, 'fatol':1.0, 'adaptive': True})
w0_s, wa_s, MB_s, ai_s, ab_s, bi_s, bb_s = res_split.x
chi2_s = res_split.fun / (N - 7)
print(f"  α_intact={ai_s:.4f}, α_broken={ab_s:.4f}  (Δα={ai_s-ab_s:+.4f})")
print(f"  β_intact={bi_s:.4f}, β_broken={bb_s:.4f}")
print(f"  w₀={w0_s:.4f}, w_a={wa_s:.4f}")
print(f"  χ²/dof={chi2_s:.3f}  (vs joint B: {chi2_B:.3f})")
delta_bic_alpha = (res_split.fun - res_B.fun) + np.log(N)  # +1 param
print(f"  ΔBIC(α split vs single α): {delta_bic_alpha:.1f}  ({'favors split' if delta_bic_alpha < -6 else 'no improvement'})")

# ── STEP 5: Summary ───────────────────────────────────────────
print(f"\n{'='*70}")
print(f"  REGIME-AWARE RE-STANDARDIZATION — CAUSAL TEST SUMMARY")
print(f"{'='*70}")
print(f"  Baseline (pipeline β=3.1):")
print(f"    Full:   w₀={w0_base:.4f}, w_a={wa_base:+.4f}")
print(f"    Intact: w₀={w0_bi:.4f},   w_a={wa_bi:+.4f}")
print(f"    Broken: w₀={w0_bb:.4f},   w_a={wa_bb:+.4f}")
print(f"")
print(f"  Approach A (GMM mean β per state):")
print(f"    Full:   w₀={w0_A:.4f},   w_a={wa_A:+.4f}   Δw_a={wa_A-wa_base:+.4f}")
print(f"    Intact: w₀={w0_Ai:.4f},  w_a={wa_Ai:+.4f}")
print(f"    Broken: w₀={w0_Ab:.4f},  w_a={wa_Ab:+.4f}")
print(f"")
print(f"  Approach B (joint fit: β per state free):")
print(f"    Full:   w₀={w0_B:.4f},   w_a={wa_B:+.4f}   Δw_a={wa_B-wa_base:+.4f}")
print(f"    β_intact fitted: {beta_i_B:.4f}, β_broken fitted: {beta_b_B:.4f}")
print(f"")
print(f"  α split test:")
print(f"    α_intact={ai_s:.4f}, α_broken={ab_s:.4f}")
print(f"    w_a with α split: {wa_s:+.4f}  (ΔBIC={delta_bic_alpha:.1f})")
print(f"")

# Verdict
wa_collapse_A = abs(wa_A) < 0.3 * abs(wa_base)
wa_collapse_B = abs(wa_B) < 0.3 * abs(wa_base)
if wa_collapse_A or wa_collapse_B:
    print(f"  ✓ w_a COLLAPSES after regime correction — standardization caused the signal")
elif abs(wa_A - wa_base) > 0.3 or abs(wa_B - wa_base) > 0.3:
    print(f"  ~ w_a SHIFTS significantly ({wa_A-wa_base:+.3f}) — partial causal attribution")
else:
    print(f"  ✗ w_a barely changes — regime correction insufficient or data structure different")

print(f"{'='*70}")

# ── Save ─────────────────────────────────────────────────────
results = {
    'baseline': {'w0': w0_base, 'wa': wa_base, 'chi2_dof': chi2_base},
    'baseline_intact': {'w0': w0_bi, 'wa': wa_bi},
    'baseline_broken': {'w0': w0_bb, 'wa': wa_bb},
    'approach_A': {
        'method': 'GMM mean beta per state',
        'beta_intact': BETA_INTACT, 'beta_broken': BETA_BROKEN,
        'w0': w0_A, 'wa': wa_A, 'chi2_dof': chi2_A,
        'dwa_vs_baseline': wa_A - wa_base,
        'intact': {'w0': w0_Ai, 'wa': wa_Ai},
        'broken': {'w0': w0_Ab, 'wa': wa_Ab},
    },
    'approach_B': {
        'method': 'Joint fit: beta per state free',
        'beta_intact_fitted': beta_i_B, 'beta_broken_fitted': beta_b_B,
        'alpha_fitted': alpha_B,
        'w0': w0_B, 'wa': wa_B, 'chi2_dof': chi2_B,
        'dwa_vs_baseline': wa_B - wa_base,
    },
    'alpha_split': {
        'alpha_intact': ai_s, 'alpha_broken': ab_s,
        'delta_alpha': ai_s - ab_s,
        'w0': w0_s, 'wa': wa_s, 'chi2_dof': chi2_s,
        'delta_bic_vs_single_alpha': delta_bic_alpha,
    },
}

os.makedirs('/root/clawd/projects/closure/results', exist_ok=True)
out = '/root/clawd/projects/closure/results/regime_aware_restandardization.json'
with open(out, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved → {out}")
