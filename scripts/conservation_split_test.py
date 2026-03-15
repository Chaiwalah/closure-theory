"""
Conservation Test — Is the Two-State Split the Skeleton?
=========================================================
After regime-aware correction collapses w_a to 0, test whether
the broken/intact split accounts for ~90% of the signal structure
or only ~70% — if ~70%, there's a third anatomical component.

Tests:
  1. Signal conservation: does w_a(broken) + w_a(intact) weighted sum
     equal the full-sample w_a? How much is conserved vs leaked?

  2. Residual structure after approach B: are the corrected residuals
     Gaussian white noise (→ 2 states sufficient) or structured
     (→ 3rd component)?

  3. 3-state GMM vs 2-state: ΔBIC tells us if a third cluster is real

  4. Explained variance: what % of the breaker/sky signal is captured
     by the GMM two-state label vs remaining in residuals?

  5. Broken-state anatomy: within the broken population, is there
     sub-structure (bimodal resid? sky clustering?)
"""

import numpy as np
import pandas as pd
from scipy import stats, optimize, interpolate
from sklearn.mixture import GaussianMixture
import json, os, warnings
warnings.filterwarnings('ignore')

rng = np.random.default_rng(42)

# ── Load data ─────────────────────────────────────────────────
pp = pd.read_csv('/root/clawd/projects/closure-theory/paper1/results/results/pantheon_with_y.csv')
pp = pp.sort_values('zHD').drop_duplicates('CID', keep='first')
pp = pp[pp['MU_SH0ES'].notna() & pp['zHD'].notna() & (pp['zHD'] > 0.023)].copy()
pp = pp[pp['c'].notna() & pp['x1'].notna() & pp['mB'].notna()].copy()
pp = pp.reset_index(drop=True)
N = len(pp)

z         = pp['zHD'].values
c         = pp['c'].values
x1        = pp['x1'].values
mB        = pp['mB'].values
mu_pipe   = pp['MU_SH0ES'].values
sig       = np.clip(pp['MU_SH0ES_ERR_DIAG'].values, 0.05, 0.5)
mu_resid  = pp['mu_resid'].values if 'mu_resid' in pp.columns else np.zeros(N)

# ── Reconstruct broken mask (same threshold as prev run) ─────
thresh = np.percentile(np.abs(mu_resid), 70)
broken_mask = np.abs(mu_resid) > thresh
n_b = broken_mask.sum(); n_i = (~broken_mask).sum()
print(f"N={N}  Broken={n_b} ({n_b/N*100:.1f}%)  Intact={n_i}")

# ── CPL fitter (reused) ───────────────────────────────────────
Z_GRID = np.linspace(0.0001, 2.5, 3000)

def mu_cpl_arr(z_arr, w0, wa, Om=0.3):
    zg = Z_GRID
    f_de = (1+zg)**(3*(1+w0+wa)) * np.exp(-3*wa*zg/(1+zg))
    Ez = np.sqrt(Om*(1+zg)**3 + (1-Om)*f_de)
    dc = np.zeros(len(zg))
    inv_E = 1/Ez
    for i in range(1, len(zg)): dc[i] = dc[i-1] + 0.5*(inv_E[i]+inv_E[i-1])*(zg[i]-zg[i-1])
    dc *= 2997.92
    dl = (1+zg)*dc; dl[0]=1e-10
    mu_th = 5*np.log10(dl)+25
    return np.interp(z_arr, zg, mu_th)

def fit_cpl_quick(z_arr, mu_arr, sig_arr):
    best = (1e12, -1.0, 0.0, -19.3)
    for w0t in np.arange(-1.8,-0.3,0.3):
        for wat in np.arange(-2.0,1.5,0.5):
            mu_th = mu_cpl_arr(z_arr,w0t,wat)
            MB = np.median(mu_arr-mu_th)
            chi2 = np.sum(((mu_arr-mu_th-MB)/sig_arr)**2)
            if chi2 < best[0]: best=(chi2,w0t,wat,MB)
    def f(p): return np.sum(((mu_arr-mu_cpl_arr(z_arr,p[0],p[1])-p[2])/sig_arr)**2)
    r = optimize.minimize(f, [best[1],best[2],best[3]], method='Nelder-Mead',
                          options={'maxiter':6000,'xatol':0.001,'fatol':0.1})
    return float(r.x[0]), float(r.x[1])

# ── TEST 1: Signal conservation ───────────────────────────────
print("\n=== TEST 1: Signal conservation across the split ===")
f_b = n_b / N
f_i = n_i / N

w0_full, wa_full   = fit_cpl_quick(z, mu_pipe, sig)
w0_i,    wa_i      = fit_cpl_quick(z[~broken_mask], mu_pipe[~broken_mask], sig[~broken_mask])
w0_b,    wa_b      = fit_cpl_quick(z[broken_mask],  mu_pipe[broken_mask],  sig[broken_mask])

wa_predicted = f_i*wa_i + f_b*wa_b   # simple weighted sum

# Conservation metric: how much of the broken-intact w_a SPREAD is
# explained by the two-state split?
# Δw_a signal = wa_broken - wa_full (how much broken state drags the full sample)
# Predicted drag = f_b * (wa_broken - wa_intact)
# Observed drag = wa_full - wa_intact (how much full differs from intact-only)
observed_drag  = wa_full - wa_i           # how much broken contaminates full
predicted_drag = f_b * (wa_b - wa_i)     # what the weighted mix predicts
conservation = 1 - abs(predicted_drag - observed_drag) / (abs(observed_drag) + 1e-6)
total_signal = abs(wa_b - wa_i)          # total broken-intact separation

print(f"  w_a full sample:             {wa_full:+.4f}")
print(f"  w_a intact ({f_i*100:.0f}%):         {wa_i:+.4f}")
print(f"  w_a broken ({f_b*100:.0f}%):          {wa_b:+.4f}")
print(f"  Broken-intact separation:    {total_signal:+.4f}")
print(f"  Observed drag on full:       {observed_drag:+.4f}")
print(f"  Predicted drag (weighted):   {predicted_drag:+.4f}")
print(f"  Conservation:                {conservation*100:.1f}%")

if conservation > 0.88:
    print(f"  → ≥88%: TWO STATES ARE THE SKELETON. Residuals are noise.")
elif conservation > 0.70:
    print(f"  → 70-88%: PARTIAL. A third component exists in the anatomy.")
else:
    print(f"  → <70%: SIGNIFICANT THIRD COMPONENT. Two states insufficient.")

# ── TEST 2: Residual structure after approach B correction ────
print("\n=== TEST 2: Residual structure post-correction ===")
BETA_INTACT = 0.881; BETA_BROKEN = 0.586; BETA_PIPE = 3.1
beta_assigned = np.where(broken_mask, BETA_BROKEN, BETA_INTACT)
mu_corrA = mu_pipe + (BETA_PIPE - beta_assigned)*c

# Fit ΛCDM to corrected distances
mu_th_lcdm = mu_cpl_arr(z, -1.0, 0.0)
MB_lcdm = np.median(mu_corrA - mu_th_lcdm)
resid_corrected = mu_corrA - mu_th_lcdm - MB_lcdm

# Normality test
stat_full, p_norm_full = stats.normaltest(resid_corrected)
stat_b,    p_norm_b    = stats.normaltest(resid_corrected[broken_mask])
stat_i,    p_norm_i    = stats.normaltest(resid_corrected[~broken_mask])

print(f"  Corrected residuals — normality test:")
print(f"    Full:   p={p_norm_full:.4f}  ({'NORMAL ✓' if p_norm_full>0.05 else 'NON-NORMAL — structure remains'})")
print(f"    Intact: p={p_norm_i:.4f}  ({'NORMAL ✓' if p_norm_i>0.05 else 'NON-NORMAL'})")
print(f"    Broken: p={p_norm_b:.4f}  ({'NORMAL ✓' if p_norm_b>0.05 else 'NON-NORMAL — 3rd component?'})")

# Skewness/kurtosis
print(f"  Skewness  — full: {stats.skew(resid_corrected):.3f}, broken: {stats.skew(resid_corrected[broken_mask]):.3f}")
print(f"  Kurtosis  — full: {stats.kurtosis(resid_corrected):.3f}, broken: {stats.kurtosis(resid_corrected[broken_mask]):.3f}")

# Is there z-dependence remaining in broken residuals?
rho_z_b, p_z_b = stats.spearmanr(z[broken_mask], resid_corrected[broken_mask])
rho_z_i, p_z_i = stats.spearmanr(z[~broken_mask], resid_corrected[~broken_mask])
print(f"  Residual vs z:  broken ρ={rho_z_b:.3f} p={p_z_b:.4f}, intact ρ={rho_z_i:.3f} p={p_z_i:.4f}")

# ── TEST 3: 3-state GMM vs 2-state ───────────────────────────
print("\n=== TEST 3: 3-state GMM vs 2-state ===")

# Features: same as original GMM — use c, x1, mu_resid
features = np.column_stack([c, x1, mu_resid])
# normalize
features_n = (features - features.mean(0)) / (features.std(0) + 1e-10)

results_gmm = {}
for n_comp in [2, 3, 4]:
    gmm = GaussianMixture(n_components=n_comp, covariance_type='full',
                          random_state=42, n_init=5, max_iter=500)
    gmm.fit(features_n)
    bic = gmm.bic(features_n)
    aic = gmm.aic(features_n)
    results_gmm[n_comp] = {'bic': bic, 'aic': aic, 'labels': gmm.predict(features_n)}
    print(f"  {n_comp}-state GMM: BIC={bic:.1f}, AIC={aic:.1f}")

delta_bic_3v2 = results_gmm[3]['bic'] - results_gmm[2]['bic']
delta_bic_4v3 = results_gmm[4]['bic'] - results_gmm[3]['bic']
print(f"  ΔBIC(3 vs 2): {delta_bic_3v2:.1f}  ({'3-state strongly favored' if delta_bic_3v2 < -10 else '3-state marginally favored' if delta_bic_3v2 < -2 else 'No improvement — 2 states sufficient'})")
print(f"  ΔBIC(4 vs 3): {delta_bic_4v3:.1f}  ({'4-state favored' if delta_bic_4v3 < -10 else 'No improvement'})")

# If 3-state favored, characterize the third state
if delta_bic_3v2 < -6:
    labels3 = results_gmm[3]['labels']
    print(f"\n  3-state characterization:")
    for s in range(3):
        mask_s = labels3 == s
        print(f"    State {s} (N={mask_s.sum()}, {mask_s.mean()*100:.0f}%):")
        print(f"      mean c={c[mask_s].mean():.3f}, mean x1={x1[mask_s].mean():.3f}")
        print(f"      mean |resid|={np.abs(mu_resid[mask_s]).mean():.4f}")
        print(f"      mean z={z[mask_s].mean():.3f}")
    # What is the 3rd state? Check overlap with broken_mask
    overlap = np.zeros(3)
    for s in range(3):
        overlap[s] = broken_mask[labels3==s].mean()
    third_state = np.argmin(np.abs(overlap - 0.5))  # most ambiguous → third
    print(f"    Broken fraction per state: {[f'{overlap[s]:.2f}' for s in range(3)]}")

# ── TEST 4: Within-broken substructure ───────────────────────
print("\n=== TEST 4: Broken-state anatomy ===")
resid_broken = mu_resid[broken_mask]
c_broken     = c[broken_mask]
z_broken     = z[broken_mask]

# Bimodality coefficient
n_b2 = len(resid_broken)
skew_b = stats.skew(resid_broken)
kurt_b = stats.kurtosis(resid_broken)
bc = (skew_b**2 + 1) / (kurt_b + 3*(n_b2-1)**2/((n_b2-2)*(n_b2-3)))
print(f"  Broken residuals bimodality coefficient: {bc:.3f}  (>0.555 = bimodal)")

# Split broken into positive/negative resid
pos_mask = resid_broken > 0
neg_mask = resid_broken <= 0
print(f"  Overcorrected (resid>0): N={pos_mask.sum()}, mean c={c_broken[pos_mask].mean():.3f}, mean z={z_broken[pos_mask].mean():.3f}")
print(f"  Undercorrected (resid<0): N={neg_mask.sum()}, mean c={c_broken[neg_mask].mean():.3f}, mean z={z_broken[neg_mask].mean():.3f}")

# Do overcorrected and undercorrected have different color profiles?
rho_sign_c, p_sign_c = stats.pointbiserialr(pos_mask.astype(int), c_broken)
rho_sign_z, p_sign_z = stats.pointbiserialr(pos_mask.astype(int), z_broken)
print(f"  Resid sign vs c: r={rho_sign_c:.3f}, p={p_sign_c:.4f}")
print(f"  Resid sign vs z: r={rho_sign_z:.3f}, p={p_sign_z:.4f}")

# GMM 2-state on broken alone
gmm_b = GaussianMixture(n_components=2, covariance_type='full', random_state=42, n_init=5)
feat_b = np.column_stack([c_broken, resid_broken, z_broken])
feat_b_n = (feat_b - feat_b.mean(0)) / (feat_b.std(0)+1e-10)
gmm_b.fit(feat_b_n)
bic_b2 = gmm_b.bic(feat_b_n)
gmm_b1 = GaussianMixture(n_components=1, covariance_type='full', random_state=42)
gmm_b1.fit(feat_b_n)
bic_b1 = gmm_b1.bic(feat_b_n)
delta_bic_broken = bic_b2 - bic_b1
print(f"  Within-broken GMM ΔBIC(2 vs 1): {delta_bic_broken:.1f}  ({'sub-states present' if delta_bic_broken < -10 else 'no sub-structure'})")

if delta_bic_broken < -10:
    sub_labels = gmm_b.predict(feat_b_n)
    for s in range(2):
        ms = sub_labels==s
        print(f"    Sub-state {s} (N={ms.sum()}): mean c={c_broken[ms].mean():.3f}, mean resid={resid_broken[ms].mean():.4f}, mean z={z_broken[ms].mean():.3f}")

# ── SUMMARY ───────────────────────────────────────────────────
print(f"\n{'='*65}")
print(f"  CONSERVATION TEST — SUMMARY")
print(f"{'='*65}")
print(f"  Signal conservation: {conservation*100:.1f}%")
print(f"  Residual normality (broken): p={p_norm_b:.4f}")
print(f"  ΔBIC 3-state vs 2-state: {delta_bic_3v2:.1f}")
print(f"  Within-broken ΔBIC: {delta_bic_broken:.1f}")
print(f"  Broken bimodality coefficient: {bc:.3f} (threshold 0.555)")
print(f"")
if conservation > 0.88 and delta_bic_3v2 > -6:
    print(f"  VERDICT: TWO STATES ARE THE SKELETON ({conservation*100:.0f}% conservation)")
    print(f"           The third thing is noise, not anatomy.")
elif conservation > 0.70 and (delta_bic_3v2 < -6 or bc > 0.555):
    print(f"  VERDICT: THIRD COMPONENT EXISTS ({conservation*100:.0f}% conservation, ΔBIC={delta_bic_3v2:.1f})")
    print(f"           Broken state has sub-structure: overcorrected vs undercorrected.")
elif conservation < 0.70 and delta_bic_3v2 < -6:
    print(f"  VERDICT: THREE-STATE ANATOMY ({conservation*100:.0f}% conservation for 2 states)")
    print(f"           Third state is real — anatomy, not noise.")
else:
    print(f"  VERDICT: AMBIGUOUS — conservation={conservation*100:.0f}%, ΔBIC={delta_bic_3v2:.1f}")
print(f"{'='*65}")

# Save
out = {
    'conservation': float(conservation),
    'total_signal_dwa': float(total_signal),
    'observed_drag': float(observed_drag),
    'predicted_drag': float(predicted_drag),
    'wa_full': float(wa_full), 'wa_intact': float(wa_i), 'wa_broken': float(wa_b),
    'wa_weighted_sum': float(wa_predicted),
    'residual_normality': {'full_p': float(p_norm_full), 'broken_p': float(p_norm_b), 'intact_p': float(p_norm_i)},
    'gmm_bic': {str(k): float(v['bic']) for k,v in results_gmm.items()},
    'delta_bic_3v2': float(delta_bic_3v2), 'delta_bic_4v3': float(delta_bic_4v3),
    'within_broken': {'bimodality_coeff': float(bc), 'delta_bic_2v1': float(delta_bic_broken)},
    'broken_anatomy': {
        'overcorrected_n': int(pos_mask.sum()), 'undercorrected_n': int(neg_mask.sum()),
        'resid_vs_c_r': float(rho_sign_c), 'resid_vs_c_p': float(p_sign_c),
        'resid_vs_z_r': float(rho_sign_z), 'resid_vs_z_p': float(p_sign_z),
    }
}
os.makedirs('/root/clawd/projects/closure/results', exist_ok=True)
with open('/root/clawd/projects/closure/results/conservation_split_test.json','w') as f:
    json.dump(out, f, indent=2)
print("Saved → conservation_split_test.json")
