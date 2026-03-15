"""
w_a Erasure Test v4 — correct version
======================================
Uses MU_SH0ES (raw distance modulus) directly, not residuals.
GMM labels from (c, mu_resid) as before.
Fits CPL w(z)=w0+wa*z/(1+z) to μ_obs vs z for full/intact/broken.
No residual inversion — correct by construction.
Jackknife uncertainty (5-fold, ~2min).
"""

import numpy as np
import pandas as pd
from scipy import optimize, interpolate
from sklearn.mixture import GaussianMixture
import json, time

# ── Fast CPL lookup table ─────────────────────────────────────
Z_GRID = np.linspace(0.0001, 2.5, 3000)
OM = 0.3

def build_mu_interp(w0, wa, Om=OM):
    zg = Z_GRID
    f_de = (1+zg)**(3*(1+w0+wa)) * np.exp(-3*wa*zg/(1+zg))
    Ez = np.sqrt(Om*(1+zg)**3 + (1-Om)*f_de)
    dc = np.zeros(len(zg))
    inv_E = 1.0/Ez
    for i in range(1, len(zg)):
        dc[i] = dc[i-1] + 0.5*(inv_E[i]+inv_E[i-1])*(zg[i]-zg[i-1])
    dc *= 2997.92  # c/H0, H0=70
    dl = (1+zg)*dc; dl[0] = 1e-10
    mu = 5*np.log10(dl) + 25
    return interpolate.interp1d(zg, mu, bounds_error=False, fill_value='extrapolate')

# ── Load & clean ─────────────────────────────────────────────
pp = pd.read_csv('/root/clawd/projects/closure-theory/paper1/results/results/pantheon_with_y.csv')
pp = pp.groupby('CID').agg({
    'zHD':'mean','c':'mean','x1':'mean',
    'mu_resid':'mean',      # for GMM labels only
    'MU_SH0ES':'mean',      # raw distance modulus — the fit target
    'MU_SH0ES_ERR_DIAG':'mean',
    'RA':'mean','DEC':'mean'
}).reset_index()

pp = pp[
    (pp['zHD'] > 0.023) &
    pp['c'].notna() &
    pp['MU_SH0ES'].notna() &
    pp['mu_resid'].notna() &
    (pp['MU_SH0ES'] > 0)   # drop calibrators / bad values
].copy()

print(f"N = {len(pp)} SNe")
print(f"z range: {pp['zHD'].min():.4f} – {pp['zHD'].max():.4f}")
print(f"MU_SH0ES range: {pp['MU_SH0ES'].min():.2f} – {pp['MU_SH0ES'].max():.2f}")

# ── GMM labels (from residuals, as before) ───────────────────
X = pp[['c','mu_resid']].values
gmm = GaussianMixture(n_components=2, covariance_type='full', random_state=42, n_init=20)
gmm.fit(X)
labels = gmm.predict(X)

# Pin: broken = higher |mu_resid| AND redder c
mu0 = pp[labels==0]['mu_resid'].abs().mean()
mu1 = pp[labels==1]['mu_resid'].abs().mean()
c0  = pp[labels==0]['c'].mean()
c1  = pp[labels==1]['c'].mean()
if (mu0 > mu1) + (c0 > c1) < (mu1 > mu0) + (c1 > c0):
    labels = 1 - labels

pp['state'] = labels
broken = pp[pp['state']==0]
intact = pp[pp['state']==1]
print(f"Broken: {len(broken)} ({len(broken)/len(pp)*100:.1f}%)  Intact: {len(intact)} ({len(intact)/len(pp)*100:.1f}%)")
print(f"  Broken: mean|mu_resid|={broken['mu_resid'].abs().mean():.4f}, mean c={broken['c'].mean():.4f}")
print(f"  Intact: mean|mu_resid|={intact['mu_resid'].abs().mean():.4f}, mean c={intact['c'].mean():.4f}")

# ── CPL fitter using MU_SH0ES directly ──────────────────────
def fit_cpl(df, label='', sigma=None):
    t0 = time.time()
    z   = df['zHD'].values
    mu  = df['MU_SH0ES'].values
    sig = df['MU_SH0ES_ERR_DIAG'].values if sigma is None else np.full(len(z), sigma)
    sig = np.clip(sig, 0.05, 2.0)  # guard against zeros

    # MB (absolute mag offset) is a free parameter — not H0-dependent for w fit
    best = {'chi2': 1e12, 'p': [-1.0, 0.0, -19.3]}
    for w0t in np.arange(-1.8, -0.2, 0.25):
        for wat in np.arange(-2.5, 1.5, 0.5):
            mu_th = build_mu_interp(w0t, wat)(z)
            MB = np.median(mu - mu_th)
            chi2 = np.sum(((mu - mu_th - MB)/sig)**2)
            if chi2 < best['chi2']:
                best = {'chi2': chi2, 'p': [w0t, wat, MB]}

    def f(p):
        mu_th = build_mu_interp(p[0], p[1])(z)
        return np.sum(((mu - mu_th - p[2])/sig)**2)

    res = optimize.minimize(f, best['p'], method='Nelder-Mead',
                            options={'maxiter':5000,'xatol':0.003,'fatol':0.1})
    w0, wa, MB = res.x
    chi2 = res.fun
    dof = len(z) - 3
    print(f"  {label} (N={len(z)}): w₀={w0:.4f}, w_a={wa:.4f}, χ²/dof={chi2/dof:.3f}  [{time.time()-t0:.1f}s]")
    return float(w0), float(wa), float(chi2), int(dof)

print("\n--- Fitting ---")
w0_f, wa_f, chi2_f, dof_f = fit_cpl(pp,     'Full sample')
w0_i, wa_i, chi2_i, dof_i = fit_cpl(intact, 'Intact only')
w0_b, wa_b, chi2_b, dof_b = fit_cpl(broken, 'Broken only')

# ── 5-fold jackknife ─────────────────────────────────────────
print("\n5-fold jackknife on Δw_a...")
np.random.seed(42)
idx_f = np.random.permutation(len(pp))
idx_i = np.random.permutation(len(intact))
intact_reset = intact.reset_index(drop=True)
pp_reset = pp.reset_index(drop=True)

dwa_jk = []
for fold in range(5):
    mf = np.ones(len(pp), dtype=bool); mf[idx_f[fold*len(pp)//5:(fold+1)*len(pp)//5]] = False
    mi = np.ones(len(intact), dtype=bool); mi[idx_i[fold*len(intact)//5:(fold+1)*len(intact)//5]] = False
    
    def quick(df, w0s, was):
        z_b = df['zHD'].values; mu_b = df['MU_SH0ES'].values
        sig = np.clip(df['MU_SH0ES_ERR_DIAG'].values, 0.05, 2.0)
        mu_th = build_mu_interp(w0s, was)(z_b)
        MB = np.median(mu_b - mu_th)
        def f(p): return np.sum(((mu_b - build_mu_interp(p[0],p[1])(z_b) - p[2])/sig)**2)
        r = optimize.minimize(f, [w0s, was, MB], method='Nelder-Mead',
                              options={'maxiter':2000,'xatol':0.008,'fatol':1.0})
        return r.x[0], r.x[1]

    _, wa_jf = quick(pp_reset[mf], w0_f, wa_f)
    _, wa_ji = quick(intact_reset[mi], w0_i, wa_i)
    dwa_jk.append(wa_ji - wa_jf)
    print(f"  Fold {fold+1}: Δw_a = {dwa_jk[-1]:.3f}")

dwa_mean = np.mean(dwa_jk)
dwa_std  = np.std(dwa_jk) * np.sqrt(4)
sigma_jk = abs(dwa_mean)/dwa_std if dwa_std > 0 else 0

# ── Summary ───────────────────────────────────────────────────
delta_w0 = w0_i - w0_f
delta_wa = wa_i - wa_f
pct_wa = (1 - abs(wa_i)/abs(wa_f))*100 if abs(wa_f) > 0.1 else None

print(f"\n{'='*62}")
print(f"  w_a ERASURE TEST v4 — FINAL RESULT")
print(f"{'='*62}")
print(f"  Full sample  (N={len(pp):4d}):  w₀={w0_f:+.4f},  w_a={wa_f:+.4f}")
print(f"  Intact only  (N={len(intact):4d}):  w₀={w0_i:+.4f},  w_a={wa_i:+.4f}")
print(f"  Broken only  (N={len(broken):4d}):  w₀={w0_b:+.4f},  w_a={wa_b:+.4f}")
print(f"")
print(f"  Δw₀ (intact − full) = {delta_w0:+.4f}")
print(f"  Δw_a (intact − full) = {delta_wa:+.4f}")
print(f"  Jackknife: {dwa_mean:+.4f} ± {dwa_std:.4f}  ({sigma_jk:.1f}σ)")
if pct_wa is not None:
    print(f"  % of w_a signal from broken regime: {-pct_wa:.1f}%")
print(f"  Intact alone shifts w_a toward 0? {'YES ✓' if abs(wa_i) < abs(wa_f) else 'NO — intact alone moves AWAY from ΛCDM'}")
print(f"  Broken alone moves OPPOSITE direction? {'YES ✓' if np.sign(wa_b) != np.sign(wa_i) else 'NO'}")
print(f"{'='*62}")

results = {
    'note': 'v4 — uses MU_SH0ES (raw distance modulus), not residuals',
    'N_full': len(pp), 'N_intact': len(intact), 'N_broken': len(broken),
    'full':   {'w0': w0_f, 'wa': wa_f, 'chi2_dof': chi2_f/dof_f},
    'intact': {'w0': w0_i, 'wa': wa_i, 'chi2_dof': chi2_i/dof_i},
    'broken': {'w0': w0_b, 'wa': wa_b, 'chi2_dof': chi2_b/dof_b},
    'delta_w0': delta_w0, 'delta_wa': delta_wa,
    'jk_mean': float(dwa_mean), 'jk_std': float(dwa_std), 'jk_sigma': float(sigma_jk),
}

import os
os.makedirs('/root/clawd/projects/closure/results', exist_ok=True)
with open('/root/clawd/projects/closure/results/wa_erasure_test_v4.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Saved → wa_erasure_test_v4.json")
