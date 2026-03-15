"""
Mock Catalog Injection — Causal Test
======================================
Generate a strict ΛCDM universe. Inject broken standardizers.
Run through a blind CPL fitter that assumes one universal β.
If w_a ≠ 0 comes out of a w_a=0 input → standardization failure
is sufficient to generate the dark energy illusion.

Key design choices vs GPT's prescription:
- β drawn from full GMM broken distribution (mean=0.586, std=0.146)
  NOT fixed at 0.58 — bidirectional overcorrection AND undercorrection
- Redshift distribution matched to real Pantheon+ sample
- canonical β_pipeline = 3.1 (SALT2 default)
- 3 runs: full mock, intact-only, broken-only
"""

import numpy as np
import pandas as pd
from scipy import optimize, interpolate
import json, os

rng = np.random.default_rng(42)

# ── 1. Load real redshift distribution ───────────────────────
pp = pd.read_csv('/root/clawd/projects/closure-theory/paper1/results/results/pantheon_with_y.csv')
pp = pp.groupby('CID').agg({'zHD':'mean','MU_SH0ES_ERR_DIAG':'mean'}).reset_index()
pp = pp[(pp['zHD'] > 0.023) & pp['MU_SH0ES_ERR_DIAG'].notna()].copy()
z_real = pp['zHD'].values
sigma_real = np.clip(pp['MU_SH0ES_ERR_DIAG'].values, 0.05, 0.5)
N = len(z_real)
print(f"Mock catalog: N={N}, z range [{z_real.min():.3f}, {z_real.max():.3f}]")

# ── 2. ΛCDM true distance moduli ─────────────────────────────
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
    dc *= 2997.92
    dl = (1+zg)*dc; dl[0] = 1e-10
    mu = 5*np.log10(dl) + 25
    return interpolate.interp1d(zg, mu, bounds_error=False, fill_value='extrapolate')

mu_true_fn = build_mu_interp(-1.0, 0.0)   # strict ΛCDM
mu_true = mu_true_fn(z_real)
MB_TRUE = -19.3

# ── 3. Assign intrinsic colors ───────────────────────────────
# Real Pantheon+ color distribution: mean~0, std~0.09 (intact), ~0.12 (broken)
c_intact = rng.normal(0.0,  0.09, N)
c_broken_draw = rng.normal(0.065, 0.12, N)   # broken state is redder (from GMM)

# ── 4. Label 30% broken ──────────────────────────────────────
broken_mask = rng.random(N) < 0.30
n_broken = broken_mask.sum()
n_intact = N - n_broken
print(f"Broken: {n_broken} ({n_broken/N*100:.1f}%)  Intact: {n_intact}")

# ── 5. Assign TRUE β per SN ──────────────────────────────────
# Intact: β_true ~ 3.1 (standard, coupling holds)
# Broken: β_true drawn from GMM broken distribution
#         mean=0.586, std=0.146 — bidirectional, not fixed
BETA_PIPELINE = 3.1    # what SALT2 assumes for everyone
beta_true = np.full(N, BETA_PIPELINE)
beta_broken_draw = rng.normal(0.586, 0.146, N)
beta_true[broken_mask] = beta_broken_draw[broken_mask]

# Colors: broken SNe draw from redder distribution
c_obs = np.where(broken_mask, c_broken_draw, c_intact)

# ── 6. Generate observed magnitudes ──────────────────────────
# True: m_b = MB + mu_true + β_true * c + noise
# Pipeline applies: mu_obs = m_b - MB - β_pipeline * c
# So: mu_obs = mu_true + (β_true - β_pipeline)*c + noise

noise = rng.normal(0, sigma_real)
delta_beta = beta_true - BETA_PIPELINE          # 0 for intact, ~-2.5 for broken mean
mu_obs = mu_true + MB_TRUE + delta_beta * c_obs + noise   # what pipeline produces
# Shift so MB absorbed (fitter will fit MB as free param)
mu_obs_corr = mu_obs - MB_TRUE    # ~= mu_true + delta_beta*c + noise

print(f"\nδβ stats:")
print(f"  Intact: mean={delta_beta[~broken_mask].mean():.3f}, std={delta_beta[~broken_mask].std():.3f}")
print(f"  Broken: mean={delta_beta[broken_mask].mean():.3f}, std={delta_beta[broken_mask].std():.3f}")
print(f"  Broken mu_obs bias: mean={( delta_beta*c_obs)[broken_mask].mean():.4f}")

# ── 7. CPL fitter ─────────────────────────────────────────────
def fit_cpl(z, mu, sig, label=''):
    best = {'chi2': 1e12, 'p': [-1.0, 0.0, -19.3]}
    for w0t in np.arange(-1.8, -0.2, 0.3):
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
                            options={'maxiter':5000,'xatol':0.002,'fatol':0.05})
    w0, wa, MB = res.x
    chi2 = res.fun
    print(f"  {label} (N={len(z)}): w₀={w0:.4f}, w_a={wa:.4f}, χ²/dof={chi2/(len(z)-3):.3f}")
    return float(w0), float(wa)

print("\n--- Mock CPL fits (input truth: w₀=-1.00, w_a=0.00) ---")
w0_full,  wa_full  = fit_cpl(z_real,             mu_obs_corr,             sigma_real, 'Full mock (all SNe)')
w0_intact,wa_intact= fit_cpl(z_real[~broken_mask],mu_obs_corr[~broken_mask],sigma_real[~broken_mask],'Intact only (truth SNe)')
w0_broken,wa_broken= fit_cpl(z_real[broken_mask], mu_obs_corr[broken_mask], sigma_real[broken_mask], 'Broken only (bad standardizers)')

# ── 8. 5-fold jackknife on Δw_a ──────────────────────────────
print("\n5-fold jackknife on recovered Δw_a (intact−full)...")
idx_all = np.random.permutation(N)
idx_int = np.where(~broken_mask)[0]; np.random.shuffle(idx_int)
dwa_jk = []
for fold in range(5):
    mf = np.ones(N, dtype=bool)
    mf[idx_all[fold*N//5:(fold+1)*N//5]] = False
    ni = len(idx_int)
    mi_idx = idx_int.copy()
    mi = np.ones(N, dtype=bool)
    mi[mi_idx[fold*ni//5:(fold+1)*ni//5]] = False
    mi_intact = mi[~broken_mask]
    
    def quick(z, mu, sig, w0s, was):
        def f(p): return np.sum(((mu - build_mu_interp(p[0],p[1])(z) - p[2])/sig)**2)
        r = optimize.minimize(f, [w0s, was, -19.3], method='Nelder-Mead',
                              options={'maxiter':2000,'xatol':0.008,'fatol':1.0})
        return r.x[1]
    
    waf = quick(z_real[mf], mu_obs_corr[mf], sigma_real[mf], w0_full, wa_full)
    wai = quick(z_real[~broken_mask][mi_intact],
                mu_obs_corr[~broken_mask][mi_intact],
                sigma_real[~broken_mask][mi_intact], w0_intact, wa_intact)
    dwa_jk.append(wai - waf)
    print(f"  Fold {fold+1}: Δw_a = {dwa_jk[-1]:.3f}")

dwa_mean = np.mean(dwa_jk)
dwa_std  = np.std(dwa_jk) * np.sqrt(4)
sigma_jk = abs(dwa_mean)/dwa_std if dwa_std > 0 else 0

# ── 9. Summary ────────────────────────────────────────────────
print(f"\n{'='*65}")
print(f"  MOCK INJECTION CAUSAL TEST — FINAL RESULT")
print(f"{'='*65}")
print(f"  TRUE cosmology injected:   w₀=-1.00, w_a= 0.00 (strict ΛCDM)")
print(f"  Broken fraction injected:  {n_broken/N*100:.1f}% with β~N(0.586,0.146)")
print(f"  Pipeline β assumed:        {BETA_PIPELINE} (SALT2 universal)")
print(f"")
print(f"  Recovered (full mock):     w₀={w0_full:+.4f}, w_a={wa_full:+.4f}")
print(f"  Recovered (intact only):   w₀={w0_intact:+.4f}, w_a={wa_intact:+.4f}")
print(f"  Recovered (broken only):   w₀={w0_broken:+.4f}, w_a={wa_broken:+.4f}")
print(f"")
print(f"  Δw_a (intact−full):        {wa_intact-wa_full:+.4f}")
print(f"  Jackknife: {dwa_mean:+.4f} ± {dwa_std:.4f}  ({sigma_jk:.1f}σ)")
print(f"")
if abs(wa_full) > 0.3:
    print(f"  ✓ ΛCDM universe → false w_a signal confirmed ({wa_full:+.3f})")
else:
    print(f"  ✗ No false signal generated (w_a={wa_full:+.3f} — check β distribution)")
print(f"{'='*65}")

results = {
    'note': 'Causal test — ΛCDM in, broken standardizers injected, CPL fitter blind',
    'truth': {'w0': -1.0, 'wa': 0.0},
    'injection': {'broken_frac': float(n_broken/N), 'beta_mean': 0.586, 'beta_std': 0.146, 'beta_pipeline': BETA_PIPELINE},
    'recovered_full':   {'w0': w0_full,   'wa': wa_full},
    'recovered_intact': {'w0': w0_intact, 'wa': wa_intact},
    'recovered_broken': {'w0': w0_broken, 'wa': wa_broken},
    'delta_wa': float(wa_intact - wa_full),
    'jk_mean': float(dwa_mean), 'jk_std': float(dwa_std), 'jk_sigma': float(sigma_jk),
}

os.makedirs('/root/clawd/projects/closure/results', exist_ok=True)
with open('/root/clawd/projects/closure/results/mock_injection_causal_test.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Saved → mock_injection_causal_test.json")
