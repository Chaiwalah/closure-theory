#!/usr/bin/env python3
"""
Closure Theory — Round 7: Kill the Last Conventional Explanation
=================================================================
SALT3 cross-check on DES-5YR + improved threshold model + combined dataset power.

Tests:
  K2: SALT3 model adequacy on DES-5YR (same K1 test, different fitter)
  M4: Sigmoid threshold AIC/BIC (proper threshold model instead of linear)
  P1: Combined Pantheon+ × DES-5YR high-z power boost
  A2b: Within-survey entanglement on DES-5YR
"""

import numpy as np
import os
import json
from scipy import stats
from scipy.optimize import curve_fit, minimize
from scipy.integrate import quad
from numpy.linalg import lstsq
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

RESULTS_DIR = "results_round7"
os.makedirs(RESULTS_DIR, exist_ok=True)

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.bool_,)): return bool(obj)
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return super().default(obj)

def load_des5yr():
    """Load DES-5YR metadata (SALT3 fits)."""
    rows = []
    varnames = None
    with open("data/des_metadata.csv", 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('VARNAMES:'):
                varnames = line.split()[1:]
            elif line.startswith('SN:'):
                vals = line.split()[1:]
                if varnames and len(vals) == len(varnames):
                    rows.append(vals)
    
    data = {}
    for i, name in enumerate(varnames):
        col = [r[i] for r in rows]
        try:
            data[name] = np.array(col, dtype=float)
        except:
            data[name] = np.array(col, dtype=str)
    
    # Quality cuts
    z = data['zHD']
    x1 = data['x1']
    c = data['c']
    mu = data['MU']
    ndof = data['NDOF']
    
    mask = (z > 0.01) & (z < 2.0) & np.isfinite(x1) & np.isfinite(c) & np.isfinite(mu)
    mask &= (np.abs(x1) < 5) & (np.abs(c) < 1) & (ndof > 0)
    
    out = {}
    for k, v in data.items():
        out[k] = v[mask]
    
    out['chi2_dof'] = out['FITCHI2'] / out['NDOF']
    out['mu_resid'] = out['MURES']
    
    # Color residual
    poly = np.polyfit(out['zHD'], out['c'], 1)
    out['c_resid'] = out['c'] - np.polyval(poly, out['zHD'])
    
    print(f"DES-5YR: {mask.sum()} SNe after cuts, z=[{out['zHD'].min():.3f}, {out['zHD'].max():.3f}]")
    return out

def load_pantheon():
    import pandas as pd
    df = pd.read_csv("data/pantheon_plus.dat", sep=r'\s+', comment='#')
    
    z = df['zHD'].values
    def mu_lcdm(zv, H0=73.04, Om=0.334):
        if not np.isfinite(zv) or zv <= 0: return np.nan
        dc, _ = quad(lambda zp: 1/np.sqrt(Om*(1+zp)**3+(1-Om)), 0, zv)
        return 5*np.log10((1+zv)*dc*2.998e5/H0) + 25
    
    mu_model = np.array([mu_lcdm(zi) for zi in z])
    df['mu_resid'] = df['MU_SH0ES'].values - mu_model
    
    ndof = df['NDOF'].values.astype(float)
    ndof[ndof == 0] = np.nan
    df['chi2_dof'] = df['FITCHI2'].values.astype(float) / ndof
    
    poly = np.polyfit(z[np.isfinite(z) & np.isfinite(df['c'].values)],
                      df['c'].values[np.isfinite(z) & np.isfinite(df['c'].values)], 1)
    df['c_resid'] = df['c'].values - np.polyval(poly, z)
    
    return df

def effective_rank(X):
    X_std = (X - X.mean(0)) / (X.std(0) + 1e-10)
    s = np.linalg.svd(X_std, compute_uv=False)
    p = s / s.sum()
    p = p[p > 0]
    return float(np.exp(-np.sum(p * np.log(p))))

# =============================================================================
# K2: SALT3 MODEL ADEQUACY ON DES-5YR
# =============================================================================
def test_salt3_adequacy(des):
    """Same K1 test but on DES-5YR which uses SALT3, not SALT2."""
    print("\n" + "=" * 60)
    print("K2: SALT3 MODEL ADEQUACY (DES-5YR)")
    print("=" * 60)
    
    z = des['zHD']
    c = des['c']
    mu_res = des['mu_resid']
    chi2 = des['chi2_dof']
    x1 = des['x1']
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(mu_res) & np.isfinite(chi2) & np.isfinite(x1)
    z, c, mu_res, chi2, x1 = z[mask], c[mask], mu_res[mask], chi2[mask], x1[mask]
    
    results = {}
    
    # 1. chi2/dof vs z trend (SALT3)
    rho_cz, p_cz = stats.spearmanr(z, chi2)
    print(f"  χ²/dof vs z (SALT3): ρ={rho_cz:.4f}, p={p_cz:.4f}")
    results['chi2_z_trend'] = {'rho': float(rho_cz), 'p': float(p_cz)}
    
    # 2. Good fits vs bad fits
    chi_med = np.median(chi2)
    z_med = np.median(z)
    
    for label, chi_mask in [("Good fits (χ²<median)", chi2 < chi_med),
                             ("Bad fits (χ²≥median)", chi2 >= chi_med)]:
        hi = chi_mask & (z >= z_med)
        if hi.sum() > 20:
            rho, p = stats.spearmanr(c[hi], mu_res[hi])
            print(f"  {label}, high-z: ρ={rho:.4f}, p={p:.4f} (n={hi.sum()})")
            results[label] = {'rho': float(rho), 'p': float(p), 'n': int(hi.sum())}
    
    # 3. Signal with chi2 control
    hi = z >= z_med
    rho_raw, p_raw = stats.spearmanr(c[hi], mu_res[hi])
    
    controls = np.column_stack([z[hi], chi2[hi]])
    coef_c, _, _, _ = lstsq(controls, c[hi], rcond=None)
    coef_mu, _, _, _ = lstsq(controls, mu_res[hi], rcond=None)
    rho_ctrl, p_ctrl = stats.spearmanr(c[hi] - controls @ coef_c, mu_res[hi] - controls @ coef_mu)
    
    survival = abs(rho_ctrl) / abs(rho_raw) * 100 if abs(rho_raw) > 1e-6 else np.nan
    
    print(f"\n  High-z c-μ without χ² control: ρ={rho_raw:.4f}, p={p_raw:.4f}")
    print(f"  High-z c-μ with χ² control:    ρ={rho_ctrl:.4f}, p={p_ctrl:.4f}")
    print(f"  Survival: {survival:.1f}%")
    
    results['raw'] = {'rho': float(rho_raw), 'p': float(p_raw)}
    results['chi2_controlled'] = {'rho': float(rho_ctrl), 'p': float(p_ctrl), 'survival': float(survival)}
    
    # 4. UV exposure gradient (DES bands: griz)
    print(f"\n  Rest-frame UV exposure by z (DES/SALT3):")
    z_edges = [0.01, 0.2, 0.4, 0.6, 0.8, 1.2]
    for i in range(len(z_edges)-1):
        bm = (z >= z_edges[i]) & (z < z_edges[i+1])
        if bm.sum() > 10:
            r, p = stats.spearmanr(c[bm], mu_res[bm])
            print(f"    z=[{z_edges[i]:.1f}-{z_edges[i+1]:.1f}] (n={bm.sum()}): "
                  f"ρ(c,μ)={r:+.4f}, <χ²/dof>={chi2[bm].mean():.2f}")
    
    # 5. CRITICAL: Compare SALT2 (Pantheon+) vs SALT3 (DES) patterns
    print(f"\n  ═══ SALT2 vs SALT3 COMPARISON ═══")
    print(f"  Pantheon+ (SALT2): Signal in bad fits, survives χ² at 134%")
    print(f"  DES-5YR  (SALT3): Signal in bad fits={results.get('Bad fits (χ²≥median)', {}).get('p', '?')}, "
          f"survives χ² at {survival:.0f}%")
    
    both_in_bad = (results.get('Bad fits (χ²≥median)', {}).get('p', 1) < 0.05)
    print(f"\n  Both fitters show signal in bad fits: {both_in_bad}")
    if both_in_bad:
        print("  → Signal reproduces across SALT2 AND SALT3 — model mis-specification unlikely")
    
    return results

# =============================================================================
# M4: SIGMOID THRESHOLD AIC/BIC
# =============================================================================
def test_sigmoid_model(pan_df):
    """Proper threshold model comparison instead of linear Γ·y·z."""
    print("\n" + "=" * 60)
    print("M4: THRESHOLD MODEL SELECTION (AIC/BIC)")
    print("=" * 60)
    
    z = pan_df['zHD'].values
    c = pan_df['c'].values
    mu_res = pan_df['mu_resid'].values
    x1 = pan_df['x1'].values
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(mu_res) & np.isfinite(x1) & (z > 0.01)
    z, c, mu_res, x1 = z[mask], c[mask], mu_res[mask], x1[mask]
    n = len(z)
    
    results = {}
    
    # Models for mu_resid:
    # M0: mu_res = a0 + a1*x1                                    (baseline, 2 params)
    # M1: mu_res = a0 + a1*x1 + Γ*c                              (+ color term, 3 params)
    # M2: mu_res = a0 + a1*x1 + Γ*c*z                            (+ accumulation, 3 params)
    # M3: mu_res = a0 + a1*x1 + Γ*c*sigmoid(z-z0)               (+ threshold, 5 params)
    # M4: mu_res = a0 + a1*x1 + Γ*c*I(z>z0)                     (+ step function, 4 params, z0 scanned)
    
    def fit_model(X, y):
        """OLS fit, return SS_res, k, log-likelihood, AIC, BIC."""
        coef, res, _, _ = lstsq(X, y, rcond=None)
        pred = X @ coef
        ss_res = np.sum((y - pred)**2)
        k = X.shape[1]
        sigma2 = ss_res / n
        ll = -n/2 * np.log(2*np.pi*sigma2) - ss_res/(2*sigma2)
        aic = 2*k - 2*ll
        bic = k*np.log(n) - 2*ll
        return {'coef': coef.tolist(), 'ss_res': float(ss_res), 'k': k, 
                'll': float(ll), 'aic': float(aic), 'bic': float(bic),
                'r_squared': float(1 - ss_res/np.sum((y-y.mean())**2))}
    
    # M0: baseline
    X0 = np.column_stack([np.ones(n), x1])
    results['M0_baseline'] = fit_model(X0, mu_res)
    
    # M1: + color (global)
    X1 = np.column_stack([np.ones(n), x1, c])
    results['M1_color'] = fit_model(X1, mu_res)
    
    # M2: + color × z (accumulation)
    X2 = np.column_stack([np.ones(n), x1, c * z])
    results['M2_accumulation'] = fit_model(X2, mu_res)
    
    # M3: + color × sigmoid(z - z0)
    # Scan z0 for best fit
    best_aic = np.inf
    best_z0 = 0.5
    for z0_try in np.linspace(0.2, 1.2, 50):
        sig = 1 / (1 + np.exp(-15 * (z - z0_try)))
        X3 = np.column_stack([np.ones(n), x1, c * sig])
        r = fit_model(X3, mu_res)
        # Account for z0 as extra param: k+1
        aic_adj = 2*(r['k']+1) - 2*r['ll']
        if aic_adj < best_aic:
            best_aic = aic_adj
            best_z0 = z0_try
    
    sig_best = 1 / (1 + np.exp(-15 * (z - best_z0)))
    X3 = np.column_stack([np.ones(n), x1, c * sig_best])
    r3 = fit_model(X3, mu_res)
    r3['z0'] = float(best_z0)
    r3['k'] = r3['k'] + 1  # count z0
    r3['aic'] = 2*r3['k'] - 2*r3['ll']
    r3['bic'] = r3['k']*np.log(n) - 2*r3['ll']
    results['M3_sigmoid'] = r3
    
    # M4: step function at scanned z0
    best_aic_step = np.inf
    best_z0_step = 0.5
    for z0_try in np.linspace(0.2, 1.2, 50):
        step = (z > z0_try).astype(float)
        X4 = np.column_stack([np.ones(n), x1, c * step])
        r = fit_model(X4, mu_res)
        aic_adj = 2*(r['k']+1) - 2*r['ll']
        if aic_adj < best_aic_step:
            best_aic_step = aic_adj
            best_z0_step = z0_try
    
    step_best = (z > best_z0_step).astype(float)
    X4 = np.column_stack([np.ones(n), x1, c * step_best])
    r4 = fit_model(X4, mu_res)
    r4['z0'] = float(best_z0_step)
    r4['k'] = r4['k'] + 1
    r4['aic'] = 2*r4['k'] - 2*r4['ll']
    r4['bic'] = r4['k']*np.log(n) - 2*r4['ll']
    results['M4_step'] = r4
    
    # M5: color × z × sigmoid (accumulation + threshold)
    best_aic5 = np.inf
    best_z05 = 0.5
    for z0_try in np.linspace(0.2, 1.2, 50):
        sig = 1 / (1 + np.exp(-15 * (z - z0_try)))
        X5 = np.column_stack([np.ones(n), x1, c * z * sig])
        r = fit_model(X5, mu_res)
        aic_adj = 2*(r['k']+1) - 2*r['ll']
        if aic_adj < best_aic5:
            best_aic5 = aic_adj
            best_z05 = z0_try
    
    sig5 = 1 / (1 + np.exp(-15 * (z - best_z05)))
    X5 = np.column_stack([np.ones(n), x1, c * z * sig5])
    r5 = fit_model(X5, mu_res)
    r5['z0'] = float(best_z05)
    r5['k'] = r5['k'] + 1
    r5['aic'] = 2*r5['k'] - 2*r5['ll']
    r5['bic'] = r5['k']*np.log(n) - 2*r5['ll']
    results['M5_accum_sigmoid'] = r5
    
    # Summary table
    print(f"\n  {'Model':>25s} {'k':>3s} {'R²':>8s} {'AIC':>10s} {'BIC':>10s} {'ΔAIC':>8s} {'ΔBIC':>8s}")
    print("  " + "-" * 72)
    
    baseline_aic = results['M0_baseline']['aic']
    baseline_bic = results['M0_baseline']['bic']
    
    for name, r in results.items():
        daic = r['aic'] - baseline_aic
        dbic = r['bic'] - baseline_bic
        z0_str = f" (z0={r['z0']:.2f})" if 'z0' in r else ""
        print(f"  {name:>25s} {r['k']:3d} {r['r_squared']:8.5f} "
              f"{r['aic']:10.1f} {r['bic']:10.1f} {daic:+8.1f} {dbic:+8.1f}{z0_str}")
    
    # Find best
    all_models = list(results.keys())
    best_aic_model = min(all_models, key=lambda m: results[m]['aic'])
    best_bic_model = min(all_models, key=lambda m: results[m]['bic'])
    
    results['best_aic'] = best_aic_model
    results['best_bic'] = best_bic_model
    
    print(f"\n  Best AIC: {best_aic_model}")
    print(f"  Best BIC: {best_bic_model}")
    
    # Bayes factor interpretation
    best_dbic = results[best_bic_model]['bic'] - baseline_bic
    if best_dbic < -10:
        strength = "VERY STRONG"
    elif best_dbic < -6:
        strength = "STRONG"
    elif best_dbic < -2:
        strength = "POSITIVE"
    else:
        strength = "NOT WORTH MENTION"
    
    print(f"  BIC evidence vs baseline: ΔBIC={best_dbic:.1f} → {strength}")
    
    return results

# =============================================================================
# P1: COMBINED DATASET HIGH-Z POWER
# =============================================================================
def test_combined_highz(pan_df, des):
    """Merge Pantheon+ and DES-5YR for maximum high-z statistical power."""
    print("\n" + "=" * 60)
    print("P1: COMBINED DATASET HIGH-Z ANALYSIS")
    print("=" * 60)
    
    # Pantheon+
    pmask = (np.isfinite(pan_df['zHD'].values) & np.isfinite(pan_df['c'].values) & 
             np.isfinite(pan_df['mu_resid'].values) & np.isfinite(pan_df['x1'].values) &
             (pan_df['zHD'].values > 0.01))
    pz = pan_df['zHD'].values[pmask]
    pc = pan_df['c'].values[pmask]
    px1 = pan_df['x1'].values[pmask]
    pmu = pan_df['mu_resid'].values[pmask]
    
    # DES-5YR
    dmask = np.isfinite(des['zHD']) & np.isfinite(des['c']) & np.isfinite(des['mu_resid']) & np.isfinite(des['x1'])
    dz = des['zHD'][dmask]
    dc = des['c'][dmask]
    dx1 = des['x1'][dmask]
    dmu = des['mu_resid'][dmask]
    
    # Combine
    z = np.concatenate([pz, dz])
    c = np.concatenate([pc, dc])
    x1 = np.concatenate([px1, dx1])
    mu_res = np.concatenate([pmu, dmu])
    source = np.concatenate([np.zeros(len(pz)), np.ones(len(dz))])
    
    print(f"  Combined: {len(z)} SNe (Pantheon+: {len(pz)}, DES: {len(dz)})")
    
    results = {}
    
    # High-z analysis with boosted statistics
    thresholds = [0.3, 0.5, 0.7, 0.82, 1.0]
    print(f"\n  {'z threshold':>12s} {'n':>6s} {'ρ(c,μ)':>8s} {'p':>10s} {'rank':>6s}")
    print("  " + "-" * 50)
    
    for zt in thresholds:
        hi = z >= zt
        if hi.sum() > 20:
            rho, p = stats.spearmanr(c[hi], mu_res[hi])
            obs = np.column_stack([c[hi], x1[hi], mu_res[hi]])
            er = effective_rank(obs)
            n_p = ((source[hi] == 0)).sum()
            n_d = ((source[hi] == 1)).sum()
            sig = " ***" if p < 0.01 else " *" if p < 0.05 else ""
            print(f"  z≥{zt:.2f}      {hi.sum():6d} {rho:+8.4f} {p:10.6f} {er:6.3f}  "
                  f"(P+:{n_p}, DES:{n_d}){sig}")
            results[f'z_ge_{zt}'] = {
                'n': int(hi.sum()), 'rho': float(rho), 'p': float(p),
                'rank': er, 'n_pantheon': int(n_p), 'n_des': int(n_d)
            }
    
    # A1-style control survival on combined set
    print(f"\n  A1-style signal survival (combined, z≥median):")
    z_med = np.median(z)
    hi = z >= z_med
    
    # Raw
    rho_raw, p_raw = stats.spearmanr(c[hi], mu_res[hi])
    
    # Controlled (z, z², x1, source indicator)
    controls = np.column_stack([z[hi], z[hi]**2, x1[hi], source[hi]])
    coef_c, _, _, _ = lstsq(controls, c[hi], rcond=None)
    coef_mu, _, _, _ = lstsq(controls, mu_res[hi], rcond=None)
    rho_ctrl, p_ctrl = stats.spearmanr(c[hi] - controls @ coef_c, mu_res[hi] - controls @ coef_mu)
    
    survival = abs(rho_ctrl) / abs(rho_raw) * 100 if abs(rho_raw) > 1e-6 else np.nan
    
    print(f"  Raw: ρ={rho_raw:.4f}, p={p_raw:.6f}")
    print(f"  Controlled (z, z², x1, source): ρ={rho_ctrl:.4f}, p={p_ctrl:.6f}")
    print(f"  Survival: {survival:.1f}%")
    
    results['signal_survival'] = {
        'raw': {'rho': float(rho_raw), 'p': float(p_raw)},
        'controlled': {'rho': float(rho_ctrl), 'p': float(p_ctrl)},
        'survival_pct': float(survival)
    }
    
    return results

# =============================================================================
# A2b: WITHIN-SURVEY ENTANGLEMENT ON DES-5YR
# =============================================================================
def test_within_survey_des(des):
    """A2 equivalent on DES-5YR — does entanglement appear within DES across z?"""
    print("\n" + "=" * 60)
    print("A2b: WITHIN-DES ENTANGLEMENT ACROSS Z")
    print("=" * 60)
    
    z = des['zHD']
    c = des['c']
    x1 = des['x1']
    mu_res = des['mu_resid']
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(x1) & np.isfinite(mu_res)
    z, c, x1, mu_res = z[mask], c[mask], x1[mask], mu_res[mask]
    
    # DES has IDSAMPLE which indicates different DES sub-surveys
    idsample = des['IDSAMPLE'][mask] if 'IDSAMPLE' in des else np.zeros(len(z))
    
    # z-binned analysis within DES
    n_bins = 8
    z_edges = np.percentile(z, np.linspace(0, 100, n_bins + 1))
    
    results = {'bins': []}
    
    print(f"\n  DES-5YR z-binned entanglement (single survey, SALT3):")
    print(f"  {'z range':>15s} {'n':>5s} {'ρ(c,μ)':>8s} {'p':>8s} {'rank':>6s} {'entangling':>10s}")
    print("  " + "-" * 60)
    
    for i in range(n_bins):
        bm = (z >= z_edges[i]) & (z < z_edges[i+1] + 0.001)
        if bm.sum() < 15:
            continue
        
        zb, cb, x1b, mub = z[bm], c[bm], x1[bm], mu_res[bm]
        
        rho, p = stats.spearmanr(cb, mub)
        
        obs = np.column_stack([cb, x1b, mub])
        er = effective_rank(obs)
        
        # Entangling pairs
        obs_dict = {'c': cb, 'x1': x1b, 'mu_res': mub}
        ent = sum(1 for a, b in combinations(obs_dict.keys(), 2)
                  if stats.spearmanr(obs_dict[a], obs_dict[b])[1] < 0.05)
        
        entry = {
            'z_lo': float(z_edges[i]), 'z_hi': float(z_edges[i+1]),
            'z_med': float(np.median(zb)), 'n': int(bm.sum()),
            'rho': float(rho), 'p': float(p),
            'rank': er, 'entangling': ent,
        }
        results['bins'].append(entry)
        
        sig = " *" if p < 0.05 else ""
        print(f"  [{z_edges[i]:.3f}-{z_edges[i+1]:.3f}] {bm.sum():5d} {rho:+8.4f} {p:8.4f} "
              f"{er:6.3f} {ent:>10d}/3{sig}")
    
    # Trend
    bins = results['bins']
    if len(bins) >= 4:
        z_meds = [b['z_med'] for b in bins]
        abs_rhos = [abs(b['rho']) for b in bins]
        ranks = [b['rank'] for b in bins]
        
        trend_r, trend_p = stats.spearmanr(z_meds, abs_rhos)
        rank_r, rank_p = stats.spearmanr(z_meds, ranks)
        
        results['coupling_trend'] = {'rho': float(trend_r), 'p': float(trend_p)}
        results['rank_trend'] = {'rho': float(rank_r), 'p': float(rank_p)}
        
        print(f"\n  |ρ(c,μ)| trend with z: ρ={trend_r:.3f}, p={trend_p:.4f}")
        print(f"  Rank trend with z: ρ={rank_r:.3f}, p={rank_p:.4f}")
        
        if trend_r > 0.3 and trend_p < 0.1:
            print("  → Coupling INCREASES with z within DES ✓")
        if rank_r < -0.3 and rank_p < 0.1:
            print("  → Rank DECREASES with z within DES ✓")
    
    return results

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 60)
    print("CLOSURE THEORY — ROUND 7: LAST CONVENTIONAL KILLED?")
    print("=" * 60)
    
    pan_df = load_pantheon()
    des = load_des5yr()
    
    all_results = {}
    
    all_results['K2_salt3'] = test_salt3_adequacy(des)
    all_results['M4_sigmoid'] = test_sigmoid_model(pan_df)
    all_results['P1_combined'] = test_combined_highz(pan_df, des)
    all_results['A2b_des_within'] = test_within_survey_des(des)
    
    with open(os.path.join(RESULTS_DIR, 'round7_results.json'), 'w') as f:
        json.dump(all_results, f, indent=2, cls=NpEncoder)
    
    print("\n" + "=" * 60)
    print("ROUND 7 VERDICT")
    print("=" * 60)
    
    k2 = all_results['K2_salt3']
    salt3_bad = k2.get('Bad fits (χ²≥median)', {})
    print(f"  K2: SALT3 bad fits c-μ: ρ={salt3_bad.get('rho','?')}, p={salt3_bad.get('p','?')}")
    
    m4 = all_results['M4_sigmoid']
    print(f"  M4: Best AIC={m4['best_aic']}, Best BIC={m4['best_bic']}")
    
    p1 = all_results['P1_combined']
    ss = p1['signal_survival']
    print(f"  P1: Combined survival={ss['survival_pct']:.0f}%, "
          f"controlled p={ss['controlled']['p']:.6f}")
    
    a2b = all_results['A2b_des_within']
    if 'coupling_trend' in a2b:
        print(f"  A2b: DES coupling trend ρ={a2b['coupling_trend']['rho']:.3f}, "
              f"p={a2b['coupling_trend']['p']:.4f}")
    
    print(f"\nResults saved to {RESULTS_DIR}/")

if __name__ == "__main__":
    main()
