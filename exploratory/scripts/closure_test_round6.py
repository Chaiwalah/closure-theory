#!/usr/bin/env python3
"""
Closure Theory — Round 6: Kill the Conventional Explanations
==============================================================
Targeted tests to distinguish "closure / missing latent variable" from
"selection + calibration + model mis-specification."

Tests:
  A2: Within-survey, within-z-bin entanglement (the decisive test)
  V1: VPEC hygiene — does peculiar velocity correction affect high-z signal?
  K1: SALT2 model adequacy — rest-frame UV proxy test
  S6: Selection/truncation test — does the high-z sample narrow in observable space?
  A3: Collider bias check — are controls inducing spurious amplification?
"""

import numpy as np
import os
import json
from scipy import stats
from scipy.integrate import quad
from numpy.linalg import lstsq
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

RESULTS_DIR = "results_round6"
os.makedirs(RESULTS_DIR, exist_ok=True)

def load_data():
    import pandas as pd
    df = pd.read_csv("data/pantheon_plus.dat", sep=r'\s+', comment='#')
    
    z = df['zHD'].values
    
    def mu_lcdm(zv, H0=73.04, Om=0.334):
        if not np.isfinite(zv) or zv <= 0: return np.nan
        dc, _ = quad(lambda zp: 1/np.sqrt(Om*(1+zp)**3+(1-Om)), 0, zv)
        return 5*np.log10((1+zv)*dc*2.998e5/H0) + 25
    
    mu_model = np.array([mu_lcdm(zi) for zi in z])
    
    df['mu_resid'] = df['MU_SH0ES'].values - mu_model if 'MU_SH0ES' in df.columns else df['m_b_corr'].values - mu_model
    
    # Color residual
    c = df['c'].values
    mask = np.isfinite(c) & np.isfinite(z)
    poly = np.polyfit(z[mask], c[mask], 1)
    df['c_resid'] = c - np.polyval(poly, z)
    
    # x1 residual
    x1 = df['x1'].values
    mask = np.isfinite(x1) & np.isfinite(z)
    poly = np.polyfit(z[mask], x1[mask], 1)
    df['x1_resid'] = x1 - np.polyval(poly, z)
    
    # chi2/dof
    if 'NDOF' in df.columns:
        ndof = df['NDOF'].values.astype(float)
        ndof[ndof == 0] = np.nan
        df['chi2_dof'] = df['FITCHI2'].values.astype(float) / ndof
    
    return df

def effective_rank(X):
    X_std = (X - X.mean(0)) / (X.std(0) + 1e-10)
    s = np.linalg.svd(X_std, compute_uv=False)
    p = s / s.sum()
    p = p[p > 0]
    return float(np.exp(-np.sum(p * np.log(p))))

# =============================================================================
# A2: WITHIN-SURVEY, WITHIN-Z-BIN ENTANGLEMENT
# =============================================================================
def test_within_survey(df):
    """The decisive test: does entanglement appear within a single survey across z?"""
    print("\n" + "=" * 60)
    print("A2: WITHIN-SURVEY ENTANGLEMENT ACROSS Z")
    print("=" * 60)
    
    z = df['zHD'].values
    c = df['c'].values
    x1 = df['x1'].values
    mu_res = df['mu_resid'].values
    survey = df['IDSURVEY'].values
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(x1) & np.isfinite(mu_res) & (z > 0.01)
    
    # Find surveys with good z coverage
    unique_surveys = np.unique(survey[mask])
    
    results = {}
    
    print(f"\n  Survey coverage:")
    for sid in unique_surveys:
        s_mask = mask & (survey == sid)
        n = s_mask.sum()
        if n > 30:
            z_range = [float(z[s_mask].min()), float(z[s_mask].max())]
            print(f"    Survey {int(sid)}: n={n}, z=[{z_range[0]:.3f}, {z_range[1]:.3f}]")
    
    # For each survey with enough range, split into z-bins and track entanglement
    print(f"\n  Within-survey entanglement profiles:")
    
    for sid in unique_surveys:
        s_mask = mask & (survey == sid)
        n = s_mask.sum()
        z_range = z[s_mask].max() - z[s_mask].min()
        
        if n < 80 or z_range < 0.2:
            continue
        
        zs = z[s_mask]
        cs = c[s_mask]
        x1s = x1[s_mask]
        mu_s = mu_res[s_mask]
        
        # Create z-bins within this survey
        n_bins = min(5, max(3, n // 40))
        z_edges = np.percentile(zs, np.linspace(0, 100, n_bins + 1))
        
        print(f"\n    Survey {int(sid)} (n={n}, z=[{zs.min():.3f}-{zs.max():.3f}]):")
        
        survey_results = {'bins': []}
        
        for i in range(n_bins):
            bin_mask = (zs >= z_edges[i]) & (zs < z_edges[i+1] + 0.001)
            if bin_mask.sum() < 15:
                continue
            
            zb = zs[bin_mask]
            cb = cs[bin_mask]
            x1b = x1s[bin_mask]
            mub = mu_s[bin_mask]
            
            # Color-distance coupling
            rho_c, p_c = stats.spearmanr(cb, mub)
            
            # Effective rank
            obs = np.column_stack([cb, x1b, mub])
            valid = np.all(np.isfinite(obs), axis=1)
            if valid.sum() > 10:
                er = effective_rank(obs[valid])
            else:
                er = np.nan
            
            # Pairwise entanglement count
            obs_dict = {'c': cb, 'x1': x1b, 'mu_res': mub}
            entangling = 0
            pairs = list(combinations(obs_dict.keys(), 2))
            for a, b in pairs:
                r, p = stats.spearmanr(obs_dict[a], obs_dict[b])
                if p < 0.05:
                    entangling += 1
            
            entry = {
                'z_lo': float(z_edges[i]), 'z_hi': float(z_edges[i+1]),
                'z_med': float(np.median(zb)), 'n': int(bin_mask.sum()),
                'rho_c_mu': float(rho_c), 'p_c_mu': float(p_c),
                'eff_rank': er,
                'entangling_pairs': entangling,
            }
            survey_results['bins'].append(entry)
            
            sig = " *" if p_c < 0.05 else ""
            print(f"      z=[{z_edges[i]:.3f}-{z_edges[i+1]:.3f}] (n={bin_mask.sum():3d}): "
                  f"ρ(c,μ)={rho_c:+.4f} p={p_c:.3f}, rank={er:.3f}, "
                  f"entangling={entangling}/3{sig}")
        
        # Trend within this survey
        bins = survey_results['bins']
        if len(bins) >= 3:
            z_meds = [b['z_med'] for b in bins]
            rhos = [abs(b['rho_c_mu']) for b in bins]
            ranks = [b['eff_rank'] for b in bins if np.isfinite(b['eff_rank'])]
            
            trend_r, trend_p = stats.spearmanr(z_meds, rhos)
            survey_results['coupling_trend'] = {'rho': float(trend_r), 'p': float(trend_p)}
            
            if len(ranks) >= 3:
                z_for_rank = [b['z_med'] for b in bins if np.isfinite(b['eff_rank'])]
                rank_r, rank_p = stats.spearmanr(z_for_rank, ranks)
                survey_results['rank_trend'] = {'rho': float(rank_r), 'p': float(rank_p)}
                print(f"      Coupling trend: ρ={trend_r:+.3f} (p={trend_p:.3f})")
                print(f"      Rank trend: ρ={rank_r:+.3f} (p={rank_p:.3f})")
            else:
                print(f"      Coupling trend: ρ={trend_r:+.3f} (p={trend_p:.3f})")
        
        results[f'survey_{int(sid)}'] = survey_results
    
    return results

# =============================================================================
# V1: VPEC HYGIENE TEST
# =============================================================================
def test_vpec(df):
    """Does peculiar velocity affect the high-z signal?"""
    print("\n" + "=" * 60)
    print("V1: PECULIAR VELOCITY HYGIENE")
    print("=" * 60)
    
    z = df['zHD'].values
    c = df['c'].values
    mu_res = df['mu_resid'].values
    
    # VPEC columns
    vpec = df['VPEC'].values if 'VPEC' in df.columns else None
    vpec_err = df['VPECERR'].values if 'VPECERR' in df.columns else None
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(mu_res) & (z > 0.01)
    
    results = {}
    
    if vpec is not None:
        mask_vpec = mask & np.isfinite(vpec)
        print(f"  SNe with VPEC data: {mask_vpec.sum()}")
        
        # 1. Does VPEC correlate with c or mu_res?
        rho_vc, p_vc = stats.spearmanr(vpec[mask_vpec], c[mask_vpec])
        rho_vm, p_vm = stats.spearmanr(vpec[mask_vpec], mu_res[mask_vpec])
        print(f"  VPEC vs c: ρ={rho_vc:.4f}, p={p_vc:.4f}")
        print(f"  VPEC vs μ_resid: ρ={rho_vm:.4f}, p={p_vm:.4f}")
        results['vpec_c'] = {'rho': float(rho_vc), 'p': float(p_vc)}
        results['vpec_mu'] = {'rho': float(rho_vm), 'p': float(p_vm)}
        
        # 2. High-z signal with VPEC as additional control
        z_med = np.median(z[mask_vpec])
        hi = mask_vpec & (z >= z_med)
        
        # Without VPEC control
        rho_raw, p_raw = stats.spearmanr(c[hi], mu_res[hi])
        
        # With VPEC control
        controls = np.column_stack([z[hi], z[hi]**2, vpec[hi]])
        coef_c, _, _, _ = lstsq(controls, c[hi], rcond=None)
        coef_mu, _, _, _ = lstsq(controls, mu_res[hi], rcond=None)
        c_ctrl = c[hi] - controls @ coef_c
        mu_ctrl = mu_res[hi] - controls @ coef_mu
        rho_ctrl, p_ctrl = stats.spearmanr(c_ctrl, mu_ctrl)
        
        print(f"\n  High-z c-μ coupling:")
        print(f"    Without VPEC control: ρ={rho_raw:.4f}, p={p_raw:.4f}")
        print(f"    With VPEC control:    ρ={rho_ctrl:.4f}, p={p_ctrl:.4f}")
        
        results['high_z_without_vpec'] = {'rho': float(rho_raw), 'p': float(p_raw)}
        results['high_z_with_vpec'] = {'rho': float(rho_ctrl), 'p': float(p_ctrl)}
        
        survival = abs(rho_ctrl) / abs(rho_raw) * 100 if abs(rho_raw) > 0 else np.nan
        print(f"    Survival: {survival:.1f}%")
        results['survival_pct'] = float(survival)
        
        # 3. VPEC contribution by z-bin
        print(f"\n  VPEC magnitude by z-bin:")
        z_bins = np.percentile(z[mask_vpec], [0, 25, 50, 75, 100])
        for i in range(4):
            bm = mask_vpec & (z >= z_bins[i]) & (z < z_bins[i+1] + 0.001)
            if bm.sum() > 10:
                v = vpec[bm]
                print(f"    z=[{z_bins[i]:.2f}-{z_bins[i+1]:.2f}] (n={bm.sum()}): "
                      f"|VPEC| mean={np.abs(v).mean():.1f}, std={v.std():.1f}")
    else:
        print("  No VPEC column found — using zHD vs zHEL difference as proxy")
        
        if 'zHEL' in df.columns:
            z_hel = df['zHEL'].values
            # VPEC proxy: difference between CMB-corrected and heliocentric z
            vpec_proxy = (z - z_hel) * 2.998e5  # km/s
            mask_vp = mask & np.isfinite(vpec_proxy)
            
            rho_vc, p_vc = stats.spearmanr(vpec_proxy[mask_vp], c[mask_vp])
            rho_vm, p_vm = stats.spearmanr(vpec_proxy[mask_vp], mu_res[mask_vp])
            print(f"  VPEC_proxy vs c: ρ={rho_vc:.4f}, p={p_vc:.4f}")
            print(f"  VPEC_proxy vs μ_resid: ρ={rho_vm:.4f}, p={p_vm:.4f}")
            results['vpec_proxy_c'] = {'rho': float(rho_vc), 'p': float(p_vc)}
            results['vpec_proxy_mu'] = {'rho': float(rho_vm), 'p': float(p_vm)}
    
    return results

# =============================================================================
# K1: SALT2 MODEL ADEQUACY — REST-FRAME UV PROXY
# =============================================================================
def test_salt2_adequacy(df):
    """Test whether SALT2 rest-frame UV handling creates the signal."""
    print("\n" + "=" * 60)
    print("K1: SALT2 MODEL ADEQUACY (REST-FRAME UV PROXY)")
    print("=" * 60)
    
    z = df['zHD'].values
    c = df['c'].values
    x1 = df['x1'].values
    mu_res = df['mu_resid'].values
    chi2_dof = df['chi2_dof'].values if 'chi2_dof' in df.columns else None
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(x1) & np.isfinite(mu_res) & (z > 0.01)
    
    results = {}
    
    # At higher z, observed-frame optical samples rest-frame UV
    # SALT2 is less constrained in UV → fits get worse → chi2 increases
    # If closure signal is just SALT2 UV mis-specification:
    #   - chi2_dof should correlate with z
    #   - chi2_dof should correlate with c-mu_resid coupling strength
    #   - removing high-chi2 fits should kill the signal
    
    if chi2_dof is not None:
        mask_chi = mask & np.isfinite(chi2_dof)
        z_m, c_m, mu_m, chi_m = z[mask_chi], c[mask_chi], mu_res[mask_chi], chi2_dof[mask_chi]
        
        # 1. chi2 vs z trend
        rho_cz, p_cz = stats.spearmanr(z_m, chi_m)
        print(f"  χ²/dof vs z: ρ={rho_cz:.4f}, p={p_cz:.4f}")
        results['chi2_z_trend'] = {'rho': float(rho_cz), 'p': float(p_cz)}
        
        # 2. Signal in good fits vs bad fits
        chi_med = np.median(chi_m)
        z_med = np.median(z_m)
        
        for label, chi_mask in [("Good fits (χ²/dof<median)", chi_m < chi_med),
                                 ("Bad fits (χ²/dof≥median)", chi_m >= chi_med)]:
            hi_z = chi_mask & (z_m >= z_med)
            if hi_z.sum() > 30:
                rho, p = stats.spearmanr(c_m[hi_z], mu_m[hi_z])
                print(f"  {label}, high-z: ρ={rho:.4f}, p={p:.4f} (n={hi_z.sum()})")
                results[label] = {'rho': float(rho), 'p': float(p), 'n': int(hi_z.sum())}
        
        # 3. Signal with chi2_dof as control
        hi = mask_chi & (z >= z_med)
        controls = np.column_stack([z_m[z_m >= z_med], chi_m[z_m >= z_med]])
        c_hi = c_m[z_m >= z_med]
        mu_hi = mu_m[z_m >= z_med]
        
        coef_c, _, _, _ = lstsq(controls, c_hi, rcond=None)
        coef_mu, _, _, _ = lstsq(controls, mu_hi, rcond=None)
        rho_ctrl, p_ctrl = stats.spearmanr(c_hi - controls @ coef_c, mu_hi - controls @ coef_mu)
        
        rho_raw, p_raw = stats.spearmanr(c_hi, mu_hi)
        print(f"\n  High-z c-μ without χ² control: ρ={rho_raw:.4f}, p={p_raw:.4f}")
        print(f"  High-z c-μ with χ² control:    ρ={rho_ctrl:.4f}, p={p_ctrl:.4f}")
        
        survival = abs(rho_ctrl) / abs(rho_raw) * 100 if abs(rho_raw) > 0 else np.nan
        print(f"  Survival: {survival:.1f}%")
        results['chi2_controlled'] = {'rho': float(rho_ctrl), 'p': float(p_ctrl), 'survival': float(survival)}
    
    # 4. Rest-frame wavelength coverage proxy
    # At z=0, observer g/r/i = rest g/r/i
    # At z=0.5, observer g = rest ~UV, r = rest ~B
    # At z=1.0, observer g = rest ~far-UV
    # SALT2 UV uncertainty increases → more model-dependent
    
    # Use FITCHI2/NDOF increase rate as UV-sensitivity proxy
    print(f"\n  Rest-frame UV exposure by z:")
    z_edges = [0.01, 0.1, 0.3, 0.5, 0.7, 1.0, 2.5]
    for i in range(len(z_edges)-1):
        bm = mask & (z >= z_edges[i]) & (z < z_edges[i+1])
        if chi2_dof is not None:
            bm_chi = bm & np.isfinite(chi2_dof)
            if bm_chi.sum() > 10:
                rest_uv_frac = min(1.0, z_edges[i] / 0.5)  # rough proxy
                r, p = stats.spearmanr(c[bm_chi], mu_res[bm_chi])
                print(f"    z=[{z_edges[i]:.1f}-{z_edges[i+1]:.1f}] (n={bm_chi.sum()}): "
                      f"ρ(c,μ)={r:+.4f}, <χ²/dof>={chi2_dof[bm_chi].mean():.2f}, UV_frac~{rest_uv_frac:.1f}")
    
    return results

# =============================================================================
# S6: SELECTION/TRUNCATION TEST
# =============================================================================
def test_selection_truncation(df):
    """Does the high-z sample narrow in observable space? (Collider/selection bias check)"""
    print("\n" + "=" * 60)
    print("S6: SELECTION & TRUNCATION BIAS CHECK")
    print("=" * 60)
    
    z = df['zHD'].values
    c = df['c'].values
    x1 = df['x1'].values
    mu_res = df['mu_resid'].values
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(x1) & np.isfinite(mu_res) & (z > 0.01)
    z, c, x1, mu_res = z[mask], c[mask], x1[mask], mu_res[mask]
    
    results = {}
    
    # 1. Observable dispersion by z-bin
    print(f"\n  Observable dispersion by z-bin:")
    z_edges = np.percentile(z, np.linspace(0, 100, 7))
    
    print(f"    {'z range':>15s} {'n':>5s} {'σ(c)':>8s} {'σ(x1)':>8s} {'σ(μ_res)':>8s} {'range(c)':>8s}")
    print("    " + "-" * 55)
    
    bin_data = []
    for i in range(len(z_edges)-1):
        bm = (z >= z_edges[i]) & (z < z_edges[i+1] + 0.001)
        if bm.sum() > 10:
            entry = {
                'z_lo': float(z_edges[i]), 'z_hi': float(z_edges[i+1]),
                'n': int(bm.sum()),
                'std_c': float(c[bm].std()),
                'std_x1': float(x1[bm].std()),
                'std_mu': float(mu_res[bm].std()),
                'range_c': float(c[bm].max() - c[bm].min()),
            }
            bin_data.append(entry)
            print(f"    [{z_edges[i]:.2f}-{z_edges[i+1]:.2f}] {bm.sum():5d} "
                  f"{c[bm].std():8.4f} {x1[bm].std():8.4f} {mu_res[bm].std():8.4f} "
                  f"{c[bm].max()-c[bm].min():8.4f}")
    
    results['dispersion_profile'] = bin_data
    
    # Does σ(c) decrease with z? (selection narrowing)
    z_meds = [(b['z_lo']+b['z_hi'])/2 for b in bin_data]
    std_cs = [b['std_c'] for b in bin_data]
    std_x1s = [b['std_x1'] for b in bin_data]
    
    rho_c_narrow, p_c_narrow = stats.spearmanr(z_meds, std_cs)
    rho_x1_narrow, p_x1_narrow = stats.spearmanr(z_meds, std_x1s)
    
    results['c_narrowing'] = {'rho': float(rho_c_narrow), 'p': float(p_c_narrow)}
    results['x1_narrowing'] = {'rho': float(rho_x1_narrow), 'p': float(p_x1_narrow)}
    
    print(f"\n  σ(c) vs z trend: ρ={rho_c_narrow:.4f}, p={p_c_narrow:.4f}")
    print(f"  σ(x1) vs z trend: ρ={rho_x1_narrow:.4f}, p={p_x1_narrow:.4f}")
    
    if rho_c_narrow < -0.5 and p_c_narrow < 0.05:
        print("  ⚠️ Color distribution NARROWS at high-z — selection truncation possible")
    else:
        print("  ✓ No significant narrowing of color distribution")
    
    # 2. Correlation STRUCTURE change (not just strength)
    # If selection creates coupling, the c-μ relationship should be LINEAR
    # If closure creates it, it should be more complex (threshold-like)
    print(f"\n  Linearity test for high-z c-μ coupling:")
    z_med = np.median(z)
    hi = z >= z_med
    
    # Fit linear vs quadratic
    c_hi, mu_hi = c[hi], mu_res[hi]
    
    from numpy.polynomial import polynomial as P
    # Linear
    lin_coef = np.polyfit(c_hi, mu_hi, 1)
    lin_pred = np.polyval(lin_coef, c_hi)
    ss_lin = np.sum((mu_hi - lin_pred)**2)
    
    # Quadratic
    quad_coef = np.polyfit(c_hi, mu_hi, 2)
    quad_pred = np.polyval(quad_coef, c_hi)
    ss_quad = np.sum((mu_hi - quad_pred)**2)
    
    n_hi = len(c_hi)
    # F-test for quadratic improvement
    f_stat = ((ss_lin - ss_quad) / 1) / (ss_quad / (n_hi - 3))
    p_quad = 1 - stats.f.cdf(f_stat, 1, n_hi - 3)
    
    results['linearity'] = {
        'ss_linear': float(ss_lin), 'ss_quadratic': float(ss_quad),
        'f_stat': float(f_stat), 'p_quadratic': float(p_quad)
    }
    print(f"  Linear SS={ss_lin:.4f}, Quadratic SS={ss_quad:.4f}")
    print(f"  F-test for nonlinearity: F={f_stat:.3f}, p={p_quad:.4f}")
    
    return results

# =============================================================================
# A3: COLLIDER BIAS CHECK
# =============================================================================
def test_collider_bias(df):
    """Does conditioning on controls CREATE the amplification (collider bias)?"""
    print("\n" + "=" * 60)
    print("A3: COLLIDER BIAS CHECK")
    print("=" * 60)
    
    z = df['zHD'].values
    c = df['c'].values
    x1 = df['x1'].values
    mu_res = df['mu_resid'].values
    host_mass = df['HOST_LOGMASS'].values if 'HOST_LOGMASS' in df.columns else np.full(len(z), 10.0)
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(x1) & np.isfinite(mu_res) & (z > 0.01)
    z, c, x1, mu_res = z[mask], c[mask], x1[mask], mu_res[mask]
    hm = host_mass[mask].copy()
    hm[~np.isfinite(hm)] = np.nanmedian(hm)
    
    z_med = np.median(z)
    hi = z >= z_med
    
    results = {}
    
    # The concern: if x1 and host_mass are EFFECTS of both c and mu_res,
    # conditioning on them opens a collider path and amplifies correlation
    
    # Test: Add controls ONE AT A TIME and track how signal changes
    print(f"\n  Sequential control addition (high-z, n={hi.sum()}):")
    
    control_sets = [
        ("Raw (no controls)", []),
        ("+ z, z²", [z[hi], z[hi]**2]),
        ("+ z, z², x1", [z[hi], z[hi]**2, x1[hi]]),
        ("+ z, z², x1, host_mass", [z[hi], z[hi]**2, x1[hi], hm[hi]]),
        ("Only x1", [x1[hi]]),
        ("Only host_mass", [hm[hi]]),
    ]
    
    for label, ctrl_list in control_sets:
        if not ctrl_list:
            rho, p = stats.spearmanr(c[hi], mu_res[hi])
        else:
            controls = np.column_stack(ctrl_list)
            coef_c, _, _, _ = lstsq(controls, c[hi], rcond=None)
            coef_mu, _, _, _ = lstsq(controls, mu_res[hi], rcond=None)
            rho, p = stats.spearmanr(c[hi] - controls @ coef_c, mu_res[hi] - controls @ coef_mu)
        
        results[label] = {'rho': float(rho), 'p': float(p)}
        arrow = ""
        if label != "Raw (no controls)":
            raw_rho = results["Raw (no controls)"]['rho']
            if abs(rho) > abs(raw_rho) * 1.1:
                arrow = " ↑ AMPLIFIED"
            elif abs(rho) < abs(raw_rho) * 0.9:
                arrow = " ↓ reduced"
        
        print(f"    {label:>35s}: ρ={rho:+.4f}, p={p:.6f}{arrow}")
    
    # The key check: if EACH control individually weakens the signal,
    # but together they amplify it → classic collider
    raw = abs(results["Raw (no controls)"]['rho'])
    only_x1 = abs(results["Only x1"]['rho'])
    only_hm = abs(results["Only host_mass"]['rho'])
    all_ctrl = abs(results["+ z, z², x1, host_mass"]['rho'])
    
    individual_weaken = only_x1 < raw and only_hm < raw
    combined_amplify = all_ctrl > raw * 1.1
    
    collider_flag = individual_weaken and combined_amplify
    results['collider_pattern'] = bool(collider_flag)
    
    print(f"\n  Individual controls weaken signal: {individual_weaken}")
    print(f"  Combined controls amplify signal: {combined_amplify}")
    print(f"  Collider pattern detected: {collider_flag}")
    
    if collider_flag:
        print("  ⚠️ COLLIDER BIAS — amplification may be artifact of conditioning")
    else:
        print("  ✓ No classic collider pattern — amplification appears genuine")
    
    return results

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 60)
    print("CLOSURE THEORY — ROUND 6: KILL THE CONVENTIONAL")
    print("=" * 60)
    
    df = load_data()
    all_results = {}
    
    all_results['A2_within_survey'] = test_within_survey(df)
    all_results['V1_vpec'] = test_vpec(df)
    all_results['K1_salt2'] = test_salt2_adequacy(df)
    all_results['S6_selection'] = test_selection_truncation(df)
    all_results['A3_collider'] = test_collider_bias(df)
    
    # Save
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.bool_,)): return bool(obj)
            if isinstance(obj, (np.integer,)): return int(obj)
            if isinstance(obj, (np.floating,)): return float(obj)
            if isinstance(obj, np.ndarray): return obj.tolist()
            return super().default(obj)
    
    with open(os.path.join(RESULTS_DIR, 'round6_results.json'), 'w') as f:
        json.dump(all_results, f, indent=2, cls=NpEncoder)
    
    print("\n" + "=" * 60)
    print("ROUND 6 VERDICT")
    print("=" * 60)
    print(f"\nResults saved to {RESULTS_DIR}/")

if __name__ == "__main__":
    main()
