#!/usr/bin/env python3
"""
Closure Theory — Round 9: Model-Independent Color Tests
=========================================================
Bypass SALT2/SALT3 color parameter entirely. Test closure predictions
using raw photometric observables where possible.

Tests:
  CI1: Multi-band magnitude spread as color proxy — does the SPREAD of 
       magnitudes across bands correlate with distance at high-z?
       (If closure is real, frequency-dependent observables lose separability,
       so the relationship between bands should change at high-z)
  CI2: Hubble residual variance by z-bin — closure predicts increasing 
       scatter at high-z specifically in color-dependent terms
  CI3: SALT2 χ²/dof vs z — model-independent: does light-curve fit quality
       degrade with z in a step-like fashion? (No SALT color param used)
  CI4: Color-independent distance test — use ONLY x1 (stretch) for 
       standardization. If closure is real, x1-only residuals should NOT 
       show the entanglement pattern (frequency fingerprint confirmation)
  CI5: Raw c distribution shape by z — does the color distribution itself
       change shape (not just mean) at high-z? Closure predicts loss of
       diversity/information, which could manifest as distribution narrowing
       or skew change.

Data: Pantheon+ SH0ES
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

BASE = Path(__file__).parent
RESULTS_DIR = BASE / "results_color_independent"
RESULTS_DIR.mkdir(exist_ok=True)


def load_pantheon():
    df = pd.read_csv(BASE / "data" / "pantheon_plus.dat", sep=r'\s+', comment='#')
    for col in ['zHD', 'MU_SH0ES', 'MU_SH0ES_ERR_DIAG', 'c', 'x1', 'mB',
                'FITPROB', 'NDOF', 'SNRMAX1', 'HOST_LOGMASS']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna(subset=['zHD', 'MU_SH0ES', 'c', 'x1'])
    df = df[df['zHD'] > 0.01].copy()
    
    # μ_resid
    from numpy.polynomial import polynomial as P
    log_z = np.log10(df['zHD'].values)
    mu = df['MU_SH0ES'].values
    m = df['zHD'].values < 0.5
    coeffs = P.polyfit(log_z[m], mu[m], 2)
    df['mu_resid'] = mu - P.polyval(log_z, coeffs)
    
    return df


def test_CI1_magnitude_spread(df):
    """CI1: Does SALT2 mB-μ relationship change character at high-z?
    
    Use mB (observed peak mag) minus the Tripp-standardized μ.
    The residual mB - μ_SH0ES isolates color+intrinsic terms.
    At high-z, if color entangles, this residual should couple differently.
    """
    if 'mB' not in df.columns or df['mB'].isna().all():
        return {'test': 'CI1: mB-μ residual coupling', 'status': 'mB not available'}
    
    m = np.isfinite(df['mB'].values) & np.isfinite(df['MU_SH0ES'].values)
    mb = df['mB'].values[m]
    mu = df['MU_SH0ES'].values[m]
    z = df['zHD'].values[m]
    
    # mB - μ isolates the color + magnitude offset terms
    delta = mb - mu
    
    # Correlation of delta with z
    rho_all, p_all = stats.spearmanr(z, delta)
    
    z_thresh = 0.82
    lo = z < z_thresh
    hi = z >= z_thresh
    
    rho_lo, p_lo = stats.spearmanr(z[lo], delta[lo]) if lo.sum() > 10 else (None, None)
    rho_hi, p_hi = stats.spearmanr(z[hi], delta[hi]) if hi.sum() > 10 else (None, None)
    
    return {
        'test': 'CI1: mB−μ residual vs z (color+intrinsic proxy)',
        'N': int(m.sum()),
        'N_low_z': int(lo.sum()),
        'N_high_z': int(hi.sum()),
        'rho_all': round(rho_all, 4),
        'p_all': round(p_all, 6),
        'rho_low_z': round(rho_lo, 4) if rho_lo else None,
        'p_low_z': round(p_lo, 6) if p_lo else None,
        'rho_high_z': round(rho_hi, 4) if rho_hi else None,
        'p_high_z': round(p_hi, 6) if p_hi else None,
    }


def test_CI2_residual_variance(df):
    """CI2: Hubble residual variance by z-bin.
    
    Closure predicts scatter should increase non-linearly at high-z
    as observables entangle and standardization degrades.
    """
    z = df['zHD'].values
    mu_r = df['mu_resid'].values
    
    # Create z bins
    bin_edges = [0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.82, 1.0, 2.5]
    bins = []
    for i in range(len(bin_edges) - 1):
        mask = (z >= bin_edges[i]) & (z < bin_edges[i+1])
        if mask.sum() > 5:
            bins.append({
                'z_range': f'{bin_edges[i]:.2f}-{bin_edges[i+1]:.2f}',
                'N': int(mask.sum()),
                'std_mu_resid': round(np.std(mu_r[mask]), 4),
                'mad_mu_resid': round(np.median(np.abs(mu_r[mask] - np.median(mu_r[mask]))), 4),
                'iqr_mu_resid': round(np.subtract(*np.percentile(mu_r[mask], [75, 25])), 4),
            })
    
    # Test: is scatter increase nonlinear? Compare linear vs step model
    z_mids = []
    scatters = []
    for b in bins:
        z_lo, z_hi = map(float, b['z_range'].split('-'))
        z_mids.append((z_lo + z_hi) / 2)
        scatters.append(b['std_mu_resid'])
    
    z_mids = np.array(z_mids)
    scatters = np.array(scatters)
    
    # Linear trend
    slope, intercept, r, p, _ = stats.linregress(z_mids, scatters)
    
    return {
        'test': 'CI2: Hubble residual scatter by z-bin',
        'bins': bins,
        'scatter_vs_z_slope': round(slope, 4),
        'scatter_vs_z_r': round(r, 4),
        'scatter_vs_z_p': round(p, 6),
        'prediction': 'Scatter increases, possibly non-linearly at z>0.82',
    }


def test_CI3_chi2_vs_z(df):
    """CI3: SALT2 fit quality (χ²/dof) vs z — step-like degradation?
    
    This is model-independent in the sense that χ²/dof measures how well
    the light-curve template fits, regardless of what color parameter comes out.
    """
    if 'FITPROB' not in df.columns and 'NDOF' not in df.columns:
        return {'test': 'CI3: χ²/dof vs z', 'status': 'fit quality columns not available'}
    
    # Use FITPROB as proxy (probability of χ²)
    m = np.isfinite(df.get('FITPROB', pd.Series(dtype=float)).values)
    if m.sum() < 50:
        return {'test': 'CI3: χ²/dof vs z', 'status': 'insufficient FITPROB data'}
    
    z = df['zHD'].values[m]
    fitprob = df['FITPROB'].values[m]
    
    rho, p = stats.spearmanr(z, fitprob)
    
    z_thresh = 0.82
    lo = z < z_thresh
    hi = z >= z_thresh
    
    mean_lo = np.mean(fitprob[lo]) if lo.sum() > 0 else None
    mean_hi = np.mean(fitprob[hi]) if hi.sum() > 0 else None
    
    # Mann-Whitney for distribution difference
    if hi.sum() > 5 and lo.sum() > 5:
        U, p_mw = stats.mannwhitneyu(fitprob[lo], fitprob[hi], alternative='two-sided')
    else:
        p_mw = None
    
    return {
        'test': 'CI3: FITPROB vs z (fit quality degradation)',
        'N': int(m.sum()),
        'rho': round(rho, 4),
        'p': round(p, 6),
        'mean_fitprob_low_z': round(mean_lo, 4) if mean_lo else None,
        'mean_fitprob_high_z': round(mean_hi, 4) if mean_hi else None,
        'mannwhitney_p': round(p_mw, 6) if p_mw else None,
        'prediction': 'Fit quality degrades at high-z (lower FITPROB)',
    }


def test_CI4_x1_only_standardization(df):
    """CI4: Standardize with x1 ONLY (no color). 
    
    If closure only affects frequency-dependent observables, then x1-only
    standardized residuals should NOT show the c-μ entanglement pattern.
    This is the frequency fingerprint test without using SALT color at all.
    """
    z = df['zHD'].values
    mu = df['MU_SH0ES'].values
    x1 = df['x1'].values
    c = df['c'].values
    
    # Fit x1-only standardization at low-z
    from numpy.polynomial import polynomial as P
    log_z = np.log10(z)
    
    # First get cosmological μ(z) — same 2nd order fit
    m_fit = z < 0.5
    coeffs_z = P.polyfit(log_z[m_fit], mu[m_fit], 2)
    mu_cosmo = P.polyval(log_z, coeffs_z)
    
    # Now fit: μ_obs - μ_cosmo = α*x1 + const (low-z only)
    delta_mu = mu - mu_cosmo
    A = np.column_stack([x1[m_fit], np.ones(m_fit.sum())])
    from numpy.linalg import lstsq
    params = lstsq(A, delta_mu[m_fit], rcond=None)[0]
    alpha_x1 = params[0]
    
    # x1-only residual
    mu_x1_only = mu - alpha_x1 * x1
    mu_x1_resid = mu_x1_only - P.polyval(log_z, P.polyfit(log_z[m_fit], mu_x1_only[m_fit], 2))
    
    # Now test: does c correlate with x1-only residuals at high-z?
    z_thresh = 0.82
    hi = z >= z_thresh
    lo = z < z_thresh
    
    rho_lo, p_lo = stats.spearmanr(c[lo], mu_x1_resid[lo])
    rho_hi, p_hi = stats.spearmanr(c[hi], mu_x1_resid[hi]) if hi.sum() > 10 else (None, None)
    
    # Compare with standard (c+x1) residuals
    c_resid_corr_lo = stats.spearmanr(c[lo], df['mu_resid'].values[lo])
    c_resid_corr_hi = stats.spearmanr(c[hi], df['mu_resid'].values[hi]) if hi.sum() > 10 else (None, None)
    
    return {
        'test': 'CI4: x1-only standardization (bypass SALT color)',
        'alpha_x1_fit': round(alpha_x1, 4),
        'N_low_z': int(lo.sum()),
        'N_high_z': int(hi.sum()),
        'c_vs_x1resid_rho_low_z': round(rho_lo, 4),
        'c_vs_x1resid_p_low_z': round(p_lo, 6),
        'c_vs_x1resid_rho_high_z': round(rho_hi, 4) if rho_hi else None,
        'c_vs_x1resid_p_high_z': round(p_hi, 6) if p_hi else None,
        'c_vs_standard_resid_rho_high_z': round(c_resid_corr_hi[0], 4) if c_resid_corr_hi else None,
        'prediction': 'c should correlate with x1-only resid at high-z (color info NOT removed)',
        'interpretation': 'If c-coupling appears even without SALT color correction, it is real physics not model artifact',
    }


def test_CI5_color_distribution_shape(df):
    """CI5: Does the color distribution change shape at high-z?
    
    Closure predicts loss of information diversity. If frequency-dependent 
    observables entangle, the color distribution might narrow, shift, or 
    change skewness at high-z.
    """
    z = df['zHD'].values
    c = df['c'].values
    
    z_thresh = 0.82
    lo = (z < z_thresh) & (z > 0.3)  # avoid very low-z selection effects
    hi = z >= z_thresh
    
    c_lo = c[lo]
    c_hi = c[hi]
    
    # Distribution statistics
    result = {
        'test': 'CI5: Color distribution shape by z regime',
        'N_low_z': int(lo.sum()),
        'N_high_z': int(hi.sum()),
        'mean_c_low_z': round(np.mean(c_lo), 4),
        'mean_c_high_z': round(np.mean(c_hi), 4),
        'std_c_low_z': round(np.std(c_lo), 4),
        'std_c_high_z': round(np.std(c_hi), 4),
        'skew_c_low_z': round(float(stats.skew(c_lo)), 4),
        'skew_c_high_z': round(float(stats.skew(c_hi)), 4),
        'kurtosis_c_low_z': round(float(stats.kurtosis(c_lo)), 4),
        'kurtosis_c_high_z': round(float(stats.kurtosis(c_hi)), 4),
    }
    
    # KS test for distribution difference
    ks_stat, ks_p = stats.ks_2samp(c_lo, c_hi)
    result['ks_stat'] = round(ks_stat, 4)
    result['ks_p'] = round(ks_p, 6)
    
    # Anderson-Darling k-sample
    ad_stat, _, ad_p = stats.anderson_ksamp([c_lo, c_hi])
    result['anderson_darling_stat'] = round(ad_stat, 4)
    result['anderson_darling_p'] = round(ad_p, 6)
    
    # Levene test for variance equality
    lev_stat, lev_p = stats.levene(c_lo, c_hi)
    result['levene_stat'] = round(lev_stat, 4)
    result['levene_p'] = round(lev_p, 6)
    
    result['prediction'] = 'Distribution change at high-z (information loss → reduced diversity)'
    
    return result


def main():
    print("=" * 70)
    print("CLOSURE THEORY — ROUND 9: MODEL-INDEPENDENT COLOR TESTS")
    print("=" * 70)
    
    print("\nLoading Pantheon+ data...")
    df = load_pantheon()
    df = df.reset_index(drop=True)
    print(f"  → {len(df)} SNe after cuts")
    
    all_results = []
    
    tests = [
        ("CI1", test_CI1_magnitude_spread),
        ("CI2", test_CI2_residual_variance),
        ("CI3", test_CI3_chi2_vs_z),
        ("CI4", test_CI4_x1_only_standardization),
        ("CI5", test_CI5_color_distribution_shape),
    ]
    
    for tid, func in tests:
        print(f"\n{'─' * 50}")
        print(f"Running {tid}...")
        result = func(df)
        all_results.append(result)
        for k, v in result.items():
            if k == 'bins':
                print(f"  {k}:")
                for b in v:
                    print(f"    {b}")
            else:
                print(f"  {k}: {v}")
    
    with open(RESULTS_DIR / "color_independent_results.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print(f"\n{'=' * 70}")
    print("SUMMARY — ROUND 9: MODEL-INDEPENDENT COLOR")
    print(f"{'=' * 70}")
    for r in all_results:
        print(f"  {r['test']}")
    
    print(f"\nResults saved to {RESULTS_DIR}/")
    print("Done.")


if __name__ == "__main__":
    main()
