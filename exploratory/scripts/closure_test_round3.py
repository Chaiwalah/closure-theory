#!/usr/bin/env python3
"""
Closure Theory — Round 3: Reverse-Deduction Tests
===================================================
Tests derived from: "If closure is true, what else must be true?"

Test D1: Rest-Frame Wavelength Stress Test
  - Does the entanglement depend on rest-frame wavelength coverage?
  - If threshold moves after controlling for coverage → conventional (K-correction)
  - If threshold persists → frequency channel effect is real

Test A1: Conditional Independence After Controls
  - Does Corr(Δμ, c) survive after conditioning on survey, z, host mass, SNR, coverage?
  - If yes → model is missing a latent "channel state"
  - If no → conventional systematics explain it

Test C1: Shared Change-Point Detection
  - Do |r(Δμ,c)|, effective rank, and MI all break at the same z*?
  - If same z* → phase transition (closure)
  - If different z* → messy systematics

Test B1: PC1 Direction Stability Across Splits
  - Does PC1 at high-z point the same direction regardless of survey/hemisphere?
  - If stable → universal boundary mode
  - If rotates → survey-coded artifact

Test F1: Predictability Asymmetry
  - Can you predict c from Δμ better than Δμ from c at high-z?
  - Asymmetry → lossy channel with directionality
  - Symmetry → just noise

Author: Humza Hafeez (Closure Theory)
Date: 2026-02-21
"""

import os
import json
import warnings
import subprocess
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit, minimize
from scipy.signal import argrelextrema

warnings.filterwarnings('ignore')

OUTPUT_DIR = Path("results_round3")
DATA_DIR = Path("data")
SEED = 42

def ensure_dir(p):
    p.mkdir(parents=True, exist_ok=True)


# ============================================================
# DATA LOADING (shared with round 1/2)
# ============================================================

def download_pantheon():
    ensure_dir(DATA_DIR)
    dat_file = DATA_DIR / "pantheon_plus.dat"
    if not dat_file.exists():
        urls = [
            "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat",
            "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2BSH0ES.dat",
        ]
        for url in urls:
            result = subprocess.run(["curl", "-fSL", "-o", str(dat_file), url], capture_output=True)
            if result.returncode == 0:
                break
    return dat_file

def load_data():
    dat_file = download_pantheon()
    df = pd.read_csv(dat_file, sep=r'\s+', comment='#')
    col_map = {}
    for col in df.columns:
        if col.lower() in ('decl', 'dec'):
            col_map[col] = 'DEC'
        elif col.lower() == 'ra':
            col_map[col] = 'RA'
    df = df.rename(columns=col_map)

    # Distance residual
    if 'MU_SH0ES' in df.columns and 'MUMODEL' in df.columns:
        df['mu_resid'] = df['MU_SH0ES'] - df['MUMODEL']
    elif 'MU_SH0ES' in df.columns and 'zHD' in df.columns:
        from scipy.integrate import quad
        H0, Om = 73.04, 0.334
        def dL(z):
            if z <= 0: return 1e-10
            dc, _ = quad(lambda zp: 1.0/np.sqrt(Om*(1+zp)**3 + (1-Om)), 0, z)
            return (1+z) * dc * (299792.458 / H0)
        mu_model = np.array([5*np.log10(dL(z))+25 for z in df['zHD']])
        df['mu_resid'] = df['MU_SH0ES'] - mu_model

    # Color residual
    z_edges = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.5]
    if 'c' in df.columns:
        df['c_resid'] = np.nan
        for i in range(len(z_edges)-1):
            mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
            if mask.sum() > 5:
                df.loc[mask, 'c_resid'] = df.loc[mask, 'c'] - df.loc[mask, 'c'].median()

    # x1 residual
    if 'x1' in df.columns:
        df['x1_resid'] = np.nan
        for i in range(len(z_edges)-1):
            mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
            if mask.sum() > 5:
                df.loc[mask, 'x1_resid'] = df.loc[mask, 'x1'] - df.loc[mask, 'x1'].median()

    # Fit quality
    if 'FITCHI2' in df.columns and 'NDOF' in df.columns:
        df['chi2_dof'] = df['FITCHI2'] / df['NDOF'].replace(0, np.nan)
    elif 'FITPROB' in df.columns:
        df['chi2_dof'] = -np.log10(df['FITPROB'].clip(1e-10, 1))

    # Survey label
    survey_map = {1: 'SDSS', 4: 'SNLS', 5: 'CSP', 10: 'DES', 15: 'PS1',
                  50: 'CfA3', 51: 'CfA4', 56: 'CfA1', 57: 'CfA2',
                  61: 'CfA3K', 62: 'CfA3S', 63: 'CfA4p1', 64: 'CfA4p2',
                  100: 'HST', 150: 'FOUNDATION'}
    if 'IDSURVEY' in df.columns:
        df['survey'] = df['IDSURVEY'].map(survey_map).fillna('Other')

    # SNR proxy (inverse of distance modulus error)
    if 'MU_SH0ES_ERR_DIAG' in df.columns:
        df['snr_proxy'] = 1.0 / df['MU_SH0ES_ERR_DIAG'].clip(0.001)
    elif 'mBERR' in df.columns:
        df['snr_proxy'] = 1.0 / df['mBERR'].clip(0.001)

    print(f"[*] Loaded {len(df)} SNe")
    return df


# ============================================================
# REST-FRAME WAVELENGTH COMPUTATION
# ============================================================

# Known observer-frame filter sets per survey (central wavelengths in Angstroms)
# These are approximate effective wavelengths for SALT2 fitting
SURVEY_FILTERS = {
    'SDSS':       {'u': 3551, 'g': 4686, 'r': 6166, 'i': 7480, 'z': 8932},
    'PS1':        {'g': 4866, 'r': 6215, 'i': 7545, 'z': 8679, 'y': 9633},
    'DES':        {'g': 4810, 'r': 6440, 'i': 7820, 'z': 9260},
    'SNLS':       {'g': 4870, 'r': 6250, 'i': 7700, 'z': 8900},
    'FOUNDATION': {'g': 4866, 'r': 6215, 'i': 7545, 'z': 8679},
    'CSP':        {'B': 4353, 'V': 5477, 'g': 4765, 'r': 6223, 'i': 7623},
    'CfA1':       {'B': 4353, 'V': 5477, 'R': 6349, 'I': 8797},
    'CfA2':       {'B': 4353, 'V': 5477, 'R': 6349, 'I': 8797},
    'CfA3':       {'B': 4353, 'V': 5477, 'r': 6223, 'i': 7623},
    'CfA3K':      {'B': 4353, 'V': 5477, 'R': 6349, 'I': 8797},
    'CfA3S':      {'B': 4353, 'V': 5477, 'r': 6223, 'i': 7623},
    'CfA4':       {'B': 4353, 'V': 5477, 'r': 6223, 'i': 7623},
    'CfA4p1':     {'B': 4353, 'V': 5477, 'r': 6223, 'i': 7623},
    'CfA4p2':     {'B': 4353, 'V': 5477, 'r': 6223, 'i': 7623},
    'HST':        {'F350LP': 5846, 'F606W': 5887, 'F775W': 7647, 'F850LP': 9033, 'F105W': 10552, 'F125W': 12486, 'F160W': 15369},
    'Other':      {'B': 4353, 'V': 5477, 'R': 6349, 'I': 8797},
}

# SALT2 rest-frame training range (where the model is well-constrained)
SALT2_BLUE_EDGE = 3000   # Angstroms (rest-frame)
SALT2_RED_EDGE = 7000    # Angstroms (rest-frame)
SALT2_UV_WEAK = 3500     # Below this, SALT2 training is sparse

def compute_restframe_coverage(df):
    """Compute rest-frame wavelength coverage metrics for each SN."""

    rf_blue = []    # Bluest rest-frame wavelength covered
    rf_red = []     # Reddest rest-frame wavelength covered
    rf_span = []    # Total rest-frame wavelength span
    uv_frac = []    # Fraction of filters probing rest-frame < 3500A (UV weak zone)
    n_rf_filters = []  # Number of filters in SALT2 training range

    for _, row in df.iterrows():
        z = row['zHD']
        survey = row.get('survey', 'Other')
        filters = SURVEY_FILTERS.get(survey, SURVEY_FILTERS['Other'])

        # Rest-frame wavelengths = observer-frame / (1+z)
        rf_wavelengths = np.array([lam / (1 + z) for lam in filters.values()])

        rf_blue.append(rf_wavelengths.min())
        rf_red.append(rf_wavelengths.max())
        rf_span.append(rf_wavelengths.max() - rf_wavelengths.min())

        # How many filters land in SALT2 training range?
        in_range = ((rf_wavelengths >= SALT2_BLUE_EDGE) &
                     (rf_wavelengths <= SALT2_RED_EDGE)).sum()
        n_rf_filters.append(in_range)

        # Fraction probing UV weak zone
        in_uv = (rf_wavelengths < SALT2_UV_WEAK).sum()
        uv_frac.append(in_uv / len(rf_wavelengths))

    df['rf_blue'] = rf_blue
    df['rf_red'] = rf_red
    df['rf_span'] = rf_span
    df['rf_uv_frac'] = uv_frac
    df['rf_n_good_filters'] = n_rf_filters

    # Coverage quality metric: how well does the filter set sample SALT2 range?
    # 1.0 = all filters in range with good UV, 0.0 = nothing useful
    df['coverage_quality'] = df['rf_n_good_filters'] / df.groupby('survey')['rf_n_good_filters'].transform('max').clip(1)

    print(f"[*] Rest-frame coverage computed")
    print(f"    rf_blue range: {df['rf_blue'].min():.0f} - {df['rf_blue'].max():.0f} A")
    print(f"    rf_red range:  {df['rf_red'].min():.0f} - {df['rf_red'].max():.0f} A")
    print(f"    UV fraction:   {df['rf_uv_frac'].mean():.2f} mean")

    return df


# ============================================================
# TEST D1: REST-FRAME WAVELENGTH STRESS TEST
# ============================================================

def test_D1_restframe_stress(df):
    """Does the color-distance entanglement depend on rest-frame wavelength coverage?

    Strategy:
    1. Split SNe by rest-frame coverage quality (good vs poor UV)
    2. Compute r(c, Δμ) vs z for each group
    3. Test if the threshold moves or disappears when coverage is controlled
    4. Partial correlation: Corr(c, Δμ | rf_blue, rf_span)
    """
    print(f"\n{'='*60}")
    print("TEST D1: Rest-Frame Wavelength Stress Test")
    print(f"{'='*60}")

    results = {}

    # ---- Panel 1: r(c, Δμ) vs z, split by UV coverage ----
    z_edges = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.5]
    z_mids = [(z_edges[i] + z_edges[i+1])/2 for i in range(len(z_edges)-1)]

    # Split by median rf_blue (proxy for UV stress)
    median_rf_blue = df['rf_blue'].median()
    good_uv = df['rf_blue'] <= median_rf_blue  # Lower = more UV coverage in rest frame
    poor_uv = df['rf_blue'] > median_rf_blue

    print(f"\n  Split by rest-frame blue edge (median={median_rf_blue:.0f} A):")
    print(f"    Good UV (rf_blue <= {median_rf_blue:.0f}): N={good_uv.sum()}")
    print(f"    Poor UV (rf_blue >  {median_rf_blue:.0f}): N={poor_uv.sum()}")

    fig, axes = plt.subplots(2, 3, figsize=(18, 11))

    for split_idx, (label, mask_split, color) in enumerate([
        ('Good UV coverage', good_uv, 'steelblue'),
        ('Poor UV coverage', poor_uv, 'crimson')
    ]):
        rs, ps, ns = [], [], []
        for i in range(len(z_edges)-1):
            zmask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
            mask = zmask & mask_split & df[['c', 'mu_resid']].notna().all(axis=1)
            n = mask.sum()
            ns.append(n)
            if n >= 15:
                r, p = stats.pearsonr(df.loc[mask, 'c'], df.loc[mask, 'mu_resid'])
                rs.append(r)
                ps.append(p)
            else:
                rs.append(np.nan)
                ps.append(np.nan)

        results[f'split_{label}'] = [
            {'z_bin': f'{z_edges[i]:.1f}-{z_edges[i+1]:.1f}', 'r': round(r, 4) if np.isfinite(r) else None,
             'p': round(p, 6) if np.isfinite(p) else None, 'N': int(n)}
            for i, (r, p, n) in enumerate(zip(rs, ps, ns))
        ]

        # Plot r(z) curve
        ax = axes[0, 0]
        valid = [np.isfinite(r) for r in rs]
        zv = [z for z, v in zip(z_mids, valid) if v]
        rv = [r for r, v in zip(rs, valid) if v]
        ax.plot(zv, [abs(r) for r in rv], 'o-', color=color, ms=8, lw=2, label=label)

        for i, (r, p, n) in enumerate(zip(rs, ps, ns)):
            sig = "***" if p and p < 0.01 else "**" if p and p < 0.05 else ""
            r_str = f"{r:.3f}" if np.isfinite(r) else "—"
            print(f"    {label} z=[{z_edges[i]:.1f},{z_edges[i+1]:.1f}): r={r_str}, N={n} {sig}")

    ax = axes[0, 0]
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('|r(c, Δμ)|', fontsize=12)
    ax.set_title('Color–Distance Coupling:\nGood vs Poor UV Coverage', fontsize=12)
    ax.legend(fontsize=10)
    ax.set_ylim(-0.02, 1.0)

    # ---- Panel 2: Scatter plot of rf_blue vs z, colored by |r| contribution ----
    ax = axes[0, 1]
    sc = ax.scatter(df['zHD'], df['rf_blue'], c=df['c'].abs(), cmap='viridis',
                    alpha=0.3, s=10, rasterized=True)
    ax.axhline(SALT2_UV_WEAK, color='red', ls='--', lw=1.5, label=f'SALT2 UV weak ({SALT2_UV_WEAK} Å)')
    ax.axhline(SALT2_BLUE_EDGE, color='orange', ls='--', lw=1.5, label=f'SALT2 blue edge ({SALT2_BLUE_EDGE} Å)')
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('Rest-frame bluest filter (Å)', fontsize=12)
    ax.set_title('Rest-Frame Coverage vs Redshift', fontsize=12)
    ax.legend(fontsize=9)
    plt.colorbar(sc, ax=ax, label='|SALT2 color c|')

    # ---- Panel 3: Does rf_blue predict the coupling strength? ----
    # Compute per-SN "coupling contribution" = |c * mu_resid| as proxy
    df['coupling_proxy'] = (df['c'] * df['mu_resid']).abs()

    ax = axes[0, 2]
    # Bin by rf_blue
    rf_bins = np.percentile(df['rf_blue'].dropna(), np.linspace(0, 100, 8))
    rf_mids, coupling_means, coupling_errs = [], [], []
    for i in range(len(rf_bins)-1):
        mask = (df['rf_blue'] >= rf_bins[i]) & (df['rf_blue'] < rf_bins[i+1])
        mask &= df['coupling_proxy'].notna()
        if mask.sum() > 10:
            rf_mids.append((rf_bins[i] + rf_bins[i+1])/2)
            coupling_means.append(df.loc[mask, 'coupling_proxy'].mean())
            coupling_errs.append(df.loc[mask, 'coupling_proxy'].std() / np.sqrt(mask.sum()))
    ax.errorbar(rf_mids, coupling_means, yerr=coupling_errs, fmt='o-', color='purple', ms=8, lw=2, capsize=3)
    ax.set_xlabel('Rest-frame blue edge (Å)', fontsize=12)
    ax.set_ylabel('Mean |c × Δμ|', fontsize=12)
    ax.set_title('Coupling Strength vs Coverage', fontsize=12)

    # ---- Panel 4: PARTIAL CORRELATION — the key test ----
    # Does r(c, Δμ) survive after partialling out rf_blue + rf_span?
    ax = axes[1, 0]

    partial_rs, raw_rs = [], []
    z_partial_mids = []

    for i in range(len(z_edges)-1):
        zmask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
        cols_needed = ['c', 'mu_resid', 'rf_blue', 'rf_span']
        mask = zmask & df[cols_needed].notna().all(axis=1)
        n = mask.sum()

        if n < 25:
            continue

        z_partial_mids.append((z_edges[i] + z_edges[i+1])/2)

        # Raw correlation
        r_raw, _ = stats.pearsonr(df.loc[mask, 'c'], df.loc[mask, 'mu_resid'])
        raw_rs.append(abs(r_raw))

        # Partial correlation: regress c and mu_resid on covariates, correlate residuals
        covariates = df.loc[mask, ['rf_blue', 'rf_span']].values
        X = np.column_stack([np.ones(n), covariates])

        # Residualize c
        coef_c = np.linalg.lstsq(X, df.loc[mask, 'c'].values, rcond=None)[0]
        c_resid_partial = df.loc[mask, 'c'].values - X @ coef_c

        # Residualize mu_resid
        coef_mu = np.linalg.lstsq(X, df.loc[mask, 'mu_resid'].values, rcond=None)[0]
        mu_resid_partial = df.loc[mask, 'mu_resid'].values - X @ coef_mu

        r_partial, p_partial = stats.pearsonr(c_resid_partial, mu_resid_partial)
        partial_rs.append(abs(r_partial))

        pct_survived = abs(r_partial) / max(abs(r_raw), 1e-10) * 100
        print(f"  z=[{z_edges[i]:.1f},{z_edges[i+1]:.1f}): raw |r|={abs(r_raw):.3f} → partial |r|={abs(r_partial):.3f} ({pct_survived:.0f}% survived)")

    ax.plot(z_partial_mids, raw_rs, 'o-', color='gray', ms=8, lw=2, label='Raw r(c, Δμ)')
    ax.plot(z_partial_mids, partial_rs, 's-', color='crimson', ms=8, lw=2, label='Partial r(c, Δμ | coverage)')
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('|Pearson r|', fontsize=12)
    ax.set_title('KEY TEST: Does Signal Survive\nAfter Controlling for Coverage?', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_ylim(-0.02, 1.0)

    results['partial_correlations'] = [
        {'z_mid': round(z, 3), 'raw_r': round(rr, 4), 'partial_r': round(pr, 4),
         'pct_survived': round(pr/max(rr, 1e-10)*100, 1)}
        for z, rr, pr in zip(z_partial_mids, raw_rs, partial_rs)
    ]

    # ---- Panel 5: FULL CONTROL — partial out rf_blue + rf_span + survey + host mass + SNR ----
    ax = axes[1, 1]

    full_partial_rs = []
    z_full_mids = []

    # Encode survey as numeric
    survey_codes = {s: i for i, s in enumerate(df['survey'].unique())}
    df['survey_code'] = df['survey'].map(survey_codes)

    for i in range(len(z_edges)-1):
        zmask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
        cols_needed = ['c', 'mu_resid', 'rf_blue', 'rf_span', 'HOST_LOGMASS', 'snr_proxy', 'survey_code']
        mask = zmask & df[cols_needed].notna().all(axis=1)
        n = mask.sum()

        if n < 30:
            continue

        z_full_mids.append((z_edges[i] + z_edges[i+1])/2)

        covariates = df.loc[mask, ['rf_blue', 'rf_span', 'HOST_LOGMASS', 'snr_proxy', 'survey_code']].values
        X = np.column_stack([np.ones(n), covariates])

        coef_c = np.linalg.lstsq(X, df.loc[mask, 'c'].values, rcond=None)[0]
        c_resid_full = df.loc[mask, 'c'].values - X @ coef_c

        coef_mu = np.linalg.lstsq(X, df.loc[mask, 'mu_resid'].values, rcond=None)[0]
        mu_resid_full = df.loc[mask, 'mu_resid'].values - X @ coef_mu

        r_full, p_full = stats.pearsonr(c_resid_full, mu_resid_full)
        full_partial_rs.append(abs(r_full))

        r_raw, _ = stats.pearsonr(df.loc[mask, 'c'], df.loc[mask, 'mu_resid'])
        pct = abs(r_full) / max(abs(r_raw), 1e-10) * 100
        print(f"  FULL CONTROL z=[{z_edges[i]:.1f},{z_edges[i+1]:.1f}): raw |r|={abs(r_raw):.3f} → full partial |r|={abs(r_full):.3f} ({pct:.0f}% survived), p={p_full:.4f}")

    # Plot all three curves
    ax.plot(z_partial_mids, raw_rs, 'o-', color='gray', ms=7, lw=1.5, label='Raw', alpha=0.7)
    ax.plot(z_partial_mids, partial_rs, 's-', color='steelblue', ms=7, lw=1.5, label='| coverage', alpha=0.8)
    ax.plot(z_full_mids, full_partial_rs, 'D-', color='crimson', ms=7, lw=2, label='| coverage+survey+mass+SNR')
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('|Pearson r|', fontsize=12)
    ax.set_title('Full Control: Signal After ALL Covariates', fontsize=12)
    ax.legend(fontsize=8)
    ax.set_ylim(-0.02, 1.0)

    results['full_control'] = [
        {'z_mid': round(z, 3), 'full_partial_r': round(r, 4)}
        for z, r in zip(z_full_mids, full_partial_rs)
    ]

    # ---- Panel 6: Sigmoid fit to controlled curve ----
    ax = axes[1, 2]

    def sigmoid(z, r_max, z0, k, r_min):
        return r_min + (r_max - r_min) / (1 + np.exp(-k * (z - z0)))

    # Fit sigmoid to raw curve
    try:
        z_arr = np.array(z_partial_mids)
        r_arr = np.array(raw_rs)
        popt_raw, _ = curve_fit(sigmoid, z_arr, r_arr, p0=[0.7, 0.5, 5, 0.1],
                                bounds=([0, 0, 0.1, 0], [1, 2, 50, 0.5]), maxfev=10000)
        results['sigmoid_raw'] = {'z_threshold': round(popt_raw[1], 3), 'steepness': round(popt_raw[2], 2)}
    except:
        popt_raw = None

    # Fit sigmoid to fully controlled curve
    try:
        z_fc = np.array(z_full_mids)
        r_fc = np.array(full_partial_rs)
        popt_ctrl, _ = curve_fit(sigmoid, z_fc, r_fc, p0=[0.7, 0.5, 5, 0.05],
                                 bounds=([0, 0, 0.1, 0], [1, 2, 50, 0.5]), maxfev=10000)
        results['sigmoid_controlled'] = {'z_threshold': round(popt_ctrl[1], 3), 'steepness': round(popt_ctrl[2], 2)}
    except:
        popt_ctrl = None

    z_fit = np.linspace(0, 1.5, 100)
    if popt_raw is not None:
        ax.plot(z_fit, sigmoid(z_fit, *popt_raw), '-', color='gray', lw=2, alpha=0.5,
                label=f'Raw threshold: z={popt_raw[1]:.2f}')
    if popt_ctrl is not None:
        ax.plot(z_fit, sigmoid(z_fit, *popt_ctrl), '-', color='crimson', lw=2,
                label=f'Controlled threshold: z={popt_ctrl[1]:.2f}')

    ax.plot(z_partial_mids, raw_rs, 'o', color='gray', ms=7, alpha=0.5)
    ax.plot(z_full_mids, full_partial_rs, 'D', color='crimson', ms=7)

    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('|r(c, Δμ)|', fontsize=12)
    ax.set_title('DECISIVE: Does Threshold Move\nAfter Controls?', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_ylim(-0.02, 1.0)

    if popt_raw is not None and popt_ctrl is not None:
        shift = popt_ctrl[1] - popt_raw[1]
        print(f"\n  *** THRESHOLD COMPARISON ***")
        print(f"  Raw sigmoid threshold:       z = {popt_raw[1]:.3f}")
        print(f"  Controlled sigmoid threshold: z = {popt_ctrl[1]:.3f}")
        print(f"  Shift: {shift:+.3f}")
        if abs(shift) < 0.1:
            print(f"  → THRESHOLD STABLE — signal survives coverage control ✅")
        else:
            print(f"  → THRESHOLD MOVED — coverage explains part of the effect ⚠️")
        results['threshold_shift'] = round(shift, 4)

    fig.suptitle('Test D1: Rest-Frame Wavelength Stress Test\n"Is this just K-correction breakdown?"',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test_D1_restframe_stress.png', dpi=150, bbox_inches='tight')
    plt.close()

    return results


# ============================================================
# TEST A1: CONDITIONAL INDEPENDENCE
# ============================================================

def test_A1_conditional_independence(df):
    """Does Corr(Δμ, c) survive after conditioning on ALL known covariates?

    If the signal collapses → conventional systematics.
    If it persists → missing latent variable (channel state).
    """
    print(f"\n{'='*60}")
    print("TEST A1: Conditional Independence After Controls")
    print(f"{'='*60}")

    results = {}

    # Define control sets (progressive)
    control_sets = {
        'raw': [],
        'z_only': ['zHD'],
        'z_survey': ['zHD', 'survey_code'],
        'z_survey_mass': ['zHD', 'survey_code', 'HOST_LOGMASS'],
        'z_survey_mass_snr': ['zHD', 'survey_code', 'HOST_LOGMASS', 'snr_proxy'],
        'full': ['zHD', 'survey_code', 'HOST_LOGMASS', 'snr_proxy', 'rf_blue', 'rf_span'],
    }

    z_edges = [0.0, 0.3, 0.5, 0.7, 2.5]  # Coarser bins for more power
    z_labels = ['low (0-0.3)', 'mid (0.3-0.5)', 'high (0.5-0.7)', 'very high (0.7+)']

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Panel 1: Bar chart — partial r at each control level, for high-z
    ax = axes[0]
    control_names = list(control_sets.keys())
    high_z_rs = []
    high_z_ps = []

    zmask_high = (df['zHD'] >= 0.5)

    for ctrl_name, ctrl_cols in control_sets.items():
        cols_needed = ['c', 'mu_resid'] + ctrl_cols
        mask = zmask_high & df[cols_needed].notna().all(axis=1)
        n = mask.sum()

        if n < 20 or not ctrl_cols:
            # Raw correlation
            r_val, p_val = stats.pearsonr(df.loc[mask, 'c'], df.loc[mask, 'mu_resid'])
        else:
            covariates = df.loc[mask, ctrl_cols].values
            X = np.column_stack([np.ones(n), covariates])
            coef_c = np.linalg.lstsq(X, df.loc[mask, 'c'].values, rcond=None)[0]
            coef_mu = np.linalg.lstsq(X, df.loc[mask, 'mu_resid'].values, rcond=None)[0]
            c_r = df.loc[mask, 'c'].values - X @ coef_c
            mu_r = df.loc[mask, 'mu_resid'].values - X @ coef_mu
            r_val, p_val = stats.pearsonr(c_r, mu_r)

        high_z_rs.append(abs(r_val))
        high_z_ps.append(p_val)
        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
        print(f"  z>0.5, controls={ctrl_name}: |r|={abs(r_val):.3f}, p={p_val:.4f} {sig}, N={n}")

    results['high_z_progressive'] = [
        {'controls': name, 'abs_r': round(r, 4), 'p': round(p, 6)}
        for name, r, p in zip(control_names, high_z_rs, high_z_ps)
    ]

    colors = ['gray', 'lightblue', 'steelblue', 'royalblue', 'darkblue', 'crimson']
    bars = ax.bar(range(len(control_names)), high_z_rs, color=colors, alpha=0.8)
    ax.set_xticks(range(len(control_names)))
    ax.set_xticklabels(control_names, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('|Partial r(c, Δμ)|', fontsize=12)
    ax.set_title('z > 0.5: Signal Strength After Progressive Controls', fontsize=12)
    ax.axhline(0.05, color='gray', ls=':', alpha=0.5, label='noise floor')
    ax.legend()

    # Annotate significance
    for i, (r, p) in enumerate(zip(high_z_rs, high_z_ps)):
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax.text(i, r + 0.01, sig, ha='center', fontsize=10, fontweight='bold')

    # Panel 2: Full control partial r across z bins
    ax = axes[1]

    for ctrl_idx, (ctrl_name, ctrl_cols) in enumerate(
        [('raw', []), ('full', ['zHD', 'survey_code', 'HOST_LOGMASS', 'snr_proxy', 'rf_blue', 'rf_span'])]
    ):
        z_fine = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.5]
        rs_ctrl = []
        zs_ctrl = []

        for i in range(len(z_fine)-1):
            zmask = (df['zHD'] >= z_fine[i]) & (df['zHD'] < z_fine[i+1])
            cols_needed = ['c', 'mu_resid'] + [c for c in ctrl_cols if c != 'zHD']  # Don't control z within z-bin
            mask = zmask & df[cols_needed].notna().all(axis=1)
            n = mask.sum()

            if n < 20:
                continue

            zs_ctrl.append((z_fine[i] + z_fine[i+1])/2)

            if not ctrl_cols or ctrl_cols == ['zHD']:
                r_val, _ = stats.pearsonr(df.loc[mask, 'c'], df.loc[mask, 'mu_resid'])
            else:
                use_cols = [c for c in ctrl_cols if c != 'zHD']
                covariates = df.loc[mask, use_cols].values
                X = np.column_stack([np.ones(n), covariates])
                coef_c = np.linalg.lstsq(X, df.loc[mask, 'c'].values, rcond=None)[0]
                coef_mu = np.linalg.lstsq(X, df.loc[mask, 'mu_resid'].values, rcond=None)[0]
                c_r = df.loc[mask, 'c'].values - X @ coef_c
                mu_r = df.loc[mask, 'mu_resid'].values - X @ coef_mu
                r_val, _ = stats.pearsonr(c_r, mu_r)

            rs_ctrl.append(abs(r_val))

        color = 'gray' if ctrl_name == 'raw' else 'crimson'
        style = 'o--' if ctrl_name == 'raw' else 'D-'
        ax.plot(zs_ctrl, rs_ctrl, style, color=color, ms=8, lw=2, label=ctrl_name)

    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('|r(c, Δμ)|', fontsize=12)
    ax.set_title('Raw vs Full-Control Across Redshift Bins', fontsize=12)
    ax.legend(fontsize=10)
    ax.set_ylim(-0.02, 1.0)

    fig.suptitle('Test A1: Conditional Independence — Is There a Missing Latent Variable?',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test_A1_conditional_independence.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Verdict
    survival_pct = high_z_rs[-1] / max(high_z_rs[0], 1e-10) * 100
    results['survival_pct'] = round(survival_pct, 1)
    print(f"\n  *** VERDICT ***")
    print(f"  Signal at z>0.5 after ALL controls: |r|={high_z_rs[-1]:.3f} ({survival_pct:.0f}% of raw)")
    if high_z_ps[-1] < 0.05:
        print(f"  → SIGNAL SURVIVES — model is missing a latent channel state ✅")
    else:
        print(f"  → SIGNAL COLLAPSES — conventional systematics explain it ❌")

    return results


# ============================================================
# TEST C1: SHARED CHANGE-POINT DETECTION
# ============================================================

def test_C1_shared_changepoint(df):
    """Do multiple independent metrics break at the same z*?

    Metrics:
    1. |r(c, Δμ)| per fine z-bin
    2. Effective rank per z-bin
    3. Mutual information proxy: MI(c; Δμ) per z-bin

    Method: For each metric, find the z that maximizes the difference between
    before/after means (CUSUM-like). If all three agree → phase transition.
    """
    print(f"\n{'='*60}")
    print("TEST C1: Shared Change-Point Detection")
    print(f"{'='*60}")

    results = {}

    # Fine z-bins
    z_edges = np.arange(0.0, 1.3, 0.05)
    z_mids = (z_edges[:-1] + z_edges[1:]) / 2

    # Metric 1: |r(c, Δμ)|
    abs_rs = []
    for i in range(len(z_edges)-1):
        mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
        mask &= df[['c', 'mu_resid']].notna().all(axis=1)
        if mask.sum() >= 15:
            r, _ = stats.pearsonr(df.loc[mask, 'c'], df.loc[mask, 'mu_resid'])
            abs_rs.append(abs(r))
        else:
            abs_rs.append(np.nan)

    # Metric 2: Effective rank
    obs_cols = [c for c in ['mu_resid', 'c', 'x1', 'chi2_dof'] if c in df.columns]
    eff_ranks = []
    for i in range(len(z_edges)-1):
        mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
        mask &= df[obs_cols].notna().all(axis=1)
        if mask.sum() >= 20:
            data = df.loc[mask, obs_cols].values
            data_std = (data - np.mean(data, axis=0)) / (np.std(data, axis=0) + 1e-10)
            eigenvalues = np.maximum(np.linalg.eigvalsh(np.cov(data_std.T)), 0)[::-1]
            ev_norm = eigenvalues / (eigenvalues.sum() + 1e-10)
            ev_norm = ev_norm[ev_norm > 1e-10]
            eff_ranks.append(np.exp(-np.sum(ev_norm * np.log(ev_norm))))
        else:
            eff_ranks.append(np.nan)

    # Metric 3: MI proxy — discretized mutual information
    mi_vals = []
    for i in range(len(z_edges)-1):
        mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
        mask &= df[['c', 'mu_resid']].notna().all(axis=1)
        n = mask.sum()
        if n >= 20:
            c_vals = df.loc[mask, 'c'].values
            mu_vals = df.loc[mask, 'mu_resid'].values
            # Discretize into quintiles and compute MI
            n_bins_mi = min(5, int(np.sqrt(n/5)))
            if n_bins_mi < 2:
                mi_vals.append(np.nan)
                continue
            c_d = pd.qcut(c_vals, n_bins_mi, labels=False, duplicates='drop')
            mu_d = pd.qcut(mu_vals, n_bins_mi, labels=False, duplicates='drop')
            # Contingency table MI
            from sklearn.metrics import mutual_info_score
            mi_vals.append(mutual_info_score(c_d, mu_d))
        else:
            mi_vals.append(np.nan)

    # Find change-point for each metric using max CUSUM
    def find_changepoint(values, z_arr):
        valid = [(z, v) for z, v in zip(z_arr, values) if np.isfinite(v)]
        if len(valid) < 6:
            return np.nan, 0
        zs = np.array([v[0] for v in valid])
        vs = np.array([v[1] for v in valid])
        best_z, best_stat = np.nan, 0
        for j in range(2, len(vs)-2):
            mean_before = np.mean(vs[:j])
            mean_after = np.mean(vs[j:])
            stat = abs(mean_after - mean_before) * np.sqrt(j * (len(vs)-j) / len(vs))
            if stat > best_stat:
                best_stat = stat
                best_z = zs[j]
        return best_z, best_stat

    cp_r, stat_r = find_changepoint(abs_rs, z_mids)
    cp_rank, stat_rank = find_changepoint([-r if np.isfinite(r) else np.nan for r in eff_ranks], z_mids)  # Negate: rank DROPS
    cp_mi, stat_mi = find_changepoint(mi_vals, z_mids)

    results['changepoints'] = {
        'abs_r': {'z_star': round(cp_r, 3), 'strength': round(stat_r, 3)},
        'eff_rank': {'z_star': round(cp_rank, 3), 'strength': round(stat_rank, 3)},
        'mutual_info': {'z_star': round(cp_mi, 3), 'strength': round(stat_mi, 3)},
    }

    print(f"\n  Change-points detected:")
    print(f"    |r(c, Δμ)|:      z* = {cp_r:.3f}  (stat = {stat_r:.3f})")
    print(f"    Effective rank:  z* = {cp_rank:.3f}  (stat = {stat_rank:.3f})")
    print(f"    Mutual info:     z* = {cp_mi:.3f}  (stat = {stat_mi:.3f})")

    cps = [cp_r, cp_rank, cp_mi]
    cps_valid = [c for c in cps if np.isfinite(c)]
    if len(cps_valid) >= 2:
        cp_spread = max(cps_valid) - min(cps_valid)
        cp_mean = np.mean(cps_valid)
        results['spread'] = round(cp_spread, 3)
        results['mean_z_star'] = round(cp_mean, 3)
        print(f"\n  Spread: {cp_spread:.3f}")
        print(f"  Mean z*: {cp_mean:.3f}")
        if cp_spread < 0.15:
            print(f"  → CONVERGED — shared phase transition ✅")
        elif cp_spread < 0.3:
            print(f"  → PARTIALLY CONVERGED ⚠️")
        else:
            print(f"  → DIVERGENT — different causes ❌")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    metrics = [
        ('|r(c, Δμ)|', abs_rs, cp_r, 'crimson'),
        ('Effective Rank', eff_ranks, cp_rank, 'steelblue'),
        ('MI(c; Δμ)', mi_vals, cp_mi, 'purple'),
    ]

    for idx, (name, vals, cp, color) in enumerate(metrics):
        ax = axes[idx]
        valid_mask = [np.isfinite(v) for v in vals]
        zv = [z for z, m in zip(z_mids, valid_mask) if m]
        vv = [v for v, m in zip(vals, valid_mask) if m]
        ax.plot(zv, vv, 'o-', color=color, ms=6, lw=1.5)
        if np.isfinite(cp):
            ax.axvline(cp, color=color, ls='--', lw=2, alpha=0.7, label=f'z* = {cp:.2f}')
        ax.set_xlabel('Redshift z', fontsize=12)
        ax.set_ylabel(name, fontsize=12)
        ax.set_title(f'{name}\nChange-point: z* = {cp:.2f}', fontsize=12)
        ax.legend()

    fig.suptitle('Test C1: Do All Metrics Break at the Same z*?\n(Closure predicts: yes — shared phase transition)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test_C1_shared_changepoint.png', dpi=150, bbox_inches='tight')
    plt.close()

    return results


# ============================================================
# TEST B1: PC1 DIRECTION STABILITY
# ============================================================

def test_B1_pc1_stability(df):
    """Does PC1 at high-z point the same direction regardless of survey/hemisphere?

    Stable PC1 → universal boundary mode.
    Rotating PC1 → survey-coded artifact.
    """
    print(f"\n{'='*60}")
    print("TEST B1: PC1 Direction Stability Across Splits")
    print(f"{'='*60}")

    results = {}
    obs_cols = [c for c in ['mu_resid', 'c', 'x1', 'chi2_dof'] if c in df.columns and df[c].notna().sum() > 100]

    # Define splits
    high_z_mask = df['zHD'] >= 0.5
    low_z_mask = df['zHD'] < 0.3

    # Survey splits (only surveys with enough high-z SNe)
    survey_counts_highz = df.loc[high_z_mask, 'survey'].value_counts()
    valid_surveys = survey_counts_highz[survey_counts_highz >= 20].index.tolist()

    # Hemisphere split
    north = df['DEC'] >= 0 if 'DEC' in df.columns else pd.Series([True]*len(df))
    south = ~north

    splits = {}
    # Survey splits at high-z
    for survey in valid_surveys:
        mask = high_z_mask & (df['survey'] == survey) & df[obs_cols].notna().all(axis=1)
        if mask.sum() >= 15:
            splits[f'highz_{survey}'] = mask

    # Hemisphere splits at high-z
    mask_hn = high_z_mask & north & df[obs_cols].notna().all(axis=1)
    mask_hs = high_z_mask & south & df[obs_cols].notna().all(axis=1)
    if mask_hn.sum() >= 15:
        splits['highz_north'] = mask_hn
    if mask_hs.sum() >= 15:
        splits['highz_south'] = mask_hs

    # All high-z and all low-z for comparison
    mask_all_high = high_z_mask & df[obs_cols].notna().all(axis=1)
    mask_all_low = low_z_mask & df[obs_cols].notna().all(axis=1)
    splits['all_highz'] = mask_all_high
    splits['all_lowz'] = mask_all_low

    # Compute PC1 for each split
    pc1_vectors = {}
    for name, mask in splits.items():
        n = mask.sum()
        if n < 15:
            continue
        data = df.loc[mask, obs_cols].values
        data_std = (data - np.mean(data, axis=0)) / (np.std(data, axis=0) + 1e-10)
        cov = np.cov(data_std.T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        pc1 = eigvecs[:, -1]  # Largest eigenvalue
        # Consistent sign convention: make largest component positive
        if pc1[np.argmax(np.abs(pc1))] < 0:
            pc1 = -pc1
        pc1_vectors[name] = pc1
        var_explained = eigvals[-1] / eigvals.sum() * 100
        print(f"  {name:25s} N={n:4d}  PC1=[{', '.join(f'{v:.3f}' for v in pc1)}]  var={var_explained:.1f}%")

    results['pc1_vectors'] = {k: [round(v, 4) for v in vec] for k, vec in pc1_vectors.items()}

    # Compute pairwise angles between high-z splits
    highz_keys = [k for k in pc1_vectors if k.startswith('highz_')]
    angles = {}
    if len(highz_keys) >= 2:
        print(f"\n  Pairwise angles between high-z PC1 vectors:")
        for k1, k2 in combinations(highz_keys, 2):
            cos_angle = np.clip(np.dot(pc1_vectors[k1], pc1_vectors[k2]), -1, 1)
            angle_deg = np.degrees(np.arccos(abs(cos_angle)))  # 0-90 range
            angles[f'{k1}_vs_{k2}'] = round(angle_deg, 1)
            print(f"    {k1} vs {k2}: {angle_deg:.1f}°")

        mean_angle = np.mean(list(angles.values()))
        results['mean_angle_highz'] = round(mean_angle, 1)
        print(f"\n  Mean angle between high-z splits: {mean_angle:.1f}°")
        if mean_angle < 20:
            print(f"  → PC1 is STABLE across splits — universal mode ✅")
        elif mean_angle < 45:
            print(f"  → PC1 is PARTIALLY STABLE ⚠️")
        else:
            print(f"  → PC1 ROTATES — survey-coded ❌")

    # Compare low-z vs high-z PC1
    if 'all_lowz' in pc1_vectors and 'all_highz' in pc1_vectors:
        cos_lh = np.clip(np.dot(pc1_vectors['all_lowz'], pc1_vectors['all_highz']), -1, 1)
        angle_lh = np.degrees(np.arccos(abs(cos_lh)))
        results['lowz_vs_highz_angle'] = round(angle_lh, 1)
        print(f"\n  Low-z vs High-z PC1 angle: {angle_lh:.1f}°")

    results['pairwise_angles'] = angles

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: PC1 components as bar chart
    ax = axes[0]
    x = np.arange(len(obs_cols))
    width = 0.8 / max(len(highz_keys), 1)
    cmap = plt.cm.Set2(np.linspace(0, 1, len(highz_keys)))

    for i, key in enumerate(highz_keys):
        if key in pc1_vectors:
            ax.bar(x + i * width - 0.4 + width/2, pc1_vectors[key], width,
                   label=key, color=cmap[i], alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(obs_cols, fontsize=10)
    ax.set_ylabel('PC1 Loading', fontsize=12)
    ax.set_title('PC1 Direction by Split (high-z)', fontsize=12)
    ax.legend(fontsize=8)
    ax.axhline(0, color='gray', ls='-', lw=0.5)

    # Panel 2: Angle matrix heatmap
    ax = axes[1]
    if angles:
        all_keys = sorted(set(k for pair in angles for k in pair.split('_vs_')))
        n_keys = len(all_keys)
        angle_matrix = np.zeros((n_keys, n_keys))
        for pair, angle in angles.items():
            k1, k2 = pair.split('_vs_')
            i1, i2 = all_keys.index(k1), all_keys.index(k2)
            angle_matrix[i1, i2] = angle
            angle_matrix[i2, i1] = angle

        im = ax.imshow(angle_matrix, cmap='RdYlGn_r', vmin=0, vmax=90)
        ax.set_xticks(range(n_keys))
        ax.set_yticks(range(n_keys))
        ax.set_xticklabels([k.replace('highz_', '') for k in all_keys], rotation=45, ha='right', fontsize=9)
        ax.set_yticklabels([k.replace('highz_', '') for k in all_keys], fontsize=9)
        for i in range(n_keys):
            for j in range(n_keys):
                if i != j:
                    ax.text(j, i, f'{angle_matrix[i,j]:.0f}°', ha='center', va='center', fontsize=9)
        plt.colorbar(im, ax=ax, label='Angle (degrees)')
        ax.set_title('PC1 Pairwise Angles\n(Green = stable, Red = rotating)', fontsize=12)
    else:
        ax.text(0.5, 0.5, 'Insufficient splits\nfor angle comparison', transform=ax.transAxes,
                ha='center', va='center', fontsize=14)

    fig.suptitle('Test B1: Is PC1 a Universal Boundary Mode or Survey Artifact?',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test_B1_pc1_stability.png', dpi=150, bbox_inches='tight')
    plt.close()

    return results


# ============================================================
# TEST F1: PREDICTABILITY ASYMMETRY
# ============================================================

def test_F1_predictability_asymmetry(df):
    """Is there a directional asymmetry in predictability at high-z?

    If the channel is lossy: predicting downstream from upstream should be
    easier than the reverse. Use simple linear regression R² as measure.
    """
    print(f"\n{'='*60}")
    print("TEST F1: Predictability Asymmetry")
    print(f"{'='*60}")

    results = {}

    z_edges = [0.0, 0.15, 0.3, 0.5, 0.7, 1.0, 2.5]

    z_mids, r2_c_from_mu, r2_mu_from_c, asymmetry = [], [], [], []

    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        mask = (df['zHD'] >= z_lo) & (df['zHD'] < z_hi)
        predictors = ['x1', 'chi2_dof', 'survey_code']
        cols = ['c', 'mu_resid'] + predictors
        mask &= df[cols].notna().all(axis=1)
        n = mask.sum()

        if n < 25:
            continue

        z_mids.append((z_lo + z_hi)/2)

        # Predict c from (Δμ, x1, χ², survey)
        X_for_c = df.loc[mask, ['mu_resid'] + predictors].values
        X_for_c = np.column_stack([np.ones(n), X_for_c])
        y_c = df.loc[mask, 'c'].values
        coef = np.linalg.lstsq(X_for_c, y_c, rcond=None)[0]
        ss_res = np.sum((y_c - X_for_c @ coef)**2)
        ss_tot = np.sum((y_c - y_c.mean())**2)
        r2_c = 1 - ss_res / max(ss_tot, 1e-10)
        r2_c_from_mu.append(max(r2_c, 0))

        # Predict Δμ from (c, x1, χ², survey)
        X_for_mu = df.loc[mask, ['c'] + predictors].values
        X_for_mu = np.column_stack([np.ones(n), X_for_mu])
        y_mu = df.loc[mask, 'mu_resid'].values
        coef = np.linalg.lstsq(X_for_mu, y_mu, rcond=None)[0]
        ss_res = np.sum((y_mu - X_for_mu @ coef)**2)
        ss_tot = np.sum((y_mu - y_mu.mean())**2)
        r2_mu = 1 - ss_res / max(ss_tot, 1e-10)
        r2_mu_from_c.append(max(r2_mu, 0))

        asym = r2_mu - r2_c  # Positive = easier to predict distance from color
        asymmetry.append(asym)

        print(f"  z=[{z_lo:.2f},{z_hi:.2f}): R²(c|Δμ,...)={r2_c:.3f}, R²(Δμ|c,...)={r2_mu:.3f}, asymmetry={asym:+.3f}")

    results['bins'] = [
        {'z_mid': round(z, 3), 'R2_c_from_mu': round(rc, 4), 'R2_mu_from_c': round(rm, 4), 'asymmetry': round(a, 4)}
        for z, rc, rm, a in zip(z_mids, r2_c_from_mu, r2_mu_from_c, asymmetry)
    ]

    # Trend in asymmetry
    if len(z_mids) >= 4:
        slope, _, r_trend, p_trend, _ = stats.linregress(z_mids, asymmetry)
        results['asymmetry_trend'] = {'slope': round(slope, 4), 'r': round(r_trend, 4), 'p': round(p_trend, 6)}
        print(f"\n  Asymmetry trend: slope={slope:.4f}, p={p_trend:.4f}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.plot(z_mids, r2_c_from_mu, 'o-', color='steelblue', ms=8, lw=2, label='R²(c | Δμ, x1, χ², survey)')
    ax.plot(z_mids, r2_mu_from_c, 's-', color='crimson', ms=8, lw=2, label='R²(Δμ | c, x1, χ², survey)')
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('R²', fontsize=12)
    ax.set_title('Predictability in Each Direction', fontsize=12)
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.plot(z_mids, asymmetry, 'D-', color='purple', ms=8, lw=2)
    ax.axhline(0, color='gray', ls='--', alpha=0.5)
    ax.fill_between(z_mids, 0, asymmetry, alpha=0.2, color='purple')
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('R²(Δμ|c) − R²(c|Δμ)', fontsize=12)
    ax.set_title('Predictability Asymmetry\n(+ = distance more predictable from color)', fontsize=12)

    fig.suptitle('Test F1: Is There Directional Information Loss?',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test_F1_predictability_asymmetry.png', dpi=150, bbox_inches='tight')
    plt.close()

    return results


# ============================================================
# MAIN
# ============================================================

def main():
    ensure_dir(OUTPUT_DIR)
    np.random.seed(SEED)

    print("=" * 60)
    print("CLOSURE THEORY — ROUND 3: REVERSE-DEDUCTION TESTS")
    print("Author: Humza Hafeez")
    print("=" * 60)

    df = load_data()
    df = compute_restframe_coverage(df)

    all_results = {}

    # D1: Rest-frame wavelength stress test
    all_results['D1_restframe_stress'] = test_D1_restframe_stress(df)

    # A1: Conditional independence
    all_results['A1_conditional_independence'] = test_A1_conditional_independence(df)

    # C1: Shared change-point
    try:
        from sklearn.metrics import mutual_info_score
        all_results['C1_shared_changepoint'] = test_C1_shared_changepoint(df)
    except ImportError:
        print("\n[!] sklearn not available — skipping C1 (add scikit-learn to requirements.txt)")
        all_results['C1_shared_changepoint'] = {'error': 'sklearn not available'}

    # B1: PC1 stability
    all_results['B1_pc1_stability'] = test_B1_pc1_stability(df)

    # F1: Predictability asymmetry
    all_results['F1_predictability_asymmetry'] = test_F1_predictability_asymmetry(df)

    # Save
    results_file = OUTPUT_DIR / 'round3_results.json'
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n[*] Results saved to {results_file}")

    print(f"\n{'='*60}")
    print("ROUND 3 SUMMARY")
    print(f"{'='*60}")
    print(f"D1: Rest-frame wavelength stress test — does signal survive coverage control?")
    print(f"A1: Conditional independence — is there a missing latent variable?")
    print(f"C1: Shared change-point — do all metrics break at the same z*?")
    print(f"B1: PC1 stability — universal mode or survey artifact?")
    print(f"F1: Predictability asymmetry — directional information loss?")
    print(f"\nOutput: {OUTPUT_DIR}/")


if __name__ == '__main__':
    main()
