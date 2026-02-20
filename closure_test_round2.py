#!/usr/bin/env python3
"""
Closure Theory — Round 2: Boundary Tests
==========================================
Tests for closure boundary signatures: factorization collapse,
threshold detection, parameter saturation, and information compression.

Based on the framework:
- Closure manifests as failure of factorization, not amplitude loss
- Effects activate beyond a threshold, not continuously
- At sufficient scale, independent observables become entangled
- Adding parameters beyond a closure threshold degrades inference

Author: Humza Hafeez (Closure Theory)
Date: 2026-02-20
"""

import os
import json
import warnings
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
from itertools import combinations

warnings.filterwarnings('ignore')

OUTPUT_DIR = Path("results_round2")
DATA_DIR = Path("data")
SEED = 42

def ensure_dir(p):
    p.mkdir(parents=True, exist_ok=True)

def download_pantheon():
    """Download Pantheon+SH0ES data."""
    ensure_dir(DATA_DIR)
    dat_file = DATA_DIR / "Pantheon+SH0ES.dat"
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
    """Load and prepare Pantheon+ data."""
    dat_file = download_pantheon()
    df = pd.read_csv(dat_file, sep=r'\s+', comment='#')
    
    # Standardize column names
    col_map = {}
    for col in df.columns:
        if col.lower() in ('decl', 'dec'):
            col_map[col] = 'DEC'
        elif col.lower() == 'ra':
            col_map[col] = 'RA'
    df = df.rename(columns=col_map)
    
    # Compute distance residual
    if 'MU_SH0ES' in df.columns and 'MUMODEL' in df.columns:
        df['mu_resid'] = df['MU_SH0ES'] - df['MUMODEL']
    elif 'biasCor_mu' in df.columns and 'mumodel' in df.columns:
        df['mu_resid'] = df['biasCor_mu'] - df['mumodel']
    
    # Compute color residual (deviation from z-bin median)
    z_edges = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.5]
    if 'c' in df.columns:
        df['c_resid'] = np.nan
        for i in range(len(z_edges)-1):
            mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
            if mask.sum() > 5:
                df.loc[mask, 'c_resid'] = df.loc[mask, 'c'] - df.loc[mask, 'c'].median()
    
    # x1 residual (same approach)
    if 'x1' in df.columns:
        df['x1_resid'] = np.nan
        for i in range(len(z_edges)-1):
            mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
            if mask.sum() > 5:
                df.loc[mask, 'x1_resid'] = df.loc[mask, 'x1'] - df.loc[mask, 'x1'].median()
    
    # Fit quality proxy
    if 'FITCHI2' in df.columns and 'NDOF' in df.columns:
        df['chi2_dof'] = df['FITCHI2'] / df['NDOF'].replace(0, np.nan)
    elif 'FITPROB' in df.columns:
        df['chi2_dof'] = -np.log10(df['FITPROB'].clip(1e-10, 1))
    
    print(f"[*] Loaded {len(df)} SNe")
    print(f"    Columns available: {sorted(df.columns.tolist())}")
    
    return df


# ============================================================
# TEST 6: FULL FACTORIZATION COLLAPSE
# ============================================================

def test6_factorization_collapse(df):
    """Test ALL pairwise correlations between observables vs redshift.
    
    If closure is real, ALL cross-correlations should strengthen with z,
    not just c vs mu. Everything entangles near the boundary.
    """
    print(f"\n{'='*60}")
    print("TEST 6: Full Factorization Collapse")
    print(f"{'='*60}")
    
    observables = {}
    for col, label in [('mu_resid', 'Δμ'), ('c', 'color c'), ('x1', 'stretch x1'), 
                        ('c_resid', 'c_resid'), ('x1_resid', 'x1_resid'), ('chi2_dof', 'χ²/dof')]:
        if col in df.columns and df[col].notna().sum() > 100:
            observables[col] = label
    
    pairs = list(combinations(observables.keys(), 2))
    print(f"  Testing {len(pairs)} observable pairs across z-bins")
    
    z_edges = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.5]
    z_mids = [(z_edges[i] + z_edges[i+1])/2 for i in range(len(z_edges)-1)]
    
    results = {}
    
    # Compute correlations for all pairs in all z-bins
    pair_curves = {}
    for (col_a, col_b) in pairs:
        pair_key = f"{col_a}_vs_{col_b}"
        rs_by_z = []
        ps_by_z = []
        ns_by_z = []
        
        for i in range(len(z_edges)-1):
            mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
            mask &= df[[col_a, col_b]].notna().all(axis=1)
            n = mask.sum()
            
            if n < 10:
                rs_by_z.append(np.nan)
                ps_by_z.append(np.nan)
                ns_by_z.append(n)
                continue
            
            r, p = stats.pearsonr(df.loc[mask, col_a], df.loc[mask, col_b])
            rs_by_z.append(r)
            ps_by_z.append(p)
            ns_by_z.append(n)
        
        pair_curves[pair_key] = {
            'r_values': rs_by_z,
            'p_values': ps_by_z,
            'n_values': ns_by_z,
            'label_a': observables[col_a],
            'label_b': observables[col_b],
        }
        
        results[pair_key] = [
            {'z_bin': f'{z_edges[i]:.1f}-{z_edges[i+1]:.1f}', 
             'r': round(r, 4) if not np.isnan(r) else None, 
             'p': round(p, 6) if not np.isnan(p) else None, 
             'N': int(n)}
            for i, (r, p, n) in enumerate(zip(rs_by_z, ps_by_z, ns_by_z))
        ]
    
    # Plot: all pair correlations vs z
    fig, axes = plt.subplots(2, 1, figsize=(14, 12))
    
    # Panel 1: r(z) for all pairs
    ax = axes[0]
    cmap = plt.cm.tab20(np.linspace(0, 1, len(pair_curves)))
    
    entanglement_scores = {}  # Track which pairs show increasing |r| with z
    
    for idx, (pair_key, curve) in enumerate(pair_curves.items()):
        valid = [not np.isnan(r) for r in curve['r_values']]
        z_valid = [z for z, v in zip(z_mids, valid) if v]
        r_valid = [r for r, v in zip(curve['r_values'], valid) if v]
        
        label = f"{curve['label_a']} vs {curve['label_b']}"
        ax.plot(z_valid, [abs(r) for r in r_valid], 'o-', color=cmap[idx], 
                label=label, alpha=0.7, ms=5, lw=1.5)
        
        # Compute trend: does |r| increase with z?
        if len(z_valid) >= 4:
            slope, _, rr, pp, _ = stats.linregress(z_valid, [abs(r) for r in r_valid])
            entanglement_scores[pair_key] = {
                'slope': round(slope, 4),
                'trend_r': round(rr, 4),
                'trend_p': round(pp, 6),
                'label': label,
            }
    
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('|Pearson r|', fontsize=12)
    ax.set_title('Test 6: Do ALL Cross-Correlations Strengthen with z?\n(Closure predicts: yes — factorization collapses near boundary)', fontsize=13)
    ax.legend(fontsize=7, ncol=2, loc='upper left')
    ax.axhline(0, color='gray', ls='--', alpha=0.3)
    ax.set_ylim(-0.02, 1.0)
    
    # Panel 2: Entanglement score summary
    ax2 = axes[1]
    if entanglement_scores:
        sorted_pairs = sorted(entanglement_scores.items(), key=lambda x: x[1]['slope'], reverse=True)
        labels = [v['label'] for _, v in sorted_pairs]
        slopes = [v['slope'] for _, v in sorted_pairs]
        colors = ['crimson' if s > 0 and v['trend_p'] < 0.1 else 'steelblue' 
                  for (_, v), s in zip(sorted_pairs, slopes)]
        
        y_pos = range(len(labels))
        ax2.barh(y_pos, slopes, color=colors, alpha=0.7)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(labels, fontsize=8)
        ax2.set_xlabel('Slope of |r| vs z (positive = entanglement increases)', fontsize=11)
        ax2.set_title('Entanglement Trend per Pair\n(Red = significant increasing entanglement, p<0.1)', fontsize=12)
        ax2.axvline(0, color='black', lw=0.5)
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test6_factorization_collapse.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    results['entanglement_scores'] = entanglement_scores
    
    # Summary
    print("\n  Entanglement trends (slope of |r| vs z):")
    for pair_key, score in sorted(entanglement_scores.items(), key=lambda x: x[1]['slope'], reverse=True):
        sig = "***" if score['trend_p'] < 0.01 else "**" if score['trend_p'] < 0.05 else "*" if score['trend_p'] < 0.1 else ""
        print(f"    {score['label']:<35} slope={score['slope']:+.4f}  (p={score['trend_p']:.4f}) {sig}")
    
    n_increasing = sum(1 for s in entanglement_scores.values() if s['slope'] > 0 and s['trend_p'] < 0.1)
    n_total = len(entanglement_scores)
    print(f"\n  {n_increasing}/{n_total} pairs show significantly increasing entanglement with z")
    
    return results


# ============================================================
# TEST 7: THRESHOLD DETECTION
# ============================================================

def test7_threshold_detection(df):
    """Find the characteristic z where factorization starts failing.
    
    Fit sigmoid and broken power law to the r(z) curve for c vs mu.
    Look for a transition zone rather than gradual drift.
    """
    print(f"\n{'='*60}")
    print("TEST 7: Threshold Detection")
    print(f"{'='*60}")
    
    # Use finer z-bins for better resolution
    z_edges = np.arange(0.0, 1.6, 0.05)
    z_mids = (z_edges[:-1] + z_edges[1:]) / 2
    
    rs, ns = [], []
    for i in range(len(z_edges)-1):
        mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
        mask &= df[['c', 'mu_resid']].notna().all(axis=1)
        n = mask.sum()
        ns.append(n)
        if n >= 15:
            r, _ = stats.pearsonr(df.loc[mask, 'c'], df.loc[mask, 'mu_resid'])
            rs.append(abs(r))
        else:
            rs.append(np.nan)
    
    rs = np.array(rs)
    ns = np.array(ns)
    valid = np.isfinite(rs)
    z_v = z_mids[valid]
    r_v = rs[valid]
    n_v = ns[valid]
    
    results = {'fine_bins': [{'z_mid': round(z, 3), 'abs_r': round(r, 4), 'N': int(n)} 
                             for z, r, n in zip(z_mids, rs, ns) if np.isfinite(r)]}
    
    # Fit models
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Model A: Linear
    try:
        slope, intercept, rr, pp, se = stats.linregress(z_v, r_v)
        z_fit = np.linspace(0, 1.5, 100)
        r_fit_linear = slope * z_fit + intercept
        residuals_linear = r_v - (slope * z_v + intercept)
        ss_res_linear = np.sum(residuals_linear**2)
        
        results['linear'] = {'slope': round(slope, 4), 'intercept': round(intercept, 4),
                            'r': round(rr, 4), 'p': round(pp, 6), 'ss_res': round(ss_res_linear, 4)}
        
        ax = axes[0]
        ax.scatter(z_v, r_v, s=n_v/5, alpha=0.6, color='steelblue')
        ax.plot(z_fit, np.clip(r_fit_linear, 0, 1), 'r-', lw=2)
        ax.set_title(f'Linear: slope={slope:.3f}\nSS_res={ss_res_linear:.3f}', fontsize=11)
        ax.set_xlabel('z')
        ax.set_ylabel('|r(c, Δμ)|')
        ax.set_ylim(-0.05, 1.0)
    except Exception as e:
        print(f"  Linear fit failed: {e}")
    
    # Model B: Sigmoid (threshold model)
    def sigmoid(z, r_max, z0, k, r_min):
        return r_min + (r_max - r_min) / (1 + np.exp(-k * (z - z0)))
    
    try:
        popt, pcov = curve_fit(sigmoid, z_v, r_v, p0=[0.7, 0.4, 5, 0.1],
                               bounds=([0, 0, 0.1, 0], [1, 2, 50, 0.5]),
                               maxfev=10000)
        r_fit_sigmoid = sigmoid(z_fit, *popt)
        residuals_sigmoid = r_v - sigmoid(z_v, *popt)
        ss_res_sigmoid = np.sum(residuals_sigmoid**2)
        
        results['sigmoid'] = {
            'r_max': round(popt[0], 4), 'z_threshold': round(popt[1], 4),
            'steepness': round(popt[2], 4), 'r_min': round(popt[3], 4),
            'ss_res': round(ss_res_sigmoid, 4)
        }
        
        ax = axes[1]
        ax.scatter(z_v, r_v, s=n_v/5, alpha=0.6, color='steelblue')
        ax.plot(z_fit, r_fit_sigmoid, 'r-', lw=2)
        ax.axvline(popt[1], color='orange', ls='--', lw=2, label=f'z_threshold = {popt[1]:.2f}')
        ax.set_title(f'Sigmoid: z₀={popt[1]:.2f}, k={popt[2]:.1f}\nSS_res={ss_res_sigmoid:.3f}', fontsize=11)
        ax.set_xlabel('z')
        ax.set_ylabel('|r(c, Δμ)|')
        ax.set_ylim(-0.05, 1.0)
        ax.legend()
    except Exception as e:
        print(f"  Sigmoid fit failed: {e}")
        results['sigmoid'] = {'error': str(e)}
    
    # Model C: Broken power law (two regimes)
    def broken_powerlaw(z, r0, z_break, slope1, slope2):
        return np.where(z < z_break, 
                       r0 + slope1 * z, 
                       r0 + slope1 * z_break + slope2 * (z - z_break))
    
    try:
        popt2, pcov2 = curve_fit(broken_powerlaw, z_v, r_v, p0=[0.1, 0.3, 0.2, 0.8],
                                 bounds=([-0.5, 0.05, -2, -2], [0.5, 1.5, 2, 2]),
                                 maxfev=10000)
        r_fit_broken = broken_powerlaw(z_fit, *popt2)
        residuals_broken = r_v - broken_powerlaw(z_v, *popt2)
        ss_res_broken = np.sum(residuals_broken**2)
        
        results['broken_powerlaw'] = {
            'r0': round(popt2[0], 4), 'z_break': round(popt2[1], 4),
            'slope_before': round(popt2[2], 4), 'slope_after': round(popt2[3], 4),
            'ss_res': round(ss_res_broken, 4)
        }
        
        ax = axes[2]
        ax.scatter(z_v, r_v, s=n_v/5, alpha=0.6, color='steelblue')
        ax.plot(z_fit, np.clip(r_fit_broken, 0, 1), 'r-', lw=2)
        ax.axvline(popt2[1], color='orange', ls='--', lw=2, label=f'z_break = {popt2[1]:.2f}')
        ax.set_title(f'Broken: z_break={popt2[1]:.2f}\nslope1={popt2[2]:.2f}, slope2={popt2[3]:.2f}\nSS_res={ss_res_broken:.3f}', fontsize=11)
        ax.set_xlabel('z')
        ax.set_ylabel('|r(c, Δμ)|')
        ax.set_ylim(-0.05, 1.0)
        ax.legend()
    except Exception as e:
        print(f"  Broken power law fit failed: {e}")
    
    fig.suptitle('Test 7: Is There a Closure Threshold?\n(Sigmoid/broken = threshold activation; Linear = gradual drift)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test7_threshold_detection.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Compare models
    print("\n  Model comparison (lower SS_res = better fit):")
    for model in ['linear', 'sigmoid', 'broken_powerlaw']:
        if model in results and 'ss_res' in results[model]:
            print(f"    {model:<20} SS_res = {results[model]['ss_res']:.4f}")
    
    if 'sigmoid' in results and 'z_threshold' in results['sigmoid']:
        print(f"\n  Sigmoid threshold: z = {results['sigmoid']['z_threshold']:.2f}")
    if 'broken_powerlaw' in results and 'z_break' in results['broken_powerlaw']:
        print(f"  Broken PL break:   z = {results['broken_powerlaw']['z_break']:.2f}")
    
    return results


# ============================================================
# TEST 8: INFORMATION COMPRESSION (EFFECTIVE RANK)
# ============================================================

def test8_information_compression(df):
    """Measure effective dimensionality of observables vs z.
    
    If closure compresses information near the boundary, the effective
    rank of the observable covariance matrix should decrease with z.
    """
    print(f"\n{'='*60}")
    print("TEST 8: Information Compression (Effective Rank)")
    print(f"{'='*60}")
    
    obs_cols = [c for c in ['mu_resid', 'c', 'x1', 'chi2_dof'] if c in df.columns and df[c].notna().sum() > 100]
    print(f"  Using observables: {obs_cols}")
    
    z_edges = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.5]
    
    results = {}
    z_mids, eff_ranks, participation_ratios, n_sne = [], [], [], []
    
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        mask = (df['zHD'] >= z_lo) & (df['zHD'] < z_hi)
        mask &= df[obs_cols].notna().all(axis=1)
        n = mask.sum()
        
        z_mid = (z_lo + z_hi) / 2
        
        if n < 20:
            print(f"  z=[{z_lo:.1f}, {z_hi:.1f}): N={n} — skipped")
            continue
        
        # Compute covariance and correlation matrix
        data = df.loc[mask, obs_cols].values
        # Standardize
        data_std = (data - np.mean(data, axis=0)) / (np.std(data, axis=0) + 1e-10)
        
        cov_mat = np.cov(data_std.T)
        eigenvalues = np.linalg.eigvalsh(cov_mat)
        eigenvalues = np.maximum(eigenvalues, 0)  # numerical stability
        eigenvalues = eigenvalues[::-1]  # descending
        
        # Effective rank (Shannon entropy of normalized eigenvalues)
        ev_norm = eigenvalues / (eigenvalues.sum() + 1e-10)
        ev_norm = ev_norm[ev_norm > 1e-10]
        shannon_entropy = -np.sum(ev_norm * np.log(ev_norm))
        eff_rank = np.exp(shannon_entropy)
        
        # Participation ratio (another measure of effective dimensionality)
        pr = (eigenvalues.sum())**2 / (np.sum(eigenvalues**2) + 1e-10)
        
        # Fraction of variance in top eigenvalue
        top_frac = eigenvalues[0] / (eigenvalues.sum() + 1e-10)
        
        z_mids.append(z_mid)
        eff_ranks.append(eff_rank)
        participation_ratios.append(pr)
        n_sne.append(n)
        
        results[f'z_{z_lo:.1f}_{z_hi:.1f}'] = {
            'N': int(n),
            'effective_rank': round(eff_rank, 3),
            'participation_ratio': round(pr, 3),
            'top_eigenvalue_fraction': round(top_frac, 3),
            'eigenvalues': [round(e, 4) for e in eigenvalues],
        }
        
        print(f"  z=[{z_lo:.1f}, {z_hi:.1f}): N={n}, eff_rank={eff_rank:.2f}, PR={pr:.2f}, top_frac={top_frac:.2f}")
    
    # Fit trend
    if len(z_mids) >= 4:
        slope_er, _, r_er, p_er, _ = stats.linregress(z_mids, eff_ranks)
        slope_pr, _, r_pr, p_pr, _ = stats.linregress(z_mids, participation_ratios)
        results['trends'] = {
            'eff_rank_slope': round(slope_er, 4), 'eff_rank_r': round(r_er, 4), 'eff_rank_p': round(p_er, 6),
            'part_ratio_slope': round(slope_pr, 4), 'part_ratio_r': round(r_pr, 4), 'part_ratio_p': round(p_pr, 6),
        }
        print(f"\n  Effective rank trend: slope={slope_er:.4f}, r={r_er:.4f}, p={p_er:.4f}")
        print(f"  Participation ratio trend: slope={slope_pr:.4f}, r={r_pr:.4f}, p={p_pr:.4f}")
        
        if slope_er < 0 and p_er < 0.1:
            print(f"  → SIGNAL: Information compression detected (rank decreases with z)")
        elif slope_er > 0 and p_er < 0.1:
            print(f"  → ANTI-SIGNAL: Dimensionality increases with z")
        else:
            print(f"  → No significant trend")
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    max_rank = len(obs_cols)
    
    ax = axes[0]
    ax.plot(z_mids, eff_ranks, 'o-', color='crimson', ms=8, lw=2)
    ax.axhline(max_rank, color='gray', ls=':', alpha=0.5, label=f'Max rank = {max_rank}')
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('Effective Rank', fontsize=12)
    ax.set_title(f'Effective Rank of Observable Space\n(Lower = more compressed/entangled)', fontsize=12)
    ax.legend()
    if len(z_mids) >= 4:
        z_fit = np.linspace(min(z_mids), max(z_mids), 50)
        ax.plot(z_fit, slope_er * z_fit + (eff_ranks[0] - slope_er * z_mids[0]), 
                '--', color='crimson', alpha=0.5)
    
    ax = axes[1]
    ax.plot(z_mids, participation_ratios, 'o-', color='steelblue', ms=8, lw=2)
    ax.axhline(max_rank, color='gray', ls=':', alpha=0.5, label=f'Max = {max_rank}')
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('Participation Ratio', fontsize=12)
    ax.set_title(f'Participation Ratio\n(Lower = variance concentrated in fewer dimensions)', fontsize=12)
    ax.legend()
    
    fig.suptitle('Test 8: Does Information Compress Near the Boundary?', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test8_information_compression.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return results


# ============================================================
# TEST 9: RECONSTRUCTION DEGRADATION
# ============================================================

def test9_reconstruction_degradation(df):
    """Test if adding correction parameters helps LESS at high-z.
    
    At low-z, SALT2 corrections (alpha*x1 - beta*c + mass step) should 
    reduce scatter. Near the closure boundary, these corrections should
    become less effective because the observables aren't separable.
    """
    print(f"\n{'='*60}")
    print("TEST 9: Reconstruction Degradation")
    print(f"{'='*60}")
    
    if not all(c in df.columns for c in ['mu_resid', 'c', 'x1', 'zHD']):
        print("  Missing required columns")
        return {}
    
    z_edges = [0.0, 0.15, 0.3, 0.5, 0.7, 1.0, 2.5]
    
    results = {}
    z_mids = []
    scatter_raw = []        # Raw mu_resid scatter
    scatter_1param = []     # After correcting with c only
    scatter_2param = []     # After correcting with c + x1
    scatter_3param = []     # After correcting with c + x1 + host mass (if available)
    improvement_1 = []      # Fractional improvement from 1st param
    improvement_2 = []      # Fractional improvement from 2nd param
    
    has_mass = 'HOST_LOGMASS' in df.columns
    
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        cols_needed = ['mu_resid', 'c', 'x1']
        mask = (df['zHD'] >= z_lo) & (df['zHD'] < z_hi)
        mask &= df[cols_needed].notna().all(axis=1)
        n = mask.sum()
        
        if n < 20:
            continue
        
        z_mid = (z_lo + z_hi) / 2
        z_mids.append(z_mid)
        
        mu = df.loc[mask, 'mu_resid'].values
        c = df.loc[mask, 'c'].values
        x1 = df.loc[mask, 'x1'].values
        
        # Raw scatter
        s0 = np.std(mu)
        scatter_raw.append(s0)
        
        # 1-param correction: mu ~ a*c
        from numpy.linalg import lstsq
        X1 = np.column_stack([np.ones(n), c])
        coef1, res1, _, _ = lstsq(X1, mu, rcond=None)
        s1 = np.std(mu - X1 @ coef1)
        scatter_1param.append(s1)
        
        # 2-param correction: mu ~ a*c + b*x1
        X2 = np.column_stack([np.ones(n), c, x1])
        coef2, res2, _, _ = lstsq(X2, mu, rcond=None)
        s2 = np.std(mu - X2 @ coef2)
        scatter_2param.append(s2)
        
        # 3-param if host mass available
        if has_mass:
            mask3 = mask & df['HOST_LOGMASS'].notna()
            if mask3.sum() >= 20:
                mu3 = df.loc[mask3, 'mu_resid'].values
                X3 = np.column_stack([np.ones(mask3.sum()), 
                                      df.loc[mask3, 'c'].values,
                                      df.loc[mask3, 'x1'].values,
                                      df.loc[mask3, 'HOST_LOGMASS'].values])
                coef3, _, _, _ = lstsq(X3, mu3, rcond=None)
                s3 = np.std(mu3 - X3 @ coef3)
                scatter_3param.append(s3)
            else:
                scatter_3param.append(np.nan)
        
        # Improvement fractions
        imp1 = (s0 - s1) / s0 if s0 > 0 else 0
        imp2 = (s1 - s2) / s1 if s1 > 0 else 0
        improvement_1.append(imp1)
        improvement_2.append(imp2)
        
        results[f'z_{z_lo:.2f}_{z_hi:.2f}'] = {
            'N': int(n),
            'scatter_raw': round(s0, 4),
            'scatter_1param': round(s1, 4),
            'scatter_2param': round(s2, 4),
            'improvement_c': round(imp1, 4),
            'improvement_x1': round(imp2, 4),
        }
        
        print(f"  z=[{z_lo:.2f}, {z_hi:.2f}): N={n}, σ_raw={s0:.3f}, σ_1p={s1:.3f} ({imp1:.1%}↓), σ_2p={s2:.3f} ({imp2:.1%}↓)")
    
    # Trend in improvement
    if len(z_mids) >= 4:
        slope1, _, r1, p1, _ = stats.linregress(z_mids, improvement_1)
        slope2, _, r2, p2, _ = stats.linregress(z_mids, improvement_2)
        results['trends'] = {
            'c_correction_slope': round(slope1, 4), 'c_correction_p': round(p1, 6),
            'x1_correction_slope': round(slope2, 4), 'x1_correction_p': round(p2, 6),
        }
        print(f"\n  Color correction effectiveness trend: slope={slope1:.4f}, p={p1:.4f}")
        print(f"  Stretch correction effectiveness trend: slope={slope2:.4f}, p={p2:.4f}")
        
        if slope1 < 0 and p1 < 0.1:
            print(f"  → SIGNAL: Color correction becomes LESS effective at high-z")
        if slope2 < 0 and p2 < 0.1:
            print(f"  → SIGNAL: Stretch correction becomes LESS effective at high-z")
    
    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    ax = axes[0]
    ax.plot(z_mids, scatter_raw, 'o-', color='gray', label='Raw', ms=8, lw=2)
    ax.plot(z_mids, scatter_1param, 's-', color='steelblue', label='+ color', ms=8, lw=2)
    ax.plot(z_mids, scatter_2param, 'D-', color='crimson', label='+ color + stretch', ms=8, lw=2)
    if scatter_3param and not all(np.isnan(s) for s in scatter_3param):
        valid_3 = [(z, s) for z, s in zip(z_mids, scatter_3param) if not np.isnan(s)]
        if valid_3:
            ax.plot([z for z, s in valid_3], [s for z, s in valid_3], 
                    '^-', color='green', label='+ c + x1 + mass', ms=8, lw=2)
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('Scatter in Δμ', fontsize=12)
    ax.set_title('Hubble Residual Scatter vs Correction Level', fontsize=12)
    ax.legend()
    
    ax = axes[1]
    ax.plot(z_mids, improvement_1, 'o-', color='steelblue', ms=8, lw=2, label='Color correction')
    ax.plot(z_mids, improvement_2, 's-', color='crimson', ms=8, lw=2, label='Stretch correction (marginal)')
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('Fractional Scatter Reduction', fontsize=12)
    ax.set_title('Correction Effectiveness vs z\n(Closure predicts: decreases at high-z)', fontsize=12)
    ax.legend()
    ax.axhline(0, color='gray', ls='--', alpha=0.3)
    
    # Panel 3: Ratio of improvement
    ax = axes[2]
    if len(improvement_1) == len(improvement_2):
        ratio = [i2/max(i1, 0.001) for i1, i2 in zip(improvement_1, improvement_2)]
        ax.plot(z_mids, ratio, 'o-', color='purple', ms=8, lw=2)
        ax.set_xlabel('Redshift z', fontsize=12)
        ax.set_ylabel('Marginal x1 improvement / c improvement', fontsize=12)
        ax.set_title('Diminishing Returns Ratio\n(Closure predicts: additional params help less at high-z)', fontsize=12)
        ax.axhline(0, color='gray', ls='--', alpha=0.3)
    
    fig.suptitle('Test 9: Do Corrections Become Less Effective Near the Boundary?', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test9_reconstruction_degradation.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return results


# ============================================================
# MAIN
# ============================================================

def main():
    ensure_dir(OUTPUT_DIR)
    np.random.seed(SEED)
    
    print("=" * 60)
    print("CLOSURE THEORY — ROUND 2: BOUNDARY TESTS")
    print("Author: Humza Hafeez")
    print("=" * 60)
    
    df = load_data()
    
    all_results = {}
    all_results['test6_factorization_collapse'] = test6_factorization_collapse(df)
    all_results['test7_threshold_detection'] = test7_threshold_detection(df)
    all_results['test8_information_compression'] = test8_information_compression(df)
    all_results['test9_reconstruction_degradation'] = test9_reconstruction_degradation(df)
    
    # Save
    results_file = OUTPUT_DIR / 'round2_results.json'
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n[*] Results saved to {results_file}")
    
    # Summary
    print(f"\n{'='*60}")
    print("ROUND 2 SUMMARY")
    print(f"{'='*60}")
    print(f"Test 6: Factorization Collapse — do ALL pairs entangle with z?")
    print(f"Test 7: Threshold Detection — is there a critical z?")
    print(f"Test 8: Information Compression — does dimensionality shrink with z?")
    print(f"Test 9: Reconstruction Degradation — do corrections fail at high-z?")
    print(f"\nOutput: {OUTPUT_DIR}/")

if __name__ == '__main__':
    main()
