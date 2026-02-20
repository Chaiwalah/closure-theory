#!/usr/bin/env python3
"""
CLOSURE THEORY — QUASAR SPECTRAL LINE TEST
============================================

Test whether quasar emission line properties show the same
information-channel degradation pattern as SNe Ia:
  - Correlations between line properties strengthen at high-z
  - Sigmoid transition near z ≈ 0.8
  - Frequency-dependent observables affected preferentially
  - Effect is NOT explained by luminosity/selection

Data: Wu & Shen (2022) — 750,414 SDSS DR16Q quasars
      with PyQSOFit emission line measurements

Key emission lines (rest-frame):
  - Hα 6563Å, Hβ 4861Å (Balmer series — frequency ratio fixed by physics)
  - MgII 2798Å (UV)
  - CIV 1549Å (high-ionization UV)
  - CIII] 1909Å (semi-forbidden)
  - Lyα 1216Å (highest frequency)

Closure prediction:
  If photon information degrades with path length, emission line
  RATIOS and WIDTHS should become more correlated at high-z,
  especially for lines separated in frequency. Lines close in
  frequency should be less affected.

Control:
  Baldwin Effect (luminosity-EW anticorrelation) is well-known.
  Must control for luminosity to isolate redshift-dependent effects.
"""

import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
from astropy.io import fits
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# CONFIGURATION
# ============================================================

DATA_DIR = Path("data")
CATALOG_FILE = DATA_DIR / "dr16q_prop_catalog.fits"
# Download from: http://quasar.astro.illinois.edu/paper_data/DR16Q/
# File: dr16q_prop-merged.fits (~4GB)
# Or the smaller main catalog

# Line property columns in Wu & Shen 2022 catalog
# Format: {LINE}_{PROPERTY} where PROPERTY is:
#   FWHM, EW, PEAK, AREA, etc.
# We use EW (equivalent width) and FWHM (velocity width)

# Redshift bins — matched to SN Ia analysis
Z_BINS = [
    (0.1, 0.3),
    (0.3, 0.5),
    (0.5, 0.7),
    (0.7, 0.9),
    (0.9, 1.2),
    (1.2, 1.6),
    (1.6, 2.2),
    (2.2, 3.0),
    (3.0, 4.5),
]

# Line pairs for correlation analysis
# Each tuple: (line1, line2, rest_freq_separation_angstrom)
# Greater wavelength separation = greater frequency separation
LINE_PAIRS_EW = [
    # Close in frequency (should be LESS affected)
    ('HALPHA', 'HBETA', 1702),    # 6563-4861 = 1702Å
    ('CIII', 'CIV', 360),          # 1909-1549 = 360Å
    # Far in frequency (should be MORE affected)
    ('HALPHA', 'CIV', 5014),       # 6563-1549 = 5014Å
    ('HBETA', 'MGII', 2063),       # 4861-2798 = 2063Å
    ('MGII', 'CIV', 1249),         # 2798-1549 = 1249Å
    ('HALPHA', 'LYA', 5347),       # 6563-1216 = 5347Å
]


def load_catalog():
    """Load the Wu & Shen 2022 quasar properties catalog."""
    
    # Try the main catalog first
    possible_files = [
        DATA_DIR / "dr16q_prop_catalog.fits",
        DATA_DIR / "dr16q_prop-merged.fits",
        DATA_DIR / "DR16Q_v4.fits",
    ]
    
    catalog_path = None
    for f in possible_files:
        if f.exists():
            catalog_path = f
            break
    
    if catalog_path is None:
        print("ERROR: No catalog file found!")
        print("Download from: http://quasar.astro.illinois.edu/paper_data/DR16Q/")
        print("  or: https://data.sdss.org/sas/dr16/eboss/qso/DR16Q/DR16Q_v4.fits")
        print(f"Place in: {DATA_DIR}/")
        sys.exit(1)
    
    print(f"  Loading {catalog_path}...")
    hdul = fits.open(catalog_path)
    
    # Print available extensions and columns for debugging
    print(f"  Extensions: {[h.name for h in hdul]}")
    
    # Try to find the right extension
    data = None
    for ext in range(len(hdul)):
        if hasattr(hdul[ext], 'columns') and hdul[ext].columns is not None:
            cols = [c.name for c in hdul[ext].columns]
            # Look for redshift column
            if 'Z' in cols or 'Z_DR16Q' in cols or 'Z_SYS' in cols:
                data = hdul[ext].data
                print(f"  Using extension {ext} ({hdul[ext].name}), {len(data)} rows")
                print(f"  Sample columns: {cols[:30]}")
                break
    
    if data is None:
        print("  Could not find data extension with redshift column")
        print("  Available columns per extension:")
        for ext in range(len(hdul)):
            if hasattr(hdul[ext], 'columns') and hdul[ext].columns is not None:
                print(f"    Ext {ext}: {[c.name for c in hdul[ext].columns][:20]}")
        sys.exit(1)
    
    return data, hdul


def extract_line_properties(data):
    """Extract emission line EWs and FWHMs into a clean DataFrame."""
    
    cols = [c.name if hasattr(c, 'name') else c for c in data.columns] if hasattr(data, 'columns') else list(data.dtype.names)
    
    print(f"\n  Available columns ({len(cols)} total)")
    
    # Find emission line columns
    line_cols = [c for c in cols if any(line in c.upper() for line in 
                 ['HALPHA', 'HBETA', 'MGII', 'CIV', 'CIII', 'LYA', 'MG_II', 'C_IV', 'C_III'])]
    print(f"  Emission line columns found: {len(line_cols)}")
    if line_cols:
        print(f"  Sample: {line_cols[:20]}")
    
    # Build DataFrame with what we can find
    df = pd.DataFrame()
    
    # Redshift
    for z_col in ['Z_SYS', 'Z_DR16Q', 'Z', 'Z_FIT']:
        if z_col in cols:
            df['z'] = np.array(data[z_col])
            print(f"  Redshift column: {z_col}")
            break
    
    # S/N
    for sn_col in ['SN_MEDIAN_ALL', 'SN_MEDIAN', 'SNMEDIAN']:
        if sn_col in cols:
            df['sn'] = np.array(data[sn_col])
            break
    
    # Luminosity (for Baldwin Effect control)
    for l_col in ['LOGL5100', 'LOGL3000', 'LOGL1350', 'LOGL_BOL', 'LOGLBOL']:
        if l_col in cols:
            df[l_col.lower()] = np.array(data[l_col])
    
    # Emission line properties — try multiple naming conventions
    line_map = {
        'HALPHA': ['HALPHA', 'HA', 'H_ALPHA'],
        'HBETA': ['HBETA', 'HB', 'H_BETA'],
        'MGII': ['MGII', 'MG_II', 'MG2'],
        'CIV': ['CIV', 'C_IV', 'C4'],
        'CIII': ['CIII', 'C_III', 'C3', 'CIII]', 'CIII_COMPLEX'],
        'LYA': ['LYA', 'LYMAN_ALPHA', 'LY_ALPHA'],
    }
    
    properties = ['EW', 'FWHM', 'PEAK', 'AREA', 'LOGL']
    
    for line_name, variants in line_map.items():
        for prop in properties:
            found = False
            for variant in variants:
                # Try patterns like: HALPHA_EW, EW_HALPHA, HALPHA_BR_EW
                patterns = [
                    f"{variant}_{prop}",
                    f"{prop}_{variant}",
                    f"{variant}_BR_{prop}",  # broad component
                    f"{prop}_{variant}_BR",
                ]
                for pat in patterns:
                    matching = [c for c in cols if c.upper() == pat.upper()]
                    if matching:
                        col_key = f"{line_name}_{prop}"
                        df[col_key] = np.array(data[matching[0]])
                        found = True
                        break
                if found:
                    break
    
    # Report what we found
    line_props_found = [c for c in df.columns if any(line in c for line in line_map.keys())]
    print(f"\n  Extracted {len(line_props_found)} line property columns:")
    for c in sorted(line_props_found):
        valid = df[c].notna() & (df[c] != 0) & (df[c] != -1)
        print(f"    {c}: {valid.sum()} valid values")
    
    return df


def test_q1_ew_correlation_vs_z(df):
    """
    Q1: Do emission line EW correlations strengthen with redshift?
    
    Closure prediction: EW correlations between lines far apart in
    frequency should increase at high-z (information entanglement).
    Lines close in frequency: weaker effect.
    """
    print("\n--- Test Q1: EW Correlation vs Redshift ---")
    print("If closure signal exists, line correlations should strengthen at high-z")
    print("especially for lines separated in frequency.\n")
    
    results = []
    
    for line1, line2, sep in LINE_PAIRS_EW:
        col1 = f"{line1}_EW"
        col2 = f"{line2}_EW"
        
        if col1 not in df.columns or col2 not in df.columns:
            print(f"  {line1}-{line2} (Δλ={sep}Å): SKIPPED (missing data)")
            continue
        
        print(f"  {line1} vs {line2} (Δλ={sep}Å):")
        
        for z_lo, z_hi in Z_BINS:
            mask = (df['z'] >= z_lo) & (df['z'] < z_hi)
            mask &= df[col1].notna() & (df[col1] > 0)
            mask &= df[col2].notna() & (df[col2] > 0)
            
            sub = df[mask]
            if len(sub) < 30:
                continue
            
            # Log-space correlation (EWs span orders of magnitude)
            r, p = stats.spearmanr(np.log10(sub[col1]), np.log10(sub[col2]))
            
            print(f"    z=[{z_lo},{z_hi}): r={r:.3f}, p={p:.2e}, N={len(sub)}")
            results.append({
                'line1': line1, 'line2': line2, 
                'freq_sep': sep,
                'z_lo': z_lo, 'z_hi': z_hi,
                'r': r, 'p': p, 'N': len(sub)
            })
    
    return results


def test_q2_fwhm_correlation_vs_z(df):
    """
    Q2: Do emission line FWHM correlations strengthen with redshift?
    
    FWHM measures velocity width — primarily a kinematic property,
    less frequency-dependent. Closure predicts WEAKER effect on FWHM
    compared to EW (analogous to stretch vs color in SNe).
    """
    print("\n--- Test Q2: FWHM Correlation vs Redshift ---")
    print("FWHM is kinematic (like SN stretch). Closure predicts WEAKER effect.\n")
    
    results = []
    
    for line1, line2, sep in LINE_PAIRS_EW:
        col1 = f"{line1}_FWHM"
        col2 = f"{line2}_FWHM"
        
        if col1 not in df.columns or col2 not in df.columns:
            continue
        
        print(f"  {line1} vs {line2} (Δλ={sep}Å):")
        
        for z_lo, z_hi in Z_BINS:
            mask = (df['z'] >= z_lo) & (df['z'] < z_hi)
            mask &= df[col1].notna() & (df[col1] > 0)
            mask &= df[col2].notna() & (df[col2] > 0)
            
            sub = df[mask]
            if len(sub) < 30:
                continue
            
            r, p = stats.spearmanr(np.log10(sub[col1]), np.log10(sub[col2]))
            
            print(f"    z=[{z_lo},{z_hi}): r={r:.3f}, p={p:.2e}, N={len(sub)}")
            results.append({
                'line1': line1, 'line2': line2,
                'freq_sep': sep,
                'z_lo': z_lo, 'z_hi': z_hi,
                'r': r, 'p': p, 'N': len(sub)
            })
    
    return results


def test_q3_rank_compression(df):
    """
    Q3: Does the rank (effective dimensionality) of the emission line
    observable space compress at high-z?
    
    This is the quasar analog of the SN rank compression test.
    Build a matrix of all line properties per quasar, compute
    eigenvalue spectrum, see if effective rank drops at high-z.
    """
    print("\n--- Test Q3: Observable Rank Compression vs Redshift ---")
    print("Closure predicts fewer independent degrees of freedom at high-z.\n")
    
    # Collect all available line properties
    line_cols = [c for c in df.columns if any(line in c for line in 
                 ['HALPHA', 'HBETA', 'MGII', 'CIV', 'CIII', 'LYA'])
                 and any(prop in c for prop in ['EW', 'FWHM', 'LOGL'])]
    
    if len(line_cols) < 4:
        print(f"  Only {len(line_cols)} line columns available, need ≥4. Skipping.")
        return []
    
    print(f"  Using {len(line_cols)} line properties: {line_cols}")
    
    results = []
    
    for z_lo, z_hi in Z_BINS:
        mask = (df['z'] >= z_lo) & (df['z'] < z_hi)
        for col in line_cols:
            mask &= df[col].notna() & (df[col] > 0)
        
        sub = df.loc[mask, line_cols]
        if len(sub) < 50:
            continue
        
        # Log-transform, standardize
        log_data = np.log10(sub.values)
        log_data = (log_data - log_data.mean(axis=0)) / (log_data.std(axis=0) + 1e-10)
        
        # PCA — compute eigenvalue spectrum
        cov = np.cov(log_data.T)
        eigenvalues = np.linalg.eigvalsh(cov)[::-1]
        
        # Effective rank (participation ratio)
        eigenvalues = eigenvalues[eigenvalues > 0]
        p = eigenvalues / eigenvalues.sum()
        eff_rank = np.exp(-np.sum(p * np.log(p + 1e-15)))
        
        # Fraction of variance in first PC
        var_pc1 = eigenvalues[0] / eigenvalues.sum()
        
        z_mid = (z_lo + z_hi) / 2
        print(f"    z=[{z_lo},{z_hi}): eff_rank={eff_rank:.2f}/{len(line_cols)}, "
              f"PC1={var_pc1:.1%}, N={len(sub)}")
        
        results.append({
            'z_lo': z_lo, 'z_hi': z_hi, 'z_mid': z_mid,
            'eff_rank': eff_rank, 'max_rank': len(line_cols),
            'var_pc1': var_pc1, 'N': len(sub),
            'eigenvalues': eigenvalues.tolist()
        })
    
    # Test for trend
    if len(results) >= 4:
        z_mids = [r['z_mid'] for r in results]
        ranks = [r['eff_rank'] for r in results]
        r_trend, p_trend = stats.spearmanr(z_mids, ranks)
        print(f"\n  Rank vs z trend: r={r_trend:.3f}, p={p_trend:.4f}")
        print(f"  Closure predicts NEGATIVE trend (rank decreases with z)")
    
    return results


def test_q4_baldwin_controlled(df):
    """
    Q4: Does the EW correlation trend survive after controlling
    for the Baldwin Effect (luminosity-EW anticorrelation)?
    
    The Baldwin Effect is a known luminosity selection effect.
    Must remove it to isolate any redshift-dependent signal.
    """
    print("\n--- Test Q4: Baldwin-Controlled EW Correlations ---")
    print("Control for luminosity to remove Baldwin Effect.\n")
    
    # Find a luminosity column
    lum_col = None
    for l in ['logl5100', 'logl3000', 'logl1350', 'logl_bol', 'loglbol']:
        if l in df.columns:
            lum_col = l
            break
    
    if lum_col is None:
        print("  No luminosity column found. Skipping Baldwin control.")
        return []
    
    print(f"  Using luminosity: {lum_col}")
    results = []
    
    for line1, line2, sep in LINE_PAIRS_EW:
        col1 = f"{line1}_EW"
        col2 = f"{line2}_EW"
        
        if col1 not in df.columns or col2 not in df.columns:
            continue
        
        print(f"\n  {line1} vs {line2} (Δλ={sep}Å):")
        
        for z_lo, z_hi in Z_BINS:
            mask = (df['z'] >= z_lo) & (df['z'] < z_hi)
            mask &= df[col1].notna() & (df[col1] > 0)
            mask &= df[col2].notna() & (df[col2] > 0)
            mask &= df[lum_col].notna() & (df[lum_col] > 0)
            
            sub = df[mask]
            if len(sub) < 50:
                continue
            
            # Partial correlation: remove luminosity dependence
            log_ew1 = np.log10(sub[col1].values)
            log_ew2 = np.log10(sub[col2].values)
            lum = sub[lum_col].values
            
            # Residualize against luminosity
            from numpy.polynomial import polynomial as P
            coef1 = P.polyfit(lum, log_ew1, 1)
            coef2 = P.polyfit(lum, log_ew2, 1)
            resid1 = log_ew1 - P.polyval(lum, coef1)
            resid2 = log_ew2 - P.polyval(lum, coef2)
            
            r, p = stats.spearmanr(resid1, resid2)
            
            print(f"    z=[{z_lo},{z_hi}): r_partial={r:.3f}, p={p:.2e}, N={len(sub)}")
            results.append({
                'line1': line1, 'line2': line2,
                'freq_sep': sep,
                'z_lo': z_lo, 'z_hi': z_hi,
                'r_partial': r, 'p': p, 'N': len(sub),
                'controlled': lum_col
            })
    
    return results


def test_q5_frequency_fingerprint(df):
    """
    Q5: Does correlation strength scale with frequency separation?
    
    THE KEY TEST. Closure predicts that lines further apart in
    frequency should show STRONGER correlation increase at high-z.
    
    Plot: Δ(correlation) vs Δ(frequency) for low-z vs high-z.
    If positive slope: frequency-dependent channel degradation.
    """
    print("\n--- Test Q5: Frequency Fingerprint ---")
    print("Closure predicts: larger frequency separation → larger correlation increase at high-z.\n")
    
    z_lo_range = (0.3, 0.7)   # "low-z" reference
    z_hi_range = (0.7, 1.5)   # "high-z" test
    
    results = []
    
    for line1, line2, sep in LINE_PAIRS_EW:
        col1 = f"{line1}_EW"
        col2 = f"{line2}_EW"
        
        if col1 not in df.columns or col2 not in df.columns:
            continue
        
        # Low-z correlation
        mask_lo = (df['z'] >= z_lo_range[0]) & (df['z'] < z_lo_range[1])
        mask_lo &= df[col1].notna() & (df[col1] > 0)
        mask_lo &= df[col2].notna() & (df[col2] > 0)
        
        # High-z correlation
        mask_hi = (df['z'] >= z_hi_range[0]) & (df['z'] < z_hi_range[1])
        mask_hi &= df[col1].notna() & (df[col1] > 0)
        mask_hi &= df[col2].notna() & (df[col2] > 0)
        
        sub_lo = df[mask_lo]
        sub_hi = df[mask_hi]
        
        if len(sub_lo) < 30 or len(sub_hi) < 30:
            continue
        
        r_lo, _ = stats.spearmanr(np.log10(sub_lo[col1]), np.log10(sub_lo[col2]))
        r_hi, _ = stats.spearmanr(np.log10(sub_hi[col1]), np.log10(sub_hi[col2]))
        
        delta_r = abs(r_hi) - abs(r_lo)
        
        print(f"  {line1}-{line2} (Δλ={sep}Å): "
              f"r_low={r_lo:.3f} (N={len(sub_lo)}), "
              f"r_high={r_hi:.3f} (N={len(sub_hi)}), "
              f"Δ|r|={delta_r:+.3f}")
        
        results.append({
            'line1': line1, 'line2': line2,
            'freq_sep': sep,
            'r_low': r_lo, 'r_high': r_hi,
            'delta_r': delta_r,
            'N_low': len(sub_lo), 'N_high': len(sub_hi)
        })
    
    # Test: does delta_r correlate with frequency separation?
    if len(results) >= 3:
        seps = [r['freq_sep'] for r in results]
        deltas = [r['delta_r'] for r in results]
        r_freq, p_freq = stats.spearmanr(seps, deltas)
        print(f"\n  Frequency fingerprint: r={r_freq:.3f}, p={p_freq:.4f}")
        print(f"  Closure predicts POSITIVE (larger sep → larger Δr)")
    
    return results


def test_q6_sigmoid_threshold(df):
    """
    Q6: Is there a sigmoid-like transition in quasar line correlations?
    
    Fit sigmoid vs linear to correlation-vs-z for key line pairs.
    Closure predicts sigmoid transition near z ≈ 0.8.
    """
    print("\n--- Test Q6: Sigmoid Threshold Detection ---")
    print("Looking for sharp transition in line correlations.\n")
    
    from scipy.optimize import curve_fit
    
    def sigmoid(z, z0, k, A, B):
        return A / (1 + np.exp(-k * (z - z0))) + B
    
    results = []
    
    for line1, line2, sep in LINE_PAIRS_EW:
        col1 = f"{line1}_EW"
        col2 = f"{line2}_EW"
        
        if col1 not in df.columns or col2 not in df.columns:
            continue
        
        # Compute correlation in fine z bins
        z_mids = []
        corrs = []
        
        fine_bins = np.arange(0.2, 3.0, 0.15)
        for z_lo in fine_bins:
            z_hi = z_lo + 0.3  # overlapping bins for smoothness
            mask = (df['z'] >= z_lo) & (df['z'] < z_hi)
            mask &= df[col1].notna() & (df[col1] > 0)
            mask &= df[col2].notna() & (df[col2] > 0)
            
            sub = df[mask]
            if len(sub) < 50:
                continue
            
            r, _ = stats.spearmanr(np.log10(sub[col1]), np.log10(sub[col2]))
            z_mids.append((z_lo + z_hi) / 2)
            corrs.append(abs(r))
        
        if len(z_mids) < 6:
            continue
        
        z_arr = np.array(z_mids)
        r_arr = np.array(corrs)
        
        # Fit sigmoid
        try:
            popt, pcov = curve_fit(sigmoid, z_arr, r_arr,
                                   p0=[0.8, 5.0, 0.2, 0.3],
                                   bounds=([0.1, 0.1, -1, -1], [3.0, 50, 1, 1]),
                                   maxfev=5000)
            
            # Compare sigmoid vs linear
            ss_sig = np.sum((r_arr - sigmoid(z_arr, *popt)) ** 2)
            lin_coef = np.polyfit(z_arr, r_arr, 1)
            ss_lin = np.sum((r_arr - np.polyval(lin_coef, z_arr)) ** 2)
            
            # F-test (sigmoid has 4 params, linear has 2)
            n = len(z_arr)
            f_stat = ((ss_lin - ss_sig) / 2) / (ss_sig / (n - 4)) if ss_sig > 0 else 0
            
            print(f"  {line1}-{line2} (Δλ={sep}Å):")
            print(f"    Sigmoid z₀={popt[0]:.2f}, k={popt[1]:.1f}")
            print(f"    SS_sigmoid={ss_sig:.4f}, SS_linear={ss_lin:.4f}")
            print(f"    Sigmoid better: {ss_sig < ss_lin}, F={f_stat:.2f}")
            
            results.append({
                'line1': line1, 'line2': line2,
                'freq_sep': sep,
                'z0_sigmoid': popt[0], 'k_sigmoid': popt[1],
                'ss_sigmoid': ss_sig, 'ss_linear': ss_lin,
                'sigmoid_better': ss_sig < ss_lin,
                'f_stat': f_stat
            })
        except Exception as e:
            print(f"  {line1}-{line2}: sigmoid fit failed ({e})")
    
    return results


def main():
    print("=" * 60)
    print("CLOSURE THEORY — QUASAR SPECTRAL LINE TEST")
    print("=" * 60)
    
    # Step 1: Load data
    print("\n[1/7] Loading SDSS DR16Q catalog...")
    data, hdul = load_catalog()
    
    # Step 2: Extract line properties
    print("\n[2/7] Extracting emission line properties...")
    df = extract_line_properties(data)
    
    print(f"\n  Total quasars: {len(df)}")
    print(f"  Redshift range: {df['z'].min():.2f} — {df['z'].max():.2f}")
    
    # Quality cuts
    if 'sn' in df.columns:
        df = df[df['sn'] > 2]
        print(f"  After S/N > 2 cut: {len(df)}")
    
    all_results = {}
    
    # Step 3-7: Run tests
    print("\n[3/7] Test Q1: EW Correlation vs Redshift...")
    all_results['Q1'] = test_q1_ew_correlation_vs_z(df)
    
    print("\n[4/7] Test Q2: FWHM Correlation vs Redshift...")
    all_results['Q2'] = test_q2_fwhm_correlation_vs_z(df)
    
    print("\n[5/7] Test Q3: Observable Rank Compression...")
    all_results['Q3'] = test_q3_rank_compression(df)
    
    print("\n[6/7] Test Q4: Baldwin-Controlled Correlations...")
    all_results['Q4'] = test_q4_baldwin_controlled(df)
    
    print("\n[7/7] Tests Q5-Q6: Frequency Fingerprint & Sigmoid...")
    all_results['Q5'] = test_q5_frequency_fingerprint(df)
    all_results['Q6'] = test_q6_sigmoid_threshold(df)
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    print("""
CLOSURE PREDICTIONS FOR QUASARS:
  ✓ Q1: EW correlations strengthen at high-z
  ✓ Q2: FWHM correlations show WEAKER trend (kinematic, like stretch)
  ✓ Q3: Effective rank of line-property space decreases with z
  ✓ Q4: Signal survives after Baldwin Effect (luminosity) control
  ✓ Q5: Correlation increase scales with frequency separation
  ✓ Q6: Sigmoid transition near z ≈ 0.8
    
NULL (conventional) PREDICTION:
  × All trends explained by luminosity selection + Baldwin Effect
  × No frequency-dependent pattern
  × No sharp threshold
""")
    
    # Save results
    import json
    output = {}
    for key, val in all_results.items():
        if isinstance(val, list):
            output[key] = val
    
    outfile = Path("results_quasar_closure.json")
    with open(outfile, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"  Results saved to {outfile}")
    
    hdul.close()


if __name__ == '__main__':
    main()
