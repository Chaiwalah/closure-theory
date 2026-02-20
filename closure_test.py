#!/usr/bin/env python3
"""
Closure Theory — Triple Coherence Test
=======================================
Tests whether SN Ia light curve observables correlate with line-of-sight
baryonic density (Planck tSZ y-map) in ways ΛCDM does not predict.

Tests:
  1. Triple Coherence: Δμ, c_resid, χ²/dof vs Planck y at each SN position
  2. Survey-Split Replication: same test per survey (SNLS, DES, PS1, CfA, CSP)
  3. Angular Scale Diagnosis: correlation vs aperture size (5', 10', 30', 60')
  4. High-z Color–Distance Coupling: r(c_resid, Δμ) in redshift bins
  5. Nested Model AIC/BIC: does adding a closure term improve the fit?

Data sources (auto-downloaded):
  - Pantheon+SH0ES: https://github.com/PantheonPlusSH0ES/DataRelease
  - Planck y-map (MILCA): ESA Planck Legacy Archive

Author: Humza Hafeez (Closure Theory)
Date: 2026-02-20
"""

import os
import sys
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
from scipy.optimize import minimize

warnings.filterwarnings('ignore')

# ============================================================
# CONFIG
# ============================================================
OUTPUT_DIR = Path("results")
DATA_DIR = Path("data")
APERTURES_ARCMIN = [5, 10, 30, 60]
Z_BINS = [0.0, 0.15, 0.3, 0.5, 0.7, 1.0, 2.5]
SEED = 42

# ============================================================
# DATA DOWNLOAD
# ============================================================

def ensure_dir(p):
    p.mkdir(parents=True, exist_ok=True)

def download_pantheon():
    """Download Pantheon+SH0ES data release."""
    url = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2BSH0ES_STAT%2BSYS.dat"
    lcparams_url = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2BSH0ES.dat"
    
    ensure_dir(DATA_DIR)
    
    dat_file = DATA_DIR / "Pantheon+SH0ES.dat"
    if not dat_file.exists():
        print("[*] Downloading Pantheon+SH0ES light curve parameters...")
        subprocess.run(["curl", "-fSL", "-o", str(dat_file), lcparams_url], check=True)
    
    return dat_file

def download_planck_ymap():
    """Download Planck MILCA y-map (full-sky HEALPix FITS).
    
    This is ~100MB. We use the MILCA reconstruction from Planck 2015.
    """
    ensure_dir(DATA_DIR)
    ymap_file = DATA_DIR / "ymap.fits"
    
    if not ymap_file.exists():
        # Planck MILCA y-map (2015 release, full resolution)
        url = "https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/COM_CompMap_YSZ_R2.00/milca_ymaps.fits"
        # Alternative: smaller HEALPix downgraded version
        alt_url = "https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_CompMap-YSZ_R2.00.fits"
        
        print("[*] Downloading Planck y-map (~100MB)...")
        print("    Trying IRSA mirror first...")
        result = subprocess.run(["curl", "-fSL", "-o", str(ymap_file), url], 
                              capture_output=True)
        if result.returncode != 0:
            print("    IRSA failed, trying ESA PLA...")
            result = subprocess.run(["curl", "-fSL", "-o", str(ymap_file), alt_url],
                                  capture_output=True)
            if result.returncode != 0:
                print("[!] Could not download y-map automatically.")
                print("    Please download manually from:")
                print("    https://pla.esac.esa.int/#maps -> Compton y-map (MILCA)")
                print(f"    Save as: {ymap_file}")
                return None
    
    return ymap_file

def load_pantheon(dat_file):
    """Load Pantheon+ data into a clean DataFrame."""
    print("[*] Loading Pantheon+ data...")
    df = pd.read_csv(dat_file, sep=r'\s+', comment='#')
    
    # Print available columns for debugging
    print(f"    Columns: {list(df.columns)}")
    print(f"    N = {len(df)} supernovae")
    
    # Expected columns (may vary by release):
    # CID, IDSURVEY, zHD, zHDERR, mB, mBERR, x1, x1ERR, c, cERR, 
    # HOST_LOGMASS, HOST_LOGMASS_ERR, PKMJD, PKMJDERR,
    # RA, DEC (or DECL), FITCHI2, NDOF, biasCor_mu, biasCorErr_mu, etc.
    
    # Standardize column names
    col_map = {}
    for col in df.columns:
        low = col.lower()
        if low in ('decl', 'dec'):
            col_map[col] = 'DEC'
        elif low == 'ra':
            col_map[col] = 'RA'
    df = df.rename(columns=col_map)
    
    return df

def extract_y_at_positions(ymap_file, ra_deg, dec_deg, apertures_arcmin):
    """Extract Planck y-map values at given sky positions for multiple apertures.
    
    Returns dict of {aperture_arcmin: array_of_y_values}.
    """
    import healpy as hp
    
    print(f"[*] Reading Planck y-map from {ymap_file}...")
    
    # Read the y-map (try different HDU indices)
    try:
        ymap = hp.read_map(ymap_file, field=0, verbose=False)
    except Exception:
        try:
            ymap = hp.read_map(ymap_file, field=0, hdu=1, verbose=False)
        except Exception as e:
            print(f"[!] Error reading y-map: {e}")
            print("    Trying to list HDUs...")
            from astropy.io import fits
            with fits.open(str(ymap_file)) as hdul:
                hdul.info()
            raise
    
    nside = hp.npix2nside(len(ymap))
    print(f"    NSIDE = {nside}, Npix = {len(ymap)}")
    
    # Convert RA/Dec to theta/phi (HEALPix convention)
    theta = np.radians(90.0 - dec_deg)  # colatitude
    phi = np.radians(ra_deg)
    
    results = {}
    
    for aperture in apertures_arcmin:
        radius_rad = np.radians(aperture / 60.0)
        y_values = np.zeros(len(ra_deg))
        
        for i in range(len(ra_deg)):
            # Get pixels within aperture
            vec = hp.ang2vec(theta[i], phi[i])
            ipix = hp.query_disc(nside, vec, radius_rad)
            
            # Mean y in aperture (exclude NaN/bad pixels)
            vals = ymap[ipix]
            good = np.isfinite(vals)
            if good.sum() > 0:
                y_values[i] = np.mean(vals[good])
            else:
                y_values[i] = np.nan
        
        n_good = np.isfinite(y_values).sum()
        print(f"    Aperture {aperture}': {n_good}/{len(y_values)} good extractions")
        results[aperture] = y_values
    
    return results

# ============================================================
# COMPUTE RESIDUALS
# ============================================================

def compute_residuals(df):
    """Compute Hubble residuals and derived quantities.
    
    Uses biasCor_mu if available, otherwise computes from Tripp formula.
    """
    print("[*] Computing residuals...")
    
    # Check what columns we have
    has_biascor = 'biasCor_mu' in df.columns or 'MU_SH0ES' in df.columns
    has_fitchi2 = 'FITCHI2' in df.columns or 'FITPROB' in df.columns
    
    # Distance modulus residual
    if 'MU_SH0ES' in df.columns and 'MUMODEL' in df.columns:
        df['mu_resid'] = df['MU_SH0ES'] - df['MUMODEL']
    elif 'biasCor_mu' in df.columns and 'mumodel' in df.columns:
        df['mu_resid'] = df['biasCor_mu'] - df['mumodel']
    elif 'mB' in df.columns and 'x1' in df.columns and 'c' in df.columns:
        # Tripp formula: μ = mB - M + α·x1 - β·c
        # Use standard values: α=0.15, β=3.1, M=-19.25 (approximate)
        alpha, beta, M0 = 0.15, 3.1, -19.25
        mu_obs = df['mB'] - M0 + alpha * df['x1'] - beta * df['c']
        
        # Simple ΛCDM prediction (flat, Ωm=0.3)
        from scipy.integrate import quad
        H0 = 73.0  # km/s/Mpc
        Om = 0.3
        
        def E(z):
            return np.sqrt(Om * (1+z)**3 + (1-Om))
        
        def dL(z):
            if z <= 0:
                return 1e-10
            dc, _ = quad(lambda zp: 1.0/E(zp), 0, z)
            return (1+z) * dc * (299792.458 / H0)  # Mpc
        
        mu_model = np.array([5 * np.log10(dL(z)) + 25 for z in df['zHD']])
        df['mu_resid'] = mu_obs - mu_model
    else:
        print("[!] Cannot compute distance residuals — missing columns")
        df['mu_resid'] = np.nan
    
    # Color residual (c deviation from population mean in z-bin)
    if 'c' in df.columns:
        df['c_resid'] = np.nan
        for i in range(len(Z_BINS) - 1):
            mask = (df['zHD'] >= Z_BINS[i]) & (df['zHD'] < Z_BINS[i+1])
            if mask.sum() > 5:
                df.loc[mask, 'c_resid'] = df.loc[mask, 'c'] - df.loc[mask, 'c'].median()
    
    # Fit quality
    if 'FITCHI2' in df.columns and 'NDOF' in df.columns:
        df['chi2_dof'] = df['FITCHI2'] / df['NDOF'].replace(0, np.nan)
    elif 'FITPROB' in df.columns:
        # Convert fit probability to approximate chi2/dof
        # FITPROB ~ 1 - CDF(chi2, dof) — not invertible without dof, use as proxy
        df['chi2_dof'] = -np.log10(df['FITPROB'].clip(1e-10, 1))  # higher = worse fit
    else:
        print("[!] No fit quality column found — chi2/dof test will be skipped")
        df['chi2_dof'] = np.nan
    
    # Survey label
    if 'IDSURVEY' in df.columns:
        survey_map = {1: 'SDSS', 4: 'SNLS', 5: 'CSP', 10: 'DES', 15: 'PS1',
                      50: 'CfA3', 51: 'CfA4', 56: 'CfA1', 57: 'CfA2',
                      61: 'CfA3K', 62: 'CfA3S', 63: 'CfA4p1', 64: 'CfA4p2',
                      100: 'HST', 150: 'FOUNDATION'}
        df['survey'] = df['IDSURVEY'].map(survey_map).fillna('Other')
    elif 'SET' in df.columns:
        df['survey'] = df['SET'].astype(str)
    else:
        df['survey'] = 'Unknown'
    
    n_valid = df[['mu_resid', 'c_resid', 'chi2_dof']].notna().sum()
    print(f"    Valid: mu_resid={n_valid.get('mu_resid',0)}, c_resid={n_valid.get('c_resid',0)}, chi2_dof={n_valid.get('chi2_dof',0)}")
    
    return df

# ============================================================
# TEST 1: TRIPLE COHERENCE
# ============================================================

def test_triple_coherence(df, y_col='y_10', output_prefix='test1'):
    """Test correlation of Δμ, c_resid, χ²/dof with Planck y."""
    print(f"\n{'='*60}")
    print(f"TEST 1: Triple Coherence (y aperture: {y_col})")
    print(f"{'='*60}")
    
    observables = {
        'mu_resid': ('Hubble Residual Δμ', 'Distance anomaly'),
        'c_resid': ('Color Residual c_resid', 'Spectral anomaly'),
        'chi2_dof': ('Fit Quality χ²/dof', 'Shape/information anomaly'),
    }
    
    results = {}
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for idx, (obs_col, (title, subtitle)) in enumerate(observables.items()):
        mask = df[[obs_col, y_col]].notna().all(axis=1)
        x = df.loc[mask, y_col].values
        y = df.loc[mask, obs_col].values
        
        if len(x) < 10:
            print(f"  [{obs_col}] Insufficient data (N={len(x)})")
            results[obs_col] = {'r': np.nan, 'p': np.nan, 'N': len(x)}
            continue
        
        # Clip outliers (3σ)
        for arr in [x, y]:
            med, mad = np.median(arr), 1.4826 * np.median(np.abs(arr - np.median(arr)))
            clip = np.abs(arr - med) < 5 * max(mad, 1e-10)
            # Don't actually remove, just note
        
        r, p = stats.pearsonr(x, y)
        rho, p_rho = stats.spearmanr(x, y)
        
        results[obs_col] = {
            'pearson_r': round(r, 4),
            'pearson_p': round(p, 6),
            'spearman_rho': round(rho, 4),
            'spearman_p': round(p_rho, 6),
            'N': int(mask.sum()),
        }
        
        print(f"  [{obs_col}] N={mask.sum()}")
        print(f"    Pearson  r = {r:.4f}  (p = {p:.2e})")
        print(f"    Spearman ρ = {rho:.4f} (p = {p_rho:.2e})")
        
        # Plot
        ax = axes[idx]
        ax.scatter(x, y, alpha=0.15, s=8, color='steelblue', rasterized=True)
        
        # Binned means
        n_bins = min(10, len(x) // 20)
        if n_bins >= 3:
            bins = np.percentile(x, np.linspace(0, 100, n_bins + 1))
            bin_centers, bin_means, bin_errs = [], [], []
            for j in range(len(bins) - 1):
                in_bin = (x >= bins[j]) & (x < bins[j+1])
                if in_bin.sum() > 3:
                    bin_centers.append(np.median(x[in_bin]))
                    bin_means.append(np.mean(y[in_bin]))
                    bin_errs.append(np.std(y[in_bin]) / np.sqrt(in_bin.sum()))
            ax.errorbar(bin_centers, bin_means, yerr=bin_errs, 
                       fmt='o-', color='crimson', ms=6, lw=2, capsize=3, zorder=5)
        
        ax.set_xlabel(f'Planck y ({y_col})', fontsize=12)
        ax.set_ylabel(title, fontsize=12)
        ax.set_title(f'{title}\nr={r:.3f} (p={p:.2e}), ρ={rho:.3f}', fontsize=11)
        ax.axhline(0, color='gray', ls='--', alpha=0.5)
    
    fig.suptitle('Test 1: Triple Coherence — SN Observables vs Baryonic Channel Density', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / f'{output_prefix}_triple_coherence.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return results

# ============================================================
# TEST 2: SURVEY-SPLIT REPLICATION
# ============================================================

def test_survey_split(df, y_col='y_10'):
    """Run triple coherence per survey to check replication."""
    print(f"\n{'='*60}")
    print("TEST 2: Survey-Split Replication")
    print(f"{'='*60}")
    
    surveys = df['survey'].value_counts()
    print(f"  Surveys: {dict(surveys)}")
    
    results = {}
    
    # Only test surveys with enough SNe
    valid_surveys = surveys[surveys >= 20].index.tolist()
    
    fig, axes = plt.subplots(len(valid_surveys), 3, 
                             figsize=(18, 4 * len(valid_surveys)),
                             squeeze=False)
    
    for i, survey in enumerate(valid_surveys):
        mask_survey = df['survey'] == survey
        df_sub = df[mask_survey].copy()
        
        print(f"\n  --- {survey} (N={mask_survey.sum()}) ---")
        
        survey_results = {}
        for j, obs_col in enumerate(['mu_resid', 'c_resid', 'chi2_dof']):
            mask = df_sub[[obs_col, y_col]].notna().all(axis=1)
            if mask.sum() < 10:
                survey_results[obs_col] = {'r': np.nan, 'p': np.nan, 'N': int(mask.sum())}
                axes[i, j].set_title(f'{survey}: {obs_col} (N={mask.sum()}, insufficient)')
                continue
            
            x = df_sub.loc[mask, y_col].values
            y = df_sub.loc[mask, obs_col].values
            r, p = stats.pearsonr(x, y)
            
            survey_results[obs_col] = {'r': round(r, 4), 'p': round(p, 6), 'N': int(mask.sum())}
            print(f"    {obs_col}: r={r:.4f}, p={p:.2e}, N={mask.sum()}")
            
            ax = axes[i, j]
            ax.scatter(x, y, alpha=0.2, s=10, color='steelblue')
            ax.set_title(f'{survey}: r={r:.3f} (p={p:.2e})', fontsize=10)
            ax.set_xlabel(f'y ({y_col})')
            ax.set_ylabel(obs_col)
            ax.axhline(0, color='gray', ls='--', alpha=0.5)
        
        results[survey] = survey_results
    
    fig.suptitle('Test 2: Survey-Split Replication', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test2_survey_split.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return results

# ============================================================
# TEST 3: ANGULAR SCALE DIAGNOSIS
# ============================================================

def test_angular_scale(df, apertures_arcmin):
    """Test which angular scale produces the strongest correlation."""
    print(f"\n{'='*60}")
    print("TEST 3: Angular Scale Diagnosis")
    print(f"{'='*60}")
    
    results = {}
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(apertures_arcmin)))
    
    for obs_idx, obs_col in enumerate(['mu_resid', 'c_resid', 'chi2_dof']):
        ax = axes[obs_idx]
        rs, ps = [], []
        
        for ap in apertures_arcmin:
            y_col = f'y_{ap}'
            if y_col not in df.columns:
                rs.append(np.nan)
                ps.append(np.nan)
                continue
            
            mask = df[[obs_col, y_col]].notna().all(axis=1)
            if mask.sum() < 10:
                rs.append(np.nan)
                ps.append(np.nan)
                continue
            
            r, p = stats.pearsonr(df.loc[mask, y_col], df.loc[mask, obs_col])
            rs.append(r)
            ps.append(p)
            print(f"  {obs_col} vs y_{ap}': r={r:.4f}, p={p:.2e}")
        
        results[obs_col] = {f'{ap}arcmin': {'r': round(r, 4), 'p': round(p, 6)} 
                           for ap, r, p in zip(apertures_arcmin, rs, ps)}
        
        ax.plot(apertures_arcmin, rs, 'o-', color='crimson', ms=8, lw=2)
        ax.set_xlabel('Aperture (arcmin)', fontsize=12)
        ax.set_ylabel('Pearson r', fontsize=12)
        ax.set_title(obs_col, fontsize=12)
        ax.axhline(0, color='gray', ls='--', alpha=0.5)
        ax.set_xscale('log')
        ax.set_xticks(apertures_arcmin)
        ax.set_xticklabels([str(a) for a in apertures_arcmin])
        
        # Add p-value annotations
        for ap, r, p in zip(apertures_arcmin, rs, ps):
            if np.isfinite(r):
                ax.annotate(f'p={p:.1e}', (ap, r), textcoords="offset points",
                          xytext=(0, 10), fontsize=8, ha='center')
    
    fig.suptitle('Test 3: Angular Scale Diagnosis — Where Does the Signal Live?', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test3_angular_scale.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return results

# ============================================================
# TEST 4: HIGH-Z COLOR–DISTANCE COUPLING
# ============================================================

def test_highz_coupling(df):
    """Test if color–distance residual correlation strengthens with z."""
    print(f"\n{'='*60}")
    print("TEST 4: High-z Color–Distance Coupling")
    print(f"{'='*60}")
    
    z_edges = [0.0, 0.15, 0.3, 0.5, 0.7, 1.0, 2.5]
    results = {}
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    bin_centers, bin_rs, bin_ps, bin_ns = [], [], [], []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        mask = (df['zHD'] >= z_lo) & (df['zHD'] < z_hi)
        mask &= df[['mu_resid', 'c']].notna().all(axis=1)
        
        n = mask.sum()
        z_mid = (z_lo + z_hi) / 2
        
        if n < 10:
            print(f"  z=[{z_lo:.2f}, {z_hi:.2f}): N={n} — skipped")
            bin_centers.append(z_mid)
            bin_rs.append(np.nan)
            bin_ps.append(np.nan)
            bin_ns.append(n)
            continue
        
        c_vals = df.loc[mask, 'c'].values
        mu_vals = df.loc[mask, 'mu_resid'].values
        
        r, p = stats.pearsonr(c_vals, mu_vals)
        rho, p_rho = stats.spearmanr(c_vals, mu_vals)
        
        print(f"  z=[{z_lo:.2f}, {z_hi:.2f}): N={n}, r={r:.4f} (p={p:.2e}), ρ={rho:.4f}")
        
        bin_centers.append(z_mid)
        bin_rs.append(r)
        bin_ps.append(p)
        bin_ns.append(n)
        
        results[f'z_{z_lo:.2f}_{z_hi:.2f}'] = {
            'pearson_r': round(r, 4), 'pearson_p': round(p, 6),
            'spearman_rho': round(rho, 4), 'N': int(n)
        }
        
        if i < len(axes):
            ax = axes[i]
            ax.scatter(c_vals, mu_vals, alpha=0.2, s=10, color='steelblue')
            ax.set_xlabel('Color c')
            ax.set_ylabel('Δμ')
            ax.set_title(f'z=[{z_lo:.2f}, {z_hi:.2f}) N={n}\nr={r:.3f} (p={p:.2e})')
            ax.axhline(0, color='gray', ls='--', alpha=0.5)
            ax.axvline(0, color='gray', ls='--', alpha=0.5)
    
    # Summary panel: r vs z
    if len(axes) > len(z_edges) - 1:
        ax = axes[-1]
        valid = [not np.isnan(r) for r in bin_rs]
        ax.plot([c for c, v in zip(bin_centers, valid) if v],
                [r for r, v in zip(bin_rs, valid) if v],
                'o-', color='crimson', ms=8, lw=2)
        ax.set_xlabel('Redshift z')
        ax.set_ylabel('Pearson r(c, Δμ)')
        ax.set_title('Color–Distance Coupling vs Redshift\n(Accumulation test)')
        ax.axhline(0, color='gray', ls='--', alpha=0.5)
    
    fig.suptitle('Test 4: Does Color–Distance Coupling Strengthen with z?', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test4_highz_coupling.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return results

# ============================================================
# TEST 5: NESTED MODEL AIC/BIC
# ============================================================

def test_nested_models(df, y_col='y_10'):
    """Compare AIC/BIC for models with and without closure term."""
    print(f"\n{'='*60}")
    print("TEST 5: Nested Model Comparison (AIC/BIC)")
    print(f"{'='*60}")
    
    mask = df[['mu_resid', y_col, 'zHD']].notna().all(axis=1)
    if 'HOST_LOGMASS' in df.columns:
        # Also require host mass for step correction
        pass  # Don't filter on it, just use where available
    
    df_clean = df[mask].copy()
    N = len(df_clean)
    
    if N < 30:
        print(f"  Insufficient data (N={N})")
        return {}
    
    y_resid = df_clean['mu_resid'].values
    y_bary = df_clean[y_col].values
    z_vals = df_clean['zHD'].values
    
    # Standardize predictors
    y_bary_std = (y_bary - np.mean(y_bary)) / (np.std(y_bary) + 1e-10)
    z_std = (z_vals - np.mean(z_vals)) / (np.std(z_vals) + 1e-10)
    
    def neg_loglik(params, X, y):
        pred = X @ params[:-1]
        sigma = np.exp(params[-1])
        resid = y - pred
        return 0.5 * N * np.log(2 * np.pi * sigma**2) + 0.5 * np.sum(resid**2) / sigma**2
    
    def fit_model(X, y, k):
        p0 = np.zeros(k + 1)  # k coefficients + 1 log-sigma
        p0[-1] = np.log(np.std(y))
        res = minimize(neg_loglik, p0, args=(X, y), method='Nelder-Mead',
                      options={'maxiter': 10000})
        nll = res.fun
        aic = 2 * (k + 1) + 2 * nll
        bic = (k + 1) * np.log(N) + 2 * nll
        return {'nll': round(nll, 2), 'aic': round(aic, 2), 'bic': round(bic, 2),
                'params': res.x[:-1].round(6).tolist(), 'sigma': round(np.exp(res.x[-1]), 4),
                'k': k}
    
    # Model 0: Baseline (intercept only)
    X0 = np.ones((N, 1))
    m0 = fit_model(X0, y_resid, 1)
    
    # Model 1: Baseline + Γ·y (closure term)
    X1 = np.column_stack([np.ones(N), y_bary_std])
    m1 = fit_model(X1, y_resid, 2)
    
    # Model 2: Baseline + Γ·y + Γ_z·(y × z) (closure + accumulation)
    X2 = np.column_stack([np.ones(N), y_bary_std, y_bary_std * z_std])
    m2 = fit_model(X2, y_resid, 3)
    
    # Model 3: Baseline + Γ_z·(y × z) only (pure accumulation)
    X3 = np.column_stack([np.ones(N), y_bary_std * z_std])
    m3 = fit_model(X3, y_resid, 2)
    
    results = {
        'M0_intercept_only': m0,
        'M1_closure_y': m1,
        'M2_closure_y_yz': m2,
        'M3_accumulation_yz': m3,
        'N': N,
    }
    
    # Compute deltas
    print(f"\n  N = {N}")
    print(f"  {'Model':<30} {'k':>3} {'AIC':>10} {'BIC':>10} {'ΔAIC':>8} {'ΔBIC':>8}")
    print(f"  {'-'*70}")
    
    models = [('M0: Intercept only', m0), 
              ('M1: + Γ·y', m1),
              ('M2: + Γ·y + Γ·y·z', m2),
              ('M3: + Γ·y·z only', m3)]
    
    for name, m in models:
        daic = m['aic'] - m0['aic']
        dbic = m['bic'] - m0['bic']
        print(f"  {name:<30} {m['k']:>3} {m['aic']:>10.1f} {m['bic']:>10.1f} {daic:>+8.1f} {dbic:>+8.1f}")
    
    print(f"\n  Interpretation:")
    best_bic = min(m['bic'] for _, m in models)
    for name, m in models:
        if m['bic'] == best_bic:
            print(f"  → Best model by BIC: {name}")
    
    dbic_closure = m1['bic'] - m0['bic']
    if dbic_closure < -6:
        print(f"  → Strong evidence for closure term (ΔBIC = {dbic_closure:.1f})")
    elif dbic_closure < -2:
        print(f"  → Positive evidence for closure term (ΔBIC = {dbic_closure:.1f})")
    elif dbic_closure < 2:
        print(f"  → No significant evidence either way (ΔBIC = {dbic_closure:.1f})")
    else:
        print(f"  → Closure term NOT favored (ΔBIC = {dbic_closure:.1f}, penalty > improvement)")
    
    return results

# ============================================================
# MAIN
# ============================================================

def main():
    ensure_dir(OUTPUT_DIR)
    np.random.seed(SEED)
    
    print("=" * 60)
    print("CLOSURE THEORY — TRIPLE COHERENCE TEST SUITE")
    print("Author: Humza Hafeez")
    print("=" * 60)
    
    # Step 1: Download data
    pantheon_file = download_pantheon()
    ymap_file = download_planck_ymap()
    
    # Step 2: Load Pantheon+
    df = load_pantheon(pantheon_file)
    
    # Step 3: Compute residuals
    df = compute_residuals(df)
    
    # Step 4: Extract y-map values
    if ymap_file and ymap_file.exists():
        if 'RA' not in df.columns or 'DEC' not in df.columns:
            print("[!] No RA/DEC columns — cannot cross-match with y-map")
            print("    Available columns:", list(df.columns))
            return
        
        y_dict = extract_y_at_positions(
            ymap_file,
            df['RA'].values, 
            df['DEC'].values,
            APERTURES_ARCMIN
        )
        
        for ap, y_vals in y_dict.items():
            df[f'y_{ap}'] = y_vals
    else:
        print("[!] No y-map available — running tests without baryonic cross-match")
        print("    Tests 1-3 will be skipped; Tests 4-5 will use z-only diagnostics")
    
    # Step 5: Run tests
    all_results = {}
    
    # Test 1: Triple coherence (use 10' as primary)
    primary_y = 'y_10' if 'y_10' in df.columns else None
    if primary_y:
        all_results['test1_triple_coherence'] = test_triple_coherence(df, primary_y)
    
    # Test 2: Survey-split replication
    if primary_y:
        all_results['test2_survey_split'] = test_survey_split(df, primary_y)
    
    # Test 3: Angular scale diagnosis
    if any(f'y_{ap}' in df.columns for ap in APERTURES_ARCMIN):
        all_results['test3_angular_scale'] = test_angular_scale(df, APERTURES_ARCMIN)
    
    # Test 4: High-z coupling (no y-map needed)
    all_results['test4_highz_coupling'] = test_highz_coupling(df)
    
    # Test 5: Nested models
    if primary_y:
        all_results['test5_nested_models'] = test_nested_models(df, primary_y)
    
    # Save results
    results_file = OUTPUT_DIR / 'closure_test_results.json'
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n[*] Results saved to {results_file}")
    
    # Save processed data
    csv_out = OUTPUT_DIR / 'pantheon_with_y.csv'
    df.to_csv(csv_out, index=False)
    print(f"[*] Processed data saved to {csv_out}")
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Output directory: {OUTPUT_DIR}/")
    print(f"Plots: test1_triple_coherence.png, test2_survey_split.png,")
    print(f"       test3_angular_scale.png, test4_highz_coupling.png")
    print(f"Data:  closure_test_results.json, pantheon_with_y.csv")
    print(f"\nSend these files back for analysis.")

if __name__ == '__main__':
    main()
