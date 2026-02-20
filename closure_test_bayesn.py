#!/usr/bin/env python3
"""
Closure Theory — BayeSN Non-SALT Replication Test
==================================================
THE KILL SHOT: If the closure signal (color-luminosity coupling strengthening
at high-z) persists when using BayeSN instead of SALT2/SALT3, then SALT model
mis-specification is ruled out as the source.

BayeSN (Mandel+2022, Grayling+2024) is a hierarchical Bayesian SED model that
fits dust (R_V) and intrinsic SN properties directly from light curves —
completely independent of SALT2/SALT3 color parameterization.

Strategy:
  1. Download Pantheon+ raw light curves (SNANA format) from GitHub
  2. Cross-match with our summary stats file for redshifts + metadata
  3. Fit a strategic subsample with BayeSN (GPU-accelerated)
  4. Extract BayeSN distance moduli + theta (intrinsic) + AV (dust)
  5. Compute residuals and test for closure signal
  6. Compare BayeSN vs SALT2 correlation structures

Author: Humza Hafeez (Closure Theory)
Date: 2026-02-20
"""

import os
import sys
import json
import warnings
import subprocess
import tempfile
import glob
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

warnings.filterwarnings('ignore')

OUTPUT_DIR = Path("results_bayesn")
DATA_DIR = Path("data")
LC_DIR = DATA_DIR / "photometry"
SEED = 42
np.random.seed(SEED)

# BayeSN needs these
try:
    from bayesn import SEDmodel
    HAS_BAYESN = True
except ImportError:
    HAS_BAYESN = False
    print("WARNING: BayeSN not installed. Will run SALT-independent analysis only.")

def ensure_dir(p):
    p.mkdir(parents=True, exist_ok=True)


# ============================================================
# STEP 1: Download Pantheon+ photometry
# ============================================================

def download_photometry():
    """Download raw SNANA light curves from Pantheon+ GitHub release."""
    ensure_dir(LC_DIR)

    # Check if we already have light curves
    existing = list(LC_DIR.glob("*.dat")) + list(LC_DIR.glob("*.DAT"))
    if len(existing) > 50:
        print(f"  Already have {len(existing)} light curve files")
        return existing

    # The photometry is organized by survey in subdirectories
    # We need to clone or download from the Pantheon+ data release
    print("  Downloading Pantheon+ photometry index...")

    # Get list of photometry files from GitHub API
    import urllib.request
    api_url = "https://api.github.com/repos/PantheonPlusSH0ES/DataRelease/git/trees/3facbb99276c7589349d8eceaac218ccd2ad0726?recursive=1"
    try:
        req = urllib.request.Request(api_url, headers={'User-Agent': 'closure-theory'})
        with urllib.request.urlopen(req, timeout=30) as resp:
            tree = json.loads(resp.read())
    except Exception as e:
        print(f"  Failed to get file tree: {e}")
        print("  Trying direct clone...")
        # Fallback: sparse clone just the photometry
        subprocess.run([
            "git", "clone", "--depth=1", "--filter=blob:none", "--sparse",
            "https://github.com/PantheonPlusSH0ES/DataRelease.git",
            str(DATA_DIR / "PantheonRelease")
        ], capture_output=True, timeout=120)
        subprocess.run([
            "git", "-C", str(DATA_DIR / "PantheonRelease"),
            "sparse-checkout", "set", "Pantheon+_Data/1_DATA/photometry"
        ], capture_output=True, timeout=60)
        # Copy files out
        src = DATA_DIR / "PantheonRelease" / "Pantheon+_Data" / "1_DATA" / "photometry"
        if src.exists():
            for f in src.rglob("*.dat"):
                (LC_DIR / f.name).write_bytes(f.read_bytes())
            for f in src.rglob("*.DAT"):
                (LC_DIR / f.name).write_bytes(f.read_bytes())
        existing = list(LC_DIR.glob("*.dat")) + list(LC_DIR.glob("*.DAT"))
        print(f"  Downloaded {len(existing)} light curve files")
        return existing

    # Download individual files (only .dat/.DAT files, skip directories)
    base_url = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon+_Data/1_DATA/photometry/"
    lc_files = [t for t in tree.get('tree', []) if t['path'].endswith(('.dat', '.DAT')) and t['type'] == 'blob']

    print(f"  Found {len(lc_files)} light curve files, downloading...")
    downloaded = 0
    for item in lc_files:
        fname = Path(item['path']).name
        outpath = LC_DIR / fname
        if outpath.exists():
            downloaded += 1
            continue
        url = base_url + item['path']
        try:
            subprocess.run(["curl", "-fSL", "-o", str(outpath), url],
                         capture_output=True, timeout=30)
            downloaded += 1
            if downloaded % 100 == 0:
                print(f"    {downloaded}/{len(lc_files)} downloaded...")
        except:
            pass

    existing = list(LC_DIR.glob("*.dat")) + list(LC_DIR.glob("*.DAT"))
    print(f"  Have {len(existing)} light curve files total")
    return existing


def download_summary():
    """Download Pantheon+SH0ES summary statistics."""
    dat_file = DATA_DIR / "pantheon_plus.dat"
    if not dat_file.exists():
        ensure_dir(DATA_DIR)
        urls = [
            "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat",
        ]
        for url in urls:
            result = subprocess.run(["curl", "-fSL", "-o", str(dat_file), url],
                                  capture_output=True, timeout=60)
            if result.returncode == 0:
                break
    return dat_file


# ============================================================
# STEP 2: Parse SNANA light curves
# ============================================================

def parse_snana_file(filepath):
    """Parse an SNANA-format light curve file. Returns dict with metadata + photometry."""
    meta = {}
    phot_rows = []
    in_phot = False

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('NOBS:'):
                meta['NOBS'] = int(line.split(':')[1].strip())
            elif line.startswith('NVAR:'):
                meta['NVAR'] = int(line.split(':')[1].strip())
            elif line.startswith('VARLIST:'):
                meta['VARLIST'] = line.split(':')[1].strip().split()
            elif line.startswith('OBS:'):
                in_phot = True
                parts = line.split()[1:]  # skip 'OBS:'
                phot_rows.append(parts)
            elif ':' in line and not in_phot:
                key, _, val = line.partition(':')
                key = key.strip()
                val = val.strip()
                # Try numeric
                try:
                    if '.' in val:
                        meta[key] = float(val)
                    else:
                        meta[key] = int(val)
                except ValueError:
                    meta[key] = val
            elif line.startswith('END:'):
                break

    if phot_rows and 'VARLIST' in meta:
        df = pd.DataFrame(phot_rows, columns=meta['VARLIST'])
        # Convert numeric columns
        for col in df.columns:
            try:
                df[col] = pd.to_numeric(df[col])
            except (ValueError, TypeError):
                pass
        meta['photometry'] = df

    return meta


# ============================================================
# STEP 3: Select strategic subsample
# ============================================================

def build_lc_index(lc_files):
    """Build an index mapping SNID -> filepath by reading headers."""
    index = {}
    for f in lc_files:
        try:
            with open(f, 'r') as fh:
                for line in fh:
                    if line.startswith('SNID:'):
                        snid = line.split(':')[1].strip().upper()
                        index[snid] = f
                        break
        except:
            pass
    # Also index by filename stem
    for f in lc_files:
        stem = Path(f).stem.upper()
        if stem not in index:
            index[stem] = f
    return index


def select_subsample(summary_df, lc_files, n_per_bin=15):
    """Select a stratified subsample across redshift bins.
    We want good coverage of the closure transition zone (z~0.5-1.0)."""

    z_bins = [0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.2, 2.5]

    print("  Building LC index from file headers...")
    lc_index = build_lc_index(lc_files)
    print(f"  Indexed {len(lc_index)} light curves")

    selected = []
    for i in range(len(z_bins)-1):
        mask = (summary_df['zHD'] >= z_bins[i]) & (summary_df['zHD'] < z_bins[i+1])
        bin_df = summary_df[mask].copy()

        # Match with available light curves using multiple ID strategies
        matched = []
        for _, row in bin_df.iterrows():
            cid = str(row.get('CID', '')).strip().upper()
            # Direct match is most reliable (works for DES, CfA, KAIT, etc.)
            if cid in lc_index:
                matched.append((row, lc_index[cid]))
                continue
            # Try without leading zeros, with SN prefix, DES format
            for name_try in [cid.lstrip('0'), f"SN{cid}",
                             f"DES_REAL_{cid.zfill(8)}"]:
                if name_try and name_try in lc_index:
                    matched.append((row, lc_index[name_try]))
                    break

        if len(matched) > n_per_bin:
            # Prioritize SNe with good fit quality
            matched.sort(key=lambda x: abs(x[0].get('FITPROB', 0.5) - 0.5))
            matched = matched[:n_per_bin]

        selected.extend(matched)
        print(f"  z=[{z_bins[i]:.1f}, {z_bins[i+1]:.1f}): {len(matched)} SNe selected")

    return selected


# ============================================================
# STEP 4: Fit with BayeSN
# ============================================================

def fit_bayesn_sample(selected_sne):
    """Fit selected SNe with BayeSN and extract parameters."""
    if not HAS_BAYESN:
        print("  BayeSN not available — skipping GPU fits")
        return None

    print(f"\n  Initializing BayeSN (T21 model)...")
    model = SEDmodel(load_model='T21_model')

    results = []
    n_total = len(selected_sne)

    for i, (summary_row, lc_path) in enumerate(selected_sne):
        snid = str(summary_row.get('CID', summary_row.get('SNID', ''))).strip()
        z = summary_row['zHD']
        print(f"  [{i+1}/{n_total}] Fitting {snid} (z={z:.3f})...")

        try:
            # Parse the light curve
            lc_data = parse_snana_file(lc_path)
            if 'photometry' not in lc_data:
                print(f"    No photometry found, skipping")
                continue

            phot = lc_data['photometry']

            # BayeSN filter mapping — survey-dependent
            survey_name = lc_data.get('SURVEY', '').strip().upper()
            survey_id = summary_row.get('IDSURVEY', 0)

            # Detect filter column
            if 'FLT' in phot.columns:
                filt_col = 'FLT'
            elif 'BAND' in phot.columns:
                filt_col = 'BAND'
            else:
                print(f"    No filter column found, skipping")
                continue

            unique_filts = [str(f).strip() for f in phot[filt_col].unique()]
            filters_header = lc_data.get('FILTERS', '')

            # Survey-specific filter maps
            # Pantheon+ uses SNANA single-letter codes defined by FILTERS: header
            # BayeSN built-in: g/r/i/z_PS1, g/r/i/z_DES, Bessell_B/V/R/I, etc.
            filt_map = {}

            if survey_name in ('DES', 'DES-SN'):
                # DES uses griz directly
                for f in unique_filts:
                    if f in ('g', 'r', 'i', 'z'):
                        filt_map[f] = f'{f}_DES'
            elif survey_name in ('PS1', 'PS1MD', 'PS1_LOW_Z', 'FOUNDATION'):
                for f in unique_filts:
                    if f in ('g', 'r', 'i', 'z'):
                        filt_map[f] = f'{f}_PS1'
            elif survey_name in ('SDSS', 'SDSS-II'):
                for f in unique_filts:
                    if f in ('g', 'r', 'i', 'z'):
                        filt_map[f] = f'{f}_SDSS'
            elif survey_name == 'SNLS':
                for f in unique_filts:
                    if f in ('g', 'r', 'i', 'z'):
                        filt_map[f] = f'{f}_MegaCam'
            elif survey_name in ('CSP', 'CSP1', 'CSP3'):
                # CSP uses B, V, g, r, i in Swope filters
                for f in unique_filts:
                    if f in ('B', 'V'):
                        filt_map[f] = f'Bessell_{f}'
                    elif f in ('g', 'r', 'i'):
                        filt_map[f] = f'{f}_PS1'  # approximate
            elif survey_name == 'HST':
                # HST uses specific filter codes; skip for now (complex)
                pass
            else:
                # KAIT, CfA, etc. use single-letter SNANA codes
                # FILTERS: ABCDEFGHIJKLMNOP — these are survey-specific
                # Map common ones: try Bessell as fallback
                # Check if filters_header can give us band info
                # For KAIT/CfA with ABCD codes, we can't reliably map
                # Try: if filters look like standard bands
                for f in unique_filts:
                    if f in ('B', 'V', 'R', 'I', 'U'):
                        filt_map[f] = f'Bessell_{f}'
                    elif f in ('g', 'r', 'i', 'z'):
                        filt_map[f] = f'{f}_PS1'

            if len(filt_map) < 2:
                print(f"    Only {len(filt_map)} recognized filters ({survey_name}, bands={unique_filts[:5]}), skipping")
                continue

            # Try fitting with fit_from_file first (handles SNANA format natively)
            try:
                samples, sn_props = model.fit_from_file(
                    str(lc_path),
                    filt_map=filt_map,
                    num_warmup=200,
                    num_samples=200,
                    num_chains=1
                )
            except Exception as e1:
                # Fallback: extract arrays and use model.fit()
                print(f"    fit_from_file failed ({e1}), trying manual fit...")
                try:
                    mjd = phot['MJD'].values.astype(float)
                    flux = phot['FLUXCAL'].values.astype(float) if 'FLUXCAL' in phot.columns else phot['FLUX'].values.astype(float)
                    flux_err = phot['FLUXCALERR'].values.astype(float) if 'FLUXCALERR' in phot.columns else phot['FLUXERR'].values.astype(float)
                    filters = phot[filt_col].values

                    # Get MW E(B-V)
                    ebv_mw = lc_data.get('MWEBV', summary_row.get('MWEBV', 0.02))
                    peak_mjd = lc_data.get('SEARCH_PEAKMJD', lc_data.get('PEAKMJD',
                                summary_row.get('PKMJD', np.median(mjd))))

                    samples, sn_props = model.fit(
                        mjd, flux, flux_err, filters, z,
                        peak_mjd=peak_mjd, ebv_mw=ebv_mw,
                        filt_map=filt_map, mag=False,
                        num_warmup=200, num_samples=200, num_chains=1
                    )
                except Exception as e2:
                    print(f"    Manual fit also failed: {e2}")
                    continue

            # Extract BayeSN parameters
            # samples contains: mu (distance modulus), theta (intrinsic), AV (dust), tmax, epsilon
            mu_bayesn = np.median(samples['mu']) if 'mu' in samples else np.nan
            mu_err = np.std(samples['mu']) if 'mu' in samples else np.nan
            theta = np.median(samples['theta']) if 'theta' in samples else np.nan
            AV = np.median(samples['AV']) if 'AV' in samples else np.nan

            results.append({
                'CID': snid,
                'zHD': z,
                'mu_bayesn': mu_bayesn,
                'mu_bayesn_err': mu_err,
                'theta_bayesn': theta,
                'AV_bayesn': AV,
                'c_salt': summary_row.get('c', np.nan),
                'x1_salt': summary_row.get('x1', np.nan),
                'mu_salt': summary_row.get('MU_SH0ES', np.nan),
                'IDSURVEY': survey,
                'n_obs': len(phot),
                'n_filters': len(filt_map),
            })
            print(f"    mu_BayeSN={mu_bayesn:.3f}±{mu_err:.3f}, theta={theta:.3f}, AV={AV:.3f}")

        except Exception as e:
            print(f"    Error: {e}")
            continue

    if results:
        return pd.DataFrame(results)
    return None


# ============================================================
# STEP 5: SALT-Independent Analysis (works without BayeSN)
# ============================================================

def salt_independent_analysis(summary_df):
    """
    Even without BayeSN, we can test whether the closure signal
    is SALT-dependent by analyzing SALT-independent observables.

    Key insight: FITPROB (chi2 goodness-of-fit) tells us how well
    SALT fits each SN. If the closure signal correlates with SALT
    fit quality, it might be SALT artifact. If it doesn't, SALT
    is not the driver.
    """
    print("\n" + "="*60)
    print("SALT-Independence Analysis")
    print("="*60)

    results = {}
    df = summary_df.copy()

    # Compute distance residuals
    from scipy.integrate import quad
    H0, Om = 73.04, 0.334
    def dL(z):
        if z <= 0: return 1e-10
        dc, _ = quad(lambda zp: 1.0/np.sqrt(Om*(1+zp)**3 + (1-Om)), 0, z)
        return (1+z) * dc * (299792.458 / H0)

    if 'MU_SH0ES' in df.columns and 'MUMODEL' in df.columns:
        df['mu_resid'] = df['MU_SH0ES'] - df['MUMODEL']
    elif 'MU_SH0ES' in df.columns:
        mu_model = np.array([5*np.log10(dL(z))+25 for z in df['zHD']])
        df['mu_resid'] = df['MU_SH0ES'] - mu_model

    # Color residual (detrended by redshift bin)
    z_edges = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.5]
    if 'c' in df.columns:
        df['c_resid'] = np.nan
        for i in range(len(z_edges)-1):
            mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
            if mask.sum() > 5:
                df.loc[mask, 'c_resid'] = df.loc[mask, 'c'] - df.loc[mask, 'c'].median()

    # Fit quality
    if 'FITCHI2' in df.columns and 'NDOF' in df.columns:
        df['chi2_dof'] = df['FITCHI2'] / df['NDOF'].replace(0, np.nan)
    elif 'FITPROB' in df.columns:
        df['chi2_dof'] = -np.log10(df['FITPROB'].clip(1e-10, 1))

    # --------------------------------------------------------
    # TEST S1: Does closure signal depend on SALT fit quality?
    # --------------------------------------------------------
    print("\n--- Test S1: SALT Fit Quality Stratification ---")
    print("If closure signal is SALT artifact, it should be STRONGER")
    print("in poorly-fit SNe and ABSENT in well-fit SNe.\n")

    if 'FITPROB' in df.columns and 'c_resid' in df.columns and 'mu_resid' in df.columns:
        # Split by fit quality
        good_fit = df[df['FITPROB'] > 0.1].copy()
        bad_fit = df[df['FITPROB'] <= 0.1].copy()

        z_bins_test = [(0.01, 0.3), (0.3, 0.6), (0.6, 2.5)]
        s1_results = []

        for label, subset in [("Good fits (FITPROB>0.1)", good_fit),
                               ("Bad fits (FITPROB≤0.1)", bad_fit)]:
            print(f"  {label}:")
            for z_lo, z_hi in z_bins_test:
                mask = (subset['zHD'] >= z_lo) & (subset['zHD'] < z_hi)
                sub = subset[mask].dropna(subset=['c_resid', 'mu_resid'])
                if len(sub) < 10:
                    print(f"    z=[{z_lo},{z_hi}): N={len(sub)} (too few)")
                    continue
                r, p = stats.spearmanr(sub['c_resid'], sub['mu_resid'])
                print(f"    z=[{z_lo},{z_hi}): r={r:.3f}, p={p:.4f}, N={len(sub)}")
                s1_results.append({
                    'subset': label, 'z_lo': z_lo, 'z_hi': z_hi,
                    'r': r, 'p': p, 'N': len(sub)
                })

        results['S1'] = s1_results

        # Verdict
        good_highz = [r for r in s1_results if 'Good' in r['subset'] and r['z_lo'] >= 0.3]
        bad_highz = [r for r in s1_results if 'Bad' in r['subset'] and r['z_lo'] >= 0.3]
        if good_highz and bad_highz:
            good_signal = np.mean([abs(r['r']) for r in good_highz])
            bad_signal = np.mean([abs(r['r']) for r in bad_highz])
            print(f"\n  Mean |r| at high-z: Good fits = {good_signal:.3f}, Bad fits = {bad_signal:.3f}")
            if good_signal >= bad_signal * 0.7:
                print("  ★ Signal persists in well-fit SNe → NOT a SALT fitting artifact")
            else:
                print("  ⚠ Signal concentrated in bad fits → SALT may contribute")

    # --------------------------------------------------------
    # TEST S2: Tripp-Free Distance Residuals
    # --------------------------------------------------------
    print("\n--- Test S2: Tripp-Free Residuals ---")
    print("Compute distances WITHOUT the Tripp c correction (β*c term).")
    print("If closure signal survives, it's not an artifact of the β*c formula.\n")

    if 'mB' in df.columns and 'x1' in df.columns:
        # Standard Tripp: mu = mB - M + alpha*x1 - beta*c
        # Tripp-free: mu_noc = mB - M + alpha*x1 (no color correction)
        alpha = 0.15  # typical Pantheon+ value
        M = -19.25    # typical

        df['mu_noc'] = df['mB'] - M + alpha * df['x1']

        # Compute Tripp-free residuals
        mu_model = np.array([5*np.log10(dL(z))+25 for z in df['zHD']])
        df['mu_noc_resid'] = df['mu_noc'] - mu_model

        # Detrend mu_noc_resid by redshift bin (remove bulk offset from missing β*c)
        for i in range(len(z_edges)-1):
            mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
            if mask.sum() > 5:
                df.loc[mask, 'mu_noc_resid'] -= df.loc[mask, 'mu_noc_resid'].median()

        # Now test: does c_resid correlate with mu_noc_resid (Tripp-free)?
        s2_results = []
        for z_lo, z_hi in z_bins_test:
            mask = (df['zHD'] >= z_lo) & (df['zHD'] < z_hi)
            sub = df[mask].dropna(subset=['c_resid', 'mu_noc_resid'])
            if len(sub) < 10:
                continue
            r, p = stats.spearmanr(sub['c_resid'], sub['mu_noc_resid'])
            print(f"  z=[{z_lo},{z_hi}): r={r:.3f}, p={p:.4f}, N={len(sub)}")
            s2_results.append({'z_lo': z_lo, 'z_hi': z_hi, 'r': r, 'p': p, 'N': len(sub)})

        results['S2'] = s2_results
    else:
        print("  Missing mB/x1 columns — skipping")

    # --------------------------------------------------------
    # TEST S3: Raw Magnitude Dispersion (Zero SALT Dependence)
    # --------------------------------------------------------
    print("\n--- Test S3: Raw Magnitude Dispersion vs Redshift ---")
    print("Completely SALT-free: just look at scatter in raw peak magnitudes.\n")

    if 'mB' in df.columns:
        # Subtract cosmological model to get raw residuals (no SALT corrections at all)
        df['mB_resid'] = df['mB'] - mu_model - (-19.25)  # subtract M and mu_model

        # Measure dispersion in sliding z windows
        z_centers = np.arange(0.05, 1.5, 0.05)
        z_width = 0.15
        dispersions = []
        for zc in z_centers:
            mask = (df['zHD'] >= zc - z_width/2) & (df['zHD'] < zc + z_width/2)
            sub = df[mask]['mB_resid'].dropna()
            if len(sub) < 10:
                dispersions.append(np.nan)
            else:
                dispersions.append(sub.std())

        dispersions = np.array(dispersions)
        valid = ~np.isnan(dispersions)

        if valid.sum() > 5:
            # Test for change in dispersion around z~0.5-0.8
            low_z = (z_centers < 0.4) & valid
            high_z = (z_centers > 0.6) & valid
            if low_z.sum() > 2 and high_z.sum() > 2:
                disp_low = np.mean(dispersions[low_z])
                disp_high = np.mean(dispersions[high_z])
                ratio = disp_high / disp_low if disp_low > 0 else np.nan
                print(f"  Mean dispersion z<0.4: {disp_low:.4f}")
                print(f"  Mean dispersion z>0.6: {disp_high:.4f}")
                print(f"  Ratio (high/low): {ratio:.3f}")
                if ratio > 1.1:
                    print("  ★ Raw dispersion increases at high-z → consistent with DOF collapse")
                else:
                    print("  Dispersion stable → no raw magnitude evidence")

        results['S3'] = {'z_centers': z_centers.tolist(),
                         'dispersions': dispersions.tolist()}

    # --------------------------------------------------------
    # TEST S4: Color-Distance Coupling Without SALT Color
    # --------------------------------------------------------
    print("\n--- Test S4: Photometric Color Index (SALT-Free) ---")
    print("Use raw photometric color (mB - mV equivalent) instead of SALT c.\n")

    # We can approximate this: SALT c ≈ (B-V) - <B-V> but let's test
    # if the coupling structure changes when we use x1 (stretch) as predictor
    # since stretch is more "temporal" and less "spectral"
    if 'x1' in df.columns and 'mu_resid' in df.columns:
        df['x1_resid'] = np.nan
        for i in range(len(z_edges)-1):
            mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
            if mask.sum() > 5:
                df.loc[mask, 'x1_resid'] = df.loc[mask, 'x1'] - df.loc[mask, 'x1'].median()

        s4_results = []
        print("  Stretch (x1) vs distance residual (should be IMMUNE per closure theory):")
        for z_lo, z_hi in z_bins_test:
            mask = (df['zHD'] >= z_lo) & (df['zHD'] < z_hi)
            sub = df[mask].dropna(subset=['x1_resid', 'mu_resid'])
            if len(sub) < 10:
                continue
            r, p = stats.spearmanr(sub['x1_resid'], sub['mu_resid'])
            print(f"  z=[{z_lo},{z_hi}): r={r:.3f}, p={p:.4f}, N={len(sub)}")
            s4_results.append({'z_lo': z_lo, 'z_hi': z_hi, 'r': r, 'p': p, 'N': len(sub)})

        results['S4'] = s4_results

        # Compare c vs x1 coupling strength
        if 'c_resid' in df.columns:
            print("\n  Color (c) vs distance residual (EXPECTED to strengthen per closure):")
            for z_lo, z_hi in z_bins_test:
                mask = (df['zHD'] >= z_lo) & (df['zHD'] < z_hi)
                sub = df[mask].dropna(subset=['c_resid', 'mu_resid'])
                if len(sub) < 10:
                    continue
                r, p = stats.spearmanr(sub['c_resid'], sub['mu_resid'])
                print(f"  z=[{z_lo},{z_hi}): r={r:.3f}, p={p:.4f}, N={len(sub)}")

    # --------------------------------------------------------
    # TEST S5: SALT2 vs SALT3 Comparison (if available)
    # --------------------------------------------------------
    # This would use SALT3 fits if available — different model, same framework

    return results, df


# ============================================================
# STEP 6: BayeSN Closure Analysis
# ============================================================

def bayesn_closure_analysis(bayesn_df):
    """Analyze BayeSN results for closure signal."""
    if bayesn_df is None or len(bayesn_df) < 20:
        print("  Insufficient BayeSN results for analysis")
        return None

    print("\n" + "="*60)
    print("BayeSN Closure Analysis")
    print("="*60)

    results = {}
    df = bayesn_df.copy()

    # Compute BayeSN residuals
    from scipy.integrate import quad
    H0, Om = 73.04, 0.334
    def dL(z):
        if z <= 0: return 1e-10
        dc, _ = quad(lambda zp: 1.0/np.sqrt(Om*(1+zp)**3 + (1-Om)), 0, z)
        return (1+z) * dc * (299792.458 / H0)

    mu_model = np.array([5*np.log10(dL(z))+25 for z in df['zHD']])
    df['mu_bayesn_resid'] = df['mu_bayesn'] - mu_model

    # Detrend
    z_edges = [0, 0.3, 0.6, 2.5]
    for i in range(len(z_edges)-1):
        mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
        if mask.sum() > 3:
            df.loc[mask, 'mu_bayesn_resid'] -= df.loc[mask, 'mu_bayesn_resid'].median()

    # TEST B1: BayeSN theta vs distance residual at high-z
    print("\n--- Test B1: BayeSN theta (intrinsic) vs BayeSN distance ---")
    z_bins = [(0.01, 0.3), (0.3, 0.6), (0.6, 2.5)]
    for z_lo, z_hi in z_bins:
        mask = (df['zHD'] >= z_lo) & (df['zHD'] < z_hi)
        sub = df[mask].dropna(subset=['theta_bayesn', 'mu_bayesn_resid'])
        if len(sub) < 5:
            continue
        r, p = stats.spearmanr(sub['theta_bayesn'], sub['mu_bayesn_resid'])
        print(f"  z=[{z_lo},{z_hi}): r={r:.3f}, p={p:.4f}, N={len(sub)}")

    # TEST B2: BayeSN AV (dust) vs distance residual at high-z
    print("\n--- Test B2: BayeSN AV (dust) vs BayeSN distance ---")
    for z_lo, z_hi in z_bins:
        mask = (df['zHD'] >= z_lo) & (df['zHD'] < z_hi)
        sub = df[mask].dropna(subset=['AV_bayesn', 'mu_bayesn_resid'])
        if len(sub) < 5:
            continue
        r, p = stats.spearmanr(sub['AV_bayesn'], sub['mu_bayesn_resid'])
        print(f"  z=[{z_lo},{z_hi}): r={r:.3f}, p={p:.4f}, N={len(sub)}")

    # TEST B3: Cross-compare BayeSN vs SALT residuals
    print("\n--- Test B3: BayeSN vs SALT2 distance comparison ---")
    if 'mu_salt' in df.columns:
        df['mu_salt_resid'] = df['mu_salt'] - mu_model
        for i in range(len(z_edges)-1):
            mask = (df['zHD'] >= z_edges[i]) & (df['zHD'] < z_edges[i+1])
            if mask.sum() > 3:
                df.loc[mask, 'mu_salt_resid'] -= df.loc[mask, 'mu_salt_resid'].median()

        sub = df.dropna(subset=['mu_bayesn_resid', 'mu_salt_resid'])
        r, p = stats.spearmanr(sub['mu_bayesn_resid'], sub['mu_salt_resid'])
        print(f"  Overall correlation: r={r:.3f}, p={p:.4f}, N={len(sub)}")

        # If BayeSN and SALT residuals correlate at high-z, the signal is real
        for z_lo, z_hi in z_bins:
            mask = (sub['zHD'] >= z_lo) & (sub['zHD'] < z_hi)
            sub2 = sub[mask]
            if len(sub2) < 5:
                continue
            r2, p2 = stats.spearmanr(sub2['mu_bayesn_resid'], sub2['mu_salt_resid'])
            print(f"  z=[{z_lo},{z_hi}): r={r2:.3f}, p={p2:.4f}, N={len(sub2)}")

    # TEST B4: SALT c vs BayeSN AV correlation evolution
    print("\n--- Test B4: SALT color vs BayeSN dust (AV) ---")
    print("If SALT c is mostly dust, these should correlate strongly at all z.")
    print("If closure is real, the relationship may change at high-z.\n")
    if 'c_salt' in df.columns:
        for z_lo, z_hi in z_bins:
            mask = (df['zHD'] >= z_lo) & (df['zHD'] < z_hi)
            sub = df[mask].dropna(subset=['c_salt', 'AV_bayesn'])
            if len(sub) < 5:
                continue
            r, p = stats.spearmanr(sub['c_salt'], sub['AV_bayesn'])
            print(f"  z=[{z_lo},{z_hi}): r={r:.3f}, p={p:.4f}, N={len(sub)}")

    return results


# ============================================================
# STEP 7: Plots
# ============================================================

def make_plots(summary_df, bayesn_df, salt_results):
    """Generate publication-quality plots."""
    ensure_dir(OUTPUT_DIR)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Plot 1: S1 — Signal in good vs bad SALT fits
    ax = axes[0, 0]
    if 'S1' in salt_results:
        s1 = pd.DataFrame(salt_results['S1'])
        for label, marker, color in [('Good fits', 'o', 'blue'), ('Bad fits', 's', 'red')]:
            sub = s1[s1['subset'].str.contains(label.split()[0])]
            if len(sub) > 0:
                z_mid = (sub['z_lo'] + sub['z_hi']) / 2
                ax.scatter(z_mid, sub['r'].abs(), marker=marker, color=color,
                          s=sub['N']*2, label=label, alpha=0.8, edgecolors='black')
        ax.set_xlabel('Redshift')
        ax.set_ylabel('|Spearman r| (Δμ vs c)')
        ax.set_title('S1: Signal in Good vs Bad SALT Fits')
        ax.legend()
        ax.axhline(0, color='gray', ls='--', alpha=0.5)

    # Plot 2: S2 — Tripp-free residuals
    ax = axes[0, 1]
    if 'S2' in salt_results:
        s2 = pd.DataFrame(salt_results['S2'])
        z_mid = (s2['z_lo'] + s2['z_hi']) / 2
        colors = ['green' if p < 0.05 else 'gray' for p in s2['p']]
        ax.bar(z_mid, s2['r'].abs(), width=0.15, color=colors, alpha=0.7, edgecolor='black')
        ax.set_xlabel('Redshift')
        ax.set_ylabel('|Spearman r| (c_resid vs μ_no_color_resid)')
        ax.set_title('S2: Tripp-Free Distance Residuals')
        ax.axhline(0.1, color='red', ls='--', alpha=0.5, label='weak threshold')
        ax.legend()

    # Plot 3: S3 — Raw dispersion
    ax = axes[1, 0]
    if 'S3' in salt_results:
        zc = salt_results['S3']['z_centers']
        disp = salt_results['S3']['dispersions']
        valid = [i for i, d in enumerate(disp) if not np.isnan(d)]
        ax.plot([zc[i] for i in valid], [disp[i] for i in valid],
                'ko-', markersize=4, alpha=0.7)
        ax.axvspan(0.5, 0.9, alpha=0.15, color='red', label='Closure zone')
        ax.set_xlabel('Redshift')
        ax.set_ylabel('Raw mB dispersion (σ)')
        ax.set_title('S3: Raw Magnitude Scatter (Zero SALT)')
        ax.legend()

    # Plot 4: BayeSN vs SALT or S4 stretch immunity
    ax = axes[1, 1]
    if 'S4' in salt_results:
        s4 = pd.DataFrame(salt_results['S4'])
        z_mid = (s4['z_lo'] + s4['z_hi']) / 2
        ax.bar(z_mid, s4['r'].abs(), width=0.15, color='purple', alpha=0.7, edgecolor='black')
        ax.set_xlabel('Redshift')
        ax.set_ylabel('|Spearman r| (x1_resid vs μ_resid)')
        ax.set_title('S4: Stretch (x1) vs Distance — Should Be Immune')

    plt.suptitle('Closure Theory — SALT-Independent Tests', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "salt_independent_tests.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved: {OUTPUT_DIR}/salt_independent_tests.png")


# ============================================================
# MAIN
# ============================================================

def main():
    ensure_dir(OUTPUT_DIR)

    print("="*60)
    print("CLOSURE THEORY — BayeSN / SALT-INDEPENDENT TEST")
    print("="*60)

    # Load summary data
    print("\n[1/5] Loading Pantheon+ summary data...")
    dat_file = download_summary()
    summary_df = pd.read_csv(dat_file, sep=r'\s+', comment='#')
    print(f"  Loaded {len(summary_df)} SNe")

    # Run SALT-independent analysis (always works)
    print("\n[2/5] Running SALT-independent analysis...")
    salt_results, enriched_df = salt_independent_analysis(summary_df)

    # Try to download light curves and run BayeSN
    bayesn_df = None
    if HAS_BAYESN:
        print("\n[3/5] Downloading Pantheon+ photometry...")
        lc_files = download_photometry()

        if lc_files:
            print(f"\n[4/5] Selecting strategic subsample...")
            selected = select_subsample(summary_df, lc_files, n_per_bin=12)
            print(f"  Total selected: {len(selected)} SNe")

            if selected:
                print(f"\n[5/5] Fitting with BayeSN (GPU)...")
                bayesn_df = fit_bayesn_sample(selected)

                if bayesn_df is not None:
                    bayesn_results = bayesn_closure_analysis(bayesn_df)
                    bayesn_df.to_csv(OUTPUT_DIR / "bayesn_fits.csv", index=False)
                    print(f"\n  Saved: {OUTPUT_DIR}/bayesn_fits.csv")
    else:
        print("\n[3-5/5] BayeSN not available — SALT-independent analysis only")

    # Generate plots
    print("\n[Plots] Generating figures...")
    make_plots(summary_df, bayesn_df, salt_results)

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print("""
Tests completed:
  S1: SALT fit quality stratification
  S2: Tripp-free (no β*c) distance residuals
  S3: Raw magnitude dispersion (zero SALT dependence)
  S4: Stretch immunity check (x1 vs Δμ)
""")
    if bayesn_df is not None:
        print(f"""  B1-B4: BayeSN non-SALT replication ({len(bayesn_df)} SNe fitted)""")

    # Save results
    with open(OUTPUT_DIR / "results.json", 'w') as f:
        # Convert numpy to python types
        def convert(obj):
            if isinstance(obj, (np.integer,)): return int(obj)
            if isinstance(obj, (np.floating,)): return float(obj)
            if isinstance(obj, np.ndarray): return obj.tolist()
            return obj
        json.dump(salt_results, f, indent=2, default=convert)

    print(f"\n  Results saved to {OUTPUT_DIR}/")
    print("  Done!")


if __name__ == '__main__':
    main()
