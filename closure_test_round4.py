#!/usr/bin/env python3
"""
Closure Theory — Round 4: Replication & FRB Cross-Match
=========================================================
Two independent validation tests:

Test R1: DES-5YR Replication
  - Run core closure tests (accumulation, factorization, threshold, rank compression)
    on DES-SN5YR — completely independent catalog, different instruments, different
    selection function, SALT3 instead of SALT2.
  - If frequency fingerprint replicates → not a Pantheon+ artifact.

Test FRB1: FRB Dispersion Measure Cross-Match
  - Build DM_excess sky field from CHIME FRB catalog
  - Sample at SN positions (Pantheon+)
  - Test: does x1 (stretch) correlate with DM? (Grok's angle)
  - Test: does color/distance entanglement correlate with DM?
  - If x1 × DM null → stretch is genuinely immune, frequency-only fingerprint holds
  - If x1 × DM signal → baryonic medium affects temporal features too

Author: Humza Hafeez (Closure Theory)
Date: 2026-02-21
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

warnings.filterwarnings('ignore')

OUTPUT_DIR = Path("results_round4")
DATA_DIR = Path("data")
SEED = 42

def ensure_dir(p):
    p.mkdir(parents=True, exist_ok=True)


# ============================================================
# DES-5YR DATA
# ============================================================

def download_des5yr():
    """Download DES-SN5YR Hubble Diagram + MetaData.
    
    DES-SN5YR uses two files:
    - DES-Dovekie_HD.csv: distances (zHD, MU, MUMODEL, MURES, etc.)
    - DES-Dovekie_Metadata.csv: SN properties (x1, c, mB, HOST_LOGMASS, etc.)
    These are merged on CID.
    """
    ensure_dir(DATA_DIR)
    merged_file = DATA_DIR / "DES-SN5YR_merged.csv"
    
    if not merged_file.exists():
        hd_file = DATA_DIR / "DES-Dovekie_HD.csv"
        meta_file = DATA_DIR / "DES-Dovekie_Metadata.csv"
        
        # Download HD file
        if not hd_file.exists():
            url_hd = "https://raw.githubusercontent.com/des-science/DES-SN5YR/main/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv"
            print(f"[*] Downloading DES-SN5YR HD file...")
            result = subprocess.run(["curl", "-fSL", "-o", str(hd_file), url_hd], capture_output=True)
            if result.returncode != 0:
                print(f"[!] HD download failed: {result.stderr.decode()}")
                return None
        
        # Download Metadata file
        if not meta_file.exists():
            url_meta = "https://raw.githubusercontent.com/des-science/DES-SN5YR/main/4_DISTANCES_COVMAT/DES-Dovekie_Metadata.csv"
            print(f"[*] Downloading DES-SN5YR Metadata file...")
            result = subprocess.run(["curl", "-fSL", "-o", str(meta_file), url_meta], capture_output=True)
            if result.returncode != 0:
                print(f"[!] Metadata download failed: {result.stderr.decode()}")
                return None
        
        # Parse and merge
        try:
            # HD file has comment lines starting with #
            df_hd = pd.read_csv(hd_file, comment='#', sep=r'\s+')
            print(f"    HD file: {len(df_hd)} rows, columns: {list(df_hd.columns)[:10]}...")
            
            # Metadata file has VARNAMES header
            df_meta = pd.read_csv(meta_file, comment='#', sep=r'\s+')
            # Drop the VARNAMES row if present
            if 'VARNAMES:' in str(df_meta.iloc[0].values):
                # The first column might be 'VARNAMES:' or similar
                cols = df_meta.columns.tolist()
                if cols[0] == 'VARNAMES:':
                    cols = cols[1:]  # shift
                df_meta.columns = cols + ['_extra'] * (len(df_meta.columns) - len(cols))
            
            # Try to find CID column for merge
            cid_col_hd = None
            for c in ['CID', 'cid', 'SNID', 'SN']:
                if c in df_hd.columns:
                    cid_col_hd = c
                    break
            
            cid_col_meta = None
            for c in ['CID', 'cid', 'SNID', 'SN']:
                if c in df_meta.columns:
                    cid_col_meta = c
                    break
            
            if cid_col_hd and cid_col_meta:
                df = pd.merge(df_hd, df_meta, left_on=cid_col_hd, right_on=cid_col_meta, 
                             how='inner', suffixes=('', '_meta'))
            else:
                # If same length, just concat
                if len(df_hd) == len(df_meta):
                    df = pd.concat([df_hd.reset_index(drop=True), df_meta.reset_index(drop=True)], axis=1)
                else:
                    df = df_meta  # Metadata has everything
            
            df.to_csv(merged_file, index=False)
            print(f"    Merged: {len(df)} rows")
        except Exception as e:
            print(f"[!] Merge failed: {e}")
            # Fall back to just metadata which has most columns
            if meta_file.exists():
                merged_file = meta_file
            else:
                return None
    
    return merged_file

def load_des5yr():
    """Load DES-5YR data and compute derived quantities.
    
    DES-SN5YR Metadata file uses SNANA format with 'VARNAMES:' header
    and 'SN:' data rows. Columns include: CID, IDSURVEY, zHD, x1, c, mB,
    MU, MUMODEL, MURES, FITCHI2, NDOF, HOST_LOGMASS, etc.
    """
    des_file = download_des5yr()
    if not des_file or not des_file.exists():
        return None

    # Try standard CSV first
    try:
        df = pd.read_csv(des_file)
        if len(df) > 10 and 'zHD' in df.columns:
            pass  # Good, standard CSV worked
        else:
            raise ValueError("Standard CSV parse didn't work")
    except:
        # SNANA format: parse VARNAMES header, then SN: rows
        try:
            with open(des_file) as f:
                lines = f.readlines()
            
            # Find VARNAMES line
            varnames = None
            data_lines = []
            for line in lines:
                line = line.strip()
                if line.startswith('VARNAMES:'):
                    varnames = line.replace('VARNAMES:', '').split()
                elif line.startswith('SN:'):
                    data_lines.append(line.replace('SN:', '').split())
            
            if varnames and data_lines:
                df = pd.DataFrame(data_lines, columns=varnames)
                # Convert numeric columns
                for col in df.columns:
                    try:
                        df[col] = pd.to_numeric(df[col])
                    except:
                        pass
            else:
                # Last resort: space-separated with comment lines
                df = pd.read_csv(des_file, comment='#', sep=r'\s+')
        except Exception as e:
            print(f"[!] Failed to parse DES data: {e}")
            return None

    print(f"[*] DES-5YR: Loaded {len(df)} SNe")
    print(f"    Columns: {sorted(df.columns.tolist())[:20]}...")
    
    if 'zHD' not in df.columns:
        print(f"[!] No zHD column found")
        return None
    
    print(f"    z range: {df['zHD'].min():.3f} - {df['zHD'].max():.3f}")

    # Distance residual
    if 'MURES' in df.columns:
        df['mu_resid'] = df['MURES']
    elif 'MU' in df.columns and 'MUMODEL' in df.columns:
        df['mu_resid'] = df['MU'] - df['MUMODEL']
    elif 'MU' in df.columns:
        from scipy.integrate import quad
        H0, Om = 67.36, 0.315
        def dL(z):
            if z <= 0: return 1e-10
            dc, _ = quad(lambda zp: 1.0/np.sqrt(Om*(1+zp)**3 + (1-Om)), 0, z)
            return (1+z) * dc * (299792.458 / H0)
        mu_model = np.array([5*np.log10(dL(z))+25 for z in df['zHD']])
        df['mu_resid'] = df['MU'] - mu_model

    # SALT3 color and stretch (DES-SN5YR uses standard SALT column names)
    # c and x1 should already be present

    # Color residual
    if 'c' in df.columns:
        z_edges = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.5]
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

    n_valid = {col: df[col].notna().sum() for col in ['mu_resid', 'c', 'x1', 'chi2_dof'] if col in df.columns}
    print(f"    Valid counts: {n_valid}")

    return df


# ============================================================
# TEST R1: DES-5YR REPLICATION
# ============================================================

def test_R1_des_replication(df_des):
    """Run core closure tests on DES-5YR independently."""
    print(f"\n{'='*60}")
    print("TEST R1: DES-5YR REPLICATION")
    print(f"{'='*60}")

    if df_des is None:
        print("  [!] No DES data available")
        return {'error': 'No DES data'}

    results = {}

    # Check available columns
    has_c = 'c' in df_des.columns and df_des['c'].notna().sum() > 50
    has_mu = 'mu_resid' in df_des.columns and df_des['mu_resid'].notna().sum() > 50
    has_x1 = 'x1' in df_des.columns and df_des['x1'].notna().sum() > 50

    if not has_c or not has_mu:
        # Try to find the right columns
        print(f"  Available columns: {sorted(df_des.columns.tolist())}")
        print(f"  [!] Missing c or mu_resid — checking alternatives...")

        # Print first few rows to debug
        print(f"\n  Sample data (first 3 rows):")
        print(df_des.head(3).to_string())
        results['error'] = 'Missing required columns (c or mu_resid)'
        results['available_columns'] = sorted(df_des.columns.tolist())
        return results

    print(f"\n  DES-5YR: N={len(df_des)}, z range [{df_des['zHD'].min():.3f}, {df_des['zHD'].max():.3f}]")

    fig, axes = plt.subplots(2, 3, figsize=(18, 11))

    # ---- Sub-test R1a: Accumulation (= Test 4) ----
    z_edges = [0.0, 0.15, 0.3, 0.5, 0.7, 1.0, 2.5]
    z_mids = [(z_edges[i]+z_edges[i+1])/2 for i in range(len(z_edges)-1)]

    rs, ps, ns = [], [], []
    for i in range(len(z_edges)-1):
        mask = (df_des['zHD'] >= z_edges[i]) & (df_des['zHD'] < z_edges[i+1])
        mask &= df_des[['c', 'mu_resid']].notna().all(axis=1)
        n = mask.sum()
        ns.append(n)
        if n >= 10:
            r, p = stats.pearsonr(df_des.loc[mask, 'c'], df_des.loc[mask, 'mu_resid'])
            rs.append(r)
            ps.append(p)
            sig = "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.1 else ""
            print(f"  R1a z=[{z_edges[i]:.2f},{z_edges[i+1]:.2f}): r={r:.4f}, p={p:.2e}, N={n} {sig}")
        else:
            rs.append(np.nan)
            ps.append(np.nan)

    results['R1a_accumulation'] = [
        {'z_bin': f'{z_edges[i]:.2f}-{z_edges[i+1]:.2f}', 'r': round(r, 4) if np.isfinite(r) else None,
         'p': round(p, 6) if np.isfinite(p) else None, 'N': int(n)}
        for i, (r, p, n) in enumerate(zip(rs, ps, ns))
    ]

    ax = axes[0, 0]
    valid = [np.isfinite(r) for r in rs]
    zv = [z for z, v in zip(z_mids, valid) if v]
    rv = [r for r, v in zip(rs, valid) if v]
    ax.plot(zv, rv, 'o-', color='crimson', ms=8, lw=2)
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('r(c, Δμ)')
    ax.set_title('R1a: Accumulation (DES-5YR)')
    ax.axhline(0, color='gray', ls='--', alpha=0.5)

    # ---- Sub-test R1b: Factorization (= Test 6 simplified) ----
    if has_x1:
        obs_pairs = [('c', 'mu_resid', 'c vs Δμ'), ('x1', 'mu_resid', 'x1 vs Δμ'),
                     ('c', 'x1', 'c vs x1')]
        if 'chi2_dof' in df_des.columns:
            obs_pairs.append(('c', 'chi2_dof', 'c vs χ²'))

        ax = axes[0, 1]
        results['R1b_factorization'] = {}

        for col_a, col_b, label in obs_pairs:
            if col_a not in df_des.columns or col_b not in df_des.columns:
                continue
            pair_rs = []
            pair_zs = []
            for i in range(len(z_edges)-1):
                mask = (df_des['zHD'] >= z_edges[i]) & (df_des['zHD'] < z_edges[i+1])
                mask &= df_des[[col_a, col_b]].notna().all(axis=1)
                if mask.sum() >= 10:
                    r, _ = stats.pearsonr(df_des.loc[mask, col_a], df_des.loc[mask, col_b])
                    pair_rs.append(abs(r))
                    pair_zs.append(z_mids[i])

            if len(pair_zs) >= 3:
                slope, _, rr, pp, _ = stats.linregress(pair_zs, pair_rs)
                ax.plot(pair_zs, pair_rs, 'o-', ms=6, lw=1.5, label=f'{label} (slope={slope:.3f})')
                results['R1b_factorization'][label] = {
                    'slope': round(slope, 4), 'trend_r': round(rr, 4), 'trend_p': round(pp, 6)
                }
                sig = "***" if pp < 0.01 else "**" if pp < 0.05 else ""
                print(f"  R1b {label}: slope={slope:.4f}, p={pp:.4f} {sig}")

        ax.set_xlabel('Redshift z')
        ax.set_ylabel('|Pearson r|')
        ax.set_title('R1b: Factorization Collapse (DES-5YR)')
        ax.legend(fontsize=8)

    # ---- Sub-test R1c: Threshold (= Test 7) ----
    z_fine = np.arange(0.0, 1.3, 0.05)
    z_fine_mids = (z_fine[:-1] + z_fine[1:]) / 2
    fine_rs = []
    for i in range(len(z_fine)-1):
        mask = (df_des['zHD'] >= z_fine[i]) & (df_des['zHD'] < z_fine[i+1])
        mask &= df_des[['c', 'mu_resid']].notna().all(axis=1)
        if mask.sum() >= 10:
            r, _ = stats.pearsonr(df_des.loc[mask, 'c'], df_des.loc[mask, 'mu_resid'])
            fine_rs.append(abs(r))
        else:
            fine_rs.append(np.nan)

    ax = axes[0, 2]
    valid_fine = [(z, r) for z, r in zip(z_fine_mids, fine_rs) if np.isfinite(r)]
    if len(valid_fine) >= 5:
        zf = np.array([v[0] for v in valid_fine])
        rf = np.array([v[1] for v in valid_fine])

        def sigmoid(z, r_max, z0, k, r_min):
            return r_min + (r_max - r_min) / (1 + np.exp(-k * (z - z0)))

        try:
            popt, _ = curve_fit(sigmoid, zf, rf, p0=[0.7, 0.5, 5, 0.1],
                                bounds=([0, 0, 0.1, 0], [1, 2, 50, 0.5]), maxfev=10000)
            z_fit = np.linspace(0, 1.3, 100)
            ax.plot(z_fit, sigmoid(z_fit, *popt), '-', color='crimson', lw=2,
                    label=f'Sigmoid: z₀={popt[1]:.2f}')
            results['R1c_threshold'] = {
                'z_threshold': round(popt[1], 3), 'steepness': round(popt[2], 2)
            }
            print(f"  R1c Threshold: z₀={popt[1]:.3f}, steepness={popt[2]:.1f}")
        except Exception as e:
            print(f"  R1c Sigmoid fit failed: {e}")
            results['R1c_threshold'] = {'error': str(e)}

        ax.plot(zf, rf, 'o', color='steelblue', ms=6)
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('|r(c, Δμ)|')
    ax.set_title('R1c: Threshold Detection (DES-5YR)')
    if 'R1c_threshold' in results and 'z_threshold' in results.get('R1c_threshold', {}):
        ax.legend()

    # ---- Sub-test R1d: Information Compression (= Test 8) ----
    obs_cols = [c for c in ['mu_resid', 'c', 'x1', 'chi2_dof']
                if c in df_des.columns and df_des[c].notna().sum() > 50]

    z_rank_edges = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.5]
    rank_zs, rank_vals = [], []

    for i in range(len(z_rank_edges)-1):
        mask = (df_des['zHD'] >= z_rank_edges[i]) & (df_des['zHD'] < z_rank_edges[i+1])
        mask &= df_des[obs_cols].notna().all(axis=1)
        n = mask.sum()
        if n < 15:
            continue

        data = df_des.loc[mask, obs_cols].values
        data_std = (data - np.mean(data, axis=0)) / (np.std(data, axis=0) + 1e-10)
        eigenvalues = np.maximum(np.linalg.eigvalsh(np.cov(data_std.T)), 0)[::-1]
        ev_norm = eigenvalues / (eigenvalues.sum() + 1e-10)
        ev_norm = ev_norm[ev_norm > 1e-10]
        eff_rank = np.exp(-np.sum(ev_norm * np.log(ev_norm)))

        rank_zs.append((z_rank_edges[i] + z_rank_edges[i+1])/2)
        rank_vals.append(eff_rank)
        print(f"  R1d z=[{z_rank_edges[i]:.1f},{z_rank_edges[i+1]:.1f}): rank={eff_rank:.2f}, N={n}")

    ax = axes[1, 0]
    if len(rank_zs) >= 3:
        ax.plot(rank_zs, rank_vals, 'o-', color='crimson', ms=8, lw=2)
        slope, _, rr, pp, _ = stats.linregress(rank_zs, rank_vals)
        results['R1d_compression'] = {
            'slope': round(slope, 4), 'r': round(rr, 4), 'p': round(pp, 6),
            'rank_range': [round(min(rank_vals), 3), round(max(rank_vals), 3)]
        }
        print(f"  R1d Rank trend: slope={slope:.4f}, r={rr:.4f}, p={pp:.4f}")
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('Effective Rank')
    ax.set_title('R1d: Information Compression (DES-5YR)')
    ax.axhline(len(obs_cols), color='gray', ls=':', alpha=0.5)

    # ---- Sub-test R1e: Conditional Independence (= Test A1) ----
    ax = axes[1, 1]
    high_z = df_des['zHD'] >= 0.3  # DES has fewer high-z, use lower threshold

    if has_c and has_mu:
        # Raw
        mask_raw = high_z & df_des[['c', 'mu_resid']].notna().all(axis=1)
        r_raw, p_raw = stats.pearsonr(df_des.loc[mask_raw, 'c'], df_des.loc[mask_raw, 'mu_resid'])

        # Controlled (what's available)
        ctrl_cols = []
        if has_x1: ctrl_cols.append('x1')
        if 'HOST_LOGMASS' in df_des.columns: ctrl_cols.append('HOST_LOGMASS')
        if 'HOSTGAL_LOGMASS' in df_des.columns: ctrl_cols.append('HOSTGAL_LOGMASS')
        if 'chi2_dof' in df_des.columns: ctrl_cols.append('chi2_dof')

        all_cols = ['c', 'mu_resid'] + ctrl_cols
        mask_ctrl = high_z & df_des[all_cols].notna().all(axis=1)
        n_ctrl = mask_ctrl.sum()

        if n_ctrl >= 30 and ctrl_cols:
            X = np.column_stack([np.ones(n_ctrl), df_des.loc[mask_ctrl, ctrl_cols].values])
            coef_c = np.linalg.lstsq(X, df_des.loc[mask_ctrl, 'c'].values, rcond=None)[0]
            coef_mu = np.linalg.lstsq(X, df_des.loc[mask_ctrl, 'mu_resid'].values, rcond=None)[0]
            c_r = df_des.loc[mask_ctrl, 'c'].values - X @ coef_c
            mu_r = df_des.loc[mask_ctrl, 'mu_resid'].values - X @ coef_mu
            r_ctrl, p_ctrl = stats.pearsonr(c_r, mu_r)
        else:
            r_ctrl, p_ctrl = r_raw, p_raw

        results['R1e_conditional'] = {
            'raw_r': round(abs(r_raw), 4), 'raw_p': round(p_raw, 6),
            'controlled_r': round(abs(r_ctrl), 4), 'controlled_p': round(p_ctrl, 6),
            'survival_pct': round(abs(r_ctrl)/max(abs(r_raw), 1e-10)*100, 1),
            'N': int(mask_raw.sum()), 'controls': ctrl_cols
        }
        print(f"  R1e DES z>0.3: raw |r|={abs(r_raw):.3f} (p={p_raw:.2e}) → controlled |r|={abs(r_ctrl):.3f} (p={p_ctrl:.2e})")

        bars = ax.bar(['Raw', 'Controlled'], [abs(r_raw), abs(r_ctrl)],
                       color=['gray', 'crimson'], alpha=0.8)
        ax.set_ylabel('|r(c, Δμ)|')
        ax.set_title(f'R1e: Conditional Independence\n(DES z>0.3, N={mask_raw.sum()})')

    # ---- Summary panel ----
    ax = axes[1, 2]
    ax.axis('off')
    summary_text = "DES-5YR REPLICATION SUMMARY\n" + "="*35 + "\n\n"

    if 'R1a_accumulation' in results:
        high_z_r = [d['r'] for d in results['R1a_accumulation'] if d['r'] is not None and d['N'] > 0]
        if high_z_r:
            summary_text += f"Accumulation: r goes from {high_z_r[0]:.3f} → {high_z_r[-1]:.3f}\n"

    if 'R1c_threshold' in results and 'z_threshold' in results.get('R1c_threshold', {}):
        summary_text += f"Threshold: z₀ = {results['R1c_threshold']['z_threshold']:.2f}\n"

    if 'R1d_compression' in results:
        d = results['R1d_compression']
        summary_text += f"Rank compression: p = {d['p']:.4f}\n"
        summary_text += f"  Range: {d['rank_range'][0]:.2f} → {d['rank_range'][1]:.2f}\n"

    if 'R1e_conditional' in results:
        d = results['R1e_conditional']
        summary_text += f"Conditional: {d['survival_pct']:.0f}% survival\n"

    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    fig.suptitle('Test R1: DES-5YR Replication — Do Closure Signatures Appear in Independent Data?',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test_R1_des_replication.png', dpi=150, bbox_inches='tight')
    plt.close()

    return results


# ============================================================
# FRB DATA
# ============================================================

def download_chime_frb():
    """Download CHIME/FRB Catalog 1 CSV."""
    ensure_dir(DATA_DIR)
    frb_file = DATA_DIR / "chimefrb_catalog1.csv"
    if not frb_file.exists():
        url = "https://raw.githubusercontent.com/chime-frb-open-data/CHIME-FRB-Catalog-1/main/catalog1.csv"
        print("[*] Downloading CHIME/FRB Catalog 1...")
        result = subprocess.run(["curl", "-fSL", "-o", str(frb_file), url], capture_output=True)
        if result.returncode != 0:
            print(f"[!] Download failed, trying alternative...")
            # Try the open data site
            url2 = "https://www.chime-frb.ca/static/data/catalog1.csv"
            result = subprocess.run(["curl", "-fSL", "-o", str(frb_file), url2], capture_output=True)
            if result.returncode != 0:
                print(f"[!] FRB download failed")
                return None
    return frb_file

def load_chime_frb():
    """Load CHIME FRB catalog."""
    frb_file = download_chime_frb()
    if not frb_file or not frb_file.exists():
        return None

    df = pd.read_csv(frb_file)
    print(f"[*] CHIME FRB Catalog: {len(df)} FRBs")
    print(f"    Columns: {sorted(df.columns.tolist())[:20]}...")

    # We need: RA, Dec, DM, DM_excess (DM - DM_MW)
    # CHIME provides dm_exc_ne2001 and dm_exc_ymw16
    for dm_col in ['dm_exc_ne2001', 'dm_exc_ymw16', 'dm_fitb', 'bonsai_dm']:
        if dm_col in df.columns:
            df['DM_excess'] = df[dm_col]
            print(f"    Using {dm_col} as DM_excess")
            break

    # RA/Dec
    if 'ra' in df.columns and 'dec' in df.columns:
        df['RA_frb'] = df['ra']
        df['DEC_frb'] = df['dec']

    n_valid = df[['RA_frb', 'DEC_frb', 'DM_excess']].notna().all(axis=1).sum()
    print(f"    Valid (RA, Dec, DM_excess): {n_valid}")

    return df


def build_dm_sky_field(df_frb, df_sn, n_neighbors=5):
    """Assign DM_excess to each SN position using inverse-distance-weighted average of nearest FRBs.

    Simple IDW interpolation — not a full GP, but fast and sufficient for a first pass.
    """
    from scipy.spatial import cKDTree

    # Filter FRBs with valid data
    frb_mask = df_frb[['RA_frb', 'DEC_frb', 'DM_excess']].notna().all(axis=1)
    frb_mask &= df_frb['DM_excess'] > 0  # Physical: positive DM excess
    frb = df_frb[frb_mask].copy()

    print(f"[*] Building DM sky field from {len(frb)} FRBs")

    # Convert to Cartesian for 3D tree (on unit sphere)
    def radec_to_xyz(ra, dec):
        ra_rad = np.radians(ra)
        dec_rad = np.radians(dec)
        x = np.cos(dec_rad) * np.cos(ra_rad)
        y = np.cos(dec_rad) * np.sin(ra_rad)
        z = np.sin(dec_rad)
        return np.column_stack([x, y, z])

    frb_xyz = radec_to_xyz(frb['RA_frb'].values, frb['DEC_frb'].values)
    tree = cKDTree(frb_xyz)

    # Query for each SN
    sn_ra = df_sn['RA'].values if 'RA' in df_sn.columns else np.zeros(len(df_sn))
    sn_dec = df_sn['DEC'].values if 'DEC' in df_sn.columns else np.zeros(len(df_sn))
    sn_xyz = radec_to_xyz(sn_ra, sn_dec)

    k = min(n_neighbors, len(frb))
    distances, indices = tree.query(sn_xyz, k=k)

    # IDW: weighted average of nearest FRB DM values
    dm_values = frb['DM_excess'].values
    weights = 1.0 / (distances + 1e-10)  # Inverse distance
    dm_weighted = np.sum(weights * dm_values[indices], axis=1) / np.sum(weights, axis=1)

    # Also track angular distance to nearest FRB (in degrees)
    nearest_dist_deg = np.degrees(2 * np.arcsin(distances[:, 0] / 2))

    df_sn['DM_excess_idw'] = dm_weighted
    df_sn['nearest_frb_deg'] = nearest_dist_deg

    print(f"    DM_excess_idw range: {dm_weighted.min():.1f} - {dm_weighted.max():.1f} pc/cm³")
    print(f"    Nearest FRB: median {np.median(nearest_dist_deg):.1f}° away")

    return df_sn


# ============================================================
# TEST FRB1: FRB DM CROSS-MATCH
# ============================================================

def test_FRB1_dm_crossmatch(df_sn, df_frb):
    """Test SN observables against FRB-derived DM excess field.

    Key tests:
    1. x1 (stretch) vs DM — Grok's angle: does baryonic medium affect temporal features?
    2. c (color) vs DM — does baryonic medium affect spectral features?
    3. mu_resid vs DM — does baryonic medium affect inferred distance?
    4. |r(c, mu_resid)| in high-DM vs low-DM environments
    """
    print(f"\n{'='*60}")
    print("TEST FRB1: FRB Dispersion Measure Cross-Match")
    print(f"{'='*60}")

    if df_frb is None:
        print("  [!] No FRB data available")
        return {'error': 'No FRB data'}

    results = {}

    # Build DM field
    df_sn = build_dm_sky_field(df_frb, df_sn)

    fig, axes = plt.subplots(2, 3, figsize=(18, 11))

    # ---- Test FRB1a: Direct correlations (like Tests 1-3 but with DM) ----
    ax = axes[0, 0]
    observables = [('mu_resid', 'Δμ'), ('c', 'color c'), ('x1', 'stretch x1')]
    if 'chi2_dof' in df_sn.columns:
        observables.append(('chi2_dof', 'χ²/dof'))

    colors_obs = ['crimson', 'steelblue', 'green', 'orange']
    results['FRB1a_correlations'] = {}

    for i, (col, label) in enumerate(observables):
        if col not in df_sn.columns:
            continue
        mask = df_sn[[col, 'DM_excess_idw']].notna().all(axis=1)
        if mask.sum() < 30:
            continue

        r, p = stats.pearsonr(df_sn.loc[mask, 'DM_excess_idw'], df_sn.loc[mask, col])
        rho, p_rho = stats.spearmanr(df_sn.loc[mask, 'DM_excess_idw'], df_sn.loc[mask, col])
        results['FRB1a_correlations'][col] = {
            'pearson_r': round(r, 4), 'pearson_p': round(p, 6),
            'spearman_rho': round(rho, 4), 'spearman_p': round(p_rho, 6),
            'N': int(mask.sum())
        }
        sig = "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.1 else ""
        print(f"  FRB1a {label:15s} vs DM: r={r:.4f} (p={p:.2e}), ρ={rho:.4f} {sig}")

    # Bar chart of correlations
    labels_bar = [label for col, label in observables if col in results.get('FRB1a_correlations', {})]
    rs_bar = [results['FRB1a_correlations'][col]['pearson_r'] for col, _ in observables if col in results.get('FRB1a_correlations', {})]
    colors_bar = [colors_obs[i] for i, (col, _) in enumerate(observables) if col in results.get('FRB1a_correlations', {})]

    if labels_bar:
        ax.bar(range(len(labels_bar)), [abs(r) for r in rs_bar], color=colors_bar, alpha=0.8)
        ax.set_xticks(range(len(labels_bar)))
        ax.set_xticklabels(labels_bar, fontsize=10)
        ax.set_ylabel('|Pearson r| with DM_excess')
        ax.set_title('FRB1a: Observable vs DM_excess')
        ax.axhline(0.05, color='gray', ls=':', alpha=0.5, label='noise floor')

    # ---- Test FRB1b: x1 vs DM scatter plot (Grok's key test) ----
    ax = axes[0, 1]
    if 'x1' in df_sn.columns:
        mask = df_sn[['x1', 'DM_excess_idw']].notna().all(axis=1)
        if mask.sum() > 20:
            ax.scatter(df_sn.loc[mask, 'DM_excess_idw'], df_sn.loc[mask, 'x1'],
                       alpha=0.2, s=10, color='green', rasterized=True)
            r_x1, p_x1 = stats.pearsonr(df_sn.loc[mask, 'DM_excess_idw'], df_sn.loc[mask, 'x1'])
            ax.set_xlabel('DM_excess (pc/cm³)')
            ax.set_ylabel('Stretch x1')
            ax.set_title(f"FRB1b: Grok's Test — x1 vs DM\nr={r_x1:.4f}, p={p_x1:.2e}")
            results['FRB1b_grok_x1_dm'] = {'r': round(r_x1, 4), 'p': round(p_x1, 6)}

    # ---- Test FRB1c: c vs DM scatter ----
    ax = axes[0, 2]
    if 'c' in df_sn.columns:
        mask = df_sn[['c', 'DM_excess_idw']].notna().all(axis=1)
        if mask.sum() > 20:
            ax.scatter(df_sn.loc[mask, 'DM_excess_idw'], df_sn.loc[mask, 'c'],
                       alpha=0.2, s=10, color='steelblue', rasterized=True)
            r_c, p_c = stats.pearsonr(df_sn.loc[mask, 'DM_excess_idw'], df_sn.loc[mask, 'c'])
            ax.set_xlabel('DM_excess (pc/cm³)')
            ax.set_ylabel('SALT2 color c')
            ax.set_title(f'FRB1c: Color vs DM\nr={r_c:.4f}, p={p_c:.2e}')

    # ---- Test FRB1d: Does entanglement differ in high-DM vs low-DM environments? ----
    ax = axes[1, 0]
    dm_median = df_sn['DM_excess_idw'].median()
    low_dm = df_sn['DM_excess_idw'] <= dm_median
    high_dm = df_sn['DM_excess_idw'] > dm_median

    z_edges_dm = [0.0, 0.15, 0.3, 0.5, 0.7, 1.0, 2.5]
    z_mids_dm = [(z_edges_dm[i]+z_edges_dm[i+1])/2 for i in range(len(z_edges_dm)-1)]

    for split_name, split_mask, color in [('Low DM', low_dm, 'steelblue'), ('High DM', high_dm, 'crimson')]:
        split_rs = []
        split_zs = []
        for i in range(len(z_edges_dm)-1):
            zmask = (df_sn['zHD'] >= z_edges_dm[i]) & (df_sn['zHD'] < z_edges_dm[i+1])
            mask = zmask & split_mask & df_sn[['c', 'mu_resid']].notna().all(axis=1)
            if mask.sum() >= 10:
                r, _ = stats.pearsonr(df_sn.loc[mask, 'c'], df_sn.loc[mask, 'mu_resid'])
                split_rs.append(abs(r))
                split_zs.append(z_mids_dm[i])

        if split_zs:
            ax.plot(split_zs, split_rs, 'o-', color=color, ms=8, lw=2, label=split_name)

    ax.set_xlabel('Redshift z')
    ax.set_ylabel('|r(c, Δμ)|')
    ax.set_title('FRB1d: Entanglement in High vs Low DM')
    ax.legend()

    # ---- Test FRB1e: x1 vs DM in z-bins (does temporal affect grow with z?) ----
    ax = axes[1, 1]
    if 'x1' in df_sn.columns:
        x1_dm_rs = []
        x1_dm_zs = []
        for i in range(len(z_edges_dm)-1):
            zmask = (df_sn['zHD'] >= z_edges_dm[i]) & (df_sn['zHD'] < z_edges_dm[i+1])
            mask = zmask & df_sn[['x1', 'DM_excess_idw']].notna().all(axis=1)
            if mask.sum() >= 15:
                r, _ = stats.pearsonr(df_sn.loc[mask, 'DM_excess_idw'], df_sn.loc[mask, 'x1'])
                x1_dm_rs.append(abs(r))
                x1_dm_zs.append(z_mids_dm[i])

        if x1_dm_zs:
            ax.plot(x1_dm_zs, x1_dm_rs, 'o-', color='green', ms=8, lw=2, label='|r(x1, DM)|')
            if len(x1_dm_zs) >= 3:
                slope, _, _, pp, _ = stats.linregress(x1_dm_zs, x1_dm_rs)
                results['FRB1e_x1_dm_trend'] = {'slope': round(slope, 4), 'p': round(pp, 6)}
                print(f"  FRB1e x1-DM trend with z: slope={slope:.4f}, p={pp:.4f}")

    ax.set_xlabel('Redshift z')
    ax.set_ylabel('|r(x1, DM)|')
    ax.set_title('FRB1e: Does x1-DM Coupling Grow with z?')
    ax.axhline(0.05, color='gray', ls=':', alpha=0.5)

    # ---- Summary ----
    ax = axes[1, 2]
    ax.axis('off')
    summary = "FRB CROSS-MATCH SUMMARY\n" + "="*30 + "\n\n"
    if 'FRB1a_correlations' in results:
        for col, data in results['FRB1a_correlations'].items():
            sig = "SIG" if data['pearson_p'] < 0.05 else "null"
            summary += f"{col:12s} vs DM: r={data['pearson_r']:.3f} [{sig}]\n"
    if 'FRB1b_grok_x1_dm' in results:
        d = results['FRB1b_grok_x1_dm']
        summary += f"\nGrok's test (x1 vs DM):\n  r={d['r']:.4f}, p={d['p']:.2e}\n"
    ax.text(0.1, 0.9, summary, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    fig.suptitle("Test FRB1: FRB Dispersion Measure Cross-Match\n"
                 '"Does the baryonic medium (traced by free electrons) affect SN observables?"',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'test_FRB1_dm_crossmatch.png', dpi=150, bbox_inches='tight')
    plt.close()

    return results


# ============================================================
# MAIN
# ============================================================

def main():
    ensure_dir(OUTPUT_DIR)
    np.random.seed(SEED)

    print("=" * 60)
    print("CLOSURE THEORY — ROUND 4: REPLICATION & FRB CROSS-MATCH")
    print("Author: Humza Hafeez")
    print("=" * 60)

    all_results = {}

    # R1: DES-5YR Replication
    df_des = load_des5yr()
    all_results['R1_des_replication'] = test_R1_des_replication(df_des)

    # FRB1: Load Pantheon+ (for SN positions) and CHIME FRB
    # CHIME FRB catalog requires browser-based download (SPA site)
    # Skip in CI if data not pre-staged
    frb_file = DATA_DIR / "chimefrb_catalog1.csv"
    if frb_file.exists():
        from closure_test_round3 import load_data as load_pantheon
        df_pantheon = load_pantheon()
        df_frb = load_chime_frb()
        all_results['FRB1_dm_crossmatch'] = test_FRB1_dm_crossmatch(df_pantheon, df_frb)
    else:
        print("\n[!] CHIME FRB data not available (requires manual download from chime-frb.ca)")
        print("    Skipping FRB1 test. Pre-stage data/chimefrb_catalog1.csv to enable.")
        all_results['FRB1_dm_crossmatch'] = {'status': 'skipped', 'reason': 'FRB data not pre-staged'}

    # Save
    results_file = OUTPUT_DIR / 'round4_results.json'
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n[*] Results saved to {results_file}")

    print(f"\n{'='*60}")
    print("ROUND 4 SUMMARY")
    print(f"{'='*60}")
    print("R1:   DES-5YR Replication — do closure signatures appear in independent data?")
    print("FRB1: FRB DM Cross-Match — does baryonic medium (electron column) affect SN observables?")
    print(f"\nOutput: {OUTPUT_DIR}/")


if __name__ == '__main__':
    main()
