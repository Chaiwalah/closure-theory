#!/usr/bin/env python3
"""
Closure Theory — FRB DM Cross-Match Test v3
=============================================
Fixes from GPT review:
1. Resolution curve: rho vs max_sep cutoff (the decisive test)
2. Galactic latitude / E(B-V) / declination controls
3. Footprint-preserving permutation (shuffle DM within dec bands)
4. Restrict to |b| > 20° to reduce MW foreground contamination
5. Partial correlation controlling for confounds
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

RESULTS_DIR = Path("results_frb_v3")
RESULTS_DIR.mkdir(exist_ok=True)


def load_pantheon():
    """Load Pantheon+ SN data with Galactic coordinates."""
    df = pd.read_csv("data/pantheon_plus.dat", sep=r'\s+', comment='#')
    
    for ra_col in ['RA', 'HOST_RA', 'ra']:
        if ra_col in df.columns:
            break
    for dec_col in ['DECL', 'DEC', 'HOST_DEC', 'dec']:
        if dec_col in df.columns:
            break
    
    df['ra_sn'] = df[ra_col].astype(float)
    df['dec_sn'] = df[dec_col].astype(float)
    
    # Compute Hubble residual
    if 'MU_SH0ES' in df.columns:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=73.04, Om0=0.3)
        mu_theory = cosmo.distmod(df['zHD'].values).value
        df['mu_resid'] = df['MU_SH0ES'] - mu_theory
    
    # Compute Galactic coordinates
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    coords = SkyCoord(ra=df['ra_sn'].values*u.deg, dec=df['dec_sn'].values*u.deg, frame='icrs')
    gal = coords.galactic
    df['gl'] = gal.l.deg
    df['gb'] = gal.b.deg
    df['abs_gb'] = np.abs(df['gb'])
    
    # E(B-V) from Pantheon if available
    if 'MWEBV' in df.columns:
        df['ebv'] = df['MWEBV'].astype(float)
    
    mask = (df['zHD'] > 0.01) & df['x1'].notna() & df['c'].notna()
    mask &= df['ra_sn'].notna() & df['dec_sn'].notna()
    
    return df[mask].reset_index(drop=True)


def load_frb():
    """Load CHIME FRB Catalog 1."""
    df = pd.read_csv("data/chimefrbcat1.csv", sep='\t')
    df = df.replace(-9999, np.nan).replace('-9999', np.nan)
    
    non_rep = df[df['repeater_name'].isna()].drop_duplicates(subset='tns_name')
    rep = df[df['repeater_name'].notna()].drop_duplicates(subset='repeater_name')
    sightlines = pd.concat([non_rep, rep]).reset_index(drop=True)
    
    sightlines = sightlines[sightlines['dm_exc_ne2001'].notna()].copy()
    sightlines['ra_frb'] = sightlines['ra'].astype(float)
    sightlines['dec_frb'] = sightlines['dec'].astype(float)
    sightlines['dm_exc'] = sightlines['dm_exc_ne2001'].astype(float)
    
    # Compute Galactic latitude for FRBs too
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    coords = SkyCoord(ra=sightlines['ra_frb'].values*u.deg,
                      dec=sightlines['dec_frb'].values*u.deg, frame='icrs')
    gal = coords.galactic
    sightlines['gb_frb'] = gal.b.deg
    sightlines['abs_gb_frb'] = np.abs(sightlines['gb_frb'])
    
    return sightlines[['tns_name', 'ra_frb', 'dec_frb', 'dm_exc',
                        'dm_fitb', 'gb_frb', 'abs_gb_frb']].reset_index(drop=True)


def ra_dec_to_xyz(ra, dec):
    ra_r, dec_r = np.radians(ra), np.radians(dec)
    x = np.cos(dec_r) * np.cos(ra_r)
    y = np.cos(dec_r) * np.sin(ra_r)
    z = np.sin(dec_r)
    return np.column_stack([x, y, z])


def assign_dm(sn_df, frb_df, n_neighbors=5, max_sep_deg=30.0):
    """Assign DM with nearest-neighbor IDW, returning sep to nearest."""
    frb_xyz = ra_dec_to_xyz(frb_df['ra_frb'].values, frb_df['dec_frb'].values)
    sn_xyz = ra_dec_to_xyz(sn_df['ra_sn'].values, sn_df['dec_sn'].values)
    
    tree = cKDTree(frb_xyz)
    dists, idxs = tree.query(sn_xyz, k=n_neighbors)
    ang_seps = np.degrees(2 * np.arcsin(np.clip(dists / 2, 0, 1)))
    
    dm_values = frb_df['dm_exc'].values
    weights = 1.0 / (ang_seps + 0.1)
    mask = ang_seps <= max_sep_deg
    weights = weights * mask
    
    weighted_dm = np.sum(weights * dm_values[idxs], axis=1)
    total_weight = np.sum(weights, axis=1)
    
    dm_assigned = np.where(total_weight > 0, weighted_dm / total_weight, np.nan)
    nearest_sep = ang_seps[:, 0]
    
    return dm_assigned, nearest_sep


def partial_spearman(x, y, covariates):
    """Spearman partial correlation controlling for covariates."""
    from scipy.stats import spearmanr
    # Rank-transform everything
    x_rank = stats.rankdata(x)
    y_rank = stats.rankdata(y)
    cov_ranks = np.column_stack([stats.rankdata(c) for c in covariates.T])
    
    # Residualize x and y against covariates via OLS
    X_design = np.column_stack([cov_ranks, np.ones(len(x_rank))])
    
    beta_x = np.linalg.lstsq(X_design, x_rank, rcond=None)[0]
    resid_x = x_rank - X_design @ beta_x
    
    beta_y = np.linalg.lstsq(X_design, y_rank, rcond=None)[0]
    resid_y = y_rank - X_design @ beta_y
    
    rho, p = stats.pearsonr(resid_x, resid_y)
    return rho, p


def footprint_preserving_permutation(dm_vals, dec_vals, n_perm=10000, n_bands=10):
    """Shuffle DM within declination bands to preserve sky geometry."""
    dec_bins = pd.qcut(dec_vals, n_bands, labels=False, duplicates='drop')
    
    permuted = np.zeros((n_perm, len(dm_vals)))
    for i in range(n_perm):
        perm = dm_vals.copy()
        for band in np.unique(dec_bins):
            band_mask = dec_bins == band
            band_indices = np.where(band_mask)[0]
            perm[band_indices] = np.random.permutation(dm_vals[band_indices])
        permuted[i] = perm
    
    return permuted


def bootstrap_correlation(x, y, n_boot=5000):
    n = len(x)
    rhos = np.zeros(n_boot)
    for i in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        rho, _ = stats.spearmanr(x[idx], y[idx])
        rhos[i] = rho
    ci_lo, ci_hi = np.percentile(rhos, [2.5, 97.5])
    rho_obs, p_obs = stats.spearmanr(x, y)
    return rho_obs, p_obs, ci_lo, ci_hi


def run_tests():
    results = {}
    
    print("=" * 70)
    print("CLOSURE THEORY — FRB DM CROSS-MATCH v3 (with GPT fixes)")
    print("=" * 70)
    
    # ---- Load data ----
    print("\n[1] Loading data...")
    sn = load_pantheon()
    frb = load_frb()
    print(f"    Pantheon+ SNe: {len(sn)}")
    print(f"    CHIME FRB sightlines: {len(frb)}")
    
    # Filter FRBs to |b| > 20° (reduce MW contamination)
    frb_clean = frb[frb['abs_gb_frb'] > 20.0].reset_index(drop=True)
    print(f"    FRBs after |b|>20° cut: {len(frb_clean)}")
    
    results['n_sn'] = len(sn)
    results['n_frb_raw'] = len(frb)
    results['n_frb_clean'] = len(frb_clean)
    
    # ---- Assign DM (using cleaned FRBs) ----
    print("\n[2] Assigning DM_excess to SNe...")
    dm_assigned, nearest_sep = assign_dm(sn, frb_clean, n_neighbors=5, max_sep_deg=30.0)
    sn['dm_exc'] = dm_assigned
    sn['nearest_frb_sep'] = nearest_sep
    
    valid = sn['dm_exc'].notna()
    sn_valid = sn[valid].copy()
    print(f"    SNe with valid DM: {len(sn_valid)} / {len(sn)}")
    print(f"    Median nearest FRB: {sn_valid['nearest_frb_sep'].median():.1f}°")
    
    # ============================================================
    # TEST R1: RESOLUTION CURVE — THE DECISIVE TEST
    # ============================================================
    print("\n" + "=" * 70)
    print("[R1] RESOLUTION CURVE: rho vs max_sep cutoff")
    print("=" * 70)
    print("  If signal appears only when DM is local, the interpolation was killing it.")
    print("  If no signal at any cutoff, DM genuinely doesn't correlate.\n")
    
    cutoffs = [2, 3, 5, 7, 10, 15, 20, 30]
    resolution_curve = {}
    
    for param in ['x1', 'c', 'mu_resid']:
        if param not in sn_valid.columns:
            continue
        print(f"  {param} vs DM:")
        curve = []
        for cutoff in cutoffs:
            subset = sn_valid[sn_valid['nearest_frb_sep'] <= cutoff]
            n = len(subset)
            if n < 30:
                print(f"    max_sep ≤ {cutoff:2d}°: N={n:4d} (too few)")
                curve.append({'cutoff': cutoff, 'n': n, 'rho': None, 'p': None})
                continue
            rho, p = stats.spearmanr(subset[param], subset['dm_exc'])
            sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"    max_sep ≤ {cutoff:2d}°: N={n:4d}, rho={rho:+.4f}, p={p:.4f} {sig}")
            curve.append({'cutoff': cutoff, 'n': n, 'rho': float(rho), 'p': float(p)})
        resolution_curve[param] = curve
        print()
    
    results['R1_resolution_curve'] = resolution_curve
    
    # ============================================================
    # TEST R2: PARTIAL CORRELATIONS (controlling for confounds)
    # ============================================================
    print("=" * 70)
    print("[R2] PARTIAL CORRELATIONS (controlling |b|, E(B-V), dec, z)")
    print("=" * 70)
    
    # Build covariate matrix
    cov_cols = ['abs_gb', 'dec_sn', 'zHD']
    if 'ebv' in sn_valid.columns:
        cov_cols.append('ebv')
    
    covariates = sn_valid[cov_cols].values
    
    for param in ['x1', 'c', 'mu_resid']:
        if param not in sn_valid.columns:
            continue
        rho_raw, p_raw = stats.spearmanr(sn_valid[param], sn_valid['dm_exc'])
        rho_partial, p_partial = partial_spearman(
            sn_valid[param].values, sn_valid['dm_exc'].values, covariates
        )
        print(f"  {param:8s} vs DM:  raw rho={rho_raw:+.4f} (p={p_raw:.4f})"
              f"  |  partial rho={rho_partial:+.4f} (p={p_partial:.4f})"
              f"  |  controls: {cov_cols}")
        results[f'R2_{param}'] = {
            'rho_raw': float(rho_raw), 'p_raw': float(p_raw),
            'rho_partial': float(rho_partial), 'p_partial': float(p_partial),
            'controls': cov_cols
        }
    
    # ============================================================
    # TEST R3: FOOTPRINT-PRESERVING PERMUTATION
    # ============================================================
    print("\n" + "=" * 70)
    print("[R3] FOOTPRINT-PRESERVING PERMUTATION (shuffle within dec bands)")
    print("=" * 70)
    
    n_perm = 10000
    for param in ['x1', 'c']:
        vals = sn_valid[param].values
        dm_vals = sn_valid['dm_exc'].values
        dec_vals = sn_valid['dec_sn'].values
        
        rho_obs, _ = stats.spearmanr(vals, dm_vals)
        
        # Footprint-preserving permutation
        permuted_dm = footprint_preserving_permutation(dm_vals, dec_vals, n_perm=n_perm)
        null_rhos = np.array([stats.spearmanr(vals, permuted_dm[i])[0] for i in range(n_perm)])
        
        p_fp = np.mean(np.abs(null_rhos) >= np.abs(rho_obs))
        
        # Compare with naive permutation
        naive_rhos = np.array([stats.spearmanr(vals, np.random.permutation(dm_vals))[0]
                              for _ in range(n_perm)])
        p_naive = np.mean(np.abs(naive_rhos) >= np.abs(rho_obs))
        
        print(f"  {param} vs DM:  rho_obs={rho_obs:+.4f}")
        print(f"    Naive permutation p     = {p_naive:.4f}")
        print(f"    Footprint-preserving p  = {p_fp:.4f}")
        if p_naive < 0.05 and p_fp > 0.05:
            print(f"    → Signal DISAPPEARS with footprint control — spatial artifact! ⚠️")
        elif p_fp < 0.05:
            print(f"    → Signal SURVIVES footprint control — may be real")
        else:
            print(f"    → No signal in either")
        
        results[f'R3_{param}'] = {
            'rho_obs': float(rho_obs),
            'p_naive': float(p_naive),
            'p_footprint': float(p_fp),
        }
    
    # ============================================================
    # TEST R4: |b| > 20° RESTRICTED ANALYSIS
    # ============================================================
    print("\n" + "=" * 70)
    print("[R4] HIGH GALACTIC LATITUDE ONLY (|b| > 20°)")
    print("=" * 70)
    
    high_b = sn_valid[sn_valid['abs_gb'] > 20.0]
    print(f"  SNe with |b|>20°: {len(high_b)} / {len(sn_valid)}")
    
    for param in ['x1', 'c', 'mu_resid']:
        if param not in high_b.columns:
            continue
        rho, p, ci_lo, ci_hi = bootstrap_correlation(
            high_b[param].values, high_b['dm_exc'].values, n_boot=5000
        )
        print(f"  {param:8s} vs DM: rho={rho:+.4f}, p={p:.4f}, CI=[{ci_lo:+.4f}, {ci_hi:+.4f}]")
        results[f'R4_{param}'] = {
            'rho': float(rho), 'p': float(p),
            'ci_lo': float(ci_lo), 'ci_hi': float(ci_hi),
            'n': len(high_b)
        }
    
    # ============================================================
    # TEST R5: REDSHIFT SPLIT (z < 0.82 vs z >= 0.82)
    # ============================================================
    print("\n" + "=" * 70)
    print("[R5] REDSHIFT SPLIT with controls")
    print("=" * 70)
    
    for z_label, z_lo, z_hi in [("z<0.82", 0.01, 0.82), ("z>=0.82", 0.82, 3.0)]:
        subset = sn_valid[(sn_valid['zHD'] >= z_lo) & (sn_valid['zHD'] < z_hi)]
        if len(subset) < 20:
            print(f"  {z_label}: N={len(subset)} (too few)")
            continue
        
        for param in ['x1', 'c']:
            rho, p = stats.spearmanr(subset[param], subset['dm_exc'])
            # Partial correlation
            cov = subset[['abs_gb', 'dec_sn']].values
            if 'ebv' in subset.columns:
                cov = np.column_stack([cov, subset['ebv'].values])
            rho_p, p_p = partial_spearman(subset[param].values, subset['dm_exc'].values, cov)
            
            print(f"  {z_label} {param}: N={len(subset)}, raw rho={rho:+.4f} (p={p:.4f})"
                  f", partial rho={rho_p:+.4f} (p={p_p:.4f})")
            results[f'R5_{z_label}_{param}'] = {
                'n': len(subset),
                'rho_raw': float(rho), 'p_raw': float(p),
                'rho_partial': float(rho_p), 'p_partial': float(p_p),
            }
    
    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    # Check if resolution curve shows emerging signal
    print("\n  Resolution curve verdict:")
    for param in ['x1', 'c']:
        if param in resolution_curve:
            significant = [c for c in resolution_curve[param]
                          if c['p'] is not None and c['p'] < 0.05]
            if significant:
                best = min(significant, key=lambda c: c['p'])
                print(f"    {param}: Signal at cutoff {best['cutoff']}° "
                      f"(rho={best['rho']:+.4f}, p={best['p']:.4f}, N={best['n']})")
            else:
                print(f"    {param}: No signal at any cutoff")
    
    print(f"\n  Footprint-preserving permutation:")
    for param in ['x1', 'c']:
        key = f'R3_{param}'
        if key in results:
            r = results[key]
            verdict = "ARTIFACT" if r['p_naive'] < 0.05 and r['p_footprint'] > 0.05 else \
                      "SURVIVES" if r['p_footprint'] < 0.05 else "NULL"
            print(f"    {param}: naive p={r['p_naive']:.4f}, footprint p={r['p_footprint']:.4f} → {verdict}")
    
    # Save
    with open(RESULTS_DIR / "frb_v3_results.json", 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n  Results saved to {RESULTS_DIR}/frb_v3_results.json")
    return results


if __name__ == "__main__":
    run_tests()
