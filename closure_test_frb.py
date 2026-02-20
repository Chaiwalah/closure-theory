#!/usr/bin/env python3
"""
Closure Theory — FRB DM Cross-Match Test
=========================================
Tests whether SN Ia light-curve parameters correlate with FRB dispersion
measure (DM_excess) along similar sightlines.

Key prediction:
- If closure is frequency-dependent: x1 (stretch, temporal-only) should NOT
  correlate with DM, while c (color, frequency-dependent) MAY correlate.
- If x1 correlates with DM: challenges the frequency-only fingerprint.
- If null for both: confirms the effect is not baryonic-matter-mediated.

Method:
1. Load Pantheon+ SN Ia catalog (z, c, x1, ra, dec, mu_residual)
2. Load CHIME FRB Catalog 1 (ra, dec, DM_excess)
3. Build GP-interpolated DM_excess sky map from FRB positions
4. Sample DM_excess at each SN position (nearest-neighbor + GP)
5. Test correlations: DM vs c, DM vs x1, DM vs mu_residual
6. Split by redshift bins (below/above z=0.82 threshold)

Data: CHIME FRB Catalog 1 (536 FRBs, Amiri et al. 2021)
       Pantheon+SH0ES (1701 SNe Ia)
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

RESULTS_DIR = Path("results_frb")
RESULTS_DIR.mkdir(exist_ok=True)


def load_pantheon():
    """Load Pantheon+ SN data."""
    df = pd.read_csv("data/pantheon_plus.dat", sep=r'\s+', comment='#')
    # Need: zHD, x1, c, MU_SH0ES, MU_SHOES_ERR_DIAG, RA, DEC (or DECL)
    # Also need HOST_RA, HOST_DEC if available
    cols_needed = ['zHD', 'x1', 'c']
    
    # Find RA/DEC columns
    for ra_col in ['RA', 'HOST_RA', 'ra']:
        if ra_col in df.columns:
            break
    for dec_col in ['DECL', 'DEC', 'HOST_DEC', 'dec']:
        if dec_col in df.columns:
            break
    
    # Compute Hubble residual
    if 'MU_SH0ES' in df.columns:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=73.04, Om0=0.3)
        mu_theory = cosmo.distmod(df['zHD'].values).value
        df['mu_resid'] = df['MU_SH0ES'] - mu_theory
    
    df['ra_sn'] = df[ra_col].astype(float)
    df['dec_sn'] = df[dec_col].astype(float)
    
    mask = (df['zHD'] > 0.01) & df['x1'].notna() & df['c'].notna()
    mask &= df['ra_sn'].notna() & df['dec_sn'].notna()
    
    return df[mask].reset_index(drop=True)


def load_frb():
    """Load CHIME FRB Catalog 1."""
    df = pd.read_csv("data/chimefrbcat1.csv", sep='\t')
    # Replace -9999 with NaN (handle both numeric and string)
    df = df.replace(-9999, np.nan).replace('-9999', np.nan)
    
    # Use one entry per unique sightline
    # Non-repeaters: repeater_name is NaN (was -9999)
    # Repeaters: take first burst per repeater source
    non_rep = df[df['repeater_name'].isna()].drop_duplicates(subset='tns_name')
    rep = df[df['repeater_name'].notna()].drop_duplicates(subset='repeater_name')
    sightlines = pd.concat([non_rep, rep]).reset_index(drop=True)
    
    # Use NE2001 excess DM (Milky Way subtracted)
    sightlines = sightlines[sightlines['dm_exc_ne2001'].notna()].copy()
    sightlines['ra_frb'] = sightlines['ra'].astype(float)
    sightlines['dec_frb'] = sightlines['dec'].astype(float)
    sightlines['dm_exc'] = sightlines['dm_exc_ne2001'].astype(float)
    
    return sightlines[['tns_name', 'ra_frb', 'dec_frb', 'dm_exc', 'dm_fitb']].reset_index(drop=True)


def angular_separation(ra1, dec1, ra2, dec2):
    """Compute angular separation in degrees between two sets of coordinates."""
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    dra = ra2 - ra1
    ddec = dec2 - dec1
    a = np.sin(ddec/2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dra/2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))))


def ra_dec_to_xyz(ra, dec):
    """Convert RA, Dec (degrees) to unit sphere XYZ for KDTree."""
    ra_r, dec_r = np.radians(ra), np.radians(dec)
    x = np.cos(dec_r) * np.cos(ra_r)
    y = np.cos(dec_r) * np.sin(ra_r)
    z = np.sin(dec_r)
    return np.column_stack([x, y, z])


def assign_dm_to_sn(sn_df, frb_df, n_neighbors=5, max_sep_deg=30.0):
    """
    Assign DM_excess to each SN using inverse-distance-weighted average
    of nearest FRB sightlines.
    """
    frb_xyz = ra_dec_to_xyz(frb_df['ra_frb'].values, frb_df['dec_frb'].values)
    sn_xyz = ra_dec_to_xyz(sn_df['ra_sn'].values, sn_df['dec_sn'].values)
    
    tree = cKDTree(frb_xyz)
    
    # Query k nearest neighbors
    dists, idxs = tree.query(sn_xyz, k=n_neighbors)
    
    # Convert chord distance to angular separation (approx)
    # chord = 2*sin(theta/2), so theta = 2*arcsin(chord/2)
    ang_seps = np.degrees(2 * np.arcsin(np.clip(dists / 2, 0, 1)))
    
    dm_values = frb_df['dm_exc'].values
    
    # Inverse-distance weighted average (with floor to avoid division by zero)
    weights = 1.0 / (ang_seps + 0.1)  # 0.1 deg floor
    
    # Mask neighbors beyond max separation
    mask = ang_seps <= max_sep_deg
    weights = weights * mask
    
    weighted_dm = np.sum(weights * dm_values[idxs], axis=1)
    total_weight = np.sum(weights, axis=1)
    
    dm_assigned = np.where(total_weight > 0, weighted_dm / total_weight, np.nan)
    nearest_sep = ang_seps[:, 0]
    
    return dm_assigned, nearest_sep


def bootstrap_correlation(x, y, n_boot=10000):
    """Bootstrap Spearman correlation with confidence interval."""
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
    print("CLOSURE THEORY — FRB DM CROSS-MATCH TEST")
    print("=" * 70)
    
    # ---- Load data ----
    print("\n[1] Loading data...")
    sn = load_pantheon()
    frb = load_frb()
    print(f"    Pantheon+ SNe: {len(sn)}")
    print(f"    CHIME FRB sightlines: {len(frb)}")
    print(f"    FRB DM_exc range: {frb['dm_exc'].min():.0f} - {frb['dm_exc'].max():.0f} pc/cm³")
    print(f"    FRB sky coverage: RA [{frb['ra_frb'].min():.0f}, {frb['ra_frb'].max():.0f}], "
          f"Dec [{frb['dec_frb'].min():.0f}, {frb['dec_frb'].max():.0f}]")
    
    results['n_sn'] = len(sn)
    results['n_frb'] = len(frb)
    
    # ---- Assign DM to SNe ----
    print("\n[2] Assigning DM_excess to SNe via nearest-neighbor interpolation...")
    dm_assigned, nearest_sep = assign_dm_to_sn(sn, frb, n_neighbors=5, max_sep_deg=30.0)
    
    sn['dm_exc'] = dm_assigned
    sn['nearest_frb_sep'] = nearest_sep
    
    valid = sn['dm_exc'].notna()
    print(f"    SNe with valid DM assignment: {valid.sum()} / {len(sn)}")
    print(f"    Median nearest FRB separation: {sn['nearest_frb_sep'].median():.1f}°")
    print(f"    Mean assigned DM_exc: {sn.loc[valid, 'dm_exc'].mean():.0f} pc/cm³")
    
    results['n_sn_with_dm'] = int(valid.sum())
    results['median_frb_sep_deg'] = float(sn['nearest_frb_sep'].median())
    
    sn_valid = sn[valid].copy()
    
    # ---- Test F1: x1 vs DM (THE KEY TEST) ----
    print("\n" + "=" * 70)
    print("[TEST F1] x1 (stretch) vs DM_excess — Frequency Fingerprint Test")
    print("=" * 70)
    print("  Prediction: NULL (x1 is temporal-only, should not correlate with baryons)")
    
    rho, p, ci_lo, ci_hi = bootstrap_correlation(
        sn_valid['x1'].values, sn_valid['dm_exc'].values
    )
    print(f"  Spearman rho = {rho:.4f}  (p = {p:.4e})")
    print(f"  95% CI: [{ci_lo:.4f}, {ci_hi:.4f}]")
    print(f"  Verdict: {'NULL ✓ (as predicted)' if p > 0.05 else 'SIGNAL ✗ (challenges theory)'}")
    
    results['F1_x1_dm'] = {
        'rho': float(rho), 'p': float(p),
        'ci_lo': float(ci_lo), 'ci_hi': float(ci_hi),
        'null': bool(p > 0.05)
    }
    
    # ---- Test F2: c (color) vs DM ----
    print("\n" + "=" * 70)
    print("[TEST F2] c (color) vs DM_excess — Color-Baryon Coupling")
    print("=" * 70)
    print("  Prediction: Possible weak signal (color is frequency-dependent)")
    
    rho, p, ci_lo, ci_hi = bootstrap_correlation(
        sn_valid['c'].values, sn_valid['dm_exc'].values
    )
    print(f"  Spearman rho = {rho:.4f}  (p = {p:.4e})")
    print(f"  95% CI: [{ci_lo:.4f}, {ci_hi:.4f}]")
    
    results['F2_c_dm'] = {
        'rho': float(rho), 'p': float(p),
        'ci_lo': float(ci_lo), 'ci_hi': float(ci_hi),
    }
    
    # ---- Test F3: mu_resid vs DM ----
    if 'mu_resid' in sn_valid.columns:
        print("\n" + "=" * 70)
        print("[TEST F3] μ_resid vs DM_excess — Hubble Residual vs Baryons")
        print("=" * 70)
        print("  Prediction: Null if closure is not matter-mediated")
        
        rho, p, ci_lo, ci_hi = bootstrap_correlation(
            sn_valid['mu_resid'].values, sn_valid['dm_exc'].values
        )
        print(f"  Spearman rho = {rho:.4f}  (p = {p:.4e})")
        print(f"  95% CI: [{ci_lo:.4f}, {ci_hi:.4f}]")
        
        results['F3_mu_dm'] = {
            'rho': float(rho), 'p': float(p),
            'ci_lo': float(ci_lo), 'ci_hi': float(ci_hi),
        }
    
    # ---- Test F4: Redshift-split analysis ----
    print("\n" + "=" * 70)
    print("[TEST F4] Redshift-split: x1 vs DM below/above z=0.82")
    print("=" * 70)
    
    z_thresh = 0.82
    lo = sn_valid[sn_valid['zHD'] < z_thresh]
    hi = sn_valid[sn_valid['zHD'] >= z_thresh]
    
    for label, subset in [("z < 0.82", lo), ("z >= 0.82", hi)]:
        if len(subset) < 20:
            print(f"  {label}: Too few SNe ({len(subset)}), skipping")
            continue
        rho, p, ci_lo, ci_hi = bootstrap_correlation(
            subset['x1'].values, subset['dm_exc'].values
        )
        print(f"  {label} (N={len(subset)}): rho={rho:.4f}, p={p:.4e}, CI=[{ci_lo:.4f}, {ci_hi:.4f}]")
        results[f'F4_x1_dm_{label.replace(" ", "").replace(".", "").replace("=","")}'] = {
            'rho': float(rho), 'p': float(p), 'n': len(subset),
            'ci_lo': float(ci_lo), 'ci_hi': float(ci_hi),
        }
    
    # ---- Test F5: c vs DM redshift-split ----
    print("\n" + "=" * 70)
    print("[TEST F5] Redshift-split: c vs DM below/above z=0.82")
    print("=" * 70)
    
    for label, subset in [("z < 0.82", lo), ("z >= 0.82", hi)]:
        if len(subset) < 20:
            print(f"  {label}: Too few SNe ({len(subset)}), skipping")
            continue
        rho, p, ci_lo, ci_hi = bootstrap_correlation(
            subset['c'].values, subset['dm_exc'].values
        )
        print(f"  {label} (N={len(subset)}): rho={rho:.4f}, p={p:.4e}, CI=[{ci_lo:.4f}, {ci_hi:.4f}]")
        results[f'F5_c_dm_{label.replace(" ", "").replace(".", "").replace("=","")}'] = {
            'rho': float(rho), 'p': float(p), 'n': len(subset),
            'ci_lo': float(ci_lo), 'ci_hi': float(ci_hi),
        }
    
    # ---- Test F6: DM gradient vs angular proximity ----
    print("\n" + "=" * 70)
    print("[TEST F6] Control: Does DM assignment quality affect results?")
    print("=" * 70)
    print("  Testing whether nearest-FRB separation correlates with SN observables")
    
    for param in ['x1', 'c']:
        rho, p = stats.spearmanr(sn_valid['nearest_frb_sep'], sn_valid[param])
        print(f"  nearest_sep vs {param}: rho={rho:.4f}, p={p:.4e}")
        results[f'F6_sep_vs_{param}'] = {'rho': float(rho), 'p': float(p)}
    
    # ---- Test F7: Permutation test (null distribution) ----
    print("\n" + "=" * 70)
    print("[TEST F7] Permutation null: shuffled DM vs x1")
    print("=" * 70)
    
    n_perm = 10000
    x1_vals = sn_valid['x1'].values
    dm_vals = sn_valid['dm_exc'].values
    rho_obs, _ = stats.spearmanr(x1_vals, dm_vals)
    
    null_rhos = np.zeros(n_perm)
    for i in range(n_perm):
        shuffled = np.random.permutation(dm_vals)
        null_rhos[i], _ = stats.spearmanr(x1_vals, shuffled)
    
    p_perm = np.mean(np.abs(null_rhos) >= np.abs(rho_obs))
    print(f"  Observed rho: {rho_obs:.4f}")
    print(f"  Permutation p-value: {p_perm:.4f}")
    print(f"  Null distribution: mean={null_rhos.mean():.4f}, std={null_rhos.std():.4f}")
    
    results['F7_permutation'] = {
        'rho_obs': float(rho_obs),
        'p_perm': float(p_perm),
        'null_mean': float(null_rhos.mean()),
        'null_std': float(null_rhos.std()),
    }
    
    # ---- Summary ----
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    f1_null = results['F1_x1_dm']['null']
    print(f"  F1 (x1 vs DM):     {'NULL ✓' if f1_null else 'SIGNAL ✗'} — "
          f"rho={results['F1_x1_dm']['rho']:.4f}, p={results['F1_x1_dm']['p']:.4e}")
    print(f"  F2 (c vs DM):      rho={results['F2_c_dm']['rho']:.4f}, p={results['F2_c_dm']['p']:.4e}")
    if 'F3_mu_dm' in results:
        print(f"  F3 (μ vs DM):      rho={results['F3_mu_dm']['rho']:.4f}, p={results['F3_mu_dm']['p']:.4e}")
    print(f"  F7 (perm null):    p_perm={results['F7_permutation']['p_perm']:.4f}")
    
    if f1_null:
        print("\n  ✅ FREQUENCY FINGERPRINT CONFIRMED:")
        print("     x1 (temporal-only) does NOT correlate with baryonic matter tracer.")
        print("     Consistent with closure affecting frequency-dependent observables only.")
    else:
        print("\n  ⚠️  FREQUENCY FINGERPRINT CHALLENGED:")
        print("     x1 correlates with DM — needs investigation.")
    
    # Save results
    with open(RESULTS_DIR / "frb_crossmatch_results.json", 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n  Results saved to {RESULTS_DIR}/frb_crossmatch_results.json")
    
    return results


if __name__ == "__main__":
    run_tests()
