#!/usr/bin/env python3
"""
Closure Theory — Anisotropy Test
==================================
Tests whether the closure signal (information entanglement/compression)
is isotropic (same in all directions) or has directional structure.

If the boundary is a perfect sphere: signal should be uniform across the sky.
If the boundary has geometry: signal should show preferred directions.

Method:
1. Compute closure metrics per SN (color residual, entanglement proxy)
2. Split sky into patches and compare signal strength
3. Compute low-ell spherical harmonic modes of signal map
4. Test for dipole/quadrupole structure
5. Compare with known anisotropy axes (CMB dipole, bulk flow)

This test challenges the assumption that the observable universe
horizon is spherically symmetric.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

RESULTS_DIR = Path("results_anisotropy")
RESULTS_DIR.mkdir(exist_ok=True)


def load_data():
    """Load Pantheon+ with sky coordinates and compute closure metrics."""
    df = pd.read_csv("data/pantheon_plus.dat", sep=r'\s+', comment='#')
    
    for ra_col in ['RA', 'HOST_RA', 'ra']:
        if ra_col in df.columns:
            break
    for dec_col in ['DECL', 'DEC', 'HOST_DEC', 'dec']:
        if dec_col in df.columns:
            break
    
    df['ra'] = df[ra_col].astype(float)
    df['dec'] = df[dec_col].astype(float)
    
    # Compute Hubble residual
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=73.04, Om0=0.3)
    mu_theory = cosmo.distmod(df['zHD'].values).value
    df['mu_resid'] = df['MU_SH0ES'] - mu_theory
    
    # Compute Galactic and ecliptic coordinates
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    coords = SkyCoord(ra=df['ra'].values*u.deg, dec=df['dec'].values*u.deg, frame='icrs')
    gal = coords.galactic
    ecl = coords.barycentricmeanecliptic
    df['gl'] = gal.l.deg
    df['gb'] = gal.b.deg
    df['ecl_lon'] = ecl.lon.deg
    df['ecl_lat'] = ecl.lat.deg
    
    # Closure metrics: residuals after standard corrections
    # The "closure signal" is what's LEFT after SALT2 does its best
    # Key metric: |c_residual| — how anomalous is the color?
    # At high-z, closure predicts these residuals carry non-random structure
    
    # Color anomaly: deviation from Tripp relation expectation
    # In SALT2, c should scatter symmetrically around 0
    # Closure predicts systematic shifts accumulating with z
    df['c_abs'] = np.abs(df['c'])
    df['mu_abs_resid'] = np.abs(df['mu_resid'])
    
    # Compute "entanglement proxy": correlation between c and mu_resid
    # in local redshift windows (rolling measure of coupling)
    
    mask = (df['zHD'] > 0.01) & df['x1'].notna() & df['c'].notna()
    mask &= df['ra'].notna() & df['dec'].notna()
    mask &= df['mu_resid'].notna()
    
    return df[mask].reset_index(drop=True)


def compute_directional_signal(df, z_min=0.3, z_max=1.5):
    """
    Compute the closure signal strength per SN for directional analysis.
    
    The "signal" is: how much does the color-distance coupling deviate
    from what's expected at that redshift?
    
    For each SN, we compute:
    - c * mu_resid product (entanglement proxy: should be ~0 if independent)
    - |c| anomaly relative to redshift trend
    """
    subset = df[(df['zHD'] >= z_min) & (df['zHD'] <= z_max)].copy()
    
    # Detrend c against z (remove known trend)
    z = subset['zHD'].values
    c = subset['c'].values
    mu_r = subset['mu_resid'].values
    
    # Fit and remove z-trend from c
    coeffs_c = np.polyfit(z, c, 2)
    c_detrended = c - np.polyval(coeffs_c, z)
    
    # Fit and remove z-trend from mu_resid
    coeffs_mu = np.polyfit(z, mu_r, 2)
    mu_detrended = mu_r - np.polyval(coeffs_mu, z)
    
    # Entanglement proxy: product of detrended residuals
    # If c and mu are entangled, this product is systematically positive or negative
    # If independent, it scatters around 0
    subset['entanglement_proxy'] = c_detrended * mu_detrended
    subset['c_detrended'] = c_detrended
    subset['mu_detrended'] = mu_detrended
    
    return subset


def sky_patch_analysis(df, n_patches_ra=6, n_patches_dec=3):
    """
    Divide sky into patches and compute signal in each.
    """
    ra_bins = np.linspace(0, 360, n_patches_ra + 1)
    dec_bins = np.linspace(-90, 90, n_patches_dec + 1)
    
    patches = []
    for i in range(n_patches_ra):
        for j in range(n_patches_dec):
            mask = ((df['ra'] >= ra_bins[i]) & (df['ra'] < ra_bins[i+1]) &
                    (df['dec'] >= dec_bins[j]) & (df['dec'] < dec_bins[j+1]))
            patch = df[mask]
            if len(patch) >= 10:
                mean_ent = patch['entanglement_proxy'].mean()
                std_ent = patch['entanglement_proxy'].std() / np.sqrt(len(patch))
                mean_c = patch['c_detrended'].mean()
                
                # Local c-mu correlation
                if len(patch) >= 20:
                    rho, p = stats.spearmanr(patch['c_detrended'], patch['mu_detrended'])
                else:
                    rho, p = np.nan, np.nan
                
                patches.append({
                    'ra_center': (ra_bins[i] + ra_bins[i+1]) / 2,
                    'dec_center': (dec_bins[j] + dec_bins[j+1]) / 2,
                    'n': len(patch),
                    'mean_entanglement': float(mean_ent),
                    'sem_entanglement': float(std_ent),
                    'mean_c_detrended': float(mean_c),
                    'local_rho': float(rho) if not np.isnan(rho) else None,
                    'local_p': float(p) if not np.isnan(p) else None,
                    'mean_z': float(patch['zHD'].mean()),
                })
    
    return patches


def dipole_fit(ra, dec, signal):
    """
    Fit a dipole (ℓ=1) to the signal on the sky.
    Returns dipole direction and amplitude.
    
    Signal = A * cos(theta) + offset
    where theta is angle from dipole axis.
    
    We scan over directions to find the best dipole axis.
    """
    ra_r = np.radians(ra)
    dec_r = np.radians(dec)
    
    # Convert to unit vectors
    x = np.cos(dec_r) * np.cos(ra_r)
    y = np.cos(dec_r) * np.sin(ra_r)
    z = np.sin(dec_r)
    
    # Fit: signal = a*x + b*y + c*z + d (linear dipole model)
    A = np.column_stack([x, y, z, np.ones_like(x)])
    result = np.linalg.lstsq(A, signal, rcond=None)
    coeffs = result[0]
    
    a, b, c, d = coeffs
    
    # Dipole amplitude and direction
    dipole_amp = np.sqrt(a**2 + b**2 + c**2)
    
    if dipole_amp > 0:
        # Direction of dipole
        dipole_dec = np.degrees(np.arcsin(c / dipole_amp))
        dipole_ra = np.degrees(np.arctan2(b, a)) % 360
    else:
        dipole_ra, dipole_dec = 0, 0
    
    # Compute explained variance
    signal_pred = A @ coeffs
    ss_res = np.sum((signal - signal_pred)**2)
    ss_tot = np.sum((signal - np.mean(signal))**2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    
    # Significance: compare to random dipole fits
    n_perm = 10000
    null_amps = np.zeros(n_perm)
    for i in range(n_perm):
        shuffled = np.random.permutation(signal)
        null_coeffs = np.linalg.lstsq(A, shuffled, rcond=None)[0]
        null_amps[i] = np.sqrt(null_coeffs[0]**2 + null_coeffs[1]**2 + null_coeffs[2]**2)
    
    p_dipole = np.mean(null_amps >= dipole_amp)
    
    return {
        'amplitude': float(dipole_amp),
        'ra': float(dipole_ra),
        'dec': float(dipole_dec),
        'r_squared': float(r_squared),
        'offset': float(d),
        'p_value': float(p_dipole),
    }


def quadrupole_fit(ra, dec, signal):
    """
    Fit ℓ=1 + ℓ=2 (dipole + quadrupole) to the signal.
    Tests for more complex directional structure.
    """
    ra_r = np.radians(ra)
    dec_r = np.radians(dec)
    
    x = np.cos(dec_r) * np.cos(ra_r)
    y = np.cos(dec_r) * np.sin(ra_r)
    z = np.sin(dec_r)
    
    # ℓ=1: x, y, z
    # ℓ=2: xy, xz, yz, x²-y², 3z²-1
    A = np.column_stack([
        x, y, z,                              # dipole
        x*y, x*z, y*z, x**2 - y**2, 3*z**2 - 1,  # quadrupole
        np.ones_like(x)                        # offset
    ])
    
    result = np.linalg.lstsq(A, signal, rcond=None)
    coeffs = result[0]
    
    signal_pred = A @ coeffs
    ss_res = np.sum((signal - signal_pred)**2)
    ss_tot = np.sum((signal - np.mean(signal))**2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    
    # F-test: does quadrupole add significant explanatory power over dipole?
    A_dipole = np.column_stack([x, y, z, np.ones_like(x)])
    coeffs_d = np.linalg.lstsq(A_dipole, signal, rcond=None)[0]
    ss_res_dipole = np.sum((signal - A_dipole @ coeffs_d)**2)
    
    n = len(signal)
    df_extra = 5  # quadrupole adds 5 params
    df_resid = n - 9  # total params in full model
    
    if ss_res > 0 and df_resid > 0:
        f_stat = ((ss_res_dipole - ss_res) / df_extra) / (ss_res / df_resid)
        from scipy.stats import f as f_dist
        p_quad = 1 - f_dist.cdf(f_stat, df_extra, df_resid)
    else:
        f_stat, p_quad = 0, 1
    
    # Permutation test for full model
    n_perm = 10000
    null_r2 = np.zeros(n_perm)
    for i in range(n_perm):
        shuffled = np.random.permutation(signal)
        null_c = np.linalg.lstsq(A, shuffled, rcond=None)[0]
        null_pred = A @ null_c
        null_ss_res = np.sum((shuffled - null_pred)**2)
        null_ss_tot = np.sum((shuffled - np.mean(shuffled))**2)
        null_r2[i] = 1 - null_ss_res / null_ss_tot if null_ss_tot > 0 else 0
    
    p_full = np.mean(null_r2 >= r_squared)
    
    return {
        'r_squared_full': float(r_squared),
        'p_full_model': float(p_full),
        'f_stat_quad': float(f_stat),
        'p_quadrupole': float(p_quad),
    }


def hemisphere_comparison(df, axes):
    """
    For each proposed axis, compare closure signal in "toward" vs "away" hemispheres.
    """
    ra_r = np.radians(df['ra'].values)
    dec_r = np.radians(df['dec'].values)
    
    x = np.cos(dec_r) * np.cos(ra_r)
    y = np.cos(dec_r) * np.sin(ra_r)
    z = np.sin(dec_r)
    
    results = {}
    for name, (ax_ra, ax_dec) in axes.items():
        ax_ra_r, ax_dec_r = np.radians(ax_ra), np.radians(ax_dec)
        ax_x = np.cos(ax_dec_r) * np.cos(ax_ra_r)
        ax_y = np.cos(ax_dec_r) * np.sin(ax_ra_r)
        ax_z = np.sin(ax_dec_r)
        
        cos_angle = x * ax_x + y * ax_y + z * ax_z
        
        toward = df[cos_angle > 0]
        away = df[cos_angle <= 0]
        
        if len(toward) < 20 or len(away) < 20:
            continue
        
        # Compare entanglement proxy
        t_stat, p_val = stats.mannwhitneyu(
            toward['entanglement_proxy'].values,
            away['entanglement_proxy'].values,
            alternative='two-sided'
        )
        
        # Compare local c-mu correlation
        rho_toward, _ = stats.spearmanr(toward['c_detrended'], toward['mu_detrended'])
        rho_away, _ = stats.spearmanr(away['c_detrended'], away['mu_detrended'])
        
        results[name] = {
            'n_toward': len(toward),
            'n_away': len(away),
            'mean_ent_toward': float(toward['entanglement_proxy'].mean()),
            'mean_ent_away': float(away['entanglement_proxy'].mean()),
            'mannwhitney_p': float(p_val),
            'rho_toward': float(rho_toward),
            'rho_away': float(rho_away),
            'rho_diff': float(rho_toward - rho_away),
        }
    
    return results


def redshift_dependent_anisotropy(df, z_bins=None):
    """
    Does the dipole direction/amplitude change with redshift?
    If the boundary isn't spherical, the anisotropy should evolve with z.
    """
    if z_bins is None:
        z_bins = [(0.01, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 1.0), (1.0, 2.5)]
    
    results = []
    for z_lo, z_hi in z_bins:
        subset = df[(df['zHD'] >= z_lo) & (df['zHD'] < z_hi)]
        if len(subset) < 50:
            results.append({
                'z_range': f"{z_lo:.1f}-{z_hi:.1f}",
                'n': len(subset),
                'note': 'too few'
            })
            continue
        
        dipole = dipole_fit(subset['ra'].values, subset['dec'].values,
                           subset['entanglement_proxy'].values)
        
        results.append({
            'z_range': f"{z_lo:.1f}-{z_hi:.1f}",
            'n': len(subset),
            'dipole_amp': dipole['amplitude'],
            'dipole_ra': dipole['ra'],
            'dipole_dec': dipole['dec'],
            'dipole_p': dipole['p_value'],
        })
    
    return results


def run_tests():
    results = {}
    
    print("=" * 70)
    print("CLOSURE THEORY — ANISOTROPY TEST")
    print("Is the closure boundary spherically symmetric, or does it have")
    print("directional structure?")
    print("=" * 70)
    
    # ---- Load and prepare ----
    print("\n[1] Loading data and computing closure metrics...")
    df = load_data()
    print(f"    Total SNe: {len(df)}")
    
    # Compute directional signal for z > 0.3 (where closure should be active)
    df_signal = compute_directional_signal(df, z_min=0.1, z_max=2.5)
    print(f"    SNe with closure metrics (z=0.1-2.5): {len(df_signal)}")
    
    # Also high-z only
    df_highz = compute_directional_signal(df, z_min=0.5, z_max=2.5)
    print(f"    High-z subset (z=0.5-2.5): {len(df_highz)}")
    
    results['n_total'] = len(df_signal)
    results['n_highz'] = len(df_highz)
    
    # ============================================================
    # TEST A1: DIPOLE FIT
    # ============================================================
    print("\n" + "=" * 70)
    print("[A1] DIPOLE FIT — Is there a preferred direction?")
    print("=" * 70)
    
    for label, subset in [("All z>0.1", df_signal), ("High-z (z>0.5)", df_highz)]:
        print(f"\n  {label} (N={len(subset)}):")
        
        # Dipole on entanglement proxy
        dipole_ent = dipole_fit(subset['ra'].values, subset['dec'].values,
                               subset['entanglement_proxy'].values)
        print(f"    Entanglement proxy dipole:")
        print(f"      Amplitude: {dipole_ent['amplitude']:.6f}")
        print(f"      Direction: RA={dipole_ent['ra']:.1f}°, Dec={dipole_ent['dec']:.1f}°")
        print(f"      R² = {dipole_ent['r_squared']:.4f}")
        print(f"      p = {dipole_ent['p_value']:.4f}")
        
        # Dipole on |c|
        dipole_c = dipole_fit(subset['ra'].values, subset['dec'].values,
                             subset['c_detrended'].values)
        print(f"    Color residual dipole:")
        print(f"      Amplitude: {dipole_c['amplitude']:.6f}")
        print(f"      Direction: RA={dipole_c['ra']:.1f}°, Dec={dipole_c['dec']:.1f}°")
        print(f"      R² = {dipole_c['r_squared']:.4f}")
        print(f"      p = {dipole_c['p_value']:.4f}")
        
        # Dipole on mu_resid
        dipole_mu = dipole_fit(subset['ra'].values, subset['dec'].values,
                              subset['mu_resid'].values)
        print(f"    Hubble residual dipole:")
        print(f"      Amplitude: {dipole_mu['amplitude']:.6f}")
        print(f"      Direction: RA={dipole_mu['ra']:.1f}°, Dec={dipole_mu['dec']:.1f}°")
        print(f"      R² = {dipole_mu['r_squared']:.4f}")
        print(f"      p = {dipole_mu['p_value']:.4f}")
        
        results[f'A1_dipole_{label.replace(" ", "_")}'] = {
            'entanglement': dipole_ent,
            'color': dipole_c,
            'hubble_resid': dipole_mu,
        }
    
    # ============================================================
    # TEST A2: QUADRUPOLE FIT
    # ============================================================
    print("\n" + "=" * 70)
    print("[A2] QUADRUPOLE FIT — Higher-order structure?")
    print("=" * 70)
    
    for label, subset in [("All z>0.1", df_signal), ("High-z (z>0.5)", df_highz)]:
        print(f"\n  {label}:")
        quad = quadrupole_fit(subset['ra'].values, subset['dec'].values,
                             subset['entanglement_proxy'].values)
        print(f"    Full model (ℓ=1+2) R² = {quad['r_squared_full']:.4f}, p = {quad['p_full_model']:.4f}")
        print(f"    Quadrupole F-test: F = {quad['f_stat_quad']:.2f}, p = {quad['p_quadrupole']:.4f}")
        
        results[f'A2_quad_{label.replace(" ", "_")}'] = quad
    
    # ============================================================
    # TEST A3: HEMISPHERE COMPARISON
    # ============================================================
    print("\n" + "=" * 70)
    print("[A3] HEMISPHERE COMPARISON — Known cosmological axes")
    print("=" * 70)
    
    # Known axes to test against
    axes = {
        'CMB_dipole': (168.0, -7.0),           # CMB dipole direction
        'CMB_cold_spot': (209.0, -57.0),        # Cold Spot direction
        'Bulk_flow': (282.0, 6.0),              # Kashlinsky bulk flow
        'Galactic_poles': (192.86, 27.13),      # North Galactic Pole
        'Ecliptic_poles': (270.0, 66.56),       # North Ecliptic Pole
    }
    
    # Also test the dipole direction we found
    if df_highz is not None and len(df_highz) > 50:
        best_dipole = dipole_fit(df_highz['ra'].values, df_highz['dec'].values,
                                df_highz['entanglement_proxy'].values)
        axes['Closure_dipole'] = (best_dipole['ra'], best_dipole['dec'])
    
    hemi_results = hemisphere_comparison(df_signal, axes)
    
    for name, res in hemi_results.items():
        sig = "**" if res['mannwhitney_p'] < 0.01 else "*" if res['mannwhitney_p'] < 0.05 else ""
        print(f"  {name:20s}: toward={res['mean_ent_toward']:+.5f}, "
              f"away={res['mean_ent_away']:+.5f}, "
              f"MW p={res['mannwhitney_p']:.4f} {sig}")
        print(f"  {'':20s}  ρ(toward)={res['rho_toward']:+.3f}, "
              f"ρ(away)={res['rho_away']:+.3f}, "
              f"Δρ={res['rho_diff']:+.3f}")
    
    results['A3_hemispheres'] = hemi_results
    
    # ============================================================
    # TEST A4: REDSHIFT-DEPENDENT ANISOTROPY
    # ============================================================
    print("\n" + "=" * 70)
    print("[A4] REDSHIFT EVOLUTION OF DIPOLE")
    print("Does the anisotropy grow/change with distance?")
    print("=" * 70)
    
    z_evo = redshift_dependent_anisotropy(df_signal)
    
    for r in z_evo:
        if 'note' in r:
            print(f"  z={r['z_range']}: N={r['n']} ({r['note']})")
        else:
            sig = "**" if r['dipole_p'] < 0.01 else "*" if r['dipole_p'] < 0.05 else ""
            print(f"  z={r['z_range']}: N={r['n']:4d}, amp={r['dipole_amp']:.5f}, "
                  f"RA={r['dipole_ra']:6.1f}°, Dec={r['dipole_dec']:+6.1f}°, "
                  f"p={r['dipole_p']:.4f} {sig}")
    
    results['A4_z_evolution'] = z_evo
    
    # ============================================================
    # TEST A5: SKY PATCH ANALYSIS
    # ============================================================
    print("\n" + "=" * 70)
    print("[A5] SKY PATCH MAP — Local signal strength")
    print("=" * 70)
    
    patches = sky_patch_analysis(df_signal)
    
    # Show patches with strongest/weakest signals
    patches_sorted = sorted([p for p in patches if p['local_rho'] is not None],
                           key=lambda p: p['local_rho'] if p['local_rho'] is not None else 0)
    
    print(f"  Total patches with enough data: {len(patches_sorted)}")
    if patches_sorted:
        print(f"\n  Strongest POSITIVE c-μ coupling (entanglement):")
        for p in patches_sorted[-3:]:
            print(f"    RA≈{p['ra_center']:5.0f}°, Dec≈{p['dec_center']:+5.0f}°: "
                  f"ρ={p['local_rho']:+.3f}, N={p['n']}, <z>={p['mean_z']:.2f}")
        
        print(f"\n  Strongest NEGATIVE c-μ coupling:")
        for p in patches_sorted[:3]:
            print(f"    RA≈{p['ra_center']:5.0f}°, Dec≈{p['dec_center']:+5.0f}°: "
                  f"ρ={p['local_rho']:+.3f}, N={p['n']}, <z>={p['mean_z']:.2f}")
        
        # Test if the variation across patches is more than expected
        rhos = [p['local_rho'] for p in patches_sorted]
        std_rhos = np.std(rhos)
        
        # Expected scatter from random: simulate
        n_per_patch = np.mean([p['n'] for p in patches_sorted])
        null_stds = []
        for _ in range(5000):
            null_rhos = []
            for _ in range(len(rhos)):
                fake_x = np.random.randn(int(n_per_patch))
                fake_y = np.random.randn(int(n_per_patch))
                null_rhos.append(stats.spearmanr(fake_x, fake_y)[0])
            null_stds.append(np.std(null_rhos))
        
        p_variation = np.mean(np.array(null_stds) >= std_rhos)
        print(f"\n  Patch-to-patch variation: σ(ρ) = {std_rhos:.4f}")
        print(f"  Expected from noise: σ = {np.mean(null_stds):.4f}")
        print(f"  p(excess variation) = {p_variation:.4f}")
        
        results['A5_patch_variation'] = {
            'std_rhos': float(std_rhos),
            'expected_std': float(np.mean(null_stds)),
            'p_excess': float(p_variation),
        }
    
    results['A5_patches'] = patches
    
    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    # Check for significant dipoles
    for label in ['All_z>0.1', 'High-z_(z>0.5)']:
        key = f'A1_dipole_{label}'
        if key in results:
            for metric in ['entanglement', 'color', 'hubble_resid']:
                d = results[key][metric]
                if d['p_value'] < 0.05:
                    print(f"  ⚠️  {label} {metric} dipole: "
                          f"RA={d['ra']:.0f}°, Dec={d['dec']:.0f}°, "
                          f"p={d['p_value']:.4f}")
    
    # Check hemisphere comparisons
    for name, res in hemi_results.items():
        if res['mannwhitney_p'] < 0.05:
            print(f"  ⚠️  Hemisphere asymmetry ({name}): p={res['mannwhitney_p']:.4f}, "
                  f"Δρ={res['rho_diff']:+.3f}")
    
    # Check z-evolution
    sig_z = [r for r in z_evo if 'dipole_p' in r and r['dipole_p'] < 0.05]
    if sig_z:
        print(f"  ⚠️  Significant dipoles in {len(sig_z)} redshift bins")
        for r in sig_z:
            print(f"      z={r['z_range']}: RA={r['dipole_ra']:.0f}°, "
                  f"Dec={r['dipole_dec']:.0f}°, p={r['dipole_p']:.4f}")
    
    # Save
    with open(RESULTS_DIR / "anisotropy_results.json", 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n  Results saved to {RESULTS_DIR}/")
    return results


if __name__ == "__main__":
    run_tests()
