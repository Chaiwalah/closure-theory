#!/usr/bin/env python3
"""
Closure Theory — Round 5: Anisotropy & Geometry Tests
======================================================
Tests whether closure effects are direction-dependent and structure-dependent.

Tests:
  S1: Sky-dependent threshold mapping — Find z* in different sky regions
  S2: Void vs filament comparison — Entanglement near voids vs dense regions
  S3: Reconstruction entropy vs z — Bits needed to reconstruct observables
  S4: DM-μ_resid deep dive — Characterize the elephant
  S5: Directional information compression — Does rank drop vary by sky direction?
"""

import numpy as np
import os
import json
from scipy import stats
from scipy.optimize import curve_fit
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

RESULTS_DIR = "results_round5"
os.makedirs(RESULTS_DIR, exist_ok=True)

# =============================================================================
# DATA LOADING
# =============================================================================
def load_pantheon():
    import pandas as pd
    df = pd.read_csv("data/pantheon_plus.dat", sep=r'\s+', comment='#')
    return df

def load_frb_crossmatch(pantheon):
    """Rebuild the FRB cross-match for DM tests."""
    import pandas as pd
    
    frb = pd.read_csv('data/chimefrbcat1.csv', sep='\t')
    frb = frb.replace('-9999', np.nan)
    for col in ['ra', 'dec', 'bonsai_dm']:
        frb[col] = pd.to_numeric(frb[col], errors='coerce')
    
    frb_valid = frb.dropna(subset=['ra', 'dec', 'bonsai_dm'])
    frb_valid = frb_valid[frb_valid['bonsai_dm'] > 0]
    
    frb_ra = frb_valid['ra'].values
    frb_dec = frb_valid['dec'].values
    frb_dm = frb_valid['bonsai_dm'].values
    
    sn_ra = pantheon['RA'].values
    sn_dec = pantheon['DEC'].values
    
    def angular_sep(ra1, dec1, ra2, dec2):
        ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
        cos_sep = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)
        return np.degrees(np.arccos(np.clip(cos_sep, -1, 1)))
    
    matched_idx = []
    matched_dm = []
    matched_sep = []
    matched_frb_dec = []
    
    for i in range(len(sn_ra)):
        if not (np.isfinite(sn_ra[i]) and np.isfinite(sn_dec[i])):
            continue
        seps = angular_sep(sn_ra[i], sn_dec[i], frb_ra, frb_dec)
        j = np.argmin(seps)
        if seps[j] <= 10.0:
            matched_idx.append(i)
            matched_dm.append(frb_dm[j])
            matched_sep.append(seps[j])
            matched_frb_dec.append(frb_dec[j])
    
    return np.array(matched_idx), np.array(matched_dm), np.array(matched_sep), np.array(matched_frb_dec)


def compute_derived(df):
    """Add derived columns."""
    z = df['zHD'].values
    
    # Color residual (detrended)
    c = df['c'].values
    mask = np.isfinite(c) & np.isfinite(z)
    poly = np.polyfit(z[mask], c[mask], 1)
    df['c_resid'] = c - np.polyval(poly, z)
    
    # x1 residual
    x1 = df['x1'].values
    mask = np.isfinite(x1) & np.isfinite(z)
    poly = np.polyfit(z[mask], x1[mask], 1)
    df['x1_resid'] = x1 - np.polyval(poly, z)
    
    # chi2/dof
    if 'NDOF' in df.columns and 'FITCHI2' in df.columns:
        df['chi2_dof'] = df['FITCHI2'] / df['NDOF'].replace(0, np.nan)
    
    # Hubble residual
    if 'MU_SH0ES' in df.columns:
        from scipy.integrate import quad
        def mu_lcdm(zv, H0=73.04, Om=0.334):
            dc, _ = quad(lambda zp: 1/np.sqrt(Om*(1+zp)**3+(1-Om)), 0, zv)
            return 5*np.log10((1+zv)*dc*2.998e5/H0) + 25
        mu_model = np.array([mu_lcdm(zi) if np.isfinite(zi) and zi > 0 else np.nan for zi in z])
        df['mu_resid'] = df['MU_SH0ES'].values - mu_model
    
    # Galactic coordinates for sky region tests
    try:
        ra_rad = np.radians(df['RA'].values)
        dec_rad = np.radians(df['DEC'].values)
        # Convert to galactic (approximate)
        # North galactic pole: RA=192.85°, DEC=27.13°
        ra_ngp = np.radians(192.85948)
        dec_ngp = np.radians(27.12825)
        l_ncp = np.radians(122.93192)
        
        sin_b = np.sin(dec_ngp)*np.sin(dec_rad) + np.cos(dec_ngp)*np.cos(dec_rad)*np.cos(ra_rad - ra_ngp)
        df['gal_b'] = np.degrees(np.arcsin(np.clip(sin_b, -1, 1)))
        
        cos_b = np.cos(np.arcsin(np.clip(sin_b, -1, 1)))
        sin_l_minus = np.cos(dec_rad)*np.sin(ra_rad - ra_ngp) / np.clip(cos_b, 1e-10, None)
        cos_l_minus = (np.sin(dec_rad) - sin_b*np.sin(dec_ngp)) / (np.clip(cos_b, 1e-10, None)*np.cos(dec_ngp))
        df['gal_l'] = np.degrees(l_ncp - np.arctan2(sin_l_minus, cos_l_minus)) % 360
    except:
        df['gal_l'] = df['RA'].values
        df['gal_b'] = df['DEC'].values
    
    return df

# =============================================================================
# TEST S1: SKY-DEPENDENT THRESHOLD MAPPING
# =============================================================================
def test_sky_threshold(df):
    """Find z* in different sky regions. If they differ → boundary not spherical."""
    print("\n" + "=" * 60)
    print("S1: SKY-DEPENDENT THRESHOLD MAPPING")
    print("=" * 60)
    
    z = df['zHD'].values
    c = df['c'].values
    mu_res = df['mu_resid'].values
    ra = df['RA'].values
    dec = df['DEC'].values
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(mu_res) & np.isfinite(ra) & np.isfinite(dec)
    z, c, mu_res, ra, dec = z[mask], c[mask], mu_res[mask], ra[mask], dec[mask]
    
    def sigmoid(x, z0, k, A, B):
        return A / (1 + np.exp(-k * (x - z0))) + B
    
    def find_threshold(z_sub, c_sub, mu_sub, label):
        """Find sigmoid threshold for a sky region."""
        z_bins = np.linspace(np.percentile(z_sub, 10), np.percentile(z_sub, 90), 20)
        window = (z_bins[-1] - z_bins[0]) / 8
        
        z_centers, abs_rhos = [], []
        for zc in z_bins:
            m = (z_sub >= zc - window) & (z_sub < zc + window)
            if m.sum() > 15:
                rho, _ = stats.spearmanr(c_sub[m], mu_sub[m])
                z_centers.append(float(zc))
                abs_rhos.append(float(abs(rho)))
        
        if len(z_centers) < 5:
            return None
        
        z_arr = np.array(z_centers)
        r_arr = np.array(abs_rhos)
        
        try:
            popt, pcov = curve_fit(sigmoid, z_arr, r_arr,
                                    p0=[0.5, 10, 0.3, 0.05],
                                    maxfev=10000,
                                    bounds=([0.05, 0.1, 0.01, -0.5], [1.5, 100, 1.0, 0.5]))
            return {'z_threshold': float(popt[0]), 
                    'z_err': float(np.sqrt(np.diag(pcov))[0]),
                    'steepness': float(popt[1]),
                    'n': int(len(z_sub)),
                    'label': label}
        except:
            return None
    
    # Split sky into regions
    # Method 1: RA quadrants
    results = {'ra_quadrants': [], 'hemisphere': [], 'galactic': []}
    
    ra_edges = [0, 90, 180, 270, 360]
    for i in range(4):
        m = (ra >= ra_edges[i]) & (ra < ra_edges[i+1])
        if m.sum() > 100:
            r = find_threshold(z[m], c[m], mu_res[m], f"RA {ra_edges[i]}-{ra_edges[i+1]}°")
            if r:
                results['ra_quadrants'].append(r)
                print(f"  RA [{ra_edges[i]:3d}°-{ra_edges[i+1]:3d}°]: z* = {r['z_threshold']:.3f} ± {r['z_err']:.3f} (n={r['n']})")
    
    # Method 2: North vs South
    for label, m in [("North (dec>0)", dec > 0), ("South (dec<0)", dec < 0)]:
        if m.sum() > 100:
            r = find_threshold(z[m], c[m], mu_res[m], label)
            if r:
                results['hemisphere'].append(r)
                print(f"  {label}: z* = {r['z_threshold']:.3f} ± {r['z_err']:.3f} (n={r['n']})")
    
    # Method 3: Galactic plane vs poles
    gal_b = df['gal_b'].values[np.isfinite(z) & np.isfinite(c) & np.isfinite(mu_res) & np.isfinite(ra) & np.isfinite(dec)]
    for label, m in [("Gal plane (|b|<30°)", np.abs(gal_b) < 30), ("Gal pole (|b|>30°)", np.abs(gal_b) >= 30)]:
        if m.sum() > 100:
            r = find_threshold(z[m], c[m], mu_res[m], label)
            if r:
                results['galactic'].append(r)
                print(f"  {label}: z* = {r['z_threshold']:.3f} ± {r['z_err']:.3f} (n={r['n']})")
    
    # Check if thresholds differ
    all_z_stars = [r['z_threshold'] for group in results.values() for r in group]
    if len(all_z_stars) >= 2:
        spread = max(all_z_stars) - min(all_z_stars)
        results['spread'] = float(spread)
        results['anisotropic'] = spread > 0.15  # threshold differs by more than 0.15
        print(f"\n  z* spread across sky: {spread:.3f}")
        print(f"  Anisotropic: {'YES' if results['anisotropic'] else 'NO'} (threshold: Δz*>0.15)")
    
    return results

# =============================================================================
# TEST S2: VOID VS FILAMENT (using DM as density proxy)
# =============================================================================
def test_void_filament(df, matched_idx, matched_dm):
    """Compare closure metrics in low-DM (void-like) vs high-DM (filament-like) sightlines."""
    print("\n" + "=" * 60)
    print("S2: VOID VS FILAMENT (DM as density proxy)")
    print("=" * 60)
    
    z = df['zHD'].values[matched_idx]
    c = df['c'].values[matched_idx]
    x1 = df['x1'].values[matched_idx]
    mu_res = df['mu_resid'].values[matched_idx]
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(x1) & np.isfinite(mu_res) & np.isfinite(matched_dm)
    z, c, x1, mu_res, dm = z[mask], c[mask], x1[mask], mu_res[mask], matched_dm[mask]
    
    # Split by DM (proxy for matter along sightline)
    dm_med = np.median(dm)
    low_dm = dm < dm_med   # void-like sightlines
    high_dm = dm >= dm_med  # filament-like sightlines
    
    results = {'dm_split': float(dm_med)}
    
    # Compare entanglement (c-mu_res coupling) in each environment
    for label, m in [("Low DM (void-like)", low_dm), ("High DM (filament-like)", high_dm)]:
        rho_c, p_c = stats.spearmanr(c[m], mu_res[m])
        rho_x1, p_x1 = stats.spearmanr(x1[m], mu_res[m])
        
        # Effective rank
        obs = np.column_stack([c[m], x1[m], mu_res[m]])
        obs_std = (obs - obs.mean(0)) / (obs.std(0) + 1e-10)
        s = np.linalg.svd(obs_std, compute_uv=False)
        p = s / s.sum()
        p = p[p > 0]
        eff_rank = float(np.exp(-np.sum(p * np.log(p))))
        
        results[label] = {
            'c_mu_rho': float(rho_c), 'c_mu_p': float(p_c),
            'x1_mu_rho': float(rho_x1), 'x1_mu_p': float(p_x1),
            'eff_rank': eff_rank,
            'n': int(m.sum()),
            'mean_dm': float(dm[m].mean()),
        }
        
        print(f"\n  {label} (n={m.sum()}, mean DM={dm[m].mean():.0f}):")
        print(f"    c→μ_res:  ρ={rho_c:.4f}, p={p_c:.4f}")
        print(f"    x1→μ_res: ρ={rho_x1:.4f}, p={p_x1:.4f}")
        print(f"    Eff rank: {eff_rank:.3f}")
    
    # Is the effect stronger in high-DM sightlines?
    r_void = abs(results["Low DM (void-like)"]['c_mu_rho'])
    r_fil = abs(results["High DM (filament-like)"]['c_mu_rho'])
    rank_void = results["Low DM (void-like)"]['eff_rank']
    rank_fil = results["High DM (filament-like)"]['eff_rank']
    
    results['stronger_in_filaments'] = r_fil > r_void
    results['more_compressed_in_filaments'] = rank_fil < rank_void
    
    print(f"\n  Coupling stronger in filaments: {results['stronger_in_filaments']}")
    print(f"  More compressed in filaments: {results['more_compressed_in_filaments']}")
    
    # Test with DM tertiles for gradient
    print("\n  DM gradient (tertiles):")
    dm_t = np.percentile(dm, [0, 33, 67, 100])
    for i in range(3):
        m = (dm >= dm_t[i]) & (dm < dm_t[i+1] + 0.01)
        if m.sum() > 20:
            rho, p = stats.spearmanr(c[m], mu_res[m])
            obs = np.column_stack([c[m], x1[m], mu_res[m]])
            obs_std = (obs - obs.mean(0)) / (obs.std(0) + 1e-10)
            s = np.linalg.svd(obs_std, compute_uv=False)
            pv = s / s.sum(); pv = pv[pv > 0]
            er = float(np.exp(-np.sum(pv * np.log(pv))))
            print(f"    DM [{dm_t[i]:.0f}-{dm_t[i+1]:.0f}] (n={m.sum()}): c→μ ρ={rho:.4f}, rank={er:.3f}")
    
    return results

# =============================================================================
# TEST S3: RECONSTRUCTION ENTROPY VS Z
# =============================================================================
def test_reconstruction_entropy(df):
    """Measure bits needed to reconstruct observables as a function of z."""
    print("\n" + "=" * 60)
    print("S3: RECONSTRUCTION ENTROPY VS REDSHIFT")
    print("=" * 60)
    
    z = df['zHD'].values
    c = df['c'].values
    x1 = df['x1'].values
    mu_res = df['mu_resid'].values
    
    # Also chi2_dof if available
    has_chi2 = 'chi2_dof' in df.columns
    if has_chi2:
        chi2 = df['chi2_dof'].values
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(x1) & np.isfinite(mu_res)
    if has_chi2:
        mask &= np.isfinite(chi2)
    
    z, c, x1, mu_res = z[mask], c[mask], x1[mask], mu_res[mask]
    if has_chi2:
        chi2 = chi2[mask]
        obs_all = np.column_stack([c, x1, mu_res, chi2])
        obs_names = ['c', 'x1', 'mu_resid', 'chi2_dof']
    else:
        obs_all = np.column_stack([c, x1, mu_res])
        obs_names = ['c', 'x1', 'mu_resid']
    
    def svd_entropy(X):
        """Shannon entropy of normalized singular values = reconstruction entropy."""
        X_std = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-10)
        s = np.linalg.svd(X_std, compute_uv=False)
        p = s / s.sum()
        p = p[p > 1e-15]
        return float(-np.sum(p * np.log2(p)))
    
    def knn_entropy(X, k=5):
        """k-NN entropy estimate (Kozachenko-Leonenko)."""
        from scipy.spatial import cKDTree
        X_std = (X - X.mean(0)) / (X.std(0) + 1e-10)
        tree = cKDTree(X_std)
        dists, _ = tree.query(X_std, k=k+1)
        eps = dists[:, -1]  # distance to k-th neighbor
        eps = eps[eps > 0]
        d = X.shape[1]
        # Kozachenko-Leonenko estimator
        n = len(eps)
        H = d * np.mean(np.log2(2 * eps)) + np.log2(n) + np.euler_gamma / np.log(2)
        return float(H)
    
    # Sliding window entropy vs z
    z_sorted = np.argsort(z)
    window = len(z) // 6
    step = window // 4
    
    results = {'svd_entropy': [], 'knn_entropy': [], 'z_centers': []}
    
    print(f"\n  Observables: {obs_names}")
    print(f"  Window size: {window}, step: {step}")
    print(f"\n  {'z_center':>10s} {'SVD H (bits)':>12s} {'kNN H':>10s} {'n':>6s}")
    print("  " + "-" * 42)
    
    for i in range(0, len(z) - window, step):
        idx = z_sorted[i:i+window]
        z_center = float(np.median(z[idx]))
        
        X = obs_all[idx]
        h_svd = svd_entropy(X)
        h_knn = knn_entropy(X)
        
        results['z_centers'].append(z_center)
        results['svd_entropy'].append(h_svd)
        results['knn_entropy'].append(h_knn)
        
        print(f"  {z_center:10.3f} {h_svd:12.4f} {h_knn:10.2f} {len(idx):6d}")
    
    # Test for trend: does entropy decrease with z?
    z_arr = np.array(results['z_centers'])
    h_svd_arr = np.array(results['svd_entropy'])
    h_knn_arr = np.array(results['knn_entropy'])
    
    rho_svd, p_svd = stats.spearmanr(z_arr, h_svd_arr)
    rho_knn, p_knn = stats.spearmanr(z_arr, h_knn_arr)
    
    results['svd_trend'] = {'rho': float(rho_svd), 'p': float(p_svd)}
    results['knn_trend'] = {'rho': float(rho_knn), 'p': float(p_knn)}
    
    print(f"\n  SVD entropy trend with z: ρ={rho_svd:.4f}, p={p_svd:.4f}")
    print(f"  kNN entropy trend with z: ρ={rho_knn:.4f}, p={p_knn:.4f}")
    
    if rho_svd < 0 and p_svd < 0.05:
        print("  → Observable space LOSES information with distance ✓")
    elif rho_knn < 0 and p_knn < 0.05:
        print("  → Observable space LOSES information with distance (kNN) ✓")
    else:
        print("  → No significant entropy decrease detected")
    
    return results

# =============================================================================
# TEST S4: DM-μ_RESID DEEP DIVE (THE ELEPHANT)
# =============================================================================
def test_dm_elephant(df, matched_idx, matched_dm, matched_sep, matched_frb_dec):
    """Deep characterization of the μ_resid-DM correlation."""
    print("\n" + "=" * 60)
    print("S4: THE ELEPHANT — DM vs HUBBLE RESIDUAL DEEP DIVE")
    print("=" * 60)
    
    z = df['zHD'].values[matched_idx]
    c = df['c'].values[matched_idx]
    x1 = df['x1'].values[matched_idx]
    mu_res = df['mu_resid'].values[matched_idx]
    host_mass = df['HOST_LOGMASS'].values[matched_idx] if 'HOST_LOGMASS' in df.columns else None
    sn_dec = df['DEC'].values[matched_idx]
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(x1) & np.isfinite(mu_res) & np.isfinite(matched_dm)
    z, c, x1, mu_res, dm = z[mask], c[mask], x1[mask], mu_res[mask], matched_dm[mask]
    sep = matched_sep[mask]
    frb_dec = matched_frb_dec[mask]
    sn_dec = sn_dec[mask]
    if host_mass is not None:
        host_mass = host_mass[mask]
    
    results = {}
    
    # 1. Raw correlation
    rho, p = stats.spearmanr(dm, mu_res)
    results['raw'] = {'rho': float(rho), 'p': float(p), 'n': int(len(z))}
    print(f"\n  Raw μ_resid vs DM: ρ={rho:.4f}, p={p:.6f}, n={len(z)}")
    
    # 2. Partial correlation controlling for EVERYTHING
    from numpy.linalg import lstsq
    controls = np.column_stack([z, z**2, frb_dec, sn_dec, sep, c, x1])
    control_names = ['z', 'z²', 'frb_dec', 'sn_dec', 'sep', 'c', 'x1']
    if host_mass is not None:
        hm = host_mass.copy()
        hm[~np.isfinite(hm)] = np.nanmedian(hm)
        controls = np.column_stack([controls, hm])
        control_names.append('host_mass')
    
    coef_mu, _, _, _ = lstsq(controls, mu_res, rcond=None)
    coef_dm, _, _, _ = lstsq(controls, dm, rcond=None)
    mu_ctrl = mu_res - controls @ coef_mu
    dm_ctrl = dm - controls @ coef_dm
    
    rho_ctrl, p_ctrl = stats.spearmanr(mu_ctrl, dm_ctrl)
    results['full_control'] = {
        'rho': float(rho_ctrl), 'p': float(p_ctrl),
        'controls': control_names
    }
    print(f"  Full control μ_resid vs DM: ρ={rho_ctrl:.4f}, p={p_ctrl:.6f}")
    print(f"    Controls: {', '.join(control_names)}")
    
    # 3. Direction of effect: positive DM → dimmer or brighter?
    dm_med = np.median(dm)
    mu_high_dm = mu_res[dm > dm_med].mean()
    mu_low_dm = mu_res[dm <= dm_med].mean()
    results['direction'] = {
        'high_dm_mean_mu': float(mu_high_dm),
        'low_dm_mean_mu': float(mu_low_dm),
        'effect': 'dimmer' if mu_high_dm > mu_low_dm else 'brighter'
    }
    print(f"\n  Direction: High-DM sightlines → SNe appear {results['direction']['effect']}")
    print(f"    <μ_res> high DM: {mu_high_dm:.4f}, low DM: {mu_low_dm:.4f}")
    
    # 4. Redshift dependence — does the DM effect change with z?
    print("\n  Redshift dependence:")
    z_bins = np.percentile(z, [0, 25, 50, 75, 100])
    results['z_bins'] = []
    for i in range(4):
        m = (z >= z_bins[i]) & (z < z_bins[i+1] + 0.001)
        if m.sum() > 20:
            rho_z, p_z = stats.spearmanr(dm[m], mu_res[m])
            entry = {'z_range': [float(z_bins[i]), float(z_bins[i+1])],
                     'rho': float(rho_z), 'p': float(p_z), 'n': int(m.sum())}
            results['z_bins'].append(entry)
            print(f"    z=[{z_bins[i]:.2f}-{z_bins[i+1]:.2f}] (n={m.sum()}): ρ={rho_z:.4f}, p={p_z:.4f}")
    
    # 5. Non-linearity — is DM effect threshold-like?
    print("\n  DM threshold scan:")
    dm_thresholds = np.percentile(dm, [20, 40, 60, 80])
    results['dm_threshold'] = []
    for thresh in dm_thresholds:
        above = dm > thresh
        if above.sum() > 20 and (~above).sum() > 20:
            t, p_t = stats.ttest_ind(mu_res[above], mu_res[~above])
            print(f"    DM>{thresh:.0f}: Δμ={mu_res[above].mean()-mu_res[~above].mean():.4f}, t={t:.2f}, p={p_t:.4f}")
            results['dm_threshold'].append({
                'dm_thresh': float(thresh),
                'delta_mu': float(mu_res[above].mean()-mu_res[~above].mean()),
                't': float(t), 'p': float(p_t)
            })
    
    # 6. Bootstrap confidence
    n_boot = 5000
    rho_boot = []
    for _ in range(n_boot):
        idx = np.random.choice(len(dm), len(dm), replace=True)
        rho_boot.append(stats.spearmanr(dm[idx], mu_res[idx])[0])
    rho_boot = np.array(rho_boot)
    results['bootstrap'] = {
        'mean_rho': float(rho_boot.mean()),
        'ci_95': [float(np.percentile(rho_boot, 2.5)), float(np.percentile(rho_boot, 97.5))],
        'frac_negative': float((rho_boot < 0).mean()),
    }
    print(f"\n  Bootstrap: ρ={rho_boot.mean():.4f} [{np.percentile(rho_boot,2.5):.4f}, {np.percentile(rho_boot,97.5):.4f}]")
    print(f"  Fraction negative: {(rho_boot < 0).mean():.1%}")
    
    return results

# =============================================================================
# TEST S5: DIRECTIONAL INFORMATION COMPRESSION
# =============================================================================
def test_directional_compression(df):
    """Does the effective rank drop vary by sky direction?"""
    print("\n" + "=" * 60)
    print("S5: DIRECTIONAL INFORMATION COMPRESSION")
    print("=" * 60)
    
    z = df['zHD'].values
    c = df['c'].values
    x1 = df['x1'].values
    mu_res = df['mu_resid'].values
    ra = df['RA'].values
    
    mask = np.isfinite(z) & np.isfinite(c) & np.isfinite(x1) & np.isfinite(mu_res) & np.isfinite(ra)
    z, c, x1, mu_res, ra = z[mask], c[mask], x1[mask], mu_res[mask], ra[mask]
    
    def effective_rank(X):
        X_std = (X - X.mean(0)) / (X.std(0) + 1e-10)
        s = np.linalg.svd(X_std, compute_uv=False)
        p = s / s.sum()
        p = p[p > 0]
        return float(np.exp(-np.sum(p * np.log(p))))
    
    z_med = np.median(z)
    
    # Split sky into 6 RA slices
    ra_edges = np.linspace(0, 360, 7)
    results = {'regions': []}
    
    print(f"\n  {'RA range':>15s} {'rank_lo':>8s} {'rank_hi':>8s} {'Δrank':>8s} {'n_lo':>6s} {'n_hi':>6s}")
    print("  " + "-" * 55)
    
    for i in range(6):
        m_ra = (ra >= ra_edges[i]) & (ra < ra_edges[i+1])
        lo = m_ra & (z < z_med)
        hi = m_ra & (z >= z_med)
        
        if lo.sum() > 30 and hi.sum() > 30:
            obs_lo = np.column_stack([c[lo], x1[lo], mu_res[lo]])
            obs_hi = np.column_stack([c[hi], x1[hi], mu_res[hi]])
            
            r_lo = effective_rank(obs_lo)
            r_hi = effective_rank(obs_hi)
            delta = r_lo - r_hi
            
            entry = {
                'ra_range': [float(ra_edges[i]), float(ra_edges[i+1])],
                'rank_low_z': r_lo, 'rank_high_z': r_hi,
                'rank_drop': delta,
                'n_low': int(lo.sum()), 'n_high': int(hi.sum())
            }
            results['regions'].append(entry)
            
            marker = " ←" if delta > 0.1 else ""
            print(f"  [{ra_edges[i]:5.0f}°-{ra_edges[i+1]:5.0f}°] {r_lo:8.3f} {r_hi:8.3f} {delta:+8.3f} {lo.sum():6d} {hi.sum():6d}{marker}")
    
    # Is there a direction with anomalously high compression?
    drops = [r['rank_drop'] for r in results['regions']]
    if drops:
        results['max_compression_direction'] = results['regions'][np.argmax(drops)]['ra_range']
        results['min_compression_direction'] = results['regions'][np.argmin(drops)]['ra_range']
        results['compression_spread'] = float(max(drops) - min(drops))
        results['anisotropic'] = results['compression_spread'] > 0.15
        
        print(f"\n  Max compression: RA {results['max_compression_direction']} (Δrank={max(drops):.3f})")
        print(f"  Min compression: RA {results['min_compression_direction']} (Δrank={min(drops):.3f})")
        print(f"  Spread: {results['compression_spread']:.3f}")
        print(f"  Anisotropic: {'YES' if results['anisotropic'] else 'NO'}")
    
    return results

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 60)
    print("CLOSURE THEORY — ROUND 5: ANISOTROPY & GEOMETRY")
    print("=" * 60)
    
    df = load_pantheon()
    print(f"Loaded {len(df)} SNe from Pantheon+")
    
    df = compute_derived(df)
    
    # Load FRB cross-match
    print("Building FRB cross-match...")
    matched_idx, matched_dm, matched_sep, matched_frb_dec = load_frb_crossmatch(df)
    print(f"Matched {len(matched_idx)} sightlines")
    
    all_results = {}
    
    # S1: Sky-dependent threshold
    all_results['S1_sky_threshold'] = test_sky_threshold(df)
    
    # S2: Void vs filament
    all_results['S2_void_filament'] = test_void_filament(df, matched_idx, matched_dm)
    
    # S3: Reconstruction entropy
    all_results['S3_entropy'] = test_reconstruction_entropy(df)
    
    # S4: DM elephant deep dive
    all_results['S4_elephant'] = test_dm_elephant(df, matched_idx, matched_dm, matched_sep, matched_frb_dec)
    
    # S5: Directional compression
    all_results['S5_directional'] = test_directional_compression(df)
    
    # Save
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.bool_,)): return bool(obj)
            if isinstance(obj, (np.integer,)): return int(obj)
            if isinstance(obj, (np.floating,)): return float(obj)
            if isinstance(obj, np.ndarray): return obj.tolist()
            return super().default(obj)
    
    with open(os.path.join(RESULTS_DIR, 'round5_results.json'), 'w') as f:
        json.dump(all_results, f, indent=2, cls=NpEncoder)
    
    # Final summary
    print("\n" + "=" * 60)
    print("ROUND 5 SUMMARY")
    print("=" * 60)
    
    s1 = all_results['S1_sky_threshold']
    if 'anisotropic' in s1:
        print(f"  S1 Sky threshold: {'ANISOTROPIC' if s1['anisotropic'] else 'Consistent'} (spread={s1.get('spread',0):.3f})")
    
    s2 = all_results['S2_void_filament']
    print(f"  S2 Void/filament: Stronger in filaments={s2.get('stronger_in_filaments')}, "
          f"More compressed={s2.get('more_compressed_in_filaments')}")
    
    s3 = all_results['S3_entropy']
    print(f"  S3 Entropy trend: SVD ρ={s3['svd_trend']['rho']:.3f} (p={s3['svd_trend']['p']:.4f}), "
          f"kNN ρ={s3['knn_trend']['rho']:.3f} (p={s3['knn_trend']['p']:.4f})")
    
    s4 = all_results['S4_elephant']
    print(f"  S4 Elephant: Raw ρ={s4['raw']['rho']:.4f} (p={s4['raw']['p']:.6f}), "
          f"Full-ctrl ρ={s4['full_control']['rho']:.4f} (p={s4['full_control']['p']:.6f})")
    
    s5 = all_results['S5_directional']
    if 'anisotropic' in s5:
        print(f"  S5 Directional: {'ANISOTROPIC' if s5['anisotropic'] else 'Isotropic'} "
              f"(spread={s5['compression_spread']:.3f})")
    
    print(f"\nResults saved to {RESULTS_DIR}/round5_results.json")

if __name__ == "__main__":
    main()
