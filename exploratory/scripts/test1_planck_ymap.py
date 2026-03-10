#!/usr/bin/env python3
"""
GPT Cornering Test #1 (v2): Planck PR4 Compton-y as foreground proxy
====================================================================
κ (total projected mass) failed. y-map traces hot electron pressure (tSZ),
which is the RIGHT baryonic proxy for our mechanism.

Tests:
1. y-quintile split: does MgII EW/FWHM ratio depend on Compton-y?
2. Sky-rotation placebos (90°, 180°, 270°): is the signal positional?
3. Partial correlation: z vs y contribution to damage
4. 2D z×y separability (GPT Test #2): grid of z-bins × y-quintiles

If y works where κ failed, it confirms: operator couples to baryonic
electron structure, not total projected mass.
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
import healpy as hp
import os
import json
from datetime import datetime

YMAP_PATH = '/tmp/planck_pr4_ymap.fits'
QUASAR_PATH = '/root/clawd/data/sdss/dr16q_prop.fits'
OUTPUT_DIR = '/root/clawd/projects/closure-theory/results_ymap'

os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_ymap():
    """Load PR4 NILC y-map (FULL column, NSIDE=2048)"""
    print("Loading Planck PR4 y-map...")
    data = hp.read_map(YMAP_PATH, field=0)  # FULL column
    print(f"  NSIDE={hp.npix2nside(len(data))}, npix={len(data)}")
    print(f"  y range: {np.nanmin(data):.2e} to {np.nanmax(data):.2e}")
    print(f"  y median: {np.nanmedian(data):.2e}")
    return data

def get_y_at_positions(ymap, ra, dec):
    """Get Compton-y values at quasar positions"""
    nside = hp.npix2nside(len(ymap))
    # Convert RA/Dec to galactic coordinates (HEALPix is in galactic)
    coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    l = coords.galactic.l.deg
    b = coords.galactic.b.deg
    # HEALPix uses colatitude (theta) and longitude (phi)
    theta = np.radians(90.0 - b)
    phi = np.radians(l)
    pix = hp.ang2pix(nside, theta, phi)
    return ymap[pix]

def load_quasars():
    """Load DR16Q quasars with MgII measurements"""
    print("Loading DR16Q quasars...")
    with fits.open(QUASAR_PATH) as hdul:
        data = hdul[1].data
        z = data['Z_DR16Q']
        ra = data['RA']
        dec = data['DEC']
        mgii = data['MGII']
        
        # MgII: index 2=EW, 4=FWHM
        ew = mgii[:, 2]
        fwhm = mgii[:, 4]
        
        # Quality cuts
        mask = (z > 0.4) & (z < 2.5) & (ew > 0) & (fwhm > 0) & (ew < 500) & (fwhm < 20000)
        
        print(f"  Total: {len(z)}, after cuts: {mask.sum()}")
        return ra[mask], dec[mask], z[mask], ew[mask], fwhm[mask]

def rotate_coordinates(ra, dec, angle_deg):
    """Rotate RA by angle for placebo test"""
    ra_rot = (ra + angle_deg) % 360.0
    return ra_rot, dec

def test_y_quintiles(y_vals, z, ew, fwhm, label="real"):
    """Split by y-quintiles, measure MgII EW/FWHM ratio in each"""
    ratio = ew / fwhm
    
    # Z-bins to control for redshift
    z_bins = [(0.5, 0.8), (0.8, 1.1), (1.1, 1.5), (1.5, 2.0)]
    
    results = {}
    for z_lo, z_hi in z_bins:
        zmask = (z >= z_lo) & (z < z_hi)
        if zmask.sum() < 500:
            continue
        
        y_z = y_vals[zmask]
        ratio_z = ratio[zmask]
        
        # Remove any NaN/inf y-values
        valid = np.isfinite(y_z) & (y_z > -1e10)
        y_z = y_z[valid]
        ratio_z = ratio_z[valid]
        
        quintiles = np.percentile(y_z, [20, 40, 60, 80])
        q_labels = ['Q1(low-y)', 'Q2', 'Q3', 'Q4', 'Q5(high-y)']
        q_edges = [-np.inf] + list(quintiles) + [np.inf]
        
        medians = []
        counts = []
        for i in range(5):
            qmask = (y_z >= q_edges[i]) & (y_z < q_edges[i+1])
            medians.append(np.median(ratio_z[qmask]))
            counts.append(qmask.sum())
        
        # Trend: Spearman of quintile index vs median ratio
        rho, p = stats.spearmanr(range(5), medians)
        
        # KS between Q1 and Q5
        q1_mask = y_z < quintiles[0]
        q5_mask = y_z >= quintiles[3]
        ks_stat, ks_p = stats.ks_2samp(ratio_z[q1_mask], ratio_z[q5_mask])
        
        key = f"z={z_lo}-{z_hi}"
        results[key] = {
            'medians': medians,
            'counts': counts,
            'trend_rho': float(rho),
            'trend_p': float(p),
            'ks_stat': float(ks_stat),
            'ks_p': float(ks_p),
            'n': int(valid.sum())
        }
        
        print(f"  {label} {key}: quintile medians = {[f'{m:.4f}' for m in medians]}")
        print(f"    trend rho={rho:.3f} p={p:.3f}, KS={ks_stat:.3f} p={ks_p:.3e}")
    
    return results

def test_partial_correlation(y_vals, z, ew, fwhm):
    """Partial correlation: does y add predictive power beyond z?"""
    ratio = np.log10(ew / fwhm)
    valid = np.isfinite(y_vals) & np.isfinite(ratio) & (y_vals > -1e10)
    
    y_v = y_vals[valid]
    z_v = z[valid]
    r_v = ratio[valid]
    
    # Simple: correlate ratio with y, controlling for z
    # Residualize both on z
    z_poly = np.column_stack([z_v, z_v**2])
    
    # Ratio residuals
    coef_r = np.linalg.lstsq(z_poly, r_v, rcond=None)[0]
    resid_r = r_v - z_poly @ coef_r
    
    # y residuals
    coef_y = np.linalg.lstsq(z_poly, y_v, rcond=None)[0]
    resid_y = y_v - z_poly @ coef_y
    
    rho, p = stats.spearmanr(resid_y, resid_r)
    
    print(f"\n  Partial correlation (y|z vs ratio|z): rho={rho:.4f} p={p:.2e}")
    
    # Also: z partial (controlling for y)
    coef_r2 = np.linalg.lstsq(y_v.reshape(-1,1), r_v, rcond=None)[0]
    resid_r2 = r_v - y_v.reshape(-1,1) @ coef_r2
    coef_z2 = np.linalg.lstsq(y_v.reshape(-1,1), z_v, rcond=None)[0]
    resid_z2 = z_v - y_v.reshape(-1,1) @ coef_z2
    
    rho_z, p_z = stats.spearmanr(resid_z2, resid_r2)
    print(f"  Partial correlation (z|y vs ratio|y): rho={rho_z:.4f} p={p_z:.2e}")
    
    return {
        'y_given_z': {'rho': float(rho), 'p': float(p)},
        'z_given_y': {'rho': float(rho_z), 'p': float(p_z)},
        'n': int(valid.sum())
    }

def test_2d_separability(y_vals, z, ew, fwhm):
    """GPT Test #2: 2D z×y grid, fit D = a·f(z) + b·y + c·f(z)·y"""
    ratio = ew / fwhm
    valid = np.isfinite(y_vals) & (y_vals > -1e10)
    
    y_v = y_vals[valid]
    z_v = z[valid]
    r_v = ratio[valid]
    
    # Grid: 5 z-bins × 5 y-quintiles
    z_edges = np.percentile(z_v, np.linspace(0, 100, 6))
    y_edges = np.percentile(y_v, np.linspace(0, 100, 6))
    
    grid = np.zeros((5, 5))
    grid_n = np.zeros((5, 5), dtype=int)
    
    for i in range(5):
        for j in range(5):
            mask = ((z_v >= z_edges[i]) & (z_v < z_edges[i+1]) &
                    (y_v >= y_edges[j]) & (y_v < y_edges[j+1]))
            if mask.sum() > 50:
                grid[i, j] = np.median(r_v[mask])
                grid_n[i, j] = mask.sum()
    
    # Test: does y-axis show variation at fixed z?
    y_variation = np.std(grid, axis=1)  # variation across y at each z
    z_variation = np.std(grid, axis=0)  # variation across z at each y
    
    print(f"\n  2D Grid (z×y) median EW/FWHM:")
    print(f"  y-axis variation per z-bin: {[f'{v:.4f}' for v in y_variation]}")
    print(f"  z-axis variation per y-bin: {[f'{v:.4f}' for v in z_variation]}")
    print(f"  Mean y-variation: {np.mean(y_variation):.4f}")
    print(f"  Mean z-variation: {np.mean(z_variation):.4f}")
    
    # Interaction test: fit additive vs multiplicative
    z_centers = [(z_edges[i]+z_edges[i+1])/2 for i in range(5)]
    y_centers = [(y_edges[j]+y_edges[j+1])/2 for j in range(5)]
    
    # Flatten for regression
    Z_flat, Y_flat, R_flat = [], [], []
    for i in range(5):
        for j in range(5):
            if grid_n[i,j] > 50:
                Z_flat.append(z_centers[i])
                Y_flat.append(y_centers[j])
                R_flat.append(grid[i,j])
    
    Z_flat = np.array(Z_flat)
    Y_flat = np.array(Y_flat)
    R_flat = np.array(R_flat)
    
    # Additive: R = a*z + b*y + c
    X_add = np.column_stack([Z_flat, Y_flat, np.ones_like(Z_flat)])
    coef_add = np.linalg.lstsq(X_add, R_flat, rcond=None)
    resid_add = R_flat - X_add @ coef_add[0]
    ss_add = np.sum(resid_add**2)
    
    # Multiplicative: R = a*z + b*y + d*z*y + c
    X_mult = np.column_stack([Z_flat, Y_flat, Z_flat*Y_flat, np.ones_like(Z_flat)])
    coef_mult = np.linalg.lstsq(X_mult, R_flat, rcond=None)
    resid_mult = R_flat - X_mult @ coef_mult[0]
    ss_mult = np.sum(resid_mult**2)
    
    print(f"  Additive SS={ss_add:.6f}, Multiplicative SS={ss_mult:.6f}")
    print(f"  Interaction coefficient: {coef_mult[0][2]:.6f}")
    
    return {
        'grid': grid.tolist(),
        'grid_n': grid_n.tolist(),
        'y_variation': y_variation.tolist(),
        'z_variation': z_variation.tolist(),
        'additive_coefs': coef_add[0].tolist(),
        'multiplicative_coefs': coef_mult[0].tolist(),
        'additive_ss': float(ss_add),
        'multiplicative_ss': float(ss_mult)
    }

def main():
    print("=" * 70)
    print("GPT Test #1 v2: Planck PR4 Compton-y as foreground proxy")
    print("=" * 70)
    
    ymap = load_ymap()
    ra, dec, z, ew, fwhm = load_quasars()
    
    # Get y at quasar positions
    print("\nExtracting y-values at quasar positions...")
    y_vals = get_y_at_positions(ymap, ra, dec)
    print(f"  y range at QSO positions: {np.nanmin(y_vals):.2e} to {np.nanmax(y_vals):.2e}")
    print(f"  y median at QSO: {np.nanmedian(y_vals):.2e}")
    
    # Save y-values for later use (small file, ~10MB)
    np.savez_compressed(os.path.join(OUTPUT_DIR, 'quasar_yvals.npz'),
                        ra=ra, dec=dec, z=z, ew=ew, fwhm=fwhm, y=y_vals)
    print("  Saved quasar y-values to quasar_yvals.npz")
    
    # TEST 1: Y-quintile split
    print("\n" + "=" * 50)
    print("TEST 1: Y-quintile split (MgII EW/FWHM ratio)")
    print("=" * 50)
    real_results = test_y_quintiles(y_vals, z, ew, fwhm, label="REAL")
    
    # TEST 2: Sky-rotation placebos
    print("\n" + "=" * 50)
    print("TEST 2: Sky-rotation placebos")
    print("=" * 50)
    placebo_results = {}
    for angle in [90, 180, 270]:
        print(f"\n  --- Rotation {angle}° ---")
        ra_rot, dec_rot = rotate_coordinates(ra, dec, angle)
        y_rot = get_y_at_positions(ymap, ra_rot, dec)
        placebo_results[f"rot{angle}"] = test_y_quintiles(y_rot, z, ew, fwhm, label=f"ROT{angle}")
    
    # TEST 3: Partial correlation
    print("\n" + "=" * 50)
    print("TEST 3: Partial correlation (y vs z)")
    print("=" * 50)
    partial_results = test_partial_correlation(y_vals, z, ew, fwhm)
    
    # TEST 4: 2D separability
    print("\n" + "=" * 50)
    print("TEST 4: 2D z×y separability")
    print("=" * 50)
    sep_results = test_2d_separability(y_vals, z, ew, fwhm)
    
    # SUMMARY
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    # Did y work where κ failed?
    any_significant = False
    for key, res in real_results.items():
        if res['ks_p'] < 0.01:
            any_significant = True
            print(f"  ✅ SIGNIFICANT: {key} KS p={res['ks_p']:.2e}")
    
    if not any_significant:
        print("  ❌ No significant y-quintile splits found")
        print("  → y-map may also be too smooth/blunt for this operator")
        print("  → Next: absorber-incidence LOS proxy (DESI catalogs)")
    
    # Check placebos
    print(f"\n  Partial correlation y|z: rho={partial_results['y_given_z']['rho']:.4f}")
    print(f"  Partial correlation z|y: rho={partial_results['z_given_y']['rho']:.4f}")
    
    # Save all results
    all_results = {
        'timestamp': datetime.now().isoformat(),
        'real_quintiles': real_results,
        'placebos': placebo_results,
        'partial_correlation': partial_results,
        'separability_2d': sep_results
    }
    
    with open(os.path.join(OUTPUT_DIR, 'ymap_test_results.json'), 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n  Results saved to {OUTPUT_DIR}/ymap_test_results.json")
    
    # Clean up y-map from /tmp
    print("\n  NOTE: Delete /tmp/planck_pr4_ymap.fits after reviewing results")

if __name__ == '__main__':
    main()
