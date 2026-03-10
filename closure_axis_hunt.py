#!/usr/bin/env python3
"""
AXIS HUNT — What's the REAL variable?

z is a "polluted angle" — it contains distance info plus systematics we 
can't separate. If we plot the same correlations against different axes,
the one that gives the SHARPEST transition is closest to the "true" variable.

Test axes:
1. z (raw spectroscopic redshift — direct measurement)
2. Proper distance d_P(z) [Mpc]
3. Comoving distance d_C(z) [Mpc]
4. Luminosity distance d_L(z) [Mpc]  
5. Lookback time t_L(z) [Gyr]
6. log(1+z) — "expansion factor"
7. LOGMBH — BH mass (model-independent observable)
8. LOGLEDD — Eddington ratio (model-independent)
9. LOGL3000 — continuum luminosity (partially model-dependent via d_L)

For each axis: fit sigmoid, measure sharpness (steepness k), R², AIC.
The axis with the best fit is the "true compass."

Also: TWEAK TESTS — vary threshold, correlation method, bin size.
"""

import numpy as np
from astropy.io import fits
from astropy.cosmology import Planck18 as cosmo
from scipy import stats, optimize
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_axis_hunt', exist_ok=True)

print("=" * 70)
print("AXIS HUNT — Finding the true variable")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']

# Extract observables
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
ciii_ew = d['CIII_BR'][:, 2]
civ_ew = d['CIV'][:, 2]

# ============================================================
# Compute distance measures
# ============================================================
print("\nComputing distance measures for z grid...")

z_grid = np.arange(0.3, 2.5, 0.05)
z_mids_grid = (z_grid[:-1] + z_grid[1:]) / 2

# Compute cosmological distances at bin midpoints
d_proper = np.array([cosmo.comoving_distance(z).value / (1+z) for z in z_mids_grid])
d_comoving = np.array([cosmo.comoving_distance(z).value for z in z_mids_grid])
d_luminosity = np.array([cosmo.luminosity_distance(z).value for z in z_mids_grid])
t_lookback = np.array([cosmo.lookback_time(z).value for z in z_mids_grid])
log1z = np.log(1 + z_mids_grid)

print(f"  z range: {z_mids_grid[0]:.2f} — {z_mids_grid[-1]:.2f}")
print(f"  d_P range: {d_proper[0]:.0f} — {d_proper[-1]:.0f} Mpc")
print(f"  d_C range: {d_comoving[0]:.0f} — {d_comoving[-1]:.0f} Mpc")
print(f"  t_L range: {t_lookback[0]:.1f} — {t_lookback[-1]:.1f} Gyr")

# ============================================================
# Compute correlations in z-bins (same for all axes)
# ============================================================
def compute_corr_series(v1, v2, z_edges):
    """Compute Spearman correlation in z-bins."""
    z_mids = (z_edges[:-1] + z_edges[1:]) / 2
    corrs, z_valid, ns = [], [], []
    bin_indices = []
    
    for i in range(len(z_edges) - 1):
        zlo, zhi = z_edges[i], z_edges[i+1]
        zmask = (z_all >= zlo) & (z_all < zhi)
        d1, d2 = v1[zmask], v2[zmask]
        valid = np.isfinite(d1) & np.isfinite(d2) & (d1 != 0) & (d2 != 0)
        if valid.sum() > 50:
            r, p = stats.spearmanr(d1[valid], d2[valid])
            corrs.append(r)
            z_valid.append(z_mids[i])
            ns.append(int(valid.sum()))
            bin_indices.append(i)
    
    return np.array(corrs), np.array(z_valid), ns, bin_indices

# Key pair: MgII↔Hβ (cleanest sigmoid)
corrs_mghb, z_valid, ns, bin_idx = compute_corr_series(mg_ew, hb_ew, z_grid)
corrs_mgciv, z_valid_civ, _, bin_idx_civ = compute_corr_series(mg_ew, civ_ew, z_grid)
corrs_mgciii, z_valid_ciii, _, bin_idx_ciii = compute_corr_series(mg_ew, ciii_ew, z_grid)

# ============================================================
# Fit sigmoid on each axis
# ============================================================
def sigmoid(x, a, b, x0, k):
    return a + b / (1 + np.exp(-k * (x - x0)))

def fit_sigmoid_on_axis(x_values, corr_values, label):
    """Fit sigmoid and return quality metrics."""
    try:
        # Normalize x to [0,1] for stable fitting
        x_norm = (x_values - x_values.min()) / (x_values.max() - x_values.min())
        
        popt, pcov = optimize.curve_fit(sigmoid, x_values, corr_values,
                                         p0=[0.8, -0.8, np.median(x_values), -5],
                                         maxfev=20000)
        
        pred = sigmoid(x_values, *popt)
        ss_res = np.sum((corr_values - pred)**2)
        ss_tot = np.sum((corr_values - np.mean(corr_values))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        n = len(corr_values)
        aic = n * np.log(ss_res / n) + 2 * 4
        
        # Residual autocorrelation (lower = better fit)
        residuals = corr_values - pred
        if len(residuals) > 3:
            autocorr = np.corrcoef(residuals[:-1], residuals[1:])[0, 1]
        else:
            autocorr = 0
        
        return {
            'label': label,
            'R2': float(r2),
            'AIC': float(aic),
            'x0': float(popt[2]),
            'k': float(popt[3]),
            'residual_autocorr': float(autocorr),
            'residual_std': float(np.std(residuals)),
        }
    except Exception as e:
        return {
            'label': label,
            'R2': -999,
            'AIC': 999,
            'x0': None,
            'k': None,
            'residual_autocorr': None,
            'residual_std': None,
            'error': str(e),
        }

# Map each axis to the MgII↔Hβ z-bins
axes = {
    'z (raw)': z_valid,
    'Proper distance [Mpc]': d_proper[bin_idx],
    'Comoving distance [Mpc]': d_comoving[bin_idx],
    'Luminosity distance [Mpc]': d_luminosity[bin_idx],
    'Lookback time [Gyr]': t_lookback[bin_idx],
    'log(1+z)': log1z[bin_idx],
}

print(f"\n{'='*70}")
print(f"MgII↔Hβ SIGMOID FIT ON DIFFERENT AXES ({len(corrs_mghb)} bins)")
print(f"{'='*70}")
print(f"\n  {'Axis':<30} {'R²':>8} {'AIC':>10} {'x₀':>12} {'k':>8} {'Resid AC':>10} {'Resid σ':>10}")
print(f"  {'-'*30} {'-'*8} {'-'*10} {'-'*12} {'-'*8} {'-'*10} {'-'*10}")

results_mghb = []
for axis_name, x_values in axes.items():
    result = fit_sigmoid_on_axis(x_values, corrs_mghb, axis_name)
    results_mghb.append(result)
    
    x0_str = f"{result['x0']:.3f}" if result['x0'] is not None else "FAIL"
    k_str = f"{result['k']:.1f}" if result['k'] is not None else "FAIL"
    ac_str = f"{result['residual_autocorr']:.3f}" if result['residual_autocorr'] is not None else "FAIL"
    rs_str = f"{result['residual_std']:.4f}" if result['residual_std'] is not None else "FAIL"
    
    print(f"  {axis_name:<30} {result['R2']:>8.4f} {result['AIC']:>10.1f} {x0_str:>12} {k_str:>8} {ac_str:>10} {rs_str:>10}")

# Find winner
valid_results = [r for r in results_mghb if r['AIC'] != 999]
if valid_results:
    winner = min(valid_results, key=lambda r: r['AIC'])
    lowest_ac = min(valid_results, key=lambda r: abs(r['residual_autocorr']) if r['residual_autocorr'] is not None else 999)
    print(f"\n  ★ BEST FIT (lowest AIC): {winner['label']}")
    print(f"  ★ LOWEST RESIDUAL AUTOCORR: {lowest_ac['label']} ({lowest_ac['residual_autocorr']:.3f})")

# ============================================================
# Now try INTERNAL (model-free) axes
# ============================================================
print(f"\n{'='*70}")
print("MODEL-FREE AXES — Binning by observable, not distance")
print(f"{'='*70}")

# Instead of z-bins, bin by BH mass and check if MgII↔Hβ varies
mbh = d['LOGMBH']
ledd = d['LOGLEDD_RATIO']
l3000 = d['LOGL3000']

internal_axes = {
    'LOGMBH': mbh,
    'LOGLEDD_RATIO': ledd,
    'LOGL3000': l3000,
}

for axis_name, axis_data in internal_axes.items():
    # Restrict to z < 1.2 where both MgII and Hβ are measurable
    z_cut = (z_all >= 0.3) & (z_all < 1.2)
    
    valid_all = z_cut & np.isfinite(mg_ew) & np.isfinite(hb_ew) & \
                (mg_ew != 0) & (hb_ew != 0) & np.isfinite(axis_data) & (axis_data != 0)
    
    if valid_all.sum() < 500:
        continue
    
    # Bin by this axis
    axis_vals = axis_data[valid_all]
    percentiles = np.percentile(axis_vals, np.arange(0, 101, 10))
    
    print(f"\n  {axis_name} (z < 1.2, N={valid_all.sum()}):")
    print(f"  {'Bin':<25} {'N':>6} {'MgII↔Hβ':>10} {'mean_z':>8}")
    
    for i in range(len(percentiles) - 1):
        lo, hi = percentiles[i], percentiles[i+1]
        bin_mask = valid_all & (axis_data >= lo) & (axis_data < hi)
        n = bin_mask.sum()
        
        if n > 30:
            r, p = stats.spearmanr(mg_ew[bin_mask], hb_ew[bin_mask])
            mean_z = np.mean(z_all[bin_mask])
            print(f"  [{lo:.2f}, {hi:.2f}){'':<{max(0, 15-len(f'[{lo:.2f}, {hi:.2f})'))}} {n:>6} {r:>+10.3f} {mean_z:>8.2f}")

# ============================================================
# TWEAK TESTS — How robust are the results?
# ============================================================
print(f"\n{'='*70}")
print("TWEAK TESTS — Robustness checks")
print(f"{'='*70}")

# Tweak 1: Different correlation methods
print(f"\n  1. CORRELATION METHOD (MgII↔Hβ at z=1.0-1.1):")
zmask_test = (z_all >= 1.0) & (z_all < 1.1)
v1 = mg_ew[zmask_test]
v2 = hb_ew[zmask_test]
valid = np.isfinite(v1) & np.isfinite(v2) & (v1 != 0) & (v2 != 0)

if valid.sum() > 30:
    r_spear, _ = stats.spearmanr(v1[valid], v2[valid])
    r_pears, _ = stats.pearsonr(v1[valid], v2[valid])
    r_kendall, _ = stats.kendalltau(v1[valid], v2[valid])
    print(f"    Spearman:  {r_spear:+.4f}")
    print(f"    Pearson:   {r_pears:+.4f}")
    print(f"    Kendall τ: {r_kendall:+.4f}")

# Tweak 2: Different bin sizes
print(f"\n  2. BIN SIZE (MgII↔Hβ sigmoid z₀):")
for dz in [0.03, 0.05, 0.08, 0.10, 0.15, 0.20]:
    z_edges_test = np.arange(0.3, 2.5, dz)
    corrs_test, z_test, _, _ = compute_corr_series(mg_ew, hb_ew, z_edges_test)
    if len(corrs_test) > 5:
        try:
            popt, _ = optimize.curve_fit(sigmoid, z_test, corrs_test,
                                          p0=[0.8, -0.8, 1.0, -10], maxfev=10000)
            pred = sigmoid(z_test, *popt)
            ss_res = np.sum((corrs_test - pred)**2)
            ss_tot = np.sum((corrs_test - np.mean(corrs_test))**2)
            r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
            print(f"    Δz={dz:.2f}: z₀={popt[2]:.3f}  k={popt[3]:.1f}  R²={r2:.4f}  ({len(corrs_test)} bins)")
        except:
            print(f"    Δz={dz:.2f}: FIT FAILED")

# Tweak 3: Different correlation thresholds for graph percolation
print(f"\n  3. THRESHOLD SENSITIVITY (giant component at z=1.05-1.15):")
zmask_trans = (z_all >= 1.05) & (z_all < 1.15)
for thresh in [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5]:
    # Quick graph for all observables
    obs_list = [mg_ew, hb_ew, ciii_ew, civ_ew, d['LOGL3000'], d['LOGMBH'], d['LOGLEDD_RATIO']]
    obs_names = ['MgII', 'Hβ', 'CIII', 'CIV', 'L3000', 'MBH', 'Ledd']
    n_obs = len(obs_list)
    adj = np.zeros((n_obs, n_obs), dtype=int)
    
    for i in range(n_obs):
        for j in range(i+1, n_obs):
            vi = obs_list[i][zmask_trans]
            vj = obs_list[j][zmask_trans]
            valid = np.isfinite(vi) & np.isfinite(vj) & (vi != 0) & (vj != 0)
            if valid.sum() > 30:
                r, _ = stats.spearmanr(vi[valid], vj[valid])
                if abs(r) > thresh:
                    adj[i, j] = 1
                    adj[j, i] = 1
    
    n_edges = adj.sum() // 2
    max_edges = n_obs * (n_obs - 1) // 2
    print(f"    thresh={thresh:.2f}: {n_edges}/{max_edges} edges ({n_edges/max_edges:.0%})")

# Tweak 4: Bootstrap confidence intervals on z₀
print(f"\n  4. BOOTSTRAP z₀ (MgII↔Hβ, 100 iterations):")
z_edges_boot = np.arange(0.3, 2.5, 0.05)
z0_samples = []

for boot in range(100):
    # Resample quasars with replacement
    np.random.seed(boot)
    n_total = len(z_all)
    idx = np.random.choice(n_total, n_total, replace=True)
    
    corrs_boot = []
    z_boot = []
    
    for i in range(len(z_edges_boot) - 1):
        zlo, zhi = z_edges_boot[i], z_edges_boot[i+1]
        zmask = (z_all[idx] >= zlo) & (z_all[idx] < zhi)
        if zmask.sum() > 50:
            v1 = mg_ew[idx[zmask]]
            v2 = hb_ew[idx[zmask]]
            valid = np.isfinite(v1) & np.isfinite(v2) & (v1 != 0) & (v2 != 0)
            if valid.sum() > 30:
                r, _ = stats.spearmanr(v1[valid], v2[valid])
                corrs_boot.append(r)
                z_boot.append((zlo + zhi) / 2)
    
    if len(corrs_boot) > 5:
        try:
            popt, _ = optimize.curve_fit(sigmoid, np.array(z_boot), np.array(corrs_boot),
                                          p0=[0.8, -0.8, 1.0, -10], maxfev=5000)
            z0_samples.append(popt[2])
        except:
            pass

if z0_samples:
    z0_arr = np.array(z0_samples)
    print(f"    z₀ = {np.mean(z0_arr):.4f} ± {np.std(z0_arr):.4f}")
    print(f"    95% CI: [{np.percentile(z0_arr, 2.5):.4f}, {np.percentile(z0_arr, 97.5):.4f}]")
    print(f"    Range: [{z0_arr.min():.4f}, {z0_arr.max():.4f}]")
    print(f"    Successful fits: {len(z0_samples)}/100")

# Save
with open('results_axis_hunt/axis_hunt_results.json', 'w') as f:
    json.dump({
        'mghb_fits': results_mghb,
    }, f, indent=2, default=str)

print(f"\n\nResults saved to results_axis_hunt/")
