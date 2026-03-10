#!/usr/bin/env python3
"""
CLOSURE MASS STEP EVOLUTION
=============================
The SN Ia "mass step" — a ~0.04 mag offset between SNe in high-mass
vs low-mass host galaxies — is one of the biggest unresolved systematics
in SN cosmology. Nobody knows what causes it.

Closure prediction:
    The mass step is a DIAGNOSTIC systematic. High-mass hosts have
    different environments (older populations, more processed ISM,
    different dust). If channel degradation affects diagnostic
    observables, the mass step should:
    
    1. EVOLVE with z (grow or change character above z₀)
    2. CORRELATE with β(z) (both driven by diagnostic degradation)
    3. Be ABSENT in locked observables (stretch less affected than color)
    4. Depend on ENVIRONMENT, not just mass (density > mass)

This test:
    - Measures the mass step in z-bins
    - Tests for z-evolution of the step
    - Decomposes step into color-driven vs stretch-driven components
    - Tests correlation between mass step and β evolution
    - Checks if the step tracks host properties (diagnostic prediction)

Data: Pantheon+ (1590 SNe with host galaxy masses)

Author: Closure Theory collaboration
Date: 2026-03-03
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.stats import spearmanr, pearsonr, mannwhitneyu, ks_2samp
from scipy.interpolate import interp1d
import json
import os

# ============================================================
# CONSTANTS & COSMOLOGY
# ============================================================
c_light = 299792.458
H0_fid = 67.4

def E(z, Om=0.315):
    return np.sqrt(Om * (1+z)**3 + (1-Om))

_mu_interp = None
def mu_theory(z_arr, Om=0.315):
    global _mu_interp
    key = round(Om, 5)
    if _mu_interp is None or _mu_interp[0] != key:
        z_grid = np.linspace(0.005, 2.5, 500)
        DH = c_light / H0_fid
        mu_grid = []
        for zi in z_grid:
            dc, _ = quad(lambda zp: 1.0/E(zp, Om), 0, zi)
            dL = (1+zi) * DH * dc
            mu_grid.append(5*np.log10(dL) + 25.0)
        _mu_interp = (key, interp1d(z_grid, mu_grid, kind='cubic', fill_value='extrapolate'))
    return _mu_interp[1](z_arr)

def sigmoid(z, z0=0.82, k=8.0):
    return 1.0 / (1.0 + np.exp(-k * (z - z0)))

# ============================================================
# LOAD DATA
# ============================================================

def load_data():
    data_path = os.path.join(os.path.dirname(__file__), 'data', 'pantheon_plus.dat')
    with open(data_path, 'r') as f:
        header = f.readline().strip().split()
        rows = [line.strip().split() for line in f if len(line.strip().split()) >= len(header)]
    
    col = {name: i for i, name in enumerate(header)}
    
    z = np.array([float(r[col['zHD']]) for r in rows])
    mB = np.array([float(r[col['mB']]) for r in rows])
    c = np.array([float(r[col['c']]) for r in rows])
    x1 = np.array([float(r[col['x1']]) for r in rows])
    c_err = np.array([float(r[col['cERR']]) for r in rows])
    x1_err = np.array([float(r[col['x1ERR']]) for r in rows])
    mB_err = np.array([float(r[col['m_b_corr_err_DIAG']]) for r in rows])
    mu_corr = np.array([float(r[col['m_b_corr']]) for r in rows])
    host_mass = np.array([float(r[col['HOST_LOGMASS']]) for r in rows])
    
    mask = (z > 0.01) & (mB_err > 0) & (mB_err < 5) & np.isfinite(mB) & (host_mass > 0)
    
    d = {k: v[mask] for k, v in {
        'z': z, 'mB': mB, 'c': c, 'x1': x1, 'c_err': c_err,
        'x1_err': x1_err, 'mB_err': mB_err, 'mu_corr': mu_corr,
        'host_mass': host_mass
    }.items()}
    
    print(f"Loaded {mask.sum()} SNe with host mass data")
    print(f"  Mass range: {d['host_mass'].min():.1f} to {d['host_mass'].max():.1f} log(M/M☉)")
    print(f"  High-mass (>10): {(d['host_mass']>10).sum()}")
    print(f"  Low-mass (≤10): {(d['host_mass']<=10).sum()}")
    return d

# ============================================================
# TRIPP FORMULA HELPERS
# ============================================================

def fit_tripp(z, mB, c, x1, mB_err, c_err, x1_err, Om=0.315):
    """Fit standard Tripp: α, β, M"""
    mu_th = mu_theory(z, Om)
    
    def chi2(params):
        a, b, M = params
        mu_obs = mB + a*x1 - b*c - M
        err = np.sqrt(mB_err**2 + (a*x1_err)**2 + (b*c_err)**2)
        return np.sum(((mu_obs - mu_th)/err)**2)
    
    res = minimize(chi2, [0.14, 3.0, -19.4],
                   bounds=[(0, 0.5), (0.5, 6), (-25, -15)],
                   method='L-BFGS-B')
    return res.x  # alpha, beta, M

def hubble_residuals(d, alpha, beta, M, Om=0.315):
    """Compute Hubble residuals after Tripp correction"""
    mu_obs = d['mB'] + alpha*d['x1'] - beta*d['c'] - M
    mu_th = mu_theory(d['z'], Om)
    return mu_obs - mu_th

# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():
    print("=" * 70)
    print("CLOSURE MASS STEP EVOLUTION")
    print("Does the host mass step evolve with redshift?")
    print("=" * 70)
    
    d = load_data()
    N = len(d['z'])
    
    # Global Tripp fit
    alpha, beta, M = fit_tripp(d['z'], d['mB'], d['c'], d['x1'],
                                d['mB_err'], d['c_err'], d['x1_err'])
    hr = hubble_residuals(d, alpha, beta, M)
    
    print(f"\nGlobal Tripp: α={alpha:.4f}, β={beta:.4f}, M={M:.4f}")
    
    high_mass = d['host_mass'] > 10.0
    low_mass = ~high_mass
    
    global_step = np.median(hr[high_mass]) - np.median(hr[low_mass])
    print(f"Global mass step: {global_step:+.4f} mag")
    
    # ============================================================
    # TEST 1: Mass step in z-bins
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 1: MASS STEP vs REDSHIFT")
    print("Closure prediction: step evolves (grows/changes) above z₀=0.82")
    print("=" * 70)
    
    z_edges = [0.01, 0.08, 0.15, 0.25, 0.4, 0.6, 0.82, 1.1, 2.3]
    
    steps = []
    step_errs = []
    z_centers = []
    n_high_arr = []
    n_low_arr = []
    p_vals = []
    
    print(f"\n  {'z-bin':>12s}  {'N_hi':>5s}  {'N_lo':>5s}  {'Step (mag)':>11s}  {'±':>7s}  {'p(MW)':>8s}  {'z>z₀':>6s}")
    print("  " + "-" * 65)
    
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        zmask = (d['z'] >= z_lo) & (d['z'] < z_hi)
        
        hi = zmask & high_mass
        lo = zmask & low_mass
        n_hi, n_lo = hi.sum(), lo.sum()
        
        if n_hi < 5 or n_lo < 5:
            continue
        
        step = np.median(hr[hi]) - np.median(hr[lo])
        # Bootstrap error
        boot_steps = []
        for _ in range(1000):
            idx_hi = np.random.choice(np.where(hi)[0], n_hi, replace=True)
            idx_lo = np.random.choice(np.where(lo)[0], n_lo, replace=True)
            boot_steps.append(np.median(hr[idx_hi]) - np.median(hr[idx_lo]))
        step_err = np.std(boot_steps)
        
        # Mann-Whitney U test
        stat, p_mw = mannwhitneyu(hr[hi], hr[lo], alternative='two-sided')
        
        z_mid = np.mean(d['z'][zmask])
        above = z_mid > 0.82
        
        steps.append(step)
        step_errs.append(step_err)
        z_centers.append(z_mid)
        n_high_arr.append(n_hi)
        n_low_arr.append(n_lo)
        p_vals.append(p_mw)
        
        sig = '***' if p_mw < 0.001 else '**' if p_mw < 0.01 else '*' if p_mw < 0.05 else ''
        print(f"  [{z_lo:.2f}, {z_hi:.2f})  {n_hi:5d}  {n_lo:5d}  {step:+11.4f}  {step_err:7.4f}  {p_mw:8.4f}{sig}  {'YES ←' if above else ''}")
    
    steps = np.array(steps)
    step_errs = np.array(step_errs)
    z_centers = np.array(z_centers)
    
    # Trend
    rho_step_z, p_step_z = spearmanr(z_centers, steps)
    print(f"\n  Mass step vs z: ρ = {rho_step_z:+.3f}, p = {p_step_z:.4f}")
    
    # Absolute step (magnitude of effect)
    rho_abs, p_abs = spearmanr(z_centers, np.abs(steps))
    print(f"  |Mass step| vs z: ρ = {rho_abs:+.3f}, p = {p_abs:.4f}")
    
    # Split at z₀
    below_z0 = z_centers < 0.82
    above_z0 = z_centers >= 0.82
    
    if below_z0.sum() > 0 and above_z0.sum() > 0:
        mean_step_low = np.mean(steps[below_z0])
        mean_step_high = np.mean(steps[above_z0])
        print(f"\n  <step> at z < 0.82: {mean_step_low:+.4f} mag")
        print(f"  <step> at z ≥ 0.82: {mean_step_high:+.4f} mag")
        print(f"  Change: {mean_step_high - mean_step_low:+.4f} mag")
    
    # ============================================================
    # TEST 2: COLOR-DRIVEN vs STRETCH-DRIVEN decomposition
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 2: DECOMPOSE MASS STEP — COLOR vs STRETCH COMPONENTS")
    print("Closure: color component should evolve, stretch should be stable")
    print("=" * 70)
    
    # Color difference between high/low mass hosts
    print(f"\n  {'z-bin':>12s}  {'Δc (hi-lo)':>11s}  {'Δx1 (hi-lo)':>12s}  {'Color step':>11s}  {'Stretch step':>13s}")
    print("  " + "-" * 70)
    
    color_steps = []
    stretch_steps = []
    z_decomp = []
    
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        zmask = (d['z'] >= z_lo) & (d['z'] < z_hi)
        hi = zmask & high_mass
        lo = zmask & low_mass
        
        if hi.sum() < 5 or lo.sum() < 5:
            continue
        
        # Mean color/stretch difference
        dc = np.mean(d['c'][hi]) - np.mean(d['c'][lo])
        dx1 = np.mean(d['x1'][hi]) - np.mean(d['x1'][lo])
        
        # Color contribution to mass step: -β·Δc
        color_step = -beta * dc
        # Stretch contribution: α·Δx1
        stretch_step = alpha * dx1
        
        z_mid = np.mean(d['z'][zmask])
        color_steps.append(color_step)
        stretch_steps.append(stretch_step)
        z_decomp.append(z_mid)
        
        above = z_mid > 0.82
        print(f"  [{z_lo:.2f}, {z_hi:.2f})  {dc:+11.4f}  {dx1:+12.4f}  {color_step:+11.4f}  {stretch_step:+13.4f}  {'← z>z₀' if above else ''}")
    
    color_steps = np.array(color_steps)
    stretch_steps = np.array(stretch_steps)
    z_decomp = np.array(z_decomp)
    
    if len(color_steps) > 3:
        rho_color_z, p_color_z = spearmanr(z_decomp, color_steps)
        rho_stretch_z, p_stretch_z = spearmanr(z_decomp, stretch_steps)
        
        print(f"\n  Color component vs z: ρ = {rho_color_z:+.3f}, p = {p_color_z:.4f}")
        print(f"  Stretch component vs z: ρ = {rho_stretch_z:+.3f}, p = {p_stretch_z:.4f}")
        
        if abs(rho_color_z) > abs(rho_stretch_z):
            print(f"  ✓ Color evolves MORE than stretch → diagnostic > locked")
        else:
            print(f"  Stretch evolves more — unexpected by simple closure model")
    
    # ============================================================
    # TEST 3: β(z) × MASS STEP CORRELATION
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 3: β(z) AND MASS STEP CORRELATION")
    print("Closure: both driven by same channel degradation → correlated")
    print("=" * 70)
    
    # Compute β in same z-bins
    betas_binned = []
    steps_matched = []
    z_matched = []
    
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        zmask = (d['z'] >= z_lo) & (d['z'] < z_hi)
        n = zmask.sum()
        if n < 15:
            continue
        
        # Fit β in this bin
        d_bin = {k: v[zmask] for k, v in d.items()}
        try:
            a_b, b_b, m_b = fit_tripp(d_bin['z'], d_bin['mB'], d_bin['c'], d_bin['x1'],
                                       d_bin['mB_err'], d_bin['c_err'], d_bin['x1_err'])
        except:
            continue
        
        # Mass step in this bin
        hi = zmask & high_mass
        lo = zmask & low_mass
        if hi.sum() < 5 or lo.sum() < 5:
            continue
        
        step = np.median(hr[hi]) - np.median(hr[lo])
        z_mid = np.mean(d['z'][zmask])
        
        betas_binned.append(b_b)
        steps_matched.append(step)
        z_matched.append(z_mid)
    
    betas_binned = np.array(betas_binned)
    steps_matched = np.array(steps_matched)
    z_matched = np.array(z_matched)
    
    if len(betas_binned) > 3:
        rho_beta_step, p_beta_step = spearmanr(betas_binned, steps_matched)
        print(f"\n  β vs mass step: ρ = {rho_beta_step:+.3f}, p = {p_beta_step:.4f}")
        
        print(f"\n  {'<z>':>6s}  {'β':>7s}  {'Step':>8s}")
        print("  " + "-" * 25)
        for z, b, s in zip(z_matched, betas_binned, steps_matched):
            print(f"  {z:6.3f}  {b:7.3f}  {s:+8.4f}")
        
        if abs(rho_beta_step) > 0.5:
            print(f"\n  ✓ β and mass step are {'positively' if rho_beta_step > 0 else 'inversely'} correlated")
            print(f"    → Consistent with shared diagnostic degradation driver")
    
    # ============================================================
    # TEST 4: Color distributions in high vs low mass hosts
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 4: COLOR DISTRIBUTIONS — HIGH vs LOW MASS HOSTS")
    print("Do color distributions diverge with z? (diagnostic degradation)")
    print("=" * 70)
    
    print(f"\n  {'z-bin':>12s}  {'<c>_hi':>8s}  {'<c>_lo':>8s}  {'σ_c_hi':>7s}  {'σ_c_lo':>7s}  {'KS p':>8s}  {'Diverge?':>9s}")
    print("  " + "-" * 70)
    
    ks_p_vals = []
    z_ks = []
    
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        zmask = (d['z'] >= z_lo) & (d['z'] < z_hi)
        hi = zmask & high_mass
        lo = zmask & low_mass
        
        if hi.sum() < 10 or lo.sum() < 10:
            continue
        
        c_hi = d['c'][hi]
        c_lo = d['c'][lo]
        
        ks_stat, ks_p = ks_2samp(c_hi, c_lo)
        
        z_mid = np.mean(d['z'][zmask])
        ks_p_vals.append(ks_p)
        z_ks.append(z_mid)
        
        diverge = 'YES' if ks_p < 0.05 else 'no'
        print(f"  [{z_lo:.2f}, {z_hi:.2f})  {np.mean(c_hi):+8.4f}  {np.mean(c_lo):+8.4f}  {np.std(c_hi):7.4f}  {np.std(c_lo):7.4f}  {ks_p:8.4f}  {diverge}")
    
    if len(ks_p_vals) > 3:
        rho_ks, p_ks_trend = spearmanr(z_ks, -np.log10(np.array(ks_p_vals) + 1e-10))
        print(f"\n  KS significance vs z: ρ = {rho_ks:+.3f}, p = {p_ks_trend:.4f}")
        if rho_ks > 0:
            print(f"  ✓ Color distributions DIVERGE more at high z")
    
    # ============================================================
    # TEST 5: MASS STEP with z-dependent β correction
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 5: DOES β(z) CORRECTION REDUCE THE MASS STEP?")
    print("If mass step is partly caused by constant-β assumption...")
    print("=" * 70)
    
    # Recompute residuals with z-dependent β
    beta0 = 2.88
    dbeta = 1.16  # From coupled model
    beta_z = beta0 - dbeta * sigmoid(d['z'], 0.82, 8.0)
    
    hr_closure = d['mB'] + alpha*d['x1'] - beta_z*d['c'] - M - mu_theory(d['z'])
    
    print(f"\n  {'z-bin':>12s}  {'Standard step':>14s}  {'Closure step':>13s}  {'Reduction':>10s}")
    print("  " + "-" * 55)
    
    reductions = []
    z_red = []
    
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        zmask = (d['z'] >= z_lo) & (d['z'] < z_hi)
        hi = zmask & high_mass
        lo = zmask & low_mass
        
        if hi.sum() < 5 or lo.sum() < 5:
            continue
        
        step_std = np.median(hr[hi]) - np.median(hr[lo])
        step_cl = np.median(hr_closure[hi]) - np.median(hr_closure[lo])
        
        z_mid = np.mean(d['z'][zmask])
        reduction = (1 - abs(step_cl)/abs(step_std))*100 if abs(step_std) > 0.001 else 0
        reductions.append(reduction)
        z_red.append(z_mid)
        
        print(f"  [{z_lo:.2f}, {z_hi:.2f})  {step_std:+14.4f}  {step_cl:+13.4f}  {reduction:+9.0f}%  {'← z>z₀' if z_mid > 0.82 else ''}")
    
    if len(reductions) > 0:
        above_z0_mask = np.array(z_red) > 0.82
        if above_z0_mask.any():
            avg_red_high = np.mean(np.array(reductions)[above_z0_mask])
            avg_red_low = np.mean(np.array(reductions)[~above_z0_mask])
            print(f"\n  Mean reduction z < 0.82: {avg_red_low:+.0f}%")
            print(f"  Mean reduction z ≥ 0.82: {avg_red_high:+.0f}%")
            if avg_red_high > avg_red_low:
                print(f"  ✓ β(z) correction reduces mass step MORE at high-z")
    
    # ============================================================
    # TEST 6: MASS-COLOR-REDSHIFT 3D STRUCTURE
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 6: 3D STRUCTURE — MASS × COLOR × REDSHIFT")
    print("The full diagnostic degradation fingerprint")
    print("=" * 70)
    
    # Split into 4 quadrants: high/low mass × red/blue
    red = d['c'] > 0
    blue = d['c'] <= 0
    
    quadrants = {
        'High-mass, Red':  high_mass & red,
        'High-mass, Blue': high_mass & blue,
        'Low-mass, Red':   low_mass & red,
        'Low-mass, Blue':  low_mass & blue,
    }
    
    print(f"\n  {'Quadrant':>20s}  {'N':>5s}  {'<z>':>6s}  {'<HR>':>8s}  {'HR vs z slope':>14s}  {'p':>8s}")
    print("  " + "-" * 70)
    
    for name, mask in quadrants.items():
        n = mask.sum()
        if n < 10:
            continue
        
        mean_z = np.mean(d['z'][mask])
        mean_hr = np.mean(hr[mask])
        
        # HR vs z slope in this quadrant
        rho_hz, p_hz = spearmanr(d['z'][mask], hr[mask])
        
        print(f"  {name:>20s}  {n:5d}  {mean_z:6.3f}  {mean_hr:+8.4f}  {rho_hz:+14.3f}  {p_hz:8.4f}")
    
    # The key closure prediction: high-mass RED should diverge most at high-z
    # (color diagnostic in dense environment = maximum degradation)
    print(f"\n  Closure prediction: 'High-mass, Red' quadrant should show")
    print(f"  strongest HR-z evolution (most diagnostic degradation)")
    
    # ============================================================
    # TEST 7: CONTINUOUS MASS-STEP FUNCTION
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 7: CONTINUOUS MASS FUNCTION — IS 10.0 THE RIGHT SPLIT?")
    print("Does the mass threshold evolve with z?")
    print("=" * 70)
    
    # Find optimal mass split in z-bins
    print(f"\n  {'z-bin':>12s}  {'Best M_split':>12s}  {'Step at split':>14s}  {'Step at 10.0':>13s}")
    print("  " + "-" * 55)
    
    optimal_splits = []
    z_splits = []
    
    for i in range(len(z_edges)-1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        zmask = (d['z'] >= z_lo) & (d['z'] < z_hi)
        n = zmask.sum()
        if n < 30:
            continue
        
        # Scan mass threshold
        best_step = 0
        best_mass = 10.0
        best_sig = 0
        
        for m_try in np.arange(9.0, 11.5, 0.1):
            hi_try = zmask & (d['host_mass'] > m_try)
            lo_try = zmask & (d['host_mass'] <= m_try)
            if hi_try.sum() < 5 or lo_try.sum() < 5:
                continue
            step_try = np.median(hr[hi_try]) - np.median(hr[lo_try])
            # Significance
            try:
                _, p_try = mannwhitneyu(hr[hi_try], hr[lo_try])
                sig_try = -np.log10(max(p_try, 1e-10))
            except:
                sig_try = 0
            
            if sig_try > best_sig:
                best_sig = sig_try
                best_step = step_try
                best_mass = m_try
        
        # Step at standard M=10.0
        hi_10 = zmask & (d['host_mass'] > 10.0)
        lo_10 = zmask & (d['host_mass'] <= 10.0)
        step_10 = np.median(hr[hi_10]) - np.median(hr[lo_10]) if hi_10.sum() > 3 and lo_10.sum() > 3 else 0
        
        z_mid = np.mean(d['z'][zmask])
        optimal_splits.append(best_mass)
        z_splits.append(z_mid)
        
        print(f"  [{z_lo:.2f}, {z_hi:.2f})  {best_mass:12.1f}  {best_step:+14.4f}  {step_10:+13.4f}")
    
    if len(optimal_splits) > 3:
        rho_split_z, p_split_z = spearmanr(z_splits, optimal_splits)
        print(f"\n  Optimal mass split vs z: ρ = {rho_split_z:+.3f}, p = {p_split_z:.4f}")
        if abs(rho_split_z) > 0.3:
            print(f"  → Mass threshold {'increases' if rho_split_z > 0 else 'decreases'} with z")
    
    # ============================================================
    # SCORECARD
    # ============================================================
    print("\n" + "=" * 70)
    print("SCORECARD")
    print("=" * 70)
    
    tests = []
    
    # T1: Step evolves
    tests.append(('Mass step evolves with z',
                  abs(rho_step_z) > 0.3 and p_step_z < 0.15,
                  f'ρ={rho_step_z:+.3f}, p={p_step_z:.4f}'))
    
    # T2: Color > stretch
    if len(color_steps) > 3:
        tests.append(('Color component evolves more than stretch',
                      abs(rho_color_z) > abs(rho_stretch_z),
                      f'|ρ_color|={abs(rho_color_z):.3f} vs |ρ_stretch|={abs(rho_stretch_z):.3f}'))
    
    # T3: β-step correlation
    if len(betas_binned) > 3:
        tests.append(('β(z) correlates with mass step',
                      abs(rho_beta_step) > 0.3,
                      f'ρ={rho_beta_step:+.3f}, p={p_beta_step:.4f}'))
    
    # T5: β(z) reduces step at high-z
    if len(reductions) > 0 and above_z0_mask.any():
        tests.append(('β(z) correction reduces step more at high-z',
                      avg_red_high > avg_red_low,
                      f'reduction: {avg_red_low:.0f}% (low-z) vs {avg_red_high:.0f}% (high-z)'))
    
    # T6: Color distributions diverge
    if len(ks_p_vals) > 3:
        tests.append(('Color distributions diverge with z',
                      rho_ks > 0,
                      f'ρ={rho_ks:+.3f}'))
    
    n_pass = sum(1 for _, p, _ in tests if p)
    n_total = len(tests)
    
    print()
    for name, passed, detail in tests:
        print(f"  {'✓' if passed else '✗'} {name}")
        print(f"    {detail}")
    
    print(f"\n  SCORE: {n_pass}/{n_total}")
    
    # ============================================================
    # SAVE
    # ============================================================
    results_dir = os.path.join(os.path.dirname(__file__), 'results_mass_step_evolution')
    os.makedirs(results_dir, exist_ok=True)
    
    results = {
        'global_mass_step': float(global_step),
        'step_z_rho': float(rho_step_z),
        'step_z_p': float(p_step_z),
        'beta_step_rho': float(rho_beta_step) if len(betas_binned) > 3 else None,
        'beta_step_p': float(p_beta_step) if len(betas_binned) > 3 else None,
        'color_z_rho': float(rho_color_z) if len(color_steps) > 3 else None,
        'stretch_z_rho': float(rho_stretch_z) if len(stretch_steps) > 3 else None,
        'tests_passed': n_pass,
        'tests_total': n_total,
        'steps_by_z': [{'z': float(z), 'step': float(s), 'err': float(e)} 
                       for z, s, e in zip(z_centers, steps, step_errs)],
        'beta_step_pairs': [{'z': float(z), 'beta': float(b), 'step': float(s)}
                           for z, b, s in zip(z_matched, betas_binned, steps_matched)] if len(betas_binned) > 0 else []
    }
    
    with open(os.path.join(results_dir, 'results.json'), 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to results_mass_step_evolution/results.json")

if __name__ == '__main__':
    main()
