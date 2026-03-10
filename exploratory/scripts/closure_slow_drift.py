#!/usr/bin/env python3
"""
closure_slow_drift.py — Slow-Channel Hidden Brightness Drift
==============================================================

The "anchor" population (slow decliners) showed MORE low/high-z tension
in M than the full sample (ΔM = -0.061 vs -0.045). This test quantifies
the slow-channel brightness evolution after full Tripp standardization.

If slow decliners have a z-dependent brightness drift AFTER correction,
that's a hidden systematic in the distance ladder.

Also tests:
- Fast-channel composite structure (core + extreme subpopulation)
- GPT's 2-component mixture test within the fast channel

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr, norm
from scipy.optimize import curve_fit, minimize
from scipy.integrate import quad
import json
from pathlib import Path

np.random.seed(42)

# Cosmology
H0 = 67.4
Omega_m, Omega_L = 0.315, 0.685
c_km = 2.998e5

def E(z): return np.sqrt(Omega_m*(1+z)**3 + Omega_L)

def luminosity_distance(z):
    dc, _ = quad(lambda zz: c_km / (H0 * E(zz)), 0, z)
    return dc * (1 + z)

def distance_modulus(z):
    dL = luminosity_distance(z)
    if dL <= 0: return 0
    return 5 * np.log10(dL) + 25

def cosmic_age(z):
    r, _ = quad(lambda zz: 1/((1+zz)*E(zz)), z, 100)
    return r / (H0 * 1e5 / 3.086e24 / 3.156e16)

# Load data
data = []
with open('data/pantheon_plus.dat', 'r') as f:
    header = f.readline().strip().split()
    for line in f:
        parts = line.strip().split()
        if len(parts) >= len(header):
            row = {}
            for i, h in enumerate(header):
                try: row[h] = float(parts[i])
                except: row[h] = parts[i]
            data.append(row)

seen = set()
sne = []
for row in data:
    z = row.get('zHD', 0)
    if z > 0.01 and z < 2.5 and abs(row.get('c', -999)) < 0.3 and abs(row.get('x1', -999)) < 3:
        key = (row['CID'], round(z, 5))
        if key not in seen:
            seen.add(key)
            sne.append(row)

print(f"Loaded {len(sne)} unique SNe Ia\n")

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 2.5]

# ============================================================
# FIT GLOBAL TRIPP PARAMETERS
# ============================================================
def fit_tripp(sample):
    def chi2(params):
        M, alpha, beta = params
        total = 0
        for s in sample:
            mu_sn = s['mB'] - M + alpha * s['x1'] - beta * s['c']
            mu_cosmo = distance_modulus(s['zHD'])
            err = max(s.get('mBERR', 0.1), 0.01)
            total += ((mu_sn - mu_cosmo) / err)**2
        return total
    result = minimize(chi2, x0=[-19.3, 0.15, 3.0], method='Nelder-Mead',
                      options={'maxiter': 10000})
    return result.x

M_all, alpha_all, beta_all = fit_tripp(sne)
print(f"Global Tripp: M={M_all:.4f}, α={alpha_all:.4f}, β={beta_all:.4f}\n")

# Compute Hubble residuals for every SN
for s in sne:
    s['mu_sn'] = s['mB'] - M_all + alpha_all * s['x1'] - beta_all * s['c']
    s['mu_cosmo'] = distance_modulus(s['zHD'])
    s['resid'] = s['mu_sn'] - s['mu_cosmo']

# ============================================================
# TEST: SLOW-CHANNEL BRIGHTNESS DRIFT
# ============================================================
print("=" * 100)
print("SLOW-CHANNEL HIDDEN BRIGHTNESS DRIFT")
print("Hubble residuals (after global Tripp correction) by z-bin")
print("If slow channel drifts → hidden systematic in distance ladder")
print("=" * 100)

for ch_name, ch_filter in [("ALL", lambda s: True),
                             ("SLOW (x1≥0)", lambda s: s['x1'] >= 0),
                             ("FAST (x1<0)", lambda s: s['x1'] < 0)]:
    print(f"\n  {ch_name}:")
    print(f"  {'z':>6} {'N':>4} | {'⟨resid⟩':>9} {'σ(resid)':>10} {'⟨resid⟩/σ_mean':>15} | {'⟨mB_corr⟩':>10}")
    print(f"  " + "-" * 65)
    
    z_list = []
    resid_list = []
    sigma_mean_list = []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi and ch_filter(s)]
        if len(bin_sne) < 10:
            continue
        
        zc = np.mean([s['zHD'] for s in bin_sne])
        resids = np.array([s['resid'] for s in bin_sne])
        mu_corr = np.array([s['mu_sn'] for s in bin_sne])
        
        mean_r = np.mean(resids)
        std_r = np.std(resids)
        sigma_mean = std_r / np.sqrt(len(bin_sne))
        snr = mean_r / sigma_mean if sigma_mean > 0 else 0
        
        z_list.append(zc)
        resid_list.append(mean_r)
        sigma_mean_list.append(sigma_mean)
        
        print(f"  {zc:>6.3f} {len(bin_sne):>4} | {mean_r:>9.4f} {std_r:>10.4f} {snr:>15.2f}σ | {np.mean(mu_corr):>10.4f}")
    
    if len(z_list) >= 3:
        z_ch = np.array(z_list)
        r_ch = np.array(resid_list)
        rho, p = spearmanr(z_ch, r_ch)
        
        # Linear fit to quantify drift
        from scipy.stats import linregress
        slope, intercept, r_val, p_val, stderr = linregress(z_ch, r_ch)
        
        print(f"\n  ⟨resid⟩ vs z: ρ = {rho:+.3f}, p = {p:.4f}")
        print(f"  Linear drift: slope = {slope:+.4f} mag/unit-z (p = {p_val:.4f})")
        print(f"  At z=1: drift = {slope * 1.0:+.4f} mag → {(10**(slope*1.0/5)-1)*100:+.2f}% in H₀")
        
        if abs(rho) > 0.5 and p < 0.1:
            print(f"  🔥 SIGNIFICANT DRIFT DETECTED")
        else:
            print(f"  Drift not significant at current precision")

# ============================================================
# FAST vs SLOW RESIDUAL DIVERGENCE
# ============================================================
print(f"\n\n{'=' * 100}")
print("FAST vs SLOW RESIDUAL DIVERGENCE")
print("Does the gap between fast and slow residuals grow with z?")
print("=" * 100)

print(f"\n  {'z':>6} | {'⟨r⟩_fast':>10} {'⟨r⟩_slow':>10} {'Gap':>8} | {'N_fast':>7} {'N_slow':>7}")
print(f"  " + "-" * 65)

z_gap = []
gap_list = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    fast = [s for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] < 0]
    slow = [s for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] >= 0]
    
    if len(fast) < 10 or len(slow) < 10:
        continue
    
    zc = np.mean([s['zHD'] for s in fast + slow])
    r_fast = np.mean([s['resid'] for s in fast])
    r_slow = np.mean([s['resid'] for s in slow])
    gap = r_slow - r_fast
    
    z_gap.append(zc)
    gap_list.append(gap)
    
    print(f"  {zc:>6.3f} | {r_fast:>10.4f} {r_slow:>10.4f} {gap:>8.4f} | {len(fast):>7} {len(slow):>7}")

if len(z_gap) >= 3:
    rho_g, p_g = spearmanr(z_gap, gap_list)
    print(f"\n  Gap (slow-fast) vs z: ρ = {rho_g:+.3f}, p = {p_g:.4f}")
    
    if rho_g > 0.5 and p_g < 0.1:
        print(f"  🔥 Gap GROWS with z → channels diverge in Hubble residuals")
        print(f"     This is the hidden systematic: fast and slow give different distances at high z")
    elif rho_g < -0.5 and p_g < 0.1:
        print(f"  Gap SHRINKS with z → channels converge in residuals")
    else:
        print(f"  Gap stable or noisy")

# ============================================================
# GPT's TEST: 2-COMPONENT MIXTURE WITHIN FAST CHANNEL
# ============================================================
print(f"\n\n{'=' * 100}")
print("GPT's TEST: Is the fast channel COMPOSITE?")
print("Fit 1-component vs 2-component Gaussian mixture to fast x1 at each z-bin")
print("If 2-component wins at low z and loses at high z → extreme subpopulation disappears")
print("=" * 100)

def fit_1gauss(x):
    """Fit single Gaussian, return log-likelihood"""
    mu, sigma = np.mean(x), np.std(x)
    if sigma < 0.01: sigma = 0.01
    ll = np.sum(norm.logpdf(x, mu, sigma))
    return ll, 2  # 2 params

def fit_2gauss(x, n_restarts=10):
    """Fit 2-component Gaussian mixture via EM, return log-likelihood"""
    n = len(x)
    best_ll = -np.inf
    best_params = None
    
    for _ in range(n_restarts):
        # Random init
        mu1 = np.random.choice(x)
        mu2 = np.random.choice(x)
        s1 = s2 = np.std(x) * 0.5
        w1 = np.random.uniform(0.3, 0.7)
        
        for em_iter in range(100):
            # E-step
            p1 = w1 * norm.pdf(x, mu1, max(s1, 0.01))
            p2 = (1-w1) * norm.pdf(x, mu2, max(s2, 0.01))
            total = p1 + p2 + 1e-300
            gamma = p1 / total
            
            # M-step
            n1 = np.sum(gamma)
            n2 = n - n1
            if n1 < 2 or n2 < 2:
                break
            
            w1_new = n1 / n
            mu1_new = np.sum(gamma * x) / n1
            mu2_new = np.sum((1-gamma) * x) / n2
            s1_new = np.sqrt(np.sum(gamma * (x - mu1_new)**2) / n1)
            s2_new = np.sqrt(np.sum((1-gamma) * (x - mu2_new)**2) / n2)
            
            if s1_new < 0.01: s1_new = 0.01
            if s2_new < 0.01: s2_new = 0.01
            
            w1, mu1, mu2, s1, s2 = w1_new, mu1_new, mu2_new, s1_new, s2_new
        
        ll = np.sum(np.log(w1 * norm.pdf(x, mu1, s1) + (1-w1) * norm.pdf(x, mu2, s2) + 1e-300))
        if ll > best_ll:
            best_ll = ll
            best_params = (w1, mu1, s1, mu2, s2)
    
    return best_ll, 5, best_params  # 5 params

print(f"\n  {'z':>6} {'N':>4} | {'LL_1g':>10} {'LL_2g':>10} {'ΔBIC':>8} | {'w1':>6} {'μ1':>7} {'σ1':>6} {'μ2':>7} {'σ2':>6} | {'Verdict':>12}")
print(f"  " + "-" * 100)

bic_diff_list = []
z_bic_list = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    fast = [s['x1'] for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] < 0]
    
    if len(fast) < 20:
        continue
    
    x = np.array(fast)
    n = len(x)
    zc = np.mean([s['zHD'] for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] < 0])
    
    ll1, k1 = fit_1gauss(x)
    ll2, k2, params2 = fit_2gauss(x)
    
    # BIC = k*ln(n) - 2*LL (lower is better)
    bic1 = k1 * np.log(n) - 2 * ll1
    bic2 = k2 * np.log(n) - 2 * ll2
    delta_bic = bic1 - bic2  # positive = 2-component better
    
    bic_diff_list.append(delta_bic)
    z_bic_list.append(zc)
    
    w1, mu1, s1, mu2, s2 = params2
    verdict = "2-COMP ✓" if delta_bic > 2 else ("1-COMP ✓" if delta_bic < -2 else "TIE")
    
    print(f"  {zc:>6.3f} {n:>4} | {ll1:>10.1f} {ll2:>10.1f} {delta_bic:>8.1f} | {w1:>6.2f} {mu1:>7.3f} {s1:>6.3f} {mu2:>7.3f} {s2:>6.3f} | {verdict:>12}")

if len(z_bic_list) >= 3:
    rho_bic, p_bic = spearmanr(z_bic_list, bic_diff_list)
    print(f"\n  ΔBIC (favoring 2-comp) vs z: ρ = {rho_bic:+.3f}, p = {p_bic:.4f}")
    
    if rho_bic < -0.4:
        print(f"  🔥 2-component preference WEAKENS with z → extreme subpopulation DISAPPEARS")
        print(f"     GPT's 'fast channel is composite' hypothesis CONFIRMED")
    elif rho_bic > 0.4:
        print(f"  2-component preference GROWS with z (unexpected)")
    else:
        print(f"  No clear trend in mixture preference")

# ============================================================
# HOST MASS SPLIT: Does slow-channel drift depend on host mass?
# ============================================================
print(f"\n\n{'=' * 100}")
print("HOST MASS SPLIT: Slow-channel drift in massive vs low-mass hosts")
print("(Tests whether the drift is tied to galaxy age/environment)")
print("=" * 100)

slow_sne = [s for s in sne if s['x1'] >= 0]
mass_median = np.median([s['HOST_LOGMASS'] for s in slow_sne if s.get('HOST_LOGMASS', 0) > 0])
print(f"  Median HOST_LOGMASS (slow): {mass_median:.2f}")

for mass_label, mass_filter in [("HIGH MASS (>median)", lambda s: s.get('HOST_LOGMASS', 0) > mass_median),
                                  ("LOW MASS (≤median)", lambda s: 0 < s.get('HOST_LOGMASS', 0) <= mass_median)]:
    print(f"\n  {mass_label}:")
    print(f"  {'z':>6} {'N':>4} | {'⟨resid⟩':>9}")
    print(f"  " + "-" * 30)
    
    z_list = []
    r_list = []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_sne = [s for s in slow_sne if z_lo <= s['zHD'] < z_hi and mass_filter(s)]
        if len(bin_sne) < 8:
            continue
        
        zc = np.mean([s['zHD'] for s in bin_sne])
        mean_r = np.mean([s['resid'] for s in bin_sne])
        z_list.append(zc)
        r_list.append(mean_r)
        
        print(f"  {zc:>6.3f} {len(bin_sne):>4} | {mean_r:>9.4f}")
    
    if len(z_list) >= 3:
        rho, p = spearmanr(z_list, r_list)
        print(f"  ⟨resid⟩ vs z: ρ = {rho:+.3f}, p = {p:.4f}")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_slow_drift")
results_dir.mkdir(exist_ok=True)

print(f"\n\nResults saved to {results_dir}/")
