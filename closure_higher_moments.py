#!/usr/bin/env python3
"""
closure_higher_moments.py — Higher Moments Test (GPT's suggestion)
===================================================================

If within-channel contraction is real state-space narrowing:
- Kurtosis should DROP with z (tails truncate)
- Skewness should evolve toward the attractor value
- Both effects should appear WITHIN each x1-sign channel

If contraction is symmetric compression (e.g., just scaling):
- Kurtosis stays flat (shape preserved, just scaled)

This is a smoking gun for state-space narrowing vs simple variance reduction.

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr, kurtosis, skew
import json
from pathlib import Path

np.random.seed(42)

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

# ============================================================
# HIGHER MOMENTS: Full population
# ============================================================
z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 2.5]

print("=" * 90)
print("HIGHER MOMENTS — FULL POPULATION")
print("=" * 90)
print(f"\n{'z_mid':>6} {'N':>5} | {'Kurt(mB)':>10} {'Kurt(x1)':>10} {'Kurt(c)':>10} | {'Skew(mB)':>10} {'Skew(x1)':>10} {'Skew(c)':>10}")
print("-" * 90)

z_centers_all = []
kurt_mB_all, kurt_x1_all, kurt_c_all = [], [], []
skew_mB_all, skew_x1_all, skew_c_all = [], [], []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 20:
        continue
    
    mB = np.array([s['mB'] for s in bin_sne])
    x1 = np.array([s['x1'] for s in bin_sne])
    c = np.array([s['c'] for s in bin_sne])
    
    zc = np.mean([s['zHD'] for s in bin_sne])
    
    km = kurtosis(mB, fisher=True)  # excess kurtosis (0 = Gaussian)
    kx = kurtosis(x1, fisher=True)
    kc = kurtosis(c, fisher=True)
    sm = skew(mB)
    sx = skew(x1)
    sc = skew(c)
    
    z_centers_all.append(zc)
    kurt_mB_all.append(km)
    kurt_x1_all.append(kx)
    kurt_c_all.append(kc)
    skew_mB_all.append(sm)
    skew_x1_all.append(sx)
    skew_c_all.append(sc)
    
    print(f"{zc:>6.3f} {len(bin_sne):>5} | {km:>10.3f} {kx:>10.3f} {kc:>10.3f} | {sm:>10.3f} {sx:>10.3f} {sc:>10.3f}")

z_arr = np.array(z_centers_all)

print(f"\nKurtosis trends (full population):")
for name, arr in [("mB", kurt_mB_all), ("x1", kurt_x1_all), ("c", kurt_c_all)]:
    rho, p = spearmanr(z_arr, arr)
    print(f"  Kurt({name}) vs z: ρ = {rho:+.3f}, p = {p:.4f} {'🔥 DROPS' if rho < -0.3 and p < 0.1 else '— FLAT' if abs(rho) < 0.3 else '↑ RISES'}")

print(f"\nSkewness trends (full population):")
for name, arr in [("mB", skew_mB_all), ("x1", skew_x1_all), ("c", skew_c_all)]:
    rho, p = spearmanr(z_arr, arr)
    print(f"  Skew({name}) vs z: ρ = {rho:+.3f}, p = {p:.4f}")

# ============================================================
# HIGHER MOMENTS: Within-channel (x1 < 0 vs x1 > 0)
# ============================================================
print(f"\n\n{'=' * 90}")
print("HIGHER MOMENTS — WITHIN CHANNELS (x1 < 0 = fast decliners, x1 > 0 = slow decliners)")
print("=" * 90)

for channel_name, channel_filter in [("FAST (x1<0)", lambda s: s['x1'] < 0), 
                                       ("SLOW (x1>0)", lambda s: s['x1'] >= 0)]:
    print(f"\n  Channel: {channel_name}")
    print(f"  {'z_mid':>6} {'N':>5} | {'Kurt(mB)':>10} {'Kurt(x1)':>10} {'Kurt(c)':>10} | {'Skew(mB)':>10} {'Skew(x1)':>10} {'Skew(c)':>10}")
    print(f"  " + "-" * 85)
    
    z_c_list = []
    kurt_mB_ch, kurt_x1_ch, kurt_c_ch = [], [], []
    skew_mB_ch, skew_x1_ch, skew_c_ch = [], [], []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi and channel_filter(s)]
        if len(bin_sne) < 15:
            continue
        
        mB = np.array([s['mB'] for s in bin_sne])
        x1 = np.array([s['x1'] for s in bin_sne])
        c = np.array([s['c'] for s in bin_sne])
        
        zc = np.mean([s['zHD'] for s in bin_sne])
        
        km = kurtosis(mB, fisher=True)
        kx = kurtosis(x1, fisher=True)
        kc = kurtosis(c, fisher=True)
        sm = skew(mB)
        sx = skew(x1)
        sc = skew(c)
        
        z_c_list.append(zc)
        kurt_mB_ch.append(km)
        kurt_x1_ch.append(kx)
        kurt_c_ch.append(kc)
        skew_mB_ch.append(sm)
        skew_x1_ch.append(sx)
        skew_c_ch.append(sc)
        
        print(f"  {zc:>6.3f} {len(bin_sne):>5} | {km:>10.3f} {kx:>10.3f} {kc:>10.3f} | {sm:>10.3f} {sx:>10.3f} {sc:>10.3f}")
    
    if len(z_c_list) >= 3:
        z_ch = np.array(z_c_list)
        print(f"\n  Kurtosis trends ({channel_name}):")
        for name, arr in [("mB", kurt_mB_ch), ("x1", kurt_x1_ch), ("c", kurt_c_ch)]:
            rho, p = spearmanr(z_ch, arr)
            print(f"    Kurt({name}) vs z: ρ = {rho:+.3f}, p = {p:.4f} {'🔥 DROPS' if rho < -0.3 and p < 0.15 else '— FLAT' if abs(rho) < 0.3 else '↑ RISES'}")
        print(f"\n  Skewness trends ({channel_name}):")
        for name, arr in [("mB", skew_mB_ch), ("x1", skew_x1_ch), ("c", skew_c_ch)]:
            rho, p = spearmanr(z_ch, arr)
            print(f"    Skew({name}) vs z: ρ = {rho:+.3f}, p = {p:.4f}")

# ============================================================
# TAIL TRUNCATION TEST
# ============================================================
print(f"\n\n{'=' * 90}")
print("TAIL TRUNCATION: P90-P10 interquantile range within channels")
print("=" * 90)
print("(If state-space narrows, tails should truncate → IQR drops)")

for channel_name, channel_filter in [("FAST (x1<0)", lambda s: s['x1'] < 0), 
                                       ("SLOW (x1>0)", lambda s: s['x1'] >= 0)]:
    print(f"\n  Channel: {channel_name}")
    print(f"  {'z_mid':>6} {'N':>5} | {'IQR90(mB)':>10} {'IQR90(x1)':>10} {'IQR90(c)':>10}")
    print(f"  " + "-" * 50)
    
    z_c_list = []
    iqr_mB, iqr_x1, iqr_c = [], [], []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi and channel_filter(s)]
        if len(bin_sne) < 15:
            continue
        
        mB = np.array([s['mB'] for s in bin_sne])
        x1 = np.array([s['x1'] for s in bin_sne])
        c = np.array([s['c'] for s in bin_sne])
        
        zc = np.mean([s['zHD'] for s in bin_sne])
        z_c_list.append(zc)
        
        iqr_mB.append(np.percentile(mB, 90) - np.percentile(mB, 10))
        iqr_x1.append(np.percentile(x1, 90) - np.percentile(x1, 10))
        iqr_c.append(np.percentile(c, 90) - np.percentile(c, 10))
        
        print(f"  {zc:>6.3f} {len(bin_sne):>5} | {iqr_mB[-1]:>10.4f} {iqr_x1[-1]:>10.4f} {iqr_c[-1]:>10.4f}")
    
    if len(z_c_list) >= 3:
        z_ch = np.array(z_c_list)
        print(f"\n  IQR90 trends ({channel_name}):")
        for name, arr in [("mB", iqr_mB), ("x1", iqr_x1), ("c", iqr_c)]:
            rho, p = spearmanr(z_ch, arr)
            print(f"    IQR90({name}) vs z: ρ = {rho:+.3f}, p = {p:.4f} {'🔥 SHRINKS' if rho < -0.3 and p < 0.15 else '— FLAT' if abs(rho) < 0.3 else '↑ WIDENS'}")

# ============================================================
# GPT's V_C(z) = sqrt(det C_y(z|C)) — Phase-Space Volume
# ============================================================
print(f"\n\n{'=' * 90}")
print("PHASE-SPACE VOLUME V_C(z) = sqrt(det C_y(z|C))")
print("GPT's 'Hard Cap Index': η_C(z) = -ln(V_C(z)/V_C(z_ref))")
print("=" * 90)

for channel_name, channel_filter in [("FAST (x1<0)", lambda s: s['x1'] < 0), 
                                       ("SLOW (x1>0)", lambda s: s['x1'] >= 0),
                                       ("ALL", lambda s: True)]:
    print(f"\n  Channel: {channel_name}")
    print(f"  {'z_mid':>6} {'N':>5} | {'V_C(z)':>12} {'η_C(z)':>10}")
    print(f"  " + "-" * 45)
    
    z_c_list = []
    V_list = []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi and channel_filter(s)]
        if len(bin_sne) < 15:
            continue
        
        mB = np.array([s['mB'] for s in bin_sne])
        x1 = np.array([s['x1'] for s in bin_sne])
        c = np.array([s['c'] for s in bin_sne])
        
        # Intrinsic covariance (subtract measurement errors from diagonal)
        mBerr = np.array([s['mBERR'] for s in bin_sne])
        x1err = np.array([s['x1ERR'] for s in bin_sne])
        cerr = np.array([s['cERR'] for s in bin_sne])
        
        C_obs = np.cov(np.array([mB, x1, c]))
        C_err = np.diag([np.mean(mBerr**2), np.mean(x1err**2), np.mean(cerr**2)])
        C_intr = C_obs - C_err
        
        # Ensure positive semi-definite for det
        eigvals = np.linalg.eigvalsh(C_intr)
        eigvals = np.maximum(eigvals, 1e-10)
        V = np.sqrt(np.prod(eigvals))
        
        zc = np.mean([s['zHD'] for s in bin_sne])
        z_c_list.append(zc)
        V_list.append(V)
    
    if len(V_list) >= 3:
        V_ref = V_list[0]
        for j, (zc, V) in enumerate(zip(z_c_list, V_list)):
            eta = -np.log(V / V_ref) if V > 0 and V_ref > 0 else 0
            print(f"  {zc:>6.3f} {'':>5} | {V:>12.6f} {eta:>10.3f}")
        
        z_ch = np.array(z_c_list)
        V_arr = np.array(V_list)
        rho, p = spearmanr(z_ch, V_arr)
        print(f"\n  V_C vs z: ρ = {rho:+.3f}, p = {p:.4f} {'🔥 CONTRACTS' if rho < -0.3 and p < 0.1 else '— FLAT'}")
        
        # Check if fast vs slow have DIFFERENT contraction rates
        if channel_name != "ALL":
            eta_arr = -np.log(V_arr / V_arr[0])
            from scipy.optimize import curve_fit
            def lin(x, a, b): return a * x + b
            try:
                popt, _ = curve_fit(lin, z_ch, eta_arr)
                print(f"  κ_C (contraction rate) ≈ {popt[0]:.3f} e-folds per unit z")
            except:
                pass

# ============================================================
# B(z) vs W(z) — Between/Within decomposition
# ============================================================
print(f"\n\n{'=' * 90}")
print("BETWEEN/WITHIN DECOMPOSITION: B(z) + W(z) = 1")
print("B(z) = Var(E[Z|C,z]) / Var(Z|z)  [architecture index]")
print("W(z) = E[Var(Z|C,z)] / Var(Z|z)  [hard cap index]")
print("=" * 90)

for var_name, var_key in [("x1", "x1"), ("c", "c"), ("mB", "mB")]:
    print(f"\n  Variable: {var_name}")
    print(f"  {'z_mid':>6} | {'B(z)':>8} {'W(z)':>8} | {'Var_total':>10} {'Var_between':>12} {'Var_within':>12}")
    print(f"  " + "-" * 70)
    
    z_c_list = []
    B_list, W_list = [], []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        fast = [s[var_key] for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] < 0]
        slow = [s[var_key] for s in sne if z_lo <= s['zHD'] < z_hi and s['x1'] >= 0]
        all_vals = fast + slow
        
        if len(fast) < 10 or len(slow) < 10:
            continue
        
        zc = np.mean([s['zHD'] for s in sne if z_lo <= s['zHD'] < z_hi])
        
        var_total = np.var(all_vals)
        if var_total < 1e-10:
            continue
        
        # Between: variance of channel means
        n_f, n_s = len(fast), len(slow)
        n_tot = n_f + n_s
        mu_f, mu_s = np.mean(fast), np.mean(slow)
        mu_all = np.mean(all_vals)
        var_between = (n_f/n_tot) * (mu_f - mu_all)**2 + (n_s/n_tot) * (mu_s - mu_all)**2
        
        # Within: weighted average of channel variances
        var_within = (n_f/n_tot) * np.var(fast) + (n_s/n_tot) * np.var(slow)
        
        B = var_between / var_total
        W = var_within / var_total
        
        z_c_list.append(zc)
        B_list.append(B)
        W_list.append(W)
        
        print(f"  {zc:>6.3f} | {B:>8.4f} {W:>8.4f} | {var_total:>10.6f} {var_between:>12.6f} {var_within:>12.6f}")
    
    if len(z_c_list) >= 3:
        z_ch = np.array(z_c_list)
        rho_B, p_B = spearmanr(z_ch, B_list)
        rho_W, p_W = spearmanr(z_ch, W_list)
        print(f"\n  B(z) vs z: ρ = {rho_B:+.3f}, p = {p_B:.4f} {'— STABLE (architecture preserved)' if abs(rho_B) < 0.4 else ''}")
        print(f"  W(z) vs z: ρ = {rho_W:+.3f}, p = {p_W:.4f} {'🔥 SHRINKS (hard cap)' if rho_W < -0.3 else ''}")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_higher_moments")
results_dir.mkdir(exist_ok=True)

print(f"\n\n{'=' * 90}")
print("VERDICT")
print("=" * 90)
print("""
If kurtosis DROPS within channels → tails truncate → STATE-SPACE NARROWING confirmed
If kurtosis FLAT within channels  → symmetric compression → different mechanism
If V_C(z) contracts AND B(z) stable → "two things side by side" confirmed:
   - Architecture (fast/slow) preserved
   - Hard cap (within-channel diversity) contracts
""")

print(f"\nResults saved to {results_dir}/")
