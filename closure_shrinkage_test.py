#!/usr/bin/env python3
"""
closure_shrinkage_test.py — Is the contraction pipeline shrinkage?
===================================================================

GPT's objection: SALT2 fitter pulls estimates toward priors as
measurement errors grow with z. This artificially contracts the
distribution even if the TRUE population is unchanged.

Kill tests:
1. Does Var(x1) track the Wiener gain formula?
2. Does x1/σ_x1 collapse toward 0 at high z?
3. Does contraction vanish in high-SNR subsample?
4. Does Var(c) behave differently (flat = different prior/shrinkage)?

If shrinkage explains it → contraction is pipeline, not astrophysics.
If it DOESN'T explain it → hard cap is real.

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr
from scipy.optimize import minimize_scalar
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
    if z > 0.01 and z < 2.5 and abs(row.get('c',-999)) < 0.3 and abs(row.get('x1',-999)) < 3:
        key = (row['CID'], round(z, 5))
        if key not in seen:
            seen.add(key)
            sne.append(row)

print(f"Loaded {len(sne)} unique SNe Ia\n")

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.82, 1.0, 2.5]

# ============================================================
# TEST 1: Wiener Gain Formula for x1
# ============================================================
print("=" * 80)
print("TEST 1: Does Var(x1) follow Wiener shrinkage?")
print("=" * 80)

# Wiener gain: α(z) = σ²_prior / (σ²_prior + ⟨σ²_meas⟩_z)
# Var(x1_hat) = α²·Var(x1_true) + (1-α)²·⟨σ²_meas⟩
# Approximately: Var(x1_hat) ≈ α²·σ²_true + noise

print(f"\n{'z':>6} {'N':>4} {'Var(x1)':>9} {'⟨σ_x1⟩':>8} {'⟨σ²_x1⟩':>9} {'x1/σ_x1':>8} {'Var_intr':>9}")
print("-" * 65)

z_centers = []
var_x1_list = []
mean_sig_x1_list = []
mean_sig2_x1_list = []
snr_x1_list = []
var_intrinsic_x1 = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 15:
        continue
    
    x1 = np.array([s['x1'] for s in bin_sne])
    x1_err = np.array([s['x1ERR'] for s in bin_sne])
    
    var_x1 = np.var(x1)
    mean_sig = np.mean(x1_err)
    mean_sig2 = np.mean(x1_err**2)
    snr = np.mean(np.abs(x1) / x1_err)
    
    # Intrinsic variance = observed - mean measurement error
    var_intr = max(0, var_x1 - mean_sig2)
    
    z_center = np.mean([s['zHD'] for s in bin_sne])
    z_centers.append(z_center)
    var_x1_list.append(var_x1)
    mean_sig_x1_list.append(mean_sig)
    mean_sig2_x1_list.append(mean_sig2)
    snr_x1_list.append(snr)
    var_intrinsic_x1.append(var_intr)
    
    print(f"{z_center:>6.2f} {len(bin_sne):>4} {var_x1:>9.4f} {mean_sig:>8.4f} {mean_sig2:>9.5f} {snr:>8.2f} {var_intr:>9.4f}")

z_arr = np.array(z_centers)
vx = np.array(var_x1_list)
sig2 = np.array(mean_sig2_x1_list)
vi = np.array(var_intrinsic_x1)

# Does σ_x1 grow with z? (Required for shrinkage argument)
rho_sig, p_sig = spearmanr(z_arr, mean_sig_x1_list)
print(f"\n  ⟨σ_x1⟩ vs z: ρ = {rho_sig:+.3f}, p = {p_sig:.4f}")

if rho_sig > 0.5:
    print(f"  ✓ Measurement errors GROW with z (shrinkage prerequisite met)")
else:
    print(f"  ✗ Errors don't grow significantly — shrinkage less plausible")

# Does Var_intrinsic still contract after subtracting measurement error?
rho_vi, p_vi = spearmanr(z_arr, vi)
print(f"\n  Var_intrinsic(x1) vs z: ρ = {rho_vi:+.3f}, p = {p_vi:.4f}")

if rho_vi < -0.3:
    print(f"  🔥 INTRINSIC variance still contracts after error subtraction")
    print(f"     → NOT just measurement shrinkage — real population effect")
elif abs(rho_vi) < 0.3:
    print(f"  ⚠️ Intrinsic variance FLAT after error subtraction")
    print(f"     → Contraction may be largely measurement shrinkage")

# Fit Wiener model: Var(x1) = σ²_prior² / (σ²_prior + σ²_meas)² × σ²_true + noise
# Simplified: find σ_prior that best predicts observed Var(x1)
def wiener_prediction(sig_prior, sig2_meas, var_true):
    alpha = sig_prior**2 / (sig_prior**2 + sig2_meas)
    return alpha**2 * var_true + (1-alpha)**2 * sig2_meas

# Assume var_true = Var at lowest z (most reliable)
var_true_est = vi[0] if vi[0] > 0 else vx[0]

def wiener_cost(log_sig_prior):
    sig_prior = np.exp(log_sig_prior)
    pred = np.array([wiener_prediction(sig_prior, s2, var_true_est) for s2 in sig2])
    return np.sum((vx - pred)**2)

result = minimize_scalar(wiener_cost, bounds=(-2, 3), method='bounded')
sig_prior_best = np.exp(result.x)
pred_wiener = np.array([wiener_prediction(sig_prior_best, s2, var_true_est) for s2 in sig2])

ss_res = np.sum((vx - pred_wiener)**2)
ss_tot = np.sum((vx - np.mean(vx))**2)
r_sq_wiener = 1 - ss_res/ss_tot if ss_tot > 0 else 0

print(f"\n  Wiener model fit:")
print(f"  σ_prior = {sig_prior_best:.3f}")
print(f"  R² = {r_sq_wiener:.3f}")
print(f"  Var_true (from z=0 bin) = {var_true_est:.4f}")

if r_sq_wiener > 0.7:
    print(f"  ❄️ Wiener shrinkage EXPLAINS the contraction (R²={r_sq_wiener:.2f})")
else:
    print(f"  🔥 Wiener shrinkage does NOT fully explain it (R²={r_sq_wiener:.2f})")

# ============================================================
# TEST 2: x1/σ_x1 SNR distribution
# ============================================================
print(f"\n\n" + "=" * 80)
print("TEST 2: x1/σ_x1 (SNR) collapse toward zero?")
print("=" * 80)

rho_snr, p_snr = spearmanr(z_arr, snr_x1_list)
print(f"\n  ⟨|x1/σ_x1|⟩ vs z: ρ = {rho_snr:+.3f}, p = {p_snr:.4f}")

if rho_snr < -0.5:
    print(f"  ✓ SNR collapses with z → consistent with shrinkage")
else:
    print(f"  ✗ SNR doesn't collapse → shrinkage not dominant")

# ============================================================
# TEST 3: High-SNR subsample — does contraction vanish?
# ============================================================
print(f"\n\n" + "=" * 80)
print("TEST 3: High-SNR Subsample — Contraction without shrinkage?")
print("=" * 80)

# Cut on σ_x1 < median (keep only well-measured SNe)
all_x1err = np.array([s['x1ERR'] for s in sne])
err_cuts = [np.percentile(all_x1err, 25), np.percentile(all_x1err, 50)]

for pct, err_cut in zip([25, 50], err_cuts):
    highsnr = [s for s in sne if s['x1ERR'] < err_cut]
    
    z_hs = []
    vx_hs = []
    vc_hs = []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_hs = [s for s in highsnr if z_lo <= s['zHD'] < z_hi]
        if len(bin_hs) < 10:
            continue
        
        x1 = np.array([s['x1'] for s in bin_hs])
        c = np.array([s['c'] for s in bin_hs])
        
        z_hs.append(np.mean([s['zHD'] for s in bin_hs]))
        vx_hs.append(np.var(x1))
        vc_hs.append(np.var(c))
    
    if len(z_hs) >= 4:
        rho_vx_hs, p_vx_hs = spearmanr(z_hs, vx_hs)
        rho_vc_hs, p_vc_hs = spearmanr(z_hs, vc_hs)
        print(f"\n  σ_x1 < P{pct} ({err_cut:.3f}): N={len(highsnr)}")
        print(f"    Var(x1) vs z: ρ = {rho_vx_hs:+.3f}, p = {p_vx_hs:.4f}")
        print(f"    Var(c) vs z:  ρ = {rho_vc_hs:+.3f}, p = {p_vc_hs:.4f}")
        
        if rho_vx_hs < -0.3:
            print(f"    🔥 Contraction PERSISTS in high-SNR subsample")
            print(f"       → NOT just shrinkage — real population effect")
        elif abs(rho_vx_hs) < 0.3:
            print(f"    ⚠️ Contraction vanishes in high-SNR subsample")
            print(f"       → Shrinkage may explain it")

# ============================================================
# TEST 4: Same analysis for color
# ============================================================
print(f"\n\n" + "=" * 80)
print("TEST 4: Color Shrinkage Check")
print("=" * 80)

print(f"\n{'z':>6} {'N':>4} {'Var(c)':>9} {'⟨σ_c⟩':>8} {'⟨σ²_c⟩':>9} {'c/σ_c':>8} {'Var_intr':>9}")
print("-" * 65)

var_c_list = []
var_intr_c_list = []
snr_c_list = []

for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 15:
        continue
    
    c = np.array([s['c'] for s in bin_sne])
    c_err = np.array([s['cERR'] for s in bin_sne])
    
    var_c = np.var(c)
    mean_sig_c = np.mean(c_err)
    mean_sig2_c = np.mean(c_err**2)
    snr_c = np.mean(np.abs(c) / c_err)
    var_intr_c = max(0, var_c - mean_sig2_c)
    
    z_center = np.mean([s['zHD'] for s in bin_sne])
    var_c_list.append(var_c)
    var_intr_c_list.append(var_intr_c)
    snr_c_list.append(snr_c)
    
    print(f"{z_center:>6.2f} {len(bin_sne):>4} {var_c:>9.6f} {mean_sig_c:>8.4f} {mean_sig2_c:>9.6f} {snr_c:>8.2f} {var_intr_c:>9.6f}")

vic = np.array(var_intr_c_list)
rho_vic, p_vic = spearmanr(z_arr[:len(vic)], vic)
print(f"\n  Var_intrinsic(c) vs z: ρ = {rho_vic:+.3f}, p = {p_vic:.4f}")

if abs(rho_vic) < 0.3:
    print(f"  Color intrinsic variance is FLAT → no shrinkage in color")
    print(f"  This matches: x1 shrinks (high σ_x1), c doesn't (low σ_c)")

# ============================================================
# SAVE & VERDICT
# ============================================================
results_dir = Path("results_shrinkage")
results_dir.mkdir(exist_ok=True)

results = {
    'wiener_r_sq': float(r_sq_wiener),
    'sig_prior': float(sig_prior_best),
    'sigma_x1_vs_z_rho': float(rho_sig),
    'intrinsic_x1_vs_z_rho': float(rho_vi),
    'snr_x1_vs_z_rho': float(rho_snr),
    'intrinsic_c_vs_z_rho': float(rho_vic),
}

with open(results_dir / "shrinkage_test.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {results_dir}/")
print("\n" + "=" * 80)
print("FINAL VERDICT: Is the contraction pipeline shrinkage?")
print("=" * 80)
print(f"  σ_x1 grows with z:           ρ = {rho_sig:+.3f} ({'YES' if rho_sig > 0.3 else 'NO'})")
print(f"  Wiener model fit:             R² = {r_sq_wiener:.3f} ({'EXPLAINS' if r_sq_wiener > 0.7 else 'PARTIAL' if r_sq_wiener > 0.4 else 'NO'})")
print(f"  Intrinsic Var(x1) contracts:  ρ = {rho_vi:+.3f} ({'REAL' if rho_vi < -0.3 else 'FLAT'})")
print(f"  SNR collapses:                ρ = {rho_snr:+.3f} ({'YES' if rho_snr < -0.3 else 'NO'})")
print(f"  Color intrinsic Var:          ρ = {rho_vic:+.3f} (FLAT = different behavior)")
