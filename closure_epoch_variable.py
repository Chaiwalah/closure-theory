#!/usr/bin/env python3
"""
closure_epoch_variable.py — Hidden Epoch Variable Test (GPT's key insight)
===========================================================================

GPT's sharpest move: the generator isn't z itself.
The true control variable η(z) is the shrinking progenitor phase-space volume.

Test: replace z with physically motivated proxies:
  1. HOST_LOGMASS (metallicity proxy — mass-metallicity relation)
  2. Prompt fraction proxy: f_prompt(z) based on DTD ~ t^{-1}
  3. Cosmic age fraction: t(z) / t_0

If any η(z) proxy gives R² > the raw-z semigroup fit (R² ≈ 0.86-0.90),
we've found the generator. That's Level 10.

Also tests GPT's "two-operator model":
  V_C(z) = V_C,0 · exp(-κ_epoch · η(z)) · exp(-κ_path · Σ(z))

Author: Closure Theory Collaboration  
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr
from scipy.optimize import curve_fit
from scipy.integrate import quad
import json
from pathlib import Path

np.random.seed(42)

# Cosmology
H0 = 67.4
Omega_m, Omega_L = 0.315, 0.685
Omega_b = 0.0493
m_p = 1.6726e-24
rho_crit = 1.878e-29
n_H0 = Omega_b * rho_crit * (H0/100)**2 / m_p * 0.76
H0_cgs = H0 * 1e5 / 3.086e24
c_cgs = 2.998e10

def E(z): return np.sqrt(Omega_m*(1+z)**3 + Omega_L)

def Sigma(z):
    if z <= 0: return 0
    r, _ = quad(lambda zz: n_H0*(1+zz)**2*c_cgs/(H0_cgs*E(zz)), 0, z)
    return r

def lookback_time(z):
    r, _ = quad(lambda zz: 1/((1+zz)*E(zz)), 0, z)
    return r / (H0 * 1e5 / 3.086e24 / 3.156e16)

def cosmic_age(z):
    """Age of universe at redshift z in Gyr"""
    age_now = lookback_time(100)  # approximate
    # Better: integrate from z to infinity
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

print(f"Loaded {len(sne)} unique SNe Ia")

# Check host mass availability
has_mass = [s for s in sne if s.get('HOST_LOGMASS', -999) > 0]
print(f"With HOST_LOGMASS: {len(has_mass)} ({100*len(has_mass)/len(sne):.0f}%)\n")

z_edges = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.82, 1.0, 2.5]

# ============================================================
# Compute bin-level statistics
# ============================================================

def compute_intrinsic_vars(bin_sne):
    """Return intrinsic variances and phase-space volume for a bin"""
    mB = np.array([s['mB'] for s in bin_sne])
    x1 = np.array([s['x1'] for s in bin_sne])
    c = np.array([s['c'] for s in bin_sne])
    mBerr = np.array([s['mBERR'] for s in bin_sne])
    x1err = np.array([s['x1ERR'] for s in bin_sne])
    cerr = np.array([s['cERR'] for s in bin_sne])
    
    vi_mB = max(1e-6, np.var(mB) - np.mean(mBerr**2))
    vi_x1 = max(1e-6, np.var(x1) - np.mean(x1err**2))
    vi_c = max(1e-6, np.var(c) - np.mean(cerr**2))
    
    # Phase-space volume
    C_obs = np.cov(np.array([mB, x1, c]))
    C_err = np.diag([np.mean(mBerr**2), np.mean(x1err**2), np.mean(cerr**2)])
    C_intr = C_obs - C_err
    eigvals = np.linalg.eigvalsh(C_intr)
    eigvals = np.maximum(eigvals, 1e-10)
    V = np.sqrt(np.prod(eigvals))
    
    return vi_mB, vi_x1, vi_c, V

# Build bins
bins = []
for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi]
    if len(bin_sne) < 15:
        continue
    
    zc = np.mean([s['zHD'] for s in bin_sne])
    vi_mB, vi_x1, vi_c, V = compute_intrinsic_vars(bin_sne)
    
    # Host mass (mean of bin, metallicity proxy)
    masses = [s['HOST_LOGMASS'] for s in bin_sne if s.get('HOST_LOGMASS', -999) > 0]
    mean_mass = np.mean(masses) if len(masses) > 5 else np.nan
    
    bins.append({
        'z': zc,
        'sigma': Sigma(zc),
        't_lookback': lookback_time(zc),
        'age': cosmic_age(zc),
        'host_mass': mean_mass,
        'N': len(bin_sne),
        'vi_mB': vi_mB, 'vi_x1': vi_x1, 'vi_c': vi_c,
        'V': V
    })

# ============================================================
# Define candidate η(z) variables
# ============================================================

# 1. Raw z
eta_z = np.array([b['z'] for b in bins])

# 2. Σ(z) — column density
eta_sigma = np.array([b['sigma'] for b in bins])

# 3. Lookback time
eta_tlb = np.array([b['t_lookback'] for b in bins])

# 4. Cosmic age at z (younger universe = less diverse)
eta_age = np.array([b['age'] for b in bins])

# 5. log(age) — logarithmic compression
eta_log_age = np.log(np.maximum(eta_age, 0.1))

# 6. Mean host mass (metallicity proxy via mass-metallicity relation)
eta_mass = np.array([b['host_mass'] for b in bins])
mass_valid = ~np.isnan(eta_mass)

# 7. Prompt fraction proxy: f_prompt ~ 1 / (1 + age/τ_delay)
tau_delay = 3.0  # Gyr, characteristic delay time
eta_prompt = 1.0 / (1.0 + eta_age / tau_delay)

# 8. DTD-truncated diversity: integral of t^{-1} from t_min to age(z)
# Represents available delay-time phase space
t_min = 0.04  # Gyr, minimum delay time
eta_dtd = np.log(np.maximum(eta_age, t_min) / t_min)

# 9. Metallicity-evolution proxy: Z(z) ~ Z_solar * 10^{-0.15z} (rough)
eta_metal = 10**(-0.15 * eta_z)

# ============================================================
# FIT: V(η) = V_0 exp(-κ η)  for each candidate
# ============================================================

print("=" * 100)
print("HIDDEN EPOCH VARIABLE: Which η(z) best predicts V_C(z)?")
print("=" * 100)

V_arr = np.array([b['V'] for b in bins])
ln_V = np.log(V_arr)

# For each intrinsic variance too
vi_mB_arr = np.array([b['vi_mB'] for b in bins])
vi_x1_arr = np.array([b['vi_x1'] for b in bins])
vi_c_arr = np.array([b['vi_c'] for b in bins])

def exp_decay(eta, V0, kappa):
    return V0 * np.exp(-kappa * eta)

def fit_and_score(eta, y, name, eta_name):
    """Fit exp decay and return R²"""
    try:
        # Normalize eta to [0, 1] for numerical stability
        eta_n = (eta - eta.min()) / (eta.max() - eta.min() + 1e-10)
        popt, _ = curve_fit(exp_decay, eta_n, y, p0=[y[0], 1.0], 
                           bounds=([0, -20], [100, 50]), maxfev=10000)
        pred = exp_decay(eta_n, *popt)
        ss_res = np.sum((y - pred)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return r_sq, popt[1]
    except:
        return -999, 0

candidates = [
    ("z (raw redshift)", eta_z, np.ones(len(bins), dtype=bool)),
    ("Σ(z) (column density)", eta_sigma, np.ones(len(bins), dtype=bool)),
    ("t_lookback (Gyr)", eta_tlb, np.ones(len(bins), dtype=bool)),
    ("Age(z) (cosmic age)", eta_age, np.ones(len(bins), dtype=bool)),
    ("log(Age)", eta_log_age, np.ones(len(bins), dtype=bool)),
    ("HOST_LOGMASS (Z proxy)", eta_mass, mass_valid),
    ("f_prompt (DTD proxy)", eta_prompt, np.ones(len(bins), dtype=bool)),
    ("DTD phase-space (log)", eta_dtd, np.ones(len(bins), dtype=bool)),
    ("Z_proxy (10^{-0.15z})", eta_metal, np.ones(len(bins), dtype=bool)),
]

print(f"\n{'Candidate η(z)':<30} | {'R²(V)':>8} {'R²(mB)':>8} {'R²(x1)':>8} {'R²(c)':>8} | {'κ_V':>8}")
print("-" * 95)

best_r2 = -999
best_name = ""

for name, eta, mask in candidates:
    eta_m = eta[mask]
    if len(eta_m) < 4 or np.std(eta_m) < 1e-10:
        print(f"{name:<30} | {'N/A':>8}")
        continue
    
    r2_V, kappa_V = fit_and_score(eta_m, V_arr[mask], "V", name)
    r2_mB, _ = fit_and_score(eta_m, vi_mB_arr[mask], "mB", name)
    r2_x1, _ = fit_and_score(eta_m, vi_x1_arr[mask], "x1", name)
    r2_c, _ = fit_and_score(eta_m, vi_c_arr[mask], "c", name)
    
    avg_r2 = np.mean([r2_V, r2_mB, r2_x1, r2_c])
    marker = " 🔥 BEST" if avg_r2 > best_r2 and r2_V > 0 else ""
    if avg_r2 > best_r2 and r2_V > 0:
        best_r2 = avg_r2
        best_name = name
    
    print(f"{name:<30} | {r2_V:>8.3f} {r2_mB:>8.3f} {r2_x1:>8.3f} {r2_c:>8.3f} | {kappa_V:>8.3f}{marker}")

print(f"\n  → Best hidden variable: {best_name} (avg R² = {best_r2:.3f})")

# ============================================================
# WITHIN-CHANNEL: Compare κ_fast vs κ_slow
# ============================================================
print(f"\n\n{'=' * 100}")
print("WITHIN-CHANNEL CONTRACTION RATES: κ_fast vs κ_slow")
print("(If they differ → tells us WHICH latent variable is shrinking)")
print("=" * 100)

for channel_name, ch_filter in [("FAST (x1<0)", lambda s: s['x1'] < 0),
                                  ("SLOW (x1>0)", lambda s: s['x1'] >= 0)]:
    ch_bins = []
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bin_sne = [s for s in sne if z_lo <= s['zHD'] < z_hi and ch_filter(s)]
        if len(bin_sne) < 10:
            continue
        zc = np.mean([s['zHD'] for s in bin_sne])
        _, _, _, V = compute_intrinsic_vars(bin_sne)
        ch_bins.append({'z': zc, 'V': V, 'N': len(bin_sne)})
    
    if len(ch_bins) < 3:
        print(f"\n  {channel_name}: insufficient bins")
        continue
    
    z_ch = np.array([b['z'] for b in ch_bins])
    V_ch = np.array([b['V'] for b in ch_bins])
    
    r2, kappa = fit_and_score(z_ch, V_ch, "V", "z")
    rho, p = spearmanr(z_ch, V_ch)
    
    print(f"\n  {channel_name}:")
    print(f"    V_C vs z: ρ = {rho:+.3f}, p = {p:.4f}")
    print(f"    Semigroup fit: R² = {r2:.3f}, κ = {kappa:.3f}")
    
    # η from V
    eta_ch = -np.log(V_ch / V_ch[0])
    print(f"    η range: {eta_ch[0]:.3f} → {eta_ch[-1]:.3f} ({eta_ch[-1] - eta_ch[0]:.2f} e-folds)")

# ============================================================
# TWO-OPERATOR MODEL: V = V_0 exp(-κ_epoch η) exp(-κ_path Σ)
# ============================================================
print(f"\n\n{'=' * 100}")
print("TWO-OPERATOR MODEL: V(z) = V_0 · exp(-κ_epoch·z) · exp(-κ_path·Σ)")
print("(Does adding path operator improve on epoch-only?)")
print("=" * 100)

def two_operator(X, V0, k_epoch, k_path):
    z, sig = X
    return V0 * np.exp(-k_epoch * z) * np.exp(-k_path * sig)

try:
    X_data = np.array([eta_z, eta_sigma])
    popt2, _ = curve_fit(two_operator, X_data, V_arr,
                         p0=[V_arr[0], 1.0, 1e-23],
                         bounds=([0, 0, 0], [100, 50, 1e-18]),
                         maxfev=10000)
    pred2 = two_operator(X_data, *popt2)
    ss_res2 = np.sum((V_arr - pred2)**2)
    ss_tot2 = np.sum((V_arr - np.mean(V_arr))**2)
    r2_two = 1 - ss_res2 / ss_tot2
    
    # Compare to epoch-only
    r2_epoch, k_epoch = fit_and_score(eta_z, V_arr, "V", "z")
    
    print(f"\n  Epoch-only:     R² = {r2_epoch:.4f}")
    print(f"  Two-operator:   R² = {r2_two:.4f}")
    print(f"  Improvement:    ΔR² = {r2_two - r2_epoch:+.4f}")
    print(f"\n  κ_epoch = {popt2[1]:.4f}")
    print(f"  κ_path  = {popt2[2]:.2e}")
    
    if r2_two > r2_epoch + 0.01:
        print(f"\n  🔥 Two-operator model is BETTER → both epoch and path contribute")
    else:
        print(f"\n  Epoch-only suffices for SN Ia → path contribution negligible here")
        print(f"  (Path operator may still matter for quasars/FRBs)")
except Exception as e:
    print(f"\n  Two-operator fit failed: {e}")

# ============================================================
# SAVE
# ============================================================
results_dir = Path("results_epoch_variable")
results_dir.mkdir(exist_ok=True)

results = {
    'bins': [{k: float(v) if isinstance(v, (int, float, np.floating)) else v 
              for k, v in b.items()} for b in bins],
    'best_variable': best_name,
    'best_avg_r2': float(best_r2),
}

with open(results_dir / "epoch_variable.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {results_dir}/")
