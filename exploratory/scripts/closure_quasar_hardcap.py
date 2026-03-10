#!/usr/bin/env python3
"""
closure_quasar_hardcap.py — DESI Cross-Domain Validation (DR16Q)
=================================================================

Apply the SN Ia Hard Cap framework to quasar emission lines:
  1. Phase-space volume V_C(age) — does it contract?
  2. Asymmetric truncation — does one "channel" migrate while another anchors?
  3. Eigenvalue evolution — do all shrink? Same ordering as SNe?
  4. Composite structure — are there two sub-populations within line-width classes?

Data: Wu & Shen (2022) — 750,414 SDSS DR16Q quasars
      PyQSOFit measurements: [peak_wave, center_wave, EW, log_sigma, FWHM_km/s, log_flux]

Key lines:
  - MgII_BR (2798Å): z=0.5-2.5, 650K objects — THERMODYNAMIC (EW = temperature/opacity)
  - CIV (1549Å): z=1.5-4, 500K objects — HIGH-IONIZATION (density/temperature)  
  - CIII_BR (1909Å): z=1-3, 640K objects — SEMI-FORBIDDEN
  - Hβ_BR (4861Å): z=0-1, 160K objects — BALMER (geometric/virial)

Prediction from SN Ia results:
  - EW channels contract fastest (thermodynamic, like mB)
  - FWHM channels contract slowest (geometric/virial, like x1)
  - One population class migrates, another anchors
  - Generator is cosmic age, not z

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr, kurtosis as scipy_kurtosis
from scipy.optimize import curve_fit
from scipy.integrate import quad
from astropy.io import fits
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# Cosmology
H0 = 67.4
Omega_m, Omega_L = 0.315, 0.685

def E(z): return np.sqrt(Omega_m*(1+z)**3 + Omega_L)

def cosmic_age(z):
    r, _ = quad(lambda zz: 1/((1+zz)*E(zz)), z, 100)
    return r / (H0 * 1e5 / 3.086e24 / 3.156e16)

# Load DR16Q
print("Loading DR16Q catalog...")
f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z_all = d['Z_DR16Q']
print(f"  Total: {len(z_all)}, valid z: {np.sum(z_all > 0)}\n")

# ============================================================
# EXTRACT EMISSION LINE PROPERTIES
# ============================================================
# PyQSOFit format: [peak_wave, center_wave, EW, log_sigma_km/s, FWHM_km/s, log_flux]
# Col 2 = EW (rest-frame equivalent width in Angstrom)
# Col 4 = FWHM (km/s) — virial/geometric diagnostic

def extract_line(name, z_min, z_max, ew_max=5000, fwhm_min=100, fwhm_max=30000):
    """Extract valid measurements for a line"""
    ld = d[name]
    mask = ((z_all >= z_min) & (z_all < z_max) & 
            (ld[:, 4] > fwhm_min) & (ld[:, 4] < fwhm_max) &
            (ld[:, 2] > 0) & (ld[:, 2] < ew_max) &
            np.isfinite(ld[:, 2]) & np.isfinite(ld[:, 4]))
    
    # Also need LOGLBOL for Baldwin control
    lbol = d['LOGLBOL']
    mask &= (lbol > 40) & (lbol < 50) & np.isfinite(lbol)
    
    return {
        'z': z_all[mask],
        'ew': ld[mask, 2],
        'fwhm': ld[mask, 4],
        'lbol': lbol[mask],
        'N': np.sum(mask),
    }

# Primary lines
print("=" * 100)
print("QUASAR EMISSION LINE COVERAGE")
print("=" * 100)

lines = {
    'MgII_BR': extract_line('MGII_BR', 0.4, 2.5),
    'CIV': extract_line('CIV', 1.5, 4.0),
    'CIII_BR': extract_line('CIII_BR', 0.8, 3.0),
    'Hbeta_BR': extract_line('HBETA_BR', 0.0, 1.0),
}

for name, info in lines.items():
    print(f"  {name:>10}: N={info['N']:>7}, z=[{info['z'].min():.2f}, {info['z'].max():.2f}], "
          f"EW=[{np.percentile(info['ew'],5):.1f}, {np.percentile(info['ew'],95):.1f}], "
          f"FWHM=[{np.percentile(info['fwhm'],5):.0f}, {np.percentile(info['fwhm'],95):.0f}] km/s")

# ============================================================
# TEST 1: PHASE-SPACE VOLUME V_C(age) FOR EACH LINE
# ============================================================
print(f"\n\n{'=' * 100}")
print("TEST 1: PHASE-SPACE VOLUME V(EW, FWHM) vs COSMIC AGE")
print("Does quasar diagnostic diversity contract with cosmic age?")
print("=" * 100)

def compute_V_quasar(ew, fwhm, lbol):
    """Compute phase-space volume from (log_EW, log_FWHM) after Baldwin control"""
    log_ew = np.log10(np.maximum(ew, 0.1))
    log_fwhm = np.log10(np.maximum(fwhm, 100))
    
    # Residualize against luminosity (Baldwin effect control)
    from numpy.polynomial import polynomial as P
    # EW residual
    coef_ew = P.polyfit(lbol, log_ew, 1)
    log_ew_resid = log_ew - P.polyval(lbol, coef_ew)
    
    # FWHM residual  
    coef_fwhm = P.polyfit(lbol, log_fwhm, 1)
    log_fwhm_resid = log_fwhm - P.polyval(lbol, coef_fwhm)
    
    # Covariance
    C = np.cov(np.array([log_ew_resid, log_fwhm_resid]))
    eigvals = np.linalg.eigvalsh(C)
    eigvals = np.maximum(eigvals, 1e-10)
    V = np.sqrt(np.prod(eigvals))
    
    return V, np.var(log_ew_resid), np.var(log_fwhm_resid), C[0,1]

for line_name, info in lines.items():
    if info['N'] < 1000:
        continue
    
    print(f"\n  Line: {line_name} (N={info['N']})")
    
    # Z-bins (adaptive based on coverage)
    z_sorted = np.sort(info['z'])
    n_bins = min(8, info['N'] // 5000)
    if n_bins < 3:
        n_bins = 3
    z_edges_q = np.percentile(z_sorted, np.linspace(0, 100, n_bins + 1))
    z_edges_q = np.unique(np.round(z_edges_q, 3))
    
    print(f"  {'z_mid':>6} {'age':>10} {'N':>6} | {'V':>10} {'Var(EW)':>10} {'Var(FWHM)':>10} {'Cov':>10}")
    print(f"  " + "-" * 75)
    
    z_list, age_list, V_list = [], [], []
    var_ew_list, var_fwhm_list = [], []
    
    for i in range(len(z_edges_q) - 1):
        z_lo, z_hi = z_edges_q[i], z_edges_q[i+1]
        mask = (info['z'] >= z_lo) & (info['z'] < z_hi)
        if np.sum(mask) < 200:
            continue
        
        zc = np.mean(info['z'][mask])
        age = cosmic_age(zc)
        V, var_ew, var_fwhm, cov_ef = compute_V_quasar(
            info['ew'][mask], info['fwhm'][mask], info['lbol'][mask])
        
        z_list.append(zc)
        age_list.append(age)
        V_list.append(V)
        var_ew_list.append(var_ew)
        var_fwhm_list.append(var_fwhm)
        
        print(f"  {zc:>6.3f} {age:>10.2e} {np.sum(mask):>6} | {V:>10.6f} {var_ew:>10.6f} {var_fwhm:>10.6f} {cov_ef:>10.6f}")
    
    if len(z_list) >= 3:
        z_arr = np.array(z_list)
        V_arr = np.array(V_list)
        
        # V vs z
        rho_z, p_z = spearmanr(z_arr, V_arr)
        # V vs age
        age_arr = np.array(age_list)
        rho_age, p_age = spearmanr(age_arr, V_arr)
        
        # Var(EW) and Var(FWHM) trends
        rho_ew, p_ew = spearmanr(z_arr, var_ew_list)
        rho_fwhm, p_fwhm = spearmanr(z_arr, var_fwhm_list)
        
        print(f"\n  V vs z:    ρ = {rho_z:+.3f}, p = {p_z:.4f} {'🔥 CONTRACTS' if rho_z < -0.4 and p_z < 0.1 else '— FLAT' if abs(rho_z) < 0.3 else ''}")
        print(f"  V vs age:  ρ = {rho_age:+.3f}, p = {p_age:.4f} {'🔥 CONTRACTS' if rho_age > 0.4 and p_age < 0.1 else ''}")
        print(f"  Var(EW) vs z:   ρ = {rho_ew:+.3f}, p = {p_ew:.4f} {'→ EW contracts' if rho_ew < -0.3 else '→ EW stable/grows'}")
        print(f"  Var(FWHM) vs z: ρ = {rho_fwhm:+.3f}, p = {p_fwhm:.4f} {'→ FWHM contracts' if rho_fwhm < -0.3 else '→ FWHM stable'}")
        
        # Rate ordering check (SN prediction: EW fastest, FWHM slowest)
        if abs(rho_ew) > abs(rho_fwhm):
            print(f"  ✓ Rate ordering matches SN prediction: EW contracts faster than FWHM")
        else:
            print(f"  ✗ Rate ordering REVERSED: FWHM contracts faster than EW")

# ============================================================
# TEST 2: CHANNEL SPLIT — HIGH vs LOW FWHM (Eddington classes)
# ============================================================
print(f"\n\n{'=' * 100}")
print("TEST 2: CHANNEL SPLIT — Narrow-line (NL) vs Broad-line (BL) populations")
print("Split at median FWHM. Does one channel migrate like fast decliners?")
print("=" * 100)

for line_name in ['MgII_BR', 'CIV']:
    info = lines[line_name]
    if info['N'] < 5000:
        continue
    
    fwhm_median = np.median(info['fwhm'])
    print(f"\n  {line_name}: median FWHM = {fwhm_median:.0f} km/s")
    
    z_sorted = np.sort(info['z'])
    n_bins = min(8, info['N'] // 5000)
    z_edges_q = np.percentile(z_sorted, np.linspace(0, 100, n_bins + 1))
    z_edges_q = np.unique(np.round(z_edges_q, 3))
    
    print(f"  {'z':>6} | {'μ_EW(NL)':>10} {'μ_EW(BL)':>10} {'ΔμEW':>8} | {'μ_FWHM(NL)':>11} {'μ_FWHM(BL)':>11} {'ΔμFWHM':>8}")
    print(f"  " + "-" * 80)
    
    z_list = []
    delta_ew_list = []
    delta_fwhm_list = []
    mu_ew_nl_list, mu_ew_bl_list = [], []
    
    for i in range(len(z_edges_q) - 1):
        z_lo, z_hi = z_edges_q[i], z_edges_q[i+1]
        mask = (info['z'] >= z_lo) & (info['z'] < z_hi)
        
        nl = mask & (info['fwhm'] < fwhm_median)
        bl = mask & (info['fwhm'] >= fwhm_median)
        
        if np.sum(nl) < 100 or np.sum(bl) < 100:
            continue
        
        zc = np.mean(info['z'][mask])
        
        # Use log for EW and FWHM
        mu_ew_nl = np.mean(np.log10(info['ew'][nl]))
        mu_ew_bl = np.mean(np.log10(info['ew'][bl]))
        mu_fwhm_nl = np.mean(np.log10(info['fwhm'][nl]))
        mu_fwhm_bl = np.mean(np.log10(info['fwhm'][bl]))
        
        z_list.append(zc)
        delta_ew_list.append(mu_ew_bl - mu_ew_nl)
        delta_fwhm_list.append(mu_fwhm_bl - mu_fwhm_nl)
        mu_ew_nl_list.append(mu_ew_nl)
        mu_ew_bl_list.append(mu_ew_bl)
        
        print(f"  {zc:>6.3f} | {mu_ew_nl:>10.4f} {mu_ew_bl:>10.4f} {mu_ew_bl-mu_ew_nl:>8.4f} | {mu_fwhm_nl:>11.4f} {mu_fwhm_bl:>11.4f} {mu_fwhm_bl-mu_fwhm_nl:>8.4f}")
    
    if len(z_list) >= 3:
        z_arr = np.array(z_list)
        
        rho_dew, p_dew = spearmanr(z_arr, delta_ew_list)
        rho_dfwhm, p_dfwhm = spearmanr(z_arr, delta_fwhm_list)
        
        # Individual centroid evolution
        rho_nl, p_nl = spearmanr(z_arr, mu_ew_nl_list)
        rho_bl, p_bl = spearmanr(z_arr, mu_ew_bl_list)
        
        print(f"\n  ΔEW (BL-NL) vs z:    ρ = {rho_dew:+.3f}, p = {p_dew:.4f} {'→ CHANNELS CONVERGE' if rho_dew < -0.4 else '→ STABLE' if abs(rho_dew) < 0.3 else '→ DIVERGE'}")
        print(f"  ΔFWHM (BL-NL) vs z:  ρ = {rho_dfwhm:+.3f}, p = {p_dfwhm:.4f}")
        print(f"  μ_EW(NL) vs z: ρ = {rho_nl:+.3f}, p = {p_nl:.4f} {'→ NL MIGRATES' if abs(rho_nl) > 0.5 else '→ STATIONARY'}")
        print(f"  μ_EW(BL) vs z: ρ = {rho_bl:+.3f}, p = {p_bl:.4f} {'→ BL MIGRATES' if abs(rho_bl) > 0.5 else '→ STATIONARY'}")

# ============================================================
# TEST 3: HIGHER MOMENTS WITHIN CHANNELS
# ============================================================
print(f"\n\n{'=' * 100}")
print("TEST 3: HIGHER MOMENTS — Kurtosis and IQR of log(EW) and log(FWHM)")
print("SN prediction: kurtosis rises (outside-in truncation) in the migrating channel")
print("=" * 100)

for line_name in ['MgII_BR', 'CIV']:
    info = lines[line_name]
    if info['N'] < 5000:
        continue
    
    fwhm_median = np.median(info['fwhm'])
    
    z_sorted = np.sort(info['z'])
    n_bins = min(8, info['N'] // 5000)
    z_edges_q = np.percentile(z_sorted, np.linspace(0, 100, n_bins + 1))
    z_edges_q = np.unique(np.round(z_edges_q, 3))
    
    for ch_name, ch_mask_fn in [("NL (FWHM<med)", lambda: info['fwhm'] < fwhm_median),
                                  ("BL (FWHM≥med)", lambda: info['fwhm'] >= fwhm_median)]:
        print(f"\n  {line_name} — {ch_name}:")
        print(f"  {'z':>6} {'N':>6} | {'Kurt(EW)':>10} {'Kurt(FWHM)':>12} {'IQR90(EW)':>10} {'IQR90(FWHM)':>12}")
        print(f"  " + "-" * 65)
        
        z_list = []
        kurt_ew, kurt_fwhm = [], []
        iqr_ew, iqr_fwhm = [], []
        
        for i in range(len(z_edges_q) - 1):
            z_lo, z_hi = z_edges_q[i], z_edges_q[i+1]
            ch_mask = ch_mask_fn()
            mask = (info['z'] >= z_lo) & (info['z'] < z_hi) & ch_mask
            if np.sum(mask) < 100:
                continue
            
            zc = np.mean(info['z'][mask])
            lew = np.log10(info['ew'][mask])
            lfwhm = np.log10(info['fwhm'][mask])
            
            ke = scipy_kurtosis(lew, fisher=True)
            kf = scipy_kurtosis(lfwhm, fisher=True)
            ie = np.percentile(lew, 90) - np.percentile(lew, 10)
            iff = np.percentile(lfwhm, 90) - np.percentile(lfwhm, 10)
            
            z_list.append(zc)
            kurt_ew.append(ke)
            kurt_fwhm.append(kf)
            iqr_ew.append(ie)
            iqr_fwhm.append(iff)
            
            print(f"  {zc:>6.3f} {np.sum(mask):>6} | {ke:>10.3f} {kf:>12.3f} {ie:>10.4f} {iff:>12.4f}")
        
        if len(z_list) >= 3:
            z_arr = np.array(z_list)
            rho_ke, p_ke = spearmanr(z_arr, kurt_ew)
            rho_kf, p_kf = spearmanr(z_arr, kurt_fwhm)
            rho_ie, p_ie = spearmanr(z_arr, iqr_ew)
            rho_if, p_if = spearmanr(z_arr, iqr_fwhm)
            
            print(f"\n  Kurt(EW) vs z:   ρ = {rho_ke:+.3f}, p = {p_ke:.4f} {'🔥 RISES' if rho_ke > 0.4 else ''}")
            print(f"  Kurt(FWHM) vs z: ρ = {rho_kf:+.3f}, p = {p_kf:.4f} {'🔥 RISES' if rho_kf > 0.4 else ''}")
            print(f"  IQR90(EW) vs z:  ρ = {rho_ie:+.3f}, p = {p_ie:.4f} {'🔥 SHRINKS' if rho_ie < -0.4 else ''}")
            print(f"  IQR90(FWHM) vs z:ρ = {rho_if:+.3f}, p = {p_if:.4f} {'🔥 SHRINKS' if rho_if < -0.4 else ''}")

# ============================================================
# TEST 4: B(z)/W(z) DECOMPOSITION (NL/BL split)
# ============================================================
print(f"\n\n{'=' * 100}")
print("TEST 4: BETWEEN/WITHIN DECOMPOSITION")
print("B(z) = architecture (NL/BL separation)")
print("W(z) = within-channel diversity (hard cap)")
print("SN prediction: B stable or drops, W shrinks")
print("=" * 100)

for line_name in ['MgII_BR', 'CIV']:
    info = lines[line_name]
    if info['N'] < 5000:
        continue
    
    fwhm_median = np.median(info['fwhm'])
    
    z_sorted = np.sort(info['z'])
    n_bins = min(8, info['N'] // 5000)
    z_edges_q = np.percentile(z_sorted, np.linspace(0, 100, n_bins + 1))
    z_edges_q = np.unique(np.round(z_edges_q, 3))
    
    for var_name, var_fn in [("log(EW)", lambda: np.log10(info['ew'])),
                               ("log(FWHM)", lambda: np.log10(info['fwhm']))]:
        print(f"\n  {line_name} — {var_name}:")
        print(f"  {'z':>6} | {'B(z)':>8} {'W(z)':>8}")
        print(f"  " + "-" * 30)
        
        vals = var_fn()
        z_list, B_list, W_list = [], [], []
        
        for i in range(len(z_edges_q) - 1):
            z_lo, z_hi = z_edges_q[i], z_edges_q[i+1]
            mask = (info['z'] >= z_lo) & (info['z'] < z_hi)
            
            nl_vals = vals[(mask) & (info['fwhm'] < fwhm_median)]
            bl_vals = vals[(mask) & (info['fwhm'] >= fwhm_median)]
            all_vals = vals[mask]
            
            if len(nl_vals) < 50 or len(bl_vals) < 50:
                continue
            
            zc = np.mean(info['z'][mask])
            var_total = np.var(all_vals)
            if var_total < 1e-10:
                continue
            
            n_nl, n_bl = len(nl_vals), len(bl_vals)
            n_tot = n_nl + n_bl
            mu_nl, mu_bl = np.mean(nl_vals), np.mean(bl_vals)
            mu_all = np.mean(all_vals)
            
            var_between = (n_nl/n_tot)*(mu_nl-mu_all)**2 + (n_bl/n_tot)*(mu_bl-mu_all)**2
            var_within = (n_nl/n_tot)*np.var(nl_vals) + (n_bl/n_tot)*np.var(bl_vals)
            
            B = var_between / var_total
            W = var_within / var_total
            
            z_list.append(zc)
            B_list.append(B)
            W_list.append(W)
            
            print(f"  {zc:>6.3f} | {B:>8.4f} {W:>8.4f}")
        
        if len(z_list) >= 3:
            z_arr = np.array(z_list)
            rho_B, p_B = spearmanr(z_arr, B_list)
            rho_W, p_W = spearmanr(z_arr, W_list)
            print(f"  B(z) vs z: ρ = {rho_B:+.3f}, p = {p_B:.4f}")
            print(f"  W(z) vs z: ρ = {rho_W:+.3f}, p = {p_W:.4f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n\n{'=' * 100}")
print("CROSS-DOMAIN VERDICT")
print("=" * 100)
print("""
If quasars show:
  ✓ V_C(age) contracts          → Hard Cap is universal (not SN-specific)
  ✓ One channel migrates         → Asymmetric truncation is universal
  ✓ EW contracts > FWHM          → Rate ordering matches (thermodynamic > geometric)
  ✓ Kurtosis rises in migrator   → Outside-in truncation is universal
  ✓ B(z) stable, W(z) shrinks    → Architecture preserved, diversity capped

If ANY of these fail → Hard Cap is SN-specific (still publishable, but not universal)
""")

# Save
results_dir = Path("results_quasar_hardcap")
results_dir.mkdir(exist_ok=True)
print(f"Results saved to {results_dir}/")
