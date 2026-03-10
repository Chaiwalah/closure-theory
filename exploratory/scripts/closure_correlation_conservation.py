#!/usr/bin/env python3
"""
closure_correlation_conservation.py — Correlation Conservation Test
=====================================================================

Hypothesis (from EPR/entanglement analogy):
  Individual variances can change with z, but the CORRELATION STRUCTURE
  between observables is conserved. Information redistributes between
  channels but the total mutual information is constant.

Tests:
  1. Mutual Information MI(EW, FWHM) vs z — is it conserved?
  2. Product (Var × peak_density) vs z — étendue conservation?
  3. Total correlation = det(C) / prod(marginals) vs z
  4. Entropy of the joint distribution vs z
  5. Cross-line MI: MI(EW_line1, EW_line2) vs z — are cross-line
     correlations conserved even as individual variances change?

For SNe Ia:
  MI(mB, x1, c) vs z — same test on the original domain

Author: Closure Theory Collaboration
Date: 2026-03-05
"""

import numpy as np
from scipy.stats import spearmanr
from astropy.io import fits
from numpy.polynomial import polynomial as P
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# MUTUAL INFORMATION ESTIMATOR (KDE-based for continuous vars)
# ============================================================
def mutual_information_kde(x, y, n_bins=30):
    """Estimate MI using binned histogram (fast, robust for large N)."""
    # Standardize
    x = (x - np.mean(x)) / (np.std(x) + 1e-10)
    y = (y - np.mean(y)) / (np.std(y) + 1e-10)
    
    # 2D histogram
    H_xy, _, _ = np.histogram2d(x, y, bins=n_bins, density=True)
    H_x, _ = np.histogram(x, bins=n_bins, density=True)
    H_y, _ = np.histogram(y, bins=n_bins, density=True)
    
    # Normalize to probabilities
    dx = (x.max() - x.min()) / n_bins
    dy = (y.max() - y.min()) / n_bins
    
    p_xy = H_xy * dx * dy
    p_x = H_x * dx
    p_y = H_y * dy
    
    # MI = sum p(x,y) log(p(x,y) / (p(x)*p(y)))
    mi = 0
    for i in range(n_bins):
        for j in range(n_bins):
            if p_xy[i, j] > 1e-10 and p_x[i] > 1e-10 and p_y[j] > 1e-10:
                mi += p_xy[i, j] * np.log(p_xy[i, j] / (p_x[i] * p_y[j]))
    
    return max(mi, 0)


def differential_entropy(x, n_bins=30):
    """Estimate differential entropy."""
    H, edges = np.histogram(x, bins=n_bins, density=True)
    dx = edges[1] - edges[0]
    p = H * dx
    ent = -np.sum(p[p > 1e-10] * np.log(p[p > 1e-10]))
    return ent


def joint_entropy(x, y, n_bins=30):
    """Estimate joint entropy."""
    H, _, _ = np.histogram2d(x, y, bins=n_bins, density=True)
    dx = (x.max() - x.min()) / n_bins
    dy = (y.max() - y.min()) / n_bins
    p = H * dx * dy
    ent = -np.sum(p[p > 1e-10] * np.log(p[p > 1e-10]))
    return ent


# ============================================================
# PART 1: QUASAR CORRELATION CONSERVATION
# ============================================================
print("Loading DR16Q catalog...")
f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z_all = d['Z_DR16Q']
lbol = d['LOGLBOL']

print(f"\n{'=' * 110}")
print("PART 1: QUASAR MI(EW, FWHM) CONSERVATION — Does correlation structure survive while marginals change?")
print("=" * 110)

for line_name, col_name, z_min, z_max in [
    ('MgII_BR', 'MGII_BR', 0.4, 2.5),
    ('CIV', 'CIV', 1.5, 4.0),
    ('CIII_BR', 'CIII_BR', 0.8, 3.0),
    ('Hβ_BR', 'HBETA_BR', 0.3, 1.0),
]:
    ld = d[col_name]
    mask = ((z_all >= z_min) & (z_all < z_max) &
            (ld[:, 4] > 100) & (ld[:, 4] < 30000) &
            (ld[:, 2] > 0) & (ld[:, 2] < 5000) &
            np.isfinite(ld[:, 2]) & np.isfinite(ld[:, 4]) &
            (lbol > 40) & (lbol < 50))
    
    N = np.sum(mask)
    if N < 5000:
        continue
    
    z_v = z_all[mask]
    log_ew = np.log10(ld[mask, 2])
    log_fwhm = np.log10(ld[mask, 4])
    
    # Residualize against L_bol
    lbol_v = lbol[mask]
    coef_e = P.polyfit(lbol_v, log_ew, 1)
    log_ew_r = log_ew - P.polyval(lbol_v, coef_e)
    coef_f = P.polyfit(lbol_v, log_fwhm, 1)
    log_fwhm_r = log_fwhm - P.polyval(lbol_v, coef_f)
    
    n_bins = min(8, N // 5000)
    if n_bins < 3:
        n_bins = 3
    z_edges = np.percentile(z_v, np.linspace(0, 100, n_bins + 1))
    z_edges = np.unique(np.round(z_edges, 3))
    
    print(f"\n  {line_name} (N={N}):")
    print(f"  {'z':>6} {'N':>6} | {'MI(EW,FWHM)':>12} {'Var(EW)':>10} {'Var(FWHM)':>10} {'H(joint)':>10} {'H(EW)':>8} {'H(FWHM)':>8} | {'Var×Peak':>10}")
    print(f"  " + "-" * 105)
    
    z_list, mi_list, var_ew_list, var_fwhm_list = [], [], [], []
    h_joint_list, etendue_list = [], []
    
    for i in range(len(z_edges) - 1):
        z_lo, z_hi = z_edges[i], z_edges[i+1]
        bm = (z_v >= z_lo) & (z_v < z_hi)
        n = np.sum(bm)
        if n < 500:
            continue
        
        zc = np.mean(z_v[bm])
        ew_bin = log_ew_r[bm]
        fwhm_bin = log_fwhm_r[bm]
        
        mi = mutual_information_kde(ew_bin, fwhm_bin)
        var_e = np.var(ew_bin)
        var_f = np.var(fwhm_bin)
        h_j = joint_entropy(ew_bin, fwhm_bin)
        h_e = differential_entropy(ew_bin)
        h_f = differential_entropy(fwhm_bin)
        
        # Étendue: variance × peak density (inverse of spread × concentration)
        from scipy.stats import gaussian_kde
        try:
            kde_e = gaussian_kde(ew_bin)
            peak_density = kde_e.evaluate(np.array([np.mean(ew_bin)]))[0]
            etendue = var_e * peak_density
        except:
            etendue = 0
        
        z_list.append(zc)
        mi_list.append(mi)
        var_ew_list.append(var_e)
        var_fwhm_list.append(var_f)
        h_joint_list.append(h_j)
        etendue_list.append(etendue)
        
        print(f"  {zc:>6.3f} {n:>6} | {mi:>12.4f} {var_e:>10.6f} {var_f:>10.6f} {h_j:>10.4f} {h_e:>8.4f} {h_f:>8.4f} | {etendue:>10.4f}")
    
    if len(z_list) >= 3:
        z_arr = np.array(z_list)
        rho_mi, p_mi = spearmanr(z_arr, mi_list)
        rho_ve, p_ve = spearmanr(z_arr, var_ew_list)
        rho_vf, p_vf = spearmanr(z_arr, var_fwhm_list)
        rho_hj, p_hj = spearmanr(z_arr, h_joint_list)
        rho_et, p_et = spearmanr(z_arr, etendue_list)
        
        print(f"\n  MI(EW,FWHM) vs z:  ρ = {rho_mi:+.3f}, p = {p_mi:.4f} {'🔥 CONSERVED' if abs(rho_mi) < 0.3 else '↗ GROWS' if rho_mi > 0.3 else '↘ DROPS'}")
        print(f"  Var(EW) vs z:      ρ = {rho_ve:+.3f}, p = {p_ve:.4f}")
        print(f"  Var(FWHM) vs z:    ρ = {rho_vf:+.3f}, p = {p_vf:.4f}")
        print(f"  H(joint) vs z:     ρ = {rho_hj:+.3f}, p = {p_hj:.4f} {'🔥 CONSERVED' if abs(rho_hj) < 0.3 else ''}")
        print(f"  Var×Peak vs z:     ρ = {rho_et:+.3f}, p = {p_et:.4f} {'🔥 CONSERVED' if abs(rho_et) < 0.3 else ''}")


# ============================================================
# PART 2: CROSS-LINE MI CONSERVATION
# ============================================================
print(f"\n\n{'=' * 110}")
print("PART 2: CROSS-LINE MI — Is MI between different lines conserved?")
print("If individual variances change but cross-line correlations don't → information redistribution")
print("=" * 110)

# MgII vs CIII (z=0.8-2.5)
mgii = d['MGII_BR']
ciii = d['CIII_BR']

mask_mc = ((z_all >= 0.8) & (z_all < 2.5) &
           (mgii[:, 2] > 0) & (mgii[:, 4] > 100) & (mgii[:, 4] < 30000) &
           (ciii[:, 2] > 0) & (ciii[:, 4] > 100) & (ciii[:, 4] < 30000) &
           np.isfinite(mgii[:, 2]) & np.isfinite(ciii[:, 2]) &
           (lbol > 40) & (lbol < 50))

z_mc = z_all[mask_mc]
log_ew_mgii = np.log10(mgii[mask_mc, 2])
log_ew_ciii = np.log10(ciii[mask_mc, 2])

# Residualize
lbol_mc = lbol[mask_mc]
coef = P.polyfit(lbol_mc, log_ew_mgii, 1)
ew_mgii_r = log_ew_mgii - P.polyval(lbol_mc, coef)
coef = P.polyfit(lbol_mc, log_ew_ciii, 1)
ew_ciii_r = log_ew_ciii - P.polyval(lbol_mc, coef)

n_bins = min(8, np.sum(mask_mc) // 5000)
z_edges = np.percentile(z_mc, np.linspace(0, 100, n_bins + 1))
z_edges = np.unique(np.round(z_edges, 3))

print(f"\n  MgII_EW vs CIII_EW (N={np.sum(mask_mc)}):")
print(f"  {'z':>6} {'N':>6} | {'MI(MgII,CIII)':>14} {'Var(MgII)':>10} {'Var(CIII)':>10} {'r(Pearson)':>11}")
print(f"  " + "-" * 65)

z_list, mi_list = [], []
for i in range(len(z_edges) - 1):
    z_lo, z_hi = z_edges[i], z_edges[i+1]
    bm = (z_mc >= z_lo) & (z_mc < z_hi)
    n = np.sum(bm)
    if n < 500:
        continue
    
    zc = np.mean(z_mc[bm])
    mi = mutual_information_kde(ew_mgii_r[bm], ew_ciii_r[bm])
    ve = np.var(ew_mgii_r[bm])
    vc = np.var(ew_ciii_r[bm])
    r = np.corrcoef(ew_mgii_r[bm], ew_ciii_r[bm])[0, 1]
    
    z_list.append(zc)
    mi_list.append(mi)
    print(f"  {zc:>6.3f} {n:>6} | {mi:>14.4f} {ve:>10.6f} {vc:>10.6f} {r:>11.4f}")

if len(z_list) >= 3:
    rho, p = spearmanr(z_list, mi_list)
    print(f"\n  MI(MgII,CIII) vs z: ρ = {rho:+.3f}, p = {p:.4f} {'🔥 CONSERVED' if abs(rho) < 0.3 else '↗ GROWS' if rho > 0.3 else '↘ DROPS'}")


# ============================================================
# PART 3: SNe Ia CORRELATION CONSERVATION
# ============================================================
print(f"\n\n{'=' * 110}")
print("PART 3: SNe Ia — MI(mB, x1, c) vs z")
print("The original domain: does the SN correlation structure redistribute while conserving total MI?")
print("=" * 110)

sn = np.loadtxt('data/pantheon_plus.dat', dtype=str)
header = sn[0]
data = sn[1:]

col_idx = {h: i for i, h in enumerate(header)}
z_sn = data[:, col_idx['zHD']].astype(float)
mb_sn = data[:, col_idx['mB']].astype(float)
x1_sn = data[:, col_idx['x1']].astype(float)
c_sn = data[:, col_idx['c']].astype(float)

# Standard cuts
mask_sn = (z_sn > 0.01) & (z_sn < 2.5) & (np.abs(c_sn) < 0.3) & (np.abs(x1_sn) < 3)
z_sn = z_sn[mask_sn]
mb_sn = mb_sn[mask_sn]
x1_sn = x1_sn[mask_sn]
c_sn = c_sn[mask_sn]

print(f"  N = {len(z_sn)}")

z_edges_sn = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 2.5]

print(f"  {'z':>6} {'N':>4} | {'MI(mB,x1)':>10} {'MI(mB,c)':>10} {'MI(x1,c)':>10} | {'Var(mB)':>8} {'Var(x1)':>8} {'Var(c)':>8} | {'MI_total':>10}")
print(f"  " + "-" * 95)

z_list, mi_total_list = [], []
mi_mb_x1_list, mi_mb_c_list, mi_x1_c_list = [], [], []

for i in range(len(z_edges_sn) - 1):
    z_lo, z_hi = z_edges_sn[i], z_edges_sn[i+1]
    bm = (z_sn >= z_lo) & (z_sn < z_hi)
    n = np.sum(bm)
    if n < 30:
        continue
    
    zc = np.mean(z_sn[bm])
    
    mi_1 = mutual_information_kde(mb_sn[bm], x1_sn[bm], n_bins=15)
    mi_2 = mutual_information_kde(mb_sn[bm], c_sn[bm], n_bins=15)
    mi_3 = mutual_information_kde(x1_sn[bm], c_sn[bm], n_bins=15)
    mi_tot = mi_1 + mi_2 + mi_3
    
    vm = np.var(mb_sn[bm])
    vx = np.var(x1_sn[bm])
    vc = np.var(c_sn[bm])
    
    z_list.append(zc)
    mi_total_list.append(mi_tot)
    mi_mb_x1_list.append(mi_1)
    mi_mb_c_list.append(mi_2)
    mi_x1_c_list.append(mi_3)
    
    print(f"  {zc:>6.3f} {n:>4} | {mi_1:>10.4f} {mi_2:>10.4f} {mi_3:>10.4f} | {vm:>8.4f} {vx:>8.4f} {vc:>8.4f} | {mi_tot:>10.4f}")

if len(z_list) >= 3:
    z_arr = np.array(z_list)
    rho_tot, p_tot = spearmanr(z_arr, mi_total_list)
    rho_1, p_1 = spearmanr(z_arr, mi_mb_x1_list)
    rho_2, p_2 = spearmanr(z_arr, mi_mb_c_list)
    rho_3, p_3 = spearmanr(z_arr, mi_x1_c_list)
    
    print(f"\n  MI(mB,x1) vs z:  ρ = {rho_1:+.3f}, p = {p_1:.4f}")
    print(f"  MI(mB,c) vs z:   ρ = {rho_2:+.3f}, p = {p_2:.4f}")
    print(f"  MI(x1,c) vs z:   ρ = {rho_3:+.3f}, p = {p_3:.4f}")
    print(f"  MI_total vs z:   ρ = {rho_tot:+.3f}, p = {p_tot:.4f} {'🔥 CONSERVED' if abs(rho_tot) < 0.3 else '↗ GROWS' if rho_tot > 0.3 else '↘ DROPS'}")


# ============================================================
# VERDICT
# ============================================================
print(f"\n\n{'=' * 110}")
print("CORRELATION CONSERVATION VERDICT")
print("=" * 110)
print("""
The EPR analogy predicts:
  - Individual variances can change (marginals redistribute)
  - But MUTUAL INFORMATION (correlation structure) is CONSERVED
  - Total MI = constant across z, even as channels shift
  
If MI is conserved: the firecracker redistributes information, doesn't create or destroy it
If MI drops: information is genuinely lost (decoherence / Hard Cap destroys correlations)  
If MI grows: information is being CONCENTRATED (compression increases density)
""")
