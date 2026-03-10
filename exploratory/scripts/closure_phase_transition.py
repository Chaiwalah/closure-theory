#!/usr/bin/env python3
"""
PHASE TRANSITION TESTS — Is the sigmoid a real phase transition?

If the lattice-to-liquid analogy is literal, not metaphorical:
1. Critical fluctuations: variance of correlations should PEAK at z₀
2. Critical exponent: correlation ~ |z - z₀|^β near the transition
3. Universality: different source classes share critical exponents
4. Correlation length: at the boundary, how far does order persist?
5. Order parameter: what IS the order parameter of this transition?
"""

import numpy as np
from astropy.io import fits
from scipy import stats
from scipy.optimize import curve_fit
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_phase_transition', exist_ok=True)

# ============================================================
# Load data
# ============================================================
print("=" * 70)
print("PHASE TRANSITION TESTS — Is the Sigmoid a Melting Point?")
print("=" * 70)

# DR16Q quasars
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']

# Extract Hβ and MgII EW (our strongest pair)
hb_ew = d['HBETA_BR'][:, 2]  # EW at index 2
mg_ew = d['MGII_BR'][:, 2]
civ_ew = d['CIV'][:, 2]

# Valid mask
mask_hb_mg = ((hb_ew > 0) & (mg_ew > 0) & np.isfinite(hb_ew) & np.isfinite(mg_ew) &
              (z_all > 0.1) & (z_all < 2.5))

z_qso = z_all[mask_hb_mg]
hb = np.log10(hb_ew[mask_hb_mg])
mg = np.log10(mg_ew[mask_hb_mg])

print(f"Quasars (Hβ-MgII): {len(z_qso)}")

# Also load Pantheon+ SNe
import csv
sn_file = 'data/pantheon_plus.csv'
if not os.path.exists(sn_file):
    sn_file = 'results/pantheon_with_y.csv'

sn_z, sn_c, sn_x1, sn_mu, sn_mB = [], [], [], [], []
if os.path.exists(sn_file):
    with open(sn_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                zz = float(row.get('zHD', row.get('z', 0)))
                cc = float(row.get('c', 0))
                xx = float(row.get('x1', 0))
                mm = float(row.get('mB', row.get('m_b_corr', 0)))
                if zz > 0.01:
                    sn_z.append(zz)
                    sn_c.append(cc)
                    sn_x1.append(xx)
                    sn_mB.append(mm)
            except:
                pass
    sn_z = np.array(sn_z)
    sn_c = np.array(sn_c)
    sn_x1 = np.array(sn_x1)
    sn_mB = np.array(sn_mB)
    print(f"SNe Ia: {len(sn_z)}")

# FRB data
frb_file = 'data/chimefrbcat1.csv'
frb_dm, frb_width, frb_si = [], [], []
if os.path.exists(frb_file):
    with open(frb_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                dm = float(row.get('bonsai_dm', row.get('dm', 0)))
                w = float(row.get('width_fitb', row.get('width', 0)))
                si = float(row.get('spectral_index', 0))
                if dm > 0 and w > 0 and si != 0:
                    frb_dm.append(dm)
                    frb_width.append(w)
                    frb_si.append(si)
            except:
                pass
    frb_dm = np.array(frb_dm)
    frb_width = np.array(frb_width)
    frb_si = np.array(frb_si)
    print(f"FRBs: {len(frb_dm)}")

# ============================================================
# TEST 1: Critical Fluctuations — Does variance peak at z₀?
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: CRITICAL FLUCTUATIONS — Variance peak at phase boundary?")
print("=" * 70)
print("In phase transitions, fluctuations diverge at the critical point.")
print("Ice: low fluctuation. Water: low fluctuation. Melting point: MAXIMUM.\n")

# Quasars: compute local correlation + variance in sliding windows
window_size = 0.15
z_centers = np.arange(0.3, 2.1, 0.05)

qso_corrs = []
qso_vars = []  # variance of residuals from the local correlation
qso_scatter = []  # bootstrap variance of the correlation itself
qso_ns = []

for zc in z_centers:
    zmask = (z_qso >= zc - window_size) & (z_qso < zc + window_size)
    n = zmask.sum()
    if n > 50:
        r, p = stats.spearmanr(hb[zmask], mg[zmask])
        qso_corrs.append(r)
        qso_ns.append(n)
        
        # Bootstrap variance of the correlation
        boot_rs = []
        for _ in range(200):
            idx = np.random.choice(n, n, replace=True)
            br, _ = stats.spearmanr(hb[zmask][idx], mg[zmask][idx])
            boot_rs.append(br)
        qso_vars.append(np.var(boot_rs))
        qso_scatter.append(np.std(boot_rs))
    else:
        qso_corrs.append(np.nan)
        qso_vars.append(np.nan)
        qso_scatter.append(np.nan)
        qso_ns.append(0)

qso_corrs = np.array(qso_corrs)
qso_vars = np.array(qso_vars)
qso_scatter = np.array(qso_scatter)
valid = np.isfinite(qso_corrs)

# Find where variance peaks
if valid.sum() > 5:
    peak_idx = np.nanargmax(qso_vars)
    z_peak_var = z_centers[peak_idx]
    
    print(f"  Quasar Hβ-MgII EW correlation:")
    print(f"    Correlation range: {np.nanmin(qso_corrs):.3f} to {np.nanmax(qso_corrs):.3f}")
    print(f"    ★ VARIANCE PEAKS at z = {z_peak_var:.2f} (bootstrap σ² = {qso_vars[peak_idx]:.6f})")
    print(f"    Known sigmoid z₀ ≈ 1.05")
    print(f"    Peak variance / z₀ distance: {abs(z_peak_var - 1.05):.2f}")
    
    # Is variance distribution peaked or flat?
    mean_var = np.nanmean(qso_vars)
    max_var = np.nanmax(qso_vars)
    print(f"    Peak/mean variance ratio: {max_var/mean_var:.2f}x")
    
    print(f"\n    z-by-z variance profile:")
    for i, zc in enumerate(z_centers):
        if valid[i]:
            bar = "█" * int(qso_vars[i] / max_var * 40)
            marker = " ◄ PEAK" if i == peak_idx else ""
            marker2 = " ◄ z₀" if abs(zc - 1.05) < 0.03 else ""
            print(f"      z={zc:.2f}: r={qso_corrs[i]:+.3f} σ²={qso_vars[i]:.6f} {bar}{marker}{marker2}")

# ============================================================
# TEST 2: Critical Exponent — Power law near transition
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: CRITICAL EXPONENT — Power law near z₀")
print("=" * 70)
print("Near a phase transition: order_parameter ~ |T - Tc|^β")
print("For us: correlation ~ |z - z₀|^β near the sigmoid midpoint\n")

z0_qso = 1.05  # Known threshold

# Use the correlation values near the transition
# Below z₀ (ordered phase): correlation is positive and drops
# Above z₀ (disordered phase): correlation is near zero

# Extract the ordered side (z < z₀)
ordered_mask = valid & (z_centers < z0_qso) & (z_centers > 0.4)
disordered_mask = valid & (z_centers > z0_qso) & (z_centers < 2.0)

if ordered_mask.sum() > 3:
    z_ord = z_centers[ordered_mask]
    r_ord = qso_corrs[ordered_mask]
    
    # The "order parameter" = correlation value
    # Near z₀: r ∝ (z₀ - z)^β
    dz = z0_qso - z_ord
    
    # Only use points close to transition
    near = dz < 0.5
    if near.sum() > 3:
        # Log-log fit for power law
        log_dz = np.log10(dz[near])
        log_r = np.log10(np.abs(r_ord[near]) + 1e-10)
        
        slope, intercept, r_fit, p_fit, _ = stats.linregress(log_dz, log_r)
        beta = slope  # critical exponent
        
        print(f"  Ordered side (z < z₀):")
        print(f"    Power law fit: r ∝ |z₀ - z|^{beta:.3f}")
        print(f"    R² = {r_fit**2:.3f}, p = {p_fit:.3e}")
        
        # Known critical exponents for comparison
        print(f"\n  Known critical exponents (β):")
        print(f"    Mean field:  β = 0.500")
        print(f"    3D Ising:    β = 0.326")
        print(f"    3D XY:       β = 0.348")
        print(f"    3D Heisenberg: β = 0.365")
        print(f"    Percolation: β = 0.418")
        print(f"    ★ OUR β = {beta:.3f}")
        
        # Check if α = 1.845 relates to known exponents
        print(f"\n  Cooperative exponent α = 1.845 comparison:")
        print(f"    Susceptibility exponent γ (mean field): 1.000")
        print(f"    Susceptibility exponent γ (3D Ising):   1.237")
        print(f"    Susceptibility exponent γ (percolation): 1.793")
        print(f"    ★ Our α = 1.845 (closest to percolation γ!)")

# ============================================================
# TEST 3: Universality — Same exponents across source classes?
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: UNIVERSALITY — Same critical behavior across sources?")
print("=" * 70)
print("If this is a universal phase transition, different source classes")
print("should share the same critical exponents (like Ising universality).\n")

# For SNe: compute correlation(c, residual) vs z
if len(sn_z) > 100:
    # Simple Tripp residual
    # μ_obs - μ_model ≈ mB - 5*log10(dL) + 19.3
    # But we just need correlation of c with Hubble residual trend
    sn_window = 0.08
    sn_zc = np.arange(0.05, 1.5, 0.03)
    sn_corrs = []
    sn_vars_boot = []
    sn_valid_zc = []
    
    for zc in sn_zc:
        zm = (sn_z >= zc - sn_window) & (sn_z < zc + sn_window)
        if zm.sum() > 30:
            r, p = stats.spearmanr(sn_c[zm], sn_mB[zm])
            sn_corrs.append(r)
            sn_valid_zc.append(zc)
            
            boot_rs = []
            for _ in range(200):
                idx = np.random.choice(zm.sum(), zm.sum(), replace=True)
                br, _ = stats.spearmanr(sn_c[zm][idx], sn_mB[zm][idx])
                boot_rs.append(br)
            sn_vars_boot.append(np.var(boot_rs))
    
    sn_corrs = np.array(sn_corrs)
    sn_vars_boot = np.array(sn_vars_boot)
    sn_valid_zc = np.array(sn_valid_zc)
    
    if len(sn_corrs) > 3:
        peak_sn = np.argmax(sn_vars_boot)
        z_peak_sn = sn_valid_zc[peak_sn]
        
        print(f"  SNe Ia (color-mB correlation):")
        print(f"    Variance peaks at z = {z_peak_sn:.2f}")
        print(f"    Known sigmoid z₀ ≈ 0.82")
        print(f"    Distance from z₀: {abs(z_peak_sn - 0.82):.2f}")
        
        # Critical exponent for SNe
        z0_sn = 0.82
        ord_sn = (sn_valid_zc < z0_sn) & (sn_valid_zc > 0.1)
        if ord_sn.sum() > 3:
            dz_sn = z0_sn - sn_valid_zc[ord_sn]
            near_sn = dz_sn < 0.5
            if near_sn.sum() > 3:
                log_dz_sn = np.log10(dz_sn[near_sn])
                log_r_sn = np.log10(np.abs(sn_corrs[ord_sn][near_sn]) + 1e-10)
                slope_sn, _, r_sn, p_sn, _ = stats.linregress(log_dz_sn, log_r_sn)
                print(f"    Critical exponent β_SN = {slope_sn:.3f} (R²={r_sn**2:.3f})")

# FRBs
if len(frb_dm) > 50:
    dm_bins = np.arange(100, 1500, 100)
    frb_corrs = []
    frb_vars_boot = []
    frb_valid_dm = []
    
    for i in range(len(dm_bins) - 1):
        dm_m = (frb_dm >= dm_bins[i]) & (frb_dm < dm_bins[i+1])
        if dm_m.sum() > 15:
            r, p = stats.spearmanr(frb_width[dm_m], frb_si[dm_m])
            frb_corrs.append(r)
            frb_valid_dm.append((dm_bins[i] + dm_bins[i+1]) / 2)
            
            boot_rs = []
            for _ in range(200):
                idx = np.random.choice(dm_m.sum(), dm_m.sum(), replace=True)
                br, _ = stats.spearmanr(frb_width[dm_m][idx], frb_si[dm_m][idx])
                boot_rs.append(br)
            frb_vars_boot.append(np.var(boot_rs))
    
    if len(frb_corrs) > 2:
        frb_corrs = np.array(frb_corrs)
        frb_vars_boot = np.array(frb_vars_boot)
        peak_frb = np.argmax(frb_vars_boot)
        
        print(f"\n  FRBs (width-spectral_index):")
        print(f"    Variance peaks at DM = {frb_valid_dm[peak_frb]:.0f}")
        print(f"    Known threshold DM ≈ 500")
        print(f"    DM bins: {[f'{d:.0f}' for d in frb_valid_dm]}")
        print(f"    Correlations: {[f'{c:.3f}' for c in frb_corrs]}")
        print(f"    Variances: {[f'{v:.5f}' for v in frb_vars_boot]}")

# ============================================================
# TEST 4: Order parameter identification
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: ORDER PARAMETER — What's the 'magnetization'?")
print("=" * 70)
print("In a magnet, order = alignment. In ice, order = lattice position.")
print("In our data, order = ...?\n")

# Candidates for order parameter:
# 1. The correlation itself (r between observables)
# 2. PC1 variance fraction (how much variance is in the leading mode)
# 3. Effective dimensionality (how many PCs needed to explain 90%)
# 4. Mutual information between observables

# Let's compute all four as function of z for quasars
from numpy.linalg import eigh

print("  Computing order parameters vs z for quasars...\n")

op_z = []
op_corr = []  # raw correlation
op_pc1frac = []  # PC1 variance fraction
op_eff_dim = []  # effective dimensionality
op_mi_proxy = []  # MI proxy (correlation-based)

z_op_bins = np.arange(0.3, 2.0, 0.1)

# Need multiple columns for PCA — use Hβ EW, MgII EW, Hβ FWHM, MgII FWHM
hb_fwhm = d['HBETA_BR'][:, 4]
mg_fwhm = d['MGII_BR'][:, 4]

full_mask = (mask_hb_mg & (hb_fwhm > 0) & (mg_fwhm > 0) & 
             np.isfinite(hb_fwhm) & np.isfinite(mg_fwhm))

z_full = z_all[full_mask]
X_full = np.column_stack([
    np.log10(hb_ew[full_mask]),
    np.log10(mg_ew[full_mask]),
    np.log10(hb_fwhm[full_mask]),
    np.log10(mg_fwhm[full_mask])
])
finite_full = np.all(np.isfinite(X_full), axis=1)
z_full = z_full[finite_full]
X_full = X_full[finite_full]

for i in range(len(z_op_bins) - 1):
    zm = (z_full >= z_op_bins[i]) & (z_full < z_op_bins[i+1])
    if zm.sum() > 100:
        zc = (z_op_bins[i] + z_op_bins[i+1]) / 2
        op_z.append(zc)
        
        X_bin = X_full[zm]
        X_std = (X_bin - X_bin.mean(axis=0)) / X_bin.std(axis=0)
        
        # 1. Correlation
        r, _ = stats.spearmanr(X_bin[:, 0], X_bin[:, 1])  # Hβ-MgII EW
        op_corr.append(r)
        
        # 2. PC1 variance fraction
        cov = np.cov(X_std.T)
        eigenvalues = np.sort(np.linalg.eigvalsh(cov))[::-1]
        pc1_frac = eigenvalues[0] / eigenvalues.sum()
        op_pc1frac.append(pc1_frac)
        
        # 3. Effective dimensionality (participation ratio)
        p = eigenvalues / eigenvalues.sum()
        eff_dim = 1.0 / np.sum(p**2)  # participation ratio
        op_eff_dim.append(eff_dim)
        
        # 4. MI proxy: -0.5 * log(1 - r^2)
        mi = -0.5 * np.log(1 - r**2) if abs(r) < 0.999 else 5.0
        op_mi_proxy.append(mi)

op_z = np.array(op_z)
op_corr = np.array(op_corr)
op_pc1frac = np.array(op_pc1frac)
op_eff_dim = np.array(op_eff_dim)
op_mi_proxy = np.array(op_mi_proxy)

print(f"  {'z':>5} {'corr':>8} {'PC1%':>8} {'eff_dim':>8} {'MI':>8}")
print(f"  {'-'*5} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
for i in range(len(op_z)):
    lattice = "█" if op_corr[i] > 0.5 else "▓" if op_corr[i] > 0.2 else "░" if op_corr[i] > 0.05 else " "
    print(f"  {op_z[i]:5.2f} {op_corr[i]:+8.3f} {op_pc1frac[i]:8.1%} {op_eff_dim[i]:8.2f} {op_mi_proxy[i]:8.3f}  {lattice}")

# Which order parameter has the sharpest transition?
for name, vals in [('correlation', op_corr), ('PC1 fraction', op_pc1frac), 
                    ('eff dimension', op_eff_dim), ('MI proxy', op_mi_proxy)]:
    if len(vals) > 4:
        # Find max gradient
        grad = np.gradient(vals, op_z)
        max_grad_idx = np.argmax(np.abs(grad))
        z_sharpest = op_z[max_grad_idx]
        print(f"\n  {name}: sharpest change at z = {z_sharpest:.2f} (gradient = {grad[max_grad_idx]:.3f})")

# ============================================================
# TEST 5: The lattice — correlation BETWEEN z-bins
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: LATTICE COHERENCE — Correlation between z-bins")
print("=" * 70)
print("In a crystal, distant atoms still correlate (long-range order).")
print("In a liquid, only nearby atoms correlate (short-range).")
print("How far does 'order' extend in z-space?\n")

# Compute the covariance structure at each z, then compare
# covariance structures between z-bins
# If structures are similar across z → long-range lattice order
# If they diverge with z-separation → short-range only

n_zbins = len(op_z)
if n_zbins > 4:
    # Recompute covariance matrices per z-bin
    cov_matrices = []
    for i in range(len(z_op_bins) - 1):
        zm = (z_full >= z_op_bins[i]) & (z_full < z_op_bins[i+1])
        if zm.sum() > 100:
            X_bin = X_full[zm]
            X_std = (X_bin - X_bin.mean(axis=0)) / X_bin.std(axis=0)
            cov_matrices.append(np.cov(X_std.T))
    
    # Compute similarity between all pairs of z-bins
    # Use Frobenius distance between covariance matrices
    n_cov = len(cov_matrices)
    print(f"  Covariance similarity matrix ({n_cov} z-bins):")
    print(f"  (Frobenius distance — smaller = more similar)\n")
    
    dist_matrix = np.zeros((n_cov, n_cov))
    for i in range(n_cov):
        for j in range(n_cov):
            dist_matrix[i, j] = np.linalg.norm(cov_matrices[i] - cov_matrices[j], 'fro')
    
    # Does distance grow with z-separation?
    pairs_dz = []
    pairs_dist = []
    for i in range(n_cov):
        for j in range(i+1, n_cov):
            pairs_dz.append(abs(op_z[i] - op_z[j]) if i < len(op_z) and j < len(op_z) else j - i)
            pairs_dist.append(dist_matrix[i, j])
    
    r_lattice, p_lattice = stats.spearmanr(pairs_dz, pairs_dist)
    print(f"  ρ(Δz, cov_distance) = {r_lattice:.3f} (p = {p_lattice:.3e})")
    print(f"  {'Long-range order BREAKS with distance' if r_lattice > 0.3 else 'Order persists across z (lattice intact)'}")
    
    # Nearby vs distant similarity
    near_dists = [d for dz, d in zip(pairs_dz, pairs_dist) if dz < 0.3]
    far_dists = [d for dz, d in zip(pairs_dz, pairs_dist) if dz > 0.5]
    if near_dists and far_dists:
        print(f"\n  Nearby z-bins (Δz < 0.3): mean cov distance = {np.mean(near_dists):.3f}")
        print(f"  Distant z-bins (Δz > 0.5): mean cov distance = {np.mean(far_dists):.3f}")
        print(f"  Ratio: {np.mean(far_dists)/np.mean(near_dists):.2f}x")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("PHASE TRANSITION SUMMARY")
print("=" * 70)

all_results = {
    'critical_fluctuations': {
        'quasar_variance_peak_z': float(z_peak_var) if 'z_peak_var' in dir() else None,
        'known_z0': 1.05,
    },
    'order_parameters': {
        'z': op_z.tolist(),
        'correlation': op_corr.tolist(),
        'pc1_fraction': op_pc1frac.tolist(),
        'effective_dim': op_eff_dim.tolist(),
        'mi_proxy': op_mi_proxy.tolist(),
    },
    'lattice_coherence': {
        'rho_dz_covdist': float(r_lattice) if 'r_lattice' in dir() else None,
    }
}

with open('results_phase_transition/phase_transition_results.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print("\nResults saved to results_phase_transition/")
