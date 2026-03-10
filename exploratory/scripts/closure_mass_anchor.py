#!/usr/bin/env python3
"""
MASS ANCHOR TESTS — Paper 2 Foundation
Does mass modulate the sigmoid? Does mass shift the threshold?

If mass moves the sigmoid: geometry alone can't explain it.
In GR, geometry IS mass. If geometry selectively degrades mass-dependent
information, it's acting on itself — demands mechanism, not just projection.

Tests:
1. BH mass vs sigmoid threshold — split quasars by mass, fit sigmoids separately
2. BH mass vs degradation rate — at fixed z, massive vs light BHs
3. Eddington ratio (thermodynamic temperature) vs degradation
4. Luminosity at fixed mass — pure energy vs mass
5. Mass-z interaction — does mass change the SLOPE?
6. The angle — PCA rotation in quasar parameter space
"""

import numpy as np
from astropy.io import fits
from scipy import stats
from scipy.optimize import curve_fit
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_mass_anchor', exist_ok=True)

# ============================================================
# Load DR16Q
# ============================================================
print("=" * 70)
print("MASS ANCHOR TESTS — Does Mass Move the Sigmoid?")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data

z_all = d['Z_DR16Q']
logmbh = d['LOGMBH']
loglbol = d['LOGLBOL']
logledd = d['LOGLEDD_RATIO']
sn_median = d['SN_MEDIAN_ALL']

# Extract line measurements
# Structure: [peak, centroid, EW, logL, FWHM, shift]
def extract_line(data, colname_pref):
    """Try _BR first, then total."""
    for suffix in [f'{colname_pref}_BR', colname_pref]:
        if suffix in data.dtype.names:
            arr = data[suffix]
            if arr.ndim == 2 and arr.shape[1] >= 5:
                ew = arr[:, 2]
                fwhm = arr[:, 4]
                logl = arr[:, 3]
                return ew, fwhm, logl
            elif arr.ndim == 1:
                # Structured differently - try accessing
                try:
                    ew = np.array([x[2] if len(x) >= 3 else 0 for x in arr])
                    fwhm = np.array([x[4] if len(x) >= 5 else 0 for x in arr])
                    logl = np.array([x[3] if len(x) >= 4 else 0 for x in arr])
                    return ew, fwhm, logl
                except:
                    pass
    return None, None, None

# Extract key lines
hb_ew, hb_fwhm, hb_logl = extract_line(d, 'HBETA')
mg_ew, mg_fwhm, mg_logl = extract_line(d, 'MGII')
civ_ew, civ_fwhm, civ_logl = extract_line(d, 'CIV')

print(f"Total quasars: {len(z_all)}")
print(f"Hβ EW valid: {np.sum((hb_ew is not None) & (hb_ew > 0)) if hb_ew is not None else 0}")
print(f"MgII EW valid: {np.sum((mg_ew is not None) & (mg_ew > 0)) if mg_ew is not None else 0}")
print(f"CIV EW valid: {np.sum((civ_ew is not None) & (civ_ew > 0)) if civ_ew is not None else 0}")
print(f"LOGMBH valid: {np.sum(logmbh > 0)}")
print(f"LOGLBOL valid: {np.sum(loglbol > 0)}")
print(f"LOGLEDD valid: {np.sum(np.isfinite(logledd) & (logledd != 0))}")

# Sigmoid function
def sigmoid(z, A, z0, k, B):
    return A / (1 + np.exp(-k * (z - z0))) + B

# ============================================================
# TEST 1: BH Mass Split — Does mass shift the sigmoid?
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: BH MASS SPLIT — Does mass shift the sigmoid threshold?")
print("=" * 70)

# Use Hβ-MgII EW correlation (our strongest quasar pair)
# Need both lines measured + valid BH mass
mask = ((hb_ew > 0) & (mg_ew > 0) & (logmbh > 0) & 
        np.isfinite(hb_ew) & np.isfinite(mg_ew) & np.isfinite(logmbh) &
        (z_all > 0.1) & (z_all < 2.5))

z_use = z_all[mask]
hb_use = np.log10(hb_ew[mask])  # Log space for EWs
mg_use = np.log10(mg_ew[mask])
mbh_use = logmbh[mask]
lbol_use = loglbol[mask]
ledd_use = logledd[mask]

print(f"\nUsable quasars (Hβ + MgII + BH mass): {np.sum(mask)}")
print(f"BH mass range: {mbh_use.min():.1f} to {mbh_use.max():.1f} (log M☉)")
print(f"Median BH mass: {np.median(mbh_use):.2f}")
print(f"z range: {z_use.min():.2f} to {z_use.max():.2f}")

# Split into mass terciles
mass_lo = np.percentile(mbh_use, 33)
mass_hi = np.percentile(mbh_use, 67)
print(f"\nMass terciles: < {mass_lo:.2f} | {mass_lo:.2f}-{mass_hi:.2f} | > {mass_hi:.2f}")

mass_labels = ['LOW mass', 'MID mass', 'HIGH mass']
mass_masks = [
    mbh_use < mass_lo,
    (mbh_use >= mass_lo) & (mbh_use <= mass_hi),
    mbh_use > mass_hi
]

# For each mass bin, compute rolling EW correlation vs z
z_bins = np.arange(0.3, 2.2, 0.15)
results_t1 = {}

for mi, (mlabel, mmask) in enumerate(zip(mass_labels, mass_masks)):
    z_m = z_use[mmask]
    hb_m = hb_use[mmask]
    mg_m = mg_use[mmask]
    
    corrs = []
    z_centers = []
    ns = []
    
    for i in range(len(z_bins) - 1):
        zmask = (z_m >= z_bins[i]) & (z_m < z_bins[i+1])
        if zmask.sum() > 30:
            r, p = stats.spearmanr(hb_m[zmask], mg_m[zmask])
            corrs.append(r)
            z_centers.append((z_bins[i] + z_bins[i+1]) / 2)
            ns.append(zmask.sum())
    
    corrs = np.array(corrs)
    z_centers = np.array(z_centers)
    
    # Fit sigmoid
    try:
        popt, _ = curve_fit(sigmoid, z_centers, corrs, 
                           p0=[-0.5, 1.0, 5.0, 0.5], maxfev=10000,
                           bounds=([-2, 0.1, 0.1, -1], [0, 3.0, 50, 1]))
        z0_fit = popt[1]
        k_fit = popt[2]
    except:
        z0_fit = np.nan
        k_fit = np.nan
    
    # Linear trend
    slope, intercept, r_trend, p_trend, _ = stats.linregress(z_centers, corrs)
    
    results_t1[mlabel] = {
        'N': int(mmask.sum()),
        'mean_logMBH': float(np.mean(mbh_use[mmask])),
        'z_centers': z_centers.tolist(),
        'correlations': corrs.tolist(),
        'n_per_bin': [int(n) for n in ns],
        'sigmoid_z0': float(z0_fit),
        'sigmoid_k': float(k_fit),
        'linear_slope': float(slope),
        'linear_r': float(r_trend),
        'linear_p': float(p_trend),
        'corr_at_low_z': float(np.mean(corrs[:3])) if len(corrs) >= 3 else None,
        'corr_at_high_z': float(np.mean(corrs[-3:])) if len(corrs) >= 3 else None,
    }
    
    print(f"\n  {mlabel} (N={mmask.sum()}, <logM>={np.mean(mbh_use[mmask]):.2f}):")
    print(f"    Low-z corr: {results_t1[mlabel]['corr_at_low_z']:.3f}")
    print(f"    High-z corr: {results_t1[mlabel]['corr_at_high_z']:.3f}")
    print(f"    Sigmoid z₀ = {z0_fit:.3f}, k = {k_fit:.2f}")
    print(f"    Linear slope = {slope:.4f} (r={r_trend:.3f}, p={p_trend:.2e})")

# Compare sigmoid thresholds across mass
z0_values = [results_t1[m]['sigmoid_z0'] for m in mass_labels]
mass_values = [results_t1[m]['mean_logMBH'] for m in mass_labels]
if all(np.isfinite(z0_values)):
    r_mass_z0, p_mass_z0 = stats.spearmanr(mass_values, z0_values)
    print(f"\n  ★ MASS vs SIGMOID THRESHOLD: ρ = {r_mass_z0:.3f} (p = {p_mass_z0:.3e})")
    print(f"    z₀ values: {[f'{v:.3f}' for v in z0_values]}")
    print(f"    Mass shifts sigmoid: {'YES' if p_mass_z0 < 0.05 else 'NO'}")

# ============================================================
# TEST 2: At fixed z, does mass modulate correlation strength?
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: MASS vs CORRELATION AT FIXED REDSHIFT")
print("=" * 70)

results_t2 = {}
z_slices = [(0.3, 0.6), (0.6, 0.9), (0.9, 1.2), (1.2, 1.6), (1.6, 2.2)]

for zlo, zhi in z_slices:
    zmask = (z_use >= zlo) & (z_use < zhi)
    if zmask.sum() < 100:
        continue
    
    z_s = z_use[zmask]
    hb_s = hb_use[zmask]
    mg_s = mg_use[zmask]
    mbh_s = mbh_use[zmask]
    
    # Split into mass quintiles within this z-slice
    mass_quints = np.percentile(mbh_s, [20, 40, 60, 80])
    quint_masks = [
        mbh_s <= mass_quints[0],
        (mbh_s > mass_quints[0]) & (mbh_s <= mass_quints[1]),
        (mbh_s > mass_quints[1]) & (mbh_s <= mass_quints[2]),
        (mbh_s > mass_quints[2]) & (mbh_s <= mass_quints[3]),
        mbh_s > mass_quints[3]
    ]
    
    corrs = []
    masses = []
    for qi, qm in enumerate(quint_masks):
        if qm.sum() > 20:
            r, p = stats.spearmanr(hb_s[qm], mg_s[qm])
            corrs.append(r)
            masses.append(np.median(mbh_s[qm]))
    
    if len(corrs) >= 3:
        r_mq, p_mq = stats.spearmanr(masses, corrs)
        zkey = f"z=[{zlo},{zhi})"
        results_t2[zkey] = {
            'N': int(zmask.sum()),
            'mass_quintile_medians': [float(m) for m in masses],
            'correlations': [float(c) for c in corrs],
            'spearman_r': float(r_mq),
            'spearman_p': float(p_mq)
        }
        direction = "MASSIVE PRESERVE" if r_mq > 0 else "MASSIVE DEGRADE"
        sig = "★" if p_mq < 0.05 else " "
        print(f"  {sig} z=[{zlo:.1f},{zhi:.1f}): ρ(mass,corr) = {r_mq:+.3f} (p={p_mq:.3e}) N={zmask.sum()} — {direction}")
        for qi, (m, c) in enumerate(zip(masses, corrs)):
            print(f"      Q{qi+1}: logM={m:.2f}, r(Hβ,MgII)={c:.3f}")

# ============================================================
# TEST 3: Eddington Ratio (Temperature Proxy) vs Degradation
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: EDDINGTON RATIO (TEMPERATURE) vs DEGRADATION")
print("=" * 70)
print("Eddington ratio = L/L_Edd — how hard the BH is accreting")
print("High Edd = hot, radiatively efficient. Low Edd = cold, inefficient.\n")

edd_valid = np.isfinite(ledd_use) & (ledd_use != 0) & (ledd_use > -5)
results_t3 = {}

for zlo, zhi in z_slices:
    zmask = (z_use >= zlo) & (z_use < zhi) & edd_valid
    if zmask.sum() < 100:
        continue
    
    hb_s = hb_use[zmask]
    mg_s = mg_use[zmask]
    ledd_s = ledd_use[zmask]
    mbh_s = mbh_use[zmask]
    
    # Split by Eddington ratio terciles
    edd_lo, edd_hi = np.percentile(ledd_s, [33, 67])
    
    edd_masks = [ledd_s < edd_lo, (ledd_s >= edd_lo) & (ledd_s <= edd_hi), ledd_s > edd_hi]
    edd_labels = ['COLD', 'WARM', 'HOT']
    
    corrs_edd = []
    edd_vals = []
    for el, em in zip(edd_labels, edd_masks):
        if em.sum() > 20:
            r, p = stats.spearmanr(hb_s[em], mg_s[em])
            corrs_edd.append(r)
            edd_vals.append(np.median(ledd_s[em]))
    
    if len(corrs_edd) >= 3:
        r_eq, p_eq = stats.spearmanr(edd_vals, corrs_edd)
        zkey = f"z=[{zlo},{zhi})"
        results_t3[zkey] = {
            'N': int(zmask.sum()),
            'edd_medians': [float(e) for e in edd_vals],
            'correlations': [float(c) for c in corrs_edd],
            'spearman_r': float(r_eq),
            'spearman_p': float(p_eq)
        }
        direction = "HOT PRESERVES" if r_eq > 0 else "HOT DEGRADES"
        sig = "★" if p_eq < 0.05 else " "
        print(f"  {sig} z=[{zlo:.1f},{zhi:.1f}): ρ(Edd,corr) = {r_eq:+.3f} (p={p_eq:.3e}) — {direction}")

# ============================================================
# TEST 4: Mass-z INTERACTION — Does mass change the SLOPE?
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: MASS-z INTERACTION — Does mass change degradation slope?")
print("=" * 70)
print("If mass shifts the sigmoid, the degradation slope should differ by mass.\n")

results_t4 = {}
for mi, (mlabel, mmask) in enumerate(zip(mass_labels, mass_masks)):
    z_m = z_use[mmask]
    hb_m = hb_use[mmask]
    mg_m = mg_use[mmask]
    
    # Rolling correlation with smaller bins for slope measurement
    z_fine = np.arange(0.3, 2.0, 0.1)
    corrs = []
    zc = []
    for i in range(len(z_fine) - 2):
        zm = (z_m >= z_fine[i]) & (z_m < z_fine[i+2])
        if zm.sum() > 50:
            r, p = stats.spearmanr(hb_m[zm], mg_m[zm])
            corrs.append(r)
            zc.append((z_fine[i] + z_fine[i+2]) / 2)
    
    if len(corrs) >= 4:
        slope, _, r_val, p_val, _ = stats.linregress(zc, corrs)
        results_t4[mlabel] = {
            'degradation_slope': float(slope),
            'slope_r': float(r_val),
            'slope_p': float(p_val)
        }
        print(f"  {mlabel}: degradation slope = {slope:.4f}/z (r={r_val:.3f}, p={p_val:.2e})")

slopes = [results_t4[m]['degradation_slope'] for m in mass_labels if m in results_t4]
if len(slopes) == 3:
    ratio = slopes[2] / slopes[0] if slopes[0] != 0 else np.nan
    print(f"\n  Slope ratio (HIGH/LOW mass): {ratio:.3f}")
    print(f"  If ≈ 1.0: mass doesn't change slope. If ≠ 1.0: mass modulates degradation rate.")

# ============================================================
# TEST 5: Luminosity at FIXED mass — Energy vs Mass
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: LUMINOSITY AT FIXED MASS — Pure energy effect")
print("=" * 70)
print("Control mass, vary luminosity. Does pure energy affect degradation?\n")

results_t5 = {}
# Fix mass in a narrow band, split by luminosity
mass_band = (mbh_use >= 8.5) & (mbh_use <= 9.5)  # ~1 dex mass window

for zlo, zhi in [(0.4, 0.8), (0.8, 1.3), (1.3, 2.0)]:
    zmask = (z_use >= zlo) & (z_use < zhi) & mass_band & (loglbol[mask][np.arange(len(z_use))] > 0)
    # Need to handle loglbol indexing properly
    lbol_valid = np.isfinite(lbol_use) & (lbol_use > 0)
    zmask = (z_use >= zlo) & (z_use < zhi) & mass_band & lbol_valid
    
    if zmask.sum() < 100:
        continue
    
    hb_s = hb_use[zmask]
    mg_s = mg_use[zmask]
    lbol_s = lbol_use[zmask]
    
    lbol_med = np.median(lbol_s)
    lo_mask = lbol_s < lbol_med
    hi_mask = lbol_s >= lbol_med
    
    r_lo, p_lo = stats.spearmanr(hb_s[lo_mask], mg_s[lo_mask]) if lo_mask.sum() > 20 else (np.nan, np.nan)
    r_hi, p_hi = stats.spearmanr(hb_s[hi_mask], mg_s[hi_mask]) if hi_mask.sum() > 20 else (np.nan, np.nan)
    
    delta = r_hi - r_lo
    zkey = f"z=[{zlo},{zhi})"
    results_t5[zkey] = {
        'N': int(zmask.sum()),
        'mass_range': '8.5-9.5',
        'r_low_lum': float(r_lo),
        'r_high_lum': float(r_hi),
        'delta_r': float(delta)
    }
    direction = "BRIGHT PRESERVES" if delta > 0 else "BRIGHT DEGRADES"
    print(f"  z=[{zlo},{zhi}): r_dim={r_lo:.3f}, r_bright={r_hi:.3f}, Δ={delta:+.3f} — {direction}")

# ============================================================
# TEST 6: PCA ROTATION IN QUASAR SPACE — The angle
# ============================================================
print("\n" + "=" * 70)
print("TEST 6: PCA ROTATION IN QUASAR PARAMETER SPACE")
print("=" * 70)
print("Does the principal component rotate with z in quasar space")
print("the same way it does in SN space (~20-25°)?\n")

# Use EW + FWHM of Hβ and MgII as 4D parameter space
# (analogous to x1, c, mB in SN space)
full_mask = ((hb_ew > 0) & (mg_ew > 0) & (hb_fwhm > 0) & (mg_fwhm > 0) &
             np.isfinite(hb_ew) & np.isfinite(mg_ew) & 
             np.isfinite(hb_fwhm) & np.isfinite(mg_fwhm) &
             (z_all > 0.1) & (z_all < 2.5))

z_pca = z_all[full_mask]
X_pca = np.column_stack([
    np.log10(hb_ew[full_mask]),
    np.log10(mg_ew[full_mask]),
    np.log10(hb_fwhm[full_mask]),
    np.log10(mg_fwhm[full_mask])
])

# Remove any remaining NaN/inf
finite_mask = np.all(np.isfinite(X_pca), axis=1)
z_pca = z_pca[finite_mask]
X_pca = X_pca[finite_mask]

print(f"Objects for PCA: {len(z_pca)}")

# Split into z-bins and compute PCA
z_pca_bins = [(0.3, 0.6), (0.6, 0.9), (0.9, 1.3), (1.3, 1.8), (1.8, 2.5)]
pc1_vectors = []
pc1_variances = []
results_t6 = {}

from numpy.linalg import eigh

for zlo, zhi in z_pca_bins:
    zm = (z_pca >= zlo) & (z_pca < zhi)
    if zm.sum() < 100:
        continue
    
    X_bin = X_pca[zm]
    X_std = (X_bin - X_bin.mean(axis=0)) / X_bin.std(axis=0)
    
    cov = np.cov(X_std.T)
    eigenvalues, eigenvectors = eigh(cov)
    
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    pc1 = eigenvectors[:, 0]
    var_explained = eigenvalues[0] / eigenvalues.sum()
    
    pc1_vectors.append(pc1)
    pc1_variances.append(var_explained)
    
    zkey = f"z=[{zlo},{zhi})"
    results_t6[zkey] = {
        'N': int(zm.sum()),
        'PC1': pc1.tolist(),
        'PC1_variance': float(var_explained),
        'eigenvalues': eigenvalues.tolist()
    }
    
    print(f"  {zkey} (N={zm.sum()}): PC1 var={var_explained:.1%}")
    print(f"    PC1 = [{pc1[0]:+.3f}, {pc1[1]:+.3f}, {pc1[2]:+.3f}, {pc1[3]:+.3f}]")
    print(f"           (Hβ_EW,  MgII_EW, Hβ_FWHM, MgII_FWHM)")

# Compute angles between consecutive PC1 vectors
if len(pc1_vectors) >= 2:
    print(f"\n  PC1 rotation between z-bins:")
    total_rotation = 0
    for i in range(len(pc1_vectors) - 1):
        cos_angle = np.abs(np.dot(pc1_vectors[i], pc1_vectors[i+1]))
        cos_angle = min(cos_angle, 1.0)  # numerical safety
        angle = np.degrees(np.arccos(cos_angle))
        total_rotation += angle
        z1 = z_pca_bins[i]
        z2 = z_pca_bins[i+1]
        print(f"    z=[{z1[0]},{z1[1]}) → z=[{z2[0]},{z2[1]}): {angle:.1f}°")
    
    # Total rotation low-z to high-z
    cos_total = np.abs(np.dot(pc1_vectors[0], pc1_vectors[-1]))
    cos_total = min(cos_total, 1.0)
    total_angle = np.degrees(np.arccos(cos_total))
    print(f"\n  ★ TOTAL PC1 ROTATION (lowest → highest z): {total_angle:.1f}°")
    print(f"    (SN space showed ~20-25° rotation)")
    
    # Variance drainage
    print(f"\n  PC1 variance explained:")
    for i, ((zlo, zhi), v) in enumerate(zip(z_pca_bins[:len(pc1_variances)], pc1_variances)):
        print(f"    z=[{zlo},{zhi}): {v:.1%}")
    
    drainage = pc1_variances[0] - pc1_variances[-1]
    print(f"    Drainage: {pc1_variances[0]:.1%} → {pc1_variances[-1]:.1%} (Δ = {drainage:+.1%})")

# ============================================================
# TEST 7: Box's M on quasar parameters — formal covariance test
# ============================================================
print("\n" + "=" * 70)
print("TEST 7: BOX'S M — Covariance non-stationarity in quasar space")
print("=" * 70)

def box_m_test(X_groups):
    """Simplified Box's M statistic."""
    p = X_groups[0].shape[1]
    k = len(X_groups)
    ns = [len(g) for g in X_groups]
    N = sum(ns)
    
    covs = [np.cov(g.T) for g in X_groups]
    S_pooled = sum((n-1)*c for n, c in zip(ns, covs)) / (N - k)
    
    M = 0
    for i, (n, c) in enumerate(zip(ns, covs)):
        sign, logdet_c = np.linalg.slogdet(c)
        sign_p, logdet_p = np.linalg.slogdet(S_pooled)
        if sign > 0 and sign_p > 0:
            M += (n - 1) * (logdet_p - logdet_c)
    
    # Correction factor
    sum_inv = sum(1/(n-1) for n in ns) - 1/(N-k)
    c1 = (2*p*p + 3*p - 1) / (6*(p+1)*(k-1)) * sum_inv
    
    df = p * (p + 1) * (k - 1) / 2
    chi2 = M * (1 - c1)
    p_val = 1 - stats.chi2.cdf(chi2, df)
    
    return chi2, df, p_val

# Split quasars into low-z and high-z groups
z_split = 1.0
lo_mask_pca = z_pca < z_split
hi_mask_pca = z_pca >= z_split

if lo_mask_pca.sum() > 100 and hi_mask_pca.sum() > 100:
    # Standardize together
    X_std_all = (X_pca - X_pca.mean(axis=0)) / X_pca.std(axis=0)
    
    X_lo = X_std_all[lo_mask_pca]
    X_hi = X_std_all[hi_mask_pca]
    
    # Subsample to avoid numerical issues with huge N
    np.random.seed(42)
    n_sub = min(10000, len(X_lo), len(X_hi))
    X_lo_sub = X_lo[np.random.choice(len(X_lo), n_sub, replace=False)]
    X_hi_sub = X_hi[np.random.choice(len(X_hi), n_sub, replace=False)]
    
    chi2, df, p_val = box_m_test([X_lo_sub, X_hi_sub])
    
    sigma = np.sqrt(2) * stats.norm.isf(p_val) if p_val > 0 else 99
    print(f"  Box's M (z < {z_split} vs z ≥ {z_split}): χ² = {chi2:.1f}, df = {df:.0f}, p = {p_val:.2e}")
    print(f"  Significance: {sigma:.1f}σ")
    print(f"  (SN space: p = 2.3×10⁻⁷ ≈ 5σ)")
    print(f"  N_low = {lo_mask_pca.sum()}, N_high = {hi_mask_pca.sum()} (subsampled to {n_sub})")
    
    results_t6['box_m'] = {
        'chi2': float(chi2),
        'df': float(df),
        'p_value': float(p_val),
        'sigma': float(sigma),
        'n_subsample': n_sub
    }

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY — Mass Anchor Results")
print("=" * 70)

all_results = {
    'test1_mass_sigmoid': results_t1,
    'test2_mass_fixed_z': results_t2,
    'test3_eddington': results_t3,
    'test4_mass_slope': results_t4,
    'test5_luminosity_fixed_mass': results_t5,
    'test6_pca_rotation': results_t6,
}

with open('results_mass_anchor/mass_anchor_results.json', 'w') as f:
    json.dump(all_results, f, indent=2, default=str)

print("\nResults saved to results_mass_anchor/mass_anchor_results.json")
print("\nKey questions answered:")
print("  1. Does BH mass shift the sigmoid threshold?")
print("  2. At fixed z, do massive BHs preserve correlations?")
print("  3. Does temperature (Eddington ratio) matter?")
print("  4. Does mass change the degradation slope?")
print("  5. Does luminosity matter at fixed mass?")
print("  6. Does quasar parameter space rotate like SN space?")
print("  7. Is quasar covariance formally non-stationary (Box's M)?")
