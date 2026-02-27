#!/usr/bin/env python3
"""
THE TUBE TEST — Do quasars sharing the same path share the same damage?
========================================================================
Find pairs that are:
- Close on the sky (< 0.5°, < 1°, < 2°)  
- Close in redshift (Δz < 0.05)
Their photons traveled through nearly the SAME tube of vacuum.

If damage residuals correlate within tubes → path matters → spatial structure exists
If they don't → the 84% is intrinsic source stochasticity

Also: SDSS plate analysis. Same plate = same instrument conditions.
If plate explains residual variance → instrumental.
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy import stats
from scipy.spatial import cKDTree
import os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_tube_test'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("THE TUBE TEST — Shared path = shared damage?")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

def fix(arr): return arr.astype(arr.dtype.newbyteorder('='))

z = fix(data['Z_DR16Q'])
ra = fix(data['RA'])
dec = fix(data['DEC'])
snr = fix(data['SN_MEDIAN_ALL'])
oiii_flux = fix(data['OIII5007'][:, 4])
hbeta_flux = fix(data['HBETA'][:, 4])
hbeta_fwhm = fix(data['HBETA'][:, 3])
logLbol = fix(data['LOGLBOL'])
logMbh = fix(data['LOGMBH'])
logLedd = fix(data['LOGLEDD_RATIO'])

# SDSS plate info
try:
    plate = fix(data['PLATE'])
    mjd = fix(data['MJD'])
    fiberid = fix(data['FIBERID'])
    HAS_PLATE = True
except:
    HAS_PLATE = False
    print("  No plate info found")

good = ((z > 0.3) & (z < 2.0) & (snr > 3) &
        (oiii_flux > 0) & (hbeta_flux > 0) &
        np.isfinite(oiii_flux) & np.isfinite(hbeta_flux) &
        np.isfinite(logLbol) & np.isfinite(logMbh) & np.isfinite(logLedd) &
        (logLbol > 0) & (logMbh > 0) &
        np.isfinite(hbeta_fwhm) & (hbeta_fwhm > 0) &
        np.isfinite(ra) & np.isfinite(dec))

idx = np.where(good)[0]
N = len(idx)
Z = z[idx]; RA = ra[idx]; DEC = dec[idx]; SNR_arr = snr[idx]
log_ratio = np.log10(oiii_flux[idx] / hbeta_flux[idx])
HB_FWHM = np.log10(np.maximum(hbeta_fwhm[idx], 1))
LEDD = logLedd[idx]; LBOL = logLbol[idx]; MBH = logMbh[idx]

print(f"  Sample: {N}")

# Damage
def compute_damage(ratio, redshift, n_bins=20):
    edges = np.percentile(redshift, np.linspace(0, 100, n_bins + 1))
    zb = np.clip(np.digitize(redshift, edges) - 1, 0, n_bins - 1)
    damage = np.full(len(ratio), np.nan)
    for b in range(n_bins):
        m = zb == b
        if np.sum(m) < 10: continue
        med = np.median(ratio[m])
        mad = max(1.4826 * np.median(np.abs(ratio[m] - med)), 1e-6)
        damage[m] = np.abs((ratio[m] - med) / mad)
    return damage

damage = compute_damage(log_ratio, Z)
valid = np.isfinite(damage)

# Source model residuals
from numpy.linalg import lstsq
v = valid
X_model = np.column_stack([Z[v], HB_FWHM[v], LEDD[v], LBOL[v], MBH[v], SNR_arr[v], np.ones(np.sum(v))])
coeffs, _, _, _ = lstsq(X_model, damage[v], rcond=None)
residuals = np.full(N, np.nan)
residuals[v] = damage[v] - X_model @ coeffs

# ============================================================
# THE TUBE TEST
# ============================================================
print("\n" + "=" * 70)
print("TUBE TEST — Angular proximity + redshift proximity")
print("=" * 70)

# Build 3D tree: (RA, Dec, z*scale)
# Scale z so that Δz=0.05 ~ 1 degree on sky
z_scale = 20.0  # 0.05 * 20 = 1 degree equivalent

coords_3d = np.column_stack([
    RA * np.cos(np.radians(DEC)),  # project RA by cos(dec)
    DEC,
    Z * z_scale
])

# Only use valid objects
valid_idx = np.where(valid)[0]
coords_valid = coords_3d[valid_idx]
resid_valid = residuals[valid_idx]
damage_valid = damage[valid_idx]
z_valid = Z[valid_idx]

print(f"  Building KDTree with {len(valid_idx)} objects...")
tree = cKDTree(coords_valid)

# For different angular separations, find pairs and test correlation
for max_ang_deg, max_dz in [(0.5, 0.03), (1.0, 0.05), (2.0, 0.05), (5.0, 0.1)]:
    max_r = np.sqrt(max_ang_deg**2 + (max_dz * z_scale)**2)
    
    # Sample: pick random subset to avoid O(N²)
    np.random.seed(42)
    n_sample = min(10000, len(valid_idx))
    sample = np.random.choice(len(valid_idx), n_sample, replace=False)
    
    pair_r1 = []
    pair_r2 = []
    pair_d1 = []
    pair_d2 = []
    
    for i in sample:
        neighbors = tree.query_ball_point(coords_valid[i], max_r)
        for j in neighbors:
            if j <= i: continue  # avoid double-counting
            # Verify angular and z separation
            dra = (RA[valid_idx[i]] - RA[valid_idx[j]]) * np.cos(np.radians(DEC[valid_idx[i]]))
            ddec = DEC[valid_idx[i]] - DEC[valid_idx[j]]
            ang_sep = np.sqrt(dra**2 + ddec**2)
            dz = abs(Z[valid_idx[i]] - Z[valid_idx[j]])
            
            if ang_sep < max_ang_deg and dz < max_dz:
                pair_r1.append(resid_valid[i])
                pair_r2.append(resid_valid[j])
                pair_d1.append(damage_valid[i])
                pair_d2.append(damage_valid[j])
    
    if len(pair_r1) > 20:
        r_resid, p_resid = stats.spearmanr(pair_r1, pair_r2)
        r_damage, p_damage = stats.spearmanr(pair_d1, pair_d2)
        print(f"\n  Tube: <{max_ang_deg}° sky, Δz<{max_dz}")
        print(f"    Pairs found: {len(pair_r1)}")
        print(f"    Residual correlation: ρ = {r_resid:+.4f}, p = {p_resid:.2e}")
        print(f"    Raw damage correlation: ρ = {r_damage:+.4f}, p = {p_damage:.2e}")
        print(f"    → {'PATH MATTERS' if p_resid < 0.01 and r_resid > 0.02 else 'No path signal'}")
    else:
        print(f"\n  Tube: <{max_ang_deg}° sky, Δz<{max_dz}: too few pairs ({len(pair_r1)})")

# ============================================================
# RANDOM PAIR CONTROL — pairs at SAME z but DIFFERENT sky position
# ============================================================
print("\n" + "=" * 70)
print("CONTROL — Same z, different sky position")
print("=" * 70)

np.random.seed(123)
ctrl_r1 = []; ctrl_r2 = []

for _ in range(50000):
    i, j = np.random.choice(len(valid_idx), 2, replace=False)
    dz = abs(z_valid[i] - z_valid[j])
    if dz < 0.05:
        ctrl_r1.append(resid_valid[i])
        ctrl_r2.append(resid_valid[j])

if len(ctrl_r1) > 100:
    r_ctrl, p_ctrl = stats.spearmanr(ctrl_r1, ctrl_r2)
    print(f"  Random pairs (Δz<0.05, any sky sep): {len(ctrl_r1)}")
    print(f"  Residual correlation: ρ = {r_ctrl:+.4f}, p = {p_ctrl:.2e}")

# ============================================================
# PLATE ANALYSIS
# ============================================================
if HAS_PLATE:
    print("\n" + "=" * 70)
    print("PLATE ANALYSIS — Instrumental vs cosmic")
    print("=" * 70)
    
    PLATE = plate[idx]
    unique_plates = np.unique(PLATE[valid])
    
    # Compute per-plate mean residual
    plate_resid = {}
    plate_n = {}
    for p in unique_plates:
        mask = valid & (PLATE == p)
        if np.sum(mask) >= 5:
            plate_resid[p] = np.mean(residuals[mask])
            plate_n[p] = np.sum(mask)
    
    print(f"  Unique plates with ≥5 objects: {len(plate_resid)}")
    
    pr_vals = np.array(list(plate_resid.values()))
    pr_counts = np.array(list(plate_n.values()))
    
    print(f"  Plate residual std: {np.std(pr_vals):.4f}")
    print(f"  Per-object residual std: {np.nanstd(residuals):.4f}")
    print(f"  Plate explains: ~{100 * np.var(pr_vals) / np.nanvar(residuals):.2f}% of variance")
    
    # F-test: do plates explain significant variance?
    # Between-plate vs within-plate variance
    grand_mean = np.nanmean(residuals[valid])
    ss_between = sum(plate_n[p] * (plate_resid[p] - grand_mean)**2 for p in plate_resid)
    ss_within = sum(sum((residuals[valid & (PLATE == p)] - plate_resid[p])**2) 
                     for p in plate_resid if np.sum(valid & (PLATE == p)) >= 5)
    
    k = len(plate_resid)
    n_total = sum(plate_n.values())
    ms_between = ss_between / (k - 1)
    ms_within = ss_within / (n_total - k)
    F_plate = ms_between / ms_within
    
    from scipy.stats import f as f_dist
    p_plate = 1 - f_dist.cdf(F_plate, k-1, n_total-k)
    print(f"  ANOVA: F = {F_plate:.2f}, p = {p_plate:.2e}")
    print(f"  → {'PLATE MATTERS' if p_plate < 0.01 else 'Plates dont explain residuals'}")
    
    # MJD analysis: observation date
    MJD = mjd[idx]
    rho_mjd, p_mjd = stats.spearmanr(MJD[valid], residuals[valid])
    print(f"\n  MJD (observation date) vs residual: ρ = {rho_mjd:+.4f}, p = {p_mjd:.2e}")
    
    # Same plate pairs vs different plate pairs
    print("\n  --- Same-plate vs different-plate pair correlation ---")
    np.random.seed(456)
    same_r1 = []; same_r2 = []
    diff_r1 = []; diff_r2 = []
    
    for _ in range(100000):
        i, j = np.random.choice(len(valid_idx), 2, replace=False)
        dz = abs(z_valid[i] - z_valid[j])
        if dz > 0.05: continue
        
        pi = PLATE[valid_idx[i]]
        pj = PLATE[valid_idx[j]]
        
        if pi == pj:
            same_r1.append(resid_valid[i])
            same_r2.append(resid_valid[j])
        else:
            diff_r1.append(resid_valid[i])
            diff_r2.append(resid_valid[j])
    
    if len(same_r1) > 50:
        r_same, p_same = stats.spearmanr(same_r1, same_r2)
        print(f"  Same plate (Δz<0.05): N={len(same_r1)}, ρ = {r_same:+.4f}, p = {p_same:.2e}")
    if len(diff_r1) > 50:
        r_diff, p_diff = stats.spearmanr(diff_r1, diff_r2)
        print(f"  Diff plate (Δz<0.05): N={len(diff_r1)}, ρ = {r_diff:+.4f}, p = {p_diff:.2e}")

# ============================================================
# THE LEVER: What unmeasured property drives the 84%?
# ============================================================
print("\n" + "=" * 70)
print("THE 84% — What drives it?")
print("=" * 70)

# Check if residuals correlate with any property we HAVEN'T used
ebv = fix(data['EBV'])[idx]  # galactic extinction
feii = fix(data['FEII_OPT_EW'])[idx]  # FeII strength (EV1)
bal = fix(data['BAL_FLAG_VI'])[idx] if 'BAL_FLAG_VI' in [c.name for c in hdu[1].columns] else None

print("  Correlations with unused properties:")

# E(B-V) galactic extinction
v_ebv = valid & np.isfinite(ebv)
rho_e, p_e = stats.spearmanr(ebv[v_ebv], residuals[v_ebv])
print(f"  E(B-V) galactic:  ρ = {rho_e:+.4f}, p = {p_e:.2e}")

# FeII (Eigenvector 1 — AGN diversity axis)
v_fe = valid & np.isfinite(feii) & (feii > 0)
if np.sum(v_fe) > 100:
    rho_f, p_f = stats.spearmanr(feii[v_fe], residuals[v_fe])
    print(f"  FeII EW (EV1):    ρ = {rho_f:+.4f}, p = {p_f:.2e}")

# BAL flag
if bal is not None:
    v_bal = valid & np.isfinite(bal)
    bal_yes = v_bal & (bal > 0)
    bal_no = v_bal & (bal == 0)
    if np.sum(bal_yes) > 50 and np.sum(bal_no) > 50:
        t_bal, p_bal = stats.ttest_ind(residuals[bal_yes], residuals[bal_no])
        print(f"  BAL quasars:      Δresid = {np.mean(residuals[bal_yes]) - np.mean(residuals[bal_no]):+.4f}, p = {p_bal:.2e}")
        print(f"    BAL: {np.mean(residuals[bal_yes]):+.4f} (N={np.sum(bal_yes)})")
        print(f"    non-BAL: {np.mean(residuals[bal_no]):+.4f} (N={np.sum(bal_no)})")

# Radio loudness
try:
    first = fix(data['FIRST_MATCHED'])[idx]
    v_radio = valid & np.isfinite(first)
    radio_loud = v_radio & (first > 0)
    radio_quiet = v_radio & (first == 0)
    if np.sum(radio_loud) > 50 and np.sum(radio_quiet) > 50:
        t_r, p_r = stats.ttest_ind(residuals[radio_loud], residuals[radio_quiet])
        print(f"  Radio loud:       Δresid = {np.mean(residuals[radio_loud]) - np.mean(residuals[radio_quiet]):+.4f}, p = {p_r:.2e}")
except:
    pass

# Continuum slope
try:
    alpha = fix(data['ALPHA_NU'])[idx]
    v_a = valid & np.isfinite(alpha)
    rho_a, p_a = stats.spearmanr(alpha[v_a], residuals[v_a])
    print(f"  Spectral index:   ρ = {rho_a:+.4f}, p = {p_a:.2e}")
except:
    pass

# Absolute magnitude
try:
    mi = fix(data['MI'])[idx]
    v_mi = valid & np.isfinite(mi)
    rho_mi, p_mi = stats.spearmanr(mi[v_mi], residuals[v_mi])
    print(f"  Abs mag (Mi):     ρ = {rho_mi:+.4f}, p = {p_mi:.2e}")
except:
    pass

print("\n✅ Done.")
