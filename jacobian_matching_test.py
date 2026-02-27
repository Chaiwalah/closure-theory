#!/usr/bin/env python3
"""
JACOBIAN MATCHING TEST — GPT's Last Door
==========================================
GPT's claim: the matched-state test doesn't control for observed-frame
diagnostic geometry. Two quasars with identical physical states can be
differentially ill-conditioned if their lines land on different observed
wavelengths (sky-line forests, telluric, continuum shape).

This test: match pairs on FWHM, Ledd, Lbol, SNR **AND** observed-frame
diagnostic geometry (line wavelength neighborhoods, continuum slope),
then test if high-z is still more damaged.

If yes → inference-layer fragility is dead.
"""

import numpy as np
from astropy.io import fits
from scipy import stats
from scipy.spatial import cKDTree
import os, warnings
warnings.filterwarnings('ignore')

OUTDIR = 'results_jacobian'
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 70)
print("JACOBIAN MATCHING TEST — Closing GPT's Last Door")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data

def fix(arr): return arr.astype(arr.dtype.newbyteorder('='))

z = fix(data['Z_DR16Q'])
snr = fix(data['SN_MEDIAN_ALL'])
oiii_flux = fix(data['OIII5007'][:, 4])
hbeta_flux = fix(data['HBETA'][:, 4])
hbeta_fwhm = fix(data['HBETA'][:, 3])
logLbol = fix(data['LOGLBOL'])
logMbh = fix(data['LOGMBH'])
logLedd = fix(data['LOGLEDD_RATIO'])

# Rest-frame line wavelengths (Angstroms)
OIII_REST = 5007.0
HB_REST = 4861.0
SII_REST = 6718.0  # [SII] 6718
OII_REST = 3728.0  # [OII] 3728
NII_REST = 6585.0  # [NII] 6585
HA_REST = 6563.0

# Quality cuts
good = ((z > 0.3) & (z < 2.0) & (snr > 3) &
        (oiii_flux > 0) & (hbeta_flux > 0) &
        np.isfinite(oiii_flux) & np.isfinite(hbeta_flux) &
        np.isfinite(logLbol) & np.isfinite(logMbh) & np.isfinite(logLedd) &
        (logLbol > 0) & (logMbh > 0) &
        np.isfinite(hbeta_fwhm) & (hbeta_fwhm > 0))

idx = np.where(good)[0]
N = len(idx)
Z = z[idx]
print(f"  Base sample: {N}")

log_ratio = np.log10(oiii_flux[idx] / hbeta_flux[idx])
HB_FWHM = np.log10(np.maximum(hbeta_fwhm[idx], 1))
LEDD = logLedd[idx]
LBOL = logLbol[idx]
SNR = snr[idx]

# ============================================================
# COMPUTE OBSERVED-FRAME DIAGNOSTIC GEOMETRY
# ============================================================
print("\n--- Computing observed-frame diagnostic geometry ---")

# 1. Observed wavelengths of key lines
oiii_obs = OIII_REST * (1 + Z)
hb_obs = HB_REST * (1 + Z)

# 2. Observed-frame line separation (proxy for continuum confusion)
line_sep_obs = oiii_obs - hb_obs  # grows with z

# 3. Sky-line density score
# Major OH sky emission bands (Angstroms, observed frame)
# These are the worst offenders in optical spectroscopy
SKY_LINES = np.array([
    5577.3,  # [OI] - strongest sky line
    5889.9, 5895.9,  # Na I D
    6300.3, 6363.8,  # [OI]
    6498.7, 6533.0, 6553.6,  # OH
    6863.9, 6912.6, 6923.2, 6939.5,  # B-band telluric
    7244.0, 7276.0, 7316.0, 7340.0,  # OH forest
    7584.0, 7632.0, 7711.0, 7750.0, 7794.0,  # A-band telluric
    7821.0, 7853.0, 7913.0, 7950.0, 7993.0,  # OH forest
    8344.0, 8399.0, 8430.0, 8465.0,  # OH
    8791.0, 8827.0, 8885.0, 8919.0, 8958.0,  # OH
])

def sky_density(obs_wave, window=30.0):
    """Count sky lines within ±window Angstroms of observed wavelength"""
    return np.array([np.sum(np.abs(SKY_LINES - w) < window) for w in obs_wave])

oiii_sky = sky_density(oiii_obs)
hb_sky = sky_density(hb_obs)
total_sky = oiii_sky + hb_sky

print(f"  Sky density range: {total_sky.min()} - {total_sky.max()}")
print(f"  Median sky density: {np.median(total_sky):.1f}")

# 4. Continuum slope proxy — ratio of observed wavelengths
# (higher z → lines move redder → different continuum regime)
cont_slope_proxy = np.log10(oiii_obs / hb_obs)  # nearly constant ~0.013, but shift matters

# 5. Jacobian proxy: local diagnostic sensitivity
# The Jacobian of the OIII/Hβ ratio w.r.t. flux perturbations
# ∂(log R) / ∂(log F_OIII) = 1, ∂(log R) / ∂(log F_Hβ) = -1
# But the NOISE on each flux depends on observed-frame position
# Proxy: SNR ratio between lines (available from data)
oiii_sigma = fix(data['OIII5007'][:, 3])[idx]  # line width as noise proxy
hb_sigma = fix(data['HBETA'][:, 3])[idx]

# Actually, better proxy: the observed-frame wavelength itself
# determines detector QE, sky background, etc.
# Use oiii_obs as the "diagnostic geometry" coordinate
# since it captures where in detector space the measurement happens

# ============================================================
# Z-BIN DAMAGE
# ============================================================
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

# ============================================================
# TEST A: ORIGINAL MATCHED PAIRS (baseline — reproduce corner_the_rat)
# ============================================================
print("\n" + "=" * 70)
print("TEST A: BASELINE MATCHED PAIRS (source only)")
print("=" * 70)

# Split into z-bins
z_med = np.median(Z)
lo_mask = valid & (Z < z_med - 0.1)
hi_mask = valid & (Z > z_med + 0.1)

# Feature vectors for matching: FWHM, Ledd, Lbol, SNR
lo_idx = np.where(lo_mask)[0]
hi_idx = np.where(hi_mask)[0]

def standardize(arr):
    return (arr - np.mean(arr)) / max(np.std(arr), 1e-10)

# Match on source properties only
feat_lo_A = np.column_stack([
    standardize(HB_FWHM[lo_idx]),
    standardize(LEDD[lo_idx]),
    standardize(LBOL[lo_idx]),
    standardize(SNR[lo_idx])
])

feat_hi_A = np.column_stack([
    standardize(HB_FWHM[hi_idx]),
    standardize(LEDD[hi_idx]),
    standardize(LBOL[hi_idx]),
    standardize(SNR[hi_idx])
])

print(f"  Low-z pool: {len(lo_idx)}, High-z pool: {len(hi_idx)}")

tree_lo = cKDTree(feat_lo_A)
dists, matches = tree_lo.query(feat_hi_A, k=1)

# Quality filter: max distance
d_thresh = np.percentile(dists, 75)
good_match = dists < d_thresh
n_pairs_A = np.sum(good_match)
print(f"  Matched pairs (source only): {n_pairs_A}")

dmg_hi_A = damage[hi_idx[good_match]]
dmg_lo_A = damage[lo_idx[matches[good_match]]]
delta_A = np.nanmean(dmg_hi_A) - np.nanmean(dmg_lo_A)
t_A, p_A = stats.ttest_rel(dmg_hi_A[np.isfinite(dmg_hi_A) & np.isfinite(dmg_lo_A)],
                             dmg_lo_A[np.isfinite(dmg_hi_A) & np.isfinite(dmg_lo_A)])

print(f"  Δdamage (hi - lo): {delta_A:+.4f}")
print(f"  t = {t_A:.2f}, p = {p_A:.2e}")

# ============================================================
# TEST B: JACOBIAN-MATCHED PAIRS (source + observed-frame geometry)
# ============================================================
print("\n" + "=" * 70)
print("TEST B: JACOBIAN-MATCHED PAIRS (source + obs-frame geometry)")
print("=" * 70)
print("  Adding: observed OIII wavelength, sky density, line separation")

# Now match on source + geometry
feat_lo_B = np.column_stack([
    standardize(HB_FWHM[lo_idx]),
    standardize(LEDD[lo_idx]),
    standardize(LBOL[lo_idx]),
    standardize(SNR[lo_idx]),
    standardize(oiii_obs[lo_idx]),
    standardize(total_sky[lo_idx].astype(float)),
    standardize(line_sep_obs[lo_idx])
])

feat_hi_B = np.column_stack([
    standardize(HB_FWHM[hi_idx]),
    standardize(LEDD[hi_idx]),
    standardize(LBOL[hi_idx]),
    standardize(SNR[hi_idx]),
    standardize(oiii_obs[hi_idx]),
    standardize(total_sky[hi_idx].astype(float)),
    standardize(line_sep_obs[hi_idx])
])

tree_lo_B = cKDTree(feat_lo_B)
dists_B, matches_B = tree_lo_B.query(feat_hi_B, k=1)

d_thresh_B = np.percentile(dists_B, 75)
good_match_B = dists_B < d_thresh_B
n_pairs_B = np.sum(good_match_B)
print(f"  Matched pairs (source + geometry): {n_pairs_B}")

dmg_hi_B = damage[hi_idx[good_match_B]]
dmg_lo_B = damage[lo_idx[matches_B[good_match_B]]]
delta_B = np.nanmean(dmg_hi_B) - np.nanmean(dmg_lo_B)
fin_B = np.isfinite(dmg_hi_B) & np.isfinite(dmg_lo_B)
t_B, p_B = stats.ttest_rel(dmg_hi_B[fin_B], dmg_lo_B[fin_B])

print(f"  Δdamage (hi - lo): {delta_B:+.4f}")
print(f"  t = {t_B:.2f}, p = {p_B:.2e}")

# ============================================================
# TEST C: EXTREME GEOMETRY MATCHING (belt AND suspenders)
# ============================================================
print("\n" + "=" * 70)
print("TEST C: EXTREME GEOMETRY MATCHING")
print("=" * 70)
print("  Matching on obs-frame wavelength BINS (force same detector region)")

# Bin by observed OIII wavelength (100 Å bins)
oiii_obs_bin = np.floor(oiii_obs / 100).astype(int)
# Bin by sky density
sky_bin = total_sky

# Within each (oiii_obs_bin, sky_density) cell, match on source properties
unique_obs_bins = np.unique(oiii_obs_bin)
pair_hi_list = []
pair_lo_list = []

for ob in unique_obs_bins:
    for sd in range(int(total_sky.max()) + 1):
        cell_lo = lo_mask & (oiii_obs_bin == ob) & (total_sky == sd)
        cell_hi = hi_mask & (oiii_obs_bin == ob) & (total_sky == sd)
        n_lo_cell = np.sum(cell_lo)
        n_hi_cell = np.sum(cell_hi)
        if n_lo_cell < 5 or n_hi_cell < 5:
            continue
        
        cl_idx = np.where(cell_lo)[0]
        ch_idx = np.where(cell_hi)[0]
        
        # Match on source within cell
        feat_cl = np.column_stack([
            standardize(HB_FWHM[cl_idx]),
            standardize(LEDD[cl_idx]),
            standardize(LBOL[cl_idx]),
            standardize(SNR[cl_idx])
        ])
        feat_ch = np.column_stack([
            standardize(HB_FWHM[ch_idx]),
            standardize(LEDD[ch_idx]),
            standardize(LBOL[ch_idx]),
            standardize(SNR[ch_idx])
        ])
        
        tree_cl = cKDTree(feat_cl)
        d_cl, m_cl = tree_cl.query(feat_ch, k=1)
        thresh_cl = np.percentile(d_cl, 75) if len(d_cl) > 4 else np.inf
        gm = d_cl < thresh_cl
        
        if np.sum(gm) > 0:
            pair_hi_list.append(damage[ch_idx[gm]])
            pair_lo_list.append(damage[cl_idx[m_cl[gm]]])

if pair_hi_list:
    all_hi_C = np.concatenate(pair_hi_list)
    all_lo_C = np.concatenate(pair_lo_list)
    fin_C = np.isfinite(all_hi_C) & np.isfinite(all_lo_C)
    n_pairs_C = np.sum(fin_C)
    delta_C = np.nanmean(all_hi_C[fin_C]) - np.nanmean(all_lo_C[fin_C])
    t_C, p_C = stats.ttest_rel(all_hi_C[fin_C], all_lo_C[fin_C])
    print(f"  Matched pairs (cell-matched): {n_pairs_C}")
    print(f"  Δdamage (hi - lo): {delta_C:+.4f}")
    print(f"  t = {t_C:.2f}, p = {p_C:.2e}")
else:
    print("  No cell pairs found")

# ============================================================
# TEST D: SKY-LINE CONTAMINATION CONTROL
# ============================================================
print("\n" + "=" * 70)
print("TEST D: SKY-LINE CONTAMINATION CONTROL")
print("=" * 70)
print("  Split by sky density: do clean sightlines show LESS damage?")

lo_sky = total_sky == 0
hi_sky = total_sky >= 3

for label, mask in [("Clean (0 sky lines)", lo_sky & valid),
                     ("Contaminated (3+ sky lines)", hi_sky & valid)]:
    if np.sum(mask) < 100:
        print(f"  {label}: too few ({np.sum(mask)})")
        continue
    sl, _, _, p, _ = stats.linregress(Z[mask], damage[mask])
    print(f"  {label}: N={np.sum(mask)}, slope={sl:+.4f}, p={p:.2e}")

# ============================================================
# TEST E: SUSCEPTIBILITY GATING WITHIN GEOMETRY-MATCHED PAIRS
# ============================================================
print("\n" + "=" * 70)
print("TEST E: SUSCEPTIBILITY × GEOMETRY INTERACTION")
print("=" * 70)
print("  Does susceptibility gating survive within fixed geometry?")

# Use Test B pairs, stratify by susceptibility
S = HB_FWHM  # susceptibility proxy (non-circular)
s_terciles = np.percentile(S, [33, 67])

for label, s_lo, s_hi in [("Low S", -np.inf, s_terciles[0]),
                            ("Mid S", s_terciles[0], s_terciles[1]),
                            ("High S", s_terciles[1], np.inf)]:
    # Get pairs from Test B where BOTH members are in this S range
    hi_s = S[hi_idx[good_match_B]]
    lo_s = S[lo_idx[matches_B[good_match_B]]]
    in_range = (hi_s >= s_lo) & (hi_s < s_hi) & (lo_s >= s_lo) & (lo_s < s_hi)
    fin = in_range & np.isfinite(dmg_hi_B) & np.isfinite(dmg_lo_B)
    
    if np.sum(fin) < 30:
        print(f"  {label}: too few ({np.sum(fin)})")
        continue
    
    d_hi = dmg_hi_B[fin]
    d_lo = dmg_lo_B[fin]
    delta = np.mean(d_hi) - np.mean(d_lo)
    t, p = stats.ttest_rel(d_hi, d_lo)
    print(f"  {label}: N={np.sum(fin)}, Δ={delta:+.4f}, t={t:.2f}, p={p:.2e}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"  Test A (source-only matching):     Δ = {delta_A:+.4f}, p = {p_A:.2e}")
print(f"  Test B (source + geometry):         Δ = {delta_B:+.4f}, p = {p_B:.2e}")
if pair_hi_list:
    print(f"  Test C (cell-matched geometry):     Δ = {delta_C:+.4f}, p = {p_C:.2e}")
print(f"\n  If Test B and C show significant Δ → inference-layer fragility is DEAD")
print(f"  If Test E shows S-gating within geometry → it's the operator, not the telescope")

# Save
with open(f'{OUTDIR}/summary.txt', 'w') as f:
    f.write("JACOBIAN MATCHING TEST RESULTS\n")
    f.write(f"Test A (source only): Δ={delta_A:+.4f}, p={p_A:.2e}\n")
    f.write(f"Test B (source+geom): Δ={delta_B:+.4f}, p={p_B:.2e}\n")
    if pair_hi_list:
        f.write(f"Test C (cell-match):  Δ={delta_C:+.4f}, p={p_C:.2e}\n")

print("\n✅ Done. Results in", OUTDIR)
