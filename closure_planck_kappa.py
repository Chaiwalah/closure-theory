#!/usr/bin/env python3
"""
PLANCK KAPPA CROSS-MATCH — The test everyone agreed on.

Does foreground lensing convergence (κ) predict early breaking?
- Positive κ = foreground mass overdensity = structure
- Negative κ = foreground mass underdensity = void
- If structure protects: breakers should have LOWER κ
- If voids damage: breakers should cluster in negative κ regions

Data: Planck 2018 CMB lensing convergence (R3.00)
      MV (minimum variance) reconstruction
"""

import numpy as np
import healpy as hp
from astropy.io import fits
from scipy import stats
import json, os, warnings, gzip
warnings.filterwarnings('ignore')

os.makedirs('results_planck_kappa', exist_ok=True)

print("=" * 70)
print("PLANCK κ CROSS-MATCH — Does foreground mass predict breaking?")
print("=" * 70)

# ============================================================
# Load Planck lensing convergence map
# ============================================================
print("\nLoading Planck 2018 lensing convergence map...")

# The data is in spherical harmonic coefficients (alm format)
# MV = minimum variance (best combined estimate)
# dat_klm = data kappa alm; mf_klm = mean field alm
# kappa_alm = dat_klm - mf_klm

dat_file = 'data/planck_lensing/COM_Lensing_4096_R3.00/MV/dat_klm.fits'
mf_file = 'data/planck_lensing/COM_Lensing_4096_R3.00/MV/mf_klm.fits'
mask_file = 'data/planck_lensing/COM_Lensing_4096_R3.00/mask.fits.gz'

# Read alm coefficients
dat_alm = hp.read_alm(dat_file)
mf_alm = hp.read_alm(mf_file)

# Subtract mean field
kappa_alm = dat_alm - mf_alm

lmax_alm = hp.Alm.getlmax(len(kappa_alm))
print(f"alm lmax: {lmax_alm}")

# Apply Wiener filter using noise power spectrum
# nlkk columns: L, C_L^kk (signal), N_L^kk (noise)  
nlkk = np.loadtxt('data/planck_lensing/COM_Lensing_4096_R3.00/MV/nlkk.dat')
L_arr = nlkk[:, 0].astype(int)
CL_signal = nlkk[:, 1]
NL_noise = nlkk[:, 2]

# Wiener filter: W_L = C_L / (C_L + N_L)
# Also cut at Lmax=400 where S/N is still reasonable
LMAX_CUT = 400
print(f"Applying Wiener filter with Lmax={LMAX_CUT}...")

for l in range(lmax_alm + 1):
    for m in range(0, min(l + 1, lmax_alm + 1)):
        idx = hp.Alm.getidx(lmax_alm, l, m)
        if idx >= len(kappa_alm):
            continue
        if l > LMAX_CUT or l >= len(CL_signal):
            kappa_alm[idx] = 0
        else:
            cl = CL_signal[l]
            nl = NL_noise[l] 
            wiener = cl / (cl + nl) if (cl + nl) > 0 else 0
            kappa_alm[idx] *= wiener

# Convert to pixel map
NSIDE = 512  # sufficient for degree-scale features
print(f"Converting alm to pixel map (NSIDE={NSIDE})...")
kappa_map = hp.alm2map(kappa_alm, NSIDE, verbose=False)

# Load and apply mask
print("Loading mask...")
mask_raw = hp.read_map(mask_file, verbose=False)
# Mask might be at different NSIDE, upgrade/downgrade
mask_nside = hp.get_nside(mask_raw)
if mask_nside != NSIDE:
    print(f"  Resampling mask from NSIDE={mask_nside} to {NSIDE}...")
    mask = hp.ud_grade(mask_raw, NSIDE)
else:
    mask = mask_raw

# Apply mask threshold (mask > 0.5 = valid)
valid_mask = mask > 0.5
kappa_map[~valid_mask] = hp.UNSEEN

print(f"Map stats (unmasked): mean={np.mean(kappa_map[valid_mask]):.6f}, "
      f"std={np.std(kappa_map[valid_mask]):.6f}, "
      f"valid pixels={valid_mask.sum()} / {len(valid_mask)}")

# ============================================================
# Load quasar data
# ============================================================
print("\nLoading DR16Q quasars...")
hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
ra_q = d['RA']
dec_q = d['DEC']

# ============================================================
# Recreate groups
# ============================================================
print("Recomputing rank agreements...")

trans_mask = (z_all >= 0.75) & (z_all < 1.15)
valid = trans_mask & np.isfinite(mg_ew) & np.isfinite(hb_ew) & (mg_ew != 0) & (hb_ew != 0)
z_v = z_all[valid]
mg_v = mg_ew[valid]
hb_v = hb_ew[valid]

window = 0.05
rank_agreement = np.full(valid.sum(), np.nan)

for i in range(valid.sum()):
    z_i = z_v[i]
    local = (z_v >= z_i - window) & (z_v < z_i + window)
    if local.sum() < 30:
        continue
    mg_pct = stats.percentileofscore(mg_v[local], mg_v[i]) / 100
    hb_pct = stats.percentileofscore(hb_v[local], hb_v[i]) / 100
    rank_agreement[i] = 1 - abs(mg_pct - hb_pct)

valid_ra = np.isfinite(rank_agreement)
valid_indices = np.where(valid)[0]

early_zone = (z_v >= 0.80) & (z_v < 0.95)
late_zone = (z_v >= 1.00) & (z_v < 1.15)

ra_early = rank_agreement[early_zone & valid_ra]
ra_late = rank_agreement[late_zone & valid_ra]

early_thresh = np.percentile(ra_early, 20)
late_thresh = np.percentile(ra_late, 80)

early_breakers = early_zone & valid_ra & (rank_agreement <= early_thresh)
early_normal = early_zone & valid_ra & (rank_agreement > np.percentile(ra_early, 40)) & (rank_agreement < np.percentile(ra_early, 60))
late_holders = late_zone & valid_ra & (rank_agreement >= late_thresh)
late_normal = late_zone & valid_ra & (rank_agreement > np.percentile(ra_late, 40)) & (rank_agreement < np.percentile(ra_late, 60))

eb_idx = valid_indices[early_breakers]
en_idx = valid_indices[early_normal]
lh_idx = valid_indices[late_holders]
ln_idx = valid_indices[late_normal]

print(f"Groups: EB={len(eb_idx)}, EN={len(en_idx)}, LH={len(lh_idx)}, LN={len(ln_idx)}")

# ============================================================
# Extract κ values at quasar positions
# ============================================================
print(f"\n{'='*70}")
print("Extracting Planck κ at quasar positions...")
print(f"{'='*70}")

def get_kappa(indices):
    """Get kappa values at quasar sky positions."""
    # Convert RA/DEC to healpy theta/phi
    # healpy: theta = colatitude (0=north pole), phi = longitude (RA)
    theta = np.radians(90.0 - dec_q[indices])  # DEC -> colatitude
    phi = np.radians(ra_q[indices])  # RA -> longitude
    
    # Get pixel indices
    pix = hp.ang2pix(NSIDE, theta, phi)
    
    # Get kappa values
    kvals = kappa_map[pix]
    
    # Check mask
    masked = ~valid_mask[pix]
    kvals[masked] = np.nan
    
    return kvals

k_eb = get_kappa(eb_idx)
k_en = get_kappa(en_idx)
k_lh = get_kappa(lh_idx)
k_ln = get_kappa(ln_idx)

# Remove NaN (masked regions)
k_eb_v = k_eb[np.isfinite(k_eb)]
k_en_v = k_en[np.isfinite(k_en)]
k_lh_v = k_lh[np.isfinite(k_lh)]
k_ln_v = k_ln[np.isfinite(k_ln)]

print(f"\n  Valid kappa values:")
print(f"  Early Breakers: {len(k_eb_v)}/{len(k_eb)} ({100*len(k_eb_v)/len(k_eb):.1f}%)")
print(f"  Early Normal:   {len(k_en_v)}/{len(k_en)} ({100*len(k_en_v)/len(k_en):.1f}%)")
print(f"  Late Holders:   {len(k_lh_v)}/{len(k_lh)} ({100*len(k_lh_v)/len(k_lh):.1f}%)")
print(f"  Late Normal:    {len(k_ln_v)}/{len(k_ln)} ({100*len(k_ln_v)/len(k_ln):.1f}%)")

# ============================================================
# TEST 1: Do breakers have lower κ than normals?
# ============================================================
print(f"\n{'='*70}")
print("TEST 1: Breaker vs Normal κ distributions")
print(f"{'='*70}")

print(f"\n  EARLY ZONE (z=0.8-0.95):")
print(f"  Breakers κ: median = {np.median(k_eb_v):.6f}, mean = {np.mean(k_eb_v):.6f}")
print(f"  Normal  κ: median = {np.median(k_en_v):.6f}, mean = {np.mean(k_en_v):.6f}")
_, p_eb = stats.mannwhitneyu(k_eb_v, k_en_v, alternative='two-sided')
_, p_eb_less = stats.mannwhitneyu(k_eb_v, k_en_v, alternative='less')
print(f"  MWU p (two-sided): {p_eb:.6f}")
print(f"  MWU p (breakers < normal): {p_eb_less:.6f}")

delta_eb = np.median(k_eb_v) - np.median(k_en_v)
print(f"  Δκ = {delta_eb:.6f}")

if np.median(k_eb_v) < np.median(k_en_v):
    print(f"  → Breakers behind LESS mass → VOID HYPOTHESIS SUPPORTED")
elif np.median(k_eb_v) > np.median(k_en_v):
    print(f"  → Breakers behind MORE mass → structure DAMAGES")
else:
    print(f"  → Same κ")

print(f"\n  LATE ZONE (z=1.0-1.15):")
print(f"  Holders κ: median = {np.median(k_lh_v):.6f}, mean = {np.mean(k_lh_v):.6f}")
print(f"  Normal  κ: median = {np.median(k_ln_v):.6f}, mean = {np.mean(k_ln_v):.6f}")
_, p_lh = stats.mannwhitneyu(k_lh_v, k_ln_v, alternative='two-sided')
_, p_lh_greater = stats.mannwhitneyu(k_lh_v, k_ln_v, alternative='greater')
print(f"  MWU p (two-sided): {p_lh:.6f}")
print(f"  MWU p (holders > normal): {p_lh_greater:.6f}")

delta_lh = np.median(k_lh_v) - np.median(k_ln_v)
print(f"  Δκ = {delta_lh:.6f}")

# ============================================================
# TEST 2: Fraction in negative vs positive κ
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: Fraction in negative κ (void) vs positive κ (structure)")
print(f"{'='*70}")

for label, kvals in [("Early Breakers", k_eb_v), ("Early Normal", k_en_v),
                      ("Late Holders", k_lh_v), ("Late Normal", k_ln_v)]:
    neg_frac = (kvals < 0).mean()
    pos_frac = (kvals > 0).mean()
    print(f"  {label:20s}: neg_κ = {100*neg_frac:.1f}%, pos_κ = {100*pos_frac:.1f}%, mean = {np.mean(kvals):.6f}")

# ============================================================
# TEST 3: Quintile analysis
# ============================================================
print(f"\n{'='*70}")
print("QUINTILE ANALYSIS — κ vs breaker fraction")
print(f"{'='*70}")

# For all early-zone objects
all_early = early_zone & valid_ra
ae_idx = valid_indices[all_early]
k_all_early = get_kappa(ae_idx)
ra_all_early = rank_agreement[all_early]

# Remove NaN
k_valid = np.isfinite(k_all_early)
k_ae = k_all_early[k_valid]
ra_ae = ra_all_early[k_valid]

quintile_edges = np.percentile(k_ae, [0, 20, 40, 60, 80, 100])
print(f"\n  κ quintile edges: {[f'{e:.5f}' for e in quintile_edges]}")

print(f"\n  {'κ range':<25} {'Mean rank agr':>15} {'Breaker frac':>15} {'N':>8}")
print(f"  {'-'*25} {'-'*15} {'-'*15} {'-'*8}")

mean_ras = []
breaker_fracs = []

for i in range(5):
    lo, hi = quintile_edges[i], quintile_edges[i+1]
    if i < 4:
        m = (k_ae >= lo) & (k_ae < hi)
    else:
        m = (k_ae >= lo) & (k_ae <= hi)
    
    ra_bin = ra_ae[m]
    mean_ra = np.mean(ra_bin)
    breaker_frac = (ra_bin <= early_thresh).mean()
    
    mean_ras.append(mean_ra)
    breaker_fracs.append(breaker_frac)
    
    label = f"κ=[{lo:+.5f},{hi:+.5f}]"
    print(f"  {label:<25} {mean_ra:>15.3f} {100*breaker_frac:>14.1f}% {m.sum():>8}")

r_mono, p_mono = stats.spearmanr(range(5), breaker_fracs)
print(f"\n  Monotonicity (κ quintile vs breaker fraction): ρ = {r_mono:+.3f} (p = {p_mono:.3f})")

if r_mono < -0.5 and p_mono < 0.1:
    print(f"  ★ HIGHER κ → FEWER breakers → MASS PROTECTS")
elif r_mono > 0.5 and p_mono < 0.1:
    print(f"  ★ HIGHER κ → MORE breakers → mass DAMAGES")
else:
    print(f"  ~ No clear monotonic relationship")

# ============================================================
# TEST 4: KS tests for distributional differences
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: KS test — distribution shape differences")
print(f"{'='*70}")

ks_eb, p_ks_eb = stats.ks_2samp(k_eb_v, k_en_v)
ks_lh, p_ks_lh = stats.ks_2samp(k_lh_v, k_ln_v)
print(f"  Breakers vs Normal: D={ks_eb:.4f}, p={p_ks_eb:.6f}")
print(f"  Holders vs Normal:  D={ks_lh:.4f}, p={p_ks_lh:.6f}")

# ============================================================
# TEST 5: Extreme tails — do breakers preferentially live in deep voids?
# ============================================================
print(f"\n{'='*70}")
print("TEST 5: Extreme κ tails")
print(f"{'='*70}")

# Bottom 10% of κ (deepest voids)
k10 = np.percentile(k_ae, 10)
# Top 10% of κ (densest structures)
k90 = np.percentile(k_ae, 90)

for label, kvals, ref_frac in [("Early Breakers", k_eb_v, None), ("Early Normal", k_en_v, None)]:
    deep_void = (kvals < k10).mean()
    dense_struct = (kvals > k90).mean()
    print(f"  {label:20s}: in deep voids (κ<{k10:.5f}) = {100*deep_void:.1f}%, in dense struct (κ>{k90:.5f}) = {100*dense_struct:.1f}%")

# ============================================================
# Save results
# ============================================================
results = {
    'kappa_breakers_median': float(np.median(k_eb_v)),
    'kappa_normals_median': float(np.median(k_en_v)),
    'kappa_holders_median': float(np.median(k_lh_v)),
    'kappa_holders_normal_median': float(np.median(k_ln_v)),
    'delta_kappa_breakers': float(delta_eb),
    'delta_kappa_holders': float(delta_lh),
    'p_breakers_vs_normal': float(p_eb),
    'p_breakers_less': float(p_eb_less),
    'p_holders_vs_normal': float(p_lh),
    'p_holders_greater': float(p_lh_greater),
    'ks_breakers': float(ks_eb),
    'p_ks_breakers': float(p_ks_eb),
    'ks_holders': float(ks_lh),
    'p_ks_holders': float(p_ks_lh),
    'monotonicity_rho': float(r_mono),
    'monotonicity_p': float(p_mono),
    'quintile_breaker_fracs': [float(x) for x in breaker_fracs],
    'n_valid': {
        'breakers': int(len(k_eb_v)),
        'normals': int(len(k_en_v)),
        'holders': int(len(k_lh_v)),
        'holders_normal': int(len(k_ln_v)),
    }
}

with open('results_planck_kappa/planck_kappa_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n{'='*70}")
print("Results saved to results_planck_kappa/")
print(f"{'='*70}")
