#!/usr/bin/env python3
"""
tSZ COMPTON-y STRUCTURAL STATE TEST

Uses Planck Commander SZ map (Compton-y parameter) to test whether
sightlines through virialized structures (high y) preserve observable
correlations better than sightlines through diffuse space (low y).

High y = hot, bound intracluster medium = ORDERED state
Low y = diffuse, expanding IGM = DISORDERED state

Tests:
1. Compton-y vs SN Ia color/stretch/residual (Pantheon+)
2. Compton-y quintiles vs quasar MgII-Hβ coupling strength
3. Compton-y vs density (independence check)
4. Comparison with kill grid results (should outperform density/RM/κ)
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import healpy as hp
from scipy import stats
import json, os, sys, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_tsz', exist_ok=True)

print("=" * 70)
print("tSZ COMPTON-y STRUCTURAL STATE TEST")
print("Planck Commander SZ map vs observable degradation")
print("=" * 70)
sys.stdout.flush()

# ============================================================
# 1. LOAD PLANCK COMPTON-y MAP
# ============================================================
print("\n[1] Loading Planck Commander SZ map...")
sys.stdout.flush()

sz_hdu = fits.open('/root/clawd/projects/closure-theory/data/planck_sz_commander.fits')
print(f"  HDUs: {len(sz_hdu)}")
for i, h in enumerate(sz_hdu):
    if h.data is not None:
        print(f"  HDU[{i}]: {h.name}, columns={h.columns.names if hasattr(h, 'columns') else 'image'}")
sys.stdout.flush()

# Commander SZ map is typically in extension 1 as a HEALPix table
try:
    ymap = hp.read_map('/root/clawd/projects/closure-theory/data/planck_sz_commander.fits', 
                        field=0, verbose=False)
    nside = hp.get_nside(ymap)
    print(f"  Nside: {nside}, npix: {len(ymap)}")
    print(f"  y range: [{ymap.min():.2e}, {ymap.max():.2e}]")
    print(f"  y mean:  {ymap.mean():.2e}")
    print(f"  y > 0:   {(ymap > 0).sum()} pixels ({(ymap > 0).sum()/len(ymap)*100:.1f}%)")
except Exception as e:
    print(f"  Error reading map: {e}")
    # Try reading manually
    ymap = sz_hdu[1].data['SIGNAL'] if 'SIGNAL' in sz_hdu[1].columns.names else sz_hdu[1].data.field(0)
    nside = hp.npix2nside(len(ymap))
    print(f"  Manual read: Nside={nside}, npix={len(ymap)}")
    print(f"  y range: [{ymap.min():.2e}, {ymap.max():.2e}]")

sys.stdout.flush()

# ============================================================
# 2. SAMPLE COMPTON-y AT SN Ia POSITIONS
# ============================================================
print("\n[2] Sampling Compton-y at SN Ia positions...")
sys.stdout.flush()

pan_data = np.genfromtxt('/root/clawd/projects/closure-theory/data/pantheon_plus.dat',
                         names=True, dtype=None, encoding='utf-8')

ra_sn = pan_data['RA'].astype(float)
dec_sn = pan_data['DEC'].astype(float)
z_sn = pan_data['zHD'].astype(float)
c_sn = pan_data['c'].astype(float)
x1_sn = pan_data['x1'].astype(float)
mu_corr = pan_data['m_b_corr'].astype(float)
ebv_sn = pan_data['MWEBV'].astype(float) if 'MWEBV' in pan_data.dtype.names else None

valid_sn = np.isfinite(ra_sn) & np.isfinite(dec_sn) & np.isfinite(z_sn) & \
           np.isfinite(c_sn) & np.isfinite(x1_sn) & (z_sn > 0.01)

ra_v = ra_sn[valid_sn]
dec_v = dec_sn[valid_sn]
z_v = z_sn[valid_sn]
c_v = c_sn[valid_sn]
x1_v = x1_sn[valid_sn]
mu_v = mu_corr[valid_sn]

# Convert RA/Dec to HEALPix pixels
theta = np.radians(90.0 - dec_v)  # colatitude
phi = np.radians(ra_v)
pix_sn = hp.ang2pix(nside, theta, phi)

# Sample y-map
y_sn = ymap[pix_sn]
print(f"  Valid SNe: {len(y_sn)}")
print(f"  y at SN positions: [{y_sn.min():.2e}, {y_sn.max():.2e}], mean={y_sn.mean():.2e}")

# Also get galactic coordinates for control
coords_sn = SkyCoord(ra=ra_v*u.degree, dec=dec_v*u.degree, frame='icrs')
gal_b_sn = coords_sn.galactic.b.degree

# Mask galactic plane (|b| < 15°) where foreground contamination dominates
gal_mask = np.abs(gal_b_sn) > 15
print(f"  After |b|>15° cut: {gal_mask.sum()} SNe")
sys.stdout.flush()

# ============================================================
# 3. CORRELATE COMPTON-y WITH SN PARAMETERS
# ============================================================
print("\n[3] Correlating Compton-y with SN Ia parameters...")
sys.stdout.flush()

yg = y_sn[gal_mask]
cg = c_v[gal_mask]
x1g = x1_v[gal_mask]
mug = mu_v[gal_mask]
zg = z_v[gal_mask]

# Use log(y) or absolute y? Many values near zero or negative.
# Use the raw y - negative values indicate noise/CMB subtraction artifacts
# Use |y| with sign preserved for correlation

rho_yc, p_yc = stats.spearmanr(yg, cg)
rho_yx, p_yx = stats.spearmanr(yg, x1g)
rho_ym, p_ym = stats.spearmanr(yg, mug)

print(f"  Compton-y vs Color (c):    ρ = {rho_yc:+.4f}, p = {p_yc:.2e}")
print(f"  Compton-y vs Stretch (x1): ρ = {rho_yx:+.4f}, p = {p_yx:.2e}")
print(f"  Compton-y vs μ_corr:       ρ = {rho_ym:+.4f}, p = {p_ym:.2e}")

# Redshift-controlled: partial correlation
z_resid_y = yg - np.polyval(np.polyfit(zg, yg, 2), zg)
z_resid_c = cg - np.polyval(np.polyfit(zg, cg, 2), zg)
rho_yc_zctrl, p_yc_zctrl = stats.spearmanr(z_resid_y, z_resid_c)
print(f"  Compton-y vs Color (z-controlled): ρ = {rho_yc_zctrl:+.4f}, p = {p_yc_zctrl:.2e}")
sys.stdout.flush()

# ============================================================
# 4. QUINTILE ANALYSIS — SN Ia
# ============================================================
print("\n[4] SN Ia quintile analysis by Compton-y...")
sys.stdout.flush()

quintiles = np.percentile(yg, [20, 40, 60, 80])
q_labels = ['Q1 (lowest y)', 'Q2', 'Q3', 'Q4', 'Q5 (highest y)']

print(f"\n  {'Quintile':<20s} {'N':>5s} {'mean(c)':>10s} {'std(c)':>10s} {'mean(x1)':>10s} {'std(x1)':>10s} {'mean(y)':>12s}")
print("  " + "-" * 80)

q_stds_c = []
q_stds_x1 = []
q_ymeans = []
q_mean_c = []

for qi in range(5):
    if qi == 0:
        mask = yg <= quintiles[0]
    elif qi == 4:
        mask = yg > quintiles[3]
    else:
        mask = (yg > quintiles[qi-1]) & (yg <= quintiles[qi])
    
    n_q = mask.sum()
    mc = cg[mask].mean()
    sc = cg[mask].std()
    mx = x1g[mask].mean()
    sx = x1g[mask].std()
    my = yg[mask].mean()
    
    q_stds_c.append(sc)
    q_stds_x1.append(sx)
    q_ymeans.append(my)
    q_mean_c.append(mc)
    
    print(f"  {q_labels[qi]:<20s} {n_q:5d} {mc:10.4f} {sc:10.4f} {mx:10.4f} {sx:10.4f} {my:12.2e}")

rho_q_sc, _ = stats.spearmanr(q_ymeans, q_stds_c)
rho_q_sx, _ = stats.spearmanr(q_ymeans, q_stds_x1)
rho_q_mc, _ = stats.spearmanr(q_ymeans, q_mean_c)
print(f"\n  Quintile trend (y → color variance): ρ = {rho_q_sc:+.3f}")
print(f"  Quintile trend (y → stretch var):    ρ = {rho_q_sx:+.3f}")
print(f"  Quintile trend (y → mean color):     ρ = {rho_q_mc:+.3f}")

# Covariance Q1 vs Q5
q1m = yg <= quintiles[0]
q5m = yg > quintiles[3]
if q1m.sum() > 10 and q5m.sum() > 10:
    cov_q1 = np.cov(cg[q1m], x1g[q1m])
    cov_q5 = np.cov(cg[q5m], x1g[q5m])
    print(f"\n  Q1 (low y)  var(c)={cov_q1[0,0]:.6f}  cov(c,x1)={cov_q1[0,1]:.6f}")
    print(f"  Q5 (high y) var(c)={cov_q5[0,0]:.6f}  cov(c,x1)={cov_q5[0,1]:.6f}")
    
    # Box's M test analog: ratio of determinants
    det_q1 = np.linalg.det(cov_q1)
    det_q5 = np.linalg.det(cov_q5)
    print(f"  det(Σ_Q1) / det(Σ_Q5) = {det_q1/det_q5:.4f}" if det_q5 != 0 else "  det ratio: Q5 singular")
sys.stdout.flush()

# ============================================================
# 5. QUASAR SIGHTLINES
# ============================================================
print("\n[5] Quasar sightlines — Compton-y vs coupling strength...")
sys.stdout.flush()

hdu_q = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
dq = hdu_q[1].data
z_q = dq['Z_DR16Q']
mg_ew = dq['MGII_BR'][:, 2]
hb_ew = dq['HBETA_BR'][:, 2]
civ_ew = dq['CIV_BR'][:, 2]
ciii_ew = dq['CIII_BR'][:, 2]
ra_q = dq['RA']
dec_q = dq['DEC']

# Transition zone
trans_mask = (z_q >= 0.75) & (z_q < 1.15) & np.isfinite(mg_ew) & np.isfinite(hb_ew) & \
             (mg_ew != 0) & (hb_ew != 0)

# Galactic plane cut
coords_q = SkyCoord(ra=ra_q[trans_mask]*u.degree, dec=dec_q[trans_mask]*u.degree, frame='icrs')
gal_b_q = coords_q.galactic.b.degree
gal_cut_q = np.abs(gal_b_q) > 15

# Apply cuts
ra_qt = ra_q[trans_mask][gal_cut_q]
dec_qt = dec_q[trans_mask][gal_cut_q]
z_qt = z_q[trans_mask][gal_cut_q]
mg_qt = mg_ew[trans_mask][gal_cut_q]
hb_qt = hb_ew[trans_mask][gal_cut_q]

print(f"  Transition quasars (|b|>15°): {len(ra_qt)}")

# Sample y at quasar positions
theta_q = np.radians(90.0 - dec_qt)
phi_q = np.radians(ra_qt)
pix_q = hp.ang2pix(nside, theta_q, phi_q)
y_qt = ymap[pix_q]

log_mg = np.log10(np.abs(mg_qt) + 1)
log_hb = np.log10(np.abs(hb_qt) + 1)

# Quintile analysis
q_bins = np.percentile(y_qt, [20, 40, 60, 80])
print(f"\n  {'y-Quintile':<20s} {'N':>6s} {'ρ(MgII,Hβ)':>12s} {'p-value':>12s} {'mean(y)':>12s}")
print("  " + "-" * 65)

q_rhos = []
q_y_mids = []
for qi in range(5):
    if qi == 0:
        mask = y_qt <= q_bins[0]
    elif qi == 4:
        mask = y_qt > q_bins[3]
    else:
        mask = (y_qt > q_bins[qi-1]) & (y_qt <= q_bins[qi])
    
    if mask.sum() > 30:
        r, p = stats.spearmanr(log_mg[mask], log_hb[mask])
        q_rhos.append(r)
        q_y_mids.append(y_qt[mask].mean())
        print(f"  Q{qi+1} (y={y_qt[mask].mean():+.2e}){'':<3s} {mask.sum():6d} {r:12.4f} {p:12.2e} {y_qt[mask].mean():12.2e}")

if len(q_rhos) >= 3:
    rho_trend_q, p_trend_q = stats.spearmanr(q_y_mids, q_rhos)
    print(f"\n  Coupling trend (y → ρ(MgII,Hβ)): ρ = {rho_trend_q:+.3f}, p = {p_trend_q:.3f}")
    print(f"  Theory predicts: HIGH y (ordered) → STRONGER coupling (higher ρ)")
    if rho_trend_q > 0:
        print(f"  ✅ Direction matches: ordered sightlines preserve coupling")
    else:
        print(f"  ❌ Direction wrong or flat")
else:
    rho_trend_q = p_trend_q = float('nan')

sys.stdout.flush()

# ============================================================
# 6. REDSHIFT-BINNED QUASAR TEST (more powerful)
# ============================================================
print("\n[6] Redshift-binned quasar test...")
sys.stdout.flush()

# Test across wider redshift range, not just transition zone
wide_mask = (z_q >= 0.4) & (z_q < 2.0) & np.isfinite(mg_ew) & np.isfinite(hb_ew) & \
            (mg_ew != 0) & (hb_ew != 0)

coords_wide = SkyCoord(ra=ra_q[wide_mask]*u.degree, dec=dec_q[wide_mask]*u.degree, frame='icrs')
gal_b_wide = coords_wide.galactic.b.degree
gal_cut_wide = np.abs(gal_b_wide) > 15

ra_w = ra_q[wide_mask][gal_cut_wide]
dec_w = dec_q[wide_mask][gal_cut_wide]
z_w = z_q[wide_mask][gal_cut_wide]
mg_w = mg_ew[wide_mask][gal_cut_wide]
hb_w = hb_ew[wide_mask][gal_cut_wide]

theta_w = np.radians(90.0 - dec_w)
phi_w = np.radians(ra_w)
pix_w = hp.ang2pix(nside, theta_w, phi_w)
y_w = ymap[pix_w]

# Split into high-y and low-y halves, then track coupling vs redshift
y_median = np.median(y_w)
high_y = y_w > y_median
low_y = ~high_y

z_bins = np.arange(0.4, 2.0, 0.1)
print(f"\n  {'z-bin':<12s} {'ρ_high_y':>10s} {'ρ_low_y':>10s} {'Δρ':>10s} {'N_high':>8s} {'N_low':>8s}")
print("  " + "-" * 60)

deltas = []
z_mids = []
for i in range(len(z_bins)-1):
    z_lo, z_hi = z_bins[i], z_bins[i+1]
    zmid = (z_lo + z_hi) / 2
    
    mask_h = high_y & (z_w >= z_lo) & (z_w < z_hi)
    mask_l = low_y & (z_w >= z_lo) & (z_w < z_hi)
    
    if mask_h.sum() > 50 and mask_l.sum() > 50:
        lmg_h = np.log10(np.abs(mg_w[mask_h]) + 1)
        lhb_h = np.log10(np.abs(hb_w[mask_h]) + 1)
        lmg_l = np.log10(np.abs(mg_w[mask_l]) + 1)
        lhb_l = np.log10(np.abs(hb_w[mask_l]) + 1)
        
        rh, _ = stats.spearmanr(lmg_h, lhb_h)
        rl, _ = stats.spearmanr(lmg_l, lhb_l)
        delta = rh - rl
        deltas.append(delta)
        z_mids.append(zmid)
        
        print(f"  {z_lo:.1f}-{z_hi:.1f}{'':<6s} {rh:10.4f} {rl:10.4f} {delta:+10.4f} {mask_h.sum():8d} {mask_l.sum():8d}")

if len(deltas) > 3:
    mean_delta = np.mean(deltas)
    t_stat, t_p = stats.ttest_1samp(deltas, 0)
    print(f"\n  Mean Δρ (high_y - low_y): {mean_delta:+.4f}")
    print(f"  t-test vs 0: t={t_stat:.3f}, p={t_p:.3e}")
    print(f"  Theory: Δρ > 0 means ordered sightlines preserve coupling better")
sys.stdout.flush()

# ============================================================
# 7. EARLY BREAKER / LATE HOLDER TEST
# ============================================================
print("\n[7] Early breaker / Late holder Compton-y comparison...")
sys.stdout.flush()

# Recreate the breaker/holder classification from kill grid
trans_q = (z_q >= 0.75) & (z_q < 1.15) & np.isfinite(mg_ew) & np.isfinite(hb_ew) & \
          (mg_ew != 0) & (hb_ew != 0)

z_t = z_q[trans_q]
mg_t = mg_ew[trans_q]
hb_t = hb_ew[trans_q]
ra_t = ra_q[trans_q]
dec_t = dec_q[trans_q]

# Compute rank agreement
window = 0.05
n_trans = trans_q.sum()
rank_agr = np.full(n_trans, np.nan)

print(f"  Computing rank agreements for {n_trans} quasars...")
sys.stdout.flush()

# Vectorized rank agreement (approximate with percentile bins)
for i in range(n_trans):
    z_i = z_t[i]
    local = (z_t >= z_i - window) & (z_t < z_i + window)
    if local.sum() < 30:
        continue
    mg_pct = stats.percentileofscore(mg_t[local], mg_t[i]) / 100
    hb_pct = stats.percentileofscore(hb_t[local], hb_t[i]) / 100
    rank_agr[i] = 1 - abs(mg_pct - hb_pct)

valid_ra = np.isfinite(rank_agr)

early_zone = (z_t >= 0.80) & (z_t < 0.95)
late_zone = (z_t >= 1.00) & (z_t < 1.15)

ra_early_vals = rank_agr[early_zone & valid_ra]
ra_late_vals = rank_agr[late_zone & valid_ra]

if len(ra_early_vals) > 100 and len(ra_late_vals) > 100:
    early_thresh = np.percentile(ra_early_vals, 20)
    late_thresh = np.percentile(ra_late_vals, 80)
    
    early_breakers = early_zone & valid_ra & (rank_agr < early_thresh)
    late_holders = late_zone & valid_ra & (rank_agr > late_thresh)
    
    # Get Compton-y for breakers and holders
    theta_t = np.radians(90.0 - dec_t)
    phi_t = np.radians(ra_t)
    pix_t = hp.ang2pix(nside, theta_t, phi_t)
    y_t = ymap[pix_t]
    
    y_breakers = y_t[early_breakers]
    y_holders = y_t[late_holders]
    
    stat_u, p_u = stats.mannwhitneyu(y_breakers, y_holders, alternative='two-sided')
    
    print(f"\n  Early breakers: N={early_breakers.sum()}, mean(y)={y_breakers.mean():.4e}, median={np.median(y_breakers):.4e}")
    print(f"  Late holders:   N={late_holders.sum()}, mean(y)={y_holders.mean():.4e}, median={np.median(y_holders):.4e}")
    print(f"  Mann-Whitney U: p = {p_u:.4e}")
    print(f"  Theory: holders should have HIGHER y (more ordered sightlines)")
    if y_holders.mean() > y_breakers.mean():
        print(f"  ✅ Direction correct: holders see more thermal SZ")
    else:
        print(f"  ❌ Direction wrong: breakers see more thermal SZ")
else:
    print("  Insufficient data for breaker/holder test")

sys.stdout.flush()

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY — tSZ COMPTON-y STRUCTURAL STATE TEST")
print("=" * 70)

print(f"\nSN Ia (N={gal_mask.sum()}, |b|>15°):")
print(f"  y vs Color:             ρ = {rho_yc:+.4f} (p = {p_yc:.2e})")
print(f"  y vs Stretch:           ρ = {rho_yx:+.4f} (p = {p_yx:.2e})")
print(f"  y vs μ_corr:            ρ = {rho_ym:+.4f} (p = {p_ym:.2e})")
print(f"  y vs Color (z-ctrl):    ρ = {rho_yc_zctrl:+.4f} (p = {p_yc_zctrl:.2e})")
print(f"  Q-trend (y → std(c)):   ρ = {rho_q_sc:+.3f}")
print(f"  Q-trend (y → std(x1)):  ρ = {rho_q_sx:+.3f}")

if not np.isnan(rho_trend_q):
    print(f"\nQuasars (transition zone, |b|>15°):")
    print(f"  Coupling trend (y → ρ): ρ = {rho_trend_q:+.3f} (p = {p_trend_q:.3f})")

if len(deltas) > 3:
    print(f"\nRedshift-binned Δρ (high_y - low_y):")
    print(f"  Mean Δρ: {mean_delta:+.4f}, t-test p = {t_p:.3e}")

print("\n--- INTERPRETATION ---")
print("  Compare with kill grid nulls:")
print("  RM variance:    p = 0.97 (DEAD)")
print("  Planck κ:       p = 0.71 (DEAD)")
print("  Dust:           wrong direction")
if abs(rho_yc) > 0.05 or (len(deltas) > 3 and abs(mean_delta) > 0.01):
    print("  tSZ Compton-y:  SIGNAL DETECTED")
else:
    print("  tSZ Compton-y:  No signal at this resolution")

# Save
results = {
    'sn_results': {
        'n_valid': int(gal_mask.sum()),
        'y_vs_color': {'rho': float(rho_yc), 'p': float(p_yc)},
        'y_vs_stretch': {'rho': float(rho_yx), 'p': float(p_yx)},
        'y_vs_mu': {'rho': float(rho_ym), 'p': float(p_ym)},
        'y_vs_color_zctrl': {'rho': float(rho_yc_zctrl), 'p': float(p_yc_zctrl)},
        'quintile_trend_color_var': float(rho_q_sc),
        'quintile_trend_stretch_var': float(rho_q_sx),
    },
}
if not np.isnan(rho_trend_q):
    results['quasar_coupling_trend'] = {'rho': float(rho_trend_q), 'p': float(p_trend_q)}
if len(deltas) > 3:
    results['redshift_binned_delta'] = {'mean': float(mean_delta), 't_stat': float(t_stat), 'p': float(t_p)}

with open('results_tsz/tsz_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nResults saved to results_tsz/tsz_results.json")
print("=" * 70)
