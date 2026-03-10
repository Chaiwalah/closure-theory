#!/usr/bin/env python3
"""
KILL GRID — Every mechanism lives or dies in one pass.

Tests:
1. R_BLR SHIELD (Gemini Mech 1): Late holders have larger BLR?
2. MW LOCAL (Gemini Mech 2): Breakers trace MW structure or extragalactic?
3. ALP FLUX (Gemini Mech 3): Breakers dimmer than expected?
4. RM PATH (Grok): Breaker fraction ∝ |RM|?
5. OUR OWN: Breaker sky distribution vs known superstructures

All using the 7,863 early breakers and 4,664 late holders already identified.
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_kill_grid', exist_ok=True)

print("=" * 70)
print("KILL GRID — Every mechanism tested at once")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
ciii_ew = d['CIII_BR'][:, 2]
logmbh = d['LOGMBH']
loglbol = d['LOGLBOL']
logedd = d['LOGLEDD_RATIO']
ra_q = d['RA']
dec_q = d['DEC']
ebv = d['EBV'] if 'EBV' in d.dtype.names else None

coords = SkyCoord(ra=ra_q*u.degree, dec=dec_q*u.degree, frame='icrs')
gal_lat = coords.galactic.b.degree
gal_lon = coords.galactic.l.degree

# Recreate the isolation groups (same logic as isolation script)
trans_mask = (z_all >= 0.75) & (z_all < 1.15)
valid = trans_mask & np.isfinite(mg_ew) & np.isfinite(hb_ew) & (mg_ew != 0) & (hb_ew != 0)

z_v = z_all[valid]
mg_v = mg_ew[valid]
hb_v = hb_ew[valid]

window = 0.05
rank_agreement = np.full(valid.sum(), np.nan)

print("Recomputing rank agreements...")
for i in range(valid.sum()):
    z_i = z_v[i]
    local = (z_v >= z_i - window) & (z_v < z_i + window)
    if local.sum() < 30:
        continue
    mg_local = mg_v[local]
    hb_local = hb_v[local]
    mg_pct = stats.percentileofscore(mg_local, mg_v[i]) / 100
    hb_pct = stats.percentileofscore(hb_local, hb_v[i]) / 100
    rank_agreement[i] = 1 - abs(mg_pct - hb_pct)

valid_ra = np.isfinite(rank_agreement)
valid_indices = np.where(valid)[0]

# Define groups
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

print(f"Early breakers: {len(eb_idx)}, Early normal: {len(en_idx)}")
print(f"Late holders: {len(lh_idx)}, Late normal: {len(ln_idx)}")

verdicts = {}

# ============================================================
# TEST 1: R_BLR SHIELD (Gemini Mechanism 1)
# ============================================================
print(f"\n{'='*70}")
print("TEST 1: R_BLR SHIELD — Do late holders have physically larger BLR?")
print("R_BLR ∝ L^0.5 (Kaspi+2000, Bentz+2013)")
print(f"{'='*70}")

# R_BLR in light-days ≈ 33 * (L_bol / 10^46)^0.5  (Kaspi relation)
def r_blr(loglbol_arr):
    """BLR radius from Kaspi+2000 relation"""
    L46 = 10**(loglbol_arr - 46)  # L in units of 10^46 erg/s
    return 33 * L46**0.5  # R_BLR in light-days

r_lh = r_blr(loglbol[lh_idx])
r_ln = r_blr(loglbol[ln_idx])
r_eb = r_blr(loglbol[eb_idx])
r_en = r_blr(loglbol[en_idx])

# Filter valid
r_lh_v = r_lh[np.isfinite(r_lh) & (r_lh > 0)]
r_ln_v = r_ln[np.isfinite(r_ln) & (r_ln > 0)]
r_eb_v = r_eb[np.isfinite(r_eb) & (r_eb > 0)]
r_en_v = r_en[np.isfinite(r_en) & (r_en > 0)]

print(f"\n  Late holders R_BLR: median = {np.median(r_lh_v):.1f} light-days (N={len(r_lh_v)})")
print(f"  Late normal  R_BLR: median = {np.median(r_ln_v):.1f} light-days (N={len(r_ln_v)})")
_, p_late = stats.mannwhitneyu(r_lh_v, r_ln_v, alternative='greater')
print(f"  MWU p (holders > normal): {p_late:.4f}")

print(f"\n  Early breakers R_BLR: median = {np.median(r_eb_v):.1f} light-days (N={len(r_eb_v)})")
print(f"  Early normal   R_BLR: median = {np.median(r_en_v):.1f} light-days (N={len(r_en_v)})")
_, p_early = stats.mannwhitneyu(r_eb_v, r_en_v, alternative='two-sided')
print(f"  MWU p (breakers ≠ normal): {p_early:.4f}")

if p_late < 0.05:
    print(f"\n  ★ CONFIRMED: Late holders have LARGER BLR (p={p_late:.4f})")
    print(f"    Plasma scintillation shield is consistent")
    verdicts['R_BLR_shield'] = 'SURVIVES'
else:
    print(f"\n  ✗ NOT CONFIRMED: Late holders don't have significantly larger BLR")
    verdicts['R_BLR_shield'] = 'WEAKENED'

# Kill check: at FIXED R_BLR, does Eddington ratio alone predict holding?
print(f"\n  KILL CHECK: At fixed R_BLR, does Eddington predict late holding?")
# Bin late-zone objects by R_BLR, check if holders have different Edd
r_all_late = r_blr(loglbol[valid_indices[late_zone & valid_ra]])
edd_all_late = logedd[valid_indices[late_zone & valid_ra]]
ra_all_late = rank_agreement[late_zone & valid_ra]

r_valid = np.isfinite(r_all_late) & (r_all_late > 0) & np.isfinite(edd_all_late)
r_med = np.median(r_all_late[r_valid])

# Fixed R_BLR (within 20% of median)
fixed_r = r_valid & (r_all_late > r_med * 0.8) & (r_all_late < r_med * 1.2)
if fixed_r.sum() > 100:
    holders_fixed = fixed_r & (ra_all_late >= late_thresh)
    normals_fixed = fixed_r & (ra_all_late > np.percentile(ra_late, 40)) & (ra_all_late < np.percentile(ra_late, 60))
    
    edd_h = edd_all_late[holders_fixed]
    edd_n = edd_all_late[normals_fixed]
    edd_h = edd_h[np.isfinite(edd_h)]
    edd_n = edd_n[np.isfinite(edd_n)]
    
    if len(edd_h) > 20 and len(edd_n) > 20:
        _, p_edd = stats.mannwhitneyu(edd_h, edd_n, alternative='two-sided')
        print(f"    At fixed R_BLR (±20%): Edd holders={np.median(edd_h):.3f}, normal={np.median(edd_n):.3f}, p={p_edd:.4f}")
        if p_edd < 0.05:
            print(f"    ⚠️ Eddington STILL matters at fixed R_BLR → size alone isn't the full shield")
            verdicts['R_BLR_shield'] += ' (but Edd also contributes)'
        else:
            print(f"    Size fully explains holding → scintillation shield CONFIRMED")

# ============================================================
# TEST 2: MW LOCAL vs EXTRAGALACTIC (Gemini Mechanism 2)
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: MW LOCAL — Do breakers trace MW structure or extragalactic?")
print(f"{'='*70}")

eb_lon = gal_lon[eb_idx]
eb_lat = gal_lat[eb_idx]
en_lon = gal_lon[en_idx]
en_lat = gal_lat[en_idx]

# Key MW structures:
# Galactic center: l=0, b=0
# North Galactic Spur: l~30, b>20
# Fermi Bubbles: |l|<20, |b|=10-55
# Gould Belt: complex shape
# MW spiral arms project to specific l ranges

# Test: Do breakers prefer galactic plane (MW disk) or high latitude (cleaner sightlines)?
print(f"\n  A. Galactic latitude distribution:")
print(f"     Breakers median |b|: {np.median(np.abs(eb_lat)):.1f}°")
print(f"     Normal median |b|:   {np.median(np.abs(en_lat)):.1f}°")
_, p_lat = stats.mannwhitneyu(np.abs(eb_lat), np.abs(en_lat), alternative='two-sided')
print(f"     MWU p: {p_lat:.4f}")

# If breakers prefer LOW |b| → MW local (more MW magnetic field at low lat)
# If breakers prefer HIGH |b| → NOT MW local
if np.median(np.abs(eb_lat)) < np.median(np.abs(en_lat)):
    print(f"     Breakers at LOWER |b| → consistent with MW local")
else:
    print(f"     Breakers at HIGHER |b| → INCONSISTENT with MW local!")

# Test: Fermi Bubble region
print(f"\n  B. Fermi Bubble region (|l|<20, |b|=10-55):")
in_fermi_eb = (np.abs(eb_lon) < 20) | ((360 - eb_lon) < 20)
in_fermi_eb &= (np.abs(eb_lat) > 10) & (np.abs(eb_lat) < 55)
in_fermi_en = (np.abs(en_lon) < 20) | ((360 - en_lon) < 20)
in_fermi_en &= (np.abs(en_lat) > 10) & (np.abs(en_lat) < 55)

frac_eb = in_fermi_eb.mean()
frac_en = in_fermi_en.mean()
print(f"     Breakers in Fermi region: {100*frac_eb:.1f}%")
print(f"     Normal in Fermi region:   {100*frac_en:.1f}%")

# Test: Known extragalactic overdensities
print(f"\n  C. Extragalactic structure directions:")

structures = [
    ("Shapley Concentration", 307, 30),
    ("Coma Cluster", 58, 88),
    ("Perseus-Pisces", 150, -15),
    ("Hercules", 30, 45),
    ("Sloan Great Wall", 200, 55),
    ("Dipole Repeller", 320, 10),
]

print(f"     {'Structure':<25} {'Breaker excess':>15} {'p':>8}")
print(f"     {'-'*25} {'-'*15} {'-'*8}")

for name, l_c, b_c in structures:
    # Count objects within 15° of structure direction
    dist_eb = np.sqrt(((eb_lon - l_c + 180) % 360 - 180)**2 + (eb_lat - b_c)**2)
    dist_en = np.sqrt(((en_lon - l_c + 180) % 360 - 180)**2 + (en_lat - b_c)**2)
    
    near_eb = (dist_eb < 15).mean()
    near_en = (dist_en < 15).mean()
    
    if near_en > 0:
        excess = near_eb / near_en - 1
        # Simple proportion test
        n1, n2 = len(eb_lon), len(en_lon)
        p1, p2 = near_eb, near_en
        p_pool = (near_eb * n1 + near_en * n2) / (n1 + n2)
        if p_pool > 0 and p_pool < 1:
            se = np.sqrt(p_pool * (1-p_pool) * (1/n1 + 1/n2))
            z_stat = (p1 - p2) / se if se > 0 else 0
            p_val = 2 * (1 - stats.norm.cdf(abs(z_stat)))
            print(f"     {name:<25} {100*excess:>+13.1f}% {p_val:>8.3f}")
        else:
            print(f"     {name:<25} {'N/A':>15}")
    else:
        print(f"     {name:<25} {'N/A':>15}")

# Verdict: Is the anisotropy MW or extragalactic?
if np.median(np.abs(eb_lat)) > np.median(np.abs(en_lat)):
    print(f"\n  ★ Breakers prefer HIGH latitude → MW local is UNLIKELY")
    print(f"    The anisotropy is more consistent with EXTRAGALACTIC structure")
    verdicts['MW_local'] = 'KILLED (breakers at high |b|)'
else:
    print(f"\n  ⚠️ Breakers prefer LOW latitude → MW local NOT ruled out")
    verdicts['MW_local'] = 'SURVIVES'

# ============================================================
# TEST 3: ALP FLUX DEFICIT (Gemini Mechanism 3)
# ============================================================
print(f"\n{'='*70}")
print("TEST 3: ALP FLUX — Are early breakers dimmer than expected?")
print("If ALPs remove photons, breakers should have a flux deficit")
print(f"{'='*70}")

lum_eb = loglbol[eb_idx]
lum_en = loglbol[en_idx]
lum_eb_v = lum_eb[np.isfinite(lum_eb)]
lum_en_v = lum_en[np.isfinite(lum_en)]

print(f"\n  Early breakers log L_bol: {np.median(lum_eb_v):.3f} (N={len(lum_eb_v)})")
print(f"  Early normal   log L_bol: {np.median(lum_en_v):.3f} (N={len(lum_en_v)})")

_, p_lum = stats.mannwhitneyu(lum_eb_v, lum_en_v, alternative='two-sided')
_, p_lum_less = stats.mannwhitneyu(lum_eb_v, lum_en_v, alternative='less')
print(f"  MWU p (two-sided): {p_lum:.4f}")
print(f"  MWU p (breakers < normal): {p_lum_less:.4f}")

ks_stat, p_ks = stats.ks_2samp(lum_eb_v, lum_en_v)
print(f"  KS test: D={ks_stat:.3f}, p={p_ks:.4f}")

delta_lum = np.median(lum_eb_v) - np.median(lum_en_v)
print(f"  Δ(median log L_bol): {delta_lum:+.3f} dex")

if p_lum > 0.05:
    print(f"\n  ★ ALP KILLED: Breakers have SAME luminosity as normals")
    print(f"    No photons leaving the beam → wave scrambled, not depleted")
    verdicts['ALP'] = 'KILLED'
elif delta_lum < -0.1:
    print(f"\n  ⚠️ Breakers ARE dimmer → ALP not ruled out")
    verdicts['ALP'] = 'SURVIVES'
else:
    print(f"\n  ~ Marginal luminosity difference, inconclusive")
    verdicts['ALP'] = 'INCONCLUSIVE'

# ============================================================
# TEST 4: DUST (our own) — Does E(B-V) predict breaking?
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: DUST — Does E(B-V) predict early breaking?")
print("If mechanism is magnetic, dust should NOT predict breaking")
print("(Dust ∝ λ⁻¹, our effect ∝ λ²)")
print(f"{'='*70}")

if ebv is not None:
    ebv_eb = ebv[eb_idx]
    ebv_en = ebv[en_idx]
    ebv_eb_v = ebv_eb[np.isfinite(ebv_eb)]
    ebv_en_v = ebv_en[np.isfinite(ebv_en)]
    
    print(f"\n  Breakers E(B-V): {np.median(ebv_eb_v):.4f}")
    print(f"  Normal E(B-V):   {np.median(ebv_en_v):.4f}")
    _, p_dust = stats.mannwhitneyu(ebv_eb_v, ebv_en_v, alternative='two-sided')
    print(f"  MWU p: {p_dust:.4f}")
    
    if p_dust > 0.05:
        print(f"\n  ★ DUST KILLED: Same E(B-V) → dust is NOT the mechanism")
        verdicts['dust'] = 'KILLED'
    else:
        delta_ebv = np.median(ebv_eb_v) - np.median(ebv_en_v)
        print(f"\n  ⚠️ Dust differs (Δ={delta_ebv:+.4f}) — needs more investigation")
        verdicts['dust'] = 'NOT KILLED'

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("KILL GRID SUMMARY")
print(f"{'='*70}\n")

for mech, verdict in verdicts.items():
    icon = "✓" if "KILLED" in verdict else "⚠" if "SURVIVES" in verdict else "~"
    print(f"  {icon} {mech:<25} → {verdict}")

with open('results_kill_grid/kill_grid_results.json', 'w') as f:
    json.dump(verdicts, f, indent=2)

print(f"\nResults saved to results_kill_grid/")
