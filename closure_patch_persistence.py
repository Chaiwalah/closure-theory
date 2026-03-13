#!/usr/bin/env python3
"""P1: Patch Persistence Test — Bootstrap + HEALPix Jackknife"""
import numpy as np
from scipy import stats
import json, os, sys, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_patch', exist_ok=True)

print("="*60)
print("P1: PATCH PERSISTENCE TEST")
print("="*60)
sys.stdout.flush()

# 1. Load
print("\n[1] Loading Pantheon+ data...")
sys.stdout.flush()
d = np.genfromtxt('data/pantheon_plus.dat', names=True, dtype=None, encoding='utf-8')
ra = d['RA'].astype(float); dec = d['DEC'].astype(float)
z = d['zHD'].astype(float); c = d['c'].astype(float)
x1 = d['x1'].astype(float)
valid = np.isfinite(ra)&np.isfinite(dec)&np.isfinite(z)&np.isfinite(c)&(z>0.01)
ra,dec,z,c,x1 = ra[valid],dec[valid],z[valid],c[valid],x1[valid]
N = len(ra)
print(f"  {N} valid SNe")

# 2. Breaker classification
print("\n[2] Classifying breakers (z-matched, |Δc|>1σ)...")
sys.stdout.flush()
is_breaker = np.zeros(N, dtype=bool)
for i in range(N):
    local = np.abs(z - z[i]) < 0.05
    if local.sum() < 10: continue
    mu_c, sig_c = c[local].mean(), c[local].std()
    if sig_c > 0 and np.abs(c[i] - mu_c) > sig_c:
        is_breaker[i] = True
print(f"  Breakers: {is_breaker.sum()}/{N} ({100*is_breaker.mean():.1f}%)")

# 3. HEALPix setup
print("\n[3] Building HEALPix maps...")
sys.stdout.flush()
try:
    import healpy as hp
except ImportError:
    os.system('pip install healpy -q')
    import healpy as hp

NSIDE = 8
NPIX = hp.nside2npix(NSIDE)
theta = np.radians(90 - dec)
phi = np.radians(ra)
pix = hp.ang2pix(NSIDE, theta, phi)

def breaker_map(is_b, px):
    count = np.zeros(NPIX)
    bcount = np.zeros(NPIX)
    for p, b in zip(px, is_b):
        count[p] += 1
        bcount[p] += b
    frac = np.full(NPIX, np.nan)
    good = count >= 3
    frac[good] = bcount[good] / count[good]
    return frac, good, count

frac0, good0, cnt0 = breaker_map(is_breaker, pix)
n_good = good0.sum()
print(f"  Nside={NSIDE}, {NPIX} pixels, {n_good} with ≥3 SNe")
print(f"  Breaker fraction range: [{np.nanmin(frac0):.3f}, {np.nanmax(frac0):.3f}]")

# Get top-10 hotspot pixels
frac_sorted = np.argsort(np.where(np.isnan(frac0), -1, frac0))[::-1]
top10_ref = set(frac_sorted[:10])
print(f"  Top-10 hotspot pixels: {sorted(top10_ref)}")

# 4. Bootstrap
NBOOT = 1000
print(f"\n[4] Bootstrap ({NBOOT} iterations)...")
sys.stdout.flush()

top10_counts = np.zeros(NPIX, dtype=int)  # how often each pixel is in top-10
centroids_ra = []
centroids_dec = []

for b in range(NBOOT):
    if b % 200 == 0:
        print(f"  ... iteration {b}/{NBOOT}")
        sys.stdout.flush()
    idx = np.random.choice(N, N, replace=True)
    fb, gb, _ = breaker_map(is_breaker[idx], pix[idx])
    # top-10 for this bootstrap
    fs = np.argsort(np.where(np.isnan(fb), -1, fb))[::-1]
    t10 = set(fs[:10])
    for p in t10:
        top10_counts[p] += 1
    # centroid of top-10
    t10_list = list(t10)
    t10_th, t10_ph = hp.pix2ang(NSIDE, t10_list)
    t10_frac = np.array([fb[p] if not np.isnan(fb[p]) else 0 for p in t10_list])
    w = t10_frac / (t10_frac.sum() + 1e-30)
    c_th = np.average(t10_th, weights=w)
    c_ph = np.average(t10_ph, weights=w)
    centroids_ra.append(np.degrees(c_ph))
    centroids_dec.append(90 - np.degrees(c_th))

print(f"  ... done")
sys.stdout.flush()

# 5. Hotspot overlap
print("\n[5] Hotspot overlap analysis...")
sys.stdout.flush()
# What fraction of the REFERENCE top-10 appear in >50% of bootstraps?
overlap_count = sum(1 for p in top10_ref if top10_counts[p] > NBOOT * 0.5)
overlap_frac = overlap_count / 10
print(f"  Reference top-10 pixels appearing in >50% of bootstraps: {overlap_count}/10")
print(f"  Overlap fraction: {overlap_frac:.2f}")

# Also: how many unique pixels EVER appear in top-10?
ever_top10 = (top10_counts > 0).sum()
print(f"  Unique pixels ever in top-10: {ever_top10}/{NPIX}")

# 6. Centroid stability
print("\n[6] Centroid stability...")
sys.stdout.flush()
c_ra = np.array(centroids_ra)
c_dec = np.array(centroids_dec)
ra_std = np.std(c_ra)
dec_std = np.std(c_dec)
print(f"  Centroid RA std:  {ra_std:.2f}°")
print(f"  Centroid DEC std: {dec_std:.2f}°")
print(f"  Mean centroid: ({np.mean(c_ra):.1f}°, {np.mean(c_dec):.1f}°)")

# 7. Jackknife
print("\n[7] HEALPix jackknife (removing one pixel at a time)...")
sys.stdout.flush()

# Global chi2: sum of (observed_breakers - expected)^2 / expected per pixel
def compute_chi2(is_b, px, good_mask):
    global_rate = is_b.sum() / len(is_b)
    chi2 = 0
    npix_used = 0
    for p in range(NPIX):
        if not good_mask[p]: continue
        mask = px == p
        n_p = mask.sum()
        if n_p < 3: continue
        obs = is_b[mask].sum()
        exp = global_rate * n_p
        if exp > 0:
            chi2 += (obs - exp)**2 / exp
            npix_used += 1
    return chi2, npix_used

chi2_full, npix_full = compute_chi2(is_breaker, pix, good0)
print(f"  Full χ² = {chi2_full:.2f} ({npix_full} pixels)")

jackknife_chi2 = []
for p in range(NPIX):
    if not good0[p]: continue
    mask = pix != p
    fj, gj, _ = breaker_map(is_breaker[mask], pix[mask])
    c2, _ = compute_chi2(is_breaker[mask], pix[mask], gj)
    jackknife_chi2.append((p, c2))

jk_vals = [x[1] for x in jackknife_chi2]
jk_min, jk_max = min(jk_vals), max(jk_vals)
max_drop = (chi2_full - jk_min) / chi2_full * 100
worst_pixel = jackknife_chi2[np.argmin(jk_vals)][0]

print(f"  Jackknife χ² range: [{jk_min:.2f}, {jk_max:.2f}]")
print(f"  Max drop: {max_drop:.1f}% (removing pixel {worst_pixel})")
print(f"  Worst pixel location: {hp.pix2ang(NSIDE, worst_pixel, lonlat=True)}")

# VERDICT
print("\n" + "="*60)
stable = overlap_frac >= 0.5 and ra_std < 20 and dec_std < 20
unstable = overlap_frac < 0.3 or ra_std > 30 or dec_std > 30
artifact = max_drop > 50

if artifact:
    verdict = "SINGLE-PIXEL ARTIFACT"
elif unstable:
    verdict = "UNSTABLE"
elif stable:
    verdict = "STABLE"
else:
    verdict = "MIXED"

print(f"VERDICT: {verdict}")
print(f"  Overlap: {overlap_frac:.2f} (threshold: >0.5=stable, <0.3=unstable)")
print(f"  Centroid RA std: {ra_std:.1f}° (threshold: <20°=stable, >30°=unstable)")
print(f"  Centroid DEC std: {dec_std:.1f}° (threshold: <20°=stable, >30°=unstable)")
print(f"  Max jackknife drop: {max_drop:.1f}% (threshold: >50%=artifact)")
print("="*60)

results = {
    'n_sne': int(N),
    'n_breakers': int(is_breaker.sum()),
    'nside': NSIDE,
    'n_good_pixels': int(n_good),
    'n_bootstraps': NBOOT,
    'hotspot_overlap_fraction': round(overlap_frac, 4),
    'centroid_ra_std_deg': round(ra_std, 2),
    'centroid_dec_std_deg': round(dec_std, 2),
    'centroid_mean_ra': round(float(np.mean(c_ra)), 2),
    'centroid_mean_dec': round(float(np.mean(c_dec)), 2),
    'chi2_full': round(chi2_full, 2),
    'jackknife_chi2_range': [round(jk_min, 2), round(jk_max, 2)],
    'jackknife_max_drop_pct': round(max_drop, 2),
    'worst_pixel': int(worst_pixel),
    'ever_in_top10': int(ever_top10),
    'verdict': verdict
}

with open('results_patch/patch_persistence.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved to results_patch/patch_persistence.json")
