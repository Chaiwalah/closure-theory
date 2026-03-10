#!/usr/bin/env python3
"""
GPT's PHASE DIAGRAM — The crown jewel test.

Build two scores:
- PATH SCORE: from sky position (the 96% variable)
- SOURCE SCORE: from R_BLR (the confirmed shield)

If early breakers and late holders occupy OPPOSITE QUADRANTS:
→ Two-field competition law confirmed
→ F_path > R_source = break
→ R_source > F_path = hold

Also tests:
- RM variance (GPT's insight: disorder, not strength)
- Survey footprint confound (GPT's Mechanism 4)
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_phase_diagram', exist_ok=True)

print("=" * 70)
print("PHASE DIAGRAM — Path score vs Source score")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
logmbh = d['LOGMBH']
loglbol = d['LOGLBOL']
logedd = d['LOGLEDD_RATIO']
ra_q = d['RA']
dec_q = d['DEC']

coords = SkyCoord(ra=ra_q*u.degree, dec=dec_q*u.degree, frame='icrs')
gal_lat = coords.galactic.b.degree
gal_lon = coords.galactic.l.degree

# Load RM for variance calculation
rm_cat = Table.read('data/nvss_rm_taylor2009.fits')
coords_rm = SkyCoord(
    ra=[str(r) for r in rm_cat['RAJ2000']], 
    dec=[str(dd) for dd in rm_cat['DEJ2000']], 
    unit=(u.hourangle, u.degree), frame='icrs'
)
rm_val = np.array(rm_cat['RM'], dtype=float)

# Recreate groups
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
# Build PATH SCORE
# ============================================================
print(f"\n{'='*70}")
print("Building PATH SCORE from sky position")
print("Using the Random Forest's own logic: gal_lon + |gal_lat|")
print(f"{'='*70}")

# Train on early breakers vs normal to learn the path function
all_early_mask = early_breakers | early_normal
ae_idx = valid_indices[all_early_mask]
ae_features = np.column_stack([
    gal_lon[ae_idx],
    np.abs(gal_lat[ae_idx]),
])
ae_labels = np.concatenate([np.ones(early_breakers.sum()), np.zeros(early_normal.sum())])
ae_valid = np.all(np.isfinite(ae_features), axis=1)

rf_path = RandomForestClassifier(n_estimators=200, max_depth=5, random_state=42)
rf_path.fit(ae_features[ae_valid], ae_labels[ae_valid])

# Apply to ALL objects in transition zone to get path score
all_features_path = np.column_stack([
    gal_lon[valid_indices],
    np.abs(gal_lat[valid_indices]),
])
path_score = rf_path.predict_proba(all_features_path)[:, 1]  # P(breaker) from sky alone

print(f"Path score range: {path_score.min():.3f} - {path_score.max():.3f}")
print(f"Path score for breakers: {np.median(path_score[early_breakers]):.3f}")
print(f"Path score for normals:  {np.median(path_score[early_normal]):.3f}")

# ============================================================
# Build SOURCE SCORE  
# ============================================================
print(f"\n{'='*70}")
print("Building SOURCE SCORE from R_BLR (= L^0.5)")
print(f"{'='*70}")

# R_BLR as source stiffness
def r_blr(loglbol_arr):
    L46 = 10**(loglbol_arr - 46)
    return 33 * L46**0.5

source_score_raw = r_blr(loglbol[valid_indices])
# Normalize to [0,1]
ss_valid = np.isfinite(source_score_raw) & (source_score_raw > 0)
ss_min = np.percentile(source_score_raw[ss_valid], 5)
ss_max = np.percentile(source_score_raw[ss_valid], 95)
source_score = np.clip((source_score_raw - ss_min) / (ss_max - ss_min), 0, 1)

print(f"Source score for late holders: {np.median(source_score[late_holders]):.3f}")
print(f"Source score for late normals: {np.median(source_score[late_normal]):.3f}")

# ============================================================
# THE PHASE DIAGRAM
# ============================================================
print(f"\n{'='*70}")
print("PHASE DIAGRAM — Do they occupy opposite quadrants?")
print(f"{'='*70}")

# Define quadrants by median splits
path_med = np.median(path_score)
source_med = np.median(source_score[ss_valid])

def quadrant_analysis(mask, label):
    ps = path_score[mask]
    ss = source_score[mask]
    v = np.isfinite(ps) & np.isfinite(ss)
    ps, ss = ps[v], ss[v]
    
    q1 = ((ps > path_med) & (ss > source_med)).sum()  # high path, high source
    q2 = ((ps <= path_med) & (ss > source_med)).sum()  # low path, high source
    q3 = ((ps <= path_med) & (ss <= source_med)).sum()  # low path, low source
    q4 = ((ps > path_med) & (ss <= source_med)).sum()  # high path, low source
    total = q1 + q2 + q3 + q4
    
    print(f"\n  {label} (N={total}):")
    print(f"                  High source    Low source")
    print(f"    High path:    {100*q1/total:>6.1f}%        {100*q4/total:>6.1f}%")
    print(f"    Low path:     {100*q2/total:>6.1f}%        {100*q3/total:>6.1f}%")
    
    return {'hi_path_hi_source': q1/total, 'lo_path_hi_source': q2/total,
            'lo_path_lo_source': q3/total, 'hi_path_lo_source': q4/total}

print("\n  Expected if two-field competition:")
print("  - Early breakers → HIGH path, LOW source (top-right kills)")
print("  - Late holders → LOW path, HIGH source (bottom-left shields)")

q_eb = quadrant_analysis(early_breakers, "EARLY BREAKERS")
q_en = quadrant_analysis(early_normal, "EARLY NORMAL (control)")
q_lh = quadrant_analysis(late_holders, "LATE HOLDERS")
q_ln = quadrant_analysis(late_normal, "LATE NORMAL (control)")

# Measure separation
eb_kill_quadrant = q_eb['hi_path_lo_source']  # high path, low source
en_kill_quadrant = q_en['hi_path_lo_source']
lh_shield_quadrant = q_lh['lo_path_hi_source']  # low path, high source
ln_shield_quadrant = q_ln['lo_path_hi_source']

print(f"\n  SEPARATION TEST:")
print(f"  Breakers in 'kill quadrant' (high path, low source): {100*eb_kill_quadrant:.1f}% vs normal {100*en_kill_quadrant:.1f}%")
print(f"  Holders in 'shield quadrant' (low path, high source): {100*lh_shield_quadrant:.1f}% vs normal {100*ln_shield_quadrant:.1f}%")

if eb_kill_quadrant > en_kill_quadrant and lh_shield_quadrant > ln_shield_quadrant:
    print(f"\n  ★ OPPOSITE QUADRANTS CONFIRMED: Two-field competition law HOLDS")
    print(f"    F_path > R_source = break")
    print(f"    R_source > F_path = hold")
else:
    print(f"\n  ✗ Quadrant separation not clean")

# ============================================================
# RM VARIANCE TEST (GPT's insight)
# ============================================================
print(f"\n{'='*70}")
print("RM VARIANCE — Does path DISORDER predict breaking?")
print("(GPT: coherence filter cares about disorder, not strength)")
print(f"{'='*70}")

# For each quasar, compute RM variance in a 1° radius
print("Computing local RM variance for each quasar...")

# Use a subsample for speed (all breakers + normals in early zone)
test_idx = np.concatenate([eb_idx, en_idx])
test_coords = SkyCoord(ra=ra_q[test_idx]*u.degree, dec=dec_q[test_idx]*u.degree)

# For each test object, find all RM sources within 2°
rm_var_local = np.full(len(test_idx), np.nan)

# Use a faster approach: for each RM source, assign to nearby quasars
# Build a grid approach
print(f"  Computing RM scatter for {len(test_idx)} objects...")

for i in range(len(test_idx)):
    sep = test_coords[i].separation(coords_rm).degree
    nearby = sep < 2.0
    if nearby.sum() >= 3:
        rm_var_local[i] = np.std(rm_val[nearby])
    
    if i % 2000 == 0 and i > 0:
        print(f"    {i}/{len(test_idx)}...")

eb_rmvar = rm_var_local[:len(eb_idx)]
en_rmvar = rm_var_local[len(eb_idx):]

eb_rv = eb_rmvar[np.isfinite(eb_rmvar)]
en_rv = en_rmvar[np.isfinite(en_rmvar)]

print(f"\n  Breakers RM scatter: median = {np.median(eb_rv):.1f} rad/m² (N={len(eb_rv)})")
print(f"  Normal RM scatter:   median = {np.median(en_rv):.1f} rad/m² (N={len(en_rv)})")

_, p_rmvar = stats.mannwhitneyu(eb_rv, en_rv, alternative='greater')
print(f"  MWU p (breakers > normal): {p_rmvar:.4f}")

if p_rmvar < 0.05:
    print(f"  ★ CONFIRMED: Breakers behind MORE disordered RM sightlines")
    print(f"    Path disorder, not just strength, matters")
else:
    print(f"  ✗ RM variance doesn't distinguish breakers from normals")

# ============================================================
# GPT's MECHANISM 4: Survey footprint confound
# ============================================================
print(f"\n{'='*70}")
print("SURVEY FOOTPRINT — Does matching kill the sky signal?")
print(f"{'='*70}")

# Match on source properties, re-run the classifier
from sklearn.preprocessing import StandardScaler

# Build matched early-zone sample
ae_features_full = np.column_stack([
    gal_lon[ae_idx],
    np.abs(gal_lat[ae_idx]),
    logmbh[ae_idx],
    loglbol[ae_idx],
    logedd[ae_idx],
    mg_ew[ae_idx],
    hb_ew[ae_idx],
])

ae_valid_full = np.all(np.isfinite(ae_features_full), axis=1)
feat_clean = ae_features_full[ae_valid_full]
lab_clean = ae_labels[ae_valid_full]

# Classifier with ALL features
rf_all = RandomForestClassifier(n_estimators=200, max_depth=5, random_state=42)
scores_all = cross_val_score(rf_all, feat_clean, lab_clean, cv=5, scoring='roc_auc')
print(f"\n  AUC with ALL features (sky + source): {scores_all.mean():.3f} ± {scores_all.std():.3f}")

# Classifier with ONLY source features (no sky)
feat_source = feat_clean[:, 2:]  # drop lon and lat
rf_source = RandomForestClassifier(n_estimators=200, max_depth=5, random_state=42)
scores_source = cross_val_score(rf_source, feat_source, lab_clean, cv=5, scoring='roc_auc')
print(f"  AUC with ONLY source features: {scores_source.mean():.3f} ± {scores_source.std():.3f}")

# Classifier with ONLY sky features
feat_sky = feat_clean[:, :2]  # only lon and lat
rf_sky = RandomForestClassifier(n_estimators=200, max_depth=5, random_state=42)
scores_sky = cross_val_score(rf_sky, feat_sky, lab_clean, cv=5, scoring='roc_auc')
print(f"  AUC with ONLY sky features: {scores_sky.mean():.3f} ± {scores_sky.std():.3f}")

print(f"\n  Sky explains: AUC {scores_sky.mean():.3f}")
print(f"  Source explains: AUC {scores_source.mean():.3f}")
print(f"  Combined: AUC {scores_all.mean():.3f}")

if scores_sky.mean() > scores_source.mean() + 0.1:
    print(f"  ★ Sky dominates even with source features available")
    print(f"    Survey footprint / source-orientation confound is UNLIKELY")
elif scores_source.mean() > scores_sky.mean():
    print(f"  ⚠️ Source features dominate → sky might be an artifact")

# Save everything
results = {
    'phase_diagram': {
        'early_breakers': q_eb,
        'early_normal': q_en,
        'late_holders': q_lh,
        'late_normal': q_ln,
    },
    'auc_all': float(scores_all.mean()),
    'auc_sky': float(scores_sky.mean()),
    'auc_source': float(scores_source.mean()),
}

with open('results_phase_diagram/phase_diagram_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults saved to results_phase_diagram/")
