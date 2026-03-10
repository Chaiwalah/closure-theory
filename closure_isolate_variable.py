#!/usr/bin/env python3
"""
ISOLATE THE VARIABLE

Stop guessing what the path variable is. Let the data tell us.

Strategy:
1. For each quasar near the transition (z=0.7-1.2), compute its 
   "local contribution" to the MgII↔Hβ correlation.
   (Does this object HELP or HURT the correlation in its z-bin?)
   
2. Identify EARLY BREAKERS (objects that decorrelate before they should)
   and LATE HOLDERS (objects that stay correlated longer than expected)
   
3. Compare EVERYTHING about these two groups:
   - Sky position (galactic lat/lon, ecliptic)
   - Source properties (mass, luminosity, Eddington)
   - RM (if in direct-match sample)
   - Emission line properties
   - Dust extinction
   
4. Whatever separates them IS the variable. 
   If nothing separates them → the effect is truly universal (no sightline component)
   If something does → we found the box
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_isolate', exist_ok=True)

print("=" * 70)
print("ISOLATE THE VARIABLE — What separates early breakers from late holders?")
print("=" * 70)

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = hdu[1].data
z_all = d['Z_DR16Q']
mg_ew = d['MGII_BR'][:, 2]
hb_ew = d['HBETA_BR'][:, 2]
ciii_ew = d['CIII_BR'][:, 2]
civ_ew = d['CIV'][:, 2]

logmbh = d['LOGMBH']
loglbol = d['LOGLBOL']
logedd = d['LOGLEDD_RATIO']
ra_q = d['RA']
dec_q = d['DEC']

# EBV if available
ebv = d['EBV'] if 'EBV' in d.dtype.names else None

coords = SkyCoord(ra=ra_q*u.degree, dec=dec_q*u.degree, frame='icrs')
gal_lat = coords.galactic.b.degree
gal_lon = coords.galactic.l.degree

# ============================================================
# Step 1: Compute per-object "correlation contribution"
# ============================================================
print("\nStep 1: Computing per-object correlation contribution...")
print("For each quasar, does it HELP or HURT the local MgII↔Hβ correlation?")

# Focus on the transition zone
trans_mask = (z_all >= 0.75) & (z_all < 1.15)
valid = trans_mask & np.isfinite(mg_ew) & np.isfinite(hb_ew) & (mg_ew != 0) & (hb_ew != 0)

print(f"Objects in transition zone (z=0.75-1.15): {valid.sum()}")

z_v = z_all[valid]
mg_v = mg_ew[valid]
hb_v = hb_ew[valid]

# Convert to ranks within narrow z-slices (removes z-dependence)
# For each object, compute its rank(MgII) and rank(Hβ) within its z-slice
# If rank(MgII) ≈ rank(Hβ), the object SUPPORTS the correlation
# If they're far apart, the object BREAKS the correlation

window = 0.05  # z-window for local ranking
rank_agreement = np.full(valid.sum(), np.nan)

for i in range(valid.sum()):
    z_i = z_v[i]
    local = (z_v >= z_i - window) & (z_v < z_i + window)
    n_local = local.sum()
    
    if n_local < 30:
        continue
    
    # Rank this object within its local z-neighborhood
    mg_local = mg_v[local]
    hb_local = hb_v[local]
    
    # Find this object's percentile in each distribution
    mg_pct = stats.percentileofscore(mg_local, mg_v[i]) / 100
    hb_pct = stats.percentileofscore(hb_local, hb_v[i]) / 100
    
    # Agreement = how close the percentiles are (1 = perfect agreement, 0 = opposite)
    rank_agreement[i] = 1 - abs(mg_pct - hb_pct)

valid_ra = np.isfinite(rank_agreement)
print(f"Objects with valid rank agreement: {valid_ra.sum()}")
print(f"Mean rank agreement: {np.nanmean(rank_agreement):.3f}")
print(f"Std: {np.nanstd(rank_agreement):.3f}")

# ============================================================
# Step 2: Identify early breakers and late holders
# ============================================================
print(f"\n{'='*70}")
print("Step 2: Identifying early breakers and late holders")
print(f"{'='*70}")

# Early breakers: low rank agreement at z < 0.95 (where correlation should still hold)
# Late holders: high rank agreement at z > 1.0 (where correlation should be dead)

early_zone = (z_v >= 0.80) & (z_v < 0.95)
late_zone = (z_v >= 1.00) & (z_v < 1.15)

# In early zone: bottom 20% of rank agreement = "early breakers"
early_valid = early_zone & valid_ra
ra_early = rank_agreement[early_valid]
early_thresh = np.percentile(ra_early, 20)

# In late zone: top 20% of rank agreement = "late holders" 
late_valid = late_zone & valid_ra
ra_late = rank_agreement[late_valid]
late_thresh = np.percentile(ra_late, 80)

early_breakers = early_zone & valid_ra & (rank_agreement <= early_thresh)
late_holders = late_zone & valid_ra & (rank_agreement >= late_thresh)
# Also get "normal" objects as control
early_normal = early_zone & valid_ra & (rank_agreement > np.percentile(ra_early, 40)) & (rank_agreement < np.percentile(ra_early, 60))
late_normal = late_zone & valid_ra & (rank_agreement > np.percentile(ra_late, 40)) & (rank_agreement < np.percentile(ra_late, 60))

print(f"Early breakers (z=0.8-0.95, bottom 20% agreement): {early_breakers.sum()}")
print(f"Early normal (z=0.8-0.95, middle 20%): {early_normal.sum()}")
print(f"Late holders (z=1.0-1.15, top 20% agreement): {late_holders.sum()}")
print(f"Late normal (z=1.0-1.15, middle 20%): {late_normal.sum()}")

# Map back to full catalog indices for property lookup
valid_indices = np.where(valid)[0]

def get_props(mask):
    """Get properties for a subset defined by mask on the valid array"""
    idx = valid_indices[mask]
    return {
        'z': z_all[idx],
        'gal_lat': gal_lat[idx],
        'gal_lon': gal_lon[idx],
        'abs_gal_lat': np.abs(gal_lat[idx]),
        'logmbh': logmbh[idx],
        'loglbol': loglbol[idx],
        'logedd': logedd[idx],
        'mg_ew': mg_ew[idx],
        'hb_ew': hb_ew[idx],
        'ciii_ew': ciii_ew[idx],
        'civ_ew': civ_ew[idx],
        'ra': ra_q[idx],
        'dec': dec_q[idx],
        'ebv': ebv[idx] if ebv is not None else None,
    }

eb = get_props(early_breakers)
en = get_props(early_normal)
lh = get_props(late_holders)
ln = get_props(late_normal)

# ============================================================
# Step 3: Compare EVERYTHING
# ============================================================
print(f"\n{'='*70}")
print("Step 3: WHAT SEPARATES THEM?")
print(f"{'='*70}")

def compare(name, arr1, arr2, label1, label2):
    """Compare two distributions, report median difference and KS test"""
    v1 = arr1[np.isfinite(arr1)]
    v2 = arr2[np.isfinite(arr2)]
    if len(v1) < 10 or len(v2) < 10:
        return None
    
    med1 = np.median(v1)
    med2 = np.median(v2)
    ks_stat, p_val = stats.ks_2samp(v1, v2)
    mwu_stat, p_mwu = stats.mannwhitneyu(v1, v2, alternative='two-sided')
    
    sig = "★★★" if p_mwu < 0.001 else "★★" if p_mwu < 0.01 else "★" if p_mwu < 0.05 else ""
    
    return {
        'name': name,
        'med1': med1, 'med2': med2,
        'delta': med2 - med1,
        'ks_p': p_val, 'mwu_p': p_mwu,
        'sig': sig,
        'n1': len(v1), 'n2': len(v2)
    }

# EARLY BREAKERS vs EARLY NORMAL (same z, what's different?)
print(f"\n  EARLY BREAKERS vs EARLY NORMAL (z=0.8-0.95)")
print(f"  What makes some objects break correlation early?")
print(f"  {'Property':<20} {'Breakers':>10} {'Normal':>10} {'Δ':>8} {'p(MWU)':>10} {'Sig':>5}")
print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*8} {'-'*10} {'-'*5}")

properties = [
    ('|gal_lat|', eb['abs_gal_lat'], en['abs_gal_lat']),
    ('gal_lon', eb['gal_lon'], en['gal_lon']),
    ('log M_BH', eb['logmbh'], en['logmbh']),
    ('log L_bol', eb['loglbol'], en['loglbol']),
    ('log L/L_Edd', eb['logedd'], en['logedd']),
    ('MgII EW', eb['mg_ew'], en['mg_ew']),
    ('Hβ EW', eb['hb_ew'], en['hb_ew']),
    ('CIII EW', eb['ciii_ew'], en['ciii_ew']),
    ('CIV EW', eb['civ_ew'], en['civ_ew']),
]

if eb['ebv'] is not None:
    properties.append(('E(B-V)', eb['ebv'], en['ebv']))

early_results = []
for name, v1, v2 in properties:
    r = compare(name, v1, v2, 'Breakers', 'Normal')
    if r:
        early_results.append(r)
        print(f"  {name:<20} {r['med1']:>10.3f} {r['med2']:>10.3f} {r['delta']:>+8.3f} {r['mwu_p']:>10.4f} {r['sig']:>5}")

# LATE HOLDERS vs LATE NORMAL (same z, what's different?)
print(f"\n  LATE HOLDERS vs LATE NORMAL (z=1.0-1.15)")
print(f"  What makes some objects hold correlation late?")
print(f"  {'Property':<20} {'Holders':>10} {'Normal':>10} {'Δ':>8} {'p(MWU)':>10} {'Sig':>5}")
print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*8} {'-'*10} {'-'*5}")

properties_late = [
    ('|gal_lat|', lh['abs_gal_lat'], ln['abs_gal_lat']),
    ('gal_lon', lh['gal_lon'], ln['gal_lon']),
    ('log M_BH', lh['logmbh'], ln['logmbh']),
    ('log L_bol', lh['loglbol'], ln['loglbol']),
    ('log L/L_Edd', lh['logedd'], ln['logedd']),
    ('MgII EW', lh['mg_ew'], ln['mg_ew']),
    ('Hβ EW', lh['hb_ew'], ln['hb_ew']),
    ('CIII EW', lh['ciii_ew'], ln['ciii_ew']),
    ('CIV EW', lh['civ_ew'], ln['civ_ew']),
]

if lh['ebv'] is not None:
    properties_late.append(('E(B-V)', lh['ebv'], ln['ebv']))

late_results = []
for name, v1, v2 in properties_late:
    r = compare(name, v1, v2, 'Holders', 'Normal')
    if r:
        late_results.append(r)
        print(f"  {name:<20} {r['med1']:>10.3f} {r['med2']:>10.3f} {r['delta']:>+8.3f} {r['mwu_p']:>10.4f} {r['sig']:>5}")

# ============================================================
# Step 4: Sky distribution — are early breakers clustered?
# ============================================================
print(f"\n{'='*70}")
print("Step 4: SKY CLUSTERING — Are early breakers from specific directions?")
print(f"{'='*70}")

# Galactic longitude distribution
for label, props in [("Early breakers", eb), ("Early normal", en)]:
    lon = props['gal_lon']
    quad_counts = [
        ((lon >= 0) & (lon < 90)).sum(),
        ((lon >= 90) & (lon < 180)).sum(),
        ((lon >= 180) & (lon < 270)).sum(),
        ((lon >= 270) & (lon < 360)).sum(),
    ]
    total = sum(quad_counts)
    pcts = [100*c/total for c in quad_counts]
    print(f"  {label}: Q1={pcts[0]:.1f}% Q2={pcts[1]:.1f}% Q3={pcts[2]:.1f}% Q4={pcts[3]:.1f}%")

# Chi-squared test for uniform distribution
eb_lon = eb['gal_lon']
en_lon = en['gal_lon']

# Bin into 8 longitude bins
bins = np.linspace(0, 360, 9)
eb_hist, _ = np.histogram(eb_lon, bins=bins)
en_hist, _ = np.histogram(en_lon, bins=bins)

# Normalize
eb_frac = eb_hist / eb_hist.sum()
en_frac = en_hist / en_hist.sum()

print(f"\n  Longitude distribution (8 bins):")
print(f"  {'l-bin':<15} {'Breakers%':>10} {'Normal%':>10} {'Δ':>8}")
for i in range(8):
    lo, hi = bins[i], bins[i+1]
    print(f"  l=[{lo:.0f},{hi:.0f})    {100*eb_frac[i]:>10.1f} {100*en_frac[i]:>10.1f} {100*(eb_frac[i]-en_frac[i]):>+8.1f}")

chi2, p_chi2 = stats.chisquare(eb_hist, f_exp=en_hist * eb_hist.sum() / en_hist.sum())
print(f"\n  χ² test (breakers vs normal longitude distribution): χ²={chi2:.1f}, p={p_chi2:.3f}")

# ============================================================  
# Step 5: The killer — COMBINED discriminant
# ============================================================
print(f"\n{'='*70}")
print("Step 5: MULTIVARIATE DISCRIMINANT")
print("Can we PREDICT early breakers from their properties?")
print(f"{'='*70}")

# Simple logistic-style: which properties best separate breakers from normal?
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

# Build feature matrix for early zone
all_early = early_breakers | early_normal
idx_all_early = valid_indices[all_early]

features = np.column_stack([
    np.abs(gal_lat[idx_all_early]),
    gal_lon[idx_all_early],
    logmbh[idx_all_early],
    loglbol[idx_all_early],
    logedd[idx_all_early],
    mg_ew[idx_all_early],
    hb_ew[idx_all_early],
])

# Handle NaNs
feature_names = ['|gal_lat|', 'gal_lon', 'log_MBH', 'log_Lbol', 'log_Edd', 'MgII_EW', 'Hβ_EW']
labels = np.concatenate([np.ones(early_breakers.sum()), np.zeros(early_normal.sum())])

# Remove rows with NaN
valid_feat = np.all(np.isfinite(features), axis=1)
features_clean = features[valid_feat]
labels_clean = labels[valid_feat]

if len(features_clean) > 50:
    try:
        rf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
        scores = cross_val_score(rf, features_clean, labels_clean, cv=5, scoring='roc_auc')
        print(f"\n  Random Forest AUC (5-fold CV): {scores.mean():.3f} ± {scores.std():.3f}")
        print(f"  (0.5 = random, 1.0 = perfect separation)")
        
        if scores.mean() > 0.55:
            # Fit on all data to get feature importances
            rf.fit(features_clean, labels_clean)
            importances = rf.feature_importances_
            
            print(f"\n  Feature importances:")
            for fname, imp in sorted(zip(feature_names, importances), key=lambda x: -x[1]):
                bar = "█" * int(imp * 50)
                print(f"    {fname:<15} {imp:.3f} {bar}")
        else:
            print(f"  → Can't separate breakers from normals using these features")
            print(f"  → The effect may be truly universal (no sightline/source predictor)")
    except Exception as e:
        print(f"  RF failed: {e}")

# Save
results = {
    'n_early_breakers': int(early_breakers.sum()),
    'n_early_normal': int(early_normal.sum()),
    'n_late_holders': int(late_holders.sum()),
    'n_late_normal': int(late_normal.sum()),
    'early_significant': [r['name'] for r in early_results if r['mwu_p'] < 0.05],
    'late_significant': [r['name'] for r in late_results if r['mwu_p'] < 0.05],
}

with open('results_isolate/isolate_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to results_isolate/")
