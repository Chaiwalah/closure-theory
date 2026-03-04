#!/usr/bin/env python3
"""
KILL B2 — Environment-Controlled Volume-Limited Test
=====================================================
B2 says: high-z quasars live in denser environments, which changes their
intrinsic line properties AND correlates with sightline density. So the
"foreground effect" is really just environment-driven source evolution.

TO KILL B2 we need:
1. Volume-limited subsample (fixed luminosity range across z)
2. Environment control (match on environment proxies)
3. Show the EW-FWHM asymmetry STILL depends on sightline properties

If the effect persists after controlling for everything B2 could hide behind,
B2 is dead.

APPROACH:
- Use LOGLBOL as luminosity control (volume-limited)
- Use E(B-V) as sightline proxy (galactic structure tracer)
- Use galactic latitude as independent sightline proxy
- Match on SNR to kill any measurement quality confound
- For each z-bin: select quasars in narrow L, SNR ranges
  Split by sightline proxy
  Check if EW coupling differs while FWHM doesn't
"""

import numpy as np
from scipy import stats
from astropy.io import fits
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = Path('results_kill_B2')
RESULTS_DIR.mkdir(exist_ok=True)

# Load data
f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z = d['Z_DR16Q']
ra = d['RA']
dec = d['DEC']
ebv = d['EBV']
snr = d['SN_MEDIAN_ALL']
lbol = d['LOGLBOL']

lines = {}
for name in ['OIII5007', 'HBETA', 'CIV', 'CIII_ALL', 'MGII', 'LYA',
             'NII6585', 'SII6718', 'OII3728', 'HALPHA']:
    col = d[name]
    if col.ndim == 2 and col.shape[1] >= 6:
        lines[name] = {'ew': col[:, 2], 'ew_err': col[:, 3], 
                       'fwhm': col[:, 4], 'fwhm_err': col[:, 5]}
f.close()

# Galactic latitude
ra_r, dec_r = np.radians(ra), np.radians(dec)
ra_ngp, dec_ngp = np.radians(192.85948), np.radians(27.12825)
sin_b = np.sin(dec_ngp)*np.sin(dec_r) + np.cos(dec_ngp)*np.cos(dec_r)*np.cos(ra_r - ra_ngp)
b_gal = np.degrees(np.arcsin(np.clip(sin_b, -1, 1)))
abs_b = np.abs(b_gal)

print("=" * 80)
print("🔪 KILL B2 — ENVIRONMENT-CONTROLLED VOLUME-LIMITED TEST")
print("=" * 80)
print("""
B2 CLAIM: "High-z quasars live in denser environments. Denser environments
change both the source AND the sightline. The foreground effect is really
just environment-driven source evolution."

KILL STRATEGY: Fix EVERYTHING B2 could hide behind (luminosity, SNR,
redshift) and show the sightline effect persists.
""")

# =====================================================================
# TEST 1: VOLUME-LIMITED LUMINOSITY-MATCHED SIGHTLINE SPLIT
# =====================================================================
print("=" * 80)
print("TEST 1: VOLUME-LIMITED SIGHTLINE SPLIT")
print("Match on: z (narrow), luminosity (narrow), SNR (narrow)")
print("Split on: E(B-V) quartile (sightline proxy)")
print("Measure: EW coupling difference and FWHM coupling difference")
print("=" * 80)

# Define diagnostic pairs
ew_pairs = [
    ('CIV', 'CIII_ALL', 'CIV-CIII'),
    ('CIV', 'MGII', 'CIV-MgII'),
    ('OIII5007', 'HBETA', 'OIII-Hβ'),
    ('OIII5007', 'NII6585', 'OIII-NII (control)'),
]

z_bins = [
    (0.4, 0.6, 'sub-threshold'),
    (0.7, 0.9, 'threshold zone'),
    (0.9, 1.2, 'above threshold'),
    (1.2, 1.8, 'well above'),
    (1.8, 2.5, 'deep field'),
]

all_results = []

for z_lo, z_hi, z_label in z_bins:
    print(f"\n--- z = {z_lo}-{z_hi} ({z_label}) ---")
    
    z_mask = (z > z_lo) & (z < z_hi) & np.isfinite(z) & np.isfinite(lbol) & np.isfinite(snr)
    
    if z_mask.sum() < 200:
        print(f"  Too few objects ({z_mask.sum()}), skipping")
        continue
    
    # Volume limit: narrow luminosity range (middle 50%)
    lbol_z = lbol[z_mask]
    lbol_p25, lbol_p75 = np.nanpercentile(lbol_z[np.isfinite(lbol_z)], [25, 75])
    lbol_mask = (lbol > lbol_p25) & (lbol < lbol_p75)
    
    # SNR control: middle 50%
    snr_z = snr[z_mask]
    snr_p25, snr_p75 = np.nanpercentile(snr_z[np.isfinite(snr_z)], [25, 75])
    snr_mask = (snr > snr_p25) & (snr < snr_p75)
    
    base_mask = z_mask & lbol_mask & snr_mask
    
    if base_mask.sum() < 100:
        print(f"  Too few after L+SNR control ({base_mask.sum()}), skipping")
        continue
    
    # Split by E(B-V) — sightline proxy
    ebv_masked = ebv[base_mask]
    ebv_med = np.nanmedian(ebv_masked)
    ebv_q1 = np.nanpercentile(ebv_masked[np.isfinite(ebv_masked)], 25)
    ebv_q3 = np.nanpercentile(ebv_masked[np.isfinite(ebv_masked)], 75)
    
    clean_sight = base_mask & (ebv < ebv_q1)   # cleanest 25% sightlines
    dirty_sight = base_mask & (ebv > ebv_q3)    # dirtiest 25% sightlines
    
    print(f"  L range: {lbol_p25:.1f}-{lbol_p75:.1f}, SNR range: {snr_p25:.1f}-{snr_p75:.1f}")
    print(f"  E(B-V) split: clean < {ebv_q1:.4f}, dirty > {ebv_q3:.4f}")
    print(f"  N(clean) = {clean_sight.sum()}, N(dirty) = {dirty_sight.sum()}")
    
    for name_a, name_b, pair_label in ew_pairs:
        if name_a not in lines or name_b not in lines:
            continue
        
        ew_a, ew_b = lines[name_a]['ew'], lines[name_b]['ew']
        fwhm_a, fwhm_b = lines[name_a]['fwhm'], lines[name_b]['fwhm']
        
        # EW coupling
        valid = np.isfinite(ew_a) & np.isfinite(ew_b) & (ew_a != 0) & (ew_b != 0)
        
        m_clean = clean_sight & valid
        m_dirty = dirty_sight & valid
        
        if m_clean.sum() < 30 or m_dirty.sum() < 30:
            continue
        
        r_ew_clean, p_ew_clean = stats.pearsonr(ew_a[m_clean], ew_b[m_clean])
        r_ew_dirty, p_ew_dirty = stats.pearsonr(ew_a[m_dirty], ew_b[m_dirty])
        
        # FWHM coupling  
        valid_f = np.isfinite(fwhm_a) & np.isfinite(fwhm_b) & (fwhm_a > 0) & (fwhm_b > 0)
        m_clean_f = clean_sight & valid_f
        m_dirty_f = dirty_sight & valid_f
        
        r_fw_clean = r_fw_dirty = np.nan
        if m_clean_f.sum() >= 30 and m_dirty_f.sum() >= 30:
            r_fw_clean, _ = stats.pearsonr(fwhm_a[m_clean_f], fwhm_b[m_clean_f])
            r_fw_dirty, _ = stats.pearsonr(fwhm_a[m_dirty_f], fwhm_b[m_dirty_f])
        
        delta_ew = r_ew_dirty - r_ew_clean
        delta_fw = r_fw_dirty - r_fw_clean if np.isfinite(r_fw_clean) else np.nan
        
        result = {
            'z_bin': f'{z_lo}-{z_hi}', 'z_label': z_label,
            'pair': pair_label,
            'r_ew_clean': float(r_ew_clean), 'r_ew_dirty': float(r_ew_dirty),
            'delta_ew': float(delta_ew),
            'r_fw_clean': float(r_fw_clean) if np.isfinite(r_fw_clean) else None,
            'r_fw_dirty': float(r_fw_dirty) if np.isfinite(r_fw_dirty) else None,
            'delta_fw': float(delta_fw) if np.isfinite(delta_fw) else None,
            'n_clean': int(m_clean.sum()), 'n_dirty': int(m_dirty.sum()),
        }
        all_results.append(result)
        
        ew_flag = "⚠️ EW DIFFERS" if abs(delta_ew) > 0.02 else "  EW stable"
        fw_flag = "✅ FWHM stable" if (np.isnan(delta_fw) or abs(delta_fw) < 0.03) else "❌ FWHM differs"
        
        print(f"  {pair_label:<15} EW: clean={r_ew_clean:+.3f} dirty={r_ew_dirty:+.3f} Δ={delta_ew:+.4f}  {ew_flag}")
        if np.isfinite(r_fw_clean):
            print(f"  {'':<15} FW: clean={r_fw_clean:+.3f} dirty={r_fw_dirty:+.3f} Δ={delta_fw:+.4f}  {fw_flag}")


# =====================================================================
# TEST 2: GALACTIC LATITUDE AS INDEPENDENT SIGHTLINE PROXY
# =====================================================================
print(f"\n\n{'='*80}")
print("TEST 2: GALACTIC LATITUDE SPLIT (independent of E(B-V))")
print("Same volume-limited approach, different sightline proxy")
print("=" * 80)

for z_lo, z_hi, z_label in z_bins:
    z_mask = (z > z_lo) & (z < z_hi) & np.isfinite(z) & np.isfinite(lbol) & np.isfinite(snr)
    
    if z_mask.sum() < 200:
        continue
    
    lbol_z = lbol[z_mask]
    lbol_p25, lbol_p75 = np.nanpercentile(lbol_z[np.isfinite(lbol_z)], [25, 75])
    snr_z = snr[z_mask]
    snr_p25, snr_p75 = np.nanpercentile(snr_z[np.isfinite(snr_z)], [25, 75])
    
    base_mask = z_mask & (lbol > lbol_p25) & (lbol < lbol_p75) & (snr > snr_p25) & (snr < snr_p75)
    
    if base_mask.sum() < 100:
        continue
    
    # Split by |b| — high lat = clean, low lat = structured
    high_lat = base_mask & (abs_b > 60)
    low_lat = base_mask & (abs_b < 30)
    
    print(f"\n--- z = {z_lo}-{z_hi} ({z_label}) | N(high-lat)={high_lat.sum()}, N(low-lat)={low_lat.sum()} ---")
    
    for name_a, name_b, pair_label in ew_pairs:
        if name_a not in lines or name_b not in lines:
            continue
        
        ew_a, ew_b = lines[name_a]['ew'], lines[name_b]['ew']
        valid = np.isfinite(ew_a) & np.isfinite(ew_b) & (ew_a != 0) & (ew_b != 0)
        
        m_hi = high_lat & valid
        m_lo = low_lat & valid
        
        if m_hi.sum() < 20 or m_lo.sum() < 20:
            continue
        
        r_hi, _ = stats.pearsonr(ew_a[m_hi], ew_b[m_hi])
        r_lo, _ = stats.pearsonr(ew_a[m_lo], ew_b[m_lo])
        delta = r_lo - r_hi
        
        print(f"  {pair_label:<15} high-lat={r_hi:+.3f} low-lat={r_lo:+.3f} Δ={delta:+.4f}")


# =====================================================================
# TEST 3: THE DOUBLE CONTROL — Match on BOTH L and environment proxy
# =====================================================================
print(f"\n\n{'='*80}")
print("TEST 3: DOUBLE CONTROL — Match z + L + SNR + E(B-V) quartile")
print("If B2 is right, WITHIN the same E(B-V) quartile, there should be")
print("NO residual z-dependent degradation. If A1 is right, there should still")
print("be degradation because E(B-V) is only a weak proxy for cosmic web density.")
print("=" * 80)

# Within the dirty sightline quartile, does degradation still increase with z?
ebv_q3_global = np.nanpercentile(ebv[np.isfinite(ebv)], 75)
dirty_global = ebv > ebv_q3_global

print(f"\nWithin DIRTY sightlines only (E(B-V) > {ebv_q3_global:.4f}):")
print(f"Does EW coupling still degrade with z?")
print(f"\n  {'Line pair':<15} {'z-bin':<12} {'r(EW)':<10} {'r(FWHM)':<10} {'N'}")
print(f"  {'-'*55}")

z_track = []
for name_a, name_b, pair_label in [('CIV', 'CIII_ALL', 'CIV-CIII'), ('OIII5007', 'HBETA', 'OIII-Hβ')]:
    if name_a not in lines or name_b not in lines:
        continue
    
    ew_a, ew_b = lines[name_a]['ew'], lines[name_b]['ew']
    fwhm_a, fwhm_b = lines[name_a]['fwhm'], lines[name_b]['fwhm']
    
    for z_lo, z_hi, z_label in z_bins:
        z_mask = (z > z_lo) & (z < z_hi) & dirty_global & np.isfinite(lbol) & np.isfinite(snr)
        valid_ew = z_mask & np.isfinite(ew_a) & np.isfinite(ew_b) & (ew_a != 0) & (ew_b != 0)
        valid_fw = z_mask & np.isfinite(fwhm_a) & np.isfinite(fwhm_b) & (fwhm_a > 0) & (fwhm_b > 0)
        
        if valid_ew.sum() < 30:
            continue
        
        r_ew, _ = stats.pearsonr(ew_a[valid_ew], ew_b[valid_ew])
        r_fw = np.nan
        if valid_fw.sum() >= 30:
            r_fw, _ = stats.pearsonr(fwhm_a[valid_fw], fwhm_b[valid_fw])
        
        z_mid = (z_lo + z_hi) / 2
        z_track.append({'pair': pair_label, 'z_mid': z_mid, 'r_ew': r_ew, 'r_fw': r_fw})
        
        fw_str = f"{r_fw:+.3f}" if np.isfinite(r_fw) else "  N/A"
        print(f"  {pair_label:<15} {z_lo:.1f}-{z_hi:.1f}    {r_ew:>+.4f}    {fw_str}   {valid_ew.sum()}")

# Check if r_ew still declines with z within dirty sightlines
for pair_label in ['CIV-CIII', 'OIII-Hβ']:
    pts = [t for t in z_track if t['pair'] == pair_label and np.isfinite(t['r_ew'])]
    if len(pts) >= 3:
        zz = [p['z_mid'] for p in pts]
        rr = [p['r_ew'] for p in pts]
        slope, _, r, p, _ = stats.linregress(zz, rr)
        print(f"\n  {pair_label}: r(EW) vs z slope = {slope:+.4f}, trend r={r:+.3f}, p={p:.4f}")
        if p < 0.05:
            print(f"  → EW STILL DEGRADES with z even within fixed dirty sightlines!")
            print(f"  → B2 cannot explain this: environment is HELD FIXED")
        else:
            print(f"  → No significant z-trend within fixed environment (B2 survives here)")


# =====================================================================
# TEST 4: PROPENSITY SCORE MATCHING
# =====================================================================
print(f"\n\n{'='*80}")
print("TEST 4: PROPENSITY-MATCHED PAIRS")
print("For each dirty-sightline quasar, find a clean-sightline quasar")
print("matched on z (±0.02), L (±0.1 dex), SNR (±20%)")
print("Then compare their EW and FWHM distributions.")
print("This is the gold standard — eliminates ALL confounders B2 can invoke.")
print("=" * 80)

# Focus on CIV (best statistics above threshold)
if 'CIV' in lines and 'CIII_ALL' in lines:
    ew_civ = lines['CIV']['ew']
    ew_ciii = lines['CIII_ALL']['ew']
    fwhm_civ = lines['CIV']['fwhm']
    fwhm_ciii = lines['CIII_ALL']['fwhm']
    
    # Select above-threshold quasars with good measurements
    above_thresh = (z > 1.0) & (z < 2.5) & np.isfinite(z) & np.isfinite(lbol) & np.isfinite(snr)
    above_thresh &= np.isfinite(ew_civ) & np.isfinite(ew_ciii) & (ew_civ != 0) & (ew_ciii != 0)
    above_thresh &= np.isfinite(fwhm_civ) & np.isfinite(fwhm_ciii) & (fwhm_civ > 0) & (fwhm_ciii > 0)
    
    # Split into dirty and clean
    ebv_q25 = np.nanpercentile(ebv[above_thresh], 25)
    ebv_q75 = np.nanpercentile(ebv[above_thresh], 75)
    
    idx_all = np.where(above_thresh)[0]
    idx_clean = idx_all[ebv[idx_all] < ebv_q25]
    idx_dirty = idx_all[ebv[idx_all] > ebv_q75]
    
    print(f"  Above threshold: {above_thresh.sum()} quasars")
    print(f"  Clean sightlines: {len(idx_clean)}, Dirty sightlines: {len(idx_dirty)}")
    
    # Match each dirty to nearest clean on (z, L, SNR)
    matched_clean = []
    matched_dirty = []
    
    # For speed, do block matching
    np.random.seed(42)
    n_match = min(5000, len(idx_dirty))
    sample_dirty = np.random.choice(idx_dirty, n_match, replace=False)
    
    for i_d in sample_dirty:
        dz = np.abs(z[idx_clean] - z[i_d])
        dl = np.abs(lbol[idx_clean] - lbol[i_d])
        ds = np.abs(snr[idx_clean] - snr[i_d]) / (snr[i_d] + 1)
        
        # Strict matching
        ok = (dz < 0.02) & (dl < 0.1) & (ds < 0.2)
        if ok.sum() == 0:
            ok = (dz < 0.05) & (dl < 0.2) & (ds < 0.3)
        if ok.sum() == 0:
            continue
        
        # Pick closest
        dist = dz + dl / 10 + ds
        dist[~ok] = 999
        best = np.argmin(dist)
        
        matched_dirty.append(i_d)
        matched_clean.append(idx_clean[best])
    
    matched_dirty = np.array(matched_dirty)
    matched_clean = np.array(matched_clean)
    
    print(f"  Matched pairs: {len(matched_dirty)}")
    
    if len(matched_dirty) >= 100:
        # Verify match quality
        dz_match = np.abs(z[matched_dirty] - z[matched_clean])
        dl_match = np.abs(lbol[matched_dirty] - lbol[matched_clean])
        print(f"  Match quality: median Δz={np.median(dz_match):.4f}, median ΔL={np.median(dl_match):.3f}")
        
        # Compare EW distributions
        ew_ratio_dirty = ew_civ[matched_dirty] / (ew_ciii[matched_dirty] + 1e-10)
        ew_ratio_clean = ew_civ[matched_clean] / (ew_ciii[matched_clean] + 1e-10)
        
        # Clip extreme ratios
        valid_ratio = (np.abs(ew_ratio_dirty) < 100) & (np.abs(ew_ratio_clean) < 100)
        ew_ratio_dirty = ew_ratio_dirty[valid_ratio]
        ew_ratio_clean = ew_ratio_clean[valid_ratio]
        
        ks_ew, p_ew = stats.ks_2samp(ew_ratio_dirty, ew_ratio_clean)
        
        # Compare FWHM distributions
        fwhm_ratio_dirty = fwhm_civ[matched_dirty] / (fwhm_ciii[matched_dirty] + 1e-10)
        fwhm_ratio_clean = fwhm_civ[matched_clean] / (fwhm_ciii[matched_clean] + 1e-10)
        
        valid_fw = (fwhm_ratio_dirty > 0) & (fwhm_ratio_dirty < 100) & (fwhm_ratio_clean > 0) & (fwhm_ratio_clean < 100)
        fwhm_ratio_dirty_v = fwhm_ratio_dirty[valid_fw]
        fwhm_ratio_clean_v = fwhm_ratio_clean[valid_fw]
        
        ks_fw, p_fw = stats.ks_2samp(fwhm_ratio_dirty_v, fwhm_ratio_clean_v)
        
        print(f"\n  PROPENSITY-MATCHED RESULTS (CIV/CIII ratio):")
        print(f"  EW ratio:   KS={ks_ew:.4f}, p={p_ew:.2e}  {'⚠️ DIFFERS' if p_ew < 0.01 else '✅ Same'}")
        print(f"  FWHM ratio: KS={ks_fw:.4f}, p={p_fw:.2e}  {'❌ DIFFERS' if p_fw < 0.01 else '✅ Same'}")
        
        if p_ew < 0.01 and p_fw > 0.01:
            print(f"\n  🔪 B2 IS DEAD: After matching on z, L, and SNR,")
            print(f"     dirty sightlines STILL show different EW ratios")
            print(f"     but SAME FWHM ratios. Environment cannot explain this.")
            print(f"     The sightline itself is doing something SELECTIVE.")
        elif p_ew < 0.01 and p_fw < 0.01:
            print(f"\n  ⚠️ BOTH differ — could be residual environment effect")
        elif p_ew > 0.01:
            print(f"\n  B2 survives: no significant EW difference after matching")
        
        # Also check raw EW and FWHM individually
        print(f"\n  Individual line checks (matched pairs):")
        for lname, arr_type in [('CIV', 'ew'), ('CIV', 'fwhm'), ('CIII_ALL', 'ew'), ('CIII_ALL', 'fwhm')]:
            vals_d = lines[lname][arr_type][matched_dirty]
            vals_c = lines[lname][arr_type][matched_clean]
            valid = np.isfinite(vals_d) & np.isfinite(vals_c) & (vals_d != 0) & (vals_c != 0)
            if valid.sum() > 50:
                ks, p = stats.ks_2samp(vals_d[valid], vals_c[valid])
                med_d = np.median(vals_d[valid])
                med_c = np.median(vals_c[valid])
                tag = "⚠️" if p < 0.01 else "✅"
                print(f"    {lname} {arr_type}: dirty median={med_d:.2f}, clean median={med_c:.2f}, KS p={p:.2e} {tag}")


# =====================================================================
# SUMMARY
# =====================================================================
print(f"\n\n{'='*80}")
print("SUMMARY — B2 KILL TEST")
print("=" * 80)
print("""
B2 (selection + environment) requires that ALL of these are true:
1. Environment drives intrinsic line properties differently for EW vs FWHM
2. Environment correlates with sightline density
3. This conspiracy reproduces the doublet ladder ordering
4. This conspiracy works across SNe, quasars, AND FRBs
5. This conspiracy survives propensity matching on z, L, SNR

If the propensity-matched test shows EW differs but FWHM doesn't,
B2 is dead. The sightline itself is selectively acting on information
content, independent of source properties.
""")

# Save results
with open(RESULTS_DIR / 'results.json', 'w') as f:
    json.dump(all_results, f, indent=2)
print(f"Results saved to {RESULTS_DIR}/")
