#!/usr/bin/env python3
"""
SHRINK THE POSSIBILITIES
=========================
Tests that distinguish:
  A1: Medium physically reweights subzones (depends on foreground properties)
  D:  Frame-dependent observables degrade, geometric ones don't (depends on observable TYPE)

KEY PREDICTIONS:
  A1 says: effect depends on MEDIUM (density, structure, wavelength)
  D says:  effect depends on OBSERVABLE TYPE (relational vs geometric)
  
  Where they AGREE: EW degrades, FWHM doesn't
  Where they DIFFER:
    - A1: wavelength should matter independently
    - D: ONLY observable type matters, wavelength is noise
    - A1: foreground density is the driver
    - D: path complexity (distance) is the driver
    - D: ALL geometric quantities survive, ALL relational ones degrade
    - A1: chromatic selectivity means wavelength-dependent, not just type-dependent

TESTS:
1. Peak wavelength stability (geometric) — does it survive like FWHM?
2. Line asymmetry stability — partially geometric
3. Velocity offsets between lines — geometric (difference of kinematics)
4. Is the driver DISTANCE or MEDIUM? (already partially done, refine here)
5. Do ALL relational observables degrade equally per unit P?
"""

import numpy as np
from scipy import stats
from astropy.io import fits
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = Path('results_shrink')
RESULTS_DIR.mkdir(exist_ok=True)

f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
z = d['Z_DR16Q']
ebv = d['EBV']
snr = d['SN_MEDIAN_ALL']

lines = {}
for name in ['OIII5007', 'HBETA', 'CIV', 'CIII_ALL', 'MGII', 'LYA',
             'NII6585', 'SII6718', 'OII3728', 'HALPHA', 'HBETA_BR',
             'HEII4687', 'HEII1640', 'NV1240']:
    col = d[name]
    if col.ndim == 2 and col.shape[1] >= 6:
        lines[name] = {
            'peak': col[:, 0],       # peak wavelength
            'peak2': col[:, 1],      # secondary peak
            'ew': col[:, 2],         # equivalent width
            'ew_err': col[:, 3],     # EW error
            'fwhm': col[:, 4],       # FWHM
            'fwhm_err': col[:, 5],   # FWHM error
        }
f.close()

print("=" * 80)
print("🔬 SHRINK THE POSSIBILITIES — A1 vs D")
print("=" * 80)

# =====================================================================
# TEST 1: PEAK WAVELENGTH STABILITY
# =====================================================================
print(f"\n{'='*80}")
print("TEST 1: PEAK WAVELENGTH CORRELATIONS (geometric observable)")
print("D predicts: peak λ correlations should SURVIVE like FWHM (geometric)")
print("A1 predicts: peak λ could shift if medium refracts chromatically")
print("=" * 80)

# For each z-bin, correlate peak wavelength between line pairs
# Peak wavelength encodes the centroid velocity — geometric/kinematic
z_bins = [(0.3, 0.7), (0.7, 1.0), (1.0, 1.5), (1.5, 2.5)]

test_pairs = [
    ('OIII5007', 'HBETA', 'OIII-Hβ'),
    ('CIV', 'CIII_ALL', 'CIV-CIII'),
    ('CIV', 'MGII', 'CIV-MgII'),
]

print(f"\n  Peak wavelength inter-line correlation:")
print(f"  {'Pair':<15} {'z-bin':<10} {'r(peak)':<12} {'r(EW)':<12} {'r(FWHM)':<12} {'Peak=FWHM?'}")
print(f"  {'-'*70}")

peak_results = []
for name_a, name_b, label in test_pairs:
    if name_a not in lines or name_b not in lines:
        continue
    
    for z_lo, z_hi in z_bins:
        z_mask = (z > z_lo) & (z < z_hi)
        
        # Peak correlation
        pk_a, pk_b = lines[name_a]['peak'], lines[name_b]['peak']
        valid_pk = z_mask & np.isfinite(pk_a) & np.isfinite(pk_b) & (pk_a > 0) & (pk_b > 0)
        
        # EW correlation
        ew_a, ew_b = lines[name_a]['ew'], lines[name_b]['ew']
        valid_ew = z_mask & np.isfinite(ew_a) & np.isfinite(ew_b) & (ew_a != 0) & (ew_b != 0)
        
        # FWHM correlation
        fw_a, fw_b = lines[name_a]['fwhm'], lines[name_b]['fwhm']
        valid_fw = z_mask & np.isfinite(fw_a) & np.isfinite(fw_b) & (fw_a > 0) & (fw_b > 0)
        
        r_pk = r_ew = r_fw = np.nan
        if valid_pk.sum() >= 50:
            r_pk, _ = stats.spearmanr(pk_a[valid_pk], pk_b[valid_pk])
        if valid_ew.sum() >= 50:
            r_ew, _ = stats.spearmanr(ew_a[valid_ew], ew_b[valid_ew])
        if valid_fw.sum() >= 50:
            r_fw, _ = stats.spearmanr(fw_a[valid_fw], fw_b[valid_fw])
        
        if np.isfinite(r_pk):
            # Does peak behave like FWHM (geometric) or like EW (relational)?
            d_to_fw = abs(r_pk - r_fw) if np.isfinite(r_fw) else np.nan
            d_to_ew = abs(r_pk - r_ew) if np.isfinite(r_ew) else np.nan
            
            if np.isfinite(d_to_fw) and np.isfinite(d_to_ew):
                tag = "≈FWHM (geometric)" if d_to_fw < d_to_ew else "≈EW (relational)"
            else:
                tag = ""
            
            pk_str = f"{r_pk:+.4f}" 
            ew_str = f"{r_ew:+.4f}" if np.isfinite(r_ew) else "N/A"
            fw_str = f"{r_fw:+.4f}" if np.isfinite(r_fw) else "N/A"
            
            print(f"  {label:<15} {z_lo:.1f}-{z_hi:.1f}  {pk_str:<12} {ew_str:<12} {fw_str:<12} {tag}")
            
            peak_results.append({
                'pair': label, 'z': f'{z_lo}-{z_hi}',
                'r_peak': float(r_pk), 'r_ew': float(r_ew) if np.isfinite(r_ew) else None,
                'r_fwhm': float(r_fw) if np.isfinite(r_fw) else None
            })


# =====================================================================
# TEST 2: VELOCITY OFFSET STABILITY
# =====================================================================
print(f"\n{'='*80}")
print("TEST 2: VELOCITY OFFSET BETWEEN LINES (geometric difference)")
print("Velocity offset = (peak_A - peak_B) / peak_B × c")
print("D predicts: offsets are geometric, should be STABLE with z")
print("A1 predicts: if medium shifts centroids chromatically, offsets could drift")
print("=" * 80)

for name_a, name_b, label in test_pairs:
    if name_a not in lines or name_b not in lines:
        continue
    
    pk_a, pk_b = lines[name_a]['peak'], lines[name_b]['peak']
    ew_a = lines[name_a]['ew']
    
    valid = np.isfinite(pk_a) & np.isfinite(pk_b) & (pk_a > 0) & (pk_b > 0)
    
    # Compute velocity offset
    v_offset = (pk_a - pk_b) / pk_b * 3e5  # km/s
    
    print(f"\n  {label} velocity offset vs z:")
    print(f"  {'z-bin':<10} {'median Δv':<15} {'std Δv':<15} {'N'}")
    
    offsets_by_z = []
    for z_lo, z_hi in [(0.3, 0.6), (0.6, 0.9), (0.9, 1.2), (1.2, 1.6), (1.6, 2.0), (2.0, 2.5), (2.5, 3.5)]:
        m = valid & (z > z_lo) & (z < z_hi)
        if m.sum() < 30:
            continue
        
        v = v_offset[m]
        v = v[np.isfinite(v) & (np.abs(v) < 1e5)]  # clip extreme
        if len(v) < 30:
            continue
        
        med = np.median(v)
        std = np.std(v)
        offsets_by_z.append({'z_mid': (z_lo+z_hi)/2, 'median': med, 'std': std, 'n': len(v)})
        print(f"  {z_lo:.1f}-{z_hi:.1f}    {med:>+8.1f} km/s   {std:>8.1f} km/s   {len(v)}")
    
    if len(offsets_by_z) >= 3:
        zz = [o['z_mid'] for o in offsets_by_z]
        stds = [o['std'] for o in offsets_by_z]
        slope, _, r, p, _ = stats.linregress(zz, stds)
        print(f"  Scatter vs z: slope={slope:+.1f}, r={r:+.3f}, p={p:.4f}")
        if p < 0.05 and slope > 0:
            print(f"  → Velocity offset scatter GROWS with z → A1 favored (medium shifts centroids)")
        elif p > 0.05:
            print(f"  → Velocity offset scatter STABLE → D favored (geometric invariant)")


# =====================================================================
# TEST 3: CLASSIFY ALL OBSERVABLES BY TYPE
# =====================================================================
print(f"\n{'='*80}")
print("TEST 3: OBSERVABLE TAXONOMY — Which ones survive?")
print("Classify each column type, measure stability with z")
print("D predicts: geometric survives, relational degrades. Clean split.")
print("A1 predicts: wavelength-dependent, not just type-dependent.")
print("=" * 80)

# For CIV (best statistics above threshold):
if 'CIV' in lines and 'CIII_ALL' in lines:
    observables = {
        'peak (geometric)':    (lines['CIV']['peak'], lines['CIII_ALL']['peak'], 
                                lambda x: x > 0, lambda x: x > 0),
        'FWHM (geometric)':    (lines['CIV']['fwhm'], lines['CIII_ALL']['fwhm'],
                                lambda x: x > 0, lambda x: x > 0),
        'EW (relational)':     (lines['CIV']['ew'], lines['CIII_ALL']['ew'],
                                lambda x: x != 0, lambda x: x != 0),
    }
    
    print(f"\n  CIV-CIII correlation by observable type across z:")
    print(f"  {'Observable':<25} {'z=1.0-1.3':<12} {'z=1.3-1.7':<12} {'z=1.7-2.2':<12} {'z=2.2-3.0':<12} {'Trend'}")
    print(f"  {'-'*80}")
    
    z_test = [(1.0, 1.3), (1.3, 1.7), (1.7, 2.2), (2.2, 3.0)]
    
    for obs_name, (arr_a, arr_b, filt_a, filt_b) in observables.items():
        row = []
        for z_lo, z_hi in z_test:
            m = (z > z_lo) & (z < z_hi) & np.isfinite(arr_a) & np.isfinite(arr_b) & filt_a(arr_a) & filt_b(arr_b)
            if m.sum() >= 50:
                r, _ = stats.spearmanr(arr_a[m], arr_b[m])
                row.append(r)
            else:
                row.append(np.nan)
        
        vals = [f"{r:+.4f}" if np.isfinite(r) else "N/A" for r in row]
        
        # Trend
        valid_r = [(i, r) for i, r in enumerate(row) if np.isfinite(r)]
        if len(valid_r) >= 3:
            idx, rr = zip(*valid_r)
            sl, _, tr, tp, _ = stats.linregress(idx, rr)
            trend = f"{'↓' if sl < -0.01 else '↑' if sl > 0.01 else '→'} (p={tp:.3f})"
        else:
            trend = ""
        
        print(f"  {obs_name:<25} {vals[0]:<12} {vals[1]:<12} {vals[2]:<12} {vals[3]:<12} {trend}")


# =====================================================================
# TEST 4: DISTANCE vs MEDIUM — Which matters MORE?
# =====================================================================
print(f"\n{'='*80}")
print("TEST 4: PARTIAL CORRELATION — z vs E(B-V) as predictors")
print("A1: foreground density (E(B-V) proxy) should add predictive power beyond z")
print("D:  distance alone should be sufficient; foreground adds nothing")
print("=" * 80)

# For MgII (pipeline-confirmed line)
if 'MGII' in lines:
    ew_mg = lines['MGII']['ew']
    fw_mg = lines['MGII']['fwhm']
    
    # Compute EW/FWHM ratio as degradation proxy
    valid = np.isfinite(ew_mg) & np.isfinite(fw_mg) & (ew_mg != 0) & (fw_mg > 0) & np.isfinite(z) & np.isfinite(ebv)
    valid &= (z > 0.3) & (z < 2.5)
    
    ew_v = np.abs(ew_mg[valid])
    fw_v = fw_mg[valid]
    z_v = z[valid]
    ebv_v = ebv[valid]
    
    # Log-transform for better scaling
    ratio = np.log10(ew_v / fw_v + 1e-10)
    
    # Partial correlation: ratio ~ z, controlling for E(B-V)
    # Partial correlation: ratio ~ E(B-V), controlling for z
    
    from numpy.linalg import lstsq
    
    def partial_corr(x, y, control):
        """Partial Spearman correlation of x,y controlling for control."""
        # Rank everything
        rx = stats.rankdata(x)
        ry = stats.rankdata(y)
        rc = stats.rankdata(control)
        
        # Residualize x and y on control
        A = np.column_stack([rc, np.ones(len(rc))])
        res_x = rx - A @ lstsq(A, rx, rcond=None)[0]
        res_y = ry - A @ lstsq(A, ry, rcond=None)[0]
        
        return stats.pearsonr(res_x, res_y)
    
    # Subsample for speed
    n_sub = min(50000, len(z_v))
    idx = np.random.RandomState(42).choice(len(z_v), n_sub, replace=False)
    
    r_z_raw, p_z_raw = stats.spearmanr(z_v[idx], ratio[idx])
    r_ebv_raw, p_ebv_raw = stats.spearmanr(ebv_v[idx], ratio[idx])
    
    r_z_partial, p_z_partial = partial_corr(z_v[idx], ratio[idx], ebv_v[idx])
    r_ebv_partial, p_ebv_partial = partial_corr(ebv_v[idx], ratio[idx], z_v[idx])
    
    print(f"\n  MgII EW/FWHM ratio predictions (N={n_sub}):")
    print(f"  Raw correlation with z:      r={r_z_raw:+.4f}, p={p_z_raw:.2e}")
    print(f"  Raw correlation with E(B-V): r={r_ebv_raw:+.4f}, p={p_ebv_raw:.2e}")
    print(f"  Partial r(z | E(B-V)):       r={r_z_partial:+.4f}, p={p_z_partial:.2e}")
    print(f"  Partial r(E(B-V) | z):       r={r_ebv_partial:+.4f}, p={p_ebv_partial:.2e}")
    
    if abs(r_ebv_partial) > 0.01 and p_ebv_partial < 0.01:
        print(f"\n  → E(B-V) adds predictive power BEYOND z → A1 favored (medium matters)")
    elif abs(r_ebv_partial) < 0.01 or p_ebv_partial > 0.05:
        print(f"\n  → E(B-V) adds NO power beyond z → D favored (distance is sufficient)")
    
    if abs(r_z_partial) > abs(r_ebv_partial):
        print(f"  → z is the STRONGER predictor → leans D (distance/perspective)")
    else:
        print(f"  → E(B-V) is STRONGER than z → leans A1 (medium properties)")


# =====================================================================
# TEST 5: WAVELENGTH vs TYPE — Which organizes the data?
# =====================================================================
print(f"\n{'='*80}")
print("TEST 5: WAVELENGTH vs OBSERVABLE TYPE as organizer")
print("A1: degradation should correlate with rest wavelength (chromatic)")
print("D:  degradation should correlate with P (susceptibility type) only")
print("=" * 80)

# Use our line rankings
line_info = {
    # name: (rest_wav_Å, empirical_P, type)
    'LYA':      (1216, 1.00, 'relational'),
    'CIV':      (1549, 0.75, 'relational'),
    'CIII_ALL': (1909, 0.60, 'relational'),
    'MGII':     (2798, 0.65, 'relational'),
    'HBETA':    (4861, 0.30, 'relational'),
    'OIII5007': (5007, 0.05, 'geometric-like'),   # branch-locked
    'NII6585':  (6585, 0.05, 'geometric-like'),   # branch-locked
    'HALPHA':   (6563, 0.25, 'relational'),
    'SII6718':  (6718, 0.70, 'relational'),
    'OII3728':  (3728, 0.40, 'relational'),
}

# Measure decorrelation rate for each line where possible
print(f"\n  Per-line EW coefficient of variation growth with z:")
print(f"  {'Line':<12} {'λ_rest':<8} {'P':<6} {'CV slope':<12} {'r':<8} {'p':<10}")
print(f"  {'-'*55}")

cv_slopes = []
for lname, (lam, P, ltype) in line_info.items():
    if lname not in lines:
        continue
    
    ew = lines[lname]['ew']
    fw = lines[lname]['fwhm']
    
    z_mids = []
    cv_ews = []
    
    for z_lo, z_hi in [(0.1, 0.4), (0.4, 0.7), (0.7, 1.0), (1.0, 1.3), (1.3, 1.8), (1.8, 2.5), (2.5, 3.5)]:
        m = (z > z_lo) & (z < z_hi) & np.isfinite(ew) & (ew != 0)
        if m.sum() < 100:
            continue
        
        cv = np.std(ew[m]) / (np.abs(np.mean(ew[m])) + 1e-10)
        z_mids.append((z_lo + z_hi) / 2)
        cv_ews.append(cv)
    
    if len(z_mids) >= 3:
        sl, _, r, p, _ = stats.linregress(z_mids, cv_ews)
        cv_slopes.append({'name': lname, 'lam': lam, 'P': P, 'slope': sl, 'r': r, 'p': p, 'type': ltype})
        print(f"  {lname:<12} {lam:<8} {P:<6.2f} {sl:>+10.3f}   {r:>+.3f}  {p:.2e}")

if len(cv_slopes) >= 5:
    lams = [c['lam'] for c in cv_slopes]
    Ps = [c['P'] for c in cv_slopes]
    slopes = [c['slope'] for c in cv_slopes]
    
    r_lam, p_lam = stats.spearmanr(lams, slopes)
    r_P, p_P = stats.spearmanr(Ps, slopes)
    
    print(f"\n  Correlation of CV slope with:")
    print(f"    Wavelength: r={r_lam:+.3f}, p={p_lam:.4f}")
    print(f"    P (susceptibility): r={r_P:+.3f}, p={p_P:.4f}")
    
    if abs(r_P) > abs(r_lam):
        print(f"\n  → P organizes better than wavelength → D favored")
        print(f"     (Observable type > chromatic dependence)")
    else:
        print(f"\n  → Wavelength organizes better → A1 favored")
        print(f"     (Chromatic process > observable type)")


# =====================================================================
# SUMMARY
# =====================================================================
print(f"\n{'='*80}")
print("SUMMARY — A1 vs D scorecard")
print("=" * 80)
print("""
A1 (medium / refractive decoherence):
  + Foreground matters at same z (z-matched test p=10⁻⁵)
  + Blue-preferential (wavelength dependence exists)  
  + Simulation matches threshold
  
D (perspective / frame-dependence):
  + P organizes better than wavelength (R²=0.23 for λ vs r=-0.726 for P)
  + Observable TYPE predicts survival (geometric vs relational)
  + The equation is dr/dz = -0.10 × P, not dr/dz = f(λ)

BOTH:
  + EW degrades, FWHM doesn't
  + Threshold exists
  + Cross-class universality

The question: Is the wavelength dependence a CONSEQUENCE of type-dependence
(relational observables happen to be bluer) or an independent variable?
""")
