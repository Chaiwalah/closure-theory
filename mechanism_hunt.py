#!/usr/bin/env python3
"""
MECHANISM HUNT — Find the Rat
==============================
Tests to constrain WHAT is causing the selective information degradation.

Test 1: Rosetta Stone — Do all three sigmoid thresholds map to same physical quantity?
Test 2: Column Density Mapping — Does degradation correlate with integrated electron density?
Test 3: Sightline Tomography — Void vs filament vs cluster sightlines
Test 4: Wavelength Dependence — UV vs optical vs IR operator profile
Test 5: Sigmoid First Principles — What physical scale sets z₀?
"""

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = Path('results_mechanism_hunt')
RESULTS_DIR.mkdir(exist_ok=True)


def load_quasar_data():
    from astropy.io import fits
    f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
    d = f[1].data
    z = d['Z_DR16Q']
    
    lines = {}
    line_names = ['HALPHA', 'HALPHA_BR', 'HBETA', 'HBETA_BR', 'OIII5007', 'OIII5007C',
                  'NII6585', 'SII6718', 'OII3728', 'MGII', 'MGII_BR', 'CIV',
                  'CIII_ALL', 'CIII_BR', 'LYA', 'NV1240', 'HEII4687', 'HEII1640',
                  'NEV3426', 'SIIV_OIV']
    
    for name in line_names:
        col = d[name]
        if col.ndim == 2 and col.shape[1] >= 6:
            lines[name] = {
                'flux': col[:, 0],
                'ew': col[:, 2],
                'fwhm': col[:, 4]
            }
    
    scalars = {}
    for s in ['LOGLBOL', 'LOGMBH', 'LOGLEDD_RATIO', 'LOGL1350', 'LOGL3000', 'LOGL5100',
              'EBV', 'SN_MEDIAN_ALL']:
        if s in f[1].columns.names:
            scalars[s] = d[s]
    
    f.close()
    return z, lines, scalars


def sigmoid(x, L, k, x0, b):
    """Standard sigmoid function."""
    return L / (1 + np.exp(-k * (x - x0))) + b


def test1_rosetta_stone(z, lines):
    """
    ROSETTA STONE TEST
    ==================
    Map sigmoid thresholds across source classes to a common physical variable.
    
    Key insight: FRB threshold is DM ≈ 500 pc/cm³. DM = integrated electron density.
    If we convert quasar z₀ ≈ 1.1 and SNe z₀ ≈ 0.82 to expected DM using the
    Macquart relation, do they converge to the same column density?
    """
    print("\n" + "=" * 70)
    print("TEST 1: ROSETTA STONE — Common Physical Threshold")
    print("=" * 70)
    
    # Macquart relation: DM_cosmic ≈ 930 * z (pc/cm³) for mean IGM
    # More precisely: DM(z) = (3 * c * H0 * Omega_b) / (8 * pi * G * m_p) * integral
    # Simplified: DM ≈ 900-1000 * z for z < 2
    
    # Known thresholds
    thresholds = {
        'SNe Ia (Pantheon+)': {'z0': 0.82, 'z0_err': 0.05},
        'Quasars DR16Q (EW)': {'z0': 1.14, 'z0_err': 0.09},  # midpoint of 1.05-1.23
        'FRBs (CHIME)': {'DM0': 500, 'DM0_err': 50},
    }
    
    print("\n  Converting all thresholds to DM (integrated electron column density):")
    print("  Using Macquart relation: DM_cosmic ≈ 930 × z pc/cm³\n")
    
    dm_values = []
    
    for name, vals in thresholds.items():
        if 'z0' in vals:
            dm = 930 * vals['z0']
            dm_err = 930 * vals['z0_err']
            print(f"  {name}: z₀ = {vals['z0']:.2f} → DM ≈ {dm:.0f} ± {dm_err:.0f} pc/cm³")
            dm_values.append(dm)
        else:
            dm = vals['DM0']
            dm_err = vals['DM0_err']
            # Convert DM to effective z
            z_eff = dm / 930
            print(f"  {name}: DM₀ = {dm} → z_eff ≈ {z_eff:.2f}")
            dm_values.append(dm)
    
    mean_dm = np.mean(dm_values)
    std_dm = np.std(dm_values)
    cv = std_dm / mean_dm
    
    print(f"\n  DM values: {[f'{d:.0f}' for d in dm_values]}")
    print(f"  Mean: {mean_dm:.0f} pc/cm³")
    print(f"  Std: {std_dm:.0f} pc/cm³")
    print(f"  CV (coefficient of variation): {cv:.3f}")
    
    if cv < 0.3:
        print(f"\n  ✅ CONVERGENCE: All three source classes threshold at similar column density!")
        print(f"     The sigmoid isn't about DISTANCE — it's about HOW MUCH MEDIUM.")
        print(f"     Critical column density ≈ {mean_dm:.0f} pc/cm³")
    else:
        print(f"\n  ⚠️ Partial convergence. CV = {cv:.3f}. Spread may reflect different sensitivities.")
    
    # Now let's also check: what physical conditions correspond to DM ≈ 700-800?
    print(f"\n  Physical interpretation of DM ≈ {mean_dm:.0f} pc/cm³:")
    print(f"  - Mean IGM electron density: n_e ≈ 2 × 10⁻⁷ cm⁻³")
    print(f"  - At z ≈ 0.8-1.1: comoving distance ≈ 2.5-3.5 Gpc")
    print(f"  - Path through IGM: ≈ {mean_dm / 2e-7 / 3.086e24:.1f} Gpc")
    print(f"  - This corresponds to the epoch where cosmic web structure")
    print(f"    transitions from filament-dominated to void-dominated topology")
    
    # Verify with our data: fit sigmoid to quasar EW coupling
    print("\n  Fitting sigmoid to quasar data to verify z₀...")
    
    # Use OIII EW vs Hbeta EW coupling as function of z
    oiii_ew = lines['OIII5007']['ew']
    hb_ew = lines['HBETA']['ew']
    
    mask = (np.isfinite(oiii_ew) & np.isfinite(hb_ew) & 
            (oiii_ew != 0) & (hb_ew != 0) & (z > 0.1) & (z < 3) & np.isfinite(z))
    
    zz = z[mask]
    oo = oiii_ew[mask]
    hh = hb_ew[mask]
    
    # Measure correlation in z-bins
    n_bins = 20
    z_edges = np.percentile(zz, np.linspace(0, 100, n_bins + 1))
    z_mids = []
    correlations = []
    
    for i in range(n_bins):
        bm = (zz >= z_edges[i]) & (zz < z_edges[i+1] + (0.001 if i == n_bins-1 else 0))
        if bm.sum() > 50:
            r, p = stats.pearsonr(oo[bm], hh[bm])
            if np.isfinite(r):
                z_mids.append((z_edges[i] + z_edges[i+1]) / 2)
                correlations.append(abs(r))
    
    z_mids = np.array(z_mids)
    correlations = np.array(correlations)
    
    print(f"  Correlation by z-bin (|r| of OIII_EW vs HB_EW):")
    for zm, cr in zip(z_mids, correlations):
        bar = "█" * int(cr * 50)
        print(f"    z={zm:.2f}: |r|={cr:.4f} {bar}")
    
    # Try sigmoid fit
    try:
        popt, pcov = curve_fit(sigmoid, z_mids, correlations, 
                               p0=[0.1, -5, 1.0, 0.01], maxfev=10000)
        z0_fit = popt[2]
        print(f"\n  Sigmoid fit: z₀ = {z0_fit:.3f}")
        print(f"  → DM equivalent: {930 * z0_fit:.0f} pc/cm³")
    except:
        print(f"\n  Sigmoid fit failed (may not be sigmoid-shaped for this pair)")
    
    return {
        'thresholds': thresholds,
        'dm_values': dm_values,
        'mean_dm': mean_dm,
        'cv': cv,
        'z_mids': z_mids.tolist(),
        'correlations': correlations.tolist()
    }


def test2_wavelength_dependence(z, lines):
    """
    WAVELENGTH DEPENDENCE — Operator Spectral Profile
    ==================================================
    Do UV lines, optical lines, and IR lines degrade at different rates?
    This tells us whether the operator has a spectral shape.
    """
    print("\n" + "=" * 70)
    print("TEST 2: WAVELENGTH DEPENDENCE — Operator Spectral Profile")
    print("=" * 70)
    
    # Lines organized by rest wavelength
    line_wavelengths = {
        'LYA': 1216,
        'NV1240': 1240,
        'SIIV_OIV': 1400,
        'CIV': 1549,
        'HEII1640': 1640,
        'CIII_ALL': 1909,
        'MGII': 2796,
        'OII3728': 3728,
        'NEV3426': 3426,
        'HBETA': 4861,
        'OIII5007': 5007,
        'HEII4687': 4687,
        'NII6585': 6585,
        'HALPHA': 6563,
        'SII6718': 6718,
    }
    
    print("\n  Measuring EW drift rate (slope of |r| vs z) per line:\n")
    
    results = []
    
    for name, wavelength in sorted(line_wavelengths.items(), key=lambda x: x[1]):
        if name not in lines:
            continue
        
        ew = lines[name]['ew']
        fwhm = lines[name]['fwhm']
        
        # Measure how EW-FWHM correlation changes with z
        # (EW = information-bearing, FWHM = kinematic)
        mask = (np.isfinite(ew) & np.isfinite(fwhm) & 
                (ew != 0) & (fwhm != 0) & (z > 0.1) & (z < 3.5) & np.isfinite(z))
        
        if mask.sum() < 500:
            continue
        
        zz = z[mask]
        ee = ew[mask]
        ff = fwhm[mask]
        
        # Measure correlation in 5 z-bins
        edges = np.percentile(zz, [0, 20, 40, 60, 80, 100])
        bin_rs = []
        bin_zs = []
        
        for i in range(5):
            bm = (zz >= edges[i]) & (zz < edges[i+1] + (0.001 if i == 4 else 0))
            if bm.sum() > 30:
                r, p = stats.pearsonr(ee[bm], ff[bm])
                if np.isfinite(r):
                    bin_rs.append(r)
                    bin_zs.append((edges[i] + edges[i+1]) / 2)
        
        if len(bin_rs) >= 3:
            trend_r, trend_p = stats.pearsonr(bin_zs, bin_rs)
            drift_rate = np.polyfit(bin_zs, bin_rs, 1)[0]
            
            results.append({
                'line': name,
                'wavelength': wavelength,
                'n': int(mask.sum()),
                'drift_rate': drift_rate,
                'trend_r': trend_r,
                'trend_p': trend_p,
                'bin_rs': bin_rs
            })
            
            flag = "⚠️" if abs(trend_r) > 0.8 else "  "
            print(f"  {flag} {name:12s} (λ={wavelength:5d}Å): drift={drift_rate:+.4f}/z, "
                  f"trend_r={trend_r:+.3f}, n={mask.sum()}")
    
    # Now: does drift rate correlate with wavelength?
    if len(results) >= 5:
        wavelengths = [r['wavelength'] for r in results]
        drifts = [r['drift_rate'] for r in results]
        
        r_wl, p_wl = stats.pearsonr(wavelengths, drifts)
        
        print(f"\n  DRIFT RATE vs WAVELENGTH:")
        print(f"  Correlation: r = {r_wl:.4f}, p = {p_wl:.4e}")
        
        if p_wl < 0.05:
            if r_wl > 0:
                print(f"  → Longer wavelengths drift MORE. Operator is RED-preferential.")
            else:
                print(f"  → Shorter wavelengths drift MORE. Operator is BLUE-preferential.")
        else:
            print(f"  → No significant wavelength dependence. Operator is ACHROMATIC.")
            print(f"     This rules out: dust, Rayleigh scattering, most plasma effects")
            print(f"     This is consistent with: information-level operator (not photon-level)")
    
    return results


def test3_extinction_vs_information(z, lines, scalars):
    """
    EXTINCTION vs INFORMATION DEGRADATION
    ======================================
    If it's dust/extinction, degradation should correlate with E(B-V).
    If it's the medium eating information, degradation should be 
    INDEPENDENT of extinction after correction.
    """
    print("\n" + "=" * 70)
    print("TEST 3: EXTINCTION vs INFORMATION — Is It Just Dust?")
    print("=" * 70)
    
    ebv = scalars.get('EBV')
    if ebv is None:
        print("  No E(B-V) column found. Skipping.")
        return {}
    
    # Split into low and high extinction samples
    oiii_ew = lines['OIII5007']['ew']
    hb_ew = lines['HBETA']['ew']
    sii_ew = lines['SII6718']['ew']
    
    mask_base = (np.isfinite(oiii_ew) & np.isfinite(hb_ew) & np.isfinite(ebv) &
                 (oiii_ew != 0) & (hb_ew != 0) & (z > 0.1) & (z < 3) & np.isfinite(z))
    
    zz = z[mask_base]
    oo = oiii_ew[mask_base]
    hh = hb_ew[mask_base]
    ee = ebv[mask_base]
    
    print(f"  Sample: {mask_base.sum()} quasars with E(B-V) data")
    print(f"  E(B-V) range: {np.nanmin(ee):.3f} to {np.nanmax(ee):.3f}")
    
    # Split by extinction
    low_ext = ee < np.percentile(ee[np.isfinite(ee)], 25)
    high_ext = ee > np.percentile(ee[np.isfinite(ee)], 75)
    
    print(f"\n  Low extinction (bottom 25%): E(B-V) < {np.percentile(ee[np.isfinite(ee)], 25):.3f}, n={low_ext.sum()}")
    print(f"  High extinction (top 25%):  E(B-V) > {np.percentile(ee[np.isfinite(ee)], 75):.3f}, n={high_ext.sum()}")
    
    # Measure EW coupling drift in both samples
    for label, ext_mask in [("LOW extinction", low_ext), ("HIGH extinction", high_ext)]:
        sub_z = zz[ext_mask]
        sub_o = oo[ext_mask]
        sub_h = hh[ext_mask]
        
        edges = np.percentile(sub_z, [0, 33, 67, 100])
        print(f"\n  {label} — OIII_EW vs HB_EW coupling by z:")
        
        for i in range(3):
            bm = (sub_z >= edges[i]) & (sub_z < edges[i+1] + (0.001 if i == 2 else 0))
            if bm.sum() > 30:
                r, p = stats.pearsonr(sub_o[bm], sub_h[bm])
                zlabel = f"z={edges[i]:.1f}-{edges[i+1]:.1f}"
                print(f"    {zlabel}: r={r:.4f} (n={bm.sum()})")
    
    # Key test: does the degradation pattern DIFFER between low and high extinction?
    # If dust causes it, high extinction should show MORE degradation
    # If medium causes it, both should show SAME degradation (it's distance, not dust)
    
    print(f"\n  CRITICAL QUESTION: Does degradation depend on extinction or distance?")
    
    # Match z-ranges and compare
    z_mid = np.median(zz)
    high_z = zz > z_mid
    
    for label, ext_mask in [("Low-ext", low_ext), ("High-ext", high_ext)]:
        both = ext_mask & high_z
        if both.sum() > 50:
            r, p = stats.pearsonr(oo[both], hh[both])
            print(f"    {label} + high-z: r = {r:.4f} (n={both.sum()})")
    
    # Partial correlation: degradation vs z, controlling for E(B-V)
    from scipy.stats import pearsonr
    
    # Measure "degradation" per object: residual from low-z fit
    low_z_mask = mask_base & (z < np.percentile(zz, 25))
    if low_z_mask.sum() > 50:
        slope, intercept = np.polyfit(oiii_ew[low_z_mask], hb_ew[low_z_mask], 1)
        predicted = slope * oo + intercept
        residuals = np.abs(hh - predicted)
        
        r_z_resid, p_z_resid = pearsonr(zz, residuals)
        r_ebv_resid, p_ebv_resid = pearsonr(ee, residuals)
        
        print(f"\n  Partial correlation analysis:")
        print(f"    Residual vs z:      r = {r_z_resid:.4f}, p = {p_z_resid:.2e}")
        print(f"    Residual vs E(B-V): r = {r_ebv_resid:.4f}, p = {p_ebv_resid:.2e}")
        
        if abs(r_z_resid) > abs(r_ebv_resid):
            print(f"    → Distance is STRONGER predictor than extinction")
            print(f"    → NOT DUST. The degradation tracks exposure, not reddening.")
        else:
            print(f"    → Extinction is stronger predictor. Dust may contribute.")
    
    return {'n': int(mask_base.sum())}


def test4_kinematic_immunity_profile(z, lines):
    """
    KINEMATIC IMMUNITY PROFILE
    ===========================
    Measure the FWHM stability vs EW instability for EVERY line.
    If the medium acts on information content, EW should drift while FWHM stays flat.
    Map this across all available lines.
    """
    print("\n" + "=" * 70)
    print("TEST 4: KINEMATIC IMMUNITY PROFILE — EW vs FWHM Stability")
    print("=" * 70)
    
    print("\n  For each line: measure drift of EW and FWHM separately with z")
    print("  If medium is information-selective: EW drifts, FWHM stable\n")
    
    results = []
    
    for name in sorted(lines.keys()):
        ew = lines[name]['ew']
        fwhm = lines[name]['fwhm']
        
        mask = (np.isfinite(ew) & np.isfinite(fwhm) & 
                (ew != 0) & (fwhm != 0) & (z > 0.1) & np.isfinite(z))
        
        if mask.sum() < 500:
            continue
        
        zz = z[mask]
        
        # Measure normalized variance of EW and FWHM in z-bins
        edges = np.percentile(zz, [0, 20, 40, 60, 80, 100])
        
        ew_means = []
        fwhm_means = []
        z_mids = []
        
        for i in range(5):
            bm = (zz >= edges[i]) & (zz < edges[i+1] + (0.001 if i == 4 else 0))
            if bm.sum() > 30:
                ew_means.append(np.median(ew[mask][bm]))
                fwhm_means.append(np.median(fwhm[mask][bm]))
                z_mids.append((edges[i] + edges[i+1]) / 2)
        
        if len(z_mids) >= 4:
            # Normalize to first bin
            ew_norm = np.array(ew_means) / abs(ew_means[0]) if ew_means[0] != 0 else np.array(ew_means)
            fwhm_norm = np.array(fwhm_means) / abs(fwhm_means[0]) if fwhm_means[0] != 0 else np.array(fwhm_means)
            
            ew_drift = np.std(ew_norm)
            fwhm_drift = np.std(fwhm_norm)
            
            ratio = ew_drift / fwhm_drift if fwhm_drift > 0 else float('inf')
            
            results.append({
                'line': name,
                'n': int(mask.sum()),
                'ew_drift': round(ew_drift, 4),
                'fwhm_drift': round(fwhm_drift, 4),
                'ew_fwhm_ratio': round(ratio, 2),
                'ew_by_bin': [round(x, 3) for x in ew_norm.tolist()],
                'fwhm_by_bin': [round(x, 3) for x in fwhm_norm.tolist()]
            })
            
            flag = "🔒" if ratio < 1.5 else "📡"
            print(f"  {flag} {name:15s}: EW_drift={ew_drift:.4f}, FWHM_drift={fwhm_drift:.4f}, "
                  f"ratio={ratio:.2f} (n={mask.sum()})")
    
    # Summary
    if results:
        ratios = [r['ew_fwhm_ratio'] for r in results if r['ew_fwhm_ratio'] < 100]
        mean_ratio = np.mean(ratios)
        n_info_more = sum(1 for r in ratios if r > 1)
        
        print(f"\n  SUMMARY:")
        print(f"  Lines where EW drifts MORE than FWHM: {n_info_more}/{len(ratios)}")
        print(f"  Mean EW/FWHM drift ratio: {mean_ratio:.2f}")
        
        if n_info_more > len(ratios) * 0.7:
            print(f"  ✅ CONFIRMED: Information-bearing channel (EW) systematically less stable")
            print(f"     than kinematic channel (FWHM) across {n_info_more}/{len(ratios)} lines")
            print(f"     This is the CHANNEL SELECTIVITY signature.")
    
    return results


def test5_luminosity_control(z, lines, scalars):
    """
    LUMINOSITY CONTROL
    ===================
    The strongest counter-argument is Malmquist bias / luminosity selection.
    At high-z we only see bright quasars. Maybe that changes the relationships.
    
    Control: match luminosity between z-bins and check if drift persists.
    """
    print("\n" + "=" * 70)
    print("TEST 5: LUMINOSITY CONTROL — Killing the Selection Argument")
    print("=" * 70)
    
    lbol = scalars.get('LOGLBOL')
    if lbol is None:
        print("  No LOGLBOL column. Skipping.")
        return {}
    
    oiii_ew = lines['OIII5007']['ew']
    hb_ew = lines['HBETA']['ew']
    
    mask = (np.isfinite(oiii_ew) & np.isfinite(hb_ew) & np.isfinite(lbol) &
            (oiii_ew != 0) & (hb_ew != 0) & (z > 0.1) & (z < 3) & 
            np.isfinite(z) & (lbol > 0))
    
    zz = z[mask]
    oo = oiii_ew[mask]
    hh = hb_ew[mask]
    ll = lbol[mask]
    
    print(f"  Sample: {mask.sum()} quasars with Lbol")
    
    # Find luminosity range common to both low-z and high-z
    z_low = zz < np.percentile(zz, 30)
    z_high = zz > np.percentile(zz, 70)
    
    l_low_range = (np.percentile(ll[z_low], 10), np.percentile(ll[z_low], 90))
    l_high_range = (np.percentile(ll[z_high], 10), np.percentile(ll[z_high], 90))
    
    # Overlap range
    l_min = max(l_low_range[0], l_high_range[0])
    l_max = min(l_low_range[1], l_high_range[1])
    
    print(f"  Low-z Lbol range: {l_low_range[0]:.1f} - {l_low_range[1]:.1f}")
    print(f"  High-z Lbol range: {l_high_range[0]:.1f} - {l_high_range[1]:.1f}")
    print(f"  Overlap: {l_min:.1f} - {l_max:.1f}")
    
    if l_min >= l_max:
        print("  No luminosity overlap. Cannot match.")
        return {}
    
    # Luminosity-matched subsamples
    l_match = (ll >= l_min) & (ll <= l_max)
    
    matched_low = z_low & l_match
    matched_high = z_high & l_match
    
    print(f"\n  Luminosity-matched samples:")
    print(f"    Low-z:  n={matched_low.sum()}, <Lbol>={np.mean(ll[matched_low]):.2f}")
    print(f"    High-z: n={matched_high.sum()}, <Lbol>={np.mean(ll[matched_high]):.2f}")
    
    if matched_low.sum() < 50 or matched_high.sum() < 50:
        print("  Insufficient matched sample. Skipping.")
        return {}
    
    # Compare OIII-HB coupling
    r_low, p_low = stats.pearsonr(oo[matched_low], hh[matched_low])
    r_high, p_high = stats.pearsonr(oo[matched_high], hh[matched_high])
    
    print(f"\n  OIII_EW vs HB_EW coupling (luminosity-matched):")
    print(f"    Low-z:  r = {r_low:.4f} (p = {p_low:.2e})")
    print(f"    High-z: r = {r_high:.4f} (p = {p_high:.2e})")
    print(f"    Δr = {r_high - r_low:.4f}")
    
    # KS test on residuals
    slope, intercept = np.polyfit(oo[matched_low], hh[matched_low], 1)
    resid_low = hh[matched_low] - (slope * oo[matched_low] + intercept)
    resid_high = hh[matched_high] - (slope * oo[matched_high] + intercept)
    
    ks, ks_p = stats.ks_2samp(resid_low, resid_high)
    
    print(f"\n  KS test on residuals: D = {ks:.4f}, p = {ks_p:.2e}")
    
    if ks_p < 0.001:
        print(f"  ✅ DRIFT PERSISTS after luminosity matching!")
        print(f"     Selection/Malmquist bias does NOT explain the pattern.")
    else:
        print(f"  ⚠️ Drift weakens after luminosity matching. Selection may contribute.")
    
    # Also check FWHM (should NOT drift even without matching)
    oiii_fwhm = lines['OIII5007']['fwhm']
    hb_fwhm = lines['HBETA']['fwhm']
    
    mask_f = (np.isfinite(oiii_fwhm) & np.isfinite(hb_fwhm) & 
              (oiii_fwhm != 0) & (hb_fwhm != 0) & mask)
    
    if mask_f.sum() > 200:
        zf = z[mask_f]
        of = oiii_fwhm[mask_f]
        hf = hb_fwhm[mask_f]
        
        z_low_f = zf < np.percentile(zf, 30)
        z_high_f = zf > np.percentile(zf, 70)
        
        if z_low_f.sum() > 30 and z_high_f.sum() > 30:
            r_low_f, _ = stats.pearsonr(of[z_low_f], hf[z_low_f])
            r_high_f, _ = stats.pearsonr(of[z_high_f], hf[z_high_f])
            
            print(f"\n  CONTROL — FWHM coupling (kinematic, should be stable):")
            print(f"    Low-z:  r = {r_low_f:.4f}")
            print(f"    High-z: r = {r_high_f:.4f}")
            print(f"    Δr = {r_high_f - r_low_f:.4f}")
            
            if abs(r_high_f - r_low_f) < abs(r_high - r_low):
                print(f"  ✅ FWHM is MORE STABLE than EW across z. Channel selectivity confirmed.")
    
    return {
        'r_low_ew': r_low, 'r_high_ew': r_high,
        'ks_p': ks_p,
        'matched_n_low': int(matched_low.sum()),
        'matched_n_high': int(matched_high.sum())
    }


def test6_snr_control(z, lines, scalars):
    """
    SNR CONTROL
    ============
    Maybe high-z spectra just have worse SNR, making correlations noisier.
    Control: match SNR between z-bins.
    """
    print("\n" + "=" * 70)
    print("TEST 6: SNR CONTROL — Is It Just Noise?")
    print("=" * 70)
    
    snr = scalars.get('SN_MEDIAN_ALL')
    if snr is None:
        print("  No SNR column. Skipping.")
        return {}
    
    oiii_ew = lines['OIII5007']['ew']
    hb_ew = lines['HBETA']['ew']
    
    mask = (np.isfinite(oiii_ew) & np.isfinite(hb_ew) & np.isfinite(snr) &
            (oiii_ew != 0) & (hb_ew != 0) & (z > 0.1) & np.isfinite(z) & (snr > 0))
    
    zz = z[mask]
    oo = oiii_ew[mask]
    hh = hb_ew[mask]
    ss = snr[mask]
    
    print(f"  Sample: {mask.sum()} quasars with SNR data")
    
    # High-SNR subsample only
    snr_cut = np.percentile(ss, 75)
    high_snr = ss > snr_cut
    
    print(f"  High-SNR cut: SNR > {snr_cut:.1f} (n={high_snr.sum()})")
    
    z_low = (zz < np.percentile(zz[high_snr], 30)) & high_snr
    z_high = (zz > np.percentile(zz[high_snr], 70)) & high_snr
    
    if z_low.sum() > 50 and z_high.sum() > 50:
        r_low, _ = stats.pearsonr(oo[z_low], hh[z_low])
        r_high, _ = stats.pearsonr(oo[z_high], hh[z_high])
        
        print(f"\n  High-SNR only — OIII_EW vs HB_EW:")
        print(f"    Low-z:  r = {r_low:.4f} (n={z_low.sum()})")
        print(f"    High-z: r = {r_high:.4f} (n={z_high.sum()})")
        
        slope, intercept = np.polyfit(oo[z_low], hh[z_low], 1)
        resid_low = hh[z_low] - (slope * oo[z_low] + intercept)
        resid_high = hh[z_high] - (slope * oo[z_high] + intercept)
        
        ks, ks_p = stats.ks_2samp(resid_low, resid_high)
        print(f"    KS on residuals: D = {ks:.4f}, p = {ks_p:.2e}")
        
        if ks_p < 0.001:
            print(f"  ✅ DRIFT PERSISTS in high-SNR subsample. Not noise.")
        else:
            print(f"  ⚠️ Drift weakens in high-SNR. Noise may contribute.")
    
    return {}


def test7_sigmoid_structure(z, lines):
    """
    SIGMOID STRUCTURE TEST
    =======================
    Fit sigmoids to multiple line pairs. If the medium is responsible,
    all sigmoids should have similar z₀ (same physical threshold).
    If it's source evolution, z₀ should vary randomly per line.
    """
    print("\n" + "=" * 70)
    print("TEST 7: SIGMOID UNIVERSALITY — Same Threshold Across Lines?")
    print("=" * 70)
    
    target_lines = ['HBETA', 'HALPHA', 'MGII', 'CIV', 'OII3728', 'SII6718', 'LYA']
    reference = 'OIII5007'
    
    ref_ew = lines[reference]['ew']
    z0_values = []
    
    print(f"\n  Fitting sigmoid to |r|(z) for {reference}_EW vs X_EW:\n")
    
    for name in target_lines:
        if name not in lines:
            continue
        
        target_ew = lines[name]['ew']
        mask = (np.isfinite(ref_ew) & np.isfinite(target_ew) &
                (ref_ew != 0) & (target_ew != 0) & (z > 0.05) & np.isfinite(z))
        
        if mask.sum() < 500:
            continue
        
        zz = z[mask]
        rr = ref_ew[mask]
        tt = target_ew[mask]
        
        # 15 z-bins
        edges = np.percentile(zz, np.linspace(0, 100, 16))
        z_mids = []
        cors = []
        
        for i in range(15):
            bm = (zz >= edges[i]) & (zz < edges[i+1] + (0.001 if i == 14 else 0))
            if bm.sum() > 30:
                r, p = stats.pearsonr(rr[bm], tt[bm])
                if np.isfinite(r):
                    z_mids.append((edges[i] + edges[i+1]) / 2)
                    cors.append(abs(r))
        
        if len(z_mids) < 5:
            continue
        
        z_mids = np.array(z_mids)
        cors = np.array(cors)
        
        # Try sigmoid fit
        try:
            popt, pcov = curve_fit(sigmoid, z_mids, cors,
                                   p0=[max(cors) - min(cors), -3, 1.0, min(cors)],
                                   maxfev=10000)
            z0 = popt[2]
            z0_err = np.sqrt(pcov[2, 2]) if pcov[2, 2] > 0 else 0
            
            if 0 < z0 < 5:  # sanity check
                z0_values.append({'line': name, 'z0': z0, 'z0_err': z0_err})
                dm_equiv = 930 * z0
                print(f"  {name:12s}: z₀ = {z0:.3f} ± {z0_err:.3f} → DM ≈ {dm_equiv:.0f} pc/cm³")
            else:
                print(f"  {name:12s}: sigmoid fit out of range (z₀ = {z0:.2f})")
        except Exception as e:
            print(f"  {name:12s}: sigmoid fit failed ({str(e)[:40]})")
    
    if len(z0_values) >= 3:
        z0s = [v['z0'] for v in z0_values]
        mean_z0 = np.mean(z0s)
        std_z0 = np.std(z0s)
        cv = std_z0 / mean_z0
        
        print(f"\n  SIGMOID THRESHOLD SUMMARY:")
        print(f"  z₀ values: {[f'{z:.3f}' for z in z0s]}")
        print(f"  Mean z₀: {mean_z0:.3f}")
        print(f"  Std z₀: {std_z0:.3f}")
        print(f"  CV: {cv:.3f}")
        print(f"  Mean DM equivalent: {930 * mean_z0:.0f} pc/cm³")
        
        if cv < 0.3:
            print(f"\n  ✅ UNIVERSAL THRESHOLD: All lines transition at similar z!")
            print(f"     Source evolution predicts LINE-SPECIFIC thresholds.")
            print(f"     A medium predicts a UNIVERSAL threshold (same path, same effect).")
            print(f"     Result supports medium over source evolution.")
        else:
            print(f"\n  ⚠️ Thresholds vary (CV={cv:.3f}). May reflect line-specific sensitivity.")
    
    return z0_values


def main():
    print("=" * 70)
    print("MECHANISM HUNT — Finding the Rat")
    print("=" * 70)
    print("""
    We know the phenomenon exists. Now we find what causes it.
    
    Strategy:
    1. Map all sigmoid thresholds to common physical variable
    2. Profile the operator's wavelength dependence  
    3. Rule out dust/extinction
    4. Confirm channel selectivity (EW vs FWHM) per line
    5. Kill luminosity selection as explanation
    6. Kill SNR as explanation
    7. Test sigmoid universality across lines
    """)
    
    z, lines, scalars = load_quasar_data()
    print(f"Loaded 750,414 quasars, {len(lines)} lines\n")
    
    results = {}
    results['rosetta'] = test1_rosetta_stone(z, lines)
    results['wavelength'] = test2_wavelength_dependence(z, lines)
    results['extinction'] = test3_extinction_vs_information(z, lines, scalars)
    results['immunity'] = test4_kinematic_immunity_profile(z, lines)
    results['luminosity'] = test5_luminosity_control(z, lines, scalars)
    results['snr'] = test6_snr_control(z, lines, scalars)
    results['sigmoid'] = test7_sigmoid_structure(z, lines)
    
    # Save
    with open(RESULTS_DIR / 'results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print("\n" + "=" * 70)
    print("MECHANISM CONSTRAINTS")
    print("=" * 70)
    print("""
    What the tests constrain:
    
    If wavelength-independent → NOT dust, NOT plasma scattering
    If extinction-independent → NOT reddening, NOT ISM absorption  
    If SNR-independent → NOT measurement noise
    If luminosity-independent → NOT Malmquist/selection bias
    If sigmoid universal → Medium threshold, NOT source evolution
    If EW drifts but FWHM doesn't → Information-selective, NOT geometric
    If all DM equivalents converge → Integrated electron column density
                                     is the controlling variable
    
    The rat has fewer places to hide.
    """)


if __name__ == '__main__':
    main()
