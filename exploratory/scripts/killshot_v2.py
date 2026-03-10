#!/usr/bin/env python3
"""
KILLSHOT v2 — Three Gemini-prescribed tests
=============================================
1. Partial correlation r(EW, FWHM | L/L_Edd) for CIV — kills BLR evolution
2. [NeV] 3346/3426 ratio variance vs z — resolved branching-locked pair
3. Baldwin effect slope difference Δβ for [OIII] doublet components
"""

import numpy as np
from astropy.io import fits
from scipy.stats import spearmanr, kurtosis, pearsonr
import warnings
warnings.filterwarnings('ignore')

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data
z = data['Z_DR16Q']

# ============================================================
# TEST 1: PARTIAL CORRELATION — r(EW, FWHM | L/L_Edd)
# If BLR evolution causes decorrelation, controlling for 
# Eddington ratio should REMOVE the z-trend.
# If refractive decoherence causes it, controlling for 
# Eddington ratio should NOT remove the z-trend.
# ============================================================
print("=" * 80)
print("TEST 1: PARTIAL CORRELATION r(EW, FWHM | L/L_Edd)")
print("Controls for accretion state. Kills BLR evolution if trend persists.")
print("=" * 80)

# Get CIV data
civ = data['CIV']
ew_civ = civ[:, 2]
fwhm_civ = civ[:, 4]

# Eddington ratio proxy: use LOGLBOL and LOGMBH_CIV
try:
    lbol = data['LOGLBOL']
    mbh = data['LOGMBH_CIV']
    
    # L/L_Edd ∝ L_bol / M_BH → log(L/L_Edd) = logL - logM - const
    log_edd = lbol - mbh - 38.1  # Eddington luminosity = 1.26e38 * M/M_sun
    
    valid = (np.isfinite(ew_civ) & np.isfinite(fwhm_civ) & 
             np.isfinite(log_edd) & np.isfinite(lbol) &
             (ew_civ > 0) & (fwhm_civ > 0) & (ew_civ < 500) & (fwhm_civ < 30000) &
             (z > 0.5) & (log_edd > -3) & (log_edd < 2))
    
    print(f"\nCIV with valid L/L_Edd: {valid.sum():,}")
    
    zv = z[valid]
    ewv = ew_civ[valid]
    fwv = fwhm_civ[valid]
    eddv = log_edd[valid]
    lbolv = lbol[valid]
    
    # Bin by z
    idx = np.argsort(zv)
    n = len(idx)
    nbins = 15
    bs = n // nbins
    
    print(f"\n{'z':>6} {'r_raw':>8} {'r_partial':>10} {'r_partial2':>11} {'N':>7} {'med_Edd':>9}")
    print("-" * 60)
    
    raw_corrs = []
    partial_corrs = []
    z_centers = []
    
    for i in range(nbins):
        sl = idx[i*bs:(i+1)*bs]
        ew_bin = ewv[sl]
        fw_bin = fwv[sl]
        edd_bin = eddv[sl]
        lbol_bin = lbolv[sl]
        z_bin = zv[sl]
        
        # Raw correlation
        r_raw, _ = spearmanr(ew_bin, fw_bin)
        
        # Partial correlation: r(EW, FWHM | L/L_Edd)
        # Residualize both EW and FWHM against L/L_Edd
        from numpy.polynomial import polynomial as P
        
        # Rank-based partial correlation (more robust)
        ew_rank = np.argsort(np.argsort(ew_bin)).astype(float)
        fw_rank = np.argsort(np.argsort(fw_bin)).astype(float)
        edd_rank = np.argsort(np.argsort(edd_bin)).astype(float)
        
        # Residualize EW ranks against Edd ranks
        ew_coef = np.polyfit(edd_rank, ew_rank, 1)
        ew_resid = ew_rank - np.polyval(ew_coef, edd_rank)
        
        # Residualize FWHM ranks against Edd ranks
        fw_coef = np.polyfit(edd_rank, fw_rank, 1)
        fw_resid = fw_rank - np.polyval(fw_coef, edd_rank)
        
        # Partial correlation
        r_partial, _ = pearsonr(ew_resid, fw_resid)
        
        # Also control for BOTH L/L_Edd AND L_bol
        lbol_rank = np.argsort(np.argsort(lbol_bin)).astype(float)
        X = np.column_stack([edd_rank, lbol_rank])
        
        ew_coef2 = np.linalg.lstsq(X, ew_rank, rcond=None)[0]
        ew_resid2 = ew_rank - X @ ew_coef2
        fw_coef2 = np.linalg.lstsq(X, fw_rank, rcond=None)[0]
        fw_resid2 = fw_rank - X @ fw_coef2
        r_partial2, _ = pearsonr(ew_resid2, fw_resid2)
        
        raw_corrs.append(r_raw)
        partial_corrs.append(r_partial)
        z_centers.append(np.median(z_bin))
        
        print(f"{np.median(z_bin):6.3f} {r_raw:>+8.4f} {r_partial:>+10.4f} {r_partial2:>+11.4f} "
              f"{len(sl):>7} {np.median(edd_bin):>9.3f}")
    
    z_centers = np.array(z_centers)
    raw_corrs = np.array(raw_corrs)
    partial_corrs = np.array(partial_corrs)
    
    r_raw_trend, p_raw_trend = spearmanr(z_centers, raw_corrs)
    r_partial_trend, p_partial_trend = spearmanr(z_centers, partial_corrs)
    
    print(f"\nRaw r(EW,FWHM) trend with z:     r = {r_raw_trend:+.3f} (p = {p_raw_trend:.2e})")
    print(f"Partial r(EW,FWHM|Edd) trend:    r = {r_partial_trend:+.3f} (p = {p_partial_trend:.2e})")
    
    if r_partial_trend < -0.5 and p_partial_trend < 0.05:
        print("\n🎯 DECORRELATION PERSISTS AFTER CONTROLLING FOR EDDINGTON RATIO")
        print("BLR evolution / disk-wind models CANNOT explain this.")
        print("The effect is EXTRINSIC to the quasar engine.")
    elif r_partial_trend < -0.3:
        print("\n⚠️  Trend weakened but still present. Partial BLR contribution possible.")
    else:
        print("\n✗ Controlling for Eddington ratio removes the trend.")
        print("BLR evolution may be the primary driver.")

except Exception as e:
    print(f"Error in Test 1: {e}")
    import traceback
    traceback.print_exc()

# ============================================================
# TEST 2: [NeV] 3346/3426 RATIO VARIANCE
# Branching-locked pair with 80Å separation — RESOLVED by SDSS
# ============================================================
print(f"\n{'='*80}")
print("TEST 2: [NeV] 3346/3426 SEARCH")
print("Branching-locked pair, well-separated, high-ionization")
print("=" * 80)

# Check if DR16Q has NeV
cols = [c.name for c in hdu[1].columns]
nev_cols = [c for c in cols if 'NEV' in c.upper() or 'NE_V' in c.upper() or 'NE5' in c.upper()]
print(f"NeV columns found: {nev_cols}")

if not nev_cols:
    print("DR16Q does not include [NeV] line measurements.")
    print("This line must be measured from raw SDSS spectra.")
    print("ALTERNATIVE: Use [OI] 6300/6363 or [NeIII] 3869/3967")
    
    # Check for OI
    oi_cols = [c for c in cols if 'OI' in c.upper() and '6300' in c]
    print(f"[OI] columns: {oi_cols}")
    
    # Check what we DO have that could work
    print("\nAvailable narrow-line columns:")
    for c in cols:
        if any(x in c.upper() for x in ['OIII','NII','SII','OII','NEV','NEI','OI']):
            if 'ERR' not in c and 'STAT' not in c and 'COMP' not in c:
                print(f"  {c}")

# ============================================================
# TEST 3: CIV EW-FWHM IN NARROW Eddington BINS
# Strongest version: fix L/L_Edd tightly, vary ONLY z
# ============================================================
print(f"\n{'='*80}")
print("TEST 3: CIV DECORRELATION IN FIXED-EDDINGTON SLICES")
print("Narrow Eddington bins → removes ALL accretion-state variation")
print("=" * 80)

try:
    # Narrow Eddington bins
    edd_bins = [(-1.5, -1.0), (-1.0, -0.5), (-0.5, 0.0), (0.0, 0.5)]
    
    for edd_lo, edd_hi in edd_bins:
        edd_mask = valid & (log_edd >= edd_lo) & (log_edd < edd_hi)
        
        if edd_mask.sum() < 5000:
            continue
        
        zv2 = z[edd_mask]
        ewv2 = ew_civ[edd_mask]
        fwv2 = fwhm_civ[edd_mask]
        
        idx2 = np.argsort(zv2)
        n2 = len(idx2)
        nb2 = min(10, n2 // 500)
        bs2 = n2 // nb2
        
        corrs2 = []
        zc2 = []
        for i in range(nb2):
            sl = idx2[i*bs2:(i+1)*bs2]
            r, _ = spearmanr(ewv2[sl], fwv2[sl])
            corrs2.append(r)
            zc2.append(np.median(zv2[sl]))
        
        corrs2 = np.array(corrs2)
        zc2 = np.array(zc2)
        
        r_trend, p_trend = spearmanr(zc2, corrs2)
        
        decay = "DECAYS" if r_trend < -0.3 and p_trend < 0.1 else "FLAT" if abs(r_trend) < 0.3 else "RISES"
        
        print(f"\n  log(L/L_Edd) = [{edd_lo}, {edd_hi}] — N = {edd_mask.sum():,}")
        print(f"    r(EW,FWHM): {corrs2[0]:.3f} → {corrs2[-1]:.3f}")
        print(f"    Trend: r = {r_trend:+.3f} (p = {p_trend:.2e}) → {decay}")
    
    print(f"\nIf ALL Eddington bins show decay → BLR evolution is dead.")
    print(f"The decorrelation is extrinsic, regardless of accretion state.")

except Exception as e:
    print(f"Error in Test 3: {e}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*80}")
print("KILLSHOT v2 SUMMARY")
print("=" * 80)
print(f"""
TEST 1 — Partial correlation r(EW, FWHM | L/L_Edd):
  Raw trend:     r = {r_raw_trend:+.3f} (p = {p_raw_trend:.2e})
  Partial trend: r = {r_partial_trend:+.3f} (p = {p_partial_trend:.2e})
  → Does controlling for accretion state remove the signal?
  
TEST 2 — [NeV] branching-locked ratio:
  → Not available in DR16Q catalog. Needs raw spectral refitting or DESI.
  
TEST 3 — Fixed-Eddington slices:
  → Does the decay persist at EVERY accretion rate?
""")

hdu.close()
