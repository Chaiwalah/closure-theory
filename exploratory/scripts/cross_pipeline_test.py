#!/usr/bin/env python3
"""
CROSS-PIPELINE TEST — C1 KILLER
=================================
Compare DR16Q (SDSS pipeline) vs Rakshit+2020 (PyQSOFit pipeline)
on the SAME quasars. If the EW-FWHM asymmetry exists in BOTH
independent pipelines, C1 (pipeline artifact) is DEAD.

Rakshit+2020: 526,265 quasars, PyQSOFit fitting (independent of SDSS pipeline)
DR16Q: 750,414 quasars, SDSS pipeline fitting

Cross-match by RA/Dec within 1 arcsec.
"""

import numpy as np
from scipy import stats
from astropy.io import fits
import gzip
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = Path('results_cross_pipeline')
RESULTS_DIR.mkdir(exist_ok=True)

print("=" * 80)
print("🔪 CROSS-PIPELINE TEST — C1 KILLER")
print("DR16Q (SDSS pipeline) vs Rakshit+2020 (PyQSOFit)")
print("=" * 80)

# =====================================================================
# LOAD DR16Q
# =====================================================================
print("\nLoading DR16Q...")
f = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
d = f[1].data
dr16_ra = d['RA']
dr16_dec = d['DEC']
dr16_z = d['Z_DR16Q']

dr16_lines = {}
for name in ['CIV', 'CIII_ALL', 'MGII', 'HBETA', 'OIII5007', 'LYA']:
    col = d[name]
    if col.ndim == 2 and col.shape[1] >= 6:
        dr16_lines[name] = {'ew': col[:, 2], 'fwhm': col[:, 4]}
f.close()
print(f"  DR16Q: {len(dr16_ra)} quasars")

# =====================================================================
# LOAD RAKSHIT+2020 (fixed-width format)
# =====================================================================
print("\nLoading Rakshit+2020 (this takes a minute)...")

# Column positions (0-indexed) from ReadMe
# RA: 39-48, Dec: 50-59, z: 95-102
# CIV EW: 2613-2624, CIV FWHM: 2561-2572
# CIII EW: 2451-2462, CIII FWHM: 2399-2410
# MgII broad EW: 2289-2300, MgII broad FWHM: 2237-2248
# Hβ broad EW: 1393-1404 (need to find), Hβ broad FWHM: ...
# OIII5007 core EW: 1341-1352 (OIII4959C) ... need OIII5007C
# Lyα EW: 2775-2786, Lyα FWHM: 2723-2734

# Let me parse more carefully from ReadMe
# Hbeta broad: 
#   FWHMHb-BR: bytes 1453-1464
#   EWHb-BR: bytes 1493-1504
# OIII5007C (core):
#   EWOIII5007C: bytes 1367-1378

col_defs = {
    'ra':        (38, 48),      # F10.6
    'dec':       (49, 59),      # F10.6
    'z':         (94, 102),     # F8.6
    'snr':       (103, 111),    # F8.3
    'ew_civ':    (2612, 2624),  # E12.6
    'fwhm_civ':  (2560, 2572),  # F12.6
    'ew_ciii':   (2450, 2462),  # E12.6
    'fwhm_ciii': (2398, 2410),  # F12.6
    'ew_mgii':   (2288, 2300),  # E12.6 (broad)
    'fwhm_mgii': (2236, 2248),  # F12.6 (broad)
    'ew_lya':    (2774, 2786),  # E12.6
    'fwhm_lya':  (2722, 2734),  # F12.6
}

rak = {k: [] for k in col_defs}

dat_path = '/root/clawd/data/sdss/rakshit2020.dat.gz'
with gzip.open(dat_path, 'rt', encoding='ascii', errors='replace') as fh:
    for i, line in enumerate(fh):
        if i % 100000 == 0:
            print(f"  Read {i} lines...")
        for col_name, (start, end) in col_defs.items():
            try:
                val = float(line[start:end].strip())
            except (ValueError, IndexError):
                val = np.nan
            rak[col_name].append(val)

for k in rak:
    rak[k] = np.array(rak[k])

# Replace -999 sentinels with NaN
for k in rak:
    if k not in ('ra', 'dec', 'z'):
        rak[k][rak[k] <= -998] = np.nan

n_rak = len(rak['ra'])
print(f"  Rakshit+2020: {n_rak} quasars")

# =====================================================================
# CROSS-MATCH by RA/Dec (1 arcsec = 1/3600 deg)
# =====================================================================
print("\nCross-matching by position (1 arcsec radius)...")
# For speed, use binning approach
from collections import defaultdict

# Bin DR16Q by rounded RA/Dec
bin_size = 0.01  # ~36 arcsec bins
dr16_bins = defaultdict(list)
for i in range(len(dr16_ra)):
    if np.isfinite(dr16_ra[i]) and np.isfinite(dr16_dec[i]):
        key = (round(dr16_ra[i] / bin_size), round(dr16_dec[i] / bin_size))
        dr16_bins[key].append(i)

match_tol = 1.0 / 3600  # 1 arcsec in degrees
matched_dr16 = []
matched_rak = []

for j in range(n_rak):
    ra_j = rak['ra'][j]
    dec_j = rak['dec'][j]
    if not (np.isfinite(ra_j) and np.isfinite(dec_j)):
        continue
    
    key = (round(ra_j / bin_size), round(dec_j / bin_size))
    
    # Check nearby bins
    best_dist = 999
    best_i = -1
    for dk_ra in (-1, 0, 1):
        for dk_dec in (-1, 0, 1):
            nkey = (key[0] + dk_ra, key[1] + dk_dec)
            for i in dr16_bins.get(nkey, []):
                dist = np.sqrt((dr16_ra[i] - ra_j)**2 + 
                              ((dr16_dec[i] - dec_j) * np.cos(np.radians(dec_j)))**2)
                if dist < best_dist:
                    best_dist = dist
                    best_i = i
    
    if best_dist < match_tol and best_i >= 0:
        matched_dr16.append(best_i)
        matched_rak.append(j)

matched_dr16 = np.array(matched_dr16)
matched_rak = np.array(matched_rak)
print(f"  Matched: {len(matched_dr16)} quasars")

# Verify match quality
dz = np.abs(dr16_z[matched_dr16] - rak['z'][matched_rak])
print(f"  Redshift agreement: median Δz = {np.nanmedian(dz):.6f}")
good_z = dz < 0.01
print(f"  Good redshift match (Δz < 0.01): {good_z.sum()}")

# Use only good matches
m_dr16 = matched_dr16[good_z]
m_rak = matched_rak[good_z]
z_matched = dr16_z[m_dr16]

print(f"\n  Final matched sample: {len(m_dr16)} quasars")

# =====================================================================
# THE TEST: EW-FWHM decorrelation in BOTH pipelines
# =====================================================================
print(f"\n{'='*80}")
print("THE TEST: EW-FWHM correlation by z-bin in BOTH pipelines")
print("If the EW≈0 / FWHM>0 asymmetry appears in BOTH → C1 is DEAD")
print("=" * 80)

z_bins = [(0.5, 0.8), (0.8, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.5)]

# Compare CIV
print(f"\n--- CIV (EW vs FWHM correlation by z-bin) ---")
print(f"  {'z-bin':<10} {'DR16Q r(EW,FW)':<18} {'Rakshit r(EW,FW)':<20} {'Agreement'}")
print(f"  {'-'*65}")

for z_lo, z_hi in z_bins:
    z_mask = (z_matched > z_lo) & (z_matched < z_hi)
    
    # DR16Q
    ew_dr16 = dr16_lines['CIV']['ew'][m_dr16]
    fw_dr16 = dr16_lines['CIV']['fwhm'][m_dr16]
    valid_dr16 = z_mask & np.isfinite(ew_dr16) & np.isfinite(fw_dr16) & (ew_dr16 != 0) & (fw_dr16 > 0)
    
    # Rakshit
    ew_rak = rak['ew_civ'][m_rak]
    fw_rak = rak['fwhm_civ'][m_rak]
    valid_rak = z_mask & np.isfinite(ew_rak) & np.isfinite(fw_rak) & (ew_rak != 0) & (fw_rak > 0)
    
    r_dr16 = r_rak = "N/A"
    if valid_dr16.sum() >= 30:
        r_d, p_d = stats.spearmanr(ew_dr16[valid_dr16], fw_dr16[valid_dr16])
        r_dr16 = f"{r_d:+.4f} (N={valid_dr16.sum()})"
    if valid_rak.sum() >= 30:
        r_r, p_r = stats.spearmanr(ew_rak[valid_rak], fw_rak[valid_rak])
        r_rak = f"{r_r:+.4f} (N={valid_rak.sum()})"
    
    agree = ""
    if valid_dr16.sum() >= 30 and valid_rak.sum() >= 30:
        if abs(r_d - r_r) < 0.1:
            agree = "✅ AGREE"
        else:
            agree = f"⚠️ Δ={abs(r_d-r_r):.3f}"
    
    print(f"  {z_lo:.1f}-{z_hi:.1f}    {r_dr16:<18} {r_rak:<20} {agree}")

# Compare MgII
print(f"\n--- MgII (EW vs FWHM correlation by z-bin) ---")
print(f"  {'z-bin':<10} {'DR16Q r(EW,FW)':<18} {'Rakshit r(EW,FW)':<20} {'Agreement'}")
print(f"  {'-'*65}")

for z_lo, z_hi in z_bins:
    z_mask = (z_matched > z_lo) & (z_matched < z_hi)
    
    ew_dr16 = dr16_lines['MGII']['ew'][m_dr16]
    fw_dr16 = dr16_lines['MGII']['fwhm'][m_dr16]
    valid_dr16 = z_mask & np.isfinite(ew_dr16) & np.isfinite(fw_dr16) & (ew_dr16 != 0) & (fw_dr16 > 0)
    
    ew_rak = rak['ew_mgii'][m_rak]
    fw_rak = rak['fwhm_mgii'][m_rak]
    valid_rak = z_mask & np.isfinite(ew_rak) & np.isfinite(fw_rak) & (ew_rak != 0) & (fw_rak > 0)
    
    r_dr16 = r_rak = "N/A"
    if valid_dr16.sum() >= 30:
        r_d, p_d = stats.spearmanr(ew_dr16[valid_dr16], fw_dr16[valid_dr16])
        r_dr16 = f"{r_d:+.4f} (N={valid_dr16.sum()})"
    if valid_rak.sum() >= 30:
        r_r, p_r = stats.spearmanr(ew_rak[valid_rak], fw_rak[valid_rak])
        r_rak = f"{r_r:+.4f} (N={valid_rak.sum()})"
    
    agree = ""
    if valid_dr16.sum() >= 30 and valid_rak.sum() >= 30:
        if abs(r_d - r_r) < 0.1:
            agree = "✅ AGREE"
        else:
            agree = f"⚠️ Δ={abs(r_d-r_r):.3f}"
    
    print(f"  {z_lo:.1f}-{z_hi:.1f}    {r_dr16:<18} {r_rak:<20} {agree}")

# =====================================================================
# THE INTER-LINE TEST: CIV-CIII EW coupling in both pipelines
# =====================================================================
print(f"\n{'='*80}")
print("INTER-LINE TEST: CIV-CIII EW coupling (both pipelines)")
print("This is the core test: does EW inter-line correlation degrade with z?")
print("=" * 80)

print(f"\n  {'z-bin':<10} {'DR16Q r(EW)':<15} {'Rakshit r(EW)':<17} {'Agreement'}")
print(f"  {'-'*55}")

for z_lo, z_hi in z_bins:
    z_mask = (z_matched > z_lo) & (z_matched < z_hi)
    
    # DR16Q: CIV-CIII EW
    ew_civ_d = dr16_lines['CIV']['ew'][m_dr16]
    ew_ciii_d = dr16_lines['CIII_ALL']['ew'][m_dr16]
    valid_d = z_mask & np.isfinite(ew_civ_d) & np.isfinite(ew_ciii_d) & (ew_civ_d != 0) & (ew_ciii_d != 0)
    
    # Rakshit: CIV-CIII EW
    ew_civ_r = rak['ew_civ'][m_rak]
    ew_ciii_r = rak['ew_ciii'][m_rak]
    valid_r = z_mask & np.isfinite(ew_civ_r) & np.isfinite(ew_ciii_r) & (ew_civ_r != 0) & (ew_ciii_r != 0)
    
    r_d_str = r_r_str = "N/A"
    if valid_d.sum() >= 30:
        r_d, _ = stats.spearmanr(ew_civ_d[valid_d], ew_ciii_d[valid_d])
        r_d_str = f"{r_d:+.4f} ({valid_d.sum()})"
    if valid_r.sum() >= 30:
        r_r, _ = stats.spearmanr(ew_civ_r[valid_r], ew_ciii_r[valid_r])
        r_r_str = f"{r_r:+.4f} ({valid_r.sum()})"
    
    agree = ""
    if valid_d.sum() >= 30 and valid_r.sum() >= 30:
        agree = "✅" if abs(r_d - r_r) < 0.1 else f"⚠️ Δ={abs(r_d-r_r):.3f}"
    
    print(f"  {z_lo:.1f}-{z_hi:.1f}    {r_d_str:<15} {r_r_str:<17} {agree}")

# FWHM inter-line (should be stable in both)
print(f"\n  CIV-CIII FWHM coupling (should be STABLE in both):")
print(f"  {'z-bin':<10} {'DR16Q r(FW)':<15} {'Rakshit r(FW)':<17}")
print(f"  {'-'*45}")

for z_lo, z_hi in z_bins:
    z_mask = (z_matched > z_lo) & (z_matched < z_hi)
    
    fw_civ_d = dr16_lines['CIV']['fwhm'][m_dr16]
    fw_ciii_d = dr16_lines['CIII_ALL']['fwhm'][m_dr16]
    valid_d = z_mask & np.isfinite(fw_civ_d) & np.isfinite(fw_ciii_d) & (fw_civ_d > 0) & (fw_ciii_d > 0)
    
    fw_civ_r = rak['fwhm_civ'][m_rak]
    fw_ciii_r = rak['fwhm_ciii'][m_rak]
    valid_r = z_mask & np.isfinite(fw_civ_r) & np.isfinite(fw_ciii_r) & (fw_civ_r > 0) & (fw_ciii_r > 0)
    
    r_d_str = r_r_str = "N/A"
    if valid_d.sum() >= 30:
        r_d, _ = stats.spearmanr(fw_civ_d[valid_d], fw_ciii_d[valid_d])
        r_d_str = f"{r_d:+.4f}"
    if valid_r.sum() >= 30:
        r_r, _ = stats.spearmanr(fw_civ_r[valid_r], fw_ciii_r[valid_r])
        r_r_str = f"{r_r:+.4f}"
    
    print(f"  {z_lo:.1f}-{z_hi:.1f}    {r_d_str:<15} {r_r_str:<17}")


# =====================================================================
# VERDICT
# =====================================================================
print(f"\n{'='*80}")
print("VERDICT")
print("=" * 80)
print("""
If BOTH pipelines show:
  - EW inter-line correlation near zero or degrading with z
  - FWHM inter-line correlation significantly positive and stable
Then C1 (pipeline artifact) is DEAD.

The SDSS pipeline and PyQSOFit use completely different:
  - Continuum fitting methods
  - Line profile models  
  - Quality criteria
  - Iron template subtraction

If both independently produce the same EW/FWHM asymmetry,
it's not the measurement. It's the universe.
""")
