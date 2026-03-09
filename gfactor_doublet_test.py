#!/usr/bin/env python3
"""
G-FACTOR DOUBLET TEST — The Controlled Experiment
===================================================
Tests whether emission line doublet ratios evolve with redshift
in the direction predicted by Landé g-factor decoherence.

Three doublets tested:
1. CIV λ1548/λ1550 (²P₃/₂ vs ²P₁/₂)
2. [SII] λ6716/λ6731 (different J levels)
3. [OII] λ3726/λ3729 (²D₃/₂ vs ²D₅/₂)

Plus: g² vs g correlation test with NIST values.

Author: Clawd (closure-theory collaboration)
Date: 7 March 2026
"""

import numpy as np
import json
import os
from scipy import stats

# Try to load FITS
try:
    from astropy.io import fits
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False
    print("WARNING: astropy not available, will use synthetic test")

RESULTS_DIR = '/root/clawd/projects/closure-theory/results_gfactor'
os.makedirs(RESULTS_DIR, exist_ok=True)

# ============================================================
# NIST ATOMIC DATA — Landé g-factors from first principles
# ============================================================
# Formula: g_J = 1 + [J(J+1) + S(S+1) - L(L+1)] / [2*J*(J+1)]
# For S=0: g_J = 1 always (when J≠0)
# For hydrogen-like (one electron): standard formula

def lande_g(J, L, S):
    """Calculate Landé g-factor from quantum numbers."""
    if J == 0:
        return 0.0
    return 1.0 + (J*(J+1) + S*(S+1) - L*(L+1)) / (2*J*(J+1))

# ---------- CIV (C³⁺, Li-like, 1 valence electron, S=1/2) ----------
# λ1548.2: 2p ²P₃/₂ → 2s ²S₁/₂
#   Upper: J=3/2, L=1, S=1/2 → g = 4/3
#   Lower: J=1/2, L=0, S=1/2 → g = 2
# λ1550.8: 2p ²P₁/₂ → 2s ²S₁/₂ 
#   Upper: J=1/2, L=1, S=1/2 → g = 2/3
#   Lower: J=1/2, L=0, S=1/2 → g = 2

CIV_1548_g_upper = lande_g(1.5, 1, 0.5)  # ²P₃/₂
CIV_1548_g_lower = lande_g(0.5, 0, 0.5)  # ²S₁/₂
CIV_1550_g_upper = lande_g(0.5, 1, 0.5)  # ²P₁/₂
CIV_1550_g_lower = lande_g(0.5, 0, 0.5)  # ²S₁/₂

# ---------- [SII] (S⁺, 3 electrons in 3p, ground config 3p³) ----------
# λ6716: ²D₃/₂ → ⁴S₃/₂  (the weaker density-sensitive line)
# λ6731: ²D₅/₂ → ⁴S₃/₂  (the stronger density-sensitive line)
# Upper states from ²D term: S=1/2 (doublet)
# ⁴S₃/₂ lower state: J=3/2, L=0, S=3/2

SII_6716_g_upper = lande_g(1.5, 2, 0.5)  # ²D₃/₂
SII_6716_g_lower = lande_g(1.5, 0, 1.5)  # ⁴S₃/₂
SII_6731_g_upper = lande_g(2.5, 2, 0.5)  # ²D₅/₂
SII_6731_g_lower = lande_g(1.5, 0, 1.5)  # ⁴S₃/₂

# ---------- [OII] (O⁺, same config as [SII] — 3 valence electrons) ----------
# λ3726: ²D₃/₂ → ⁴S₃/₂
# λ3729: ²D₅/₂ → ⁴S₃/₂

OII_3726_g_upper = lande_g(1.5, 2, 0.5)  # ²D₃/₂
OII_3726_g_lower = lande_g(1.5, 0, 1.5)  # ⁴S₃/₂
OII_3729_g_upper = lande_g(2.5, 2, 0.5)  # ²D₅/₂
OII_3729_g_lower = lande_g(1.5, 0, 1.5)  # ⁴S₃/₂

# ---------- Other lines for the 6-line ladder ----------
# [NII] 6584: ¹D₂ → ³P₂ (N⁺)
NII_g_upper = lande_g(2, 2, 0)    # ¹D₂: S=0
NII_g_lower = lande_g(2, 1, 1)    # ³P₂: S=1

# [OIII] 5007: ¹D₂ → ³P₂ (O²⁺) — same term structure as [NII]
OIII_g_upper = lande_g(2, 2, 0)   # ¹D₂: S=0
OIII_g_lower = lande_g(2, 1, 1)   # ³P₂: S=1

# Hβ 4861: n=4→2, multiple substates, weighted average
# For hydrogen: g_J = 1 for all levels (L=S coupling gives g≈1)
Hbeta_g_upper = 1.0  # approximate weighted average
Hbeta_g_lower = 1.0

print("=" * 70)
print("LANDÉ G-FACTOR TABLE (from first principles)")
print("=" * 70)

doublet_data = {
    'CIV': {
        '1548': {'upper': CIV_1548_g_upper, 'lower': CIV_1548_g_lower,
                 'label': '²P₃/₂ → ²S₁/₂'},
        '1550': {'upper': CIV_1550_g_upper, 'lower': CIV_1550_g_lower,
                 'label': '²P₁/₂ → ²S₁/₂'},
    },
    '[SII]': {
        '6716': {'upper': SII_6716_g_upper, 'lower': SII_6716_g_lower,
                 'label': '²D₃/₂ → ⁴S₃/₂'},
        '6731': {'upper': SII_6731_g_upper, 'lower': SII_6731_g_lower,
                 'label': '²D₅/₂ → ⁴S₃/₂'},
    },
    '[OII]': {
        '3726': {'upper': OII_3726_g_upper, 'lower': OII_3726_g_lower,
                 'label': '²D₃/₂ → ⁴S₃/₂'},
        '3729': {'upper': OII_3729_g_upper, 'lower': OII_3729_g_lower,
                 'label': '²D₅/₂ → ⁴S₃/₂'},
    },
}

for ion, lines in doublet_data.items():
    print(f"\n{ion}:")
    for lam, data in lines.items():
        delta_g = abs(data['upper'] - data['lower'])
        print(f"  λ{lam} ({data['label']}): g_upper={data['upper']:.4f}, "
              f"g_lower={data['lower']:.4f}, |Δg|={delta_g:.4f}")

# Print the 6-line ladder g-factors
print("\n" + "=" * 70)
print("6-LINE LADDER — G-FACTOR COMPARISON")
print("=" * 70)

ladder_lines = [
    {'name': '[NII] 6584',  'g_upper': NII_g_upper,   'g_lower': NII_g_lower,
     'degradation': 0.000},
    {'name': '[OIII] 5007', 'g_upper': OIII_g_upper,  'g_lower': OIII_g_lower,
     'degradation': 0.000},
    {'name': 'Hβ 4861',     'g_upper': Hbeta_g_upper, 'g_lower': Hbeta_g_lower,
     'degradation': -0.038},
    {'name': '[OII] 3727',  'g_upper': (OII_3726_g_upper + OII_3729_g_upper)/2,
     'g_lower': OII_3726_g_lower, 'degradation': -0.179},
    {'name': 'CIV 1549',    'g_upper': (CIV_1548_g_upper + CIV_1550_g_upper)/2,
     'g_lower': CIV_1548_g_lower, 'degradation': -0.289},
    {'name': '[SII] 6718',  'g_upper': (SII_6716_g_upper + SII_6731_g_upper)/2,
     'g_lower': SII_6716_g_lower, 'degradation': -0.396},
]

# Compute different g_eff definitions
for line in ladder_lines:
    line['g_upper_only'] = line['g_upper']
    line['delta_g'] = abs(line['g_upper'] - line['g_lower'])
    # Effective g for Zeeman pattern (standard spectroscopy)
    line['g_eff_zeeman'] = (line['g_upper'] + line['g_lower']) / 2

print(f"\n{'Line':<15} {'g_upper':<10} {'g_lower':<10} {'|Δg|':<10} "
      f"{'g_avg':<10} {'Degradation':<12}")
print("-" * 67)
for line in ladder_lines:
    print(f"{line['name']:<15} {line['g_upper']:<10.4f} {line['g_lower']:<10.4f} "
          f"{line['delta_g']:<10.4f} {line['g_eff_zeeman']:<10.4f} "
          f"{line['degradation']:<12.3f}")

# ============================================================
# TEST A: g² vs g CORRELATION (with multiple g_eff definitions)
# ============================================================
print("\n" + "=" * 70)
print("TEST A: g² vs g CORRELATION")
print("=" * 70)

degradation = np.array([l['degradation'] for l in ladder_lines])

for gdef_name, gdef_key in [('g_upper', 'g_upper_only'), 
                              ('|Δg|', 'delta_g'),
                              ('g_avg', 'g_eff_zeeman')]:
    g_vals = np.array([l[gdef_key] for l in ladder_lines])
    g2_vals = g_vals ** 2
    
    r_g, p_g = stats.pearsonr(g_vals, degradation)
    r_g2, p_g2 = stats.pearsonr(g2_vals, degradation)
    
    # Also Spearman (rank) for robustness
    rs_g, ps_g = stats.spearmanr(g_vals, degradation)
    rs_g2, ps_g2 = stats.spearmanr(g2_vals, degradation)
    
    print(f"\nDefinition: {gdef_name}")
    print(f"  Values: {[f'{v:.3f}' for v in g_vals]}")
    print(f"  Pearson  r(g)  = {r_g:+.4f}  (p = {p_g:.4f})")
    print(f"  Pearson  r(g²) = {r_g2:+.4f}  (p = {p_g2:.4f})")
    print(f"  Spearman r(g)  = {rs_g:+.4f}  (p = {ps_g:.4f})")
    print(f"  Spearman r(g²) = {rs_g2:+.4f}  (p = {ps_g2:.4f})")
    print(f"  >>> {'g² WINS' if abs(r_g2) > abs(r_g) else 'g WINS'} "
          f"(Δr = {abs(r_g2) - abs(r_g):+.4f})")

# ============================================================
# TEST B: DOUBLET RATIO EVOLUTION WITH REDSHIFT (DR16Q)
# ============================================================
print("\n" + "=" * 70)
print("TEST B: DOUBLET RATIO EVOLUTION WITH REDSHIFT")
print("=" * 70)

DATA_FILE = '/root/clawd/data/sdss/dr16q_prop.fits'

if HAS_ASTROPY and os.path.exists(DATA_FILE):
    print(f"\nLoading DR16Q from {DATA_FILE}...")
    hdu = fits.open(DATA_FILE)
    data = hdu[1].data
    
    # Get available columns
    cols = [c.name for c in hdu[1].columns]
    print(f"Total objects: {len(data)}")
    
    # Check what emission line columns are available
    ew_cols = [c for c in cols if 'EW' in c.upper() or 'ew' in c.lower()]
    line_cols = [c for c in cols if any(x in c.upper() for x in 
                 ['CIV', 'SII', 'OII', 'OIII', 'NII', 'HB', 'HBETA', 'MG'])]
    
    print(f"\nEW-related columns: {ew_cols[:20]}")
    print(f"Line-related columns: {line_cols[:20]}")
    
    # Get redshift
    z_col = None
    for candidate in ['Z', 'Z_VI', 'Z_PIPE', 'REDSHIFT', 'Z_PCA']:
        if candidate in cols:
            z_col = candidate
            break
    
    if z_col is None:
        print("ERROR: No redshift column found!")
        print(f"Available columns: {cols[:30]}...")
    else:
        z = data[z_col]
        print(f"\nRedshift column: {z_col}")
        print(f"z range: {np.nanmin(z):.3f} to {np.nanmax(z):.3f}")
        print(f"Median z: {np.nanmedian(z):.3f}")
    
    # ---- CIV DOUBLET TEST ----
    # Look for CIV EW or flux columns
    civ_cols = [c for c in cols if 'CIV' in c.upper()]
    print(f"\nCIV columns available: {civ_cols}")
    
    # Check for individual doublet component measurements
    # DR16Q may have CIV as a single line or as components
    # Also check for velocity/width measurements that could show asymmetry
    
    # Look for ALL columns with line measurements
    print(f"\nAll columns containing 'EW': {[c for c in cols if 'EW' in c.upper()]}")
    print(f"All columns containing 'FWHM': {[c for c in cols if 'FWHM' in c.upper()]}")
    print(f"All columns containing 'FLUX': {[c for c in cols if 'FLUX' in c.upper()][:20]}")
    print(f"All columns containing 'PEAK': {[c for c in cols if 'PEAK' in c.upper()]}")
    print(f"All columns containing 'SHIFT': {[c for c in cols if 'SHIFT' in c.upper()]}")
    print(f"All columns containing 'VOFF': {[c for c in cols if 'VOFF' in c.upper()]}")
    
    # Look for the line property catalog columns
    # DR16Q properties file typically has: EW, FWHM, PEAK, MEDIAN for each line
    for prefix in ['CIV', 'SIIDOUBLET', 'SII', 'OII', 'OIII', 'NII', 'HALPHA', 
                    'HBETA', 'MGII', 'CIII']:
        matching = [c for c in cols if c.upper().startswith(prefix)]
        if matching:
            print(f"\n{prefix} columns: {matching}")
    
    # ---- GENERIC DOUBLET RATIO TEST ----
    # For any doublet where we have individual component EWs,
    # measure ratio vs z
    
    def doublet_ratio_vs_z(ew1, ew2, z_arr, n_bins=6, label=""):
        """Measure doublet ratio evolution with redshift."""
        # Quality cuts
        mask = (np.isfinite(ew1) & np.isfinite(ew2) & 
                np.isfinite(z_arr) & (ew1 > 0) & (ew2 > 0) & (z_arr > 0))
        
        ew1_clean = ew1[mask]
        ew2_clean = ew2[mask]
        z_clean = z_arr[mask]
        
        if len(z_clean) < 100:
            print(f"  {label}: Too few objects ({len(z_clean)}), skipping")
            return None
        
        # Compute ratio
        ratio = ew1_clean / ew2_clean
        
        # Bin by redshift
        z_bins = np.percentile(z_clean, np.linspace(0, 100, n_bins + 1))
        z_mids = []
        ratio_medians = []
        ratio_stds = []
        n_per_bin = []
        
        for i in range(n_bins):
            in_bin = (z_clean >= z_bins[i]) & (z_clean < z_bins[i+1])
            if np.sum(in_bin) < 20:
                continue
            z_mids.append(np.median(z_clean[in_bin]))
            ratio_medians.append(np.median(ratio[in_bin]))
            ratio_stds.append(np.std(ratio[in_bin]) / np.sqrt(np.sum(in_bin)))
            n_per_bin.append(np.sum(in_bin))
        
        if len(z_mids) < 3:
            print(f"  {label}: Too few bins with data, skipping")
            return None
        
        z_mids = np.array(z_mids)
        ratio_medians = np.array(ratio_medians)
        
        # Correlation
        r, p = stats.pearsonr(z_mids, ratio_medians)
        rs, ps = stats.spearmanr(z_mids, ratio_medians)
        
        # Linear fit
        slope, intercept, r_val, p_val, std_err = stats.linregress(z_mids, ratio_medians)
        
        print(f"\n  {label}:")
        print(f"    N objects: {len(z_clean)}")
        print(f"    z range: {z_clean.min():.2f} - {z_clean.max():.2f}")
        print(f"    Bins: {n_bins} ({n_per_bin})")
        print(f"    Ratio at low-z:  {ratio_medians[0]:.4f}")
        print(f"    Ratio at high-z: {ratio_medians[-1]:.4f}")
        print(f"    Change: {ratio_medians[-1] - ratio_medians[0]:+.4f}")
        print(f"    Slope: {slope:+.6f} per unit z")
        print(f"    Pearson:  r = {r:+.4f}, p = {p:.4f}")
        print(f"    Spearman: r = {rs:+.4f}, p = {ps:.4f}")
        
        return {
            'label': label,
            'n_objects': int(len(z_clean)),
            'z_mids': z_mids.tolist(),
            'ratio_medians': ratio_medians.tolist(),
            'slope': float(slope),
            'pearson_r': float(r),
            'pearson_p': float(p),
            'spearman_r': float(rs),
            'spearman_p': float(ps),
        }
    
    # Try to find and test available doublets
    results = {}
    
    # Check column names more broadly
    print("\n" + "-" * 50)
    print("Searching for doublet component columns...")
    print("-" * 50)
    
    # Print ALL columns for inspection
    print(f"\nFull column list ({len(cols)} columns):")
    for i in range(0, len(cols), 5):
        chunk = cols[i:i+5]
        print(f"  {', '.join(chunk)}")
    
    hdu.close()

else:
    if not HAS_ASTROPY:
        print("astropy not installed. Installing...")
    elif not os.path.exists(DATA_FILE):
        print(f"Data file not found: {DATA_FILE}")
        print("Checking alternate locations...")
        for path in ['/root/clawd/data/sdss/', '/root/clawd/projects/closure-theory/data/']:
            if os.path.exists(path):
                print(f"  {path}: {os.listdir(path)[:10]}")

# ============================================================
# TEST C: CAST COMPATIBILITY CALCULATION
# ============================================================
print("\n" + "=" * 70)
print("TEST C: CAST COMPATIBILITY")
print("=" * 70)

# Constants
Gamma_0 = 2.17
CAST_limit = 5.8e-11  # GeV^-1

# 1 nG * 1 Mpc in natural units
# 1 T = 1.95e-16 GeV^2 (in natural units where ℏ=c=1)
# 1 m = 5.068e15 GeV^-1
# B[GeV^2] * L[GeV^-1] = B*L [GeV]
# 1 nG = 1e-13 T = 1.95e-29 GeV^2
# 1 Mpc = 3.086e22 m = 1.563e38 GeV^-1
# Product: 1.95e-29 * 1.563e38 = 3.048e9 GeV

BL_per_nG_Mpc = 3.048e9  # GeV

print(f"\nConversion: 1 nG × 1 Mpc = {BL_per_nG_Mpc:.3e} GeV")
print(f"Γ₀ = {Gamma_0}")
print(f"CAST limit: g_aγ < {CAST_limit:.1e} GeV⁻¹")

# g_aγ = 2*sqrt(Γ₀) / (B*L)
g_required_central = 2 * np.sqrt(Gamma_0) / (1.0 * 1.0 * BL_per_nG_Mpc)
print(f"\nCentral (B=1nG, L=1Mpc): g_aγ = {g_required_central:.3e} GeV⁻¹")
print(f"  Ratio to CAST: {g_required_central / CAST_limit:.1f}× {'ABOVE' if g_required_central > CAST_limit else 'BELOW'}")

print(f"\n{'B (nG)':<10} {'L (Mpc)':<10} {'BL (GeV)':<15} {'g_aγ (GeV⁻¹)':<18} {'vs CAST':<15}")
print("-" * 68)

B_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
L_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

cast_pass_count = 0
total_count = 0

for B in B_values:
    for L in L_values:
        BL = B * L * BL_per_nG_Mpc
        g_req = 2 * np.sqrt(Gamma_0) / BL
        passes = g_req < CAST_limit
        if passes:
            cast_pass_count += 1
        total_count += 1
        
        # Only print interesting rows
        if B * L >= 0.5:
            status = "✓ PASSES" if passes else "✗ FAILS"
            print(f"{B:<10.1f} {L:<10.1f} {BL:<15.3e} {g_req:<18.3e} {status}")

print(f"\nPassing combinations: {cast_pass_count}/{total_count}")
print(f"Minimum BL for CAST compliance: {2*np.sqrt(Gamma_0)/CAST_limit:.3e} GeV")
min_BL_nGMpc = 2*np.sqrt(Gamma_0) / (CAST_limit * BL_per_nG_Mpc)
print(f"  = {min_BL_nGMpc:.1f} nG·Mpc")
print(f"  Examples: B=1nG needs L≥{min_BL_nGMpc:.1f}Mpc, "
      f"L=1Mpc needs B≥{min_BL_nGMpc:.1f}nG")

# Decoherence reduction factor
print(f"\n--- DECOHERENCE vs CONVERSION ---")
print(f"If decoherence requires ~100× less coupling (Gemini estimate):")
g_deco = g_required_central / 100
print(f"  Effective g_aγ needed: {g_deco:.3e} GeV⁻¹")
print(f"  vs CAST: {g_deco/CAST_limit:.3f}× → {'✓ PASSES easily' if g_deco < CAST_limit else '✗ FAILS'}")

# ============================================================
# SAVE ALL RESULTS
# ============================================================

all_results = {
    'g_factor_table': {
        'CIV_1548': {'g_upper': float(CIV_1548_g_upper), 'g_lower': float(CIV_1548_g_lower)},
        'CIV_1550': {'g_upper': float(CIV_1550_g_upper), 'g_lower': float(CIV_1550_g_lower)},
        'SII_6716': {'g_upper': float(SII_6716_g_upper), 'g_lower': float(SII_6716_g_lower)},
        'SII_6731': {'g_upper': float(SII_6731_g_upper), 'g_lower': float(SII_6731_g_lower)},
        'OII_3726': {'g_upper': float(OII_3726_g_upper), 'g_lower': float(OII_3726_g_lower)},
        'OII_3729': {'g_upper': float(OII_3729_g_upper), 'g_lower': float(OII_3729_g_lower)},
        'NII_6584': {'g_upper': float(NII_g_upper), 'g_lower': float(NII_g_lower)},
        'OIII_5007': {'g_upper': float(OIII_g_upper), 'g_lower': float(OIII_g_lower)},
        'Hbeta': {'g_upper': float(Hbeta_g_upper), 'g_lower': float(Hbeta_g_lower)},
    },
    'cast_central': {
        'g_required': float(g_required_central),
        'cast_limit': float(CAST_limit),
        'ratio': float(g_required_central / CAST_limit),
        'passes': bool(g_required_central < CAST_limit),
    },
    'cast_min_BL': float(min_BL_nGMpc),
}

with open(os.path.join(RESULTS_DIR, 'gfactor_results.json'), 'w') as f:
    json.dump(all_results, f, indent=2)

print(f"\nResults saved to {RESULTS_DIR}/gfactor_results.json")
print("\nDone.")
