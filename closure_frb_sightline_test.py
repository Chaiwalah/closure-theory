#!/usr/bin/env python3
"""
FRB SIGHTLINE TEST — Grok's proposal
Cross-match FRB scattering excess with SN Ia breaker map.
If turbulent scattering correlates with breaker patches, we have a new observable class.
"""
import numpy as np
from scipy import stats
from scipy.spatial import cKDTree
import healpy as hp
import json, os, sys, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_frb_sightline', exist_ok=True)

print("="*60)
print("FRB SIGHTLINE TEST")
print("="*60)
sys.stdout.flush()

# 1. Load CHIME FRB catalog
print("\n[1] Loading CHIME/FRB Catalog 1...")
sys.stdout.flush()

# Parse TSV manually (has multi-line header)
lines = []
with open('data/chime_full.tsv') as f:
    for line in f:
        if line.startswith('#') or line.strip() == '': continue
        lines.append(line.strip())

# First 3 lines are header/units/separator
header = lines[0].split('\t')
data_lines = lines[3:]  # skip header, units, separator

# Parse into arrays
frb_data = {}
for col in header:
    frb_data[col] = []

for line in data_lines:
    fields = line.split('\t')
    for i, col in enumerate(header):
        if i < len(fields):
            frb_data[col].append(fields[i].strip())
        else:
            frb_data[col].append('')

# Extract key columns
ra_frb = np.array([float(x) if x and x not in ('', ' ') else np.nan for x in frb_data['RAJ2000']])
dec_frb = np.array([float(x) if x and x not in ('', ' ') else np.nan for x in frb_data['DEJ2000']])
dm_frb = np.array([float(x) if x and x not in ('', ' ') else np.nan for x in frb_data['DM']])
dm_ne2001 = np.array([float(x) if x and x not in ('', ' ') else np.nan for x in frb_data['DMeNE2001']])
scat_frb = np.array([float(x) if x and x not in ('', ' ') else np.nan for x in frb_data['Scat']])
sp_idx = np.array([float(x) if x and x not in ('', ' ') else np.nan for x in frb_data['SpInd']])
sp_run = np.array([float(x) if x and x not in ('', ' ') else np.nan for x in frb_data['spRun']])
glat = np.array([float(x) if x and x not in ('', ' ') else np.nan for x in frb_data['GLAT']])
fluence = np.array([float(x) if x and x not in ('', ' ') else np.nan for x in frb_data['Fluence']])

# DM excess = DM_observed - DM_galactic (NE2001 model)
dm_excess = dm_frb - dm_ne2001

# Valid FRBs: have position, DM, scattering, galactic cut
valid = np.isfinite(ra_frb) & np.isfinite(dec_frb) & np.isfinite(dm_frb) & \
        np.isfinite(scat_frb) & np.isfinite(dm_ne2001) & (np.abs(glat) > 10)

ra_v = ra_frb[valid]; dec_v = dec_frb[valid]
dm_v = dm_frb[valid]; dme_v = dm_excess[valid]
scat_v = scat_frb[valid]; sp_v = sp_idx[valid]
sp_run_v = sp_run[valid]; flu_v = fluence[valid]

print(f"  Total FRBs: {len(ra_frb)}")
print(f"  Valid (pos+DM+scat+|b|>10°): {valid.sum()}")
print(f"  DM range: [{dm_v.min():.1f}, {dm_v.max():.1f}] pc/cm³")
print(f"  DM excess range: [{dme_v.min():.1f}, {dme_v.max():.1f}] pc/cm³")
print(f"  Scattering range: [{scat_v.min():.6f}, {scat_v.max():.4f}] s")

# Scattering excess: residual after DM-scattering relation
# τ ∝ DM^α typically. Compute residual.
log_scat = np.log10(scat_v + 1e-10)
log_dm = np.log10(dme_v.clip(1))
finite_both = np.isfinite(log_scat) & np.isfinite(log_dm) & (dme_v > 10)
if finite_both.sum() > 20:
    slope, intercept, r, p, se = stats.linregress(log_dm[finite_both], log_scat[finite_both])
    scat_excess = log_scat - (slope * log_dm + intercept)
    print(f"\n  DM-scattering relation: log(τ) = {slope:.2f}*log(DM_ex) + {intercept:.2f}, r={r:.3f}")
    print(f"  Scattering excess = residual from this fit")
else:
    scat_excess = log_scat
    print("  Not enough for DM-scat fit, using raw scattering")

sys.stdout.flush()

# 2. Load SN breaker map
print("\n[2] Building SN Ia breaker HEALPix map...")
sys.stdout.flush()

d = np.genfromtxt('data/pantheon_plus.dat', names=True, dtype=None, encoding='utf-8')
ra_sn = d['RA'].astype(float); dec_sn = d['DEC'].astype(float)
z_sn = d['zHD'].astype(float); c_sn = d['c'].astype(float)
val_sn = np.isfinite(ra_sn)&np.isfinite(dec_sn)&np.isfinite(z_sn)&np.isfinite(c_sn)&(z_sn>0.01)
ra_sn,dec_sn,z_sn,c_sn = ra_sn[val_sn],dec_sn[val_sn],z_sn[val_sn],c_sn[val_sn]

is_breaker = np.zeros(len(ra_sn), dtype=bool)
for i in range(len(ra_sn)):
    local = np.abs(z_sn - z_sn[i]) < 0.05
    if local.sum() < 10: continue
    mu_c, sig_c = c_sn[local].mean(), c_sn[local].std()
    if sig_c > 0 and np.abs(c_sn[i] - mu_c) > sig_c:
        is_breaker[i] = True

NSIDE = 8
theta_sn = np.radians(90 - dec_sn); phi_sn = np.radians(ra_sn)
pix_sn = hp.ang2pix(NSIDE, theta_sn, phi_sn)

count_map = np.zeros(hp.nside2npix(NSIDE))
break_map = np.zeros(hp.nside2npix(NSIDE))
for i in range(len(ra_sn)):
    count_map[pix_sn[i]] += 1
    break_map[pix_sn[i]] += is_breaker[i]

frac_map = np.full(hp.nside2npix(NSIDE), np.nan)
good = count_map >= 3
frac_map[good] = break_map[good] / count_map[good]
print(f"  Breaker map: {good.sum()} populated pixels")

# 3. Cross-match FRBs with breaker map
print("\n[3] Cross-matching FRB positions with breaker map...")
sys.stdout.flush()

theta_frb = np.radians(90 - dec_v); phi_frb = np.radians(ra_v)
pix_frb = hp.ang2pix(NSIDE, theta_frb, phi_frb)

frb_breaker_frac = frac_map[pix_frb]
has_sn_data = np.isfinite(frb_breaker_frac) & finite_both[:valid.sum()] if len(finite_both) >= valid.sum() else np.isfinite(frb_breaker_frac)

# Limit to FRBs that land on populated SN pixels
n_matched = has_sn_data.sum()
print(f"  FRBs on populated SN pixels: {n_matched}/{valid.sum()}")

if n_matched > 20:
    # Split FRBs by breaker fraction
    bf = frb_breaker_frac[has_sn_data]
    se = scat_excess[has_sn_data] if len(scat_excess) >= valid.sum() else scat_excess[:valid.sum()][has_sn_data]
    
    # High vs low breaker fraction
    median_bf = np.nanmedian(bf)
    high_bf = bf > median_bf
    low_bf = bf <= median_bf
    
    if high_bf.sum() > 5 and low_bf.sum() > 5:
        mean_scat_high = se[high_bf].mean()
        mean_scat_low = se[low_bf].mean()
        u_stat, u_p = stats.mannwhitneyu(se[high_bf], se[low_bf], alternative='two-sided')
        
        print(f"\n  High breaker-fraction sightlines ({high_bf.sum()}):")
        print(f"    Mean scattering excess: {mean_scat_high:.4f}")
        print(f"  Low breaker-fraction sightlines ({low_bf.sum()}):")
        print(f"    Mean scattering excess: {mean_scat_low:.4f}")
        print(f"  Mann-Whitney U p = {u_p:.4e}")
        print(f"  Theory: HIGH breaker → MORE scattering excess (disordered medium)")
        
        if mean_scat_high > mean_scat_low:
            print(f"  ✅ Correct direction")
        else:
            print(f"  ❌ Wrong direction")
    
    # Correlation
    rho_bf_scat, p_bf_scat = stats.spearmanr(bf, se)
    print(f"\n  Breaker fraction vs scattering excess: ρ = {rho_bf_scat:+.4f}, p = {p_bf_scat:.4e}")

# 4. FRB self-consistency: DM vs scattering decorrelation
print("\n[4] FRB internal consistency checks...")
sys.stdout.flush()

# Does scattering become unpredictable at high DM (like coupling breaks at high z)?
dm_bins = np.percentile(dme_v[finite_both[:valid.sum()] if len(finite_both) >= valid.sum() else np.ones(valid.sum(), dtype=bool)], [25, 50, 75])
print(f"  DM excess quartile boundaries: {dm_bins}")

fb = finite_both[:valid.sum()] if len(finite_both) >= valid.sum() else np.ones(valid.sum(), dtype=bool)
dme_fb = dme_v[fb]; scat_fb = scat_v[fb]; se_fb = scat_excess[fb] if len(scat_excess) >= valid.sum() else scat_excess[:valid.sum()][fb]

bin_edges = [0] + list(dm_bins) + [99999]
print(f"\n  {'DM_ex bin':<20s} {'N':>5s} {'std(τ_excess)':>15s} {'ρ(DM,τ)':>10s}")
print("  " + "-"*55)
variances = []
for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
    mask = (dme_fb >= lo) & (dme_fb < hi)
    if mask.sum() > 10:
        std_excess = se_fb[mask].std()
        rho_local, _ = stats.spearmanr(dme_fb[mask], scat_fb[mask])
        variances.append(std_excess)
        print(f"  {lo:.0f}-{hi:.0f}{'':>8s} {mask.sum():5d} {std_excess:15.4f} {rho_local:10.4f}")

if len(variances) >= 3:
    rho_var, _ = stats.spearmanr(range(len(variances)), variances)
    print(f"\n  Variance trend across DM bins: ρ = {rho_var:+.3f}")
    print(f"  Theory: variance should INCREASE with DM (coupling degrades)")

sys.stdout.flush()

# 5. Spectral index as path variable
print("\n[5] Spectral index analysis...")
sys.stdout.flush()

sp_valid = np.isfinite(sp_v) & np.isfinite(frb_breaker_frac)
if sp_valid.sum() > 20:
    rho_sp, p_sp = stats.spearmanr(frb_breaker_frac[sp_valid], sp_v[sp_valid])
    print(f"  Spectral index vs breaker fraction: ρ = {rho_sp:+.4f}, p = {p_sp:.4e}")
    print(f"  Theory: breaker sightlines may show steeper/different spectral index")

# SUMMARY
print("\n" + "="*60)
print("SUMMARY — FRB SIGHTLINE TEST")
print("="*60)
print(f"  FRBs matched to SN breaker map: {n_matched}")
if n_matched > 20:
    print(f"  Breaker frac vs scat excess: ρ = {rho_bf_scat:+.4f}, p = {p_bf_scat:.4e}")
print(f"  DM-scat relation: log(τ) = {slope:.2f}*log(DM) + {intercept:.2f}")

results = {
    'n_frbs_valid': int(valid.sum()),
    'n_matched_to_sn_map': int(n_matched),
    'dm_scat_slope': round(slope, 4) if 'slope' in dir() else None,
    'dm_scat_r': round(r, 4) if 'r' in dir() else None,
}
if n_matched > 20:
    results['bf_vs_scat_rho'] = round(rho_bf_scat, 4)
    results['bf_vs_scat_p'] = round(p_bf_scat, 6)

with open('results_frb_sightline/frb_sightline_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nSaved to results_frb_sightline/frb_sightline_results.json")
print("="*60)
