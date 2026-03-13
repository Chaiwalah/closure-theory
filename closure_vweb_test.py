#!/usr/bin/env python3
"""
V-WEB STRUCTURAL STATE TEST (Vectorized)

Uses Cosmicflows-4 reconstructed velocity field to classify cosmic web
as void/sheet/filament/knot via velocity shear tensor eigenvalues.
Ray-traces SN Ia and quasar sightlines to compute Fractional Disordered
Path (FDP).

Data: CF4 grouped 64^3 velocity grid (Courtois et al. 2023, A&A 670, L15)
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
import json, os, sys, warnings
warnings.filterwarnings('ignore')

os.makedirs('results_vweb', exist_ok=True)

print("=" * 70)
print("V-WEB STRUCTURAL STATE TEST (Vectorized)")
print("=" * 70)
sys.stdout.flush()

# ============================================================
# 1. LOAD CF4 VELOCITY GRID
# ============================================================
print("\n[1] Loading Cosmicflows-4 velocity grid...")
sys.stdout.flush()

vel_hdu = fits.open('/root/clawd/projects/closure-theory/data/CF4gp_velocity.fits')
vel_raw = vel_hdu[0].data  # (3, 64, 64, 64)
vel = vel_raw * 52.0  # km/s

den_hdu = fits.open('/root/clawd/projects/closure-theory/data/CF4gp_delta.fits')
delta = den_hdu[0].data  # (64, 64, 64)

N = 64
L = 500.0  # half-box Mpc/h
dx = 2 * L / N
H0 = 75.0

print(f"  Grid: {N}^3, box ±{L} Mpc/h, voxel {dx:.1f} Mpc/h")
print(f"  Velocity range: [{vel.min():.0f}, {vel.max():.0f}] km/s")
print(f"  Density range: [{delta.min():.2f}, {delta.max():.2f}]")
sys.stdout.flush()

# ============================================================
# 2. COMPUTE V-WEB (vectorized eigenvalues)
# ============================================================
print("\n[2] Computing velocity shear tensor eigenvalues...")
sys.stdout.flush()

# Partial derivatives via np.gradient
dvdx = np.zeros((3, 3, N, N, N), dtype=np.float32)
for i in range(3):
    for j in range(3):
        dvdx[i][j] = np.gradient(vel[i], dx, axis=j)

# Symmetric shear tensor
sigma = np.zeros((3, 3, N, N, N), dtype=np.float32)
for a in range(3):
    for b in range(3):
        sigma[a][b] = -1.0 / (2.0 * H0) * (dvdx[a][b] + dvdx[b][a])

# Reshape to (N^3, 3, 3) for batch eigenvalue computation
sigma_batch = np.zeros((N**3, 3, 3), dtype=np.float32)
for a in range(3):
    for b in range(3):
        sigma_batch[:, a, b] = sigma[a][b].ravel()

print("  Batch eigenvalue computation...")
sys.stdout.flush()
eigenvalues_flat = np.linalg.eigvalsh(sigma_batch)  # (N^3, 3)
eigenvalues = eigenvalues_flat.reshape(N, N, N, 3)

# Sort descending at each voxel
eigenvalues = np.sort(eigenvalues, axis=-1)[:,:,:,::-1]

for i in range(3):
    print(f"  λ_{i+1}: [{eigenvalues[:,:,:,i].min():.4f}, {eigenvalues[:,:,:,i].max():.4f}]")

del dvdx, sigma, sigma_batch  # free memory
sys.stdout.flush()

# ============================================================
# 3. V-WEB CLASSIFICATION
# ============================================================
print("\n[3] V-web classification...")
sys.stdout.flush()

n_positive = np.sum(eigenvalues > 0, axis=-1)  # (N,N,N)

vweb_class = n_positive.astype(np.int8)
labels = {0: 'Void', 1: 'Sheet', 2: 'Filament', 3: 'Knot'}
for c in range(4):
    count = int(np.sum(vweb_class == c))
    frac = count / N**3 * 100
    print(f"  {labels[c]:10s}: {count:6d} voxels ({frac:5.1f}%)")

frac_disordered = float(np.sum(vweb_class <= 1)) / N**3 * 100
print(f"\n  Disordered (Void+Sheet): {frac_disordered:.1f}%")
print(f"  Ordered (Fil+Knot):      {100-frac_disordered:.1f}%")
sys.stdout.flush()

# ============================================================
# 4. VECTORIZED RAY-TRACING (SN Ia)
# ============================================================
print("\n[4] Ray-tracing Pantheon+ SN Ia sightlines (vectorized)...")
sys.stdout.flush()

pan_data = np.genfromtxt('/root/clawd/projects/closure-theory/data/pantheon_plus.dat',
                         names=True, dtype=None, encoding='utf-8')

ra_sn = pan_data['RA'].astype(float)
dec_sn = pan_data['DEC'].astype(float)
z_sn = pan_data['zHD'].astype(float)
c_sn = pan_data['c'].astype(float)
x1_sn = pan_data['x1'].astype(float)
mu_corr = pan_data['m_b_corr'].astype(float)

valid_sn = np.isfinite(ra_sn) & np.isfinite(dec_sn) & np.isfinite(z_sn) & \
           np.isfinite(c_sn) & np.isfinite(x1_sn) & (z_sn > 0.01)

ra_v = ra_sn[valid_sn]
dec_v = dec_sn[valid_sn]
z_v = z_sn[valid_sn]
c_v = c_sn[valid_sn]
x1_v = x1_sn[valid_sn]
mu_v = mu_corr[valid_sn]
print(f"  Valid SNe: {len(ra_v)}")

c_light = 299792.458

def batch_fdp(ra_arr, dec_arr, z_arr, n_steps=100):
    """Vectorized FDP computation for arrays of sightlines."""
    n_obj = len(ra_arr)
    
    # Convert all positions to supergalactic
    coords = SkyCoord(ra=ra_arr*u.degree, dec=dec_arr*u.degree, frame='icrs')
    sg = coords.supergalactic
    sgl_rad = sg.sgl.rad
    sgb_rad = sg.sgb.rad
    
    # For each step in z, compute grid position for ALL objects at once
    z_max_trace = np.minimum(z_arr, 0.12)  # grid edge
    
    fdp_arr = np.full(n_obj, np.nan)
    den_arr = np.full(n_obj, np.nan)
    
    disordered_counts = np.zeros(n_obj, dtype=int)
    total_counts = np.zeros(n_obj, dtype=int)
    den_sums = np.zeros(n_obj, dtype=float)
    
    for step in range(n_steps):
        frac = (step + 1) / n_steps
        z_step = 0.005 + frac * (z_max_trace - 0.005)
        
        # Skip objects where z_max_trace < 0.005
        active = z_max_trace > 0.005
        if not active.any():
            continue
        
        d_mpc = z_step * c_light / H0
        sgx = d_mpc * np.cos(sgb_rad) * np.cos(sgl_rad)
        sgy = d_mpc * np.cos(sgb_rad) * np.sin(sgl_rad)
        sgz = d_mpc * np.sin(sgb_rad)
        
        # Check bounds
        in_grid = active & (np.abs(sgx) < L) & (np.abs(sgy) < L) & (np.abs(sgz) < L)
        if not in_grid.any():
            continue
        
        # Convert to grid indices
        ix = ((sgx[in_grid] + L) / (2 * L) * N).astype(int)
        iy = ((sgy[in_grid] + L) / (2 * L) * N).astype(int)
        iz = ((sgz[in_grid] + L) / (2 * L) * N).astype(int)
        ix = np.clip(ix, 0, N-1)
        iy = np.clip(iy, 0, N-1)
        iz = np.clip(iz, 0, N-1)
        
        # Look up V-web class (axis order: iz, iy, ix for FITS)
        web_vals = vweb_class[iz, iy, ix]
        den_vals = delta[iz, iy, ix]
        
        # Update counts
        disordered_counts[in_grid] += (web_vals <= 1).astype(int)
        total_counts[in_grid] += 1
        den_sums[in_grid] += den_vals
    
    valid = total_counts >= 10
    fdp_arr[valid] = disordered_counts[valid] / total_counts[valid]
    den_arr[valid] = den_sums[valid] / total_counts[valid]
    
    return fdp_arr, den_arr

print("  Computing FDP...")
sys.stdout.flush()
fdp_sn, den_sn_path = batch_fdp(ra_v, dec_v, z_v)

good = np.isfinite(fdp_sn)
print(f"  SNe with valid FDP: {good.sum()}")
if good.sum() > 0:
    print(f"  FDP range: [{fdp_sn[good].min():.3f}, {fdp_sn[good].max():.3f}]")
    print(f"  FDP mean:  {fdp_sn[good].mean():.3f}")
sys.stdout.flush()

# ============================================================
# 5. CORRELATE FDP WITH SN PARAMETERS
# ============================================================
print("\n[5] Correlating FDP with SN Ia parameters...")
sys.stdout.flush()

fdp_g = fdp_sn[good]
c_g = c_v[good]
x1_g = x1_v[good]
mu_g = mu_v[good]
z_g = z_v[good]
den_g = den_sn_path[good]

rho_c, p_c = stats.spearmanr(fdp_g, c_g)
rho_x1, p_x1 = stats.spearmanr(fdp_g, x1_g)
rho_mu, p_mu = stats.spearmanr(fdp_g, mu_g)
rho_den_c, p_den_c = stats.spearmanr(den_g, c_g)

# Partial correlation: FDP vs color controlling for density
fdp_resid = fdp_g - np.polyval(np.polyfit(den_g, fdp_g, 1), den_g)
c_resid = c_g - np.polyval(np.polyfit(den_g, c_g, 1), den_g)
rho_partial, p_partial = stats.spearmanr(fdp_resid, c_resid)

print(f"  FDP vs Color (c):     ρ = {rho_c:+.4f}, p = {p_c:.2e}")
print(f"  FDP vs Stretch (x1):  ρ = {rho_x1:+.4f}, p = {p_x1:.2e}")
print(f"  FDP vs μ_corr:        ρ = {rho_mu:+.4f}, p = {p_mu:.2e}")
print(f"  Density vs Color:     ρ = {rho_den_c:+.4f}, p = {p_den_c:.2e}")
print(f"  FDP vs Color (dens-controlled): ρ = {rho_partial:+.4f}, p = {p_partial:.2e}")
sys.stdout.flush()

# ============================================================
# 6. QUINTILE ANALYSIS
# ============================================================
print("\n[6] Quintile analysis...")
sys.stdout.flush()

quintiles = np.percentile(fdp_g, [20, 40, 60, 80])
q_labels_list = ['Q1 (ordered)', 'Q2', 'Q3', 'Q4', 'Q5 (disordered)']

print(f"\n  {'Quintile':<20s} {'N':>5s} {'mean(c)':>10s} {'std(c)':>10s} {'mean(x1)':>10s} {'std(x1)':>10s} {'FDP':>8s}")
print("  " + "-" * 75)

q_stds_c = []
q_stds_x1 = []
q_fdps = []

for qi in range(5):
    if qi == 0:
        mask = fdp_g <= quintiles[0]
    elif qi == 4:
        mask = fdp_g > quintiles[3]
    else:
        mask = (fdp_g > quintiles[qi-1]) & (fdp_g <= quintiles[qi])
    
    n_q = mask.sum()
    if n_q < 5:
        continue
    mc = c_g[mask].mean()
    sc = c_g[mask].std()
    mx = x1_g[mask].mean()
    sx = x1_g[mask].std()
    mf = fdp_g[mask].mean()
    
    q_stds_c.append(sc)
    q_stds_x1.append(sx)
    q_fdps.append(mf)
    
    print(f"  {q_labels_list[qi]:<20s} {n_q:5d} {mc:10.4f} {sc:10.4f} {mx:10.4f} {sx:10.4f} {mf:8.4f}")

if len(q_fdps) >= 3:
    rho_q_c, _ = stats.spearmanr(q_fdps, q_stds_c)
    rho_q_x1, _ = stats.spearmanr(q_fdps, q_stds_x1)
    print(f"\n  Quintile FDP vs color variance:   ρ = {rho_q_c:+.3f}")
    print(f"  Quintile FDP vs stretch variance: ρ = {rho_q_x1:+.3f}")
else:
    rho_q_c = rho_q_x1 = float('nan')

# Covariance comparison Q1 vs Q5
q1_mask = fdp_g <= quintiles[0]
q5_mask = fdp_g > quintiles[3]
if q1_mask.sum() > 5 and q5_mask.sum() > 5:
    cov_q1 = np.cov(c_g[q1_mask], x1_g[q1_mask])
    cov_q5 = np.cov(c_g[q5_mask], x1_g[q5_mask])
    print(f"\n  Q1 cov(c,x1): {cov_q1[0,1]:.6f}  |  Q5 cov(c,x1): {cov_q5[0,1]:.6f}")
    print(f"  Q1 var(c):    {cov_q1[0,0]:.6f}  |  Q5 var(c):    {cov_q5[0,0]:.6f}")

sys.stdout.flush()

# ============================================================
# 7. QUASAR SIGHTLINES
# ============================================================
print("\n[7] Quasar sightlines through local V-web volume...")
sys.stdout.flush()

hdu_q = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
dq = hdu_q[1].data
z_q = dq['Z_DR16Q']
mg_ew = dq['MGII_BR'][:, 2]
hb_ew = dq['HBETA_BR'][:, 2]
ra_q = dq['RA']
dec_q = dq['DEC']

trans_mask = (z_q >= 0.75) & (z_q < 1.15) & np.isfinite(mg_ew) & np.isfinite(hb_ew) & \
             (mg_ew != 0) & (hb_ew != 0)

np.random.seed(42)
q_indices = np.where(trans_mask)[0]
n_sample = min(30000, len(q_indices))
q_sample = np.random.choice(q_indices, n_sample, replace=False)
print(f"  Sampling {n_sample} transition-zone quasars")
sys.stdout.flush()

# All quasars are at z>0.75, so their sightlines traverse the ENTIRE grid
# Use z=0.12 as effective trace depth for all
fdp_q, den_q = batch_fdp(ra_q[q_sample], dec_q[q_sample], 
                          np.full(n_sample, 0.12))

good_q = np.isfinite(fdp_q)
print(f"  Quasars with valid FDP: {good_q.sum()}")
sys.stdout.flush()

if good_q.sum() > 100:
    mg_q = mg_ew[q_sample[good_q]]
    hb_q = hb_ew[q_sample[good_q]]
    fdp_qq = fdp_q[good_q]
    
    log_mg = np.log10(np.abs(mg_q) + 1)
    log_hb = np.log10(np.abs(hb_q) + 1)
    
    q_fdp_bins = np.percentile(fdp_qq, [20, 40, 60, 80])
    print(f"\n  {'FDP Quintile':<20s} {'N':>6s} {'ρ(MgII,Hβ)':>12s} {'p-value':>12s}")
    print("  " + "-" * 55)
    
    q_rhos = []
    q_fdp_mids = []
    for qi in range(5):
        if qi == 0:
            mask = fdp_qq <= q_fdp_bins[0]
        elif qi == 4:
            mask = fdp_qq > q_fdp_bins[3]
        else:
            mask = (fdp_qq > q_fdp_bins[qi-1]) & (fdp_qq <= q_fdp_bins[qi])
        
        if mask.sum() > 30:
            r, p = stats.spearmanr(log_mg[mask], log_hb[mask])
            q_rhos.append(r)
            q_fdp_mids.append(fdp_qq[mask].mean())
            print(f"  Q{qi+1} (FDP={fdp_qq[mask].mean():.3f}){'':<5s} {mask.sum():6d} {r:12.4f} {p:12.2e}")
    
    if len(q_rhos) >= 3:
        rho_trend, p_trend = stats.spearmanr(q_fdp_mids, q_rhos)
        print(f"\n  Coupling trend (FDP vs ρ): ρ = {rho_trend:+.3f}, p = {p_trend:.3f}")
    
    # Overall FDP-density independence
    rho_fd_q, p_fd_q = stats.spearmanr(fdp_qq, den_q[good_q])
    print(f"  FDP vs Density (quasars): ρ = {rho_fd_q:+.4f}")
else:
    rho_trend = p_trend = rho_fd_q = float('nan')
    q_rhos = []

sys.stdout.flush()

# ============================================================
# 8. FDP vs DENSITY INDEPENDENCE
# ============================================================
print("\n[8] FDP vs Density independence...")
rho_fd, p_fd = stats.spearmanr(fdp_g, den_g)
print(f"  FDP vs Density (SN Ia): ρ = {rho_fd:+.4f}, p = {p_fd:.2e}")
print(f"  If FDP ≠ density AND FDP predicts degradation → STATE wins over DENSITY")
sys.stdout.flush()

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY — V-WEB STRUCTURAL STATE TEST")
print("=" * 70)

print(f"\nCosmic Web Composition (V-web, λ_th=0):")
for cv in range(4):
    count = int(np.sum(vweb_class == cv))
    print(f"  {labels[cv]:10s}: {count/N**3*100:.1f}%")

print(f"\nSN Ia Results (N={good.sum()}):")
print(f"  FDP vs Color:             ρ = {rho_c:+.4f} (p = {p_c:.2e})")
print(f"  FDP vs Stretch:           ρ = {rho_x1:+.4f} (p = {p_x1:.2e})")
print(f"  FDP vs μ_corr:            ρ = {rho_mu:+.4f} (p = {p_mu:.2e})")
print(f"  Density vs Color:         ρ = {rho_den_c:+.4f} (p = {p_den_c:.2e})")
print(f"  FDP vs Color (dens-ctrl): ρ = {rho_partial:+.4f} (p = {p_partial:.2e})")
if not np.isnan(rho_q_c):
    print(f"  Q-trend (FDP→color var):  ρ = {rho_q_c:+.3f}")
    print(f"  Q-trend (FDP→x1 var):     ρ = {rho_q_x1:+.3f}")

if good_q.sum() > 100 and len(q_rhos) >= 3:
    print(f"\nQuasar Results (N={good_q.sum()}):")
    print(f"  FDP coupling trend:       ρ = {rho_trend:+.3f} (p = {p_trend:.3f})")
    print(f"  FDP vs Density:           ρ = {rho_fd_q:+.4f}")

print(f"\nFDP-Density Independence: ρ = {rho_fd:+.4f}")

print("\n--- INTERPRETATION ---")
if abs(rho_c) > abs(rho_den_c) and abs(rho_c) > 0.05:
    print("  ✅ FDP predicts color BETTER than density → State hypothesis supported")
elif abs(rho_den_c) > abs(rho_c):
    print("  ❌ Density predicts color better than FDP → State hypothesis not supported")
else:
    print("  ~ Neither strongly predictive at this resolution")
    print("    (64^3 grid = 15.6 Mpc/h voxels, may be too coarse)")

# Save
results = {
    'vweb_fractions': {labels[cv]: float(np.sum(vweb_class == cv) / N**3) for cv in range(4)},
    'sn_results': {
        'n_valid': int(good.sum()),
        'fdp_vs_color': {'rho': float(rho_c), 'p': float(p_c)},
        'fdp_vs_stretch': {'rho': float(rho_x1), 'p': float(p_x1)},
        'fdp_vs_mu': {'rho': float(rho_mu), 'p': float(p_mu)},
        'density_vs_color': {'rho': float(rho_den_c), 'p': float(p_den_c)},
        'fdp_vs_color_ctrl': {'rho': float(rho_partial), 'p': float(p_partial)},
    },
    'fdp_density_independence': {'rho': float(rho_fd), 'p': float(p_fd)},
}

with open('results_vweb/vweb_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nResults saved to results_vweb/vweb_results.json")
print("=" * 70)
