#!/usr/bin/env python3
"""Closure Theory Data Extraction Pipeline.

Cross-matches Pantheon+ SNe Ia with Planck CMB lensing (κ) and tSZ (y) maps,
computes Hubble residuals, and produces diagnostic plots.
"""

import os, sys, json, subprocess, tarfile, glob, gc
import numpy as np
import pandas as pd
from pathlib import Path

# Paths
BASE = Path("/root/clawd/projects/closure-theory")
DATA_DIR = BASE / "data"
MAP_DIR = BASE / "maps"
PLOT_DIR = BASE / "plots"
OUT_DIR = BASE / "output"
for d in [MAP_DIR, PLOT_DIR, OUT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── Step 0: Load SN data ──────────────────────────────────────────────────
print("Loading Pantheon+ data...")
df = pd.read_csv(DATA_DIR / "pantheon_plus.dat", delim_whitespace=True, comment='#')
print(f"  Loaded {len(df)} rows, {df['CID'].nunique()} unique SNe")

# Average duplicate observations per unique SN
agg_cols = {
    'zHD': 'first', 'zHDERR': 'first', 'MU_SH0ES': 'mean', 'MU_SH0ES_ERR_DIAG': 'first',
    'RA': 'first', 'DEC': 'first', 'MWEBV': 'first', 'HOST_LOGMASS': 'first',
    'IDSURVEY': 'first', 'zCMB': 'first',
}
sn = df.groupby('CID', sort=False).agg(agg_cols).reset_index()
sn.rename(columns={'MU_SH0ES': 'MU', 'MU_SH0ES_ERR_DIAG': 'MUERR'}, inplace=True)
print(f"  Unique SNe: {len(sn)}")

# ── Step 1: Convert RA/DEC to Galactic ────────────────────────────────────
print("Converting to Galactic coordinates...")
from astropy.coordinates import SkyCoord
import astropy.units as u
coords = SkyCoord(ra=sn['RA'].values*u.deg, dec=sn['DEC'].values*u.deg, frame='icrs')
gal = coords.galactic
sn['GLON'] = gal.l.deg
sn['GLAT'] = gal.b.deg

# ── Step 2: Download and process maps ─────────────────────────────────────
import healpy as hp

NSIDE_WORK = 512  # Safe for 2GB RAM (~12MB per map at float32)

def download_file(url, dest):
    """Download with wget."""
    print(f"  Downloading {url}...")
    r = subprocess.run(['wget', '-q', '--no-check-certificate', '-O', str(dest), url],
                       capture_output=True, timeout=600)
    if r.returncode != 0:
        print(f"  FAILED: {r.stderr.decode()[:200]}")
        return False
    print(f"  Downloaded {dest.stat().st_size / 1e6:.1f} MB")
    return True

def extract_map_from_tgz(tgz_path, map_dir, keywords):
    """Extract first FITS file matching keywords from tarball."""
    with tarfile.open(tgz_path, 'r:gz') as tar:
        for member in tar.getmembers():
            name_lower = member.name.lower()
            if any(k in name_lower for k in keywords) and name_lower.endswith('.fits'):
                print(f"  Extracting {member.name}...")
                tar.extract(member, path=map_dir)
                return map_dir / member.name
        # If no keyword match, list what's inside
        print("  Available files in tarball:")
        for m in tar.getmembers():
            if m.name.endswith('.fits'):
                print(f"    {m.name} ({m.size/1e6:.1f}MB)")
    return None

def load_and_downgrade(fits_path, field=0):
    """Load a HEALPix map and downgrade to NSIDE_WORK."""
    print(f"  Loading {fits_path}...")
    try:
        m = hp.read_map(str(fits_path), field=field, dtype=np.float64)
    except Exception as e:
        print(f"  Error reading field {field}: {e}")
        # Try reading all and taking first
        m = hp.read_map(str(fits_path), dtype=np.float64)
    nside_orig = hp.npix2nside(len(m))
    print(f"  Original NSIDE={nside_orig}, downgrading to {NSIDE_WORK}")
    if nside_orig > NSIDE_WORK:
        m = hp.ud_grade(m, NSIDE_WORK)
    elif nside_orig < NSIDE_WORK:
        NSIDE_WORK_actual = nside_orig  # use native
        m = m  # keep as is
    return m

def try_download_kappa_alm(map_dir):
    """Download lensing tgz, convert alm->map."""
    tgz = map_dir / "lensing.tgz"
    url = "https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/lensing/COM_Lensing_4096_R3.00.tgz"
    if not download_file(url, tgz):
        return None, None

    # Extract - look for klm or convergence
    with tarfile.open(tgz, 'r:gz') as tar:
        fits_files = [(m.name, m.size) for m in tar.getmembers() if m.name.endswith('.fits')]
        print(f"  Found {len(fits_files)} FITS files in tarball:")
        for name, size in fits_files:
            print(f"    {name} ({size/1e6:.1f}MB)")

        # Look for the convergence/kappa map or klm
        target = None
        for name, size in fits_files:
            nl = name.lower()
            if 'klm' in nl or 'kappa' in nl or 'convergence' in nl:
                target = name
                break
        if target is None:
            # Take the map file
            for name, size in fits_files:
                nl = name.lower()
                if 'map' in nl or 'dat' in nl:
                    target = name
                    break
        if target is None and fits_files:
            target = fits_files[0][0]

        if target:
            print(f"  Extracting {target}...")
            tar.extract(tar.getmember(target), path=map_dir)
            # Also extract mean-field if exists
            mf_name = target.replace('dat_klm', 'mf_klm')
            try:
                tar.extract(tar.getmember(mf_name), path=map_dir)
                print(f"  Also extracted {mf_name}")
            except KeyError:
                pass

    # Clean up tarball
    tgz.unlink()
    gc.collect()

    if target is None:
        return None, None

    extracted = map_dir / target
    nl = target.lower()

    if 'klm' in nl or 'alm' in nl:
        # It's alms of the lensing potential φ - convert to convergence κ
        # κ_lm = l(l+1)/2 * φ_lm
        # Also need to subtract mean-field: dat - mf
        print("  Converting alm to convergence map...")
        alm_dat = hp.read_alm(str(extracted))
        lmax = hp.Alm.getlmax(len(alm_dat))
        
        # Try to load mean-field
        mf_path = str(extracted).replace('dat_klm', 'mf_klm')
        if os.path.exists(mf_path):
            print("  Subtracting mean-field...")
            alm_mf = hp.read_alm(mf_path)
            alm_dat = alm_dat - alm_mf
            del alm_mf
            os.unlink(mf_path)
        
        # dat_klm stores lensing potential φ_lm (mean-field subtracted)
        # Convert to convergence κ_lm = l(l+1)/2 * φ_lm with Wiener filter
        # to suppress reconstruction noise (SNR<<1 per pixel)
        print(f"  Building Wiener-filtered κ map (lmax={lmax})...")
        ls = np.arange(lmax + 1, dtype=np.float64)
        
        # Measure total C_l^φφ (signal + noise)
        cl_total = hp.alm2cl(alm_dat)
        
        # Build signal model: smooth the measured C_l to get signal estimate
        # Use binned average as signal proxy (noise is roughly flat)
        # Better: use Planck 2018 best-fit C_l^φφ approximation
        cl_signal = np.zeros(lmax + 1)
        for l in range(2, lmax + 1):
            # Planck 2018 fiducial C_phi(l) approximation
            cl_signal[l] = 2.1e-7 * (l/60.0)**(-0.3) / (1 + (l/60.0)**3.8)
        
        # Wiener filter in phi space
        wiener = np.where(cl_total > 0, np.minimum(1.0, cl_signal / cl_total), 0.0)
        
        # Combined: Wiener * phi-to-kappa
        fl_kappa = ls * (ls + 1) / 2.0
        fl_kappa[0] = 0
        combined = fl_kappa * wiener
        
        hp.almxfl(alm_dat, combined, inplace=True)
        kappa = hp.alm2map(alm_dat, NSIDE_WORK)
        # Normalize to zero mean, unit std for correlation analysis
        kappa_std = kappa.std()
        kappa = (kappa - kappa.mean()) / kappa_std if kappa_std > 0 else kappa
        print(f"  κ (Wiener-filtered, normalized): std={kappa.std():.4f}, range=[{kappa.min():.4f}, {kappa.max():.4f}]")
        del alm_dat; gc.collect()
        extracted.unlink()
        print(f"  κ range: [{kappa.min():.6f}, {kappa.max():.6f}]")
        return kappa, None
    else:
        # It's a map - try to load
        try:
            # Try field 0 (kappa) and field 1 (mask)
            kappa = load_and_downgrade(extracted, field=0)
            try:
                mask = load_and_downgrade(extracted, field=1)
            except:
                mask = None
            extracted.unlink()
            return kappa, mask
        except Exception as e:
            print(f"  Error: {e}")
            extracted.unlink()
            return None, None

def try_download_ysz(map_dir):
    """Download tSZ y-map."""
    tgz = map_dir / "ysz.tgz"
    url = "https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_YSZ_R2.01.fits.tgz"
    if not download_file(url, tgz):
        return None, None

    with tarfile.open(tgz, 'r:gz') as tar:
        fits_files = [(m.name, m.size) for m in tar.getmembers() if m.name.endswith('.fits')]
        print(f"  Found FITS: {[f[0] for f in fits_files]}")

        # Look for the y-map (MILCA or NILC)
        target = None
        for name, size in fits_files:
            nl = name.lower()
            if 'milca' in nl or 'nilc' in nl or 'ysz' in nl or 'compton' in nl:
                target = name
                break
        if target is None and fits_files:
            target = fits_files[0][0]

        if target:
            print(f"  Extracting {target}...")
            tar.extract(tar.getmember(target), path=map_dir)

    tgz.unlink()
    gc.collect()

    if target is None:
        return None, None

    extracted = map_dir / target
    try:
        ymap = load_and_downgrade(extracted, field=0)
        try:
            ymask = load_and_downgrade(extracted, field=1)
        except:
            ymask = None
        extracted.unlink()
        return ymap, ymask
    except Exception as e:
        print(f"  Error loading y-map: {e}")
        extracted.unlink()
        return None, None

# Download kappa (skip if already extracted)
print("\n=== Downloading Planck Lensing (κ) map ===")
cached_klm = MAP_DIR / "COM_Lensing_4096_R3.00" / "MV" / "dat_klm.fits"
if cached_klm.exists():
    print("  Using cached alm files...")
    import healpy as hp_tmp
    alm_dat = hp.read_alm(str(cached_klm))
    mf_path = str(cached_klm).replace('dat_klm', 'mf_klm')
    if os.path.exists(mf_path):
        alm_mf = hp.read_alm(mf_path)
        alm_dat = alm_dat - alm_mf
        del alm_mf
    lmax = hp.Alm.getlmax(len(alm_dat))
    
    # Wiener-filtered kappa
    ls = np.arange(lmax + 1, dtype=np.float64)
    cl_total = hp.alm2cl(alm_dat)
    cl_signal = np.zeros(lmax + 1)
    for l in range(2, lmax + 1):
        cl_signal[l] = 2.1e-7 * (l/60.0)**(-0.3) / (1 + (l/60.0)**3.8)
    wiener = np.where(cl_total > 0, np.minimum(1.0, cl_signal / cl_total), 0.0)
    fl_kappa = ls * (ls + 1) / 2.0; fl_kappa[0] = 0
    hp.almxfl(alm_dat, fl_kappa * wiener, inplace=True)
    kappa_map = hp.alm2map(alm_dat, NSIDE_WORK)
    kappa_std = kappa_map.std()
    kappa_map = (kappa_map - kappa_map.mean()) / kappa_std if kappa_std > 0 else kappa_map
    kappa_mask = None
    del alm_dat; gc.collect()
    print(f"  κ map: std={kappa_map.std():.4f}, range=[{kappa_map.min():.4f}, {kappa_map.max():.4f}]")
else:
    kappa_map, kappa_mask = try_download_kappa_alm(MAP_DIR)
have_kappa = kappa_map is not None
if have_kappa:
    nside_kappa = hp.npix2nside(len(kappa_map))
    print(f"  κ map loaded, NSIDE={nside_kappa}, range=[{kappa_map.min():.4f}, {kappa_map.max():.4f}]")
else:
    print("  WARNING: Failed to get κ map")

# Download y
print("\n=== Downloading Planck tSZ (y) map ===")
ymap, ymask = try_download_ysz(MAP_DIR)
have_y = ymap is not None
if have_y:
    nside_y = hp.npix2nside(len(ymap))
    print(f"  y map loaded, NSIDE={nside_y}, range=[{ymap.min():.6f}, {ymap.max():.6f}]")
else:
    print("  WARNING: Failed to get y map")

# Clean up any extracted directories
for p in MAP_DIR.iterdir():
    if p.is_dir():
        import shutil
        shutil.rmtree(p)
gc.collect()

# ── Step 3: Cross-match SNe with maps ─────────────────────────────────────
print("\n=== Cross-matching SNe with maps ===")

RADII_ARCMIN = [5, 10, 20, 30]

def disc_average(hpmap, theta, phi, radius_arcmin, nside):
    """Average map value in a disc of given radius."""
    vec = hp.ang2vec(theta, phi)
    radius_rad = np.radians(radius_arcmin / 60.0)
    pixels = hp.query_disc(nside, vec, radius_rad)
    if len(pixels) == 0:
        return np.nan, 0
    vals = hpmap[pixels]
    # Exclude UNSEEN
    good = vals != hp.UNSEEN
    if good.sum() == 0:
        return np.nan, 0
    return np.mean(vals[good]), int(good.sum())

def compensated_aperture(hpmap, theta, phi, r_inner, r_outer, nside):
    """Compensated aperture: mean(inner) - mean(annulus)."""
    vec = hp.ang2vec(theta, phi)
    r_in_rad = np.radians(r_inner / 60.0)
    r_out_rad = np.radians(r_outer / 60.0)
    pix_inner = set(hp.query_disc(nside, vec, r_in_rad))
    pix_outer = set(hp.query_disc(nside, vec, r_out_rad))
    pix_annulus = pix_outer - pix_inner

    def mean_good(pixels):
        if not pixels:
            return np.nan
        vals = hpmap[list(pixels)]
        good = vals != hp.UNSEEN
        return np.mean(vals[good]) if good.sum() > 0 else np.nan

    return mean_good(pix_inner) - mean_good(pix_annulus)

# Convert Galactic coords to theta, phi for healpy
theta = np.radians(90.0 - sn['GLAT'].values)  # colatitude
phi = np.radians(sn['GLON'].values)

# Initialize columns
for r in RADII_ARCMIN:
    sn[f'kappa_{r}'] = np.nan
    sn[f'y_{r}'] = np.nan
for pair in [(5,10), (10,20), (20,40)]:
    sn[f'kappa_comp_{pair[0]}_{pair[1]}'] = np.nan
    sn[f'y_comp_{pair[0]}_{pair[1]}'] = np.nan
sn['mask_ok_kappa'] = False
sn['mask_ok_y'] = False
sn['n_pix_kappa'] = 0
sn['n_pix_y'] = 0

for i in range(len(sn)):
    if i % 200 == 0:
        print(f"  Processing SN {i}/{len(sn)}...")
    th, ph = theta[i], phi[i]

    if have_kappa:
        ns = hp.npix2nside(len(kappa_map))
        for r in RADII_ARCMIN:
            val, npix = disc_average(kappa_map, th, ph, r, ns)
            sn.loc[sn.index[i], f'kappa_{r}'] = val
            if r == 10:
                sn.loc[sn.index[i], 'n_pix_kappa'] = npix
        for r_in, r_out in [(5,10), (10,20), (20,40)]:
            sn.loc[sn.index[i], f'kappa_comp_{r_in}_{r_out}'] = compensated_aperture(
                kappa_map, th, ph, r_in, r_out, ns)
        sn.loc[sn.index[i], 'mask_ok_kappa'] = not np.isnan(sn.loc[sn.index[i], 'kappa_10'])

    if have_y:
        ns = hp.npix2nside(len(ymap))
        for r in RADII_ARCMIN:
            val, npix = disc_average(ymap, th, ph, r, ns)
            sn.loc[sn.index[i], f'y_{r}'] = val
            if r == 10:
                sn.loc[sn.index[i], 'n_pix_y'] = npix
        for r_in, r_out in [(5,10), (10,20), (20,40)]:
            sn.loc[sn.index[i], f'y_comp_{r_in}_{r_out}'] = compensated_aperture(
                ymap, th, ph, r_in, r_out, ns)
        sn.loc[sn.index[i], 'mask_ok_y'] = not np.isnan(sn.loc[sn.index[i], 'y_10'])

print("  Cross-match complete.")

# ── Step 4: Compute Hubble residuals ──────────────────────────────────────
print("\n=== Computing Hubble residuals ===")
from scipy.integrate import quad

H0 = 73.04  # km/s/Mpc
Om = 0.334
c_km_s = 299792.458

def E(z):
    return np.sqrt(Om * (1+z)**3 + (1-Om))

def luminosity_distance(z):
    """d_L in Mpc for flat ΛCDM."""
    if z <= 0:
        return np.nan
    integral, _ = quad(lambda zp: 1.0/E(zp), 0, z)
    return (c_km_s / H0) * (1+z) * integral

def mu_LCDM(z):
    dL = luminosity_distance(z)
    if np.isnan(dL) or dL <= 0:
        return np.nan
    return 5 * np.log10(dL * 1e6 / 10)  # dL in pc / 10

print("  Computing ΛCDM distance moduli...")
sn['mu_LCDM'] = sn['zHD'].apply(mu_LCDM)
sn['delta_mu'] = sn['MU'] - sn['mu_LCDM']

# Deconfounded residual: regress out MWEBV, host mass step, survey ID
print("  Computing deconfounded residuals...")
valid = sn['delta_mu'].notna() & sn['zHD'].notna() & (sn['zHD'] > 0)
sn_valid = sn[valid].copy()

# Build design matrix
X_cols = []
sn_valid['const'] = 1.0
X_cols.append('const')
sn_valid['MWEBV_reg'] = sn_valid['MWEBV']
X_cols.append('MWEBV_reg')
# Host mass step
sn_valid['mass_step'] = (sn_valid['HOST_LOGMASS'] > 10).astype(float)
sn_valid.loc[sn_valid['HOST_LOGMASS'] < 0, 'mass_step'] = 0  # missing
X_cols.append('mass_step')
# Survey dummies (top surveys)
surveys = sn_valid['IDSURVEY'].value_counts()
for sid in surveys.index[1:min(8, len(surveys))]:  # skip most common (reference)
    col = f'survey_{int(sid)}'
    sn_valid[col] = (sn_valid['IDSURVEY'] == sid).astype(float)
    X_cols.append(col)

X = sn_valid[X_cols].values
y_resid = sn_valid['delta_mu'].values
# OLS
try:
    beta = np.linalg.lstsq(X, y_resid, rcond=None)[0]
    predicted = X @ beta
    sn_valid['delta_mu_deconf'] = y_resid - predicted
    # Map back
    sn['delta_mu_deconf'] = np.nan
    sn.loc[valid, 'delta_mu_deconf'] = sn_valid['delta_mu_deconf'].values
    print(f"  Deconfounding coefficients: MWEBV={beta[1]:.3f}, mass_step={beta[2]:.3f}")
except Exception as e:
    print(f"  Deconfounding failed: {e}")
    sn['delta_mu_deconf'] = sn['delta_mu']

# Clean sample flag
sn['clean_sample'] = (
    sn['mask_ok_kappa'] &
    (sn['MWEBV'] < 0.15) &
    sn['delta_mu'].notna() &
    (sn['zHD'] > 0.01)
)
print(f"  Clean sample: {sn['clean_sample'].sum()} SNe")

# ── Step 5: Save output ──────────────────────────────────────────────────
print("\n=== Saving output ===")

out_cols = ['CID', 'zHD', 'zHDERR', 'MU', 'MUERR', 'RA', 'DEC', 'MWEBV', 'HOST_LOGMASS', 'IDSURVEY',
            'kappa_5', 'kappa_10', 'kappa_20', 'kappa_30',
            'kappa_comp_5_10', 'kappa_comp_10_20', 'kappa_comp_20_40',
            'y_5', 'y_10', 'y_20', 'y_30',
            'y_comp_5_10', 'y_comp_10_20', 'y_comp_20_40',
            'mask_ok_kappa', 'mask_ok_y', 'n_pix_kappa', 'n_pix_y',
            'delta_mu', 'delta_mu_deconf', 'clean_sample']
sn[out_cols].to_csv(OUT_DIR / "sn_closure_table.csv", index=False)
print(f"  Saved sn_closure_table.csv ({len(sn)} rows)")

# Manifest
manifest = {
    "kappa_map": "Planck PR3 Lensing COM_Lensing_4096_R3.00" if have_kappa else "NOT AVAILABLE",
    "y_map": "Planck PR2 tSZ COM_CompMap_YSZ_R2.01" if have_y else "NOT AVAILABLE",
    "NSIDE_kappa": int(hp.npix2nside(len(kappa_map))) if have_kappa else None,
    "NSIDE_y": int(hp.npix2nside(len(ymap))) if have_y else None,
    "radii_arcmin": RADII_ARCMIN,
    "compensated_pairs": [[5,10],[10,20],[20,40]],
    "H0": H0, "Omega_m": Om,
    "n_sne_total": len(sn),
    "n_clean": int(sn['clean_sample'].sum()),
}
with open(OUT_DIR / "maps_used_manifest.json", 'w') as f:
    json.dump(manifest, f, indent=2)

# Selection report
with open(OUT_DIR / "selection_report.md", 'w') as f:
    f.write("# Closure Theory Selection Report\n\n")
    f.write(f"- Total SNe (unique): {len(sn)}\n")
    f.write(f"- With valid κ: {sn['mask_ok_kappa'].sum()}\n")
    f.write(f"- With valid y: {sn['mask_ok_y'].sum()}\n")
    f.write(f"- MWEBV < 0.15: {(sn['MWEBV'] < 0.15).sum()}\n")
    f.write(f"- z > 0.01: {(sn['zHD'] > 0.01).sum()}\n")
    f.write(f"- Clean sample (all cuts): {sn['clean_sample'].sum()}\n")
    f.write(f"\n## Maps\n")
    f.write(f"- κ: {'loaded' if have_kappa else 'MISSING'}\n")
    f.write(f"- y: {'loaded' if have_y else 'MISSING'}\n")

# ── Step 6: Diagnostic plots ─────────────────────────────────────────────
print("\n=== Making diagnostic plots ===")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

clean = sn[sn['clean_sample']].copy()
print(f"  Plotting with {len(clean)} clean SNe")

def plot_binned(ax, x, y, nbins=10, xlabel='', ylabel='', title=''):
    """Bin y by deciles of x and plot."""
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 20:
        ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes, ha='center')
        return
    bins = np.percentile(x, np.linspace(0, 100, nbins+1))
    bin_centers, bin_means, bin_errs = [], [], []
    for j in range(nbins):
        mask = (x >= bins[j]) & (x < bins[j+1]) if j < nbins-1 else (x >= bins[j])
        if mask.sum() > 2:
            bin_centers.append(np.mean(x[mask]))
            bin_means.append(np.mean(y[mask]))
            bin_errs.append(np.std(y[mask]) / np.sqrt(mask.sum()))
    ax.errorbar(bin_centers, bin_means, yerr=bin_errs, fmt='o-', capsize=3, color='navy')
    ax.axhline(0, ls='--', color='gray', alpha=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    # Linear fit
    if len(bin_centers) >= 3:
        from numpy.polynomial.polynomial import polyfit
        c = polyfit(x, y, 1)
        ax.plot([min(x), max(x)], [c[0]+c[1]*min(x), c[0]+c[1]*max(x)],
                'r--', alpha=0.5, label=f'slope={c[1]:.4f}')
        ax.legend(fontsize=8)

# 1. Elephant kappa
if have_kappa:
    fig, ax = plt.subplots(figsize=(8,5))
    plot_binned(ax, clean['kappa_comp_10_20'].values, clean['delta_mu_deconf'].values,
                xlabel='κ_comp (10-20 arcmin)', ylabel='Δμ⊥', title='Hubble residual vs CMB lensing κ')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'elephant_kappa.png', dpi=150)
    plt.close()
    print("  Saved elephant_kappa.png")

# 2. Elephant y
if have_y:
    fig, ax = plt.subplots(figsize=(8,5))
    plot_binned(ax, clean['y_comp_10_20'].values, clean['delta_mu_deconf'].values,
                xlabel='y_comp (10-20 arcmin)', ylabel='Δμ⊥', title='Hubble residual vs tSZ y')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'elephant_y.png', dpi=150)
    plt.close()
    print("  Saved elephant_y.png")

# 3. Scale test
if have_kappa:
    fig, ax = plt.subplots(figsize=(8,5))
    slopes = []
    for r in RADII_ARCMIN:
        x = clean[f'kappa_{r}'].values
        y = clean['delta_mu_deconf'].values
        valid = np.isfinite(x) & np.isfinite(y)
        if valid.sum() > 20:
            c = np.polyfit(x[valid], y[valid], 1)
            slopes.append((r, c[0]))
    if slopes:
        rs, ss = zip(*slopes)
        ax.plot(rs, ss, 'o-', color='darkred', markersize=8)
        ax.set_xlabel('Aperture radius (arcmin)')
        ax.set_ylabel('Slope dΔμ⊥/dκ')
        ax.set_title('Scale dependence of κ-Δμ correlation')
        ax.axhline(0, ls='--', color='gray')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'scale_test.png', dpi=150)
    plt.close()
    print("  Saved scale_test.png")

# 4. Redshift split
if have_kappa:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    z_bins = [(0.01, 0.1, 'z < 0.1'), (0.1, 0.5, '0.1 < z < 0.5'), (0.5, 3.0, 'z > 0.5')]
    for ax, (zlo, zhi, label) in zip(axes, z_bins):
        sub = clean[(clean['zHD'] >= zlo) & (clean['zHD'] < zhi)]
        plot_binned(ax, sub['kappa_comp_10_20'].values, sub['delta_mu_deconf'].values,
                    xlabel='κ_comp (10-20\')', ylabel='Δμ⊥',
                    title=f'{label} (N={len(sub)})')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'redshift_split.png', dpi=150)
    plt.close()
    print("  Saved redshift_split.png")

# ── Summary ───────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("PIPELINE COMPLETE")
print("="*60)
print(f"Output: {OUT_DIR / 'sn_closure_table.csv'}")
print(f"Manifest: {OUT_DIR / 'maps_used_manifest.json'}")
print(f"Report: {OUT_DIR / 'selection_report.md'}")
print(f"Plots: {PLOT_DIR}")
print(f"\nSummary:")
print(f"  Total unique SNe: {len(sn)}")
print(f"  Clean sample: {sn['clean_sample'].sum()}")
if have_kappa:
    c = clean[['kappa_comp_10_20', 'delta_mu_deconf']].dropna()
    if len(c) > 10:
        corr = np.corrcoef(c['kappa_comp_10_20'], c['delta_mu_deconf'])[0,1]
        slope = np.polyfit(c['kappa_comp_10_20'], c['delta_mu_deconf'], 1)[0]
        print(f"  κ_comp vs Δμ⊥: r={corr:.4f}, slope={slope:.4f}")
if have_y:
    c = clean[['y_comp_10_20', 'delta_mu_deconf']].dropna()
    if len(c) > 10:
        corr = np.corrcoef(c['y_comp_10_20'], c['delta_mu_deconf'])[0,1]
        slope = np.polyfit(c['y_comp_10_20'], c['delta_mu_deconf'], 1)[0]
        print(f"  y_comp vs Δμ⊥: r={corr:.4f}, slope={slope:.4f}")
