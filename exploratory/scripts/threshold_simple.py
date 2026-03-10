#!/usr/bin/env python3
"""
THRESHOLD FINDER — Simple. Let data speak.
"""
import numpy as np
from astropy.io import fits
from scipy.stats import spearmanr, kurtosis
import warnings
warnings.filterwarnings('ignore')

hdu = fits.open('/root/clawd/data/sdss/dr16q_prop.fits')
data = hdu[1].data
z = data['Z_DR16Q']

lines = [
    ('LYA', 'LYA', 1216), ('CIV', 'CIV', 1549), ('CIII', 'CIII_ALL', 1909),
    ('MGII', 'MGII', 2798), ('OII', 'OII3728', 3728), ('HBETA', 'HBETA', 4861),
    ('OIII', 'OIII5007', 5007), ('HALPHA', 'HALPHA', 6563),
    ('NII', 'NII6585', 6585), ('SII', 'SII6718', 6718),
]

def get(col_name):
    try:
        c = data[col_name]
        if len(c.shape) < 2 or c.shape[1] < 5: return None, None, None
        ew, fw = c[:,2], c[:,4]
        v = np.isfinite(ew) & np.isfinite(fw) & (ew>0) & (fw>0) & (ew<500) & (fw<30000) & (z>0.05)
        return ew, fw, v
    except: return None, None, None

# METHOD 1: EW/FWHM ratio gradient
print("=" * 70)
print("METHOD 1: EW/FWHM RATIO — where does it change fastest?")
print("=" * 70)

for name, col, wav in lines:
    ew, fw, v = get(col)
    if ew is None or v.sum() < 3000: continue
    zv, rv = z[v], ew[v]/fw[v]
    
    # 20 equal-count bins
    idx = np.argsort(zv)
    n = len(idx)
    bs = n // 20
    
    zmeds, rmeds = [], []
    for i in range(20):
        sl = idx[i*bs:(i+1)*bs]
        zmeds.append(np.median(zv[sl]))
        rmeds.append(np.median(rv[sl]))
    
    zmeds, rmeds = np.array(zmeds), np.array(rmeds)
    grad = np.gradient(rmeds, zmeds)
    
    r_t, p_t = spearmanr(zmeds, rmeds)
    drop_z = zmeds[np.argmin(grad)]
    
    print(f"\n{name:>8} ({wav}Å) N={v.sum():>7,} z=[{zv.min():.2f}–{zv.max():.2f}]")
    print(f"  ratio: {rmeds[0]:.4f} → {rmeds[-1]:.4f}  trend r={r_t:+.3f} p={p_t:.1e}")
    print(f"  steepest change at z = {drop_z:.3f}")

# METHOD 2: Running kurtosis
print("\n" + "=" * 70)
print("METHOD 2: KURTOSIS — tail fattening = stochastic damage")
print("=" * 70)

for name, col, wav in lines:
    ew, fw, v = get(col)
    if ew is None or v.sum() < 5000: continue
    zv, ewv, fwv = z[v], ew[v], fw[v]
    
    idx = np.argsort(zv)
    n = len(idx)
    bs = n // 15
    
    zmeds, kew, kfw = [], [], []
    for i in range(15):
        sl = idx[i*bs:(i+1)*bs]
        zmeds.append(np.median(zv[sl]))
        kew.append(kurtosis(ewv[sl]))
        kfw.append(kurtosis(fwv[sl]))
    
    zmeds = np.array(zmeds)
    kew, kfw = np.array(kew), np.array(kfw)
    
    r_ew, p_ew = spearmanr(zmeds, kew)
    r_fw, p_fw = spearmanr(zmeds, kfw)
    pk_ew = zmeds[np.argmax(kew)]
    
    print(f"\n{name:>8} ({wav}Å)")
    print(f"  EW kurt:   r={r_ew:+.3f} p={p_ew:.1e} peak={pk_ew:.3f}")
    print(f"  FWHM kurt: r={r_fw:+.3f} p={p_fw:.1e}")
    if r_ew > 0.3 and p_ew < 0.05 and (r_fw < 0.3 or p_fw > 0.05):
        print(f"  → SELECTIVE: EW tails fatten, FWHM doesn't ✓")

# METHOD 3: Inter-line coherence
print("\n" + "=" * 70)
print("METHOD 3: INTER-LINE EW COHERENCE")
print("=" * 70)

pairs = [
    ('MGII', 'CIII_ALL', 'MgII-CIII'),
    ('CIV', 'CIII_ALL', 'CIV-CIII'),
    ('HBETA', 'OIII5007', 'Hβ-OIII'),
    ('CIV', 'LYA', 'CIV-Lyα'),
]

for c1, c2, pname in pairs:
    try:
        d1, d2 = data[c1], data[c2]
        if len(d1.shape)<2 or len(d2.shape)<2: continue
        ew1, ew2 = d1[:,2], d2[:,2]
        v = np.isfinite(ew1) & np.isfinite(ew2) & (ew1>0) & (ew2>0) & (ew1<500) & (ew2<500) & (z>0.05)
        if v.sum() < 3000: continue
        
        zv, e1, e2 = z[v], ew1[v], ew2[v]
        idx = np.argsort(zv)
        n = len(idx)
        bs = n // 12
        
        zmeds, corrs = [], []
        for i in range(12):
            sl = idx[i*bs:(i+1)*bs]
            r, _ = spearmanr(e1[sl], e2[sl])
            zmeds.append(np.median(zv[sl]))
            corrs.append(r)
        
        zmeds, corrs = np.array(zmeds), np.array(corrs)
        r_t, p_t = spearmanr(zmeds, corrs)
        
        grad = np.gradient(corrs, zmeds)
        steepest = zmeds[np.argmin(grad)]
        
        print(f"\n{pname}: N={v.sum():,}")
        print(f"  coherence: {corrs[0]:.3f} → {corrs[-1]:.3f}  trend r={r_t:+.3f} p={p_t:.1e}")
        print(f"  steepest drop at z = {steepest:.3f}")
    except Exception as e:
        pass

# METHOD 4: RAW TABLE — MgII
print("\n" + "=" * 70)
print("METHOD 4: RAW TABLE — MGII")
print("=" * 70)

ew, fw, v = get('MGII')
zv, ewv, fwv = z[v], ew[v], fw[v]
idx = np.argsort(zv)
n = len(idx)
bs = n // 20

print(f"{'z':>6} {'N':>6} {'EW_med':>7} {'FW_med':>7} {'ratio':>7} {'EW_std':>7} {'FW_std':>7} {'EW_k':>6} {'FW_k':>6} {'r_EW_FW':>8}")
print("-" * 80)

for i in range(20):
    sl = idx[i*bs:(i+1)*bs]
    zbin = zv[sl]; ewbin = ewv[sl]; fwbin = fwv[sl]
    r_ef, _ = spearmanr(ewbin, fwbin)
    print(f"{np.median(zbin):6.3f} {len(sl):6d} {np.median(ewbin):7.2f} {np.median(fwbin):7.1f} "
          f"{np.median(ewbin)/np.median(fwbin):7.4f} {np.std(ewbin):7.2f} {np.std(fwbin):7.1f} "
          f"{kurtosis(ewbin):6.2f} {kurtosis(fwbin):6.2f} {r_ef:8.3f}")

# METHOD 5: RAW TABLE — CIV
print(f"\nCIV:")
ew, fw, v = get('CIV')
zv, ewv, fwv = z[v], ew[v], fw[v]
idx = np.argsort(zv)
n = len(idx)
bs = n // 15

print(f"{'z':>6} {'N':>6} {'EW_med':>7} {'FW_med':>7} {'ratio':>7} {'EW_std':>7} {'FW_std':>7} {'EW_k':>6} {'FW_k':>6} {'r_EW_FW':>8}")
print("-" * 80)

for i in range(15):
    sl = idx[i*bs:(i+1)*bs]
    zbin = zv[sl]; ewbin = ewv[sl]; fwbin = fwv[sl]
    r_ef, _ = spearmanr(ewbin, fwbin)
    print(f"{np.median(zbin):6.3f} {len(sl):6d} {np.median(ewbin):7.2f} {np.median(fwbin):7.1f} "
          f"{np.median(ewbin)/np.median(fwbin):7.4f} {np.std(ewbin):7.2f} {np.std(fwbin):7.1f} "
          f"{kurtosis(ewbin):6.2f} {kurtosis(fwbin):6.2f} {r_ef:8.3f}")

hdu.close()
print("\nDone.")
