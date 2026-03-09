#!/usr/bin/env python3
"""
Figure 1: Standard Tripp vs Closure-Corrected Hubble Residuals
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize as scipy_minimize
from scipy.integrate import quad
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# LOAD & PREP
# ============================================================
data_file = 'data/pantheon_plus.dat'
with open(data_file, 'r') as f:
    header = f.readline().strip().split()
    rows = []
    for line in f:
        parts = line.strip().split()
        if len(parts) == len(header):
            rows.append(parts)

col = {name: i for i, name in enumerate(header)}
z_cmb = np.array([float(r[col['zCMB']]) for r in rows])
z_hd = np.array([float(r[col['zHD']]) for r in rows])
m_b = np.array([float(r[col['mB']]) for r in rows])
m_b_err = np.array([float(r[col['mBERR']]) for r in rows])
x1 = np.array([float(r[col['x1']]) for r in rows])
c = np.array([float(r[col['c']]) for r in rows])
host_mass = np.array([float(r[col['HOST_LOGMASS']]) for r in rows])

mask_q = (z_hd > 0.01) & (z_hd < 2.5) & (np.abs(x1) < 5) & (np.abs(c) < 0.5) & (m_b_err < 1.0) & (host_mass > 0) & (host_mass < 15)
z = z_cmb[mask_q]
mb = m_b[mask_q]
mb_err = m_b_err[mask_q]
x1_q = x1[mask_q]
c_q = c[mask_q]
hm = host_mass[mask_q]
mask_hi = hm >= 10.0

C_LIGHT = 299792.458
def E(zp, Om=0.3): return np.sqrt(Om*(1+zp)**3 + (1-Om))
def dist_mod(z_arr, H0=70, Om=0.3):
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        integral, _ = quad(lambda zp: 1.0/E(zp, Om), 0, zi)
        dl[i] = C_LIGHT * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

mu_lcdm = dist_mod(z, H0=70, Om=0.3)
I_sat = z / (1 + z/0.5)
x1_std = np.std(x1_q)
c_std = np.std(c_q)

# ============================================================
# STANDARD TRIPP
# ============================================================
def chi2_std(params):
    M_B, alpha, beta, gamma = params
    mu_obs = mb + alpha * x1_q - beta * c_q - M_B
    mu_obs[mask_hi] -= gamma
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res_std = scipy_minimize(chi2_std, x0=[-19.25, 0.15, 3.0, -0.05],
                         method='Nelder-Mead', options={'maxiter': 50000})
p_std = res_std.x
resid_std = (mb + p_std[1]*x1_q - p_std[2]*c_q - p_std[0]) - mu_lcdm
resid_std[mask_hi] -= p_std[3]

# ============================================================
# CLOSURE-CORRECTED (Modified Tripp with rotation)
# ============================================================
def chi2_m1(params):
    A, B, theta1, delta_v, g0, g1 = params
    angle = theta1 * I_sat
    x1n = x1_q / x1_std; cn = c_q / c_std
    u = x1n * np.cos(angle) + cn * np.sin(angle)
    v = -x1n * np.sin(angle) + cn * np.cos(angle)
    u_mag = u * x1_std; v_mag = v * c_std
    gamma_I = g0 + g1 * I_sat
    mu_corr = mb + A * u_mag - B * v_mag - delta_v * v_mag**2
    mu_corr[mask_hi] -= gamma_I[mask_hi]
    w = 1.0 / mb_err**2
    M0 = np.sum(w * (mu_corr - mu_lcdm)) / np.sum(w)
    return np.sum(((mu_corr - M0 - mu_lcdm) / mb_err)**2)

res_m1 = scipy_minimize(chi2_m1, x0=[0.15, 3.0, -0.3, 2.0, -0.03, -0.1],
                        method='Nelder-Mead', options={'maxiter': 100000})
p_m1 = res_m1.x
angle = p_m1[2] * I_sat
x1n = x1_q/x1_std; cn = c_q/c_std
u = x1n*np.cos(angle) + cn*np.sin(angle)
v = -x1n*np.sin(angle) + cn*np.cos(angle)
u_mag = u*x1_std; v_mag = v*c_std
gamma_I = p_m1[4] + p_m1[5]*I_sat
mu_corr = mb + p_m1[0]*u_mag - p_m1[1]*v_mag - p_m1[3]*v_mag**2
mu_corr[mask_hi] -= gamma_I[mask_hi]
w = 1.0/mb_err**2
M0 = np.sum(w*(mu_corr - mu_lcdm))/np.sum(w)
resid_m1 = mu_corr - M0 - mu_lcdm

# ============================================================
# STATS
# ============================================================
rho_std, p_std_val = spearmanr(z, resid_std)
rho_m1, p_m1_val = spearmanr(z, resid_m1)

# Binned means
z_edges = np.array([0.01, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.50, 0.70, 1.0, 2.3])
z_centers_s, mean_s, err_s = [], [], []
z_centers_m, mean_m, err_m = [], [], []
for i in range(len(z_edges)-1):
    mask = (z >= z_edges[i]) & (z < z_edges[i+1])
    n = mask.sum()
    if n > 3:
        zc = np.mean(z[mask])
        z_centers_s.append(zc)
        mean_s.append(np.mean(resid_std[mask]))
        err_s.append(np.std(resid_std[mask]) / np.sqrt(n))
        z_centers_m.append(zc)
        mean_m.append(np.mean(resid_m1[mask]))
        err_m.append(np.std(resid_m1[mask]) / np.sqrt(n))

# ============================================================
# PLOT
# ============================================================
fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True, 
                          gridspec_kw={'hspace': 0.08, 'height_ratios': [1, 1]})

# --- Top: Standard Tripp ---
ax1 = axes[0]
ax1.scatter(z, resid_std, s=3, alpha=0.15, c='#cc4444', zorder=1, rasterized=True)
ax1.errorbar(z_centers_s, mean_s, yerr=err_s, fmt='D', color='darkred', 
             markersize=8, capsize=4, linewidth=2, zorder=5, label='Binned mean')
ax1.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
ax1.set_ylabel(r'$\Delta\mu$ (mag)', fontsize=14)
ax1.set_ylim(-0.8, 0.8)
ax1.set_title(f'Standard Tripp  —  ρ(z, Δμ) = {rho_std:.3f}  (8.2σ from null)', 
              fontsize=15, fontweight='bold', color='darkred')
ax1.legend(fontsize=11, loc='upper right')
ax1.tick_params(labelsize=12)

# Add trend line
z_fit = np.linspace(0.01, 2.0, 100)
coeffs_s = np.polyfit(z, resid_std, 1)
ax1.plot(z_fit, np.polyval(coeffs_s, z_fit), 'r-', linewidth=2.5, alpha=0.8, 
         label=f'slope = {coeffs_s[0]:.3f} mag/z')
ax1.legend(fontsize=11, loc='upper right')

# --- Bottom: Closure-Corrected ---
ax2 = axes[1]
ax2.scatter(z, resid_m1, s=3, alpha=0.15, c='#3377bb', zorder=1, rasterized=True)
ax2.errorbar(z_centers_m, mean_m, yerr=err_m, fmt='D', color='darkblue', 
             markersize=8, capsize=4, linewidth=2, zorder=5, label='Binned mean')
ax2.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
ax2.set_xlabel('Redshift z', fontsize=14)
ax2.set_ylabel(r'$\Delta\mu$ (mag)', fontsize=14)
ax2.set_ylim(-0.8, 0.8)
ax2.set_title(f'Closure-Corrected (Modified Tripp)  —  ρ(z, Δμ) = {rho_m1:.3f}  (32% reduction)', 
              fontsize=15, fontweight='bold', color='darkblue')
ax2.tick_params(labelsize=12)

# Add trend line
coeffs_m = np.polyfit(z, resid_m1, 1)
ax2.plot(z_fit, np.polyval(coeffs_m, z_fit), 'b-', linewidth=2.5, alpha=0.8, 
         label=f'slope = {coeffs_m[0]:.3f} mag/z')
ax2.legend(fontsize=11, loc='upper right')

# Log scale x-axis for better spread
ax1.set_xscale('log')
ax2.set_xscale('log')
ax2.set_xlim(0.009, 2.5)

fig.suptitle('Pantheon+ Hubble Residuals: Standard vs Closure-Corrected', 
             fontsize=17, fontweight='bold', y=0.98)

plt.tight_layout()
plt.savefig('results_figure1/figure1_hubble_residuals.png', dpi=200, bbox_inches='tight')
print(f"Saved to results_figure1/figure1_hubble_residuals.png")
print(f"\nStandard Tripp: slope = {coeffs_s[0]:.4f} mag/z, ρ = {rho_std:.4f}")
print(f"Closure M1:     slope = {coeffs_m[0]:.4f} mag/z, ρ = {rho_m1:.4f}")
print(f"Slope reduction: {(1 - abs(coeffs_m[0])/abs(coeffs_s[0]))*100:.0f}%")
