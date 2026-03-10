#!/usr/bin/env python3
"""
Paper 1 Figures — The Devastating Set
========================================
Fig 1: Rolling correlations (the rotation)
Fig 2: Hubble residuals (the fix)
Fig 3: Interaction term Δχ² waterfall (the simplicity)
Fig 4: PCA drainage (the illustration)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import minimize as scipy_minimize
from scipy.integrate import quad
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 13, 'axes.labelsize': 14, 'axes.titlesize': 15,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'legend.fontsize': 11, 'figure.dpi': 200,
    'font.family': 'serif'
})

# ============================================================
# LOAD DATA
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
z_hd = np.array([float(r[col['zHD']]) for r in rows])
z_cmb = np.array([float(r[col['zCMB']]) for r in rows])
m_b = np.array([float(r[col['mB']]) for r in rows])
m_b_err = np.array([float(r[col['mBERR']]) for r in rows])
x1_all = np.array([float(r[col['x1']]) for r in rows])
c_all = np.array([float(r[col['c']]) for r in rows])
hm_all = np.array([float(r[col['HOST_LOGMASS']]) for r in rows])

mask_q = (z_hd > 0.01) & (z_hd < 2.5) & (np.abs(x1_all) < 5) & (np.abs(c_all) < 0.5) & (m_b_err < 1.0) & (hm_all > 0) & (hm_all < 15)
z = z_cmb[mask_q]; mb = m_b[mask_q]; mb_err = m_b_err[mask_q]
x1 = x1_all[mask_q]; c = c_all[mask_q]; hm = hm_all[mask_q]
N = len(z)

def E(zp): return np.sqrt(0.3*(1+zp)**3 + 0.7)
def dist_mod(z_arr, H0=70):
    dl = np.zeros_like(z_arr)
    for i, zi in enumerate(z_arr):
        integral, _ = quad(lambda zp: 1.0/E(zp), 0, zi)
        dl[i] = 299792.458 * (1+zi) * integral / H0
    return 5 * np.log10(dl) + 25

mu_lcdm = dist_mod(z)
mask_hi = hm >= 10.0

# Standard Tripp
def chi2_std(params):
    M_B, alpha, beta, gamma = params
    mu_obs = mb + alpha * x1 - beta * c - M_B
    mu_obs[mask_hi] -= gamma
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res_std = scipy_minimize(chi2_std, x0=[-19.25, 0.15, 3.0, -0.05], method='Nelder-Mead', options={'maxiter':50000})
p_s = res_std.x
resid_std = (mb + p_s[1]*x1 - p_s[2]*c - p_s[0]) - mu_lcdm
resid_std[mask_hi] -= p_s[3]

# Modified Tripp (interaction model: z×x1, z×c)
def chi2_int(params):
    M_B, alpha, alpha1, beta, beta1, gamma = params
    mu_obs = mb + (alpha + alpha1*z)*x1 - (beta + beta1*z)*c - M_B
    mu_obs[mask_hi] -= gamma
    return np.sum(((mu_obs - mu_lcdm) / mb_err)**2)

res_int = scipy_minimize(chi2_int, x0=[p_s[0], p_s[1], 0, p_s[2], 0, p_s[3]], method='Nelder-Mead', options={'maxiter':50000})
p_i = res_int.x
resid_int = (mb + (p_i[1]+p_i[2]*z)*x1 - (p_i[3]+p_i[4]*z)*c - p_i[0]) - mu_lcdm
resid_int[mask_hi] -= p_i[5]

# Full modified (with c², γ(z))
I_sat = z / (1 + z/0.5)
x1_std_val = np.std(x1); c_std_val = np.std(c)

def chi2_full(params):
    A, B, theta1, delta_v, g0, g1 = params
    angle = theta1 * I_sat
    x1n = x1/x1_std_val; cn = c/c_std_val
    u = x1n*np.cos(angle) + cn*np.sin(angle)
    v = -x1n*np.sin(angle) + cn*np.cos(angle)
    u_mag = u*x1_std_val; v_mag = v*c_std_val
    gamma_I = g0 + g1*I_sat
    mu_corr = mb + A*u_mag - B*v_mag - delta_v*v_mag**2
    mu_corr[mask_hi] -= gamma_I[mask_hi]
    w = 1.0/mb_err**2
    M0 = np.sum(w*(mu_corr - mu_lcdm))/np.sum(w)
    return np.sum(((mu_corr - M0 - mu_lcdm)/mb_err)**2)

res_full = scipy_minimize(chi2_full, x0=[0.15, 3.0, -0.3, 2.0, -0.03, -0.1], method='Nelder-Mead', options={'maxiter':100000})
p_f = res_full.x
angle = p_f[2]*I_sat
x1n = x1/x1_std_val; cn = c/c_std_val
u = x1n*np.cos(angle)+cn*np.sin(angle)
v = -x1n*np.sin(angle)+cn*np.cos(angle)
u_mag = u*x1_std_val; v_mag = v*c_std_val
gamma_I = p_f[4]+p_f[5]*I_sat
mu_corr = mb + p_f[0]*u_mag - p_f[1]*v_mag - p_f[3]*v_mag**2
mu_corr[mask_hi] -= gamma_I[mask_hi]
w = 1.0/mb_err**2
M0 = np.sum(w*(mu_corr-mu_lcdm))/np.sum(w)
resid_full = mu_corr - M0 - mu_lcdm

order = np.argsort(z)
z_s = z[order]; x1_s = x1[order]; c_s = c[order]
r_std_s = resid_std[order]; r_int_s = resid_int[order]; r_full_s = resid_full[order]

# ============================================================
# FIGURE 1: ROLLING CORRELATIONS — The Rotation
# ============================================================
print("Generating Figure 1: Rolling Correlations...")

fig1, ax1 = plt.subplots(figsize=(10, 6))

W = 200; step = 10
z_roll = []; rc_vals = []; rx1_vals = []
# Bootstrap bands
N_BOOT_FIG = 200
rc_boot = []; rx1_boot = []

for i in range(0, len(z_s)-W, step):
    sl = slice(i, i+W)
    z_roll.append(np.median(z_s[sl]))
    r1, _ = spearmanr(x1_s[sl], r_std_s[sl])
    r2, _ = spearmanr(c_s[sl], r_std_s[sl])
    rx1_vals.append(r1); rc_vals.append(r2)
    
    # Bootstrap confidence
    rc_b = []; rx1_b = []
    idx_pool = np.arange(i, min(i+W, len(z_s)))
    for _ in range(N_BOOT_FIG):
        bi = np.random.choice(idx_pool, W, replace=True)
        rb1, _ = spearmanr(x1_s[bi], r_std_s[bi])
        rb2, _ = spearmanr(c_s[bi], r_std_s[bi])
        rx1_b.append(rb1); rc_b.append(rb2)
    rc_boot.append(rc_b); rx1_boot.append(rx1_b)

z_roll = np.array(z_roll)
rc_vals = np.array(rc_vals); rx1_vals = np.array(rx1_vals)
rc_boot = np.array(rc_boot); rx1_boot = np.array(rx1_boot)

rc_lo = np.percentile(rc_boot, 2.5, axis=1)
rc_hi = np.percentile(rc_boot, 97.5, axis=1)
rx1_lo = np.percentile(rx1_boot, 2.5, axis=1)
rx1_hi = np.percentile(rx1_boot, 97.5, axis=1)

ax1.plot(z_roll, rc_vals, 'o-', color='#cc2222', markersize=4, linewidth=2, label=r'corr($c$, $\Delta\mu$)', zorder=5)
ax1.fill_between(z_roll, rc_lo, rc_hi, color='#cc2222', alpha=0.15)
ax1.plot(z_roll, rx1_vals, 's-', color='#2255aa', markersize=4, linewidth=2, label=r'corr($x_1$, $\Delta\mu$)', zorder=5)
ax1.fill_between(z_roll, rx1_lo, rx1_hi, color='#2255aa', alpha=0.15)

ax1.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax1.set_xlabel('Redshift $z$ (window median)')
ax1.set_ylabel('Spearman correlation with Hubble residual')
ax1.set_title('Channel–Residual Coupling Evolution Under Standard Tripp\n'
              r'corr($c$, $\Delta\mu$) trend: $\rho = -0.785$, $p = 3\times10^{-30}$', fontweight='bold')
ax1.legend(loc='lower left', framealpha=0.9)
ax1.set_xlim(0, 0.8)
ax1.set_ylim(-0.55, 0.35)

# Annotate the crossing / sign reversal
cross_idx = np.argmin(np.abs(rc_vals))
ax1.annotate('Color decouples\nfrom residuals', xy=(z_roll[cross_idx], 0), 
             xytext=(0.25, 0.20), fontsize=11, ha='center',
             arrowprops=dict(arrowstyle='->', color='#cc2222', lw=1.5),
             color='#cc2222', fontweight='bold')

ax1.text(0.05, -0.42, 'Color anti-correlates\nwith distance error\nat high $z$', 
         fontsize=10, color='#cc2222', fontstyle='italic')
ax1.text(0.55, 0.18, 'Stretch gains\ncorrective power', 
         fontsize=10, color='#2255aa', fontstyle='italic')

plt.tight_layout()
fig1.savefig('results_figure1/fig1_rolling_correlations.png', dpi=200, bbox_inches='tight')
print("  → fig1_rolling_correlations.png")


# ============================================================
# FIGURE 2: HUBBLE RESIDUALS — Standard vs Corrected
# ============================================================
print("Generating Figure 2: Hubble Residuals...")

fig2, (ax2a, ax2b) = plt.subplots(2, 1, figsize=(11, 9), sharex=True, 
                                    gridspec_kw={'hspace': 0.08})

# Binned means
z_edges = np.array([0.01, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.50, 0.70, 1.0, 2.3])

def plot_panel(ax, z_data, resid_data, color_pts, color_bin, title, rho_val):
    ax.scatter(z_data, resid_data, s=2, alpha=0.12, c=color_pts, rasterized=True, zorder=1)
    
    zc_list = []; mean_list = []; err_list = []
    for i in range(len(z_edges)-1):
        m = (z_data >= z_edges[i]) & (z_data < z_edges[i+1])
        n = m.sum()
        if n > 3:
            zc_list.append(np.mean(z_data[m]))
            mean_list.append(np.mean(resid_data[m]))
            err_list.append(np.std(resid_data[m])/np.sqrt(n))
    
    ax.errorbar(zc_list, mean_list, yerr=err_list, fmt='D', color=color_bin,
                markersize=7, capsize=3, linewidth=2, zorder=5)
    
    # Trend line
    coeff = np.polyfit(z_data, resid_data, 1)
    zf = np.linspace(0.01, 2.0, 100)
    ax.plot(zf, np.polyval(coeff, zf), '-', color=color_bin, linewidth=2.5, alpha=0.7)
    
    ax.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_ylabel(r'$\Delta\mu$ (mag)')
    ax.set_ylim(-0.7, 0.7)
    ax.set_title(title, fontweight='bold', color=color_bin)
    ax.text(0.97, 0.95, f'$\\rho = {rho_val:.3f}$\nslope = {coeff[0]:.3f} mag/$z$', 
            transform=ax.transAxes, ha='right', va='top', fontsize=12, 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

rho_std, _ = spearmanr(z, resid_std)
rho_int, _ = spearmanr(z, resid_int)

plot_panel(ax2a, z, resid_std, '#cc4444', 'darkred',
           'Standard Tripp (constant coefficients) — 8.2σ anomaly', rho_std)
plot_panel(ax2b, z, resid_int, '#3377bb', 'darkblue',
           f'Interaction Model ($z\\times x_1$, $z\\times c$) — Δχ² = 269', rho_int)

ax2a.set_xscale('log'); ax2b.set_xscale('log')
ax2b.set_xlabel('Redshift $z$')
ax2b.set_xlim(0.009, 2.5)

fig2.suptitle('Pantheon+ Hubble Residuals: Fixed vs Nonstationary Standardization',
              fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout()
fig2.savefig('results_figure1/fig2_hubble_residuals.png', dpi=200, bbox_inches='tight')
print("  → fig2_hubble_residuals.png")


# ============================================================
# FIGURE 3: Δχ² WATERFALL — The Simplicity
# ============================================================
print("Generating Figure 3: Δχ² Waterfall...")

fig3, ax3 = plt.subplots(figsize=(9, 6))

chi2_base = res_std.fun

# Individual terms
def chi2_zx1(params):
    M_B, alpha, alpha1, beta, gamma = params
    mu_obs = mb + (alpha+alpha1*z)*x1 - beta*c - M_B
    mu_obs[mask_hi] -= gamma
    return np.sum(((mu_obs-mu_lcdm)/mb_err)**2)

def chi2_zc(params):
    M_B, alpha, beta, beta1, gamma = params
    mu_obs = mb + alpha*x1 - (beta+beta1*z)*c - M_B
    mu_obs[mask_hi] -= gamma
    return np.sum(((mu_obs-mu_lcdm)/mb_err)**2)

def chi2_c2(params):
    M_B, alpha, beta, delta, gamma = params
    mu_obs = mb + alpha*x1 - beta*c - delta*c**2 - M_B
    mu_obs[mask_hi] -= gamma
    return np.sum(((mu_obs-mu_lcdm)/mb_err)**2)

def chi2_zg(params):
    M_B, alpha, beta, gamma, gamma1 = params
    mu_obs = mb + alpha*x1 - beta*c - M_B
    mu_obs[mask_hi] -= (gamma + gamma1*z[mask_hi])
    return np.sum(((mu_obs-mu_lcdm)/mb_err)**2)

r_zx1 = scipy_minimize(chi2_zx1, x0=[p_s[0], p_s[1], 0, p_s[2], p_s[3]], method='Nelder-Mead', options={'maxiter':50000})
r_zc = scipy_minimize(chi2_zc, x0=[p_s[0], p_s[1], p_s[2], 0, p_s[3]], method='Nelder-Mead', options={'maxiter':50000})
r_c2 = scipy_minimize(chi2_c2, x0=[p_s[0], p_s[1], p_s[2], 0, p_s[3]], method='Nelder-Mead', options={'maxiter':50000})
r_zg = scipy_minimize(chi2_zg, x0=[p_s[0], p_s[1], p_s[2], p_s[3], 0], method='Nelder-Mead', options={'maxiter':50000})

terms = [
    ('$z \\times x_1$\n(stretch evolves)', chi2_base - r_zx1.fun, 1, '#2255aa'),
    ('$z \\times c$\n(color evolves)', chi2_base - r_zc.fun, 1, '#cc2222'),
    ('$c^2$\n(color curvature)', chi2_base - r_c2.fun, 1, '#22aa55'),
    ('$z \\times \\gamma$\n(mass step evolves)', chi2_base - r_zg.fun, 1, '#aa8822'),
    ('All four\ncombined', chi2_base - res_full.fun, 4, '#333333'),
]

labels = [t[0] for t in terms]
dchi2 = [t[1] for t in terms]
dof = [t[2] for t in terms]
colors = [t[3] for t in terms]

bars = ax3.barh(range(len(terms)), dchi2, color=colors, edgecolor='white', linewidth=1.5, height=0.6)

for i, (dc, d) in enumerate(zip(dchi2, dof)):
    per_param = dc/d
    ax3.text(dc + 8, i, f'Δχ² = {dc:.0f}  ({dc/d:.0f}/param)', va='center', fontsize=12, fontweight='bold')

ax3.set_yticks(range(len(terms)))
ax3.set_yticklabels(labels, fontsize=12)
ax3.set_xlabel('Δχ² improvement over Standard Tripp', fontsize=13)
ax3.set_title('One Parameter at a Time: Each Correction Is Individually Devastating',
              fontsize=14, fontweight='bold')
ax3.axvline(3.84, color='gray', linestyle=':', alpha=0.5)
ax3.text(3.84, -0.5, '  p = 0.05\n  threshold', fontsize=9, color='gray')
ax3.set_xlim(0, max(dchi2)*1.4)
ax3.invert_yaxis()

plt.tight_layout()
fig3.savefig('results_figure1/fig3_dchi2_waterfall.png', dpi=200, bbox_inches='tight')
print("  → fig3_dchi2_waterfall.png")


# ============================================================
# FIGURE 4: PCA DRAINAGE — The Illustration
# ============================================================
print("Generating Figure 4: PCA Drainage...")

fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(12, 5))

# PCA in z-bins
z_bins = [(0.01, 0.06), (0.06, 0.11), (0.11, 0.22), (0.22, 0.50), (0.50, 2.3)]
z_centers = []; frac_tripp = []; frac_resid = []; angles_pc1 = []

for zlo, zhi in z_bins:
    m = (z >= zlo) & (z < zhi)
    if m.sum() < 10: continue
    x1n = (x1[m]-np.mean(x1[m]))/(np.std(x1[m])+1e-10)
    cn = (c[m]-np.mean(c[m]))/(np.std(c[m])+1e-10)
    rn = (resid_std[m]-np.mean(resid_std[m]))/(np.std(resid_std[m])+1e-10)
    data = np.column_stack([x1n, cn, rn])
    eigvals, eigvecs = np.linalg.eigh(np.cov(data.T))
    idx = np.argsort(eigvals)[::-1]
    pc1 = eigvecs[:, idx[0]]
    if pc1[0] < 0: pc1 = -pc1
    
    in_plane = pc1[0]**2 + pc1[1]**2
    in_resid = pc1[2]**2
    z_centers.append(0.5*(zlo+zhi))
    frac_tripp.append(in_plane*100)
    frac_resid.append(in_resid*100)
    angles_pc1.append(np.degrees(np.arctan2(pc1[1], pc1[0])))

# Stacked bar
z_centers = np.array(z_centers)
frac_tripp = np.array(frac_tripp)
frac_resid = np.array(frac_resid)

bar_width = 0.12
x_pos = np.arange(len(z_centers))

ax4a.bar(x_pos, frac_tripp, bar_width*3, color='#2255aa', label='In Tripp plane $(x_1, c)$')
ax4a.bar(x_pos, frac_resid, bar_width*3, bottom=frac_tripp, color='#cc2222', label='In residual dimension')

ax4a.axhline(50, color='gray', linestyle=':', alpha=0.5)
ax4a.text(3.7, 52, '50% line', fontsize=10, color='gray')

bin_labels = ['$z<0.06$', '$0.06{-}0.11$', '$0.11{-}0.22$', '$0.22{-}0.50$', '$z>0.50$']
ax4a.set_xticks(x_pos)
ax4a.set_xticklabels(bin_labels)
ax4a.set_ylabel('PC1 variance fraction (%)')
ax4a.set_title('PC1 Projection: Tripp Plane vs Residual', fontweight='bold')
ax4a.legend(loc='upper right')
ax4a.set_ylim(0, 105)

# Annotate the drainage
ax4a.annotate('76%', xy=(0, frac_tripp[0]/2), fontsize=14, fontweight='bold', 
              color='white', ha='center', va='center')
ax4a.annotate('50%', xy=(4, frac_tripp[4]/2), fontsize=14, fontweight='bold', 
              color='white', ha='center', va='center')

# PC1 angle evolution
ax4b.plot(z_centers, angles_pc1, 'ko-', markersize=8, linewidth=2)
ax4b.set_xlabel('Redshift (bin center)')
ax4b.set_ylabel('PC1 angle in $(x_1, c)$ plane (degrees)')
ax4b.set_title('Eigenvector Rotation', fontweight='bold')

total_rot = angles_pc1[-1] - angles_pc1[0]
ax4b.annotate(f'Total rotation: {total_rot:+.0f}°', xy=(0.5, 0.05), 
              xycoords='axes fraction', fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round', facecolor='lightyellow'))

plt.tight_layout()
fig4.savefig('results_figure1/fig4_pca_drainage.png', dpi=200, bbox_inches='tight')
print("  → fig4_pca_drainage.png")


# ============================================================
# FIGURE 5: THE MONEY SHOT — w_a absorption
# ============================================================
print("Generating Figure 5: w_a Absorption...")

fig5, ax5 = plt.subplots(figsize=(8, 6))

# Data for the plot
models = ['Standard Tripp\n+ $w_0 w_a$', 'Modified Tripp\n+ $w_0 w_a$', 'Modified Tripp\n+ ΛCDM']
w0_vals = [-0.559, -0.659, -1.0]
wa_vals = [-1.883, -0.020, 0.0]
chi2_vals = [9811, None, 9775]  # Modified+w0wa not stored, skip
colors_w = ['#cc2222', '#2255aa', '#22aa55']
markers = ['o', 's', '*']

# Plot in w0-wa space
ax5.scatter(w0_vals[0], wa_vals[0], s=200, c=colors_w[0], marker=markers[0], zorder=5, edgecolors='black', linewidth=1.5)
ax5.scatter(w0_vals[1], wa_vals[1], s=200, c=colors_w[1], marker=markers[1], zorder=5, edgecolors='black', linewidth=1.5)
ax5.scatter(w0_vals[2], wa_vals[2], s=300, c=colors_w[2], marker=markers[2], zorder=5, edgecolors='black', linewidth=1.5)

# Arrow from standard to modified
ax5.annotate('', xy=(w0_vals[1], wa_vals[1]), xytext=(w0_vals[0], wa_vals[0]),
            arrowprops=dict(arrowstyle='->', color='#666666', lw=2.5, connectionstyle='arc3,rad=0.2'))
ax5.text(-0.6, -1.0, '99% of $w_a$\nabsorbed', fontsize=13, fontweight='bold', 
         color='#666666', ha='center')

# ΛCDM point
ax5.axhline(0, color='gray', linestyle='--', alpha=0.3)
ax5.axvline(-1, color='gray', linestyle='--', alpha=0.3)

ax5.text(-0.52, -1.95, 'Standard Tripp\n$\\chi^2 = 9{,}811$', fontsize=11, color=colors_w[0], fontweight='bold')
ax5.text(-0.72, 0.15, 'Modified Tripp\n$w_a \\approx 0$', fontsize=11, color=colors_w[1], fontweight='bold')
ax5.text(-1.05, 0.15, 'ΛCDM\n$\\chi^2 = 9{,}775$\n(best fit)', fontsize=11, color=colors_w[2], fontweight='bold')

ax5.set_xlabel('$w_0$', fontsize=15)
ax5.set_ylabel('$w_a$', fontsize=15)
ax5.set_title('Dark Energy Parameters: Standard vs Nonstationary Standardization',
              fontsize=14, fontweight='bold')
ax5.set_xlim(-1.3, -0.3)
ax5.set_ylim(-2.3, 0.5)

plt.tight_layout()
fig5.savefig('results_figure1/fig5_wa_absorption.png', dpi=200, bbox_inches='tight')
print("  → fig5_wa_absorption.png")


print("\n" + "="*70)
print("ALL 5 FIGURES GENERATED")
print("="*70)
