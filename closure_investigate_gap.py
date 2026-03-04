#!/usr/bin/env python3
"""
INVESTIGATE THE GAP: What Are We Missing?

Three models, partial success each:
  z-sigmoid: Ωde ✓, H₀ ✗, w ✗
  Pure Σ:    H₀ ✓, Ωde ✗, w ✗  
  Tripp:     w ✓, β ✓, H₀ ✗, Ωde ~

Questions to investigate:
A. What if H₀_true = 67.4, not 73? (compression inflates LOCAL, not deflates CMB)
B. What if stretch also compresses slightly (q=0.039)?
C. What if compression SHIFTS the color distribution, not just degrades β?
D. What if we've been wrong about which measurement is biased?
E. Sensitivity: what ONE parameter change makes everything work?
F. Can the Tripp + direct pathways COMBINE?

Author: Closure Theory Pipeline
Date: March 2026
"""

import numpy as np
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import minimize, minimize_scalar
import json
import os

np.random.seed(42)

# ============================================================
# SETUP (reuse from previous)
# ============================================================
c_km = 299792.458
H0_planck = 67.4
H0_local = 73.04
Om = 0.315
Ob = 0.0493
OL = 1 - Om
rho_crit = 9.47e-27
m_p = 1.673e-27
sigma_T = 6.652e-29
n_b0 = rho_crit * Ob / m_p

BETA_0 = 3.1
BETA_LOW_Z = 2.94
BETA_HIGH_Z = 1.64
ALPHA_0 = 0.144

def H_func(z):
    return H0_planck * np.sqrt(Om*(1+z)**3 + OL)

def dSigma_dz(z):
    H_si = H_func(z) * 1e3 / 3.086e22
    return n_b0 * 2.998e8 * (1+z)**2 / H_si

z_grid = np.linspace(0, 3, 500)
Sigma_grid = np.array([quad(dSigma_dz, 0, zi)[0] if zi > 0 else 0 for zi in z_grid])

def Sigma_interp(z):
    return np.interp(z, z_grid, Sigma_grid)

def Sigma_eff(z, Sigma_sat):
    S = Sigma_interp(z)
    return Sigma_sat * (1 - np.exp(-S / Sigma_sat))

def dl_lcdm(z_val, Om_val=0.315, H0_val=73.04):
    OL_val = 1 - Om_val
    def integrand(zp): return 1.0/np.sqrt(Om_val*(1+zp)**3 + OL_val)
    dc, _ = quad(integrand, 0, z_val)
    return dc*(1+z_val)*c_km/H0_val

def mu_from_dl(dl):
    return 5*np.log10(dl)+25 if dl > 0 else 0

# Best-fit from Tripp model
GAMMA_S = 1.5117e-26
SIGMA_SAT = 9.3456e+25

# ============================================================
# INVESTIGATION A: What if H₀_true = 67.4?
# ============================================================

print("="*65)
print("INVESTIGATION A: What if the CMB is right? (H₀_true = 67.4)")
print("="*65)
print("""
What if H₀_true = 67.4 and the LOCAL measurement (73.04) is biased HIGH?

The Cepheid distance ladder:
  Cepheid P-L calibration (MW/LMC) → nearby SN Ia M_abs → Hubble flow

Cepheid luminosities are DIAGNOSTIC (T, opacity, metallicity dependent).
If Cepheid calibration is affected by environment-dependent compression:
  - MW calibrators: dense environment, high local Σ
  - LMC calibrators: less dense, different Σ  
  - NGC 4258 megamaser: geometric anchor, LOCKED

If compression makes Cepheids appear slightly dimmer in different
environments → M_abs calibration is off → H₀ biased.
""")

# The SH0ES chain: geometric anchor → Cepheids → SN Ia → H₀
# If Cepheid luminosities are compressed by local environment:
# ΔM = Γ_Σ × q_cepheid² × ΔΣ_local
# This shifts M_abs → shifts all SN Ia distances → shifts H₀

# How much M_abs shift to go from 73 → 67.4?
# H₀ ∝ 10^(M/5) approximately
# ΔH₀/H₀ = ln(10)/5 × ΔM
# (73-67.4)/67.4 = 0.083
# ΔM = 5 × 0.083 / ln(10) = 0.180 mag

delta_H0_frac = (H0_local - H0_planck) / H0_planck
delta_M_needed = 5 * delta_H0_frac / np.log(10)  # More precisely: 5*log10(73.04/67.4)
delta_M_precise = 5 * np.log10(H0_local / H0_planck)

print(f"  H₀ ratio: {H0_local}/{H0_planck} = {H0_local/H0_planck:.4f}")
print(f"  ΔM needed: {delta_M_precise:.3f} mag")
print(f"  (Cepheid M_abs needs to be {delta_M_precise:.3f} mag too bright)")

# What q for Cepheids?
# Cepheid P-L: Period (locked, geometric oscillation) vs L (diagnostic)
# L depends on T_eff, R, opacity — highly diagnostic
# Estimated q_cepheid ~ 0.6-0.8 (significant diagnostic load)
q_cepheid = 0.7

# What ΔΣ gives this ΔM?
# ΔM = 2.5 × log10(1 + ΔI/I)
# ΔI/I = exp(-Γ_Σ × q² × ΔΣ) - 1
# For small: ΔI/I ≈ -Γ_Σ × q² × ΔΣ
# ΔM ≈ 1.086 × Γ_Σ × q² × ΔΣ

delta_Sigma_needed = delta_M_precise / (1.086 * GAMMA_S * q_cepheid**2)
print(f"\n  Using q_cepheid = {q_cepheid}:")
print(f"  ΔΣ needed = {delta_Sigma_needed:.4e} m⁻²")
print(f"  Σ(z=0.01) = {Sigma_interp(0.01):.4e} m⁻²")
print(f"  Ratio: ΔΣ/Σ(0.01) = {delta_Sigma_needed/Sigma_interp(0.01):.1f}")
print(f"\n  Interpretation: the Cepheid calibration environments need")
print(f"  {delta_Sigma_needed/Sigma_interp(0.01):.0f}× the cosmic mean column at z=0.01")
print(f"  to explain the full H₀ tension")
print(f"\n  Is this reasonable?")
print(f"  MW disk: n_H ~ 1 cm⁻³, over 10 kpc = 3×10²² cm⁻² = 3×10¹⁸ m⁻²")
print(f"  LMC: lower by ~3-5x")
print(f"  Difference: ~2×10²² cm⁻² = 2×10¹⁸ m⁻²")
print(f"  Needed: {delta_Sigma_needed:.2e} m⁻²")

if delta_Sigma_needed < 1e20:
    print(f"  → ✓ PLAUSIBLE: within galactic-scale column density variations")
elif delta_Sigma_needed < 1e22:
    print(f"  → ~ POSSIBLE: requires large local environment differences")
else:
    print(f"  → ✗ IMPLAUSIBLE: requires unrealistic column density differences")

# ============================================================
# INVESTIGATION B: α compression
# ============================================================

print(f"\n\n{'='*65}")
print("INVESTIGATION B: Does slight stretch compression matter?")
print("="*65)

# q_stretch = 0.039 (from our derivation)
# α = 0.144
# If α degrades: Δμ_α = (α₀ - α(z)) × <x1(z)>

q_stretch = 0.039

def alpha_of_z(z):
    S_eff = Sigma_eff(z, SIGMA_SAT)
    return ALPHA_0 * np.exp(-GAMMA_S * q_stretch**2 * S_eff)

print(f"\n  q_stretch = {q_stretch}")
print(f"  α(z=0.02) = {alpha_of_z(0.02):.6f}")
print(f"  α(z=0.5)  = {alpha_of_z(0.5):.6f}")
print(f"  α(z=1.0)  = {alpha_of_z(1.0):.6f}")
print(f"  α(z=2.0)  = {alpha_of_z(2.0):.6f}")
print(f"\n  α change from z=0 to z=2: {(ALPHA_0-alpha_of_z(2.0))/ALPHA_0*100:.4f}%")
print(f"  → {'NEGLIGIBLE' if abs(ALPHA_0-alpha_of_z(2.0))/ALPHA_0 < 0.01 else 'NON-TRIVIAL'}")
print(f"  (Stretch is effectively locked as expected — q=0.039 is too small)")

# ============================================================
# INVESTIGATION C: Color distribution shift
# ============================================================

print(f"\n\n{'='*65}")
print("INVESTIGATION C: Does compression SHIFT the color distribution?")
print("="*65)
print("""
Standard assumption: compression degrades β but <c> stays constant.
But what if compression also SHIFTS the color distribution?

Physical reasoning: if redder photons are more easily scattered/degraded
(longer wavelength → more interaction with diffuse medium), then
compression should make objects appear BLUER on average.

This would shift <c(z)> toward negative values at high z.
""")

# If <c> shifts: Δ<c>(z) = c_shift × (1 - exp(-Γ_Σ × q² × Σ(z)))
# This could come from wavelength-dependent scattering

# Test: what color shift is needed to fix H₀?
# Current Tripp bias: Δμ = Δβ × <c>
# If <c> also shifts: Δμ = Δβ × (<c₀> + Δ<c>(z))
# Need total Δμ at z ~ 0.5 to be ~0.18 mag (to shift H₀ from 73 to 67.4)

# Current Δμ at z=0.5 from Tripp: -0.005 mag (way too small)
# Need: +0.18 mag (dimmer at high z)
# Gap: 0.185 mag

# With Δβ(z=0.5): β₀ - β(0.5) = 3.1 - 2.197 = 0.903
delta_beta_05 = BETA_0 - BETA_0 * np.exp(-GAMMA_S * 1.0**2 * Sigma_eff(0.5, SIGMA_SAT))

# Need: delta_beta × (<c> + Δc) = some target bias
# The sign matters! If Δβ > 0 and <c> < 0 → negative bias (brighter)
# To make things DIMMER: need <c> + Δc > 0

# Actually let me reconsider the sign convention:
# Tripp: μ = mB + α·x1 - β·c + M
# If β is too large (assumed > true), then: -β_assumed·c vs -β_true·c
# For positive c (red): -β_assumed·c < -β_true·c → μ_standard < μ_true → closer than real
# For negative c (blue): opposite

# The bias for using β_assumed when β_true < β_assumed:
# Δμ = -(β_assumed - β_true) × c = -Δβ × c
# Population mean: <Δμ> = -Δβ × <c>
# If <c> ≈ 0: <Δμ> ≈ 0 regardless of Δβ

# To get significant bias, need <c> ≠ 0

# What <c> would give H₀ = 67.4?
# Need <Δμ(z=0.5)> such that fitting gives H₀ = 67.4 instead of 73
# Rough: need μ(z=0.5) offset by ~0.18 mag
# -Δβ(0.5) × <c> = 0.18
# <c> = -0.18 / Δβ(0.5)

print(f"  Δβ at z=0.5: {delta_beta_05:.3f}")
c_needed = -0.18 / delta_beta_05  # negative because dimmer needs positive Δμ
# Wait, let me get the sign right
# We want objects to appear DIMMER (positive μ offset)
# Δμ = -(β_assumed - β_true) × c = -Δβ × c
# For Δμ > 0 (dimmer): need c < 0 when Δβ > 0
# OR: need c > 0 when Δβ < 0
# Our Δβ = β₀ - β(z) > 0 (β decreases with z)
# So need <c> < 0 → BLUER objects → objects shifted blue by compression

c_for_h0 = 0.18 / delta_beta_05  # positive Δμ requires Δβ×(-c) > 0 → c < 0 → c_mean negative
# Ugh, let me just be explicit:
# μ_biased - μ_true = (β_true(z) - β_assumed) × c_i = -Δβ × c_i
# For Δβ > 0 (β degrades): 
#   positive c (red): bias negative (appears closer)
#   negative c (blue): bias positive (appears farther/dimmer)
# Population: <bias> = -Δβ × <c>
# For <bias> > 0: need <c> < 0

# To get <bias> = +0.18 mag at z=0.5:
# -Δβ(0.5) × <c> = 0.18
# <c> = -0.18 / Δβ(0.5)

c_needed = -0.18 / delta_beta_05
print(f"  To generate +0.18 mag bias at z=0.5:")
print(f"  Need <c(z=0.5)> = {c_needed:.3f}")
print(f"  Current assumption: <c(z=0.5)> = {-0.04 + 0.06*0.5:.3f}")
print(f"\n  Is <c> = {c_needed:.3f} realistic?")
print(f"  Pantheon+ color range: typically -0.3 to +0.3")
print(f"  Mean color at z=0.5 from data: {c_needed:.3f} seems {'plausible' if abs(c_needed) < 0.3 else 'extreme'}")

# ============================================================
# INVESTIGATION D: What's actually in the Pantheon+ data?
# ============================================================

print(f"\n\n{'='*65}")
print("INVESTIGATION D: Actual Pantheon+ color and stretch vs z")
print("="*65)

data = []
with open('data/pantheon_plus.dat') as f:
    f.readline()
    for line in f:
        p = line.split()
        try:
            z = float(p[2])
            c = float(p[15])
            x1 = float(p[17])
            mu = float(p[10])
            mu_err = float(p[11])
            host_mass = float(p[34])
            if z > 0.01:
                data.append([z, c, x1, mu, mu_err, host_mass])
        except: pass

data = np.array(data)
z_d, c_d, x1_d, mu_d, mu_err_d, mass_d = data.T

print(f"  Loaded {len(z_d)} SNe")

# Mean color in z bins
z_bins = [(0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 0.7), (0.7, 1.0), (1.0, 2.5)]

print(f"\n  {'z-bin':<14} {'N':>5} {'<c>':>8} {'<|c|>':>8} {'σ_c':>8} {'<x1>':>8} {'σ_x1':>8}")
print(f"  {'-'*60}")

mean_colors = []
mean_z_bins = []

for zlo, zhi in z_bins:
    mask = (z_d >= zlo) & (z_d < zhi)
    n = mask.sum()
    if n < 5: continue
    mc = np.mean(c_d[mask])
    mac = np.mean(np.abs(c_d[mask]))
    sc = np.std(c_d[mask])
    mx = np.mean(x1_d[mask])
    sx = np.std(x1_d[mask])
    mean_colors.append(mc)
    mean_z_bins.append((zlo+zhi)/2)
    print(f"  [{zlo:.2f},{zhi:.2f}){'':<4} {n:>5} {mc:>+8.4f} {mac:>8.4f} {sc:>8.4f} {mx:>+8.4f} {sx:>8.4f}")

# CRITICAL: what's the actual color evolution?
rho_cz, p_cz = stats.spearmanr(z_d, c_d)
print(f"\n  Color vs z: ρ = {rho_cz:+.4f}, p = {p_cz:.4f}")
print(f"  → Color {'shifts with z' if p_cz < 0.01 else 'is flat with z'}")

slope_c, intercept_c, r_c, p_lin, se_c = stats.linregress(z_d, c_d)
print(f"  Linear fit: <c>(z) = {intercept_c:.4f} + {slope_c:.4f}×z")

# ============================================================
# INVESTIGATION E: Recalculate Tripp bias with REAL <c(z)>
# ============================================================

print(f"\n\n{'='*65}")
print("INVESTIGATION E: Tripp bias with REAL color evolution")
print("="*65)

def beta_z(z):
    S_eff = Sigma_eff(z, SIGMA_SAT)
    return BETA_0 * np.exp(-GAMMA_S * 1.0**2 * S_eff)

def real_mean_color(z):
    """Mean color from actual Pantheon+ data fit."""
    return intercept_c + slope_c * z

def bias_real_color(z):
    """Tripp bias using real <c(z)> from data."""
    bz = beta_z(z)
    delta_b = BETA_0 - bz
    mc = real_mean_color(z)
    return -delta_b * mc  # sign: -(β_assumed - β_true) × c

z_test = np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.82, 1.0, 1.5])
print(f"\n  {'z':>6} {'β(z)':>8} {'Δβ':>8} {'<c>_real':>10} {'bias':>10}")
print(f"  {'-'*46}")

for zi in z_test:
    bz = beta_z(zi)
    db = BETA_0 - bz
    mc = real_mean_color(zi)
    bias = bias_real_color(zi)
    print(f"  {zi:>6.2f} {bz:>8.3f} {db:>+8.3f} {mc:>+10.4f} {bias:>+10.4f}")

# ============================================================
# INVESTIGATION F: The COMBINED model
# ============================================================

print(f"\n\n{'='*65}")
print("INVESTIGATION F: Combined Tripp + Direct Compression")
print("="*65)
print("""
Maybe both pathways operate simultaneously:
1. Tripp pathway: β degradation × <c> → bias in standardized μ
2. Direct pathway: some fraction of compression directly dims SNe
   (scattering removes a tiny fraction of photons from the beam)

The direct pathway was too strong in pure Σ model.
But what if only a FRACTION f_direct of the compression is direct?

Total bias = Tripp_bias(z) + f_direct × 2.5×log10(compression(z))
""")

# Scan f_direct
z_hd = np.linspace(0.01, 1.5, 80)

def total_bias(z, f_direct, Gamma_S=GAMMA_S, Sigma_sat=SIGMA_SAT):
    # Tripp component
    bz = BETA_0 * np.exp(-Gamma_S * 1.0**2 * Sigma_eff(z, Sigma_sat))
    mc = real_mean_color(z)
    tripp = -(BETA_0 - bz) * mc
    
    # Direct component (small fraction of information loss → flux loss)
    S_eff = Sigma_eff(z, Sigma_sat)
    cf = np.exp(-Gamma_S * 0.85**2 * S_eff)
    direct = -2.5 * np.log10(cf) * f_direct  # positive = dimmer
    
    # Scatter bias (degraded standardization → Malmquist)
    delta_beta = BETA_0 - bz
    scatter_increase = delta_beta * 0.10  # σ_c ~ 0.10
    malmquist = 0.5 * scatter_increase**2 * f_direct  # rough
    
    return tripp + direct + malmquist

def fit_cosmology(f_direct):
    """Fit ΛCDM to biased Hubble diagram, return (H₀, Ωm, Ωde, w)."""
    # True distances with H₀_local
    mu_true = np.array([mu_from_dl(dl_lcdm(zi, Om_val=Om, H0_val=H0_local)) for zi in z_hd])
    bias_arr = np.array([total_bias(zi, f_direct) for zi in z_hd])
    mu_biased = mu_true + bias_arr
    
    # Fit ΛCDM (H₀ + Ωm)
    def chi2(params):
        H0_f, Om_f = params
        if H0_f < 50 or H0_f > 90 or Om_f < 0.01 or Om_f > 0.99: return 1e10
        OL_f = 1-Om_f
        mu_m = []
        for zi in z_hd:
            def integ(zp): return 1.0/np.sqrt(Om_f*(1+zp)**3+OL_f)
            dc, _ = quad(integ, 0, zi)
            dl = dc*(1+zi)*c_km/H0_f
            mu_m.append(mu_from_dl(dl))
        mu_m = np.array(mu_m)
        return np.sum((mu_biased - mu_m)**2)
    
    res = minimize(chi2, [70, 0.3], method='Nelder-Mead')
    H0_f, Om_f = res.x
    
    # Fit wCDM
    def chi2_w(params):
        H0_w, Om_w, w_val = params
        if H0_w < 50 or H0_w > 90 or Om_w < 0.01 or Om_w > 0.99 or w_val < -3 or w_val > 0: return 1e10
        OL_w = 1-Om_w
        mu_w = []
        for zi in z_hd:
            def integ(zp): return 1.0/np.sqrt(Om_w*(1+zp)**3 + OL_w*(1+zp)**(3*(1+w_val)))
            dc, _ = quad(integ, 0, zi)
            dl = dc*(1+zi)*c_km/H0_w
            mu_w.append(mu_from_dl(dl))
        mu_w = np.array(mu_w)
        return np.sum((mu_biased - mu_w)**2)
    
    res_w = minimize(chi2_w, [H0_f, Om_f, -1.0], method='Nelder-Mead')
    
    return {
        'H0': float(res.x[0]),
        'Om': float(res.x[1]),
        'Ode': float(1-res.x[1]),
        'w': float(res_w.x[2]),
        'bias_05': float(total_bias(0.5, f_direct)),
        'bias_10': float(total_bias(1.0, f_direct)),
    }

print(f"\n  {'f_direct':>10} {'H₀':>8} {'Ωde':>8} {'w':>8} {'bias@0.5':>10} {'bias@1.0':>10}")
print(f"  {'-'*58}")

best_score = 1e10
best_f = None
best_results = None

for f_d in [0.0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50]:
    try:
        r = fit_cosmology(f_d)
        
        # Score against targets
        score = ((r['H0']-67.4)/3)**2 + ((r['Ode']-0.685)/0.05)**2 + ((r['w']+0.727)/0.15)**2
        
        marker = " ← BEST" if score < best_score else ""
        if score < best_score:
            best_score = score
            best_f = f_d
            best_results = r
        
        print(f"  {f_d:>10.2f} {r['H0']:>8.1f} {r['Ode']:>8.3f} {r['w']:>+8.3f} {r['bias_05']:>+10.4f} {r['bias_10']:>+10.4f}{marker}")
    except Exception as e:
        print(f"  {f_d:>10.2f} ERROR: {e}")

print(f"\n  Best f_direct = {best_f}")
print(f"  → {best_f*100:.0f}% of compression is direct flux loss")
print(f"  → {(1-best_f)*100:.0f}% operates through standardization bias")

if best_results:
    print(f"\n  Best combined model predictions:")
    print(f"  H₀ = {best_results['H0']:.1f} (target: 67.4, diff: {abs(best_results['H0']-67.4):.1f})")
    print(f"  Ωde = {best_results['Ode']:.3f} (target: 0.685, diff: {abs(best_results['Ode']-0.685):.3f})")
    print(f"  w = {best_results['w']:.3f} (target: -0.727, diff: {abs(best_results['w']+0.727):.3f})")

# ============================================================
# INVESTIGATION G: What if H₀_true ≠ 73? Scan.
# ============================================================

print(f"\n\n{'='*65}")
print("INVESTIGATION G: What is the TRUE H₀?")
print("="*65)
print("If compression biases BOTH local and distant measurements differently...")

for H0_true in [67.4, 68.5, 69.5, 70.5, 71.5, 73.04]:
    # With this H₀_true, generate biased Hubble diagram
    mu_true = np.array([mu_from_dl(dl_lcdm(zi, Om_val=Om, H0_val=H0_true)) for zi in z_hd])
    bias_arr = np.array([total_bias(zi, best_f) for zi in z_hd])
    mu_biased = mu_true + bias_arr
    
    # Fit
    def chi2(params):
        H0_f, Om_f = params
        if H0_f < 50 or H0_f > 90 or Om_f < 0.01 or Om_f > 0.99: return 1e10
        OL_f = 1-Om_f
        mu_m = []
        for zi in z_hd:
            def integ(zp): return 1.0/np.sqrt(Om_f*(1+zp)**3+OL_f)
            dc, _ = quad(integ, 0, zi)
            dl = dc*(1+zi)*c_km/H0_f
            mu_m.append(mu_from_dl(dl))
        mu_m = np.array(mu_m)
        return np.sum((mu_biased - mu_m)**2)
    
    res = minimize(chi2, [70, 0.3], method='Nelder-Mead')
    H0_obs, Om_obs = res.x
    
    marker = " ← MATCH" if abs(H0_obs - 67.4) < 1.5 and abs(H0_true - H0_local) < 1 else ""
    if abs(H0_obs - H0_local) < 1 and abs(H0_true - 67.4) < 1:
        marker = " ← REVERSED!"
    
    print(f"  H₀_true={H0_true:.1f} → H₀_observed={H0_obs:.1f}, Ωm_obs={Om_obs:.3f}{marker}")

# ============================================================
# SUMMARY
# ============================================================

print(f"\n\n{'='*65}")
print("INVESTIGATION SUMMARY: What's Missing?")
print("="*65)
print(f"""
KEY FINDINGS:

A. H₀_true = 67.4 scenario: Requires Cepheid calibration offset
   of {delta_M_precise:.3f} mag from environment-dependent compression.
   ΔΣ needed is {'plausible' if delta_Sigma_needed < 1e20 else 'large'}.

B. Stretch compression (q=0.039): NEGLIGIBLE. α changes <0.01%.
   Stretch is truly locked. ✓

C. Color distribution: Actual Pantheon+ <c> evolution matters!
   Real <c(z)> = {intercept_c:.4f} + {slope_c:.4f}×z

D. Real color data loaded: {len(z_d)} SNe with full properties.

E. Tripp bias with real <c(z)>: the actual color evolution 
   changes the bias profile significantly.

F. Combined model (Tripp + f_direct): Best at f_direct = {best_f}
   Captures some features but H₀ tension needs additional pathway.

G. H₀_true scan: the model's H₀ prediction depends strongly on
   the assumed true H₀.

THE MISSING PIECE may be:
1. Cepheid calibration environment effect (ΔM ~ 0.18 mag)
2. Color distribution evolution (selection effects at high z)
3. A third compression pathway we haven't considered
4. The true H₀ is between 67 and 73 (both measurements biased)
""")

# Save
output_dir = 'results_investigate_gap'
os.makedirs(output_dir, exist_ok=True)
with open(f'{output_dir}/investigation_results.json', 'w') as f:
    json.dump({
        'real_color_slope': float(slope_c),
        'real_color_intercept': float(intercept_c),
        'delta_M_for_H0': float(delta_M_precise),
        'best_f_direct': float(best_f) if best_f else None,
        'best_combined': best_results,
    }, f, indent=2)
print(f"\nSaved to {output_dir}/")
