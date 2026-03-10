#!/usr/bin/env python3
"""
THE Σ FRAMEWORK: Column Density as the Compression Variable

Key hypothesis: the z-sigmoid (z₀=0.82, k=8.0) was an empirical 
approximation. The REAL compression variable is Σ(z) — the integrated
matter column density along the sightline.

If this is correct:
1. Σ(z) should be naturally sigmoid-shaped (matter density evolution)
2. The inflection point of Σ(z) should fall near z ≈ 0.82
3. Using Σ directly should give STRONGER compression at low z
4. z₀ becomes a DERIVED quantity, not a free parameter

The compression law becomes:
   I(z) = I₀ · exp(-Γ_Σ · q² · Σ(z))

where Σ(z) = ∫₀ᶻ n_baryon(z') · c · dt/dz' · dz'
            = ∫₀ᶻ n_b0·(1+z')³ · c/(H(z')·(1+z')) · dz'

Author: Closure Theory Pipeline
Date: March 2026
"""

import numpy as np
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import minimize, minimize_scalar, brentq
import json
import os

np.random.seed(42)

# ============================================================
# PHYSICAL CONSTANTS
# ============================================================

c_km = 299792.458        # km/s
H0 = 67.4                # km/s/Mpc (Planck)
H0_local = 73.04         # km/s/Mpc (SH0ES)
Om = 0.315               # Ωm
Ob = 0.0493              # Ωb
OL = 1 - Om              # ΩΛ
rho_crit = 9.47e-27      # kg/m³
m_p = 1.673e-27          # kg (proton mass)
Mpc_to_cm = 3.086e24     # cm per Mpc
c_cm = 2.998e10          # cm/s

# Mean baryon number density today
n_b0 = rho_crit * Ob / m_p  # baryons per m³
print(f"Mean baryon density today: n_b0 = {n_b0:.4f} m⁻³ = {n_b0*1e-6:.4f} cm⁻³")

# ============================================================
# Σ(z): INTEGRATED BARYON COLUMN DENSITY
# ============================================================

def H(z):
    """Hubble parameter at z (km/s/Mpc)."""
    return H0 * np.sqrt(Om * (1+z)**3 + OL)

def dSigma_dz(z):
    """
    Column density integrand: n_b(z) · c · |dt/dz|
    
    n_b(z) = n_b0 · (1+z)³  (baryon density scales with volume)
    |dt/dz| = 1 / (H(z) · (1+z))
    
    Σ(z) = ∫₀ᶻ n_b0·(1+z')³ · c / (H(z')·(1+z')) dz'
         = n_b0 · c · ∫₀ᶻ (1+z')² / H(z') dz'
    
    Returns in units of [baryons / m²] = [m⁻²]
    """
    # c in m/s, H in km/s/Mpc → convert H to s⁻¹
    H_si = H(z) * 1e3 / (3.086e22)  # km/s/Mpc → s⁻¹
    c_ms = 2.998e8  # m/s
    return n_b0 * c_ms * (1+z)**2 / H_si

def Sigma(z):
    """Integrated baryon column density from 0 to z [m⁻²]."""
    if z <= 0:
        return 0.0
    result, _ = quad(dSigma_dz, 0, z)
    return result

# Compute Σ(z) curve
z_arr = np.linspace(0, 3, 200)
sigma_arr = np.array([Sigma(zi) for zi in z_arr])

# Normalize to Σ(z=2) for convenient units
sigma_norm = sigma_arr / Sigma(2.0)

print(f"\nΣ(z) values (normalized to Σ(z=2)):")
for zi in [0.1, 0.3, 0.5, 0.82, 1.0, 1.5, 2.0, 3.0]:
    si = Sigma(zi) / Sigma(2.0)
    print(f"  z = {zi:.2f}: Σ_norm = {si:.4f}")

# ============================================================
# IS Σ(z) SIGMOID-SHAPED?
# ============================================================

print(f"\n{'='*65}")
print("QUESTION: Is Σ(z) naturally sigmoid-shaped?")
print("="*65)

# Compute dΣ/dz to find the inflection point
dsigma = np.gradient(sigma_norm, z_arr)
d2sigma = np.gradient(dsigma, z_arr)

# Inflection point: where d²Σ/dz² changes sign
sign_changes = np.where(np.diff(np.sign(d2sigma)))[0]
if len(sign_changes) > 0:
    z_inflection = z_arr[sign_changes[0]]
    print(f"\n  Σ(z) inflection point: z = {z_inflection:.3f}")
    print(f"  Our empirical z₀: 0.82")
    print(f"  Match: {'✓ YES' if abs(z_inflection - 0.82) < 0.2 else '✗ NO'}")
    print(f"  Discrepancy: {abs(z_inflection - 0.82):.3f}")
else:
    print(f"  No inflection point found — Σ(z) is monotonically accelerating")
    # Actually check: dΣ/dz should decrease because H(z) grows
    # while (1+z)² grows. The competition determines shape.
    print(f"  dΣ/dz at z=0: {dsigma[1]:.4f}")
    print(f"  dΣ/dz at z=1: {dsigma[100]:.4f}")
    print(f"  dΣ/dz at z=2: {dsigma[133]:.4f}")
    z_inflection = None

# Check shape by fitting a sigmoid
def sigmoid(z, z0, k, A):
    return A / (1 + np.exp(-k * (z - z0)))

from scipy.optimize import curve_fit
try:
    popt, pcov = curve_fit(sigmoid, z_arr, sigma_norm, p0=[1.0, 2.0, 1.0], maxfev=10000)
    z0_fit, k_fit, A_fit = popt
    sigma_sigmoid_fit = sigmoid(z_arr, *popt)
    residual_rms = np.sqrt(np.mean((sigma_norm - sigma_sigmoid_fit)**2))
    
    print(f"\n  Sigmoid fit to Σ(z):")
    print(f"  z₀ = {z0_fit:.3f}  (our empirical: 0.82)")
    print(f"  k  = {k_fit:.3f}   (our empirical: 8.0)")
    print(f"  A  = {A_fit:.3f}")
    print(f"  RMS residual: {residual_rms:.6f}")
    
    # The KEY question: does the physics give z₀ ≈ 0.82?
    print(f"\n  Does physics derive z₀ ≈ 0.82?")
    if abs(z0_fit - 0.82) < 0.3:
        print(f"  → ✓ YES! z₀ = {z0_fit:.3f} from baryon density evolution alone")
        print(f"  → z₀ is NOT a free parameter — it's derivable from Ωb and H(z)")
    else:
        print(f"  → ✗ No, z₀ = {z0_fit:.3f} differs from 0.82")
except Exception as e:
    print(f"  Sigmoid fit failed: {e}")
    z0_fit = None

# ============================================================
# Σ-BASED COMPRESSION LAW
# ============================================================

print(f"\n\n{'='*65}")
print("Σ-BASED COMPRESSION LAW")
print("="*65)

# The law: I(z) = I₀ · exp(-Γ_Σ · q² · Σ(z))
# Where Σ(z) is the physical baryon column density
# Γ_Σ has units of [m²] — an effective cross-section!

# From our empirical Γ₀ = 0.533 (dimensionless, per unit sigmoid):
# Γ_Σ = Γ₀ / Σ_scale
# where Σ_scale normalizes the sigmoid ∫σ(z')dz' to physical Σ

# Let's work backwards: what Γ_Σ reproduces our observed β degradation?
# β drops from 2.94 to 1.64 between z=0.02 and z=0.9
# That's a factor of 0.56, meaning compression ≈ 0.44 information loss

# I(z=0.9)/I(z=0.02) = exp(-Γ_Σ · q² · [Σ(0.9) - Σ(0.02)])
# 0.56 = exp(-Γ_Σ · 1.0² · ΔΣ)

delta_Sigma_phys = Sigma(0.9) - Sigma(0.02)
ln_ratio = np.log(0.56)  # = -0.58
Gamma_Sigma = -ln_ratio / delta_Sigma_phys
print(f"\n  β degradation: 2.94 → 1.64 (ratio = 0.56)")
print(f"  Σ(0.9) - Σ(0.02) = {delta_Sigma_phys:.4e} m⁻²")
print(f"  Γ_Σ = {Gamma_Sigma:.4e} m²")
print(f"\n  This is an EFFECTIVE CROSS-SECTION for information loss!")
print(f"  Thomson cross-section: σ_T = 6.65 × 10⁻²⁹ m²")
print(f"  Ratio Γ_Σ / σ_T: {Gamma_Sigma / 6.652e-29:.2e}")

# If Γ_Σ is much smaller than σ_T, compression isn't from scattering
# If comparable, it could be related to Thomson scattering

# ============================================================
# REDERIVE COSMOLOGICAL MYSTERIES WITH Σ
# ============================================================

Q_SN_EFF = 0.85

def compression_Sigma(z, q=Q_SN_EFF):
    """Information survival fraction using physical Σ."""
    S = Sigma(z)
    return np.exp(-Gamma_Sigma * q**2 * S)

def distance_inflation_Sigma(z, q=Q_SN_EFF):
    """Distance inflation from Σ-based compression."""
    cf = compression_Sigma(z, q)
    # Compression reduces diagnostic flux → SN appears dimmer → distance overestimated
    # For SN Ia: dimmer by factor cf → distance inflated by 1/sqrt(cf)
    return 1.0 / np.sqrt(cf)

print(f"\n\n{'='*65}")
print("REDERIVATION 1: H₀ TENSION (Σ-based)")
print("="*65)

# Distance inflation at various z
print(f"\n  {'z':>6} {'Σ(z)/Σ(2)':>10} {'compression':>12} {'d_inflation':>12}")
print(f"  {'-'*44}")
for zi in [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.82, 1.0, 1.5, 2.0]:
    cf = compression_Sigma(zi)
    di = distance_inflation_Sigma(zi)
    sn = Sigma(zi)/Sigma(2.0)
    print(f"  {zi:>6.2f} {sn:>10.4f} {cf:>12.6f} {di:>12.4f}x")

# H₀ tension: at z ~ 0.5 (median Pantheon+ for H₀ analysis)
z_pivot = 0.5
infl = distance_inflation_Sigma(z_pivot)
h0_biased = H0_local / infl

print(f"\n  At z = {z_pivot}:")
print(f"  Distance inflation: {infl:.4f}x")
print(f"  H₀_true = {H0_local} → H₀_apparent = {h0_biased:.1f}")
print(f"  Observed H₀_CMB = {H0}")
print(f"  Match: {'✓' if abs(h0_biased - H0) < 3 else '✗'} (diff = {abs(h0_biased - H0):.1f})")

# What z_pivot gives exact match?
try:
    z_exact = brentq(lambda z: H0_local/distance_inflation_Sigma(z) - H0, 0.01, 10.0)
    print(f"\n  z where H₀_apparent exactly = {H0}: z = {z_exact:.3f}")
except:
    print(f"\n  Could not find exact match z")

# ============================================================
# REDERIVATION 2: APPARENT Ωde (Σ-based)
# ============================================================

print(f"\n\n{'='*65}")
print("REDERIVATION 2: APPARENT Ωde (Σ-based)")
print("="*65)

def dl_model(z_val, Om_true, H0_val=73.04):
    OL_val = 1 - Om_true
    def integrand(zp):
        return 1.0 / np.sqrt(Om_true*(1+zp)**3 + OL_val)
    dc, _ = quad(integrand, 0, z_val)
    return dc * (1+z_val) * c_km / H0_val

def mu_from_dl(dl):
    return 5*np.log10(dl) + 25 if dl > 0 else 0

z_mock = np.linspace(0.01, 1.5, 100)

# For each true Ωm, compute compressed Hubble diagram and fit apparent Ωm
def get_apparent_Ode(Om_true):
    # True distances
    mu_true = np.array([mu_from_dl(dl_model(zi, Om_true)) for zi in z_mock])
    # Add Σ-based compression
    inflation = np.array([distance_inflation_Sigma(zi) for zi in z_mock])
    mu_compressed = mu_true + 5 * np.log10(inflation)
    
    # Fit ΛCDM
    def chi2(params):
        Om_fit = params[0]
        if Om_fit < 0.01 or Om_fit > 0.99: return 1e10
        OL_fit = 1 - Om_fit
        mu_model = []
        for zi in z_mock:
            def integ(zp):
                return 1.0/np.sqrt(Om_fit*(1+zp)**3 + OL_fit)
            dc, _ = quad(integ, 0, zi)
            dl = dc*(1+zi)*c_km/73.04
            mu_model.append(mu_from_dl(dl))
        mu_model = np.array(mu_model)
        offset = np.mean(mu_compressed - mu_model)
        return np.sum((mu_compressed - mu_model - offset)**2)
    
    res = minimize(chi2, [0.5], method='Nelder-Mead')
    return 1 - res.x[0]

print(f"\n  {'True Ωm':>10} {'True Ωde':>10} {'Apparent Ωde':>14} {'Match?':>8}")
print(f"  {'-'*46}")

for Om_t in np.arange(0.30, 0.85, 0.05):
    Ode_app = get_apparent_Ode(Om_t)
    match = " ← MATCH" if abs(Ode_app - 0.685) < 0.03 else ""
    print(f"  {Om_t:>10.2f} {1-Om_t:>10.2f} {Ode_app:>14.3f} {match:>8}")

# ============================================================
# REDERIVATION 3: w EVOLUTION (Σ-based)
# ============================================================

print(f"\n\n{'='*65}")
print("REDERIVATION 3: w EVOLUTION (Σ-based)")
print("="*65)

def dl_wcdm(z_val, Om=0.315, w=-1.0, H0_val=73.04):
    OL_val = 1 - Om
    def integrand(zp):
        return 1.0/np.sqrt(Om*(1+zp)**3 + OL_val*(1+zp)**(3*(1+w)))
    dc, _ = quad(integrand, 0, z_val)
    return dc*(1+z_val)*c_km/H0_val

# True ΛCDM distances + compression → fit w at each z
z_w = np.linspace(0.1, 1.5, 15)
print(f"\n  {'z':>6} {'inflation':>10} {'w_eff':>8}")
print(f"  {'-'*28}")

w_eff_arr = []
for zi in z_w:
    dl_t = dl_model(zi, 0.315, 73.04)
    infl_i = distance_inflation_Sigma(zi)
    dl_c = dl_t * infl_i
    
    def w_obj(w):
        try:
            dl_w = dl_wcdm(zi, w=w)
            return (dl_w - dl_c)**2
        except: return 1e10
    
    res = minimize_scalar(w_obj, bounds=(-3, 0), method='bounded')
    w_eff = res.x
    w_eff_arr.append(w_eff)
    print(f"  {zi:>6.2f} {infl_i:>10.4f} {w_eff:>+8.3f}")

w_eff_arr = np.array(w_eff_arr)
desi_w0 = -0.727
desi_mask = (z_w >= 0.3) & (z_w <= 0.8)
our_w = np.mean(w_eff_arr[desi_mask])
print(f"\n  Mean w_eff (z=0.3-0.8): {our_w:.3f}")
print(f"  DESI w₀: {desi_w0}")
print(f"  Match: {'✓' if abs(our_w - desi_w0) < 0.15 else '✗'}")

# ============================================================
# THE CROSS-SECTION
# ============================================================

print(f"\n\n{'='*65}")
print("THE INFORMATION CROSS-SECTION")
print("="*65)

sigma_T = 6.652e-29  # Thomson
sigma_info = Gamma_Sigma

print(f"""
  Γ_Σ = {sigma_info:.4e} m²
  σ_Thomson = {sigma_T:.4e} m²
  
  Ratio: Γ_Σ / σ_T = {sigma_info/sigma_T:.4e}
""")

if sigma_info < sigma_T:
    print(f"  Γ_Σ << σ_Thomson")
    print(f"  → Compression is NOT standard Thomson scattering")
    print(f"  → It's a MUCH weaker interaction")
    print(f"  → Consistent with: information degradation through")
    print(f"     cumulative weak interactions, not direct scattering")
elif sigma_info < sigma_T * 100:
    print(f"  Γ_Σ ~ σ_Thomson (within 2 orders)")
    print(f"  → Compression could be related to electromagnetic interaction")
    print(f"     with intervening plasma (IGM, CGM)")
else:
    print(f"  Γ_Σ >> σ_Thomson")
    print(f"  → Compression is NOT simple scattering")
    print(f"  → May involve collective/resonant effects")

# ============================================================
# SUMMARY
# ============================================================

print(f"\n\n{'='*65}")
print("SUMMARY: Σ FRAMEWORK RESULTS")
print("="*65)

print(f"""
The compression law reparameterized with physical Σ(z):

   I(z) = I₀ · exp(-Γ_Σ · q² · Σ(z))
   
   Γ_Σ = {Gamma_Sigma:.4e} m² (effective information cross-section)
   Σ(z) = ∫₀ᶻ n_b(z') · c / H(z') · (1+z')² dz'

Key findings:
1. Σ(z) shape and whether z₀ emerges naturally
2. Distance inflation profile (linear in Σ vs sigmoid in z)
3. H₀ tension prediction
4. Apparent Ωde prediction  
5. w evolution prediction
6. The information cross-section value and its physics
""")

# Save
output_dir = 'results_sigma_framework'
os.makedirs(output_dir, exist_ok=True)
results = {
    'Gamma_Sigma_m2': float(Gamma_Sigma),
    'Gamma_Sigma_over_Thomson': float(Gamma_Sigma/sigma_T),
    'Sigma_z082_over_z2': float(Sigma(0.82)/Sigma(2.0)),
    'z0_from_sigmoid_fit': float(z0_fit) if z0_fit else None,
}
with open(f'{output_dir}/sigma_framework_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"Results saved to {output_dir}/")
