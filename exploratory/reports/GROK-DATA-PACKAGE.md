# Data Package for Grok — Evolving Coefficient Section

## Full Modified Tripp Coefficients (from closure_eigenvector_rotation.py, N=530 subsample)

```
μ = m_B + α(z)×x1 − β(z)×c − δ×c² − γ(z)×δ_M − M_B

α(z) = 0.1484 + (−0.0321) × z
β(z) = 3.108 + (−0.615) × z  
γ(z) = −0.0277 + (−0.0701) × z
δ(c²) = +3.638
M_B = −19.253 (reference)

Δχ² = 600.9 vs standard Tripp (5 extra parameters, p ≈ 0)
```

## Coefficient Values at Key Redshifts

| z | α(z) | β(z) | γ(z) |
|------|--------|--------|---------|
| 0.01 | 0.1481 | 3.102 | −0.0284 |
| 0.05 | 0.1468 | 3.078 | −0.0312 |
| 0.10 | 0.1452 | 3.047 | −0.0348 |
| 0.30 | 0.1388 | 2.924 | −0.0488 |
| 0.50 | 0.1324 | 2.801 | −0.0628 |
| 1.00 | 0.1163 | 2.493 | −0.0978 |

## α(z) Standalone Fit (from closure_beta_z_reanalysis.py, N=1590)

```
α(z) = 0.1592 − 0.0668 × z
H₀ = 68.79 km/s/Mpc (tension = 1.43)
Δχ² = 178.13 vs constant α
```

## β(z) Binned Values (from closure_beta_z_reanalysis.py)

| z-bin | N | β | α |
|------------|------|-------|-------|
| [0.01,0.05]| 176 | 3.049 | 0.151 |
| [0.05,0.10]| 31 | 3.161 | 0.125 |
| [0.10,0.25]| 116 | 3.183 | 0.129 |
| [0.25,0.50]| 134 | 2.825 | 0.111 |
| [0.50,1.00]| 62 | 1.968 | 0.098 |

## Mass Step — Binned by Redshift (from closure_mass_step_impedance.py)

| z-bin | N_hi | N_lo | ⟨resid⟩_hi | ⟨resid⟩_lo | Step | p-value |
|------------|------|------|------------|------------|---------|---------|
| [0.01,0.05]| 97 | 79 | −0.0053 | +0.0554 | −0.0606 | 0.0262 |
| [0.05,0.10]| 16 | 15 | −0.0115 | −0.0057 | −0.0058 | 0.8986 |
| [0.10,0.25]| 60 | 56 | −0.0668 | +0.0331 | −0.0999 | 0.0007 |
| [0.25,0.50]| 75 | 59 | −0.0569 | +0.0012 | −0.0581 | 0.0302 |
| [0.50,1.00]| 26 | 36 | −0.0622 | −0.0699 | +0.0077 | 0.8796 |

Mass step evolution with z: Spearman ρ = +0.500

## Mass Step Decomposition

| Standardization | Mass step (γ) |
|----------------|--------------|
| Full Tripp | −0.0461 mag |
| Stretch-only | −0.0135 mag |
| Color-driven | 0.0326 mag (71%) |

## Quadratic Color Term (c²)

```
Standard Tripp: χ² = 20385.8
Tripp + c²: χ² = 20047.1 (Δχ² = 338.7, p ≈ 0)
δ = +2.045

Tripp + c² + c×z: χ² = 19992.1 (Δχ² = 393.7)
δ = +1.805, η(c×z) = −0.385
```

## w₀wₐ Results (from closure_w0wa_refit.py, N=795)

| Model | w₀ | wₐ | χ² | k |
|-------|------|-------|--------|---|
| Standard Tripp + ΛCDM | −1.0 | 0.0 | 10126.7 | 4 |
| Standard Tripp + w₀wₐ | −0.559 | −1.883 | 9811.4 | 6 |
| α(z) + w₀wₐ | −0.614 | −1.387 | 9623.0 | 7 |
| γ(z) + w₀wₐ | −0.538 | −1.754 | 9792.7 | 7 |
| Full Modified + w₀wₐ | −0.659 | −0.020 | 9277.1 | 10 |
| Full Modified + ΛCDM | −1.0 | 0.0 | 9774.7 | 8 |

### Key Result
- wₐ: −1.883 → −0.020 (99% absorbed by standardization)
- w₀: −0.559 → −0.659 (23% absorbed, shifts toward ΛCDM)
- Full Modified + ΛCDM BEATS Standard + w₀wₐ by Δχ² = 37

## Impedance Framework Constants

```
D_i(z) = a × N_modes_i^α × Z_g(z)

α = 1.845 (cooperative exponent)
a = 0.0204 (eating law amplitude)
Γ₀ = 2.17 (coupling constant)
β_sightline = 1.43 (Z_g power law)
```

## N_modes Assignments

| Observable | N_modes | Basis |
|-----------|---------|-------|
| CMB/BAO | 0 | geometric |
| GW sirens | 0 | geometric |
| TRGB | 2 | Z, age |
| SN stretch | ~1 | Ni56 mass |
| SN color | ~3 | T, Z, extinction |
| SN Tripp | 4 | stretch + color |
| Cepheid | 4-5 | T_eff, Z, extinction, pulsation |

## Data Source
- Pantheon+ SN Ia sample: 1,590 SNe after quality cuts
- File: `data/pantheon_plus.dat` (1,701 entries, 47 columns)
- SALT2 parameters: x1, c, mB, HOST_LOGMASS, zCMB, zHD
