# Pre-Registered Predictions — Environment-Dependent Closure
# Timestamped: 2026-02-22
# These predictions are committed to git BEFORE running the tests.

## Context
We found that CIV/MgII flux Spearman correlation is higher in dense vs sparse 
environments (6/6 z-bins, surviving SNR+lum matching, forward model, shuffle null).
Mean Δρ ≈ +0.01 (dense - sparse).

## Predictions for NEW line pairs (not yet tested with environment split):

### Prediction 1: HeII(4686) / CIV(1549)
- **Coupling**: Both high-ionization, tight physical coupling
- **Expected environment Δρ**: SMALLER than CIV/MgII (~0.005 or less)
- **Reasoning**: Tightly coupled pairs have less room to degrade differentially.
  If closure degrades loosely-coupled pairs more than tightly-coupled ones,
  the environment modulation should also be weaker for tight pairs.

### Prediction 2: MgII(2796) / Hβ_broad
- **Coupling**: Cross-species (low-ionization Mg vs Balmer), moderate coupling
- **Expected environment Δρ**: LARGER than CIV/MgII (~0.015 or more)
- **Reasoning**: More loosely coupled → more susceptible to environment modulation

### Prediction 3: CIV(1549) / NeV(3426)
- **Coupling**: Very different ionization zones (mid vs very high ionization)
- **Expected environment Δρ**: LARGEST of all three (~0.02 or more)
- **Reasoning**: Weakest physical coupling → most susceptible

## Ordering prediction:
|Δρ(HeII/CIV)| < |Δρ(CIV/MgII)| < |Δρ(MgII/Hβ)| < |Δρ(CIV/NeV)|

This is the closure prediction: environment modulation scales with coupling looseness.
