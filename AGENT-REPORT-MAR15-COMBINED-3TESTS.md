# AGENT REPORT — Three New Tests
**Date:** 2026-03-15  
**Status:** Day 2 — Results for agent review

These are the three tests run today that agents haven't seen yet. Results only — no interpretation beyond what the numbers say directly.

---

## TEST 1 — Two-State GMM: Is Pantheon+ one population or two?

**Method:** Gaussian Mixture Model (2-component) fit jointly on (c, mu_resid) across 1,473 SNe.

**Results:**

| State | N | % | Mean β (residual) | Mean \|mu_resid\| | Mean c | Breaker rate |
|-------|---|---|--------------------|-------------------|--------|--------------|
| Intact | 1,027 | 70% | 0.881 | 0.094 | −0.059 | 14.6% |
| Broken | 446 | 30% | 0.586 | 0.211 | +0.065 | 54.0% |

- **ΔBIC = 240.0** (two-state vs one-state). ΔBIC > 10 = very strong.
- **ΔAIC = 271.8** — same conclusion independently.
- **Chi2 = 245.9, p = 2×10⁻⁵⁵** — GMM state assignment and breaker status are not independent.
- Residuals non-Gaussian: p = 2.5×10⁻¹⁴ (they are a mixture, not a deformed single distribution).
- **Zero-shot transfer to DES5YR (440 SNe, independent survey): 91.4% agreement.** DES internal ΔBIC = 35.8 (two-state wins independently).

**Cross-survey consistency:**

| Dataset | ΔBIC | Broken % | Broken mean c | Broken mean \|MURES\| |
|---------|------|----------|---------------|----------------------|
| Pantheon+ | 240.0 | 30% | +0.065 | 0.211 |
| DES5YR | 35.8 | 32% | +0.085 | 0.144 |

The two states reproduce in an independent survey with a frozen model — no refitting.

**One anomaly:** Bootstrap of sky-position mapping (Test B) was not significant (p=0.234, NSIDE=4). The GMM state labels do not cleanly map onto great circle distance bins. The breaker map is a geometric object; the GMM labels are a standardization-regime object. They overlap but are not identical.

---

## TEST 2 — Mass Step: What drives it?

**Method:** Extended Tripp correction with γ(z) = γ₀ + γ₁·z term (host mass × redshift interaction). N=1,473 SNe.

**Results:**

- Raw mass×z correlation: ρ = −0.091, **p = 0.00044, 3.5σ**
- After standard Tripp correction: ρ = 0.068, **2.6σ** (residual remains)
- Full α(z)+β(z)+γ(z) model vs standard Tripp: **ΔBIC = 71** (decisive — extended model wins)
- γ₁ = +0.057 (slope of mass correction with z)
- β₁ = −1.08 (color evolution term dominant)
- Standard mass step: −0.0152 mag
- Extended model mass step: −0.0154 mag
- **γ(z) absorbs 1.5% of the mass step** (not 71% — agents may have inflated this number)

**What the numbers say:** The extended model fits better (ΔBIC=71), driven by β evolving with z (β₁ = −1.08), not by the mass step dissolving. The mass step itself is nearly unchanged by adding γ(z).

---

## TEST 3 — w_a Erasure: What happens to dark energy fits when broken SNe are removed?

**Data note:** Uses `MU_SH0ES` column (SH0ES distance modulus, diagonal errors only). Not the full Pantheon+ covariance matrix. Absolute values differ from published DESI/Pantheon+ results — the comparison is internal only.

**Method:** CPL fit w(z) = w₀ + w_a·z/(1+z) to MU_SH0ES vs zHD. Free parameters: w₀, w_a, M_B. Three fits: full sample, intact-only, broken-only. N=1,309 SNe after quality cuts (zHD > 0.023, valid MU_SH0ES).

**Results:**

| Sample | N | w₀ | w_a |
|--------|---|----|-----|
| Full | 1,309 | −0.800 | −0.515 |
| Intact only | 970 | −1.020 | +0.590 |
| Broken only | 339 | −0.317 | −1.994 |

- **Δw₀ (intact − full) = −0.22** (intact alone gives w₀ = −1.02, consistent with ΛCDM)
- **Δw_a (intact − full) = +1.10**
- **Jackknife: +1.14 ± 0.02** — stable across all 5 folds, not a fluke
- **Broken and intact prefer opposite signs of w_a** (broken: −1.99, intact: +0.59)
- The full-sample w_a ≈ −0.5 sits between them — a weighted mixture

**What the numbers don't say:** We have not shown the broken regime causes DESI's published w_a ≈ −1.3. Our full-sample w_a ≈ −0.5 differs from DESI's already (diagonal-only errors). The finding is that the two sub-populations pull in opposite directions.

---

## Questions for agents

1. **Test 1 — GMM geometry gap:** The two states don't map cleanly onto great circle distance. Broken SNe are at lower mean z (median 0.159 vs 0.222). Is this a redshift confound, or does the z-distribution difference itself encode the geometry (walls at lower z are more populated)?

2. **Test 2 — Mass step:** β₁ = −1.08 (β evolves with z) is the dominant signal, not γ(z). Does β evolving with redshift have a natural explanation in the closure framework, or does this need its own mechanism?

3. **Test 3 — w_a sign flip:** The broken-only fit gives w₀ = −0.32 (far from −1). What does a sub-population of SNe with w₀ = −0.32 mean physically, or is this just what a badly standardized sample looks like when forced through a CPL fit?

4. **Across all three:** Is there a single test that would let us make a causal claim — that standardization failure in the broken regime is responsible for dark energy evolution signals — rather than just showing correlation?

---

*Scripts: `wa_erasure_test_v4.py`, `mass_z_interaction_test.py` (closure-theory/scripts/)*  
*Results: `wa_erasure_test_v4.json`, `two_state_mixture_tests.json`, `mass_z_interaction_test.json`*
