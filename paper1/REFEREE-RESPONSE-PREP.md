# Pre-emptive Referee Response Guide
## Paper: "Geometric Non-Stationarity of the Type Ia Supernova Standardization Manifold and Its Impact on Dark Energy Inference"

*Prepared 9 March 2026 — based on adversarial review by GPT, Grok, and Gemini*

---

## Attack 1: "Over-fitting"

**Expected complaint:** "The improvement in χ² and the reduction in ρ are simply the result of adding five additional free parameters to a large dataset."

**Defense:**

1. **Significance-per-parameter:** Δχ² = 601 for 5 parameters = ~120 per parameter. The interaction-only model (2 params) gives Δχ² = 269 — that's 134.5 per parameter, far exceeding any conventional threshold for nuisance curve-fitting.

2. **Box's M is parameter-free:** Section 3.1 proves the manifold is nonstationary *before* a single extra parameter is added. Box's M (p = 2.3 × 10⁻⁷) is a standard test of the raw covariance structure, not a model fit.

3. **Cross-domain anchor:** Over-fitting in a supernova pipeline cannot explain identical sigmoid-structured coupling collapses in 750,000 quasars (SDSS DR16Q) or 535 FRBs (CHIME).

---

## Attack 2: "Selection Bias / Malmquist"

**Expected complaint:** "The evolution of β(z) is likely an artifact of Malmquist bias. At higher redshifts, we only see the brightest, bluest supernovae, creating a false correlation between color and distance."

**Defense:**

1. **Sign reversal:** Malmquist bias might change the *magnitude* of a correlation, but it rarely causes a *sign reversal*. corr(c, Δμ) flips from +0.09 to −0.44 — this is a hallmark of basis rotation, not a selection cutoff.

2. **Kinematic immunity:** If the effect were purely selection bias, it should affect stretch (x₁) and color (c) similarly. Instead, x₁ remains relatively locked while c degrades — consistent with a channel-selective mechanism, not a brightness threshold.

3. **Quality cuts already applied:** |x₁| < 5, |c| < 0.5, σ(mB) < 1.0 — standard Pantheon+ cuts that mitigate the most obvious selection effects.

---

## Attack 3: "Progenitor Evolution"

**Expected complaint:** "Supernovae at z ~ 1.5 are younger and come from different stellar populations. You are measuring stellar evolution, not geometric non-stationarity."

**Defense:**

1. **Sigmoid vs SFH:** Progenitor evolution follows Star Formation History — a gradual, smooth curve. The effects we observe follow path-length-dependent sigmoids with sharp thresholds.

2. **Cross-domain kills it:** It is astrophysically impossible for SN progenitors, quasar accretion disks, and FRB engines to all evolve with the same geometric signature across five orders of magnitude in frequency. Three independent source classes, three independent instruments, one pattern.

3. **Kinematic invariance:** If progenitor evolution drove the effect, stretch should be affected too (different ⁵⁶Ni yields at different metallicities). Instead, stretch is relatively immune — consistent with a channel-selective propagation effect, not a source population drift.

---

## Attack 4: "H₀–wₐ Degeneracy"

**Expected complaint:** "By allowing α and β to evolve, you've just moved the uncertainty into the standardization. You haven't solved anything."

**Defense:**

1. **Geometric necessity:** We are not "moving" uncertainty — we are identifying its source. The covariance structure is *measurably* different at low vs high z (Box's M, p = 10⁻⁷). A fixed basis applied to a nonstationary manifold *necessarily* produces structured residuals. The correction is geometric, not ad hoc.

2. **Model comparison:** ΛCDM + nonstationary Tripp (fewer cosmological free parameters) outperforms w₀wₐ + standard Tripp by Δχ² = 37. The standardization correction is *more parsimonious* than the dark energy extension.

3. **Direction is physical:** The H₀ shift toward Planck (73.0 → 68.8) is not a manual adjustment — it's the natural result of removing the color-channel bias. The direction is independently predicted by the channel-count hierarchy.

---

## Attack 5: "Residual Floor"

**Expected complaint:** "Your corrected residuals still show ρ = −0.135 at 7σ. The correction is incomplete. How do you know the remaining structure isn't evidence against your model?"

**Defense:**

1. **Acknowledged and predicted:** Section 8 explicitly states the linear repair is incomplete. The remaining structure is consistent with manifold curvature and saturation at z > 1, where the [1.0, 1.5) bin shows over-correction under linear terms.

2. **Mock null confirms partial correction:** The mock null injection test shows pre-correction 8.2σ → post-correction 7.0σ, with the correction itself significant at 4.6σ. The floor is real but the correction is also real.

3. **Higher-order terms exist:** The quadratic color term (Δχ² = 339) and evolving δ(z) capture some of this curvature. A full nonlinear manifold correction (polynomial or spline basis) is reserved for follow-up work.

---

## Attack 6: "Covariance Mismatch" (GPT's primary concern)

**Expected complaint:** "You use diagonal uncertainties while the Pantheon+ analysis uses a full covariance matrix. Your significance levels may be inflated."

**Defense:**

1. **Explicitly acknowledged:** Limitation #2 states this clearly. The public covariance is in post-standardization μ-space; our corrections operate upstream in (mB, x₁, c) space.

2. **Jacobian propagation planned:** A formal Jacobian propagation from the correction space into μ-space covariance is planned for follow-up using the SALT2 parameter-level covariance.

3. **Box's M is independent:** The covariance nonstationarity test operates on the raw (x₁, c, Δμ) vectors, not on goodness-of-fit. Even with a full covariance, the *structure* of the covariance changing with z remains.

4. **Conservative direction:** Keeley et al. (2023) showed the Pantheon+ covariance may be *overestimated* by ~7%. If anything, proper covariance treatment would strengthen our significance levels.

---

## Attack 7: "Nature of Nonstationarity" (GPT's geometry caveat)

**Expected complaint:** "You call this 'rotation' but haven't proven it's purely rotation. It could be scale change, curvature, or something else."

**Defense:**

1. **Acknowledged:** Limitation #7 explicitly states: "the rolling-window and covariance-based tests are descriptive of nonstationarity [and] do not by themselves uniquely determine whether the underlying geometry is best described as pure rotation, curvature, scale change, or a combination thereof."

2. **The claim is nonstationarity, not rotation:** The paper's title and core claim are about *nonstationarity*, not a specific geometric transformation. PCA drainage is presented as illustration (Section 4), not proof.

3. **The practical consequence is identical:** Whether the manifold rotates, curves, or scales, the fixed Tripp basis fails. The correction is the same regardless of the geometric classification.

---

## Summary: Lines of Defense

| Attack | Primary Shield | Secondary Shield | Tertiary Shield |
|--------|---------------|-----------------|-----------------|
| Over-fitting | Box's M (no params) | 134.5 Δχ²/param | Cross-domain |
| Selection bias | Sign reversal | Kinematic immunity | Quality cuts |
| Progenitor evolution | Cross-domain (3 classes) | Sigmoid vs SFH | Stretch immunity |
| H₀–wₐ degeneracy | Model comparison (Δχ²=37) | Geometric necessity | Direction physical |
| Residual floor | Acknowledged + predicted | Mock null (4.6σ) | Higher-order planned |
| Covariance mismatch | Explicitly acknowledged | Jacobian planned | Keeley et al. direction |
| Nature of geometry | Claim is nonstationarity | PCA is illustrative | Consequence identical |

---

*This document is for internal preparation. Do not submit with the manuscript.*
