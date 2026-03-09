# CALCULATION ASSIGNMENT — The Final Wall
## Modular Tasks for Grok and Gemini
### 7 March 2026

---

## OBJECTIVE

GPT says one thing stands between us and a closed mechanism: the actual numbers. No more narrative. No more qualitative arguments. Just math. Each task below is independent — if one fails, the others still stand. Complete as many as you can.

---

## TASK 1: THE CAST NUMBER (Critical — this decides everything)

**The equation:**
Γ₀ = (g_aγ × B_⊥ × L_coh)² / 4

**Known values:**
- Γ₀ = 2.17 (measured from 5 independent density tests, obs/pred = 1.03 ± 0.02)
- B_⊥ (IGM magnetic field): use range 0.1 – 10 nG (observational constraints from Faraday rotation, CMB, blazar observations)
- L_coh (magnetic coherence length): use range 0.1 – 10 Mpc (from IGM simulations)

**Calculate:**
1. Solve for g_aγ at the central values (B = 1 nG, L_coh = 1 Mpc)
2. Solve for g_aγ across the full parameter space (B × L_coh grid)
3. Compare every solution to CAST bound: g_aγ < 5.8 × 10⁻¹¹ GeV⁻¹ (95% CL, low mass)
4. Also compare to IAXO projected sensitivity and astrophysical bounds (globular cluster, SN1987A)
5. **Verdict:** Is there ANY region of realistic (B, L_coh) parameter space where the required g_aγ is BELOW the CAST bound?

**Unit conversion reminder:** Make sure B is in Tesla (1 nG = 10⁻¹³ T), L_coh is in meters (1 Mpc ≈ 3.086 × 10²² m), and Γ₀ is dimensionless (it's an optical depth parameter). Show your unit analysis.

**If the number passes CAST:** GPT's last quantitative objection dies.
**If the number fails CAST:** We need to understand why, and whether the decoherence-vs-conversion argument saves it (decoherence needs less coupling than conversion by how many orders of magnitude?).

---

## TASK 2: THE MICROSCOPIC LINDBLAD DERIVATION (The theoretical core)

GPT's specific objection: "L_k ∝ σ_z with γ_k ∝ g_k² is an ansatz, not derived from the Hamiltonian."

**The assignment:**

Start from the full Hamiltonian:
H = H_photon + H_ALP + H_IGM + H_int

where H_int = (1/4) g_aγ F_μν F̃^μν_ALP (the Primakoff coupling)

In the presence of a stochastic external magnetic field B_IGM(x) with correlation function:
⟨B_i(x) B_j(x')⟩ = (B²/3) δ_ij × C(|x-x'|/L_coh)

Use the Nakajima-Zwanzig projection operator technique OR the Born-Markov master equation derivation to:

1. Write the reduced dynamics of the photon subsystem after tracing over the ALP field and stochastic B-field
2. Show that the resulting Lindblad operators naturally couple to the polarization coherences
3. Show that the dephasing rate for a photon with initial coherence d_i gives γ_i ∝ |d_i|²
4. Combined with Deliverable A (d_i ∝ g_i), this gives γ_i ∝ g_i⁴ (or g_i² — clarify the exact scaling)

**The key step GPT wants:** The Lindblad operators must EMERGE from the Hamiltonian + environment correlator. They cannot be inserted.

---

## TASK 3: SOURCE-SIDE DECOHERENCE SURVIVAL (The weakest link)

GPT keeps asking: why don't the magnetic coherences decohere at the source before the photon escapes?

**Calculate:**

1. **Emission timescale** for each line type:
   - Permitted lines (CIV, Hβ): τ_emission ~ 10⁻⁸ s (spontaneous, A-coefficient)
   - Forbidden lines ([OIII], [NII], [SII]): τ_emission ~ 1–100 s
   - Semi-forbidden ([OII]): τ_emission ~ 10⁻²–10⁻¹ s

2. **Source-side magnetic decoherence timescale:**
   - BLR magnetic field: ~1–100 Gauss
   - NLR magnetic field: ~10⁻⁴–10⁻² Gauss
   - Estimate τ_decoherence = 1/(g_i² × μ_B² × B_source² × τ_correlation)
   - where τ_correlation is the correlation time of B-field fluctuations at the source

3. **Compare:** If τ_emission < τ_decoherence for each line, the coherences are "frozen in" before source-side decoherence can act. If τ_emission > τ_decoherence for some lines, those lines lose their coherences at the source.

4. **Critical question:** Even if partial source-side decoherence occurs, does the RESIDUAL coherence still scale with g_i? (Answer should be yes — partial decoherence reduces the magnitude but preserves the g_i dependence.)

---

## TASK 4: THE g² vs g TEST (Empirical — can we do this NOW?)

The derivation predicts γ ∝ g_i² (quadratic), not γ ∝ g_i (linear).

**Data we have (6 emission lines):**

| Line | Degradation Rate | Diagnostic Sensitivity (q) |
|------|-----------------|---------------------------|
| [NII] 6585 | 0.000 | 0.0 |
| [OIII] 5007 | 0.000 | 0.0 |
| Hβ 4861 | -0.038 | 0.3 |
| [OII] 3727 | -0.179 | 0.4 |
| CIV 1549 | -0.289 | ~0.6 |
| [SII] 6718 | -0.396 | 0.7 |

**What we need:** The actual Landé g-factors for these transitions from NIST atomic data.

**Calculate:**
1. Look up the effective Landé g-factor for each transition (upper and lower levels, compute g_eff)
2. Compute correlation: degradation rate vs g
3. Compute correlation: degradation rate vs g²
4. **If |r(g²)| > |r(g)|:** The quadratic scaling is confirmed → derivation predicts the data
5. **If |r(g)| > |r(g²)|:** Linear scaling fits better → derivation needs correction

This is the smoking gun empirical test of the mechanism.

---

## TASK 5: CIV DOUBLET NUMERICAL PREDICTION (Bonus — the kill shot)

The CIV doublet (λ1548.2, λ1550.8) has two components with DIFFERENT g-factors.

**Calculate:**
1. Look up the exact Landé g-factors for CIV λ1548 and CIV λ1550 transitions
2. From the mechanism: γ_1548/γ_1550 = (g_1548/g_1550)²
3. Predict the expected asymmetry drift: Δv(z) from the differential decoherence
4. Compare to the observed drift rate (r = +0.995 across 4 z-bins from 311K quasars)
5. Does the predicted magnitude match? Does the predicted sign match?

**If this works:** The mechanism makes a NUMERICAL prediction for a doublet that matches observations. That's not a qualitative story. That's a quantitative mechanism.

---

## DELIVERY

Each task is self-contained. Do whichever ones you can. Number your responses by task. If a task requires assumptions, STATE THEM CLEARLY so we can evaluate.

The goal: give GPT numbers, not words. He's asked for numbers four times now. Let's give him numbers.

---

*"You do not need a giant institution to make progress."* — GPT-4, March 7, 2026
*"Then let's make progress."* — Us, right now.
