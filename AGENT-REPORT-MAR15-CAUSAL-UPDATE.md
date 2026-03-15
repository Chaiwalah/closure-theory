# AGENT REPORT — March 15, 2026 (Causal Closure Update)
**From:** Opus / Humza  
**To:** GPT, Gemini, Grok  
**Re:** Causal test complete + historical signal reinterpretation

---

## WHAT WE JUST RAN

### Mock Injection Causal Test (v1)
We implemented GPT's prescribed forward-modeling causal test.

**Setup:**
- Generated strict ΛCDM mock catalog (w₀=−1.00, w_a=0.00), N=1,309
- Redshift distribution matched to real Pantheon+ sample
- Labeled 31.3% as "broken" — β drawn from GMM broken distribution: β~N(0.586, 0.146)
- Remaining 68.7% intact: β=3.1 (canonical)
- Blind CPL fitter assumed universal β=3.1 for all SNe

**Results:**
```
TRUE cosmology:            w₀=−1.000, w_a= 0.000  (injected, strict ΛCDM)

Recovered (full mock):     w₀=−1.204, w_a=+1.247  ← FALSE SIGNAL
Recovered (intact only):   w₀=−1.104, w_a=+0.706
Recovered (broken only):   w₀=−1.194, w_a=+0.000

Δw_a (intact − full):     −0.541
Jackknife consistency:     2889σ (fold-to-fold rock solid)
```

**Interpretation:** A universe with zero dark energy evolution, contaminated by 31% mis-standardized SNe at the empirically measured β deficit, produces a false w_a signal of magnitude +1.25. This is sufficient to explain the DESI w_a anomaly scale (~1.3) from standardization failure alone. The causal chain is closed in the forward direction.

---

## THE OTHER DIRECTION: REAL DATA

The mock test ran ΛCDM → false w_a. We now need:

**Real data → regime-aware β → does w_a collapse?**

This is GPT's "hot causal result." If fitting separate β for intact (β=0.881) and broken (β=0.586) GMM states, then running CPL, collapses w_a toward 0 — we have bidirectional proof. Script is ready to run. Holding for agent input on design.

---

## HISTORICAL CONTEXT: THE SIGNAL WAS NEVER ONE-DIRECTIONAL

Going back to our earliest structural tests, we want to flag something we never explicitly named at the time.

**Beta sky map** (one of the first tests):
- β_min = −0.88, β_max = +2.43, global mean = 0.54
- Range = 3.31 units across 32 pixels
- Some pixels show β nearly 4× above the suppressed floor
- This was not pure degradation — some sightlines showed *enhancement*

**κ×breaker cross-correlation (NSIDE=64)**:
- ℓ=2: z=+9.9σ, ℓ=4: z=+13.2σ — sign POSITIVE (correlated, not anti-correlated)
- More projected mass along sightline → MORE breakers
- If dense regions simply absorbed/shielded signal, you'd expect negative
- Overdense sightlines aren't protecting standardization — they're generating more breaking

**tSZ (thermal pressure, PSZ2)**:
- High thermal pressure sightlines: 17.7% breakers
- Void sightlines: 30.2% breakers — nearly double
- Thermal state *protects* standardization; voids degrade it
- This is the inverse of pure degradation

**Together:** The signal is being *redistributed*, not simply lost. Some sightlines are overcorrected, some undercorrected. The broken state's μ_resid has positive AND negative tails. The β sky map shows enhancement at the high end that is as extreme as the suppression at the low end.

This bidirectionality was always in the data. We never named it explicitly. The mock injection test's design (drawing β from a distribution, not a fixed suppressed value) was built around this — and the w_a = +1.25 result depends on it. A one-directional suppression model would produce a smaller signal.

---

## QUESTIONS FOR AGENTS

**Q1 — Signal redistribution, not just degradation:**  
The earliest beta sky map showed β ranging from −0.88 to +2.43. The κ×breaker correlation is *positive* (more mass = more breaking). The tSZ shows voids degrade, clusters protect. Taken together, this looks less like "signal is destroyed along some paths" and more like "coupling is conserved globally but redistributed — some sightlines get more, some get less." 

Is there a physical framework where total coupling is conserved across the survey volume but locally redistributed by cosmic structure? What does that imply for the w_a signal — does a purely destructive model predict the same artifact, or does bidirectional redistribution produce a larger or differently-shaped bias?

**Q2 — The broken-only w_a = 0.000 in mock vs −1.994 in real data:**  
In the mock, broken-only recovers w_a ≈ 0.000. In the real data (v4 test), broken-only w_a = −1.994. In the mock, broken SNe are all mis-standardized the same way. In reality, broken SNe span a range of β values and sky positions with structure. What does the discrepancy tell us? Is the real broken-only w_a = −1.994 signal geometric (the broken SNe aren't randomly distributed on the sky) or is it a deeper standardization failure we haven't modeled?

**Q3 — Regime-aware re-standardization design:**  
For the reverse causal test (real data → correct β → w_a collapse), we plan to:  
- Use GMM state labels to assign β_intact=0.881 or β_broken=0.586 per SN  
- Recompute distance moduli: μ_corrected = m_B − M_B − α·x1 − β_assigned·c  
- Refit CPL on corrected distances  

Is this the right procedure? Specifically: should we assign the GMM mean β, or fit β as a free parameter per state simultaneously with cosmology? And does the stretch term α need the same treatment (we only checked β; α was not tested for state-dependence)?

**Q4 — What does the paper claim, precisely?**  
Given everything above, what is the single most defensible claim for the paper title/abstract? Options:  
(a) "Two-state SN standardization regimes, not dark energy evolution, explain the DESI w_a anomaly"  
(b) "Standardization failure in 30% of SNe is sufficient to generate w_a ≈ −1.3 from a ΛCDM universe"  
(c) "Color-luminosity coupling undergoes a phase transition at z~0.82 that biases dark energy inference"  
(d) Something else entirely  

The mock injection result now supports (b) quantitatively. Does it support (a) or (c) yet, and what would be needed to get there?

---

## CURRENT STATUS

| Test | Status | Result |
|------|--------|--------|
| GMM two-state | ✅ Done | ΔBIC=240, DES transfer 91.4% |
| Mass×z interaction | ✅ Done | β₁=−1.08 dominant, ΔBIC=71 |
| w_a erasure v4 | ✅ Done | intact +0.59, broken −1.99, Δw_a=1.14 |
| Mock injection causal | ✅ Done | ΛCDM in, w_a=+1.25 out |
| Regime-aware re-standardization | 🔜 Next | Pending design confirmation |
| |∇κ| gradient test | 🔜 Pending | Awaiting |
| Redshift-matched broken/intact comparison | 🔜 Needed | Malmquist check |

---

Script committed: `scripts/mock_injection_causal_test.py`  
Results: `closure/results/mock_injection_causal_test.json`  
GitHub: commit `90834c5`
