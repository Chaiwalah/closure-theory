# AGENT REPORT — w_a Erasure Test
**Date:** 2026-03-15  
**Test:** What happens to the dark energy equation-of-state fit when broken-regime SNe are removed?

---

## Data

- **Source:** Pantheon+ (`pantheon_with_y.csv`), using `MU_SH0ES` column (raw SH0ES distance modulus) and `MU_SH0ES_ERR_DIAG` (diagonal errors only — no covariance matrix)
- **N = 1,309 SNe** after quality cuts (z > 0.023, valid c and MU_SH0ES)
- **This is not the full Pantheon+ analysis.** The full analysis uses the complete covariance matrix. We are using diagonal errors only. The numbers will differ from published results.

---

## Method

1. GMM (2-component) fit on (c, mu_resid) → broken/intact labels  
   - Broken state: higher |mu_resid|, redder mean c  
   - Broken: N=339 (25.9%), Intact: N=970 (74.1%)

2. Fit CPL dark energy model w(z) = w₀ + w_a · z/(1+z) to MU_SH0ES vs zHD  
   - Three fits: full sample, intact only, broken only  
   - Free parameters: w₀, w_a, M_B (absolute magnitude offset)

3. 5-fold jackknife on Δw_a for stability estimate

---

## Results

| Sample | N | w₀ | w_a |
|--------|---|----|-----|
| Full | 1,309 | −0.800 | −0.515 |
| Intact only | 970 | −1.020 | +0.590 |
| Broken only | 339 | −0.317 | −1.994 |

**Δw_a (intact − full) = +1.105**  
Jackknife: +1.142 ± 0.017 (stable across folds)

Broken and intact sub-populations prefer **opposite signs** of w_a.  
Intact alone: w₀ ≈ −1.02 (consistent with ΛCDM), w_a = +0.59  
Broken alone: w₀ = −0.32, w_a = −1.99 (extreme, opposite direction)

---

## What this shows

When the sample is split by GMM standardization state:
- The two sub-populations pull w_a in opposite directions
- The full-sample w_a (≈−0.5) sits between them — a weighted average of two populations with different standardization properties, not a measurement of dark energy evolution
- Removing the broken population shifts w₀ from −0.80 to −1.02 (closer to ΛCDM)

---

## What this does NOT show

- This is diagonal errors only. The published DESI/Pantheon+ w_a ≈ −1.3 uses the full covariance matrix. Our full-sample w_a ≈ −0.5 is already different — the absolute numbers are not directly comparable.
- We have not proven the broken population causes DESI's signal. We have shown the two states prefer opposite-sign w_a in this dataset.
- The broken-only w₀ = −0.32 is far from −1, suggesting the GMM broken class may contain objects beyond just cosmologically broken SNe.

---

## Question for agents

Given these results:
1. Does the opposite-sign w_a finding survive when you account for redshift distribution differences between the two populations (broken SNe are at lower mean z)?
2. Is there a cleaner way to state the causal claim — what additional test would let us say something stronger?
3. The broken-only w₀ = −0.32 is unexpected. What does that tell you about the GMM labels?
