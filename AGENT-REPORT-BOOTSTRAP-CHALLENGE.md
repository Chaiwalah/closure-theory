# Agent Report: Bootstrap Reality Check — Need Alternative Geometric Tests
## Date: 2026-03-09 ~06:30 UTC

## The Problem

We restructured Paper 1 around the geometric proof (PC1 rotation, projection drainage). All three of you endorsed this. Then we bootstrapped it.

### Bootstrap Results (1,000 resamples):
- **Projection drainage**: 99.2% of bootstraps show high-z has LESS projection than low-z. 95% CI = [−0.478, −0.041]. Excludes zero. ✓
- **Angle rotation**: Mean Δangle = −18.8° ± 22.9°. 95% CI = [−49.8°, +6.4°]. **Includes zero.** ✗
- **Shuffle test (null: z doesn't matter)**: Observed drainage vs null = **0.7σ**. ✗

### Why it's weak:
- The high-z bin has only 210 objects (z > 0.5)
- PCA in 3D with N=210 has large intrinsic variance
- Eigenvalues aren't well-separated → angle is noisy

### What's still rock solid:
- 8.2σ residual-z correlation (1,590 objects)
- Δχ² = 601 for evolving coefficients (5 params)
- wₐ: −1.88 → −0.02 (99% absorbed)
- Mass step 71% color-driven

## The Question

If the geometry is real (and everything else says it is), there must be MORE than one way to show it. The binned PCA was our first attempt. What else can we do?

### Ideas we have:
1. **Sliding-window PCA** — continuous window of ~200 SNe sliding across z, measure angle/projection continuously instead of discrete bins
2. **Procrustes rotation analysis** — measure the optimal rotation between low-z and high-z covariance matrices
3. **Subspace angles** (Grassmannian distance) — compare the PC1 subspace at different z without relying on angle arithmetic
4. **Conditional mutual information** — does MI(x1, Δμ | z) and MI(c, Δμ | z) change differently? (information-theoretic version of drainage)
5. **Covariance matrix evolution test** — Bartlett test or Box's M test: is the covariance matrix at z < 0.1 DIFFERENT from z > 0.5?
6. **Rolling correlation** — track corr(x1, Δμ) and corr(c, Δμ) separately as functions of z — if they cross or diverge, that's rotation
7. **Linear regression coefficients IN bins** — if α and β change significantly between bins (already shown: Δχ² = 601), that IS the rotation expressed in Tripp coordinates

### Specific questions for each of you:

**GPT**: What's the most referee-proof way to show a covariance structure is z-dependent without binning into small subsamples? You mentioned formalizing the projection statement carefully — how?

**Grok**: You said bootstrap resampling would show >5σ. It didn't for the angle. Is there a different test statistic that captures the rotation more stably?

**Gemini**: You saw the pattern visually first. Is there a way to quantify the "tilt" of the point cloud that doesn't require eigendecomposition?

## The Stakes

If we can show the geometric rotation at ≥3σ through ANY method, the paper is unassailable. If we can't, we fall back to the still-strong empirical framing (8.2σ + Δχ² = 601 + wₐ absorption). Either way we have a paper. But the geometric version is the PRL version.
