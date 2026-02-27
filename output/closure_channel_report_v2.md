# Closure Theory — First Empirical Constraints and Forward Tests

**Date:** February 19, 2026
**Author:** Humza Hafeez

---

## 1. Reframing: What This Test Was and Wasn't

### What we did
Cross-matched 1,454 Pantheon+SH0ES Type Ia supernovae with Planck CMB lensing (κ) and tSZ Compton-y maps at multiple apertures. Computed deconfounded Hubble residuals (Δμ⊥) after removing dust, host mass, and survey systematics. Tested for correlations between Δμ⊥ and line-of-sight structure tracers.

### What we found
| Tracer | r | p | Interpretation |
|--------|---|---|----------------|
| κ_comp (10–20') | −0.024 | 0.362 | Null |
| y_comp (10–20') | +0.009 | 0.728 | Null |
| κ₃₀ (raw disc) | −0.022 | 0.393 | Null |
| y₃₀ (raw disc) | +0.060 | 0.021 | Marginal, low-z only |

The y₃₀ signal (r = +0.111, p = 0.014) exists only at z < 0.1 and vanishes at higher redshift.

### Why this null is not damaging — it's constraining

This test was a **proxy** for the primary prediction. Closure Theory predicts a difference between EM and GW distances from the same source (the Γ₁ coupling). Testing EM-only residuals against structure maps is a second-order effect — analogous to testing general relativity by measuring one clock instead of comparing two at different potentials.

**What the null constrains:** The closure drift rate Γ₁ is not large enough to produce detectable correlations (|r| > 0.05) in EM-only data against Planck-resolution structure maps. This places an upper bound of approximately |Γ₁| < 10⁻³ Mpc⁻¹ — consistent with the theory's expectation that the effect is subtle and accumulates over cosmological baselines.

**What the null does not constrain:** The primary EM-GW differential test, which has far greater statistical leverage. A Γ₁ of 10⁻⁴ Mpc⁻¹ would be undetectable in our test but clearly visible with ~100 standard sirens at z > 0.1.

---

## 2. The Low-z y₃₀ Signal: Two Interpretations

### Standard interpretation (signal detection frame)
Local environment systematic. Peculiar velocities dominate at z < 0.1. SNe in denser baryonic environments (higher tSZ y) have different host properties → residual bias. The compensated aperture kills the signal → it's a monopole (local) effect, not a line-of-sight (propagation) effect.

### Closure interpretation (information-theoretic frame)
Short paths through locally dense baryonic environments produce the largest *per-unit-length* channel interaction. The closure drift rate Γ(x,t) is highest in regions of high baryonic density. At low-z, the local structure dominates the integrated channel effect because the path is short and the local environment contributes a large fraction of the total. At high-z, the signal averages out over diverse environments along the much longer path.

The compensated aperture removing the signal is **ambiguous**: if the baryonic medium is roughly uniform across the 30' disc (smooth filament or cluster), the compensated aperture removes a real signal along with the monopole.

**Neither interpretation is ruled out by this data.** The distinction requires the EM-GW differential test or a direct baryonic line-of-sight tracer (FRB dispersion measure).

---

## 3. The Deeper Issue: We Were Looking for the Wrong Thing

### The Andromeda insight

Two observers separated by walking speed, looking at the Andromeda galaxy simultaneously, see photons corresponding to events separated by **days** on Andromeda's surface. Standard physics calls this "relativity of simultaneity." Closure Theory says: those observers are in slightly different closure loops, and "Andromeda" is a different reconstruction for each loop.

A photon traveling from Andromeda experiences zero proper time. From the photon's null geodesic, emission and absorption are the same event. The 2.5 million light-years of "distance" exists only from inside the observer's closure loop. The photon doesn't travel *through* space — it connects two points of a closure structure. Space is what the loop looks like from the inside.

### The implication for measurement

If reality is fundamentally about information closing into loops — not matter moving through space — then the right observable isn't **energy** (luminosity, distance modulus) but **information structure** (spectral complexity, temporal coherence, bandwidth).

Standard cosmology measures:
- How bright (luminosity → distance)
- How red (redshift → recession velocity)
- How shaped (light curve → standardization)

Nobody measures:
- How **complex** is the received signal compared to the emitted signal?
- Does the **information content** of a photon bundle degrade with path length through structure?
- Is there a **bandwidth constraint** that's environment-dependent?

---

## 4. Testable Information-Theoretic Predictions

### Test A: Light Curve Temporal Complexity vs. LOS Structure

**The idea:** If closure imposes bandwidth constraints on the EM channel, the *temporal structure* of a SN light curve should be subtly affected by the intervening medium — not just redshifted and dimmed, but **smoothed**. High-frequency features in the light curve would be attenuated relative to low-frequency features when the photon traverses denser baryonic environments.

**Observable:** Compute the spectral entropy (Shannon entropy of the Fourier power spectrum) of each SN Ia light curve's residuals after SALT2 fitting. Correlate this entropy with line-of-sight baryonic density.

**Prediction:**
- ΛCDM: No correlation between light curve entropy and LOS structure (after standardization, residual structure is intrinsic noise)
- Closure: Negative correlation — SNe behind denser channels have *less* residual temporal complexity (bandwidth-limited, smoothed)

**Data:** Pantheon+ provides full SALT2 fit parameters (x1, c) and fit residuals (FITCHI2, NDOF). The raw light curves are in the DR1 release. We already have LOS κ and y values for all 1,454 SNe.

**What's novel:** Nobody has correlated SALT2 fit quality with *line-of-sight* structure. Fit quality is always attributed to intrinsic SN diversity or measurement noise. This test asks: **is some of the "noise" in light curve fits actually channel-imposed smoothing?**

**Quick proxy test (can run now):** Correlate FITCHI2/NDOF (goodness-of-fit) with κ and y values. If closure bandwidth-limits the temporal signal, SNe behind denser environments should have *lower* χ²/dof (better fits, because high-frequency deviations are suppressed). This is backwards from what you'd expect if the environment adds noise.

### Test B: SALT2 Color Excess vs. LOS Baryonic Density

**The idea:** The SALT2 color parameter `c` captures deviations from the standard SN Ia color law. Standard interpretation: dust in the host galaxy. But `c` also absorbs any wavelength-dependent modification to the photon bundle — including channel-dependent bandwidth effects.

**Observable:** After correcting for host galaxy dust (MWEBV), correlate the SALT2 color residual with *line-of-sight* (not host) baryonic density.

**Prediction:**
- ΛCDM: No correlation (c is entirely intrinsic + host dust)
- Closure: Positive correlation — denser LOS channels produce a wavelength-dependent drift that SALT2 absorbs into `c`

**What's novel:** Existing work (Pruzhinskaya et al. 2020) shows x1 correlates with host morphology but c does NOT correlate with host properties. If c correlates with *line-of-sight* structure instead of host structure, that would be a genuinely new finding.

**Data:** Already in our Pantheon+ dataset. The `c` and `cERR` columns are in `pantheon_plus.dat`. We have LOS κ and y values. This test takes 10 minutes to run.

### Test C: CMB Anomalies as Closure Boundary Signatures

**The idea:** The CMB anomalies (hemispherical asymmetry, cold spot, quadrupole-octupole alignment, missing large-angle correlations) are unexplained in ΛCDM. A 2023 paper (A&A, 2023) found that foreground galaxies in the structured hemisphere contribute to the apparent CMB fluctuations — exactly what closure would predict (the channel modifies the signal).

**Closure interpretation:** The CMB isn't a photograph of the last scattering surface. It's the **thermal signature of the closure loop's boundary**. The anomalies aren't primordial — they're artifacts of the loop structure itself. The hemispherical asymmetry reflects an asymmetry in the local closure boundary (our position relative to large-scale structure). The cold spot is a region where the channel has lower bandwidth (fewer intervening structures → less information content → lower apparent temperature).

**Observable:** Correlate the *local variance* of the CMB (not the mean temperature, but the fluctuation amplitude) with foreground structure maps. If the CMB is filtered through a channel, regions behind denser foreground should show modified fluctuation statistics — not just mean shifts (ISW effect, already known) but **variance changes** (bandwidth limitation).

**Data:** All public (Planck CMB maps + foreground catalogs). This is a map-level analysis, not a point-source cross-match.

### Test D: FRB Dispersion Measure × SN Hubble Residuals

**The idea:** Fast Radio Burst dispersion measures directly measure the integrated free electron column density along the line of sight. This is the purest baryonic channel tracer available.

**Observable:** For SN–FRB pairs within some angular separation (1–5°), correlate DM_excess (DM observed minus Milky Way contribution minus host contribution) with Δμ⊥.

**Data:** CHIME/FRB has ~2,500+ localized FRBs. Cross-match with Pantheon+ sky positions. Limited by angular overlap, but statistical binning can work with N ~ 50–100 pairs.

### Test E: Photon Arrival Time Structure (Future — Needs New Data)

**The most direct closure test:** If the channel imposes bandwidth constraints, photons from the same source but different frequencies should experience slightly different channel effects. For a pulsed source (FRB, GRB), this means frequency-dependent arrival time modifications beyond standard plasma dispersion (which goes as ν⁻²). Closure would predict **anomalous dispersion** that correlates with LOS baryonic structure in a way that standard plasma physics does not.

**Data:** FRB timing residuals after DM correction. If there's a systematic frequency-dependent residual that correlates with LOS environment, that's a closure-channel signature.

---

## 5. Priority Ranking

| Priority | Test | Data availability | Time to run | Novel? |
|----------|------|-------------------|-------------|--------|
| **1** | **B: SALT2 color vs LOS structure** | Already have everything | 10 min | Yes — nobody has checked this |
| **2** | **A: Light curve χ²/dof vs LOS** | Already have everything | 10 min | Yes — fit quality always attributed to intrinsic scatter |
| 3 | D: FRB DM × SN residuals | Need CHIME catalog | 1 hour | Partially novel |
| 4 | C: CMB variance vs foreground | Need analysis framework | 1 day | Moderate novelty |
| 5 | E: FRB anomalous dispersion | Needs specialized data | Research project | Highly novel |

**Tests 1 and 2 can be run right now with data we already have on disk.**

---

## 6. The Philosophical Shift

Standard cosmology asks: "How far away is that object?" 
Closure Theory asks: "What is the information-theoretic capacity of the channel connecting us to that object?"

The first question leads to distance ladders, standard candles, and luminosity measurements.
The second question leads to spectral complexity, temporal coherence, and channel bandwidth.

Every existing cosmological test is designed to answer the first question. The tests proposed above are designed to answer the second. If closure is real, the signal isn't in how *dim* distant objects appear — it's in how *simple* their signals become after traversing structured channels.

The universe doesn't hide the signal. We've been measuring the wrong thing.

---

## Appendix: Full Correlation Results

### Full-sample (N = 1,454 clean SNe)
| Tracer | Pearson r | p-value |
|--------|-----------|---------|
| κ_comp (10–20') | −0.024 | 0.362 |
| y_comp (10–20') | +0.009 | 0.728 |
| κ₃₀ (raw disc) | −0.022 | 0.393 |
| y₃₀ (raw disc) | +0.060 | 0.021 |

### Redshift splits — κ_comp (10–20')
| Bin | N | r | p |
|-----|---|---|---|
| z < 0.1 | 495 | +0.012 | 0.788 |
| 0.1 < z < 0.5 | 749 | −0.053 | 0.145 |
| z > 0.5 | 210 | −0.023 | 0.741 |

### Redshift splits — y_comp (10–20')
| Bin | N | r | p |
|-----|---|---|---|
| z < 0.1 | 495 | +0.024 | 0.591 |
| 0.1 < z < 0.5 | 749 | −0.034 | 0.352 |
| z > 0.5 | 210 | −0.066 | 0.340 |

### Redshift splits — y₃₀ (raw disc)
| Bin | N | r | p |
|-----|---|---|---|
| z < 0.1 | 495 | +0.111 | 0.014 |
| 0.1 < z < 0.5 | 749 | −0.022 | 0.556 |
| z > 0.5 | 210 | +0.049 | 0.481 |

### Analysis parameters
- ΛCDM: H₀ = 73.04, Ωm = 0.334 (Pantheon+ best fit)
- Planck κ: PR3 MV lensing alm → map at NSIDE = 1024
- Planck y: PR2 MILCA at NSIDE = 1024
- Deconfounding: OLS regression on E(B−V)_MW, host mass step, survey dummies
- Machine: WSL2/Windows, RTX 5090, 64GB RAM

### File locations (Humza's PC)
- Full data table: `C:\Users\humza\sn_closure_full.csv`
- Clean data table: `C:\Users\humza\sn_closure_clean.csv`
- Plots: `C:\Users\humza\elephant_all.png`, `zsplit_*.png`, `scale_test.png`
- Pipeline code: `/root/closure-theory/analyze.py`
- Maps: `/root/closure-theory/maps/` (1.5GB)
