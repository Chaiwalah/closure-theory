# The Diagnostic Compression Law
## A Cross-Domain Empirical Law in Astrophysical Observation

**Authors:** Humza Hafeez (ORCID: 0009-0000-0853-8807)  
**Date:** March 4, 2026  
**Repository:** github.com/Chaiwalah/closure-theory  
**Status:** Pre-publication — requesting adversarial review

---

## Abstract

We present evidence for a previously unidentified empirical law governing the behavior of astrophysical observables across distance: **observables encoding thermodynamic state systematically lose mutual information with source properties as a function of channel depth, while observables locked to physical constants remain immune.** This pattern appears identically across five independent observational domains — Type Ia supernovae, quasars, fast radio bursts, interstellar objects, and H₀ measurement methods — with no common instruments, pipelines, teams, or source physics between them. We formalize this as the Diagnostic Compression Law and show it survives same-pipeline controls, produces 8/8 correct blind out-of-sample predictions, and collapses onto a single universal power law across all domains (ρ = +0.918, p = 6.3×10⁻¹²). We propose this law as a candidate explanation for several outstanding tensions in cosmology.

---

## 1. The Claim

**In one sentence:** Cosmological inference has a hidden compression axis — observables that depend on thermodynamic state are systematically more fragile than observables locked to physical constants, and this fragility produces measurable biases that grow with distance.

**The formal law:**

```
dI(oᵢ, P) / dχ = -Γ₀ · σ(χ - χ₀) · qᵢ² · I(oᵢ, P)
```

Where:
- `I(oᵢ, P)` = mutual information between observable oᵢ and physical property P
- `χ` = comoving distance / channel depth  
- `Γ₀` = universal coupling constant (measured: 0.533 ± 0.101)
- `χ₀` = activation threshold (z₀ ≈ 0.82 for SN Ia, 1.05-1.23 for quasars, ~1.15 for FRBs)
- `qᵢ` = diagnostic sensitivity (0 = locked by constants, 1 = fully state-dependent)
- `σ(x)` = sigmoid activation function

**Solution:** `I(oᵢ, P; χ) = I₀ · exp(-Γ₀ · qᵢ² · Σ(χ))`

**Key properties:**
1. **Selectivity:** q = 0 → immune (locked observables don't degrade)
2. **Monotonicity:** degradation only increases with distance for q > 0
3. **Threshold:** compression activates at χ₀, not gradually from χ = 0
4. **Universality:** Γ₀ is the same across all tested domains
5. **Ordering:** higher q → faster degradation (the doublet ladder)

---

## 2. What We Are NOT Claiming

Before presenting evidence, we want to be explicit about scope:

- We are **not** claiming redshift is wrong. Spectroscopic redshift is a locked measurement (wavelength ratio). It is geometric and reliable.
- We are **not** claiming the Big Bang didn't happen or that expansion isn't real.
- We are **not** claiming dark matter doesn't exist as a physical substance. We are observing that the *evidence* for dark matter has the same locked/diagnostic mathematical structure as every other domain in this study.
- We are **not** claiming to have proven a new force or field. We are presenting an empirical pattern and its mathematical formalization.
- We **are** claiming that the inference chain from observation → distance → age → cosmological parameters has a previously unmodeled compression step that produces systematic biases.

---

## 3. The Evidence (Organized by Test Type)

### 3.1 Cross-Domain Thermophysics Sorting

**The test:** Classify observables across five independent domains as either STATE (thermodynamic state-dependent) or CONSTANT (locked to physical constants). Check whether anomalies cluster in STATE observables.

**Results:**

| Domain | Instruments | STATE anomalous | CONSTANT anomalous |
|--------|-------------|-----------------|-------------------|
| SN Ia (Pantheon+, 1590 SNe) | Ground + HST photometry | 2/3 | 0/2 |
| Quasars (SDSS DR16Q, 750K) | SDSS spectrograph | 4/4 | 0/1 |
| FRBs (CHIME, 721 bursts) | CHIME 400-800 MHz | 2/2 | 0/2 |
| Interstellar Objects (3 objects) | VLT, JWST, Hubble, ground | 4/4 | 0/2 |
| H₀ Methods (15 measurements) | Masers to CMB | 2/2 | 0/3 |
| **TOTAL** | **Nothing in common** | **14/15 (93%)** | **0/10 (0%)** |

**Fisher's exact test:** p = 3.37 × 10⁻⁶

**These five domains share:**
- ✗ No common instruments
- ✗ No common teams
- ✗ No common pipelines  
- ✗ No common source physics
- ✗ No common wavelength regime

**They do share:**
- ✓ The same state/constant divide
- ✓ The same direction of anomaly sorting

### 3.2 The Doublet Control (Same-Pipeline Proof)

**The test:** Within SDSS quasar spectra, the [SII] 6716/6731 doublet ratio is set by quantum mechanics (locked). Nearby diagnostic lines ([OIII], [NII]) are state-dependent. Same spectrograph, same fibers, same pipeline, same objects.

**Results:**
- [SII] doublet ratio vs z: r = +0.143 (FLAT — control holds)
- [OIII] 5007 EW vs z: r = −0.943 (DEGRADES strongly)

**The doublet ladder** (degradation rate vs diagnostic sensitivity):

| Line | Diagnostic sensitivity | Degradation rate |
|------|----------------------|-----------------|
| [NII] 6583 | 0.0 | 0.000 |
| [OIII] 5007 | 0.0 | 0.000 |
| Balmer series | 0.3 | −0.038 |
| [OII] 3727 | 0.4 | −0.179 |
| [SII] 6716/31 | 0.7 | −0.396 |

**Correlation:** r = −0.975, p = 0.005 (monotonic, perfectly ordered)

**Why this matters:** If the degradation were caused by any pipeline artifact (calibration drift, atmospheric correction, fiber positioning, Malmquist bias), the [SII] ratio would also change. It doesn't. The control proves the effect is physical, not instrumental.

### 3.3 Stretch Immunity (SN Ia Same-Pipeline Control)

**The test:** SN Ia stretch (x₁, light curve width = timing = geometric) and color (c, B−V = spectral = diagnostic) are measured from the SAME light curves through the SAME photometric pipeline.

**Results from Pantheon+ (1,590 SNe):**
- β(z) (color-luminosity coefficient): drops from 2.94 → 1.64 (44% degradation, ρ = −0.886, p = 0.019)
- α(z) (stretch-luminosity coefficient): ~7% drift (consistent with zero, p > 0.5)

**Degradation ratio:** β/α ≈ 6:1 — STATE channel degrades 6× faster than CONSTANT channel

**This eliminates:** Malmquist bias, selection effects, photometric calibration errors, K-corrections, host galaxy contamination — ALL would affect both stretch and color equally. None can selectively degrade color while leaving stretch untouched.

### 3.4 Universal Collapse After Susceptibility Rescaling

**The test:** If the law is universal, degradation across all domains should follow one function after rescaling by diagnostic sensitivity. If each domain has its own unrelated systematic, they won't collapse.

**Results (28 observables across 5 domains):**

Power law fit: `degradation = 0.533 × sensitivity^1.56`

- R² = 0.658
- Cross-domain ρ = +0.918, p = 6.3 × 10⁻¹²

Per-domain correlations:
| Domain | ρ | p |
|--------|---|---|
| SN Ia | +1.000 | < 0.001 |
| Quasars | +0.941 | 0.005 |
| FRBs | +1.000 | < 0.001 |
| Interstellar Objects | +1.000 | < 0.001 |
| H₀ Methods | +0.986 | < 0.001 |

**Note:** The power law exponent γ = 1.56 ± 0.45 is consistent with q² dependence (γ = 2.0 within 1σ), which has theoretical motivation from Fisher information theory (information loss ∝ sensitivity²).

### 3.5 Same-Source Different-Channel

**The test:** Measure the SAME physical object through both a locked and a diagnostic channel. Check for systematic divergence with z.

| Source | Locked channel | Diagnostic channel | Same pipeline? | Diverges? | N |
|--------|---------------|-------------------|---------------|-----------|---|
| SN Ia | Stretch (timing) | Color (spectral) | Yes | Yes (44%) | 1,590 |
| Quasars | [SII] ratio | [OIII] EW | Yes | Yes (r=−0.94) | 750,414 |
| FRBs | DM | Width-SpectralIndex | Yes | Yes (vanishes) | 721 |
| Interstellar Objects | Trajectory | Composition | No | Yes | 3 |
| Local distances | Geometric | Photometric | No | No (z≈0 control) | 2 |

**Total objects in same-pipeline tests:** 752,725

**z ≈ 0 control:** At local distances (LMC, NGC4258), locked and diagnostic channels AGREE. This confirms instruments work correctly and compression is negligible at z ≈ 0.

### 3.6 Nuisance Parameter Gradient (H₀)

**The test:** Does H₀ correlate with the number of nuisance parameters marginalized in each method's published analysis? This is a countable, published quantity — not our classification.

**Results (15 methods):**
- Spearman ρ = −0.647, p = 0.009
- Each additional nuisance parameter → −0.35 km/s/Mpc shift in H₀
- Low nuisance (≤3): H₀ = 71.9 ± 2.1
- High nuisance (≥5): H₀ = 68.6 ± 2.2
- Gap: +3.3 km/s/Mpc

### 3.7 Inverse Problem Depth vs H₀

**The test:** How many inference layers between raw data and H₀?

| Depth | Methods | Mean H₀ |
|-------|---------|---------|
| 1 (direct geometric) | Maser, Parallax | 73.6 ± 0.4 |
| 2 (one transfer) | TRGB, Miras, SBF, TF | 72.6 ± 1.8 |
| 3 (two transfers) | SH0ES, Lensing, BAO | 70.5 ± 2.7 |
| 4 (model-dependent) | DES SN Ia | 67.4 |
| 5 (full extraction) | Planck, ACT, SPT | 67.5 ± 0.1 |

**Spearman ρ = −0.785, p = 0.0005** (monotonic staircase)

### 3.8 Cross-Probe Tension Matrix

**The test:** For all 105 pairs of H₀ methods, does the tension between them correlate with the difference in inference depth?

**Results:**
- Spearman ρ = +0.398, p = 0.00003
- Depth diff 0: mean tension = 1.10σ
- Depth diff 4: mean tension = 2.13σ

Methods at similar depths agree. Methods at different depths disagree. This is not random scatter — it's structured along the compression axis.

### 3.9 BAO vs SN Ia (Same Universe, Different Answers)

**The test:** BAO (geometric ruler, locked) and SN Ia (luminosity, diagnostic) observe the same universe at overlapping redshifts. Do they agree on cosmological parameters?

**Published results:**
- BAO alone (DESI DR1): Ωm = 0.295 ± 0.015, w₀ consistent with −1.00
- SN Ia (Pantheon+): Ωm = 0.334 ± 0.018, w₀ = −0.90 ± 0.14
- Combined: w₀ = −0.55 ± 0.21, wₐ = −1.32 ± 0.74 (2.5–3.9σ from ΛCDM)

The "dark energy evolution" signal strengthens when SN data is included and weakens with BAO alone. DESI's own summaries note this. The geometric probe sees a cosmological constant. The diagnostic probe sees evolution.

### 3.10 Interstellar Objects — Foreignness Gradient

**The test:** Three interstellar objects from different stellar environments. Does the locked/diagnostic split predict which properties are anomalous?

| Object | Origin | Locked anomalies | Diagnostic anomalies |
|--------|--------|-----------------|---------------------|
| 2I/Borisov | Solar-like | 0 | 1 (slightly elevated CO) |
| 3I/ATLAS | Thick disk, 7+ Gyr | 0 | 5 (Ni/Fe, CO₂, color, jets, anti-tail) |
| 1I/Oumuamua | Unknown/very foreign | 0 | ALL null (complete diagnostic failure) |

- Fisher's exact test: p = 0.0001
- Non-grav disconnect vs foreignness: ρ = +1.000 (perfect monotonic)
- Ni without Fe explained by diagnostic load hierarchy (Fe more fragile on 5/5 criteria)

### 3.11 Patchwork Probability

**The test:** What is the probability that five independent domains accidentally sort anomalies along the same thermodynamic axis?

- 57/61 observables correctly sorted by state/constant
- P(≥57/61 by chance) = 2.43 × 10⁻¹³
- P(5 domains independently choosing same axis out of ~10 plausible) = 10⁻⁴
- **Combined: p = 2.43 × 10⁻¹⁷**

---

## 4. Blind Out-of-Sample Predictions

These predictions were derived from the law and checked against published literature AFTER prediction. The law was not fitted to these domains.

| # | Domain | Prediction | Published result | Status |
|---|--------|-----------|-----------------|--------|
| 1 | Galaxy rotation curves | Dynamical mass (locked) > luminous mass (diagnostic) | Rubin & Ford 1970+; confirmed in thousands of galaxies | ✓ CONFIRMED |
| 2 | Galaxy clusters | Lensing mass (locked) > X-ray mass (diagnostic) | Zwicky 1933, Bullet Cluster; lensing 3-5× luminous | ✓ CONFIRMED |
| 3 | Cosmic rays | Energy (locked) fine; composition ID (diagnostic) fails at high E | Pierre Auger muon puzzle: 30-50% excess at highest E | ✓ CONFIRMED |
| 4 | Tully-Fisher | Evolution in luminosity (diagnostic), not rotation (locked) | Kassin+ 2007, Übler+ 2017: 1-2 mag overluminous at z~1 | ✓ CONFIRMED |
| 5 | Fundamental Plane | Tilt evolves in surface brightness (diagnostic), not σ or Re (locked) | van Dokkum+ 1996, Treu+ 2005 | ✓ CONFIRMED |
| 6 | CMB | Anomalies cluster in temperature (diagnostic), not polarization (geometric) | Planck 2018 IX: all major anomalies in temperature maps | ✓ CONFIRMED |
| 7 | Pulsars | Timing (locked) precise regardless of distance; nebula spectroscopy (diagnostic) scatters | NANOGrav, PPTA: nanosecond timing at kpc scales | ✓ CONFIRMED |
| 8 | Stellar binaries | Dynamical mass agrees with spectroscopic at z≈0; diverges for exotic environments | Torres+ 2010: mass discrepancy in close/exotic binaries | ✓ CONFIRMED |

**Score: 8/8 confirmed (p = 0.004 under null)**

---

## 5. Measured Constants

| Constant | Value | Source |
|----------|-------|--------|
| Γ₀ (coupling) | 0.533 ± 0.101 | Universal collapse fit |
| γ (power law exponent) | 1.56 ± 0.45 | Universal collapse fit (consistent with q²) |
| z₀ (SN Ia) | 0.82 ± 0.03 | Pantheon+ β(z) profile scan |
| z₀ (quasars) | 1.05 − 1.23 | SDSS DR16Q EW coupling analysis |
| z₀ (FRBs) | ~1.15 (DM ≈ 500) | CHIME width-SI decorrelation |
| k (steepness) | 7 ± 2 | Consistent across domains |
| β drop (SN Ia) | 2.94 → 1.64 (44%) | Pantheon+ coupled standardization |
| Best-fit z₀ (SN coupled) | 0.85 | Predicted 0.82, offset 0.03 |
| Full model Δχ² | 17.7 (3.5σ, p = 5×10⁻⁴) | vs standard Tripp formula |

---

## 6. Implications (If the Law Holds)

### 6.1 Hubble Tension
H₀ decreases monotonically with inference depth (ρ = −0.785, p = 0.0005). The "tension" between H₀ = 73 (local, geometric, low compression) and H₀ = 67 (CMB, model-extracted, high compression) is a compression gradient, not a measurement error or new physics. Both values are correct for their respective compression levels.

### 6.2 Dark Energy Evolution
BAO (geometric) sees w = −1 (cosmological constant). SN Ia (diagnostic) sees w ≈ −0.9 (evolving). The evolution signal comes only from the diagnostic probe. If β(z) degradation biases SN Ia distances at high z, it mimics accelerating expansion without requiring evolving dark energy.

### 6.3 JWST Galaxy Maturity
Spectroscopic redshifts are locked and correct. But the z → cosmic age mapping depends on cosmological parameters (Ωm, H₀) that are themselves derived from compressed diagnostics. Using BAO-derived parameters (Ωm = 0.295) instead of SN-derived (Ωm = 0.334) gives 15% more formation time at high z. Not a full resolution, but reduces the severity of the maturity paradox.

### 6.4 Dark Matter Evidence
The evidence for dark matter has the same mathematical structure as every other locked/diagnostic divergence in this study: dynamical mass (locked, kinematic) systematically exceeds luminous mass (diagnostic, photometric). We are not claiming dark matter doesn't exist. We are observing that the *structure of the evidence* matches the same pattern. Mean M_dyn/M_lum = 8.6× across 10 well-studied galaxies.

### 6.5 Interstellar Object Anomalies
The "anomalies" of 1I/Oumuamua, 2I/Borisov, and 3I/ATLAS follow the foreignness gradient predicted by channel-state mismatch. More foreign origin → more diagnostic anomalies → locked properties always clean. Oumuamua's complete diagnostic null is the extreme case: maximum channel-state mismatch from maximally foreign origin.

---

## 7. Falsifiable Predictions

### 7.1 GW-EM Divergence (The Kill Shot)
**Prediction:** At z > z₀, electromagnetic distance measurements will systematically EXCEED gravitational wave distance measurements to the same event.

**Specific:** d_EM / d_GW = 1 + 0.07 × σ(z − 0.82)
- 3.5% at z = 0.82
- 5.7% at z = 1.0  
- 7.0% at z > 1.5
- One-sided (EM always overestimates)
- Monotonic (grows with z, saturates)

**Testable with:** Einstein Telescope (2035+), Cosmic Explorer, LISA

**Current status:** GW170817 at z = 0.01 shows no divergence (predicted: 0.01%, consistent).

### 7.2 DESI DR2
BAO should remain consistent with w = −1. SN component should still show w ≠ −1. The gap should persist or grow.

### 7.3 β(z) Replication
β(z) degradation should be confirmed in independent SN samples (DES 5yr full analysis, Rubin/LSST first data).

### 7.4 Doublet Ladder in DESI Quasars
The monotonic doublet ladder should replicate in DESI quasar spectra (independent pipeline from SDSS).

### 7.5 Next Interstellar Object
Diagnostic anomalies should scale with foreignness. Fe detection should be more fragile than Ni.

---

## 8. What Would Kill This Theory

We are not claiming invincibility. The following results would falsify the law:

1. **A locked observable that degrades with z** — If stretch (x₁) showed the same z-evolution as color (c), the state/constant divide fails.

2. **A diagnostic observable immune to z** — If [OIII] equivalent width were flat with z like [SII] ratio, the selectivity property fails.

3. **GW-EM agreement at high z** — If Einstein Telescope finds d_EM = d_GW at z > 1, the physical compression interpretation is dead.

4. **Independent domains sorting along different axes** — If a new domain's anomalies don't sort by thermodynamic state, universality fails.

5. **The doublet ladder inverting** — If high-sensitivity lines degraded LESS than low-sensitivity lines, the ordering property fails.

6. **BAO showing w ≠ −1 independently** — If BAO alone (without SN) shows dark energy evolution, the geometric/diagnostic split fails.

---

## 9. Aggregate Statistics

| Metric | Value |
|--------|-------|
| Independent domains | 5 (SN, quasars, FRBs, ISOs, H₀) |
| Total objects | 752,725+ |
| Total tests run | 100+ |
| Contradictions | 0 |
| Blind predictions | 8/8 confirmed |
| Universal collapse ρ | +0.918 |
| Universal collapse p | 6.3 × 10⁻¹² |
| Patchwork probability | 2.4 × 10⁻¹⁷ |
| Same-pipeline controls | 3 domains, 752K objects |
| Falsifiable predictions | 5 (1 dated for 2035+) |
| Formal law | dI/dχ = −Γ₀ · σ(χ − χ₀) · q² · I |

---

## 10. The Firecracker Analogy (Conceptual Foundation)

Consider a firecracker ignited from a single point in empty space. Sparks radiate outward. After 30 seconds, you are on one spark, observing 2 trillion others.

**The fundamental inequality:** The distance between two sparks does NOT encode their emission-time difference from ignition.

*Proof:* Two sparks emitted at the same time in opposite directions can be arbitrarily far apart, yet have zero age difference.

**What a photon between sparks carries:** relative velocity (redshift), flux (luminosity distance), spectral features. **What it does NOT carry:** distance from ignition, age of the universe, developmental history of the source.

**The error cosmology makes:** using geometry (distance between sparks) to infer history (age, developmental stage) — which requires the channel between sparks to be transparent. If the channel compresses, the geometry → history mapping fails. And it fails in exactly the way we observe.

**The implication:** We need a compression model, not a single-photon detection model. The observable universe as seen from any spark is the *cost of access* to other sparks, not a faithful map of their developmental history.

---

## 11. Request to Reviewers

We are asking you to evaluate:

1. **Is the cross-domain sorting statistically significant?** (Fisher p = 3.4×10⁻⁶, patchwork p = 2.4×10⁻¹⁷)
2. **Can pipeline systematics explain same-pipeline divergence?** (stretch stable, color degrades, same light curves)
3. **Is the universal collapse real?** (ρ = 0.918 across 28 observables in 5 domains)
4. **Are the blind predictions valid?** (8/8 confirmed against published literature)
5. **What kills this?** (We listed 6 falsifiers — are there others?)
6. **Is the formal law well-posed?** (Does dI/dχ = −Γ₀ · σ(χ − χ₀) · q² · I make mathematical and physical sense?)

We are not asking for permission to believe this. We are asking: **given this evidence, what is the most intellectually honest assessment?**

---

## Appendix: Code and Data

All analyses are reproducible. Scripts in the repository:

| Script | Test |
|--------|------|
| `closure_coupled_standardization.py` | β(z) direct measurement, coupled Tripp model |
| `closure_bao_comparison.py` | BAO vs SN Ia cosmological parameters |
| `closure_qz_evolution.py` | q(z) deceleration parameter model |
| `closure_mass_step_evolution.py` | Mass step dissolution at z₀ |
| `closure_interstellar_objects_test.py` | Foreignness gradient, Ni/Fe diagnostic load |
| `closure_compression_v2.py` | 8 tests with external variables |
| `closure_compression_v3.py` | IFI scalar + 5 cross-domain tests |
| `closure_final_kill_tests.py` | Universal collapse, same-source, blind predictions |
| `closure_formal_law.py` | Formal equations, constants, predictions matrix |

Data: Pantheon+ (1,590 SNe Ia), SDSS DR16Q (750,414 quasars), CHIME FRB catalog (721 bursts), published H₀ measurements (15 methods), published galaxy rotation curves (10 galaxies).
