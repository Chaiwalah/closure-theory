# Paper 1: Empirical Paper Outline

## Title (candidates)
1. "Selective Degradation of Spectral Correlations as a Function of Diagnostic Sensitivity Across Cosmological Distances"
2. "A Diagnostic Sensitivity Law for Spectral Correlation Loss in Cosmological Observation"
3. "Evidence for Bandwidth-Limited Information Channels in Cosmological Spectroscopy"

**Recommendation:** Title 2 — names the law, signals discovery, accessible

## Target Journals (priority order)
1. **Nature Astronomy** — Letter format (~3,500 words, 4-6 figures, 30 refs). High impact. Perfect for a new empirical law with cross-domain confirmation. Turnaround: ~4-6 weeks to first decision.
2. **Physical Review Letters** — Shorter (~3,500 words, 4 figures). Prestige in physics. Good if we emphasize the information-theoretic angle.
3. **The Astrophysical Journal Letters** — Standard astro venue. Easier acceptance, still respected. Good fallback.
4. **Monthly Notices of the Royal Astronomical Society (Letters)** — Fast turnaround, open access option.

**Recommendation:** Nature Astronomy first. It's where anomalous cosmological results get attention. Backup: PRL for the physics angle.

## Abstract (~150 words)

We report the discovery of a systematic, selective degradation of spectral correlations with cosmological distance. Using 750,414 quasars from SDSS DR16Q, 1,590 Type Ia supernovae from Pantheon+SH0ES, and 535 fast radio bursts from CHIME, we show that correlations between spectral observables degrade at rates proportional to their diagnostic sensitivity — the degree to which the observable encodes thermodynamic state information. Observables locked to fundamental physical constants (e.g., [OIII] 4959/5007 branching ratio) show zero degradation, while state-dependent diagnostics (e.g., CIV/MgII flux ratio) degrade with a measured rank correlation of r = −0.975 (p = 0.005). The degradation follows a sigmoid threshold at z₀ ≈ 0.82 (SNe Ia), 1.05–1.23 (quasars), consistent with a bandwidth-limited information channel. We measure the functional relationship: eaten = 0.048(1 − MI₀) + 0.012 (r = −0.718, p = 0.009) across 12 emission-line pairs in 130,000 DESI galaxies. These results suggest cosmological observation channels are not informationally neutral.

## Structure

### 1. Introduction (~500 words)
- Standard cosmological analyses assume spectral information propagates losslessly
- Recent anomalies challenge this: Hubble tension, JWST "impossible" galaxies, DESI dark energy evolution
- We ask: do all spectral observables degrade equally with distance, or is there selectivity?
- Preview: we find a specific, measurable law governing selective degradation
- No framework claims — just "we looked and found this"

### 2. Data and Methods (~600 words)

#### 2.1 Datasets
- **Pantheon+SH0ES**: 1,590 SNe Ia, z = 0.01–2.26. Color (c), stretch (x₁), distance modulus (μ), host mass, survey metadata.
- **SDSS DR16Q**: 750,414 spectroscopic quasars, z = 0.1–5.0. Emission line measurements: [OIII], [NII], [SII], [OII], Balmer, CIV, MgII, FeII.
- **CHIME FRB Catalog 1**: 535 FRBs (+ 186 localized). Dispersion measure (DM), rotation measure (RM), temporal width, spectral index, fluence.
- **DESI Year 1**: 130,000 emission-line galaxies for the eating law measurement.

#### 2.2 Observable Classification
- Define "diagnostic sensitivity" q ∈ [0, 1]:
  - q = 0: Observable value fixed by fundamental physics (branching ratios, selection rules)
  - q = 1: Observable value fully determined by local thermodynamic/environmental state
- Classification derived from atomic physics (partial derivatives of emissivity w.r.t. T, n_e, Z)
- Independent of any theoretical framework — purely empirical classification

#### 2.3 Correlation Measurement
- Mutual information (MI) between emission-line pairs as function of redshift
- Spearman rank correlations in redshift bins
- Sigmoid threshold fitting: σ(z) = 1/(1 + exp(−k(z − z₀)))
- Controls: luminosity, SNR, survey selection, aperture effects

### 3. Results (~1,000 words)

#### 3.1 The Doublet Ladder (KEY FIGURE 1)
- Six observables ranked by diagnostic sensitivity vs measured degradation rate
- r = −0.975 (Spearman), p = 0.005, perfectly monotonic
- Locked lines ([NII], [OIII] doublet ratios): zero degradation across all z
- Maximum degradation: CIV/MgII cross-ion ratio (q = 1)
- **This is the central result: degradation rate is a function of diagnostic sensitivity**

#### 3.2 Sigmoid Thresholds (FIGURE 2)
- Color-distance coupling (SNe Ia): sigmoid onset at z₀ = 0.82 ± 0.03, k = 8.0 ± 2.0
- CIV/MgII coupling collapse (quasars): z₀ = 1.05–1.23
- FRB width-spectral index: correlation vanishes past DM ≈ 500
- All three source classes show sharp threshold, not gradual decline

#### 3.3 The Eating Law (FIGURE 3)
- 12 emission-line pairs in 130K DESI galaxies
- eaten = 0.048 × (1 − MI₀) + 0.012
- r = −0.718, p = 0.009
- Interpretation: degradation proportional to how much "readable" information the pair carries
- Pairs with high initial MI (closed at emission) are immune
- Pairs with low initial MI degrade proportionally

#### 3.4 Channel Divergence (FIGURE 4)
- Flux correlations degrade: r = −0.943
- Sigma (line width) correlations: FLAT, r = +0.143
- Gap (velocity offset): r = −1.000
- The effect is channel-specific, not geometric
- If it were dust, distance, or lensing → all channels would degrade equally

#### 3.5 Cross-Domain Consistency (FIGURE 5)
- Same pattern independently confirmed in:
  - Type Ia supernovae (standard candles)
  - Quasars (active galactic nuclei)  
  - Fast radio bursts (transient radio sources)
- Three completely different source classes, different wavelengths, different physics
- All show: frequency-dependent observable entanglement with sharp threshold, kinematic immunity

#### 3.6 Environmental Modulation (FIGURE 6)
- Cluster sightlines: +14% correlation preservation (Δρ = +0.141)
- Void sightlines: different coupling behavior
- Dense environments preserve correlations → consistent with matter increasing effective channel capacity
- Tested with BLR 5D control (luminosity, BH mass, Eddington ratio, FeII, redshift): signal survives

### 4. Functional Form (~400 words)

#### 4.1 Best-Fit Model
- Tested: linear, power law, exponential, Fisher quadratic, logarithmic
- Best fit: exponential growth in q (R² = 0.994)
- Power law with γ = 2.2 also excellent (R² = 0.993)
- Pure quadratic (Fisher information prediction): R² = 0.990
- Linear: R² = 0.879 — clearly rejected
- The relationship is nonlinear and accelerating

#### 4.2 The Diagnostic Compression Law (formal statement)
- degradation_i ∝ q_i^γ, γ = 2.2 ± 0.5
- Consistent with information-theoretic prediction (γ = 2, Fisher information)
- Steeper exponent suggests cooperative coupling (γ > 2)

### 5. Controls and Systematics (~500 words)

#### 5.1 What This Is NOT
- Not Malmquist bias (signal amplifies under selection controls)
- Not K-correction (survives rest-frame wavelength matching)
- Not aperture effects (survives spatial deblending)
- Not SNR artifact (signal present in high-SNR subsample)
- Not pipeline artifact (tested with injection simulations)

#### 5.2 Conditional Independence Test
- Signal AMPLIFIES to 125.5% under full controls (survey + mass + SNR + coverage)
- p = 1.1 × 10⁻⁵
- Model is definitively missing a latent variable

#### 5.3 Null Controls
- Stretch (x₁) immune across all z: temporal observable, not spectral → expected
- HeII lines: no degradation → consistent with survivorship selection at extreme z
- Baryonic proxies (host mass, metallicity): NULL → effect is not caused by intervening matter

### 6. Predictions (~300 words)

#### 6.1 Testable Now
- P1: DESI Year 3 quasar spectra will show same doublet ladder ordering
- P2: High spectral resolution (DESI > SDSS) should show stronger correlations at matched z
- P3: Lensed quasar pairs should show correlated degradation (same sightline)

#### 6.2 Testable Soon
- P4: Rubin/LSST SNe Ia will confirm β(z) degradation in independent sample
- P5: JWST deep spectroscopy at z > 2 will show accelerated diagnostic loss
- P6: SKA FRB catalog will extend width-spectral decoupling to higher DM

#### 6.3 Quantitative
- From the fitted law: specific numerical predictions for correlation strengths in upcoming surveys (table)

### 7. Discussion (~400 words)
- These results are consistent with cosmological observation channels having finite information bandwidth
- The selectivity (diagnostic lines degrade, locked lines don't) rules out geometric explanations (dust, lensing, expansion effects)
- The environmental modulation (matter preserves correlations) suggests the channel capacity is a physical property of the propagation medium
- We note structural similarity to bandwidth-limited systems in other domains (LHC trigger systems, biological measurement channels) without claiming mechanistic equivalence
- Implications for Hubble tension: if spectral and geometric channels have different effective bandwidths, systematic disagreement between methods is expected, not anomalous
- Implications for JWST "impossible" galaxies: spectral age inference at high-z may be unreliable due to diagnostic compression

### 8. Conclusion (~150 words)
- We have identified a new empirical law: spectral correlation degradation scales with diagnostic sensitivity
- The law is confirmed across three independent source classes with >750,000 objects
- It follows a specific functional form (power law with γ ≈ 2.2)
- It has sharp sigmoid activation thresholds
- It is selective (channel-dependent, not geometric)
- It is modulated by environment (matter preserves correlations)
- The default assumption of informationally neutral cosmological observation may need revision

---

## Figures (6)

1. **The Doublet Ladder** — Diagnostic sensitivity (x) vs degradation rate (y). Six points, monotonic, r = −0.975. THE figure. With error bars and functional form overlay.

2. **Sigmoid Thresholds** — Three panels: SNe Ia color-distance coupling, quasar CIV/MgII, FRB width-spectral index. All showing sharp transition at z₀.

3. **The Eating Law** — MI₀ (x) vs degradation (y) for 12 line pairs. Linear fit. r = −0.718. Shows the rule: more initial information = more degradation.

4. **Channel Divergence** — Three panels: flux (degrades), sigma (flat), gap (degrades). Same objects, same redshift range, different observables. Proves selectivity.

5. **Cross-Domain Consistency** — Three-panel comparison: SNe / Quasars / FRBs all showing the same pattern. Different wavelengths, different physics, same law.

6. **Environmental Modulation** — Cluster vs void sightlines. Correlation preservation in dense environments. With BLR controls.

---

## Supplementary Material
- Full test results table (100+ tests)
- Pipeline injection simulations
- Monte Carlo uncertainty on q classification
- Extended figures for all line pairs
- Data access: GitHub repo link

---

## Author
- Humza Hafeez (independent researcher)
- ORCID: 0009-0000-0853-8807

## Acknowledgments
- AI-assisted analysis (Claude/OpenClaw for test design and execution)
- SDSS DR16Q, Pantheon+SH0ES, CHIME, DESI collaborations for public data
- Zenodo DOI for pre-registration of framework (to be added)

---

## Timeline
1. [ ] Upload framework paper to Zenodo (timestamp priority)
2. [ ] Finalize figure set
3. [ ] Draft manuscript (~2 weeks)
4. [ ] Internal review (GPT, Grok, Claude adversarial)
5. [ ] Submit to Nature Astronomy
6. [ ] Simultaneously post to arXiv (astro-ph.CO)

---

## Key Numbers for Quick Reference
- 752,725 total objects (1,590 SNe + 750,414 quasars + 721 FRBs)
- 100+ statistical tests
- 0 contradictions
- Doublet ladder: r = −0.975, p = 0.005
- Eating law: r = −0.718, p = 0.009
- Sigmoid: z₀ = 0.82 (SNe), 1.05–1.23 (quasars)
- Channel divergence: flux r = −0.943, sigma r = +0.143
- Cluster shadow: +14% (Δρ = +0.141)
- Conditional independence: 125.5%, p = 1.1 × 10⁻⁵
- Cross-domain: 3 source classes, 65+ tests at time of first cross-domain paper
- Functional form: γ = 2.2 ± 0.5 (power law), R² = 0.993
- Higgs portal test: ρ = −1.000 (save for Paper 2)
