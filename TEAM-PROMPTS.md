# TEAM PROMPTS — Closure Theory Next Phase
*Feb 27, 2026 — Copy-paste ready*

---

## 🔬 GEMINI — Plasma Lens Simulation

```
I need you to simulate photon propagation through an intermittent plasma lens network to test whether selective spectral degradation emerges naturally.

SETUP:
- Generate 10,000 "spectral lines," each with two properties:
  1. Intensity structure: a Gaussian with random substructure (2-5 subcomponents with random amplitudes) — this represents EW-like information (complex, spatially structured)
  2. Velocity centroid: a single number drawn from N(0, σ) — this represents FWHM-like information (simple, scalar)

- Define a "plasma lens" as a phase screen that applies wavelength-dependent phase shifts:
  φ(λ) = (e²/2πmₑc) × DM_lens / (c/λ)
  where DM_lens is drawn from LogNormal(mean=10, σ=5) pc/cm³

- Each "sightline" passes through N lenses, where N represents path length. Vary N from 10 to 200 (representing DM_total from ~100 to ~2000 pc/cm³)

- Between lenses, there's a gap. The gap has probability p_gap of being "empty" (no lens). Vary p_gap from 0 (fully connected network) to 0.8 (fragmented network).

MEASURE:
For each sightline:
1. Cross-correlate the output intensity structure with the input → get "EW fidelity"
2. Compare output velocity centroid to input → get "FWHM fidelity"
3. Compute Spearman correlation between EW fidelity and FWHM fidelity across all 10,000 lines

QUESTIONS:
1. Does the EW-FWHM correlation decay with increasing N (total DM)?
2. Is there a THRESHOLD in N where decorrelation accelerates?
3. Does the threshold depend on p_gap (network connectivity)?
4. Is the degradation wavelength-dependent? (Run at λ = 1200, 2800, 5000, 6700 Å)
5. If you introduce a "susceptibility" parameter S for each line (S = number of subcomponents / total complexity), does degradation rate scale with S?

USE: Standard IGM parameters. n_e = 2×10⁻⁷ cm⁻³, T = 10⁴ K.

This tests whether a physically motivated plasma lens network naturally produces:
- Selective degradation (EW damaged, FWHM preserved)
- A DM threshold
- Complexity-dependent vulnerability
- Wavelength dependence

Plot all results. If a threshold emerges, report the DM value and compare to 774 pc/cm³.
```

---

## 🧠 GPT — Derive Susceptibility P From Atomic Physics

```
I have an empirical result from 750,000 SDSS quasars that I can't explain from first principles. I need your help identifying the physical property that determines a spectral line's vulnerability to information degradation.

THE RESULT:
The correlation between emission line Equivalent Width (EW) and Full Width at Half Maximum (FWHM) decays with redshift. The decay rate is:

dr/dz = -0.10 × P

where P is an empirically measured "susceptibility." Here is the measured ranking:

| Line      | λ_rest (Å) | Type | Ionization | dr/dz    | P_empirical |
|-----------|-----------|------|------------|----------|-------------|
| Lyα       | 1216      | BLR  | H I        | -0.1042  | 1.000       |
| [SII]     | 6718      | NLR  | S II       | -0.0921  | 0.983       |
| CIV       | 1549      | BLR  | C IV       | -0.0554  | 0.930       |
| MgII      | 2798      | BLR  | Mg II      | -0.0170  | 0.875       |
| [OIII]    | 5007      | NLR  | O III      | -0.0094  | 0.864       |
| [OII]     | 3728      | NLR  | O II       | +0.0043  | 0.845       |
| CIII]     | 1909      | BLR  | C III      | +0.0417  | 0.791       |
| Hβ        | 4861      | BLR  | H I        | +0.2338  | 0.517       |
| Hα        | 6563      | BLR  | H I        | +0.5084  | 0.124       |
| [NII]     | 6585      | NLR  | N II       | +0.5951  | 0.000       |

WHAT DOESN'T EXPLAIN THE ORDERING:
- NOT wavelength: SII (6718Å, red) is #2 most vulnerable while Hα (6563Å, also red) is immune
- NOT ionization state: SII (low ionization) is vulnerable, OIII (high ionization) is not
- NOT BLR vs NLR: Both types appear at both ends
- NOT element: Both hydrogen lines appear at both ends (Lyα most vulnerable, Hα immune)

WHAT I NEED FROM YOU:
1. What PHYSICAL PROPERTY of an emission line could produce this specific ordering?
2. Consider: critical density, number of fine-structure sublevels, sensitivity to electron temperature vs density, transition probability (Einstein A coefficient), optical depth effects, line formation physics
3. The key constraint: Lyα and [SII] must be MOST vulnerable. [NII] and Hα must be LEAST vulnerable. CIV and MgII in between.
4. Can you find a single atomic/plasma parameter that, when plotted against P_empirical, gives a monotonic relationship?
5. If the answer involves multiple parameters, what combination works? 

CONTEXT: This appears in a study of whether the intergalactic medium selectively degrades spectral information. The degradation is NOT dust, NOT noise, NOT selection — it survives all controls. It depends on foreground structure (z-matched test, p=10⁻⁵). It follows a sigmoid with threshold at z≈0.82.

Think step by step. Consider every line's formation physics individually. What makes Lyα special? What makes [NII] immune?
```

---

## 🔍 GROK — Recent Observational Evidence

```
Search for recent papers and observations (2024-2026) in these specific categories. I need actual citations, not summaries of general knowledge:

1. ANOMALOUS LINE RATIOS AT HIGH-z:
Any papers reporting that emission line ratios (BPT diagram, [OIII]/Hβ, [NII]/Hα, CIV/CIII, etc.) behave unexpectedly at z > 0.8 compared to local calibrations. Especially anything about line ratio "evolution" that doesn't match photoionization models.

2. LENSED QUASAR SPECTRAL ANOMALIES:
Papers on multiply-imaged quasars where different images show different equivalent widths, line profiles, or flux ratios that CAN'T be explained by microlensing. Especially: HE1104-1805, RXJ1131-1231, HE0435-1223, SDSS1004+4112, Q0957+561.

3. DESI YEAR 1 DARK ENERGY:
The specific w₀-wa results from DESI BAO. What is their reported w₀ and wa? At what redshift does their dark energy equation of state deviate most from w=-1? Does the deviation have a sigmoid-like shape?

4. JWST "IMPOSSIBLE" GALAXIES:
Papers about high-z galaxies (z>6) whose emission line properties (line widths, equivalent widths, ratios) don't match expectations. Any claims that spectroscopic measurements disagree with photometric estimates.

5. HUBBLE TENSION — SPECTROSCOPIC ANGLE:
Any 2024-2026 papers suggesting the Hubble tension could involve a systematic in spectroscopic distance indicators (Cepheid metallicity effects, SNe Ia spectral evolution, TRGB calibration issues).

6. COSMIC WEB EFFECT ON SPECTRA:
Any papers studying how intervening large-scale structure affects background source spectra beyond simple absorption. Especially anything about foreground filaments correlating with background quasar properties.

7. FRB DM ANOMALIES:
Papers reporting FRBs whose DM doesn't match expectations from the Macquart relation, especially at DM > 500 pc/cm³. Any studies of how FRB spectral properties change with DM.

For each paper found: give the title, authors, year, journal, and the specific result that's relevant. I need real citations I can look up.
```

---

## HOW TO USE THESE

1. **Gemini**: Paste the simulation prompt. If it can run code, great. If not, ask it to write Python and run locally.
2. **GPT**: Paste the susceptibility prompt. Use GPT-4/o1 in a fresh conversation for maximum reasoning.
3. **Grok**: Paste the search prompt. Grok has real-time X/web access — it should find recent preprints.

**When you get results back, send them to me and I'll integrate.**
