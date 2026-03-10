# Paper 2: The Diagnostic Sensitivity Law
## Working Title
"A Diagnostic Sensitivity Law for Spectral Correlation Loss Across Cosmological Distances"

## Target: Nature Astronomy (Letter, ~3,500 words) or ApJ Letters

## Core Result
Observable-pair correlations degrade at rates proportional to diagnostic sensitivity:
- Doublet ladder: r = -0.975 (p = 0.005) across emission lines
- Quasar EW coupling: r = 0.82 → 0.03 (collapse with z)
- FRB width-spectral index: r = 0.27 → 0.00 past DM ≈ 500
- Eating law: eaten = 0.048(1-MI₀) + 0.012 across 130K DESI galaxies
- Three source classes, one pattern, no conventional mechanism

## Key Scripts (existing)
- closure_test_quasar.py
- closure_test_frb_v3.py
- closure_void_galaxy_test.py
- closure_sfr_split_test.py
- closure_crossover_map.py
- closure_eating_threshold.py
- closure_bridge_tests.py
- closure_doublet_ladder.py

## Datasets
- [x] SDSS DR16Q (750,414 quasars)
- [x] CHIME FRB Catalog 1 (535 + 186 localized)
- [x] DESI Year 1 (130K emission-line galaxies)
- [x] Pantheon+SH0ES (1,590 SNe Ia)
- [ ] LAMOST DR9? (millions of stellar spectra)
- [ ] GALAH DR3? (galactic archaeology survey)
- [ ] 4MOST? (upcoming)

## Structure
1. Introduction — spectral information assumed lossless, anomalies suggest otherwise
2. Diagnostic Sensitivity Classification — q ∈ [0,1] from atomic physics
3. Quasar Analysis — EW coupling, FWHM, Baldwin control
4. FRB Analysis — width-spectral index, DM-RM decoupling
5. Galaxy Analysis — eating law, void test, SFR split
6. Cross-Domain Consistency — same pattern, different physics
7. Discussion — no conventional mechanism, implications

## Conventional Alternatives to Kill
- [ ] Selection effects (Malmquist, aperture)
- [ ] Evolution (intrinsic source evolution with z)
- [ ] Dust/extinction (wavelength-dependent)
- [ ] IGM absorption (Lyman-alpha forest effects)
- [ ] Instrumental systematics (survey-specific)
- [ ] Statistical artifacts (look-elsewhere effect)

## Figures (planned)
1. Doublet ladder — degradation vs diagnostic sensitivity
2. Quasar EW coupling collapse with z
3. FRB width-spectral index vanishing
4. Eating law — linear relationship
5. Cross-domain sigmoid comparison
6. Null controls (branching ratios, baryonic proxies)

## Status: OUTLINE STAGE
