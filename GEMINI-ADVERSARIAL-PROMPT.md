# GEMINI — Adversarial Audit + Gap Analysis

```
I need you to play hostile referee AND discovery engine simultaneously. I have an empirical result with a proposed mechanism, and I need you to:

A) Find every hole in my case
B) Tell me what I'm not accounting for that could explain the data WITHOUT my mechanism
C) If my mechanism survives, tell me what one test would make it undeniable

THE EMPIRICAL RESULTS (these are raw measurements, not interpretations):

1. Across 750,414 SDSS DR16Q quasars, the Spearman correlation between emission line EW and FWHM decays monotonically with redshift for CIV (r drops from 0.248 to 0.135 across z=1.3-3.5, trend p=10⁻⁹) and Lyα (r drops from 0.268 to 0.050 across z=1.9-3.7, p=4×10⁻⁴).

2. Meanwhile, [OIII] 5007 (z=0.3-1.1) and [NII] 6583 (z=0.2-0.6) show NO decorrelation or even positive trends. [SII] 6718 shows decorrelation.

3. EW kurtosis rises with z for CIV (r=+0.979, p=10⁻¹⁰), CIII (r=+0.975), MgII (r=+0.982). FWHM kurtosis does NOT rise for these same lines.

4. Z-matched quasar pairs (same z, different foreground density estimated from projected galaxy counts): EW differs (p=5.67×10⁻⁵ for [OIII]), FWHM is stable (p=0.49).

5. 40/40 local EW predictions fail when applied globally (KS p=0 for all pairs).

6. The decorrelation rate dr/dz scales with our empirically-assigned "susceptibility" P at r=-0.726, p=0.017.

7. Three source classes (1,590 SNe Ia, 750,414 quasars, 535 FRBs) all show sigmoid-like transitions at consistent DM thresholds.

THE PROPOSED MECHANISM:
Refractive decoherence in the cosmic web. Plasma lenses differentially magnify/demagnify unresolved source subregions. Lines whose emissivity varies across subregions (high ∂ln j/∂ln n_e or ∂ln j/∂ln T) get their observed EW scrambled. Lines locked by atomic branching ratios (same upper level) are immune.

Susceptibility formula: P = ‖∇ ln j‖ × E_trap × (1 - B_lock)

Threshold at z≈0.8 attributed to cosmic web fragmentation from dark energy transition.

WHAT I NEED FROM YOU:

1. CONVENTIONAL EXPLANATIONS I HAVEN'T KILLED:
   - Could the Baldwin effect (EW-luminosity anticorrelation) produce ALL of these signatures? Including the z-matched foreground test? Including the kurtosis asymmetry?
   - Could spectral fitting pipeline systematics (template mismatch at high z, lower SNR) produce the monotonic CIV decorrelation?
   - Could quasar evolution (BLR structure changes with cosmic time, different accretion modes at high z) explain why EW-FWHM decouples? This is NOT the same as "source evolution" — this is about whether the INTRINSIC EW-FWHM relationship itself evolves.
   - Could selection effects (Malmquist bias, flux limits that select different quasar populations at different z) create an apparent decorrelation?
   
2. SPECIFIC WEAKNESSES IN MY MECHANISM:
   - The P correlation (r=0.721) used trapping leverage values (E_trap) that I assigned semi-subjectively. If you recompute P using ONLY published atomic data (Einstein A coefficients, critical densities, oscillator strengths), does the correlation hold?
   - The z-matched test used projected galaxy counts as a foreground density proxy. How reliable is this? Could foreground galaxy density correlate with the BACKGROUND quasar properties through some selection channel?
   - My "threshold matches dark energy transition" claim: is this a coincidence? How many other astrophysical transitions happen at z≈0.6-1.0?

3. THE THING I'M PROBABLY MISSING:
   - Is there a known effect in AGN physics that predicts EW-FWHM decorrelation with z that I should be testing against?
   - Is there published work on quasar emission line evolution with z that already explains my observations without invoking IGM effects?
   - What does the reverberation mapping literature say about EW-FWHM relationships and how they change with luminosity/redshift?

4. IF MY MECHANISM IS RIGHT, WHAT'S THE KILLSHOT TEST?
   - What single observation would make this undeniable?
   - Not "get more data" — what SPECIFIC measurement with EXISTING data would definitively separate IGM refractive decoherence from all conventional explanations?

Be brutal. I'd rather know I'm wrong now than after submitting to Nature.
```
