# Closure Theory — Humza Hafeez
*Started: 2026-02-20*

## Core Idea
Photons don't just lose **energy** over cosmological distances — they lose **information**. The EM propagation channel (baryonic medium: free electrons, plasma, IGM) degrades signal fidelity in a way that accumulates with distance and baryonic column density.

This is a **signal processing** framing of cosmology, not a geometry/energy one.

## Key Predictions
1. SN light curve fit quality (χ²/dof) should correlate with line-of-sight baryonic density (Planck tSZ y-map)
2. Effect strengthens with redshift (longer path = more accumulation)
3. ΛCDM predicts NO such correlation — signal quality should be independent of foreground baryons after extinction correction
4. High-z SNe should show stronger color–distance coupling than low-z (Wolf-test style anomaly)
5. FRB dispersion measure (DM) as baryon-channel ruler could reveal channel-specific drift

## The Make-or-Break Test
**Pantheon+ χ²/dof vs Planck y-map cross-match**
- Download Pantheon+ SN catalog (public, GitHub)
- Download Planck y-map (public, HEALPix)
- For each SN: extract tSZ y-value at its sky position
- Plot χ²/dof vs y, binned by redshift
- Positive correlation = anomaly paper. No correlation = rethink.

## Theoretical Framework
- Closure drift parameter: Γ₁ (EM channel), predicted ~10⁻³ Mpc⁻¹
- Deconfounded residual: Δμ⊥
- Nested model test: Baseline ΛCDM vs Baseline + Γ·(baryon LOS metric) vs Baseline + Γz·(baryon × z interaction)
- Compare via AIC/BIC

## Information-First Framing
Instead of "did the photon get dimmer?" ask:
- Did the light curve lose high-frequency temporal content?
- Did the spectrum lose fine structure?
- Did fit residuals become more template-like?
These are **channel/information** signatures, not cosmological ones.

## Why This Might Be Hidden
- Physicists think energy/geometry, not signal/channel
- Nobody cross-matches SN fit quality with baryonic foreground maps
- FRB DM catalogs and SN catalogs have never been systematically cross-correlated
- "Template-fit smoothing behind baryons" isn't a category anyone looks for

## Connections to Known Tensions
| Tension | Closure Interpretation |
|---------|----------------------|
| Hubble tension (5σ) | Local vs CMB chains traverse different baryonic paths → channel-specific disagreement |
| Lithium problem (3x deficit) | Li-7 measured via narrow absorption lines → fine spectral structure degrades over cosmological paths → systematic undercount |
| S8/σ₈ tension (~10% low) | Lensing measured via galaxy shape distortion → high-spatial-frequency image content degrades → underestimate shear |

## Data Sources
- **Pantheon+SH0ES**: GitHub (public, .dat file)
- **Planck y-map**: ESA archive (public, HEALPix FITS)
- **CHIME/FRB catalog**: ~2,500 localized FRBs with DM (public CSV)

## Priority Stack
1. ⬜ Run χ²/dof vs y-map test (one-afternoon calculation)
2. ⬜ Document y_10 result robustness (χ² improvement with baryonic conditioning)
3. ⬜ High-z color–distance coupling analysis (Wolf-test anomaly)
4. ⬜ Build full-sky DM_excess field from FRBs (Gaussian process interpolation)
5. ⬜ Nested model comparison (AIC/BIC) on Pantheon+
6. ⬜ Write anomaly paper (NOT theory paper first)

## IP Protection
- [ ] Private GitHub repo with timestamped commits
- [ ] Email PDF summary to self
- [ ] All agent transcripts saved with dates
