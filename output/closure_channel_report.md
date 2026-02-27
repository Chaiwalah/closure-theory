# Closure Channel Test — First-Pass Results

**Date:** February 19, 2026
**Author:** Humza Hafeez
**Analysis pipeline:** VPS (NSIDE=512) + PC (NSIDE=1024, RTX 5090)

---

## 1. Executive Summary

We performed the first empirical test of the Closure Theory channel-drift hypothesis using 1,454 Type Ia supernovae from the Pantheon+SH0ES compilation cross-matched with Planck PR3 CMB lensing (κ) and PR2 thermal Sunyaev–Zel'dovich (tSZ Compton-y) all-sky maps.

**Primary finding:** No statistically significant closure-channel signal was detected in compensated aperture measurements at any scale or redshift bin. One marginal detection (tSZ y at 30 arcmin, p = 0.021) was identified but fails the redshift-gradient test required by the theory, indicating a local-environment systematic rather than a propagation effect.

**Status:** The closure-channel hypothesis is not falsified by this analysis — the test lacks the statistical power and tracer specificity required for a definitive conclusion — but no supporting evidence was found.

---

## 2. Data

### 2.1 Supernova Sample

| Property | Value |
|----------|-------|
| Source catalog | Pantheon+SH0ES (Scolnic et al. 2022) |
| Raw light curves | 1,701 |
| Unique SNe Ia | 1,543 |
| Clean sample (after cuts) | 1,454 |
| Redshift range | 0.010 – 2.261 |
| Cuts applied | z > 0.01, E(B−V)_MW < 0.15, valid κ, valid Δμ |

Duplicate observations of the same SN (from different surveys) were averaged. Distance moduli (MU_SH0ES) and associated diagonal errors were taken directly from the Pantheon+ release.

### 2.2 Maps

| Map | Source | Resolution | NSIDE (working) |
|-----|--------|-----------|-----------------|
| CMB lensing convergence κ | Planck PR3 COM_Lensing_4096_R3.00 MV (alm→map) | ~2 arcmin native | 1024 |
| Compton-y (tSZ) | Planck PR2 COM_CompMap_YSZ_R2.02 MILCA | 10 arcmin FWHM | 1024 |

### 2.3 ΛCDM Reference Cosmology

Flat ΛCDM with H₀ = 73.04 km/s/Mpc, Ωm = 0.334 (Pantheon+SH0ES best fit). Luminosity distances computed via numerical integration.

---

## 3. Method

### 3.1 Map Extraction

For each SN position (converted from ICRS to Galactic coordinates):

- **Disc averages** at R = 5, 10, 20, 30 arcmin
- **Compensated apertures**: inner disc minus surrounding annulus at (5,10), (10,20), (20,40) arcmin pairs
- Mask validity flags and pixel counts recorded per SN

Compensated apertures are the primary observable. They suppress large-scale gradients and isolate local structure around each SN, which is the physically relevant quantity for a propagation-dependent effect.

### 3.2 Hubble Residuals

**Raw residual:**

Δμᵢ = μ_obs,i − μ_ΛCDM(zᵢ)

**Deconfounded residual (Δμ⊥):**

Regressed out via OLS:
- Milky Way dust extinction (E(B−V)_MW)
- Host galaxy mass step (split at log M = 10)
- Survey indicator variables (top 7 surveys as dummies)

The deconfounded residual removes known systematics before testing for structure-dependent correlations. All results below use Δμ⊥.

---

## 4. Results

### 4.1 Full-Sample Correlations

| Tracer | Pearson r | p-value | N |
|--------|-----------|---------|---|
| κ_comp (10–20 arcmin) | −0.024 | 0.362 | 1,454 |
| y_comp (10–20 arcmin) | +0.009 | 0.728 | 1,454 |
| κ₃₀ (raw disc) | −0.022 | 0.393 | 1,454 |
| **y₃₀ (raw disc)** | **+0.060** | **0.021** | **1,454** |

The compensated apertures — which are the theoretically preferred estimator — show no correlation. The one detection (y₃₀, p = 0.021) occurs only in the raw disc average, which is dominated by the local environment rather than the integrated line-of-sight.

### 4.2 Redshift Splits

**κ_comp (10–20 arcmin) vs Δμ⊥:**

| Redshift bin | N | r | p |
|-------------|---|---|---|
| z < 0.1 | 495 | +0.012 | 0.788 |
| 0.1 < z < 0.5 | 749 | −0.053 | 0.145 |
| z > 0.5 | 210 | −0.023 | 0.741 |

No redshift trend. The mid-z hint (p = 0.145) has the wrong sign (negative) for closure drift.

**y_comp (10–20 arcmin) vs Δμ⊥:**

| Redshift bin | N | r | p |
|-------------|---|---|---|
| z < 0.1 | 495 | +0.024 | 0.591 |
| 0.1 < z < 0.5 | 749 | −0.034 | 0.352 |
| z > 0.5 | 210 | −0.066 | 0.340 |

No redshift trend. If anything, the correlation is weakly *negative* at high-z — opposite to the closure prediction.

**y₃₀ (raw disc) vs Δμ⊥ — the one detection:**

| Redshift bin | N | r | p |
|-------------|---|---|---|
| **z < 0.1** | **495** | **+0.111** | **0.014** |
| 0.1 < z < 0.5 | 749 | −0.022 | 0.556 |
| z > 0.5 | 210 | +0.049 | 0.481 |

The signal is concentrated entirely in the low-z bin. This is the **opposite** of the closure prediction (signal should strengthen with redshift as channel drift accumulates). The low-z concentration is characteristic of local-environment contamination — peculiar velocities, host-cluster membership, and local density effects that dominate at z < 0.1.

### 4.3 Scale Dependence

The slope dΔμ⊥/dκ was computed at each aperture radius. For κ, the slope is weakly negative at all scales (5–30 arcmin) with no preferred physical scale. For y, a positive slope appears only at the largest aperture (30 arcmin), consistent with the local-environment explanation (larger aperture captures more cluster/filament signal).

No scale-dependent saturation or turnover was observed — a feature that would be expected if a physical propagation mechanism were at work.

---

## 5. Interpretation

### 5.1 What the Data Says

The Planck CMB lensing convergence shows no correlation with Hubble residuals at any scale, aperture configuration, or redshift bin. This is the cleanest projected-mass tracer available and the most natural proxy for total line-of-sight structure. Its null result places a constraint on any theory that predicts EM distance modulation by intervening mass.

The Planck tSZ Compton-y parameter shows a marginal correlation at 30 arcmin (p = 0.021), but this signal:
1. Vanishes in compensated apertures
2. Is concentrated at z < 0.1
3. Does not strengthen with redshift

All three properties point to a local-environment systematic, not a propagation effect.

### 5.2 What This Does NOT Rule Out

This analysis has important limitations:

1. **Wrong primary observable.** The closure-channel test (§3 of the theory paper) predicts a difference between EM and GW distances from the *same source*. We tested EM distances alone against structure tracers — a second-order proxy, not the primary signal.

2. **κ traces total mass, not baryonic medium.** Closure drift is hypothesized to depend on the baryonic channel (free electrons, plasma), not dark matter. CMB lensing κ is ~85% dark matter dominated. The tSZ y-parameter is a better baryonic tracer, but the MILCA y-map has significant foreground contamination at the scales we probe.

3. **No matched line-of-sight galaxy overdensity.** The most closure-specific tracer — galaxy counts weighted by proximity to the SN redshift — was not computed in this run due to API rate limits on the DESI Legacy Survey. This tracer isolates the baryonic medium *along the actual photon path*, rather than the full projected column.

4. **Statistical power.** With 1,454 SNe and Planck resolution, a correlation of |r| < 0.05 is indistinguishable from noise. A Γ₁ coupling of order 10⁻⁴ Mpc⁻¹ would produce correlations at the |r| ~ 0.01–0.02 level — well below our detection threshold.

5. **No GW distance data.** Only one standard siren with EM counterpart exists (GW170817, z = 0.01). The primary closure test requires O4/O5 LIGO runs (2025–2028) and ultimately Einstein Telescope / Cosmic Explorer (2030s).

### 5.3 Falsification Status

Per the falsification criteria defined in the theory paper (§3.9):

- **Criterion 1** (𝔼[Δ] = 0 and Cov(Δ,X) = 0): Partially met for κ and y_comp. The test cannot be applied to the primary observable (EM vs GW distance ratio) with current data.
- **Criteria 2–5**: Not applicable to this proxy analysis.

**Verdict:** The closure-channel hypothesis is **not falsified** by this analysis, but receives **no empirical support** from Planck×Pantheon+ cross-correlation. The theory remains viable pending:
- Standard siren catalogs from O4/O5
- Matched LOS galaxy overdensity analysis
- Higher-resolution baryonic tracers (FRB dispersion measures, kinetic SZ)

---

## 6. Next Steps

| Priority | Action | Data source | Timeline |
|----------|--------|-------------|----------|
| 1 | Matched LOS galaxy overdensity (δg) | DESI Legacy Survey DR10 photo-z | Days (API queries) |
| 2 | FRB dispersion measure cross-match | CHIME/FRB catalog | Available now |
| 3 | Kinetic SZ velocity reconstruction | ACT/SPT + spectroscopic surveys | Public data |
| 4 | Standard siren EM-GW comparison | LIGO O4 events | As events are published |
| 5 | Full covariance matrix analysis | Pantheon+ STAT+SYS cov (32MB, on disk) | Requires MCMC framework |

The matched LOS δg (Priority 1) is the most promising near-term test. It is the only tracer that isolates baryonic structure specifically along the photon propagation path, weighted by proximity to the SN redshift. No published analysis has performed this exact cross-correlation with Pantheon+ data.

---

## 7. Reproducibility

All code, data, and maps are archived at:

| Asset | Location |
|-------|----------|
| Pantheon+ raw data | `/root/closure-theory/data/pantheon_plus.dat` |
| Planck lensing alm | `/root/closure-theory/maps/COM_Lensing_4096_R3.00/MV/dat_klm.fits` |
| Planck MILCA y-map | `/root/closure-theory/maps/COM_CompMap_YSZ_R2.02/milca_ymaps.fits` |
| Analysis script | `/root/closure-theory/analyze.py` |
| Output table | `/root/closure-theory/output/sn_closure_full.csv` |
| Diagnostic plots | `/root/closure-theory/plots/` |

Machine: WSL2 on Windows 11, RTX 5090, 64GB RAM. Python 3.12, healpy 1.19.0, astropy 7.2.0.

---

## Appendix: Deconfounding Coefficients

| Parameter | Coefficient |
|-----------|-------------|
| E(B−V)_MW | (fitted, removes dust correlation) |
| Host mass step (log M > 10) | (fitted, removes ~0.04 mag step) |
| Survey dummies (7) | (fitted, removes inter-survey calibration offsets) |

The deconfounded residual Δμ⊥ has mean ≈ 0 and no residual correlation with any control variable by construction.
