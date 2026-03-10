# Closure Theory — Geometric Non-Stationarity of Observable Correlations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18928542.svg)](https://doi.org/10.5281/zenodo.18928542)

**Humza Hafeez** · Independent Researcher · [ORCID 0009-0000-0853-8807](https://orcid.org/0009-0000-0853-8807)

---

## Paper 1 — Geometric Non-Stationarity of the Type Ia Supernova Standardization Manifold

**Status**: Submitted to Physical Review Letters (March 2026)  
**Preprint**: [Zenodo DOI 10.5281/zenodo.18928542](https://doi.org/10.5281/zenodo.18928542)

### Abstract

We demonstrate that the Tripp standardization basis for Type Ia supernovae is geometrically nonstationary — the covariance structure of the (x₁, c, Δμ) space is strongly redshift-dependent (Box's M test, p = 2.3 × 10⁻⁷). Correcting this nonstationarity absorbs 99% of the reported wₐ dark energy evolution signal from Pantheon+.

### Key Results

| Result | Value |
|--------|-------|
| Residual–color coupling reversal | ρ = −0.785, p = 3 × 10⁻³⁰ |
| Interaction term improvement | Δχ² = 269 (1,590 SNe) |
| wₐ shift under correction | −1.88 → −0.02 |
| ΛCDM+nonstationary vs w₀wₐ+standard | Δχ² = 37 in favor |
| Host mass step via color channel | 71% mediated |

### Cross-Domain Confirmation

The same pattern — sigmoidal collapse of observable correlations with cosmological distance — appears independently in:

- **750,414 quasars** (SDSS DR16Q) — emission line equivalent width coupling collapses
- **535 fast radio bursts** (CHIME Catalog 1) — width-spectral index correlation vanishes past DM ≈ 500
- **130,000 galaxies** (DESI) — environmental density coupling follows eating law

## Data Sources

All analyses use publicly available data:

| Dataset | Source | Objects |
|---------|--------|---------|
| Pantheon+SH0ES | [GitHub](https://github.com/PantheonPlusSH0ES/DataRelease) | 1,590 SNe Ia |
| SDSS DR16Q | [SDSS](https://www.sdss.org/dr16/algorithms/qso_catalog/) | 750,414 quasars |
| CHIME/FRB Cat 1 | [CHIME](https://www.chime-frb.ca/catalog) | 535 FRBs |
| DESI Galaxy Catalog | [DESI](https://data.desi.lbl.gov/) | ~130,000 galaxies |

## Reproducing Results

```bash
# Install dependencies
pip install numpy scipy pandas matplotlib astropy scikit-learn

# Paper 1 — Core SN analysis
python paper1/scripts/closure_test.py

# Paper 1 — Quasar cross-domain test
python paper1/scripts/closure_test_quasar.py

# Paper 1 — FRB cross-domain test
python paper1/scripts/closure_test_frb_v3.py

# Paper 1 — Unified cross-domain analysis
python paper1/scripts/closure_cross_domain.py

# Paper 1 — N_modes mechanism
python paper1/scripts/closure_multipole_mechanism.py

# Paper 1 — Generate figures
python paper1/scripts/plot_paper_figures.py
```

## Repository Structure

```
paper1/
├── scripts/          # Analysis scripts that generate Paper 1 results
├── results/          # Output from each analysis
├── main.tex          # LaTeX source
├── references.bib    # Bibliography
└── *.png             # Figures

data/                 # Local data files (Pantheon+, CHIME FRBs)
zenodo/               # Preprint PDF
references/           # Supporting literature
exploratory/          # Development scripts, agent reports, exploratory analyses
requirements.txt      # Python dependencies
```

## License

Code: MIT · Paper: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Citation

```bibtex
@article{hafeez2026geometric,
  title={Geometric Non-Stationarity of the Type Ia Supernova 
         Standardization Manifold and Its Impact on Dark Energy Inference},
  author={Hafeez, Humza},
  year={2026},
  doi={10.5281/zenodo.18928542},
  publisher={Zenodo}
}
```
