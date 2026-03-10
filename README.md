# Geometric Non-Stationarity of SN Ia Standardization

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18928542.svg)](https://doi.org/10.5281/zenodo.18928542)

**Humza Hafeez** · Independent Researcher · [ORCID 0009-0000-0853-8807](https://orcid.org/0009-0000-0853-8807)

## Summary

We demonstrate that the Tripp standardization basis for Type Ia supernovae is geometrically nonstationary — the covariance structure of the (x₁, c, Δμ) space is strongly redshift-dependent (Box's M test, p = 2.3 × 10⁻⁷). Correcting this nonstationarity absorbs 99% of the reported wₐ dark energy evolution signal from Pantheon+.

**Key results:**
- Residual–color coupling reverses sign across the Pantheon+ redshift range (ρ = −0.785, p = 3 × 10⁻³⁰)
- Adding interaction terms (z × x₁, z × c) to the Tripp equation yields Δχ² = 269 for 1,590 SNe
- Under corrected standardization, wₐ shifts from −1.88 to −0.02
- ΛCDM + nonstationary Tripp outperforms w₀wₐ + standard Tripp by Δχ² = 37
- 71% of the host-galaxy mass step is mediated through the color channel

## Paper

- **Preprint**: [Zenodo DOI 10.5281/zenodo.18928542](https://doi.org/10.5281/zenodo.18928542)
- **Submitted to**: Physical Review Letters (March 2026)
- **LaTeX source**: [`paper1/`](paper1/)

## Data

All analyses use publicly available data:
- [Pantheon+SH0ES](https://github.com/PantheonPlusSH0ES/DataRelease) (1,590 SNe Ia)
- [SDSS DR16Q](https://www.sdss.org/dr16/algorithms/qso_catalog/) (750,414 quasars)
- [CHIME/FRB Catalog 1](https://www.chime-frb.ca/catalog) (535 FRBs)

## Reproducing Results

```bash
# Install dependencies
pip install numpy scipy pandas matplotlib astropy scikit-learn

# Core analysis (Paper 1)
python closure_paper1_analysis.py

# Cross-domain tests
python closure_test_quasar.py
python closure_test_frb_v3.py
python closure_void_galaxy_test.py
```

## Repository Structure

```
paper1/              # LaTeX source for Paper 1
results_*/           # Output directories for each analysis
closure_*.py         # Analysis scripts
zenodo/              # Preprint PDF
```

## License

Code: MIT · Paper: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Citation

```bibtex
@article{hafeez2026geometric,
  title={Geometric Non-Stationarity of the Type Ia Supernova Standardization Manifold and Its Impact on Dark Energy Inference},
  author={Hafeez, Humza},
  year={2026},
  doi={10.5281/zenodo.18928542},
  publisher={Zenodo}
}
```
