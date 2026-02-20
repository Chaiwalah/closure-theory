# Closure Theory Test Suite — Setup

## Quick Start
```bash
pip install -r requirements.txt
python closure_test.py
```

## What It Does
Downloads Pantheon+SH0ES SN catalog and Planck y-map automatically, then runs 5 tests:

1. **Triple Coherence** — Do Δμ, color residual, and χ²/dof ALL correlate with LOS baryonic density?
2. **Survey-Split** — Does the signal replicate across independent surveys (DES, SNLS, PS1, etc)?
3. **Angular Scale** — At which aperture (5', 10', 30', 60') is the signal strongest? Physical vs artifact?
4. **High-z Coupling** — Does color–distance correlation strengthen with redshift? (accumulation test)
5. **Nested Models** — Does adding a closure drift term improve AIC/BIC over baseline?

## Output
All results go to `results/`:
- `closure_test_results.json` — all numbers
- `test1_triple_coherence.png` — main result plot
- `test2_survey_split.png` — replication check
- `test3_angular_scale.png` — artifact diagnosis
- `test4_highz_coupling.png` — accumulation test
- `pantheon_with_y.csv` — processed data with y-values

## Data Sources
- Pantheon+SH0ES: https://github.com/PantheonPlusSH0ES/DataRelease
- Planck y-map (MILCA): ~100MB FITS file, auto-downloaded

## What to Look For
- **Triple coherence**: All three observables correlating with y in the same direction = channel effect
- **Survey replication**: Signal in multiple surveys = not calibration
- **Angular scale**: Peak at physical scale (30'+) = real; peak at beam size (10') = maybe artifact
- **z-accumulation**: r(c, Δμ) increasing with z = path-length dependent = channel
- **BIC improvement**: ΔBIC < -6 for closure model = strong evidence
