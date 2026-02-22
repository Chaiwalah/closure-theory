#!/usr/bin/env python3
"""
β(z) × BAO Tension Overlay
============================
Plot the coincidence between:
1. SN color standardization breakdown (β(z) from our Hubble v2 analysis)
2. DESI BAO residuals vs ΛCDM (from NotebookLM extraction)
3. Environment-dependent correlation degradation (from void tests)

No new data — just synthesizing what we already have into one view.
"""
import numpy as np
import json, os

RESULTS_DIR = 'results_bao_overlay'
os.makedirs(RESULTS_DIR, exist_ok=True)

# ============================================================
# DATA: All from our previous analyses + DESI BAO paper
# ============================================================

# β(z) from closure_hubble_v2.py results
# β drops from ~2.42 at low z to ~1.85 at high z, breakpoint z≈0.6
beta_z = {
    'description': 'SN Ia color standardization coefficient vs redshift',
    'source': 'closure_hubble_v2.py on Pantheon+ sample',
    'bins': [
        {'z_mid': 0.025, 'beta': 2.42, 'beta_err': 0.15},
        {'z_mid': 0.075, 'beta': 2.38, 'beta_err': 0.12},
        {'z_mid': 0.15, 'beta': 2.35, 'beta_err': 0.10},
        {'z_mid': 0.30, 'beta': 2.28, 'beta_err': 0.11},
        {'z_mid': 0.50, 'beta': 2.10, 'beta_err': 0.14},
        {'z_mid': 0.70, 'beta': 1.92, 'beta_err': 0.16},
        {'z_mid': 0.90, 'beta': 1.85, 'beta_err': 0.20},
    ],
    'breakpoint_z': 0.60,
    'delta_chi2': -51.76,
    'note': 'Δχ²=-51.76 for 2 extra params vs constant β'
}

# DESI BAO residuals (from NotebookLM extraction of DESI 2024 VI)
bao_residuals = {
    'description': 'DESI DR1 BAO distance residuals vs Planck ΛCDM',
    'source': 'DESI 2024 VI (arXiv:2404.03002)',
    'bins': [
        {'z_eff': 0.30, 'tracer': 'BGS', 'DV_resid': -0.14, 'DV_err': 0.15,
         'DM_resid': None, 'DH_resid': None, 'type': '1D'},
        {'z_eff': 0.51, 'tracer': 'LRG1', 
         'DM_resid': +0.13, 'DM_err': 0.25, 
         'DH_resid': -1.78, 'DH_err': 0.61, 'type': '2D',
         'note': 'LARGEST anomaly: D_H ~2.9σ below ΛCDM'},
        {'z_eff': 0.71, 'tracer': 'LRG2',
         'DM_resid': -0.84, 'DM_err': 0.32,
         'DH_resid': -0.09, 'DH_err': 0.60, 'type': '2D',
         'note': 'D_M ~2.6σ below ΛCDM'},
        {'z_eff': 0.93, 'tracer': 'LRG3+ELG1',
         'DM_resid': -0.21, 'DM_err': 0.28,
         'DH_resid': +0.26, 'DH_err': 0.35, 'type': '2D',
         'note': 'Consistent with ΛCDM'},
        {'z_eff': 1.32, 'tracer': 'ELG2',
         'DM_resid': -0.22, 'DM_err': 0.69,
         'DH_resid': -0.29, 'DH_err': 0.42, 'type': '2D',
         'note': 'Consistent with ΛCDM'},
        {'z_eff': 1.49, 'tracer': 'QSO', 'DV_resid': +0.04, 'DV_err': 0.67,
         'DM_resid': None, 'DH_resid': None, 'type': '1D',
         'note': 'Consistent with ΛCDM'},
        {'z_eff': 2.33, 'tracer': 'Ly-α', 
         'DM_resid': 0, 'DH_resid': 0, 'type': 'approx',
         'note': 'α_∥=0.993, consistent with ΛCDM'},
    ]
}

# Environment effect from our void tests
env_effect = {
    'description': 'CIV/MgII Spearman Δρ (dense-sparse) after full photometry matching',
    'source': 'closure_photometry_control.py',
    'bins': [
        {'z_mid': 1.45, 'delta_rho': +0.0123},
        {'z_mid': 1.75, 'delta_rho': +0.0160},
        {'z_mid': 2.05, 'delta_rho': +0.0073},
        {'z_mid': 2.40, 'delta_rho': +0.0258},
    ]
}

# Closure breakpoints from various tests
breakpoints = {
    'description': 'Measured closure breakpoints from different probes',
    'points': [
        {'probe': '[SII] doublet', 'z0': 0.35, 'source': 'DESI emfit'},
        {'probe': 'β(z) SN color', 'z0': 0.60, 'source': 'Pantheon+'},
        {'probe': 'MgII doublet', 'z0': 0.82, 'source': 'SDSS DR16Q'},
        {'probe': 'Hβ-MgII cross', 'z0': 1.05, 'source': 'SDSS DR16Q'},
        {'probe': 'CIII-CIV cross', 'z0': 1.21, 'source': 'SDSS DR16Q'},
    ]
}

# ============================================================
# ANALYSIS: Quantify the BAO-β(z) coincidence
# ============================================================
print("=" * 70, flush=True)
print("β(z) × BAO TENSION OVERLAY ANALYSIS", flush=True)
print("=" * 70, flush=True)

# BAO tension magnitude (normalized)
print("\n[1] BAO tension by redshift:", flush=True)
for b in bao_residuals['bins']:
    if b['type'] == '2D' and b.get('DH_resid') is not None:
        # Combined tension: sqrt(DM² + DH²) / sqrt(DM_err² + DH_err²)
        tension = np.sqrt((b.get('DM_resid',0))**2 + (b.get('DH_resid',0))**2)
        err = np.sqrt((b.get('DM_err',1))**2 + (b.get('DH_err',1))**2)
        sigma = tension / err if err > 0 else 0
        print(f"  z={b['z_eff']:.2f} ({b['tracer']}): |residual|={tension:.2f}, ~{sigma:.1f}σ", flush=True)
    elif b['type'] == '1D':
        if b.get('DV_resid') is not None:
            sigma = abs(b['DV_resid']) / b.get('DV_err', 1)
            print(f"  z={b['z_eff']:.2f} ({b['tracer']}): |DV resid|={abs(b['DV_resid']):.2f}, ~{sigma:.1f}σ", flush=True)

# β(z) values in the BAO tension window
print("\n[2] β(z) in the BAO tension window (z=0.4-0.8):", flush=True)
for b in beta_z['bins']:
    if 0.3 <= b['z_mid'] <= 0.9:
        print(f"  z={b['z_mid']:.2f}: β={b['beta']:.2f} ± {b['beta_err']:.2f}", flush=True)
print(f"  β drops from 2.28 (z=0.3) to 1.92 (z=0.7) — a 16% decline", flush=True)
print(f"  Breakpoint at z≈0.6, Δχ²=-51.76", flush=True)

# Key coincidence
print("\n[3] THE COINCIDENCE:", flush=True)
print(f"  BAO max tension:  z = 0.51 (D_H, ~2.9σ) and z = 0.71 (D_M, ~2.6σ)", flush=True)
print(f"  β(z) breakpoint:  z ≈ 0.60", flush=True)
print(f"  [SII] breakpoint: z ≈ 0.35", flush=True)
print(f"  BAO clean above:  z > 0.8", flush=True)
print(f"  β(z) flattens at: z > 0.7", flush=True)

# Window overlap
print(f"\n  Overlap window: z = 0.35 — 0.80", flush=True)
print(f"  Both probes show maximum anomaly in this window.", flush=True)
print(f"  Both probes are clean/consistent outside it.", flush=True)

# Statistical framing
print(f"\n[4] Independence argument:", flush=True)
print(f"  β(z): Measured from Type Ia SN color-luminosity relation (Pantheon+)", flush=True)
print(f"  BAO:  Measured from galaxy clustering acoustic scale (DESI LRG/ELG)", flush=True)
print(f"  Different objects, different physics, different pipelines.", flush=True)
print(f"  Shared: photons traversing the same z=0.4-0.8 volume of spacetime.", flush=True)

# Directional scrambling argument
print(f"\n[5] Directional inconsistency in BAO:", flush=True)
print(f"  z=0.51: D_H too SMALL (line-of-sight compressed)", flush=True)
print(f"  z=0.71: D_M too SMALL (transverse compressed)", flush=True)
print(f"  → Not a coherent systematic bias (would push both same direction)", flush=True)
print(f"  → Consistent with 'scrambling': distances unreliable in both directions", flush=True)
print(f"  → Different geometric projections disagree at neighboring redshifts", flush=True)

# Environment connection
print(f"\n[6] Environment connection:", flush=True)
print(f"  Our void test shows spectral correlations degrade FASTER in sparse environments", flush=True)
print(f"  BAO is measured from galaxy clustering — which is SPARSE in voids, DENSE in walls", flush=True)
print(f"  If spectral info degrades faster in voids, BAO measured THROUGH voids is less reliable", flush=True)
print(f"  BAO tension could be partially driven by void-traversing sightlines", flush=True)

# ============================================================
# ASSESSMENT: Can we make a serious claim?
# ============================================================
print(f"\n{'='*70}", flush=True)
print("CLAIM ASSESSMENT", flush=True)
print(f"{'='*70}", flush=True)

print(f"""
WHAT WE CAN CLAIM (empirically supported):
1. Spectral coupling (CIV/MgII) degrades with redshift in a way that survives
   all known systematic controls (SNR, luminosity, E(B-V), psf depth, nobs,
   CIV wind proxy, forward model, shuffle null). [200k quasars, DESI DR1]

2. This degradation is environment-dependent: denser 3D environments preserve
   correlations better than sparse ones, at fixed z, after full photometry
   matching. 4/4 bins, mean Δρ=+0.015. [120k quasars with photometry]

3. The effect scales with coupling looseness: CIV/NeV (weakest coupling) shows
   largest environment modulation (+0.022), CIV/MgII intermediate (+0.016),
   MgII/Hβ smallest (+0.006). [Pre-registered, git-timestamped prediction]

4. The redshift window where SN spectral standardization breaks down (z≈0.5-0.7,
   Δχ²=-51.76) coincides with where DESI BAO shows maximum tension with ΛCDM
   (z=0.51 at 2.9σ, z=0.71 at 2.6σ). Both probes are clean outside this window.

WHAT WE CANNOT YET CLAIM:
- That BAO tension is CAUSED by spectral degradation (correlation ≠ causation)
- That the environment effect is due to "information channel" physics (no mechanism proven)
- That dark energy is an artifact (we've shown β(z) affects ΩΛ by 9-30%, not that it's fake)
- That the Hubble tension is explained (suggestive, not demonstrated)

STRENGTH OF EVIDENCE:
- Environment effect: ROBUST (survives 8+ adversarial controls, pre-registered predictions)
- β(z)-BAO coincidence: SUGGESTIVE (two independent probes, same window, but no causal link)
- Coupling ordering: CONFIRMED for 3 pairs (need more, and at matched z-ranges)
- Galaxy emission lines: SUPPORTING only (z-range insufficient for breakpoints)

PUBLICATION READINESS:
- Paper 1 (environment-dependent spectral coupling): ~80% ready
  Need: independent pipeline cross-check (Wu & Shen), 1-2 more pair predictions at matched z
- Paper 2 (β(z) + cosmological implications): ~60% ready
  Need: red/blue ΩΛ convergence test, explicit BAO connection analysis
- Paper 3 (multi-channel information conservation): ~20% ready
  Need: gravity channel proxy at same z, total MI measurement
""", flush=True)

# Save everything
output = {
    'beta_z': beta_z,
    'bao_residuals': bao_residuals,
    'env_effect': env_effect,
    'breakpoints': breakpoints,
    'overlap_window': {'z_min': 0.35, 'z_max': 0.80},
    'claim_level': {
        'environment_effect': 'ROBUST',
        'bao_coincidence': 'SUGGESTIVE', 
        'coupling_ordering': 'CONFIRMED_3_PAIRS',
        'paper1_readiness': '80%',
        'paper2_readiness': '60%',
        'paper3_readiness': '20%',
    }
}

with open(f'{RESULTS_DIR}/bao_overlay_analysis.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"Saved to {RESULTS_DIR}/bao_overlay_analysis.json", flush=True)
