# Closure Theory — Vacuum Membrane Transport Framework
## Working Draft — March 7, 2026

### Core Premise
The vacuum is not passive empty space. It is a **selective information membrane** with transport properties analogous to biological cell membranes. The doublet ladder, sigmoid thresholds, and channel divergence are all natural consequences of membrane transport kinetics.

---

## 1. The Biological Membrane Model

In biological membranes:
- **Passive diffusion**: small, simple molecules pass freely
- **Facilitated/active transport**: complex molecules are selectively filtered
- **Selectivity**: proportional to molecular complexity/information content
- **Saturation**: follows Michaelis-Menten kinetics (sigmoid)
- **Flux vs energy**: membrane affects HOW MANY molecules pass, not their individual kinetic energy

## 2. Translation to Vacuum Membrane

| Biological Membrane | Vacuum Membrane |
|---|---|
| Molecule | Spectral line's encoded information |
| Molecular complexity | Diagnostic sensitivity q |
| Membrane thickness | Comoving path length (∝ z, DM) |
| Transport rate | Information degradation rate |
| Flux (particle count) | Photon flux in line |
| Kinetic energy | Line width (velocity dispersion σ) |
| Partition coefficient K | Coupling strength Γ × q |
| Non-specific leak | Thermal vacuum noise floor |

## 3. The Transport Equation

For spectral line i with diagnostic sensitivity q_i, the correlation remaining after path length z:

**C_i(z) = C₀_i / (1 + exp(Γ × q_i × (z - z₀_i) / w))**

Where:
- C₀_i = intrinsic correlation at source (z=0)
- Γ = universal membrane coupling constant
- q_i = diagnostic sensitivity (from atomic physics: CHIANTI/NIST)
- z₀_i = midpoint distance (where C = C₀/2)
- w = transition width

### Michaelis-Menten form:
The rate of information loss for channel i:

**V_i = V_max × q_i × z / (K_m + z)**

- V_max = max degradation rate (membrane property)
- K_m = half-saturation path length
- z₀ = K_m / q_i (higher q → earlier threshold)

## 4. Deriving the Numbers

### 4.1 The Doublet Ladder: r = -0.975

**Biological parallel**: Membrane selectivity coefficient.

In a perfect membrane, degradation D_i ∝ q_i exactly → r = -1.000.

Departure comes from **non-specific permeability** (the membrane's thermal noise floor):

D_i = Γ × q_i + ε_i

where ε_i is non-specific leak with variance σ²_ε.

The correlation:
r = -1 / √(1 + σ²_ε / (Γ² × σ²_q))

For r = -0.975:
σ²_ε / (Γ² × σ²_q) = 1/0.975² - 1 = 0.0519

**Signal-to-noise ratio of the membrane selectivity = 4.39:1**
**Non-specific leak = 5.2% of total variance**

In biological membranes, non-specific leak current is typically **2-5%** of total conductance. The vacuum membrane has the SAME order of magnitude. This is either a deep universality or needs explanation.

### 4.2 Flux vs Sigma: r = -0.943 vs r = +0.143

**Biological parallel**: Membranes affect particle COUNT, not particle ENERGY.

A K⁺ ion that passes through a potassium channel has the same kinetic energy on both sides. But the NUMBER of ions passing is controlled by the channel.

Translation:
- **Flux** = number of information-carrying photons → controlled by membrane → DEGRADES (r = -0.943)
- **Sigma** = velocity structure of emitting gas → property of the SOURCE, not the path → PRESERVED (r = +0.143 ≈ 0, within noise)

The membrane operates on the **statistical ensemble** of photons (how many carry state information), not on individual photon properties (energy, wavelength). This is why wavelength doesn't shift — the membrane doesn't interact with E = hν. It interacts with the CORRELATION between photon populations.

### 4.3 Sigmoid Thresholds: z₀ = 0.82 (SNe), 1.05-1.23 (Quasars), DM≈500 (FRBs)

**Biological parallel**: Different substrates have different K_m values.

The midpoint z₀ depends on:
- **Source baseline correlation C₀**: stronger initial correlation → takes more path to degrade → higher z₀
- **Effective q of the source class**: higher mean q → lower z₀

SNe Ia color-distance coupling has moderate initial correlation and moderate q → z₀ = 0.82
Quasar EW coupling has strong initial correlation → z₀ = 1.05-1.23 (needs more path to break)
FRBs have weak initial coupling → DM threshold is lower (in equivalent path units)

From Michaelis-Menten: **z₀ = K_m × C₀ / (V_max × q_eff)**

If K_m and V_max are universal (membrane properties), then:
z₀_quasar / z₀_SN = (C₀_quasar × q_SN) / (C₀_SN × q_quasar)

This ratio should be calculable from the intrinsic properties of each source class.

### 4.4 The -0.873 Constant (Quasar IP Correlation)

**Biological parallel**: The Hill coefficient.

In cooperative membrane transport, binding follows the Hill equation:
θ = [L]^n / (K_d^n + [L]^n)

where n = Hill coefficient measures cooperativity.

For n = 1: independent binding (Michaelis-Menten)
For n > 1: positive cooperativity
For n < 1: negative cooperativity

The -0.873 across 7 quasar ionization lines...

If the 7 lines probe the membrane at different "depths" (ionization potential = depth), the cross-correlation between adjacent layers follows:

r_cross = -cos(π × (N-1)/N) for N stratified layers

For N = 7: r = -cos(6π/7) = -cos(154.3°) = **-(-0.901)** = needs adjustment

Alternative: if the membrane has **3 independent transport channels** (T, n_e, Z — the three thermodynamic axes that q responds to), and each line samples a weighted combination:

The maximum anti-correlation in 3D state space between the "most sensitive" and "least sensitive" directions:

r = -cos(arctan(√2)) = -cos(54.7°) = **-0.578** — that's the magic angle again but wrong value.

Hmm. Let me try: for a 3D system with cooperative binding across channels:

r = -(1 - 1/N²) where N = number of independent membrane channels?

For N = 3.5 (fractional cooperativity): r = -(1 - 1/12.25) = -0.918 — closer but not exact.

**THIS NEEDS MORE WORK. The -0.873 derivation is not yet clean.**

### 4.5 The 0.533 Universal Coupling Constant

If Γ × V_max = 0.533, this sets the membrane's intrinsic transport rate.

In biological terms, this would be the membrane's **permeability coefficient P**:

P = D × K / Δx

where D = diffusion coefficient, K = partition coefficient, Δx = membrane thickness.

For the vacuum membrane:
- D = related to the speed of information propagation (c?)
- K = partition between "encoded" and "free" information
- Δx = "thickness" of the vacuum membrane per unit comoving distance

**0.533 ≈ 1/√(2π × e^(1/2))** ??? No obvious fundamental constant decomposition yet.

---

## 5. Predictions from the Membrane Model

### 5.1 Kill-shot test: GW vs EM distances
Gravitational waves don't carry atomic transition information → they don't couple to the membrane's selectivity gradient → GW distances should EXCEED EM distances beyond z₀. The discrepancy should scale as:

Δd/d = Γ × q_eff × sigmoid(z, z₀)

### 5.2 Laboratory null test
The membrane effect requires cosmological path lengths. Laboratory atomic spectroscopy should show ZERO correlation degradation regardless of diagnostic sensitivity. (Trivially satisfied but important as a control.)

### 5.3 Gravitational lensing split test
Two images of the same lensed quasar traverse different path lengths. The longer-path image should show slightly MORE correlation degradation in diagnostic lines, with the difference proportional to Δz_path × q.

### 5.4 Line-by-line prediction
Given the q values from CHIANTI for any NEW emission line not in the original 6, the model predicts its degradation rate: D_new = Γ × q_new. This is testable with [SIII] 6312, [NII] 5755, [ArIII], etc.

### 5.5 Biological membrane constant comparison
If the vacuum membrane's non-specific leak (5.2%) is related to a fundamental constant, it should appear in other information-processing boundaries. Check: is there a 2-5% noise floor in other cosmological measurements that has been attributed to "intrinsic scatter"?

---

## 6. What IS the Membrane?

The membrane is not a physical barrier in space. It is a **property of how information degrades across expanding spacetime**.

Possible identifications:
1. **The expansion itself** — as space expands, the phase coherence between correlated photon populations degrades. The rate of decoherence depends on how much "state information" the correlation encodes (= q).

2. **Vacuum energy as membrane** — the cosmological constant Λ isn't just accelerating expansion; it's the "lipid bilayer" of the cosmic membrane. Its energy density sets the transport rate.

3. **Holographic boundary** — if information in a volume is encoded on its boundary (holographic principle), then the "membrane" is literal: it's the boundary of the observable universe, and the degradation is the cost of encoding diagnostic information on a surface that's stretching.

---

## Status
- Doublet ladder (r = -0.975): DERIVABLE from membrane noise model ✓
- Flux vs sigma divergence: EXPLAINED by count-vs-energy separation ✓  
- Sigmoid thresholds: EXPLAINED by Michaelis-Menten saturation ✓
- Different z₀ across sources: EXPLAINED by different C₀/q_eff ✓
- Isotropic: EXPLAINED if membrane = property of spacetime ✓
- -0.873 constant: NEEDS MORE WORK
- 0.533 coupling: NEEDS MORE WORK
- Specific z₀ values: NEED CALCULATION from source properties
