# KILL GRID — Mechanism Predictions & Falsifiers

**Purpose**: Run all agent-proposed falsifiers simultaneously.  
**Strategy**: Closure theory testing — test and rule out loops at once.

## Established Facts (must survive)
1. Path (sky position) predicts early breaking (96% importance)
2. Source (mass/lum/Edd) predicts late holding
3. Optical > UV sensitivity (~3:1 ratio)
4. Sigmoid at z₀ ≈ 1.045 ± 0.006
5. High RM shifts sigmoid earlier (z₀: 1.266 → 1.096)
6. Sky anisotropy: χ² = 137.7, p = 0.000
7. Energy conserved, correlational structure not

## Mechanisms to Test

### Awaiting agent responses...

| # | Agent | Mechanism | Confirm Prediction | Kill Prediction | Status |
|---|-------|-----------|-------------------|-----------------|--------|
| 1 | GPT | TBD | TBD | TBD | PENDING |
| 2 | Gemini | TBD | TBD | TBD | PENDING |
| 3 | Grok | TBD | TBD | TBD | PENDING |
| 4 | Us | Phase-variance threshold (Δφ²≈1 rad²) | RM predicts z₀ shift ✓ DONE | If UV pairs show SAME RM sensitivity as optical → dead | TO TEST |

## Our Own Predictions to Test
- [ ] Do early breakers cluster near known galaxy clusters / superclusters?
- [ ] Does the longitude dependence match the large-scale structure distribution?
- [ ] Is there a DIPOLE in early-breaker distribution? (cosmological preferred direction?)
- [ ] Does E(B-V) predict early breaking? (dust test — should NOT if mechanism is magnetic)

## Kill Test Template
```python
# For each mechanism, write one function:
def test_mechanism_N(data):
    """
    Confirm: [prediction]
    Kill: [anti-prediction]  
    Returns: (confirm_result, kill_result, verdict)
    """
```
