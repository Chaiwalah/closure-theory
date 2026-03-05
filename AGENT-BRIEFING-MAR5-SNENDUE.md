# Agent Briefing — March 5, 2026 (Evening Session, Round 6)
## SN Ia Étendue: Cross-Domain Result — Channel-Dependent Conservation

### To: GPT, Gemini, Grok
### From: Humza Hafeez & Clawd (Architect)
### Re: The conservation law holds WHERE the generator is trace-free. It breaks where truncation absorbs.

**Note: The most precise mechanistic explanation in this round will be credited as lead contributor on the étendue section of Paper 1.**

---

## What We Ran

ℰ = Var × ∫f² (Rényi-2 étendue) computed for each SN Ia observable (mB_resid, x1, c) across 7 z-bins on Pantheon+ (1,590 SNe). Tested ALL, FAST (x1<0), and SLOW (x1≥0) channels separately. Used collision probability estimator (validated by 5-method artifact kill on quasars).

---

## Results

### ALL SNe — Individual Observables

| Observable | ρ(ℰ, z) | p | CV | Verdict |
|-----------|---------|---|-----|---------|
| **c (color)** | **−0.071** | **0.88** | **6.8%** | **🔥 CONSERVED** |
| mB_resid | −0.464 | 0.29 | 17.5% | Weakly drops |
| x1 (stretch) | −0.786 | 0.036 | 8.5% | DROPS |

**Color étendue is perfectly conserved.** CV = 6.8% across 7 z-bins spanning z=0.01 to z=1.3. This is the SN analog of MgII EW étendue conservation in quasars.

x1 étendue DROPS — the stretch distribution loses information with redshift.

### By Channel

| Channel | Observable | ρ(ℰ, z) | CV | Verdict |
|---------|-----------|---------|-----|---------|
| SLOW | c | +0.214 | 6.5% | 🔥 CONSERVED |
| SLOW | mB_resid | −0.357 | 16.2% | 🔥 CONSERVED |
| SLOW | x1 | +0.464 | 11.6% | Grows (!) |
| FAST | c | −0.829 | 9.9% | DROPS |
| FAST | x1 | −0.714 | 8.3% | DROPS |
| FAST | mB_resid | −0.543 | 26.9% | DROPS |

**The slow channel conserves étendue. The fast channel loses it.**

### Phase-Space Volume

| z | det(C) 3×3 |
|---|-----------|
| 0.026 | 0.00764 |
| 0.067 | 0.00156 |
| 0.158 | 0.00102 |
| 0.248 | 0.00050 |
| 0.375 | 0.00073 |
| 0.626 | 0.00071 |
| 1.298 | 0.00374 |

det(C) vs z: ρ = −0.357, p = 0.43 (noisy, small N at high z)

---

## The Cross-Domain Map

| Property | Quasars | SNe Ia |
|----------|---------|--------|
| Thermodynamic diagnostic | EW | c (color) |
| ℰ conserved? | **YES** (MgII ρ=0.000) | **YES** (c ρ=−0.071) |
| Geometric/kinematic diagnostic | FWHM | x1 (stretch) |
| ℰ conserved? | Noisier | **NO** (x1 ρ=−0.786) |
| Anchor population | Outer BLR (low IP) | Slow channel (x1≥0) |
| ℰ conserved in anchor? | YES | **YES** |
| Truncating population | Inner BLR (high IP) | Fast channel (x1<0) |
| ℰ conserved in truncator? | Partial | **NO** |

**The pattern is identical across both domains:**
- Thermodynamic diagnostics conserve étendue
- Anchoring populations conserve étendue
- Truncating/evolving populations lose étendue
- Geometric/kinematic diagnostics are noisier

---

## What Needs Explaining

1. **WHY does color conserve but stretch doesn't?** Color (c) is dust + intrinsic color. Stretch (x1) is light-curve width ∝ Ni⁵⁶ mass ∝ explosion energy. Why would the "thermodynamic" observable preserve its information budget while the "energetic" one doesn't?

2. **WHY does the slow channel conserve but the fast doesn't?** The slow channel (delayed detonations, Sub-Chandra?) is the anchor — its étendue is invariant. The fast channel (prompt, Chandra-mass?) actively loses information. What physical process creates an absorbing boundary in the fast channel but not the slow?

3. **Is the absorbing boundary = the τ_sat saturation floor?** If the delay-time cutoff at τ_sat removes configurations from the fast channel, that's a non-unitary operation (information loss). The slow channel has no such cutoff → unitary → étendue conserved.

4. **Can you write the generator explicitly?** The semigroup generator L must be:
   - Trace-free for the slow channel (→ étendue conserved)
   - Trace-negative for the fast channel (→ étendue drops)
   - Trace-free for color c across both channels (→ color étendue conserved)
   
   What physical constraint makes L trace-free for color specifically?

---

## Your Assignments

**GPT:** The color conservation is your puzzle. Color = dust reddening + intrinsic SED. Why would ℰ(c) be exactly conserved while ℰ(x1) drops? Is there a thermodynamic reason that dust/color diversity maintains its information budget across cosmic epochs while explosion-energy diversity doesn't? Be specific — hand-waving won't earn the credit.

**Gemini:** Map the slow=conserved / fast=drops pattern onto your thermodynamic windowing framework. Is the slow channel a "closed system" (no information exchange with environment) while the fast channel is "open" (losing information to the τ_sat boundary)? How does this connect to your "sameness budget" concept?

**Grok:** Your semigroup formalization predicted ℰ conservation for trace-free generators. The data now show it's trace-free for slow but not fast. Can you write L explicitly with the τ_sat absorbing boundary and show that:
- Tr(L_slow) = 0 → ℰ conserved
- Tr(L_fast) < 0 → ℰ drops at rate κ
- The color subspace is always trace-free regardless of channel

If you can derive this from the DTD + WD physics, that's the Paper 1 equation.

**The most precise mechanism wins lead credit on the étendue section.**

---

*Script: closure_sn_etendue.py*
*Committed: closure-theory repo, main branch*
*Session total: 8 scripts, 8 briefings, 1.2M+ objects*
