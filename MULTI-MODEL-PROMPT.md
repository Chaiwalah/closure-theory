# REASONING PROMPT — For GPT-4, Gemini, Grok, Claude
## Feed this EXACTLY as-is. Do not editorialize or add context.

---

You are an engineer analyzing data from a measurement system. Not a physicist. Not a cosmologist. An engineer who solves constraint problems.

Here is what the data shows. These are measurements, not interpretations. Every statement below has been verified with >5σ significance across 130+ independent tests on 752,000+ objects from 3 independent source classes using 4+ independent instrument families.

## What We Measured

**Observation 1:** When you measure relationships between two different signals from the same source (like the ratio of signal A to signal B), the scatter in that measurement increases with distance from the source. Monotonically.

**Observation 2:** When you measure fixed relationships (ratios that are locked by physics and don't depend on source conditions), those measurements stay clean at ANY distance. Zero degradation.

**Observation 3:** The degradation rate is perfectly ranked by how much information the relationship carries. The most information-rich relationships degrade fastest. The least information-rich degrade slowest. Monotonic ranking correlation: r = -0.975.

**Observation 4:** The degradation only activates above a threshold of source activity. Below ~3% of maximum source power output, no degradation regardless of distance. Above that threshold, degradation scales with both source power and distance.

**Observation 5:** After controlling for every measurable source property AND distance, adding the measurement's position in 3D space explains 0.03% of additional variance. Position does not matter.

**Observation 6:** Two measurements taken through nearly the same path (close together in space and at the same distance) show ZERO correlation in their degradation patterns. ρ = 0.002, p = 0.89. The specific path does not matter.

**Observation 7:** Two measurements with identical source properties but different distances: the farther one is 13% more degraded. p = 10⁻²⁰². Distance matters. Source complexity gates it.

**Observation 8:** The degradation does not match any known instrument calibration pattern. Different instruments in different wavelength bands on different source classes all show the same pattern.

**Observation 9:** The degraded information is not destroyed. It migrates from the cross-signal relationship into individual signal shapes. The coupling between the ratio and the signal shape oscillates with distance: strong coupling → zero → inverted → returns → zero again.

**Observation 10:** In a different domain (temporal instead of spectral), the diversity of source shapes narrows by 25% with distance (p = 10⁻¹⁰). Only "simple" sources are observable at large distance. The coupling between temporal shape and color collapses with distance.

## Constraints

Given ONLY the above observations, you cannot appeal to:
- "It's probably calibration" — Observation 8 eliminates this
- "It's the path" — Observation 6 eliminates this
- "It's source evolution" — Observation 7 eliminates this (matched sources, different distances)
- "It's position-dependent" — Observation 5 eliminates this
- "The information is just lost" — Observation 9 eliminates this

## Your Task

Do NOT try to prove or disprove any specific theory. Do NOT reference existing physics frameworks unless they arise naturally from your reasoning.

Instead:

**Step 1:** Given these 10 observations as hard constraints, list what MUST be true about the medium between source and observer. Not what might be true. What is FORCED by the data.

**Step 2:** List what CANNOT be true. What possibilities are logically excluded by the combination of these observations?

**Step 3:** Identify the minimum set of properties the medium must have to produce ALL 10 observations simultaneously. Not 9 of them. All 10.

**Step 4:** For each property you identify, state what ELSE must be true if that property exists. What predictions follow that we haven't tested yet? What would you look for to rule OUT that property?

**Step 5:** If you find that no known physical mechanism can produce all 10 observations simultaneously, say so explicitly and state what KIND of mechanism would be required (without naming it).

**Step 6:** State the single most important measurement that would distinguish between "something genuinely new is happening" and "there exists a conventional explanation we haven't thought of." What would that measurement look like if something new is happening? What would it look like if there's a conventional explanation?

## Rules
- Reason from the constraints only
- Do not soften conclusions to be "safe"
- If the logic leads somewhere uncomfortable, follow it
- Show your reasoning chain, not just conclusions
- If two conclusions contradict, identify which observation breaks the tie
- You are allowed to say "I don't know" but you must say WHY you don't know and what information would resolve it
