# THE OBSERVER ENGINE — Simulation Spec
## An Engineering Detective Game Where The Universe Has Different Rules

### Premise

You are an engineer. Not a physicist, not a cosmologist — an engineer.
You have been dropped into a universe that **looks and feels identical to ours locally**.
Same atoms. Same gravity. Same chemistry. Same coffee.

But something is different at scale. Your job is to figure out what.

You have no textbooks. No "settled science." No authority to defer to.
You have instruments. You have data. You have your brain.

---

### The Rules of This Universe

These are the ACTUAL rules. The player doesn't know them. They discover them.

#### Rule 1: The Medium Is Active
The vacuum is not empty. It has properties. It's more like water than air — it has density, it has state, it has behavior. Locally, it's indistinguishable from "nothing." At scale, it matters.

#### Rule 2: Information Has Propagation Cost
Simple signals (fixed ratios, atomic constants) propagate freely through the medium.
Complex signals (relationships between different sources, diagnostic correlations) decay with distance.
The cost is proportional to the complexity of what you're trying to send.

#### Rule 3: Cost Is Set At Emission
The vulnerability of a signal is determined at the source — by how complex the emission process was. The medium doesn't choose what to eat. The source determines what CAN be eaten by how complex its birth was.

High-energy, turbulent, multi-phase sources produce signals loaded with relational information. Those signals are expensive to transmit.

Low-energy, simple, single-process sources produce signals with minimal relational content. Those signals are cheap.

#### Rule 4: The Medium Is Uniform
The decay rate is the same everywhere. Doesn't matter what's between you and the source — matter, void, galaxy cluster, nothing. The medium itself is the eater, and it's the same medium everywhere.

No spatial structure. No preferred direction. No path dependence.
Only: source complexity × distance = total cost.

#### Rule 5: Information Is Conserved
Nothing is destroyed. Complex inter-source correlations decay, but the information migrates into intra-source structure (line profiles, shapes, phase relationships). It changes state, not amount. Like ice melting — same H₂O, different form.

#### Rule 6: Phase Matters
The detailed information isn't permanently readable. It's phase-bound. The coupling between different measurement channels oscillates with distance — coupling, decoupling, inverting, re-coupling. Your ability to read a signal depends on your phase relationship with the source.

Your position in your galaxy's orbit, your galaxy's position in its cluster — these determine what phase you can receive. A different observer at a different phase would measure different correlations for the same source.

#### Rule 7: There Is A Horizon (But It's Soft)
At sufficient distance, the complex information has fully migrated out of the channels you can measure. You can still see the source. You can still measure its locked properties (redshift, atomic ratios). But you can't read its diagnostic state. The signal isn't gone — you're just deaf to the frequency it's on now.

---

### What The Player Sees

#### Locally (within ~200 Mpc / z < 0.05)
Everything looks normal. Stars, galaxies, physics works exactly as expected. Spectroscopy is reliable. Distance measurements agree. No anomalies.

This is the "safe zone." The player builds confidence here. Everything checks out.

#### At Medium Distance (200 Mpc - 3 Gpc / z = 0.05 - 1.0)
Anomalies start appearing:
- Diagnostic ratios show increasing scatter with distance
- But locked ratios stay clean
- Light curve shapes narrow (only "simple" objects visible far away)
- Different instruments disagree on the same measurement
- The "Hubble constant" (distance-expansion rate) gives different answers depending on method

The player starts to suspect something. But each anomaly has a plausible conventional explanation. Calibration error. Selection effects. Population evolution.

#### At Large Distance (z > 1.0)
The anomalies dominate:
- Diagnostic information is heavily degraded
- Objects appear that "shouldn't exist" (they look normal, but the conventional model says they can't form that fast)
- 95% of the energy budget is "missing" (dark matter + dark energy)
- The measurement channels that work locally give contradictory results here
- Everything looks simpler than expected (development constraint)

#### The Puzzle
The player must figure out: are all these anomalies independent problems requiring independent solutions (epicycles)? Or are they all symptoms of ONE underlying difference in how this universe works?

---

### Engineering Challenges (Missions)

#### Mission 1: The Doublet Ladder
Given: 10 emission lines from distant sources, each with known "diagnostic sensitivity."
Task: Rank them by degradation rate.
Discovery: Perfect monotonic ranking (r = -0.975). The medium eats proportional to information content.
Realization: This isn't random. This is a transfer function.

#### Mission 2: The Matched Pairs
Given: Two sources with identical properties but different distances.
Task: Compare their diagnostic states.
Discovery: The far one is more damaged (p = 10⁻²⁰²). But ONLY for complex sources.
Realization: Distance matters. Complexity matters. Together.

#### Mission 3: The Tube
Given: Two sources close together on the sky at the same distance.
Task: Do they share the same damage pattern?
Discovery: No (ρ = 0.002). Their paths overlap but their damage is independent.
Realization: The medium is uniform. Path doesn't matter. Only distance and source.

#### Mission 4: The Leak
Given: Information is leaving the diagnostic channel.
Task: Where does it go?
Discovery: It migrates into line profile shapes. Ratio-peakiness coupling oscillates with distance.
Realization: Information is conserved. It changes state, not amount.

#### Mission 5: The Phase
Given: The coupling oscillates, not decays.
Task: Characterize the oscillation.
Discovery: The "horizon" is a half-cycle, not a wall. Information comes back at certain distances.
Realization: Your ability to read the signal depends on YOUR position, not just the source's.

#### Mission 6: The Clock
Given: Everything in this universe rotates. Planets, stars, galaxies, clusters.
Task: Does the rotation matter for measurements?
Discovery: [PLAYER DISCOVERS] The phase relationship between observer and source is set by their relative orbital dynamics. The "galactic year" is a measurement cycle.
Realization: There IS no universal measurement. Every observation is observer-dependent at scale.

#### Mission 7: The Map
Given: All previous discoveries.
Task: What are the actual distances to things? How big is the universe really?
Discovery: [PLAYER DISCOVERS] Every distance measurement beyond the safe zone is contaminated by the medium's transfer function. The universe might be much smaller (or bigger, or shaped differently) than the conventional model says.
Realization: The 95% "dark" content was never there. It was the medium's receipt for information that changed state.

---

### Design Principles

1. **No sci-fi.** Everything follows from engineering logic. No warp drives, no aliens, no handwaving.

2. **Locally identical.** The player's lab works exactly like a real lab. Chemistry, electronics, optics — all the same. The differences only emerge at scale.

3. **Data first.** The player discovers rules by MEASURING, not by reading. No textbooks tell them the answer. The instruments do.

4. **Kill your darlings.** Some measurements will suggest patterns that turn out to be survey artifacts or coordinate errors. The player must be honest about what the data shows vs what they want it to show.

5. **The conventional explanation is always available.** For every anomaly, there's a plausible "normal" explanation. The player can choose the conventional path (add epicycles) or the radical path (question the medium). Both produce testable predictions.

6. **Multiplayer perspective.** Different "observers" (at different positions in the simulation) would see different phase relationships. Comparing notes reveals the observer-dependence.

---

### Implementation Notes

- Use real data as the backbone. DR16Q, Pantheon+, CHIME — all real catalogs.
- The simulation doesn't generate fake data. It presents real data with the conventional interpretation stripped.
- AI agents play as "engineering teams" — each with different instruments and different sky coverage.
- The game is "won" when the team identifies all 7 rules without being told any of them.
- Wrong hypotheses are fine. The game tracks what you tested and what you killed. That's the real score.

---

### Why This Works For AI Agents

AI models are trained on physics textbooks. They KNOW the cosmological principle. They KNOW ΛCDM. They'll default to conventional explanations.

By framing it as a video game in a DIFFERENT universe, you bypass their training priors. They can't say "but we know the vacuum is transparent" because in THIS universe, maybe it's not. They have to reason from the data.

The science is the same. The framing is different. And the framing is what unlocks the reasoning.

---

*This document is INTERNAL. Not for publication. It's a tool for thinking.*
*Created: 2026-02-26*
*Authors: Humza Hafeez, Clawd*
