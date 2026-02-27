# Literature Sweep: Grokking Experiments for Closure Theory Validation

*Compiled 2026-02-15*

---

## 1. Power et al. (2022) — Original Grokking Paper

**Title:** "Grokking: Generalization Beyond Overfitting on Small Algorithmic Datasets"  
**ArXiv:** https://arxiv.org/abs/2201.02177  
**GitHub:** No official code release (OpenAI internal)

### Key Results
- **Architecture:** Small decoder-only transformer (2 layers, 4 heads, d_model=128)
- **Tasks:** Binary operations on Z_p (modular arithmetic), trained on various operations (addition, subtraction, division, permutation compositions)
- **Grokking observed:** Training accuracy hits ~100% quickly, test accuracy stays at chance for thousands of epochs, then suddenly jumps to ~100%
- **Hyperparameter sweeps:** Varied data fraction (percentage of operation table used for training). Grokking is most prominent at intermediate data fractions (~30-50%). Below a critical fraction: no generalization. Above: fast generalization without delay.
- **Weight decay critical:** Grokking requires regularization (weight decay). Without it, the model memorizes and never generalizes
- **No width sweep reported.** They used a single architecture. Author comment on Reddit: "I did not experiment on a full sweep of network sizes, but I trained a much smaller network (1 layer, 1 head, d_model=32) which never overfit."

### Relevance to Closure Theory
- ❌ No systematic width/capacity variation
- ❌ No internal representation metrics logged
- ✅ Data fraction sweeps could be reinterpreted through closure lens (more data → lower effective dimensionality of task relative to model → changes closure ratio)
- ⚠️ The observation that d_model=32 "never overfit" is interesting — may be a capacity threshold effect interpretable via closure

---

## 2. Liu et al. (2023) — Omnigrok (ICLR 2023 Spotlight)

**Title:** "Omnigrok: Grokking Beyond Algorithmic Data"  
**ArXiv:** https://arxiv.org/abs/2210.01117  
**OpenReview:** https://openreview.net/forum?id=zDiHoIWa0q1  
**GitHub:** https://github.com/KindXiaoming/Omnigrok ✅ Full code available

### Key Results
- **Core insight:** Grokking is caused by mismatch between training and test loss landscapes as a function of weight norm. There exists a "Goldilocks zone" of weight norms where generalization occurs.
- **Datasets:** Modular addition, teacher-student, MNIST, IMDb, QM9 molecules
- **Two experiment types per dataset:** (1) Reduced landscape analysis (fixed weight norm), (2) Standard grokking experiments
- **Weight initialization scale (α) is critical:** Large α → overfitting first, then slow weight decay brings weights to Goldilocks zone → grokking. Small α → direct generalization.
- **Data size effects:** Larger datasets make the Goldilocks zone broader, reducing or eliminating grokking
- **Critical data size:** Below a threshold, grokking is impossible regardless of weight norm

### Width/Capacity Experiments
- The modular addition setup uses an MLP decoder: `Dec(E_i + E_j)` with embeddings. The MLP architecture is configurable.
- **They do NOT systematically vary network width** in the main paper. The focus is on weight norm, not architecture capacity.
- The code repo has separate folders per dataset with configurable architectures — could be used for width sweeps.

### Relevance to Closure Theory
- ✅ **HIGHLY RELEVANT.** Their "Goldilocks zone" concept maps directly to closure: the weight norm regime where the network's effective representation capacity matches the task complexity.
- ✅ Code is modular and minimal — easily extensible for width sweep experiments
- ✅ The reduced landscape analysis (fixing weight norm) isolates the representation quality from optimization dynamics — could directly measure closure ratio at different scales
- ⚠️ Missing: they don't track effective rank, activation covariance, or singular values during training
- 🔑 **Key reanalysis opportunity:** Run their code with varying widths while logging activation covariance spectrum → test if grokking time correlates with effective bandwidth rather than parameter count

---

## 3. Nanda et al. (2023) — Progress Measures for Grokking (ICLR 2023)

**Title:** "Progress measures for grokking via mechanistic interpretability"  
**ArXiv:** https://arxiv.org/abs/2301.05217  
**Website:** https://www.neelnanda.io/grokking-paper  
**Code:** Available on the paper website  
**Replication walkthrough:** https://www.neelnanda.io/mechanistic-interpretability/modular-addition-walkthrough

### Key Results
- **Architecture:** 1-layer transformer on modular addition (mod 113)
- **Mechanistic finding:** The model learns a Fourier-based algorithm: embeds inputs as rotations, composes them, reads off logits via cos(ω(a+b-c))
- **Progress measures tracked:**
  - Fourier component magnitudes in embeddings and weights
  - "Excluded loss" (loss on held-out data)
  - Weight norms
  - Logit structure evolution
- **Key finding:** Progress measures increase "relatively smoothly" before the phase transition, but they "lack a general notion of criticality that would allow prediction of when the phase transition will happen ex ante"
- **No width sweep.** Single architecture studied.

### Internal Representation Metrics
- ✅ Fourier component evolution (could be reinterpreted as spectral structure of representations)
- ✅ Weight norm tracking
- ❌ No activation covariance or effective rank
- ❌ No mutual information measures
- The code allows training with checkpoints — could add custom metrics

### Relevance to Closure Theory
- ✅ **The Fourier decomposition IS effectively a bandwidth measure.** The number of active Fourier components at any training step is analogous to the effective rank of the representation
- ✅ Their admission that they can't predict grokking timing is exactly the gap closure theory aims to fill
- 🔑 **Key reanalysis:** Track the effective rank of embedding/activation covariance during their training runs. Hypothesis: grokking occurs when effective rank of representation matches the intrinsic dimensionality of the task (which for mod-p addition is related to the number of required Fourier modes)
- ⚠️ Only 1-layer transformer — may be too simple for rich closure dynamics

---

## 4. Liu, Kitouni & Tegmark (2022) — Effective Theory of Representation Learning (NeurIPS 2022)

**Title:** "Towards Understanding Grokking: An Effective Theory of Representation Learning"  
**NeurIPS:** https://proceedings.neurips.cc/paper/2022/hash/dfc310e81992d2e4cedc09ac47eff13e-Abstract-Conference.html  
**ArXiv:** https://arxiv.org/abs/2205.10343

### Key Results
- **Develops an effective theory** predicting grokking dynamics in a toy setting
- **Key prediction:** Grokking time t ~ 1/λ₃, where λ₃ is a representation quality eigenvalue that depends on data fraction
- **Phase diagrams:** Maps (data fraction, weight decay, learning rate) → {memorization, grokking, generalization, confusion}
- **Critical data fraction ~0.4** for modular addition (above this, grokking time decreases with more data)
- **Representation Quality Index (RQI):** They define a metric for how well representations support the task

### Relevance to Closure Theory
- ✅ **VERY RELEVANT.** Their λ₃ eigenvalue is essentially a measure of representation bandwidth
- ✅ t ~ 1/λ₃ is a scaling law for grokking time that could be rewritten in closure terms
- ✅ Phase diagrams across hyperparameters are directly testable against closure predictions
- 🔑 **Key connection:** If λ₃ ∝ effective bandwidth / task complexity, then their scaling law IS a closure scaling law. Need to verify if λ₃ scales with width in the predicted way.
- ⚠️ Only validated in toy setting

---

## 5. Thilak et al. (2022) — Slingshot Mechanism

**Title:** "The Slingshot Mechanism: An Empirical Study of Adaptive Optimizers and the Grokking Phenomenon"  
**ArXiv:** https://arxiv.org/abs/2206.04817  
**Apple ML Research:** https://machinelearning.apple.com/research/slingshot-effect

### Key Results
- **Discovery:** Cyclic phase transitions between stable/unstable training regimes in Adam optimizer
- **The "slingshot":** Periodic oscillations in last-layer weight norms, caused by adaptive optimizer dynamics
- **Key finding:** "Without explicit regularization, grokking almost exclusively happens at the onset of slingshots"
- **Implication:** Grokking may be an optimization artifact of adaptive gradient methods, not purely a representation learning phenomenon

### Relevance to Closure Theory
- ⚠️ **COMPLICATING FACTOR.** If grokking is tied to optimizer dynamics (slingshots), then closure-based predictions need to account for this
- ✅ However, the slingshot mechanism could be the *mechanism* by which the network explores weight space to find the closure point — the slingshot pushes the network from a memorizing configuration to a generalizing one
- ❌ No width/capacity sweeps
- ❌ No representation metrics beyond weight norms
- 🔑 **Test:** Does the slingshot timing correlate with effective bandwidth metrics? If closure theory is correct, slingshots that lead to grokking should coincide with representations reaching sufficient effective rank

---

## 6. Davies, Langosco & Krueger (2023) — Unifying Grokking and Double Descent

**Title:** "Unifying Grokking and Double Descent"  
**ArXiv:** https://arxiv.org/abs/2303.06173 (NeurIPS 2022 ML Safety Workshop)

### Key Results
- **Framework:** Both grokking and double descent arise from different "pattern learning speeds" — the model learns some patterns fast (memorization) and others slow (generalization)
- **Three types of double descent patterns identified**
- **Hypothesis extends to model capacity:** When varying model capacity instead of training time, the same framework applies
- **"Predictiveness" measure:** Tracks how much each pattern type contributes to predictions over training

### Relevance to Closure Theory
- ✅ **VERY RELEVANT.** Their "pattern learning speeds" framework is compatible with closure theory — faster patterns have higher bandwidth relative to their complexity, slower patterns have lower bandwidth
- ✅ The capacity variation prediction is directly testable: closure theory would predict specific scaling of the "crossing point" with width
- ⚠️ Mostly conceptual — limited experimental validation
- 🔑 **Key connection:** If "learning speed of pattern i" ∝ bandwidth allocated to pattern i / complexity of pattern i, then their framework becomes a special case of closure theory

---

## 7. Merrill, Tsilivis & Shukla (2023) — A Tale of Two Circuits

**Title:** "A Tale of Two Circuits: Grokking as Competition of Sparse and Dense Subnetworks"  
**ArXiv:** https://arxiv.org/abs/2303.11873  
**GitHub:** https://github.com/Tsili42/parity-nn (minimal parity learning codebase)

### Key Results
- **Model:** MLPs on sparse parity learning tasks
- **Core finding:** Two competing circuits during grokking:
  - **Dense (memorizing) circuit:** Uses many parameters, learns fast, poor generalization
  - **Sparse (generalizing) circuit:** Uses few parameters, learns slow, good generalization
- **Weight decay** gradually favors the sparse circuit by penalizing the dense one
- **Sparsity metrics tracked during training**

### Relevance to Closure Theory
- ✅ **STRONGLY SUPPORTIVE.** The sparse vs. dense circuit competition maps to closure: the dense circuit has high bandwidth but low efficiency (bandwidth >> task complexity), while the sparse circuit achieves bandwidth ≈ task complexity (closure)
- ✅ Weight decay drives the system toward closure by penalizing excess bandwidth
- ✅ Code available for parity tasks
- 🔑 **Key reanalysis:** Measure effective rank of each circuit's activations. Hypothesis: grokking = transition point where sparse circuit's effective bandwidth first exceeds task dimensionality
- ⚠️ Parity tasks, not modular arithmetic — different but complementary

---

## 8. Varma et al. (2023) — Explaining Grokking Through Circuit Efficiency

**Title:** "Explaining grokking through circuit efficiency"  
**ArXiv:** https://arxiv.org/abs/2309.02390

### Key Results
- **Framework:** Grokking occurs when the generalizing circuit is more weight-efficient than the memorizing circuit, and weight decay selects for efficiency
- **Novel phenomena discovered:** "ungrokking" (regression from perfect to low test accuracy) and "semi-grokking" (delayed partial generalization)
- **Three necessary conditions for grokking:**
  1. Generalizing solution must be more parameter-efficient
  2. Memorizing solution must be learned faster
  3. Regularization must favor efficiency

### Relevance to Closure Theory
- ✅ "Weight efficiency" ≈ closure ratio. The generalizing circuit achieves closure (bandwidth matches task) while the memorizing circuit is over-parameterized for the task
- ✅ Ungrokking and semi-grokking provide additional test cases for closure predictions
- 🔑 **Prediction from closure theory:** Ungrokking should occur when perturbations push effective bandwidth below the task requirement

---

## 9. Tian (2025) — Provable Scaling Laws of Feature Emergence

**Title:** "Provable Scaling Laws of Feature Emergence from Learning Dynamics of Grokking"  
**ArXiv:** https://arxiv.org/abs/2509.21519  
**GitHub:** https://github.com/yuandong-tian/understanding/tree/main/ssl/real-dataset/cogo ✅

### Key Results
- **Li2 framework:** Three stages: (I) Lazy learning, (II) Independent feature learning, (III) Interactive feature learning
- **Energy function ℰ:** The dynamics of independent feature learning follows gradient ascent of ℰ, whose local maxima are the emerging features
- **Provable scaling laws:** Relates feature emergence, memorization, and generalization to weight decay, learning rate, sample sizes
- **Key hyperparameter roles characterized mathematically**
- **2-layer nonlinear networks on group arithmetic tasks**

### Relevance to Closure Theory
- ✅ **EXTREMELY RELEVANT.** Their energy function ℰ essentially measures the "nonlinear CCA" between input and target — this is closely related to mutual information and thus to closure
- ✅ Their scaling laws could potentially be re-derived from closure theory
- ✅ Three-stage framework maps to closure: Stage I (memorization, bandwidth >> task), Stage II (features emerge, bandwidth approaches task needs), Stage III (refinement, closure achieved)
- ✅ Code available with reproducible experiments
- 🔑 **Critical test:** Their ℰ function's relationship to effective bandwidth needs to be formalized. If ℰ ∝ effective bandwidth measure, then their scaling laws ARE closure scaling laws.

---

## 10. Clauw, Stramaglia & Marinazzo (2024) — Information-Theoretic Progress Measures

**Title:** "Information-Theoretic Progress Measures reveal Grokking is an Emergent Phase Transition"  
**ArXiv:** https://arxiv.org/abs/2408.08944

### Key Results
- **Uses higher-order mutual information** (O-Information) to decompose neuron interactions into synergy and redundancy
- **Three phases identified before grokking:**
  1. Feature learning (low-level patterns)
  2. Emergence (generalizing sub-network forms via synergy)
  3. Decoupling (compression)
- **Synergy peaks predict grokking onset**
- **Weight decay and initialization affect the emergence phase timing**

### Relevance to Closure Theory
- ✅ **DIRECTLY RELEVANT.** This is the closest existing work to measuring what closure theory predicts
- ✅ Synergy = collective bandwidth that emerges from neuron interaction. Redundancy = wasted bandwidth. The transition from redundant to synergistic representations IS the approach to closure
- ✅ Three phases map perfectly: Feature learning (building bandwidth), Emergence (bandwidth reaches task complexity = closure), Decoupling (pruning excess = tightening closure)
- 🔑 **Key connection:** Their synergy metric may be mathematically related to the closure ratio. If synergy/(synergy+redundancy) ≈ closure ratio, this provides direct empirical validation
- ⚠️ Fully connected networks on modular arithmetic only

---

## 11. Yunis et al. (2024) — Spectral Dynamics of Weights

**Title:** "Approaching Deep Learning through the Spectral Dynamics of Weights"  
**ArXiv:** https://arxiv.org/abs/2408.11804

### Key Results
- **Tracks singular value evolution of weight matrices during training**
- **Key finding:** "Grokking is closely tied to rank minimization" — validation loss drops coincide with discovery of low-rank weight solutions
- **Singular value plots during grokking** show progressive rank reduction
- **Alignment metrics** between different weight matrices tracked

### Relevance to Closure Theory
- ✅ **STRONG EVIDENCE FOR CLOSURE.** Rank minimization during grokking = the network reducing its effective bandwidth to match the task. This IS closure.
- ✅ They have actual singular value evolution plots — these could be reanalyzed to compute effective rank trajectories and test closure predictions
- 🔑 **Key test:** Does the effective rank at grokking onset match the intrinsic task dimensionality? Closure theory predicts yes.

---

## 12. Shwartz-Ziv & Tishby (2017) — Information Bottleneck in Deep Learning

**Title:** "Opening the Black Box of Deep Neural Networks via Information"  
**ArXiv:** https://arxiv.org/abs/1703.00810

### Key Results
- **Two training phases:**
  1. **Fitting phase:** I(X;T) and I(T;Y) both increase (network learns representations)
  2. **Compression phase:** I(X;T) decreases while I(T;Y) remains high (network discards irrelevant info)
- **Information plane:** Plots I(X;T) vs I(T;Y) for each layer during training
- **Controversial:** Later work (Saxe et al. 2018, "On the Information Bottleneck Theory of Deep Learning") argued compression doesn't occur with ReLU activations, only with saturating nonlinearities

### Connection to Closure Theory
- ✅ The compression phase IS the approach to closure: the network reduces its effective bandwidth (I(X;T) decreases) while maintaining task-relevant capacity (I(T;Y) stays high)
- ✅ Closure point = the information plane location where I(T;Y) is sufficient for the task AND I(X;T) is minimal
- ⚠️ The controversy about ReLU is important — closure theory needs to work for all activation functions, not just saturating ones
- 🔑 **Reconciliation hypothesis:** Closure doesn't require information compression in the Tishby sense. Instead, it requires that the effective bandwidth (measured by activation covariance rank, not MI) matches task complexity. This may occur even when I(X;T) doesn't decrease (as in ReLU networks where effective rank decreases without MI decrease)

---

## 13. Related: Explaining Grokking and Information Bottleneck through Neural Collapse (2025)

**Title:** "Explaining Grokking and Information Bottleneck through Neural Collapse Emergence"  
**ArXiv:** https://arxiv.org/abs/2509.20829

### Key Results
- **Unifies grokking and information bottleneck** through neural collapse geometry
- **Key measure:** Contraction of population within-class variance drives both grokking and information bottleneck
- **Neural collapse metric** relates to representation geometry

### Relevance to Closure Theory
- ✅ Within-class variance contraction = bandwidth focusing on task-relevant directions = approach to closure
- ✅ Provides geometric interpretation compatible with closure

---

## Summary: Reanalysis Opportunities

### Tier 1: Highest Value (code available, directly testable)

| Paper | What to measure | Code |
|-------|----------------|------|
| **Omnigrok (Liu 2023)** | Run width sweeps, log activation covariance eigenvalues, measure effective rank vs grokking time | https://github.com/KindXiaoming/Omnigrok |
| **Tian (2025)** | Compute effective bandwidth from their energy function ℰ, test if t_g ~ f(bandwidth/task_dim) | https://github.com/yuandong-tian/understanding |
| **Merrill (2023)** | Measure effective rank of sparse vs dense circuits during competition | https://github.com/Tsili42/parity-nn |
| **Nanda (2023)** | Add effective rank tracking to their transformer training, correlate with Fourier component count | Code on neelnanda.io |

### Tier 2: Reinterpret Published Results

| Paper | What's reusable |
|-------|----------------|
| **Liu, Kitouni & Tegmark (2022)** | λ₃ eigenvalue scaling with data fraction — reinterpret as bandwidth measure |
| **Yunis (2024)** | Singular value evolution plots — extract effective rank trajectories |
| **Clauw (2024)** | Synergy/redundancy measures — test if synergy/(synergy+redundancy) ≈ closure ratio |
| **Davies (2023)** | Pattern learning speed framework — formalize as closure rate |

### Tier 3: Conceptual Support

| Paper | Connection |
|-------|-----------|
| **Thilak (2022)** | Slingshot = mechanism for reaching closure point |
| **Varma (2023)** | Circuit efficiency ≈ closure ratio |
| **Shwartz-Ziv & Tishby (2017)** | Compression phase ≈ approach to closure |
| **Neural Collapse paper (2025)** | Within-class variance contraction = bandwidth focusing |

---

## Key Gap in Literature

**No existing paper systematically varies network width while tracking:**
1. Grokking onset time
2. Effective rank of activation covariances
3. Singular value spectrum of weight matrices
4. Any bandwidth-like measure

**This is the primary experimental contribution closure theory can make.** The experiment would be:
- Fix task (e.g., modular addition mod 113)
- Vary MLP width: [32, 64, 128, 256, 512, 1024]
- For each width, run multiple seeds
- Log: activation covariance eigenvalues per layer per epoch, weight SVD per epoch, grokking time
- Test: does t_g scale with effective bandwidth (computed from activation covariance) rather than raw parameter count?

### Specific Predictions to Test
1. **Closure ratio prediction:** t_g ∝ 1/(effective_bandwidth - task_dimensionality) near the critical point
2. **Width scaling:** Wider networks should grok faster, but with diminishing returns once effective bandwidth >> task dimensionality
3. **Effective rank at grokking:** At the moment of grokking, the effective rank of layer activations should approximately equal the intrinsic task dimensionality (≈ number of required Fourier modes for modular addition)
4. **Bandwidth-not-parameters:** Two architectures with the same parameter count but different effective bandwidth (e.g., deep-narrow vs shallow-wide) should have different grokking times predicted by bandwidth, not parameter count

---

## Existing Results That May Support/Contradict Closure

### Supporting
- Yunis (2024): Grokking coincides with rank minimization ✅
- Nanda (2023): Fourier components (= bandwidth components) smoothly increase before grokking ✅
- Merrill (2023): Generalizing circuit is sparse (= low effective bandwidth) ✅
- Liu et al. (2022): t ~ 1/λ₃ where λ₃ is a representation quality eigenvalue ✅
- Clauw (2024): Synergy (= cooperative bandwidth) peaks predict grokking ✅

### Potentially Contradicting
- Thilak (2022): If grokking is purely an optimizer artifact (slingshots), then bandwidth-based predictions may not capture the full picture
- Shwartz-Ziv controversy: If information compression doesn't occur with ReLU, closure-as-compression needs refinement

### Neither (need more data)
- Power (2022): d_model=32 never overfitting could mean bandwidth was never sufficient, OR that the architecture was below a capacity threshold
- Davies (2023): Framework is compatible but not specific enough to distinguish closure from other explanations
