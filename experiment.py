#!/usr/bin/env python3
"""
Closure Theory × Grokking Experiment
=====================================

Tests the core prediction of Closure Theory:
    Grokking onset time t_g ∝ B_accessible / İ_internal

where:
    B_accessible = total representational bandwidth (bits available)
    İ_internal   = rate of internal correlation production

Task: modular multiplication (a * b mod 97) learned by a small MLP.
Sweep over bottleneck widths to vary B_accessible while keeping the
task fixed, then check whether t_g scales linearly with bandwidth.

Author: Generated for Closure Theory research
License: MIT
"""

import csv
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

# ──────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────

PRIME = 97                          # modulus for a*b mod p
TRAIN_FRAC = 0.30                   # 30 % train, 70 % test
MAX_STEPS = 50_000
LR = 1e-3
WEIGHT_DECAY = 1.0
LOG_EVERY = 100                     # log closure metrics every N steps
BOTTLENECK_WIDTHS = [32, 64, 128, 256, 512]
EMBED_DIM = 128                     # token embedding dimension
HIDDEN_DIM = 256                    # fixed outer hidden width
BITS_PER_PARAM = 32                 # float32 precision
ETA = 1e-6                          # regulariser for log-det covariance
SEED = 42
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Grokking detection parameters
GROK_JUMP = 0.15                    # minimum accuracy jump
GROK_WINDOW = 500 // LOG_EVERY      # window in logged steps (5 entries)
GROK_SUSTAIN_TOL = 0.03             # sustain tolerance

OUT_DIR = Path(__file__).parent / "results"

# ──────────────────────────────────────────────────────────────────────
# Dataset: modular multiplication
# ──────────────────────────────────────────────────────────────────────

def make_dataset(prime: int = PRIME, train_frac: float = TRAIN_FRAC, seed: int = SEED):
    """Create (a, b, a*b mod p) dataset with a fixed train/test split.

    Returns tensors on DEVICE:
        train_a, train_b, train_y, test_a, test_b, test_y
    """
    rng = np.random.RandomState(seed)
    pairs = [(a, b) for a in range(prime) for b in range(prime)]  # p^2 pairs
    rng.shuffle(pairs)
    n_train = int(len(pairs) * train_frac)

    def to_tensors(subset):
        a = torch.tensor([p[0] for p in subset], dtype=torch.long, device=DEVICE)
        b = torch.tensor([p[1] for p in subset], dtype=torch.long, device=DEVICE)
        y = (a * b) % prime
        return a, b, y

    return to_tensors(pairs[:n_train]), to_tensors(pairs[n_train:])


# ──────────────────────────────────────────────────────────────────────
# Model: Embedding + MLP with configurable bottleneck
# ──────────────────────────────────────────────────────────────────────

class GrokMLP(nn.Module):
    """
    Architecture:
        Embed(a) ‖ Embed(b)  →  Linear(2·embed, hidden)  →  ReLU
        →  Linear(hidden, bottleneck)  →  ReLU           ← bandwidth knob
        →  Linear(bottleneck, hidden)  →  ReLU
        →  Linear(hidden, prime)

    The bottleneck layer controls B_accessible: narrower bottleneck =
    less representational bandwidth = (prediction) later grokking.
    """

    def __init__(self, prime: int, embed_dim: int, hidden_dim: int, bottleneck_dim: int):
        super().__init__()
        self.embed_a = nn.Embedding(prime, embed_dim)
        self.embed_b = nn.Embedding(prime, embed_dim)

        self.layers = nn.Sequential(
            nn.Linear(2 * embed_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, bottleneck_dim),   # ← bottleneck
            nn.ReLU(),
            nn.Linear(bottleneck_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, prime),
        )

    def forward(self, a: torch.Tensor, b: torch.Tensor):
        x = torch.cat([self.embed_a(a), self.embed_b(b)], dim=-1)
        return self.layers(x)

    def get_layer_activations(self, a: torch.Tensor, b: torch.Tensor):
        """Return list of post-activation tensors for each hidden layer."""
        x = torch.cat([self.embed_a(a), self.embed_b(b)], dim=-1)
        activations = []
        for module in self.layers:
            x = module(x)
            if isinstance(module, nn.ReLU):
                activations.append(x.detach())
        return activations

    def layer_widths(self) -> list[int]:
        """Return widths of each linear layer (input side)."""
        return [m.out_features for m in self.layers if isinstance(m, nn.Linear)]


# ──────────────────────────────────────────────────────────────────────
# Closure metrics
# ──────────────────────────────────────────────────────────────────────

def compute_I_internal(activations: list[torch.Tensor], eta: float = ETA) -> float:
    """
    I_internal(t) = Σ_ℓ  log det(Σ_ℓ + η·I)

    where Σ_ℓ is the empirical covariance of activations at layer ℓ
    over the probe batch (training data).

    This measures how much internal correlation structure the network
    has built.  Higher values → richer internal representations.
    """
    total = 0.0
    for act in activations:
        # act shape: (N, d_ℓ)
        act_centered = act - act.mean(dim=0, keepdim=True)
        cov = (act_centered.T @ act_centered) / (act.shape[0] - 1)
        cov += eta * torch.eye(cov.shape[0], device=cov.device)
        # log-det via Cholesky for numerical stability
        try:
            L = torch.linalg.cholesky(cov)
            logdet = 2.0 * L.diagonal().log().sum().item()
        except RuntimeError:
            # Fallback: eigenvalue decomposition
            eigvals = torch.linalg.eigvalsh(cov)
            logdet = eigvals.clamp(min=eta).log().sum().item()
        total += logdet
    return total


def compute_B_accessible(layer_widths: list[int], bits: int = BITS_PER_PARAM) -> float:
    """
    B_accessible = Σ_ℓ  d_ℓ · b_ℓ

    where d_ℓ is layer width and b_ℓ is precision in bits.
    This is the total representational bandwidth — the maximum number
    of bits the network *could* use to encode information.
    """
    return float(sum(w * bits for w in layer_widths))


def compute_effective_rank(activations: list[torch.Tensor], threshold: float = 0.01) -> list[int]:
    """
    Effective rank per layer: count of singular values > 1% of maximum.

    Tracks the intrinsic dimensionality of each layer's representation.
    During grokking, effective rank often drops then recovers as the
    network finds a compact generalising solution.
    """
    ranks = []
    for act in activations:
        act_centered = act - act.mean(dim=0, keepdim=True)
        sv = torch.linalg.svdvals(act_centered)
        ranks.append(int((sv > threshold * sv[0]).sum().item()))
    return ranks


# ──────────────────────────────────────────────────────────────────────
# Grokking detection
# ──────────────────────────────────────────────────────────────────────

def detect_grokking(test_accs: list[float], steps: list[int]) -> Optional[int]:
    """
    Detect grokking onset: find the first step where test accuracy
    jumps by >= GROK_JUMP over a GROK_WINDOW-sized window and stays
    within GROK_SUSTAIN_TOL of the peak for at least GROK_WINDOW more steps.

    Returns the step index of grokking onset, or None if not detected.
    """
    if len(test_accs) < 2 * GROK_WINDOW:
        return None

    for i in range(GROK_WINDOW, len(test_accs) - GROK_WINDOW):
        jump = test_accs[i] - test_accs[i - GROK_WINDOW]
        if jump >= GROK_JUMP:
            # Check sustain: accuracy stays within tolerance of test_accs[i]
            peak = test_accs[i]
            sustained = all(
                abs(test_accs[j] - peak) < GROK_SUSTAIN_TOL
                for j in range(i, min(i + GROK_WINDOW, len(test_accs)))
            )
            if sustained:
                return steps[i]
    return None


# ──────────────────────────────────────────────────────────────────────
# Single run
# ──────────────────────────────────────────────────────────────────────

@dataclass
class RunResult:
    bottleneck: int
    B_accessible: float
    t_grok: Optional[int]
    C_at_grok: Optional[float]
    steps: list[int] = field(default_factory=list)
    train_accs: list[float] = field(default_factory=list)
    test_accs: list[float] = field(default_factory=list)
    I_internals: list[float] = field(default_factory=list)
    closures: list[float] = field(default_factory=list)
    I_dots: list[float] = field(default_factory=list)
    eff_ranks: list[list[int]] = field(default_factory=list)


def run_experiment(bottleneck: int, seed: int = SEED) -> RunResult:
    """Train one model and collect all metrics."""
    torch.manual_seed(seed)
    np.random.seed(seed)

    (train_a, train_b, train_y), (test_a, test_b, test_y) = make_dataset()
    model = GrokMLP(PRIME, EMBED_DIM, HIDDEN_DIM, bottleneck).to(DEVICE)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

    widths = model.layer_widths()
    B = compute_B_accessible(widths)
    result = RunResult(bottleneck=bottleneck, B_accessible=B, t_grok=None, C_at_grok=None)

    prev_I = None
    print(f"\n{'='*60}")
    print(f"  Bottleneck width = {bottleneck}  |  B_accessible = {B:.0f}")
    print(f"{'='*60}")

    for step in range(1, MAX_STEPS + 1):
        # ── Training step ──
        model.train()
        logits = model(train_a, train_b)
        loss = F.cross_entropy(logits, train_y)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if step % LOG_EVERY != 0:
            continue

        # ── Evaluation ──
        model.eval()
        with torch.no_grad():
            train_pred = model(train_a, train_b).argmax(dim=-1)
            train_acc = (train_pred == train_y).float().mean().item()

            test_pred = model(test_a, test_b).argmax(dim=-1)
            test_acc = (test_pred == test_y).float().mean().item()

        # ── Closure metrics ──
        with torch.no_grad():
            acts = model.get_layer_activations(train_a, train_b)

        I = compute_I_internal(acts)
        C = I / B
        I_dot = (I - prev_I) / LOG_EVERY if prev_I is not None else 0.0
        prev_I = I
        eff_rank = compute_effective_rank(acts)

        # ── Store ──
        result.steps.append(step)
        result.train_accs.append(train_acc)
        result.test_accs.append(test_acc)
        result.I_internals.append(I)
        result.closures.append(C)
        result.I_dots.append(I_dot)
        result.eff_ranks.append(eff_rank)

        if step % (LOG_EVERY * 50) == 0:
            print(f"  step {step:>6d}  |  train {train_acc:.3f}  test {test_acc:.3f}  "
                  f"|  I={I:.1f}  C={C:.4f}  İ={I_dot:.4f}  rank={eff_rank}")

    # ── Detect grokking ──
    result.t_grok = detect_grokking(result.test_accs, result.steps)
    if result.t_grok is not None:
        idx = result.steps.index(result.t_grok)
        result.C_at_grok = result.closures[idx]
        print(f"  ✓ Grokking detected at step {result.t_grok}  (C = {result.C_at_grok:.4f})")
    else:
        print(f"  ✗ No grokking detected within {MAX_STEPS} steps")

    return result


# ──────────────────────────────────────────────────────────────────────
# CSV output
# ──────────────────────────────────────────────────────────────────────

def save_csv(result: RunResult, out_dir: Path):
    """Save per-step metrics to CSV."""
    path = out_dir / f"metrics_bottleneck_{result.bottleneck}.csv"
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["step", "train_acc", "test_acc", "I_internal", "B_accessible",
                  "closure_ratio", "I_dot"] + [f"eff_rank_L{i}" for i in range(len(result.eff_ranks[0]))]
        writer.writerow(header)
        for i, step in enumerate(result.steps):
            row = [step, f"{result.train_accs[i]:.6f}", f"{result.test_accs[i]:.6f}",
                   f"{result.I_internals[i]:.6f}", f"{result.B_accessible:.0f}",
                   f"{result.closures[i]:.6f}", f"{result.I_dots[i]:.6f}"]
            row += [str(r) for r in result.eff_ranks[i]]
            writer.writerow(row)
    print(f"  Saved {path}")


# ──────────────────────────────────────────────────────────────────────
# Plotting
# ──────────────────────────────────────────────────────────────────────

def generate_plots(results: list[RunResult], out_dir: Path):
    """Generate the four key diagnostic plots."""
    fig_dir = out_dir
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(results)))

    # ── Plot A: Test accuracy vs step (all widths overlaid) ──
    fig, ax = plt.subplots(figsize=(10, 6))
    for res, c in zip(results, colors):
        ax.plot(res.steps, res.test_accs, color=c, label=f"w={res.bottleneck}", alpha=0.85)
        if res.t_grok is not None:
            ax.axvline(res.t_grok, color=c, linestyle="--", alpha=0.4)
    ax.set_xlabel("Training Step")
    ax.set_ylabel("Test Accuracy")
    ax.set_title("Grokking: Test Accuracy vs Training Step")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "plot_a_test_accuracy.png", dpi=150)
    plt.close(fig)
    print(f"  Saved plot_a_test_accuracy.png")

    # ── Plot B: Closure ratio C(t) vs step ──
    fig, ax = plt.subplots(figsize=(10, 6))
    for res, c in zip(results, colors):
        ax.plot(res.steps, res.closures, color=c, label=f"w={res.bottleneck}", alpha=0.85)
        if res.t_grok is not None:
            idx = res.steps.index(res.t_grok)
            ax.plot(res.t_grok, res.closures[idx], "o", color=c, markersize=10,
                    markeredgecolor="black", zorder=5)
    ax.set_xlabel("Training Step")
    ax.set_ylabel("C(t) = I_internal / B_accessible")
    ax.set_title("Closure Ratio Over Training (● = grokking onset)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "plot_b_closure_ratio.png", dpi=150)
    plt.close(fig)
    print(f"  Saved plot_b_closure_ratio.png")

    # ── Plot C: t_g vs B_accessible (key prediction) ──
    grokked = [r for r in results if r.t_grok is not None]
    fig, ax = plt.subplots(figsize=(8, 6))
    if len(grokked) >= 2:
        Bs = [r.B_accessible for r in grokked]
        tgs = [r.t_grok for r in grokked]
        ax.scatter(Bs, tgs, s=100, zorder=5, edgecolors="black")
        for r in grokked:
            ax.annotate(f"w={r.bottleneck}", (r.B_accessible, r.t_grok),
                       textcoords="offset points", xytext=(8, 4), fontsize=9)
        # Linear fit
        coeffs = np.polyfit(Bs, tgs, 1)
        fit_x = np.linspace(min(Bs) * 0.9, max(Bs) * 1.1, 100)
        ax.plot(fit_x, np.polyval(coeffs, fit_x), "--", color="red", alpha=0.7,
                label=f"Linear fit: t_g = {coeffs[0]:.3f}·B + {coeffs[1]:.0f}")
        # R² calculation
        predicted = np.polyval(coeffs, Bs)
        ss_res = np.sum((np.array(tgs) - predicted) ** 2)
        ss_tot = np.sum((np.array(tgs) - np.mean(tgs)) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        ax.set_title(f"KEY TEST: t_g vs B_accessible  (R² = {r_squared:.3f})")
        ax.legend()
    else:
        ax.set_title("KEY TEST: t_g vs B_accessible  (insufficient grokking events)")
    ax.set_xlabel("B_accessible (bits)")
    ax.set_ylabel("Grokking onset step t_g")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "plot_c_tg_vs_bandwidth.png", dpi=150)
    plt.close(fig)
    print(f"  Saved plot_c_tg_vs_bandwidth.png")

    # ── Plot D: C at grokking onset vs bottleneck width ──
    fig, ax = plt.subplots(figsize=(8, 6))
    if len(grokked) >= 2:
        widths = [r.bottleneck for r in grokked]
        Cs = [r.C_at_grok for r in grokked]
        ax.scatter(widths, Cs, s=100, zorder=5, edgecolors="black")
        ax.axhline(np.mean(Cs), color="red", linestyle="--", alpha=0.7,
                   label=f"Mean C = {np.mean(Cs):.4f} ± {np.std(Cs):.4f}")
        for r in grokked:
            ax.annotate(f"w={r.bottleneck}", (r.bottleneck, r.C_at_grok),
                       textcoords="offset points", xytext=(8, 4), fontsize=9)
        cv = np.std(Cs) / np.mean(Cs) if np.mean(Cs) != 0 else float("inf")
        ax.set_title(f"Closure Ratio at Grokking Onset  (CV = {cv:.3f})")
        ax.legend()
    else:
        ax.set_title("Closure Ratio at Grokking Onset  (insufficient data)")
    ax.set_xlabel("Bottleneck Width")
    ax.set_ylabel("C(t_g)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "plot_d_closure_at_grok.png", dpi=150)
    plt.close(fig)
    print(f"  Saved plot_d_closure_at_grok.png")


# ──────────────────────────────────────────────────────────────────────
# Summary table
# ──────────────────────────────────────────────────────────────────────

def print_summary(results: list[RunResult]):
    """Print a summary table of all runs."""
    print(f"\n{'='*80}")
    print(f"  SUMMARY TABLE")
    print(f"{'='*80}")
    print(f"  {'Width':>6s}  {'B_access':>10s}  {'t_grok':>8s}  {'C(t_g)':>10s}  "
          f"{'Final Train':>12s}  {'Final Test':>11s}")
    print(f"  {'-'*6}  {'-'*10}  {'-'*8}  {'-'*10}  {'-'*12}  {'-'*11}")
    for r in results:
        tg = str(r.t_grok) if r.t_grok else "—"
        cg = f"{r.C_at_grok:.4f}" if r.C_at_grok else "—"
        print(f"  {r.bottleneck:>6d}  {r.B_accessible:>10.0f}  {tg:>8s}  {cg:>10s}  "
              f"{r.train_accs[-1]:>12.4f}  {r.test_accs[-1]:>11.4f}")

    grokked = [r for r in results if r.t_grok is not None]
    print(f"\n  Grokked: {len(grokked)}/{len(results)} runs")

    if len(grokked) >= 2:
        Bs = np.array([r.B_accessible for r in grokked])
        tgs = np.array([r.t_grok for r in grokked])
        corr = np.corrcoef(Bs, tgs)[0, 1]
        print(f"  Pearson correlation(B_accessible, t_g) = {corr:.4f}")

        Cs = np.array([r.C_at_grok for r in grokked])
        cv = np.std(Cs) / np.mean(Cs) if np.mean(Cs) != 0 else float("inf")
        print(f"  C at grokking: mean={np.mean(Cs):.4f}, std={np.std(Cs):.4f}, CV={cv:.3f}")

        print(f"\n  VERDICT:")
        if corr > 0.85:
            print(f"  ✓ Strong linear relationship (r={corr:.3f}): SUPPORTS Closure Theory")
        elif corr > 0.5:
            print(f"  ~ Moderate relationship (r={corr:.3f}): INCONCLUSIVE")
        else:
            print(f"  ✗ Weak/no relationship (r={corr:.3f}): CHALLENGES Closure Theory")

        if cv < 0.15:
            print(f"  ✓ C(t_g) approximately constant (CV={cv:.3f}): SUPPORTS universal threshold")
        elif cv < 0.30:
            print(f"  ~ C(t_g) moderately variable (CV={cv:.3f}): INCONCLUSIVE")
        else:
            print(f"  ✗ C(t_g) highly variable (CV={cv:.3f}): CHALLENGES universal threshold")

    print(f"{'='*80}\n")


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────
# Phase 2: Intervention Regimes
# ──────────────────────────────────────────────────────────────────────

# We train a single model, checkpoint it at multiple points before and
# after grokking, then apply perturbations at each checkpoint and
# measure whether global behavior changes.

INTERVENTION_WIDTH = 128            # use mid-range bottleneck
INTERVENTION_CHECKPOINTS = 20       # number of checkpoints to save
PERTURBATION_NORMS = [0.01, 0.05, 0.1, 0.5, 1.0]  # L2 norms for weight perturbation
N_PERTURBATION_TRIALS = 5          # random trials per norm level
REPRESENTATION_INTERVENTIONS = ["rotation", "projection", "sign_flip"]


@dataclass
class InterventionResult:
    """Result of a single intervention at a checkpoint."""
    step: int
    phase: str                      # "pre_closure" or "post_closure"
    intervention_type: str          # "weight_perturb" or "repr_rotation" etc.
    param: float                    # L2 norm or intervention strength
    test_acc_before: float
    test_acc_after: float
    train_acc_before: float
    train_acc_after: float
    delta_test: float
    delta_train: float
    # Representation-level metrics
    repr_similarity: Optional[float] = None  # cosine sim of activations before/after
    effective_rank_before: Optional[list] = None
    effective_rank_after: Optional[list] = None


def get_test_metrics(model, test_a, test_b, test_y, train_a, train_b, train_y):
    """Get accuracy and activations."""
    model.eval()
    with torch.no_grad():
        test_pred = model(test_a, test_b).argmax(dim=-1)
        test_acc = (test_pred == test_y).float().mean().item()
        train_pred = model(train_a, train_b).argmax(dim=-1)
        train_acc = (train_pred == train_y).float().mean().item()
        acts = model.get_layer_activations(train_a, train_b)
    return test_acc, train_acc, acts


def perturb_weights(model, l2_norm: float, seed: int = 0):
    """Add random perturbation to ALL weights with bounded L2 norm."""
    rng = torch.Generator(device=DEVICE).manual_seed(seed)
    # Collect all parameter perturbations
    perturbations = []
    total_numel = sum(p.numel() for p in model.parameters())
    for p in model.parameters():
        delta = torch.randn(p.shape, device=DEVICE, generator=rng)
        perturbations.append(delta)
    # Normalize to desired L2 norm
    total_norm = torch.sqrt(sum((d**2).sum() for d in perturbations))
    scale = l2_norm / (total_norm + 1e-12)
    for p, delta in zip(model.parameters(), perturbations):
        p.data.add_(delta * scale)


def intervene_representation(model, method: str, layer_idx: int = 2,
                              train_a=None, train_b=None):
    """
    Apply intervention at the REPRESENTATION level (post-bottleneck).

    layer_idx=2 targets the bottleneck→hidden weight matrix (layers[4] in Sequential,
    which is Linear(bottleneck, hidden)).

    Methods:
    - rotation: Apply random orthogonal rotation to bottleneck output space
    - projection: Project out the top principal component of activations
    - sign_flip: Flip sign of the dominant feature direction
    """
    # Find the bottleneck output layer (layer after bottleneck ReLU)
    linear_layers = [m for m in model.layers if isinstance(m, nn.Linear)]
    # linear_layers: [in→hidden, hidden→bottleneck, bottleneck→hidden, hidden→out]
    # We target index 2: bottleneck→hidden (reads from bottleneck representation)
    target_layer = linear_layers[layer_idx]
    W = target_layer.weight.data  # shape: (hidden, bottleneck)

    if method == "rotation":
        # Random orthogonal rotation of the input space (bottleneck dims)
        d = W.shape[1]
        # Generate random orthogonal matrix via QR decomposition
        random_mat = torch.randn(d, d, device=DEVICE)
        Q, _ = torch.linalg.qr(random_mat)
        # Apply: new_W = W @ Q (rotate the input representation)
        target_layer.weight.data = W @ Q

    elif method == "projection":
        # Find dominant direction of bottleneck activations, project it out
        if train_a is not None:
            with torch.no_grad():
                acts = model.get_layer_activations(train_a, train_b)
                # Bottleneck is activation index 1 (after 2nd ReLU)
                bottleneck_acts = acts[1]  # shape: (N, bottleneck_dim)
                # SVD to find top direction
                U, S, Vh = torch.linalg.svd(bottleneck_acts - bottleneck_acts.mean(0), full_matrices=False)
                top_dir = Vh[0]  # shape: (bottleneck_dim,)
                # Project out: W' = W @ (I - v v^T)
                proj = torch.eye(W.shape[1], device=DEVICE) - torch.outer(top_dir, top_dir)
                target_layer.weight.data = W @ proj

    elif method == "sign_flip":
        # Flip sign of the dominant feature direction
        if train_a is not None:
            with torch.no_grad():
                acts = model.get_layer_activations(train_a, train_b)
                bottleneck_acts = acts[1]
                U, S, Vh = torch.linalg.svd(bottleneck_acts - bottleneck_acts.mean(0), full_matrices=False)
                top_dir = Vh[0]
                # Reflect: W' = W @ (I - 2 v v^T)
                reflect = torch.eye(W.shape[1], device=DEVICE) - 2 * torch.outer(top_dir, top_dir)
                target_layer.weight.data = W @ reflect


def cosine_similarity_activations(acts_a: list[torch.Tensor], acts_b: list[torch.Tensor]) -> float:
    """Mean cosine similarity between corresponding layer activations."""
    sims = []
    for a, b in zip(acts_a, acts_b):
        # Flatten and compute cosine sim
        a_flat = a.reshape(-1)
        b_flat = b.reshape(-1)
        sim = F.cosine_similarity(a_flat.unsqueeze(0), b_flat.unsqueeze(0)).item()
        sims.append(sim)
    return float(np.mean(sims))


def run_intervention_experiment(seed: int = SEED):
    """
    Phase 2: Train one model, checkpoint before/after grokking,
    then test both weight-level and representation-level interventions.
    """
    print(f"\n{'#'*80}")
    print(f"  PHASE 2: INTERVENTION REGIMES")
    print(f"{'#'*80}")

    torch.manual_seed(seed)
    np.random.seed(seed)

    (train_a, train_b, train_y), (test_a, test_b, test_y) = make_dataset()
    model = GrokMLP(PRIME, EMBED_DIM, HIDDEN_DIM, INTERVENTION_WIDTH).to(DEVICE)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

    # ── Phase 2a: Train and save checkpoints ──
    checkpoints = {}   # step → state_dict
    checkpoint_interval = MAX_STEPS // INTERVENTION_CHECKPOINTS
    test_accs_log = []
    steps_log = []

    print(f"\n  Training with checkpoints (width={INTERVENTION_WIDTH})...")
    for step in range(1, MAX_STEPS + 1):
        model.train()
        logits = model(train_a, train_b)
        loss = F.cross_entropy(logits, train_y)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if step % LOG_EVERY == 0:
            model.eval()
            with torch.no_grad():
                test_acc = (model(test_a, test_b).argmax(-1) == test_y).float().mean().item()
            test_accs_log.append(test_acc)
            steps_log.append(step)

        if step % checkpoint_interval == 0:
            checkpoints[step] = {k: v.clone() for k, v in model.state_dict().items()}
            test_acc_now = test_accs_log[-1] if test_accs_log else 0
            print(f"    Checkpoint at step {step:>6d}  test_acc={test_acc_now:.3f}")

    # Detect grokking
    t_grok = detect_grokking(test_accs_log, steps_log)
    if t_grok is None:
        print("  ⚠ No grokking detected — using midpoint as boundary")
        t_grok = MAX_STEPS // 2

    print(f"\n  Grokking boundary: step {t_grok}")
    print(f"  Pre-closure checkpoints: {[s for s in sorted(checkpoints) if s < t_grok]}")
    print(f"  Post-closure checkpoints: {[s for s in sorted(checkpoints) if s >= t_grok]}")

    # ── Phase 2b: Weight-level perturbations (pre vs post closure) ──
    all_results = []

    print(f"\n  {'─'*60}")
    print(f"  REGIME 1: Weight-level perturbations")
    print(f"  {'─'*60}")

    for ckpt_step in sorted(checkpoints.keys()):
        phase = "pre_closure" if ckpt_step < t_grok else "post_closure"

        for norm in PERTURBATION_NORMS:
            deltas_test = []
            deltas_train = []

            for trial in range(N_PERTURBATION_TRIALS):
                # Restore checkpoint
                model.load_state_dict({k: v.clone() for k, v in checkpoints[ckpt_step].items()})
                test_before, train_before, acts_before = get_test_metrics(
                    model, test_a, test_b, test_y, train_a, train_b, train_y)

                # Perturb
                perturb_weights(model, norm, seed=trial)
                test_after, train_after, acts_after = get_test_metrics(
                    model, test_a, test_b, test_y, train_a, train_b, train_y)

                delta_test = abs(test_after - test_before)
                delta_train = abs(train_after - train_before)
                deltas_test.append(delta_test)
                deltas_train.append(delta_train)

                repr_sim = cosine_similarity_activations(acts_before, acts_after)
                eff_rank_before = compute_effective_rank(acts_before)
                eff_rank_after = compute_effective_rank(acts_after)

                all_results.append(InterventionResult(
                    step=ckpt_step, phase=phase,
                    intervention_type="weight_perturb", param=norm,
                    test_acc_before=test_before, test_acc_after=test_after,
                    train_acc_before=train_before, train_acc_after=train_after,
                    delta_test=delta_test, delta_train=delta_train,
                    repr_similarity=repr_sim,
                    effective_rank_before=eff_rank_before,
                    effective_rank_after=eff_rank_after,
                ))

            mean_dt = np.mean(deltas_test)
            mean_dtr = np.mean(deltas_train)
            if norm in [0.01, 0.1, 1.0]:  # print subset to avoid spam
                print(f"    step={ckpt_step:>6d} [{phase:>13s}]  L2={norm:.2f}  "
                      f"Δtest={mean_dt:.4f}  Δtrain={mean_dtr:.4f}")

    # ── Phase 2c: Representation-level interventions (post-closure only) ──
    print(f"\n  {'─'*60}")
    print(f"  REGIME 2: Representation-level interventions (post-closure)")
    print(f"  {'─'*60}")

    post_checkpoints = {s: sd for s, sd in checkpoints.items() if s >= t_grok}
    if not post_checkpoints:
        # Use the last checkpoint
        last_step = max(checkpoints.keys())
        post_checkpoints = {last_step: checkpoints[last_step]}

    for ckpt_step in sorted(post_checkpoints.keys()):
        for method in REPRESENTATION_INTERVENTIONS:
            # Restore
            model.load_state_dict({k: v.clone() for k, v in checkpoints[ckpt_step].items()})
            test_before, train_before, acts_before = get_test_metrics(
                model, test_a, test_b, test_y, train_a, train_b, train_y)

            # Intervene at representation level
            intervene_representation(model, method, layer_idx=2,
                                      train_a=train_a, train_b=train_b)
            test_after, train_after, acts_after = get_test_metrics(
                model, test_a, test_b, test_y, train_a, train_b, train_y)

            repr_sim = cosine_similarity_activations(acts_before, acts_after)
            eff_rank_before = compute_effective_rank(acts_before)
            eff_rank_after = compute_effective_rank(acts_after)

            delta_test = test_after - test_before  # signed — we want to see direction
            delta_train = train_after - train_before

            all_results.append(InterventionResult(
                step=ckpt_step, phase="post_closure",
                intervention_type=f"repr_{method}", param=1.0,
                test_acc_before=test_before, test_acc_after=test_after,
                train_acc_before=train_before, train_acc_after=train_after,
                delta_test=delta_test, delta_train=delta_train,
                repr_similarity=repr_sim,
                effective_rank_before=eff_rank_before,
                effective_rank_after=eff_rank_after,
            ))

            print(f"    step={ckpt_step:>6d}  {method:>12s}  "
                  f"test: {test_before:.3f}→{test_after:.3f} (Δ={delta_test:+.4f})  "
                  f"repr_sim={repr_sim:.4f}")

    return all_results, t_grok, test_accs_log, steps_log, checkpoints


# ──────────────────────────────────────────────────────────────────────
# Phase 2 Plotting
# ──────────────────────────────────────────────────────────────────────

def generate_intervention_plots(results: list[InterventionResult], t_grok: int,
                                  test_accs: list, steps: list, out_dir: Path):
    """Generate intervention diagnostic plots."""

    # ── Plot E: Weight perturbation sensitivity (pre vs post closure) ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, phase, title in zip(axes, ["pre_closure", "post_closure"],
                                 ["Pre-Closure", "Post-Closure"]):
        phase_results = [r for r in results
                         if r.phase == phase and r.intervention_type == "weight_perturb"]
        if not phase_results:
            ax.set_title(f"{title}: No data")
            continue

        # Group by step
        step_groups = {}
        for r in phase_results:
            step_groups.setdefault(r.step, []).append(r)

        colors_e = plt.cm.coolwarm(np.linspace(0, 1, len(step_groups)))
        for (step, group), c in zip(sorted(step_groups.items()), colors_e):
            # Group by norm
            norm_to_deltas = {}
            for r in group:
                norm_to_deltas.setdefault(r.param, []).append(abs(r.delta_test))
            norms = sorted(norm_to_deltas.keys())
            means = [np.mean(norm_to_deltas[n]) for n in norms]
            ax.plot(norms, means, 'o-', color=c, label=f"step {step}", alpha=0.8)

        ax.set_xlabel("Perturbation L2 Norm")
        ax.set_ylabel("|Δ Test Accuracy|")
        ax.set_title(f"{title}: Weight Perturbation Sensitivity")
        ax.legend(fontsize=7, ncol=2)
        ax.grid(True, alpha=0.3)
        ax.set_xscale("log")

    fig.tight_layout()
    fig.savefig(out_dir / "plot_e_weight_perturbation.png", dpi=150)
    plt.close(fig)
    print(f"  Saved plot_e_weight_perturbation.png")

    # ── Plot F: Representation intervention impact ──
    repr_results = [r for r in results if r.intervention_type.startswith("repr_")]
    if repr_results:
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

        # F1: Delta test accuracy by intervention type
        ax = axes[0]
        for method in REPRESENTATION_INTERVENTIONS:
            method_results = [r for r in repr_results if r.intervention_type == f"repr_{method}"]
            if method_results:
                steps_m = [r.step for r in method_results]
                deltas = [r.delta_test for r in method_results]
                ax.plot(steps_m, deltas, 'o-', label=method, markersize=8)
        ax.axhline(0, color='black', linestyle='-', alpha=0.3)
        ax.axvline(t_grok, color='red', linestyle='--', alpha=0.5, label=f"grok @ {t_grok}")
        ax.set_xlabel("Checkpoint Step")
        ax.set_ylabel("Δ Test Accuracy (signed)")
        ax.set_title("Representation Interventions: Global Effect")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # F2: Representation similarity after intervention
        ax = axes[1]
        for method in REPRESENTATION_INTERVENTIONS:
            method_results = [r for r in repr_results if r.intervention_type == f"repr_{method}"]
            if method_results:
                steps_m = [r.step for r in method_results]
                sims = [r.repr_similarity for r in method_results]
                ax.plot(steps_m, sims, 's-', label=method, markersize=8)
        ax.axvline(t_grok, color='red', linestyle='--', alpha=0.5, label=f"grok @ {t_grok}")
        ax.set_xlabel("Checkpoint Step")
        ax.set_ylabel("Cosine Similarity (before/after)")
        ax.set_title("Representation Similarity After Intervention")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # F3: Effective rank change
        ax = axes[2]
        for method in REPRESENTATION_INTERVENTIONS:
            method_results = [r for r in repr_results if r.intervention_type == f"repr_{method}"]
            if method_results:
                steps_m = [r.step for r in method_results]
                rank_deltas = [
                    np.mean(r.effective_rank_after) - np.mean(r.effective_rank_before)
                    for r in method_results
                    if r.effective_rank_before and r.effective_rank_after
                ]
                if rank_deltas:
                    ax.plot(steps_m[:len(rank_deltas)], rank_deltas, 'D-', label=method, markersize=8)
        ax.axhline(0, color='black', linestyle='-', alpha=0.3)
        ax.axvline(t_grok, color='red', linestyle='--', alpha=0.5, label=f"grok @ {t_grok}")
        ax.set_xlabel("Checkpoint Step")
        ax.set_ylabel("Δ Mean Effective Rank")
        ax.set_title("Effective Rank Reorganization")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(out_dir / "plot_f_representation_interventions.png", dpi=150)
        plt.close(fig)
        print(f"  Saved plot_f_representation_interventions.png")

    # ── Plot G: Pre vs Post closure sensitivity comparison ──
    fig, ax = plt.subplots(figsize=(10, 6))
    weight_results = [r for r in results if r.intervention_type == "weight_perturb"]

    for phase, color, marker in [("pre_closure", "blue", "o"), ("post_closure", "red", "s")]:
        phase_r = [r for r in weight_results if r.phase == phase]
        if not phase_r:
            continue
        # Aggregate: mean delta per norm across all checkpoints
        norm_to_deltas = {}
        for r in phase_r:
            norm_to_deltas.setdefault(r.param, []).append(abs(r.delta_test))
        norms = sorted(norm_to_deltas.keys())
        means = [np.mean(norm_to_deltas[n]) for n in norms]
        stds = [np.std(norm_to_deltas[n]) for n in norms]
        ax.errorbar(norms, means, yerr=stds, fmt=f'{marker}-', color=color,
                     label=f"{phase.replace('_', ' ').title()}", capsize=4, markersize=8)

    ax.set_xlabel("Perturbation L2 Norm")
    ax.set_ylabel("|Δ Test Accuracy| (mean ± std)")
    ax.set_title("KEY TEST: Pre vs Post Closure Weight Perturbation Sensitivity")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xscale("log")
    fig.tight_layout()
    fig.savefig(out_dir / "plot_g_pre_vs_post_sensitivity.png", dpi=150)
    plt.close(fig)
    print(f"  Saved plot_g_pre_vs_post_sensitivity.png")


def print_intervention_summary(results: list[InterventionResult], t_grok: int):
    """Print Phase 2 summary with verdicts."""
    print(f"\n{'#'*80}")
    print(f"  PHASE 2 SUMMARY: INTERVENTION REGIMES")
    print(f"{'#'*80}")

    weight_results = [r for r in results if r.intervention_type == "weight_perturb"]
    repr_results = [r for r in results if r.intervention_type.startswith("repr_")]

    # ── Weight perturbation analysis ──
    print(f"\n  REGIME 1: Weight Perturbations (bounded L2)")
    print(f"  {'─'*50}")

    pre_deltas = [abs(r.delta_test) for r in weight_results if r.phase == "pre_closure"]
    post_deltas = [abs(r.delta_test) for r in weight_results if r.phase == "post_closure"]

    if pre_deltas and post_deltas:
        pre_mean = np.mean(pre_deltas)
        post_mean = np.mean(post_deltas)
        print(f"  Pre-closure  mean |Δtest|: {pre_mean:.4f}")
        print(f"  Post-closure mean |Δtest|: {post_mean:.4f}")
        print(f"  Ratio (post/pre): {post_mean/pre_mean:.2f}x")

        # At high perturbation norms specifically
        pre_high = [abs(r.delta_test) for r in weight_results
                     if r.phase == "pre_closure" and r.param >= 0.5]
        post_high = [abs(r.delta_test) for r in weight_results
                      if r.phase == "post_closure" and r.param >= 0.5]
        if pre_high and post_high:
            print(f"\n  At high perturbation (L2≥0.5):")
            print(f"    Pre-closure  mean |Δtest|: {np.mean(pre_high):.4f}")
            print(f"    Post-closure mean |Δtest|: {np.mean(post_high):.4f}")

        print(f"\n  PREDICTION: Pre-closure perturbations should FAIL to induce generalization")
        if pre_mean < 0.1:
            print(f"  ✓ Pre-closure weight perturbations had minimal effect ({pre_mean:.4f})")
        else:
            print(f"  ✗ Pre-closure weight perturbations had significant effect ({pre_mean:.4f})")

    # ── Representation intervention analysis ──
    print(f"\n  REGIME 2: Representation-Level Interventions (post-closure)")
    print(f"  {'─'*50}")

    if repr_results:
        for method in REPRESENTATION_INTERVENTIONS:
            method_r = [r for r in repr_results if r.intervention_type == f"repr_{method}"]
            if method_r:
                deltas = [r.delta_test for r in method_r]
                sims = [r.repr_similarity for r in method_r if r.repr_similarity is not None]
                print(f"\n  {method.upper()}:")
                print(f"    Δtest: {np.mean(deltas):+.4f} (mean), range [{min(deltas):+.4f}, {max(deltas):+.4f}]")
                if sims:
                    print(f"    Repr similarity: {np.mean(sims):.4f}")
                # Check for global reorganization
                if abs(np.mean(deltas)) > 0.1:
                    print(f"    ✓ Global reorganization detected (|Δ|>{0.1})")
                else:
                    print(f"    ~ Modest effect (|Δ|<{0.1})")

        print(f"\n  PREDICTION: Single representation-level intervention should cause GLOBAL reorganization")
        all_repr_deltas = [abs(r.delta_test) for r in repr_results]
        mean_repr_delta = np.mean(all_repr_deltas) if all_repr_deltas else 0
        mean_weight_delta = np.mean(post_deltas) if post_deltas else 0

        if mean_repr_delta > mean_weight_delta * 1.5:
            print(f"  ✓ Representation interventions ({mean_repr_delta:.4f}) >> "
                  f"weight perturbations ({mean_weight_delta:.4f})")
            print(f"    SUPPORTS: Post-closure representations are causally structured")
        elif mean_repr_delta > mean_weight_delta:
            print(f"  ~ Representation interventions ({mean_repr_delta:.4f}) > "
                  f"weight perturbations ({mean_weight_delta:.4f})")
            print(f"    INCONCLUSIVE: Trend in right direction but not definitive")
        else:
            print(f"  ✗ Representation interventions ({mean_repr_delta:.4f}) ≤ "
                  f"weight perturbations ({mean_weight_delta:.4f})")
            print(f"    CHALLENGES: No evidence of special causal structure at representation level")

    print(f"\n{'#'*80}\n")


# ──────────────────────────────────────────────────────────────────────
# Phase 2 CSV
# ──────────────────────────────────────────────────────────────────────

def save_intervention_csv(results: list[InterventionResult], out_dir: Path):
    """Save intervention results to CSV."""
    path = out_dir / "intervention_results.csv"
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["step", "phase", "intervention_type", "param",
                         "test_acc_before", "test_acc_after", "train_acc_before", "train_acc_after",
                         "delta_test", "delta_train", "repr_similarity"])
        for r in results:
            writer.writerow([r.step, r.phase, r.intervention_type, f"{r.param:.4f}",
                             f"{r.test_acc_before:.6f}", f"{r.test_acc_after:.6f}",
                             f"{r.train_acc_before:.6f}", f"{r.train_acc_after:.6f}",
                             f"{r.delta_test:.6f}", f"{r.delta_train:.6f}",
                             f"{r.repr_similarity:.6f}" if r.repr_similarity is not None else ""])
    print(f"  Saved {path}")


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Device: {DEVICE}")
    print(f"Output: {OUT_DIR}")
    print(f"Bottleneck widths: {BOTTLENECK_WIDTHS}")
    print(f"Max steps: {MAX_STEPS}, LR: {LR}, WD: {WEIGHT_DECAY}")
    print(f"Dataset: {PRIME}² = {PRIME**2} pairs, {TRAIN_FRAC*100:.0f}% train")

    t0 = time.time()

    # ── Phase 1: Observation (bandwidth sweep) ──
    print(f"\n{'#'*80}")
    print(f"  PHASE 1: OBSERVATION (Bandwidth Sweep)")
    print(f"{'#'*80}")
    results = []
    for width in BOTTLENECK_WIDTHS:
        result = run_experiment(width)
        save_csv(result, OUT_DIR)
        results.append(result)

    print(f"\nGenerating Phase 1 plots...")
    generate_plots(results, OUT_DIR)
    print_summary(results)

    # ── Phase 2: Intervention regimes ──
    intervention_results, t_grok, test_accs, steps, ckpts = run_intervention_experiment()
    save_intervention_csv(intervention_results, OUT_DIR)
    generate_intervention_plots(intervention_results, t_grok, test_accs, steps, OUT_DIR)
    print_intervention_summary(intervention_results, t_grok)

    print(f"Total time: {time.time() - t0:.1f}s")
    print(f"All results in: {OUT_DIR}/")


if __name__ == "__main__":
    main()
