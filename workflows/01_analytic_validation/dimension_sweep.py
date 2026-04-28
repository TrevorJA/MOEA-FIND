"""Diagnostic: shell-versus-interior coverage of the current MOEA-FIND formulation.

The method critic and the user both raised the concern that the L1 auxiliary
objective might produce only shell coverage of the feasible drought region,
that is, solutions concentrated at the maximum achievable L1 distance from
the anti-ideal rather than filling the interior. If the method only delivers
shell coverage, the central claim of MOEA-FIND fails for any realistic
hydrologic application where the feasible region is a proper subset of the
bounding box.

This diagnostic tests the claim empirically on a k-dimensional analytic
problem with a non-trivial constrained feasible region, using the legacy
objective path that the current analytic and Kirsch experiments actually
invoke (src.objectives.analytic_objectives and src.objectives.drought_objectives
with target=None):

    J_i(X) = X_i            for i = 1, ..., k
    J_{k+1}(X) = ||X - D*||_1     (Manhattan norm from anti-ideal)

Decision space: X in [DV_LO, DV_HI]^k.
Anti-ideal:    D* = (DV_HI, ..., DV_HI).
Feasible set:  a k-ball of radius R centered at the origin, so sum X_i^2 <= R^2.

The k-ball has a non-trivial interior distinct from its shell, the
anti-ideal sits outside it, and its volume ratio inside the bounding box
shrinks rapidly with k. The latter is deliberate: it makes the higher
dimensional runs an informative stress test of whether Borg's exploration
can still find and tile the feasible interior when rejection cost is high.

Reference samples: uniform, Latin hypercube, and Sobol samples inside the
k-ball by rejection, matched to the archive size.

Diagnostic metrics:
    1. Distance-from-anti-ideal distribution (the direct shell-only probe).
    2. Distance-from-ball-boundary histogram (interior mass versus shell mass).
    3. Orthant occupancy. A k-ball contains 2^k signed orthants relative to
       the origin; how many of them does each sampler populate?
    4. Cell occupancy on a coarse regular partition of the bounding box,
       restricted to cells whose centers lie inside the ball.
    5. A 2x3 visual panel for k <= 3; a metrics-only panel for k >= 4.

Outputs to outputs/diag_shell_vs_interior/k{K}/:
    - results.json
    - samples.npz
    - fig_shell_vs_interior_k{K}.pdf (or metrics-only for higher k)

Run:
    python scripts/diag_shell_vs_interior.py --k 3 --nfe 40000 --seed 42

Per the active plan file, this diagnostic is blocking on any further
manuscript or figure redesign work, and the sweep across k informs the
where-does-it-break-down question.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

_stub = MagicMock()
sys.modules.setdefault("synhydro", _stub)
sys.modules.setdefault("synhydro.droughts", _stub.droughts)
sys.modules.setdefault("synhydro.droughts.ssi", _stub.droughts.ssi)

from platypus import EpsNSGAII, Problem, Real  # noqa: E402

from src.analysis import (  # noqa: E402
    generate_lhs_samples,
    generate_sobol_samples,
)
from src.objectives import analytic_objectives  # noqa: E402


# Global problem constants. DV_LO, DV_HI, and FEASIBLE_RADIUS are kept
# fixed across k so that results from different dimensionalities are
# directly comparable; the bounding box side is 6 and the feasible-region
# half-width / radius is 2.5. The anti-ideal sits at the positive corner
# of the bounding box, so it is outside both the ball (radius 2.5) and
# the cube (half-width 2.5).
DV_LO, DV_HI = -3.0, 3.0
FEASIBLE_RADIUS = 2.5
VALID_SHAPES = ("ball", "cube")
OUTPUT_SLUG = "diag_shell_vs_interior"


def _in_feasible(
    x: np.ndarray,
    shape: str = "cube",
    radius: float = FEASIBLE_RADIUS,
) -> np.ndarray:
    """Boolean mask indicating which rows of x lie inside the feasible set.

    Args:
        x: Array of shape (..., k).
        shape: ``"ball"`` for ``sum(X_i^2) <= R^2`` or ``"cube"`` for
            ``max(|X_i|) <= R`` (an L-inf ball centered at origin).
        radius: Ball radius / cube half-width.
    """
    if shape == "ball":
        return (x ** 2).sum(axis=-1) <= radius ** 2
    if shape == "cube":
        return np.max(np.abs(x), axis=-1) <= radius
    raise ValueError(f"unknown feasible shape {shape!r} (choose from {VALID_SHAPES})")


def _reject_into_feasible(
    samples: np.ndarray,
    shape: str = "cube",
    radius: float = FEASIBLE_RADIUS,
) -> np.ndarray:
    """Return only rows of samples that fall inside the feasible set."""
    return samples[_in_feasible(samples, shape, radius)]


def _sample_in_feasible(
    n_target: int,
    method: str,
    seed: int,
    k: int,
    shape: str = "cube",
    radius: float = FEASIBLE_RADIUS,
    oversample: int = 20,
) -> np.ndarray:
    """Draw n_target samples inside the feasible set by rejection.

    Uses method in {"uniform", "lhs", "sobol"}. The rejection rate is
    (bounding-box volume) / (feasible-set volume). For a ball this
    grows quickly with k; for a cube it grows as ``((DV_HI - DV_LO) / (2R))^k``
    (much milder). The oversample factor is bumped automatically if the
    first pass falls short.
    """
    lb = np.full(k, DV_LO)
    ub = np.full(k, DV_HI)
    n_draw = max(n_target * oversample, 1024)

    if method == "uniform":
        rng = np.random.default_rng(seed)
        raw = lb + rng.random((n_draw, k)) * (ub - lb)
    elif method == "lhs":
        raw = generate_lhs_samples(n_draw, k, lb, ub, seed=seed)
    elif method == "sobol":
        raw = generate_sobol_samples(n_draw, k, lb, ub, seed=seed)
    else:
        raise ValueError(f"unknown method {method!r}")

    accepted = _reject_into_feasible(raw, shape, radius)
    if len(accepted) < n_target:
        return _sample_in_feasible(
            n_target, method, seed + 1, k, shape, radius,
            oversample=oversample * 4,
        )
    return accepted[:n_target]


def _anti_ideal(k: int) -> np.ndarray:
    return np.full(k, DV_HI)


def _distance_from_antiideal_l1(points: np.ndarray) -> np.ndarray:
    if len(points) == 0:
        return np.zeros(0)
    return np.abs(points - _anti_ideal(points.shape[1])).sum(axis=1)


def _distance_from_boundary(
    points: np.ndarray,
    shape: str = "cube",
    radius: float = FEASIBLE_RADIUS,
) -> np.ndarray:
    """Positive inside, zero on boundary, negative outside.

    For the ball, uses the Euclidean distance from the surface; for the
    cube, uses the L-infinity distance from the nearest face.
    """
    if shape == "ball":
        return radius - np.linalg.norm(points, axis=1)
    if shape == "cube":
        return radius - np.max(np.abs(points), axis=1)
    raise ValueError(f"unknown feasible shape {shape!r}")


def _orthant_occupancy(points: np.ndarray) -> dict:
    """How many signed-sign orthants (relative to origin) the sample hits.

    There are 2^k orthants. For each point we compute a binary sign
    vector and convert it to an integer. The orthant occupancy is the
    number of distinct integers represented in the sample divided by
    2^k. For any interior-filling sampler in a ball centered at the
    origin, this should approach 1.
    """
    if len(points) == 0:
        return {"occupied": 0, "total": 0, "fraction": 0.0}
    k = points.shape[1]
    signs = (points > 0).astype(int)
    ids = np.zeros(len(points), dtype=int)
    for j in range(k):
        ids += signs[:, j] << j
    occupied = len(np.unique(ids))
    total = 1 << k
    return {
        "occupied": int(occupied),
        "total": int(total),
        "fraction": float(occupied) / total,
    }


def _grid_occupancy(
    points: np.ndarray,
    k: int,
    n_bins: int,
    shape: str = "cube",
    radius: float = FEASIBLE_RADIUS,
) -> dict:
    """Feasible-cell occupancy on a coarse regular partition.

    The full grid has n_bins^k cells, of which only those whose center
    lies inside the feasible set are counted as feasible. Occupancy is
    the number of distinct feasible cells that contain at least one
    sample point, divided by the number of feasible cells. For a cube
    feasibility every cell is feasible (the grid exactly tiles the
    cube). For a ball a fraction of cells are outside the ball and are
    excluded.

    For k >= 6 with n_bins = 5 this is 15,625 cells; with n_bins = 4
    it is 4,096. The function caps n_bins^k at 10^5 to bound memory.
    """
    max_cells = 100_000
    while n_bins ** k > max_cells and n_bins > 2:
        n_bins -= 1

    edges = np.linspace(-radius, radius, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    mesh = np.meshgrid(*([centers] * k), indexing="ij")
    cells = np.stack([m.ravel() for m in mesh], axis=1)
    feasible_mask = _in_feasible(cells, shape, radius)
    n_feasible = int(feasible_mask.sum())

    in_feasible_points = (
        points[_in_feasible(points, shape, radius)] if len(points) else points
    )
    if len(in_feasible_points) == 0:
        return {"occupied": 0, "feasible": n_feasible, "fraction": 0.0, "n_bins": n_bins}

    ix = np.clip(
        np.digitize(in_feasible_points, edges) - 1, 0, n_bins - 1,
    )  # shape (n, k)
    flat = np.zeros(len(in_feasible_points), dtype=int)
    for j in range(k):
        flat = flat * n_bins + ix[:, j]

    occupied_ids = set(int(f) for f in flat)
    feasible_ids = set(int(i) for i, m in enumerate(feasible_mask) if m)
    occupied = len(occupied_ids & feasible_ids)

    return {
        "occupied": int(occupied),
        "feasible": n_feasible,
        "fraction": float(occupied) / max(1, n_feasible),
        "n_bins": int(n_bins),
    }


def _interior_fraction(
    points: np.ndarray,
    interior_eps: float,
    shape: str = "cube",
    radius: float = FEASIBLE_RADIUS,
) -> float:
    if len(points) == 0:
        return 0.0
    inside = points[_in_feasible(points, shape, radius)]
    if len(inside) == 0:
        return 0.0
    return float((_distance_from_boundary(inside, shape, radius) > interior_eps).mean())


def run_borg(
    k: int,
    nfe: int,
    seed: int,
    epsilon: float,
    shape: str = "cube",
    radius: float = FEASIBLE_RADIUS,
) -> np.ndarray:
    """Run EpsNSGAII on the k-dim constrained analytic problem.

    Penalty strictly dominates every feasible objective vector: feasible
    ``J_i in [-R, R]`` and ``J_{k+1} in [0, 2k*DV_HI]``, so penalty at
    ``[10, ..., 10, 30k]`` is always worse.
    """
    np.random.seed(seed)
    anti_ideal = _anti_ideal(k)
    penalty = np.array([10.0] * k + [30.0 * k])

    def evaluate(variables):
        dvs = np.array([float(v) for v in variables])
        if not _in_feasible(dvs[None, :], shape, radius)[0]:
            return penalty.tolist()
        return analytic_objectives(dvs, anti_ideal).tolist()

    problem = Problem(k, k + 1)
    for i in range(k):
        problem.types[i] = Real(DV_LO, DV_HI)
    problem.function = evaluate

    algo = EpsNSGAII(problem, epsilons=[epsilon] * (k + 1))
    t0 = time.perf_counter()
    algo.run(nfe)
    elapsed = time.perf_counter() - t0

    all_dvs = np.array([
        [float(v) for v in s.variables[:]] for s in algo.result
    ]) if len(algo.result) else np.zeros((0, k))
    dvs = (
        all_dvs[_in_feasible(all_dvs, shape, radius)]
        if len(all_dvs) else all_dvs
    )

    print(
        f"  borg(k={k}, shape={shape}): nfe={nfe} seed={seed} eps={epsilon} "
        f"archive={len(algo.result)} feasible={len(dvs)} "
        f"elapsed={elapsed:.1f}s"
    )
    return dvs


def _boundary_2d(shape: str, radius: float) -> tuple:
    """Closed (x, y) polyline tracing the feasibility boundary in 2D."""
    if shape == "ball":
        theta = np.linspace(0, 2 * np.pi, 200)
        return radius * np.cos(theta), radius * np.sin(theta)
    # cube: walk the four edges
    x = np.array([-radius, radius, radius, -radius, -radius])
    y = np.array([-radius, -radius, radius, radius, -radius])
    return x, y


def plot_low_d(
    k: int,
    borg: np.ndarray,
    unif: np.ndarray,
    lhs: np.ndarray,
    sobol: np.ndarray,
    out_path: Path,
    shape: str = "cube",
    radius: float = FEASIBLE_RADIUS,
) -> None:
    """2x3 panel for k in {2, 3}: scatter row + metrics row."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(12, 8))
    shape_label = shape  # "ball" or "cube"

    if k == 2:
        axes_top = [fig.add_subplot(2, 3, i + 1) for i in range(3)]
        bx, by = _boundary_2d(shape, radius)
        for ax, pts, color, name in zip(
            axes_top,
            [borg, lhs, sobol],
            ["tab:blue", "tab:orange", "tab:green"],
            ["MOEA-FIND archive",
             f"LHS in {shape_label}",
             f"Sobol in {shape_label}"],
        ):
            ax.plot(bx, by, "k--", lw=1.0)
            ax.scatter(*_anti_ideal(2), marker="X", color="red", s=60)
            ax.scatter(pts[:, 0], pts[:, 1], s=10, color=color)
            ax.set_xlim(-3.2, 3.5)
            ax.set_ylim(-3.2, 3.5)
            ax.set_aspect("equal")
            ax.set_xlabel(r"$X_1$")
            ax.set_ylabel(r"$X_2$")
            ax.set_title(f"{name} (n={len(pts)})")
    elif k == 3:
        for i, (pts, color, name) in enumerate(zip(
            [borg, lhs, sobol],
            ["tab:blue", "tab:orange", "tab:green"],
            ["MOEA-FIND archive",
             f"LHS in {shape_label}",
             f"Sobol in {shape_label}"],
        )):
            ax = fig.add_subplot(2, 3, i + 1, projection="3d")
            ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=6, color=color)
            ax.scatter(*_anti_ideal(3), marker="X", color="red", s=60)
            ax.set_xlim(-3, 3)
            ax.set_ylim(-3, 3)
            ax.set_zlim(-3, 3)
            ax.set_xlabel(r"$X_1$")
            ax.set_ylabel(r"$X_2$")
            ax.set_zlabel(r"$X_3$")
            ax.set_title(f"{name} (n={len(pts)})")

    # Metrics row.
    ax_d = fig.add_subplot(2, 3, 4)
    bins = np.linspace(0, 4 * k, 25)
    d_borg = _distance_from_antiideal_l1(borg)
    d_lhs = _distance_from_antiideal_l1(lhs)
    d_sobol = _distance_from_antiideal_l1(sobol)
    ax_d.hist(d_borg, bins=bins, alpha=0.6, color="tab:blue", label="MOEA-FIND")
    ax_d.hist(d_lhs, bins=bins, alpha=0.4, color="tab:orange", label="LHS")
    ax_d.hist(d_sobol, bins=bins, alpha=0.4, color="tab:green", label="Sobol")
    ax_d.set_xlabel(r"Manhattan distance from $D^*$")
    ax_d.set_ylabel("count")
    ax_d.set_title("Distance from anti-ideal")
    ax_d.legend(fontsize=8)

    ax_b = fig.add_subplot(2, 3, 5)
    bins_b = np.linspace(0, radius, 20)
    b_borg = _distance_from_boundary(
        borg[_in_feasible(borg, shape, radius)], shape, radius)
    b_lhs = _distance_from_boundary(
        lhs[_in_feasible(lhs, shape, radius)], shape, radius)
    b_sobol = _distance_from_boundary(
        sobol[_in_feasible(sobol, shape, radius)], shape, radius)
    ax_b.hist(b_borg, bins=bins_b, alpha=0.6, color="tab:blue", label="MOEA-FIND")
    ax_b.hist(b_lhs, bins=bins_b, alpha=0.4, color="tab:orange", label="LHS")
    ax_b.hist(b_sobol, bins=bins_b, alpha=0.4, color="tab:green", label="Sobol")
    ax_b.set_xlabel(f"distance from {shape_label} boundary")
    ax_b.set_ylabel("count")
    ax_b.set_title("Interior versus shell mass")
    ax_b.legend(fontsize=8)

    ax_g = fig.add_subplot(2, 3, 6)
    occ = {
        "MOEA-FIND": _grid_occupancy(borg, k, 12 if k <= 3 else 6, shape, radius),
        "uniform":   _grid_occupancy(unif, k, 12 if k <= 3 else 6, shape, radius),
        "LHS":       _grid_occupancy(lhs, k, 12 if k <= 3 else 6, shape, radius),
        "Sobol":     _grid_occupancy(sobol, k, 12 if k <= 3 else 6, shape, radius),
    }
    labels = list(occ.keys())
    values = [occ[x]["fraction"] for x in labels]
    colors = ["tab:blue", "tab:gray", "tab:orange", "tab:green"]
    ax_g.bar(labels, values, color=colors)
    ax_g.set_ylabel("occupied / feasible cells")
    ax_g.set_ylim(0, 1.05)
    ax_g.set_title("Grid cell coverage")
    for i, v in enumerate(values):
        ax_g.text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)

    fig.suptitle(
        f"Shell versus interior diagnostic, k={k} "
        f"({shape_label} half-width {radius}, anti-ideal at the corner)"
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_high_d(
    k: int,
    borg: np.ndarray,
    unif: np.ndarray,
    lhs: np.ndarray,
    sobol: np.ndarray,
    out_path: Path,
    shape: str = "cube",
    radius: float = FEASIBLE_RADIUS,
) -> None:
    """Metrics-only 1x4 panel for k >= 4, since no honest 4+D scatter exists."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))

    # Distance from anti-ideal.
    bins = np.linspace(0, 4 * k, 25)
    d_borg = _distance_from_antiideal_l1(borg)
    d_lhs = _distance_from_antiideal_l1(lhs)
    d_sobol = _distance_from_antiideal_l1(sobol)
    axes[0].hist(d_borg, bins=bins, alpha=0.6, color="tab:blue", label="MOEA-FIND")
    axes[0].hist(d_lhs, bins=bins, alpha=0.4, color="tab:orange", label="LHS")
    axes[0].hist(d_sobol, bins=bins, alpha=0.4, color="tab:green", label="Sobol")
    axes[0].set_xlabel(r"Manhattan distance from $D^*$")
    axes[0].set_ylabel("count")
    axes[0].set_title("Distance from anti-ideal")
    axes[0].legend(fontsize=8)

    # Distance from feasibility boundary.
    bins_b = np.linspace(0, radius, 20)
    b_borg = _distance_from_boundary(
        borg[_in_feasible(borg, shape, radius)], shape, radius)
    b_lhs = _distance_from_boundary(
        lhs[_in_feasible(lhs, shape, radius)], shape, radius)
    b_sobol = _distance_from_boundary(
        sobol[_in_feasible(sobol, shape, radius)], shape, radius)
    axes[1].hist(b_borg, bins=bins_b, alpha=0.6, color="tab:blue", label="MOEA-FIND")
    axes[1].hist(b_lhs, bins=bins_b, alpha=0.4, color="tab:orange", label="LHS")
    axes[1].hist(b_sobol, bins=bins_b, alpha=0.4, color="tab:green", label="Sobol")
    axes[1].set_xlabel(f"distance from {shape} boundary")
    axes[1].set_ylabel("count")
    axes[1].set_title("Interior versus shell mass")
    axes[1].legend(fontsize=8)

    # Orthant occupancy.
    orth = {
        "MOEA-FIND": _orthant_occupancy(borg),
        "uniform": _orthant_occupancy(unif),
        "LHS": _orthant_occupancy(lhs),
        "Sobol": _orthant_occupancy(sobol),
    }
    labels = list(orth.keys())
    values = [orth[x]["fraction"] for x in labels]
    colors = ["tab:blue", "tab:gray", "tab:orange", "tab:green"]
    axes[2].bar(labels, values, color=colors)
    axes[2].set_ylabel(f"occupied / {1 << k} orthants")
    axes[2].set_ylim(0, 1.05)
    axes[2].set_title("Orthant occupancy")
    for i, v in enumerate(values):
        axes[2].text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)

    # Grid cell occupancy on a coarser partition.
    n_bins = 6 if k <= 4 else 4
    occ = {
        "MOEA-FIND": _grid_occupancy(borg, k, n_bins, shape, radius),
        "uniform":   _grid_occupancy(unif, k, n_bins, shape, radius),
        "LHS":       _grid_occupancy(lhs, k, n_bins, shape, radius),
        "Sobol":     _grid_occupancy(sobol, k, n_bins, shape, radius),
    }
    labels = list(occ.keys())
    values = [occ[x]["fraction"] for x in labels]
    axes[3].bar(labels, values, color=colors)
    axes[3].set_ylabel(f"occupied / feasible cells ({n_bins}$^{k}$)")
    axes[3].set_ylim(0, 1.05)
    axes[3].set_title("Grid cell coverage")
    for i, v in enumerate(values):
        axes[3].text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)

    fig.suptitle(
        f"Shell versus interior diagnostic, k={k} "
        f"({shape} half-width {radius}, anti-ideal at the corner)"
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def summarize(
    k: int,
    borg: np.ndarray,
    unif: np.ndarray,
    lhs: np.ndarray,
    sobol: np.ndarray,
    interior_eps: float,
    shape: str = "cube",
    radius: float = FEASIBLE_RADIUS,
) -> dict:
    def pack(arr):
        d = _distance_from_antiideal_l1(arr)
        return {
            "n": int(len(arr)),
            "dist_from_D*": {
                "mean": float(d.mean()) if len(d) else None,
                "std": float(d.std()) if len(d) else None,
                "min": float(d.min()) if len(d) else None,
                "max": float(d.max()) if len(d) else None,
            },
            "interior_fraction": _interior_fraction(
                arr, interior_eps, shape, radius),
            "orthant_occupancy": _orthant_occupancy(arr),
            "grid_occupancy": _grid_occupancy(
                arr, k, 6 if k >= 4 else 12, shape, radius),
        }

    # Reference sampler labels track the chosen feasibility shape for
    # downstream plotting and the summary table.
    ref_key = f"in_{shape}"
    return {
        "k": k,
        "dv_range": [DV_LO, DV_HI],
        "feasible_shape": shape,
        "feasible_radius": radius,
        "MOEA-FIND": pack(borg),
        f"uniform_{ref_key}": pack(unif),
        f"lhs_{ref_key}": pack(lhs),
        f"sobol_{ref_key}": pack(sobol),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--k", type=int, default=2, help="decision dimensionality")
    p.add_argument("--nfe", type=int, default=30_000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--epsilon", type=float, default=None,
                   help="epsilon-box side (default 0.10 for k<=3, else 0.15+)")
    p.add_argument("--interior-eps", type=float, default=0.25)
    p.add_argument("--feasible-shape", choices=VALID_SHAPES, default="cube",
                   help="Feasible-region shape (default: cube — K-dim square)")
    p.add_argument("--feasible-radius", type=float, default=FEASIBLE_RADIUS,
                   help="Ball radius or cube half-width")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    p.add_argument("--figure-dir", type=Path,
                   default=PROJECT_ROOT / "figures")
    args = p.parse_args()

    shape = args.feasible_shape
    radius = args.feasible_radius

    eps = args.epsilon
    if eps is None:
        eps = {2: 0.10, 3: 0.15, 4: 0.20, 5: 0.25, 6: 0.30, 7: 0.35}.get(
            args.k, 0.4,
        )

    k_dir = args.output_dir / f"k{args.k}"
    k_dir.mkdir(parents=True, exist_ok=True)
    (k_dir / "config.json").write_text(json.dumps({
        "script": "diag_shell_vs_interior.py",
        "k": args.k,
        "dv_range": [DV_LO, DV_HI],
        "anti_ideal": _anti_ideal(args.k).tolist(),
        "feasible_shape": shape,
        "feasible_radius": radius,
        "nfe": args.nfe,
        "seed": args.seed,
        "epsilon": eps,
        "interior_eps": args.interior_eps,
    }, indent=2))

    print(
        f"[diag k={args.k}] shape={shape} radius={radius} "
        f"anti-ideal={_anti_ideal(args.k).tolist()}"
    )

    borg = run_borg(args.k, args.nfe, args.seed, eps, shape, radius)
    if len(borg) == 0:
        print(f"[diag k={args.k}] ERROR: no feasible Borg solutions.")
        return

    n = len(borg)
    unif = _sample_in_feasible(n, "uniform", args.seed, args.k, shape, radius)
    lhs = _sample_in_feasible(n, "lhs", args.seed, args.k, shape, radius)
    sobol = _sample_in_feasible(n, "sobol", args.seed, args.k, shape, radius)

    stats = summarize(
        args.k, borg, unif, lhs, sobol, args.interior_eps, shape, radius
    )
    (k_dir / "results.json").write_text(
        json.dumps(stats, indent=2, default=float),
    )
    np.savez(k_dir / "samples.npz",
             borg=borg, unif=unif, lhs=lhs, sobol=sobol)

    fig_path = args.figure_dir / f"figSI_shell_interior_k{args.k}.pdf"
    if args.k <= 3:
        plot_low_d(args.k, borg, unif, lhs, sobol, fig_path, shape, radius)
    else:
        plot_high_d(args.k, borg, unif, lhs, sobol, fig_path, shape, radius)
    print(f"[diag k={args.k}] figure: {fig_path}")
    print(json.dumps(stats, indent=2, default=float))


if __name__ == "__main__":
    main()
