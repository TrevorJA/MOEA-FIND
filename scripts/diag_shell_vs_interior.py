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

PROJECT_ROOT = Path(__file__).resolve().parents[1]
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


# Global problem constants. DV_LO, DV_HI, and DISK_RADIUS are kept fixed
# across k so that results from different dimensionalities are directly
# comparable; the bounding box side is 6, the k-ball radius is 2.5, and
# the anti-ideal sits at the positive corner.
DV_LO, DV_HI = -3.0, 3.0
DISK_RADIUS = 2.5
OUTPUT_SLUG = "diag_shell_vs_interior"


def _in_ball(x: np.ndarray, radius: float = DISK_RADIUS) -> np.ndarray:
    """Boolean mask indicating which rows of x lie inside the k-ball."""
    return (x ** 2).sum(axis=-1) <= radius ** 2


def _reject_into_ball(
    samples: np.ndarray,
    radius: float = DISK_RADIUS,
) -> np.ndarray:
    """Return only rows of samples that fall inside the k-ball."""
    return samples[_in_ball(samples, radius)]


def _sample_in_ball(
    n_target: int,
    method: str,
    seed: int,
    k: int,
    radius: float = DISK_RADIUS,
    oversample: int = 20,
) -> np.ndarray:
    """Draw n_target samples inside the k-ball by rejection.

    Uses method in {"uniform", "lhs", "sobol"}. The rejection rate is
    (cube volume) / (ball volume), which grows quickly with k, so the
    oversample factor is automatically increased if the first pass
    falls short.
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

    accepted = _reject_into_ball(raw, radius)
    if len(accepted) < n_target:
        return _sample_in_ball(
            n_target, method, seed + 1, k, radius,
            oversample=oversample * 4,
        )
    return accepted[:n_target]


def _anti_ideal(k: int) -> np.ndarray:
    return np.full(k, DV_HI)


def _distance_from_antiideal_l1(points: np.ndarray) -> np.ndarray:
    if len(points) == 0:
        return np.zeros(0)
    return np.abs(points - _anti_ideal(points.shape[1])).sum(axis=1)


def _distance_from_ball_boundary(
    points: np.ndarray,
    radius: float = DISK_RADIUS,
) -> np.ndarray:
    """Positive inside, zero on boundary, negative outside."""
    r = np.linalg.norm(points, axis=1)
    return radius - r


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
    radius: float = DISK_RADIUS,
) -> dict:
    """Feasible-cell occupancy on a coarse regular partition.

    The full grid has n_bins^k cells, of which only those whose center
    lies inside the ball are counted as feasible. Occupancy is the
    number of distinct feasible cells that contain at least one
    sample point, divided by the number of feasible cells.

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
    feasible_mask = _in_ball(cells, radius)
    n_feasible = int(feasible_mask.sum())

    in_ball_points = points[_in_ball(points, radius)] if len(points) else points
    if len(in_ball_points) == 0:
        return {"occupied": 0, "feasible": n_feasible, "fraction": 0.0, "n_bins": n_bins}

    ix = np.clip(
        np.digitize(in_ball_points, edges) - 1, 0, n_bins - 1,
    )  # shape (n, k)
    flat = np.zeros(len(in_ball_points), dtype=int)
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


def _interior_fraction(points: np.ndarray, interior_eps: float) -> float:
    if len(points) == 0:
        return 0.0
    in_ball = points[_in_ball(points)]
    if len(in_ball) == 0:
        return 0.0
    return float((_distance_from_ball_boundary(in_ball) > interior_eps).mean())


def run_borg(k: int, nfe: int, seed: int, epsilon: float) -> np.ndarray:
    """Run EpsNSGAII on the k-dimensional constrained analytic problem."""
    np.random.seed(seed)
    anti_ideal = _anti_ideal(k)
    # Penalty strictly dominates every feasible objective vector. Feasible
    # bounds: J_i in [-R, R] = [-2.5, 2.5], J_{k+1} in [0, 2k*DV_HI] bounded
    # above by 2k*3 = 6k. Penalty at 10 and 30k always strictly worse.
    penalty = np.array([10.0] * k + [30.0 * k])

    def evaluate(variables):
        dvs = np.array([float(v) for v in variables])
        if (dvs ** 2).sum() > DISK_RADIUS ** 2:
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

    feasible_dvs = []
    for s in algo.result:
        dvs = np.array([float(v) for v in s.variables[:]])
        if (dvs ** 2).sum() <= DISK_RADIUS ** 2:
            feasible_dvs.append(dvs)
    dvs = np.array(feasible_dvs) if feasible_dvs else np.zeros((0, k))

    print(
        f"  borg(k={k}): nfe={nfe} seed={seed} eps={epsilon} "
        f"archive={len(algo.result)} feasible={len(dvs)} "
        f"elapsed={elapsed:.1f}s"
    )
    return dvs


def plot_low_d(
    k: int,
    borg: np.ndarray,
    unif: np.ndarray,
    lhs: np.ndarray,
    sobol: np.ndarray,
    out_path: Path,
) -> None:
    """2x3 panel for k in {2, 3}: scatter row + metrics row."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(12, 8))

    if k == 2:
        axes_top = [fig.add_subplot(2, 3, i + 1) for i in range(3)]
        theta = np.linspace(0, 2 * np.pi, 200)
        cx = DISK_RADIUS * np.cos(theta)
        cy = DISK_RADIUS * np.sin(theta)
        for ax, pts, color, name in zip(
            axes_top,
            [borg, lhs, sobol],
            ["tab:blue", "tab:orange", "tab:green"],
            ["MOEA-FIND archive", "LHS in ball", "Sobol in ball"],
        ):
            ax.plot(cx, cy, "k--", lw=1.0)
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
            ["MOEA-FIND archive", "LHS in ball", "Sobol in ball"],
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
    bins_b = np.linspace(0, DISK_RADIUS, 20)
    b_borg = _distance_from_ball_boundary(borg[_in_ball(borg)])
    b_lhs = _distance_from_ball_boundary(lhs[_in_ball(lhs)])
    b_sobol = _distance_from_ball_boundary(sobol[_in_ball(sobol)])
    ax_b.hist(b_borg, bins=bins_b, alpha=0.6, color="tab:blue", label="MOEA-FIND")
    ax_b.hist(b_lhs, bins=bins_b, alpha=0.4, color="tab:orange", label="LHS")
    ax_b.hist(b_sobol, bins=bins_b, alpha=0.4, color="tab:green", label="Sobol")
    ax_b.set_xlabel("distance from ball boundary")
    ax_b.set_ylabel("count")
    ax_b.set_title("Interior versus shell mass")
    ax_b.legend(fontsize=8)

    ax_g = fig.add_subplot(2, 3, 6)
    occ = {
        "MOEA-FIND": _grid_occupancy(borg, k, 12 if k <= 3 else 6),
        "uniform": _grid_occupancy(unif, k, 12 if k <= 3 else 6),
        "LHS": _grid_occupancy(lhs, k, 12 if k <= 3 else 6),
        "Sobol": _grid_occupancy(sobol, k, 12 if k <= 3 else 6),
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
        f"(ball radius {DISK_RADIUS}, anti-ideal at the corner)"
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

    # Distance from ball boundary.
    bins_b = np.linspace(0, DISK_RADIUS, 20)
    b_borg = _distance_from_ball_boundary(borg[_in_ball(borg)])
    b_lhs = _distance_from_ball_boundary(lhs[_in_ball(lhs)])
    b_sobol = _distance_from_ball_boundary(sobol[_in_ball(sobol)])
    axes[1].hist(b_borg, bins=bins_b, alpha=0.6, color="tab:blue", label="MOEA-FIND")
    axes[1].hist(b_lhs, bins=bins_b, alpha=0.4, color="tab:orange", label="LHS")
    axes[1].hist(b_sobol, bins=bins_b, alpha=0.4, color="tab:green", label="Sobol")
    axes[1].set_xlabel("distance from ball boundary")
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

    # Grid cell occupancy on a courser partition.
    n_bins = 6 if k <= 4 else 4
    occ = {
        "MOEA-FIND": _grid_occupancy(borg, k, n_bins),
        "uniform": _grid_occupancy(unif, k, n_bins),
        "LHS": _grid_occupancy(lhs, k, n_bins),
        "Sobol": _grid_occupancy(sobol, k, n_bins),
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
        f"(ball radius {DISK_RADIUS}, anti-ideal at the corner)"
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
            "interior_fraction": _interior_fraction(arr, interior_eps),
            "orthant_occupancy": _orthant_occupancy(arr),
            "grid_occupancy": _grid_occupancy(arr, k, 6 if k >= 4 else 12),
        }

    return {
        "k": k,
        "dv_range": [DV_LO, DV_HI],
        "ball_radius": DISK_RADIUS,
        "MOEA-FIND": pack(borg),
        "uniform_in_ball": pack(unif),
        "lhs_in_ball": pack(lhs),
        "sobol_in_ball": pack(sobol),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--k", type=int, default=2, help="decision dimensionality")
    p.add_argument("--nfe", type=int, default=30_000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--epsilon", type=float, default=None,
                   help="epsilon-box side (default 0.10 for k<=3, else 0.15+)")
    p.add_argument("--interior-eps", type=float, default=0.25)
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    p.add_argument("--figure-dir", type=Path,
                   default=PROJECT_ROOT / "figures")
    args = p.parse_args()

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
        "ball_radius": DISK_RADIUS,
        "nfe": args.nfe,
        "seed": args.seed,
        "epsilon": eps,
        "interior_eps": args.interior_eps,
    }, indent=2))

    print(f"[diag k={args.k}] ball radius={DISK_RADIUS} anti-ideal={_anti_ideal(args.k).tolist()}")

    borg = run_borg(args.k, args.nfe, args.seed, eps)
    if len(borg) == 0:
        print(f"[diag k={args.k}] ERROR: no feasible Borg solutions.")
        return

    n = len(borg)
    unif = _sample_in_ball(n, "uniform", args.seed, args.k)
    lhs = _sample_in_ball(n, "lhs", args.seed, args.k)
    sobol = _sample_in_ball(n, "sobol", args.seed, args.k)

    stats = summarize(args.k, borg, unif, lhs, sobol, args.interior_eps)
    (k_dir / "results.json").write_text(
        json.dumps(stats, indent=2, default=float),
    )
    np.savez(k_dir / "samples.npz",
             borg=borg, unif=unif, lhs=lhs, sobol=sobol)

    fig_path = args.figure_dir / f"figSI_shell_interior_k{args.k}.pdf"
    if args.k <= 3:
        plot_low_d(args.k, borg, unif, lhs, sobol, fig_path)
    else:
        plot_high_d(args.k, borg, unif, lhs, sobol, fig_path)
    print(f"[diag k={args.k}] figure: {fig_path}")
    print(json.dumps(stats, indent=2, default=float))


if __name__ == "__main__":
    main()
