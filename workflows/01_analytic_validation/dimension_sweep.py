"""K-dimensional analytic dimension sweep (manuscript Fig 4).

Tests interior-filling coverage of the L1 + epsilon-tile MOEA-FIND device
on a K-dimensional analytic problem with a constrained feasible region.
The default feasible region is the **hypercube** ``max(|X_i|) <= R``;
the K-ball variant is preserved via ``--feasible-shape ball`` for
historical comparison but is no longer the manuscript default.

Optimization is MM Borg MOEA (launched via MPI). Per-K outputs under
``outputs/01_analytic_validation/dimension_sweep/k{K}/``:
    - config.json
    - results.json (summary statistics for borg + uniform/lhs/sobol references)
    - samples.npz  (raw point clouds for downstream plotting)

The companion plotting driver
``workflows/99_manuscript_figures/dimension_sweep.py`` reads
every ``k{K}`` directory and produces the K=1..6 sweep figure.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

_stub = MagicMock()
sys.modules.setdefault("synhydro", _stub)
sys.modules.setdefault("synhydro.droughts", _stub.droughts)
sys.modules.setdefault("synhydro.droughts.ssi", _stub.droughts.ssi)

from src.discovery.analysis import generate_lhs_samples, generate_sobol_samples  # noqa: E402
from src.optimization.borg_runner import run_optimization  # noqa: E402
from src.metrics.objectives import analytic_objectives  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "01_analytic_validation"
DRIVER = "dimension_sweep"

DV_LO, DV_HI = -3.0, 3.0
FEASIBLE_RADIUS = 2.5
VALID_SHAPES = ("cube", "ball")


def _in_feasible(x, shape="cube", radius=FEASIBLE_RADIUS):
    if shape == "cube":
        return np.max(np.abs(x), axis=-1) <= radius
    if shape == "ball":
        return (x ** 2).sum(axis=-1) <= radius ** 2
    raise ValueError(f"unknown feasible shape {shape!r}")


def _reject_into_feasible(samples, shape="cube", radius=FEASIBLE_RADIUS):
    return samples[_in_feasible(samples, shape, radius)]


def _sample_in_feasible(n_target, method, seed, k, shape="cube",
                        radius=FEASIBLE_RADIUS, oversample=20):
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
        return _sample_in_feasible(n_target, method, seed + 1, k, shape, radius,
                                   oversample=oversample * 4)
    return accepted[:n_target]


def _anti_ideal(k):
    return np.full(k, DV_HI)


def _distance_from_antiideal_l1(points):
    if len(points) == 0:
        return np.zeros(0)
    return np.abs(points - _anti_ideal(points.shape[1])).sum(axis=1)


def _distance_from_boundary(points, shape="cube", radius=FEASIBLE_RADIUS):
    if shape == "cube":
        return radius - np.max(np.abs(points), axis=1)
    if shape == "ball":
        return radius - np.linalg.norm(points, axis=1)
    raise ValueError(shape)


def _orthant_occupancy(points):
    if len(points) == 0:
        return {"occupied": 0, "total": 0, "fraction": 0.0}
    k = points.shape[1]
    signs = (points > 0).astype(int)
    ids = np.zeros(len(points), dtype=int)
    for j in range(k):
        ids += signs[:, j] << j
    occupied = len(np.unique(ids))
    total = 1 << k
    return {"occupied": int(occupied), "total": int(total),
            "fraction": float(occupied) / total}


def _grid_occupancy(points, k, n_bins, shape="cube", radius=FEASIBLE_RADIUS):
    max_cells = 100_000
    while n_bins ** k > max_cells and n_bins > 2:
        n_bins -= 1
    edges = np.linspace(-radius, radius, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mesh = np.meshgrid(*([centers] * k), indexing="ij")
    cells = np.stack([m.ravel() for m in mesh], axis=1)
    feasible_mask = _in_feasible(cells, shape, radius)
    n_feasible = int(feasible_mask.sum())
    in_pts = points[_in_feasible(points, shape, radius)] if len(points) else points
    if len(in_pts) == 0:
        return {"occupied": 0, "feasible": n_feasible, "fraction": 0.0, "n_bins": n_bins}
    ix = np.clip(np.digitize(in_pts, edges) - 1, 0, n_bins - 1)
    flat = np.zeros(len(in_pts), dtype=int)
    for j in range(k):
        flat = flat * n_bins + ix[:, j]
    occupied_ids = set(int(f) for f in flat)
    feasible_ids = set(int(i) for i, m in enumerate(feasible_mask) if m)
    return {"occupied": len(occupied_ids & feasible_ids),
            "feasible": n_feasible,
            "fraction": float(len(occupied_ids & feasible_ids)) / max(1, n_feasible),
            "n_bins": int(n_bins)}


def _interior_fraction(points, interior_eps, shape, radius):
    if len(points) == 0:
        return 0.0
    inside = points[_in_feasible(points, shape, radius)]
    if len(inside) == 0:
        return 0.0
    return float((_distance_from_boundary(inside, shape, radius) > interior_eps).mean())


def run_borg(k, nfe, seed, epsilon, shape, radius, output_dir):
    """MM Borg run with infeasibility encoded as a hard penalty objective."""
    np.random.seed(seed)
    anti_ideal = _anti_ideal(k)
    penalty = np.array([10.0] * k + [30.0 * k])

    def evaluate(unit_dvs):
        # Borg DVs in [0,1]; scale to physical [-3,3] domain.
        dvs = DV_LO + unit_dvs * (DV_HI - DV_LO)
        if not _in_feasible(dvs[None, :], shape, radius)[0]:
            return penalty.tolist(), []
        return analytic_objectives(dvs, anti_ideal).tolist(), []

    opt = run_optimization(
        algorithm="borg_mm",
        evaluate=evaluate,
        n_dvs=k,
        n_objs=k + 1,
        n_constrs=0,
        epsilons=[epsilon] * (k + 1),
        nfe=nfe,
        seed=seed,
        output_dir=output_dir,
        # n_islands auto-picked by borg_runner._auto_islands
    )

    if opt.pareto_dvs.shape[0] == 0:
        return np.zeros((0, k))

    physical = DV_LO + opt.pareto_dvs * (DV_HI - DV_LO)
    feasible = physical[_in_feasible(physical, shape, radius)]
    print(f"  borg(k={k}, shape={shape}): nfe={nfe} seed={seed} eps={epsilon} "
          f"archive={len(physical)} feasible={len(feasible)} "
          f"elapsed={opt.elapsed_s:.1f}s")
    return feasible


def summarize(k, borg, unif, lhs, sobol, interior_eps, shape, radius):
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
            "interior_fraction": _interior_fraction(arr, interior_eps, shape, radius),
            "orthant_occupancy": _orthant_occupancy(arr),
            "grid_occupancy": _grid_occupancy(arr, k, 6 if k >= 4 else 12, shape, radius),
        }
    ref_key = f"in_{shape}"
    return {
        "k": k,
        "algorithm": "borg_mm",
        "dv_range": [DV_LO, DV_HI],
        "feasible_shape": shape,
        "feasible_radius": radius,
        "MOEA-FIND": pack(borg),
        f"uniform_{ref_key}": pack(unif),
        f"lhs_{ref_key}": pack(lhs),
        f"sobol_{ref_key}": pack(sobol),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--k", type=int, required=True)
    p.add_argument("--nfe", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--epsilon", type=float, required=True)
    p.add_argument("--interior-eps", type=float, default=0.25)
    p.add_argument("--feasible-shape", choices=VALID_SHAPES, default="cube")
    p.add_argument("--feasible-radius", type=float, default=FEASIBLE_RADIUS)
    args = p.parse_args()

    shape = args.feasible_shape
    radius = args.feasible_radius
    k_dir = stage_output_dir(STAGE, DRIVER, f"k{args.k}")
    (k_dir / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "algorithm": "borg_mm",
        "k": args.k, "dv_range": [DV_LO, DV_HI],
        "anti_ideal": _anti_ideal(args.k).tolist(),
        "feasible_shape": shape, "feasible_radius": radius,
        "nfe": args.nfe, "seed": args.seed, "epsilon": args.epsilon,
        "interior_eps": args.interior_eps,
    }, indent=2))

    print(f"[01/dimension_sweep k={args.k}] shape={shape} radius={radius}")
    borg = run_borg(args.k, args.nfe, args.seed, args.epsilon, shape, radius, k_dir)
    if len(borg) == 0:
        print(f"[01/dimension_sweep k={args.k}] ERROR: no feasible Borg solutions.")
        return
    n = len(borg)
    unif = _sample_in_feasible(n, "uniform", args.seed, args.k, shape, radius)
    lhs = _sample_in_feasible(n, "lhs", args.seed, args.k, shape, radius)
    sobol = _sample_in_feasible(n, "sobol", args.seed, args.k, shape, radius)
    stats = summarize(args.k, borg, unif, lhs, sobol, args.interior_eps, shape, radius)
    (k_dir / "results.json").write_text(json.dumps(stats, indent=2, default=float))
    np.savez(k_dir / "samples.npz", borg=borg, unif=unif, lhs=lhs, sobol=sobol)
    print(f"  outputs: {k_dir}")


if __name__ == "__main__":
    main()
