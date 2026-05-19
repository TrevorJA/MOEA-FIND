"""Kirsch ensemble convergence analysis (compute only; SI-D).

Subsamples a pre-generated Kirsch library at varying sizes and measures
convergence of drought-characteristic range and coverage metrics; if a
MOEA-FIND results.json is provided, the Pareto front is added as a
reference point.

Outputs under outputs/02_calibration/kirsch_convergence/:
    - config.json
    - convergence.json (metrics at each ensemble size x seed)

Plotting lives in workflows/99_supporting_info_figures/kirsch_convergence.py.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.discovery.analysis import coverage_metrics  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "02_calibration"
DRIVER = "kirsch_convergence"

ENSEMBLE_SIZES = [50, 100, 200, 500, 1000, 2000, 5000, 10000]
N_SEEDS = 10


def load_kirsch_library(path: Path) -> np.ndarray:
    """Load Kirsch library objectives from NPZ or JSON.

    Returns array of shape (n_traces, 2) with [mean_duration, mean_avg_severity].
    """
    path = Path(path)
    if path.suffix == ".npz":
        data = np.load(path, allow_pickle=True)
        return data["objectives"]
    elif path.suffix == ".json":
        chars = json.loads(path.read_text())
        return np.array([
            [c["mean_duration"], c["mean_avg_severity"]] for c in chars
        ])
    raise ValueError(f"Unknown format: {path.suffix}")


def load_moea_front(path: Path) -> Optional[np.ndarray]:
    """Load MOEA-FIND Pareto front objectives from results.json.

    Returns array of shape (n_pareto, 2) or None if no solutions.
    """
    result = json.loads(Path(path).read_text())
    dm = result.get("drought_metrics", [])
    if not dm:
        return None
    return np.array(dm)[:, :2]  # first 2 columns: mean_duration, mean_avg_severity


def compute_range_coverage(subset: np.ndarray, reference: np.ndarray) -> Dict:
    """Fraction of reference range covered by subset, per objective."""
    ref_ranges = reference.max(axis=0) - reference.min(axis=0)
    sub_ranges = subset.max(axis=0) - subset.min(axis=0)
    # Avoid division by zero
    safe_ref = np.where(ref_ranges > 0, ref_ranges, 1.0)
    frac = np.clip(sub_ranges / safe_ref, 0.0, 1.0)
    return {
        "mean_duration_frac": float(frac[0]),
        "mean_avg_severity_frac": float(frac[1]),
        "overall_frac": float(frac.mean()),
        "subset_range_duration": [float(subset[:, 0].min()), float(subset[:, 0].max())],
        "subset_range_severity": [float(subset[:, 1].min()), float(subset[:, 1].max())],
    }


def run_convergence(
    library: np.ndarray,
    moea_front: Optional[np.ndarray],
    ensemble_sizes: List[int],
    n_seeds: int,
    base_seed: int,
) -> List[Dict]:
    """Compute convergence metrics at each ensemble size × seed."""
    rng = np.random.default_rng(base_seed)
    n_library = len(library)

    # Bounds for coverage metrics (from full library extent)
    lb = library.min(axis=0)
    ub = library.max(axis=0)

    results = []
    for N in ensemble_sizes:
        if N > n_library:
            print(f"  skipping N={N} (library has only {n_library} traces)")
            continue
        for s in range(n_seeds):
            idx = rng.choice(n_library, size=N, replace=False)
            subset = library[idx]

            cm = coverage_metrics(subset, lb, ub)
            entry = {
                "ensemble_size": N,
                "seed": s,
                "n_points": N,
                "L2_star": cm.get("L2_star_discrepancy", float("nan")),
                "nn_cv": cm.get("nn_cv", float("nan")),
                "nn_mean": cm.get("nn_mean", float("nan")),
                "range_duration": [float(subset[:, 0].min()), float(subset[:, 0].max())],
                "range_severity": [float(subset[:, 1].min()), float(subset[:, 1].max())],
            }
            if moea_front is not None:
                rc = compute_range_coverage(subset, moea_front)
                entry["range_coverage"] = rc
            results.append(entry)
        print(f"  N={N}: done ({n_seeds} seeds)")

    # Add MOEA-FIND as a reference point
    if moea_front is not None and len(moea_front) > 0:
        cm = coverage_metrics(moea_front, lb, ub)
        results.append({
            "ensemble_size": len(moea_front),
            "seed": -1,
            "method": "moea_find",
            "n_points": len(moea_front),
            "L2_star": cm.get("L2_star_discrepancy", float("nan")),
            "nn_cv": cm.get("nn_cv", float("nan")),
            "nn_mean": cm.get("nn_mean", float("nan")),
            "range_duration": [float(moea_front[:, 0].min()), float(moea_front[:, 0].max())],
            "range_severity": [float(moea_front[:, 1].min()), float(moea_front[:, 1].max())],
            "range_coverage": {
                "overall_frac": 1.0,
                "mean_duration_frac": 1.0,
                "mean_avg_severity_frac": 1.0,
            },
        })
        print(f"  MOEA-FIND: {len(moea_front)} solutions, L2*={cm.get('L2_star_discrepancy', 0):.4f}")

    return results


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--library", type=Path, required=True,
                   help="Path to characteristics.npz or characteristics.json from stage 03 build_library")
    p.add_argument("--moea-front", type=Path, default=None,
                   help="Path to results.json from a stage 04 run_moea_find run")
    p.add_argument("--n-seeds", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    args = p.parse_args()

    out = stage_output_dir(STAGE, DRIVER)
    (out / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "library": str(args.library),
        "moea_front": str(args.moea_front) if args.moea_front else None,
        "ensemble_sizes": ENSEMBLE_SIZES,
        "n_seeds": args.n_seeds, "seed": args.seed,
    }, indent=2))

    print(f"[02/kirsch_convergence] loading library from {args.library}")
    library = load_kirsch_library(args.library)
    print(f"  loaded {len(library)} traces")

    moea_front = None
    if args.moea_front and args.moea_front.exists():
        moea_front = load_moea_front(args.moea_front)
        if moea_front is not None:
            print(f"  loaded {len(moea_front)} Pareto solutions")
        else:
            print(f"  MOEA-FIND results have 0 solutions; skipping comparison")

    sizes = [s for s in ENSEMBLE_SIZES if s <= len(library)]
    results = run_convergence(library, moea_front, sizes, args.n_seeds, args.seed)
    (out / "convergence.json").write_text(json.dumps(results, indent=2))
    print(f"  outputs: {out}")


if __name__ == "__main__":
    main()
