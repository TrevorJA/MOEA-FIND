"""Diagnostic: Kirsch ensemble convergence analysis.

Measures how drought characteristic range and coverage metrics converge as
Kirsch ensemble size increases. Subsamples a pre-generated library at varying
sizes and compares to a MOEA-FIND Pareto front (if provided).

This analysis supports manuscript H1 (range convergence) and H2 (coverage
efficiency) — see manuscript/experiment_plan_baseline_comparison.md.

Outputs under outputs/diag_kirsch_convergence/:
    - convergence.json (metrics at each ensemble size × seed)
    - convergence_range.pdf (range coverage vs ensemble size)
    - convergence_coverage.pdf (L2* and NN_CV vs ensemble size)

Usage:
    python workflows/02_calibration/kirsch_convergence.py \\
        --library outputs/exp05_kirsch_library/characteristics.npz \\
        --moea-front outputs/exp04_kirsch_single_site/{variant}/results.json
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

from src.analysis import coverage_metrics  # noqa: E402

OUTPUT_SLUG = "diag_kirsch_convergence"

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


def plot_convergence(results: List[Dict], moea_front: Optional[np.ndarray],
                     output_dir: Path):
    """Generate convergence diagnostic plots."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Separate Kirsch and MOEA results
    kirsch = [r for r in results if r.get("method") != "moea_find"]
    moea = [r for r in results if r.get("method") == "moea_find"]

    sizes = sorted(set(r["ensemble_size"] for r in kirsch))

    # --- Range coverage ---
    if any("range_coverage" in r for r in kirsch):
        fig, ax = plt.subplots(figsize=(7, 4))
        medians, lo, hi = [], [], []
        for N in sizes:
            fracs = [r["range_coverage"]["overall_frac"]
                     for r in kirsch if r["ensemble_size"] == N]
            medians.append(np.median(fracs))
            lo.append(np.percentile(fracs, 10))
            hi.append(np.percentile(fracs, 90))

        ax.plot(sizes, medians, "o-", color="#2b6cb0", label="Kirsch random (median)")
        ax.fill_between(sizes, lo, hi, alpha=0.2, color="#2b6cb0", label="10-90th pct")
        ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5, label="100% coverage")

        if moea:
            ax.axvline(moea[0]["n_points"], color="#d62728", linestyle=":",
                       label=f"MOEA-FIND ({moea[0]['n_points']} solutions)")

        ax.set_xscale("log")
        ax.set_xlabel("Ensemble size")
        ax.set_ylabel("Fraction of MOEA-FIND range covered")
        ax.set_title("Range convergence: Kirsch vs MOEA-FIND")
        ax.legend(fontsize=8)
        ax.set_ylim(0, 1.05)
        fig.tight_layout()
        fig.savefig(fig_dir / "convergence_range.pdf", dpi=300)
        plt.close(fig)
        print(f"[conv] wrote {fig_dir / 'convergence_range.pdf'}")

    # --- Coverage metrics (L2*, NN_CV) ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    for metric, ax, label in [
        ("L2_star", ax1, "L2* discrepancy"),
        ("nn_cv", ax2, "NN coefficient of variation"),
    ]:
        medians, lo, hi = [], [], []
        for N in sizes:
            vals = [r[metric] for r in kirsch
                    if r["ensemble_size"] == N and np.isfinite(r[metric])]
            if vals:
                medians.append(np.median(vals))
                lo.append(np.percentile(vals, 10))
                hi.append(np.percentile(vals, 90))
            else:
                medians.append(float("nan"))
                lo.append(float("nan"))
                hi.append(float("nan"))

        ax.plot(sizes, medians, "o-", color="#2b6cb0", label="Kirsch random")
        ax.fill_between(sizes, lo, hi, alpha=0.2, color="#2b6cb0")

        if moea and np.isfinite(moea[0][metric]):
            ax.axhline(moea[0][metric], color="#d62728", linestyle="--",
                       label=f"MOEA-FIND ({moea[0]['n_points']})")

        ax.set_xscale("log")
        ax.set_xlabel("Ensemble size")
        ax.set_ylabel(label)
        ax.legend(fontsize=8)

    ax1.set_title("Uniformity: L2* discrepancy (lower = better)")
    ax2.set_title("Regularity: NN_CV (lower = more uniform spacing)")
    fig.tight_layout()
    fig.savefig(fig_dir / "convergence_coverage.pdf", dpi=300)
    plt.close(fig)
    print(f"[conv] wrote {fig_dir / 'convergence_coverage.pdf'}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--library", type=Path, required=True,
                   help="Path to characteristics.npz or characteristics.json from script 05")
    p.add_argument("--moea-front", type=Path, default=None,
                   help="Path to results.json from script 04 MOEA-FIND run")
    p.add_argument("--n-seeds", type=int, default=N_SEEDS)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    print(f"[conv] loading Kirsch library from {args.library} ...")
    library = load_kirsch_library(args.library)
    print(f"[conv] loaded {len(library)} traces")

    moea_front = None
    if args.moea_front and args.moea_front.exists():
        print(f"[conv] loading MOEA-FIND front from {args.moea_front} ...")
        moea_front = load_moea_front(args.moea_front)
        if moea_front is not None:
            print(f"[conv] loaded {len(moea_front)} Pareto solutions")
        else:
            print(f"[conv] MOEA-FIND results have 0 solutions; skipping comparison")

    # Filter ensemble sizes to library size
    sizes = [s for s in ENSEMBLE_SIZES if s <= len(library)]

    print(f"[conv] running convergence analysis: sizes={sizes}, seeds={args.n_seeds}")
    results = run_convergence(library, moea_front, sizes, args.n_seeds, args.seed)

    (out / "convergence.json").write_text(json.dumps(results, indent=2))
    print(f"[conv] wrote {out / 'convergence.json'}")

    print(f"[conv] generating plots ...")
    plot_convergence(results, moea_front, out)


if __name__ == "__main__":
    main()
