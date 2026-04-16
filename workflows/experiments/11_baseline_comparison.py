"""Script 11 — Kirsch vs MOEA-FIND baseline comparison (manuscript §6.3, Fig 7).

Produces the manuscript comparison figure showing:
  (a) Kirsch random ensemble scatter in drought space
  (b) MOEA-FIND Pareto front in same space
  (c) Range convergence curves (ensemble size vs coverage fraction)

Reads outputs from:
  - Script 05 (Kirsch library): outputs/exp05_kirsch_library/
  - Script 04 (MOEA-FIND): outputs/exp04_kirsch_single_site/{variant}/
  - Convergence diagnostic: outputs/diag_kirsch_convergence/

Usage:
    python workflows/experiments/11_baseline_comparison.py \\
        --kirsch-library outputs/exp05_kirsch_library/characteristics.npz \\
        --moea-front outputs/exp04_kirsch_single_site/{variant}/results.json \\
        --convergence outputs/diag_kirsch_convergence/convergence.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis import coverage_metrics  # noqa: E402


def load_kirsch(path: Path) -> np.ndarray:
    """Load Kirsch objectives (n, 2)."""
    path = Path(path)
    if path.suffix == ".npz":
        return np.load(path, allow_pickle=True)["objectives"]
    chars = json.loads(path.read_text())
    return np.array([[c["mean_duration"], c["mean_avg_severity"]] for c in chars])


def load_moea(path: Path):
    """Load MOEA-FIND results. Returns (objectives (n,2), result_dict)."""
    r = json.loads(Path(path).read_text())
    dm = r.get("drought_metrics", [])
    if not dm:
        return None, r
    return np.array(dm)[:, :2], r


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--kirsch-library", type=Path, required=True)
    p.add_argument("--moea-front", type=Path, required=True)
    p.add_argument("--convergence", type=Path, default=None,
                   help="convergence.json from diag_kirsch_convergence.py")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / "exp11_baseline_comparison")
    args = p.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    fig_dir = out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # --- Load data ---
    kirsch = load_kirsch(args.kirsch_library)
    moea, moea_result = load_moea(args.moea_front)

    print(f"[11] Kirsch library: {len(kirsch)} traces")
    if moea is not None:
        print(f"[11] MOEA-FIND: {len(moea)} Pareto solutions")
    else:
        print(f"[11] MOEA-FIND: no solutions")
        return

    # Coverage metrics
    lb = kirsch.min(axis=0)
    ub = kirsch.max(axis=0)
    kirsch_cm = coverage_metrics(kirsch, lb, ub)
    moea_cm = coverage_metrics(moea, lb, ub)

    print(f"[11] Kirsch L2*={kirsch_cm.get('L2_star_discrepancy', 0):.4f}, "
          f"NN_CV={kirsch_cm.get('nn_cv', 0):.3f}")
    print(f"[11] MOEA   L2*={moea_cm.get('L2_star_discrepancy', 0):.4f}, "
          f"NN_CV={moea_cm.get('nn_cv', 0):.3f}")

    # --- Figure: scatter comparison ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # Panel (a): Kirsch
    ax1.scatter(kirsch[:, 0], kirsch[:, 1], s=3, alpha=0.15, c="gray",
                label=f"Kirsch ({len(kirsch)} traces)")
    ax1.set_xlabel("Mean duration (months)")
    ax1.set_ylabel("Mean avg. severity")
    ax1.set_title(f"(a) Kirsch random ensemble\n"
                  f"L2*={kirsch_cm.get('L2_star_discrepancy', 0):.3f}")
    ax1.legend(fontsize=8, loc="upper right")

    # Panel (b): MOEA-FIND
    ax2.scatter(moea[:, 0], moea[:, 1], s=15, alpha=0.7, c="#d62728",
                label=f"MOEA-FIND ({len(moea)} solutions)")
    # Show Kirsch cloud as background
    ax2.scatter(kirsch[:, 0], kirsch[:, 1], s=1, alpha=0.05, c="gray", zorder=0)
    ax2.set_xlabel("Mean duration (months)")
    ax2.set_title(f"(b) MOEA-FIND Pareto front\n"
                  f"L2*={moea_cm.get('L2_star_discrepancy', 0):.3f}")
    ax2.legend(fontsize=8, loc="upper right")

    fig.suptitle("Drought Characteristic Space: Kirsch Random vs MOEA-FIND", fontsize=13)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig07_scatter_comparison.pdf", dpi=300)
    plt.close(fig)
    print(f"[11] wrote {fig_dir / 'fig07_scatter_comparison.pdf'}")

    # --- Figure: convergence (if data available) ---
    if args.convergence and args.convergence.exists():
        conv_data = json.loads(args.convergence.read_text())
        kirsch_data = [r for r in conv_data if r.get("method") != "moea_find"]
        moea_data = [r for r in conv_data if r.get("method") == "moea_find"]

        sizes = sorted(set(r["ensemble_size"] for r in kirsch_data))

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Panel (c): Range coverage
        if any("range_coverage" in r for r in kirsch_data):
            medians, lo, hi = [], [], []
            for N in sizes:
                fracs = [r["range_coverage"]["overall_frac"]
                         for r in kirsch_data if r["ensemble_size"] == N]
                medians.append(np.median(fracs))
                lo.append(np.percentile(fracs, 10))
                hi.append(np.percentile(fracs, 90))

            ax1.plot(sizes, medians, "o-", color="#2b6cb0", label="Kirsch random")
            ax1.fill_between(sizes, lo, hi, alpha=0.2, color="#2b6cb0",
                             label="10-90th pct")
            ax1.axhline(1.0, color="gray", linestyle="--", alpha=0.5)
            if moea_data:
                ax1.axvline(moea_data[0]["n_points"], color="#d62728", linestyle=":",
                            label=f"MOEA-FIND n={moea_data[0]['n_points']}")
            ax1.set_xscale("log")
            ax1.set_xlabel("Ensemble size")
            ax1.set_ylabel("Fraction of MOEA-FIND range covered")
            ax1.set_title("(c) Range convergence")
            ax1.legend(fontsize=8)
            ax1.set_ylim(0, 1.05)

        # Panel (d): L2* comparison
        l2_medians = []
        for N in sizes:
            vals = [r["L2_star"] for r in kirsch_data
                    if r["ensemble_size"] == N and np.isfinite(r["L2_star"])]
            l2_medians.append(np.median(vals) if vals else float("nan"))

        ax2.plot(sizes, l2_medians, "o-", color="#2b6cb0", label="Kirsch random")
        if moea_data:
            ax2.axhline(moea_data[0]["L2_star"], color="#d62728", linestyle="--",
                        label=f"MOEA-FIND (n={moea_data[0]['n_points']})")
        ax2.set_xscale("log")
        ax2.set_xlabel("Ensemble size")
        ax2.set_ylabel("L2* discrepancy (lower = better)")
        ax2.set_title("(d) Coverage uniformity")
        ax2.legend(fontsize=8)

        fig.suptitle("Convergence: How Many Kirsch Traces to Match MOEA-FIND?",
                     fontsize=13)
        fig.tight_layout()
        fig.savefig(fig_dir / "fig07_convergence.pdf", dpi=300)
        plt.close(fig)
        print(f"[11] wrote {fig_dir / 'fig07_convergence.pdf'}")

    # --- Save summary ---
    summary = {
        "kirsch_n": len(kirsch),
        "moea_n": len(moea),
        "kirsch_L2_star": kirsch_cm.get("L2_star_discrepancy"),
        "moea_L2_star": moea_cm.get("L2_star_discrepancy"),
        "kirsch_nn_cv": kirsch_cm.get("nn_cv"),
        "moea_nn_cv": moea_cm.get("nn_cv"),
        "kirsch_range_duration": [float(kirsch[:, 0].min()), float(kirsch[:, 0].max())],
        "kirsch_range_severity": [float(kirsch[:, 1].min()), float(kirsch[:, 1].max())],
        "moea_range_duration": [float(moea[:, 0].min()), float(moea[:, 0].max())],
        "moea_range_severity": [float(moea[:, 1].min()), float(moea[:, 1].max())],
    }
    (out / "comparison_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[11] wrote {out / 'comparison_summary.json'}")


if __name__ == "__main__":
    main()
