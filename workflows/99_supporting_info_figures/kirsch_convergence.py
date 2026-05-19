"""Render Kirsch ensemble convergence figures from cached compute outputs.

Reads ``outputs/02_calibration/kirsch_convergence/convergence.json`` and
writes range / coverage convergence panels to
``figures/02_calibration/kirsch_convergence/``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "02_calibration"
DRIVER = "kirsch_convergence"


def main():
    in_dir = stage_output_dir(STAGE, DRIVER, create=False)
    json_path = in_dir / "convergence.json"
    if not json_path.exists():
        sys.exit(f"missing {json_path} -- run the compute driver first")
    results = json.loads(json_path.read_text())

    fig_dir = stage_figure_dir(STAGE, DRIVER)
    kirsch = [r for r in results if r.get("method") != "moea_find"]
    moea = [r for r in results if r.get("method") == "moea_find"]
    sizes = sorted(set(r["ensemble_size"] for r in kirsch))

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
    print(f"[plots/02/kirsch_convergence] fig_dir={fig_dir}")


if __name__ == "__main__":
    main()
