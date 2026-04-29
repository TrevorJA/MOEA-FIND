"""Render the Kirsch-library vs MOEA-FIND coverage figure (main Fig 7).

Reads ``outputs/04_moea_find_single_site/baseline_comparison/pooled.npz``
and writes the comparison panels to
``figures/04_moea_find_single_site/baseline_comparison/``.

Plotting only -- never re-runs the comparison.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "baseline_comparison"


def _fig_scatter(kirsch, moea, anti_ideal, fig_path):
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.scatter(kirsch[:, 0], kirsch[:, 1], s=6, alpha=0.35,
               color="#7f7f7f", label=f"Kirsch library (n={len(kirsch)})")
    ax.scatter(moea[:, 0], moea[:, 1], s=14, alpha=0.85,
               color="#2b6cb0", label=f"MOEA-FIND Pareto (n={len(moea)})")
    if anti_ideal.size >= 2:
        ax.scatter(anti_ideal[0], anti_ideal[1], marker="X", s=80,
                   color="#d62728", label="anti-ideal D*")
    ax.set_xlabel("Mean duration (months)")
    ax.set_ylabel("Mean avg. severity")
    ax.set_title("Drought-space coverage: Kirsch library vs MOEA-FIND")
    ax.legend(fontsize=9, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_summary_bars(summary, fig_path):
    metrics = [("L2_star", "L2* discrepancy (lower better)"),
               ("nn_cv",   "NN-CV (lower better)")]
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 4.0))
    for ax, (key, title) in zip(axes, metrics):
        kv = summary.get(f"kirsch_{key}", float("nan"))
        mv = summary.get(f"moea_{key}", float("nan"))
        ax.bar(["Kirsch library", "MOEA-FIND"], [kv, mv],
               color=["#7f7f7f", "#2b6cb0"])
        for i, v in enumerate([kv, mv]):
            if v is not None and np.isfinite(v):
                ax.text(i, v, f"{v:.4g}", ha="center", va="bottom", fontsize=9)
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    in_dir = stage_output_dir(STAGE, DRIVER, create=False)
    pooled_path = in_dir / "pooled.npz"
    summary_path = in_dir / "comparison_summary.json"
    if not pooled_path.exists():
        sys.exit(f"missing {pooled_path} -- run the compute driver first")
    npz = np.load(pooled_path)
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}

    fig_dir = stage_figure_dir(STAGE, DRIVER)
    print(f"[plots/04/baseline_comparison] fig_dir={fig_dir}")
    _fig_scatter(npz["kirsch"], npz["moea"], npz["anti_ideal"],
                 fig_dir / "fig07_library_vs_moea.pdf")
    if summary:
        _fig_summary_bars(summary, fig_dir / "fig07_coverage_bars.pdf")


if __name__ == "__main__":
    main()
