"""Render Cannonsville Pareto + plausibility figures (main Figs 5, 6).

Reads ``outputs/04_moea_find_single_site/run_moea_find/<slug>/`` and
writes the drought-space scatter, plausibility panels, and trace
diagnostics into ``figures/04_moea_find_single_site/run_moea_find/<slug>/``.

Plotting only -- never re-runs MOEA-FIND. Slug must match a completed
compute run.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "run_moea_find"


def _fig_drought_space(metrics, anti_ideal, hist_chars, fig_path):
    fig, ax = plt.subplots(figsize=(7.0, 5.5))
    if hist_chars is not None and hist_chars.shape[0] > 0:
        ax.scatter(hist_chars[:, 0], hist_chars[:, 1],
                   s=18, marker="x", color="#7f7f7f",
                   label=f"historical T-blocks (n={len(hist_chars)})")
    ax.scatter(metrics[:, 0], metrics[:, 1], s=12, alpha=0.85,
               color="#2b6cb0", label=f"MOEA-FIND Pareto (n={len(metrics)})")
    if anti_ideal.size >= 2:
        ax.scatter(anti_ideal[0], anti_ideal[1], marker="X", s=80,
                   color="#d62728", label=r"anti-ideal $D^*$")
    ax.set_xlabel("Mean duration (months)")
    ax.set_ylabel("Mean avg. severity")
    ax.set_title("Cannonsville drought-space coverage")
    ax.legend(fontsize=9, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_drought_space_3d(metrics, anti_ideal, hist_chars, objective_keys, fig_path):
    if metrics.shape[1] < 3:
        return
    fig = plt.figure(figsize=(7.0, 6.0))
    ax = fig.add_subplot(111, projection="3d")
    if hist_chars is not None and hist_chars.shape[0] > 0 and hist_chars.shape[1] >= 3:
        ax.scatter(hist_chars[:, 0], hist_chars[:, 1], hist_chars[:, 2],
                   s=14, marker="x", color="#7f7f7f", label="historical")
    ax.scatter(metrics[:, 0], metrics[:, 1], metrics[:, 2],
               s=8, alpha=0.7, color="#2b6cb0", label="MOEA-FIND")
    if anti_ideal.size >= 3:
        ax.scatter([anti_ideal[0]], [anti_ideal[1]], [anti_ideal[2]],
                   marker="X", s=80, color="#d62728", label=r"$D^*$")
    if len(objective_keys) >= 3:
        ax.set_xlabel(str(objective_keys[0]))
        ax.set_ylabel(str(objective_keys[1]))
        ax.set_zlabel(str(objective_keys[2]))
    ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True,
                   help="Variant slug under outputs/04_moea_find_single_site/run_moea_find/")
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    results_path = in_dir / "results.json"
    if not results_path.exists():
        sys.exit(f"missing {results_path} -- run the compute driver first")
    result = json.loads(results_path.read_text())
    metrics = np.array(result.get("drought_metrics", []), dtype=float)
    anti_ideal = np.asarray(result.get("anti_ideal", []), dtype=float)
    objective_keys = result.get("objective_keys", [])

    hist_chars = None
    hist_npz = in_dir / "historical_block_chars.npz"
    if hist_npz.exists():
        hist_chars = np.load(hist_npz)["chars"]

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    print(f"[plots/04/run_moea_find] fig_dir={fig_dir}")
    if metrics.size == 0:
        print("  no Pareto solutions; nothing to render.")
        return
    _fig_drought_space(metrics, anti_ideal, hist_chars,
                       fig_dir / "fig05_drought_space_2d.pdf")
    _fig_drought_space_3d(metrics, anti_ideal, hist_chars, objective_keys,
                          fig_dir / "fig06_drought_space_3d.pdf")


if __name__ == "__main__":
    main()
