"""Render per-preset 3D scatter figures from cached metric_blocks output.

Reads ``outputs/02_calibration/metric_blocks/{block_chars.csv,
full_hist_chars.pkl}`` and writes one 3D scatter per preset under
``figures/02_calibration/metric_blocks/``.
"""

from __future__ import annotations

import pickle
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.metrics.drought_metrics import PRESETS, resolve_metric_set  # noqa: E402
from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "02_calibration"
DRIVER = "metric_blocks"


def _plot_metric_set_3d(chars_df, preset_name, metric_set, fig_path, *,
                        historical_aggregate=None):
    if len(metric_set) < 3:
        return
    pts = np.array([[m.extract(row) for m in metric_set]
                    for row in chars_df.to_dict("records")])
    fig = plt.figure(figsize=(7.0, 6.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2],
               c=np.arange(len(pts)), cmap="viridis",
               s=24, alpha=0.85, edgecolors="black", linewidths=0.3)
    if historical_aggregate is not None:
        ha = np.array([m.extract(historical_aggregate) for m in metric_set])
        ax.scatter([ha[0]], [ha[1]], [ha[2]],
                   marker="*", s=220, c="red",
                   edgecolors="black", linewidths=0.6,
                   label="full historical record")
        ax.legend(loc="upper left", fontsize=8)
    m0, m1, m2 = metric_set[0], metric_set[1], metric_set[2]
    ax.set_xlabel(f"{m0.label}\n({m0.units})", fontsize=9)
    ax.set_ylabel(f"{m1.label}\n({m1.units})", fontsize=9)
    ax.set_zlabel(f"{m2.label}\n({m2.units})", fontsize=9)
    ax.set_title(
        f"Preset: {preset_name}\n"
        f"{len(pts)} historical T-year blocks (colour = chronological order)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(fig_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(fig_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)


def main():
    in_dir = stage_output_dir(STAGE, DRIVER, create=False)
    chars_path = in_dir / "block_chars.csv"
    pkl_path = in_dir / "full_hist_chars.pkl"
    if not chars_path.exists():
        sys.exit(f"missing {chars_path} -- run the compute driver first")
    chars_df = pd.read_csv(chars_path)
    full_hist_chars = pickle.loads(pkl_path.read_bytes()) if pkl_path.exists() else None

    fig_dir = stage_figure_dir(STAGE, DRIVER)
    print(f"[plots/02/metric_blocks] fig_dir={fig_dir}")
    for preset_name in PRESETS:
        ms = resolve_metric_set(preset_name)
        fig_path = fig_dir / f"fig_{preset_name}_3d"
        _plot_metric_set_3d(chars_df, preset_name, ms, fig_path,
                            historical_aggregate=full_hist_chars)


if __name__ == "__main__":
    main()
