"""Render wrapper-geometry mapping figures from cached compute outputs.

Reads ``outputs/02_calibration/wrapper_geometry/sweep.npz`` and writes
mapping smoothness + cartoon panels to
``figures/02_calibration/wrapper_geometry/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "02_calibration"
DRIVER = "wrapper_geometry"

MODE_COLORS = {"index": "#2b6cb0", "residual": "#c05621"}


def _fig_mapping_smoothness(npz, fig_path):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), sharex=True)
    for mode_idx, mode in enumerate(("index", "residual")):
        p_grid = npz[f"{mode}_p_grid"]
        flow_at_k = npz[f"{mode}_flow_at_k"]
        annual_totals = npz[f"{mode}_annual_totals"]
        axes[0].plot(p_grid, flow_at_k, color=MODE_COLORS[mode], lw=1.6, label=mode)
        med = np.median(annual_totals, axis=1)
        lo = np.percentile(annual_totals, 10, axis=1)
        hi = np.percentile(annual_totals, 90, axis=1)
        axes[1].plot(p_grid, med, color=MODE_COLORS[mode], lw=1.6, label=f"{mode} median")
        axes[1].fill_between(p_grid, lo, hi, color=MODE_COLORS[mode], alpha=0.18)
    axes[0].set_xlabel("DV value (sweep coordinate)")
    axes[0].set_ylabel("Flow at the swept month (cfs)")
    axes[0].set_title("(a) Per-month flow response")
    axes[1].set_xlabel("DV value (sweep coordinate)")
    axes[1].set_ylabel("Annual total (cfs); 10/50/90 envelope")
    axes[1].set_title("(b) Annual total response")
    for ax in axes:
        ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_rep_profiles(npz, fig_path):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), sharex=True, sharey=True)
    months = np.arange(12)
    rep_dvs = npz["rep_dvs"]  # shared across modes
    for ax, mode in zip(axes, ("index", "residual")):
        rep_profiles = npz[f"{mode}_rep_profiles"]
        for i, dv in enumerate(rep_dvs):
            ax.plot(months, rep_profiles[i], lw=1.5, label=f"DV={float(dv):.2f}")
        ax.set_xticks(months)
        ax.set_xlabel("Month")
        ax.set_ylabel("Median monthly flow (cfs)")
        ax.set_title(f"{mode} wrapper")
        ax.legend(fontsize=8, framealpha=0.9)
    fig.suptitle("Representative seasonal profiles per DV value", fontsize=11)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    in_dir = stage_output_dir(STAGE, DRIVER, create=False)
    sweep_path = in_dir / "sweep.npz"
    if not sweep_path.exists():
        sys.exit(f"missing {sweep_path} -- run the compute driver first")
    npz = np.load(sweep_path)

    fig_dir = stage_figure_dir(STAGE, DRIVER)
    print(f"[plots/02/wrapper_geometry] fig_dir={fig_dir}")
    _fig_mapping_smoothness(npz, fig_dir / "fig_mapping_smoothness.pdf")
    _fig_rep_profiles(npz, fig_dir / "fig_mapping_cartoon.pdf")


if __name__ == "__main__":
    main()
