"""Render the plausibility-constraint bootstrap calibration histograms (SI A).

Reads ``outputs/02_calibration/constraint_calibration/{calibrated_tolerances.json,
bootstrap_samples.npz}`` and writes
``figures/02_calibration/constraint_calibration/constraint_calibration.pdf``.

Plotting-only driver -- never re-runs the bootstrap. Pair with the
compute driver ``workflows/02_calibration/constraint_calibration.py``.
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
from src.plotting.style import apply_style  # noqa: E402

STAGE = "02_calibration"
DRIVER = "constraint_calibration"

STAT_NAMES = (
    "annual_mean",
    "annual_cv",
    "lag1_ac_monthly",
    "non_drought_mean",
    "seasonal_cycle_max_dev",
)


def _make_pdf(summaries, hist_map, kirsch_map, output_path):
    n = len(STAT_NAMES)
    fig, axes = plt.subplots(n, 1, figsize=(8, 2.4 * n), constrained_layout=True)
    if n == 1:
        axes = [axes]

    for ax, name in zip(axes, STAT_NAMES):
        s = summaries[name]
        ref = s["reference"]
        mode = s["mode"]
        h = hist_map[name]
        k = kirsch_map[name]
        if mode == "fractional":
            h_dev = h / ref - 1.0 if ref else h * 0.0
            k_dev = k / ref - 1.0 if ref else k * 0.0
            xlabel = f"{name} fractional deviation"
        else:
            h_dev = h - ref
            k_dev = k - ref
            xlabel = f"{name} absolute deviation"

        h_dev = h_dev[np.isfinite(h_dev)]
        k_dev = k_dev[np.isfinite(k_dev)]

        bins = 40
        ax.hist(h_dev, bins=bins, alpha=0.5, color="#2b6cb0",
                label=f"Historical (n={len(h_dev)})", density=True)
        ax.hist(k_dev, bins=bins, alpha=0.5, color="#c05621",
                label=f"Kirsch (n={len(k_dev)})", density=True)
        tol = s["tolerance"]
        ax.axvline(-tol, color="black", linestyle="--", lw=1)
        ax.axvline(+tol, color="black", linestyle="--", lw=1,
                   label=f"tol = {tol:.3f}")
        ax.set_title(name)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("density")
        ax.legend(fontsize=8, loc="upper right")

    fig.suptitle("Plausibility-constraint bootstrap calibration",
                 fontsize=12)
    fig.savefig(output_path, format="pdf")
    plt.close(fig)


def main():
    apply_style()

    in_dir = stage_output_dir(STAGE, DRIVER, create=False)
    json_path = in_dir / "calibrated_tolerances.json"
    npz_path = in_dir / "bootstrap_samples.npz"
    if not json_path.exists() or not npz_path.exists():
        sys.exit(f"missing inputs in {in_dir} -- run the compute driver first")

    payload = json.loads(json_path.read_text())
    # The JSON has a single "{site}_T{T}" key; pull its details dict.
    block = next(iter(payload.values()))
    summaries = block["details"]

    npz = np.load(npz_path)
    hist_map = {n: npz[f"h_{n}"] for n in STAT_NAMES}
    kirsch_map = {n: npz[f"k_{n}"] for n in STAT_NAMES}

    fig_dir = stage_figure_dir(STAGE, DRIVER)
    out_pdf = fig_dir / "constraint_calibration.pdf"
    _make_pdf(summaries, hist_map, kirsch_map, out_pdf)
    print(f"[plots/02/constraint_calibration] {out_pdf}")


if __name__ == "__main__":
    main()
