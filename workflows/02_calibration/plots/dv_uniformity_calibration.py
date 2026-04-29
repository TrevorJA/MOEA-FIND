"""Render the DV-uniformity bootstrap calibration histograms (SI B).

Reads ``outputs/02_calibration/dv_uniformity_calibration/<dv_mode>/{
calibrated_dv_tolerances.json, bootstrap_samples.npz}`` and writes
``figures/02_calibration/dv_uniformity_calibration/<dv_mode>/dv_uniformity_calibration.pdf``.

Plotting-only driver -- never re-runs the bootstrap. Pair with the
compute driver ``workflows/02_calibration/dv_uniformity_calibration.py``.
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

from src.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import apply_style  # noqa: E402

STAGE = "02_calibration"
DRIVER = "dv_uniformity_calibration"


def _make_pdf(summaries, sample_map, output_path):
    names = list(summaries.keys())
    n = len(names)
    fig, axes = plt.subplots(n, 1, figsize=(8, 2.6 * n), constrained_layout=True)
    if n == 1:
        axes = [axes]

    for ax, name in zip(axes, names):
        s = summaries[name]
        x = sample_map[name]
        x = x[np.isfinite(x)]
        ax.hist(x, bins=40, alpha=0.55, color="#c05621",
                label=f"Random U[0,1] bootstrap (n={len(x)})",
                density=True)
        tol = s["tolerance"]
        ax.axvline(tol, color="black", linestyle="--", lw=1,
                   label=f"tol (97.5%) = {tol:.4g}")
        ax.axvline(s["q025"], color="gray", linestyle=":", lw=1,
                   label=f"2.5% = {s['q025']:.4g}")
        ax.axvline(s["q975"], color="gray", linestyle=":", lw=1)
        ax.set_title(f"{name} on random U[0,1] DV draws (N={s['n_dvs']})")
        ax.set_xlabel(name)
        ax.set_ylabel("density")
        ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(
        "DV-uniformity bootstrap calibration (random U[0,1] DV draws)",
        fontsize=12,
    )
    fig.savefig(output_path, format="pdf")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dv-mode", choices=["residual", "index"], required=True,
                   help="Which calibration arm to render.")
    args = p.parse_args()

    apply_style()

    in_dir = stage_output_dir(STAGE, DRIVER, slug=args.dv_mode, create=False)
    json_path = in_dir / "calibrated_dv_tolerances.json"
    npz_path = in_dir / "bootstrap_samples.npz"
    if not json_path.exists() or not npz_path.exists():
        sys.exit(f"missing inputs in {in_dir} -- run the compute driver first")

    payload = json.loads(json_path.read_text())
    # Each value carries {statistic, tolerance, details}; reconstruct the
    # summary dict the original plotter consumed.
    summaries = {entry["statistic"]: entry["details"] for entry in payload.values()}
    npz = np.load(npz_path)
    sample_map = {name: npz[name] for name in summaries}

    fig_dir = stage_figure_dir(STAGE, DRIVER, slug=args.dv_mode)
    out_pdf = fig_dir / "dv_uniformity_calibration.pdf"
    _make_pdf(summaries, sample_map, out_pdf)
    print(f"[plots/02/dv_uniformity_calibration:{args.dv_mode}] {out_pdf}")


if __name__ == "__main__":
    main()
