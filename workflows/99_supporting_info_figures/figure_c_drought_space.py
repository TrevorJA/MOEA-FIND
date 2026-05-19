"""Figure C — drought-space coverage at the recommended (K*, T*).

Pairwise 2-D scatter matrix in the K* drought-metric subspace, with
three superimposed point clouds:

* historical T*-block points (large markers),
* baseline-Kirsch ensemble (small grey dots),
* MOEA-FIND Pareto archive (coloured by archive rank or rank index).

Side panel: bar chart of L2-star discrepancy and NN-CV for each cloud.

Output: ``figures/02_calibration/figure_c_drought_space/figure_c.{png,pdf}``.
"""

from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.discovery.analysis import coverage_metrics  # noqa: E402
from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import COLORS, apply_style  # noqa: E402

STAGE = "02_calibration"
HIST_DRIVER = "t_sensitivity_historical"
KIRSCH_STAGE = "03_kirsch_library"
KIRSCH_DRIVER = "build_library_extended"
DECISION_DRIVER = "decision_matrix"
FIGURE_DRIVER = "figure_c_drought_space"


def _load_recommendation():
    p = stage_output_dir(STAGE, DECISION_DRIVER, create=False) \
        / "pareto_front_KxT.json"
    if not p.exists():
        return None
    return json.loads(p.read_text())


def _hist_block_chars(T: int) -> Optional[pd.DataFrame]:
    p = stage_output_dir(STAGE, HIST_DRIVER, slug=f"T{T:02d}", create=False) \
        / "block_chars_extended.csv"
    return pd.read_csv(p) if p.exists() else None


def _kirsch_chars(T: int, n_traces: int = 10_000, seed: int = 42) -> Optional[pd.DataFrame]:
    p = stage_output_dir(KIRSCH_STAGE, KIRSCH_DRIVER,
                         slug=f"n{n_traces}_t{T}_ssi3-12_s{seed}",
                         create=False) / "characteristics_extended.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    return pd.DataFrame(z["values"], columns=[str(n) for n in z["metric_names"]])


def _load_pareto_metrics(pareto_archive: Optional[Path]) -> Optional[pd.DataFrame]:
    if pareto_archive is None or not pareto_archive.exists():
        return None
    raw = json.loads(pareto_archive.read_text())
    chars = raw.get("pareto_chars") or raw.get("drought_metrics") or []
    if not chars:
        return None
    return pd.DataFrame(chars)


def _coverage(points: np.ndarray, lb: np.ndarray, ub: np.ndarray) -> dict:
    pts = points[np.isfinite(points).all(axis=1)]
    if pts.shape[0] < 2:
        return {"L2_star_discrepancy": float("nan"), "nn_cv": float("nan")}
    return coverage_metrics(pts, lb, ub)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--metrics", type=str, default=None,
                   help="Comma-separated K* metrics. Default: read from "
                        "decision_matrix recommendation.")
    p.add_argument("--T-star", type=int, default=None)
    p.add_argument("--pareto-archive", type=Path, default=None)
    args = p.parse_args()

    apply_style()
    fig_dir = stage_figure_dir(STAGE, FIGURE_DRIVER)

    rec = _load_recommendation()
    if args.metrics:
        cols = [m.strip() for m in args.metrics.split(",") if m.strip()]
    elif rec is not None:
        cols = list(rec["recommended"]["metrics"])
    else:
        cols = ["mean_severity", "mean_magnitude", "time_in_drought_fraction"]
    T_star = (args.T_star if args.T_star is not None
              else (int(rec["recommended"]["T_years"]) if rec is not None else 20))
    print(f"[figure_c] T*={T_star}, metrics={cols}")

    hist_df = _hist_block_chars(T_star)
    kir_df = _kirsch_chars(T_star)
    pareto_df = _load_pareto_metrics(args.pareto_archive)

    if hist_df is None or kir_df is None:
        raise SystemExit(f"[figure_c] missing artifacts at T={T_star}")

    available_cols = [c for c in cols if c in hist_df.columns and c in kir_df.columns]
    if pareto_df is not None:
        available_cols = [c for c in available_cols if c in pareto_df.columns]
    if len(available_cols) < 2:
        raise SystemExit("Fewer than 2 K* metrics available across sources.")

    pairs = list(combinations(available_cols, 2))
    n_pairs = len(pairs)
    n_cols_grid = min(3, n_pairs)
    n_rows_grid = int(np.ceil(n_pairs / n_cols_grid))

    fig = plt.figure(figsize=(3.6 * (n_cols_grid + 1), 3.4 * n_rows_grid))
    gs = fig.add_gridspec(
        n_rows_grid, n_cols_grid + 1,
        width_ratios=[1.0] * n_cols_grid + [0.7],
        hspace=0.35, wspace=0.35,
    )

    # Determine common bounds per metric across sources.
    bounds = {}
    for c in available_cols:
        vals = np.concatenate([
            hist_df[c].dropna().values,
            kir_df[c].dropna().values,
            (pareto_df[c].dropna().values if pareto_df is not None else np.empty(0)),
        ])
        bounds[c] = (float(np.nanmin(vals)), float(np.nanmax(vals)))

    cov_records: List[dict] = []
    for idx, (a, b) in enumerate(pairs):
        r, c = divmod(idx, n_cols_grid)
        ax = fig.add_subplot(gs[r, c])
        ax.scatter(kir_df[a], kir_df[b], s=2, alpha=0.10,
                   color=COLORS["muted"], rasterized=True, label="Kirsch")
        ax.scatter(hist_df[a], hist_df[b], s=24, alpha=0.85,
                   color=COLORS["historical"], edgecolor="white",
                   linewidth=0.4, label="Historical T-blocks")
        if pareto_df is not None:
            ax.scatter(pareto_df[a], pareto_df[b], s=18, alpha=0.85,
                       color=COLORS["highlight"], edgecolor="white",
                       linewidth=0.4, label="MOEA-FIND")
        ax.set_xlabel(a, fontsize=8)
        ax.set_ylabel(b, fontsize=8)
        if idx == 0:
            ax.legend(fontsize=7, loc="best")

    # Coverage panel
    lb = np.array([bounds[c][0] for c in available_cols])
    ub = np.array([bounds[c][1] for c in available_cols])
    ub = np.where(ub > lb, ub, lb + 1e-9)
    sources = {
        "Historical": hist_df[available_cols].dropna().values,
        "Kirsch": kir_df[available_cols].dropna().values,
    }
    if pareto_df is not None:
        sources["MOEA-FIND"] = pareto_df[available_cols].dropna().values
    for label, pts in sources.items():
        cov = _coverage(pts, lb, ub)
        cov_records.append({"source": label,
                            "L2_star_discrepancy": cov["L2_star_discrepancy"],
                            "nn_cv": cov.get("nn_cv", float("nan")),
                            "n_points": int(pts.shape[0])})

    cov_df = pd.DataFrame(cov_records)
    ax_cov = fig.add_subplot(gs[:, -1])
    x = np.arange(len(cov_df))
    width = 0.35
    ax_cov.bar(x - width / 2, cov_df["L2_star_discrepancy"], width,
               label="L2*-discrepancy",
               color=COLORS["empirical"], alpha=0.85)
    ax_cov.bar(x + width / 2, cov_df["nn_cv"], width,
               label="NN CV",
               color=COLORS["highlight"], alpha=0.85)
    ax_cov.set_xticks(x)
    ax_cov.set_xticklabels(cov_df["source"], rotation=45, ha="right")
    ax_cov.set_ylabel("coverage metric (lower = better)")
    ax_cov.legend(fontsize=8)
    ax_cov.set_title("Coverage diagnostics", fontsize=10)

    fig.suptitle(
        f"Figure C — drought-space coverage in K={len(available_cols)} space at T* = {T_star} yr",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(fig_dir / "figure_c.png", dpi=200, bbox_inches="tight")
    fig.savefig(fig_dir / "figure_c.pdf", bbox_inches="tight")
    cov_df.to_csv(fig_dir / "coverage_summary.csv", index=False)
    print(f"[figure_c] wrote {fig_dir / 'figure_c.png'}")
    print(cov_df.to_string(index=False))


if __name__ == "__main__":
    main()
