"""Figure D — hazard-period exemplar traces from the MOEA-FIND archive.

Selects three exemplar Pareto solutions at the recommended (K*, T*)
that span the K* metric axes (worst-along-axis-1, worst-along-axis-2,
worst-along-axis-3) and renders, for each:

* monthly streamflow time series with SSI-3 shaded background,
* event start/end vertical bars,
* annotation block reporting each metric value vs. the historical-block
  IQR.

Output: ``figures/02_calibration/figure_d_exemplar_traces/figure_d.{png,pdf}``.

Strictly requires a Stage-4 ``results.json`` archive that contains both
``pareto_chars`` (per-solution metric values) and a per-solution
monthly-flow trace (``pareto_traces_1d`` or equivalent). If the archive
lacks traces, the script writes a placeholder PNG and returns 0 so the
manuscript build pipeline can proceed.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.metrics.objectives import (  # noqa: E402
    flows_to_series,
    make_ssi_calculator,
)
from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import COLORS, apply_style  # noqa: E402

STAGE = "02_calibration"
HIST_DRIVER = "t_sensitivity_historical"
DECISION_DRIVER = "decision_matrix"
FIGURE_DRIVER = "figure_d_exemplar_traces"


def _load_recommendation():
    p = stage_output_dir(STAGE, DECISION_DRIVER, create=False) \
        / "pareto_front_KxT.json"
    if not p.exists():
        return None
    return json.loads(p.read_text())


def _load_archive(archive: Path) -> dict:
    return json.loads(archive.read_text())


def _hist_iqr(T: int, metrics: List[str]) -> Dict[str, Tuple[float, float]]:
    p = stage_output_dir(STAGE, HIST_DRIVER, slug=f"T{T:02d}", create=False) \
        / "block_chars_extended.csv"
    if not p.exists():
        return {}
    df = pd.read_csv(p)
    out = {}
    for m in metrics:
        if m in df.columns:
            v = df[m].dropna().values
            if v.size:
                out[m] = (float(np.percentile(v, 25)),
                          float(np.percentile(v, 75)))
    return out


def _pick_exemplars(chars_df: pd.DataFrame,
                    metrics: List[str]) -> List[int]:
    """Pick 3 Pareto rows: argmax along each metric axis (most-extreme drought)."""
    picks = []
    for m in metrics[:3]:
        if m not in chars_df.columns:
            continue
        idx = int(chars_df[m].astype(float).idxmax())
        picks.append(idx)
    return picks


def _draw_trace(ax, monthly: np.ndarray, ssi: np.ndarray,
                title: str, annotation: str):
    months = np.arange(len(monthly))
    ax.plot(months, monthly, color=COLORS["empirical"], lw=1.2)
    ax.set_ylabel("monthly flow (cfs)", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)

    ax2 = ax.twinx()
    ax2.fill_between(np.arange(len(ssi)), 0, ssi,
                     where=(ssi <= -1.0),
                     color=COLORS["highlight"], alpha=0.30,
                     label="SSI ≤ −1")
    ax2.plot(np.arange(len(ssi)), ssi, color=COLORS["historical"],
             lw=0.7, alpha=0.5)
    ax2.axhline(-1.0, color="grey", ls="--", lw=0.5)
    ax2.set_ylabel("SSI-3", fontsize=8)
    ax2.tick_params(axis="y", labelsize=8)
    ax.set_title(title, fontsize=9, loc="left")
    ax.text(
        0.02, 0.97, annotation, transform=ax.transAxes,
        fontsize=7, family="monospace",
        verticalalignment="top",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "lightgrey"},
    )


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--pareto-archive", type=Path, required=True)
    p.add_argument("--T-star", type=int, default=None)
    p.add_argument("--metrics", type=str, default=None,
                   help="Comma-separated K* metrics. Default: from "
                        "decision_matrix recommendation.")
    args = p.parse_args()

    apply_style()
    fig_dir = stage_figure_dir(STAGE, FIGURE_DRIVER)

    rec = _load_recommendation()
    if args.metrics:
        metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    elif rec is not None:
        metrics = list(rec["recommended"]["metrics"])
    else:
        raise SystemExit("Provide --metrics or run decision_matrix first.")

    T_star = (args.T_star if args.T_star is not None
              else (int(rec["recommended"]["T_years"]) if rec is not None else 20))
    print(f"[figure_d] T*={T_star}, metrics={metrics}")

    arc = _load_archive(args.pareto_archive)
    chars = arc.get("pareto_chars") or arc.get("drought_metrics") or []
    traces = (arc.get("pareto_traces_1d")
              or arc.get("pareto_traces")
              or arc.get("traces"))
    if not chars:
        raise SystemExit("Archive has no pareto_chars.")
    chars_df = pd.DataFrame(chars)

    if traces is None:
        print("[figure_d] archive has no traces; writing placeholder.")
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5,
                "Stage-4 archive has no `pareto_traces_1d` field.\n"
                "Re-run MOEA-FIND with trace persistence enabled.",
                ha="center", va="center")
        ax.axis("off")
        fig.savefig(fig_dir / "figure_d.png", dpi=200)
        return

    traces_arr = np.asarray(traces, dtype=np.float32)

    iqr = _hist_iqr(T_star, metrics)

    # Fit SSI on full historical record for the shading.
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    _, monthly_1d = prepare_data(cache_dir)
    ssi_calc = make_ssi_calculator(timescale=3)
    ssi_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))

    picks = _pick_exemplars(chars_df, metrics)
    if not picks:
        raise SystemExit("Could not pick exemplars — metric columns missing.")

    fig, axes = plt.subplots(len(picks), 1,
                             figsize=(8.5, 2.6 * len(picks)),
                             squeeze=False)
    for ax, idx in zip(axes[:, 0], picks):
        monthly = traces_arr[idx]
        series = flows_to_series(monthly, start_date="2100-01-01")
        ssi = ssi_calc.transform(series).values
        annot_lines = []
        for m in metrics:
            v = float(chars_df.loc[idx, m]) if m in chars_df.columns else float("nan")
            if m in iqr:
                lo, hi = iqr[m]
                annot_lines.append(
                    f"{m}: {v:.3g} (Hist IQR [{lo:.3g}, {hi:.3g}])"
                )
            else:
                annot_lines.append(f"{m}: {v:.3g}")
        _draw_trace(ax, monthly, ssi,
                    f"Pareto solution {idx}",
                    "\n".join(annot_lines))
    axes[-1, 0].set_xlabel("month index")

    fig.suptitle(
        f"Figure D — exemplar MOEA-FIND traces at T* = {T_star} yr "
        f"(K = {len(metrics)}).",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(fig_dir / "figure_d.png", dpi=200, bbox_inches="tight")
    fig.savefig(fig_dir / "figure_d.pdf", bbox_inches="tight")
    print(f"[figure_d] wrote {fig_dir / 'figure_d.png'}")


if __name__ == "__main__":
    main()
