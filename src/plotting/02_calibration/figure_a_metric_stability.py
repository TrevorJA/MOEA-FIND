"""Figure A — drought-metric distributional stability across T.

Multipanel grid: rows = each metric in the recommended K* set + 4-5
runner-up metrics (one per concept axis); columns = the surviving T
values from the joint K x T sweep. Each cell overlays:

* the historical T-block distribution as a violin (Stage 1 output),
* the baseline-Kirsch ensemble as a translucent ridge / KDE histogram
  (Stage 2 output),
* MOEA-FIND Pareto values at T* as scatter dots (Stage 4 output;
  optional — falls back to Hist+Kirsch only if absent).

Drives the manuscript narrative for the metric-and-T joint choice:
the reader can see at a glance which metrics are stable across T and
which axes the Kirsch generator faithfully reproduces.

Output: ``figures/02_calibration/figure_a_metric_stability/figure_a.{png,pdf}``.
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

from src.metrics.extended import CONCEPT_MAP  # noqa: E402
from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import COLORS, apply_style  # noqa: E402

STAGE = "02_calibration"
HIST_DRIVER = "t_sensitivity_historical"
KIRSCH_STAGE = "03_kirsch_library"
KIRSCH_DRIVER = "build_library_extended"
DECISION_DRIVER = "decision_matrix"
FIGURE_DRIVER = "figure_a_metric_stability"


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
    """Load Pareto archive metric values; archive must contain a
    ``pareto_chars`` field shaped as a list of dicts (per-solution
    drought characteristics keyed by metric name)."""
    if pareto_archive is None or not pareto_archive.exists():
        return None
    raw = json.loads(pareto_archive.read_text())
    chars = raw.get("pareto_chars") or raw.get("drought_metrics") or []
    if not chars:
        return None
    return pd.DataFrame(chars)


def _metrics_to_plot(
    decision_path: Path,
    n_runners_up: int = 4,
) -> Tuple[List[str], int]:
    """Pick K* metrics + N runner-ups (one per orthogonal concept).

    Falls back to a default Tier-A trio if ``pareto_front_KxT.json``
    is absent (allows Figure A to be rendered before Stage 3 finishes
    for a preview).
    """
    if decision_path.exists():
        rec = json.loads(decision_path.read_text())
        chosen = list(rec["recommended"]["metrics"])
        T_star = int(rec["recommended"]["T_years"])
    else:
        chosen = ["mean_severity", "mean_magnitude", "time_in_drought_fraction"]
        T_star = 20
    chosen_concepts = {CONCEPT_MAP.get(m, "?") for m in chosen}
    runner_pool = [m for m, c in CONCEPT_MAP.items()
                   if c not in chosen_concepts]
    runners: List[str] = []
    seen_concepts = set()
    for m in runner_pool:
        c = CONCEPT_MAP.get(m, "?")
        if c in seen_concepts:
            continue
        runners.append(m)
        seen_concepts.add(c)
        if len(runners) >= n_runners_up:
            break
    return chosen + runners, T_star


def _draw_cell(
    ax,
    metric: str,
    T: int,
    hist_vals: np.ndarray,
    kir_vals: np.ndarray,
    pareto_vals: Optional[np.ndarray],
    is_chosen: bool,
):
    parts: List[float] = []
    if hist_vals.size:
        ax.violinplot(
            hist_vals, positions=[0.0], widths=0.7, vert=True,
            showextrema=False, showmedians=True,
        )
        parts.extend(hist_vals.tolist())
    if kir_vals.size:
        # Kirsch as horizontal jittered scatter behind the violin.
        rng = np.random.default_rng(0)
        jitter = rng.uniform(0.45, 1.05, size=kir_vals.size)
        ax.scatter(jitter, kir_vals, s=2, alpha=0.05,
                   color=COLORS["kde"], rasterized=True, label="Kirsch")
        parts.extend(kir_vals.tolist())
    if pareto_vals is not None and pareto_vals.size:
        rng = np.random.default_rng(1)
        jitter = rng.uniform(1.4, 1.7, size=pareto_vals.size)
        ax.scatter(jitter, pareto_vals, s=10, alpha=0.7,
                   color=COLORS["highlight"], edgecolor="none",
                   label="MOEA-FIND")
        parts.extend(pareto_vals.tolist())

    ax.set_xlim(-0.5, 2.0)
    ax.set_xticks([])
    if not is_chosen:
        for spine in ("left", "bottom"):
            ax.spines[spine].set_alpha(0.4)
    ax.tick_params(axis="y", labelsize=8)
    if parts:
        lo, hi = float(np.percentile(parts, 1)), float(np.percentile(parts, 99))
        if hi > lo:
            margin = 0.05 * (hi - lo)
            ax.set_ylim(lo - margin, hi + margin)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-grid", type=int, nargs="+", default=[5, 10, 20, 30])
    p.add_argument("--pareto-archive", type=Path, default=None,
                   help="Stage-4 results.json with `pareto_chars` (optional).")
    p.add_argument("--n-runners-up", type=int, default=4)
    args = p.parse_args()

    apply_style()
    fig_dir = stage_figure_dir(STAGE, FIGURE_DRIVER)
    print(f"[figure_a] → {fig_dir}")

    decision_path = stage_output_dir(STAGE, DECISION_DRIVER, create=False) \
        / "pareto_front_KxT.json"
    metrics, T_star = _metrics_to_plot(decision_path, args.n_runners_up)
    chosen_set = set(metrics[: max(0, len(metrics) - args.n_runners_up)])
    print(f"[figure_a] T*={T_star}, metrics={metrics}")

    hist_T: Dict[int, pd.DataFrame] = {}
    kirsch_T: Dict[int, pd.DataFrame] = {}
    for T in args.T_grid:
        h = _hist_block_chars(T)
        k = _kirsch_chars(T)
        if h is None:
            print(f"[warn] T={T}: missing historical chars")
            continue
        if k is None:
            print(f"[warn] T={T}: missing Kirsch chars")
        hist_T[T] = h
        if k is not None:
            kirsch_T[T] = k

    pareto_df = _load_pareto_metrics(args.pareto_archive)
    if pareto_df is not None:
        print(f"[figure_a] pareto archive: {len(pareto_df)} solutions")

    n_rows = len(metrics)
    n_cols = len(args.T_grid)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(2.4 * n_cols, 1.6 * n_rows),
        squeeze=False,
        sharex="col",
    )
    for i, metric in enumerate(metrics):
        for j, T in enumerate(args.T_grid):
            ax = axes[i, j]
            hist_vals = hist_T.get(T, pd.DataFrame()).get(metric, pd.Series([])).dropna().values
            kir_vals = kirsch_T.get(T, pd.DataFrame()).get(metric, pd.Series([])).dropna().values
            pareto_vals = (
                pareto_df[metric].dropna().values
                if (pareto_df is not None and metric in pareto_df.columns and T == T_star)
                else None
            )
            _draw_cell(ax, metric, T, hist_vals, kir_vals, pareto_vals,
                       is_chosen=metric in chosen_set)
            if i == 0:
                ax.set_title(f"T = {T} yr")
            if j == 0:
                tag = "★ " if metric in chosen_set else ""
                ax.set_ylabel(f"{tag}{metric}\n[{CONCEPT_MAP.get(metric,'?')}]",
                              fontsize=8)

    fig.suptitle(
        "Figure A — historical T-blocks, Kirsch ensembles, "
        f"and MOEA-FIND Pareto at T* = {T_star} yr.",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(fig_dir / "figure_a.png", dpi=200)
    fig.savefig(fig_dir / "figure_a.pdf")
    print(f"[figure_a] wrote figure_a.{{png,pdf}}")


if __name__ == "__main__":
    main()
