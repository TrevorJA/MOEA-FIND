"""Figure B — correlation structure preservation (Hist | Kirsch | MOEA).

Three Spearman clustermaps side-by-side at the recommended T*:
historical T-blocks, baseline Kirsch ensemble, MOEA-FIND Pareto. The
metric ordering and the dendrogram leaf order are locked by the
historical-blocks dendrogram so the three heatmaps are directly
comparable. Demonstrates that MOEA-FIND extends but does not distort
the inter-metric structure of the drought hazard space.

Output: ``figures/02_calibration/figure_b_corr_clustermaps/figure_b.{png,pdf}``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import apply_style  # noqa: E402

STAGE = "02_calibration"
HIST_DRIVER = "t_sensitivity_historical"
KIRSCH_STAGE = "03_kirsch_library"
KIRSCH_DRIVER = "build_library_extended"
DECISION_DRIVER = "decision_matrix"
FIGURE_DRIVER = "figure_b_corr_clustermaps"


def _hist_block_chars(T: int) -> Optional[pd.DataFrame]:
    p = stage_output_dir(STAGE, HIST_DRIVER, slug=f"T{T:02d}", create=False) \
        / "block_chars_extended.csv"
    return pd.read_csv(p) if p.exists() else None


def _hist_spread(T: int) -> Optional[pd.DataFrame]:
    p = stage_output_dir(STAGE, HIST_DRIVER, slug=f"T{T:02d}", create=False) \
        / "per_metric_spread.csv"
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


def _load_T_star() -> int:
    p = stage_output_dir(STAGE, DECISION_DRIVER, create=False) \
        / "pareto_front_KxT.json"
    if p.exists():
        return int(json.loads(p.read_text())["recommended"]["T_years"])
    return 20


def _spearman(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    sub = df[cols].astype(float)
    return sub.corr(method="spearman")


def _leaf_order(rho: pd.DataFrame) -> List[str]:
    if rho.shape[0] < 2:
        return list(rho.columns)
    dist = 1.0 - rho.abs().values
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0
    Z = linkage(squareform(np.clip(dist, 0, None), checks=False), method="average")
    order = leaves_list(Z)
    return [rho.columns[i] for i in order]


def _draw_heatmap(ax, rho: pd.DataFrame, order: List[str], title: str):
    arr = rho.loc[order, order].values
    im = ax.imshow(arr, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_yticklabels(order, fontsize=6)
    ax.set_title(title, fontsize=10)
    return im


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-star", type=int, default=None,
                   help="Override T* (default: read from decision_matrix output).")
    p.add_argument("--pareto-archive", type=Path, default=None,
                   help="Stage-4 results.json with `pareto_chars` (optional).")
    args = p.parse_args()

    apply_style()
    fig_dir = stage_figure_dir(STAGE, FIGURE_DRIVER)
    T_star = args.T_star if args.T_star is not None else _load_T_star()
    print(f"[figure_b] T*={T_star}")

    hist_df = _hist_block_chars(T_star)
    spread = _hist_spread(T_star)
    kir_df = _kirsch_chars(T_star)
    pareto_df = _load_pareto_metrics(args.pareto_archive)

    if hist_df is None or spread is None:
        raise SystemExit(f"[figure_b] missing Stage-1 artifacts at T={T_star}")

    surviving = spread[spread["passes_screen"]]["metric"].tolist()
    cols = [m for m in surviving if m in hist_df.columns]
    if kir_df is not None:
        cols = [m for m in cols if m in kir_df.columns]
    if not cols:
        raise SystemExit("No common metrics across Hist/Kirsch.")

    rho_hist = _spearman(hist_df, cols)
    order = _leaf_order(rho_hist)
    rho_kir = _spearman(kir_df, cols) if kir_df is not None else None
    rho_par = (_spearman(pareto_df, [c for c in cols if c in pareto_df.columns])
               if pareto_df is not None else None)

    panels = [("Historical T-blocks", rho_hist)]
    if rho_kir is not None:
        panels.append(("Baseline Kirsch ensemble", rho_kir))
    if rho_par is not None:
        # Reindex to the historical column set so panels share an axis
        # — missing columns will appear as NaN cells (light grey).
        rho_par = rho_par.reindex(index=cols, columns=cols)
        panels.append(("MOEA-FIND Pareto", rho_par))

    fig, axes = plt.subplots(1, len(panels),
                             figsize=(4.0 * len(panels), 4.6),
                             squeeze=False)
    for ax, (title, rho) in zip(axes[0], panels):
        rho2 = rho.copy()
        rho2.fillna(0.0, inplace=True)
        im = _draw_heatmap(ax, rho2, order, title)
    fig.colorbar(im, ax=axes[0].tolist(), shrink=0.8, label="Spearman ρ")
    fig.suptitle(f"Figure B — Spearman correlation skeleton at T* = {T_star} yr",
                 fontsize=11)
    fig.savefig(fig_dir / "figure_b.png", dpi=200, bbox_inches="tight")
    fig.savefig(fig_dir / "figure_b.pdf", bbox_inches="tight")
    print(f"[figure_b] wrote {fig_dir / 'figure_b.png'}")


if __name__ == "__main__":
    main()
