"""Render diagnostic figures for the drought-metric selection workflow.

Reads ``outputs/02_calibration/metric_explorer/`` and writes:

* ``spearman_clustermap.{pdf,png}`` — Spearman ρ heatmap with hierarchical
  dendrogram (the load-bearing figure for the "which pairs are redundant"
  question).
* ``pearson_clustermap.{pdf,png}`` — same on Pearson correlation.
* ``spread_bar.{pdf,png}`` — robust spread (IQR/(|median|+σ)) bar, ranked,
  color-coded by hydrologic concept; cluster-representative status is
  shown as bold edge.
* ``recommended_k3_3d.{pdf,png}`` — 3D scatter of T-blocks on the
  recommended K=3 set.
* ``recommended_k3_scatter_matrix.{pdf,png}`` — pair plot of the
  recommended K=3 set.
* ``recommended_k4_scatter_matrix.{pdf,png}`` — pair plot of the
  recommended K=4 set.

PCA-based figures and metrics are intentionally not produced: the
recommendation is purely in terms of single, named hazard characteristics.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers projection)
import numpy as np
import pandas as pd
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "02_calibration"
DRIVER = "metric_explorer"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Clustermap (Spearman / Pearson)
# ---------------------------------------------------------------------------


def plot_corr_clustermap(
    corr: pd.DataFrame,
    title: str,
    fig_path: Path,
    *,
    cmap: str = "RdBu_r",
) -> None:
    """Seaborn clustermap with dendrogram on a correlation matrix."""
    if corr.shape[0] < 2:
        return
    g = sns.clustermap(
        corr, cmap=cmap, vmin=-1.0, vmax=1.0, center=0.0,
        figsize=(0.32 * len(corr) + 4.0, 0.32 * len(corr) + 4.0),
        annot=False, cbar_kws={"label": "ρ"},
        method="average", metric="euclidean",
    )
    g.fig.suptitle(title, y=1.02, fontsize=11)
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(), rotation=70, fontsize=8, ha="right"
    )
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8
    )
    g.fig.savefig(fig_path.with_suffix(".pdf"), bbox_inches="tight")
    g.fig.savefig(fig_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(g.fig)


# ---------------------------------------------------------------------------
# Spread bar chart
# ---------------------------------------------------------------------------


def plot_spread_bar(
    spread_df: pd.DataFrame,
    cluster_df: pd.DataFrame,
    fig_path: Path,
) -> None:
    """Robust spread bar chart, ranked, color-coded by hydrologic concept.

    Cluster-representative status is shown as bold edge. Metrics that
    failed the screen are drawn in light grey on the right.
    """
    sub = spread_df.copy()
    cluster_lookup = cluster_df.set_index("metric")["cluster_id"].to_dict()
    rep_lookup = cluster_df.set_index("metric")["is_representative"].to_dict()
    concept_lookup = cluster_df.set_index("metric")["concept"].to_dict()
    sub["cluster_id"] = sub["metric"].map(cluster_lookup)
    sub["is_representative"] = sub["metric"].map(rep_lookup).fillna(False)
    if "concept" not in sub.columns:
        sub["concept"] = sub["metric"].map(concept_lookup).fillna("unknown")

    sub = sub.sort_values(
        ["passes_screen", "spread_score"], ascending=[False, False]
    ).reset_index(drop=True)

    concepts = sorted([c for c in sub["concept"].unique() if isinstance(c, str)])
    cmap = plt.get_cmap("tab20")
    concept_color = {c: cmap(i % cmap.N) for i, c in enumerate(concepts)}

    n = len(sub)
    fig, ax = plt.subplots(figsize=(max(8, 0.32 * n), 5.5))
    for i, row in sub.iterrows():
        if not row["passes_screen"]:
            color = "lightgrey"
        else:
            color = concept_color.get(row["concept"], "lightgrey")
        edgecolor = "black" if row["is_representative"] else "none"
        lw = 1.4 if row["is_representative"] else 0.0
        score = row["spread_score"] if np.isfinite(row["spread_score"]) else 0.0
        ax.bar(i, score, color=color, edgecolor=edgecolor, linewidth=lw)
    ax.set_xticks(range(n))
    ax.set_xticklabels(sub["metric"], rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("Robust spread score: IQR / (|median| + σ)")
    ax.set_title("Per-metric spread across historical T-blocks\n"
                 "color = hydrologic concept; bold edge = cluster representative; "
                 "grey = failed screen")

    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color=concept_color[c], label=c)
        for c in concepts
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=7,
              ncol=2, framealpha=0.9, title="concept")
    fig.tight_layout()
    _save(fig, fig_path)


# ---------------------------------------------------------------------------
# PCA biplot
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Recommended set 3D scatter and pair plot
# ---------------------------------------------------------------------------


def plot_recommended_3d(
    chars_df: pd.DataFrame,
    metrics: List[str],
    fig_path: Path,
) -> None:
    """3D scatter of T-blocks projected onto the recommended K=3 set."""
    if len(metrics) < 3:
        return
    pts = chars_df[metrics].astype(float).values
    fig = plt.figure(figsize=(7.0, 6.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        pts[:, 0], pts[:, 1], pts[:, 2],
        c=np.arange(len(pts)), cmap="viridis",
        s=24, alpha=0.85, edgecolors="black", linewidths=0.3,
    )
    ax.set_xlabel(metrics[0], fontsize=9)
    ax.set_ylabel(metrics[1], fontsize=9)
    ax.set_zlabel(metrics[2], fontsize=9)
    ax.set_title(
        "Recommended K=3 metric set\n"
        f"{len(pts)} historical T-year blocks (color = chronological)",
        fontsize=10,
    )
    fig.tight_layout()
    _save(fig, fig_path)


def plot_scatter_matrix(
    chars_df: pd.DataFrame,
    metrics: List[str],
    fig_path: Path,
) -> None:
    """Pair plot of the recommended metric set."""
    if len(metrics) < 2:
        return
    sub = chars_df[metrics].astype(float).copy()
    sub["block_idx"] = np.arange(len(sub))
    g = sns.pairplot(
        sub, vars=metrics, diag_kind="hist",
        plot_kws={"s": 18, "alpha": 0.7, "edgecolor": "black", "linewidth": 0.2},
        corner=True,
    )
    g.fig.suptitle("Recommended set — pairwise relationships across T-blocks",
                   y=1.02, fontsize=11)
    g.fig.savefig(fig_path.with_suffix(".pdf"), bbox_inches="tight")
    g.fig.savefig(fig_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(g.fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _load_recommended(json_path: Path) -> List[str]:
    """Pull the top-1 metric tuple from a {k3,k4}_alternatives.json file."""
    if not json_path.exists():
        return []
    payload = json.loads(json_path.read_text())
    alts = payload.get("alternatives", [])
    if not alts:
        return []
    return list(alts[0].get("metrics", []))


def main():
    in_dir = stage_output_dir(STAGE, DRIVER, create=False)
    fig_dir = stage_figure_dir(STAGE, DRIVER)
    print(f"[plots/02/metric_explorer] in_dir={in_dir} fig_dir={fig_dir}")

    chars_path = in_dir / "block_chars_extended.csv"
    spread_path = in_dir / "per_metric_spread.csv"
    spearman_path = in_dir / "spearman_corr.csv"
    pearson_path = in_dir / "pearson_corr.csv"
    clusters_path = in_dir / "clusters.json"

    missing = [p for p in [chars_path, spread_path, spearman_path,
                           pearson_path, clusters_path]
               if not p.exists()]
    if missing:
        sys.exit(f"missing inputs: {missing} — run the compute driver first")

    chars_df = pd.read_csv(chars_path)
    spread_df = pd.read_csv(spread_path)
    spearman = pd.read_csv(spearman_path, index_col=0)
    pearson = pd.read_csv(pearson_path, index_col=0)
    cluster_df = pd.DataFrame(json.loads(clusters_path.read_text()))

    k3 = _load_recommended(in_dir / "k3_alternatives.json")
    k4 = _load_recommended(in_dir / "k4_alternatives.json")
    print(f"[plots] recommended K=3 = {k3}")
    print(f"[plots] recommended K=4 = {k4}")

    plot_corr_clustermap(
        spearman, "Spearman ρ — surviving candidate metrics",
        fig_dir / "spearman_clustermap",
    )
    plot_corr_clustermap(
        pearson, "Pearson r — surviving candidate metrics",
        fig_dir / "pearson_clustermap",
    )
    plot_spread_bar(spread_df, cluster_df, fig_dir / "spread_bar")

    if len(k3) >= 3:
        plot_recommended_3d(chars_df, k3, fig_dir / "recommended_k3_3d")
        plot_scatter_matrix(chars_df, k3, fig_dir / "recommended_k3_scatter_matrix")
    if len(k4) >= 2:
        plot_scatter_matrix(chars_df, k4, fig_dir / "recommended_k4_scatter_matrix")

    # Remove any stale PCA-era artifacts from previous runs.
    for stale in ["pca_biplot.pdf", "pca_biplot.png"]:
        sp = fig_dir / stale
        if sp.exists():
            sp.unlink()

    print("[plots] done.")


if __name__ == "__main__":
    main()
