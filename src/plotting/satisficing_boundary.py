"""Decision-boundary and manifest-summary plots for the satisficing sweep.

Consumes the artifacts produced by :mod:`src.satisficing_labels`:
    - A per-definition ``model.joblib`` (sklearn classifier)
    - A per-definition ``cv_predictions.csv`` (``y_true``, ``y_prob``)
    - The manifest summary DataFrame (``classifier_summary.csv``)

Figures:
    - :func:`plot_gbt_boundary_2d`  — one-definition boundary in (D1, D2)
      with D3 aggregated to a fixed slice.
    - :func:`plot_manifest_summary` — AUC +/- std bars across definitions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def _satisficing_cmap():
    """Red-white-green diverging cmap for P(satisficing) fills."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(
        "satisficing", ["#d62728", "#ffffff", "#2ca02c"]
    )


def _apply_style():
    try:
        from src.plotting.style import apply_style
        apply_style()
    except Exception:
        pass


def plot_gbt_boundary_2d(
    model,
    X: pd.DataFrame,
    y: np.ndarray,
    x_col: str = "mean_duration",
    y_col: str = "mean_avg_severity",
    agg_col: Optional[str] = "peak_severity_month",
    agg_fn: str = "median",
    output_path: Optional[Path] = None,
    title: str = "",
    n_grid: int = 80,
):
    """Two-panel boundary figure: decision contour + probability heatmap.

    Args:
        model: Fitted sklearn classifier with ``predict_proba``.
        X: Feature DataFrame. Must contain ``x_col``, ``y_col``, and
            ``agg_col`` if not None. Extra columns are allowed and held
            at their aggregate value when generating the grid.
        y: Binary labels for the training points (for scatter overlay).
        x_col / y_col: Feature names for the two axes.
        agg_col: Third feature held at its aggregate on the grid. If the
            model was trained on fewer than three features, pass None.
        agg_fn: Aggregator for columns not on the axes ("median" or "mean").
        output_path: If provided, save PDF here; otherwise return the figure.
        n_grid: Resolution of the decision grid along each axis.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _apply_style()

    aggregator = np.median if agg_fn == "median" else np.mean

    x_lo, x_hi = float(X[x_col].min()), float(X[x_col].max())
    y_lo, y_hi = float(X[y_col].min()), float(X[y_col].max())
    xx, yy = np.meshgrid(
        np.linspace(x_lo, x_hi, n_grid),
        np.linspace(y_lo, y_hi, n_grid),
    )
    grid = pd.DataFrame({x_col: xx.ravel(), y_col: yy.ravel()})
    for col in X.columns:
        if col in (x_col, y_col):
            continue
        grid[col] = float(aggregator(X[col]))
    # Align to training column order
    grid = grid[list(X.columns)]

    probs = model.predict_proba(grid.values)[:, 1].reshape(xx.shape)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel A — decision boundary at p=0.5 with scatter overlay
    ax = axes[0]
    boundary_cmap = _satisficing_cmap()
    ax.contourf(xx, yy, probs, levels=np.linspace(0, 1, 11),
                cmap=boundary_cmap, alpha=0.45)
    cs = ax.contour(xx, yy, probs, levels=[0.5], colors="black", linewidths=1.2)
    ax.clabel(cs, fmt={0.5: "p=0.5"}, inline=True, fontsize=8)
    ax.scatter(X.loc[y == 1, x_col], X.loc[y == 1, y_col],
               s=18, c="#2ca02c", alpha=0.55, edgecolors="none",
               label=f"Satisficing (n={int((y == 1).sum())})")
    ax.scatter(X.loc[y == 0, x_col], X.loc[y == 0, y_col],
               s=18, c="#d62728", alpha=0.55, edgecolors="none",
               label=f"Failure (n={int((y == 0).sum())})")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    if agg_col:
        ax.set_title(f"Decision boundary ({agg_col} at {agg_fn})")
    else:
        ax.set_title("Decision boundary")
    ax.legend(fontsize=8, loc="best")

    # Panel B — probability heatmap
    ax = axes[1]
    im = ax.pcolormesh(xx, yy, probs, cmap="RdYlGn", shading="auto",
                       vmin=0.0, vmax=1.0)
    fig.colorbar(im, ax=ax, label="P(satisficing)")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title("Classifier probability")

    if title:
        fig.suptitle(title, fontsize=12)
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return output_path
    return fig


def plot_manifest_summary(
    summary_df: pd.DataFrame,
    output_path: Optional[Path] = None,
    sort_by: str = "auc_mean",
    title: str = "Satisficing sweep — classifier AUC per definition",
):
    """Horizontal bar chart of AUC mean +/- std per definition."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _apply_style()

    df = summary_df.copy()
    if sort_by in df.columns:
        df = df.sort_values(sort_by, ascending=True)

    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.45 * len(df) + 1.5)))
    ys = np.arange(len(df))
    auc_mean = df["auc_mean"].astype(float).values
    auc_std = df["auc_std"].astype(float).fillna(0).values

    colors = [
        "#2ca02c" if s == "ok" else
        "#ff7f0e" if s == "cv_reduced_folds" else
        "#d62728"
        for s in df["status"].values
    ]

    ax.barh(ys, auc_mean, xerr=auc_std, color=colors, alpha=0.8,
            edgecolor="black", linewidth=0.6)
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=0.7,
               label="random guessing")
    ax.set_yticks(ys)
    ax.set_yticklabels(df["definition_id"].values, fontsize=9)
    ax.set_xlabel("AUC (stratified K-fold CV)")
    ax.set_xlim(0, 1)
    ax.set_title(title, fontsize=11)

    # Annotate n_pos / n_neg
    for y_pos, (_, row) in zip(ys, df.iterrows()):
        ax.text(0.02, y_pos,
                f"n+={int(row['n_pos'])}, n-={int(row['n_neg'])}",
                va="center", fontsize=7, color="black")

    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return output_path
    return fig
