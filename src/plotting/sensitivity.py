"""Plotting helpers for Stage-08 NYC sensitivity analysis.

Tornado, factor-by-outcome heatmap, sample-size convergence, method ×
method rank-correlation, and outcome × outcome rank-correlation. Each
function accepts an Axes (creates a fresh figure if ``ax`` is None) and
either returns the Axes (for composition) or writes a file when
``output_path`` is supplied.

All figures use :func:`src.plotting.style.apply_style` for consistent
fonts, DPI, and spines.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .style import apply_style


# ---------------------------------------------------------------------------
# Tornado: per-(outcome, method) bar of indices with CI whiskers.
# ---------------------------------------------------------------------------


def plot_tornado(
    indices_df: pd.DataFrame,
    *,
    headline_col: str,
    ci_lo_col: str = "ci_lo",
    ci_hi_col: str = "ci_hi",
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
    color: str = "#1f77b4",
    output_path: Optional[Path] = None,
) -> plt.Axes:
    """Tornado plot of factor sensitivity indices with CI whiskers.

    Args:
        indices_df: Per-factor DataFrame (index = factor names) with
            ``headline_col`` and the two CI columns.
        headline_col: Column carrying the index magnitude (e.g.,
            ``delta``, ``median``, ``S1``).
        ci_lo_col, ci_hi_col: Column names for the lower and upper CI.
        ax: Optional matplotlib Axes; a fresh figure is created if None.
        title: Optional axes title.
        color: Bar fill color.
        output_path: If provided, ``fig.savefig`` is called and the path
            is returned via ``ax.figure``.

    Returns:
        The matplotlib Axes.
    """
    apply_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(5.0, 0.45 * len(indices_df) + 1.0))
    else:
        fig = ax.figure

    df = indices_df.copy().sort_values(headline_col, ascending=True)
    y = np.arange(len(df))
    vals = df[headline_col].values
    lo = df[ci_lo_col].values
    hi = df[ci_hi_col].values
    err_lo = np.maximum(0.0, vals - lo)
    err_hi = np.maximum(0.0, hi - vals)

    ax.barh(y, vals, color=color, edgecolor="black", linewidth=0.6,
            alpha=0.85)
    ax.errorbar(vals, y, xerr=[err_lo, err_hi], fmt="none",
                ecolor="black", capsize=2.5, linewidth=0.8)
    ax.axvline(0.0, color="black", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(df.index)
    ax.set_xlabel(headline_col)
    if title:
        ax.set_title(title)

    if output_path is not None:
        fig.savefig(Path(output_path))
        plt.close(fig)
    return ax


# ---------------------------------------------------------------------------
# Heatmap: factor × outcome for one method.
# ---------------------------------------------------------------------------


def plot_index_heatmap(
    method_outcome_indices: Mapping[str, pd.DataFrame],
    *,
    method_name: str,
    headline_col: str,
    factor_order: Optional[Sequence[str]] = None,
    output_path: Optional[Path] = None,
    title: Optional[str] = None,
    cmap: str = "viridis",
) -> plt.Axes:
    """Factor × outcome heatmap of sensitivity indices for one method.

    Args:
        method_outcome_indices: ``{outcome_name: DataFrame}`` where each
            DataFrame is the per-factor index table for the same method.
        method_name: Used in the default title only.
        headline_col: Column carrying the index magnitude.
        factor_order: Optional explicit factor ordering for the y-axis.
            Defaults to the index of the first DataFrame.
        output_path: If provided, ``fig.savefig`` is called.
        title: Optional axes title.
        cmap: Matplotlib colormap.

    Returns:
        The matplotlib Axes.
    """
    apply_style()
    outcomes = list(method_outcome_indices.keys())
    if factor_order is None:
        factor_order = list(next(iter(method_outcome_indices.values())).index)

    matrix = np.full((len(factor_order), len(outcomes)), np.nan, dtype=float)
    for j, oc in enumerate(outcomes):
        df = method_outcome_indices[oc]
        for i, f in enumerate(factor_order):
            if f in df.index:
                matrix[i, j] = float(df.loc[f, headline_col])

    fig, ax = plt.subplots(
        figsize=(0.9 * len(outcomes) + 2.5, 0.4 * len(factor_order) + 1.5)
    )
    im = ax.imshow(matrix, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(len(outcomes)))
    ax.set_xticklabels(outcomes, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(factor_order)))
    ax.set_yticklabels(factor_order)
    for i in range(len(factor_order)):
        for j in range(len(outcomes)):
            v = matrix[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="white" if v > matrix[np.isfinite(matrix)].mean()
                        else "black", fontsize=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(headline_col)
    ax.set_title(title or f"{method_name}: {headline_col} per factor × outcome")

    if output_path is not None:
        fig.savefig(Path(output_path))
        plt.close(fig)
    return ax


# ---------------------------------------------------------------------------
# Convergence: index vs n with error band.
# ---------------------------------------------------------------------------


def plot_convergence(
    curves_df: pd.DataFrame,
    *,
    title: Optional[str] = None,
    output_path: Optional[Path] = None,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Sample-size convergence: index vs n per factor with 5–95% band.

    Args:
        curves_df: Long-form DataFrame from
            :func:`src.sensitivity.convergence_curve` with columns
            ``factor``, ``n``, ``mean``, ``p05``, ``p95``.
        title: Optional axes title.
        output_path: If provided, ``fig.savefig`` is called.
        ax: Optional matplotlib Axes.

    Returns:
        The matplotlib Axes.
    """
    apply_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 4.0))
    else:
        fig = ax.figure

    factors = list(curves_df["factor"].unique())
    cmap = plt.cm.tab10
    for k, f in enumerate(factors):
        sub = curves_df[curves_df["factor"] == f].sort_values("n")
        color = cmap(k % 10)
        ax.plot(sub["n"], sub["mean"], color=color, marker="o",
                linewidth=1.2, markersize=4, label=f)
        ax.fill_between(sub["n"], sub["p05"], sub["p95"],
                        color=color, alpha=0.18, linewidth=0)
    ax.set_xlabel("Sample size n")
    ax.set_ylabel("Sensitivity index")
    ax.set_xscale("log")
    ax.legend(frameon=False, fontsize=8, loc="best", ncol=2)
    if title:
        ax.set_title(title)

    if output_path is not None:
        fig.savefig(Path(output_path))
        plt.close(fig)
    return ax


# ---------------------------------------------------------------------------
# Rank-correlation heatmap (method × method or outcome × outcome).
# ---------------------------------------------------------------------------


def plot_rank_correlation(
    rho_df: pd.DataFrame,
    *,
    title: Optional[str] = None,
    output_path: Optional[Path] = None,
    ax: Optional[plt.Axes] = None,
    cmap: str = "RdBu_r",
    vmin: float = -1.0,
    vmax: float = 1.0,
) -> plt.Axes:
    """Spearman rank-correlation heatmap with annotated cells.

    Args:
        rho_df: Square DataFrame; row and column labels are used.
        title: Optional axes title.
        output_path: If provided, ``fig.savefig`` is called.
        ax: Optional matplotlib Axes.
        cmap, vmin, vmax: Matplotlib colorbar parameters.

    Returns:
        The matplotlib Axes.
    """
    apply_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(0.6 * len(rho_df) + 2.0,
                                         0.6 * len(rho_df) + 1.5))
    else:
        fig = ax.figure

    arr = rho_df.values.astype(float)
    im = ax.imshow(arr, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(rho_df.columns)))
    ax.set_xticklabels(rho_df.columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(rho_df.index)))
    ax.set_yticklabels(rho_df.index)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=8,
                        color="white" if abs(v) > 0.55 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Spearman rho")
    if title:
        ax.set_title(title)

    if output_path is not None:
        fig.savefig(Path(output_path))
        plt.close(fig)
    return ax
