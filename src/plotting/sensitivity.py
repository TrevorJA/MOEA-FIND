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
            :func:`src.sensitivity.sensitivity.convergence_curve` with columns
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


# ---------------------------------------------------------------------------
# Tornado matrix: factor x outcome grouped horizontal bar with CI whiskers.
#
# This is the headline literature-style figure for SI-style SA reports
# (Pianosi & Wagener 2015; Reed group e.g. Quinn et al. 2018; Saltelli
# 2008). For each factor we render one row of bars, one bar per outcome,
# sorted with the most-sensitive outcome on top per factor. Bootstrap
# CIs are drawn as whiskers; an optional reference line marks the
# Delta noise floor (~0.05 for moment-independent at typical n).
# ---------------------------------------------------------------------------


def plot_tornado_matrix(
    indices_long: pd.DataFrame,
    bootstrap_long: Optional[pd.DataFrame] = None,
    *,
    headline_col: str,
    factor_order: Optional[Sequence[str]] = None,
    outcome_order: Optional[Sequence[str]] = None,
    outcome_labels: Optional[Mapping[str, str]] = None,
    factor_labels: Optional[Mapping[str, str]] = None,
    noise_floor: Optional[float] = None,
    title: Optional[str] = None,
    cmap_name: str = "viridis",
    output_path: Optional[Path] = None,
) -> plt.Figure:
    """Headline SA summary: factor-stack of grouped horizontal bars.

    Each factor occupies one row of the figure. Within that row the
    bars are the per-outcome sensitivity index (with bootstrap CI
    whiskers if ``bootstrap_long`` is provided), one outcome per bar.
    A consistent outcome colour is used down the stack so the reader
    can track which outcome each bar refers to across factors.

    Args:
        indices_long: Long-form table with columns ``outcome``,
            ``factor``, and ``headline_col`` (the index value).
            Optional ``ci_lo``, ``ci_hi`` are honoured if present.
        bootstrap_long: Optional long-form bootstrap table
            (``outcome``, ``factor``, ``ci_lo``, ``ci_hi``). Overrides
            ``ci_lo``/``ci_hi`` from ``indices_long``.
        headline_col: Column carrying the index magnitude.
        factor_order: Optional explicit ordering for the factor rows.
            Defaults to factor average ranking (most-important first).
        outcome_order: Optional ordering for outcomes within each
            factor row. Defaults to alphabetical.
        outcome_labels: Display label per raw outcome name.
        factor_labels: Display label per raw factor name.
        noise_floor: If given, plot a vertical dashed reference line
            on every panel.
        title: Suptitle.
        cmap_name: Matplotlib colormap used to assign one colour per
            outcome.
        output_path: If provided, ``fig.savefig`` is called.

    Returns:
        The matplotlib Figure.
    """
    apply_style()

    df = indices_long.copy()
    if bootstrap_long is not None and not bootstrap_long.empty:
        boot = bootstrap_long.set_index(["outcome", "factor"])
        df = df.set_index(["outcome", "factor"])
        for col in ("ci_lo", "ci_hi"):
            if col in boot.columns:
                df[col] = boot[col]
        df = df.reset_index()

    if factor_order is None:
        # Rank factors by mean headline across outcomes, descending.
        factor_order = (
            df.groupby("factor")[headline_col].mean()
              .sort_values(ascending=False)
              .index.tolist()
        )
    if outcome_order is None:
        outcome_order = sorted(df["outcome"].unique())

    n_factors = len(factor_order)
    n_outcomes = len(outcome_order)
    cmap = plt.get_cmap(cmap_name, n_outcomes)
    outcome_colors = {oc: cmap(i) for i, oc in enumerate(outcome_order)}

    fig, axes = plt.subplots(
        n_factors, 1,
        figsize=(7.5, 0.45 * n_outcomes * n_factors + 0.6),
        sharex=True,
    )
    if n_factors == 1:
        axes = [axes]

    xmax = float(df[headline_col].max() * 1.15)

    for ax, factor in zip(axes, factor_order):
        sub = (df[df["factor"] == factor]
               .set_index("outcome")
               .reindex(outcome_order))
        y = np.arange(n_outcomes)
        vals = sub[headline_col].values
        colors = [outcome_colors[oc] for oc in outcome_order]
        ax.barh(y, vals, color=colors, edgecolor="black",
                linewidth=0.5, alpha=0.9)
        if {"ci_lo", "ci_hi"}.issubset(sub.columns):
            lo = sub["ci_lo"].values
            hi = sub["ci_hi"].values
            err_lo = np.maximum(0.0, vals - lo)
            err_hi = np.maximum(0.0, hi - vals)
            ax.errorbar(vals, y, xerr=[err_lo, err_hi], fmt="none",
                        ecolor="black", capsize=2.5, linewidth=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels(
            [outcome_labels.get(oc, oc) if outcome_labels else oc
             for oc in outcome_order],
            fontsize=8,
        )
        ax.invert_yaxis()
        ax.axvline(0.0, color="black", linewidth=0.5)
        if noise_floor is not None:
            ax.axvline(noise_floor, color="grey",
                       linewidth=0.8, linestyle="--",
                       label=f"noise floor ~{noise_floor:.2f}"
                       if ax is axes[0] else None)
        ax.set_xlim(0, xmax)
        flab = factor_labels.get(factor, factor) if factor_labels else factor
        ax.set_ylabel(flab, fontsize=10, fontweight="bold")
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    axes[-1].set_xlabel(headline_col)
    if noise_floor is not None and len(axes) > 0:
        axes[0].legend(loc="lower right", fontsize=7, frameon=False)
    if title:
        fig.suptitle(title, fontsize=11, y=0.995)
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(Path(output_path))
        plt.close(fig)
    return fig


# ---------------------------------------------------------------------------
# Rank-consistency dot plot: factor rank within each outcome.
# ---------------------------------------------------------------------------


def plot_rank_consistency(
    indices_long: pd.DataFrame,
    *,
    headline_col: str,
    factor_order: Optional[Sequence[str]] = None,
    outcome_order: Optional[Sequence[str]] = None,
    outcome_labels: Optional[Mapping[str, str]] = None,
    factor_labels: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    output_path: Optional[Path] = None,
) -> plt.Figure:
    """Dot plot of factor rank (1 = most sensitive) per outcome.

    For each outcome we rank the factors by ``headline_col`` and plot
    a dot at the rank position. Useful for "is the dominant factor
    stable across outcomes?" -- a horizontal line of 1s for one factor
    means it dominates everywhere; vertical scatter means rankings
    flip between outcomes.

    Args:
        indices_long: Long-form table (``outcome``, ``factor``,
            ``headline_col``).
        headline_col: Column to rank by (descending = lowest rank
            number).
        factor_order: Y-axis ordering of factors. Default: alphabetical.
        outcome_order: X-axis ordering of outcomes. Default: alphabetical.
        outcome_labels, factor_labels: Display-name mappings.
        title: Suptitle.
        output_path: If provided, ``fig.savefig`` is called.

    Returns:
        The matplotlib Figure.
    """
    apply_style()
    df = indices_long.copy()
    if factor_order is None:
        factor_order = sorted(df["factor"].unique())
    if outcome_order is None:
        outcome_order = sorted(df["outcome"].unique())

    rank_long: list[tuple[str, str, int]] = []
    for oc in outcome_order:
        sub = (df[df["outcome"] == oc]
               .set_index("factor")
               .reindex(factor_order))
        ranks = sub[headline_col].rank(ascending=False, method="min")
        for f in factor_order:
            r = ranks.get(f, np.nan)
            if np.isfinite(r):
                rank_long.append((oc, f, int(r)))

    n_factors = len(factor_order)
    n_outcomes = len(outcome_order)

    fig, ax = plt.subplots(
        figsize=(0.9 * n_outcomes + 2.5, 0.5 * n_factors + 1.5)
    )

    rank_palette = {1: "#d62728", 2: "#ff9f1c", 3: "#1f77b4"}
    rank_size = {1: 220, 2: 130, 3: 80}

    factor_to_y = {f: i for i, f in enumerate(factor_order)}
    outcome_to_x = {o: j for j, o in enumerate(outcome_order)}

    for oc, f, r in rank_long:
        ax.scatter(
            outcome_to_x[oc], factor_to_y[f],
            s=rank_size.get(r, 60),
            color=rank_palette.get(r, "#888"),
            edgecolor="black", linewidth=0.6, zorder=3,
        )
        ax.text(outcome_to_x[oc], factor_to_y[f], str(r),
                ha="center", va="center", fontsize=8,
                color="white", zorder=4, fontweight="bold")

    ax.set_xticks(np.arange(n_outcomes))
    ax.set_xticklabels(
        [outcome_labels.get(o, o) if outcome_labels else o
         for o in outcome_order],
        rotation=30, ha="right", fontsize=9,
    )
    ax.set_yticks(np.arange(n_factors))
    ax.set_yticklabels(
        [factor_labels.get(f, f) if factor_labels else f
         for f in factor_order],
        fontsize=9,
    )
    ax.invert_yaxis()
    ax.set_xlim(-0.5, n_outcomes - 0.5)
    ax.set_ylim(n_factors - 0.5, -0.5)
    ax.grid(axis="x", linestyle=":", linewidth=0.4, alpha=0.7)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    if title:
        ax.set_title(title, fontsize=11)

    # Custom legend for rank dots.
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", linestyle="",
               markerfacecolor=rank_palette[r],
               markeredgecolor="black", markersize=np.sqrt(rank_size[r]),
               label=f"Rank {r}")
        for r in (1, 2, 3)
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=8)

    fig.tight_layout()
    if output_path is not None:
        fig.savefig(Path(output_path))
        plt.close(fig)
    return fig
