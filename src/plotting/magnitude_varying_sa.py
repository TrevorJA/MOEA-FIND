"""Plotting helpers for Stage-09 magnitude-varying sensitivity analysis.

Three figure types matching the manuscript layout:

- :func:`plot_mv_sa_stacked_area` -- normalized factor share vs.
  magnitude percentile (the Hadjimichael Fig. 4 / Fig. 5 layout).
  This is the *headline* MV-SA figure: one row of stacked colour
  bands per method, x-axis = magnitude percentile, y-axis = factor's
  share of total non-negative sensitivity.

- :func:`plot_mv_sa_lines_with_ci` -- one line per factor with a
  bootstrap CI ribbon, vs. magnitude percentile. Better than the
  stacked area when only 2-3 factors dominate and the reader needs
  to see absolute index magnitudes (and the noise-floor band of the
  control factor).

- :func:`plot_mv_sa_method_panel` -- multi-method comparison: one
  stacked-area panel per method, sharing x-axis. Drives the
  three-method triangulation defence (Delta vs PAWN vs RBD-FAST
  agreement at each percentile).

All functions accept a long-form DataFrame as produced by
:func:`src.sensitivity.magnitude_varying_sa.compute_mv_sa`. The headline-index
column is always ``"headline_index"`` in MV-SA output.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .style import apply_style
from src.sensitivity.magnitude_varying_sa import CONTROL_FACTOR_NAME, stacked_share


# ---------------------------------------------------------------------------
# Stacked area: factor share vs magnitude percentile.
# ---------------------------------------------------------------------------


def plot_mv_sa_stacked_area(
    df: pd.DataFrame,
    *,
    factor_order: Optional[Sequence[str]] = None,
    factor_labels: Optional[Mapping[str, str]] = None,
    cmap_name: str = "tab10",
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
    magnitude_label: str = "Magnitude percentile",
    output_path: Optional[Path] = None,
) -> plt.Axes:
    """Stacked-area plot of normalized factor share vs percentile.

    The control factor (if present) is rendered as a hatched grey band
    so the reader can see the empirical noise floor at every percentile.

    Args:
        df: Long-form output of :func:`src.sensitivity.magnitude_varying_sa.compute_mv_sa`
            (one method).
        factor_order: Order of factors from bottom to top of the stack.
            Default: alphabetical from the data, with the control
            factor (if present) pinned to the bottom so it cannot
            visually mask other factors.
        factor_labels: Display labels per raw factor name.
        cmap_name: Matplotlib colormap for non-control factors.
        ax: Optional Axes; new figure created if None.
        title: Optional axes title.
        magnitude_label: X-axis label (caller specifies what M was).
        output_path: If provided, ``fig.savefig`` is called.

    Returns:
        The matplotlib Axes.
    """
    apply_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.0, 3.6))
    else:
        fig = ax.figure

    factors_in_data = sorted(df["factor"].unique())
    if factor_order is None:
        non_ctrl = [f for f in factors_in_data if f != CONTROL_FACTOR_NAME]
        factor_order = (
            ([CONTROL_FACTOR_NAME] if CONTROL_FACTOR_NAME in factors_in_data
             else []) + non_ctrl
        )

    shares = stacked_share(df, factor_order=factor_order)
    if shares.empty:
        ax.text(0.5, 0.5, "no MV-SA data",
                transform=ax.transAxes, ha="center", va="center")
        return ax

    x = shares.index.values.astype(float)

    # Colour assignment: control is grey-hatched, others draw from cmap.
    cmap = plt.get_cmap(cmap_name, max(1, len(factor_order) - (
        1 if CONTROL_FACTOR_NAME in factor_order else 0
    )))
    color_map = {}
    cidx = 0
    for f in factor_order:
        if f == CONTROL_FACTOR_NAME:
            color_map[f] = "#cccccc"
        else:
            color_map[f] = cmap(cidx)
            cidx += 1

    # Stack: bottom-up cumulative shares.
    cum = np.zeros(len(x))
    for f in factor_order:
        if f not in shares.columns:
            continue
        vals = shares[f].fillna(0.0).values
        ax.fill_between(
            x, cum, cum + vals,
            facecolor=color_map[f],
            edgecolor="white", linewidth=0.4,
            label=factor_labels.get(f, f) if factor_labels else f,
            hatch="//" if f == CONTROL_FACTOR_NAME else None,
            alpha=0.55 if f == CONTROL_FACTOR_NAME else 0.9,
        )
        cum = cum + vals

    ax.set_xlim(float(x.min()), float(x.max()))
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(magnitude_label)
    ax.set_ylabel("Share of total sensitivity")
    if title:
        ax.set_title(title)
    ax.legend(
        loc="center left", bbox_to_anchor=(1.01, 0.5),
        frameon=False, fontsize=8,
    )

    if output_path is not None:
        fig.savefig(Path(output_path))
        plt.close(fig)
    return ax


# ---------------------------------------------------------------------------
# Lines with CI ribbon: per-factor index vs percentile.
# ---------------------------------------------------------------------------


def plot_mv_sa_lines_with_ci(
    df: pd.DataFrame,
    *,
    factor_order: Optional[Sequence[str]] = None,
    factor_labels: Optional[Mapping[str, str]] = None,
    cmap_name: str = "tab10",
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
    magnitude_label: str = "Magnitude percentile",
    headline_label: str = "Sensitivity index",
    output_path: Optional[Path] = None,
) -> plt.Axes:
    """Line plot of per-factor sensitivity vs percentile, with CI ribbons.

    Args:
        df: Long-form MV-SA output (one method).
        factor_order: Z-order of factors (front to back). Default:
            alphabetical, with the control factor last.
        factor_labels: Display labels.
        cmap_name: Matplotlib colormap.
        ax: Optional Axes.
        title: Optional axes title.
        magnitude_label, headline_label: Axis labels.
        output_path: If provided, ``fig.savefig`` is called.

    Returns:
        The matplotlib Axes.
    """
    apply_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.0, 3.8))
    else:
        fig = ax.figure

    factors_in_data = sorted(df["factor"].unique())
    if factor_order is None:
        non_ctrl = [f for f in factors_in_data if f != CONTROL_FACTOR_NAME]
        factor_order = non_ctrl + (
            [CONTROL_FACTOR_NAME] if CONTROL_FACTOR_NAME in factors_in_data
            else []
        )

    cmap = plt.get_cmap(cmap_name, max(1, len(factor_order) - (
        1 if CONTROL_FACTOR_NAME in factor_order else 0
    )))
    color_idx = 0
    for f in factor_order:
        sub = (df[df["factor"] == f]
               .sort_values("percentile"))
        if sub.empty:
            continue
        x = sub["percentile"].values.astype(float)
        y = sub["headline_index"].values.astype(float)
        if f == CONTROL_FACTOR_NAME:
            color = "#888888"
            ls = "--"
        else:
            color = cmap(color_idx)
            color_idx += 1
            ls = "-"
        label = factor_labels.get(f, f) if factor_labels else f
        ax.plot(x, y, color=color, linestyle=ls, linewidth=1.6,
                marker="o", markersize=3.5, label=label)
        if "ci_lo" in sub.columns and sub["ci_lo"].notna().any():
            lo = sub["ci_lo"].values.astype(float)
            hi = sub["ci_hi"].values.astype(float)
            ax.fill_between(x, lo, hi, color=color, alpha=0.15, linewidth=0)

    ax.set_xlabel(magnitude_label)
    ax.set_ylabel(headline_label)
    if title:
        ax.set_title(title)
    ax.legend(loc="best", frameon=False, fontsize=8, ncol=2)

    if output_path is not None:
        fig.savefig(Path(output_path))
        plt.close(fig)
    return ax


# ---------------------------------------------------------------------------
# Multi-method panel: one row per method, shared x-axis.
# ---------------------------------------------------------------------------


def plot_mv_sa_method_panel(
    dfs_by_method: Mapping[str, pd.DataFrame],
    *,
    factor_order: Optional[Sequence[str]] = None,
    factor_labels: Optional[Mapping[str, str]] = None,
    method_labels: Optional[Mapping[str, str]] = None,
    cmap_name: str = "tab10",
    title: Optional[str] = None,
    magnitude_label: str = "Magnitude percentile",
    output_path: Optional[Path] = None,
) -> plt.Figure:
    """Stack of stacked-area panels, one per method, shared x-axis.

    Mirrors Hadjimichael Fig. 4: rows index method (Delta / PAWN /
    RBD-FAST), shared x-axis is magnitude percentile, each panel shows
    the same factor stack so cross-method agreement is visually
    obvious (similar stack profile = method consensus).

    Args:
        dfs_by_method: ``{method_name: long-form DataFrame}``.
        factor_order: As in :func:`plot_mv_sa_stacked_area`.
        factor_labels: Display labels per factor.
        method_labels: Display labels per method (default: uppercase).
        cmap_name: Matplotlib colormap for factors.
        title: Optional figure suptitle.
        magnitude_label: X-axis label.
        output_path: If provided, ``fig.savefig`` is called.

    Returns:
        The matplotlib Figure.
    """
    apply_style()
    methods = list(dfs_by_method.keys())
    n = len(methods)
    if n == 0:
        raise ValueError("dfs_by_method is empty")

    fig, axes = plt.subplots(
        n, 1, figsize=(7.5, 2.6 * n + 0.6), sharex=True,
    )
    if n == 1:
        axes = [axes]

    # Build factor order from the union of factors across methods.
    union = set()
    for df in dfs_by_method.values():
        union.update(df["factor"].unique())
    if factor_order is None:
        non_ctrl = sorted(f for f in union if f != CONTROL_FACTOR_NAME)
        factor_order = (
            ([CONTROL_FACTOR_NAME] if CONTROL_FACTOR_NAME in union else [])
            + non_ctrl
        )

    for ax, m in zip(axes, methods):
        plot_mv_sa_stacked_area(
            dfs_by_method[m],
            factor_order=factor_order,
            factor_labels=factor_labels,
            cmap_name=cmap_name,
            ax=ax,
            title=method_labels.get(m, m.upper()) if method_labels
            else m.upper(),
            magnitude_label=magnitude_label if ax is axes[-1] else "",
        )
        # Single legend on the top panel only; remove from others.
        if ax is not axes[0]:
            leg = ax.get_legend()
            if leg is not None:
                leg.remove()

    if title:
        fig.suptitle(title, y=1.00)
    fig.tight_layout()
    if output_path is not None:
        fig.savefig(Path(output_path))
        plt.close(fig)
    return fig
