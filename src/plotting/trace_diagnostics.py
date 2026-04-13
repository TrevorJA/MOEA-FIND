"""Synthetic trace plausibility diagnostics.

Reusable diagnostic plots comparing properties of generated synthetic
traces against the historical record. Designed to work with any generator
mode (empirical, KDE, parametric) and any number of traces.

All functions accept data arrays and return (fig, axes) tuples for
flexible composition into diagnostic panels.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt

from src.plotting.style import COLORS, WATER_YEAR_MONTHS, apply_style


def plot_hydrograph(
    ax: plt.Axes,
    monthly_flows: np.ndarray,
    threshold: float,
    title: str = "",
    color: str = COLORS["empirical"],
    highlight_droughts: bool = True,
):
    """Plot a monthly flow hydrograph with drought events highlighted.

    Args:
        ax: Matplotlib axes to plot on.
        monthly_flows: 1D array of monthly flows.
        threshold: Drought threshold for highlighting.
        title: Axes title.
        color: Line color.
        highlight_droughts: Whether to shade below-threshold periods.
    """
    t = np.arange(len(monthly_flows))
    ax.plot(t, monthly_flows, color=color, linewidth=0.8, alpha=0.8)
    ax.axhline(threshold, color=COLORS["muted"], linestyle="--",
               linewidth=0.7, alpha=0.6, label="P20 threshold")

    if highlight_droughts:
        below = monthly_flows < threshold
        ax.fill_between(t, monthly_flows, threshold,
                        where=below, alpha=0.2, color=COLORS["anti_ideal"],
                        label="Drought")

    ax.set_ylabel("Flow (cfs)")
    if title:
        ax.set_title(title, fontsize=10)


def plot_hydrograph_panel(
    synthetic_traces: List[np.ndarray],
    historical_1d: np.ndarray,
    threshold: float,
    trace_labels: Optional[List[str]] = None,
    trace_chars: Optional[List[dict]] = None,
    figsize: Tuple[float, float] = (14, 10),
) -> Tuple[plt.Figure, np.ndarray]:
    """Panel of synthetic hydrographs from different drought regions.

    Args:
        synthetic_traces: List of 1D monthly flow arrays.
        historical_1d: Historical monthly flows (1D) for reference.
        threshold: Drought threshold.
        trace_labels: Optional label for each trace.
        trace_chars: Optional drought characteristics dict per trace.
        figsize: Figure size.

    Returns:
        (fig, axes) tuple.
    """
    apply_style()
    n = len(synthetic_traces)
    ncols = min(n, 2)
    nrows = (n + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=False, sharey=True)
    if n == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for i, trace in enumerate(synthetic_traces):
        ax = axes[i]
        label = trace_labels[i] if trace_labels else f"Trace {i+1}"
        plot_hydrograph(ax, trace, threshold, title=label,
                        color=COLORS["empirical"])

        if trace_chars and i < len(trace_chars):
            chars = trace_chars[i]
            info = (f"Dur: {chars['mean_duration']:.1f} mo\n"
                    f"Int: {chars['mean_intensity']:.1f} cfs\n"
                    f"Events: {chars['n_events']}")
            ax.text(0.02, 0.98, info, transform=ax.transAxes, va="top",
                    fontsize=8, bbox=dict(boxstyle="round,pad=0.3",
                    facecolor="white", alpha=0.9))

        ax.set_xlabel("Month index")

    # Hide extra axes
    for i in range(n, len(axes)):
        axes[i].set_visible(False)

    fig.tight_layout()
    return fig, axes


def plot_autocorrelation_comparison(
    synthetic_traces: List[np.ndarray],
    historical_1d: np.ndarray,
    max_lag: int = 24,
    trace_labels: Optional[List[str]] = None,
    figsize: Tuple[float, float] = (8, 5),
) -> Tuple[plt.Figure, plt.Axes]:
    """Compare autocorrelation functions of synthetic traces vs historical.

    Args:
        synthetic_traces: List of 1D monthly flow arrays.
        historical_1d: Historical monthly flows (1D).
        max_lag: Maximum lag to compute.
        trace_labels: Optional labels for each trace.
        figsize: Figure size.

    Returns:
        (fig, ax) tuple.
    """
    apply_style()

    def acf(x, max_lag):
        n = len(x)
        mean = np.mean(x)
        var = np.var(x)
        if var < 1e-12:
            return np.zeros(max_lag + 1)
        result = np.zeros(max_lag + 1)
        for lag in range(max_lag + 1):
            cov = np.mean((x[:n-lag] - mean) * (x[lag:] - mean))
            result[lag] = cov / var
        return result

    fig, ax = plt.subplots(figsize=figsize)
    lags = np.arange(max_lag + 1)

    # Historical
    hist_acf = acf(historical_1d, max_lag)
    ax.plot(lags, hist_acf, color=COLORS["historical"], linewidth=2,
            label="Historical", zorder=10)

    # Synthetic traces
    synth_acfs = []
    for i, trace in enumerate(synthetic_traces):
        trace_acf = acf(trace, max_lag)
        synth_acfs.append(trace_acf)

    synth_acfs = np.array(synth_acfs)

    if len(synth_acfs) > 5:
        # Show envelope (10th-90th percentile) + median
        p10 = np.percentile(synth_acfs, 10, axis=0)
        p50 = np.percentile(synth_acfs, 50, axis=0)
        p90 = np.percentile(synth_acfs, 90, axis=0)
        ax.fill_between(lags, p10, p90, alpha=0.2, color=COLORS["empirical"],
                         label="Synthetic (10-90th pctl)")
        ax.plot(lags, p50, color=COLORS["empirical"], linewidth=1.5,
                label="Synthetic median", linestyle="--")
    else:
        for i, trace_acf in enumerate(synth_acfs):
            lbl = trace_labels[i] if trace_labels else f"Trace {i+1}"
            ax.plot(lags, trace_acf, alpha=0.6, linewidth=1, label=lbl)

    ax.set_xlabel("Lag (months)")
    ax.set_ylabel("Autocorrelation")
    ax.set_title("Autocorrelation Comparison")
    ax.legend(fontsize=8)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle=":")
    fig.tight_layout()
    return fig, ax


def plot_seasonal_cycle_comparison(
    synthetic_traces_2d: List[np.ndarray],
    historical_2d: np.ndarray,
    figsize: Tuple[float, float] = (8, 5),
) -> Tuple[plt.Figure, plt.Axes]:
    """Compare monthly mean and std of synthetic traces vs historical.

    Args:
        synthetic_traces_2d: List of (n_years, 12) arrays.
        historical_2d: Historical monthly flows (n_years, 12).
        figsize: Figure size.

    Returns:
        (fig, axes) tuple with mean and std panels.
    """
    apply_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    months = np.arange(12)

    hist_mean = np.mean(historical_2d, axis=0)
    hist_std = np.std(historical_2d, axis=0, ddof=1)

    ax1.plot(months, hist_mean, "o-", color=COLORS["historical"],
             linewidth=2, markersize=5, label="Historical", zorder=10)
    ax2.plot(months, hist_std, "o-", color=COLORS["historical"],
             linewidth=2, markersize=5, label="Historical", zorder=10)

    if len(synthetic_traces_2d) > 5:
        syn_means = np.array([np.mean(t, axis=0) for t in synthetic_traces_2d])
        syn_stds = np.array([np.std(t, axis=0, ddof=1) for t in synthetic_traces_2d])

        for ax, syn_stat in [(ax1, syn_means), (ax2, syn_stds)]:
            p10 = np.percentile(syn_stat, 10, axis=0)
            p50 = np.percentile(syn_stat, 50, axis=0)
            p90 = np.percentile(syn_stat, 90, axis=0)
            ax.fill_between(months, p10, p90, alpha=0.2, color=COLORS["empirical"])
            ax.plot(months, p50, "--", color=COLORS["empirical"],
                    linewidth=1.5, label="Synthetic median")
    else:
        for i, trace_2d in enumerate(synthetic_traces_2d):
            lbl = f"Trace {i+1}"
            ax1.plot(months, np.mean(trace_2d, axis=0), alpha=0.5, linewidth=1, label=lbl)
            ax2.plot(months, np.std(trace_2d, axis=0, ddof=1), alpha=0.5, linewidth=1)

    for ax in [ax1, ax2]:
        ax.set_xticks(months)
        ax.set_xticklabels(WATER_YEAR_MONTHS, rotation=45, fontsize=8)
        ax.legend(fontsize=8)

    ax1.set_ylabel("Mean Monthly Flow (cfs)")
    ax1.set_title("Seasonal Mean")
    ax2.set_ylabel("Std Monthly Flow (cfs)")
    ax2.set_title("Seasonal Variability")

    fig.tight_layout()
    return fig, (ax1, ax2)


def plot_flow_duration_curve(
    synthetic_traces: List[np.ndarray],
    historical_1d: np.ndarray,
    figsize: Tuple[float, float] = (8, 5),
) -> Tuple[plt.Figure, plt.Axes]:
    """Compare flow duration curves of synthetic traces vs historical.

    Args:
        synthetic_traces: List of 1D monthly flow arrays.
        historical_1d: Historical monthly flows.
        figsize: Figure size.

    Returns:
        (fig, ax) tuple.
    """
    apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    def fdc(flows):
        sorted_flows = np.sort(flows)[::-1]
        exceedance = np.arange(1, len(sorted_flows) + 1) / (len(sorted_flows) + 1)
        return exceedance, sorted_flows

    exc_h, flows_h = fdc(historical_1d)
    ax.semilogy(exc_h * 100, flows_h, color=COLORS["historical"],
                linewidth=2, label="Historical", zorder=10)

    if len(synthetic_traces) > 5:
        # Envelope
        all_exc = np.linspace(0.01, 0.99, 200)
        interp_flows = []
        for trace in synthetic_traces:
            exc_s, flows_s = fdc(trace)
            interp_flows.append(np.interp(all_exc, exc_s, flows_s))
        interp_flows = np.array(interp_flows)
        p10 = np.percentile(interp_flows, 10, axis=0)
        p50 = np.percentile(interp_flows, 50, axis=0)
        p90 = np.percentile(interp_flows, 90, axis=0)
        ax.fill_between(all_exc * 100, p10, p90, alpha=0.2, color=COLORS["empirical"])
        ax.semilogy(all_exc * 100, p50, "--", color=COLORS["empirical"],
                     linewidth=1.5, label="Synthetic median")
    else:
        for i, trace in enumerate(synthetic_traces):
            exc_s, flows_s = fdc(trace)
            ax.semilogy(exc_s * 100, flows_s, alpha=0.4, linewidth=0.8)

    ax.set_xlabel("Exceedance Probability (%)")
    ax.set_ylabel("Flow (cfs)")
    ax.set_title("Flow Duration Curve")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig, ax
