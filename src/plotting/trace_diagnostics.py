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


# -----------------------------------------------------------------------------
# Figure 5 — SSI-3 drought metric definitions example
# -----------------------------------------------------------------------------
#
# Manuscript §4.3 / Figure 5. This is a didactic figure: it shows how SSI-3
# is computed from a historical monthly streamflow series and how mean
# severity, mean duration, and peak severity month are extracted from a
# representative drought event. The figure is rendered from a cached
# daily USGS CSV under ``outputs/data_cache/`` — not from SynHydro — so
# it is self-contained and renders without the SynHydro dependency.
# -----------------------------------------------------------------------------


def _daily_csv_to_monthly(csv_path) -> "pd.Series":
    """Load a cached daily USGS CSV and aggregate to monthly mean flow."""
    import pandas as pd
    df = pd.read_csv(csv_path, index_col=0, parse_dates=True)
    col = df.columns[0] if df.shape[1] >= 1 else "flow_cfs"
    monthly = df[col].resample("MS").mean()
    monthly.name = "flow_cfs"
    return monthly


def _ssi3_from_monthly(monthly: "pd.Series") -> "pd.Series":
    """Compute SSI-3 from a monthly flow series using a per-month gamma fit.

    Light-weight replacement for SynHydro's SSI calculator sufficient for
    a didactic figure. Steps:
        1. Take a centered 3-month rolling mean (SSI-3 accumulation).
        2. For each calendar month, fit a two-parameter gamma to the
           non-zero historical values.
        3. Transform the 3-month accumulated flow to a standard normal
           via the gamma CDF and the inverse normal CDF.
    """
    import pandas as pd
    from scipy.stats import gamma, norm

    acc = monthly.rolling(3, min_periods=3).mean()
    ssi = pd.Series(index=acc.index, dtype=float)
    for m in range(1, 13):
        mask = acc.index.month == m
        vals = acc[mask].dropna()
        positive = vals[vals > 0]
        if len(positive) < 10:
            ssi[mask] = np.nan
            continue
        shape, loc, scale = gamma.fit(positive, floc=0.0)
        p = gamma.cdf(acc[mask].values, shape, loc=loc, scale=scale)
        # Clip to (0,1) to avoid inf under inverse normal
        p = np.clip(p, 1e-6, 1 - 1e-6)
        ssi[mask] = norm.ppf(p)
    return ssi


def _detect_drought_events(ssi: "pd.Series") -> list:
    """Detect SSI-based drought events using the SynHydro convention.

    Rules:
        - Onset at SSI < 0.
        - Event becomes "critical" once SSI <= -1.
        - Termination after three consecutive SSI > 0 months.
        - Only critical events are returned.

    Returns a list of dicts with keys: ``start``, ``end``, ``duration``,
    ``severity`` (min SSI), ``avg_severity`` (mean SSI), and
    ``max_severity_date`` (the timestamp of the minimum SSI).
    """
    events = []
    drought_idx: list = []
    in_critical = False
    pos = 0
    vals = ssi.values
    dates = ssi.index
    for i, v in enumerate(vals):
        if np.isnan(v):
            continue
        if v < 0:
            drought_idx.append(i)
            pos = 0
            if v <= -1:
                in_critical = True
        else:
            if in_critical:
                pos += 1
                if pos >= 3:
                    arr = vals[drought_idx]
                    peak_i = drought_idx[int(np.argmin(arr))]
                    events.append({
                        "start": dates[drought_idx[0]],
                        "end": dates[drought_idx[-1]],
                        "duration": len(drought_idx),
                        "severity": float(arr.min()),
                        "avg_severity": float(arr.mean()),
                        "max_severity_date": dates[peak_i],
                    })
                    in_critical = False
                    drought_idx = []
                    pos = 0
            else:
                drought_idx = []
                pos = 0
    return events


def fig5_ssi_definitions_example(
    monthly_flow_csv,
    figsize: Tuple[float, float] = (7.0, 4.8),
    event_choice: str = "most_severe",
):
    """Manuscript §4.3 / Figure 5 — SSI-3 drought metric definitions example.

    Three-panel figure:

    (a) Full historical SSI-3 time series with the ``SSI = -1`` critical
        threshold line and all detected critical drought events shaded.
    (b) Zoom on a representative drought event with onset, trough (peak
        severity month), and recovery annotated.
    (c) Annotation box defining mean severity, mean duration, and mean
        peak severity month as used by the MOEA-FIND objectives.

    Parameters
    ----------
    monthly_flow_csv : str or Path
        Path to a cached daily USGS CSV under ``outputs/data_cache/``.
        The file must have a datetime index and a single flow column.
    figsize : tuple of float
        Figure size in inches.
    event_choice : {"most_severe", "longest", "first"}
        How to pick the representative event drawn in panel (b).

    Returns
    -------
    matplotlib.figure.Figure
    """
    from matplotlib.gridspec import GridSpec

    apply_style()

    monthly = _daily_csv_to_monthly(monthly_flow_csv)
    ssi = _ssi3_from_monthly(monthly)
    events = _detect_drought_events(ssi)
    if not events:
        raise RuntimeError("No critical drought events detected in the "
                           "cached historical record; Figure 5 needs real "
                           "drought events to illustrate the definitions.")

    if event_choice == "most_severe":
        ev = min(events, key=lambda e: e["severity"])
    elif event_choice == "longest":
        ev = max(events, key=lambda e: e["duration"])
    else:
        ev = events[0]

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.2],
                  width_ratios=[1.4, 1.0], hspace=0.45, wspace=0.3)

    # -- Panel (a): full SSI-3 timeline --
    ax_a = fig.add_subplot(gs[0, :])
    ax_a.plot(ssi.index, ssi.values, color=COLORS["historical"], linewidth=0.8)
    ax_a.axhline(-1.0, color=COLORS["anti_ideal"], linestyle="--",
                 linewidth=0.8, label="SSI = -1 (critical)")
    for e in events:
        ax_a.axvspan(e["start"], e["end"], alpha=0.18,
                     color=COLORS["anti_ideal"])
    ax_a.axvspan(ev["start"], ev["end"], alpha=0.35,
                 color=COLORS["highlight"], label="representative event")
    ax_a.set_ylabel("SSI-3")
    ax_a.set_title(f"(a) Historical SSI-3 with {len(events)} critical events")
    ax_a.legend(fontsize=7, loc="lower left", framealpha=0.9)

    # -- Panel (b): zoom on the representative event --
    ax_b = fig.add_subplot(gs[1, 0])
    pad = max(6, ev["duration"])
    zoom_start = ev["start"] - pd_offset_months(pad)
    zoom_end = ev["end"] + pd_offset_months(pad)
    mask = (ssi.index >= zoom_start) & (ssi.index <= zoom_end)
    ax_b.plot(ssi.index[mask], ssi.values[mask],
              color=COLORS["historical"], linewidth=1.1)
    ax_b.axhline(-1.0, color=COLORS["anti_ideal"], linestyle="--",
                 linewidth=0.8)
    ax_b.axhline(0.0, color=COLORS["muted"], linestyle=":", linewidth=0.6)
    ax_b.axvspan(ev["start"], ev["end"], alpha=0.25,
                 color=COLORS["highlight"])
    ax_b.plot(ev["max_severity_date"], ev["severity"], "v",
              color=COLORS["anti_ideal"], markersize=9,
              markeredgecolor="black", markeredgewidth=0.6)
    ax_b.annotate(
        f"peak\nmonth {ev['max_severity_date'].month}",
        xy=(ev["max_severity_date"], ev["severity"]),
        xytext=(8, -10), textcoords="offset points", fontsize=7,
        arrowprops=dict(arrowstyle="->", lw=0.6),
    )
    ax_b.set_ylabel("SSI-3")
    ax_b.set_title("(b) Representative event")
    ax_b.tick_params(axis="x", labelsize=7)
    # Light rotation of x tick labels
    for lbl in ax_b.get_xticklabels():
        lbl.set_rotation(25)
        lbl.set_ha("right")

    # -- Panel (c): definitions box --
    ax_c = fig.add_subplot(gs[1, 1])
    ax_c.axis("off")
    ax_c.set_title("(c) Drought characteristic definitions")
    lines = [
        r"$d_1$  mean severity   $= \overline{\min\,\mathrm{SSI}_{\mathrm{event}}}$",
        r"$d_2$  mean duration   $=$ mean months below onset",
        r"$d_3$  peak severity   $=$ cyclic mean of months of",
        r"        month               per-event SSI minimum",
        "",
        r"L1 distance to target $D^*$:",
        r"  $f_j = |d_j - D^*_j|$  (severity, duration)",
        r"  $f_j = \min(|m - m^*|,\,12 - |m - m^*|)$ (cyclic month)",
        r"  $f_{k+1} = \sum_j f_j$   (Manhattan norm)",
        "",
        f"Representative event: dur = {ev['duration']} mo, "
        f"min SSI = {ev['severity']:.2f}, "
        f"peak month = {ev['max_severity_date'].month}",
    ]
    for i, line in enumerate(lines):
        ax_c.text(0.02, 0.95 - 0.08 * i, line,
                  transform=ax_c.transAxes, fontsize=8,
                  family="serif", va="top")

    fig.tight_layout()
    return fig


def pd_offset_months(n: int):
    """Return a pandas DateOffset of n months. Defined at module level
    so `fig5_ssi_definitions_example` stays importable when the caller
    has not pre-imported pandas."""
    import pandas as pd
    return pd.DateOffset(months=int(n))
