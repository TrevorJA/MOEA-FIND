"""Diagnostic: Kirsch-wrapper fidelity under uniform random DVs.

Compares four ensembles of 20-year traces at the trace level:

    historical        — overlapping T-year blocks of the historical record.
    kirsch_baseline   — 1000 traces from SynHydro's KirschGenerator.generate()
                        (no DV wrapper). Provides the "natural" Kirsch
                        distribution the wrappers should match under uniform
                        DVs.
    index             — KirschBorgWrapper(mode="index") under U[0,1] DVs.
    residual          — KirschBorgWrapper(mode="residual") under U[0,1] DVs.

Primary questions answered:

    Q1 (invariant / sanity). Under uniform DVs, does each wrapper mode's
       flow distribution match kirsch_baseline? If NO, that is a wrapper-
       level bias and invalidates any downstream mode comparison.
    Q2 (cross-month dependence). Does residual mode (which samples
       quantiles independently per month × site) preserve the within-year
       cross-month correlation that index mode preserves by construction?

A second section (Phase B) converts each trace into drought characteristics
and compares per-mode coverage against the historical T-block cloud, an
in-process 1000-trace kirsch_baseline cloud, and the 10 000-trace exp05
Kirsch library.

Outputs under ``outputs/diag_kirsch_wrapper_fidelity/``:

    figures/
      fig_wrapper_baseline_check.pdf     — wrappers-vs-baseline under U[0,1]
      fig_monthly_marginal_boxplots.pdf  — per-month marginals
      fig_fdc_envelopes.pdf              — FDC 10-90% bands
      fig_acf_envelopes.pdf              — ACF lags 0..24 with 1/12 insets
      fig_cross_month_corr_heatmap.pdf   — 12x12 correlations + diffs
      fig_seasonal_cycle.pdf             — monthly mean and std envelopes
      fig_annual_tails.pdf               — annual-total + annual-min ECDFs
      fig_drought_space_random_2d.pdf    — mode clouds vs historical + Kirsch
      fig_drought_space_random_3d.pdf    — as above in 3D
      fig_coverage_vs_ensemble_size.pdf  — convergence of range / NN-CV
    fidelity_summary.json
    coverage_summary.json

Run (SLURM only — not on the login node):
    sbatch workflows/02_calibration/slurm/wrapper_fidelity.slurm
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.experiment_utils import (  # noqa: E402
    prepare_data,
    compute_historical_ssi_chars,
)
from src.historical_blocks import (  # noqa: E402
    resample_historical_blocks,
    resample_historical_blocks_2d,
    compute_historical_block_chars,
)
from src.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.objectives import (  # noqa: E402
    flows_to_series,
    compute_ssi_drought_characteristics,
)


OUTPUT_SLUG = "diag_kirsch_wrapper_fidelity"

# Ensemble colors/labels. Residual/index picked to match the default
# ablation_dv palette (terracotta and a distinct blue) so downstream figures
# are visually consistent with exp14's SI set.
ENS_COLORS: Dict[str, str] = {
    "historical": "#000000",
    "kirsch_baseline": "#7f7f7f",
    "index": "#2b6cb0",
    "residual": "#c05621",
}
ENS_LABELS: Dict[str, str] = {
    "historical": "historical T-blocks",
    "kirsch_baseline": "Kirsch baseline (no DV)",
    "index": "index wrapper (U[0,1] DVs)",
    "residual": "residual wrapper (U[0,1] DVs)",
}

# Order controls legend order and cloud-overlap z-order in figures.
ENS_ORDER: Tuple[str, ...] = ("historical", "kirsch_baseline", "index", "residual")


# ----------------------------------------------------------------------------
# Generator helpers
# ----------------------------------------------------------------------------

def _baseline_traces(
    kirsch_gen,
    n_traces: int,
    n_years: int,
    seed: int,
) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Generate kirsch_baseline traces via SynHydro's native generate().

    Returns (list of 1D traces, list of 2D traces), water-year aligned to
    match the KirschBorgWrapper output convention.
    """
    ensemble = kirsch_gen.generate(
        n_realizations=n_traces,
        n_years=n_years,
        seed=seed,
    )
    traces_1d: List[np.ndarray] = []
    traces_2d: List[np.ndarray] = []
    for rid in sorted(ensemble.data_by_realization.keys()):
        trace_df = ensemble.data_by_realization[rid]
        trace_1d = trace_df.values.flatten()
        n_yrs = len(trace_1d) // 12
        trace_2d = trace_1d[: n_yrs * 12].reshape(n_yrs, 12)
        # SynHydro emits calendar-year order; roll to water-year (Oct first).
        trace_2d = np.roll(trace_2d, 3, axis=1)
        traces_1d.append(trace_2d.flatten())
        traces_2d.append(trace_2d.copy())
    return traces_1d, traces_2d


def _wrapper_chunk_worker(
    mode: str,
    n_traces_chunk: int,
    n_years: int,
    seed_chunk: int,
    monthly_2d: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Worker for parallel wrapper-trace generation.

    Each worker fits its own KirschGenerator (cheap, < 1s) and then
    generates ``n_traces_chunk`` traces under uniform random DVs. Results
    are returned as stacked numpy arrays so pickling is trivial.
    """
    kirsch_gen = build_kirsch_generator(monthly_2d)
    wrapper = KirschBorgWrapper(kirsch_gen, mode=mode, n_years_out=n_years)
    rng = np.random.default_rng(seed_chunk)
    dvs = rng.uniform(0.0, 1.0, size=(n_traces_chunk, wrapper.n_dvs))
    traces = np.empty((n_traces_chunk, n_years, 12), dtype=float)
    for i in range(n_traces_chunk):
        traces[i] = np.asarray(wrapper.generate(dvs[i]), dtype=float)
    return traces, dvs


def _wrapper_traces(
    kirsch_gen,
    mode: str,
    n_traces: int,
    n_years: int,
    seed: int,
    monthly_2d: np.ndarray,
    workers: int = 1,
) -> Tuple[List[np.ndarray], List[np.ndarray], np.ndarray]:
    """Generate n_traces traces via KirschBorgWrapper under uniform DVs.

    When ``workers > 1`` the N traces are split across a
    :class:`concurrent.futures.ProcessPoolExecutor`; each worker fits its
    own KirschGenerator and generates its share, with distinct seeds
    derived from ``seed`` so the total ensemble is reproducible given
    (seed, workers, n_traces).

    Returns (1D traces list, 2D traces list, DV matrix of shape
    (n_traces, n_dvs)).
    """
    if workers <= 1:
        # Serial path — no pool overhead.
        wrapper = KirschBorgWrapper(kirsch_gen, mode=mode, n_years_out=n_years)
        rng = np.random.default_rng(seed)
        dvs = rng.uniform(0.0, 1.0, size=(n_traces, wrapper.n_dvs))
        traces_1d: List[np.ndarray] = []
        traces_2d: List[np.ndarray] = []
        for i in range(n_traces):
            trace_2d = np.asarray(wrapper.generate(dvs[i]), dtype=float)
            traces_2d.append(trace_2d)
            traces_1d.append(trace_2d.flatten())
        return traces_1d, traces_2d, dvs

    # Parallel path — chunk trace generation.
    chunk_sizes = [n_traces // workers] * workers
    for i in range(n_traces % workers):
        chunk_sizes[i] += 1
    # Distinct per-worker seeds so the full draw is reproducible.
    seeds = [seed + 1_000_003 * i for i in range(workers)]

    traces_chunks: List[np.ndarray] = [None] * workers  # type: ignore
    dvs_chunks: List[np.ndarray] = [None] * workers  # type: ignore
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(
                _wrapper_chunk_worker, mode, chunk_sizes[i], n_years,
                seeds[i], monthly_2d,
            ): i
            for i in range(workers) if chunk_sizes[i] > 0
        }
        for fut in as_completed(futures):
            i = futures[fut]
            tr, dv = fut.result()
            traces_chunks[i] = tr
            dvs_chunks[i] = dv

    traces_all = np.concatenate([c for c in traces_chunks if c is not None], axis=0)
    dvs_all = np.concatenate([c for c in dvs_chunks if c is not None], axis=0)
    traces_1d = [traces_all[i].flatten() for i in range(traces_all.shape[0])]
    traces_2d = [traces_all[i] for i in range(traces_all.shape[0])]
    return traces_1d, traces_2d, dvs_all


# ----------------------------------------------------------------------------
# Statistical primitives
# ----------------------------------------------------------------------------

def _acf(x: np.ndarray, max_lag: int) -> np.ndarray:
    n = len(x)
    mean = np.mean(x)
    var = np.var(x)
    if var < 1e-12:
        return np.zeros(max_lag + 1)
    out = np.zeros(max_lag + 1)
    for lag in range(max_lag + 1):
        out[lag] = np.mean((x[: n - lag] - mean) * (x[lag:] - mean)) / var
    return out


def _envelope(stack: np.ndarray, pcts=(10, 50, 90)) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    lo, mid, hi = (np.percentile(stack, p, axis=0) for p in pcts)
    return lo, mid, hi


def _fdc(flows: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    sorted_flows = np.sort(flows)[::-1]
    exceedance = np.arange(1, len(sorted_flows) + 1) / (len(sorted_flows) + 1)
    return exceedance, sorted_flows


def _stack_fdc(traces: List[np.ndarray], grid: np.ndarray) -> np.ndarray:
    stacked = np.empty((len(traces), len(grid)))
    for i, tr in enumerate(traces):
        e, f = _fdc(tr)
        stacked[i] = np.interp(grid, e, f)
    return stacked


def _cross_month_corr(traces_2d: List[np.ndarray]) -> np.ndarray:
    """Average 12x12 correlation matrix across the ensemble.

    For each (T_years, 12) trace we compute the 12x12 Pearson correlation
    across years (i.e. correlation between calendar months). The result is
    averaged over the ensemble; this is the standard "cross-month
    correlation" summary used in Kirsch-Nowak evaluations.
    """
    mats = []
    for t in traces_2d:
        if t.shape[0] < 2:
            continue
        c = np.corrcoef(t, rowvar=False)
        mats.append(c)
    return np.mean(np.array(mats), axis=0) if mats else np.zeros((12, 12))


def _annual_totals(trace_2d: np.ndarray) -> np.ndarray:
    return trace_2d.sum(axis=1)


def _annual_min_monthly(trace_2d: np.ndarray) -> np.ndarray:
    return trace_2d.min(axis=1)


def _ks_two_sample(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sample KS statistic (no p-value, just the D statistic)."""
    from scipy.stats import ks_2samp
    return float(ks_2samp(a, b).statistic)


def _ad_two_sample(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sample Anderson-Darling statistic."""
    from scipy.stats import anderson_ksamp
    try:
        # anderson_ksamp is sensitive to small samples; guard it.
        return float(anderson_ksamp([a, b]).statistic)
    except Exception:
        return float("nan")


# ----------------------------------------------------------------------------
# Figure builders
# ----------------------------------------------------------------------------

def _fig_baseline_check(
    ensembles_1d: Dict[str, List[np.ndarray]],
    ensembles_2d: Dict[str, List[np.ndarray]],
    fig_path: Path,
) -> None:
    """fig_wrapper_baseline_check.pdf: wrappers vs baseline under U[0,1] DVs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.plotting.style import apply_style
    apply_style()

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2))
    ax_fdc, ax_lag1, ax_ann = axes

    # FDC bands (log y)
    any_trace = next(iter(ensembles_1d.values()))[0]
    grid = np.arange(1, len(any_trace) + 1) / (len(any_trace) + 1)
    grid_pct = grid * 100
    for ens in ENS_ORDER:
        traces = ensembles_1d[ens]
        if not traces:
            continue
        stack = _stack_fdc(traces, grid)
        lo, mid, hi = _envelope(stack)
        ax_fdc.fill_between(
            grid_pct, lo, hi, alpha=0.18, color=ENS_COLORS[ens],
            label=f"{ENS_LABELS[ens]} (n={len(traces)})",
        )
        ax_fdc.semilogy(grid_pct, mid, color=ENS_COLORS[ens], linewidth=1.6)
    ax_fdc.set_xlabel("Exceedance probability (%)")
    ax_fdc.set_ylabel("Flow (cfs, log)")
    ax_fdc.set_title("(a) FDC 10–90% envelope")
    ax_fdc.legend(fontsize=7, loc="upper right", framealpha=0.85)

    # Lag-1 ACF distribution
    for ens in ENS_ORDER:
        traces = ensembles_1d[ens]
        if not traces:
            continue
        lag1 = np.array([_acf(t, 1)[1] for t in traces])
        ax_lag1.hist(
            lag1, bins=40, density=True, alpha=0.45,
            color=ENS_COLORS[ens], label=f"{ENS_LABELS[ens]} (med={np.median(lag1):.2f})",
        )
    ax_lag1.set_xlabel("Lag-1 autocorrelation")
    ax_lag1.set_ylabel("Density")
    ax_lag1.set_title("(b) Lag-1 ACF distribution")
    ax_lag1.legend(fontsize=7, framealpha=0.85)

    # Annual total distribution
    for ens in ENS_ORDER:
        traces = ensembles_2d[ens]
        if not traces:
            continue
        totals = np.concatenate([_annual_totals(t) for t in traces])
        ax_ann.hist(
            totals, bins=60, density=True, alpha=0.4,
            color=ENS_COLORS[ens], label=ENS_LABELS[ens],
        )
    ax_ann.set_xlabel("Annual total flow (cfs)")
    ax_ann.set_ylabel("Density")
    ax_ann.set_title("(c) Annual total distribution")
    ax_ann.legend(fontsize=7, framealpha=0.85)

    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_monthly_marginal_boxplots(
    ensembles_2d: Dict[str, List[np.ndarray]],
    fig_path: Path,
) -> None:
    """fig_monthly_marginal_boxplots.pdf."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.plotting.style import apply_style, WATER_YEAR_MONTHS
    apply_style()

    fig, ax = plt.subplots(figsize=(12.0, 4.5))
    months = np.arange(12)
    width = 0.2  # box width per ensemble
    active = [e for e in ENS_ORDER if ensembles_2d.get(e)]
    n_active = len(active)
    start = -(n_active - 1) / 2.0

    for i, ens in enumerate(active):
        data = np.concatenate(ensembles_2d[ens], axis=0)  # (n_years_total, 12)
        positions = months + (start + i) * width
        box_data = [data[:, m] for m in range(12)]
        bp = ax.boxplot(
            box_data, positions=positions, widths=width * 0.9,
            patch_artist=True, showfliers=False,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(ENS_COLORS[ens])
            patch.set_alpha(0.45)
            patch.set_edgecolor(ENS_COLORS[ens])
        for whisker in bp["whiskers"]:
            whisker.set_color(ENS_COLORS[ens])
        for cap in bp["caps"]:
            cap.set_color(ENS_COLORS[ens])
        for median in bp["medians"]:
            median.set_color("black")
            median.set_linewidth(1.2)

    ax.set_xticks(months)
    ax.set_xticklabels(WATER_YEAR_MONTHS)
    ax.set_ylabel("Monthly flow (cfs)")
    ax.set_yscale("log")
    ax.set_title("Per-month marginals by ensemble (log scale)")
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=ENS_COLORS[e], alpha=0.45,
                              edgecolor=ENS_COLORS[e], label=ENS_LABELS[e])
               for e in active]
    ax.legend(handles=handles, fontsize=8, loc="upper right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_fdc_envelopes(
    ensembles_1d: Dict[str, List[np.ndarray]],
    fig_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.plotting.style import apply_style
    apply_style()

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    any_trace = next(iter(ensembles_1d.values()))[0]
    grid = np.arange(1, len(any_trace) + 1) / (len(any_trace) + 1)
    grid_pct = grid * 100
    for ens in ENS_ORDER:
        traces = ensembles_1d[ens]
        if not traces:
            continue
        stack = _stack_fdc(traces, grid)
        lo, mid, hi = _envelope(stack)
        ax.fill_between(
            grid_pct, lo, hi, alpha=0.2, color=ENS_COLORS[ens],
            label=f"{ENS_LABELS[ens]} 10–90% (n={len(traces)})",
        )
        ax.semilogy(grid_pct, mid, color=ENS_COLORS[ens], linewidth=1.6)
    ax.set_xlabel("Exceedance probability (%)")
    ax.set_ylabel("Flow (cfs, log)")
    ax.set_title("Flow duration curve envelopes")
    ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_acf_envelopes(
    ensembles_1d: Dict[str, List[np.ndarray]],
    fig_path: Path,
    max_lag: int = 24,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.plotting.style import apply_style
    apply_style()

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    lags = np.arange(max_lag + 1)
    per_ens_lag_arrays: Dict[str, np.ndarray] = {}
    for ens in ENS_ORDER:
        traces = ensembles_1d[ens]
        if not traces:
            continue
        acfs = np.array([_acf(t, max_lag) for t in traces])
        per_ens_lag_arrays[ens] = acfs
        lo, mid, hi = _envelope(acfs)
        ax.fill_between(
            lags, lo, hi, alpha=0.18, color=ENS_COLORS[ens],
            label=f"{ENS_LABELS[ens]} 10–90%",
        )
        ax.plot(lags, mid, color=ENS_COLORS[ens], linewidth=1.6)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle=":")
    ax.set_xlabel("Lag (months)")
    ax.set_ylabel("Autocorrelation")
    ax.set_title(f"ACF envelopes (lags 0–{max_lag})")
    ax.legend(fontsize=8, framealpha=0.9)

    # Inset: lag-1 and lag-12 bar comparison.
    ins = ax.inset_axes([0.55, 0.50, 0.42, 0.42])
    positions = np.arange(2)
    width = 0.18
    active = [e for e in ENS_ORDER if e in per_ens_lag_arrays]
    for i, ens in enumerate(active):
        arr = per_ens_lag_arrays[ens]
        medians = [np.median(arr[:, 1]), np.median(arr[:, 12])]
        lows = [np.percentile(arr[:, 1], 10), np.percentile(arr[:, 12], 10)]
        his = [np.percentile(arr[:, 1], 90), np.percentile(arr[:, 12], 90)]
        offsets = (i - (len(active) - 1) / 2.0) * width
        ins.bar(
            positions + offsets, medians, width=width * 0.95,
            color=ENS_COLORS[ens], alpha=0.8,
        )
        ins.errorbar(
            positions + offsets, medians,
            yerr=[np.array(medians) - np.array(lows), np.array(his) - np.array(medians)],
            fmt="none", ecolor="black", linewidth=0.8, capsize=2,
        )
    ins.set_xticks(positions)
    ins.set_xticklabels(["lag 1", "lag 12"], fontsize=7)
    ins.tick_params(axis="y", labelsize=7)
    ins.set_title("median ± 10/90%", fontsize=7)

    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_cross_month_corr_heatmap(
    ensembles_2d: Dict[str, List[np.ndarray]],
    fig_path: Path,
) -> None:
    """One heatmap per ensemble; second row shows (ens - historical) diffs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.plotting.style import apply_style, WATER_YEAR_MONTHS
    apply_style()

    active = [e for e in ENS_ORDER if ensembles_2d.get(e)]
    mats = {e: _cross_month_corr(ensembles_2d[e]) for e in active}

    n = len(active)
    fig, axes = plt.subplots(2, n, figsize=(3.0 * n + 1.0, 6.3))
    if n == 1:
        axes = axes.reshape(2, 1)

    # Top row: raw 12x12 correlation
    for j, ens in enumerate(active):
        ax = axes[0, j]
        im = ax.imshow(mats[ens], vmin=-0.5, vmax=1.0, cmap="RdBu_r", aspect="equal")
        ax.set_title(ENS_LABELS[ens], fontsize=9)
        ax.set_xticks(range(12))
        ax.set_xticklabels(WATER_YEAR_MONTHS, fontsize=6, rotation=90)
        ax.set_yticks(range(12))
        ax.set_yticklabels(WATER_YEAR_MONTHS, fontsize=6)
    cbar = fig.colorbar(im, ax=axes[0, :].tolist(), shrink=0.75)
    cbar.set_label("cross-month Pearson r", fontsize=8)

    # Bottom row: difference vs historical (first active)
    ref_key = "historical" if "historical" in mats else active[0]
    ref = mats[ref_key]
    for j, ens in enumerate(active):
        ax = axes[1, j]
        diff = mats[ens] - ref
        im2 = ax.imshow(diff, vmin=-0.4, vmax=0.4, cmap="RdBu_r", aspect="equal")
        frob = np.linalg.norm(diff, ord="fro")
        ax.set_title(f"{ens} − {ref_key}  (‖·‖_F = {frob:.2f})", fontsize=8)
        ax.set_xticks(range(12))
        ax.set_xticklabels(WATER_YEAR_MONTHS, fontsize=6, rotation=90)
        ax.set_yticks(range(12))
        ax.set_yticklabels(WATER_YEAR_MONTHS, fontsize=6)
    cbar2 = fig.colorbar(im2, ax=axes[1, :].tolist(), shrink=0.75)
    cbar2.set_label("Δ cross-month r", fontsize=8)

    fig.suptitle("Cross-month correlation (ensemble-averaged)", fontsize=11)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_seasonal_cycle(
    ensembles_2d: Dict[str, List[np.ndarray]],
    fig_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.plotting.style import apply_style, WATER_YEAR_MONTHS
    apply_style()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.5))
    months = np.arange(12)
    for ens in ENS_ORDER:
        traces = ensembles_2d[ens]
        if not traces:
            continue
        means = np.array([t.mean(axis=0) for t in traces])
        stds = np.array([t.std(axis=0, ddof=1) for t in traces])
        for ax, stat in [(ax1, means), (ax2, stds)]:
            lo, mid, hi = _envelope(stat)
            ax.fill_between(months, lo, hi, alpha=0.18, color=ENS_COLORS[ens])
            ax.plot(months, mid, color=ENS_COLORS[ens], linewidth=1.6,
                    label=ENS_LABELS[ens])
    for ax, ttl, ylbl in [(ax1, "Seasonal mean", "Mean monthly flow (cfs)"),
                          (ax2, "Seasonal std", "Std monthly flow (cfs)")]:
        ax.set_xticks(months)
        ax.set_xticklabels(WATER_YEAR_MONTHS, rotation=45)
        ax.set_ylabel(ylbl)
        ax.set_title(ttl)
    ax1.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_annual_tails(
    ensembles_2d: Dict[str, List[np.ndarray]],
    fig_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.plotting.style import apply_style
    apply_style()

    fig, (ax_tot, ax_min) = plt.subplots(1, 2, figsize=(11.0, 4.5))
    for ens in ENS_ORDER:
        traces = ensembles_2d[ens]
        if not traces:
            continue
        totals = np.concatenate([_annual_totals(t) for t in traces])
        mins = np.concatenate([_annual_min_monthly(t) for t in traces])
        for ax, data in [(ax_tot, totals), (ax_min, mins)]:
            srt = np.sort(data)
            cdf = np.arange(1, len(srt) + 1) / (len(srt) + 1)
            ax.plot(srt, cdf, color=ENS_COLORS[ens], linewidth=1.6,
                    label=ENS_LABELS[ens])
    ax_tot.set_xlabel("Annual total flow (cfs)")
    ax_tot.set_title("Annual total ECDF")
    ax_tot.set_ylabel("Cumulative probability")
    ax_min.set_xlabel("Annual minimum monthly flow (cfs)")
    ax_min.set_title("Annual min-month ECDF  (drought tail)")
    ax_min.set_xscale("log")
    for ax in (ax_tot, ax_min):
        ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Phase B — drought-space coverage
# ----------------------------------------------------------------------------

def _drought_chars_chunk_worker(
    traces_1d: List[np.ndarray],
    ssi_calc,
    objective_keys: Tuple[str, ...],
) -> np.ndarray:
    """Worker: compute drought-char vectors for a list of traces.

    Passes the raw monthly flows through so trace-level extras like
    ``q10_flow_neg`` are available in the chars dict for metric
    extraction.
    """
    out = np.zeros((len(traces_1d), len(objective_keys)))
    for i, t in enumerate(traces_1d):
        series = flows_to_series(t, start_date="2100-01-01")
        ssi = ssi_calc.transform(series)
        chars = compute_ssi_drought_characteristics(ssi, monthly_flows=t)
        for j, k in enumerate(objective_keys):
            out[i, j] = float(chars.get(k, np.nan))
    return out


def _drought_chars_for_traces(
    traces_1d: List[np.ndarray],
    ssi_calc,
    objective_keys: Tuple[str, ...],
    workers: int = 1,
) -> np.ndarray:
    """Compute drought-characteristic vectors per trace (n_traces, K).

    When ``workers > 1`` traces are chunked and distributed across a
    ProcessPoolExecutor. The fitted ``ssi_calc`` is pickled into each
    worker (the SynHydroSSI object is picklable because its monthly-
    distribution parameters are plain numpy/scalar state).
    """
    n = len(traces_1d)
    if workers <= 1 or n < workers * 4:
        return _drought_chars_chunk_worker(traces_1d, ssi_calc, objective_keys)

    chunk_size = n // workers
    chunks = [
        traces_1d[i * chunk_size : (i + 1) * chunk_size]
        for i in range(workers - 1)
    ]
    chunks.append(traces_1d[(workers - 1) * chunk_size :])

    parts: List[np.ndarray] = [None] * workers  # type: ignore
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(
                _drought_chars_chunk_worker, chunks[i], ssi_calc, objective_keys,
            ): i
            for i in range(workers) if len(chunks[i]) > 0
        }
        for fut in as_completed(futures):
            parts[futures[fut]] = fut.result()
    return np.concatenate([p for p in parts if p is not None], axis=0)


def _hull_volume(points: np.ndarray) -> float:
    """Convex hull volume in R^d (d in {2, 3}). Returns NaN on degenerate input."""
    try:
        from scipy.spatial import ConvexHull
    except ImportError:
        return float("nan")
    pts = points[~np.isnan(points).any(axis=1)]
    if len(pts) < pts.shape[1] + 1:
        return float("nan")
    try:
        h = ConvexHull(pts)
        return float(h.volume)
    except Exception:
        return float("nan")


def _fraction_inside_hull(points: np.ndarray, hull_points: np.ndarray) -> float:
    """Fraction of points that lie inside the convex hull of hull_points."""
    try:
        from scipy.spatial import Delaunay
    except ImportError:
        return float("nan")
    pts = points[~np.isnan(points).any(axis=1)]
    hpts = hull_points[~np.isnan(hull_points).any(axis=1)]
    if len(hpts) < hpts.shape[1] + 1 or len(pts) == 0:
        return float("nan")
    try:
        tri = Delaunay(hpts)
        inside = tri.find_simplex(pts) >= 0
        return float(np.mean(inside))
    except Exception:
        return float("nan")


def _fig_drought_space_2d(
    mode_clouds: Dict[str, np.ndarray],
    hist_cloud: np.ndarray,
    in_proc_baseline_cloud: np.ndarray,
    exp05_cloud: np.ndarray,
    objective_labels: Tuple[str, str],
    fig_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.plotting.style import apply_style
    apply_style()

    fig, ax = plt.subplots(figsize=(7.5, 6.2))

    if exp05_cloud is not None and exp05_cloud.shape[0] > 0:
        ax.scatter(
            exp05_cloud[:, 0], exp05_cloud[:, 1],
            color="#d9d9d9", s=5, alpha=0.35,
            label=f"exp05 Kirsch library (n={len(exp05_cloud)})",
        )
    if hist_cloud is not None and hist_cloud.shape[0] > 0:
        ax.scatter(
            hist_cloud[:, 0], hist_cloud[:, 1],
            color=ENS_COLORS["historical"], s=30, marker="x",
            label=f"historical T-blocks (n={len(hist_cloud)})",
        )
    if in_proc_baseline_cloud is not None and in_proc_baseline_cloud.shape[0] > 0:
        ax.scatter(
            in_proc_baseline_cloud[:, 0], in_proc_baseline_cloud[:, 1],
            color=ENS_COLORS["kirsch_baseline"], s=8, alpha=0.45,
            label=f"Kirsch baseline (n={len(in_proc_baseline_cloud)})",
        )
    for mode in ("index", "residual"):
        cloud = mode_clouds.get(mode)
        if cloud is None or cloud.shape[0] == 0:
            continue
        ax.scatter(
            cloud[:, 0], cloud[:, 1],
            color=ENS_COLORS[mode], s=10, alpha=0.55, marker="o",
            edgecolor="none",
            label=f"{ENS_LABELS[mode]} (n={len(cloud)})",
        )
    ax.set_xlabel(objective_labels[0])
    ax.set_ylabel(objective_labels[1])
    ax.set_title("Drought-space coverage under random DVs")
    ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_drought_space_3d(
    mode_clouds: Dict[str, np.ndarray],
    hist_cloud: np.ndarray,
    in_proc_baseline_cloud: np.ndarray,
    exp05_cloud: np.ndarray,
    objective_labels: Tuple[str, str, str],
    fig_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.plotting.style import apply_style
    apply_style()

    fig = plt.figure(figsize=(8.0, 6.5))
    ax = fig.add_subplot(111, projection="3d")

    if exp05_cloud is not None and exp05_cloud.shape[0] > 0 and exp05_cloud.shape[1] >= 3:
        ax.scatter(
            exp05_cloud[:, 0], exp05_cloud[:, 1], exp05_cloud[:, 2],
            color="#d9d9d9", s=3, alpha=0.3, label=f"exp05 Kirsch library (n={len(exp05_cloud)})",
        )
    if hist_cloud is not None and hist_cloud.shape[0] > 0:
        ax.scatter(
            hist_cloud[:, 0], hist_cloud[:, 1], hist_cloud[:, 2],
            color=ENS_COLORS["historical"], s=20, marker="x",
            label=f"historical T-blocks (n={len(hist_cloud)})",
        )
    if in_proc_baseline_cloud is not None and in_proc_baseline_cloud.shape[0] > 0:
        ax.scatter(
            in_proc_baseline_cloud[:, 0], in_proc_baseline_cloud[:, 1], in_proc_baseline_cloud[:, 2],
            color=ENS_COLORS["kirsch_baseline"], s=6, alpha=0.4,
            label=f"Kirsch baseline (n={len(in_proc_baseline_cloud)})",
        )
    for mode in ("index", "residual"):
        cloud = mode_clouds.get(mode)
        if cloud is None or cloud.shape[0] == 0 or cloud.shape[1] < 3:
            continue
        ax.scatter(
            cloud[:, 0], cloud[:, 1], cloud[:, 2],
            color=ENS_COLORS[mode], s=8, alpha=0.5,
            label=f"{ENS_LABELS[mode]} (n={len(cloud)})",
        )
    ax.set_xlabel(objective_labels[0])
    ax.set_ylabel(objective_labels[1])
    ax.set_zlabel(objective_labels[2])
    ax.set_title("Drought-space 3D coverage (random DVs)")
    ax.legend(fontsize=7, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_coverage_vs_n(
    mode_clouds: Dict[str, np.ndarray],
    exp05_cloud: np.ndarray,
    objective_labels: Tuple[str, ...],
    fig_path: Path,
    ns: Tuple[int, ...] = (50, 100, 200, 500, 1000),
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.plotting.style import apply_style
    apply_style()

    n_obj = len(objective_labels)
    fig, axes = plt.subplots(1, n_obj, figsize=(4.0 * n_obj, 4.2), sharex=True)
    if n_obj == 1:
        axes = [axes]

    for j, lbl in enumerate(objective_labels):
        ax = axes[j]
        for mode in ("index", "residual"):
            cloud = mode_clouds.get(mode)
            if cloud is None or cloud.shape[0] == 0 or cloud.shape[1] <= j:
                continue
            max_curve = []
            for n in ns:
                n_use = min(n, cloud.shape[0])
                max_curve.append(np.max(cloud[:n_use, j]))
            ax.plot(ns, max_curve, "o-", color=ENS_COLORS[mode],
                    label=ENS_LABELS[mode], linewidth=1.6, markersize=5)
        if exp05_cloud is not None and exp05_cloud.shape[0] > 0 and exp05_cloud.shape[1] > j:
            ax.axhline(
                np.max(exp05_cloud[:, j]), color="#7f7f7f",
                linestyle="--", linewidth=1.2,
                label=f"exp05 library max (n={len(exp05_cloud)})",
            )
        ax.set_xlabel("Ensemble size N")
        ax.set_ylabel(f"max {lbl}")
        ax.set_xscale("log")
        ax.set_title(lbl)
        ax.legend(fontsize=7, framealpha=0.9)
    fig.suptitle("Per-axis range coverage vs random-DV ensemble size", fontsize=11)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-traces", type=int, default=1000,
                   help="Traces per ensemble (wrappers and kirsch_baseline).")
    p.add_argument("--n-years", type=int, default=20,
                   help="Trace length in water years.")
    p.add_argument("--ssi", type=int, default=3, choices=[1, 3, 6, 12])
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    p.add_argument("--kirsch-library-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / "exp05_kirsch_library",
                   help="Location of the exp05 Kirsch library for the Phase B "
                        "extended reference cloud.")
    p.add_argument("--skip-phase-b", action="store_true",
                   help="Skip the drought-space coverage section (Phase B).")
    p.add_argument("--metric-set", default="primary",
                   help="Drought metric set: a preset name from "
                        "src.drought_metrics.PRESETS or a single metric name.")
    p.add_argument("--workers", type=int,
                   default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")),
                   help="Parallel worker count for trace generation and SSI "
                        "characterisation. Defaults to $SLURM_CPUS_PER_TASK.")
    args = p.parse_args()

    out_dir = args.output_dir
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    print(f"[diag_wrapper_fidelity] output_dir={out_dir}")

    # ------------------------------------------------------------------
    # Data
    # ------------------------------------------------------------------
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)
    print(f"[diag_wrapper_fidelity] historical: {monthly_2d.shape[0]} water years")

    hist_blocks_1d = resample_historical_blocks(monthly_1d, T_years=args.n_years, stride=1)
    hist_blocks_2d = resample_historical_blocks_2d(monthly_2d, T_years=args.n_years, stride=1)
    print(f"[diag_wrapper_fidelity] historical T-blocks: n={len(hist_blocks_1d)}")

    # ------------------------------------------------------------------
    # Generator + ensembles
    # ------------------------------------------------------------------
    print(f"[diag_wrapper_fidelity] fitting Kirsch generator ...")
    kirsch_gen = build_kirsch_generator(monthly_2d)

    print(f"[diag_wrapper_fidelity] generating kirsch_baseline (n={args.n_traces}) ...")
    baseline_1d, baseline_2d = _baseline_traces(
        kirsch_gen, args.n_traces, args.n_years, args.seed,
    )

    workers = max(1, int(args.workers))
    print(f"[diag_wrapper_fidelity] workers={workers}")

    print(f"[diag_wrapper_fidelity] generating index-wrapper traces ...")
    index_1d, index_2d, index_dvs = _wrapper_traces(
        kirsch_gen, "index", args.n_traces, args.n_years, args.seed + 1,
        monthly_2d=monthly_2d, workers=workers,
    )

    print(f"[diag_wrapper_fidelity] generating residual-wrapper traces ...")
    resid_1d, resid_2d, resid_dvs = _wrapper_traces(
        kirsch_gen, "residual", args.n_traces, args.n_years, args.seed + 2,
        monthly_2d=monthly_2d, workers=workers,
    )

    ensembles_1d: Dict[str, List[np.ndarray]] = {
        "historical": hist_blocks_1d,
        "kirsch_baseline": baseline_1d,
        "index": index_1d,
        "residual": resid_1d,
    }
    ensembles_2d: Dict[str, List[np.ndarray]] = {
        "historical": hist_blocks_2d,
        "kirsch_baseline": baseline_2d,
        "index": index_2d,
        "residual": resid_2d,
    }

    # ------------------------------------------------------------------
    # Phase A figures
    # ------------------------------------------------------------------
    print(f"[diag_wrapper_fidelity] Phase A figures ...")
    _fig_baseline_check(ensembles_1d, ensembles_2d, fig_dir / "fig_wrapper_baseline_check.pdf")
    _fig_monthly_marginal_boxplots(ensembles_2d, fig_dir / "fig_monthly_marginal_boxplots.pdf")
    _fig_fdc_envelopes(ensembles_1d, fig_dir / "fig_fdc_envelopes.pdf")
    _fig_acf_envelopes(ensembles_1d, fig_dir / "fig_acf_envelopes.pdf")
    _fig_cross_month_corr_heatmap(ensembles_2d, fig_dir / "fig_cross_month_corr_heatmap.pdf")
    _fig_seasonal_cycle(ensembles_2d, fig_dir / "fig_seasonal_cycle.pdf")
    _fig_annual_tails(ensembles_2d, fig_dir / "fig_annual_tails.pdf")

    # ------------------------------------------------------------------
    # Fidelity summary JSON
    # ------------------------------------------------------------------
    hist_flat_by_month = [
        np.concatenate([b[:, m] for b in hist_blocks_2d]) for m in range(12)
    ]
    summary: Dict = {
        "config": {
            "n_traces": args.n_traces, "n_years": args.n_years,
            "seed": args.seed,
        },
        "n_hist_blocks": len(hist_blocks_1d),
    }
    summary_by_ens: Dict[str, Dict] = {}
    for ens in ENS_ORDER:
        traces_1d = ensembles_1d[ens]
        traces_2d = ensembles_2d[ens]
        if not traces_1d:
            continue
        annual_totals = np.concatenate([_annual_totals(t) for t in traces_2d])
        annual_mins = np.concatenate([_annual_min_monthly(t) for t in traces_2d])
        acfs = np.array([_acf(t, 12) for t in traces_1d])
        cxm = _cross_month_corr(traces_2d)
        flat_by_month = [np.concatenate([t[:, m] for t in traces_2d]) for m in range(12)]
        per_month_ks = [_ks_two_sample(flat_by_month[m], hist_flat_by_month[m])
                        for m in range(12)]
        per_month_ad = [_ad_two_sample(flat_by_month[m], hist_flat_by_month[m])
                        for m in range(12)]
        summary_by_ens[ens] = {
            "n_traces": len(traces_1d),
            "annual_mean": float(annual_totals.mean()),
            "annual_cv": float(annual_totals.std(ddof=1) / annual_totals.mean()),
            "annual_min_median": float(np.median(annual_mins)),
            "lag1_ac_median": float(np.median(acfs[:, 1])),
            "lag12_ac_median": float(np.median(acfs[:, 12])),
            "cross_month_frob_vs_hist": float(np.linalg.norm(
                cxm - _cross_month_corr(hist_blocks_2d), ord="fro",
            )),
            "per_month_ks_vs_hist": per_month_ks,
            "per_month_ad_vs_hist": per_month_ad,
            "per_month_ks_max": float(max(per_month_ks)),
            "per_month_ad_max": float(max(per_month_ad)),
        }
    summary["ensembles"] = summary_by_ens
    (out_dir / "fidelity_summary.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )
    print(f"[diag_wrapper_fidelity] wrote {out_dir / 'fidelity_summary.json'}")

    # ------------------------------------------------------------------
    # Phase B — drought-space coverage
    # ------------------------------------------------------------------
    if args.skip_phase_b:
        print(f"[diag_wrapper_fidelity] Phase B skipped.")
        return

    print(f"[diag_wrapper_fidelity] Phase B: drought characterisation ...")
    from src.drought_metrics import metric_labels, metric_names, resolve_metric_set
    metric_set = resolve_metric_set(args.metric_set)
    objective_keys = metric_names(metric_set)
    objective_labels = tuple(
        f"{m.label} ({m.units})" for m in metric_set
    )
    print(f"[diag_wrapper_fidelity] metric set: {args.metric_set} → {objective_keys}")
    _, ssi_calc, _ = compute_historical_ssi_chars(monthly_1d, args.ssi)

    mode_clouds: Dict[str, np.ndarray] = {
        "index": _drought_chars_for_traces(index_1d, ssi_calc, objective_keys,
                                           workers=workers),
        "residual": _drought_chars_for_traces(resid_1d, ssi_calc, objective_keys,
                                              workers=workers),
    }
    in_proc_baseline_cloud = _drought_chars_for_traces(
        baseline_1d, ssi_calc, objective_keys, workers=workers,
    )
    hist_cloud = compute_historical_block_chars(
        monthly_1d, T_years=args.n_years, ssi_calc=ssi_calc,
        objective_keys=objective_keys, stride=1,
    )

    exp05_cloud: np.ndarray = np.zeros((0, len(objective_keys)))
    exp05_path = args.kirsch_library_dir / "characteristics.npz"
    if exp05_path.exists():
        try:
            klib = np.load(exp05_path, allow_pickle=True)
            k_all_keys = [str(k) for k in klib["all_keys"]]
            k_all_vals = klib["all_values"]
            obj_indices: List[int] = []
            for ok in objective_keys:
                if ok in k_all_keys:
                    obj_indices.append(k_all_keys.index(ok))
                else:
                    obj_indices.append(-1)
            if all(j >= 0 for j in obj_indices):
                exp05_cloud = k_all_vals[:, obj_indices].astype(float)
                print(f"[diag_wrapper_fidelity] loaded exp05 library: "
                      f"n={exp05_cloud.shape[0]}")
            else:
                missing = [ok for ok, j in zip(objective_keys, obj_indices) if j < 0]
                print(f"[diag_wrapper_fidelity] exp05 library missing keys {missing}; "
                      f"proceeding without exp05 reference.")
        except Exception as exc:
            print(f"[diag_wrapper_fidelity] WARNING loading exp05 library: {exc}")
    else:
        print(f"[diag_wrapper_fidelity] exp05 library not found at {exp05_path}; "
              "drought-space plots will omit exp05 reference cloud.")

    print(f"[diag_wrapper_fidelity] Phase B figures ...")
    _fig_drought_space_2d(
        mode_clouds={m: c[:, :2] for m, c in mode_clouds.items()},
        hist_cloud=hist_cloud[:, :2] if hist_cloud.shape[0] > 0 else hist_cloud,
        in_proc_baseline_cloud=in_proc_baseline_cloud[:, :2],
        exp05_cloud=exp05_cloud[:, :2] if exp05_cloud.shape[0] > 0 else exp05_cloud,
        objective_labels=objective_labels[:2],
        fig_path=fig_dir / "fig_drought_space_random_2d.pdf",
    )
    _fig_drought_space_3d(
        mode_clouds=mode_clouds,
        hist_cloud=hist_cloud,
        in_proc_baseline_cloud=in_proc_baseline_cloud,
        exp05_cloud=exp05_cloud,
        objective_labels=objective_labels,
        fig_path=fig_dir / "fig_drought_space_random_3d.pdf",
    )
    _fig_coverage_vs_n(
        mode_clouds=mode_clouds,
        exp05_cloud=exp05_cloud,
        objective_labels=objective_labels,
        fig_path=fig_dir / "fig_coverage_vs_ensemble_size.pdf",
    )

    # ------------------------------------------------------------------
    # Coverage summary
    # ------------------------------------------------------------------
    cov: Dict = {"objective_keys": list(objective_keys)}
    per_mode_cov: Dict[str, Dict] = {}
    ref_hull_2d = exp05_cloud[:, :2] if exp05_cloud.shape[0] > 0 else hist_cloud[:, :2]
    ref_hull_3d = exp05_cloud if exp05_cloud.shape[0] > 0 else hist_cloud
    for mode in ("index", "residual"):
        cloud = mode_clouds[mode]
        if cloud.size == 0:
            continue
        per_mode_cov[mode] = {
            "n": int(cloud.shape[0]),
            "per_axis_min": [float(x) for x in np.nanmin(cloud, axis=0)],
            "per_axis_max": [float(x) for x in np.nanmax(cloud, axis=0)],
            "hull_volume_2d": _hull_volume(cloud[:, :2]),
            "hull_volume_3d": _hull_volume(cloud),
            "frac_inside_exp05_hull_2d": _fraction_inside_hull(cloud[:, :2], ref_hull_2d),
            "frac_inside_exp05_hull_3d": _fraction_inside_hull(cloud, ref_hull_3d),
        }
    cov["per_mode"] = per_mode_cov
    cov["in_proc_baseline_n"] = int(in_proc_baseline_cloud.shape[0])
    cov["exp05_library_n"] = int(exp05_cloud.shape[0])
    cov["historical_block_n"] = int(hist_cloud.shape[0])
    (out_dir / "coverage_summary.json").write_text(
        json.dumps(cov, indent=2, default=str)
    )
    print(f"[diag_wrapper_fidelity] wrote {out_dir / 'coverage_summary.json'}")
    print(f"[diag_wrapper_fidelity] done.")


if __name__ == "__main__":
    main()
