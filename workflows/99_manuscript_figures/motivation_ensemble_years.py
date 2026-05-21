"""Motivation figure: how many baseline Kirsch ensemble years are needed
to reach 1960s-level drought severity and magnitude.

Two panels, x-axis = number of scenario years (log scale):
  left  — maximum drought duration  (months)             vs ensemble years
  right — maximum drought magnitude (|SSI|-month deficit) vs ensemble years

Each ensemble size is bootstrap-resampled from the existing baseline Kirsch
library (``outputs/03_kirsch_library/build_library/library__N=10000__T=20__s=42``,
10,000 traces x 20 water years, SSI-3, seed 42). For an ensemble of ``m``
traces the panel value is the maximum per-trace ``max_duration`` /
``max_magnitude`` across the ``m`` drawn traces; ``m`` traces = ``20 m``
scenario years. A horizontal line marks the historical 1964 drought of
record (Jun 1964 - Apr 1967), the worst event in the 1950-2023 USGS
record at the Cannonsville inflow site.

Duration is used rather than |peak SSI| severity: per-trace
``worst_severity`` saturates at the gamma->normal transform tail ceiling
(~6.36, hit by >50 % of traces) and yields a flat, uninformative curve,
whereas ``max_duration`` grows with ensemble size.

The point: even ~100k Kirsch years do not reach the 1964 duration or
magnitude, so untargeted ensemble generation cannot economically
surface historically extreme conditions.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.data import load_usgs_daily, daily_to_monthly  # noqa: E402
from src.metrics.objectives import make_ssi_calculator  # noqa: E402
from src.plotting.style import apply_style, COLORS  # noqa: E402
from synhydro.droughts.ssi import get_drought_metrics  # noqa: E402

LIBRARY_NPZ = (
    PROJECT_ROOT
    / "outputs/03_kirsch_library/build_library/library__N=10000__T=20__s=42"
    / "characteristics.npz"
)
TRACE_YEARS = 20  # each library trace is 20 water years
SSI_TIMESCALE = 3


def historical_1964_reference() -> tuple[float, float]:
    """Return (duration_months, |magnitude|) of the historical drought of record.

    Fits the SSI-3 calculator on the full water-year-aligned USGS record
    (same gamma/monthly configuration the library used) and selects the
    single worst critical drought event, which is the 1964-1967 drought
    of record.
    """
    cache = PROJECT_ROOT / "outputs" / "data_cache"
    cache.mkdir(parents=True, exist_ok=True)
    monthly = daily_to_monthly(load_usgs_daily(cache_dir=cache))
    first_oct = monthly.index[monthly.index.month == 10][0]
    last_sep = monthly.index[monthly.index.month == 9][-1]
    monthly = monthly[first_oct:last_sep]

    calc = make_ssi_calculator(timescale=SSI_TIMESCALE)
    ssi = calc.fit_transform(monthly)
    dm = get_drought_metrics(ssi)

    worst = dm.loc[dm["severity"].astype(float).idxmin()]
    dur = float(worst["duration"])
    mag = abs(float(worst["magnitude"]))
    print(
        f"[motivation] 1964 reference: {worst['start'].date()} -> "
        f"{worst['end'].date()}  duration={dur:.0f}mo  |magnitude|={mag:.2f}"
    )
    return dur, mag


def load_library() -> tuple[np.ndarray, np.ndarray]:
    """Per-trace (max_duration, max_magnitude) from the baseline library."""
    d = np.load(LIBRARY_NPZ, allow_pickle=True)
    keys = [str(k) for k in d["all_keys"]]
    vals = d["all_values"]
    dur = vals[:, keys.index("max_duration")].astype(float)
    mag = vals[:, keys.index("max_magnitude")].astype(float)
    return dur, mag


def bootstrap_curve(
    per_trace: np.ndarray,
    trace_counts: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Bootstrap ensemble maxima.

    Returns an ``(len(trace_counts), n_bootstrap)`` array: for each
    ensemble size, ``n_bootstrap`` independent resampled maxima of the
    per-trace metric (sampling traces with replacement from the pool).
    """
    n_pool = per_trace.size
    out = np.empty((trace_counts.size, n_bootstrap), dtype=float)
    for i, m in enumerate(trace_counts):
        for b in range(n_bootstrap):
            idx = rng.integers(0, n_pool, size=int(m))
            out[i, b] = per_trace[idx].max()
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--max-years", type=int, default=100_000,
                   help="Largest ensemble size on the x-axis (<= 100k).")
    p.add_argument("--n-bootstrap", type=int, default=10,
                   help="Bootstrap resamples per ensemble size.")
    p.add_argument("--n-grid", type=int, default=30,
                   help="Number of log-spaced ensemble sizes.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=Path,
                   default=PROJECT_ROOT / "figures" / "main"
                   / "fig_motivation_ensemble_years")
    args = p.parse_args()

    dur_pool, mag_pool = load_library()
    n_pool = dur_pool.size
    dur_1964, mag_1964 = historical_1964_reference()

    # Ensemble sizes: log-spaced years -> integer trace counts (>=1),
    # capped so ensemble years <= max-years and <= available pool.
    max_traces = min(args.max_years // TRACE_YEARS, n_pool)
    grid_years = np.logspace(
        np.log10(TRACE_YEARS), np.log10(args.max_years), args.n_grid
    )
    trace_counts = np.clip(
        np.round(grid_years / TRACE_YEARS).astype(int), 1, max_traces
    )
    trace_counts = np.unique(trace_counts)
    ensemble_years = trace_counts * TRACE_YEARS

    rng = np.random.default_rng(args.seed)
    dur_bs = bootstrap_curve(dur_pool, trace_counts, args.n_bootstrap, rng)
    mag_bs = bootstrap_curve(mag_pool, trace_counts, args.n_bootstrap, rng)

    apply_style()
    import matplotlib.pyplot as plt

    fig, (ax_d, ax_m) = plt.subplots(1, 2, figsize=(7.0, 3.2))

    for ax, bs, ref, ylabel in (
        (ax_d, dur_bs, dur_1964, "Maximum drought duration  (months)"),
        (ax_m, mag_bs, mag_1964,
         "Maximum drought magnitude  (|SSI|·month)"),
    ):
        # Bootstrap samples (faint) + median curve.
        x_rep = np.repeat(ensemble_years, bs.shape[1])
        ax.scatter(
            x_rep, bs.ravel(), s=8, color=COLORS["empirical"],
            alpha=0.25, edgecolors="none", zorder=2,
        )
        ax.plot(
            ensemble_years, np.median(bs, axis=1),
            color=COLORS["empirical"], lw=1.6, zorder=3,
            label="Kirsch ensemble",
        )
        ax.axhline(
            ref, color=COLORS["historical"], ls="--", lw=1.4,
            zorder=4, label="1964 drought",
        )
        ax.set_xscale("log")
        ax.set_xlabel("Number of Scenario Years")
        ax.set_ylabel(ylabel)
        ax.set_xlim(ensemble_years.min(), ensemble_years.max())

    ax_d.legend(frameon=False, loc="upper left")

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out.with_suffix(".pdf"))
    fig.savefig(args.out.with_suffix(".png"))
    print(f"[motivation] wrote {args.out.with_suffix('.pdf')}")
    print(f"[motivation] wrote {args.out.with_suffix('.png')}")


if __name__ == "__main__":
    main()
