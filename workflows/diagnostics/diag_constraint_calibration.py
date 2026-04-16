"""Diagnostic: bootstrap calibration of plausibility-constraint tolerances.

Produces calibrated tolerances for every statistic used in MOEA-FIND plausibility
constraints, at the target synthetic trace length (T=20 years). The output JSON
feeds `src.constraints.ConstraintConfig` and drives all downstream HPC runs.

The calibration compares two bootstrap distributions for each statistic:

    (a) Historical block bootstrap. Draw 2000 contiguous T-year blocks from the
        historical record with wrap-around, compute the statistic on each block.
        This is the natural finite-sample noise floor of a hydrologically honest
        20-year sample.

    (b) Kirsch-generator bootstrap. Draw 2000 standard Kirsch traces of length T
        from the fitted generator (no optimizer), compute the same statistic.
        This is how much an unperturbed Kirsch draw moves away from the historical
        reference purely from generator sampling.

The constraint tolerance for a statistic is
    max( half-width of historical 95% CI, half-width of Kirsch 95% CI )
so that neither an honest historical window nor an honest Kirsch draw would be
flagged infeasible.

Statistics:
    - annual_mean            (fractional dev. from historical annual mean)
    - annual_cv              (fractional dev. from historical annual CV)
    - lag1_ac_monthly        (absolute dev. from historical lag-1 AC)
    - non_drought_mean       (fractional dev. from historical non-drought mean)
    - seasonal_cycle_max_dev (max across months of fractional dev. per month)

Note: the Hurst-exponent constraint mentioned in constraints_spec.md §4.2 is
intentionally excluded here per the scoping decision for this plan.

Usage:
    python scripts/diag_constraint_calibration.py --n-boot 2000 --T 20

Outputs under outputs/diag_constraint_calibration/:
    - calibrated_tolerances.json
    - constraint_calibration.pdf  (per-statistic histogram panels)
    - config.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment_utils import prepare_data  # noqa: E402
from src.objectives import (  # noqa: E402
    flows_to_series,
    make_ssi_calculator,
)

OUTPUT_SLUG = "diag_constraint_calibration"


# ---------------------------------------------------------------------------
# Statistic functions (all operate on a 1D monthly flow vector of length T*12)
# ---------------------------------------------------------------------------

def _annual_totals(monthly_1d: np.ndarray) -> np.ndarray:
    """Sum monthly flows into annual totals. Assumes len divisible by 12."""
    n_years = len(monthly_1d) // 12
    return monthly_1d[: n_years * 12].reshape(n_years, 12).sum(axis=1)


def stat_annual_mean(monthly_1d: np.ndarray) -> float:
    return float(np.mean(_annual_totals(monthly_1d)))


def stat_annual_cv(monthly_1d: np.ndarray) -> float:
    totals = _annual_totals(monthly_1d)
    mean = totals.mean()
    if mean <= 0:
        return 0.0
    return float(totals.std(ddof=1) / mean)


def stat_lag1_ac_monthly(monthly_1d: np.ndarray) -> float:
    x = np.asarray(monthly_1d, dtype=float)
    if len(x) < 3:
        return 0.0
    return float(np.corrcoef(x[:-1], x[1:])[0, 1])


def _seasonal_cycle(monthly_1d: np.ndarray) -> np.ndarray:
    """Per-calendar-month mean flow over the trace."""
    n_years = len(monthly_1d) // 12
    return monthly_1d[: n_years * 12].reshape(n_years, 12).mean(axis=0)


def stat_seasonal_cycle_max_fractional_dev(
    monthly_1d: np.ndarray,
    historical_monthly_means: np.ndarray,
) -> float:
    """Max over months of |sample_month_mean / hist_month_mean - 1|."""
    cycle = _seasonal_cycle(monthly_1d)
    hist = np.asarray(historical_monthly_means, dtype=float)
    safe_hist = np.where(hist > 0, hist, 1.0)
    rel = np.abs(cycle / safe_hist - 1.0)
    return float(np.max(rel))


def stat_non_drought_mean(
    monthly_1d: np.ndarray,
    ssi_calc,
    start_date: str = "2100-01-01",
) -> float:
    """Mean flow across months whose SSI is strictly positive.

    Uses a pre-fitted SSI calculator (fitted on the full historical record) so
    that the drought/non-drought split is defined against the same reference
    distribution for every sample.
    """
    series = flows_to_series(monthly_1d, start_date=start_date)
    ssi = ssi_calc.transform(series)
    ssi_arr = np.asarray(ssi.values, dtype=float)
    # Align flows to SSI output (SSI drops leading months for accumulation > 1)
    flows_arr = np.asarray(series.loc[ssi.index].values, dtype=float)
    mask = np.isfinite(ssi_arr) & (ssi_arr > 0)
    if not mask.any():
        return float(np.nan)
    return float(flows_arr[mask].mean())


# ---------------------------------------------------------------------------
# Bootstrap drivers
# ---------------------------------------------------------------------------

def historical_block_bootstrap(
    monthly_1d_hist: np.ndarray,
    T_years: int,
    n_boot: int,
    seed: int,
) -> np.ndarray:
    """Return a (n_boot, T_years*12) array of contiguous wrap-around blocks.

    monthly_1d_hist must have length divisible by 12 (water-year aligned).
    """
    n_years_hist = len(monthly_1d_hist) // 12
    if n_years_hist < T_years:
        raise ValueError(
            f"Historical record has {n_years_hist} years, need at least {T_years}."
        )

    yearly_grid = monthly_1d_hist[: n_years_hist * 12].reshape(n_years_hist, 12)
    rng = np.random.default_rng(seed)
    start_idx = rng.integers(0, n_years_hist, size=n_boot)

    samples = np.empty((n_boot, T_years * 12), dtype=float)
    for i, s in enumerate(start_idx):
        idx = (np.arange(s, s + T_years) % n_years_hist)
        samples[i] = yearly_grid[idx].flatten()
    return samples


def kirsch_bootstrap(
    kirsch_gen,
    T_years: int,
    n_boot: int,
    seed: int,
) -> np.ndarray:
    """Return a (n_boot, T_years*12) array of standard Kirsch traces."""
    ensemble = kirsch_gen.generate(
        n_realizations=n_boot,
        n_years=T_years,
        seed=seed,
    )
    n_months = T_years * 12
    traces = np.empty((n_boot, n_months), dtype=float)
    for i, rid in enumerate(sorted(ensemble.data_by_realization.keys())):
        arr = ensemble.data_by_realization[rid].values.flatten()
        traces[i] = arr[:n_months]
    return traces


# ---------------------------------------------------------------------------
# Tolerance derivation
# ---------------------------------------------------------------------------

def _ci_half_width(values: np.ndarray, reference: float, mode: str) -> Tuple[float, float, float]:
    """Return (lo, hi, half_width) of the 95% empirical interval.

    mode = "fractional" reports (sample / reference - 1).
    mode = "absolute"   reports (sample - reference).
    Half-width is max(|lo|, |hi|), i.e. the symmetric tolerance that covers 95%.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if mode == "fractional":
        if reference == 0:
            dev = np.zeros_like(values)
        else:
            dev = values / reference - 1.0
    elif mode == "absolute":
        dev = values - reference
    else:
        raise ValueError(f"Unknown mode {mode!r}")

    lo, hi = np.quantile(dev, [0.025, 0.975])
    return float(lo), float(hi), float(max(abs(lo), abs(hi)))


def summarize_statistic(
    name: str,
    hist_samples: np.ndarray,
    kirsch_samples: np.ndarray,
    reference: float,
    mode: str,
) -> Dict:
    """Compute bootstrap CIs for both distributions and combine into tolerance.

    hist_samples / kirsch_samples are 1D arrays of statistic values per bootstrap
    draw. The final calibrated tolerance is
        max( historical half-width, kirsch half-width )
    so that neither distribution is systematically flagged infeasible.
    """
    h_lo, h_hi, h_hw = _ci_half_width(hist_samples, reference, mode)
    k_lo, k_hi, k_hw = _ci_half_width(kirsch_samples, reference, mode)
    tolerance = float(max(h_hw, k_hw))
    return {
        "name": name,
        "mode": mode,
        "reference": float(reference),
        "historical_ci_95": [h_lo, h_hi],
        "historical_half_width": h_hw,
        "kirsch_ci_95": [k_lo, k_hi],
        "kirsch_half_width": k_hw,
        "tolerance": tolerance,
        "n_hist": int(np.isfinite(hist_samples).sum()),
        "n_kirsch": int(np.isfinite(kirsch_samples).sum()),
    }


# ---------------------------------------------------------------------------
# Diagnostic plotting
# ---------------------------------------------------------------------------

def make_diagnostic_pdf(
    summaries: Dict[str, Dict],
    hist_stat_samples: Dict[str, np.ndarray],
    kirsch_stat_samples: Dict[str, np.ndarray],
    output_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = list(summaries.keys())
    n = len(names)
    fig, axes = plt.subplots(n, 1, figsize=(8, 2.4 * n), constrained_layout=True)
    if n == 1:
        axes = [axes]

    for ax, name in zip(axes, names):
        s = summaries[name]
        ref = s["reference"]
        mode = s["mode"]
        h = hist_stat_samples[name]
        k = kirsch_stat_samples[name]
        if mode == "fractional":
            h_dev = h / ref - 1.0 if ref else h * 0.0
            k_dev = k / ref - 1.0 if ref else k * 0.0
            xlabel = f"{name} fractional deviation"
        else:
            h_dev = h - ref
            k_dev = k - ref
            xlabel = f"{name} absolute deviation"

        h_dev = h_dev[np.isfinite(h_dev)]
        k_dev = k_dev[np.isfinite(k_dev)]

        bins = 40
        ax.hist(h_dev, bins=bins, alpha=0.5, color="#2b6cb0",
                label=f"Historical (n={len(h_dev)})", density=True)
        ax.hist(k_dev, bins=bins, alpha=0.5, color="#c05621",
                label=f"Kirsch (n={len(k_dev)})", density=True)
        tol = s["tolerance"]
        ax.axvline(-tol, color="black", linestyle="--", lw=1)
        ax.axvline(+tol, color="black", linestyle="--", lw=1,
                   label=f"tol = {tol:.3f}")
        ax.set_title(name)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("density")
        ax.legend(fontsize=8, loc="upper right")

    fig.suptitle("Plausibility-constraint bootstrap calibration (T=20y)",
                 fontsize=12)
    fig.savefig(output_path, format="pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _build_kirsch_generator(monthly_2d: np.ndarray):
    import pandas as pd
    from synhydro.methods.generation.nonparametric.kirsch import KirschGenerator
    gen = KirschGenerator(generate_using_log_flow=True)
    dates = pd.date_range(start="1950-10-01", periods=monthly_2d.size, freq="MS")
    gen.fit(pd.DataFrame({"flow_cfs": monthly_2d.flatten()}, index=dates))
    return gen


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T", type=int, default=20,
                   help="Synthetic trace length in years (target for calibration).")
    p.add_argument("--n-boot", type=int, default=2000,
                   help="Number of bootstrap draws per distribution.")
    p.add_argument("--ssi-timescale", type=int, default=3, choices=[1, 3, 6, 12])
    p.add_argument("--seed", type=int, default=20260415)
    p.add_argument("--site-label", default="cannonsville",
                   help="Label used as the top-level key in the output JSON.")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    (out / "config.json").write_text(json.dumps({
        "script": "diag_constraint_calibration.py",
        "T_years": args.T,
        "n_boot": args.n_boot,
        "ssi_timescale": args.ssi_timescale,
        "seed": args.seed,
        "site_label": args.site_label,
    }, indent=2))

    print(f"[calib] loading historical data ...")
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)
    n_years_hist = monthly_2d.shape[0]
    print(f"[calib] historical: {n_years_hist} water years, "
          f"mean={monthly_1d.mean():.1f} cfs")

    # Precompute reference statistics on the full historical record
    hist_annual_totals = _annual_totals(monthly_1d)
    ref_annual_mean = float(hist_annual_totals.mean())
    ref_annual_cv = float(hist_annual_totals.std(ddof=1) / ref_annual_mean)
    ref_lag1_ac = stat_lag1_ac_monthly(monthly_1d)
    ref_monthly_means = monthly_2d.mean(axis=0)

    # Fit SSI once on full history; reused for non-drought classification of every sample.
    ssi_calc = make_ssi_calculator(timescale=args.ssi_timescale)
    hist_series = flows_to_series(monthly_1d, start_date="1950-10-01")
    ssi_calc.fit(hist_series)
    ref_non_drought_mean = stat_non_drought_mean(
        monthly_1d, ssi_calc, start_date="1950-10-01"
    )
    # Seasonal cycle reference deviation is always 0 by construction; we calibrate
    # on the *absolute* max-fractional-deviation statistic itself (not a difference).
    ref_seasonal_max_dev = 0.0

    print(f"[calib] historical ref annual_mean = {ref_annual_mean:.1f}")
    print(f"[calib] historical ref annual_cv   = {ref_annual_cv:.3f}")
    print(f"[calib] historical ref lag1_ac     = {ref_lag1_ac:.3f}")
    print(f"[calib] historical ref non_drought_mean = {ref_non_drought_mean:.1f}")

    # ----- Historical bootstrap -----
    print(f"[calib] drawing {args.n_boot} historical T={args.T}y blocks ...")
    h_samples = historical_block_bootstrap(
        monthly_1d, T_years=args.T, n_boot=args.n_boot, seed=args.seed,
    )

    # ----- Kirsch bootstrap -----
    print(f"[calib] building Kirsch generator and drawing {args.n_boot} traces ...")
    kirsch_gen = _build_kirsch_generator(monthly_2d)
    k_samples = kirsch_bootstrap(
        kirsch_gen, T_years=args.T, n_boot=args.n_boot, seed=args.seed + 1,
    )

    # ----- Compute statistic distributions -----
    def _batch(statfn, samples, desc):
        out_arr = np.empty(samples.shape[0], dtype=float)
        for i, row in enumerate(samples):
            out_arr[i] = statfn(row)
        print(f"[calib]   {desc}: mean={np.nanmean(out_arr):.4f} "
              f"std={np.nanstd(out_arr):.4f}")
        return out_arr

    print("[calib] computing statistics on historical bootstrap ...")
    h_annual_mean = _batch(stat_annual_mean, h_samples, "h annual_mean")
    h_annual_cv = _batch(stat_annual_cv, h_samples, "h annual_cv")
    h_lag1 = _batch(stat_lag1_ac_monthly, h_samples, "h lag1_ac")
    h_seasonal = _batch(
        lambda x: stat_seasonal_cycle_max_fractional_dev(x, ref_monthly_means),
        h_samples, "h seasonal_max_dev",
    )
    h_non_drought = _batch(
        lambda x: stat_non_drought_mean(x, ssi_calc),
        h_samples, "h non_drought_mean",
    )

    print("[calib] computing statistics on Kirsch bootstrap ...")
    k_annual_mean = _batch(stat_annual_mean, k_samples, "k annual_mean")
    k_annual_cv = _batch(stat_annual_cv, k_samples, "k annual_cv")
    k_lag1 = _batch(stat_lag1_ac_monthly, k_samples, "k lag1_ac")
    k_seasonal = _batch(
        lambda x: stat_seasonal_cycle_max_fractional_dev(x, ref_monthly_means),
        k_samples, "k seasonal_max_dev",
    )
    k_non_drought = _batch(
        lambda x: stat_non_drought_mean(x, ssi_calc),
        k_samples, "k non_drought_mean",
    )

    # ----- Summaries and tolerance derivation -----
    summaries = {
        "annual_mean": summarize_statistic(
            "annual_mean", h_annual_mean, k_annual_mean,
            ref_annual_mean, "fractional",
        ),
        "annual_cv": summarize_statistic(
            "annual_cv", h_annual_cv, k_annual_cv,
            ref_annual_cv, "fractional",
        ),
        "lag1_ac_monthly": summarize_statistic(
            "lag1_ac_monthly", h_lag1, k_lag1,
            ref_lag1_ac, "absolute",
        ),
        "non_drought_mean": summarize_statistic(
            "non_drought_mean", h_non_drought, k_non_drought,
            ref_non_drought_mean, "fractional",
        ),
        # For seasonal cycle, both distributions measure the same "max fractional
        # deviation" random variable; reference is 0 and mode is absolute, so
        # tolerance becomes the 97.5 percentile of the worse of the two.
        "seasonal_cycle_max_dev": summarize_statistic(
            "seasonal_cycle_max_dev", h_seasonal, k_seasonal,
            ref_seasonal_max_dev, "absolute",
        ),
    }

    # ----- Output JSON -----
    payload = {
        f"{args.site_label}_T{args.T}": {
            "T_years": args.T,
            "n_boot": args.n_boot,
            "historical": {
                "annual_mean": ref_annual_mean,
                "annual_cv": ref_annual_cv,
                "lag1_ac_monthly": ref_lag1_ac,
                "non_drought_mean": ref_non_drought_mean,
                "monthly_means": ref_monthly_means.tolist(),
            },
            "tolerances": {
                name: s["tolerance"] for name, s in summaries.items()
            },
            "details": summaries,
        }
    }
    json_path = out / "calibrated_tolerances.json"
    json_path.write_text(json.dumps(payload, indent=2))
    print(f"[calib] wrote {json_path}")

    hist_map = {
        "annual_mean": h_annual_mean,
        "annual_cv": h_annual_cv,
        "lag1_ac_monthly": h_lag1,
        "non_drought_mean": h_non_drought,
        "seasonal_cycle_max_dev": h_seasonal,
    }
    kirsch_map = {
        "annual_mean": k_annual_mean,
        "annual_cv": k_annual_cv,
        "lag1_ac_monthly": k_lag1,
        "non_drought_mean": k_non_drought,
        "seasonal_cycle_max_dev": k_seasonal,
    }
    pdf_path = out / "constraint_calibration.pdf"
    make_diagnostic_pdf(summaries, hist_map, kirsch_map, pdf_path)
    print(f"[calib] wrote {pdf_path}")

    print("[calib] Calibrated tolerances:")
    for name, s in summaries.items():
        print(f"  {name:<24} tol = {s['tolerance']:.4f}  "
              f"(hist={s['historical_half_width']:.4f}, "
              f"kirsch={s['kirsch_half_width']:.4f}, mode={s['mode']})")


if __name__ == "__main__":
    main()
