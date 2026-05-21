"""All-events historical envelope + random-DV coverage probe (Phase B+C).

Combines the two diagnostic phases for the T=10y first-SSI3-drought pivot:

    Phase B: per-event characteristics of every critical SSI-3 drought in
             the 73-year historical record. Output is the historical
             envelope (5-axis distribution + correlation matrix) the
             MOEA targets — replaces the per-block historical comparator.

    Phase C: random-DV coverage probe at T=10y. Generates N synthetic
             sequences with uniform-random DVs in [0, 1], extracts their
             first-event metrics, and reports the fraction with no
             critical SSI-3 event (the at-least-one-event constraint
             binding rate at random DVs). Also computes per-calendar-
             month Anderson-Darling k-sample plausibility statistics
             (synthetic vs historical monthly flows) so the AD signal
             is visible alongside the constraint behavior.

Outputs under the canonical stage paths
``outputs/02_calibration/diagnose_first_event_envelope/<slug>/`` and
``figures/02_calibration/diagnose_first_event_envelope/<slug>/``, where
slug = ``<dv_mode>_T<n_years>``:

    outputs/02_calibration/diagnose_first_event_envelope/<slug>/summary.json
    figures/02_calibration/diagnose_first_event_envelope/<slug>/envelope.png
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.hydrology.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.hydrology.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.metrics.objectives import (  # noqa: E402
    compute_ssi,
    compute_ssi_drought_characteristics,
    flows_to_series,
    get_drought_metrics,
    make_ssi_calculator,
)

STAGE = "02_calibration"
DRIVER = "diagnose_first_event_envelope"


def historical_event_records(monthly_1d: np.ndarray, timescale: int = 3) -> pd.DataFrame:
    """All critical SSI-N events in the historical record, with 5 first-event axes."""
    ssi_series, _calc = compute_ssi(monthly_1d, timescale=timescale)
    dm = get_drought_metrics(ssi_series, end_drought_threshold_months=3)
    if len(dm) == 0:
        return pd.DataFrame(
            columns=["duration", "severity", "magnitude",
                     "start_month", "peak_month"]
        )
    out = pd.DataFrame({
        "duration": dm["duration"].astype(float).values,
        "severity": dm["severity"].abs().astype(float).values,    # peak depth
        "magnitude": dm["magnitude"].abs().astype(float).values,  # cum. deficit
        "start_month": pd.to_datetime(dm["start"]).dt.month.values.astype(float),
        "peak_month": pd.to_datetime(
            dm["max_severity_date"]
        ).dt.month.values.astype(float),
    })
    return out


def random_dv_first_event_records(
    monthly_2d: np.ndarray,
    monthly_1d: np.ndarray,
    n_samples: int,
    n_years_out: int,
    dv_mode: str,
    seed: int,
    timescale: int = 3,
) -> tuple[pd.DataFrame, list[np.ndarray]]:
    """Sample uniform-random DVs, generate synthetic T-year sequences,
    and return per-sample first-event metrics + the raw monthly-flow
    matrix (for the AD plausibility test)."""
    rng = np.random.default_rng(seed)
    kirsch = build_kirsch_generator(monthly_2d)
    gen = KirschBorgWrapper(kirsch, mode=dv_mode, n_years_out=n_years_out)
    n_dvs = gen.n_dvs

    # Pre-fit SSI on historical so per-sample transform is calibrated to history.
    ssi_calc = make_ssi_calculator(timescale=timescale)
    hist_series = flows_to_series(monthly_1d, start_date="1950-10-01")
    ssi_calc.fit(hist_series)

    rows = []
    synthetic_monthly_2d = []
    for i in range(n_samples):
        dvs = rng.random(n_dvs)
        synthetic_2d = gen.generate(dvs)
        if synthetic_2d.ndim == 1:
            synthetic_2d = synthetic_2d.reshape(n_years_out, 12)
        synthetic_1d = synthetic_2d.flatten()
        synthetic_monthly_2d.append(synthetic_2d)

        syn_series = flows_to_series(synthetic_1d, start_date="2100-01-01")
        ssi_syn = ssi_calc.transform(syn_series)
        chars = compute_ssi_drought_characteristics(ssi_syn, monthly_flows=synthetic_1d)

        rows.append({
            "first_event_present": int(chars["first_event_present"]),
            "duration": float(chars["first_event_duration"]),
            "severity": float(chars["first_event_severity"]),
            "magnitude": float(chars["first_event_magnitude"]),
            "start_month": float(chars["first_event_start_month"]),
            "peak_month": float(chars["first_event_peak_month"]),
        })
    return pd.DataFrame(rows), synthetic_monthly_2d


def ad_plausibility_per_month(
    historical_2d: np.ndarray,
    synthetic_2d_list: list[np.ndarray],
) -> dict[str, float]:
    """Per-calendar-month Anderson-Darling k-sample statistic comparing
    synthetic vs historical monthly flow distributions. The AD A² is
    a tail-weighted goodness-of-fit measure: small A² (≲ 1) is consistent
    with the null of identical distributions; A² ≳ 4 is strong rejection.
    """
    from scipy.stats import anderson_ksamp

    # historical_2d: (n_hist_years, 12); synthetic stack: (n_samples * n_years, 12)
    syn_stack = np.vstack(synthetic_2d_list)
    out: dict[str, float] = {}
    for m in range(12):
        hist_col = historical_2d[:, m]
        syn_col = syn_stack[:, m]
        try:
            res = anderson_ksamp([hist_col, syn_col])
            out[f"month_{m + 1:02d}"] = float(res.statistic)
        except Exception:
            out[f"month_{m + 1:02d}"] = float("nan")
    out["mean_A2"] = float(np.nanmean(list(out.values())))
    out["max_A2"] = float(np.nanmax([v for k, v in out.items() if k != "mean_A2"]))
    return out


def plot_envelope(
    hist_df: pd.DataFrame,
    syn_df: pd.DataFrame,
    out_path: Path,
) -> None:
    axes_names = ["duration", "severity", "magnitude", "start_month", "peak_month"]
    n = len(axes_names)
    fig, axes = plt.subplots(n, n, figsize=(2.4 * n, 2.4 * n))
    syn_present = syn_df[syn_df["first_event_present"] == 1]
    for i, yname in enumerate(axes_names):
        for j, xname in enumerate(axes_names):
            ax = axes[i, j]
            if i == j:
                ax.hist(hist_df[xname], bins=20, alpha=0.6, color="C0",
                        label=f"historical (n={len(hist_df)})", density=True)
                if len(syn_present) > 0:
                    ax.hist(syn_present[xname], bins=20, alpha=0.4, color="C1",
                            label=f"random-DV (n={len(syn_present)})", density=True)
                if i == 0:
                    ax.legend(fontsize=6)
            else:
                ax.scatter(hist_df[xname], hist_df[yname], s=18, alpha=0.7,
                           color="C0", label="historical")
                if len(syn_present) > 0:
                    ax.scatter(syn_present[xname], syn_present[yname], s=4,
                               alpha=0.25, color="C1", label="random-DV")
            if i == n - 1:
                ax.set_xlabel(xname, fontsize=8)
            if j == 0:
                ax.set_ylabel(yname, fontsize=8)
            ax.tick_params(labelsize=6)
    fig.suptitle("First-event SSI-3 envelope: historical (all events) vs random-DV synthetic")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-years", type=int, default=10,
                   help="Synthetic candidate trace length (years).")
    p.add_argument("--n-samples", type=int, default=2000,
                   help="Number of random-DV synthetic sequences.")
    p.add_argument("--dv-mode", choices=["index", "residual"], default="residual",
                   help="DV mapping mode for the Kirsch wrapper.")
    p.add_argument("--ssi-timescale", type=int, default=3,
                   help="SSI accumulation period (default 3 = SSI-3).")
    p.add_argument("--seed", type=int, default=20260506)
    p.add_argument("--cache-dir", type=Path,
                   default=PROJECT_ROOT / "data" / "usgs_cache")
    args = p.parse_args()

    slug = f"{args.dv_mode}_T{args.n_years}"
    out_dir = stage_output_dir(STAGE, DRIVER, slug)
    fig_dir = stage_figure_dir(STAGE, DRIVER, slug)

    print(f"[diagnose_first_event] T={args.n_years}y, mode={args.dv_mode}, "
          f"N={args.n_samples}, seed={args.seed}")

    monthly_2d, monthly_1d = prepare_data(args.cache_dir)

    print(f"[diagnose_first_event] Historical record: {monthly_2d.shape[0]} years")
    hist_df = historical_event_records(monthly_1d, timescale=args.ssi_timescale)
    print(f"[diagnose_first_event] Historical critical SSI-{args.ssi_timescale} events: "
          f"n={len(hist_df)}")

    syn_df, synthetic_2d_list = random_dv_first_event_records(
        monthly_2d, monthly_1d,
        n_samples=args.n_samples, n_years_out=args.n_years,
        dv_mode=args.dv_mode, seed=args.seed,
        timescale=args.ssi_timescale,
    )
    n_present = int(syn_df["first_event_present"].sum())
    rate_empty = 1.0 - n_present / max(len(syn_df), 1)
    print(f"[diagnose_first_event] Random-DV present rate: {n_present}/{len(syn_df)} "
          f"({100 * (1 - rate_empty):.1f}% feasible; "
          f"{100 * rate_empty:.1f}% empty-event)")

    ad_stats = ad_plausibility_per_month(monthly_2d, synthetic_2d_list)
    print(f"[diagnose_first_event] AD k-sample (synthetic vs historical, by month): "
          f"mean A²={ad_stats['mean_A2']:.3f}, max A²={ad_stats['max_A2']:.3f}")

    syn_present = syn_df[syn_df["first_event_present"] == 1]
    corr_hist = hist_df.corr().to_dict() if len(hist_df) else {}
    corr_syn = syn_present.corr().to_dict() if len(syn_present) else {}

    summary = {
        "config": {
            "n_years_out": args.n_years,
            "dv_mode": args.dv_mode,
            "ssi_timescale": args.ssi_timescale,
            "n_samples": args.n_samples,
            "seed": args.seed,
        },
        "historical": {
            "n_events": int(len(hist_df)),
            "stats": {
                col: {
                    "min": float(hist_df[col].min()),
                    "q25": float(hist_df[col].quantile(0.25)),
                    "median": float(hist_df[col].median()),
                    "q75": float(hist_df[col].quantile(0.75)),
                    "max": float(hist_df[col].max()),
                } for col in ["duration", "severity", "magnitude",
                              "start_month", "peak_month"]
            } if len(hist_df) else {},
            "correlation": corr_hist,
        },
        "random_dv": {
            "n_samples": int(len(syn_df)),
            "n_present": n_present,
            "constraint_binding_rate": float(rate_empty),
            "stats_present_only": {
                col: {
                    "min": float(syn_present[col].min()),
                    "q25": float(syn_present[col].quantile(0.25)),
                    "median": float(syn_present[col].median()),
                    "q75": float(syn_present[col].quantile(0.75)),
                    "max": float(syn_present[col].max()),
                } for col in ["duration", "severity", "magnitude",
                              "start_month", "peak_month"]
            } if len(syn_present) else {},
            "correlation": corr_syn,
        },
        "ad_plausibility": ad_stats,
    }

    json_path = out_dir / "summary.json"
    json_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"[diagnose_first_event] wrote {json_path}")

    fig_path = fig_dir / "envelope.png"
    plot_envelope(hist_df, syn_df, fig_path)
    print(f"[diagnose_first_event] wrote {fig_path}")


if __name__ == "__main__":
    main()
