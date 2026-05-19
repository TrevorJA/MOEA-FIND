"""Render Stage-08 SA diagnostic figures from cached compute outputs.

Reads the parquet artifacts written by
``workflows/08_nyc_sensitivity/run_sa.py`` under
``outputs/08_nyc_sensitivity/run_sa/<slug>/results/`` and emits PDFs to
``figures/08_nyc_sensitivity/run_sa/<slug>/``.

This is a plotting-only driver — it never re-runs SA. The slug must
match the compute-driver invocation. Figures emitted:

    tornado_<outcome>_<method>.pdf
    heatmap_indices_<method>.pdf            (factor x outcome)
    convergence_<outcome>.pdf               (method-faceted)
    cross_method_rank_corr_<outcome>.pdf
    cross_outcome_rank_corr_<method>.pdf
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.sensitivity.sensitivity import HEADLINE_INDEX  # noqa: E402

STAGE = "08_nyc_sensitivity"
DRIVER = "run_sa"


def _wide_indices(idx_long: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """Reshape long-form indices parquet into ``{outcome: df_indexed_by_factor}``."""
    out: Dict[str, pd.DataFrame] = {}
    for outcome, sub in idx_long.groupby("outcome"):
        df = sub.drop(columns=["outcome"]).set_index("factor")
        out[str(outcome)] = df
    return out


def _wide_bootstrap(boot_long: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    out: Dict[str, pd.DataFrame] = {}
    for outcome, sub in boot_long.groupby("outcome"):
        df = sub.drop(columns=["outcome"]).set_index("factor")
        out[str(outcome)] = df
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True,
                   help="Variant slug under outputs/08_nyc_sensitivity/run_sa/")
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    results_dir = in_dir / "results"
    if not results_dir.exists():
        sys.exit(f"missing {results_dir} — run the compute driver first")

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)

    # Discover methods from indices_<method>.parquet files present.
    methods: List[str] = sorted(
        p.stem.removeprefix("indices_")
        for p in results_dir.glob("indices_*.parquet")
    )
    if not methods:
        sys.exit(f"no indices_*.parquet files in {results_dir}")
    print(f"[plots/08/run_sa] methods: {methods}")

    # Load per-method long-form tables and reshape to per-outcome dicts.
    method_results: Dict[str, Dict[str, pd.DataFrame]] = {}
    bootstrap_results: Dict[str, Dict[str, pd.DataFrame]] = {}
    convergence_results: Dict[str, Dict[str, pd.DataFrame]] = {}
    for m in methods:
        method_results[m] = _wide_indices(
            pd.read_parquet(results_dir / f"indices_{m}.parquet")
        )
        bootstrap_results[m] = _wide_bootstrap(
            pd.read_parquet(results_dir / f"bootstrap_{m}.parquet")
        )
        conv_long = pd.read_parquet(results_dir / f"convergence_{m}.parquet")
        convergence_results[m] = {
            str(oc): sub.drop(columns=["outcome", "method"], errors="ignore")
            for oc, sub in conv_long.groupby("outcome")
        }

    # Outcome list = union across methods, sorted.
    outcomes = sorted({
        oc for m in methods for oc in method_results[m].keys()
    })
    print(f"[plots/08/run_sa] outcomes: {outcomes}")

    # Optional: read selection_log.json to log the anchor choice.
    sel_path = results_dir / "selection_log.json"
    if sel_path.exists():
        sel_log = json.loads(sel_path.read_text())
        anchors = {oc: sel_log.get(oc, {}).get("_anchor", "none")
                   for oc in outcomes}
        print(f"[plots/08/run_sa] anchors: {anchors}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from src.plotting.sensitivity import (
        plot_convergence,
        plot_index_heatmap,
        plot_rank_consistency,
        plot_rank_correlation,
        plot_tornado,
        plot_tornado_matrix,
    )
    from src.sensitivity.sensitivity import cross_method_rank_corr, cross_outcome_rank_corr

    # Tornado per (outcome, method) — uses the bootstrap CI for whiskers.
    for method_name in methods:
        headline = HEADLINE_INDEX[method_name]
        for outcome_label in outcomes:
            if outcome_label not in method_results[method_name]:
                continue
            indices_df = method_results[method_name][outcome_label].copy()
            boot_df = bootstrap_results[method_name].get(outcome_label)
            if boot_df is not None:
                # Replace per-method CI columns with diagnostic-layer
                # bootstrap CI so all tornadoes use a uniform definition.
                indices_df["ci_lo"] = boot_df["ci_lo"]
                indices_df["ci_hi"] = boot_df["ci_hi"]
            plot_tornado(
                indices_df, headline_col=headline,
                title=f"{method_name}: {outcome_label}",
                output_path=fig_dir / f"tornado_{outcome_label}_{method_name}.pdf",
            )

    # Heatmap factor x outcome per method.
    for method_name in methods:
        headline = HEADLINE_INDEX[method_name]
        plot_index_heatmap(
            method_results[method_name],
            method_name=method_name,
            headline_col=headline,
            output_path=fig_dir / f"heatmap_indices_{method_name}.pdf",
        )

    # Convergence per outcome (method-faceted via stacked subplots).
    for outcome_label in outcomes:
        present_methods = [m for m in methods
                           if outcome_label in convergence_results[m]]
        if not present_methods:
            continue
        fig, axes = plt.subplots(
            len(present_methods), 1,
            figsize=(6.5, 3.0 * len(present_methods)),
            sharex=True,
        )
        if len(present_methods) == 1:
            axes = [axes]
        for ax, method_name in zip(axes, present_methods):
            plot_convergence(
                convergence_results[method_name][outcome_label],
                title=f"{method_name}: {outcome_label}",
                ax=ax,
            )
        fig.tight_layout()
        fig.savefig(fig_dir / f"convergence_{outcome_label}.pdf")
        plt.close(fig)

    # Cross-method rank correlation per outcome.
    for outcome_label in outcomes:
        per_method = {
            m: method_results[m][outcome_label]
            for m in methods if outcome_label in method_results[m]
        }
        if len(per_method) < 2:
            continue
        rho = cross_method_rank_corr(per_method)
        plot_rank_correlation(
            rho, title=f"Cross-method rank correlation: {outcome_label}",
            output_path=fig_dir / f"cross_method_rank_corr_{outcome_label}.pdf",
        )

    # Cross-outcome rank correlation per method.
    for method_name in methods:
        if len(method_results[method_name]) < 2:
            continue
        rho = cross_outcome_rank_corr(
            method_results[method_name], method=method_name,
        )
        plot_rank_correlation(
            rho, title=f"Cross-outcome rank correlation: {method_name}",
            output_path=fig_dir / f"cross_outcome_rank_corr_{method_name}.pdf",
        )

    # ------------------------------------------------------------------
    # Headline literature-style summary figures (per method).
    # ------------------------------------------------------------------

    # Display labels: drop verbose `_log1p` suffix for axis legibility.
    OUTCOME_DISPLAY = {
        "nyc_min_storage_frac":              "NYC min storage frac",
        "nyc_drawdown_days_below_0.25":      "NYC drawdown days <25%",
        "nyc_drawdown_days_below_0.25_log1p":"NYC drawdown days <25% (log1p)",
        "montague_flow_reliability":         "Montague flow reliability",
        "montague_flow_vulnerability":       "Montague flow vulnerability",
        "montague_flow_vulnerability_log1p": "Montague flow vulnerability (log1p)",
    }
    FACTOR_DISPLAY = {
        "mean_severity":            "Mean severity",
        "mean_magnitude":           "Mean magnitude",
        "time_in_drought_fraction": "Time in drought (frac.)",
    }

    for method_name in methods:
        headline = HEADLINE_INDEX[method_name]
        idx_long = pd.read_parquet(results_dir / f"indices_{method_name}.parquet")
        boot_long = pd.read_parquet(results_dir / f"bootstrap_{method_name}.parquet")

        # Tornado matrix: factor stack of horizontal bars across outcomes.
        plot_tornado_matrix(
            idx_long, boot_long, headline_col=headline,
            outcome_labels=OUTCOME_DISPLAY,
            factor_labels=FACTOR_DISPLAY,
            noise_floor=0.05,  # Delta moment-independent noise floor
            title=f"Sensitivity of NYC outcomes to drought-hazard "
                  f"characteristics ({method_name.upper()} index)",
            output_path=fig_dir / f"summary_tornado_matrix_{method_name}.pdf",
        )

        # Rank-consistency dot plot.
        plot_rank_consistency(
            idx_long, headline_col=headline,
            outcome_labels=OUTCOME_DISPLAY,
            factor_labels=FACTOR_DISPLAY,
            title=f"Factor rank by outcome ({method_name.upper()})",
            output_path=fig_dir / f"summary_rank_consistency_{method_name}.pdf",
        )

    print(f"[plots/08/run_sa] figures -> {fig_dir}")


if __name__ == "__main__":
    main()
