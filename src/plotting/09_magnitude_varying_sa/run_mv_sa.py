"""Render Stage-09 MV-SA figures from cached compute outputs.

Reads ``mv_sa_<method>.parquet`` from
``outputs/09_magnitude_varying_sa/run_mv_sa/<slug>/results/`` and
emits PDFs to ``figures/09_magnitude_varying_sa/run_mv_sa/<slug>/``:

    stacked_area_<method>.pdf       per-method headline (factor share vs tau)
    lines_with_ci_<method>.pdf      per-factor lines with bootstrap CIs
    method_panel.pdf                multi-method comparison (one row per method)

This is a plotting-only driver; it never re-runs SA.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "09_magnitude_varying_sa"
DRIVER = "run_mv_sa"


# Display-name overrides matched to Stage-08 plotting driver so
# headline figures across stages use the same vocabulary.
FACTOR_DISPLAY = {
    "mean_severity":            "Mean severity",
    "mean_magnitude":           "Mean magnitude",
    "time_in_drought_fraction": "Time in drought (frac.)",
    "control_uniform":          "Control (uniform)",
}

AXIS_DISPLAY = {
    "nyc_min_storage_frac":         "NYC min storage frac",
    "nyc_drawdown_days_below_0.25": "NYC drawdown days <25%",
    "montague_flow_reliability":    "Montague flow reliability",
    "montague_flow_vulnerability":  "Montague flow vulnerability",
}

METHOD_DISPLAY = {
    "delta":    "Delta moment-independent",
    "pawn":     "PAWN density-based",
    "rbd_fast": "RBD-FAST",
}


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True,
                   help="Variant slug under "
                        "outputs/09_magnitude_varying_sa/run_mv_sa/")
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    results_dir = in_dir / "results"
    if not results_dir.exists():
        sys.exit(f"missing {results_dir} -- run the compute driver first")

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)

    # Pull the magnitude axis label from the saved config so the
    # x-axis label matches what the user actually ran.
    cfg_path = in_dir / "config.json"
    cfg = json.loads(cfg_path.read_text()) if cfg_path.exists() else {}
    axis_col = cfg.get("magnitude_axis", "")
    axis_label = AXIS_DISPLAY.get(axis_col, axis_col or "Magnitude")
    magnitude_label = f"Percentile of {axis_label}"

    method_paths = sorted(results_dir.glob("mv_sa_*.parquet"))
    method_paths = [p for p in method_paths if p.stem != "mv_sa_combined"]
    if not method_paths:
        sys.exit(f"no mv_sa_*.parquet files in {results_dir}")

    methods = [p.stem.removeprefix("mv_sa_") for p in method_paths]
    print(f"[plots/09/run_mv_sa] methods: {methods}")

    import matplotlib
    matplotlib.use("Agg")

    from src.plotting.magnitude_varying_sa import (
        plot_mv_sa_lines_with_ci,
        plot_mv_sa_method_panel,
        plot_mv_sa_stacked_area,
    )

    dfs = {m: pd.read_parquet(results_dir / f"mv_sa_{m}.parquet")
           for m in methods}

    for m, df in dfs.items():
        plot_mv_sa_stacked_area(
            df,
            factor_labels=FACTOR_DISPLAY,
            title=f"{METHOD_DISPLAY.get(m, m)}: factor share vs "
                  f"{axis_label} percentile",
            magnitude_label=magnitude_label,
            output_path=fig_dir / f"stacked_area_{m}.pdf",
        )
        plot_mv_sa_lines_with_ci(
            df,
            factor_labels=FACTOR_DISPLAY,
            title=f"{METHOD_DISPLAY.get(m, m)}: per-factor "
                  f"sensitivity vs {axis_label} percentile",
            magnitude_label=magnitude_label,
            output_path=fig_dir / f"lines_with_ci_{m}.pdf",
        )

    if len(dfs) >= 2:
        plot_mv_sa_method_panel(
            dfs,
            factor_labels=FACTOR_DISPLAY,
            method_labels={m: METHOD_DISPLAY.get(m, m) for m in dfs},
            title=f"Magnitude-varying sensitivity of drought-hazard "
                  f"characteristics to {axis_label} (method comparison)",
            magnitude_label=magnitude_label,
            output_path=fig_dir / "method_panel.pdf",
        )

    print(f"[plots/09/run_mv_sa] figures -> {fig_dir}")


if __name__ == "__main__":
    main()
