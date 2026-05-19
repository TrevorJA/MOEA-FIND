"""Stage 07 / scenario_discovery_plots -- final SD figure generator.

Reads the persisted outputs of stage 06 (metric bank, drought levels)
and the Pareto archive from stage 04 (drought characteristics) and
produces every scenario-discovery figure: the legacy single-classifier
satisficing map, plus per-classifier multi-panel maps and per-model
AUC/Brier summary bars for ``gbt`` and ``logreg``.

No Pywr-DRB, no MPI. Safe to rerun whenever the metric bank changes,
or to iterate on the manifest / classifier recipes.

Numerical derivatives (re-fit classifiers, refreshed satisficing
tables) land under
``outputs/07_scenario_discovery/scenario_discovery_plots/<slug>/``.
Figures land under
``figures/07_scenario_discovery/scenario_discovery_plots/<slug>/``.

Inputs:
    --pareto-results
        outputs/04_moea_find_single_site/run_moea_find/<slug>/results.json
    Stage 06 metric bank derived from the Pareto slug (or --src-slug):
        outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/results/metric_bank.parquet
        outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/results/drought_levels.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.discovery.scenario_discovery import (  # noqa: E402
    build_satisficing_table,
    plot_satisficing_map,
    plot_satisficing_map_multi,
    save_results,
)
from src.discovery.satisficing_labels import (  # noqa: E402
    apply_labels,
    load_manifest,
    sweep_manifest,
)
from src.plotting.satisficing_boundary import plot_manifest_summary  # noqa: E402
from src.io_paths.paths import stage_output_dir, stage_figure_dir  # noqa: E402

STAGE = "07_scenario_discovery"
DRIVER = "scenario_discovery_plots"
STAGE06 = "06_pywrdrb_reeval"
STAGE06_DRIVER = "policy_reeval"

DEFAULT_MANIFEST = (
    PROJECT_ROOT / "workflows" / "07_scenario_discovery"
    / "satisficing_manifest.yaml"
)

CLASSIFIER_DESCRIPTIONS = {
    "gbt": "sklearn GradientBoostingClassifier, default hyperparams, "
           "5-fold stratified CV",
    "logreg": "sklearn LogisticRegression with degree-2 polynomial "
              "features + StandardScaler, 5-fold stratified CV",
}


def _read_metric_bank(results_dir: Path) -> pd.DataFrame:
    """Locate ``metric_bank.{parquet,csv}`` under ``results_dir``."""
    from src.io_paths.io import load_metric_bank
    return load_metric_bank(results_dir / "metric_bank.parquet")


def main():
    p = argparse.ArgumentParser(
        description="Generate DRB scenario-discovery figures from stage 06 outputs."
    )
    p.add_argument(
        "--pareto-results", type=Path, required=True,
        help="Path to stage 04's results.json (Pareto archive that fed "
             "the stage 06 re-evaluation).",
    )
    p.add_argument(
        "--src-slug", type=str, default=None,
        help="Stage 06 source slug. Default: parent dir name of "
             "--pareto-results.",
    )
    p.add_argument(
        "--stage06-dir", type=Path, default=None,
        help="Override the full stage 06 variant dir. Default: "
             f"outputs/{STAGE06}/{STAGE06_DRIVER}/<src_slug>.",
    )
    p.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST,
        help=f"Satisficing manifest YAML. Default: {DEFAULT_MANIFEST}.",
    )
    p.add_argument(
        "--classifier-models", nargs="+", default=["gbt", "logreg"],
        choices=["gbt", "logreg"],
        help="Which classifiers to fit and plot (default: both).",
    )
    p.add_argument(
        "--slug", type=str, default=None,
        help="Output slug. Default: parent dir name of --pareto-results.",
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Shared seed for classifier CV and fold generation.",
    )
    args = p.parse_args()

    # --- Load Pareto archive --------------------------------------------------
    results = json.loads(args.pareto_results.read_text())
    pareto_chars = results.get("pareto_chars") or []
    drought_metrics = np.asarray(
        results.get("drought_metrics") or [], dtype=float
    )
    objective_keys = list(results.get("objective_keys") or [])
    anti_ideal = np.asarray(results.get("anti_ideal") or [], dtype=float)
    n_pareto = int(results.get("n_pareto", len(pareto_chars)))
    realization_ids = [str(i) for i in range(n_pareto)]

    # --- Resolve stage 06 dir -------------------------------------------------
    src_slug = args.src_slug or args.pareto_results.parent.name
    if args.stage06_dir is not None:
        stage06_dir = args.stage06_dir
    else:
        stage06_dir = (PROJECT_ROOT / "outputs" / STAGE06
                       / STAGE06_DRIVER / src_slug)
    if not stage06_dir.exists():
        raise SystemExit(f"[07/scenario_discovery_plots] stage 06 dir not "
                         f"found: {stage06_dir}")

    src_results_dir = stage06_dir / "results"

    # --- Resolve our own output / figure dirs --------------------------------
    slug = args.slug or args.pareto_results.parent.name
    out_dir = stage_output_dir(STAGE, DRIVER, slug)
    fig_dir = stage_figure_dir(STAGE, DRIVER, slug)

    print(f"[07/scenario_discovery_plots] pareto-results: {args.pareto_results}")
    print(f"[07/scenario_discovery_plots] stage 06 dir:   "
          f"{stage06_dir.relative_to(PROJECT_ROOT)}")
    print(f"[07/scenario_discovery_plots] outputs:        "
          f"{out_dir.relative_to(PROJECT_ROOT)}")
    print(f"[07/scenario_discovery_plots] figures:        "
          f"{fig_dir.relative_to(PROJECT_ROOT)}")
    print(f"[07/scenario_discovery_plots] n_pareto:       {n_pareto}")
    print(f"[07/scenario_discovery_plots] objective keys: {objective_keys}")
    print(f"[07/scenario_discovery_plots] manifest:       {args.manifest}")
    print(f"[07/scenario_discovery_plots] classifiers:    {args.classifier_models}")

    # --- Load persisted artifacts --------------------------------------------
    bank = _read_metric_bank(src_results_dir)
    print(f"[07/scenario_discovery_plots] metric bank: {bank.shape}")

    drought_levels_path = src_results_dir / "drought_levels.json"
    if not drought_levels_path.exists():
        raise SystemExit(
            f"[07/scenario_discovery_plots] missing {drought_levels_path} -- "
            f"run stage 06 first."
        )
    drought_levels = json.loads(drought_levels_path.read_text())

    # --- Legacy single-panel map ---------------------------------------------
    df = build_satisficing_table(
        drought_levels, pareto_chars, drought_metrics, objective_keys,
    )
    # Persist the refreshed legacy classification table under our own
    # outputs dir (not under stage 06).
    save_results(df, drought_levels, out_dir)

    legacy_fig = fig_dir / "fig09_satisficing_map.pdf"
    plot_satisficing_map(df, anti_ideal=anti_ideal, output_path=legacy_fig)
    print(f"[07/scenario_discovery_plots] wrote "
          f"{legacy_fig.relative_to(PROJECT_ROOT)}")

    # --- Manifest sweep per classifier ---------------------------------------
    manifest = load_manifest(args.manifest)
    chars_records = [dict(c) for c in pareto_chars]
    chars_df = pd.DataFrame(chars_records)
    chars_df.index = realization_ids
    chars_df.index.name = "realization_id"
    feature_cols = tuple(objective_keys)
    labels_long = apply_labels(bank, manifest)

    for model_name in args.classifier_models:
        print(f"[07/scenario_discovery_plots] --- classifier: {model_name} ---")
        model_dir = out_dir / model_name
        summary_df = sweep_manifest(
            bank_df=bank, chars_df=chars_df, manifest=manifest,
            feature_cols=list(feature_cols),
            output_dir=model_dir,
            seed=args.seed,
            model_name=model_name,
        )
        summary_path = model_dir / "classifier_summary.csv"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(summary_path, index=False)
        print(f"[07/scenario_discovery_plots]   wrote "
              f"{summary_path.relative_to(PROJECT_ROOT)}")

        if summary_df.empty:
            print(f"[07/scenario_discovery_plots]   {model_name} summary "
                  f"empty; skipping plots")
            continue

        multi_fig = fig_dir / f"fig09_satisficing_map_{model_name}.pdf"
        plot_satisficing_map_multi(
            labels_long=labels_long,
            chars_df=chars_df,
            manifest=manifest,
            classifiers_dir=model_dir / "classifiers",
            feature_cols=list(feature_cols),
            anti_ideal=anti_ideal,
            classifier_label=CLASSIFIER_DESCRIPTIONS.get(
                model_name, model_name
            ),
            output_path=multi_fig,
        )
        print(f"[07/scenario_discovery_plots]   wrote "
              f"{multi_fig.relative_to(PROJECT_ROOT)}")

        auc_fig = fig_dir / f"fig_manifest_summary_{model_name}.pdf"
        plot_manifest_summary(
            summary_df,
            output_path=auc_fig,
            title=f"Satisficing sweep ({model_name}) -- {slug}",
        )
        print(f"[07/scenario_discovery_plots]   wrote "
              f"{auc_fig.relative_to(PROJECT_ROOT)}")

    print(f"\n[07/scenario_discovery_plots] complete.")
    print(f"  outputs: {out_dir.relative_to(PROJECT_ROOT)}")
    print(f"  figures: {fig_dir.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
