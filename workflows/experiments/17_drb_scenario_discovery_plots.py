"""Script 17 — DRB scenario-discovery figure generator.

Reads the persisted outputs of Script 09 (metric bank, drought levels,
Pareto archive) and produces every scenario-discovery figure we know
how to make: the legacy single-classifier satisficing map, plus per-
classifier multi-panel maps and per-model AUC/Brier summary bars for
``gbt`` and ``logreg``.

No Pywr-DRB, no MPI.  Safe to rerun whenever the metric bank changes,
or to iterate on the manifest / classifier recipes without touching
the simulation stack.

Usage:
    python workflows/experiments/17_drb_scenario_discovery_plots.py \\
        --pareto-results outputs/exp04_kirsch_single_site/{variant}/results.json

Assumes exp09 has already written
``outputs/exp09_drb_policy_reeval/{variant}/results/{metric_bank.{parquet,csv},
drought_levels.json}``.
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

from src.scenario_discovery import (  # noqa: E402
    build_satisficing_table,
    plot_satisficing_map,
    plot_satisficing_map_multi,
    save_results,
)
from src.satisficing_labels import (  # noqa: E402
    apply_labels,
    load_manifest,
    sweep_manifest,
)
from src.plotting.satisficing_boundary import plot_manifest_summary  # noqa: E402

EXP09_SLUG_ROOT = "exp09_drb_policy_reeval"
DEFAULT_MANIFEST = (
    PROJECT_ROOT / "workflows" / "experiments" / "satisficing_manifest.yaml"
)

CLASSIFIER_DESCRIPTIONS = {
    "gbt": "sklearn GradientBoostingClassifier, default hyperparams, "
           "5-fold stratified CV",
    "logreg": "sklearn LogisticRegression with degree-2 polynomial "
              "features + StandardScaler, 5-fold stratified CV",
}


def _read_metric_bank(results_dir: Path) -> pd.DataFrame:
    """Prefer parquet; fall back to CSV (what the bridge writes when
    pyarrow isn't installed)."""
    parquet = results_dir / "metric_bank.parquet"
    if parquet.exists():
        return pd.read_parquet(parquet)
    csv = results_dir / "metric_bank.csv"
    if csv.exists():
        return pd.read_csv(csv, index_col=0)
    raise FileNotFoundError(
        f"No metric_bank.{{parquet,csv}} under {results_dir}. "
        "Run exp09 Stage 4 classify first."
    )


def main():
    p = argparse.ArgumentParser(
        description="Generate DRB scenario-discovery figures from exp09 outputs."
    )
    p.add_argument(
        "--pareto-results", type=Path, required=True,
        help="Path to exp04's results.json — the Pareto archive that "
             "fed the exp09 re-evaluation.",
    )
    p.add_argument(
        "--exp09-dir", type=Path, default=None,
        help="Override exp09 variant dir.  Default: "
             f"outputs/{EXP09_SLUG_ROOT}/<slug-from-pareto-results-parent>.",
    )
    p.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST,
        help=f"Satisficing manifest YAML.  Default: {DEFAULT_MANIFEST}.",
    )
    p.add_argument(
        "--classifier-models", nargs="+", default=["gbt", "logreg"],
        choices=["gbt", "logreg"],
        help="Which classifiers to fit and plot (default: both).",
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

    # --- Resolve exp09 dir ----------------------------------------------------
    # Derive the slug from the Pareto results path rather than from the
    # JSON contents: the on-disk dir name is the authoritative slug used
    # by exp09/10 too, and parsing it avoids rebuilding the slug from
    # make_variant_slug which takes fields not always echoed in results.json.
    if args.exp09_dir is not None:
        exp09_dir = args.exp09_dir
    else:
        slug = args.pareto_results.parent.name
        exp09_dir = PROJECT_ROOT / "outputs" / EXP09_SLUG_ROOT / slug
    if not exp09_dir.exists():
        raise SystemExit(f"[17] exp09 dir not found: {exp09_dir}")

    results_dir = exp09_dir / "results"
    fig_dir = exp09_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"[17] pareto-results: {args.pareto_results}")
    print(f"[17] exp09 dir:      {exp09_dir.relative_to(PROJECT_ROOT)}")
    print(f"[17] n_pareto:       {n_pareto}")
    print(f"[17] objective keys: {objective_keys}")
    print(f"[17] manifest:       {args.manifest}")
    print(f"[17] classifiers:    {args.classifier_models}")

    # --- Load persisted artifacts --------------------------------------------
    bank = _read_metric_bank(results_dir)
    print(f"[17] metric bank: {bank.shape}")

    drought_levels_path = results_dir / "drought_levels.json"
    if not drought_levels_path.exists():
        raise SystemExit(
            f"[17] missing {drought_levels_path} — run exp09 Stage 4 first."
        )
    drought_levels = json.loads(drought_levels_path.read_text())

    # --- Legacy single-panel map ---------------------------------------------
    df = build_satisficing_table(
        drought_levels, pareto_chars, drought_metrics, objective_keys,
    )
    # Re-persist the legacy classification table — cheap, keeps
    # results_dir self-consistent if pareto_chars changed since exp09 ran.
    save_results(df, drought_levels, results_dir)

    legacy_fig = fig_dir / "fig09_satisficing_map.pdf"
    plot_satisficing_map(df, anti_ideal=anti_ideal, output_path=legacy_fig)
    print(f"[17] wrote {legacy_fig.relative_to(PROJECT_ROOT)}")

    # --- Manifest sweep per classifier ---------------------------------------
    manifest = load_manifest(args.manifest)
    chars_records = [dict(c) for c in pareto_chars]
    chars_df = pd.DataFrame(chars_records)
    chars_df.index = realization_ids
    chars_df.index.name = "realization_id"
    feature_cols = tuple(objective_keys)
    labels_long = apply_labels(bank, manifest)

    slug_label = exp09_dir.name

    for model_name in args.classifier_models:
        print(f"[17] --- classifier: {model_name} ---")
        model_dir = results_dir / model_name
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
        print(f"[17]   wrote {summary_path.relative_to(PROJECT_ROOT)}")

        if summary_df.empty:
            print(f"[17]   {model_name} summary empty; skipping plots")
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
        print(f"[17]   wrote {multi_fig.relative_to(PROJECT_ROOT)}")

        auc_fig = fig_dir / f"fig_manifest_summary_{model_name}.pdf"
        plot_manifest_summary(
            summary_df,
            output_path=auc_fig,
            title=f"Satisficing sweep ({model_name}) — {slug_label}",
        )
        print(f"[17]   wrote {auc_fig.relative_to(PROJECT_ROOT)}")

    print(f"\n[17] scenario-discovery figures complete in "
          f"{fig_dir.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
