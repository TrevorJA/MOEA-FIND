"""Stage 07 / satisficing_sweep -- cheap-to-rerun manifest sweep.

Reads a pre-computed metric bank (produced by stage 06's policy
re-evaluation) plus the Pareto archive's drought characteristics (from
stage 04), runs every entry in a satisficing manifest through the GBT
classifier, and writes a summary table plus per-definition boundary
PDFs.

Does not re-run Pywr-DRB. To add new satisficing definitions, edit
``workflows/07_scenario_discovery/satisficing_manifest.yaml`` and
re-invoke this script against the same metric bank.

Outputs under ``outputs/07_scenario_discovery/satisficing_sweep/<slug>/``:
    - <primary_label>/classifiers/...
    - <primary_label>/figures/boundary_<def>.pdf
    - [<baseline_label>/...] (if --baseline-bank given)
    - classifier_summary.csv
    - figures/manifest_summary.pdf

The default slug mirrors the Pareto archive's variant (the directory
name of ``--pareto-results``'s parent).

Usage:
    python workflows/07_scenario_discovery/satisficing_sweep.py \\
        --pareto-results outputs/04_moea_find_single_site/run_moea_find/<slug>/results.json \\
        --bank outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/results/metric_bank.parquet \\
        --manifest workflows/07_scenario_discovery/satisficing_manifest.yaml
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Tuple

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.satisficing_labels import (  # noqa: E402
    load_manifest,
    sweep_manifest,
)
from src.plotting.satisficing_boundary import (  # noqa: E402
    plot_gbt_boundary_2d,
    plot_manifest_summary,
)
from src.io import load_metric_bank as _read_bank  # noqa: E402
from src.io import load_pareto_chars as _load_chars  # noqa: E402
from src.paths import stage_output_dir  # noqa: E402

STAGE = "07_scenario_discovery"
DRIVER = "satisficing_sweep"


def _run_one(
    bank: pd.DataFrame,
    chars: pd.DataFrame,
    manifest,
    feature_cols,
    output_dir: Path,
    seed: int,
    source_label: str,
) -> Tuple[pd.DataFrame, Path]:
    """Run the sweep for one (bank, chars) pair and emit figures."""
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = sweep_manifest(
        bank_df=bank, chars_df=chars, manifest=manifest,
        feature_cols=list(feature_cols),
        output_dir=output_dir, seed=seed,
    )
    summary.insert(0, "source", source_label)

    import joblib
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    for _, row in summary.iterrows():
        if row["status"] in ("no_labels", "degenerate_labels"):
            continue
        def_dir = output_dir / "classifiers" / row["definition_id"]
        model_path = def_dir / "model.joblib"
        if not model_path.exists():
            continue
        model = joblib.load(model_path)
        real_ids = bank.index.astype(str)
        keep_ids = [r for r in real_ids if r in chars.index]
        X = chars.loc[keep_ids, list(feature_cols)].astype(float)

        metric_series = pd.to_numeric(bank.loc[keep_ids, row["metric"]],
                                      errors="coerce")
        direction = row["direction"]
        threshold = float(row["threshold"])
        if direction == "le":
            y = (metric_series <= threshold)
        elif direction == "lt":
            y = (metric_series < threshold)
        elif direction == "ge":
            y = (metric_series >= threshold)
        else:
            y = (metric_series > threshold)
        mask = metric_series.notna()
        X_fit = X.loc[mask]
        y_fit = y.loc[mask].astype(int).values

        try:
            plot_gbt_boundary_2d(
                model=model, X=X_fit, y=y_fit,
                x_col=feature_cols[0],
                y_col=feature_cols[1] if len(feature_cols) > 1 else feature_cols[0],
                agg_col=feature_cols[2] if len(feature_cols) > 2 else None,
                output_path=fig_dir / f"boundary_{row['definition_id']}.pdf",
                title=f"{source_label} -- {row['definition_id']} (AUC={row['auc_mean']:.3f})",
            )
        except Exception as exc:
            print(f"[07/satisficing_sweep] boundary plot failed for "
                  f"{row['definition_id']}: {exc}")

    return summary, fig_dir


def _resolve_bank(args) -> Path:
    if args.bank is not None:
        return args.bank
    src_slug = args.src_slug or args.pareto_results.parent.name
    return (PROJECT_ROOT / "outputs" / "06_pywrdrb_reeval" / "policy_reeval"
            / src_slug / "results" / "metric_bank.parquet")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--pareto-results", type=Path, required=True,
                   help="Path to stage 04's results.json (Pareto archive).")
    p.add_argument("--bank", type=Path, default=None,
                   help="Primary metric_bank.parquet. Default: derived from "
                        "stage 06 outputs using --src-slug or pareto parent.")
    p.add_argument("--src-slug", type=str, default=None,
                   help="Override stage 06 source slug (when stage 06 was "
                        "re-evaluated under a different slug than stage 04).")
    p.add_argument("--manifest", type=Path,
                   default=PROJECT_ROOT / "workflows" / "07_scenario_discovery"
                   / "satisficing_manifest.yaml",
                   help="YAML satisficing manifest.")
    p.add_argument("--slug", type=str, default=None,
                   help="Output slug. Default: pareto archive parent dir name.")
    p.add_argument("--feature-cols", nargs="+",
                   default=("mean_duration", "mean_avg_severity",
                            "peak_severity_month"),
                   help="Drought-characteristic feature columns.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--baseline-bank", type=Path, default=None,
                   help="Optional second metric bank for side-by-side "
                        "comparison (e.g. library baseline).")
    p.add_argument("--baseline-chars", type=Path, default=None,
                   help="Drought characteristics JSON for the baseline.")
    p.add_argument("--baseline-label", default="baseline")
    p.add_argument("--primary-label", default="moea_find")
    args = p.parse_args()

    slug = args.slug or args.pareto_results.parent.name
    out = stage_output_dir(STAGE, DRIVER, slug)

    manifest = load_manifest(args.manifest)
    print(f"[07/satisficing_sweep] manifest: {len(manifest)} definitions")

    bank_path = _resolve_bank(args)
    primary_bank = _read_bank(bank_path)
    primary_chars = _load_chars(args.pareto_results, args.feature_cols)
    print(f"[07/satisficing_sweep] bank: {bank_path}")
    print(f"[07/satisficing_sweep] primary bank: {len(primary_bank)} "
          f"realizations, chars: {len(primary_chars)} records")

    primary_out = out / args.primary_label
    primary_summary, _ = _run_one(
        primary_bank, primary_chars, manifest,
        args.feature_cols, primary_out, args.seed, args.primary_label,
    )

    combined = [primary_summary]

    if args.baseline_bank is not None:
        if args.baseline_chars is None:
            raise SystemExit("--baseline-bank requires --baseline-chars")
        baseline_bank = _read_bank(args.baseline_bank)
        baseline_chars = _load_chars(args.baseline_chars, args.feature_cols)
        print(f"[07/satisficing_sweep] baseline bank: {len(baseline_bank)} "
              f"realizations, chars: {len(baseline_chars)} records")

        baseline_out = out / args.baseline_label
        baseline_summary, _ = _run_one(
            baseline_bank, baseline_chars, manifest,
            args.feature_cols, baseline_out, args.seed, args.baseline_label,
        )
        combined.append(baseline_summary)

    summary_all = pd.concat(combined, ignore_index=True)
    summary_path = out / "classifier_summary.csv"
    summary_all.to_csv(summary_path, index=False)
    print(f"[07/satisficing_sweep] wrote {summary_path}")

    plot_manifest_summary(
        summary_all.assign(
            definition_id=lambda d: d["source"] + ":" + d["definition_id"]
        ),
        output_path=out / "figures" / "manifest_summary.pdf",
        title="Satisficing sweep -- AUC per (source, definition)",
    )
    print(f"[07/satisficing_sweep] outputs: {out}")
    print(summary_all.to_string(index=False))


if __name__ == "__main__":
    main()
