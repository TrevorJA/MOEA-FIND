"""Script 12 — Cheap-to-rerun satisficing sweep.

Reads a pre-computed metric bank (produced by script 09 Stage 4) plus
the Pareto archive's drought characteristics (from script 04), runs
every entry in a satisficing manifest through the GBT classifier, and
writes a summary table plus per-definition boundary PDFs.

Does not re-run Pywr-DRB. To add new satisficing definitions, edit
``workflows/experiments/satisficing_manifest.yaml`` and re-invoke this
script against the same metric bank.

Usage:
    python workflows/experiments/12_satisficing_sweep.py \\
        --bank outputs/exp09_drb_policy_reeval/<slug>/results/metric_bank.parquet \\
        --chars outputs/exp04_kirsch_single_site/<slug>/results.json \\
        --manifest workflows/experiments/satisficing_manifest.yaml \\
        --output-dir outputs/exp09_drb_policy_reeval/<slug>/satisficing

    # Side-by-side comparison of two banks (e.g. MOEA-FIND vs library baseline)
    python workflows/experiments/12_satisficing_sweep.py \\
        --bank outputs/exp09_drb_policy_reeval/<moea_slug>/results/metric_bank.parquet \\
        --chars outputs/exp04_kirsch_single_site/<moea_slug>/results.json \\
        --baseline-bank outputs/exp09_drb_policy_reeval/<lhs_slug>/results/metric_bank.parquet \\
        --baseline-chars outputs/exp06_library_subsample/<lhs_slug>/selection.json \\
        --manifest workflows/experiments/satisficing_manifest.yaml \\
        --output-dir outputs/exp09_drb_policy_reeval/<moea_slug>/satisficing
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
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


def _read_bank(path: Path) -> pd.DataFrame:
    """Read a metric bank, handling parquet and CSV transparently."""
    path = Path(path)
    if path.suffix == ".parquet":
        try:
            return pd.read_parquet(path).assign(
                realization_id=lambda d: d.index.astype(str)
            ).set_index("realization_id", drop=True)
        except Exception:
            csv_alt = path.with_suffix(".csv")
            if csv_alt.exists():
                print(f"[12] parquet read failed; falling back to {csv_alt}")
                return pd.read_csv(csv_alt, index_col="realization_id",
                                   dtype={"realization_id": str})
            raise
    return pd.read_csv(path, index_col="realization_id",
                       dtype={"realization_id": str})


def _load_chars(chars_path: Path, feature_cols) -> pd.DataFrame:
    """Load per-realization drought characteristics from a script-04
    ``results.json`` or a script-06 selection JSON with ``pareto_chars``
    or ``selected_chars`` field."""
    chars_path = Path(chars_path)
    payload = json.loads(chars_path.read_text())
    rows = (
        payload.get("pareto_chars")
        or payload.get("selected_chars")
        or []
    )
    if not rows:
        raise ValueError(
            f"{chars_path} has neither 'pareto_chars' nor 'selected_chars'."
        )
    df = pd.DataFrame(rows)
    df.index = df.index.astype(str)
    df.index.name = "realization_id"
    missing = [c for c in feature_cols if c not in df.columns]
    if missing:
        raise KeyError(
            f"drought characteristics in {chars_path} missing columns {missing}. "
            f"Available: {list(df.columns)}"
        )
    return df


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

    # Per-definition boundary figure
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
        # Reconstruct the X/y used for this definition by re-applying the label
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
                title=f"{source_label} — {row['definition_id']} (AUC={row['auc_mean']:.3f})",
            )
        except Exception as exc:
            print(f"[12] boundary plot failed for {row['definition_id']}: {exc}")

    return summary, fig_dir


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--bank", type=Path, required=True,
                   help="Primary metric_bank.parquet (or .csv).")
    p.add_argument("--chars", type=Path, required=True,
                   help="Pareto archive results.json (has 'pareto_chars').")
    p.add_argument("--manifest", type=Path, required=True,
                   help="YAML satisficing manifest.")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Output directory for summary + figures.")
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

    manifest = load_manifest(args.manifest)
    print(f"[12] manifest: {len(manifest)} definitions")

    primary_bank = _read_bank(args.bank)
    primary_chars = _load_chars(args.chars, args.feature_cols)
    print(f"[12] primary bank: {len(primary_bank)} realizations, "
          f"chars: {len(primary_chars)} records")

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

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
        print(f"[12] baseline bank: {len(baseline_bank)} realizations, "
              f"chars: {len(baseline_chars)} records")

        baseline_out = out / args.baseline_label
        baseline_summary, _ = _run_one(
            baseline_bank, baseline_chars, manifest,
            args.feature_cols, baseline_out, args.seed, args.baseline_label,
        )
        combined.append(baseline_summary)

    summary_all = pd.concat(combined, ignore_index=True)
    summary_path = out / "classifier_summary.csv"
    summary_all.to_csv(summary_path, index=False)
    print(f"[12] wrote {summary_path}")

    # Aggregate AUC bar chart covering all sources
    plot_manifest_summary(
        summary_all.assign(
            definition_id=lambda d: d["source"] + ":" + d["definition_id"]
        ),
        output_path=out / "figures" / "manifest_summary.pdf",
        title="Satisficing sweep — AUC per (source, definition)",
    )
    print(f"[12] summary head:")
    print(summary_all.to_string(index=False))


if __name__ == "__main__":
    main()
