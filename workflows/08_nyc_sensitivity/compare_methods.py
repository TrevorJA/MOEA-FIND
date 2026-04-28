"""Stage 08 — compare SA across runs.

Reads multiple ``run_sa.py`` output directories (different upstream
MOEA-FIND archives, different metric sets, or different SA-method
hyperparameters) and produces side-by-side comparisons.

Use cases:
- Sensitivity-of-method: did factor rankings change with the upstream
  metric set (e.g., ``primary`` vs ``extreme_event`` MOEA-FIND ablation)?
- Robustness-of-SA: did indices stabilize at the production sample size
  vs a smaller pilot run on the same upstream archive?

The script does not re-run SA — every upstream directory must already
contain ``results/indices_<method>.parquet``.

Usage:
    python workflows/08_nyc_sensitivity/compare_methods.py \\
        --runs outputs/exp10_nyc_sensitivity/<slug_a> \\
                outputs/exp10_nyc_sensitivity/<slug_b> \\
        --labels primary extreme_event \\
        --method delta \\
        --output-dir figures/nyc_sensitivity/comparisons
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.plotting.style import apply_style  # noqa: E402
from src.sensitivity import HEADLINE_INDEX  # noqa: E402


def _load_indices(run_dir: Path, method: str) -> pd.DataFrame:
    path = run_dir / "results" / f"indices_{method}.parquet"
    if not path.exists():
        raise FileNotFoundError(f"missing {path}")
    return pd.read_parquet(path)


def _common_outcomes(frames: Sequence[pd.DataFrame]) -> List[str]:
    sets = [set(df["outcome"].unique()) for df in frames]
    return sorted(set.intersection(*sets))


def _common_factors(frames: Sequence[pd.DataFrame]) -> List[str]:
    sets = [set(df["factor"].unique()) for df in frames]
    return sorted(set.intersection(*sets))


def _plot_side_by_side_tornado(
    *,
    run_indices: Dict[str, pd.DataFrame],
    outcome: str,
    method: str,
    output_path: Path,
) -> None:
    """One panel per run, sharing y-axis (factors)."""
    import matplotlib.pyplot as plt

    apply_style()
    headline = HEADLINE_INDEX[method]
    labels = list(run_indices.keys())
    factors = _common_factors([run_indices[l] for l in labels])

    fig, axes = plt.subplots(
        1, len(labels),
        figsize=(2.6 * len(labels) + 1.5, 0.4 * len(factors) + 1.5),
        sharey=True,
    )
    if len(labels) == 1:
        axes = [axes]

    # Order factors by the first run's headline index so the comparison
    # panels share a stable y-axis.
    first_df = run_indices[labels[0]]
    first_df_oc = first_df.query("outcome == @outcome").set_index("factor")
    if first_df_oc.empty:
        print(f"[compare] outcome {outcome!r} not in first run; skipping")
        plt.close(fig)
        return
    factor_order = (
        first_df_oc.loc[factors, "headline_index"]
        .sort_values(ascending=True).index.tolist()
    )
    y = np.arange(len(factor_order))

    cmap = plt.cm.tab10
    for k, (label, ax) in enumerate(zip(labels, axes)):
        df = run_indices[label]
        sub = df.query("outcome == @outcome").set_index("factor")
        if sub.empty:
            ax.set_title(f"{label}\n(missing)")
            continue
        vals = sub.reindex(factor_order)["headline_index"].values
        ci_lo = sub.reindex(factor_order).get("ci_lo")
        ci_hi = sub.reindex(factor_order).get("ci_hi")
        ax.barh(y, vals, color=cmap(k % 10), edgecolor="black",
                linewidth=0.5, alpha=0.85)
        if ci_lo is not None and ci_hi is not None:
            err_lo = np.maximum(0.0, vals - ci_lo.values)
            err_hi = np.maximum(0.0, ci_hi.values - vals)
            ax.errorbar(vals, y, xerr=[err_lo, err_hi], fmt="none",
                        ecolor="black", capsize=2.0, linewidth=0.7)
        ax.axvline(0.0, color="black", linewidth=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels(factor_order)
        ax.set_xlabel(headline)
        ax.set_title(label)

    fig.suptitle(f"{method}: {outcome}", fontsize=11)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def _plot_run_rank_correlation(
    *,
    run_indices: Dict[str, pd.DataFrame],
    outcome: str,
    method: str,
    output_path: Path,
) -> None:
    """Run × run Spearman correlation of factor rankings on one outcome."""
    from scipy.stats import spearmanr

    from src.plotting.sensitivity import plot_rank_correlation

    apply_style()
    labels = list(run_indices.keys())
    factors = _common_factors([run_indices[l] for l in labels])
    rho = np.full((len(labels), len(labels)), np.nan, dtype=float)

    rankings: Dict[str, np.ndarray] = {}
    for label in labels:
        df = run_indices[label].query("outcome == @outcome").set_index("factor")
        if df.empty:
            rankings[label] = np.full(len(factors), np.nan)
            continue
        rankings[label] = df.reindex(factors)["headline_index"].values

    for i, li in enumerate(labels):
        for j, lj in enumerate(labels):
            if i == j:
                rho[i, j] = 1.0
                continue
            r, _ = spearmanr(rankings[li], rankings[lj])
            rho[i, j] = r if np.isfinite(r) else np.nan

    rho_df = pd.DataFrame(rho, index=labels, columns=labels)
    plot_rank_correlation(
        rho_df,
        title=f"Run rank correlation: {method} / {outcome}",
        output_path=output_path,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--runs", nargs="+", type=Path, required=True,
                   help="Stage-08 output directories to compare.")
    p.add_argument("--labels", nargs="+", default=None,
                   help="Display labels for each run (default: directory name).")
    p.add_argument("--method", default="delta",
                   choices=sorted(HEADLINE_INDEX.keys()),
                   help="SA method to compare across runs.")
    p.add_argument("--outcomes", nargs="+", default=None,
                   help="Outcomes to compare. Defaults to the intersection "
                        "of outcomes present in every run.")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "figures" / "nyc_sensitivity" / "comparisons")
    args = p.parse_args()

    if args.labels is None:
        labels = [r.name for r in args.runs]
    else:
        if len(args.labels) != len(args.runs):
            raise SystemExit("--labels must have the same length as --runs")
        labels = list(args.labels)

    print(f"[compare] loading {len(args.runs)} runs for method={args.method}")
    run_indices: Dict[str, pd.DataFrame] = {}
    for label, run_dir in zip(labels, args.runs):
        run_indices[label] = _load_indices(Path(run_dir), args.method)
        print(f"[compare]   {label}: {len(run_indices[label])} rows from {run_dir}")

    outcomes = args.outcomes or _common_outcomes(list(run_indices.values()))
    print(f"[compare] comparing on outcomes: {outcomes}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for outcome in outcomes:
        _plot_side_by_side_tornado(
            run_indices=run_indices, outcome=outcome, method=args.method,
            output_path=args.output_dir / f"compare_tornado_{args.method}_{outcome}.pdf",
        )
        if len(run_indices) >= 2:
            _plot_run_rank_correlation(
                run_indices=run_indices, outcome=outcome, method=args.method,
                output_path=args.output_dir / f"compare_rho_{args.method}_{outcome}.pdf",
            )

    print(f"[compare] wrote figures to {args.output_dir}")


if __name__ == "__main__":
    main()
