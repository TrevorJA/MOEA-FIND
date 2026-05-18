"""Exploratory drought-metric selection for MOEA-FIND.

Computes a ~28-metric candidate library on the historical T-year blocks,
runs spread/correlation diagnostics, and recommends K=3 and K=4 metric sets
to replace the current ``"primary"`` preset (where ``mean_severity`` and
``mean_magnitude`` waste a Borg dimension by being highly correlated).

**Selection is hazard-interpretable**: every recommendation is a tuple of
named, single drought characteristics (no PCA combinations). The driver
ranks K-sets by **min robust spread**, subject to three hard constraints:

1. Each metric is from a different |ρ_S|≥0.7 redundancy cluster.
2. Each metric is from a different hydrologic concept (severity, frequency,
   duration, recovery, deficit_volume, flow_tail, ...) per
   :data:`src.metrics.extended.CONCEPT_MAP`.
3. All pairwise |ρ_S| < ``--correlation-cap`` (default 0.6).

**Robust spread score**: ``IQR / (|median| + σ)``. Uses σ as a stable
fallback when |median| is small, so sign-crossing metrics like
``sen_slope_annual_min_neg`` no longer get artificially inflated rankings.

**Statistical framing**: stride-1 T-blocks share up to 19/20 years and are
NOT independent (effective n ≈ 3.7). The 54-block matrix is treated as a
deterministic descriptor of the historical record's behaviour under a
20-year lens; rank-order of correlations and cluster structure are
reported, absolute p-values are not.

Output: ``outputs/02_calibration/metric_explorer/`` —

* ``block_chars_extended.csv`` — per-block candidate matrix.
* ``per_metric_spread.csv`` — descriptors + robust spread + screen flag.
* ``pearson_corr.csv``, ``spearman_corr.csv`` — pairwise matrices.
* ``clusters.json`` — cluster ID + concept tag + representative per metric.
* ``k3_alternatives.json``, ``k4_alternatives.json`` — top-N enumerated sets.
* ``recommendation.md`` — top set as the primary recommendation, plus
  the next four alternatives with their trade-offs.
* ``config.json`` — reproducibility metadata.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.metrics.extended import CANDIDATE_METRIC_NAMES  # noqa: E402
from src.metrics.screening import (  # noqa: E402
    CLUSTER_DISTANCE_CUT,
    N_ALTERNATIVES,
    cluster_metrics,
    compute_block_matrix,
    correlation_matrices,
    per_metric_spread,
    relax_until_nonempty,
)
from src.metrics.objectives import flows_to_series, make_ssi_calculator  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "02_calibration"
DRIVER = "metric_explorer"


# ---------------------------------------------------------------------------
# Markdown helpers
# ---------------------------------------------------------------------------


def _df_to_markdown(df: pd.DataFrame) -> str:
    """Render a small DataFrame as a markdown table (no ``tabulate`` dep)."""
    cols = list(df.columns)
    header = "| " + " | ".join([""] + [str(c) for c in cols]) + " |"
    divider = "|" + "|".join(["---"] * (len(cols) + 1)) + "|"
    rows = [header, divider]
    for idx, row in df.iterrows():
        cells = [str(idx)] + [
            f"{v:.3f}" if isinstance(v, (int, float)) else str(v) for v in row
        ]
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join(rows)


def _alternatives_table(alts: List[Dict]) -> str:
    """Render top-N candidate K-sets as a markdown table."""
    if not alts:
        return "_No feasible K-set under the active constraints._"
    lines = [
        "| rank | metrics | concepts | min_spread | sum_spread | max \\|ρ_S\\| |",
        "|------|---------|----------|-----------:|-----------:|-----------:|",
    ]
    for i, a in enumerate(alts, 1):
        metrics = ", ".join(f"`{m}`" for m in a["metrics"])
        concepts = ", ".join(a["concepts"])
        lines.append(
            f"| {i} | {metrics} | {concepts} | "
            f"{a['min_spread']:.3f} | {a['sum_spread']:.3f} | "
            f"{a['max_pairwise_rho']:.3f} |"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Recommendation
# ---------------------------------------------------------------------------


def first_without_concept(
    alts: List[Dict], excluded_concept: str,
) -> Dict | None:
    """Return the first K-set in ``alts`` that does NOT include any metric
    tagged with ``excluded_concept``. ``None`` if no such set exists."""
    for a in alts:
        if excluded_concept not in a["concepts"]:
            return a
    return None


def best_balanced_without_concept(
    alts: List[Dict], excluded_concept: str,
) -> Dict | None:
    """Return the K-set (excluding ``excluded_concept``) maximising the
    composite score ``min_spread × (1 − max_pairwise_rho)``.

    Pure "lowest max |ρ_S|" picks degenerate sets where every metric is
    near-constant across T-blocks (their correlations are trivially low
    because their variance is tiny). The composite penalises both
    dimensions: a K-set with low spread or with high correlation both
    lose. Approximates the Pareto-optimal trade-off between axis
    informativeness and axis independence.
    """
    pool = [a for a in alts if excluded_concept not in a["concepts"]]
    if not pool:
        return None
    return max(
        pool,
        key=lambda a: a["min_spread"] * (1.0 - a["max_pairwise_rho"]),
    )


def write_recommendation_md(
    out_path: Path,
    *,
    config: Dict,
    spread_df: pd.DataFrame,
    cluster_df: pd.DataFrame,
    k3_alts: List[Dict],
    k3_constraint_rung: str,
    k4_alts: List[Dict],
    k4_constraint_rung: str,
    spearman_corr: pd.DataFrame,
) -> None:
    """Render the recommendation markdown report (interpretable hazards only)."""
    surviving = spread_df[spread_df["passes_screen"]]["metric"].tolist()
    n_total = len(spread_df)
    n_pass = len(surviving)
    n_clusters = int(cluster_df["cluster_id"].nunique())
    n_concepts = int(cluster_df["concept"].nunique())

    lines: List[str] = []
    lines.append("# MOEA-FIND drought-metric selection — recommendation")
    lines.append("")
    lines.append("This file is auto-generated by "
                 "`workflows/02_calibration/metric_explorer.py`.")
    lines.append("")
    lines.append(
        f"- T-block length: **{config['T_years']} yr**, stride **{config['stride']}**, "
        f"**{config['n_blocks']}** blocks.")
    lines.append("- SSI timescales: SSI-3 (production) and SSI-12 (long-term hydrologic).")
    lines.append(
        f"- Q80 threshold (fixed, full-record): **{config['q80_threshold']:.2f} cfs**.")
    lines.append(
        f"- Candidates: **{n_total}** total → **{n_pass}** survive spread/skew screen "
        f"→ **{n_clusters}** redundancy clusters at |ρ_S| ≥ 0.7 spanning "
        f"**{n_concepts}** distinct hydrologic concepts.")
    lines.append("")
    lines.append("Selection is **hazard-interpretable**: every metric in every "
                 "recommendation is a single, named drought characteristic — no "
                 "PCA combinations, no synthetic axes. Sets are scored by ")
    lines.append("**min(spread score)** across members (so no weak axis), "
                 "subject to three hard constraints: (1) distinct redundancy "
                 "clusters, (2) distinct hydrologic concepts (severity vs "
                 "frequency vs duration vs recovery vs deficit_volume vs "
                 "flow_tail vs ...), (3) pairwise |ρ_S| < "
                 f"**{config['correlation_cap']}**.")
    lines.append("")
    lines.append(
        "Spread score is **IQR / (|median| + σ)** — robust to sign-crossing "
        "metrics where the original ``IQR / |median|`` would explode.")
    lines.append("")

    # ---------------- K=3 ----------------
    lines.append("## Recommended K=3 set")
    lines.append("")
    if k3_alts:
        primary = k3_alts[0]
        hazard_only = first_without_concept(k3_alts, "trend")
        lines.append("```")
        lines.append("(" + ", ".join(f'"{m}"' for m in primary["metrics"]) + ")")
        lines.append("```")
        lines.append("")
        if k3_constraint_rung != "distinct_clusters_and_concepts":
            lines.append(
                f"> ⚠️  Constraint relaxation: no K=3 set satisfied both "
                f"distinct clusters AND distinct concepts; falling back to "
                f"`{k3_constraint_rung}`. Inspect the alternatives table "
                f"and consider whether two members share enough hydrologic "
                f"meaning to be redundant in practice.")
            lines.append("")
        lines.append(
            f"- min spread = **{primary['min_spread']:.3f}**, "
            f"sum spread = **{primary['sum_spread']:.3f}**, "
            f"max pairwise |ρ_S| = **{primary['max_pairwise_rho']:.3f}**.")
        lines.append(
            f"- concepts: {', '.join(primary['concepts'])}.")
        lines.append("")
        sub = spearman_corr.loc[primary["metrics"], primary["metrics"]].round(3)
        lines.append("Pairwise Spearman ρ:")
        lines.append("")
        lines.append(_df_to_markdown(sub))
        lines.append("")

        # Hazard-only alternative: best K=3 without any trend-concept metric.
        if hazard_only is not None and hazard_only["metrics"] != primary["metrics"]:
            lines.append("### Hazard-only alternatives (no trend metrics)")
            lines.append("")
            lines.append(
                "Sen's slope is a within-block trend, not a steady-state "
                "drought hazard. The two hazard-only options below differ in "
                "what they optimise:")
            lines.append("")
            lines.append("**Option A — highest min spread** "
                         "(every axis is informative across T-blocks):")
            lines.append("")
            lines.append("```")
            lines.append("(" + ", ".join(f'"{m}"' for m in hazard_only["metrics"]) + ")")
            lines.append("```")
            lines.append("")
            lines.append(
                f"- min spread = **{hazard_only['min_spread']:.3f}**, "
                f"sum spread = **{hazard_only['sum_spread']:.3f}**, "
                f"max pairwise |ρ_S| = **{hazard_only['max_pairwise_rho']:.3f}**.")
            lines.append(
                f"- concepts: {', '.join(hazard_only['concepts'])}.")
            lines.append("")
            sub_h = spearman_corr.loc[
                hazard_only["metrics"], hazard_only["metrics"]
            ].round(3)
            lines.append(_df_to_markdown(sub_h))
            lines.append("")

            balanced = best_balanced_without_concept(k3_alts, "trend")
            if balanced is not None and balanced["metrics"] != hazard_only["metrics"]:
                bal_score = balanced["min_spread"] * (1.0 - balanced["max_pairwise_rho"])
                lines.append(
                    "**Option B — best balanced** (maximises "
                    "``min_spread × (1 − max |ρ_S|)`` — penalises both "
                    "low informativeness and high correlation; "
                    "**recommended for MOEA spread of the Pareto front**):")
                lines.append("")
                lines.append("```")
                lines.append("(" + ", ".join(f'"{m}"' for m in balanced["metrics"]) + ")")
                lines.append("```")
                lines.append("")
                lines.append(
                    f"- min spread = **{balanced['min_spread']:.3f}**, "
                    f"sum spread = **{balanced['sum_spread']:.3f}**, "
                    f"max pairwise |ρ_S| = **{balanced['max_pairwise_rho']:.3f}**, "
                    f"composite score = **{bal_score:.3f}**.")
                lines.append(
                    f"- concepts: {', '.join(balanced['concepts'])}.")
                lines.append("")
                sub_o = spearman_corr.loc[
                    balanced["metrics"], balanced["metrics"]
                ].round(3)
                lines.append(_df_to_markdown(sub_o))
                lines.append("")

        lines.append(f"### Top {min(N_ALTERNATIVES, len(k3_alts))} alternatives")
        lines.append("")
        lines.append(_alternatives_table(k3_alts[:N_ALTERNATIVES]))
        lines.append("")
    else:
        lines.append("_No feasible K=3 set under any constraint rung._")
        lines.append("")

    # ---------------- K=4 ----------------
    lines.append("## Recommended K=4 set")
    lines.append("")
    if k4_alts:
        primary = k4_alts[0]
        hazard_only = first_without_concept(k4_alts, "trend")
        lines.append("```")
        lines.append("(" + ", ".join(f'"{m}"' for m in primary["metrics"]) + ")")
        lines.append("```")
        lines.append("")
        if k4_constraint_rung != "distinct_clusters_and_concepts":
            lines.append(
                f"> ⚠️  Constraint relaxation: no K=4 set satisfied both "
                f"distinct clusters AND distinct concepts; falling back to "
                f"`{k4_constraint_rung}`.")
            lines.append("")
        lines.append(
            f"- min spread = **{primary['min_spread']:.3f}**, "
            f"sum spread = **{primary['sum_spread']:.3f}**, "
            f"max pairwise |ρ_S| = **{primary['max_pairwise_rho']:.3f}**.")
        lines.append(
            f"- concepts: {', '.join(primary['concepts'])}.")
        lines.append("")
        sub = spearman_corr.loc[primary["metrics"], primary["metrics"]].round(3)
        lines.append("Pairwise Spearman ρ:")
        lines.append("")
        lines.append(_df_to_markdown(sub))
        lines.append("")

        if hazard_only is not None and hazard_only["metrics"] != primary["metrics"]:
            lines.append("### Hazard-only K=4 alternative (no trend metrics)")
            lines.append("")
            lines.append("```")
            lines.append("(" + ", ".join(f'"{m}"' for m in hazard_only["metrics"]) + ")")
            lines.append("```")
            lines.append("")
            lines.append(
                f"- min spread = **{hazard_only['min_spread']:.3f}**, "
                f"sum spread = **{hazard_only['sum_spread']:.3f}**, "
                f"max pairwise |ρ_S| = **{hazard_only['max_pairwise_rho']:.3f}**.")
            lines.append(
                f"- concepts: {', '.join(hazard_only['concepts'])}.")
            lines.append("")
            sub_h = spearman_corr.loc[
                hazard_only["metrics"], hazard_only["metrics"]
            ].round(3)
            lines.append(_df_to_markdown(sub_h))
            lines.append("")

        lines.append(f"### Top {min(N_ALTERNATIVES, len(k4_alts))} alternatives")
        lines.append("")
        lines.append(_alternatives_table(k4_alts[:N_ALTERNATIVES]))
        lines.append("")
    else:
        lines.append("_No feasible K=4 set under any constraint rung._")
        lines.append("")

    # ---------------- How to apply ----------------
    lines.append("## How to apply")
    lines.append("")
    lines.append("This file is a recommendation only; PRESETS / REGISTRY are "
                 "not auto-modified. To apply:")
    lines.append("")
    lines.append("1. For any metric in the recommended set that is not already "
                 "in `src.metrics.drought_metrics.REGISTRY`, add a `DroughtMetric` "
                 "instance with an appropriate `AntiIdealRule` "
                 "(HEADROOM_TIMES_MAX for unbounded-above metrics; CONSTANT "
                 "for fractions in [0, 1]; CYCLIC_HEADROOM only for cyclic-"
                 "month metrics — none in this recommendation).")
    lines.append("2. Define a new preset (suggest `'primary_v2'`) in "
                 "`src.metrics.drought_metrics.PRESETS` with the recommended tuple. "
                 "Keep the existing `'primary'` preset for reproducibility of "
                 "prior MOEA runs.")
    lines.append("3. Re-run `workflows/02_calibration/metric_blocks.py` to "
                 "produce the 3D scatter for the new preset.")
    lines.append("4. Update manuscript governance (DD-04) citing this report.")
    lines.append("")

    out_path.write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-years", type=int, default=20,
                   help="Length of each T-year block (default 20).")
    p.add_argument("--stride", type=int, default=1,
                   help="Block stride in years (default 1).")
    p.add_argument("--cache-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / "data_cache")
    p.add_argument("--correlation-cap", type=float, default=0.6,
                   help="Max |ρ_S| allowed between any two metrics in a "
                        "recommended set.")
    p.add_argument("--cluster-cut", type=float, default=CLUSTER_DISTANCE_CUT,
                   help="Average-linkage cut on 1−|ρ_S| (default 0.30).")
    p.add_argument("--metrics-subset", type=str, default=None,
                   help="Comma-separated subset to force-include in the pool.")
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER)
    print(f"[02/metric_explorer] output_dir={out_dir}")

    # --- Data ---
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(args.cache_dir)
    n_hist_years = monthly_2d.shape[0]
    n_blocks = (n_hist_years - args.T_years) // args.stride + 1
    print(f"[diag] historical: {n_hist_years} water years, "
          f"{n_blocks} blocks of T={args.T_years} yr (stride={args.stride})")

    # --- SSI calculators (pre-fit on full record) ---
    ssi3_calc = make_ssi_calculator(timescale=3)
    ssi3_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    ssi12_calc = make_ssi_calculator(timescale=12)
    ssi12_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))

    # --- Q80 threshold from full record ---
    q80 = float(np.percentile(monthly_1d, 20.0))
    print(f"[diag] Q80 (low-flow threshold, fixed): {q80:.2f} cfs")

    # --- Per-block candidate matrix ---
    print("[diag] computing per-block candidate matrix ...")
    chars_df = compute_block_matrix(
        monthly_1d, args.T_years, args.stride, ssi3_calc, ssi12_calc, q80
    )
    chars_df.to_csv(out_dir / "block_chars_extended.csv", index=False)
    print(f"[diag] wrote block_chars_extended.csv "
          f"({len(chars_df)} rows × {len(chars_df.columns)} cols)")

    # --- Per-metric spread + screening (robust spread score) ---
    spread_df = per_metric_spread(chars_df)
    spread_df.to_csv(out_dir / "per_metric_spread.csv", index=False)
    surviving = spread_df[spread_df["passes_screen"]]["metric"].tolist()
    print(f"[diag] {len(surviving)}/{len(spread_df)} metrics pass screen")

    if args.metrics_subset:
        forced = [m.strip() for m in args.metrics_subset.split(",") if m.strip()]
        for m in forced:
            if m in chars_df.columns and m not in surviving:
                surviving.append(m)

    # --- Correlations on surviving metrics ---
    pearson, spearman = correlation_matrices(chars_df, surviving)
    pearson.to_csv(out_dir / "pearson_corr.csv")
    spearman.to_csv(out_dir / "spearman_corr.csv")

    # --- Hierarchical clustering (representative = max robust spread) ---
    spread_lookup = spread_df.set_index("metric")["spread_score"]
    cluster_df = cluster_metrics(spearman, spread_lookup, args.cluster_cut)
    cluster_df.to_json(out_dir / "clusters.json", orient="records", indent=2)
    n_clusters = int(cluster_df["cluster_id"].nunique())
    n_concepts = int(cluster_df["concept"].nunique())
    print(f"[diag] {n_clusters} clusters / {n_concepts} distinct concepts on "
          f"{len(surviving)} surviving metrics")

    # --- Brute-force K-set enumeration with strict-then-relaxed ladder ---
    k3_alts, k3_rung = relax_until_nonempty(
        surviving, spearman, cluster_df, spread_df, K=3,
        correlation_cap=args.correlation_cap,
    )
    k4_alts, k4_rung = relax_until_nonempty(
        surviving, spearman, cluster_df, spread_df, K=4,
        correlation_cap=args.correlation_cap,
    )

    (out_dir / "k3_alternatives.json").write_text(
        json.dumps({"constraint_rung": k3_rung, "alternatives": k3_alts[:N_ALTERNATIVES]},
                   indent=2)
    )
    (out_dir / "k4_alternatives.json").write_text(
        json.dumps({"constraint_rung": k4_rung, "alternatives": k4_alts[:N_ALTERNATIVES]},
                   indent=2)
    )

    if k3_alts:
        print(f"[rec] K=3 ({k3_rung}): {k3_alts[0]['metrics']}")
        print(f"        min_spread={k3_alts[0]['min_spread']:.3f}, "
              f"max|ρ_S|={k3_alts[0]['max_pairwise_rho']:.3f}")
        print(f"[rec] K=3 considered {len(k3_alts)} feasible sets.")
    else:
        print("[rec] K=3: NO FEASIBLE SET")
    if k4_alts:
        print(f"[rec] K=4 ({k4_rung}): {k4_alts[0]['metrics']}")
        print(f"        min_spread={k4_alts[0]['min_spread']:.3f}, "
              f"max|ρ_S|={k4_alts[0]['max_pairwise_rho']:.3f}")

    # --- Reproducibility config ---
    cfg = {
        "script": "workflows/02_calibration/metric_explorer.py",
        "T_years": int(args.T_years),
        "stride": int(args.stride),
        "n_hist_years": int(n_hist_years),
        "n_blocks": int(len(chars_df)),
        "q80_threshold": q80,
        "correlation_cap": float(args.correlation_cap),
        "cluster_distance_cut": float(args.cluster_cut),
        "candidate_metric_names": list(CANDIDATE_METRIC_NAMES),
        "n_metrics_total": int(len(spread_df)),
        "n_metrics_pass_screen": int(len(surviving)),
        "n_clusters": n_clusters,
        "n_concepts": n_concepts,
        "k3_constraint_rung": k3_rung,
        "k4_constraint_rung": k4_rung,
        "n_k3_feasible": len(k3_alts),
        "n_k4_feasible": len(k4_alts),
    }
    (out_dir / "config.json").write_text(json.dumps(cfg, indent=2))

    # --- Recommendation markdown ---
    write_recommendation_md(
        out_dir / "recommendation.md",
        config=cfg,
        spread_df=spread_df,
        cluster_df=cluster_df,
        k3_alts=k3_alts,
        k3_constraint_rung=k3_rung,
        k4_alts=k4_alts,
        k4_constraint_rung=k4_rung,
        spearman_corr=spearman,
    )
    print(f"[diag] wrote {out_dir / 'recommendation.md'}")
    print("[diag] done.")


if __name__ == "__main__":
    main()
