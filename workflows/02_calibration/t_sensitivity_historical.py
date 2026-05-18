"""Stage 1 — T-sensitivity on historical T-year blocks (per-T driver).

For one value of ``T``, run the full 28-metric screening pipeline on
the overlapping stride-1 historical blocks of the Cannonsville record
and write per-T diagnostic tables. SSI-3, SSI-12, and Q80 calibrators
are fit **once** on the full historical record so every T sees the
same calibrated distributions (DD-11 lock-in).

Outputs (under ``outputs/02_calibration/t_sensitivity_historical/T{T:02d}/``):

* ``block_chars_extended.csv`` — per-block 28-metric matrix.
* ``per_metric_spread.csv`` — descriptors + robust spread + screen flag.
* ``degeneracy.csv`` — per-metric zero/saturated/NaN fractions plus
  event-count diagnostics. **Stage-1 headline table** for assessing
  small-T viability: a metric whose value is zero in >25% of blocks at
  T=5 cannot serve as a Borg objective there because Borg cannot move
  along an axis with a degenerate cluster.
* ``pearson_corr.csv``, ``spearman_corr.csv`` — correlation matrices on
  surviving metrics.
* ``clusters.json`` — hierarchical clusters at ``1−|ρ_S| ≤ 0.30``.
* ``k3_alternatives.json``, ``k4_alternatives.json`` — top-N feasible
  K-sets; surfaces stage-3 candidate rows.
* ``config.json`` — reproducibility metadata.

Driven by :file:`slurm/t_sensitivity_historical.slurm` array job
(one task per T value). Aggregation across T is the responsibility of
:mod:`workflows.02_calibration.t_sensitivity_aggregate`.
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
from src.metrics.extended import (  # noqa: E402
    CANDIDATE_METRIC_NAMES,
    CONCEPT_MAP,
    FullRecordRefs,
)
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
DRIVER = "t_sensitivity_historical"

#: Saturation tolerance: a metric value is "saturated at max" if it
#: equals the per-T max within this relative tolerance.
SATURATION_RTOL = 1e-6


def degeneracy_table(chars_df: pd.DataFrame) -> pd.DataFrame:
    """Per-metric degeneracy diagnostics + per-block event counts.

    Computes for every candidate metric the fraction of blocks where:

    * the value is exactly zero (proxy for "no events captured at all"
      on count- and event-derived metrics),
    * the value equals the per-T maximum (saturated — flag for
      bounded-above metrics like ``max_duration`` clustering at 9 mo
      across most 20-yr blocks),
    * the value is NaN (definitional failure, e.g. ``cv_annual_min``
      on a degenerate block).

    Also embeds the SSI-3 and SSI-12 event counts (``n_events``,
    ``n_events_ssi12``) — kept as separate columns alongside the metric
    rows so downstream scripts can read both without joining tables.
    """
    rows: List[Dict[str, object]] = []
    n_blocks = len(chars_df)

    for name in CANDIDATE_METRIC_NAMES:
        if name not in chars_df.columns:
            continue
        v = chars_df[name].astype(float).values
        finite = np.isfinite(v)
        n_finite = int(finite.sum())
        n_nan = int(n_blocks - n_finite)
        if n_finite == 0:
            rows.append({
                "metric": name,
                "concept": CONCEPT_MAP.get(name, "unknown"),
                "n_blocks": n_blocks,
                "n_nan": n_nan,
                "frac_nan": 1.0,
                "frac_zero": float("nan"),
                "frac_saturated_at_max": float("nan"),
                "n_distinct_values": 0,
                "max_value": float("nan"),
            })
            continue
        v_finite = v[finite]
        max_val = float(np.max(v_finite))
        n_zero = int(np.sum(v_finite == 0.0))
        if max_val > 0:
            n_saturated = int(np.sum(np.isclose(v_finite, max_val,
                                                rtol=SATURATION_RTOL)))
        else:
            n_saturated = int(np.sum(v_finite == max_val))
        rows.append({
            "metric": name,
            "concept": CONCEPT_MAP.get(name, "unknown"),
            "n_blocks": n_blocks,
            "n_nan": n_nan,
            "frac_nan": n_nan / n_blocks,
            "frac_zero": n_zero / n_finite,
            "frac_saturated_at_max": n_saturated / n_finite,
            "n_distinct_values": int(np.unique(v_finite).size),
            "max_value": max_val,
        })

    deg_df = pd.DataFrame(rows)

    # Event-count diagnostics: per-block counts from the candidate
    # extractor, surfaced as the *aggregate* mean / median / min / max
    # across blocks.
    diag = {}
    for col in ("n_events", "n_events_ssi12"):
        if col in chars_df.columns:
            v = chars_df[col].astype(float).values
            v = v[np.isfinite(v)]
            if v.size > 0:
                diag[col] = {
                    "mean": float(np.mean(v)),
                    "median": float(np.median(v)),
                    "min": float(np.min(v)),
                    "max": float(np.max(v)),
                    "frac_blocks_zero": float(np.mean(v == 0)),
                }
    deg_df.attrs["event_count_diagnostics"] = diag
    return deg_df


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-years", type=int, required=True,
                   help="Length of each T-year block (one of the sweep grid).")
    p.add_argument("--stride", type=int, default=1,
                   help="Block stride in years (default 1; max overlap).")
    p.add_argument("--cache-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / "data_cache")
    p.add_argument("--correlation-cap", type=float, default=0.6)
    p.add_argument("--cluster-cut", type=float, default=CLUSTER_DISTANCE_CUT)
    args = p.parse_args()

    slug = f"T{args.T_years:02d}"
    out_dir = stage_output_dir(STAGE, DRIVER, slug=slug)
    print(f"[02/t_sensitivity_historical] T={args.T_years} → {out_dir}")

    # --- Data ---
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(args.cache_dir)
    n_hist_years = monthly_2d.shape[0]
    if n_hist_years < args.T_years:
        raise SystemExit(
            f"Historical record has {n_hist_years} years; cannot block at T={args.T_years}"
        )
    n_blocks = (n_hist_years - args.T_years) // args.stride + 1
    print(f"[diag] historical: {n_hist_years} water years, {n_blocks} blocks")

    # --- SSI / Q80 calibration on FULL record (DD-11 lock-in) ---
    ssi3_calc = make_ssi_calculator(timescale=3)
    ssi3_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    ssi12_calc = make_ssi_calculator(timescale=12)
    ssi12_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    q80 = float(np.percentile(monthly_1d, 20.0))
    refs = FullRecordRefs.from_full_record(monthly_1d)
    print(f"[diag] Q80 (fixed full-record): {q80:.2f} cfs")
    print(f"[diag] full-record refs: monthly μ={refs.monthly_mean:.1f} σ={refs.monthly_std:.1f}; "
          f"annual μ={refs.annual_mean_mean:.1f} σ={refs.annual_mean_std:.1f}")

    # --- Per-block candidate matrix (28 legacy + 6 Tier-G hazard-clean) ---
    print("[diag] computing per-block candidate matrix ...")
    chars_df = compute_block_matrix(
        monthly_1d, args.T_years, args.stride, ssi3_calc, ssi12_calc, q80,
        full_record_refs=refs,
    )
    chars_df.to_csv(out_dir / "block_chars_extended.csv", index=False)
    print(f"[diag] wrote block_chars_extended.csv "
          f"({len(chars_df)} blocks × {len(chars_df.columns)} cols)")

    # --- Per-metric spread + screen ---
    spread_df = per_metric_spread(chars_df)
    spread_df.to_csv(out_dir / "per_metric_spread.csv", index=False)
    surviving = spread_df[spread_df["passes_screen"]]["metric"].tolist()
    print(f"[diag] {len(surviving)}/{len(spread_df)} metrics pass screen")

    # --- Degeneracy table ---
    deg_df = degeneracy_table(chars_df)
    deg_df.to_csv(out_dir / "degeneracy.csv", index=False)
    event_diag = deg_df.attrs.get("event_count_diagnostics", {})
    (out_dir / "event_count_diagnostics.json").write_text(
        json.dumps(event_diag, indent=2)
    )
    if event_diag:
        for k, v in event_diag.items():
            print(f"[diag] {k}: median={v['median']:.1f}, "
                  f"frac_blocks_zero={v['frac_blocks_zero']:.3f}")

    # --- Correlations on surviving metrics ---
    if len(surviving) >= 2:
        pearson, spearman = correlation_matrices(chars_df, surviving)
        pearson.to_csv(out_dir / "pearson_corr.csv")
        spearman.to_csv(out_dir / "spearman_corr.csv")
    else:
        print("[warn] fewer than 2 surviving metrics; skipping correlations")
        spearman = pd.DataFrame()
        pearson = pd.DataFrame()

    # --- Clusters ---
    if len(surviving) >= 2:
        spread_lookup = spread_df.set_index("metric")["spread_score"]
        cluster_df = cluster_metrics(spearman, spread_lookup, args.cluster_cut)
        cluster_df.to_json(out_dir / "clusters.json", orient="records", indent=2)
        n_clusters = int(cluster_df["cluster_id"].nunique())
        n_concepts = int(cluster_df["concept"].nunique())
        print(f"[diag] {n_clusters} clusters / {n_concepts} concepts")
    else:
        cluster_df = pd.DataFrame(columns=["metric", "cluster_id", "concept",
                                           "is_representative"])
        n_clusters = n_concepts = 0

    # --- K-set enumeration (K=3 and K=4) ---
    k3_alts, k3_rung = ([], "none")
    k4_alts, k4_rung = ([], "none")
    if len(surviving) >= 3:
        k3_alts, k3_rung = relax_until_nonempty(
            surviving, spearman, cluster_df, spread_df, K=3,
            correlation_cap=args.correlation_cap,
        )
    if len(surviving) >= 4:
        k4_alts, k4_rung = relax_until_nonempty(
            surviving, spearman, cluster_df, spread_df, K=4,
            correlation_cap=args.correlation_cap,
        )

    (out_dir / "k3_alternatives.json").write_text(
        json.dumps(
            {"constraint_rung": k3_rung,
             "alternatives": k3_alts[:N_ALTERNATIVES]},
            indent=2,
        )
    )
    (out_dir / "k4_alternatives.json").write_text(
        json.dumps(
            {"constraint_rung": k4_rung,
             "alternatives": k4_alts[:N_ALTERNATIVES]},
            indent=2,
        )
    )
    if k3_alts:
        print(f"[rec] K=3 ({k3_rung}): {k3_alts[0]['metrics']} "
              f"min_spread={k3_alts[0]['min_spread']:.3f}")
    if k4_alts:
        print(f"[rec] K=4 ({k4_rung}): {k4_alts[0]['metrics']} "
              f"min_spread={k4_alts[0]['min_spread']:.3f}")

    # --- Config ---
    cfg = {
        "script": "workflows/02_calibration/t_sensitivity_historical.py",
        "T_years": int(args.T_years),
        "stride": int(args.stride),
        "n_hist_years": int(n_hist_years),
        "n_blocks": int(n_blocks),
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
    print(f"[diag] done. wrote {out_dir / 'config.json'}")


if __name__ == "__main__":
    main()
