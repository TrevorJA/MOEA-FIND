"""Stage 1 — cross-T aggregation of historical-block sensitivity diagnostics.

Reads the per-T outputs of
:mod:`workflows.02_calibration.t_sensitivity_historical` (one slug per
T value) and produces:

* ``t_sensitivity_summary.csv`` — long-form ``(metric, T)`` →
  spread/range/skew/degeneracy. Easy to plot.
* ``cluster_stability.csv`` — Adjusted Rand Index between cluster
  memberships at each pair of surviving T values, computed over the
  metrics that pass screen at *both* Ts. Quantifies whether the
  redundancy structure is preserved across T.
* ``kset_stability.csv`` — for every (T, K) pair, the Jaccard score of
  the top-rank K-set against the T=20 reference set, plus the metric
  list. Identifies *T-portable* K-sets — those that emerge at multiple
  Ts and therefore survive the joint K×T decision.
* ``surviving_T_grid.json`` — list of T values that satisfy the
  Stage-1 viability gate (≤25% degenerate metrics on ≤2 candidates;
  non-empty K=3 strict-rung). Stage 2 reads this list to know which Ts
  to validate Kirsch fidelity for.

Run after the per-T array job completes.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "02_calibration"
PER_T_DRIVER = "t_sensitivity_historical"
DRIVER = "t_sensitivity_aggregate"

#: Reference T against which K-set stability is benchmarked.
REFERENCE_T = 20

#: A metric is considered hard-degenerate at a given T if its zero or
#: NaN fraction exceeds this threshold. (Saturation-at-max is excluded:
#: it does not violate the spread-screen IQR>0 rule on its own and is
#: instead surfaced as a per-metric quality flag in the figures.)
DEGENERACY_THRESHOLD = 0.25

#: Stage-1 viability gate: drop a T only if MORE than this fraction of
#: the 28 candidate metrics are hard-degenerate (zero/NaN > threshold)
#: AND the strict clusters-and-concepts K=3 rung is empty. Lenient by
#: design — we want to keep T values where a sensible K=3 axis triple
#: still emerges even after excluding the SSI-12 family at small T.
MAX_HARD_DEGENERATE_METRIC_FRAC = 0.5


def _per_t_dir(T: int) -> Path:
    return stage_output_dir(STAGE, PER_T_DRIVER, slug=f"T{T:02d}", create=False)


def _load_per_t_artifacts(T: int) -> Optional[Dict[str, object]]:
    """Load the artifact bundle for one T value, or None if missing."""
    d = _per_t_dir(T)
    if not d.exists():
        return None
    try:
        spread = pd.read_csv(d / "per_metric_spread.csv")
        deg = pd.read_csv(d / "degeneracy.csv")
        cfg = json.loads((d / "config.json").read_text())
        clusters_path = d / "clusters.json"
        clusters = (
            pd.read_json(clusters_path) if clusters_path.exists() else pd.DataFrame()
        )
        k3 = json.loads((d / "k3_alternatives.json").read_text())
        k4 = json.loads((d / "k4_alternatives.json").read_text())
    except FileNotFoundError as e:
        print(f"[warn] T={T}: missing artifact ({e}); skipping")
        return None
    return {
        "T": T,
        "spread": spread,
        "deg": deg,
        "cfg": cfg,
        "clusters": clusters,
        "k3": k3,
        "k4": k4,
    }


def long_form_summary(bundles: List[Dict[str, object]]) -> pd.DataFrame:
    """Stack per-T spread + degeneracy into a single long-form table."""
    rows: List[Dict[str, object]] = []
    for b in bundles:
        T = b["T"]
        spread = b["spread"].set_index("metric")
        deg = b["deg"].set_index("metric")
        for metric in spread.index.union(deg.index):
            row: Dict[str, object] = {"metric": metric, "T_years": T}
            if metric in spread.index:
                s = spread.loc[metric]
                row.update({
                    "concept": s.get("concept"),
                    "median": s.get("median"),
                    "iqr": s.get("iqr"),
                    "range": s.get("range"),
                    "spread_score": s.get("spread_score"),
                    "skew": s.get("skew"),
                    "passes_screen": bool(s.get("passes_screen", False)),
                })
            if metric in deg.index:
                d = deg.loc[metric]
                row.update({
                    "frac_zero": d.get("frac_zero"),
                    "frac_nan": d.get("frac_nan"),
                    "frac_saturated_at_max": d.get("frac_saturated_at_max"),
                    "n_distinct_values": d.get("n_distinct_values"),
                })
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["metric", "T_years"]).reset_index(drop=True)


def cluster_stability(bundles: List[Dict[str, object]]) -> pd.DataFrame:
    """Pairwise Adjusted Rand Index between cluster memberships."""
    rows: List[Dict[str, object]] = []
    for i, a in enumerate(bundles):
        for b in bundles[i:]:
            ca, cb = a["clusters"], b["clusters"]
            if ca.empty or cb.empty:
                continue
            ca_map = ca.set_index("metric")["cluster_id"].to_dict()
            cb_map = cb.set_index("metric")["cluster_id"].to_dict()
            shared = sorted(set(ca_map) & set(cb_map))
            if len(shared) < 2:
                continue
            la = [ca_map[m] for m in shared]
            lb = [cb_map[m] for m in shared]
            ari = float(adjusted_rand_score(la, lb))
            rows.append({
                "T_a": a["T"],
                "T_b": b["T"],
                "n_shared_metrics": len(shared),
                "adjusted_rand_index": ari,
            })
    return pd.DataFrame(rows)


def _top_kset(alt_bundle: Dict[str, object]) -> Optional[List[str]]:
    alts = alt_bundle.get("alternatives", [])
    if not alts:
        return None
    return list(alts[0]["metrics"])


def _jaccard(a: Optional[List[str]], b: Optional[List[str]]) -> float:
    if not a or not b:
        return float("nan")
    sa, sb = set(a), set(b)
    union = sa | sb
    if not union:
        return float("nan")
    return len(sa & sb) / len(union)


def kset_stability(bundles: List[Dict[str, object]]) -> pd.DataFrame:
    """Top-K set Jaccard vs the reference T=20 K-set."""
    by_T = {b["T"]: b for b in bundles}
    ref = by_T.get(REFERENCE_T)
    if ref is None:
        print(f"[warn] no reference T={REFERENCE_T} bundle; jaccard against "
              f"smallest T instead.")
        ref = bundles[0] if bundles else None
        ref_T = ref["T"] if ref else None
    else:
        ref_T = REFERENCE_T

    rows: List[Dict[str, object]] = []
    for K_label, key in (("K=3", "k3"), ("K=4", "k4")):
        ref_set = _top_kset(ref[key]) if ref else None
        for b in bundles:
            top = _top_kset(b[key])
            rows.append({
                "K": K_label,
                "T_years": b["T"],
                "constraint_rung": b[key].get("constraint_rung"),
                "top_set": top,
                "reference_T": ref_T,
                "reference_set": ref_set,
                "jaccard_vs_reference": _jaccard(top, ref_set),
            })
    return pd.DataFrame(rows)


def surviving_T_grid(bundles: List[Dict[str, object]]) -> List[int]:
    """Apply the Stage-1 viability gate to identify usable T values.

    Lenient by design: we keep a T if its strict clusters-and-concepts
    K=3 rung is non-empty AND fewer than half of the 28 candidate
    metrics are hard-degenerate (zero or NaN > 25% of blocks). At small
    T the SSI-12 family will hard-degenerate; that does not kill the
    T, it just shrinks the viable candidate pool — which is fine
    because we still find hundreds of feasible Tier-A/C/D/E/F K=3 sets.
    """
    survivors: List[int] = []
    for b in bundles:
        deg = b["deg"]
        if deg.empty:
            continue
        hard_mask = (
            (deg["frac_zero"].fillna(0) > DEGENERACY_THRESHOLD)
            | (deg["frac_nan"].fillna(0) > DEGENERACY_THRESHOLD)
        )
        n_total = len(deg)
        hard_frac = float(hard_mask.sum()) / max(n_total, 1)
        sat_frac = float(
            (deg["frac_saturated_at_max"].fillna(0)
             > DEGENERACY_THRESHOLD).sum()
        ) / max(n_total, 1)
        k3_alts = b["k3"].get("alternatives", [])
        k3_rung = b["k3"].get("constraint_rung")
        ok_hard = hard_frac <= MAX_HARD_DEGENERATE_METRIC_FRAC
        ok_strict_k3 = (k3_rung == "distinct_clusters_and_concepts"
                        and len(k3_alts) > 0)
        verdict = "PASS" if (ok_hard and ok_strict_k3) else "DROP"
        print(f"[gate] T={b['T']}: hard_degenerate_frac={hard_frac:.3f} "
              f"(threshold {MAX_HARD_DEGENERATE_METRIC_FRAC:.2f}), "
              f"saturation_frac={sat_frac:.3f} (informational), "
              f"k3_rung={k3_rung} k3_alts={len(k3_alts)} → {verdict}")
        if ok_hard and ok_strict_k3:
            survivors.append(int(b["T"]))
    return sorted(set(survivors))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-grid", type=int, nargs="+",
                   default=[5, 10, 20, 30],
                   help="T values to aggregate (default: 5 10 20 30).")
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER)
    print(f"[02/{DRIVER}] T-grid={args.T_grid} → {out_dir}")

    bundles: List[Dict[str, object]] = []
    for T in args.T_grid:
        b = _load_per_t_artifacts(T)
        if b is None:
            print(f"[warn] T={T} artifacts missing; run per-T driver first.")
            continue
        bundles.append(b)

    if not bundles:
        raise SystemExit("No per-T artifacts found. Run "
                         "t_sensitivity_historical.py first.")

    summary_df = long_form_summary(bundles)
    summary_df.to_csv(out_dir / "t_sensitivity_summary.csv", index=False)
    print(f"[diag] wrote t_sensitivity_summary.csv "
          f"({len(summary_df)} rows)")

    ari_df = cluster_stability(bundles)
    ari_df.to_csv(out_dir / "cluster_stability.csv", index=False)
    print(f"[diag] wrote cluster_stability.csv "
          f"({len(ari_df)} pairwise rows)")

    kset_df = kset_stability(bundles)
    kset_df.to_csv(out_dir / "kset_stability.csv", index=False)
    print(f"[diag] wrote kset_stability.csv ({len(kset_df)} rows)")

    survivors = surviving_T_grid(bundles)
    (out_dir / "surviving_T_grid.json").write_text(
        json.dumps({"surviving_T_grid": survivors,
                    "input_T_grid": list(args.T_grid),
                    "criteria": {
                        "max_hard_degenerate_metric_frac":
                            MAX_HARD_DEGENERATE_METRIC_FRAC,
                        "degeneracy_threshold": DEGENERACY_THRESHOLD,
                        "saturation_excluded_from_gate": True,
                        "require_strict_k3_rung": True,
                    }}, indent=2)
    )
    print(f"[gate] surviving T grid: {survivors}")
    print("[diag] done.")


if __name__ == "__main__":
    main()
