"""Stage 2 — compare baseline-Kirsch metric distributions to historical T-blocks.

Joins the per-T outputs of
:mod:`workflows.02_calibration.t_sensitivity_historical` (historical
block-metric matrices) and
:mod:`workflows.03_kirsch_library.build_library_extended` (Kirsch
ensemble metric matrices) and produces fidelity diagnostics across the
full T-grid.

Outputs (under ``outputs/02_calibration/t_sensitivity_kirsch_compare/``):

* ``ks_vs_historical.csv`` — long-form (metric, T) → KS statistic and
  p-value comparing baseline-Kirsch ensemble distribution to the
  historical T-block distribution. Threshold ``ks < 0.10`` is the
  acceptable-fidelity gate.
* ``frobenius_corrshift.csv`` — per T, ``‖ρ_Kirsch − ρ_hist‖_F`` with
  the metrics restricted to the intersection of surviving subsets at
  that T. Identifies the T below which the Kirsch generator fails to
  preserve the inter-metric correlation skeleton.
* ``coverage.csv`` — per (T, K-set) coverage diagnostics
  (L2-star discrepancy, NN_CV) of the Kirsch ensemble in each candidate
  K-set's sub-space. Provides the coverage *baseline* against which
  MOEA-FIND must improve in Stage 5 figures.
* ``surviving_T_grid_kirsch.json`` — list of T values where the
  Kirsch generator passes the KS gate (max metric KS ≤ 0.10) AND the
  correlation-preservation gate. Stage 3 reads this list to know which
  Ts to keep as decision-matrix columns.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.discovery.analysis import coverage_metrics  # noqa: E402
from src.metrics.extended import CANDIDATE_METRIC_NAMES  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "02_calibration"
HIST_DRIVER = "t_sensitivity_historical"
KIRSCH_STAGE = "03_kirsch_library"
KIRSCH_DRIVER = "build_library_extended"
DRIVER = "t_sensitivity_kirsch_compare"

#: KS-statistic threshold; per-metric KS above this is a fidelity flag.
KS_GATE = 0.10

#: Frobenius-norm gate on the elementwise correlation difference. Above
#: this, the Kirsch correlation skeleton has diverged enough from
#: historical that downstream metric-pair selections may not transfer.
FROBENIUS_GATE = 1.5


def _hist_dir(T: int) -> Path:
    return stage_output_dir(STAGE, HIST_DRIVER, slug=f"T{T:02d}", create=False)


def _kirsch_slug(n_traces: int, T: int, seed: int) -> str:
    return f"n{n_traces}_t{T}_ssi3-12_s{seed}"


def _kirsch_dir(n_traces: int, T: int, seed: int) -> Path:
    return stage_output_dir(
        KIRSCH_STAGE, KIRSCH_DRIVER,
        slug=_kirsch_slug(n_traces, T, seed),
        create=False,
    )


def _load_hist_block_chars(T: int) -> Optional[pd.DataFrame]:
    p = _hist_dir(T) / "block_chars_extended.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def _load_hist_spread(T: int) -> Optional[pd.DataFrame]:
    p = _hist_dir(T) / "per_metric_spread.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def _load_kirsch_chars(n_traces: int, T: int, seed: int) -> Optional[pd.DataFrame]:
    p = _kirsch_dir(n_traces, T, seed) / "characteristics_extended.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    metric_names = [str(n) for n in z["metric_names"]]
    return pd.DataFrame(z["values"], columns=metric_names)


def _load_kset_alternatives(T: int, K: int) -> List[List[str]]:
    """Top-N alternative metric tuples at this T from the Stage-1 driver."""
    p = _hist_dir(T) / f"k{K}_alternatives.json"
    if not p.exists():
        return []
    raw = json.loads(p.read_text())
    return [list(a["metrics"]) for a in raw.get("alternatives", [])]


def ks_long_form(
    hist_T: Dict[int, pd.DataFrame],
    kirsch_T: Dict[int, pd.DataFrame],
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for T in sorted(set(hist_T) & set(kirsch_T)):
        hist_df = hist_T[T]
        kir_df = kirsch_T[T]
        for m in CANDIDATE_METRIC_NAMES:
            if m not in hist_df.columns or m not in kir_df.columns:
                continue
            h = hist_df[m].astype(float).values
            k = kir_df[m].astype(float).values
            h = h[np.isfinite(h)]
            k = k[np.isfinite(k)]
            if h.size < 2 or k.size < 2:
                continue
            stat, pval = ks_2samp(h, k, alternative="two-sided", mode="auto")
            rows.append({
                "metric": m,
                "T_years": T,
                "n_hist": int(h.size),
                "n_kirsch": int(k.size),
                "ks_statistic": float(stat),
                "ks_pvalue": float(pval),
                "passes_ks_gate": bool(stat <= KS_GATE),
            })
    return pd.DataFrame(rows)


def frobenius_corr_shift(
    hist_T: Dict[int, pd.DataFrame],
    kirsch_T: Dict[int, pd.DataFrame],
    hist_spread_T: Dict[int, pd.DataFrame],
) -> pd.DataFrame:
    """Frobenius norm of the Spearman-correlation difference per T.

    Restricted to the metrics that pass spread/skew screening on the
    historical record at that T (so we don't compare correlations on
    metrics that were already excluded as MOEA-unviable axes).
    """
    rows: List[Dict[str, object]] = []
    for T in sorted(set(hist_T) & set(kirsch_T) & set(hist_spread_T)):
        hist_df = hist_T[T]
        kir_df = kirsch_T[T]
        spread = hist_spread_T[T]
        surviving = spread[spread["passes_screen"]]["metric"].tolist()
        cols = [m for m in surviving
                if m in hist_df.columns and m in kir_df.columns]
        if len(cols) < 3:
            continue
        rho_hist = hist_df[cols].astype(float).corr(method="spearman").values
        rho_kir = kir_df[cols].astype(float).corr(method="spearman").values
        diff = rho_kir - rho_hist
        frob = float(np.linalg.norm(diff, ord="fro"))
        rows.append({
            "T_years": T,
            "n_metrics_compared": len(cols),
            "frobenius_corr_shift": frob,
            "passes_correlation_gate": bool(frob <= FROBENIUS_GATE),
        })
    return pd.DataFrame(rows)


def coverage_long(
    kirsch_T: Dict[int, pd.DataFrame],
    hist_T: Dict[int, pd.DataFrame],
) -> pd.DataFrame:
    """L2-star and NN-CV of the Kirsch ensemble in each candidate K-set's space.

    Bounds for normalization come from the union of historical-block
    and Kirsch values per metric, so that historical and synthetic
    extents share the same coverage frame.
    """
    rows: List[Dict[str, object]] = []
    for T in sorted(set(hist_T) & set(kirsch_T)):
        hist_df = hist_T[T]
        kir_df = kirsch_T[T]
        for K in (3, 4):
            for rank, kset in enumerate(_load_kset_alternatives(T, K)[:5], 1):
                cols = [m for m in kset
                        if m in hist_df.columns and m in kir_df.columns]
                if len(cols) != K:
                    continue
                pts = kir_df[cols].astype(float).values
                pts = pts[np.isfinite(pts).all(axis=1)]
                if pts.shape[0] < 2:
                    continue
                hist_pts = hist_df[cols].astype(float).values
                hist_pts = hist_pts[np.isfinite(hist_pts).all(axis=1)]
                combined = np.vstack([pts, hist_pts]) if hist_pts.size else pts
                lb = combined.min(axis=0)
                ub = combined.max(axis=0)
                ub = np.where(ub > lb, ub, lb + 1e-9)
                cov = coverage_metrics(pts, lb, ub)
                rows.append({
                    "T_years": T,
                    "K": K,
                    "rank": rank,
                    "metrics": "|".join(cols),
                    "n_traces": int(pts.shape[0]),
                    "L2_star_discrepancy": cov["L2_star_discrepancy"],
                    "nn_mean": cov.get("nn_mean", float("nan")),
                    "nn_cv": cov.get("nn_cv", float("nan")),
                })
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-grid", type=int, nargs="+",
                   default=[5, 10, 20, 30])
    p.add_argument("--n-traces", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER)
    print(f"[02/{DRIVER}] T-grid={args.T_grid} → {out_dir}")

    hist_T: Dict[int, pd.DataFrame] = {}
    kirsch_T: Dict[int, pd.DataFrame] = {}
    hist_spread_T: Dict[int, pd.DataFrame] = {}
    for T in args.T_grid:
        h = _load_hist_block_chars(T)
        k = _load_kirsch_chars(args.n_traces, T, args.seed)
        s = _load_hist_spread(T)
        if h is None:
            print(f"[warn] T={T}: missing historical block_chars; skipping")
            continue
        if k is None:
            print(f"[warn] T={T}: missing Kirsch characteristics_extended.npz; "
                  f"skipping")
            continue
        if s is None:
            print(f"[warn] T={T}: missing per_metric_spread; skipping")
            continue
        hist_T[T] = h
        kirsch_T[T] = k
        hist_spread_T[T] = s
        print(f"[diag] T={T}: hist={len(h)} blocks, kirsch={len(k)} traces")

    if not hist_T:
        raise SystemExit("No (T) bundles available — run Stage 1 + "
                         "build_library_extended first.")

    ks_df = ks_long_form(hist_T, kirsch_T)
    ks_df.to_csv(out_dir / "ks_vs_historical.csv", index=False)
    print(f"[diag] wrote ks_vs_historical.csv ({len(ks_df)} rows)")

    frob_df = frobenius_corr_shift(hist_T, kirsch_T, hist_spread_T)
    frob_df.to_csv(out_dir / "frobenius_corrshift.csv", index=False)
    print(f"[diag] wrote frobenius_corrshift.csv ({len(frob_df)} rows)")

    cov_df = coverage_long(kirsch_T, hist_T)
    cov_df.to_csv(out_dir / "coverage.csv", index=False)
    print(f"[diag] wrote coverage.csv ({len(cov_df)} rows)")

    # --- Apply Kirsch-fidelity gate per T ---
    survivors: List[int] = []
    by_T = {}
    for T in sorted(hist_T):
        ks_max = ks_df[ks_df["T_years"] == T]["ks_statistic"].max() \
            if not ks_df.empty else float("nan")
        frob = frob_df[frob_df["T_years"] == T]["frobenius_corr_shift"].max() \
            if not frob_df.empty else float("nan")
        ks_pass = bool(np.isfinite(ks_max) and ks_max <= KS_GATE)
        frob_pass = bool(np.isfinite(frob) and frob <= FROBENIUS_GATE)
        by_T[T] = {"ks_max": float(ks_max),
                   "frobenius": float(frob),
                   "ks_pass": ks_pass,
                   "frobenius_pass": frob_pass}
        if ks_pass and frob_pass:
            survivors.append(T)
        else:
            print(f"[gate] T={T}: ks_max={ks_max:.3f} (gate {KS_GATE}) "
                  f"frob={frob:.3f} (gate {FROBENIUS_GATE}) → DROP")

    (out_dir / "surviving_T_grid_kirsch.json").write_text(
        json.dumps({
            "surviving_T_grid": survivors,
            "input_T_grid": list(args.T_grid),
            "per_T": by_T,
            "ks_gate": KS_GATE,
            "frobenius_gate": FROBENIUS_GATE,
        }, indent=2)
    )
    print(f"[gate] surviving T (Kirsch fidelity): {survivors}")
    print("[diag] done.")


if __name__ == "__main__":
    main()
