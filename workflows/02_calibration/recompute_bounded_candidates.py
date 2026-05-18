"""DD-15c — recompute the 48 bounded candidate metrics on existing data.

Phase 2 of the bounded T=1 reformulation. Loads existing artifacts and
emits parquet outputs that drive the diagnostic-driven (mapping × K=4
windows) selection in :mod:`select_bounded_kset` and the upgraded
calibration figures.

Produces:

* ``historical.parquet`` — per-water-year metric matrix on the historical
  record. Each row is a single T=1 water year; columns include the 48
  bounded candidates plus a ``year`` index.
* ``{archive_slug}.parquet`` — per-Pareto-member metric matrix. One row
  per Pareto solution; columns include the 48 bounded candidates, the
  legacy ``short_block_drb`` objective values (for backward comparison),
  and per-trace flow diagnostics (``max_monthly_flow``,
  ``frac_monthly_above_10k_cfs``, ``min_monthly_flow``) for flood-corner
  scoring in Phase 3.
* ``family_refs.npz`` — frozen :class:`BoundedFamilyRefs` for the
  figure-rendering stage (mu/sigma/sorted historical sample per window).
* ``config.json`` — reproducibility metadata.

The four "before" Pareto archives processed by default:

1. ``residual_T1_nfe200000_s42_constrained_cmdv_uniform_stad`` (residual+AD baseline)
2. ``residual_T1_nfe200000_s42_constrained_cmdv_uniform_sfxiter2_stad`` (residual+AD-iter2)
3. ``residual_T1_nfe200000_s42_constrained_cmdv_uniform_stks`` (residual+KS)
4. ``index_T1_nfe200000_s42_constrained_cmdv_uniform_stad`` (index+AD)

Light Python compute (~1–2 minutes for ~40k Pareto members × 48 candidates).
Submitted via SLURM per project convention (see ``slurm/recompute_bounded_candidates.slurm``).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.metrics.extended import (  # noqa: E402
    BoundedFamilyRefs,
    FullRecordRefs,
)
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.metrics.short_block import (  # noqa: E402
    CANDIDATE_BOUNDED_METRIC_NAMES,
    CANDIDATE_BOUNDED_CONCEPT_MAP,
    SHORT_BLOCK_METRIC_NAMES,
    WINDOW_SPECS,
    compute_candidate_bounded_metrics,
    compute_short_block_metrics,
)


_DEFAULT_ARCHIVES: Tuple[str, ...] = (
    "residual_T1_nfe200000_s42_constrained_cmdv_uniform_stad",
    "residual_T1_nfe200000_s42_constrained_cmdv_uniform_sfxiter2_stad",
    "residual_T1_nfe200000_s42_constrained_cmdv_uniform_stks",
    "index_T1_nfe200000_s42_constrained_cmdv_uniform_stad",
)

_FLOOD_THRESHOLD_CFS = 10_000.0


def _archive_results_path(slug: str) -> Path:
    return (
        PROJECT_ROOT
        / "outputs"
        / "04_moea_find_single_site"
        / "run_moea_find"
        / slug
        / "results.json"
    )


def _save_family_refs_npz(refs: BoundedFamilyRefs, path: Path) -> None:
    """Persist the family refs as a flat ``.npz`` (one array per field)."""
    arrays = {}
    for window, payload in refs.per_window.items():
        arrays[f"{window}__mu"] = np.asarray(payload["mu"], dtype=float)
        arrays[f"{window}__sigma"] = np.asarray(payload["sigma"], dtype=float)
        arrays[f"{window}__sorted"] = np.asarray(payload["sorted"], dtype=float)
    np.savez(path, **arrays)


def _historical_per_year_table(
    monthly_2d: np.ndarray,
    family_refs: BoundedFamilyRefs,
    full_refs: FullRecordRefs,
) -> pd.DataFrame:
    """Compute candidate + legacy metrics on each historical water year."""
    n_years = monthly_2d.shape[0]
    rows: List[Dict[str, float]] = []
    for y in range(n_years):
        yr = monthly_2d[y]  # shape (12,)
        bounded = compute_candidate_bounded_metrics(yr, family_refs)
        legacy = compute_short_block_metrics(yr, ssi3_full_window_series=None, refs=full_refs)
        row: Dict[str, float] = {"year_index": int(y)}
        row["max_monthly_flow"] = float(yr.max())
        row["min_monthly_flow"] = float(yr.min())
        row["frac_monthly_above_10k"] = float((yr > _FLOOD_THRESHOLD_CFS).mean())
        for k in CANDIDATE_BOUNDED_METRIC_NAMES:
            row[k] = bounded[k]
        for k in SHORT_BLOCK_METRIC_NAMES:
            row[f"legacy__{k}"] = legacy.get(k, 0.0)
        rows.append(row)
    return pd.DataFrame(rows)


def _archive_pareto_table(
    archive_slug: str,
    family_refs: BoundedFamilyRefs,
    full_refs: FullRecordRefs,
) -> pd.DataFrame:
    """Compute candidate + legacy metrics + flow diagnostics on a Pareto archive."""
    results_path = _archive_results_path(archive_slug)
    if not results_path.exists():
        raise FileNotFoundError(
            f"Archive results.json not found for {archive_slug!r}: {results_path}"
        )
    payload = json.loads(results_path.read_text())
    traces = payload.get("pareto_traces_1d", [])
    drought_metrics = payload.get("drought_metrics", [])
    objective_keys = payload.get("objective_keys", [])
    n_pareto = len(traces)
    if n_pareto == 0:
        raise ValueError(f"{archive_slug}: no Pareto traces in results.json")

    rows: List[Dict[str, float]] = []
    for i, trace in enumerate(traces):
        flows = np.asarray(trace, dtype=float).ravel()
        if flows.size != 12:
            raise ValueError(
                f"{archive_slug} Pareto member {i}: expected 12 monthly flows "
                f"(T=1), got {flows.size}"
            )
        bounded = compute_candidate_bounded_metrics(flows, family_refs)
        legacy = compute_short_block_metrics(flows, ssi3_full_window_series=None, refs=full_refs)
        row: Dict[str, float] = {"pareto_index": int(i)}
        row["max_monthly_flow"] = float(flows.max())
        row["min_monthly_flow"] = float(flows.min())
        row["frac_monthly_above_10k"] = float((flows > _FLOOD_THRESHOLD_CFS).mean())
        for k in CANDIDATE_BOUNDED_METRIC_NAMES:
            row[k] = bounded[k]
        for k in SHORT_BLOCK_METRIC_NAMES:
            row[f"legacy__{k}"] = legacy.get(k, 0.0)
        # Persist the live optimisation objective values so plots can compare
        # bounded candidates against the actual optimised metric values.
        if drought_metrics and i < len(drought_metrics):
            obj_row = drought_metrics[i]
            for j, name in enumerate(objective_keys):
                if j < len(obj_row):
                    row[f"objective__{name}"] = float(obj_row[j])
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archives",
        nargs="*",
        default=list(_DEFAULT_ARCHIVES),
        help="Archive slugs to recompute (default: 4 DD-15c 'before' archives).",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=PROJECT_ROOT / "outputs" / "data_cache",
        help="USGS cache directory.",
    )
    args = parser.parse_args()

    out_dir = stage_output_dir("02_calibration", "recompute_bounded_candidates")
    print(f"[recompute_bounded_candidates] output dir = {out_dir}")

    t0 = time.perf_counter()

    print("[recompute_bounded_candidates] loading historical record ...")
    monthly_2d, monthly_1d = prepare_data(args.cache_dir)
    n_full_years = monthly_2d.shape[0]
    print(f"  historical: {n_full_years} water years, "
          f"{monthly_1d.size} months total")

    print("[recompute_bounded_candidates] building BoundedFamilyRefs ...")
    family_refs = BoundedFamilyRefs.from_full_record(monthly_1d)
    full_refs = FullRecordRefs.from_full_record(monthly_1d)

    print(f"[recompute_bounded_candidates] {len(WINDOW_SPECS)} windows × 2 "
          f"mappings = {len(CANDIDATE_BOUNDED_METRIC_NAMES)} candidate metrics")

    print("[recompute_bounded_candidates] computing historical per-year table ...")
    hist_df = _historical_per_year_table(monthly_2d, family_refs, full_refs)
    hist_path = out_dir / "historical.parquet"
    hist_df.to_parquet(hist_path, index=False)
    print(f"  wrote {hist_path} ({len(hist_df)} rows × {len(hist_df.columns)} cols)")

    refs_path = out_dir / "family_refs.npz"
    _save_family_refs_npz(family_refs, refs_path)
    print(f"  wrote {refs_path}")

    archive_summary: List[Dict[str, object]] = []
    for slug in args.archives:
        print(f"[recompute_bounded_candidates] archive: {slug}")
        try:
            df = _archive_pareto_table(slug, family_refs, full_refs)
        except (FileNotFoundError, ValueError) as exc:
            print(f"  SKIP: {exc}")
            archive_summary.append({"slug": slug, "n_pareto": 0, "skipped": True, "reason": str(exc)})
            continue
        path = out_dir / f"{slug}.parquet"
        df.to_parquet(path, index=False)
        # Summary stats for sanity checks
        cand_max = float(df[list(CANDIDATE_BOUNDED_METRIC_NAMES)].max().max())
        cand_min = float(df[list(CANDIDATE_BOUNDED_METRIC_NAMES)].min().min())
        flood_frac = float((df["max_monthly_flow"] > _FLOOD_THRESHOLD_CFS).mean())
        archive_summary.append({
            "slug": slug,
            "n_pareto": int(len(df)),
            "skipped": False,
            "candidate_d_min": cand_min,
            "candidate_d_max": cand_max,
            "frac_pareto_with_max_monthly_above_10k": flood_frac,
        })
        print(f"  wrote {path} ({len(df)} rows × {len(df.columns)} cols); "
              f"D ∈ [{cand_min:.4f}, {cand_max:.4f}]; "
              f"frac flood-corner Pareto = {flood_frac:.3%}")

    config = {
        "stage": "02_calibration",
        "driver": "recompute_bounded_candidates",
        "n_historical_years": int(n_full_years),
        "n_windows": len(WINDOW_SPECS),
        "n_candidate_metrics": len(CANDIDATE_BOUNDED_METRIC_NAMES),
        "candidate_metric_names": list(CANDIDATE_BOUNDED_METRIC_NAMES),
        "candidate_concept_map": CANDIDATE_BOUNDED_CONCEPT_MAP,
        "flood_threshold_cfs": _FLOOD_THRESHOLD_CFS,
        "archives": archive_summary,
        "elapsed_s": float(time.perf_counter() - t0),
    }
    config_path = out_dir / "config.json"
    config_path.write_text(json.dumps(config, indent=2, default=str))
    print(f"  wrote {config_path}")
    print(f"[recompute_bounded_candidates] DONE in "
          f"{config['elapsed_s']:.1f}s")


if __name__ == "__main__":
    main()
