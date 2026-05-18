"""Recharacterize an already-generated Kirsch library with the current
``compute_all_candidates`` definitions, without regenerating traces.

Use this whenever the metric library
(:mod:`src.metrics.extended`) gains new metrics and the existing
``library.npy`` traces are still valid. The driver reads ``library.npy``
from
``outputs/03_kirsch_library/build_library_extended/<slug>/`` and writes
a refreshed ``characteristics_extended.npz`` next to it. Generation
parameters (n_traces, T, seed) are inferred from ``config.json``.

Cheaper than re-running :mod:`build_library_extended` when the heavy
work is the Kirsch generation itself, not the metric computation.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.metrics.extended import (  # noqa: E402
    CANDIDATE_METRIC_NAMES,
    FullRecordRefs,
    compute_all_candidates,
)
from src.metrics.objectives import flows_to_series, make_ssi_calculator  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "03_kirsch_library"
DRIVER = "build_library_extended"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-years", type=int, required=True)
    p.add_argument("--n-traces", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    slug = f"n{args.n_traces}_t{args.T_years}_ssi3-12_s{args.seed}"
    out = stage_output_dir(STAGE, DRIVER, slug=slug, create=False)
    if not out.exists():
        raise SystemExit(f"library directory missing: {out}")
    library_path = out / "library.npy"
    if not library_path.exists():
        raise SystemExit(f"library.npy missing at {library_path}")

    print(f"[recharacterize] T={args.T_years}, slug={slug}")
    library = np.load(library_path)
    print(f"[recharacterize] loaded library shape={library.shape}")

    # --- Full-record references (must match the historical pipeline) ---
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d_full, monthly_1d_full = prepare_data(cache_dir)

    ssi3_calc = make_ssi_calculator(timescale=3)
    ssi3_calc.fit(flows_to_series(monthly_1d_full, start_date="1950-10-01"))
    ssi12_calc = make_ssi_calculator(timescale=12)
    ssi12_calc.fit(flows_to_series(monthly_1d_full, start_date="1950-10-01"))
    q80 = float(np.percentile(monthly_1d_full, 20.0))
    refs = FullRecordRefs.from_full_record(monthly_1d_full)
    print(f"[recharacterize] Q80 (full record): {q80:.2f} cfs")

    # --- Recompute metrics on every saved trace ---
    print(f"[recharacterize] characterising {library.shape[0]} traces ...")
    t0 = time.perf_counter()
    rows: List[Dict[str, float]] = []
    n_yrs = args.T_years
    for rid in range(library.shape[0]):
        trace_1d = library[rid].astype(np.float32)
        trace_2d = trace_1d.reshape(n_yrs, 12)
        chars = compute_all_candidates(
            trace_1d, trace_2d, ssi3_calc, ssi12_calc, q80,
            full_record_refs=refs,
        )
        chars["trace_id"] = int(rid)
        rows.append(chars)
    wall = time.perf_counter() - t0
    print(f"[recharacterize] characterisation done in {wall:.1f}s")

    chars_df = pd.DataFrame(rows)
    metric_cols = [m for m in CANDIDATE_METRIC_NAMES if m in chars_df.columns]
    values = chars_df[metric_cols].astype(np.float32).values
    n_events = chars_df["n_events"].astype(np.int32).values \
        if "n_events" in chars_df.columns else np.array([], dtype=np.int32)
    n_events_ssi12 = chars_df["n_events_ssi12"].astype(np.int32).values \
        if "n_events_ssi12" in chars_df.columns else np.array([], dtype=np.int32)

    np.savez_compressed(
        out / "characteristics_extended.npz",
        metric_names=np.array(metric_cols),
        values=values,
        n_events=n_events,
        n_events_ssi12=n_events_ssi12,
        trace_ids=chars_df["trace_id"].astype(np.int32).values,
    )
    print(f"[recharacterize] wrote characteristics_extended.npz "
          f"({values.shape[0]} traces × {values.shape[1]} metrics)")


if __name__ == "__main__":
    main()
