"""Stage 2 — Kirsch baseline library with full 28-metric characterization.

Variant of :mod:`workflows.03_kirsch_library.build_library` that computes
the full :data:`src.metrics.extended.CANDIDATE_METRIC_NAMES`
candidate set on every Kirsch realization, instead of only the
production ``compute_ssi_drought_characteristics`` keys. Used by Stage-2
of the metric/T joint-justification workflow to compare baseline-Kirsch
metric distributions against historical T-blocks.

Outputs under ``outputs/03_kirsch_library/build_library_extended/<slug>/``:

* ``config.json``               — invocation record.
* ``library.npy``               — ``(n_traces, n_years*12)`` water-year
  monthly flow matrix.
* ``characteristics_extended.npz`` — compact arrays:

    - ``metric_names`` (list[str], length 28)
    - ``values`` (n_traces, 28) ndarray of metric values
    - ``n_events``, ``n_events_ssi12`` (n_traces,) integer event counts
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
from src.hydrology.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.metrics.objectives import (  # noqa: E402
    flows_to_series,
    make_ssi_calculator,
)
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "03_kirsch_library"
DRIVER = "build_library_extended"


def slug(n_traces: int, n_years: int, seed: int) -> str:
    return f"n{n_traces}_t{n_years}_ssi3-12_s{seed}"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-traces", type=int, default=10_000)
    p.add_argument("--n-years", type=int, required=True)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    out = stage_output_dir(STAGE, DRIVER,
                           slug(args.n_traces, args.n_years, args.seed))
    (out / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "n_traces": args.n_traces,
        "n_years": args.n_years,
        "ssi_timescales": [3, 12],
        "seed": args.seed,
    }, indent=2))
    print(f"[03/build_library_extended] T={args.n_years}, "
          f"n_traces={args.n_traces}, seed={args.seed} → {out}")

    # --- Load data and fit generator ---
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    print("[03/build_library_extended] fitting Kirsch generator ...")
    generator = build_kirsch_generator(monthly_2d)

    # --- Pre-fit SSI-3, SSI-12 and Q80 on full historical record (DD-11
    # lock-in: identical calibration to Stage-1 historical-block matrix). ---
    ssi3_calc = make_ssi_calculator(timescale=3)
    ssi3_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    ssi12_calc = make_ssi_calculator(timescale=12)
    ssi12_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    q80 = float(np.percentile(monthly_1d, 20.0))
    refs = FullRecordRefs.from_full_record(monthly_1d)
    print(f"[03/build_library_extended] Q80 (full record): {q80:.2f} cfs")

    # --- Generate ensemble ---
    print(f"[03/build_library_extended] generating "
          f"{args.n_traces} traces × {args.n_years} years ...")
    t0 = time.perf_counter()
    ensemble = generator.generate(
        n_realizations=args.n_traces,
        n_years=args.n_years,
        seed=args.seed,
    )
    gen_wall = time.perf_counter() - t0
    print(f"[03/build_library_extended] generation done in {gen_wall:.1f}s")

    # --- Characterize each trace with the 28-metric library ---
    print("[03/build_library_extended] computing 28-metric characteristics ...")
    t1 = time.perf_counter()
    rows: List[Dict[str, float]] = []
    library_rows: List[np.ndarray] = []
    for rid in sorted(ensemble.data_by_realization.keys()):
        trace_df = ensemble.data_by_realization[rid]
        trace_1d = trace_df.values.flatten()

        # SynHydro outputs calendar-year order; convert to water-year
        # ordering (Oct-start) per the project convention used in
        # build_library.py — this preserves cross-stage comparability.
        n_yrs = len(trace_1d) // 12
        trace_2d = trace_1d[: n_yrs * 12].reshape(n_yrs, 12)
        trace_2d = np.roll(trace_2d, 3, axis=1)
        trace_1d = trace_2d.flatten()
        library_rows.append(trace_1d.astype(np.float32))

        chars = compute_all_candidates(
            trace_1d, trace_2d, ssi3_calc, ssi12_calc, q80,
            full_record_refs=refs,
        )
        chars["trace_id"] = int(rid)
        rows.append(chars)
    char_wall = time.perf_counter() - t1
    print(f"[03/build_library_extended] characterization done in {char_wall:.1f}s")

    chars_df = pd.DataFrame(rows)

    # --- Save library + extended characteristics ---
    library = np.vstack(library_rows)
    np.save(out / "library.npy", library)

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

    print(f"[03/build_library_extended] wrote {len(rows)} traces to {out}")
    print(f"  library.npy shape={library.shape}")
    print(f"  characteristics_extended.npz: "
          f"{values.shape[0]} × {values.shape[1]} metrics")
    print(f"[03/build_library_extended] total wall: "
          f"{gen_wall + char_wall:.1f}s "
          f"({gen_wall:.1f}s gen + {char_wall:.1f}s chars)")


if __name__ == "__main__":
    main()
