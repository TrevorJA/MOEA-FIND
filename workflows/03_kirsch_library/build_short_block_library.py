"""DD-15b — Kirsch ensemble for short-block screening at T ∈ {1, 2}.

Generates one set of n_years=3 Kirsch traces and extracts metrics for
both T=1 and T=2 from the same library, with 3-month SSI-3 burn-in
prepended to each evaluation block. This saves one generation pass.

Layout per 36-month trace (water-year ordered after roll-by-3):

  months  0-11 = year 1 (used as full burn-in for T=2 trio)
  months 12-35 = T=2 evaluation block (24 months)
  months  9-35 = T=2 SSI-3 input window (27 months: 3 burnin + 24 eval)
  months 24-35 = T=1 evaluation block (12 months)
  months 21-35 = T=1 SSI-3 input window (15 months: 3 burnin + 12 eval)

Outputs (under
``outputs/03_kirsch_library/build_short_block_library/n{n_traces}_s{seed}/``):

* ``library_36mo.npy``                — (n_traces, 36) monthly flow matrix
* ``characteristics_short_T1.npz``    — Tier H + Tier I metrics on T=1 eval block
* ``characteristics_short_T2.npz``    — same on T=2 eval block
* ``config.json``
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
from src.metrics.extended import FullRecordRefs  # noqa: E402
from src.hydrology.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.metrics.objectives import flows_to_series, make_ssi_calculator  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.metrics.short_block import (  # noqa: E402
    SHORT_BLOCK_METRIC_NAMES,
    compute_short_block_metrics,
)

STAGE = "03_kirsch_library"
DRIVER = "build_short_block_library"
GEN_N_YEARS = 3
BURNIN_MONTHS = 3


def _slug(n_traces: int, seed: int) -> str:
    return f"n{n_traces}_s{seed}"


def _save_chars_npz(
    rows: List[Dict[str, float]], out_path: Path,
):
    df = pd.DataFrame(rows)
    metric_cols = [m for m in SHORT_BLOCK_METRIC_NAMES if m in df.columns]
    values = df[metric_cols].astype(np.float32).values
    np.savez_compressed(
        out_path,
        metric_names=np.array(metric_cols),
        values=values,
        trace_ids=np.arange(values.shape[0], dtype=np.int32),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-traces", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    out = stage_output_dir(STAGE, DRIVER, slug=_slug(args.n_traces, args.seed))
    print(f"[03/{DRIVER}] n_traces={args.n_traces}, seed={args.seed} → {out}")

    cache = PROJECT_ROOT / "outputs" / "data_cache"
    cache.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache)

    print("[diag] fitting Kirsch generator ...")
    generator = build_kirsch_generator(monthly_2d)

    ssi3_calc = make_ssi_calculator(timescale=3)
    ssi3_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    refs = FullRecordRefs.from_full_record(monthly_1d)

    print(f"[diag] generating {args.n_traces} traces × {GEN_N_YEARS} years ...")
    t0 = time.perf_counter()
    ensemble = generator.generate(
        n_realizations=args.n_traces,
        n_years=GEN_N_YEARS,
        seed=args.seed,
    )
    gen_wall = time.perf_counter() - t0
    print(f"[diag] generation done in {gen_wall:.1f}s")

    # --- Convert to water-year ordering and stack into one matrix ---
    print("[diag] reshaping to water-year ordering ...")
    library_rows: List[np.ndarray] = []
    for rid in sorted(ensemble.data_by_realization.keys()):
        trace_df = ensemble.data_by_realization[rid]
        trace_1d = np.asarray(trace_df.values).flatten()
        n_yrs = len(trace_1d) // 12
        trace_2d = trace_1d[: n_yrs * 12].reshape(n_yrs, 12)
        # Calendar Jan-Dec → water-year Oct-Sep (np.roll on month axis)
        trace_2d = np.roll(trace_2d, 3, axis=1)
        library_rows.append(trace_2d.flatten().astype(np.float32))
    library = np.vstack(library_rows)
    np.save(out / "library_36mo.npy", library)
    print(f"[diag] library_36mo.npy shape={library.shape}")

    # --- Characterise each trace at T=1 and T=2 ---
    print("[diag] characterising at T=1 (months 24..35) "
          "and T=2 (months 12..35) ...")
    t1 = time.perf_counter()
    rows_T1: List[Dict[str, float]] = []
    rows_T2: List[Dict[str, float]] = []
    for rid in range(library.shape[0]):
        trace = library[rid]

        # T=2 evaluation: months 12..35; SSI-3 input: months 9..35
        eval_T2 = trace[12:36]
        ssi_input_T2 = trace[9:36]
        ssi_series_T2 = ssi3_calc.transform(
            flows_to_series(ssi_input_T2, start_date="2100-01-01")
        )
        chars_T2 = compute_short_block_metrics(
            eval_T2, ssi_series_T2, refs,
            eval_first_idx_in_ssi=BURNIN_MONTHS,
        )
        chars_T2["trace_id"] = int(rid)
        rows_T2.append(chars_T2)

        # T=1 evaluation: months 24..35; SSI-3 input: months 21..35
        eval_T1 = trace[24:36]
        ssi_input_T1 = trace[21:36]
        ssi_series_T1 = ssi3_calc.transform(
            flows_to_series(ssi_input_T1, start_date="2100-01-01")
        )
        chars_T1 = compute_short_block_metrics(
            eval_T1, ssi_series_T1, refs,
            eval_first_idx_in_ssi=BURNIN_MONTHS,
        )
        chars_T1["trace_id"] = int(rid)
        rows_T1.append(chars_T1)

    char_wall = time.perf_counter() - t1
    print(f"[diag] characterisation done in {char_wall:.1f}s")

    _save_chars_npz(rows_T1, out / "characteristics_short_T1.npz")
    _save_chars_npz(rows_T2, out / "characteristics_short_T2.npz")
    print(f"[diag] wrote characteristics_short_T1.npz "
          f"and characteristics_short_T2.npz")

    cfg = {
        "n_traces": int(args.n_traces),
        "seed": int(args.seed),
        "gen_n_years": GEN_N_YEARS,
        "burnin_months": BURNIN_MONTHS,
        "metric_names": list(SHORT_BLOCK_METRIC_NAMES),
        "T_evaluation_layout": {
            "T1": {"eval_months": [24, 36], "ssi_input_months": [21, 36]},
            "T2": {"eval_months": [12, 36], "ssi_input_months": [9, 36]},
        },
        "wall_seconds": {
            "generation": gen_wall, "characterization": char_wall,
        },
    }
    (out / "config.json").write_text(json.dumps(cfg, indent=2))
    print(f"[diag] total wall: {gen_wall + char_wall:.1f}s")


if __name__ == "__main__":
    main()
