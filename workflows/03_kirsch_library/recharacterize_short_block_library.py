"""DD-15c — recharacterize the existing short-block Kirsch library
without regenerating traces.

Reads the 36-month water-year-aligned monthly-flow library
(``library_36mo.npy``) saved by
:mod:`workflows.03_kirsch_library.build_short_block_library` and
re-extracts the (expanded) short-block metric pool with the current
:func:`src.metrics.short_block.compute_short_block_metrics`. Writes
refreshed ``characteristics_short_T1.npz`` and
``characteristics_short_T2.npz`` next to the input library.

Cheap (~5 min for 10 000 traces) compared to regeneration (~20 min).
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
from src.metrics.objectives import flows_to_series, make_ssi_calculator  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.metrics.short_block import (  # noqa: E402
    SHORT_BLOCK_METRIC_NAMES,
    compute_short_block_metrics,
)

STAGE = "03_kirsch_library"
DRIVER = "build_short_block_library"
BURNIN_MONTHS = 3


def _save_chars_npz(rows: List[Dict[str, float]], out_path: Path):
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

    out = stage_output_dir(
        STAGE, DRIVER, slug=f"n{args.n_traces}_s{args.seed}", create=False,
    )
    library_path = out / "library_36mo.npy"
    if not library_path.exists():
        raise SystemExit(f"library_36mo.npy missing at {library_path}")
    print(f"[recharacterize_short] reading {library_path}")
    library = np.load(library_path)
    print(f"[recharacterize_short] library shape={library.shape}")

    cache = PROJECT_ROOT / "outputs" / "data_cache"
    cache.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache)
    ssi3_calc = make_ssi_calculator(timescale=3)
    ssi3_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    refs = FullRecordRefs.from_full_record(monthly_1d)

    print(f"[recharacterize_short] characterising {library.shape[0]} traces"
          f" with {len(SHORT_BLOCK_METRIC_NAMES)} metrics ...")
    t0 = time.perf_counter()
    rows_T1: List[Dict[str, float]] = []
    rows_T2: List[Dict[str, float]] = []
    for rid in range(library.shape[0]):
        trace = library[rid]
        # T=2: months 12..35 evaluation, 9..35 SSI window
        eval_T2 = trace[12:36]
        ssi_input_T2 = trace[9:36]
        ssi_T2 = ssi3_calc.transform(
            flows_to_series(ssi_input_T2, start_date="2100-01-01")
        )
        chars_T2 = compute_short_block_metrics(
            eval_T2, ssi_T2, refs, eval_first_idx_in_ssi=BURNIN_MONTHS,
        )
        chars_T2["trace_id"] = int(rid)
        rows_T2.append(chars_T2)

        # T=1: months 24..35 evaluation, 21..35 SSI window
        eval_T1 = trace[24:36]
        ssi_input_T1 = trace[21:36]
        ssi_T1 = ssi3_calc.transform(
            flows_to_series(ssi_input_T1, start_date="2100-01-01")
        )
        chars_T1 = compute_short_block_metrics(
            eval_T1, ssi_T1, refs, eval_first_idx_in_ssi=BURNIN_MONTHS,
        )
        chars_T1["trace_id"] = int(rid)
        rows_T1.append(chars_T1)
    wall = time.perf_counter() - t0
    print(f"[recharacterize_short] characterisation done in {wall:.1f}s")

    _save_chars_npz(rows_T1, out / "characteristics_short_T1.npz")
    _save_chars_npz(rows_T2, out / "characteristics_short_T2.npz")
    print(f"[recharacterize_short] wrote characteristics_short_T1.npz "
          f"and characteristics_short_T2.npz "
          f"({len(SHORT_BLOCK_METRIC_NAMES)} metrics)")


if __name__ == "__main__":
    main()
