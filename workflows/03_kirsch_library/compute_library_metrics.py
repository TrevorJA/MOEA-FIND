"""Compute drought-metric characterization for an existing Kirsch library.

Reads ``library.npy`` from
``outputs/03_kirsch_library/build_library/<slug>/`` and writes
``characteristics.json`` + ``characteristics.npz`` to the same
directory. The SSI calculator is fit on the full historical record.

Splits cleanly from ``build_library.py``: the heavy work is the Kirsch
generation, so recomputing metrics (different ``--ssi`` timescale, or
when the metric library gains new entries) does not need a regenerated
ensemble.
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
from src.metrics.objectives import (  # noqa: E402
    flows_to_series,
    compute_ssi_drought_characteristics,
    make_ssi_calculator,
)
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.io_paths.slugs import library_slug  # noqa: E402

STAGE = "03_kirsch_library"
SOURCE_DRIVER = "build_library"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-traces", type=int, required=True)
    p.add_argument("--n-years", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--ssi", type=int, required=True, choices=[1, 3, 6, 12])
    args = p.parse_args()

    slug = library_slug(n_traces=args.n_traces, n_years=args.n_years, seed=args.seed)
    out = stage_output_dir(STAGE, SOURCE_DRIVER, slug, create=False)
    library_path = out / "library.npy"
    if not library_path.exists():
        raise SystemExit(
            f"library.npy missing at {library_path}\n"
            f"Run build_library.py with matching --n-traces/--n-years/--seed first."
        )

    print(f"[03/compute_metrics] loading {library_path}")
    library = np.load(library_path)
    print(f"[03/compute_metrics] library shape={library.shape}")

    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    _monthly_2d, monthly_1d = prepare_data(cache_dir)

    ssi_calc = make_ssi_calculator(timescale=args.ssi)
    hist_series = flows_to_series(monthly_1d, start_date="1950-10-01")
    ssi_calc.fit(hist_series)
    print(f"[03/compute_metrics] SSI-{args.ssi} calculator fit on full record")

    print(f"[03/compute_metrics] characterising {library.shape[0]} traces ...")
    t0 = time.perf_counter()
    characteristics: List[Dict] = []
    for rid in range(library.shape[0]):
        trace_1d = library[rid]
        series = flows_to_series(trace_1d, start_date="2100-01-01")
        ssi = ssi_calc.transform(series)
        chars = compute_ssi_drought_characteristics(ssi)
        chars["trace_id"] = int(rid)
        characteristics.append(chars)
    char_wall = time.perf_counter() - t0
    print(f"[03/compute_metrics] done in {char_wall:.1f}s")

    chars_df = pd.DataFrame(characteristics)
    objective_keys = ["mean_duration", "mean_avg_severity"]
    (out / "characteristics.json").write_text(json.dumps(characteristics, indent=2))
    np.savez_compressed(
        out / "characteristics.npz",
        objectives=chars_df[objective_keys].values,
        objective_keys=objective_keys,
        all_keys=list(chars_df.columns),
        all_values=chars_df.values,
        ssi_timescale=np.int32(args.ssi),
    )
    print(f"[03/compute_metrics] wrote characteristics.{{json,npz}} to {out}")


if __name__ == "__main__":
    main()
