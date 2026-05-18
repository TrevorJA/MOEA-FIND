"""Kirsch library generation (manuscript Fig 7 baseline).

Builds a large library of standard Kirsch-Nowak synthetic traces using
SynHydro's native ensemble API. Each trace is characterized by SSI-based
drought metrics and stored for downstream convergence analysis and
subsampling by ``subsample_baseline.py``.

This uses ``KirschGenerator.generate()`` directly -- no MOEA wrapper, no
decision-variable injection. The resulting ensemble represents the
natural sampling distribution of the Kirsch bootstrap generator,
which serves as the "random baseline" for comparison with MOEA-FIND.

Outputs under ``outputs/03_kirsch_library/build_library/<slug>/``:
    - config.json          (invocation record)
    - library.npy          (n_traces x n_months monthly flow matrix)
    - characteristics.json (per-trace drought metrics, full)
    - characteristics.npz  (compact arrays for fast loading)
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.hydrology.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.metrics.objectives import (  # noqa: E402
    flows_to_series,
    compute_ssi_drought_characteristics,
    make_ssi_calculator,
)
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "03_kirsch_library"
DRIVER = "build_library"


def slug(n_traces: int, n_years: int, ssi: int, seed: int) -> str:
    return f"n{n_traces}_t{n_years}_ssi{ssi}_s{seed}"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-traces", type=int, required=True)
    p.add_argument("--n-years", type=int, required=True)
    p.add_argument("--ssi", type=int, required=True, choices=[1, 3, 6, 12])
    p.add_argument("--seed", type=int, required=True)
    args = p.parse_args()

    out = stage_output_dir(STAGE, DRIVER, slug(args.n_traces, args.n_years, args.ssi, args.seed))
    (out / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "n_traces": args.n_traces,
        "n_years": args.n_years,
        "ssi": args.ssi,
        "seed": args.seed,
    }, indent=2))

    # --- Load data and fit generator ---
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    print(f"[03/build_library] fitting Kirsch generator ...")
    generator = build_kirsch_generator(monthly_2d)

    # --- Pre-fit SSI calculator on historical data ---
    ssi_calc = make_ssi_calculator(timescale=args.ssi)
    hist_series = flows_to_series(monthly_1d, start_date="1950-10-01")
    ssi_calc.fit(hist_series)

    # --- Generate ensemble using SynHydro's native API ---
    print(f"[03/build_library] generating {args.n_traces} traces x {args.n_years} years ...")
    t0 = time.perf_counter()
    ensemble = generator.generate(
        n_realizations=args.n_traces,
        n_years=args.n_years,
        seed=args.seed,
    )
    gen_wall = time.perf_counter() - t0
    print(f"[03/build_library] generation done in {gen_wall:.1f}s")

    # --- Characterize each trace and collect monthly flow matrix ---
    print(f"[03/build_library] computing drought characteristics ...")
    t1 = time.perf_counter()
    characteristics: List[Dict] = []
    library_rows: List[np.ndarray] = []
    for rid in sorted(ensemble.data_by_realization.keys()):
        trace_df = ensemble.data_by_realization[rid]
        trace_1d = trace_df.values.flatten()

        # SynHydro outputs calendar-year order; convert to water-year
        n_yrs = len(trace_1d) // 12
        trace_2d = trace_1d[:n_yrs * 12].reshape(n_yrs, 12)
        trace_2d = np.roll(trace_2d, 3, axis=1)
        trace_1d = trace_2d.flatten()
        library_rows.append(trace_1d.astype(np.float32))

        series = flows_to_series(trace_1d, start_date="2100-01-01")
        ssi = ssi_calc.transform(series)
        chars = compute_ssi_drought_characteristics(ssi)
        chars["trace_id"] = int(rid)
        characteristics.append(chars)

    char_wall = time.perf_counter() - t1
    print(f"[03/build_library] characterization done in {char_wall:.1f}s")

    # --- Save monthly flow library as npy ---
    library = np.vstack(library_rows)
    np.save(out / "library.npy", library)

    # --- Save JSON (full characteristics) ---
    (out / "characteristics.json").write_text(
        json.dumps(characteristics, indent=2)
    )

    # --- Save NPZ (compact arrays for fast loading in convergence analysis) ---
    import pandas as pd
    chars_df = pd.DataFrame(characteristics)
    objective_keys = ["mean_duration", "mean_avg_severity"]
    np.savez_compressed(
        out / "characteristics.npz",
        objectives=chars_df[objective_keys].values,
        objective_keys=objective_keys,
        all_keys=list(chars_df.columns),
        all_values=chars_df.values,
    )

    print(f"[03/build_library] wrote {len(characteristics)} traces to {out}")
    print(f"  library.npy shape={library.shape}")
    print(f"  characteristics.json (full) + characteristics.npz (compact)")
    print(f"[03/build_library] total wall time: {gen_wall + char_wall:.1f}s "
          f"({gen_wall:.1f}s generation + {char_wall:.1f}s characterization)")


if __name__ == "__main__":
    main()
