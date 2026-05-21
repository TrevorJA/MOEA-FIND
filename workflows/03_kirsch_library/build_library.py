"""Kirsch library generation (manuscript Fig 7 baseline).

Builds a large library of standard Kirsch-Nowak synthetic traces using
SynHydro's native ensemble API. Generation only -- drought-metric
characterization lives in the separate ``compute_library_metrics.py``
driver. The split lets metrics be (re)computed without regenerating
traces.

This uses ``KirschGenerator.generate()`` directly -- no MOEA wrapper, no
decision-variable injection. The resulting ensemble represents the
natural sampling distribution of the Kirsch bootstrap generator, which
serves as the "random baseline" for comparison with MOEA-FIND.

Outputs under ``outputs/03_kirsch_library/build_library/<slug>/``:
    - config.json   (invocation record)
    - library.npy   (n_traces x n_months monthly flow matrix, water-year)
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import List

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.hydrology.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.io_paths.slugs import library_slug  # noqa: E402

STAGE = "03_kirsch_library"
DRIVER = "build_library"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-traces", type=int, required=True)
    p.add_argument("--n-years", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    args = p.parse_args()

    slug = library_slug(n_traces=args.n_traces, n_years=args.n_years, seed=args.seed)
    out = stage_output_dir(STAGE, DRIVER, slug)
    (out / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "n_traces": args.n_traces,
        "n_years": args.n_years,
        "seed": args.seed,
    }, indent=2))

    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, _monthly_1d = prepare_data(cache_dir)

    print(f"[03/build_library] fitting Kirsch generator ...")
    generator = build_kirsch_generator(monthly_2d)

    print(f"[03/build_library] generating {args.n_traces} traces x {args.n_years} years ...")
    t0 = time.perf_counter()
    ensemble = generator.generate(
        n_realizations=args.n_traces,
        n_years=args.n_years,
        seed=args.seed,
    )
    gen_wall = time.perf_counter() - t0
    print(f"[03/build_library] generation done in {gen_wall:.1f}s")

    # SynHydro outputs calendar-year order; convert to water-year by
    # rolling each year's monthly axis 3 positions (Jan->Oct start).
    library_rows: List[np.ndarray] = []
    for rid in sorted(ensemble.data_by_realization.keys()):
        trace_df = ensemble.data_by_realization[rid]
        trace_1d = trace_df.values.flatten()
        n_yrs = len(trace_1d) // 12
        trace_2d = trace_1d[:n_yrs * 12].reshape(n_yrs, 12)
        trace_2d = np.roll(trace_2d, 3, axis=1)
        library_rows.append(trace_2d.flatten().astype(np.float32))

    library = np.vstack(library_rows)
    np.save(out / "library.npy", library)

    print(f"[03/build_library] wrote {library.shape[0]} traces to {out}")
    print(f"  library.npy shape={library.shape}")
    print(f"[03/build_library] total wall time: {gen_wall:.1f}s")


if __name__ == "__main__":
    main()
