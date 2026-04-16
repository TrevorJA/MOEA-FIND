"""Script 05 — Kirsch library generation (manuscript §6.3, Fig 7 baseline).

Builds a large library of standard Kirsch-Nowak synthetic traces using
SynHydro's native ensemble API. Each trace is characterized by SSI-based
drought metrics and stored for downstream convergence analysis and
subsampling by script 06.

This uses KirschGenerator.generate() directly — no MOEA wrapper, no
decision-variable injection. The resulting ensemble represents the
natural sampling distribution of the Kirsch bootstrap generator,
which serves as the "random baseline" for comparison with MOEA-FIND.

Outputs under outputs/exp05_kirsch_library/:
    - characteristics.json (per-trace drought metrics)
    - config.json

Run locally (small test):
    python workflows/experiments/05_kirsch_library_build.py --n-traces 200

Run on HPC:
    sbatch workflows/slurm/05_kirsch_library_build.slurm
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

from src.experiment_utils import prepare_data, compute_historical_ssi_chars  # noqa: E402
from src.objectives import (  # noqa: E402
    flows_to_series,
    compute_ssi_drought_characteristics,
)

OUTPUT_SLUG = "exp05_kirsch_library"


def _build_generator(monthly_2d: np.ndarray):
    """Fit SynHydro KirschGenerator on historical monthly flows."""
    import pandas as pd
    from synhydro.methods.generation.nonparametric.kirsch import KirschGenerator
    gen = KirschGenerator(generate_using_log_flow=True)
    dates = pd.date_range(start="1950-10-01", periods=monthly_2d.size, freq="MS")
    gen.fit(pd.DataFrame({"flow_cfs": monthly_2d.flatten()}, index=dates))
    return gen


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-traces", type=int, default=10_000)
    p.add_argument("--n-years", type=int, default=20)
    p.add_argument("--ssi", type=int, default=3, choices=[1, 3, 6, 12])
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    (out / "config.json").write_text(json.dumps({
        "script": "05_kirsch_library_build.py",
        "n_traces": args.n_traces,
        "n_years": args.n_years,
        "ssi": args.ssi,
        "seed": args.seed,
    }, indent=2))

    # --- Load data and fit generator ---
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    print(f"[05] fitting Kirsch generator ...")
    generator = _build_generator(monthly_2d)

    # --- Pre-fit SSI calculator on historical data ---
    from src.objectives import make_ssi_calculator
    ssi_calc = make_ssi_calculator(timescale=args.ssi)
    hist_series = flows_to_series(monthly_1d, start_date="1950-10-01")
    ssi_calc.fit(hist_series)

    # --- Generate ensemble using SynHydro's native API ---
    print(f"[05] generating {args.n_traces} traces x {args.n_years} years ...")
    t0 = time.perf_counter()
    ensemble = generator.generate(
        n_realizations=args.n_traces,
        n_years=args.n_years,
        seed=args.seed,
    )
    gen_wall = time.perf_counter() - t0
    print(f"[05] generation done in {gen_wall:.1f}s")

    # --- Characterize each trace ---
    print(f"[05] computing drought characteristics ...")
    t1 = time.perf_counter()
    characteristics: List[Dict] = []
    for rid in sorted(ensemble.data_by_realization.keys()):
        trace_df = ensemble.data_by_realization[rid]
        trace_1d = trace_df.values.flatten()

        # SynHydro outputs calendar-year order; convert to water-year
        n_yrs = len(trace_1d) // 12
        trace_2d = trace_1d[:n_yrs * 12].reshape(n_yrs, 12)
        trace_2d = np.roll(trace_2d, 3, axis=1)
        trace_1d = trace_2d.flatten()

        series = flows_to_series(trace_1d, start_date="2100-01-01")
        ssi = ssi_calc.transform(series)
        chars = compute_ssi_drought_characteristics(ssi)
        chars["trace_id"] = int(rid)
        characteristics.append(chars)

    char_wall = time.perf_counter() - t1
    print(f"[05] characterization done in {char_wall:.1f}s")

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

    print(f"[05] wrote {len(characteristics)} traces to {out}")
    print(f"     characteristics.json (full) + characteristics.npz (compact)")
    print(f"[05] total wall time: {gen_wall + char_wall:.1f}s "
          f"({gen_wall:.1f}s generation + {char_wall:.1f}s characterization)")


if __name__ == "__main__":
    main()
