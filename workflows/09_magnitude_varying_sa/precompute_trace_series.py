"""Stage 09 -- precompute per-trace annual operational-outcome timeseries.

Reads the Stage-06 ``simulations/pywrdrb_output.hdf5`` and writes a
parquet of per-trace annual NYC min storage fraction (rows = calendar
year, columns = realization id). This is the input the within-trace-
percentile MV-SA response form (Hadjimichael Variant 2) consumes.

One-shot driver: cheap (single hdf5 read, daily-to-annual reduction),
runs in seconds. Submit via SLURM only if the hdf5 is on slow storage.

Output:

    outputs/09_magnitude_varying_sa/precompute_trace_series/<src_slug>/
        annual_nyc_min_storage_frac.parquet

Usage:
    python workflows/09_magnitude_varying_sa/precompute_trace_series.py \\
        --pywrdrb-output outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/simulations/pywrdrb_output.hdf5 \\
        --src-slug <src_slug>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.io import save_experiment_config  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.hydrology.precompute_trace_series import (  # noqa: E402
    compute_annual_nyc_min_storage_frac,
)

STAGE = "09_magnitude_varying_sa"
DRIVER = "precompute_trace_series"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--pywrdrb-output", type=Path, required=True,
                   help="Path to Stage-06 simulations/pywrdrb_output.hdf5.")
    p.add_argument("--src-slug", required=True,
                   help="Upstream archive slug; used to name the output "
                        "directory (so the SA driver can find it).")
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER, args.src_slug)
    out_path = out_dir / "annual_nyc_min_storage_frac.parquet"

    print(f"[09/precompute] reading {args.pywrdrb_output}", flush=True)
    df = compute_annual_nyc_min_storage_frac(args.pywrdrb_output)
    print(f"[09/precompute] computed shape={df.shape} "
          f"(years x traces)", flush=True)
    print(f"[09/precompute] storage frac range "
          f"min={df.values.min():.4f} max={df.values.max():.4f} "
          f"mean={df.values.mean():.4f}", flush=True)

    df.to_parquet(out_path)
    save_experiment_config(out_dir, {
        "script": "workflows/09_magnitude_varying_sa/precompute_trace_series.py",
        "stage": STAGE,
        "driver": DRIVER,
        "src_slug": args.src_slug,
        "pywrdrb_output": str(args.pywrdrb_output),
        "n_years": int(df.shape[0]),
        "n_traces": int(df.shape[1]),
        "outcome": "annual_nyc_min_storage_frac",
    })
    print(f"[09/precompute] wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
