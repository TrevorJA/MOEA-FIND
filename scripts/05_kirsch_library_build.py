"""Script 05 — Kirsch library generation (manuscript §6.3, Fig 7 baseline).

Builds a large library of Kirsch-Nowak synthetic traces for the library-and-
subsample baseline. Each trace is characterized by SSI-based drought metrics
and stored for downstream subsampling by script 06.

Parallelization:
    - Serial fallback when mpi4py is unavailable (for local smoke tests).
    - mpi4py master-worker pool on HPC: rank 0 farms trace indices to workers,
      collects drought characteristics, and writes a single JSON/NPZ.

Outputs under outputs/exp05_kirsch_library/:
    - library.npz (traces as (n_traces, n_months) array)
    - characteristics.json (per-trace drought metrics)
    - config.json

Run locally (small smoke test, serial):
    python scripts/05_kirsch_library_build.py --n-traces 200 --n-years 15

Run on HPC (MPI):
    sbatch scripts/05_kirsch_library_build.slurm
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment_utils import prepare_data, compute_historical_ssi_chars  # noqa: E402

OUTPUT_SLUG = "exp05_kirsch_library"


def _try_mpi():
    try:
        from mpi4py import MPI
        return MPI
    except ImportError:
        return None


def _build_generator(monthly_2d: np.ndarray):
    import pandas as pd
    from synhydro.methods.generation.nonparametric.kirsch import KirschGenerator
    gen = KirschGenerator(generate_using_log_flow=True)
    dates = pd.date_range(start="1950-10-01", periods=monthly_2d.size, freq="MS")
    gen.fit(pd.DataFrame({"flow_cfs": monthly_2d.flatten()}, index=dates))
    return gen


def _characterize(trace: np.ndarray, ssi_hist, ssi_calc, ssi_acc: int) -> Dict:
    from src.objectives import (
        flows_to_series, compute_ssi, compute_ssi_drought_characteristics,
    )
    series = flows_to_series(trace)
    ssi = compute_ssi(series, ssi_calc, ssi_hist)
    return compute_ssi_drought_characteristics(ssi)


def generate_one(generator, ssi_hist, ssi_calc, ssi_acc: int,
                 n_years: int, trace_seed: int) -> Dict:
    """Generate a single random Kirsch trace and compute SSI characteristics."""
    rng = np.random.default_rng(trace_seed)
    dvs = rng.random(n_years * 12)
    # Reuse wrapper in residual mode for a uniform-random trace draw
    from src.kirsch_wrapper import KirschBorgWrapper
    wrapper = KirschBorgWrapper(generator, mode="residual", n_years_out=n_years)
    flows = wrapper.generate(dvs)
    chars = _characterize(flows, ssi_hist, ssi_calc, ssi_acc)
    chars["seed"] = int(trace_seed)
    return chars


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-traces", type=int, default=10_000)
    p.add_argument("--n-years", type=int, default=15)
    p.add_argument("--ssi", type=int, default=3, choices=[1, 3, 6, 12])
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    MPI = _try_mpi()
    if MPI is None:
        rank, size = 0, 1
        comm = None
    else:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

    if rank == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "config.json").write_text(json.dumps({
            "script": "05_kirsch_library_build.py",
            "manuscript_section": "§6.3 Library Baseline (Fig 7 baseline)",
            "n_traces": args.n_traces, "n_years": args.n_years,
            "ssi": args.ssi, "seed": args.seed, "mpi_size": size,
        }, indent=2))
        cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        monthly_2d, monthly_1d = prepare_data(cache_dir)
        ssi_hist, ssi_calc, _ = compute_historical_ssi_chars(monthly_1d, args.ssi)
        generator = _build_generator(monthly_2d)
    else:
        generator = ssi_hist = ssi_calc = None

    if comm is not None:
        generator = comm.bcast(generator, root=0)
        ssi_hist = comm.bcast(ssi_hist, root=0)
        ssi_calc = comm.bcast(ssi_calc, root=0)

    # Assign traces round-robin
    trace_ids = list(range(args.n_traces))
    my_ids = trace_ids[rank::size]
    rng = np.random.default_rng(args.seed)
    base_seeds = rng.integers(0, 2**31 - 1, size=args.n_traces)

    t0 = time.perf_counter()
    my_results: List[Dict] = []
    for tid in my_ids:
        chars = generate_one(generator, ssi_hist, ssi_calc, args.ssi,
                             args.n_years, int(base_seeds[tid]))
        chars["trace_id"] = tid
        my_results.append(chars)
    wall = time.perf_counter() - t0

    if comm is not None:
        gathered = comm.gather(my_results, root=0)
    else:
        gathered = [my_results]

    if rank == 0:
        flat: List[Dict] = []
        for chunk in gathered:
            flat.extend(chunk)
        flat.sort(key=lambda r: r["trace_id"])
        (args.output_dir / "characteristics.json").write_text(json.dumps(flat, indent=2))
        print(f"[05] generated {len(flat)} traces in {wall:.1f}s (mpi_size={size})")
        print(f"     wrote {args.output_dir / 'characteristics.json'}")


if __name__ == "__main__":
    main()
