"""Stage 3 — per-T MOEA evaluation-cost timing.

For each T in the sweep grid, time N sequential calls to the
production MOEA inner loop (`KirschBorgWrapper.generate()` →
flatten → SSI transform → ``compute_ssi_drought_characteristics``)
on a single rank with random DVs in [0, 1]. This mirrors the
per-evaluation cost the MM Borg run will incur at that T.

Output (``outputs/02_calibration/eval_cost_timing/T{T:02d}/``):
``timing.json`` with median, mean, std, min, max wall-time per
evaluation and the parameters of the fit.

The Stage-3 aggregator :mod:`workflows.02_calibration.decision_matrix`
loads these to fit ``t_eval(T) = a + b·T`` and project total MOEA
wall-time at NFE=200,000 × 120 ranks for each candidate T.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.hydrology.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.hydrology.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.metrics.objectives import (  # noqa: E402
    compute_ssi_drought_characteristics,
    flows_to_series,
    make_ssi_calculator,
)
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "02_calibration"
DRIVER = "eval_cost_timing"


def time_evaluations(
    wrapper: KirschBorgWrapper,
    ssi_calc,
    n_warmup: int,
    n_eval: int,
    seed: int,
) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    n_dvs = wrapper.n_dvs

    # Warmup (JIT, allocator, OS page-cache).
    for _ in range(n_warmup):
        dvs = rng.uniform(0.0, 1.0, n_dvs)
        synth_2d = wrapper.generate(dvs)
        synth_1d = synth_2d.flatten() if synth_2d.ndim > 1 else synth_2d
        series = flows_to_series(synth_1d, start_date="2100-01-01")
        ssi = ssi_calc.transform(series)
        compute_ssi_drought_characteristics(ssi, monthly_flows=synth_1d)

    # Timed loop.
    walls = np.empty(n_eval, dtype=np.float64)
    for i in range(n_eval):
        dvs = rng.uniform(0.0, 1.0, n_dvs)
        t0 = time.perf_counter()
        synth_2d = wrapper.generate(dvs)
        synth_1d = synth_2d.flatten() if synth_2d.ndim > 1 else synth_2d
        series = flows_to_series(synth_1d, start_date="2100-01-01")
        ssi = ssi_calc.transform(series)
        compute_ssi_drought_characteristics(ssi, monthly_flows=synth_1d)
        walls[i] = time.perf_counter() - t0

    return {
        "n_eval": int(n_eval),
        "mean": float(np.mean(walls)),
        "median": float(np.median(walls)),
        "std": float(np.std(walls, ddof=1)),
        "min": float(np.min(walls)),
        "max": float(np.max(walls)),
        "p10": float(np.percentile(walls, 10)),
        "p90": float(np.percentile(walls, 90)),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-years", type=int, required=True)
    p.add_argument("--mode", choices=["index", "residual"], default="residual")
    p.add_argument("--n-warmup", type=int, default=20)
    p.add_argument("--n-eval", type=int, default=200)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    slug = f"T{args.T_years:02d}"
    out_dir = stage_output_dir(STAGE, DRIVER, slug=slug)
    print(f"[02/{DRIVER}] T={args.T_years} mode={args.mode} → {out_dir}")

    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    print("[diag] fitting Kirsch generator ...")
    generator = build_kirsch_generator(monthly_2d)
    wrapper = KirschBorgWrapper(generator, mode=args.mode,
                                n_years_out=args.T_years)
    print(f"[diag] wrapper n_dvs={wrapper.n_dvs}")

    print("[diag] fitting SSI-3 on full record ...")
    ssi_calc = make_ssi_calculator(timescale=3)
    ssi_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))

    print(f"[diag] timing: warmup={args.n_warmup}, n_eval={args.n_eval}")
    stats = time_evaluations(
        wrapper, ssi_calc,
        n_warmup=args.n_warmup,
        n_eval=args.n_eval,
        seed=args.seed,
    )
    payload = {
        "T_years": int(args.T_years),
        "mode": args.mode,
        "n_dvs": int(wrapper.n_dvs),
        **stats,
        "n_warmup": int(args.n_warmup),
        "seed": int(args.seed),
    }
    (out_dir / "timing.json").write_text(json.dumps(payload, indent=2))
    print(f"[diag] median={stats['median']*1000:.2f} ms, "
          f"mean={stats['mean']*1000:.2f} ms ± "
          f"{stats['std']*1000:.2f} ms")
    print(f"[diag] wrote {out_dir / 'timing.json'}")


if __name__ == "__main__":
    main()
