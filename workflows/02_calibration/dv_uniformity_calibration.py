"""Bootstrap calibration of the DV-uniformity constraint (SI B).

Parallel to :mod:`constraint_calibration` but in DV space rather than
hydrologic space. Random Kirsch DV draws are analytically ``U[0, 1]^N``,
and there is no historical-DV bootstrap. The calibration is a single
``n_boot`` draw from ``U[0, 1]^N``, followed by the 95% envelope of the
chosen scalar uniformity statistic (L2-star discrepancy and/or KS vs
``U[0, 1]``).

Outputs, written under
``outputs/02_calibration/dv_uniformity_calibration/<dv_mode>/``:

    - ``calibrated_dv_tolerances.json`` -- consumed by
      :meth:`src.optimization.constraints_dv.DVUniformityConfig.from_calibration_json`.
      Keyed as ``{site_label}_T{T_years}_{statistic}`` so L2-star and KS
      entries coexist.
    - ``bootstrap_samples.npz`` -- one array per statistic name.
    - ``config.json`` -- arguments the script was invoked with.

Plots are produced separately by
``workflows/99_supporting_info_figures/dv_uniformity_calibration.py``.

Usage:
    python workflows/02_calibration/dv_uniformity_calibration.py \\
        --T 20 --n-boot 2000 --dv-mode residual --statistic all
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.optimization.constraints_dv import (  # noqa: E402
    VALID_STATISTICS,
    statistic_fn,
)
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "02_calibration"
DRIVER = "dv_uniformity_calibration"


def random_dv_bootstrap(
    n_dvs: int,
    n_boot: int,
    seed: int,
) -> np.ndarray:
    """Return a ``(n_boot, n_dvs)`` array of i.i.d. ``U[0, 1]`` DV vectors."""
    rng = np.random.default_rng(seed)
    return rng.uniform(0.0, 1.0, size=(n_boot, n_dvs))


def batch_statistic(
    samples: np.ndarray,
    statistic: str,
) -> np.ndarray:
    fn = statistic_fn(statistic)
    out = np.empty(samples.shape[0], dtype=float)
    for i, row in enumerate(samples):
        out[i] = fn(row)
    return out


def ninety_five_percent_envelope(values: np.ndarray) -> Tuple[float, float, float]:
    """Return ``(q025, q975, tolerance)`` for a nonnegative statistic."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    q025, q975 = np.quantile(values, [0.025, 0.975])
    return float(q025), float(q975), float(q975)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T", type=int, default=20,
                   help="Synthetic trace length in years.")
    p.add_argument("--n-boot", type=int, default=2000,
                   help="Number of bootstrap draws.")
    p.add_argument("--dv-mode", choices=["residual", "index"], default="residual",
                   help="DV injection mode. Controls n_dvs: residual => "
                        "12*T DVs; index => 12*(T+1) DVs.")
    p.add_argument("--statistic", default="all",
                   choices=("all", "both") + VALID_STATISTICS,
                   help="Which statistic(s) to calibrate. 'all' (default) "
                        "writes entries for every statistic in VALID_STATISTICS.")
    p.add_argument("--seed", type=int, default=20260420)
    p.add_argument("--site-label", default="cannonsville",
                   help="Top-level label in the output JSON keys.")
    args = p.parse_args()

    n_dvs = 12 * args.T if args.dv_mode == "residual" else 12 * (args.T + 1)
    out = stage_output_dir(STAGE, DRIVER, slug=args.dv_mode)
    (out / "config.json").write_text(json.dumps({
        "stage": STAGE,
        "driver": DRIVER,
        "T_years": args.T,
        "n_dvs": n_dvs,
        "dv_mode": args.dv_mode,
        "n_boot": args.n_boot,
        "statistic": args.statistic,
        "seed": args.seed,
        "site_label": args.site_label,
    }, indent=2))

    which: List[str] = (
        list(VALID_STATISTICS)
        if args.statistic in ("all", "both")
        else [args.statistic]
    )
    print(f"[02/dv_uniformity_calibration] n_dvs={n_dvs} (T={args.T}, "
          f"mode={args.dv_mode}), n_boot={args.n_boot}, statistics={which}")

    print(f"  drawing {args.n_boot} random U[0,1] DV vectors ...")
    samples = random_dv_bootstrap(n_dvs, args.n_boot, args.seed)

    summaries: Dict[str, Dict] = {}
    sample_map: Dict[str, np.ndarray] = {}
    payload: Dict[str, Dict] = {}

    for stat_name in which:
        print(f"  computing {stat_name} ...")
        values = batch_statistic(samples, stat_name)
        q025, q975, tolerance = ninety_five_percent_envelope(values)
        summary = {
            "name": stat_name,
            "n_dvs": n_dvs,
            "T_years": args.T,
            "n_boot": args.n_boot,
            "mean": float(np.nanmean(values)),
            "std": float(np.nanstd(values)),
            "q025": q025,
            "q975": q975,
            "tolerance": tolerance,
        }
        summaries[stat_name] = summary
        sample_map[stat_name] = values
        key = f"{args.site_label}_T{args.T}_{stat_name}"
        payload[key] = {
            "T_years": args.T,
            "n_dvs": n_dvs,
            "dv_mode": args.dv_mode,
            "statistic": stat_name,
            "tolerance": tolerance,
            "details": summary,
        }
        print(f"    mean={summary['mean']:.4g} std={summary['std']:.4g}  "
              f"95%CI=[{q025:.4g}, {q975:.4g}]  tol={tolerance:.4g}")

    json_path = out / "calibrated_dv_tolerances.json"
    json_path.write_text(json.dumps(payload, indent=2))
    print(f"  wrote {json_path}")

    npz_path = out / "bootstrap_samples.npz"
    np.savez(npz_path, **{k: v for k, v in sample_map.items()})
    print(f"  wrote {npz_path}")


if __name__ == "__main__":
    main()
