"""Bootstrap calibration for the DV-space uniformity constraint (SI ablation).

Parallel to ``diag_constraint_calibration.py`` but in DV space instead of
hydrologic space. Produces calibrated tolerances for both available
uniformity statistics (L2-star discrepancy, KS vs U[0, 1]) at the target
DV-vector length ``N = n_years_out * 12`` (residual mode) so the DV-uniform
ablation arm can be configured against an empirical 95% envelope.

Calibration asymmetry vs the hydrologic script. Random Kirsch DV draws are
analytically ``U[0, 1]^N``, and there is no "historical DV" bootstrap (no DVs
exist for historical observations). The calibration is therefore a single
``n_boot`` draw from ``U[0, 1]^N``, followed by the 95% envelope of the chosen
scalar statistic — no min/max over two distributions like the hydrologic
script. The script still generates a diagnostic PDF comparable in structure
to ``constraint_calibration.pdf`` so reviewers can eyeball both calibrations
the same way.

Outputs, under ``outputs/diag_dv_uniformity_calibration/``:

    - ``calibrated_dv_tolerances.json`` — consumed by
      :meth:`src.constraints_dv.DVUniformityConfig.from_calibration_json`.
      Keyed as ``{site_label}_T{T_years}_{statistic}`` so L2-star and KS
      entries coexist.
    - ``dv_calibration.pdf`` — one histogram panel per statistic.
    - ``config.json`` — the arguments the script was invoked with.

Usage:
    python workflows/diagnostics/diag_dv_uniformity_calibration.py \\
        --T 20 --n-boot 2000 --statistic l2_star
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

from src.constraints_dv import (  # noqa: E402
    VALID_STATISTICS,
    statistic_fn,
)

OUTPUT_SLUG = "diag_dv_uniformity_calibration"


# ---------------------------------------------------------------------------
# Bootstrap driver
# ---------------------------------------------------------------------------

def random_dv_bootstrap(
    n_dvs: int,
    n_boot: int,
    seed: int,
) -> np.ndarray:
    """Return a ``(n_boot, n_dvs)`` array of i.i.d. ``U[0, 1]`` DV vectors.

    Mirrors what the default (non-MOEA) Kirsch generator would produce:
    every DV is an independent uniform draw, which is exactly the null
    distribution the DV-uniformity constraint should accept.
    """
    rng = np.random.default_rng(seed)
    return rng.uniform(0.0, 1.0, size=(n_boot, n_dvs))


def batch_statistic(
    samples: np.ndarray,
    statistic: str,
) -> np.ndarray:
    """Apply ``statistic`` to each row of ``samples`` and return a 1D array."""
    fn = statistic_fn(statistic)
    out = np.empty(samples.shape[0], dtype=float)
    for i, row in enumerate(samples):
        out[i] = fn(row)
    return out


def ninety_five_percent_envelope(values: np.ndarray) -> Tuple[float, float, float]:
    """Return ``(q025, q975, tolerance)`` for a nonnegative statistic.

    For the DV-uniformity statistics, the reference is the 95% upper bound
    (the statistic cannot go below 0, so the 2.5% quantile is not a
    meaningful tolerance edge). ``tolerance`` here is the 97.5th percentile.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    q025, q975 = np.quantile(values, [0.025, 0.975])
    return float(q025), float(q975), float(q975)


# ---------------------------------------------------------------------------
# Diagnostic plotting
# ---------------------------------------------------------------------------

def make_diagnostic_pdf(
    summaries: Dict[str, Dict],
    sample_map: Dict[str, np.ndarray],
    output_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = list(summaries.keys())
    n = len(names)
    fig, axes = plt.subplots(n, 1, figsize=(8, 2.6 * n), constrained_layout=True)
    if n == 1:
        axes = [axes]

    for ax, name in zip(axes, names):
        s = summaries[name]
        x = sample_map[name]
        x = x[np.isfinite(x)]
        ax.hist(x, bins=40, alpha=0.55, color="#c05621",
                label=f"Random U[0,1] bootstrap (n={len(x)})",
                density=True)
        tol = s["tolerance"]
        ax.axvline(tol, color="black", linestyle="--", lw=1,
                   label=f"tol (97.5%) = {tol:.4g}")
        ax.axvline(s["q025"], color="gray", linestyle=":", lw=1,
                   label=f"2.5% = {s['q025']:.4g}")
        ax.axvline(s["q975"], color="gray", linestyle=":", lw=1)
        ax.set_title(f"{name} on random U[0,1] DV draws (N={s['n_dvs']})")
        ax.set_xlabel(name)
        ax.set_ylabel("density")
        ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(
        "DV-uniformity bootstrap calibration (random U[0,1] DV draws)",
        fontsize=12,
    )
    fig.savefig(output_path, format="pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

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
                        "writes entries for every statistic in "
                        "VALID_STATISTICS. 'both' is an alias kept for "
                        "back-compat with earlier scripts.")
    p.add_argument("--seed", type=int, default=20260420)
    p.add_argument("--site-label", default="cannonsville",
                   help="Top-level label in the output JSON keys.")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    n_dvs = 12 * args.T if args.dv_mode == "residual" else 12 * (args.T + 1)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    (out / "config.json").write_text(json.dumps({
        "script": "diag_dv_uniformity_calibration.py",
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
    print(f"[dv-calib] n_dvs={n_dvs} (T={args.T}, mode={args.dv_mode}), "
          f"n_boot={args.n_boot}, statistics={which}")

    print(f"[dv-calib] drawing {args.n_boot} random U[0,1] DV vectors ...")
    samples = random_dv_bootstrap(n_dvs, args.n_boot, args.seed)

    summaries: Dict[str, Dict] = {}
    sample_map: Dict[str, np.ndarray] = {}
    payload: Dict[str, Dict] = {}

    for stat_name in which:
        print(f"[dv-calib] computing {stat_name} statistic ...")
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
        print(f"[dv-calib]   mean={summary['mean']:.4g} "
              f"std={summary['std']:.4g}  "
              f"95%CI=[{q025:.4g}, {q975:.4g}]  "
              f"tol={tolerance:.4g}")

    json_path = out / "calibrated_dv_tolerances.json"
    json_path.write_text(json.dumps(payload, indent=2))
    print(f"[dv-calib] wrote {json_path}")

    pdf_path = out / "dv_calibration.pdf"
    make_diagnostic_pdf(summaries, sample_map, pdf_path)
    print(f"[dv-calib] wrote {pdf_path}")


if __name__ == "__main__":
    main()
