"""Stage 04 / baseline_comparison -- Kirsch library vs MOEA-FIND (Fig 7).

Computes coverage metrics for the stage 03 Kirsch random library against
the stage 04 MOEA-FIND Pareto front and writes a numerical summary.
Figures are produced separately by
``workflows/04_moea_find_single_site/plots/baseline_comparison.py``.

Outputs under ``outputs/04_moea_find_single_site/baseline_comparison/``:
    config.json
    comparison_summary.json
    pooled.npz                (kirsch + moea objective points, anti_ideal)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis import coverage_metrics  # noqa: E402
from src.paths import stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "baseline_comparison"


def load_kirsch(path: Path) -> np.ndarray:
    """Load Kirsch objectives (n, 2)."""
    path = Path(path)
    if path.suffix == ".npz":
        data = np.load(path, allow_pickle=True)
        if "objectives" in data.files:
            return data["objectives"]
        # Stage 03 build_library writes (all_keys, all_values); pull the
        # first two columns as the headline 2D drought-space.
        return data["all_values"][:, :2]
    chars = json.loads(path.read_text())
    return np.array([[c["mean_duration"], c["mean_avg_severity"]] for c in chars])


def load_moea(path: Path):
    """Load MOEA-FIND results. Returns (objectives (n,2), result_dict)."""
    r = json.loads(Path(path).read_text())
    dm = r.get("drought_metrics", [])
    if not dm:
        return None, r
    return np.array(dm)[:, :2], r


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--kirsch-library", type=Path, required=True,
        help="Path to stage 03 Kirsch library characteristics.npz.",
    )
    p.add_argument(
        "--moea-front", type=Path, required=True,
        help="Path to stage 04 run_moea_find results.json for the production slug.",
    )
    p.add_argument(
        "--convergence", type=Path, default=None,
        help="Optional convergence.json from a stage 03 convergence diagnostic.",
    )
    args = p.parse_args()

    for label, path in [("--kirsch-library", args.kirsch_library),
                        ("--moea-front", args.moea_front)]:
        if not path.exists():
            raise SystemExit(
                f"[04/baseline_comparison] {label}={path} does not exist."
            )

    out = stage_output_dir(STAGE, DRIVER)
    (out / "config.json").write_text(json.dumps({
        "stage": STAGE,
        "driver": DRIVER,
        "kirsch_library": str(args.kirsch_library),
        "moea_front": str(args.moea_front),
        "convergence": str(args.convergence) if args.convergence else None,
    }, indent=2))

    # --- Load data ---
    kirsch = load_kirsch(args.kirsch_library)
    moea, moea_result = load_moea(args.moea_front)

    print(f"[04/baseline_comparison] Kirsch library: {len(kirsch)} traces")
    if moea is None:
        print(f"[04/baseline_comparison] MOEA-FIND: no solutions; aborting.")
        return
    print(f"[04/baseline_comparison] MOEA-FIND: {len(moea)} Pareto solutions")

    # Coverage metrics on the union bounding box
    lb = kirsch.min(axis=0)
    ub = kirsch.max(axis=0)
    kirsch_cm = coverage_metrics(kirsch, lb, ub)
    moea_cm = coverage_metrics(moea, lb, ub)

    print(f"[04/baseline_comparison] Kirsch L2*={kirsch_cm.get('L2_star_discrepancy', 0):.4f}, "
          f"NN_CV={kirsch_cm.get('nn_cv', 0):.3f}")
    print(f"[04/baseline_comparison] MOEA   L2*={moea_cm.get('L2_star_discrepancy', 0):.4f}, "
          f"NN_CV={moea_cm.get('nn_cv', 0):.3f}")

    summary = {
        "kirsch_n": int(len(kirsch)),
        "moea_n": int(len(moea)),
        "kirsch_L2_star": kirsch_cm.get("L2_star_discrepancy"),
        "moea_L2_star": moea_cm.get("L2_star_discrepancy"),
        "kirsch_nn_cv": kirsch_cm.get("nn_cv"),
        "moea_nn_cv": moea_cm.get("nn_cv"),
        "kirsch_range_duration": [float(kirsch[:, 0].min()), float(kirsch[:, 0].max())],
        "kirsch_range_severity": [float(kirsch[:, 1].min()), float(kirsch[:, 1].max())],
        "moea_range_duration": [float(moea[:, 0].min()), float(moea[:, 0].max())],
        "moea_range_severity": [float(moea[:, 1].min()), float(moea[:, 1].max())],
    }
    (out / "comparison_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[04/baseline_comparison] wrote {out / 'comparison_summary.json'}")

    # Persist the pooled point clouds so the plotting driver does not have
    # to re-load the full Kirsch library .npz.
    np.savez(
        out / "pooled.npz",
        kirsch=kirsch,
        moea=moea,
        anti_ideal=np.asarray(moea_result.get("anti_ideal", []), dtype=float),
    )
    print(f"[04/baseline_comparison] wrote {out / 'pooled.npz'}")

    # Optional convergence pass-through: copy the convergence JSON into our
    # output dir so the plot driver has a deterministic local input path.
    if args.convergence and args.convergence.exists():
        conv_data = json.loads(args.convergence.read_text())
        (out / "convergence.json").write_text(json.dumps(conv_data, indent=2))
        print(f"[04/baseline_comparison] copied convergence -> {out / 'convergence.json'}")


if __name__ == "__main__":
    main()
