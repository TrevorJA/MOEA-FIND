"""Script 06 — Library subsample baseline (manuscript §6.3, Fig 7).

Consumes the library built by script 05 and produces LHS / Sobol / random
subsamples in drought-characteristic space via nearest-neighbor matching.
Outputs the coverage-metric comparison used in Figure 7.

Run:
    python scripts/06_library_subsample_baseline.py \
        --library outputs/exp05_kirsch_library/characteristics.json \
        --n-select 200 --method lhs
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis import (  # noqa: E402
    coverage_metrics,
    generate_lhs_samples,
    generate_sobol_samples,
)

OUTPUT_SLUG = "exp06_library_subsample"
FEATURE_KEYS = ("mean_duration", "mean_avg_severity")


def load_features(library_path: Path) -> tuple[np.ndarray, list]:
    records = json.loads(library_path.read_text())
    feats = np.array([[r[k] for k in FEATURE_KEYS] for r in records], dtype=float)
    return feats, records


def nearest_neighbor_select(design: np.ndarray, library: np.ndarray) -> np.ndarray:
    """For each design point, pick the library index closest in feature space."""
    from scipy.spatial import cKDTree
    tree = cKDTree(library)
    _, idx = tree.query(design, k=1)
    return np.unique(idx)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--library", type=Path, required=True,
                   help="Path to characteristics.json produced by script 05")
    p.add_argument("--n-select", type=int, default=200)
    p.add_argument("--method", choices=["lhs", "sobol", "random"], default="lhs")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "config.json").write_text(json.dumps({
        "script": "06_library_subsample_baseline.py",
        "manuscript_section": "§6.3 Library Baseline (Fig 7)",
        "library": str(args.library),
        "n_select": args.n_select, "method": args.method, "seed": args.seed,
    }, indent=2))

    feats, records = load_features(args.library)
    print(f"[06] library: {len(feats)} traces, features={FEATURE_KEYS}")

    lb = feats.min(axis=0)
    ub = feats.max(axis=0)
    d = feats.shape[1]

    if args.method == "lhs":
        design = generate_lhs_samples(args.n_select, d, lb, ub, seed=args.seed)
    elif args.method == "sobol":
        design = generate_sobol_samples(args.n_select, d, lb, ub, seed=args.seed)[:args.n_select]
    else:
        rng = np.random.default_rng(args.seed)
        idx = rng.choice(len(feats), args.n_select, replace=False)
        design = feats[idx]

    selected_idx = nearest_neighbor_select(design, feats)
    selected_feats = feats[selected_idx]
    print(f"[06] {args.method}: {len(selected_idx)} unique library traces selected")

    metrics_subset = coverage_metrics(selected_feats, lb, ub)
    metrics_full = coverage_metrics(feats, lb, ub)

    summary = {
        "method": args.method,
        "n_library": int(len(feats)),
        "n_selected": int(len(selected_idx)),
        "feature_keys": list(FEATURE_KEYS),
        "feature_bounds": {"lb": lb.tolist(), "ub": ub.tolist()},
        "subset_coverage": metrics_subset,
        "full_library_coverage": metrics_full,
        "selected_trace_ids": [int(records[i]["trace_id"]) for i in selected_idx],
    }
    out = args.output_dir / f"subsample_{args.method}.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"[06] wrote {out}")


if __name__ == "__main__":
    main()
