"""Library subsample baseline (manuscript Fig 7).

Consumes the library built by ``build_library.py`` and produces
LHS / Sobol / random subsamples in drought-characteristic space via
nearest-neighbor matching. Outputs the per-method coverage summaries
consumed by stage 04's ``baseline_comparison.py`` driver, which
renders Fig 7.

Outputs under ``outputs/03_kirsch_library/subsample_baseline/<method>/<slug>/``:
    - config.json
    - subsample_<method>.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.discovery.analysis import (  # noqa: E402
    coverage_metrics,
    generate_lhs_samples,
    generate_sobol_samples,
)
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.io_paths.slugs import subsample_slug  # noqa: E402

STAGE = "03_kirsch_library"
DRIVER = "subsample_baseline"
FEATURE_KEYS = ("mean_duration", "mean_avg_severity")


def _src_id(library: Path) -> str:
    """Compact id of the upstream library run, safe to embed in another slug.

    ``library.parent.name`` is the library's own slug
    (``library__N=10000__T=20__s=42``); strip the ``__`` / ``=``
    separators so the resulting subsample slug round-trips through
    :func:`src.io_paths.slugs.parse_slug` unambiguously.
    """
    return library.parent.name.replace("__", "_").replace("=", "")


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
                   help="Path to characteristics.json produced by build_library.py")
    p.add_argument("--n-select", type=int, required=True)
    p.add_argument("--method", choices=["lhs", "sobol", "random"], required=True)
    p.add_argument("--seed", type=int, required=True)
    args = p.parse_args()

    method_driver = f"{DRIVER}/{args.method}"
    slug = subsample_slug(src=_src_id(args.library), method=args.method,
                          n_select=args.n_select, seed=args.seed)
    out_dir = stage_output_dir(STAGE, method_driver, slug)
    (out_dir / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "library": str(args.library),
        "n_select": args.n_select,
        "method": args.method,
        "seed": args.seed,
    }, indent=2))

    feats, records = load_features(args.library)
    print(f"[03/subsample_baseline] library: {len(feats)} traces, features={FEATURE_KEYS}")

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
    print(f"[03/subsample_baseline] {args.method}: "
          f"{len(selected_idx)} unique library traces selected")

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
    out = out_dir / f"subsample_{args.method}.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"[03/subsample_baseline] wrote {out}")


if __name__ == "__main__":
    main()
