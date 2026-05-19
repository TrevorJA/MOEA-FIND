"""Pairwise scatter of N random-DV Kirsch-generated years vs the historical
record in a bounded-metric space.

Companion diagnostic to ``bounded_archive_extras.py``: same figure shape
(K×K pairwise grid + diagonal marginals) but with **uniformly random-DV
Kirsch synthetics** in place of the Pareto archive. This isolates the
generator's natural distribution in metric space — i.e., what the
Kirsch-Nowak engine produces *before* the MOEA reshapes it. If the
random-DV cloud already covers the historical envelope tightly, the
MOEA only needs to push outward into the unfilled corners.

Reads the slug's ``config.json`` for ``metric_set``, ``n_years``,
``mode``, and ``ssi`` so the generator settings exactly mirror what
the production run consumed. No MOEA / Pareto data is loaded.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "run_moea_find"


def historical_bounded_points(
    monthly_1d: np.ndarray, objective_keys: list[str]
) -> np.ndarray:
    from src.metrics.extended import BoundedFamilyRefs
    from src.metrics.short_block import compute_candidate_bounded_metrics

    refs = BoundedFamilyRefs.from_full_record(monthly_1d)
    flows = np.asarray(monthly_1d, dtype=float).ravel()
    n_full_yrs = flows.size // 12
    out = np.empty((n_full_yrs, len(objective_keys)), dtype=float)
    for y in range(n_full_yrs):
        block = flows[y * 12 : (y + 1) * 12]
        chars = compute_candidate_bounded_metrics(block, refs)
        for j, k in enumerate(objective_keys):
            out[y, j] = float(chars.get(k, np.nan))
    finite = np.all(np.isfinite(out), axis=1)
    return out[finite]


def kirsch_random_dv_points(
    monthly_2d: np.ndarray,
    monthly_1d: np.ndarray,
    objective_keys: list[str],
    n_samples: int,
    n_years_out: int,
    dv_mode: str,
    seed: int,
) -> np.ndarray:
    from src.metrics.extended import BoundedFamilyRefs
    from src.hydrology.kirsch_utils import build_kirsch_generator
    from src.hydrology.kirsch_wrapper import KirschBorgWrapper
    from src.metrics.short_block import compute_candidate_bounded_metrics

    refs = BoundedFamilyRefs.from_full_record(monthly_1d)
    kirsch = build_kirsch_generator(monthly_2d)
    gen = KirschBorgWrapper(kirsch, mode=dv_mode, n_years_out=n_years_out)
    rng = np.random.default_rng(seed)

    out = np.empty((n_samples, len(objective_keys)), dtype=float)
    for i in range(n_samples):
        dvs = rng.random(gen.n_dvs)
        synthetic_2d = gen.generate(dvs)
        if synthetic_2d.ndim == 1:
            synthetic_2d = synthetic_2d.reshape(n_years_out, 12)
        synthetic_1d = synthetic_2d.flatten()
        chars = compute_candidate_bounded_metrics(synthetic_1d, refs)
        for j, k in enumerate(objective_keys):
            out[i, j] = float(chars.get(k, np.nan))
    finite = np.all(np.isfinite(out), axis=1)
    return out[finite]


def _short(name: str, max_len: int = 14) -> str:
    if len(name) <= max_len:
        return name
    for prefix in ("first_event_", "mean_", "max_", "worst_"):
        if name.startswith(prefix) and len(name) - len(prefix) <= max_len:
            return name[len(prefix):]
    return name[:max_len]


def render(
    hist: np.ndarray,
    kirsch: np.ndarray,
    obj_names: list[str],
    out_path: Path,
    title: str,
) -> None:
    K = len(obj_names)
    fig, axes = plt.subplots(K, K, figsize=(2.6 * K, 2.6 * K))
    short = [_short(n) for n in obj_names]
    for i in range(K):
        for j in range(K):
            ax = axes[i, j]
            if i == j:
                ax.hist(kirsch[:, i], bins=20, color="#16a085", alpha=0.55,
                        density=True, label=f"random-DV Kirsch (n={len(kirsch)})")
                ax.hist(hist[:, i], bins=20, color="#d62728", alpha=0.55,
                        density=True, label=f"historical (n={len(hist)})")
                ax.set_title(short[i], fontsize=10)
                ax.set_yticks([])
                if i == 0:
                    ax.legend(fontsize=6, loc="upper left")
            elif i > j:
                ax.scatter(kirsch[:, j], kirsch[:, i], s=18, alpha=0.7,
                           color="#16a085", marker="^",
                           edgecolor="white", linewidth=0.4, zorder=3,
                           label="random-DV Kirsch")
                ax.scatter(hist[:, j], hist[:, i], s=42, alpha=0.95,
                           color="#d62728", marker="o",
                           edgecolor="white", linewidth=0.6, zorder=4,
                           label="historical years")
                if i == K - 1:
                    ax.set_xlabel(short[j], fontsize=9)
                if j == 0:
                    ax.set_ylabel(short[i], fontsize=9)
            else:
                ax.axis("off")
    fig.suptitle(title, fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True,
                   help="Slug whose config.json drives metric_set, n_years, mode.")
    p.add_argument("--n-samples", type=int, default=100,
                   help="Number of random-DV Kirsch years to generate.")
    p.add_argument("--seed", type=int, default=20260506)
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    cfg = json.loads((in_dir / "config.json").read_text())
    metric_set_name = cfg["metric_set"]
    n_years = int(cfg["n_years"])
    dv_mode = cfg["mode"]

    from src.metrics.drought_metrics import metric_names, resolve_metric_set
    metric_set = resolve_metric_set(metric_set_name)
    obj_names = metric_names(metric_set)

    print(f"[kirsch_history_pairwise] slug={args.slug}")
    print(f"  preset={metric_set_name}  n_years={n_years}  mode={dv_mode}")
    print(f"  metrics: {obj_names}")
    print(f"  n_samples={args.n_samples}  seed={args.seed}")

    from src.experiment import prepare_data
    cache_dir = PROJECT_ROOT / "data" / "usgs_cache"
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    hist = historical_bounded_points(monthly_1d, obj_names)
    print(f"  historical points: n={len(hist)}")
    kirsch = kirsch_random_dv_points(
        monthly_2d, monthly_1d, obj_names,
        n_samples=args.n_samples,
        n_years_out=n_years,
        dv_mode=dv_mode,
        seed=args.seed,
    )
    print(f"  random-DV Kirsch points: n={len(kirsch)}")

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    out_png = fig_dir / "kirsch_history_pairwise.png"
    title = (
        f"Random-DV Kirsch ({args.n_samples}× T={n_years} {dv_mode}) vs historical — "
        f"{metric_set_name} preset, slug {args.slug}"
    )
    render(hist, kirsch, obj_names, out_png, title)
    print(f"  wrote {out_png}")
    print(f"  wrote {out_png.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
