"""Space-filling comparison: historical / random-DV Kirsch / Pareto in
the bounded metric hypercube ``[0, 1]^K``.

Three point clouds, all interpreted as samples in the same hypercube:

  1. **Historical** — n_hist = 73, one point per water year via
     ``compute_candidate_bounded_metrics`` on the historical record.
  2. **Random-DV Kirsch** — n_kirsch = 10,000 by default, generated with
     uniform DVs in [0,1]^n_dvs through the same KirschBorgWrapper /
     mode the production run consumed.
  3. **Pareto** — ``results.json`` ``drought_metrics`` columns matching
     ``objective_keys`` (Pareto archive of the production MOEA run).

Per source we report two complementary kinds of metric:

  * **Discrepancies** (``scipy.stats.qmc.discrepancy``) — L2-star,
    wrap-around (WD), centered (CD), and mixture (MD). Lower = more
    uniform space-filling.

  * **Geometric coverage** — fill distance ``h_X`` (covering radius
    against a Halton 4D grid) and minimum pairwise distance ``d_min``.
    These pick up on gaps and clumping that integrated discrepancies
    can miss.

To make the three sources comparable despite very different sample
sizes (73 vs 10,000 vs ~95,000), every per-source value is reported in
two forms:

  * **Full-set** — computed once on the full point cloud.
  * **Subsampled-to-73 bootstrap** — n_boot replicates of size 73 drawn
    with replacement; the median and 95% interval form a fair side-by-
    side comparison anchored at the historical sample size.
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
from scipy.spatial import cKDTree
from scipy.stats import qmc

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "run_moea_find"

DISCREPANCY_METHODS = ("L2-star", "WD", "CD", "MD")
DISCREPANCY_LABELS = {
    "L2-star": "L²-star",
    "WD":      "wrap-around",
    "CD":      "centered",
    "MD":      "mixture",
}


def historical_bounded_points(monthly_1d, objective_keys: list[str]) -> np.ndarray:
    from src.metrics.extended import BoundedFamilyRefs
    from src.metrics.short_block import compute_candidate_bounded_metrics
    refs = BoundedFamilyRefs.from_full_record(monthly_1d)
    flows = np.asarray(monthly_1d, dtype=float).ravel()
    n = flows.size // 12
    out = np.empty((n, len(objective_keys)), dtype=float)
    for y in range(n):
        chars = compute_candidate_bounded_metrics(flows[y * 12 : (y + 1) * 12], refs)
        for j, k in enumerate(objective_keys):
            out[y, j] = float(chars.get(k, np.nan))
    return out[np.all(np.isfinite(out), axis=1)]


def kirsch_random_dv_points(
    monthly_2d, monthly_1d, objective_keys, n_samples, n_years_out, dv_mode, seed
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
        synthetic_2d = gen.generate(rng.random(gen.n_dvs))
        if synthetic_2d.ndim == 1:
            synthetic_2d = synthetic_2d.reshape(n_years_out, 12)
        chars = compute_candidate_bounded_metrics(synthetic_2d.flatten(), refs)
        for j, k in enumerate(objective_keys):
            out[i, j] = float(chars.get(k, np.nan))
    return out[np.all(np.isfinite(out), axis=1)]


def pareto_points(result: dict, objective_keys: list[str]) -> np.ndarray:
    dm = np.asarray(result["drought_metrics"], dtype=float)
    keys = list(result["objective_keys"])
    cols = [keys.index(k) for k in objective_keys]
    pts = dm[:, cols]
    return pts[np.all(np.isfinite(pts), axis=1)]


def _safe_clip(pts: np.ndarray) -> np.ndarray:
    """Clip to (0,1) by a hair so qmc.discrepancy doesn't reject the input.
    Mapping G values are in [0, 1), but tiny FP / saturation can overshoot."""
    eps = 1e-9
    return np.clip(pts, eps, 1.0 - eps)


def discrepancies(pts: np.ndarray) -> dict[str, float]:
    p = _safe_clip(pts)
    out = {}
    for m in DISCREPANCY_METHODS:
        try:
            out[m] = float(qmc.discrepancy(p, method=m))
        except Exception as e:
            out[m] = float("nan")
    return out


def fill_distance(pts: np.ndarray, n_grid: int, dim: int) -> float:
    """Halton-grid covering radius: max nearest-sample distance over a
    quasi-random 4D grid of size n_grid."""
    halton = qmc.Halton(d=dim, scramble=False).random(n_grid)
    tree = cKDTree(pts)
    d, _ = tree.query(halton, k=1)
    return float(d.max())


def min_pairwise_distance(pts: np.ndarray) -> float:
    if len(pts) < 2:
        return float("nan")
    tree = cKDTree(pts)
    d, _ = tree.query(pts, k=2)
    return float(d[:, 1].min())


def bootstrap_metrics(
    pts: np.ndarray, n_match: int, n_boot: int, seed: int, n_grid_fd: int, dim: int
) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    n = len(pts)
    keys = list(DISCREPANCY_METHODS) + ["fill_distance", "min_pairwise"]
    out = {k: np.empty(n_boot) for k in keys}
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n_match)
        sub = pts[idx]
        d = discrepancies(sub)
        for m in DISCREPANCY_METHODS:
            out[m][b] = d[m]
        out["fill_distance"][b] = fill_distance(sub, n_grid_fd, dim)
        out["min_pairwise"][b] = min_pairwise_distance(sub)
    return out


def render(
    sources: dict[str, dict],
    out_path: Path,
    title: str,
    n_match: int,
) -> None:
    metrics = list(DISCREPANCY_METHODS) + ["fill_distance", "min_pairwise"]
    metric_titles = [
        f"L²-star\n(lower = more uniform)",
        f"wrap-around\n(lower = more uniform)",
        f"centered\n(lower = more uniform)",
        f"mixture\n(lower = more uniform)",
        f"fill distance\n(lower = fewer gaps)",
        f"min pairwise dist\n(higher = less clumping)",
    ]
    n_panels = len(metrics)
    n_cols = 3
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.6 * n_cols, 3.4 * n_rows))
    axes = axes.flatten()

    src_names = list(sources.keys())
    src_colors = {"historical": "#d62728",
                  "random-DV Kirsch": "#16a085",
                  "Pareto": "#2b6cb0"}

    for k, m in enumerate(metrics):
        ax = axes[k]
        # Boxplots of bootstrap distributions
        data = [sources[s]["boot"][m] for s in src_names]
        bp = ax.boxplot(data, labels=src_names, patch_artist=True, widths=0.55,
                        showfliers=False, medianprops={"color": "black", "lw": 1.6})
        for patch, sname in zip(bp["boxes"], src_names):
            patch.set_facecolor(src_colors.get(sname, "#888"))
            patch.set_alpha(0.65)
        # Full-set values as red dashes
        for i, s in enumerate(src_names, start=1):
            full = sources[s]["full"].get(m, np.nan)
            if np.isfinite(full):
                ax.axhline(full, xmin=(i - 1) / len(src_names) + 0.05,
                           xmax=i / len(src_names) - 0.05,
                           color="black", lw=1.4, ls="--")
                ax.scatter([i], [full], s=44, color="black",
                           marker="D", zorder=5,
                           label="full set" if (k == 0 and i == 1) else None)
        ax.set_title(metric_titles[k], fontsize=10)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.25)
        if k == 0:
            ax.legend(loc="upper right", fontsize=7)

    for k in range(len(metrics), len(axes)):
        axes[k].axis("off")

    fig.suptitle(title + f"  (boxes: bootstrap n={n_match}; diamonds: full-set)",
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True)
    p.add_argument("--n-kirsch", type=int, default=10000)
    p.add_argument("--n-boot", type=int, default=500,
                   help="Number of bootstrap replicates per source.")
    p.add_argument("--n-grid-fd", type=int, default=4096,
                   help="Halton grid size for fill-distance coverage radius.")
    p.add_argument("--seed", type=int, default=20260506)
    p.add_argument("--full-set-cap", type=int, default=5000,
                   help="Cap the full-set discrepancy/coverage computation at "
                        "this many uniformly subsampled points per source. "
                        "Bootstrap (matched to historical n) is unaffected.")
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    cfg = json.loads((in_dir / "config.json").read_text())
    result = json.loads((in_dir / "results.json").read_text())
    metric_set_name = cfg["metric_set"]
    n_years = int(cfg["n_years"])
    dv_mode = cfg["mode"]

    from src.metrics.drought_metrics import metric_names, resolve_metric_set
    metric_set = resolve_metric_set(metric_set_name)
    objective_keys = list(metric_names(metric_set))
    K = len(objective_keys)

    print(f"[space_filling_compare] slug={args.slug}")
    print(f"  preset={metric_set_name}  K={K}  n_years={n_years}  mode={dv_mode}")

    from src.experiment import prepare_data
    cache_dir = PROJECT_ROOT / "data" / "usgs_cache"
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    hist_pts = historical_bounded_points(monthly_1d, objective_keys)
    print(f"  historical points: n={len(hist_pts)}")
    kirsch_pts = kirsch_random_dv_points(
        monthly_2d, monthly_1d, objective_keys,
        n_samples=args.n_kirsch, n_years_out=n_years,
        dv_mode=dv_mode, seed=args.seed,
    )
    print(f"  random-DV Kirsch points: n={len(kirsch_pts)}")
    pareto_pts = pareto_points(result, objective_keys)
    print(f"  Pareto points: n={len(pareto_pts)}")

    n_match = int(len(hist_pts))  # bootstrap each source down to historical n
    # qmc.discrepancy is O(n^2 d) (O(n^2 d^2) for MD) — quadratic blow-up on
    # the n≈95k Pareto archive dominates wall time. The full-set values are
    # only used as cosmetic anchors on the comparison plot, so cap their
    # input at a uniform random subsample of `full_set_cap` points per source.
    rng_full = np.random.default_rng(args.seed + 7)
    sources = {}
    for name, pts in [
        ("historical", hist_pts),
        ("random-DV Kirsch", kirsch_pts),
        ("Pareto", pareto_pts),
    ]:
        if len(pts) > args.full_set_cap:
            sub_idx = rng_full.choice(len(pts), size=args.full_set_cap, replace=False)
            pts_full = pts[sub_idx]
            full_n = args.full_set_cap
        else:
            pts_full = pts
            full_n = len(pts)
        full = discrepancies(pts_full)
        full["fill_distance"] = fill_distance(pts_full, args.n_grid_fd, K)
        full["min_pairwise"] = min_pairwise_distance(pts_full)
        boot = bootstrap_metrics(
            pts, n_match=n_match, n_boot=args.n_boot,
            seed=args.seed + hash(name) % 10_000,
            n_grid_fd=args.n_grid_fd, dim=K,
        )
        sources[name] = {"n": int(len(pts)), "n_full_set": int(full_n),
                         "full": full, "boot": boot}
        print(f"  {name:18s} full(n={full_n}) L²-star={full['L2-star']:.4f}  "
              f"fill_dist={full['fill_distance']:.3f}  "
              f"min_pair={full['min_pairwise']:.4f}")

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    out_png = fig_dir / "space_filling_compare.png"
    title = (
        f"Space-filling comparison — {metric_set_name} preset, slug {args.slug}\n"
        f"K={K} hypercube  |  hist n={len(hist_pts)}, "
        f"Kirsch n={len(kirsch_pts)}, Pareto n={len(pareto_pts)}"
    )
    render(sources, out_png, title, n_match=n_match)
    print(f"  wrote {out_png}")
    print(f"  wrote {out_png.with_suffix('.pdf')}")

    # Summary JSON
    summary = {
        "slug": args.slug,
        "metric_set": metric_set_name,
        "K": K,
        "n_match_bootstrap": n_match,
        "n_boot": args.n_boot,
        "full_set_cap": args.full_set_cap,
        "sources": {
            name: {
                "n": s["n"],
                "n_full_set": s["n_full_set"],
                "full": s["full"],
                "bootstrap": {
                    m: {
                        "median": float(np.median(s["boot"][m])),
                        "q025": float(np.quantile(s["boot"][m], 0.025)),
                        "q975": float(np.quantile(s["boot"][m], 0.975)),
                        "mean": float(np.mean(s["boot"][m])),
                        "std": float(np.std(s["boot"][m])),
                    }
                    for m in list(DISCREPANCY_METHODS) + ["fill_distance", "min_pairwise"]
                },
            }
            for name, s in sources.items()
        },
    }
    json_path = fig_dir / "space_filling_compare.json"
    json_path.write_text(json.dumps(summary, indent=2))
    print(f"  wrote {json_path}")


if __name__ == "__main__":
    main()
