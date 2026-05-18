"""Extra diagnostics for any bounded-metric MOEA-FIND archive (e.g.
``short_block_drb_v2``). Two outputs:

  1. 2D pairwise Pareto scatter with the **historical 1-year envelope**
     overlaid in the same bounded metric space. The historical points
     are computed by running the registered preset's bounded-metric
     extractor on every historical water year, using the same
     ``BoundedFamilyRefs`` the production run consumed.

  2. 3D scatter montage covering the first 4 distinct 3-axis projections
     of the preset (``C(K, 3)`` capped at 4 panels). Anti-ideal corner
     marked. Historical overlay shown in the same projections.

Generic across any preset whose objective keys are produced by
``compute_candidate_bounded_metrics`` — i.e., any ``<window>_g`` /
``<window>_e`` family in :data:`src.metrics.short_block.WINDOW_SPECS`.
Reads the run's ``results.json`` / ``pareto.npz`` directly so it works
without re-running the MOEA.
"""

from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "run_moea_find"


def historical_bounded_points(
    monthly_1d: np.ndarray, objective_keys: list[str]
) -> np.ndarray:
    """Bounded-metric vector for every historical water year.

    Returns shape ``(n_hist_years, K)`` aligned with ``objective_keys``.
    Skips any historical year for which an extractor returns NaN.
    """
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


def _short(name: str, max_len: int = 14) -> str:
    if len(name) <= max_len:
        return name
    for prefix in ("first_event_", "mean_", "max_", "worst_"):
        if name.startswith(prefix) and len(name) - len(prefix) <= max_len:
            return name[len(prefix):]
    return name[:max_len]


def render_pairwise_with_history(
    pareto: np.ndarray,
    hist: np.ndarray,
    anti_ideal: np.ndarray,
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
                ax.hist(pareto[:, i], bins=60, color="#2b6cb0", alpha=0.8,
                        density=True, label=f"Pareto (n={len(pareto)})")
                ax.hist(hist[:, i], bins=20, color="#d62728", alpha=0.55,
                        density=True, label=f"historical (n={len(hist)})")
                ax.axvline(anti_ideal[i], color="black", linewidth=1.2,
                           linestyle="--", label=f"D*={anti_ideal[i]:.2g}")
                ax.set_title(short[i], fontsize=10)
                ax.set_yticks([])
                if i == 0:
                    ax.legend(fontsize=6, loc="upper left")
            elif i > j:
                ax.scatter(pareto[:, j], pareto[:, i], s=2, alpha=0.25,
                           color="#2b6cb0", label="Pareto")
                ax.scatter(hist[:, j], hist[:, i], s=42, alpha=0.95,
                           color="#d62728", marker="o",
                           edgecolor="white", linewidth=0.6, zorder=4,
                           label="historical years")
                ax.scatter(anti_ideal[j], anti_ideal[i],
                           marker="X", s=80, color="black", zorder=5,
                           edgecolor="white", linewidth=0.6, label="D*")
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
    plt.close(fig)


def render_3d_montage(
    pareto: np.ndarray,
    hist: np.ndarray,
    anti_ideal: np.ndarray,
    obj_names: list[str],
    out_path: Path,
    title: str,
    max_panels: int = 4,
) -> None:
    K = len(obj_names)
    if K < 3:
        return
    triplets = list(combinations(range(K), 3))[:max_panels]
    n = len(triplets)
    short = [_short(n_, 14) for n_ in obj_names]
    n_cols = 2 if n > 1 else 1
    n_rows = (n + n_cols - 1) // n_cols
    fig = plt.figure(figsize=(6.2 * n_cols, 5.4 * n_rows))
    # Subsample Pareto for legibility (3D is expensive to render at 95k pts)
    rng = np.random.default_rng(0)
    if len(pareto) > 8000:
        idx = rng.choice(len(pareto), size=8000, replace=False)
        pareto_plot = pareto[idx]
    else:
        pareto_plot = pareto
    for k, (a, b, c) in enumerate(triplets):
        ax = fig.add_subplot(n_rows, n_cols, k + 1, projection="3d")
        ax.scatter(pareto_plot[:, a], pareto_plot[:, b], pareto_plot[:, c],
                   s=2, alpha=0.25, color="#2b6cb0",
                   label=f"Pareto subsample (n={len(pareto_plot)})")
        ax.scatter(hist[:, a], hist[:, b], hist[:, c],
                   s=44, alpha=0.95, color="#d62728", marker="o",
                   edgecolor="white", linewidth=0.7,
                   label=f"historical years (n={len(hist)})")
        ax.scatter([anti_ideal[a]], [anti_ideal[b]], [anti_ideal[c]],
                   marker="X", s=110, color="black",
                   edgecolor="white", linewidth=0.6, label="D*")
        ax.set_xlabel(short[a], fontsize=8, labelpad=2)
        ax.set_ylabel(short[b], fontsize=8, labelpad=2)
        ax.set_zlabel(short[c], fontsize=8, labelpad=2)
        ax.tick_params(labelsize=6)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
        ax.set_title(f"{short[a]} × {short[b]} × {short[c]}", fontsize=9)
        if k == 0:
            ax.legend(fontsize=7, loc="upper left")
    fig.suptitle(title, fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True)
    p.add_argument("--max-3d-panels", type=int, default=4,
                   help="Cap the 3D montage at this many distinct triplets.")
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    result = json.loads((in_dir / "results.json").read_text())
    pareto = np.array(result["drought_metrics"], dtype=float)
    anti_ideal = np.asarray(result["anti_ideal"], dtype=float)
    obj_names = list(result["objective_keys"])
    print(f"[bounded_extras] slug={args.slug}")
    print(f"  K={len(obj_names)} objectives: {obj_names}")
    print(f"  Pareto n={len(pareto)} anti_ideal={anti_ideal}")

    from src.experiment import prepare_data
    cache_dir = PROJECT_ROOT / "data" / "usgs_cache"
    _, monthly_1d = prepare_data(cache_dir)
    hist = historical_bounded_points(monthly_1d, obj_names)
    print(f"  historical 1-year-block points: n={len(hist)}")

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    cfg = json.loads((in_dir / "config.json").read_text())
    base_title = (
        f"{cfg.get('metric_set', '?')} archive — {args.slug}\n"
        f"K={len(obj_names)} T={cfg.get('n_years','?')} "
        f"NFE={cfg['nfe']:,} seed={cfg['seed']} "
        f"n_pareto={result['n_pareto']}"
    )
    render_pairwise_with_history(
        pareto, hist, anti_ideal, obj_names,
        fig_dir / "bounded_pairwise_with_history.png",
        title=base_title + "\nPareto + historical 1-year blocks",
    )
    print(f"  wrote {fig_dir / 'bounded_pairwise_with_history.png'}")

    render_3d_montage(
        pareto, hist, anti_ideal, obj_names,
        fig_dir / "bounded_3d_montage.png",
        title=base_title + "\n3D projections (subsample)",
        max_panels=args.max_3d_panels,
    )
    print(f"  wrote {fig_dir / 'bounded_3d_montage.png'}")


if __name__ == "__main__":
    main()
