"""SI bootstrap-uncertainty diagnostic for the Mapping G transformation.

Mapping G uses ``μ̂``, ``σ̂`` estimated from n=73 historical water years.
Sampling variability in those estimates propagates into ``D_g``: a
candidate's bounded objective is *not* a fixed function of its
hydrology — it depends on which 73 years happened to be observed.

This script quantifies that propagation by **resampling the historical
record with replacement** (n_boot replicates, default 2000), refitting
``(μ̂*, σ̂*)`` per replicate, and evaluating ``D_g`` on a dense grid of
``s_j`` values. For each grid point we record the bootstrap distribution
and report:

  - 95% / 50% bootstrap interval around the point estimate
  - Standard error ``SE(D_g)``
  - Implied "percentile uncertainty" — i.e., how many percentile points
    the metric could shift if the same 73-year record were redrawn.

Two output panels per window:

  * **Band view**: ``D_g(s_j)`` point estimate with 50% and 95% bootstrap
    bands; empirical CDF underlay for context.
  * **SE profile**: ``SE(D_g)`` as a function of the point-estimate
    ``D_g`` — quantifies which percentile bands are most stable under
    historical resampling.

Generic across any bounded preset whose objective keys are
``<window>_g``. Reads only ``results.json`` for window discovery; does
not require re-running the MOEA.
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
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "run_moea_find"


def _window_name_from_metric(metric_key: str) -> str | None:
    for suffix in ("_g", "_e"):
        if metric_key.endswith(suffix):
            return metric_key[: -len(suffix)]
    return None


def historical_window_samples(monthly_1d: np.ndarray, window_name: str) -> np.ndarray:
    from src.metrics.short_block import WINDOW_SPECS, _compute_window_summary

    if window_name not in WINDOW_SPECS:
        raise KeyError(f"unknown window {window_name!r}")
    kind, params = WINDOW_SPECS[window_name]
    flows = np.asarray(monthly_1d, dtype=float).ravel()
    n = flows.size // 12
    flows_2d = flows[: n * 12].reshape(n, 12)
    out = np.array([_compute_window_summary(flows_2d[y : y + 1], kind, params)
                    for y in range(n)])
    return out[np.isfinite(out)]


def bootstrap_dg(
    samples: np.ndarray,
    s_grid: np.ndarray,
    n_boot: int,
    seed: int,
) -> np.ndarray:
    """Bootstrap distribution of D_g(s_j) at each s_j on s_grid.

    Returns shape ``(n_boot, len(s_grid))``. Each row is one bootstrap
    replicate's ``D_g`` curve. Samples drawn with replacement from
    ``samples``; ``μ̂*``, ``σ̂*`` recomputed; ``D_g* = Φ((s − μ̂*)/σ̂*)``.
    """
    rng = np.random.default_rng(seed)
    n = samples.size
    n_grid = s_grid.size
    out = np.empty((n_boot, n_grid), dtype=float)
    # Vectorised bootstrap: draw all indices at once.
    idx = rng.integers(0, n, size=(n_boot, n))
    resamples = samples[idx]                    # shape (n_boot, n)
    mu = resamples.mean(axis=1)                 # shape (n_boot,)
    sigma = resamples.std(axis=1, ddof=1)       # shape (n_boot,)
    sigma = np.where(sigma > 1e-12, sigma, np.nan)
    z = (s_grid[None, :] - mu[:, None]) / sigma[:, None]   # (n_boot, n_grid)
    out = stats.norm.cdf(z)
    return out


def render_bootstrap_grid(
    window_data: list[dict],
    out_path: Path,
    title: str,
) -> None:
    n_w = len(window_data)
    fig, axes = plt.subplots(2, n_w, figsize=(3.6 * n_w, 6.4))
    if n_w == 1:
        axes = axes[:, None]

    for j, w in enumerate(window_data):
        s_grid = w["s_grid"]
        boot = w["boot"]
        d_point = w["d_point"]
        samples = w["samples"]
        wname = w["window"]

        # Bootstrap quantile bands
        q025 = np.nanpercentile(boot, 2.5, axis=0)
        q25 = np.nanpercentile(boot, 25.0, axis=0)
        q75 = np.nanpercentile(boot, 75.0, axis=0)
        q975 = np.nanpercentile(boot, 97.5, axis=0)
        se = np.nanstd(boot, axis=0)

        # Row 0: band view of D_g(s_j)
        ax = axes[0, j]
        ax.fill_between(s_grid, q025, q975, color="#9fb6d6", alpha=0.55,
                        label="95% bootstrap")
        ax.fill_between(s_grid, q25, q75, color="#5b8ac0", alpha=0.7,
                        label="50% bootstrap")
        ax.plot(s_grid, d_point, color="#1f3a5f", lw=1.6,
                label=r"point estimate $D_g$")
        # empirical CDF for visual anchor
        s_sorted = np.sort(samples)
        cdf_emp = np.arange(1, samples.size + 1) / samples.size
        ax.plot(s_sorted, cdf_emp, drawstyle="steps-post",
                color="#d62728", lw=1.0, alpha=0.85,
                label="empirical CDF")
        ax.set_xlabel(r"$s_j$ (drought-positive)", fontsize=8)
        if j == 0:
            ax.set_ylabel(r"$D_g$", fontsize=9)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(f"{wname}\nmedian SE={np.nanmedian(se):.3f}, max SE={np.nanmax(se):.3f}",
                     fontsize=9)
        if j == 0:
            ax.legend(fontsize=7, loc="lower right")
        ax.tick_params(labelsize=7)

        # Row 1: SE(D_g) vs D_g — uncertainty profile
        ax = axes[1, j]
        order = np.argsort(d_point)
        ax.plot(d_point[order], se[order], color="#16a085", lw=1.6)
        ax.axhline(np.nanmedian(se), color="0.5", lw=0.8, ls=":",
                   label=f"median={np.nanmedian(se):.3f}")
        # Annotate worst-case
        i_max = int(np.nanargmax(se))
        ax.scatter([d_point[i_max]], [se[i_max]], s=30, color="#d62728",
                   zorder=4, label=f"max={se[i_max]:.3f} @ D_g={d_point[i_max]:.2f}")
        ax.set_xlabel(r"$D_g$ point estimate", fontsize=8)
        if j == 0:
            ax.set_ylabel(r"$SE(D_g)$", fontsize=9)
        ax.set_xlim(0, 1)
        ax.legend(fontsize=7, loc="upper right")
        ax.set_title("uncertainty profile", fontsize=9)
        ax.tick_params(labelsize=7)

        # Stash summary
        # Probe percentile uncertainty at canonical D_g points (10/50/90).
        probe_targets = np.array([0.10, 0.50, 0.90])
        # Find s_grid points whose point-estimate D_g is closest to each target.
        probe_se = []
        for t in probe_targets:
            i = int(np.argmin(np.abs(d_point - t)))
            probe_se.append(float(se[i]))
        w["summary"] = {
            "n_boot": int(boot.shape[0]),
            "median_se": float(np.nanmedian(se)),
            "max_se": float(np.nanmax(se)),
            "max_se_at_dg": float(d_point[i_max]),
            "se_at_dg_0.10": probe_se[0],
            "se_at_dg_0.50": probe_se[1],
            "se_at_dg_0.90": probe_se[2],
        }

    fig.suptitle(title, fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True)
    p.add_argument("--n-boot", type=int, default=2000,
                   help="Number of bootstrap replicates per window.")
    p.add_argument("--n-grid", type=int, default=200,
                   help="Number of grid points across the historical envelope.")
    p.add_argument("--grid-extrap", type=float, default=0.5,
                   help="Extra width as a fraction of σ̂ on each side of "
                        "the historical [s_min, s_max] envelope.")
    p.add_argument("--seed", type=int, default=20260506)
    p.add_argument("--windows", nargs="*", default=None)
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    result = json.loads((in_dir / "results.json").read_text())
    obj_keys = list(result["objective_keys"])

    if args.windows:
        windows = list(args.windows)
    else:
        windows = []
        for k in obj_keys:
            w = _window_name_from_metric(k)
            if w is not None and w not in windows:
                windows.append(w)
        if not windows:
            sys.exit(f"no _g/_e suffixed metrics in objective_keys: {obj_keys}")

    print(f"[mapping_g_bootstrap] slug={args.slug}")
    print(f"  windows: {windows}  n_boot={args.n_boot}  n_grid={args.n_grid}")

    from src.experiment import prepare_data
    cache_dir = PROJECT_ROOT / "data" / "usgs_cache"
    _, monthly_1d = prepare_data(cache_dir)

    window_data = []
    for w_name in windows:
        samples = historical_window_samples(monthly_1d, w_name)
        mu = float(samples.mean())
        sigma = float(samples.std(ddof=1))
        # Grid covering historical envelope plus modest extrapolation
        pad = args.grid_extrap * sigma
        s_grid = np.linspace(samples.min() - pad, samples.max() + pad, args.n_grid)
        d_point = stats.norm.cdf((s_grid - mu) / sigma)
        boot = bootstrap_dg(samples, s_grid, args.n_boot, args.seed)
        window_data.append({
            "window": w_name,
            "samples": samples,
            "mu": mu,
            "sigma": sigma,
            "s_grid": s_grid,
            "d_point": d_point,
            "boot": boot,
        })
        med_se = float(np.nanmedian(np.nanstd(boot, axis=0)))
        max_se = float(np.nanmax(np.nanstd(boot, axis=0)))
        print(f"  {w_name:24s} n={samples.size:3d}  μ={mu:.3f}  σ={sigma:.3f}  "
              f"median SE={med_se:.3f}  max SE={max_se:.3f}")

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    cfg = json.loads((in_dir / "config.json").read_text())
    title = (
        f"Mapping G bootstrap uncertainty (n_boot={args.n_boot}) — "
        f"{cfg.get('metric_set','?')} preset, slug {args.slug}\n"
        f"Top: 50%/95% bootstrap bands on D_g(s_j).  "
        f"Bottom: SE(D_g) profile."
    )
    out_png = fig_dir / "mapping_g_bootstrap.png"
    render_bootstrap_grid(window_data, out_png, title)
    print(f"  wrote {out_png}")
    print(f"  wrote {out_png.with_suffix('.pdf')}")

    # JSON for the SI table
    summary = {
        "slug": args.slug,
        "metric_set": cfg.get("metric_set"),
        "n_years": cfg.get("n_years"),
        "n_boot": args.n_boot,
        "n_grid": args.n_grid,
        "windows": [
            {
                "window": w["window"],
                "n": int(w["samples"].size),
                "mu": w["mu"],
                "sigma": w["sigma"],
                **w["summary"],
            }
            for w in window_data
        ],
    }
    json_path = fig_dir / "mapping_g_bootstrap.json"
    json_path.write_text(json.dumps(summary, indent=2))
    print(f"  wrote {json_path}")


if __name__ == "__main__":
    main()
