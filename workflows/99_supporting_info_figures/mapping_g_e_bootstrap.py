"""SI bootstrap-uncertainty comparison: Mapping G vs Mapping E.

Mirrors ``mapping_g_bootstrap.py`` but evaluates **both** mappings on
each bootstrap replicate, so that body-vs-tail SE behavior is directly
comparable. Per replicate:

  1. Resample n=73 historical samples with replacement.
  2. Refit ``μ̂*``, ``σ̂*`` for Mapping G; ``sorted_hist*`` for Mapping E.
  3. Evaluate ``D_g*(s)`` and ``D_e*(s)`` on the same ``s_grid``.

Two-row figure per window:

  * **Top — band view**: 95% bootstrap bands for G (blue) and E
    (orange) overlaid on a single axis with the point estimate of each.
  * **Bottom — SE profile**: ``SE(D_g)`` and ``SE(D_e)`` as functions
    of the point-estimate ``D``.

For SI: the body should show similar SE for G and E (both are fit to
the same n=73 sample); the **tails** are the discriminator. Mapping E's
polynomial tail is anchored to the top/bottom k=5 samples, so its
tail SE inherits high sampling variance there. Mapping G saturates
toward 0 / 1 in the tails, so its SE collapses by construction.
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

PROJECT_ROOT = Path(__file__).resolve().parents[2]
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


def bootstrap_g_and_e(
    samples: np.ndarray,
    s_grid: np.ndarray,
    n_boot: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Joint bootstrap distribution of D_g(s_j) and D_e(s_j).

    Returns ``(boot_g, boot_e)`` each of shape ``(n_boot, len(s_grid))``.
    """
    from src.metrics.short_block import _apply_mapping_e

    rng = np.random.default_rng(seed)
    n = samples.size
    n_grid = s_grid.size
    boot_g = np.empty((n_boot, n_grid), dtype=float)
    boot_e = np.empty((n_boot, n_grid), dtype=float)
    idx = rng.integers(0, n, size=(n_boot, n))
    for b in range(n_boot):
        resample = samples[idx[b]]
        sorted_re = np.sort(resample)
        mu = float(resample.mean())
        sigma = float(resample.std(ddof=1))
        if sigma <= 1e-12:
            boot_g[b] = np.nan
        else:
            boot_g[b] = stats.norm.cdf((s_grid - mu) / sigma)
        for i, s in enumerate(s_grid):
            boot_e[b, i] = _apply_mapping_e(float(s), sorted_re)
    return boot_g, boot_e


def render_compare(
    window_data: list[dict], out_path: Path, title: str
) -> None:
    n_w = len(window_data)
    fig, axes = plt.subplots(2, n_w, figsize=(3.6 * n_w, 6.6))
    if n_w == 1:
        axes = axes[:, None]

    for j, w in enumerate(window_data):
        s_grid = w["s_grid"]
        bg, be = w["boot_g"], w["boot_e"]
        dg, de = w["d_point_g"], w["d_point_e"]
        s = w["samples"]
        wname = w["window"]

        gq025 = np.nanpercentile(bg, 2.5, axis=0)
        gq975 = np.nanpercentile(bg, 97.5, axis=0)
        eq025 = np.nanpercentile(be, 2.5, axis=0)
        eq975 = np.nanpercentile(be, 97.5, axis=0)
        se_g = np.nanstd(bg, axis=0)
        se_e = np.nanstd(be, axis=0)

        # Top: band overlay
        ax = axes[0, j]
        ax.fill_between(s_grid, gq025, gq975, color="#9fb6d6", alpha=0.45,
                        label="G 95% band")
        ax.fill_between(s_grid, eq025, eq975, color="#f6c9a3", alpha=0.45,
                        label="E 95% band")
        ax.plot(s_grid, dg, color="#1f3a5f", lw=1.6, label=r"$D_g$ point")
        ax.plot(s_grid, de, color="#a04000", lw=1.6, label=r"$D_e$ point")
        # Empirical envelope ticks
        s_min, s_max = float(s.min()), float(s.max())
        ax.axvline(s_min, color="0.5", lw=0.6, ls=":")
        ax.axvline(s_max, color="0.5", lw=0.6, ls=":")
        ax.set_xlabel(r"$s_j$", fontsize=8)
        if j == 0:
            ax.set_ylabel(r"$D$", fontsize=9)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(
            f"{wname}\n"
            f"med SE_G={np.nanmedian(se_g):.3f}, "
            f"med SE_E={np.nanmedian(se_e):.3f}",
            fontsize=9,
        )
        if j == 0:
            ax.legend(fontsize=7, loc="lower right")
        ax.tick_params(labelsize=7)

        # Bottom: SE profile vs point-estimate D
        ax = axes[1, j]
        order_g = np.argsort(dg)
        order_e = np.argsort(de)
        ax.plot(dg[order_g], se_g[order_g], color="#1f3a5f", lw=1.6,
                label="G")
        ax.plot(de[order_e], se_e[order_e], color="#a04000", lw=1.6,
                label="E")
        ax.set_xlim(0, 1)
        ax.set_xlabel(r"$D$ point estimate", fontsize=8)
        if j == 0:
            ax.set_ylabel(r"$SE(D)$", fontsize=9)
        ax.set_title(f"max SE_G={np.nanmax(se_g):.3f}, max SE_E={np.nanmax(se_e):.3f}",
                     fontsize=9)
        ax.legend(fontsize=7, loc="upper right")
        ax.tick_params(labelsize=7)

        # Stash for JSON
        # Probe SEs at canonical D points (interpret each mapping on its own grid)
        probe_targets = np.array([0.10, 0.50, 0.90])
        probe_g, probe_e = [], []
        for t in probe_targets:
            ig = int(np.argmin(np.abs(dg - t)))
            ie = int(np.argmin(np.abs(de - t)))
            probe_g.append(float(se_g[ig]))
            probe_e.append(float(se_e[ie]))
        w["summary"] = {
            "G": {
                "median_se": float(np.nanmedian(se_g)),
                "max_se": float(np.nanmax(se_g)),
                "se_at_d_0.10": probe_g[0],
                "se_at_d_0.50": probe_g[1],
                "se_at_d_0.90": probe_g[2],
            },
            "E": {
                "median_se": float(np.nanmedian(se_e)),
                "max_se": float(np.nanmax(se_e)),
                "se_at_d_0.10": probe_e[0],
                "se_at_d_0.50": probe_e[1],
                "se_at_d_0.90": probe_e[2],
            },
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
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--n-grid", type=int, default=200)
    p.add_argument("--grid-extrap", type=float, default=0.5)
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

    print(f"[mapping_g_e_bootstrap] slug={args.slug}")
    print(f"  windows: {windows}  n_boot={args.n_boot}  n_grid={args.n_grid}")

    from src.experiment import prepare_data
    from src.metrics.short_block import _apply_mapping_e
    cache_dir = PROJECT_ROOT / "data" / "usgs_cache"
    _, monthly_1d = prepare_data(cache_dir)

    window_data = []
    for w_name in windows:
        samples = historical_window_samples(monthly_1d, w_name)
        mu = float(samples.mean())
        sigma = float(samples.std(ddof=1))
        pad = args.grid_extrap * sigma
        s_grid = np.linspace(samples.min() - pad, samples.max() + pad, args.n_grid)
        d_point_g = stats.norm.cdf((s_grid - mu) / sigma)
        d_point_e = np.array(
            [_apply_mapping_e(float(x), np.sort(samples)) for x in s_grid]
        )
        boot_g, boot_e = bootstrap_g_and_e(
            samples, s_grid, args.n_boot, args.seed
        )
        window_data.append({
            "window": w_name,
            "samples": samples,
            "s_grid": s_grid,
            "d_point_g": d_point_g,
            "d_point_e": d_point_e,
            "boot_g": boot_g,
            "boot_e": boot_e,
        })
        print(
            f"  {w_name:24s} "
            f"med SE_G={np.nanmedian(np.nanstd(boot_g, axis=0)):.3f}  "
            f"max SE_G={np.nanmax(np.nanstd(boot_g, axis=0)):.3f}  "
            f"med SE_E={np.nanmedian(np.nanstd(boot_e, axis=0)):.3f}  "
            f"max SE_E={np.nanmax(np.nanstd(boot_e, axis=0)):.3f}"
        )

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    cfg = json.loads((in_dir / "config.json").read_text())
    title = (
        f"Mapping G vs E bootstrap uncertainty (n_boot={args.n_boot}) — "
        f"{cfg.get('metric_set','?')} preset, slug {args.slug}"
    )
    out_png = fig_dir / "mapping_g_e_bootstrap.png"
    render_compare(window_data, out_png, title)
    print(f"  wrote {out_png}")
    print(f"  wrote {out_png.with_suffix('.pdf')}")

    summary = {
        "slug": args.slug,
        "metric_set": cfg.get("metric_set"),
        "n_boot": args.n_boot,
        "n_grid": args.n_grid,
        "windows": [
            {"window": w["window"], "n": int(w["samples"].size), **w["summary"]}
            for w in window_data
        ],
    }
    json_path = fig_dir / "mapping_g_e_bootstrap.json"
    json_path.write_text(json.dumps(summary, indent=2))
    print(f"  wrote {json_path}")


if __name__ == "__main__":
    main()
