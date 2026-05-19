"""Comprehensive diagnostics figure for a MOEA-FIND archive.

K- and T-agnostic: reads ``objective_keys`` and ``n_years`` from the
archive's ``results.json`` / ``config.json`` and adapts the figure
layout (K×K pairwise grid, K+2 montage panels, diagnostic panels) to
whatever metric set and trace length were used. Tested on T=1 short-
block bounded archives (DD-15c) and on T=10 first-event SSI-3 archives.

Renders a single multi-panel PDF showing:
  - Lower-triangle pairwise Pareto scatter, diagonal histograms
  - Trace montage: extreme corner of each objective + closest-to-D* + median
  - DV-uniformity constraint distribution (any of l2_star / ks / ad)
  - Hyperplane sum diagnostic (DD-11 K+1 device)
  - Per-month flow CDF vs historical envelope

Reads outputs/04_moea_find_single_site/run_moea_find/<slug>/results.json
and writes figures/04_moea_find_single_site/run_moea_find/<slug>/archive_diagnostics.{pdf,png}.
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

WY_MONTHS = ["Oct", "Nov", "Dec", "Jan", "Feb", "Mar",
             "Apr", "May", "Jun", "Jul", "Aug", "Sep"]


def _short_label(name: str, max_len: int = 11) -> str:
    """Shorten a metric name for tight axis labels (registry name → ≤max_len)."""
    if len(name) <= max_len:
        return name
    # Drop common suffixes/prefixes that add no information for axis labels
    s = name.replace("_total_neg", "_tot").replace("_zscore", "_z")
    s = s.replace("summer_recession", "sum_rec").replace("_total", "_tot")
    # Common prefix stripping (do not let "first_event_*" all collapse to
    # the same prefix; keep the discriminating tail).
    for prefix in ("first_event_", "mean_", "max_", "worst_"):
        if s.startswith(prefix) and len(s) - len(prefix) <= max_len:
            return s[len(prefix):]
    return s[:max_len]


def _select_extreme_indices(dm: np.ndarray, anti_ideal: np.ndarray,
                            obj_names: list) -> dict:
    """Pick K+2 representative Pareto solutions (one per objective max +
    closest-to-D* + Pareto median). Returns ordered dict-like list of
    (label, index) tuples in display order."""
    out = []
    for j, name in enumerate(obj_names):
        out.append((f"max_{name}", int(np.argmax(dm[:, j]))))
    d_to_dstar = np.sum(np.abs(dm - anti_ideal[None, :]), axis=1)
    out.append(("closest_to_D*", int(np.argmin(d_to_dstar))))
    centroid = np.median(dm, axis=0)
    out.append(("Pareto median", int(np.argmin(np.sum(np.abs(dm - centroid), axis=1)))))
    return out


def _historical_monthly_bounds(monthly_2d: np.ndarray) -> tuple:
    return monthly_2d.min(axis=0), monthly_2d.max(axis=0), monthly_2d.mean(axis=0)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True)
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    result = json.loads((in_dir / "results.json").read_text())

    dm = np.array(result["drought_metrics"], dtype=float)
    anti_ideal = np.asarray(result["anti_ideal"], dtype=float)
    obj_names = list(result["objective_keys"])
    K = len(obj_names)
    short_names = [_short_label(n) for n in obj_names]

    if K < 1:
        sys.exit("results.json has no objective_keys")

    # Reconstruct pareto traces (1D per Pareto solution)
    if "pareto_traces_1d" in result and result["pareto_traces_1d"]:
        traces_1d = np.array(result["pareto_traces_1d"], dtype=float)
    else:
        t2d = np.array(result.get("pareto_traces_2d", []), dtype=float)
        traces_1d = t2d.reshape(t2d.shape[0], -1) if t2d.ndim == 3 else t2d

    # Constraint behavior (graceful if absent — e.g. constraint_mode=none)
    diag = result.get("constraint_diagnostics", [])
    have_constraint = bool(diag) and "violations" in diag[0] and diag[0]["violations"]
    if have_constraint:
        stat_name = result["dv_constraint_config"]["statistic"]
        tol = diag[0]["violations"][0]["tolerance"]
        hard_cut = diag[0]["violations"][0]["hard_cutoff"]
        stat_vals = np.array([d["violations"][0]["deviation"] for d in diag])
    else:
        stat_name, tol, hard_cut, stat_vals = None, None, None, None

    # Historical envelope: per-water-year-month bounds tiled to match the
    # archive's trace length. Works for any T, including T=10 first-event
    # archives (n_months_per_block = 120, tile = 10).
    from src.experiment import prepare_data
    cache_dir = PROJECT_ROOT / "data" / "usgs_cache"
    monthly_2d, _ = prepare_data(cache_dir)
    n_months_per_block = traces_1d.shape[1]
    n_blocks = max(1, n_months_per_block // 12)
    hmin, hmax, hmean = _historical_monthly_bounds(monthly_2d)
    if n_months_per_block == 12:
        hist_min, hist_max, hist_mean = hmin, hmax, hmean
    else:
        hist_min = np.tile(hmin, n_blocks)
        hist_max = np.tile(hmax, n_blocks)
        hist_mean = np.tile(hmean, n_blocks)

    extreme = _select_extreme_indices(dm, anti_ideal, obj_names)
    n_montage = len(extreme)  # K + 2

    # ---- Adaptive figure size ----
    pairwise_h = max(2.0 * K, 6.0)
    montage_h = 3.0
    diag_h = 3.0
    fig_w = max(2.5 * K, 12.0)
    fig_h = pairwise_h + montage_h + diag_h + 1.0
    fig = plt.figure(figsize=(fig_w, fig_h), constrained_layout=False)

    # Vertical fractions for the three row-blocks. Generous spacing between
    # blocks (0.06) to avoid collision between pairwise x-tick labels and
    # the montage row titles (which are multi-line stat strings).
    h_total = pairwise_h + montage_h + diag_h
    usable = 0.88
    top_pair = 0.96
    bot_pair = top_pair - (pairwise_h / h_total) * usable
    top_mont = bot_pair - 0.06
    bot_mont = top_mont - (montage_h / h_total) * usable
    top_diag = bot_mont - 0.06
    bot_diag = top_diag - (diag_h / h_total) * usable

    # ---- Block 1: K×K pairwise grid ----
    pair_gs = fig.add_gridspec(K, K, top=top_pair, bottom=bot_pair,
                               left=0.06, right=0.97, hspace=0.40, wspace=0.40)
    corner_palette = plt.colormaps.get_cmap("tab10")
    for i in range(K):
        for j in range(K):
            ax = fig.add_subplot(pair_gs[i, j])
            if i == j:
                ax.hist(dm[:, i], bins=60, color="#2b6cb0", alpha=0.85)
                ax.axvline(anti_ideal[i], color="#d62728", linewidth=2,
                           label=f"D*={anti_ideal[i]:.2g}")
                ax.set_title(short_names[i], fontsize=10)
                ax.legend(fontsize=7)
                ax.set_yticks([])
            elif i > j:
                ax.scatter(dm[:, j], dm[:, i], s=2, alpha=0.4, color="#2b6cb0")
                ax.scatter(anti_ideal[j], anti_ideal[i],
                           marker="X", s=70, color="#d62728", zorder=5,
                           edgecolor="black", linewidth=0.5)
                # Highlight extreme corners
                for k_idx, (label, idx) in enumerate(extreme):
                    color = corner_palette(k_idx % 10)
                    ax.scatter(dm[idx, j], dm[idx, i],
                               marker="o", s=35, color=color,
                               edgecolor="white", linewidth=0.6, zorder=4)
                if i == K - 1:
                    ax.set_xlabel(short_names[j], fontsize=9)
                if j == 0:
                    ax.set_ylabel(short_names[i], fontsize=9)
            else:
                ax.axis("off")

    # ---- Block 2: trace montage (K+2 panels) ----
    mont_gs = fig.add_gridspec(1, n_montage,
                               top=top_mont, bottom=bot_mont,
                               left=0.06, right=0.97, wspace=0.30)
    n_months = traces_1d.shape[1]
    months_axis = np.arange(n_months)
    if n_months == 12:
        xticks_pos = np.arange(12)
        xticks_lab = WY_MONTHS
    else:
        xticks_pos = np.arange(0, n_months, 12)
        xticks_lab = [f"yr{i+1}" for i in range(len(xticks_pos))]

    for k, (label, idx) in enumerate(extreme):
        ax = fig.add_subplot(mont_gs[0, k])
        ax.fill_between(months_axis, hist_min, hist_max,
                        color="#bbb", alpha=0.5, label="hist range")
        ax.plot(months_axis, hist_mean, color="#444",
                linestyle="--", linewidth=1, label="hist mean")
        ax.plot(months_axis, traces_1d[idx], color=corner_palette(k % 10),
                linewidth=1.6, marker="o", markersize=3)
        ax.set_xticks(xticks_pos)
        ax.set_xticklabels(xticks_lab, fontsize=7, rotation=45)
        ax.set_yscale("log")
        m = dm[idx]
        title = label.replace("max_", "max ")
        l1 = float(np.sum(np.abs(m - anti_ideal)))
        ax.set_title(f"{title}\nL1(D-D*)={l1:.1f}", fontsize=8)
        if k == 0:
            ax.set_ylabel("monthly flow (cfs, log)", fontsize=9)
            ax.legend(fontsize=7, loc="lower left")

    # ---- Block 3: diagnostic panels ----
    diag_gs = fig.add_gridspec(1, 3, top=top_diag, bottom=bot_diag,
                               left=0.06, right=0.97, wspace=0.30)
    # (a) Constraint statistic distribution
    ax = fig.add_subplot(diag_gs[0, 0])
    if have_constraint:
        ax.hist(stat_vals, bins=60, color="#9467bd", alpha=0.8)
        ax.axvline(tol, color="#e67e22", linestyle="--", linewidth=2,
                   label=f"soft tol = {tol:.2f}")
        ax.axvline(hard_cut, color="#d62728", linestyle="--", linewidth=2,
                   label=f"hard cut = {hard_cut:.2f}")
        ax.set_xlabel(f"{stat_name.upper()} statistic")
        ax.set_ylabel("Pareto count")
        ax.set_title(f"DV-uniformity ({stat_name}): "
                     f"frac>tol={(stat_vals>tol).mean():.2%}")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no DV-uniformity constraint\n(constraint_mode=none/hydrologic)",
                ha="center", va="center", transform=ax.transAxes, fontsize=10)
        ax.set_title("Constraint diagnostic")
        ax.set_xticks([]); ax.set_yticks([])

    # (b) Hyperplane diagnostic — DD-11 K+1 device: D_j + ||D-D*||_1 = ΣD*_j
    ax = fig.add_subplot(diag_gs[0, 1])
    manhattan = np.sum(np.abs(dm - anti_ideal[None, :]), axis=1)
    sums = dm.sum(axis=1) + manhattan  # all K+1 objectives
    expected = float(anti_ideal.sum())
    # Resilient binning: when DD-15c bounded metrics make the identity
    # hold so tightly that the data range collapses below numerical
    # tolerance, expand the histogram range around the expected value
    # so np.histogram doesn't error on zero-width bins.
    sums_min, sums_max = float(sums.min()), float(sums.max())
    range_pad = max(sums_max - sums_min, 1e-3) * 1.05
    hist_range = (
        min(sums_min, expected) - 0.1 * range_pad,
        max(sums_max, expected) + 0.1 * range_pad,
    )
    n_bins = max(10, min(60, int(round(np.sqrt(sums.size)))))
    ax.hist(sums, bins=n_bins, range=hist_range, color="#16a085", alpha=0.8)
    ax.axvline(expected, color="#d62728", linewidth=2,
               label=f"ΣD* = {expected:.3f}")
    ax.set_xlabel(r"$\Sigma_j D_j(x) + \|D-D^*\|_1$")
    ax.set_ylabel("Pareto count")
    ax.set_title(
        f"Hyperplane: mean={sums.mean():.4f}, std={sums.std():.2e}"
    )
    ax.legend(fontsize=8)

    # (c) Flow plausibility CDF
    ax = fig.add_subplot(diag_gs[0, 2])
    hist_flat = monthly_2d.ravel()
    pareto_flat = traces_1d.ravel()
    h_sorted = np.sort(hist_flat)
    p_sorted = np.sort(pareto_flat)
    ax.plot(h_sorted, np.linspace(0, 1, len(h_sorted)),
            color="#444", linewidth=2, label=f"historical (n={len(h_sorted)})")
    ax.plot(p_sorted, np.linspace(0, 1, len(p_sorted)),
            color="#2b6cb0", linewidth=2,
            label=f"Pareto (n={len(p_sorted)})")
    ax.set_xscale("log")
    ax.set_xlabel("monthly flow (cfs)")
    ax.set_ylabel("CDF")
    ax.set_title(f"Flow plausibility — Pareto max={pareto_flat.max():.0f}, "
                 f"hist max={hist_flat.max():.0f}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # ---- Header ----
    cfg = json.loads((in_dir / "config.json").read_text())
    header = (f"{cfg.get('metric_set', '?')} archive — {args.slug}\n"
              f"K={K}  T={cfg.get('n_years', '?')}  "
              f"NFE={cfg['nfe']:,}  seed={cfg['seed']}  algo={cfg['algorithm']}  "
              f"mode={cfg['mode']}  constraint={cfg['constraint_mode']}/"
              f"{cfg.get('statistic','-')}  "
              f"n_pareto={result['n_pareto']}  elapsed={result['elapsed_s']:.1f}s")
    fig.suptitle(header, y=0.995, fontsize=10, fontfamily="monospace")

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    out_pdf = fig_dir / "archive_diagnostics.pdf"
    out_png = fig_dir / "archive_diagnostics.png"
    fig.savefig(out_pdf, dpi=180, bbox_inches="tight")
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_pdf}")
    print(f"  wrote {out_png}")


if __name__ == "__main__":
    main()
