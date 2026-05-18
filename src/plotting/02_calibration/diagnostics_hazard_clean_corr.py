"""Diagnostic figure for the Tier-G hazard-clean pool correlation structure.

Renders a side-by-side comparison of pairwise Spearman correlations
between the seven hazard-clean candidate metrics, computed on:

* historical T-blocks (the 69 stride-1 5-yr blocks of the 73-WY record),
* baseline Kirsch ensemble (10 000 independent realisations at T=5),

plus a third panel showing the 7×7 pairwise scatter matrix on the
historical blocks so the dependence structure is visible.

Reads:
  outputs/02_calibration/t_sensitivity_historical/T05/block_chars_extended.csv
  outputs/03_kirsch_library/build_library_extended/n10000_t5_ssi3-12_s42/characteristics_extended.npz

Writes:
  figures/02_calibration/diagnostics_for_decision/df5_hazard_clean_corr.{png,pdf}
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.metrics.extended import HAZARD_CLEAN_METRICS  # noqa: E402
from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import COLORS, apply_style  # noqa: E402


def main():
    apply_style()
    fig_dir = stage_figure_dir("02_calibration", "diagnostics_for_decision")

    # --- Historical T=5 block-chars ---
    hist = pd.read_csv(stage_output_dir(
        "02_calibration", "t_sensitivity_historical", slug="T05", create=False,
    ) / "block_chars_extended.csv")
    cols = [m for m in HAZARD_CLEAN_METRICS if m in hist.columns]
    if len(cols) != len(HAZARD_CLEAN_METRICS):
        print(f"[df5] missing hazard-clean metrics in hist: "
              f"{set(HAZARD_CLEAN_METRICS) - set(cols)}")
    hist_sub = hist[cols].astype(float)

    # --- Kirsch T=5 ensemble ---
    kp = stage_output_dir(
        "03_kirsch_library", "build_library_extended",
        slug="n10000_t5_ssi3-12_s42", create=False,
    ) / "characteristics_extended.npz"
    if kp.exists():
        z = np.load(kp, allow_pickle=False)
        kir_full = pd.DataFrame(z["values"],
                                columns=[str(n) for n in z["metric_names"]])
        kir_cols = [m for m in cols if m in kir_full.columns]
        kir_sub = kir_full[kir_cols].astype(float)
    else:
        kir_sub = pd.DataFrame(columns=cols)

    rho_hist = hist_sub.corr(method="spearman")
    rho_kir = (kir_sub.corr(method="spearman")
               if not kir_sub.empty else pd.DataFrame())

    # --- Layout: 2 panels (heatmaps) on top + 1 panel (scatter matrix) below
    fig = plt.figure(figsize=(13.0, 11.5))
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[1.0, 1.7],
        hspace=0.45, wspace=0.30,
    )

    def _draw_heat(ax, rho, title):
        if rho.empty:
            ax.axis("off")
            ax.set_title(f"{title}: no data", fontsize=10)
            return None
        arr = rho.values
        im = ax.imshow(np.abs(arr), cmap="RdYlBu_r", vmin=0.5, vmax=1.0,
                       aspect="auto")
        ax.set_xticks(range(len(rho.columns)))
        ax.set_yticks(range(len(rho.columns)))
        ax.set_xticklabels(rho.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels(rho.columns, fontsize=7)
        for i in range(len(rho.columns)):
            for j in range(len(rho.columns)):
                v = arr[i, j]
                col = "white" if abs(v) >= 0.85 else "black"
                ax.text(j, i, f"{abs(v):.2f}", ha="center", va="center",
                        fontsize=7, color=col)
        ax.set_title(title, fontsize=10)
        return im

    ax0 = fig.add_subplot(gs[0, 0])
    im0 = _draw_heat(ax0, rho_hist,
                     "Hazard-clean pool — |Spearman ρ| on 69 historical 5-yr blocks")
    ax1 = fig.add_subplot(gs[0, 1])
    im1 = _draw_heat(ax1, rho_kir,
                     "Hazard-clean pool — |Spearman ρ| on 10 000 Kirsch traces (T=5)")

    if im0 is not None:
        cb = fig.colorbar(im0, ax=[ax0, ax1], shrink=0.8, pad=0.02,
                          location="right")
        cb.set_label("|Spearman ρ|", fontsize=9)

    # --- Pairwise scatter matrix (historical blocks only) ---
    n = len(cols)
    inner_gs = gs[1, :].subgridspec(n, n, wspace=0.08, hspace=0.08)
    for i, mi in enumerate(cols):
        for j, mj in enumerate(cols):
            ax = fig.add_subplot(inner_gs[i, j])
            if i == j:
                # Diagonal — histogram of metric values, hist + kirsch
                vals_h = hist_sub[mi].dropna().values
                vals_k = kir_sub[mi].dropna().values \
                    if mi in kir_sub.columns else np.array([])
                if vals_k.size:
                    ax.hist(vals_k, bins=40, density=True, alpha=0.35,
                            color=COLORS["kde"], edgecolor="none")
                if vals_h.size:
                    ax.hist(vals_h, bins=20, density=True, alpha=0.7,
                            color=COLORS["historical"], edgecolor="none")
                ax.set_yticks([])
            else:
                # Off-diagonal — scatter; lower triangle = hist, upper = kirsch
                if i > j:
                    src_x, src_y = hist_sub[mj], hist_sub[mi]
                    color = COLORS["historical"]
                    sz = 6
                    alpha = 0.7
                else:
                    if mi in kir_sub.columns and mj in kir_sub.columns:
                        src_x, src_y = kir_sub[mj], kir_sub[mi]
                    else:
                        src_x, src_y = pd.Series([]), pd.Series([])
                    color = COLORS["kde"]
                    sz = 1
                    alpha = 0.05
                ax.scatter(src_x, src_y, s=sz, c=color, alpha=alpha,
                           rasterized=True, edgecolors="none")
                ax.set_xticks([])
                ax.set_yticks([])

            if i == n - 1:
                ax.set_xlabel(mj, fontsize=6, rotation=30, ha="right")
            if j == 0:
                ax.set_ylabel(mi, fontsize=6, rotation=0, ha="right",
                              va="center", labelpad=30)
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.tick_params(axis="both", which="both",
                           bottom=False, left=False,
                           labelbottom=(i == n - 1),
                           labelleft=False)

    fig.suptitle(
        "DF5 — Tier-G hazard-clean pool correlation structure at T = 5\n"
        "Top: pairwise |Spearman ρ| on historical blocks (left) and Kirsch ensemble (right).\n"
        "Bottom: pairwise scatter matrix — lower triangle = historical "
        "(black), upper triangle = Kirsch (green), diagonal = marginal histograms.",
        fontsize=10,
    )
    fig.savefig(fig_dir / "df5_hazard_clean_corr.png", dpi=180,
                bbox_inches="tight")
    fig.savefig(fig_dir / "df5_hazard_clean_corr.pdf", bbox_inches="tight")
    print(f"[df5] wrote {fig_dir / 'df5_hazard_clean_corr.png'}")

    # Also write a small CSV summary of the off-diagonal correlation
    # statistics, for quick textual reference.
    summary_rows = []
    for src, mat in (("historical", rho_hist), ("kirsch", rho_kir)):
        if mat.empty:
            continue
        off = mat.where(~np.eye(len(mat), dtype=bool)).abs().stack()
        summary_rows.append({
            "source": src,
            "n_pairs": int(off.size / 2),  # symmetric
            "min_abs_rho": float(off.min()),
            "median_abs_rho": float(off.median()),
            "max_abs_rho": float(off.max()),
        })
    pd.DataFrame(summary_rows).to_csv(
        fig_dir / "df5_hazard_clean_corr_summary.csv", index=False,
    )


if __name__ == "__main__":
    main()
