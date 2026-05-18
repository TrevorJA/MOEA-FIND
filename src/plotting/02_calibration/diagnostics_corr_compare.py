"""DF6 — three-way comparison of Spearman correlations on the Tier-G
hazard-clean pool at T = 5: stride-1 historical (overlap), disjoint
historical tilings (each year used once per tiling), and the 10 000-trace
Kirsch ensemble.

Demonstrates that the high pairwise correlations on the stride-1
historical reference are inflated by overlapping windows of a 73-WY
record. The disjoint-tiling correlations are statistically unbiased
(each year used once per tiling) but noisy at N = 14 per tiling. The
Kirsch ensemble has 10 000 independent observations and produces the
tightest, most reliable estimate of the population correlation.

Output:
  figures/02_calibration/diagnostics_for_decision/df6_corr_compare.{png,pdf}
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.metrics.extended import (  # noqa: E402
    HAZARD_CLEAN_METRICS,
    FullRecordRefs,
    compute_all_candidates,
)
from src.hydrology.historical_blocks import (  # noqa: E402
    resample_disjoint_tilings,
    resample_historical_blocks,
)
from src.metrics.objectives import (  # noqa: E402
    flows_to_series,
    make_ssi_calculator,
)
from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import apply_style  # noqa: E402

T = 5


def _block_chars_2d(blocks_1d: List[np.ndarray], ssi3, ssi12, q80, refs):
    rows = []
    for blk in blocks_1d:
        blk_2d = blk.reshape(T, 12)
        rows.append(compute_all_candidates(blk, blk_2d, ssi3, ssi12, q80,
                                           full_record_refs=refs))
    return pd.DataFrame(rows)


def main():
    apply_style()
    fig_dir = stage_figure_dir("02_calibration", "diagnostics_for_decision")

    # --- Historical record + calibrators ---
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    monthly_2d, monthly_1d = prepare_data(cache_dir)
    ssi3 = make_ssi_calculator(timescale=3)
    ssi3.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    ssi12 = make_ssi_calculator(timescale=12)
    ssi12.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    q80 = float(np.percentile(monthly_1d, 20.0))
    refs = FullRecordRefs.from_full_record(monthly_1d)

    cols = list(HAZARD_CLEAN_METRICS)

    # --- Hist overlap (stride 1) ---
    hist_overlap = _block_chars_2d(
        resample_historical_blocks(monthly_1d, T_years=T, stride=1),
        ssi3, ssi12, q80, refs,
    )[cols]
    rho_hist_ovl = hist_overlap.corr(method="spearman")

    # --- Hist disjoint (T tilings, each year ≤1× per tiling) ---
    tilings = resample_disjoint_tilings(monthly_1d, T_years=T)
    per_tiling_chars: List[pd.DataFrame] = []
    per_tiling_rho: List[pd.DataFrame] = []
    union_chars: List[pd.DataFrame] = []
    for blocks in tilings:
        chars = _block_chars_2d(blocks, ssi3, ssi12, q80, refs)[cols]
        per_tiling_chars.append(chars)
        per_tiling_rho.append(chars.corr(method="spearman"))
        union_chars.append(chars)
    union_chars_df = pd.concat(union_chars, ignore_index=True)
    rho_hist_disj_pooled = union_chars_df.corr(method="spearman")
    # Tiling-averaged correlation (Fisher-z averaging avoided for simplicity;
    # arithmetic mean is fine at small magnitudes around |ρ|≈0.7).
    rho_hist_disj_avg = sum([r.fillna(0) for r in per_tiling_rho]) \
        / len(per_tiling_rho)

    # --- Kirsch ---
    kp = stage_output_dir(
        "03_kirsch_library", "build_library_extended",
        slug="n10000_t5_ssi3-12_s42", create=False,
    ) / "characteristics_extended.npz"
    z = np.load(kp, allow_pickle=False)
    kir_full = pd.DataFrame(z["values"],
                            columns=[str(n) for n in z["metric_names"]])
    rho_kir = kir_full[cols].astype(float).corr(method="spearman")

    # --- Layout: 3 heatmaps side-by-side, plus a bar chart of off-diagonal stats
    fig = plt.figure(figsize=(14.5, 7.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.4],
                          hspace=0.5, wspace=0.30)

    def _draw(ax, rho, title, subtitle):
        if rho.empty:
            ax.axis("off"); return None
        arr = rho.values
        im = ax.imshow(np.abs(arr), cmap="RdYlBu_r", vmin=0.3, vmax=1.0,
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
                        fontsize=6.5, color=col)
        ax.set_title(f"{title}\n{subtitle}", fontsize=9)
        return im

    ax0 = fig.add_subplot(gs[0, 0])
    _draw(ax0, rho_hist_ovl,
          "Historical (stride 1, overlap)",
          f"n_blocks = {len(hist_overlap)}, each year used up to T times")
    ax1 = fig.add_subplot(gs[0, 1])
    _draw(ax1, rho_hist_disj_avg,
          f"Historical (disjoint, mean of {len(tilings)} tilings)",
          f"n_blocks = {len(tilings[0])} per tiling × {len(tilings)} tilings")
    ax2 = fig.add_subplot(gs[0, 2])
    im2 = _draw(ax2, rho_kir,
                "Kirsch ensemble",
                f"n_traces = {len(kir_full)} (independent draws)")

    if im2 is not None:
        cb = fig.colorbar(im2, ax=[ax0, ax1, ax2], shrink=0.6,
                          location="right", pad=0.02)
        cb.set_label("|Spearman ρ|", fontsize=8)

    # --- Off-diagonal summary bar chart ---
    def _offdiag_stats(rho):
        if rho.empty:
            return {"min": np.nan, "median": np.nan, "max": np.nan}
        off = rho.where(~np.eye(len(rho), dtype=bool)).abs().stack()
        return {"min": float(off.min()),
                "median": float(off.median()),
                "max": float(off.max())}

    rows = [
        ("Historical (overlap)", _offdiag_stats(rho_hist_ovl), len(hist_overlap)),
        ("Historical (disjoint, avg)", _offdiag_stats(rho_hist_disj_avg),
         len(tilings[0])),
        ("Historical (disjoint, pooled)",
         _offdiag_stats(rho_hist_disj_pooled), len(union_chars_df)),
        ("Kirsch ensemble", _offdiag_stats(rho_kir), len(kir_full)),
    ]
    summary = pd.DataFrame([
        {"source": s, "n": n, **stats} for s, stats, n in rows
    ])
    summary.to_csv(fig_dir / "df6_corr_compare_summary.csv", index=False)

    ax3 = fig.add_subplot(gs[1, :])
    x = np.arange(len(summary))
    width = 0.25
    ax3.bar(x - width, summary["min"], width, label="min", color="#5fa3d6")
    ax3.bar(x, summary["median"], width, label="median", color="#7f7f7f")
    ax3.bar(x + width, summary["max"], width, label="max", color="#d65a5a")
    ax3.set_xticks(x)
    ax3.set_xticklabels([f"{s}\nN = {summary['n'].iloc[i]}"
                         for i, s in enumerate(summary["source"])],
                        fontsize=8)
    ax3.set_ylabel("|Spearman ρ| over off-diagonal pairs", fontsize=9)
    ax3.set_ylim(0, 1.0)
    ax3.axhline(0.6, color="grey", ls="--", lw=0.6, alpha=0.6)
    ax3.text(3.4, 0.61, "ρ=0.6 cap", fontsize=7, color="grey")
    ax3.axhline(0.75, color="grey", ls=":", lw=0.6, alpha=0.6)
    ax3.text(3.4, 0.76, "ρ=0.75 cap", fontsize=7, color="grey")
    ax3.legend(loc="upper right", fontsize=8, frameon=False)
    ax3.set_title(
        "DF6 — pairwise |Spearman ρ| summary across reference sets "
        "(7 hazard-clean metrics × 21 unique pairs)",
        fontsize=9,
    )

    fig.suptitle(
        "DF6 — Tier-G hazard-clean correlation structure: historical "
        "overlap vs disjoint tilings vs Kirsch ensemble (T = 5)",
        fontsize=10,
    )
    fig.savefig(fig_dir / "df6_corr_compare.png", dpi=180, bbox_inches="tight")
    fig.savefig(fig_dir / "df6_corr_compare.pdf", bbox_inches="tight")
    print(f"[df6] wrote {fig_dir / 'df6_corr_compare.png'}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
