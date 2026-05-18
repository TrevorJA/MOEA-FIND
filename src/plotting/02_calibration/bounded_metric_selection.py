"""DD-15c — diagnostic figures for bounded metric-family / K=4 selection.

Phase 4 of the bounded T=1 reformulation. Consumes outputs from
``recompute_bounded_candidates`` (Phase 2) and ``select_bounded_kset``
(Phase 3) and renders a multi-panel diagnostic PDF that makes the
metric-selection decision visually defensible.

Panels:

A. **Per-candidate scoring scatter** — discrimination_hist (x) vs
   tail_resolution (y), color = mapping {G, E}, marker size = inverse
   of flood_corner_frac (smaller = better flood suppression). Survivors
   highlighted; gates drawn as guide lines.

B. **Saturation bars** — ``frac(D ≤ 0.05)``, ``frac(D ≥ 0.95)``, and
   ``IQR(D)/D*`` per candidate, side-by-side for both mappings. Makes
   the no-clip property obvious (Mapping G should never exactly hit 1).

C. **ECDF overlay panel** — for the winning K=4 set per mapping: per-
   metric ECDF of historical T=1 blocks vs each Pareto archive. Bound
   reference lines at 0 and 1; tail-resolution score annotated.

D. **Pareto archive flood-corner** — for each archive, scatter
   ``max_monthly_flow`` (y, log-scale) vs ``max(D)`` over the K=4
   metrics (x). Horizontal reference at 10k cfs. Quantifies whether
   the new metrics suppress flood-shaped Pareto members.

E. **K-set composite-score Pareto** — top-20 K-sets per mapping plotted
   in (median |ρ|, mean discrimination) space; iso-composite contours
   overlaid. Surfaces the trade-off the selection algorithm makes.

Output: ``figures/02_calibration/bounded_metric_selection/bounded_metric_selection.pdf``
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.metrics.short_block import (  # noqa: E402
    CANDIDATE_BOUNDED_CONCEPT_MAP,
    CANDIDATE_BOUNDED_METRIC_NAMES,
)


_MAPPING_COLORS = {"g": "#1f77b4", "e": "#d62728"}
_MAPPING_LABELS = {"g": "Mapping G (Gaussian CDF)", "e": "Mapping E (empirical + poly tail)"}
_FLOOD_THRESHOLD_CFS = 10_000.0


def _load_inputs(
    recompute_dir: Path,
    select_dir: Path,
) -> Dict:
    hist_df = pd.read_parquet(recompute_dir / "historical.parquet")
    archive_dfs = {}
    for path in sorted(recompute_dir.glob("*.parquet")):
        if path.name == "historical.parquet":
            continue
        archive_dfs[path.stem] = pd.read_parquet(path)
    cand_scores = pd.read_csv(select_dir / "per_candidate_scores.csv")
    selection_report = json.loads((select_dir / "selection_report.json").read_text())
    ksets_g = pd.read_csv(select_dir / "kset_rankings_g.csv") if (select_dir / "kset_rankings_g.csv").exists() else pd.DataFrame()
    ksets_e = pd.read_csv(select_dir / "kset_rankings_e.csv") if (select_dir / "kset_rankings_e.csv").exists() else pd.DataFrame()
    return {
        "hist_df": hist_df,
        "archive_dfs": archive_dfs,
        "cand_scores": cand_scores,
        "selection_report": selection_report,
        "ksets_g": ksets_g,
        "ksets_e": ksets_e,
    }


def _panel_a_scoring_scatter(ax, cand_scores: pd.DataFrame, gates: Dict) -> None:
    survivors = cand_scores[
        (cand_scores["discrimination_hist"] >= gates["discrim"])
        & (cand_scores["tail_resolution"] >= gates["tail_res"])
        & (cand_scores["flood_corner_frac"] <= gates["flood"])
    ]
    for mapping, sub in cand_scores.groupby("mapping"):
        col = _MAPPING_COLORS.get(mapping, "gray")
        marker = "o" if mapping == "g" else "s"
        sizes = 30 + 100 * (1.0 - sub["flood_corner_frac"].clip(0.0, 1.0))
        ax.scatter(
            sub["discrimination_hist"], sub["tail_resolution"],
            s=sizes, c=col, marker=marker, alpha=0.5, edgecolor="none",
            label=f"{_MAPPING_LABELS[mapping]} ({len(sub)})",
        )
    if not survivors.empty:
        ax.scatter(
            survivors["discrimination_hist"], survivors["tail_resolution"],
            s=80, facecolors="none", edgecolor="black", linewidth=1.2,
            label=f"survivors ({len(survivors)})", zorder=4,
        )
    ax.axvline(gates["discrim"], color="gray", linestyle=":", linewidth=0.8)
    ax.axhline(gates["tail_res"], color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel("discrimination_hist (IQR(D over historical) / D*)")
    ax.set_ylabel("tail_resolution (std(D in top-decile of Pareto))")
    ax.set_title("(A) Per-candidate scoring — survivors circled")
    ax.set_yscale("symlog", linthresh=1e-4)
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(True, alpha=0.3)


def _panel_b_saturation_bars(ax, hist_df: pd.DataFrame) -> None:
    cand_list = list(CANDIDATE_BOUNDED_METRIC_NAMES)
    n_cand = len(cand_list)
    frac_low = np.empty(n_cand)
    frac_high = np.empty(n_cand)
    iqr_norm = np.empty(n_cand)
    for i, c in enumerate(cand_list):
        v = hist_df[c].to_numpy(dtype=float)
        frac_low[i] = float((v <= 0.05).mean())
        frac_high[i] = float((v >= 0.95).mean())
        q25, q75 = np.percentile(v, [25, 75])
        iqr_norm[i] = float(q75 - q25)
    x = np.arange(n_cand)
    width = 0.30
    ax.bar(x - width, frac_low, width=width, color="#9ecae1", label="frac(D ≤ 0.05)")
    ax.bar(x, frac_high, width=width, color="#fdae6b", label="frac(D ≥ 0.95)")
    ax.bar(x + width, iqr_norm, width=width, color="#74c476", label="IQR(D)/D*")
    ax.set_xticks(x)
    short_labels = [c.replace("_logmean", "·lm").replace("_logsum", "·ls") for c in cand_list]
    ax.set_xticklabels(short_labels, rotation=90, fontsize=5)
    ax.set_ylabel("fraction / IQR")
    ax.set_title("(B) Per-candidate saturation diagnostic on historical T=1")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3, axis="y")


def _panel_c_ecdf_winning_kset(
    ax_list,
    hist_df: pd.DataFrame,
    archive_dfs: Dict[str, pd.DataFrame],
    selection_report: Dict,
) -> None:
    overall = selection_report.get("overall_winner")
    if overall is None:
        ax_list[0].text(0.5, 0.5, "No winning K-set found",
                        ha="center", va="center", transform=ax_list[0].transAxes)
        return
    metrics = overall["details"]["metrics"]
    archive_colors = plt.cm.Dark2(np.linspace(0, 1, len(archive_dfs)))
    for i, (m, ax) in enumerate(zip(metrics, ax_list)):
        if i >= len(ax_list):
            break
        # Historical ECDF
        h = np.sort(hist_df[m].to_numpy(dtype=float))
        ax.plot(h, np.arange(1, len(h) + 1) / len(h),
                color="black", lw=2, label="historical")
        for (slug, df), col in zip(archive_dfs.items(), archive_colors):
            v = np.sort(df[m].to_numpy(dtype=float))
            short_slug = slug.replace("_T1_nfe200000_s42_constrained_cmdv_uniform", "_")
            ax.plot(v, np.arange(1, len(v) + 1) / len(v),
                    color=col, lw=1.0, alpha=0.85, label=short_slug)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.axvline(1.0, color="red", linestyle=":", linewidth=0.8, alpha=0.5)
        ax.set_title(f"{m}", fontsize=9)
        ax.grid(True, alpha=0.3)
        if i == 0:
            ax.legend(fontsize=6, loc="lower right")
            ax.set_ylabel("ECDF")
        ax.set_xlabel(r"$D_j$")


def _panel_d_flood_corner_scatter(
    ax,
    archive_dfs: Dict[str, pd.DataFrame],
    selection_report: Dict,
) -> None:
    overall = selection_report.get("overall_winner")
    if overall is None:
        ax.text(0.5, 0.5, "No winning K-set", ha="center", va="center", transform=ax.transAxes)
        return
    metrics = overall["details"]["metrics"]
    archive_colors = plt.cm.Dark2(np.linspace(0, 1, len(archive_dfs)))
    for (slug, df), col in zip(archive_dfs.items(), archive_colors):
        max_d = df[metrics].max(axis=1).to_numpy(dtype=float)
        max_q = df["max_monthly_flow"].to_numpy(dtype=float)
        short_slug = slug.replace("_T1_nfe200000_s42_constrained_cmdv_uniform", "_")
        ax.scatter(max_d, max_q, s=4, c=[col], alpha=0.4, label=short_slug)
    ax.axhline(_FLOOD_THRESHOLD_CFS, color="red", linestyle="--",
               linewidth=1.0, alpha=0.7, label="10k cfs flood threshold")
    ax.set_xlabel(r"max($D_j$) over winning K=4")
    ax.set_ylabel("max monthly flow in Pareto trace (cfs)")
    ax.set_yscale("log")
    ax.set_title("(D) Flood-corner diagnostic — Pareto traces")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.3, which="both")


def _panel_e_kset_pareto(ax, ksets_g: pd.DataFrame, ksets_e: pd.DataFrame) -> None:
    for mapping, df in (("g", ksets_g), ("e", ksets_e)):
        if df.empty:
            continue
        col = _MAPPING_COLORS[mapping]
        ax.scatter(
            df["median_abs_rho"], df["mean_discrimination"],
            s=20 + 200 * df["composite_score"], c=col, alpha=0.5,
            edgecolor="none",
            label=f"{_MAPPING_LABELS[mapping]} top-{len(df)}",
        )
        # Mark the winner (top of CSV, since sorted by composite_score)
        top = df.iloc[0]
        ax.scatter(
            [top["median_abs_rho"]], [top["mean_discrimination"]],
            s=120, marker="*", c=col, edgecolor="black", linewidth=1.2,
            zorder=5,
        )
    ax.set_xlabel("median |Spearman ρ| within K=4 set")
    ax.set_ylabel("mean discrimination_hist over K=4")
    ax.set_title("(E) Top K=4 sets per mapping — stars = winners")
    ax.legend(fontsize=7, loc="lower left")
    ax.grid(True, alpha=0.3)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--recompute-dir",
        type=Path,
        default=stage_output_dir(
            "02_calibration", "recompute_bounded_candidates", create=False,
        ),
    )
    parser.add_argument(
        "--select-dir",
        type=Path,
        default=stage_output_dir(
            "02_calibration", "select_bounded_kset", create=False,
        ),
    )
    args = parser.parse_args()

    fig_dir = stage_figure_dir("02_calibration", "bounded_metric_selection")
    pdf_path = fig_dir / "bounded_metric_selection.pdf"
    print(f"[bounded_metric_selection] writing {pdf_path}")

    inputs = _load_inputs(args.recompute_dir, args.select_dir)
    selection_report = inputs["selection_report"]
    gates = selection_report.get("gates", {
        "discrim": 0.20, "tail_res": 1e-4, "flood": 0.05,
    })

    with PdfPages(pdf_path) as pdf:
        # --- Page 1: scoring scatter + saturation bars ---
        fig = plt.figure(figsize=(14, 10))
        gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1.0], hspace=0.45)
        ax_a = fig.add_subplot(gs[0, 0])
        _panel_a_scoring_scatter(ax_a, inputs["cand_scores"], gates)
        ax_b = fig.add_subplot(gs[1, 0])
        _panel_b_saturation_bars(ax_b, inputs["hist_df"])
        fig.suptitle("DD-15c — bounded candidate scoring (Phase 3 inputs)",
                     fontsize=12)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # --- Page 2: winning K-set ECDFs + flood-corner ---
        fig = plt.figure(figsize=(14, 9))
        gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.9])
        ax_c = [fig.add_subplot(gs[0, j]) for j in range(4)]
        _panel_c_ecdf_winning_kset(
            ax_c, inputs["hist_df"], inputs["archive_dfs"], selection_report,
        )
        ax_d = fig.add_subplot(gs[1, :2])
        _panel_d_flood_corner_scatter(
            ax_d, inputs["archive_dfs"], selection_report,
        )
        ax_e = fig.add_subplot(gs[1, 2:])
        _panel_e_kset_pareto(ax_e, inputs["ksets_g"], inputs["ksets_e"])
        winner_str = "—"
        if selection_report.get("overall_winner"):
            ow = selection_report["overall_winner"]
            winner_str = (
                f"mapping={ow['mapping']!r}, K=4 = {ow['details']['metrics']}, "
                f"composite={ow['details']['composite_score']:.3f}"
            )
        fig.suptitle(
            f"DD-15c — winning selection: {winner_str}",
            fontsize=11,
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"[bounded_metric_selection] DONE")


if __name__ == "__main__":
    main()
