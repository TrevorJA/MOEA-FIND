"""Diagnostic figures supporting the (K*, T*) decision in DD-15.

Produces four PNG/PDF figures under
``figures/02_calibration/diagnostics_for_decision/``:

* ``df1_degeneracy_heatmap.{png,pdf}`` — 28-metric × 4-T heatmap of
  ``max(frac_zero, frac_saturated_at_max)`` with cell annotations
  marking which mode dominates. Visual answer to "which metrics are
  usable as Borg axes at each T?".
* ``df2_metric_distributions.{png,pdf}`` — per-metric historical-block
  violins + Kirsch-ensemble KDE overlay across T. Rows: the union of
  metrics appearing in the top-10 K=3 sets at any T (≈ 12 metrics).
  Cells where the metric fails the viability gate are greyed out.
* ``df3_decision_matrix_heatmap.{png,pdf}`` — composite-score heatmap
  for the top-10 K-sets per (K, T) cell, with the recommended pick
  starred.
* ``df4_drought_space_candidates.{png,pdf}`` — 3 × 3 grid of pairwise
  scatters: rows = the recommended K* triple plus the two alternates,
  columns = the 3 metric pairs in each tuple. Historical T-blocks
  overlaid on Kirsch ensemble (sub-sampled).

All four read from existing Stage-1 / Stage-2 / Stage-3 outputs only;
no MOEA archive needed.
"""

from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.metrics.extended import CANDIDATE_METRIC_NAMES, CONCEPT_MAP  # noqa: E402
from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import COLORS, apply_style  # noqa: E402

STAGE = "02_calibration"
HIST_DRIVER = "t_sensitivity_historical"
KIRSCH_STAGE = "03_kirsch_library"
KIRSCH_DRIVER = "build_library_extended"
DECISION_DRIVER = "decision_matrix"
FIGURE_DRIVER = "diagnostics_for_decision"

T_GRID = (5, 10, 20, 30)

#: Tier ordering used to sort metrics on figure rows.
TIER_OF_METRIC = {
    # Tier A — SSI-3 events
    **{m: ("A", i) for i, m in enumerate([
        "frequency", "mean_duration", "max_duration",
        "mean_magnitude", "max_magnitude",
        "mean_severity", "worst_severity", "mean_avg_severity",
        "time_in_drought_fraction",
    ])},
    # Tier B — SSI-12 events
    **{m: ("B", i) for i, m in enumerate([
        "frequency_ssi12", "mean_duration_ssi12", "max_duration_ssi12",
        "mean_magnitude_ssi12", "max_magnitude_ssi12",
        "mean_severity_ssi12", "worst_severity_ssi12",
        "mean_avg_severity_ssi12", "time_in_drought_fraction_ssi12",
    ])},
    # Tier C — recovery / inter-arrival
    **{m: ("C", i) for i, m in enumerate([
        "mean_recovery_time", "mean_drought_free_spell_neg",
    ])},
    # Tier D — FDC tail
    **{m: ("D", i) for i, m in enumerate([
        "q10_flow_neg", "q25_flow_neg",
        "mean_annual_min_neg", "cv_annual_min",
    ])},
    # Tier E — Q80 deficit
    **{m: ("E", i) for i, m in enumerate([
        "q80_events_per_decade", "q80_mean_deficit", "q80_max_deficit",
    ])},
    # Tier F — trend
    "sen_slope_annual_min_neg": ("F", 0),
}

TIER_COLORS = {
    "A": "#1f77b4", "B": "#9ecae1", "C": "#2ca02c",
    "D": "#ff7f0e", "E": "#d62728", "F": "#7f7f7f",
}


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def _hist_dir(T: int) -> Path:
    return stage_output_dir(STAGE, HIST_DRIVER, slug=f"T{T:02d}", create=False)


def _kirsch_dir(T: int, n_traces: int = 10_000, seed: int = 42) -> Path:
    return stage_output_dir(
        KIRSCH_STAGE, KIRSCH_DRIVER,
        slug=f"n{n_traces}_t{T}_ssi3-12_s{seed}", create=False,
    )


def _load_block_chars(T: int) -> pd.DataFrame:
    return pd.read_csv(_hist_dir(T) / "block_chars_extended.csv")


def _load_kirsch_chars(T: int) -> Optional[pd.DataFrame]:
    p = _kirsch_dir(T) / "characteristics_extended.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    return pd.DataFrame(z["values"], columns=[str(n) for n in z["metric_names"]])


def _load_degeneracy(T: int) -> pd.DataFrame:
    return pd.read_csv(_hist_dir(T) / "degeneracy.csv")


def _load_spread(T: int) -> pd.DataFrame:
    return pd.read_csv(_hist_dir(T) / "per_metric_spread.csv")


def _load_decision_matrix() -> pd.DataFrame:
    p = stage_output_dir(STAGE, DECISION_DRIVER, create=False) \
        / "decision_matrix.csv"
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def _load_recommendation() -> Optional[Dict]:
    p = stage_output_dir(STAGE, DECISION_DRIVER, create=False) \
        / "pareto_front_KxT.json"
    return json.loads(p.read_text()) if p.exists() else None


def _viability_flag(spread_row, deg_row) -> str:
    """One-letter viability flag for a (metric, T) pair."""
    fz = float(deg_row.get("frac_zero") or 0.0)
    fn = float(deg_row.get("frac_nan") or 0.0)
    fs = float(deg_row.get("frac_saturated_at_max") or 0.0)
    passes = bool(spread_row["passes_screen"])
    if fz > 0.25:
        return "Z"
    if fn > 0.25:
        return "N"
    if fs > 0.25:
        return "S"
    if not passes:
        return "F"
    return ""


# ---------------------------------------------------------------------------
# Figure DF1 — degeneracy heatmap
# ---------------------------------------------------------------------------


def figure_df1(out_dir: Path):
    apply_style()
    metrics = sorted(CANDIDATE_METRIC_NAMES,
                     key=lambda m: TIER_OF_METRIC.get(m, ("Z", 99)))

    deg_per_T: Dict[int, pd.DataFrame] = {T: _load_degeneracy(T) for T in T_GRID}
    spread_per_T: Dict[int, pd.DataFrame] = {T: _load_spread(T) for T in T_GRID}

    score = np.full((len(metrics), len(T_GRID)), np.nan)
    flag = np.empty_like(score, dtype=object); flag[:] = ""
    for i, m in enumerate(metrics):
        for j, T in enumerate(T_GRID):
            deg = deg_per_T[T].set_index("metric")
            spr = spread_per_T[T].set_index("metric")
            if m not in deg.index:
                continue
            row_deg = deg.loc[m]
            fz = float(row_deg.get("frac_zero") or 0.0)
            fs = float(row_deg.get("frac_saturated_at_max") or 0.0)
            score[i, j] = max(fz, fs)
            row_spr = spr.loc[m] if m in spr.index else {"passes_screen": True}
            flag[i, j] = _viability_flag(row_spr, row_deg)

    fig, ax = plt.subplots(figsize=(5.5, 9.5))
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad("#eeeeee")
    norm = mcolors.Normalize(vmin=0.0, vmax=0.5)
    masked = np.ma.array(score, mask=np.isnan(score))
    im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(range(len(T_GRID)))
    ax.set_xticklabels([f"T={T}" for T in T_GRID])
    ax.set_yticks(range(len(metrics)))

    # Colour metric labels by tier
    yticklabels = []
    for m in metrics:
        tier = TIER_OF_METRIC.get(m, ("?", 0))[0]
        yticklabels.append(m)
    ax.set_yticklabels(yticklabels, fontsize=7)
    for tick, m in zip(ax.get_yticklabels(), metrics):
        tier = TIER_OF_METRIC.get(m, ("?", 0))[0]
        tick.set_color(TIER_COLORS.get(tier, "black"))

    for i in range(len(metrics)):
        for j in range(len(T_GRID)):
            v = score[i, j]
            if np.isnan(v):
                continue
            txt = flag[i, j]
            color = "white" if v > 0.30 else "black"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=8, color=color, fontweight="bold")

    ax.set_xlabel("trace length T (water years)")
    ax.set_title("DF1 — per-metric degeneracy across T\n"
                 "color = max(frac_zero, frac_saturated_at_max); "
                 "Z = zero, S = saturated, N = NaN, F = fail-screen",
                 fontsize=9)

    cb = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.8)
    cb.set_label("max degenerate fraction", fontsize=8)

    legend_patches = [Patch(facecolor=c, edgecolor="black", label=f"Tier {t}")
                      for t, c in TIER_COLORS.items()]
    legend_patches.insert(0, Patch(facecolor="white", edgecolor="black",
                                    label="(metrics colored by tier)"))
    ax.legend(handles=legend_patches, loc="lower left",
              bbox_to_anchor=(1.18, -0.05), fontsize=7, frameon=False)

    fig.tight_layout()
    fig.savefig(out_dir / "df1_degeneracy_heatmap.png", dpi=200,
                bbox_inches="tight")
    fig.savefig(out_dir / "df1_degeneracy_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[DF1] wrote {out_dir / 'df1_degeneracy_heatmap.png'}")


# ---------------------------------------------------------------------------
# Figure DF2 — per-metric distributions across T
# ---------------------------------------------------------------------------


def _gather_top_metrics(decision_df: pd.DataFrame, n_top: int = 10) -> List[str]:
    """Union of metrics appearing in the top-n composite-ranked rows."""
    if decision_df.empty:
        return [
            "cv_annual_min", "worst_severity", "mean_recovery_time",
            "mean_annual_min_neg", "max_duration_ssi12", "frequency_ssi12",
            "mean_duration_ssi12", "q80_mean_deficit", "q10_flow_neg",
            "time_in_drought_fraction", "q25_flow_neg",
            "time_in_drought_fraction_ssi12",
        ]
    top = decision_df.sort_values("composite_score", ascending=False).head(n_top)
    seen: List[str] = []
    for s in top["metrics"]:
        for m in str(s).split("|"):
            if m and m not in seen:
                seen.append(m)
    return seen


def figure_df2(out_dir: Path, decision_df: pd.DataFrame):
    apply_style()
    metrics = _gather_top_metrics(decision_df, n_top=10)
    metrics = sorted(metrics, key=lambda m: TIER_OF_METRIC.get(m, ("Z", 99)))

    hist = {T: _load_block_chars(T) for T in T_GRID}
    kirsch = {T: _load_kirsch_chars(T) for T in T_GRID}
    deg = {T: _load_degeneracy(T).set_index("metric") for T in T_GRID}
    spread = {T: _load_spread(T).set_index("metric") for T in T_GRID}

    n_rows = len(metrics)
    n_cols = len(T_GRID)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(2.2 * n_cols, 1.3 * n_rows),
        squeeze=False,
        sharex=False, sharey=False,
    )

    for i, m in enumerate(metrics):
        all_vals = np.concatenate([
            hist[T][m].dropna().astype(float).values
            for T in T_GRID if m in hist[T].columns
        ] + [
            kirsch[T][m].dropna().astype(float).values
            for T in T_GRID
            if kirsch[T] is not None and m in kirsch[T].columns
        ])
        if all_vals.size == 0:
            continue
        lo, hi = float(np.nanpercentile(all_vals, 1)), \
                 float(np.nanpercentile(all_vals, 99))
        if hi <= lo:
            hi = lo + 1.0
        margin = 0.05 * (hi - lo)
        ymin, ymax = lo - margin, hi + margin

        for j, T in enumerate(T_GRID):
            ax = axes[i, j]
            hv = hist[T][m].dropna().astype(float).values \
                if m in hist[T].columns else np.array([])
            kv = (kirsch[T][m].dropna().astype(float).values
                  if (kirsch[T] is not None and m in kirsch[T].columns)
                  else np.array([]))

            if kv.size:
                rng = np.random.default_rng(0)
                jitter = rng.uniform(0.45, 1.05, size=min(kv.size, 2000))
                kv_sub = rng.choice(kv, size=min(kv.size, 2000), replace=False)
                ax.scatter(jitter, kv_sub, s=2, alpha=0.10,
                           color=COLORS["kde"], rasterized=True)
            if hv.size:
                ax.violinplot(hv, positions=[0.0], widths=0.7,
                              showextrema=False, showmedians=True)

            row_deg = deg[T].loc[m] if m in deg[T].index \
                else pd.Series({"frac_zero": 0, "frac_nan": 0,
                                "frac_saturated_at_max": 0})
            row_spr = spread[T].loc[m] if m in spread[T].index \
                else pd.Series({"passes_screen": True})
            flag = _viability_flag(row_spr, row_deg)

            if flag in ("Z", "N", "S", "F"):
                ax.set_facecolor("#fff0f0")
                ax.text(0.98, 0.95, flag, transform=ax.transAxes,
                        ha="right", va="top", fontsize=10,
                        color="#a30000", fontweight="bold")

            ax.set_xlim(-0.5, 1.5)
            ax.set_xticks([])
            ax.set_ylim(ymin, ymax)
            ax.tick_params(axis="y", labelsize=7)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if i == 0:
                ax.set_title(f"T = {T}", fontsize=9)
            if j == 0:
                tier = TIER_OF_METRIC.get(m, ("?", 0))[0]
                ax.set_ylabel(f"[{tier}] {m}", fontsize=7,
                              color=TIER_COLORS.get(tier, "black"))

    fig.suptitle("DF2 — historical T-block violins (black) and "
                 "Kirsch ensemble jitter (green) across T\n"
                 "Red overlays mark degenerate (Z/N/S) or screen-fail (F) cells",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_dir / "df2_metric_distributions.png", dpi=200,
                bbox_inches="tight")
    fig.savefig(out_dir / "df2_metric_distributions.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[DF2] wrote {out_dir / 'df2_metric_distributions.png'}")


# ---------------------------------------------------------------------------
# Figure DF3 — decision matrix heatmap
# ---------------------------------------------------------------------------


def _shorten(name: str, width: int = 14) -> str:
    """Compact name for use inside a heat-map cell."""
    repl = {
        "time_in_drought_fraction": "tidf",
        "time_in_drought_fraction_ssi12": "tidf12",
        "mean_recovery_time": "recov",
        "mean_drought_free_spell_neg": "ifree-",
        "mean_annual_min_neg": "amin-",
        "cv_annual_min": "cvAmin",
        "worst_severity": "wSev",
        "mean_severity": "mSev",
        "mean_severity_ssi12": "mSev12",
        "worst_severity_ssi12": "wSev12",
        "mean_avg_severity": "mAvg",
        "mean_avg_severity_ssi12": "mAvg12",
        "mean_magnitude": "mMag",
        "max_magnitude": "xMag",
        "mean_magnitude_ssi12": "mMag12",
        "max_magnitude_ssi12": "xMag12",
        "mean_duration": "mDur",
        "max_duration": "xDur",
        "mean_duration_ssi12": "mDur12",
        "max_duration_ssi12": "xDur12",
        "frequency": "freq",
        "frequency_ssi12": "freq12",
        "q10_flow_neg": "q10-",
        "q25_flow_neg": "q25-",
        "q80_mean_deficit": "q80mDef",
        "q80_max_deficit": "q80xDef",
        "q80_events_per_decade": "q80fr",
        "sen_slope_annual_min_neg": "senSlp",
    }
    return repl.get(name, name[:width])


def figure_df3(out_dir: Path, decision_df: pd.DataFrame, rec: Optional[Dict]):
    apply_style()
    if decision_df.empty:
        print("[DF3] decision_matrix.csv missing; skipping")
        return

    fig, axes = plt.subplots(2, 1, figsize=(11.0, 9.0), squeeze=False)
    cmap = plt.get_cmap("viridis")

    rec_K = rec_T = None
    rec_metrics = ""
    if rec is not None:
        rec_K = int(rec["recommended"]["K"])
        rec_T = int(rec["recommended"]["T_years"])
        rec_metrics = "|".join(rec["recommended"]["metrics"])

    for ax, K in zip([axes[0, 0], axes[1, 0]], (3, 4)):
        sub = decision_df[decision_df["K"] == K].copy()
        if sub.empty:
            ax.axis("off"); ax.set_title(f"K = {K}: no entries"); continue
        sub = sub.sort_values(["T_years", "composite_score"],
                              ascending=[True, False])
        # Re-rank within each T column from 1 to N
        sub["display_rank"] = sub.groupby("T_years")["composite_score"] \
            .rank(method="dense", ascending=False).astype(int)
        sub = sub[sub["display_rank"] <= 10]

        wide = sub.pivot_table(
            index="display_rank", columns="T_years",
            values="composite_score", aggfunc="first",
        )
        wide = wide.sort_index()
        im = ax.imshow(wide.values, aspect="auto", cmap=cmap,
                       vmin=0.0, vmax=1.0)
        ax.set_xticks(range(len(wide.columns)))
        ax.set_xticklabels([f"T={T}" for T in wide.columns])
        ax.set_yticks(range(len(wide.index)))
        ax.set_yticklabels([f"#{r}" for r in wide.index], fontsize=8)
        ax.set_title(f"K = {K}", fontsize=11)

        # Annotate metric labels
        labels = sub.pivot_table(
            index="display_rank", columns="T_years",
            values="metrics", aggfunc="first",
        ).reindex_like(wide)
        for i, r in enumerate(wide.index):
            for j, T in enumerate(wide.columns):
                v = wide.iat[i, j]
                if pd.isna(v):
                    continue
                metr = str(labels.iat[i, j]).split("|")
                short = "\n".join([_shorten(m) for m in metr])
                txt_color = "white" if v < 0.45 else "black"
                ax.text(j, i, short, ha="center", va="center",
                        fontsize=6.2, color=txt_color, linespacing=0.95)
                # Star the recommended pick
                if (rec_K == K and rec_T == T
                        and labels.iat[i, j] == rec_metrics):
                    ax.scatter(j, i, marker="*", s=180,
                               edgecolor="white", linewidth=1.4,
                               facecolor="gold", zorder=10)

        ax.set_xlabel("T (water years)")
        ax.set_ylabel("rank within (T, K) cell")

    fig.colorbar(im, ax=axes[0].tolist(), shrink=0.8,
                 label="composite score")
    fig.suptitle("DF3 — joint K × T decision matrix "
                 "(top-10 by composite score per cell; ★ = recommended)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_dir / "df3_decision_matrix_heatmap.png", dpi=200,
                bbox_inches="tight")
    fig.savefig(out_dir / "df3_decision_matrix_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[DF3] wrote {out_dir / 'df3_decision_matrix_heatmap.png'}")


# ---------------------------------------------------------------------------
# Figure DF4 — drought-space coverage for the top candidates
# ---------------------------------------------------------------------------


def figure_df4(out_dir: Path, rec: Optional[Dict]):
    apply_style()
    if rec is None:
        print("[DF4] no recommendation; skipping")
        return

    candidates: List[Dict] = []
    candidates.append({
        "label": "RECOMMENDED",
        "K": int(rec["recommended"]["K"]),
        "T": int(rec["recommended"]["T_years"]),
        "metrics": list(rec["recommended"]["metrics"]),
        "score": float(rec["recommended"]["composite_score"]),
    })
    for alt in rec.get("alternates", []):
        candidates.append({
            "label": "ALTERNATE",
            "K": int(alt["K"]),
            "T": int(alt["T_years"]),
            "metrics": list(alt["metrics"]),
            "score": float(alt["composite_score"]),
        })
    n_rows = len(candidates)
    if n_rows == 0:
        print("[DF4] no candidates; skipping")
        return

    fig, axes = plt.subplots(
        n_rows, 3,
        figsize=(11.5, 3.2 * n_rows),
        squeeze=False,
    )

    for i, cand in enumerate(candidates):
        T = cand["T"]
        metrics = cand["metrics"][:3]  # 3 axes
        if len(metrics) < 3:
            for ax in axes[i]:
                ax.axis("off")
            continue
        hist = _load_block_chars(T)
        kir = _load_kirsch_chars(T)

        pairs = list(combinations(metrics, 2))
        for j, (a, b) in enumerate(pairs):
            ax = axes[i, j]
            if kir is not None and a in kir.columns and b in kir.columns:
                ax.scatter(kir[a], kir[b], s=2, alpha=0.06,
                           color=COLORS["muted"], rasterized=True)
            if a in hist.columns and b in hist.columns:
                ax.scatter(hist[a], hist[b], s=22, alpha=0.85,
                           color=COLORS["historical"], edgecolor="white",
                           linewidth=0.4)
            ax.set_xlabel(a, fontsize=8)
            ax.set_ylabel(b, fontsize=8)
            ax.tick_params(axis="both", labelsize=7)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)

        title = (f"{cand['label']}: K={cand['K']}, T={cand['T']} yr  |  "
                 f"composite = {cand['score']:.3f}")
        axes[i, 0].annotate(
            title,
            xy=(0.0, 1.10), xycoords="axes fraction",
            ha="left", va="bottom", fontsize=9,
            fontweight="bold" if cand["label"] == "RECOMMENDED" else "normal",
            color="black" if cand["label"] == "RECOMMENDED" else "#555",
        )

    legend_handles = [
        Patch(facecolor=COLORS["historical"], edgecolor="white",
              label="Historical T-blocks"),
        Patch(facecolor=COLORS["muted"], edgecolor="none",
              alpha=0.4, label="Baseline Kirsch ensemble (10 k)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=2, fontsize=9, frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("DF4 — drought-space coverage for the recommended K-set "
                 "and two alternates",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    fig.savefig(out_dir / "df4_drought_space_candidates.png", dpi=200,
                bbox_inches="tight")
    fig.savefig(out_dir / "df4_drought_space_candidates.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"[DF4] wrote {out_dir / 'df4_drought_space_candidates.png'}")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--figures", type=str, nargs="*",
                   default=["df1", "df2", "df3", "df4"],
                   help="Subset of figures to render (default: all).")
    args = p.parse_args()

    out_dir = stage_figure_dir(STAGE, FIGURE_DRIVER)
    print(f"[diagnostics_for_decision] → {out_dir}")

    decision_df = _load_decision_matrix()
    rec = _load_recommendation()

    if "df1" in args.figures:
        figure_df1(out_dir)
    if "df2" in args.figures:
        figure_df2(out_dir, decision_df)
    if "df3" in args.figures:
        figure_df3(out_dir, decision_df, rec)
    if "df4" in args.figures:
        figure_df4(out_dir, rec)


if __name__ == "__main__":
    main()
