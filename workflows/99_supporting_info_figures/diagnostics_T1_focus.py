"""DD-15b — comprehensive T=1 / T=2 diagnostic figure pack.

Produces seven PNGs under
``figures/02_calibration/diagnostics_for_decision/`` (no PDFs):

* ``fig01_distributions.png`` — per-metric Hist violin + Kirsch density,
  T=1 left column, T=2 right column. 12 metrics × 2 T = 24 panels.
* ``fig02_degeneracy_spread.png`` — three-panel summary: frac_zero,
  frac_saturated_at_max, robust spread score; bars per metric × T.
* ``fig03_correlations.png`` — 2×2 grid of |Spearman ρ| heatmaps with
  cell numbers; rows = T={1,2}, cols = {Hist, Kirsch}.
* ``fig04_scatter_T1.png`` — 12×12 lower-triangle pairwise scatter at
  T=1 with Hist and Kirsch overlaid; diagonal = marginal histograms.
* ``fig05_kset_viability.png`` — strict-rung K-set count by correlation
  cap, K ∈ {3, 4, 5, 6}, T ∈ {1, 2} — the "K-set headroom" view.
* ``fig06_top_ksets.png`` — top-10 K=4 / K=5 / K=6 candidate tuples at
  T=1, ranked by composite score; annotated with concepts and max |ρ|.
* ``fig07_exemplar_traces.png`` — hand-picked T=1 hydrographs from
  history and Kirsch chosen to span extremes of the leading metrics,
  with metric values annotated.
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import COLORS, apply_style  # noqa: E402
from src.metrics.short_block import (  # noqa: E402
    SHORT_BLOCK_CONCEPT_MAP,
    SHORT_BLOCK_METRIC_NAMES,
)


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def _hist_chars(T: int) -> Optional[pd.DataFrame]:
    p = stage_output_dir(
        "02_calibration", "short_block_screening", slug=f"T{T}",
        create=False,
    ) / "block_chars_short.csv"
    return pd.read_csv(p) if p.exists() else None


def _kirsch_chars(T: int, n_traces: int = 10_000, seed: int = 42) -> Optional[pd.DataFrame]:
    p = stage_output_dir(
        "03_kirsch_library", "build_short_block_library",
        slug=f"n{n_traces}_s{seed}", create=False,
    ) / f"characteristics_short_T{T}.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    return pd.DataFrame(z["values"], columns=[str(n) for n in z["metric_names"]])


def _kirsch_library_36mo(n_traces: int = 10_000, seed: int = 42) -> Optional[np.ndarray]:
    p = stage_output_dir(
        "03_kirsch_library", "build_short_block_library",
        slug=f"n{n_traces}_s{seed}", create=False,
    ) / "library_36mo.npy"
    return np.load(p) if p.exists() else None


# ---------------------------------------------------------------------------
# Fig 01 — distributions per metric, T=1 vs T=2
# ---------------------------------------------------------------------------


def fig01_distributions(out_dir: Path):
    apply_style()
    metrics = list(SHORT_BLOCK_METRIC_NAMES)
    n = len(metrics)
    h1, k1 = _hist_chars(1), _kirsch_chars(1)
    h2, k2 = _hist_chars(2), _kirsch_chars(2)
    if h1 is None or k1 is None:
        print("[fig01] missing artifacts"); return

    # Dynamic layout: 4 columns (2 metrics × 2 T = 4 cells per row); rows
    # = ceil(n_metrics / 2). Each metric occupies 2 adjacent columns
    # (T=1 then T=2), so n_rows = ceil(n / 2).
    n_rows = (n + 1) // 2
    fig, axes = plt.subplots(n_rows, 4, figsize=(13.0, 2.0 * n_rows),
                             squeeze=False)

    rng = np.random.default_rng(0)
    for i, m in enumerate(metrics):
        row = i // 2
        col_offset = (i % 2) * 2
        for j, (T_label, hist_df, kir_df) in enumerate([
            ("T=1", h1, k1), ("T=2", h2, k2),
        ]):
            ax = axes[row, col_offset + j]
            if hist_df is None or m not in hist_df.columns:
                ax.axis("off"); continue
            hv = hist_df[m].dropna().astype(float).values
            kv = kir_df[m].dropna().astype(float).values \
                if (kir_df is not None and m in kir_df.columns) else np.array([])

            if kv.size:
                samp = rng.choice(kv, size=min(kv.size, 2000), replace=False)
                jit = rng.uniform(0.45, 1.05, size=samp.size)
                ax.scatter(jit, samp, s=2, alpha=0.10,
                           color=COLORS["kde"], rasterized=True)
            if hv.size:
                ax.violinplot(hv, positions=[0.0], widths=0.7,
                              showextrema=False, showmedians=True)

            ax.set_xlim(-0.5, 1.4)
            ax.set_xticks([])
            ax.tick_params(axis="y", labelsize=7)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if j == 0:
                ax.set_ylabel(
                    f"{m}\n[{SHORT_BLOCK_CONCEPT_MAP.get(m, '?')}]",
                    fontsize=7,
                )
            ax.set_title(T_label, fontsize=8)

    fig.suptitle(
        "Fig 01 — per-metric distributions at T = 1 vs T = 2\n"
        "blue violin = historical blocks  |  green jitter = Kirsch (10 000 traces)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_dir / "fig01_distributions.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig01] wrote {out_dir / 'fig01_distributions.png'}")


# ---------------------------------------------------------------------------
# Fig 02 — degeneracy and spread bars
# ---------------------------------------------------------------------------


def _spread_score(values: np.ndarray) -> float:
    if values.size < 2:
        return 0.0
    med = float(np.median(values))
    std = float(np.std(values, ddof=1))
    q1, q3 = np.percentile(values, [25.0, 75.0])
    iqr = float(q3 - q1)
    denom = abs(med) + std
    return float(iqr / denom) if denom > 1e-12 else 0.0


def _degeneracy_stats(values: np.ndarray) -> Dict[str, float]:
    finite = np.isfinite(values)
    n = values.size
    if n == 0 or not finite.any():
        return {"frac_zero": np.nan, "frac_nan": 1.0,
                "frac_saturated_at_max": np.nan}
    v = values[finite]
    n_finite = v.size
    n_nan = n - n_finite
    n_zero = int((v == 0.0).sum())
    max_val = float(np.max(v))
    if max_val > 0:
        n_sat = int(np.isclose(v, max_val, rtol=1e-6).sum())
    else:
        n_sat = int((v == max_val).sum())
    return {
        "frac_zero": n_zero / n_finite,
        "frac_nan": n_nan / n,
        "frac_saturated_at_max": n_sat / n_finite,
    }


def fig02_degeneracy_spread(out_dir: Path):
    apply_style()
    metrics = list(SHORT_BLOCK_METRIC_NAMES)
    rows = []
    for T in (1, 2):
        h = _hist_chars(T)
        k = _kirsch_chars(T)
        if h is None:
            continue
        for m in metrics:
            if m not in h.columns:
                continue
            hv = h[m].astype(float).values
            kv = k[m].astype(float).values if (k is not None and m in k.columns) \
                else np.array([])
            d = _degeneracy_stats(hv)
            rows.append({
                "metric": m, "T": f"T={T}", "source": "Hist",
                "spread": _spread_score(hv[np.isfinite(hv)]),
                **d,
            })
            d_k = _degeneracy_stats(kv) if kv.size else {
                "frac_zero": np.nan, "frac_nan": np.nan,
                "frac_saturated_at_max": np.nan,
            }
            rows.append({
                "metric": m, "T": f"T={T}", "source": "Kirsch",
                "spread": _spread_score(kv[np.isfinite(kv)]) if kv.size else np.nan,
                **d_k,
            })
    df = pd.DataFrame(rows)

    fig, axes = plt.subplots(3, 1, figsize=(11.5, 10.0))
    pivot_label = lambda src, T: f"{src} {T}"

    def _bar(ax, col, title, ymax=None):
        # Show 4 bars per metric: Hist T=1, Kirsch T=1, Hist T=2, Kirsch T=2
        x = np.arange(len(metrics))
        width = 0.21
        offsets = {
            ("Hist", "T=1"): -1.5 * width,
            ("Kirsch", "T=1"): -0.5 * width,
            ("Hist", "T=2"): 0.5 * width,
            ("Kirsch", "T=2"): 1.5 * width,
        }
        colors = {
            ("Hist", "T=1"): "#1f77b4", ("Kirsch", "T=1"): "#9ecae1",
            ("Hist", "T=2"): "#2ca02c", ("Kirsch", "T=2"): "#a1d99b",
        }
        for (src, Tlbl), off in offsets.items():
            sub = df[(df["source"] == src) & (df["T"] == Tlbl)]
            sub = sub.set_index("metric").reindex(metrics)
            ax.bar(x + off, sub[col].values, width,
                   label=pivot_label(src, Tlbl), color=colors[(src, Tlbl)])
        ax.set_xticks(x)
        ax.set_xticklabels(metrics, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel(title)
        if ymax is not None:
            ax.set_ylim(0, ymax)
        ax.axhline(0.25, color="grey", ls="--", lw=0.6, alpha=0.6)

    _bar(axes[0], "frac_zero", "frac of blocks with metric = 0", 1.0)
    axes[0].set_title(
        "Fraction of blocks where the metric is exactly zero "
        "(Borg-axis viability gate at 0.25)",
        fontsize=10,
    )
    _bar(axes[1], "frac_saturated_at_max", "frac of blocks at the per-T max", 1.0)
    axes[1].set_title(
        "Fraction of blocks pinned at the per-T max value "
        "(saturation gate at 0.25)",
        fontsize=10,
    )
    _bar(axes[2], "spread", "robust spread = IQR / (|median| + σ)")
    axes[2].set_title("Robust spread score (higher = more discriminating)",
                      fontsize=10)
    axes[0].legend(loc="upper right", fontsize=8, frameon=False)

    fig.suptitle(
        "Fig 02 — degeneracy and spread per metric at T = 1, T = 2 "
        "(Hist + Kirsch)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_dir / "fig02_degeneracy_spread.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig02] wrote {out_dir / 'fig02_degeneracy_spread.png'}")


# ---------------------------------------------------------------------------
# Fig 03 — correlation matrices Hist × Kirsch × T
# ---------------------------------------------------------------------------


def fig03_correlations(out_dir: Path):
    apply_style()
    cols = list(SHORT_BLOCK_METRIC_NAMES)
    sources = []
    for T in (1, 2):
        h, k = _hist_chars(T), _kirsch_chars(T)
        if h is not None and k is not None:
            sources.append((T, h, k))

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 10.0))
    last_im = None
    for r, (T, hdf, kdf) in enumerate(sources):
        for c, (label, df) in enumerate([("Historical", hdf), ("Kirsch", kdf)]):
            ax = axes[r, c]
            avail = [m for m in cols if m in df.columns]
            sub = df[avail].astype(float)
            rho = sub.corr(method="spearman")
            arr = rho.values
            im = ax.imshow(np.abs(arr), cmap="RdYlBu_r",
                           vmin=0.0, vmax=1.0, aspect="auto")
            last_im = im
            ax.set_xticks(range(len(rho.columns)))
            ax.set_yticks(range(len(rho.columns)))
            ax.set_xticklabels(rho.columns, rotation=45, ha="right", fontsize=7)
            ax.set_yticklabels(rho.columns, fontsize=7)
            for i in range(len(rho.columns)):
                for j in range(len(rho.columns)):
                    v = arr[i, j]
                    txt_color = "white" if abs(v) >= 0.85 else "black"
                    ax.text(j, i, f"{abs(v):.2f}",
                            ha="center", va="center",
                            fontsize=6.0, color=txt_color)
            off = rho.where(~np.eye(len(rho), dtype=bool)).abs().stack()
            n = len(df)
            ax.set_title(
                f"T = {T} | {label}\n"
                f"n = {n} | median |ρ| = {off.median():.3f} | "
                f"min |ρ| = {off.min():.3f}",
                fontsize=9,
            )

    if last_im is not None:
        fig.colorbar(last_im, ax=axes.ravel().tolist(),
                     shrink=0.6, pad=0.02, location="right",
                     label="|Spearman ρ|")
    fig.suptitle(
        "Fig 03 — pairwise |Spearman ρ| heatmaps "
        "(rows = T = 1, T = 2; columns = Historical, Kirsch ensemble)",
        fontsize=11,
    )
    fig.savefig(out_dir / "fig03_correlations.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig03] wrote {out_dir / 'fig03_correlations.png'}")


# ---------------------------------------------------------------------------
# Fig 04 — pairwise scatter matrix at T=1
# ---------------------------------------------------------------------------


def _scatter_matrix_metric_subset(rho: pd.DataFrame, hist_chars: pd.DataFrame,
                                   T: int) -> List[str]:
    """Curate the 8 highest-spread, lowest-redundancy metrics for fig04.

    Goals: include all 5 Tier-J additions, retain the most informative
    Tier-H / Tier-I metrics, exclude perfectly-correlated and saturated
    members. Keeps the 17-metric pool's scatter matrix readable
    (8×8 panels rather than 17×17).
    """
    # Mandatory inclusions: Tier-J novelties so user can see them.
    must_include = ["amj_total_neg", "ond_total_neg", "min_6mo_rolling_neg",
                    "aug_zscore", "q90_zscore"]
    # Strong Tier-H complements: distinct concepts, not perfectly
    # redundant, not saturated.
    candidates = ["djf_total_neg", "jja_total_neg", "min_monthly_flow_neg",
                  "summer_recession", "min_annual_zscore",
                  "total_deficit_ssi3_block", "time_in_drought_fraction_block",
                  "q10_zscore"]
    excl = _kset_exclude(T)
    out = [m for m in must_include
           if m in hist_chars.columns and m not in excl]
    for m in candidates:
        if len(out) >= 8:
            break
        if m in hist_chars.columns and m not in excl and m not in out:
            out.append(m)
    return out


def fig04_scatter_T1(out_dir: Path):
    apply_style()
    h = _hist_chars(1)
    k = _kirsch_chars(1)
    if h is None or k is None:
        print("[fig04] missing artifacts"); return
    cols_full = [m for m in SHORT_BLOCK_METRIC_NAMES if m in k.columns]
    rho = k[cols_full].astype(float).corr(method="spearman")
    cols = _scatter_matrix_metric_subset(rho, h, T=1)

    n = len(cols)
    fig, axes = plt.subplots(n, n, figsize=(2.0 * n, 2.0 * n), squeeze=False)
    rng = np.random.default_rng(1)
    kir_idx = rng.choice(len(k), size=min(len(k), 1500), replace=False)

    for i, mi in enumerate(cols):
        for j, mj in enumerate(cols):
            ax = axes[i, j]
            if i == j:
                # marginal histograms — Hist and Kirsch overlaid
                if mi in h.columns:
                    hv = h[mi].dropna().astype(float).values
                    kv = k[mi].dropna().astype(float).values
                    ax.hist(kv, bins=40, density=True, alpha=0.4,
                            color=COLORS["kde"], edgecolor="none")
                    ax.hist(hv, bins=20, density=True, alpha=0.7,
                            color=COLORS["historical"], edgecolor="none")
                    ax.set_yticks([])
            elif i > j:
                # lower triangle: Hist (large black) + Kirsch (small grey)
                if mj in k.columns and mi in k.columns:
                    ax.scatter(k[mj].iloc[kir_idx], k[mi].iloc[kir_idx],
                               s=2, alpha=0.10, color=COLORS["muted"],
                               rasterized=True)
                if mj in h.columns and mi in h.columns:
                    ax.scatter(h[mj], h[mi], s=14, alpha=0.85,
                               color=COLORS["historical"],
                               edgecolor="white", linewidth=0.3)
                ax.set_xticks([]); ax.set_yticks([])
            else:
                # upper triangle: |ρ| value
                if mi in h.columns and mj in h.columns:
                    rho_h = h[[mi, mj]].astype(float).corr(method="spearman").iat[0, 1]
                    rho_k = k[[mi, mj]].astype(float).corr(method="spearman").iat[0, 1]
                    ax.text(0.5, 0.65, f"H |ρ| = {abs(rho_h):.2f}",
                            transform=ax.transAxes, ha="center",
                            fontsize=8,
                            color="black",
                            weight="bold" if abs(rho_h) >= 0.6 else "normal")
                    ax.text(0.5, 0.35, f"K |ρ| = {abs(rho_k):.2f}",
                            transform=ax.transAxes, ha="center",
                            fontsize=8,
                            color=COLORS["kde"],
                            weight="bold" if abs(rho_k) >= 0.6 else "normal")
                ax.set_xticks([]); ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_alpha(0.2)

            if i == n - 1:
                ax.set_xlabel(mj, fontsize=6, rotation=30, ha="right")
            if j == 0:
                ax.set_ylabel(mi, fontsize=6, rotation=30, ha="right")
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)

    fig.suptitle(
        f"Fig 04 — pairwise scatter matrix at T = 1, curated {n}-metric "
        "subset (Tier-J additions + strong Tier-H complements)\n"
        "lower triangle: scatter (black = Hist, grey = Kirsch)  |  "
        "upper triangle: |Spearman ρ| (H = Hist, K = Kirsch)  |  "
        "diagonal: marginal histograms",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_dir / "fig04_scatter_T1.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig04] wrote {out_dir / 'fig04_scatter_T1.png'}")


# ---------------------------------------------------------------------------
# Fig 05 — K-set viability bars, K=3..6, T=1 vs T=2
# ---------------------------------------------------------------------------


def _count_strict_ksets(rho: pd.DataFrame, K: int, cap: float,
                         exclude_metrics: Optional[set] = None) -> int:
    metrics = [m for m in rho.columns
               if (exclude_metrics is None or m not in exclude_metrics)]
    if len(metrics) < K:
        return 0
    sub_rho = rho.loc[metrics, metrics].abs().values.copy()
    np.fill_diagonal(sub_rho, 0.0)
    concept_lookup = {m: SHORT_BLOCK_CONCEPT_MAP.get(m, "?") for m in metrics}
    name_to_idx = {m: i for i, m in enumerate(metrics)}
    count = 0
    for combo in combinations(metrics, K):
        if len({concept_lookup[m] for m in combo}) < K:
            continue
        idxs = [name_to_idx[m] for m in combo]
        max_off = 0.0
        for i in range(K):
            for j in range(i + 1, K):
                v = sub_rho[idxs[i], idxs[j]]
                if v > max_off: max_off = v
        if max_off < cap:
            count += 1
    return count


#: Metrics excluded from K-set enumeration (kept in distributions/correlations
#: figures for transparency, but never recommended as Borg axes):
#:
#: * ``total_flow_neg`` at T=1 — deterministic transform of
#:   ``min_annual_zscore`` (perfect ρ=1.0 at T=1).
#: * ``min_ssi3_neg`` at all T — saturates at the SynHydro SSI-3
#:   calibration floor (~6.36) in 11–69% of blocks. The same-value
#:   pin makes it a non-discriminating Borg axis.
ALWAYS_EXCLUDE = {"min_ssi3_neg"}


def _kset_exclude(T: int) -> set:
    excl = set(ALWAYS_EXCLUDE)
    if T == 1:
        excl.add("total_flow_neg")
    return excl


def fig05_kset_viability(out_dir: Path):
    apply_style()
    caps = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    Ks = [3, 4, 5, 6]
    settings = []
    for T in (1, 2):
        k = _kirsch_chars(T)
        if k is None:
            continue
        cols = [m for m in SHORT_BLOCK_METRIC_NAMES if m in k.columns]
        rho = k[cols].astype(float).corr(method="spearman")
        excl = _kset_exclude(T)
        settings.append((f"T={T}", rho, excl))

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5))
    axes = axes.flatten()
    colors = {"T=1": "#1f77b4", "T=2": "#2ca02c"}
    width = 0.4
    for ax_idx, K in enumerate(Ks):
        ax = axes[ax_idx]
        x = np.arange(len(caps))
        for s_idx, (label, rho, excl) in enumerate(settings):
            counts = [_count_strict_ksets(rho, K=K, cap=c, exclude_metrics=excl)
                      for c in caps]
            offset = (s_idx - 0.5) * width
            ax.bar(x + offset, counts, width, label=label,
                   color=colors[label])
            for i, c in enumerate(counts):
                ax.text(x[i] + offset, c + max(counts) * 0.01 + 0.5,
                        str(c), ha="center", va="bottom", fontsize=7)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{c:.1f}" for c in caps])
        ax.set_xlabel("Kirsch correlation cap |ρ|")
        ax.set_ylabel(f"# strict-rung K = {K} sets passing cap")
        ax.set_title(f"K = {K}", fontsize=10)
        ax.legend(fontsize=8, frameon=False)

    n_pool = len(SHORT_BLOCK_METRIC_NAMES)
    n_T1 = n_pool - len(_kset_exclude(1))
    n_T2 = n_pool - len(_kset_exclude(2))
    fig.suptitle(
        "Fig 05 — K-set viability (strict-rung: distinct concepts AND "
        "pairwise |ρ_Kirsch| < cap)\n"
        f"T = 1 effective pool: {n_T1}/{n_pool} metrics (excl "
        f"total_flow_neg redundancy + min_ssi3_neg saturation)  |  "
        f"T = 2 effective pool: {n_T2}/{n_pool} metrics (excl "
        f"min_ssi3_neg saturation)",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_dir / "fig05_kset_viability.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig05] wrote {out_dir / 'fig05_kset_viability.png'}")


# ---------------------------------------------------------------------------
# Fig 06 — top K=4, K=5, K=6 candidate tuples at T=1
# ---------------------------------------------------------------------------


def _enumerate_strict(rho: pd.DataFrame, K: int, cap: float,
                       exclude_metrics: Optional[set] = None,
                       hist_chars: Optional[pd.DataFrame] = None) -> List[Dict]:
    metrics = [m for m in rho.columns
               if (exclude_metrics is None or m not in exclude_metrics)]
    if len(metrics) < K:
        return []
    sub_rho = rho.loc[metrics, metrics]
    abs_vals = sub_rho.abs().values.copy()
    np.fill_diagonal(abs_vals, 0.0)
    name_to_idx = {m: i for i, m in enumerate(metrics)}
    concept_lookup = {m: SHORT_BLOCK_CONCEPT_MAP.get(m, "?") for m in metrics}
    spread_lookup: Dict[str, float] = {}
    if hist_chars is not None:
        for m in metrics:
            if m in hist_chars.columns:
                v = hist_chars[m].astype(float).values
                v = v[np.isfinite(v)]
                spread_lookup[m] = _spread_score(v)

    cands = []
    for combo in combinations(metrics, K):
        if len({concept_lookup[m] for m in combo}) < K:
            continue
        idxs = [name_to_idx[m] for m in combo]
        max_off = 0.0
        for i in range(K):
            for j in range(i + 1, K):
                v = abs_vals[idxs[i], idxs[j]]
                if v > max_off: max_off = v
        if max_off >= cap:
            continue
        spreads = [spread_lookup.get(m, 0.0) for m in combo]
        cands.append({
            "metrics": list(combo),
            "concepts": [concept_lookup[m] for m in combo],
            "min_spread": float(min(spreads)),
            "sum_spread": float(sum(spreads)),
            "max_pairwise_rho": float(max_off),
        })

    if cands:
        max_ms = max(c["min_spread"] for c in cands) or 1.0
        for c in cands:
            score_spread = max(c["min_spread"] / max_ms, 1e-3)
            score_indep = max(1.0 - c["max_pairwise_rho"], 1e-3)
            c["composite"] = float(np.exp(np.log([score_spread, score_indep]).mean()))
    cands.sort(key=lambda c: -c["composite"])
    return cands


def fig06_top_ksets(out_dir: Path):
    """One figure per K cardinality at T=1 (and T=2 for K=5, 6 since
    T=1 cannot admit them). Each figure is large enough to read."""
    apply_style()
    h1, k1 = _hist_chars(1), _kirsch_chars(1)
    h2, k2 = _hist_chars(2), _kirsch_chars(2)

    settings: List[Tuple[str, int, pd.DataFrame, pd.DataFrame, set]] = []
    if h1 is not None and k1 is not None:
        cols1 = [m for m in SHORT_BLOCK_METRIC_NAMES if m in k1.columns]
        rho1 = k1[cols1].astype(float).corr(method="spearman")
        excl1 = _kset_exclude(1)
        for K in (3, 4, 5, 6):
            settings.append((f"T=1, K={K}", K, rho1, h1, excl1))

    n = len(settings)
    fig, axes = plt.subplots(n, 1, figsize=(13.0, 4.5 * n),
                             squeeze=False)
    for s_idx, (label, K, rho, hist_df, excl) in enumerate(settings):
        ax = axes[s_idx, 0]
        cands_per_cap: Dict[float, List[Dict]] = {}
        chosen_cap: Optional[float] = None
        for cap in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
            cands = _enumerate_strict(rho, K, cap, excl, hist_chars=hist_df)
            cands_per_cap[cap] = cands
            if chosen_cap is None and len(cands) >= 6:
                chosen_cap = cap
        if chosen_cap is None:
            chosen_cap = max(cands_per_cap,
                             key=lambda c: len(cands_per_cap[c]))
        cands = cands_per_cap[chosen_cap][:10]

        if not cands:
            ax.axis("off")
            ax.set_title(
                f"{label}: NO strict-rung set up to cap = 0.8",
                fontsize=11,
            )
            continue

        y = np.arange(len(cands))
        scores = [c["composite"] for c in cands]
        rhos = [c["max_pairwise_rho"] for c in cands]
        ax.barh(y, scores, color="#1f77b4", alpha=0.85)
        ax.set_yticks(y)
        labels_y = [
            f"#{i+1}: {' | '.join(c['metrics'])}\n"
            f"      concepts: {', '.join(c['concepts'])}"
            for i, c in enumerate(cands)
        ]
        ax.set_yticklabels(labels_y, fontsize=8.5)
        ax.invert_yaxis()
        ax.set_xlabel("composite score = √(spread_score × (1 − max|ρ|))",
                      fontsize=9)
        ax.set_xlim(0, 1.0)
        for i, (s, r) in enumerate(zip(scores, rhos)):
            ax.text(min(s + 0.01, 0.97), i,
                    f"max|ρ| = {r:.2f}", va="center",
                    fontsize=8, color="grey")
        ax.set_title(
            f"{label} — top {len(cands)} strict-rung candidates "
            f"(cap = {chosen_cap:.1f}, "
            f"{len(cands_per_cap[chosen_cap])} total)",
            fontsize=11, loc="left",
        )
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    fig.suptitle(
        "Fig 06 — top candidate K-set tuples ranked by composite score\n"
        "(Kirsch correlations; strict rung = distinct concepts AND "
        "pairwise |ρ| < cap; min_ssi3_neg excluded for SSI-3 calibration "
        "saturation; total_flow_neg excluded at T=1 for redundancy)",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_dir / "fig06_top_ksets.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig06] wrote {out_dir / 'fig06_top_ksets.png'}")


# ---------------------------------------------------------------------------
# Fig 07 — exemplar T=1 hydrographs spanning metric extremes
# ---------------------------------------------------------------------------


WY_MONTH_LABELS = ["Oct", "Nov", "Dec", "Jan", "Feb", "Mar",
                   "Apr", "May", "Jun", "Jul", "Aug", "Sep"]


def fig07_exemplar_traces(out_dir: Path):
    apply_style()
    h = _hist_chars(1)
    library = _kirsch_library_36mo()
    k = _kirsch_chars(1)
    if h is None or library is None or k is None:
        print("[fig07] missing"); return

    # Reload historical monthly record to extract specific water years.
    cache = PROJECT_ROOT / "outputs" / "data_cache"
    monthly_2d_full, monthly_1d_full = prepare_data(cache)

    # Pick exemplar metrics for axis extremes
    pick_metrics = [
        "djf_total_neg", "min_monthly_flow_neg",
        "summer_recession", "jja_total_neg",
    ]
    rng = np.random.default_rng(7)

    # Hist exemplars: most-extreme on each metric (top-2 of each axis)
    hist_picks: List[Tuple[str, int]] = []
    seen = set()
    for m in pick_metrics:
        order = h.sort_values(m, ascending=False)
        for _, row in order.head(2).iterrows():
            wy = int(row["water_year_index"])
            if wy not in seen:
                hist_picks.append((m, wy))
                seen.add(wy)

    # Kirsch exemplars: similar — top-1 trace by each axis
    kir_picks: List[Tuple[str, int]] = []
    for m in pick_metrics:
        if m in k.columns:
            tid = int(k[m].astype(float).idxmax())
            if tid not in {p[1] for p in kir_picks}:
                kir_picks.append((m, tid))

    n_rows = max(len(hist_picks), len(kir_picks))
    fig, axes = plt.subplots(n_rows, 2, figsize=(13.0, 1.7 * n_rows),
                             squeeze=False)

    def _draw(ax, monthly: np.ndarray, title: str, metric_values: Dict[str, float]):
        x = np.arange(12)
        ax.fill_between(x, 0, monthly, color=COLORS["empirical"],
                        alpha=0.35, step="mid")
        ax.plot(x, monthly, color=COLORS["empirical"], lw=1.4)
        ax.set_xticks(x)
        ax.set_xticklabels(WY_MONTH_LABELS, fontsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.set_ylabel("monthly flow (cfs)", fontsize=7)
        ax.set_title(title, fontsize=8, loc="left")
        annot = "\n".join([
            f"{name}: {val:.2f}" for name, val in metric_values.items()
        ])
        ax.text(0.985, 0.95, annot, transform=ax.transAxes,
                ha="right", va="top", fontsize=6.5,
                family="monospace",
                bbox={"facecolor": "white", "alpha": 0.85,
                      "edgecolor": "lightgrey"})

    # Historical
    for i, (m_extr, wy) in enumerate(hist_picks):
        if i >= n_rows: break
        ax = axes[i, 0]
        # Pull WY j: months [j*12, j*12+12) in monthly_1d
        block = monthly_1d_full[wy * 12: (wy + 1) * 12]
        # Recover labels from the historical chars dataframe
        row = h[h["water_year_index"] == wy].iloc[0]
        mvals = {m: float(row[m]) for m in pick_metrics if m in row.index}
        _draw(ax, block,
              f"Historical WY index {wy} (extreme on {m_extr})", mvals)

    # Kirsch
    for i, (m_extr, tid) in enumerate(kir_picks):
        if i >= n_rows: break
        ax = axes[i, 1]
        # T=1 evaluation block = months 24..35 of the 36-month trace
        block = library[tid, 24:36]
        row = k.iloc[tid]
        mvals = {m: float(row[m]) for m in pick_metrics if m in row.index}
        _draw(ax, block,
              f"Kirsch trace #{tid} (extreme on {m_extr})", mvals)

    axes[0, 0].set_title("Historical (n = 72 water years)",
                         fontsize=10, loc="center", fontweight="bold")
    axes[0, 1].set_title("Kirsch ensemble (n = 10 000)",
                        fontsize=10, loc="center", fontweight="bold")

    fig.suptitle(
        "Fig 07 — exemplar T = 1 yr hydrographs (water-year ordered, "
        "Oct → Sep) selected to span extremes of djf_total, "
        "min_monthly, summer_recession, jja_total",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_dir / "fig07_exemplar_traces.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig07] wrote {out_dir / 'fig07_exemplar_traces.png'}")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main():
    out_dir = stage_figure_dir("02_calibration", "diagnostics_for_decision")
    fig01_distributions(out_dir)
    fig02_degeneracy_spread(out_dir)
    fig03_correlations(out_dir)
    fig04_scatter_T1(out_dir)
    fig05_kset_viability(out_dir)
    fig06_top_ksets(out_dir)
    fig07_exemplar_traces(out_dir)


if __name__ == "__main__":
    main()
