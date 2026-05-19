"""Diagnostic figures supporting the DD-15b short-block reframe.

Produces three PNGs under
``figures/02_calibration/diagnostics_for_decision/``:

* ``df7_short_block_distributions.{png,pdf}`` — per-metric distribution
  comparison: violins of historical 1-yr and 2-yr blocks and
  jittered-strip Kirsch ensembles, side-by-side per metric. Lets the
  reader see the marginal shapes that the K-set selection sees.

* ``df8_short_block_corr.{png,pdf}`` — pairwise |Spearman ρ|
  heatmaps. Three columns (T=1, T=2, plus the existing T=5 hazard-clean
  result for context). Two rows: historical reference and Kirsch
  ensemble. Demonstrates the per-T evolution of metric independence.

* ``df9_short_block_kset_summary.{png,pdf}`` — bar chart of the
  number of K-sets passing the strict-rung filter at varying
  correlation caps for T = 1, 2, 5. Headline: how much room does
  each T provide for a low-correlation K-set?
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.metrics.extended import HAZARD_CLEAN_METRICS  # noqa: E402
from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import COLORS, apply_style  # noqa: E402
from src.metrics.short_block import (  # noqa: E402
    SHORT_BLOCK_CONCEPT_MAP,
    SHORT_BLOCK_METRIC_NAMES,
)


def _hist_short_chars(T: int) -> Optional[pd.DataFrame]:
    p = stage_output_dir(
        "02_calibration", "short_block_screening", slug=f"T{T}",
        create=False,
    ) / "block_chars_short.csv"
    return pd.read_csv(p) if p.exists() else None


def _kirsch_short_chars(T: int, n_traces: int = 10_000, seed: int = 42) -> Optional[pd.DataFrame]:
    p = stage_output_dir(
        "03_kirsch_library", "build_short_block_library",
        slug=f"n{n_traces}_s{seed}", create=False,
    ) / f"characteristics_short_T{T}.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    return pd.DataFrame(z["values"], columns=[str(n) for n in z["metric_names"]])


def _hist_T5_chars() -> Optional[pd.DataFrame]:
    p = stage_output_dir(
        "02_calibration", "t_sensitivity_historical", slug="T05",
        create=False,
    ) / "block_chars_extended.csv"
    return pd.read_csv(p) if p.exists() else None


def _kirsch_T5_chars() -> Optional[pd.DataFrame]:
    p = stage_output_dir(
        "03_kirsch_library", "build_library_extended",
        slug="n10000_t5_ssi3-12_s42", create=False,
    ) / "characteristics_extended.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    return pd.DataFrame(z["values"], columns=[str(n) for n in z["metric_names"]])


# ---------------------------------------------------------------------------
# DF7 — per-metric distributions across T
# ---------------------------------------------------------------------------


def figure_df7(out_dir: Path):
    apply_style()
    metrics = list(SHORT_BLOCK_METRIC_NAMES)

    hist_T1 = _hist_short_chars(1)
    hist_T2 = _hist_short_chars(2)
    kir_T1 = _kirsch_short_chars(1)
    kir_T2 = _kirsch_short_chars(2)

    if hist_T1 is None:
        print("[df7] missing historical T=1 chars; skipping")
        return

    n_rows = len(metrics)
    fig, axes = plt.subplots(
        n_rows, 2, figsize=(7.5, 1.6 * n_rows),
        squeeze=False, sharex=False,
    )

    for i, m in enumerate(metrics):
        for j, (label, hist_df, kir_df) in enumerate([
            ("T = 1 yr", hist_T1, kir_T1),
            ("T = 2 yr", hist_T2, kir_T2),
        ]):
            ax = axes[i, j]
            if hist_df is None or m not in hist_df.columns:
                ax.set_visible(False)
                continue
            hv = hist_df[m].dropna().astype(float).values
            ax.violinplot(hv, positions=[0.0], widths=0.7,
                          showextrema=False, showmedians=True)
            if kir_df is not None and m in kir_df.columns:
                rng = np.random.default_rng(0)
                kv = kir_df[m].dropna().astype(float).values
                jit = rng.uniform(0.45, 1.05, size=min(kv.size, 1500))
                samp = rng.choice(kv, size=min(kv.size, 1500), replace=False)
                ax.scatter(jit, samp, s=2, alpha=0.10,
                           color=COLORS["kde"], rasterized=True)
            ax.set_xlim(-0.5, 1.4)
            ax.set_xticks([])
            ax.tick_params(axis="y", labelsize=7)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if i == 0:
                ax.set_title(label, fontsize=9)
            if j == 0:
                ax.set_ylabel(
                    f"{m}\n[{SHORT_BLOCK_CONCEPT_MAP.get(m, '?')}]",
                    fontsize=7,
                )

    fig.suptitle(
        "DF7 — per-metric distributions at T = 1 yr and T = 2 yr\n"
        "blue violin = historical blocks (Hist) | green jitter = Kirsch ensemble",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_dir / "df7_short_block_distributions.png",
                dpi=180, bbox_inches="tight")
    fig.savefig(out_dir / "df7_short_block_distributions.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"[df7] wrote {out_dir / 'df7_short_block_distributions.png'}")


# ---------------------------------------------------------------------------
# DF8 — correlation matrices across T
# ---------------------------------------------------------------------------


def figure_df8(out_dir: Path):
    apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(15.0, 9.0), squeeze=False)

    panels: List[tuple] = []
    # Column 1: T=1 short metrics
    h1 = _hist_short_chars(1)
    k1 = _kirsch_short_chars(1)
    cols1 = list(SHORT_BLOCK_METRIC_NAMES) if h1 is not None else None
    panels.append((0, "T=1 yr — short-block pool", h1, k1, cols1))
    # Column 2: T=2 short metrics
    h2 = _hist_short_chars(2)
    k2 = _kirsch_short_chars(2)
    cols2 = list(SHORT_BLOCK_METRIC_NAMES) if h2 is not None else None
    panels.append((1, "T=2 yr — short-block pool", h2, k2, cols2))
    # Column 3: T=5 hazard-clean for context
    h5 = _hist_T5_chars()
    k5 = _kirsch_T5_chars()
    cols5 = list(HAZARD_CLEAN_METRICS) if h5 is not None else None
    panels.append((2, "T=5 yr — hazard-clean (context)", h5, k5, cols5))

    summary_rows: List[Dict] = []

    for col_idx, title, hdf, kdf, cols in panels:
        for row, (src, df) in enumerate([("Historical", hdf), ("Kirsch", kdf)]):
            ax = axes[row, col_idx]
            if df is None or cols is None:
                ax.axis("off")
                ax.set_title(f"{title}\n{src}: missing", fontsize=9)
                continue
            avail = [c for c in cols if c in df.columns]
            if len(avail) < 2:
                ax.axis("off"); continue
            sub = df[avail].astype(float)
            rho = sub.corr(method="spearman")
            arr = rho.values
            im = ax.imshow(np.abs(arr), cmap="RdYlBu_r",
                           vmin=0.0, vmax=1.0, aspect="auto")
            ax.set_xticks(range(len(rho.columns)))
            ax.set_yticks(range(len(rho.columns)))
            ax.set_xticklabels(rho.columns, rotation=45, ha="right",
                               fontsize=6)
            ax.set_yticklabels(rho.columns, fontsize=6)
            for i in range(len(rho.columns)):
                for j in range(len(rho.columns)):
                    v = arr[i, j]
                    c = "white" if abs(v) >= 0.85 else "black"
                    ax.text(j, i, f"{abs(v):.2f}", ha="center", va="center",
                            fontsize=5.5, color=c)
            if row == 0:
                ax.set_title(f"{title}\n{src}", fontsize=9)
            else:
                ax.set_title(f"{src}", fontsize=9)
            off = rho.where(~np.eye(len(rho), dtype=bool)).abs().stack()
            summary_rows.append({
                "T": title.split(" ")[0], "source": src,
                "n_samples": int(len(df)),
                "n_metrics": len(rho.columns),
                "min_abs_rho": float(off.min()),
                "median_abs_rho": float(off.median()),
                "max_abs_rho": float(off.max()),
            })

    fig.suptitle(
        "DF8 — pairwise |Spearman ρ| across T (short-block pool at T=1,2; "
        "hazard-clean pool at T=5).\n"
        "Top row = historical blocks, bottom row = Kirsch ensemble",
        fontsize=10,
    )
    if im is not None:
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, pad=0.02,
                     location="right", label="|Spearman ρ|")
    fig.savefig(out_dir / "df8_short_block_corr.png", dpi=180,
                bbox_inches="tight")
    fig.savefig(out_dir / "df8_short_block_corr.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[df8] wrote {out_dir / 'df8_short_block_corr.png'}")

    pd.DataFrame(summary_rows).to_csv(
        out_dir / "df8_corr_summary.csv", index=False,
    )


# ---------------------------------------------------------------------------
# DF9 — K-set viability by correlation cap
# ---------------------------------------------------------------------------


def _count_ksets(rho: pd.DataFrame, K: int, cap: float) -> int:
    metrics = list(rho.columns)
    n = len(metrics)
    if n < K:
        return 0
    abs_rho = rho.abs().values.copy()
    np.fill_diagonal(abs_rho, 0.0)
    count = 0
    for combo in combinations(range(n), K):
        max_off = 0.0
        for i in range(K):
            for j in range(i + 1, K):
                v = abs_rho[combo[i], combo[j]]
                if v > max_off:
                    max_off = v
        if max_off < cap:
            count += 1
    return count


def figure_df9(out_dir: Path):
    apply_style()
    caps = [0.4, 0.5, 0.6, 0.7, 0.8]
    Ks = [3, 4]

    settings = []
    h1 = _hist_short_chars(1); k1 = _kirsch_short_chars(1)
    h2 = _hist_short_chars(2); k2 = _kirsch_short_chars(2)
    h5 = _hist_T5_chars(); k5 = _kirsch_T5_chars()
    if h1 is not None and k1 is not None:
        settings.append(("T=1", k1, list(SHORT_BLOCK_METRIC_NAMES)))
    if h2 is not None and k2 is not None:
        settings.append(("T=2", k2, list(SHORT_BLOCK_METRIC_NAMES)))
    if h5 is not None and k5 is not None:
        settings.append(("T=5 (HC)", k5, list(HAZARD_CLEAN_METRICS)))

    fig, axes = plt.subplots(1, len(Ks), figsize=(4.6 * len(Ks), 4.0),
                             squeeze=False)
    for ax_idx, K in enumerate(Ks):
        ax = axes[0, ax_idx]
        x = np.arange(len(caps))
        bw = 0.85 / len(settings)
        for s_idx, (label, df, cols) in enumerate(settings):
            avail = [c for c in cols if c in df.columns]
            sub = df[avail].astype(float)
            rho = sub.corr(method="spearman")
            counts = [_count_ksets(rho, K=K, cap=cap) for cap in caps]
            ax.bar(x + (s_idx - len(settings) / 2) * bw + bw / 2,
                   counts, bw, label=label)
            for i, c in enumerate(counts):
                ax.text(x[i] + (s_idx - len(settings) / 2) * bw + bw / 2,
                        c + 0.5, str(c), ha="center", va="bottom",
                        fontsize=7)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{c:.2f}" for c in caps])
        ax.set_xlabel("correlation cap (Kirsch |ρ|)")
        ax.set_ylabel(f"# strict-rung K={K} sets passing cap")
        ax.set_title(f"K = {K}", fontsize=10)
        ax.legend(fontsize=8, frameon=False)

    fig.suptitle(
        "DF9 — K-set viability by Kirsch correlation cap and T\n"
        "Higher bars = more independent K-set choices available",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_dir / "df9_short_block_kset_summary.png",
                dpi=180, bbox_inches="tight")
    fig.savefig(out_dir / "df9_short_block_kset_summary.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"[df9] wrote {out_dir / 'df9_short_block_kset_summary.png'}")


def main():
    out_dir = stage_figure_dir("02_calibration", "diagnostics_for_decision")
    figure_df7(out_dir)
    figure_df8(out_dir)
    figure_df9(out_dir)


if __name__ == "__main__":
    main()
