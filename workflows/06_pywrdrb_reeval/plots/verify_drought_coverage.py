"""Render the five Stage-A verification figures from cached compute outputs.

Reads the JSON + npz artifacts written by
``workflows/06_pywrdrb_reeval/verify_drought_coverage.py`` and emits one
PDF per criterion under
``figures/06_pywrdrb_reeval/verify_drought_coverage/<src_slug>/``.

Plotting-only driver -- it never re-runs the verification logic.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Sequence

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from src.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "06_pywrdrb_reeval"
DRIVER = "verify_drought_coverage"


def _save_fig(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[plots/06/verify] wrote {path}")


def plot_c1_coverage(
    pareto_drought_metrics: np.ndarray,
    hist_block_chars: np.ndarray,
    objective_keys: Sequence[str],
    anti_ideal: np.ndarray,
    output_path: Path,
) -> None:
    K = pareto_drought_metrics.shape[1]
    n_panels = K * (K - 1) // 2 if K >= 2 else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 4.5), squeeze=False)
    axes = axes.flatten()

    p = 0
    for i in range(K):
        for j in range(i + 1, K):
            ax = axes[p]
            ax.scatter(hist_block_chars[:, i], hist_block_chars[:, j],
                       s=40, c="#1f77b4", alpha=0.6, label="Historical blocks")
            ax.scatter(pareto_drought_metrics[:, i], pareto_drought_metrics[:, j],
                       s=8, c="#ff7f0e", alpha=0.5, label="Pareto archive")
            if anti_ideal.size >= K:
                ax.scatter([anti_ideal[i]], [anti_ideal[j]],
                           marker="X", s=180, c="red", zorder=5, label="D*")
            ax.set_xlabel(objective_keys[i])
            ax.set_ylabel(objective_keys[j])
            ax.legend(fontsize=8, loc="best")
            p += 1
    fig.suptitle("Criterion 1 -- Drought-space coverage", fontsize=12)
    fig.tight_layout()
    _save_fig(fig, output_path)


def plot_c2_low_flow(report: Dict, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    if report.get("status") == "insufficient_data":
        ax.text(0.5, 0.5,
                f"insufficient_data\nn_drought_subset={report.get('n_drought_subset', 0)}",
                ha="center", va="center", transform=ax.transAxes)
        _save_fig(fig, output_path)
        return

    xs = [r["exceedance"] * 100 for r in report["per_exceedance"]]
    syn_p50 = [r["syn_p50"] for r in report["per_exceedance"]]
    hist_p50 = [r["hist_p50"] for r in report["per_exceedance"]]
    syn_p10 = [r["syn_p10"] for r in report["per_exceedance"]]
    hist_p10 = [r["hist_p10"] for r in report["per_exceedance"]]

    ax.plot(xs, hist_p50, "-o", color="#1f77b4", label="Historical p50")
    ax.plot(xs, hist_p10, "--", color="#1f77b4", label="Historical p10")
    ax.plot(xs, syn_p50, "-o", color="#ff7f0e", label="Drought-subset p50")
    ax.plot(xs, syn_p10, "--", color="#ff7f0e", label="Drought-subset p10")

    for r in report["per_exceedance"]:
        if not r["pass"]:
            ax.axvspan(r["exceedance"] * 100 - 1, r["exceedance"] * 100 + 1,
                       alpha=0.12, color="red")

    ax.set_yscale("log")
    ax.set_xlabel("Exceedance (%)")
    ax.set_ylabel("Flow (cfs)")
    ax.set_title(f"Criterion 2 -- Low-flow directionality ({report['status']})")
    ax.legend(fontsize=9)
    fig.tight_layout()
    _save_fig(fig, output_path)


def plot_c3_nominal(report: Dict, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    status = report.get("status", "?")

    if status == "insufficient_data":
        for ax in axes:
            ax.text(0.5, 0.5, f"insufficient_data\nn={report.get('n_nominal_subset', 0)}",
                    ha="center", va="center", transform=ax.transAxes)
        fig.suptitle(f"Criterion 3 -- Nominal-subset fidelity ({status})")
        fig.tight_layout()
        _save_fig(fig, output_path)
        return

    fdc = report["fdc_high_flow"]
    axes[0].bar([str(e) for e in fdc["exceedances"]],
                fdc["pct_inside_per_exceedance"], color="#2ca02c")
    axes[0].axhline(90, color="red", linestyle="--", label="90% target")
    axes[0].set_ylabel("% of nominal traces inside envelope")
    axes[0].set_title(f"High-flow FDC containment ({'pass' if fdc['pass'] else 'fail'})")
    axes[0].legend()

    ac = report["lag1_ac"]
    axes[1].barh(["Synthetic", "Historical"],
                 [ac["syn_5_95"][1] - ac["syn_5_95"][0],
                  ac["hist_5_95"][1] - ac["hist_5_95"][0]],
                 left=[ac["syn_5_95"][0], ac["hist_5_95"][0]],
                 color=["#ff7f0e", "#1f77b4"], alpha=0.7)
    axes[1].set_xlabel("lag-1 AC")
    axes[1].set_title(f"Lag-1 AC 5-95 interval ({'pass' if ac['intervals_overlap'] else 'fail'})")

    sc = report["seasonal_cycle"]
    axes[2].bar(["Monthly mean", "Monthly std"],
                [sc["min_month_pct_inside_mean"], sc["min_month_pct_inside_std"]],
                color=["#2ca02c", "#2ca02c"])
    axes[2].axhline(90, color="red", linestyle="--", label="90% target")
    axes[2].set_ylabel("Worst-month % inside envelope")
    axes[2].set_title(f"Seasonal cycle ({'pass' if sc['pass'] else 'fail'})")
    axes[2].legend()

    fig.suptitle(f"Criterion 3 -- Nominal-subset fidelity ({status})")
    fig.tight_layout()
    _save_fig(fig, output_path)


def plot_c4_spread(
    A_syn: np.ndarray,
    A_hist: np.ndarray,
    drought_metrics: np.ndarray,
    drought_mask: np.ndarray,
    nominal_mask: np.ndarray,
    report: Dict,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    bins = 30

    axes[0].hist(A_hist, bins=bins, alpha=0.5, color="#1f77b4",
                 label=f"Historical blocks (std={report['std_hist']:.0f})")
    axes[0].hist(A_syn, bins=bins, alpha=0.5, color="#ff7f0e",
                 label=f"Pareto archive (std={report['std_syn']:.0f})")
    axes[0].set_xlabel("Annual mean flow (cfs-month)")
    axes[0].set_ylabel("Count")
    spread_pass_str = "pass" if report["spread_pass"] else "fail"
    axes[0].set_title(
        f"Spread ratio = {report['spread_ratio']:.2f} [0.9, 1.1] -- {spread_pass_str}"
    )
    axes[0].legend()

    axes[1].scatter(drought_metrics[:, 0], A_syn,
                    s=6, c="#888888", alpha=0.3, label="All Pareto")
    if nominal_mask.any():
        axes[1].scatter(drought_metrics[nominal_mask, 0], A_syn[nominal_mask],
                        s=12, c="#2ca02c", alpha=0.6,
                        label=f"Nominal corner (n={int(nominal_mask.sum())})")
    if drought_mask.any():
        axes[1].scatter(drought_metrics[drought_mask, 0], A_syn[drought_mask],
                        s=12, c="#d62728", alpha=0.5,
                        label=f"Drought corner (n={int(drought_mask.sum())})")
    hist_med = report["corner_bifurcation"]["hist_median_annual_mean"]
    axes[1].axhline(hist_med, color="black", linestyle="--",
                    label=f"Historical median = {hist_med:.0f}")
    axes[1].set_xlabel("mean_duration")
    axes[1].set_ylabel("Annual mean flow")
    axes[1].set_title("Corner bifurcation")
    axes[1].legend(fontsize=8)

    fig.suptitle(f"Criterion 4 -- Flow-space spread ({report['status']})")
    fig.tight_layout()
    _save_fig(fig, output_path)


def plot_c5_onset(report: Dict, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    if report.get("status") == "insufficient_data":
        ax.text(0.5, 0.5, "insufficient_data",
                ha="center", va="center", transform=ax.transAxes)
        _save_fig(fig, output_path)
        return

    months = np.arange(1, 13)
    syn = np.array(report["synthetic_histogram"], dtype=float)
    hist = np.array(report["historical_histogram"], dtype=float)
    syn_n = syn / max(syn.sum(), 1.0)
    hist_n = hist / max(hist.sum(), 1.0)

    width = 0.4
    ax.bar(months - width / 2, hist_n, width, color="#1f77b4",
           label=f"Historical (n={int(hist.sum())})")
    ax.bar(months + width / 2, syn_n, width, color="#ff7f0e",
           label=f"Synthetic (n={int(syn.sum())})")

    ax.set_xticks(months)
    ax.set_xlabel("Calendar month")
    ax.set_ylabel("Fraction of drought events")
    ax.set_title(
        f"Criterion 5 -- Drought onset seasonality "
        f"(chi2={report['chi2_statistic']:.1f}, p={report['p_value']:.3f} "
        f"-> {report['status']})"
    )
    ax.legend()
    fig.tight_layout()
    _save_fig(fig, output_path)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True,
                   help="Source slug under outputs/06_pywrdrb_reeval/"
                        "verify_drought_coverage/")
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    report_json = in_dir / "verification_report.json"
    if not report_json.exists():
        sys.exit(f"missing {report_json} -- run the compute driver first")

    report = json.loads(report_json.read_text())
    objective_keys = tuple(report["objective_keys"])

    subsets = np.load(in_dir / "subsets.npz", allow_pickle=True)
    drought_mask = subsets["drought_mask"].astype(bool)
    nominal_mask = subsets["nominal_mask"].astype(bool)
    drought_metrics = subsets["drought_metrics"]
    anti_ideal = subsets["anti_ideal"]

    chars_npz = np.load(in_dir / "hist_block_chars.npz", allow_pickle=True)
    hist_block_chars = chars_npz["chars"]

    annual = np.load(in_dir / "annual_means.npz")
    A_syn = annual["A_syn"]
    A_hist = annual["A_hist"]

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    crit = report["criteria"]
    plot_c1_coverage(drought_metrics, hist_block_chars,
                     objective_keys, anti_ideal,
                     fig_dir / "c1_drought_coverage.pdf")
    plot_c2_low_flow(crit["c2_low_flow"], fig_dir / "c2_drought_subset_fdc.pdf")
    plot_c3_nominal(crit["c3_nominal_fidelity"], fig_dir / "c3_nominal_subset.pdf")
    plot_c4_spread(A_syn, A_hist, drought_metrics,
                   drought_mask, nominal_mask, crit["c4_spread"],
                   fig_dir / "c4_annual_spread.pdf")
    plot_c5_onset(crit["c5_drought_onset"], fig_dir / "c5_drought_onset.pdf")
    print(f"[plots/06/verify] all figures written under {fig_dir}")


if __name__ == "__main__":
    main()
