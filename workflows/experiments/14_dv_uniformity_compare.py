"""Script 14 — DV-uniformity ablation comparison and SI figure set.

Loads every ``results.json`` written by ``13_dv_uniformity_ablation.py`` under
``outputs/exp13_dv_uniformity_ablation/{hydrologic,dv_uniform}/...``, pools
the Pareto archive across seeds within each arm, and renders the SI figures:

    figSI_ablation_pareto_2d.pdf         — side-by-side drought-space Pareto
    figSI_ablation_pareto_3d.pdf         — same in 3D (if >=3 objectives)
    figSI_ablation_manhattan_dist.pdf    — Manhattan-objective distribution
    figSI_ablation_hydrology_panels.pdf  — ACF / FDC / seasonal panel
    figSI_ablation_per_arm_timeseries_stats.pdf  — per-trace stat boxplots
    figSI_ablation_dv_distributions.pdf  — DV QQ + histogram per arm
    figSI_ablation_infeasibility_bar.pdf — final infeasibility per arm

Also writes ``comparison_summary.json`` with the falsification checklist
values (max-duration ratio, max-severity ratio, median Manhattan distance,
hyperplane identity, infeasibility rates) so the check is machine-readable.

Both arms must produce at least one Pareto archive under
``outputs/exp13_dv_uniformity_ablation/`` or the script warns and emits
whatever figures are possible.

Run serially (no MPI, no SLURM required for local execution):
    python workflows/experiments/14_dv_uniformity_compare.py

Or via SLURM:
    sbatch workflows/slurm/14_dv_uniformity_compare.slurm
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment_utils import (  # noqa: E402
    prepare_data,
    compute_historical_ssi_chars,
    compute_ssi_anti_ideal,
)
from src.experiment_config import DEFAULT_EXPERIMENT  # noqa: E402

INPUT_SLUG = "exp13_dv_uniformity_ablation"
OUTPUT_SLUG = "exp14_dv_uniformity_compare"
ARMS = ("hydrologic", "dv_uniform")


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _iter_arm_result_files(root: Path, arm: str) -> List[Path]:
    arm_dir = root / arm
    if not arm_dir.exists():
        return []
    return sorted(arm_dir.glob("*/results.json"))


def _load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text())


# ---------------------------------------------------------------------------
# Per-trace statistic helpers (parallel to src/constraints.py primitives but
# returning raw magnitudes rather than deviations, since we want to compare
# each arm to the historical distribution directly).
# ---------------------------------------------------------------------------

def _annual_mean(trace: np.ndarray) -> float:
    n = (len(trace) // 12) * 12
    if n == 0:
        return float("nan")
    return float(trace[:n].reshape(-1, 12).sum(axis=1).mean())


def _annual_cv(trace: np.ndarray) -> float:
    n = (len(trace) // 12) * 12
    if n == 0:
        return float("nan")
    totals = trace[:n].reshape(-1, 12).sum(axis=1)
    mean = totals.mean()
    if mean <= 0:
        return float("nan")
    return float(totals.std(ddof=1) / mean)


def _lag1_ac(trace: np.ndarray) -> float:
    if len(trace) < 3:
        return float("nan")
    return float(np.corrcoef(trace[:-1], trace[1:])[0, 1])


def _seasonal_max_frac_dev(trace_2d: np.ndarray,
                            hist_monthly_means: np.ndarray) -> float:
    cycle = np.asarray(trace_2d).mean(axis=0)
    hist = np.asarray(hist_monthly_means, dtype=float)
    safe_hist = np.where(hist > 0, hist, 1.0)
    return float(np.max(np.abs(cycle / safe_hist - 1.0)))


def _compute_per_trace_stats(
    traces_1d: List[np.ndarray],
    traces_2d: List[np.ndarray],
    hist_monthly_means: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Compute the summary statistics used for the boxplot panel."""
    ann_mean = np.array([_annual_mean(t) for t in traces_1d])
    ann_cv = np.array([_annual_cv(t) for t in traces_1d])
    lag1 = np.array([_lag1_ac(t) for t in traces_1d])
    seas = np.array([
        _seasonal_max_frac_dev(t2, hist_monthly_means) for t2 in traces_2d
    ])
    return {
        "annual_mean": ann_mean,
        "annual_cv": ann_cv,
        "lag1_ac": lag1,
        "seasonal_max_frac_dev": seas,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--input-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / INPUT_SLUG)
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    p.add_argument("--max-dv-rows-per-arm", type=int, default=200,
                   help="Cap on number of Pareto-member DV vectors sampled "
                        "for the QQ/histogram figure.")
    p.add_argument("--min-nfe", type=int, default=0,
                   help="Ignore any results.json whose 'nfe' is below this "
                        "threshold. Use to exclude smoke-test runs from the "
                        "production comparison pool.")
    args = p.parse_args()

    out_dir = args.output_dir
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load per-arm archives
    # ------------------------------------------------------------------
    arm_results: Dict[str, List[Dict[str, Any]]] = {a: [] for a in ARMS}
    for arm in ARMS:
        files = _iter_arm_result_files(args.input_dir, arm)
        print(f"[14] arm={arm} found {len(files)} results.json files")
        for f in files:
            try:
                r = _load_json(f)
            except Exception as exc:
                print(f"[14] WARNING: failed to load {f}: {exc}")
                continue
            nfe = int(r.get("nfe", 0))
            if nfe < args.min_nfe:
                print(f"[14]   skipping {f} (nfe={nfe} < min_nfe={args.min_nfe})")
                continue
            arm_results[arm].append(r)

    # ------------------------------------------------------------------
    # Pool per arm
    # ------------------------------------------------------------------
    def _pool(arm: str, key: str) -> List:
        out: List = []
        for r in arm_results[arm]:
            out.extend(r.get(key, []))
        return out

    pooled_traces_1d = {a: [np.asarray(t, dtype=float)
                             for t in _pool(a, "pareto_traces_1d")]
                        for a in ARMS}
    pooled_traces_2d = {a: [np.asarray(t, dtype=float)
                             for t in _pool(a, "pareto_traces_2d")]
                        for a in ARMS}
    pooled_metrics = {a: np.array(_pool(a, "drought_metrics"))
                      if _pool(a, "drought_metrics") else np.zeros((0, 0))
                      for a in ARMS}
    pooled_dvs = {a: np.array(_pool(a, "pareto_dvs"))
                  if _pool(a, "pareto_dvs") else np.zeros((0, 0))
                  for a in ARMS}
    # Manhattan objective = last column of the full Pareto objective vector
    # (objectives = drought_metrics + Manhattan). We don't have the raw
    # objective vector persisted directly, so recompute from drought_metrics
    # using the anti-ideal stored in results.json.
    def _manhattan_from_result(r: Dict[str, Any]) -> np.ndarray:
        dm = np.asarray(r.get("drought_metrics", []), dtype=float)
        ai = np.asarray(r.get("anti_ideal", []), dtype=float)
        if dm.size == 0 or ai.size == 0:
            return np.array([])
        return np.sum(ai[None, :] - dm, axis=1)

    pooled_manhattan: Dict[str, np.ndarray] = {a: np.array([]) for a in ARMS}
    for a in ARMS:
        parts = [_manhattan_from_result(r) for r in arm_results[a]]
        parts = [p for p in parts if p.size > 0]
        if parts:
            pooled_manhattan[a] = np.concatenate(parts)

    print(f"[14] pool summary:")
    for a in ARMS:
        print(f"     {a}: {len(pooled_traces_1d[a])} Pareto traces "
              f"(across {len(arm_results[a])} seeds)")

    # ------------------------------------------------------------------
    # Historical baseline (same for both arms)
    # ------------------------------------------------------------------
    cfg = DEFAULT_EXPERIMENT
    n_years = cfg.n_years_out
    ssi_timescale = cfg.ssi_timescale
    objective_keys = cfg.objective_keys

    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)
    _, ssi_calc, hist_chars = compute_historical_ssi_chars(
        monthly_1d, ssi_timescale,
    )
    anti_ideal = compute_ssi_anti_ideal(
        hist_chars, objective_keys, headroom=cfg.anti_ideal_headroom,
    )

    from src.historical_blocks import (
        resample_historical_blocks,
        resample_historical_blocks_2d,
        compute_historical_block_chars,
    )
    hist_blocks_1d = resample_historical_blocks(
        monthly_1d, T_years=n_years, stride=1,
    )
    hist_blocks_2d = resample_historical_blocks_2d(
        monthly_2d, T_years=n_years, stride=1,
    )
    hist_monthly_means = monthly_2d.mean(axis=0)
    hist_per_trace_stats = _compute_per_trace_stats(
        hist_blocks_1d, hist_blocks_2d, hist_monthly_means,
    )
    hist_block_chars = compute_historical_block_chars(
        monthly_1d, T_years=n_years, ssi_calc=ssi_calc,
        objective_keys=objective_keys, stride=1,
    )

    # ------------------------------------------------------------------
    # Per-arm stats
    # ------------------------------------------------------------------
    stats_by_arm: Dict[str, Dict[str, np.ndarray]] = {}
    for a in ARMS:
        if pooled_traces_1d[a]:
            stats_by_arm[a] = _compute_per_trace_stats(
                pooled_traces_1d[a], pooled_traces_2d[a], hist_monthly_means,
            )
        else:
            stats_by_arm[a] = {k: np.array([])
                               for k in hist_per_trace_stats}

    # ------------------------------------------------------------------
    # Falsification checklist
    # ------------------------------------------------------------------
    def _ratio(arm_val, ref_val):
        if ref_val == 0 or not np.isfinite(ref_val):
            return float("nan")
        return float(arm_val / ref_val)

    checklist: Dict[str, Any] = {
        "arms_seen": {a: len(arm_results[a]) for a in ARMS},
        "n_pareto_pooled": {a: len(pooled_traces_1d[a]) for a in ARMS},
    }

    def _range(a: str, j: int) -> Dict[str, float]:
        if pooled_metrics[a].size == 0:
            return {"min": float("nan"), "max": float("nan")}
        return {
            "min": float(pooled_metrics[a][:, j].min()),
            "max": float(pooled_metrics[a][:, j].max()),
        }

    if pooled_metrics["hydrologic"].size and pooled_metrics["dv_uniform"].size:
        ranges_hy = [_range("hydrologic", j) for j in range(pooled_metrics["hydrologic"].shape[1])]
        ranges_dv = [_range("dv_uniform", j) for j in range(pooled_metrics["dv_uniform"].shape[1])]
        checklist["ranges"] = {
            "hydrologic": {k: r for k, r in zip(objective_keys, ranges_hy)},
            "dv_uniform": {k: r for k, r in zip(objective_keys, ranges_dv)},
        }
        checklist["max_ratio_dv_over_hydrologic"] = {
            objective_keys[j]: _ratio(ranges_dv[j]["max"], ranges_hy[j]["max"])
            for j in range(len(ranges_hy))
        }
        checklist["falsified_by_max_duration_lt_80pct"] = bool(
            checklist["max_ratio_dv_over_hydrologic"].get("mean_duration", 1.0) < 0.80
        )
        if "mean_avg_severity" in checklist["max_ratio_dv_over_hydrologic"]:
            checklist["falsified_by_max_avg_severity_lt_80pct"] = bool(
                checklist["max_ratio_dv_over_hydrologic"]["mean_avg_severity"] < 0.80
            )

    if pooled_manhattan["hydrologic"].size and pooled_manhattan["dv_uniform"].size:
        med_hy = float(np.median(pooled_manhattan["hydrologic"]))
        med_dv = float(np.median(pooled_manhattan["dv_uniform"]))
        checklist["manhattan_median"] = {"hydrologic": med_hy, "dv_uniform": med_dv}
        checklist["manhattan_median_delta"] = med_dv - med_hy
        checklist["falsified_by_manhattan_median"] = bool(med_dv > med_hy)

    # Infeasibility rates per arm (pooled across seeds)
    infeasibility_rates: Dict[str, float] = {}
    for a in ARMS:
        infeas = sum(r.get("n_infeasible", 0) for r in arm_results[a])
        total = sum(r.get("n_evals_total", 0) for r in arm_results[a])
        infeasibility_rates[a] = (infeas / total) if total > 0 else 0.0
    checklist["infeasibility_rate"] = infeasibility_rates

    # Hyperplane identity (one value per seed)
    hp: Dict[str, List[Dict[str, float]]] = {}
    for a in ARMS:
        hp[a] = [
            {
                "expected_sum": float(r["hyperplane"]["expected_sum"]),
                "actual_mean": float(r["hyperplane"]["actual_mean"]),
                "actual_std": float(r["hyperplane"]["actual_std"]),
            }
            for r in arm_results[a] if "hyperplane" in r
        ]
    checklist["hyperplane"] = hp

    (out_dir / "comparison_summary.json").write_text(
        json.dumps(checklist, indent=2, default=str)
    )
    print(f"[14] wrote {out_dir / 'comparison_summary.json'}")

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from src.plotting.drought_space import (
            plot_scatter_with_marginals,
            plot_drought_space_3d,
        )
        from src.plotting.ablation_dv import (
            ARM_COLORS,
            plot_hydrology_panels_two_arms,
            plot_manhattan_distribution,
            plot_dv_distributions,
            plot_per_trace_stats,
            plot_infeasibility_bar,
        )
    except ImportError as exc:
        print(f"[14] ABORT: required plotting modules missing: {exc}")
        return

    # ---- figSI_ablation_pareto_2d.pdf ----
    if any(pooled_metrics[a].size for a in ARMS):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
        hist_point = (
            float(hist_chars.get("mean_duration", 0.0)),
            float(hist_chars.get("mean_avg_severity", 0.0)),
        )
        for ax, arm in zip(axes, ARMS):
            dm = pooled_metrics[arm]
            if dm.size == 0:
                ax.text(0.5, 0.5, f"{arm}: no solutions",
                        ha="center", va="center", transform=ax.transAxes)
                ax.set_title(arm)
                continue
            ax.scatter(hist_block_chars[:, 0], hist_block_chars[:, 1],
                       s=15, alpha=0.3, color="gray",
                       label=f"historical blocks (n={len(hist_block_chars)})")
            ax.scatter(dm[:, 0], dm[:, 1], s=25, alpha=0.7,
                       color=ARM_COLORS[arm],
                       label=f"{arm} Pareto (n={dm.shape[0]})")
            ax.scatter(*hist_point, marker="*", s=220, color="black",
                       zorder=6, label="historical")
            ax.set_xlabel("mean duration (months)")
            ax.set_ylabel("mean avg severity")
            ax.set_title(arm)
            ax.legend(fontsize=7, loc="lower right")
        fig.suptitle("Pareto drought-space comparison (pooled across seeds)")
        fig.tight_layout()
        fig.savefig(fig_dir / "figSI_ablation_pareto_2d.pdf", dpi=200,
                    bbox_inches="tight")
        plt.close(fig)
        print(f"[14] wrote {fig_dir / 'figSI_ablation_pareto_2d.pdf'}")

    # ---- figSI_ablation_pareto_3d.pdf ----
    if (all(pooled_metrics[a].size for a in ARMS)
            and pooled_metrics["hydrologic"].shape[1] >= 3
            and len(objective_keys) >= 3):
        for arm in ARMS:
            fig3 = plot_drought_space_3d(
                pooled_metrics[arm][:, :3],
                anti_ideal=anti_ideal[:3],
                objective_labels=tuple(objective_keys[:3]),
                historical_point=np.array([
                    float(hist_chars.get(k, 0.0)) for k in objective_keys[:3]
                ]),
                historical_cloud=hist_block_chars[:, :3],
                title=f"{arm} Pareto (pooled)",
            )
            fig3.savefig(fig_dir / f"figSI_ablation_pareto_3d_{arm}.pdf",
                         dpi=180, bbox_inches="tight")
            plt.close(fig3)
            print(f"[14] wrote {fig_dir / f'figSI_ablation_pareto_3d_{arm}.pdf'}")

    # ---- figSI_ablation_manhattan_dist.pdf ----
    if any(pooled_manhattan[a].size for a in ARMS):
        fig_m, _ = plot_manhattan_distribution(pooled_manhattan)
        fig_m.savefig(fig_dir / "figSI_ablation_manhattan_dist.pdf", dpi=200,
                      bbox_inches="tight")
        plt.close(fig_m)
        print(f"[14] wrote {fig_dir / 'figSI_ablation_manhattan_dist.pdf'}")

    # ---- figSI_ablation_hydrology_panels.pdf ----
    if any(pooled_traces_1d[a] for a in ARMS):
        fig_h, _ = plot_hydrology_panels_two_arms(
            traces_1d_by_arm={a: pooled_traces_1d[a] for a in ARMS},
            traces_2d_by_arm={a: pooled_traces_2d[a] for a in ARMS},
            historical_blocks_1d=hist_blocks_1d,
            historical_blocks_2d=hist_blocks_2d,
        )
        fig_h.savefig(fig_dir / "figSI_ablation_hydrology_panels.pdf",
                      dpi=200, bbox_inches="tight")
        plt.close(fig_h)
        print(f"[14] wrote {fig_dir / 'figSI_ablation_hydrology_panels.pdf'}")

    # ---- figSI_ablation_per_arm_timeseries_stats.pdf ----
    fig_b, _ = plot_per_trace_stats(
        stats_by_arm=stats_by_arm,
        hist_stats=hist_per_trace_stats,
        stat_order=["annual_mean", "annual_cv", "lag1_ac",
                    "seasonal_max_frac_dev"],
    )
    fig_b.savefig(fig_dir / "figSI_ablation_per_arm_timeseries_stats.pdf",
                  dpi=200, bbox_inches="tight")
    plt.close(fig_b)
    print(f"[14] wrote "
          f"{fig_dir / 'figSI_ablation_per_arm_timeseries_stats.pdf'}")

    # ---- figSI_ablation_dv_distributions.pdf ----
    dvs_sampled: Dict[str, np.ndarray] = {}
    for a in ARMS:
        if pooled_dvs[a].size:
            arr = pooled_dvs[a]
            if arr.ndim == 2 and arr.shape[0] > args.max_dv_rows_per_arm:
                rng = np.random.default_rng(0)
                pick = rng.choice(arr.shape[0],
                                   size=args.max_dv_rows_per_arm,
                                   replace=False)
                arr = arr[pick]
            dvs_sampled[a] = arr
    if dvs_sampled:
        fig_d, _ = plot_dv_distributions(dvs_sampled)
        fig_d.savefig(fig_dir / "figSI_ablation_dv_distributions.pdf",
                      dpi=200, bbox_inches="tight")
        plt.close(fig_d)
        print(f"[14] wrote {fig_dir / 'figSI_ablation_dv_distributions.pdf'}")

    # ---- figSI_ablation_infeasibility_bar.pdf ----
    fig_ib, _ = plot_infeasibility_bar(infeasibility_rates)
    fig_ib.savefig(fig_dir / "figSI_ablation_infeasibility_bar.pdf",
                   dpi=200, bbox_inches="tight")
    plt.close(fig_ib)
    print(f"[14] wrote {fig_dir / 'figSI_ablation_infeasibility_bar.pdf'}")

    print("[14] done.")


if __name__ == "__main__":
    main()
