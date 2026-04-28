"""Script 16 — Kirsch wrapper-mode comparison and SI figure set.

Loads every ``results.json`` written by ``15_wrapper_mode_ablation.py`` under
``outputs/exp15_wrapper_mode_ablation/{index,residual}/...``, pools the Pareto
archive across seeds within each mode, and renders the SI figures:

    figSI_wrapper_pareto_2d.pdf         — side-by-side drought-space Pareto
    figSI_wrapper_pareto_3d.pdf         — same in 3D (if >=3 objectives)
    figSI_wrapper_manhattan_dist.pdf    — Manhattan-objective distribution
    figSI_wrapper_hydrology_panels.pdf  — ACF / FDC / seasonal panel
    figSI_wrapper_per_trace_stats.pdf   — per-trace stat boxplots
    figSI_wrapper_dv_distributions.pdf  — DV QQ + histogram per mode
    figSI_wrapper_dv_tail_mass.pdf      — DV tail-mass comparison

Also writes ``wrapper_comparison_summary.json`` with the comparison checklist
values (max-duration ratio, max-severity ratio, median Manhattan distance,
hyperplane identity, infeasibility rates) so the check is machine-readable.

Both modes must produce at least one Pareto archive under
``outputs/exp15_wrapper_mode_ablation/`` or the script warns and emits
whatever figures are possible.

Run serially (no MPI, no SLURM required for local execution):
    python workflows/04_moea_find_single_site/wrapper_mode_compare.py

Or via SLURM:
    sbatch workflows/04_moea_find_single_site/slurm/wrapper_mode_compare.slurm
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment_utils import (  # noqa: E402
    prepare_data,
    compute_historical_ssi_chars,
    compute_ssi_anti_ideal,
)
from src.experiment_config import DEFAULT_EXPERIMENT  # noqa: E402

INPUT_SLUG = "exp15_wrapper_mode_ablation"
OUTPUT_SLUG = "exp16_wrapper_mode_compare"

DEFAULT_ARMS: Tuple[str, ...] = ("index", "residual")

_LOGICAL_ARM_SPEC: Dict[str, Dict[str, Any]] = {
    "index":    {"disk_arm": "index",    "require_statistic": None},
    "residual": {"disk_arm": "residual", "require_statistic": None},
}


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _iter_arm_result_files(root: Path, disk_arm: str) -> List[Path]:
    arm_dir = root / disk_arm
    if not arm_dir.exists():
        return []
    return sorted(arm_dir.glob("*/results.json"))


def _load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text())


# ---------------------------------------------------------------------------
# Per-trace statistic helpers (parallel to src/constraints.py primitives but
# returning raw magnitudes rather than deviations, since we want to compare
# each mode to the historical distribution directly).
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
    p.add_argument("--arms", nargs="+", default=list(DEFAULT_ARMS),
                   help="Which wrapper modes to compare. Options: "
                        "index, residual. Default is both.")
    p.add_argument("--kirsch-library-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / "exp05_kirsch_library",
                   help="Directory containing exp05 Kirsch library "
                        "characteristics.npz, used for the reference panel.")
    args = p.parse_args()

    unknown = [a for a in args.arms if a not in _LOGICAL_ARM_SPEC]
    if unknown:
        raise SystemExit(
            f"Unknown arm(s): {unknown}. "
            f"Valid: {list(_LOGICAL_ARM_SPEC)}"
        )
    ARMS: Tuple[str, ...] = tuple(args.arms)

    out_dir = args.output_dir
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load per-arm archives
    # ------------------------------------------------------------------
    arm_results: Dict[str, List[Dict[str, Any]]] = {a: [] for a in ARMS}
    for arm in ARMS:
        spec = _LOGICAL_ARM_SPEC[arm]
        disk_arm = spec["disk_arm"]
        files = _iter_arm_result_files(args.input_dir, disk_arm)
        print(f"[16] arm={arm} (disk_arm={disk_arm}) found "
              f"{len(files)} results.json files")
        for f in files:
            try:
                r = _load_json(f)
            except Exception as exc:
                print(f"[16] WARNING: failed to load {f}: {exc}")
                continue
            nfe = int(r.get("nfe", 0))
            if nfe < args.min_nfe:
                print(f"[16]   skipping {f} "
                      f"(nfe={nfe} < min_nfe={args.min_nfe})")
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

    print(f"[16] pool summary:")
    for a in ARMS:
        print(f"     {a}: {len(pooled_traces_1d[a])} Pareto traces "
              f"(across {len(arm_results[a])} seeds)")

    # ------------------------------------------------------------------
    # Historical baseline (same for both modes)
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

    # Median T-block point for plotting (T-block scale, not full-record).
    hist_block_median = np.median(hist_block_chars, axis=0)

    # Kirsch library reference cloud from exp05 (random unconstrained traces)
    kirsch_cloud_2d: Optional[np.ndarray] = None
    kirsch_cloud_3d: Optional[np.ndarray] = None
    kirsch_lib_path = args.kirsch_library_dir / "characteristics.npz"
    if kirsch_lib_path.exists():
        try:
            klib = np.load(kirsch_lib_path, allow_pickle=True)
            k_all_keys = [str(k) for k in klib["all_keys"]]
            k_all_vals = klib["all_values"]
            obj_indices = []
            for ok in objective_keys:
                if ok in k_all_keys:
                    obj_indices.append(k_all_keys.index(ok))
                else:
                    obj_indices.append(None)
            if all(j is not None for j in obj_indices):
                kirsch_chars = k_all_vals[:, obj_indices]
                kirsch_cloud_2d = kirsch_chars[:, :2]
                kirsch_cloud_3d = kirsch_chars[:, :3] if len(obj_indices) >= 3 else None
                print(f"[16] loaded Kirsch library: {k_all_vals.shape[0]} traces")
            else:
                missing = [k for k, j in zip(objective_keys, obj_indices) if j is None]
                print(f"[16] WARNING: Kirsch library missing keys {missing}, "
                      "skipping reference panel")
        except Exception as exc:
            print(f"[16] WARNING: could not load Kirsch library: {exc}")
    else:
        print(f"[16] Kirsch library not found at {kirsch_lib_path}, "
              "reference panel will show only historical T-blocks")

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
    # Comparison checklist
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

    # Per-arm range summary for every arm that produced Pareto solutions.
    arm_ranges: Dict[str, Dict[str, Dict[str, float]]] = {}
    for a in ARMS:
        if pooled_metrics[a].size == 0:
            continue
        arm_ranges[a] = {
            k: _range(a, j) for j, k in enumerate(objective_keys)
        }
    checklist["ranges"] = arm_ranges

    # Pairwise max-ratio: residual / index (if both modes have Pareto archives).
    if "index" in arm_ranges and "residual" in arm_ranges:
        idx_ranges = arm_ranges["index"]
        ratios: Dict[str, float] = {
            k: _ratio(arm_ranges["residual"][k]["max"], idx_ranges[k]["max"])
            for k in objective_keys
        }
        checklist["max_ratio_residual_over_index"] = ratios

    # Median Manhattan per arm + delta residual - index.
    manhattan_median: Dict[str, float] = {}
    for a in ARMS:
        if pooled_manhattan[a].size:
            manhattan_median[a] = float(np.median(pooled_manhattan[a]))
    checklist["manhattan_median"] = manhattan_median
    if "index" in manhattan_median and "residual" in manhattan_median:
        checklist["manhattan_median_delta_residual_over_index"] = (
            manhattan_median["residual"] - manhattan_median["index"]
        )

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

    (out_dir / "wrapper_comparison_summary.json").write_text(
        json.dumps(checklist, indent=2, default=str)
    )
    print(f"[16] wrote {out_dir / 'wrapper_comparison_summary.json'}")

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from src.plotting.ablation_dv import (
            ARM_COLORS,
            ARM_LABELS,
            plot_hydrology_panels_two_arms,
            plot_manhattan_distribution,
            plot_dv_distributions,
            plot_dv_tail_mass,
            plot_pareto_2d_subpanels,
            plot_pareto_3d_subpanels,
            plot_per_trace_stats,
        )
    except ImportError as exc:
        print(f"[16] ABORT: required plotting modules missing: {exc}")
        return

    # ---- figSI_wrapper_pareto_2d.pdf (multi-panel, shared axes) ----
    if any(pooled_metrics[a].size for a in ARMS):
        pareto_by_arm_2d = {
            a: pooled_metrics[a][:, :2] if pooled_metrics[a].size else np.zeros((0, 2))
            for a in ARMS
        }
        fig, _ = plot_pareto_2d_subpanels(
            pareto_by_arm=pareto_by_arm_2d,
            historical_cloud=hist_block_chars[:, :2],
            historical_median=tuple(hist_block_median[:2]),
            kirsch_cloud=kirsch_cloud_2d,
            objective_labels=("mean duration (months)", "mean avg severity"),
            arm_order=list(ARMS),
        )
        fig.savefig(fig_dir / "figSI_wrapper_pareto_2d.pdf", dpi=200,
                    bbox_inches="tight")
        plt.close(fig)
        print(f"[16] wrote {fig_dir / 'figSI_wrapper_pareto_2d.pdf'}")

    # ---- figSI_wrapper_pareto_3d.pdf (multi-panel 3D, shared limits) ----
    if (any(pooled_metrics[a].size for a in ARMS)
            and len(objective_keys) >= 3):
        k3 = tuple(objective_keys[:3])
        any_arm = next(a for a in ARMS if pooled_metrics[a].size)
        if pooled_metrics[any_arm].shape[1] >= 3:
            pareto_by_arm_3d = {
                a: (pooled_metrics[a][:, :3] if pooled_metrics[a].size
                    else np.zeros((0, 3)))
                for a in ARMS
            }
            fig3, _ = plot_pareto_3d_subpanels(
                pareto_by_arm=pareto_by_arm_3d,
                historical_cloud=hist_block_chars[:, :3],
                historical_median=hist_block_median[:3],
                kirsch_cloud=kirsch_cloud_3d,
                objective_labels=k3,
                arm_order=list(ARMS),
            )
            fig3.savefig(fig_dir / "figSI_wrapper_pareto_3d.pdf",
                         dpi=180, bbox_inches="tight")
            plt.close(fig3)
            print(f"[16] wrote {fig_dir / 'figSI_wrapper_pareto_3d.pdf'}")

    # ---- figSI_wrapper_manhattan_dist.pdf ----
    if any(pooled_manhattan[a].size for a in ARMS):
        fig_m, _ = plot_manhattan_distribution(pooled_manhattan)
        fig_m.savefig(fig_dir / "figSI_wrapper_manhattan_dist.pdf", dpi=200,
                      bbox_inches="tight")
        plt.close(fig_m)
        print(f"[16] wrote {fig_dir / 'figSI_wrapper_manhattan_dist.pdf'}")

    # ---- figSI_wrapper_hydrology_panels.pdf ----
    if any(pooled_traces_1d[a] for a in ARMS):
        fig_h, _ = plot_hydrology_panels_two_arms(
            traces_1d_by_arm={a: pooled_traces_1d[a] for a in ARMS},
            traces_2d_by_arm={a: pooled_traces_2d[a] for a in ARMS},
            historical_blocks_1d=hist_blocks_1d,
            historical_blocks_2d=hist_blocks_2d,
        )
        fig_h.savefig(fig_dir / "figSI_wrapper_hydrology_panels.pdf",
                      dpi=200, bbox_inches="tight")
        plt.close(fig_h)
        print(f"[16] wrote {fig_dir / 'figSI_wrapper_hydrology_panels.pdf'}")

    # ---- figSI_wrapper_per_trace_stats.pdf ----
    fig_b, _ = plot_per_trace_stats(
        stats_by_arm=stats_by_arm,
        hist_stats=hist_per_trace_stats,
        stat_order=["annual_mean", "annual_cv", "lag1_ac",
                    "seasonal_max_frac_dev"],
    )
    fig_b.savefig(fig_dir / "figSI_wrapper_per_trace_stats.pdf",
                  dpi=200, bbox_inches="tight")
    plt.close(fig_b)
    print(f"[16] wrote {fig_dir / 'figSI_wrapper_per_trace_stats.pdf'}")

    # ---- figSI_wrapper_dv_distributions.pdf ----
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
        fig_d.savefig(fig_dir / "figSI_wrapper_dv_distributions.pdf",
                      dpi=200, bbox_inches="tight")
        plt.close(fig_d)
        print(f"[16] wrote {fig_dir / 'figSI_wrapper_dv_distributions.pdf'}")

        fig_tm, _ = plot_dv_tail_mass(dvs_sampled, tail_bounds=(0.05, 0.95))
        fig_tm.savefig(fig_dir / "figSI_wrapper_dv_tail_mass.pdf",
                       dpi=200, bbox_inches="tight")
        plt.close(fig_tm)
        print(f"[16] wrote {fig_dir / 'figSI_wrapper_dv_tail_mass.pdf'}")

    print("[16] done.")


if __name__ == "__main__":
    main()
