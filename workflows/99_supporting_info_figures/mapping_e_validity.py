"""SI validity diagnostic for Mapping E (empirical CDF with polynomial-
tail extrapolation), structured to mirror ``mapping_g_validity.py``
but with panels appropriate to E's design.

E makes no parametric distribution claim, so the four panels test
different properties than the G version:

  1. **Sample histogram + KDE**: what E is implicitly approximating
     (the empirical density). No fit overlay — E *is* the empirical
     CDF inside the envelope.

  2. **Leave-one-out predictive check**: for each historical year y,
     refit ``BoundedFamilyRefs`` on the remaining n−1 years and predict
     ``D_e(s_y)``. Plot LOO-predicted vs that year's true Hazen
     percentile in the full sample. Bias / RMSE quantify how stable E
     is to dropping a single year — a direct test of leverage.

  3. **Body vs tail behavior**: ``D_e(s)`` over a grid extending past
     the envelope ``[s_min, s_max]``, with the empirical envelope
     marked. Inside the envelope D_e is the linear interpolation of
     Hazen plotting positions; outside it is the Cauchy-shaped
     polynomial tail anchored to the local empirical slope.

  4. **E vs G mapping agreement**: per-historical-year ``D_e`` vs
     ``D_g`` on a 1:1 reference — same panel as G's validity figure
     row 4 with axes swapped, kept for visual parallelism.

Reported statistics per window:
  - LOO bias, RMSE, max |residual|
  - In-envelope max ``|D_e − F_Hazen|`` (zero by construction; reported
    to confirm)
  - Out-of-envelope tail slope at top/bottom k_tail samples
  - max ``|D_g − D_e|`` across historical years
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
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "run_moea_find"


def _window_name_from_metric(metric_key: str) -> str | None:
    for suffix in ("_g", "_e"):
        if metric_key.endswith(suffix):
            return metric_key[: -len(suffix)]
    return None


def historical_window_samples(monthly_1d: np.ndarray, window_name: str) -> np.ndarray:
    from src.metrics.short_block import WINDOW_SPECS, _compute_window_summary

    if window_name not in WINDOW_SPECS:
        raise KeyError(f"unknown window {window_name!r}")
    kind, params = WINDOW_SPECS[window_name]
    flows = np.asarray(monthly_1d, dtype=float).ravel()
    n = flows.size // 12
    flows_2d = flows[: n * 12].reshape(n, 12)
    out = np.array([_compute_window_summary(flows_2d[y : y + 1], kind, params)
                    for y in range(n)])
    return out[np.isfinite(out)]


def _hazen_position(arr_sorted: np.ndarray, s: float) -> float:
    """Hazen plotting position of s in arr_sorted by linear interpolation."""
    n = arr_sorted.size
    ranks = np.arange(1, n + 1, dtype=float) / (n + 1.0)
    if s <= arr_sorted[0]:
        return float(ranks[0])
    if s >= arr_sorted[-1]:
        return float(ranks[-1])
    return float(np.interp(s, arr_sorted, ranks))


def loo_predict(samples: np.ndarray) -> np.ndarray:
    """Leave-one-out D_e prediction for each sample.

    For each ``y``, refit ``sorted_hist`` on the remaining ``n−1`` years
    and apply Mapping E to get the predicted ``D_e(s_y)``.
    """
    from src.metrics.short_block import _apply_mapping_e

    n = samples.size
    out = np.empty(n, dtype=float)
    for i in range(n):
        held = np.delete(samples, i)
        out[i] = _apply_mapping_e(float(samples[i]), np.sort(held))
    return out


def per_window_stats(samples: np.ndarray, loo_pred: np.ndarray) -> dict:
    n = samples.size
    arr_sorted = np.sort(samples)
    ranks = np.arange(1, n + 1, dtype=float) / (n + 1.0)
    # In-envelope sanity: D_e on the training points equals Hazen rank by
    # construction — report the max residual to confirm.
    from src.metrics.short_block import _apply_mapping_e
    in_envelope = np.array([_apply_mapping_e(float(s), arr_sorted) for s in samples])
    rank_of = np.array([_hazen_position(arr_sorted, float(s)) for s in samples])
    in_max = float(np.max(np.abs(in_envelope - rank_of)))

    # LOO residuals against the year's full-sample Hazen rank
    loo_resid = loo_pred - rank_of
    loo_bias = float(loo_resid.mean())
    loo_rmse = float(np.sqrt(np.mean(loo_resid ** 2)))
    loo_max = float(np.max(np.abs(loo_resid)))

    # Tail slopes (per Mapping E's k_tail=5 anchor)
    k_tail = min(5, n)
    if k_tail >= 2:
        slope_low = float((ranks[k_tail - 1] - ranks[0])
                          / max(arr_sorted[k_tail - 1] - arr_sorted[0], 1e-12))
        slope_high = float((ranks[-1] - ranks[-k_tail])
                           / max(arr_sorted[-1] - arr_sorted[-k_tail], 1e-12))
    else:
        slope_low = slope_high = float("nan")

    return {
        "n": int(n),
        "mu": float(samples.mean()),
        "sigma": float(samples.std(ddof=1)),
        "in_envelope_max_resid": in_max,
        "loo_bias": loo_bias,
        "loo_rmse": loo_rmse,
        "loo_max_abs": loo_max,
        "tail_slope_low": slope_low,
        "tail_slope_high": slope_high,
    }


def _g_e_for_each_sample(samples: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    from src.metrics.short_block import _apply_mapping_e, _apply_mapping_g
    arr_sorted = np.sort(samples)
    mu = float(samples.mean())
    sigma = float(samples.std(ddof=1))
    g = np.array([_apply_mapping_g(float(s), mu, sigma) for s in samples])
    e = np.array([_apply_mapping_e(float(s), arr_sorted) for s in samples])
    return g, e


def render_e_validity_grid(
    window_data: list[dict], out_path: Path, title: str
) -> None:
    from src.metrics.short_block import _apply_mapping_e

    n_w = len(window_data)
    fig, axes = plt.subplots(4, n_w, figsize=(3.6 * n_w, 11.0))
    if n_w == 1:
        axes = axes[:, None]

    for j, w in enumerate(window_data):
        s = w["samples"]
        st = w["stats"]
        wname = w["window"]
        arr_sorted = np.sort(s)
        ranks = np.arange(1, st["n"] + 1, dtype=float) / (st["n"] + 1.0)

        # Row 0: histogram + KDE of s_j
        ax = axes[0, j]
        ax.hist(s, bins=14, density=True, color="#bbb", edgecolor="white",
                alpha=0.85, label=f"hist (n={st['n']})")
        try:
            kde = stats.gaussian_kde(s)
            x_grid = np.linspace(s.min() - 0.5 * st["sigma"],
                                 s.max() + 0.5 * st["sigma"], 200)
            ax.plot(x_grid, kde(x_grid), color="#9467bd", lw=2,
                    label="empirical KDE")
        except Exception:
            pass
        ax.set_title(f"{wname}\nμ={st['mu']:.2f}, σ={st['sigma']:.2f}",
                     fontsize=9)
        ax.set_xlabel(r"$s_j$ (drought-positive)", fontsize=8)
        if j == 0:
            ax.set_ylabel("density", fontsize=9)
        ax.legend(fontsize=7)
        ax.tick_params(labelsize=7)

        # Row 1: LOO predicted D_e vs full-sample Hazen rank
        ax = axes[1, j]
        loo = w["loo_pred"]
        rank_of = np.array([_hazen_position(arr_sorted, float(v)) for v in s])
        ax.plot([0, 1], [0, 1], color="0.6", lw=1, ls="--")
        ax.scatter(rank_of, loo, s=18, color="#2b6cb0", alpha=0.85,
                   edgecolor="white", linewidth=0.4)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_title(f"LOO RMSE={st['loo_rmse']:.3f}, max|res|={st['loo_max_abs']:.3f}",
                     fontsize=9)
        ax.set_xlabel("full-sample Hazen rank", fontsize=8)
        if j == 0:
            ax.set_ylabel(r"$D_e$ (LOO predict)", fontsize=9)
        ax.tick_params(labelsize=7)

        # Row 2: D_e(s) over a grid extending past the empirical envelope
        ax = axes[2, j]
        s_min, s_max = float(arr_sorted[0]), float(arr_sorted[-1])
        pad = 0.6 * st["sigma"]
        s_grid = np.linspace(s_min - pad, s_max + pad, 240)
        d_e = np.array([_apply_mapping_e(float(x), arr_sorted) for x in s_grid])
        # Inside-envelope segment vs outside-envelope tails
        in_env = (s_grid >= s_min) & (s_grid <= s_max)
        ax.plot(s_grid[in_env], d_e[in_env], color="#2b6cb0", lw=1.8,
                label="inside envelope")
        out_low = s_grid < s_min
        out_hi = s_grid > s_max
        if out_low.any():
            ax.plot(s_grid[out_low], d_e[out_low], color="#d62728", lw=1.8,
                    label="lower tail (poly)")
        if out_hi.any():
            ax.plot(s_grid[out_hi], d_e[out_hi], color="#d62728", lw=1.8,
                    label="upper tail (poly)" if not out_low.any() else None)
        ax.axvline(s_min, color="0.4", lw=0.8, ls=":")
        ax.axvline(s_max, color="0.4", lw=0.8, ls=":")
        # Empirical CDF underlay for context
        cdf_emp = np.arange(1, st["n"] + 1) / st["n"]
        ax.plot(arr_sorted, cdf_emp, drawstyle="steps-post",
                color="0.55", lw=1.0, alpha=0.85, label="F_emp")
        ax.set_xlabel(r"$s_j$", fontsize=8)
        if j == 0:
            ax.set_ylabel(r"$D_e$", fontsize=9)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(f"slope_low={st['tail_slope_low']:.2g}, "
                     f"slope_high={st['tail_slope_high']:.2g}",
                     fontsize=9)
        ax.legend(fontsize=7, loc="lower right")
        ax.tick_params(labelsize=7)

        # Row 3: D_e vs D_g agreement (mirror of G version with swapped axes)
        ax = axes[3, j]
        g, e = _g_e_for_each_sample(s)
        ax.plot([0, 1], [0, 1], color="0.6", lw=1, ls="--")
        ax.scatter(e, g, s=22, color="#9467bd", alpha=0.85,
                   edgecolor="white", linewidth=0.4)
        max_div = float(np.max(np.abs(g - e)))
        ax.set_title(f"E vs G agreement\nmax|D_e − D_g|={max_div:.3f}",
                     fontsize=9)
        ax.set_xlabel(r"$D_e$ (empirical CDF)", fontsize=8)
        if j == 0:
            ax.set_ylabel(r"$D_g = \Phi(z)$", fontsize=9)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.tick_params(labelsize=7)
        w["max_e_minus_g"] = max_div

    fig.suptitle(title, fontsize=12, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True)
    p.add_argument("--windows", nargs="*", default=None)
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    result = json.loads((in_dir / "results.json").read_text())
    obj_keys = list(result["objective_keys"])
    if args.windows:
        windows = list(args.windows)
    else:
        windows = []
        for k in obj_keys:
            w = _window_name_from_metric(k)
            if w is not None and w not in windows:
                windows.append(w)
        if not windows:
            sys.exit(f"no _g/_e suffixed metrics in objective_keys: {obj_keys}")

    print(f"[mapping_e_validity] slug={args.slug}")
    print(f"  windows: {windows}")

    from src.experiment import prepare_data
    cache_dir = PROJECT_ROOT / "data" / "usgs_cache"
    _, monthly_1d = prepare_data(cache_dir)

    window_data = []
    for w in windows:
        samples = historical_window_samples(monthly_1d, w)
        loo_pred = loo_predict(samples)
        st = per_window_stats(samples, loo_pred)
        window_data.append({"window": w, "samples": samples,
                            "loo_pred": loo_pred, "stats": st})
        print(f"  {w:24s} n={st['n']:3d}  LOO RMSE={st['loo_rmse']:.3f}  "
              f"max|res|={st['loo_max_abs']:.3f}  in-env={st['in_envelope_max_resid']:.1e}  "
              f"slope_low={st['tail_slope_low']:.2g} slope_high={st['tail_slope_high']:.2g}")

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    cfg = json.loads((in_dir / "config.json").read_text())
    title = (
        f"Mapping E validity — {cfg.get('metric_set','?')} preset, slug {args.slug}\n"
        f"K={len(obj_keys)} windows × 4 panels: KDE · LOO predict · body vs tail · E vs G"
    )
    out_png = fig_dir / "mapping_e_validity.png"
    render_e_validity_grid(window_data, out_png, title)
    print(f"  wrote {out_png}")
    print(f"  wrote {out_png.with_suffix('.pdf')}")

    summary = {
        "slug": args.slug,
        "metric_set": cfg.get("metric_set"),
        "n_years": cfg.get("n_years"),
        "windows": [
            {"window": w["window"], **w["stats"],
             "max_e_minus_g": float(w.get("max_e_minus_g", float("nan")))}
            for w in window_data
        ],
    }
    json_path = fig_dir / "mapping_e_validity.json"
    json_path.write_text(json.dumps(summary, indent=2))
    print(f"  wrote {json_path}")


if __name__ == "__main__":
    main()
