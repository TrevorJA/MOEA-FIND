"""SI-quality validity diagnostics for the Mapping G transformation.

Mapping G converts a hydrologic window summary ``s_j`` into a bounded
drought score ``D_g = Φ((s_j − μ_log) / σ_log)`` using historical
``μ_log``, ``σ_log``. The percentile interpretation of ``D_g`` rests
on ``s_j`` being approximately Gaussian under the historical record.
This script tests that assumption per window and writes a single
multi-panel figure suitable for a manuscript supporting-information
section, plus a JSON of test statistics.

Per window, four diagnostic panels:

  1. **Sample histogram + fitted Gaussian PDF**: visual fit check.
  2. **QQ plot vs normal**: tail-departure diagnostic.
  3. **Empirical vs fitted Gaussian CDF**: calibration plot. The
     supremum gap is the one-sample KS statistic against ``N(μ, σ)``.
  4. **G vs E mapping agreement**: scatter of ``D_g`` against ``D_e``
     for each historical year on a 1:1 reference. Tail divergence is
     where Mapping G saturates relative to the empirical-CDF mapping.

Reported statistics per window:
  - Shapiro-Wilk ``W`` and ``p_W``
  - Anderson-Darling ``A²`` (one-sample, normal) and critical value at α=5%
  - One-sample KS distance ``D_KS`` (empirical vs fitted Gaussian)
  - Sample mean / std / skewness / excess kurtosis
  - max ``|D_g − D_e|`` across historical years

Generic across any preset whose objective keys are bounded windows
produced by ``compute_candidate_bounded_metrics``. The script auto-
discovers the four (or however many) ``<window>_g`` keys in the
archive's ``results.json`` and runs the validity checks per window.
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

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "run_moea_find"


def _window_name_from_metric(metric_key: str) -> str | None:
    """Strip the ``_g`` / ``_e`` suffix to recover the window name."""
    for suffix in ("_g", "_e"):
        if metric_key.endswith(suffix):
            return metric_key[: -len(suffix)]
    return None


def historical_window_samples(
    monthly_1d: np.ndarray, window_name: str
) -> np.ndarray:
    """Per-water-year ``s_j`` samples for one window.

    Reproduces the construction in
    :class:`src.metrics.extended.BoundedFamilyRefs.from_full_record`
    so the samples returned here are the exact distribution Mapping G
    is fit on.
    """
    from src.metrics.short_block import WINDOW_SPECS, _compute_window_summary

    if window_name not in WINDOW_SPECS:
        raise KeyError(f"unknown window {window_name!r}")
    kind, params = WINDOW_SPECS[window_name]
    flows = np.asarray(monthly_1d, dtype=float).ravel()
    n_full_yrs = flows.size // 12
    flows_2d = flows[: n_full_yrs * 12].reshape(n_full_yrs, 12)
    out = np.empty(n_full_yrs, dtype=float)
    for y in range(n_full_yrs):
        out[y] = _compute_window_summary(flows_2d[y : y + 1], kind, params)
    return out[np.isfinite(out)]


def per_window_stats(samples: np.ndarray) -> dict:
    """Return all summary stats and normality tests for one window."""
    s = np.asarray(samples, dtype=float)
    n = s.size
    mu = float(s.mean())
    sigma = float(s.std(ddof=1)) if n > 1 else 0.0

    # Shapiro-Wilk
    sw_w, sw_p = stats.shapiro(s) if 3 <= n <= 5000 else (np.nan, np.nan)

    # Anderson-Darling (one-sample, normal)
    ad = stats.anderson(s, dist="norm")
    a2 = float(ad.statistic)
    crit_5 = float(ad.critical_values[2])  # 5% alpha is index 2

    # KS distance vs fitted Gaussian (parameters estimated from sample;
    # this is the Lilliefors variant — exact null is non-Gaussian, but
    # the statistic itself is interpretable as a calibration gap).
    cdf_emp = np.arange(1, n + 1) / n
    s_sorted = np.sort(s)
    cdf_fit = stats.norm.cdf(s_sorted, loc=mu, scale=sigma) if sigma > 0 else np.full(n, 0.5)
    ks = float(np.max(np.abs(cdf_emp - cdf_fit)))

    return {
        "n": int(n),
        "mu": mu,
        "sigma": sigma,
        "skew": float(stats.skew(s, bias=False)) if n > 2 else float("nan"),
        "excess_kurtosis": float(stats.kurtosis(s, bias=False, fisher=True)) if n > 3 else float("nan"),
        "shapiro_W": float(sw_w),
        "shapiro_p": float(sw_p),
        "anderson_A2": a2,
        "anderson_crit_5pct": crit_5,
        "ks_lilliefors": ks,
    }


def _g_e_for_each_sample(samples: np.ndarray, mu: float, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    from src.metrics.short_block import _apply_mapping_e, _apply_mapping_g
    arr_sorted = np.sort(samples)
    g = np.array([_apply_mapping_g(float(s_j), mu, sigma) for s_j in samples])
    e = np.array([_apply_mapping_e(float(s_j), arr_sorted) for s_j in samples])
    return g, e


def render_validity_grid(
    window_data: list[dict],
    out_path: Path,
    title: str,
) -> None:
    """4 panels (rows) × N_windows (cols) validity grid."""
    n_w = len(window_data)
    fig, axes = plt.subplots(4, n_w, figsize=(3.6 * n_w, 11.0))
    if n_w == 1:
        axes = axes[:, None]

    for j, w in enumerate(window_data):
        s = w["samples"]
        mu, sigma = w["stats"]["mu"], w["stats"]["sigma"]
        st = w["stats"]
        wname = w["window"]

        # Row 0: histogram + Gaussian PDF
        ax = axes[0, j]
        ax.hist(s, bins=14, density=True, color="#bbb", edgecolor="white",
                alpha=0.85, label=f"hist (n={st['n']})")
        x_grid = np.linspace(s.min() - 0.5 * sigma, s.max() + 0.5 * sigma, 200)
        ax.plot(x_grid, stats.norm.pdf(x_grid, loc=mu, scale=sigma),
                color="#d62728", lw=2, label=f"N(μ={mu:.2f}, σ={sigma:.2f})")
        ax.set_title(f"{wname}\nSW W={st['shapiro_W']:.3f} p={st['shapiro_p']:.3f}",
                     fontsize=9)
        ax.set_xlabel(r"$s_j$ (drought-positive)", fontsize=8)
        if j == 0:
            ax.set_ylabel("density", fontsize=9)
        ax.legend(fontsize=7)
        ax.tick_params(labelsize=7)

        # Row 1: QQ plot vs Normal
        ax = axes[1, j]
        (osm, osr), (slope, intercept, r) = stats.probplot(s, dist="norm", plot=None)
        ax.scatter(osm, osr, s=18, color="#2b6cb0", alpha=0.85)
        line_y = slope * osm + intercept
        ax.plot(osm, line_y, color="#d62728", lw=1.4,
                label=f"r²={r**2:.3f}")
        ax.set_title(f"AD A²={st['anderson_A2']:.2f} (crit5%={st['anderson_crit_5pct']:.2f})",
                     fontsize=9)
        ax.set_xlabel("theoretical quantile", fontsize=8)
        if j == 0:
            ax.set_ylabel("sample quantile", fontsize=9)
        ax.legend(fontsize=7)
        ax.tick_params(labelsize=7)

        # Row 2: empirical vs fitted Gaussian CDF
        ax = axes[2, j]
        cdf_emp = np.arange(1, st["n"] + 1) / st["n"]
        s_sorted = np.sort(s)
        cdf_fit = stats.norm.cdf(s_sorted, loc=mu, scale=sigma)
        ax.plot(s_sorted, cdf_emp, drawstyle="steps-post",
                color="#2b6cb0", lw=1.6, label="empirical")
        ax.plot(s_sorted, cdf_fit, color="#d62728", lw=1.4, ls="--",
                label="N(μ, σ)")
        # KS gap
        gap_idx = int(np.argmax(np.abs(cdf_emp - cdf_fit)))
        ax.vlines(s_sorted[gap_idx], cdf_emp[gap_idx], cdf_fit[gap_idx],
                  color="#16a085", lw=2.4, label=f"KS={st['ks_lilliefors']:.3f}")
        ax.set_title("calibration: F_emp vs F_G", fontsize=9)
        ax.set_xlabel(r"$s_j$", fontsize=8)
        if j == 0:
            ax.set_ylabel("CDF", fontsize=9)
        ax.legend(fontsize=7, loc="lower right")
        ax.tick_params(labelsize=7)

        # Row 3: G vs E agreement
        ax = axes[3, j]
        g, e = _g_e_for_each_sample(s, mu, sigma)
        ax.plot([0, 1], [0, 1], color="0.6", lw=1, ls="--")
        ax.scatter(g, e, s=22, color="#9467bd", alpha=0.85,
                   edgecolor="white", linewidth=0.4)
        max_div = float(np.max(np.abs(g - e)))
        ax.set_title(f"G vs E agreement\nmax|D_g − D_e|={max_div:.3f}",
                     fontsize=9)
        ax.set_xlabel(r"$D_g = \Phi(z)$", fontsize=8)
        if j == 0:
            ax.set_ylabel(r"$D_e$ (empirical CDF)", fontsize=9)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.tick_params(labelsize=7)
        # Stash for JSON output
        w["max_g_minus_e"] = max_div

    fig.suptitle(title, fontsize=12, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True,
                   help="Variant slug under outputs/04_moea_find_single_site/run_moea_find/")
    p.add_argument("--windows", nargs="*", default=None,
                   help="Optional explicit window names (without _g suffix). "
                        "Default: auto-discover from the archive's objective_keys.")
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

    print(f"[mapping_g_validity] slug={args.slug}")
    print(f"  windows: {windows}")

    from src.experiment import prepare_data
    cache_dir = PROJECT_ROOT / "data" / "usgs_cache"
    _, monthly_1d = prepare_data(cache_dir)

    window_data = []
    for w in windows:
        samples = historical_window_samples(monthly_1d, w)
        st = per_window_stats(samples)
        window_data.append({"window": w, "samples": samples, "stats": st})
        print(f"  {w:24s} n={st['n']:3d}  μ={st['mu']:.3f}  σ={st['sigma']:.3f}  "
              f"SW p={st['shapiro_p']:.3f}  AD A²={st['anderson_A2']:.2f}  "
              f"KS={st['ks_lilliefors']:.3f}")

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    cfg = json.loads((in_dir / "config.json").read_text())
    title = (
        f"Mapping G validity — {cfg.get('metric_set','?')} preset, slug {args.slug}\n"
        f"K={len(obj_keys)} windows × 4 panels: PDF fit · QQ · CDF calibration · G vs E"
    )
    render_validity_grid(window_data, fig_dir / "mapping_g_validity.png", title)
    print(f"  wrote {fig_dir / 'mapping_g_validity.png'}")
    print(f"  wrote {fig_dir / 'mapping_g_validity.pdf'}")

    # Statistics JSON for the SI table
    summary = {
        "slug": args.slug,
        "metric_set": cfg.get("metric_set"),
        "n_years": cfg.get("n_years"),
        "windows": [
            {
                "window": w["window"],
                **w["stats"],
                "max_g_minus_e": float(w.get("max_g_minus_e", float("nan"))),
            }
            for w in window_data
        ],
    }
    json_path = fig_dir / "mapping_g_validity.json"
    json_path.write_text(json.dumps(summary, indent=2))
    print(f"  wrote {json_path}")


if __name__ == "__main__":
    main()
