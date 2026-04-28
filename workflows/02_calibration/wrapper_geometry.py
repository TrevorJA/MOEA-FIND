"""Diagnostic: Kirsch-wrapper mapping geometry.

Illustrates the discrete-vs-continuous nature of the two decision-variable
mappings in src/kirsch_wrapper.py (index vs residual). Sweeps a single DV
coordinate over [0, 1] with all other DVs held fixed, and plots the resulting
flow response for each mode.

Outputs:
    outputs/diag_kirsch_wrapper_geometry/
      ├── figures/
      │   ├── fig_mapping_smoothness.pdf
      │   └── fig_mapping_cartoon.pdf
      └── geometry_summary.json   (step counts, sweep stats)

Run (SLURM only, not login-node):
    sbatch workflows/02_calibration/slurm/wrapper_geometry.slurm
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.plotting.style import COLORS, WATER_YEAR_MONTHS, apply_style  # noqa: E402

# ---------------------------------------------------------------------------
# Mode display config
# ---------------------------------------------------------------------------
MODE_COLORS = {"index": "#2b6cb0", "residual": "#c05621"}
MODE_LABELS = {"index": "Index wrapper", "residual": "Residual wrapper"}

# The swept coordinate: we pick coordinate k = n_dvs_residual // 2, which
# falls in the middle of the DV vector.  For residual mode n_dvs = n_years * 12
# and for index mode n_dvs = (n_years + 1) * 12.  We sweep the SAME positional
# index k in both wrappers.  Because residual DVs are laid out as
# dvs.reshape(n_years_out, 12), coordinate k maps to month = k % 12.
SWEEP_K_FRACTION = 0.5   # place k at the midpoint of the residual DV vector


# ---------------------------------------------------------------------------
# Data / generator helpers
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Sweep logic
# ---------------------------------------------------------------------------

def run_sweep(wrapper, k: int, n_sweep: int, rng: np.random.Generator) -> dict:
    """Sweep DV coordinate k over [0,1] and record flow responses.

    Args:
        wrapper: KirschBorgWrapper instance (already built and fitted).
        k: Index of the DV coordinate to sweep.
        n_sweep: Number of evenly-spaced points across [0, 1].
        rng: NumPy random generator (unused here; baseline is fixed at 0.5).

    Returns:
        dict with keys:
          - p_grid: sweep values (length n_sweep)
          - flow_at_k: flow for the month that coordinate k controls,
                       averaged over the years generated (length n_sweep)
          - annual_totals: 2-D array (n_sweep, n_years_out) of annual totals
          - annual_mean: per-sweep mean annual total (length n_sweep)
          - distinct_count: number of distinct flow_at_k values
          - k: the swept coordinate
          - month_idx: calendar month that k controls (k % 12)
    """
    from src.kirsch_wrapper import KirschBorgWrapper  # noqa: F401 — already imported in caller

    n_dvs = wrapper.n_dvs
    n_years = wrapper.n_years_out

    # Neutral baseline: all DVs fixed at 0.5
    v0 = np.full(n_dvs, 0.5)

    # The month that coordinate k controls (calendar-year order inside wrapper)
    month_idx = k % 12

    p_grid = np.linspace(0.0, 1.0, n_sweep)
    flow_at_k = np.zeros(n_sweep)
    annual_totals = np.zeros((n_sweep, n_years))

    for i, p in enumerate(p_grid):
        v = v0.copy()
        v[k] = p
        synthetic = wrapper.generate(v)  # shape (n_years, 12) water-year order
        # synthetic columns are water-year order (Oct=0 … Sep=11)
        # but month_idx is in terms of the DV layout (calendar: Jan=0..Dec=11)
        # The generate() method rolls by +3, so calendar month m_cal -> water col:
        #   water_col = (m_cal + 3) % 12   (Jan->3, Feb->4, … Oct->0, Nov->1, Dec->2)
        water_col = (month_idx + 3) % 12
        # Record average flow across all years for that month column
        flow_at_k[i] = synthetic[:, water_col].mean()
        annual_totals[i, :] = synthetic.sum(axis=1)

    distinct_vals = len(np.unique(np.round(flow_at_k, decimals=6)))
    return {
        "p_grid": p_grid,
        "flow_at_k": flow_at_k,
        "annual_totals": annual_totals,
        "annual_mean": annual_totals.mean(axis=1),
        "distinct_count": distinct_vals,
        "k": k,
        "month_idx": month_idx,
    }


# ---------------------------------------------------------------------------
# Smoothness metric helpers
# ---------------------------------------------------------------------------

def _flow_change_stats(p_grid: np.ndarray, flow: np.ndarray) -> dict:
    """Compute mean and max flow change per unit DV step."""
    dp = np.diff(p_grid)
    df = np.abs(np.diff(flow))
    rate = df / np.where(dp > 0, dp, np.nan)
    return {
        "mean_change_per_unit_dv": float(np.nanmean(rate)),
        "max_change_per_unit_dv": float(np.nanmax(rate)),
    }


# ---------------------------------------------------------------------------
# Figure 1: mapping smoothness
# ---------------------------------------------------------------------------

def plot_mapping_smoothness(
    results: dict,
    output_path: Path,
    n_years_hist: int,
) -> None:
    """Four-panel figure showing the DV-sweep flow response for each mode.

    Panels:
      TL: swept DV → monthly flow (both modes)
      TR: swept DV → annual total (both modes)
      BL: 12-month flow profile for three representative DV values
      BR: bar chart of distinct-value counts per mode
    """
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax_tl, ax_tr, ax_bl, ax_br = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    rep_dvs = [0.1, 0.5, 0.9]
    rep_styles = ["-", "--", ":"]

    # --- Top-left: monthly flow vs swept DV ---
    for mode, res in results.items():
        ax_tl.plot(
            res["p_grid"],
            res["flow_at_k"],
            color=MODE_COLORS[mode],
            label=MODE_LABELS[mode],
            lw=1.4,
            alpha=0.85,
        )
    ax_tl.set_xlabel("Decision variable value $p_k$")
    ax_tl.set_ylabel("Mean monthly flow (cfs)")
    ax_tl.set_title("(a) Monthly flow response to swept DV")
    ax_tl.legend(frameon=False)
    ax_tl.annotate(
        f"Hist. years = {n_years_hist}",
        xy=(0.97, 0.05), xycoords="axes fraction",
        ha="right", fontsize=8, color=COLORS["muted"],
    )

    # --- Top-right: annual total vs swept DV ---
    for mode, res in results.items():
        ax_tr.plot(
            res["p_grid"],
            res["annual_mean"],
            color=MODE_COLORS[mode],
            label=MODE_LABELS[mode],
            lw=1.4,
            alpha=0.85,
        )
    ax_tr.set_xlabel("Decision variable value $p_k$")
    ax_tr.set_ylabel("Mean annual total (cfs·mo)")
    ax_tr.set_title("(b) Annual total response to swept DV")
    ax_tr.legend(frameon=False)

    # --- Bottom-left: 12-month profile for three DV values ---
    # Use the index-mode wrapper to get full 12-month profiles for representative
    # DV values; show both modes with the same linestyle scheme but different colors.
    month_lbls = WATER_YEAR_MONTHS  # Oct … Sep
    x_months = np.arange(12)

    for s_idx, (p_rep, ls) in enumerate(zip(rep_dvs, rep_styles)):
        for mode, wrapper_info in results.items():
            # Find the nearest sweep point
            nearest_i = int(np.argmin(np.abs(wrapper_info["p_grid"] - p_rep)))
            # Annual totals already recorded; we need per-month profile.
            # The representative-year profiles were stored during sweep.
            # We stored only summary stats; re-access via the stored annual_totals.
            # We pick the first year of the trace as the representative year.
            # (The yearly variation is captured by annual_totals; the 12-month
            # seasonal shape comes from wrapper_info["profiles"] if stored.)
            # Since we did not cache per-year profiles, use annual_totals directly
            # as a proxy (single-value proxy: mean annual / 12 per month) — the
            # cartoon already shows the full shape. For the actual month-by-month
            # pattern we surface the raw sweep record we stored in profiles_by_rep.
            pass  # handled below via stored profiles

    # Re-run for the three representative points to get full 12-month profiles
    # This is fast (only 3×2 = 6 generate() calls total).
    for mode, wrapper_ref in [("index", results["index"]["_wrapper"]),
                               ("residual", results["residual"]["_wrapper"])]:
        n_dvs = wrapper_ref.n_dvs
        n_years = wrapper_ref.n_years_out
        k = results[mode]["k"]
        v0 = np.full(n_dvs, 0.5)
        for s_idx, (p_rep, ls) in enumerate(zip(rep_dvs, rep_styles)):
            v = v0.copy()
            v[k] = p_rep
            syn = wrapper_ref.generate(v)  # (n_years, 12) water-year order
            # Use median year profile
            year_totals = syn.sum(axis=1)
            med_yr = int(np.argsort(year_totals)[n_years // 2])
            profile = syn[med_yr, :]  # 12 water-year months

            label = f"$p_k={p_rep}$" if mode == "index" else None
            ax_bl.plot(
                x_months,
                profile,
                color=MODE_COLORS[mode],
                ls=ls,
                lw=1.3,
                label=label,
                alpha=0.85,
            )

    ax_bl.set_xticks(x_months)
    ax_bl.set_xticklabels(month_lbls, rotation=45, ha="right", fontsize=8)
    ax_bl.set_ylabel("Monthly flow (cfs)")
    ax_bl.set_title("(c) Seasonal profile at representative DV values")

    # Build a combined legend: mode patches + linestyle for DV value
    mode_patches = [
        mpatches.Patch(color=MODE_COLORS[m], label=MODE_LABELS[m])
        for m in ("index", "residual")
    ]
    style_lines = [
        plt.Line2D([0], [0], color="gray", ls=ls, lw=1.2, label=f"$p_k={p}$")
        for p, ls in zip(rep_dvs, rep_styles)
    ]
    ax_bl.legend(
        handles=mode_patches + style_lines,
        frameon=False,
        fontsize=8,
        ncol=2,
    )

    # --- Bottom-right: distinct-value bar chart ---
    modes = ["index", "residual"]
    distinct_counts = [results[m]["distinct_count"] for m in modes]
    colors_bar = [MODE_COLORS[m] for m in modes]
    x_bar = np.arange(len(modes))
    bars = ax_br.bar(x_bar, distinct_counts, color=colors_bar, width=0.5, alpha=0.85)
    ax_br.set_xticks(x_bar)
    ax_br.set_xticklabels([MODE_LABELS[m] for m in modes], fontsize=9)
    ax_br.set_ylabel("Distinct output values (400-point sweep)")
    ax_br.set_title("(d) Step count: distinct flow values")
    ax_br.axhline(n_years_hist, color=COLORS["muted"], ls="--", lw=1,
                  label=f"n_years_hist = {n_years_hist}")
    ax_br.legend(frameon=False, fontsize=8)
    for bar, cnt in zip(bars, distinct_counts):
        ax_br.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.5,
            str(cnt),
            ha="center", va="bottom", fontsize=9,
        )

    fig.suptitle(
        "Kirsch-wrapper DV mapping geometry: index vs residual mode",
        fontsize=12, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[geometry] saved {output_path}")


# ---------------------------------------------------------------------------
# Figure 2: mapping cartoon
# ---------------------------------------------------------------------------

def plot_mapping_cartoon(output_path: Path) -> None:
    """Conceptual schematic of index vs residual DV mapping.

    Two panels side-by-side drawn with matplotlib patches and text.
    Left:  index mode — 12 DVs → one year index → one historical year block.
    Right: residual mode — 12 DVs → 12 quantile lookups → Cholesky transform.
    """
    fig, (ax_idx, ax_res) = plt.subplots(1, 2, figsize=(13, 6))

    _draw_index_cartoon(ax_idx)
    _draw_residual_cartoon(ax_res)

    fig.suptitle(
        "Kirsch-wrapper decision-variable mapping: schematic",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[geometry] saved {output_path}")


def _box(ax, x, y, w, h, text, fc="white", ec="#444", fontsize=9, bold=False):
    """Draw a labeled rectangle."""
    rect = mpatches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02",
        linewidth=1.2,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(rect)
    weight = "bold" if bold else "normal"
    ax.text(
        x + w / 2, y + h / 2, text,
        ha="center", va="center", fontsize=fontsize,
        fontweight=weight, wrap=True,
        color="#1a1a1a",
    )


def _arrow(ax, x0, y0, x1, y1, color="#555"):
    ax.annotate(
        "", xy=(x1, y1), xytext=(x0, y0),
        arrowprops=dict(arrowstyle="->", color=color, lw=1.2),
    )


def _draw_index_cartoon(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Index mode", fontsize=11, fontweight="bold",
                 color=MODE_COLORS["index"], pad=8)

    # --- Row of 12 DV boxes ---
    dv_y = 8.5
    dv_h = 0.7
    dv_w = 0.55
    gap = 0.18
    start_x = 0.4
    for i in range(12):
        xpos = start_x + i * (dv_w + gap)
        fc = "#d6e4f7" if i < 6 else "#eaf3fb"
        _box(ax, xpos, dv_y, dv_w, dv_h, f"$p_{{{i+1}}}$",
             fc=fc, ec=MODE_COLORS["index"], fontsize=7)

    # Brace / arrow down to floor() node
    mid_x = start_x + 5.5 * (dv_w + gap)
    _arrow(ax, mid_x, dv_y - 0.02, mid_x, 7.2)

    # floor() node
    floor_x, floor_y, floor_w, floor_h = 3.8, 6.3, 2.4, 0.75
    _box(ax, floor_x, floor_y, floor_w, floor_h,
         r"$\lfloor p_k \times N_{hist} \rfloor$  → year index $i$",
         fc="#c7dcf5", ec=MODE_COLORS["index"], fontsize=8, bold=True)
    ax.text(
        floor_x + floor_w / 2, floor_y - 0.35,
        "One index shared by all 12 months",
        ha="center", va="top", fontsize=7.5, color=COLORS["muted"], style="italic",
    )

    _arrow(ax, floor_x + floor_w / 2, floor_y - 0.02,
           floor_x + floor_w / 2, 5.2)

    # Historical year block
    yr_x, yr_y, yr_w, yr_h = 1.5, 3.6, 7.0, 1.45
    _box(ax, yr_x, yr_y, yr_w, yr_h,
         "Historical year $i$ block\n"
         "[Oct$_i$  Nov$_i$  Dec$_i$  Jan$_i$  Feb$_i$  Mar$_i$  "
         "Apr$_i$  May$_i$  Jun$_i$  Jul$_i$  Aug$_i$  Sep$_i$]",
         fc="#e8f4e8", ec="#2e7d32", fontsize=8)

    _arrow(ax, yr_x + yr_w / 2, yr_y - 0.02, yr_x + yr_w / 2, 2.4)

    # Output label
    _box(ax, 3.0, 1.4, 4.0, 0.8,
         "Synthetic monthly flows (12 months)",
         fc="#f0f4f8", ec="#555", fontsize=8)

    # Key annotation
    ax.text(
        5.0, 0.6,
        "All 12 months come from the same historical year\n"
        "→ within-year correlation preserved by construction.",
        ha="center", va="center", fontsize=8,
        color="#1a472a",
        bbox=dict(boxstyle="round,pad=0.3", fc="#e8f4e8", ec="#2e7d32", alpha=0.7),
    )

    # DV-space label
    ax.text(
        5.0, 9.45,
        "Decision variables   $\\mathbf{p} \\in [0,1]^{12}$",
        ha="center", va="center", fontsize=9, fontweight="bold",
        color=MODE_COLORS["index"],
    )


def _draw_residual_cartoon(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Residual mode", fontsize=11, fontweight="bold",
                 color=MODE_COLORS["residual"], pad=8)

    # --- Row of 12 DV boxes ---
    dv_y = 8.5
    dv_h = 0.7
    dv_w = 0.55
    gap = 0.18
    start_x = 0.4
    for i in range(12):
        xpos = start_x + i * (dv_w + gap)
        fc = "#fde8d8" if i % 2 == 0 else "#fef3ec"
        _box(ax, xpos, dv_y, dv_w, dv_h, f"$p_{{{i+1}}}$",
             fc=fc, ec=MODE_COLORS["residual"], fontsize=7)

    ax.text(
        5.0, 9.45,
        "Decision variables   $\\mathbf{p} \\in [0,1]^{12}$",
        ha="center", va="center", fontsize=9, fontweight="bold",
        color=MODE_COLORS["residual"],
    )

    # 12 independent arrows pointing down to 12 quantile lookup boxes
    lookup_y = 6.15
    lookup_h = 0.72
    lookup_w = 0.62
    for i in range(12):
        xpos = start_x + i * (dv_w + gap)
        cx = xpos + dv_w / 2
        _arrow(ax, cx, dv_y - 0.02, cx, lookup_y + lookup_h + 0.02)
        _box(ax, xpos - 0.04, lookup_y, lookup_w, lookup_h,
             f"Q$_{{{i+1}}}$\n↕",
             fc="#fde8d8", ec=MODE_COLORS["residual"], fontsize=6.5)

    # Label
    ax.text(
        5.0, lookup_y - 0.35,
        "12 independent quantile lookups\n"
        "(sorted residuals $\\mathbf{Z}_h$ per month $\\times$ site)",
        ha="center", va="top", fontsize=7.5, color=COLORS["muted"], style="italic",
    )

    # Arrow down to Cholesky / normal-score box
    mid_cx = start_x + 5.5 * (dv_w + gap)
    _arrow(ax, mid_cx, lookup_y - 0.02, mid_cx, 4.7)

    chol_x, chol_y, chol_w, chol_h = 2.2, 3.7, 5.6, 0.85
    _box(ax, chol_x, chol_y, chol_w, chol_h,
         "Cholesky / normal-score transform  (SynHydro KirschGenerator)",
         fc="#fde8d8", ec=MODE_COLORS["residual"], fontsize=8, bold=True)

    _arrow(ax, chol_x + chol_w / 2, chol_y - 0.02,
           chol_x + chol_w / 2, 2.8)

    # Output box
    _box(ax, 3.0, 1.85, 4.0, 0.8,
         "Synthetic monthly flows (12 months)",
         fc="#f0f4f8", ec="#555", fontsize=8)

    # Key annotation
    ax.text(
        5.0, 0.6,
        "Each month's quantile is chosen independently;\n"
        "month-to-month dependence relies on the\n"
        "downstream SynHydro Cholesky transform.",
        ha="center", va="center", fontsize=8,
        color="#7b2504",
        bbox=dict(boxstyle="round,pad=0.3", fc="#fde8d8", ec=MODE_COLORS["residual"],
                  alpha=0.7),
    )


# ---------------------------------------------------------------------------
# Summary JSON
# ---------------------------------------------------------------------------

def write_summary(
    results: dict,
    output_dir: Path,
    n_years_hist: int,
) -> None:
    summary = {"n_years_hist": n_years_hist, "modes": {}}
    for mode, res in results.items():
        stats = _flow_change_stats(res["p_grid"], res["flow_at_k"])
        summary["modes"][mode] = {
            "distinct_flow_values": res["distinct_count"],
            "swept_coordinate_k": int(res["k"]),
            "controlled_month_idx": int(res["month_idx"]),
            **stats,
        }
    out = output_dir / "geometry_summary.json"
    with open(out, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[geometry] saved {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Diagnostic: Kirsch-wrapper mapping geometry."
    )
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed (default: 42).")
    p.add_argument("--n-sweep-points", type=int, default=400,
                   help="Number of DV sweep points (default: 400).")
    p.add_argument("--n-years", type=int, default=20,
                   help="Number of output years to generate (default: 20).")
    p.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "outputs" / "diag_kirsch_wrapper_geometry",
        help="Root output directory.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()

    rng = np.random.default_rng(args.seed)

    # ---- Setup output dirs ----
    fig_dir = args.output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # ---- Load data ----
    cache_dir = PROJECT_ROOT / "data" / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    from src.experiment_utils import prepare_data  # noqa: E402
    monthly_2d, monthly_1d = prepare_data(cache_dir)
    n_years_hist = monthly_2d.shape[0]
    print(f"[geometry] n_years_hist = {n_years_hist}")

    # ---- Fit KirschGenerator ----
    from src.kirsch_utils import build_kirsch_generator  # noqa: E402
    gen = build_kirsch_generator(monthly_2d)

    # ---- Build wrappers ----
    from src.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
    wrapper_index = KirschBorgWrapper(gen, mode="index", n_years_out=args.n_years)
    wrapper_residual = KirschBorgWrapper(gen, mode="residual", n_years_out=args.n_years)

    # ---- Choose swept coordinate k ----
    # k is chosen as the midpoint of the residual DV vector so both wrappers
    # sweep a coordinate in the interior of the sequence.
    # residual n_dvs = n_years * 12  → k ≈ n_years * 6
    k_residual = int(SWEEP_K_FRACTION * wrapper_residual.n_dvs)
    # Clamp k for index wrapper (its n_dvs = (n_years+1)*12; use same k value
    # as long as k < index.n_dvs, which it always is since k < n_years*12).
    k_index = min(k_residual, wrapper_index.n_dvs - 1)
    print(f"[geometry] sweeping k={k_residual} (month {k_residual % 12}) "
          f"for residual, k={k_index} (month {k_index % 12}) for index")

    # ---- Run sweeps ----
    print("[geometry] running index sweep ...")
    res_index = run_sweep(wrapper_index, k_index, args.n_sweep_points, rng)
    res_index["_wrapper"] = wrapper_index

    print("[geometry] running residual sweep ...")
    res_residual = run_sweep(wrapper_residual, k_residual, args.n_sweep_points, rng)
    res_residual["_wrapper"] = wrapper_residual

    results = {"index": res_index, "residual": res_residual}

    # ---- Figure 1: smoothness ----
    fig1_path = fig_dir / "fig_mapping_smoothness.pdf"
    plot_mapping_smoothness(results, fig1_path, n_years_hist)

    # ---- Figure 2: cartoon ----
    fig2_path = fig_dir / "fig_mapping_cartoon.pdf"
    plot_mapping_cartoon(fig2_path)

    # ---- Summary JSON ----
    write_summary(results, args.output_dir, n_years_hist)

    print("[geometry] done.")


if __name__ == "__main__":
    main()
