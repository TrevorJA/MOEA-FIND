"""Analytic proof-of-concept figure building blocks.

Each function here is a *reusable panel* targeted at a specific manuscript
figure. Callers (scripts/01, scripts/02, scripts/03, scripts/10) compose
them into final figures with layout conventions pinned to the manuscript
display-item specification (single-column ~3.5 in, 1.5-column ~5.5 in,
double-column ~7.0 in at WRR).

Manuscript cross-reference:
    - fig1_manhattan_concept       → Main §3.1 / Fig 1 (conceptual, no data)
    - fig2_2d_coverage_comparison  → Main §5.1 / Fig 2a
    - fig2_3d_projections          → Main §5.2 / Fig 2b-c
    - fig3_eps_nfe_heatmap         → Main §5.3 / Fig 3
    - fig_si_hyperplane_check      → SI-1 hyperplane-residual histogram
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.analysis import (
    coverage_metrics,
    generate_lhs_samples,
    generate_sobol_samples,
)
from src.plotting.style import COLORS, apply_style


# -----------------------------------------------------------------------------
# Figure 1 — Manhattan-norm concept (schematic, no experimental data)
# -----------------------------------------------------------------------------
def fig1_manhattan_concept(
    anti_ideal: Tuple[float, float] = (3.0, 3.0),
    bounds: Tuple[float, float] = (-3.0, 3.0),
    n_schematic: int = 60,
    figsize: Tuple[float, float] = (7.0, 3.2),
) -> plt.Figure:
    """Manuscript Figure 1 — the Manhattan-norm trick (WRR 1.5-column width).

    Two panels: (a) the drought objective plane with the anti-ideal D* and
    a schematic Pareto front on the Manhattan line; (b) the three-objective
    simplex viewed in 3D showing epsilon-box tiling.
    """
    apply_style()
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1, 1.1], wspace=0.28)

    # -- Panel (a): 2D concept --
    ax_a = fig.add_subplot(gs[0, 0])
    xs = np.linspace(bounds[0], anti_ideal[0], n_schematic)
    ys = anti_ideal[0] + anti_ideal[1] - xs - (anti_ideal[0] - bounds[0])
    ax_a.plot(xs, ys, "-", color=COLORS["empirical"], linewidth=1.6,
              label=r"Pareto front ($J_1 + J_2 + J_3 = C$)")
    ax_a.scatter(xs[::3], ys[::3], s=18, c=COLORS["empirical"], zorder=5,
                 label="Epsilon-box centers")
    ax_a.plot(anti_ideal[0], anti_ideal[1], "x", color=COLORS["anti_ideal"],
              markersize=10, mew=2, label=r"Anti-ideal $D^*$")
    ax_a.plot(bounds[0], bounds[0], ".", color=COLORS["historical"],
              markersize=10, label="Ideal point")
    ax_a.set_xlim(bounds); ax_a.set_ylim(bounds)
    ax_a.set_aspect("equal")
    ax_a.set_xlabel(r"$J_1 = D_1$")
    ax_a.set_ylabel(r"$J_2 = D_2$")
    ax_a.set_title("(a) Manhattan-norm trick in 2D")
    ax_a.legend(fontsize=7, loc="lower left", framealpha=0.9)

    # -- Panel (b): 3D simplex --
    ax_b = fig.add_subplot(gs[0, 1], projection="3d")
    vertices = np.array([
        [anti_ideal[0], bounds[0], bounds[0]],
        [bounds[0], anti_ideal[1], bounds[0]],
        [bounds[0], bounds[0], anti_ideal[0]],
    ])
    tri = np.vstack([vertices, vertices[0]])
    ax_b.plot(tri[:, 0], tri[:, 1], tri[:, 2], color=COLORS["empirical"], lw=1.3)
    rng = np.random.default_rng(0)
    bary = rng.dirichlet(alpha=(1.0, 1.0, 1.0), size=80)
    pts = bary @ vertices
    ax_b.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=8,
                 color=COLORS["empirical"], alpha=0.8)
    ax_b.scatter(*anti_ideal, bounds[0], color=COLORS["anti_ideal"],
                 marker="x", s=40)
    ax_b.set_xlabel(r"$J_1$"); ax_b.set_ylabel(r"$J_2$"); ax_b.set_zlabel(r"$J_3$")
    ax_b.set_title("(b) 2-simplex in 3D")
    ax_b.view_init(elev=22, azim=35)

    return fig


# -----------------------------------------------------------------------------
# Figure 2 — Analytic validation (2D and 3D)
# -----------------------------------------------------------------------------
def fig2_2d_coverage_comparison(
    pareto_dvs: np.ndarray,
    seed: int = 42,
    figsize: Tuple[float, float] = (7.0, 2.4),
) -> plt.Figure:
    """Manuscript Figure 2 panel a — 2D Pareto vs LHS vs Sobol vs Random.

    Four-panel row. Each panel annotates L2-star and NN_CV coverage metrics.
    Sized for WRR 1.5-column width.
    """
    apply_style()
    n = len(pareto_dvs)
    lb = np.array([-3.0, -3.0])
    ub = np.array([3.0, 3.0])

    lhs = generate_lhs_samples(n, 2, lb, ub, seed=seed)
    sobol = generate_sobol_samples(n, 2, lb, ub, seed=seed)[:n]
    rng = np.random.default_rng(seed)
    rand = lb + rng.random((n, 2)) * (ub - lb)

    datasets = [
        ("MOEA-FIND", pareto_dvs, COLORS["empirical"]),
        ("LHS", lhs, COLORS["lhs"]),
        ("Sobol", sobol, COLORS["sobol"]),
        ("Random", rand, COLORS["random"]),
    ]

    fig, axes = plt.subplots(1, 4, figsize=figsize, sharex=True, sharey=True)
    for ax, (label, pts, color) in zip(axes, datasets):
        ax.scatter(pts[:, 0], pts[:, 1], s=1.5, alpha=0.6, c=color, rasterized=True)
        m = coverage_metrics(pts, lb, ub)
        ax.text(0.04, 0.96,
                f"L2*={m['L2_star_discrepancy']:.4f}\nNN_CV={m['nn_cv']:.3f}",
                transform=ax.transAxes, fontsize=7, va="top",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85,
                          edgecolor="none"))
        ax.set_title(label, fontsize=9)
        ax.set_aspect("equal")
        ax.set_xlim(-3.2, 3.2); ax.set_ylim(-3.2, 3.2)

    axes[0].set_ylabel(r"$X_2$")
    for ax in axes:
        ax.set_xlabel(r"$X_1$")
    fig.tight_layout()
    return fig


def fig2_3d_projections(
    pareto_dvs: np.ndarray,
    figsize: Tuple[float, float] = (7.0, 4.8),
) -> plt.Figure:
    """Manuscript Figure 2 panel b — 3D Pareto front projections + 3D scatter.

    Top row: three pairwise projections (X1X2, X1X3, X2X3).
    Bottom row: 3D scatter + coverage-quality bar chart (L2*, NN_CV).
    """
    apply_style()
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.32)

    pairs = [(0, 1, r"$X_1$", r"$X_2$"),
             (0, 2, r"$X_1$", r"$X_3$"),
             (1, 2, r"$X_2$", r"$X_3$")]

    for idx, (i, j, xi, xj) in enumerate(pairs):
        ax = fig.add_subplot(gs[0, idx])
        ax.scatter(pareto_dvs[:, i], pareto_dvs[:, j], s=1.2, alpha=0.4,
                   c=COLORS["empirical"], rasterized=True)
        ax.set_xlim(-3.2, 3.2); ax.set_ylim(-3.2, 3.2)
        ax.set_aspect("equal")
        ax.set_xlabel(xi); ax.set_ylabel(xj)
        ax.set_title(f"{xi} vs {xj}", fontsize=9)

    ax3d = fig.add_subplot(gs[1, 0:2], projection="3d")
    ax3d.scatter(pareto_dvs[:, 0], pareto_dvs[:, 1], pareto_dvs[:, 2],
                 s=1, alpha=0.35, c=COLORS["empirical"], rasterized=True)
    ax3d.set_xlabel(r"$X_1$"); ax3d.set_ylabel(r"$X_2$"); ax3d.set_zlabel(r"$X_3$")
    ax3d.set_title(f"3D Pareto front (n={len(pareto_dvs)})", fontsize=9)
    ax3d.view_init(elev=20, azim=40)

    ax_bar = fig.add_subplot(gs[1, 2])
    lb = np.full(3, -3.0); ub = np.full(3, 3.0)
    n = len(pareto_dvs)
    lhs = generate_lhs_samples(n, 3, lb, ub)
    sobol = generate_sobol_samples(n, 3, lb, ub)[:n]
    methods = ["MOEA-FIND", "LHS", "Sobol"]
    l2 = [coverage_metrics(pareto_dvs, lb, ub)["L2_star_discrepancy"],
          coverage_metrics(lhs, lb, ub)["L2_star_discrepancy"],
          coverage_metrics(sobol, lb, ub)["L2_star_discrepancy"]]
    colors = [COLORS["empirical"], COLORS["lhs"], COLORS["sobol"]]
    bars = ax_bar.bar(methods, l2, color=colors, alpha=0.85)
    ax_bar.set_ylabel("L2* discrepancy")
    ax_bar.set_title("3D coverage quality", fontsize=9)
    for bar, val in zip(bars, l2):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f"{val:.4f}", ha="center", va="bottom", fontsize=7)
    ax_bar.tick_params(axis="x", labelsize=7)

    return fig


# -----------------------------------------------------------------------------
# Figure 3 — Epsilon × NFE sensitivity heatmap
# -----------------------------------------------------------------------------
def fig3_eps_nfe_heatmap(
    aggregated: List[Dict],
    figsize: Tuple[float, float] = (7.0, 2.6),
) -> plt.Figure:
    """Manuscript Figure 3 — coverage quality vs (epsilon, NFE).

    Three panels: (a) archive size, (b) L2-star, (c) NN_CV.
    Input is the aggregated list produced by scripts/03_eps_nfe_sweep.py's
    aggregate() function (one row per (eps, nfe) cell, mean across seeds).
    """
    apply_style()
    eps_vals = sorted({r["epsilon"] for r in aggregated})
    nfe_vals = sorted({r["nfe"] for r in aggregated})
    n_eps, n_nfe = len(eps_vals), len(nfe_vals)

    def grid(field: str) -> np.ndarray:
        g = np.full((n_eps, n_nfe), np.nan)
        for r in aggregated:
            i = eps_vals.index(r["epsilon"])
            j = nfe_vals.index(r["nfe"])
            g[i, j] = r[field]
        return g

    n_sol = grid("n_solutions_mean")
    l2 = grid("L2_star_mean")
    nncv = grid("nn_cv_mean")

    fig, axes = plt.subplots(1, 3, figsize=figsize, sharey=True)
    panels = [
        (n_sol, "(a) Archive size", "viridis", "{:.0f}"),
        (l2, "(b) L2* discrepancy", "magma_r", "{:.3f}"),
        (nncv, "(c) NN spacing CV", "cividis_r", "{:.2f}"),
    ]
    for ax, (g, title, cmap, fmt) in zip(axes, panels):
        im = ax.imshow(g, origin="lower", aspect="auto", cmap=cmap)
        ax.set_xticks(range(n_nfe))
        ax.set_xticklabels([f"{n:,}" for n in nfe_vals], rotation=30, ha="right", fontsize=7)
        ax.set_yticks(range(n_eps))
        ax.set_yticklabels([f"{e:.3f}" for e in eps_vals], fontsize=7)
        ax.set_xlabel("NFE")
        ax.set_title(title, fontsize=9)
        for i in range(n_eps):
            for j in range(n_nfe):
                if not np.isnan(g[i, j]):
                    ax.text(j, i, fmt.format(g[i, j]),
                            ha="center", va="center", fontsize=6,
                            color="white")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    axes[0].set_ylabel(r"Epsilon $\varepsilon$")
    fig.tight_layout()
    return fig


# -----------------------------------------------------------------------------
# SI figure — hyperplane residual histogram
# -----------------------------------------------------------------------------
def fig_si_hyperplane_check(
    objs_2d: Optional[np.ndarray],
    objs_3d: Optional[np.ndarray],
    figsize: Tuple[float, float] = (7.0, 2.6),
) -> plt.Figure:
    """SI-1 hyperplane verification — histograms of Σ J_i − Σ D*_i.

    Shows that Pareto solutions satisfy the Manhattan-norm identity to
    machine precision. Used as Figure SI-1 in the supporting information.
    """
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    def _panel(ax, objs, expected, label, color):
        if objs is None:
            ax.text(0.5, 0.5, "(no data)", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_title(label)
            return
        sums = objs.sum(axis=1)
        dev = sums - expected
        if np.std(dev) < 1e-10:
            ax.bar(["all"], [len(sums)], color=color, alpha=0.85)
            ax.text(0.5, 0.72,
                    f"$\\Sigma J_i$ = {np.mean(sums):.6f}\n"
                    f"(machine precision,\nstd={np.std(dev):.1e})",
                    transform=ax.transAxes, ha="center", fontsize=8)
        else:
            ax.hist(dev, bins=40, color=color, alpha=0.85,
                    edgecolor="black", linewidth=0.4)
            ax.set_xlabel(r"$\Sigma J_i - \Sigma D^*_i$")
            ax.set_ylabel("count")
        ax.set_title(label, fontsize=9)

    _panel(axes[0], objs_2d, expected=6.0,
           label="(a) 2D: expected $\\Sigma$=6", color=COLORS["empirical"])
    _panel(axes[1], objs_3d, expected=9.0,
           label="(b) 3D: expected $\\Sigma$=9", color=COLORS["parametric"])

    fig.tight_layout()
    return fig
