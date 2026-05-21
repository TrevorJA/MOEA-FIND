"""Analytic proof-of-concept figure building blocks.

Each function here is a *reusable panel* targeted at a specific manuscript
figure. Callers (workflows/01_analytic_validation/ and 99_manuscript_figures/) compose
them into final figures with layout conventions pinned to the manuscript
display-item specification (single-column ~3.5 in, 1.5-column ~5.5 in,
double-column ~7.0 in at WRR).

Manuscript cross-reference (restructured 2026-04-14 for the 7-figure
narrative sequence):
    - fig1_param_vs_hazard_space   → Main §1.3 / Figure 1
    - fig3_manhattan_construction  → Main §2.4 / Figure 3
    - fig4_dimension_sweep         → Main §3.1 / Figure 4
    - fig_si_hyperplane_check      → SI-1 residual histogram
    - fig3_eps_nfe_heatmap         → SI-3 epsilon x NFE heatmap

Legacy (deprecated, used by archived scripts):
    - fig1_manhattan_concept, fig2_2d_coverage_comparison,
      fig2_3d_projections, fig2_metrics_bar
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.discovery.analysis import (
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
# Figure 2 panel c — coverage metrics bar chart
# -----------------------------------------------------------------------------
def fig2_metrics_bar(
    pareto_points: np.ndarray,
    lhs_points: np.ndarray,
    sobol_points: np.ndarray,
    random_points: Optional[np.ndarray] = None,
    lb: Optional[np.ndarray] = None,
    ub: Optional[np.ndarray] = None,
    figsize: Tuple[float, float] = (7.0, 2.4),
) -> plt.Figure:
    """Manuscript Figure 2 panel c — NN_CV, L2*, centered discrepancy.

    Three-panel bar chart comparing coverage metrics across MOEA-FIND,
    LHS, Sobol, and (optionally) uniform-random samples of the same size
    on the same bounding box. Sized as a single row so it can be stacked
    below the 2D/3D scatter panels (panels a and b) of Figure 2.

    Parameters
    ----------
    pareto_points, lhs_points, sobol_points : np.ndarray
        (n, d) arrays of equal-size samples. All three must share the
        same dimensionality d.
    random_points : np.ndarray, optional
        Optional fourth sampler to include as a "Random" baseline.
    lb, ub : np.ndarray, optional
        Bounding-box lower and upper bounds. If None, inferred from the
        min/max of the union of all samples.
    """
    apply_style()
    samples = [("MOEA-FIND", pareto_points, COLORS["empirical"]),
               ("LHS",       lhs_points,    COLORS["lhs"]),
               ("Sobol",     sobol_points,  COLORS["sobol"])]
    if random_points is not None:
        samples.append(("Random", random_points, COLORS["random"]))

    if lb is None or ub is None:
        stack = np.vstack([s[1] for s in samples])
        lb = stack.min(axis=0)
        ub = stack.max(axis=0)

    rows = [coverage_metrics(s[1], lb, ub) for s in samples]
    labels = [s[0] for s in samples]
    colors = [s[2] for s in samples]

    l2 = [r["L2_star_discrepancy"] for r in rows]
    nn_cv = [r["nn_cv"] for r in rows]
    nn_mean = [r["nn_mean"] for r in rows]

    fig, axes = plt.subplots(1, 3, figsize=figsize, sharey=False)
    panels = [
        (axes[0], l2, "L2* discrepancy", "{:.4f}"),
        (axes[1], nn_cv, "NN spacing CV", "{:.3f}"),
        (axes[2], nn_mean, "NN mean distance", "{:.3f}"),
    ]
    x = np.arange(len(labels))
    for ax, vals, title, fmt in panels:
        bars = ax.bar(x, vals, color=colors, alpha=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=7)
        ax.set_title(title, fontsize=9)
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    fmt.format(v), ha="center", va="bottom", fontsize=7)
        ax.margins(y=0.15)

    fig.tight_layout()
    return fig


# -----------------------------------------------------------------------------
# Figure 3 (SI-2 in new numbering) — Epsilon × NFE sensitivity heatmap
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


# =============================================================================
# Main-text Figure 3 — Manhattan-distance auxiliary objective construction
# =============================================================================
def fig3_manhattan_construction(
    diag_root=None,
    figsize: Tuple[float, float] = (11.0, 3.6),
) -> plt.Figure:
    """Manuscript Figure 3 — space-filling behavior in the K-dim objective space.

    Uses the real Borg archive from the K=2 shell-vs-interior diagnostic
    to make three points:
      (a) the archive fills the K-dim drought-objective cube — it is not
          confined to a simplex or any other codimension-1 manifold;
      (b) lifting each archive point into (K+1) dimensions by the
          auxiliary J_{K+1} = ||D - D^*||_1 shows that the auxiliary is a
          deterministic function of (D_1, ..., D_K) — a piecewise-linear
          tent over the cube, not an independent axis;
      (c) side-by-side with an equal-size LHS sample inside the same
          feasible cube, the MOEA-FIND archive matches the space-filling
          density quality of a designed space-filling sampler.

    Parameters
    ----------
    diag_root : path-like, optional
        Directory containing k2/samples.npz and k2/results.json from
        ``workflows/01_analytic_validation/dimension_sweep.py``. Defaults to
        ``outputs/01_analytic_validation/dimension_sweep`` relative to the project root.
    """
    from pathlib import Path

    apply_style()

    # --- Load K=2 Borg + LHS samples from the dimension-sweep compute ---
    if diag_root is None:
        here = Path(__file__).resolve()
        diag_root = (
            here.parents[2] / "outputs" / "01_analytic_validation" / "dimension_sweep"
        )
    diag_root = Path(diag_root)
    k2_dir = diag_root / "k2"
    samples_path = k2_dir / "samples.npz"
    if not samples_path.exists():
        raise FileNotFoundError(
            f"Fig 3 requires Borg samples at {samples_path}. "
            "Run workflows/01_analytic_validation/dimension_sweep.py --k 2."
        )
    samples = dict(np.load(samples_path))
    borg = samples["borg"]     # (n, 2) — J_1 = x_1, J_2 = x_2
    lhs = samples["lhs"]
    D_star = np.array([3.0, 3.0])
    half_width = 2.5  # feasible cube half-width from the diagnostic

    # Auxiliary objective: J_{K+1} = Manhattan distance to anti-ideal D*
    def aux(pts):
        return np.abs(pts - D_star).sum(axis=1)

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 1.2, 1.0], wspace=0.32)

    # -- Panel (a): Borg archive fills the K-dim drought objective cube --
    ax_a = fig.add_subplot(gs[0, 0])
    feas = plt.matplotlib.patches.Rectangle(
        (-half_width, -half_width), 2 * half_width, 2 * half_width,
        fill=False, lw=0.8, ls="--", edgecolor="black",
        label=r"feasible cube",
    )
    ax_a.add_patch(feas)
    ax_a.scatter(borg[:, 0], borg[:, 1], s=10,
                 color=COLORS["parametric"], alpha=0.6,
                 edgecolor="none", label=f"MOEA-FIND archive (n={len(borg)})")
    ax_a.plot(*D_star, "x", color=COLORS["anti_ideal"], mew=2.2, ms=12,
              label=r"anti-ideal $D^\star$")
    ax_a.set_xlim(-3.3, 3.5); ax_a.set_ylim(-3.3, 3.5)
    ax_a.set_aspect("equal")
    ax_a.set_xlabel(r"$J_1 = D_1$")
    ax_a.set_ylabel(r"$J_2 = D_2$")
    ax_a.set_title(r"(a) archive fills $K$-dim objective cube",
                   fontsize=10)
    ax_a.legend(fontsize=7, loc="lower left", framealpha=0.9)

    # -- Panel (b): K+1 lift — auxiliary is a deterministic tent surface --
    ax_b = fig.add_subplot(gs[0, 1], projection="3d")
    # Render the piecewise-linear tent J_{K+1}(J_1, J_2) = |J_1-3| + |J_2-3|
    gg = np.linspace(-half_width, half_width, 30)
    J1, J2 = np.meshgrid(gg, gg)
    Jaux = np.abs(J1 - D_star[0]) + np.abs(J2 - D_star[1])
    ax_b.plot_surface(J1, J2, Jaux, cmap="Greys", alpha=0.25,
                      edgecolor="none", rstride=1, cstride=1)
    # Scatter Borg archive lifted to (J1, J2, J_{K+1})
    ax_b.scatter(borg[:, 0], borg[:, 1], aux(borg), s=6,
                 color=COLORS["parametric"], alpha=0.65,
                 label="MOEA-FIND archive")
    # Anti-ideal marker at J_{K+1}=0
    ax_b.scatter([D_star[0]], [D_star[1]], [0],
                 color=COLORS["anti_ideal"], marker="x", s=60, linewidth=2.2,
                 label=r"$D^\star$")
    ax_b.set_xlabel(r"$J_1$", labelpad=-2)
    ax_b.set_ylabel(r"$J_2$", labelpad=-2)
    ax_b.set_zlabel(r"$J_{K+1}=\|D-D^\star\|_1$", labelpad=-2)
    ax_b.set_title(r"(b) auxiliary $J_{K+1}$ is a function of $J_{1:K}$",
                   fontsize=10)
    ax_b.view_init(elev=22, azim=-135)
    for lbl in (ax_b.get_xticklabels() + ax_b.get_yticklabels()
                + ax_b.get_zticklabels()):
        lbl.set_fontsize(7)
    ax_b.legend(fontsize=7, loc="upper right", framealpha=0.9)

    # -- Panel (c): MOEA-FIND vs equal-size LHS in the feasible cube --
    ax_c = fig.add_subplot(gs[0, 2])
    feas2 = plt.matplotlib.patches.Rectangle(
        (-half_width, -half_width), 2 * half_width, 2 * half_width,
        fill=False, lw=0.8, ls="--", edgecolor="black",
    )
    ax_c.add_patch(feas2)
    # Match n for a fair overlay
    n_show = min(len(borg), len(lhs))
    rng = np.random.default_rng(0)
    idx_b = rng.choice(len(borg), size=n_show, replace=False)
    idx_l = rng.choice(len(lhs), size=n_show, replace=False)
    ax_c.scatter(lhs[idx_l, 0], lhs[idx_l, 1], s=10,
                 color=COLORS["muted"], alpha=0.5,
                 edgecolor="none", label="LHS reference")
    ax_c.scatter(borg[idx_b, 0], borg[idx_b, 1], s=10,
                 color=COLORS["parametric"], alpha=0.55,
                 edgecolor="none", label="MOEA-FIND")
    ax_c.plot(*D_star, "x", color=COLORS["anti_ideal"], mew=2.2, ms=12)
    ax_c.set_xlim(-3.3, 3.5); ax_c.set_ylim(-3.3, 3.5)
    ax_c.set_aspect("equal")
    ax_c.set_xlabel(r"$J_1 = D_1$")
    ax_c.set_ylabel(r"$J_2 = D_2$")
    ax_c.set_title("(c) matches space-filling reference", fontsize=10)
    ax_c.legend(fontsize=7, loc="lower left", framealpha=0.9)

    return fig


# =============================================================================
# Main-text Figure 4 — Dimension sweep at K = 2 through K = 6
# =============================================================================
def _load_dimension_sweep(diag_root) -> Dict[int, Dict[str, object]]:
    """Load all per-K diagnostic outputs from the dimension sweep.

    Slugs follow ``analytic_slug(k=K, nfe=N, seed=S)`` ->
    ``analytic__k=K__nfe=...__s=...``. If multiple slugs exist for the
    same K (e.g. different NFE budgets), the one with the highest NFE
    wins (parsed back via :func:`src.io_paths.slugs.parse_slug`).
    """
    import json
    from pathlib import Path
    from src.io_paths.slugs import parse_slug
    diag_root = Path(diag_root)
    out: Dict[int, Dict[str, object]] = {}
    if not diag_root.exists():
        return out
    by_k: Dict[int, list] = {}
    for sub in diag_root.iterdir():
        if not sub.is_dir():
            continue
        try:
            parsed = parse_slug(sub.name)
        except Exception:  # noqa: BLE001
            continue
        if parsed.get("_stage") != "analytic" or "k" not in parsed:
            continue
        try:
            k = int(parsed["k"])
        except ValueError:
            continue
        if not ((sub / "results.json").exists() and (sub / "samples.npz").exists()):
            continue
        nfe_token = parsed.get("nfe", 0)
        nfe = nfe_token if isinstance(nfe_token, int) else 0
        by_k.setdefault(k, []).append((nfe, sub))
    for k, candidates in by_k.items():
        nfe, kd = max(candidates, key=lambda t: t[0])
        out[k] = {
            "results": json.loads((kd / "results.json").read_text()),
            "samples": dict(np.load(kd / "samples.npz")),
        }
    return out


def fig4_dimension_sweep(
    diag_root,
    figsize: Tuple[float, float] = (12.0, 7.4),
) -> plt.Figure:
    """Manuscript Figure 4 — interior-filling coverage across K = 2..6.

    Five panels (top row: 3D scatter; bottom row: per-K metric lines):
      (a) K=3 MOEA-FIND archive rendered as a 3D scatter inside the
          feasible cube, anti-ideal marked at the corner.
      (b) K=3 LHS and Sobol reference samples overlaid in 3D, same cube.
      (c) Mean Manhattan distance from the anti-ideal versus K.
      (d) Interior mass fraction versus K.
      (e) Signed orthant occupancy fraction versus K.

    The four samplers (MOEA-FIND, uniform, LHS, Sobol) use distinct
    line styles in panels (c)-(e) so that overlapping curves can still
    be traced individually.

    Parameters
    ----------
    diag_root : path-like
        Directory containing k{2..6}/results.json and samples.npz, e.g.
        ``outputs/01_analytic_validation/dimension_sweep``.
    """
    apply_style()
    data = _load_dimension_sweep(diag_root)
    if not data:
        raise FileNotFoundError(
            f"No k{{2..6}} outputs found under {diag_root}"
        )

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 6, figure=fig, wspace=1.1, hspace=0.35,
                  height_ratios=[1.1, 1.0])
    ax_a = fig.add_subplot(gs[0, 0:3], projection="3d")
    ax_b = fig.add_subplot(gs[0, 3:6], projection="3d")
    ax_c = fig.add_subplot(gs[1, 0:2])
    ax_d = fig.add_subplot(gs[1, 2:4])
    ax_e = fig.add_subplot(gs[1, 4:6])

    sampler_colors = {
        "MOEA-FIND": COLORS["parametric"],
        "uniform_in_cube": COLORS["muted"],
        "lhs_in_cube": COLORS["empirical"],
        "sobol_in_cube": COLORS["sobol"],
    }
    sampler_linestyles = {
        "MOEA-FIND": "-",
        "uniform_in_cube": "--",
        "lhs_in_cube": "-.",
        "sobol_in_cube": ":",
    }
    sampler_markers = {
        "MOEA-FIND": "o",
        "uniform_in_cube": "s",
        "lhs_in_cube": "^",
        "sobol_in_cube": "D",
    }
    sampler_labels = {
        "MOEA-FIND": "MOEA-FIND",
        "uniform_in_cube": "uniform",
        "lhs_in_cube": "LHS",
        "sobol_in_cube": "Sobol",
    }
    sampler_keys = ["MOEA-FIND", "uniform_in_cube", "lhs_in_cube",
                    "sobol_in_cube"]

    def _cube_wireframe(ax, r: float) -> None:
        """Draw a dashed wireframe of the feasible cube [-r, r]^3."""
        corners = np.array([
            [sx, sy, sz]
            for sx in (-r, r) for sy in (-r, r) for sz in (-r, r)
        ])
        edges = [
            (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
            (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7),
        ]
        for i, j in edges:
            ax.plot(
                [corners[i, 0], corners[j, 0]],
                [corners[i, 1], corners[j, 1]],
                [corners[i, 2], corners[j, 2]],
                color="black", lw=0.6, ls="--", alpha=0.55,
            )

    # --- Panel (a): K=3 MOEA-FIND archive in 3D ---
    if 3 in data:
        k3 = data[3]
        samples = k3["samples"]
        results = k3["results"]
        radius = results.get("feasible_radius", 2.5)
        _cube_wireframe(ax_a, radius)
        borg = samples["borg"]
        # Downsample for visual clarity
        n_show = min(800, len(borg))
        rng = np.random.default_rng(0)
        idx_b = rng.choice(len(borg), size=n_show, replace=False)
        ax_a.scatter(borg[idx_b, 0], borg[idx_b, 1], borg[idx_b, 2],
                     s=6, color=COLORS["parametric"], alpha=0.6,
                     label="MOEA-FIND")
        ax_a.scatter([3], [3], [3], marker="X",
                     color=COLORS["anti_ideal"], s=60,
                     label=r"anti-ideal $D^\star$")
        ax_a.set_xlim(-3.2, 3.2); ax_a.set_ylim(-3.2, 3.2); ax_a.set_zlim(-3.2, 3.2)
        ax_a.set_xlabel(r"$x_1$", labelpad=-2)
        ax_a.set_ylabel(r"$x_2$", labelpad=-2)
        ax_a.set_zlabel(r"$x_3$", labelpad=-2)
        ax_a.view_init(elev=22, azim=42)
        for lbl in (ax_a.get_xticklabels() + ax_a.get_yticklabels()
                    + ax_a.get_zticklabels()):
            lbl.set_fontsize(7)
        ax_a.legend(fontsize=7, loc="upper left", framealpha=0.9)
    ax_a.set_title(r"(a) $K=3$ MOEA-FIND archive", fontsize=10)

    # --- Panel (b): K=3 LHS + Sobol reference samples in 3D ---
    if 3 in data:
        samples = data[3]["samples"]
        radius = data[3]["results"].get("feasible_radius", 2.5)
        _cube_wireframe(ax_b, radius)
        lhs = samples["lhs"]
        sobol = samples["sobol"]
        n_show = min(800, len(lhs), len(sobol))
        rng = np.random.default_rng(1)
        idx_l = rng.choice(len(lhs), size=n_show, replace=False)
        idx_s = rng.choice(len(sobol), size=n_show, replace=False)
        ax_b.scatter(lhs[idx_l, 0], lhs[idx_l, 1], lhs[idx_l, 2],
                     s=6, color=COLORS["empirical"], alpha=0.45,
                     label="LHS")
        ax_b.scatter(sobol[idx_s, 0], sobol[idx_s, 1], sobol[idx_s, 2],
                     s=6, color=COLORS["sobol"], alpha=0.45, marker="^",
                     label="Sobol")
        ax_b.scatter([3], [3], [3], marker="X",
                     color=COLORS["anti_ideal"], s=60,
                     label=r"anti-ideal $D^\star$")
        ax_b.set_xlim(-3.2, 3.2); ax_b.set_ylim(-3.2, 3.2); ax_b.set_zlim(-3.2, 3.2)
        ax_b.set_xlabel(r"$x_1$", labelpad=-2)
        ax_b.set_ylabel(r"$x_2$", labelpad=-2)
        ax_b.set_zlabel(r"$x_3$", labelpad=-2)
        ax_b.view_init(elev=22, azim=42)
        for lbl in (ax_b.get_xticklabels() + ax_b.get_yticklabels()
                    + ax_b.get_zticklabels()):
            lbl.set_fontsize(7)
        ax_b.legend(fontsize=7, loc="upper left", framealpha=0.9)
    ax_b.set_title(r"(b) $K=3$ LHS + Sobol references", fontsize=10)

    # Helper to pull per-K metric series
    def _series(metric_key: str, subkey: str = None):
        Ks = sorted(data.keys())
        series: Dict[str, List[float]] = {k: [] for k in sampler_keys}
        for K in Ks:
            r = data[K]["results"]
            for s in sampler_keys:
                if s in r:
                    val = r[s][metric_key]
                    if subkey is not None:
                        val = val[subkey]
                    series[s].append(val)
                else:
                    series[s].append(np.nan)
        return Ks, series

    def _plot_series(ax, Ks, series, ylabel, title):
        for s in sampler_keys:
            ax.plot(Ks, series[s],
                    linestyle=sampler_linestyles[s],
                    marker=sampler_markers[s],
                    color=sampler_colors[s],
                    label=sampler_labels[s], markersize=5, lw=1.4)
        ax.set_xticks(Ks)
        ax.set_xlabel(r"target dimensionality $K$")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, loc="best", framealpha=0.9)

    # --- Panel (c): mean Manhattan distance from anti-ideal vs K ---
    Ks, ser_mean = _series("dist_from_D*", subkey="mean")
    _plot_series(ax_c, Ks, ser_mean,
                 r"mean $L^1$ distance from $D^\star$",
                 r"(c) Manhattan distance to anti-ideal")

    # --- Panel (d): interior mass fraction vs K ---
    Ks, ser_int = _series("interior_fraction")
    _plot_series(ax_d, Ks, ser_int,
                 "interior mass fraction",
                 "(d) interior mass fraction")
    ax_d.set_ylim(0.0, 1.0)

    # --- Panel (e): orthant occupancy fraction vs K ---
    Ks, ser_orth = _series("orthant_occupancy", subkey="fraction")
    _plot_series(ax_e, Ks, ser_orth,
                 "signed orthant occupancy",
                 r"(e) signed orthant occupancy ($2^K$ orthants)")
    ax_e.set_ylim(0.0, 1.05)

    return fig


# =============================================================================
# Main-text Figure 1 — Parameter space versus drought hazard space
# =============================================================================
def fig1_param_vs_hazard_space(
    param_pts: np.ndarray,
    hazard_pts: np.ndarray,
    historical: Optional[np.ndarray] = None,
    param_labels: Tuple[str, str] = (r"$x_1$", r"$x_2$"),
    hazard_labels: Tuple[str, str] = (
        r"mean drought duration $D_2$ (months)",
        r"mean drought severity $D_1$ (SSI units)",
    ),
    figsize: Tuple[float, float] = (10.5, 3.4),
) -> plt.Figure:
    """Manuscript Figure 1 — parameter space vs drought hazard space contrast.

    Three panels: (a) space-filling LHS design in two generator decision
    variables, colored by position in the design; (b) schematic arrow
    labeled g denoting the generator map; (c) the same points projected
    into the SSI-3 drought hazard space under g, coloured identically
    to panel a, with the historical reference drought marked.

    Parameters
    ----------
    param_pts : (n, 2) array
        Latin hypercube design in two generator decision variables.
    hazard_pts : (n, 2) array
        Same n points projected into (duration, severity) drought
        hazard space under the generator map. Rows must correspond to
        ``param_pts`` rows.
    historical : (2,) array, optional
        Historical reference drought in (duration, severity) for panel c.
    """
    apply_style()
    if param_pts.shape[0] != hazard_pts.shape[0]:
        raise ValueError("param_pts and hazard_pts must have same length")

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 0.45, 1.0], wspace=0.28)

    # Colour by position in parameter space (left-to-right gradient)
    order = np.argsort(param_pts[:, 0] + 0.5 * param_pts[:, 1])
    rank = np.empty_like(order)
    rank[order] = np.arange(len(order))
    cmap = plt.get_cmap("viridis")
    colors = cmap(rank / max(1, len(order) - 1))

    # -- Panel (a): parameter space --
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.scatter(param_pts[:, 0], param_pts[:, 1], s=14, c=colors,
                 edgecolor="black", linewidth=0.2)
    ax_a.set_xlabel(param_labels[0])
    ax_a.set_ylabel(param_labels[1])
    ax_a.set_title("(a) parameter space: LHS design", fontsize=10)
    ax_a.set_aspect("equal", adjustable="box")

    # -- Panel (b): generator map arrow --
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.set_xlim(0, 1); ax_b.set_ylim(0, 1)
    ax_b.axis("off")
    ax_b.annotate(
        "", xy=(0.95, 0.5), xytext=(0.05, 0.5),
        arrowprops=dict(arrowstyle="->", lw=2.0, color=COLORS["historical"]),
    )
    ax_b.text(0.5, 0.62, r"generator map $g$",
              ha="center", va="bottom", fontsize=11)
    ax_b.text(
        0.5, 0.30,
        "Kirsch-Nowak\n+ SSI-3 extraction",
        ha="center", va="top", fontsize=8, color=COLORS["muted"],
    )

    # -- Panel (c): hazard space projection --
    ax_c = fig.add_subplot(gs[0, 2])
    ax_c.scatter(hazard_pts[:, 0], hazard_pts[:, 1], s=14, c=colors,
                 edgecolor="black", linewidth=0.2)
    if historical is not None:
        ax_c.plot(historical[0], historical[1], "x",
                  color=COLORS["anti_ideal"], mew=2, ms=12,
                  label="historical drought")
        ax_c.legend(fontsize=8, loc="best", framealpha=0.9)
    ax_c.set_xlabel(hazard_labels[0])
    ax_c.set_ylabel(hazard_labels[1])
    ax_c.set_title("(c) drought hazard space: projection under $g$",
                   fontsize=10)

    return fig
