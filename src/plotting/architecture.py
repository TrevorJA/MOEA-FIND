"""Architecture / algorithmic schematics for the MOEA-FIND manuscript.

Reproducible matplotlib-only diagrams (no Inkscape) so that every main-
text figure regenerates from Python source. Current contents:

    fig3_wrapper_schematic → Main §3.4 / Figure 3 (KirschBorgWrapper
                             control flow)

Functions return a matplotlib Figure. They take no data arguments because
architecture schematics do not depend on experimental output.
"""

from __future__ import annotations

from typing import Tuple

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from src.plotting.style import COLORS, apply_style


def _box(
    ax,
    xy: Tuple[float, float],
    width: float,
    height: float,
    text: str,
    fill: str,
    edge: str = "#333333",
    fontsize: int = 8,
    weight: str = "normal",
) -> FancyBboxPatch:
    box = FancyBboxPatch(
        (xy[0] - width / 2, xy[1] - height / 2),
        width, height,
        boxstyle="round,pad=0.02,rounding_size=0.05",
        linewidth=1.0, edgecolor=edge, facecolor=fill, alpha=0.92,
    )
    ax.add_patch(box)
    ax.text(xy[0], xy[1], text, ha="center", va="center",
            fontsize=fontsize, weight=weight)
    return box


def _arrow(ax, start: Tuple[float, float], end: Tuple[float, float],
           label: str = "", color: str = "#333333",
           rad: float = 0.0, fontsize: int = 7) -> None:
    arr = FancyArrowPatch(
        start, end,
        arrowstyle="-|>", mutation_scale=12,
        connectionstyle=f"arc3,rad={rad}",
        color=color, lw=1.1,
    )
    ax.add_patch(arr)
    if label:
        mx = 0.5 * (start[0] + end[0])
        my = 0.5 * (start[1] + end[1]) + 0.03
        ax.text(mx, my, label, ha="center", va="bottom",
                fontsize=fontsize, style="italic", color="#444")


def fig3_wrapper_schematic(
    figsize: Tuple[float, float] = (7.0, 3.6),
) -> plt.Figure:
    """Manuscript Figure 3 — KirschBorgWrapper control flow.

    Pure matplotlib schematic. Boxes for:
        Borg MOEA → wrapper DV decoder → Kirsch-Nowak generator →
        SSI-3 / drought feature extractor → L1 accumulator → Borg archive

    Side callout listing DV injection modes (residual main, index SI)
    and the default epsilon vector used in §§6-7.

    Returns
    -------
    matplotlib.figure.Figure
    """
    apply_style()
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.6)
    ax.set_aspect("equal")
    ax.axis("off")

    # -- Main pipeline boxes (left to right, top row) --
    c_borg = COLORS["empirical"]
    c_wrap = "#ffd8a8"
    c_gen = "#b2df8a"
    c_feat = "#fdbf6f"
    c_obj = COLORS["anti_ideal"]

    _box(ax, (1.0, 4.3), 1.6, 0.7, "Borg MOEA\n(ε-box archive)",
         fill="#cfe2f3", edge=c_borg, weight="bold", fontsize=8)
    _box(ax, (3.2, 4.3), 1.6, 0.7, "KirschBorgWrapper\n(DV decoder)",
         fill=c_wrap, fontsize=8)
    _box(ax, (5.4, 4.3), 1.8, 0.7, "Kirsch-Nowak\nbootstrap generator",
         fill=c_gen, fontsize=8)
    _box(ax, (7.8, 4.3), 1.8, 0.7, "SSI-3 +\ndrought feature extractor",
         fill=c_feat, fontsize=8)

    # -- Bottom row --
    _box(ax, (7.8, 2.5), 1.8, 0.7,
         r"L1 accumulator" "\n" r"$f_{k+1} = \sum_j f_j$",
         fill=c_obj, fontsize=8)
    _box(ax, (3.2, 2.5), 1.6, 0.7,
         r"Objective vector" "\n" r"$(f_1, \ldots, f_{k+1})$",
         fill="#f4cccc", fontsize=8)

    # -- Arrows (counterclockwise loop) --
    _arrow(ax, (1.8, 4.3), (2.4, 4.3), label=r"$\mathbf{x}\in[0,1]^d$")
    _arrow(ax, (4.0, 4.3), (4.5, 4.3), label="residuals / indices")
    _arrow(ax, (6.3, 4.3), (6.9, 4.3), label="monthly trace")
    _arrow(ax, (8.7, 4.3), (8.7, 2.85), label=r"$\mathbf{D}(\mathbf{x})$")
    _arrow(ax, (6.9, 2.5), (4.0, 2.5), label=r"$f_j = |D_j - D^*_j|$")
    _arrow(ax, (2.4, 2.5), (1.0, 2.5), label="return")
    _arrow(ax, (1.0, 2.85), (1.0, 3.95), label="archive update", rad=0.0)

    # -- Side callout: DV modes + epsilon vector --
    callout = FancyBboxPatch(
        (0.15, 0.2), 9.7, 1.5,
        boxstyle="round,pad=0.08,rounding_size=0.1",
        linewidth=0.8, edgecolor="#888888",
        facecolor="#f7f7f7", alpha=0.95,
    )
    ax.add_patch(callout)
    ax.text(0.3, 1.4,
            "Decision-variable injection modes",
            fontsize=9, weight="bold", va="top")
    ax.text(0.3, 1.1,
            "  • Residual (main text): "
            r"$p \mapsto \Phi^{-1}(p)$ standardized residual passed through Cholesky + normal-score pipeline",
            fontsize=7, va="top")
    ax.text(0.3, 0.88,
            "  • Index (SI-4): "
            r"$p \mapsto \lfloor p \cdot N_{\mathrm{years}} \rfloor$ historical-year index injected into the resampler",
            fontsize=7, va="top")
    ax.text(0.3, 0.62,
            "Default ε vector: "
            r"$\varepsilon_{\mathrm{severity}}=0.05,\;$"
            r"$\varepsilon_{\mathrm{duration}}=0.05,\;$"
            r"$\varepsilon_{\mathrm{month}}=0.5,\;$"
            r"$\varepsilon_{L_1}=0.05$",
            fontsize=7, va="top")
    ax.text(0.3, 0.35,
            "Constraints (soft, Deb-style): "
            r"$|\rho_1^{\mathrm{syn}} - \rho_1^{\mathrm{hist}}| < 0.05$, "
            "non-drought annual mean within ±15 % of historical",
            fontsize=7, va="top")

    fig.tight_layout(pad=0.6)
    return fig
