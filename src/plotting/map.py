"""Study-area map plots for the MOEA-FIND manuscript.

Produces Figure 4 — the Delaware River Basin study area — in pure
matplotlib. The basin polygon, reservoirs, gauges, and flow targets
are hard-coded as lat/lon coordinates so that the figure renders on
any machine without a working shapely/geopandas installation.

The coordinates are approximate and intended for a schematic main-text
display item, not for quantitative spatial analysis. When a vetted
shapefile from the Pywr-DRB data package becomes available, this module
can be upgraded to a geopandas-based implementation without changing
the function signatures or the script 10 wiring.

Manuscript cross-reference:
    fig4_drb_map → Main §4.1 / Figure 4 (Delaware River Basin study area)
"""

from __future__ import annotations

from typing import Dict, Iterable, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from src.plotting.style import COLORS, apply_style


# -----------------------------------------------------------------------------
# Hard-coded DRB geometry (approximate; schematic only)
# -----------------------------------------------------------------------------
# Outer basin boundary — anti-clockwise polygon in (lon, lat). Digitized by
# hand from the Pywr-DRB v2 basin map; values are rounded to 0.05 degrees.
DRB_BASIN_OUTLINE: Tuple[Tuple[float, float], ...] = (
    (-75.05, 42.55), (-74.70, 42.45), (-74.55, 42.20), (-74.40, 41.95),
    (-74.45, 41.65), (-74.65, 41.40), (-74.85, 41.15), (-75.00, 40.95),
    (-75.20, 40.75), (-75.15, 40.45), (-75.05, 40.15), (-75.20, 39.90),
    (-75.40, 39.65), (-75.55, 39.50), (-75.55, 39.30), (-75.40, 39.20),
    (-75.25, 39.10), (-75.10, 39.05), (-74.95, 39.10), (-74.80, 39.20),
    (-74.70, 39.40), (-74.75, 39.65), (-74.85, 39.90), (-75.00, 40.15),
    (-75.05, 40.35), (-75.00, 40.55), (-74.90, 40.70), (-74.75, 40.85),
    (-74.65, 41.00), (-74.60, 41.20), (-74.70, 41.45), (-74.80, 41.70),
    (-74.95, 41.95), (-75.10, 42.25), (-75.05, 42.55),
)

# NYC supply reservoirs (lon, lat, label).
DRB_NYC_RESERVOIRS: Tuple[Tuple[float, float, str], ...] = (
    (-75.00, 42.05, "Cannonsville"),
    (-74.90, 42.10, "Pepacton"),
    (-74.60, 41.85, "Neversink"),
)

# USGS inflow gauges actually used by MOEA-FIND (lon, lat, gauge_id).
DRB_GAUGES: Tuple[Tuple[float, float, str], ...] = (
    (-75.00, 42.05, "USGS 01425000"),  # Cannonsville
    (-74.90, 42.10, "USGS 01417000"),  # Pepacton
    (-74.60, 41.85, "USGS 01435000"),  # Neversink
)

# FFMP policy-relevant flow targets.
DRB_FLOW_TARGETS: Tuple[Tuple[float, float, str], ...] = (
    (-74.80, 41.30, "Montague"),
    (-74.75, 40.20, "Trenton"),
)

# Approximate mainstem reach — used as a stylized river line.
DRB_MAINSTEM: Tuple[Tuple[float, float], ...] = (
    (-75.05, 42.55), (-74.95, 42.15), (-74.80, 41.85), (-74.75, 41.50),
    (-74.85, 41.15), (-74.90, 40.80), (-74.95, 40.45), (-75.10, 40.10),
    (-75.20, 39.85), (-75.35, 39.60), (-75.55, 39.40),
)


def _draw_basin(ax) -> None:
    poly = Polygon(
        DRB_BASIN_OUTLINE,
        closed=True,
        facecolor="#e8eef5",
        edgecolor=COLORS["historical"],
        linewidth=1.2,
        alpha=0.9,
    )
    ax.add_patch(poly)


def _draw_mainstem(ax) -> None:
    lon = [p[0] for p in DRB_MAINSTEM]
    lat = [p[1] for p in DRB_MAINSTEM]
    ax.plot(lon, lat, color=COLORS["empirical"], linewidth=1.6, zorder=3,
            label="Delaware mainstem (schematic)")


def _draw_reservoirs(ax) -> None:
    for lon, lat, label in DRB_NYC_RESERVOIRS:
        ax.plot(lon, lat, marker="s", markersize=8,
                markerfacecolor=COLORS["anti_ideal"],
                markeredgecolor="black", markeredgewidth=0.8,
                zorder=5)
        ax.annotate(label, (lon, lat), textcoords="offset points",
                    xytext=(6, 3), fontsize=8)


def _draw_gauges(ax) -> None:
    for lon, lat, gid in DRB_GAUGES:
        ax.plot(lon, lat, marker="^", markersize=7,
                markerfacecolor=COLORS["lhs"],
                markeredgecolor="black", markeredgewidth=0.6,
                zorder=4)


def _draw_flow_targets(ax) -> None:
    for lon, lat, label in DRB_FLOW_TARGETS:
        ax.plot(lon, lat, marker="X", markersize=10,
                markerfacecolor=COLORS["highlight"],
                markeredgecolor="black", markeredgewidth=0.8,
                zorder=5)
        ax.annotate(label, (lon, lat), textcoords="offset points",
                    xytext=(6, -3), fontsize=8)


def fig4_drb_map(
    figsize: Tuple[float, float] = (5.0, 7.0),
    show_inset: bool = True,
) -> plt.Figure:
    """Manuscript Figure 4 — Delaware River Basin study area.

    Pure-matplotlib schematic. A (hand-digitized) DRB boundary polygon,
    a stylized Delaware mainstem river, the three NYC supply reservoirs
    (Cannonsville, Pepacton, Neversink), the USGS inflow gauges used by
    MOEA-FIND, and the Montague/Trenton FFMP flow targets. An optional
    small regional inset gives geographic context.

    Parameters
    ----------
    figsize : tuple of float
        Figure size in inches. Default is a single-column WRR figure.
    show_inset : bool
        If True, draw a small regional inset (east-coast US) in the
        upper-right showing the DRB's location.

    Returns
    -------
    matplotlib.figure.Figure
    """
    apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    _draw_basin(ax)
    _draw_mainstem(ax)
    _draw_gauges(ax)
    _draw_reservoirs(ax)
    _draw_flow_targets(ax)

    ax.set_xlim(-75.8, -74.1)
    ax.set_ylim(39.0, 42.7)
    ax.set_aspect(1.3)  # rough Mercator correction at ~41 N
    ax.set_xlabel("Longitude (°W)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_title("Delaware River Basin study area")

    # Manual legend (mixed marker types)
    legend_handles = [
        plt.Line2D([0], [0], marker="s", linestyle="",
                   markerfacecolor=COLORS["anti_ideal"],
                   markeredgecolor="black", markersize=8,
                   label="NYC supply reservoir"),
        plt.Line2D([0], [0], marker="^", linestyle="",
                   markerfacecolor=COLORS["lhs"],
                   markeredgecolor="black", markersize=7,
                   label="USGS inflow gauge"),
        plt.Line2D([0], [0], marker="X", linestyle="",
                   markerfacecolor=COLORS["highlight"],
                   markeredgecolor="black", markersize=9,
                   label="FFMP flow target"),
        plt.Line2D([0], [0], color=COLORS["empirical"], linewidth=1.6,
                   label="Delaware mainstem"),
    ]
    ax.legend(handles=legend_handles, loc="lower left", fontsize=7,
              framealpha=0.95)

    if show_inset:
        # Small regional inset: east-coast US rectangle with DRB bbox.
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        iax = inset_axes(ax, width="32%", height="22%", loc="upper right",
                         borderpad=0.5)
        # Rough US east-coast rectangle (for orientation only).
        iax.add_patch(Polygon(
            [(-82, 35), (-69, 35), (-69, 45), (-82, 45)],
            facecolor="#f0f0f0", edgecolor="#888", linewidth=0.6,
        ))
        # DRB bounding box highlight
        iax.add_patch(Polygon(
            [(-75.7, 39.1), (-74.3, 39.1), (-74.3, 42.6), (-75.7, 42.6)],
            facecolor=COLORS["anti_ideal"], alpha=0.3,
            edgecolor=COLORS["anti_ideal"], linewidth=1.0,
        ))
        iax.set_xlim(-83, -68)
        iax.set_ylim(34, 46)
        iax.set_xticks([]); iax.set_yticks([])
        iax.set_aspect(1.3)
        for spine in iax.spines.values():
            spine.set_linewidth(0.6); spine.set_edgecolor("#888")

    fig.tight_layout()
    return fig
