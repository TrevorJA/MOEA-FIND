"""Shared matplotlib style configuration for MOEA-FIND figures.

Provides consistent styling across all diagnostic and publication figures.
Import and call `apply_style()` at the top of any plotting script.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

# Color palette
COLORS = {
    "empirical": "#1f77b4",   # blue
    "kde": "#2ca02c",         # green
    "parametric": "#d62728",  # red
    "lhs": "#ff7f0e",         # orange
    "sobol": "#9467bd",       # purple
    "random": "#8c564b",      # brown
    "historical": "#000000",  # black
    "anti_ideal": "#d62728",  # red
    "highlight": "#e377c2",   # pink
    "muted": "#7f7f7f",       # gray
}

# Month abbreviations for water year (Oct-Sep)
WATER_YEAR_MONTHS = [
    "Oct", "Nov", "Dec", "Jan", "Feb", "Mar",
    "Apr", "May", "Jun", "Jul", "Aug", "Sep",
]

CALENDAR_MONTHS = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
]


def apply_style():
    """Apply consistent matplotlib rcParams for MOEA-FIND figures."""
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.titlesize": 12,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
