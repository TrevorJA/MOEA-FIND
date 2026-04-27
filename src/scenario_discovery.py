"""Scenario discovery for the DRB policy re-evaluation pipeline (Stage 4).

Extracts FFMP drought levels from Pywr-DRB simulation outputs, classifies
each Pareto solution as satisficing (no Drought Emergency) or non-satisficing
(Level 6 triggered), and produces the scenario discovery scatter plot in
drought characteristic space.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import json
import math
import numpy as np
import pandas as pd


# ===================================================================
# FFMP Level 6 extraction
# ===================================================================


def extract_drought_levels(
    output_dir: Path,
    realization_ids: List[str],
) -> Dict[str, dict]:
    """Extract FFMP drought levels from Pywr-DRB batch output files.

    Scans all ``batch_*.hdf5`` files in *output_dir* and loads the
    ``drought_level_agg_nyc`` parameter for each realization.

    FFMP Level 6 (Drought Emergency) is indicated by
    ``drought_level_agg_nyc == 6`` in the Pywr-DRB output.

    Args:
        output_dir: Directory containing ``batch_*.hdf5`` files.
        realization_ids: List of realization ID strings to extract.

    Returns:
        ``{real_id: {max_level, level6_triggered, n_days_level6,
        first_level6_date}}`` dict.
    """
    import pywrdrb

    output_dir = Path(output_dir)
    # run_pywrdrb_batch combines and deletes the per-rank batch files,
    # so we look for the combined output first and fall back to any
    # remaining ``*.hdf5`` files for backwards compatibility.
    combined = output_dir / "pywrdrb_output.hdf5"
    if combined.exists():
        output_files = [combined]
    else:
        output_files = sorted(output_dir.glob("*.hdf5"))
    if not output_files:
        raise FileNotFoundError(
            f"No Pywr-DRB output HDF5 found in {output_dir}")

    data = pywrdrb.Data()
    data.load_output(
        output_filenames=[str(f) for f in output_files],
        results_sets=["res_level"],
    )

    results: Dict[str, dict] = {}
    requested = set(realization_ids)

    for dataset_name, scenarios in data.res_level.items():
        for scenario_idx, scenario_df in scenarios.items():
            # scenario_idx may be int or str depending on pywrdrb version
            real_id = str(scenario_idx)
            if real_id not in requested:
                continue

            # Look for drought_level_agg_nyc column
            if "drought_level_agg_nyc" in scenario_df.columns:
                levels = scenario_df["drought_level_agg_nyc"]
            elif "nyc" in scenario_df.columns:
                levels = scenario_df["nyc"]
            else:
                # Try to find any column with drought level data
                level_cols = [c for c in scenario_df.columns
                              if "drought" in c.lower() or "level" in c.lower()]
                if level_cols:
                    levels = scenario_df[level_cols[0]]
                else:
                    print(f"[scenario_discovery] WARNING: no drought level "
                          f"column found for realization {real_id}")
                    results[real_id] = {
                        "max_level": -1,
                        "level6_triggered": False,
                        "n_days_level6": 0,
                        "first_level6_date": None,
                        "warning": "no drought level column found",
                    }
                    continue

            max_level = int(levels.max())
            level6_mask = levels >= 6
            n_days_level6 = int(level6_mask.sum())
            triggered = n_days_level6 > 0
            first_date = (str(levels.index[level6_mask][0])
                          if triggered else None)

            results[real_id] = {
                "max_level": max_level,
                "level6_triggered": triggered,
                "n_days_level6": n_days_level6,
                "first_level6_date": first_date,
            }

    # Flag any missing realizations
    found = set(results.keys())
    missing = requested - found
    if missing:
        print(f"[scenario_discovery] WARNING: {len(missing)} realizations "
              f"not found in output: {sorted(missing)[:5]}...")

    print(f"[scenario_discovery] Extracted drought levels for "
          f"{len(results)} realizations")
    return results


# ===================================================================
# Satisficing classification
# ===================================================================


def build_satisficing_table(
    drought_levels: Dict[str, dict],
    pareto_chars: List[dict],
    drought_metrics: np.ndarray,
    objective_keys: Tuple[str, ...],
) -> pd.DataFrame:
    """Map satisficing classification back to drought characteristic space.

    Satisficing criterion: a scenario is satisficing if and only if it
    does NOT trigger FFMP Drought Emergency (Level 6) at any timestep.

    Args:
        drought_levels: Output from :func:`extract_drought_levels`.
        pareto_chars: List of drought characteristic dicts from Script 04
            results (``results["pareto_chars"]``).
        drought_metrics: Array of shape ``(n_pareto, n_objectives)`` from
            Script 04 results (``results["drought_metrics"]``).
        objective_keys: Tuple of objective names
            (e.g., ``("mean_duration", "mean_avg_severity")``).

    Returns:
        DataFrame with columns: ``pareto_idx``, ``satisficing``,
        each drought characteristic, each objective, and drought level
        details.
    """
    records = []
    for i, (chars, metrics) in enumerate(
        zip(pareto_chars, drought_metrics)
    ):
        real_id = str(i)
        level_info = drought_levels.get(real_id, {})

        record = {
            "pareto_idx": i,
            "satisficing": not level_info.get("level6_triggered", True),
            "max_level": level_info.get("max_level", -1),
            "n_days_level6": level_info.get("n_days_level6", 0),
            "first_level6_date": level_info.get("first_level6_date"),
        }

        # Add all drought characteristics from SSI analysis
        if isinstance(chars, dict):
            for key, val in chars.items():
                record[key] = val

        # Add objective values
        for j, key in enumerate(objective_keys):
            record[f"obj_{key}"] = float(metrics[j])

        records.append(record)

    df = pd.DataFrame(records)

    n_satisficing = df["satisficing"].sum()
    n_total = len(df)
    print(f"[scenario_discovery] Classification: {n_satisficing}/{n_total} "
          f"satisficing ({100 * n_satisficing / max(n_total, 1):.1f}%)")

    return df


# ===================================================================
# Plotting
# ===================================================================


def _data_axis_limits(
    values: np.ndarray, pad_frac: float = 0.05,
) -> Tuple[float, float]:
    """Axis limits from data min/max with symmetric fractional padding."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return (0.0, 1.0)
    lo, hi = float(v.min()), float(v.max())
    span = hi - lo if hi > lo else max(abs(hi), 1.0)
    pad = pad_frac * span
    return (lo - pad, hi + pad)


def _anti_ideal_offplot_corner(
    anti_ideal: Optional[np.ndarray],
    xlim: Tuple[float, float],
    ylim: Tuple[float, float],
) -> Optional[str]:
    """Return 'upper right' / 'lower left' / ... if anti-ideal is outside the
    axes, else None."""
    if anti_ideal is None or len(anti_ideal) < 2:
        return None
    ax_x, ax_y = float(anti_ideal[0]), float(anti_ideal[1])
    if xlim[0] <= ax_x <= xlim[1] and ylim[0] <= ax_y <= ylim[1]:
        return None
    horiz = "right" if ax_x > xlim[1] else ("left" if ax_x < xlim[0] else "right")
    vert = "upper" if ax_y > ylim[1] else ("lower" if ax_y < ylim[0] else "upper")
    return f"{vert} {horiz}"


def plot_satisficing_map(
    df: pd.DataFrame,
    hist_chars: Optional[dict] = None,
    anti_ideal: Optional[np.ndarray] = None,
    output_path: Optional[Path] = None,
    x_col: str = "mean_severity",
    y_col: str = "mean_magnitude",
    x_label: str = "Mean drought severity (|SSI|)",
    y_label: str = "Mean event cumulative deficit (|SSI|·month)",
):
    """Legacy single-panel scatter used when no manifest is supplied.

    Axes are cropped to the data cloud (+5% pad), so the anti-ideal
    point is deliberately off-plot when it sits outside the realized
    drought-characteristic cloud — a corner annotation flags this.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    satisficing = df[df["satisficing"]]
    non_satisficing = df[~df["satisficing"]]

    ax.scatter(
        satisficing[x_col], satisficing[y_col],
        c="#2ca02c", s=8, alpha=0.6, edgecolors="none",
        label=f"Satisficing (n={len(satisficing)})",
    )
    ax.scatter(
        non_satisficing[x_col], non_satisficing[y_col],
        c="#d62728", s=8, alpha=0.6, edgecolors="none",
        label=f"Non-satisficing (n={len(non_satisficing)})",
    )

    if hist_chars is not None and x_col in hist_chars and y_col in hist_chars:
        ax.scatter(
            hist_chars[x_col], hist_chars[y_col],
            marker="*", s=250, c="black", zorder=5,
            label="Historical",
        )

    xlim = _data_axis_limits(df[x_col].values)
    ylim = _data_axis_limits(df[y_col].values)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    corner = _anti_ideal_offplot_corner(anti_ideal, xlim, ylim)
    if corner is not None:
        vert, horiz = corner.split()
        ax.text(
            0.98 if horiz == "right" else 0.02,
            0.98 if vert == "upper" else 0.02,
            "D* anti-ideal off-plot",
            transform=ax.transAxes,
            ha=horiz, va=("top" if vert == "upper" else "bottom"),
            fontsize=8, color="gray", style="italic",
        )

    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title("Scenario Discovery: DRB Policy Re-evaluation", fontsize=13)
    ax.legend(fontsize=10, loc="best")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[scenario_discovery] Saved {output_path}")

    return fig


def plot_satisficing_map_multi(
    labels_long: pd.DataFrame,
    chars_df: pd.DataFrame,
    manifest: Sequence,
    classifiers_dir: Path,
    feature_cols: Optional[Sequence[str]] = None,
    x_col: str = "mean_severity",
    y_col: str = "mean_magnitude",
    agg_col: Optional[str] = "time_in_drought_fraction",
    agg_fn: str = "median",
    x_label: str = "Mean drought severity (|SSI|)",
    y_label: str = "Mean event cumulative deficit (|SSI|·month)",
    hist_chars: Optional[dict] = None,
    anti_ideal: Optional[np.ndarray] = None,
    output_path: Optional[Path] = None,
    n_grid: int = 80,
    n_cols: int = 3,
    classifier_label: str = (
        "sklearn GradientBoostingClassifier, default hyperparams, 5-fold stratified CV"
    ),
):
    """Manuscript Figure 9 (multi-panel): one subplot per satisficing
    criterion, each with a classifier-derived P(satisficing) background.

    Args:
        labels_long: Output of
            :func:`src.satisficing_labels.apply_labels(bank_df, manifest)`.
            Columns include ``definition_id``, ``realization_id``, ``y``.
        chars_df: Drought-characteristic features indexed by
            ``realization_id``. Must contain every column in
            *feature_cols*.
        manifest: Sequence of :class:`~src.satisficing_labels.SatisficingDefinition`.
        classifiers_dir: Directory holding one subdir per definition
            (``<definition_id>/model.joblib``) as written by
            :func:`src.satisficing_labels.sweep_manifest`.
        feature_cols: Feature ordering used when the classifiers were
            trained. The prediction grid preserves this ordering.
        x_col / y_col: Drought characteristics to place on each panel's
            axes.
        agg_col: Third feature held at its aggregate (median/mean) across
            the grid. Set to None if classifiers were trained on two
            features only.
        agg_fn: "median" (default) or "mean" for the held-fixed feature.
        hist_chars: Historical drought-characteristic dict (plotted as
            black star per panel if both axes are present).
        anti_ideal: Anti-ideal point. Deliberately clipped off-plot by
            the data-driven axis limits; a corner note flags this.
        output_path: PDF destination. If None, the figure is returned
            without saving.
        n_grid: Resolution of the P(satisficing) evaluation grid.
        n_cols: Number of panel columns.

    Returns:
        matplotlib Figure object.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import joblib

    n_defs = len(manifest)
    if n_defs == 0:
        raise ValueError("plot_satisficing_map_multi: empty manifest")

    aggregator = np.median if agg_fn == "median" else np.mean
    # Hadjimichael et al. 2020 aesthetics: RdBu diverging map (blue =
    # high P(satisficing), red = low), with matching scatter colors.
    boundary_cmap = plt.get_cmap("RdBu")
    sat_color = "#4575b4"     # muted blue (RdBu high end)
    fail_color = "#c0392b"    # muted red (RdBu low end)

    # Feature columns we need to evaluate classifiers. Grid values are
    # filled from chars_df aggregates for columns other than x/y. Default
    # to the primary metric preset.
    if feature_cols is None:
        from src.drought_metrics import metric_names, resolve_metric_set
        feature_cols = list(metric_names(resolve_metric_set("primary")))
    else:
        feature_cols = list(feature_cols)
    missing = [c for c in feature_cols if c not in chars_df.columns]
    if missing:
        raise KeyError(
            f"chars_df missing feature columns {missing}; "
            f"needed to evaluate classifiers"
        )

    # Global axis limits so every panel shares the same frame of reference
    xlim = _data_axis_limits(chars_df[x_col].values)
    ylim = _data_axis_limits(chars_df[y_col].values)

    xx, yy = np.meshgrid(
        np.linspace(xlim[0], xlim[1], n_grid),
        np.linspace(ylim[0], ylim[1], n_grid),
    )
    grid = pd.DataFrame({x_col: xx.ravel(), y_col: yy.ravel()})
    for col in feature_cols:
        if col in (x_col, y_col):
            continue
        grid[col] = float(aggregator(pd.to_numeric(chars_df[col], errors="coerce").dropna()))
    grid = grid[feature_cols]

    # Figure layout: panels + bottom legend row + right colorbar column
    n_rows = math.ceil(n_defs / n_cols)
    panel_w, panel_h = 4.2, 3.6
    fig_w = panel_w * n_cols + 1.3  # +colorbar strip
    fig_h = panel_h * n_rows + 1.1 + 0.25 * n_defs  # +legend strip
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        nrows=n_rows + 1, ncols=n_cols + 1,
        width_ratios=[1.0] * n_cols + [0.06],
        height_ratios=[1.0] * n_rows + [0.12 + 0.04 * n_defs],
        wspace=0.28, hspace=0.45,
    )

    axes_for_mesh = []
    for idx, d in enumerate(manifest):
        row, col = divmod(idx, n_cols)
        ax = fig.add_subplot(gs[row, col])

        # Background heatmap from the trained classifier
        model_path = Path(classifiers_dir) / d.definition_id / "model.joblib"
        mesh = None
        if model_path.exists():
            try:
                model = joblib.load(model_path)
                probs = model.predict_proba(grid.values)[:, 1].reshape(xx.shape)
                mesh = ax.contourf(
                    xx, yy, probs, levels=np.linspace(0, 1, 11),
                    cmap=boundary_cmap, alpha=0.85,
                )
            except Exception as exc:
                ax.text(
                    0.5, 0.5, f"model load failed\n({exc})",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="gray",
                )
        else:
            ax.text(
                0.5, 0.5, "no classifier",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color="gray",
            )
        if mesh is not None:
            axes_for_mesh.append(mesh)

        # Realized outcomes — scatter coloured by y
        sub = labels_long[labels_long["definition_id"] == d.definition_id]
        if not sub.empty:
            real_ids = sub["realization_id"].astype(str).values
            keep = np.array([r in chars_df.index for r in real_ids])
            real_ids = real_ids[keep]
            y_vec = sub["y"].values.astype(int)[keep]
            x_vals = chars_df.loc[real_ids, x_col].astype(float).values
            y_vals = chars_df.loc[real_ids, y_col].astype(float).values
            ax.scatter(
                x_vals[y_vec == 1], y_vals[y_vec == 1],
                c=sat_color, s=8, alpha=0.8, edgecolors="none",
            )
            ax.scatter(
                x_vals[y_vec == 0], y_vals[y_vec == 0],
                c=fail_color, s=8, alpha=0.8, edgecolors="none",
            )
            n_pos = int((y_vec == 1).sum())
            n_neg = int((y_vec == 0).sum())
        else:
            n_pos = n_neg = 0

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        corner = _anti_ideal_offplot_corner(anti_ideal, xlim, ylim)
        if corner is not None:
            vert, horiz = corner.split()
            ax.text(
                0.98 if horiz == "right" else 0.02,
                0.98 if vert == "upper" else 0.02,
                "D* off-plot",
                transform=ax.transAxes,
                ha=horiz, va=("top" if vert == "upper" else "bottom"),
                fontsize=7, color="gray", style="italic",
            )

        if row == n_rows - 1 or idx >= n_defs - n_cols:
            ax.set_xlabel(x_label, fontsize=10)
        if col == 0:
            ax.set_ylabel(y_label, fontsize=10)
        ax.tick_params(axis="both", labelsize=8)

        # Panel letter + definition_id in the title (linked to the
        # extended legend below). n_sat/n_fail counts go inside the axes
        # in the corner opposite the anti-ideal flag so they don't
        # overlap the title.
        panel_letter = chr(ord("a") + idx)
        ax.set_title(
            f"({panel_letter}) {d.definition_id}",
            fontsize=10, pad=4,
        )
        # Place counts in the lower-left when anti-ideal is in the
        # upper-right (typical case); otherwise upper-left.
        count_vert = "bottom" if (corner is None or "upper" in corner) else "top"
        ax.text(
            0.02, 0.02 if count_vert == "bottom" else 0.98,
            f"n={n_pos}/{n_pos + n_neg} sat.",
            transform=ax.transAxes, ha="left", va=count_vert,
            fontsize=7, color="black",
            bbox=dict(facecolor="white", edgecolor="none",
                      alpha=0.7, pad=1.5),
        )

    # Shared colorbar on the right, spanning all panel rows
    cax = fig.add_subplot(gs[:n_rows, -1])
    if axes_for_mesh:
        cbar = fig.colorbar(axes_for_mesh[0], cax=cax)
        cbar.set_label("Probability of satisficing", fontsize=10)
        cbar.ax.tick_params(labelsize=8)
    else:
        cax.set_axis_off()

    # Extended legend spanning the bottom row. Each entry is keyed by
    # the matching panel letter so a reader can map (a) → first subplot.
    legend_ax = fig.add_subplot(gs[-1, :])
    legend_ax.set_axis_off()
    lines = [
        f"Satisficing criteria (classifier: {classifier_label}; "
        f"features: {', '.join(feature_cols)}):",
    ]
    for idx, d in enumerate(manifest):
        letter = chr(ord("a") + idx)
        direction_symbol = {"le": "≤", "lt": "<", "ge": "≥", "gt": ">"}[
            d.direction
        ]
        desc = d.description or "(no description)"
        lines.append(
            f"  ({letter}) {d.definition_id}: {desc}  "
            f"[{d.metric} {direction_symbol} {d.threshold}]"
        )
    legend_ax.text(
        0.0, 1.0, "\n".join(lines),
        transform=legend_ax.transAxes, ha="left", va="top",
        fontsize=8, family="monospace",
    )

    fig.suptitle(
        "Scenario Discovery: DRB Policy Re-evaluation",
        fontsize=13, y=0.995,
    )

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[scenario_discovery] Saved {output_path}")

    return fig


def save_results(
    df: pd.DataFrame,
    drought_levels: Dict[str, dict],
    output_dir: Path,
) -> None:
    """Save classification results to disk.

    Writes:
        - ``satisficing_classification.csv``
        - ``drought_levels.json``
        - ``satisficing_summary.json``

    Args:
        df: Output from :func:`build_satisficing_table`.
        drought_levels: Output from :func:`extract_drought_levels`.
        output_dir: Results directory.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_dir / "satisficing_classification.csv", index=False)

    with open(output_dir / "drought_levels.json", "w") as f:
        json.dump(drought_levels, f, indent=2, default=str)

    summary = {
        "n_total": len(df),
        "n_satisficing": int(df["satisficing"].sum()),
        "n_non_satisficing": int((~df["satisficing"]).sum()),
        "pct_satisficing": float(
            100 * df["satisficing"].sum() / max(len(df), 1)
        ),
        "max_level_distribution": (
            df["max_level"].value_counts().sort_index().to_dict()
        ),
    }
    with open(output_dir / "satisficing_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[scenario_discovery] Saved results to {output_dir}")
