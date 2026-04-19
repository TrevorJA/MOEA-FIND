"""Scenario discovery for the DRB policy re-evaluation pipeline (Stage 4).

Extracts FFMP drought levels from Pywr-DRB simulation outputs, classifies
each Pareto solution as satisficing (no Drought Emergency) or non-satisficing
(Level 6 triggered), and produces the scenario discovery scatter plot in
drought characteristic space.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import json
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
    batch_files = sorted(output_dir.glob("batch_*.hdf5"))

    if not batch_files:
        raise FileNotFoundError(
            f"No batch_*.hdf5 files found in {output_dir}")

    # Load all batch outputs
    data = pywrdrb.Data()
    data.load_output(
        output_filenames=[str(f) for f in batch_files],
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


def plot_satisficing_map(
    df: pd.DataFrame,
    hist_chars: Optional[dict] = None,
    anti_ideal: Optional[np.ndarray] = None,
    output_path: Optional[Path] = None,
    x_col: str = "mean_duration",
    y_col: str = "mean_avg_severity",
    x_label: str = "Mean drought duration (months)",
    y_label: str = "Mean avg. severity (|SSI|)",
):
    """Manuscript Figure 9: Drought characteristic space colored by satisficing.

    Green markers = satisficing (no Level 6 triggered).
    Red markers = non-satisficing (Level 6 triggered at least once).

    Args:
        df: Output from :func:`build_satisficing_table`.
        hist_chars: Historical drought characteristics dict (optional,
            plotted as a star).
        anti_ideal: Anti-ideal point (optional, plotted as an X).
        output_path: Path to save figure (PDF). If None, figure is
            returned without saving.
        x_col: Column name for x-axis drought characteristic.
        y_col: Column name for y-axis drought characteristic.
        x_label: Display label for x-axis.
        y_label: Display label for y-axis.

    Returns:
        matplotlib Figure object.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    satisficing = df[df["satisficing"]]
    non_satisficing = df[~df["satisficing"]]

    ax.scatter(
        satisficing[x_col], satisficing[y_col],
        c="#2ca02c", s=30, alpha=0.7, edgecolors="none",
        label=f"Satisficing (n={len(satisficing)})",
    )
    ax.scatter(
        non_satisficing[x_col], non_satisficing[y_col],
        c="#d62728", s=30, alpha=0.7, edgecolors="none",
        label=f"Non-satisficing (n={len(non_satisficing)})",
    )

    if hist_chars is not None and x_col in hist_chars and y_col in hist_chars:
        ax.scatter(
            hist_chars[x_col], hist_chars[y_col],
            marker="*", s=250, c="black", zorder=5,
            label="Historical",
        )

    if anti_ideal is not None and len(anti_ideal) >= 2:
        ax.scatter(
            anti_ideal[0], anti_ideal[1],
            marker="X", s=200, c="orange", zorder=5, edgecolors="black",
            label="Anti-ideal D*",
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
