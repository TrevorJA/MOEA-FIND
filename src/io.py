"""I/O helpers for MOEA-FIND experiment artifacts.

Loaders for the file types produced and consumed across workflow stages:
- Pareto ``results.json`` from MOEA runs (stage 04, 05)
- ``metric_bank.parquet`` / ``.csv`` from Pywr-DRB re-evaluation (stage 06)
- Per-realization drought characteristics (``pareto_chars`` /
  ``selected_chars``) consumed by stage 07 scenario discovery

Also provides ``save_experiment_config`` for the ``config.json`` invocation
record that every driver writes alongside its outputs.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

import numpy as np


def load_pareto_results(results_json_path: Path) -> dict:
    """Load a MOEA results.json into a dict.

    Args:
        results_json_path: Path to a ``results.json`` written by stage 04
            (``run_moea_find``) or stage 05 (``run_drb_multisite``).

    Returns:
        Parsed JSON payload.
    """
    return json.loads(Path(results_json_path).read_text())


def load_pareto_metrics(results_json_path: Path) -> np.ndarray:
    """Return the ``drought_metrics`` array from a Pareto results.json.

    Returns:
        ``(n_pareto, n_objectives)`` float array. Empty array if no Pareto
        solutions were found.
    """
    payload = load_pareto_results(results_json_path)
    return np.asarray(payload.get("drought_metrics", []), dtype=float)


def load_pareto_chars(
    chars_path: Path,
    feature_cols: Optional[Iterable[str]] = None,
):
    """Load per-realization drought characteristics as a DataFrame.

    Reads a stage-04 ``results.json`` (which has ``pareto_chars``) or a
    stage-06 selection JSON (which has ``selected_chars``). The DataFrame
    is indexed by string ``realization_id``.

    Args:
        chars_path: Path to a JSON file with ``pareto_chars`` or
            ``selected_chars``.
        feature_cols: If provided, raise ``KeyError`` when any of these
            columns is missing from the loaded characteristics.

    Returns:
        DataFrame indexed by ``realization_id``.

    Raises:
        ValueError: Neither ``pareto_chars`` nor ``selected_chars`` is set.
        KeyError: A requested feature column is not in the loaded data.
    """
    import pandas as pd

    chars_path = Path(chars_path)
    payload = json.loads(chars_path.read_text())
    rows = payload.get("pareto_chars") or payload.get("selected_chars") or []
    if not rows:
        raise ValueError(
            f"{chars_path} has neither 'pareto_chars' nor 'selected_chars'."
        )
    df = pd.DataFrame(rows)
    df.index = df.index.astype(str)
    df.index.name = "realization_id"
    if feature_cols is not None:
        missing = [c for c in feature_cols if c not in df.columns]
        if missing:
            raise KeyError(
                f"drought characteristics in {chars_path} missing columns "
                f"{missing}. Available: {list(df.columns)}"
            )
    return df


def load_metric_bank(path: Path):
    """Load a ``metric_bank.parquet`` (or ``.csv``) DataFrame.

    Falls back from parquet to CSV if a parquet read fails (e.g. pyarrow
    not installed). The returned DataFrame is indexed by string
    ``realization_id``.

    Args:
        path: Path to ``metric_bank.parquet`` or ``metric_bank.csv``. If the
            given suffix is missing on disk, the sibling with the other
            suffix is tried.

    Returns:
        DataFrame indexed by ``realization_id``.
    """
    import pandas as pd

    path = Path(path)
    if path.suffix == ".parquet":
        try:
            df = pd.read_parquet(path)
            return df.assign(
                realization_id=lambda d: d.index.astype(str)
            ).set_index("realization_id", drop=True)
        except Exception:
            csv_alt = path.with_suffix(".csv")
            if csv_alt.exists():
                return pd.read_csv(
                    csv_alt, index_col="realization_id",
                    dtype={"realization_id": str},
                )
            raise
    if path.suffix == ".csv":
        return pd.read_csv(
            path, index_col="realization_id", dtype={"realization_id": str},
        )
    parquet = path.with_suffix(".parquet")
    csv = path.with_suffix(".csv")
    if parquet.exists():
        return load_metric_bank(parquet)
    if csv.exists():
        return load_metric_bank(csv)
    raise FileNotFoundError(
        f"No metric_bank.{{parquet,csv}} at or near {path}"
    )


def save_experiment_config(
    output_dir: Path,
    config: Mapping[str, Any],
    *,
    filename: str = "config.json",
) -> Path:
    """Write a script's invocation record to ``output_dir/config.json``.

    Path objects in ``config`` are converted to strings; numpy scalars are
    converted to Python natives. The output directory is created if needed.

    Args:
        output_dir: Directory to write into.
        config: Mapping of invocation parameters.
        filename: Filename inside ``output_dir`` (defaults to
            ``config.json``).

    Returns:
        The path the config was written to.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / filename
    out_path.write_text(json.dumps(config, indent=2, default=str))
    return out_path
