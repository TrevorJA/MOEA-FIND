"""Satisficing metric bank for Pywr-DRB ensemble outputs.

Computes a fixed panel of candidate satisficing metrics per realization
from a Pywr-DRB combined HDF5 output file. The *metric bank* is the
reproducibility boundary described in the project plan: the expensive
Pywr-DRB simulation runs once, producing this bank; downstream label
manifests and GBT classifiers (see :mod:`src.discovery.satisficing_labels`) are
re-runnable cheaply without re-simulating.

Metric families
---------------

1. FFMP drought-level exposure (from ``res_level``)
2. NYC combined-storage minima + drawdown (from ``res_storage``)
3. Hashimoto reliability / resilience / vulnerability on NYC delivery
   shortages (from ``ibt_diversions`` + ``ibt_demands``)
4. Flow-target reliability at Montague and Trenton (from ``major_flow``
   + ``mrf_target``)

The extractors tolerate missing results_sets: a missing set emits NaN
for its columns plus one warning per set, so the bank remains meaningful
even when a run is missing a variable. Column names and conventions
match ``../StochasticExploratoryExperiment/methods/metrics/shortfall.py``
so the two projects can compare results directly.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


# Default NYC reservoir capacities (million gallons). Used to normalise
# min-storage metrics. Matches
# ``../StochasticExploratoryExperiment/methods/config.py::NYC_STORAGE_CAPACITIES``.
NYC_STORAGE_CAPACITIES = {
    "cannonsville": 95700.0,
    "pepacton": 140200.0,
    "neversink": 34900.0,
}
NYC_TOTAL_CAPACITY = float(sum(NYC_STORAGE_CAPACITIES.values()))

# FFMP drought levels tracked in the bank. 0 = normal; higher = deeper.
_FFMP_LEVELS = (1, 2, 3, 6)

# Flow-target nodes
_FLOW_TARGET_NODES = ("delMontague", "delTrenton")

# Default shortage tolerance (MGD); sub-tolerance deficits are zeroed.
_DEFAULT_SHORTAGE_TOLERANCE_MGD = 1.0

# Shortfall break length (days) for Hashimoto event detection.
_SHORTFALL_BREAK_LENGTH_DAYS = 7


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _safe_first_column(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[pd.Series]:
    """Return the first existing column in *df* from *candidates*, else None."""
    for col in candidates:
        if col in df.columns:
            return df[col]
    return None


def _longest_run(mask: np.ndarray) -> int:
    """Longest run of ``True`` in a 1D bool array."""
    m = np.asarray(mask, dtype=bool)
    if not m.any():
        return 0
    # Pad with False at both ends and diff
    padded = np.concatenate(([False], m, [False]))
    diffs = np.diff(padded.astype(np.int8))
    starts = np.where(diffs == 1)[0]
    ends = np.where(diffs == -1)[0]
    return int((ends - starts).max()) if len(starts) else 0


def _hashimoto(
    flow: pd.Series,
    target: pd.Series,
    tolerance: float = _DEFAULT_SHORTAGE_TOLERANCE_MGD,
    break_length: int = _SHORTFALL_BREAK_LENGTH_DAYS,
) -> Dict[str, float]:
    """Hashimoto reliability / resilience / vulnerability on a deficit series.

    Mirrors ``methods/metrics/shortfall.py::calculate_hashimoto_metrics`` in
    the reference project. Returns reliability, resiliency (probability of
    recovery given current deficit), vulnerability (maximum daily deficit
    normalised by mean target during events), max_event_duration, and
    event_count.
    """
    if flow is None or target is None:
        return {
            "reliability": np.nan,
            "resiliency": np.nan,
            "vulnerability": np.nan,
            "max_event_duration_days": np.nan,
            "event_count": np.nan,
        }

    f = np.asarray(flow.values, dtype=float)
    t = np.asarray(target.values, dtype=float)
    deficits = np.maximum(t - f, 0.0)
    is_deficit = deficits >= tolerance

    reliability = float((~is_deficit).mean())

    if reliability >= 1.0 - 1e-9:
        return {
            "reliability": reliability,
            "resiliency": 1.0,
            "vulnerability": 0.0,
            "max_event_duration_days": 0.0,
            "event_count": 0.0,
        }

    # Resiliency — P(recover | currently in deficit). Defined as in
    # Hashimoto et al. (1982): transitions from deficit to non-deficit
    # divided by total deficit-day count.
    transitions = np.logical_and(is_deficit[:-1], ~is_deficit[1:]).sum()
    resiliency = float(transitions) / float(is_deficit[:-1].sum()) \
        if is_deficit[:-1].any() else np.nan

    # Event detection with break_length tolerance
    durations: List[int] = []
    vulnerabilities: List[float] = []
    targets_sum = 0.0
    duration = 0
    vulnerability = 0.0
    in_event = False
    target_sum_event = 0.0
    n = len(f)
    for i in range(n):
        if in_event or is_deficit[i]:
            if not in_event:
                in_event = True
                duration = 0
                vulnerability = 0.0
                target_sum_event = 0.0
            duration += 1
            target_sum_event += t[i]
            vulnerability = max(vulnerability, float(deficits[i]))
            look_ahead = is_deficit[i + 1 : i + 1 + break_length]
            if not np.any(look_ahead):
                durations.append(duration)
                targets_sum += target_sum_event
                if target_sum_event > 0:
                    vulnerabilities.append(
                        float(vulnerability) / float(target_sum_event / duration)
                    )
                else:
                    vulnerabilities.append(float(vulnerability))
                in_event = False
                duration = 0
                vulnerability = 0.0

    max_vul = float(max(vulnerabilities)) if vulnerabilities else 0.0
    max_dur = float(max(durations)) if durations else 0.0
    return {
        "reliability": reliability,
        "resiliency": resiliency,
        "vulnerability": max_vul,
        "max_event_duration_days": max_dur,
        "event_count": float(len(durations)),
    }


# ---------------------------------------------------------------------------
# Per-realization extractors
# ---------------------------------------------------------------------------


def _ffmp_exposure(res_level_df: Optional[pd.DataFrame]) -> Dict[str, float]:
    """FFMP drought-level exposure statistics from ``res_level``.

    ``res_level_df`` is the per-realization DataFrame from
    ``data.res_level[dataset_id][realization_id]``. The aggregated NYC
    level is tried under both common column names.
    """
    if res_level_df is None:
        return {
            "max_level": np.nan,
            "first_L6_day": np.nan,
            **{f"days_L{l}": np.nan for l in _FFMP_LEVELS},
            "days_drought_any": np.nan,
        }

    levels = _safe_first_column(res_level_df, ("drought_level_agg_nyc", "nyc"))
    if levels is None:
        return {
            "max_level": np.nan,
            "first_L6_day": np.nan,
            **{f"days_L{l}": np.nan for l in _FFMP_LEVELS},
            "days_drought_any": np.nan,
        }

    lv = np.asarray(levels.values, dtype=int)
    stats = {
        "max_level": int(lv.max()),
        "days_drought_any": int((lv > 0).sum()),
    }
    for level in _FFMP_LEVELS:
        stats[f"days_L{level}"] = int((lv >= level).sum())

    first_l6 = np.argmax(lv >= 6) if (lv >= 6).any() else -1
    stats["first_L6_day"] = int(first_l6) if first_l6 >= 0 else np.nan
    return stats


def _storage_minima(res_storage_df: Optional[pd.DataFrame]) -> Dict[str, float]:
    """Min storage fraction + max drawdown from ``res_storage``."""
    base_keys = [f"{r}_min_storage_frac" for r in NYC_STORAGE_CAPACITIES]
    out = {k: np.nan for k in base_keys}
    out["nyc_min_storage_frac"] = np.nan
    out["nyc_drawdown_days_below_0.5"] = np.nan
    out["nyc_drawdown_days_below_0.25"] = np.nan

    if res_storage_df is None:
        return out

    # Per-reservoir minima
    nyc_total = np.zeros(len(res_storage_df), dtype=float)
    any_series = False
    for res, cap in NYC_STORAGE_CAPACITIES.items():
        if res in res_storage_df.columns:
            s = res_storage_df[res].values.astype(float)
            out[f"{res}_min_storage_frac"] = float(s.min() / cap)
            nyc_total += s
            any_series = True
    if any_series:
        frac = nyc_total / NYC_TOTAL_CAPACITY
        out["nyc_min_storage_frac"] = float(frac.min())
        out["nyc_drawdown_days_below_0.5"] = int(_longest_run(frac < 0.5))
        out["nyc_drawdown_days_below_0.25"] = int(_longest_run(frac < 0.25))
    return out


def _nyc_delivery_shortage(
    ibt_div_df: Optional[pd.DataFrame],
    ibt_dem_df: Optional[pd.DataFrame],
) -> Dict[str, float]:
    """Hashimoto on NYC delivery vs demand."""
    prefix = "nyc_delivery"
    if ibt_div_df is None or ibt_dem_df is None:
        return {f"{prefix}_{k}": np.nan
                for k in ("reliability", "resiliency", "vulnerability",
                           "max_event_duration_days", "event_count")}

    flow = _safe_first_column(ibt_div_df, ("delivery_nyc",))
    target = _safe_first_column(ibt_dem_df, ("demand_nyc",))
    h = _hashimoto(flow, target)
    return {f"{prefix}_{k}": v for k, v in h.items()}


def _flow_target_reliability(
    major_flow_df: Optional[pd.DataFrame],
    mrf_target_df: Optional[pd.DataFrame],
) -> Dict[str, float]:
    """Reliability at Montague and Trenton flow targets."""
    out = {}
    for node in _FLOW_TARGET_NODES:
        prefix = f"{node.lower().replace('del', '')}_flow"
        if major_flow_df is None or mrf_target_df is None:
            out[f"{prefix}_reliability"] = np.nan
            out[f"{prefix}_vulnerability"] = np.nan
            continue
        flow = major_flow_df.get(node)
        if node == "delTrenton" and "delTrenton_equiv" in major_flow_df.columns:
            flow = major_flow_df["delTrenton_equiv"]
        target = mrf_target_df.get(node)
        if flow is None or target is None:
            out[f"{prefix}_reliability"] = np.nan
            out[f"{prefix}_vulnerability"] = np.nan
            continue
        h = _hashimoto(flow, target)
        out[f"{prefix}_reliability"] = h["reliability"]
        out[f"{prefix}_vulnerability"] = h["vulnerability"]
    return out


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def _load_data(
    output_file: Path,
    results_sets: Sequence[str],
) -> "pywrdrb.Data":
    """Load a pywrdrb.Data object for the requested results sets.

    Silent to a missing results set: pywrdrb raises, which we catch and
    warn about so downstream extractors emit NaN rather than crash.
    """
    import pywrdrb
    data = pywrdrb.Data()
    try:
        data.load_output(
            output_filenames=[str(output_file)],
            results_sets=list(results_sets),
        )
    except Exception as exc:
        warnings.warn(
            f"[satisficing_metrics] pywrdrb.Data.load_output raised {exc!r} "
            f"for {output_file}; trying per-set loads",
            stacklevel=2,
        )
        for rs in results_sets:
            try:
                data.load_output(
                    output_filenames=[str(output_file)],
                    results_sets=[rs],
                )
            except Exception as inner:
                warnings.warn(
                    f"[satisficing_metrics] skipping results_set={rs!r} "
                    f"({inner!r})",
                    stacklevel=2,
                )
    return data


def _get_realization_frame(
    data_attr: Dict,
    realization_id: str,
) -> Optional[pd.DataFrame]:
    """Pull ``data.<results_set>[dataset_id][realization_id]`` safely.

    ``data_attr`` is the attribute dict keyed by dataset name (one key
    per registered flow_type). We collapse across dataset keys because
    the MOEA-FIND pipeline has exactly one flow_type per invocation.
    """
    if not data_attr:
        return None
    for _, scenarios in data_attr.items():
        if realization_id in scenarios:
            return scenarios[realization_id]
        # Try int and stringified-int as fallback
        try:
            as_int = int(realization_id)
            if as_int in scenarios:
                return scenarios[as_int]
        except (ValueError, TypeError):
            pass
    return None


def compute_metric_bank(
    output_file: Path,
    realization_ids: Sequence[str],
    results_sets: Sequence[str] = (
        "res_level",
        "res_storage",
        "major_flow",
        "mrf_target",
        "ibt_diversions",
        "ibt_demands",
    ),
    comm=None,
) -> pd.DataFrame:
    """Compute the satisficing metric bank for a combined Pywr-DRB output.

    Args:
        output_file: Combined HDF5 produced by :func:`src.pywrdrb.bridge.run_pywrdrb_batch`.
        realization_ids: String IDs for each realization in the output.
        results_sets: pywrdrb results_sets to load. The defaults cover
            every extractor in this module; subset if memory is tight.
        comm: Optional MPI communicator. If provided and ``size>1``,
            realizations are split across ranks and rank 0 returns the
            gathered DataFrame. Other ranks return an empty DataFrame.

    Returns:
        DataFrame indexed by ``realization_id`` with one column per metric.
    """
    output_file = Path(output_file)
    rank, size = 0, 1
    if comm is not None and comm.Get_size() > 1:
        rank, size = comm.Get_rank(), comm.Get_size()

    # Load Data object once per rank; its internal dicts are shared across
    # realizations. For very large ensembles a streaming load would be
    # lower memory, but at typical scale (~2 k realizations) this is fine.
    data = _load_data(output_file, results_sets)

    all_ids = list(realization_ids)
    rank_ids = np.array_split(all_ids, size)[rank] if size > 1 else all_ids

    rows: List[Dict] = []
    for real_id in rank_ids:
        res_level_df = _get_realization_frame(
            getattr(data, "res_level", {}) or {}, str(real_id),
        )
        res_storage_df = _get_realization_frame(
            getattr(data, "res_storage", {}) or {}, str(real_id),
        )
        major_flow_df = _get_realization_frame(
            getattr(data, "major_flow", {}) or {}, str(real_id),
        )
        mrf_target_df = _get_realization_frame(
            getattr(data, "mrf_target", {}) or {}, str(real_id),
        )
        ibt_div_df = _get_realization_frame(
            getattr(data, "ibt_diversions", {}) or {}, str(real_id),
        )
        ibt_dem_df = _get_realization_frame(
            getattr(data, "ibt_demands", {}) or {}, str(real_id),
        )

        record: Dict = {"realization_id": str(real_id)}
        record.update(_ffmp_exposure(res_level_df))
        record.update(_storage_minima(res_storage_df))
        record.update(_nyc_delivery_shortage(ibt_div_df, ibt_dem_df))
        record.update(_flow_target_reliability(major_flow_df, mrf_target_df))
        rows.append(record)

    local_df = pd.DataFrame(rows)

    if comm is None or size <= 1:
        return local_df.set_index("realization_id", drop=True)

    # Gather to rank 0 via point-to-point (avoids collective barriers on
    # Hopper libfabric; matches ../StochasticExploratoryExperiment
    # /methods/mpi_utils.py::global_point_to_point_gather).
    if rank == 0:
        gathered: List[pd.DataFrame] = [local_df]
        for src in range(1, size):
            piece = comm.recv(source=src, tag=900 + src)
            gathered.append(piece)
        out = pd.concat(gathered, ignore_index=True)
        return out.set_index("realization_id", drop=True)
    else:
        comm.send(local_df, dest=0, tag=900 + rank)
        return pd.DataFrame()


def write_metric_bank(
    df: pd.DataFrame,
    output_path: Path,
) -> Path:
    """Persist the metric bank to parquet (preferred) or CSV fallback.

    Requires ``pyarrow`` for parquet; falls back to CSV with a warning
    if not installed. The reproducibility contract is that this file
    plus the Pareto archive's ``results.json`` is enough to re-run the
    label/classifier layer without any Pywr-DRB step.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        df.to_parquet(output_path)
    except Exception as exc:
        csv_path = output_path.with_suffix(".csv")
        warnings.warn(
            f"[satisficing_metrics] parquet write failed ({exc!r}); "
            f"falling back to {csv_path}",
            stacklevel=2,
        )
        df.to_csv(csv_path)
        return csv_path
    return output_path
