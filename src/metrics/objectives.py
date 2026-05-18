"""Objective functions for MOEA-FIND.

Includes analytic test problems (POC), SSI (Standardized Streamflow Index)
via SynHydro, and drought metric calculations.

SSI calculation is delegated to synhydro.droughts.ssi.SSI, which:
  1. Accumulates monthly flows over a rolling window (timescale parameter)
  2. Fits a distribution (Gamma by default) per calendar month
  3. Transforms to standard normal via the fitted CDF

Drought events are identified by synhydro.droughts.ssi.get_drought_metrics,
which uses SSI < 0 for mild onset, SSI <= -1 for critical drought, and a
configurable recovery hysteresis (default 3 consecutive positive months).

References:
  - Vicente-Serrano et al. (2012) for standardized drought index methodology
  - Zaniolo et al. (2023, FIND) for SSI-based drought control
  - McKee et al. (1993) for SPI (same procedure applied to precipitation)
"""

from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from scipy import stats


def manhattan_norm(x: np.ndarray, anti_ideal: np.ndarray) -> float:
    """L1 (Manhattan) distance from anti-ideal point.

    Args:
        x: Point in objective space (k-dimensional).
        anti_ideal: Anti-ideal point D* (k-dimensional).

    Returns:
        L1 norm ||x - anti_ideal||_1.
    """
    return float(np.sum(np.abs(x - anti_ideal)))


def cyclic_mean(values: np.ndarray, period: float = 12.0) -> float:
    """Arithmetic mean of a cyclic coordinate via sin/cos embedding.

    Handles the December/January wraparound correctly: `cyclic_mean`
    of (1, 12) is ~12.5 mod period, not 6.5. Used for the per-trace
    aggregation of peak-severity months across events.

    Args:
        values: 1D array of cyclic coordinates (e.g. months 1..12).
        period: Period length. Default 12.

    Returns:
        Cyclic mean in [0, period).
    """
    if len(values) == 0:
        return 0.0
    theta = 2.0 * np.pi * np.asarray(values, dtype=float) / period
    mean_sin = np.mean(np.sin(theta))
    mean_cos = np.mean(np.cos(theta))
    angle = np.arctan2(mean_sin, mean_cos)
    if angle < 0:
        angle += 2.0 * np.pi
    return float(angle * period / (2.0 * np.pi))


def analytic_objectives(
    dvs: np.ndarray,
    anti_ideal: np.ndarray,
) -> np.ndarray:
    """Analytic test problem: J_i = X_i for i=1..k, J_{k+1} = ||X - X*||_1.

    Works for any dimensionality k. Decision variables ARE the first k objectives.

    Args:
        dvs: Decision variable vector (k-dimensional).
        anti_ideal: Anti-ideal point X* (k-dimensional).

    Returns:
        Objective vector (k+1 dimensional): [X_1, ..., X_k, manhattan_norm].
    """
    objs = np.empty(len(dvs) + 1)
    objs[:-1] = dvs
    objs[-1] = manhattan_norm(dvs, anti_ideal)
    return objs


def compute_discrepancy(
    points: np.ndarray,
    method: str = "L2-star",
) -> float:
    """Compute discrepancy of a point set in [0,1]^d.

    Args:
        points: Array of shape (n, d), values in [0, 1].
        method: Discrepancy type ("L2-star", "CD", "MD", "WD").

    Returns:
        Discrepancy value (lower = more uniform).
    """
    from scipy.stats.qmc import discrepancy
    return float(discrepancy(points, method=method))


###############################################################################
# SSI via SynHydro
###############################################################################

# Re-export SynHydro's SSI and get_drought_metrics for convenience
from synhydro.droughts.ssi import SSI as SynHydroSSI
from synhydro.droughts.ssi import get_drought_metrics


def make_ssi_calculator(
    timescale: int = 3,
    dist: str = "gamma",
    fit_freq: str = "M",
) -> SynHydroSSI:
    """Create a configured SynHydro SSI calculator.

    Thin factory to set MOEA-FIND defaults (monthly fitting, gamma dist)
    while exposing SynHydro's full configuration.

    Args:
        timescale: Accumulation period in months (1, 3, 6, 12).
        dist: Distribution name ("gamma", "lognorm", "pearson3", etc.).
        fit_freq: Frequency for seasonal fitting ("M" for monthly).

    Returns:
        Configured (unfitted) SynHydroSSI instance.
    """
    return SynHydroSSI(
        dist=dist,
        timescale=timescale,
        fit_freq=fit_freq,
    )


def flows_to_series(
    monthly_flows: np.ndarray,
    start_date: str = "1950-10-01",
    freq: str = "MS",
) -> pd.Series:
    """Convert a numpy array of monthly flows to a pd.Series with DatetimeIndex.

    SynHydro's SSI requires a pd.Series with DatetimeIndex. This helper
    wraps numpy arrays (from generators) with a synthetic date index.

    Args:
        monthly_flows: 1D or 2D (n_years, 12) array of monthly flows.
        start_date: Start date for the index.
        freq: Pandas frequency string (default "MS" for month start).

    Returns:
        pd.Series with DatetimeIndex.
    """
    if monthly_flows.ndim == 2:
        monthly_flows = monthly_flows.flatten()
    index = pd.date_range(start=start_date, periods=len(monthly_flows), freq=freq)
    return pd.Series(monthly_flows, index=index, name="flow_cfs")


def compute_ssi(
    monthly_flows: np.ndarray,
    timescale: int = 3,
    dist: str = "gamma",
    reference_flows: Optional[np.ndarray] = None,
    start_date: str = "1950-10-01",
) -> Tuple[pd.Series, SynHydroSSI]:
    """Compute SSI values using SynHydro.

    If reference_flows is provided, the distribution is fit on the reference
    data and applied to monthly_flows (train/test split for synthetic traces).
    Otherwise, fit and transform on the same data.

    Args:
        monthly_flows: Monthly flow array (1D or 2D).
        timescale: SSI accumulation period (1, 3, 6, 12).
        dist: Distribution for fitting.
        reference_flows: Historical flows for fitting (optional).
        start_date: Start date for pd.Series index.

    Returns:
        Tuple of (ssi_series, fitted_calculator).
    """
    calc = make_ssi_calculator(timescale=timescale, dist=dist)

    if reference_flows is not None:
        ref_series = flows_to_series(reference_flows, start_date=start_date)
        calc.fit(ref_series)
        # Synthetic flows get a different start date to avoid index collision
        syn_series = flows_to_series(monthly_flows, start_date="2100-01-01")
        ssi = calc.transform(syn_series)
    else:
        series = flows_to_series(monthly_flows, start_date=start_date)
        ssi = calc.fit_transform(series)

    return ssi, calc


#: SSI value below which a month is counted as drought-stressed for
#: ``time_in_drought_fraction``. Matches SynHydro's "critical drought"
#: convention. A trace's time-in-drought fraction is independent of the
#: event-merging logic that aggregates consecutive months into events.
TIME_IN_DROUGHT_THRESHOLD: float = -1.0

#: Interim cap on first-event drought severity (|min SSI| within the event).
#: SSI values beyond ~4.5 in this pipeline are distribution-fitting tail
#: artifacts (small-sample gamma->normal transform), being fixed properly in
#: SynHydro. Until then we clip ``first_event_severity`` here so the artifact
#: cannot dominate the MOEA objective; the clip also flows into
#: ``first_event_onset_intensification_rate`` (= severity / rising_months),
#: bounding it at the same ceiling. Remove once the SynHydro fix lands.
FIRST_EVENT_SEVERITY_CLIP: float = 4.5


def compute_ssi_drought_characteristics(
    ssi_values: pd.Series,
    end_drought_threshold_months: int = 3,
    monthly_flows: Optional[np.ndarray] = None,
) -> Dict:
    """Compute aggregate drought characteristics from SSI series.

    Wraps SynHydro's get_drought_metrics to return a flat dict of aggregate
    statistics suitable for MOEA objectives.

    SynHydro drought identification rules:
      - Onset: SSI < 0 (mild), becomes critical when SSI <= -1
      - Termination: N consecutive months with SSI > 0 (default N=3)
      - Only critical droughts are recorded as events

    Args:
        ssi_values: SSI pd.Series from SynHydro.
        end_drought_threshold_months: Consecutive positive months to end
            a critical drought.
        monthly_flows: Optional 1D array of raw monthly flows. When
            provided, ``q10_flow_neg`` (negated 10th-percentile flow) is
            included so it can be picked up by the ``q10_flow`` metric.
            When ``None``, the key is absent and the metric extractor
            falls back to ``0.0``.

    Returns:
        Dict with the standard event-level metrics (frequency,
        mean_duration, mean_magnitude, mean_severity, mean_avg_severity,
        peak_severity_month, max_duration, max_magnitude, worst_severity,
        n_events) plus trace-level extras (time_in_drought_fraction, and
        q10_flow_neg when ``monthly_flows`` is supplied).
    """
    dm = get_drought_metrics(
        ssi_values,
        end_drought_threshold_months=end_drought_threshold_months,
    )

    valid_months = ssi_values.dropna()
    n_valid = len(valid_months)
    n_years = n_valid / 12.0 if n_valid > 0 else 0.0

    if n_valid > 0:
        time_in_drought_fraction = float(
            (valid_months <= TIME_IN_DROUGHT_THRESHOLD).sum() / n_valid
        )
    else:
        time_in_drought_fraction = 0.0

    base: Dict[str, float] = {
        "time_in_drought_fraction": time_in_drought_fraction,
    }
    if monthly_flows is not None:
        flows = np.asarray(monthly_flows, dtype=float).flatten()
        if flows.size > 0:
            # Negate so that "more drought" → larger objective, matching
            # the L1 device's sign convention.
            base["q10_flow_neg"] = float(-np.percentile(flows, 10.0))

    # The 10 Tier-A event stats are shared verbatim with
    # extended.compute_ssi_event_metrics; reuse the single implementation
    # (pass the already-computed dm so get_drought_metrics isn't re-run).
    from src.metrics.ssi_common import compute_ssi_tier_a
    core = compute_ssi_tier_a(ssi_values, dm=dm)

    if len(dm) == 0:
        return {
            **base,
            **core,
            "peak_severity_month": 0.0,
            "first_event_present": 0,
            "first_event_duration": 0.0,
            "first_event_severity": 0.0,
            "first_event_magnitude": 0.0,
            "first_event_start_month": 0.0,
            "first_event_peak_month": 0.0,
            "first_event_onset_intensification_rate": 0.0,
            "first_event_rising_limb_fraction": 0.0,
            "max_intensification_rate": 0.0,
        }

    # Cyclic mean of peak-severity calendar month across all events.
    # max_severity_date is a pandas Timestamp; extract .month (1..12) per
    # event and combine with cyclic_mean() to handle the Dec/Jan wraparound.
    if "max_severity_date" in dm.columns:
        try:
            peak_months = pd.to_datetime(dm["max_severity_date"]).dt.month.values
            peak_severity_month = cyclic_mean(peak_months, period=12.0)
        except Exception:
            peak_severity_month = 0.0
    else:
        peak_severity_month = 0.0

    # First-event characteristics (chronologically first critical SSI-3 drought).
    # Used by the first_event_* metric family at T=10y. The L1 device anchors
    # against the all-events historical maxes (worst_severity, max_magnitude,
    # max_duration) via DroughtMetric.max_partner.
    dm_sorted = dm.sort_values("start")
    first_row = dm_sorted.iloc[0]
    first_start_month = float(pd.to_datetime(first_row["start"]).month)
    if "max_severity_date" in dm_sorted.columns:
        first_peak_month = float(pd.to_datetime(first_row["max_severity_date"]).month)
    else:
        first_peak_month = 0.0

    # Shape descriptors. Two month-counts per event:
    #   peak_offset_months — calendar months from start to peak (rising span − 1)
    #   span_months        — calendar months from start to end (end − start + 1)
    # SynHydro's ``duration`` is the count of negative-SSI months only; it
    # excludes any brief positive gaps that didn't trigger termination.
    # ``span_months ≥ duration`` always, and using span as the shape
    # denominator keeps rising_limb_fraction bounded in (0, 1] (using
    # duration breaks the bound for events with mid-event positive gaps).
    def _peak_offset_months(row) -> int:
        if "max_severity_date" not in dm_sorted.columns:
            return 0
        start_ts = pd.to_datetime(row["start"])
        peak_ts = pd.to_datetime(row["max_severity_date"])
        return (peak_ts.year - start_ts.year) * 12 + (peak_ts.month - start_ts.month)

    def _span_months(row) -> int:
        start_ts = pd.to_datetime(row["start"])
        end_ts = pd.to_datetime(row["end"])
        return (end_ts.year - start_ts.year) * 12 + (end_ts.month - start_ts.month) + 1

    first_peak_offset = _peak_offset_months(first_row)
    first_rising_months = first_peak_offset + 1
    first_span = _span_months(first_row)
    first_duration = float(first_row["duration"])
    # Clip the SSI tail artifact (see FIRST_EVENT_SEVERITY_CLIP). Applied
    # before deriving the onset rate so both objectives share the ceiling.
    first_severity_abs = min(
        float(abs(first_row["severity"])), FIRST_EVENT_SEVERITY_CLIP
    )
    first_onset_rate = first_severity_abs / first_rising_months
    first_rising_fraction = (
        first_rising_months / first_span if first_span > 0 else 0.0
    )

    # Anti-ideal anchor for first_event_onset_intensification_rate: max
    # intensification rate across ALL events in the trace. Mirrors the
    # max_partner pattern used by first_event_severity / _magnitude / _duration.
    all_peak_offsets = dm_sorted.apply(_peak_offset_months, axis=1).to_numpy()
    all_rising_months = all_peak_offsets + 1
    all_severity_abs = dm_sorted["severity"].abs().to_numpy()
    all_rates = all_severity_abs / all_rising_months
    max_rate = float(all_rates.max()) if len(all_rates) > 0 else 0.0

    return {
        **base,
        **core,
        "peak_severity_month": float(peak_severity_month),
        "max_intensification_rate": max_rate,
        "first_event_present": 1,
        "first_event_duration": first_duration,
        "first_event_severity": first_severity_abs,
        "first_event_magnitude": float(abs(first_row["magnitude"])),
        "first_event_start_month": first_start_month,
        "first_event_peak_month": first_peak_month,
        "first_event_onset_intensification_rate": float(first_onset_rate),
        "first_event_rising_limb_fraction": float(first_rising_fraction),
    }


###############################################################################
# Threshold-based drought calculations (legacy, kept for comparison)
###############################################################################


def compute_drought_events(
    monthly_flows: np.ndarray,
    threshold: float,
) -> list:
    """Identify drought events from monthly flow time series.

    A drought event is a consecutive sequence of months where flow is
    below the threshold. Returns list of events, each a dict with
    start index, duration (months), and total deficit.

    Args:
        monthly_flows: 1D array of monthly flows.
        threshold: Flow threshold (e.g., P20 of historical).

    Returns:
        List of dicts with keys: 'start', 'duration', 'deficit', 'intensity'.
    """
    below = monthly_flows < threshold
    events = []
    in_event = False
    start = 0
    deficit = 0.0

    for i, is_below in enumerate(below):
        if is_below and not in_event:
            in_event = True
            start = i
            deficit = threshold - monthly_flows[i]
        elif is_below and in_event:
            deficit += threshold - monthly_flows[i]
        elif not is_below and in_event:
            duration = i - start
            events.append({
                "start": start,
                "duration": duration,
                "deficit": float(deficit),
                "intensity": float(deficit / duration),
            })
            in_event = False
            deficit = 0.0

    # Handle event at end of series
    if in_event:
        duration = len(monthly_flows) - start
        events.append({
            "start": start,
            "duration": duration,
            "deficit": float(deficit),
            "intensity": float(deficit / duration),
        })

    return events


def compute_drought_characteristics(
    monthly_flows: np.ndarray,
    threshold: float,
    n_years: float,
) -> dict:
    """Compute aggregate drought characteristics from monthly flows.

    Args:
        monthly_flows: 1D array of monthly flows.
        threshold: Drought threshold (e.g., P20 of historical).
        n_years: Number of years in the trace (for frequency calc).

    Returns:
        Dict with 'frequency' (events/decade), 'mean_duration' (months),
        'mean_intensity' (deficit/month), 'max_duration', 'max_deficit'.
    """
    events = compute_drought_events(monthly_flows, threshold)

    if len(events) == 0:
        return {
            "frequency": 0.0,
            "mean_duration": 0.0,
            "mean_intensity": 0.0,
            "max_duration": 0.0,
            "max_deficit": 0.0,
            "n_events": 0,
        }

    durations = [e["duration"] for e in events]
    intensities = [e["intensity"] for e in events]
    deficits = [e["deficit"] for e in events]

    return {
        "frequency": float(len(events) / n_years * 10),  # events per decade
        "mean_duration": float(np.mean(durations)),
        "mean_intensity": float(np.mean(intensities)),
        "max_duration": float(np.max(durations)),
        "max_deficit": float(np.max(deficits)),
        "n_events": len(events),
    }


#: Drought characteristic names whose natural metric is a calendar month
#: (1..12).
#:
#: Deprecated 2026-04-27: the cyclic-vs-non-cyclic distinction now lives
#: on :class:`src.metrics.drought_metrics.DroughtMetric` instances via the
#: ``is_cyclic`` flag and the ``CYCLIC_HEADROOM`` anti-ideal rule.
#: Retained here so that legacy code paths that still pass tuples of
#: string keys (rather than ``DroughtMetric`` objects) continue to
#: work; new code should consume metric instances and not consult this
#: set.
CYCLIC_MONTH_KEYS = frozenset({
    "peak_severity_month",
    "onset_month",
})


def drought_objectives(
    synthetic_chars: dict,
    anti_ideal: np.ndarray,
    objective_keys=("mean_duration", "mean_avg_severity"),
) -> np.ndarray:
    """MOEA-FIND L1 Device — the authoritative DD-11 formulation.

    See ``manuscript/design_decisions.md §DD-11`` for the full derivation.
    This function is the load-bearing piece of the MOEA-FIND method; do
    not introduce alternative formulations here.

    Produces the (K+1)-objective vector Borg minimises:

        f_j(x)     = D_j(x)                    for j = 1..K
        f_{K+1}(x) = ||D(x) - D*||_1           (Manhattan L1 distance)

    where ``D*`` is the anti-ideal point placed **outside the feasible
    region** by :func:`src.experiment.compute_ssi_anti_ideal`
    (``1.5 × max_hist`` for non-cyclic metrics, ``12 × headroom`` for
    cyclic-month metrics so that ``D_j <= D*_j`` for every feasible
    point).

    Under the ``D_j <= D*_j`` assumption the hyperplane identity
    ``sum_{j=1}^{K+1} f_j = sum_j D*_j`` is constant across the feasible
    set, so every feasible point is Pareto non-dominated: a decrease in
    any ``f_j`` (smaller ``D_j``) is exactly offset by an increase in
    ``f_{K+1}`` (farther from ``D*``). The feasible set maps injectively
    onto a hyperplane in (K+1)-objective space, and Borg's epsilon-box
    archive tiles that hyperplane — delivering interior coverage in
    drought-characteristic space.

    .. warning::
       Do **not** change ``f_j`` to ``|D_j - D*_j|``. That reformulation
       makes ``f_{K+1}`` linearly redundant with the first K objectives
       and collapses the Pareto front to the feasible point closest to
       ``D*``. DD-11 calls this reading "degenerate" and rejects it.

    Args:
        synthetic_chars: Drought characteristics dict from
            :func:`compute_ssi_drought_characteristics`.
        anti_ideal: Anti-ideal point ``D*`` (K-dim). Must satisfy
            ``D_j <= D*_j`` for every feasible ``D``; produced by
            :func:`src.metrics.drought_metrics.compute_anti_ideal` (or the
            legacy :func:`src.experiment.compute_ssi_anti_ideal`).
        objective_keys: Either a tuple of metric names (legacy) or a
            tuple of :class:`src.metrics.drought_metrics.DroughtMetric` instances
            (new). Order determines the objective vector layout.

    Returns:
        Array of length ``K+1``. The first ``K`` entries are the raw
        drought characteristics ``D_j``; the last entry is the Manhattan
        distance ``||D - D*||_1``.
    """
    # Local import to avoid a top-level cycle (drought_metrics may import
    # from objectives in the future).
    from src.metrics.drought_metrics import DroughtMetric

    target = np.asarray(anti_ideal, dtype=float)
    keys = tuple(objective_keys)
    k = len(keys)
    if target.shape != (k,):
        raise ValueError(
            f"anti_ideal must be length {k} to match objective_keys, "
            f"got shape {target.shape}"
        )

    metrics = np.zeros(k)
    if k > 0 and isinstance(keys[0], DroughtMetric):
        for i, m in enumerate(keys):
            metrics[i] = float(m.extract(synthetic_chars))
    else:
        for i, key in enumerate(keys):
            metrics[i] = float(synthetic_chars.get(key, 0.0))

    objs = np.empty(k + 1)
    objs[:k] = metrics
    objs[k] = manhattan_norm(metrics, target)
    return objs


def normalize_to_unit_cube(
    points: np.ndarray,
    lb: np.ndarray,
    ub: np.ndarray,
) -> np.ndarray:
    """Normalize points from [lb, ub] to [0, 1]^d.

    Args:
        points: Array of shape (n, d).
        lb: Lower bounds (d,).
        ub: Upper bounds (d,).

    Returns:
        Normalized array in [0, 1]^d.
    """
    return (points - lb) / (ub - lb)
