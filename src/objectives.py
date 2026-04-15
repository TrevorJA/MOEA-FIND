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


def cyclic_l1_distance(a: float, b: float, period: float = 12.0) -> float:
    """Wrapped L1 distance for cyclic (e.g. month-of-year) features.

    For a cyclic coordinate with period P (e.g. 12 months), the shortest
    L1 distance between two values is min(|a - b|, P - |a - b|). This
    preserves the hyperplane identity in the MOEA-FIND L1-simplex theorem
    (§3.2) because the objective is still a sum of nonnegative terms.

    Args:
        a, b: Cyclic coordinates in the same unit (e.g. both months).
        period: Period length. Default 12 for month-of-year.

    Returns:
        Wrapped distance in [0, period / 2].
    """
    d = abs(float(a) - float(b))
    return float(min(d, period - d))


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


def compute_ssi_drought_characteristics(
    ssi_values: pd.Series,
    end_drought_threshold_months: int = 3,
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

    Returns:
        Dict with frequency (events/decade), mean_duration, mean_magnitude,
        mean_severity (peak SSI), mean_avg_severity, max_duration,
        max_magnitude, worst_severity, n_events.
    """
    dm = get_drought_metrics(
        ssi_values,
        end_drought_threshold_months=end_drought_threshold_months,
    )

    valid_months = ssi_values.dropna()
    n_years = len(valid_months) / 12.0

    if len(dm) == 0:
        return {
            "frequency": 0.0,
            "mean_duration": 0.0,
            "mean_magnitude": 0.0,
            "mean_severity": 0.0,
            "mean_avg_severity": 0.0,
            "max_duration": 0.0,
            "max_magnitude": 0.0,
            "worst_severity": 0.0,
            "n_events": 0,
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

    return {
        "frequency": float(len(dm) / n_years * 10),
        "mean_duration": float(dm["duration"].mean()),
        "mean_magnitude": float(dm["magnitude"].abs().mean()),
        "mean_severity": float(dm["severity"].abs().mean()),
        "mean_avg_severity": float(dm["avg_severity"].abs().mean()),
        "peak_severity_month": float(peak_severity_month),
        "max_duration": float(dm["duration"].max()),
        "max_magnitude": float(dm["magnitude"].abs().max()),
        "worst_severity": float(dm["severity"].abs().max()),
        "n_events": len(dm),
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


#: Drought characteristic names whose natural metric is a cyclic calendar
#: month (1..12). The MOEA-FIND L1-simplex theorem supports cyclic axes
#: through the wrapped distance cyclic_l1_distance (§3.3 of the manuscript).
CYCLIC_MONTH_KEYS = frozenset({
    "peak_severity_month",
    "onset_month",
})


def drought_objectives(
    synthetic_chars: dict,
    anti_ideal: np.ndarray,
    objective_keys: tuple = ("mean_duration", "mean_avg_severity"),
    target: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Compute MOEA-FIND objectives from drought characteristics.

    Produces the (k+1)-objective vector used by Borg MOEA:

        f_j(x) = |D_j(x) - D*_j|              for j = 1..k
        f_{k+1}(x) = sum_j f_j(x)             (L1 Manhattan norm)

    For any objective whose name appears in :data:`CYCLIC_MONTH_KEYS`,
    the per-component distance |D_j - D*_j| is replaced by the wrapped
    cyclic L1 distance with period 12 months. This preserves the
    L1-simplex theorem because the objective is still a sum of
    nonnegative absolute-value terms.

    The function supports two target conventions for backwards
    compatibility:

    - If ``target`` is None, the historical anti-ideal pattern is used:
      f_j = D_j (i.e. the raw metric value). The simplex identity still
      holds with J_{k+1} = sum_j |D_j - 0|. This matches the pre-2026-04
      behavior of the MOEA-FIND Phase 1 experiments.
    - If ``target`` is supplied, it is treated as the user-specified
      drought characteristic target D* and f_j = d(D_j, D*_j), where d
      is the cyclic L1 distance for cyclic axes and the absolute-value
      distance otherwise. This is the convention described in §3 of
      the manuscript.

    Args:
        synthetic_chars: Drought characteristics of synthetic trace, e.g.
            from :func:`compute_ssi_drought_characteristics`.
        anti_ideal: Anti-ideal point in metric space; used for the
            target=None path and to set dimensionality. (k-dim.)
        objective_keys: Which characteristics to use as objectives.
        target: Optional drought characteristic target D* in the same
            order as ``objective_keys``. If provided, all f_j become
            cyclic-aware absolute deviations from D*_j.

    Returns:
        Objective vector (k+1 dimensional). Element k is the L1 sum of
        elements 0..k-1 (the Manhattan-norm auxiliary objective).
    """
    k = len(objective_keys)
    metrics = np.zeros(k)
    for i, key in enumerate(objective_keys):
        metrics[i] = synthetic_chars.get(key, 0.0)

    if target is None:
        # Legacy path: f_j = D_j, Manhattan norm from anti_ideal.
        objs = np.empty(k + 1)
        objs[:k] = metrics
        objs[k] = manhattan_norm(metrics, anti_ideal)
        return objs

    # Target-aware path with cyclic-aware per-component distance.
    target = np.asarray(target, dtype=float)
    assert target.shape == (k,), \
        f"target must be length {k}, got {target.shape}"

    deviations = np.zeros(k)
    for i, key in enumerate(objective_keys):
        if key in CYCLIC_MONTH_KEYS:
            deviations[i] = cyclic_l1_distance(metrics[i], target[i], period=12.0)
        else:
            deviations[i] = abs(metrics[i] - target[i])

    objs = np.empty(k + 1)
    objs[:k] = deviations
    objs[k] = float(np.sum(deviations))  # L1 simplex identity by construction
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
