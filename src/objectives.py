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

    return {
        "frequency": float(len(dm) / n_years * 10),
        "mean_duration": float(dm["duration"].mean()),
        "mean_magnitude": float(dm["magnitude"].abs().mean()),
        "mean_severity": float(dm["severity"].abs().mean()),
        "mean_avg_severity": float(dm["avg_severity"].abs().mean()),
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


def drought_objectives(
    synthetic_chars: dict,
    anti_ideal: np.ndarray,
    objective_keys: tuple = ("mean_duration", "mean_avg_severity"),
) -> np.ndarray:
    """Compute MOEA-FIND objectives from drought characteristics.

    Objectives are the actual drought metrics (not ratios). The Pareto
    front maps directly to physically interpretable drought space:
    duration in months, intensity in cfs deficit, etc.

    J_i = D_i (actual metric value)
    J_{k+1} = ||D - D*||_1 (Manhattan norm from anti-ideal)

    All objectives are MINIMIZED. The anti-ideal is placed at the
    maximum plausible values of each metric, so the Manhattan norm
    forces solutions toward the anti-ideal (more severe droughts)
    while individual objectives push toward mild droughts.

    Args:
        synthetic_chars: Drought characteristics of synthetic trace.
        anti_ideal: Anti-ideal point in actual metric space (k-dim).
            E.g., [max_duration_months, max_intensity_cfs].
        objective_keys: Which characteristics to use as objectives.

    Returns:
        Objective vector (k+1 dimensional).
    """
    k = len(objective_keys)
    metrics = np.zeros(k)
    for i, key in enumerate(objective_keys):
        metrics[i] = synthetic_chars.get(key, 0.0)

    objs = np.empty(k + 1)
    objs[:k] = metrics
    objs[k] = manhattan_norm(metrics, anti_ideal)
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
