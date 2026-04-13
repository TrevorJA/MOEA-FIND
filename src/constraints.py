"""Plausibility constraints for MOEA-FIND generated traces.

Constraints ensure that generated traces are physically plausible,
preventing Borg from finding degenerate traces that hit target
drought metrics through unrealistic flow patterns.
"""

import numpy as np


def autocorrelation_constraint(
    synthetic_monthly: np.ndarray,
    historical_monthly: np.ndarray,
    tolerance: float = 0.3,
) -> float:
    """Lag-1 autocorrelation constraint.

    Returns violation magnitude (0 = feasible, >0 = infeasible).

    Args:
        synthetic_monthly: 1D array of synthetic monthly flows.
        historical_monthly: 1D array of historical monthly flows.
        tolerance: Maximum allowed absolute difference in lag-1 autocorrelation.

    Returns:
        Constraint violation (0 if feasible).
    """
    def lag1_autocorr(x):
        if len(x) < 3:
            return 0.0
        return float(np.corrcoef(x[:-1], x[1:])[0, 1])

    rho_s = lag1_autocorr(synthetic_monthly)
    rho_h = lag1_autocorr(historical_monthly)
    violation = max(0.0, abs(rho_s - rho_h) - tolerance)
    return violation


def annual_stats_constraint(
    synthetic_monthly: np.ndarray,
    historical_monthly: np.ndarray,
    mean_tolerance: float = 0.5,
    cv_tolerance: float = 0.5,
) -> float:
    """Annual mean and CV constraint.

    Ensures the synthetic trace has plausible annual statistics
    (within tolerance fraction of historical).

    Args:
        synthetic_monthly: Synthetic flows reshaped or as 1D.
        historical_monthly: Historical flows as 1D.
        mean_tolerance: Fractional tolerance on annual mean ratio.
        cv_tolerance: Fractional tolerance on CV ratio.

    Returns:
        Sum of constraint violations (0 if feasible).
    """
    s_mean = np.mean(synthetic_monthly)
    h_mean = np.mean(historical_monthly)

    if h_mean > 0:
        mean_ratio = s_mean / h_mean
        mean_violation = max(0.0, abs(mean_ratio - 1.0) - mean_tolerance)
    else:
        mean_violation = 0.0

    s_cv = np.std(synthetic_monthly) / max(s_mean, 1e-6)
    h_cv = np.std(historical_monthly) / max(h_mean, 1e-6)

    if h_cv > 0:
        cv_ratio = s_cv / h_cv
        cv_violation = max(0.0, abs(cv_ratio - 1.0) - cv_tolerance)
    else:
        cv_violation = 0.0

    return mean_violation + cv_violation


def seasonal_cycle_constraint(
    synthetic_monthly_2d: np.ndarray,
    historical_monthly_2d: np.ndarray,
    tolerance: float = 0.5,
) -> float:
    """Seasonal cycle preservation constraint.

    Checks that monthly means of synthetic trace are within tolerance
    (as fraction) of historical monthly means.

    Args:
        synthetic_monthly_2d: Shape (n_years, 12).
        historical_monthly_2d: Shape (n_years, 12).
        tolerance: Fractional tolerance per month.

    Returns:
        Sum of violations across months (0 if feasible).
    """
    s_monthly_means = np.mean(synthetic_monthly_2d, axis=0)
    h_monthly_means = np.mean(historical_monthly_2d, axis=0)

    violation = 0.0
    for m in range(12):
        if h_monthly_means[m] > 0:
            ratio = s_monthly_means[m] / h_monthly_means[m]
            violation += max(0.0, abs(ratio - 1.0) - tolerance)

    return violation


def compute_all_constraints(
    synthetic_1d: np.ndarray,
    synthetic_2d: np.ndarray,
    historical_1d: np.ndarray,
    historical_2d: np.ndarray,
    autocorr_tol: float = 0.3,
    mean_tol: float = 0.5,
    cv_tol: float = 0.5,
    seasonal_tol: float = 0.5,
) -> list:
    """Compute all plausibility constraints.

    Args:
        synthetic_1d: Synthetic flows as 1D array.
        synthetic_2d: Synthetic flows as (n_years, 12).
        historical_1d: Historical flows as 1D.
        historical_2d: Historical flows as (n_years, 12).
        autocorr_tol: Tolerance for autocorrelation.
        mean_tol: Tolerance for annual mean.
        cv_tol: Tolerance for CV.
        seasonal_tol: Tolerance for seasonal cycle.

    Returns:
        List of constraint violation values (all 0 = feasible).
    """
    c1 = autocorrelation_constraint(
        synthetic_1d, historical_1d, autocorr_tol,
    )
    c2 = annual_stats_constraint(
        synthetic_1d, historical_1d, mean_tol, cv_tol,
    )
    c3 = seasonal_cycle_constraint(
        synthetic_2d, historical_2d, seasonal_tol,
    )
    return [c1, c2, c3]
