"""Plausibility constraints for MOEA-FIND generated traces.

Constraints ensure that synthetic traces produced by Borg are hydrologically
plausible, preventing the optimizer from finding degenerate traces that reach
target drought metrics through unrealistic flow patterns. The constraint set
is deliberately simple (no Hurst exponent, no lag-2 AC) and relies on
**bootstrap-calibrated tolerances** rather than hand-picked ones: tolerances
are loaded from ``ConstraintConfig``, which in turn is produced by
``workflows/02_calibration/constraint_calibration.py``.

Every constraint function returns a :class:`Violation` that carries both a
hard-infeasibility magnitude and a soft quadratic penalty:

    - hard_violation > 0 when the sample statistic moves outside
      ``2 * calibrated_tolerance`` from the historical reference. This is
      what Borg's constraint-domination rule ranks by.
    - soft_penalty is quadratic in the fractional excess inside the
      1×–2× tolerance band and zero elsewhere. The Manhattan-norm
      auxiliary objective receives it with a configurable weight so the
      optimizer has gradient guidance away from the constraint boundary.

Statistics covered:

    - annual_mean         (fractional deviation from historical annual mean)
    - annual_cv           (fractional deviation from historical annual CV)
    - lag1_ac_monthly     (absolute deviation from historical lag-1 AC)
    - non_drought_mean    (fractional deviation from historical non-drought mean)
    - seasonal_cycle      (max fractional deviation across calendar months)

``Hurst exponent``, lag-2 AC, skewness, and annual quantile constraints are
intentionally out of scope (see plan file and constraints_spec.md scoping note).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class ConstraintConfig:
    """Calibrated constraint tolerances and precomputed historical references.

    Populated from ``calibrated_tolerances.json`` written by
    ``workflows/02_calibration/constraint_calibration.py``. All ``*_tol`` fields are the
    calibrated *half-width* of the 95% bootstrap CI of the corresponding
    statistic at the target trace length. The hard-infeasibility cut-off is
    ``hard_multiplier * tol`` for every statistic.

    Attributes:
        T_years: Target synthetic trace length in years.
        hist_annual_mean: Historical annual-total mean flow.
        hist_annual_cv: Historical annual-total CV.
        hist_lag1_ac_monthly: Historical lag-1 AC of flattened monthly series.
        hist_non_drought_mean: Mean monthly flow over non-drought months
            (SSI > 0) of the full historical record.
        hist_monthly_means: Per-calendar-month mean flows (length 12).
        annual_mean_tol: Fractional-deviation tolerance on annual mean.
        annual_cv_tol: Fractional-deviation tolerance on annual CV.
        lag1_ac_tol: Absolute-deviation tolerance on lag-1 AC.
        non_drought_mean_tol: Fractional-deviation tolerance on non-drought mean.
        seasonal_cycle_tol: Tolerance on max per-month fractional deviation.
        hard_multiplier: Hard-infeasibility activates beyond
            ``hard_multiplier * tol``. Default 2.0.
        soft_weight: Weight applied to the aggregated soft penalty when
            folded into the Manhattan-norm objective.
        enabled: Per-statistic on/off flags.
    """

    T_years: int
    hist_annual_mean: float
    hist_annual_cv: float
    hist_lag1_ac_monthly: float
    hist_non_drought_mean: float
    hist_monthly_means: List[float]

    annual_mean_tol: float
    annual_cv_tol: float
    lag1_ac_tol: float
    non_drought_mean_tol: float
    seasonal_cycle_tol: float

    hard_multiplier: float = 2.0
    soft_weight: float = 1.0
    enabled: Dict[str, bool] = field(default_factory=lambda: {
        "annual_mean": True,
        "annual_cv": True,
        "lag1_ac_monthly": True,
        "non_drought_mean": True,
        "seasonal_cycle": True,
    })

    def to_dict(self) -> Dict:
        d = asdict(self)
        d["hist_monthly_means"] = list(self.hist_monthly_means)
        return d

    def to_json(self, path: Union[str, Path]) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2))

    @classmethod
    def from_dict(cls, d: Dict) -> "ConstraintConfig":
        return cls(
            T_years=int(d["T_years"]),
            hist_annual_mean=float(d["hist_annual_mean"]),
            hist_annual_cv=float(d["hist_annual_cv"]),
            hist_lag1_ac_monthly=float(d["hist_lag1_ac_monthly"]),
            hist_non_drought_mean=float(d["hist_non_drought_mean"]),
            hist_monthly_means=[float(x) for x in d["hist_monthly_means"]],
            annual_mean_tol=float(d["annual_mean_tol"]),
            annual_cv_tol=float(d["annual_cv_tol"]),
            lag1_ac_tol=float(d["lag1_ac_tol"]),
            non_drought_mean_tol=float(d["non_drought_mean_tol"]),
            seasonal_cycle_tol=float(d["seasonal_cycle_tol"]),
            hard_multiplier=float(d.get("hard_multiplier", 2.0)),
            soft_weight=float(d.get("soft_weight", 1.0)),
            enabled=dict(d.get("enabled", {
                "annual_mean": True, "annual_cv": True,
                "lag1_ac_monthly": True, "non_drought_mean": True,
                "seasonal_cycle": True,
            })),
        )

    @classmethod
    def from_calibration_json(
        cls,
        path: Union[str, Path],
        site_label: str = "cannonsville",
        T_years: int = 20,
    ) -> "ConstraintConfig":
        """Build a ConstraintConfig from the calibration script output.

        The calibration JSON is keyed as ``{site_label}_T{T_years}``.
        """
        payload = json.loads(Path(path).read_text())
        key = f"{site_label}_T{T_years}"
        if key not in payload:
            raise KeyError(
                f"No entry {key!r} in {path}. Available: {list(payload.keys())}"
            )
        entry = payload[key]
        hist = entry["historical"]
        tol = entry["tolerances"]
        return cls(
            T_years=int(entry.get("T_years", T_years)),
            hist_annual_mean=float(hist["annual_mean"]),
            hist_annual_cv=float(hist["annual_cv"]),
            hist_lag1_ac_monthly=float(hist["lag1_ac_monthly"]),
            hist_non_drought_mean=float(hist["non_drought_mean"]),
            hist_monthly_means=[float(x) for x in hist["monthly_means"]],
            annual_mean_tol=float(tol["annual_mean"]),
            annual_cv_tol=float(tol["annual_cv"]),
            lag1_ac_tol=float(tol["lag1_ac_monthly"]),
            non_drought_mean_tol=float(tol["non_drought_mean"]),
            seasonal_cycle_tol=float(tol["seasonal_cycle_max_dev"]),
        )


# ---------------------------------------------------------------------------
# Violation primitive
# ---------------------------------------------------------------------------

@dataclass
class Violation:
    """Single-constraint evaluation result.

    hard_violation: 0 if the sample is inside the hard band
        (``|deviation| <= hard_multiplier * tol``). Otherwise, the excess
        magnitude beyond the hard cut-off. This is what Borg's constraint
        vector receives.
    soft_penalty: 0 outside the soft band. Inside ``[tol, hard_mult*tol]``
        it is a quadratic penalty normalised so that it reaches 1.0 at the
        hard cut-off. Folded into the Manhattan-norm objective with a weight.
    deviation: signed deviation value actually computed, for diagnostics.
    """
    name: str
    deviation: float
    tolerance: float
    hard_cutoff: float
    hard_violation: float
    soft_penalty: float


def _two_tier(
    name: str,
    deviation: float,
    tolerance: float,
    hard_multiplier: float,
) -> Violation:
    """Return a :class:`Violation` for a one-sided deviation magnitude.

    deviation must already be a nonnegative scalar (e.g. ``abs(...)``).
    Soft penalty is quadratic in the fractional progress through the soft
    band, so it is 0 at ``|deviation| == tolerance`` and 1 at
    ``|deviation| == hard_multiplier * tolerance``.
    """
    dev = float(abs(deviation))
    tol = float(tolerance)
    hard_cut = float(hard_multiplier) * tol

    if tol <= 0 or not np.isfinite(dev):
        return Violation(name, dev, tol, hard_cut, 0.0, 0.0)

    if dev <= tol:
        return Violation(name, dev, tol, hard_cut, 0.0, 0.0)

    hard_violation = max(0.0, dev - hard_cut)

    if hard_cut > tol:
        frac = (min(dev, hard_cut) - tol) / (hard_cut - tol)
    else:
        frac = 1.0 if dev > tol else 0.0
    soft_penalty = float(frac * frac)

    return Violation(name, dev, tol, hard_cut, float(hard_violation), soft_penalty)


# ---------------------------------------------------------------------------
# Statistic helpers
# ---------------------------------------------------------------------------

def _annual_totals(monthly_1d: np.ndarray) -> np.ndarray:
    x = np.asarray(monthly_1d, dtype=float)
    n = (len(x) // 12) * 12
    if n == 0:
        return np.array([0.0])
    return x[:n].reshape(-1, 12).sum(axis=1)


def _lag1_ac(monthly_1d: np.ndarray) -> float:
    x = np.asarray(monthly_1d, dtype=float)
    if len(x) < 3:
        return 0.0
    return float(np.corrcoef(x[:-1], x[1:])[0, 1])


def _monthly_cycle(monthly_2d: np.ndarray) -> np.ndarray:
    return np.asarray(monthly_2d, dtype=float).mean(axis=0)


# ---------------------------------------------------------------------------
# Individual constraint functions — each returns a Violation
# ---------------------------------------------------------------------------

def annual_mean_violation(
    synthetic_1d: np.ndarray,
    cfg: ConstraintConfig,
) -> Violation:
    """Fractional deviation of annual mean flow from historical."""
    if cfg.hist_annual_mean <= 0:
        return Violation("annual_mean", 0.0, cfg.annual_mean_tol,
                         cfg.hard_multiplier * cfg.annual_mean_tol, 0.0, 0.0)
    syn_mean = float(_annual_totals(synthetic_1d).mean())
    dev = syn_mean / cfg.hist_annual_mean - 1.0
    return _two_tier("annual_mean", dev, cfg.annual_mean_tol, cfg.hard_multiplier)


def annual_cv_violation(
    synthetic_1d: np.ndarray,
    cfg: ConstraintConfig,
) -> Violation:
    """Fractional deviation of annual CV from historical."""
    if cfg.hist_annual_cv <= 0:
        return Violation("annual_cv", 0.0, cfg.annual_cv_tol,
                         cfg.hard_multiplier * cfg.annual_cv_tol, 0.0, 0.0)
    totals = _annual_totals(synthetic_1d)
    mean = totals.mean()
    if mean <= 0:
        return _two_tier("annual_cv", 1.0, cfg.annual_cv_tol, cfg.hard_multiplier)
    syn_cv = float(totals.std(ddof=1) / mean)
    dev = syn_cv / cfg.hist_annual_cv - 1.0
    return _two_tier("annual_cv", dev, cfg.annual_cv_tol, cfg.hard_multiplier)


def lag1_ac_violation(
    synthetic_1d: np.ndarray,
    cfg: ConstraintConfig,
) -> Violation:
    """Absolute deviation of lag-1 AC from historical."""
    syn_rho = _lag1_ac(synthetic_1d)
    dev = syn_rho - cfg.hist_lag1_ac_monthly
    return _two_tier("lag1_ac_monthly", dev, cfg.lag1_ac_tol, cfg.hard_multiplier)


def non_drought_mean_violation(
    synthetic_1d: np.ndarray,
    cfg: ConstraintConfig,
    ssi_calc,
    start_date: str = "2100-01-01",
) -> Violation:
    """Fractional deviation of non-drought mean (SSI > 0 months).

    Requires a pre-fitted SynHydro SSI calculator; the classification of a
    month as "drought" vs "non-drought" is made by transforming the synthetic
    series through the calibrated distribution and thresholding at SSI > 0.
    """
    if cfg.hist_non_drought_mean <= 0:
        return Violation("non_drought_mean", 0.0, cfg.non_drought_mean_tol,
                         cfg.hard_multiplier * cfg.non_drought_mean_tol,
                         0.0, 0.0)

    from src.objectives import flows_to_series

    series = flows_to_series(synthetic_1d, start_date=start_date)
    ssi = ssi_calc.transform(series)
    ssi_arr = np.asarray(ssi.values, dtype=float)
    # Align flows to SSI output (SSI drops leading months for accumulation > 1)
    flows_arr = np.asarray(series.loc[ssi.index].values, dtype=float)
    mask = np.isfinite(ssi_arr) & (ssi_arr > 0)
    if not mask.any():
        return _two_tier("non_drought_mean", 1.0, cfg.non_drought_mean_tol,
                         cfg.hard_multiplier)
    syn_mean = float(flows_arr[mask].mean())
    dev = syn_mean / cfg.hist_non_drought_mean - 1.0
    return _two_tier(
        "non_drought_mean", dev, cfg.non_drought_mean_tol, cfg.hard_multiplier
    )


def seasonal_cycle_violation(
    synthetic_2d: np.ndarray,
    cfg: ConstraintConfig,
) -> Violation:
    """Max over calendar months of |syn_month_mean / hist_month_mean - 1|."""
    syn_cycle = _monthly_cycle(synthetic_2d)
    hist = np.asarray(cfg.hist_monthly_means, dtype=float)
    safe_hist = np.where(hist > 0, hist, 1.0)
    rel = np.abs(syn_cycle / safe_hist - 1.0)
    max_dev = float(np.max(rel))
    return _two_tier(
        "seasonal_cycle", max_dev, cfg.seasonal_cycle_tol, cfg.hard_multiplier
    )


# ---------------------------------------------------------------------------
# Aggregated result
# ---------------------------------------------------------------------------

@dataclass
class ConstraintResult:
    """Aggregated constraint evaluation for a single trace.

    hard_violations is the vector consumed by Borg's constraint-domination
    rule: one entry per enabled constraint. A trace is feasible iff every
    entry is 0. soft_penalty_total is added (with ``cfg.soft_weight``) to the
    Manhattan-norm auxiliary objective so the optimizer gets a gradient away
    from the constraint boundary instead of a cliff.
    """
    hard_violations: List[float]
    soft_penalty_total: float
    soft_penalty_weighted: float
    violations: List[Violation]
    feasible: bool

    def to_dict(self) -> Dict:
        return {
            "hard_violations": list(self.hard_violations),
            "soft_penalty_total": self.soft_penalty_total,
            "soft_penalty_weighted": self.soft_penalty_weighted,
            "feasible": self.feasible,
            "violations": [asdict(v) for v in self.violations],
        }


def compute_all_constraints(
    synthetic_1d: np.ndarray,
    synthetic_2d: np.ndarray,
    cfg: ConstraintConfig,
    ssi_calc=None,
) -> ConstraintResult:
    """Evaluate every enabled plausibility constraint on one synthetic trace.

    Args:
        synthetic_1d: Flattened monthly flows (length T*12).
        synthetic_2d: Same flows reshaped as (T, 12).
        cfg: ConstraintConfig with calibrated tolerances.
        ssi_calc: Pre-fitted SynHydro SSI calculator, required if the
            non-drought-mean constraint is enabled.

    Returns:
        ConstraintResult with the per-constraint hard-violation vector,
        total soft penalty, and a list of individual :class:`Violation`
        objects for diagnostics.
    """
    violations: List[Violation] = []

    if cfg.enabled.get("annual_mean", True):
        violations.append(annual_mean_violation(synthetic_1d, cfg))

    if cfg.enabled.get("annual_cv", True):
        violations.append(annual_cv_violation(synthetic_1d, cfg))

    if cfg.enabled.get("lag1_ac_monthly", True):
        violations.append(lag1_ac_violation(synthetic_1d, cfg))

    if cfg.enabled.get("non_drought_mean", True):
        if ssi_calc is None:
            raise ValueError(
                "non_drought_mean constraint is enabled but no ssi_calc was "
                "supplied to compute_all_constraints()."
            )
        violations.append(
            non_drought_mean_violation(synthetic_1d, cfg, ssi_calc)
        )

    if cfg.enabled.get("seasonal_cycle", True):
        violations.append(seasonal_cycle_violation(synthetic_2d, cfg))

    hard = [v.hard_violation for v in violations]
    soft = float(sum(v.soft_penalty for v in violations))
    return ConstraintResult(
        hard_violations=hard,
        soft_penalty_total=soft,
        soft_penalty_weighted=float(cfg.soft_weight * soft),
        violations=violations,
        feasible=all(h <= 0.0 for h in hard),
    )
