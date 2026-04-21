"""DV-space plausibility constraint for MOEA-FIND (ablation for manuscript SI).

Alternative to the five hydrologic-space constraints in :mod:`src.constraints`.
Constrains the flattened Kirsch decision-variable vector to stay close to the
``U[0, 1]`` distribution from which it would be drawn under random (non-MOEA)
Kirsch bootstrap.

Two scalar statistics are available, both bounded in ``[0, 1]``:

    - ``"l2_star"``: L2-star discrepancy of the flattened DVs treated as a
      point set in ``[0, 1]`` (shape ``(N, 1)``). Larger values indicate DVs
      cluster away from uniform coverage — severe-drought traces concentrate
      DVs near 0 which inflates this statistic.
    - ``"ks"``: one-sample Kolmogorov-Smirnov statistic of the flattened DVs
      against ``U[0, 1]``.

Both are bootstrap-calibrated against random ``U[0, 1]`` draws of the same
length as the DV vector by
``workflows/diagnostics/diag_dv_uniformity_calibration.py``, which produces
``calibrated_dv_tolerances.json``.

Unlike the hydrologic calibration, the DV calibration only bootstraps the
Kirsch branch: random Kirsch DV draws are analytically ``U[0, 1]^{N}``, and
there is no "historical DV" distribution to compare against (historical
observations don't have DVs).

This module exists to *falsify* the hypothesis that DV-uniformity alone
suffices to enforce hydrologic plausibility during drought search. It is
*expected* to yield a collapsed Pareto front in the ablation study, because
severe droughts require non-uniform DVs. See the plan file for full rationale.

``Violation``, ``_two_tier``, and ``ConstraintResult`` are imported from
:mod:`src.constraints` so that downstream code (``run_experiment``, diagnostic
serialization) can treat DV and hydrologic constraints through the same
interface.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Union

import numpy as np
from scipy import stats
from scipy.stats import qmc

from src.constraints import ConstraintResult, Violation, _two_tier


VALID_STATISTICS = ("l2_star", "ks", "ad")


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class DVUniformityConfig:
    """Calibrated tolerance for the DV-uniformity constraint.

    Attributes:
        T_years: Target synthetic trace length in years (metadata).
        n_dvs: Length of the flattened DV vector the tolerance was calibrated
            for. :func:`compute_dv_constraint` checks this matches the DV
            vector it receives.
        statistic: Which scalar uniformity statistic to use
            (``"l2_star"`` or ``"ks"``).
        tolerance: Half-width of the 95% envelope of the chosen statistic
            under random ``U[0, 1]`` DV draws of length ``n_dvs``. The
            statistic is nonnegative, so this is the upper edge of the
            soft band rather than a symmetric half-width.
        hard_multiplier: Hard-infeasibility activates beyond
            ``hard_multiplier * tolerance``. Default 2.0, matching the
            hydrologic default.
        soft_weight: Weight applied to the soft penalty when folded into
            the Manhattan-norm objective. Default 0.5 — one-half the
            hydrologic default because there is only one constraint here
            versus five in the hydrologic arm, and the combined weighted
            penalty should stay in the same order of magnitude.
        enabled: Single on/off flag; mirrors ``ConstraintConfig.enabled``
            API at lower cardinality.
    """

    T_years: int
    n_dvs: int
    statistic: str
    tolerance: float
    hard_multiplier: float = 2.0
    soft_weight: float = 0.5
    enabled: bool = True

    def __post_init__(self) -> None:
        if self.statistic not in VALID_STATISTICS:
            raise ValueError(
                f"statistic must be one of {VALID_STATISTICS}, "
                f"got {self.statistic!r}"
            )
        if self.tolerance < 0:
            raise ValueError(
                f"tolerance must be nonnegative, got {self.tolerance!r}"
            )
        if self.n_dvs <= 0:
            raise ValueError(
                f"n_dvs must be positive, got {self.n_dvs!r}"
            )

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "DVUniformityConfig":
        return cls(
            T_years=int(d["T_years"]),
            n_dvs=int(d["n_dvs"]),
            statistic=str(d["statistic"]),
            tolerance=float(d["tolerance"]),
            hard_multiplier=float(d.get("hard_multiplier", 2.0)),
            soft_weight=float(d.get("soft_weight", 0.5)),
            enabled=bool(d.get("enabled", True)),
        )

    @classmethod
    def from_calibration_json(
        cls,
        path: Union[str, Path],
        site_label: str = "cannonsville",
        T_years: int = 20,
        statistic: str = "l2_star",
    ) -> "DVUniformityConfig":
        """Build a config from ``diag_dv_uniformity_calibration.py`` output.

        The calibration JSON is keyed as
        ``{site_label}_T{T_years}_{statistic}`` so the two statistic
        variants (l2_star, ks) can coexist in one file.
        """
        payload = json.loads(Path(path).read_text())
        key = f"{site_label}_T{T_years}_{statistic}"
        if key not in payload:
            raise KeyError(
                f"No entry {key!r} in {path}. "
                f"Available: {list(payload.keys())}"
            )
        entry = payload[key]
        return cls(
            T_years=int(entry.get("T_years", T_years)),
            n_dvs=int(entry["n_dvs"]),
            statistic=statistic,
            tolerance=float(entry["tolerance"]),
        )


# ---------------------------------------------------------------------------
# Statistic functions
# ---------------------------------------------------------------------------

def dv_l2_star_statistic(dvs: np.ndarray) -> float:
    """L2-star discrepancy of the flattened DV vector as a 1D point set.

    Returns the star-discrepancy of the DVs viewed as points in ``[0, 1]^1``
    via ``scipy.stats.qmc.discrepancy``. Values near 0 indicate near-uniform
    coverage; larger values indicate clustering.
    """
    x = np.asarray(dvs, dtype=float).reshape(-1)
    if x.size == 0:
        return 0.0
    x = np.clip(x, 0.0, 1.0)
    return float(qmc.discrepancy(x.reshape(-1, 1), method="L2-star"))


def dv_ks_statistic(dvs: np.ndarray) -> float:
    """One-sample KS statistic of the flattened DVs vs ``U[0, 1]``.

    Range: ``[0, 1]``. Equals 0 when the empirical CDF matches the uniform
    CDF exactly at every order-statistic. Measures the maximum gap between
    empirical and theoretical CDFs — relatively insensitive to local tail
    clustering because tails contribute to a single gap rather than
    multiple ones.
    """
    x = np.asarray(dvs, dtype=float).reshape(-1)
    if x.size == 0:
        return 0.0
    x = np.clip(x, 0.0, 1.0)
    return float(stats.kstest(x, "uniform").statistic)


def dv_anderson_darling_statistic(dvs: np.ndarray) -> float:
    """Anderson-Darling ``A²`` statistic of flattened DVs vs ``U[0, 1]``.

    The AD weighting ``1 / [F(x) · (1 - F(x))]`` diverges at the
    boundaries, so this statistic is explicitly *tail-weighted*: small
    numbers of DVs clustered near 0 or 1 inflate A² more than comparable
    clustering in the middle of [0, 1]. Complementary to L2-star
    (integrates squared discrepancy uniformly) and KS (max single gap).

    Closed form for ``F = U[0, 1]``:

        A² = -n - (1/n) · Σ_{i=1..n} (2i − 1) · [ln(x_{(i)}) + ln(1 − x_{(n+1−i)})]

    where ``x_{(1)} ≤ ... ≤ x_{(n)}`` are the sorted DVs. Under the null
    (samples truly ``U[0, 1]``), ``A²`` is ≈ 1 on average with a 95th
    percentile around 2.5 — bootstrap calibration picks up the
    finite-sample envelope.
    """
    x = np.asarray(dvs, dtype=float).reshape(-1)
    n = x.size
    if n == 0:
        return 0.0
    # Clip inward to avoid log(0) / log(negative) at the sorted extremes.
    eps = 1e-12
    x_sorted = np.clip(np.sort(x), eps, 1.0 - eps)
    i = np.arange(1, n + 1, dtype=float)
    log_terms = np.log(x_sorted) + np.log(1.0 - x_sorted[::-1])
    a2 = -n - np.sum((2.0 * i - 1.0) * log_terms) / n
    return float(a2)


def statistic_fn(name: str):
    """Return the statistic function for ``name``; raise on unknown names."""
    if name == "l2_star":
        return dv_l2_star_statistic
    if name == "ks":
        return dv_ks_statistic
    if name == "ad":
        return dv_anderson_darling_statistic
    raise ValueError(
        f"Unknown DV-uniformity statistic {name!r}; "
        f"expected one of {VALID_STATISTICS}"
    )


# ---------------------------------------------------------------------------
# Constraint evaluation
# ---------------------------------------------------------------------------

def compute_dv_constraint(
    dvs: np.ndarray,
    cfg: DVUniformityConfig,
) -> ConstraintResult:
    """Evaluate the DV-uniformity constraint on one decision-variable vector.

    Returns a :class:`src.constraints.ConstraintResult` with a single-element
    ``hard_violations`` list, matching the shape expected by
    :func:`src.experiment_utils.run_experiment` when ``n_constrs = 1``.
    If ``cfg.enabled`` is False, returns a trivially feasible result.
    """
    name = f"dv_uniformity_{cfg.statistic}"
    hard_cut = cfg.hard_multiplier * cfg.tolerance

    if not cfg.enabled:
        v = Violation(name, 0.0, cfg.tolerance, hard_cut, 0.0, 0.0)
        return ConstraintResult(
            hard_violations=[0.0],
            soft_penalty_total=0.0,
            soft_penalty_weighted=0.0,
            violations=[v],
            feasible=True,
        )

    x = np.asarray(dvs, dtype=float).reshape(-1)
    if len(x) != cfg.n_dvs:
        raise ValueError(
            f"DV length mismatch: calibration used n_dvs={cfg.n_dvs} "
            f"but compute_dv_constraint received {len(x)}."
        )

    stat_value = statistic_fn(cfg.statistic)(x)
    v = _two_tier(name, stat_value, cfg.tolerance, cfg.hard_multiplier)
    soft = float(v.soft_penalty)
    return ConstraintResult(
        hard_violations=[v.hard_violation],
        soft_penalty_total=soft,
        soft_penalty_weighted=float(cfg.soft_weight * soft),
        violations=[v],
        feasible=(v.hard_violation <= 0.0),
    )
