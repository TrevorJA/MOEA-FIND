"""Loaders for calibrated constraint configurations.

Wraps :class:`src.constraints.ConstraintConfig` and
:class:`src.constraints_dv.DVUniformityConfig` with the print-and-fall-back
behavior previously inlined in ``04_kirsch_single_site.py``: a missing or
malformed JSON downgrades to ``None`` (i.e. unconstrained) with a warning,
rather than raising. This keeps the driver scripts free of constraint
plumbing.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from src.constraints import ConstraintConfig
from src.constraints_dv import DVUniformityConfig


def load_hydrologic_constraints(
    json_path: Optional[Path],
    site_label: str,
    T_years: int,
) -> Optional[ConstraintConfig]:
    """Load a calibrated five-statistic ConstraintConfig.

    Returns ``None`` (and prints a warning) if the JSON is missing or fails
    to parse — callers treat this as "run unconstrained."

    Args:
        json_path: Path to a calibrated_tolerances.json. If ``None`` is
            passed, returns ``None`` silently.
        site_label: Site key inside the calibration JSON (e.g.
            ``"cannonsville"``).
        T_years: Synthetic trace length in years; calibration tolerances
            are keyed by ``(site_label, T_years)``.

    Returns:
        Loaded :class:`ConstraintConfig`, or ``None`` if unavailable.
    """
    if json_path is None:
        return None
    json_path = Path(json_path)
    if not json_path.exists():
        print(f"[constraints] WARNING: hydrologic JSON not found at "
              f"{json_path}; running unconstrained.")
        return None
    try:
        cfg = ConstraintConfig.from_calibration_json(
            json_path, site_label=site_label, T_years=T_years,
        )
    except Exception as exc:
        print(f"[constraints] WARNING: failed to load hydrologic JSON: "
              f"{exc}; running unconstrained.")
        return None
    print(f"[constraints] Loaded hydrologic constraints from {json_path}")
    print(f"   annual_mean_tol={cfg.annual_mean_tol:.3f}, "
          f"annual_cv_tol={cfg.annual_cv_tol:.3f}, "
          f"lag1_ac_tol={cfg.lag1_ac_tol:.3f}, "
          f"non_drought_mean_tol={cfg.non_drought_mean_tol:.3f}, "
          f"seasonal_cycle_tol={cfg.seasonal_cycle_tol:.3f}")
    return cfg


def load_dv_uniformity_constraints(
    json_path: Optional[Path],
    site_label: str,
    T_years: int,
    statistic: str,
) -> Optional[DVUniformityConfig]:
    """Load a calibrated DVUniformityConfig.

    Returns ``None`` (and prints a warning) if the JSON is unset, missing,
    or fails to parse.

    Args:
        json_path: Path to a calibrated_dv_tolerances.json.
        site_label: Site key inside the calibration JSON.
        T_years: Synthetic trace length in years.
        statistic: Which statistic the calibration entry to use
            (``"l2_star"``, ``"ks"``, or ``"ad"``).

    Returns:
        Loaded :class:`DVUniformityConfig`, or ``None`` if unavailable.
    """
    if json_path is None:
        print("[constraints] WARNING: dv_uniform constraint requested but "
              "no JSON path supplied; running unconstrained.")
        return None
    json_path = Path(json_path)
    if not json_path.exists():
        print(f"[constraints] WARNING: DV-uniformity JSON not found at "
              f"{json_path}; running unconstrained.")
        return None
    cfg = DVUniformityConfig.from_calibration_json(
        json_path, site_label=site_label, T_years=T_years, statistic=statistic,
    )
    print(f"[constraints] Loaded DV-uniformity constraint from {json_path}")
    print(f"   statistic={cfg.statistic}, tol={cfg.tolerance:.4g}, "
          f"n_dvs={cfg.n_dvs}")
    return cfg
