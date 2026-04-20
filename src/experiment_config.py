"""Centralized MOEA-FIND experiment configuration.

Single source of truth for the knobs that define a MOEA-FIND run:
objective set, anti-ideal placement, epsilons, NFE, trace length, DV
injection mode, algorithm, and constraint source. Experiment driver
scripts (e.g. ``workflows/experiments/04_kirsch_single_site.py``) read
their defaults from :data:`DEFAULT_EXPERIMENT`, which can be either (a)
used as-is, (b) overridden by CLI args, or (c) replaced with a
site-specific :class:`ExperimentConfig` instance built at call time.

All runtime data (historical flows, SSI calculators, generator
instances) stay out of this module; it is purely declarative config.

Changing a default here affects every experiment script that reads from
``DEFAULT_EXPERIMENT``. CLI args that specify the same knob override the
default for that particular invocation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple

# Module-level project root so config paths resolve consistently.
PROJECT_ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# Per-objective epsilon defaults
# ---------------------------------------------------------------------------
#
# Used to build the epsilon list for epsilon-dominance MOEAs. The
# trailing epsilon for the Manhattan-norm auxiliary objective is the sum
# of per-objective epsilons (see :func:`src.experiment_utils.build_epsilons`).
#
# Keys not listed here fall back to 0.5.

EPSILON_MAP: Dict[str, float] = {
    "mean_duration": 0.3,
    "mean_magnitude": 0.3,
    "mean_severity": 0.05,
    "mean_avg_severity": 0.05,
    "peak_severity_month": 0.5,
    "max_duration": 1.0,
    "max_magnitude": 1.0,
    "worst_severity": 0.1,
    "frequency": 0.5,
}


# ---------------------------------------------------------------------------
# Experiment configuration dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ExperimentConfig:
    """All knobs for a single MOEA-FIND experiment.

    Attributes
    ----------
    ssi_timescale:
        SSI accumulation period in months. 3 is the manuscript default
        (short-term drought).
    objective_keys:
        Drought-characteristic keys optimised by the MOEA, in order.
        Must be keys handled by :func:`src.experiment_utils.compute_ssi_anti_ideal`
        — either a non-cyclic metric with a ``max_*`` entry there, or a
        cyclic calendar-month metric listed in
        :data:`src.objectives.CYCLIC_MONTH_KEYS`.
    anti_ideal_headroom:
        Multiplier applied to place the anti-ideal ``D*`` outside the
        feasible region. ``1.5`` gives a 50 % headroom for non-cyclic
        metrics (``1.5 × max_hist``) and ``12 × 1.5 = 18`` for cyclic
        calendar-month metrics. DD-11 requires ``D_j <= D*_j`` for every
        feasible point; do not reduce headroom below the value that
        guarantees this. See :func:`src.experiment_utils.compute_ssi_anti_ideal`.
    epsilon_map:
        Per-objective epsilon defaults; keyed by objective name. The
        Manhattan-norm epsilon is derived as the sum of the selected
        entries.
    nfe:
        Maximum function evaluations for the MOEA run.
    seed:
        Random seed used by both the generator and the MOEA backend.
    n_years_out:
        Synthetic trace length in water years.
    dv_mode:
        Kirsch-Borg DV injection mode. ``"residual"`` (default) perturbs
        the generator's standardised residuals; ``"index"`` uses the
        shuffled index space.
    algorithm:
        MOEA backend name. Overridden at runtime by
        ``MOEA_FIND_ALGORITHM`` env var when set (see
        :mod:`src.borg_runner`).
    n_islands:
        Borg multi-master island count. 1 = master-slave.
    checkpoint_freq:
        NFE between checkpoint dumps.
    constraints_json:
        Path to calibrated plausibility tolerances produced by
        ``workflows/diagnostics/diag_constraint_calibration.py``. None →
        unconstrained run.
    site_label:
        Site identifier used to key the constraint calibration JSON and
        (eventually) site-specific historical data loaders.
    """

    # SSI
    ssi_timescale: int = 3

    # Objectives
    objective_keys: Tuple[str, ...] = (
        "mean_duration",
        "mean_avg_severity",
        "peak_severity_month",
    )
    anti_ideal_headroom: float = 1.5

    # Epsilons — read through EPSILON_MAP by default. Override by
    # constructing an ExperimentConfig with a custom epsilon_map.
    epsilon_map: Dict[str, float] = field(
        default_factory=lambda: dict(EPSILON_MAP)
    )
    # Epsilon for the Manhattan-norm auxiliary objective f_{K+1}. The
    # earlier "sum rule" (``manhattan_eps = sum(per_axis_eps)``) was
    # coarser than any per-axis epsilon and suppressed archive spread
    # along the L1 axis. A small fixed value matches the per-axis
    # resolution and lets Borg keep solutions that differ in
    # ``sum(D_j)`` even when per-axis-epsilon-tied.
    manhattan_eps: float = 0.1

    # MOEA defaults.
    # NFE=100 is a plumbing smoke-test size (~seconds to a few minutes);
    # the resulting Pareto will be tiny and noisy. Bump to 50_000+ for
    # serious convergence checks and 500_000+ for production diagnostics.
    nfe: int = 100
    seed: int = 42
    n_years_out: int = 20
    dv_mode: str = "residual"
    algorithm: str = "borg_mm"
    n_islands: int = 1
    checkpoint_freq: int = 10_000

    # Constraints
    constraints_json: Optional[Path] = field(
        default_factory=lambda: (
            PROJECT_ROOT
            / "outputs"
            / "diag_constraint_calibration"
            / "calibrated_tolerances.json"
        )
    )

    # Anti-ideal D* placement. If ``anti_ideal_reference_json`` is set and
    # exists, its Pareto drought_metrics max is used as the reference_max
    # for NON-cyclic objectives (multiplied by ``anti_ideal_headroom`` for
    # the final D*). Cyclic calendar-month metrics always use
    # ``12 × headroom`` regardless, to preserve the DD-11 hyperplane
    # identity guarantee ``D_j <= D*_j``. When None (default), D* falls
    # back to ``historical_max × headroom`` for non-cyclic objectives.
    anti_ideal_reference_json: Optional[Path] = None

    # DV-space uniformity constraint (ablation arm; see src/constraints_dv.py)
    # Default constraint_mode is "hydrologic" so every existing run is
    # unchanged. "dv_uniform" swaps in a single DV-space constraint and
    # ignores ``constraints_json``. "none" disables all constraints.
    constraint_mode: str = "hydrologic"
    dv_uniformity_json: Optional[Path] = field(
        default_factory=lambda: (
            PROJECT_ROOT
            / "outputs"
            / "diag_dv_uniformity_calibration"
            / "calibrated_dv_tolerances.json"
        )
    )
    dv_uniformity_statistic: str = "l2_star"

    # Site
    site_label: str = "cannonsville"

    def with_overrides(self, **overrides) -> "ExperimentConfig":
        """Return a copy with the given attributes replaced.

        Use when a CLI arg or env var should adjust one knob without
        rebuilding the full config from scratch.
        """
        from dataclasses import replace

        return replace(self, **overrides)


#: Default experiment configuration. Scripts should treat this as the
#: baseline and apply CLI/env overrides on top via ``with_overrides``.
DEFAULT_EXPERIMENT: ExperimentConfig = ExperimentConfig()
