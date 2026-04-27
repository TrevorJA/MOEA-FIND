"""Centralized MOEA-FIND experiment configuration.

Single source of truth for the knobs that define a MOEA-FIND run:
metric set, anti-ideal placement, epsilons, NFE, trace length, DV
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
from typing import Optional, Tuple

# Module-level project root so config paths resolve consistently.
PROJECT_ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# Experiment configuration dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ExperimentConfig:
    """All knobs for a single MOEA-FIND experiment.

    The drought characteristics being optimised are selected by
    ``metric_set`` — a preset name from
    :data:`src.drought_metrics.PRESETS` (``"primary"``,
    ``"extreme_event"``, ``"trace_fdc"``, ``"legacy"``) or any single
    metric name from :data:`src.drought_metrics.REGISTRY`. Per-axis
    epsilons are read from the resolved metric instances. The cyclic
    handling, anti-ideal placement rule, and human labels for each
    metric live on the :class:`src.drought_metrics.DroughtMetric`
    instance, not on this config.

    Attributes
    ----------
    ssi_timescale:
        SSI accumulation period in months.
    metric_set:
        Preset name or single-metric name from
        :mod:`src.drought_metrics`. Default ``"primary"``: mean event
        severity, mean event cumulative deficit, time-in-drought
        fraction.
    anti_ideal_headroom:
        Multiplier applied when placing the anti-ideal ``D*`` for
        ``HEADROOM_TIMES_MAX`` metrics. ``1.5`` is the default;
        ``CYCLIC_HEADROOM`` and ``CONSTANT`` rules are unaffected. DD-11
        requires ``D_j <= D*_j``; do not reduce headroom below the value
        that guarantees this.
    manhattan_eps:
        Epsilon for the ``f_{K+1}`` Manhattan-norm auxiliary objective.
    nfe, seed, n_years_out, dv_mode, algorithm, n_islands,
    checkpoint_freq:
        MOEA backend knobs. ``nfe=100`` is a smoke-test size; bump to
        50 000+ for serious convergence checks.
    constraints_json:
        Path to calibrated plausibility tolerances. Used when
        ``constraint_mode = "hydrologic"``.
    constraint_mode:
        ``"hydrologic"``, ``"dv_uniform"``, or ``"none"``.
    dv_uniformity_json, dv_uniformity_statistic:
        Used when ``constraint_mode = "dv_uniform"``.
    anti_ideal_reference_json:
        Path to a prior ``results.json`` whose Pareto-archive maxima
        override the historical maxima for ``HEADROOM_TIMES_MAX``
        metrics. Cyclic and constant metrics ignore this.
    site_label:
        Site identifier used to key calibration JSONs and (eventually)
        site-specific historical data loaders.
    """

    # SSI
    ssi_timescale: int = 3

    # Metric set — preset name or single-metric name from drought_metrics.
    metric_set: str = "primary"
    anti_ideal_headroom: float = 1.5

    # Epsilon for the Manhattan-norm auxiliary objective f_{K+1}. The
    # earlier "sum rule" (``manhattan_eps = sum(per_axis_eps)``) was
    # coarser than any per-axis epsilon and suppressed archive spread
    # along the L1 axis. A small fixed value matches the per-axis
    # resolution and lets Borg keep solutions that differ in
    # ``sum(D_j)`` even when per-axis-epsilon-tied.
    manhattan_eps: float = 0.1

    # MOEA defaults.
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

    # Anti-ideal D* placement reference (optional; see DD-11 and
    # extract_pareto_maxes).
    anti_ideal_reference_json: Optional[Path] = None

    # DV-space uniformity constraint (production default per DD-14).
    constraint_mode: str = "dv_uniform"
    dv_uniformity_json: Optional[Path] = field(
        default_factory=lambda: (
            PROJECT_ROOT
            / "outputs"
            / "diag_dv_uniformity_calibration"
            / "calibrated_dv_tolerances.json"
        )
    )
    dv_uniformity_statistic: str = "ad"

    # Site
    site_label: str = "cannonsville"

    # ------------------------------------------------------------------
    # Derived properties (cheap; recomputed on each access).
    # ------------------------------------------------------------------

    @property
    def metrics(self) -> Tuple:
        """Resolve ``metric_set`` to a tuple of :class:`DroughtMetric`.

        Returns the tuple of metric instances that the MOEA optimises
        plus the metadata required for anti-ideal placement, epsilons,
        and labels.
        """
        from src.drought_metrics import resolve_metric_set
        return resolve_metric_set(self.metric_set)

    @property
    def objective_keys(self) -> Tuple[str, ...]:
        """Return the metric names from the resolved metric set.

        Backwards-compatible alias for code that still expects a tuple
        of string keys (e.g., for serialisation as ``objective_keys``
        in ``results.json``).
        """
        return tuple(m.name for m in self.metrics)

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
