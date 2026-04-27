"""Drought metric registry and configuration for MOEA-FIND.

Single source of truth for drought characteristic metadata: identifier,
human label, units, anti-ideal placement rule, default epsilon, and the
function that pulls the metric value out of the chars dict produced by
:func:`src.objectives.compute_ssi_drought_characteristics`.

Design intent: experiments specify a metric set by preset name (e.g.,
``"primary"``) or by a list of metric names. Anti-ideal placement,
epsilon values, and cyclic-vs-non-cyclic handling all flow from the
metric instances rather than from string-keyed lookup tables scattered
across modules.

See ``manuscript/governance/design_decisions.md`` §DD-04 for the open
question this module operationalises and §DD-11 for the L1 device that
consumes the metric set.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple, Union

import numpy as np


class AntiIdealRule(Enum):
    """How to place the anti-ideal coordinate :math:`D^*_j` for a given metric.

    HEADROOM_TIMES_MAX
        :math:`D^*_j = \\text{headroom} \\times \\text{max\\_historical}`. The
        default for non-cyclic unbounded-above metrics (severity, magnitude).
        Headroom (typically 1.5) is supplied by the caller of
        :func:`compute_anti_ideal`.
    CYCLIC_HEADROOM
        :math:`D^*_j = \\text{period} \\times \\text{headroom}`. Used for
        cyclic-month metrics so :math:`D^*` lies strictly outside the
        ``[1, 12]`` calendar.
    CONSTANT
        :math:`D^*_j = \\text{anti\\_ideal\\_constant}` (fixed on the metric
        instance). Used for metrics with a natural upper bound (e.g., a
        fraction in ``[0, 1]`` whose anti-ideal is ``1.0``).
    """

    HEADROOM_TIMES_MAX = "headroom_times_max"
    CYCLIC_HEADROOM = "cyclic_headroom"
    CONSTANT = "constant"


CharsDict = Mapping[str, float]


@dataclass(frozen=True)
class DroughtMetric:
    """A drought characteristic with the metadata required to use it as a
    MOEA-FIND objective.

    Each metric carries its own metadata so that anti-ideal placement,
    epsilon defaults, and labelling can dispatch on the metric instance
    rather than on string-keyed lookup tables in disconnected modules.

    Attributes:
        name: Stable identifier (used by presets, CLI flags, JSON output).
        label: Human-readable label for plots and logs.
        units: Units string for axis labels.
        is_cyclic: True for cyclic-calendar metrics. Affects anti-ideal
            placement. Continuous metrics in the primary set have
            ``is_cyclic = False``.
        epsilon: Default per-axis epsilon for the Borg archive.
        anti_ideal_rule: How to place :math:`D^*_j` for this metric.
        extract: Pulls this metric's value out of the chars dict produced
            by :func:`src.objectives.compute_ssi_drought_characteristics`.
            Returns ``0.0`` when the key is missing.
        max_partner: Optional name of the chars-dict key whose value is
            used as the historical-maximum reference for
            ``HEADROOM_TIMES_MAX`` placement. ``None`` means use the
            metric's own value in ``hist_chars``.
        anti_ideal_constant: Used when ``anti_ideal_rule = CONSTANT``;
            ignored otherwise.
    """

    name: str
    label: str
    units: str
    is_cyclic: bool
    epsilon: float
    anti_ideal_rule: AntiIdealRule
    extract: Callable[[CharsDict], float]
    max_partner: Optional[str] = None
    anti_ideal_constant: float = 1.0


def _extractor(key: str) -> Callable[[CharsDict], float]:
    """Build a chars-dict extractor for a single key.

    Returns ``0.0`` if the key is absent so that empty-Pareto edge cases
    don't raise during evaluation.
    """

    def _fn(chars: CharsDict) -> float:
        return float(chars.get(key, 0.0))

    return _fn


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
#
# Each metric instance is keyed by its stable ``name``. Add a new metric
# here, give it metadata, and it becomes available everywhere
# ``resolve_metric_set`` is consulted. Do not key dispatch on the metric
# string elsewhere — dispatch on the ``DroughtMetric`` instance.

REGISTRY: Dict[str, DroughtMetric] = {
    # -- Primary continuous metrics --------------------------------------
    "mean_severity": DroughtMetric(
        name="mean_severity",
        label="Mean drought severity",
        units="|SSI|",
        is_cyclic=False,
        epsilon=0.05,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("mean_severity"),
        max_partner="worst_severity",
    ),
    "mean_magnitude": DroughtMetric(
        name="mean_magnitude",
        label="Mean event cumulative deficit",
        units="|SSI|·month",
        is_cyclic=False,
        epsilon=0.3,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("mean_magnitude"),
        max_partner="max_magnitude",
    ),
    "time_in_drought_fraction": DroughtMetric(
        name="time_in_drought_fraction",
        label="Time-in-drought fraction",
        units="dimensionless",
        is_cyclic=False,
        epsilon=0.02,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=1.0,
        extract=_extractor("time_in_drought_fraction"),
    ),
    # -- Extreme-event variants -----------------------------------------
    "worst_severity": DroughtMetric(
        name="worst_severity",
        label="Worst event severity",
        units="|SSI|",
        is_cyclic=False,
        epsilon=0.1,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("worst_severity"),
    ),
    "max_magnitude": DroughtMetric(
        name="max_magnitude",
        label="Worst event cumulative deficit",
        units="|SSI|·month",
        is_cyclic=False,
        epsilon=1.0,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("max_magnitude"),
    ),
    # -- Trace-level FDC alternative -----------------------------------
    "q10_flow": DroughtMetric(
        # NB: q10_flow is the 10th percentile of monthly flows. Larger
        # values correspond to *less* drought-stressed traces. The L1
        # construction expects "more severe = larger objective", so the
        # extractor returns ``-q10`` (negated) to match. Anti-ideal is
        # placed at headroom × |min historical Q10|.
        name="q10_flow",
        label="Negated Q10 (10th-pct) monthly flow",
        units="cfs",
        is_cyclic=False,
        epsilon=10.0,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("q10_flow_neg"),
    ),
    # -- Other dict-resident metrics (kept for completeness) ----------
    "mean_avg_severity": DroughtMetric(
        name="mean_avg_severity",
        label="Mean within-event severity",
        units="|SSI|",
        is_cyclic=False,
        epsilon=0.05,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("mean_avg_severity"),
        max_partner="worst_severity",
    ),
    "frequency": DroughtMetric(
        name="frequency",
        label="Drought frequency",
        units="events/decade",
        is_cyclic=False,
        epsilon=0.5,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("frequency"),
    ),
    # -- Legacy (clustered or cyclic; retained for reproducibility) --
    "mean_duration": DroughtMetric(
        name="mean_duration",
        label="Mean drought duration",
        units="months",
        is_cyclic=False,
        epsilon=0.3,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("mean_duration"),
        max_partner="max_duration",
    ),
    "peak_severity_month": DroughtMetric(
        name="peak_severity_month",
        label="Peak severity month",
        units="month",
        is_cyclic=True,
        epsilon=0.5,
        anti_ideal_rule=AntiIdealRule.CYCLIC_HEADROOM,
        extract=_extractor("peak_severity_month"),
    ),
}


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------
#
# Named metric sets that experiments can refer to by name. Edit here to
# add a new ablation arm or to redefine the production set.

PRESETS: Dict[str, Tuple[str, ...]] = {
    "primary": (
        "mean_severity",
        "mean_magnitude",
        "time_in_drought_fraction",
    ),
    "extreme_event": (
        "worst_severity",
        "max_magnitude",
        "time_in_drought_fraction",
    ),
    "trace_fdc": (
        "q10_flow",
        "time_in_drought_fraction",
        "mean_magnitude",
    ),
    "legacy": (
        "mean_duration",
        "mean_avg_severity",
        "peak_severity_month",
    ),
}


# ---------------------------------------------------------------------------
# Resolution helpers
# ---------------------------------------------------------------------------


def resolve_metric_set(
    spec: Union[str, Sequence[str], Sequence[DroughtMetric]],
) -> Tuple[DroughtMetric, ...]:
    """Resolve a preset name, list of metric names, or list of metrics.

    Args:
        spec: Either (a) a preset name — key of :data:`PRESETS`, (b) a
            single metric name in :data:`REGISTRY`, (c) a sequence of
            metric names, or (d) a sequence of :class:`DroughtMetric`
            instances. The last form is a no-op pass-through; the other
            forms look up the registry.

    Returns:
        Tuple of :class:`DroughtMetric` in the supplied order.

    Raises:
        KeyError: Unknown preset or metric name.
    """
    if isinstance(spec, str):
        if spec in PRESETS:
            names: Sequence[str] = PRESETS[spec]
        elif spec in REGISTRY:
            names = (spec,)
        else:
            raise KeyError(
                f"Unknown metric set or metric name: {spec!r}. "
                f"Known presets: {sorted(PRESETS.keys())}. "
                f"Known metrics: {sorted(REGISTRY.keys())}."
            )
        return tuple(REGISTRY[n] for n in names)

    spec_tuple = tuple(spec)
    if len(spec_tuple) == 0:
        raise KeyError("Empty metric-set specification.")

    if isinstance(spec_tuple[0], DroughtMetric):
        return spec_tuple  # type: ignore[return-value]

    # Sequence of strings.
    out = []
    for name in spec_tuple:
        if not isinstance(name, str):
            raise TypeError(
                f"Mixed metric-set spec contains a non-string element: {name!r}"
            )
        if name not in REGISTRY:
            raise KeyError(
                f"Metric {name!r} not in registry. "
                f"Known metrics: {sorted(REGISTRY.keys())}."
            )
        out.append(REGISTRY[name])
    return tuple(out)


def metric_names(metric_set: Sequence[DroughtMetric]) -> Tuple[str, ...]:
    """Return the stable names from a tuple of metrics."""
    return tuple(m.name for m in metric_set)


def metric_labels(metric_set: Sequence[DroughtMetric]) -> Tuple[str, ...]:
    """Return the human labels from a tuple of metrics."""
    return tuple(m.label for m in metric_set)


def metric_epsilons(metric_set: Sequence[DroughtMetric]) -> Tuple[float, ...]:
    """Return the per-axis default epsilons from a tuple of metrics."""
    return tuple(m.epsilon for m in metric_set)


# ---------------------------------------------------------------------------
# Anti-ideal placement
# ---------------------------------------------------------------------------


def compute_anti_ideal(
    metric_set: Sequence[DroughtMetric],
    hist_chars: CharsDict,
    headroom: float = 1.5,
    feasible_maxes: Optional[Mapping[str, float]] = None,
) -> np.ndarray:
    """Build the DD-11 anti-ideal vector :math:`D^*` for a metric set.

    Each metric contributes one coordinate of :math:`D^*` using its
    :class:`AntiIdealRule`:

    - ``HEADROOM_TIMES_MAX``: ``headroom × reference_max``. The
      ``reference_max`` is ``feasible_maxes[name]`` if provided, else the
      historical value of the metric's ``max_partner`` (or its own value
      if no partner is declared) in ``hist_chars``.
    - ``CYCLIC_HEADROOM``: ``12 × headroom`` (period × headroom).
    - ``CONSTANT``: the metric's ``anti_ideal_constant``.

    DD-11 requires :math:`D_j(x) \\le D^*_j` for every feasible ``x`` and
    every ``j``; otherwise the affine-subspace identity fails silently.
    The cyclic and constant rules guarantee this by construction; for
    ``HEADROOM_TIMES_MAX`` the headroom factor (default 1.5) provides the
    safety margin against an unobserved larger value during search.

    Args:
        metric_set: Ordered tuple of :class:`DroughtMetric`.
        hist_chars: Drought-characteristics dict computed on the
            historical record.
        headroom: Safety factor :math:`\\ge 1` applied to the maximum.
        feasible_maxes: Optional override mapping metric name → empirical
            maximum. Used when an external Pareto archive has discovered
            characteristic values exceeding the historical record. Cyclic
            and constant metrics ignore this argument by design.

    Returns:
        ``np.ndarray`` of length ``len(metric_set)``.
    """
    out = np.zeros(len(metric_set), dtype=float)
    for i, metric in enumerate(metric_set):
        rule = metric.anti_ideal_rule
        if rule is AntiIdealRule.CYCLIC_HEADROOM:
            out[i] = 12.0 * headroom
            continue
        if rule is AntiIdealRule.CONSTANT:
            out[i] = float(metric.anti_ideal_constant)
            continue
        # HEADROOM_TIMES_MAX
        if feasible_maxes is not None and metric.name in feasible_maxes:
            ref = float(feasible_maxes[metric.name])
        else:
            partner = metric.max_partner or metric.name
            ref = float(hist_chars.get(partner, 0.0))
        if ref <= 0:
            ref = 10.0  # fallback: zero-event historical case
        out[i] = ref * headroom
    return out
