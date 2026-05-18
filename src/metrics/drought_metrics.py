"""Drought metric registry and configuration for MOEA-FIND.

Single source of truth for drought characteristic metadata: identifier,
human label, units, anti-ideal placement rule, default epsilon, and the
function that pulls the metric value out of the chars dict produced by
:func:`src.metrics.objectives.compute_ssi_drought_characteristics`.

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
            by :func:`src.metrics.objectives.compute_ssi_drought_characteristics`.
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
    # -- Short-block raw-flow metrics (T=1 event-discovery, DD-15) ----------
    "djf_total_neg": DroughtMetric(
        # Negated DJF (Dec-Jan-Feb) flow total: values are always ≤ 0 (flow is
        # non-negative), so D* = 0 is the physical bound (zero DJF flow).
        # HEADROOM_TIMES_MAX would use the negative historical max and fall back
        # to the emergency 10.0 sentinel; CONSTANT is the correct rule here.
        name="djf_total_neg",
        label="DJF total flow deficit",
        units="cfs·month",
        is_cyclic=False,
        epsilon=500.0,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=0.0,
        extract=_extractor("djf_total_neg"),
    ),
    "summer_recession": DroughtMetric(
        name="summer_recession",
        label="Summer recession slope",
        units="cfs/month",
        is_cyclic=False,
        epsilon=50.0,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("summer_recession"),
    ),
    "aug_zscore": DroughtMetric(
        name="aug_zscore",
        label="August flow z-score",
        units="σ",
        is_cyclic=False,
        epsilon=0.1,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("aug_zscore"),
    ),
    "ond_total_neg": DroughtMetric(
        # Negated OND (Oct-Nov-Dec) flow total: same sign convention as
        # djf_total_neg — always ≤ 0, so D* = 0 is the physical bound.
        name="ond_total_neg",
        label="OND total flow deficit",
        units="cfs·month",
        is_cyclic=False,
        epsilon=500.0,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=0.0,
        extract=_extractor("ond_total_neg"),
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
    # -- First-event family (T=10y first-SSI3-drought objective set) ---------
    # Anti-ideal anchored against the all-events historical maxes via
    # max_partner. The "≥1 event" hard constraint in the evaluator gates
    # candidates with first_event_present == 0 as infeasible, so these
    # extractors are only meaningful on feasible candidates.
    "first_event_duration": DroughtMetric(
        name="first_event_duration",
        label="First event duration",
        units="months",
        is_cyclic=False,
        epsilon=0.5,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("first_event_duration"),
        max_partner="max_duration",
    ),
    "first_event_severity": DroughtMetric(
        name="first_event_severity",
        label="First event peak depth",
        units="|SSI|",
        is_cyclic=False,
        epsilon=0.1,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("first_event_severity"),
        max_partner="worst_severity",
    ),
    "first_event_magnitude": DroughtMetric(
        name="first_event_magnitude",
        label="First event cumulative deficit",
        units="|SSI|·month",
        is_cyclic=False,
        epsilon=0.5,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("first_event_magnitude"),
        max_partner="max_magnitude",
    ),
    "first_event_start_month": DroughtMetric(
        name="first_event_start_month",
        label="First event start month",
        units="month",
        is_cyclic=True,
        epsilon=0.5,
        anti_ideal_rule=AntiIdealRule.CYCLIC_HEADROOM,
        extract=_extractor("first_event_start_month"),
    ),
    "first_event_peak_month": DroughtMetric(
        name="first_event_peak_month",
        label="First event peak month",
        units="month",
        is_cyclic=True,
        epsilon=0.5,
        anti_ideal_rule=AntiIdealRule.CYCLIC_HEADROOM,
        extract=_extractor("first_event_peak_month"),
    ),
    # Shape descriptors for the first event (added 2026-05-14):
    #   onset_intensification_rate captures flash-drought vs. slow-onset
    #   character (Otkin et al. 2018, Pendergrass et al. 2020); rising_limb
    #   fraction is the pure temporal-asymmetry shape factor (Lloyd-Hughes
    #   2014). Both are orthogonal to peak depth / cumulative deficit and
    #   replace first_event_start_month in the *_shape preset.
    "first_event_onset_intensification_rate": DroughtMetric(
        name="first_event_onset_intensification_rate",
        label="First event onset intensification rate",
        units="|SSI|/month",
        is_cyclic=False,
        epsilon=0.05,
        anti_ideal_rule=AntiIdealRule.HEADROOM_TIMES_MAX,
        extract=_extractor("first_event_onset_intensification_rate"),
        max_partner="max_intensification_rate",
    ),
    "first_event_rising_limb_fraction": DroughtMetric(
        # Bounded in (0, 1]: rising_months / duration. CONSTANT anti-ideal at
        # D* = 1.0 (DD-11 contract D_j ≤ D*_j holds by construction). Smaller
        # values = front-loaded (rapid onset); larger = back-loaded.
        name="first_event_rising_limb_fraction",
        label="First event rising-limb fraction",
        units="dimensionless",
        is_cyclic=False,
        epsilon=0.02,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=1.0,
        extract=_extractor("first_event_rising_limb_fraction"),
    ),
    # -- First-event fixed-anti-ideal variants (2026-05-14, user-specified) ---
    # Same chars-dict extractors as the first_event_* metrics above, but with
    # CONSTANT anti-ideal placement at user-chosen D*_j values and per-axis
    # epsilons. Used by the ``first_event_ssi3_t10_fixed`` preset. The ``_fc``
    # suffix ("fixed constant") keeps these distinct from the headroom/cyclic
    # anchored originals so both framings stay runnable (no clobber).
    #
    # DD-11 contract D_j(x) ≤ D*_j for the unbounded-above axes (severity,
    # magnitude, onset rate): the constants below are the user's tight
    # targets (4.5 / 80 / 4.5) inflated by the DD-11 default 1.5× headroom
    # to 6.75 / 120 / 6.75, so search exceeding the tight target no longer
    # breaks the affine-subspace identity. peak_month ≤ 12 < 13 and
    # rising_limb_fraction ≤ 1 < 1.1 are bounded-safe and unchanged.
    "first_event_severity_fc": DroughtMetric(
        name="first_event_severity_fc",
        label="First event peak depth",
        units="|SSI|",
        is_cyclic=False,
        epsilon=0.1,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=6.75,
        extract=_extractor("first_event_severity"),
    ),
    "first_event_magnitude_fc": DroughtMetric(
        name="first_event_magnitude_fc",
        label="First event cumulative deficit",
        units="|SSI|·month",
        is_cyclic=False,
        epsilon=1.0,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=120.0,
        extract=_extractor("first_event_magnitude"),
    ),
    "first_event_peak_month_fc": DroughtMetric(
        name="first_event_peak_month_fc",
        label="First event peak month",
        units="month",
        is_cyclic=True,
        epsilon=1.0,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=13.0,
        extract=_extractor("first_event_peak_month"),
    ),
    "first_event_onset_intensification_rate_fc": DroughtMetric(
        name="first_event_onset_intensification_rate_fc",
        label="First event onset intensification rate",
        units="|SSI|/month",
        is_cyclic=False,
        epsilon=0.1,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=6.75,
        extract=_extractor("first_event_onset_intensification_rate"),
    ),
    "first_event_rising_limb_fraction_fc": DroughtMetric(
        name="first_event_rising_limb_fraction_fc",
        label="First event rising-limb fraction",
        units="dimensionless",
        is_cyclic=False,
        epsilon=0.1,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=1.1,
        extract=_extractor("first_event_rising_limb_fraction"),
    ),
    # -- Coarse first-event variants (2026-05-15, windowed-pipeline dev) -----
    # Same chars extractors as _fc, but with COARSE epsilons (~1/10 of the
    # observed Pareto range from job 229763) for a small, fast dev ensemble,
    # and D* decreased for the severity-derived axes to reflect the interim
    # severity clip at 4.5 (see objectives.FIRST_EVENT_SEVERITY_CLIP):
    #   severity: clip 4.5 → D*=5.0 (4.5 + 0.5 strict-outside margin),
    #     ε=(4.5−1.0)/10=0.35
    #   onset_rate: severity-derived, same 4.5 ceiling → D*=5.0, ε≈0.44
    #   magnitude/peak_month/rising_limb: not severity-derived, D* unchanged,
    #     ε = obsRange/10 (mag 55.8→5.6, peak 11→1.1, rise 0.94→0.094).
    # Used by the ``first_event_ssi3_t10_coarse`` preset. ``_fcc`` keeps these
    # distinct from _fc / _fixed so all framings stay runnable.
    "first_event_severity_fcc": DroughtMetric(
        name="first_event_severity_fcc",
        label="First event peak depth",
        units="|SSI|",
        is_cyclic=False,
        epsilon=0.35,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=5.0,
        extract=_extractor("first_event_severity"),
    ),
    "first_event_magnitude_fcc": DroughtMetric(
        name="first_event_magnitude_fcc",
        label="First event cumulative deficit",
        units="|SSI|·month",
        is_cyclic=False,
        epsilon=5.6,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=120.0,
        extract=_extractor("first_event_magnitude"),
    ),
    "first_event_peak_month_fcc": DroughtMetric(
        name="first_event_peak_month_fcc",
        label="First event peak month",
        units="month",
        is_cyclic=True,
        epsilon=1.1,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=13.0,
        extract=_extractor("first_event_peak_month"),
    ),
    "first_event_onset_intensification_rate_fcc": DroughtMetric(
        name="first_event_onset_intensification_rate_fcc",
        label="First event onset intensification rate",
        units="|SSI|/month",
        is_cyclic=False,
        epsilon=0.44,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=5.0,
        extract=_extractor("first_event_onset_intensification_rate"),
    ),
    "first_event_rising_limb_fraction_fcc": DroughtMetric(
        name="first_event_rising_limb_fraction_fcc",
        label="First event rising-limb fraction",
        units="dimensionless",
        is_cyclic=False,
        epsilon=0.094,
        anti_ideal_rule=AntiIdealRule.CONSTANT,
        anti_ideal_constant=1.1,
        extract=_extractor("first_event_rising_limb_fraction"),
    ),
}


# ---------------------------------------------------------------------------
# Bounded T=1 candidate pool (DD-15c reformulation, 2026-05-01)
# ---------------------------------------------------------------------------
#
# 24 sub-annual log-space windows × 2 bounded mappings = 48 candidate
# metrics. The candidate pool feeds the diagnostic-driven (mapping ×
# K=4 windows) selection in workflows/02_calibration; the winning K=4
# subset is promoted to the production preset ``short_block_drb_v2``.
#
# All candidates use ``CONSTANT`` anti-ideal placement at ``D* = 1.0``:
# both Mapping G (Gaussian CDF) and Mapping E (empirical CDF + monotone
# tail extrapolation) produce ``D ∈ [0, 1)`` by construction, so the
# DD-11 contract ``D_j(x) ≤ D*_j`` holds in zero re-anchor iterations.
# Per-axis epsilon = 0.01 (≈100 archive bins per axis).

_BOUNDED_MAPPING_LABELS: Dict[str, str] = {
    "g": "Φ(z) of log-flow",
    "e": "empirical CDF + tail extrap",
}


def _register_bounded_candidates() -> None:
    """Add the 48 bounded-candidate metrics to :data:`REGISTRY`.

    Lazy import of :data:`WINDOW_SPECS` avoids a circular import at
    ``src.metrics.drought_metrics`` load time.
    """
    from src.metrics.short_block import WINDOW_SPECS

    for window_name in WINDOW_SPECS:
        for mapping_suffix, mapping_label in _BOUNDED_MAPPING_LABELS.items():
            metric_name = f"{window_name}_{mapping_suffix}"
            human_label = f"{window_name.replace('_', ' ')} ({mapping_label})"
            REGISTRY[metric_name] = DroughtMetric(
                name=metric_name,
                label=human_label,
                units="dimensionless",
                is_cyclic=False,
                epsilon=0.01,
                anti_ideal_rule=AntiIdealRule.CONSTANT,
                anti_ideal_constant=1.0,
                extract=_extractor(metric_name),
            )


_register_bounded_candidates()


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
    "short_block_drb": (
        "djf_total_neg",
        "summer_recession",
        "aug_zscore",
        "ond_total_neg",
    ),
    # DD-15c reformulation (2026-05-01): bounded log-space windows × Mapping G
    # (Gaussian CDF). Selected by diagnostic-driven K=4 search over a 48-
    # candidate pool with gates discrim ≥ 0.2, tail_res ≥ 1e-4, flood_corner
    # ≤ 0.05, median |ρ| ≤ 0.6, max |ρ| ≤ 0.55, ≥ 4 distinct concepts.
    # Winner: composite = 0.674, median ρ = 0.148, max ρ = 0.319.
    # Concepts: early_winter_month / winter_month / early_summer_month /
    # drawdown_rate. All four metrics use AntiIdealRule.CONSTANT with D* = 1.0.
    "short_block_drb_v2": (
        "oct_logmean_g",
        "dec_logmean_g",
        "jun_logmean_g",
        "recession_amjjas_g",
    ),
    # T=10y pivot (2026-05-06): characterise the FIRST critical SSI-3 drought
    # of each 10-year candidate sequence across five axes. Hard "≥1 event"
    # constraint enforced in the evaluator gates empty-event candidates.
    "first_event_ssi3_t10": (
        "first_event_duration",
        "first_event_severity",
        "first_event_magnitude",
        "first_event_start_month",
        "first_event_peak_month",
    ),
    # 2026-05-14 shape variant: drops start_month (near-redundant with
    # peak_month) and adds onset_intensification_rate + rising_limb_fraction
    # to characterize event shape. 6 objectives total. The original
    # first_event_ssi3_t10 preset stays runnable for back-compat.
    "first_event_ssi3_t10_shape": (
        "first_event_duration",
        "first_event_severity",
        "first_event_magnitude",
        "first_event_peak_month",
        "first_event_onset_intensification_rate",
        "first_event_rising_limb_fraction",
    ),
    # 2026-05-14 user-specified 5-objective first-event set with fixed
    # anti-ideal constants and per-axis epsilons. No duration axis. All
    # metrics use AntiIdealRule.CONSTANT (_fc variants). The original
    # first_event_ssi3_t10 / _shape presets stay runnable for back-compat.
    #   severity      D*=6.75  ε=0.1   (1.5× headroom over tight 4.5)
    #   magnitude     D*=120   ε=1     (1.5× headroom over tight 80)
    #   peak_month    D*=13    ε=1     (bounded-safe, cyclic)
    #   onset rate    D*=6.75  ε=0.1   (1.5× headroom over tight 4.5)
    #   rising limb   D*=1.1   ε=0.1   (bounded-safe)
    "first_event_ssi3_t10_fixed": (
        "first_event_severity_fc",
        "first_event_magnitude_fc",
        "first_event_peak_month_fc",
        "first_event_onset_intensification_rate_fc",
        "first_event_rising_limb_fraction_fc",
    ),
    # 2026-05-15 coarse variant for the windowed-pipeline build: same 5
    # first-event objectives, CONSTANT D*, but ~1/10-range epsilons →
    # small fast ensemble. Severity/onset D* lowered to 5.0 (interim 4.5
    # clip). Runs alongside _fixed. Slug suffix: 'coarse'.
    "first_event_ssi3_t10_coarse": (
        "first_event_severity_fcc",
        "first_event_magnitude_fcc",
        "first_event_peak_month_fcc",
        "first_event_onset_intensification_rate_fcc",
        "first_event_rising_limb_fraction_fcc",
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
