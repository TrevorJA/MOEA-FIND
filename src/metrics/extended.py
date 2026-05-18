"""Extended candidate drought metrics for MOEA-FIND objective selection.

Builds an exploratory metric library used by
``workflows/02_calibration/metric_explorer.py`` to choose a low-correlation,
high-spread set of MOEA objectives. **Production REGISTRY in
:mod:`src.metrics.drought_metrics` is unchanged**; this module is read-only with
respect to the production hot-path.

Tiers (~28 candidates total before screening):

* **A** — SSI-3 event metrics (9): same definitions as
  :func:`src.metrics.objectives.compute_ssi_drought_characteristics`, recomputed via
  the shared event extractor below.
* **B** — SSI-12 event metrics (9): same 9 metrics at the longer
  hydrologic-drought timescale; suffixed ``_ssi12``.
* **C** — Recovery and inter-arrival (2): ``mean_recovery_time``,
  ``mean_drought_free_spell``.
* **D** — FDC-based low-flow (4): ``q10_flow_neg``, ``q25_flow_neg``,
  ``mean_annual_min``, ``cv_annual_min``.
* **E** — Run-theory threshold deficit at fixed Q80 (3): ``q80_events_per_decade``,
  ``q80_mean_deficit``, ``q80_max_deficit``. Q80 is computed once on the full
  historical record so values are comparable across blocks (Hisdal & Tallaksen
  2004 convention).
* **F** — Trend (1): ``sen_slope_annual_min``.

Sign convention: every metric is constructed so that **larger value = more
drought-stressed** (matching the L1 device). FDC-based metrics carry a
``_neg`` suffix because their underlying flow percentile decreases under
drought.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from scipy import stats

from src.metrics.objectives import (
    TIME_IN_DROUGHT_THRESHOLD,
    cyclic_mean,
    compute_drought_events,
)
from synhydro.droughts.ssi import get_drought_metrics


@dataclass
class FullRecordRefs:
    """Full-record reference statistics for cross-block-comparable metrics.

    Computed **once** on the historical record and reused on every T-block
    (historical or synthetic) so that the new "hazard-clean" metrics are
    anchored to a fixed reference rather than to per-block aggregates.
    Construct with :meth:`from_full_record`.
    """

    monthly_mean: float           #: μ of all historical monthly flows
    monthly_std: float            #: σ of all historical monthly flows
    annual_mean_mean: float       #: μ of full-record annual-mean flows
    annual_mean_std: float        #: σ of full-record annual-mean flows
    full_record_q10: float        #: 10th-percentile of full-record monthly flows
    full_record_q25: float        #: 25th-percentile of full-record monthly flows
    full_record_q80_threshold: float  #: 20th-percentile (Q80 deficit threshold)
    per_month_mean: np.ndarray = field(
        default_factory=lambda: np.zeros(12)
    )  #: per-calendar-month μ in water-year ordering (Oct..Sep)
    per_month_std: np.ndarray = field(
        default_factory=lambda: np.zeros(12)
    )  #: per-calendar-month σ in water-year ordering

    @classmethod
    def from_full_record(cls, monthly_1d: np.ndarray) -> "FullRecordRefs":
        """Compute reference stats from a 1D full-record monthly-flow array.

        ``monthly_1d`` is assumed water-year-aligned (first month = October).
        """
        flows = np.asarray(monthly_1d, dtype=float).ravel()
        n_full_yrs = flows.size // 12
        monthly_2d = flows[: n_full_yrs * 12].reshape(n_full_yrs, 12)
        annual_means = monthly_2d.mean(axis=1)
        per_month_mean = monthly_2d.mean(axis=0)
        per_month_std = (monthly_2d.std(axis=0, ddof=1)
                         if n_full_yrs > 1 else np.zeros(12))
        return cls(
            monthly_mean=float(flows.mean()),
            monthly_std=float(flows.std(ddof=1)) if flows.size > 1 else 0.0,
            annual_mean_mean=float(annual_means.mean()),
            annual_mean_std=float(annual_means.std(ddof=1))
                if annual_means.size > 1 else 0.0,
            full_record_q10=float(np.percentile(flows, 10.0)),
            full_record_q25=float(np.percentile(flows, 25.0)),
            full_record_q80_threshold=float(np.percentile(flows, 20.0)),
            per_month_mean=per_month_mean.astype(float),
            per_month_std=per_month_std.astype(float),
        )


@dataclass
class BoundedFamilyRefs:
    """Per-window historical reference for the bounded T=1 metric pool (DD-15c).

    Built once from the sliding-window T=1 historical record. For each
    candidate window in
    :data:`src.metrics.short_block.WINDOW_SPECS` we compute the
    drought-positive log-space summary on every historical water year,
    then store:

    * ``mu`` — mean of the historical sample (Mapping G centering).
    * ``sigma`` — std of the historical sample (Mapping G scaling).
    * ``sorted`` — sorted historical sample (Mapping E body + monotone
      tail extrapolation).

    These are loaded into the ``compute_candidate_bounded_metrics``
    evaluator once per experiment, then reused on every Pareto candidate.
    """

    per_window: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    @classmethod
    def from_full_record(cls, monthly_1d: np.ndarray) -> "BoundedFamilyRefs":
        """Build per-window references from a 1D water-year-ordered record.

        Each historical water year (12 months, Oct..Sep) becomes a single
        T=1 block; the candidate window summary is evaluated on that
        block; the resulting ``n_years``-sample distribution is stored
        per window.
        """
        # Imported here to avoid a circular import at module load time.
        from src.metrics.short_block import (
            WINDOW_SPECS,
            _compute_window_summary,
        )

        flows = np.asarray(monthly_1d, dtype=float).ravel()
        n_full_yrs = flows.size // 12
        if n_full_yrs == 0:
            return cls(per_window={})
        flows_2d = flows[: n_full_yrs * 12].reshape(n_full_yrs, 12)

        per_window: Dict[str, Dict[str, Any]] = {}
        for name, (kind, params) in WINDOW_SPECS.items():
            samples = np.empty(n_full_yrs, dtype=float)
            for y in range(n_full_yrs):
                samples[y] = _compute_window_summary(flows_2d[y:y + 1], kind, params)
            arr_sorted = np.sort(samples)
            mu = float(samples.mean())
            sigma = float(samples.std(ddof=1)) if samples.size > 1 else 0.0
            per_window[name] = {
                "mu": mu,
                "sigma": sigma,
                "sorted": arr_sorted,
            }
        return cls(per_window=per_window)


# ---------------------------------------------------------------------------
# Tier A/B — generic SSI event extractor (parameterised by suffix)
# ---------------------------------------------------------------------------


def compute_ssi_event_metrics(
    ssi_series: pd.Series,
    suffix: str = "",
    end_drought_threshold_months: int = 3,
) -> Dict[str, float]:
    """Compute the 9 standard SSI event-level metrics from an SSI series.

    Identical event definition to
    :func:`src.metrics.objectives.compute_ssi_drought_characteristics` (delegates
    to SynHydro's :func:`get_drought_metrics`); the only change is that
    every key is suffixed so SSI-3 and SSI-12 metrics can coexist in the
    same flat dict.

    Args:
        ssi_series: SSI :class:`pandas.Series` (DatetimeIndex required).
        suffix: Appended to every output key. Empty string reproduces the
            production SSI-3 names; ``"_ssi12"`` is used for the long-
            timescale candidate set.
        end_drought_threshold_months: Recovery hysteresis threshold (months
            of consecutive SSI > 0 to terminate a critical drought). Mirrors
            the production default.

    Returns:
        Dict with keys (each with ``suffix`` appended):
        ``frequency``, ``mean_duration``, ``max_duration``,
        ``mean_magnitude``, ``max_magnitude``, ``mean_severity``,
        ``worst_severity``, ``mean_avg_severity``,
        ``time_in_drought_fraction``, plus ``n_events`` (integer count,
        not a candidate metric — used by recovery/spell aggregations).
    """
    # Tier-A logic is shared with
    # objectives.compute_ssi_drought_characteristics; see src.metrics.ssi_common.
    from src.metrics.ssi_common import compute_ssi_tier_a

    return compute_ssi_tier_a(
        ssi_series,
        suffix=suffix,
        end_drought_threshold_months=end_drought_threshold_months,
    )


# ---------------------------------------------------------------------------
# Tier C — recovery & inter-arrival
# ---------------------------------------------------------------------------


def compute_recovery_metrics(
    ssi_series: pd.Series,
    end_drought_threshold_months: int = 3,
) -> Dict[str, float]:
    """Mean recovery time and mean drought-free spell length.

    * ``mean_recovery_time``: mean of SynHydro's ``recovery_period`` across
      events (months from peak severity to event end). Larger = slower
      recovery → more drought-stressed (matches sign convention).
    * ``mean_drought_free_spell``: mean gap in months between the end of
      one critical drought and the start of the next. Smaller spells →
      more frequent/persistent stress, so we **negate** the spell value
      (``mean_drought_free_spell_neg``) to keep the "larger = worse"
      convention. With <2 events the metric is ``0.0`` (no inter-arrival
      data); the screening step in the explorer will drop SSI-12 cases
      where this happens too often.
    """
    dm = get_drought_metrics(
        ssi_series,
        end_drought_threshold_months=end_drought_threshold_months,
    )
    if len(dm) == 0:
        return {
            "mean_recovery_time": 0.0,
            "mean_drought_free_spell_neg": 0.0,
        }

    mean_recovery = float(dm["recovery_period"].mean())

    if len(dm) < 2:
        spell_neg = 0.0
    else:
        starts = pd.to_datetime(dm["start"]).reset_index(drop=True)
        ends = pd.to_datetime(dm["end"]).reset_index(drop=True)
        # Months between event i's end and event (i+1)'s start.
        gaps_months = []
        for i in range(len(dm) - 1):
            gap = (
                (starts.iloc[i + 1].year - ends.iloc[i].year) * 12
                + (starts.iloc[i + 1].month - ends.iloc[i].month)
                - 1
            )
            gaps_months.append(max(gap, 0))
        spell_neg = -float(np.mean(gaps_months))

    return {
        "mean_recovery_time": mean_recovery,
        "mean_drought_free_spell_neg": spell_neg,
    }


# ---------------------------------------------------------------------------
# Tier D — FDC-based low-flow
# ---------------------------------------------------------------------------


def compute_fdc_metrics(
    monthly_1d: np.ndarray,
    monthly_2d: np.ndarray,
) -> Dict[str, float]:
    """FDC tail and annual-minimum statistics.

    Sign convention: percentile-based metrics are negated (``_neg`` suffix)
    so larger value = lower flow = more drought-stressed. ``mean_annual_min``
    is also negated for the same reason; ``cv_annual_min`` is dimensionless
    and reported directly (a higher CV indicates more variable annual
    minima, weakly correlated with drought stress on its own).
    """
    flows = np.asarray(monthly_1d, dtype=float).ravel()
    if flows.size == 0:
        return {
            "q10_flow_neg": 0.0,
            "q25_flow_neg": 0.0,
            "mean_annual_min_neg": 0.0,
            "cv_annual_min": 0.0,
        }

    q10 = float(np.percentile(flows, 10.0))
    q25 = float(np.percentile(flows, 25.0))

    annual_min = monthly_2d.min(axis=1)  # shape (n_years,)
    if annual_min.size == 0:
        mean_min = 0.0
        cv_min = 0.0
    else:
        mean_min = float(annual_min.mean())
        std_min = float(annual_min.std(ddof=1)) if annual_min.size > 1 else 0.0
        cv_min = float(std_min / abs(mean_min)) if abs(mean_min) > 1e-12 else 0.0

    return {
        "q10_flow_neg": -q10,
        "q25_flow_neg": -q25,
        "mean_annual_min_neg": -mean_min,
        "cv_annual_min": cv_min,
    }


# ---------------------------------------------------------------------------
# Tier E — run-theory threshold deficit at fixed Q80
# ---------------------------------------------------------------------------


def compute_threshold_deficit_metrics(
    monthly_1d: np.ndarray,
    q80_threshold: float,
) -> Dict[str, float]:
    """Hisdal & Tallaksen (2004) threshold-deficit metrics at fixed Q80.

    Q80 is supplied by the caller (computed **once on the full historical
    record** so values are comparable across blocks). Run-theory event
    extraction reuses :func:`src.metrics.objectives.compute_drought_events`.

    Returns:
        Dict with ``q80_events_per_decade``, ``q80_mean_deficit``,
        ``q80_max_deficit``. Deficit values are positive (sum of
        ``threshold − flow`` over the event); larger = more drought.
    """
    flows = np.asarray(monthly_1d, dtype=float).ravel()
    n_months = flows.size
    n_years = n_months / 12.0 if n_months > 0 else 0.0
    events = compute_drought_events(flows, q80_threshold)
    if len(events) == 0 or n_years <= 0:
        return {
            "q80_events_per_decade": 0.0,
            "q80_mean_deficit": 0.0,
            "q80_max_deficit": 0.0,
        }

    deficits = np.array([e["deficit"] for e in events], dtype=float)
    return {
        "q80_events_per_decade": float(len(events) / n_years * 10),
        "q80_mean_deficit": float(deficits.mean()),
        "q80_max_deficit": float(deficits.max()),
    }


# ---------------------------------------------------------------------------
# Tier G — hazard-clean continuous integrals (DD-15 user-driven, 2026-04-29)
# ---------------------------------------------------------------------------
#
# The Tier G metrics are added to address four methodological concerns:
#
# 1. **Block-normalisation.** Tier D's ``cv_annual_min`` divides by the
#    block's own mean — a wet and a dry block can collide on CV. Replaced
#    here by ``min_annual_zscore``, which scores the driest annual mean
#    in a block against the *full-record* annual-mean distribution.
# 2. **Event averaging.** Tier A's ``mean_severity`` and friends average
#    across all events inside the block, so a 5-yr block with one extreme
#    plus one minor event lands at the midpoint and the resulting MOEA
#    axis loses its operational interpretation. ``total_deficit_ssi3``
#    sums |SSI-3| over all deficit-months without event extraction;
#    every deficit month contributes its full depth additively.
# 3. **Single-month extremum saturation.** Tier A's ``worst_severity``
#    saturates because most overlapping blocks happen to contain the
#    same record-worst monthly SSI value. ``min_36mo_ssi3`` uses a
#    36-month rolling mean (smaller than the 60-month T=5 block) so the
#    block-min of the rolling series varies smoothly even when the
#    single worst month is shared.
# 4. **Block-internal flow quantiles.** Tier D's ``q10_flow_neg`` and
#    ``q25_flow_neg`` are absolute flow values but lose their reference
#    frame across blocks; ``q10_zscore`` and ``q25_zscore`` rescale to
#    the full-record monthly-flow distribution.
#
# All Tier G metrics are continuous, full-record-anchored, non-averaged,
# and non-stratified. Sign convention: larger value = more drought-stressed.


def compute_hazard_clean_metrics(
    monthly_1d: np.ndarray,
    monthly_2d: np.ndarray,
    ssi3_series: pd.Series,
    refs: FullRecordRefs,
    *,
    rolling_window_months: int = 36,
) -> Dict[str, float]:
    """Tier G — hazard-clean continuous metrics anchored to the full record."""
    flows = np.asarray(monthly_1d, dtype=float).ravel()
    block_2d = np.asarray(monthly_2d, dtype=float)
    ssi3 = np.asarray(ssi3_series.values, dtype=float)

    # 1) Total integrated SSI-3 deficit (sum of |SSI| where SSI < 0).
    deficit_mask = np.isfinite(ssi3) & (ssi3 < 0)
    total_deficit_ssi3 = float(-ssi3[deficit_mask].sum()) if deficit_mask.any() else 0.0

    # 2) Block-min of 36-month rolling SSI-3 mean (negated for sign convention).
    if ssi3.size >= rolling_window_months:
        rolling = pd.Series(ssi3).rolling(
            window=rolling_window_months, min_periods=rolling_window_months,
        ).mean().dropna().values
        if rolling.size > 0:
            min_36mo_ssi3 = float(-np.nanmin(rolling))
        else:
            min_36mo_ssi3 = 0.0
    else:
        min_36mo_ssi3 = 0.0

    # 3) Z-score of the driest annual-mean year in the block, against the
    #    full-record annual-mean distribution.
    if block_2d.shape[0] > 0:
        block_annual_means = block_2d.mean(axis=1)
        block_min_annual = float(np.min(block_annual_means))
        if refs.annual_mean_std > 1e-12:
            min_annual_zscore = float(
                (refs.annual_mean_mean - block_min_annual) / refs.annual_mean_std
            )
        else:
            min_annual_zscore = 0.0
    else:
        min_annual_zscore = 0.0

    # 4) Z-scores of block flow tail percentiles against the full-record
    #    monthly-flow distribution.
    if flows.size > 0 and refs.monthly_std > 1e-12:
        block_q10 = float(np.percentile(flows, 10.0))
        block_q25 = float(np.percentile(flows, 25.0))
        q10_zscore = float((refs.monthly_mean - block_q10) / refs.monthly_std)
        q25_zscore = float((refs.monthly_mean - block_q25) / refs.monthly_std)
    else:
        q10_zscore = 0.0
        q25_zscore = 0.0

    # 5) Total Q80 deficit (sum of (Q80 - flow) where flow < Q80) using the
    #    full-record Q80 threshold. Continuous integral, no event extraction.
    threshold = refs.full_record_q80_threshold
    deficit_below = np.maximum(threshold - flows, 0.0)
    q80_total_deficit = float(np.sum(deficit_below))

    return {
        "total_deficit_ssi3": total_deficit_ssi3,
        "min_36mo_ssi3": min_36mo_ssi3,
        "min_annual_zscore": min_annual_zscore,
        "q10_zscore": q10_zscore,
        "q25_zscore": q25_zscore,
        "q80_total_deficit": q80_total_deficit,
    }


# ---------------------------------------------------------------------------
# Tier F — trend
# ---------------------------------------------------------------------------


def compute_trend_metrics(monthly_2d: np.ndarray) -> Dict[str, float]:
    """Sen's slope of within-block annual minima.

    Captures within-block trend, orthogonal to mean/spread axes. Sign
    convention: a more negative slope (annual minima decreasing through
    the block) indicates worsening drought, so we **negate** the slope
    so larger = more drought-stressed.

    Falls back to 0.0 if the block has fewer than 3 annual minima or if
    Theil-Sen returns a non-finite slope (constant input).
    """
    annual_min = np.asarray(monthly_2d, dtype=float).min(axis=1)
    if annual_min.size < 3:
        return {"sen_slope_annual_min_neg": 0.0}
    slope, _, _, _ = stats.theilslopes(annual_min, np.arange(annual_min.size))
    if not np.isfinite(slope):
        return {"sen_slope_annual_min_neg": 0.0}
    return {"sen_slope_annual_min_neg": float(-slope)}


# ---------------------------------------------------------------------------
# Single dispatch — full candidate library on one block
# ---------------------------------------------------------------------------

#: Concept tag per candidate metric. Used by the explorer to enforce
#: concept-diversity in K-set selection: different concepts = genuinely
#: distinct hydrologic dimensions, irrespective of correlation algebra.
#: SSI-3 and SSI-12 versions of the same characteristic carry different
#: concept tags (``severity`` vs ``severity_long``) so a recommendation
#: can include both timescales when their numeric correlation is low.
CONCEPT_MAP = {
    # Tier A — SSI-3 events
    "frequency": "frequency",
    "mean_duration": "duration",
    "max_duration": "duration",
    "mean_magnitude": "magnitude",
    "max_magnitude": "magnitude",
    "mean_severity": "severity",
    "worst_severity": "severity",
    "mean_avg_severity": "severity",
    "time_in_drought_fraction": "persistence_share",
    # Tier B — SSI-12 events (separate concept namespace by suffix)
    "frequency_ssi12": "frequency_long",
    "mean_duration_ssi12": "duration_long",
    "max_duration_ssi12": "duration_long",
    "mean_magnitude_ssi12": "magnitude_long",
    "max_magnitude_ssi12": "magnitude_long",
    "mean_severity_ssi12": "severity_long",
    "worst_severity_ssi12": "severity_long",
    "mean_avg_severity_ssi12": "severity_long",
    "time_in_drought_fraction_ssi12": "persistence_share_long",
    # Tier C
    "mean_recovery_time": "recovery",
    "mean_drought_free_spell_neg": "inter_arrival",
    # Tier D — FDC
    "q10_flow_neg": "flow_tail",
    "q25_flow_neg": "flow_tail",
    "mean_annual_min_neg": "annual_min",
    "cv_annual_min": "annual_min",
    # Tier E — Q80 run-theory deficit
    "q80_events_per_decade": "frequency",  # identical concept to SSI-3 frequency
    "q80_mean_deficit": "deficit_volume",
    "q80_max_deficit": "deficit_volume",
    # Tier F — trend
    "sen_slope_annual_min_neg": "trend",
    # Tier G — hazard-clean continuous integrals (DD-15)
    "total_deficit_ssi3": "volume_integral",
    "min_36mo_ssi3": "sustained_severity",
    "min_annual_zscore": "annual_min",
    "q10_zscore": "flow_tail",
    "q25_zscore": "flow_tail",
    "q80_total_deficit": "volume_integral_q80",
}

#: Suffix used on every SSI-12 event metric to namespace it away from SSI-3.
SSI12_SUFFIX = "_ssi12"

#: Names of all candidate metrics produced by :func:`compute_all_candidates`,
#: in stable order. Used by the explorer driver to assemble the per-block
#: DataFrame and by downstream selection.
CANDIDATE_METRIC_NAMES = (
    # Tier A — SSI-3 (no suffix; matches existing REGISTRY names)
    "frequency",
    "mean_duration",
    "max_duration",
    "mean_magnitude",
    "max_magnitude",
    "mean_severity",
    "worst_severity",
    "mean_avg_severity",
    "time_in_drought_fraction",
    # Tier B — SSI-12
    "frequency_ssi12",
    "mean_duration_ssi12",
    "max_duration_ssi12",
    "mean_magnitude_ssi12",
    "max_magnitude_ssi12",
    "mean_severity_ssi12",
    "worst_severity_ssi12",
    "mean_avg_severity_ssi12",
    "time_in_drought_fraction_ssi12",
    # Tier C — recovery & inter-arrival
    "mean_recovery_time",
    "mean_drought_free_spell_neg",
    # Tier D — FDC
    "q10_flow_neg",
    "q25_flow_neg",
    "mean_annual_min_neg",
    "cv_annual_min",
    # Tier E — Q80 threshold deficit
    "q80_events_per_decade",
    "q80_mean_deficit",
    "q80_max_deficit",
    # Tier F — trend
    "sen_slope_annual_min_neg",
    # Tier G — hazard-clean continuous integrals (DD-15)
    "total_deficit_ssi3",
    "min_36mo_ssi3",
    "min_annual_zscore",
    "q10_zscore",
    "q25_zscore",
    "q80_total_deficit",
)

#: Hazard-clean subset (Tier G + the SSI-3 time-in-drought fraction). These
#: are the metrics recommended for K-set enumeration at small T per
#: 2026-04-29 user review of DD-15 diagnostics: every metric is continuous,
#: full-record-anchored, non-averaged across events, and not stratified at
#: integer values.
HAZARD_CLEAN_METRICS = (
    "total_deficit_ssi3",
    "min_36mo_ssi3",
    "min_annual_zscore",
    "q10_zscore",
    "q25_zscore",
    "q80_total_deficit",
    "time_in_drought_fraction",
)


def compute_all_candidates(
    monthly_1d: np.ndarray,
    monthly_2d: np.ndarray,
    ssi3_calc,
    ssi12_calc,
    q80_threshold: float,
    *,
    full_record_refs: Optional[FullRecordRefs] = None,
    start_date: str = "2100-01-01",
) -> Dict[str, float]:
    """Run every candidate-metric extractor on one T-year block.

    ``ssi3_calc`` and ``ssi12_calc`` are pre-fitted SynHydro
    :class:`SSI` calculators (fit on the **full** historical record so
    every block sees the same calibrated distribution — DD-11 lock-in).
    ``q80_threshold`` is also computed once on the full record.

    Args:
        monthly_1d: 1D monthly flows for the block, length ``T_years*12``.
        monthly_2d: Same flows reshaped to ``(T_years, 12)``.
        ssi3_calc: Pre-fitted SSI-3 calculator.
        ssi12_calc: Pre-fitted SSI-12 calculator.
        q80_threshold: Q80 of the **full** historical record.
        start_date: Dummy date for the synthetic pandas index passed to
            the SSI calculators (SSI is a stationary transform on the
            fitted distributions, so any well-formed date works).

    Returns:
        Flat dict containing every name in
        :data:`CANDIDATE_METRIC_NAMES` plus the integer event counts
        ``n_events`` and ``n_events_ssi12`` (kept for diagnostics; not
        used as candidates).
    """
    from src.metrics.objectives import flows_to_series

    series = flows_to_series(monthly_1d, start_date=start_date)
    ssi3 = ssi3_calc.transform(series)
    ssi12 = ssi12_calc.transform(series)

    out: Dict[str, float] = {}
    out.update(compute_ssi_event_metrics(ssi3, suffix=""))
    out.update(compute_ssi_event_metrics(ssi12, suffix=SSI12_SUFFIX))
    out.update(compute_recovery_metrics(ssi3))
    out.update(compute_fdc_metrics(monthly_1d, monthly_2d))
    out.update(compute_threshold_deficit_metrics(monthly_1d, q80_threshold))
    out.update(compute_trend_metrics(monthly_2d))
    if full_record_refs is not None:
        out.update(compute_hazard_clean_metrics(
            monthly_1d, monthly_2d, ssi3, full_record_refs,
        ))
    else:
        # Tier G unavailable without full-record refs; emit zeros so the
        # output shape is stable for downstream consumers that expect
        # CANDIDATE_METRIC_NAMES coverage.
        for name in (
            "total_deficit_ssi3", "min_36mo_ssi3", "min_annual_zscore",
            "q10_zscore", "q25_zscore", "q80_total_deficit",
        ):
            out[name] = 0.0
    return out
