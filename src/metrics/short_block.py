"""Sub-annual / single-year hazard metrics for the T=1 / T=2 reframe (DD-15b).

Per-block hazard metrics designed for short trace lengths (1–2 water years)
where multi-event aggregation is unnecessary because each block typically
contains a single drought event (or none). Two metric tiers:

* **Tier H — raw-flow metrics.** Concrete physical flow integrals and
  extrema computed directly on the monthly trace. No SSI fitting needed;
  no event extraction. Stakeholder-interpretable in cfs / cfs·month.
* **Tier I — SSI-3 single-block metrics with burn-in.** SSI-3 calculator
  fitted on the **full** historical record; applied to a series that
  includes 3 burn-in months ahead of the T-year evaluation block so that
  every evaluation month has a valid 3-month rolling reference. Metrics
  are then aggregated only on the evaluation window.

Sign convention: every metric is constructed so that **larger value =
more drought-stressed**. Flow-based metrics carry a ``_neg`` suffix
because their underlying flow magnitude decreases under drought.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy.stats import norm

from src.metrics.extended import FullRecordRefs


#: Water-year month indices (block months 0..11 = Oct, Nov, Dec, Jan, ...,
#: Sep). Indices used for sub-annual aggregations.
WY_DJF_IDX = (2, 3, 4)            # Dec, Jan, Feb — wet season
WY_JJA_IDX = (8, 9, 10)           # Jun, Jul, Aug — peak demand period
WY_RECESSION_IDX = (8, 9, 10, 11) # Jun, Jul, Aug, Sep — summer recession window
WY_AMJ_IDX = (6, 7, 8)            # Apr, May, Jun — spring refill window
WY_OND_IDX = (0, 1, 2)            # Oct, Nov, Dec — fall recharge window
WY_AUG_IDX = 10                   # August (peak demand month)

#: Floor added before log-transform to keep ``log(0)`` finite without
#: distorting flows above ~1 cfs (typical Cannonsville monthly flow is
#: 100s–1000s of cfs).
_LOG_FLOW_EPS = 1.0


def _slope(values: np.ndarray) -> float:
    """Linear-regression slope (cfs/month)."""
    if values.size < 2:
        return 0.0
    x = np.arange(values.size, dtype=float)
    a = np.vstack([x, np.ones_like(x)]).T
    sol, *_ = np.linalg.lstsq(a, values, rcond=None)
    return float(sol[0])


def _rolling_min_3mo(monthly: np.ndarray) -> float:
    """Minimum 3-month rolling sum across the block."""
    s = pd.Series(monthly).rolling(window=3, min_periods=3).sum()
    if s.dropna().empty:
        return 0.0
    return float(s.dropna().min())


def _rolling_min_window(monthly: np.ndarray, window: int) -> float:
    """Minimum N-month rolling sum across the block."""
    s = pd.Series(monthly).rolling(window=window, min_periods=window).sum()
    if s.dropna().empty:
        return 0.0
    return float(s.dropna().min())


def compute_short_block_metrics(
    monthly_eval: np.ndarray,
    ssi3_full_window_series: Optional[pd.Series],
    refs: FullRecordRefs,
    *,
    eval_first_idx_in_ssi: int = 3,
) -> Dict[str, float]:
    """Compute Tier-H raw-flow + Tier-I SSI-3 short-block metrics.

    Args:
        monthly_eval: Evaluation block, 1D water-year-ordered monthly
            flows. Length must be ``T_years * 12`` with T in {1, 2}.
            Index 0 = October of the first evaluation water year.
        ssi3_full_window_series: SSI-3 series spanning *burn-in plus
            evaluation*, e.g. ``ssi3_full_window_series.values`` has
            length ``3 + T_years * 12``. The first ``eval_first_idx_in_ssi``
            entries are burn-in and are not used for metric aggregation.
            Pass ``None`` to skip Tier I.
        refs: Full-record reference statistics.
        eval_first_idx_in_ssi: Index of the first evaluation month
            inside ``ssi3_full_window_series``. Default 3 corresponds
            to a 3-month burn-in.

    Returns:
        Flat dict of metric values. Keys defined in
        :data:`SHORT_BLOCK_METRIC_NAMES`.
    """
    flows = np.asarray(monthly_eval, dtype=float).ravel()
    n_months = flows.size
    if n_months < 12 or n_months % 12 != 0:
        raise ValueError(
            f"monthly_eval must be a multiple of 12 months, got {n_months}"
        )
    n_years = n_months // 12
    flows_2d = flows.reshape(n_years, 12)

    out: Dict[str, float] = {}

    # ---- Tier H — raw-flow metrics ----------------------------------
    out["total_flow_neg"] = float(-flows.sum())

    # Per-water-year aggregations: take the worst (most-drought) year
    # if T_years > 1; for T=1 these collapse to that year's value.
    jja_per_year = flows_2d[:, list(WY_JJA_IDX)].sum(axis=1)
    djf_per_year = flows_2d[:, list(WY_DJF_IDX)].sum(axis=1)
    out["jja_total_neg"] = float(-jja_per_year.min())
    out["djf_total_neg"] = float(-djf_per_year.min())

    out["min_monthly_flow_neg"] = float(-flows.min())
    out["min_3mo_rolling_neg"] = float(-_rolling_min_3mo(flows))

    # Summer recession: linear slope across Jun-Sep, taken from the
    # year with the worst recession (most negative slope = steepest
    # drawdown). Negate so larger = more drought-stressed.
    recess_per_year = np.array([
        _slope(flows_2d[y, list(WY_RECESSION_IDX)])
        for y in range(n_years)
    ])
    out["summer_recession"] = float(-recess_per_year.min())

    # Annual-mean z-score: driest year's z-score against the full-record
    # annual-mean distribution. For T=1, this is just the year's z-score.
    annual_means = flows_2d.mean(axis=1)
    block_min_annual = float(annual_means.min())
    if refs.annual_mean_std > 1e-12:
        out["min_annual_zscore"] = float(
            (refs.annual_mean_mean - block_min_annual) / refs.annual_mean_std
        )
    else:
        out["min_annual_zscore"] = 0.0

    # Block flow-tail z-scores against the full-record monthly distribution
    if refs.monthly_std > 1e-12:
        block_q10 = float(np.percentile(flows, 10.0))
        block_q25 = float(np.percentile(flows, 25.0))
        block_q90 = float(np.percentile(flows, 90.0))
        out["q10_zscore"] = float(
            (refs.monthly_mean - block_q10) / refs.monthly_std
        )
        out["q25_zscore"] = float(
            (refs.monthly_mean - block_q25) / refs.monthly_std
        )
        # q90 z-score: positive when block's 90th percentile is below the
        # full-record monthly mean — captures whether peak flows arrived.
        out["q90_zscore"] = float(
            (refs.monthly_mean - block_q90) / refs.monthly_std
        )
    else:
        out["q10_zscore"] = 0.0
        out["q25_zscore"] = 0.0
        out["q90_zscore"] = 0.0

    # ---- Tier J — additional physically-motivated metrics (DRB-aware) ----
    # Spring refill window (Apr-May-Jun) volume — critical for NYC reservoir
    # filling per Pizzuto (2009). Take worst year's AMJ if T>1.
    amj_per_year = flows_2d[:, list(WY_AMJ_IDX)].sum(axis=1)
    out["amj_total_neg"] = float(-amj_per_year.min())

    # Fall recharge window (Oct-Nov-Dec) volume — early-winter antecedent.
    ond_per_year = flows_2d[:, list(WY_OND_IDX)].sum(axis=1)
    out["ond_total_neg"] = float(-ond_per_year.min())

    # Worst sustained 6-month period — half-year-scale low-flow envelope
    # (Smakhtin 2001 family of low-flow indices, generalised from 1-mo).
    out["min_6mo_rolling_neg"] = float(-_rolling_min_window(flows, 6))

    # August z-score — single-month anchor against full-record August
    # distribution. Late-summer month is peak demand for NYC.
    aug_flows = flows_2d[:, WY_AUG_IDX]
    block_aug_min = float(aug_flows.min())
    if (refs.per_month_std.size == 12
            and refs.per_month_std[WY_AUG_IDX] > 1e-12):
        mu_aug = float(refs.per_month_mean[WY_AUG_IDX])
        sd_aug = float(refs.per_month_std[WY_AUG_IDX])
        out["aug_zscore"] = float((mu_aug - block_aug_min) / sd_aug)
    else:
        out["aug_zscore"] = 0.0

    # ---- Tier I — SSI-3 single-block metrics with burn-in -----------
    if ssi3_full_window_series is None:
        out["min_ssi3_neg"] = 0.0
        out["total_deficit_ssi3_block"] = 0.0
        out["time_in_drought_fraction_block"] = 0.0
    else:
        ssi_arr = np.asarray(ssi3_full_window_series.values, dtype=float)
        # Discard burn-in: keep only the evaluation window's SSI values.
        eval_ssi = ssi_arr[eval_first_idx_in_ssi : eval_first_idx_in_ssi + n_months]
        finite = np.isfinite(eval_ssi)
        ssi_eval = eval_ssi[finite]
        if ssi_eval.size > 0:
            out["min_ssi3_neg"] = float(-ssi_eval.min())
            deficit_mask = ssi_eval < 0
            out["total_deficit_ssi3_block"] = float(
                -ssi_eval[deficit_mask].sum()
            ) if deficit_mask.any() else 0.0
            out["time_in_drought_fraction_block"] = float(
                (ssi_eval <= -1.0).sum() / ssi_eval.size
            )
        else:
            out["min_ssi3_neg"] = 0.0
            out["total_deficit_ssi3_block"] = 0.0
            out["time_in_drought_fraction_block"] = 0.0

    return out


#: Stable order for metric column layout in DataFrames / npz.
SHORT_BLOCK_METRIC_NAMES = (
    # Tier H — raw-flow (annual / seasonal volumes)
    "total_flow_neg",
    "jja_total_neg",
    "djf_total_neg",
    "min_monthly_flow_neg",
    "min_3mo_rolling_neg",
    "summer_recession",
    "min_annual_zscore",
    "q10_zscore",
    "q25_zscore",
    # Tier I — SSI-3 single-block (burn-in 3 months)
    "min_ssi3_neg",
    "total_deficit_ssi3_block",
    "time_in_drought_fraction_block",
    # Tier J — DRB-aware additions (DD-15c, 2026-04-30)
    "amj_total_neg",
    "ond_total_neg",
    "min_6mo_rolling_neg",
    "aug_zscore",
    "q90_zscore",
)

#: Concept tags for K-set diversity filtering.
SHORT_BLOCK_CONCEPT_MAP = {
    # Tier H
    "total_flow_neg": "annual_volume",
    "jja_total_neg": "summer_volume",
    "djf_total_neg": "winter_volume",
    "min_monthly_flow_neg": "month_extremum",
    "min_3mo_rolling_neg": "season_extremum",
    "summer_recession": "drawdown_rate",
    "min_annual_zscore": "annual_min",
    "q10_zscore": "flow_tail",
    "q25_zscore": "flow_tail",
    # Tier I
    "min_ssi3_neg": "ssi3_extremum",
    "total_deficit_ssi3_block": "ssi3_deficit_integral",
    "time_in_drought_fraction_block": "ssi3_persistence",
    # Tier J
    "amj_total_neg": "spring_volume",
    "ond_total_neg": "fall_volume",
    "min_6mo_rolling_neg": "half_year_extremum",
    "aug_zscore": "late_summer_month",
    "q90_zscore": "flow_high_tail",
}


@dataclass
class BlockExtractionPlan:
    """Metadata for extracting a single block + burn-in from a series."""
    T_years: int
    burnin_months: int = 3

    @property
    def total_months(self) -> int:
        return self.T_years * 12 + self.burnin_months

    def slice_indices(self, generator_n_years: int):
        """Return (burnin_start, eval_start, eval_end) indices into a
        flat monthly array of length ``generator_n_years*12`` such that:

        * ``[burnin_start, eval_end)`` is the SSI computation window
          (length ``total_months``).
        * ``[eval_start, eval_end)`` is the evaluation window
          (length ``T_years*12``).
        """
        eval_end = generator_n_years * 12
        eval_start = eval_end - self.T_years * 12
        burnin_start = eval_start - self.burnin_months
        if burnin_start < 0:
            raise ValueError(
                f"generator_n_years={generator_n_years} too small for "
                f"T_years={self.T_years} + burn-in={self.burnin_months}"
            )
        return burnin_start, eval_start, eval_end


# ---------------------------------------------------------------------------
# Bounded T=1 metric reformulation (DD-15c, 2026-05-01)
# ---------------------------------------------------------------------------
#
# The legacy short-block metrics (djf_total_neg, summer_recession,
# aug_zscore, ond_total_neg) failed DD-15c on four counts: unbounded on
# the non-drought side, scale-heterogeneous, anti-ideal not converging,
# and Pareto admitted flood-shaped corners (33,020 cfs vs historical
# max ≈ 7,700 cfs).
#
# This block adds a *candidate pool* of 24 sub-annual log-space window
# summaries, each emitted via two bounded mappings:
#
#   * Mapping G — Gaussian CDF of standardized log-flow:
#     ``D_g = Phi((s_j - mu_log_hist) / sigma_log_hist)``. Closed-form,
#     smooth, monotone, never saturates at any finite ``s_j``.
#
#   * Mapping E — empirical CDF with monotone bounded tail extrapolation:
#     inside the empirical envelope ``D_e = F_hist(s_j)``; outside,
#     extrapolate using a saturating exponential anchored to the local
#     empirical slope at the top/bottom k_tail samples. Two extreme
#     ``s_j`` values give *distinct* ``D_e`` values — no clipping, no
#     saturation flatten (DD-15c hard requirement: tail resolution must
#     be preserved so two DV sets with different extreme low-flows do
#     not collapse to the same metric).
#
# Both mappings produce ``D in [0, 1)`` so ``D*_j = 1.0`` (CONSTANT
# anti-ideal rule) preserves the L1 device contract ``D_j(x) <= D*_j``
# in zero re-anchor iterations.
#
# The candidate pool is *not* a production preset — it is the search
# pool for the diagnostic-driven (mapping x K=4 windows) selection in
# workflows/02_calibration. After Phase 3 selection, the winning K=4
# subset becomes a new preset ``short_block_drb_v2``.


#: Catalog of candidate window summaries for the bounded T=1 reformulation.
#: ``(kind, params)`` where ``kind in {"logmean", "logslope", "logsummin"}``:
#:   * ``logmean``: mean of ``log(Q + eps)`` over the index tuple. For T>1
#:     the per-year window mean is taken; the worst-year value is returned
#:     as drought-positive ``-min(per_year_mean)``. For T=1 collapses to
#:     ``-(window mean)``.
#:   * ``logslope``: least-squares slope of ``log(Q + eps)`` over the index
#:     tuple, taken on the worst (steepest) year and returned drought-
#:     positive (``-min(slopes)``). Captures relative drawdown rate.
#:   * ``logsummin``: minimum N-month rolling sum of ``log(Q + eps)`` across
#:     the entire block, returned drought-positive (``-min``).
#:
#: Indices are water-year ordering: 0=Oct, 1=Nov, ..., 11=Sep.
WINDOW_SPECS: Dict[str, Tuple[str, Tuple[int, ...]]] = {
    # 12 single calendar months (rank-equivalent to raw flow but kept in
    # log-space for numerical stability and figure clarity)
    "oct_logmean": ("logmean", (0,)),
    "nov_logmean": ("logmean", (1,)),
    "dec_logmean": ("logmean", (2,)),
    "jan_logmean": ("logmean", (3,)),
    "feb_logmean": ("logmean", (4,)),
    "mar_logmean": ("logmean", (5,)),
    "apr_logmean": ("logmean", (6,)),
    "may_logmean": ("logmean", (7,)),
    "jun_logmean": ("logmean", (8,)),
    "jul_logmean": ("logmean", (9,)),
    "aug_logmean": ("logmean", (10,)),
    "sep_logmean": ("logmean", (11,)),
    # 7 multi-month windows (DRB-physical seasons)
    "djf_logmean":   ("logmean", (2, 3, 4)),         # winter
    "mam_logmean":   ("logmean", (5, 6, 7)),         # spring
    "jja_logmean":   ("logmean", (8, 9, 10)),        # peak demand
    "amj_logmean":   ("logmean", (6, 7, 8)),         # spring refill
    "jjas_logmean":  ("logmean", (8, 9, 10, 11)),    # summer recession volume
    "mjjas_logmean": ("logmean", (7, 8, 9, 10, 11)), # extended summer
    "ond_logmean":   ("logmean", (0, 1, 2)),         # fall recharge
    # 3 recession candidates (log-slope captures relative drawdown,
    # which is the hydrologically meaningful drought signal in low-flow
    # regimes; linear-slope is dominated by absolute flow magnitude)
    "recession_jjas":   ("logslope", (8, 9, 10, 11)),
    "recession_mjjas":  ("logslope", (7, 8, 9, 10, 11)),
    "recession_amjjas": ("logslope", (6, 7, 8, 9, 10, 11)),
    # 2 sliding-window minima (sustained low-flow envelope)
    "min_3mo_logsum": ("logsummin", (3,)),
    "min_6mo_logsum": ("logsummin", (6,)),
}

#: Concept tags for K-set diversity (Phase 3 selection enforces ≥4
#: distinct concepts). Both mappings of the same window share the tag.
_WINDOW_CONCEPTS: Dict[str, str] = {
    "oct_logmean": "early_winter_month",
    "nov_logmean": "early_winter_month",
    "dec_logmean": "winter_month",
    "jan_logmean": "winter_month",
    "feb_logmean": "winter_month",
    "mar_logmean": "spring_month",
    "apr_logmean": "spring_month",
    "may_logmean": "spring_month",
    "jun_logmean": "early_summer_month",
    "jul_logmean": "summer_month",
    "aug_logmean": "summer_month",
    "sep_logmean": "late_summer_month",
    "djf_logmean": "winter_season",
    "mam_logmean": "spring_season",
    "jja_logmean": "summer_season",
    "amj_logmean": "spring_refill",
    "jjas_logmean": "summer_recession_vol",
    "mjjas_logmean": "extended_summer",
    "ond_logmean": "fall_recharge",
    "recession_jjas": "drawdown_rate",
    "recession_mjjas": "drawdown_rate",
    "recession_amjjas": "drawdown_rate",
    "min_3mo_logsum": "rolling_min_short",
    "min_6mo_logsum": "rolling_min_long",
}

#: Stable order of candidate metric names. 24 windows × 2 mappings = 48.
CANDIDATE_BOUNDED_METRIC_NAMES: Tuple[str, ...] = tuple(
    f"{w}_{m}"
    for w in WINDOW_SPECS
    for m in ("g", "e")
)

#: Concept tag per candidate metric (both mappings inherit the window's
#: concept; Phase 3 selection enforces concept diversity within K-sets).
CANDIDATE_BOUNDED_CONCEPT_MAP: Dict[str, str] = {
    f"{w}_{m}": _WINDOW_CONCEPTS[w]
    for w in WINDOW_SPECS
    for m in ("g", "e")
}


def _log_flows(flows: np.ndarray) -> np.ndarray:
    """``log(max(Q, 0) + eps)`` with eps = 1 cfs.

    Floor avoids ``log(0)`` while leaving flows above ~1 cfs (i.e., all
    realistic Cannonsville monthly flows) effectively undistorted in
    rank space.
    """
    return np.log(np.maximum(np.asarray(flows, dtype=float), 0.0) + _LOG_FLOW_EPS)


def _compute_window_summary(
    flows_2d: np.ndarray,
    kind: str,
    params: Tuple[int, ...],
) -> float:
    """Drought-positive log-space window summary on a (n_years, 12) block.

    Sign convention: returns a value that grows as drought severity grows.
    For T=1 the per-year worst-of reduces to that single year's value.
    """
    log_2d = _log_flows(flows_2d)
    n_years = log_2d.shape[0]
    if kind == "logmean":
        idx = list(params)
        per_year_mean = log_2d[:, idx].mean(axis=1)  # shape (n_years,)
        return float(-per_year_mean.min())
    if kind == "logslope":
        idx = list(params)
        slopes = np.array([_slope(log_2d[y, idx]) for y in range(n_years)])
        return float(-slopes.min())
    if kind == "logsummin":
        window_size = int(params[0])
        log_1d = log_2d.flatten()
        s = pd.Series(log_1d).rolling(window=window_size, min_periods=window_size).sum()
        s = s.dropna()
        if s.empty:
            return 0.0
        return float(-s.min())
    raise ValueError(f"Unknown window kind: {kind!r}")


def _apply_mapping_g(s_j: float, mu_log: float, sigma_log: float) -> float:
    """Mapping G: Gaussian CDF of standardized log-space summary.

    ``D = Phi((s_j - mu) / sigma) in (0, 1)``. Smooth, monotone, never
    saturates at any finite ``s_j``. Resolution at large ``|z|`` drops
    by design (Phi(4) ≈ 1−3.2e-5); diagnostic figures expose this so
    the user can compare against Mapping E.
    """
    if not np.isfinite(s_j):
        return 0.5
    if sigma_log <= 1e-12:
        return 0.5
    z = (s_j - mu_log) / sigma_log
    return float(norm.cdf(z))


def _apply_mapping_e(
    s_j: float,
    sorted_hist: np.ndarray,
    *,
    k_tail: int = 5,
) -> float:
    """Mapping E: empirical CDF with polynomial-tailed monotone extrapolation.

    Inside the empirical envelope ``[s_min, s_max]``: linear interpolation
    of the Hazen plotting positions ``i / (n+1)``.

    Above ``s_max`` (drought direction for drought-positive ``s``):
    Cauchy-shaped tail anchored to the empirical slope of the top
    ``k_tail`` samples. Closed form:

        rank = y_n + head * (slope_high * delta) / (head + slope_high * delta)

    where ``head = 1 - y_n``, ``y_n = n/(n+1)``, ``delta = s_j - s_max``.

    Properties:
    - At ``delta=0``: ``rank = y_n`` (continuity).
    - Local slope at ``delta=0`` is ``slope_high`` (matches the empirical
      tail slope, so curvature is hydrologically meaningful).
    - As ``delta → ∞``: ``rank → 1``, but the *gap* ``1 − rank`` decays
      polynomially (``∼ head²/(slope_high · delta)``) rather than
      exponentially. This is the key DD-15c property: two extreme
      drought values with delta=10 vs delta=100 give *distinct* D
      values that persist past the Gaussian-CDF float64 saturation
      point. No clipping; hydrologic signature in the extreme tail is
      preserved.

    Below ``s_min`` (flood direction): mirrored polynomial tail toward 0.

    Returns ``D ∈ (0, 1)``.
    """
    if not np.isfinite(s_j):
        return 0.5
    arr = np.asarray(sorted_hist, dtype=float)
    n = arr.size
    if n == 0:
        return 0.5
    ranks = np.arange(1, n + 1, dtype=float) / (n + 1.0)
    s_min = float(arr[0])
    s_max = float(arr[-1])

    if s_j <= s_min:
        kt = min(k_tail, n)
        denom = float(arr[kt - 1] - s_min)
        if kt < 2 or denom <= 1e-12:
            return float(ranks[0]) if s_j == s_min else float(ranks[0] * 0.5)
        slope_low = (ranks[kt - 1] - ranks[0]) / denom
        delta = s_min - s_j
        if slope_low <= 0 or ranks[0] <= 0:
            return float(max(1e-300, ranks[0] - slope_low * delta))
        # Polynomial tail toward 0 (Cauchy-shaped, mirrored):
        #   rank = ranks[0] * ranks[0] / (ranks[0] + slope_low * delta)
        # At delta=0: rank = ranks[0]; local rate -slope_low; as delta→∞,
        # rank → 0 with polynomial decay (preserves resolution).
        return float(ranks[0] * ranks[0] / (ranks[0] + slope_low * delta))

    if s_j >= s_max:
        kt = min(k_tail, n)
        denom = float(s_max - arr[-kt])
        head = 1.0 - ranks[-1]
        if kt < 2 or denom <= 1e-12 or head <= 0:
            base = float(ranks[-1])
            return min(base + 1e-9, 1.0 - 1e-15) if s_j > s_max else base
        slope_high = (ranks[-1] - ranks[-kt]) / denom
        delta = s_j - s_max
        if slope_high <= 0:
            return min(float(ranks[-1] + 1e-9), 1.0 - 1e-15)
        # Polynomial tail toward 1 (Cauchy-shaped). Local rate at
        # delta=0 = slope_high; gap (1 - rank) ~ head²/(slope_high · delta)
        # for large delta — preserves resolution past Gaussian saturation.
        return float(ranks[-1] + head * (slope_high * delta) / (head + slope_high * delta))

    return float(np.interp(s_j, arr, ranks))


def compute_candidate_bounded_metrics(
    monthly_eval: np.ndarray,
    family_refs: "BoundedFamilyRefs",
) -> Dict[str, float]:
    """Compute the 48 candidate bounded metrics on one block.

    Each :data:`WINDOW_SPECS` entry emits two values:

    * ``<window>_g`` — Mapping G (Gaussian CDF)
    * ``<window>_e`` — Mapping E (empirical CDF + tail extrapolation)

    Both ``D in [0, 1)``. The full candidate pool is used by the
    diagnostic-driven (mapping × K=4 windows) selection in
    ``workflows/02_calibration``; the winning subset becomes a new
    production preset.

    Args:
        monthly_eval: Water-year-ordered monthly flows, length T_years*12.
        family_refs: Per-window historical reference distributions.

    Returns:
        Flat dict with keys in :data:`CANDIDATE_BOUNDED_METRIC_NAMES`.
    """
    flows = np.asarray(monthly_eval, dtype=float).ravel()
    n_months = flows.size
    if n_months < 12 or n_months % 12 != 0:
        raise ValueError(
            f"monthly_eval must be a multiple of 12 months, got {n_months}"
        )
    flows_2d = flows.reshape(n_months // 12, 12)

    out: Dict[str, float] = {}
    for name, (kind, params) in WINDOW_SPECS.items():
        s_j = _compute_window_summary(flows_2d, kind, params)
        ref = family_refs.per_window.get(name) if family_refs is not None else None
        if ref is None:
            out[f"{name}_g"] = 0.5
            out[f"{name}_e"] = 0.5
            continue
        out[f"{name}_g"] = _apply_mapping_g(s_j, ref["mu"], ref["sigma"])
        out[f"{name}_e"] = _apply_mapping_e(s_j, ref["sorted"])
    return out


# Forward-declared for static analysis; the real class lives in
# :mod:`src.metrics.extended`. Imported here to avoid a circular
# import at module load time.
if False:  # pragma: no cover
    from src.metrics.extended import BoundedFamilyRefs  # noqa: F401
