"""Magnitude-Varying Sensitivity Analysis (MV-SA) in drought-hazard space.

Adaptation of Hadjimichael et al. (2020) MV-SA to the MOEA-FIND
diagnostic question: instead of asking which *parameters* drive
shortage at each magnitude (parameter-space MV-SA), we ask which
*drought-hazard characteristics* (D_1,...,D_K) drive operational
stress at each percentile of an operational hazard outcome (e.g.,
NYC minimum reservoir storage). Factor space and magnitude axis are
both inhabitants of the empirical drought-hazard / outcome space the
project structures coverage over.

Three response forms are supported:

    exceedance              -- I(tau) = 1 if M <= M_quantile(tau)
                               else 0. Cheap binary indicator. Known
                               to be degenerate when the archive has
                               no inner stochastic ensemble: the
                               binary Y has small variance and SALib
                               delta returns a noise-floor result
                               near 0.13 for all factors. Kept for
                               diagnostic comparison only; not used
                               for production figures.

    conditional             -- For a secondary outcome R, restrict
                               the sample to a window of size
                               ``window_frac * N`` realizations
                               centered (by rank) on M's tau-quantile,
                               then run SA of D on R within that
                               subset. Diagnostic for "given a year
                               sits at this severity of the primary
                               hazard, which D drives the secondary
                               outcome?"

    within_trace_percentile -- Y_i(tau) = np.percentile(
                               trace_series[i, :], 100 * tau).
                               Direct port of Hadjimichael Variant 2
                               (the magnitude-of-shortage version,
                               implemented in the project's tutorial
                               repo). Inputs require an (N, T) per-
                               trace timeseries -- e.g., the 20-year
                               annual NYC min storage series produced
                               by ``src.hydrology.precompute_trace_series``. Y
                               is continuous, varies smoothly with
                               tau, and avoids the small-N binary
                               degeneracy of ``exceedance``. This is
                               the production headline form.

A *control* uniform-random factor is optionally appended (per
Hadjimichael's recommendation) so the magnitude-varying noise floor
is visible alongside the real factors.

Public surface (consumed by ``workflows/09_magnitude_varying_sa/``):

- :func:`compute_mv_sa`
- :func:`stacked_share`

The module reuses :mod:`src.sensitivity.sensitivity` per-method ``compute_*`` and
``bootstrap_*`` helpers per percentile slice; nothing in
:mod:`SALib` is invoked directly here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from src.sensitivity.sensitivity import (
    HEADLINE_INDEX,
    METHODS,
    SAMethod,
    bootstrap_indices,
    bootstrap_rank_stability,
    resolve_method,
)

CONTROL_FACTOR_NAME = "control_uniform"


# ---------------------------------------------------------------------------
# Helpers: control factor, response builders, percentile windows
# ---------------------------------------------------------------------------


def _augment_with_control(
    X: np.ndarray,
    factor_names: Sequence[str],
    *,
    seed: int,
) -> Tuple[np.ndarray, List[str]]:
    """Append a uniform-random control column to X.

    Hadjimichael et al. (2020) include a placebo factor so the
    magnitude-varying noise floor is empirically visible: any factor
    whose index sits below the control band cannot be distinguished
    from sampling artefact.
    """
    rng = np.random.default_rng(seed)
    ctrl = rng.uniform(0.0, 1.0, size=X.shape[0]).reshape(-1, 1)
    X_aug = np.concatenate([X, ctrl], axis=1)
    names_aug = list(factor_names) + [CONTROL_FACTOR_NAME]
    return X_aug, names_aug


def _exceedance_indicator(M: np.ndarray, tau: float) -> np.ndarray:
    """Binary I(tau) = 1{M <= quantile(M, tau)} as float."""
    threshold = float(np.quantile(M, tau))
    return (M <= threshold).astype(float)


def _conditional_window_idx(
    M: np.ndarray,
    tau: float,
    *,
    window_frac: float,
) -> np.ndarray:
    """Index mask of the ``window_frac * N`` rows whose M-rank is closest to tau.

    Returns a boolean array of length len(M).
    """
    n = len(M)
    target_rank = tau * (n - 1)
    window = max(2, int(round(window_frac * n)))
    half = window // 2
    # rank realization by M, then keep ranks in [target_rank - half, target_rank + half].
    order = np.argsort(M, kind="mergesort")
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(n)
    lo = max(0, int(round(target_rank - half)))
    hi = min(n, lo + window)
    lo = max(0, hi - window)
    mask = (ranks >= lo) & (ranks < hi)
    return mask


# ---------------------------------------------------------------------------
# Core engine
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MVSAConfig:
    """Frozen configuration for a single MV-SA invocation.

    Stored alongside results so a downstream consumer can reconstruct
    the percentile sweep without re-reading the driver's CLI.
    """
    method_name: str
    response_form: str
    percentiles: Tuple[float, ...]
    n_bootstrap: int
    seed: int
    include_control: bool
    window_frac: Optional[float]
    secondary: Optional[str]


def _per_percentile_sa(
    method: SAMethod,
    X: np.ndarray,
    Y: np.ndarray,
    factor_names: Sequence[str],
    *,
    n_bootstrap: int,
    seed: int,
    method_kwargs: Dict,
) -> Dict[str, pd.DataFrame]:
    """Run base + bootstrap + rank-stability for a single percentile slice."""
    base = method.compute(X, Y, factor_names, seed=seed, **method_kwargs)
    if n_bootstrap > 0:
        boot = bootstrap_indices(
            method, X, Y, factor_names,
            n_bootstrap=n_bootstrap, seed=seed,
            method_kwargs=method_kwargs,
        )
        rank = bootstrap_rank_stability(
            method, X, Y, factor_names,
            n_bootstrap=n_bootstrap, seed=seed,
            method_kwargs=method_kwargs,
        )
    else:
        boot = pd.DataFrame(
            {"mean": base[method.headline_col].values,
             "ci_lo": np.full(len(factor_names), np.nan),
             "ci_hi": np.full(len(factor_names), np.nan)},
            index=pd.Index(list(factor_names), name="factor"),
        )
        rank = pd.DataFrame(
            {"full_rank": pd.Series(base[method.headline_col].values,
                                    index=factor_names).rank(
                ascending=False, method="average"
             ).values,
             "median_rank": np.full(len(factor_names), np.nan),
             "rank_iqr_lo": np.full(len(factor_names), np.nan),
             "rank_iqr_hi": np.full(len(factor_names), np.nan)},
            index=pd.Index(list(factor_names), name="factor"),
        )
        rank.attrs["rank_spearman_median"] = np.nan
        rank.attrs["rank_spearman_iqr_lo"] = np.nan
        rank.attrs["rank_spearman_iqr_hi"] = np.nan
    return {"base": base, "bootstrap": boot, "rank_stability": rank}


def compute_mv_sa(
    X: np.ndarray,
    M: Optional[np.ndarray],
    factor_names: Sequence[str],
    *,
    method: str,
    response_form: str = "within_trace_percentile",
    percentiles: Sequence[float] = tuple(np.round(np.linspace(0.05, 0.95, 19), 3)),
    secondary: Optional[np.ndarray] = None,
    trace_series: Optional[np.ndarray] = None,
    window_frac: float = 0.30,
    n_bootstrap: int = 200,
    seed: int = 42,
    include_control: bool = True,
    method_kwargs: Optional[Dict] = None,
) -> pd.DataFrame:
    """Run magnitude-varying sensitivity analysis along the M-percentile axis.

    Args:
        X: ``(N, K)`` factor matrix; rows are realizations, columns are
            drought-hazard characteristics.
        M: ``(N,)`` magnitude axis values (e.g., NYC min storage frac).
        factor_names: Names of the K factor columns of X. The control
            factor (if requested) is appended after these.
        method: Method key from :data:`src.sensitivity.sensitivity.METHODS`
            (``"delta"``, ``"pawn"``, ``"rbd_fast"``).
        response_form: ``"exceedance"`` (default) builds a binary
            indicator ``I(tau) = 1{M <= M_quantile(tau)}`` per percentile;
            ``"conditional"`` restricts to a window centered on
            percentile ``tau`` and runs SA on ``secondary``.
        percentiles: Iterable of tau values in (0, 1). Default is
            ``np.linspace(0.05, 0.95, 19)`` (19 points; matches the
            stacked-area visualisation grain in Hadjimichael Fig. 4).
        secondary: Required when ``response_form == "conditional"``.
            Length-N array of the secondary outcome (e.g., FFMP Level 6
            indicator, Montague vulnerability).
        window_frac: Fraction of the sample retained around each tau
            in conditional mode. Default 0.30 keeps 150 of 500
            realizations per slice — large enough for SA to be stable,
            narrow enough to be severity-conditioned.
        n_bootstrap: Bootstrap replicates for CIs and rank stability.
            Set to 0 to disable (test/diagnostic use only).
        seed: RNG seed; pinned for reproducibility.
        include_control: Append a uniform-random control factor whose
            sensitivity defines the empirical noise floor.
        method_kwargs: Extra kwargs forwarded to
            ``method.compute`` / ``bootstrap_*``.

    Returns:
        Long-form DataFrame with columns
        ``percentile, method, factor, headline_index, ci_lo, ci_hi,
        full_rank, median_rank, rank_iqr_lo, rank_iqr_hi,
        rank_spearman_median, n_used, threshold``.

        ``threshold`` is the M-quantile at this percentile (always
        recorded so the figure can show what value of the magnitude
        axis each tau corresponds to). ``n_used`` is the number of
        realizations actually entering the SA at that slice
        (full N for exceedance; window subset for conditional).
    """
    if response_form not in {"exceedance", "conditional",
                             "within_trace_percentile"}:
        raise ValueError(
            f"response_form must be one of {{'exceedance', "
            f"'conditional', 'within_trace_percentile'}}; "
            f"got {response_form!r}"
        )
    if response_form == "conditional" and secondary is None:
        raise ValueError(
            "response_form='conditional' requires a 'secondary' array."
        )
    if response_form == "within_trace_percentile" and trace_series is None:
        raise ValueError(
            "response_form='within_trace_percentile' requires a "
            "'trace_series' (N, T) array."
        )
    if method not in METHODS:
        raise KeyError(
            f"Unknown method {method!r}; known: {sorted(METHODS.keys())}"
        )

    method_kwargs = dict(method_kwargs or {})
    method_obj = resolve_method(method)
    headline_col = HEADLINE_INDEX[method]

    X = np.asarray(X, dtype=float)
    if response_form == "within_trace_percentile":
        # M is unused; trace_series is the source of all responses. Build
        # a row-length proxy for the row count check, drop rows non-finite
        # anywhere in the trace series (any year missing for a trace
        # disqualifies that trace from the SA).
        trace_series = np.asarray(trace_series, dtype=float)
        if trace_series.ndim != 2:
            raise ValueError(
                f"trace_series must be 2-D (N, T); got shape {trace_series.shape}"
            )
        if X.shape[0] != trace_series.shape[0]:
            raise ValueError(
                f"X has {X.shape[0]} rows but trace_series has "
                f"{trace_series.shape[0]}; must match."
            )
        finite = (np.all(np.isfinite(trace_series), axis=1)
                  & np.all(np.isfinite(X), axis=1))
        X = X[finite]
        trace_series = trace_series[finite]
        # M is the trace-min for legacy threshold/n_used reporting only.
        M = trace_series.min(axis=1)
    else:
        M = np.asarray(M, dtype=float)
        if X.shape[0] != M.shape[0]:
            raise ValueError(
                f"X has {X.shape[0]} rows but M has {M.shape[0]}; must match."
            )
        if secondary is not None:
            secondary = np.asarray(secondary, dtype=float)
            if secondary.shape[0] != M.shape[0]:
                raise ValueError(
                    f"secondary has {secondary.shape[0]} rows but M has "
                    f"{M.shape[0]}; must match."
                )
        finite = np.isfinite(M) & np.all(np.isfinite(X), axis=1)
        if secondary is not None:
            finite = finite & np.isfinite(secondary)
        X = X[finite]
        M = M[finite]
        if secondary is not None:
            secondary = secondary[finite]

    if include_control:
        X, factor_names_full = _augment_with_control(
            X, factor_names, seed=seed,
        )
    else:
        factor_names_full = list(factor_names)

    rows = []
    for tau in percentiles:
        tau = float(tau)
        threshold = float(np.quantile(M, tau))
        if response_form == "exceedance":
            Y = _exceedance_indicator(M, tau)
            X_used = X
            n_used = int(X_used.shape[0])
            # Degenerate slices (tau=0 or 1) collapse Y to constant.
            if Y.std() == 0.0:
                # Skip this percentile cleanly.
                for f in factor_names_full:
                    rows.append({
                        "percentile": tau,
                        "method": method,
                        "factor": f,
                        "headline_index": np.nan,
                        "ci_lo": np.nan, "ci_hi": np.nan,
                        "full_rank": np.nan, "median_rank": np.nan,
                        "rank_iqr_lo": np.nan, "rank_iqr_hi": np.nan,
                        "rank_spearman_median": np.nan,
                        "n_used": n_used,
                        "threshold": threshold,
                    })
                continue
        elif response_form == "within_trace_percentile":
            # Y_i(tau) = the tau-th percentile of trace i's own series.
            # Hadjimichael Variant 2 (magnitude-of-shortage) with the
            # trace's 20-year annual sequence playing the role of the
            # SOW's inner ensemble.
            Y = np.percentile(trace_series, 100.0 * tau, axis=1)
            X_used = X
            n_used = int(X_used.shape[0])
            # `threshold` reported here is the cross-trace median of
            # Y(tau): a useful axis annotation (e.g. "the typical 20-th
            # percentile annual storage across traces is X").
            threshold = float(np.median(Y))
            if Y.std() == 0.0:
                for f in factor_names_full:
                    rows.append({
                        "percentile": tau,
                        "method": method,
                        "factor": f,
                        "headline_index": np.nan,
                        "ci_lo": np.nan, "ci_hi": np.nan,
                        "full_rank": np.nan, "median_rank": np.nan,
                        "rank_iqr_lo": np.nan, "rank_iqr_hi": np.nan,
                        "rank_spearman_median": np.nan,
                        "n_used": n_used,
                        "threshold": threshold,
                    })
                continue
        else:  # conditional
            mask = _conditional_window_idx(
                M, tau, window_frac=window_frac,
            )
            X_used = X[mask]
            Y = secondary[mask]
            n_used = int(X_used.shape[0])
            if n_used < 50 or np.unique(Y).size < 2:
                for f in factor_names_full:
                    rows.append({
                        "percentile": tau,
                        "method": method,
                        "factor": f,
                        "headline_index": np.nan,
                        "ci_lo": np.nan, "ci_hi": np.nan,
                        "full_rank": np.nan, "median_rank": np.nan,
                        "rank_iqr_lo": np.nan, "rank_iqr_hi": np.nan,
                        "rank_spearman_median": np.nan,
                        "n_used": n_used,
                        "threshold": threshold,
                    })
                continue

        try:
            slice_out = _per_percentile_sa(
                method_obj, X_used, Y, factor_names_full,
                n_bootstrap=n_bootstrap, seed=seed,
                method_kwargs=method_kwargs,
            )
        except Exception as exc:
            for f in factor_names_full:
                rows.append({
                    "percentile": tau,
                    "method": method,
                    "factor": f,
                    "headline_index": np.nan,
                    "ci_lo": np.nan, "ci_hi": np.nan,
                    "full_rank": np.nan, "median_rank": np.nan,
                    "rank_iqr_lo": np.nan, "rank_iqr_hi": np.nan,
                    "rank_spearman_median": np.nan,
                    "n_used": n_used,
                    "threshold": threshold,
                    "error": str(exc),
                })
            continue

        base = slice_out["base"]
        boot = slice_out["bootstrap"]
        rank = slice_out["rank_stability"]
        rank_rho = float(rank.attrs.get("rank_spearman_median", np.nan))
        for f in factor_names_full:
            rows.append({
                "percentile": tau,
                "method": method,
                "factor": f,
                "headline_index": float(base.loc[f, headline_col]),
                "ci_lo": float(boot.loc[f, "ci_lo"])
                if "ci_lo" in boot.columns else np.nan,
                "ci_hi": float(boot.loc[f, "ci_hi"])
                if "ci_hi" in boot.columns else np.nan,
                "full_rank": float(rank.loc[f, "full_rank"]),
                "median_rank": float(rank.loc[f, "median_rank"]),
                "rank_iqr_lo": float(rank.loc[f, "rank_iqr_lo"]),
                "rank_iqr_hi": float(rank.loc[f, "rank_iqr_hi"]),
                "rank_spearman_median": rank_rho,
                "n_used": n_used,
                "threshold": threshold,
            })

    df = pd.DataFrame(rows)
    df.attrs["mv_sa_config"] = MVSAConfig(
        method_name=method,
        response_form=response_form,
        percentiles=tuple(float(t) for t in percentiles),
        n_bootstrap=int(n_bootstrap),
        seed=int(seed),
        include_control=bool(include_control),
        window_frac=float(window_frac) if response_form == "conditional"
        else None,
        secondary="<provided>" if secondary is not None else None,
    ).__dict__
    return df


# ---------------------------------------------------------------------------
# Post-processing for stacked-area visualisation
# ---------------------------------------------------------------------------


def stacked_share(
    df: pd.DataFrame,
    *,
    headline_col: str = "headline_index",
    factor_order: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Convert per-(percentile, factor) indices into normalized shares.

    For each percentile, the share of factor i is

        share_i(tau) = max(0, idx_i(tau)) / sum_j max(0, idx_j(tau))

    Negative indices (which can occur in finite-sample RBD-FAST) are
    clipped to zero before normalization. Rows whose total is zero
    receive equal shares so the stack still tiles unit area.

    Args:
        df: Long-form output of :func:`compute_mv_sa`.
        headline_col: Column carrying the index magnitude.
        factor_order: Optional explicit factor order for the returned
            wide table. Defaults to alphabetical from ``df``.

    Returns:
        Wide DataFrame indexed by ``percentile`` with one column per
        factor; rows sum to 1.0 (or NaN for percentiles that failed).
    """
    sub = df.copy()
    sub["_clipped"] = sub[headline_col].clip(lower=0.0)
    pivot = sub.pivot_table(
        index="percentile", columns="factor", values="_clipped",
        aggfunc="first",
    )
    if factor_order is not None:
        pivot = pivot.reindex(columns=list(factor_order))
    totals = pivot.sum(axis=1)
    out = pivot.div(totals.where(totals > 0, np.nan), axis=0)
    # When totals are all zero (degenerate slice), assign equal share so
    # the stacked area still tiles. NaN rows propagate when the slice
    # was skipped entirely.
    zero_mask = (totals == 0) & pivot.notna().all(axis=1)
    if zero_mask.any():
        equal = 1.0 / pivot.shape[1]
        out.loc[zero_mask, :] = equal
    return out
