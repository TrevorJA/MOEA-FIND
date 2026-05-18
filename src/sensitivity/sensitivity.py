"""Global sensitivity analysis on MOEA-FIND drought-hazard realizations.

Wrappers around :mod:`SALib` for the three sample-free SA methods compatible
with a fixed structured realization sample: Delta moment-independent
(Borgonovo, 2007), PAWN density-based (Pianosi & Wagener, 2015), and
RBD-FAST (Tarantola et al., 2006). Sobol / Morris / FAST require dedicated
Saltelli / radial / frequency-coded designs and are excluded by Stage-08
scope.

The module is intentionally narrow: each method has one entry point that
takes a paired ``(X, Y)`` and returns a tidy DataFrame indexed by factor.
Bootstrap CIs, rank-stability, cross-method rank correlation, and
sample-size convergence helpers all consume the same per-method index
function so additional methods can be slotted in later without touching
the diagnostic layer.

Public surface (consumed by ``workflows/08_nyc_sensitivity/run_sa.py``):

- :func:`compute_delta`
- :func:`compute_pawn`
- :func:`compute_rbd_fast`
- :func:`bootstrap_indices`
- :func:`bootstrap_rank_stability`
- :func:`cross_method_rank_corr`
- :func:`convergence_curve`
- :func:`drop_low_variance_factors`
- :func:`apply_method_selection_criterion`

Notes
-----
The Manhattan-norm auxiliary objective (:math:`f_{K+1}` in DD-11) is
algorithmic, not a hazard attribute, and must be excluded by the caller
before the X array reaches this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# Per-method headline-index column. The value of this column is what the
# diagnostic layer (rank stability, cross-method correlation, convergence,
# selection criterion) consumes; anything else in the per-method DataFrame
# is decoration. Update only when SALib's output schema changes.
HEADLINE_INDEX = {
    "delta": "delta",
    "pawn": "median",
    "rbd_fast": "S1",
}


# ---------------------------------------------------------------------------
# SALib problem builder
# ---------------------------------------------------------------------------


def _build_problem(factor_names: Sequence[str], X: np.ndarray) -> Dict:
    """Build the SALib ``problem`` dict from the empirical X matrix.

    SALib only uses ``num_vars`` and ``names`` for the analyze step;
    ``bounds`` are required by the schema but unused by Delta / PAWN /
    RBD-FAST when X is supplied directly. We pass empirical min/max so
    that any debug-print is interpretable.
    """
    if X.ndim != 2:
        raise ValueError(f"X must be 2-D; got shape {X.shape}")
    if X.shape[1] != len(factor_names):
        raise ValueError(
            f"X has {X.shape[1]} factor columns but {len(factor_names)} "
            f"factor_names were supplied: {list(factor_names)}"
        )
    bounds = [[float(X[:, i].min()), float(X[:, i].max())]
              for i in range(X.shape[1])]
    return {"num_vars": len(factor_names), "names": list(factor_names),
            "bounds": bounds}


def _validate_xy(X: np.ndarray, Y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Drop rows with non-finite X or Y; raise if too few rows remain."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    if Y.ndim != 1:
        raise ValueError(f"Y must be 1-D; got shape {Y.shape}")
    if X.shape[0] != Y.shape[0]:
        raise ValueError(
            f"X and Y must have matching row counts; got {X.shape[0]} vs "
            f"{Y.shape[0]}"
        )
    finite = np.isfinite(Y) & np.all(np.isfinite(X), axis=1)
    if not finite.any():
        raise ValueError("All rows have non-finite X or Y after dropping NaN.")
    return X[finite], Y[finite]


# ---------------------------------------------------------------------------
# Per-method analyze wrappers
# ---------------------------------------------------------------------------


def compute_delta(
    X: np.ndarray,
    Y: np.ndarray,
    factor_names: Sequence[str],
    *,
    num_resamples: int = 100,
    conf_level: float = 0.95,
    seed: int = 42,
) -> pd.DataFrame:
    """Delta moment-independent SA via :mod:`SALib.analyze.delta`.

    Returns a DataFrame indexed by factor with columns
    ``delta``, ``delta_conf`` (half-width of the conf_level CI),
    ``S1`` (Plischke-estimator first-order Sobol), ``S1_conf``,
    ``ci_lo``, ``ci_hi`` (Delta CI bounds).
    """
    from SALib.analyze import delta as sa_delta

    X, Y = _validate_xy(X, Y)
    problem = _build_problem(factor_names, X)
    result = sa_delta.analyze(
        problem, X, Y,
        num_resamples=num_resamples,
        conf_level=conf_level,
        seed=seed,
    )
    df = pd.DataFrame({
        "delta": np.asarray(result["delta"], dtype=float),
        "delta_conf": np.asarray(result["delta_conf"], dtype=float),
        "S1": np.asarray(result["S1"], dtype=float),
        "S1_conf": np.asarray(result["S1_conf"], dtype=float),
    }, index=pd.Index(list(factor_names), name="factor"))
    df["ci_lo"] = df["delta"] - df["delta_conf"]
    df["ci_hi"] = df["delta"] + df["delta_conf"]
    return df


def compute_pawn(
    X: np.ndarray,
    Y: np.ndarray,
    factor_names: Sequence[str],
    *,
    S: int = 10,
    seed: int = 42,
    n_bootstrap: int = 1000,
    conf_level: float = 0.95,
) -> pd.DataFrame:
    """PAWN density-based SA via :mod:`SALib.analyze.pawn`.

    SALib's PAWN does not return CIs; the wrapper computes a percentile
    bootstrap CI on the median (the headline index) by resampling rows
    with replacement.

    Returns a DataFrame indexed by factor with columns ``minimum``,
    ``mean``, ``median``, ``maximum``, ``CV``, ``stdev``, ``ci_lo``,
    ``ci_hi`` (CIs on ``median``).
    """
    from SALib.analyze import pawn as sa_pawn

    X, Y = _validate_xy(X, Y)
    problem = _build_problem(factor_names, X)

    point = sa_pawn.analyze(problem, X, Y, S=S, seed=seed)
    df = pd.DataFrame({
        "minimum": np.asarray(point["minimum"], dtype=float),
        "mean": np.asarray(point["mean"], dtype=float),
        "median": np.asarray(point["median"], dtype=float),
        "maximum": np.asarray(point["maximum"], dtype=float),
        "CV": np.asarray(point["CV"], dtype=float),
        "stdev": np.asarray(point["stdev"], dtype=float),
    }, index=pd.Index(list(factor_names), name="factor"))

    if n_bootstrap > 0:
        rng = np.random.default_rng(seed)
        n_rows = X.shape[0]
        boots = np.empty((n_bootstrap, len(factor_names)), dtype=float)
        for b in range(n_bootstrap):
            idx = rng.integers(0, n_rows, size=n_rows)
            try:
                rb = sa_pawn.analyze(
                    problem, X[idx], Y[idx], S=S, seed=int(rng.integers(2**31)),
                )
                boots[b] = np.asarray(rb["median"], dtype=float)
            except Exception:
                boots[b] = np.nan
        alpha = (1.0 - conf_level) / 2.0
        df["ci_lo"] = np.nanpercentile(boots, 100.0 * alpha, axis=0)
        df["ci_hi"] = np.nanpercentile(boots, 100.0 * (1.0 - alpha), axis=0)
    else:
        df["ci_lo"] = np.nan
        df["ci_hi"] = np.nan
    return df


def compute_rbd_fast(
    X: np.ndarray,
    Y: np.ndarray,
    factor_names: Sequence[str],
    *,
    M: int = 10,
    num_resamples: int = 100,
    conf_level: float = 0.95,
    seed: int = 42,
) -> pd.DataFrame:
    """RBD-FAST first-order SA via :mod:`SALib.analyze.rbd_fast`.

    SALib's RBD-FAST returns ``S1`` and ``S1_conf`` (half-width). Returns
    a DataFrame indexed by factor with columns ``S1``, ``S1_conf``,
    ``ci_lo``, ``ci_hi``.
    """
    from SALib.analyze import rbd_fast as sa_rbd

    X, Y = _validate_xy(X, Y)
    problem = _build_problem(factor_names, X)

    result = sa_rbd.analyze(
        problem, X, Y,
        M=M,
        num_resamples=num_resamples,
        conf_level=conf_level,
        seed=seed,
    )
    s1 = np.asarray(result["S1"], dtype=float)
    s1_conf = np.asarray(result.get("S1_conf", np.full_like(s1, np.nan)),
                         dtype=float)
    df = pd.DataFrame({
        "S1": s1,
        "S1_conf": s1_conf,
        "ci_lo": s1 - s1_conf,
        "ci_hi": s1 + s1_conf,
    }, index=pd.Index(list(factor_names), name="factor"))
    return df


# ---------------------------------------------------------------------------
# Method registry — keeps the diagnostic layer method-agnostic
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SAMethod:
    """Metadata for a sensitivity-analysis method.

    The ``compute`` callable has signature
    ``compute(X, Y, factor_names, **kwargs) -> pd.DataFrame`` indexed by
    factor with at least ``ci_lo`` / ``ci_hi`` and ``HEADLINE_INDEX[name]``.
    """
    name: str
    label: str
    headline_col: str
    compute: Callable[..., pd.DataFrame]


METHODS: Dict[str, SAMethod] = {
    "delta": SAMethod(
        name="delta", label="Delta moment-independent",
        headline_col=HEADLINE_INDEX["delta"], compute=compute_delta,
    ),
    "pawn": SAMethod(
        name="pawn", label="PAWN density-based",
        headline_col=HEADLINE_INDEX["pawn"], compute=compute_pawn,
    ),
    "rbd_fast": SAMethod(
        name="rbd_fast", label="RBD-FAST",
        headline_col=HEADLINE_INDEX["rbd_fast"], compute=compute_rbd_fast,
    ),
}


def resolve_method(name: str) -> SAMethod:
    """Look up a registered method by name. Raises KeyError on miss."""
    if name not in METHODS:
        raise KeyError(
            f"Unknown SA method {name!r}. Known: {sorted(METHODS.keys())}"
        )
    return METHODS[name]


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------


def drop_low_variance_factors(
    X: pd.DataFrame,
    cv_threshold: float = 0.05,
) -> Tuple[pd.DataFrame, list]:
    """Drop factor columns whose coefficient of variation is below threshold.

    A near-constant factor produces meaningless SA indices; the upstream
    optimization may have collapsed an objective onto a single value.

    Args:
        X: Per-realization factor DataFrame (rows = realizations, columns
            = factor names).
        cv_threshold: Minimum |std / mean| required to keep a factor. The
            CV is computed in absolute value to handle negative-mean
            metrics (e.g., negated Q10).

    Returns:
        Tuple of ``(X_kept, dropped_names)``. ``dropped_names`` is empty
        if every factor cleared the threshold.
    """
    means = X.mean(axis=0)
    stds = X.std(axis=0)
    cv = stds.abs() / np.where(means.abs() > 1e-12, means.abs(), 1.0)
    keep = cv >= cv_threshold
    dropped = [c for c, k in zip(X.columns, keep) if not k]
    return X.loc[:, keep], dropped


def bootstrap_indices(
    method: SAMethod,
    X: np.ndarray,
    Y: np.ndarray,
    factor_names: Sequence[str],
    *,
    n_bootstrap: int = 1000,
    seed: int = 42,
    conf_level: float = 0.95,
    method_kwargs: Optional[Dict] = None,
) -> pd.DataFrame:
    """Bootstrap the headline index of *method* by resampling rows.

    Returns a DataFrame indexed by factor with columns ``mean``,
    ``ci_lo``, ``ci_hi``, and the per-replicate matrix dropped (use
    :func:`bootstrap_rank_stability` for rank-level diagnostics).

    Methods whose own ``compute`` already bootstraps internally (e.g.,
    ``compute_pawn`` with ``n_bootstrap > 0``) are still re-bootstrapped
    here so the diagnostic layer has a uniform CI definition across
    methods.
    """
    method_kwargs = dict(method_kwargs or {})
    rng = np.random.default_rng(seed)
    n_rows = X.shape[0]
    boots = np.empty((n_bootstrap, len(factor_names)), dtype=float)
    for b in range(n_bootstrap):
        idx = rng.integers(0, n_rows, size=n_rows)
        sub_seed = int(rng.integers(2**31))
        kw = {**method_kwargs, "seed": sub_seed}
        try:
            df_b = method.compute(X[idx], Y[idx], factor_names, **kw)
            boots[b] = df_b[method.headline_col].values
        except Exception:
            boots[b] = np.nan
    alpha = (1.0 - conf_level) / 2.0
    out = pd.DataFrame({
        "mean": np.nanmean(boots, axis=0),
        "ci_lo": np.nanpercentile(boots, 100.0 * alpha, axis=0),
        "ci_hi": np.nanpercentile(boots, 100.0 * (1.0 - alpha), axis=0),
    }, index=pd.Index(list(factor_names), name="factor"))
    return out


def bootstrap_rank_stability(
    method: SAMethod,
    X: np.ndarray,
    Y: np.ndarray,
    factor_names: Sequence[str],
    *,
    n_bootstrap: int = 1000,
    seed: int = 42,
    method_kwargs: Optional[Dict] = None,
) -> pd.DataFrame:
    """Bootstrap-rank Spearman of factor rankings.

    For each bootstrap replicate, compute the headline index per factor,
    rank factors by index magnitude, and Spearman-correlate the ranking
    against the full-sample ranking. Returns a per-factor median rank +
    rank-IQR plus a top-line ``rank_spearman`` distribution summary.

    Returns:
        DataFrame indexed by factor with columns ``full_rank``,
        ``median_rank``, ``rank_iqr_lo`` (25th percentile),
        ``rank_iqr_hi`` (75th percentile). Plus a single-row "summary"
        index containing the median + IQR of the per-replicate Spearman
        correlation against the full-sample ranking.
    """
    from scipy.stats import spearmanr

    method_kwargs = dict(method_kwargs or {})
    full_df = method.compute(X, Y, factor_names, **method_kwargs)
    full_idx = full_df[method.headline_col].values
    # rankdata-style: largest index → rank 1
    full_rank = pd.Series(full_idx, index=factor_names).rank(
        ascending=False, method="average"
    )

    rng = np.random.default_rng(seed)
    n_rows = X.shape[0]
    rank_replicates = np.empty((n_bootstrap, len(factor_names)), dtype=float)
    spearman_replicates = np.empty(n_bootstrap, dtype=float)
    for b in range(n_bootstrap):
        idx = rng.integers(0, n_rows, size=n_rows)
        sub_seed = int(rng.integers(2**31))
        kw = {**method_kwargs, "seed": sub_seed}
        try:
            df_b = method.compute(X[idx], Y[idx], factor_names, **kw)
            ranks = pd.Series(df_b[method.headline_col].values,
                              index=factor_names).rank(
                ascending=False, method="average"
            )
            rank_replicates[b] = ranks.values
            rho, _ = spearmanr(full_rank.values, ranks.values)
            spearman_replicates[b] = rho if np.isfinite(rho) else np.nan
        except Exception:
            rank_replicates[b] = np.nan
            spearman_replicates[b] = np.nan

    out = pd.DataFrame({
        "full_rank": full_rank.values,
        "median_rank": np.nanmedian(rank_replicates, axis=0),
        "rank_iqr_lo": np.nanpercentile(rank_replicates, 25, axis=0),
        "rank_iqr_hi": np.nanpercentile(rank_replicates, 75, axis=0),
    }, index=pd.Index(list(factor_names), name="factor"))
    # Append a one-row summary of the rank-Spearman distribution. This
    # is the diagnostic the selection criterion consumes.
    out.attrs["rank_spearman_median"] = float(
        np.nanmedian(spearman_replicates)
    )
    out.attrs["rank_spearman_iqr_lo"] = float(
        np.nanpercentile(spearman_replicates, 25)
    )
    out.attrs["rank_spearman_iqr_hi"] = float(
        np.nanpercentile(spearman_replicates, 75)
    )
    return out


def cross_method_rank_corr(
    method_results: Mapping[str, pd.DataFrame],
) -> pd.DataFrame:
    """Spearman rank correlation between methods.

    Args:
        method_results: ``{method_name: per-method DataFrame from compute_*}``.
            Each DataFrame must include the method's headline-index
            column (looked up via :data:`HEADLINE_INDEX`).

    Returns:
        Square DataFrame (method × method) of Spearman rho. The diagonal
        is 1.0 by construction.
    """
    from scipy.stats import spearmanr

    methods = list(method_results.keys())
    n = len(methods)
    out = np.full((n, n), np.nan, dtype=float)
    for i, mi in enumerate(methods):
        ci = HEADLINE_INDEX[mi]
        ri = method_results[mi][ci].values
        for j, mj in enumerate(methods):
            cj = HEADLINE_INDEX[mj]
            rj = method_results[mj][cj].values
            if i == j:
                out[i, j] = 1.0
                continue
            rho, _ = spearmanr(ri, rj)
            out[i, j] = rho if np.isfinite(rho) else np.nan
    return pd.DataFrame(out, index=methods, columns=methods)


def cross_outcome_rank_corr(
    per_outcome_indices: Mapping[str, pd.DataFrame],
    method: str,
) -> pd.DataFrame:
    """Spearman rank correlation of factor rankings between outcomes.

    Sensitivity-of-sensitivity diagnostic: for one method, do different
    NYC outcomes attribute variance to the same drought characteristics?

    Args:
        per_outcome_indices: ``{outcome_name: per-method DataFrame}``.
            Each DataFrame must include the method's headline column.
        method: Method name (must be a key of :data:`HEADLINE_INDEX`).

    Returns:
        Square DataFrame (outcome × outcome) of Spearman rho.
    """
    from scipy.stats import spearmanr

    headline = HEADLINE_INDEX[method]
    outcomes = list(per_outcome_indices.keys())
    n = len(outcomes)
    out = np.full((n, n), np.nan, dtype=float)
    for i, oi in enumerate(outcomes):
        ri = per_outcome_indices[oi][headline].values
        for j, oj in enumerate(outcomes):
            rj = per_outcome_indices[oj][headline].values
            if i == j:
                out[i, j] = 1.0
                continue
            rho, _ = spearmanr(ri, rj)
            out[i, j] = rho if np.isfinite(rho) else np.nan
    return pd.DataFrame(out, index=outcomes, columns=outcomes)


def convergence_curve(
    method: SAMethod,
    X: np.ndarray,
    Y: np.ndarray,
    factor_names: Sequence[str],
    *,
    sizes: Sequence[int],
    n_replicates: int = 200,
    seed: int = 42,
    method_kwargs: Optional[Dict] = None,
) -> pd.DataFrame:
    """Sample-size convergence curve for *method*'s headline index.

    For each ``n`` in ``sizes`` (and the full sample), draw
    ``n_replicates`` random subsets of size ``n`` without replacement,
    compute the index, and report mean + 5th/95th percentile per factor.

    Args:
        method: SA method.
        X, Y: Full sample.
        factor_names: Factor names.
        sizes: Subset sizes to scan. The full sample size is appended
            automatically and is computed only once (no resampling).
        n_replicates: Sub-sample replicates per ``n``.
        seed: RNG seed.
        method_kwargs: Extra kwargs for ``method.compute``.

    Returns:
        Long-form DataFrame with columns ``factor``, ``n``, ``mean``,
        ``p05``, ``p95``, ``replicates``.
    """
    method_kwargs = dict(method_kwargs or {})
    rng = np.random.default_rng(seed)
    full_n = X.shape[0]
    sizes = sorted({int(n) for n in sizes if 0 < int(n) < full_n}) + [full_n]

    rows = []
    for n in sizes:
        if n == full_n:
            df_full = method.compute(X, Y, factor_names, **method_kwargs)
            vals = df_full[method.headline_col].values
            for f, v in zip(factor_names, vals):
                rows.append({
                    "factor": f, "n": n,
                    "mean": float(v), "p05": float(v), "p95": float(v),
                    "replicates": 1,
                })
            continue
        replicates = np.empty((n_replicates, len(factor_names)), dtype=float)
        for b in range(n_replicates):
            idx = rng.choice(full_n, size=n, replace=False)
            sub_seed = int(rng.integers(2**31))
            kw = {**method_kwargs, "seed": sub_seed}
            try:
                df_b = method.compute(X[idx], Y[idx], factor_names, **kw)
                replicates[b] = df_b[method.headline_col].values
            except Exception:
                replicates[b] = np.nan
        means = np.nanmean(replicates, axis=0)
        p05 = np.nanpercentile(replicates, 5, axis=0)
        p95 = np.nanpercentile(replicates, 95, axis=0)
        for f, m, lo, hi in zip(factor_names, means, p05, p95):
            rows.append({
                "factor": f, "n": n, "mean": float(m),
                "p05": float(lo), "p95": float(hi),
                "replicates": n_replicates,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Method-selection criterion (per outcome)
# ---------------------------------------------------------------------------


def apply_method_selection_criterion(
    *,
    bootstrap_df: pd.DataFrame,
    rank_stability_df: pd.DataFrame,
    cross_method_rho: pd.Series,
    rank_spearman_threshold: float = 0.8,
    cross_method_threshold: float = 0.7,
) -> Dict[str, object]:
    """Evaluate the three-condition method-selection criterion.

    The criterion (fixed before HPC, see plan §"Diagnostics + manuscript-
    method selection criterion") requires *simultaneously*:

    (a) Bootstrap CI not spanning zero for the top-ranked factor.
    (b) Bootstrap-rank Spearman ≥ ``rank_spearman_threshold`` (median).
    (c) Cross-method rank correlation ≥ ``cross_method_threshold`` with
        at least one comparator on the same outcome.

    Args:
        bootstrap_df: Output of :func:`bootstrap_indices` (per factor;
            must contain ``mean``, ``ci_lo``, ``ci_hi``).
        rank_stability_df: Output of :func:`bootstrap_rank_stability`
            (must contain ``full_rank`` column and the
            ``rank_spearman_median`` attr).
        cross_method_rho: Series indexed by *other-method-name* with
            Spearman rho of this method against each comparator.
        rank_spearman_threshold: Threshold for condition (b).
        cross_method_threshold: Threshold for condition (c).

    Returns:
        Dict with keys ``passes`` (bool), ``top_factor`` (str),
        ``top_factor_ci_excludes_zero`` (bool), ``rank_spearman_median``
        (float), ``rank_spearman_passes`` (bool),
        ``best_cross_method_rho`` (float), ``cross_method_passes`` (bool),
        and ``failing_conditions`` (list of strings).
    """
    full_rank = rank_stability_df["full_rank"]
    top_factor = full_rank.idxmin()  # rank 1 = top
    top_lo = float(bootstrap_df.loc[top_factor, "ci_lo"])
    top_hi = float(bootstrap_df.loc[top_factor, "ci_hi"])
    ci_excludes_zero = (top_lo > 0.0) or (top_hi < 0.0)

    rank_rho = float(rank_stability_df.attrs.get(
        "rank_spearman_median", np.nan
    ))
    rank_passes = np.isfinite(rank_rho) and rank_rho >= rank_spearman_threshold

    rho_vals = cross_method_rho.dropna()
    best_rho = float(rho_vals.max()) if len(rho_vals) else np.nan
    cross_passes = (
        np.isfinite(best_rho) and best_rho >= cross_method_threshold
    )

    failing = []
    if not ci_excludes_zero:
        failing.append("ci_spans_zero")
    if not rank_passes:
        failing.append("rank_unstable")
    if not cross_passes:
        failing.append("cross_method_disagreement")

    return {
        "passes": (ci_excludes_zero and rank_passes and cross_passes),
        "top_factor": str(top_factor),
        "top_factor_ci_excludes_zero": bool(ci_excludes_zero),
        "rank_spearman_median": rank_rho,
        "rank_spearman_passes": bool(rank_passes),
        "best_cross_method_rho": best_rho,
        "cross_method_passes": bool(cross_passes),
        "failing_conditions": failing,
    }
