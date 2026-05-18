"""Unit tests for src.sensitivity.magnitude_varying_sa.

Constructs a synthetic regime where the dominant factor *changes*
across percentiles of the magnitude axis, then checks that the MV-SA
engine recovers that severity-dependent ranking.

Synthetic design::

    X1, X2 ~ U[0, 1] (independent)
    M = X1 if X1 > 0.6 else 0.5 * X2

The two regimes produce non-overlapping M ranges: 60% of realizations
have M = 0.5*X2 in [0, 0.5] and 40% have M = X1 in (0.6, 1]. The
indicator I(tau) = 1{M <= q(tau)} therefore decomposes cleanly:

  - At low tau (e.g. tau=0.20, q approx 0.167): only the X2-regime
    can contribute I=1; X1-regime is uniformly above the threshold so
    I=0 there. Sensitivity is overwhelmingly to X2.
  - At high tau (e.g. tau=0.85, q approx 0.85): every X2-regime row
    contributes I=1 deterministically; the X1-regime contributes
    I = 1{X1 <= 0.85}. Sensitivity is overwhelmingly to X1.

A uniform-random control factor must sit at the noise floor at every
percentile.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("SALib")

from src.sensitivity.magnitude_varying_sa import (  # noqa: E402
    CONTROL_FACTOR_NAME,
    compute_mv_sa,
    stacked_share,
)


def _synthetic_split_regime(n: int = 800, seed: int = 0):
    """Two-factor magnitude-varying regime with a clean low/high switch."""
    rng = np.random.default_rng(seed)
    X1 = rng.uniform(0.0, 1.0, size=n)
    X2 = rng.uniform(0.0, 1.0, size=n)
    M = np.where(X1 > 0.6, X1, 0.5 * X2)
    X = np.column_stack([X1, X2])
    return X, M, ["X1", "X2"]


def _factor_at_quantile(df: pd.DataFrame, tau: float, factor: str) -> float:
    """Headline-index value for one factor at one percentile."""
    sub = df[(df["percentile"] == tau) & (df["factor"] == factor)]
    if sub.empty:
        return float("nan")
    return float(sub["headline_index"].iloc[0])


# ---------------------------------------------------------------------------
# Severity-dependent ranking on the synthetic regime
# ---------------------------------------------------------------------------


class TestSeverityDependentRanking:
    """Engine recovers the magnitude-dependent dominant factor."""

    def test_low_percentile_x2_dominates(self):
        X, M, names = _synthetic_split_regime(n=800, seed=0)
        df = compute_mv_sa(
            X, M, names,
            method="delta",
            percentiles=(0.10, 0.20),
            n_bootstrap=0,  # base indices only — fast
            include_control=False,
            method_kwargs={"num_resamples": 20},
        )
        for tau in (0.10, 0.20):
            x1 = _factor_at_quantile(df, tau, "X1")
            x2 = _factor_at_quantile(df, tau, "X2")
            assert np.isfinite(x1) and np.isfinite(x2)
            assert x2 > x1, (
                f"At low tau={tau} expected X2 > X1; got X1={x1}, X2={x2}"
            )

    def test_high_percentile_x1_dominates(self):
        X, M, names = _synthetic_split_regime(n=800, seed=0)
        df = compute_mv_sa(
            X, M, names,
            method="delta",
            percentiles=(0.80, 0.90),
            n_bootstrap=0,
            include_control=False,
            method_kwargs={"num_resamples": 20},
        )
        for tau in (0.80, 0.90):
            x1 = _factor_at_quantile(df, tau, "X1")
            x2 = _factor_at_quantile(df, tau, "X2")
            assert np.isfinite(x1) and np.isfinite(x2)
            assert x1 > x2, (
                f"At high tau={tau} expected X1 > X2; got X1={x1}, X2={x2}"
            )


# ---------------------------------------------------------------------------
# Control factor must stay at the noise floor everywhere
# ---------------------------------------------------------------------------


def test_control_factor_below_real_factors():
    X, M, names = _synthetic_split_regime(n=800, seed=1)
    df = compute_mv_sa(
        X, M, names,
        method="delta",
        percentiles=(0.10, 0.30, 0.50, 0.70, 0.90),
        n_bootstrap=0,
        include_control=True,
        method_kwargs={"num_resamples": 20},
    )
    assert CONTROL_FACTOR_NAME in df["factor"].unique()
    # The dominant real factor at each tau should beat the control.
    for tau in (0.10, 0.30, 0.50, 0.70, 0.90):
        ctrl = _factor_at_quantile(df, tau, CONTROL_FACTOR_NAME)
        real_max = max(
            _factor_at_quantile(df, tau, f) for f in ("X1", "X2")
        )
        assert real_max > ctrl, (
            f"At tau={tau} dominant real factor ({real_max}) does not "
            f"exceed control ({ctrl})"
        )


# ---------------------------------------------------------------------------
# Long-form schema and stacked-share post-processing
# ---------------------------------------------------------------------------


def test_long_form_schema():
    X, M, names = _synthetic_split_regime(n=400, seed=2)
    df = compute_mv_sa(
        X, M, names,
        method="delta",
        percentiles=(0.25, 0.50, 0.75),
        n_bootstrap=10,
        include_control=True,
        method_kwargs={"num_resamples": 10},
    )
    expected = {
        "percentile", "method", "factor", "headline_index",
        "ci_lo", "ci_hi", "full_rank", "median_rank",
        "rank_iqr_lo", "rank_iqr_hi", "rank_spearman_median",
        "n_used", "threshold",
    }
    assert expected.issubset(set(df.columns))
    # Each (percentile, factor) pair appears exactly once.
    counts = df.groupby(["percentile", "factor"]).size()
    assert (counts == 1).all()
    # 3 percentiles x 3 factors (X1, X2, control) = 9 rows.
    assert len(df) == 9


def test_stacked_share_rows_sum_to_one():
    X, M, names = _synthetic_split_regime(n=400, seed=3)
    df = compute_mv_sa(
        X, M, names,
        method="delta",
        percentiles=(0.20, 0.50, 0.80),
        n_bootstrap=0,
        include_control=False,
        method_kwargs={"num_resamples": 10},
    )
    shares = stacked_share(df)
    # Each non-degenerate row sums to 1.
    row_sums = shares.sum(axis=1, skipna=True)
    np.testing.assert_allclose(row_sums.values, 1.0, atol=1e-9)


# ---------------------------------------------------------------------------
# Conditional response form runs end-to-end
# ---------------------------------------------------------------------------


def test_conditional_response_form_runs():
    X, M, names = _synthetic_split_regime(n=600, seed=4)
    # Construct a secondary outcome that depends mainly on X2 in the
    # M-tail (low M = X1<=0.5 regime). At tau=0.20 we expect X2 to
    # dominate the conditional SA; we only test that the call runs and
    # returns a sensible-shape table.
    rng = np.random.default_rng(5)
    secondary = X[:, 1] + 0.05 * rng.normal(size=X.shape[0])
    df = compute_mv_sa(
        X, M, names,
        method="delta",
        response_form="conditional",
        percentiles=(0.20, 0.50, 0.80),
        secondary=secondary,
        window_frac=0.40,
        n_bootstrap=0,
        include_control=True,
        method_kwargs={"num_resamples": 10},
    )
    assert df["n_used"].max() > 100  # at least the window subset
    assert df["n_used"].max() <= 600
    # Each percentile reports the same n_used (deterministic window).
    n_per_tau = df.groupby("percentile")["n_used"].nunique()
    assert (n_per_tau == 1).all()


# ---------------------------------------------------------------------------
# Bootstrap CIs are returned ordered when requested
# ---------------------------------------------------------------------------


def test_bootstrap_ci_ordered():
    X, M, names = _synthetic_split_regime(n=400, seed=6)
    df = compute_mv_sa(
        X, M, names,
        method="delta",
        percentiles=(0.30, 0.70),
        n_bootstrap=20,
        include_control=False,
        method_kwargs={"num_resamples": 10},
    )
    finite = df.dropna(subset=["ci_lo", "ci_hi"])
    assert len(finite) > 0
    assert (finite["ci_lo"] <= finite["ci_hi"] + 1e-12).all()


# ---------------------------------------------------------------------------
# Degenerate slice (constant indicator) returns NaN, not a crash
# ---------------------------------------------------------------------------


def test_degenerate_percentile_returns_nan():
    rng = np.random.default_rng(7)
    n = 200
    X = np.column_stack([rng.uniform(0, 1, n), rng.uniform(0, 1, n)])
    M = np.full(n, 0.5)  # constant magnitude axis -> indicator is constant
    df = compute_mv_sa(
        X, M, ["X1", "X2"],
        method="delta",
        percentiles=(0.25, 0.75),
        n_bootstrap=0,
        include_control=False,
        method_kwargs={"num_resamples": 5},
    )
    # All headline_index entries should be NaN (degenerate slice).
    assert df["headline_index"].isna().all()
