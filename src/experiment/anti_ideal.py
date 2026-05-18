"""Anti-ideal placement and epsilon construction for MOEA-FIND.

Split out of the former ``src.experiment_utils`` god module: pure
objective-space helpers, no data IO and no MOEA orchestration, so they
can be reused without pulling in the experiment harness.
"""

from typing import Dict, Optional

import numpy as np


def compute_ssi_anti_ideal(
    hist_chars: dict,
    objective_keys,
    headroom: float = 1.5,
    feasible_maxes: Optional[Dict[str, float]] = None,
) -> np.ndarray:
    """Build the DD-11 anti-ideal point ``D*`` from historical drought chars.

    Thin shim that resolves ``objective_keys`` to a tuple of
    :class:`src.metrics.drought_metrics.DroughtMetric` and delegates to
    :func:`src.metrics.drought_metrics.compute_anti_ideal`. Each metric's
    :class:`AntiIdealRule` decides placement: ``HEADROOM_TIMES_MAX`` for
    unbounded-above non-cyclic metrics, ``CYCLIC_HEADROOM`` for cyclic
    calendar metrics (``12 × headroom``), and ``CONSTANT`` for metrics
    with a natural upper bound (e.g. fractions in ``[0, 1]``).

    DD-11 requires ``D_j(x) <= D*_j`` for every feasible ``x``; the
    rules above guarantee this provided ``headroom > 1``.

    Args:
        hist_chars: Historical drought characteristics dict.
        objective_keys: Either a tuple of metric names or a tuple of
            :class:`DroughtMetric` instances.
        headroom: Safety factor applied to the historical maximum.
        feasible_maxes: Optional ``{metric_name: observed_max}`` override
            for ``HEADROOM_TIMES_MAX`` metrics. Cyclic and constant
            metrics ignore this argument.

    Returns:
        Anti-ideal array ``D*`` of length ``len(objective_keys)``.
    """
    from src.metrics.drought_metrics import compute_anti_ideal, resolve_metric_set

    metric_set = resolve_metric_set(objective_keys)
    return compute_anti_ideal(
        metric_set,
        hist_chars,
        headroom=headroom,
        feasible_maxes=feasible_maxes,
    )


def build_epsilons(
    objective_keys,
    epsilon_map: Optional[Dict[str, float]] = None,
    manhattan_eps: Optional[float] = None,
) -> list:
    """Return epsilon list for objectives plus the Manhattan-norm auxiliary.

    Args:
        objective_keys: Either a tuple of metric names or a tuple of
            :class:`src.metrics.drought_metrics.DroughtMetric` instances. When
            metric instances are passed, the per-axis epsilon is read
            from each metric's ``epsilon`` field. When string names are
            passed, ``epsilon_map`` is consulted (falls back to ``0.5``).
        epsilon_map: Optional override mapping metric name → epsilon.
            Used only when ``objective_keys`` is a tuple of strings.
        manhattan_eps: Epsilon for the ``f_{K+1}`` auxiliary objective.
            Defaults to
            :data:`src.experiment.config.DEFAULT_EXPERIMENT.manhattan_eps`.
    """
    from src.metrics.drought_metrics import DroughtMetric, resolve_metric_set

    keys = tuple(objective_keys)
    if len(keys) > 0 and isinstance(keys[0], DroughtMetric):
        eps = [m.epsilon for m in keys]
    elif epsilon_map is not None:
        eps = [epsilon_map.get(k, 0.5) for k in keys]
    else:
        # Resolve string names to metrics so each carries its own epsilon.
        metric_set = resolve_metric_set(keys)
        eps = [m.epsilon for m in metric_set]

    if manhattan_eps is None:
        from src.experiment.config import DEFAULT_EXPERIMENT
        manhattan_eps = DEFAULT_EXPERIMENT.manhattan_eps
    eps.append(float(manhattan_eps))
    return eps
