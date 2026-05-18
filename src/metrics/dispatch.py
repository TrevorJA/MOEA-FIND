"""Metric-family dispatcher: which metric family applies at a given T.

The T -> family decision was duplicated across ``run_moea_find.py`` and
``src.experiment.runner.run_experiment`` (short-block raw-flow path for
very short traces vs the SSI drought-characteristic path otherwise).
This module centralises that single decision so future call sites route
through one place instead of re-deriving the threshold.

It is a thin lazy-import router — importing it pulls in nothing heavy,
and it changes no existing behavior on its own.
"""

from typing import Literal, Tuple

import numpy as np

# Traces of <= this many water years use the sub-annual raw-flow
# short-block family (DD-15); longer traces use the SSI path.
SHORT_BLOCK_MAX_YEARS = 2

Family = Literal["short_block", "extended", "registry"]


def resolve_metric_family(
    timescale: int,
    *,
    mode: str = "objectives",
) -> Tuple[Family, object]:
    """Return ``(family_name, handle)`` for an SSI timescale ``T``.

    * ``T <= SHORT_BLOCK_MAX_YEARS``            -> ``("short_block", module)``
    * ``mode == "candidates"`` (longer T)       -> ``("extended", module)``
    * otherwise (production objectives)         -> ``("registry", REGISTRY)``

    The handle is returned lazily so callers that only need the family
    label don't import the heavy metric modules.
    """
    if timescale <= SHORT_BLOCK_MAX_YEARS:
        from src.metrics import short_block
        return "short_block", short_block
    if mode == "candidates":
        from src.metrics import extended
        return "extended", extended
    from src.metrics.drought_metrics import REGISTRY
    return "registry", REGISTRY


def chars_fn_for(timescale: int, monthly_1d: np.ndarray):
    """Return an ``evaluate``-compatible ``chars_fn(syn_1d, syn_2d) -> dict``.

    For ``T <= SHORT_BLOCK_MAX_YEARS`` this is the short-block factory
    (:func:`src.experiment.data.make_short_block_chars_fn`). For longer
    traces it returns ``None``, signalling callers to use the default
    pre-fitted-SSI path inside
    :func:`src.experiment.runner.run_experiment` (kept inline there so
    the SSI calculator is fitted once per run, not per evaluation).
    """
    if timescale <= SHORT_BLOCK_MAX_YEARS:
        from src.experiment.data import make_short_block_chars_fn
        return make_short_block_chars_fn(monthly_1d)
    return None
