"""Helpers for fitting SynHydro's Kirsch generator on historical flows.

Centralizes the boilerplate that previously appeared inline in every
workflow script that needs a fitted ``KirschGenerator``: building a
DatetimeIndex from a flat ``(n_years, 12)`` array of monthly flows and
calling ``KirschGenerator.fit`` with the expected schema.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from synhydro.methods.generation.nonparametric.kirsch import KirschGenerator


_DEFAULT_START_DATE = "1950-10-01"


def build_kirsch_generator(
    monthly_2d: np.ndarray,
    *,
    generate_using_log_flow: bool = True,
    start_date: str = _DEFAULT_START_DATE,
) -> "KirschGenerator":
    """Fit SynHydro's KirschGenerator on a (n_years, 12) monthly array.

    Args:
        monthly_2d: Historical monthly flows shaped ``(n_years, 12)``.
        generate_using_log_flow: Forwarded to ``KirschGenerator``.
        start_date: First month of the synthetic DatetimeIndex used to fit
            the generator. The actual dates are immaterial to Kirsch fitting
            but the index must exist.

    Returns:
        A fitted ``KirschGenerator`` instance.
    """
    import pandas as pd
    from synhydro.methods.generation.nonparametric.kirsch import KirschGenerator

    gen = KirschGenerator(generate_using_log_flow=generate_using_log_flow)
    dates = pd.date_range(start=start_date, periods=monthly_2d.size, freq="MS")
    fit_df = pd.DataFrame({"flow_cfs": monthly_2d.flatten()}, index=dates)
    gen.fit(fit_df)
    return gen
