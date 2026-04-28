"""MOEA-FIND plotting package.

Modular, manuscript-specific figure building blocks. Each submodule carries
a header comment mapping its functions to the manuscript section / figure
they produce. Scripts under ``workflows/0N_<stage>/`` import from here rather than
defining their own plotting code.

Submodules
----------
style              shared rcParams and color palette
analytic           Figs 1-3 (Manhattan concept, 2D/3D PoC, epsilon-NFE sweep)
drought_space      reusable drought-space scatter/marginals (Figs 5, 8)
coverage           Fig 7 library-vs-MOEA-FIND coverage comparison
trace_diagnostics  Fig 6 plausibility panels (acf, FDC, seasonal cycle)
convergence        SI Borg convergence diagnostics
"""

from src.plotting import (  # noqa: F401
    analytic,
    coverage,
    convergence,
    drought_space,
    style,
    trace_diagnostics,
)
