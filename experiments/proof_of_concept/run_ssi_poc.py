"""DEPRECATED: SSI POC with BlockBootstrap and ParametricGenerator has been removed.

This script previously compared BlockBootstrap (empirical, KDE) and ParametricGenerator
with SSI-based objectives. BlockBootstrap has been removed in favor of the Kirsch
nonparametric method.

Use run_kirsch_poc.py instead for SSI-based experiments with Kirsch:
    python experiments/proof_of_concept/run_kirsch_poc.py --ssi 3 --nfe 5000
    python experiments/proof_of_concept/run_kirsch_poc.py --compare --ssi 3

For parametric generator experiments, see run_parametric_poc.py.
"""

raise RuntimeError(
    "run_ssi_poc.py has been removed. "
    "Use run_kirsch_poc.py for Kirsch-based SSI experiments or "
    "run_parametric_poc.py for parametric generator experiments."
)
