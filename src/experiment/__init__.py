"""src.experiment — MOEA-FIND experiment harness.

Successor to the former flat ``src.experiment_utils`` module, now split
by concern:

* :mod:`src.experiment.data`      — USGS prep & historical characterization
* :mod:`src.experiment.anti_ideal`— anti-ideal placement & epsilons
* :mod:`src.experiment.runner`    — the MOEA orchestrator
* :mod:`src.experiment.plots`     — cross-generator comparison figure
* :mod:`src.experiment.config`    — experiment configuration dataclass

The public callables are re-exported here so callers can simply
``from src.experiment import run_experiment, prepare_data, ...``.
"""

from src.experiment.data import (
    prepare_data,
    compute_historical_ssi_chars,
    extract_pareto_maxes,
    compute_historical_short_block_chars,
    make_short_block_chars_fn,
)
from src.experiment.anti_ideal import compute_ssi_anti_ideal, build_epsilons
from src.experiment.runner import run_experiment
from src.experiment.plots import plot_comparison

__all__ = [
    "prepare_data",
    "compute_historical_ssi_chars",
    "extract_pareto_maxes",
    "compute_historical_short_block_chars",
    "make_short_block_chars_fn",
    "compute_ssi_anti_ideal",
    "build_epsilons",
    "run_experiment",
    "plot_comparison",
]
