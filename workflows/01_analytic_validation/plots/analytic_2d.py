"""Render the analytic-2D coverage figure from cached compute outputs.

Reads ``outputs/01_analytic_validation/analytic_2d/<slug>/pareto.npz``
and writes the PDF to
``figures/01_analytic_validation/analytic_2d/<slug>/analytic_2d.pdf``.

This is a plotting-only driver — it never re-runs the MOEA. Pair with
``analytic_2d.py`` (the compute driver) which produces the input npz.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.analytic import fig2_2d_coverage_comparison  # noqa: E402

STAGE = "01_analytic_validation"
DRIVER = "analytic_2d"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True,
                   help="Variant slug under outputs/01_analytic_validation/analytic_2d/")
    p.add_argument("--seed", type=int, required=True,
                   help="Seed used by the LHS reference panel; must match the compute run")
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    npz = in_dir / "pareto.npz"
    if not npz.exists():
        sys.exit(f"missing {npz} — run the compute driver first")
    pareto_dvs = np.load(npz)["dvs"]

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    out_pdf = fig_dir / "analytic_2d.pdf"
    fig = fig2_2d_coverage_comparison(pareto_dvs, seed=args.seed)
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)
    print(f"[plots/01/analytic_2d] {out_pdf}")


if __name__ == "__main__":
    main()
