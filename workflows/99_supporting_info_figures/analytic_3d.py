"""Render the analytic-3D projections figure from cached compute outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.analytic import fig2_3d_projections  # noqa: E402

STAGE = "01_analytic_validation"
DRIVER = "analytic_3d"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--slug", required=True)
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, args.slug, create=False)
    npz = in_dir / "pareto.npz"
    if not npz.exists():
        sys.exit(f"missing {npz} — run the compute driver first")
    pareto_dvs = np.load(npz)["dvs"]

    fig_dir = stage_figure_dir(STAGE, DRIVER, args.slug)
    out_pdf = fig_dir / "analytic_3d.pdf"
    fig = fig2_3d_projections(pareto_dvs)
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)
    print(f"[plots/01/analytic_3d] {out_pdf}")


if __name__ == "__main__":
    main()
