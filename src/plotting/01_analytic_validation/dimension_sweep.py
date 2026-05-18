"""Render the K=2..6 dimension-sweep figure from cached compute outputs.

Reads every ``outputs/01_analytic_validation/dimension_sweep/k{K}/``
directory and emits one combined figure under
``figures/01_analytic_validation/dimension_sweep/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.analytic import fig4_dimension_sweep  # noqa: E402

STAGE = "01_analytic_validation"
DRIVER = "dimension_sweep"


def main():
    diag_root = stage_output_dir(STAGE, DRIVER, create=False)
    if not diag_root.exists():
        sys.exit(f"missing {diag_root} — run the compute driver for K=2..6 first")

    fig_dir = stage_figure_dir(STAGE, DRIVER)
    out_pdf = fig_dir / "dimension_sweep.pdf"
    fig = fig4_dimension_sweep(diag_root)
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)
    print(f"[plots/01/dimension_sweep] {out_pdf}")


if __name__ == "__main__":
    main()
