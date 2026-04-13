"""Script 08 — DRB multi-site MOEA-FIND (manuscript §7.1-7.2, Fig 8).

Four DRB inflow sites (Cannonsville, Pepacton, Neversink, Delaware lateral)
with shared Kirsch bootstrap indices to preserve cross-site correlation.
Three drought objectives (frequency, mean duration, mean intensity) plus
the Manhattan-norm auxiliary objective. Designed for MM Borg under MPI.

STATUS (2026-04-13): SCAFFOLDED. Multi-site wrapper and DRB data loader
are not yet implemented. This driver documents the intended interface,
parses CLI, and runs end-to-end under --dry-run.

Run:
    sbatch scripts/08_drb_multisite_moea.slurm        # on HPC (Phase C)
    python scripts/08_drb_multisite_moea.py --dry-run # local pipeline check
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

OUTPUT_SLUG = "exp08_drb_multisite"

DEFAULT_SITES = ("cannonsville", "pepacton", "neversink", "delaware_lateral")
DEFAULT_OBJECTIVES = ("frequency", "mean_duration", "mean_avg_severity")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--nfe", type=int, default=500_000)
    p.add_argument("--n-years", type=int, default=30)
    p.add_argument("--seeds", type=int, nargs="+", default=[42, 17, 91, 5, 13, 31, 77, 109])
    p.add_argument("--sites", nargs="+", default=list(DEFAULT_SITES))
    p.add_argument("--objectives", nargs="+", default=list(DEFAULT_OBJECTIVES))
    p.add_argument("--ssi", type=int, default=3, choices=[1, 3, 6, 12])
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "config.json").write_text(json.dumps({
        "script": "08_drb_multisite_moea.py",
        "manuscript_section": "§7.1-7.2 DRB Case Study (Fig 8)",
        "status": "SCAFFOLDED",
        "nfe": args.nfe, "n_years": args.n_years,
        "seeds": args.seeds, "sites": args.sites,
        "objectives": args.objectives, "ssi": args.ssi,
    }, indent=2))

    if args.dry_run:
        print(f"[08] dry-run OK. Sites={args.sites} objectives={args.objectives} "
              f"seeds={len(args.seeds)} nfe={args.nfe}")
        return

    raise NotImplementedError(
        "Multi-site Kirsch wrapper and DRB inflow loader pending. "
        "See notes/publication_plan.md Phase C3/C4."
    )


if __name__ == "__main__":
    main()
