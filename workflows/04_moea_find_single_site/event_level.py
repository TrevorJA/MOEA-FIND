"""Event-level Kirsch MOEA-FIND scaffold (SI-H).

Short-trace (5-10 year) Kirsch formulation with event-level objectives
(event duration, peak intensity, cumulative severity, onset month) rather
than the trace-level aggregates used by run_moea_find.py.

Status: SCAFFOLDED. The event-level objective helpers in
``src/objectives.py`` are not yet implemented; this driver runs end to
end under ``--dry-run`` for pipeline validation.

Outputs to ``outputs/04_moea_find_single_site/event_level/<slug>/``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.io_paths.slugs import moea_slug  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "event_level"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--nfe", type=int, required=True)
    p.add_argument("--n-years", type=int, required=True,
                   help="Short-trace length (event-focused).")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--ssi", type=int, required=True, choices=[1, 3, 6, 12])
    p.add_argument("--dry-run", action="store_true",
                   help="Validate pipeline without calling unimplemented helpers.")
    args = p.parse_args()

    out = stage_output_dir(STAGE, DRIVER, moea_slug(
        mode="event", n_years=args.n_years, nfe=args.nfe,
        seed=args.seed, ssi=args.ssi,
    ))
    (out / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER, "status": "SCAFFOLDED",
        "nfe": args.nfe, "n_years": args.n_years, "seed": args.seed, "ssi": args.ssi,
    }, indent=2))

    if args.dry_run:
        print(f"[04/event_level] dry-run OK; output_dir={out}")
        print("  waiting on src.metrics.objectives event-level helpers "
              "(duration, peak_intensity, cumulative_severity, onset_month).")
        return

    raise NotImplementedError(
        "Event-level objectives not yet implemented; see DD-01."
    )


if __name__ == "__main__":
    main()
