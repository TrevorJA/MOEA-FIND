"""Script 07 — Event-level Kirsch MOEA-FIND (manuscript §6.4).

Short-trace (5-10 year) Kirsch formulation with event-level objectives
(event duration, peak intensity, cumulative severity, onset month) rather
than the trace-level aggregates used in script 04. Addresses DD-01 Option B.

STATUS (2026-04-13): SCAFFOLDED. The short-trace objective helpers in
src/objectives.py are not yet implemented; this driver will import them
once available. Runs end-to-end under --dry-run for pipeline validation.

Run:
    python scripts/07_event_level_kirsch.py --nfe 20000 --n-years 8 --dry-run
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

OUTPUT_SLUG = "exp07_event_level_kirsch"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--nfe", type=int, default=20_000)
    p.add_argument("--n-years", type=int, default=8,
                   help="Short-trace length (event-focused).")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--ssi", type=int, default=3, choices=[1, 3, 6, 12])
    p.add_argument("--dry-run", action="store_true",
                   help="Validate pipeline without calling unimplemented helpers.")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "config.json").write_text(json.dumps({
        "script": "07_event_level_kirsch.py",
        "manuscript_section": "§6.4 Event-level Formulation",
        "status": "SCAFFOLDED",
        "nfe": args.nfe, "n_years": args.n_years, "seed": args.seed, "ssi": args.ssi,
    }, indent=2))

    if args.dry_run:
        print("[07] dry-run OK. Waiting on src.objectives event-level helpers "
              "(duration, peak_intensity, cumulative_severity, onset_month).")
        return

    raise NotImplementedError(
        "Event-level objectives not yet implemented. "
        "See notes/design_decisions.md DD-01 and notes/publication_plan.md Phase B5."
    )


if __name__ == "__main__":
    main()
