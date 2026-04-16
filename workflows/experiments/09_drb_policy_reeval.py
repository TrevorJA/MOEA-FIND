"""Script 09 — DRB policy re-evaluation hand-off (manuscript §7.3, Fig 9).

Consumes the MOEA-FIND ensemble produced by script 08 and exports it in the
format Pywr-DRB / NYCOptimization expect for policy re-evaluation. Runs BART
scenario discovery on the re-evaluation results.

STATUS (2026-04-13): SCAFFOLDED. Pending script 08 output and the
NYCOptimization interface specification.

Run:
    python scripts/09_drb_policy_reeval.py \
        --ensemble outputs/exp08_drb_multisite/ensemble.npz --dry-run
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

OUTPUT_SLUG = "exp09_drb_policy_reeval"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--ensemble", type=Path, required=False,
                   help="Path to ensemble produced by script 08.")
    p.add_argument("--policies", type=Path, required=False,
                   help="Path to NYCOptimization Pareto policies.")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "config.json").write_text(json.dumps({
        "script": "09_drb_policy_reeval.py",
        "manuscript_section": "§7.3 Policy Re-evaluation (Fig 9)",
        "status": "SCAFFOLDED",
        "ensemble": str(args.ensemble) if args.ensemble else None,
        "policies": str(args.policies) if args.policies else None,
    }, indent=2))

    if args.dry_run:
        print("[09] dry-run OK. Waiting on script 08 ensemble and "
              "NYCOptimization handoff spec.")
        return

    raise NotImplementedError(
        "Policy re-evaluation pipeline pending. "
        "See notes/publication_plan.md Phase C5."
    )


if __name__ == "__main__":
    main()
