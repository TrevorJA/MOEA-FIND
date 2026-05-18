"""Behavior gate for the Phase-6 Tier-A dedup (src.metrics.ssi_common).

`conftest.py` mocks synhydro during pytest, so the golden equality check
cannot run as a normal test. This standalone script uses the REAL
synhydro install:

    # before the refactor
    venv/bin/python tests/golden/ssi_common_gate.py capture
    # after the refactor
    venv/bin/python tests/golden/ssi_common_gate.py verify

`verify` exits non-zero if any value of
``compute_ssi_drought_characteristics`` or ``compute_ssi_event_metrics``
changed vs the captured golden, proving the dedup is behavior-preserving.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

GOLDEN = Path(__file__).with_name("ssi_common_golden.json")


def _fixture_ssi_series():
    """Deterministic monthly-flow series with embedded multi-year droughts."""
    rng = np.random.default_rng(20260518)
    n_years = 40
    base = rng.lognormal(mean=6.0, sigma=0.4, size=n_years * 12).astype(float)
    # Carve two distinct drought spells (scale flows down hard).
    base[60:78] *= 0.18    # ~1.5 yr severe
    base[300:312] *= 0.30  # ~1 yr moderate, later in the record
    from src.metrics.objectives import compute_ssi
    ssi, _ = compute_ssi(base, timescale=3)
    return ssi, base


def _snapshot():
    from src.metrics.objectives import compute_ssi_drought_characteristics
    from src.metrics.extended import compute_ssi_event_metrics

    ssi, flows = _fixture_ssi_series()
    chars = compute_ssi_drought_characteristics(ssi, monthly_flows=flows)
    chars_no_flow = compute_ssi_drought_characteristics(ssi)
    ev = compute_ssi_event_metrics(ssi, suffix="")
    ev_sfx = compute_ssi_event_metrics(ssi, suffix="_ssi12")

    def norm(d):
        # JSON round-trip normalises numpy scalars / float repr.
        return json.loads(json.dumps({k: (None if v is None else float(v))
                                      for k, v in d.items()}, sort_keys=True))

    return {
        "chars_with_flow": norm(chars),
        "chars_no_flow": norm(chars_no_flow),
        "event_metrics": norm(ev),
        "event_metrics_suffixed": norm(ev_sfx),
    }


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "verify"
    snap = _snapshot()
    if mode == "capture":
        GOLDEN.write_text(json.dumps(snap, indent=2, sort_keys=True))
        print(f"captured golden -> {GOLDEN}")
        for k, v in snap.items():
            print(f"  {k}: {len(v)} keys")
        return 0
    if mode == "verify":
        if not GOLDEN.exists():
            print("ERROR: no golden file; run `capture` first", file=sys.stderr)
            return 2
        want = json.loads(GOLDEN.read_text())
        if snap == want:
            print("GOLDEN GATE PASSED — Tier-A dedup is behavior-preserving")
            return 0
        print("GOLDEN GATE FAILED — outputs changed:", file=sys.stderr)
        for sect in want:
            for key in sorted(set(want[sect]) | set(snap.get(sect, {}))):
                a, b = want[sect].get(key), snap.get(sect, {}).get(key)
                if a != b:
                    print(f"  [{sect}] {key}: golden={a!r} now={b!r}",
                          file=sys.stderr)
        return 1
    print(f"unknown mode {mode!r} (use capture|verify)", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
