"""Inline smoke test for workflows/08_nyc_sensitivity/run_sa.py.

Runs the driver against synthetic chars + bank fixtures with tiny
bootstrap settings to confirm:
- the (X, Y) wiring works,
- every output file lands at the expected path,
- the long-form indices DataFrame schema is correct,
- selection_log.json renders.

Not part of the standard pytest suite — invoke explicitly:
    python tests/smoke_run_sa.py
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    rng = np.random.default_rng(0)
    N = 200
    factor_names = ["mean_severity", "mean_magnitude", "time_in_drought_fraction"]

    chars_rows = [
        {
            "mean_severity": float(rng.uniform(0.5, 3.0)),
            "mean_magnitude": float(rng.uniform(1.0, 30.0)),
            "time_in_drought_fraction": float(rng.uniform(0.0, 0.9)),
        }
        for _ in range(N)
    ]
    xs = np.array([
        [r["mean_severity"], r["mean_magnitude"], r["time_in_drought_fraction"]]
        for r in chars_rows
    ])

    bank = pd.DataFrame({
        "nyc_min_storage_frac": (
            0.7 - 0.1 * xs[:, 0] - 0.005 * xs[:, 1] - 0.3 * xs[:, 2]
            + 0.05 * rng.standard_normal(N)
        ),
        "nyc_drawdown_days_below_0.25": np.clip(
            50.0 + 30.0 * xs[:, 2] + 5.0 * xs[:, 0] + 10.0 * rng.standard_normal(N),
            0, None,
        ).astype(int),
        "montague_flow_reliability": np.clip(
            1.0 - 0.05 * xs[:, 2] - 0.01 * xs[:, 0] + 0.02 * rng.standard_normal(N),
            0, 1,
        ),
        "montague_flow_vulnerability": np.clip(
            0.05 + 0.02 * xs[:, 0] + 0.01 * xs[:, 2] + 0.005 * rng.standard_normal(N),
            0, None,
        ),
    }, index=pd.Index([str(i) for i in range(N)], name="realization_id"))

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        chars_json = td / "chars.json"
        chars_json.write_text(json.dumps({
            "objective_keys": factor_names,
            "pareto_chars": chars_rows,
            "metric_set": "primary",
            "variant": "smoke_test",
        }))
        bank_path = td / "metric_bank.parquet"
        bank.to_parquet(bank_path)
        out_dir = td / "sa_out"

        cmd = [
            sys.executable,
            str(Path(__file__).resolve().parents[1]
                / "workflows" / "08_nyc_sensitivity" / "run_sa.py"),
            "--bank", str(bank_path),
            "--chars", str(chars_json),
            "--methods", "delta",
            "--n-bootstrap", "10",
            "--n-convergence-replicates", "3",
            "--convergence-sizes", "100",
            "--seed", "0",
            "--output-dir", str(out_dir),
        ]
        print("running:", " ".join(cmd))
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if res.returncode != 0:
            print("=== STDOUT ===\n" + res.stdout)
            print("=== STDERR ===\n" + res.stderr)
            sys.exit(res.returncode)

        runs = list(out_dir.glob("sa__*"))
        assert runs, f"no run dirs created under {out_dir}"
        slug_dir = runs[0]
        print(f"slug: {slug_dir.name}")

        results_dir = slug_dir / "results"
        for f in (
            "indices_delta.parquet",
            "bootstrap_delta.parquet",
            "rank_stability_delta.parquet",
            "rank_stability_summary_delta.parquet",
            "convergence_delta.parquet",
            "cross_method_rank_corr.parquet",
            "cross_outcome_rank_corr.parquet",
            "selection_log.json",
        ):
            assert (results_dir / f).exists(), f"missing {f}"
            size = (results_dir / f).stat().st_size
            print(f"  {f}: {size} bytes")

        idx = pd.read_parquet(results_dir / "indices_delta.parquet")
        print(f"indices_delta schema: {list(idx.columns)}")
        assert {"outcome", "factor", "headline_index", "delta", "ci_lo", "ci_hi"}.issubset(
            set(idx.columns)
        )
        assert set(idx["outcome"].unique()) >= {
            "nyc_min_storage_frac",
            "nyc_drawdown_days_below_0.25",
            "nyc_drawdown_days_below_0.25_log1p",
            "montague_flow_reliability",
            "montague_flow_vulnerability",
            "montague_flow_vulnerability_log1p",
        }

        sel = json.loads((results_dir / "selection_log.json").read_text())
        for oc, decisions in sel.items():
            if oc == "_meta":
                continue
            print(f"  {oc}: anchor={decisions['_anchor']}")

    print("smoke test PASSED")


if __name__ == "__main__":
    main()
