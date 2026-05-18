"""Stage 06 (windowed) -- aggregate per-scenario outputs into the standard
results/ bank that stages 07/08/09 consume.

Runs ONCE after all policy_reeval_windowed array tasks finish. Scans
scenarios/<id>/ for successful runs, reuses the existing
src.pywrdrb.satisficing_metrics / src.discovery.scenario_discovery helpers per scenario,
and writes:

  results/metric_bank.parquet      (rows = successful global scenario ids)
  results/satisficing_*.{csv,json} (build_satisficing_table + save_results)
  results/window_manifest.parquet  (per-scenario window + short_spinup)
  results/aggregate_config.json    (n_total / n_ok / n_failed, window rule)

Only successful scenarios are included; pareto_chars / drought_metrics
are sliced to the same ordered subset so the satisficing table aligns.
Serial, cheap, single core (R4-safe to run on a small SLURM job).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.pywrdrb.satisficing_metrics import compute_metric_bank, write_metric_bank  # noqa: E402
from src.discovery.scenario_discovery import (  # noqa: E402
    extract_drought_levels,
    build_satisficing_table,
    save_results,
)
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "06_pywrdrb_reeval"
DRIVER = "policy_reeval"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pareto-results", type=Path, required=True)
    args = ap.parse_args()

    results = json.loads(Path(args.pareto_results).read_text())
    pareto_chars = results.get("pareto_chars", []) or []
    drought_metrics = np.asarray(results.get("drought_metrics", []), dtype=float)
    objective_keys = tuple(results.get("objective_keys", ()))
    src_slug = args.pareto_results.parent.name
    out = stage_output_dir(STAGE, DRIVER, src_slug)
    scen_root = out / "scenarios"
    results_dir = out / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    if not scen_root.exists():
        raise SystemExit(f"[06w-agg] no scenarios dir at {scen_root}")

    ok_ids, banks, drought_levels, manifest = [], [], {}, []
    n_total = sum(1 for _ in scen_root.iterdir() if _.is_dir())
    for sdir in sorted(scen_root.iterdir(), key=lambda p: int(p.name)
                       if p.name.isdigit() else 1 << 30):
        if not sdir.is_dir() or not sdir.name.isdigit():
            continue
        sid = int(sdir.name)
        status = sdir / "status.json"
        out_h5 = sdir / "simulations" / "batch_0.hdf5"
        if not status.exists() or not out_h5.exists():
            continue
        st = json.loads(status.read_text())
        if not st.get("ok"):
            continue
        try:
            b = compute_metric_bank(out_h5, ["0"])
            b = b.copy()
            # Preserve compute_metric_bank's "realization_id"-named index
            # (stages 08/09 require it); just rekey the single row to the
            # global scenario id as a string.
            b.index = pd.Index([str(sid)] * len(b), name="realization_id")
            banks.append(b)
            dl = extract_drought_levels(sdir / "simulations", ["0"])
            # remap the single realization key "0" -> str(sid)
            for _k, v in list(dl.items()):
                drought_levels[str(sid)] = v
            wj = sdir / "window.json"
            if wj.exists():
                manifest.append(json.loads(wj.read_text()))
            ok_ids.append(sid)
        except Exception as e:  # noqa: BLE001
            print(f"[06w-agg] scenario {sid} aggregation error: {e!r}",
                  flush=True)

    if not ok_ids:
        raise SystemExit("[06w-agg] no successful scenarios to aggregate")

    bank = pd.concat(banks, axis=0)
    bank_path = write_metric_bank(bank, results_dir / "metric_bank.parquet")
    # CSV companion for stage 09 (reads metric_bank.csv).
    try:
        bank.to_csv(results_dir / "metric_bank.csv", index=True)
    except Exception as e:  # noqa: BLE001
        print(f"[06w-agg] csv companion skipped: {e!r}", flush=True)

    sub_chars = [pareto_chars[i] for i in ok_ids] if pareto_chars else []
    sub_dm = drought_metrics[ok_ids] if drought_metrics.size else drought_metrics
    df = build_satisficing_table(drought_levels, sub_chars, sub_dm,
                                 objective_keys)
    save_results(df, drought_levels, results_dir)

    # Stages 08/09 treat objective_keys as literal pareto_chars columns,
    # but the coarse preset's metric NAMES carry a "_fcc" suffix while the
    # chars dict uses the raw extraction keys. Emit a SA-facing copy of
    # results.json whose objective_keys are mapped (via the registry's
    # extractor) to those raw chars columns. The _coarse 08/09 SLURMs
    # point --chars at this file.
    try:
        from src.metrics.drought_metrics import REGISTRY
        def _chars_key(name):
            m = REGISTRY.get(name)
            if m is None:
                return name
            try:
                return m.extract.__closure__[0].cell_contents
            except Exception:  # noqa: BLE001
                return name
        sa_results = dict(results)
        sa_results["objective_keys"] = [_chars_key(k) for k in objective_keys]
        sa_json = json.dumps(sa_results)
        # (a) stage-06 results/ — consumed by 08/09 via --chars.
        (results_dir / "results_sa.json").write_text(sa_json)
        # (b) stage-04 coarse dir — consumed by S07b, whose driver derives
        # src_slug from pareto_results.parent.name (must be the coarse
        # slug, not "results"). Write a sibling of the stage-04 results.json.
        (Path(args.pareto_results).parent / "results_sa.json").write_text(sa_json)
        print(f"[06w-agg] wrote results_sa.json x2 (objective_keys -> "
              f"{sa_results['objective_keys']})", flush=True)
    except Exception as e:  # noqa: BLE001
        print(f"[06w-agg] results_sa.json skipped: {e!r}", flush=True)

    if manifest:
        pd.DataFrame(manifest).to_parquet(
            results_dir / "window_manifest.parquet", index=False)
    n_short = sum(1 for m in manifest if m.get("short_spinup"))
    (results_dir / "aggregate_config.json").write_text(json.dumps({
        "src_slug": src_slug, "n_total_scenarios": n_total,
        "n_ok": len(ok_ids), "n_failed": n_total - len(ok_ids),
        "n_short_spinup": n_short,
        "window_rule": "first critical SSI-3 event: [onset-36mo, recovery+12mo] clamped",
        "metric_bank_rows": int(len(bank)),
    }, indent=2))
    print(f"[06w-agg] ok={len(ok_ids)}/{n_total} bank={bank.shape} "
          f"short_spinup={n_short} -> {bank_path}", flush=True)


if __name__ == "__main__":
    main()
