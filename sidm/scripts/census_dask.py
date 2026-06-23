#!/usr/bin/env python3
"""Run the file census at scale on LPC Condor workers via dask (lpcjobqueue).

Enumerates every file of a location YAML (including commented ones), batches them, and
maps the shallow probe over the batches on dask workers spun up by
sidm.tools.scaleout.make_lpc_client, then folds the results into one manifest. Best for
the cheap shallow census; the DEEP (full-SidmProcessor) census is better run via the
Condor scripts (make_census_args.py / submit_census.sub), where each job is isolated.

    python sidm/scripts/census_dask.py --location-cfg backgrounds.yaml \\
        --version skimmed_llpNanoAOD_v2 --batch-size 50 --max-workers 40 \\
        --out census_backgrounds.json --run-id backgrounds_dask
"""

import argparse
import json
import os
import sys
from dataclasses import asdict
from pathlib import Path

repo = Path(__file__).resolve().parent.parent.parent
if str(repo) not in sys.path:
    sys.path.insert(0, str(repo))

from sidm.tools import file_census as fc
from sidm.tools import scaleout


def probe_batch(entries, redirector):
    """Runs on a dask worker (sidm/ shipped via UploadDirectory). entries are plain dicts."""
    from sidm.tools import file_census as _fc
    out = []
    for e in entries:
        pr = _fc.probe_file(e["url"], redirector=redirector)
        row = dict(e)
        row.update(status=pr.status, genEventSumw=pr.genEventSumw, genEventCount=pr.genEventCount,
                   events_entries=pr.events_entries, n_runs_entries=pr.n_runs_entries,
                   has_runs=pr.has_runs, n_attempts=pr.n_attempts, last_error=pr.error,
                   process_status=pr.process_status, process_error=pr.process_error,
                   probed_utc=_fc._utc())
        out.append(row)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--location-cfg", required=True)
    p.add_argument("--version", default=None)
    p.add_argument("--batch-size", type=int, default=50)
    p.add_argument("--max-workers", type=int, default=40)
    p.add_argument("--redirector", default="root://cmseos.fnal.gov")
    p.add_argument("--out", required=True)
    p.add_argument("--run-id", default="manual")
    args = p.parse_args()

    cfg = args.location_cfg
    if not os.path.isfile(cfg):
        cfg = os.path.join(repo, "sidm", "configs", "ntuples", cfg)
    started = fc._utc()
    entries = [asdict(e) for e in fc.enumerate_location(cfg, version=args.version)]
    batches = [entries[i:i + args.batch_size] for i in range(0, len(entries), args.batch_size)]
    print(f"{len(entries)} files in {len(batches)} batches; spinning up dask workers...")

    cluster, client = scaleout.make_lpc_client(min_workers=1, max_workers=args.max_workers)
    rows = []
    failed_batches = 0
    try:
        from dask.distributed import as_completed
        # retries>0 re-runs a batch whose worker was evicted/OOM'd; raise_errors=False keeps a
        # batch that still fails after retries from crashing the WHOLE run and discarding every
        # row already gathered (over thousands of preemptible LPC workers, a death is near-certain).
        futures = client.map(probe_batch, batches, redirector=args.redirector, retries=2)
        done = 0
        for fut, res in as_completed(futures, with_results=True, raise_errors=False):
            done += 1
            if isinstance(res, Exception):
                failed_batches += 1
                print(f"  batch FAILED ({type(res).__name__}: {str(res)[:120]}) -- its files "
                      f"will be missing from the manifest", file=sys.stderr)
            else:
                rows.extend(res)
            if done % 10 == 0:
                print(f"  {done}/{len(batches)} batches done "
                      f"({len(rows)} files, {failed_batches} failed)")
    finally:
        client.close()
        cluster.close()

    complete = (len(rows) == len(entries))
    meta = fc._provenance(cfg, args.version, "ALL", None, args.redirector, "shallow",
                          started, fc._utc(), None)
    meta["backend"] = f"dask (max_workers={args.max_workers})"
    meta["n_enumerated"] = len(entries)
    meta["n_probed"] = len(rows)
    meta["complete"] = complete
    manifest = {
        "schema_version": 2, "run_id": args.run_id,
        "meta": meta, "files": rows, "samples": fc._rollup(rows),
    }
    with open(args.out, "w") as f:
        json.dump(manifest, f, indent=2)
    print(fc.render_census(manifest))
    print(f"\nWrote {args.out} ({len(rows)} files, complete={complete})")
    if not complete:
        print(f"\n*** PARTIAL CENSUS: probed {len(rows)}/{len(entries)} files "
              f"({failed_batches} failed batches). Manifest is marked complete=false; "
              f"re-run before publishing. ***", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
