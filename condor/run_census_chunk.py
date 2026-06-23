#!/usr/bin/env python3
"""Condor worker for the file census: probe one shard of files and write a JSON shard.

Reads a per-job filelist (one ROOT URL per line; `#` lines kept as a sentinel encode the
sample + commented-source so the merge can reconstruct the manifest rows), probes each
file with sidm.tools.file_census (shallow by default, or deep = full SidmProcessor when
--deep), and writes the result rows as a small JSON shard. run_census_job.sh xrdcp's the
shard to EOS, and merge_census_shards.py folds all shards into one manifest. Mirrors the
run_sidm_chunk.py / run_job.sh worker contract.

Filelist line format (written by make_census_args.py):
    <url>\t<sample>\t<source: live|commented>\t<skim_factor>\t<is_data>\t<year>[\t<commented_reason>]
"""

import argparse
import json
import os
import sys
from pathlib import Path

repo_parent = Path.cwd()
if str(repo_parent) not in sys.path:
    sys.path.insert(0, str(repo_parent))

from sidm.tools import file_census as fc


def read_shard(path):
    rows = []
    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("##"):
                continue
            parts = line.split("\t")
            # version is FIRST (make_census_args writes it there); _rollup requires row["version"].
            version = parts[0]
            url = parts[1] if len(parts) > 1 else ""
            sample = parts[2] if len(parts) > 2 else "UNKNOWN"
            source = parts[3] if len(parts) > 3 else "live"
            skim = float(parts[4]) if len(parts) > 4 and parts[4] else 1.0
            is_data = (parts[5].lower() == "true") if len(parts) > 5 else False
            year = parts[6] if len(parts) > 6 else "2018"
            reason = parts[7] if len(parts) > 7 else ""
            rows.append(dict(version=version, url=url, sample=sample, source=source,
                             skim_factor=skim, is_data=is_data, year=year,
                             commented_reason=reason))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--filelist", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--redirector", default="root://cmseos.fnal.gov")
    p.add_argument("--deep", action="store_true")
    p.add_argument("--channels", default="base")
    p.add_argument("--hist-collections", default="muon_base")
    p.add_argument("--chunksize", type=int, default=50000)
    args = p.parse_args()

    shard = read_shard(args.filelist)
    print(f"Probing {len(shard)} files (mode={'deep' if args.deep else 'shallow'})")
    channels = [x for x in args.channels.split(",") if x.strip()]
    hcoll = [x for x in args.hist_collections.split(",") if x.strip()]

    out_rows = []
    for i, e in enumerate(shard):
        print(f"  [{i+1}/{len(shard)}] {e['sample']}/{os.path.basename(e['url'])}")
        if args.deep:
            md = {"is_data": e["is_data"], "skim_factor": e["skim_factor"], "year": e["year"]}
            pr = fc.probe_file_deep(e["url"], e["sample"], md, channels, hcoll,
                                    redirector=args.redirector, chunksize=args.chunksize)
        else:
            pr = fc.probe_file(e["url"], redirector=args.redirector)
        row = dict(e)
        row["filename"] = os.path.basename(e["url"])
        row.update(status=pr.status, genEventSumw=pr.genEventSumw, genEventCount=pr.genEventCount,
                   events_entries=pr.events_entries, n_runs_entries=pr.n_runs_entries,
                   has_runs=pr.has_runs, n_attempts=pr.n_attempts, last_error=pr.error,
                   runs_anomaly=pr.runs_anomaly,
                   process_status=pr.process_status, process_error=pr.process_error,
                   probed_utc=fc._utc())
        out_rows.append(row)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    # n_expected lets merge_census_shards detect an under-filled shard (a job that died
    # mid-probe) distinct from a legitimately small final chunk.
    with open(args.output, "w") as f:
        json.dump({"n_expected": len(shard), "files": out_rows}, f)
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
