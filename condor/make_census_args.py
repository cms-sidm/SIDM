#!/usr/bin/env python3
"""Enumerate every file of a location YAML (INCLUDING commented-out ones) and chunk them
into Condor census jobs. Mirrors make_job_args.py, but uses file_census.enumerate_location
(so the commented files are tested too) and writes a richer per-line shard format that
run_census_chunk.py reads and merge_census_shards.py reconstructs into a manifest.

    python condor/make_census_args.py --location-cfg backgrounds.yaml \\
        --version skimmed_llpNanoAOD_v2 --files-per-job 200
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


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--location-cfg", required=True,
                   help="path or bare name under sidm/configs/ntuples/")
    p.add_argument("--version", default=None)
    p.add_argument("--files-per-job", type=int, default=200,
                   help="census probes are cheap, so use big chunks (default 200)")
    p.add_argument("--max-files", type=int, default=None,
                   help="cap the total files enumerated (smoke tests / sampling)")
    p.add_argument("--outdir", default="condor/filelists_census")
    p.add_argument("--job-args", default="condor/census_job_args.txt")
    p.add_argument("--redirector", default="root://cmseos.fnal.gov")
    args = p.parse_args()

    cfg = args.location_cfg
    if not os.path.isfile(cfg):
        cfg = os.path.join(repo_parent, "sidm", "configs", "ntuples", cfg)
    entries = fc.enumerate_location(cfg, version=args.version)
    if args.max_files:
        entries = entries[:args.max_files]
    os.makedirs(args.outdir, exist_ok=True)

    n_jobs = 0
    with open(args.job_args, "w") as ja:
        for i in range(0, len(entries), args.files_per_job):
            chunk = entries[i:i + args.files_per_job]
            name = f"census_{i // args.files_per_job}.txt"
            fpath = os.path.join(args.outdir, name)
            with open(fpath, "w") as f:
                for e in chunk:
                    # version is FIRST: run_census_chunk reconstructs it into the row, and
                    # file_census._rollup requires r["version"]. Omitting it makes the whole
                    # Condor merge crash with KeyError('version').
                    f.write("\t".join([e.version, e.url, e.sample, e.source, str(e.skim_factor),
                                       str(e.is_data), e.year,
                                       (e.commented_reason or "").replace("\t", " ")]) + "\n")
            # path relative to condor/ (condor_submit runs from there)
            rel = os.path.relpath(fpath, "condor")
            ja.write(f"{i // args.files_per_job} {rel}\n")
            n_jobs += 1

    # Expected-count sidecar so merge_census_shards can detect a missing shard / missing job
    # and refuse to publish a PARTIAL census as complete. Named after the job-args file so two
    # campaigns (e.g. backgrounds + data_skimmed) can be chunked concurrently without clobbering.
    exp_path = os.path.splitext(args.job_args)[0] + "_expected.json"
    with open(exp_path, "w") as f:
        json.dump({"n_files": len(entries), "n_jobs": n_jobs,
                   "files_per_job": args.files_per_job,
                   "version": args.version or "ALL",
                   "source_yaml": os.path.abspath(cfg)}, f, indent=2)

    print(f"Enumerated {len(entries)} files ({sum(1 for e in entries if e.source=='commented')} "
          f"commented) -> {n_jobs} census jobs")
    print(f"Wrote {args.job_args}, {args.outdir}/*.txt, and {exp_path}")


if __name__ == "__main__":
    main()
