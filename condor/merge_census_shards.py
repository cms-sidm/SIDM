#!/usr/bin/env python3
"""Fold the per-job census JSON shards (on EOS) into one manifest.

xrdfs-ls the shard directory, xrdcp each shard down, concatenate their `files` rows,
rebuild the per-sample rollup with sidm.tools.file_census, attach a provenance meta
block, and write the merged manifest. Mirrors merge_coffea_chunks_eos.py's EOS I/O.

    python condor/merge_census_shards.py \\
        --shard-eos-dir /store/user/$USER/sidm_census/backgrounds_v2 \\
        --out census_backgrounds_v2.json --run-id backgrounds_v2 --mode shallow
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

repo_parent = Path.cwd()
if str(repo_parent) not in sys.path:
    sys.path.insert(0, str(repo_parent))

from sidm.tools import file_census as fc


def xrdfs_ls(redirector, eos_dir):
    out = subprocess.run(["xrdfs", redirector, "ls", eos_dir],
                         capture_output=True, text=True)
    return [l for l in out.stdout.split() if l.endswith(".json")]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shard-eos-dir", required=True)
    p.add_argument("--redirector", default="root://cmseos.fnal.gov")
    p.add_argument("--out", required=True)
    p.add_argument("--run-id", default="manual")
    p.add_argument("--mode", default="shallow")
    p.add_argument("--source-yaml", default="", help="recorded in the manifest meta")
    args = p.parse_args()

    shards = xrdfs_ls(args.redirector, args.shard_eos_dir)
    print(f"Found {len(shards)} shards under {args.shard_eos_dir}")
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for i, eos_path in enumerate(sorted(shards)):
            local = os.path.join(tmp, os.path.basename(eos_path))
            subprocess.run(["xrdcp", "-f", "-s", f"{args.redirector}/{eos_path}", local],
                           check=True)
            with open(local) as f:
                rows.extend(json.load(f).get("files", []))
            if (i + 1) % 50 == 0:
                print(f"  merged {i+1}/{len(shards)} shards")

    started = rows[0].get("probed_utc", fc._utc()) if rows else fc._utc()
    ended = fc._utc()
    manifest = {
        "schema_version": 2,
        "run_id": args.run_id,
        "meta": fc._provenance(args.source_yaml or args.out, None, "ALL", None,
                               args.redirector, args.mode, started, ended, None),
        "files": rows,
        "samples": fc._rollup(rows),
    }
    with open(args.out, "w") as f:
        json.dump(manifest, f, indent=2)
    print(fc.render_census(manifest))
    print(f"\nWrote {args.out} ({len(rows)} files)")


if __name__ == "__main__":
    main()
