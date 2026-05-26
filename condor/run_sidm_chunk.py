#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path

from coffea import processor
import coffea.util

# Make sure parent of sidm/ is importable
repo_parent = Path.cwd()
if str(repo_parent) not in sys.path:
    sys.path.insert(0, str(repo_parent))

from sidm.tools import sidm_processor
from sidm.tools import llpnanoaodschema


def read_filelist(path):
    with open(path) as f:
        return [
            line.strip()
            for line in f
            if line.strip() and not line.strip().startswith("#")
        ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--filelist", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--treename", default="Events")
    parser.add_argument("--chunksize", type=int, default=50_000)
    parser.add_argument("--workers", type=int, default=1)

    parser.add_argument("--channels", default="base")
    parser.add_argument("--hist-collections", default="muon_base")
    parser.add_argument("--unweighted-hist", action="store_true")

    args = parser.parse_args()

    files = read_filelist(args.filelist)
    # SidmProcessor.process reads events.metadata["is_data"] / ["skim_factor"] /
    # ["year"] without .get() defaults, so a bare {sample: files} fileset crashes
    # with KeyError. Supply MC-sample defaults; data submissions should override.
    fileset = {
        args.sample: {
            "files": files,
            "metadata": {
                "is_data": False,
                "skim_factor": 1.0,
                "year": "2018",
            },
        }
    }

    print("Repo parent:", repo_parent)
    print("Sample:", args.sample)
    print("Number of files:", len(files))
    print("Files:")
    for f in files:
        print("  ", f)

    runner = processor.Runner(
        executor=processor.FuturesExecutor(workers=args.workers),
        schema=llpnanoaodschema.LLPNanoAODSchema,
        skipbadfiles=True,
        chunksize=args.chunksize,
    )

    channels = [x.strip() for x in args.channels.split(",") if x.strip()]
    hist_collections = [x.strip() for x in args.hist_collections.split(",") if x.strip()]

    p = sidm_processor.SidmProcessor(
        channels,
        hist_collections,
        unweighted_hist=args.unweighted_hist,
    )

    output = runner.run(
        fileset,
        treename=args.treename,
        processor_instance=p,
    )

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    coffea.util.save(output, args.output)

    print("Saved:", args.output)


if __name__ == "__main__":
    main()