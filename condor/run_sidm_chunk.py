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

from sidm import BASE_DIR
from sidm.tools import sidm_processor
from sidm.tools import llpnanoaodschema

NTUPLE_DIR = Path(BASE_DIR) / "configs" / "ntuples"


def read_filelist(path):
    with open(path) as f:
        return [
            line.strip()
            for line in f
            if line.strip() and not line.strip().startswith("#")
        ]


def store_path(url):
    """/store/... part of a file URL, so xcache and cmseos URLs compare equal"""
    return url[url.index("/store/"):] if "/store/" in url else url


def sample_metadata(sample, files):
    """The fileset metadata of sample, from the location YAML its files come from

    Every YAML under sidm/configs/ntuples is searched for an entry named sample
    whose directory holds the first file, the way make_job_args.py built the
    filelist. Returns (metadata, "<yaml>:<ntuple version>"), or (None, None)
    when no entry matches.
    """
    import yaml
    first = store_path(files[0])
    for cfg in sorted(NTUPLE_DIR.glob("*.yaml")):
        with open(cfg) as f:
            menu = yaml.load(f, Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader)) or {}
        for version, locations in menu.items():
            entry = (locations.get("samples") or {}).get(sample)
            if not isinstance(entry, dict):
                continue
            if first.startswith(store_path(locations["path"] + entry["path"])):
                metadata = {
                    "skim_factor": entry.get("skim_factor", 1.0),
                    "is_data": entry.get("is_data", False) or cfg.name.startswith("data"),
                    "year": entry.get("year", "2018"),
                }
                return metadata, f"{cfg.name}:{version}"
    return None, None


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
    parser.add_argument(
        "--corrections", default="",
        help="standard corrections to apply (sidm/configs/corrections.yaml): 'default' for "
             "the config's default list, or a comma-separated list of names; empty for none")
    parser.add_argument("--syst-hists", action="store_true",
                        help="add a syst axis with the weight variations to every histogram")
    parser.add_argument(
        "--monitor-corrections", default="",
        help="corrections to evaluate without applying them (comma-separated); they add "
             "with_<name> weight categories")
    parser.add_argument(
        "--allow-data", action="store_true",
        help="process collision data. Off by default because the analysis is blinded: "
             "a chunk whose files are data stops before anything is read.")

    args = parser.parse_args()

    files = read_filelist(args.filelist)
    # SidmProcessor.process reads events.metadata["is_data"] / ["skim_factor"] /
    # ["year"] without .get() defaults, so a bare {sample: files} fileset crashes
    # with KeyError. Take them from the location YAML the filelist was built
    # from: the skim factor sets the normalisation of the skimmed backgrounds.
    metadata, source = sample_metadata(args.sample, files)
    if metadata is None:
        print(f"WARNING: no location YAML lists {args.sample} with these files; "
              "using MC defaults (skim_factor 1, year 2018)")
        metadata = {"is_data": False, "skim_factor": 1.0, "year": "2018"}
    looks_like_data = metadata["is_data"] or any("/SIDM/Data/" in f for f in files)
    if looks_like_data and not args.allow_data:
        sys.exit(f"ERROR: {args.sample} is collision data ({source or 'a /SIDM/Data/ path'}). "
                 "The analysis is blinded; pass --allow-data to process it on purpose.")
    print("Sample metadata:", metadata, "from", source)
    fileset = {args.sample: {"files": files, "metadata": metadata}}

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
    if args.corrections.strip() == "default":
        corrections = True
    else:
        corrections = [x.strip() for x in args.corrections.split(",") if x.strip()]
    monitored = [x.strip() for x in args.monitor_corrections.split(",") if x.strip()]

    p = sidm_processor.SidmProcessor(
        channels,
        hist_collections,
        unweighted_hist=args.unweighted_hist,
        corrections=corrections,
        syst_hists=args.syst_hists,
        monitor_corrections=monitored or None,
    )

    try:
        output = runner.run(
            fileset,
            treename=args.treename,
            processor_instance=p,
        )
    except TypeError as e:
        # coffea's Runner unpacks the executor result without guarding the empty case:
        # when every file in this chunk has 0 events (legitimate after a skim), the
        # preprocessing produces zero work items, the executor returns None, and
        # `wrapped_out, e = executor(...)` raises "cannot unpack non-iterable NoneType
        # object". That is not a failure of this chunk -- it simply has nothing to
        # process -- so emit a valid EMPTY output (the additive identity for the merge's
        # processor.accumulate) and exit cleanly instead of crashing the whole job.
        if "unpack non-iterable NoneType" not in str(e):
            raise
        print("All files in this chunk have 0 events; saving an empty output.")
        output = {"out": {}, "processed": set(), "exception": 0}

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    coffea.util.save(output, args.output)

    print("Saved:", args.output)


if __name__ == "__main__":
    main()