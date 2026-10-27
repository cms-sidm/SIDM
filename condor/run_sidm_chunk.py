#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path

from coffea import processor
import coffea.util
import yaml

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


NTUPLE_DIR = Path("sidm/configs/ntuples")


def _store_path(path):
    """The /store/... part of a file path or URL, so that one file reached through different
    xrootd servers, or through the xcache placeholder of the ntuple configs, compares equal.
    A path without /store/ is returned unchanged."""
    text = str(path)
    index = text.find("/store/")
    return text[index:] if index >= 0 else text


def find_sample_cfg(sample, files, primary_yaml, primary_section, ntuple_dir=NTUPLE_DIR):
    """The metadata entry (is_data, skim_factor, year) of the ntuple section the chunk's files
    come from, with the config file and section it was read from.

    One sample name can be listed in several sections: a background in its skimmed and its
    unskimmed section, whose skim factors differ (0.0067 against 1.0 for DYJetsToLL_M10to50),
    a signal in several signal configs and in the ffNtuple config.  The name alone therefore
    does not identify the metadata.  Every section of `primary_yaml` and of the configs in
    `ntuple_dir` that lists the sample is a candidate, and the one whose directory -- the
    section path followed by the sample path -- holds every file of the chunk is used.  No
    such section, or more than one, raises with every candidate named.  `primary_section`
    only orders the candidates in that message.
    """
    primary_yaml = Path(primary_yaml)
    yamls = [primary_yaml] + sorted(path for path in Path(ntuple_dir).glob("*.yaml")
                                    if path.resolve() != primary_yaml.resolve())
    candidates = []
    for yaml_path in yamls:
        if not yaml_path.exists():
            continue
        with open(yaml_path) as handle:
            cfg = yaml.safe_load(handle) or {}
        sections = sorted(cfg, key=lambda name: name != primary_section)
        for section in sections:
            block = cfg[section]
            if not isinstance(block, dict):
                continue
            entry = (block.get("samples") or {}).get(sample)
            if entry is None:
                continue
            entry = entry if isinstance(entry, dict) else {}
            directory = _store_path(f"{block.get('path', '')}{entry.get('path', '')}")
            candidates.append((entry, str(yaml_path), section, directory.rstrip("/") + "/"))

    matches = [candidate for candidate in candidates
               if all(_store_path(name).startswith(candidate[3]) for name in files)]
    if len(matches) == 1:
        entry, yaml_path, section, _ = matches[0]
        return entry, yaml_path, section
    listed = "; ".join(f"{yaml_path}:{section} -> {directory}"
                       for _, yaml_path, section, directory in candidates) or "none"
    problem = "no section" if not matches else f"{len(matches)} sections"
    raise KeyError(
        f"{problem} of the ntuple configs listing sample={sample!r} holds every file of this "
        f"chunk (first file: {files[0] if files else None!r}). Candidates: {listed}"
    )


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
    parser.add_argument("--backgrounds-yaml", default="sidm/configs/ntuples/backgrounds.yaml")
    parser.add_argument("--backgrounds-section", default="skimmed_llpNanoAOD_v2")
    # channel definitions, relative to sidm/ as SidmProcessor expects
    parser.add_argument("--selections-cfg", default="configs/selections.yaml")

    args = parser.parse_args()

    files = read_filelist(args.filelist)

    sample_cfg, sample_cfg_path, sample_cfg_section = find_sample_cfg(
        args.sample,
        files,
        args.backgrounds_yaml,
        args.backgrounds_section,
    )

    metadata = {
        "is_data": bool(sample_cfg.get("is_data", False)),
        "skim_factor": float(sample_cfg.get("skim_factor", 1.0)),
        "year": str(sample_cfg.get("year", "2018")),
    }

    fileset = {
        args.sample: {
            "files": files,
            "metadata": metadata,
        }
    }

    print("Repo parent:", repo_parent)
    print("Sample:", args.sample)
    print("Number of files:", len(files))
    print("Metadata source:", f"{sample_cfg_path}:{sample_cfg_section}")
    print("Metadata:", metadata)
    print("Selections:", args.selections_cfg)
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
        selections_cfg=args.selections_cfg,
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
