#!/usr/bin/env python3

import argparse
import os
import glob
import shutil
import subprocess
import sys
from collections import defaultdict

from coffea.util import load, save
from coffea import processor

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from sidm.tools.metadata import write_run_metadata


def run_cmd(cmd, check=True):
    print("+", " ".join(cmd))
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.stdout.strip():
        print(result.stdout.strip())

    if result.stderr.strip():
        print(result.stderr.strip())

    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join(cmd)}")

    return result


def eos_ls(eos_dir):
    """
    List files in an EOS directory using xrdfs.
    eos_dir should look like:
      /store/user/scampbel/sidm_condor/BackgroundChunks_v1
    """
    result = run_cmd(
        ["xrdfs", "root://cmseos.fnal.gov", "ls", eos_dir],
        check=True,
    )

    files = [
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip().endswith(".coffea")
    ]

    return files


def eos_mkdir(eos_dir):
    run_cmd(
        ["xrdfs", "root://cmseos.fnal.gov", "mkdir", "-p", eos_dir],
        check=True,
    )


def xrdcp_from_eos(eos_path, local_path):
    """
    eos_path should look like:
      /store/user/scampbel/...
    """
    url = f"root://cmseos.fnal.gov/{eos_path}"
    run_cmd(["xrdcp", "-f", url, local_path], check=True)


def xrdcp_to_eos(local_path, eos_dir):
    """
    Copy local merged file to EOS output directory.
    """
    basename = os.path.basename(local_path)
    url = f"root://cmseos.fnal.gov/{eos_dir}/{basename}"
    run_cmd(["xrdcp", "-f", local_path, url], check=True)


def parse_sample_chunk(path):
    """
    Expected chunk filename:
      SAMPLE_CHUNK.coffea

    Examples:
      QCD_Pt800To1000_407.coffea
      2Mu2E_200GeV_0p25GeV_0p1mm_3.coffea
    """
    fname = os.path.basename(path)
    stem = fname.replace(".coffea", "")

    try:
        sample, chunk = stem.rsplit("_", 1)
        int(chunk)
    except Exception:
        raise RuntimeError(f"Could not parse sample/chunk from filename: {fname}")

    return sample, chunk


def main():
    parser = argparse.ArgumentParser(
        description="Merge SIDM Coffea chunk files from EOS and copy merged outputs back to EOS."
    )

    parser.add_argument(
        "--input-eos-dir",
        required=True,
        help="EOS input directory, e.g. /store/user/scampbel/sidm_condor/BackgroundChunks_v1",
    )

    parser.add_argument(
        "--output-eos-dir",
        required=True,
        help="EOS output directory, e.g. /store/user/scampbel/sidm_condor/BackgroundMerged_v1",
    )

    parser.add_argument(
        "--workdir",
        default="merge_tmp",
        help="Local temporary working directory",
    )

    parser.add_argument(
        "--keep-local",
        action="store_true",
        help="Keep local temporary files after merging",
    )

    # Metadata sidecar options. If --filelists-dir is supplied a .meta.yaml is
    # written next to each merged .coffea (and xrdcp'd to EOS) recording the
    # full list of ROOT files, selections, hist collections, schema, etc.
    parser.add_argument(
        "--filelists-dir",
        default=None,
        help=(
            "Directory containing per-chunk filelists (e.g. condor/filelists)."
            " Required to write the .meta.yaml sidecar."
        ),
    )
    parser.add_argument(
        "--selections",
        default=None,
        help="Comma-separated selection names (e.g. 'base'). Saved into sidecar.",
    )
    parser.add_argument(
        "--hist-collections",
        default=None,
        help="Comma-separated hist-collection names (e.g. 'muon_base'). Saved into sidecar.",
    )
    parser.add_argument("--schema", default="LLPNanoAODSchema", help="Schema name saved into sidecar.")
    parser.add_argument("--chunksize", type=int, default=50_000, help="Chunksize saved into sidecar.")
    parser.add_argument(
        "--unweighted-hist",
        action="store_true",
        help="Records `unweighted_hist: true` in the sidecar (must match how chunks were run).",
    )

    parser.add_argument(
        "--delete-chunks-after-merge",
        action="store_true",
        help=(
            "After a successful merge, delete the per-chunk .coffea files from the input EOS"
            " directory. Off by default since deletion is irreversible; opt in once a merge"
            " has been verified."
        ),
    )

    args = parser.parse_args()

    input_eos_dir = args.input_eos_dir.rstrip("/")
    output_eos_dir = args.output_eos_dir.rstrip("/")

    local_chunk_dir = os.path.join(args.workdir, "chunks")
    local_merged_dir = os.path.join(args.workdir, "merged")

    os.makedirs(local_chunk_dir, exist_ok=True)
    os.makedirs(local_merged_dir, exist_ok=True)

    print("=" * 80)
    print("Input EOS directory:")
    print(input_eos_dir)
    print("Output EOS directory:")
    print(output_eos_dir)
    print("Local workdir:")
    print(args.workdir)
    print("=" * 80)

    eos_mkdir(output_eos_dir)

    eos_files = eos_ls(input_eos_dir)
    print(f"Found {len(eos_files)} .coffea chunk files on EOS")

    if len(eos_files) == 0:
        raise RuntimeError(f"No .coffea files found in {input_eos_dir}")

    groups = defaultdict(list)

    for eos_path in eos_files:
        sample, chunk = parse_sample_chunk(eos_path)
        groups[sample].append(eos_path)

    print(f"Found {len(groups)} samples")

    for sample, eos_paths in sorted(groups.items()):
        print("=" * 80)
        print(f"Merging sample: {sample}")
        print(f"Chunks: {len(eos_paths)}")

        local_paths = []

        for i, eos_path in enumerate(sorted(eos_paths)):
            basename = os.path.basename(eos_path)
            local_path = os.path.join(local_chunk_dir, basename)

            if not os.path.exists(local_path):
                print(f"[{i+1}/{len(eos_paths)}] Copying {basename}")
                xrdcp_from_eos(eos_path, local_path)
            else:
                print(f"[{i+1}/{len(eos_paths)}] Already local: {basename}")

            local_paths.append(local_path)

        print(f"Loading and accumulating {len(local_paths)} chunks for {sample}")

        outputs = []

        for i, path in enumerate(sorted(local_paths)):
            if i % 50 == 0:
                print(f"  loading chunk {i+1}/{len(local_paths)}")
            outputs.append(load(path))

        merged = processor.accumulate(outputs)

        local_merged_path = os.path.join(local_merged_dir, f"{sample}.coffea")
        save(merged, local_merged_path)

        size_mb = os.path.getsize(local_merged_path) / 1024 / 1024
        print(f"Wrote local merged file: {local_merged_path} ({size_mb:.2f} MB)")

        print("Copying merged file to EOS")
        xrdcp_to_eos(local_merged_path, output_eos_dir)

        if args.filelists_dir:
            sample_files = []
            for fl in sorted(glob.glob(os.path.join(args.filelists_dir, f"{sample}_*.txt"))):
                with open(fl) as f:
                    sample_files.extend(line.strip() for line in f if line.strip())

            selections = [s.strip() for s in (args.selections or "").split(",") if s.strip()]
            hist_collections = [s.strip() for s in (args.hist_collections or "").split(",") if s.strip()]

            sidecar = write_run_metadata(
                local_merged_path,
                fileset={sample: {"files": sample_files, "metadata": {}}},
                selections=selections,
                hist_collections=hist_collections,
                schema=args.schema,
                chunksize=args.chunksize,
                unweighted_hist=args.unweighted_hist,
                extra={"n_chunks": len(local_paths)},
                sidm_root=_REPO_ROOT,
            )
            print(f"Wrote sidecar: {sidecar}")
            xrdcp_to_eos(sidecar, output_eos_dir)
        else:
            print("--filelists-dir not provided; skipping .meta.yaml sidecar.")

        if args.delete_chunks_after_merge:
            for eos_path in eos_paths:
                run_cmd(["xrdfs", "root://cmseos.fnal.gov", "rm", eos_path], check=True)
            print(f"Deleted {len(eos_paths)} chunks for {sample} from {input_eos_dir}")

    print("=" * 80)
    print("All samples merged and copied to EOS")

    if not args.keep_local:
        print(f"Removing local workdir: {args.workdir}")
        shutil.rmtree(args.workdir, ignore_errors=True)

    print("Done")


if __name__ == "__main__":
    main()