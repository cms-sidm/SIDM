#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path

repo_parent = Path.cwd()
if str(repo_parent) not in sys.path:
    sys.path.insert(0, str(repo_parent))

from sidm.tools import utilities


def read_samples(path):
    with open(path) as f:
        return [
            line.strip()
            for line in f
            if line.strip() and not line.strip().startswith("#")
        ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-file", required=True)
    parser.add_argument("--dataset", default="llpNanoAOD_v2")
    parser.add_argument("--location-cfg", default="signal_2mu2e_v10.yaml")
    parser.add_argument("--max-files", type=int, default=None)
    parser.add_argument("--files-per-job", type=int, default=1)
    parser.add_argument("--outdir", default="condor/filelists")
    parser.add_argument("--job-args", default="condor/job_args.txt")
    parser.add_argument("--replace-xcache", action="store_true",
                        help="Rewrite root://xcache// to root://cmseos.fnal.gov// in the "
                             "fileset before chunking (use on LPC where xcache does not resolve).")
    parser.add_argument("--census-skip", default=None,
                        help="census skip-list (path, or a name under sidm/configs/census/); its "
                             "flagged files are dropped in addition to the YAML's commented ones.")
    args = parser.parse_args()

    samples = read_samples(args.samples_file)

    fileset = utilities.make_fileset(
        samples,
        args.dataset,
        location_cfg=args.location_cfg,
        max_files=args.max_files,
        replace_xcache=args.replace_xcache,
        census_skip=args.census_skip,
    )

    os.makedirs(args.outdir, exist_ok=True)

    n_jobs = 0

    with open(args.job_args, "w") as job_args:
        for sample, sample_info in fileset.items():

            if isinstance(sample_info, dict) and "files" in sample_info:
                files = list(sample_info["files"])
            else:
                files = list(sample_info)

            for i in range(0, len(files), args.files_per_job):
                chunk = i // args.files_per_job
                chunk_files = files[i:i + args.files_per_job]

                filelist_name = f"{sample}_{chunk}.txt"
                filelist_path = Path(args.outdir) / filelist_name

                with open(filelist_path, "w") as f:
                    for root_file in chunk_files:
                        f.write(root_file + "\n")

                # job_args.txt is consumed by `condor_submit` running from inside condor/,
                # so write the path relative to that directory (e.g. filelists/foo_0.txt,
                # not condor/filelists/foo_0.txt).
                path_for_submit = os.path.relpath(filelist_path, "condor")
                job_args.write(f"{sample} {chunk} {path_for_submit}\n")
                n_jobs += 1

        print(f"Made {n_jobs} jobs")
        print(f"Wrote {args.job_args}")


if __name__ == "__main__":
    main()