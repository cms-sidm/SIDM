#!/usr/bin/env python3
"""Check the unweighted_hist field of .meta.yaml sidecars on EOS against their outputs.

A merged output records in its own metadata whether its hists are unweighted. The sidecar
beside it can say otherwise when the merge that wrote it was run with an --unweighted-hist
flag that did not match the chunks. This script compares the two for every sidecar under
an EOS directory and reports the ones that disagree. With --apply it rewrites that one
field in those sidecars; without it nothing is written. Sidecars that state no mode
(a null field), and outputs that record none, are listed and left alone.

    python sidm/scripts/fix_sidecar_hist_mode.py \\
      --eos-dir /store/group/lpcmetx/SIDM/coffea_outputs/$USER/my_study
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile

import yaml
from coffea.util import load

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from sidm.scripts.merge_coffea_chunks_eos import (  # pylint: disable=wrong-import-position
    hists_unweighted, sample_outputs, xrdcp_from_eos, xrdcp_to_eos)

SIDECAR_SUFFIX = ".meta.yaml"
FIELD = re.compile(r"^unweighted_hist:.*$", flags=re.MULTILINE)


def eos_find(eos_dir):
    """All file paths under an EOS directory, e.g. /store/user/<user>/sidm_condor/Merged_v1."""
    result = subprocess.run(["xrdfs", "root://cmseos.fnal.gov", "ls", "-R", eos_dir],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"cannot list {eos_dir}: {result.stderr.strip()}")
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def recorded_hist_mode(output):
    """The unweighted_hist setting an output records: one bool, or None when no sample
    records it or the samples disagree."""
    values = set()
    for sout in sample_outputs(output).values():
        if isinstance(sout, dict) and "metadata" in sout:
            values.add(hists_unweighted(sout["metadata"]))
    values.discard(None)
    return values.pop() if len(values) == 1 else None


def hist_mode_status(text, output):
    """Compare a sidecar's text with its output: the status and the mode the output records"""
    recorded = recorded_hist_mode(output)
    stated = yaml.safe_load(text).get("unweighted_hist")
    if recorded is None:
        return "output records no mode", None
    if stated is None:
        return "sidecar states no mode", recorded
    return ("consistent" if stated is recorded else "wrong"), recorded


def with_hist_mode(text, unweighted):
    """The sidecar text with its top-level unweighted_hist field set, the rest untouched."""
    new, n = FIELD.subn(f"unweighted_hist: {'true' if unweighted else 'false'}", text, count=1)
    if n != 1:
        raise ValueError("sidecar has no top-level unweighted_hist field")
    return new


def main():
    """Command-line entry point"""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", maxsplit=1)[0])
    parser.add_argument("--eos-dir", required=True,
                        help="EOS directory to search, subdirectories included, "
                             "e.g. /store/group/lpcmetx/SIDM/coffea_outputs/$USER/my_study")
    parser.add_argument("--apply", action="store_true",
                        help="Rewrite the field in the sidecars that disagree with their "
                             "output. Default: only report them.")
    args = parser.parse_args()

    paths = set(eos_find(args.eos_dir.rstrip("/")))
    sidecars = sorted(path for path in paths if path.endswith(SIDECAR_SUFFIX))
    print(f"Found {len(sidecars)} sidecars under {args.eos_dir}")

    counts = {"consistent": 0, "wrong": 0, "sidecar states no mode": 0,
              "output records no mode": 0, "no output beside it": 0}
    with tempfile.TemporaryDirectory() as tmp:
        for sidecar in sidecars:
            output_path = sidecar[:-len(SIDECAR_SUFFIX)] + ".coffea"
            if output_path not in paths:
                counts["no output beside it"] += 1
                print(f"NO OUTPUT BESIDE IT: {sidecar}")
                continue
            local_sidecar = os.path.join(tmp, os.path.basename(sidecar))
            local_output = os.path.join(tmp, os.path.basename(output_path))
            xrdcp_from_eos(sidecar, local_sidecar)
            xrdcp_from_eos(output_path, local_output)
            with open(local_sidecar, encoding="utf8") as f:
                text = f.read()
            status, recorded = hist_mode_status(text, load(local_output))
            os.remove(local_output)
            counts[status] += 1
            if status == "wrong":
                print(f"WRONG: {sidecar} (the output records unweighted_hist = {recorded})")
                if args.apply:
                    with open(local_sidecar, "w", encoding="utf8") as f:
                        f.write(with_hist_mode(text, recorded))
                    xrdcp_to_eos(local_sidecar, os.path.dirname(sidecar))
            elif status != "consistent":
                print(f"{status.upper()}: {sidecar}")
            os.remove(local_sidecar)

    print("=" * 80)
    print("\n".join(f"{label:>24}: {n}" for label, n in counts.items()))
    if counts["wrong"]:
        print("The wrong sidecars were rewritten." if args.apply
              else "Nothing was written. Run again with --apply to rewrite the wrong sidecars.")


if __name__ == "__main__":
    main()
