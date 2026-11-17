#!/usr/bin/env python3
"""Likelihood scan of one parameter of one datacard, with Combine's MultiDimFit

    python sidm/scripts/run_combine_scan.py CARD --outdir DIR --poi POI --range=LO,HI
        [--expect-signal R] [--points 40] [--freeze PARAMS] [--name TAG]

Scans -2 Delta ln L against POI, which is r or a nuisance (r then floats), on the
Asimov dataset with the signal at strength R (default 0, background only). Every
other parameter is profiled except those in --freeze: a comma-separated list, or
allConstrainedNuisances for the statistical uncertainty alone. Writes
DIR/scan_TAG.json (TAG defaults to POI): the scanned values, -2 Delta ln L at each,
and the best fit. r is kept at or above 0, as in the limits.

The scan runs in the Combine release, dispatched as run_combine_limits.py does.
"""

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_combine_limits as rcl  # noqa: E402  pylint: disable=wrong-import-position


def extract(root_file, poi, json_path):
    """Write the scan of a MultiDimFit grid output as JSON (run in the Combine release)"""
    # ROOT exists only in the Combine release, where this function runs
    import ROOT  # pylint: disable=import-outside-toplevel,import-error
    f = ROOT.TFile.Open(root_file)  # pylint: disable=no-member
    # the first entry is the best fit, the others are the grid points
    points = [(getattr(entry, poi), 2.0 * entry.deltaNLL) for entry in f.Get("limit")]
    x, dnll2 = zip(*sorted(points))
    out = {"poi": poi, "best": points[0][0], "x": list(x), "dnll2": list(dnll2)}
    Path(json_path).write_text(json.dumps(out, indent=1), encoding="utf-8")


def main(argv=None):
    """Run the scan of one card in the Combine release and write scan_TAG.json"""
    p = argparse.ArgumentParser(description=__doc__.split("\n", maxsplit=1)[0])
    p.add_argument("card")
    p.add_argument("--outdir", required=True)
    p.add_argument("--poi", default="r", help="r, or a nuisance (r then floats)")
    p.add_argument("--range", required=True,
                   help="LO,HI of the scanned parameter; write --range=LO,HI if LO is negative")
    p.add_argument("--points", type=int, default=40)
    p.add_argument("--expect-signal", type=float, default=0.0,
                   help="signal strength r of the Asimov dataset (default 0, background only)")
    p.add_argument("--freeze", default="",
                   help="parameters held at their pre-fit values, or allConstrainedNuisances")
    p.add_argument("--name", default="", help="TAG of the output file (default: POI)")
    rcl.add_release_options(p)
    p.add_argument("--extract", metavar="ROOT_FILE", help=argparse.SUPPRESS)
    args = p.parse_args(argv)

    outdir = Path(args.outdir).resolve()
    tag = args.name or args.poi
    if args.extract:
        extract(args.extract, args.poi, outdir / f"scan_{tag}.json")
        return 0

    outdir.mkdir(parents=True, exist_ok=True)
    card = str(Path(args.card).resolve())
    rmax = rcl.estimate_rmax(card) or 20.0
    ranges = f"{args.poi}={args.range}" + ("" if args.poi == "r" else f":r=0,{rmax:.6g}")
    scan = ["combine", "-M", "MultiDimFit", card, "-n", f"_scan_{tag}", "-m", "125",
            "--algo", "grid", "--points", str(args.points), "-P", args.poi,
            "--floatOtherPOIs", "1", "--setParameterRanges", ranges,
            "-t", "-1", "--expectSignal", f"{args.expect_signal:g}"]
    if args.freeze:
        scan += ["--freezeParameters", args.freeze]
    this = shlex.quote(str(Path(__file__).resolve()))
    inner = (f"cd {shlex.quote(str(outdir))} && {shlex.join(scan)} > scan_{tag}.log 2>&1 && "
             f"python3 {this} {shlex.quote(card)} --outdir {shlex.quote(str(outdir))} "
             f"--poi {shlex.quote(args.poi)} --range={shlex.quote(args.range)} "
             f"--name {shlex.quote(tag)} --extract higgsCombine_scan_{tag}.MultiDimFit.mH125.root")
    proc = subprocess.run(rcl.in_combine_release(inner, args), check=False)
    if proc.returncode or not (outdir / f"scan_{tag}.json").exists():
        sys.exit(f"the scan failed; see {outdir / f'scan_{tag}.log'}")
    print(f"wrote {outdir / f'scan_{tag}.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
