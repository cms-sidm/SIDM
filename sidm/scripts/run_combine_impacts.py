#!/usr/bin/env python3
"""Nuisance impacts of one datacard, with combineTool.py -M Impacts

    python sidm/scripts/run_combine_impacts.py CARD --outdir DIR --expect-signal R [-j 4]

Fits the Asimov dataset with the signal injected at strength R and, for every
constrained nuisance, refits with that nuisance fixed one post-fit standard
deviation above and below its best-fit value: the shift of r is the nuisance's
impact. Writes DIR/impacts.json in combineTool's format: for every nuisance its best
fit with its -1 and +1 sigma values ("fit"), the same before the fit ("prefit"), and
r at those three values ("r"), plus the fit of r ("POIs"). The signal has to be
injected: at r = 0 a nuisance that acts on the signal alone moves r by nothing. The
free background parameters of the ABCD model (bkg_*) are left out, as in
run_combine_pulls.py. r is kept at or above 0, as in the limits.

The fits run in the Combine release, dispatched as run_combine_limits.py does.
"""

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_combine_limits as rcl  # noqa: E402  pylint: disable=wrong-import-position


def main(argv=None):
    """Run the impact fits of one card in the Combine release and write impacts.json"""
    p = argparse.ArgumentParser(description=__doc__.split("\n", maxsplit=1)[0])
    p.add_argument("card")
    p.add_argument("--outdir", required=True)
    p.add_argument("--expect-signal", type=float, required=True,
                   help="signal strength r injected into the Asimov dataset")
    p.add_argument("-j", "--jobs", type=int, default=4, help="fits run in parallel")
    rcl.add_release_options(p)
    args = p.parse_args(argv)

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    card = Path(args.card).resolve()
    rmax = rcl.estimate_rmax(card) or 20.0
    # combineTool needs the same options in its three steps. The tolerances are 100
    # times tighter than Combine's defaults, which leave r uncertain by about 1% of its
    # own uncertainty, as large as the impact of a 1% signal uncertainty.
    common = ["-M", "Impacts", "-d", "workspace.root", "-m", "125", "-t", "-1",
              "--expectSignal", f"{args.expect_signal:g}", "--rMin", "0", "--rMax", f"{rmax:.6g}",
              "--robustFit", "1", "--cminDefaultMinimizerTolerance", "0.001",
              "--setRobustFitTolerance", "0.001", "--exclude", "rgx{bkg_.*}"]
    steps = [["text2workspace.py", str(card), "-o", "workspace.root", "-m", "125"],
             ["combineTool.py", *common, "--doInitialFit"],
             ["combineTool.py", *common, "--doFits", "--parallel", str(args.jobs)],
             ["combineTool.py", *common, "-o", "impacts.json"]]
    inner = f"cd {shlex.quote(str(outdir))} && " + " && ".join(
        f"{shlex.join(step)} >> impacts.log 2>&1" for step in steps)
    proc = subprocess.run(rcl.in_combine_release(inner, args), check=False)
    if proc.returncode or not (outdir / "impacts.json").exists():
        sys.exit(f"the impact fits failed; see {outdir / 'impacts.log'}")
    print(f"wrote {outdir / 'impacts.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
