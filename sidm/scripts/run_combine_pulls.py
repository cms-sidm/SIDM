#!/usr/bin/env python3
"""Nuisance pulls of one datacard, from Combine's FitDiagnostics

    python sidm/scripts/run_combine_pulls.py CARD --outdir DIR [--expect-signal R] [--fit-data-obs]

Fits the expected background-only (Asimov) dataset, as a blinded analysis does, or
with --expect-signal the Asimov dataset with the signal at strength R, or with
--fit-data-obs the card's data_obs, and writes DIR/pulls.json: for every
constrained nuisance, its pre-fit value and uncertainty and its post-fit value and
uncertainty in the background-only fit (fit_b) and the signal-plus-background fit
(fit_s), plus the fitted r and the correlation matrix of every floating parameter
of each fit. On the Asimov dataset the pulls are zero by construction and the
post-fit uncertainties show how much the fit constrains each nuisance. The free
background parameters of the ABCD model (flatParam, named bkg_*) are left out of
the pulls, since they have no pre-fit constraint to be pulled against, but not out
of the correlations.

Every card these scripts write holds simulation as data_obs, so a fit to it is a
closure test of the background model and never a measurement. The fit runs in
the Combine release, dispatched as run_combine_limits.py does.
"""

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_combine_limits as rcl  # noqa: E402  pylint: disable=wrong-import-position


def extract(root_file, json_path):
    """Write the pulls and correlations of a fitDiagnostics file as JSON (run in the
    Combine release)"""
    # ROOT exists only in the Combine release, where this function runs
    import ROOT  # pylint: disable=import-outside-toplevel,import-error
    f = ROOT.TFile.Open(root_file)  # pylint: disable=no-member
    prefit = f.Get("nuisances_prefit")
    out = {"r": {}, "status": {}, "nuisances": {}, "correlation": {}}
    for fit in ("fit_b", "fit_s"):
        result = f.Get(fit)
        if not result:
            continue
        out["status"][fit] = result.status()
        names = [par.GetName() for par in result.floatParsFinal()]
        matrix = result.correlationMatrix()
        out["correlation"][fit] = {"names": names,
                                   "matrix": [[matrix[i][j] for j in range(len(names))]
                                              for i in range(len(names))]}
        for par in result.floatParsFinal():
            name = par.GetName()
            if name == "r":
                out["r"][fit] = [par.getVal(), par.getError()]
            elif not name.startswith("bkg_"):
                # nuisances_prefit keeps the values only; every constrained nuisance
                # these cards write (lnN, param) is a unit Gaussian
                before = prefit.find(name) if prefit else None
                entry = out["nuisances"].setdefault(
                    name, {"prefit": [before.getVal() if before else 0.0, 1.0]})
                entry[fit] = [par.getVal(), par.getError()]
    Path(json_path).write_text(json.dumps(out, indent=1), encoding="utf-8")


def main(argv=None):
    """Run FitDiagnostics on one card in the Combine release and write its pulls"""
    p = argparse.ArgumentParser(description=__doc__.split("\n", maxsplit=1)[0])
    p.add_argument("card")
    p.add_argument("--outdir", required=True)
    p.add_argument("--expect-signal", type=float, default=0.0,
                   help="signal strength r of the Asimov dataset (default 0, background only)")
    p.add_argument("--fit-data-obs", action="store_true",
                   help="fit the card's data_obs (simulation) instead of the Asimov dataset")
    rcl.add_release_options(p)
    p.add_argument("--extract", metavar="ROOT_FILE", help=argparse.SUPPRESS)
    args = p.parse_args(argv)

    outdir = Path(args.outdir).resolve()
    if args.extract:
        extract(args.extract, outdir / "pulls.json")
        return 0

    outdir.mkdir(parents=True, exist_ok=True)
    rmax = rcl.estimate_rmax(args.card) or 20.0
    fit = ["combine", "-M", "FitDiagnostics", str(Path(args.card).resolve()), "-n", "_pulls",
           "-m", "125", "--rMin", "0", "--rMax", f"{rmax:.6g}"]
    if not args.fit_data_obs:
        fit += ["-t", "-1", "--expectSignal", f"{args.expect_signal:g}"]
    this = shlex.quote(str(Path(__file__).resolve()))
    inner = (f"cd {shlex.quote(str(outdir))} && {shlex.join(fit)} > fit.log 2>&1 && "
             f"python3 {this} {shlex.quote(str(args.card))} --outdir {shlex.quote(str(outdir))} "
             f"--extract fitDiagnostics_pulls.root")
    proc = subprocess.run(rcl.in_combine_release(inner, args), check=False)
    if proc.returncode or not (outdir / "pulls.json").exists():
        sys.exit(f"FitDiagnostics failed; see {outdir / 'fit.log'}")
    print(f"wrote {outdir / 'pulls.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
