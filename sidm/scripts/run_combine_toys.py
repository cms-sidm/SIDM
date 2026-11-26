#!/usr/bin/env python3
"""Signal-injection test of one datacard: fits of toy datasets with a known signal

    python sidm/scripts/run_combine_toys.py CARD --outdir DIR --expect-signal R [--toys 500] [-j 8]

Generates toy datasets from the model with the signal at strength R and fits each
with FitDiagnostics. The toys are frequentist ones made from the pre-fit model
(--toysFrequentist --bypassFrequentistFit): the global observables of the nuisances
are randomised with each toy, and the card's data_obs is not used. Writes
DIR/toys.json: the injected R and, for every toy, the fitted r with its lower and
upper (MINOS) uncertainty, the fit status, and the fitted value of every constrained
nuisance. For an unbiased fit with the right uncertainties, (r - R) divided by the
uncertainty on the side of R is distributed as a unit Gaussian, and so is every
nuisance that the data do not constrain. r is kept at or above 0, as in the limits,
so toys that fluctuate below 0 stop there.

The toys run in the Combine release, in parallel jobs with their own seeds,
dispatched as run_combine_limits.py does.
"""

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_combine_limits as rcl  # noqa: E402  pylint: disable=wrong-import-position


def extract(files, injected, json_path):
    """Write the fits of FitDiagnostics toy outputs as JSON (run in the Combine release)"""
    # ROOT exists only in the Combine release, where this function runs
    import ROOT  # pylint: disable=import-outside-toplevel,import-error
    out = {"injected": injected, "r": [], "r_lo": [], "r_hi": [], "status": [], "nuisances": {}}
    for path in files:
        f = ROOT.TFile.Open(path)  # pylint: disable=no-member
        tree = f.Get("tree_fit_sb") if f and not f.IsZombie() else None
        if not tree:
            continue
        branches = {b.GetName() for b in tree.GetListOfBranches()}
        # the constrained nuisances are the parameters with a global observable, NAME_In
        nuisances = sorted(n for n in branches if f"{n}_In" in branches)
        for entry in tree:
            out["r"].append(entry.r)
            out["r_lo"].append(entry.rLoErr)
            out["r_hi"].append(entry.rHiErr)
            out["status"].append(entry.fit_status)
            for name in nuisances:
                out["nuisances"].setdefault(name, []).append(getattr(entry, name))
    Path(json_path).write_text(json.dumps(out), encoding="utf-8")


def main(argv=None):
    """Fit toys of one card in the Combine release and write toys.json"""
    p = argparse.ArgumentParser(description=__doc__.split("\n", maxsplit=1)[0])
    p.add_argument("card")
    p.add_argument("--outdir", required=True)
    p.add_argument("--expect-signal", type=float, required=True,
                   help="signal strength r the toys are generated with")
    p.add_argument("--toys", type=int, default=500)
    p.add_argument("-j", "--jobs", type=int, default=8, help="parallel jobs, one seed each")
    rcl.add_release_options(p)
    p.add_argument("--extract", nargs="+", metavar="ROOT_FILE", help=argparse.SUPPRESS)
    args = p.parse_args(argv)

    outdir = Path(args.outdir).resolve()
    if args.extract:
        extract(args.extract, args.expect_signal, outdir / "toys.json")
        return 0

    outdir.mkdir(parents=True, exist_ok=True)
    card = str(Path(args.card).resolve())
    rmax = rcl.estimate_rmax(card) or 20.0
    per_job = -(-args.toys // args.jobs)
    seeds = range(1, args.jobs + 1)
    jobs = []
    for seed in seeds:
        fit = ["combine", "-M", "FitDiagnostics", card, "-n", f"_toys_s{seed}", "-m", "125",
               "-s", str(seed), "-t", str(per_job), "--expectSignal", f"{args.expect_signal:g}",
               "--toysFrequentist", "--bypassFrequentistFit", "--skipBOnlyFit",
               "--rMin", "0", "--rMax", f"{rmax:.6g}"]
        jobs.append(f"{shlex.join(fit)} > toys_s{seed}.log 2>&1 &")
    files = " ".join(f"fitDiagnostics_toys_s{seed}.root" for seed in seeds)
    this = shlex.quote(str(Path(__file__).resolve()))
    # "cd DIR && job &" would put the cd, and the release set-up before it, into the
    # first job's background shell only; ending the chain with ";" keeps both here
    inner = (f"cd {shlex.quote(str(outdir))} || exit 1; {' '.join(jobs)} wait; "
             f"python3 {this} {shlex.quote(card)} --outdir {shlex.quote(str(outdir))} "
             f"--expect-signal {args.expect_signal:g} --extract {files}")
    proc = subprocess.run(rcl.in_combine_release(inner, args), check=False)
    if proc.returncode or not (outdir / "toys.json").exists():
        sys.exit(f"the toy fits failed; see the toys_s*.log files in {outdir}")
    print(f"wrote {outdir / 'toys.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
