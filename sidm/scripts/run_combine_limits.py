#!/usr/bin/env python3
"""Run Combine's AsymptoticLimits over a directory of SIDM counting datacards.

The datacards are produced by ``sidm/studies/limit_plotting/make_datacards.ipynb``
(see ``sidm/studies/limit_plotting/datacard_tools.py``).  Each is a single-bin
counting experiment in one ABCD signal region for one signal point, with the
signal normalised to a 1 fb reference cross section -- so Combine's signal
strength ``r`` is directly the limit on the signal cross section in fb.

Combine lives in its own CMSSW release, which is generally not the release this
repository is checked out in, so each ``combine`` call is dispatched through a
sub-shell that sources ``cmsset_default.sh`` and runs ``scram runtime`` in the
Combine release.  Point ``--combine-cmssw`` (or ``$COMBINE_CMSSW_BASE``) at that
release.  If ``combine`` is already on ``$PATH``, pass ``--no-cmsenv`` to call it
directly.

Results are collected into ``limits.csv`` and ``limits.json`` under
``sidm/studies/limit_plotting/limits/`` (override with ``--outdir``), with one
row per datacard holding the observed limit and the -2/-1/median/+1/+2 sigma
expected limits, plus the signal point's grid coordinates parsed from its name.

Examples
--------
Expected (blinded) limits for every datacard, 8 at a time::

    python sidm/scripts/run_combine_limits.py -j 8

Re-run a subset without redoing the rest::

    python sidm/scripts/run_combine_limits.py --pattern 'datacard_SR_4mu_*' --force
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import getpass
import socket
import json
import math
import os
import re
import shlex
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
STUDY_DIR = REPO_ROOT / "sidm" / "studies" / "limit_plotting"
DEFAULT_DATACARD_DIR = STUDY_DIR / "datacards"
DEFAULT_OUTDIR = STUDY_DIR / "limits"
DEFAULT_COMBINE_CMSSW = "/uscms_data/d3/scampbel/CMSSW_14_1_0_pre4"
CMSSET = "/cvmfs/cms.cern.ch/cmsset_default.sh"

# Combine prints one line per quantile, e.g. "Expected 50.0%: r < 14.3125".
QUANTILE_LINE = re.compile(r"Expected\s+([\d.]+)%:\s*r\s*<\s*([\deE.+-]+)")
OBSERVED_LINE = re.compile(r"Observed Limit:\s*r\s*<\s*([\deE.+-]+)")
COMBINE_VERSION = re.compile(r"<<<\s*(v[\d.]+)\s*>>>")

# Map the quantiles Combine reports onto the column names used downstream.
QUANTILE_COLUMNS = {
    "2.5": "exp_m2",
    "16.0": "exp_m1",
    "50.0": "exp",
    "84.0": "exp_p1",
    "97.5": "exp_p2",
}

# The grid coordinates are recovered from the filename rather than imported from
# datacard_tools, so this script stays runnable in the Combine release, where
# coffea (which that module imports) is not installed.
# Combine version strings seen during this run, for the provenance sidecar.
COMBINE_VERSION_SEEN = set()

# datacard_<channel>_<signal>.txt, e.g. datacard_SR_2mu2e_2Mu2E_1000GeV_1p2GeV_0p96mm.txt
# datacard_<method>_<channel>_<signal>.txt, or datacard_comb_<signal>.txt for the
# combined cards, which span both channels and so carry no channel in the name.
CARD_NAME = re.compile(
    r"^datacard_(?:(?P<method>abcd|shape)_)?(?P<channel>SR_\w+?)_"
    r"(?P<signal>(?:2Mu2E|4Mu)_[\dp]+GeV_[\dp]+GeV_[\dp]+mm)$"
)
COMBINED_CARD_NAME = re.compile(
    r"^datacard_comb_(?:(?P<method>abcd)_)?"
    r"(?P<signal>Comb_[\dp]+GeV_[\dp]+GeV_[\dp]+mm)$"
)
SIGNAL_NAME = re.compile(
    r"^(?P<final_state>2Mu2E|4Mu|Comb)_(?P<mzd>[\dp]+)GeV_(?P<mdp>[\dp]+)GeV_(?P<ctau>[\dp]+)mm$"
)

CSV_COLUMNS = [
    "datacard", "method", "normalisation", "channel", "signal", "final_state",
    "m_mediator", "m_darkphoton", "ctau",
    "exp_m2", "exp_m1", "exp", "exp_p1", "exp_p2", "obs",
]


def parse_card_name(stem):
    """Pull the channel and signal point out of a datacard filename stem."""
    m = CARD_NAME.match(stem)
    combined = None
    if not m:
        combined = COMBINED_CARD_NAME.match(stem)
        if not combined:
            return {"channel": "", "signal": stem, "method": ""}
    if combined:
        # the combined cards carry both signal regions in one fit
        info = {"channel": "combined", "signal": combined["signal"],
                "method": "combined_abcd" if combined["method"] else "combined"}
        m = combined
    else:
        info = {"channel": m["channel"], "signal": m["signal"],
                "method": m["method"] or "counting"}
    s = SIGNAL_NAME.match(m["signal"])
    if s:
        num = lambda x: float(x.replace("p", "."))
        info.update(
            final_state=s["final_state"],
            m_mediator=num(s["mzd"]),
            m_darkphoton=num(s["mdp"]),
            ctau=num(s["ctau"]),
        )
    return info


def read_normalisation(card):
    """Whether a card normalises signal to 1 fb or to its theory cross section.

    Taken from the header comment the writer emits, so the results table says
    what the card actually did rather than what the output directory implies.
    """
    for line in Path(card).read_text().splitlines():
        if not line.startswith("#"):
            break
        if "theory cross section" in line:
            return "theory"
        if "normalised to" in line:
            return "1fb"
    return ""


def read_rates(card):
    """Return ``(signal_rate, background_rate)`` for a datacard's signal region.

    For the single-bin counting cards this is just the signal rate and the sum
    of the background rates.  The four-bin ABCD cards need different handling:
    their background ``rate`` column is 1 in every bin (the normalisation lives
    in the rateParams), so the background level has to come from the signal
    region's ``observation`` instead, and the signal rate has to be picked out
    of the ``_A`` bin rather than summed over all four.
    """
    obs_bins, observations, proc_bins, process_ids, rates = None, None, None, None, None
    header = {}
    for line in Path(card).read_text().splitlines():
        # The RooParametricHist shape cards have no yields in the rate column --
        # the signal rate is -1 (taken from its shape) and the background rate
        # is 1 (its normalisation lives in the workspace) -- so the writer
        # records the signal region levels in the header for this to read.
        if line.startswith("#") and "_rate_A" in line:
            key, _, value = line.lstrip("# ").partition("=")
            try:
                header[key.strip()] = float(value)
            except ValueError:
                pass
        fields = line.split()
        if not fields:
            continue
        head, rest = fields[0], fields[1:]
        if head == "bin":
            if obs_bins is None:
                obs_bins = rest          # first `bin` line pairs with observation
            else:
                proc_bins = rest         # second pairs with the process columns
        elif head == "observation":
            observations = [float(f) for f in rest]
        elif head == "process" and all(f.lstrip("-").isdigit() for f in rest):
            process_ids = [int(f) for f in rest]
        elif head == "rate":
            rates = [float(f) for f in rest]

    if "signal_rate_A" in header:
        return header["signal_rate_A"], header.get("background_rate_A", 0.0)

    if process_ids is None or rates is None or len(process_ids) != len(rates):
        return None, None

    # Single-bin card: no per-bin disambiguation needed.
    if not obs_bins or len(obs_bins) == 1:
        signal = sum(r for i, r in zip(process_ids, rates) if i <= 0)
        background = sum(r for i, r in zip(process_ids, rates) if i > 0)
        return signal, background

    # Multi-bin ABCD card: work with the signal-region bin only.
    sr_bin = next((b for b in obs_bins if b.endswith("_A")), obs_bins[0])
    background = None
    if observations and len(observations) == len(obs_bins):
        background = observations[obs_bins.index(sr_bin)]
    signal = sum(r for b, i, r in zip(proc_bins or [], process_ids, rates)
                 if b == sr_bin and i <= 0)
    return signal, background


def estimate_rmax(card, headroom=50.0, floor=20.0):
    """Guess an upper end for Combine's signal-strength scan.

    Combine's asymptotic scan only brackets the CLs crossing inside
    ``[rMin, rMax]``, and its default ``rMax`` of 20 is far too small for the
    faintest signal points here, where a few times 1e-4 expected events against
    ~30 background events puts the limit at ``r`` of order 1e4.  The 95% CL
    upper limit on the event count is roughly ``3 + 2*sqrt(B)``, so dividing
    that by the signal rate gives the right order of magnitude; ``headroom``
    then leaves a wide margin around it.
    """
    signal, background = read_rates(card)
    if not signal or signal <= 0:
        return None
    n_up = 3.0 + 2.0 * math.sqrt(max(background, 0.0))
    return max(floor, headroom * n_up / signal)


def build_command(card, workdir, args, rmax=None):
    """Assemble the shell command that runs Combine on one datacard.

    Combine writes ``higgsCombine*.root`` into its working directory and has no
    option to redirect that, so every card is run in its own scratch directory
    to keep parallel jobs from colliding.
    """
    combine = [
        args.combine, "-M", "AsymptoticLimits", str(Path(card).resolve()),
        "-n", f"_{Path(card).stem}", "-m", "125",
    ]
    if args.blind:
        combine += ["--run", "blind"]
    if rmax:
        combine += ["--rMin", "0", "--rMax", f"{rmax:.6g}"]
    if args.extra:
        combine += shlex.split(args.extra)

    inner = f"cd {shlex.quote(str(workdir))} && {shlex.join(combine)}"
    if args.no_cmsenv:
        return ["bash", "-c", inner]

    # `scram runtime` must run from inside the release, and any CMSSW state
    # inherited from the calling shell has to be dropped first or the wrong
    # release wins.
    setup = (
        "unset CMSSW_BASE CMSSW_RELEASE_BASE SCRAM_ARCH PYTHONPATH PYTHONHOME "
        "PYTHONUSERBASE; "
        f"source {shlex.quote(CMSSET)} >/dev/null 2>&1 && "
        f"cd {shlex.quote(str(Path(args.combine_cmssw) / 'src'))} && "
        'eval "$(scramv1 runtime -sh)" && '
    )
    return ["bash", "-c", setup + inner]


def parse_output(text):
    """Extract the limit quantiles from Combine's stdout."""
    result = {}
    for quantile, value in QUANTILE_LINE.findall(text):
        column = QUANTILE_COLUMNS.get(quantile)
        if column:
            result[column] = float(value)
    observed = OBSERVED_LINE.search(text)
    if observed:
        result["obs"] = float(observed.group(1))
    return result


def run_one(card, args):
    """Run Combine on a single datacard and return its parsed limits.

    A scan that fails to bracket the CLs crossing exits cleanly but prints no
    limit, so a miss is retried once with a hundredfold wider range before it
    is reported as a failure.

    Returns ``(row, error)``; exactly one of the two is ``None``.
    """
    card = Path(card)
    if args.rmax:
        rmax = args.rmax
    elif args.auto_rmax:
        rmax = estimate_rmax(card)
    else:
        rmax = None

    output, limits = "", {}
    for attempt in range(2):
        with tempfile.TemporaryDirectory(prefix="combine_", dir=args.scratch) as workdir:
            cmd = build_command(card, workdir, args, rmax=rmax)
            try:
                proc = subprocess.run(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, timeout=args.timeout,
                )
            except subprocess.TimeoutExpired:
                return None, f"{card.name}: timed out after {args.timeout}s"
            output = proc.stdout

        if proc.returncode != 0:
            break
        limits = parse_output(output)
        if "exp" in limits:
            m = COMBINE_VERSION.search(output)
            if m:
                COMBINE_VERSION_SEEN.add(m.group(1))
            break
        if rmax is None or args.rmax:
            break  # a fixed range was requested; do not second-guess it
        rmax *= 100

    if args.logdir:
        Path(args.logdir).mkdir(parents=True, exist_ok=True)
        (Path(args.logdir) / f"{card.stem}.log").write_text(output)

    if proc.returncode != 0:
        tail = "\n".join(output.strip().splitlines()[-15:])
        return None, f"{card.name}: combine exited {proc.returncode}\n{tail}"
    # A card with no expected median is a scan that never found the crossing.
    if "exp" not in limits:
        tail = "\n".join(output.strip().splitlines()[-15:])
        return None, f"{card.name}: no limit found in combine output\n{tail}"

    row = {"datacard": card.name, "normalisation": read_normalisation(card)}
    row.update(parse_card_name(card.stem))
    row.update(limits)
    return row, None


def load_existing(csv_path):
    """Read an existing limits.csv so completed points can be skipped."""
    if not csv_path.exists():
        return {}
    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    out = {}
    for row in rows:
        for key, value in list(row.items()):
            if key in {"datacard", "method", "normalisation", "channel",
                       "signal", "final_state"}:
                continue
            row[key] = float(value) if value not in ("", None) else None
        out[row["datacard"]] = row
    return out


def write_limit_metadata(outdir, args, rows, errors, cards):
    """Write ``limits.meta.yaml`` next to ``limits.csv``.

    Mirrors the ``.meta.yaml`` convention the merged coffea inputs use, so a set
    of limits records the conditions it was produced under and can be sorted by
    them. The datacard directory's own sidecar is embedded whole, which in turn
    carries the upstream coffea provenance -- so the chain
    coffea -> datacards -> limits is recoverable from this one file.
    """
    datacard_meta = None
    sidecar_in = Path(args.datacards) / "datacards.meta.yaml"
    if sidecar_in.exists():
        try:
            with open(sidecar_in) as f:
                datacard_meta = yaml.safe_load(f)
        except Exception as exc:
            datacard_meta = {"error": f"could not read {sidecar_in}: {exc}"}

    exp = [r["exp"] for r in rows if r.get("exp") is not None]
    obs = [r["obs"] for r in rows if r.get("obs") is not None]
    record = {
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "produced_by": "sidm/scripts/run_combine_limits.py",
        "producer": getpass.getuser(),
        "host": socket.gethostname(),
        "combine": {
            "method": "AsymptoticLimits",
            "version": sorted(COMBINE_VERSION_SEEN) or None,
            "cmssw_release": None if args.no_cmsenv else str(args.combine_cmssw),
            "blinded": bool(args.blind),
            "expected_kind": "a-priori (Asimov, observation ignored)" if args.blind
                             else "a-posteriori (fit to the observation)",
            "observed_limits_present": bool(obs),
            "rmax_strategy": ("fixed" if args.rmax else
                              "per-card estimate from S and B, one 100x retry"
                              if args.auto_rmax else "combine default"),
            "rmax_fixed": args.rmax,
            "extra_args": args.extra or None,
            "timeout_s": args.timeout,
            "jobs": args.jobs,
        },
        "run": {
            "datacard_dir": str(args.datacards),
            "pattern": args.pattern,
            "n_datacards_matched": len(cards),
            "n_limits": len(rows),
            "n_failed": len(errors),
            "failures": [e.splitlines()[0] for e in errors] or None,
        },
        "results_summary": {
            "median_expected_r": round(sorted(exp)[len(exp) // 2], 6) if exp else None,
            "min_expected_r": round(min(exp), 6) if exp else None,
            "max_expected_r": round(max(exp), 6) if exp else None,
            "median_observed_r": round(sorted(obs)[len(obs) // 2], 6) if obs else None,
        },
        "warning": ("Observed limits here are the fit run against simulation standing in "
                    "for data -- a closure test, not a measurement.") if obs else None,
        "datacards": datacard_meta,
    }
    path = Path(outdir) / "limits.meta.yaml"
    with open(path, "w") as f:
        yaml.safe_dump(record, f, sort_keys=False, default_flow_style=False)
    return path


def write_results(rows, outdir):
    """Write limits.csv and limits.json, sorted for stable diffs."""
    outdir.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda r: (r.get("normalisation", ""), r.get("method", ""),
                                       r.get("channel", ""), r.get("signal", "")))

    csv_path = outdir / "limits.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    json_path = outdir / "limits.json"
    json_path.write_text(json.dumps(rows, indent=2) + "\n")
    return csv_path, json_path


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--datacards", default=str(DEFAULT_DATACARD_DIR),
                   help="directory of datacard .txt files (default: %(default)s)")
    p.add_argument("--pattern", default="datacard_*.txt",
                   help="glob selecting datacards within that directory "
                        "(default: %(default)s)")
    p.add_argument("--outdir", default=str(DEFAULT_OUTDIR),
                   help="where limits.csv / limits.json are written "
                        "(default: %(default)s)")
    p.add_argument("--logdir", default=None,
                   help="if set, save each combine log here")
    p.add_argument("--combine-cmssw", default=os.environ.get(
                       "COMBINE_CMSSW_BASE", DEFAULT_COMBINE_CMSSW),
                   help="CMSSW release holding HiggsAnalysis/CombinedLimit "
                        "(default: $COMBINE_CMSSW_BASE or %(default)s)")
    p.add_argument("--combine", default="combine",
                   help="combine executable name (default: %(default)s)")
    p.add_argument("--no-cmsenv", action="store_true",
                   help="call combine directly instead of sourcing a CMSSW env")
    p.add_argument("--blind", action="store_true", default=True,
                   help="expected limits only, from the Asimov dataset (default)")
    p.add_argument("--unblind", dest="blind", action="store_false",
                   help="also compute the observed limit from the card's "
                        "observation line")
    p.add_argument("--rmax", type=float, default=None,
                   help="fixed --rMax for every combine call; by default the "
                        "range is estimated per datacard from its S and B")
    p.add_argument("--no-auto-rmax", dest="auto_rmax", action="store_false",
                   default=True,
                   help="leave the signal-strength range at combine's default")
    p.add_argument("--extra", default="",
                   help="extra arguments appended to every combine call")
    p.add_argument("-j", "--jobs", type=int, default=4,
                   help="datacards to process in parallel (default: %(default)s)")
    p.add_argument("--timeout", type=int, default=600,
                   help="per-datacard timeout in seconds (default: %(default)s)")
    p.add_argument("--scratch", default=None,
                   help="parent directory for per-job scratch dirs "
                        "(default: system temp)")
    p.add_argument("--force", action="store_true",
                   help="re-run datacards already present in limits.csv")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    datacard_dir = Path(args.datacards)
    cards = sorted(datacard_dir.glob(args.pattern))
    if not cards:
        sys.exit(f"no datacards matching {args.pattern!r} under {datacard_dir}")

    if not args.no_cmsenv:
        release = Path(args.combine_cmssw)
        if not (release / "src" / "HiggsAnalysis" / "CombinedLimit").is_dir():
            sys.exit(
                f"no HiggsAnalysis/CombinedLimit under {release}/src.\n"
                "Set --combine-cmssw (or $COMBINE_CMSSW_BASE) to a release with "
                "Combine installed, or pass --no-cmsenv if combine is already "
                "on $PATH."
            )

    outdir = Path(args.outdir)
    # Existing rows are always loaded, so that running with a --pattern that
    # selects part of the grid rewrites those points without discarding the
    # rest of limits.csv.  --force only drops the rows being re-run.
    existing = load_existing(outdir / "limits.csv")
    if args.force:
        selected = {c.name for c in cards}
        existing = {name: row for name, row in existing.items() if name not in selected}
    todo = [c for c in cards if c.name not in existing]
    if existing:
        print(f"keeping {len(existing)} limits already in limits.csv, "
              f"{len(todo)} to run" + ("" if args.force else " (use --force to redo them)"))

    print(f"running combine on {len(todo)} datacards with {args.jobs} jobs")
    rows, errors = list(existing.values()), []
    if todo:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(run_one, c, args): c for c in todo}
            for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
                row, error = future.result()
                if error:
                    errors.append(error)
                    print(f"[{i}/{len(todo)}] FAILED {error.splitlines()[0]}")
                else:
                    rows.append(row)
                    print(f"[{i}/{len(todo)}] {row['datacard']}: "
                          f"expected r < {row['exp']:.4g}")

    if not rows:
        sys.exit("no limits were produced")

    csv_path, json_path = write_results(rows, outdir)
    meta_path = write_limit_metadata(outdir, args, rows, errors, cards)
    print(f"\nwrote {len(rows)} limits to {csv_path} and {json_path}")
    print(f"provenance sidecar: {meta_path}")
    if errors:
        print(f"\n{len(errors)} datacards failed:")
        for error in errors:
            print(f"  {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
