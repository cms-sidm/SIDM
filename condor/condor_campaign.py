#!/usr/bin/env python3
"""SIDM Condor campaign orchestrator.

Phase 0 (this file): a read-only reconciler for a Condor campaign you submitted
with the existing condor/submit.sub. It compares the chunks you asked for
(job_args.txt) against the per-job Condor logs and the .coffea outputs on EOS,
classifies every failure by root cause, and writes one human report plus a
machine-readable report.json. It submits and resubmits nothing.

Run it on the LPC submitter node (the sidm_venv is fine, but a bare python3
works too -- this tool is standard-library only):

    python condor/condor_campaign.py reconcile \\
        --job-args     condor/job_args.txt \\
        --logs-dir     condor/logs \\
        --eos-chunk-dir /store/user/$USER/sidm_condor/SignalChunks_v1 \\
        --run-id       signals_v1

    python condor/condor_campaign.py status --report condor/campaigns/signals_v1/report.json

Later phases add init/submit/resubmit/run on top of campaign_lib; the read-only
verbs here never change.
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import campaign_lib as cl


def _default_outdir(run_id):
    return os.path.join("condor", "campaigns", run_id)


def _exit_code(report):
    """Exit contract (so 'reconcile && merge' is safe):
       0 = every chunk DONE; 1 = not all DONE (failures / in-flight / unknown);
       2 = EOS could not be listed; 3 = both.
    """
    t = report["totals"]
    code = 0
    if t["done"] != t["expected"]:
        code |= 1
    if not report["eos_listing_ok"]:
        code |= 2
    return code


def cmd_reconcile(args):
    redirector = args.redirector
    results, eos_ok, eos_note = cl.reconcile(
        job_args_path=args.job_args,
        logs_dir=args.logs_dir,
        eos_chunk_dir=args.eos_chunk_dir,
        redirector=redirector,
        condor_dir=args.condor_dir,
        stat_missing=not args.no_stat_missing,
    )
    inputs = {
        "job_args": os.path.abspath(args.job_args),
        "logs_dir": os.path.abspath(args.logs_dir),
        "eos_chunk_dir": cl.strip_redirector(args.eos_chunk_dir),
        "redirector": redirector,
    }
    report = cl.build_report(results, eos_ok, inputs, run_id=args.run_id, eos_note=eos_note)

    out_dir = args.out_dir or _default_outdir(args.run_id)
    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, "report.json")
    md_path = os.path.join(out_dir, "report.md")
    rendered = cl.render_report(report)
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2)
    with open(md_path, "w") as f:
        f.write(rendered + "\n")

    print(rendered)
    print(f"\nWrote {json_path}")
    print(f"Wrote {md_path}")
    return _exit_code(report)


def cmd_report(args):
    with open(args.report) as f:
        report = json.load(f)
    print(cl.render_report(report))
    return _exit_code(report)


def cmd_status(args):
    with open(args.report) as f:
        report = json.load(f)
    print(cl.render_status_line(report))
    return _exit_code(report)


def build_parser():
    p = argparse.ArgumentParser(
        prog="condor_campaign.py",
        description="Reconcile and report a SIDM Condor campaign (Phase 0: read-only).",
    )
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("reconcile",
                       help="classify every chunk and write report.{json,md}")
    r.add_argument("--job-args", required=True,
                   help="condor/job_args.txt that was submitted")
    r.add_argument("--logs-dir", required=True,
                   help="directory of per-job Condor .log/.out/.err files")
    r.add_argument("--eos-chunk-dir", required=True,
                   help="EOS dir the chunk .coffea outputs were written to "
                        "(a bare /store/... path; a full root:// URL is also accepted)")
    r.add_argument("--redirector", default="root://cmseos.fnal.gov",
                   help="EOS xrootd redirector (default: %(default)s)")
    r.add_argument("--run-id", default="manual",
                   help="label for this campaign (names the output dir)")
    r.add_argument("--out-dir", default=None,
                   help="where to write report.{json,md} (default: condor/campaigns/<run-id>)")
    r.add_argument("--condor-dir", default=None,
                   help="dir the job_args filelist paths are relative to "
                        "(default: the directory of --job-args)")
    r.add_argument("--no-stat-missing", action="store_true",
                   help="skip the per-file xrdfs stat that splits MISSING_ROOT_FILE "
                        "from EOS_TIMEOUT (faster; everything open-failed becomes "
                        "EOS_TIMEOUT, i.e. nothing is ever blamed on a missing file)")
    r.set_defaults(func=cmd_reconcile)

    rep = sub.add_parser("report", help="re-render a report.md from an existing report.json")
    rep.add_argument("--report", required=True, help="path to report.json")
    rep.set_defaults(func=cmd_report)

    st = sub.add_parser("status", help="print a one-line rollup from a report.json")
    st.add_argument("--report", required=True, help="path to report.json")
    st.set_defaults(func=cmd_status)

    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
