#!/usr/bin/env python3
"""Drive the corrections-study Condor campaign to completion

The campaign (sidm/studies/corrections/README.md, "The Condor campaign") runs
SidmProcessor over every v10 signal point and every skimmed background, simulation
only, with condor/submit_corrections.sub. This script does the bookkeeping between
submissions; run it from the repository root:

    python sidm/scripts/corrections_campaign.py status
    python sidm/scripts/corrections_campaign.py retry [--submit]
    python sidm/scripts/corrections_campaign.py skipped [--submit]
    python sidm/scripts/corrections_campaign.py merge [--samples a,b,...]

status   reconciles every job list against the Condor logs and the chunk outputs on
         EOS (condor/condor_campaign.py), and lists the samples that are complete.
retry    writes a job list of the chunks that failed, with the input files that are
         confirmed missing on EOS dropped, and submits it with --submit. Chunks
         that need a human (a deterministic python error) are listed, not retried.
skipped  finds the input files that finished chunks did not process (coffea runs
         with skipbadfiles=True, so an unreadable file is skipped, not fatal) and
         writes extra chunks of just those files, so a file that failed for a
         transient reason is processed in a later round. A file skipped in two
         rounds is reported as bad and left out.
merge    merges every complete sample that is not merged yet (or the --samples
         given) into MERGED_DIR, with a .meta.yaml sidecar per sample.

Repeat retry and skipped until nothing is left to resubmit, then merge.
"""

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CONDOR = os.path.join(REPO, "condor")
sys.path.insert(0, CONDOR)
import campaign_lib as cl  # noqa: E402  pylint: disable=wrong-import-position

REDIRECTOR = "root://cmseos.fnal.gov"
CHUNK_DIR = "/store/user/scampbel/sidm_condor/CorrectionsChunks_v1"
MERGED_DIR = "/store/user/scampbel/sidm_coffea_outputs/corrections_v1"
JOB_LISTS = ["job_args_corrections_2mu2e.txt", "job_args_corrections_4mu.txt",
             "job_args_corrections_background.txt"]
EXTRA_LIST = "job_args_corrections_extra.txt"   # chunks of skipped files, written by `skipped`
LOGS = os.path.join(CONDOR, "logs_corrections")
FILELISTS = os.path.join(CONDOR, "filelists_corrections")
STATE = os.path.join(CONDOR, "campaigns", "corrections_v1")
CHANNELS = "base,2mu2e,4mu,2mu2e_sr,4mu_sr,Z_to_mumu_control_region"
HIST_COLLECTIONS = "corrections_study"
EXTRA_CHUNK_OFFSET = 9000   # chunk numbers of the extra chunks: 9000, 9001, ...

# A sized listing of a few thousand chunks can take minutes when EOS is busy, and
# reconcile lists the chunk directory once per job list: list it once per command,
# with a generous timeout, and share the result.
_LISTINGS = {}
_xrdfs_ls_sizes = cl.xrdfs_ls_sizes


def _shared_listing(redirector, eos_dir, timeout=900):
    key = (redirector, eos_dir)
    if key not in _LISTINGS:
        _LISTINGS[key] = _xrdfs_ls_sizes(redirector, eos_dir, timeout=max(timeout, 900))
    return _LISTINGS[key]


cl.xrdfs_ls_sizes = _shared_listing


def job_lists():
    """The job lists of the campaign, with the extra chunks once there are any"""
    lists = [os.path.join(CONDOR, j) for j in JOB_LISTS]
    extra = os.path.join(CONDOR, EXTRA_LIST)
    return lists + ([extra] if os.path.exists(extra) else [])


def reconcile_all(stat_missing=False):
    """ChunkResults of every chunk of every job list"""
    results = []
    for path in job_lists():
        res, eos_ok, note = cl.reconcile(path, LOGS, CHUNK_DIR, REDIRECTOR, condor_dir=CONDOR,
                                         stat_missing=stat_missing)
        if not eos_ok:
            sys.exit(f"cannot list {CHUNK_DIR}: {note}")
        results += res
    return results


def per_sample(results):
    """{sample: counts of expected, done and in-progress chunks}"""
    samples = {}
    for r in results:
        s = samples.setdefault(r.sample, {"expected": 0, "done": 0, "in_progress": 0})
        s["expected"] += 1
        s["done"] += r.state == cl.DONE
        s["in_progress"] += r.state == cl.IN_PROGRESS
    return samples


def merged_samples():
    """Samples that already have a merged output in MERGED_DIR"""
    out = subprocess.run(["xrdfs", REDIRECTOR, "ls", MERGED_DIR], capture_output=True, text=True,
                         check=False)   # the directory does not exist before the first merge
    return {os.path.basename(p)[:-len(".coffea")] for p in out.stdout.split() if p.endswith(".coffea")}


def cmd_status(args):
    """Totals, completeness and every chunk that is neither done nor in progress"""
    results = reconcile_all(stat_missing=args.stat_missing)
    by_state = {}
    for r in results:
        by_state[r.state] = by_state.get(r.state, 0) + 1
    print(f"{len(results)} chunks: " + ", ".join(f"{k} {v}" for k, v in sorted(by_state.items())))
    samples = per_sample(results)
    complete = sorted(s for s, v in samples.items() if v["done"] == v["expected"])
    merged = merged_samples()
    print(f"{len(complete)}/{len(samples)} samples complete, {len(merged)} merged, "
          f"{len(set(complete) - merged)} complete but not merged")
    for r in results:
        if r.state not in (cl.DONE, cl.IN_PROGRESS):
            print(f"  {r.key:45s} {r.state:18s} {r.reason[:110]}")


def write_job_list(rows, name):
    """Write (sample, chunk, filelist) rows as a job list in condor/"""
    path = os.path.join(CONDOR, name)
    with open(path, "w", encoding="utf8") as f:
        for sample, chunk, filelist in rows:
            f.write(f"{sample} {chunk} {os.path.relpath(filelist, CONDOR)}\n")
    return path


def submit(job_list, memory=None):
    """condor_submit a job list with submit_corrections.sub"""
    # LPC's condor_submit is a wrapper script without a shebang: run it with sh
    cmd = ["sh", shutil.which("condor_submit") or "condor_submit", f"JOB_ARGS={os.path.basename(job_list)}"]
    if memory:
        cmd += ["-append", f"request_memory = {memory}"]
    cmd.append("submit_corrections.sub")
    print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=CONDOR, check=True)


def next_round():
    """Number of the next resubmission round"""
    os.makedirs(STATE, exist_ok=True)
    return len(glob.glob(os.path.join(CONDOR, "job_args_corrections_retry*.txt"))) + 1


FAILED_FILE = re.compile(r"Failed processing file: WorkItem\(dataset='[^']*', filename='([^']+)'")


def failing_file(result):
    """The input file a chunk failed on, if it failed on one file in every execution

    coffea names the file of the work item that raised. When Condor's own retries
    (max_retries in submit_corrections.sub) failed the same way, at least twice, the
    file is unreadable rather than unlucky: an OSError is skipped by skipbadfiles, so
    what remains is e.g. a corrupt basket (uproot DeserializationError).
    """
    log, err = result.log_paths.get("log"), result.log_paths.get("err")
    if not log or not err or not os.path.exists(err):
        return None
    with open(os.path.join(REPO, log) if not os.path.isabs(log) else log, encoding="utf8") as f:
        n_failed = len(re.findall(r"Normal termination \(return value [1-9]", f.read()))
    with open(os.path.join(REPO, err) if not os.path.isabs(err) else err, encoding="utf8",
              errors="replace") as f:
        files = set(FAILED_FILE.findall(f.read()))
    return files.pop() if n_failed >= 2 and len(files) == 1 else None


def recorded_bad_files():
    """/store paths of the files recorded as unreadable in bad_files.txt"""
    path = os.path.join(STATE, "bad_files.txt")
    if not os.path.exists(path):
        return set()
    with open(path, encoding="utf8") as f:
        return {cl.strip_redirector(line.split()[0]) for line in f if line.strip()}


def drop_files(filelist, bad):
    """Remove the bad files from a filelist in place; False if nothing is left"""
    files = [f for f in cl.read_filelist(filelist) if cl.strip_redirector(f) not in bad]
    with open(filelist, "w", encoding="utf8") as f:
        f.write("\n".join(files) + "\n")
    return bool(files)


def cmd_retry(args):
    """Job list of the failed chunks, with the inputs confirmed missing dropped"""
    results = reconcile_all(stat_missing=True)
    failed = [r for r in results if r.state not in (cl.DONE, cl.IN_PROGRESS)]
    specs = {}
    for path in job_lists():
        specs.update({s.key: s for s in cl.parse_job_args(path)})
    n = next_round()
    retry_dir = os.path.join(FILELISTS, f"retry{n}")
    rows, rows_oom, human = [], [], []
    for r in failed:
        if not r.auto_retryable:
            human.append(r)
            continue
        files = cl.read_filelist(os.path.join(CONDOR, specs[r.key].filelist))
        bad = set(r.offending_root_files) if r.state == cl.MISSING_ROOT_FILE else set()
        unreadable = failing_file(r)
        if unreadable:
            bad.add(unreadable)
            os.makedirs(STATE, exist_ok=True)
            with open(os.path.join(STATE, "bad_files.txt"), "a", encoding="utf8") as f:
                f.write(f"{unreadable}  # {r.key}: failed every attempt\n")
            print(f"  {r.key}: dropping {os.path.basename(unreadable)}, which failed every attempt")
        files = [f for f in files if f not in bad and cl.strip_redirector(f) not in
                 {cl.strip_redirector(b) for b in bad}]
        if not files:
            print(f"  {r.key}: every input file is missing on EOS; nothing to retry")
            continue
        os.makedirs(retry_dir, exist_ok=True)
        filelist = os.path.join(retry_dir, f"{r.key}.txt")
        with open(filelist, "w", encoding="utf8") as f:
            f.write("\n".join(files) + "\n")
        (rows_oom if r.state == cl.OOM else rows).append((r.sample, r.chunk, filelist))
        if bad - {unreadable}:
            print(f"  {r.key}: dropped {len(bad - {unreadable})} file(s) missing on EOS")
    for r in human:
        print(f"  NEEDS A HUMAN {r.key} {r.state}: {r.reason[:150]}")
    # every file dropped in an earlier round stays dropped
    known_bad = recorded_bad_files()
    rows = [(s, c, f) for s, c, f in rows if drop_files(f, known_bad)]
    rows_oom = [(s, c, f) for s, c, f in rows_oom if drop_files(f, known_bad)]
    if not rows and not rows_oom:
        print("nothing to retry")
        return
    for chunk_rows, suffix, memory in ((rows, "", None), (rows_oom, "_oom", "8GB")):
        if chunk_rows:
            path = write_job_list(chunk_rows, f"job_args_corrections_retry{n}{suffix}.txt")
            print(f"wrote {path} ({len(chunk_rows)} chunks)")
            if args.submit:
                submit(path, memory)


def processed_files(eos_path):
    """Input files a chunk output says it processed"""
    # only this command needs coffea
    from coffea.util import load   # pylint: disable=import-outside-toplevel
    out = load(cl_local(eos_path))
    return {item.filename for item in out.get("processed", ())}


def cl_local(eos_path):
    """A chunk output, through the EOS FUSE mount when there is one, else copied locally"""
    fuse = "/eos/uscms" + eos_path
    if os.path.exists(fuse):
        return fuse
    local = os.path.join(STATE, "tmp", os.path.basename(eos_path))
    os.makedirs(os.path.dirname(local), exist_ok=True)
    subprocess.run(["xrdcp", "-s", "-f", f"{REDIRECTOR}/{eos_path}", local], check=True)
    return local


def n_entries(url):
    """Events in a file, or None if it cannot be opened"""
    import uproot   # pylint: disable=import-outside-toplevel
    fuse = "/eos/uscms" + cl.strip_redirector(url)
    try:
        with uproot.open(fuse if os.path.exists(fuse) else url) as f:
            return f["Events"].num_entries
    except Exception:   # pylint: disable=broad-except
        return None


def cmd_skipped(args):
    """Extra chunks of the files that finished chunks skipped"""
    results = reconcile_all()
    specs = {}
    for path in job_lists():
        specs.update({s.key: s for s in cl.parse_job_args(path)})
    checked_path = os.path.join(STATE, "skipped_checked.json")
    checked = {}
    if os.path.exists(checked_path):
        with open(checked_path, encoding="utf8") as f:
            checked = json.load(f)
    for r in results:
        if r.state != cl.DONE or r.key in checked:
            continue
        files = cl.read_filelist(os.path.join(CONDOR, specs[r.key].filelist))
        done = processed_files(f"{CHUNK_DIR}/{r.key}.coffea")
        # compare /store paths: the filelists may name a different redirector
        done = {cl.strip_redirector(f) for f in done}
        checked[r.key] = [f for f in files if cl.strip_redirector(f) not in done]
        with open(checked_path, "w", encoding="utf8") as f:
            json.dump(checked, f, indent=1)
    # files dropped by `retry` as unreadable are not skipped, they are gone
    dropped = set()
    bad_path = os.path.join(STATE, "bad_files.txt")
    if os.path.exists(bad_path):
        with open(bad_path, encoding="utf8") as f:
            dropped = {cl.strip_redirector(line.split()[0]) for line in f if line.strip()}
    skipped = {k: [f for f in v if cl.strip_redirector(f) not in dropped] for k, v in checked.items()}
    skipped = {k: v for k, v in skipped.items() if v}
    n_skipped = sum(len(v) for v in skipped.values())
    print(f"{len(checked)} finished chunks checked; {n_skipped} input file(s) skipped in {len(skipped)} chunk(s)")

    # a file is retried once in an extra chunk; skipped again, it is bad
    extra_path = os.path.join(CONDOR, EXTRA_LIST)
    extra_specs = cl.parse_job_args(extra_path) if os.path.exists(extra_path) else []
    retried = set()
    for spec in extra_specs:
        retried.update(cl.read_filelist(os.path.join(CONDOR, spec.filelist)))
    bad = sorted({f for k, v in skipped.items() if int(k.rsplit("_", 1)[1]) >= EXTRA_CHUNK_OFFSET
                  for f in v})
    if bad:
        with open(os.path.join(STATE, "bad_files_skipped.txt"), "w", encoding="utf8") as f:
            f.write("\n".join(bad) + "\n")
        print(f"{len(bad)} file(s) skipped twice: listed as bad in {STATE}/bad_files_skipped.txt")
    todo = {}
    for key, files in skipped.items():
        sample = key.rsplit("_", 1)[0]
        todo.setdefault(sample, []).extend(f for f in files if f not in retried)
    # a file without events has no work item, so it never shows up as processed
    empty = [f for files in todo.values() for f in files if n_entries(f) == 0]
    if empty:
        with open(os.path.join(STATE, "empty_files.txt"), "w", encoding="utf8") as f:
            f.write("\n".join(sorted(empty)) + "\n")
        print(f"{len(empty)} of them have no events (listed in {STATE}/empty_files.txt)")
    todo = {s: [f for f in files if f not in empty] for s, files in todo.items()}
    todo = {s: f for s, f in todo.items() if f}
    if not todo:
        print("no skipped files left to retry")
        return
    next_chunk = {}
    for spec in extra_specs:
        next_chunk[spec.sample] = max(next_chunk.get(spec.sample, EXTRA_CHUNK_OFFSET), spec.chunk + 1)
    if not args.submit:
        for sample, files in sorted(todo.items()):
            print(f"  {sample}: {len(files)} skipped file(s) to retry")
        print("rerun with --submit to submit them as extra chunks")
        return
    rows = []
    os.makedirs(os.path.join(FILELISTS, "extra"), exist_ok=True)
    for sample, files in sorted(todo.items()):
        chunk = next_chunk.get(sample, EXTRA_CHUNK_OFFSET)
        filelist = os.path.join(FILELISTS, "extra", f"{sample}_{chunk}.txt")
        with open(filelist, "w", encoding="utf8") as f:
            f.write("\n".join(files) + "\n")
        rows.append((sample, chunk, filelist))
    with open(extra_path, "a", encoding="utf8") as f:
        for sample, chunk, filelist in rows:
            f.write(f"{sample} {chunk} {os.path.relpath(filelist, CONDOR)}\n")
    new_list = write_job_list(rows, f"job_args_corrections_extra_round{next_round()}.txt")
    print(f"{len(rows)} extra chunk(s) appended to {extra_path}; this round: {new_list}")
    submit(new_list)


def cmd_merge(args):
    """Merge the complete samples into MERGED_DIR"""
    results = reconcile_all()
    samples = per_sample(results)
    if args.samples:
        wanted = [s.strip() for s in args.samples.split(",") if s.strip()]
    else:
        merged = set() if args.remerge else merged_samples()
        wanted = sorted(s for s, v in samples.items()
                        if v["done"] == v["expected"] and s not in merged)
    incomplete = [s for s in wanted if samples[s]["done"] != samples[s]["expected"]]
    if incomplete:
        sys.exit(f"not every chunk of {incomplete} is done; refusing to merge them")
    if not wanted:
        print("nothing to merge")
        return
    os.makedirs(STATE, exist_ok=True)
    sample_file = os.path.join(STATE, "merge_samples.txt")
    with open(sample_file, "w", encoding="utf8") as f:
        f.write("\n".join(wanted) + "\n")
    # the sidecar lists the input files of every chunk of a sample, extra chunks
    # included, without the files dropped as unreadable
    bad = set()
    for name in ("bad_files.txt", "bad_files_skipped.txt"):
        path = os.path.join(STATE, name)
        if os.path.exists(path):
            with open(path, encoding="utf8") as f:
                bad.update(cl.strip_redirector(line.split()[0]) for line in f if line.strip())
    filelists = os.path.join(STATE, "filelists_all")
    os.makedirs(filelists, exist_ok=True)
    for path in glob.glob(os.path.join(FILELISTS, "**", "*.txt"), recursive=True):
        if "/retry" in path:
            continue    # a retry replaces its chunk's filelist, it does not add one
        files = [f for f in cl.read_filelist(path) if cl.strip_redirector(f) not in bad]
        with open(os.path.join(filelists, os.path.basename(path)), "w", encoding="utf8") as f:
            f.write("\n".join(files) + "\n")
    cmd = [sys.executable, os.path.join(REPO, "sidm", "scripts", "merge_coffea_chunks_eos.py"),
           "--input-eos-dir", CHUNK_DIR, "--output-eos-dir", MERGED_DIR,
           "--workdir", os.path.join(STATE, "merge_tmp"), "--samples", sample_file,
           "--filelists-dir", filelists, "--selections", CHANNELS,
           "--hist-collections", HIST_COLLECTIONS]
    print(f"merging {len(wanted)} sample(s)")
    subprocess.run(cmd, cwd=REPO, check=True)


def main():
    """Command-line entry point"""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", maxsplit=1)[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    st = sub.add_parser("status")
    st.add_argument("--stat-missing", action="store_true",
                    help="stat the inputs of failed chunks to tell missing files from timeouts")
    st.set_defaults(func=cmd_status)
    rt = sub.add_parser("retry")
    rt.add_argument("--submit", action="store_true")
    rt.set_defaults(func=cmd_retry)
    sk = sub.add_parser("skipped")
    sk.add_argument("--submit", action="store_true")
    sk.set_defaults(func=cmd_skipped)
    mg = sub.add_parser("merge")
    mg.add_argument("--samples", default=None)
    mg.add_argument("--remerge", action="store_true", help="also merge samples merged before")
    mg.set_defaults(func=cmd_merge)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
