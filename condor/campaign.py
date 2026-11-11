#!/usr/bin/env python3
"""Drive a SIDM Condor campaign from one YAML config

A campaign is one file in condor/campaign_configs/: the processor options, the job
lists and how to make them, where the jobs write their chunk outputs and where the
merged samples go (template.yaml documents every key). A new study writes a config
instead of a script. Run from the repository root:

    python condor/campaign.py CONFIG jobs [--force]
    python condor/campaign.py CONFIG submit [--job-list FILE]
    python condor/campaign.py CONFIG status [--stat-missing]
    python condor/campaign.py CONFIG retry [--submit]
    python condor/campaign.py CONFIG skipped [--submit]
    python condor/campaign.py CONFIG merge [--samples a,b,...] [--remerge]

CONFIG is a path, or the name of a file in condor/campaign_configs/ without .yaml.

jobs     writes one job list per entry of job_lists with condor/make_job_args.py, with
         at most max_files_per_job input files per job (100 unless the config says
         otherwise: a job of a few hundred tiny files spends hours on the fixed cost
         per file). Existing job lists are kept unless --force is given.
submit   submits the job lists, or the one given, with condor/submit_campaign.sub. A
         job that runs longer than stall_hours is held, released onto another node up
         to stall_restarts times, and then removed, which status reports as STALLED
         and retry resubmits.
status   reconciles every job list against the Condor logs and the chunk outputs on
         EOS (condor/campaign_lib.py, as condor_campaign.py does), and lists the
         samples that are complete.
retry    writes a job list of the chunks that failed, with the input files that are
         confirmed missing on EOS dropped, and submits it with --submit; chunks killed
         for memory get oom_memory. A file that failed in every execution of its job,
         at least twice, is dropped and recorded in bad_files.txt. Chunks that need a
         human (a deterministic python error, a manual condor_rm) are listed, not
         retried.
skipped  finds the input files that finished chunks did not process (coffea runs
         with skipbadfiles=True, so an unreadable file is skipped, not fatal) and
         writes extra chunks of just those files, so a file that failed for a
         transient reason is processed in a later round. A file skipped in two
         rounds is reported as bad and left out.
merge    merges every complete sample that is not merged yet (or the --samples
         given) into merged_dir, with a .meta.yaml sidecar per sample.

Repeat retry and skipped until nothing is left to resubmit, then merge. The working
state (the bad-file lists, the files skipped per chunk, merge scratch space) is kept
in condor/campaigns/<name>/, which is not committed.
"""

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys

import yaml

CONDOR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(CONDOR)
CONFIG_DIR = os.path.join(CONDOR, "campaign_configs")
SUBMIT_FILE = "submit_campaign.sub"
sys.path.insert(0, CONDOR)
import campaign_lib as cl  # noqa: E402  pylint: disable=wrong-import-position

HIST_MODES = ("unweighted", "weighted", "syst")
EXTRA_CHUNK_OFFSET = 9000   # chunk numbers of the extra chunks made by `skipped`: 9000, 9001, ...
FAILURES_TO_DROP = 2        # a file that failed this many executions of its job is dropped

REQUIRED = ("name", "chunk_dir", "merged_dir", "channels", "hist_collections", "job_lists")
DEFAULTS = {
    "redirector": "root://cmseos.fnal.gov",
    "corrections": "default",
    "hist_mode": "unweighted",
    "monitor_corrections": [],
    "logs_dir": None,          # logs_<name>
    "filelists_dir": None,     # filelists_<name>
    "job_list_prefix": None,   # job_args_<name>
    "max_files_per_job": 100,
    "request_memory": "2GB",
    "oom_memory": "8GB",
    "max_retries": 2,
    "stall_hours": 12,
    "stall_restarts": 1,
}
JOB_LIST_REQUIRED = ("job_args", "samples_file", "dataset", "location_cfg")
JOB_LIST_DEFAULTS = {"files_per_job": 1, "files_per_job_map": None, "replace_xcache": False,
                     "census_skip": None, "max_files": None}


def _names(value, key):
    """A list of names from a YAML list or a comma-separated string"""
    if isinstance(value, str):
        value = [v.strip() for v in value.split(",")]
    if not isinstance(value, list) or not all(isinstance(v, str) and v for v in value):
        raise ValueError(f"{key} must be a list of names, got {value!r}")
    return value


class Campaign:
    """One campaign config, with its defaults filled in and its paths resolved"""

    def __init__(self, cfg, path="<config>"):
        if not isinstance(cfg, dict):
            raise ValueError(f"{path}: a campaign config is a mapping of keys to values")
        missing = [k for k in REQUIRED if k not in cfg]
        unknown = sorted(set(cfg) - set(REQUIRED) - set(DEFAULTS))
        if missing or unknown:
            raise ValueError(f"{path}: missing keys {missing}, unknown keys {unknown}")
        c = {**DEFAULTS, **cfg}
        self.path = path
        self.name = str(c["name"])
        self.redirector = c["redirector"].rstrip("/")
        # $USER and other environment variables are expanded, so one config serves everyone
        self.chunk_dir = cl.strip_redirector(os.path.expandvars(c["chunk_dir"])).rstrip("/")
        self.merged_dir = cl.strip_redirector(os.path.expandvars(c["merged_dir"])).rstrip("/")
        for key, value in (("chunk_dir", self.chunk_dir), ("merged_dir", self.merged_dir)):
            if not value.startswith("/store/") or "$" in value:
                raise ValueError(f"{path}: {key} must be an EOS path under /store/, got {value!r}")
        self.channels = _names(c["channels"], "channels")
        self.hist_collections = _names(c["hist_collections"], "hist_collections")
        if c["corrections"] in ("default", "none", None):
            self.corrections = c["corrections"] or "none"
        else:
            self.corrections = ",".join(_names(c["corrections"], "corrections"))
        self.hist_mode = c["hist_mode"]
        if self.hist_mode not in HIST_MODES:
            raise ValueError(f"{path}: hist_mode must be one of {HIST_MODES}, got {self.hist_mode!r}")
        self.monitor = _names(c["monitor_corrections"], "monitor_corrections") if c["monitor_corrections"] else []
        if self.corrections == "none" and (self.monitor or self.hist_mode == "syst"):
            raise ValueError(f"{path}: monitor_corrections and hist_mode syst need corrections")
        self.logs_dir = os.path.join(CONDOR, c["logs_dir"] or f"logs_{self.name}")
        self.filelists_dir = os.path.join(CONDOR, c["filelists_dir"] or f"filelists_{self.name}")
        self.prefix = c["job_list_prefix"] or f"job_args_{self.name}"
        self.state_dir = os.path.join(CONDOR, "campaigns", self.name)
        self.max_files_per_job = c["max_files_per_job"]
        if self.max_files_per_job is not None and (not isinstance(self.max_files_per_job, int)
                                                   or self.max_files_per_job < 1):
            raise ValueError(f"{path}: max_files_per_job must be a positive integer or null")
        self.request_memory, self.oom_memory = str(c["request_memory"]), str(c["oom_memory"])
        self.max_retries, self.stall_hours = int(c["max_retries"]), float(c["stall_hours"])
        self.stall_restarts = int(c["stall_restarts"])
        if self.max_retries < 0 or self.stall_hours <= 0 or self.stall_restarts < 0:
            raise ValueError(f"{path}: max_retries and stall_restarts must be >= 0, stall_hours > 0")
        if not isinstance(c["job_lists"], list) or not c["job_lists"]:
            raise ValueError(f"{path}: job_lists must be a non-empty list")
        self.job_list_specs = []
        for i, spec in enumerate(c["job_lists"]):
            missing = [k for k in JOB_LIST_REQUIRED if k not in (spec or {})]
            unknown = sorted(set(spec or {}) - set(JOB_LIST_REQUIRED) - set(JOB_LIST_DEFAULTS))
            if missing or unknown:
                raise ValueError(f"{path}: job_lists[{i}]: missing keys {missing}, unknown keys {unknown}")
            self.job_list_specs.append({**JOB_LIST_DEFAULTS, **spec})
        names = [s["job_args"] for s in self.job_list_specs]
        if len(set(names)) != len(names):
            raise ValueError(f"{path}: two job lists have the same job_args file")

    @property
    def chunk_url(self):
        """Where the jobs write their chunk outputs, as an xrootd URL"""
        return f"{self.redirector}/{self.chunk_dir}"

    @property
    def extra_list(self):
        """Job list of the extra chunks of skipped files, which `skipped` appends to"""
        return os.path.join(CONDOR, f"{self.prefix}_extra.txt")

    def job_lists(self):
        """Paths of the campaign's job lists, with the extra chunks once there are any"""
        lists = [os.path.join(CONDOR, s["job_args"]) for s in self.job_list_specs]
        return lists + ([self.extra_list] if os.path.exists(self.extra_list) else [])

    def submit_macros(self, job_list, memory=None):
        """The values condor/submit_campaign.sub takes for one job list"""
        return {
            "JOB_ARGS": os.path.relpath(job_list, CONDOR),
            "CHUNK_URL": self.chunk_url,
            "CHANNELS": ",".join(self.channels),
            "HIST_COLLECTIONS": ",".join(self.hist_collections),
            "CORRECTIONS": self.corrections,
            "HIST_MODE": self.hist_mode,
            "MONITOR": ",".join(self.monitor),
            "LOGS": os.path.relpath(self.logs_dir, CONDOR),
            "MEMORY": memory or self.request_memory,
            "MAX_RETRIES": str(self.max_retries),
            "STALL_HOURS": f"{self.stall_hours:g}",
            "STALL_RESTARTS": str(self.stall_restarts),
        }

    def submit_command(self, job_list, memory=None):
        """condor_submit command line for one job list, run from condor/"""
        # LPC's condor_submit is a wrapper script without a shebang: run it with sh
        cmd = ["sh", shutil.which("condor_submit") or "condor_submit"]
        cmd += [f"{k}={v}" for k, v in self.submit_macros(job_list, memory).items()]
        return cmd + [SUBMIT_FILE]

    def make_job_args_command(self, spec):
        """condor/make_job_args.py command line for one entry of job_lists, run from the repo"""
        cmd = [sys.executable, os.path.join("condor", "make_job_args.py"),
               "--samples-file", os.path.join("condor", spec["samples_file"]),
               "--dataset", spec["dataset"], "--location-cfg", spec["location_cfg"],
               "--files-per-job", str(spec["files_per_job"]),
               "--outdir", os.path.relpath(self.filelists_dir, REPO),
               "--job-args", os.path.join("condor", spec["job_args"])]
        if spec["files_per_job_map"]:
            cmd += ["--files-per-job-map", os.path.join("condor", spec["files_per_job_map"])]
        if self.max_files_per_job:
            cmd += ["--max-files-per-job", str(self.max_files_per_job)]
        if spec["replace_xcache"]:
            cmd.append("--replace-xcache")
        if spec["census_skip"]:
            cmd += ["--census-skip", spec["census_skip"]]
        if spec["max_files"]:
            cmd += ["--max-files", str(spec["max_files"])]
        return cmd


def load_campaign(config):
    """The Campaign of a config path, or of a name in condor/campaign_configs/"""
    path = config
    if not os.path.exists(path):
        path = os.path.join(CONFIG_DIR, config if config.endswith(".yaml") else f"{config}.yaml")
    if not os.path.exists(path):
        sys.exit(f"no campaign config {config!r} (looked for it as a path and in {CONFIG_DIR})")
    with open(path, encoding="utf8") as f:
        return Campaign(yaml.safe_load(f), path)


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


def reconcile_all(camp, stat_missing=False):
    """ChunkResults of every chunk of every job list"""
    results = []
    for path in camp.job_lists():
        res, eos_ok, note = cl.reconcile(path, camp.logs_dir, camp.chunk_dir, camp.redirector,
                                         condor_dir=CONDOR, stat_missing=stat_missing)
        if not eos_ok:
            sys.exit(f"cannot list {camp.chunk_dir}: {note}")
        results += res
    return results


def job_specs(camp):
    """{chunk key: ChunkSpec} of every job list"""
    specs = {}
    for path in camp.job_lists():
        specs.update({s.key: s for s in cl.parse_job_args(path)})
    return specs


def per_sample(results):
    """{sample: counts of expected, done and in-progress chunks}"""
    samples = {}
    for r in results:
        s = samples.setdefault(r.sample, {"expected": 0, "done": 0, "in_progress": 0})
        s["expected"] += 1
        s["done"] += r.state == cl.DONE
        s["in_progress"] += r.state == cl.IN_PROGRESS
    return samples


def merged_samples(camp):
    """Samples that already have a merged output in merged_dir"""
    out = subprocess.run(["xrdfs", camp.redirector, "ls", camp.merged_dir], capture_output=True,
                         text=True, check=False)   # the directory does not exist before the first merge
    return {os.path.basename(p)[:-len(".coffea")] for p in out.stdout.split() if p.endswith(".coffea")}


def write_job_list(rows, name):
    """Write (sample, chunk, filelist) rows as a job list in condor/"""
    path = os.path.join(CONDOR, name)
    with open(path, "w", encoding="utf8") as f:
        for sample, chunk, filelist in rows:
            f.write(f"{sample} {chunk} {os.path.relpath(filelist, CONDOR)}\n")
    return path


def submit(camp, job_list, memory=None):
    """condor_submit one job list with condor/submit_campaign.sub"""
    cmd = camp.submit_command(job_list, memory)
    print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=CONDOR, check=True)


def next_round(camp):
    """Number of the next resubmission round: one more than any retry or extra round so far"""
    pattern = re.compile(re.escape(camp.prefix) + r"_(?:retry|extra_round)(\d+)(?:_oom)?\.txt$")
    rounds = [int(m.group(1)) for p in os.listdir(CONDOR) for m in [pattern.match(p)] if m]
    return max(rounds, default=0) + 1


def read_file_list(path):
    """/store paths of the files listed in a bad-file list, if it exists"""
    if not os.path.exists(path):
        return set()
    with open(path, encoding="utf8") as f:
        return {cl.strip_redirector(line.split()[0]) for line in f if line.strip()}


FAILED_FILE = re.compile(r"Failed processing file: WorkItem\(dataset='[^']*', filename='([^']+)'")


def failing_file(result):
    """The input file a chunk failed on, if it failed on that one file in every execution

    coffea names the file of the work item that raised. When Condor's own retries
    (max_retries) failed the same way, at least FAILURES_TO_DROP times in all, the file
    is unreadable rather than unlucky: an OSError is skipped by skipbadfiles, so what
    remains is e.g. a corrupt basket (uproot DeserializationError).
    """
    log, err = result.log_paths.get("log"), result.log_paths.get("err")
    if not log or not err:
        return None
    log, err = (p if os.path.isabs(p) else os.path.join(REPO, p) for p in (log, err))
    if not os.path.exists(log) or not os.path.exists(err):
        return None
    with open(log, encoding="utf8", errors="replace") as f:
        n_failed = len(re.findall(r"Normal termination \(return value [1-9]", f.read()))
    with open(err, encoding="utf8", errors="replace") as f:
        files = set(FAILED_FILE.findall(f.read()))
    return files.pop() if n_failed >= FAILURES_TO_DROP and len(files) == 1 else None


def drop_files(filelist, bad):
    """Remove the bad files from a filelist in place; False if nothing is left"""
    files = [f for f in cl.read_filelist(filelist) if cl.strip_redirector(f) not in bad]
    with open(filelist, "w", encoding="utf8") as f:
        f.write("\n".join(files) + "\n")
    return bool(files)


def cmd_jobs(camp, args):
    """Make the job lists of the campaign with condor/make_job_args.py"""
    for spec in camp.job_list_specs:
        path = os.path.join(CONDOR, spec["job_args"])
        if os.path.exists(path) and not args.force:
            print(f"{os.path.relpath(path, REPO)} exists; keeping it (--force remakes it)")
            continue
        cmd = camp.make_job_args_command(spec)
        print("+", " ".join(cmd))
        subprocess.run(cmd, cwd=REPO, check=True)


def cmd_submit(camp, args):
    """Submit the job lists of the campaign, or the one given"""
    lists = [os.path.join(CONDOR, args.job_list)] if args.job_list else camp.job_lists()
    missing = [p for p in lists if not os.path.exists(p)]
    if missing:
        sys.exit(f"no job list {[os.path.relpath(p, REPO) for p in missing]}: run `jobs` first")
    if not os.path.exists(os.path.join(CONDOR, "sidm_code.tar.gz")):
        sys.exit("condor/sidm_code.tar.gz is missing: package the code first (condor/README.md, section 8)")
    if not os.environ.get("X509_USER_PROXY"):
        print("warning: X509_USER_PROXY is not set; the jobs need a proxy to read from EOS "
              "(condor/README.md, section 11)")
    os.makedirs(camp.logs_dir, exist_ok=True)
    subprocess.run(["xrdfs", camp.redirector, "mkdir", "-p", camp.chunk_dir], check=True)
    for path in lists:
        submit(camp, path)


def cmd_status(camp, args):
    """Totals, completeness and every chunk that is neither done nor in progress"""
    results = reconcile_all(camp, stat_missing=args.stat_missing)
    by_state = {}
    for r in results:
        by_state[r.state] = by_state.get(r.state, 0) + 1
    print(f"{len(results)} chunks: " + ", ".join(f"{k} {v}" for k, v in sorted(by_state.items())))
    samples = per_sample(results)
    complete = sorted(s for s, v in samples.items() if v["done"] == v["expected"])
    merged = merged_samples(camp)
    print(f"{len(complete)}/{len(samples)} samples complete, {len(merged)} merged, "
          f"{len(set(complete) - merged)} complete but not merged")
    for r in results:
        if r.state not in (cl.DONE, cl.IN_PROGRESS):
            print(f"  {r.key:45s} {r.state:18s} {r.reason[:110]}")


def cmd_retry(camp, args):
    """Job list of the failed chunks, with the inputs confirmed missing dropped"""
    results = reconcile_all(camp, stat_missing=True)
    failed = [r for r in results if r.state not in (cl.DONE, cl.IN_PROGRESS)]
    specs = job_specs(camp)
    n = next_round(camp)
    retry_dir = os.path.join(camp.filelists_dir, f"retry{n}")
    bad_path = os.path.join(camp.state_dir, "bad_files.txt")
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
            os.makedirs(camp.state_dir, exist_ok=True)
            with open(bad_path, "a", encoding="utf8") as f:
                f.write(f"{unreadable}  # {r.key}: failed every attempt\n")
            print(f"  {r.key}: dropping {os.path.basename(unreadable)}, which failed every attempt")
        bad_paths = {cl.strip_redirector(b) for b in bad}
        files = [f for f in files if cl.strip_redirector(f) not in bad_paths]
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
    known_bad = read_file_list(bad_path)
    rows = [(s, c, f) for s, c, f in rows if drop_files(f, known_bad)]
    rows_oom = [(s, c, f) for s, c, f in rows_oom if drop_files(f, known_bad)]
    if not rows and not rows_oom:
        print("nothing to retry")
        return
    for chunk_rows, suffix, memory in ((rows, "", None), (rows_oom, "_oom", camp.oom_memory)):
        if chunk_rows:
            path = write_job_list(chunk_rows, f"{camp.prefix}_retry{n}{suffix}.txt")
            print(f"wrote {path} ({len(chunk_rows)} chunks)")
            if args.submit:
                submit(camp, path, memory)


def local_copy(camp, eos_path):
    """A chunk output, through the EOS FUSE mount when there is one, else copied locally"""
    fuse = "/eos/uscms" + eos_path
    if os.path.exists(fuse):
        return fuse
    local = os.path.join(camp.state_dir, "tmp", os.path.basename(eos_path))
    os.makedirs(os.path.dirname(local), exist_ok=True)
    subprocess.run(["xrdcp", "-s", "-f", f"{camp.redirector}/{eos_path}", local], check=True)
    return local


def processed_files(camp, eos_path):
    """Input files a chunk output says it processed"""
    # only this command needs coffea
    from coffea.util import load   # pylint: disable=import-outside-toplevel
    out = load(local_copy(camp, eos_path))
    return {item.filename for item in out.get("processed", ())}


def n_entries(url):
    """Events in a file, or None if it cannot be opened"""
    import uproot   # pylint: disable=import-outside-toplevel
    fuse = "/eos/uscms" + cl.strip_redirector(url)
    try:
        with uproot.open(fuse if os.path.exists(fuse) else url) as f:
            return f["Events"].num_entries
    except Exception:   # pylint: disable=broad-except
        return None


def cmd_skipped(camp, args):
    """Extra chunks of the files that finished chunks skipped"""
    results = reconcile_all(camp)
    specs = job_specs(camp)
    checked_path = os.path.join(camp.state_dir, "skipped_checked.json")
    checked = {}
    if os.path.exists(checked_path):
        with open(checked_path, encoding="utf8") as f:
            checked = json.load(f)
    os.makedirs(camp.state_dir, exist_ok=True)
    for r in results:
        if r.state != cl.DONE or r.key in checked:
            continue
        files = cl.read_filelist(os.path.join(CONDOR, specs[r.key].filelist))
        # compare /store paths: the filelists may name a different redirector
        done = {cl.strip_redirector(f) for f in processed_files(camp, f"{camp.chunk_dir}/{r.key}.coffea")}
        checked[r.key] = [f for f in files if cl.strip_redirector(f) not in done]
        with open(checked_path, "w", encoding="utf8") as f:
            json.dump(checked, f, indent=1)
    # files dropped by `retry` as unreadable are not skipped, they are gone
    dropped = read_file_list(os.path.join(camp.state_dir, "bad_files.txt"))
    skipped = {k: [f for f in v if cl.strip_redirector(f) not in dropped] for k, v in checked.items()}
    skipped = {k: v for k, v in skipped.items() if v}
    n_skipped = sum(len(v) for v in skipped.values())
    print(f"{len(checked)} finished chunks checked; {n_skipped} input file(s) skipped in {len(skipped)} chunk(s)")

    # a file is retried once in an extra chunk; skipped again, it is bad
    extra_specs = cl.parse_job_args(camp.extra_list) if os.path.exists(camp.extra_list) else []
    retried = set()
    for spec in extra_specs:
        retried.update(cl.read_filelist(os.path.join(CONDOR, spec.filelist)))
    bad = sorted({f for k, v in skipped.items() if int(k.rsplit("_", 1)[1]) >= EXTRA_CHUNK_OFFSET
                  for f in v})
    if bad:
        with open(os.path.join(camp.state_dir, "bad_files_skipped.txt"), "w", encoding="utf8") as f:
            f.write("\n".join(bad) + "\n")
        print(f"{len(bad)} file(s) skipped twice: listed as bad in {camp.state_dir}/bad_files_skipped.txt")
    todo = {}
    for key, files in skipped.items():
        sample = key.rsplit("_", 1)[0]
        todo.setdefault(sample, []).extend(f for f in files if f not in retried)
    # a file without events has no work item, so it never shows up as processed
    empty = [f for files in todo.values() for f in files if n_entries(f) == 0]
    if empty:
        with open(os.path.join(camp.state_dir, "empty_files.txt"), "w", encoding="utf8") as f:
            f.write("\n".join(sorted(empty)) + "\n")
        print(f"{len(empty)} of them have no events (listed in {camp.state_dir}/empty_files.txt)")
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
    os.makedirs(os.path.join(camp.filelists_dir, "extra"), exist_ok=True)
    for sample, files in sorted(todo.items()):
        chunk = next_chunk.get(sample, EXTRA_CHUNK_OFFSET)
        filelist = os.path.join(camp.filelists_dir, "extra", f"{sample}_{chunk}.txt")
        with open(filelist, "w", encoding="utf8") as f:
            f.write("\n".join(files) + "\n")
        rows.append((sample, chunk, filelist))
    with open(camp.extra_list, "a", encoding="utf8") as f:
        for sample, chunk, filelist in rows:
            f.write(f"{sample} {chunk} {os.path.relpath(filelist, CONDOR)}\n")
    new_list = write_job_list(rows, f"{camp.prefix}_extra_round{next_round(camp)}.txt")
    print(f"{len(rows)} extra chunk(s) appended to {camp.extra_list}; this round: {new_list}")
    submit(camp, new_list)


def cmd_merge(camp, args):
    """Merge the complete samples into merged_dir"""
    results = reconcile_all(camp)
    samples = per_sample(results)
    if args.samples:
        wanted = [s.strip() for s in args.samples.split(",") if s.strip()]
    else:
        merged = set() if args.remerge else merged_samples(camp)
        wanted = sorted(s for s, v in samples.items() if v["done"] == v["expected"] and s not in merged)
    unknown = [s for s in wanted if s not in samples]
    if unknown:
        sys.exit(f"{unknown} are not in this campaign's job lists")
    incomplete = [s for s in wanted if samples[s]["done"] != samples[s]["expected"]]
    if incomplete:
        sys.exit(f"not every chunk of {incomplete} is done; refusing to merge them")
    if not wanted:
        print("nothing to merge")
        return
    os.makedirs(camp.state_dir, exist_ok=True)
    sample_file = os.path.join(camp.state_dir, "merge_samples.txt")
    with open(sample_file, "w", encoding="utf8") as f:
        f.write("\n".join(wanted) + "\n")
    # the sidecar lists the input files of every chunk of a sample, extra chunks
    # included, without the files dropped as unreadable
    bad = set()
    for name in ("bad_files.txt", "bad_files_skipped.txt"):
        bad |= read_file_list(os.path.join(camp.state_dir, name))
    filelists = os.path.join(camp.state_dir, "filelists_all")
    os.makedirs(filelists, exist_ok=True)
    for path in glob.glob(os.path.join(camp.filelists_dir, "**", "*.txt"), recursive=True):
        if "/retry" in path:
            continue    # a retry replaces its chunk's filelist, it does not add one
        files = [f for f in cl.read_filelist(path) if cl.strip_redirector(f) not in bad]
        with open(os.path.join(filelists, os.path.basename(path)), "w", encoding="utf8") as f:
            f.write("\n".join(files) + "\n")
    cmd = [sys.executable, os.path.join(REPO, "sidm", "scripts", "merge_coffea_chunks_eos.py"),
           "--input-eos-dir", camp.chunk_dir, "--output-eos-dir", camp.merged_dir,
           "--workdir", os.path.join(camp.state_dir, "merge_tmp"), "--samples", sample_file,
           "--filelists-dir", filelists, "--selections", ",".join(camp.channels),
           "--hist-collections", ",".join(camp.hist_collections)]
    print(f"merging {len(wanted)} sample(s)")
    subprocess.run(cmd, cwd=REPO, check=True)


def build_parser():
    """The command-line interface"""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", maxsplit=1)[0])
    parser.add_argument("config", help="campaign config: a path, or a name in condor/campaign_configs/")
    sub = parser.add_subparsers(dest="cmd", required=True)
    jb = sub.add_parser("jobs", help="make the job lists")
    jb.add_argument("--force", action="store_true", help="remake job lists that exist")
    jb.set_defaults(func=cmd_jobs)
    sm = sub.add_parser("submit", help="submit the job lists")
    sm.add_argument("--job-list", default=None, help="submit only this job list (relative to condor/)")
    sm.set_defaults(func=cmd_submit)
    st = sub.add_parser("status", help="reconcile the job lists against the logs and EOS")
    st.add_argument("--stat-missing", action="store_true",
                    help="stat the inputs of failed chunks to tell missing files from timeouts")
    st.set_defaults(func=cmd_status)
    rt = sub.add_parser("retry", help="job list of the failed chunks")
    rt.add_argument("--submit", action="store_true")
    rt.set_defaults(func=cmd_retry)
    sk = sub.add_parser("skipped", help="extra chunks of the input files that were skipped")
    sk.add_argument("--submit", action="store_true")
    sk.set_defaults(func=cmd_skipped)
    mg = sub.add_parser("merge", help="merge the complete samples")
    mg.add_argument("--samples", default=None, help="comma-separated samples to merge")
    mg.add_argument("--remerge", action="store_true", help="also merge samples merged before")
    mg.set_defaults(func=cmd_merge)
    return parser


def main(argv=None):
    """Command-line entry point"""
    args = build_parser().parse_args(argv)
    args.func(load_campaign(args.config), args)


if __name__ == "__main__":
    main()
