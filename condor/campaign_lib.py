"""Helpers for the SIDM Condor campaign orchestrator.

Phase 0 of the orchestrator is a *read-only* reconciler: given the job arguments
that were submitted (``job_args.txt``), the per-job Condor logs, and the chunk
outputs that actually landed on EOS, it works out the state of every chunk,
classifies each failure by its root cause, and renders one human report plus a
machine-readable JSON record. It submits and resubmits nothing.

To keep that reconcile fast and robust even when the analysis environment is not
available, this module uses only the Python standard library. It shells out to
``xrdfs`` (EOS) but imports no coffea / awkward / sidm. Do not add heavy imports here.

The output filename contract is ``{sample}_{chunk}.coffea`` and the Condor log
filename contract is ``{sample}_{chunk}_{cluster}_{proc}.{log,out,err}`` -- both
fixed by ``condor/submit.sub``. The ``{sample}_{chunk}`` stem is the chunk key
used everywhere below; it matches the stem parsed by
``sidm/scripts/merge_coffea_chunks_eos.py``. Keep the two parsers consistent.
"""

import os
import re
import subprocess
from collections import OrderedDict
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone


# --------------------------------------------------------------------------
# Failure taxonomy. Each chunk that is not DONE gets exactly one of these.
# The strings are the public vocabulary of the report -- keep them stable.
# --------------------------------------------------------------------------
DONE = "DONE"                          # output present on EOS, size > 0
ZERO_BYTE_OUTPUT = "ZERO_BYTE_OUTPUT"  # output on EOS but 0 bytes (partial copy)
MISSING_ROOT_FILE = "MISSING_ROOT_FILE"  # an input ROOT file is confirmed absent now
EOS_TIMEOUT = "EOS_TIMEOUT"            # transient EOS read failure; inputs still present/unprobed
OOM = "OOM"                            # killed by signal / memory exhausted
HELD_PROXY = "HELD_PROXY"              # held because the x509 proxy was unreadable
HELD_TRANSFER = "HELD_TRANSFER"        # held on a non-proxy file-transfer failure
STALLED = "STALLED"                    # removed by our periodic_remove stall guard
REMOVED = "REMOVED"                    # removed manually (condor_rm) -- needs a human
CODE_BUG = "CODE_BUG"                  # deterministic python exception; a retry will re-fail
XRDCP_FAILED = "XRDCP_FAILED"          # processor ran but the output copy to EOS failed
MISSING_OUTPUT = "MISSING_OUTPUT"      # job returned 0 but no output on EOS (look at the dir)
UNKNOWN = "UNKNOWN"                    # terminal, no output, nothing above matched
LOST = "LOST"                          # no terminal log event found for this chunk at all
IN_PROGRESS = "IN_PROGRESS"            # submitted/running, no terminal event yet (not a failure)
EOS_UNLISTABLE = "EOS_UNLISTABLE"      # EOS dir could not be listed; output presence unknown

# Classes a capped resubmit (a later phase) should NOT auto-retry, because the
# next attempt would deterministically reproduce the same failure.
DO_NOT_AUTO_RETRY = {CODE_BUG, REMOVED}


# --------------------------------------------------------------------------
# Small dataclasses
# --------------------------------------------------------------------------
@dataclass
class ChunkSpec:
    """One line of job_args.txt: the intent of a single chunk."""
    sample: str
    chunk: int
    filelist: str           # path as written in job_args.txt, relative to condor/

    @property
    def key(self):
        return f"{self.sample}_{self.chunk}"


@dataclass
class AttemptLog:
    """Summary of one Condor .log (one (cluster, proc) attempt of one chunk)."""
    sample: str
    chunk: int
    cluster: int
    proc: int
    log_path: str
    err_path: str = ""
    out_path: str = ""
    submitted: bool = False
    executed: bool = False
    terminated: bool = False
    return_value: int = None      # exit code if "Normal termination"
    signal: int = None            # signal number if "Abnormal termination"
    held: bool = False
    hold_reason: str = ""
    aborted: bool = False
    abort_reason: str = ""
    evicted: bool = False
    released: bool = False        # a 013 "Job was released" event followed a hold
    error_event: str = ""         # text of a 021 "Error from starter" event, if any
    peak_mem_mb: int = 0
    request_mem_mb: int = 0
    last_event_local: str = ""    # timestamp of the last event seen (schedd-LOCAL time, not UTC)


@dataclass
class ChunkResult:
    """The reconciled verdict for one chunk."""
    sample: str
    chunk: int
    key: str
    state: str                    # DONE | one of the failure classes
    reason: str = ""
    offending_root_files: list = field(default_factory=list)
    n_files: int = 0              # number of input ROOT files in the chunk
    n_attempts: int = 0           # distinct (cluster, proc) attempts seen in the logs
    eos_size: int = None          # bytes of the EOS output, or None if absent
    auto_retryable: bool = True
    log_paths: dict = field(default_factory=dict)  # newest attempt's {log,err,out}


# --------------------------------------------------------------------------
# Parsing job_args.txt and filelists
# --------------------------------------------------------------------------
def parse_job_args(path):
    """Read job_args.txt -> list[ChunkSpec]. Lines are 'sample chunk filelist'."""
    specs = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            sample, chunk, filelist = parts[0], parts[1], parts[2]
            try:
                chunk_i = int(chunk)
            except ValueError:
                continue
            specs.append(ChunkSpec(sample=sample, chunk=chunk_i, filelist=filelist))
    return specs


def read_filelist(path):
    """Return the ROOT-file URLs listed in a chunk filelist (skipping comments)."""
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [ln.strip() for ln in f
                if ln.strip() and not ln.strip().startswith("#")]


def parse_log_name(filename):
    """'{sample}_{chunk}_{cluster}_{proc}.log' -> (sample, chunk, cluster, proc).

    sample may itself contain underscores, so we peel the last three integer
    fields off the right. Returns None if the tail is not three integers.
    """
    stem = re.sub(r"\.(log|out|err)$", "", os.path.basename(filename))
    bits = stem.rsplit("_", 3)
    if len(bits) != 4:
        return None
    sample, chunk, cluster, proc = bits
    try:
        return sample, int(chunk), int(cluster), int(proc)
    except ValueError:
        return None


def parse_coffea_stem(basename):
    """'{sample}_{chunk}.coffea' -> (sample, chunk). Mirrors the merge script."""
    stem = re.sub(r"\.coffea$", "", os.path.basename(basename))
    sample, _, chunk = stem.rpartition("_")
    if not sample or not chunk.isdigit():
        return None
    return sample, int(chunk)


# --------------------------------------------------------------------------
# Parsing a Condor user log (.log) -- the durable, local, per-job event record.
# Event codes seen on LPC: 000 submitted, 001 executing, 004 evicted,
# 005 terminated, 006 image-size update, 009 aborted, 012 held, 013 released,
# 021 error.
# Each event is a header line "NNN (cl.pr.sub) DATE TIME text" followed by
# indented body lines and closed by a line that is exactly "...".
# --------------------------------------------------------------------------
_EVENT_HEADER = re.compile(
    r"^(?P<code>\d{3}) \((?P<cl>\d+)\.(?P<pr>\d+)\.\d+\) "
    r"(?P<date>\d{4}-\d{2}-\d{2}) (?P<time>\d{2}:\d{2}:\d{2}) (?P<text>.*)$"
)
_RETURN_VALUE = re.compile(r"Normal termination \(return value (\d+)\)")
_SIGNAL = re.compile(r"Abnormal termination \(signal (\d+)\)")
_MEM_USAGE = re.compile(r"(\d+)\s+-\s+MemoryUsage of job \(MB\)")
_REQUEST_MEM = re.compile(r"^\s*Memory = (\d+)\s*$", re.M)


def _iter_events(lines):
    """Yield (header_match, body_lines) for each ULOG event.

    Collect a body until the "..." terminator OR the next event header. Resyncing
    on the next header means a missing/truncated terminator on one event cannot
    swallow every later event (and the real 005/009 terminal) into its body.
    """
    i, n = 0, len(lines)
    while i < n:
        m = _EVENT_HEADER.match(lines[i])
        if not m:
            i += 1
            continue
        body = []
        j = i + 1
        while j < n and lines[j].rstrip() != "..." and not _EVENT_HEADER.match(lines[j]):
            body.append(lines[j])
            j += 1
        yield m, body
        if j < n and lines[j].rstrip() == "...":
            i = j + 1      # consume the terminator
        else:
            i = j          # resync on the next header (or EOF)


def parse_event_log(path):
    """Parse one .log file into an AttemptLog (without err/out paths filled in)."""
    info = parse_log_name(path)
    if info is None:
        return None
    sample, chunk, cluster, proc = info
    a = AttemptLog(sample=sample, chunk=chunk, cluster=cluster, proc=proc,
                   log_path=path)
    try:
        with open(path, errors="replace") as f:
            lines = f.read().splitlines()
    except OSError:
        return a

    for m, body in _iter_events(lines):
        code = m.group("code")
        a.last_event_local = f"{m.group('date')} {m.group('time')}"
        text = m.group("text")
        body_text = "\n".join(body)
        if code == "000":
            a.submitted = True
        elif code == "001":
            a.executed = True
            rm = _REQUEST_MEM.search(body_text)
            if rm:
                a.request_mem_mb = int(rm.group(1))
        elif code == "005":
            a.terminated = True
            rv = _RETURN_VALUE.search(body_text)
            sg = _SIGNAL.search(body_text)
            if rv:
                a.return_value = int(rv.group(1))
            elif sg:
                a.signal = int(sg.group(1))
        elif code == "006":
            for mm in _MEM_USAGE.finditer(body_text):
                a.peak_mem_mb = max(a.peak_mem_mb, int(mm.group(1)))
        elif code == "009":
            a.aborted = True
            a.abort_reason = _first_nonempty(body) or text
        elif code == "012":
            a.held = True
            a.hold_reason = _first_nonempty(body) or text
        elif code == "021":
            a.error_event = (text + " " + _first_nonempty(body)).strip()
        elif code == "004":
            a.evicted = True
        elif code == "013":
            a.released = True
    return a


def _first_nonempty(body_lines):
    for ln in body_lines:
        s = ln.strip()
        if s:
            return s
    return ""


# --------------------------------------------------------------------------
# Reading the .err tail for python-level error signatures.
# --------------------------------------------------------------------------
# Ordered most-specific first. The first signature that matches the .err tail
# picks the class for a job that ran but exited non-zero.
_PROXY_HOLD = re.compile(r"x509up|/tmp/x509up|Transfer input files failure", re.I)
_OPEN_FAIL = re.compile(r"Failed to open file|No such file or directory|"
                        r"\[ERROR\].*(?:Operation|open)", re.I)
_EOS_TIMEOUT_SIG = re.compile(r"Operation expired|Socket timeout|"
                              r"Empty list provided to reduction", re.I)
_OOM_SIG = re.compile(r"MemoryError|Killed|Out of memory|std::bad_alloc", re.I)
_XRDCP_FAIL = re.compile(r"xrdcp failed", re.I)
# Deterministic python exceptions: a retry re-fails. The final traceback line.
_CODE_BUG_EXC = re.compile(
    r"^\s*(KeyError|AttributeError|ImportError|ModuleNotFoundError|TypeError|"
    r"IndexError|NameError|SyntaxError|FileNotFoundError):", re.M)
_ROOT_URL = re.compile(r"root://\S+?\.root")


def read_err_tail(path, max_bytes=20000):
    """Return the last ~max_bytes of a .err/.out file as text (or '')."""
    if not path or not os.path.isfile(path):
        return ""
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as f:
            if size > max_bytes:
                f.seek(-max_bytes, os.SEEK_END)
            data = f.read()
        return data.decode("utf-8", errors="replace")
    except OSError:
        return ""


def last_traceback_line(err_text):
    """The exception line of the last python traceback in err_text, if any."""
    matches = list(_CODE_BUG_EXC.finditer(err_text))
    if matches:
        return matches[-1].group(0).strip()
    # Fall back to the last non-empty line if a bare 'Traceback' is present.
    if "Traceback" in err_text:
        for ln in reversed(err_text.splitlines()):
            if ln.strip():
                return ln.strip()
    return ""


# --------------------------------------------------------------------------
# EOS reconcile via xrdfs (read-only; works with just a Kerberos ticket).
# --------------------------------------------------------------------------
def _run(cmd, timeout=120):
    """Run a subprocess; return (returncode, stdout, stderr). -1 on timeout."""
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "timeout"
    except FileNotFoundError as e:
        return -2, "", str(e)


_DATE_TOK = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def xrdfs_ls_sizes(redirector, eos_dir, timeout=180):
    """List an EOS directory once. Return (ok, {basename: size_bytes}, note).

    ok is False if the directory could not be listed; ``note`` then says why:
    "not-found" for a positively-absent dir (usually a typo'd --eos-chunk-dir),
    vs "error: ..." for a transient/degraded EOS or a missing xrdfs binary. Uses
    'xrdfs <redir> ls -l <dir>' so we get sizes in the same call.
    """
    rc, out, err = _run(["xrdfs", redirector, "ls", "-l", eos_dir], timeout=timeout)
    if rc == 0:
        return True, parse_ls_l_sizes(out), ""
    blob = (out + " " + err).lower()
    if "no such file" in blob or "[3011]" in blob:
        return False, {}, "not-found"        # dir positively absent -- usually a typo
    return False, {}, "error: " + (err.strip() or out.strip() or f"rc={rc}")[:200]


def parse_ls_l_sizes(out_text):
    """Parse 'xrdfs ls -l' stdout into {basename: size_bytes}. Pure / testable.

    FNAL EOS format: "<flags> <owner> <group> <size> <YYYY-MM-DD> <HH:MM:SS> <path>".
    The column count varies by xrootd build, so anchor on the date token and take
    the size as the integer immediately before it; the path is the last token.
    """
    sizes = {}
    for line in out_text.splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        name = os.path.basename(parts[-1])
        # Size is the integer immediately before the YYYY-MM-DD date column. We
        # deliberately do NOT guess a size when the date anchor is absent: a wrong
        # non-zero size would be a false "DONE" (size>0). No size -> not recorded
        # -> the chunk is treated as not-yet-on-EOS, which is the safe direction.
        for i, tok in enumerate(parts):
            if i >= 1 and _DATE_TOK.match(tok) and parts[i - 1].isdigit():
                sizes[name] = int(parts[i - 1])
                break
    return sizes


# Probe outcomes for a single file's existence check.
ABSENT = "absent"            # positively confirmed not present
PRESENT = "present"          # confirmed present
INCONCLUSIVE = "inconclusive"  # could not determine (EOS error / timeout)


def url_to_eos_path(root_url):
    """root://host//store/... -> /store/...  (the path xrdfs stat expects)."""
    m = re.search(r"root://[^/]+/+(/store/\S+)", root_url)
    if m:
        return m.group(1)
    m = re.search(r"(/store/\S+)", root_url)
    return m.group(1) if m else root_url


def strip_redirector(path):
    """root://host//store/... -> /store/... ; a bare path is returned unchanged.

    'xrdfs <redir> ls <dir>' wants <dir> WITHOUT the redirector, so a full root://
    URL pasted out of submit.sub would otherwise fail as a missing directory.
    """
    m = re.match(r"root://[^/]+/+(/.*)$", path)
    return m.group(1) if m else path


def xrdfs_stat_exists(redirector, root_url, timeout=60):
    """Probe one ROOT file. Returns ABSENT / PRESENT / INCONCLUSIVE.

    The asymmetry here is the whole safety story for a later prune step: we only
    ever report ABSENT when the server explicitly says 'No such file'. A timeout
    or any other error is INCONCLUSIVE, never ABSENT -- so an EOS outage can
    never be mistaken for a genuinely missing input.
    """
    path = url_to_eos_path(root_url)
    rc, out, err = _run(["xrdfs", redirector, "stat", path], timeout=timeout)
    if rc == 0:
        return PRESENT
    if rc < 0:
        # xrdfs did not run (not installed) or timed out -> we learned nothing.
        # Never read Python's own "[Errno 2] ... 'xrdfs'" as a server "No such
        # file": ABSENT must rest only on a positive answer from EOS itself.
        return INCONCLUSIVE
    blob = (out + " " + err).lower()
    if "no such file" in blob or "[3011]" in blob:
        return ABSENT
    return INCONCLUSIVE


# --------------------------------------------------------------------------
# Classification: combine the newest attempt log, the .err tail, and EOS
# reality into one verdict per chunk.
# --------------------------------------------------------------------------
def classify_chunk(spec, attempts, eos_size, root_files, redirector,
                   stat_missing=True, eos_listable=True):
    """Return a ChunkResult for one chunk.

    spec        : the ChunkSpec (intent).
    attempts    : list[AttemptLog] for this chunk, any order.
    eos_size    : size in bytes of the EOS output, or None if absent/unknown.
    root_files  : the input ROOT URLs of this chunk (for prune diagnosis).
    redirector  : EOS redirector for the (read-only) per-file existence probe.
    stat_missing: if True, xrdfs-stat the implicated ROOT files to split
                  MISSING_ROOT_FILE from EOS_TIMEOUT. If False, default to
                  EOS_TIMEOUT (never prune on an unprobed open failure).
    eos_listable: whether the EOS output dir could be listed. If False, eos_size
                  is None for every chunk and "no output" cannot be concluded, so
                  an otherwise-clean job is reported EOS_UNLISTABLE, not missing.
    """
    # Order attempts chronologically: ClusterIds are NOT monotonic in time on LPC
    # (jobs span multiple schedds, e.g. 2.8M vs 84M), so the last event timestamp
    # is a better "newest" key than the cluster number. (Schedd-local time, so
    # cross-schedd ordering can still skew under a clock offset -- acceptable; it
    # is strictly better than cluster id and correct within a schedd.)
    attempts = sorted(attempts, key=lambda a: (a.last_event_local, a.cluster, a.proc))
    newest = attempts[-1] if attempts else None
    res = ChunkResult(
        sample=spec.sample, chunk=spec.chunk, key=spec.key,
        n_files=len(root_files), n_attempts=len(attempts), eos_size=eos_size,
        log_paths=({"log": newest.log_path, "err": newest.err_path,
                    "out": newest.out_path} if newest else {}),
        state=UNKNOWN,
    )

    # 1. EOS reality wins. A present, non-empty output is DONE no matter what
    #    any earlier attempt's logs say (a late success or a manual rerun).
    if eos_size is not None and eos_size > 0:
        res.state = DONE
        res.reason = "output present on EOS"
        return res
    if eos_size == 0:
        res.state = ZERO_BYTE_OUTPUT
        res.reason = "output on EOS is 0 bytes (partial copy); rerun"
        res.auto_retryable = True
        return res

    # 2. No usable log at all -> the job never ran or its log is gone.
    if newest is None:
        res.state = LOST
        res.reason = "no Condor log found for this chunk"
        return res

    err_text = read_err_tail(newest.err_path) + "\n" + read_err_tail(newest.out_path)
    any_proxy_hold = any(_PROXY_HOLD.search(a.hold_reason or "") for a in attempts)

    # 3. The newest attempt TERMINATED -> its exit is the definitive outcome of
    #    this chunk, regardless of holds/evictions on earlier attempts (or earlier
    #    in this same attempt). A terminal exit settles the verdict.
    if newest.terminated:
        rv = newest.return_value
        # 3a. Exited clean. If EOS could not be listed we cannot conclude "missing".
        if rv == 0:
            if not eos_listable:
                res.state = EOS_UNLISTABLE
                res.reason = "exit 0 but the EOS dir could not be listed -- output presence unknown"
                return res
            res.state = MISSING_OUTPUT
            res.reason = ("job returned 0 but no .coffea on EOS -- wrong output dir, "
                          "or the file was deleted")
            return res
        # 3b. Killed by OOM: signal 9, exit 137 (128+SIGKILL from the OOM killer),
        #     or an explicit OOM message.
        if newest.signal == 9 or newest.return_value == 137 or _OOM_SIG.search(err_text):
            res.state = OOM
            res.reason = (f"out of memory (signal 9 / exit 137 / OOM message; peak "
                          f"{newest.peak_mem_mb} MB of {newest.request_mem_mb} MB requested); "
                          f"resubmit at higher memory")
            return res
        # 3c. Ran and exited non-zero -> read the python error to subclassify.
        if rv not in (None, 0):
            # Deterministic python exception -> a retry would re-fail.
            if _CODE_BUG_EXC.search(err_text):
                res.state = CODE_BUG
                res.reason = last_traceback_line(err_text) or "deterministic python exception"
                res.auto_retryable = False
                return res
            # Input-open failure -> confirm (read-only) whether a file is really gone.
            if _OPEN_FAIL.search(err_text):
                # Only consider URLs that are actually THIS chunk's declared inputs;
                # a redirector/pileup echo of some other absent file must never get
                # a real input pruned later. If the filelist could not be read
                # (root_files empty), nothing is declared -> nothing is ever pruned.
                declared = {url_to_eos_path(u) for u in root_files}
                implicated = [u for u in dict.fromkeys(_ROOT_URL.findall(err_text))
                              if url_to_eos_path(u) in declared]
                confirmed_absent = []
                if stat_missing:
                    for url in implicated:
                        if xrdfs_stat_exists(redirector, url) == ABSENT:
                            confirmed_absent.append(url)
                if confirmed_absent:
                    res.state = MISSING_ROOT_FILE
                    res.offending_root_files = confirmed_absent
                    res.reason = (f"{len(confirmed_absent)} input ROOT file(s) confirmed "
                                  f"absent on EOS")
                    return res
                # Open failure but no declared input confirmed gone -> transient (never prune).
                res.state = EOS_TIMEOUT
                res.reason = ("input-open failure but no declared input confirmed missing "
                              "(probe inconclusive, files present, or non-input URL) -> transient")
                return res
            # Explicit transient EOS markers.
            if _EOS_TIMEOUT_SIG.search(err_text):
                res.state = EOS_TIMEOUT
                res.reason = "transient EOS read failure (Operation expired / empty reduction)"
                return res
            # The processor ran but the output copy to EOS failed.
            if _XRDCP_FAIL.search(err_text):
                res.state = XRDCP_FAILED
                res.reason = "processor finished but xrdcp of the output to EOS failed"
                return res
            # No specific signature. Note memory pressure if RSS was near the
            # request, but do NOT assert OOM without corroboration -- a job that
            # merely drifted to ~95% and then died for an unrelated reason should
            # not be (mis)routed to "bump the memory".
            mem_note = ""
            if (newest.request_mem_mb and newest.peak_mem_mb
                    and newest.peak_mem_mb >= 0.95 * newest.request_mem_mb):
                mem_note = (f" [peak {newest.peak_mem_mb} MB of {newest.request_mem_mb} MB "
                            f"requested -- possibly memory]")
            res.state = UNKNOWN
            res.reason = "non-zero exit; see .err: " + _tail_snippet(err_text) + mem_note
            return res
        # 3d. Abnormal termination by a signal other than 9.
        res.state = UNKNOWN
        res.reason = (f"abnormal termination (signal {newest.signal}); see .err: "
                      + _tail_snippet(err_text))
        return res

    # 4. The newest attempt was REMOVED (aborted) without terminating.
    if newest.aborted:
        if "SIDM_STALL" in newest.abort_reason:
            res.state = STALLED
            res.reason = "removed by stall guard: " + newest.abort_reason
            return res
        if any_proxy_hold:
            res.state = HELD_PROXY
            res.reason = "removed after an x509 proxy hold; renew proxy and resubmit"
            return res
        if "condor_rm" in newest.abort_reason:
            res.state = REMOVED
            res.reason = "removed manually (" + newest.abort_reason + "); needs a human"
            res.auto_retryable = False
            return res
        # A non-manual removal (e.g. a SYSTEM_PERIODIC_REMOVE policy) is a stall-like,
        # retryable removal rather than a human action.
        res.state = STALLED
        res.reason = "removed by a system/periodic policy: " + (newest.abort_reason or "unknown")
        return res

    # 5. The newest attempt is still HELD (never terminated, aborted, or released).
    if newest.held and not newest.released:
        reason_text = newest.hold_reason
        if _PROXY_HOLD.search(reason_text) and ("x509" in reason_text.lower()
                                                or "/tmp/x509" in reason_text.lower()):
            res.state = HELD_PROXY
            res.reason = ("held: x509 proxy unreadable by the schedd "
                          "(copy proxy to shared NFS; see README §11)")
            return res
        res.state = HELD_TRANSFER
        res.reason = "held: " + (reason_text[:160] or "input file-transfer failure")
        return res

    # 6. Submitted/executed, no terminal event, not currently held -> still in
    #    flight (queued or running). Not a failure; a live condor_q (a later phase)
    #    settles it. This keeps a mid-campaign reconcile honest.
    if newest.submitted or newest.executed:
        res.state = IN_PROGRESS
        res.reason = "no terminal Condor event yet; job appears queued or still running"
        return res

    # 7. A log exists but shows no recognizable lifecycle -> treat as lost.
    res.state = LOST
    res.reason = "no recognizable Condor lifecycle events in the log"
    return res


def _tail_snippet(err_text, n=3):
    lines = [ln for ln in err_text.splitlines() if ln.strip()]
    return " | ".join(lines[-n:])[:300]


# --------------------------------------------------------------------------
# Top-level reconcile: glue parsing + EOS + classification together.
# --------------------------------------------------------------------------
def collect_attempts(logs_dir):
    """Index condor/logs into {(sample, chunk): [AttemptLog, ...]}."""
    by_chunk = {}
    if not os.path.isdir(logs_dir):
        return by_chunk
    for name in os.listdir(logs_dir):
        if not name.endswith(".log"):
            continue
        a = parse_event_log(os.path.join(logs_dir, name))
        if a is None:
            continue
        base = name[:-4]  # strip .log
        for ext, attr in ((".err", "err_path"), (".out", "out_path")):
            cand = os.path.join(logs_dir, base + ext)
            if os.path.isfile(cand):
                setattr(a, attr, cand)
        by_chunk.setdefault((a.sample, a.chunk), []).append(a)
    return by_chunk


def reconcile(job_args_path, logs_dir, eos_chunk_dir, redirector,
              condor_dir=None, stat_missing=True):
    """Reconcile a campaign. Returns (results, eos_ok, eos_note).

    results : list[ChunkResult], one per chunk in job_args.txt.
    eos_ok  : False if the EOS chunk dir could not be listed (report says so).
    eos_note: "" on success, else "not-found" / "error: ..." explaining why.

    condor_dir: directory the filelist paths in job_args.txt are relative to
                (normally the 'condor/' dir). Defaults to the dir of job_args.
    """
    specs = parse_job_args(job_args_path)
    attempts = collect_attempts(logs_dir)
    eos_ok, eos_sizes, eos_note = xrdfs_ls_sizes(redirector, strip_redirector(eos_chunk_dir))
    if condor_dir is None:
        condor_dir = os.path.dirname(os.path.abspath(job_args_path))

    results = []
    for spec in specs:
        root_files = read_filelist(os.path.join(condor_dir, spec.filelist))
        eos_size = eos_sizes.get(spec.key + ".coffea") if eos_ok else None
        res = classify_chunk(spec, attempts.get((spec.sample, spec.chunk), []),
                             eos_size, root_files, redirector,
                             stat_missing=stat_missing, eos_listable=eos_ok)
        results.append(res)
    return results, eos_ok, eos_note


# --------------------------------------------------------------------------
# Report building + rendering.
# --------------------------------------------------------------------------
def _utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_report(results, eos_ok, inputs, run_id="manual", eos_note=""):
    """Project the per-chunk results into the report.json structure."""
    by_class = OrderedDict()
    for r in results:
        by_class[r.state] = by_class.get(r.state, 0) + 1

    done = by_class.get(DONE, 0)
    in_progress = by_class.get(IN_PROGRESS, 0)
    unlistable = by_class.get(EOS_UNLISTABLE, 0)
    expected = len(results)

    # "failed" excludes DONE, in-flight jobs, and chunks whose presence is unknown
    # because EOS could not be listed -- none of those are confirmed failures.
    not_failure = {DONE, IN_PROGRESS, EOS_UNLISTABLE}

    per_sample = OrderedDict()
    for r in results:
        s = per_sample.setdefault(r.sample, {"sample": r.sample, "expected": 0,
                                             "done": 0, "missing_chunks": []})
        s["expected"] += 1
        if r.state == DONE:
            s["done"] += 1
        else:
            s["missing_chunks"].append(r.chunk)
    for s in per_sample.values():
        s["complete"] = (s["done"] == s["expected"])
        s["missing_chunks"].sort()

    failed = [r for r in results if r.state not in not_failure]
    terminal = [r for r in failed if not r.auto_retryable]

    report = OrderedDict()
    report["schema_version"] = 1
    report["run_id"] = run_id
    report["generated_utc"] = _utc_now()
    report["inputs"] = inputs
    report["eos_listing_ok"] = eos_ok
    report["eos_note"] = eos_note
    report["totals"] = {
        "expected": expected,
        "done": done,
        "failed": len(failed),
        "terminal_failures": len(terminal),
        "in_progress": in_progress,
        "eos_unlistable": unlistable,
        "pct_done": round(100.0 * done / expected, 1) if expected else 0.0,
    }
    report["by_class"] = dict(by_class)
    report["per_sample"] = list(per_sample.values())
    report["chunks"] = [asdict(r) for r in results]
    report["terminal_failures"] = [
        {"key": r.key, "state": r.state, "reason": r.reason} for r in terminal]
    report["next_steps"] = _next_steps(by_class, eos_ok, eos_note, per_sample, inputs)
    return report


def _next_steps(by_class, eos_ok, eos_note, per_sample, inputs):
    """Deterministic 'what to do next' hints from the class histogram."""
    steps = []
    eos_dir = inputs.get("eos_chunk_dir", "the EOS chunk dir")
    if not eos_ok:
        if eos_note == "not-found":
            steps.append(f"{eos_dir} does not exist -- check --eos-chunk-dir (it is a "
                         "bare /store/... path, no root:// prefix). 'done' counts are "
                         "meaningless until it is correct.")
        else:
            steps.append("The EOS chunk dir could not be listed (" + (eos_note or "error") +
                         ") -- confirm EOS health with 'xrdfs root://cmseos.fnal.gov stat "
                         "/store/group/lpcmetx/SIDM' and re-run. 'done' counts are unreliable "
                         "until then.")
    if by_class.get(OOM):
        steps.append(f"{by_class[OOM]} OOM chunk(s): resubmit at higher request_memory.")
    if by_class.get(CODE_BUG):
        steps.append(f"{by_class[CODE_BUG]} CODE_BUG chunk(s): a deterministic python "
                     "error -- fix the processor and rebuild sidm_code.tar.gz; a plain "
                     "retry will re-fail.")
    if by_class.get(HELD_PROXY):
        steps.append(f"{by_class[HELD_PROXY]} HELD_PROXY chunk(s): renew the VOMS proxy "
                     "and copy it to shared NFS (README §11), then resubmit unchanged.")
    if by_class.get(MISSING_ROOT_FILE):
        steps.append(f"{by_class[MISSING_ROOT_FILE]} chunk(s) have input ROOT files "
                     "confirmed absent on EOS -- these prune (a later phase) or shrink "
                     "the dataset; verify against the location YAML.")
    if by_class.get(EOS_TIMEOUT):
        steps.append(f"{by_class[EOS_TIMEOUT]} EOS_TIMEOUT chunk(s): transient. Confirm "
                     "EOS is healthy, then resubmit the whole chunk unchanged.")
    incomplete = [s["sample"] for s in per_sample.values() if not s["complete"]]
    if incomplete:
        steps.append(f"{len(incomplete)} sample(s) incomplete -- DO NOT merge them yet. "
                     "Complete samples are ready to merge.")
    else:
        steps.append("All samples complete -- ready to merge.")
    return steps


# Short, stable labels for the class lines in the rendered report.
_CLASS_NOTE = {
    DONE: "on EOS, size>0",
    ZERO_BYTE_OUTPUT: "0-byte output; rerun",
    MISSING_ROOT_FILE: "input ROOT file confirmed absent",
    EOS_TIMEOUT: "transient EOS read failure",
    OOM: "out of memory",
    HELD_PROXY: "x509 proxy unreadable",
    HELD_TRANSFER: "input transfer failure",
    STALLED: "removed by stall guard",
    REMOVED: "removed manually (condor_rm)",
    CODE_BUG: "deterministic python error -- NOT a retry",
    XRDCP_FAILED: "output copy to EOS failed",
    MISSING_OUTPUT: "exit 0 but no output on EOS",
    UNKNOWN: "unrecognized non-zero exit",
    LOST: "no terminal log event",
    IN_PROGRESS: "queued or still running",
    EOS_UNLISTABLE: "EOS dir unlistable; presence unknown",
}


def render_report(report):
    """Render the human-facing report (the text printed and saved as report.md)."""
    t = report["totals"]
    inp = report["inputs"]
    L = []
    bar = "=" * 64
    L.append(bar)
    L.append(f" SIDM Condor campaign reconcile -- run_id: {report['run_id']}")
    L.append(f" generated {report['generated_utc']}")
    L.append(f" logs: {inp.get('logs_dir')}")
    L.append(f" EOS : {inp.get('eos_chunk_dir')}")
    if not report["eos_listing_ok"]:
        note = report.get("eos_note") or "error"
        L.append(f" !! EOS dir could not be listed ({note}) -- 'done' counts are unreliable")
    L.append(bar)
    L.append("CAMPAIGN TOTALS                                  chunks")
    L.append(f"  expected (from job_args) .................... {t['expected']}")
    L.append(f"  DONE (on EOS, size>0) ....................... {t['done']}   ({t['pct_done']}%)")
    L.append(f"  failed (confirmed) ......................... {t['failed']}")
    L.append(f"  terminal (need a human) .................... {t['terminal_failures']}")
    if t.get("in_progress"):
        L.append(f"  in progress (queued/running) .............. {t['in_progress']}")
    if t.get("eos_unlistable"):
        L.append(f"  presence unknown (EOS unlistable) ......... {t['eos_unlistable']}")
    L.append("")
    L.append("FAILURE CLASSES (confirmed failures)")
    skip = {DONE, IN_PROGRESS, EOS_UNLISTABLE}
    shown = 0
    for cls, n in sorted(report["by_class"].items(), key=lambda kv: -kv[1]):
        if cls in skip:
            continue
        L.append(f"  {cls:<18} {n:>4}  ({_CLASS_NOTE.get(cls, '')})")
        shown += 1
    if shown == 0:
        L.append("  (none)")
    L.append("")
    L.append("PER-SAMPLE COMPLETENESS (incomplete only)")
    any_incomplete = False
    for s in report["per_sample"]:
        if s["complete"]:
            continue
        any_incomplete = True
        miss = ", ".join(str(c) for c in s["missing_chunks"][:12])
        more = "" if len(s["missing_chunks"]) <= 12 else f" (+{len(s['missing_chunks'])-12} more)"
        L.append(f"  {s['sample']:<32} {s['done']}/{s['expected']}  "
                 f"INCOMPLETE  missing chunks: {miss}{more}")
    n_complete = sum(1 for s in report["per_sample"] if s["complete"])
    if not any_incomplete:
        L.append(f"  (all {n_complete} sample(s) complete)")
    else:
        L.append(f"  ({n_complete} sample(s) complete, collapsed)")
    if report["terminal_failures"]:
        L.append("")
        L.append("TERMINAL FAILURES (manual action required)")
        for r in report["terminal_failures"][:40]:
            L.append(f"  {r['key']:<36} {r['state']:<16} {r['reason'][:60]}")
    L.append("")
    L.append("WHAT TO DO NEXT")
    for step in report["next_steps"]:
        L.append(f"  - {step}")
    L.append(bar)
    return "\n".join(L)


def render_status_line(report):
    t = report["totals"]
    return (f"[{report['run_id']}] {t['done']}/{t['expected']} done "
            f"({t['pct_done']}%), {t['failed']} failed, "
            f"{t.get('in_progress', 0)} in-progress, "
            f"{t['terminal_failures']} terminal")
