#!/usr/bin/env python3
"""Standalone regression tests for campaign_lib (no pytest needed).

Run on the LPC submitter (or anywhere with python3):

    python condor/test_campaign_lib.py

Covers the filename/log/EOS parsers and every branch of the failure classifier
with synthetic fixtures, so the classification logic is checked without needing a
live Condor pool or EOS. A few extra checks run against the real condor/logs if
they happen to be present; they are skipped otherwise so the suite stays portable.
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import campaign_lib as cl

_PASS = 0
_FAIL = 0


def check(name, cond):
    global _PASS, _FAIL
    if cond:
        _PASS += 1
    else:
        _FAIL += 1
        print("FAIL:", name)


# --------------------------------------------------------------------------
# Filename + EOS-listing parsers
# --------------------------------------------------------------------------
def test_parsers():
    check("log_name ok",
          cl.parse_log_name("2Mu2E_500GeV_1p2GeV_1p9mm_3_59157219_3.log")
          == ("2Mu2E_500GeV_1p2GeV_1p9mm", 3, 59157219, 3))
    check("log_name underscore-sample",
          cl.parse_log_name("A_B_C_12_345_6.err") == ("A_B_C", 12, 345, 6))
    check("log_name bad -> None",
          cl.parse_log_name("not_a_real_logfile.log") is None)
    check("coffea_stem ok",
          cl.parse_coffea_stem("2Mu2E_500GeV_1p2GeV_1p9mm_7.coffea")
          == ("2Mu2E_500GeV_1p2GeV_1p9mm", 7))
    check("coffea_stem bad -> None",
          cl.parse_coffea_stem("merged_summary.coffea") is None)

    ls = ("-rw-r--r-- murtazas us_cms 21604 2026-05-26 16:48:04 "
          "/store/user/murtazas/sidm_condor/SignalSmokeTest_v3/2Mu2E_500GeV_1p2GeV_1p9mm_0.coffea\n"
          "drwxr-xr-x murtazas us_cms 0 2026-05-26 16:00:00 /store/user/murtazas/sub\n"
          "garbage line with no path\n")
    sizes = cl.parse_ls_l_sizes(ls)
    check("ls size parsed", sizes.get("2Mu2E_500GeV_1p2GeV_1p9mm_0.coffea") == 21604)
    # No date column -> we must NOT guess a digit as the size (avoids a false DONE).
    check("ls no-date -> no size",
          cl.parse_ls_l_sizes("-rw 1700000000 /store/x/foo_0.coffea") == {})

    check("url_to_eos_path",
          cl.url_to_eos_path("root://cmseos.fnal.gov//store/x/y.root") == "/store/x/y.root")
    check("strip_redirector full url",
          cl.strip_redirector("root://cmseos.fnal.gov//store/u/x") == "/store/u/x")
    check("strip_redirector bare path",
          cl.strip_redirector("/store/u/x") == "/store/u/x")


# --------------------------------------------------------------------------
# EOS per-file probe: the ABSENT-only-on-positive-answer invariant (H1)
# --------------------------------------------------------------------------
def test_eos_probe():
    orig = cl._run
    try:
        cl._run = lambda cmd, **k: (-2, "", "[Errno 2] No such file or directory: 'xrdfs'")
        check("xrdfs binary missing -> INCONCLUSIVE (never ABSENT)",
              cl.xrdfs_stat_exists("root://x", "root://x//store/a.root") == cl.INCONCLUSIVE)
        cl._run = lambda cmd, **k: (-1, "", "timeout")
        check("xrdfs timeout -> INCONCLUSIVE",
              cl.xrdfs_stat_exists("root://x", "root://x//store/a.root") == cl.INCONCLUSIVE)
        cl._run = lambda cmd, **k: (54, "", "[ERROR] Server responded: [3011] No such file or directory")
        check("server no-such-file -> ABSENT",
              cl.xrdfs_stat_exists("root://x", "root://x//store/a.root") == cl.ABSENT)
        cl._run = lambda cmd, **k: (0, "", "")
        check("rc 0 -> PRESENT",
              cl.xrdfs_stat_exists("root://x", "root://x//store/a.root") == cl.PRESENT)
    finally:
        cl._run = orig


# --------------------------------------------------------------------------
# ULOG (.log) event parsing against synthetic fixtures
# --------------------------------------------------------------------------
_TERM_LOG = """\
000 (123.000.000) 2026-05-28 09:00:00 Job submitted from host: <1.2.3.4:9618>
...
001 (123.000.000) 2026-05-28 09:00:05 Job executing on host: <1.2.3.4:9618>
\tSlotName: slot1@wn
\tCpus = 1
\tMemory = 4096
...
006 (123.000.000) 2026-05-28 09:01:00 Image size of job updated: 500
\t500  -  MemoryUsage of job (MB)
\t512000  -  ResidentSetSize of job (KB)
...
005 (123.000.000) 2026-05-28 09:05:00 Job terminated.
\t(1) Normal termination (return value 1)
\t\tUsr 0 00:00:00, Sys 0 00:00:00  -  Run Remote Usage
...
"""

_HELD_LOG = """\
000 (456.000.000) 2026-05-28 09:33:00 Job submitted from host: <1.2.3.4:9618>
...
012 (456.000.000) 2026-05-28 09:33:06 Job was held.
\tTransfer input files failure at access point lpcschedd6 ... reading from file /tmp/x509up_u58548: (errno 2) No such file or directory
\tCode 13 Subcode 2
...
009 (456.000.000) 2026-05-28 09:33:44 Job was aborted.
\tvia condor_rm (by user murtazas)
...
"""

# Two intermediate events (001, 006) are MISSING their "..." terminators; only the
# real 005 terminal is well-formed. The parser must resync on each header (M3).
_UNTERMINATED_LOG = """\
000 (789.000.000) 2026-05-28 10:00:00 Job submitted from host: <1.2.3.4:9618>
...
001 (789.000.000) 2026-05-28 10:00:05 Job executing on host: <1.2.3.4:9618>
\tMemory = 4096
006 (789.000.000) 2026-05-28 10:01:00 Image size of job updated: 700
\t700  -  MemoryUsage of job (MB)
005 (789.000.000) 2026-05-28 10:05:00 Job terminated.
\t(1) Normal termination (return value 0)
...
"""


# A job that ran, made progress until 12:04, then hung for ~2 days until LPC's
# site policy reaped it (real shape from a student's log).
_STALL_LOG = """\
000 (3124017.069.000) 2026-06-12 11:43:24 Job submitted from host: <1.2.3.4:9618>
...
001 (3124017.069.000) 2026-06-12 11:44:18 Job executing on host: <1.2.3.4:9618>
\tMemory = 4096
...
006 (3124017.069.000) 2026-06-12 12:04:35 Image size of job updated: 983
\t983  -  MemoryUsage of job (MB)
...
009 (3124017.069.000) 2026-06-14 11:44:49 Job was aborted.
\tLPC job removed by SYSTEM_PERIODIC_REMOVE due to job running for more than 2 days.
...
"""


def _write_log(tmpdir, name, text):
    path = os.path.join(tmpdir, name)
    with open(path, "w") as f:
        f.write(text)
    return path


def test_event_log():
    with tempfile.TemporaryDirectory() as d:
        a = cl.parse_event_log(_write_log(d, "Sig_2_123_0.log", _TERM_LOG))
        check("term sample", a.sample == "Sig" and a.chunk == 2)
        check("term executed", a.executed)
        check("term return_value", a.return_value == 1)
        check("term request_mem", a.request_mem_mb == 4096)
        check("term peak_mem", a.peak_mem_mb == 500)
        check("term terminated", a.terminated and not a.held and not a.aborted)

        h = cl.parse_event_log(_write_log(d, "Sig_0_456_0.log", _HELD_LOG))
        check("held flag", h.held)
        check("held reason x509", "x509up" in h.hold_reason)
        check("held aborted", h.aborted)
        check("held abort_reason", "condor_rm" in h.abort_reason)

        u = cl.parse_event_log(_write_log(d, "Sig_5_789_0.log", _UNTERMINATED_LOG))
        check("resync recovers terminal", u.terminated and u.return_value == 0)
        check("resync recovers exec+mem",
              u.executed and u.request_mem_mb == 4096 and u.peak_mem_mb == 700)

        s = cl.parse_event_log(_write_log(d, "QCD_63_3124017_69.log", _STALL_LOG))
        check("stall parsed aborted", s.aborted and "SYSTEM_PERIODIC_REMOVE" in s.abort_reason)
        check("stall last progress = last 006", s.last_progress_local == "2026-06-12 12:04:35")
        check("stall not OOM (mem under request)", s.peak_mem_mb == 983 and s.request_mem_mb == 4096)


def test_stall_diagnostics():
    check("hours_between ~47.7h",
          abs(cl._hours_between("2026-06-12 12:04:35", "2026-06-14 11:44:49") - 47.67) < 0.1)
    check("hours_between bad -> None", cl._hours_between("nope", "2026-06-14 11:44:49") is None)
    check("idle_note empty for short gap",
          cl._idle_note(cl.AttemptLog("S", 0, 1, 0, "x.log",
                        last_progress_local="2026-06-12 12:00:00",
                        last_event_local="2026-06-12 12:10:00")) == "")
    check("output_tail keeps last n lines", cl._output_tail("a\nb\nc\nd\ne", n=3) == "c\nd\ne")
    check("output_tail skips blank lines", cl._output_tail("a\n\n\nb\n", n=5) == "a\nb")


# --------------------------------------------------------------------------
# Classifier: one synthetic case per branch
# --------------------------------------------------------------------------
SPEC = cl.ChunkSpec("S", 0, "filelists/S_0.txt")
_ORIG_ERR = cl.read_err_tail
_ORIG_STAT = cl.xrdfs_stat_exists


def att(**kw):
    base = dict(sample="S", chunk=0, cluster=1, proc=0, log_path="x.log")
    base.update(kw)
    return cl.AttemptLog(**base)


def classify(attempts, eos_size=None, err="", stat=cl.INCONCLUSIVE,
             stat_missing=True, root_files=None, eos_listable=True):
    # Monkeypatch the err reader and the per-file EOS probe so the classifier is
    # exercised without touching the filesystem or the network.
    cl.read_err_tail = lambda *a, **k: err
    cl.xrdfs_stat_exists = lambda redirector, url, **k: stat
    try:
        return cl.classify_chunk(SPEC, attempts, eos_size, root_files or [],
                                 "root://cmseos.fnal.gov", stat_missing=stat_missing,
                                 eos_listable=eos_listable)
    finally:
        cl.read_err_tail = _ORIG_ERR
        cl.xrdfs_stat_exists = _ORIG_STAT


BAD_URL = "root://cmseos.fnal.gov//store/x/bad_42.root"
OPEN_ERR = ("OSError: XRootD error: [3011]\nFailed to open file "
            + BAD_URL + ": No such file or directory\n")


def test_classifier():
    r = classify([att(terminated=True, return_value=0)], eos_size=21604)
    check("DONE", r.state == cl.DONE)

    r = classify([att(terminated=True, return_value=0)], eos_size=0)
    check("ZERO_BYTE", r.state == cl.ZERO_BYTE_OUTPUT)

    r = classify([att(terminated=True, return_value=0)], eos_size=None)
    check("MISSING_OUTPUT", r.state == cl.MISSING_OUTPUT)

    r = classify([att(terminated=True, return_value=0)], eos_size=None, eos_listable=False)
    check("EOS_UNLISTABLE (exit 0, dir unlistable)", r.state == cl.EOS_UNLISTABLE)

    r = classify([att(terminated=True, return_value=1)],
                 err="Traceback (most recent call last):\nKeyError: 'LJ_pt'\n")
    check("CODE_BUG (reason keeps the key)",
          r.state == cl.CODE_BUG and not r.auto_retryable and "LJ_pt" in r.reason)

    r = classify([att(terminated=True, return_value=1)],
                 err="step 1 of 5\nopening file X\nTraceback (most recent call last):\nKeyError: 'LJ_pt'\n")
    check("CODE_BUG carries multi-line output tail",
          r.state == cl.CODE_BUG and "\n" in r.last_output and "opening file X" in r.last_output)

    r = classify([att(terminated=True, return_value=1)], err=OPEN_ERR, stat=cl.ABSENT,
                 root_files=[BAD_URL])
    check("MISSING_ROOT_FILE (declared input, confirmed absent)",
          r.state == cl.MISSING_ROOT_FILE and r.offending_root_files)

    # The absent URL is NOT one of this chunk's declared inputs -> must not blame it.
    r = classify([att(terminated=True, return_value=1)], err=OPEN_ERR, stat=cl.ABSENT,
                 root_files=["root://cmseos.fnal.gov//store/x/other.root"])
    check("non-input absent -> EOS_TIMEOUT (not MISSING)", r.state == cl.EOS_TIMEOUT)

    r = classify([att(terminated=True, return_value=1)], err=OPEN_ERR,
                 stat=cl.INCONCLUSIVE, root_files=[BAD_URL])
    check("open-fail inconclusive -> EOS_TIMEOUT", r.state == cl.EOS_TIMEOUT)

    r = classify([att(terminated=True, return_value=1)], err=OPEN_ERR,
                 stat=cl.ABSENT, stat_missing=False, root_files=[BAD_URL])
    check("no-stat -> never MISSING", r.state == cl.EOS_TIMEOUT)

    r = classify([att(terminated=True, return_value=1)],
                 err="XRootD error: Operation expired\nValueError: Empty list provided to reduction\n")
    check("EOS_TIMEOUT explicit", r.state == cl.EOS_TIMEOUT)

    r = classify([att(terminated=True, signal=9)])
    check("OOM signal 9", r.state == cl.OOM)

    r = classify([att(terminated=True, return_value=137)])
    check("OOM exit 137", r.state == cl.OOM)

    r = classify([att(terminated=True, return_value=1, request_mem_mb=4096,
                      peak_mem_mb=4000)], err="some unrecognized failure")
    check("mem-pressure -> UNKNOWN (not OOM), with note",
          r.state == cl.UNKNOWN and "possibly memory" in r.reason)

    r = classify([att(terminated=True, return_value=54)],
                 err="xrdcp failed with exit code 54\n")
    check("XRDCP_FAILED", r.state == cl.XRDCP_FAILED)

    r = classify([att(terminated=True, return_value=3)], err="mysterious nonzero exit")
    check("UNKNOWN nonzero", r.state == cl.UNKNOWN)

    r = classify([att(held=True, hold_reason="Transfer input files failure ... /tmp/x509up_u58548 ...")])
    check("HELD_PROXY still-held", r.state == cl.HELD_PROXY)

    r = classify([att(held=True, hold_reason="Transfer input files failure for sidm_code.tar.gz")])
    check("HELD_TRANSFER", r.state == cl.HELD_TRANSFER)

    r = classify([att(aborted=True, abort_reason="SIDM_STALL:WALLCLOCK")])
    check("STALLED (our stall guard)", r.state == cl.STALLED)

    r = classify([att(aborted=True, abort_reason="SYSTEM_PERIODIC_REMOVE evaluated to TRUE")])
    check("system periodic remove -> STALLED retryable",
          r.state == cl.STALLED and r.auto_retryable)

    r = classify([att(aborted=True,
                      abort_reason="LPC job removed by SYSTEM_PERIODIC_REMOVE due to job running for more than 2 days.",
                      last_progress_local="2026-06-12 12:04:35",
                      last_event_local="2026-06-14 11:44:49")])
    check("site-policy stall -> STALLED with idle note",
          r.state == cl.STALLED and "no progress for" in r.reason and r.auto_retryable)

    r = classify([att(aborted=True, abort_reason="via condor_rm (by user murtazas)")])
    check("REMOVED manual", r.state == cl.REMOVED and not r.auto_retryable)

    r = classify([att(cluster=1, held=True, aborted=True,
                      hold_reason="/tmp/x509up_u58548 No such file",
                      abort_reason="via condor_rm (by user murtazas)"),
                  att(cluster=2, aborted=True,
                      abort_reason="via condor_rm (by user murtazas)")])
    check("HELD_PROXY aborted-after-hold", r.state == cl.HELD_PROXY)

    r = classify([], eos_size=None)
    check("LOST no-attempts", r.state == cl.LOST)

    r = classify([att()])  # a log with no recognizable lifecycle
    check("LOST no-lifecycle", r.state == cl.LOST)

    r = classify([att(executed=True)])  # ran, no terminal yet
    check("IN_PROGRESS executed-no-terminal", r.state == cl.IN_PROGRESS)

    # Regression (H2): the chronologically-latest attempt wins even when its
    # ClusterId is numerically SMALLER (cross-schedd retries are not monotonic).
    r = classify([att(cluster=84000000, last_event_local="2026-05-28 09:33:44",
                      held=True, aborted=True,
                      hold_reason="/tmp/x509up_u1 No such file", abort_reason="via condor_rm"),
                  att(cluster=2826673, last_event_local="2026-05-28 18:19:00",
                      terminated=True, return_value=0)], eos_size=None)
    check("H2 newest-by-timestamp wins (lower cluster, later time) -> MISSING_OUTPUT",
          r.state == cl.MISSING_OUTPUT)

    # Regression: EOS present wins even if an earlier attempt failed hard.
    r = classify([att(cluster=1, last_event_local="2026-05-28 09:00:00",
                      terminated=True, return_value=1),
                  att(cluster=2, last_event_local="2026-05-28 10:00:00",
                      terminated=True, return_value=0)], eos_size=999)
    check("EOS-present-wins", r.state == cl.DONE)


# --------------------------------------------------------------------------
# Report building / rendering
# --------------------------------------------------------------------------
def test_report():
    results = [
        cl.ChunkResult("S", 0, "S_0", cl.DONE, n_attempts=1, eos_size=100),
        cl.ChunkResult("S", 1, "S_1", cl.CODE_BUG, reason="KeyError",
                       auto_retryable=False, n_attempts=2),
        cl.ChunkResult("S", 2, "S_2", cl.EOS_TIMEOUT, reason="transient", n_attempts=3),
        cl.ChunkResult("S", 3, "S_3", cl.IN_PROGRESS, reason="running", n_attempts=1),
    ]
    rep = cl.build_report(results, True, {"logs_dir": "x", "eos_chunk_dir": "y"},
                          run_id="t")
    check("report totals", rep["totals"]["expected"] == 4 and rep["totals"]["done"] == 1)
    # IN_PROGRESS is not a failure: 2 confirmed failures (CODE_BUG, EOS_TIMEOUT).
    check("report failed excludes in-progress", rep["totals"]["failed"] == 2)
    check("report in_progress", rep["totals"]["in_progress"] == 1)
    check("report terminal", rep["totals"]["terminal_failures"] == 1)
    check("report per-sample incomplete", rep["per_sample"][0]["complete"] is False)
    text = cl.render_report(rep)
    check("render mentions CODE_BUG", "CODE_BUG" in text)
    check("render hides in_progress from failure list",
          "IN_PROGRESS" not in text.split("FAILURE CLASSES")[1].split("PER-SAMPLE")[0])
    check("status line", "1/4 done" in cl.render_status_line(rep))


def test_real_logs_if_present():
    """Opportunistic checks against the real logs, skipped if absent."""
    real = "/uscms_data/d3/murtazas/SIDM/condor/logs"
    held = os.path.join(real, "2Mu2E_500GeV_1p2GeV_1p9mm_0_84267545_0.log")
    if not os.path.isfile(held):
        return
    a = cl.parse_event_log(held)
    check("real held log: held", a.held)
    check("real held log: x509 reason", "x509up" in a.hold_reason)


if __name__ == "__main__":
    test_parsers()
    test_eos_probe()
    test_event_log()
    test_stall_diagnostics()
    test_classifier()
    test_report()
    test_real_logs_if_present()
    print(f"\n{_PASS} passed, {_FAIL} failed")
    sys.exit(1 if _FAIL else 0)
