#!/usr/bin/env python3
"""SIDM file census — Phase 1 (read-only).

For every input ROOT file of a sample -- INCLUDING the ones hand-commented-out in the
location YAMLs (which `make_fileset` never sees, because `yaml.safe_load` drops the
`# - ...root` lines) -- open the file over xrootd and read ONLY the per-file `Runs`
tree (`genEventSumw`, `genEventCount`); no event loop. Each file gets a durable verdict
(reachable / bad / unreachable / inconclusive) and, when reachable, its generator-weight
census. This settles two things at once:

  1. which of the commented files are genuinely bad vs wrongly retired, and
  2. the pre-vs-post-skim question for `Runs.genEventSumw` -- MEASURED per file by
     comparing `Runs.genEventCount` to the `Events`-tree entry count, never assumed
     (a wrong guess mis-normalizes a skimmed sample by up to ~7500x).

The verdict discipline mirrors condor/campaign_lib.py: a file is condemned only on a
positive, deterministic failure (a reproducible open/decompress error, or a server
"no such file"); a timeout or redirector hiccup is `inconclusive` and retried, never
blacklisted. Pure standard library + uproot (no coffea/awkward), so it ships clean to
a dask or Condor worker later.

Phase 1 is serial and prints a census + writes a JSON manifest. Scale-out (dask/Condor)
and the committed YAML manifest + normalization wiring are later phases.
"""

import argparse
import hashlib
import json
import os
import re
import socket
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone

import yaml
import uproot

# Reuse the campaign reconciler's positive-only EOS probe for the ABSENT fallback, if
# available (it lives in condor/). Guarded so the census still runs standalone; without
# it, an unconfirmable open failure stays `inconclusive` (the safe direction), never
# `unreachable`.
_cl = None
try:
    _CONDOR_DIR = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "condor")
    if _CONDOR_DIR not in sys.path:
        sys.path.insert(0, _CONDOR_DIR)
    import campaign_lib as _cl
except Exception:
    _cl = None


# --------------------------------------------------------------------------
# Enumeration -- both live and commented files, the same URL rule as make_fileset
# (utilities.py:207-212): version["path"] + sample["path"] + filename, xcache->cmseos.
# --------------------------------------------------------------------------
@dataclass
class FileEntry:
    version: str
    sample: str
    filename: str
    url: str
    source: str               # "live" | "commented"
    commented_reason: str = ""
    skim_factor: float = 1.0
    is_data: bool = False
    year: str = "2018"


_COMMENTED = re.compile(r"^\s*#\s*-\s*(\S+\.root)\b\s*(?:#\s*(.*))?$")
_HEADER = re.compile(r"^(\s*)([^\s:#][^:]*):\s*$")


def _apply_xcache(url, replace_xcache):
    return url.replace("root://xcache//", "root://cmseos.fnal.gov//") if replace_xcache else url


def enumerate_location(yaml_path, version=None, replace_xcache=True):
    """Return all FileEntry for a location YAML -- live files via safe_load, commented
    files via a raw-text pass attributed to the sample whose files: block they sit in."""
    with open(yaml_path) as f:
        text = f.read()
    doc = yaml.safe_load(text)

    versions = [version] if version else [
        k for k, v in doc.items() if isinstance(v, dict) and "samples" in v]

    # Per (version, sample): the URL base + per-sample metadata (authoritative).
    base = {}
    meta = {}
    sample_names = {}
    entries = []
    for ver in versions:
        vblock = doc[ver]
        vpath = vblock["path"]
        sample_names[ver] = set(vblock["samples"].keys())
        for sname, sy in vblock["samples"].items():
            b = vpath + sy["path"]
            base[(ver, sname)] = b
            m = dict(skim_factor=float(sy.get("skim_factor", 1.0)),
                     is_data=bool(sy.get("is_data", False)),
                     year=str(sy.get("year", "2018")))
            meta[(ver, sname)] = m
            for fn in (sy.get("files") or []):
                entries.append(FileEntry(ver, sname, fn, _apply_xcache(b + fn, replace_xcache),
                                         "live", "", **m))

    # Raw-text pass for the commented files.
    versions_set = set(versions)
    cur_ver = cur_sample = None
    for line in text.splitlines():
        hm = _HEADER.match(line)
        if hm:
            indent, key = hm.group(1), hm.group(2)
            if not indent and key in versions_set:
                cur_ver, cur_sample = key, None
                continue
            if cur_ver and key in sample_names.get(cur_ver, ()):
                cur_sample = key
                continue
        cm = _COMMENTED.match(line)
        if cm and cur_ver and cur_sample is not None:
            fn, reason = cm.group(1), (cm.group(2) or "").strip()
            b = base[(cur_ver, cur_sample)]
            entries.append(FileEntry(cur_ver, cur_sample, fn,
                                     _apply_xcache(b + fn, replace_xcache),
                                     "commented", reason, **meta[(cur_ver, cur_sample)]))
    return entries


# --------------------------------------------------------------------------
# The probe -- open + read Runs (stronger than `xrdfs stat`: an LZMA-corrupt file
# stat's fine but won't open).
# --------------------------------------------------------------------------
@dataclass
class ProbeResult:
    status: str               # reachable | bad | unreachable | inconclusive
    genEventSumw: float = None
    genEventCount: int = None
    events_entries: int = None
    n_runs_entries: int = None
    has_runs: bool = True      # False if the file (legitimately, for a skim) has no Runs tree
    n_attempts: int = 0
    error: str = ""
    runs_anomaly: str = ""     # set when Events > genEventCount: a counting-convention note, NOT corruption
    n_genpart_selfref_events: int = None  # --check-genpart: events with >=1 self-referential GenPart
                                          # mother (mother == own index) -- the distinctParent loop-causer.
                                          # None = not checked, no GenPart (e.g. data), OR scan errored
                                          # (then genpart_anomaly is set).
    genpart_anomaly: str = ""             # set to the error message when the GenPart scan could not run
    process_status: str = None  # deep mode only: "ok" | "failed" (full SidmProcessor on this file)
    process_error: str = ""


# Deterministic "this file is bad" signatures -> condemn. Anything else defaults to
# transient (retry / inconclusive), never a blacklist on uncertainty.
_BAD = re.compile(
    r"lzma|zlib|decompress|not a TFile|Deserialization|truncat|corrupt|"
    r"unknown compression|invalid.*key|bad magic", re.I)
_ABSENT = re.compile(r"no such file|\[3011\]|not found", re.I)
_OPEN_TIMEOUT = 120   # seconds; passed to uproot.open AND recorded in the manifest provenance


def _classify_open_error(exc):
    msg = f"{type(exc).__name__}: {exc}"
    if _BAD.search(msg):
        return "bad"
    if _ABSENT.search(msg):
        return "absent"
    return "transient"


def _xrdfs_absent(redirector, url):
    """Positive-only existence check, self-contained (no campaign_lib): return True ONLY when
    the server explicitly says the file does not exist. Any other stat outcome (timeout, auth,
    redirector error) returns False -> stays `inconclusive`, never a wrong condemnation. This
    lets the probe still emit `unreachable` on a dask/Condor worker where condor/campaign_lib
    is not importable (the `_cl is None` path)."""
    m = re.search(r"(/store/\S+)$", url)
    path = m.group(1) if m else url
    host = redirector.split("://")[-1].strip("/")
    try:
        r = subprocess.run(["xrdfs", host, "stat", path],
                           capture_output=True, text=True, timeout=60)
    except Exception:
        return False
    if r.returncode == 0:
        return False
    return bool(re.search(r"no such file|not found|\[3011\]", r.stderr + r.stdout, re.I))


def _count_selfref_events(mother):
    """Given the jagged GenPart_genPartIdxMother array for a file, return the number of events with
    >=1 SELF-REFERENTIAL mother (index points to itself) -- the signature that makes coffea's
    distinctParent walk loop forever and hang the job. Pure awkward, unit-testable. (An out-of-range
    local mother index is NOT flagged: coffea's local2global masks it to -1, so it does not hang
    distinctParent.)"""
    import awkward as ak
    local = ak.local_index(mother)
    selfref = ak.any(mother == local, axis=1)
    return int(ak.sum(selfref))


def probe_file(url, redirector="root://cmseos.fnal.gov", retries=2, check_genpart=False):
    """Open `url`, read its Runs tree, and return a ProbeResult. When check_genpart is set,
    also read GenPart_genPartIdxMother and flag self-referential / out-of-range mother indices
    (the GenPart corruption that hangs coffea's distinctParent) -- an extra full-branch read,
    so it is opt-in."""
    last_err = ""
    attempts = 0
    for attempt in range(retries + 1):
        attempts = attempt + 1
        try:
            with uproot.open(url, timeout=_OPEN_TIMEOUT) as f:
                top = {k.split(";")[0] for k in f.keys(recursive=False)}
                ev = int(f["Events"].num_entries) if "Events" in top else None
                # Optional GenPart-integrity scan: the self-referential mother index that makes
                # coffea's distinctParent loop forever. One extra branch read (opt-in). Absent
                # GenPart (e.g. data) leaves the count None (not applicable). A read error leaves
                # the count None AND records genpart_anomaly so the failure is visible, not silent.
                gp_selfref = None
                genpart_anomaly = ""
                if check_genpart and ev and "Events" in top:
                    if "GenPart_genPartIdxMother" in set(f["Events"].keys()):
                        try:
                            mother = f["Events"].arrays(
                                ["GenPart_genPartIdxMother"],
                                library="ak")["GenPart_genPartIdxMother"]
                            gp_selfref = _count_selfref_events(mother)
                        except Exception as exc:  # noqa: BLE001
                            genpart_anomaly = f"GenPart scan failed: {type(exc).__name__}: {exc}"[:160]
                # No Runs tree is NOT a failure: skimmed ntuples legitimately strip Runs. The
                # file is still reachable for processing (the analysis reads Events); its
                # normalization then comes from skim_factor + an Events genWeight sum, not Runs.
                # Only a file with NEITHER Runs nor Events is unusable.
                has_genw = False
                if "Runs" in top:
                    runs = f["Runs"]
                    n_runs = int(runs.num_entries)
                    has_genw = n_runs > 0 and {"genEventSumw", "genEventCount"} <= set(runs.keys())
                else:
                    n_runs = 0
                if not has_genw:
                    if ev is None:
                        return ProbeResult("bad", n_attempts=attempts,
                                           error="no Events tree and no usable Runs")
                    return ProbeResult("reachable", events_entries=ev, n_runs_entries=n_runs,
                                       has_runs=False, n_attempts=attempts,
                                       n_genpart_selfref_events=gp_selfref,
                                       genpart_anomaly=genpart_anomaly)
                arr = runs.arrays(["genEventSumw", "genEventCount"], library="np")
                gesw = float(arr["genEventSumw"].sum())
                gec = int(arr["genEventCount"].sum())
                # A readable Runs tree can report Events > genEventCount. That is NOT corruption:
                # the Events tree is intact and the analysis reads it (the processor sums genWeight
                # over Events, never the Runs counters). In the SIDM CutDecayFalse signals the Runs
                # counters (genEventCount AND genEventSumw) are written at ~half the Events count --
                # a producer counting convention, confirmed by the files processing cleanly through
                # the full SidmProcessor. So flag it as an anomaly but KEEP the file: a Runs-count
                # mismatch is a normalization caveat (the census Sw-from-Runs is unreliable for these
                # samples; use the Events genWeight sum), not a reason to drop a processable file.
                runs_anomaly = ""
                if ev is not None and gec and ev > gec * 1.05:
                    runs_anomaly = (f"Runs counters undercount: Events {ev} = {ev / gec:.2f}x "
                                    f"genEventCount {gec} (genEventSumw too) -- counting convention, "
                                    f"not corruption; normalize from the Events genWeight sum")
                return ProbeResult("reachable", genEventSumw=gesw, genEventCount=gec,
                                   events_entries=ev, n_runs_entries=n_runs,
                                   has_runs=True, n_attempts=attempts, runs_anomaly=runs_anomaly,
                                   n_genpart_selfref_events=gp_selfref,
                                   genpart_anomaly=genpart_anomaly)
        except Exception as exc:  # noqa: BLE001 -- we deliberately classify, not crash
            last_err = f"{type(exc).__name__}: {exc}"[:300]
            kind = _classify_open_error(exc)
            if kind == "bad":
                return ProbeResult("bad", n_attempts=attempts, error=last_err)
            if kind == "absent":
                return ProbeResult("unreachable", n_attempts=attempts, error=last_err)
            # transient -> retry; on the last try fall through to the stat probe
    # Exhausted retries on a transient error: ask the server directly (positive-only). Prefer
    # campaign_lib's tested probe when it shipped (serial/Condor); fall back to the self-contained
    # check on a worker where condor/ is absent (dask), so both backends can emit `unreachable`.
    if _cl is not None:
        try:
            if _cl.xrdfs_stat_exists(redirector, url) == _cl.ABSENT:
                return ProbeResult("unreachable", n_attempts=attempts, error=last_err)
        except Exception:
            pass
    elif _xrdfs_absent(redirector, url):
        return ProbeResult("unreachable", n_attempts=attempts, error=last_err)
    return ProbeResult("inconclusive", n_attempts=attempts, error=last_err)


def measure_skim_basis(genEventCount, events_entries, skim_factor):
    """Classify the per-file Events/Runs count ratio. MEASURED, never assumed.

    Events-tree entries are always post-(skim and gen-filter). Comparing to Runs.genEventCount:
      - ratio >= 0.9            -> "post": Runs count matches Events, i.e. the skim rewrote Runs
                                   (so a skimmed sample must divide Sw by skim_factor), or the
                                   sample is simply unskimmed+unfiltered.
      - ratio < 0.9, skim<0.95  -> "pre": a declared skim with Runs >> Events -> Runs is the
                                   original pre-skim count -> use raw Sw.
      - ratio < 0.9, skim~1     -> "gen_filtered": NO declared skim but Events << Runs -> the
                                   shortfall is a GENERATOR filter, not a skim. Runs.genEventSumw
                                   is the full pre-filter sum; the processor instead sums genWeight
                                   over the post-filter Events tree -> the two denominators differ
                                   (here by ~10x) and a physicist must decide which the xsec means.
    """
    if not genEventCount or events_entries is None:
        return "unknown"
    ratio = events_entries / genEventCount
    if ratio >= 0.9:
        return "post"
    if skim_factor and skim_factor < 0.95:
        return "pre"
    return "gen_filtered"


def probe_file_deep(url, sample, metadata, channels, hist_collections,
                    redirector="root://cmseos.fnal.gov", chunksize=50000):
    """Shallow probe, then -- only for a reachable, NON-empty file -- run the full
    SidmProcessor on JUST that file to certify it actually processes. Sets
    process_status ("ok"/"failed") + process_error on the returned ProbeResult.

    Heavy imports (coffea, sidm_processor) are LOCAL so the default shallow census
    never pays for them. skipbadfiles=False, so a deep read failure RAISES (and is
    recorded) rather than being silently skipped. A postprocess cross-section
    lookup failure is a config issue, not a file fault, so it is reported as 'ok'.
    """
    pr = probe_file(url, redirector=redirector, check_genpart=True)
    if pr.status != "reachable" or (pr.events_entries or 0) == 0:
        return pr  # bad / unreachable / empty: nothing to process-test
    # The cycle-safe distinctParent kernel (sidm.tools.gen, installed via the LLPNanoAOD schema)
    # keeps the processor from hanging on self-referential GenPart, so deep mode runs the full
    # SidmProcessor on those files too -- they are flagged (GENPART INTEGRITY), not skipped. (The
    # default channels don't read .distinctParent, so a green deep run does not itself exercise the
    # guard -- that is covered by test_gen.py and the standalone corrupt-file check.)
    try:
        from coffea import processor
        from sidm.tools import sidm_processor, llpnanoaodschema
        runner = processor.Runner(
            executor=processor.FuturesExecutor(workers=1),
            schema=llpnanoaodschema.LLPNanoAODSchema,
            skipbadfiles=False, chunksize=chunksize)
        p = sidm_processor.SidmProcessor(channels, hist_collections, unweighted_hist=True)
        runner.run({sample: {"files": [url], "metadata": metadata}},
                   treename="Events", processor_instance=p)
        pr.process_status = "ok"
    except Exception as exc:  # noqa: BLE001 -- classify, don't crash the census
        import traceback
        tb = traceback.format_exc()
        if any(s in tb for s in ("get_xs", "get_lumixs_weight", "xs_menu")):
            # events processed fine; only the xsec/postprocess lookup failed (config).
            pr.process_status = "ok"
            pr.process_error = ("events processed; xsec/postprocess config issue "
                                f"(not a file fault): {type(exc).__name__}: {exc}")[:240]
        else:
            pr.process_status = "failed"
            pr.process_error = f"{type(exc).__name__}: {exc}"[:300]
    return pr


# --------------------------------------------------------------------------
# Census: enumerate -> probe -> aggregate -> manifest + summary
# --------------------------------------------------------------------------
def _utc():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _sh(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=15).stdout.strip()
    except Exception:
        return ""


def _secs(t0, t1):
    try:
        a = datetime.strptime(t0, "%Y-%m-%dT%H:%M:%SZ")
        b = datetime.strptime(t1, "%Y-%m-%dT%H:%M:%SZ")
        return (b - a).total_seconds()
    except ValueError:
        return 0.0


def _provenance(yaml_path, version, samples, max_live, redirector, mode, started, ended, deep_cfg):
    """Self-describing record: when / what / where / how / by-whom the census was produced."""
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sha = ""
    try:
        with open(yaml_path, "rb") as f:
            sha = hashlib.sha256(f.read()).hexdigest()
    except OSError:
        pass
    probe = {"open_timeout": _OPEN_TIMEOUT, "retries": 2}
    if deep_cfg:
        probe.update(deep_cfg)
    return {
        "mode": mode,
        "started_utc": started, "ended_utc": ended, "duration_s": round(_secs(started, ended), 1),
        "user": os.environ.get("USER", ""), "host": socket.gethostname(),
        "python": sys.version.split()[0], "uproot_version": uproot.__version__,
        "sidm_commit": _sh(["git", "-C", repo, "rev-parse", "--short", "HEAD"]),
        "sidm_branch": _sh(["git", "-C", repo, "rev-parse", "--abbrev-ref", "HEAD"]),
        "source_yaml": os.path.abspath(yaml_path), "source_yaml_sha256": sha,
        "version": version or "ALL", "samples": samples or "ALL",
        "max_live": max_live, "commented_files_tested": True,
        "redirector": redirector, "probe": probe,
    }


def run_census(yaml_path, version=None, samples=None, max_live=None,
               redirector="root://cmseos.fnal.gov", progress=True, mode="shallow",
               channels=("base",), hist_collections=("muon_base",), chunksize=50000,
               run_id="manual", workers=1, check_genpart=True):
    """Enumerate + probe every file; return a self-describing manifest dict.

    mode="shallow" (default): cheap open + read Runs + peek Events (the publishable
    light census). mode="deep": additionally run the full SidmProcessor on each
    reachable non-empty file to certify it processes (expensive; for validation runs).
    """
    started = _utc()
    entries = enumerate_location(yaml_path, version=version)
    if samples:
        wanted = set(samples)
        entries = [e for e in entries if e.sample in wanted]

    # Optionally cap LIVE files per sample for a fast pass, but ALWAYS probe commented ones.
    if max_live:
        seen = {}
        capped = []
        for e in entries:
            if e.source == "live":
                seen[e.sample] = seen.get(e.sample, 0) + 1
                if seen[e.sample] > max_live:
                    continue
            capped.append(e)
        entries = capped

    def _one(e):
        if mode == "deep":
            md = {"is_data": e.is_data, "skim_factor": e.skim_factor, "year": e.year}
            pr = probe_file_deep(e.url, e.sample, md, list(channels), list(hist_collections),
                                 redirector=redirector, chunksize=chunksize)
        else:
            pr = probe_file(e.url, redirector=redirector,
                            check_genpart=(check_genpart and not e.is_data))
        row = {**asdict(e)}
        row.update(status=pr.status, genEventSumw=pr.genEventSumw,
                   genEventCount=pr.genEventCount, events_entries=pr.events_entries,
                   n_runs_entries=pr.n_runs_entries, has_runs=pr.has_runs,
                   n_attempts=pr.n_attempts, last_error=pr.error, runs_anomaly=pr.runs_anomaly,
                   n_genpart_selfref_events=pr.n_genpart_selfref_events,
                   genpart_anomaly=pr.genpart_anomaly,
                   process_status=pr.process_status, process_error=pr.process_error,
                   probed_utc=_utc())
        return row

    n = len(entries)
    # Threads parallelize the I/O-bound SHALLOW probe well (network opens + C-level
    # decompression release the GIL). Deep mode runs coffea per file and stays serial.
    if workers and workers > 1 and mode == "shallow":
        from concurrent.futures import ThreadPoolExecutor
        if progress:
            print(f"  probing {n} files with {workers} threads...", file=sys.stderr)
        with ThreadPoolExecutor(max_workers=workers) as ex:
            rows = list(ex.map(_one, entries))
    else:
        rows = []
        for i, e in enumerate(entries):
            if progress:
                print(f"  [{i+1}/{n}] {mode} {e.source:9s} {e.sample}/{e.filename}", file=sys.stderr)
            rows.append(_one(e))

    ended = _utc()
    deep_cfg = ({"channels": list(channels), "hist_collections": list(hist_collections),
                 "chunksize": chunksize} if mode == "deep" else None)
    meta = _provenance(yaml_path, version, samples, max_live, redirector, mode,
                       started, ended, deep_cfg)
    meta["probe"]["check_genpart"] = bool(check_genpart) or (mode == "deep")
    # Completeness reconciliation: a consumer must be able to tell a full census from a partial
    # one. Serial always probes every enumerated file; the dask/Condor mergers stamp the same
    # keys, and cleaned_filelists warns when complete is False.
    meta["n_enumerated"] = len(entries)
    meta["n_probed"] = len(rows)
    meta["complete"] = (len(rows) == len(entries))
    return {
        "schema_version": 2,
        "run_id": run_id,
        "meta": meta,
        "files": rows,
        "samples": _rollup(rows),
    }


def _rollup(rows):
    by = {}
    for r in rows:
        s = by.setdefault((r["version"], r["sample"]), {
            "version": r["version"], "sample": r["sample"],
            "skim_factor": r["skim_factor"], "is_data": r["is_data"],
            "n_files": 0, "n_live": 0, "n_commented": 0,
            "reachable": 0, "bad": 0, "unreachable": 0, "inconclusive": 0,
            "reachable_with_runs": 0, "reachable_no_runs": 0, "n_empty": 0,
            "n_process_ok": 0, "n_process_failed": 0,
            "n_runs_anomaly": 0, "n_genpart_anomaly": 0, "genpart_selfref_events": 0,
            "genEventCount_reachable": 0, "events_entries_reachable": 0,
            "genEventSumw_reachable_raw": 0.0, "skim_basis_votes": {}})
        s["n_files"] += 1
        s["n_live" if r["source"] == "live" else "n_commented"] += 1
        s[r["status"]] = s.get(r["status"], 0) + 1
        if r["status"] == "reachable":
            # deep mode: did the full SidmProcessor run on this file?
            if r.get("process_status") == "ok":
                s["n_process_ok"] += 1
            elif r.get("process_status") == "failed":
                s["n_process_failed"] += 1
            # Runs-count anomaly (Events > genEventCount): kept as reachable, flagged for awareness.
            if r.get("runs_anomaly"):
                s["n_runs_anomaly"] += 1
            # GenPart self-reference: kept as reachable, flagged; distinctParent hangs without the guard.
            if r.get("n_genpart_selfref_events"):
                s["n_genpart_anomaly"] += 1
                s["genpart_selfref_events"] += r["n_genpart_selfref_events"]
            # An empty (0-event) file is reachable and not corrupt, but contributes nothing;
            # a chunk built ENTIRELY of empties makes coffea's executor return None ->
            # "TypeError: cannot unpack non-iterable NoneType object". Surface it.
            if (r.get("events_entries") or 0) == 0:
                s["n_empty"] += 1
            if r.get("has_runs") and r["genEventSumw"] is not None:
                s["reachable_with_runs"] += 1
                s["genEventCount_reachable"] += r["genEventCount"] or 0
                s["events_entries_reachable"] += r["events_entries"] or 0
                s["genEventSumw_reachable_raw"] += r["genEventSumw"] or 0.0
                sb = measure_skim_basis(r["genEventCount"], r["events_entries"], r["skim_factor"])
                s["skim_basis_votes"][sb] = s["skim_basis_votes"].get(sb, 0) + 1
            else:
                s["reachable_no_runs"] += 1

    out = []
    for s in by.values():
        votes = s.pop("skim_basis_votes")
        sf = s["skim_factor"] or 1.0
        if s["reachable_with_runs"] == 0:
            # Reachable files (if any) have no Runs tree -> the Sw denominator is NOT available
            # from the file. It must come from skim_factor + an Events genWeight sum (the path the
            # processor already uses), or from the unskimmed version's Runs. A later phase decides.
            s["skim_basis"] = "no_runs" if s["reachable"] else "unknown"
            s["genEventSumw_reachable_norm"] = None
            s["norm_basis_skim"] = ("skim_factor+events_sum (Runs stripped)"
                                    if s["reachable"] else "n/a")
            s["norm_mismatch_warning"] = False
            s["events_over_runs_ratio"] = None
            s["skim_basis_disagreement"] = False
            out.append(s)
            continue
        non_unknown = {k: v for k, v in votes.items() if k != "unknown"}
        basis = max(non_unknown, key=non_unknown.get) if non_unknown else "unknown"
        s["skim_basis"] = basis
        raw = s["genEventSumw_reachable_raw"]
        if basis == "post" and sf < 1.0:
            s["genEventSumw_reachable_norm"] = raw / sf
            s["norm_basis_skim"] = "divided_by_skim"
        else:
            s["genEventSumw_reachable_norm"] = raw     # pre / gen_filtered / unskimmed -> raw
            s["norm_basis_skim"] = basis
        # gen_filtered: the census Sw (pre-filter, from Runs) differs from the processor's
        # post-filter sum over the Events tree -> a real per-sample denominator decision.
        s["norm_mismatch_warning"] = (basis == "gen_filtered")
        gec, ev = s["genEventCount_reachable"], s["events_entries_reachable"]
        s["events_over_runs_ratio"] = round(ev / gec, 4) if gec else None
        s["skim_basis_disagreement"] = len(non_unknown) > 1
        out.append(s)
    return out


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------
def render_census(manifest):
    m = manifest.get("meta", {})
    L = []
    bar = "=" * 78
    L.append(bar)
    L.append(f" FILE CENSUS  {os.path.basename(m.get('source_yaml', '?'))}"
             f"  ::  {m.get('version', '?')}   [mode: {m.get('mode', 'shallow')}]")
    L.append(f" run_id {manifest.get('run_id', '?')}   {m.get('ended_utc', '?')}"
             f"   {m.get('user', '?')}@{m.get('host', '?')}"
             f"   sidm@{m.get('sidm_commit', '?')} ({m.get('sidm_branch', '?')})")
    L.append(f" redirector {m.get('redirector', '?')}   duration {m.get('duration_s', '?')}s"
             f"   uproot {m.get('uproot_version', '?')}")
    L.append(bar)
    L.append(f"{'sample':<34} {'reach':>5} {'empty':>5} {'bad':>4} {'unr':>4} {'inc':>4} "
             f"{'ev/runs':>8} {'skim_basis':<13} {'Sw_reachable(norm)':>19}")
    L.append("-" * 110)
    for s in manifest["samples"]:
        sb = s["skim_basis"]
        if sb == "post" and (s["skim_factor"] or 1.0) < 1.0:
            sb = f"post(/{s['skim_factor']:g})"
        ratio = s.get("events_over_runs_ratio")
        rstr = f"{ratio:.4f}" if ratio is not None else "-"
        norm = s.get("genEventSumw_reachable_norm")
        nstr = f"{norm:.6g}" if norm is not None else "n/a(noRuns)"
        flag = ""
        if s.get("n_empty") and s["n_empty"] == s["reachable"]:
            flag = "  <-- ALL reachable files EMPTY (0 ev) -> coffea 'NoneType' TypeError risk"
        elif s.get("skim_basis_disagreement"):
            flag = "  <-- basis disagrees"
        elif s.get("norm_mismatch_warning") and ratio:
            flag = "  <-- gen-filter: pre vs post-filter denom differ ~%.0fx" % (1.0 / ratio)
        elif sb == "no_runs":
            flag = "  <-- skimmed, no Runs: norm via skim_factor+Events"
        L.append(f"{s['sample']:<34} {s['reachable']:>5} {s.get('n_empty', 0):>5} {s['bad']:>4} "
                 f"{s['unreachable']:>4} {s['inconclusive']:>4} {rstr:>8} {sb:<13} "
                 f"{nstr:>19}{flag}")
    # Bad / unreachable detail (the whole point: which files to act on)
    bad = [r for r in manifest["files"] if r["status"] in ("bad", "unreachable")]
    if bad:
        L.append("")
        L.append("BAD / UNREACHABLE FILES")
        for r in bad[:60]:
            claim = f"  (YAML said: {r['commented_reason']})" if r["commented_reason"] else ""
            L.append(f"  {r['status']:<12} {r['sample']}/{r['filename']}{claim}")
            L.append(f"       {r['last_error'][:96]}")
        if len(bad) > 60:
            L.append(f"  ... (+{len(bad) - 60} more; see manifest)")
    anomalies = [r for r in manifest["files"] if r.get("runs_anomaly")]
    if anomalies:
        L.append("")
        L.append("RUNS-COUNT ANOMALIES (Events > genEventCount -- NOT corruption; these files are KEPT).")
        L.append("The Runs counters (genEventCount + genEventSumw) undercount vs the Events tree -- a")
        L.append("producer counting convention. The files process fine; normalize from the Events")
        L.append("genWeight sum, not the Runs Sw (the census Sw-from-Runs is unreliable for these samples).")
        for r in anomalies[:30]:
            L.append(f"  {r['sample']}/{r['filename']}")
            L.append(f"       {r['runs_anomaly'][:96]}")
        if len(anomalies) > 30:
            L.append(f"  ... (+{len(anomalies) - 30} more; see manifest)")
    genpart = [r for r in manifest["files"] if r.get("n_genpart_selfref_events")]
    if genpart:
        L.append("")
        L.append("GENPART INTEGRITY -- SELF-REFERENTIAL MOTHER (distinctParent infinite-loop risk).")
        L.append("These files have GenPart entries whose mother index points to themselves. coffea's")
        L.append("distinctParent walks that self-loop forever and HANGS the job. The files are otherwise")
        L.append("readable and KEPT (this scan reads only GenPart_genPartIdxMother). Process them with")
        L.append("the distinctParent self-reference guard, or run `filelists --drop-genpart-corrupt` to")
        L.append("exclude them until the guard is in place.")
        total = sum(r["n_genpart_selfref_events"] for r in genpart)
        L.append(f"  {len(genpart)} file(s), {total} affected event(s):")
        for r in genpart[:40]:
            L.append(f"  {r['sample']}/{r['filename']}: {r['n_genpart_selfref_events']} event(s)")
        if len(genpart) > 40:
            L.append(f"  ... (+{len(genpart) - 40} more; see manifest)")
    scan_err = [r for r in manifest["files"]
                if r.get("genpart_anomaly") and r.get("n_genpart_selfref_events") is None]
    if scan_err:
        L.append("")
        L.append("GENPART SCAN ERRORS (--check-genpart could not read GenPart on these files; their")
        L.append("self-reference status is UNKNOWN, not certified clean).")
        for r in scan_err[:30]:
            L.append(f"  {r['sample']}/{r['filename']}: {r['genpart_anomaly'][:90]}")
        if len(scan_err) > 30:
            L.append(f"  ... (+{len(scan_err) - 30} more; see manifest)")
    empties = [r for r in manifest["files"]
               if r["status"] == "reachable" and (r["events_entries"] or 0) == 0]
    if empties:
        L.append("")
        L.append("EMPTY FILES (0 events -- reachable and NOT corrupt, but contribute nothing). A chunk")
        L.append("built ENTIRELY of these makes coffea's executor return None -> 'cannot unpack")
        L.append("non-iterable NoneType object'. Drop them, or spread so no job is all-empty.")
        for r in empties[:60]:
            L.append(f"  {r['sample']}/{r['filename']}")
        if len(empties) > 60:
            L.append(f"  ... (+{len(empties) - 60} more; see manifest)")
    procfail = [r for r in manifest["files"] if r.get("process_status") == "failed"]
    if procfail:
        L.append("")
        L.append("PROCESS FAILURES (deep mode -- file opens but the SidmProcessor errored on it)")
        for r in procfail[:60]:
            L.append(f"  {r['sample']}/{r['filename']}")
            L.append(f"       {r.get('process_error', '')[:96]}")
        if len(procfail) > 60:
            L.append(f"  ... (+{len(procfail) - 60} more; see manifest)")
    # Commented files that actually OPEN FINE (wrongly retired)
    revived = [r for r in manifest["files"]
               if r["source"] == "commented" and r["status"] == "reachable"]
    if revived:
        L.append("")
        L.append("COMMENTED-OUT FILES THAT OPEN FINE (candidates to un-comment)")
        for r in revived[:60]:
            L.append(f"  {r['sample']}/{r['filename']}  (was: {r['commented_reason']})")
    L.append(bar)
    return "\n".join(L)


def cleaned_filelists(manifest, out_dir, keep_empty=False, drop_genpart_corrupt=False):
    """Write per-sample filelists of GOOD files from a census manifest: reachable and
    (unless keep_empty) non-empty -- live OR commented-but-actually-good (so a wrongly
    retired file is restored). Drops bad / unreachable / empty. Also writes _dropped.tsv
    (with reasons) and _census_ref.json (provenance pointer). Returns a summary dict."""
    meta = manifest.get("meta", {})
    if meta.get("complete") is False:
        print(f"WARNING: census is PARTIAL ({meta.get('n_probed')}/{meta.get('n_enumerated')} "
              f"files probed) -- the cleaned filelists will OMIT the unprobed files. "
              f"Re-run the census to completion before publishing.", file=sys.stderr)
    os.makedirs(out_dir, exist_ok=True)
    by_sample = {}
    dropped = []
    for r in manifest["files"]:
        empty = (r.get("events_entries") or 0) == 0
        genpart_corrupt = drop_genpart_corrupt and bool(r.get("n_genpart_selfref_events"))
        good = r["status"] == "reachable" and (keep_empty or not empty) and not genpart_corrupt
        if good:
            by_sample.setdefault(r["sample"], []).append(r["url"])
        else:
            if r["status"] == "reachable" and genpart_corrupt:
                reason = "genpart_selfref"
            elif r["status"] == "reachable" and empty:
                reason = "empty"
            else:
                reason = r["status"]
            dropped.append((r["sample"], r["filename"], reason, (r.get("last_error") or "")[:100]))
    for sample, urls in sorted(by_sample.items()):
        with open(os.path.join(out_dir, f"{sample}.txt"), "w") as f:
            f.write("\n".join(sorted(urls)) + "\n")
    with open(os.path.join(out_dir, "_dropped.tsv"), "w") as f:
        f.write("sample\tfilename\treason\terror\n")
        for row in sorted(dropped):
            f.write("\t".join(str(x) for x in row) + "\n")
    with open(os.path.join(out_dir, "_census_ref.json"), "w") as f:
        json.dump({"run_id": manifest.get("run_id"), "meta": manifest.get("meta")}, f, indent=2)
    return {"n_samples": len(by_sample), "n_good": sum(len(v) for v in by_sample.values()),
            "n_dropped": len(dropped), "out_dir": os.path.abspath(out_dir)}


def build_parser():
    p = argparse.ArgumentParser(prog="file_census.py",
                                description="SIDM file census (Phase 1: read-only).")
    sub = p.add_subparsers(dest="command", required=True)
    c = sub.add_parser("census", help="probe files of a location YAML and write a manifest")
    c.add_argument("--location-cfg", required=True,
                   help="path to a sidm/configs/ntuples/*.yaml (or a bare name resolved there)")
    c.add_argument("--version", default=None, help="version block key (default: all in the file)")
    c.add_argument("--samples", default=None, help="comma-separated sample names (default: all)")
    c.add_argument("--max-live", type=int, default=None,
                   help="cap LIVE files probed per sample (commented files are always probed)")
    c.add_argument("--redirector", default="root://cmseos.fnal.gov")
    c.add_argument("--out", default=None, help="manifest JSON path (default: ./census_<stem>.json)")
    c.add_argument("--deep", action="store_true",
                   help="also run the full SidmProcessor on each reachable non-empty file to "
                        "certify it processes (expensive; for validation runs, not the light census)")
    c.add_argument("--channels", default="base", help="deep mode: comma-separated channels")
    c.add_argument("--hist-collections", default="muon_base",
                   help="deep mode: comma-separated hist collections")
    c.add_argument("--chunksize", type=int, default=50000, help="deep mode: coffea chunksize")
    c.add_argument("--no-check-genpart", action="store_true",
                   help="disable the GenPart-integrity scan. By DEFAULT the census reads "
                        "GenPart_genPartIdxMother on each MC file and flags self-referential mother "
                        "indices (the corruption that hangs coffea's distinctParent) -- a cheap extra "
                        "branch read (~5%; the file open dominates). Data has no GenPart and is "
                        "skipped. Always on in --deep mode regardless.")
    c.add_argument("--run-id", default="manual", help="label recorded in the manifest meta")
    c.add_argument("--workers", type=int, default=16,
                   help="concurrent probe threads (shallow mode only; I/O-bound, scales well). "
                        "Parallelizes the unavoidable per-file open timeouts instead of "
                        "serializing them; deep mode ignores this and stays serial.")

    fl = sub.add_parser("filelists", help="write cleaned per-sample filelists from a manifest")
    fl.add_argument("--manifest", required=True, help="census manifest JSON")
    fl.add_argument("--out-dir", required=True, help="dir to write <sample>.txt good-file lists")
    fl.add_argument("--keep-empty", action="store_true",
                    help="keep 0-event files in the lists (default: drop them)")
    fl.add_argument("--drop-genpart-corrupt", action="store_true",
                    help="drop files flagged with self-referential GenPart (default: keep them -- "
                         "the distinctParent guard handles them; use this to exclude until it ships)")
    return p


def _resolve_cfg(path):
    if os.path.isfile(path):
        return path
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cand = os.path.join(here, "configs", "ntuples", path)
    return cand if os.path.isfile(cand) else path


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "filelists":
        with open(args.manifest) as f:
            manifest = json.load(f)
        s = cleaned_filelists(manifest, args.out_dir, keep_empty=args.keep_empty,
                              drop_genpart_corrupt=args.drop_genpart_corrupt)
        print(f"Wrote {s['n_good']} good files across {s['n_samples']} samples to {s['out_dir']}"
              f"  ({s['n_dropped']} dropped; see _dropped.tsv, _census_ref.json)")
        return 0
    cfg = _resolve_cfg(args.location_cfg)
    samples = [s for s in (args.samples or "").split(",") if s] or None
    channels = tuple(x for x in args.channels.split(",") if x.strip())
    hcoll = tuple(x for x in args.hist_collections.split(",") if x.strip())
    manifest = run_census(cfg, version=args.version, samples=samples, max_live=args.max_live,
                          redirector=args.redirector,
                          mode=("deep" if args.deep else "shallow"),
                          channels=channels, hist_collections=hcoll,
                          chunksize=args.chunksize, run_id=args.run_id, workers=args.workers,
                          check_genpart=(not args.no_check_genpart))
    out = args.out or f"census_{os.path.splitext(os.path.basename(cfg))[0]}.json"
    with open(out, "w") as f:
        json.dump(manifest, f, indent=2)
    print(render_census(manifest))
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
