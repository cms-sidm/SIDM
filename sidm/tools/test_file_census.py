#!/usr/bin/env python3
"""Standalone tests for file_census pure logic (no network).

    python sidm/tools/test_file_census.py

Covers the YAML enumerator (incl. commented-file attribution + URL building), the
skim-basis measurement, the open-error classifier, and the per-sample rollup
(norm denominator, gen-filter mismatch, no-Runs handling). The actual xrootd probe
needs a live EOS and is exercised by hand, not here.
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import file_census as fc

_P = _F = 0


def check(name, cond):
    global _P, _F
    if cond:
        _P += 1
    else:
        _F += 1
        print("FAIL:", name)


_YAML = """\
verA:
  path: root://xcache//store/base/
  samples:
    SampA:
      files:
      - a_part-0.root
      # - a_part-1.root # LZMA Corrupt input data
      - a_part-2.root
      path: subA/
      skim_factor: 0.5
    SampB:
      files:
      - b_part-0.root
      path: subB/
"""


def test_enumerate():
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as fh:
        fh.write(_YAML)
        path = fh.name
    try:
        e = fc.enumerate_location(path)
        check("total entries", len(e) == 4)
        live = [x for x in e if x.source == "live"]
        comm = [x for x in e if x.source == "commented"]
        check("3 live / 1 commented", len(live) == 3 and len(comm) == 1)
        c = comm[0]
        check("commented attributed to SampA", c.sample == "SampA")
        check("commented filename", c.filename == "a_part-1.root")
        check("commented reason", c.commented_reason == "LZMA Corrupt input data")
        check("xcache rewritten + path joined",
              c.url == "root://cmseos.fnal.gov//store/base/subA/a_part-1.root")
        a0 = next(x for x in live if x.filename == "a_part-0.root")
        check("live url", a0.url == "root://cmseos.fnal.gov//store/base/subA/a_part-0.root")
        check("per-sample skim_factor", a0.skim_factor == 0.5)
        b0 = next(x for x in live if x.sample == "SampB")
        check("SampB default skim_factor", b0.skim_factor == 1.0 and b0.url.endswith("subB/b_part-0.root"))
        check("no-xcache option keeps xcache",
              fc.enumerate_location(path, replace_xcache=False)[0].url.startswith("root://xcache//"))
    finally:
        os.unlink(path)


def test_skim_basis():
    check("post: Events ~ Runs", fc.measure_skim_basis(10000, 9500, 1.0) == "post")
    check("pre: skimmed, Runs >> Events", fc.measure_skim_basis(10000, 130, 0.0132) == "pre")
    check("gen_filtered: no skim, Events << Runs", fc.measure_skim_basis(10000, 940, 1.0) == "gen_filtered")
    check("unknown: no count", fc.measure_skim_basis(0, 5, 1.0) == "unknown")
    check("unknown: no events", fc.measure_skim_basis(10000, None, 1.0) == "unknown")


def test_classify_error():
    check("lzma -> bad", fc._classify_open_error(OSError("LZMA: corrupt input data")) == "bad")
    check("not-a-TFile -> bad", fc._classify_open_error(ValueError("not a TFile")) == "bad")
    check("3011 -> absent", fc._classify_open_error(RuntimeError("[3011] No such file or directory")) == "absent")
    check("timeout -> transient", fc._classify_open_error(TimeoutError("Operation expired")) == "transient")
    check("unknown -> transient (safe default)",
          fc._classify_open_error(ValueError("something weird")) == "transient")


def _row(sample, status, has_runs=True, sw=None, gec=None, ev=None, skim=1.0, source="live"):
    return dict(version="v", sample=sample, skim_factor=skim, is_data=False, source=source,
                status=status, has_runs=has_runs, genEventSumw=sw, genEventCount=gec,
                events_entries=ev)


def test_rollup():
    # gen-filtered sample: 2 reachable-with-Runs + 1 bad
    rows = [_row("S", "reachable", True, 10000.0, 10000, 940),
            _row("S", "reachable", True, 10000.0, 10000, 940),
            _row("S", "bad")]
    s = {x["sample"]: x for x in fc._rollup(rows)}["S"]
    check("rollup reachable count", s["reachable"] == 2 and s["bad"] == 1)
    check("rollup with_runs", s["reachable_with_runs"] == 2)
    check("rollup raw Sw", s["genEventSumw_reachable_raw"] == 20000.0)
    check("rollup gen_filtered basis", s["skim_basis"] == "gen_filtered")
    check("rollup norm = raw (gen_filtered)", s["genEventSumw_reachable_norm"] == 20000.0)
    check("rollup mismatch flag", s["norm_mismatch_warning"] is True)

    # skimmed POST sample -> divide by skim_factor
    rows = [_row("P", "reachable", True, 1000.0, 10000, 9800, skim=0.0132)]
    p = {x["sample"]: x for x in fc._rollup(rows)}["P"]
    check("rollup post basis", p["skim_basis"] == "post")
    check("rollup norm divided by skim", abs(p["genEventSumw_reachable_norm"] - 1000.0 / 0.0132) < 1e-6)
    check("rollup norm_basis_skim", p["norm_basis_skim"] == "divided_by_skim")

    # skimmed sample with NO Runs -> reachable, norm unavailable
    rows = [_row("N", "reachable", has_runs=False, ev=170, skim=0.0132),
            _row("N", "reachable", has_runs=False, ev=160, skim=0.0132)]
    n = {x["sample"]: x for x in fc._rollup(rows)}["N"]
    check("rollup no_runs basis", n["skim_basis"] == "no_runs")
    check("rollup no_runs reachable", n["reachable"] == 2 and n["reachable_no_runs"] == 2)
    check("rollup no_runs norm None", n["genEventSumw_reachable_norm"] is None)

    # empty (0-event) files: reachable, not corrupt, counted as n_empty
    rows = [_row("E", "reachable", has_runs=False, ev=0),
            _row("E", "reachable", has_runs=False, ev=0),
            _row("E", "reachable", has_runs=False, ev=5)]
    e = {x["sample"]: x for x in fc._rollup(rows)}["E"]
    check("rollup n_empty", e["n_empty"] == 2 and e["reachable"] == 3)

    # deep mode: per-file process_status is counted
    rows = [_row("D", "reachable", True, 100.0, 100, 50),
            _row("D", "reachable", True, 100.0, 100, 50), _row("D", "bad")]
    rows[0]["process_status"] = "ok"
    rows[1]["process_status"] = "failed"
    d = {x["sample"]: x for x in fc._rollup(rows)}["D"]
    check("rollup process counts", d["n_process_ok"] == 1 and d["n_process_failed"] == 1)

    # Runs-count anomaly (Events > genEventCount): KEPT as reachable, flagged, never 'bad'
    rows = [_row("A", "reachable", True, 141.0, 141, 284),
            _row("A", "reachable", True, 154.0, 154, 297)]
    for r in rows:
        r["runs_anomaly"] = "Runs counters undercount: Events 2x genEventCount"
    a = {x["sample"]: x for x in fc._rollup(rows)}["A"]
    check("rollup anomaly stays reachable", a["reachable"] == 2 and a["bad"] == 0)
    check("rollup n_runs_anomaly counted", a["n_runs_anomaly"] == 2)


def test_genpart_selfref():
    import awkward as ak
    # event0 clean; event1 has a self-referential mother (idx1 -> 1); event2 a normal record.
    # Pure-array, no network.
    mother = ak.Array([[-1, 0, 1], [-1, 1, 1], [-1, 0, 0]])
    check("genpart self-ref event count", fc._count_selfref_events(mother) == 1)
    clean = ak.Array([[-1, 0, 1], [-1, 0, 0]])
    check("genpart clean -> zero", fc._count_selfref_events(clean) == 0)
    empty = ak.Array([[], []])
    check("genpart empty events -> zero", fc._count_selfref_events(empty) == 0)

    # rollup: a self-ref file stays reachable, is flagged, never 'bad'
    rows = [_row("G", "reachable", True, 100.0, 100, 100),
            _row("G", "reachable", True, 100.0, 100, 100)]
    rows[0]["n_genpart_selfref_events"] = 3
    rows[1]["n_genpart_selfref_events"] = None
    g = {x["sample"]: x for x in fc._rollup(rows)}["G"]
    check("rollup genpart stays reachable", g["reachable"] == 2 and g["bad"] == 0)
    check("rollup n_genpart_anomaly", g["n_genpart_anomaly"] == 1)
    check("rollup genpart_selfref_events total", g["genpart_selfref_events"] == 3)

    # filelists: KEEP by default, DROP only with drop_genpart_corrupt
    manifest = {"meta": {"complete": True}, "files": [
        {"sample": "G", "filename": "f0.root", "url": "root://x//f0.root", "status": "reachable",
         "events_entries": 100, "n_genpart_selfref_events": 3, "commented_reason": "",
         "last_error": ""},
        {"sample": "G", "filename": "f1.root", "url": "root://x//f1.root", "status": "reachable",
         "events_entries": 100, "n_genpart_selfref_events": None, "commented_reason": "",
         "last_error": ""}]}
    d_keep = tempfile.mkdtemp()
    s_keep = fc.cleaned_filelists(manifest, d_keep)
    check("filelists keep genpart by default", s_keep["n_good"] == 2 and s_keep["n_dropped"] == 0)
    d_drop = tempfile.mkdtemp()
    s_drop = fc.cleaned_filelists(manifest, d_drop, drop_genpart_corrupt=True)
    check("filelists drop genpart on request", s_drop["n_good"] == 1 and s_drop["n_dropped"] == 1)
    # _skip.json (the runner veto) lists exactly the dropped files, per sample
    import json as _json
    skip = _json.load(open(os.path.join(d_drop, "_skip.json")))["skip"]
    check("skip.json lists the dropped genpart file", skip == {"G": ["f0.root"]})
    skip_keep = _json.load(open(os.path.join(d_keep, "_skip.json")))["skip"]
    check("skip.json empty when nothing dropped", skip_keep == {})


if __name__ == "__main__":
    test_enumerate()
    test_skim_basis()
    test_classify_error()
    test_rollup()
    test_genpart_selfref()
    print(f"\n{_P} passed, {_F} failed")
    sys.exit(1 if _F else 0)
