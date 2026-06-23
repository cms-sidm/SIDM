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


if __name__ == "__main__":
    test_enumerate()
    test_skim_basis()
    test_classify_error()
    test_rollup()
    print(f"\n{_P} passed, {_F} failed")
    sys.exit(1 if _F else 0)
