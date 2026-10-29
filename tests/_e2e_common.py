"""Shared helpers for the per-PR chain report (chain_report.py).

Runs the SidmProcessor over a small committed fixture (a 200-event 2Mu2E signal
file) and exposes the pieces the report compares between `main` and the PR:
  * broken hist collections    -- static: collections referencing undefined hists
  * processor warnings         -- the "Unable to apply ... Skipping" lines, captured
  * per-channel cutflow counts -- raw + weighted cumulative counts (selection regression)
  * channel health             -- which channels are regression-insensitive on the fixture

There is no committed baseline: chain_report.py recomputes both sides from the
checkout it runs against, so there is nothing for these helpers to keep in sync.
"""
import os
import re
import sys
import io
import contextlib

# Heavy deps (coffea, sidm.*) are imported lazily inside the functions that use
# them, so importing this module stays cheap for the static-only helpers.
# pylint: disable=import-outside-toplevel

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

FIXTURE = os.path.join(_HERE, "data", "events_2mu2e_500GeV_200ev.root")
# A signal-shaped dataset name so the processor takes its "assume 1 fb" path
# instead of failing the cross-section lookup.
DATASET = "2Mu2E_500GeV_1p2GeV_1p9mm"

# Real-data fixture: a 200-event slice of a DoubleMuon 2018C skim
# (skimmed_output_583.root, entries 11649-11849) chosen to SPAN a golden-JSON
# dead-run boundary. Events 0-99 are run 320007, which is absent from the 2018
# golden JSON entirely; events 100-199 are run 319993, which is fully golden.
# With CHUNKSIZE=100 the first chunk is 100% non-golden, so any code that
# filters data by the golden JSON produces a ZERO-EVENT chunk here -- the
# configuration that crashed real 2018C processing (PR #42) and that MC
# fixtures can never probe, because the mask only runs when is_data is true.
DATA_FIXTURE = os.path.join(_HERE, "data", "events_DoubleMuon2018C_200ev.root")
DATA_DATASET = "DoubleMuon_2018C_slice"

# The processor applies the standard corrections by default, and their files live
# on CVMFS, which GitHub's runners do not mount. data/corrections holds copies of
# the files the default list and the corrections_* histograms read, in the CVMFS
# layout, so the chain runs with the corrections everywhere. A correction file
# that the processor or a histogram starts reading needs copying here too.
os.environ.setdefault("SIDM_CORRECTIONS_ROOT", os.path.join(_HERE, "data", "corrections"))

# Two chunks per 200-event fixture. For the MC fixture this only splits the
# accumulation (sums are chunk-order independent, so all counts are unchanged);
# for the data fixture it is load-bearing -- see DATA_FIXTURE above.
CHUNKSIZE = 100
# Deliberately non-unit: postprocess scales the weighted column by
# lumi*xs/(sum_w/skim_factor), so at skim_factor=1.0 (with the fixture's
# genWeights identically 1.0) the /skim_factor path is a numerical no-op and a
# PR that broke or dropped the skim_factor handling would leave every weighted
# number unchanged; at 0.5 such a bug moves every weighted row.
SKIM_FACTOR = 0.5

# One channel per selection family, for the histogram-coverage / plotting pass.
# chain_report.py intersects this list with all_channels() at run time, so a
# channel renamed or removed on one side is skipped (and reported) on that side
# rather than crashing it.
REPRESENTATIVE_CHANNELS = [
    "base", "baseNoLj", "2mu2e", "4mu", "base_sr", "2mu2e_sr",
    "gen_leptons_final", "gen_leptons_born", "base_ljObjCut_ljIso",
    "baseNoLj_displacedA", "baseNoLj_A_mumu_matched_lj", "baseNoLj_A_ee_matched_lj",
    "baseNoLj_Two_Ormore_Muon", "data_control_region_1muLj",
]


def _hist_defs():
    from sidm.definitions.hists import hist_defs
    return hist_defs


def _collection_menu():
    from sidm import BASE_DIR
    from sidm.tools import utilities
    return utilities.load_yaml(f"{BASE_DIR}/configs/hist_collections.yaml")


def all_channels():
    """Every runnable selection (has an obj_cuts block); excludes pure cut-fragment
    anchors (e.g. trigger, pv_cuts) that are not standalone channels."""
    from sidm import BASE_DIR
    from sidm.tools import utilities
    sel = utilities.load_yaml(f"{BASE_DIR}/configs/selections.yaml")
    return sorted(k for k, v in sel.items() if isinstance(v, dict) and "obj_cuts" in v)


def broken_collections():
    """{collection: [missing hist names]} for collections referencing hists that are
    not defined in hist_defs (dangling references)."""
    from sidm.tools import utilities
    hist_defs, menu = _hist_defs(), _collection_menu()
    out = {}
    for coll in menu:
        missing = sorted({n for n in utilities.flatten(menu[coll]) if n not in hist_defs})
        if missing:
            out[coll] = missing
    return out


def valid_collections():
    """Collections whose every referenced hist is defined (safe to fill)."""
    from sidm.tools import utilities
    hist_defs, menu = _hist_defs(), _collection_menu()
    return sorted(c for c in menu
                  if all(n in hist_defs for n in utilities.flatten(menu[c])))


def _normalize_warning(line):
    """Strip event/array-size-specific bits so warnings compare stably across runs."""
    line = re.sub(r"length \d+", "length N", line)
    line = re.sub(r"with \[\[.*", "with [...]", line)
    return line.strip()


def run_chain(channels, collections, data=False):
    """Run the processor over a fixture. Returns (per_sample_output, warning_set).
    data=False runs the MC signal fixture; data=True runs the real-data fixture
    (is_data metadata, so data-only code paths like golden-JSON filtering run)."""
    unknown = [c for c in channels if c not in set(all_channels())]
    if unknown:
        raise ValueError(f"Unknown channel(s) not defined as selections with obj_cuts: {unknown}")
    from coffea import processor
    from sidm.tools import sidm_processor, llpnanoaodschema
    if data:
        fileset = {DATA_DATASET: {"files": [DATA_FIXTURE],
                                  "metadata": {"is_data": True, "year": "2018",
                                               "skim_factor": 1.0}}}
    else:
        fileset = {DATASET: {"files": [FIXTURE],
                             "metadata": {"is_data": False, "year": "2018",
                                          "skim_factor": SKIM_FACTOR}}}
    runner = processor.Runner(
        executor=processor.IterativeExecutor(),
        schema=llpnanoaodschema.LLPNanoAODSchema,
        skipbadfiles=False,
        chunksize=CHUNKSIZE,
    )
    proc = sidm_processor.SidmProcessor(channels, collections, verbose=False)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = runner.run(fileset, treename="Events", processor_instance=proc)
    warnings = {_normalize_warning(l) for l in buf.getvalue().splitlines()
                if l.lstrip().startswith("Warning:")}
    return out["out"][DATA_DATASET if data else DATASET], warnings


def cutflow_counts(per_sample_output):
    """{channel: {cut_name: raw_cumulative_count}} for each channel's EVENT-LEVEL cut
    sequence (object-level cuts slim collections but add no cutflow rows). Raw integer
    counts -- deterministic and independent of the weighting."""
    cf = per_sample_output["cutflow"]
    return {ch: {cut: int(cf[ch].rows[cut]["raw"]) for cut in cf[ch].rows}
            for ch in cf}


def cutflow_weighted(per_sample_output):
    """{channel: {cut_name: weighted_cumulative_count}} -- the weighted twin of
    cutflow_counts, rounded to 6 decimals for a stable diff. This is the column that
    moves when weights/normalization change without any raw count changing."""
    cf = per_sample_output["cutflow"]
    return {ch: {cut: round(float(cf[ch].rows[cut]["weighted"]), 6) for cut in cf[ch].rows}
            for ch in cf}


def scaled_sum_weights(per_sample_output):
    """The sample-level sum(genWeight)/skim_factor the processor uses for lumi*xs/Sw
    normalization, rounded for a stable diff. None if absent."""
    val = per_sample_output.get("metadata", {}).get("scaled_sum_weights")
    return None if val is None else round(float(val), 6)


def channel_health(per_sample_output):
    """Which channels' cutflows the fixture cannot regression-test:
    {"zero_final": [...], "no_cut_rows": [...]}. zero_final = the final event-level
    cut leaves 0 events, so downstream count regressions are invisible. no_cut_rows
    = the cutflow holds only the initial "None" total (the channel applies no
    event-level cuts -- object-cut-only and cut-fragment-style channels), so the
    count check sees nothing beyond the event total. Note this is about the CUTFLOW
    only: a channel whose distinguishing cuts are all object-level (e.g. the
    electron-ID family) has generic event rows and is NOT flagged here, but its ID
    content is still invisible to the count check. Uses the in-memory row order
    (the persisted JSON is key-sorted, so the final cut is only known here)."""
    cf = per_sample_output["cutflow"]
    zero_final, no_cut_rows = [], []
    for ch in cf:
        rows = cf[ch].rows
        if set(rows) <= {"None"}:
            no_cut_rows.append(ch)
        elif int(list(rows.values())[-1]["raw"]) == 0:
            zero_final.append(ch)
    return {"zero_final": sorted(zero_final), "no_cut_rows": sorted(no_cut_rows)}
