"""Shared helpers for the per-PR chain report (chain_report.py).

Runs the SidmProcessor over a small committed fixture (a 200-event 2Mu2E signal
file) and exposes the pieces the report compares between `main` and the PR:
  * broken hist collections    -- static: collections referencing undefined hists
  * processor warnings         -- the "Unable to apply ... Skipping" lines, captured
  * per-channel cutflow counts -- raw cumulative event counts (selection regression)

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

# One channel per selection family, for the histogram-coverage / plotting pass.
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


def run_chain(channels, collections):
    """Run the processor over the fixture. Returns (per_sample_output, warning_set)."""
    unknown = [c for c in channels if c not in set(all_channels())]
    if unknown:
        raise ValueError(f"Unknown channel(s) not defined as selections with obj_cuts: {unknown}")
    from coffea import processor
    from sidm.tools import sidm_processor, llpnanoaodschema
    fileset = {DATASET: {"files": [FIXTURE],
                         "metadata": {"is_data": False, "year": "2018", "skim_factor": 1.0}}}
    runner = processor.Runner(
        executor=processor.IterativeExecutor(),
        schema=llpnanoaodschema.LLPNanoAODSchema,
        skipbadfiles=False,
    )
    proc = sidm_processor.SidmProcessor(channels, collections, verbose=False)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = runner.run(fileset, treename="Events", processor_instance=proc)
    warnings = {_normalize_warning(l) for l in buf.getvalue().splitlines()
                if l.lstrip().startswith("Warning:")}
    return out["out"][DATASET], warnings


def cutflow_counts(per_sample_output):
    """{channel: {cut_name: raw_cumulative_count}} -- raw integer counts, which are
    independent of the (1 fb placeholder) cross-section weighting."""
    cf = per_sample_output["cutflow"]
    return {ch: {cut: int(cf[ch].rows[cut]["raw"]) for cut in cf[ch].rows}
            for ch in cf}
