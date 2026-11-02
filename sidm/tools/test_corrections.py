#!/usr/bin/env python3
"""Standalone test for sidm.tools.corrections and the weight-variation plumbing

    python sidm/tools/test_corrections.py

The first part builds its own correctionlib JSON and config in a temporary
directory, so it needs no network and no CVMFS: the event-weight arithmetic, the
object SFs and their working-point logic, the cutflow variations and the
histogram syst axis. The second part reads the pinned POG files and is skipped
when /cvmfs is not mounted.
"""
import os
import sys
import tempfile

import awkward as ak
import hist
import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from sidm.tools import corrections  # noqa: E402
from sidm.tools.cutflow import SimpleCutflow  # noqa: E402
from sidm.tools.histogram import Axis, Histogram  # noqa: E402

_P = _F = 0


def check(name, cond):
    """Count a passed or failed check, printing the failures"""
    global _P, _F
    if cond:
        _P += 1
    else:
        _F += 1
        print("FAIL:", name)


def test_event_weights():
    """Products, shifts, subsets and input checks of EventWeights"""
    w = corrections.EventWeights(3)
    w.add("genWeight", [1.0, -2.0, 1.0])
    w.add("pileup", [1.0, 0.0, 2.0], up=[1.1, 0.0, 2.2], down=[0.9, 0.0, 1.8])
    w.add("sf", ak.Array([0.9, 1.0, 1.0]), up=[1.0, 1.1, 1.0], down=[0.8, 0.9, 1.0])
    check("nominal is the product", np.allclose(w.weight(), [0.9, -0.0, 2.0]))
    check("up shifts one factor", np.allclose(w.weight("pileupUp"), [0.99, 0.0, 2.2]))
    check("down shifts one factor", np.allclose(w.weight("sfDown"), [0.8, 0.0, 2.0]))
    # a zero nominal factor must not make the other variations NaN
    check("no NaN from a zero factor", np.all(np.isfinite(w.weight("sfUp"))))
    check("include", np.allclose(w.weight(include=["genWeight"]), [1.0, -2.0, 1.0]))
    check("exclude", np.allclose(w.weight(exclude=["pileup"]), [0.9, -2.0, 1.0]))
    check("variations listed", w.variations == ["pileupUp", "pileupDown", "sfUp", "sfDown"])
    c = w.copy()
    c.add("extra", [2.0, 2.0, 2.0])
    check("copy is independent", "extra" not in w.names and np.allclose(c.weight(), 2 * w.weight()))
    for bad, kwargs in [("up only", {"up": [1, 1, 1]}), ("wrong length", {})]:
        try:
            w.add(f"bad_{bad.replace(' ', '_')}", [1.0] * (3 if bad == "up only" else 2), **kwargs)
            check(f"rejects {bad}", False)
        except ValueError:
            check(f"rejects {bad}", True)
    try:
        w.weight("pileupSideways")
        check("rejects unknown variation", False)
    except KeyError:
        check("rejects unknown variation", True)


def _write_fixture(tmp):
    """A two-bin muon SF and a matching config, in the POG file layout"""
    corrections.import_correctionlib()
    import correctionlib.schemav2 as schema  # pylint: disable=import-outside-toplevel

    def sf_node(nominal):
        # SF = nominal, up/down = nominal +- 0.1, binned in |eta| and pt
        return schema.Binning(nodetype="binning", input="eta", edges=[0.0, 1.2, 2.4], flow="error",
                              content=[schema.Binning(nodetype="binning", input="pt",
                                                      edges=[20.0, 50.0, "inf"], flow="error",
                                                      content=[schema.Category(
                                                          nodetype="category", input="scale_factors",
                                                          content=[{"key": "nominal", "value": v},
                                                                   {"key": "systup", "value": v + 0.1},
                                                                   {"key": "systdown", "value": v - 0.1}])
                                                          for v in row])
                                       for row in nominal])

    corr = schema.Correction(
        name="TOY_ID", version=1, output={"name": "weight", "type": "real"},
        inputs=[{"name": "eta", "type": "real"}, {"name": "pt", "type": "real"},
                {"name": "scale_factors", "type": "string"}],
        data=schema.Transform(nodetype="transform", input="eta",
                              rule=schema.Formula(nodetype="formula", expression="abs(x)",
                                                  parser="TFormula", variables=["eta"]),
                              content=sf_node([[0.95, 0.97], [0.90, 0.92]])))
    cset = schema.CorrectionSet(schema_version=2, corrections=[corr])
    with open(os.path.join(tmp, "toy.json"), "w", encoding="utf-8") as f:
        f.write(cset.model_dump_json(exclude_unset=True))

    cfg = {
        "json_root": tmp,
        "default": ["toy_id"],
        "2018": {
            "toy_id": {
                "kind": "object", "collection": "pfMuons", "cut_object": "muons",
                "working_points": {"looseID": [{"file": "toy.json", "correction": "TOY_ID",
                                                "pt": [20.0, float("inf")]}]},
                "variations": {"nominal": "nominal", "up": "systup", "down": "systdown"},
                "max_abs_eta": 2.4,
            },
        },
    }
    path = os.path.join(tmp, "corrections.yaml")
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f)
    return path


def test_object_sfs():
    """Object SFs, their pt range, working points and per-event products"""
    with tempfile.TemporaryDirectory() as tmp:
        corr = corrections.Corrections("2018", True, cfg=_write_fixture(tmp))
        muons = ak.Array([
            [{"pt": 30.0, "eta": 0.5}, {"pt": 60.0, "eta": -2.4}],   # |eta| = 2.4 exactly: clipped
            [{"pt": 10.0, "eta": 0.1}],                              # below the measured pt: SF 1
            [],
        ])
        sf = corr.object_sf("toy_id", muons, "looseID")
        check("object SF values", ak.to_list(sf) == [[0.95, 0.92], [1.0], []])
        check("object SF up", np.allclose(ak.flatten(corr.object_sf("toy_id", muons, "looseID", "up")),
                                          [1.05, 1.02, 1.0]))

        check("wp from cuts", corr.working_point("toy_id", {"muons": ["looseID", "pT > 5 GeV"]}) == "looseID")
        check("no wp without the cut", corr.working_point("toy_id", {"muons": ["pT > 5 GeV"]}) is None)

        # lepton jets: events -> LJs -> constituents, as built by the processor
        ljs = ak.Array([
            [{"pfMuons": [{"pt": 30.0, "eta": 0.5}]}, {"pfMuons": [{"pt": 60.0, "eta": 2.0}]}],
            [{"pfMuons": []}],
        ])
        w = corrections.EventWeights(2)
        w.add("genWeight", [1.0, 1.0])
        applied = corr.add_object_weights(w, ljs, {"muons": ["looseID"]})
        check("applied reported", applied == {"toy_id": "looseID"})
        check("per-event product", np.allclose(w.weight(), [0.95 * 0.92, 1.0]))
        check("correlated up", np.allclose(w.weight("toy_idUp"), [1.05 * 1.02, 1.0]))
        w2 = corrections.EventWeights(2)
        check("not applied without the cut", not corr.add_object_weights(w2, ljs, {"muons": []}))
        check("provenance lists the file", "toy.json" in corr.provenance()["corrections"]["toy_id"]["files"])

        # selected objects instead of the constituents, as with sf_objects
        selected = {"muons": ak.Array([[{"pt": 30.0, "eta": 0.5}], [{"pt": 30.0, "eta": 1.5}]])}
        w3 = corrections.EventWeights(2)
        applied = corr.add_object_weights(w3, ljs, {"muons": ["looseID"]}, selected)
        check("SFs of the selected objects", applied == {"toy_id": "looseID"}
              and np.allclose(w3.weight(), [0.95, 0.90]))
        w4 = corrections.EventWeights(2)
        check("other collections skipped",
              not corr.add_object_weights(w4, ljs, {"muons": ["looseID"]}, {"electrons": selected["muons"]}))


def test_evaluate_structure():
    """evaluate() keeps flat and jagged shapes"""
    with tempfile.TemporaryDirectory() as tmp:
        corr = corrections.Corrections("2018", True, cfg=_write_fixture(tmp))
        c = corr.correction("toy.json", "TOY_ID")
        jagged = corrections.evaluate(c, ak.Array([[0.5, 1.5], []]), ak.Array([[30.0, 30.0], []]), "nominal")
        check("jagged in, jagged out", ak.to_list(jagged) == [[0.95, 0.90], []])
        flat = corrections.evaluate(c, ak.Array([0.5]), ak.Array([60.0]), "systup")
        check("flat in, flat out", ak.to_list(flat) == [1.07])
        check("numpy passthrough", np.allclose(c.evaluate(np.array([0.5]), np.array([30.0]), "nominal"), [0.95]))


def test_cutflow_variations():
    """SimpleCutflow carries, adds and scales the variations"""
    a = SimpleCutflow()
    a.add_row("None", 10, 10.0, {"pileupUp": 11.0})
    a.add_row("cut", 5, 5.0, {"pileupUp": 6.0})
    b = a.identity()
    check("identity keeps variation names", b.rows["cut"]["variations"] == {"pileupUp": 0.0})
    total = a + a
    check("add sums variations", total.rows["cut"]["variations"]["pileupUp"] == 12.0)
    total.scale(0.5)
    check("scale scales variations", total.variation_yields() == {"nominal": 5.0, "pileupUp": 6.0})
    plain = SimpleCutflow()
    plain.add_row("None", 3, 3.0)
    check("no variations: nominal only", plain.variation_yields() == {"nominal": 3.0})
    check("old rows unchanged", plain.rows["None"] == {"raw": 3, "weighted": 3.0})


def test_hist_syst_axis():
    """Histogram fills one syst category per variation given"""
    h = Histogram([Axis(hist.axis.Regular(4, 0, 4, name="x"), lambda objs, mask: objs["x"])])
    h.make_hist("x", channels=["ch"], systematics=["nominal", "sfUp", "sfDown"])
    check("syst axis after channel", [a.name for a in h.hist.axes] == ["channel", "syst", "x"])
    objs = {"ch": "ch", "x": ak.Array([[0.5, 1.5], [2.5]])}
    w = ak.Array([1.0, 2.0])
    h.fill(objs, {"nominal": w, "sfUp": 2 * w})
    view = h.hist[{"channel": "ch"}]
    check("nominal filled", np.allclose(view[{"syst": "nominal"}].values(), [1, 1, 2, 0]))
    check("variation filled", np.allclose(view[{"syst": "sfUp"}].values(), [2, 2, 4, 0]))
    check("missing category stays empty", view[{"syst": "sfDown"}].sum().value == 0)

    plain = Histogram([Axis(hist.axis.Regular(4, 0, 4, name="x"), lambda objs, mask: objs["x"])])
    plain.make_hist("x", channels=["ch"])
    plain.fill(objs, {"nominal": w, "sfUp": 2 * w})
    check("no syst axis: nominal used", np.allclose(plain.hist[{"channel": "ch"}].values(), [1, 1, 2, 0]))


def test_processor_config():
    """The processor applies the default list unless told otherwise; sf_objects is read per channel"""
    from sidm.tools.sidm_processor import SidmProcessor  # pylint: disable=import-outside-toplevel
    p = SidmProcessor(["Z_to_mumu_control_region", "base"], [])
    check("default list by default", p.corrections == corrections.resolve_names(True))
    check("None turns them off", SidmProcessor(["base"], [], corrections=None).corrections == [])
    cuts = p.build_cuts()
    check("sf_objects read", cuts["Z_to_mumu_control_region"]["sf_objects"] == ["muons"])
    check("lepton jets without sf_objects", cuts["base"]["sf_objects"] is None)


def test_pog_files():
    """Reads the pinned POG files; skipped without CVMFS"""
    corr = corrections.Corrections("2018")
    if not os.path.isdir(corr.json_root):
        print(f"skipping POG checks: {corr.json_root} is not mounted")
        return
    for name in corr.config:
        if isinstance(corr.config[name], dict):
            for f in corr.files(name):
                check(f"{name}: {f} exists", os.path.exists(corr.path(f)))
    pu = corr.correction(corr.config["pileup"]["file"], corr.config["pileup"]["correction"])
    nominal = pu.evaluate(np.array([10.0, 30.0, 50.0]), "nominal")
    check("pileup weights positive", np.all(nominal > 0))
    check("pileup up and down differ", not np.allclose(pu.evaluate(30.0, "up"), pu.evaluate(30.0, "down")))
    muons = ak.Array([[{"pt": 5.0, "eta": 0.3}, {"pt": 25.0, "eta": -1.5}, {"pt": 500.0, "eta": 2.4}]])
    sf = ak.flatten(corr.object_sf("muon_id", muons, "looseID"))
    check("muon SFs near 1 across the J/psi-Z split", np.all(np.abs(np.asarray(sf) - 1) < 0.05))
    electrons = ak.Array([[{"pt": 15.0, "eta": 0.3, "deltaEtaSC": 0.01}, {"pt": 45.0, "eta": -2.0, "deltaEtaSC": 0.0}]])
    sf = ak.flatten(corr.object_sf("electron_reco", electrons, "always"))
    check("electron reco SFs near 1", np.all(np.abs(np.asarray(sf) - 1) < 0.1))
    photons = ak.Array([[{"pt": 40.0, "eta": 0.3, "isScEtaEB": True}, {"pt": 40.0, "eta": 2.0, "isScEtaEB": False}]])
    sf = np.asarray(ak.flatten(corr.object_sf("photon_pixel_veto", photons, "pixelSeed")))
    check("pixel veto SF differs between barrel and endcap", sf[0] != sf[1])


if __name__ == "__main__":
    test_event_weights()
    test_object_sfs()
    test_evaluate_structure()
    test_cutflow_variations()
    test_hist_syst_axis()
    test_processor_config()
    test_pog_files()
    print(f"{_P} passed, {_F} failed")
    sys.exit(1 if _F else 0)
