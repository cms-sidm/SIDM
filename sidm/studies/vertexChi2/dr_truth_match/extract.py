"""Per-LJ truth records of signal mu-type LJs, input of analyze.py (signal MC only).

usage: python extract.py SAMPLE FILELIST OUT.parquet [--chunksize N] [--workers N]
                         [--no-skip-bad-files]

Runs SidmProcessor over the ROOT files listed in FILELIST with the channel
base_ljDisplacementIso_cosAlpha of configs/selections_vertexChi2_as_run.yaml, no corrections,
unweighted histograms and the hist collection mu_lj_vtx_outcome_lxy_dr (which makes the
processor compute ljs.vtx_chi2), and with one debug branch, lj_rows, that returns one record per
mu-type LJ of the selected events. Writes the records to OUT.parquet and prints the event and
record counts and the vertex-outcome counts of the mu_lj_vtxChi2_muonN_mumu_dR histogram of the
same run, for comparison with the records. Files that cannot be opened are skipped unless
--no-skip-bad-files is given.

Fields of a record: run, lumi, event; lj_pt, lj_eta, lj_phi; n_pf, n_dsa (constituent counts);
vtx_chi2 (as the vertex cut reads it, < 0 for no stored vertex); mumu_dr_ref and gen_lxy_ref
(the values of the histogram helpers mu_lj_mumu_dr and mu_lj_gen_a_lxy); gen_matched (a
genAs_toMu within dR < 0.4 of the LJ) and the pt, eta, phi, Lxy (cm) and dR to the LJ of the
nearest genAs_toMu; n_gen_dau, d1_*, d2_* (pt, eta, phi, charge of the generated daughter
muons of that nearest genAs_toMu, highest pT first, -1 / -99 / 0 when absent); dau_parent_pt
(pt of the dark photon d1 came from); n_other_dau, o1_*, o2_* (the same for the Zd daughter
muons of the event that are not d1, d2); o_a_lxy, o_a_phi (Lxy and phi of the first other
genAs_toMu, -1 / -99 when none); n_genAs_toMu; muons (the PF and DSA constituents,
highest pT first: type 0 = PF and 1 = DSA, pt, eta, phi, charge, idx and, for DSA muons,
nSegments and the PF-muon segment matches muonMatch1..5 / muonMatch1..5idx, -1 for PF
muons; for DSA muons the track reference point vx, vy, vz (cm) and the track-quality branches of
DSA_TRACK, for PF muons the inner-track reference point and momentum and the muon-type flags of
PF_TRACK, -999 for the other type); pv_n, pv_dca, pv_dcaStatus, pv_isValid, pv_normChi2 (the
stored dimuon-vertex rows of the two leading muons m0, m1 in the DSAMuonVertex, PatMuonVertex and
PatDSAMuonVertex tables, legs matched by index and type in either order: their number and the
values of the first, -999 when there is none); sample.
"""
import argparse
import sys

import awkward as ak
import numpy as np
from coffea import processor

from sidm.definitions import hists as H
from sidm.tools import llpnanoaodschema, sidm_processor

CHANNEL = "base_ljDisplacementIso_cosAlpha"
SEL_CFG = "configs/selections_vertexChi2_as_run.yaml"
DSA_MATCH = [f"muonMatch{i}" for i in range(1, 6)] + [f"muonMatch{i}idx" for i in range(1, 6)]
DSA_TRACK = ["vx", "vy", "vz", "normChi2", "chi2", "ndof", "ptErr", "etaErr", "phiErr",
             "nDTSegments", "nCSCSegments", "trkNumHits", "trkNumDTHits", "trkNumCSCHits",
             "trkNumPlanes", "dxy", "dz", "dxyBS", "dxyBSErr", "dzBS", "dzBSErr", "dzPV",
             "dzPVErr", "displacedID", "outerEta", "outerPhi"]
PF_TRACK = ["innerVx", "innerVy", "innerVz", "innerPt", "innerEta", "innerPhi", "isGlobal",
            "isStandalone", "isTracker"]
# (NanoAOD collection, leg type of both legs: 1 = DSA, 0 = PF, None = from isDSAMuon1/2)
VTX_TABLES = (("DSAMuonVertex", 1), ("PatMuonVertex", 0), ("PatDSAMuonVertex", None))
VTX_FIELDS = ("dca", "dcaStatus", "isValid", "normChi2")


def pair_vertex(events, entry, mus):
    """Stored vertex rows of each LJ's two leading muons (m0, m1), from the three dimuon-vertex
    tables, legs matched by muon index and type in either order: the number of rows (pv_n) and
    dca, dcaStatus, isValid, normChi2 of the first (pv_*, -999 when there is none). The tables
    are read from the unselected events of the chunk at the selected events' positions (entry):
    the processor drops them from the selected objects."""
    lead = ak.pad_none(mus, 2, axis=2)
    m0, m1 = lead[:, :, 0], lead[:, :, 1]
    found = {f: [] for f in VTX_FIELDS}
    for name, kind in VTX_TABLES:
        v = events[name][entry]
        i1, i2 = v.originalMuonIdx1[:, None, :], v.originalMuonIdx2[:, None, :]
        if kind is None:
            k1, k2 = v.isDSAMuon1[:, None, :], v.isDSAMuon2[:, None, :]
        else:
            k1 = k2 = kind

        def leg(i, k, m):
            return (i == m.idx[:, :, None]) & (k == m["type"][:, :, None])

        ok = ak.fill_none((leg(i1, k1, m0) & leg(i2, k2, m1))
                          | (leg(i1, k1, m1) & leg(i2, k2, m0)), False)
        for f in VTX_FIELDS:
            found[f].append(ak.broadcast_arrays(v[f][:, None, :], ok)[0][ok])
    out = {}
    for f in VTX_FIELDS:
        allv = ak.concatenate(found[f], axis=2)
        if f == "dca":
            out["pv_n"] = ak.values_astype(ak.num(allv, axis=2), np.int32)
        out[f"pv_{f}"] = ak.fill_none(ak.firsts(allv, axis=2), -999.0)
    return out


def lj_rows(sel_objs, events):
    """One record per mu-type LJ, flattened over events."""
    ljs = sel_objs["mu_ljs"]
    ids = sel_objs["evt_ids"]
    gens = sel_objs["gens"]

    # constituents, PF and DSA pooled, with a type flag (0 = PF, 1 = DSA)
    pf, dsa = ljs.pfMuons, ljs.dsaMuons

    def minus_one(m):
        return ak.values_astype(ak.ones_like(m.pt, dtype=np.int32) * -1, np.int32)

    def fill_value(m):
        return ak.values_astype(ak.ones_like(m.pt, dtype=np.float32) * -999, np.float32)

    def muon_record(m, is_dsa):
        d = {"type": ak.values_astype(ak.ones_like(m.pt, dtype=np.int32) * is_dsa, np.int32),
             "pt": m.pt, "eta": m.eta, "phi": m.phi,
             "charge": ak.values_astype(m.charge, np.int32),
             "idx": ak.values_astype(m.idx, np.int32)}
        for f in DSA_MATCH + ["nSegments"]:
            d[f] = ak.values_astype(m[f], np.int32) if is_dsa else minus_one(m)
        for f in DSA_TRACK:
            d[f] = ak.values_astype(m[f], np.float32) if is_dsa else fill_value(m)
        for f in PF_TRACK:
            d[f] = fill_value(m) if is_dsa else ak.values_astype(m[f], np.float32)
        return ak.zip(d, depth_limit=3)

    mus = ak.concatenate([muon_record(pf, 0), muon_record(dsa, 1)], axis=-1)
    mus = mus[ak.argsort(mus.pt, axis=-1, ascending=False)]

    # matched gen dark photon, as mu_lj_gen_a_lxy finds it
    gen_a, dr_gen = ljs.nearest(sel_objs["genAs_toMu"], return_metric=True)
    matched = ak.fill_none(dr_gen < 0.4, False)

    # final-state daughter muons of a Zd: status-1 muons whose first ancestor with a different
    # pdgId (distinctParent, which skips the muon's own radiation copies) is a Zd (pdgId 32)
    gm = gens[(abs(gens.pdgId) == 13) & (gens.status == 1)]
    gm = gm[ak.fill_none(abs(gm.distinctParent.pdgId) == 32, False)]
    gm = gm[ak.argsort(gm.pt, axis=-1, ascending=False)]
    # per LJ: the daughters of the matched dark photon (their distinctParent is that genAs_toMu
    # entry, compared by pt, eta, phi) and those of the other dark photons of the event
    par = gm.distinctParent
    gmr = ak.zip({"pt": gm.pt, "eta": gm.eta, "phi": gm.phi, "pdgId": gm.pdgId,
                  "ppt": par.pt, "peta": par.eta, "pphi": par.phi})
    gar = ak.zip({"pt": gen_a.pt, "eta": gen_a.eta, "phi": gen_a.phi})
    c = ak.cartesian({"a": gar, "m": gmr}, axis=1, nested=True)
    own = ak.fill_none((c.m.ppt == c.a.pt) & (c.m.peta == c.a.eta) & (c.m.pphi == c.a.phi), False)
    n_dau, n_oth = ak.num(c.m[own], axis=2), ak.num(c.m[~own], axis=2)
    gmp, gmo = ak.pad_none(c.m[own], 2, axis=2), ak.pad_none(c.m[~own], 2, axis=2)
    d1, d2, o1, o2 = gmp[:, :, 0], gmp[:, :, 1], gmo[:, :, 0], gmo[:, :, 1]
    # the first other genAs_toMu entry of the event (in 4Mu the other dark photon)
    ga_all = sel_objs["genAs_toMu"]
    c2 = ak.cartesian({"a": gar, "g": ak.zip({"pt": ga_all.pt, "eta": ga_all.eta,
                                              "phi": ga_all.phi, "lxy": ga_all.lxy})},
                      axis=1, nested=True)
    same = ak.fill_none((c2.g.pt == c2.a.pt) & (c2.g.eta == c2.a.eta) & (c2.g.phi == c2.a.phi),
                        False)
    oth_a = ak.pad_none(c2.g[~same], 1, axis=2)[:, :, 0]

    def bc(x):
        """Broadcast an event-level array to the LJs."""
        return ak.broadcast_arrays(x, ljs.pt)[0]

    out = {
        "run": bc(ids.run), "lumi": bc(ids.lumi), "event": bc(ids.event),
        "lj_pt": ljs.pt, "lj_eta": ljs.eta, "lj_phi": ljs.phi,
        "n_pf": ak.values_astype(ljs.pfMu_n, np.int32),
        "n_dsa": ak.values_astype(ljs.dsaMu_n, np.int32),
        "vtx_chi2": ljs.vtx_chi2,
        "mumu_dr_ref": H.mu_lj_mumu_dr(sel_objs),
        "gen_lxy_ref": H.mu_lj_gen_a_lxy(sel_objs),
        "gen_matched": matched,
        "gen_a_dr_lj": ak.fill_none(dr_gen, -1.0),
        "gen_a_pt": ak.fill_none(gen_a.pt, -1.0),
        "gen_a_eta": ak.fill_none(gen_a.eta, -99.0),
        "gen_a_phi": ak.fill_none(gen_a.phi, -99.0),
        "gen_a_lxy": ak.fill_none(gen_a.lxy, -1.0),
        "n_gen_dau": ak.values_astype(n_dau, np.int32),
        "dau_parent_pt": ak.fill_none(d1.ppt, -1.0),
        "n_genAs_toMu": bc(ak.values_astype(ak.num(sel_objs["genAs_toMu"]), np.int32)),
        "n_other_dau": ak.values_astype(n_oth, np.int32),
        "o_a_lxy": ak.fill_none(oth_a.lxy, -1.0),
        "o_a_phi": ak.fill_none(oth_a.phi, -99.0),
    }
    for k, d in (("d1", d1), ("d2", d2), ("o1", o1), ("o2", o2)):
        out[f"{k}_pt"] = ak.fill_none(d.pt, -1.0)
        out[f"{k}_eta"] = ak.fill_none(d.eta, -99.0)
        out[f"{k}_phi"] = ak.fill_none(d.phi, -99.0)
        out[f"{k}_charge"] = ak.fill_none(ak.values_astype(-np.sign(d.pdgId), np.int32), 0)
    out.update(pair_vertex(events, ids.entry, mus))
    out["muons"] = mus
    rows = ak.zip(out, depth_limit=2)
    return ak.flatten(rows, axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("sample")
    ap.add_argument("filelist")
    ap.add_argument("out")
    ap.add_argument("--chunksize", type=int, default=50_000)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--no-skip-bad-files", dest="skip_bad_files", action="store_false",
                    help="fail on a file that cannot be opened instead of skipping it")
    a = ap.parse_args()
    with open(a.filelist, encoding="utf-8") as fh:
        files = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    if any("/SIDM/Data/" in f for f in files):
        sys.exit("FILELIST holds collision data; this script reads signal MC only")
    fileset = {a.sample: {"files": files, "metadata": {"is_data": False, "skim_factor": 1.0,
                                                        "year": "2018"}}}
    p = sidm_processor.SidmProcessor(
        [CHANNEL], ["mu_lj_vtx_outcome_lxy_dr"], unweighted_hist=True, corrections=None,
        selections_cfg=SEL_CFG, debug=True, debug_branches={"rows": lj_rows},
        include_default_debug_branches=False, debug_suppress_failures=False)
    # run/lumi/event and the position in the chunk ride along with the other event-level
    # objects through the event cuts
    p.obj_defs["evt_ids"] = lambda e: ak.zip({"run": e.run, "lumi": e.luminosityBlock,
                                              "event": e.event,
                                              "entry": ak.local_index(e.run, axis=0)})
    runner = processor.Runner(executor=processor.FuturesExecutor(workers=a.workers),
                              schema=llpnanoaodschema.LLPNanoAODSchema,
                              skipbadfiles=a.skip_bad_files, chunksize=a.chunksize)
    out = runner.run(fileset, treename="Events", processor_instance=p)
    if out.get("exception"):
        sys.exit(f"processor exception: {out['exception']}")
    o = out["out"][a.sample]
    rows = o["debug"]["0.4"][CHANNEL]["rows"].value
    rows = ak.with_field(rows, a.sample, "sample")
    ak.to_parquet(rows, a.out)
    # histogram-side outcome counts of LJs with >= 2 muons, for comparison with the records
    h = o["hists"]["mu_lj_vtxChi2_muonN_mumu_dR"][{"channel": CHANNEL}].values(flow=True)
    print("n_evts", o["metadata"]["n_evts"], "n_rows", len(rows),
          "hist_total", h.sum(), "hist_novtx_ge2", h[[1, 2]][:, 2:].sum(),
          "hist_pass_ge2", h[[3, 4, 5]][:, 2:].sum(), "hist_fail_ge2", h[6:][:, 2:].sum())


if __name__ == "__main__":
    main()
