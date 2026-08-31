"""Best within-LJ dimuon-vertex normChi2, per lepton jet.

For each mu-type LJ, find every stored dimuon vertex (DSAMuonVertex,
PatMuonVertex, PatDSAMuonVertex) whose BOTH legs are constituents of that LJ,
and return the smallest normChi2 among them; -1.0 if no such vertex exists.

SidmProcessor fills ljs["vtx_chi2"] from the three vertex tables registered in
sidm/definitions/objects.py whenever a requested channel cuts on vtx_chi2 or a
requested histogram reads it; otherwise the LJs carry no vtx_chi2 field.

Standalone usage on NanoEvents (LLPNanoAOD; the vertex tables are present in
the 2018 signal files and in the data/background skims):

    ljs["vtx_chi2"] = lj_best_vertex_chi2(
        ljs.pfMuons.idx, ljs.dsaMuons.idx,
        [(evts.DSAMuonVertex, "dsa"), (evts.PatMuonVertex, "pat"),
         (evts.PatDSAMuonVertex, "mix")])

    # keep-cut: rejects LJs whose muons never form a converged dimuon vertex,
    # the signature of fake or cosmic muon pairings
    keep = (ljs.vtx_chi2 >= 0) & (ljs.vtx_chi2 < 5)

Notes:
- Only defined for LJs with >= 2 muons. A single-muon LJ can never have a
  within-LJ vertex, so it is always -1.0; pair any vtx_chi2 cut or histogram
  with "Mu >= 2".
- -1.0 means "no stored within-LJ vertex". A stored normChi2 is returned as
  it is: a real vertex with a huge normChi2 keeps its value (values above 1e7
  occur) and lands in the histogram overflow, and a negative one, should the
  producer write one, comes back negative; the vtx_chi2 cuts treat every
  negative value as "no vertex". Being the minimum, a negative row also masks
  any positive normChi2 of another vertex in the same LJ.
- Unlike the *Spread_* variables, which fill 0 (= passes) when a muon category
  is absent, -1.0 FAILS the keep-cuts by design, so the veto rejects both
  bad-fit LJs and no-vertex LJs.
- If the input sample lacks the vertex tables entirely, SidmProcessor sets
  vtx_chi2 to NaN and every vtx_chi2 cut (keep AND inverse) then selects
  nothing; an exactly empty channel is the symptom, so check the worker logs
  for "dsaMuonVertex not found". The samples in sidm/configs/ntuples/
  (signal_2mu2e_v10, signal_4mu_v10, data_skimmed, backgrounds) all carry the
  tables.
- The producer (cms-sidm/LLPNanoAOD, MuonVertexTableProducer) only stores
  pairs whose Kalman fit converged and with DCA < 15 cm, so "no within-LJ
  vertex" is itself informative: fakes rarely have one. The isValid filter
  below is a defensive re-check of the same requirement.
- Vertex-table columns are float-typed (including indices and flags); the
  membership test compares float indices, exact for these values.
- The same muon appears in many vertex rows (all pairings are fitted); the
  min over matched rows makes that harmless.
- If your LJ objects do not carry constituent .idx, pass any two per-LJ
  index lists (PF-muon and DSA-muon indices of the LJ constituents).
- Validated by the brute-force closure test in
  tests/test_lj_vertex_chi2.py, which also checks the PF-leg index convention
  on the committed NanoAOD test files.
"""
import numpy as np
import awkward as ak


def _member(vidx, pool):
    # vidx: (evt, nvtx) vertex-leg indices; pool: (evt, nlj, nmu) LJ constituents
    # -> (evt, nlj, nvtx): is this leg one of this LJ's constituents?
    return ak.any(vidx[:, None, :, None] == pool[:, :, None, :], axis=-1)


def lj_best_vertex_chi2(pf_pool, dsa_pool, tables):
    """pf_pool/dsa_pool: per-LJ constituent index lists, e.g. ljs.pfMuons.idx.
    tables: list of (vertex_collection, kind) with kind in {"dsa","pat","mix"}."""
    best = None
    for vtx, kind in tables:
        i1, i2 = vtx.originalMuonIdx1, vtx.originalMuonIdx2
        if kind == "dsa":
            ok = _member(i1, dsa_pool) & _member(i2, dsa_pool)
        elif kind == "pat":
            ok = _member(i1, pf_pool) & _member(i2, pf_pool)
        else:  # mixed table: leg type from isDSAMuon1/2
            ok1 = ak.where(vtx.isDSAMuon1[:, None, :] == 1,
                           _member(i1, dsa_pool), _member(i1, pf_pool))
            ok2 = ak.where(vtx.isDSAMuon2[:, None, :] == 1,
                           _member(i2, dsa_pool), _member(i2, pf_pool))
            ok = ok1 & ok2
        ok = ok & (vtx.isValid[:, None, :] == 1)
        chi2 = ak.where(ok, vtx.normChi2[:, None, :], np.inf)
        b = ak.fill_none(ak.min(chi2, axis=-1), np.inf)
        best = b if best is None else np.minimum(best, b)
    return ak.where(np.isinf(best), -1.0, best)
