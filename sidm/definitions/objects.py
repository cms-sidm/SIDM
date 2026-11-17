"""Define all commonly used objects"""

import awkward as ak
import numpy as np
from sidm.tools.utilities import matched, get_pairs, dR

# define helper functions
def pid(part, val):
    return part[abs(part.pdgId) == val]

def toPid(part, val):
    return part[ak.all(abs(part.children.pdgId) == val, axis=-1)]

def fromPid(part, val):
    return part[abs(part.parent.pdgId) == val]

def yesMu(lj):
    return lj[lj.muon_n > 0]

def noMu(lj):
    return lj[lj.muon_n == 0]

def noDsa(lj):
    return lj[lj.dsaMu_n == 0]

def noPf(lj):
    return lj[lj.pfMu_n == 0]

def noE(lj):
    return lj[lj.electron_n == 0]

def noPhoton(lj):
    return lj[lj.photon_n == 0]

def nE(lj, n):
    return lj[lj.electron_n == n]

def nPhoton(lj, n):
    return lj[lj.photon_n == n]

# These are diagnostic views, not an exhaustive partition: mixed PF/DSA and
# mixed electron/photon jets belong to their parent but neither pure child.
lj_parents = {
    "mu_ljs": "ljs",
    "egm_ljs": "ljs",
    "pfmu_ljs": "mu_ljs",
    "dsamu_ljs": "mu_ljs",
    "electron_ljs": "egm_ljs",
    "photon_ljs": "egm_ljs",
}
lj_candidate_names = ("candidate_ljs", "candidate_mu_ljs", "candidate_egm_ljs")

def lj_membership(ljs, selected):
    """Per-event membership by the stable reconstruction index, never rank/pT."""
    return ak.fill_none(ak.any(
        ljs.idx[:, :, None] == selected.idx[:, None, :], axis=-1
    ), False)

def refresh_lj_children(objs, parent):
    """Prune existing descendants after a parent cut, preserving child cuts.

    Child views do not select their parent. In particular, their union must not
    replace it, since the pure-child categories omit mixed-constituent jets.
    Missing indices are a coherence error and must propagate to the caller.
    """
    for child, source in lj_parents.items():
        if source == parent and child in objs:
            objs[child] = objs[child][lj_membership(objs[child], objs[parent])]
            refresh_lj_children(objs, child)

def add_lj_candidates(objs):
    """Add the leading pair of eligible LJs after all typed selections.

    ``ljs`` is the pT-ordered baseline after generic LJ cuts; it remains intact
    for reconstruction diagnostics. Keep jets surviving their typed selection
    (with >=2 muons for muon candidates), then take the leading two in baseline
    order. Later jets can replace rejected jets. Return a complete pair or an
    empty list; extra jets are allowed. Call after all subtype cuts.
    """
    ljs = objs["ljs"]
    mu_candidates = objs["mu_ljs"][objs["mu_ljs"].muon_n >= 2]
    egm_candidates = objs["egm_ljs"][objs["egm_ljs"].muon_n == 0]
    surviving = (lj_membership(ljs, mu_candidates)
                 | lj_membership(ljs, egm_candidates))
    pair = ljs[surviving][:, :2]
    pair = pair[ak.broadcast_arrays(ak.num(pair) == 2, pair.idx)[0]]
    selected = objs.copy()
    selected["candidate_ljs"] = pair
    selected["candidate_mu_ljs"] = pair[pair.muon_n >= 2]
    selected["candidate_egm_ljs"] = pair[pair.muon_n == 0]
    return selected

def withMass(part, mass):
    return ak.zip(
        {
            "pt": part.pt,
            "eta": part.eta,
            "phi": part.phi,
            "mass": ak.full_like(part.pt, mass),
        },
        with_name="PtEtaPhiMLorentzVector",
    )

# define objects whose definitions don't depend on LJs
preLj_objs = {}
preLj_objs["pvs"]        = lambda evts: evts.PV
preLj_objs["bs"]         = lambda evts: evts.BS
preLj_objs["met"]        = lambda evts: evts.MET
preLj_objs["hlt"]        = lambda evts: evts.HLT
preLj_objs["electrons"]  = lambda evts: evts.Electron
preLj_objs["photons"]    = lambda evts: evts.Photon
preLj_objs["muons"]      = lambda evts: evts.Muon
preLj_objs["dsaMuons"]   = lambda evts: evts.DSAMuon
# uncut copy of the DSA muon collection. No selection lists it under obj_cuts, so it survives
# the LJ-source cross-cleaning that removes DSA muons sharing a PF muon, and histograms can
# still reach that discarded population at fill time.
preLj_objs["dsaMuonsRaw"] = lambda evts: evts.DSAMuon
preLj_objs["weight"]     = lambda evts: evts.genWeight
preLj_objs["gens"]       = lambda evts: evts.GenPart
preLj_objs["genMus"]     = lambda evts: pid(preLj_objs["gens"](evts), 13)
preLj_objs["genEs"]      = lambda evts: pid(preLj_objs["gens"](evts), 11)
preLj_objs["genAs"]      = lambda evts: pid(preLj_objs["gens"](evts), 32)
preLj_objs["genAs_toMu"] = lambda evts: toPid(preLj_objs["genAs"](evts), 13)
preLj_objs["genAs_toE"]  = lambda evts: toPid(preLj_objs["genAs"](evts), 11)
preLj_objs["rho_PFIso"]  = lambda evts: evts.fixedGridRhoFastjetAll
preLj_objs["jets"]       = lambda evts: evts.Jet
preLj_objs["flags"]       = lambda evts: evts.Flag
preLj_objs["bjets"]      = lambda evts: evts.Jet
# simulation only (absent from data): the pileup truth that the pileup weights use
preLj_objs["pileup"]     = lambda evts: evts.Pileup
# the event number seeds the reproducible JER smearing
preLj_objs["event_number"] = lambda evts: evts.event
# dimuon-vertex tables. They carry vx/vy/vxy but no pt/px/py, so the
# processor's automatic beamspot-dxy step skips them (it requires pt)
preLj_objs["dsaMuonVertex"]    = lambda evts: evts.DSAMuonVertex
preLj_objs["patMuonVertex"]    = lambda evts: evts.PatMuonVertex
preLj_objs["patDsaMuonVertex"] = lambda evts: evts.PatDSAMuonVertex
# loaded only when a requested channel or histogram needs them (see SidmProcessor.__init__)
on_demand_preLj_objs = ("dsaMuonsRaw", "dsaMuonVertex", "patMuonVertex", "patDsaMuonVertex")
# uncut copies of the photon and muon collections for the lepton-jet control samples
# (sidm/studies/systematics/lj_reco_efficiency.md): their probes are found among all
# reconstructed objects, while the lepton jets are built from the cut collections. Loaded only
# when a requested channel or histogram needs them (SidmProcessor.needs_probe_inputs)
preLj_objs["photonsRaw"] = lambda evts: evts.Photon
preLj_objs["muonsRaw"]   = lambda evts: evts.Muon
on_demand_probe_objs = ("photonsRaw", "muonsRaw")
# define objects whose that will be added to objs by the sidm_processor after LJs are clustered
# and LJ cuts are applied. postLj_obj cuts can be applied to these
postLj_objs = {}
postLj_objs_MC = {}
postLj_objs["mu_ljs"]       = lambda objs: yesMu(objs["ljs"])
postLj_objs["egm_ljs"]      = lambda objs: noMu(objs["ljs"])
postLj_objs["pfmu_ljs"]     = lambda objs: noDsa(objs["mu_ljs"])
postLj_objs["dsamu_ljs"]    = lambda objs: noPf(objs["mu_ljs"])
postLj_objs["electron_ljs"] = lambda objs: noPhoton(objs["egm_ljs"])
postLj_objs["photon_ljs"]   = lambda objs: noE(objs["egm_ljs"])
postLj_objs["allMuons"] = lambda objs: ak.concatenate([ ak.with_name(objs["muons"],"PtEtaPhiMLorentzVector"),ak.with_name(objs["dsaMuons"],"PtEtaPhiMLorentzVector")],axis=1)
# Adding the following here since I want the cuts on genMus and genEs to be applied
postLj_objs_MC["genMus_fromA"] = lambda objs: fromPid(objs["genMus"], 32)
postLj_objs_MC["genEs_fromA"]  = lambda objs: fromPid(objs["genEs"],  32)

# define objects that depend on extra parameters determined in hist or cut definitions
derived_objs = {}
derived_objs["n_electron_ljs"] = lambda objs, n: nE(objs["electron_ljs"], n)
derived_objs["n_photon_ljs"]   = lambda objs, n: nPhoton(objs["photon_ljs"], n)
derived_objs["genAs_matched_lj"]        = lambda objs, r: matched(objs["genAs"], objs["ljs"], r)
derived_objs["genAs_toMu_matched_lj"]   = lambda objs, r: matched(objs["genAs_toMu"], objs["ljs"], r)
derived_objs["genAs_toE_matched_lj"]    = lambda objs, r: matched(objs["genAs_toE"], objs["ljs"], r)
derived_objs["genAs_matched_muLj"]      = lambda objs, r: matched(objs["genAs"], objs["mu_ljs"], r)
derived_objs["genAs_toMu_matched_muLj"] = lambda objs, r: matched(objs["genAs_toMu"], objs["mu_ljs"], r)
derived_objs["genAs_matched_egmLj"]     = lambda objs, r: matched(objs["genAs"], objs["egm_ljs"], r)
derived_objs["genAs_toE_matched_egmLj"] = lambda objs, r: matched(objs["genAs_toE"], objs["egm_ljs"], r)
derived_objs["mu_lj_matched_genAs_toMu"]   = lambda objs, r: matched(objs["mu_ljs"], objs["genAs_toMu"], r)
derived_objs["back_to_back_dsa_pairs"]   = lambda objs: (lambda pairs, v1, v2: pairs[np.cos(v1.deltaangle(v2)) <= -0.99])(get_pairs(objs["dsaMuons"]),*ak.unzip(get_pairs(objs["dsaMuons"])))
derived_objs["parallel_dsa_pairs"]   = lambda objs: (lambda pairs, v1, v2: pairs[np.cos(v1.deltaangle(v2)) >= 0.99])(get_pairs(objs["dsaMuons"]),*ak.unzip(get_pairs(objs["dsaMuons"])))

# Gen-level objects that depend on PIDs not present in all samples (signal-only).
# Defined as derived_objs so they're only evaluated when explicitly referenced by a histogram or cut.
# Derive from objs["gens"] (uncut) rather than objs["genMus"]/objs["genEs"] (which are channel-filtered
# by status, breaking the fromPid(., 32) link to the Zd parent).
derived_objs["genA_from_genMus"] = lambda objs: withMass(fromPid(pid(objs["gens"], 13), 32), 0.105658).sum()
derived_objs["genA_from_genEs"]  = lambda objs: withMass(fromPid(pid(objs["gens"], 11), 32), 0.000511).sum()
derived_objs["genBSs"]           = lambda objs: pid(objs["gens"], 35)
derived_objs["genBSs_toA"]       = lambda objs: toPid(derived_objs["genBSs"](objs), 32)
derived_objs["genBS_from_genAs"] = lambda objs: pid(objs["gens"], 32).sum()
# final-state generator photons above 5 GeV, to tell real photons from electrons
derived_objs["genPhotons"] = lambda objs: (lambda g: g[(g.status == 1) & (g.pt > 5)])(pid(objs["gens"], 22))


# Probes of the lepton-jet control samples (sidm/studies/systematics/lj_reco_efficiency.md).
# Each probe carries "passing": whether a selected lepton jet of the type it tests lies within
# dR < 0.2. The probe definitions do not depend on the lepton jets, so the efficiency is the
# fraction of passing probes.
def _alone(probes, objs):
    """The probes with no lepton jet of the signal selection's kinds (in ljs, so pT > 30 GeV and
    |eta| < 2.4, with two or more muons or none) farther than dR = 0.4 from them. An event of the
    signal region holds two such lepton jets, at least one of them away from any probe, so none
    of its probes is kept and region A is never read; the condition is the same for passing and
    failing probes."""
    ljs = objs["ljs"]
    kinds = ljs[(ljs.muon_n >= 2) | (ljs.muon_n == 0)]
    pairs = ak.cartesian({"p": probes, "lj": kinds}, nested=True)
    return probes[~ak.any(pairs.p.delta_r(pairs.lj) > 0.4, axis=-1)]


def zmmg_probe_photons(objs, alone=False):
    """Photons from Z -> mumu gamma, the probes of the egamma-type lepton jets.

    The event has exactly two selected PF muons, of opposite charge, with 35 < m(mumu) < 85 GeV
    (final-state radiation lowers it). A probe is any reconstructed photon with pT > 20 GeV and
    |eta| < 2.4, more than 0.5 in dR from both muons, so that the lepton-jet clustering
    (R = 0.4) keeps it apart from them, and with 60 < m(mumu gamma) < 120 GeV. The field m_mmg
    holds m(mumu gamma). With alone, only the probes without another lepton jet (_alone).
    """
    muons = objs["muons"]
    pairs = ak.combinations(muons, 2, fields=["mu0", "mu1"])
    # one pair in events with exactly two muons, none otherwise
    pairs = pairs[ak.broadcast_arrays(ak.num(muons) == 2, pairs.mu0.pt)[0]]
    trip = ak.cartesian({"mm": pairs, "g": objs["photonsRaw"]})
    dimuon = trip.mm.mu0 + trip.mm.mu1
    m_mmg = (dimuon + trip.g).mass
    keep = ((trip.mm.mu0.charge != trip.mm.mu1.charge)
            & (dimuon.mass > 35) & (dimuon.mass < 85)
            & (trip.g.pt > 20) & (abs(trip.g.eta) < 2.4)
            & (trip.g.delta_r(trip.mm.mu0) > 0.5) & (trip.g.delta_r(trip.mm.mu1) > 0.5)
            & (m_mmg > 60) & (m_mmg < 120))
    probes = trip.g[keep]
    probes["m_mmg"] = m_mmg[keep]
    probes["passing"] = dR(probes, objs["egm_ljs"]) < 0.2
    return _alone(probes, objs) if alone else probes


def jpsi_probe_dimuons(objs, alone=False, pt_max=None):
    """Muon pairs near the J/psi mass, the probes of the muon-type lepton jets.

    Any two tracker muons of the uncut collection, of opposite charge and with |eta| < 2.4,
    with 2.6 < m < 3.6 GeV (the peak and its sidebands), dR < 0.4 (one lepton-jet cone) and
    pT(mumu) > 20 GeV, below pt_max if given. No other mass is kept. The signal's dark photon
    masses (0.25, 1.2 and 5 GeV) are outside the window, but mismeasured pairs of displaced
    5 GeV dark photons reach it at high pT, which the safety check of the study measures.
    Fields: the pair's four-vector, dR, and min_pt, the pT of the softer muon. With alone, only
    the probes without another lepton jet (_alone).
    """
    muons = objs["muonsRaw"]
    muons = muons[muons.isTracker & (abs(muons.eta) < 2.4)]
    pairs = ak.combinations(muons, 2, fields=["mu0", "mu1"])
    dimuon = pairs.mu0 + pairs.mu1
    dr = pairs.mu0.delta_r(pairs.mu1)
    keep = ((pairs.mu0.charge != pairs.mu1.charge) & (dimuon.mass > 2.6) & (dimuon.mass < 3.6)
            & (dr < 0.4) & (dimuon.pt > 20))
    if pt_max is not None:
        keep = keep & (dimuon.pt < pt_max)
    probes = dimuon[keep]
    probes["dR"] = dr[keep]
    probes["min_pt"] = np.minimum(pairs.mu0.pt, pairs.mu1.pt)[keep]
    probes["passing"] = dR(probes, objs["mu_ljs"]) < 0.2
    return _alone(probes, objs) if alone else probes


derived_objs["zmmg_probe_photons"] = zmmg_probe_photons
derived_objs["jpsi_probe_dimuons"] = jpsi_probe_dimuons
