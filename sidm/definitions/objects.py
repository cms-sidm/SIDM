"""Define all commonly used objects"""

import awkward as ak
from sidm.tools.utilities import matched
# define helper functions
def pid(part, val):
    return part[abs(part.pdgId) == val]

def toPid(part, val):
    return part[ak.all(abs(part.children.pdgId) == val, axis=-1)]

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

def nEnPhoton(lj, n_e, n_photon):
    noMu = (lj.muon_n == 0)
    nE = (lj.electron_n == n_e)
    nP = (lj.photon_n == n_photon)

    return lj[noMu & nE & nP]





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
preLj_objs["weight"]     = lambda evts: evts.genWeight
preLj_objs["gens"]       = lambda evts: evts.GenPart
preLj_objs["genMus"]     = lambda evts: pid(preLj_objs["gens"](evts), 13)
preLj_objs["genEs"]      = lambda evts: pid(preLj_objs["gens"](evts), 11)
preLj_objs["genAs"]      = lambda evts: pid(preLj_objs["gens"](evts), 32)
preLj_objs["genAs_toMu"] = lambda evts: toPid(preLj_objs["genAs"](evts), 13)
preLj_objs["genAs_toE"]  = lambda evts: toPid(preLj_objs["genAs"](evts), 11)
preLj_objs["rho_PFIso"]  = lambda evts: evts.fixedGridRhoFastjetAll
preLj_objs["jets"]       = lambda evts: evts.Jet



# define objects whose that will be added to objs by the sidm_processor after LJs are clustered
# and LJ cuts are applied. postLj_obj cuts can be applied to these
postLj_objs = {}
postLj_objs["mu_ljs"]       = lambda objs: yesMu(objs["ljs"])
postLj_objs["egm_ljs"]      = lambda objs: noMu(objs["ljs"])
postLj_objs["pfmu_ljs"]     = lambda objs: noDsa(objs["mu_ljs"])
postLj_objs["dsamu_ljs"]    = lambda objs: noPf(objs["mu_ljs"])
postLj_objs["electron_ljs"] = lambda objs: noPhoton(objs["egm_ljs"])
postLj_objs["photon_ljs"]   = lambda objs: noE(objs["egm_ljs"])
postLj_objs["one_photon_ljs"] = lambda objs: nEnPhoton(objs["ljs"], 0, 1)
postLj_objs["two_photon_ljs"] = lambda objs: nEnPhoton(objs["ljs"], 0, 2)
postLj_objs["one_electron_ljs"] = lambda objs: nEnPhoton(objs["ljs"], 1, 0)
postLj_objs["two_electron_ljs"] = lambda objs: nEnPhoton(objs["ljs"], 2, 0)
postLj_objs["one_e_one_p_ljs"] = lambda objs: nEnPhoton(objs["ljs"], 1, 1)


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
derived_objs['genAs_toE_matched_one_photon_ljs'] = lambda objs, r: matched(objs["genAs_toE"], objs["one_photon_ljs"], r)
derived_objs['genAs_toE_matched_two_photon_ljs'] = lambda objs, r: matched(objs["genAs_toE"], objs["two_photon_ljs"], r)
derived_objs['genAs_toE_matched_one_electron_ljs'] = lambda objs, r: matched(objs["genAs_toE"], objs["one_electron_ljs"], r)
derived_objs['genAs_toE_matched_two_electron_ljs'] = lambda objs, r: matched(objs["genAs_toE"], objs["one_electron_ljs"], r)
derived_objs['genAs_toE_matched_one_e_one_p_ljs'] = lambda objs, r: matched(objs["genAs_toE"], objs["one_e_one_p_ljs"], r)
derived_objs['genAs_toE_matched_photons'] = lambda objs, r: matched(objs['genAs_toE'], objs['photons'], r)
derived_objs['genAs_toE_matched_electrons'] = lambda objs, r: matched(objs['genAs_toE'], objs['electrons'], r)
derived_objs['genA_egmLj_oneEoneP_ptRatio_PS'] = lambda objs: objs["one_e_one_p_ljs"].pt / objs["one_e_one_p_ljs"].nearest(objs["genAs_toE"], threshold=0.4).pt


