"""Define all available cuts"""

# columnar analysis
import awkward as ak
import numpy as np
# local
from sidm.definitions.objects import derived_objs
from sidm.tools.utilities import dR, lxy, rho, check_bits, returnBitMapTArrayPhoton, dR_outer, cosAlpha, nearest_lj_index

obj_cut_defs = {
    "pvs": {
        "ndof > 4": lambda objs: objs["pvs"].ndof > 4,
        "|z| < 24 cm": lambda objs: abs(objs["pvs"].z) < 24,
        "|rho| < 0.02 cm": lambda objs: rho(objs["pvs"], ref=objs["bs"]) < 0.02,
    },
    "ljs": {
        "pT > 10 GeV": lambda objs: objs["ljs"].pt > 10,
        "pT > 20 GeV": lambda objs: objs["ljs"].pt > 20,
        "pT > 30 GeV": lambda objs: objs["ljs"].pt > 30,
        "pT > 40 GeV": lambda objs: objs["ljs"].pt > 40,
        "pT > 50 GeV": lambda objs: objs["ljs"].pt > 50,
        "pT > 60 GeV": lambda objs: objs["ljs"].pt > 60,
        "pT < 150 GeV": lambda objs: objs["ljs"].pt < 150,
        "|eta| < 2.4": lambda objs: abs(objs["ljs"].eta) < 2.4,
        "barrel": lambda objs: abs(objs["ljs"].eta) < 1.479,
        "endcap": lambda objs: ((abs(objs["ljs"].eta) > 1.479)& (abs(objs["ljs"].eta) < 2.4)),
        "mu_charge == 0": lambda objs: ak.sum(objs["ljs"].muons.charge, axis= -1) == 0,
        "dR(LJ, A) < 0.2": lambda objs: dR(objs["ljs"], objs["genAs"]) < 0.2,
        "egmLj": lambda objs: ak.num(objs["ljs"].muons) == 0,
        "eLj": lambda objs: (ak.num(objs["ljs"].muons) == 0) & (ak.num(objs["ljs"].electrons) > 0) & (ak.num(objs["ljs"].photons) == 0),
        "gLj": lambda objs: (ak.num(objs["ljs"].muons) == 0) & (ak.num(objs["ljs"].electrons) == 0) & (ak.num(objs["ljs"].photons) > 0),
        "eAndGLj": lambda objs: (ak.num(objs["ljs"].muons) == 0) & (ak.num(objs["ljs"].electrons) > 0) & (ak.num(objs["ljs"].photons) > 0),
        "muLj": lambda objs: ak.num(objs["ljs"].muons) > 0,
        "pfMuLj": lambda objs: (ak.num(objs["ljs"].pfMuons) > 0) & (ak.num(objs["ljs"].dsaMuons) == 0),
        "dsaMuLj": lambda objs: ak.num(objs["ljs"].dsaMuons) > 0,
        "2dsaMuLj": lambda objs: ak.num(objs["ljs"].dsaMuons) > 1,
        "pfDsaMuLj": lambda objs: (ak.num(objs["ljs"].pfMuons) > 0) & (ak.num(objs["ljs"].dsaMuons) > 0),
        "2muLj": lambda objs: objs["ljs"].muon_n >= 2,
        "lj_iso < 0.1": lambda objs: objs["ljs"].isolation < 0.1,
        "lj_iso < 0.2": lambda objs: objs["ljs"].isolation < 0.2,
    },
    "egm_ljs": {
        "eLj": lambda objs: (objs["egm_ljs"].electron_n > 0) & (objs["egm_ljs"].photon_n == 0),
        "gLj": lambda objs: (objs["egm_ljs"].electron_n == 0) & (objs["egm_ljs"].photon_n > 0),
        "egLj": lambda objs: (objs["egm_ljs"].electron_n > 0) & (objs["egm_ljs"].photon_n > 0),
        "1eLj": lambda objs: (objs["egm_ljs"].electron_n == 1) & (objs["egm_ljs"].photon_n == 0),
        "2eLj": lambda objs: (objs["egm_ljs"].electron_n == 2) & (objs["egm_ljs"].photon_n == 0),
        "1gLj": lambda objs: (objs["egm_ljs"].electron_n == 0) & (objs["egm_ljs"].photon_n == 1),
        "2gLj": lambda objs: (objs["egm_ljs"].electron_n == 0) & (objs["egm_ljs"].photon_n == 2),
        "egm_lj_iso < 0.2": lambda objs: objs["egm_ljs"].isolation < 0.2,
        "reverse_iso": lambda objs: objs["egm_ljs"].isolation >= 0.2,
        "lostHits >= 1": lambda objs: ak.min(objs["egm_ljs"].electrons.trkNumPixelHits, axis=-1) >= 1,
        "displaced": lambda objs: (ak.min(objs["egm_ljs"].egamma.lostHits, axis=-1) >= 1),
        "reverse_displaced": lambda objs: (ak.max(objs["egm_ljs"].egamma.lostHits, axis=-1) <1),
    },
    "mu_ljs": {
        "pfMuLj": lambda objs: (objs["mu_ljs"].pfMu_n > 0) & (objs["mu_ljs"].dsaMu_n == 0),
        "dsaMuLj": lambda objs: (objs["mu_ljs"].pfMu_n == 0) & (objs["mu_ljs"].dsaMu_n > 0),
        "pf_dsa_muLj": lambda objs: (objs["mu_ljs"].pfMu_n > 0) & (objs["mu_ljs"].dsaMu_n > 0),
        "1dsaMuLj": lambda objs: objs["mu_ljs"].dsaMu_n > 0,
        "mu_lj_iso < 0.1": lambda objs: objs["mu_ljs"].isolation < 0.1,
        "reverse_iso": lambda objs: objs["mu_ljs"].isolation >= 0.1,
        "pf_pixelhits <= 2": lambda objs: ak.max(objs["mu_ljs"].pfMuons.trkNumPixelHits, axis=-1) <= 2,
        "displaced": lambda objs: (ak.max(objs["mu_ljs"].muons.trkNumPixelHits, axis=-1) <= 2),
        "Mu == 1": lambda objs: objs["mu_ljs"].muon_n == 1,
        "Mu == 2": lambda objs: objs["mu_ljs"].muon_n == 2,
        "Mu == 3": lambda objs: objs["mu_ljs"].muon_n == 3,
        "Mu >= 4": lambda objs: objs["mu_ljs"].muon_n >= 4,
        "Mu >= 2": lambda objs: objs["mu_ljs"].muon_n >= 2,
        "reverse_displaced": lambda objs: (ak.min(objs["mu_ljs"].muons.trkNumPixelHits, axis=-1) > 2),
        "dzSpread_dsa <= 50": lambda objs: objs["mu_ljs"].dzSpread_dsa <= 50,
        "vxySpread_dsa <= 50": lambda objs: objs["mu_ljs"].vxySpread_dsa <= 50,
        "dzSpread_dsa <= 5": lambda objs: objs["mu_ljs"].dzSpread_dsa <= 5,
        "vxySpread_dsa <= 5": lambda objs: objs["mu_ljs"].vxySpread_dsa <= 5,
        "dzSpread_pf <= 5": lambda objs: objs["mu_ljs"].dzSpread_pf <= 5,
        "vxySpread_pf <= 5": lambda objs: objs["mu_ljs"].vxySpread_pf <= 5,
        "dzSpread_pf <= 10": lambda objs: objs["mu_ljs"].dzSpread_pf <= 10,
        "vxySpread_pf <= 10": lambda objs: objs["mu_ljs"].vxySpread_pf <= 10,
        "dxySpread_pf <= 10": lambda objs: objs["mu_ljs"].dxySpread_pf <= 10,
        "dzSpread_dsa <= 500000": lambda objs: objs["mu_ljs"].dzSpread_dsa <= 500000,
        "vxySpread_dsa <= 500000": lambda objs: objs["mu_ljs"].vxySpread_dsa <= 500000,
        "dzSpread_mu <= 500000": lambda objs: objs["mu_ljs"].dzSpread_mu <= 500000,
        "vxySpread_mu <= 500000": lambda objs: objs["mu_ljs"].vxySpread_mu <= 500000,
        "dzSpread_mu <= 50": lambda objs: objs["mu_ljs"].dzSpread_mu <= 50,
        "vxySpread_mu <= 50": lambda objs: objs["mu_ljs"].vxySpread_mu <= 50,
        "dxySpread_mu <= 50": lambda objs: objs["mu_ljs"].dxySpread_mu <= 50,
        "dzSpread_mu > 50": lambda objs: objs["mu_ljs"].dzSpread_mu > 50,
        "vxySpread_mu > 50": lambda objs: objs["mu_ljs"].vxySpread_mu > 50,
        "dzSpread_pf > 10": lambda objs: objs["mu_ljs"].dzSpread_pf > 10,
        "vxySpread_pf > 10": lambda objs: objs["mu_ljs"].vxySpread_pf > 10,
        "vzSpread_pf <= 10": lambda objs: objs["mu_ljs"].vzSpread_pf <= 10,

        "inverse_spread" : lambda objs: ((objs["mu_ljs"].vzSpread_pf > 10)
                                         |(objs["mu_ljs"].vxySpread_pf > 10)
                                         |(objs["mu_ljs"].vxySpread_mu > 50)
                                         |(objs["mu_ljs"].vxySpread_mu > 50) )
    },
    "genMus":{
        "pT >= 10 GeV": lambda objs: objs["genMus"].pt >= 10,
        "status 1": lambda objs: objs["genMus"].status == 1,
        "status 23": lambda objs: objs["genMus"].status == 23,
        "fromGenA": lambda objs: abs(objs["genMus"].distinctParent.pdgId) == 32,
        "parent == A": lambda objs: abs(objs["genMus"].parent.pdgId) == 32,
    },
    "genEs":{
        "status 1": lambda objs: objs["genEs"].status == 1,
        "status 23": lambda objs: objs["genEs"].status == 23,
        "fromGenA": lambda objs: abs(objs["genEs"].distinctParent.pdgId) == 32,
        "parent == A": lambda objs: abs(objs["genEs"].parent.pdgId) == 32,
    },
    "genAs": {
        "dR(A, LJ) < 0.2": lambda objs: dR(objs["genAs"], objs["ljs"]) < 0.2,
        "dR(A, LJ) < 0.4": lambda objs: dR(objs["genAs"], objs["ljs"]) < 0.4,
        "lxy < 10 cm": lambda objs: lxy(objs["genAs"]) < 10,
        "lxy < 40 cm": lambda objs: lxy(objs["genAs"]) < 40,
        "10 cm <= lxy < 100 cm": lambda objs: ((lxy(objs["genAs"]) >= 10)
                                               & (lxy(objs["genAs"]) < 100)),
        "100 cm <= lxy < 135 cm": lambda objs: ((lxy(objs["genAs"]) >= 100)
                                               & (lxy(objs["genAs"]) < 135)),
        "lxy >= 100 cm": lambda objs: lxy(objs["genAs"]) >= 100,
        "lxy <= 100 cm": lambda objs: lxy(objs["genAs"]) <= 100,
        "lxy <= 150 cm": lambda objs: lxy(objs["genAs"]) <= 150,
        "lxy <= 250 cm": lambda objs: lxy(objs["genAs"]) <= 250,
        "lxy <= 400 cm": lambda objs: lxy(objs["genAs"]) <= 400,
        "pT > 30 GeV": lambda objs: objs["genAs"].pt > 30,
        "pT < 300 GeV": lambda objs: objs["genAs"].pt < 300,
    },
    "genAs_toMu": {
        "dR(A, LJ) < 0.2": lambda objs: dR(objs["genAs_toMu"], objs["ljs"]) < 0.2,
        "dR(A, LJ) < 0.4": lambda objs: dR(objs["genAs_toMu"], objs["ljs"]) < 0.4,
        "lxy < 10 cm": lambda objs: lxy(objs["genAs_toMu"]) < 10,
        "10 cm <= lxy < 100 cm": lambda objs: (lxy(objs["genAs_toMu"]) >= 10) & (lxy(objs["genAs_toMu"]) < 100),
        "lxy <= 60 cm": lambda objs: lxy(objs["genAs_toMu"]) <= 60,
        "lxy >= 100 cm": lambda objs: lxy(objs["genAs_toMu"]) >= 100,
        "lxy <= 100 cm": lambda objs: lxy(objs["genAs_toMu"]) <= 100,
        "lxy <= 150 cm": lambda objs: lxy(objs["genAs_toMu"]) <= 150,
        "lxy <= 250 cm": lambda objs: lxy(objs["genAs_toMu"]) <= 250,
        "lxy <= 400 cm": lambda objs: lxy(objs["genAs_toMu"]) <= 400,
        "pT > 30 GeV": lambda objs: objs["genAs_toMu"].pt > 30,
        "pT < 300 GeV": lambda objs: objs["genAs_toMu"].pt < 300,
    },
    "genAs_toE": {
        "dR(A, LJ) < 0.2": lambda objs: dR(objs["genAs_toE"], objs["ljs"]) < 0.2,
        "dR(A, LJ) < 0.4": lambda objs: dR(objs["genAs_toE"], objs["ljs"]) < 0.4,
        "lxy <= 5 cm": lambda objs: lxy(objs["genAs_toE"]) <= 5,
        "lxy <= 2.5 cm": lambda objs: lxy(objs["genAs_toE"]) <= 2.5,
        "lxy < 10 cm": lambda objs: lxy(objs["genAs_toE"]) < 10,
        "lxy < 40 cm": lambda objs: lxy(objs["genAs_toE"]) < 40,
        "10 cm <= lxy < 100 cm": lambda objs: (lxy(objs["genAs_toE"]) >= 10) & (lxy(objs["genAs_toE"]) < 100),
        "40 cm <= lxy < 77 cm": lambda objs: (lxy(objs["genAs_toE"]) >= 40) & (lxy(objs["genAs_toE"]) < 77),
        "100 cm <= lxy < 135 cm": lambda objs: (lxy(objs["genAs_toE"]) >= 100) & (lxy(objs["genAs_toE"]) < 135),
        "125 cm <= lxy <= 135 cm": lambda objs: (lxy(objs["genAs_toE"]) >= 125) & (lxy(objs["genAs_toE"]) < 135),
        "lxy >= 100 cm": lambda objs: lxy(objs["genAs_toE"]) >= 100,
        "lxy <= 100 cm": lambda objs: lxy(objs["genAs_toE"]) <= 100,
        "lxy <= 150 cm": lambda objs: lxy(objs["genAs_toE"]) <= 150,
        "lxy <= 250 cm": lambda objs: lxy(objs["genAs_toE"]) <= 250,
        "lxy <= 400 cm": lambda objs: lxy(objs["genAs_toE"]) <= 400,
        "pT > 30 GeV": lambda objs: objs["genAs_toE"].pt > 30,
        "pT < 300 GeV": lambda objs: objs["genAs_toE"].pt < 300,
    },
    "electrons": {
        "pT > 10 GeV": lambda objs: objs["electrons"].pt > 10,
        "pT > 30 GeV": lambda objs: objs["electrons"].pt > 30,
        "|eta| < 1.479": lambda objs: abs(objs["electrons"].eta) < 1.479,
        "eta-phi_veto" : lambda objs: ~ ((objs["electrons"].eta <= 1.4)
                                         & (objs["electrons"].eta >= -0.1)
                                         & (objs["electrons"].phi >= 0.4)
                                         & (objs["electrons"].phi <=0.8)),
        "1.479 < |eta| < 2.4": lambda objs: ((abs(objs["electrons"].eta) > 1.479)
                                             & (abs(objs["electrons"].eta) < 2.4)),
        "|eta| < 2.4": lambda objs: abs(objs["electrons"].eta) < 2.4,
        "dR(e, A) < 0.5": lambda objs: dR(objs["electrons"], objs["genAs_toE"]) < 0.5,
        "looseID": lambda objs: objs["electrons"].cutBased >= 2,
        "photonIdx": lambda objs: objs["electrons"].photonIdx != -1,
        "barrel": lambda objs: abs(objs["electrons"].eta + objs["electrons"].deltaEtaSC) <= 1.479 ,
        "endcap": lambda objs: ((abs(objs["electrons"].eta + objs["electrons"].deltaEtaSC) > 1.479) & (abs(objs["electrons"].eta + objs["electrons"].deltaEtaSC) < 2.4)),
        "barrel_sieie": lambda objs: objs["electrons"].sieie <= 0.0112,
        "barrel_iso": lambda objs: objs["electrons"].pfRelIso03_all < (.112 + .506/objs["electrons"].pt),
        "barrel_hoe": lambda objs: objs["electrons"].hoe < 0.05 + 1.16/((1 + objs["electrons"].scEtOverPt) * objs["electrons"].pt) + 0.0324*objs["rho_PFIso"]/((1 + objs["electrons"].scEtOverPt) * objs["electrons"].pt),
        "barrel_eInvMinusPInv": lambda objs: objs["electrons"].eInvMinusPInv < 0.193,
        "lostHits": lambda objs: objs["electrons"].lostHits <= 1,
        "convVeto": lambda objs: objs["electrons"].convVeto == True,
        "endcap_sieie": lambda objs: objs["electrons"].sieie <= 0.0425,
        "endcap_iso": lambda objs: objs["electrons"].pfRelIso03_all < (.108 + .963/objs["electrons"].pt),
        "endcap_hoe": lambda objs: objs["electrons"].hoe <  0.0441 + 2.54/((1 + objs["electrons"].scEtOverPt) * objs["electrons"].pt) + 0.183*objs["rho_PFIso"]/((1 + objs["electrons"].scEtOverPt) * objs["electrons"].pt),
        "endcap_eInvMinusPInv": lambda objs: objs["electrons"].eInvMinusPInv < 0.111,
        'MVANonIsoWPL': lambda objs: objs['electrons'].mvaFall17V2noIso_WPL,
        "missing_hits == 0" : lambda objs: objs["electrons"].lostHits == 0,
        "missing_hits == 1" : lambda objs: objs["electrons"].lostHits == 1,
       "missing_hits >= 1" : lambda objs: objs["electrons"].lostHits >= 1,
        "dxy >= 0.05" : lambda objs: objs["electrons"].dxy >= 0.05,
    },
    "jets":{
        "btagDeepFlavB_tight": lambda objs: objs["jets"].btagDeepFlavB >= 0.7100,
        "tightID": lambda objs: objs["jets"].jetId >= 6,
        "|eta| < 2.4": lambda objs: abs(objs["jets"].eta) < 2.4,
        "pT > 30 GeV": lambda objs: objs["jets"].pt > 30,
        "dR(jets, mu) > 0.4": lambda objs: dR(objs["jets"], objs["muons"]) > 0.4,
    },
    "bjets":{
        "tightID": lambda objs: objs["bjets"].jetId >= 6,
        "|eta| < 2.4": lambda objs: abs(objs["bjets"].eta) < 2.4,
        "pT > 30 GeV": lambda objs: objs["bjets"].pt > 30,
        "dR(bjets, mu) > 0.4": lambda objs: dR(objs["bjets"], objs["muons"]) > 0.4,
    },
    "muons": {
        #Tested the following to try to enable us to apply these cuts to muons *and* matched_muons associated to dsas
        "looseID": lambda objs, muons: muons.looseId,
        "tightID": lambda objs, muons: muons.tightId,
        "pT > 5 GeV": lambda objs, muons: muons.pt > 5,
        "pT > 30 GeV": lambda objs, muons: muons.pt > 30,
        "eta-phi_veto" : lambda objs, muons: ~ ((muons.eta <= 1.4)
                                              & (muons.eta >= -0.1)
                                              & (muons.phi >= 0.4)
                                              & (muons.phi <=0.8)),
        "|eta| < 2.4": lambda objs, muons: abs(muons.eta) < 2.4,
         "|eta| < 2.1": lambda objs, muons: abs(muons.eta) < 2.1,
        "dR(mu, A) < 0.5": lambda objs, muons: dR(muons, objs["genAs_toMu"]) < 0.5,
        "dR(mu, A) < 0.5 nested": lambda objs, muons: dR(muons, objs["genAs_toMu"][:,:,None]) < 0.5,
        "dxy >= 0.008" : lambda objs, muons: muons.dxy >= 0.008,
        "trkNumPixelHits <= 2" : lambda objs, muons: muons.trkNumPixelHits <= 2,
        "barrel": lambda objs, muons: abs(muons.eta) <= 1.479,
        "endcap": lambda objs, muons: ((abs(muons.eta) > 1.479) & (abs(muons.eta) < 2.4)),
        "pfRelIso04_all <= 0.15": lambda objs, muons: muons.pfRelIso04_all <= 0.15,
        "sip3d < 4": lambda objs, muons: abs(muons.sip3d) < 4,
    },
    "photons":{
        "pT > 20 GeV": lambda objs: objs["photons"].pt > 20,
        "pT > 30 GeV": lambda objs: objs["photons"].pt > 30,
        "eta-phi_veto" : lambda objs: ~ ((objs["photons"].eta <= 1.4)
                                         & (objs["photons"].eta >= -0.1)
                                         & (objs["photons"].phi >= 0.4)
                                         & (objs["photons"].phi <=0.8)),
        "|eta| < 2.5": lambda objs: abs(objs["photons"].eta) < 2.5, # fixme: do we want eta or scEta
        "eta": lambda objs: objs["photons"].isScEtaEB | objs["photons"].isScEtaEE,
        "barrel": lambda objs: objs["photons"].isScEtaEB,
        "endcap": lambda objs: objs["photons"].isScEtaEE,
        "looseID": lambda objs: objs["photons"].cutBased >= 1,
        "pixelSeed": lambda objs: objs["photons"].pixelSeed == False ,
        "electronVeto": lambda objs: objs["photons"].electronVeto,
        "dR(gm, A) < 0.5": lambda objs: dR(objs["photons"], objs["genAs_toE"]) < 0.5,
        "barrel_hoe": lambda objs: objs["photons"].hoe <= .04596,
        "barrel_sieie": lambda objs: objs["photons"].sieie <=  0.0106 ,
        "|eta| < 1.479": lambda objs: abs(objs["photons"].eta) < 1.479,
        'Custom Cutbased': lambda objs: returnBitMapTArrayPhoton(objs['photons'].vidNestedWPBitmap, 'PhoFull5x5SigmaIEtaIEtaCut', 'PhoIsoWithEALinScalingCut'),
        'Photon DR Veto 0p025': lambda objs: (dR(objs['photons'],objs['electrons']) > 0.025),
        'Photon DR Veto 0p035': lambda objs: (dR(objs['photons'],objs['electrons']) > 0.035),
        'Photon DR Veto 0p045': lambda objs: (dR(objs['photons'],objs['electrons']) > 0.045),
    },
    "dsaMuons": {
        "pT > 10 GeV": lambda objs, dsa: dsa.pt > 10,
        "pT > 100 GeV": lambda objs, dsa: dsa.pt > 100,
        "eta-phi_veto" : lambda objs, dsa: ~ ((dsa.eta <= 1.4)
                                              & (dsa.eta >= -0.1)
                                              & (dsa.phi >= 0.4)
                                              & (dsa.phi <=0.8)),
        "|dxy| <= 40": lambda objs, dsa: abs(dsa.dxy) <= 40,
        "|dz| <= 60": lambda objs, dsa: abs(dsa.dz) <= 60,
        "|eta| < 2.4": lambda objs, dsa: abs(dsa.eta) < 2.4,
        "dR(dsa, A) < 0.4": lambda objs, dsa: dR(dsa, objs["genAs_toMu"]) < 0.4,
        # displaced ID as a single flag and as individual cuts
        "displaced ID" : lambda objs, dsa: dsa.displacedID > 0,
        "DT + CSC hits > 12": lambda objs, dsa: (dsa.trkNumDTHits
                                            + dsa.trkNumCSCHits) > 12,
        "ifcsczero": lambda objs, dsa: ak.where(((dsa.trkNumCSCHits == 0)
                                           & (dsa.trkNumDTHits <= 18)), False, True),
        "normChi2 < 2.5": lambda objs, dsa: dsa.normChi2 < 2.5,
        "ptErrorOverPT < 1": lambda objs, dsa: (dsa.ptErr / dsa.pt) < 1.0,
        "barrel": lambda objs, dsa: abs(dsa.eta) <= 1.479,
        "endcap": lambda objs, dsa: ((abs(dsa.eta) > 1.479) & (abs(dsa.eta) < 2.4)),
        # just use segment-based matching
       # "no PF match" : lambda objs, dsa: dsa.muonMatch1/dsa.nSegments < 0.667,
        "dR(mu, A) < 0.5": lambda objs, dsa: dR(dsa, objs["genAs_toMu"]) < 0.5,
        "dR(mu, A) < 0.5 nested": lambda objs, dsa: dR(dsa, objs["genAs_toMu"][:,:,None]) < 0.5,
        "dR(dsa, pf) > 0.2": lambda objs, dsa: dR(dsa, objs["muons"]) > 0.2,
        "segmatch veto": lambda objs, dsa: ak.all(dsa.good_matched_muons.numMatch < 1, axis=2),
        "segmatch veto + charge": lambda objs, dsa: ak.all((dsa.good_matched_muons.numMatch < 1) | (dsa.charge[:,:,None] != dsa.good_matched_muons.charge), axis=2),
        "segmatch veto + dR outer": lambda objs, dsa: ak.all((dsa.good_matched_muons.numMatch < 1) | (dR_outer(dsa[:,:,None], dsa.good_matched_muons) > 0.1), axis=2),
        "segmatch veto + num seg": lambda objs, dsa: ak.all((dsa.good_matched_muons.numMatch < 1) | (dsa.good_matched_muons.numMatch/(dsa.nSegments[:,:,None]) < 0.34), axis=2),
        "old_all": lambda objs, dsa: ak.all((dsa.good_matched_muons.numMatch < 1) | (dR_outer(dsa[:,:,None], dsa.good_matched_muons) > 0.1) | (dsa.good_matched_muons.numMatch/(dsa.nSegments[:,:,None]) < 0.34), axis=2),
        "all": lambda objs, dsa: ak.all(
            ((dsa.good_matched_muons.numMatch < 1) | (dR_outer(dsa[:,:,None], dsa.good_matched_muons) > 0.1) | (dsa.good_matched_muons.numMatch/(dsa.nSegments[:,:,None]) < 0.34))
            & ((dsa.good_matched_muons.numMatch < 1) | (dsa.good_matched_muons.numMatch/(dsa.nSegments[:,:,None]) < 1)),
            axis=2,
        ),
        "all + charge": lambda objs, dsa: ak.all((dsa.good_matched_muons.numMatch < 1) | (dsa.charge[:,:,None] != dsa.good_matched_muons.charge) | (dR_outer(dsa[:,:,None], dsa.good_matched_muons) > 0.1) | (dsa.good_matched_muons.numMatch/(dsa.nSegments[:,:,None]) < 0.34), axis=2),
    },
    "dsaMuonPairs":{
        "back_to_back_pairs": lambda objs: (lambda v1, v2: np.cos(v1.deltaangle(v2)))
                                         (*ak.unzip(objs["dsaMuonPairs"])) <= -0.95,
        "parallel_pairs": lambda objs: (lambda v1, v2: np.cos(v1.deltaangle(v2)))
                                         (*ak.unzip(objs["dsaMuonPairs"])) >= 0.95,
        "lj_index_diff = 0": lambda objs: (lambda v1, v2, objs: abs(nearest_lj_index(v1, objs)  - nearest_lj_index(v2, objs)))(*ak.unzip(objs["dsaMuonPairs"]), objs["ljs"]) == 0,
}
}

def _gen_filter_count(leps, flag_bits, status=None):
    # Count leptons passing the RunIIFall18 production gen filter: pT > 5 GeV,
    # |eta| < 2.4, production vertex rho < 740 cm and |z| < 960 cm, plus the
    # statusFlags bits required by the filter's cut string
    # (0=isPrompt, 7=isHardProcess, 8=fromHardProcess, 13=isLastCopy).
    if status is not None:
        leps = leps[leps.status == status]
    passing = (
        (leps.pt > 5.0)
        & (abs(leps.eta) < 2.4)
        & (rho(leps, use_v=True) < 740.0)
        & (abs(leps.vz) < 960.0)
        & check_bits(leps.statusFlags, flag_bits)
    )
    return ak.num(leps[passing])

evt_cut_defs = {
    # This following will be True for every event. There's probably a more intuitive way to do this
    "Keep all evts": lambda objs: objs["pvs"].npvs >= 0,
    # Re-application of the central-production gen filter. The v10 signal samples
    # were PRODUCED with this filter (their GenFilterInfo records the tried/passed
    # counts), so on them the cutflow of these cuts is a re-pass/validation rate
    # (0.84-1.00, mass-dependent) -- NOT a filter efficiency. Per-point production
    # efficiencies are recorded in
    # sidm/studies/truth_kinematics_forAN/central_genFilterEfficiencies.yml.
    # "pass gen filter" follows the McM request's cut string (isLastCopy &&
    # isPromptFinalState && fromHardProcessFinalState); the isHardProcess variant
    # follows the piedpiper template. Both count electrons and muons together, >= 4.
    "pass gen filter": lambda objs: (
        _gen_filter_count(objs["genMus"], [0, 8, 13], status=1)
        + _gen_filter_count(objs["genEs"], [0, 8, 13], status=1)
    ) >= 4,
    "pass gen filter isHardProcess": lambda objs: (
        _gen_filter_count(objs["genMus"], [7])
        + _gen_filter_count(objs["genEs"], [7])
    ) >= 4,
    "pass triggers": lambda objs: (
          objs["hlt"].DoubleL2Mu23NoVtx_2Cha
        | objs["hlt"].DoubleL2Mu23NoVtx_2Cha_CosmicSeed
        | objs["hlt"].DoubleL2Mu25NoVtx_2Cha_Eta2p4
        | objs["hlt"].DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4
    ),
    "pass flags": lambda objs: (
          objs["flags"].goodVertices
        & objs["flags"].globalSuperTightHalo2016Filter
        & objs["flags"].HBHENoiseFilter
        & objs["flags"].HBHENoiseIsoFilter
        & objs["flags"].EcalDeadCellTriggerPrimitiveFilter
        & objs["flags"].BadPFMuonFilter
        & objs["flags"].BadPFMuonDzFilter
        & objs["flags"].eeBadScFilter
        & objs["flags"].ecalBadCalibFilter
        & objs["flags"].hfNoisyHitsFilter
    ),
    ">=1 muon": lambda objs: ak.num(objs["muons"]) >= 1,
    ">=2 muons (pf or dsa) pT > 26 GeV": lambda objs: (
        ak.sum(objs["muons"].pt > 26, axis=-1) + ak.sum(objs["dsaMuons"].pt > 26, axis=-1) >= 2
    ),

    "PV filter": lambda objs: ak.flatten(objs["pvs"].npvsGood) >= 1,
    #"Cosmic veto": lambda objs: objs["cosmicveto"].result,
    ">=2 LJs": lambda objs: ak.num(objs["ljs"]) >= 2,
    ">=1 LJs": lambda objs: ak.num(objs["ljs"]) >= 1,
    "< 2 LJs": lambda objs: ak.num(objs["ljs"]) < 2,
    ">=1 egm_ljs": lambda objs: ak.num(objs["egm_ljs"]) >= 1,
    ">=1 mu_ljs": lambda objs: ak.num(objs["mu_ljs"]) >= 1,
    "1 egm_ljs": lambda objs: ak.num(objs["egm_ljs"]) == 1,
    "1 mu_ljs": lambda objs: ak.num(objs["mu_ljs"]) == 1,
    "= 2 LJs": lambda objs: ak.num(objs["ljs"]) == 2,
    "= 1 LJs": lambda objs: ak.num(objs["ljs"]) == 1,
    ">=2 matched As": lambda objs: ak.num(derived_objs["genAs_matched_lj"](objs, 0.2)) >= 2,
    # 4mu: leading two LJs are both mu-type
    "4mu": lambda objs: ak.count_nonzero(objs["ljs"][:, :2].muon_n >= 2, axis=-1) == 2,
    # 2mu2e: leading two LJs contain exactly 1 mu-type and exactly 1 egm-type
    "2mu2e": lambda objs: ((ak.count_nonzero(objs["ljs"][:, :2].muon_n >= 2, axis=-1) == 1)
                           & (ak.count_nonzero(objs["ljs"][:, :2].muon_n == 0, axis=-1) == 1)),
    "genAs_toE_matched_egmLj": lambda objs: ak.num(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4)) >= 1,
    "genAs_toMu_matched_muLj": lambda objs: ak.num(derived_objs["genAs_toMu_matched_muLj"](objs, 0.4)) >= 1,
    "genAs_toE": lambda objs: ak.num(objs["genAs_toE"]) >= 1,
    "genAs_toMu": lambda objs: ak.num(objs["genAs_toMu"]) >= 1,
    "ljs": lambda objs: ak.num(objs["ljs"]) >= 1,
    "50 GeV <= GenMu0_pT <= 60 GeV": lambda objs : (objs["genMus"][:, 0].pt >=50) & (objs["genMus"][:, 0].pt <=60),
    "genMus": lambda objs: ak.num(objs["genMus"]) > 1,
    "dR(Mu_0, Mu_1) > 0.03": lambda objs: objs["genMus"][:,0].delta_r(objs["genMus"][:,1]) > 0.03,
    "LJ-LJ dPhi > 2": lambda objs: abs(objs["ljs"][:, 0].delta_phi(objs["ljs"][:, 1])) > 2.0,
    "= 1 LJs": lambda objs: ak.num(objs["ljs"]) == 1,
    "= 1 muLj or 1 egmLJ": lambda objs:((ak.num(objs["mu_ljs"]) == 1) & (ak.num(objs["egm_ljs"]) == 0)) | ((ak.num(objs["mu_ljs"]) == 0) & (ak.num(objs["egm_ljs"]) == 1)),
    "= 1 muLJs": lambda objs: (ak.num(objs["mu_ljs"]) == 1) & (ak.num(objs["egm_ljs"]) == 0),
    "= 1 egmLJs": lambda objs: (ak.num(objs["mu_ljs"]) == 0) & (ak.num(objs["egm_ljs"]) == 1),
    "= 2 LJs": lambda objs: ak.num(objs["ljs"]) == 2,
    "= 1 muLJs and = 1 egmLJs": lambda objs: (ak.num(objs["mu_ljs"]) == 1) & (ak.num(objs["egm_ljs"]) == 1),
    "= 2 muLJs": lambda objs: (ak.num(objs["mu_ljs"]) == 2) & (ak.num(objs["egm_ljs"]) == 0),
    "60 <= inv(Mu_0, Mu_1) <= 120": lambda objs : ((objs["muons"][:,:2].sum().mass) <= 120) &  ((objs["muons"][:,:2].sum().mass) >= 60),
    "n_mu == 2": lambda objs : ak.num(objs["muons"]) == 2,
    "n_dsa >= 2": lambda objs : ak.num(objs["dsaMuons"]) >= 2,
    "n_bjet == 2": lambda objs : ak.num(objs["bjets"]) == 2,
    "n_jet == 2": lambda objs : ak.num(objs["jets"]) == 2,
    "2j2btag": lambda objs : (ak.num(objs["jets"]) == 2) &  (ak.num(objs["bjets"]) == 2),
    "4j2btag": lambda objs : (ak.num(objs["jets"]) >= 4) &  (ak.num(objs["bjets"]) == 2),
    "dPhi(Mu_0, Mu_1) > 2.4": lambda objs: ak.fill_none(abs((ak.pad_none(objs["muons"], 2)[:,0].delta_phi(ak.pad_none(objs["muons"], 2)[:,1]))) > 2.4, False),
    "pv_ndof >=4" : lambda objs :  ak.flatten(objs["pvs"].ndof) >=4.0,
    "pv_z <= 24" : lambda objs :  ak.flatten(objs["pvs"].z) <= 24,
    "inv(lj1, lj2) <= 100": lambda objs : (objs["ljs"][:,:2].sum().mass) <= 100,
    "all cos_alpha(dsa, dsa) > -0.95" : lambda objs : ak.all(cosAlpha(objs["dsaMuons"]) > -0.95, axis =1),
    "all cos_alpha(mu, mu) > -0.95" : lambda objs : ak.all(cosAlpha(objs["muons"]) > -0.95, axis =1),
    "any cos_alpha(dsa, dsa) <= -0.95" : lambda objs : ak.any(cosAlpha(objs["dsaMuons"]) <=-0.95, axis =1),
    "1 back_to_back_dsa": lambda objs: (lambda v1, v2: ak.sum(np.cos(v1.deltaangle(v2)) <= -0.95, axis=1))(*ak.unzip(objs["dsaMuonPairs"]))  == 1,
    "2 back_to_back_dsa": lambda objs: (lambda v1, v2: ak.sum(np.cos(v1.deltaangle(v2)) <= -0.95, axis=1))(*ak.unzip(objs["dsaMuonPairs"]))  == 2,
    "pass two missing triggers": lambda objs: (
    (
        objs["hlt"].DoubleL2Mu23NoVtx_2Cha_NoL2Matched
      | objs["hlt"].DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched
    )
    & ~(objs["hlt"].DoubleL2Mu23NoVtx_2Cha)
    & ~(objs["hlt"].DoubleL2Mu23NoVtx_2Cha_CosmicSeed)
    & ~(objs["hlt"].DoubleL2Mu25NoVtx_2Cha_Eta2p4)
    & ~(objs["hlt"].DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4)
),
    "all cosAlpha(mu, mu) > -0.99"  : lambda objs: ak.all(cosAlpha(objs["allMuons"]) > -0.99, axis=1),
    "any cosAlpha(mu, mu) <= -0.99" : lambda objs: ak.any(cosAlpha(objs["allMuons"]) <= -0.99, axis=1),
}
