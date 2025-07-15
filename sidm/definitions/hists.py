"""Define all available histograms

All hists are defined as Histogram objects whose axes are given as a list of Axis objects, which
bundle a hist.axis with a function that defines how the axis will be filled. The underlying
hist.Hists storage is weight unless otherwise specified.
"""

# python
import math
import importlib
# columnar analysis
import hist
import awkward as ak
# local
from sidm.tools import histogram as h
from sidm.tools.utilities import dR, lxy, matched, dxy
from sidm.definitions.objects import derived_objs
# always reload local modules to pick up changes during development
importlib.reload(h)


# define counters
counter_defs = {
    "Total LJs": lambda objs: ak.count(objs["ljs"].pt),
    "Gen As to muons": lambda objs: ak.count(objs["genAs_toMu"].pt),
    "Gen As to electrons": lambda objs: ak.count(objs["genAs_toE"].pt),
    "Matched gen As to muons": lambda objs: ak.count(derived_objs["genAs_toMu_matched_lj"](objs, 0.4).pt),
    "Matched gen As to electrons": lambda objs: ak.count(derived_objs["genAs_toE_matched_lj"](objs, 0.4).pt),
}
    
# define default labels and binnings
obj_labels = {
    "electrons": "Electron",
    "photons": "Photon",
    "muons": "PF Muon",
    "dsaMuons": "DSA Muon",
    "ljs": "Lepton Jet",
    "mu_ljs": r"$\mu$-type Lepton Jet",
    "egm_ljs": r"$e\gamma$-type Lepton Jet",
    "genAs": r"$Z_d$",
    "genAs_toMu": r"$Z_d\rightarrow \mu\mu$",
    "genAs_toE": r"$Z_d\rightarrow ee$",
    "pvs": "PV",
}
attr_labels = {
    "pt": r"$p_T$ (GeV)",
    "eta": r"$\eta$",
    "phi": r"$\phi$",
    "lxy": r"$L_{{xy}}$ (cm) ",
    "dxy": r"$d_0$",
}
default_binnings = {
    "n":  (10, 0, 10),
    "pt":  (100, 0, 100),
    "eta": (50, -3, 3),
    "phi": (50, -1*math.pi, math.pi),
    "lxy": (100, 0, 100),
}

# define convenience functions to simplify creating basic hists
def make_label(obj, attr, absval):
    obj_label = obj_labels.get(obj, obj)
    if attr == "n":
        return f"Number of {obj_label}s"
    attr = attr_labels.get(attr, attr)
    if absval:
        attr = f"|{attr}|"
    return f"{obj_label} {attr}"

def obj_attr(obj, attr, absval=False, nbins=None, xmin=None, xmax=None, label=None):
    (_nbins, _xmin, _xmax) = default_binnings.get(attr, (100, 0, 100))
    nbins = _nbins if nbins is None else nbins
    xmin = _xmin if xmin is None else xmin
    xmax = _xmax if xmax is None else xmax
    label = make_label(obj, attr, absval) if label is None else label
    return h.Histogram.simple_hist(obj, attr, absval, nbins, xmin, xmax, label)

def make_2d(h1, h2):
    return h.Histogram([h1.axes[-1], h2.axes[-1]])

def obj_eta_phi(obj, nbins_x=None, xmin=None, xmax=None, nbins_y=None, ymin=None, ymax=None):
    return make_2d(
        obj_attr(obj, "eta", nbins_x, xmin, xmax),
        obj_attr(obj, "phi", nbins_y, ymin, ymax),
    )


# define histograms
hist_defs = {
    # pv
    "pv_n": obj_attr("pvs", "npvs", nbins=50, label="Number of PVs"),
    "pv_ndof": obj_attr("pvs", "ndof", nbins=25, xmax=100),
    "pv_z": obj_attr("pvs", "z", xmin=-50, xmax=50),
    "pv_rho": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, -0.5, 0.5, name="pv_rho"),
                   lambda objs, mask: objs["pvs"].pos.rho),
        ],
    ),
    # GSFelectron: Plottting electron ID varaiables and plotting 2D hists of the leading electron
    # ID variables in barrel within Delta R < .5 of a dark photon vs the lxy of the dark photon
    "electron_GsfEleDEtaInSeedCut": h.Histogram(
        [
            h.Axis(hist.axis.Regular(35, 0, .0070, name="electron_GsfEleDEtaInSeedCut"),
                   lambda objs, mask: objs["electrons"].GsfEleDEtaInSeedCut_0),
        ],
    ),
    "electron_GsfEleDEtaInSeedCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]),
            h.Axis(hist.axis.Regular(35, 0, .0070, name="electron_GsfEleDEtaInSeedCut"),
                   lambda objs, mask: matched(objs["electrons"], objs["genAs_toE"], 0.5)[mask, 0:1].GsfEleDEtaInSeedCut_0)
        ],
        evt_mask = lambda objs: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5)) > 0,
    ),
    "electron_GsfEleDPhiInCut": h.Histogram(
        [
            h.Axis(hist.axis.Regular(45, 0, .0450, name="electron_GsfEleDPhiInCut"),
                   lambda objs, mask: objs["electrons"].GsfEleDPhiInCut_0),
        ],
    ),
    "electron_GsfEleDPhiInCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]),
            h.Axis(hist.axis.Regular(45, 0, .09, name="electron_GsfEleDPhiInCut"),
                   lambda objs, mask: matched(objs["electrons"], objs["genAs_toE"], 0.5)[mask, 0:1].GsfEleDPhiInCut_0)
        ],
        evt_mask = lambda objs: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5)) > 0,
    ),
    "electron_GsfEleEInverseMinusPInverseCut": h.Histogram(
        [
            h.Axis(hist.axis.Regular(60, 0, .3, name="electron_GsfEleEInverseMinusPInverseCut"),
                   lambda objs, mask: objs["electrons"].GsfEleEInverseMinusPInverseCut_0),
        ],
    ),
    "electron_GsfEleEInverseMinusPInverseCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]),
            h.Axis(hist.axis.Regular(60, 0, .3, name="electron_GsfEleEInverseMinusPInverseCut"),
                   lambda objs, mask: matched(objs["electrons"], objs["genAs_toE"], 0.5)[mask, 0:1].GsfEleEInverseMinusPInverseCut_0)
        ],
        evt_mask = lambda objs: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5)) > 0,
    ),
    "electron_GsfEleRelPFIsoScaledCut": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, .2, name="electron_GsfEleRelPFIsoScaledCut"),
                   lambda objs, mask: (objs["electrons"].GsfEleRelPFIsoScaledCut_0
                                       - .506/objs["electrons"].pt)),
        ],
    ),
    "electron_GsfEleRelPFIsoScaledCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
           # added the alegbra relIso has in the analysis note
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]),
            h.Axis(hist.axis.Regular(40, 0, .2, name="electron_GsfEleRelPFIsoScaledCut"),
                   lambda objs, mask: (matched(objs["electrons"], objs["genAs_toE"], 0.5)[mask, 0:1].GsfEleRelPFIsoScaledCut_0
                                       - .506/(matched(objs["electrons"], objs["genAs_toE"], 0.5)[mask, 0:1].pt))),
        ],
        evt_mask = lambda objs: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5)) > 0,
    ),
    "electron_GsfEleFull5x5SigmaIEtaIEtaCut": h.Histogram(
        [
            h.Axis(hist.axis.Regular(45, 0, .045, name="electron_GsfEleFull5x5SigmaIEtaIEtaCut"),
                   lambda objs, mask: objs["electrons"].GsfEleFull5x5SigmaIEtaIEtaCut_0),
        ],
    ),
    "electron_GsfEleFull5x5SigmaIEtaIEtaCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]),
            h.Axis(hist.axis.Regular(45, 0, .045, name="electron_GsfEleFull5x5SigmaIEtaIEtaCut"),
                   lambda objs, mask: matched(objs["electrons"], objs["genAs_toE"], 0.5)[mask, 0:1].GsfEleFull5x5SigmaIEtaIEtaCut_0)
        ],
        evt_mask = lambda objs: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5)) > 0,
    ),
    "electron_GsfEleConversionVetoCut": h.Histogram(
        [
            h.Axis(hist.axis.Regular(2, 0, 2, name="electron_GsfEleConversionVetoCut"),
                   lambda objs, mask: objs["electrons"].GsfEleConversionVetoCut_0),
        ],
    ),
    "electron_GsfEleConversionVetoCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]),
            h.Axis(hist.axis.Regular(2, 0, 2, name="electron_GsfEleConversionVetoCut"),
                   lambda objs, mask: matched(objs["electrons"], objs["genAs_toE"], 0.5)[mask, 0:1].GsfEleConversionVetoCut_0)
        ],
        evt_mask = lambda objs: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5)) > 0,
    ),
    "electron_GsfEleHadronicOverEMEnergyScaledCut": h.Histogram(
         [
             h.Axis(hist.axis.Regular(30, 0, .15, name="electron_GsfEleHadronicOverEMEnergyScaledCut"),
                    lambda objs, mask: objs["electrons"].GsfEleHadronicOverEMEnergyScaledCut_0),
         ],
     ),
    "electron_GsfEleHadronicOverEMEnergyScaledCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]),
            h.Axis(hist.axis.Regular(30, 0, .15, name="electron_GsfEleHadronicOverEMEnergyScaledCut"),
                   lambda objs, mask: matched(objs["electrons"], objs["genAs_toE"], 0.5)[mask, 0:1].GsfEleHadronicOverEMEnergyScaledCut_0)
        ],
        evt_mask = lambda objs: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5)) > 0,
    ),
    "electron_GsfEleMissingHitsCut": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="electron_GsfEleMissingHitsCut"),
                   lambda objs, mask: objs["electrons"].GsfEleMissingHitsCut_0),
        ],
    ),
    "electron_GsfEleMissingHitsCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]),
            h.Axis(hist.axis.Regular(10, 0, 10, name="electron_GsfEleMissingHitsCut"),
                   lambda objs, mask: matched(objs["electrons"], objs["genAs_toE"], 0.5)[mask, 0:1].GsfEleMissingHitsCut_0)
        ],
        evt_mask = lambda objs: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5)) > 0,
    ),
    # pfelectron
    "electron_n": obj_attr("electrons", "n", nbins=10),
    "electron_pt": obj_attr("electrons", "pt", xmax=500),
    "electron_dxy": obj_attr("electrons", "dxy",),
    "electron_eta_phi": obj_eta_phi("electrons"),
    "electron_photonIdx": obj_attr("electrons", "photonIdx", xmin=-1, xmax=10, nbins=10),
    "electron_dxy": obj_attr("electrons", "dxy", xmax=0.2),
    "electron_pfRelIso03_all": obj_attr("electrons", "pfRelIso03_all"),
    "electron_pfRelIso03_all_lowRange": obj_attr("electrons", "pfRelIso03_all", xmax=5),
    "electron_r9": obj_attr("electrons", "r9"),
    "electron_scEtOverPt": obj_attr("electrons", "scEtOverPt"),
    "electron_sieie": obj_attr("electrons", "sieie", xmax=.05),
    "electron_hoe": obj_attr("electrons", "hoe", xmax=1),
    "electron_eInvMinusPInv": obj_attr("electrons", "eInvMinusPInv", xmax=0.5),
    "electron_lostHits": obj_attr("electrons", "lostHits", xmax=10, nbins=10),
    "electron_deltaEtaSC": obj_attr("electrons", "deltaEtaSC", xmin=-0.1, xmax=0.1),
    "electron_nearGenA_n": h.Histogram(
        [
            # number of electrons within dR=0.5 of a genA that decays to electrons
            h.Axis(hist.axis.Integer(0, 10, name="electron_nearGenA_n"),
                   lambda objs, mask: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5))),
        ],
    ),

    "electron_nearGenE_n": h.Histogram(
        [
            # number of electrons within dR=0.5 of a genA that decays to electrons
            h.Axis(hist.axis.Integer(0, 10, name="electron_nearGenE_n"),
                   lambda objs, mask: ak.num(matched(objs["electrons"], objs["genEs"], 0.5))),
        ],
    ),


    # pfelectron-genA
    "electron_nearGenA_n_genA_lxy": h.Histogram(
        [
            # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])),
            # number of electrons within dR=0.5 of a genA that decays to electrons
            h.Axis(hist.axis.Integer(0, 4, name="electron_nearGenA_n", label="$N_{e}$"),
                   lambda objs, mask: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5))),
        ],
    ),
    "electron_genA_dR": h.Histogram(
        [
            # dR(e, nearest gen A)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="electron_genA_dR"),
                   lambda objs, mask: dR(objs["electrons"], objs["genAs"]))
        ],
    ),
    # pfelectron-genElectron
    "electron_genE_dR": h.Histogram(
        [
            # dR(e, nearest gen e)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="electron_genE_dR"),
                   lambda objs, mask: dR(objs["electrons"], objs["genEs"]))
        ],
    ),
    #pfelectron-pfphoton
    "electron_photon_dR": h.Histogram(
        [
            # dR(e, nearest gen e)
            h.Axis(hist.axis.Regular(50, 0, 0.1, name="electron_photon_dR"),
                   lambda objs, mask: dR(objs["electrons"], objs["photons"]))
        ],
    ),
    # pfphoton
    "photon_n": obj_attr("photons", "n"),
    "photon_pt":obj_attr("photons", "pt", xmax=500),
    "photon_hoe": obj_attr("photons", "hoe", xmax=1),
    "photon_sieie": obj_attr("photons", "sieie", xmax=.05),
    "photon_pfRelIso03_all": obj_attr("photons", "pfRelIso03_all", xmax=5),
    "photon_pfRelIso03_chg": obj_attr("photons", "pfRelIso03_chg", xmax=5),
    "photon_eta_phi": obj_eta_phi("photons"),
    "photon_nearGenA_n": h.Histogram(
        [
            # number of photons within dR=0.5 of a genA that decays to electrons
            h.Axis(hist.axis.Integer(0, 10, name="photon_nearGenA_n"),
                   lambda objs, mask: ak.num(matched(objs["photons"], objs["genAs_toE"], 0.5))),
        ],
    ),
    #electron-photon
    "electron_photon_dR": h.Histogram(
        [
            # dR(e, nearest gen e)
            h.Axis(hist.axis.Regular(50, 0, .3, name="electron_photon_dR"),
                   lambda objs, mask: dR(objs["electrons"], objs["photons"]))
        ],
    ),
    # pfphoton-genA
    "photon_nearGenA_n_genA_lxy": h.Histogram(
        [
            # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])),
            # number of photons within dR=0.5 of a genA that decays to electrons
            h.Axis(hist.axis.Integer(0, 4, name="photon_nearGenA_n", label=r"$N_{\gamma}$"),
                   lambda objs, mask: ak.num(matched(objs["photons"], objs["genAs_toE"], 0.5))),
        ],
    ),
    "photon_genA_dR": h.Histogram(
        [
            # dR(photon, nearest gen A)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="photon_genA_dR"),
                   lambda objs, mask: dR(objs["photons"], objs["genAs"]))
        ],
    ),
    # pfphoton-genElectron
    "photon_genE_dR": h.Histogram(
        [
            # dR(photon, nearest gen e)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="photon_genE_dR"),
                   lambda objs, mask: dR(objs["photons"], objs["genEs"]))
        ],
    ),
    # pfmuon
    "muon_n": obj_attr("muons", "n"),
    "muon_pt":obj_attr("muons", "pt", xmax=500),
    "muon_dxy":obj_attr("muons", "dxy"),
    "muon_eta_phi": obj_eta_phi("muons"),
    "muon_absD0": obj_attr("muons", "dxy", absval=True, xmax=500),
    "muon_absD0_lowRange": obj_attr("muons", "dxy", absval=True, xmax=10),
    "lj_muonChg_nearGenA": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -4, 4, name="lj_muonChg_nearGenA",
                                     label=r"LJ Muon Chg Near GenA"),
                   lambda objs, mask: (matched(objs["mu_ljs"], objs["genAs_toMu"], 0.5).muon_chg)),
        ],
    ),
    "muon_nearGenA_n": h.Histogram(
        [
            # number of muons within dR=0.5 of a genA that decays to muons
            h.Axis(hist.axis.Integer(0, 10, name="muon_nearGenA_n"),
                   lambda objs, mask: ak.num(matched(objs["muons"], objs["genAs_toMu"], 0.5))),
        ],
    ),
    # pfmuon-genA
    "muon_nearGenA_n_genA_lxy": h.Histogram(
        [
            # lxy of dark photon that decays to muons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toMu"])),
            # number of muons within dR=0.5 of a genA that decays to muons
            h.Axis(hist.axis.Integer(0, 4, name="muon_nearGenA_n", label=r"$N_{\mu^{PF}}$"),
                   lambda objs, mask: ak.num(matched(objs["muons"], objs["genAs_toMu"], 0.5))),
        ],
    ),
    "muon_genA_dR": h.Histogram(
        [
            # dR(mu, nearest gen A)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="muon_genA_dR"),
                   lambda objs, mask: dR(objs["muons"], objs["genAs"]))
        ],
    ),
    # pfmuon-genMuon
    "muon_genMu_dR": h.Histogram(
        [
            # dR(mu, nearest gen mu)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="muon_genMu_dR"),
                   lambda objs, mask: dR(objs["muons"], objs["genMus"]))
        ],
    ),
    # dsamuon
    "dsaMuon_n": obj_attr("dsaMuons", "n"),
    "dsaMuon_pt":obj_attr("dsaMuons", "pt", xmax=500),
    "dsaMuon_dxy":obj_attr("dsaMuons", "dxy"),
    "dsaMuon_eta_phi": obj_eta_phi("dsaMuons"),
    "dsaMuon_absD0": obj_attr("dsaMuons", "dxy", absval=True, xmax=500),
    "dsaMuon_absD0_lowRange": obj_attr("dsaMuons", "dxy", absval=True, xmax=10),
    "dsaMuon_nearGenA_n": h.Histogram(
        [
            # number of muons within dR=0.5 of a genA that decays to muons
            h.Axis(hist.axis.Integer(0, 10, name="dsaMuon_nearGenA_n"),
                   lambda objs, mask: ak.num(matched(objs["dsaMuons"], objs["genAs_toMu"], 0.5))),
        ],
    ),
    # dsamuon-genA
    "dsaMuon_nearGenA_n_genA_lxy": h.Histogram(
        [
            # lxy of dark photon that decays to dsaMuons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toMu"])),
            # number of dsaMuons within dR=0.5 of a genA that decays to muons
            h.Axis(hist.axis.Integer(0, 4, name="dsaMuon_nearGenA_n", label=r"$N_{\mu^{DSA}}$"),
                   lambda objs, mask: ak.num(matched(objs["dsaMuons"], objs["genAs_toMu"], 0.5))),
        ],
    ),
    "dsaMuon_genA_dR": h.Histogram(
        [
            # dR(dsa mu, nearest gen A)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="dsaMuon_genA_dR"),
                   lambda objs, mask: dR(objs["dsaMuons"], objs["genAs"]))
        ],
    ),
    # dsaMuon-genMuon
    "dsaMuon_genMu_dR": h.Histogram(
        [
            # dR(dsa mu, nearest gen mu)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="dsaMuon_genMu_dR"),
                   lambda objs, mask: dR(objs["dsaMuons"], objs["genMus"]))
        ],
    ),
    # lj
    "lj_n": obj_attr("ljs", "n"),
    "egm_lj_n": obj_attr("egm_ljs", "n"),
    "mu_lj_n": obj_attr("mu_ljs", "n"),
    "lj_pt": obj_attr("ljs", "pt", xmax=700),
    "lj0_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 400, name="lj0_pt",
                                     label="Leading lepton jet pT [GeV]"),
                   lambda objs, mask: objs["ljs"][mask, 0].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 0,
    ),
    "lj1_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 400, name="lj1_pt",
                                     label="Subleading lepton jet pT [GeV]"),
                   lambda objs, mask: objs["ljs"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj0_e": h.Histogram(
        [
            h.Axis(hist.axis.Regular(350, 0, 700, name="lj_e",
                                     label="Leading lepton jet E [GeV]"),
                   lambda objs, mask: objs["ljs"][mask, 0].energy),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 0,
    ),
    "lj1_e": h.Histogram(
        [
            h.Axis(hist.axis.Regular(350, 0, 700, name="lj_e",
                                     label="Subleading lepton jet E [GeV]"),
                   lambda objs, mask: objs["ljs"][mask, 1].energy),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj0_dRSpread": h.Histogram(
        [
            h.Axis(hist.axis.Regular(250, 0, 1.0, name="lj0_dRSpread",
                                     label="Leading lepton jet dRSpread"),
                   lambda objs, mask: objs["ljs"][mask, 0].dRSpread),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 0,
    ),
    "lj1_dRSpread": h.Histogram(
        [
            h.Axis(hist.axis.Regular(250, 0, 1.0, name="lj1_dRSpread",
                                     label="Subleading lepton jet dRSpread"),
                   lambda objs, mask: objs["ljs"][mask, 1].dRSpread),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_eta_phi": obj_eta_phi("ljs"),
    "lj_electronN": obj_attr("ljs", "electron_n", xmax=10, nbins=10),
    "lj_photonN": obj_attr("ljs", "photon_n", xmax=10, nbins=10),
    "lj_muonN": obj_attr("ljs", "muon_n", xmax=10, nbins=10),
    "lj_dsaMuN": obj_attr("ljs", "dsaMu_n", xmax=10, nbins=10),
    "lj_pfMuN": obj_attr("ljs", "pfMu_n", xmax=10, nbins=10),
    "mu_lj_pt": obj_attr("mu_ljs", "pt", xmax=700),
    "mu_lj_muonN": obj_attr("mu_ljs", "muon_n", xmax=10, nbins=10),
    "mu_lj_pfMu_n": obj_attr("mu_ljs", "pfMu_n", xmax=10, nbins=10),
    "mu_lj_dsaMu_n": obj_attr("mu_ljs", "dsaMu_n", xmax=10, nbins=10),
    "egm_lj_pt": obj_attr("egm_ljs", "pt", xmax=700),
    "egm_lj_electronN": obj_attr("egm_ljs", "electron_n", xmax=10, nbins=10),
    "egm_lj_photonN": obj_attr("egm_ljs", "photon_n", xmax=10, nbins=10),
    "egm_lj_electron_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name=r"egm- type LJ e pT (GeV)"),
                   lambda objs, mask: objs["egm_ljs"].electrons.pt),
        ],
    ),
    "egm_lj_photon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name=r"egm- type LJ $\gamma$ pT (GeV)"),
                   lambda objs, mask: objs["egm_ljs"].photons.pt),
        ],
    ),
    "egm_lj_electron_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2, name=r"egm- type LJ e dxy (cm)"),
                   lambda objs, mask: abs(objs["egm_ljs"].electrons.dxy)),
        ],
    ),
    "egm_lj_electron_min_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2, name=r"egm- type LJ e min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_max_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2, name=r"egm- type LJ e max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_dxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1, name=r"egm- type LJ e dxy (cm)"),
                   lambda objs, mask: abs(objs["egm_ljs"].electrons.dxy)),
        ],
    ),
    "egm_lj_electron_min_dxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1, name=r"egm- type LJ e min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_max_dxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1, name=r"egm- type LJ e max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, .4, name=r"egm- type LJ e dxy (cm)"),
                   lambda objs, mask: abs(objs["egm_ljs"].electrons.dxy)),
        ],
    ),
    "egm_lj_electron_min_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, .4, name=r"egm- type LJ e min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_max_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, .4, name=r"egm- type LJ e max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_dxy_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, .2, name=r"egm- type LJ e dxy (cm)"),
                   lambda objs, mask: abs(objs["egm_ljs"].electrons.dxy)),
        ],
    ),
    "egm_lj_electron_min_dxy_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, .2, name=r"egm- type LJ e min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_max_dxy_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, .2, name=r"egm- type LJ e max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_dxy_XXXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, .1, name=r"egm- type LJ e dxy (cm)"),
                   lambda objs, mask: abs(objs["egm_ljs"].electrons.dxy)),
        ],
    ),
    "egm_lj_electron_min_dxy_XXXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, .1, name=r"egm- type LJ e min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_max_dxy_XXXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, .1, name=r"egm- type LJ e max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_dxy_XXXXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, .01, name=r"egm- type LJ e dxy (cm)"),
                   lambda objs, mask: abs(objs["egm_ljs"].electrons.dxy)),
        ],
    ),
    "egm_lj_electron_min_dxy_XXXXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, .01, name=r"egm- type LJ e min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_max_dxy_XXXXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, .01, name=r"egm- type LJ e max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["egm_ljs"].electrons.dxy), axis=-1)),
        ],
    ),
    "egm_lj_electron_lostHits": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name=r"egm- type LJ e lostHits"),
                   lambda objs, mask: objs["egm_ljs"].electrons.lostHits),
        ],
    ),
    "mu_lj_muon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name=r"$\mu$- type LJ $\mu$ pT (GeV)"),
                   lambda objs, mask: objs["mu_ljs"].muons.pt),
        ],
    ),
    "mu_lj_pfMuon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name=r"$\mu$- type LJ PF $\mu$ pT (GeV)"),
                   lambda objs, mask: objs["mu_ljs"].pfMuons.pt),
        ],
    ),
    "mu_lj_dsaMuon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name=r"$\mu$- type LJ DSA $\mu$ pT (GeV)"),
                   lambda objs, mask: objs["mu_ljs"].dsaMuons.pt),
        ],
    ),
    "mu_lj_muon_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 50, name=r"$\mu$- type LJ $\mu$ dxy (cm)"),
                   lambda objs, mask: abs(objs["mu_ljs"].muons.dxy)),
        ],
    ),
    "mu_lj_pfMuon_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 20, name=r"$\mu$- type LJ PF $\mu$ dxy (cm)"),
                   lambda objs, mask: abs(objs["mu_ljs"].pfMuons.dxy)),
        ],
    ),
    "mu_lj_dsaMuon_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 50, name=r"$\mu$- type LJ DSA $\mu$ dxy (cm)"),
                   lambda objs, mask: abs(objs["mu_ljs"].dsaMuons.dxy)),
        ],
    ),
    "mu_lj_muon_dxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 10, name=r"$\mu$- type LJ $\mu$ dxy (cm)"),
                   lambda objs, mask: abs(objs["mu_ljs"].muons.dxy)),
        ],
    ),
    "mu_lj_pfMuon_dxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 5, name=r"$\mu$- type LJ PF $\mu$ dxy (cm)"),
                   lambda objs, mask: abs(objs["mu_ljs"].pfMuons.dxy)),
        ],
    ),
    "mu_lj_dsaMuon_dxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 10, name=r"$\mu$- type LJ DSA $\mu$ dxy (cm)"),
                   lambda objs, mask: abs(objs["mu_ljs"].dsaMuons.dxy)),
        ],
    ),
    "mu_lj_muon_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1, name=r"$\mu$- type LJ $\mu$ dxy (cm)"),
                   lambda objs, mask: abs(objs["mu_ljs"].muons.dxy)),
        ],
    ),
    "mu_lj_pfMuon_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1, name=r"$\mu$- type LJ PF $\mu$ dxy (cm)"),
                   lambda objs, mask: abs(objs["mu_ljs"].pfMuons.dxy)),
        ],
    ),
    "mu_lj_dsaMuon_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1, name=r"$\mu$- type LJ DSA $\mu$ dxy (cm)"),
                   lambda objs, mask: abs(objs["mu_ljs"].dsaMuons.dxy)),
        ],
    ),
    "mu_lj_muon_min_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 50, name=r"$\mu$- type LJ $mu$ min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["mu_ljs"].muons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_dsaMuon_min_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 50, name=r"$\mu$- type LJ DSA $\mu$ min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["mu_ljs"].dsaMuons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_pfMuon_min_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 5, name=r"$\mu$- type LJ PF $\mu$ min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["mu_ljs"].pfMuons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_muon_min_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1, name=r"$\mu$- type LJ $mu$ min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["mu_ljs"].muons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_dsaMuon_min_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1, name=r"$\mu$- type LJ DSA $\mu$ min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["mu_ljs"].dsaMuons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_pfMuon_min_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 0.1, name=r"$\mu$- type LJ PF $\mu$ min dxy (cm)"),
                   lambda objs, mask: ak.min(abs(objs["mu_ljs"].pfMuons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_muon_max_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 50, name=r"$\mu$- type LJ $\mu$ max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["mu_ljs"].muons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_dsaMuon_max_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 50, name=r"$\mu$- type LJ DSA $\mu$ max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["mu_ljs"].dsaMuons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_pfMuon_max_dxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 5, name=r"$\mu$- type LJ PF $\mu$ max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["mu_ljs"].pfMuons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_muon_max_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1, name=r"$\mu$- type LJ $\mu$ max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["mu_ljs"].muons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_dsaMuon_max_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1, name=r"$\mu$- type LJ DSA $\mu$ max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["mu_ljs"].dsaMuons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_pfMuon_max_dxy_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 0.1, name=r"$\mu$- type LJ PF $\mu$ max dxy (cm)"),
                   lambda objs, mask: ak.max(abs(objs["mu_ljs"].pfMuons.dxy), axis=-1)),
        ],
    ),
    "mu_lj_pfMu_nTrackerLayers": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ PF $\mu$ nTrackerLayers"),
                   lambda objs, mask: objs["mu_ljs"].pfMuons.nTrackerLayers),
        ],
    ),
    "mu_lj_dsaMu_trkNumPlanes": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ DSA $\mu$ trkNumPlanes"),
                   lambda objs, mask: objs["mu_ljs"].dsaMuons.trkNumPlanes),
        ],
    ),
    "mu_lj_dsaMu_trkNumHits": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ DSA $\mu$ trkNumHits"),
                   lambda objs, mask: objs["mu_ljs"].dsaMuons.trkNumHits),
        ],
    ),
    "mu_lj_dsaMu_trkNumDTHits": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ DSA $\mu$ trkDTHHits"),
                   lambda objs, mask: objs["mu_ljs"].dsaMuons.trkNumDTHits),
        ],
    ),
    "mu_lj_dsaMu_trkNumCSCHits": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ DSA $\mu$ trkNumCSCHits"),
                   lambda objs, mask: objs["mu_ljs"].dsaMuons.trkNumCSCHits),
        ],
    ),
    "mu_lj_dsaMu_nSegments": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ DSA $\mu$ nSegments"),
                   lambda objs, mask: objs["mu_ljs"].dsaMuons.nSegments),
        ],
    ),
    "mu_lj_dsaMu_nDTSegments": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ DSA $\mu$ nDTSegments"),
                   lambda objs, mask: objs["mu_ljs"].dsaMuons.nDTSegments),
        ],
    ),
    "mu_lj_dsaMu_nCSCSegments": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ DSA $\mu$ nCSCSegments"),
                   lambda objs, mask: objs["mu_ljs"].dsaMuons.nCSCSegments),
        ],
    ),
    "mu_lj_pfMu_trkNumPlanes": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ PF $\mu$ trkNumPlanes"),
                   lambda objs, mask: objs["mu_ljs"].pfMuons.trkNumPlanes),
        ],
    ),
    "mu_lj_pfMu_trkNumHits": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ PF $\mu$ trkNumHits"),
                   lambda objs, mask: objs["mu_ljs"].pfMuons.trkNumHits),
        ],
    ),
    "mu_lj_pfMu_trkNumDTHits": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ PF $\mu$ trkNumDTHHits"),
                   lambda objs, mask: objs["mu_ljs"].pfMuons.trkNumDTHits),
        ],
    ),
    "mu_lj_pfMu_trkNumCSCHits": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ PF $\mu$ trkNumCSCits"),
                   lambda objs, mask: objs["mu_ljs"].pfMuons.trkNumCSCHits),
        ],
    ),
    "mu_lj_pfMu_trkNumPixelHits": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ PF $\mu$ trkNumPixelHits"),
                   lambda objs, mask: objs["mu_ljs"].pfMuons.trkNumPixelHits),
        ],
    ),
    "mu_lj_pfMu_trkNumTrkLayers": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, 40, name=r"$\mu$- type LJ PF $\mu$ trkNumTrkLayers"),
                   lambda objs, mask: objs["mu_ljs"].pfMuons.trkNumTrkLayers),
        ],
    ),
    "mu_lj_muon_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="mu_lj_mu_eta"),
                   lambda objs, mask: objs["mu_ljs"].muons.eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="mu_lj_mu_phi"),
                   lambda objs, mask: objs["mu_ljs"].muons.phi),
        ],
    ),
    "lj_electronPhotonN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_electronPhotonN"),
                   lambda objs, mask: objs["ljs"].electron_n + objs["ljs"].photon_n),
        ],
    ),
    # pfelectron-lj
    "electron_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="electron_lj_dR"),
                   lambda objs, mask: dR(objs["electrons"], objs["ljs"]))
        ],
    ),
    
    "electron_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 1.0, name="electron_lj_dR_lowRange"),
                   lambda objs, mask: dR(objs["electrons"], objs["ljs"]))
        ],
    ),
    # pfelectron-egm lj
    "electron_egmLj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 6, name="electron_egmLj_dR"),
                   lambda objs, mask: dR(objs["electrons"], objs["egm_ljs"]))
        ],
    ),
    # pfphoton-lj
    "photon_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="photon_lj_dR"),
                   lambda objs, mask: dR(objs["photons"], objs["ljs"]))
        ],
    ),
    "photon_lj_dR_lowRange": h.Histogram(
        [
            # dR(photon, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 1.0, name="photon_lj_dR_lowRange"),
                   lambda objs, mask: dR(objs["photons"], objs["ljs"]))
        ],
    ),
    "photon_lj_dR_reallyLowRange": h.Histogram(
        [
            # dR(photon, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 0.1, name="photon_lj_dR_reallyLowRange"),
                   lambda objs, mask: dR(objs["photons"], objs["ljs"]))
        ],
    ),
    # pfmuon-lj
    "muon_lj_dR": h.Histogram(
        [
            # dR(mu, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="muon_lj_dR"),
                   lambda objs, mask: dR(objs["muons"], objs["ljs"]))
        ],
    ),
    "muon_lj_dR_lowRange": h.Histogram(
        [
            # dR(mu, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 1.0, name="muon_lj_dR_lowRange"),
                   lambda objs, mask: dR(objs["muons"], objs["ljs"]))
        ],
    ),
    # dsamuon-lj
    "dsaMuon_lj_dR": h.Histogram(
        [
            # dR(dsa mu, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="dsaMuon_lj_dR"),
                   lambda objs, mask: dR(objs["dsaMuons"], objs["ljs"]))
        ],
    ),
    "dsaMuon_lj_dR_lowRange": h.Histogram(
        [
            # dR(dsa mu, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="dsaMuon_lj_dR_lowRange"),
                   lambda objs, mask: dR(objs["dsaMuons"], objs["ljs"]))
        ],
    ),
    # lj-lj
    "lj_lj_absdphi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="|$\Delta\phi$| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].phi - objs["ljs"][mask, 0].phi)),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_absdR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 6, name="|$\Delta$R| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, 1].delta_r(objs["ljs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_absdeta": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 6, name="|$\Delta\eta$| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].eta - objs["ljs"][mask, 0].eta)),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_invmass": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 1200, name="ljlj_mass",
                                     label=r"Invariant Mass ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, :2].sum().mass),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_invmass_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="ljlj_mass",
                                     label=r"InvMass($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, :2].sum().mass),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 1.0, 2.0, name="lj_lj_ptRatio",
                   label="Leading LJ PT / Subleading LJ PT"),
                   lambda objs, mask: objs["ljs"][mask, 0].pt / objs["ljs"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    # matchedjet
    # matchedjet
    "matched_jet_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 800, name="matched_jet_pt",
                   label="Matched Jet PT [GeV]"),
                   lambda objs, mask:  objs["ljs"].matched_jet.pt),
        ],
    ),
    "mu_matched_jet_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 800, name="mu_matched_jet_pt",
                   label="Mu Matched Jet PT [GeV]"),
                   lambda objs, mask:  objs["mu_ljs"].matched_jet.pt),
        ],
    ),
    "pfmu_matched_jet_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 800, name="pfmu_matched_jet_pt",
                   label="PF Mu Matched Jet PT [GeV]"),
                   lambda objs, mask:  objs["pfmu_ljs"].matched_jet.pt),
        ],
    ),
    "dsamu_matched_jet_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 800, name="dsamu_matched_jet_pt",
                   label="DSA Mu Matched Jet PT [GeV]"),
                   lambda objs, mask:  objs["dsamu_ljs"].matched_jet.pt),
        ],
    ),
    "egm_matched_jet_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 800, name="egm_matched_jet_pt",
                   label="EGM Matched Jet PT [GeV]"),
                   lambda objs, mask:  objs["egm_ljs"].matched_jet.pt),
        ],
    ),
    "electron_matched_jet_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 800, name="electron_matched_jet_pt",
                   label="Electron Matched Jet PT [GeV]"),
                   lambda objs, mask:  objs["electron_ljs"].matched_jet.pt),
        ],
    ),
    "photon_matched_jet_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 800, name="photon_matched_jet_pt",
                   label="Photon Matched Jet PT [GeV]"),
                   lambda objs, mask:  objs["photon_ljs"].matched_jet.pt),
        ],
    ),
    "matched_jet_lepfraction": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 1, name="matched_jet_lepfraction",
                   label="Matched Jet Lepton Fraction"),
                   lambda objs, mask:  objs["ljs"].lepton_fraction),
        ],
    ),
    "mu_matched_jet_lepfraction": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 1, name="mu_matched_jet_lepfraction",
                   label="Mu Matched Jet Lepton Fraction"),
                   lambda objs, mask:  objs["mu_ljs"].lepton_fraction),
        ],
    ),
    "pfmu_matched_jet_lepfraction": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 1, name="pfmu_matched_jet_lepfraction",
                   label="PF Mu Matched Jet Lepton Fraction"),
                   lambda objs, mask:  objs["pfmu_ljs"].lepton_fraction),
        ],
    ),
    "dsamu_matched_jet_lepfraction": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 1, name="dsamu_matched_jet_lepfraction",
                   label="DSA Mu Matched Jet Lepton Fraction"),
                   lambda objs, mask:  objs["dsamu_ljs"].lepton_fraction),
        ],
    ),
    "egm_matched_jet_lepfraction": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 1, name="egm_matched_jet_lepfraction",
                   label="EGM Matched Jet Lepton Fraction"),
                   lambda objs, mask:  objs["egm_ljs"].lepton_fraction),
        ],
    ),
    "electron_matched_jet_lepfraction": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 1, name="electron_matched_jet_lepfraction",
                   label="Electron Matched Jet Lepton Fraction"),
                   lambda objs, mask:  objs["electron_ljs"].lepton_fraction),
        ],
    ),
    "photon_matched_jet_lepfraction": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 1, name="photon_matched_jet_lepfraction",
                   label="Photon Matched Jet Lepton Fraction"),
                   lambda objs, mask:  objs["photon_ljs"].lepton_fraction),
        ],
    ),
    # matchedjet-lj
    "matched_jet_lj_dR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.4, name="matched_jet_dR",
                   label="dR(LJ, Matched Jet)"),
                   lambda objs, mask: objs["ljs"].dR_matched_jet),
        ],
    ),
    "mu_matched_jet_lj_dR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.4, name="mu_matched_jet_lj_dR",
                   label="dR(Mu-LJ, Mu-Matched Jet)"),
                   lambda objs, mask: objs["mu_ljs"].dR_matched_jet),
        ],
    ),
    "pfmu_matched_jet_lj_dR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.4, name="pfmu_matched_jet_lj_dR",
                   label="dR(PF Mu-LJ, Mu-Matched Jet)"),
                   lambda objs, mask: objs["pfmu_ljs"].dR_matched_jet),
        ],
    ),
    "dsamu_matched_jet_lj_dR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.4, name="dsamu_matched_jet_lj_dR",
                   label="dR(DSA Mu-LJ, DSA Mu-Matched Jet)"),
                   lambda objs, mask: objs["dsamu_ljs"].dR_matched_jet),
        ],
    ),
    "egm_matched_jet_lj_dR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.4, name="egm_matched_jet_lj_dR",
                   label="dR(EGM-LJ, EGM-Matched Jet)"),
                   lambda objs, mask: objs["egm_ljs"].dR_matched_jet),
        ],
    ), 
    "electron_matched_jet_lj_dR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.4, name="electron_matched_jet_lj_dR",
                   label="dR(Electron-LJ, Electron EGM-Matched Jet)"),
                   lambda objs, mask: objs["electron_ljs"].dR_matched_jet),
        ],
    ), 
    "photon_matched_jet_lj_dR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.4, name="photon_matched_jet_lj_dR",
                   label="dR(Photon-LJ, Photon EGM-Matched Jet)"),
                   lambda objs, mask: objs["photon_ljs"].dR_matched_jet),
        ],
    ),
    # lj isolation
    "lj_isolation": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2, name="lj_isolation",
                   label="LJ Isolation"),
                   lambda objs, mask:  objs["ljs"].isolation),
        ],
    ),
    "lj_isolation_zoom": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.2, name="lj_isolation",
                   label="LJ Isolation"),
                   lambda objs, mask:  objs["ljs"].isolation),
        ],
    ),
    "mu_lj_isolation": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2, name="mu_lj_isolation",
                   label="Mu-LJ Isolation"),
                   lambda objs, mask:  objs["mu_ljs"].isolation),
        ],
    ),
    "mu_lj_isolation_zoom": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.2, name="mu_lj_isolation",
                   label="Mu-LJ Isolation"),
                   lambda objs, mask:  objs["mu_ljs"].isolation),
        ],
    ),
    "pfmu_lj_isolation": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2, name="pfmu_lj_isolation",
                   label="PF Mu-LJ Isolation"),
                   lambda objs, mask:  objs["pfmu_ljs"].isolation),
        ],
    ),
    "pfmu_lj_isolation_zoom": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.2, name="pfmu_lj_isolation",
                   label="PF Mu-LJ Isolation"),
                   lambda objs, mask:  objs["pfmu_ljs"].isolation),
        ],
    ),
    "dsamu_lj_isolation": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2, name="dsamu_lj_isolation",
                   label="DSA Mu-LJ Isolation"),
                   lambda objs, mask:  objs["dsamu_ljs"].isolation),
        ],
    ),
    "dsamu_lj_isolation_zoom": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.2, name="dsamu_lj_isolation",
                   label="DSA Mu-LJ Isolation"),
                   lambda objs, mask:  objs["dsamu_ljs"].isolation),
        ],
    ),
    "egm_lj_isolation": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2, name="egm_lj_isolation",
                   label="EGM-LJ Isolation"),
                   lambda objs, mask:  objs["egm_ljs"].isolation),
        ],
    ),
    "egm_lj_isolation_zoom": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.2, name="egm_lj_isolation",
                   label="EGM-LJ Isolation"),
                   lambda objs, mask:  objs["egm_ljs"].isolation),
        ],
    ),
    "electron_lj_isolation": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2, name="electron_lj_isolation",
                   label="Electron-LJ Isolation"),
                   lambda objs, mask:  objs["electron_ljs"].isolation),
        ],
    ),
    "electron_lj_isolation_zoom": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.2, name="electron_lj_isolation",
                   label="Electron-LJ Isolation"),
                   lambda objs, mask:  objs["electron_ljs"].isolation),
        ],
    ),
    "photon_lj_isolation": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2, name="photon_lj_isolation",
                   label="Photon-LJ Isolation"),
                   lambda objs, mask:  objs["photon_ljs"].isolation),
        ],
    ),
    "photon_lj_isolation_zoom": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 0.2, name="photon_lj_isolation",
                   label="Photon-LJ Isolation"),
                   lambda objs, mask:  objs["photon_ljs"].isolation),
        ],
    ),
    "lj_isolation_sum": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2, name="lj_isolation_sum",
                   label="LJ1 Isolation + LJ2 Isolation"),
                   lambda objs, mask:  (objs["ljs"][mask, 1].isolation + objs["ljs"][mask, 0].isolation)),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "dlj_isolation": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2, name="lj_isolation",
                   label="|LJ1 Isolation - LJ2 Isolation|"),
                   lambda objs, mask:  abs(objs["ljs"][mask, 1].isolation - objs["ljs"][mask, 0].isolation)),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_isolation_ratio1": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 10, name="lj_isolation_ratio1",
                   label="LJ1 Isolation / LJ2 Isolation"),
                   lambda objs, mask:  (objs["ljs"][mask, 0].isolation / objs["ljs"][mask, 1].isolation)),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_isolation_ratio2": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 10, name="lj_isolation_ratio2",
                   label="LJ2 Isolation / LJ1 Isolation"),
                   lambda objs, mask:  (objs["ljs"][mask, 1].isolation / objs["ljs"][mask, 0].isolation)),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    # ABCD plane
    "lj_lj_absdphi_invmass": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="|$\Delta\phi$| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].phi - objs["ljs"][mask, 0].phi)),
            h.Axis(hist.axis.Regular(100, 0, 1200, name="ljlj_mass",
                                     label=r"Invariant Mass ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, :2].sum().mass),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_absdphi_absdR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="|$\Delta\phi$| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].phi - objs["ljs"][mask, 0].phi)),
            h.Axis(hist.axis.Regular(100, 0, 6, name="|$\Delta$R| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, 1].delta_r(objs["ljs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_absdphi_absdeta": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="|$\Delta\phi$| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].phi - objs["ljs"][mask, 0].phi)),
            h.Axis(hist.axis.Regular(100, 0, 6, name="|$\Delta\eta$| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].eta - objs["ljs"][mask, 0].eta)),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_absdphi_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="|$\Delta\phi$| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].phi - objs["ljs"][mask, 0].phi)),
            h.Axis(hist.axis.Regular(100, 1.0, 2.0, name="lj_lj_ptRatio",
                                     label="Leading LJ PT / Subleading LJ PT"),
                   lambda objs, mask: objs["ljs"][mask, 0].pt / objs["ljs"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_absdR_absdeta": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 6, name="|$\Delta$R| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, 1].delta_r(objs["ljs"][mask, 0])),
            h.Axis(hist.axis.Regular(100, 0, 6, name="|$\Delta\eta$| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].eta - objs["ljs"][mask, 0].eta)),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_absdR_invmass": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 6, name="|$\Delta$R| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, 1].delta_r(objs["ljs"][mask, 0])),
            h.Axis(hist.axis.Regular(100, 0, 1200, name="ljlj_mass",
                                     label=r"Invariant Mass ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, :2].sum().mass),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_absdR_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 6, name="|$\Delta$R| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, 1].delta_r(objs["ljs"][mask, 0])),
            h.Axis(hist.axis.Regular(100, 1.0, 2.0, name="lj_lj_ptRatio",
                   label="Leading LJ PT / Subleading LJ PT"),
                   lambda objs, mask: objs["ljs"][mask, 0].pt / objs["ljs"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_absdeta_invmass": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 6, name="|$\Delta\eta$| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].eta - objs["ljs"][mask, 0].eta)),
            h.Axis(hist.axis.Regular(100, 0, 1200, name="ljlj_mass",
                                     label=r"Invariant Mass ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, :2].sum().mass),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_absdeta_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 6, name="|$\Delta\eta$| ($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].eta - objs["ljs"][mask, 0].eta)),
            h.Axis(hist.axis.Regular(100, 1.0, 2.0, name="lj_lj_ptRatio",
                   label="Leading LJ PT / Subleading LJ PT"),
                   lambda objs, mask: objs["ljs"][mask, 0].pt / objs["ljs"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    # gen
    "gen_abspid": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 40, name="gen_abspid"),
                   lambda objs, mask: abs(objs["gens"].pdgId)),
        ],
    ),
    # genelectron
    "genE_n": obj_attr("genEs", "n"),
    "genE_pt": obj_attr("genEs", "pt"),
    "genE_pt_highRange": obj_attr("genEs", "pt", xmax=700),
    "genE_dxy": obj_attr("genEs", "dxy", absval=True, xmax=10, nbins=100),
    "genE_dxy_lowRange": obj_attr("genEs", "dxy", absval=True, xmax=1, nbins=100),
    "genE_dxy_XLowRange": obj_attr("genEs", "dxy", absval=True, xmax=0.1, nbins=100),
    "genE_dxy_XXLowRange": obj_attr("genEs", "dxy", absval=True, xmax=0.01, nbins=100),
    "genE0_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE0_pt",
                                     label=r"Leading gen-level electron $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genEs"][mask, 0].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 0,
    ),
    "genE0_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(70, 0, 700, name="genE_pt",
                                     label=r"Leading gen-level electron $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genEs"][mask, 0].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 0,
    ),
    "genE1_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE1_pt",
                                     label=r"Sub-leading gen-level electron $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genEs"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE1_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(70, 0, 700, name="genE_pt",
                                     label=r"Sub-leading gen-level electron $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genEs"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE_eta_phi": obj_eta_phi("genEs"),
    "genE_parent_absPdgId": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 50, name="genE_parent_absPdgId"),
                   lambda objs, mask: abs(objs["genEs"].parent.pdgId)),
        ],
    ),
    # genelectron-genelectron
    "genE_genE_dR": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(100, 0, 1.0, name="genE_genE_dR",
                                     label=r"$\Delta R$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE_genE_dR_lowRange": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(75, 0, 0.5, name="genE_genE_dR_lowRange",
                                     label=r"$\Delta R$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE_genE_dR_matched_Photon": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(50, 0, 0.5, name="genE_genE_dR",
                                     label=r"$\Delta R$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: (ak.num(objs["genEs"]) > 1) & (ak.num(derived_objs['genAs_toE_matched_photons'](objs, 0.4)) > 0),
    ),
    "genE_genE_dR_matched_Photon_lowRange": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(50, 0, 0.1, name="genE_genE_dR",
                                     label=r"$\Delta R$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: (ak.num(objs["genEs"]) > 1) & (ak.num(derived_objs['genAs_toE_matched_photons'](objs, 0.4)) > 0),
    ),
    "genE_genE_dR_matched_Photon_XLowRange": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(50, 0, 0.04, name="genE_genE_dR",
                                     label=r"$\Delta R$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: (ak.num(objs["genEs"]) > 1) & (ak.num(derived_objs['genAs_toE_matched_photons'](objs, 0.4)) > 0),
    ),
    "genE_genE_dR_XLowRange": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(50, 0, 0.1, name="genE_genE_dR_lowRange",
                                     label=r"$\Delta R$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE_genE_dR_XXLowRange": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(50, 0, 0.04, name="genE_genE_dR_lowRange",
                                     label=r"$\Delta R$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    
    "genE_genE_dEta": h.Histogram(
        [
            # abs(dEta(subleading gen E, leading gen E))
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genE_genE_dEta",
                                     label=r"$\Delta\, \eta$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: abs(objs["genEs"][mask, 1].eta
                                          - objs["genEs"][mask, 0].eta)),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE_genE_dPhi": h.Histogram(
        [
            # abs(dEta(subleading gen E, leading gen E))
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genE_genE_dPhi",
                                     label=r"$\Delta\, \phi$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: abs(objs["genEs"][mask, 1].phi
                                          - objs["genEs"][mask, 0].phi)),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE_genE_dPt": h.Histogram(
        [
            # abs(dEta(subleading gen E, leading gen E))
            h.Axis(hist.axis.Regular(50, 0, 100, name="genE_genE_dPt",
                                     label=r"$\Delta\, p_T$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: abs(objs["genEs"][mask, 1].pt
                                          - objs["genEs"][mask, 0].pt)),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE_genE_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE_genE_pt"),
                   lambda objs, mask: objs["genEs"][mask, :2].sum().pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    # genmuon
    "genMu_n": obj_attr("genMus", "n"),
    "genMu_pt": obj_attr("genMus", "pt"),
    "genMu_pt_highRange": obj_attr("genMus", "pt", xmax=700),
    "genMu_dxy": obj_attr("genMus", "dxy", absval=True, xmax=10, nbins=100),
    "genMu_dxy_lowRange": obj_attr("genMus", "dxy", absval=True, xmax=1, nbins=100),
    "genMu_dxy_XLowRange": obj_attr("genMus", "dxy", absval=True, xmax=0.1, nbins=100),
    "genMu_dxy_XXLowRange": obj_attr("genMus", "dxy", absval=True, xmax=0.01, nbins=100),
    "genMu0_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 200, name="genMu0_pt",
                                     label=r"Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 0].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 0,
    ),
    "genMu0_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(70, 0, 700, name="genMu0_pt",
                                     label=r"Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 0].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 0,
    ),
    "genMu1_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genMu1_pt",
                                     label=r"Sub-leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu1_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(70, 0, 700, name="genMu1_pt",
                                     label=r"Sub-leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_eta_phi": obj_eta_phi("genMus"),
    "genMu_parent_absPdgId": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 50, name="genMu_parent_absPdgId"),
                   lambda objs, mask: abs(objs["genMus"].parent.pdgId)),
        ],
    ),
    # genmuon-genmuon
    "genMu_genMu_dR": h.Histogram(
        [
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(
                       objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_genMu_dR_lowRange": h.Histogram(
        [
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(50, 0, 0.5, name="genMu_genMu_dR_lowRange",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(
                       objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_genMu_dR_XLowRange": h.Histogram(
        [
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(50, 0, 0.1, name="genMu_genMu_dR_lowRange",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(
                       objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_genMu_dR_XXLowRange": h.Histogram(
        [
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(50, 0, 0.04, name="genMu_genMu_dR_lowRange",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(
                       objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_genMu_dEta": h.Histogram(
        [
            # abs(dEta(subleading gen Mu, leading gen Mu))
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genMu_genMu_dEta",
                                     label=r"$\Delta\, \eta$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: abs(objs["genMus"][mask, 1].eta
                                          - objs["genMus"][mask, 0].eta)),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_genMu_dPhi": h.Histogram(
        [
            # abs(dEta(subleading gen Mu, leading gen Mu))
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genMu_genMu_dPhi",
                                     label=r"$\Delta\, \phi$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: abs(objs["genMus"][mask, 1].phi
                                          - objs["genMus"][mask, 0].phi)),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_genMu_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genMu_genMu_pt"),
                   lambda objs, mask: objs["genMus"][mask, :2].sum().pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    #dsamuon-genAs_toMu
    "dsamuon_absd0_genAs_toMu_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 200, name="dsamuon_absd0",
                                     label=r"dsa muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"].d0)),
            h.Axis(hist.axis.Regular(25, 0,400, name="genA_lxy"),
                   #Added the function ak.ones_like to match delta R array with the d0 array.
                   lambda objs, mask: lxy(objs["genAs_toMu"])[:,0]*ak.ones_like(objs["dsaMuons"].d0)),
        ],
        evt_mask=lambda objs: ak.num(objs["genAs_toMu"]) > 0,
    ),
    #dsamuon-genmuon
    "dsaMuon_absD0_genMus_dR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 200, name="dsaMuon_absD0",
                                     label=r"DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"].d0)),
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(25, 0, 0.4, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   #Added the function ak.ones_like to match delta R array with the d0 array.
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(
                       objs["genMus"][mask, 0])*ak.ones_like(objs["dsaMuons"].d0)),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1
    ),
    "leadingDsaMuon_absD0_genMus_dR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 200, name="dsaMuon_absD0",
                                     label=r"Leading DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"][mask, 0].d0)),
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(25, 0, 0.4, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ((ak.num(objs["genMus"]) > 1) & (ak.num(objs["dsaMuons"]) > 0))
    ),
    "subLeadingDsaMuon_absD0_genMus_dR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 200, name="dsaMuon_absD0",
                                     label=r"Sub Leading DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"][mask, 1].d0)),
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(25, 0, 0.4, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ((ak.num(objs["genMus"]) > 1) & (ak.num(objs["dsaMuons"]) > 1))
    ),
    "dsaMuon_absD0_genMus_dR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 200, name="dsaMuon_absD0",
                                     label=r"DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"].d0)),
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(25, 0, 0.1, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   #Added the function ak.ones_like to match delta R array with the d0 array.
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(
                       objs["genMus"][mask, 0])*ak.ones_like(objs["dsaMuons"].d0)),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1
    ),
    "leadingDsaMuon_absD0_genMus_dR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 200, name="dsaMuon_absD0",
                                     label=r"Leading DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"][mask, 0].d0)),
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(25, 0, 0.1, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ((ak.num(objs["genMus"]) > 1) & (ak.num(objs["dsaMuons"]) > 0))
    ),
    "subLeadingDsaMuon_absD0_genMus_dR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 200, name="dsaMuon_absD0",
                                     label=r"Sub Leading DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"][mask, 1].d0)),
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(25, 0, 0.1, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ((ak.num(objs["genMus"]) > 1) & (ak.num(objs["dsaMuons"]) > 1))
    ),
    "dsaMuon_absD0_genMus_dR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 200, name="dsaMuon_absD0",
                                     label=r"DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"].d0)),
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(25, 0, 0.03, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   #Added the function ak.ones_like to match delta R array with the d0 array.
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(
                       objs["genMus"][mask, 0])*ak.ones_like(objs["dsaMuons"].d0)),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1
    ),
    "leadingDsaMuon_absD0_genMus_dR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 200, name="dsaMuon_absD0",
                                     label=r"Leading DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"][mask, 0].d0)),
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(25, 0, 0.03, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ((ak.num(objs["genMus"]) > 1) & (ak.num(objs["dsaMuons"]) > 0))
    ),
    "subLeadingDsaMuon_absD0_genMus_dR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 200, name="dsaMuon_absD0",
                                     label=r"Sub Leading DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"][mask, 1].d0)),
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(25, 0, 0.03, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ((ak.num(objs["genMus"]) > 1) & (ak.num(objs["dsaMuons"]) > 1))
    ),
    # dsamuon-genmuon, dR 0.4 window
    "dsaMuon_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="dsaMuon_genMu_ptRatio"),
                   lambda objs, mask: objs["dsaMuons"].pt
                       / objs["dsaMuons"].nearest(objs["genMus"], threshold=0.4).pt),
        ],
    ),
    "dsaMuon0_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="dsaMuon0_genMu_ptRatio"),
                   lambda objs, mask: (objs["dsaMuons"][mask, 0:1].pt
                       / objs["dsaMuons"][mask, 0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ak.num(matched(objs["dsaMuons"] ,objs["genMus"], 0.4)) > 0,
    ),
    "dsaMuon1_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="dsaMuon1_genMu_ptRatio"),
                   lambda objs, mask: (objs["dsaMuons"][mask, 1:2].pt
                       / objs["dsaMuons"][mask, 1:2].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ak.num(matched(objs["dsaMuons"], objs["genMus"], 0.4)) > 1,
    ),
    # pfmuon-genmuon, dR 0.4 window
    "pfMuon_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="pfMuon_genMu_ptRatio"),
                   lambda objs, mask: objs["muons"].pt
                       / objs["muons"].nearest(objs["genMus"], threshold=0.4).pt),
        ],
    ),
    "pfMuon0_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="pfMuon0_genMu_ptRatio"),
                   lambda objs, mask: (objs["muons"][mask,0:1].pt
                       / objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ak.num(matched(objs["muons"], objs["genMus"], 0.4)) > 0,
    ),
    "pfMuon1_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="pfMuon1_genMu_ptRatio"),
                   lambda objs, mask: (objs["muons"][mask,1:2].pt
                       / objs["muons"][mask,1:2].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ak.num(matched(objs["muons"], objs["genMus"], 0.4)) > 1,
    ),
    # gen dark photons (A)
    "genAs_n": obj_attr("genAs", "n"),
    "genAs_toMu_n": obj_attr("genAs_toMu", "n"),
    "genAs_toE_n": obj_attr("genAs_toE", "n"),
    "genAs_pt": obj_attr("genAs", "pt", xmax=200),
    "genAs_pt_highRange": obj_attr("genAs", "pt", xmax=700),
    "egm_lj_pt": obj_attr("genAs", "pt", xmax=700),
    
    "genAs_eta_phi": obj_eta_phi("genAs"),
    "genAs_toMu_matched_muLj_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="genA_n"),
                   lambda objs, mask: ak.num(derived_objs["genAs_toMu_matched_muLj"](objs, 0.4)))
        ],
    ),
    "genAs_toE_matched_egmLj_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="genA_n"),
                   lambda objs, mask: ak.num(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4)))
        ],
    ),
    "genAs_x_y": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0.000, 0.025, name="genAs_x"),
                   lambda objs, mask: objs["genAs"].vx),
            h.Axis(hist.axis.Regular(100, 0.025, 0.050, name="genAs_y"),
                   lambda objs, mask: objs["genAs"].vy),
        ],
    ),
    "genAs_children_x_y": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, -200, 200, name="genAs_children_x"),
                   lambda objs, mask: objs["genAs"].children.vx),
            h.Axis(hist.axis.Regular(100, -200, 200, name="genAs_children_y"),
                   lambda objs, mask: objs["genAs"].children.vy),
        ],
    ),
    "genAs_lxy": obj_attr("genAs", "lxy", xmax=500),
    "genAs_lxy_lowRange": obj_attr("genAs", "lxy", xmax=10),
    "genAs_children_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="genAs_children_n"),
                   lambda objs, mask: ak.num(objs["genAs"].children)),
        ],
    ),
    "genAs_children_absPdgId": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 50, name="genAs_children_absPdgId"),
                   lambda objs, mask: abs(objs["genAs"].children.pdgId)),
        ],
    ),
    "genAs_toMu_lxy": obj_attr("genAs_toMu", "lxy", xmax=500, nbins=100),
    "genAs_toMu_pt": obj_attr("genAs_toMu", "pt", xmax=200, nbins=50),
    "genAs_toMu_pt_highRange": obj_attr("genAs_toMu", "pt", xmax=700, nbins=200),
    "genAs_toMu_eta": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name=r"$Z_d$ $\eta$"),
                   lambda objs, mask: objs["genAs_toMu"].eta ),
        ],
    ),
    "genAs_toE_lxy": obj_attr("genAs_toE", "lxy", xmax=150, nbins=50),
    "genAs_toE_lxy_lowRange": obj_attr("genAs_toE", "lxy", xmax=20),
    "genAs_toE_lxy_highRange": obj_attr("genAs_toE", "lxy", xmax=500),
    "genAs_toE_lxy_midRange": obj_attr("genAs_toE", "lxy", xmin=40, xmax=80),
    "genAs_toE_lxy_ecal": obj_attr("genAs_toE", "lxy", xmin=125, xmax=135),
    "genAs_toE_pt": obj_attr("genAs_toE", "pt", xmax=200, nbins=50),
    "genAs_toE_pt_highRange": obj_attr("genAs_toE", "pt", xmax=700, nbins=200),
    "genAs_toE_eta": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name=r"$Z_d$ $\eta$"),
                   lambda objs, mask: objs["genAs_toE"].eta ),
        ],
    ),
    "genAs_matched_lj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_matched_lj"](objs, 0.4)) ),
        ],
    ),
    "genAs_toMu_matched_lj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toMu_matched_lj"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_lj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_lj"](objs, 0.4)) ),
        ],
    ),
    "genAs_matched_muLj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_matched_muLj"](objs, 0.4)) ),
        ],
    ),
    "genAs_toMu_matched_muLj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toMu_matched_muLj"](objs, 0.4)) ),
        ],
    ),
    "genAs_toMu_matched_muLj_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name=r"$Z_d$ $p_T$ $(GeV)$"),
                   lambda objs, mask: abs(derived_objs["genAs_toMu_matched_muLj"](objs, 0.4).pt) ),
        ],
    ),
    "genAs_toMu_matched_muLj_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(140, 0, 700, name=r"$Z_d$ $p_T$ $(GeV)$"),
                   lambda objs, mask: abs(derived_objs["genAs_toMu_matched_muLj"](objs, 0.4).pt) ),
        ],
    ),
    "genAs_toMu_matched_muLj_eta": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name=r"$Z_d$ $\eta$"),
                   lambda objs, mask: derived_objs["genAs_toMu_matched_muLj"](objs, 0.4).eta ),
        ],
    ),
    "genAs_matched_egmLj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 150, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_matched_egmLj"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_egmLj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 150, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_egmLj_lxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 20, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_egmLj_lxy_midRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 40, 80, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_egmLj_lxy_ecal": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 125, 135, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4)) ),
        ],
    ),
    "genAs_matched_lj_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="genAs_matched_lj_n"),
                   lambda objs, mask: ak.num(derived_objs["genAs_matched_lj"](objs, 0.4)) ),
        ],
    ),
    "genAs_matched_lj_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name=r"$Z_d$ $p_T$ $(GeV)$"),
                   lambda objs, mask: abs(derived_objs["genAs_matched_lj"](objs, 0.4).pt) ),
        ],
    ),
    "genAs_matched_lj_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(140, 0, 700, name=r"$Z_d$ $p_T$ $(GeV)$"),
                   lambda objs, mask: abs(derived_objs["genAs_matched_lj"](objs, 0.4).pt) ),
        ],
    ),
    "genAs_matched_lj_eta": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name=r"$Z_d$ $\eta$"),
                   lambda objs, mask: derived_objs["genAs_matched_lj"](objs, 0.4).eta ),
        ],
    ),
    "genAs_toE_matched_egmLj_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name=r"$Z_d$ $p_T$ $(GeV)$"),
                   lambda objs, mask: abs(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4).pt) ),
        ],
    ),
    "genAs_toE_matched_egmLj_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(140, 0, 700, name=r"$Z_d$ $p_T$ $(GeV)$"),
                   lambda objs, mask: abs(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4).pt) ),
        ],
    ),
    "genAs_toE_matched_egmLj_eta": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name=r"$Z_d$ $\eta$"),
                   lambda objs, mask: derived_objs["genAs_toE_matched_egmLj"](objs, 0.4).eta ),
        ],
    ),
    "genAs_pt_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs"].pt)),
            h.Axis(hist.axis.Regular(250, 0, 500, name="genAs_lxy"),
                   lambda objs, mask: lxy(objs["genAs"])),
        ],
    ),
    "genMu0_pt_MuMudR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 200, name="genMu0_pt",
                                     label=r"Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 0].pt),
            h.Axis(hist.axis.Regular(25, 0, 0.4, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu0_pt_MuMudR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 200, name="genMu0_pt",
                                     label=r"Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 0].pt),
            h.Axis(hist.axis.Regular(25, 0, 0.1, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu0_pt_MuMudR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 200, name="genMu0_pt",
                                     label=r"Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 0].pt),
            h.Axis(hist.axis.Regular(25, 0, 0.04, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu0_pt_highRange_MuMudR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 700, name="genMu0_pt",
                                     label=r"Sub-Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 0].pt),
            h.Axis(hist.axis.Regular(20, 0, 0.25, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu0_pt_highRange_MuMudR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 700, name="genMu0_pt",
                                     label=r"Sub-Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 0].pt),
            h.Axis(hist.axis.Regular(20, 0, 0.06, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu0_pt_highRange_MuMudR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 700, name="genMu0_pt",
                                     label=r"Sub-Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 0].pt),
            h.Axis(hist.axis.Regular(10, 0, 0.01, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu1_pt_MuMudR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 200, name="genMu0_pt",
                                     label=r"Sub-Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 1].pt),
            h.Axis(hist.axis.Regular(25, 0, 0.4, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu1_pt_MuMudR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 200, name="genMu0_pt",
                                     label=r"Sub-Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 1].pt),
            h.Axis(hist.axis.Regular(25, 0, 0.1, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu1_pt_MuMudR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 200, name="genMu0_pt",
                                     label=r"Sub-Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 1].pt),
            h.Axis(hist.axis.Regular(25, 0, 0.04, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu1_pt_highRange_MuMudR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 700, name="genMu0_pt",
                                     label=r"Sub-Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 1].pt),
            h.Axis(hist.axis.Regular(20, 0, 0.25, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu1_pt_highRange_MuMudR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 700, name="genMu0_pt",
                                     label=r"Sub-Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 1].pt),
            h.Axis(hist.axis.Regular(20, 0, 0.06, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu1_pt_highRange_MuMudR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 700, name="genMu0_pt",
                                     label=r"Sub-Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 1].pt),
            h.Axis(hist.axis.Regular(10, 0, 0.01, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genAs_toMu_pt_MuMudR_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,200, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toMu"].pt)),
            h.Axis(hist.axis.Regular(25, 0, 0.4, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genAs_toMu_pt_MuMudR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,200, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toMu"].pt)),
            h.Axis(hist.axis.Regular(25, 0, 0.1, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genAs_toMu_pt_MuMudR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,200, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toMu"].pt)),
            h.Axis(hist.axis.Regular(25, 0, 0.04, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genAs_toMu_lxy_MuMudR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,400, name="genAs_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toMu"])),
            h.Axis(hist.axis.Regular(25, 0, 0.4, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genAs_toMu_lxy_MuMudR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,400, name="genAs_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toMu"])),
            h.Axis(hist.axis.Regular(25, 0, 0.1, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genAs_toMu_lxy_MuMudR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,400, name="genAs_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toMu"])),
            h.Axis(hist.axis.Regular(25, 0, 0.04, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genAs_toMu_lxy_pt_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 400, name="genAs_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toMu"])),
            h.Axis(hist.axis.Regular(25, 0,200, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toMu"].pt)),
        ],
    ),
    "genAs_toE_pt_EEdR_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,200, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toE"].pt)),
            h.Axis(hist.axis.Regular(25, 0, 0.4, name="genE_genE_dR_lowRange"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genAs_toE_pt_EEdR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,200, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toE"].pt)),
            h.Axis(hist.axis.Regular(25, 0, 0.1, name="genE_genE_dR_lowRange"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genAs_toE_pt_EEdR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,200, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toE"].pt)),
            h.Axis(hist.axis.Regular(25, 0, 0.04, name="genE_genE_dR_lowRange"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genAs_toE_lxy_EEdR": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,150, name="genAs_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toE"])),
            h.Axis(hist.axis.Regular(25, 0, 0.4, name="genE_genE_dR_lowRange"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genAs_toE_lxy_EEdR_XLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,150, name="genAs_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toE"])),
            h.Axis(hist.axis.Regular(25, 0, 0.1, name="genE_genE_dR_lowRange"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genAs_toE_lxy_EEdR_XXLowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,150, name="genAs_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toE"])),
            h.Axis(hist.axis.Regular(25, 0, 0.04, name="genE_genE_dR_lowRange"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genAs_toE_lxy_pt_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 150, name="genAs_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toE"])),
            h.Axis(hist.axis.Regular(25, 0,200, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toE"].pt)),
        ],
    ),
    "genAs_toE_pt_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(140, 0,700, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toE"].pt)),
            h.Axis(hist.axis.Regular(50, 0, 200, name="genAs_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toE"])),
        ],
    ),
    "genAs_toE_matched_egmLj_pt_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(140, 0, 700, name="genAs_pt"),
                   lambda objs, mask: abs(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4).pt)),
            h.Axis(hist.axis.Regular(50, 0, 200, name="genAs_lxy"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4))),
        ],
    ),
    "genAs_toMu_pt_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 500, name="genAs_lxy", label =r"$Z_d$ $L_{xy}$"),
                   lambda objs, mask: lxy(objs["genAs_toMu"])),
            h.Axis(hist.axis.Regular(50, 0,700, name="genAs_pt", label =r"$Z_d$ $p_T$"),
                   lambda objs, mask: abs(objs["genAs_toMu"].pt)),
        ],
    ),
    "genAs_toMu_matched_muLj_pt_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 500, name="genAs_lxy"),
                   lambda objs, mask: lxy(derived_objs["genAs_toMu_matched_muLj"](objs, 0.4))),
            h.Axis(hist.axis.Regular(140, 0, 700, name="genAs_pt", label =r"$Z_d$ $p_T$"),
                   lambda objs, mask: abs(derived_objs["genAs_toMu_matched_muLj"](objs, 0.4).pt)),
        ],
    ),
    "genAs_toMu_pt_MuMudR_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,700, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toMu"].pt)),
            h.Axis(hist.axis.Regular(25, 0, 0.3, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genAs_toMu_pt_highRange_MuMudR_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0,700, name="genAs_pt"),
                   lambda objs, mask: abs(objs["genAs_toMu"].pt)),
            h.Axis(hist.axis.Regular(25, 0, 0.04, name="genMu_genMu_dR_lowRange"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    # genA-genA
    "genAs_genAs_absdphi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, math.pi, name="genAs_genAs_absdphi",
                                     label=r"|$\Delta\phi$ between $Z_d$|"),
                   lambda objs, mask: abs(objs["genAs"][mask, 1].delta_phi(objs["genAs"][mask, 0]))),
        ],
        evt_mask=lambda objs: ak.num(objs["genAs"]) > 1,
    ),
    # genA-LJ
    "genAs_lj_dR": h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(200, 0, 2*math.pi, name="genAs_lj_dR"),
                   lambda objs, mask: dR(objs["genAs"], objs["ljs"]))
        ],
    ),
    "genAs_toE_lj_dR": h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(200, 0, 2*math.pi, name="genAs_toE_lj_dR"),
                   lambda objs, mask: dR(objs["genAs_toE"], objs["ljs"]))
        ],
    ),
    
    "genAs_Photon_dR": h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(200, 0, 1.0, name="dR(genAs, RecoPhotons)"),
                   lambda objs, mask: dR(objs["genAs"], objs["photons"]))
        ],
    ),
    # genA - LJ 0.4 matching radius, pT Ratios
    "genA_lj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 4.0, name="genA_lj_ptRatio",
                   label=r"Lepton Jet pT / (closest) $Z_d$ pT"),
                   lambda objs, mask: objs["ljs"].pt
                       / objs["ljs"].nearest(objs["genAs"], threshold=0.4).pt),
        ],
    ),
    "genA_egmLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 4.0, name="genA_egmLj_ptRatio",
                   label=r"EGM Lepton Jet pT / (closest) $Z_d$ pT"),
                   lambda objs, mask: objs["egm_ljs"].pt
                       / objs["egm_ljs"].nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_oneElectronLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_oneElectronLj_ptRatio",
                   label=r"(1) Electron Lepton Jet / (closest) $Z_d$ pT"),
                   lambda objs, mask: derived_objs["n_electron_ljs"](objs, 1).pt
                       / derived_objs["n_electron_ljs"](objs, 1).nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_twoElectronLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_twoElectronLj_ptRatio",
                   label=r"(2) Electron Lepton Jet / (closest) $Z_d$ pT"),
                   lambda objs, mask: derived_objs["n_electron_ljs"](objs, 2).pt
                       / derived_objs["n_electron_ljs"](objs, 2).nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_onePhotonLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_onePhotonLj_ptRatio",
                   label=r"(1) Photon Lepton Jet / (closest) $Z_d$ pT"),
                   lambda objs, mask: derived_objs["n_photon_ljs"](objs, 1).pt
                       / derived_objs["n_photon_ljs"](objs, 1).nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_twoPhotonLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_twoPhotonLj_ptRatio",
                   label=r"(2) Photon Lepton Jet / (closest) $Z_d$ pT"),
                   lambda objs, mask: derived_objs["n_photon_ljs"](objs, 2).pt
                       / derived_objs["n_photon_ljs"](objs, 2).nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_muLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_muLj_ptRatio",
                   label=r"Muon Lepton Jet pT / (closest) $Z_d$ pT"),
                   lambda objs, mask: objs["mu_ljs"].pt
                       / objs["mu_ljs"].nearest(objs["genAs_toMu"], threshold=0.4).pt),
        ],
    ),
    "genA_dsaMuonLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_dsaMuonLj_ptRatio",
                   label=r"DSA Muon Lepton Jet pT / (closest) $Z_d$ pT"),
                   lambda objs, mask: (objs["dsaMuons"][mask]).nearest(objs["ljs"][mask], threshold=0.4).pt
                       / (objs["dsaMuons"][mask]).nearest(objs["genAs_toMu"][mask], threshold=0.4).pt),
        ],
    ),
    "genA_pfMuonLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_pfMuonLj_ptRatio",
                   label=r"PF Muon Lepton Jet pT / (closest) $Z_d$ pT"),
                   lambda objs, mask: (objs["muons"][mask]).nearest(objs["ljs"][mask], threshold=0.4).pt
                       / (objs["muons"][mask]).nearest(objs["genAs_toMu"][mask], threshold=0.4).pt),
        ],
    ),
    "genA_dsaMuon0Lj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_dsaMuonLj_ptRatio",
                   label=r"Lead DSA Muon Lepton Jet / (closest) $Z_d$ pT"),
                   lambda objs, mask: ((objs["dsaMuons"][mask, 0:1]).nearest(objs["ljs"][mask], threshold=0.4).pt
                       / (objs["dsaMuons"][mask, 0:1]).nearest(objs["genAs_toMu"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["dsaMuons"], objs["genAs_toMu"], 0.4)) > 0)
                               & (ak.num(matched(objs["dsaMuons"], objs["ljs"], 0.4)) > 0)),
    ),
    "genA_pfMuon0Lj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_pfMuonLj_ptRatio",
                   label=r"Lead PF Muon Lepton Jet / (closest) $Z_d$ pT"),
                   lambda objs, mask: ((objs["muons"][mask, 0:1]).nearest(objs["ljs"][mask], threshold=0.4).pt
                       / (objs["muons"][mask, 0:1]).nearest(objs["genAs_toMu"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["muons"], objs["genAs_toMu"], 0.4)) > 0)
                               & (ak.num(matched(objs["muons"], objs["ljs"], 0.4)) > 0)),
    ),
    # genA - LJ 0.4 matching radius, LJ Reco Lxy / True Lxy
    "genA_muLj_lxyRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_muLj_lxyRatio",
                                    label=r"Muon Lepton Jet Reco L$_{xy}$ / (closest) $Z_d$ L$_{xy}$"),
                   lambda objs, mask: objs["mu_ljs"].kinvtx.lxy
                       / lxy(objs["mu_ljs"].nearest(objs["genAs"], threshold=0.4))),
        ],
    ),
    "genA_egmLj_lxyRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_egmLj_lxyRatio",
                                    label=r"EGM Lepton Jet Reco L$_{xy}$ / (closest) $Z_d$ L$_{xy}$"),
                   lambda objs, mask: objs["egm_ljs"].kinvtx.lxy
                       / lxy(objs["egm_ljs"].nearest(objs["genAs"], threshold=0.4))),
        ],
    ),
    # LJ Res vs Reco Lxy
    "mu_lj_genA_ptRatio_vs_recolxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="mu_lj_genA_ptRatio"),
                   lambda objs, mask: objs["mu_ljs"].pt
                       / objs["mu_ljs"].nearest(objs["genAs"]).pt),
            h.Axis(hist.axis.Regular(100, 0, 300, name="mu_lj_recolxy"),
                   lambda objs, mask: objs["mu_ljs"].kinvtx.lxy),
        ],
    ),
    "egm_lj_genA_ptRatio_vs_recolxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="egm_lj_genA_ptRatio"),
                   lambda objs, mask: objs["egm_ljs"].pt
                       / objs["egm_ljs"].nearest(objs["genAs"]).pt),
            h.Axis(hist.axis.Regular(100, 0, 300, name="egm_lj_recolxy"),
                   lambda objs, mask: objs["egm_ljs"].kinvtx.lxy),
        ],
    ),
    # LJ Res vs True Lxy, 0.4 thresholds on dR matching
    "egm_lj_genA_ptRatio_vs_truelxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="egm_lj_genA_ptRatio"),
                   lambda objs, mask: objs["egm_ljs"].pt
                       / objs["egm_ljs"].nearest(objs["genAs"], threshold=0.4).pt),
            h.Axis(hist.axis.Regular(100, 0, 300, name="egm_lj_truelxy"),
                   lambda objs, mask: lxy(objs["egm_ljs"].nearest(objs["genAs"], threshold=0.4))),
        ],
    ),
    "mu_lj_genA_ptRatio_vs_truelxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="mu_lj_genA_ptRatio"),
                   lambda objs, mask: objs["mu_ljs"].pt
                       / objs["mu_ljs"].nearest(objs["genAs"], threshold=0.4).pt),
            h.Axis(hist.axis.Regular(100, 0, 300, name="mu_lj_truelxy"),
                   lambda objs, mask: lxy(objs["mu_ljs"].nearest(objs["genAs"], threshold=0.4))),
        ],
    ),
    "dsaMuon0_genMu0_ptRatio_vs_truelxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0., 2.0, name="dsaMuon0_genMu0_ptRatio"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].pt
                       / objs["dsaMuons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(100, 0, 300, name="dsaMuon0_lj_truelxy"),
                   lambda objs, mask: lxy(objs["dsaMuons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4))),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["genMus"], objs["dsaMuons"], 0.4)) > 0)
                               & (ak.num(matched(objs["genAs"], objs["dsaMuons"], 0.4)) > 0)),
    ),
    "muon0_genMu0_ptRatio_vs_truelxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0., 2.0, name="muon0_genMu0_ptRatio"),
                   lambda objs, mask: (objs["muons"][mask,0:1].pt
                       / objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(100, 0, 300, name="pfMuon0_lj_truelxy"),
                   lambda objs, mask: lxy(objs["muons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4))),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["genMus"], objs["muons"], 0.4)) > 0)
                               & (ak.num(matched(objs["genAs"], objs["muons"], 0.4)) > 0)),
    ),
    # LJ Res vs True pT, dR 0.4 matching window
    "dsaMuon0_genMu0_ptRatio_vs_truept": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2.0, name="dsaMuon0_genMu0_ptRatio"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].pt
                       / objs["dsaMuons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: (ak.num(objs["dsaMuons"]) > 0),
    ),
    "muon0_genMu0_ptRatio_vs_truept": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2.0, name="muon0_genMu0_ptRatio"),
                   lambda objs, mask: (objs["muons"][mask,0:1].pt
                       / objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: (ak.num(objs["muons"]) > 0),
    ),
    "dsaMuon0_muLj_ptRatio_vs_truept": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2.0, name="dsaMuon0_genMu0_ptRatio"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].nearest(objs["ljs"][mask], threshold=0.4).pt
                       / objs["dsaMuons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: (ak.num(objs["dsaMuons"]) > 0),
    ),
    "muon0_muLj_ptRatio_vs_truept": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0., 2.0, name="muon0_muLj_ptRatio"),
                   lambda objs, mask: (objs["muons"][mask,0:1].nearest(objs["ljs"][mask], threshold=0.4).pt
                       / objs["muons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["ljs"], objs["muons"], 0.4)) > 0)
                               & (ak.num(matched(objs["genMus"], objs["muons"], 0.4)) > 0)),
    ),
    "egmLj_ptRatio_vs_egm_truept": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0., 2.0, name="egm_lj_genA_ptRatio"),
                   lambda objs, mask: objs["egm_ljs"].pt
                       / objs["egm_ljs"].nearest(objs["genAs"], threshold=0.4).pt),
            h.Axis(hist.axis.Regular(100, 0, 1000, name="genE_pt"),
                   lambda objs, mask: (objs["egm_ljs"].nearest(objs["genEs"], threshold=0.4).pt)[mask,0]),
        ],
    ),
    # LJ True pT vs True Lxy, dR 0.4 matching window
    "genMu0_truept_vs_dsaMuon0_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(100, 0, 300, name="dsaMuon0_lj_truelxy"),
                   lambda objs, mask: lxy(objs["dsaMuons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4))),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["genMus"], objs["dsaMuons"], 0.4)) > 0)
                               & (ak.num(matched(objs["genAs"], objs["dsaMuons"], 0.4)) > 0)),
    ),
    "genMu0_truept_vs_muon0_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(100, 0, 300, name="pfMuon0_lj_truelxy"),
                   lambda objs, mask: lxy(objs["muons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4))),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["genMus"], objs["muons"], 0.4)) > 0)
                               & (ak.num(matched(objs["genAs"], objs["muons"], 0.4)) > 0)),
    ),
    "genAs_toE_matched_one_photon_ljs_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 150, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_one_photon_ljs"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_two_photon_ljs_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 150, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_two_photon_ljs"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_one_electron_ljs_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 150, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_one_electron_ljs"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_two_electron_ljs_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 150, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_two_electron_ljs"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_one_e_one_p_ljs_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 150, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs['genAs_toE_matched_one_e_one_p_ljs'](objs, 0.4)) ),
        ],
    ),
    ###
    "genAs_toE_matched_one_photon_ljs_lxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 5, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_one_photon_ljs"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_two_photon_ljs_lxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 5, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_two_photon_ljs"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_one_electron_ljs_lxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 5, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_one_electron_ljs"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_two_electron_ljs_lxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 5, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_two_electron_ljs"](objs, 0.4)) ),
        ],
    ),
    "genAs_toE_matched_one_e_one_p_ljs_lxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 5, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs['genAs_toE_matched_one_e_one_p_ljs'](objs, 0.4)) ),
        ],
    ),
    "genA_egmLj_ptRatio_onePhoton": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 4.0, name="genA_egmLj_ptRatio_onePhoton",
                   label=r"EGM (1 Photon) Lepton Jet pT / (closest) $Z_d$ pT"),
                   lambda objs, mask: objs["one_photon_ljs"].pt / objs["one_photon_ljs"].nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_egmLj_ptRatio_twoPhoton": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 4.0, name="genA_egmLj_ptRatio_twoPhoton",
                   label=r"EGM (2 Photon) Lepton Jet pT / (closest) $Z_d$ pT"),
                   lambda objs, mask: objs["two_photon_ljs"].pt / objs["two_photon_ljs"].nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_egmLj_ptRatio_oneElectron": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 4.0, name="genA_egmLj_ptRatio_oneElectron",
                   label=r"EGM (1 Electron) Lepton Jet pT / (closest) $Z_d$ pT"),
                   lambda objs, mask: objs["one_electron_ljs"].pt / objs["one_electron_ljs"].nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_egmLj_ptRatio_twoElectron": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 4.0, name="genA_egmLj_ptRatio_twoElectron",
                   label=r"EGM (2 Electron) Lepton Jet pT / (closest) $Z_d$ pT"),
                   lambda objs, mask: objs["two_electron_ljs"].pt / objs["two_electron_ljs"].nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_egmLj_ptRatio_oneEoneP": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 4.0, name="genA_egmLj_ptRatio_oneEoneP",
                   label=r"EGM (1E 1P) Lepton Jet pT / (closest) $Z_d$ pT"),
                   lambda objs, mask: objs["one_e_one_p_ljs"].pt / objs["one_e_one_p_ljs"].nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),

    "photon_electron_dR_from_oneEonePLJ": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(30, 0, 0.2, name="photon_electron_dR_from_oneEonePLJ",
                                     label=r"$\Delta R$($p$, $e$)"),
                   lambda objs, mask: dR(objs["one_e_one_p_ljs"].photons, objs["one_e_one_p_ljs"].electrons))
        ],
    ),
    
    "photon_electron_dR_from_oneEonePLJ_lowThreshold": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(50, 0, 0.2, name="photon_electron_dR_from_oneEonePLJ",
                                     label=r"$\Delta R$($p$, $e$)"),
                   lambda objs, mask: dR(objs["one_e_one_p_ljs"][derived_objs['genA_egmLj_oneEoneP_ptRatio_PS'](objs) < 1.2].photons, objs["one_e_one_p_ljs"][derived_objs['genA_egmLj_oneEoneP_ptRatio_PS'](objs) < 1.2].electrons))
        ],
    ),
    
    "photon_electron_dR_from_oneEonePLJ_highThreshold": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(50, 0, 0.2, name="photon_electron_dR_from_oneEonePLJ",
                                     label=r"$\Delta R$($p$, $e$)"),
                  lambda objs, mask: dR(objs["one_e_one_p_ljs"][derived_objs['genA_egmLj_oneEoneP_ptRatio_PS'](objs) > 1.8].photons, objs["one_e_one_p_ljs"][derived_objs['genA_egmLj_oneEoneP_ptRatio_PS'](objs) > 1.8].electrons))
        ],
    ),
    'genA_toE_matched_photons_lxyLowRange': h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 20, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs['genAs_toE_matched_photons'](objs,0.4)))
        ],
    ),
     'genA_toE_matched_photons_lxy': h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 150, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs['genAs_toE_matched_photons'](objs,0.4)))
        ],
    ),
    'genA_toE_matched_electrons_lxyLowRange': h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 20, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs['genAs_toE_matched_electrons'](objs,0.4)))
        ],
    ),
     'genA_toE_matched_electrons_lxy': h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 150, name=r"$Z_d$ $L_{xy}$ $(cm)$"),
                   lambda objs, mask: lxy(derived_objs['genAs_toE_matched_electrons'](objs,0.4)))
        ],
    ),
    "mJJ_2mu2e": h.Histogram(
        [
            h.Axis(hist.axis.Regular(30, 0, 1200.0, name=r"M_{jj}",
                   label=r"$M_{jj}$ [GeV]"),
                   lambda objs, mask: (objs["egm_ljs"][mask, 0] + objs["mu_ljs"][mask, 0]).mass)
        ],
        evt_mask=lambda objs: (ak.num(objs["egm_ljs"].mass, axis=1) > 0) & (ak.num(objs["mu_ljs"].mass, axis=1) > 0),
    ),

}
