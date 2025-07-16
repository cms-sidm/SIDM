"""Module to define the base SIDM processor"""

# python
import copy
import numpy as np
# columnar analysis
from coffea import processor
from coffea.nanoevents.methods import nanoaod
from coffea.nanoevents.methods import vector as cvec
import awkward as ak
import fastjet
import vector
#local
from sidm import BASE_DIR
from sidm.tools import selection, cutflow, utilities
from sidm.definitions.hists import hist_defs, counter_defs
from sidm.definitions.objects import preLj_objs, postLj_objs

class SidmProcessor(processor.ProcessorABC):
    """Class to apply selections, make histograms, and make cutflows

    Accepts NanoEvents records that are assumed to have been produced by FFSchema. Selections are
    chosen by supplying a list of selection names (as defined in selections.yaml), and histograms
    are chosen by providing a list of histogram collection names (as definined in
    hist_collections.yaml).
    """

    def __init__(
        self,
        channel_names,
        hist_collection_names,
        lj_reco_choices=["0.4"],
        selections_cfg="configs/selections.yaml",
        histograms_cfg="configs/hist_collections.yaml",
        unweighted_hist=False,
        verbose=False,
    ):
        self.channel_names = channel_names
        self.hist_collection_names = hist_collection_names
        self.lj_reco_choices = lj_reco_choices
        self.selections_cfg = selections_cfg
        self.histograms_cfg = histograms_cfg
        self.unweighted_hist = unweighted_hist
        self.obj_defs = preLj_objs
        self.verbose = verbose
        self.year = "2018" # fixme: may be better to store as event metadata

    def process(self, events):
        """Apply selections, make histograms and cutflow"""
        # create object collections
        # fixme: only include objs used in cuts or hists
        objs = {}
        for obj_name, obj_def in self.obj_defs.items():
            try:
                obj = obj_def(events)
            except AttributeError:
                print(f"Warning: {obj_name} not found in this sample. Skipping.")
                continue
            objs[obj_name] = obj

            # pt order
            objs[obj_name] = self.order(objs[obj_name])

            # use nanoevents.Muon behaviors for dsa muons
            if obj_name == "dsaMuons":
                forms = {f: objs[obj_name][f] for f in objs[obj_name].fields}
                objs[obj_name] = ak.zip(forms, with_name="Muon", behavior=nanoaod.behavior)

            # add lxy attribute to particles with children
            if hasattr(obj, "children"):
                objs[obj_name]["lxy"] = utilities.lxy(objs[obj_name])

            # add dxy wrt beamspot for all objs that don't already have it
            if hasattr(obj, "vx") and not hasattr(obj, "dxy") and "bs" in objs:
                objs[obj_name]["dxy"] = utilities.dxy(objs[obj_name], ref=objs["bs"])

            # add dimension to one-per-event objects to allow independent obj and evt cuts
            # skip objects with no fields
            if objs[obj_name].ndim == 1 and "x" in obj.fields:
                counts = ak.ones_like(objs[obj_name].x, dtype=np.int32)
                objs[obj_name] = ak.unflatten(objs[obj_name], counts)

        cutflows = {}
        counters = {}

        # define histograms
        hists = self.build_histograms()

        ### define pre-lj object, lj, post-lj obj, and event cuts per channel
        ch_cuts = self.build_cuts()

        # loop through lj reco choices and channels, treating each lj+channel pair as a unique Selection
        for channel, cuts in ch_cuts.items():
            obj_selection = selection.JaggedSelection(cuts["obj"], self.verbose)

            for lj_reco in self.lj_reco_choices:
                # apply pre-LJ object selection
                sel_objs = obj_selection.apply_obj_cuts(objs)

                # reconstruct lepton jets
                sel_objs["ljs"] = self.build_lepton_jets(sel_objs, float(lj_reco))

                # apply obj selection to ljs
                lj_selection = selection.JaggedSelection(cuts["lj"], self.verbose)
                sel_objs = lj_selection.apply_obj_cuts(sel_objs)

                # add post-lj objects to sel_objs
                for obj in postLj_objs:
                    sel_objs[obj] = postLj_objs[obj](sel_objs)

                # apply post-lj obj selection
                postLj_selection = selection.JaggedSelection(cuts["postLj_obj"], self.verbose)
                sel_objs = postLj_selection.apply_obj_cuts(sel_objs)

                # build Selection objects and apply event selection
                evt_selection = selection.Selection(cuts["evt"], self.verbose)
                sel_objs = evt_selection.apply_evt_cuts(sel_objs)

                # fill all hists
                sel_objs["ch"] = channel
                sel_objs["lj_reco"] = lj_reco

                # define event weights
                evt_weights =  self.obj_defs["weight"](events)*events.metadata["skim_factor"]

                # make cutflow
                if lj_reco not in cutflows:
                    cutflows[str(lj_reco)] = {}
                cutflows[str(lj_reco)][channel] = cutflow.Cutflow(evt_selection.all_evt_cuts, evt_selection.evt_cuts, evt_weights)

                # fill histograms for this channel+lj_reco pair
                hist_weights = evt_weights[evt_selection.all_evt_cuts.all(*evt_selection.evt_cuts)]
                if self.unweighted_hist:
                    hist_weights =  ak.ones_like(hist_weights)
                for h in hists.values():
                    h.fill(sel_objs, hist_weights)

                # Fill counters
                if lj_reco not in counters:
                    counters[lj_reco] = {}
                counters[lj_reco][channel] = {}

                for name, counter in counter_defs.items():
                    try:
                        counters[lj_reco][channel][name] = counter(sel_objs)
                    except (KeyError, AttributeError) as e:
                        print(f"Warning: cannot fill counter {name}. Skipping.")

        # lose lj_reco dimension to cutflows if only one reco was run
        if len(self.lj_reco_choices) == 1:
            cutflows = cutflows[self.lj_reco_choices[0]]

        out = {
            "cutflow": cutflows,
            "hists": {n: h.hist for n, h in hists.items()}, # output hist.Hists, not Histograms
            "counters": counters,
            "metadata": {
                "n_evts": events.metadata["entrystop"] - events.metadata["entrystart"],
            },
        }

        return {events.metadata["dataset"]: out}

    def make_vector(self, objs, collection, fields, type_id=None, mass=None):
        shape = ak.ones_like(objs[collection].pt)
        # all objects must have the same fields to later concatenate and cluster them
        # set fields that aren't available for a given object to be None
        # these additional fields will be removed after clustering
        nan = ak.full_like(shape, None)
        forms = {f: objs[collection][f] if f in objs[collection].fields else nan for f in fields}
        forms["part_type"] = objs[collection]["type"] if type_id is None else type_id*shape
        forms["mass"] = objs[collection]["mass"] if mass is None else mass*shape
        return vector.zip(forms)

    def make_constituent(self, consts, type_ids, name, fields):
        """Return array of particles of given type_ids, name, and only specified fields"""
        relevant_consts = consts[ak.any((consts.part_type == x for x in type_ids), axis=0)]
        forms = {f: relevant_consts.__getattr__(f) for f in fields}
        return ak.zip(forms, with_name=name, behavior=nanoaod.behavior)

    def build_lepton_jets(self, objs, lj_reco):
        """Reconstruct lepton jets according to defintion given by lj_reco"""

        # Use electron/muon/photon/dsamuon collections with a custom distance parameter
        collections = ["muons", "dsaMuons", "electrons", "photons"]
        fields = [objs[c].fields for c in collections]
        all_fields = list(set().union(*fields))
        muon_inputs = self.make_vector(objs, "muons", all_fields,  type_id=3)
        dsa_inputs = self.make_vector(objs, "dsaMuons", all_fields, type_id=8, mass=0.106)
        ele_inputs = self.make_vector(objs, "electrons", all_fields, type_id=2)
        photon_inputs = self.make_vector(objs, "photons", all_fields, type_id=4)
        lj_inputs = ak.concatenate([muon_inputs, dsa_inputs, ele_inputs, photon_inputs], axis=-1)

        distance_param = abs(lj_reco)
        jet_def = fastjet.JetDefinition(fastjet.antikt_algorithm, distance_param)
        cluster = fastjet.ClusterSequence(lj_inputs, jet_def)
        jets = cluster.inclusive_jets()

        # turn lepton jets back into LorentzVectors that match existing structures
        ljs = ak.zip(
            {"x": jets.x,
             "y": jets.y,
             "z": jets.z,
             "t": jets.t},
            with_name="LorentzVector",
            behavior=nanoaod.behavior
        )

        # add fields to access LJ constituents
        consts = cluster.constituents()
        common_fields = list(set(fields[0]).intersection(*fields[1:]))
        muon_fields = list(set(objs["muons"].fields).intersection(objs["dsaMuons"].fields))
        ljs["constituents"] = self.make_constituent(consts, [2, 3, 4, 8], "PtEtaPhiMCollection", common_fields)
        ljs["muons"] = self.make_constituent(consts, [3, 8], "Muon", muon_fields)
        ljs["pfMuons"] = self.make_constituent(consts, [3], "Muon", objs["muons"].fields)
        ljs["dsaMuons"] = self.make_constituent(consts, [8], "Muon", objs["dsaMuons"].fields)
        ljs["electrons"] = self.make_constituent(consts, [2], "Electron", objs["electrons"].fields)
        ljs["photons"] = self.make_constituent(consts, [4], "Photon", objs["photons"].fields)

        # define LJ-level quantities

        # number of constituents
        ljs["pfMu_n"] = ak.num(ljs.pfMuons, axis=-1)
        ljs["dsaMu_n"] = ak.num(ljs.dsaMuons, axis=-1)
        ljs["muon_n"] = ak.num(ljs.muons, axis=-1)
        ljs["electron_n"] = ak.num(ljs.electrons, axis=-1)
        ljs["photon_n"] = ak.num(ljs.photons, axis=-1)

        # dRSpread (the maximum dR betwen any pair of constituents in each lepton jet)
        # a) for each constituent, find the dR between it and all other constituents in the same LJ
        # b) flatten that into a list of dRs per LJ
        # c) and then take the maximum dR per LJ, leaving us with a single value per LJ
        ljs["dRSpread"] = ak.max(ak.flatten(
            ljs["constituents"].metric_table(ljs["constituents"], axis=2), axis=-1), axis=-1)

        # LJ isolation
        ljs["matched_jet"] = ljs.nearest(objs["jets"], threshold=0.4)
        ljs["lepton_fraction"] =  ljs["matched_jet"].chEmEF + ljs["matched_jet"].neEmEF + ljs["matched_jet"].muEF
        ljs["isolation"] = ak.fill_none((ljs["matched_jet"].energy / ljs.energy) * (1 - (ljs["lepton_fraction"])), 0)
        ljs["dR_matched_jet"] = ljs.delta_r(ljs["matched_jet"])
        
        # todo: add LJ displacement

        # pt order the new LJs
        ljs = self.order(ljs)

        # return the new LJ collection
        return ljs

    def build_cuts(self):
        """ Make list of pre-lj object, lj, post-lj obj, and event cuts per channel"""

        selection_menu = utilities.load_yaml(f"{BASE_DIR}/{self.selections_cfg}")

        ch_cuts = {}

        for channel in self.channel_names:
            ch_cuts[channel] = {}
            ch_cuts[channel]["obj"] = {}
            ch_cuts[channel]["lj"] = {}
            ch_cuts[channel]["postLj_obj"] = {}
            ch_cuts[channel]["evt"] = {}

            cuts = selection_menu[channel]
            for obj, obj_cuts in cuts["obj_cuts"].items():
                if obj not in ch_cuts[channel]["obj"]:
                    ch_cuts[channel]["obj"][obj] = []
                ch_cuts[channel]["obj"][obj] = utilities.flatten(obj_cuts)

            if "postLj_obj_cuts" in cuts:
                for obj, obj_cuts in cuts["postLj_obj_cuts"].items():
                    if obj == "ljs":
                        ch_cuts[channel]["lj"][obj] = utilities.flatten(obj_cuts)
                    else:
                        ch_cuts[channel]["postLj_obj"][obj] = utilities.flatten(obj_cuts)

            if "evt_cuts" in cuts:
                ch_cuts[channel]["evt"] = utilities.flatten(cuts["evt_cuts"])

        return ch_cuts

    def build_histograms(self):
        """Create dictionary of Histogram objects"""
        hist_menu = utilities.load_yaml(f"{BASE_DIR}/{self.histograms_cfg}")
        # build dictionary and create hist.Hist objects
        hists = {}
        for collection in self.hist_collection_names:
            collection = utilities.flatten(hist_menu[collection])
            for hist_name in collection:
                hists[hist_name] = copy.deepcopy(hist_defs[hist_name])
                # Add lj_reco axis only when more than one reco is run
                lj_reco_names = self.lj_reco_choices if len(self.lj_reco_choices) > 1 else None
                hists[hist_name].make_hist(hist_name, self.channel_names, lj_reco_names)
        return hists

    def order(self, obj):
        """Explicitly order objects"""
        # pt order objects with a pt attribute
        if hasattr(obj, "pt"):
            obj = obj[ak.argsort(obj.pt, ascending=False)]
        # fixme: would be good to explicitly order other objects as well
        return obj

    def postprocess(self, accumulator):
        """Modify accumulator after process has run on all chunks"""
        # scale cutflow and hists according to lumi*xs
        for sample, output in accumulator.items():
            n_evts = output["metadata"]["n_evts"]
            lumixs_weight = utilities.get_lumixs_weight(sample, self.year, n_evts)
            for name in output["cutflow"]:
                if len(self.lj_reco_choices) > 1:
                    for key in output['cutflow'][name]:
                         accumulator[sample]["cutflow"][name][key].scale(lumixs_weight)
                else:
                    accumulator[sample]["cutflow"][name].scale(lumixs_weight)
            if not self.unweighted_hist:
                for name in output["hists"]:
                    accumulator[sample]["hists"][name] *= lumixs_weight
