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
from sidm.definitions.objects import preLj_objs, postLj_objs, postLj_objs_MC
import coffea.nanoevents.transforms as tr

def _patched_local2global(stack):
    """
    Original: index,target_offsets,!local2global
    Turn jagged local index into global index
    """
    target_offsets = ak.Array(stack.pop())
    index = ak.Array(stack.pop())
    index = index.mask[index >= 0] + target_offsets[:-1]
    index = index.mask[index < target_offsets[1:]]

    out = ak.flatten(ak.fill_none(index, -1), axis=None)
    out = ak.values_astype(out, np.int64)

    stack.append(out)
tr.local2global = _patched_local2global

class _AwkwardColumnAccumulator(processor.AccumulatorABC):
    """Accumulator that concatenates per-chunk awkward arrays.

    Used only for the optional debug output. Each chunk contributes one awkward
    array; chunks are stored by reference and concatenated lazily when `.value`
    is accessed, which keeps accumulation memory-efficient and avoids converting
    columns to Python lists.
    """

    def __init__(self, chunks=None):
        self._chunks = list(chunks) if chunks is not None else []

    def identity(self):
        return _AwkwardColumnAccumulator()

    def add(self, other):
        self._chunks.extend(other._chunks)
        return self

    @property
    def value(self):
        if not self._chunks:
            return ak.Array([])
        if len(self._chunks) == 1:
            return self._chunks[0]
        return ak.concatenate(self._chunks)

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
        debug=False,
        debug_branches=None,
        include_default_debug_branches=True,
        debug_suppress_failures=True,
    ):
        self.channel_names = channel_names
        self.hist_collection_names = hist_collection_names
        self.lj_reco_choices = lj_reco_choices
        self.selections_cfg = selections_cfg
        self.histograms_cfg = histograms_cfg
        self.unweighted_hist = unweighted_hist
        self.obj_defs = preLj_objs
        self.postLj_objs = postLj_objs
        self.postLj_objs_MC = postLj_objs_MC
        self.verbose = verbose

        # Optional debug output. Disabled by default, no effect on the standard
        # output. When debug=True, an additional out["debug"] dict is added:
        #
        #     out["debug"][lj_reco][channel][branch_name]  -> _AwkwardColumnAccumulator
        #     read with: output[dataset]["debug"]["0.4"][channel]["name"].value
        #
        # apply_evt_cuts already filters sel_objs to the events passing the event
        # selection, so branches that read sel_objs are over selected events;
        # branches that read raw `events` (e.g. gen weights) stay full-chunk.
        #
        # Supply analysis-specific branches via debug_branches, e.g.:
        #     SidmProcessor(..., debug=True, debug_branches=sidm_debug_branches(),
        #                   include_default_debug_branches=False)
        self.debug = debug
        self.debug_suppress_failures = debug_suppress_failures

        self.debug_branches = {}
        if include_default_debug_branches:
            self.debug_branches.update(self.default_debug_branches())
        if debug_branches is not None:
            self.debug_branches.update(debug_branches)

    def process(self, events):
        """Apply selections, make histograms and cutflow"""
        is_data = events.metadata["is_data"]
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
        debug_output = {} if self.debug else None

        # define histograms
        hists = self.build_histograms()

        # define pre-lj object, lj, post-lj obj, and event cuts per channel
        ch_cuts = self.build_cuts()

        # define event weights
        if not is_data:
            evt_weights = self.obj_defs["weight"](events)
        else:
            evt_weights = ak.broadcast_arrays(1.0, self.obj_defs["met"](events))[0]

        # loop through lj reco choices and channels, treating each lj+channel pair as a unique Selection
        for channel, cuts in ch_cuts.items():
            obj_selection = selection.JaggedSelection(cuts["obj"], self.verbose)
            nested_selection = selection.NestedSelection(cuts["obj"], self.verbose)

            for lj_reco in self.lj_reco_choices:
                sel_objs = objs.copy()

                # apply selections on matched_muons within the DSA muons and matched_dsa_muons within the PF muons
                # remove None entries from matched PF or DSA muons before applying cuts
                sel_objs["dsaMuons"]["good_matched_muons"] = nested_selection.apply_obj_cuts(sel_objs, ak.drop_none(sel_objs["dsaMuons"].matched_muons), "muons")
                sel_objs["muons"]["good_matched_dsa_muons"] = nested_selection.apply_obj_cuts(sel_objs, ak.drop_none(sel_objs["muons"].matched_dsa_muons), "dsaMuons")

                # apply pre-LJ object selection
                sel_objs = obj_selection.apply_obj_cuts(sel_objs)

                # reconstruct lepton jets
                sel_objs["ljs"] = self.build_lepton_jets(sel_objs, float(lj_reco))

                # apply obj selection to ljs
                lj_selection = selection.JaggedSelection(cuts["lj"], self.verbose)
                sel_objs = lj_selection.apply_obj_cuts(sel_objs)

                # add post-lj objects to sel_objs
                if not is_data:
                    self.postLj_objs = {**self.postLj_objs, **self.postLj_objs_MC}
                for obj in self.postLj_objs:
                    sel_objs[obj] = self.postLj_objs[obj](sel_objs)

                # apply post-lj obj selection
                postLj_selection = selection.JaggedSelection(cuts["postLj_obj"], self.verbose)
                sel_objs = postLj_selection.apply_obj_cuts(sel_objs)
 
                # build Selection objects and apply event selection
                sel_objs["evt_weights"] = evt_weights
                evt_selection = selection.Selection(cuts["evt"], self.verbose)
                sel_objs = evt_selection.apply_evt_cuts(sel_objs)

                # optional debug output (no effect unless debug=True).
                # apply_evt_cuts has already trimmed sel_objs to events passing the
                # event selection, so branches that read sel_objs are over selected
                # events; branches that read raw `events` stay full-chunk.
                if self.debug:
                    lj_reco_key = str(lj_reco)
                    if lj_reco_key not in debug_output:
                        debug_output[lj_reco_key] = {}
                    debug_output[lj_reco_key][channel] = self.fill_debug_branches(
                        sel_objs, events
                    )

                # fill all hists

                # fixme: disable cutflows due to sequential event cut implementation
                # store cutflow in separate dict
                if lj_reco not in cutflows:
                    cutflows[str(lj_reco)] = {}
                cutflows[str(lj_reco)][channel] = evt_selection.cutflow

                # fill histograms for this channel+lj_reco pair
                sel_objs["ch"] = channel
                sel_objs["lj_reco"] = lj_reco
                hist_weights = sel_objs["evt_weights"]
                if self.unweighted_hist:
                    hist_weights =  ak.ones_like(hist_weights)
                for h in hists.values():
                    h.fill(sel_objs, hist_weights, self.verbose)

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
        # fixme: disable cutflows due to sequential event cut implemention
        if len(self.lj_reco_choices) == 1:
            cutflows = cutflows[self.lj_reco_choices[0]]

        out = {
            "cutflow": cutflows,
            "hists": {n: h.hist for n, h in hists.items()}, # output hist.Hists, not Histograms
            "counters": counters,
            "metadata": {
                "n_evts": events.metadata["entrystop"] - events.metadata["entrystart"],
                "scaled_sum_weights": ak.sum(evt_weights)/events.metadata["skim_factor"],
                # add sample metadata as set_accumulator to only keep unique values during accumulation
                "year": processor.set_accumulator([events.metadata["year"]]),
                "is_data": processor.set_accumulator([events.metadata["is_data"]]),
                "unweighted_hist": processor.set_accumulator([self.unweighted_hist]),
            },
        }

        # optional debug output: keep the standard return shape, add one key
        if self.debug:
            out["debug"] = debug_output

        return {events.metadata["dataset"]: out}

    @staticmethod
    def default_debug_branches():
        """Minimal, analysis-agnostic debug branches.

        passing_weights reads sel_objs, so it is over selected events.
        gen_weights reads raw events, so it stays full-chunk (useful for scaling);
        it is skipped on data (no Generator) when debug_suppress_failures=True.
        Anything analysis-specific should be supplied via debug_branches.
        """
        return {
            "passing_weights": lambda sel_objs, events: sel_objs["evt_weights"],
            "gen_weights": lambda sel_objs, events: (
                events.Generator.weight if hasattr(events, "Generator") else ak.Array([])
            ),
        }

    def fill_debug_branches(self, sel_objs, events):
        """Evaluate all registered debug branches for one channel + lj_reco pair.

        apply_evt_cuts has already trimmed sel_objs to the events passing the full
        event selection, so branches reading sel_objs are automatically over
        selected events; branches reading raw `events` (e.g. generator weights)
        remain full-chunk. Each branch is wrapped in its own try/except so one
        failing branch does not abort the job (unless debug_suppress_failures is
        False). Results are stored as awkward arrays in an accumulator that
        concatenates across chunks.
        """
        debug = {}
        for name, branch_func in self.debug_branches.items():
            try:
                value = self._to_debug_array(branch_func(sel_objs, events))
                debug[name] = _AwkwardColumnAccumulator([value])
            except Exception as e:
                if not self.debug_suppress_failures:
                    raise
                print(f"Warning: cannot fill debug branch {name}. Skipping. Error: {e}")
                debug[name] = _AwkwardColumnAccumulator()
        return debug

    @staticmethod
    def _to_debug_array(value):
        """Coerce a branch result into a plain, picklable awkward array.

        NanoAOD / `vector` LorentzVector behaviors are implemented with lambdas,
        which the distributed (condor) pickler cannot serialize when the output is
        shipped back from the workers. Any array that still carries that behavior
        (e.g. a saved `ljs`/object collection, or a weight that kept its behavior)
        breaks pickling. So we strip behavior and record/array names and store plain
        data. numpy arrays and Python lists are wrapped directly; a bare scalar
        becomes a length-1 array.
        """
        if isinstance(value, ak.Array):
            arr = value
        else:
            try:
                arr = ak.Array(value)
            except (ValueError, TypeError):
                arr = ak.Array([value])

        # drop behavior (the unpicklable lambdas live here) and metadata names
        return ak.Array(ak.without_parameters(arr).layout)

    def make_vector(self, objs, collection, fields, type_id=None, mass=None):
        shape = ak.ones_like(objs[collection].pt, dtype=np.dtype(int))
        # all objects must have the same fields to later concatenate and cluster them
        # set fields that aren't available for a given object to be -1
        # these additional fields will be removed after clustering anyway
        forms = {f: objs[collection][f] if f in objs[collection].fields else -1*shape for f in fields}
        forms["part_type"] = objs[collection]["type"] if type_id is None else type_id*shape
        forms["mass"] = objs[collection]["mass"] if mass is None else mass*shape
        if type_id == 8:
            forms["trkNumPixelHits"] = 0*shape
            forms["trkNumTrkLayers"] = 0*shape
        if type_id == 4:
            forms["lostHits"] = 999*shape
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

        unsafe_fields = ['muonIdxG','dsaIdxG','matched_muons','matched_dsa_muons','good_matched_muons','good_matched_dsa_muons']

        all_fields = list(set().union(*fields))
        for field in unsafe_fields:
            try:
                all_fields.remove(field)
            except ValueError:
                continue

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
        ljs["constituents"] = self.make_constituent(consts, [2, 3, 4, 8], "PtEtaPhiMCollection", common_fields)


    ######
        ## FIX ME! Won't be able to access the dsaMuon matches from the LJ constituent muon, and vice versa
        ## (can only access it from the original muon collection in objects)

        objs["dsaMuons"]["mass"] = ak.full_like(objs["dsaMuons"].pt, 0.105712890625)

        safe_pf_fields = list(objs["muons"].fields)
        safe_dsa_fields = list(objs["dsaMuons"].fields) +  ["trkNumPixelHits","trkNumTrkLayers" ]

        for field in unsafe_fields:
            if field in safe_pf_fields:
                safe_pf_fields.remove(field)
            if field in safe_dsa_fields:
                safe_dsa_fields.remove(field)

        extra_muon_fields =  ["trkNumPixelHits","trkNumTrkLayers" ]
        muon_fields = list(set(safe_pf_fields).intersection(safe_dsa_fields)) + extra_muon_fields
       

        ljs["muons"] = self.make_constituent(consts, [3, 8], "Muon", muon_fields)
        ljs["pfMuons"] = self.make_constituent(consts, [3], "Muon", safe_pf_fields)
        ljs["dsaMuons"] = self.make_constituent(consts, [8], "DSAMuon", safe_dsa_fields)
    ######
        extra_egamma_fields  = ["lostHits"]
        safe_electron_fields = list(objs["electrons"].fields)
        safe_photon_fields = list(objs["photons"].fields)
        egamma_fields  =  list(set(safe_electron_fields).intersection(safe_photon_fields)) + extra_egamma_fields
        ljs ["egamma"]  = self.make_constituent(consts, [2, 4], "Egamma", egamma_fields)
        ljs["electrons"] = self.make_constituent(consts, [2], "Electron",safe_electron_fields )
        ljs["photons"] = self.make_constituent(consts, [4], "Photon", safe_photon_fields)

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
            ch_cuts[channel]["preLj_obj"] = {}
            ch_cuts[channel]["lj"] = {}
            ch_cuts[channel]["postLj_obj"] = {}
            ch_cuts[channel]["evt"] = {}

            cuts = selection_menu[channel]
            for obj, obj_cuts in cuts["obj_cuts"].items():
                if obj not in ch_cuts[channel]["obj"]:
                    ch_cuts[channel]["obj"][obj] = []
                ch_cuts[channel]["obj"][obj] = utilities.flatten(obj_cuts)

            if "preLj_obj_cuts" in cuts:
                for obj, obj_cuts in cuts["preLj_obj_cuts"].items():
                    ch_cuts[channel]["preLj_obj"][obj] = utilities.flatten(obj_cuts)

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
            if len(output["metadata"]["is_data"]) != 1 or len(output["metadata"]["year"]) != 1:
                print(f"WARNING: {sample} has more than one value for is_data or year. Not scaling histograms or cutflows.")
                continue

            if output["metadata"]["is_data"].pop():
                print(f"{sample} is data. Not scaling histograms or cutflows.")
                continue

            print(f"{sample} is simulation. Scaling histograms or cutflows according to lumi*xs.")
            year = output["metadata"]["year"].pop()
            sum_weights = output["metadata"]["scaled_sum_weights"]
            lumixs_weight = utilities.get_lumixs_weight(sample, year, sum_weights)
            for name in output["cutflow"]:
                accumulator[sample]["cutflow"][name].scale(lumixs_weight)
            if not self.unweighted_hist:
                for name in output["hists"]:
                    accumulator[sample]["hists"][name] *= lumixs_weight