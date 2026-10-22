"""Standard CMS corrections through correctionlib

configs/corrections.yaml says, per year, which JSON file holds each correction,
what selection it corrects, and whether it was measured for objects like ours.
This module turns that into event weights with up/down variations.

From the processor:

    SidmProcessor(channels, hists, corrections=True)             # the config's default list
    SidmProcessor(channels, hists, corrections=["pileup"])       # a subset

From a notebook:

    from sidm.tools import corrections
    corr = corrections.Corrections("2018")
    weights = corr.event_weights(events)             # genWeight, pileup, prefiring
    corr.add_object_weights(weights, ljs, obj_cuts)  # SFs for the objects in the lepton jets
    weights.weight(), weights.weight("pileupUp")

Only event-weight corrections are applied by the processor. The jet and
b-tagging helpers at the bottom of Corrections (jec_factor, jer_smearing_factor,
...) are for studies; the README in sidm/studies/corrections says why they are
not applied yet.
"""

import ctypes
import functools
import hashlib
import importlib
import importlib.util
import os
import sys

import awkward as ak
import numpy as np

from sidm import BASE_DIR
from sidm.tools.utilities import load_yaml

CONFIG = "configs/corrections.yaml"
REPO_ROOT = os.path.dirname(BASE_DIR)
SHIFTS = ("Up", "Down")


def import_correctionlib():
    """Import correctionlib so that its C++ core binds the library it ships with

    Under cmsenv, LD_LIBRARY_PATH lists CMSSW's external lib directory, which
    has its own libcorrectionlib.so (2.7.0 in CMSSW_17_0_0_pre4, older in the
    CMSSW_14 used by the Condor jobs). The PyPI wheel's _core extension finds
    its bundled copy through a RUNPATH, which the loader searches after
    LD_LIBRARY_PATH, so _core binds the CMSSW copy and the import fails with
    "undefined symbol: __cxa_call_terminate". Loading the bundled copy by path
    first makes the loader reuse it. This is harmless when there is no clash.
    """
    if "correctionlib" in sys.modules:
        return sys.modules["correctionlib"]
    spec = importlib.util.find_spec("correctionlib")
    if spec is None:
        raise ImportError("correctionlib is not installed")
    bundled = os.path.join(os.path.dirname(spec.origin), "lib", "libcorrectionlib.so")
    if os.path.exists(bundled):
        ctypes.CDLL(bundled, mode=ctypes.RTLD_GLOBAL)
    return importlib.import_module("correctionlib")


@functools.lru_cache(maxsize=None)
def load_correction_set(path):
    """Load a correctionlib JSON (plain or gzipped) once per process"""
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"correction file {path} not found. The CMS corrections live on CVMFS; "
            "set SIDM_CORRECTIONS_ROOT to use a local copy of the same tree.")
    return import_correctionlib().CorrectionSet.from_file(path)


@functools.lru_cache(maxsize=None)
def md5sum(path):
    """Checksum of a correction file, recorded as provenance"""
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def load_config(cfg=CONFIG):
    """The corrections config; a relative path is relative to the sidm package"""
    return load_yaml(cfg if os.path.isabs(cfg) else f"{BASE_DIR}/{cfg}")


def resolve_names(names, cfg=CONFIG):
    """Turn the processor's `corrections` argument into a list of names

    None, False or [] mean no corrections; True means the config's default
    list; a string is a single correction.
    """
    if not names:
        return []
    if names is True:
        return list(load_config(cfg)["default"])
    if isinstance(names, str):
        return [names]
    return list(names)


def variation_names(names):
    """Names of the up/down variations the given corrections produce"""
    return [f"{name}{shift}" for name in names for shift in SHIFTS]


def weight_categories(names, monitored=()):
    """Every weight category the processor fills besides "nominal"

    "uncorrected" is the plain genWeight and "only_<name>" the genWeight times
    one correction, as references for what each correction does; then the
    up/down variations of the applied corrections. A monitored correction is
    evaluated but not applied: "with_<name>" is the nominal weight times it,
    and "with_<name>Up"/"Down" the same with it shifted.
    """
    out = ["uncorrected"] + [f"only_{name}" for name in names] + variation_names(names)
    for name in monitored:
        out += [f"with_{name}"] + [f"with_{name}{shift}" for shift in SHIFTS]
    return out


def evaluate(correction, *args):
    """correction.evaluate for scalars, numpy arrays, or awkward arrays

    Awkward arguments are flattened, evaluated in one vectorised call and given
    back the structure of the first awkward argument, which may be flat (one
    value per event) or jagged (one value per object). Other arguments, such as
    the string that picks a variation, are broadcast as they are. Missing values
    must be filled first (ak.fill_none): flattening would drop them.
    """
    template = next((a for a in args if isinstance(a, ak.Array)), None)
    if template is None:
        return correction.evaluate(*args)
    if template.ndim > 2:
        raise ValueError("evaluate() handles one level of nesting; flatten the array first")
    flat = [np.asarray(ak.flatten(a, axis=None)) if isinstance(a, ak.Array) else a
            for a in args]
    out = correction.evaluate(*flat)
    if template.ndim == 1:
        return ak.Array(out)
    return ak.unflatten(out, ak.num(template, axis=1))


def _evaluate_numpy(correction, args):
    """correction.evaluate on numpy arguments, some of which may be string arrays

    correctionlib takes arrays only for int and real inputs, so string inputs
    that vary per object (e.g. the photon barrel/endcap bin) are split into one
    call per distinct value.
    """
    categorical = [i for i, a in enumerate(args)
                   if isinstance(a, np.ndarray) and a.dtype.kind in "US"]
    if not categorical:
        return correction.evaluate(*args)
    size = len(args[categorical[0]])
    out = np.empty(size)
    for combo in set(zip(*(args[i] for i in categorical))):
        sel = np.ones(size, dtype=bool)
        for i, value in zip(categorical, combo):
            sel &= args[i] == value
        values = dict(zip(categorical, combo))
        out[sel] = correction.evaluate(*[
            values[i] if i in values else (a[sel] if isinstance(a, np.ndarray) else a)
            for i, a in enumerate(args)])
    return out


def _branch(events, name):
    """events.Pileup.nTrueInt from "Pileup_nTrueInt", the NanoAOD branch name"""
    if name in events.fields:
        return events[name]
    prefix, _, field = name.partition("_")
    return events[prefix][field]


class EventWeights:
    """The event weight as a product of named factors with optional up/down shifts

    Follows the interface of coffea.analysis_tools.Weights (add, weight,
    variations), but builds each variation as a product of factors instead of
    a ratio to the nominal weight, so a factor that is zero for some events
    cannot turn a variation into NaN.
    """

    def __init__(self, size, factors=None):
        self.size = size
        self._factors = dict(factors or {})

    def add(self, name, nominal, up=None, down=None):
        """Add a factor; up and down are its shifted values, not ratios"""
        if name in self._factors:
            raise ValueError(f"weight {name!r} already added")
        if name.endswith(SHIFTS):
            raise ValueError(f"weight names must not end in Up or Down: {name!r}")
        if (up is None) != (down is None):
            raise ValueError(f"give both up and down shifts for {name!r}, or neither")
        values = []
        for v in (nominal, up, down):
            if v is not None:
                v = np.asarray(ak.to_numpy(v) if isinstance(v, ak.Array) else v, dtype=np.float64)
                if v.shape != (self.size,):
                    raise ValueError(f"{name!r} has shape {v.shape}, expected ({self.size},)")
            values.append(v)
        self._factors[name] = tuple(values)

    def copy(self):
        """A copy that shares the factor arrays, for adding channel-specific factors"""
        return EventWeights(self.size, self._factors)

    @property
    def names(self):
        """Names of the factors, in the order they were added"""
        return list(self._factors)

    @property
    def variations(self):
        """Names of the available variations, e.g. pileupUp and pileupDown"""
        return [f"{name}{shift}" for name, (_, up, _) in self._factors.items()
                if up is not None for shift in SHIFTS]

    def factor(self, name, shift=None):
        """One factor: nominal, or shifted with shift="Up" or "Down" """
        nominal, up, down = self._factors[name]
        return {None: nominal, "Up": up, "Down": down}[shift]

    def weight(self, variation=None, include=None, exclude=()):
        """Product of the factors, with at most one of them shifted

        include/exclude restrict the product to a subset of factors, e.g.
        weight(include=["genWeight"]) is the uncorrected weight.
        """
        shifted, shift = None, None
        if variation is not None:
            if variation not in self.variations:
                raise KeyError(f"unknown variation {variation!r}; have {self.variations}")
            shift = "Up" if variation.endswith("Up") else "Down"
            shifted = variation[:-len(shift)]
        out = np.ones(self.size)
        for name in (self.names if include is None else include):
            if name not in exclude:
                out = out * self.factor(name, shift if name == shifted else None)
        return out

    def category(self, category, monitored=()):
        """The weight of "nominal" or of one category of weight_categories()

        The factors named in monitored enter only their own with_<name>
        categories. A category whose correction was not added here (an object
        SF the channel does not select) equals the weight without that factor,
        e.g. muon_idUp is the nominal weight.
        """
        applied = [name for name in self.names if name not in monitored]
        if category == "uncorrected":
            return self.weight(include=["genWeight"])
        if category.startswith("only_"):
            name = category[len("only_"):]
            return self.weight(include=["genWeight"] + [n for n in self.names if n == name])
        if category.startswith("with_"):
            name, variation = category[len("with_"):], None
            for shift in SHIFTS:
                if name.endswith(shift) and name[:-len(shift)] in monitored:
                    variation, name = name, name[:-len(shift)]
            if name not in self.names:
                return self.weight(include=applied)
            return self.weight(variation, include=applied + [name])
        if category in self.variations:
            return self.weight(category, include=applied)
        return self.weight(include=applied)

    def factor_columns(self, names):
        """{name: nominal, nameUp: up, nameDown: down} for each of names

        A factor that was not added is one, so that every channel gets the same
        columns; genWeight has no shifts.
        """
        ones = np.ones(self.size)
        out = {}
        for name in names:
            nominal, up, down = self._factors.get(name, (ones, None, None))
            out[name] = nominal
            if name != "genWeight":
                out[f"{name}Up"] = nominal if up is None else up
                out[f"{name}Down"] = nominal if down is None else down
        return out


class Corrections:
    """The corrections for one year, as configured in configs/corrections.yaml

    Holds only configuration, so it pickles cheaply to dask and Condor workers;
    the JSON files are loaded on first use in each process.
    """

    def __init__(self, year, names=True, cfg=CONFIG):
        menu = load_config(cfg)
        self.year = str(year)
        if self.year not in menu:
            raise KeyError(f"no corrections configured for {self.year}; add them to {cfg}")
        self.config = menu[self.year]
        self.json_root = os.environ.get("SIDM_CORRECTIONS_ROOT", menu["json_root"])
        self.names = resolve_names(names, cfg)
        unknown = [n for n in self.names if n not in self.config]
        if unknown:
            raise KeyError(f"corrections {unknown} are not configured for {self.year} in {cfg}")
        helpers = [n for n in self.names if self.config[n]["kind"] not in ("event", "object")]
        if helpers:
            raise ValueError(f"{helpers} are helpers, not event weights, and cannot be applied "
                             "by the processor")

    # --- files -------------------------------------------------------------

    def path(self, file):
        """Absolute path of a correction file named in the config"""
        if os.path.isabs(file):
            return file
        if file.startswith("sidm/"):
            return os.path.join(REPO_ROOT, file)
        return os.path.join(self.json_root, file)

    def correction(self, file, name):
        """One correction (or compound correction) from a file named in the config"""
        cset = load_correction_set(self.path(file))
        # `name in cset` raises IndexError for a missing key, so test the list
        return cset[name] if name in list(cset) else cset.compound[name]

    def entry(self, name):
        """The config entry of one correction"""
        return self.config[name]

    def files(self, name):
        """Every file a correction reads"""
        entry = self.config[name]
        paths = [entry[k] for k in ("file", "smear_file") if k in entry]
        for pieces in entry.get("working_points", {}).values():
            paths += [p["file"] for p in pieces]
        return sorted(set(paths))

    def provenance(self, names=None):
        """What was applied, from which file, with checksums, for metadata sidecars"""
        out = {"year": self.year, "json_root": self.json_root, "corrections": {}}
        for name in (self.names if names is None else names):
            entry = self.config[name]
            out["corrections"][name] = {
                "kind": entry["kind"],
                "validity": entry.get("validity"),
                "files": {f: md5sum(self.path(f)) for f in self.files(name)},
            }
        return out

    # --- event-level corrections -------------------------------------------

    def event_weights(self, events):
        """EventWeights with genWeight and every configured event-level correction"""
        weights = EventWeights(len(events))
        weights.add("genWeight", events.genWeight)
        for name in self.names:
            entry = self.config[name]
            if entry["kind"] == "event":
                weights.add(name, *self.event_factor(name, events))
        return weights

    def event_factor(self, name, events):
        """(nominal, up, down) of one event-level correction"""
        entry = self.config[name]
        if "branches" in entry:
            return tuple(_branch(events, entry["branches"][k]) for k in ("nominal", "up", "down"))
        corr = self.correction(entry["file"], entry["correction"])
        inputs = [_branch(events, b) for b in entry["inputs"]]
        return tuple(evaluate(corr, *inputs, entry["variations"][k])
                     for k in ("nominal", "up", "down"))

    # --- object-level corrections ------------------------------------------

    def working_point(self, name, obj_cuts):
        """The working point of an object correction that a channel selects, or None

        obj_cuts is the channel's {object: [cut names]}, as built by
        SidmProcessor.build_cuts.
        """
        entry = self.config[name]
        cuts = obj_cuts.get(entry["cut_object"], [])
        for wp in entry["working_points"]:
            if wp == "always" or wp in cuts:
                return wp
        return None

    def coverage(self, obj_cuts):
        """{object correction: working point or None} for one channel"""
        return {name: self.working_point(name, obj_cuts)
                for name in self.names if self.config[name]["kind"] == "object"}

    def add_object_weights(self, weights, ljs, obj_cuts):
        """Multiply in the object SFs that apply to one channel

        The SFs are evaluated for the constituents of the channel's lepton jets
        after the lepton-jet cuts, i.e. the objects the event selection is built
        from, and multiplied over each event. Up and down move every object's SF
        together (fully correlated). Returns {correction: working point} for the
        corrections that were applied.
        """
        applied = {}
        for name, wp in self.coverage(obj_cuts).items():
            if wp is None:
                continue
            objs = ak.flatten(ljs[self.config[name]["collection"]], axis=2)
            per_event = [ak.prod(self.object_sf(name, objs, wp, which), axis=1)
                         for which in ("nominal", "up", "down")]
            weights.add(name, *per_event)
            applied[name] = wp
        return applied

    def object_sf(self, name, objs, working_point, which="nominal"):
        """Per-object SFs, jagged like objs; objects outside the measured range get 1"""
        entry = self.config[name]
        counts = ak.num(objs, axis=1)
        values = self._object_inputs(entry, ak.flatten(objs, axis=1))
        values["scale_factors"] = values["ValType"] = entry["variations"][which]
        sf = np.ones(len(values["pt"]))
        for piece in entry["working_points"][working_point]:
            lo, hi = piece.get("pt", (-np.inf, np.inf))
            in_range = (values["pt"] >= lo) & (values["pt"] < hi)
            if not np.any(in_range):
                continue
            corr = self.correction(piece["file"], piece["correction"])
            # anything else in the piece is a fixed input, e.g. WorkingPoint
            fixed = {k: v for k, v in piece.items() if k not in ("file", "correction", "pt")}
            args = []
            for var in corr.inputs:
                value = fixed.get(var.name, values.get(var.name))
                if value is None:
                    raise KeyError(f"{name}: no value for input {var.name!r} of {piece['correction']}")
                args.append(value[in_range] if isinstance(value, np.ndarray) else value)
            sf[in_range] = _evaluate_numpy(corr, args)
        return ak.unflatten(sf, counts)

    def _object_inputs(self, entry, objs):
        """numpy inputs for object SFs, keyed by the input names the POG files use"""
        pt = np.asarray(objs.pt, dtype=np.float64)
        eta = np.asarray(objs.eta, dtype=np.float64)
        if entry["collection"] == "electrons":
            # EGM SFs are binned in supercluster eta
            eta = eta + np.asarray(objs.deltaEtaSC, dtype=np.float64)
        if "max_abs_eta" in entry:
            # the last eta bin is closed on the right with flow=error
            edge = np.nextafter(entry["max_abs_eta"], 0)
            eta = np.clip(eta, -edge, edge)
        values = {"pt": pt, "eta": eta, "abseta": np.abs(eta), "year": self.year}
        if entry["collection"] == "photons":
            barrel = np.asarray(ak.values_astype(objs.isScEtaEB, bool))
            values["HasPixBin"] = values["CSEVBin"] = np.where(barrel, "EBInc", "EEInc")
        return values

    # --- jet and b-tagging helpers (not applied by the processor) -----------

    def _jet(self, key):
        entry = self.config["jet_energy"]
        return self.correction(entry["file"], entry[key])

    def jec_factor(self, jets, rho):
        """L1L2L3 JEC recomputed from the raw jet pt (corrected over raw)

        NanoAODv9 jets are already corrected; comparing this with
        1/(1 - rawFactor) checks that they carry the current JEC tag.
        """
        raw_pt = jets.pt * (1 - jets.rawFactor)
        corr = self._jet("jec")
        inputs = {"JetA": jets.area, "JetEta": jets.eta, "JetPt": raw_pt,
                  "Rho": ak.broadcast_arrays(rho, jets.pt)[0]}
        return evaluate(corr, *[inputs[v.name] for v in corr.inputs])

    def jes_uncertainty(self, jets):
        """Relative total JES uncertainty per jet; shifted pt is pt * (1 +- unc)"""
        return evaluate(self._jet("jes_uncertainty"), jets.eta, jets.pt)

    def jer_smearing_factor(self, jets, rho, event_id, variation="nominal"):
        """Hybrid JER smearing factor per jet for MC (multiply pt and mass by it)

        Jets matched to a generator jet (dR < 0.2, |pt - pt_gen| < 3 sigma) are
        scaled towards it, others are smeared stochastically with a seed built
        from the jet and the event, so the result is reproducible. variation is
        "nominal", "up" or "down", which shift the JER SF by its uncertainty.
        """
        gen = jets.matched_gen
        columns = [jets.pt, jets.eta, ak.fill_none(gen.pt, -1.0),
                   ak.fill_none(jets.delta_r(gen), np.inf),
                   ak.broadcast_arrays(rho, jets.pt)[0],
                   ak.broadcast_arrays(event_id, jets.pt)[0]]
        flat = [np.asarray(ak.flatten(c, axis=None)) for c in columns]
        return ak.unflatten(self.jer_smearing_values(*flat, variation=variation),
                            ak.num(jets.pt, axis=1))

    def jer_smearing_values(self, pt, eta, gen_pt, gen_dr, rho, event_id, variation="nominal"):
        """jer_smearing_factor for flat numpy arrays with one entry per jet

        gen_pt and gen_dr are the pt of the generator jet NanoAOD matched to the
        jet (Jet_genJetIdx) and their distance; give gen_pt < 0 for a jet
        without one. Useful when the jets are optional, e.g. the jet matched to
        each lepton jet.
        """
        pt, eta, rho = (np.asarray(a, dtype=np.float64) for a in (pt, eta, rho))
        resolution = self._jet("jer_resolution").evaluate(eta, pt, rho)
        sf = self._jet("jer_sf").evaluate(eta, pt)
        if variation != "nominal":
            unc = self._jet("jer_sf_uncertainty").evaluate(eta, pt)
            sf = sf + unc if variation == "up" else sf - unc
        matched = (gen_pt >= 0) & (gen_dr < 0.2) & (np.abs(pt - gen_pt) < 3 * resolution * pt)
        entry = self.config["jet_energy"]
        smear = self.correction(entry["smear_file"], entry["smear"])
        inputs = {"JetPt": pt, "JetEta": eta, "GenPt": np.where(matched, gen_pt, -1.0),
                  "Rho": rho, "EventID": np.asarray(event_id, dtype=np.int64),
                  "JER": resolution, "JERSF": sf}
        return smear.evaluate(*[inputs[v.name] for v in smear.inputs])

    def jet_veto_mask(self, jets):
        """True for jets inside the jet veto map (the HEM 15/16 and hot regions)"""
        entry = self.config["jet_veto"]
        corr = self.correction(entry["file"], entry["correction"])
        # outside the map the flow value is 0, i.e. not vetoed
        return evaluate(corr, entry["type"], jets.eta, jets.phi) > 0

    def btag_working_point(self, working_point="T"):
        """DeepJet discriminant threshold of a fixed working point (L, M, T)"""
        entry = self.config["btag"]
        return self.correction(entry["file"], entry["wp_values"]).evaluate(working_point)

    def btag_sf(self, jets, working_point="T", systematic="central"):
        """Fixed-working-point DeepJet SF per jet (b and c from comb, light from incl)

        Jets outside the measured range (pt 20-1000 GeV, |eta| < 2.5) get 1.
        """
        entry = self.config["btag"]
        heavy = self.correction(entry["file"], entry["sf_heavy"])
        light = self.correction(entry["file"], entry["sf_light"])
        counts = ak.num(jets, axis=1)
        flat = ak.flatten(jets, axis=1)
        flavor = np.asarray(flat.hadronFlavour, dtype=np.int64)
        pt = np.asarray(flat.pt, dtype=np.float64)
        abseta = np.abs(np.asarray(flat.eta, dtype=np.float64))
        sf = np.ones(len(pt))
        in_range = (pt >= 20) & (pt < 1000) & (abseta < 2.5)
        for corr, is_flavor in ((heavy, flavor > 0), (light, flavor == 0)):
            sel = in_range & is_flavor
            if np.any(sel):
                sf[sel] = corr.evaluate(systematic, working_point, flavor[sel], abseta[sel], pt[sel])
        return ak.unflatten(sf, counts)

    # --- summaries ----------------------------------------------------------

    def summary(self, names=None):
        """[(name, kind, validity, note)] for the configured corrections"""
        names = [n for n, e in self.config.items() if isinstance(e, dict)] if names is None else names
        return [(n, self.config[n]["kind"], self.config[n].get("validity"),
                 " ".join(self.config[n].get("note", "").split())) for n in names]
