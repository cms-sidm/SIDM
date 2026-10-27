"""Synthetic LJ-selection regressions; no NanoEvents, files or dataset access."""

import contextlib
import copy
import io
import itertools
import sys
import unittest
from unittest.mock import patch

import awkward as ak
import numpy as np
from coffea.nanoevents.methods import vector

import sidm
from sidm.definitions.cuts import evt_cut_defs, obj_cut_defs
from sidm.definitions.hists import hist_defs
from sidm.definitions.objects import add_lj_candidates, lj_parents
from sidm.tools.selection import JaggedSelection, Selection
from sidm.tools.sidm_processor import SidmProcessor


def jet(idx, mu=2, pf=None, ele=0, photon=0, pt=None, phi=0.0,
        isolation=0.01, pixels=1, lost=1):
    """A reconstructed LJ with literal constituents and stable identity."""
    pf = mu if pf is None else pf
    pt = 100.0 - idx if pt is None else pt
    constituent = dict(pt=30.0, eta=0.1, phi=phi, mass=0.1,
                       trkNumPixelHits=pixels, trkNumTrkLayers=8, lostHits=lost)
    pf_muons = [dict(constituent) for _ in range(pf)]
    dsa_muons = [dict(constituent, trkNumPixelHits=0) for _ in range(mu - pf)]
    electrons = [dict(constituent) for _ in range(ele)]
    photons = [dict(constituent, lostHits=999) for _ in range(photon)]
    return dict(idx=idx, pt=pt, eta=0.1, phi=phi, mass=3.0,
                muon_n=mu, pfMu_n=pf, dsaMu_n=mu - pf,
                electron_n=ele, photon_n=photon, isolation=isolation,
                muons=pf_muons + dsa_muons, pfMuons=pf_muons,
                dsaMuons=dsa_muons, electrons=electrons, photons=photons,
                egamma=electrons + photons, dRSpread=0.05,
                dzSpread_mu=0.0, dxySpread_mu=0.0,
                dzSpread_pf=0.0, dxySpread_pf=0.0)


def objects(events):
    """Keep a concrete schema even for empty events and zero-event batches."""
    schema_event = [jet(0, mu=2, pf=1, ele=1, photon=1)]
    ljs = ak.with_name(ak.Array(list(events) + [schema_event])[:-1],
                       "PtEtaPhiMLorentzVector", behavior=vector.behavior)
    return dict(ljs=ljs, muons=ak.flatten(ljs.pfMuons, axis=2),
                dsaMuons=ak.flatten(ljs.dsaMuons, axis=2),
                evt_weights=ak.Array(np.ones(len(events))),
                event_id=ak.Array(np.arange(len(events))))


def selected(events, channel="4mu_sr"):
    """Exercise the real processor stage and expanded YAML subtype cuts."""
    processor = SidmProcessor([channel], [])
    cuts = processor.build_cuts()[channel]
    # is_data controls only whether signal-only post-LJ views are requested in
    # this stage. All records are synthetic, with no data/MC input opened.
    return processor.select_lepton_jets(objects(events), cuts, is_data=True)


def ids(objs, name):
    return ak.to_list(objs[name].idx)


PAIR_HISTS = tuple(name for name in hist_defs if name.startswith("lj_lj_"))
MIXED_HISTS = (
    "mulj_egmlj_invmass", "mulj_egmlj_invmass_pixelHits_lostHits",
    "mulj_egmlj_invmass_trkLayers_lostHits", "pf_mulj_egmlj_invmass",
    "1pf_mulj_egmlj_invmass", "pf_mulj_egmlj_invmass_pixelHits",
    "1pf_mulj_egmlj_invmass_pixelHits2", "1pf_mulj_egmlj_invmass_pixelHits1",
    "1pf_mulj_egmlj_invmass_pixelHits0", "1pf_mulj_egmlj_invmass_pixelHits3",
    "1pf_mulj_egmlj_invmass_pixelHits4", "pf_dsa_mulj_egmlj_invmass_pixelHits",
    "pf_mulj_egmlj_invmass_trkLayers", "1pf_mulj_egmlj_invmass_trkLayers10",
    "1pf_mulj_egmlj_invmass_trkLayers12", "1pf_mulj_egmlj_invmass_trkLayers11",
    "1pf_mulj_egmlj_invmass_trkLayers9", "1pf_mulj_egmlj_invmass_trkLayers8",
    "1pf_mulj_egmlj_invmass_trkLayers7", "pf_dsa_mulj_egmlj_invmass_trkLayers",
    "pf_dsa_mulj_egmlj_invmass", "dsa_mulj_egmlj_invmass",
    "mulj_e_lj_invmass", "mulj_1e_lj_invmass", "mulj_e_lj_invmass_lostHits",
    "mulj_1e_lj_invmass_lostHits1", "mulj_1e_lj_invmass_lostHits0",
    "mulj_1e_lj_invmass_lostHits2", "mulj_eg_lj_invmass",
    "mulj_eg_lj_invmass_lostHits", "mulj_g_lj_invmass",
)


class LJCollectionsTest(unittest.TestCase):
    """Selection membership, identity, channel and histogram invariants."""

    def setUp(self):
        self.output = io.StringIO()
        self.redirect = contextlib.redirect_stdout(self.output)
        self.redirect.__enter__()

    def tearDown(self):
        self.redirect.__exit__(None, None, None)
        self.assertNotIn("Warning:", self.output.getvalue())

    def test_missing_or_failed_original_pair_never_promotes(self):
        events = [[], [jet(0)],
                  [jet(0, isolation=0.2), jet(1), jet(2)],
                  [jet(0), jet(1, isolation=0.2), jet(2)],
                  [jet(0, pixels=3), jet(1), jet(2)],
                  [jet(0), jet(1, pixels=3), jet(2), jet(3)],
                  [jet(0, mu=1), jet(1), jet(2)],
                  [jet(0), jet(1), jet(2, isolation=0.2)],
                  [jet(0), jet(1), jet(2), jet(3)]]
        objs = selected(events)
        self.assertEqual(ids(objs, "candidate_ljs"), [[]] * 7 + [[0, 1], [0, 1]])
        self.assertEqual(ak.to_list(evt_cut_defs["4mu"](objs)), [False] * 7 + [True, True])
        self.assertEqual(ids(objs, "ljs"), [[j["idx"] for j in e] for e in events])
        # Diagnostic multiplicity remains distinct from candidate acceptance.
        self.assertTrue(ak.to_list(evt_cut_defs[">=2 LJs"](objs))[2])

    def test_egamma_failure_and_extra_jets(self):
        electron = jet(1, mu=0, ele=2, phi=2.5)
        events = [[jet(0), electron, jet(2, mu=0, ele=1)],
                  [jet(0), dict(electron, isolation=0.3), jet(2, mu=0, ele=1)],
                  [jet(0), jet(1, mu=0, ele=1, lost=0), jet(2, mu=0, ele=1)],
                  [jet(0), jet(1), jet(2, mu=0, ele=1)],
                  [jet(0), electron, jet(2), jet(3, mu=0, ele=1)]]
        objs = selected(events, "2mu2e_sr")
        self.assertEqual(ak.to_list(evt_cut_defs["2mu2e"](objs)), [True, False, False, False, True])
        self.assertEqual(ak.to_list(evt_cut_defs["4mu"](objs)), [False, False, False, True, False])
        self.assertEqual(ids(objs, "candidate_ljs"), [[0, 1], [], [], [0, 1], [0, 1]])

    def test_channel_exclusivity_for_three_and_four_jets(self):
        events = []
        expected_4mu, expected_mixed = [], []
        for count in (2, 3, 4):
            for kinds in itertools.product((0, 1, 2), repeat=count):
                events.append([jet(i, mu=kind, ele=int(kind == 0)) for i, kind in enumerate(kinds)])
                expected_4mu.append(kinds[:2] == (2, 2))
                expected_mixed.append(kinds[:2] in ((0, 2), (2, 0)))
        objs = selected(events)
        four, mixed = evt_cut_defs["4mu"](objs), evt_cut_defs["2mu2e"](objs)
        self.assertEqual(ak.to_list(four), expected_4mu)
        self.assertEqual(ak.to_list(mixed), expected_mixed)
        self.assertFalse(bool(ak.any(four & mixed)))

    def test_original_indices_after_generic_cuts_and_equal_pt(self):
        objs = selected([[jet(0, pt=20), jet(3, pt=80), jet(8, pt=80), jet(11, pt=70)]])
        self.assertEqual(ids(objs, "candidate_ljs"), [[3, 8]])
        self.assertEqual(ids(objs, "mu_ljs"), [[3, 8, 11]])
        failed = selected([[jet(0, pt=20), jet(3, pt=80, isolation=0.4),
                            jet(8, pt=80), jet(11, pt=70)]])
        self.assertEqual(ids(failed, "candidate_ljs"), [[]])

    def test_parent_cuts_refresh_pure_children_keep_mixed(self):
        objs = selected([[jet(0), jet(1, pf=0), jet(2, pf=1),
                          jet(3, isolation=0.4), jet(4, pf=0, isolation=0.4),
                          jet(5, mu=0, ele=1), jet(6, mu=0, photon=1),
                          jet(7, mu=0, ele=1, photon=1),
                          jet(8, mu=0, ele=1, isolation=0.4),
                          jet(9, mu=0, photon=1, isolation=0.4)]])
        self.assertEqual(ids(objs, "mu_ljs"), [[0, 1, 2]])
        self.assertEqual(ids(objs, "pfmu_ljs"), [[0]])
        self.assertEqual(ids(objs, "dsamu_ljs"), [[1]])
        self.assertEqual(ids(objs, "egm_ljs"), [[5, 6, 7]])
        self.assertEqual(ids(objs, "electron_ljs"), [[5]])
        self.assertEqual(ids(objs, "photon_ljs"), [[6]])
        for child, parent in lj_parents.items():
            self.assertLessEqual(set(ids(objs, child)[0]), set(ids(objs, parent)[0]))

    def test_child_cuts_are_preserved_in_either_order(self):
        processor = SidmProcessor([], [])
        raw = objects([[jet(0), jet(1), jet(2, isolation=0.4), jet(3, pf=1)]])
        custom = {"drop first": lambda objs: objs["pfmu_ljs"].idx != 0}
        for order in (0, 1):
            cuts = [("pfmu_ljs", ["drop first"]), ("mu_ljs", ["mu_lj_iso < 0.1"])]
            if order:
                cuts.reverse()
            with patch.dict(obj_cut_defs, {"pfmu_ljs": custom}):
                objs = processor.select_lepton_jets(raw, {"lj": {}, "postLj_obj": dict(cuts)}, True)
            self.assertEqual(ids(objs, "pfmu_ljs"), [[1]])
            self.assertEqual(ids(objs, "mu_ljs"), [[0, 1, 3]])
            self.assertEqual(ids(objs, "candidate_ljs"), [[0, 1]])
        self.assertEqual(ids(raw, "ljs"), [[0, 1, 2, 3]])

    def test_generic_parent_cut_prunes_existing_descendants(self):
        objs = selected([[jet(0, pt=45), jet(1, pt=80), jet(2, pf=0, pt=70)]])
        objs = JaggedSelection({"ljs": ["pT > 50 GeV"]}).apply_obj_cuts(objs)
        self.assertEqual(ids(objs, "mu_ljs"), [[1, 2]])
        self.assertEqual(ids(objs, "pfmu_ljs"), [[1]])
        self.assertEqual(ids(objs, "dsamu_ljs"), [[2]])
        self.assertEqual(ids(objs, "candidate_ljs"), [[1, 2]])

    def test_later_subtype_cut_refreshes_pair_without_promoting(self):
        objs = selected([[jet(0, isolation=0.15), jet(1), jet(2)]], "4mu")
        self.assertEqual(ids(objs, "candidate_ljs"), [[0, 1]])
        objs = JaggedSelection({"mu_ljs": ["mu_lj_iso < 0.1"]}).apply_obj_cuts(objs)
        self.assertEqual(ids(objs, "mu_ljs"), [[1, 2]])
        self.assertEqual(ids(objs, "candidate_ljs"), [[]])

    def test_event_cuts_require_candidate_setup(self):
        with self.assertRaises(KeyError):
            Selection(["4mu"]).apply_evt_cuts(objects([[jet(0), jet(1)]]))

    def test_invalid_identity_cannot_silently_skip_refresh(self):
        objs = selected([[jet(0), jet(1)]])
        objs["mu_ljs"] = ak.without_field(objs["mu_ljs"], "idx")
        with self.assertRaises(AttributeError):
            JaggedSelection({"mu_ljs": ["Mu >= 2"]}).apply_obj_cuts(objs)
        with self.assertRaises(AttributeError):
            add_lj_candidates(objs)

    def test_declared_custom_child_is_not_rebuilt(self):
        processor = SidmProcessor([], [])
        calls = []

        def custom_view(objs):
            calls.append(1)
            return objs["mu_ljs"]

        processor.postLj_objs = {**processor.postLj_objs, "custom_ljs": custom_view}
        cuts = {"lj": {}, "postLj_obj": {
            "custom_ljs": ["drop first"], "mu_ljs": ["mu_lj_iso < 0.1"]}}
        custom_cuts = {"drop first": lambda objs: objs["custom_ljs"].idx != 0}
        with patch.dict(lj_parents, {"custom_ljs": "mu_ljs"}), \
                patch.dict(obj_cut_defs, {"custom_ljs": custom_cuts}):
            objs = processor.select_lepton_jets(
                objects([[jet(0), jet(1), jet(2, isolation=0.4)]]), cuts, True)
        self.assertEqual(calls, [1])
        self.assertEqual(ids(objs, "custom_ljs"), [[1]])

    def test_mc_extra_views_do_not_leak_into_data_stage(self):
        processor = SidmProcessor([], [])
        raw = objects([[jet(0), jet(1)]])
        leptons = ak.Array([[dict(pt=25.0, parent=dict(pdgId=32)),
                             dict(pt=20.0, parent=dict(pdgId=23))]])
        raw.update(genMus=leptons, genEs=leptons)
        cuts = {"lj": {}, "postLj_obj": {}}
        mc = processor.select_lepton_jets(raw, cuts, False)
        self.assertEqual(ak.to_list(ak.num(mc["genMus_fromA"])), [1])
        self.assertEqual(ak.to_list(ak.num(mc["genEs_fromA"])), [1])
        data_stage = processor.select_lepton_jets(objects([[jet(0), jet(1)]]), cuts, True)
        self.assertNotIn("genMus_fromA", data_stage)

    def test_none_criterion_rejects_instead_of_placeholder(self):
        raw = objects([[jet(0), jet(1), jet(2)]])
        processor = SidmProcessor([], [])
        custom = {"defined": lambda objs: ak.Array([[True, None, True]])}
        with patch.dict(obj_cut_defs["mu_ljs"], custom):
            objs = processor.select_lepton_jets(raw, {"lj": {}, "postLj_obj": {"mu_ljs": ["defined"]}}, True)
        self.assertEqual(ids(objs, "mu_ljs"), [[0, 2]])
        self.assertEqual(ids(objs, "pfmu_ljs"), [[0, 2]])
        self.assertEqual(ids(objs, "candidate_ljs"), [[]])

    def test_event_cuts_mask_every_view_and_cutflow(self):
        objs = selected([[jet(0, phi=-1), jet(1, phi=1.5)],
                         [jet(0, isolation=0.5), jet(1, phi=2.5), jet(2)],
                         [jet(0), jet(1, phi=1)]])
        objs["evt_weights"] = ak.Array([2.0, 3.0, 5.0])
        selection = Selection([">=2 LJs", "LJ-LJ dPhi > 2", "4mu"])
        passing = selection.apply_evt_cuts(objs)
        self.assertEqual(ak.to_list(passing["event_id"]), [0])
        for name in ("ljs", "mu_ljs", "pfmu_ljs", "candidate_ljs", "candidate_mu_ljs"):
            self.assertEqual(ids(passing, name), [[0, 1]])
        self.assertEqual(selection.cutflow.rows["4mu"], {"raw": 1, "weighted": 2.0})

    def test_empty_batches_and_incomplete_pairs_are_false(self):
        for events in ([], [[]], [[jet(0)]]):
            objs = selected(events)
            for cut in ("4mu", "2mu2e", "LJ-LJ dPhi > 2", "inv(lj1, lj2) <= 100"):
                self.assertEqual(ak.to_list(evt_cut_defs[cut](objs)), [False] * len(events))
            for name in PAIR_HISTS + MIXED_HISTS + ("mu_lj_sum_pt_score",):
                histogram = copy.deepcopy(hist_defs[name])
                histogram.make_hist(name)
                histogram.fill(objs, objs["evt_weights"])
                self.assertEqual(histogram.hist.sum(flow=True).value, 0.0, name)

    def test_pair_histograms_and_cuts_use_identical_objects(self):
        events = [[jet(0, phi=-1), jet(1, phi=1.5), jet(2, mu=0, ele=1)],
                  [jet(0, phi=0.2), jet(1, mu=0, ele=1, phi=2.7), jet(2)],
                  [jet(0, isolation=0.4), jet(1), jet(2)],
                  [jet(0), jet(1, phi=0.1)]]
        objs = selected(events)
        objs["evt_weights"] = ak.Array([2.0, 3.0, 5.0, 7.0])
        expected_pair_mask = [True, True, False, True]
        pair = objs["candidate_ljs"]
        self.assertEqual(ak.to_list(evt_cut_defs["inv(lj1, lj2) <= 100"](objs)),
                         [False, False, False, True])
        dphi_hist = hist_defs["lj_lj_absdphi"]
        dphi = dphi_hist.axes[0].fill_func(objs, dphi_hist.evt_mask(objs))
        self.assertEqual(ak.to_list(evt_cut_defs["LJ-LJ dPhi > 2"](objs)), [True, True, False, False])
        self.assertTrue(np.allclose(ak.to_numpy(dphi), [2.5, 2.5, 0.1]))
        for name in PAIR_HISTS:
            with self.subTest(histogram=name):
                histogram = copy.deepcopy(hist_defs[name])
                self.assertEqual(ak.to_list(histogram.evt_mask(objs)), expected_pair_mask)
                histogram.make_hist(name)
                histogram.fill(objs, objs["evt_weights"])
                total = histogram.hist.sum(flow=True)
                self.assertEqual(total.value, 12.0)
                self.assertEqual(total.variance, 62.0)
                for axis in histogram.axes:
                    if axis.name == "ljlj_mass":
                        values = axis.fill_func(objs, histogram.evt_mask(objs))
                        self.assertTrue(np.allclose(ak.to_numpy(values), ak.to_numpy(pair[expected_pair_mask].sum().mass)))

    def test_mixed_mass_variants_do_not_pick_extra_jet(self):
        events = [[jet(0), jet(1), jet(2, mu=0, ele=1)],
                  [jet(0), jet(1, mu=0, ele=1, phi=2.5), jet(2)],
                  [jet(0, pf=1), jet(1, mu=0, ele=1, photon=1, phi=2.5)],
                  [jet(0, pf=0), jet(1, mu=0, photon=1, phi=2.5)],
                  [jet(0, isolation=0.5), jet(1, mu=0, ele=1), jet(2)]]
        objs = selected(events)
        for name in MIXED_HISTS:
            with self.subTest(histogram=name):
                histogram = copy.deepcopy(hist_defs[name])
                mask = histogram.evt_mask(objs)
                self.assertFalse(ak.to_list(mask)[0])
                self.assertFalse(ak.to_list(mask)[4])
                values = histogram.axes[0].fill_func(objs, mask)
                self.assertTrue(np.allclose(ak.to_numpy(values), ak.to_numpy(objs["candidate_ljs"][mask].sum().mass)))
                histogram.make_hist(name)
                histogram.fill(objs, objs["evt_weights"])
                self.assertEqual(histogram.hist.sum(flow=True).value, int(ak.sum(mask)))
        self.assertEqual(ak.to_list(hist_defs["mulj_egmlj_invmass"].evt_mask(objs)),
                         [False, True, True, True, False])
        self.assertEqual(ak.to_list(hist_defs["mu_lj_sum_pt_score"].evt_mask(objs)),
                         [True, False, False, False, False])

    def test_inclusive_diagnostics_retain_extras(self):
        objs = selected([[jet(0), jet(1), jet(2), jet(3, isolation=0.4)]])
        for name, expected in (("lj_pt", 4), ("mu_lj_pt", 3), ("pfmu_lj_pt", 3),
                               ("mu_lj_sum_pt_score", 1)):
            histogram = copy.deepcopy(hist_defs[name])
            histogram.make_hist(name)
            histogram.fill(objs, objs["evt_weights"])
            self.assertEqual(histogram.hist.sum(flow=True).value, expected)


if __name__ == "__main__":
    print(f"Python: {sys.executable}", flush=True)
    print(f"SIDM: {sidm.__file__}", flush=True)
    unittest.main(verbosity=2)
