"""Synthetic angular regressions; no event files or processor execution.

Run from the repository root with its analysis dependencies:
    python -m unittest discover -s tests -p 'test_delta_phi.py' -v
"""

import copy
import io
import math
import unittest
from contextlib import redirect_stdout

import awkward as ak
import numpy as np
from coffea.nanoevents.methods import vector

from sidm.definitions.cuts import evt_cut_defs, obj_cut_defs
from sidm.definitions.hists import hist_defs
from sidm.tools.utilities import dR_outer, delta_phi, get_pairs


LJ_HISTS = (
    "lj_lj_absdphi", "lj_lj_absdphi_invmass", "lj_lj_absdphi_absdR",
    "lj_lj_absdphi_absdeta", "lj_lj_absdphi_ptRatio",
)
GEN_HISTS = {"genE_genE_dPhi": "genEs", "genMu_genMu_dPhi": "genMus"}
DSA_HISTS = (
    "dsa_pair_phi_diff", "back_to_back_dsa_pair_phi_diff",
    "parallel_dsa_pair_phi_diff",
)


def vectors(phis, eta_step=0.0):
    """Small in-memory vectors, with decreasing pt to preserve pair ordering."""
    phi = ak.Array(phis)
    index = ak.local_index(phi, axis=1)
    return ak.zip(
        {"pt": 60.0 - 10.0 * index, "eta": eta_step * index,
         "phi": phi, "mass": ak.ones_like(phi)},
        with_name="PtEtaPhiMLorentzVector", behavior=vector.behavior,
    )


def angular_distance(a, b):
    """Independent trigonometric reference, including antipodal directions."""
    return np.abs(np.arctan2(np.sin(a - b), np.cos(a - b)))


class DeltaPhiTests(unittest.TestCase):
    def assert_array_close(self, actual, expected):
        np.testing.assert_allclose(ak.to_numpy(actual), expected, atol=2e-14, rtol=1e-13)

    def fill_hist(self, name, objs, weights):
        definition = copy.deepcopy(hist_defs[name])
        definition.make_hist(name)
        output = io.StringIO()
        with redirect_stdout(output):
            definition.fill(objs, ak.Array(weights), verbose=True)
        self.assertEqual(output.getvalue(), "", "Histogram filling must not silently skip")
        return definition.hist

    def test_scalar_boundaries(self):
        self.assertEqual(delta_phi(0.0, 0.0), 0.0)
        self.assertEqual(delta_phi(math.pi, 0.0), -math.pi)
        self.assertEqual(delta_phi(-math.pi, 0.0), -math.pi)
        self.assertEqual(delta_phi(math.pi, -math.pi), 0.0)
        self.assertAlmostEqual(delta_phi(math.pi - .02, -math.pi + .02), -.04)
        self.assertAlmostEqual(delta_phi(-math.pi + .02, math.pi - .02), .04)

    def test_numpy_reference_swap_and_periodicity(self):
        first = np.linspace(-math.pi, math.pi, 73)[:, None]
        second = np.linspace(-math.pi, math.pi, 41)[None, :]
        result = delta_phi(first, second)
        self.assertTrue(np.all((-math.pi <= result) & (result < math.pi)))
        np.testing.assert_allclose(abs(result), angular_distance(first, second), atol=2e-14)
        np.testing.assert_allclose(abs(result), abs(delta_phi(second, first)), atol=2e-14)
        # At the antipodal tie only the absolute value is invariant to floating-point shifts.
        for turns in (-3, -1, 1, 4):
            np.testing.assert_allclose(abs(delta_phi(first + turns * 2 * math.pi, second)),
                                       abs(result), atol=2e-14)
            np.testing.assert_allclose(abs(delta_phi(first, second + turns * 2 * math.pi)),
                                       abs(result), atol=2e-14)
        away_from_tie = abs(result) < math.pi - 1e-12
        np.testing.assert_allclose(result[away_from_tie],
                                   -delta_phi(second, first)[away_from_tie], atol=2e-14)

    def test_jagged_missing_and_empty(self):
        first = ak.Array([[math.pi - .02, None], [], [-math.pi + .02, 0.0]])
        second = ak.Array([[-math.pi + .02, 1.0], [], [math.pi - .02, 0.0]])
        result = delta_phi(first, second)
        self.assertEqual(ak.to_list(ak.num(result)), [2, 0, 2])
        self.assertEqual(ak.to_list(ak.is_none(result, axis=-1)), [[False, True], [], [False, False]])
        self.assert_array_close(ak.flatten(result, axis=None), [-.04, .04, 0.0])
        self.assertEqual(ak.to_list(delta_phi(first[:0], second[:0])), [])
        self.assertEqual(ak.to_list(delta_phi(first[:, :0], 0.0)), [[], [], []])

    def test_float32_and_float64_histogram_endpoints(self):
        for dtype in (np.float32, np.float64):
            with self.subTest(dtype=dtype):
                first = np.array([math.pi, -math.pi, 0.0, math.pi - .02], dtype=dtype)
                second = np.array([0.0, 0.0, 0.0, -math.pi + .02], dtype=dtype)
                result = delta_phi(first, second)
                self.assertEqual(result.dtype, np.dtype(np.float64))
                self.assertTrue(np.all(abs(result) <= math.pi))
                np.testing.assert_allclose(abs(result), angular_distance(first.astype(np.float64),
                                                                        second.astype(np.float64)), atol=2e-14)
                objects = vectors([[0.0, math.pi], [0.0, -math.pi], [0.0, 0.0],
                                   [math.pi - .02, -math.pi + .02]])
                objects = ak.with_field(objects, ak.values_astype(objects.phi, dtype), "phi")
                for name in LJ_HISTS + DSA_HISTS:
                    key = "ljs" if name in LJ_HISTS else "dsaMuons"
                    definition = hist_defs[name]
                    objs = {key: objects}
                    actual = definition.axes[0].fill_func(objs, definition.evt_mask(objs))
                    flat = ak.to_numpy(ak.flatten(actual, axis=None))
                    self.assertEqual(flat.dtype, np.dtype(np.float64))
                    self.assertTrue(np.all(abs(flat) <= math.pi))
                    filled = self.fill_hist(name, objs, [1, 1, 1, 1])
                    self.assertEqual(filled.sum().value, len(flat))
                    self.assertEqual(filled.sum(flow=True).value, len(flat))
                    if name != "parallel_dsa_pair_phi_diff":
                        self.assertEqual(filled.project(0).values()[-1], 2)

    def lj_fixture(self):
        return vectors([[math.pi - .02, -math.pi + .02], [0.0, 2.5], [], [.2],
                        [0.0, math.pi], [1.0, 1.0],
                        [math.pi - .03, -math.pi + .04, 0.0]], eta_step=.15)

    def test_lj_histograms_fill_and_agree_with_cut(self):
        ljs = self.lj_fixture()
        objs = {"ljs": ljs}
        expected = [.04, 2.5, math.pi, 0.0, .07]
        weights = [2.0, 3.0, 101.0, 103.0, 5.0, 7.0, 11.0]
        selected = ak.num(ljs) > 1
        self.assert_array_close(abs(ljs[selected, 1].delta_phi(ljs[selected, 0])), expected)
        self.assertEqual(ak.to_list(evt_cut_defs["LJ-LJ dPhi > 2"]({"ljs": ljs[selected]})),
                         [False, True, True, False, False])
        for name in LJ_HISTS:
            with self.subTest(histogram=name):
                definition = hist_defs[name]
                self.assert_array_close(definition.axes[0].fill_func(objs, definition.evt_mask(objs)), expected)
                filled = self.fill_hist(name, objs, weights)
                axis = filled.axes[0]
                self.assertEqual(axis.edges[0], 0.0)
                self.assertAlmostEqual(axis.edges[-1], math.pi)
                self.assertEqual(axis.index(math.pi), axis.size - 1)
                self.assertAlmostEqual(filled.sum().value, 28.0)
                self.assertAlmostEqual(filled.sum().variance, 208.0)
                self.assertAlmostEqual(filled.sum(flow=True).value, 28.0)
                projected = filled.project(0)
                self.assertAlmostEqual(projected.values()[0], 7.0)
                self.assertAlmostEqual(projected.values()[-1], 5.0)

    def test_lj_histogram_swap_periodicity_and_empty_masks(self):
        ljs = self.lj_fixture()
        for name in LJ_HISTS:
            definition = hist_defs[name]
            with self.subTest(histogram=name):
                for changed in (
                    ak.concatenate([ljs[:, 1:2], ljs[:, :1], ljs[:, 2:]], axis=1),
                    ak.with_field(ljs, ljs.phi + 4 * math.pi * ak.local_index(ljs, axis=1), "phi"),
                ):
                    objs = {"ljs": changed}
                    self.assert_array_close(definition.axes[0].fill_func(objs, definition.evt_mask(objs)),
                                            [.04, 2.5, math.pi, 0.0, .07])
                for empty in (ljs[:0], ljs[:, :1], ljs[:, :0]):
                    objs = {"ljs": empty}
                    self.assertEqual(ak.to_list(definition.axes[0].fill_func(objs, definition.evt_mask(objs))), [])
                    self.assertEqual(self.fill_hist(name, objs, np.ones(len(empty))).sum(flow=True).value, 0)

    def test_generator_zoom_histograms(self):
        for name, obj in GEN_HISTS.items():
            with self.subTest(histogram=name):
                collection = self.lj_fixture()
                objs = {obj: collection}
                definition = hist_defs[name]
                self.assert_array_close(definition.axes[0].fill_func(objs, definition.evt_mask(objs)),
                                        [.04, 2.5, math.pi, 0.0, .07])
                filled = self.fill_hist(name, objs, [2, 3, 101, 103, 5, 7, 11])
                self.assertEqual(filled.axes[0].edges[-1], 1.0)
                self.assertEqual(filled.sum().value, 20)
                self.assertEqual(filled.sum(flow=True).value, 28)
                # Boundary-crossing close pairs belong inside the zoom, not overflow.
                self.assertEqual(filled.values(flow=True)[-1], 8)
                self.assertEqual(self.fill_hist(name, {obj: collection[:, :1]}, np.ones(len(collection))).sum(flow=True).value, 0)

    def test_dsa_pair_histograms_and_existing_pair_filters(self):
        dsa = vectors([[math.pi - .02, -math.pi + .02], [.02, -math.pi + .01],
                       [], [.3], [0.0, math.pi], [0.0, .1, math.pi]])
        expected = {
            DSA_HISTS[0]: ([[.04], [math.pi - .01], [], [], [math.pi], [.1, math.pi, math.pi - .1]], 31),
            DSA_HISTS[1]: ([[], [math.pi - .01], [], [], [math.pi], [math.pi, math.pi - .1]], 22),
            DSA_HISTS[2]: ([[.04], [], [], [], [], [.1]], 9),
        }
        for name, (values, total) in expected.items():
            with self.subTest(histogram=name):
                actual = hist_defs[name].axes[0].fill_func({"dsaMuons": dsa}, slice(None))
                self.assertEqual(ak.to_list(ak.num(actual)), [len(row) for row in values])
                self.assert_array_close(ak.flatten(actual), ak.flatten(ak.Array(values)))
                filled = self.fill_hist(name, {"dsaMuons": dsa}, [2, 3, 101, 103, 5, 7])
                self.assertEqual(filled.sum().value, total)
                self.assertEqual(filled.sum(flow=True).value, total)
                self.assertEqual(filled.axes[0].index(math.pi), filled.axes[0].size - 1)
                for empty in (dsa[:0], dsa[:, :1]):
                    self.assertEqual(self.fill_hist(name, {"dsaMuons": empty}, np.ones(len(empty))).sum(flow=True).value, 0)

    def test_outer_track_distance_and_real_cross_cleaning_cuts(self):
        def dsa(match_phi, segments=5, charge=1):
            return {"outerEta": 0.0, "outerPhi": math.pi - .02,
                    "nSegments": 10, "charge": 1,
                    "good_matched_muons": [{"outerEta": 0.0, "outerPhi": match_phi,
                                            "numMatch": segments, "charge": charge}]}

        close = -math.pi + .02
        far = math.pi - .22
        unmatched = dsa(close)
        unmatched["good_matched_muons"] = []
        collection = ak.Array([[dsa(close), dsa(far), dsa(close, segments=0),
                                dsa(close, segments=2), dsa(close, charge=-1),
                                unmatched, dsa(None), dsa(far, segments=10)], []])
        distances = dR_outer(collection[:, :, None], collection.good_matched_muons)
        self.assert_array_close(ak.flatten(distances, axis=None), [.04, .2, .04, .04, .04, np.inf, .2])
        self.assertEqual(ak.to_list(ak.num(distances, axis=-1)), [[1, 1, 1, 1, 1, 0, 1, 1], []])
        self.assertEqual(ak.to_list(dR_outer(collection[:0, :, None], collection[:0].good_matched_muons)), [])
        expected = {
            "segmatch veto + dR outer": [False, True, True, False, False, True, True, True],
            "old_all": [False, True, True, True, False, True, True, True],
            "all": [False, True, True, True, False, True, True, False],
            "all + charge": [False, True, True, True, True, True, True, True],
        }
        for name, decisions in expected.items():
            with self.subTest(cut=name):
                self.assertEqual(ak.to_list(obj_cut_defs["dsaMuons"][name]({}, collection)), [decisions, []])
        # Eta is still an ordinary, unwrapped distance component.
        first = ak.Array([{"outerEta": .3, "outerPhi": math.pi - .02}])
        second = ak.Array([{"outerEta": -.1, "outerPhi": close}])
        self.assert_array_close(dR_outer(first, second), [math.hypot(.4, .04)])
        second = ak.with_field(second, second.outerPhi + 6 * math.pi, "outerPhi")
        self.assert_array_close(dR_outer(first, second), [math.hypot(.4, .04)])

    def test_existing_wrapped_observables_and_cosmic_invariance(self):
        objects = self.lj_fixture()
        mask = ak.num(objects) > 1
        definition = hist_defs["genAs_genAs_absdphi"]
        self.assert_array_close(definition.axes[0].fill_func({"genAs": objects}, mask),
                                [.04, 2.5, math.pi, 0.0, .07])
        self.assertEqual(ak.to_list(evt_cut_defs["dPhi(Mu_0, Mu_1) > 2.4"]({"muons": objects})),
                         [False, True, False, False, True, False, False])
        first, second = ak.unzip(get_pairs(objects))
        wrapped_cosmic = np.sqrt((first.eta + second.eta)**2
                                + (math.pi - angular_distance(first.phi, second.phi))**2)
        for name, key in (("dsaMu_dsaMu_deltaR_cosmic", "dsaMuons"),
                          ("muon_muon_deltaR_cosmic", "muons"),
                          ("Allmuon_Allmuon_deltaR_cosmic", "allMuons")):
            with self.subTest(histogram=name):
                actual = hist_defs[name].axes[0].fill_func({key: objects}, slice(None))
                self.assert_array_close(ak.flatten(actual), ak.flatten(wrapped_cosmic))


if __name__ == "__main__":
    unittest.main()
