"""Synthetic muon-pair DeltaR naming and filling regressions; no event files."""

import copy
import io
import math
import unittest
import warnings
from contextlib import redirect_stdout

import awkward as ak
import hist
import numpy as np
from coffea.nanoevents.methods import vector

from sidm.definitions.hists import hist_defs
from sidm.tools.sidm_processor import SidmProcessor


NAME = "muon_muon_dr"
OLD_NAME = "muon_muon_dphi"


def muons(etas, phis):
    """Vectors with increasing pt, so input-order tests also detect sorting."""
    phi = ak.Array(phis)
    return ak.zip(
        {"pt": 10.0 + 10.0 * ak.local_index(phi, axis=1),
         "eta": ak.Array(etas), "phi": phi, "mass": ak.ones_like(phi)},
        with_name="PtEtaPhiMLorentzVector", behavior=vector.behavior,
    )


class MuonMuonDrTests(unittest.TestCase):
    def fill_hist(self, objects, weights, definition=None):
        """Use production booking/filling and fail on its skip-warning path."""
        if definition is None:
            definition = copy.deepcopy(hist_defs[NAME])
            definition.make_hist(NAME)
        output = io.StringIO()
        with warnings.catch_warnings(record=True) as caught, redirect_stdout(output):
            warnings.simplefilter("always")
            definition.fill(objects, ak.Array(weights), verbose=True)
        self.assertEqual(output.getvalue(), "", "Histogram filling must not skip")
        self.assertEqual(caught, [], "Histogram filling must not warn")
        return definition.hist

    def fixture(self):
        return muons(
            [[0.0, .8, .8], [0.0, .8], [0.0, .03], [0.0, 7.0], [], [.1]],
            [[0.0, .6, .6], [0.0, .6],
             [math.pi - .02, -math.pi + .02], [0.0, 0.0], [], [.2]],
        )

    def test_delta_r_wrapped_phi_first_two_and_mask(self):
        collection = self.fixture()
        objects = {"muons": collection}
        definition = hist_defs[NAME]
        mask = definition.evt_mask(objects)
        self.assertEqual(ak.to_list(mask), [True, True, True, True, False, False])
        actual = definition.axes[0].fill_func(objects, mask)
        # Independent references: hypot(.8, .6) = 1 and hypot(.03, .04) = .05.
        np.testing.assert_allclose(ak.to_numpy(actual), [1.0, 1.0, .05, 7.0],
                                   atol=2e-14, rtol=1e-13)
        self.assertNotAlmostEqual(float(actual[0]), .6)  # abs(DeltaPhi)
        # The highest-pt third muon coincides with the second: sorting or choosing
        # the closest/all pairs would change the first event's result/count.
        self.assertEqual(len(actual), 4)
        sorted_objects = {"muons": collection[ak.argsort(collection.pt, ascending=False)]}
        sorted_values = definition.axes[0].fill_func(sorted_objects, mask)
        self.assertAlmostEqual(float(sorted_values[0]), 0.0)

    def test_weighted_fill_named_projection_and_variance(self):
        filled = self.fill_hist({"muons": self.fixture()}, [2, -1, -3, 5, 101, 103])
        self.assertIs(filled.storage_type, hist.storage.Weight)
        projected = filled.project(NAME)
        expected_values = np.zeros(50)
        expected_variances = np.zeros(50)
        expected_values[projected.axes[0].index(1.0)] = 1.0
        expected_variances[projected.axes[0].index(1.0)] = 5.0
        expected_values[projected.axes[0].index(.05)] = -3.0
        expected_variances[projected.axes[0].index(.05)] = 9.0
        np.testing.assert_array_equal(projected.values(), expected_values)
        np.testing.assert_array_equal(projected.variances(), expected_variances)
        self.assertEqual(projected.sum().value, -2.0)
        self.assertEqual(projected.sum().variance, 14.0)
        self.assertEqual(projected.sum(flow=True).value, 3.0)
        self.assertEqual(projected.sum(flow=True).variance, 39.0)
        self.assertEqual(projected.values(flow=True)[-1], 5.0)

    def test_empty_events_and_zero_event_input(self):
        collection = self.fixture()
        for empty in (collection[:0], collection[:, :0], collection[:, :1]):
            with self.subTest(type=str(ak.type(empty))):
                objects = {"muons": empty}
                definition = hist_defs[NAME]
                values = definition.axes[0].fill_func(objects, definition.evt_mask(objects))
                self.assertEqual(ak.to_list(values), [])
                filled = self.fill_hist(objects, np.ones(len(empty)))
                self.assertEqual(filled.sum(flow=True).value, 0.0)
                self.assertEqual(filled.sum(flow=True).variance, 0.0)

    def test_preserved_regular_axis_range_and_flow(self):
        collection = muons([[0.0, value] for value in (0.0, 6.279, 6.28, 7.0)],
                           [[0.0, 0.0]] * 4)
        filled = self.fill_hist({"muons": collection}, [1, 2, -3, 4])
        axis = filled.axes[NAME]
        self.assertEqual(axis.size, 50)
        np.testing.assert_allclose(axis.edges, np.linspace(0.0, 6.28, 51))
        self.assertTrue(axis.traits.underflow)
        self.assertTrue(axis.traits.overflow)
        self.assertEqual(axis.index(0.0), 0)
        self.assertEqual(axis.index(6.28), 50)
        self.assertEqual(filled.values()[0], 1.0)
        self.assertEqual(filled.values()[-1], 2.0)
        self.assertEqual(filled.values(flow=True)[0], 0.0)
        self.assertEqual(filled.values(flow=True)[-1], 1.0)
        self.assertEqual(filled.variances(flow=True)[-1], 25.0)

    def test_ttjets_config_resolves_and_books_new_name(self):
        self.assertIn(NAME, hist_defs)
        self.assertNotIn(OLD_NAME, hist_defs)
        self.assertEqual(hist_defs[NAME].axes[0].name, NAME)
        self.assertEqual(hist_defs[NAME].axes[0].axis.name, NAME)
        self.assertEqual(hist_defs[NAME].axes[0].axis.label,
                         r"$\Delta R$($\mu_0$, $\mu_1$)")
        # This resolves the real YAML collection and books all its members;
        # no event-processing entry point or dataset is used.
        booked = SidmProcessor(["base"], ["ttjets_base"]).build_histograms()
        self.assertIn(NAME, booked)
        self.assertNotIn(OLD_NAME, booked)
        definition = booked[NAME]
        self.assertEqual(tuple(definition.hist.axes.name), ("channel", NAME))
        filled = self.fill_hist(
            {"muons": muons([[0.0, .8]], [[0.0, .6]]), "ch": "base"},
            [2.5], definition=definition,
        )
        projected = filled[{"channel": "base"}].project(NAME)
        self.assertEqual(projected.values()[projected.axes[0].index(1.0)], 2.5)
        self.assertEqual(projected.sum().value, 2.5)
        self.assertEqual(projected.sum().variance, 6.25)


if __name__ == "__main__":
    unittest.main()
