"""Synthetic filling and output-schema checks for the back-to-back DSA phi axis."""

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
from sidm.tools.histogram import Axis


PHI_KEY = "back_to_back_dsa_pair_phi_diff"
PHI_AXIS = "back to back dsa_pair_phi_diff"
ETA_KEY = "back_to_back_dsa_pair_eta_diff"
ETA_AXIS = "back to back dsa_pair_eta_diff"
PHI_LABEL = r"Back-to-back $|\Delta\phi|$(DSA $\mu_1$, DSA $\mu_2$)"


def dsa_fixture():
    """Pairs with wrapped phi, both pi signs, 3D rejection and empty events."""
    phis = ak.Array([
        [math.pi - .02, -math.pi + .02],  # close across the branch cut: reject
        [.12, -math.pi + .02],           # raw difference > pi: keep, wrap
        [0.0, math.pi],
        [0.0, .1, math.pi],             # two selected pairs, one event weight each
        [0.0, math.pi],                 # opposite phi but same positive eta: reject
        [0.0, math.pi],                 # opposite phi and opposite eta: keep
        [0.0, math.pi - .2],            # fails cos(3D angle) <= -0.99
        [], [.3], [0.0, -math.pi],
    ])
    etas = ak.Array([
        [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0, 0.0],
        [.5, .5], [.4, -.4], [0.0, 0.0], [], [0.0], [0.0, 0.0],
    ])
    dsa = ak.zip(
        {"pt": 60.0 - 10.0 * ak.local_index(phis, axis=1), "eta": etas,
         "phi": phis, "mass": ak.ones_like(phis)},
        with_name="PtEtaPhiMLorentzVector", behavior=vector.behavior,
    )
    return dsa, ak.Array([101, 2, -3, 5, 103, 7, 107, 109, 113, 11])


class DsaPhiAxisTests(unittest.TestCase):
    def fill_hist(self, definition, key, dsa, weights):
        """Exercise keyword filling with the runtime category axes present."""
        definition.make_hist(key, channels=["synthetic"], lj_reco_choices=["test"])
        output = io.StringIO()
        with redirect_stdout(output), warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            definition.fill({"dsaMuons": dsa, "ch": "synthetic", "lj_reco": "test"},
                            weights, verbose=True)
        self.assertEqual(output.getvalue(), "", "Histogram filling must not silently skip")
        self.assertEqual(caught, [], "Histogram filling must not emit warnings")
        return definition.hist

    def test_named_projection_weighted_bins_and_exact_pi(self):
        dsa, weights = dsa_fixture()
        definition = copy.deepcopy(hist_defs[PHI_KEY])
        self.assertEqual(definition.axes[0].name, PHI_AXIS)
        self.assertEqual(definition.storage, "weight")
        filled = self.fill_hist(definition, PHI_KEY, dsa, weights)
        self.assertEqual(definition.name, PHI_KEY)
        self.assertEqual(filled.axes.name, ("lj_reco", "channel", PHI_AXIS))
        axis = filled.axes[PHI_AXIS]
        self.assertEqual(axis.label, PHI_LABEL)
        self.assertEqual(axis.size, 100)
        self.assertEqual(axis.edges[0], 0.0)
        self.assertEqual(axis.edges[-1], np.nextafter(math.pi, math.inf))
        self.assertTrue(axis.traits.underflow)
        self.assertTrue(axis.traits.overflow)
        self.assertEqual(axis.index(math.pi), 99)

        projected = filled.project(PHI_AXIS)
        # The six selected pairs carry weights 2, -3, 5, 5, 7 and 11.
        values, variances = np.zeros(102), np.zeros(102)
        values[97], variances[97] = 7, 29      # pi - 0.1, including flow offset
        values[100], variances[100] = 20, 204  # exact pi in the final normal bin
        np.testing.assert_array_equal(projected.values(flow=True), values)
        np.testing.assert_array_equal(projected.variances(flow=True), variances)
        self.assertEqual(projected.sum().value, 27)
        self.assertEqual(projected.sum().variance, 233)
        self.assertEqual(filled[{PHI_AXIS: hist.loc(math.pi)}].sum().value, 20)
        with self.assertRaises(ValueError):
            filled.project(ETA_AXIS)

    def test_pair_selection_and_empty_fills(self):
        dsa, weights = dsa_fixture()
        definition = hist_defs[PHI_KEY]
        objs = {"dsaMuons": dsa}
        values = definition.axes[0].fill_func(objs, definition.evt_mask(objs))
        self.assertEqual(ak.to_list(ak.num(values)), [0, 1, 1, 2, 0, 1, 0, 0, 0, 1])
        np.testing.assert_allclose(ak.to_numpy(ak.flatten(values)),
                                   [math.pi - .1, math.pi, math.pi, math.pi - .1,
                                    math.pi, math.pi], atol=2e-14, rtol=0)
        for empty, empty_weights in ((dsa[:0], weights[:0]),
                                     (dsa[:, :0], weights), (dsa[:, :1], weights)):
            with self.subTest(counts=ak.to_list(ak.num(empty))):
                filled = self.fill_hist(copy.deepcopy(definition), PHI_KEY,
                                        empty, empty_weights)
                self.assertEqual(filled.sum(flow=True).value, 0)
                self.assertEqual(filled.sum(flow=True).variance, 0)

    def test_axis_rename_preserves_contents_and_variances(self):
        dsa, weights = dsa_fixture()
        current = copy.deepcopy(hist_defs[PHI_KEY])
        legacy = copy.deepcopy(hist_defs[PHI_KEY])
        # Recreate the former schema only; the value function and storage are shared.
        legacy.axes[0] = Axis(hist.axis.Regular(
            100, 0, np.nextafter(math.pi, math.inf), name=ETA_AXIS, label=PHI_LABEL),
            legacy.axes[0].fill_func)
        before = self.fill_hist(legacy, PHI_KEY, dsa, weights).project(ETA_AXIS)
        after = self.fill_hist(current, PHI_KEY, dsa, weights).project(PHI_AXIS)
        np.testing.assert_array_equal(before.axes[0].edges, after.axes[0].edges)
        self.assertEqual(before.axes[0].traits, after.axes[0].traits)
        self.assertEqual(before.axes[0].label, after.axes[0].label)
        np.testing.assert_array_equal(before.values(flow=True), after.values(flow=True))
        np.testing.assert_array_equal(before.variances(flow=True), after.variances(flow=True))

    def test_real_eta_histogram_retains_its_name_and_values(self):
        dsa, weights = dsa_fixture()
        definition = copy.deepcopy(hist_defs[ETA_KEY])
        self.assertEqual(definition.axes[0].name, ETA_AXIS)
        self.assertEqual(definition.storage, "weight")
        filled = self.fill_hist(definition, ETA_KEY, dsa, weights).project(ETA_AXIS)
        self.assertEqual(filled.axes[0].size, 100)
        self.assertEqual(filled.axes[0].edges[0], 0.0)
        self.assertEqual(filled.axes[0].edges[-1], 6.0)
        self.assertEqual(filled.axes[0].label,
                         r" back to back abs(DSA $\mu$1 $\eta$ - DSA $\mu$1 $\eta$ ) ")
        values, variances = np.zeros(102), np.zeros(102)
        values[1], variances[1] = 20, 184      # eta separation zero
        values[14], variances[14] = 7, 49     # eta separation 0.8
        np.testing.assert_array_equal(filled.values(flow=True), values)
        np.testing.assert_array_equal(filled.variances(flow=True), variances)


if __name__ == "__main__":
    unittest.main()
