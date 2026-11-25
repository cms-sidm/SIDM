"""Synthetic checks of the step-3 helpers of the systematics study (sidm/studies/systematics):
the passing fraction with scale factors, the merging of sparse bins of the DSA scale factor, and
the size of a per-lepton-jet systematic in each channel; no files or dataset access."""

import unittest

import hist
import numpy as np

from sidm.studies.systematics import systematics_study as sysstudy

SPEC = {"mass": "m", "peak": (2.9, 3.3), "sidebands": ((2.6, 2.9), (3.3, 3.6))}


def probe_hist():
    """A probe histogram with a mass, a passing and an sf axis: in the peak, 80 passing and 20
    failing probes; the passing ones carry a scale factor of 0.5 in the category nominal"""
    h = hist.Hist(hist.axis.Variable([2.6, 2.9, 3.3, 3.6], name="m"),
                  hist.axis.Integer(0, 2, name="passing", underflow=False, overflow=False),
                  hist.axis.StrCategory(["none", "nominal"], name="sf"), storage=hist.storage.Weight())
    h.fill(m=np.full(80, 3.1), passing=np.ones(80, dtype=int), sf="none")
    h.fill(m=np.full(20, 3.1), passing=np.zeros(20, dtype=int), sf="none")
    h.fill(m=np.full(80, 3.1), passing=np.ones(80, dtype=int), sf="nominal", weight=0.5)
    h.fill(m=np.full(20, 3.1), passing=np.zeros(20, dtype=int), sf="nominal")
    return h


class SystematicsHelpersTest(unittest.TestCase):
    """corrected_efficiency, merge_groups, flat_factor and acts_in"""

    def test_scale_factor_changes_passing_not_probes(self):
        """With the corrections the passing fraction is 0.8 times the scale factor, 0.4, not the
        0.67 of the weighted passing over the weighted total"""
        r = sysstudy.corrected_efficiency(probe_hist(), SPEC, None, "nominal")
        np.testing.assert_allclose(r["eff"], [0.4])
        np.testing.assert_allclose(r["probes"], [100.0])
        r = sysstudy.corrected_efficiency(probe_hist(), SPEC, None, "none")
        np.testing.assert_allclose(r["eff"], [0.8])

    def test_sparse_bins_join_their_neighbours(self):
        """A bin with fewer probes than the minimum joins the bins above it; a short group at
        the top joins the group below"""
        self.assertEqual(sysstudy.merge_groups(np.array([50, 5000, 3000, 100]), 2000),
                         [[0, 1], [2, 3]])
        self.assertEqual(sysstudy.merge_groups(np.array([10, 20]), 2000), [[0, 1]])

    def test_per_lj_sizes_follow_the_lepton_jets(self):
        """A per-lepton-jet entry counts once in 2mu2e and twice in 4mu for the muon type, and
        not at all in 4mu for the egamma type"""
        mu = {"per_lj": 0.1, "lj_type": "muon"}
        egm = {"per_lj": 0.1, "lj_type": "egamma"}
        self.assertAlmostEqual(sysstudy.flat_factor(mu, "SR_2mu2e"), 1.1)
        self.assertAlmostEqual(sysstudy.flat_factor(mu, "SR_4mu"), 1.21)
        self.assertAlmostEqual(sysstudy.flat_factor(egm, "SR_4mu"), 1.0)
        self.assertAlmostEqual(sysstudy.flat_factor({"flat": 0.0083}, "SR_4mu"), 1.0083)

    def test_per_lj_entry_without_lepton_jets_does_not_act(self):
        """The egamma-type entry cannot change the 4mu signal, so the fit checks skip it there;
        every other entry acts in both channels"""
        egm = {"per_lj": 0.044, "lj_type": "egamma"}
        self.assertFalse(sysstudy.acts_in(egm, "SR_4mu"))
        self.assertTrue(sysstudy.acts_in(egm, "SR_2mu2e"))
        self.assertTrue(sysstudy.acts_in({"per_lj": 0.017, "lj_type": "muon"}, "SR_4mu"))
        self.assertTrue(sysstudy.acts_in({"syst": "pileup"}, "SR_4mu"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
