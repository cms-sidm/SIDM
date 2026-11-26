"""Synthetic checks of the step-3 helpers of the systematics study (sidm/studies/systematics):
the passing fraction with scale factors, the merging of sparse bins of the DSA scale factor, the
size of a per-lepton-jet systematic in each channel, and the pixel-veto factor per photon against
an efficiency reweighting of the veto; no files or dataset access."""

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


def seed_hist():
    """A seed histogram (hists.zmmg_probe_seed_hist) of generator-matched FSR probes at 50 GeV that
    pass the other photon cuts: in each of the barrel and the endcap, 90 without a pixel seed whose
    lepton jet holds the photon, and 10 with a seed, 4 of which have a lepton jet of one electron;
    plus, in the barrel, 10 probes of the first kind in each sideband and 10 more in the peak,
    which the sideband subtraction removes"""
    def flag(name):
        return hist.axis.Integer(0, 2, name=name, underflow=False, overflow=False)
    h = hist.Hist(flag("gen_fsr"), hist.axis.Variable([20, 30, 40, 60, 120], name="probe_pt"),
                  hist.axis.Variable([0, 1.4442, 1.566, 2.4], name="probe_abseta"),
                  hist.axis.Variable([0, 0.94, 1.5], name="probe_r9"), flag("probe_pixelSeed"),
                  flag("probe_photon_id"), flag("probe_electron_near"),
                  hist.axis.Variable([75, 85, 95, 105], name="m_mmg"),
                  hist.axis.Integer(-1, 3, name="lj_photons", underflow=False),
                  hist.axis.Integer(-1, 3, name="lj_electrons", underflow=False),
                  storage=hist.storage.Weight())

    def fill(n, abseta, seed, lepton_jet, mass=90.0):
        """n probes; lepton_jet: (photons, electrons) of the matched lepton jet, -1 for none"""
        h.fill(gen_fsr=1, probe_pt=50.0, probe_abseta=abseta, probe_r9=1.0, probe_pixelSeed=seed,
               probe_photon_id=1, probe_electron_near=0, m_mmg=np.full(n, mass),
               lj_photons=lepton_jet[0], lj_electrons=lepton_jet[1])
    for abseta in (0.5, 2.0):
        fill(90, abseta, 0, (1, 0))
        fill(4, abseta, 1, (0, 1))
        fill(6, abseta, 1, (-1, -1))
    for mass in (80.0, 90.0, 100.0):
        fill(10, 0.5, 0, (1, 0), mass)
    return h


class SystematicsHelpersTest(unittest.TestCase):
    """corrected_efficiency, merge_groups, flat_factor, acts_in, seed_outcomes and
    pixel_veto_factors"""

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

    def test_seeded_photons_pass_through_electrons(self):
        """Without a seed every probe has a lepton jet with its photon; with one, 8 of 20 have a
        lepton jet of electrons only and the rest none"""
        h = seed_hist()
        np.testing.assert_allclose(sysstudy.seed_outcomes(h, 0), (1.0, 0.0, 0.0))
        np.testing.assert_allclose(sysstudy.seed_outcomes(h, 1), (0.0, 0.4, 0.6))

    def test_period_spread(self):
        """Two periods that agree give chi2 0; one 2 sigma above the other, each with sigma 1,
        gives pulls of 1 about their mean, so chi2 2 for one degree of freedom per bin"""
        same = sysstudy.period_spread(np.array([[0.85, 0.80], [0.85, 0.80]]), np.full((2, 2), 0.01))
        self.assertEqual(same, (0.0, 2, 0.0))
        chi2, ndf, pull = sysstudy.period_spread(np.array([[1.0], [3.0]]), np.ones((2, 1)))
        self.assertEqual((chi2, ndf, pull), (2.0, 1, 1.0))
        self.assertTrue(sysstudy.in_period("B")("DoubleMuon_2018B_13"))
        self.assertFalse(sysstudy.in_period("B")("DoubleMuon_2018D_3"))

    def test_veto_reweighting_keeps_the_photons_that_pass_through_electrons(self):
        """Barrel, after the subtraction: 90 probes without a seed pass with their photon, 4 of 10
        with a seed pass through an electron. Per photon only the 90 carry the factor 0.8; the
        reweighting also raises the 4 by (1 - 0.8 * 0.9) / (1 - 0.9) = 2.8"""
        f = sysstudy.pixel_veto_factors(seed_hist(), {"barrel": 0.8, "endcap": 0.9})
        np.testing.assert_allclose(f["barrel"], (0.9, 76 / 94, 83.2 / 94))
        np.testing.assert_allclose(f["endcap"], (0.9, 85 / 94, 88.6 / 94))
        np.testing.assert_allclose(f["both"][1:], (161 / 188, 171.8 / 188))


if __name__ == "__main__":
    unittest.main(verbosity=2)
