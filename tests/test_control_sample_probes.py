"""Synthetic checks of the probe helpers of the lepton-jet control samples, version 3, and of the
DSA-muon study (sidm/studies/systematics/lj_reco_efficiency.md); no files or dataset access."""

import unittest

import awkward as ak
import numpy as np
from coffea.nanoevents.methods import vector

from sidm.definitions.hists import combine_sfs, jpsi_dsa_count, lj_composition, nearest_factors
from sidm.definitions.objects import dsa_probe_muons, jpsi_probe_vertex_lxy


def candidates(events):
    """Lorentz-vector records of literal fields, one list per event"""
    return ak.with_name(ak.Array(events), "PtEtaPhiMLorentzVector", behavior=vector.behavior)


def muon(idx, pt=45.0, phi=0.0, charge=1, iso=0.05):
    """A reconstructed tight, prompt muon at eta = 0 with the fields the probe definitions read"""
    return {"idx": idx, "pt": pt, "eta": 0.0, "phi": phi, "mass": 0.105658, "charge": charge,
            "tightId": True, "isTracker": True, "pfRelIso04_all": iso, "dxy": 0.001, "dz": 0.01}


def point(idx, phi):
    """A massive candidate at eta = 0 with an index"""
    return {"idx": idx, "pt": 10.0, "eta": 0.0, "phi": phi, "mass": 0.1}


class ControlSampleProbesTest(unittest.TestCase):
    """The DSA count of J/psi probes, the dimuon-vertex lookup, the Z -> mumu DSA probes, the scale
    factors of passing probes and the make-up of their lepton jets"""

    def test_dsa_count_counts_distinct_matches(self):
        """Two DSA muons near the two muons give 2, one DSA muon nearest to both gives 1"""
        pair = {"mu0": point(0, 0.0), "mu1": point(1, 0.05)}
        probes = ak.Array([[pair], [pair], [pair]])
        probes = ak.zip({"mu0": candidates(probes.mu0), "mu1": candidates(probes.mu1)},
                        depth_limit=2)
        dsa = candidates([[point(0, 0.0), point(1, 0.05)], [point(0, 0.02)], []])
        self.assertEqual(ak.to_list(jpsi_dsa_count(probes, dsa)), [[2], [1], [0]])

    def test_vertex_lxy_matches_both_legs_in_either_order(self):
        """The PF-PF vertex of the two probe muons gives its distance from the PV; a mixed-leg
        or invalid vertex, or none, gives -1"""
        probes = ak.zip({"mu0": ak.Array([[{"idx": 3}], [{"idx": 1}]]),
                         "mu1": ak.Array([[{"idx": 5}], [{"idx": 2}]])}, depth_limit=2)
        row = {"originalMuonIdx1": 5.0, "originalMuonIdx2": 3.0, "isDSAMuon1": 0.0,
               "isDSAMuon2": 0.0, "isValid": 1.0, "vx": 0.3, "vy": 0.4}
        vertices = ak.Array([[{**row, "vx": 9.0, "isDSAMuon1": 1.0}, row],
                             [{**row, "originalMuonIdx1": 1.0, "originalMuonIdx2": 2.0,
                               "isValid": 0.0}]])
        pv = ak.Array([[{"x": 0.0, "y": 0.0}], [{"x": 0.0, "y": 0.0}]])  # one per event, as in the processor
        np.testing.assert_allclose(ak.to_list(ak.flatten(jpsi_probe_vertex_lxy(probes, vertices, pv))),
                                   [0.5, -1.0])

    def test_z_probes_need_a_tag_of_opposite_charge(self):
        """Two isolated tight muons of opposite charge at m = 90 GeV are each a tag and a probe;
        a non-isolated muon is neither"""
        muons = candidates([[muon(0), muon(1, phi=np.pi, charge=-1),
                             muon(2, phi=np.pi, charge=-1, iso=0.5)]])
        probes = dsa_probe_muons({"muonsRaw": muons})
        self.assertEqual(sorted(ak.to_list(probes.idx[0])), [0, 1])
        np.testing.assert_allclose(ak.to_list(probes.m_tp[0]), [90.0, 90.0], rtol=1e-4)

    def test_combine_sfs_shifts_one_correction_at_a_time(self):
        """nominal multiplies every correction's nominal factor, <name>Up shifts that one only"""
        ones = ak.Array([[1.0, 1.0]])
        factors = combine_sfs({"a": {"nominal": ak.Array([[0.9, 0.8]]), "up": ak.Array([[1.0, 0.9]]),
                                     "down": ak.Array([[0.8, 0.7]])},
                               "b": {"nominal": ak.Array([[0.5, 1.0]]), "up": ak.Array([[0.6, 1.0]]),
                                     "down": ak.Array([[0.4, 1.0]])}}, ones)
        self.assertEqual(list(factors), ["none", "nominal", "aUp", "aDown", "bUp", "bDown"])
        np.testing.assert_allclose(ak.to_list(factors["nominal"][0]), [0.45, 0.8])
        np.testing.assert_allclose(ak.to_list(factors["aUp"][0]), [0.5, 0.9])
        np.testing.assert_allclose(ak.to_list(factors["bDown"][0]), [0.36, 0.8])
        np.testing.assert_allclose(ak.to_list(factors["none"][0]), [1.0, 1.0])

    def test_nearest_factors_go_to_passing_probes_only(self):
        """A probe within dR < 0.2 of an object takes its factor, a probe without one gets 1"""
        probes = candidates([[point(0, 0.0), point(1, 1.0)], [point(0, 0.0)]])
        others = candidates([[point(5, 0.05), point(6, 2.0)], []])
        factors = nearest_factors(probes, others, {"nominal": ak.Array([[0.8, 0.5], []])})
        self.assertEqual(ak.to_list(factors["nominal"]), [[0.8, 1.0], [1.0]])

    def test_lj_composition_counts_the_matched_lepton_jet(self):
        """A probe takes the photon and electron counts of the lepton jet within dR < 0.2, a
        probe without one gets -1"""
        def lepton_jet(phi, n_photons, n_electrons):
            return {**point(0, phi), "photons": [{"pt": 30.0}] * n_photons,
                    "electrons": [{"pt": 30.0}] * n_electrons}
        probes = candidates([[point(0, 0.0), point(1, 1.0)], [point(0, 0.0)]])
        ljs = candidates([[lepton_jet(0.05, 0, 1), lepton_jet(2.0, 1, 0)], []])
        photons, electrons = lj_composition(probes, ljs)
        self.assertEqual(ak.to_list(photons), [[0, -1], [-1]])
        self.assertEqual(ak.to_list(electrons), [[1, -1], [-1]])


if __name__ == "__main__":
    unittest.main(verbosity=2)
