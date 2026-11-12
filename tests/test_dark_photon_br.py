"""Dark photon branching fractions and the 4Mu/2Mu2E final-state factors.

Run from the repository root:
    python -m pytest tests/test_dark_photon_br.py -v
"""

import numpy as np
import pytest

from sidm.tools.dark_photon_br import (CHANNELS, M_E, M_MU, M_TAU, _r_table, branching_fractions,
                                       final_state_factor, lepton_width_factor, r_hadronic,
                                       r_perturbative)


def test_closed_form_below_hadron_threshold():
    """At 0.25 GeV only ee and mumu are open, so BR(mumu) follows from the lepton factors."""
    br = branching_fractions(0.25)
    assert br["mumu"] == pytest.approx(0.420, abs=0.005)
    assert br["hadrons"] == 0.0
    assert br["tautau"] == 0.0


def test_r_table_pinned():
    """The shipped table is the full PDG compilation: its ends, the phi peak, R(1.2) and the
    continuity with the perturbative continuum at the last point."""
    sqrt_s, r = _r_table()
    assert (sqrt_s[0], r[0]) == pytest.approx((0.300, 0.01996))
    assert (sqrt_s[-1], r[-1]) == pytest.approx((1.990, 2.17542))
    assert abs(sqrt_s[np.argmax(r)] - 1.0195) < 0.002
    assert r.max() > 20
    assert r_hadronic(1.2) == pytest.approx(0.927, abs=0.02)
    assert abs(r_hadronic(1.99) / r_perturbative(1.99) - 1) < 0.05


def test_r_table_read_only():
    """Callers cannot corrupt the cached table."""
    for column in _r_table():
        with pytest.raises(ValueError):
            column[0] = 0.0


def test_resonance_region():
    """At 1.2 GeV R is taken from the PDG table, above the rho/omega/phi peaks."""
    assert 0.33 <= branching_fractions(1.2)["mumu"] <= 0.41


def test_perturbative_region():
    """At 5 GeV u, d, s, c and tau are open and R is about 3.55."""
    assert r_hadronic(5.0) == pytest.approx(3.55, abs=0.01)
    assert 0.150 <= branching_fractions(5.0)["mumu"] <= 0.160
    assert branching_fractions(5.0)["tautau"] > 0.1


@pytest.mark.parametrize("m", np.linspace(0.21, 10.0, 20))
def test_fractions_sum_to_one(m):
    """Every channel is non-negative and the four sum to one."""
    br = branching_fractions(m)
    assert set(br) == set(CHANNELS)
    assert sum(br.values()) == pytest.approx(1.0, abs=1e-12)
    assert all(v >= 0 for v in br.values())


def test_vectorised_matches_scalar():
    """The array path gives the scalar results element by element."""
    masses = np.linspace(0.21, 10.0, 20)
    vector = branching_fractions(masses)
    for i, m in enumerate(masses):
        for ch, value in branching_fractions(m).items():
            assert vector[ch][i] == pytest.approx(value, rel=1e-12)


@pytest.mark.parametrize("m_lepton", [M_E, M_MU, M_TAU])
def test_lepton_factor_monotonic(m_lepton):
    """Zero up to 2 m_l, then rising towards 1 (the massless limit)."""
    masses = np.linspace(0.5 * m_lepton, 50 * m_lepton, 2000)
    factor = lepton_width_factor(masses, m_lepton)
    assert np.all(factor[masses <= 2 * m_lepton] == 0)
    assert np.all(np.diff(factor[masses > 2 * m_lepton]) > 0)
    assert np.all(factor <= 1)
    assert lepton_width_factor(1e4 * m_lepton, m_lepton) == pytest.approx(1.0, abs=1e-6)


def test_final_state_closure():
    """Below the hadron threshold the two dark photons give 4mu, 2mu2e or 4e and nothing else."""
    closure = final_state_factor(0.25, "4Mu") + final_state_factor(0.25, "2Mu2E")
    assert closure + branching_fractions(0.25)["ee"] ** 2 == pytest.approx(1.0)


def test_edge_masses():
    """Only ee is open at 0.20 GeV; nothing is open at m = 0; NaN propagates."""
    assert branching_fractions(0.20)["ee"] == 1.0
    assert all(v == 0.0 for v in branching_fractions(0.0).values())
    assert all(np.isnan(v) for v in branching_fractions(float("nan")).values())
    vector = branching_fractions(np.array([0.0, 0.20, np.nan]))
    assert vector["ee"][0] == 0.0 and vector["ee"][1] == 1.0 and np.isnan(vector["ee"][2])


def test_final_state_factors():
    """4Mu = BR(mumu)^2 and 2Mu2E = 2 BR(mumu) BR(ee), both dark photons decaying."""
    br = branching_fractions(1.2)
    assert final_state_factor(1.2, "4Mu") == pytest.approx(br["mumu"] ** 2)
    assert final_state_factor(1.2, "2Mu2E") == pytest.approx(2 * br["mumu"] * br["ee"])
    assert 1 / 45 <= final_state_factor(5.0, "4Mu") <= 1 / 38
    with pytest.raises(ValueError):
        final_state_factor(1.2, "4E")
