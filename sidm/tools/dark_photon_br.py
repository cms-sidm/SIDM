"""Branching fractions of a kinetically mixed dark photon and the per-final-state factors of the
signal samples.

A dark photon Zd of mass m couples to the electromagnetic current with strength eps*e, so each
partial width is the photon's coupling times eps^2. Up to the common factor alpha_EM eps^2 m / 3,
which cancels in every branching fraction:

- leptons, l = e, mu, tau: Gamma_ll = sqrt(1 - 4 ml^2/m^2) (1 + 2 ml^2/m^2), zero below 2 ml;
- hadrons: Gamma_had = Gamma_mumu(m) R(m), with R = sigma(e+e- -> hadrons)/sigma(e+e- -> mu+mu-)
  at sqrt(s) = m, zero below the two-pion threshold.

See Batell, Pospelov, Ritz, Phys. Rev. D 79, 115008 (2009), arXiv:0903.0363, and the data-driven
treatment of Ilten, Soreq, Williams, Xue, JHEP 06 (2018) 004, arXiv:1801.04847 (Darkcast).

R(m) comes from ``sidm/configs/r_hadronic.dat``: the full PDG 2024 compilation of the sum of
exclusive channels, 874 points over 0.300-1.990 GeV (provenance in the file header), interpolated
linearly in sqrt(s). Between the two-pion threshold and the first data point R rises linearly from
zero; pi0 gamma below the two-pion threshold is neglected (R of order 1e-3 or less there). The PDG R
is normalised to the pointlike lowest-order sigma(ee -> mumu), while Gamma_had here multiplies it by
Gamma_mumu(m), which carries the muon phase-space factor, as Darkcast does; the difference is that
factor, 1.3% at 0.5 GeV, 4e-4 at 1.2 GeV and nil at 5 GeV. Above the last
data point R is perturbative, R = 3 sum_q e_q^2 (1 + alpha_s/pi) with alpha_s = 0.2 fixed, quarks
open by hadron threshold: u, d, s always, c above 2 m_D0, b above 2 m_B0. That continuum carries
no charmonium resonance structure (psi(3770), psi(4040), ...) and no narrow J/psi or psi(2S).

The signal samples force both dark photons of the bound state to decay, to muons (4Mu) or one to
muons and one to electrons (2Mu2E), so their yield per production cross section carries
BR(mumu)^2 and 2 BR(mumu) BR(ee) respectively.

    python -m sidm.tools.dark_photon_br     # table at 0.25, 1.2, 5 GeV
"""
import functools
import os

import numpy as np

from sidm import BASE_DIR

# PDG 2024 masses [GeV]
M_E = 0.51099895e-3
M_MU = 0.1056583755
M_TAU = 1.77693
M_PI_CHARGED = 0.13957039
M_D0 = 1.86484
M_B0 = 5.27966

LEPTON_MASSES = {"ee": M_E, "mumu": M_MU, "tautau": M_TAU}
HADRON_THRESHOLD = 2 * M_PI_CHARGED
ALPHA_S = 0.2
R_TABLE = os.path.join(BASE_DIR, "configs", "r_hadronic.dat")

# (threshold [GeV], charge) per quark for the perturbative continuum
QUARKS = ((0.0, 2 / 3), (0.0, -1 / 3), (0.0, -1 / 3), (2 * M_D0, 2 / 3), (2 * M_B0, -1 / 3))
CHANNELS = ("ee", "mumu", "tautau", "hadrons")
FINAL_STATES = ("4Mu", "2Mu2E")


def _as_output(values, m):
    """Return a float for scalar mass input, the array otherwise."""
    return float(values) if np.ndim(m) == 0 else values


@functools.lru_cache(maxsize=None)
def _r_table():
    """The (sqrt(s), R) table, read once."""
    sqrt_s, r = np.loadtxt(R_TABLE, comments="#", unpack=True)
    sqrt_s.setflags(write=False)
    r.setflags(write=False)
    return sqrt_s, r


def lepton_width_factor(m, m_lepton):
    """Gamma(Zd -> l l) / (alpha_EM eps^2 m / 3); zero at and below threshold, NaN for NaN m."""
    m = np.asarray(m, dtype=float)
    open_ = m > 2 * m_lepton
    x = np.divide(m_lepton**2, m**2, out=np.zeros_like(m), where=open_)
    factor = np.where(open_, np.sqrt(1 - 4 * x) * (1 + 2 * x), 0.0)
    return _as_output(np.where(np.isnan(m), np.nan, factor), m)


def r_perturbative(m):
    """R = 3 sum_q e_q^2 (1 + alpha_s/pi) over the quarks open at m."""
    m = np.asarray(m, dtype=float)
    charge2 = sum(np.where(m > threshold, q**2, 0.0) for threshold, q in QUARKS)
    return _as_output(3 * charge2 * (1 + ALPHA_S / np.pi), m)


def r_hadronic(m):
    """R(m): PDG data within the table, perturbative above it, zero below 2 m_pi, NaN for NaN m."""
    m = np.asarray(m, dtype=float)
    sqrt_s, r = _r_table()
    data = np.interp(m, np.concatenate(([HADRON_THRESHOLD], sqrt_s)), np.concatenate(([0.0], r)))
    out = np.where(m > sqrt_s[-1], r_perturbative(m), data)
    out = np.where(m > HADRON_THRESHOLD, out, 0.0)
    return _as_output(np.where(np.isnan(m), np.nan, out), m)


def branching_fractions(m):
    """Dark photon branching fractions at mass m [GeV]: dict over ee, mumu, tautau, hadrons.

    All fractions are zero where no channel is open (m <= 2 m_e) and NaN for NaN m.
    """
    widths = {ch: np.asarray(lepton_width_factor(m, ml)) for ch, ml in LEPTON_MASSES.items()}
    widths["hadrons"] = widths["mumu"] * np.asarray(r_hadronic(m))
    total = sum(widths.values())
    out = {}
    for ch, w in widths.items():
        frac = np.divide(w, total, out=np.zeros_like(w), where=total > 0)
        out[ch] = _as_output(np.where(np.isnan(total), np.nan, frac), m)
    return out


def final_state_factor(m, final_state):
    """Fraction of bound-state decays (two dark photons) giving the sample's final state."""
    br = branching_fractions(m)
    if final_state == "4Mu":
        return br["mumu"] ** 2
    if final_state == "2Mu2E":
        return 2 * br["mumu"] * br["ee"]
    raise ValueError(f"unknown final state {final_state!r}; expected one of {FINAL_STATES}")


def main():
    """Print the branching fractions and final-state factors at the sample dark photon masses."""
    print(f"{'m [GeV]':>8} {'R':>6} " + " ".join(f"{c:>8}" for c in CHANNELS + FINAL_STATES))
    for m in (0.25, 1.2, 5.0):
        br = branching_fractions(m)
        cols = [br[c] for c in CHANNELS] + [final_state_factor(m, fs) for fs in FINAL_STATES]
        print(f"{m:8.2f} {r_hadronic(m):6.3f} " + " ".join(f"{v:8.4f}" for v in cols))


if __name__ == "__main__":
    main()
