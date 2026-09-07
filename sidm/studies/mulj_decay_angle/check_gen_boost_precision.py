"""How the single-precision boost of the generated histograms reads one simulated file.

Reads, over xrootd, one LLPNanoAOD file (part-0) of the 4mu signal at m_Zd = 0.25 GeV for
each bound-state mass given (default 500 and 1000 GeV), applies no selection, and builds
the generator-level muons of each dark photon with the processor's own object definitions.
For each file it prints:

* the number of dark photons with two muons, and the 1 and 99 percent points of their
  dimuon mass, computed in double precision from the muons' px, py and pz;
* the median Lorentz factor of those dark photons, their energy over the nominal m_Zd, the
  definition the toy of section b uses;
* |cos theta*| as the generated histogram computes it (cos_theta_in_parent_frame, single
  precision): alpha of 1 + alpha cos^2 theta* fitted to 25 bins, hi/lo, and the number of
  muons whose cosine comes out as not a number or at one or above;
* the same muons boosted in double precision (toy_abs_cos_theta of _mulj_lib);
* the boost-free reading |cos theta*| = |E1 - E2| / (beta* |p1 + p2|), one entry per dark
  photon, with beta* = sqrt(1 - 4 m_mu^2 / m^2) taken from each dark photon's own
  generated mass m, and the alpha that mixture of beta* values predicts:
  mean(a / (1 + a/3)) / mean(1 / (1 + a/3)) with a = beta*^2 / (2 - beta*^2).

The histogram fills both muons of a dark photon, which carry the same |cos theta*|, so the
alpha errors of the two boosted readings are those of one entry per dark photon: the fit
error over all muons times sqrt(2).

Run from this directory with the repository on PYTHONPATH:
    python check_gen_boost_precision.py [--m-xx 500 1000]
"""
import argparse

import awkward as ak
import numpy as np
from coffea.nanoevents import NanoAODSchema, NanoEventsFactory

import _mulj_lib as lib
from sidm.definitions.hists import cos_theta_in_parent_frame
from sidm.definitions.objects import postLj_objs_MC, preLj_objs

MU_MASS = 0.105658
M_ZD = 0.25
BASE = ("root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/ULSignalSamples/2018_v10/"
        "BsTo2DpTo4Mu/")
CTAU = {500: "0p004", 1000: "0p002"}
# the toy's fit and hi/lo, so that these numbers are computed exactly as the toy table's
fit_alpha, hi_lo = lib._toy_fit_alpha, lib._toy_hi_lo  # pylint: disable=protected-access


def file_url(m_xx):
    """The part-0 file of the 4mu sample at this bound-state mass and m_Zd = 0.25 GeV."""
    name = f"CutDecayFalse_SIDM_BsTo2DpTo4Mu_MBs-{m_xx}_MDp-0p25_ctau-{CTAU[m_xx]}_v3"
    return f"{BASE}{name}/LLPnanoAODv2/{name}_part-0.root"


def four_vectors(records, mass=None):
    """(E, px, py, pz) rows in float64 from the px, py, pz (and mass) of flat records."""
    px, py, pz = (ak.to_numpy(ak.flatten(records[k])).astype(np.float64)
                  for k in ("px", "py", "pz"))
    m = (np.full(len(px), mass) if mass is not None
         else ak.to_numpy(ak.flatten(records["mass"])).astype(np.float64))
    return np.stack([np.sqrt(px * px + py * py + pz * pz + m * m), px, py, pz], axis=-1)


def check(m_xx):
    """Print every number of one file."""
    events = NanoEventsFactory.from_root({file_url(m_xx): "Events"},
                                         schemaclass=NanoAODSchema, mode="virtual").events()
    objs = {"gens": events.GenPart}
    objs["genMus"] = preLj_objs["genMus"](events)
    objs["genAs"] = preLj_objs["genAs"](events)
    objs["genMus_fromA"] = postLj_objs_MC["genMus_fromA"](objs)
    mask = ak.num(objs["genMus_fromA"]) >= 2
    single = ak.to_numpy(ak.flatten(abs(cos_theta_in_parent_frame(
        objs, mask, "genMus_fromA", mass=MU_MASS)))).astype(np.float64)
    muons = objs["genMus_fromA"][mask]
    parents = muons.parent
    double = lib.toy_abs_cos_theta(four_vectors(muons, MU_MASS), four_vectors(parents))

    # the two muons of each dark photon, for the boost-free reading and the dimuon mass
    mother = ak.to_numpy(ak.flatten(muons.genPartIdxMotherG))
    order = np.argsort(mother, kind="stable")
    mu_p4 = four_vectors(muons, MU_MASS)[order]
    first = np.r_[True, mother[order][1:] != mother[order][:-1]]
    starts = np.flatnonzero(first)
    pairs = starts[np.diff(np.r_[starts, len(order)]) == 2]
    mu1, mu2 = mu_p4[pairs], mu_p4[pairs + 1]
    zd_p4 = four_vectors(parents)[order][pairs]
    total = mu1 + mu2
    dimuon_mass = np.sqrt(np.maximum(total[:, 0] ** 2 - (total[:, 1:] ** 2).sum(axis=1), 0.0))
    zd_mass = ak.to_numpy(ak.flatten(parents.mass)).astype(np.float64)[order][pairs]
    beta_star = np.sqrt(np.clip(1.0 - 4.0 * MU_MASS ** 2 / zd_mass ** 2, 1e-12, None))
    energy = np.abs(mu1[:, 0] - mu2[:, 0]) / (beta_star * np.linalg.norm(total[:, 1:], axis=1))
    a = beta_star ** 2 / (2.0 - beta_star ** 2)
    expected = float(np.mean(a / (1.0 + a / 3.0)) / np.mean(1.0 / (1.0 + a / 3.0)))

    print(f"m_XX = {m_xx} GeV, m_Zd = {M_ZD} GeV, no selection: {file_url(m_xx)}")
    print(f"  dark photons with two muons {len(pairs):,}; dimuon mass 1 to 99 percent "
          f"{np.percentile(dimuon_mass, 1):.3f} to {np.percentile(dimuon_mass, 99):.3f} GeV; "
          f"median Lorentz factor E / m_Zd {np.median(zd_p4[:, 0]) / M_ZD:,.0f}")
    for label, values in (("single precision (as the histogram)", single),
                          ("double precision, same muons", double)):
        alpha, error = fit_alpha(values)
        finite = values[np.isfinite(values)]
        print(f"  {label:36s}: alpha {alpha:+.3f} +- {error * np.sqrt(2.0):.3f}, hi/lo "
              f"{hi_lo(values):.2f}, above 0.96 {np.mean(finite > 0.96):.3f}, "
              f"not a number {int(np.sum(~np.isfinite(values))):,} and at one or above "
              f"{int(np.sum(finite >= 1.0)):,} of {len(values):,} muons")
    alpha, error = fit_alpha(energy)
    print(f"  {'boost-free |E1 - E2| / (beta* |p1 + p2|)':36s}: alpha {alpha:+.3f} +- {error:.3f} "
          f"with beta* from each dark photon's generated mass; that mass distribution "
          f"predicts {expected:.3f}")


def main():
    """Parse the bound-state masses and check one file of each."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--m-xx", type=int, nargs="+", default=[500, 1000],
                        choices=sorted(CTAU), help="bound-state masses to check [GeV]")
    for m_xx in parser.parse_args().m_xx:
        check(m_xx)


if __name__ == "__main__":
    main()
