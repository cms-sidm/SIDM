"""Shared analysis helpers for the dark-photon lifetime / proper-cτ study.

Loads the full-grid coffea and provides the three cτ estimators used by the
companion notebooks:
  - the histogram mean of the proper-lxyz distribution (faithful regime),
  - the Option-B core-slope fit (exponential slope on the un-truncated core),
  - the Option-C acceptance-corrected fit (folds in the lab-decay-length cap).

The grid was produced by _fullgrid_run.py and lives at
root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/coffea_outputs/<user>/lifetime_study/.
"""
import os
import subprocess
from collections import defaultdict

import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import coffea.util

CHANNEL = "baseNoLj_noTrigger"   # gen-level channel, no HLT-trigger sculpting
EOS_COFFEA = (
    "root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/"
    "coffea_outputs/murtazas/lifetime_study/lifetime_fullgrid.coffea"
)


def load_grid(local_path="lifetime_fullgrid.coffea", eos_url=EOS_COFFEA):
    """Load the full-grid coffea; xrdcp it from the shared EOS area if not cached locally."""
    if not os.path.exists(local_path):
        subprocess.run(["xrdcp", "-f", eos_url, local_path], check=True)
    return coffea.util.load(local_path)


def ctau_cm(name):
    """Nominal proper cτ in cm, parsed from the sample-name token (in mm)."""
    return float(name.split("_")[-1].replace("mm", "").replace("p", ".")) / 10.0


def mass_point(name):
    """(channel, BS-mass, Dp-mass) key, e.g. ('4Mu', '500GeV', '5p0GeV')."""
    p = name.split("_")
    return (p[0], p[1], p[2])


def proper_density(output, sample, channel=CHANNEL):
    """Proper-lxyz dN/dx, effective counts (value**2/variance), and the histogram-mean cτ."""
    h = output[sample]["hists"]["genAs_lxyz_proper"][{"channel": channel}]
    edges = h.axes[-1].edges
    centers = h.axes[-1].centers
    vals = h.values()
    var = h.variances()
    effN = np.where(var > 0, vals**2 / var, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        dens = np.where(vals > 0, vals / np.diff(edges), np.nan)
    mean = (vals * centers).sum() / max(vals.sum(), 1)
    return centers, dens, effN, mean


def core_slope_fit(output, sample, channel=CHANNEL, k=3.0):
    """Option B: fit log(dN/dx) = a - x/cτ on the un-truncated low-x core (cτ = -1/slope),
    iteratively trimming the systematic high-x deficit (the lab-cap truncation)."""
    centers, dens, effN, mean = proper_density(output, sample, channel)
    ly = np.where(dens > 0, np.log(dens), np.nan)
    mask = effN >= 5
    slope = np.nan
    if mask.sum() >= 4:
        for _ in range(30):
            slope, ic = np.polyfit(centers[mask], ly[mask], 1, w=np.sqrt(effN[mask]))
            sig = np.where(effN > 0, 1 / np.sqrt(effN), np.inf)
            trim = mask & (ly < ic + slope * centers - k * sig) & (centers > np.median(centers[mask]))
            if not trim.any():
                break
            mask = mask & ~trim
    ctau = -1 / slope if (np.isfinite(slope) and slope < 0) else np.nan
    return ctau, mean


def betagamma_cdf(output, anchor, channel=CHANNEL):
    """Empirical F_βγ and ⟨βγ⟩ from an (un-truncated) anchor sample's βγ histogram."""
    h = output[anchor]["hists"]["genAs_betagamma"][{"channel": channel}]
    bg = h.axes[-1].centers
    w = h.values()
    cdf = np.cumsum(w) / max(w.sum(), 1)
    return interp1d(bg, cdf, bounds_error=False, fill_value=(0.0, 1.0)), (bg * w).sum() / max(w.sum(), 1)


def acceptance_fit(output, sample, F_bg, mean, channel=CHANNEL):
    """Option C: fit log(dN/dx) = logA - x/cτ + log F_βγ(R_max/x) for (cτ, R_max).

    ε(x) = F_βγ(R_max/x) is the probability a decay with proper length x survives the
    lab-decay-length cap R_max, given the (lifetime-independent) boost distribution F_βγ."""
    centers, dens, effN, _ = proper_density(output, sample, channel)
    good = effN >= 5
    if good.sum() < 4:
        return np.nan, np.nan
    x = centers[good]
    y = np.log(dens[good])
    yerr = 1 / np.sqrt(effN[good])

    def logmodel(x, logA, ctau, logR):
        return logA - x / ctau + np.log(np.clip(F_bg(np.exp(logR) / x), 1e-12, 1))

    try:
        popt, _ = curve_fit(
            logmodel, x, y, p0=[y.max(), mean, np.log(800.0)], sigma=yerr,
            bounds=([-50, 1e-6, np.log(10)], [50, 1e5, np.log(1e6)]), maxfev=40000,
        )
        return popt[1], np.exp(popt[2])
    except Exception:
        return np.nan, np.nan


def group_by_mass_point(output):
    """{(channel, BS, Dp): [samples sorted by cτ]}."""
    g = defaultdict(list)
    for s in output:
        g[mass_point(s)].append(s)
    for k in g:
        g[k].sort(key=ctau_cm)
    return dict(g)


def compute_grid(output, channel=CHANNEL):
    """Per-sample nominal / mean / Option-B / Option-C cτ + R_max, and ⟨βγ⟩ per mass point.

    F_βγ for each mass point is taken from its shortest-cτ (faithful) anchor, since βγ
    depends on the production kinematics (masses), not the lifetime."""
    groups = group_by_mass_point(output)
    rows = {}
    bg_mean = {}
    for key, samples in groups.items():
        F_bg, bgm = betagamma_cdf(output, samples[0], channel)
        bg_mean[key] = bgm
        for s in samples:
            ctB, mean = core_slope_fit(output, s, channel)
            ctC, Rmax = acceptance_fit(output, s, F_bg, mean, channel)
            rows[s] = dict(mass_point=key, nominal=ctau_cm(s), mean=mean,
                           optionB=ctB, optionC=ctC, Rmax=Rmax)
    return rows, groups, bg_mean
