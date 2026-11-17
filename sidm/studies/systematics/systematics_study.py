"""Helpers shared by the notebooks of the systematics study

Notebook-side code only, like corrections_study, which this builds on for loading the
merged outputs of the corrections campaign and for the figure conventions. It adds the
working-point cuts of the limit planes, the signal points the limits use, the
isolation efficiencies read from the per-lepton-jet histograms, the statistical
precision of a pileup variation, and, for the limits notebooks, the reading of the
ABCD inputs and their systematic variations through the limit_plotting tools and the
running of the Combine scripts in sidm/scripts.
"""

import argparse
import collections
import functools
import json
import math
import operator
import os
import shutil
import subprocess
import sys
from pathlib import Path

import fsspec
import hist
import matplotlib.pyplot as plt
import numpy as np

from sidm import BASE_DIR
from sidm.studies.corrections import corrections_study as study
from sidm.tools.utilities import get_eff_hist, load_yaml

# datacard_tools and shape_tools import each other by module name, as the limit
# scripts use them, so their directory goes on the path, as does that of the Combine
# scripts, which the notebooks run
LIMIT_PLOTTING = os.path.join(BASE_DIR, "studies", "limit_plotting")
SCRIPTS = os.path.join(BASE_DIR, "scripts")
for _path in (LIMIT_PLOTTING, SCRIPTS):
    if _path not in sys.path:
        sys.path.insert(0, _path)
# pylint: disable=wrong-import-position,wrong-import-order,import-error
import datacard_tools as dt  # noqa: E402
import shape_tools as st  # noqa: E402
import run_combine_limits as rcl  # noqa: E402
# pylint: enable=wrong-import-position,wrong-import-order,import-error

REGISTRY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "systematics.yaml")

# The lepton-jet isolation cuts of region A in the limit planes
# (sidm/studies/limit_plotting): the muon-type lepton jet below 0.3 and the other
# lepton jet below 0.5, which is the egamma-type one in 2mu2e and the subleading
# muon-type one in 4mu. Each entry is (is_mu, cut) for the per-lepton-jet histograms.
ISO_CUTS = {
    "2Mu2E": {"muon LJ": (True, 0.3), "other LJ": (False, 0.5)},
    "4Mu": {"muon LJ": (True, 0.3), "other LJ": (True, 0.5)},
}

# The campaign channel used for every number here: the category with |dphi| > 2.
# Its isolation cuts act on mu_ljs/egm_ljs, not on the lepton jets the category is
# built from (corrections README, "Found along the way"), so it is inclusive in the
# isolation, as a study of migrations across the isolation cut needs.
CHANNEL = dict(study.SR)


def lifetime_rank(samples):
    """{sample: 1 for the shortest lifetime of its (final state, m_phi, m_zd), 2, ...}"""
    groups = collections.defaultdict(list)
    for s in samples:
        p = study.parse_point(s)
        groups[(p["final_state"], p["m_phi"], p["m_zd"])].append((p["ctau"], s))
    return {s: i + 1 for points in groups.values() for i, (_, s) in enumerate(sorted(points))}


def limit_grid(samples):
    """The signal points the limits use: m_phi of 200 GeV and above (the masses with a
    model cross section), without the shortest lifetime (median L_xy about 0.2 cm),
    which the group no longer considers. 48 points per final state."""
    rank = lifetime_rank(samples)
    return [s for s in samples if study.parse_point(s)["m_phi"] >= 200 and rank[s] >= 2]


def an_hypotheses(samples):
    """The signal hypotheses of the AN's limits (Sec. 8): m_phi of 200 GeV and above, every
    lifetime, the shortest included. 60 points per final state; limit_grid is a subset."""
    return [s for s in samples if study.parse_point(s)["m_phi"] >= 200]


def iso_efficiency(h, is_mu, cut, scale=1.0):
    """Fraction of the lepton jets of one type with isolation below cut

    h is a per-lepton-jet isolation histogram (is_mu x iso on [0, 0.5) in 0.01 bins).
    scale multiplies every isolation, which moves the cut to cut/scale; inside a bin
    the entries are taken as uniform. Values of 0.5 and above are in the overflow, so a
    cut that moves above 0.5 cannot gain anything: only scale >= 1 is exact there.
    """
    values = h[{"is_mu": is_mu}].values(flow=True)  # underflow, 50 bins, overflow
    total = values.sum()
    if total <= 0:
        return np.nan
    edges = h.axes["iso"].edges
    cut = cut / scale
    inside = values[1:-1]
    passed = values[0] + inside[edges[1:] <= cut + 1e-12].sum()
    partial = (edges[:-1] < cut) & (edges[1:] > cut + 1e-12)
    if partial.any():
        i = int(np.argmax(partial))
        passed += inside[i] * (cut - edges[i]) / (edges[i + 1] - edges[i])
    return passed / total


def event_iso_shift(nominal, varied, final_state, scale=1.0):
    """Relative change of the two-lepton-jet isolation efficiency of region A

    Sum over the two lepton jets of the relative change of the per-lepton-jet
    efficiency at its cut: exact to first order in the change. nominal and varied are
    per-lepton-jet isolation histograms; scale is passed to iso_efficiency for the
    varied one.
    """
    total = 0.0
    for is_mu, cut in ISO_CUTS[final_state].values():
        e_nom = iso_efficiency(nominal, is_mu, cut)
        e_var = iso_efficiency(varied, is_mu, cut, scale)
        total += e_var / e_nom - 1 if e_nom > 0 else 0.0
    return total


def ordered_points(samples, final_state):
    """The signal points of one final state, ordered by m_phi, then m_zd, then ctau"""
    points = [s for s in samples if study.parse_point(s)["final_state"] == final_state]
    return sorted(points, key=lambda s: tuple(study.parse_point(s)[k] for k in ("m_phi", "m_zd", "ctau")))


def mass_ticks(ax, points, label=True):
    """Separate and label the m_phi groups on an axis with one entry per point of ordered_points

    label=False leaves out the axis title, for the upper panels of a stack.
    """
    masses = [study.parse_point(s)["m_phi"] for s in points]
    starts = [0] + [i for i in range(1, len(masses)) if masses[i] != masses[i - 1]]
    ends = starts[1:] + [len(masses)]
    for start in starts[1:]:
        ax.axvline(start - 0.5, color=study.GRAY, lw=0.8, ls=":")
    ax.set_xticks([(a + b - 1) / 2 for a, b in zip(starts, ends)])
    ax.set_xticklabels([f"{masses[a]:g}" for a in starts])
    if label:
        ax.set_xlabel(r"$m_\Phi$ (GeV); within each mass, $m_{Z_d}$ and then $c\tau$ increase")


def headroom(ax, fraction=0.35):
    """Extend the y range upwards by a fraction of itself, to keep a legend clear of the data"""
    low, high = ax.get_ylim()
    ax.set_ylim(low, high + fraction * (high - low))


def _ratios(h_nominal, h_varied):
    """Nominal yields, their variances, and the varied-over-nominal ratio per bin"""
    y = h_nominal.values(flow=True)
    ratio = np.divide(h_varied.values(flow=True), y, out=np.ones_like(y), where=y > 0)
    return y, h_nominal.variances(flow=True), ratio


def _weighted_mean(y, v, d):
    """Mean of a per-event quantity d that is constant within each bin, and its uncertainty"""
    total = y.sum()
    mean = np.sum(y * d) / total
    return float(mean), float(np.sqrt(np.sum(v * (d - mean) ** 2)) / total)


def pileup_shift(h_nominal, h_varied):
    """(relative shift of the yield, its statistical uncertainty) from nTrueInt histograms

    The pileup weight is a function of nTrueInt alone, so in each bin b every event
    changes by the same factor r_b = Y_varied(b) / Y_nominal(b). The shift is
    sum_b Y(b) (r_b - 1) / sum_b Y(b) and its variance sum_b V(b) (r_b - 1 - shift)^2 /
    (sum_b Y(b))^2, with V(b) the sum of squared nominal weights in the bin. In the
    campaign's binning three bins hold two integers (corrections README, "Integer
    quantities on bin edges"); there r_b is their average, a negligible approximation.
    """
    y, v, ratio = _ratios(h_nominal, h_varied)
    return _weighted_mean(y, v, ratio - 1)


def pileup_halfwidth(h_nominal, h_up, h_down):
    """((up - down) / 2 of the relative yield change, its statistical uncertainty)

    The up and down shifts come from the same events and are anti-correlated, so the
    uncertainty is computed for the per-event half-difference directly rather than by
    adding the two shifts' uncertainties.
    """
    y, v, r_up = _ratios(h_nominal, h_up)
    _, _, r_down = _ratios(h_nominal, h_down)
    return _weighted_mean(y, v, (r_up - r_down) / 2)


# --- pileup in data: the good primary vertices in the Z -> mumu control region ---------
# Approved 2026-11-16 (the advisor, via the user): data in Z_to_mumu_control_region, the
# good-PV count only (campaign npv_data_v1, every 100th file of each run period). The
# simulation is the corrections_v1 campaign, which has the same channel and histogram.
NPV_CHANNEL = "Z_to_mumu_control_region"
NPV_DATA_DIR = "/store/user/scampbel/sidm_coffea_outputs/npv_data_v1"
NPV_RESULTS_DIR = "/store/user/scampbel/sidm_coffea_outputs/npv_data_v1_results"
# only_pileup: the pileup weight alone; pileupUp/Down are relative to nominal (all corrections)
NPV_CATEGORIES = ("only_pileup", "uncorrected", "nominal", "pileupUp", "pileupDown")
# the background groups of the limits: DYJetsToLL, which holds the muon pairs of DYJetsToMuMu
NPV_BACKGROUNDS = {name: study.BACKGROUNDS[name]
                   for name in (r"DY$\to\ell\ell$", "QCD", r"$t\bar{t}$", "diboson")}


def npv_data(directory=NPV_DATA_DIR):
    """(values, variances, edges) of the good-PV count in data in NPV_CHANNEL, summed over
    the data samples of directory"""
    total = None
    for url in list_coffea(directory):
        for sample, out in dt.read_coffea(url).items():
            if dt.is_simulation(sample, out):
                raise ValueError(f"{url}: {sample} is simulation, not data")
            h = out["hists"]["pv_npvsGood"][{"channel": NPV_CHANNEL}]
            total = h if total is None else total + h
    return total.values(), total.variances(), total.axes[0].edges


def npv_simulation(outputs):
    """{category: {group: (values, variances)}} of the good-PV count of the simulated
    background in NPV_CHANNEL, for each of NPV_CATEGORIES, from corrections_study.load"""
    out = {category: {} for category in NPV_CATEGORIES}
    for group, samples in NPV_BACKGROUNDS.items():
        for category in NPV_CATEGORIES:
            values = variances = 0.0
            for sample in samples:
                h = study.syst_hist(outputs[sample]["hists"]["pv_npvsGood"], category, NPV_CHANNEL)
                h = h[{"channel": NPV_CHANNEL}]
                values, variances = values + h.values(), variances + h.variances()
            out[category][group] = (values, variances)
    return out


def plot_npv(data, simulation, edges, xmax=70):
    """The good-PV count in data against the simulated background, stacked by group, every
    simulated distribution normalised to the data: with the pileup weight (stack), with the
    genWeight alone (dashed), and the change under the pileup variations (band in the ratio).
    Returns (figure, {"data_mean", "mc_mean", "mc_mean_uncorrected", "chi2", "ndf"})."""
    import plot_tools as pt   # pylint: disable=import-outside-toplevel,import-error
    y, v = data[0], data[1]
    totals = {c: sum(vals for vals, _ in simulation[c].values()) for c in NPV_CATEGORIES}
    scale = {c: y.sum() / totals[c].sum() for c in NPV_CATEGORIES}
    weighted = totals["only_pileup"] * scale["only_pileup"]
    weighted_var = (sum(var for _, var in simulation["only_pileup"].values())
                    * scale["only_pileup"] ** 2)
    nominal = totals["nominal"] * scale["nominal"]
    centers = 0.5 * (edges[1:] + edges[:-1])

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    bottom = np.zeros_like(y)
    colors = (study.BLUE, study.ORANGE, study.AQUA, study.GRAY)
    for (group, (vals, _)), color in zip(simulation["only_pileup"].items(), colors):
        ax.stairs(bottom + vals * scale["only_pileup"], edges, baseline=bottom, fill=True,
                  color=color, label=group)
        bottom = bottom + vals * scale["only_pileup"]
    ax.stairs(totals["uncorrected"] * scale["uncorrected"], edges, color=study.BLACK, ls="--",
              lw=1.5, label="simulation, genWeight only")
    ax.errorbar(centers, y, yerr=np.sqrt(v), fmt="o", color="k", ms=4, label="data")
    ax.set_ylabel("events (simulation normalised to data)")
    ax.set_xlim(edges[0], xmax)
    ax.legend(fontsize=13)

    ok = (weighted > 0) & (nominal > 0)
    ratio = np.divide(y, weighted, out=np.full_like(y, np.nan), where=ok)
    rax.errorbar(centers, ratio, yerr=np.divide(np.sqrt(v), weighted, out=np.zeros_like(y), where=ok),
                 fmt="o", color="k", ms=4)
    up, down = (np.divide(totals[c] * scale[c], nominal, out=np.ones_like(y), where=ok)
                for c in ("pileupUp", "pileupDown"))
    rax.fill_between(edges[:-1], np.minimum(up, down), np.maximum(up, down), step="post",
                     color=study.ORANGE, alpha=0.4, lw=0, label="pileup weight variations")
    rax.stairs(np.divide(totals["uncorrected"] * scale["uncorrected"], weighted,
                         out=np.ones_like(y), where=ok), edges, color=study.BLACK, ls="--", lw=1.5)
    rax.axhline(1.0, color="k", lw=0.8)
    rax.set_ylim(0.5, 1.5)
    rax.set_xlabel("number of good primary vertices")
    rax.set_ylabel("data / simulation")
    rax.legend(fontsize=12, loc="upper left")
    pt.cms_label(ax, lumi=None, data=True, text="Work in progress")

    # shape agreement where both have entries; the simulation is normalised to the data,
    # so one degree of freedom goes into that
    keep = ok & (y > 0)
    chi2 = float(np.sum((y[keep] - weighted[keep]) ** 2 / (v[keep] + weighted_var[keep])))

    def mean(values):
        """mean number of good vertices of a distribution"""
        return float(np.sum(values * centers) / np.sum(values))
    summary = {"data_mean": mean(y), "mc_mean": mean(weighted),
               "mc_mean_uncorrected": mean(totals["uncorrected"]),
               "mc_mean_up": mean(totals["pileupUp"]), "mc_mean_down": mean(totals["pileupDown"]),
               "chi2": chi2, "ndf": int(keep.sum()) - 1}
    return fig, summary


def npv_ratio(data, simulation):
    """Data over the pileup-weighted simulation of the good-PV count in the Z -> mumu control
    region, both normalised to unit area, per bin (one where the simulation is empty)"""
    y = data[0] / data[0].sum()
    total = sum(values for values, _ in simulation["only_pileup"].values())
    mc = total / total.sum()
    return np.divide(y, mc, out=np.ones_like(y), where=mc > 0)


def npv_reweighted_shift(sample_out, channel, ratio):
    """Relative change of a sample's yield in channel when its good-PV distribution (with the
    pileup weight alone) is weighted by ratio (npv_ratio): how much the signal yield would move
    if the simulation were made to match the vertex count of the data"""
    h = study.syst_hist(sample_out["hists"]["pv_npvsGood"], "only_pileup", channel)
    values = h[{"channel": channel}].values()
    return float(np.sum(values * ratio) / values.sum() - 1) if values.sum() > 0 else float("nan")


# --- the ABCD inputs of the limits ---------------------------------------------------

# Region A of the working points (iso x |dphi| plane, m_LJLJ >= 150 GeV). 2mu2e: muon-type
# LJ 0 <= iso < 0.3, |dphi| >= 2.6 against a control side below 2.4, egamma-type LJ
# iso < 0.5; the default displacement cuts, at most two pixel hits on the muon-type LJ and
# at least one missing hit on the egamma-type LJ. 4mu: leading muon-type LJ 0 <= iso < 0.3,
# |dphi| >= 2.0, subleading iso < 0.5; at most three pixel hits on both LJs (the group's
# loosened cut; pixel_hit_max is the upper edge of the 3-hit bin of the abcd6d_* axes).
PLANE = "iso_dphi"
WORKING_POINTS = {
    "SR_2mu2e": st.Selection("2mu2e working point", iso_cut=0.3, dphi_cut=2.6, dphi_cr=2.4,
                             other_iso_cut=0.5, iso_floor=0.0),
    "SR_4mu": st.Selection("4mu working point", iso_cut=0.3, dphi_cut=2.0,
                           other_iso_cut=0.5, iso_floor=0.0, pixel_hit_max=3.5),
}
FINAL_STATE_CHANNEL = {"2Mu2E": "SR_2mu2e", "4Mu": "SR_4mu"}
# DYJetsToMuMu holds the muon pairs of DYJetsToLL, so the two must never be added
# (corrections_study.BACKGROUNDS); the limits use DYJetsToLL, as the earlier
# productions did (sixd_inclusive_v1 has no DYJetsToMuMu)
EXCLUDED_BACKGROUNDS = ("DYJetsToMuMu_M10to50", "DYJetsToMuMu_M50")


def load_registry(path=REGISTRY):
    """{nuisance: entry} of systematics.yaml, each entry checked to be of exactly one kind"""
    registry = load_yaml(path)
    for name, entry in registry.items():
        kinds = [k for k in ("flat", "syst", "variant") if k in entry]
        if len(kinds) != 1:
            raise ValueError(f"{path}: {name} needs exactly one of flat, syst, variant; has {kinds}")
    return registry


def region_totals(regions):
    """{region index: summed yield} of {region index: Shape}"""
    return {r: float(shape.sumw.sum()) for r, shape in regions.items()}


def variations(out, sample, channel, entry, nominal):
    """{"up": {region index: Shape}, "down": ...} of one systematics.yaml entry for one
    signal sample in its channel, nominal being its {region index: Shape}

    A flat entry scales the nominal shapes, a syst entry takes the Up and Down weight
    categories of the syst axis, and a variant has an up variation only: the histogram
    filled with the varied input.
    """
    if "flat" in entry:
        factor = 1 + entry["flat"]
        return {"up": {r: s.scaled(factor) for r, s in nominal.items()},
                "down": {r: s.scaled(1 / factor) for r, s in nominal.items()}}

    def regions(**varied):
        return st.region_shapes_6d(out, channel, PLANE, sample,
                                   selection=WORKING_POINTS[channel], **varied)

    if "syst" in entry:
        return {"up": regions(syst=f"{entry['syst']}Up"),
                "down": regions(syst=f"{entry['syst']}Down")}
    return {"up": regions(variant=entry["variant"])}


def signal_kappas(out, sample, channel, nominal, registry):
    """{nuisance: {region name: kappa}} of one signal point, a kappa being a number or [down, up]

    nominal is the point's {region index: Shape} in its channel; the varied yields are
    sliced the same way (variations). A flat entry gives its value as it is.
    """
    nom = region_totals(nominal)
    kappas = {}
    for name, entry in registry.items():
        if "flat" in entry:
            values = {r: 1 + entry["flat"] for r in nom}
        else:
            ratio = {}
            for direction, shapes in variations(out, sample, channel, entry, nominal).items():
                totals = region_totals(shapes)
                ratio[direction] = {r: totals[r] / nom[r] if nom[r] > 0 else 1.0 for r in nom}
            values = ({r: [ratio["down"][r], ratio["up"][r]] for r in nom} if "down" in ratio
                      else ratio["up"])
        kappas[name] = {dt.ABCD_REGIONS[r]: v for r, v in values.items()}
    return kappas


def kappa_size(kappa):
    """Size of a log-normal kappa as a fraction: |kappa - 1|, or half the difference of
    [down, up]"""
    return abs(kappa[1] - kappa[0]) / 2 if isinstance(kappa, list) else abs(kappa - 1)


def signal_mc_stat(payload, point, channel, region="0"):
    """Relative statistical uncertainty of a signal point's yield in one region (its index,
    as a string), from a shapes payload: the mcstat row abcd_workspace.py writes"""
    shape = payload["signal"][point][channel][region]
    total = sum(shape["sumw"])
    return math.sqrt(sum(shape["sumw2"])) / total if total > 0 else float("nan")


def freeze_variants(registry):
    """{label: what --freezeParameters fixes} of the limit breakdown: each systematic of the
    registry left out in turn, all of them, and every constrained nuisance, the statistical
    uncertainty of the simulated signal included, which leaves the statistics of the
    (Asimov) data alone"""
    variants = {f"without {name}": name for name in registry}
    variants["without any systematic"] = ",".join(registry)
    variants["statistics only"] = "allConstrainedNuisances"
    return variants


def abcd_inputs(directory, registry, weights=("nominal",)):
    """Background shapes, signal shapes and signal kappas from one campaign's merged outputs

    One pass over the .coffea files of directory, which holds signal and background
    (told apart by name, as abcd_datacards.py does for such a campaign), leaving out
    EXCLUDED_BACKGROUNDS. The shapes are made with each event weight in weights, a
    category of the syst axis ("nominal" is with the corrections, "uncorrected" the
    genWeight alone), and with "nominal" in any case, since the kappas are relative to it.
    Returns {weights: (background, signal)}, background being
    {sample: {channel: {region: Shape}}} and signal the same for the signal in the
    channel of its final state, and {point: {nuisance: {region name: kappa}}}.
    """
    categories = ["nominal"] + [w for w in weights if w != "nominal"]
    shapes = {w: ({}, {}) for w in categories}
    kappas = {}
    files = [url for url in list_coffea(directory)
             if url.rsplit("/", 1)[-1][:-len(".coffea")] not in EXCLUDED_BACKGROUNDS]
    if not files:
        raise FileNotFoundError(f"no .coffea files in {directory}")
    for url in files:
        for sample, out in dt.read_coffea(url).items():
            info = dt.parse_signal_name(sample)
            if info is None:
                for w in categories:
                    shapes[w][0][sample] = sample_shapes(out, sample, list(WORKING_POINTS), w)
                continue
            channel = FINAL_STATE_CHANNEL[info["final_state"]]
            per_weight = {w: sample_shapes(out, sample, [channel], w) for w in categories}
            if per_weight["nominal"]:
                for w in categories:
                    shapes[w][1][sample] = per_weight[w]
                kappas[sample] = signal_kappas(out, sample, channel,
                                               per_weight["nominal"][channel], registry)
    return shapes, kappas


def sample_shapes(out, sample, channels, syst="nominal"):
    """{channel: {region: Shape}} of one sample at the working points with the event weight
    syst (a category of the syst axis), channels it lacks left out"""
    shapes = {ch: st.region_shapes_6d(out, ch, PLANE, sample, selection=WORKING_POINTS[ch],
                                      syst=syst)
              for ch in channels}
    return {ch: regions for ch, regions in shapes.items() if regions}


# --- collision data in the control regions ---------------------------------------------
# Approved 2026-11-16 (the advisor, via the user): regions B, C and D of the signal plane in
# data, their m_LJLJ distributions included, for blinded data-driven expected limits. Region
# A of data is never sliced: shape_tools.region_shapes_6d leaves it out of anything not
# positively identified as simulation, and data_control_regions stops if it gets it anyway.
DATA_APPROVAL = ("signal-plane regions B, C and D in data, m_LJLJ included (approved "
                 "2026-11-16); region A never read")
# Where 03_limits.ipynb writes its data-driven limits (next to the campaign's merged
# directory), and the estimate of the region-A background that 04 and 05 check: the
# integrated one (one m_LJLJ bin), the likelihood of the AN's limits (Sec. 8).
LIMITS_SUFFIX = "_limits_data"
LIMITS_METHOD = "integrated"
SR_BLIND = "blinded: never read"   # the mark abcd_workspace.py --sr-data-blinded checks


def data_control_regions(directory):
    """{channel: {region: Shape}} of regions B, C and D in data at the working points,
    summed over the data samples of directory (one campaign's merged outputs)"""
    files = list_coffea(directory)
    if not files:
        raise FileNotFoundError(f"no .coffea files in {directory}")
    total = {}
    for url in files:
        for sample, out in dt.read_coffea(url).items():
            if dt.is_simulation(sample, out):
                raise ValueError(f"{url}: {sample} is simulation, not data")
            for channel, selection in WORKING_POINTS.items():
                regions = st.region_shapes_6d(out, channel, PLANE, sample, selection=selection)
                if not regions:
                    continue
                if dt.SR_ABCD_REGION in regions:
                    raise RuntimeError(f"{sample}: region A of data was sliced: BLINDING")
                slot = total.setdefault(channel, {})
                for region, shape in regions.items():
                    slot[region] = slot[region] + shape if region in slot else shape
    return total


def write_payloads(outdir, background, signal, weights="nominal", data=None):
    """shapes_<channel>.json for abcd_workspace.py, one per channel with its working point,
    recording the event weight (syst category) the shapes were made with. With data
    (data_control_regions), each payload also carries regions B, C and D of data, region A
    marked as never read, for abcd_workspace.py --sr-data-blinded."""
    grouped = st.group_shapes(background)
    paths = {}
    for channel, selection in WORKING_POINTS.items():
        observable = st.Observable(name=f"ljlj_mass_{PLANE}", hists={}, axis="ljlj_mass",
                                   label=r"$m_{LJLJ}$ [GeV]",
                                   note=f"6D slice, plane {PLANE}, {selection.describe()}")
        points = {s: per_channel for s, per_channel in signal.items() if channel in per_channel}
        paths[channel] = st.export_shapes(
            Path(outdir) / f"shapes_{channel}.json", grouped, points, observable,
            channels=[channel],
            meta={"plane": PLANE, "selection": selection.describe(), "weights": weights,
                  "normalisation": "signal normalised to 1 fb at 59.83 /fb (production "
                                   "cross section divided out), so r = sigma / 1 fb"})
        if data is not None:
            payload = json.loads(paths[channel].read_text())
            payload["mc_only"] = False
            payload["data"] = {
                "A": SR_BLIND, "approved": DATA_APPROVAL,
                "channels": {channel: {name: data[channel][region].to_dict()
                                       for region, name in ((1, "B"), (2, "C"), (3, "D"))}}}
            paths[channel].write_text(json.dumps(payload, indent=1))
    return paths


# --- EOS through xrootd ----------------------------------------------------------------
# The limits notebook reads its inputs from EOS and writes its outputs there through
# xrootd (fsspec-xrootd), so it runs the same where EOS is mounted (the LPC) and where it
# is not (CoffeaCasa).

EOS_REDIRECTOR = "root://cmseos.fnal.gov"


def eos_url(path):
    """The xrootd URL of an EOS path (/eos/uscms/store/... or /store/...); URLs and other
    local paths are returned as they are"""
    path = str(path)
    if path.startswith("/eos/uscms/store/"):
        path = path[len("/eos/uscms"):]
    return f"{EOS_REDIRECTOR}/{path}" if path.startswith("/store/") else path


def list_coffea(directory):
    """URLs of the .coffea files of a directory on EOS (or local), in name order"""
    fs, root = fsspec.core.url_to_fs(eos_url(directory))
    if not fs.exists(root):
        return []
    return sorted(fs.unstrip_protocol(p) for p in fs.glob(f"{root.rstrip('/')}/*.coffea"))


def upload(local_dir, directory):
    """Copy every file under local_dir to directory on EOS (or local), keeping the tree

    Written through fs.open, since fsspec-xrootd does not implement put_file, and made
    readable by everyone afterwards, as xrdcp leaves them (an xrootd open makes them private).
    """
    fs, root = fsspec.core.url_to_fs(eos_url(directory))
    for path in sorted(Path(local_dir).rglob("*")):
        if path.is_file():
            target = f"{root.rstrip('/')}/{path.relative_to(local_dir)}"
            fs.makedirs(os.path.dirname(target), exist_ok=True)
            with open(path, "rb") as source, fs.open(target, "wb") as copy:
                shutil.copyfileobj(source, copy)
            if hasattr(fs, "chmod"):
                fs.chmod(target, 0o644)
    # The directory listings cached while writing are stale now, and with a listing of a
    # parent in its cache fsspec-xrootd fails to open the files below it ("_ls_from_cache()
    # failed to function"), so reading the outputs back would fail.
    fs.invalidate_cache()


def download(directory, local_dir):
    """Copy the files of directory on EOS (or local) into local_dir, for Combine, which
    reads local files only"""
    fs, root = fsspec.core.url_to_fs(eos_url(directory))
    os.makedirs(local_dir, exist_ok=True)
    for entry in fs.ls(root, detail=True):
        if entry["type"] == "file":
            target = os.path.join(local_dir, os.path.basename(entry["name"]))
            with fs.open(entry["name"], "rb") as source, open(target, "wb") as copy:
                shutil.copyfileobj(source, copy)
    fs.invalidate_cache()   # as in upload


def limits_outputs(limits_dir, registry):
    """The signal systematics 03_limits.ipynb wrote to limits_dir, {"sources", "kappas",
    "registry"}, checked to come from the current registry, so that a notebook built on
    the outputs of 03_limits stops when 03_limits needs a rerun"""
    try:
        with fsspec.open(f"{eos_url(limits_dir)}/inputs/with/systematics.json") as f:
            systematics = json.load(f)
    except OSError as err:
        raise FileNotFoundError(f"no outputs of 03_limits.ipynb in {limits_dir}: "
                                "run it first") from err
    if systematics.get("registry") != registry:
        raise RuntimeError(f"the outputs in {limits_dir} do not come from the current "
                           "systematics.yaml: rerun 03_limits.ipynb first")
    return systematics


# --- Combine, from the notebooks -------------------------------------------------------

def find_combine():
    """Where Combine is, as the scripts in sidm/scripts take it: an argparse.Namespace
    with no_cmsenv and combine_cmssw, or None. Combine is taken from the PATH if it is
    there (e.g. a conda installation), else from the release run_combine_limits.py uses
    at the LPC, which COMBINE_CMSSW_BASE overrides."""
    release = os.environ.get("COMBINE_CMSSW_BASE", rcl.DEFAULT_COMBINE_CMSSW)
    if shutil.which("combine"):
        return argparse.Namespace(no_cmsenv=True, combine_cmssw=None)
    if os.path.isdir(release):
        return argparse.Namespace(no_cmsenv=False, combine_cmssw=release)
    return None


def run_logged(command, log, check=True):
    """Run a command with its output in the file log, and return its exit code; with
    check, a failure raises with the end of the log"""
    with open(log, "w", encoding="utf-8") as f:
        proc = subprocess.run([str(c) for c in command], stdout=f, stderr=subprocess.STDOUT,
                              check=False)
    if proc.returncode and check:
        raise RuntimeError(f"failed, see {log}:\n" + Path(log).read_text(encoding="utf-8")[-2000:])
    return proc.returncode


def run_script(script, *args, combine, log, check=True):
    """Run a script of sidm/scripts with this Python, told where Combine is (find_combine),
    with its output in the file log (run_logged)"""
    where = ["--no-cmsenv"] if combine.no_cmsenv else ["--combine-cmssw", combine.combine_cmssw]
    return run_logged([sys.executable, os.path.join(SCRIPTS, script), *args, *where], log, check)


# --- The lepton-jet reconstruction efficiency (lj_reco_efficiency.md, step 1) --------------

# The merged outputs of condor/campaign_configs/lj_reco_eff_v1.yaml, written out as
# corrections_study.CAMPAIGN_DIR is, so the notebook reads them wherever it runs
LJ_RECO_DIR = "/store/user/scampbel/sidm_coffea_outputs/lj_reco_eff_v1"

# Per lepton-jet type, the channels of the efficiency chain of sidm/studies/ljRecoEffi_studies
# and its per-dark-photon histograms: L_xy, pT, and the map of L_xy against the opening angle
# of the two leptons. The denominator holds every generator dark photon decaying to that
# lepton, in events that pass the base selection; the numerator those within Delta R < 0.4 of
# a lepton jet, in events with one matched to a selected lepton jet of that type.
LJ_EFFICIENCY = {
    "muon-type": {"num": "baseNoLj_A_mumu_matched_lj", "den": "baseNoLj_A_mumu",
                  "lxy": "genAs_toMu_lxy", "pt": "genAs_toMu_pt",
                  "map": "genA_toMu_lxy_vs_daughters_dR"},
    "egamma-type": {"num": "baseNoLj_A_ee_matched_lj", "den": "baseNoLj_A_ee",
                    "lxy": "genAs_toE_lxy", "pt": "genAs_toE_pt",
                    "map": "genA_toE_lxy_vs_daughters_dR"},
}

# Where the dark photon decays: inside the pixel layers, in the strip tracker, beyond the
# tracker; the bounds are bin edges of the L_xy histograms (5 cm bins)
LXY_REGIMES = {"Lxy < 10 cm": (0.0, 10.0), "10 cm <= Lxy < 100 cm": (10.0, 100.0),
               "Lxy >= 100 cm": (100.0, np.inf)}


def lj_efficiency_inputs(directory=LJ_RECO_DIR):
    """{sample: {lepton-jet type: {variable: (numerator, denominator)}}} from the merged
    outputs of the lj_reco_eff_v1 campaign, variable being "lxy", "pt" or "map"; a type
    without dark photons in a sample (egamma-type in 4mu) is left out"""
    inputs = {}
    for url in list_coffea(directory):
        for sample, out in dt.read_coffea(url).items():
            hists = out["hists"]
            per_type = {}
            for lj_type, spec in LJ_EFFICIENCY.items():
                pairs = {v: (hists[spec[v]][{"channel": spec["num"]}],
                             hists[spec[v]][{"channel": spec["den"]}]) for v in ("lxy", "pt", "map")}
                if pairs["lxy"][1].sum(flow=True).value > 0:
                    per_type[lj_type] = pairs
            inputs[sample] = per_type
    return inputs


def counts(hist_):
    """A histogram of one sample as counts of simulated events. Every event of a signal sample
    carries the same weight, so (sum of weights)^2 / (sum of squared weights) per bin is the
    number of events in it."""
    out = hist_.copy()
    view = out.view(flow=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        n = np.where(view.variance > 0, view.value ** 2 / view.variance, 0.0)
    view.value, view.variance = n, n
    return out


def summed(inputs, samples, lj_type, variable):
    """(numerator, denominator) of one variable, summed over the samples that have it as
    counts of simulated dark photons, so each counts once whatever the cross section of its
    sample. Weighted sums would let the samples of the largest cross sections dominate, and
    their effective counts could put more events in a numerator bin than in its denominator."""
    pairs = [inputs[s][lj_type][variable] for s in samples if lj_type in inputs[s]]
    return (functools.reduce(operator.add, (counts(num) for num, _ in pairs)),
            functools.reduce(operator.add, (counts(den) for _, den in pairs)))


def efficiency(num, den):
    """utilities.get_eff_hist (efficiency histogram, Clopper-Pearson errors), quiet about
    empty bins"""
    with np.errstate(divide="ignore", invalid="ignore"):
        return get_eff_hist(num, den)


def regime_fractions(num, den):
    """The efficiency of one sample and, per LXY_REGIMES entry, the fraction of its matched
    dark photons decaying there, from its L_xy histograms (flow bins included)"""
    edges = num.axes[0].edges
    lower = np.concatenate([[-np.inf], edges[:-1], [edges[-1]]])   # of each flow bin
    matched = num.values(flow=True)
    total = matched.sum()
    fractions = {name: float(matched[(lower >= lo) & (lower < hi)].sum() / total)
                 if total > 0 else np.nan for name, (lo, hi) in LXY_REGIMES.items()}
    return float(total / den.values(flow=True).sum()), fractions


def efficiency_points(inputs, samples, lj_type, variable):
    """Bin centres, efficiencies and their errors of one lepton-jet type against "lxy", "pt" or
    "dR" (the map projected over L_xy, in logarithmic bins), summed over samples"""
    if variable == "dR":
        num, den = summed(inputs, samples, lj_type, "map")
        num, den = num.project(num.axes[1].name), den.project(den.axes[1].name)
        edges = num.axes[0].edges
        centers = np.sqrt(edges[:-1] * edges[1:])
    else:
        num, den = summed(inputs, samples, lj_type, variable)
        centers = num.axes[0].centers
    eff, errors = efficiency(num, den)
    return centers, eff.values(), errors


def plot_lj_efficiency(inputs, groups, lj_type, title):
    """The efficiency of one lepton-jet type against L_xy, the dark photon pT and the opening
    angle of its leptons, one curve per group of samples ({label: [samples]}); the opening
    angle comes from the map, projected over L_xy"""
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.2), sharey=True)
    for (label, samples), colour in zip(groups.items(), (study.BLUE, study.ORANGE, study.AQUA,
                                                          study.GRAY)):
        if not any(lj_type in inputs[s] for s in samples):
            continue
        for ax, variable in zip(axes, ("lxy", "pt", "dR")):
            x, y, errors = efficiency_points(inputs, samples, lj_type, variable)
            ax.errorbar(x, y, yerr=errors, fmt="o", ms=3.5, color=colour, label=label)
    for ax, xlabel in zip(axes, (r"dark photon $L_{xy}$ (cm)", r"dark photon $p_T$ (GeV)",
                                 r"$\Delta R(\ell, \ell)$ of the dark photon")):
        ax.set_xlabel(xlabel)
    axes[2].set_xscale("log")
    axes[0].set_ylabel(f"{lj_type} lepton-jet efficiency")
    axes[0].set_ylim(0, 1.15)
    axes[0].legend(fontsize=12, title=title, title_fontsize=12, loc="upper right")
    study.cms_label(axes)
    return fig


def plot_lj_efficiency_map(inputs, samples, lj_type, ax, min_events=2.0):
    """The efficiency of one lepton-jet type in the plane of L_xy and the opening angle, summed
    over samples (summed); bins with fewer than min_events simulated dark photons are blank"""
    num, den = summed(inputs, samples, lj_type, "map")
    eff, _ = efficiency(num, den)
    values = np.where(den.values() >= min_events, eff.values(), np.nan)
    mesh = ax.pcolormesh(num.axes[0].edges, num.axes[1].edges, values.T, vmin=0, vmax=1,
                         cmap="viridis")
    ax.set_yscale("log")
    ax.set_xlabel(r"dark photon $L_{xy}$ (cm)")
    ax.set_ylabel(r"$\Delta R(\ell, \ell)$")
    return mesh


def plot_regimes(inputs, points, lj_type, title):
    """Per signal point (ordered_points), the efficiency of one lepton-jet type, and where its
    matched dark photons decay (LXY_REGIMES)"""
    results = [regime_fractions(*inputs[s][lj_type]["lxy"]) for s in points]
    x = np.arange(len(points))
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(15, 7.5), sharex=True,
                                      gridspec_kw={"height_ratios": [1, 1.3], "hspace": 0.06})
    top.plot(x, [eff for eff, _ in results], "o", ms=4, color=study.BLACK)
    top.set_ylim(0, 1.05)
    top.set_ylabel("efficiency")
    top.set_title(title, fontsize=14)
    stacked = np.zeros(len(points))
    for name, colour in zip(LXY_REGIMES, (study.BLUE, study.ORANGE, study.AQUA)):
        values = np.array([fractions[name] for _, fractions in results])
        bottom.bar(x, values, bottom=stacked, width=0.85, color=colour, label=name)
        stacked += values
    bottom.set_ylim(0, 1.25)
    bottom.set_ylabel("fraction of the matched\ndark photons")
    bottom.legend(ncol=3, fontsize=12, loc="upper left")
    mass_ticks(top, points, label=False)
    mass_ticks(bottom, points)
    study.cms_label(top)
    return fig


# ------------------------------------------------------------------------------------------
# Lepton-jet control samples (lj_reco_efficiency.md, step 2)

# The merged simulation outputs of condor/campaign_configs/lj_reco_cr_v2.yaml
LJ_CR_DIR = "/store/user/scampbel/sidm_coffea_outputs/lj_reco_cr_v2"
# Per control sample: its channel and probe histogram (sidm/configs/selections.yaml,
# sidm/definitions/hists.py) in the version a data campaign reads, whose probes have no other
# lepton jet of the signal selection's kinds (objects._alone), so no event of the signal region
# enters; the version without that condition (*_all), for the simulation only; the lepton-jet
# type it tests; and its mass axis with the peak and the sidebands that subtract the background
# under the peak (GeV). For the Z, the photons that do not come from the Z decay make a
# background that rises steeply with m(mumu gamma), which a straight line between distant
# sidebands overestimates under the peak; a narrow peak with adjacent sidebands keeps the
# subtraction close to linear (first run of 07: Drell-Yan efficiencies above 1 with 80-100 GeV
# and sidebands of 20 GeV, 0.86 to 0.99 with these windows, as a quadratic fit gives)
LJ_CR = {
    "Z -> mumu gamma": {"channel": "lj_cr_Z_to_mumugamma_alone", "hist": "lj_cr_zmmg_probes_alone",
                        "channel_all": "lj_cr_Z_to_mumugamma", "hist_all": "lj_cr_zmmg_probes",
                        "label": r"$Z \to \mu\mu\gamma$",
                        "lj_type": "egamma-type", "mass": "m_mmg", "peak": (85.0, 95.0),
                        "sidebands": ((75.0, 85.0), (95.0, 105.0))},
    "J/psi -> mumu": {"channel": "lj_cr_JPsi_to_mumu_alone", "hist": "lj_cr_jpsi_probes_alone",
                      "channel_all": "lj_cr_JPsi_to_mumu", "hist_all": "lj_cr_jpsi_probes",
                      "label": r"$J/\psi \to \mu\mu$",
                      "lj_type": "muon-type", "mass": "m_mumu", "peak": (2.9, 3.3),
                      "sidebands": ((2.6, 2.9), (3.3, 3.6))},
}

# The group's criteria for reading a region of data, agreed in the review of the first VR-safety
# deck (2026-11-04) and applied to every working point since (for example the 4mu working point
# with |dphi| >= 2.4): judged in simulation, every quantity at its worst over the 96 points of
# the limit grid (4mu and 2mu2e, L_xy >= 3 cm) at their model cross sections. Safe if the region
# expects under one signal event, or if S/(S+B) < 5% and the Asimov significance Z < 1; not safe
# if S/(S+B) > 10% or Z > 2; caution in between. A shape use is judged with the binned Z. The
# same thresholds are in the local slides/vr_safety.py of sidm/studies/limit_plotting.
SAFE_EVENTS = 1.0
SAFE_CONTAM, UNSAFE_CONTAM = 0.05, 0.10
SAFE_Z, UNSAFE_Z = 1.0, 2.0
# One rule beyond the decks: the bins of a region that expect under one background event are
# judged together as one more region with the same criteria. The decks report the signal in
# those bins but leave it out of the verdict; a fine probe histogram can hold a signal-dominated
# tail there while the region as a whole passes (the first J/psi test: about 5 signal events in
# bins without simulated background). Used with the user's agreement (2026-11-16) until the
# advisor rules on it; set to False to judge exactly as the decks do.
STRICT_SPARSE_BINS = True


def asimov(s, b):
    """Asimov discovery significance of s over b, per element, with no background uncertainty:
    leaving it out makes Z larger, the cautious direction for a safety call. b > 0."""
    s, b = np.asarray(s, dtype=float), np.asarray(b, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        z2 = 2.0 * ((s + b) * np.log1p(s / b) - s)
    return np.sqrt(np.clip(np.nan_to_num(z2), 0.0, None))


def safety_verdict(s, contam, z):
    """("safe" | "caution" | "not safe", branch) of the group's criteria; branch is "events"
    when the S < 1 rule decided it, "rel" otherwise"""
    if s < SAFE_EVENTS:
        return "safe", ("rel" if contam < SAFE_CONTAM and z < SAFE_Z else "events")
    if contam < SAFE_CONTAM and z < SAFE_Z:
        return "safe", "rel"
    if contam > UNSAFE_CONTAM or z > UNSAFE_Z:
        return "not safe", "rel"
    return "caution", "rel"


def region_metrics(s, b, b_w2):
    """Contamination and visibility of the signal s over the background b, both flat arrays of
    bins, b_w2 the background's sums of squared weights. z_binned combines the bins expecting at
    least one background event, z_binned_all every bin with background; s_sparse is the signal
    in the bins expecting under one background event."""
    total_s, total_b = float(s.sum()), float(b.sum())
    meas, pos = b >= 1.0, b > 0
    return {"S": total_s, "B": total_b,
            "neff": total_b ** 2 / float(b_w2.sum()) if b_w2.sum() > 0 else 0.0,
            "contam": total_s / (total_s + total_b) if total_s + total_b > 0 else 0.0,
            "z": float(asimov(total_s, total_b)) if total_b > 0 else (math.inf if total_s > 0 else 0.0),
            "z_binned": float(np.sqrt(np.sum(asimov(s[meas], b[meas]) ** 2))),
            "z_binned_all": float(np.sqrt(np.sum(asimov(s[pos], b[pos]) ** 2))),
            "s_sparse": float(s[~meas].sum())}


def lj_cr_inputs(directory=LJ_CR_DIR, alone=True):
    """{"signal": {point: {control sample: hist}}, "background": {sample: {...}}} of the probe
    histograms (LJ_CR) from the merged simulation outputs, the signal at the 1 fb reference
    (datacard_tools.read_coffea): the versions a data campaign reads, or with alone=False those
    without the condition on other lepton jets. Stops at a data sample: this reads simulation."""
    channel, name = ("channel", "hist") if alone else ("channel_all", "hist_all")
    inputs = {"signal": {}, "background": {}}
    for url in list_coffea(directory):
        for sample, out in dt.read_coffea(url).items():
            if not dt.is_simulation(sample, out):
                raise ValueError(f"{url}: {sample} is data; the safety check reads simulation only")
            kind = "signal" if dt.parse_signal_name(sample) else "background"
            inputs[kind][sample] = {cr: out["hists"][spec[name]][{"channel": spec[channel]}]
                                    for cr, spec in LJ_CR.items()}
    return inputs


def _slices(axis):
    """(label, index) of every bin of a histogram axis in the flow=True layout, flow bins
    included, so that nothing a data histogram holds escapes the check"""
    flow = axis.traits.underflow, axis.traits.overflow
    edges = axis.edges
    out = [(f"< {edges[0]:g}", 0)] if flow[0] else []
    out += [(f"{lo:g} to {hi:g}", i + flow[0]) for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:]))]
    if flow[1]:
        out.append((f">= {edges[-1]:g}", len(edges) - 1 + flow[0]))
    return out


def _heaviest_event(background, take):
    """(weight, sample, flat bin) of the heaviest simulated background event in a region,
    estimated per sample and bin as sumw2/sumw: exact for a bin of one event or of events of one
    weight, a lower bound otherwise"""
    best = (0.0, None, None)
    for sample, h in background.items():
        w, w2 = take(h.values(flow=True)), take(h.variances(flow=True))
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(w > 0, w2 / w, 0.0)
        i = int(np.argmax(ratio))
        if ratio[i] > best[0]:
            best = (float(ratio[i]), sample, i)
    return best


def _judge(signal, b, b_w2):
    """The metrics of each point's signal (signal: {point: flat bins}) over the background b, the
    worst point, and the verdicts of the group's criteria (counting, shape, sparse bins)"""
    metrics = {p: region_metrics(s, b, b_w2) for p, s in signal.items()}
    worst = max(metrics, key=lambda p: metrics[p]["S"])
    contam = max(m["contam"] for m in metrics.values())
    # the bins expecting under one background event, judged together as one more region
    sparse = b < 1.0
    s_sp = max(float(s[sparse].sum()) for s in signal.values())
    b_sp = float(b[sparse].sum())
    z_sp = float(asimov(s_sp, b_sp)) if b_sp > 0 else (math.inf if s_sp > 0 else 0.0)
    verdicts = {
        "counting": safety_verdict(metrics[worst]["S"], contam,
                                   max(m["z"] for m in metrics.values())),
        "shape": safety_verdict(metrics[worst]["S"], contam,
                                max(m["z_binned_all"] for m in metrics.values())),
        "sparse bins": safety_verdict(s_sp, s_sp / (s_sp + b_sp) if s_sp + b_sp > 0 else 0.0, z_sp),
        "s_sparse": s_sp, "b_sparse": b_sp}
    return metrics, worst, verdicts


def _region_row(cr, region, take, ctx):
    """One row of lj_cr_safety: the metrics of a region (take maps a flow=True array of the probe
    histogram to the region's flat bins) and its verdicts, worst case over the points; ctx holds
    the background histograms, the signal histograms and the model cross sections"""
    total = sum(ctx["background"].values())
    b, b_w2 = take(total.values(flow=True)), take(total.variances(flow=True))
    signal = {p: take(h.values(flow=True)) * ctx["xs"][p] for p, h in ctx["signal"].items()}
    metrics, worst, verdicts = _judge(signal, b, b_w2)
    worst_z = max(metrics, key=lambda p: metrics[p]["z_binned_all"])
    row = {"control sample": cr, "region": region, "points": len(signal),
           "B": metrics[worst]["B"], "neff": metrics[worst]["neff"], "worst point": worst,
           "S": metrics[worst]["S"], "contam": max(m["contam"] for m in metrics.values()),
           "z": max(m["z"] for m in metrics.values()),
           "z_binned": max(m["z_binned"] for m in metrics.values()),
           "z_binned_all": metrics[worst_z]["z_binned_all"], "worst point binned": worst_z,
           "s_sparse": verdicts["s_sparse"], "b_sparse": verdicts["b_sparse"]}
    for use in ("counting", "shape", "sparse bins"):
        row[use], row[f"{use} branch"] = verdicts[use]
    # stress test: the heaviest simulated background event of the region taken out of its bin
    weight, sample, i = _heaviest_event(ctx["background"], take)
    b_cut = b.copy()
    if sample is not None:
        b_cut[i] = max(b_cut[i] - weight, 0.0)
    cut = _judge(signal, b_cut, b_w2)[2]
    row["heaviest event"] = f"{sample}, weight {weight:.3g}" if sample else "none"
    row["B without it"] = float(b_cut.sum())
    for use in ("counting", "shape", "sparse bins"):
        row[f"{use}, stress test"] = cut[use][0]
    return row


def _taker(indices, axis):
    """A function mapping a flow=True array of a probe histogram to the flat bins of the given
    indices of one axis (None: every bin)"""
    def take(arr):
        return (arr if indices is None else np.take(arr, indices, axis=axis)).ravel()
    return take


VERDICT_KEYS = ("counting", "shape", "counting, stress test", "shape, stress test")
SPARSE_VERDICT_KEYS = ("sparse bins", "sparse bins, stress test")


def passes(row):
    """Whether every verdict of a row of lj_cr_safety is safe, the sparse-bin ones only with
    STRICT_SPARSE_BINS"""
    keys = VERDICT_KEYS + (SPARSE_VERDICT_KEYS if STRICT_SPARSE_BINS else ())
    return all(row[key] == "safe" for key in keys)


def lj_cr_safety(inputs, split_axis="probe_pt"):
    """The group's criteria applied to each control sample before any data, one row per region:
    the whole probe histogram, each bin of split_axis with the other axes kept (the regions a data
    step reads), and the union of the split-axis regions that are safe only through the S < 1
    rule, which per region does not bound their sum (the union check of the working-point decks).

    Each row has the background B and its effective number of simulated events, the worst signal
    point and its S, the worst S/(S+B), Z and binned Z over the 96 points of the limit grid at
    their model cross sections, and three verdicts: counting (integrated Z), shape (binned Z over
    the bins with background) and sparse bins (the bins expecting under one background event,
    judged together as one region; it decides only with STRICT_SPARSE_BINS), each also with the
    heaviest simulated background event removed from its bin (the stress test)."""
    # pylint: disable-next=import-outside-toplevel,import-error
    import plot_tools as pt
    xs = pt.theory_xs_fb()
    # the points with a model cross section and L_xy >= 3 cm; not limit_grid, which ranks the
    # lifetimes among the samples it is given, and the campaign holds only these points
    points = sorted(p for p in inputs["signal"] if p in xs and pt.parse_point(p)["lxy"] >= 3.0)
    if len(points) != 96:
        raise ValueError(f"{len(points)} signal points with a model cross section and "
                         "L_xy >= 3 cm, not the 96 of the limit grid")
    rows = []
    for cr in LJ_CR:
        ctx = {"background": {s: per[cr] for s, per in inputs["background"].items()},
               "signal": {p: inputs["signal"][p][cr] for p in points}, "xs": xs}
        axes = next(iter(ctx["background"].values())).axes
        k = [a.name for a in axes].index(split_axis)
        slices = _slices(axes[split_axis])
        mine = [_region_row(cr, "whole histogram", _taker(None, k), ctx)]
        mine += [_region_row(cr, f"{split_axis} {label}", _taker([index], k), ctx)
                 for label, index in slices]
        by_events = [index for (_, index), row in zip(slices, mine[1:])
                     if passes(row) and "events" in (row["counting branch"], row["shape branch"])]
        if len(by_events) > 1:
            mine.append(_region_row(cr, "union of the S < 1 regions", _taker(by_events, k), ctx))
        rows += mine
    return rows


def lj_cr_readable(rows, cr):
    """The split-axis regions of a control sample that may be read in data: those whose verdicts
    are all safe (passes), without the ones safe only through S < 1 when their union is not"""
    mine = [r for r in rows if r["control sample"] == cr]
    union = next((r for r in mine if r["region"] == "union of the S < 1 regions"), None)
    out = []
    for r in mine:
        if r["region"] in ("whole histogram", "union of the S < 1 regions") or not passes(r):
            continue
        if (union is not None and not passes(union)
                and "events" in (r["counting branch"], r["shape branch"])):
            continue
        out.append(r["region"])
    return out


def require_safe(rows, cr, regions):
    """Stop unless every region a data step reads is safe (lj_cr_readable): the data step of a
    control sample calls this before it opens anything derived from data"""
    allowed = set(lj_cr_readable(rows, cr))
    refused = [r for r in regions if r not in allowed]
    if refused or not regions:
        raise RuntimeError(f"{cr}: not safe to read in data: {refused or 'nothing requested'}; "
                           f"safe regions: {sorted(allowed)}")


# where 07_lj_reco_control_samples.ipynb writes the safety table and its figures; the data step
# reads the table back (lj_cr_verdicts) and refuses every region it does not mark safe
LJ_CR_RESULTS_DIR = "/store/user/scampbel/sidm_coffea_outputs/lj_reco_cr_v2_results"


def write_verdicts(rows, directory=LJ_CR_RESULTS_DIR):
    """Write the safety table of lj_cr_safety to directory as safety.json, with the switch it was
    judged with and the readable regions per control sample"""
    table = {"strict_sparse_bins": STRICT_SPARSE_BINS, "rows": rows,
             "readable": {cr: lj_cr_readable(rows, cr) for cr in LJ_CR}}
    fs, root = fsspec.core.url_to_fs(eos_url(directory))
    fs.makedirs(root, exist_ok=True)
    with fs.open(f"{root.rstrip('/')}/safety.json", "w") as f:
        json.dump(table, f, indent=1, default=float)
    fs.invalidate_cache()


def lj_cr_verdicts(directory=LJ_CR_RESULTS_DIR):
    """The safety table written by write_verdicts, checked to be judged with the current
    STRICT_SPARSE_BINS, so that a change of the rule forces the check to be redone"""
    with fsspec.open(f"{eos_url(directory)}/safety.json", "rt") as f:
        table = json.load(f)
    if table["strict_sparse_bins"] != STRICT_SPARSE_BINS:
        raise RuntimeError("the safety table was judged with STRICT_SPARSE_BINS = "
                           f"{table['strict_sparse_bins']}; rerun the safety check")
    return table["rows"]


def mass_window(h, axis, lo, hi):
    """h summed over the bins of axis between the bin edges lo and hi"""
    return h[{axis: slice(hist.loc(lo), hist.loc(hi), sum)}]


def probe_efficiency(h, spec, axis):
    """The passing fraction of the probes of one control sample (h: its probe histogram, summed
    over samples; spec: its LJ_CR entry) against one axis, with the background under the mass
    peak subtracted using the sidebands, scaled by the ratio of the widths (a linear background).
    Returns {"edges", "eff", "err", "probes", "passing"} over the bins of axis, its flow bins
    included at either end; err follows from the sums of squared weights, which are the counts
    in data."""
    (p_lo, p_hi), sidebands = spec["peak"], spec["sidebands"]
    k = (p_hi - p_lo) / sum(hi - lo for lo, hi in sidebands)
    peak = mass_window(h, spec["mass"], p_lo, p_hi).project(axis, "passing")
    side = sum(mass_window(h, spec["mass"], lo, hi).project(axis, "passing")
               for lo, hi in sidebands)
    values = peak.values(flow=True) - k * side.values(flow=True)
    variances = peak.variances(flow=True) + k * k * side.variances(flow=True)
    passing, probes = values[:, 1], values.sum(axis=1)
    var_pass, var_fail = variances[:, 1], variances[:, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        eff = np.where(probes > 0, passing / probes, np.nan)
        # passing and failing probes are independent: d(eff) = ((1-eff) dP, -eff dF) / N
        err = np.sqrt((1 - eff) ** 2 * var_pass + eff ** 2 * var_fail) / probes
    return {"axis": peak.axes[0], "edges": peak.axes[0].edges, "eff": eff, "err": err,
            "probes": probes, "passing": passing}


def efficiency_ratio(data, simulation):
    """The data/simulation ratio of two probe_efficiency results, with its uncertainty"""
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = data["eff"] / simulation["eff"]
        err = ratio * np.sqrt((data["err"] / data["eff"]) ** 2
                              + (simulation["err"] / simulation["eff"]) ** 2)
    return ratio, err


# the background groups of the control samples: those of corrections_study without DYJetsToMuMu,
# whose muon pairs DYJetsToLL holds
LJ_CR_GROUPS = {name: samples for name, samples in study.BACKGROUNDS.items()
                if not set(samples) & set(EXCLUDED_BACKGROUNDS)}


def in_range(axis, values):
    """values over the flow=True bins of a histogram axis, cut to its in-range bins"""
    start = 1 if axis.traits.underflow else 0
    return values[start:start + len(axis)]


def plot_probe_distribution(inputs, cr, axis, points=(), ax=None):
    """The probes of one control sample in simulation against one axis of its probe histogram,
    on a log scale: the background stacked by group (LJ_CR_GROUPS) and the given signal points at
    their model cross sections, in-range bins only"""
    # pylint: disable-next=import-outside-toplevel,import-error
    import plot_tools as pt
    xs = pt.theory_xs_fb()
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    bottom, edges = None, None
    colours = (study.BLUE, study.ORANGE, study.AQUA, study.GRAY)
    for (group, samples), colour in zip(LJ_CR_GROUPS.items(), colours):
        hists = [inputs["background"][s][cr] for s in samples if s in inputs["background"]]
        if not hists:
            continue
        proj = sum(hists).project(axis)
        edges, values = proj.axes[0].edges, proj.values()
        base = np.zeros_like(values) if bottom is None else bottom
        ax.stairs(base + values, edges, baseline=base, fill=True, color=colour, label=group)
        bottom = base + values
    for point, style in zip(points, ("-", "--", ":", "-.")):
        proj = inputs["signal"][point][cr].project(axis)
        ax.stairs(proj.values() * xs[point], proj.axes[0].edges, color=study.BLACK,
                  linestyle=style, linewidth=2, label=f"{point} (model)")
    ax.set_yscale("log")
    ax.set_ylim(bottom=1e-2)
    if edges is not None:
        ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel(inputs["background"][next(iter(inputs["background"]))][cr].axes[axis].label)
    ax.set_ylabel("probes")
    ax.legend(fontsize=11)
    return ax


def plot_probe_efficiency(results, ax=None):
    """Passing fractions against one axis ({label: probe_efficiency result}), in-range bins"""
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    for (label, r), colour in zip(results.items(), (study.BLACK, study.BLUE, study.ORANGE)):
        edges = r["edges"]
        centres = 0.5 * (edges[1:] + edges[:-1])
        ax.errorbar(centres, in_range(r["axis"], r["eff"]), yerr=in_range(r["axis"], r["err"]),
                    xerr=0.5 * np.diff(edges), fmt="o", color=colour, label=label)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel(next(iter(results.values()))["axis"].label)
    ax.set_ylabel("passing fraction")
    ax.legend(fontsize=11)
    return ax


# The merged data outputs of the control samples (campaign lj_reco_cr_data_v1, not yet run): the
# alone versions only, filled for the regions the safety table marks safe
LJ_CR_DATA_DIR = "/store/user/scampbel/sidm_coffea_outputs/lj_reco_cr_data_v1"


def readable_bins(rows, cr, axis):
    """Indices, in the flow=True layout of axis, of the bins of the split axis (probe pT) that the
    safety table marks readable for one control sample"""
    allowed = set(lj_cr_readable(rows, cr))
    return [index for label, index in _slices(axis)
            if f"{axis.name} {label}" in allowed]


def lj_cr_data(rows, directory=LJ_CR_DATA_DIR, split_axis="probe_pt"):
    """{control sample: probe histogram summed over the data samples} of the alone versions, with
    every bin of split_axis that the safety table (rows, from lj_cr_verdicts) does not mark
    readable set to zero before anything else is done with it; stops through require_safe unless
    at least one region is readable, and at any simulated sample"""
    out = {}
    for url in list_coffea(directory):
        for sample, sample_out in dt.read_coffea(url).items():
            if dt.is_simulation(sample, sample_out):
                raise ValueError(f"{url}: {sample} is simulation, not data")
            for cr, spec in LJ_CR.items():
                h = sample_out["hists"][spec["hist"]][{"channel": spec["channel"]}]
                out[cr] = h if cr not in out else out[cr] + h
    for cr, h in out.items():
        axis = h.axes[split_axis]
        keep = readable_bins(rows, cr, axis)
        require_safe(rows, cr, [f"{split_axis} {label}" for label, index in _slices(axis)
                                if index in keep])
        k = [a.name for a in h.axes].index(split_axis)
        view = h.view(flow=True)
        index = [slice(None)] * view.ndim
        index[k] = [i for i in range(view.shape[k]) if i not in keep]
        for field in ("value", "variance"):     # the weighted view takes no scalar as a whole
            view[field][tuple(index)] = 0.0
    return out
