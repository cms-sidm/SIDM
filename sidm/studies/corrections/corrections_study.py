"""Helpers shared by the notebooks of the corrections study

Only notebook-side code lives here: loading the merged outputs of the Condor
campaign from EOS (with the known generator-weight outliers taken out of the
normalisation), the sample lists, a quiet context for the tutorial's processor
run, and the figure conventions (colours, the CMS Simulation header, saving into
slides/figures for the deck). Nothing here runs inside a dask task, so it does
not need to be under sidm/tools (see README section 7 of the repo).
"""

import contextlib
import functools
import io
import os
import subprocess
import warnings

import hist
import mplhep as hep
import numpy as np
import rich
from coffea import util as coffea_util

from sidm import BASE_DIR
from sidm.tools import utilities

STUDY_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(STUDY_DIR, "outputs")              # *.coffea is gitignored
FIG_DIR = os.path.join(STUDY_DIR, "slides", "figures")      # *.pdf is gitignored

# One colour per role in every figure: the first three slots of the validated
# categorical palette, black for data, gray for uncorrected references.
BLUE, ORANGE, AQUA, GRAY, BLACK = "#2a78d6", "#eb6834", "#1baf7a", "#8c8c8c", "#0b0b0b"

LUMI_2018 = 59.8   # /fb, run_periods.yaml

REDIRECTOR = "root://cmseos.fnal.gov"

# The merged outputs of the Condor campaign (README, "The Condor campaign"): one
# <sample>.coffea and <sample>.meta.yaml per sample, readable by everyone. To use
# another production, set CAMPAIGN_DIR before the first load().
CAMPAIGN_DIR = "/store/user/scampbel/sidm_coffea_outputs/corrections_v1"


# Events with an unphysical generator weight that dominates a sample's sum of
# weights (README, "Found along the way"). In the DYJetsToMuMu_M50 skim one event
# has genWeight 1.45e14 (the LHE weight too), against ~2400 for the others: it
# alone is 99.8% of the sample's sum of weights, so lumi*xs/sumw, and with it every
# yield of the sample, comes out ~630 times too small. The event fails the trigger,
# so it is in no selected yield, and taking it out of the normalisation is exact.
GENWEIGHT_OUTLIERS = {
    "DYJetsToMuMu_M50": [{
        "file": "/store/group/lpcmetx/SIDM/Backgrounds/2018_v2/Skims/DYJetsToMuMu_M-50_massWgtFix/"
                "skimmed_output_6129.root",
        "event": 147047969,
        "genWeight": 145123002286080.0,
        # its event-level factors, for the corrected sums of weights
        "factors": {"pileup": 0.9187099917298068, "l1_prefiring": 0.998309850692749},
    }],
}


def drop_genweight_outliers(sample, out, skim_factor):
    """Take the GENWEIGHT_OUTLIERS of sample out of its normalisation, in place

    Every cutflow and histogram is scaled by sumw / (sumw - outliers), the
    normalisation without them, and the sums of weights in the metadata lose them.
    Valid only for outliers that fail every selection; the cutflow's "None" row,
    which contains them, is scaled too and should not be used.
    """
    outliers = GENWEIGHT_OUTLIERS.get(sample, [])
    if not outliers:
        return out
    meta = out["metadata"]
    removed = sum(o["genWeight"] for o in outliers) / skim_factor
    scale = meta["scaled_sum_weights"] / (meta["scaled_sum_weights"] - removed)
    for cf in out["cutflow"].values():
        cf.scale(scale)
    for name in out["hists"]:
        out["hists"][name] = out["hists"][name] * scale
    meta["scaled_sum_weights"] -= removed
    corrected = meta.get("scaled_sum_weights_corrected", {})
    for name in list(corrected):
        for o in outliers:
            factors = o["factors"]
            f = (np.prod(list(factors.values())) if name == "all_event_level"
                 else factors.get(name, 1.0))
            corrected[name] -= o["genWeight"] * f / skim_factor
    meta["genweight_outliers_dropped"] = len(outliers)
    return out


def load(samples, campaign_dir=None):
    """{sample: processor output} from the merged campaign outputs on EOS

    Each file is copied once with xrdcp into outputs/<campaign name>/ (gitignored)
    and read from there afterwards; delete that directory to fetch them again.
    A Kerberos ticket is enough to read them on LPC. Known generator-weight
    outliers are taken out of the normalisation (GENWEIGHT_OUTLIERS).
    """
    campaign_dir = (campaign_dir or CAMPAIGN_DIR).rstrip("/")
    local = os.path.join(CACHE_DIR, os.path.basename(campaign_dir))
    os.makedirs(local, exist_ok=True)
    out = {}
    for sample in samples:
        path = os.path.join(local, f"{sample}.coffea")
        if not os.path.exists(path):
            for ext in (".coffea", ".meta.yaml"):
                subprocess.run(["xrdcp", "-s", "-f", f"{REDIRECTOR}/{campaign_dir}/{sample}{ext}",
                                os.path.join(local, f"{sample}{ext}")], check=ext == ".coffea")
        merged = coffea_util.load(path)
        out[sample] = merged.get("out", merged)[sample]
        if sample in GENWEIGHT_OUTLIERS:
            drop_genweight_outliers(sample, out[sample], background_skim_factor(sample))
    return out


@functools.lru_cache(maxsize=None)
def background_skim_factor(sample):
    """skim_factor of a skimmed background, from configs/ntuples/backgrounds.yaml"""
    menu = utilities.load_yaml(os.path.join(BASE_DIR, "configs", "ntuples", "backgrounds.yaml"))
    return menu["skimmed_llpNanoAOD_v2"]["samples"][sample]["skim_factor"]


FINAL_STATES = {"2Mu2E": "signal_2mu2e_v10.yaml", "4Mu": "signal_4mu_v10.yaml"}
SR = {"2Mu2E": "2mu2e_sr", "4Mu": "4mu_sr"}
CATEGORY = {"2Mu2E": "2mu2e", "4Mu": "4mu"}

# background groups as in configs/composite_samples.yaml. The two DY samples
# overlap (DYJetsToLL contains the muon pairs of DYJetsToMuMu) and must not be
# added; which one the analysis uses is still open, so totals are made with each.
BACKGROUNDS = {
    r"DY$\to\mu\mu$": ["DYJetsToMuMu_M10to50", "DYJetsToMuMu_M50"],
    r"DY$\to\ell\ell$": ["DYJetsToLL_M10to50", "DYJetsToLL_M50"],
    "QCD": ["QCD_Pt15To20", "QCD_Pt20To30", "QCD_Pt30To50", "QCD_Pt50To80", "QCD_Pt80To120",
            "QCD_Pt120To170", "QCD_Pt170To300", "QCD_Pt300To470", "QCD_Pt470To600",
            "QCD_Pt600To800", "QCD_Pt800To1000", "QCD_Pt1000"],
    r"$t\bar{t}$": ["TTJets"],
    "diboson": ["WW", "WZ", "ZZ"],
}
TOTAL_BACKGROUND = {
    r"total, DY$\to\mu\mu$": [r"DY$\to\mu\mu$", "QCD", r"$t\bar{t}$", "diboson"],
    r"total, DY$\to\ell\ell$": [r"DY$\to\ell\ell$", "QCD", r"$t\bar{t}$", "diboson"],
}


def signal_points():
    """Every v10 signal point, 2Mu2E then 4Mu, in the order of the location YAMLs"""
    points = []
    for final_state, cfg in FINAL_STATES.items():
        menu = utilities.load_yaml(os.path.join(BASE_DIR, "configs", "ntuples", cfg))
        points += list(menu["llpNanoAOD_v2"]["samples"])
    return points


def parse_point(sample):
    """{"final_state", "m_phi", "m_zd", "ctau"} of a signal sample name such as 4Mu_500GeV_1p2GeV_1p9mm"""
    final_state, m_phi, m_zd, ctau = sample.split("_")
    number = lambda s, unit: float(s.replace(unit, "").replace("p", "."))
    return {"final_state": final_state, "m_phi": number(m_phi, "GeV"),
            "m_zd": number(m_zd, "GeV"), "ctau": number(ctau, "mm")}


def median_lxy():
    """{signal sample: median generator-level dark-photon Lxy in cm}, from the trigger study"""
    path = os.path.join(BASE_DIR, "studies", "truth_kinematics_forAN", "v3_trigger_efficiency.yml")
    return {s: v["median_genA_lxy_cm"] for s, v in utilities.load_yaml(path).items()}


def yields(out, channel):
    """Weighted yield after the last cut of channel, for nominal and every weight category"""
    return out["cutflow"][channel].variation_yields()


def raw(out, channel):
    """Raw number of simulated events after the last cut of channel"""
    rows = out["cutflow"][channel].rows
    return rows[list(rows)[-1]]["raw"]


def binned_mean(h):
    """Mean of a 1D histogram from its bin centres (flow bins excluded)"""
    values = h.values()
    return float(np.sum(h.axes[0].centers * values) / np.sum(values))


def campaign_samples(campaign_dir=None):
    """Names of the samples with a merged output in the campaign directory"""
    campaign_dir = (campaign_dir or CAMPAIGN_DIR).rstrip("/")
    listing = subprocess.run(["xrdfs", REDIRECTOR, "ls", campaign_dir], check=True,
                             capture_output=True, text=True).stdout.split()
    return sorted(os.path.basename(p)[:-len(".coffea")] for p in listing if p.endswith(".coffea"))


def cms_label(axes, lumi=None, fontsize=None):
    """CMS Simulation header on every panel of a finished figure

    Call once, after every panel and colorbar exists: the canvas is drawn so the
    layout settles, the layout is then frozen, and only then are the labels
    placed, as in truth_kinematics_forAN/an_style.py. lumi (in /fb) is shown for
    histograms normalised to the 2018 luminosity; otherwise only the year.
    """
    axes = list(np.ravel(np.asarray(axes, dtype=object)))
    fig = axes[0].figure
    if fontsize is None and len(axes) > 2:
        fontsize = 17
    fig.canvas.draw()
    if fig.get_layout_engine() is not None:
        fig.set_layout_engine("none")
    rlabel = rf"{lumi} fb$^{{-1}}$ (13 TeV)" if lumi else "2018 (13 TeV)"
    for ax in axes:
        hep.cms.label(ax=ax, data=False, rlabel=rlabel, fontsize=fontsize)


def save(fig, name):
    """Write a figure as vector PDF into slides/figures, where the deck reads it"""
    os.makedirs(FIG_DIR, exist_ok=True)
    path = os.path.join(FIG_DIR, f"{name}.pdf")
    fig.savefig(path, bbox_inches="tight")
    return path


def h1(values, axis, weights=None):
    """A weighted 1D hist.Hist filled with values"""
    h = hist.Hist(axis, storage=hist.storage.Weight())
    h.fill(np.asarray(values), weight=None if weights is None else np.asarray(weights))
    return h


def fractions(h):
    """h divided by its total, flow included: the fraction of entries per bin

    Use this instead of density=True on log-binned axes, where a density divides
    by the linear bin width and distorts the comparison of shapes.
    """
    return h / h.sum(flow=True).value


def band(ax, upper, lower, color, *, label=None, alpha=0.25):
    """Shade between two histograms with the same binning, e.g. up and down"""
    up, down = upper.values(), lower.values()
    ax.stairs(np.maximum(up, down), upper.axes[0].edges, baseline=np.minimum(up, down),
              fill=True, color=color, alpha=alpha, lw=0, label=label)


@contextlib.contextmanager
def quiet():
    """Keep the processor's progress bars and pre-existing warnings out of a notebook

    coffea draws its progress bar with rich, which in Jupyter installs its own stdout
    proxy, so a plain redirect_stdout does not catch the processor's prints. Pointing
    rich at a non-Jupyter console for the duration fixes that.
    """
    buf = io.StringIO()
    rich.reconfigure(force_jupyter=False)
    try:
        with contextlib.redirect_stdout(buf), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            yield buf
    finally:
        rich.reconfigure()
