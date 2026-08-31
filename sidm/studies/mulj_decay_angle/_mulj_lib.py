"""Helpers for the muon lepton-jet decay-angle study.

`mulj_decay_angle.ipynb` is a thin driver: every figure and every number in it comes from a
function defined here, so a single figure can be regenerated on its own after the inputs
change.

Inputs are the per-sample `.coffea` files written by the study's Condor campaign, one
merged file per sample.  Point `INPUT_DIR` at the directory holding them -- a local path,
or an EOS `root://` URL, which is copied to a local cache on first use -- and the sample
lists follow from what is actually there.  Nothing else in the notebook needs editing.

Normalisation, in one place because it is easy to get wrong: the campaign ran with
`--unweighted-hist`, so the histograms are RAW COUNTS with no cross-section or generator
weight applied.  The cutflows stored in the same files ARE luminosity x cross-section
scaled.  The two must never be combined.  `lumi_xs_weight()` implements the per-sample MC
rescaling (lumi x sigma / N_generated) for when a normalised MC prediction is wanted; the
shape comparisons in this study are unit-area and do not use it.

Figure conventions implemented by `Canvas`:
  * figures are designed for placement at 7.5 in text width, so every piece of in-figure
    text is sized as `base_pt / print_scale` with `print_scale = 7.5 / figure_width`;
  * a grid is never more than 1.15x taller than it is wide;
  * the CMS label is drawn exactly once, by `Canvas.finish()`, after all panels, legends
    and colour bars exist -- that call also freezes the constrained layout;
  * legends take an explicit location, and get an opaque backing when they sit over data;
  * efficiency points are drawn only where the denominator holds at least 20 entries;
  * log axes get explicit minor ticks.
"""

from __future__ import annotations

import glob
import os
import re
import subprocess

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter
import mplhep as hep
import coffea.util
import hist
from hist.intervals import clopper_pearson_interval, poisson_interval

try:  # the framework is only needed for the cross-section / luminosity lookup
    from sidm.tools import utilities as sidm_utilities
except ImportError:  # pragma: no cover - the study is still readable without it
    sidm_utilities = None

import yaml


# --------------------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------------------

#: width, in inches, at which every figure is placed in the analysis note.
TEXT_WIDTH_IN = 7.5

#: smallest denominator for which an efficiency point is drawn.
MIN_DENOMINATOR = 20

#: smallest number of entries for which a shape is drawn at all.  A curve can clear the
#: per-bin rule below simply by being drawn in few enough bins, so a shape also has to
#: rest on this many entries however coarsely it is binned.
MIN_SHAPE_ENTRIES = 50

#: smallest mean number of effective entries per bin for which a shape is drawn.  Below
#: about five the Poisson error on a bin is close to half its content and the curve
#: carries no shape information, only noise.  The count is the effective one, so a
#: cross-section-weighted sum of simulated samples is judged on the statistics it really
#: has rather than on the number of generated events behind it.
MIN_ENTRIES_PER_BIN = 5

_HERE = os.path.dirname(os.path.abspath(__file__))
_CONFIG_DIR = os.path.normpath(os.path.join(_HERE, "..", "..", "configs"))

INPUT_DIR = "/uscms_data/d3/murtazas/mulj_smoke"
LUMI_FB = None
DATASET_TAG = ""
YEAR = "2018"
FIG_DIR = os.path.join(_HERE, "figures")
CACHE_DIR = os.path.join(os.path.expanduser("~"), ".cache", "mulj_decay_angle")

_LOADED: dict = {}
_NOTES: list = []
_UNSET = object()


def configure(input_dir=None, lumi_fb=_UNSET, dataset_tag=None, year=None, fig_dir=None,
              cache_dir=None):
    """Set the module-level configuration.

    `input_dir` is the only setting that must change when the smoke-test inputs are
    swapped for the merged campaign outputs.  `lumi_fb` controls the luminosity printed in
    the CMS label of data figures; leave it as None while running on a partial sample, so
    that no figure claims a luminosity it does not have.
    """
    global INPUT_DIR, LUMI_FB, DATASET_TAG, YEAR, FIG_DIR, CACHE_DIR
    if input_dir is not None:
        INPUT_DIR = input_dir.rstrip("/")
        _LOADED.clear()
    if lumi_fb is not _UNSET:
        LUMI_FB = None if not lumi_fb else float(lumi_fb)
    if dataset_tag is not None:
        DATASET_TAG = dataset_tag
    if year is not None:
        YEAR = str(year)
    if fig_dir is not None:
        FIG_DIR = fig_dir
    if cache_dir is not None:
        CACHE_DIR = cache_dir
    os.makedirs(FIG_DIR, exist_ok=True)
    return {"input_dir": INPUT_DIR, "lumi_fb": LUMI_FB, "dataset_tag": DATASET_TAG,
            "year": YEAR, "fig_dir": FIG_DIR}


def note(msg):
    """Record and print a note about something that could not be drawn.

    A message that has already been recorded is counted but not printed again, so that a
    grid of panels sharing one empty input does not bury the rest of the output.
    """
    if msg in _NOTES:
        return
    _NOTES.append(msg)
    print(f"  [note] {msg}")


def notes():
    """Every note recorded so far, in order."""
    return list(_NOTES)


# --------------------------------------------------------------------------------------
# regions
# --------------------------------------------------------------------------------------

#: the five regions this study compares, in the order used across the notebook.
REGIONS = [
    "data_control_region_1muLj",
    "data_control_region_1muLj_cosmic_veto",
    "data_control_region_1muLj_spread_cosAlpha_mu_veto",
    "test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto",
    "test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto",
]

#: compact two-line names used as legend titles inside panels.
REGION_SHORT = {
    "data_control_region_1muLj": "1 $\\mu$-LJ CR\nno cleaning",
    "data_control_region_1muLj_cosmic_veto": "1 $\\mu$-LJ CR\ncosmic veto",
    "data_control_region_1muLj_spread_cosAlpha_mu_veto": "1 $\\mu$-LJ CR\nspread + $\\cos\\alpha$ veto",
    "test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto": "2$\\mu$2e VR, inv. displaced\nPF-muon tagged",
    "test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto": "4$\\mu$ VR, inv. displaced\nPF-muon tagged",
}

#: single-line names for tables and prints.
REGION_LABEL = {
    "data_control_region_1muLj": "1 mu-LJ CR (no cleaning)",
    "data_control_region_1muLj_cosmic_veto": "1 mu-LJ CR + cosmic veto",
    "data_control_region_1muLj_spread_cosAlpha_mu_veto": "1 mu-LJ CR + spread and cos-alpha veto",
    "test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto": "2mu2e inverted-displacement VR",
    "test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto": "4mu inverted-displacement VR",
}

#: regions whose inverted-displacement requirement can only be satisfied by a PF muon, so
#: that pairs containing displaced-standalone muons are absent by construction.
PF_TAGGED_REGIONS = {
    "test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto",
    "test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto",
}

#: the only channels collision data may be read in: this study's control and validation
#: regions.  Every other channel is the signal region or a variant of it and stays blind.
DATA_CHANNELS = frozenset({
    "data_control_region_1muLj",
    "data_control_region_1muLj_cosmic_veto",
    "data_control_region_1muLj_spread_cosAlpha_mu_veto",
    "data_control_region_1egmLj_cosmic_veto",
    "cr_2lj_invDisplaced_spread_cosAlpha_mu_veto",
    "test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto",
    "test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto",
})


def refuse_blinded_data(samples, channel):
    """Raise when collision-data samples are asked for in a channel outside DATA_CHANNELS.

    Called before anything is read, so a signal-region channel is never opened on data,
    whatever the caller passes.
    """
    if channel in DATA_CHANNELS:
        return
    opened = [name for name in samples if parse_data(name) is not None]
    if opened:
        raise ValueError(
            f"channel {channel!r} is not one of this study's control or validation regions "
            f"and {len(opened)} collision-data sample(s) were passed to it. Data is read "
            "only in: " + ", ".join(sorted(DATA_CHANNELS)))


#: region used as the source of the fake sample in the correlation and cut-scan sections.
UNCUT_CR = "data_control_region_1muLj"
CLEANED_CR = "data_control_region_1muLj_spread_cosAlpha_mu_veto"

#: signal channels, one per signal final state.
SIGNAL_CHANNEL = {"4Mu": "4mu", "2Mu2E": "2mu2e"}

MUON_MASS = 0.1056583745  # GeV


# --------------------------------------------------------------------------------------
# sample bookkeeping
# --------------------------------------------------------------------------------------

_SIGNAL_RE = re.compile(r"^(4Mu|2Mu2E)_(\d+)GeV_([0-9p]+)GeV_([0-9p]+)mm$")
_DATA_RE = re.compile(r"^(DoubleMuon|SingleMuon|EGamma)_(\d{4})([A-Z])(?:_(\d+))?$")


def _p2f(text):
    """Convert the `1p2` style number used in sample names to a float."""
    return float(text.replace("p", "."))


def parse_signal(name):
    """Decompose a signal sample name, or return None if it is not a signal sample."""
    match = _SIGNAL_RE.match(name)
    if match is None:
        return None
    return {
        "name": name,
        "final_state": match.group(1),
        "m_xx": float(match.group(2)),
        "m_zd": _p2f(match.group(3)),
        "ctau_mm": _p2f(match.group(4)),
        "channel": SIGNAL_CHANNEL[match.group(1)],
    }


def parse_data(name):
    """Decompose a collision-data sample name, or return None."""
    match = _DATA_RE.match(name)
    if match is None:
        return None
    return {"name": name, "primary_dataset": match.group(1), "year": match.group(2),
            "era": match.group(3), "part": match.group(4)}


def _remote_listing(url):
    server, _, path = url[len("root://"):].partition("/")
    path = "/" + path.lstrip("/")
    out = subprocess.run(["xrdfs", server, "ls", path], capture_output=True, text=True,
                         check=True)
    return [os.path.basename(line) for line in out.stdout.split() if line.endswith(".coffea")]


def list_samples():
    """Sample names available in `INPUT_DIR`, sorted."""
    if INPUT_DIR.startswith("root://"):
        files = _remote_listing(INPUT_DIR)
    else:
        files = [os.path.basename(p) for p in glob.glob(os.path.join(INPUT_DIR, "*.coffea"))]
    return sorted(f[: -len(".coffea")] for f in files)


def signal_samples(final_state=None, m_zd=None, m_xx=None):
    """Signal samples present in `INPUT_DIR`, optionally filtered.

    All lifetimes of the selected mass point are returned; summing them is the caller's
    choice (see `sum_over` in `get_hist`).
    """
    picked = []
    for name in list_samples():
        info = parse_signal(name)
        if info is None:
            continue
        if final_state is not None and info["final_state"] != final_state:
            continue
        if m_zd is not None and not np.isclose(info["m_zd"], m_zd):
            continue
        if m_xx is not None and not np.isclose(info["m_xx"], m_xx):
            continue
        picked.append(name)
    return picked


def data_samples(era=None):
    """Collision-data samples present in `INPUT_DIR`; all eras unless one is named."""
    picked = []
    for name in list_samples():
        info = parse_data(name)
        if info is None:
            continue
        if era is not None and info["era"] != era:
            continue
        picked.append(name)
    return picked


def background_samples():
    """Simulated-background samples present in `INPUT_DIR`.

    A file only counts as a background if a cross section is configured for it: that is
    what a normalised background prediction needs, and it keeps stray files in the input
    directory out of the sum.  Anything skipped is reported once.
    """
    table = load_cross_sections()
    picked, ignored = [], []
    for name in list_samples():
        if parse_signal(name) is not None or parse_data(name) is not None:
            continue
        if name in table:
            picked.append(name)
        else:
            ignored.append(name)
    if ignored:
        note("no cross section configured for " + ", ".join(sorted(ignored)) +
             "; not treated as background")
    return picked


def signal_grid():
    """Mass points present in `INPUT_DIR`, as a sorted list of (final_state, m_zd, m_xx)."""
    points = set()
    for name in signal_samples():
        info = parse_signal(name)
        points.add((info["final_state"], info["m_zd"], info["m_xx"]))
    return sorted(points)


# --------------------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------------------

def _local_path(name):
    """Local path of a sample file, fetching it from EOS on first use if needed."""
    if not INPUT_DIR.startswith("root://"):
        return os.path.join(INPUT_DIR, f"{name}.coffea")
    os.makedirs(CACHE_DIR, exist_ok=True)
    local = os.path.join(CACHE_DIR, f"{name}.coffea")
    if not os.path.exists(local):
        subprocess.run(["xrdcp", "-f", f"{INPUT_DIR}/{name}.coffea", local], check=True)
    return local


def load_sample(name):
    """Load one sample and return `{hists, cutflow, metadata, is_data, sum_weights, year}`.

    Returns None, with a note, when the file is missing -- so that a notebook run against
    a partial set of inputs degrades instead of raising.
    """
    if name in _LOADED:
        return _LOADED[name]
    try:
        path = _local_path(name)
        blob = coffea.util.load(path)
    except (FileNotFoundError, OSError, subprocess.CalledProcessError) as exc:
        note(f"sample {name} could not be read ({type(exc).__name__}); skipping it")
        _LOADED[name] = None
        return None
    out = blob["out"] if "out" in blob else blob
    if name in out:
        payload = out[name]
    elif len(out) == 1:
        payload = next(iter(out.values()))
    else:
        note(f"file for {name} holds {len(out)} samples and none matches the name; skipping")
        _LOADED[name] = None
        return None
    meta = payload["metadata"]
    record = {
        "hists": payload["hists"],
        "cutflow": payload.get("cutflow", {}),
        "metadata": meta,
        "is_data": bool(_scalar(meta.get("is_data"), False)),
        "year": str(_scalar(meta.get("year"), YEAR)),
        "sum_weights": float(meta.get("scaled_sum_weights", meta.get("n_evts", 0.0))),
        "n_evts": float(meta.get("n_evts", 0.0)),
        "unweighted_hist": bool(_scalar(meta.get("unweighted_hist"), True)),
    }
    _LOADED[name] = record
    return record


def _scalar(value, default):
    """Unwrap the single-element set accumulators coffea stores metadata in."""
    if value is None:
        return default
    if isinstance(value, (set, frozenset)):
        return next(iter(value)) if value else default
    try:
        if hasattr(value, "__iter__") and not isinstance(value, (str, bytes)):
            items = list(value)
            return items[0] if items else default
    except TypeError:
        pass
    return value


def load_cross_sections():
    """The cross-section table used for MC normalisation, in pb."""
    with open(os.path.join(_CONFIG_DIR, "cross_sections.yaml"), encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_run_periods():
    """The run-period table, which carries the integrated luminosity per year."""
    with open(os.path.join(_CONFIG_DIR, "run_periods.yaml"), encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def lumi_xs_weight(sample, year=None):
    """Per-sample MC weight lumi x sigma / N_generated.

    IMPORTANT: with `--unweighted-hist` the stored histograms are raw counts, so this
    factor is what turns them into an expected yield.  The denominator taken here is the
    `scaled_sum_weights` entry of the sample metadata; for signal samples that quantity
    counts the events that survived the skim, not the events that were generated, so a
    signal yield built from it is too large by the skim acceptance and the generated count
    from the file census has to be substituted before any yield is quoted.  Shapes -- what
    this study actually plots -- are unaffected, because the factor cancels in a unit-area
    normalisation.
    """
    record = load_sample(sample)
    if record is None:
        return None
    if record["is_data"]:
        return 1.0
    year = year or record["year"] or YEAR
    if sidm_utilities is not None:
        try:
            return sidm_utilities.get_lumixs_weight(sample, year, record["sum_weights"])
        except (KeyError, FileNotFoundError):
            pass
    table = load_cross_sections()
    if sample in table:
        xsec = table[sample]
    elif sample.startswith(("4Mu", "2Mu2E")):
        xsec = 0.001  # 1 fb is the convention for the signal grid
    else:
        note(f"no cross section for {sample}; MC weight set to 1")
        return 1.0
    lumi = load_run_periods()[year]["lumi"]
    return lumi * xsec / record["sum_weights"] if record["sum_weights"] else 1.0


# --------------------------------------------------------------------------------------
# histogram access
# --------------------------------------------------------------------------------------

def _scale_hist(histogram, factor):
    scaled = histogram.copy()
    view = scaled.view(flow=True)
    if view.dtype.names and "value" in view.dtype.names:
        view["value"] *= factor
        view["variance"] *= factor * factor
    else:
        scaled *= factor
    return scaled


def get_hist(samples, hist_name, channel, weighted=False, quiet=False):
    """Sum one histogram over samples, with the `channel` axis already selected.

    `samples` is a name or a list of names -- passing every lifetime of a mass point sums
    over lifetimes, passing every data era sums over eras.  `weighted=True` applies the
    per-sample lumi x sigma / N factor before summing, which is what a background stack
    needs; the default leaves raw counts, which is what every statistical statement in
    this study is built from.  Returns None, with a note, when nothing is available.
    Collision data is refused outside the control and validation regions
    (`refuse_blinded_data`).
    """
    if isinstance(samples, str):
        samples = [samples]
    refuse_blinded_data(samples, channel)
    total = None
    for name in samples:
        record = load_sample(name)
        if record is None:
            continue
        histogram = record["hists"].get(hist_name)
        if histogram is None:
            if not quiet:
                note(f"{name} has no histogram {hist_name}; skipping it")
            continue
        available = list(histogram.axes["channel"])
        if channel not in available:
            if not quiet:
                note(f"{name} has no channel {channel} for {hist_name}; skipping it")
            continue
        selected = histogram[{"channel": channel}]
        if weighted:
            factor = lumi_xs_weight(name)
            if factor is None:
                continue
            selected = _scale_hist(selected, factor)
        else:
            selected = selected.copy()
        total = selected if total is None else total + selected
    if total is None and not quiet:
        note(f"no input for {hist_name} in channel {channel} from {len(samples)} sample(s)")
    return total


def raw_counts(histogram, flow=False):
    """Unweighted entry count per bin, recovered as value^2 / variance.

    For the raw-count histograms this study reads, value and variance are equal and this
    returns the counts unchanged; it stays correct if a globally rescaled histogram is
    ever passed in.
    """
    values = np.asarray(histogram.values(flow=flow), dtype=float)
    variances = histogram.variances(flow=flow)
    if variances is None:
        return values
    variances = np.asarray(variances, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        counts = np.where(variances > 0, values * values / variances, 0.0)
    return np.nan_to_num(counts)


def total_entries(histogram):
    """Total unweighted entries in a histogram, overflow included."""
    if histogram is None:
        return 0.0
    return float(raw_counts(histogram, flow=True).sum())


def unit_area(histogram, axis=0):
    """Bin centres, unit-area density and its uncertainty for a 1D histogram.

    The density integrates to one over the axis range: dividing by the bin width as well
    as by the total keeps distributions with different binnings directly comparable.
    Also returns the bin widths and the total, which a caller needs to put an
    independently computed error bar on the same scale.
    """
    edges = histogram.axes[axis].edges
    widths = np.diff(edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    values = np.asarray(histogram.values(), dtype=float)
    variances = histogram.variances()
    variances = values if variances is None else np.asarray(variances, dtype=float)
    total = float(values.sum())
    if total <= 0:
        return centres, np.zeros_like(values), np.zeros_like(values), widths, 0.0
    density = values / widths / total
    errors = np.sqrt(variances) / widths / total
    return centres, density, errors, widths, total


def poisson_errors(counts):
    """Asymmetric Poisson (Garwood) error bars for observed counts."""
    counts = np.asarray(counts, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        interval = poisson_interval(counts, counts)
    lower = np.nan_to_num(interval[0], nan=0.0)
    upper = np.nan_to_num(interval[1], nan=0.0)
    lower = np.clip(counts - lower, 0.0, None)
    upper = np.clip(upper - counts, 0.0, None)
    # a zero-count bin has no lower error and an upper edge at 1.84 events
    upper = np.where(counts == 0, 1.841, upper)
    return np.vstack([lower, upper])


def clopper_pearson(passed, total, coverage=0.6827):
    """Clopper-Pearson interval, returned as (efficiency, lower error, upper error)."""
    passed = np.atleast_1d(np.asarray(passed, dtype=float))
    total = np.atleast_1d(np.asarray(total, dtype=float))
    efficiency = np.divide(passed, total, out=np.zeros_like(passed),
                           where=total > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        interval = clopper_pearson_interval(passed, np.where(total > 0, total, 1.0),
                                            coverage=coverage)
    lower = np.nan_to_num(interval[0], nan=0.0)
    upper = np.nan_to_num(interval[1], nan=0.0)
    return efficiency, np.clip(efficiency - lower, 0, None), np.clip(upper - efficiency, 0, None)


def efficiency_series(passed, total, min_denominator=MIN_DENOMINATOR):
    """Efficiency with Clopper-Pearson errors, masked where the denominator is too small.

    Points with fewer than `min_denominator` entries in the denominator are dropped rather
    than drawn: a Clopper-Pearson interval on one or two events spans essentially the whole
    range and only adds a picket fence to the figure.
    """
    efficiency, low, high = clopper_pearson(passed, total)
    keep = np.asarray(total, dtype=float) >= min_denominator
    return efficiency, low, high, keep


def profile_mean(histogram, value_axis=0, slice_axis=1, min_entries=MIN_DENOMINATOR):
    """Mean of the value axis in each bin of the companion axis, with its standard error.

    Returns companion bin centres, the mean, its uncertainty and a mask of the bins that
    hold at least `min_entries` entries.
    """
    counts = raw_counts(histogram)
    if value_axis == 1:
        counts = counts.T
    value_centres = histogram.axes[value_axis].centers
    companion = histogram.axes[slice_axis]
    per_bin = counts.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = (counts * value_centres[:, None]).sum(axis=0) / np.where(per_bin > 0, per_bin, 1)
        second = (counts * (value_centres[:, None] ** 2)).sum(axis=0) / np.where(per_bin > 0, per_bin, 1)
    variance = np.clip(second - mean ** 2, 0, None)
    error = np.sqrt(variance / np.where(per_bin > 0, per_bin, 1))
    keep = per_bin >= min_entries
    return _axis_centres(companion), np.nan_to_num(mean), np.nan_to_num(error), keep, per_bin


def _axis_centres(axis):
    """Bin centres of an axis, geometric for a wide variable binning, index for categories."""
    try:
        edges = np.asarray(axis.edges, dtype=float)
    except (AttributeError, TypeError):
        return np.arange(len(axis), dtype=float)
    lo, hi = edges[:-1], edges[1:]
    use_geometric = np.all(lo >= 0) and (hi[-1] / max(hi[0], 1e-12) > 50)
    if use_geometric:
        # a bin that starts at zero has no geometric centre; its arithmetic centre still
        # sits inside the bin and is representable on a log axis
        centres = np.sqrt(np.where(lo > 0, lo, hi) * hi)
        return np.where(lo > 0, centres, 0.5 * hi)
    return 0.5 * (lo + hi)


def beta_star(pair_mass, muon_mass=MUON_MASS):
    """Muon speed in the rest frame of a pair of the given invariant mass."""
    ratio = (2.0 * muon_mass / pair_mass) ** 2
    return float(np.sqrt(max(0.0, 1.0 - ratio)))


def ridge_ratio(cos_theta, pair_mass, muon_mass=MUON_MASS):
    """Transverse-momentum ratio expected for a boosted two-body decay.

    For a pair produced with a large boost, the two muons carry lab momenta in the ratio
    (1 - beta* cos theta*) / (1 + beta* cos theta*), which is where the diagonal band in
    the (|cos theta*|, pT ratio) plane comes from.
    """
    speed = beta_star(pair_mass, muon_mass)
    cos_theta = np.asarray(cos_theta, dtype=float)
    return (1.0 - speed * cos_theta) / (1.0 + speed * cos_theta)


def spin_one_density(cos_theta):
    """Unit-area |cos theta*| density of a transversely unpolarised 1 + cos^2 decay."""
    return 0.75 * (1.0 + np.asarray(cos_theta, dtype=float) ** 2)


# --------------------------------------------------------------------------------------
# figure scaffolding
# --------------------------------------------------------------------------------------

def set_style(dpi=90):
    """Apply the CMS plotting style used across the analysis."""
    plt.style.use(hep.style.CMS)
    plt.rcParams["figure.dpi"] = dpi
    plt.rcParams["savefig.dpi"] = dpi
    plt.rcParams["figure.constrained_layout.use"] = False


#: colours used consistently across the study
COLOURS = {
    "4Mu": "#3f6fd6",
    "2Mu2E": "#d55e00",
    "signal": "#3f6fd6",
    "background": "#009e73",
    "data": "black",
    "gen": "#7b3294",
    "reco": "#3f6fd6",
    "reference": "#999999",
}

_SLICE_COLOURS = ["#264653", "#2a9d8f", "#e9c46a", "#e76f51", "#8e7dbe"]


class Canvas:
    """A figure laid out to the analysis-note figure standards.

    Construct it with the panel grid and the size of one panel; the class derives the
    print scale from the resulting figure width and exposes `font`, the design font size
    that prints at `base_pt`.  Call `finish()` once, after every panel is complete, and
    then `save()`.
    """

    def __init__(self, nrows=1, ncols=1, panel_w=5.0, panel_h=4.0, base_pt=9.0,
                 sharex=False, sharey=False, width_ratios=None, height_ratios=None):
        self.width = panel_w * ncols
        self.height = panel_h * nrows
        self.nrows = nrows
        self.ncols = ncols
        self.print_scale = TEXT_WIDTH_IN / self.width
        self.font = base_pt / self.print_scale
        #: a colour bar's ticks may be smaller than the rest, but never below eight
        #: points once the figure is placed at text width
        self.tick_font = max(self.font * 0.9, 8.0 / self.print_scale)
        self.fig, axes = plt.subplots(
            nrows, ncols, figsize=(self.width, self.height), squeeze=False,
            sharex=sharex, sharey=sharey, layout="constrained",
            gridspec_kw={"width_ratios": width_ratios, "height_ratios": height_ratios},
        )
        self.axes = axes
        self._finished = False
        self._colourbar_axes = []
        for ax in self.axes.flat:
            ax.tick_params(labelsize=self.font, which="major")
            ax.tick_params(which="minor", labelsize=self.font * 0.85)
        # reserve one text line above the top row for the CMS label
        for ax in self.axes[0]:
            ax.set_title(" ", fontsize=self.font)
        if self.height > 1.15 * self.width:
            note(f"figure is {self.height:.1f} x {self.width:.1f} in, taller than 1.15x its "
                 "width; split it along a physics axis instead of shrinking the panels")

    # -- panel helpers ------------------------------------------------------------------

    def ax(self, row=0, col=0):
        """The axes at a grid position."""
        return self.axes[row][col]

    def labels(self, ax, xlabel=None, ylabel=None):
        """Set axis titles at the size that prints legibly at text width."""
        if xlabel is not None:
            ax.set_xlabel(xlabel, fontsize=self.font * 1.1)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=self.font * 1.1)

    def figure_ylabel(self, text):
        """One y-axis title for the whole grid.

        A per-panel title longer than its own panel overruns into the rows above and
        below it and lands on their tick labels, so a quantity shared by every row of a
        tall grid gets a figure-level title instead.  Call it before `finish()`.
        """
        self.fig.supylabel(text, fontsize=self.font * 1.1)

    def legend(self, ax, loc, handles=None, labels=None, title=None, ncol=1,
               over_data=False, **kwargs):
        """Add a legend at an explicit location.

        Matplotlib's automatic placement routinely lands a long legend in the middle of
        the data, so `loc` is required.  `over_data=True` gives the legend an opaque
        backing, which the CMS style does not do by default because it turns the frame
        off entirely.
        """
        style = dict(loc=loc, ncol=ncol, fontsize=self.font,
                     title_fontsize=self.font, borderpad=0.4, labelspacing=0.35,
                     handlelength=1.4, handletextpad=0.5, columnspacing=1.0)
        if over_data:
            style.update(frameon=True, facecolor="white", edgecolor="none", framealpha=0.92)
        style.update(kwargs)
        if handles is not None:
            legend = ax.legend(handles, labels, title=title, **style)
        else:
            legend = ax.legend(title=title, **style)
        if title is not None and legend.get_title() is not None:
            legend.get_title().set_fontsize(self.font)
        return legend

    def stamp(self, ax, text, loc="upper left", **kwargs):
        """A short opaque-backed annotation inside a panel."""
        # the box is padded, and the ticks point inwards, so the anchor has to leave
        # room for both: at 0.92 with a small pad the box clears the innermost tick
        positions = {"upper left": (0.03, 0.92, "left", "top"),
                     "upper right": (0.97, 0.92, "right", "top"),
                     "lower left": (0.03, 0.08, "left", "bottom"),
                     "lower right": (0.97, 0.08, "right", "bottom")}
        x, y, ha, va = positions[loc]
        return ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va,
                       fontsize=self.font,
                       bbox=dict(facecolor="white", edgecolor="none", alpha=0.92,
                                 boxstyle="round,pad=0.15"), **kwargs)

    def headroom(self, ax, factor=1.6, bottom=0.0):
        """Leave empty space above the drawn curves so a legend never sits on the data."""
        ax.autoscale_view()
        top = ax.dataLim.y1
        if not np.isfinite(top) or top <= 0:
            return
        ax.set_ylim(bottom, top * factor)

    def clear_legends(self, entries, pad=0.05, most=2.0):
        """Raise the top of each axis until its legend no longer covers the curves.

        `headroom` guesses how much empty space a legend needs; this measures it.  Every
        legend is anchored to the top of its own axes, so its height as a fraction of
        that axes does not change when the limit does, and one draw is enough to measure
        the whole grid.  `entries` holds `(axes, legend)`, or `(axes, legend, height)`
        where `height` is the height of the curves under that legend alone: a key placed
        in a corner the curves leave empty need not clear the tallest curve in the panel,
        and forcing it to would set the scale of every panel sharing the axis.
        """
        self.fig.canvas.draw()
        for entry in entries:
            ax, legend = entry[0], entry[1]
            reference = entry[2] if len(entry) > 2 else None
            box = legend.get_window_extent()
            low = ax.transAxes.inverted().transform((box.x0, box.y0))[1]
            bottom, top = ax.get_ylim()
            data_top = float(ax.dataLim.y1 if reference is None else reference)
            if not np.isfinite(data_top) or data_top <= bottom:
                continue
            room = max(low - pad, 0.05)
            # capped, so that a tall key cannot squeeze the curves into a strip at the
            # bottom of the panel; a key that still does not fit wants fewer lines
            needed = min(bottom + (data_top - bottom) / room,
                         bottom + (data_top - bottom) * most)
            if needed > top:
                ax.set_ylim(bottom, needed)

    def share_ylimits(self, entries, axes=None):
        """Put every panel of a grid on one vertical scale.

        The panels of a correlation figure all show the same normalised quantity, so a
        per-panel scale invites the reader to compare shapes that are drawn at different
        magnifications.  `entries` may be a list of axes or of `(axes, legend)` pairs.
        """
        panels = list(axes) if axes is not None else [
            item[0] if isinstance(item, tuple) else item for item in entries]
        if not panels:
            return
        bottom = min(a.get_ylim()[0] for a in panels)
        top = max(a.get_ylim()[1] for a in panels)
        for a in panels:
            a.set_ylim(bottom, top)

    def log_axis(self, ax, which="x"):
        """Turn on a log scale with explicit minor ticks."""
        locator = LogLocator(base=10.0, subs=tuple(np.arange(2, 10) * 0.1 * 10), numticks=100)
        if "x" in which:
            ax.set_xscale("log")
            ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
            ax.xaxis.set_minor_formatter(NullFormatter())
        if "y" in which:
            ax.set_yscale("log")
            ax.yaxis.set_minor_locator(locator)
            ax.yaxis.set_minor_formatter(NullFormatter())
        return ax

    def colorbar(self, ax, mappable, label=None):
        """A colour bar sized to print legibly; add it before `finish()` freezes the layout."""
        cbar = self.fig.colorbar(mappable, ax=ax, pad=0.02, fraction=0.05)
        cbar.ax.tick_params(labelsize=self.tick_font)
        if label:
            cbar.set_label(label, fontsize=self.font)
        self._colourbar_axes.append(ax)
        return cbar

    def shared_colourbar(self, axes, meshes, label=None):
        """One colour bar, on one scale, for a row of maps.

        Side-by-side maps on silently different colour scales cannot be compared cell by
        cell, which is the comparison the reader will make, so the row is put on the
        common scale of its largest cell and gets a single bar.
        """
        axes = list(axes)
        meshes = [m for m in meshes if m is not None]
        if not meshes:
            return None
        top = max(float(np.nanmax(np.asarray(m.get_array(), dtype=float)))
                  for m in meshes)
        for mesh in meshes:
            mesh.set_clim(0.0, top)
        cbar = self.fig.colorbar(meshes[-1], ax=axes, pad=0.02, fraction=0.05)
        cbar.ax.tick_params(labelsize=self.tick_font)
        if label:
            cbar.set_label(label, fontsize=self.font)
        self._colourbar_axes.extend(axes)
        return cbar

    def hide(self, ax, message=None):
        """Blank a panel that has no input, optionally saying why."""
        if message:
            ax.text(0.5, 0.5, message, transform=ax.transAxes, ha="center", va="center",
                    fontsize=self.font, color="#666666", wrap=True)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            ax.set_visible(False)

    # -- finishing ----------------------------------------------------------------------

    def finish(self, data=False, lumi_fb=None, text=None, com=13):
        """Draw the CMS label once, then freeze the layout.

        Call this after every panel's axis titles, legends, colour bars and annotations
        exist.  Anything added afterwards gets no space reserved for it and will be drawn
        on top of a neighbouring panel.
        """
        if self._finished:
            note("finish() called twice on the same figure; the second call was ignored")
            return self
        visible = [ax for ax in self.axes[0] if ax.get_visible()]
        if not visible:
            visible = [ax for ax in self.axes.flat if ax.get_visible()]
        left = visible[0]
        right = visible[-1]
        if text is None:
            text = "Preliminary" if data else ""
        hep.cms.label(ax=left, text=text, data=data, rlabel="", loc=0,
                      fontsize=self.font * 1.05)
        pieces = []
        if data and lumi_fb:
            pieces.append(f"{lumi_fb:.1f} fb$^{{-1}}$")
        if DATASET_TAG:
            pieces.append(DATASET_TAG)
        pieces.append(f"({com} TeV)")
        right.text(1.0, 1.005, " ".join(pieces), transform=right.transAxes, ha="right",
                   va="bottom", fontsize=self.font * 1.05)
        self.fig.canvas.draw()
        self.fig.set_layout_engine("none")
        # a colour bar shrinks its parent after the y title was placed, which leaves the
        # title riding at the top of the panel instead of centred on it
        for ax in self._colourbar_axes:
            if not ax.get_ylabel():
                continue
            boxes = [t.get_window_extent() for t in ax.get_yticklabels() if t.get_text()]
            if boxes:
                left = min(box.x0 for box in boxes)
                x = ax.transAxes.inverted().transform((left, 0.0))[0] - 0.03
            else:
                x = -0.12
            ax.yaxis.set_label_coords(x, 0.5)
        self._finished = True
        return self

    def save(self, stem):
        """Write the figure as a vector PDF into the study's (untracked) figures directory."""
        os.makedirs(FIG_DIR, exist_ok=True)
        path = os.path.join(FIG_DIR, f"{stem}.pdf")
        self.fig.savefig(path, bbox_inches="tight")
        print(f"  wrote {path}")
        return path


# --------------------------------------------------------------------------------------
# drawing primitives
# --------------------------------------------------------------------------------------

def draw_shape(ax, histogram, label, colour, canvas, linestyle="-", as_points=False,
               poisson=False, context=""):
    """Draw a 1D histogram as a unit-area shape.

    `poisson=True` puts asymmetric Poisson error bars on the points, which is what an
    observed-count distribution needs; everything else gets the propagated statistical
    error of the raw counts.
    """
    if histogram is None:
        return None
    entries = total_entries(histogram)
    if entries < MIN_SHAPE_ENTRIES:
        where = f" in {context}" if context else ""
        plural = "entry" if entries == 1 else "entries"
        note(f"{_plain(label)}{where}: {entries:.0f} {plural} over "
             f"{len(histogram.axes[0].centers)} bins, below the {MIN_SHAPE_ENTRIES} "
             f"needed to draw a shape at all; not drawn")
        return None
    n_bins = len(histogram.axes[0].centers)
    if entries < MIN_ENTRIES_PER_BIN * n_bins:
        where = f" in {context}" if context else ""
        note(f"{_plain(label)}{where}: {entries:.0f} entries over {n_bins} bins, fewer "
             f"than the {MIN_ENTRIES_PER_BIN} per bin a shape is drawn at; not drawn")
        return None
    centres, density, errors, widths, total = unit_area(histogram)
    # the gate above judges the curve as a whole; a curve that passes it can still hold
    # individual bins with one or two effective entries, which are noise at the width
    # they are drawn at, so those bins are left out rather than drawn
    per_bin = raw_counts(histogram)
    keep = per_bin >= MIN_ENTRIES_PER_BIN
    if not keep.any():
        where = f" in {context}" if context else ""
        note(f"{_plain(label)}{where}: no bin holds {MIN_ENTRIES_PER_BIN} effective "
             f"entries; not drawn")
        return None
    dropped = int((~keep).sum())
    if dropped > 0.2 * len(keep):
        where = f" in {context}" if context else ""
        note(f"{_plain(label)}{where}: {dropped} of {len(keep)} bins hold fewer than "
             f"{MIN_ENTRIES_PER_BIN} effective entries and are not drawn")
    if poisson:
        counts = raw_counts(histogram)
        scale = 1.0 / widths / total if total > 0 else np.zeros_like(widths)
        asym = poisson_errors(counts) * scale
        ax.errorbar(centres[keep], density[keep], yerr=asym[:, keep], fmt="o",
                    color=colour, label=label,
                    markersize=canvas.font * 0.22, elinewidth=1.2, capsize=0)
        return entries
    if as_points:
        ax.errorbar(centres[keep], density[keep], yerr=errors[keep], fmt="s",
                    color=colour, label=label,
                    markersize=canvas.font * 0.20, elinewidth=1.2, capsize=0)
        return entries
    drawn = np.where(keep, density, np.nan)
    ax.stairs(drawn, histogram.axes[0].edges, color=colour, linewidth=2.0,
              linestyle=linestyle, label=label, baseline=None)
    ax.stairs(np.where(keep, density + errors, np.nan), histogram.axes[0].edges,
              baseline=np.where(keep, density - errors, np.nan),
              color=colour, alpha=0.18, fill=True)
    return entries


def draw_map(ax, histogram, canvas, cmap="viridis", column_normalise=True, log_x=False,
             y_index=False, y_tick_every=1):
    """Draw a 2D histogram, by default with each column of the first axis normalised.

    Column normalisation is what makes a correlation visible when the first axis spans
    several orders of magnitude in population.  `y_index=True` puts the companion axis on
    a bin index and labels the ticks with the bin edges, which is the only way to show a
    variable binning whose first edge is zero without dropping that first bin.
    """
    values = np.asarray(histogram.values(), dtype=float)
    if column_normalise:
        totals = values.sum(axis=1, keepdims=True)
        values = np.divide(values, totals, out=np.zeros_like(values), where=totals > 0)
    x_edges = histogram.axes[0].edges
    y_edges = np.asarray(histogram.axes[1].edges, dtype=float)
    if y_index:
        drawn_y = np.arange(len(y_edges), dtype=float)
    else:
        drawn_y = y_edges
    mesh = ax.pcolormesh(x_edges, drawn_y, values.T, cmap=cmap, shading="flat")
    if y_index:
        ticks = np.arange(0, len(y_edges), y_tick_every, dtype=float)
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{y_edges[int(t)]:g}" for t in ticks])
        ax.set_ylim(0, len(y_edges) - 1)
    if log_x:
        canvas.log_axis(ax, "x")
    return mesh


def edge_index(histogram, value, axis=1):
    """Index of the bin edge closest to `value` on the named axis."""
    edges = np.asarray(histogram.axes[axis].edges, dtype=float)
    return int(np.argmin(np.abs(edges - value)))


def slice_groups(axis, boundaries):
    """Group the bins of an axis into slices at the requested boundaries.

    `boundaries` are values, not indices, so a change of binning does not silently move a
    slice.  Returns a list of `(low_index, high_index, label)`.
    """
    edges = np.asarray(axis.edges, dtype=float)
    groups = []
    for low, high in zip(boundaries[:-1], boundaries[1:]):
        i_low = int(np.argmin(np.abs(edges - low)))
        i_high = int(np.argmin(np.abs(edges - high)))
        if i_high <= i_low:
            continue
        groups.append((i_low, i_high, _slice_label(edges[i_low], edges[i_high])))
    return groups


def _slice_label(low, high):
    def fmt(value):
        if value >= 1000:
            return f"{value:.0f}"
        if value >= 10:
            return f"{value:.0f}"
        if value >= 1:
            return f"{value:g}"
        return f"{value:g}"
    if low <= 0:
        return f"< {fmt(high)}"
    return f"{fmt(low)}-{fmt(high)}"


# --------------------------------------------------------------------------------------
# companion variables used in the correlation section
# --------------------------------------------------------------------------------------

#: (histogram, companion axis, axis title, slice boundaries, log x) for every companion the
#: |cos theta*| shape is examined against.  Slice boundaries are values, not bin indices.
COMPANIONS = [
    ("mu_lj_mumu_absCosTheta_vs_iso", "lj_iso", "$\\mu$-LJ isolation",
     [0, 0.1, 0.3, 1.0, 10.0], True),
    ("mu_lj_mumu_absCosTheta_vs_minDxy", "minDxy", "min $|d_{xy}|$ of the pair [cm]",
     [0, 0.02, 0.1, 1.0, 1000.0], True),
    ("mu_lj_mumu_absCosTheta_vs_maxDxy", "maxDxy", "max $|d_{xy}|$ of the pair [cm]",
     [0, 0.02, 0.1, 1.0, 1000.0], True),
    ("mu_lj_mumu_absCosTheta_vs_vxySpread_mu", "vxySpread_mu",
     "LJ transverse-impact-parameter spread [cm]", [0, 0.1, 1.0, 10.0, 1000.0], True),
    ("mu_lj_mumu_absCosTheta_vs_dzSpread_mu", "dzSpread_mu",
     "LJ longitudinal-impact-parameter spread [cm]", [0, 0.1, 1.0, 10.0, 1000.0], True),
    ("mu_lj_mumu_absCosTheta_vs_pt", "lj_pt", "$\\mu$-LJ $p_{T}$ [GeV]",
     [0, 50, 100, 200, 1000], False),
    ("mu_lj_mumu_absCosTheta_vs_mass", "mumu_mass", "$m_{\\mu\\mu}$ [GeV]",
     [0, 1.0, 2.0, 6.0, 200.0], True),
    ("mu_lj_mumu_absCosTheta_vs_nDsa", "nDsa", "displaced-standalone muons in the pair",
     None, False),
]

NDSA_TICKS = {0: "0 (PF-PF)", 1: "1 (PF-DSA)", 2: "2 (DSA-DSA)"}

#: the companions in two groups of four.  Eight panels placed at text width leave each
#: one under two inches across, which is too little for a panel that has to carry a
#: four-entry key; four panels give each one nearly four inches.
COMPANION_GROUPS = {
    "displacement": ["lj_iso", "minDxy", "maxDxy", "vxySpread_mu"],
    "kinematics": ["lj_pt", "mumu_mass", "nDsa", "dzSpread_mu"],
}


def companion_group(name):
    """The entries of `COMPANIONS` belonging to one group, in the group's order."""
    wanted = COMPANION_GROUPS[name]
    by_axis = {entry[1]: entry for entry in COMPANIONS}
    return [by_axis[axis] for axis in wanted]


def _companion_slices(histogram, boundaries):
    """Slices of the companion axis, either at the requested boundaries or one per category."""
    axis = histogram.axes[1]
    if boundaries is None:
        return [(i, i + 1, NDSA_TICKS.get(i, str(i))) for i in range(len(axis))]
    return slice_groups(axis, boundaries)


def _project_slice(histogram, low, high):
    """The |cos theta*| projection of one companion slice, as a 1D histogram."""
    return histogram[:, low:high:sum]


# --------------------------------------------------------------------------------------
# section b -- signal reference
# --------------------------------------------------------------------------------------

def _hi_lo_ratio(histogram, low_edge=0.8, high_edge=0.2):
    """Mean density above `low_edge` divided by the mean density below `high_edge`."""
    centres, density, _, _, _ = unit_area(histogram)
    high = density[centres >= low_edge]
    low = density[centres <= high_edge]
    if len(high) == 0 or len(low) == 0 or low.mean() <= 0:
        return float("nan")
    return float(high.mean() / low.mean())


def fig_gen_vs_reco(m_zd_values=(0.25, 1.2, 5.0), m_xx=None, stem="b1_gen_vs_reco"):
    """Generator-level and reconstructed |cos theta*|, one row per dark-photon mass.

    Columns are the two signal final states.  Every lifetime and mediator mass of a mass
    point is summed, which is legitimate because the angular distribution of the decay
    does not depend on how far the dark photon travelled before decaying.  The masses are
    shown together because the generator-level shape changes with mass and the
    reconstructed one does not follow it everywhere.
    """
    canvas = Canvas(len(m_zd_values), 2, panel_w=6.0, panel_h=4.0)
    drew = False
    panels = []
    grid = [(m, f) for m in m_zd_values for f in ["4Mu", "2Mu2E"]]
    for index, (m_zd, final_state) in enumerate(grid):
        row, col = divmod(index, 2)
        ax = canvas.ax(row, col)
        samples = signal_samples(final_state, m_zd=m_zd, m_xx=m_xx)
        channel = SIGNAL_CHANNEL[final_state]
        if not samples:
            canvas.hide(ax, f"no {final_state} samples\nwith $m_{{Z_d}}$ = {m_zd} GeV")
            note(f"gen/reco panel: no {final_state} sample at m_Zd = {m_zd} GeV")
            continue
        gen = get_hist(samples, "genMu_AFrame_absCosTheta", channel, quiet=True)
        reco = get_hist(samples, "mu_lj_mumu_absCosTheta", channel, quiet=True)
        reference_x = np.linspace(0, 1, 200)
        ax.plot(reference_x, spin_one_density(reference_x), color=COLOURS["reference"],
                linestyle="--", linewidth=2.0, label="$1+\\cos^{2}\\theta^{*}$")
        if gen is not None and total_entries(gen) >= MIN_SHAPE_ENTRIES:
            draw_shape(ax, gen, f"generated (hi/lo = {_hi_lo_ratio(gen):.2f})",
                       COLOURS["gen"], canvas)
            drew = True
        else:
            note(f"gen/reco panel {final_state}: no generator-level histogram available")
        if reco is not None and total_entries(reco) >= MIN_SHAPE_ENTRIES:
            draw_shape(ax, reco, f"reconstructed (hi/lo = {_hi_lo_ratio(reco):.2f})",
                       COLOURS["reco"], canvas, as_points=True)
            drew = True
        canvas.labels(ax, "$|\\cos\\theta^{*}|$" if row == len(m_zd_values) - 1 else None,
                      "normalised pairs / bin width" if col == 0 else None)
        ax.set_xlim(0, 1)
        canvas.headroom(ax, 1.9)
        title = f"{_fs_label(final_state)}, $m_{{Z_d}}$ = {m_zd} GeV\n{_sample_tally(samples)}"
        canvas.legend(ax, loc="upper left", title=title, over_data=True, ncol=1)
        panels.append(ax)
    # every panel carries the same quantity under one axis title, so one scale
    canvas.share_ylimits(panels)
    canvas.finish(data=False)
    if drew:
        canvas.save(stem)
    return canvas


def _sample_tally(samples):
    """How many samples a curve sums, and along which axes of the grid.

    "15 lifetimes summed" was wrong: the fifteen are three mediator masses times five
    lifetimes, and a reader who takes them for one mediator mass will misread the spread.
    """
    if not samples:
        return "no samples"
    parsed = [parse_signal(name) for name in samples]
    by_mass = {}
    for entry in parsed:
        by_mass.setdefault(entry["m_xx"], set()).add(entry["ctau_mm"])
    counts = {len(v) for v in by_mass.values()}
    # the lifetime values are not the same from one mediator mass to the next, so the
    # union of them is not the factor; the count per mediator mass is
    if len(by_mass) > 1 and len(counts) == 1:
        per = counts.pop()
        return (f"{len(samples)} samples ({len(by_mass)} mediator masses "
                f"$\\times$ {per} lifetime{'s' if per > 1 else ''} each)")
    if len(by_mass) > 1:
        return f"{len(samples)} samples over {len(by_mass)} mediator masses"
    per = counts.pop() if counts else 0
    return f"{len(samples)} samples ({per} lifetime{'s' if per > 1 else ''})"


def _fs_label(final_state):
    return "$4\\mu$" if final_state == "4Mu" else "$2\\mu 2e$"


def fig_ridge_map(m_zd=1.2, m_xx=None, stem="b2_ridge_map"):
    """The (|cos theta*|, pT ratio) plane, with the two-body-decay band drawn on top.

    Left: reconstructed pairs.  Right: the same plane at generator level.  Each column of
    the map is normalised to unit sum so the band is visible at every angle.
    """
    canvas = Canvas(1, 2, panel_w=6.6, panel_h=5.4)
    samples = signal_samples("4Mu", m_zd=m_zd, m_xx=m_xx)
    if not samples:
        for col in range(2):
            canvas.hide(canvas.ax(0, col), f"no 4Mu samples\nwith $m_{{Z_d}}$ = {m_zd} GeV")
        note(f"ridge map: no 4Mu sample at m_Zd = {m_zd} GeV")
        canvas.finish(data=False)
        return canvas
    channel = SIGNAL_CHANNEL["4Mu"]
    panels = [("reconstructed", "mu_lj_mumu_ptRatio_vs_absCosTheta"),
              ("generated", "genMu_ptRatio_vs_absCosTheta")]
    speed = beta_star(m_zd)
    meshes, map_axes = [], []
    curve_x = np.linspace(0, 1, 200)
    curve_y = ridge_ratio(curve_x, m_zd)
    for col, (title, name) in enumerate(panels):
        ax = canvas.ax(0, col)
        histogram = get_hist(samples, name, channel, quiet=(col == 1))
        if histogram is None or total_entries(histogram) < MIN_SHAPE_ENTRIES:
            canvas.hide(ax, f"no {title} map available")
            note(f"ridge map: {title} histogram {name} unavailable or empty")
            continue
        mesh = draw_map(ax, histogram, canvas)
        meshes.append(mesh)
        map_axes.append(ax)
        ax.plot(curve_x, curve_y, color="white", linewidth=3.4)
        ax.plot(curve_x, curve_y, color="#d62728", linewidth=2.0,
                label="$(1-\\beta^{*}\\cos\\theta^{*})"
                      "/(1+\\beta^{*}\\cos\\theta^{*})$")
        canvas.labels(ax, "$|\\cos\\theta^{*}|$",
                      "$p_{T}^{sub}/p_{T}^{lead}$" if col == 0 else None)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        canvas.legend(ax, loc="upper right", over_data=True,
                      title=f"{title}, $4\\mu$\n$m_{{Z_d}}$ = {m_zd} GeV, "
                            f"$\\beta^{{*}}$ = {speed:.3f}")
    canvas.shared_colourbar(map_axes, meshes, "fraction of pairs in the column")
    if map_axes:
        canvas.fig.legend(
            handles=[Line2D([], [], color="none")],
            # wrapped: a single line of this length is wider than the figure, and a
            # figure saved wider than it was laid out prints all of its text smaller
            labels=["the panels share one colour scale; the generator-level panel is "
                    "fainter\nonly because the common maximum comes from the "
                    "reconstructed panel"],
            loc="outside lower center", fontsize=canvas.font, frameon=False,
            handlelength=0)
    canvas.finish(data=False)
    canvas.save(stem)
    reconstructed = get_hist(samples, panels[0][1], channel, quiet=True)
    if reconstructed is not None:
        report_ridge_agreement(reconstructed, m_zd)
    return canvas


def report_ridge_agreement(histogram, m_zd, tolerance=0.04):
    """Print how far the mean momentum ratio of each angle column sits from the band.

    The figure draws the analytic curve over the map, which shows that the two agree but
    not by how much; this puts the number on the page next to the claim.
    """
    counts = raw_counts(histogram)
    angle = np.asarray(histogram.axes[0].centers, dtype=float)
    ratio = np.asarray(histogram.axes[1].centers, dtype=float)
    predicted = ridge_ratio(angle, m_zd)
    column = counts.sum(axis=1)
    keep = column >= MIN_DENOMINATOR
    if not keep.any():
        note("ridge agreement: no column holds enough pairs to compare")
        return None
    measured = (counts * ratio).sum(axis=1) / np.maximum(column, 1)
    deviation = np.abs(measured - predicted) / np.maximum(predicted, 1e-12)
    inside = int((deviation[keep] <= tolerance).sum())
    total = int(keep.sum())
    worst = int(np.argmax(np.where(keep, deviation, -1.0)))
    print(f"  reconstructed column means against the band, {total} columns with at least "
          f"{MIN_DENOMINATOR} pairs:")
    print(f"    median deviation {np.median(deviation[keep]):.1%}, "
          f"{inside} of {total} columns within {tolerance:.0%}")
    print(f"    worst column |cos theta*| = {angle[worst]:.2f}: measured "
          f"{measured[worst]:.3f} against {predicted[worst]:.3f}, "
          f"{(measured[worst] - predicted[worst]) / predicted[worst]:+.0%}")
    return deviation


# --------------------------------------------------------------------------------------
# section c -- acceptance
# --------------------------------------------------------------------------------------

TRIGGER_THRESHOLDS = (20.0, 26.0)


def _fraction_below(histogram, threshold, axis_name="subMuPt"):
    """Per |cos theta*| bin: pairs whose companion value is below `threshold`, and the total."""
    axis = histogram.axes[axis_name]
    edges = np.asarray(axis.edges, dtype=float)
    cut = int(np.argmin(np.abs(edges - threshold)))
    counts = raw_counts(histogram)
    return counts[:, :cut].sum(axis=1), counts.sum(axis=1)


def fig_acceptance(m_zd=1.2, stem="c1_acceptance"):
    """Reconstructed |cos theta*| in the two signal final states, and where the
    sub-leading muon momentum sits as a function of the angle."""
    canvas = Canvas(2, 2, panel_w=6.5, panel_h=5.2)
    shapes_ax = canvas.ax(0, 0)
    frac_ax = canvas.ax(0, 1)
    reference_x = np.linspace(0, 1, 200)
    shapes_ax.plot(reference_x, spin_one_density(reference_x), color=COLOURS["reference"],
                   linestyle="--", linewidth=2.0, label="$1+\\cos^{2}\\theta^{*}$")
    any_shape = False
    for final_state in ["4Mu", "2Mu2E"]:
        samples = signal_samples(final_state, m_zd=m_zd)
        if not samples:
            note(f"acceptance: no {final_state} sample at m_Zd = {m_zd} GeV")
            continue
        histogram = get_hist(samples, "mu_lj_mumu_absCosTheta", SIGNAL_CHANNEL[final_state])
        if histogram is None:
            continue
        if draw_shape(shapes_ax, histogram,
                      f"{_fs_label(final_state)} (hi/lo = {_hi_lo_ratio(histogram):.2f})",
                      COLOURS[final_state], canvas) is not None:
            any_shape = True
    canvas.labels(shapes_ax, "$|\\cos\\theta^{*}|$",
                  "normalised pairs / bin width")
    shapes_ax.set_xlim(0, 1)
    canvas.headroom(shapes_ax, 1.75)
    canvas.legend(shapes_ax, loc="upper left", over_data=True,
                  title=f"$m_{{Z_d}}$ = {m_zd} GeV\n"
                        f"{_sample_tally(signal_samples('4Mu', m_zd=m_zd))}")

    drew_fraction = False
    for final_state in ["4Mu", "2Mu2E"]:
        samples = signal_samples(final_state, m_zd=m_zd)
        if not samples:
            continue
        histogram = get_hist(samples, "mu_lj_mumu_absCosTheta_vs_subMuPt",
                             SIGNAL_CHANNEL[final_state])
        if histogram is None:
            continue
        centres = histogram.axes[0].centers
        for threshold, style in zip(TRIGGER_THRESHOLDS, ["-", "--"]):
            passed, total = _fraction_below(histogram, threshold)
            eff, low, high, keep = efficiency_series(passed, total)
            if not keep.any():
                note(f"acceptance: fewer than {MIN_DENOMINATOR} pairs per bin for "
                     f"{final_state}; no points drawn")
                continue
            frac_ax.errorbar(centres[keep], eff[keep], yerr=[low[keep], high[keep]],
                             fmt="o" if style == "-" else "s", linestyle=style,
                             color=COLOURS[final_state], markersize=canvas.font * 0.20,
                             elinewidth=1.2, capsize=0,
                             label=f"{_fs_label(final_state)}, $p_{{T}}^{{sub}} <$ {threshold:.0f} GeV")
            drew_fraction = True
    canvas.labels(frac_ax, "$|\\cos\\theta^{*}|$", "fraction of pairs")
    frac_ax.set_xlim(0, 1)
    frac_ax.set_ylim(0, 1.05)
    if drew_fraction:
        canvas.legend(frac_ax, loc="upper left", over_data=True, ncol=1,
                      title="denominator $\\geq$ 20 pairs per bin")
    else:
        canvas.hide(frac_ax, "no pairs available")

    map_meshes, map_axes = [], []
    for col, final_state in enumerate(["4Mu", "2Mu2E"]):
        ax = canvas.ax(1, col)
        samples = signal_samples(final_state, m_zd=m_zd)
        histogram = get_hist(samples, "mu_lj_mumu_absCosTheta_vs_subMuPt",
                             SIGNAL_CHANNEL[final_state]) if samples else None
        if histogram is None or total_entries(histogram) < MIN_SHAPE_ENTRIES:
            canvas.hide(ax, f"no {final_state} pairs available")
            continue
        mesh = draw_map(ax, histogram, canvas, y_index=True)
        map_meshes.append(mesh)
        map_axes.append(ax)
        canvas.labels(ax, "$|\\cos\\theta^{*}|$",
                      "$p_{T}^{sub}$ bin [GeV]" if col == 0 else None)
        ax.set_xlim(0, 1)
        for threshold in TRIGGER_THRESHOLDS:
            position = edge_index(histogram, threshold)
            ax.axhline(position, color="white", linewidth=2.8)
            ax.axhline(position, color="#d62728", linewidth=1.6, linestyle="--")
        dashed = dict(color="#d62728", ls="--", lw=1.6)
        # the busiest cells of both maps sit along the top; the band below the two lines
        # is empty in 2mu2e and near-empty in 4mu, so the key goes there
        canvas.legend(ax, loc="lower left", over_data=True,
                      handles=[Line2D([], [], **dashed), Line2D([], [], **dashed)],
                      labels=["sub-leading muon $p_{T}$ = 20 GeV",
                              "= 26 GeV (offline requirement)"],
                      title=f"{_fs_label(final_state)}, $m_{{Z_d}}$ = {m_zd} GeV")
    canvas.shared_colourbar(map_axes, map_meshes, "fraction of pairs in the column")
    if map_axes:
        note_text = ("lower row: the vertical axis is the bin index of the sub-leading "
                     "muon momentum, so the bins are equally wide\nand the ticks label "
                     "their edges; the two maps share one colour scale")
        canvas.fig.legend(handles=[Line2D([], [], color="none")], labels=[note_text],
                          loc="outside lower center", fontsize=canvas.font,
                          frameon=False, handlelength=0)
    canvas.finish(data=False)
    if any_shape:
        canvas.save(stem)
    return canvas


def _record_ratio(store, where, name, mass, histogram):
    """Keep the high-to-low ratio of a drawn curve, with its statistical error."""
    counts = raw_counts(histogram)
    edges = np.asarray(histogram.axes[0].edges, dtype=float)
    high = counts[edges[:-1] >= 0.8 - 1e-9].sum()
    low = counts[edges[1:] <= 0.2 + 1e-9].sum()
    value = _hi_lo_ratio(histogram)
    if not np.isfinite(value) or high <= 0 or low <= 0:
        return
    store.append((where, name, mass, value,
                  value * float(np.sqrt(1.0 / high + 1.0 / low))))


def rebin(histogram, factor):
    """Merge `factor` neighbouring bins of a 1D histogram."""
    if histogram is None or factor <= 1:
        return histogram
    return histogram[:: hist.rebin(factor)]


#: short names used where a companion has to fit inside a legend title.
COMPANION_SHORT = {
    "lj_iso": "$\\mu$-LJ isolation",
    "minDxy": "min $|d_{xy}|$ [cm]",
    "maxDxy": "max $|d_{xy}|$ [cm]",
    "vxySpread_mu": "$d_{xy}$ spread [cm]",
    "dzSpread_mu": "$d_{z}$ spread [cm]",
    "lj_pt": "$\\mu$-LJ $p_{T}$ [GeV]",
    "mumu_mass": "$m_{\\mu\\mu}$ [GeV]",
    "nDsa": "DSA muons in pair",
}


# --------------------------------------------------------------------------------------
# section d -- signal, simulated background and data in every region
# --------------------------------------------------------------------------------------

#: (histogram, axis title, bin merging, bin merging in the low-statistics regions).
#: The last factor is applied to the histogram as stored, not on top of the third; the
#: validation regions hold a few hundred effective simulated entries, so at the binning
#: the control regions use, most of their bins would hold one or two.
OBSERVABLES = {
    "absCosTheta": ("mu_lj_mumu_absCosTheta", "$|\\cos\\theta^{*}|$", 1, 5),
    "ptRatio": ("mu_lj_mumu_ptRatio", "$p_{T}^{sub}/p_{T}^{lead}$", 2, 5),
}

#: regions whose simulated statistics need the coarser binning.
LOW_STAT_REGIONS = frozenset(PF_TAGGED_REGIONS)


def fig_region_shapes(observable="absCosTheta", m_zd_values=(0.25, 1.2, 5.0),
                      stem=None, log_y=False):
    """Unit-area signal, simulated-background and observed shapes, region by region.

    Rows are the five regions, columns the dark-photon mass; every lifetime and every
    mediator mass of a mass point is summed.  Simulated backgrounds are combined with
    their cross-section weights before the shape is taken; the observed distribution
    carries Poisson error bars.  Both signal and background are drawn under the same
    selection as the data in the same panel.
    """
    hist_name, x_label, rebin_factor, coarse_factor = OBSERVABLES[observable]
    stem = stem or f"d1_region_shapes_{observable}"
    canvas = Canvas(len(REGIONS), len(m_zd_values), panel_w=4.6, panel_h=2.95)
    backgrounds = background_samples()
    observed = data_samples()
    drew_any = False
    drawn_series = set()
    ratios = []
    for row, region in enumerate(REGIONS):
        row_axes = []
        factor = coarse_factor if region in LOW_STAT_REGIONS else rebin_factor
        for col, m_zd in enumerate(m_zd_values):
            ax = canvas.ax(row, col)
            where = f"{REGION_LABEL[region]}, m_Zd = {m_zd} GeV, {observable}"
            drew_panel = False
            for final_state in ["4Mu", "2Mu2E"]:
                samples = signal_samples(final_state, m_zd=m_zd)
                if not samples:
                    continue
                histogram = rebin(get_hist(samples, hist_name, region, quiet=True),
                                  factor)
                if draw_shape(ax, histogram,
                              f"signal {_fs_label(final_state)}", COLOURS[final_state],
                              canvas, context=where) is not None:
                    drew_panel = True
                    drawn_series.add(f"signal {final_state}")
                    _record_ratio(ratios, REGION_LABEL[region],
                                  f"signal {_plain(_fs_label(final_state))}", m_zd,
                                  histogram)
            if backgrounds:
                histogram = rebin(get_hist(backgrounds, hist_name, region, weighted=True,
                                           quiet=True), factor)
                if draw_shape(ax, histogram, "simulated background",
                              COLOURS["background"], canvas, linestyle=":",
                              context=where) is not None:
                    drew_panel = True
                    drawn_series.add("background")
                    if col == 0:
                        _record_ratio(ratios, REGION_LABEL[region],
                                      "simulated background", None, histogram)
            if observed:
                histogram = rebin(get_hist(observed, hist_name, region, quiet=True),
                                  factor)
                if draw_shape(ax, histogram, "data", COLOURS["data"], canvas,
                              poisson=True, context=where) is not None:
                    drew_panel = True
                    drawn_series.add("data")
                    if col == 0:
                        # the ratio itself goes in the table printed under the figure:
                        # no corner of these panels is free of curves at every row, and
                        # an opaque stamp in the lower right covers the high-angle tail
                        # the ratio is about
                        _record_ratio(ratios, REGION_LABEL[region], "data", None,
                                      histogram)
            drew_any = drew_any or drew_panel
            if col == 0:
                canvas.stamp(ax, REGION_SHORT[region], loc="upper left")
            if row == 0:
                canvas.stamp(ax, f"$m_{{Z_d}}$ = {m_zd} GeV", loc="upper right")
            if row == len(REGIONS) - 1:
                canvas.labels(ax, x_label, None)
            if not drew_panel:
                canvas.hide(ax, "no entries")
                continue
            ax.set_xlim(0, 1)
            if log_y:
                canvas.log_axis(ax, "y")
            else:
                canvas.headroom(ax, 1.55)
                row_axes.append(ax)
        # the observed and the simulated-background curves are the same in every column
        # of a row, so the columns share one vertical scale and stay comparable by eye
        if row_axes:
            top = max(a.get_ylim()[1] for a in row_axes)
            for a in row_axes:
                a.set_ylim(0.0, top)
    handles = [h for key, h in [
        ("signal 4Mu", Line2D([], [], color=COLOURS["4Mu"], lw=2.4,
                              label="signal $4\\mu$")),
        ("signal 2Mu2E", Line2D([], [], color=COLOURS["2Mu2E"], lw=2.4,
                                label="signal $2\\mu 2e$")),
        ("background", Line2D([], [], color=COLOURS["background"], lw=2.4, ls=":",
                              label="simulated background")),
        ("data", Line2D([], [], color=COLOURS["data"], marker="o", ls="none",
                        label="data")),
    ] if key in drawn_series]
    if handles:
        canvas.fig.legend(handles=handles, labels=[h.get_label() for h in handles],
                          loc="outside lower center", ncol=len(handles),
                          fontsize=canvas.font, frameon=False)
    if ratios:
        print("  high-to-low density ratio of every curve drawn "
              "(mean density above 0.8 over mean density below 0.2):")
        width = max(len(r[0]) for r in ratios)
        print(f"    {'region':<{width}} | {'series':<22} | {'m_Zd':>5} | {'hi/lo':>14}")
        for where, name, mass, value, error in ratios:
            mass_text = "-" if mass is None else f"{mass:.2f}"
            print(f"    {where:<{width}} | {name:<22} | {mass_text:>5} | "
                  f"{value:>6.3f} +- {error:<5.3f}")
    canvas.figure_ylabel("normalised pairs / bin width")
    canvas.finish(data=bool(observed), lumi_fb=LUMI_FB)
    if drew_any:
        canvas.save(stem)
    return canvas


# --------------------------------------------------------------------------------------
# section e -- what the angle correlates with
# --------------------------------------------------------------------------------------

def fig_correlation_overlays(samples, region, set_label, stem, is_data=False,
                            companions=None):
    """|cos theta*| shapes in slices of each companion variable.

    One panel per companion; within a panel each curve is the unit-area angular
    distribution of one slice, so a change of shape between curves is a correlation.
    """
    companions = list(companions if companions is not None else COMPANIONS)
    n_cols = 2 if len(companions) <= 4 else 4
    n_rows = -(-len(companions) // n_cols)
    canvas = Canvas(n_rows, n_cols, panel_w=4.1,
                    panel_h=4.6 if n_cols == 4 else 3.6)
    drew_any = False
    legends = []
    populations = []
    for index, (hist_name, axis_name, _label, boundaries, _log) in enumerate(companions):
        ax = canvas.ax(index // n_cols, index % n_cols)
        histogram = get_hist(samples, hist_name, region, quiet=True)
        if histogram is None or total_entries(histogram) < MIN_SHAPE_ENTRIES:
            canvas.hide(ax, f"no entries for\n{COMPANION_SHORT[axis_name]}")
            note(f"correlation overlay ({set_label}): {hist_name} empty in {region}")
            continue
        drawn = 0
        left_top = right_top = 0.0
        for colour, (low, high, label) in zip(_SLICE_COLOURS,
                                              _companion_slices(histogram, boundaries)):
            projection = _project_slice(histogram, low, high)
            entries = total_entries(projection)
            if entries < MIN_SHAPE_ENTRIES:
                continue
            if draw_shape(ax, projection, label, colour, canvas) is None:
                continue
            _, density, errors, _, _ = unit_area(projection)
            top = np.asarray(density + errors, dtype=float)
            half = len(top) // 2
            left_top = max(left_top, float(np.nanmax(top[:half])))
            right_top = max(right_top, float(np.nanmax(top[half:])))
            populations.append((COMPANION_SHORT[axis_name], label, entries))
            drawn += 1
        if drawn == 0:
            canvas.hide(ax, f"every slice of\n{COMPANION_SHORT[axis_name]}\nis below "
                            f"{MIN_SHAPE_ENTRIES} entries")
            continue
        drew_any = True
        ax.set_xlim(0, 1)
        # a modest starting margin: the key is measured afterwards and the axis raised
        # only as far as it actually needs, which at four panels is not far
        canvas.headroom(ax, 1.15)
        canvas.labels(ax, "$|\\cos\\theta^{*}|$",
                      "normalised pairs / bin width" if index % n_cols == 0 else None)
        # one narrow column, so that the key cannot reach past the spine and on to the
        # tick labels; the companion name is its title rather than a second box, which
        # at this panel width would be drawn on top of it.  It goes to whichever top
        # corner the curves leave emptier, and is measured against that side only
        corner = "upper left" if right_top >= left_top else "upper right"
        reference = left_top if corner == "upper left" else right_top
        legends.append((ax, canvas.legend(ax, loc=corner, over_data=True, ncol=1,
                                          title=COMPANION_SHORT[axis_name],
                                          borderpad=0.25, labelspacing=0.22),
                        reference))
    canvas.clear_legends(legends, most=2.3)
    canvas.share_ylimits(legends)
    canvas.fig.legend(handles=[Line2D([], [], color="none", label=set_label)],
                      labels=[set_label], loc="outside lower center",
                      fontsize=canvas.font, frameon=False, handlelength=0)
    canvas.finish(data=is_data, lumi_fb=LUMI_FB)
    if drew_any:
        canvas.save(stem)
        print(f"  pairs behind each slice of {_plain(set_label)}:")
        width = max((len(_plain(c)) for c, _, _ in populations), default=10)
        for companion, label, entries in populations:
            print(f"    {_plain(companion):<{width}} | {_plain(label):>10} | {entries:>12,.0f}")
    return canvas


#: the three data regions the fake shape is measured in, and the colour each gets when
#: they are drawn on one set of axes.
#: Two of the three regions differ by a veto that removes few pairs, so their points land
#: on top of each other; the marker shape is what tells them apart when they do.
PROFILE_REGIONS = [
    ("data_control_region_1muLj", "no cleaning", "#1f2d3d", "o"),
    ("data_control_region_1muLj_cosmic_veto", "cosmic veto", "#2a9d8f", "s"),
    ("data_control_region_1muLj_spread_cosAlpha_mu_veto", "spread + cos alpha veto",
     "#e07a5f", "^"),
]


def fig_correlation_profiles(samples, region, set_label, stem, is_data=False,
                             mean_label=None, regions=None):
    """Mean |cos theta*| as a function of each companion variable.

    Points are drawn only where the companion bin holds at least 20 pairs; the error bar
    is the standard error of the mean in that bin.  `regions` draws one series per region
    on the same axes, which is how a dependence that is flat in one region and not in
    another becomes visible; with a single region the mean over it is drawn as a dashed
    line instead.
    """
    mean_label = mean_label or ("region mean" if is_data else "sample mean")
    series = regions or [(region, None,
                          COLOURS["data"] if is_data else COLOURS["signal"], "o")]
    canvas = Canvas(2, 4, panel_w=4.1, panel_h=3.5)
    drew_any = False
    drawn_axes, seen = [], []
    key_axes, key_entries = None, []
    for index, (hist_name, axis_name, _label, _bounds, log_x) in enumerate(COMPANIONS):
        ax = canvas.ax(index // 4, index % 4)
        drew_panel = False
        for where, region_label, colour, marker in series:
            histogram = get_hist(samples, hist_name, where, quiet=True)
            if histogram is None or total_entries(histogram) < MIN_DENOMINATOR:
                note(f"correlation profile ({set_label}): {hist_name} empty in {where}")
                continue
            centres, mean, error, keep, _counts = profile_mean(histogram)
            if axis_name == "nDsa":
                centres = np.arange(len(histogram.axes[1]), dtype=float)
            if not keep.any():
                continue
            drew_panel = True
            inclusive = float((raw_counts(histogram).sum(axis=1) *
                               histogram.axes[0].centers).sum()
                              / max(total_entries(histogram), 1))
            if len(series) == 1:
                ax.axhline(inclusive, color=COLOURS["reference"], linestyle="--",
                           linewidth=2.0)
            ax.errorbar(centres[keep], mean[keep], yerr=error[keep], fmt=marker,
                        color=colour, markersize=canvas.font * 0.20, elinewidth=1.4,
                        capsize=0, markerfacecolor="none" if marker != "o" else colour)
            seen.append(mean[keep])
            if index == 0:
                if len(series) == 1:
                    key_entries.append((Line2D([], [], color=COLOURS["reference"], ls="--",
                                               lw=2.0),
                                        f"{mean_label} = {inclusive:.3f}"))
                else:
                    key_entries.append((Line2D([], [], color=colour, marker=marker,
                                               ls="none",
                                               markersize=canvas.font * 0.20,
                                               markerfacecolor="none" if marker != "o"
                                               else colour),
                                        f"{region_label}, mean {inclusive:.3f}"))
        if not drew_panel:
            canvas.hide(ax, f"fewer than {MIN_DENOMINATOR} pairs\nin every bin")
            continue
        drew_any = True
        if log_x:
            canvas.log_axis(ax, "x")
        if axis_name == "nDsa":
            ax.set_xticks(list(NDSA_TICKS))
            ax.set_xticklabels([text.split(" ", maxsplit=1)[0]
                                for text in NDSA_TICKS.values()])
            ax.set_xlim(-0.5, len(NDSA_TICKS) - 0.5)
        canvas.labels(ax, COMPANION_SHORT[axis_name],
                      "$\\langle|\\cos\\theta^{*}|\\rangle$" if index % 4 == 0 else None)
        drawn_axes.append(ax)
        if index == 0:
            key_axes = ax
    if drawn_axes:
        # one window for the whole grid: the panels show the same quantity against
        # different companions, and only the left column carries the axis title
        values = np.concatenate(seen)
        low = float(values.min())
        high = float(values.max())
        pad = max(0.08 * (high - low), 0.02)
        low_limit = max(0.0, low - pad)
        # the first panel carries the mean key, so the window keeps its top clear
        # a single-region figure keeps its key inside the first panel and needs the room;
        # a three-region key is wider than a panel and goes under the figure instead
        room = 0.74 if len(series) == 1 else 0.92
        high_limit = min(1.0, low_limit + (high + pad - low_limit) / room)
        for ax in drawn_axes:
            ax.set_ylim(low_limit, high_limit)
    if key_entries and len(series) == 1 and key_axes is not None:
        canvas.legend(key_axes, loc="upper left", over_data=True,
                      handles=[h for h, _ in key_entries],
                      labels=[t for _, t in key_entries])
    elif key_entries:
        # one block under the figure: a second legend at the same anchor is drawn on top
        # of the first rather than beside it
        legend = canvas.fig.legend(handles=[h for h, _ in key_entries],
                                   labels=[t for _, t in key_entries],
                                   loc="outside lower center", ncol=len(key_entries),
                                   fontsize=canvas.font, frameon=False, title=set_label)
        legend.get_title().set_fontsize(canvas.font)
    else:
        canvas.fig.legend(handles=[Line2D([], [], color="none")], labels=[set_label],
                          loc="outside lower center", fontsize=canvas.font,
                          frameon=False, handlelength=0)
    canvas.finish(data=is_data, lumi_fb=LUMI_FB)
    if drew_any:
        canvas.save(stem)
    return canvas


# --------------------------------------------------------------------------------------
# section f -- what a straight cut on the angle would buy
# --------------------------------------------------------------------------------------

def _cumulative_efficiency(histogram, direction="below"):
    """Efficiency of |cos theta*| < x (or > x) at every bin edge, with its interval."""
    counts = raw_counts(histogram)
    total = counts.sum()
    edges = np.asarray(histogram.axes[0].edges, dtype=float)
    if direction == "below":
        passed = np.concatenate([[0.0], np.cumsum(counts)])
    else:
        passed = total - np.concatenate([[0.0], np.cumsum(counts)])
    efficiency, low, high = clopper_pearson(passed, np.full_like(passed, total))
    return edges, efficiency, low, high, total


def cut_scan(m_zd_values=(0.25, 1.2, 5.0), fake_regions=None):
    """Signal efficiency and fake rejection for a straight cut on |cos theta*|.

    Returns a nested dictionary keyed by direction, then by curve name, holding the
    threshold grid, the efficiency and its Clopper-Pearson errors, and the number of pairs
    the efficiency was measured from.
    """
    fake_regions = fake_regions or [UNCUT_CR, CLEANED_CR]
    observed = data_samples()
    result = {}
    for direction in ["below", "above"]:
        curves = {}
        for m_zd in m_zd_values:
            for final_state in ["4Mu", "2Mu2E"]:
                samples = signal_samples(final_state, m_zd=m_zd)
                if not samples:
                    continue
                histogram = get_hist(samples, "mu_lj_mumu_absCosTheta",
                                     SIGNAL_CHANNEL[final_state], quiet=True)
                if histogram is None or total_entries(histogram) < MIN_DENOMINATOR:
                    note(f"cut scan: {final_state} at m_Zd = {m_zd} GeV has too few pairs")
                    continue
                key = (m_zd, f"signal {_fs_label(final_state)}")
                curves[key] = _cumulative_efficiency(histogram, direction)
        for region in fake_regions:
            if not observed:
                continue
            histogram = get_hist(observed, "mu_lj_mumu_absCosTheta", region, quiet=True)
            if histogram is None or total_entries(histogram) < MIN_DENOMINATOR:
                note(f"cut scan: data region {REGION_LABEL[region]} has too few pairs")
                continue
            curves[(None, REGION_LABEL[region])] = _cumulative_efficiency(histogram, direction)
        result[direction] = curves
    return result


DIRECTION_LABEL = {"below": "$|\\cos\\theta^{*}| <$ threshold",
                   "above": "$|\\cos\\theta^{*}| >$ threshold"}


def fig_cut_scan(scan, m_zd_values=(0.25, 1.2, 5.0), stem="f1_cut_scan"):
    """Signal efficiency and fake rejection against the cut threshold.

    Solid curves are the fraction of signal pairs kept; dashed curves are the fraction of
    control-region pairs rejected.  Bands are Clopper-Pearson intervals.
    """
    canvas = Canvas(len(m_zd_values), 2, panel_w=6.4, panel_h=4.6)
    drew_any = False
    drawn_series = {}
    for row, m_zd in enumerate(m_zd_values):
        for col, direction in enumerate(["below", "above"]):
            ax = canvas.ax(row, col)
            curves = scan.get(direction, {})
            drew_panel = False
            for (key_zd, name), (edges, eff, low, high, _n) in curves.items():
                if key_zd is not None and not np.isclose(key_zd, m_zd):
                    continue
                if key_zd is None:
                    is_uncut = _region_of(name) == UNCUT_CR
                    colour = COLOURS["data"] if is_uncut else "#8c564b"
                    ax.plot(edges, 1 - eff, color=colour, linestyle="--", linewidth=2.2)
                    ax.fill_between(edges, 1 - eff - high, 1 - eff + low, color=colour,
                                    alpha=0.18, linewidth=0)
                    drawn_series["uncut" if is_uncut else "cleaned"] = colour
                else:
                    final = "4Mu" if "4" in name else "2Mu2E"
                    colour = COLOURS[final]
                    ax.plot(edges, eff, color=colour, linewidth=2.4)
                    ax.fill_between(edges, eff - low, eff + high, color=colour, alpha=0.18,
                                    linewidth=0)
                    drawn_series[final] = colour
                drew_panel = True
            if not drew_panel:
                canvas.hide(ax, f"no curves for\n$m_{{Z_d}}$ = {m_zd} GeV")
                continue
            drew_any = True
            ax.set_xlim(0, 1)
            # every curve is a fraction, so the band above one is free for the stamp
            ax.set_ylim(0, 1.25)
            canvas.labels(ax, DIRECTION_LABEL[direction] if row == len(m_zd_values) - 1 else None,
                          "fraction of pairs kept (signal)\nor rejected (control region)"
                          if col == 0 else None)
            canvas.stamp(ax, f"$m_{{Z_d}}$ = {m_zd} GeV",
                         loc="upper right" if direction == "below" else "upper left")
    handles = []
    for key, label, style in [
            ("4Mu", "keep signal $4\\mu$", "-"),
            ("2Mu2E", "keep signal $2\\mu 2e$", "-"),
            ("uncut", "reject 1 $\\mu$-LJ CR", "--"),
            ("cleaned", "reject 1 $\\mu$-LJ CR, spread + $\\cos\\alpha$ veto", "--")]:
        if key in drawn_series:
            handles.append(Line2D([], [], color=drawn_series[key], ls=style,
                                  lw=2.4, label=label))
    if handles:
        canvas.fig.legend(handles=handles, labels=[h.get_label() for h in handles],
                          loc="outside lower center", ncol=min(len(handles), 2),
                          fontsize=canvas.font, frameon=False)
    canvas.finish(data=True, lumi_fb=LUMI_FB)
    if drew_any:
        canvas.save(stem)
    return canvas


def _region_of(label):
    for key, value in REGION_LABEL.items():
        if value == label:
            return key
    return ""


#: plain-text forms of the scan direction; a printed table renders no maths markup.
_TABLE_DIRECTION = {"below": "|cos theta*| < threshold",
                    "above": "|cos theta*| > threshold"}


def cut_scan_table(scan, thresholds=(0.2, 0.4, 0.6, 0.8, 0.9), direction="below"):
    """Tabulate the scan at a few thresholds and print it; returns the rows."""
    curves = scan.get(direction, {})
    rows = []
    for (m_zd, name), (edges, eff, low, high, total) in curves.items():
        for threshold in thresholds:
            index = int(np.argmin(np.abs(edges - threshold)))
            rows.append({
                "m_zd": m_zd, "curve": name, "threshold": float(edges[index]),
                "efficiency": float(eff[index]),
                "err_low": float(low[index]), "err_high": float(high[index]),
                "rejection": float(1 - eff[index]), "pairs": float(total),
            })
    if not rows:
        note("cut-scan table: no curves available")
        return rows
    header = (f"{'m_Zd':>6} | {'curve':<38} | {'cut':>5} | {'kept':>18} | {'rejected':>8} | "
              f"{'pairs':>8}")
    print(f"  {_TABLE_DIRECTION[direction]}")
    print("  " + header)
    print("  " + "-" * len(header))
    for row in rows:
        mass = "-" if row["m_zd"] is None else f"{row['m_zd']:.2f}"
        kept = f"{row['efficiency']:.3f} -{row['err_low']:.3f}/+{row['err_high']:.3f}"
        print(f"  {mass:>6} | {_plain(row['curve']):<38} | {row['threshold']:>5.2f} | "
              f"{kept:>18} | {row['rejection']:>8.3f} | {row['pairs']:>8.0f}")
    return rows


def _kept_fraction(histogram, threshold, direction="above"):
    """Fraction of pairs a cut at `threshold` keeps, with its Clopper-Pearson interval."""
    edges, efficiency, low, high, total = _cumulative_efficiency(histogram, direction)
    index = int(np.argmin(np.abs(edges - threshold)))
    return float(efficiency[index]), float(low[index]), float(high[index]), float(total)


def efficiency_by_mediator_mass(m_zd_values=(0.25, 1.2, 5.0),
                                thresholds=(0.2, 0.4, 0.6, 0.8), direction="above"):
    """Print the signal efficiency of the cut separately for each mediator mass.

    The figures sum the three mediator masses of a mass point with raw-count weights, so
    the curve they draw is dominated by whichever mediator mass contributed most pairs.
    This shows what that sum hides.
    """
    rows = []
    print(f"  signal pairs kept by {_TABLE_DIRECTION[direction]}, one row per mediator "
          f"mass:")
    header = (f"    {'channel':<7} | {'m_Zd':>5} | {'m_XX':>6} | {'pairs':>9} | "
              + " | ".join(f"{t:>7.2f}" for t in thresholds))
    print(header)
    print("    " + "-" * (len(header) - 4))
    for final_state in ["4Mu", "2Mu2E"]:
        for m_zd in m_zd_values:
            for m_xx in sorted({parse_signal(n)["m_xx"]
                                for n in signal_samples(final_state, m_zd=m_zd)}):
                samples = signal_samples(final_state, m_zd=m_zd, m_xx=m_xx)
                histogram = get_hist(samples, "mu_lj_mumu_absCosTheta",
                                     SIGNAL_CHANNEL[final_state], quiet=True)
                if histogram is None or total_entries(histogram) < MIN_DENOMINATOR:
                    note(f"per-mass efficiency: {final_state} m_Zd = {m_zd}, "
                         f"m_XX = {m_xx} has too few pairs")
                    continue
                kept = [_kept_fraction(histogram, t, direction) for t in thresholds]
                rows.append({"final_state": final_state, "m_zd": m_zd, "m_xx": m_xx,
                             "pairs": kept[0][3],
                             "kept": [k[0] for k in kept]})
                print(f"    {_plain(_fs_label(final_state)):<7} | {m_zd:>5.2f} | "
                      f"{m_xx:>6.0f} | {kept[0][3]:>9,.0f} | "
                      + " | ".join(f"{k[0]:>7.3f}" for k in kept))
    if not rows:
        note("per-mass efficiency: no signal sample had enough pairs")
    return rows


def validation_region_rejection(thresholds=(0.2, 0.4, 0.6, 0.8), direction="above"):
    """Print what the same cut removes in the two validation regions, in data.

    The cut-scan figure draws only the two single-lepton-jet control regions, so the
    validation-region numbers quoted elsewhere have nowhere on the page to be checked.
    """
    observed = data_samples()
    rows = []
    print(f"  data pairs rejected by {_TABLE_DIRECTION[direction]}, validation regions:")
    header = (f"    {'region':<32} | {'pairs':>8} | "
              + " | ".join(f"{t:>7.2f}" for t in thresholds))
    print(header)
    print("    " + "-" * (len(header) - 4))
    for region in sorted(PF_TAGGED_REGIONS):
        histogram = get_hist(observed, "mu_lj_mumu_absCosTheta", region, quiet=True) \
            if observed else None
        if histogram is None or total_entries(histogram) < MIN_DENOMINATOR:
            note(f"validation-region rejection: {REGION_LABEL[region]} has too few pairs")
            continue
        kept = [_kept_fraction(histogram, t, direction) for t in thresholds]
        rows.append({"region": region, "pairs": kept[0][3],
                     "rejected": [1.0 - k[0] for k in kept]})
        print(f"    {REGION_LABEL[region]:<32} | {kept[0][3]:>8,.0f} | "
              + " | ".join(f"{1.0 - k[0]:>7.3f}" for k in kept))
    if not rows:
        note("validation-region rejection: neither region had enough pairs")
    return rows


def _plain(text):
    """Strip the maths markup from a label so it can go in a plain-text table."""
    return (text.replace("$", "").replace("\\mu", "mu").replace("\\", "")
            .replace("{", "").replace("}", "")
            .replace("4mu", "4mu").replace("2mu 2e", "2mu2e"))


# --------------------------------------------------------------------------------------
# inventory
# --------------------------------------------------------------------------------------

def inventory(regions=None):
    """Print how many pairs each sample group contributes to each region."""
    regions = regions or REGIONS
    groups = [("data", data_samples()), ("background", background_samples())]
    for final_state in ["4Mu", "2Mu2E"]:
        for m_zd in sorted({parse_signal(n)["m_zd"] for n in signal_samples(final_state)}):
            groups.append((f"signal {final_state} m_Zd={m_zd}",
                           signal_samples(final_state, m_zd=m_zd)))
    width = max(len(name) for name, _ in groups) + 2
    header = f"{'group':<{width}}" + "".join(
        f"{r.rsplit('_region_', maxsplit=1)[-1][:22]:>24}" for r in regions)
    print(header)
    print("-" * len(header))
    table = {}
    for name, samples in groups:
        counts = []
        for region in regions:
            histogram = get_hist(samples, "mu_lj_mumu_absCosTheta", region, quiet=True) \
                if samples else None
            counts.append(total_entries(histogram))
        table[name] = counts
        print(f"{name:<{width}}" + "".join(f"{c:>24.0f}" for c in counts))
    return table
