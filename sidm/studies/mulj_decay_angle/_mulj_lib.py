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
from matplotlib.ticker import LogLocator, NullFormatter, NullLocator
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
    "test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto": "2$\\mu$2e VR, inv. displaced\n$\\geq$1 $\\mu$ with >2 pixel hits",
    "test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto": "4$\\mu$ VR, inv. displaced\n$\\geq$1 $\\mu$ with >2 pixel hits",
}

#: single-line names for tables and prints.
REGION_LABEL = {
    "data_control_region_1muLj": "1 mu-LJ CR (no cleaning)",
    "data_control_region_1muLj_cosmic_veto": "1 mu-LJ CR + cosmic veto",
    "data_control_region_1muLj_spread_cosAlpha_mu_veto": "1 mu-LJ CR + spread and cos-alpha veto",
    "test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto": "2mu2e inverted-displacement VR",
    "test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto": "4mu inverted-displacement VR",
}

#: the two inverted-displacement validation regions.  Their requirement asks each muon
#: lepton jet for at least one constituent muon with more than two pixel hits.  A
#: displaced-standalone muon is given zero pixel hits, so that muon is always a PF muon,
#: but the two leading muons the angle is computed from can still include
#: displaced-standalone ones.
VALIDATION_REGIONS = {
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
    """Per-sample MC weight lumi x sigma / sum of generator weights.

    IMPORTANT: with `--unweighted-hist` the stored histograms are raw counts, so this
    factor is what turns them into an expected yield.  The denominator is the
    `scaled_sum_weights` entry of the sample metadata, the sum of genWeight over the events
    stored in the input files, which is the normalisation SidmProcessor itself applies.  The
    figure of about one event in ten kept by the signal skim, which would call for a
    generated count here instead, came from summing a cumulative counter of the Runs tree
    across files, an artefact found on 2026-08-10; no generator filter is
    recorded in these files, and neither signal ntuple location file the campaigns read
    (`signal_4mu_v10.yaml`, `signal_2mu2e_v10.yaml`) sets a `skim_factor`.  Signal samples
    carry 1 fb unless the SIDM utilities are not importable, in which case the fallback
    below takes the value of `cross_sections.yaml` when the sample is listed there.
    Shapes, which are what this study plots, do not depend on the factor at all.
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

    def room_above(self, ax, share=0.45):
        """Leave the top of a panel free for a legend, keeping the curves' own scale.

        `headroom` scales from zero, which suits a spectrum that starts there.  A
        normalised angular distribution varies about a non-zero mean instead, and scaling
        it from zero either buries it in a strip at the top or, once a legend is anchored
        to the bottom, hides its minima underneath.  This keeps the drawn range and gives
        the curves `share` of the panel height.
        """
        ax.autoscale_view()
        low, high = float(ax.dataLim.y0), float(ax.dataLim.y1)
        if not np.isfinite(low) or not np.isfinite(high) or high <= low:
            return
        span = high - low
        ax.set_ylim(low - 0.10 * span, low + span / max(share, 0.05))

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

    def keep_label_position(self, *axes):
        """Exempt a panel from the y-title nudge `finish` applies after a colour bar.

        That nudge re-anchors the title at the middle of the axis.  The plotting style
        anchors axis titles at their top instead, so a title longer than half the panel is
        then drawn from the middle downwards and runs off the foot of the figure.  A panel
        whose tick labels are short numbers gains nothing from the nudge, so it opts out
        and keeps the style's own anchor.
        """
        for ax in axes:
            while ax in self._colourbar_axes:
                self._colourbar_axes.remove(ax)
        return self

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
# section b -- what the stored precision and the arithmetic can resolve
# --------------------------------------------------------------------------------------

TOY_MU_MASS = 0.105658     # GeV, the muon mass the processor assigns to every muon
TOY_SEED = 20261001
TOY_ETA_MAX = 2.4          # both dark photons are kept inside this |eta|
TOY_MUONS = ("lep0", "bar0", "lep1", "bar1")


def round_mantissa(values, bits):
    """values rounded to `bits` mantissa bits, the way NanoAOD stores reduced-precision floats.

    The float32 bit pattern is rounded half up on the dropped mantissa bits, a carry moving
    into the exponent. This reproduces the stored values of the generator record (8 bits for
    GenPart pt, eta, phi and mass) and of the reconstructed muon eta and phi (12 bits).
    """
    values = np.asarray(values, dtype=np.float32)
    if bits >= 23:
        return values
    shift = 23 - bits
    pattern = values.view(np.uint32).astype(np.uint64)
    pattern = (pattern + (1 << (shift - 1))) & ~np.uint64((1 << shift) - 1)
    return pattern.astype(np.uint32).view(np.float32)


def _toy_p4(pt, eta, phi, mass):
    """(E, px, py, pz) rows, in float64, of pt, eta, phi and mass."""
    pt, eta, phi = (np.asarray(v, dtype=np.float64) for v in (pt, eta, phi))
    mass = np.broadcast_to(np.asarray(mass, dtype=np.float64), pt.shape)
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    return np.stack([np.sqrt(px * px + py * py + pz * pz + mass * mass), px, py, pz], axis=-1)


def _toy_boost(p4, beta):
    """p4 boosted by the velocity beta, in the form the vector library uses for .boost()."""
    beta2 = (beta * beta).sum(axis=-1)
    gamma = 1.0 / np.sqrt(1.0 - beta2)
    projection = (beta * p4[..., 1:]).sum(axis=-1)
    space = p4[..., 1:] + ((gamma * gamma / (1.0 + gamma) * projection
                            + gamma * p4[..., 0])[..., None] * beta)
    return np.concatenate([(gamma * (p4[..., 0] + projection))[..., None], space], axis=-1)


def _toy_unit(vectors):
    """Spatial parts of (E, px, py, pz) rows, or plain 3-vectors, scaled to unit length."""
    space = vectors[..., -3:]
    return space / np.sqrt((space * space).sum(axis=-1))[..., None]


def _toy_pt(p4):
    """Transverse momentum of (E, px, py, pz) rows."""
    return np.hypot(p4[..., 1], p4[..., 2])


def toy_abs_cos_theta(daughter, parent):
    """|cos theta*| as the generated histogram of section b computes it, in float64: the
    daughter boosted into the parent's rest frame, its direction measured from the parent's
    lab momentum.

    Both arguments are (E, px, py, pz) rows. cos_theta_in_parent_frame in
    sidm/definitions/hists.py does the same through boost_to_frame and deltaangle; the
    reconstructed histogram of section c is the same construction with the leading muon and
    the two-muon sum as the parent.
    """
    boosted = _toy_boost(daughter, -parent[..., 1:] / parent[..., :1])
    cosine = (_toy_unit(boosted) * _toy_unit(parent)).sum(axis=-1)
    return np.abs(np.clip(cosine, -1.0, 1.0))


def toy_plane_phi(pair0, pair1, ref0, ref1):
    """Phi between the two decay planes as _a_frame_observables computes it, in float64.

    pair0 and pair1 are the two pairs (pair0 the higher-pT one), ref0 and ref1 their reference
    leptons, all (E, px, py, pz) rows. Everything is boosted into the rest frame of the two
    pairs' sum; each plane's normal is unit(P_i x l_i) and Phi is the signed angle between the
    normals about unit(P_0).
    """
    total = pair0 + pair1
    to_a = -total[..., 1:] / total[..., :1]
    p0a, p1a, l0a, l1a = (_toy_boost(v, to_a)[..., 1:] for v in (pair0, pair1, ref0, ref1))
    n0 = _toy_unit(np.cross(p0a, l0a))
    n1 = _toy_unit(np.cross(p1a, l1a))
    return np.arctan2((np.cross(n0, n1) * _toy_unit(p0a)).sum(axis=-1),
                      (n0 * n1).sum(axis=-1))


def toy_helicity_cos(pair, ref, total):
    """Signed helicity angle cos theta_A of one pair as _a_frame_observables computes it, in
    float64: pair and reference lepton boosted into the rest frame of `total` (the sum of the
    two pairs), the lepton then into its pair's rest frame and measured from the pair's
    direction in the `total` frame. All arguments are (E, px, py, pz) rows."""
    to_a = -total[..., 1:] / total[..., :1]
    pair_a, ref_a = _toy_boost(pair, to_a), _toy_boost(ref, to_a)
    ref_star = _toy_boost(ref_a, -pair_a[..., 1:] / pair_a[..., :1])
    return (_toy_unit(ref_star) * _toy_unit(pair_a)).sum(axis=-1)


def _toy_correlation(first, second):
    """Unbinned linear correlation of two angles over the events where both are finite, with
    the error (1 - r^2) / sqrt(N - 1)."""
    first, second = (np.asarray(v, dtype=np.float64) for v in (first, second))
    finite = np.isfinite(first) & np.isfinite(second)
    r = float(np.corrcoef(first[finite], second[finite])[0, 1])
    return r, float((1.0 - r * r) / np.sqrt(finite.sum() - 1.0))


def _toy_decay_angles(beta2, size, rng):
    """cos(theta*) of both decays and the azimuth difference of the two leptons, accepted
    from the joint distribution of a spin-0 decay into two transversely polarised vectors,

        (2 - b + b c0^2)(2 - b + b c1^2) + b^2 (1 - c0^2)(1 - c1^2) cos(2 dphi),  b = beta*^2,

    whose marginals are 1 + alpha cos^2 theta* (alpha = b / (2 - b)) for each decay and
    1 + A cos(2 Phi) (A = [b / (3 - b)]^2) for the plane angle. The density never exceeds 4.
    """
    c0 = rng.uniform(-1.0, 1.0, size)
    c1 = rng.uniform(-1.0, 1.0, size)
    dphi = rng.uniform(0.0, 2.0 * np.pi, size)
    density = ((2.0 - beta2 + beta2 * c0 * c0) * (2.0 - beta2 + beta2 * c1 * c1)
               + beta2 * beta2 * (1.0 - c0 * c0) * (1.0 - c1 * c1) * np.cos(2.0 * dphi))
    keep = rng.uniform(0.0, 4.0, size) < density
    return c0[keep], c1[keep], dphi[keep]


def _toy_generate(m_xx, m_zd, n_events, rng, max_rapidity):
    """Lab four-vectors of n_events XX -> Zd Zd -> 4 mu decays inside the dark-photon |eta|
    acceptance: a dict of (n, 4) float64 arrays zd0/zd1 (the dark photons), lep0/lep1 (the
    reference muon of each) and bar0/bar1 (the other muon)."""
    beta2 = 1.0 - 4.0 * TOY_MU_MASS ** 2 / m_zd ** 2
    q_star = 0.5 * m_zd * np.sqrt(beta2)
    p_zd = np.sqrt(0.25 * m_xx * m_xx - m_zd * m_zd)
    chunks, kept = [], 0
    while kept < n_events:
        c0, c1, dphi = _toy_decay_angles(beta2, 3 * n_events, rng)
        size = len(c0)
        cos_t = rng.uniform(-1.0, 1.0, size)
        sin_t = np.sqrt(1.0 - cos_t * cos_t)
        phi_z = rng.uniform(0.0, 2.0 * np.pi, size)
        phi0 = rng.uniform(0.0, 2.0 * np.pi, size)
        axis = np.stack([sin_t * np.cos(phi_z), sin_t * np.sin(phi_z), cos_t], axis=-1)
        e1 = np.stack([-np.sin(phi_z), np.cos(phi_z), np.zeros(size)], axis=-1)
        e2 = np.cross(axis, e1)
        # the plane azimuths of the two decays share one transverse basis, so their
        # difference is the angle between the planes up to a sign and a half turn, neither
        # of which changes cos(2 Phi)
        out = {}
        for tag, cosine, azimuth, sign in (("0", c0, phi0, 1.0), ("1", c1, phi0 + dphi, -1.0)):
            sine = np.sqrt(1.0 - cosine * cosine)
            direction = (sign * cosine[:, None] * axis
                         + sine[:, None] * (np.cos(azimuth)[:, None] * e1
                                            + np.sin(azimuth)[:, None] * e2))
            zd = np.concatenate([np.full((size, 1), 0.5 * m_xx), sign * p_zd * axis], axis=-1)
            to_lab = zd[:, 1:] / zd[:, :1]
            lep = np.concatenate([np.full((size, 1), 0.5 * m_zd), q_star * direction], axis=-1)
            bar = np.concatenate([np.full((size, 1), 0.5 * m_zd), -q_star * direction], axis=-1)
            out["zd" + tag] = zd
            out["lep" + tag] = _toy_boost(lep, to_lab)
            out["bar" + tag] = _toy_boost(bar, to_lab)
        rapidity = rng.uniform(-max_rapidity, max_rapidity, size)
        longitudinal = np.zeros((size, 3))
        longitudinal[:, 2] = np.tanh(rapidity)
        out = {key: _toy_boost(value, longitudinal) for key, value in out.items()}
        accepted = np.ones(size, dtype=bool)
        for key in ("zd0", "zd1"):
            accepted &= np.abs(np.arcsinh(out[key][:, 3] / _toy_pt(out[key]))) < TOY_ETA_MAX
        chunks.append({key: value[accepted] for key, value in out.items()})
        kept += int(accepted.sum())
    return {key: np.concatenate([chunk[key] for chunk in chunks])[:n_events]
            for key in chunks[0]}


def _toy_fit_alpha(abs_cos):
    """alpha of a 1 + alpha cos^2 density fitted to 25 bins of |cos theta*| over [0, 1], with
    its error: a weighted linear fit of a + b c^2 to the bin contents, alpha = b / a."""
    abs_cos = np.asarray(abs_cos, dtype=np.float64)
    counts, edges = np.histogram(abs_cos[np.isfinite(abs_cos)], bins=25, range=(0.0, 1.0))
    width = edges[1] - edges[0]
    design = np.vstack([np.full(len(counts), width), (edges[1:] ** 3 - edges[:-1] ** 3) / 3]).T
    weight = 1.0 / np.maximum(counts, 1)
    covariance = np.linalg.inv(design.T @ (design * weight[:, None]))
    const, slope = covariance @ (design.T @ (weight * counts))
    jacobian = np.array([-slope / const ** 2, 1.0 / const])
    return float(slope / const), float(np.sqrt(jacobian @ covariance @ jacobian))


def _toy_hi_lo(abs_cos):
    """Mean density of the 25-bin |cos theta*| histogram above 0.8 over that below 0.2, the
    ratio section b quotes."""
    abs_cos = np.asarray(abs_cos, dtype=np.float64)
    counts, _ = np.histogram(abs_cos[np.isfinite(abs_cos)], bins=25, range=(0.0, 1.0))
    centres = np.linspace(0.02, 0.98, 25)
    return float(counts[centres >= 0.8].mean() / counts[centres <= 0.2].mean())


def _toy_plane_coefficient(phi):
    """A of 1 + A cos(2 Phi), as twice the mean of cos(2 Phi), with its error."""
    phi = np.asarray(phi, dtype=np.float64)
    harmonic = np.cos(2.0 * phi[np.isfinite(phi)])
    return 2.0 * float(harmonic.mean()), 2.0 * float(harmonic.std() / np.sqrt(len(harmonic)))


def _toy_summary(abs_cos, phi):
    """alpha, hi/lo, fraction above 0.96 and A of one set of toy angles."""
    abs_cos = np.asarray(abs_cos, dtype=np.float64)
    return {"alpha": _toy_fit_alpha(abs_cos), "hi_lo": _toy_hi_lo(abs_cos),
            "above": float(np.mean(abs_cos[np.isfinite(abs_cos)] > 0.96)),
            "A": _toy_plane_coefficient(phi)}


def _toy_by_pt(first, second, pt_first, pt_second):
    """first and second (equal-length tuples of arrays, events along the first axis)
    reordered event by event so that the one with the higher pT comes first."""
    flip = np.asarray(pt_second > pt_first)

    def choose(keep, swap):
        return np.where(flip.reshape(flip.shape + (1,) * (np.ndim(keep) - 1)), swap, keep)

    lead = tuple(choose(a, b) for a, b in zip(first, second))
    trail = tuple(choose(b, a) for a, b in zip(first, second))
    return lead, trail


def _toy_exact(event):
    """The angles computed exactly (float64) from the toy's four-vectors."""
    abs_cos = np.concatenate([toy_abs_cos_theta(event[mu], event[zd])
                              for mu, zd in (("lep0", "zd0"), ("lep1", "zd1"))])
    (pair0, ref0), (pair1, ref1) = _toy_by_pt((event["zd0"], event["lep0"]),
                                              (event["zd1"], event["lep1"]),
                                              _toy_pt(event["zd0"]), _toy_pt(event["zd1"]))
    summary = _toy_summary(abs_cos, toy_plane_phi(pair0, pair1, ref0, ref1))
    summary["r"] = _toy_correlation(toy_helicity_cos(pair0, ref0, pair0 + pair1),
                                    toy_helicity_cos(pair1, ref1, pair0 + pair1))
    return summary


def _toy_processor(event, m_zd, precision):
    """The angles as the processor computes them, through its own functions and float32
    arithmetic, from the toy's four-vectors stored at one of two precisions.

    precision "generated"     -- a GenPart record as coffea reads it: px, py and pz in
                                 float32 (the momentum is built from those, not from the
                                 8-bit pt, eta and phi branches) and the 8-bit mass branch
    precision "reconstructed" -- muons with pt in float32 and eta and phi at 12 bits, each
                                 pair the sum of its two muons, as for the lepton jets
    """
    import awkward as ak  # pylint: disable=import-outside-toplevel
    from coffea.nanoevents.methods import vector  # pylint: disable=import-outside-toplevel
    from sidm.definitions.hists import (  # pylint: disable=import-outside-toplevel
        _a_frame_observables, boost_to_frame)

    def lorentz(fields):
        return ak.zip(fields, with_name="PtEtaPhiMLorentzVector", behavior=vector.behavior)

    def as_float32(value, size):
        return np.broadcast_to(np.asarray(value, dtype=np.float32), (size,))

    size = len(event["zd0"])
    if precision == "generated":
        mus = {key: lorentz({"px": event[key][:, 1].astype(np.float32),
                             "py": event[key][:, 2].astype(np.float32),
                             "pz": event[key][:, 3].astype(np.float32),
                             "mass": as_float32(TOY_MU_MASS, size)})
               for key in TOY_MUONS}
        pairs = {key: lorentz({"px": event[key][:, 1].astype(np.float32),
                               "py": event[key][:, 2].astype(np.float32),
                               "pz": event[key][:, 3].astype(np.float32),
                               "mass": round_mantissa(np.full(size, m_zd), 8)})
                 for key in ("zd0", "zd1")}
        parts = [(mus["lep0"], pairs["zd0"]), (mus["lep1"], pairs["zd1"])]
        refs = {"zd0": mus["lep0"], "zd1": mus["lep1"]}
    else:
        mus = {key: lorentz({"pt": _toy_pt(event[key]).astype(np.float32),
                             "eta": round_mantissa(np.arcsinh(event[key][:, 3]
                                                              / _toy_pt(event[key])), 12),
                             "phi": round_mantissa(np.arctan2(event[key][:, 2],
                                                              event[key][:, 1]), 12),
                             "mass": as_float32(TOY_MU_MASS, size)})
               for key in TOY_MUONS}
        pairs = {"zd0": mus["lep0"] + mus["bar0"], "zd1": mus["lep1"] + mus["bar1"]}
        parts = []
        for lep, bar, zd in (("lep0", "bar0", "zd0"), ("lep1", "bar1", "zd1")):
            lead = ak.where(mus[lep].pt >= mus[bar].pt, mus[lep], mus[bar])
            parts.append((lead, pairs[zd]))
        refs = {"zd0": mus["lep0"], "zd1": mus["lep1"]}
    abs_cos = []
    for daughter, parent in parts:
        parent_p4 = lorentz({"pt": parent.pt, "eta": parent.eta, "phi": parent.phi,
                             "mass": parent.mass})
        boosted = boost_to_frame(daughter, parent_p4, mass=TOY_MU_MASS)
        abs_cos.append(ak.to_numpy(abs(np.cos(boosted.deltaangle(parent_p4)))))
    columns = {}
    for key in ("zd0", "zd1"):
        for name, source in (("pair", pairs[key]), ("ref", refs[key])):
            columns[key + name] = (ak.to_numpy(source.pt), ak.to_numpy(source.eta),
                                   ak.to_numpy(source.phi), ak.to_numpy(source.mass))
    lead, trail = _toy_by_pt(columns["zd0pair"] + columns["zd0ref"],
                             columns["zd1pair"] + columns["zd1ref"],
                             columns["zd0pair"][0], columns["zd1pair"][0])
    pair0, ref0, pair1, ref1 = lead[:4], lead[4:], trail[:4], trail[4:]
    recs = ak.zip({"pairPt": np.stack([pair0[0], pair1[0]], axis=1),
                   "pairEta": np.stack([pair0[1], pair1[1]], axis=1),
                   "pairPhi": np.stack([pair0[2], pair1[2]], axis=1),
                   "pairMass": np.stack([pair0[3], pair1[3]], axis=1),
                   "refPt": np.stack([ref0[0], ref1[0]], axis=1),
                   "refEta": np.stack([ref0[1], ref1[1]], axis=1),
                   "refPhi": np.stack([ref0[2], ref1[2]], axis=1),
                   "refMass": np.stack([ref0[3], ref1[3]], axis=1)},
                  behavior=vector.behavior)
    frame = _a_frame_observables(recs)
    summary = _toy_summary(np.concatenate(abs_cos), ak.to_numpy(frame["Phi"]))
    summary["r"] = _toy_correlation(ak.to_numpy(frame["cosThetaA0"]),
                                    ak.to_numpy(frame["cosThetaA1"]))
    if precision == "generated":
        # |cos theta*| = |E1 - E2| / (beta* |p1 + p2|) from the float32 energies and momenta,
        # which needs no boost
        beta_star = np.float32(np.sqrt(1.0 - 4.0 * TOY_MU_MASS ** 2 / m_zd ** 2))
        energy = []
        for lep, bar in (("lep0", "bar0"), ("lep1", "bar1")):
            total = mus[lep] + mus[bar]
            energy.append(ak.to_numpy(abs(mus[lep].t - mus[bar].t) / (beta_star * total.p)))
        summary["alpha_energy"] = _toy_fit_alpha(np.concatenate(energy))
    return summary


def _toy_point(m_xx, m_zd, n_events, rng, max_rapidity):
    """Every toy observable of one mass point."""
    beta2 = 1.0 - 4.0 * TOY_MU_MASS ** 2 / m_zd ** 2
    event = _toy_generate(m_xx, m_zd, n_events, rng, max_rapidity)
    opening = (_toy_unit(event["lep0"]) * _toy_unit(event["bar0"])).sum(axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return {"m_xx": m_xx, "m_zd": m_zd, "n": n_events,
                "gamma": float(np.median(event["zd0"][:, 0]) / m_zd),
                "opening": float(np.median(np.arccos(np.clip(opening, -1.0, 1.0)))),
                "alpha_in": beta2 / (2.0 - beta2), "A_in": (beta2 / (3.0 - beta2)) ** 2,
                "exact": _toy_exact(event),
                "generated": _toy_processor(event, m_zd, "generated"),
                "reconstructed": _toy_processor(event, m_zd, "reconstructed")}


def storage_precision_toy(m_xx_values=(100.0, 500.0, 1000.0), m_zd_values=(0.25, 1.2, 5.0),
                          n_events=200_000, max_rapidity=1.0, seed=TOY_SEED):
    """How much of the decay angle and of the plane angle survives the precision with which
    the simulated events are stored and the histograms computed.

    The generated histograms of sections b and j take the generator record's momenta from
    the GenPart px, py and pz branches, which the files keep in full float32 (the pt, eta, phi
    and mass branches keep 8 mantissa bits, and of those only the mass is read), and the
    processor computes in float32. For a dark photon of mass m_Zd and energy E,
    1 - beta^2 = (m_Zd / E)^2 is about 1e-6 at m_Zd = 0.25 GeV and E = 250 GeV, close to the
    float32 resolution of 6e-8, so the boost into its rest frame is not computed reliably
    there. The reconstructed muons keep 12 mantissa bits of eta and phi.

    A toy with no file input. XX decays at rest into two dark photons, back to back along an
    isotropic axis; the event is then given a longitudinal boost of rapidity uniform in
    +-max_rapidity and no transverse motion, and kept when both dark photons lie inside
    |eta| < 2.4. Each dark photon has its nominal mass (no line shape) and decays to two
    muons with no radiation. The angles follow the spin-0 to two transverse vectors
    distribution of _toy_decay_angles: 1 + alpha cos^2 theta* per decay with
    alpha = beta*^2 / (2 - beta*^2), and 1 + A cos(2 Phi) between the planes with
    A = [beta*^2 / (3 - beta*^2)]^2. There is no generator filter, no selection and no
    detector resolution.

    The angles are computed three ways: exactly, in float64 (toy_abs_cos_theta and
    toy_plane_phi, which reproduce the processor's functions to better than 1e-7 on the same
    float64 input); as the generated histograms are, with the processor's own functions on the
    generator record as coffea reads it; and the same at the reconstructed precision, with
    each pair built from its two muons. |cos theta*| follows section b's generated histogram
    (a muon in its dark-photon frame; the histogram fills both muons, which carry the same
    angle, so one per dark photon is used here) and, at the reconstructed precision,
    section c (the leading muon in the two-muon frame); Phi, and the correlation r of the two
    signed helicity angles of section j, follow _a_frame_observables.
    "energy" is alpha from |E1 - E2| / (beta* |p1 + p2|) on the float32 generated muons, a
    reading of the decay angle that needs no boost.
    """
    rng = np.random.default_rng(seed)
    rows = [_toy_point(m_xx, m_zd, n_events, rng, max_rapidity)
            for m_zd in m_zd_values for m_xx in m_xx_values]
    print(f"{n_events} toy events per point, XX rapidity uniform in +-{max_rapidity}, both dark "
          f"photons in |eta| < {TOY_ETA_MAX}; opening: median angle between the two muons [rad]")
    print("\n|cos theta*|: alpha of 1 + alpha cos^2 theta* fitted to 25 bins, hi/lo as in "
          "section b, >0.96: fraction above 0.96")
    print(f"{'':25s} | {'':22s} | {'as the generated histograms':^40s} | "
          f"{'reconstructed precision':^19s}")
    print(f"{'m_XX/m_Zd':>9s} {'gamma':>6s} {'opening':>8s} | {'alpha in':>8s} {'exact':>13s} | "
          f"{'alpha':>13s} {'hi/lo':>5s} {'>0.96':>5s} {'energy':>13s} | {'alpha':>13s} "
          f"{'hi/lo':>5s}")
    print("-" * 113)
    for row in rows:
        gen, reco = row["generated"], row["reconstructed"]
        print(f"{row['m_xx']:5.0f}/{row['m_zd']:<4} {row['gamma']:6.0f} {row['opening']:8.1e} | "
              f"{row['alpha_in']:8.3f} {row['exact']['alpha'][0]:+6.3f}+-"
              f"{row['exact']['alpha'][1]:.3f} | {gen['alpha'][0]:+6.3f}+-{gen['alpha'][1]:.3f} "
              f"{gen['hi_lo']:5.2f} {gen['above']:5.3f} {gen['alpha_energy'][0]:+6.3f}+-"
              f"{gen['alpha_energy'][1]:.3f} | {reco['alpha'][0]:+6.3f}+-"
              f"{reco['alpha'][1]:.3f} {reco['hi_lo']:5.2f}")
    print("\nPhi: A of 1 + A cos(2 Phi), as twice the mean of cos(2 Phi); r: linear correlation "
          "of the two signed helicity angles cos theta_A")
    print(f"{'m_XX/m_Zd':>9s} | {'A in':>6s} {'exact':>15s} | {'generated':>15s} | "
          f"{'reconstructed':>15s} {'/ exact':>7s} | {'r exact':>15s} {'r generated':>15s}")
    print("-" * 112)
    for row in rows:
        exact, gen, reco = (row[key]["A"] for key in ("exact", "generated", "reconstructed"))
        r_exact, r_gen = row["exact"]["r"], row["generated"]["r"]
        print(f"{row['m_xx']:5.0f}/{row['m_zd']:<4} | {row['A_in']:6.3f} "
              f"{exact[0]:+7.4f}+-{exact[1]:.4f} | {gen[0]:+7.4f}+-{gen[1]:.4f} | "
              f"{reco[0]:+7.4f}+-{reco[1]:.4f} {reco[0] / exact[0]:7.2f} | "
              f"{r_exact[0]:+7.4f}+-{r_exact[1]:.4f} {r_gen[0]:+7.4f}+-{r_gen[1]:.4f}")
    return rows


def _overflow_fraction(histogram, axis_name):
    """Raw entries in the overflow bin of one axis, and all raw entries, flow bins included."""
    if histogram is None:
        return 0.0, 0.0
    counts = raw_counts(histogram, flow=True)
    index = [axis.name for axis in histogram.axes].index(axis_name)
    return float(np.take(counts, -1, axis=index).sum()), float(counts.sum())


def overflow_table(m_xx_values=(100.0, 500.0, 1000.0), m_zd_values=(0.25, 1.2, 5.0)):
    """What the overflow bins of the stored |cos theta*| histograms hold, from raw counts.

    Data is read in the uncut single-lepton-jet control region only: the overflow of the
    reconstructed angle and of the momentum ratio, and, from the angle-against-dimuon-mass
    histogram, the dimuon mass of the angle's overflow entries.  Signal is read in its
    inclusive channel, all lifetimes of a mass point summed: the overflow fraction of the
    reconstructed and of the generated angle, and the generated hi/lo of section b from the
    in-range bins.
    """
    data = data_samples()
    cos_data = get_hist(data, "mu_lj_mumu_absCosTheta", UNCUT_CR, quiet=True)
    ratio_data = get_hist(data, "mu_lj_mumu_ptRatio", UNCUT_CR, quiet=True)
    print("overflow of the stored histograms, raw entries, flow bins included")
    over, total = _overflow_fraction(cos_data, "cosTheta")
    ratio_over, ratio_total = _overflow_fraction(ratio_data, "ptRatio")
    print(f"2018 data, uncut single-lepton-jet control region: |cos theta*| {over:,.0f} of "
          f"{total:,.0f} ({over / max(total, 1.0):.1e}); momentum ratio {ratio_over:,.0f} of "
          f"{ratio_total:,.0f}")
    by_mass = get_hist(data, "mu_lj_mumu_absCosTheta_vs_mass", UNCUT_CR, quiet=True)
    if by_mass is not None:
        names = [axis.name for axis in by_mass.axes]
        counts = np.take(raw_counts(by_mass, flow=True), -1, axis=names.index("cosTheta"))
        edges = np.asarray(by_mass.axes["mumu_mass"].edges, dtype=float)
        upper = np.concatenate([[edges[0]], edges[1:], [np.inf]])  # underflow, bins, overflow
        below = counts[upper <= 0.2].sum()
        above = counts[upper > 0.8].sum()
        print(f"  dimuon mass of those |cos theta*| overflow entries: {below:,.0f} below 0.2 GeV "
              f"(negative masses included), {counts.sum() - below - above:,.0f} from 0.2 to "
              f"0.8 GeV, {above:,.0f} above 0.8 GeV")
    print("\nsignal, inclusive channel, all lifetimes summed: overflow fraction of the "
          "reconstructed and generated |cos theta*|, generated hi/lo from the in-range bins")
    print(f"{'final state':>11s} {'m_XX/m_Zd':>10s} | {'reconstructed':>13s} {'entries':>10s} | "
          f"{'generated':>10s} {'entries':>10s} | {'gen hi/lo':>9s}")
    print("-" * 82)
    rows = []
    for final_state in ("4Mu", "2Mu2E"):
        for m_zd in m_zd_values:
            for m_xx in m_xx_values:
                samples = signal_samples(final_state, m_zd=m_zd, m_xx=m_xx)
                if not samples:
                    continue
                channel = SIGNAL_CHANNEL[final_state]
                reco = _overflow_fraction(get_hist(samples, "mu_lj_mumu_absCosTheta", channel,
                                                   quiet=True), "cosTheta")
                gen_hist = get_hist(samples, "genMu_AFrame_absCosTheta", channel, quiet=True)
                gen = _overflow_fraction(gen_hist, "cosTheta")
                hi_lo = _hi_lo_ratio(gen_hist) if gen_hist is not None else float("nan")
                rows.append({"final_state": final_state, "m_xx": m_xx, "m_zd": m_zd,
                             "reco": reco, "gen": gen, "gen_hi_lo": hi_lo})
                print(f"{channel:>11s} {m_xx:5.0f}/{m_zd:<4} | "
                      f"{reco[0] / max(reco[1], 1.0):13.4f} {reco[1]:10,.0f} | "
                      f"{gen[0] / max(gen[1], 1.0):10.4f} {gen[1]:10,.0f} | {hi_lo:9.2f}")
    return rows


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
LOW_STAT_REGIONS = frozenset(VALIDATION_REGIONS)


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
    for region in sorted(VALIDATION_REGIONS):
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


# ======================================================================================
# the channel definitions the inputs were produced with
# ======================================================================================

#: the channel definitions every campaign of this study ran with, fully expanded and under
#: the names stored in the histograms' channel axis.
AS_RUN_SELECTIONS = os.path.join(_CONFIG_DIR, "selections_mulj_decay_angle_as_run.yaml")

#: the analysis-wide channel definitions of the checkout this module is read from.
REFERENCE_SELECTIONS = os.path.join(_CONFIG_DIR, "selections.yaml")

#: how the tables below group the channels.
SELECTION_GROUPS = [
    ("signal-like (simulation only)",
     ["4mu", "4mu_sr", "4mu_noMuPtCut", "2mu2e", "2mu2e_sr", "2mu2e_noMuPtCut"]),
    ("control",
     ["data_control_region_1muLj", "data_control_region_1muLj_cosmic_veto",
      "data_control_region_1muLj_spread_cosAlpha_mu_veto",
      "data_control_region_1egmLj_cosmic_veto",
      "cr_2lj_invDisplaced_spread_cosAlpha_mu_veto"]),
    ("validation",
     ["test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto",
      "test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto"]),
]

#: generator-level collections; their cuts select truth particles, not reconstructed
#: objects, and are left out of the tables.
_GEN_OBJECTS = ("genMus", "genEs", "genAs")

_OBJECT_KEYS = ("obj_cuts", "preLj_obj_cuts")
_LJ_KEYS = ("postLj_obj_cuts",)


def _flatten_cuts(cuts):
    """A cut list flattened the way the processor flattens it: dictionaries contribute
    their values, nested lists are expanded in order."""
    flat = []
    for item in (cuts.values() if isinstance(cuts, dict) else (cuts or [])):
        if isinstance(item, (dict, list)):
            flat.extend(_flatten_cuts(item))
        else:
            flat.append(item)
    return flat


def load_channel_definitions(path):
    """Every channel of a selections file as `{name: {"obj": {...}, "lj": {...}, "evt": [...]}}`.

    YAML anchors and merge keys are resolved by the loader and the cut lists are
    flattened, so two channels compare equal exactly when the processor would apply the
    same cuts in the same order.
    """
    with open(path, encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    channels = {}
    for name, definition in raw.items():
        if not isinstance(definition, dict):
            continue
        entry = {"obj": {}, "lj": {}, "evt": _flatten_cuts(definition.get("evt_cuts") or [])}
        for slot, keys in (("obj", _OBJECT_KEYS), ("lj", _LJ_KEYS)):
            for key in keys:
                for obj, cuts in (definition.get(key) or {}).items():
                    entry[slot].setdefault(obj, []).extend(_flatten_cuts(cuts))
        channels[name] = entry
    return channels


def _md(text):
    """Escape the characters a markdown table cell cannot hold as they are."""
    return str(text).replace("|", "\\|")


def _cell(per_object):
    """One table cell from `{object: [cuts]}`; objects without cuts are left out."""
    parts = [f"{obj}: " + ", ".join(_md(c) for c in cuts)
             for obj, cuts in per_object.items() if cuts]
    return "<br>".join(parts) if parts else "none beyond the baseline"


def _grouped_channels(channels):
    """The channels of `channels` in table order, with any not in a group appended."""
    grouped, seen = [], set()
    for title, names in SELECTION_GROUPS:
        present = [n for n in names if n in channels]
        seen.update(present)
        grouped.append((title, present))
    rest = [n for n in channels if n not in seen]
    if rest:
        grouped.append(("other", rest))
    return grouped


def selections_table(path=AS_RUN_SELECTIONS):
    """A markdown table of the channel definitions in `path`.

    The first row is the baseline every channel shares; each channel row lists only the
    cuts it adds to it, in the order they are applied.  Generator-level selections are left
    out.
    """
    channels = load_channel_definitions(path)
    names = list(channels)

    def common(slot, obj):
        lists = [channels[n][slot].get(obj, []) for n in names]
        return [c for c in lists[0] if all(c in other for other in lists[1:])]

    base = {}
    for slot in ("obj", "lj"):
        objects = []
        for n in names:
            objects += [o for o in channels[n][slot] if o not in objects and o not in _GEN_OBJECTS]
        base[slot] = {o: common(slot, o) for o in objects}
    base_evt = [c for c in channels[names[0]]["evt"]
                if all(c in channels[n]["evt"] for n in names)]

    lines = ["| channel | object cuts | lepton-jet cuts | event cuts |",
             "| --- | --- | --- | --- |",
             f"| **baseline, every channel** | {_cell(base['obj'])} | {_cell(base['lj'])} | "
             f"{', '.join(_md(c) for c in base_evt)} |"]
    for title, members in _grouped_channels(channels):
        if not members:
            continue
        lines.append(f"| ***{title}*** | | | |")
        for name in members:
            entry = channels[name]
            extra = {}
            for slot in ("obj", "lj"):
                extra[slot] = {o: [c for c in cuts if c not in base[slot].get(o, [])]
                               for o, cuts in entry[slot].items() if o not in _GEN_OBJECTS}
            evt = [c for c in entry["evt"] if c not in base_evt]
            lines.append(f"| `{name}` | {_cell(extra['obj'])} | {_cell(extra['lj'])} | "
                         f"{', '.join(_md(c) for c in evt) or 'none beyond the baseline'} |")
    return "\n".join(lines)


#: the parts of the cosmic-veto baseline the analysis-wide definitions adopted after these
#: channels were produced, as (label, test on a loaded channel).
_BASELINE_FEATURES = [
    ("eta-phi veto",
     lambda ch: any("eta-phi" in c for cuts in ch["obj"].values() for c in cuts)),
    (">=2 muons (pf or dsa) pT > 26 GeV",
     lambda ch: ">=2 muons (pf or dsa) pT > 26 GeV" in ch["evt"]),
    ("all cosAlpha(mu, mu) > -0.99",
     lambda ch: "all cosAlpha(mu, mu) > -0.99" in ch["evt"]),
]


def _spread_cuts(channel):
    """The muon lepton-jet spread variables a channel cuts on, as a short string."""
    cuts = [c.split()[0] for c in channel["lj"].get("mu_ljs", []) if "Spread" in c]
    return ", ".join(cuts) if cuts else "none"


def selections_vs_reference(path=AS_RUN_SELECTIONS, reference=REFERENCE_SELECTIONS):
    """A markdown table comparing each channel of `path` with the channel of the same name
    in `reference`: whether it exists there, whether it is defined identically, and which
    of the later baseline cuts each version contains, as "as run / reference" in each cell.
    """
    ours = load_channel_definitions(path)
    theirs = load_channel_definitions(reference)
    head = (["channel", "in reference", "identical"]
            + [_md(label) for label, _ in _BASELINE_FEATURES] + ["muon-LJ spread cuts"])
    lines = ["| " + " | ".join(head) + " |", "|" + " --- |" * len(head)]

    def yes(flag):
        return "yes" if flag else "no"

    for title, members in _grouped_channels(ours):
        if not members:
            continue
        lines.append(f"| ***{title}*** |" + " |" * (len(head) - 1))
        for name in members:
            mine, other = ours[name], theirs.get(name)
            row = [f"`{name}`", yes(other is not None)]
            if other is None:
                row.append("--")
                row += [f"{yes(test(mine))} / --" for _, test in _BASELINE_FEATURES]
                row.append(f"{_spread_cuts(mine)} / --")
            else:
                row.append(yes(mine == other))
                row += [f"{yes(test(mine))} / {yes(test(other))}"
                        for _, test in _BASELINE_FEATURES]
                row.append(f"{_spread_cuts(mine)} / {_spread_cuts(other)}")
            lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


# ======================================================================================
# displaced-standalone muons in the data pairs of each region
# ======================================================================================

#: the regions the composition table covers, in order, with the campaign whose merged
#: outputs each is read from.  The two-lepton-jet control region was written by the
#: decay-plane campaign only.
COMPOSITION_REGIONS = [
    (UNCUT_CR, "decay angle"),
    (CLEANED_CR, "decay angle"),
    ("cr_2lj_invDisplaced_spread_cosAlpha_mu_veto", "decay planes"),
    ("test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto", "decay angle"),
    ("test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto", "decay angle"),
]

_COMPOSITION_LABEL = {
    UNCUT_CR: "1 mu-LJ CR, no cleaning",
    CLEANED_CR: "1 mu-LJ CR, spread + cos-alpha",
    "cr_2lj_invDisplaced_spread_cosAlpha_mu_veto": "2-LJ inv.-displaced CR",
    "test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto": "2mu2e inv.-displaced VR",
    "test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto": "4mu inv.-displaced VR",
}

_PAIR_CLASS = {0: "PF-PF", 1: "PF-DSA", 2: "DSA-DSA"}

_COMPOSITION_HIST = "mu_lj_mumu_absCosTheta_vs_nDsa"


def _composition_counts(region, campaign):
    """Data counts of the decay angle against the number of displaced-standalone muons in
    the pair, and the angle's bin edges, in one control or validation region.

    In-range bins only.  Returns (None, None) when the region holds no data pairs; a
    region that is neither a control nor a validation region raises before anything is
    read.
    """
    if region not in set(REGIONS) | set(PLANES_DATA_CHANNELS):
        raise ValueError(f"{region!r} is not a control or validation region; "
                         "collision data is not read there")
    if campaign == "decay planes":
        counts = planes_counts(planes_data_samples(), _COMPOSITION_HIST, region)
        edges = planes_edges(_COMPOSITION_HIST, "cosTheta")
    else:
        histogram = get_hist(data_samples(), _COMPOSITION_HIST, region, quiet=True)
        counts = None if histogram is None else raw_counts(histogram)
        edges = None if histogram is None else histogram.axes["cosTheta"].edges
    if counts is None or edges is None or float(np.sum(counts)) == 0.0:
        return None, None
    return np.asarray(counts, dtype=float), np.asarray(edges, dtype=float)


def dsa_composition_table(threshold=0.4):
    """Data pairs of each control and validation region split by the number of
    displaced-standalone muons in the pair, and the fraction of each class with
    |cos theta*| above `threshold`, with its binomial uncertainty sqrt(p (1 - p) / n).

    Collision data only, read from the stored histogram of the angle against that number;
    in-range bins only, as everywhere else in this study.  Every region read here is a
    control or validation region; anything else raises.
    """
    rows = []
    print("  data muon pairs by number of displaced-standalone (DSA) muons in the pair; "
          f"f = fraction with |cos theta*| > {threshold}")
    header = (f"    {'region':<32} | {'pair':<7} | {'pairs':>9} | {'share':>6} | "
              f"{'f':>6} | {'+-':>6}")
    print(header)
    print("    " + "-" * (len(header) - 4))
    for region, campaign in COMPOSITION_REGIONS:
        counts, edges = _composition_counts(region, campaign)
        if counts is None:
            note(f"DSA composition: no data pairs in {_COMPOSITION_LABEL[region]}")
            continue
        above = edges[:-1] >= threshold - 1e-9
        total = float(counts.sum())
        for n_dsa, label in _PAIR_CLASS.items():
            n = float(counts[:, n_dsa].sum())
            row = {"region": region, "class": label, "pairs": n, "share": n / total}
            first = _COMPOSITION_LABEL[region] if n_dsa == 0 else ""
            if n > 0:
                p = float(counts[above, n_dsa].sum()) / n
                row.update(fraction=p, error=float(np.sqrt(p * (1.0 - p) / n)))
                flag = "" if n >= MIN_DENOMINATOR else "   (fewer than 20 pairs)"
                print(f"    {first:<32} | {label:<7} | {n:>9,.0f} | {n / total:>6.3f} | "
                      f"{p:>6.3f} | {row['error']:>6.3f}{flag}")
            else:
                print(f"    {first:<32} | {label:<7} | {n:>9,.0f} | {0.0:>6.3f} | "
                      f"{'--':>6} | {'--':>6}")
            rows.append(row)
    return rows


#: the regions the composition bound is computed in: the two whose particle-flow and
#: mixed pairs have unlike shapes in the composition table.
COMPOSITION_BOUND_REGIONS = [
    ("test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto", "decay angle"),
    ("cr_2lj_invDisplaced_spread_cosAlpha_mu_veto", "decay planes"),
]


def _shared_type_excess(above, total):
    """sum_k p_k f_k^2 / (sum_k p_k f_k)^2 - 1 from per-type counts above a threshold and
    in total, along the last axis; types with no pairs drop out."""
    above = np.asarray(above, dtype=float)
    total = np.asarray(total, dtype=float)
    squares = np.divide(above ** 2, total, out=np.zeros_like(above), where=total > 0)
    return total.sum(axis=-1) * squares.sum(axis=-1) / above.sum(axis=-1) ** 2 - 1.0


def composition_excess_bound(thresholds=(0.4, 0.8), observed=None, replicas=2000,
                             seed=20261001):
    """The largest excess over the squared per-jet fraction that the mixture of pair types
    alone can make.

    Let p_k be the share of muon pairs of type k (particle-flow, mixed, two
    displaced-standalone) and f_k the fraction of them above a threshold.  If the two muon
    lepton jets of an event always carried the same pair type and were independent within
    it, the fraction of events with both above the threshold would be sum_k p_k f_k^2,
    against the squared per-jet fraction (sum_k p_k f_k)^2.  The ratio minus one is the
    excess this composition makes with no angular correlation at fixed type; jets that
    share their type less often make less, so for these shares it is an upper bound.

    p_k and f_k are per-jet quantities, from the stored data histogram of the angle against
    the number of displaced-standalone muons in the pair over every qualifying muon lepton
    jet of the region, in-range bins; each threshold is moved to the nearest bin edge, which
    is printed.  The uncertainty is the standard deviation over replicas that resample
    events, not jets: the histogram does not record which two jets share an event, so each
    replica draws half as many events as there are jets from the type-by-side table and gives
    both jets of an event the same cell.  That is the limit of complete correlation within an
    event; independent jets would give a spread smaller by a factor of about sqrt(2).
    `observed`, the rows `pair_factorisation_table` returns, puts the
    excess measured in the four-muon validation region beside the bound; no other region
    has a measured one.
    """
    measured = {round(row["threshold"], 6): row for row in (observed or [])}
    rows = []
    print("excess over the squared per-jet fraction that the mixture of pair types alone makes")
    print("if both muon LJs of an event always share a pair type and are independent within it")
    print("2018 data; shares p and fractions f per pair type PF-PF / PF-DSA / DSA-DSA; percent")
    header = (f"    {'region':<24} | {'cut at':>6} | {'p':^17} | {'f':^17} | "
              f"{'bound':>12} | {'observed incl.':>14} | {'observed ord.':>14}")
    print(header)
    print("    " + "-" * (len(header) - 4))
    for region, campaign in COMPOSITION_BOUND_REGIONS:
        counts, edges = _composition_counts(region, campaign)
        if counts is None:
            note(f"composition bound: no data pairs in {_COMPOSITION_LABEL[region]}")
            continue
        total = counts.sum(axis=0)
        for threshold in thresholds:
            index = int(np.argmin(np.abs(edges - threshold)))
            above = counts[index:, :].sum(axis=0)
            bound = float(_shared_type_excess(above, total))
            cells = np.concatenate([above, total - above])
            events = int(round(cells.sum() / 2.0))
            draws = 2 * np.random.default_rng(seed).multinomial(
                events, cells / cells.sum(), size=replicas)
            sigma = float(np.std(_shared_type_excess(draws[:, :3], draws[:, :3] + draws[:, 3:]),
                                 ddof=1))
            shares = total / total.sum()
            fractions = np.divide(above, total, out=np.full_like(above, np.nan),
                                  where=total > 0)
            row = {"region": region, "threshold": threshold, "edge": float(edges[index]),
                   "shares": shares, "fractions": fractions, "bound": bound, "sigma": sigma}
            seen = measured.get(round(threshold, 6)) if region == PAIR_VR else None
            if seen is not None:
                row.update(observed_inclusive=seen["ratio_inclusive"] - 1.0,
                           sigma_inclusive=seen["sigma_inclusive"],
                           observed_ordered=seen["ratio_ordered"] - 1.0,
                           sigma_ordered=seen["sigma_ordered"])
                incl = (f"{100 * row['observed_inclusive']:5.1f} +- "
                        f"{100 * row['sigma_inclusive']:4.1f}")
                ordered = (f"{100 * row['observed_ordered']:5.1f} +- "
                           f"{100 * row['sigma_ordered']:4.1f}")
            else:
                incl = ordered = "not measured"
            label = _COMPOSITION_LABEL[region] if threshold == thresholds[0] else ""
            print(f"    {label:<24} | {row['edge']:>6.2f} | "
                  + " ".join(f"{value:5.3f}" for value in shares) + " | "
                  + " ".join(f"{value:5.3f}" for value in fractions)
                  + f" | {100 * bound:5.1f} +- {100 * sigma:4.1f} | {incl:>14} | {ordered:>14}")
            rows.append(row)
    return rows


def pair_type_separation(threshold=0.4, m_xx=500.0, m_zd=1.2,
                         region="cr_2lj_invDisplaced_spread_cosAlpha_mu_veto"):
    """Signal efficiency over root fake efficiency of a threshold on the decay angle,
    against all the data pairs of a region and against each pair type alone.

    The signal efficiency is that of the four-muon point (m_xx, m_zd) in its own channel,
    read as `planes_mu_separation_table` reads it, so the all-pairs row reproduces that
    table's entry; the fake fractions come from the data histogram of the angle against
    the number of displaced-standalone muons in the pair.
    """
    samples = planes_points("4Mu").get((m_xx, m_zd))
    signal = None if not samples else planes_counts(samples, PLANES_MU_HIST, "4mu")
    counts, edges = _composition_counts(region, dict(COMPOSITION_REGIONS)[region])
    if signal is None or counts is None:
        note(f"pair-type separation: no input for 4Mu {m_xx:.0f}/{m_zd} or "
             f"{_COMPOSITION_LABEL[region]}")
        return []
    eff_s = kept_fraction(signal.sum(axis=1), threshold)
    index = int(np.argmin(np.abs(edges - threshold)))
    print(f"signal 4mu, m_XX = {m_xx:.0f} GeV, m_Zd = {m_zd} GeV: "
          f"eff(|cos theta*| > {eff_s['edge']:.2f}) = {eff_s['fraction']:.4f}")
    print(f"against the data pairs of the {_COMPOSITION_LABEL[region]}, "
          f"cut at {edges[index]:.2f}")
    print(f"    {'fake pairs':<10} | {'pairs':>9} | {'eff_b':>6} | {'eff_s/sqrt(eff_b)':>17}")
    rows = []
    for label, column in [("all", slice(None))] + [(name, n_dsa)
                                                   for n_dsa, name in _PAIR_CLASS.items()]:
        selected = counts[:, column]
        n = float(selected.sum())
        if n <= 0:
            continue
        eff_b = float(selected[index:].sum()) / n
        gain = eff_s["fraction"] / np.sqrt(eff_b) if eff_b > 0 else float("nan")
        flag = "" if n >= MIN_DENOMINATOR else "   (fewer than 20 pairs)"
        print(f"    {label:<10} | {n:>9,.0f} | {eff_b:6.3f} | {gain:17.3f}{flag}")
        rows.append({"pairs": label, "n": n, "eff_s": eff_s["fraction"], "eff_b": eff_b,
                     "gain": gain})
    return rows


# ======================================================================================
# the two muon lepton jets of a four-muon event
# ======================================================================================
#
# The pair campaign fills two histograms the per-lepton-jet ones above cannot:
#
#   mulj_pair_absCosTheta0_vs_absCosTheta1
#       one entry per event holding two qualifying muon lepton jets, carrying the decay
#       angle of each, ordered by the transverse momentum of the parent lepton jet;
#   mulj_pair_absCosTheta_vs_mjj
#       one entry per qualifying muon lepton jet, carrying its decay angle and the
#       muon-lepton-jet-pair invariant mass of its event, on a binning fine enough to cut
#       the four-muon mass axis into windows.
#
# They live in their own merged directory, so the loader below is separate from the one
# above and leaves its cache alone.

#: merged outputs of the pair campaign.
PAIR_INPUT = "/uscms_data/d3/murtazas/mulj_pair_merged"

#: the only region in which collision data was written for the pair histograms.
PAIR_VR = "test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto"

PAIR_PLANE = "mulj_pair_absCosTheta0_vs_absCosTheta1"
PAIR_MJJ = "mulj_pair_absCosTheta_vs_mjj"

#: the search uses bound-state masses at or above this, so every number in this section is
#: quoted for the points at or above it.
PAIR_MIN_M_XX = 200.0

#: the muon-lepton-jet-pair mass range each bound-state mass populates.  The four-muon
#: validation region cuts on that mass at 150 GeV only, so these are the bins of the
#: histogram axis that hold the signal, not an existing analysis cut.
PAIR_WINDOW = {500.0: (380.0, 600.0), 1000.0: (600.0, 1.0e6)}

#: fake sample used where a signal's own window holds too few data pairs to measure a
#: shape in.
PAIR_POOLED_WINDOW = (300.0, 1.0e6)

_PAIR_LOADED: dict = {}
#: stored sum of generator weights per pair-campaign sample, from its own metadata.
_PAIR_SUM_WEIGHTS: dict = {}


def configure_pair(input_dir=None):
    """Point the pair loader at a directory of merged pair-campaign outputs."""
    global PAIR_INPUT
    if input_dir is not None and input_dir.rstrip("/") != PAIR_INPUT:
        PAIR_INPUT = input_dir.rstrip("/")
        _PAIR_LOADED.clear()
        _PAIR_SUM_WEIGHTS.clear()
    return PAIR_INPUT


def _pair_listing():
    """Sample names available in `PAIR_INPUT`, sorted."""
    if PAIR_INPUT.startswith("root://"):
        files = _remote_listing(PAIR_INPUT)
    else:
        files = [os.path.basename(p) for p in glob.glob(os.path.join(PAIR_INPUT, "*.coffea"))]
    return sorted(f[: -len(".coffea")] for f in files)


def _pair_load(name):
    """The histogram dictionary of one pair-campaign sample, or None with a note."""
    if name in _PAIR_LOADED:
        return _PAIR_LOADED[name]
    try:
        if PAIR_INPUT.startswith("root://"):
            os.makedirs(CACHE_DIR, exist_ok=True)
            path = os.path.join(CACHE_DIR, f"pair_{name}.coffea")
            if not os.path.exists(path):
                subprocess.run(["xrdcp", "-f", f"{PAIR_INPUT}/{name}.coffea", path], check=True)
        else:
            path = os.path.join(PAIR_INPUT, f"{name}.coffea")
        blob = coffea.util.load(path)
    except (FileNotFoundError, OSError, subprocess.CalledProcessError) as exc:
        note(f"pair sample {name} could not be read ({type(exc).__name__}); skipping it")
        _PAIR_LOADED[name] = None
        return None
    out = blob["out"] if "out" in blob else blob
    if name in out:
        payload = out[name]
    elif len(out) == 1:
        payload = next(iter(out.values()))
    else:
        note(f"pair file for {name} holds {len(out)} samples and none matches its name; skipping")
        _PAIR_LOADED[name] = None
        return None
    meta = payload.get("metadata", {})
    _PAIR_SUM_WEIGHTS[name] = float(
        meta.get("scaled_sum_weights", meta.get("n_evts", 0.0)) or 0.0)
    _PAIR_LOADED[name] = payload["hists"]
    return _PAIR_LOADED[name]


def pair_data_samples():
    """Collision-data samples of the pair campaign, all eras."""
    return [n for n in _pair_listing() if parse_data(n) is not None]


def pair_points(min_m_xx=PAIR_MIN_M_XX, final_state="4Mu"):
    """Signal mass points of the pair campaign as `{(m_xx, m_zd): [samples]}`.

    Every lifetime of a point is grouped together: none of the quantities in this section
    changes with the proper lifetime once the selection has been applied, and summing
    them is what makes the mass dependence measurable.
    """
    points: dict = {}
    for name in _pair_listing():
        info = parse_signal(name)
        if info is None or info["final_state"] != final_state:
            continue
        if info["m_xx"] < min_m_xx:
            continue
        points.setdefault((info["m_xx"], info["m_zd"]), []).append(name)
    return {key: sorted(value) for key, value in sorted(points.items())}


def pair_counts(samples, hist_name, channel):
    """Raw bin counts of one pair histogram, summed over samples, channel selected.
    Collision data is refused outside the control and validation regions
    (`refuse_blinded_data`)."""
    if isinstance(samples, str):
        samples = [samples]
    refuse_blinded_data(samples, channel)
    total = None
    for name in samples:
        hists = _pair_load(name)
        if hists is None or hist_name not in hists:
            continue
        histogram = hists[hist_name]
        if channel not in list(histogram.axes["channel"]):
            continue
        values = np.asarray(histogram[{"channel": channel}].values(), dtype=float)
        total = values if total is None else total + values
    if total is None:
        note(f"no pair input for {hist_name} in channel {channel}")
    return total


def pair_edges(hist_name, axis):
    """Bin edges of one axis of a pair histogram, read from the first sample that has it."""
    for name in _pair_listing():
        hists = _pair_load(name)
        if hists is not None and hist_name in hists:
            return np.asarray(hists[hist_name].axes[axis].edges, dtype=float)
    return None


def rebin2d(counts, factor=2):
    """Merge `factor` adjacent bins on both axes of a two-dimensional count array."""
    rows = counts.shape[0] // factor
    cols = counts.shape[1] // factor
    return counts[: rows * factor, : cols * factor].reshape(
        rows, factor, cols, factor).sum(axis=(1, 3))


def rebin1d(counts, factor=2):
    """Merge `factor` adjacent bins of a one-dimensional count array."""
    rows = len(counts) // factor
    return np.asarray(counts, dtype=float)[: rows * factor].reshape(rows, factor).sum(axis=1)


def _uniform_centres(n_bins):
    """Centres of `n_bins` equal bins spanning zero to one."""
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    return 0.5 * (edges[:-1] + edges[1:])


# --------------------------------------------------------------------------------------
# is the pair correlated
# --------------------------------------------------------------------------------------

def _pearson(counts):
    """Pearson correlation of a two-dimensional count array evaluated at bin centres."""
    total = float(counts.sum())
    if total <= 1:
        return float("nan")
    x = _uniform_centres(counts.shape[0])[:, None]
    y = _uniform_centres(counts.shape[1])[None, :]
    mean_x = float((counts * x).sum() / total)
    mean_y = float((counts * y).sum() / total)
    sigma_x = float(np.sqrt((counts * (x - mean_x) ** 2).sum() / total))
    sigma_y = float(np.sqrt((counts * (y - mean_y) ** 2).sum() / total))
    if sigma_x <= 0 or sigma_y <= 0:
        return float("nan")
    return float((counts * (x - mean_x) * (y - mean_y)).sum() / total / sigma_x / sigma_y)


def plane_correlation(counts, replicas=2000, seed=20260901):
    """Correlation of the two decay angles of an event, with a bootstrap uncertainty.

    The coefficient is the ordinary linear one evaluated at the bin centres, which is what
    a binned distribution allows.  Its uncertainty comes from drawing `replicas` samples
    of the observed size from the observed distribution, which assumes nothing about the
    shape; the analytic (1 - r^2) / sqrt(N - 1) is returned beside it and agrees with the
    bootstrap to a few percent on every set here.
    """
    counts = np.asarray(counts, dtype=float)
    total = float(counts.sum())
    r = _pearson(counts)
    if not np.isfinite(r) or total <= 2:
        return {"r": float("nan"), "sigma": float("nan"), "sigma_analytic": float("nan"),
                "n": total}
    rng = np.random.default_rng(seed)
    probabilities = (counts / total).ravel()
    draws = rng.multinomial(int(round(total)), probabilities, size=replicas)
    spread = np.std([_pearson(d.reshape(counts.shape).astype(float)) for d in draws])
    return {"r": r, "sigma": float(spread),
            "sigma_analytic": float((1.0 - r * r) / np.sqrt(total - 1.0)), "n": total}


def plane_independence(counts):
    """Ratio of the observed plane to the product of its own marginals.

    A ratio of one in every cell is what independence of the two angles looks like.  The
    uncertainty quoted per cell is the Poisson one of its own content, which is the
    dominant term; the chi-square is the standard test of independence of a contingency
    table, on (rows - 1) x (columns - 1) degrees of freedom.
    """
    from scipy.stats import chi2 as _chi2_distribution  # the rest of the study needs no scipy

    counts = np.asarray(counts, dtype=float)
    total = float(counts.sum())
    rows = counts.sum(axis=1)
    columns = counts.sum(axis=0)
    expected = np.outer(rows, columns) / total if total > 0 else np.zeros_like(counts)
    filled = expected > 0
    ratio = np.divide(counts, expected, out=np.full_like(counts, np.nan), where=filled)
    error = np.divide(ratio, np.sqrt(np.maximum(counts, 1.0)),
                      out=np.full_like(counts, np.nan), where=filled & (counts > 0))
    chi2 = float((((counts - expected) ** 2 / np.where(filled, expected, 1.0))[filled]).sum())
    ndof = (counts.shape[0] - 1) * (counts.shape[1] - 1)
    return {"ratio": ratio, "error": error, "chi2": chi2, "ndof": ndof,
            "p_value": float(_chi2_distribution.sf(chi2, ndof)), "n": total,
            "counts": counts, "expected": expected}


def joint_vs_product(counts, thresholds=(0.2, 0.4, 0.6, 0.8), replicas=1500,
                     seed=20260905):
    """How far a cut on both lepton jets misses the square of a one-jet rate.

    "The one-jet rate" is two different numbers, and at a tight threshold they disagree,
    so both are returned rather than one of them silently standing for the other:

      * `ratio_ordered` divides the observed both-pass fraction by the product of the two
        ORDERED marginals, the leading jet's pass rate times the subleading jet's.  That
        is the comparison for a prediction that keeps the two jets apart, and at a
        threshold of 0.8 those two marginals are far apart, because the harder jet of a
        fake pair keeps its asymmetric muon pairs more often than the softer one.
      * `ratio_inclusive` divides it by the square of ONE pooled per-jet efficiency,
        measured on both jets of every event together.  That is the comparison for the
        usual shorthand, a single per-jet number squared.

    The two agree to three decimals up to a threshold of 0.6 and part company at 0.8,
    where the joint bin holds a few tens of events.  Both carry a bootstrap uncertainty
    drawn from the observed plane at the observed sample size, because at that threshold
    the uncertainty is the size of the effect.
    """
    counts = np.asarray(counts, dtype=float)
    total = float(counts.sum())
    edges = np.linspace(0.0, 1.0, counts.shape[0] + 1)
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(int(round(total)), (counts / total).ravel(),
                            size=replicas).reshape(-1, *counts.shape).astype(float)
    rows = []
    for threshold in thresholds:
        index = int(np.argmin(np.abs(edges - threshold)))

        def measure(table, index=index):
            weight = float(table.sum())
            joint = float(table[index:, index:].sum()) / weight
            lead = float(table[index:, :].sum()) / weight
            sub = float(table[:, index:].sum()) / weight
            pooled = 0.5 * (lead + sub)
            return (joint, lead, sub, pooled,
                    joint / (lead * sub) if lead * sub > 0 else float("nan"),
                    joint / pooled ** 2 if pooled > 0 else float("nan"))

        joint, lead, sub, pooled, ratio_ordered, ratio_inclusive = measure(counts)
        spread_ordered = float(np.std([measure(d)[4] for d in draws]))
        spread_inclusive = float(np.std([measure(d)[5] for d in draws]))
        rows.append({"threshold": float(edges[index]), "joint": joint,
                     "eff_lead": lead, "eff_sub": sub, "eff_inclusive": pooled,
                     "product": lead * sub, "product_inclusive": pooled ** 2,
                     "ratio": ratio_ordered, "ratio_ordered": ratio_ordered,
                     "sigma_ordered": spread_ordered,
                     "ratio_inclusive": ratio_inclusive,
                     "sigma_inclusive": spread_inclusive,
                     "n_joint": float(counts[index:, index:].sum())})
    return rows


def pair_correlation_table(channels=("4mu", "4mu_sr"), min_m_xx=PAIR_MIN_M_XX):
    """Print the correlation of the two decay angles for data and each signal point."""
    rows = []
    short = {PAIR_VR: "4mu VR", "4mu": "4mu", "4mu_sr": "4mu SR",
             "4mu_noMuPtCut": "4mu no mu pT"}
    sets = [("2018 data", pair_data_samples(), PAIR_VR)]
    for (m_xx, m_zd), samples in pair_points(min_m_xx).items():
        for channel in channels:
            sets.append((f"signal m_XX={m_xx:.0f} m_Zd={m_zd}", samples, channel))
    print(f"{'set':30s} {'channel':>8s} {'N events':>10s} {'r':>9s} {'sigma(r)':>9s} "
          f"{'chi2/ndof':>13s} {'p':>10s}")
    print("-" * 92)
    for label, samples, channel in sets:
        fine = pair_counts(samples, PAIR_PLANE, channel)
        if fine is None:
            continue
        correlation = plane_correlation(fine)
        independence = plane_independence(rebin2d(fine, 2))
        print(f"{label:30s} {short.get(channel, channel):>8s} "
              f"{correlation['n']:10.0f} {correlation['r']:+9.4f} "
              f"{correlation['sigma']:9.4f} "
              f"{independence['chi2']:7.1f}/{independence['ndof']:<4d} "
              f"{independence['p_value']:10.2e}")
        rows.append({"label": label, "channel": channel, **correlation,
                     "chi2": independence["chi2"], "ndof": independence["ndof"],
                     "p_value": independence["p_value"]})
    return rows


def pair_factorisation_table(thresholds=(0.2, 0.4, 0.6, 0.8)):
    """Print how much a both-lepton-jet cut misses the square of a one-jet rate.

    Both definitions of "the one-jet rate" are printed side by side -- see
    `joint_vs_product` -- because they answer different questions and differ by more than
    either uncertainty only at the tightest threshold.
    """
    counts = pair_counts(pair_data_samples(), PAIR_PLANE, PAIR_VR)
    if counts is None:
        return []
    rows = joint_vs_product(counts, thresholds)
    print("2018 data, four-muon inverted-displacement validation region: cut on both mu-LJs")
    print("ordered   = joint / (leading-jet efficiency x subleading-jet efficiency)")
    print("inclusive = joint / (one pooled per-jet efficiency)^2")
    print(f"{'threshold':>10s} {'both pass':>10s} {'eff lead':>9s} {'eff sub':>8s} "
          f"{'ordered':>17s} {'eff pooled':>11s} {'inclusive':>17s} {'events kept':>12s}")
    print("-" * 104)
    for row in rows:
        ordered = f"{row['ratio_ordered']:.3f} +- {row['sigma_ordered']:.3f}"
        inclusive = f"{row['ratio_inclusive']:.3f} +- {row['sigma_inclusive']:.3f}"
        print(f"  |cos| > {row['threshold']:.1f} {row['joint']:10.4f} {row['eff_lead']:9.4f} "
              f"{row['eff_sub']:8.4f} {ordered:>17s} {row['eff_inclusive']:11.4f} "
              f"{inclusive:>17s} {row['n_joint']:12.0f}")
    return rows


# --------------------------------------------------------------------------------------
# does it separate inside a mass window
# --------------------------------------------------------------------------------------

def window_slice(edges, low, high):
    """Indices of the bins of `edges` that lie inside [low, high)."""
    return [i for i in range(len(edges) - 1) if edges[i] >= low and edges[i + 1] <= high]


def window_counts(samples, channel, low, high):
    """Decay-angle counts of the pairs whose event mass falls in a window."""
    counts = pair_counts(samples, PAIR_MJJ, channel)
    if counts is None:
        return None
    edges = pair_edges(PAIR_MJJ, "mjj")
    columns = window_slice(edges, low, high)
    if not columns:
        return None
    return counts[:, columns].sum(axis=1)


def window_label(low, high):
    """A short label for a mass window."""
    if high >= 1.0e5:
        return f"> {low:.0f} GeV"
    return f"{low:.0f}-{high:.0f} GeV"


def cut_efficiency(counts, index):
    """Fraction of entries at or above a bin index, with its Clopper-Pearson interval."""
    total = float(np.asarray(counts, dtype=float).sum())
    passed = float(np.asarray(counts, dtype=float)[index:].sum())
    efficiency, low, high = clopper_pearson(passed, max(total, 1.0))
    return float(efficiency[0]), float(low[0]), float(high[0]), total


def cut_gain(signal, background, index):
    """Gain factor of a one-sided cut on the angle, signal efficiency over root background.

    The band comes from the background interval alone: the simulated signal samples hold
    three to four orders of magnitude more entries than the data window, so their own
    interval is invisible beside it.
    """
    eff_s, _, _, n_s = cut_efficiency(signal, index)
    eff_b, low_b, high_b, n_b = cut_efficiency(background, index)
    if eff_b <= 0:
        return {"gain": float("nan"), "low": float("nan"), "high": float("nan"),
                "eff_s": eff_s, "eff_b": eff_b, "n_s": n_s, "n_b": n_b}
    return {"gain": eff_s / np.sqrt(eff_b),
            "low": eff_s / np.sqrt(min(eff_b + high_b, 1.0)),
            "high": eff_s / np.sqrt(max(eff_b - low_b, 1.0e-9)),
            "eff_s": eff_s, "eff_b": eff_b, "n_s": n_s, "n_b": n_b}


def _shape_gain(signal, background):
    """sqrt(sum_i f_s,i^2 / f_b,i) for two count arrays."""
    signal = np.asarray(signal, dtype=float)
    background = np.asarray(background, dtype=float)
    if signal.sum() <= 0 or background.sum() <= 0:
        return float("nan")
    f_s = signal / signal.sum()
    f_b = background / background.sum()
    filled = f_b > 0
    return float(np.sqrt((f_s[filled] ** 2 / f_b[filled]).sum()))


def shape_gain(signal, background, replicas=1500, seed=20260902, max_signal=2_000_000):
    """How much the whole shape could add over the inclusive count, with its bias.

    The figure of merit is sqrt(sum_i f_s,i^2 / f_b,i): the ratio of the signal over root
    background a binned likelihood over the shape reaches to the one the inclusive count
    reaches, in the small-signal Gaussian limit.  It is bounded below by one, so a value
    of 1.10 means the variable could add ten percent inside the window and 1.00 means it
    adds nothing.

    The estimator is biased upwards when the background sample is small -- noise in a
    background bin is read as separation -- so both samples are resampled and the
    bootstrap mean is returned beside the point estimate.  A bootstrap mean well above the
    point estimate is a warning that the point estimate is itself inflated; the
    bias-corrected value is twice the estimate minus that mean.
    """
    signal = np.asarray(signal, dtype=float)
    background = np.asarray(background, dtype=float)
    value = _shape_gain(signal, background)
    if not np.isfinite(value):
        return {"gain": value, "sigma": float("nan"), "bootstrap_mean": float("nan"),
                "corrected": float("nan"), "n_b": float(background.sum())}
    rng = np.random.default_rng(seed)
    n_b = int(round(background.sum()))
    n_s = int(min(round(signal.sum()), max_signal))
    p_b = (background / background.sum()).ravel()
    p_s = (signal / signal.sum()).ravel()
    replica_b = rng.multinomial(n_b, p_b, size=replicas)
    replica_s = rng.multinomial(n_s, p_s, size=replicas)
    values = [_shape_gain(s.reshape(signal.shape).astype(float),
                          b.reshape(background.shape).astype(float))
              for s, b in zip(replica_s, replica_b)]
    values = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    mean = float(values.mean()) if values.size else float("nan")
    return {"gain": value, "sigma": float(values.std()) if values.size else float("nan"),
            "bootstrap_mean": mean, "corrected": 2.0 * value - mean,
            "n_b": float(background.sum())}


def pair_window_table(min_entries=MIN_DENOMINATOR):
    """Print the data shape of the angle in each window of the event mass."""
    counts = pair_counts(pair_data_samples(), PAIR_MJJ, PAIR_VR)
    if counts is None:
        return []
    edges = pair_edges(PAIR_MJJ, "mjj")
    centres = _uniform_centres(counts.shape[0])
    pooled = counts[:, window_slice(edges, 150.0, 1.0e6)].sum(axis=1)
    pooled_shape = pooled / pooled.sum()
    rows = []
    print("2018 data, four-muon inverted-displacement validation region, one entry per mu-LJ")
    print(f"{'m(muLJ,muLJ)':>16s} {'pairs':>7s}   {'shape of |cos(theta*)| in 5 bins':<34s} "
          f"{'mean':>13s} {'chi2/4 vs all':>14s}")
    print("-" * 96)
    for index in range(len(edges) - 1):
        if edges[index] < 150.0:
            continue
        column = counts[:, index]
        total = float(column.sum())
        if total < min_entries:
            print(f"{window_label(edges[index], edges[index + 1]):>16s} {total:7.0f}   "
                  "too few pairs to measure a shape")
            rows.append({"low": edges[index], "high": edges[index + 1], "n": total,
                         "shape": None})
            continue
        shape = column / total
        mean = float((shape * centres).sum())
        error = float(np.sqrt(((centres - mean) ** 2 * shape).sum() / total))
        expected = pooled_shape * total
        chi2 = float((((column - expected) ** 2 / np.where(expected > 0, expected, 1))
                      [expected > 0]).sum())
        print(f"{window_label(edges[index], edges[index + 1]):>16s} {total:7.0f}   "
              f"{' '.join(f'{v:.3f}' for v in shape):<34s} {mean:.3f} +- {error:.3f} "
              f"{chi2:14.1f}")
        rows.append({"low": edges[index], "high": edges[index + 1], "n": total,
                     "shape": shape, "mean": mean, "error": error, "chi2": chi2})
    print()
    print(f"{'all m > 150 GeV':>16s} {pooled.sum():7.0f}   "
          f"{' '.join(f'{v:.3f}' for v in pooled_shape):<34s} "
          f"{float((pooled_shape * centres).sum()):.3f}")
    return rows


def pair_gain_table(channel="4mu_sr", thresholds=(0.2, 0.4, 0.6, 0.8),
                    min_m_xx=PAIR_MIN_M_XX):
    """Print, per signal point, what the angle could add inside its own mass window."""
    data = pair_data_samples()
    rows = []
    print(f"signal channel {channel}; fakes from the four-muon inverted-displacement VR")
    print(f"{'m_XX':>6s} {'m_Zd':>5s} {'window':>12s} {'fake pairs':>11s} "
          f"{'best cut':>22s} {'shape gain':>22s} {'vs pooled fakes':>18s}")
    print("-" * 106)
    for (m_xx, m_zd), samples in pair_points(min_m_xx).items():
        low, high = PAIR_WINDOW.get(m_xx, (150.0, 1.0e6))
        signal = window_counts(samples, channel, low, high)
        background = window_counts(data, PAIR_VR, low, high)
        pooled = window_counts(data, PAIR_VR, *PAIR_POOLED_WINDOW)
        if signal is None or background is None:
            continue
        scan = []
        for index, threshold in enumerate(thresholds, start=1):
            gain = cut_gain(signal, background, index)
            gain["threshold"] = threshold
            scan.append(gain)
        best = max(scan, key=lambda g: g["gain"] if np.isfinite(g["gain"]) else -np.inf)
        matched = shape_gain(signal, background)
        pooled_gain = shape_gain(signal, pooled)
        best_text = f"|cos| > {best['threshold']:.1f}: {best['gain']:.2f}"
        matched_text = f"{matched['gain']:.3f} +- {matched['sigma']:.3f}"
        pooled_text = f"{pooled_gain['gain']:.3f} +- {pooled_gain['sigma']:.3f}"
        print(f"{m_xx:6.0f} {m_zd:5.2f} {window_label(low, high):>12s} "
              f"{background.sum():11.0f} {best_text:>22s} {matched_text:>22s} "
              f"{pooled_text:>18s}")
        rows.append({"m_xx": m_xx, "m_zd": m_zd, "low": low, "high": high,
                     "signal": signal, "background": background, "best_cut": best,
                     "matched": matched, "pooled": pooled_gain})
    print()
    print("shape gain is sqrt(sum f_s^2 / f_b) over the five bins of the angle: the factor by")
    print("which a shape fit beats the inclusive count in the same window.  One means nothing")
    print("is added.  The bias of the estimator at these sample sizes is printed below.")
    for row in rows:
        matched = row["matched"]
        print(f"  m_XX={row['m_xx']:.0f} m_Zd={row['m_zd']}: point {matched['gain']:.3f}, "
              f"bootstrap mean {matched['bootstrap_mean']:.3f}, "
              f"bias-corrected {matched['corrected']:.3f} "
              f"(on {matched['n_b']:.0f} fake pairs)")
    return rows


def pair_plane_gain_table(channel="4mu_sr", min_m_xx=PAIR_MIN_M_XX):
    """Print the shape gain of the whole two-jet plane against the one-jet gain squared."""
    background = rebin2d(pair_counts(pair_data_samples(), PAIR_PLANE, PAIR_VR), 2)
    print("event level, all event masses: the two-dimensional plane against two one-jet cuts")
    print(f"{'m_XX':>6s} {'m_Zd':>5s} {'plane gain':>22s} {'one jet':>10s} "
          f"{'one jet squared':>16s} {'plane / squared':>16s}")
    print("-" * 82)
    rows = []
    for (m_xx, m_zd), samples in pair_points(min_m_xx).items():
        signal = rebin2d(pair_counts(samples, PAIR_PLANE, channel), 2)
        plane = shape_gain(signal, background)
        one_jet = _shape_gain(signal.sum(axis=1) + signal.sum(axis=0),
                              background.sum(axis=1) + background.sum(axis=0))
        plane_text = f"{plane['gain']:.3f} +- {plane['sigma']:.3f}"
        print(f"{m_xx:6.0f} {m_zd:5.2f} {plane_text:>22s} {one_jet:10.3f} "
              f"{one_jet ** 2:16.3f} {plane['gain'] / one_jet ** 2:16.3f}")
        rows.append({"m_xx": m_xx, "m_zd": m_zd, "plane": plane, "one_jet": one_jet})
    return rows


# --------------------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------------------

def _draw_ratio_map(ax, result, canvas, span=0.45, headroom=1.42, edges=None):
    """Draw a ratio-to-marginals map with the value written into every filled cell.

    The vertical range is extended past the map so that the annotation naming the set sits
    on empty canvas: a legend box over a map of numbers hides exactly the cells a reader
    came for.  `edges` defaults to zero to one, which is the range of an absolute cosine;
    a signed angle passes its own edges in.
    """
    ratio = result["ratio"]
    counts = result["counts"]
    if edges is None:
        edges = np.linspace(0.0, 1.0, ratio.shape[0] + 1)
    edges = np.asarray(edges, dtype=float)
    mesh = ax.pcolormesh(edges, edges, ratio.T, cmap="RdBu_r", shading="flat",
                         vmin=1.0 - span, vmax=1.0 + span)
    centres = 0.5 * (edges[:-1] + edges[1:])
    for i, x in enumerate(centres):
        for j, y in enumerate(centres):
            if not np.isfinite(ratio[i, j]) or counts[i, j] <= 0:
                continue
            ax.text(x, y, f"{ratio[i, j]:.2f}", ha="center", va="center",
                    fontsize=canvas.font * 0.90, color="#111111")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(edges[0], edges[0] + (edges[-1] - edges[0]) * headroom)
    ax.set_yticks(np.round(edges, 2))
    return mesh


def fig_pair_plane(channel="4mu_sr", point=(500.0, 1.2), stem="i1_pair_plane"):
    """The two decay angles of an event, as a ratio to the product of their marginals.

    One mass point is drawn, not a sum over the grid: summing points whose one-jet shapes
    differ builds a correlation into the sum that neither point has, and the table beside
    the figure carries every point separately.
    """
    data = pair_data_samples()
    signal_samples_list = pair_points(0.0).get(point, [])
    if not signal_samples_list:
        note(f"no pair samples for m_XX = {point[0]} GeV, m_Zd = {point[1]} GeV")
    canvas = Canvas(1, 2, panel_w=5.7, panel_h=5.0, base_pt=9.0)
    meshes = []
    sets = [
        (0, data, PAIR_VR, "2018 data, 4$\\mu$ VR, inv. displaced"),
        (1, signal_samples_list, channel,
         f"signal 4$\\mu$, $m_{{XX}}$ = {point[0]:.0f} GeV,\n"
         f"$m_{{Z_d}}$ = {point[1]} GeV, signal-region channel"),
    ]
    for column, samples, region, label in sets:
        ax = canvas.ax(0, column)
        fine = pair_counts(samples, PAIR_PLANE, region)
        if fine is None:
            canvas.hide(ax, "no input")
            continue
        result = plane_independence(rebin2d(fine, 2))
        correlation = plane_correlation(fine)
        meshes.append(_draw_ratio_map(ax, result, canvas))
        canvas.labels(ax, r"$|\cos\theta^*|$, leading $\mu$-LJ",
                      r"$|\cos\theta^*|$, subleading $\mu$-LJ" if column == 0 else None)
        ax.text(0.02, 0.985,
                f"{label}\n"
                f"$r$ = {correlation['r']:+.3f} $\pm$ {correlation['sigma']:.3f}, "
                f"{result['n']:.0f} events\n"
                f"independence $\chi^2$/dof = {result['chi2']:.0f}/{result['ndof']}",
                transform=ax.transAxes, ha="left", va="top", fontsize=canvas.font,
                linespacing=1.35)
    if meshes:
        colour_bar = canvas.fig.colorbar(meshes[-1], ax=list(canvas.axes[0]), pad=0.02,
                                         fraction=0.05)
        colour_bar.ax.tick_params(labelsize=canvas.tick_font)
        colour_bar.set_label("observed / product of marginals", fontsize=canvas.font)
    canvas.finish(data=True, lumi_fb=LUMI_FB)
    canvas.save(stem)
    return canvas


def _draw_counts_shape(ax, counts, label, colour, canvas, marker="o", linestyle="-",
                       as_points=True, min_per_bin=MIN_ENTRIES_PER_BIN, offset=0.0):
    """Draw a count array as a unit-area density with Poisson errors."""
    counts = np.asarray(counts, dtype=float)
    total = float(counts.sum())
    n_bins = counts.shape[0]
    if total < MIN_SHAPE_ENTRIES or total < min_per_bin * n_bins:
        note(f"{_plain(label)}: {total:.0f} entries over {n_bins} bins, below what this "
             "study draws a shape at; not drawn")
        return None
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    widths = np.diff(edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    density = counts / widths / total
    errors = poisson_errors(counts) / widths / total
    if as_points:
        # several series share one binning in these panels, so a small horizontal offset
        # keeps their error bars from being drawn on top of each other
        ax.errorbar(centres + offset, density, yerr=errors, fmt=marker, color=colour,
                    label=label, markersize=canvas.font * 0.20, elinewidth=1.2, capsize=0)
    else:
        ax.stairs(density, edges, color=colour, linewidth=2.0, linestyle=linestyle,
                  label=label, baseline=None)
    return total


def fig_pair_conditionals(channel="4mu_sr", point=(500.0, 1.2),
                          stem="i2_pair_conditionals"):
    """The subleading jet's angle in slices of the leading jet's, for data and signal."""
    data = pair_data_samples()
    signal_samples_list = pair_points(0.0).get(point, [])
    slices = [(0, 3, "0.0 - 0.3"), (3, 7, "0.3 - 0.7"), (7, 10, "0.7 - 1.0")]
    canvas = Canvas(1, 2, panel_w=5.7, panel_h=4.8, base_pt=9.0)
    entries = []
    for column, (samples, region, title) in enumerate([
            (data, PAIR_VR, "2018 data, 4$\\mu$ VR, inv. displaced"),
            (signal_samples_list, channel,
             f"signal 4$\\mu$, $m_{{XX}}$ = {point[0]:.0f}, "
             f"$m_{{Z_d}}$ = {point[1]} GeV, SR")]):
        ax = canvas.ax(0, column)
        fine = pair_counts(samples, PAIR_PLANE, region)
        if fine is None:
            canvas.hide(ax, "no input")
            continue
        for index, (low, high, label) in enumerate(slices):
            band = rebin1d(fine[low:high, :].sum(axis=0), 2)
            _draw_counts_shape(ax, band, label, _SLICE_COLOURS[index], canvas,
                               offset=0.02 * (index - 1))
        canvas.labels(ax, r"$|\cos\theta^*|$, subleading $\mu$-LJ",
                      "events / unit angle (normalised)" if column == 0 else None)
        legend = canvas.legend(
            ax, "upper left", ncol=3, over_data=True,
            title=title + "\n" + r"leading $\mu$-LJ $|\cos\theta^*|$ in")
        entries.append((ax, legend))
    canvas.clear_legends(entries)
    canvas.share_ylimits(entries)
    canvas.finish(data=True, lumi_fb=LUMI_FB)
    canvas.save(stem)
    return canvas


def fig_window_shapes(channel="4mu_sr", min_m_xx=PAIR_MIN_M_XX, stem="i3_window_shapes"):
    """The angle's shape window by window in data, and against signal where signal lives."""
    data = pair_data_samples()
    counts = pair_counts(data, PAIR_MJJ, PAIR_VR)
    edges = pair_edges(PAIR_MJJ, "mjj")
    canvas = Canvas(1, 2, panel_w=5.7, panel_h=4.8, base_pt=9.0)
    entries = []

    ax = canvas.ax(0, 0)
    drawn = 0
    for index in range(len(edges) - 1):
        if edges[index] < 150.0:
            continue
        column = counts[:, index]
        if column.sum() < MIN_SHAPE_ENTRIES:
            continue
        _draw_counts_shape(ax, column, window_label(edges[index], edges[index + 1]),
                           _SLICE_COLOURS[drawn % len(_SLICE_COLOURS)], canvas,
                           offset=0.016 * (drawn - 2))
        drawn += 1
    canvas.labels(ax, r"$|\cos\theta^*|$ (per $\mu$-LJ)", "entries / unit angle (normalised)")
    legend = canvas.legend(ax, "upper left", ncol=2, over_data=True,
                           title="2018 data, 4$\\mu$ VR: $m(\\mu$LJ, $\\mu$LJ)")
    entries.append((ax, legend))

    # windows the table shows but this panel does not: name them and say why, so the
    # figure and the table above it cannot be read as disagreeing
    skipped = [(edges[i], edges[i + 1], float(counts[:, i].sum()))
               for i in range(len(edges) - 1)
               if edges[i] >= 150.0 and 0 < counts[:, i].sum() < MIN_SHAPE_ENTRIES]
    if skipped:
        canvas.stamp(canvas.ax(0, 0),
                     f"windows above {skipped[0][0]:.0f} GeV hold\n"
                     + ", ".join(f"{n:.0f}" for _, _, n in skipped)
                     + f" pairs, under the {MIN_SHAPE_ENTRIES}-entry\n"
                       "floor: see the table above",
                     loc="lower left")

    ax = canvas.ax(0, 1)
    for index, (m_xx, low_high) in enumerate(sorted(PAIR_WINDOW.items())):
        low, high = low_high
        samples = [name for (mass, m_zd), names in pair_points(min_m_xx).items()
                   if mass == m_xx and m_zd >= 1.0 for name in names]
        colour = COLOURS["4Mu"] if index == 0 else COLOURS["gen"]
        _draw_counts_shape(ax, window_counts(samples, channel, low, high),
                           f"signal $m_{{XX}}$ = {m_xx:.0f} GeV", colour, canvas,
                           as_points=False, linestyle="-" if index == 0 else "--")
        fakes = window_counts(data, PAIR_VR, low, high)
        _draw_counts_shape(ax, fakes, f"data, {window_label(low, high)}", colour, canvas,
                           marker="s" if index == 0 else "^")
        if fakes is not None and fakes.sum() < MIN_SHAPE_ENTRIES:
            canvas.stamp(ax, f"data at {window_label(low, high)}: {fakes.sum():.0f} pairs,\n"
                             "too few to draw a shape", loc="lower left")
    canvas.labels(ax, r"$|\cos\theta^*|$ (per $\mu$-LJ)", None)
    legend = canvas.legend(ax, "upper left", title="inside the window the signal fills",
                           over_data=True)
    entries.append((ax, legend))

    canvas.clear_legends(entries)
    canvas.share_ylimits(entries)
    canvas.finish(data=True, lumi_fb=LUMI_FB)
    canvas.save(stem)
    return canvas


def fig_window_gain(rows, thresholds=(0.2, 0.4, 0.6, 0.8), stem="i4_window_gain"):
    """What a cut, and what the whole shape, would add inside the signal's own window."""
    canvas = Canvas(1, 2, panel_w=5.7, panel_h=4.8, base_pt=9.0)
    entries = []

    ax = canvas.ax(0, 0)
    for index, row in enumerate(rows):
        colour = _SLICE_COLOURS[index % len(_SLICE_COLOURS)]
        values, lows, highs = [], [], []
        for position in range(1, len(thresholds) + 1):
            gain = cut_gain(row["signal"], row["background"], position)
            values.append(gain["gain"])
            lows.append(gain["gain"] - gain["low"])
            highs.append(gain["high"] - gain["gain"])
        ax.errorbar(np.asarray(thresholds) + 0.004 * (index - 2.5), values,
                    yerr=[lows, highs], fmt="o-", color=colour,
                    markersize=canvas.font * 0.18, elinewidth=1.2, capsize=0,
                    label=f"{row['m_xx']:.0f} / {row['m_zd']}")
    ax.axhline(1.0, color="#666666", linestyle=":", linewidth=1.5)
    canvas.labels(ax, r"cut at $|\cos\theta^*| >$",
                  r"$\epsilon_{s} / \sqrt{\epsilon_{b}}$ per $\mu$-LJ")
    legend = canvas.legend(ax, "upper left", ncol=2, over_data=True,
                           title="one cut, inside the signal's own window\n"
                                 "$m_{XX}$ / $m_{Z_d}$ [GeV]")
    entries.append((ax, legend))

    ax = canvas.ax(0, 1)
    positions = np.arange(len(rows), dtype=float)
    for offset, key, label, colour, marker in [
            (-0.16, "matched", "fakes from the same window", COLOURS["data"], "o"),
            (0.16, "pooled", "fakes from $m(\\mu$LJ, $\\mu$LJ$) > 300$ GeV",
             COLOURS["background"], "s")]:
        values = [row[key]["gain"] for row in rows]
        errors = [row[key]["sigma"] for row in rows]
        ax.errorbar(positions + offset, values, yerr=errors, fmt=marker, color=colour,
                    markersize=canvas.font * 0.20, elinewidth=1.2, capsize=0, label=label)
    # the point estimate is biased upwards on a small fake sample, and on the 1000 GeV
    # points that bias is larger than the distance of the point from one, so the corrected
    # value is drawn beside it: without it the figure reads as six points above one while
    # the verdict says the heavy points are consistent with adding nothing
    ax.errorbar(positions, [row["matched"]["corrected"] for row in rows],
                yerr=[row["matched"]["sigma"] for row in rows], marker="D",
                markerfacecolor="white", markeredgecolor=COLOURS["data"],
                markeredgewidth=1.6, color=COLOURS["data"], linestyle="none",
                markersize=canvas.font * 0.18, elinewidth=1.2, capsize=0,
                label="same-window fakes, bias corrected")
    ax.axhline(1.0, color="#666666", linestyle=":", linewidth=1.5)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{row['m_xx']:.0f}\n{row['m_zd']}" for row in rows])
    ax.set_xlim(-0.5, len(rows) - 0.5)
    canvas.labels(ax, r"$m_{XX}$ [GeV] and $m_{Z_d}$ [GeV]",
                  r"shape gain $\sqrt{\sum_i f_{s,i}^2 / f_{b,i}}$")
    legend = canvas.legend(ax, "upper left", title="the whole shape, inside the window",
                           over_data=True)
    entries.append((ax, legend))

    sizes = {}
    for row in rows:
        sizes.setdefault(row["m_xx"], row["matched"]["n_b"])
    canvas.stamp(canvas.ax(0, 1),
                 "fake pairs behind each column:\n"
                 + ", ".join(f"$m_{{XX}}$ = {mass:.0f} GeV: {n:.0f}"
                             for mass, n in sorted(sizes.items())),
                 loc="lower left")

    canvas.clear_legends(entries)
    canvas.finish(data=True, lumi_fb=LUMI_FB)
    canvas.save(stem)
    return canvas


# ======================================================================================
# section j -- the lepton-jet decay planes
#
# Inputs are a third campaign, LjAnglePlanes_Run2018_v1, merged one file per sample.  It
# adds the signed decay angle, the electron pair of an egamma-type lepton jet, and the
# event-level observables of the two leading qualifying lepton jets: the two helicity
# angles, the azimuthal angle between the two decay planes, the transverse balance and
# the missing momentum.  Its channel list and its sample set are both wider than the
# earlier campaigns', and collision data was written in control and validation regions
# only, which `planes_counts` enforces rather than assumes.
# ======================================================================================

PLANES_INPUT = "/uscms_data/d3/murtazas/lj_angle_planes_merged"

#: campaign tag, for the provenance line of the notebook.
PLANES_CAMPAIGN = "LjAnglePlanes_Run2018_v1"

PLANES_CR_1MU = "data_control_region_1muLj_spread_cosAlpha_mu_veto"
PLANES_CR_1EGM = "data_control_region_1egmLj_cosmic_veto"
PLANES_CR_2LJ = "cr_2lj_invDisplaced_spread_cosAlpha_mu_veto"
PLANES_VR_4MU = "test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto"
PLANES_VR_2MU2E = "test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto"

#: the ONLY channels collision data may be read in.  Everything else this campaign wrote
#: is the signal region or a variant of it and stays blind; `planes_counts` refuses to
#: open a data file outside this set rather than trusting the caller.
PLANES_DATA_CHANNELS = frozenset({PLANES_CR_1MU, PLANES_CR_1EGM, PLANES_CR_2LJ,
                                  PLANES_VR_4MU, PLANES_VR_2MU2E})

PLANES_REGION_LABEL = {
    PLANES_CR_1MU: "one $\\mu$-LJ CR",
    PLANES_CR_1EGM: "one $e\\gamma$-LJ CR",
    PLANES_CR_2LJ: "two-LJ inv. displaced CR",
    PLANES_VR_4MU: "4$\\mu$ inv. displaced VR",
    PLANES_VR_2MU2E: "2$\\mu$2e inv. displaced VR",
    "4mu": "4$\\mu$ channel",
    "4mu_sr": "4$\\mu$ signal region",
    "2mu2e": "2$\\mu$2e channel",
    "2mu2e_sr": "2$\\mu$2e signal region",
}

#: the search uses bound-state masses at or above this.  The grid this campaign ran holds
#: 100, 500 and 1000 GeV only, so "at or above 200" means 500 and 1000 here and the 200 to
#: 250 GeV range at the low edge of the search is not in it.
PLANES_MIN_M_XX = 200.0

#: the coefficient of cos(2 Phi) an unsculpted CP-even parent gives, from the transverse
#: polarisation of the two dark photons: dN/dPhi is proportional to 1 + A cos(2 Phi) with
#: A = [beta*^2 / (3 - beta*^2)]^2, the product of the two decays' analyzing powers, which is
#: +1/4 for massless leptons (the value drawn), 0.228 for muons at m_Zd = 1.2 GeV and 0.011
#: at 0.25 GeV.  A CP-odd parent of the same mass would give -A.  Neither value is what
#: a reconstructed distribution has to match: the reference is the generator-level curve
#: from the same files, which carries the production filter and the selection with it.
PLANES_A_CP_EVEN = 0.25

#: threshold the missing-momentum numbers are quoted at, in GeV.  Chosen because it is a
#: bin edge of both the forty-bin one-dimensional axes and the twenty-bin plane, so the
#: one-dimensional and two-dimensional fractions are quoted at the same place.
PLANES_MET_THRESHOLD = 40.0

#: short names for the printed tables.  The channel names are long enough that a column
#: sized for them pushes the numbers off the page.
PLANES_SHORT = {
    PLANES_CR_1MU: "1 muLJ CR",
    PLANES_CR_1EGM: "1 egmLJ CR",
    PLANES_CR_2LJ: "2 LJ invDisp CR",
    PLANES_VR_4MU: "4mu invDisp VR",
    PLANES_VR_2MU2E: "2mu2e invDisp VR",
}

_PLANES_LOADED: dict = {}


def planes_short(channel):
    """A short name for a channel, for the printed tables."""
    return PLANES_SHORT.get(channel, channel)


def configure_planes(input_dir=None):
    """Point the decay-plane loader at a directory of merged plane-campaign outputs."""
    global PLANES_INPUT
    if input_dir is not None and input_dir.rstrip("/") != PLANES_INPUT:
        PLANES_INPUT = input_dir.rstrip("/")
        _PLANES_LOADED.clear()
    return PLANES_INPUT


def _planes_listing():
    """Sample names available in `PLANES_INPUT`, sorted."""
    if PLANES_INPUT.startswith("root://"):
        files = _remote_listing(PLANES_INPUT)
    else:
        files = [os.path.basename(p) for p in glob.glob(os.path.join(PLANES_INPUT, "*.coffea"))]
    return sorted(f[: -len(".coffea")] for f in files)


def _planes_load(name):
    """One plane-campaign sample as `{hists, metadata, is_data, sum_weights, year}`."""
    if name in _PLANES_LOADED:
        return _PLANES_LOADED[name]
    try:
        if PLANES_INPUT.startswith("root://"):
            os.makedirs(CACHE_DIR, exist_ok=True)
            path = os.path.join(CACHE_DIR, f"planes_{name}.coffea")
            if not os.path.exists(path):
                subprocess.run(["xrdcp", "-f", f"{PLANES_INPUT}/{name}.coffea", path], check=True)
        else:
            path = os.path.join(PLANES_INPUT, f"{name}.coffea")
        blob = coffea.util.load(path)
    except (FileNotFoundError, OSError, subprocess.CalledProcessError) as exc:
        note(f"plane sample {name} could not be read ({type(exc).__name__}); skipping it")
        _PLANES_LOADED[name] = None
        return None
    out = blob["out"] if "out" in blob else blob
    if name in out:
        payload = out[name]
    elif len(out) == 1:
        payload = next(iter(out.values()))
    else:
        note(f"plane file for {name} holds {len(out)} samples and none matches its name; skipping")
        _PLANES_LOADED[name] = None
        return None
    meta = payload.get("metadata", {})
    _PLANES_LOADED[name] = {
        "hists": payload["hists"],
        "metadata": meta,
        "is_data": bool(_scalar(meta.get("is_data"), parse_data(name) is not None)),
        "year": str(_scalar(meta.get("year"), YEAR)),
        "sum_weights": float(meta.get("scaled_sum_weights", meta.get("n_evts", 0.0)) or 0.0),
    }
    return _PLANES_LOADED[name]


def planes_data_samples():
    """Collision-data samples of the plane campaign, all eras."""
    return [n for n in _planes_listing() if parse_data(n) is not None]


def planes_background_samples():
    """Simulated background samples of the plane campaign."""
    return [n for n in _planes_listing()
            if parse_data(n) is None and parse_signal(n) is None]


def planes_points(final_state="4Mu", min_m_xx=PLANES_MIN_M_XX):
    """Signal mass points of the plane campaign as `{(m_xx, m_zd): [samples]}`.

    Every proper lifetime of a point is grouped together, as in section i: none of the
    angular quantities changes with the lifetime once a channel has been applied, and
    summing them is what makes the mass dependence measurable.
    """
    points: dict = {}
    for name in _planes_listing():
        info = parse_signal(name)
        if info is None or info["final_state"] != final_state:
            continue
        if info["m_xx"] < min_m_xx:
            continue
        points.setdefault((info["m_xx"], info["m_zd"]), []).append(name)
    return {key: sorted(value) for key, value in sorted(points.items())}


def planes_weight(name):
    """Per-sample lumi x sigma / N for a plane-campaign simulated sample.

    The same normalisation as `lumi_xs_weight`: the denominator is the sum of genWeight over
    the stored events of this campaign's file, as SidmProcessor applies it, with no
    generator filter recorded in the files.  Everything drawn in this section is unit-area,
    where the factor cancels within a sample and only reweights one simulated sample
    against another, which is what a background sum needs.
    """
    record = _planes_load(name)
    if record is None:
        return None
    if record["is_data"]:
        return 1.0
    if sidm_utilities is not None:
        try:
            return sidm_utilities.get_lumixs_weight(name, record["year"], record["sum_weights"])
        except (KeyError, FileNotFoundError):
            pass
    table = load_cross_sections()
    if name in table:
        cross_section = table[name]
    elif name.startswith(("4Mu", "2Mu2E")):
        cross_section = 0.001
    else:
        note(f"no cross section for {name}; plane weight set to 1")
        return 1.0
    lumi = load_run_periods()[record["year"]]["lumi"]
    return lumi * cross_section / record["sum_weights"] if record["sum_weights"] else 1.0


def planes_counts(samples, hist_name, channel, weighted=False):
    """Bin counts of one plane-campaign histogram, summed over samples.

    Raw counts by default, which is what every statistical statement here is built from;
    `weighted=True` applies the per-sample cross-section factor first, which is what a
    sum over simulated background samples needs before its shape means anything.

    BLINDING: this study never reads collision data outside a control or validation
    region.  Passing a data sample with any other channel raises rather than returning
    something a figure could quietly draw.
    """
    if isinstance(samples, str):
        samples = [samples]
    if channel not in PLANES_DATA_CHANNELS:
        opened = [n for n in samples if parse_data(n) is not None]
        if opened:
            raise ValueError(
                f"channel {channel!r} is not one of this study's control or validation "
                f"regions and {len(opened)} collision-data sample(s) were passed to it. "
                "Data is read only in: " + ", ".join(sorted(PLANES_DATA_CHANNELS)))
    total = None
    for name in samples:
        record = _planes_load(name)
        if record is None or hist_name not in record["hists"]:
            continue
        histogram = record["hists"][hist_name]
        if channel not in list(histogram.axes["channel"]):
            continue
        values = np.asarray(histogram[{"channel": channel}].values(), dtype=float)
        if weighted:
            factor = planes_weight(name)
            if factor is None:
                continue
            values = values * factor
        total = values if total is None else total + values
    if total is None:
        note(f"no plane input for {hist_name} in channel {channel}")
    return total


def planes_edges(hist_name, axis):
    """Bin edges of one axis of a plane-campaign histogram."""
    for name in _planes_listing():
        record = _planes_load(name)
        if record is not None and hist_name in record["hists"]:
            return np.asarray(record["hists"][hist_name].axes[axis].edges, dtype=float)
    return None


# --------------------------------------------------------------------------------------
# statistics on plain count arrays
# --------------------------------------------------------------------------------------

def binned_mean(counts, low=0.0, high=1.0):
    """Mean of a uniformly binned distribution and the standard error of that mean."""
    counts = np.asarray(counts, dtype=float)
    total = float(counts.sum())
    if total <= 1:
        return float("nan"), float("nan"), total
    edges = np.linspace(low, high, len(counts) + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    mean = float((counts * centres).sum() / total)
    variance = float((counts * (centres - mean) ** 2).sum() / total)
    return mean, float(np.sqrt(variance / total)), total


def harmonic_moment(counts, low, high, harmonic=2):
    """<cos(n x)> of a uniformly binned angular distribution, with its statistical error.

    For a distribution proportional to 1 + A cos(2 x) over a full period, <cos 2x> = A / 2,
    so the fitted coefficient is twice this moment.  Taking the moment rather than fitting
    uses every entry with no binning-dependent minimisation, and its error is the ordinary
    error of a mean.
    """
    counts = np.asarray(counts, dtype=float)
    total = float(counts.sum())
    if total <= 1:
        return {"moment": float("nan"), "sigma": float("nan"), "coefficient": float("nan"),
                "sigma_coefficient": float("nan"), "n": total}
    edges = np.linspace(low, high, len(counts) + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    weights = np.cos(harmonic * centres)
    mean = float((counts * weights).sum() / total)
    variance = float((counts * weights * weights).sum() / total) - mean * mean
    sigma = float(np.sqrt(max(variance, 0.0) / total))
    return {"moment": mean, "sigma": sigma, "coefficient": 2.0 * mean,
            "sigma_coefficient": 2.0 * sigma, "n": total}


def kept_fraction(counts, threshold, low=0.0, high=1.0, direction="above"):
    """Fraction of entries on one side of a threshold, with a Clopper-Pearson interval."""
    counts = np.asarray(counts, dtype=float)
    total = float(counts.sum())
    if total <= 0:
        return {"fraction": float("nan"), "low": float("nan"), "high": float("nan"), "n": 0.0}
    edges = np.linspace(low, high, len(counts) + 1)
    index = int(np.argmin(np.abs(edges - threshold)))
    passed = float(counts[index:].sum() if direction == "above" else counts[:index].sum())
    value, lower, upper = clopper_pearson(passed, total)
    return {"fraction": float(value[0]), "low": float(lower[0]), "high": float(upper[0]),
            "n": total, "passed": passed, "edge": float(edges[index])}


def separation_curve(signal, background, low=0.0, high=1.0, direction="above"):
    """Signal efficiency and background efficiency of every threshold on a shared axis.

    Returns the two efficiency arrays and the bin edges they are evaluated at, which is
    what a threshold-by-threshold comparison of two shapes needs and all a binned
    distribution supports.
    """
    signal = np.asarray(signal, dtype=float)
    background = np.asarray(background, dtype=float)
    edges = np.linspace(low, high, len(signal) + 1)
    if direction == "above":
        eff_s = np.array([signal[i:].sum() for i in range(len(signal) + 1)]) / max(signal.sum(), 1)
        eff_b = np.array([background[i:].sum()
                          for i in range(len(background) + 1)]) / max(background.sum(), 1)
    else:
        eff_s = np.array([signal[:i].sum() for i in range(len(signal) + 1)]) / max(signal.sum(), 1)
        eff_b = np.array([background[:i].sum()
                          for i in range(len(background) + 1)]) / max(background.sum(), 1)
    return edges, eff_s, eff_b


def _first_key(mapping, *preferred):
    """The first of `preferred` present in `mapping`, else its first key, else None."""
    for key in preferred:
        if key in mapping:
            return key
    return next(iter(mapping), None)


# --------------------------------------------------------------------------------------
# drawing primitives for plain count arrays
# --------------------------------------------------------------------------------------

def draw_counts(ax, counts, edges, label, colour, canvas, as_points=True, linestyle="-",
                marker="o", offset=0.0, raw=None, min_per_bin=MIN_ENTRIES_PER_BIN,
                context=""):
    """Draw a count array over arbitrary edges as a unit-area density.

    `raw` is the unweighted entry count behind `counts`, which is what the drawing gates
    are judged on: a cross-section-weighted sum of simulated samples can look smooth while
    resting on a handful of high-weight events.  Bins below the per-bin floor are dropped
    rather than drawn, and a curve with too few entries anywhere is not drawn at all.
    """
    if counts is None:
        return None
    counts = np.asarray(counts, dtype=float)
    entries = np.asarray(raw if raw is not None else counts, dtype=float)
    total_raw = float(entries.sum())
    n_bins = len(counts)
    where = f" in {context}" if context else ""
    if total_raw < MIN_SHAPE_ENTRIES:
        note(f"{_plain(label)}{where}: {total_raw:.0f} entries over {n_bins} bins, below "
             f"the {MIN_SHAPE_ENTRIES} needed to draw a shape at all; not drawn")
        return None
    if total_raw < min_per_bin * n_bins:
        note(f"{_plain(label)}{where}: {total_raw:.0f} entries over {n_bins} bins, fewer "
             f"than the {min_per_bin} per bin a shape is drawn at; not drawn")
        return None
    edges = np.asarray(edges, dtype=float)
    widths = np.diff(edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    total = float(counts.sum())
    if total <= 0:
        return None
    density = counts / widths / total
    with np.errstate(divide="ignore", invalid="ignore"):
        relative = np.where(entries > 0, 1.0 / np.sqrt(np.maximum(entries, 1.0)), 0.0)
    errors = density * relative
    keep = entries >= min_per_bin
    if not keep.any():
        note(f"{_plain(label)}{where}: no bin holds {min_per_bin} entries; not drawn")
        return None
    span = edges[-1] - edges[0]
    if as_points:
        ax.errorbar(centres[keep] + offset * span, density[keep], yerr=errors[keep],
                    fmt=marker, color=colour, label=label,
                    markersize=canvas.font * 0.20, elinewidth=1.2, capsize=0)
    else:
        ax.stairs(np.where(keep, density, np.nan), edges, color=colour, linewidth=2.0,
                  linestyle=linestyle, label=label, baseline=None)
        ax.stairs(np.where(keep, density + errors, np.nan), edges,
                  baseline=np.where(keep, density - errors, np.nan), color=colour,
                  alpha=0.18, fill=True)
    return total_raw


def draw_row_normalised_map(ax, counts, y_edges, canvas, min_row=MIN_DENOMINATOR,
                            tick_every=2, cmap="viridis", gate=None, trim=False):
    """A (companion, angle) map with every companion row normalised to unit sum.

    The angle axis is the first axis of `counts` and spans zero to one; the companion is
    the second and carries the variable binning, drawn on a bin index so that a first edge
    at zero is not lost to a log scale.  Rows holding fewer than `min_row` entries are
    left blank instead of drawn, since a row of two entries is noise at full contrast.
    """
    counts = np.asarray(counts, dtype=float)
    rows = counts.T
    # the population test is made on unweighted entries: a cross-section-weighted row can
    # be large and still rest on two events
    entries = np.asarray(gate if gate is not None else counts, dtype=float).T
    totals = rows.sum(axis=1, keepdims=True)
    populations = entries.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        shape = np.divide(rows, totals, out=np.full_like(rows, np.nan), where=totals > 0)
    shape[np.repeat(populations < min_row, rows.shape[1], axis=1)] = np.nan
    x_edges = np.linspace(0.0, 1.0, counts.shape[0] + 1)
    y_index = np.arange(len(y_edges), dtype=float)
    mesh = ax.pcolormesh(x_edges, y_index, np.ma.masked_invalid(shape), cmap=cmap,
                         shading="flat")
    ticks = np.arange(0, len(y_edges), tick_every, dtype=float)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{y_edges[int(t)]:g}" for t in ticks], fontsize=canvas.font)
    ax.set_ylim(0, len(y_edges) - 1)
    if trim:
        # a quantity that only exists over part of its axis -- a resolved electron pair
        # never has a mass below about a GeV -- otherwise spends most of the panel blank
        filled = np.flatnonzero(populations.ravel() >= min_row)
        if filled.size:
            ax.set_ylim(max(int(filled[0]) - 1, 0),
                        min(int(filled[-1]) + 2, len(y_edges) - 1))
    ax.set_xlim(0, 1)
    return mesh


# --------------------------------------------------------------------------------------
# j1 / j2 -- the decay angle of a muon lepton jet, mass by mass
# --------------------------------------------------------------------------------------

PLANES_MU_HIST = "mu_lj_mumu_absCosTheta_vs_mass"


def _planes_mu_sets(m_xx=500.0):
    """The (label, samples, channel) sets the muon-LJ figures compare."""
    points = planes_points("4Mu")
    sets = []
    for m_zd in (0.25, 1.2, 5.0):
        samples = points.get((m_xx, m_zd))
        if samples:
            sets.append((f"signal 4$\\mu$\n$m_{{XX}}$ = {m_xx:.0f}, $m_{{Z_d}}$ = {m_zd} GeV",
                         samples, "4mu", False))
    sets.append(("simulated background\ntwo-LJ inv. displaced CR",
                 planes_background_samples(), PLANES_CR_2LJ, True))
    sets.append(("2018 data\none $\\mu$-LJ CR", planes_data_samples(), PLANES_CR_1MU, False))
    sets.append(("2018 data\ntwo-LJ inv. displaced CR", planes_data_samples(),
                 PLANES_CR_2LJ, False))
    return sets


def fig_mulj_mass_angle(m_xx=500.0, stem="j1_mulj_mass_angle"):
    """The (dimuon mass, decay angle) plane of a muon lepton jet, signal against fakes.

    Every row of every panel is normalised on its own, so the question the figure asks is
    "at this dimuon mass, how is the angle distributed", which is the comparison that does
    not depend on the wildly different mass spectra of signal and fakes.
    """
    canvas = Canvas(2, 3, panel_w=5.6, panel_h=4.5, base_pt=9.0)
    y_edges = planes_edges(PLANES_MU_HIST, "mumu_mass")
    sets = _planes_mu_sets(m_xx)
    meshes, map_axes = [], []
    for index, (label, samples, channel, weighted) in enumerate(sets[:6]):
        ax = canvas.ax(index // 3, index % 3)
        counts = planes_counts(samples, PLANES_MU_HIST, channel, weighted=weighted)
        raw = counts if not weighted else planes_counts(samples, PLANES_MU_HIST, channel)
        if counts is None or raw is None or counts.sum() <= 0 or y_edges is None:
            canvas.hide(ax, "no input")
            continue
        meshes.append(draw_row_normalised_map(ax, counts, y_edges, canvas, gate=raw))
        map_axes.append(ax)
        canvas.labels(ax, r"$|\cos\theta^{*}|$ of the leading $\mu\mu$ pair",
                      r"$m_{\mu\mu}$ [GeV]" if index % 3 == 0 else None)
        canvas.stamp(ax, f"{label}\n{raw.sum():,.0f} $\\mu$-LJ pairs", loc="upper left")
    colour_bar = canvas.shared_colourbar(map_axes, meshes,
                                        "fraction of the mass row (99th pct scale)")
    # a single sparsely populated row, in the weighted simulation especially, otherwise
    # sets the scale for all six panels and flattens the structure everywhere else
    if meshes and colour_bar is not None:
        pooled = np.concatenate([np.asarray(mesh.get_array(), dtype=float).ravel()
                                 for mesh in meshes])
        pooled = pooled[np.isfinite(pooled)]
        if pooled.size:
            top = max(float(np.percentile(pooled, 99.0)), 1.0e-6)
            for mesh in meshes:
                mesh.set_clim(0.0, top)
    canvas.finish(data=True, lumi_fb=LUMI_FB)
    canvas.save(stem)
    return canvas



def fig_mulj_angle_overlays(m_xx=500.0, stem="j2_mulj_angle_overlays"):
    """The muon-LJ decay angle as a shape, and what a threshold on it would separate.

    Left: unit-area shapes of signal and of the three fake samples under their own
    selections.  Right: the pair of efficiencies a one-sided cut reaches, signal against
    each fake sample; what a single "separation" number hides is that the three fake
    samples do not agree with each other.
    """
    canvas = Canvas(1, 2, panel_w=6.4, panel_h=5.0, base_pt=9.0)
    points = planes_points("4Mu")
    edges = np.linspace(0.0, 1.0, 25 + 1)

    fakes = [("data, 1 $\\mu$-LJ CR", planes_data_samples(), PLANES_CR_1MU,
              False, COLOURS["data"], "o"),
             ("data, 2-LJ CR", planes_data_samples(), PLANES_CR_2LJ,
              False, "#8e7dbe", "s"),
             ("sim bkg, 2-LJ CR", planes_background_samples(),
              PLANES_CR_2LJ, True, COLOURS["background"], "^")]

    ax = canvas.ax(0, 0)
    for index, m_zd in enumerate((0.25, 1.2, 5.0)):
        samples = points.get((m_xx, m_zd))
        if not samples:
            continue
        draw_counts(ax, planes_counts(samples, PLANES_MU_HIST, "4mu").sum(axis=1), edges,
                    f"signal {m_zd}", _SLICE_COLOURS[index], canvas,
                    as_points=False, linestyle=["-", "--", "-."][index],
                    context="4mu channel")
    for offset, (label, samples, channel, weighted, colour, marker) in enumerate(fakes):
        counts = planes_counts(samples, PLANES_MU_HIST, channel, weighted=weighted)
        raw = counts if not weighted else planes_counts(samples, PLANES_MU_HIST, channel)
        if counts is None:
            continue
        draw_counts(ax, counts.sum(axis=1), edges, label, colour, canvas, marker=marker,
                    raw=raw.sum(axis=1), offset=0.006 * (offset - 1), context=channel)
    canvas.labels(ax, r"$|\cos\theta^{*}|$ of the leading $\mu\mu$ pair",
                  "$\\mu$-LJ pairs / unit angle (normalised)")
    canvas.room_above(ax, share=0.52)
    canvas.legend(ax, "upper center", ncol=2, over_data=True,
                  title=f"signal 4$\\mu$, $m_{{XX}}$ = {m_xx:.0f} GeV, four-muon channel;\n"
                        "signal curves labelled by $m_{Z_d}$ [GeV]")

    ax = canvas.ax(0, 1)
    reference = points.get((m_xx, 1.2))
    signal = (planes_counts(reference, PLANES_MU_HIST, "4mu").sum(axis=1)
              if reference else None)
    if signal is not None:
        for label, samples, channel, weighted, colour, marker in fakes:
            counts = planes_counts(samples, PLANES_MU_HIST, channel, weighted=weighted)
            if counts is None:
                continue
            thresholds, eff_s, eff_b = separation_curve(signal, counts.sum(axis=1))
            ax.plot(eff_b, eff_s, "-", color=colour, linewidth=2.0, label=label)
            marks = [int(round(t * (len(thresholds) - 1))) for t in (0.2, 0.4, 0.6, 0.8)]
            ax.plot(eff_b[marks], eff_s[marks], marker, color=colour,
                    markersize=canvas.font * 0.20, linestyle="none")
        ax.plot([0, 1], [0, 1], ":", color="#666666", linewidth=1.5,
                label="no separation")
    else:
        canvas.hide(ax, "no signal input")
    canvas.labels(ax, "fraction of fake $\\mu$-LJ pairs kept",
                  "fraction of signal $\\mu$-LJ pairs kept")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    canvas.legend(ax, "lower right", over_data=True,
                  title="cut $|\\cos\\theta^{*}| >$ threshold, markers at\n"
                        "0.2, 0.4, 0.6, 0.8; signal $m_{Z_d}$ = 1.2 GeV")

    canvas.finish(data=True, lumi_fb=LUMI_FB)
    canvas.save(stem)
    return canvas


def planes_mu_separation_table(thresholds=(0.2, 0.4, 0.6, 0.8),
                               min_m_xx=PLANES_MIN_M_XX):
    """Print what a threshold on the muon-LJ decay angle keeps and rejects.

    Signal efficiencies come from each point's own signal channel and fake efficiencies
    from the control regions, so a row is "what this threshold would do here", not "what
    it would do to one sample".
    """
    rows = []
    fakes = [("data, one mu-LJ CR", planes_data_samples(), PLANES_CR_1MU, False),
             ("data, two-LJ invDisp CR", planes_data_samples(), PLANES_CR_2LJ, False),
             ("data, 4mu invDisp VR", planes_data_samples(), PLANES_VR_4MU, False),
             ("sim bkg, two-LJ invDisp CR", planes_background_samples(), PLANES_CR_2LJ, True)]
    sets = []
    for final_state, channel in (("4Mu", "4mu"), ("2Mu2E", "2mu2e")):
        for (m_xx, m_zd), samples in planes_points(final_state, min_m_xx).items():
            sets.append((f"{final_state} {m_xx:.0f}/{m_zd}", samples, channel, False))
    sets.extend(fakes)
    header = "  ".join(f"eff(>{t:.1f})" for t in thresholds)
    print("muon-LJ |cos(theta*)|: one entry per mu-LJ with two or more muons")
    print(f"{'set':28s} {'channel':>17s} {'pairs':>10s} {'<|cos|>':>18s}   {header}")
    print("-" * 118)
    for label, samples, channel, weighted in sets:
        counts = planes_counts(samples, PLANES_MU_HIST, channel, weighted=weighted)
        if counts is None:
            continue
        raw = (counts if not weighted
               else planes_counts(samples, PLANES_MU_HIST, channel)).sum(axis=1)
        projection = counts.sum(axis=1)
        mean, error, _ = binned_mean(raw)
        cells = "  ".join(f"{kept_fraction(projection, t)['fraction']:9.4f}"
                          for t in thresholds)
        print(f"{label:28s} {planes_short(channel):>17s} {raw.sum():10.0f} "
              f"{mean:.4f} +- {error:.4f}   {cells}")
        rows.append({"label": label, "channel": channel, "counts": projection, "raw": raw,
                     "mean": mean, "error": error})
    print()
    print("signal efficiency over root fake efficiency, per mu-LJ, at |cos(theta*)| > 0.4")
    fake_rows = {row["label"]: row for row in rows if row["label"].startswith(("data", "sim"))}
    print(f"{'signal point':28s} " + "  ".join(f"{name:>26s}" for name in fake_rows))
    print("-" * 120)
    for row in rows:
        if row["label"].startswith(("data", "sim")):
            continue
        eff_s = kept_fraction(row["counts"], 0.4)["fraction"]
        cells = []
        for name, fake in fake_rows.items():
            eff_b = kept_fraction(fake["counts"], 0.4)["fraction"]
            cells.append(f"{eff_s / np.sqrt(eff_b) if eff_b > 0 else float('nan'):26.3f}")
        print(f"{row['label']:28s} " + "  ".join(cells))
    return rows


# --------------------------------------------------------------------------------------
# j3 -- the electron pair of an egamma-type lepton jet
# --------------------------------------------------------------------------------------

def planes_egm_resolved_table(min_m_xx=0.0):
    """Print how often an egamma-type lepton jet holds two reconstructed electrons.

    The denominator is every egamma-type LJ in the channel, from the electron-multiplicity
    histogram; the numerator is those with two or more, which are the only ones any
    `egm_lj_ee_*` histogram is filled from.  At a low dark-photon mass the two electrons
    of the decay merge into one supercluster plus photons, so this fraction is the first
    thing to quote about any electron-pair figure.
    """
    rows = []
    sets = []
    for (m_xx, m_zd), samples in planes_points("2Mu2E", min_m_xx).items():
        sets.append((f"2Mu2E {m_xx:.0f}/{m_zd}", samples, "2mu2e", False))
    sets += [("data, one egm-LJ CR", planes_data_samples(), PLANES_CR_1EGM, False),
             ("data, two-LJ invDisp CR", planes_data_samples(), PLANES_CR_2LJ, False),
             ("data, 2mu2e invDisp VR", planes_data_samples(), PLANES_VR_2MU2E, False),
             ("sim bkg, two-LJ invDisp CR", planes_background_samples(), PLANES_CR_2LJ, True)]
    print("egamma-type LJs holding two or more reconstructed electrons")
    print(f"{'set':28s} {'channel':>17s} {'egm LJs':>10s} {'with 2+ e':>10s} "
          f"{'fraction':>10s} {'ee entries':>11s}")
    print("-" * 90)
    for label, samples, channel, _weighted in sets:
        # counts, not yields: the resolved fraction is a per-object reconstruction
        # efficiency, so it is measured on unweighted entries in every sample
        multiplicity = planes_counts(samples, "egm_lj_electronN", channel)
        pairs = planes_counts(samples, "egm_lj_ee_absCosTheta", channel)
        if multiplicity is None:
            continue
        total = float(multiplicity.sum())
        resolved = float(multiplicity[2:].sum())
        fraction = resolved / total if total > 0 else float("nan")
        print(f"{label:28s} {planes_short(channel):>17s} {total:10.0f} "
              f"{resolved:10.0f} {fraction:10.5f} "
              f"{0.0 if pairs is None else pairs.sum():11.0f}")
        rows.append({"label": label, "channel": channel, "n_egm": total,
                     "n_resolved": resolved, "fraction": fraction})
    print()
    print("A pair histogram of an egamma-type LJ is filled only by the resolved fraction")
    print("above, so every electron-pair number below is a statement about that subset and")
    print("not about the egamma lepton jets of the channel.")
    return rows



def fig_egm_ee(stem="j3_egm_ee"):
    """The resolved electron pair of an egamma-type lepton jet: how often, and how it looks."""
    canvas = Canvas(2, 3, panel_w=5.6, panel_h=4.5, base_pt=9.0)
    points = planes_points("2Mu2E", 0.0)
    masses = sorted({m_zd for _, m_zd in points})

    # (0,0) how often the pair is resolved at all, against the dark-photon mass
    ax = canvas.ax(0, 0)
    for index, m_xx in enumerate(sorted({m_xx for m_xx, _ in points})):
        values, drawn = [], []
        for m_zd in masses:
            samples = points.get((m_xx, m_zd))
            if not samples:
                continue
            multiplicity = planes_counts(samples, "egm_lj_electronN", "2mu2e")
            if multiplicity is None or multiplicity.sum() <= 0:
                continue
            values.append(float(multiplicity[2:].sum()) / float(multiplicity.sum()))
            drawn.append(m_zd)
        if drawn:
            ax.plot(drawn, values, "o-", color=_SLICE_COLOURS[index],
                    markersize=canvas.font * 0.22, linewidth=2.0,
                    label=f"{m_xx:.0f}")
    measured = []
    for label, samples, channel, colour, style in [
            ("1 $e\\gamma$-LJ CR", planes_data_samples(), PLANES_CR_1EGM,
             COLOURS["data"], "--"),
            ("2-LJ CR", planes_data_samples(), PLANES_CR_2LJ, "#8e7dbe", ":")]:
        multiplicity = planes_counts(samples, "egm_lj_electronN", channel)
        if multiplicity is None or multiplicity.sum() <= 0:
            continue
        fraction = float(multiplicity[2:].sum()) / float(multiplicity.sum())
        # a key entry for each of these sat on the rising signal curves, and the number
        # itself is more use than a line style: it goes in a note in the empty corner
        ax.axhline(fraction, color=colour, linestyle=style, linewidth=2.0)
        measured.append(f"{label}: {100 * fraction:.2f}%")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(masses)
    ax.set_xticklabels([f"{m:g}" for m in masses])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_minor_locator(NullLocator())
    canvas.labels(ax, r"$m_{Z_d}$ [GeV]", "fraction with two electrons")
    canvas.legend(ax, "lower right", over_data=True,
                  title="2$\\mu$2e; $m_{XX}$ [GeV]")
    if measured:
        canvas.stamp(ax, "2018 data (horizontal lines)\n" + "\n".join(measured),
                     loc="upper left")

    # (0,1) and (0,2) the (mass, angle) plane where there are enough pairs to draw one
    y_edges = planes_edges("egm_lj_ee_absCosTheta_vs_mass", "ee_mass")
    plane_sets = [(_first_key(points, (500.0, 5.0), (1000.0, 5.0)), "2mu2e"),
                  (None, PLANES_CR_2LJ)]
    meshes, map_axes = [], []
    for column, (key, channel) in enumerate(plane_sets, start=1):
        ax = canvas.ax(0, column)
        if key is None:
            samples, label = planes_data_samples(), "2018 data, 2-LJ CR"
        else:
            samples = points.get(key, [])
            label = f"signal {key[0]:.0f} / {key[1]} GeV"
        counts = planes_counts(samples, "egm_lj_ee_absCosTheta_vs_mass", channel)
        if counts is None or counts.sum() < MIN_SHAPE_ENTRIES or y_edges is None:
            canvas.hide(ax, "too few resolved\nelectron pairs")
            continue
        meshes.append(draw_row_normalised_map(ax, counts, y_edges, canvas, trim=True))
        map_axes.append(ax)
        canvas.labels(ax, r"$|\cos\theta^{*}|$ of the $ee$ pair",
                      r"$m_{ee}$ [GeV]" if column == 1 else None)
        canvas.stamp(ax, f"{label}, {counts.sum():,.0f} pairs", loc="upper left")
    if map_axes:
        canvas.shared_colourbar(map_axes, meshes, "fraction of the mass row")
        canvas.keep_label_position(*map_axes)

    # (1,*) the three one-dimensional shapes, keyed once
    shape_sets = []
    for key in ((500.0, 5.0), (1000.0, 5.0), (500.0, 1.2)):
        if key in points:
            shape_sets.append((f"signal {key[0]:.0f} / {key[1]}", points[key], "2mu2e",
                               False, _SLICE_COLOURS[len(shape_sets)], "-", None))
    shape_sets.append(("data, 1 $e\\gamma$-LJ CR", planes_data_samples(),
                       PLANES_CR_1EGM, False, COLOURS["data"], None, "o"))
    shape_sets.append(("data, 2-LJ CR", planes_data_samples(),
                      PLANES_CR_2LJ, False, "#8e7dbe", None, "s"))
    shape_sets.append(("sim bkg, 2-LJ CR", planes_background_samples(),
                       PLANES_CR_2LJ, True, COLOURS["background"], None, "^"))
    panels = [("egm_lj_ee_absCosTheta", 5, 0.0, 1.0,
               r"$|\cos\theta^{*}|$ of the $ee$ pair", "pairs / unit (normalised)"),
              ("egm_lj_ee_ptRatio", 5, 0.0, 1.0,
               r"$p_T^{sub}/p_T^{lead}$ of the $ee$ pair", None),
              ("egm_lj_ee_dR", 5, 0.0, 0.5, r"$\Delta R(e^{lead}, e^{sub})$", None)]
    for column, (name, factor, low, high, xlabel, ylabel) in enumerate(panels):
        ax = canvas.ax(1, column)
        for index, (label, samples, channel, weighted, colour, style, marker) in \
                enumerate(shape_sets):
            counts = planes_counts(samples, name, channel, weighted=weighted)
            if counts is None:
                continue
            raw = counts if not weighted else planes_counts(samples, name, channel)
            rebinned = rebin1d(counts, factor)
            rebinned_raw = rebin1d(raw, factor)
            edges = np.linspace(low, high, len(rebinned) + 1)
            draw_counts(ax, rebinned, edges, label, colour, canvas,
                        as_points=marker is not None, marker=marker or "o",
                        linestyle=style or "-", raw=rebinned_raw,
                        offset=0.004 * (index - 2), context=channel)
        canvas.labels(ax, xlabel, ylabel)
        if column == 0:
            canvas.room_above(ax, share=0.44)
            legend = canvas.legend(ax, "upper center", ncol=2, over_data=True,
                                   title="signal $m_{XX}$ / $m_{Z_d}$",
                                   fontsize=canvas.font * 0.90)
            # the signal shapes rise towards one, so there is no empty corner here and the
            # room has to be measured from the key rather than guessed at
            canvas.clear_legends([(ax, legend)], most=6.0)
        else:
            # the low-mass signal puts nearly every pair in one bin while the fake
            # samples are spread over all of them, a range no linear scale shows at once
            canvas.log_axis(ax, "y")
            top = float(ax.dataLim.y1)
            if np.isfinite(top) and top > 0:
                ax.set_ylim(max(top * 1.0e-3, 1.0e-4), top * 4.0)
    canvas.finish(data=True, lumi_fb=LUMI_FB)
    canvas.save(stem)
    return canvas


def planes_egm_shape_table(min_m_xx=PLANES_MIN_M_XX):
    """Print the mean angle, momentum ratio and opening angle of resolved electron pairs."""
    sets = []
    for (m_xx, m_zd), samples in planes_points("2Mu2E", min_m_xx).items():
        sets.append((f"2Mu2E {m_xx:.0f}/{m_zd}", samples, "2mu2e", False))
    sets += [("data, one egm-LJ CR", planes_data_samples(), PLANES_CR_1EGM, False),
             ("data, two-LJ invDisp CR", planes_data_samples(), PLANES_CR_2LJ, False),
             ("sim bkg, two-LJ invDisp CR", planes_background_samples(), PLANES_CR_2LJ, True)]
    print("resolved electron pairs of egamma-type LJs")
    print(f"{'set':28s} {'channel':>17s} {'pairs':>8s} {'<|cos|>':>18s} "
          f"{'<ptRatio>':>18s} {'<dR>':>18s}")
    print("-" * 112)
    rows = []
    for label, samples, channel, weighted in sets:
        counts = planes_counts(samples, "egm_lj_ee_absCosTheta", channel, weighted=weighted)
        raw = (counts if not weighted
               else planes_counts(samples, "egm_lj_ee_absCosTheta", channel))
        if counts is None or raw is None or raw.sum() < 1:
            continue
        ratio = planes_counts(samples, "egm_lj_ee_ptRatio", channel, weighted=weighted)
        radius = planes_counts(samples, "egm_lj_ee_dR", channel, weighted=weighted)
        # the value carries the cross-section weights, the uncertainty the entries behind it
        angle_mean = binned_mean(counts)[0]
        ratio_mean = binned_mean(ratio, 0.0, 1.0)[0]
        radius_mean = binned_mean(radius, 0.0, 0.5)[0]
        total = float(raw.sum())
        angle_error = binned_mean(raw)[1]
        ratio_error = binned_mean(planes_counts(samples, "egm_lj_ee_ptRatio", channel),
                                  0.0, 1.0)[1]
        radius_error = binned_mean(planes_counts(samples, "egm_lj_ee_dR", channel),
                                   0.0, 0.5)[1]
        flag = "" if total >= MIN_SHAPE_ENTRIES else "   (below the shape floor)"
        print(f"{label:28s} {planes_short(channel):>17s} {total:8.0f} "
              f"{angle_mean:.4f} +- {angle_error:.4f} "
              f"{ratio_mean:.4f} +- {ratio_error:.4f} "
              f"{radius_mean:.4f} +- {radius_error:.4f}{flag}")
        rows.append({"label": label, "n": total, "mean": angle_mean, "error": angle_error,
                     "ptRatio": ratio_mean, "dR": radius_mean})
    return rows


# --------------------------------------------------------------------------------------
# j4 -- the two helicity angles of the lepton-jet pair
# --------------------------------------------------------------------------------------

PLANES_PAIR_PLANE = "lj_pair_cosThetaA0_vs_cosThetaA1"
PLANES_GEN_PLANE = "gen_pair_cosThetaA0_vs_cosThetaA1"


def fig_pair_helicity_plane(point=(500.0, 1.2), stem="j4_pair_helicity_plane"):
    """The two signed helicity angles of an event, as a ratio to the product of marginals.

    The same construction as section i, on the signed angle of the lepton-jet pair rather
    than the absolute angle of the muon pair inside one lepton jet: a plane of ones is
    independence, and a cell above one is a combination that occurs more often than two
    independent draws would give.
    """
    canvas = Canvas(1, 2, panel_w=5.7, panel_h=5.0, base_pt=9.0)
    samples = planes_points("4Mu", 0.0).get(point, [])
    if not samples:
        note(f"no plane samples for m_XX = {point[0]} GeV, m_Zd = {point[1]} GeV")
    sets = [(0, planes_data_samples(), PLANES_CR_2LJ,
             "2018 data, two-LJ inv. displaced CR"),
            (1, samples, "4mu",
             f"signal 4$\\mu$, $m_{{XX}}$ = {point[0]:.0f} GeV,\n"
             f"$m_{{Z_d}}$ = {point[1]} GeV, four-muon channel")]
    # five bins a side, the binning the printed table also tests independence on: on the
    # 8,400 control-region events that is several hundred entries a cell, where ten bins a
    # side would put eighty in a cell and fill the map with noise
    edges = np.linspace(-1.0, 1.0, 6)
    meshes = []
    for column, group, channel, label in sets:
        ax = canvas.ax(0, column)
        fine = planes_counts(group, PLANES_PAIR_PLANE, channel)
        if fine is None or fine.sum() < MIN_SHAPE_ENTRIES:
            canvas.hide(ax, "no input")
            continue
        result = plane_independence(rebin2d(fine, 4))
        correlation = plane_correlation(fine, replicas=400)
        meshes.append(_draw_ratio_map(ax, result, canvas, edges=edges, headroom=1.42))
        canvas.labels(ax, r"$\cos\theta_A$, leading LJ",
                      r"$\cos\theta_A$, subleading LJ" if column == 0 else None)
        ax.text(0.02, 0.985,
                f"{label}\n"
                f"$r$ = {correlation['r']:+.3f} $\\pm$ {correlation['sigma']:.3f}, "
                f"{result['n']:,.0f} events\n"
                f"independence $\\chi^2$/dof = {result['chi2']:.0f}/{result['ndof']}",
                transform=ax.transAxes, ha="left", va="top", fontsize=canvas.font,
                linespacing=1.35)
    if meshes:
        colour_bar = canvas.fig.colorbar(meshes[-1], ax=list(canvas.axes[0]), pad=0.02,
                                         fraction=0.05)
        colour_bar.ax.tick_params(labelsize=canvas.tick_font)
        colour_bar.set_label("observed / product of marginals", fontsize=canvas.font)
    canvas.finish(data=True, lumi_fb=LUMI_FB)
    canvas.save(stem)
    return canvas


def planes_pair_correlation_table(min_m_xx=PLANES_MIN_M_XX):
    """Print the correlation of the two lepton jets' helicity angles, event by event."""
    sets = [("2018 data", planes_data_samples(), PLANES_CR_2LJ),
            ("2018 data", planes_data_samples(), PLANES_VR_4MU),
            ("2018 data", planes_data_samples(), PLANES_VR_2MU2E),
            ("simulated background", planes_background_samples(), PLANES_CR_2LJ)]
    for final_state, channel in (("4Mu", "4mu"), ("2Mu2E", "2mu2e")):
        for (m_xx, m_zd), samples in planes_points(final_state, min_m_xx).items():
            sets.append((f"{final_state} {m_xx:.0f}/{m_zd}", samples, channel))
    print("linear correlation of cos(theta_A) between the two leading lepton jets")
    print(f"{'set':24s} {'channel':>17s} {'events':>9s} {'r':>9s} {'sigma(r)':>9s} "
          f"{'chi2/ndof':>13s} {'p':>9s}")
    print("-" * 95)
    rows = []
    for label, samples, channel in sets:
        fine = planes_counts(samples, PLANES_PAIR_PLANE, channel)
        if fine is None or fine.sum() < MIN_DENOMINATOR:
            print(f"{label:24s} {planes_short(channel):>17s} "
                  f"{0.0 if fine is None else fine.sum():9.0f}   too few events")
            continue
        correlation = plane_correlation(fine, replicas=400)
        independence = plane_independence(rebin2d(fine, 4))
        print(f"{label:24s} {planes_short(channel):>17s} {correlation['n']:9.0f} "
              f"{correlation['r']:+9.4f} {correlation['sigma']:9.4f} "
              f"{independence['chi2']:8.1f}/{independence['ndof']:<4d} "
              f"{independence['p_value']:9.3f}")
        rows.append({"label": label, "channel": channel, **correlation,
                     "p_value": independence["p_value"]})
    print()
    print("the same at generator level, for the signal points that carry a gen record")
    for final_state, channel in (("4Mu", "4mu"), ("2Mu2E", "2mu2e")):
        for (m_xx, m_zd), samples in planes_points(final_state, min_m_xx).items():
            fine = planes_counts(samples, PLANES_GEN_PLANE, channel)
            if fine is None or fine.sum() < MIN_DENOMINATOR:
                continue
            correlation = plane_correlation(fine, replicas=400)
            print(f"  gen {final_state} {m_xx:.0f}/{m_zd}: r = {correlation['r']:+.4f} "
                  f"+- {correlation['sigma']:.4f} on {correlation['n']:,.0f} events")
    return rows


# --------------------------------------------------------------------------------------
# j5 -- the angle between the two decay planes
# --------------------------------------------------------------------------------------

def _phi_reference(edges, coefficient=PLANES_A_CP_EVEN, harmonic=2):
    """The unit-area density 1 + A cos(n x) evaluated at the centres of `edges`."""
    centres = 0.5 * (np.asarray(edges)[:-1] + np.asarray(edges)[1:])
    span = float(edges[-1] - edges[0])
    return centres, (1.0 + coefficient * np.cos(harmonic * centres)) / span



def fig_pair_planes_phi(points=((500.0, 1.2), (500.0, 5.0)), stem="j5_pair_phi"):
    """The azimuthal angle between the two lepton-jet decay planes, reconstructed and generated.

    The generator-level curve comes from the same files and the same channel as the
    reconstructed one, so it carries the production filter and the selection with it; the
    analytic curve is drawn only as the unsculpted reference it is.

    Every distribution here varies about a non-zero mean by about a quarter of that mean,
    so the panels take explicit room for their keys rather than scaling from zero: a key
    anchored to the bottom of an autoscaled panel sits exactly on the minima of
    1 + A cos(2 Phi), which are the half of the curve that carries the coefficient.
    """
    canvas = Canvas(2, 2, panel_w=6.0, panel_h=4.8, base_pt=9.0)
    grid = planes_points("4Mu", 0.0)
    phi_edges = np.linspace(-np.pi, np.pi, 13)
    angle_edges = np.linspace(0.0, np.pi / 2, 10)

    # (0,0) signal, reconstructed against generated
    ax = canvas.ax(0, 0)
    for index, key in enumerate(points):
        samples = grid.get(key, [])
        if not samples:
            continue
        colour = _SLICE_COLOURS[index]
        draw_counts(ax, rebin1d(planes_counts(samples, "gen_pair_Phi", "4mu"), 2),
                    phi_edges, f"gen, {key[1]}", colour, canvas,
                    as_points=False, linestyle="-", context="4mu channel")
        draw_counts(ax, rebin1d(planes_counts(samples, "lj_pair_Phi", "4mu"), 2),
                    phi_edges, f"reco, {key[1]}", colour, canvas,
                    marker=["o", "s"][index % 2], offset=0.004 * (index - 0.5),
                    context="4mu channel")
    centres, reference = _phi_reference(phi_edges)
    ax.plot(centres, reference, ":", color="#666666", linewidth=2.0,
            label=r"$1 + \frac{1}{4}\cos 2\Phi$")
    canvas.labels(ax, r"$\Phi$ between the two LJ decay planes [rad]",
                  "events / rad (normalised)")
    canvas.room_above(ax, share=0.44)
    canvas.legend(ax, "upper center", ncol=3, over_data=True,
                  title=f"signal 4$\\mu$, $m_{{XX}}$ = {points[0][0]:.0f} GeV, four-muon "
                        "channel;\ncurves labelled by $m_{Z_d}$ [GeV]")

    # (0,1) fakes
    ax = canvas.ax(0, 1)
    for index, (label, samples, channel, weighted, colour, marker) in enumerate([
            ("data, 2-LJ CR", planes_data_samples(), PLANES_CR_2LJ,
             False, COLOURS["data"], "o"),
            ("data, 4$\\mu$ VR", planes_data_samples(), PLANES_VR_4MU,
             False, "#8e7dbe", "s"),
            ("sim bkg, 2-LJ CR", planes_background_samples(),
             PLANES_CR_2LJ, True, COLOURS["background"], "^")]):
        counts = planes_counts(samples, "lj_pair_Phi", channel, weighted=weighted)
        if counts is None:
            continue
        raw = counts if not weighted else planes_counts(samples, "lj_pair_Phi", channel)
        draw_counts(ax, rebin1d(counts, 2), phi_edges, label, colour, canvas,
                    marker=marker, raw=rebin1d(raw, 2), offset=0.005 * (index - 1),
                    context=channel)
    ax.plot(centres, reference, ":", color="#666666", linewidth=2.0,
            label=r"$1 + \frac{1}{4}\cos 2\Phi$")
    canvas.labels(ax, r"$\Phi$ between the two LJ decay planes [rad]", None)
    canvas.room_above(ax, share=0.46)
    canvas.legend(ax, "upper center", ncol=2, over_data=True, title="fake LJ pairs")

    # (1,0) the unoriented angle between the planes
    ax = canvas.ax(1, 0)
    for index, key in enumerate(points):
        samples = grid.get(key, [])
        if not samples:
            continue
        colour = _SLICE_COLOURS[index]
        draw_counts(ax, rebin1d(planes_counts(samples, "gen_pair_planeAngle", "4mu"), 2),
                    angle_edges, f"gen, {key[1]}", colour, canvas,
                    as_points=False, context="4mu channel")
        draw_counts(ax, rebin1d(planes_counts(samples, "lj_pair_planeAngle", "4mu"), 2),
                    angle_edges, f"reco, {key[1]}", colour, canvas,
                    marker=["o", "s"][index % 2], offset=0.004 * (index - 0.5),
                    context="4mu channel")
    counts = planes_counts(planes_data_samples(), "lj_pair_planeAngle", PLANES_CR_2LJ)
    draw_counts(ax, rebin1d(counts, 2), angle_edges, "data, 2-LJ CR",
                COLOURS["data"], canvas, marker="v", context=PLANES_CR_2LJ)
    centres_angle, reference_angle = _phi_reference(angle_edges)
    ax.plot(centres_angle, reference_angle, ":", color="#666666", linewidth=2.0,
            label=r"$1 + \frac{1}{4}\cos 2\alpha$")
    canvas.labels(ax, r"$\alpha$ between the two LJ decay planes [rad]",
                  "events / rad (normalised)")
    canvas.room_above(ax, share=0.44)
    canvas.legend(ax, "upper center", ncol=3, over_data=True,
                  title="signal curves labelled by $m_{Z_d}$ [GeV]")

    # (1,1) the fitted coefficient, point by point
    ax = canvas.ax(1, 1)
    masses = [0.25, 1.2, 5.0]
    series = []
    for final_state, channel in (("4Mu", "4mu"), ("2Mu2E", "2mu2e")):
        for m_xx in sorted({key[0] for key in planes_points(final_state)}):
            series.append((final_state, channel, m_xx))
    for index, (final_state, channel, m_xx) in enumerate(series):
        colour = _SLICE_COLOURS[index % len(_SLICE_COLOURS)]
        grid_fs = planes_points(final_state, 0.0)
        # gen and reco are told apart by the marker fill and the line style, so the key
        # carries one entry per mass point instead of two
        for name, hist_name, style in (("generated", "gen_pair_Phi", "--"),
                                       ("reconstructed", "lj_pair_Phi", "-")):
            x_values, y_values, y_errors = [], [], []
            for m_zd in masses:
                samples = grid_fs.get((m_xx, m_zd))
                if not samples:
                    continue
                counts = planes_counts(samples, hist_name, channel)
                if counts is None or counts.sum() < MIN_SHAPE_ENTRIES:
                    continue
                moment = harmonic_moment(counts, -np.pi, np.pi)
                x_values.append(m_zd)
                y_values.append(moment["coefficient"])
                y_errors.append(moment["sigma_coefficient"])
            if not x_values:
                continue
            ax.errorbar(x_values, y_values, yerr=y_errors, fmt="o", linestyle=style,
                        color=colour, markersize=canvas.font * 0.20, elinewidth=1.2,
                        capsize=0, linewidth=1.8,
                        markerfacecolor="white" if name == "generated" else colour,
                        markeredgecolor=colour, markeredgewidth=1.5,
                        label=(f"{final_state} {m_xx:.0f}"
                               if name == "reconstructed" else None))
    fakes = planes_counts(planes_data_samples(), "lj_pair_Phi", PLANES_CR_2LJ)
    if fakes is not None and fakes.sum() >= MIN_SHAPE_ENTRIES:
        moment = harmonic_moment(fakes, -np.pi, np.pi)
        ax.axhspan(moment["coefficient"] - moment["sigma_coefficient"],
                   moment["coefficient"] + moment["sigma_coefficient"],
                   color=COLOURS["data"], alpha=0.25)
        canvas.stamp(ax, "2018 data, 2-LJ CR (band):\n"
                         f"$A$ = {moment['coefficient']:+.3f} $\\pm$ "
                         f"{moment['sigma_coefficient']:.3f}", loc="lower left")
    ax.axhline(PLANES_A_CP_EVEN, color="#666666", linestyle=":", linewidth=2.0)
    ax.axhline(-PLANES_A_CP_EVEN, color="#b03060", linestyle=":", linewidth=2.0)
    style_handles = [
        Line2D([], [], color="#444444", linestyle="-", marker="o",
               markerfacecolor="#444444", label="reco (filled)"),
        Line2D([], [], color="#444444", linestyle="--", marker="o",
               markerfacecolor="white", label="gen (open)"),
        Line2D([], [], color="#666666", linestyle=":", label="$0^{+}$: $A$ = +1/4"),
        Line2D([], [], color="#b03060", linestyle=":", label="$0^{-}$: $A$ = -1/4")]
    ax.set_xscale("log")
    ax.set_xticks(masses)
    ax.set_xticklabels([f"{m:g}" for m in masses])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_minor_locator(NullLocator())
    canvas.labels(ax, r"$m_{Z_d}$ [GeV]", r"$A$ in $1 + A\cos 2\Phi$")
    handles, labels = ax.get_legend_handles_labels()
    canvas.room_above(ax, share=0.55)
    canvas.legend(ax, "upper left", handles=handles + style_handles,
                  labels=labels + [h.get_label() for h in style_handles],
                  ncol=2, over_data=True,
                  title="$m_{XX}$ [GeV] by final state", fontsize=canvas.font * 0.90)

    canvas.finish(data=True, lumi_fb=LUMI_FB)
    canvas.save(stem)
    return canvas


def planes_phi_table(min_m_xx=0.0):
    """Print the fitted cos(2 Phi) coefficient, generated and reconstructed, per point."""
    print("A in 1 + A cos(2 Phi), as twice the mean of cos(2 Phi) over the distribution")
    print("gen and reco are the same files and the same channel; for the four-muon points")
    print("they are also the same events, and for 2mu2e they are not, because a")
    print("reconstructed pair needs a resolved electron pair the generated one does not.")
    print(f"{'set':22s} {'channel':>10s} {'N reco':>9s} {'A reco':>18s} {'N gen':>9s} "
          f"{'A gen':>18s} {'reco/gen':>9s}")
    print("-" * 100)
    rows = []
    for final_state, channel in (("4Mu", "4mu"), ("2Mu2E", "2mu2e")):
        for (m_xx, m_zd), samples in planes_points(final_state, min_m_xx).items():
            reco = planes_counts(samples, "lj_pair_Phi", channel)
            gen = planes_counts(samples, "gen_pair_Phi", channel)
            if reco is None or reco.sum() < 1:
                continue
            reco_moment = harmonic_moment(reco, -np.pi, np.pi)
            gen_moment = (harmonic_moment(gen, -np.pi, np.pi) if gen is not None
                          and gen.sum() > 0 else None)
            share = (reco_moment["coefficient"] / gen_moment["coefficient"]
                     if gen_moment and abs(gen_moment["coefficient"]) > 0.02 else float("nan"))
            gen_text = ("        no gen record" if gen_moment is None else
                        f"{gen_moment['coefficient']:+.3f} +- "
                        f"{gen_moment['sigma_coefficient']:.3f}")
            print(f"{final_state + ' ' + f'{m_xx:.0f}/{m_zd}':22s} {channel:>10s} "
                  f"{reco_moment['n']:9.0f} "
                  f"{reco_moment['coefficient']:+.3f} +- "
                  f"{reco_moment['sigma_coefficient']:.3f} "
                  f"{0 if gen_moment is None else gen_moment['n']:9.0f} {gen_text:>18s} "
                  f"{share:9.2f}")
            rows.append({"final_state": final_state, "m_xx": m_xx, "m_zd": m_zd,
                         "reco": reco_moment, "gen": gen_moment})
    print()
    print(f"{'fake sample':22s} {'channel':>17s} {'N':>9s} {'A (unweighted)':>18s} "
          f"{'A (xs weighted)':>16s}")
    print("-" * 86)
    for label, samples, channel, weighted in [
            ("2018 data", planes_data_samples(), PLANES_CR_2LJ, False),
            ("2018 data", planes_data_samples(), PLANES_VR_4MU, False),
            ("2018 data", planes_data_samples(), PLANES_VR_2MU2E, False),
            ("simulated background", planes_background_samples(), PLANES_CR_2LJ, True)]:
        counts = planes_counts(samples, "lj_pair_Phi", channel, weighted=weighted)
        if counts is None or counts.sum() < 1:
            continue
        raw = counts if not weighted else planes_counts(samples, "lj_pair_Phi", channel)
        # the error belongs to the unweighted entry count; the weighted value is the shape
        # a background sum actually predicts, so both are printed and neither stands alone
        moment = harmonic_moment(raw, -np.pi, np.pi)
        weighted_moment = harmonic_moment(counts, -np.pi, np.pi)
        print(f"{label:22s} {planes_short(channel):>17s} {moment['n']:9.0f} "
              f"{moment['coefficient']:+.3f} +- {moment['sigma_coefficient']:.3f} "
              f"{weighted_moment['coefficient']:+16.3f}")
        rows.append({"final_state": label, "channel": channel, "reco": moment, "gen": None})
    return rows


# --------------------------------------------------------------------------------------
# j6 -- transverse balance and missing momentum
# --------------------------------------------------------------------------------------

PLANES_BALANCE = [("lj_pair_ptAsym", 2, 0.0, 1.0,
                   r"$|p_T^{0} - p_T^{1}| / (p_T^{0} + p_T^{1})$"),
                  ("lj_pair_ptVecSum", 4, 0.0, 400.0,
                   r"$p_T$ of the two LJs' vector sum [GeV]")]



def _planes_balance_sets():
    """The (label, samples, channel, weighted, colour, marker) sets the balance figure draws.

    Labels are short because a five-entry key at nine printed points has to fit inside a
    5.6-inch panel; the section's prose and the printed table carry the full region names.
    """
    grid = planes_points("4Mu", 0.0)
    sets = []
    for index, key in enumerate(((500.0, 1.2), (1000.0, 1.2))):
        if key in grid:
            sets.append((f"signal {key[0]:.0f} / {key[1]}", grid[key], "4mu", False,
                         _SLICE_COLOURS[index], None))
    sets.append(("data, 2-LJ CR", planes_data_samples(), PLANES_CR_2LJ, False,
                 COLOURS["data"], "o"))
    sets.append(("data, 4$\\mu$ VR", planes_data_samples(), PLANES_VR_4MU, False,
                 "#8e7dbe", "s"))
    sets.append(("sim bkg, 2-LJ CR", planes_background_samples(), PLANES_CR_2LJ, True,
                 COLOURS["background"], "^"))
    return sets




def fig_pair_balance(stem="j6_pair_balance"):
    """Whether the two lepton jets balance, and what the missing momentum adds.

    One key, in the first panel, for the five samples every panel draws in the same
    colours; the missing-momentum panel adds a second key for its two line styles only.
    Six panels each carrying the same five labels would leave no room for the curves, and
    a key wider than a 2.5-inch printed panel is drawn straight over its neighbour.
    """
    canvas = Canvas(2, 3, panel_w=5.6, panel_h=4.5, base_pt=9.0)
    sets = _planes_balance_sets()

    for column, (name, factor, low, high, xlabel) in enumerate(PLANES_BALANCE):
        ax = canvas.ax(0, column)
        for index, (label, samples, channel, weighted, colour, marker) in enumerate(sets):
            counts = planes_counts(samples, name, channel, weighted=weighted)
            if counts is None:
                continue
            raw = counts if not weighted else planes_counts(samples, name, channel)
            rebinned, rebinned_raw = rebin1d(counts, factor), rebin1d(raw, factor)
            edges = np.linspace(low, high, len(rebinned) + 1)
            draw_counts(ax, rebinned, edges, label, colour, canvas,
                        as_points=marker is not None, marker=marker or "o",
                        raw=rebinned_raw, offset=0.004 * (index - 2), context=channel)
        canvas.labels(ax, xlabel,
                      "LJ pairs / unit (normalised)" if column == 0 else None)
        if column == 0:
            canvas.room_above(ax, share=0.44)
            canvas.legend(ax, "upper right", ncol=1, over_data=True,
                          title="signal $m_{XX}$ / $m_{Z_d}$")
        else:
            canvas.headroom(ax, factor=1.2)

    # (0,2) raw and DSA-corrected missing momentum, two samples, two line styles
    ax = canvas.ax(0, 2)
    for label, samples, channel, weighted, colour, _marker in (sets[0], sets[2]):
        for name, dash in (("lj_pair_met", "-"), ("lj_pair_metCorr", "--")):
            counts = planes_counts(samples, name, channel, weighted=weighted)
            if counts is None:
                continue
            raw = counts if not weighted else planes_counts(samples, name, channel)
            rebinned, rebinned_raw = rebin1d(counts, 4), rebin1d(raw, 4)
            edges = np.linspace(0.0, 400.0, len(rebinned) + 1)
            draw_counts(ax, rebinned, edges, "_nolegend_", colour, canvas, as_points=False,
                        linestyle=dash, raw=rebinned_raw, context=channel)
    canvas.labels(ax, r"$p_T^{miss}$ [GeV]", None)
    canvas.headroom(ax, factor=1.6)
    canvas.legend(ax, "upper right", over_data=True,
                  handles=[Line2D([], [], color=sets[0][4], linestyle="-"),
                           Line2D([], [], color=sets[2][4], linestyle="-"),
                           Line2D([], [], color="#444444", linestyle="-"),
                           Line2D([], [], color="#444444", linestyle="--")],
                  labels=[sets[0][0], sets[2][0], "particle-flow", "DSA corrected"],
                  title="colour = sample, style = MET")

    # (1,0) and (1,1) the plane of the two
    meshes, map_axes = [], []
    joint = "lj_pair_metCorr_vs_ptVecSum"
    lead = sets[0] if sets else ("signal", [], "4mu", False, COLOURS["signal"], None)
    for column, (label, samples, channel, weighted) in enumerate([
            (lead[0], lead[1], lead[2], lead[3]),
            ("data, 2-LJ CR", planes_data_samples(), PLANES_CR_2LJ, False)]):
        ax = canvas.ax(1, column)
        counts = planes_counts(samples, joint, channel, weighted=weighted)
        if counts is None or counts.sum() < MIN_SHAPE_ENTRIES:
            canvas.hide(ax, "no input")
            continue
        edges = np.linspace(0.0, 400.0, counts.shape[0] + 1)
        meshes.append(ax.pcolormesh(edges, edges, (counts / counts.sum()).T,
                                    cmap="viridis", shading="flat"))
        map_axes.append(ax)
        canvas.labels(ax, r"DSA-corrected $p_T^{miss}$ [GeV]",
                      r"vector sum $p_T$ [GeV]" if column == 0 else None)
        canvas.stamp(ax, f"{label}\n{counts.sum():,.0f} pairs", loc="upper right")
    if map_axes:
        canvas.shared_colourbar(map_axes, meshes, "fraction of LJ pairs")
        canvas.keep_label_position(*map_axes)

    # (1,2) what a joint requirement on both would keep
    ax = canvas.ax(1, 2)
    for label, samples, channel, weighted, colour, marker in sets:
        counts = planes_counts(samples, joint, channel, weighted=weighted)
        if counts is None or counts.sum() <= 0:
            continue
        edges = np.linspace(0.0, 400.0, counts.shape[0] + 1)
        cumulative = np.array([counts[:i, :i].sum() for i in range(counts.shape[0] + 1)])
        ax.plot(edges, cumulative / counts.sum(), "-" if marker is None else "--",
                color=colour, linewidth=2.0, marker=marker,
                markersize=canvas.font * 0.16, markevery=2)
    canvas.labels(ax, "threshold $X$ [GeV]", "fraction with both below $X$")
    ax.set_xlim(0, 400)
    ax.set_ylim(0, 1.42)
    canvas.stamp(ax, "DSA-corrected $p_T^{miss}$ and\n"
                     "the vector sum both below $X$;\n"
                     "colours as in the first panel",
                 loc="upper left")

    canvas.finish(data=True, lumi_fb=LUMI_FB)
    canvas.save(stem)
    return canvas


def planes_balance_table(min_m_xx=PLANES_MIN_M_XX):
    """Print the balance and missing-momentum numbers, signal against fakes."""
    sets = []
    for final_state, channel in (("4Mu", "4mu"), ("2Mu2E", "2mu2e")):
        for (m_xx, m_zd), samples in planes_points(final_state, min_m_xx).items():
            sets.append((f"{final_state} {m_xx:.0f}/{m_zd}", samples, channel, False))
    sets += [("data, two-LJ invDisp CR", planes_data_samples(), PLANES_CR_2LJ, False),
             ("data, 4mu invDisp VR", planes_data_samples(), PLANES_VR_4MU, False),
             ("sim bkg, two-LJ invDisp CR", planes_background_samples(), PLANES_CR_2LJ, True)]
    print("transverse balance and missing momentum of the two leading lepton jets")
    cut = f"{PLANES_MET_THRESHOLD:.0f}"
    print(f"{'set':26s} {'pairs':>9s} {'ptAsym<0.2':>11s} {'<ptVecSum>':>11s} "
          f"{'<PF MET>':>9s} {'<corr MET>':>11s} {'corr<' + cut:>9s} "
          f"{'both<' + cut:>9s}")
    print("-" * 104)
    rows = []
    for label, samples, channel, weighted in sets:
        # values carry the cross-section weights of the simulated samples; the pair count
        # stays the unweighted one, so a weighted row cannot look better supported than it is
        asymmetry = planes_counts(samples, "lj_pair_ptAsym", channel, weighted=weighted)
        vector_sum = planes_counts(samples, "lj_pair_ptVecSum", channel, weighted=weighted)
        met = planes_counts(samples, "lj_pair_met", channel, weighted=weighted)
        corrected = planes_counts(samples, "lj_pair_metCorr", channel, weighted=weighted)
        joint = planes_counts(samples, "lj_pair_metCorr_vs_ptVecSum", channel,
                              weighted=weighted)
        entries = (asymmetry if not weighted
                   else planes_counts(samples, "lj_pair_ptAsym", channel))
        if asymmetry is None or entries is None or entries.sum() < 1:
            continue
        edge = int(round(PLANES_MET_THRESHOLD / 400.0 * joint.shape[0]))
        both = float(joint[:edge, :edge].sum()) / float(joint.sum())
        below = kept_fraction(corrected, PLANES_MET_THRESHOLD, 0, 400,
                              direction="below")["fraction"]
        print(f"{label:26s} {entries.sum():9.0f} "
              f"{kept_fraction(asymmetry, 0.2, direction='below')['fraction']:11.4f} "
              f"{binned_mean(vector_sum, 0, 400)[0]:11.1f} "
              f"{binned_mean(met, 0, 400)[0]:9.1f} "
              f"{binned_mean(corrected, 0, 400)[0]:11.1f} "
              f"{below:9.4f} {both:9.4f}")
        rows.append({
            "label": label, "channel": channel, "n": float(entries.sum()),
            "ptAsym": kept_fraction(asymmetry, 0.2, direction="below")["fraction"],
            "met": binned_mean(met, 0, 400)[0],
            "metCorr": binned_mean(corrected, 0, 400)[0],
            "metCorr_below": below, "both_below": both})
    print()
    print("The DSA correction subtracts the displaced-standalone constituents of the two")
    print("leading lepton jets back out of the missing-momentum vector, so on a signal pair")
    print("it should shrink it.  In every set above it does the opposite: read the two mean")
    print("columns side by side before using either.")
    return rows


# --------------------------------------------------------------------------------------
# the closing summary
# --------------------------------------------------------------------------------------

def planes_summary_table(min_m_xx=PLANES_MIN_M_XX):
    """One compact table of every separation number this section measures."""
    print(f"summary, bound-state masses at or above {min_m_xx:.0f} GeV first")
    print("(this campaign's grid holds 100, 500 and 1000 GeV only, so the 200 to 250 GeV")
    print("range at the low edge of the search is not represented in any row below)")
    print()
    fake = planes_counts(planes_data_samples(), PLANES_MU_HIST, PLANES_CR_2LJ).sum(axis=1)
    fake_eff = kept_fraction(fake, 0.4)["fraction"]
    fake_1mu = planes_counts(planes_data_samples(), PLANES_MU_HIST,
                             PLANES_CR_1MU).sum(axis=1)
    fake_1mu_eff = kept_fraction(fake_1mu, 0.4)["fraction"]
    print(f"{'point':16s} {'mu-LJ eff(|cos|>0.4)':>21s} {'vs 1muLJ CR':>12s} "
          f"{'vs 2LJ CR':>10s} {'egm 2e frac':>12s} {'A gen':>9s} {'A reco':>9s} "
          f"{'ptAsym<0.2':>11s}")
    print("-" * 108)
    order = []
    for final_state in ("4Mu", "2Mu2E"):
        for key in planes_points(final_state, min_m_xx):
            order.append((final_state, key))
        for key in planes_points(final_state, 0.0):
            if key[0] < min_m_xx:
                order.append((final_state, key))
    for final_state, (m_xx, m_zd) in order:
        channel = "4mu" if final_state == "4Mu" else "2mu2e"
        samples = planes_points(final_state, 0.0)[(m_xx, m_zd)]
        angle = planes_counts(samples, PLANES_MU_HIST, channel)
        if angle is None:
            continue
        eff_s = kept_fraction(angle.sum(axis=1), 0.4)["fraction"]
        # the resolved-electron fraction is a statement about the 2mu2e final state; a 4Mu
        # sample has almost no egamma lepton jets and its value would be noise on a handful
        multiplicity = (planes_counts(samples, "egm_lj_electronN", channel)
                        if final_state == "2Mu2E" else None)
        resolved = (float(multiplicity[2:].sum()) / float(multiplicity.sum())
                    if multiplicity is not None and multiplicity.sum() > 0 else float("nan"))
        reco = planes_counts(samples, "lj_pair_Phi", channel)
        gen = planes_counts(samples, "gen_pair_Phi", channel)
        a_reco = (harmonic_moment(reco, -np.pi, np.pi)["coefficient"]
                  if reco is not None and reco.sum() >= MIN_SHAPE_ENTRIES else float("nan"))
        a_gen = (harmonic_moment(gen, -np.pi, np.pi)["coefficient"]
                 if gen is not None and gen.sum() >= MIN_SHAPE_ENTRIES else float("nan"))
        asymmetry = planes_counts(samples, "lj_pair_ptAsym", channel)
        balance = (kept_fraction(asymmetry, 0.2, direction="below")["fraction"]
                   if asymmetry is not None and asymmetry.sum() >= MIN_SHAPE_ENTRIES
                   else float("nan"))
        def cell(value, width, digits=3, signed=True):
            """A blank rather than a nan, so an unmeasured entry cannot be misread."""
            if not np.isfinite(value):
                return f"{'--':>{width}s}"
            return f"{value:+{width}.{digits}f}" if signed else f"{value:{width}.{digits}f}"

        print(f"{final_state + ' ' + str(int(m_xx)) + '/' + str(m_zd):16s} {eff_s:21.4f} "
              f"{eff_s / np.sqrt(fake_1mu_eff):12.2f} {eff_s / np.sqrt(fake_eff):10.2f} "
              f"{cell(resolved, 12, 5, False)} {cell(a_gen, 9)} {cell(a_reco, 9)} "
              f"{cell(balance, 11, 4, False)}")
    print()
    print(f"fake reference, |cos(theta*)| > 0.4: one mu-LJ CR keeps {fake_1mu_eff:.4f}, "
          f"two-LJ inverted-displacement CR keeps {fake_eff:.4f}")
    for label, samples, channel, weighted in [
            ("2018 data, two-LJ invDisp CR", planes_data_samples(), PLANES_CR_2LJ, False),
            ("2018 data, 4mu invDisp VR", planes_data_samples(), PLANES_VR_4MU, False),
            ("simulated bkg, two-LJ CR (weighted)", planes_background_samples(),
             PLANES_CR_2LJ, True)]:
        phi = planes_counts(samples, "lj_pair_Phi", channel, weighted=weighted)
        raw_phi = (phi if not weighted
                   else planes_counts(samples, "lj_pair_Phi", channel))
        asymmetry = planes_counts(samples, "lj_pair_ptAsym", channel, weighted=weighted)
        plane = planes_counts(samples, PLANES_PAIR_PLANE, channel)
        if phi is None or plane is None:
            continue
        moment = harmonic_moment(phi, -np.pi, np.pi)
        moment["sigma_coefficient"] = harmonic_moment(raw_phi, -np.pi,
                                                      np.pi)["sigma_coefficient"]
        correlation = plane_correlation(plane, replicas=400)
        print(f"  {label}: A = {moment['coefficient']:+.3f} "
              f"+- {moment['sigma_coefficient']:.3f}, "
              f"r(cos theta_A) = {correlation['r']:+.4f} +- {correlation['sigma']:.4f}, "
              f"ptAsym < 0.2 keeps "
              f"{kept_fraction(asymmetry, 0.2, direction='below')['fraction']:.4f}")


def planes_inventory():
    """Print what the plane campaign holds, sample class by sample class."""
    print(f"campaign      : {PLANES_CAMPAIGN}")
    print(f"merged inputs : {PLANES_INPUT}")
    print(f"signal points : 4Mu {list(planes_points('4Mu', 0.0))}")
    print(f"                2Mu2E {list(planes_points('2Mu2E', 0.0))}")
    print(f"data samples  : {len(planes_data_samples())}")
    print(f"bkg samples   : {len(planes_background_samples())}")
    print(f"data channels : {', '.join(sorted(PLANES_DATA_CHANNELS))}")
    return None


# ======================================================================================
# signal in the data regions
# ======================================================================================
#
# Every data region of this study is a control or validation region, but none of them is
# free of signal: the two validation regions and the two-lepton-jet control region are the
# signal-region selection with every lepton jet required to be prompt, so a signal with a
# short proper lifetime lands in them.  The table below counts how much, region by region
# and campaign by campaign, with each campaign's own histogram and its own normalisation.

#: the histogram each campaign's region counts are read from.  The first campaign and the
#: pair campaign fill one entry per muon lepton jet with two or more muons; the decay-plane
#: histogram is filled once per event with two qualifying lepton jets, and is the one the
#: cos(2 Phi) coefficient of section j is measured on.
CONTAMINATION_HIST = "mu_lj_mumu_absCosTheta"
CONTAMINATION_PLANES_HIST = "lj_pair_Phi"

#: the cross section every yield in the table assumes, in pb (1 fb, the processor's
#: default for the signal grid).
CONTAMINATION_XS_PB = 0.001


def _one_fb_weight(year, sum_weights):
    """Luminosity x 1 fb / stored sum of generator weights: the factor `lumi_xs_weight` and
    `planes_weight` return for a signal sample, written out so that it does not depend on
    which cross-section lookup is importable."""
    if not sum_weights:
        return None
    return load_run_periods()[year]["lumi"] * CONTAMINATION_XS_PB / sum_weights


def _contamination_regions():
    """Each data region as each campaign draws it:
    `(campaign, unit, region, final states, data total, counter, weight)`."""
    def first(sample, channel):
        histogram = get_hist(sample, CONTAMINATION_HIST, channel, quiet=True)
        return None if histogram is None else total_entries(histogram)

    def first_weight(sample):
        record = load_sample(sample)
        return None if record is None else _one_fb_weight(record["year"], record["sum_weights"])

    def pair(sample, channel):
        counts = pair_counts(sample, PAIR_MJJ, channel)
        return None if counts is None else float(np.asarray(counts).sum())

    def pair_weight(sample):
        if _pair_load(sample) is None:
            return None
        return _one_fb_weight(YEAR, _PAIR_SUM_WEIGHTS.get(sample, 0.0))

    def planes(sample, channel):
        counts = planes_counts(sample, CONTAMINATION_PLANES_HIST, channel)
        return None if counts is None else float(np.asarray(counts).sum())

    def planes_one_fb(sample):
        record = _planes_load(sample)
        return None if record is None else _one_fb_weight(record["year"], record["sum_weights"])

    first_data = data_samples()
    pair_data = pair_data_samples()
    planes_data = planes_data_samples()
    both = ("4Mu", "2Mu2E")
    regions = []
    for region, final_states in (
            ("test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto", ("4Mu",)),
            ("test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto", ("2Mu2E",)),
            (CLEANED_CR, both), ("data_control_region_1muLj_cosmic_veto", both),
            (UNCUT_CR, both)):
        regions.append(("first campaign", "muon-LJ entries", region, final_states,
                        first(first_data, region), first, first_weight, signal_samples()))
    regions.append(("pair campaign", "muon-LJ entries", PAIR_VR, ("4Mu",),
                    pair(pair_data, PAIR_VR), pair, pair_weight,
                    [n for n in _pair_listing() if parse_signal(n) is not None]))
    planes_signal = [n for n in _planes_listing() if parse_signal(n) is not None]
    for region, final_states in ((PLANES_CR_2LJ, both), (PLANES_VR_4MU, ("4Mu",)),
                                 (PLANES_VR_2MU2E, ("2Mu2E",))):
        regions.append(("decay-plane campaign", "events", region, final_states,
                        planes(planes_data, region), planes, planes_one_fb, planes_signal))
    return regions


def _hi_lo_with_error(histogram):
    """The high-to-low density ratio of a histogram and its statistical error."""
    counts = raw_counts(histogram)
    edges = np.asarray(histogram.axes[0].edges, dtype=float)
    high = counts[edges[:-1] >= 0.8 - 1e-9].sum()
    low = counts[edges[1:] <= 0.2 + 1e-9].sum()
    value = _hi_lo_ratio(histogram)
    if not np.isfinite(value) or high <= 0 or low <= 0:
        return value, float("nan")
    return value, value * float(np.sqrt(1.0 / high + 1.0 / low))


def signal_contamination_table(min_m_xx=200.0, full=False):
    """Print how much signal each data region of the study would hold, and return the rows.

    One row is returned, and printed only with `full=True` (about 290 lines), for every
    signal sample with m_XX >= `min_m_xx` and every data region it can enter, as each
    campaign draws that region: the entries in the sample's signal channel (`4mu` or
    `2mu2e`: trigger, two lepton jets and the final state, before the displacement and
    isolation requirements of `4mu_sr` and `2mu2e_sr`), the entries in the region, their
    ratio, the expected count at a cross section of 1 fb, and that count as a percentage of
    the data entries in the region.  The first and the pair campaign count muon lepton jets,
    so their expected counts are muon-LJ entries times the per-event weight, not events; the
    decay-plane campaign counts events.  The weight is luminosity x 1 fb divided by the
    stored sum of generator weights of the campaign's own file.

    Always printed: the largest share per region, the factor by which it scales at the
    cross section `configs/cross_sections.yaml` lists for that sample, why the two-muon-
    two-electron validation region holds so few data events in the decay-plane campaign,
    and the effect of the contamination on the two numbers the study quotes from the
    validation regions: the high-to-low ratio of the four-muon one and the cos(2 Phi)
    coefficient A of the two-muon-two-electron one.
    """
    regions = _contamination_regions()
    short = {**PLANES_SHORT, UNCUT_CR: "1 muLJ CR uncut", CLEANED_CR: "1 muLJ CR cleaned",
             "data_control_region_1muLj_cosmic_veto": "1 muLJ CR cosmic veto"}
    rows = []
    for campaign, unit, region, final_states, data_total, counter, weight, samples in regions:
        picked = []
        for name in samples:
            info = parse_signal(name)
            if info["final_state"] not in final_states or info["m_xx"] < min_m_xx:
                continue
            picked.append(info)
        picked.sort(key=lambda i: (i["final_state"] != "4Mu", i["m_xx"], i["m_zd"],
                                   i["ctau_mm"]))
        if full:
            print(f"{campaign}, {short[region]} -- data: "
                  f"{0 if data_total is None else data_total:,.0f} {unit}")
            print(f"  {'sample':31s} {'ctau mm':>8s} {'channel':>9s} {'region':>9s} "
                  f"{'reg/ch':>7s} {'N at 1 fb':>10s} {'% of data':>10s}")
        for info in picked:
            name = info["name"]
            in_channel = counter(name, info["channel"])
            in_region = counter(name, region)
            factor = weight(name)
            if in_region is None or factor is None:
                continue
            expected = in_region * factor
            percent = 100.0 * expected / data_total if data_total else float("nan")
            ratio = in_region / in_channel if in_channel else float("nan")
            if full:
                print(f"  {name:31s} {info['ctau_mm']:8.4g} {in_channel or 0:9.0f} "
                      f"{in_region:9.0f} {ratio:7.3f} {expected:10.2f} {percent:10.4f}")
            rows.append({"campaign": campaign, "unit": unit, "region": region,
                         "sample": name, "final_state": info["final_state"],
                         "m_xx": info["m_xx"], "m_zd": info["m_zd"],
                         "ctau_mm": info["ctau_mm"], "channel_entries": in_channel,
                         "region_entries": in_region, "ratio": ratio,
                         "n_1fb": expected, "data_entries": data_total,
                         "percent": percent})
        if full:
            print()

    with open(os.path.join(_CONFIG_DIR, "cross_sections.yaml"), encoding="utf-8") as handle:
        listed = yaml.safe_load(handle)
    print(f"largest share of the data per region, m_XX >= {min_m_xx:.0f} GeV, at 1 fb, and at")
    print("the cross section configs/cross_sections.yaml lists for the same sample")
    print(f"  {'campaign, region':44s} {'sample':31s} {'% at 1 fb':>9s} {'>= 1 %':>7s} "
          f"{'listed fb':>9s} {'factor':>7s} {'% listed':>9s}")
    worst = {}
    for campaign, _unit, region, *_rest in regions:
        mine = [r for r in rows if r["campaign"] == campaign and r["region"] == region
                and np.isfinite(r["percent"])]
        if not mine:
            continue
        top = max(mine, key=lambda r: r["percent"])
        worst[(campaign, region)] = top
        listed_pb = listed.get(top["sample"])
        factor = listed_pb / CONTAMINATION_XS_PB if listed_pb else float("nan")
        above = sum(1 for r in mine if r["percent"] >= 1.0)
        where = f"{campaign}, {short[region]}"
        print(f"  {where:44s} {top['sample']:31s} {top['percent']:9.4f} "
              f"{above:3d}/{len(mine):<3d} "
              f"{(listed_pb or float('nan')) * 1000.0:9.3f} {factor:7.3f} "
              f"{top['percent'] * factor:9.4f}")
    print("  (the expected count is proportional to the cross section, so every percentage")
    print("   above scales by the factor in the same row when the listed value is assumed)")
    by_mass = {}
    for row in rows:
        if row["sample"] in listed:
            by_mass.setdefault(row["m_xx"], set()).add(listed[row["sample"]])
    print("  values listed for the samples in the rows above: " + ", ".join(
        f"m_XX = {m_xx:.0f} GeV: " + " / ".join(f"{v * 1000.0:.3f} fb" for v in sorted(values))
        for m_xx, values in sorted(by_mass.items())))
    print()
    _contamination_effects(regions, rows)
    return rows


def _contamination_effects(regions, rows):
    """The second half of `signal_contamination_table`: why the two-muon-two-electron
    validation region is so small in the decay-plane campaign, and what the expected signal
    does to the high-to-low ratio and the cos(2 Phi) coefficient quoted from the two
    validation regions."""
    data_in = {(r[0], r[2]): r[4] for r in regions}
    four_first = data_in[("first campaign", PLANES_VR_4MU)]
    four_planes = data_in[("decay-plane campaign", PLANES_VR_4MU)]
    two_first = data_in[("first campaign", PLANES_VR_2MU2E)]
    two_planes = data_in[("decay-plane campaign", PLANES_VR_2MU2E)]
    multiplicity = planes_counts(planes_data_samples(), "egm_lj_electronN", PLANES_VR_2MU2E)
    print("data in the validation regions, per muon lepton jet (first campaign) and per event")
    print("with two qualifying lepton jets (decay-plane campaign):")
    print(f"  4mu VR  : {four_first:,.0f} muon-LJ entries, {four_planes:,.0f} events, "
          f"ratio {four_first / four_planes:.2f}")
    print(f"  2mu2e VR: {two_first:,.0f} muon-LJ entries, {two_planes:,.0f} events, "
          f"ratio {two_first / two_planes:.1f}")
    if multiplicity is not None:
        print(f"  egamma LJs in the 2mu2e VR data: {multiplicity.sum():,.0f}, of which "
              f"{multiplicity[2:].sum():,.0f} hold two or more electrons")
    print()

    # effect on the four-muon validation region's high-to-low ratio
    observed = get_hist(data_samples(), CONTAMINATION_HIST, PLANES_VR_4MU, quiet=True)
    if observed is not None:
        value, error = _hi_lo_with_error(observed)
        print(f"4mu VR high-to-low ratio of the data: {value:.3f} +- {error:.3f}; with each")
        print("sample's expected entries at 1 fb subtracted (samples above 1 % of the data):")
        shifts = []
        for row in rows:
            if (row["campaign"] != "first campaign" or row["region"] != PLANES_VR_4MU
                    or row["percent"] < 1.0):
                continue
            signal = get_hist(row["sample"], CONTAMINATION_HIST, PLANES_VR_4MU, quiet=True)
            remainder = observed + _scale_hist(signal, -row["n_1fb"] / row["region_entries"])
            corrected = _hi_lo_ratio(remainder)
            shifts.append((row["sample"], row["percent"], _hi_lo_ratio(signal), corrected,
                           value - corrected))
        for name, percent, signal_ratio, corrected, shift in sorted(
                shifts, key=lambda item: -abs(item[4])):
            print(f"  {name:31s} {percent:6.2f} %  signal {signal_ratio:5.2f}  "
                  f"data minus signal {corrected:.3f}  shift {shift:+.3f}")
        print()

    # effect on the cos(2 Phi) coefficient of the two-muon-two-electron validation region
    region_phi = planes_counts(planes_data_samples(), CONTAMINATION_PLANES_HIST,
                               PLANES_VR_2MU2E)
    if region_phi is not None and region_phi.sum() > 1:
        data_moment = harmonic_moment(region_phi, -np.pi, np.pi)
        print(f"2mu2e VR cos(2 Phi) coefficient of the data: A = "
              f"{data_moment['coefficient']:+.3f} +- {data_moment['sigma_coefficient']:.3f} "
              f"on {data_moment['n']:.0f} events;")
        print("a signal making up a fraction f of them, with coefficient A_s, moves A by")
        print("f x (A_s - A) (samples above 1 % of the data, at 1 fb):")
        shifts = []
        for row in rows:
            if (row["campaign"] != "decay-plane campaign" or row["region"] != PLANES_VR_2MU2E
                    or row["percent"] < 1.0):
                continue
            signal = harmonic_moment(planes_counts(row["sample"], CONTAMINATION_PLANES_HIST,
                                                   PLANES_VR_2MU2E), -np.pi, np.pi)
            fraction = row["percent"] / 100.0
            shifts.append((row["sample"], row["percent"], signal,
                           fraction * (signal["coefficient"] - data_moment["coefficient"])))
        for name, percent, signal, shift in sorted(shifts, key=lambda item: -abs(item[3])):
            print(f"  {name:31s} {percent:6.2f} %  A_s = {signal['coefficient']:+.3f} +- "
                  f"{signal['sigma_coefficient']:.3f}  shift {shift:+.3f}")
