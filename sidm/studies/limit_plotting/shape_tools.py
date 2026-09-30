"""Per-ABCD-region *shape* extraction, for both coffea histogram layouts.

The counting and four-bin ABCD cards built by :mod:`datacard_tools` use one
number per ABCD region.  This module extracts a binned *distribution* per
region instead, which is what a ``RooParametricHist`` ABCD fit needs: the
control regions supply one free parameter per bin and the signal region's
background is a per-bin formula in them.

The intended observable is the lepton-jet--lepton-jet invariant mass,
``m_ljlj``.  **No histogram in the current productions carries that observable
on the ``abcd_region`` axis** -- see :data:`OBSERVABLES` and
:func:`missing_observable_help` -- so until a production provides one, two
stand-ins are available:

``mulj_pt`` and friends
    Real multi-bin distributions per ABCD region, from the histograms that do
    exist.  Right for exercising the machinery; not the physics observable.

:func:`mljlj_two_bin`
    A genuine two-bin ``m_ljlj`` distribution, split at the 150 GeV signal
    region cut, recovered by pairing the SR channel (``m_ljlj >= 150``) with
    the ``VR_invMass`` channel (``m_ljlj < 150``).  Those two selections are
    identical in every other cut, so together they tile the plane.  Coarse, but
    it is the observable itself and needs no new production.

MC only
-------
Every entry point here reads the background and signal merges.  The data
directory is refused outright by :func:`_reject_data_dir`, and the fail-closed
blinding guard in :mod:`datacard_tools` still applies on top of that, so a
mislabelled sample has its signal region withheld rather than returned.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import hist
import numpy as np

import datacard_tools as dt
from datacard_tools import ABCD_REGIONS, SR_ABCD_REGION, BlindingError, Yield


class ShapeError(RuntimeError):
    """A shape was requested that the inputs cannot supply."""


# --------------------------------------------------------------------------- #
# Observables
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Observable:
    """A distribution that exists (or should exist) per ABCD region.

    ``hists`` maps a channel name onto the coffea histogram holding this
    observable for that channel; ``axis`` is the name of its value axis.
    ``available`` is False for observables a production has not produced yet --
    they are declared here so the request for one gives a useful error instead
    of a ``KeyError``.
    """

    name: str
    hists: dict
    axis: str
    label: str
    available: bool = True
    note: str = ""


OBSERVABLES = {
    # --- the observable we actually want, and do not yet have ---------------
    "ljlj_invmass": Observable(
        name="ljlj_invmass",
        hists={"SR_2mu2e": "abcd_2mu2e_ljlj_invmass",
               "SR_4mu": "abcd_4mu_ljlj_invmass"},
        axis="lj_lj_invmass",
        label=r"$m_{LJLJ}$ [GeV]",
        available=False,
        note="not produced yet; the abcd_* collection carries only per-LJ "
             "kinematics. See missing_observable_help().",
    ),
    # --- what the current productions do carry per ABCD region --------------
    "mulj_pt": Observable(
        name="mulj_pt",
        hists={"SR_2mu2e": "abcd_2mu2e_mulj_pt", "SR_4mu": "abcd_4mu_mulj0_pt"},
        axis="mu_ljs_pt",
        label=r"$\mu$-LJ $p_{T}$ [GeV]",
        note="the same histogram the counting cards integrate over",
    ),
    "mulj_mass": Observable(
        name="mulj_mass",
        hists={"SR_2mu2e": "abcd_2mu2e_mulj_mass", "SR_4mu": "abcd_4mu_mulj0_mass"},
        axis="mu_ljs_mass",
        label=r"$\mu$-LJ mass [GeV]",
    ),
    "egmlj_pt": Observable(
        name="egmlj_pt",
        hists={"SR_2mu2e": "abcd_2mu2e_egmlj_pt"},
        axis="egm_ljs_pt",
        label=r"$e\gamma$-LJ $p_{T}$ [GeV]",
        note="2mu2e only; there is no egm LJ in the 4mu final state",
    ),
}


def missing_observable_help(name="ljlj_invmass"):
    """What a production would have to add to make an observable available."""
    obs = OBSERVABLES[name]
    wanted = "\n".join(f"      {ch}: {h}" for ch, h in obs.hists.items())
    return (
        f"Observable {name!r} is not in the current productions.\n"
        f"  The `abcd_*` histogram collection carries per-LJ kinematics only\n"
        f"  (pt, eta, phi, mass of the mu-LJ and egm-LJ), each with axes\n"
        f"  (channel, <observable>, abcd_region).\n"
        f"  To get it, the production needs a histogram filled once per event\n"
        f"  with the LJ-LJ invariant mass and the same (channel, value,\n"
        f"  abcd_region) axes, named:\n{wanted}\n"
        f"  The SR cut is m_ljlj >= 150 GeV, so the axis should start at 150\n"
        f"  for an in-SR shape fit -- or lower, if the fit is meant to span the\n"
        f"  cut. Until then use mljlj_two_bin() for the observable itself, or\n"
        f"  observable='mulj_pt' to exercise the machinery."
    )


# The SR selection and the selection with only the m_ljlj cut inverted.  They
# are identical in every other cut, so the pair tiles m_ljlj at the 150 GeV
# boundary and gives a two-bin distribution per ABCD region.
MLJLJ_SPLIT_GEV = 150.0

MLJLJ_CHANNELS = {
    "SR_2mu2e": {"low": "test_VR_2mu2e_invMass_spread_cosAlpha_mu_veto",
                 "high": "test_SR_2mu2e_spread_cosAlpha_mu_veto"},
    "SR_4mu": {"low": "test_VR_4mu_invMass_spread_cosAlpha_mu_veto",
               "high": "test_SR_4mu_spread_cosAlpha_mu_veto"},
}


# --------------------------------------------------------------------------- #
# Binned yields
# --------------------------------------------------------------------------- #
@dataclass
class Shape:
    """A binned distribution with its MC statistical variances.

    ``sumw`` and ``sumw2`` are per-bin; ``edges`` has one more entry.  ``raw``
    is the effective number of simulated events per bin, ``sumw^2/sumw2``,
    which is what the statistical treatment in a datacard has to be built from.
    """

    edges: np.ndarray
    sumw: np.ndarray
    sumw2: np.ndarray

    def __post_init__(self):
        self.edges = np.asarray(self.edges, dtype=float)
        self.sumw = np.asarray(self.sumw, dtype=float)
        self.sumw2 = np.asarray(self.sumw2, dtype=float)
        if len(self.edges) != len(self.sumw) + 1:
            raise ValueError(
                f"{len(self.edges)} edges for {len(self.sumw)} bins")

    @property
    def nbins(self):
        return len(self.sumw)

    @property
    def total(self):
        return Yield(float(self.sumw.sum()), float(self.sumw2.sum()))

    @property
    def n_eff(self):
        """Effective simulated events per bin, ``sumw^2 / sumw2``."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(self.sumw2 > 0, self.sumw ** 2 / self.sumw2, 0.0)

    def __add__(self, other):
        if not np.allclose(self.edges, other.edges):
            raise ValueError("cannot add shapes with different binning")
        return Shape(self.edges, self.sumw + other.sumw, self.sumw2 + other.sumw2)

    def scaled(self, factor):
        return Shape(self.edges, self.sumw * factor, self.sumw2 * factor ** 2)

    def rebinned(self, edges):
        """Merge into coarser ``edges``, which must be a subset of the current.

        Content below the first edge and above the last is folded into the
        first and last bins, so no events are dropped -- the same convention
        the counting cards use with ``flow=True``.
        """
        edges = np.asarray(edges, dtype=float)
        idx = np.searchsorted(self.edges, edges)
        if not np.allclose(self.edges[idx], edges):
            raise ValueError(
                "rebin edges must fall on existing bin edges; "
                f"{edges} is not a subset of the current binning")
        sumw = np.zeros(len(edges) - 1)
        sumw2 = np.zeros(len(edges) - 1)
        for i in range(len(edges) - 1):
            lo = idx[i] if i else 0                       # fold in the low tail
            hi = idx[i + 1] if i + 1 < len(edges) - 1 else len(self.sumw)
            sumw[i] = self.sumw[lo:hi].sum()
            sumw2[i] = self.sumw2[lo:hi].sum()
        return Shape(edges, sumw, sumw2)

    def to_dict(self):
        return {"edges": self.edges.tolist(), "sumw": self.sumw.tolist(),
                "sumw2": self.sumw2.tolist()}

    @classmethod
    def from_dict(cls, d):
        return cls(d["edges"], d["sumw"], d["sumw2"])


def _reject_data_dir(directory):
    """Refuse to read the data merge from this module.

    This study is MC-only by construction.  ``datacard_tools``'s blinding guard
    would already withhold region A of a data sample, but refusing the
    directory outright means a data path cannot be read at all, not even for
    the control regions.
    """
    resolved = str(Path(directory).resolve())
    for campaign in dt.CAMPAIGNS.values():
        if resolved == str(Path(campaign.data_dir).resolve()):
            raise BlindingError(
                f"shape_tools is MC-only and refuses to read the data merge "
                f"({campaign.name}). Use the background or signal directory."
            )


def region_shapes(sample_out, channel, observable, sample_name="",
                  selection=None, flow=True, allow_data_sr=False):
    """Binned distributions in each ABCD region for one sample in one channel.

    Returns ``{region_index: Shape}``, or ``None`` if this sample has no
    histogram for the observable in this channel.  ``selection`` overrides the
    channel's own selection, which is what :func:`mljlj_two_bin` uses to reach
    the inverted-mass channel.

    With ``flow=True`` the underflow and overflow are folded into the first and
    last bins, so the bin sum reproduces the counting cards' region total.

    Blinding: region A is omitted unless the sample is positively identified as
    simulation, exactly as in :func:`datacard_tools.region_yields`.
    """
    obs = OBSERVABLES[observable] if isinstance(observable, str) else observable
    if not obs.available:
        raise KeyError(missing_observable_help(obs.name))
    hist_name = obs.hists.get(channel.name)
    if hist_name is None:
        return None
    h = sample_out["hists"].get(hist_name)
    if h is None:
        return None

    selection = selection or channel.selection
    if selection not in list(h.axes["channel"]):
        return None

    axis = h.axes[obs.axis]
    edges = np.asarray(axis.edges, dtype=float)
    blind_sr = not allow_data_sr and not dt.is_simulation(sample_name, sample_out)
    return _slice_regions(h, selection, obs, edges, flow, blind_sr)


def _slice_regions(h, selection, obs, edges, flow, blind_sr):
    out = {}
    for region in ABCD_REGIONS:
        if region == SR_ABCD_REGION and blind_sr:
            continue
        sliced = h[{"channel": selection, "abcd_region": hist.loc(region)}]
        view = sliced.view(flow=True)
        sumw = np.asarray([v["value"] for v in view], dtype=float)
        sumw2 = np.asarray([v["variance"] for v in view], dtype=float)
        if flow:
            # index 0 is underflow and -1 overflow; fold them into the edges
            sumw[1] += sumw[0]
            sumw2[1] += sumw2[0]
            sumw[-2] += sumw[-1]
            sumw2[-2] += sumw2[-1]
        out[region] = Shape(edges, sumw[1:-1], sumw2[1:-1])
    return out


def collect_shapes(directory, observable, channels=None, flow=True,
                   progress=None, allow_data_sr=False):
    """Per-region shapes for every ``.coffea`` file in ``directory``.

    Returns ``{sample: {channel_name: {region_index: Shape}}}``.
    """
    _reject_data_dir(directory)
    channels = channels or dt.CHANNELS
    files = sorted(Path(directory).glob("*.coffea"))
    if not files:
        raise FileNotFoundError(f"no .coffea files under {directory}")

    out = {}
    for i, path in enumerate(files):
        if progress is not None:
            progress(i, len(files), path.name)
        for sample, sample_out in dt.read_coffea(path).items():
            per_channel = {}
            for ch_name, channel in channels.items():
                obs = OBSERVABLES[observable] if isinstance(observable, str) else observable
                if not obs.available:
                    raise KeyError(missing_observable_help(obs.name))
                hist_name = obs.hists.get(ch_name)
                if hist_name is None:
                    continue
                h = sample_out["hists"].get(hist_name)
                if h is None or channel.selection not in list(h.axes["channel"]):
                    continue
                blind_sr = not allow_data_sr and not dt.is_simulation(sample, sample_out)
                edges = np.asarray(h.axes[obs.axis].edges, dtype=float)
                per_channel[ch_name] = _slice_regions(
                    h, channel.selection, obs, edges, flow, blind_sr)
            if per_channel:
                out[sample] = per_channel
    return out


def mljlj_two_bin(directory, channels=None, flow=True, progress=None,
                  allow_data_sr=False):
    """A two-bin ``m_ljlj`` distribution per ABCD region, split at 150 GeV.

    Built by pairing each SR selection with the selection that inverts only its
    ``m_ljlj`` cut, so bin 0 is ``m_ljlj < 150`` and bin 1 is ``m_ljlj >= 150``.
    Any per-region histogram integrates to the same thing in each channel, so
    ``mulj_pt`` is used as the counter.

    Returns ``{sample: {channel_name: {region_index: Shape}}}`` with two bins.
    """
    _reject_data_dir(directory)
    channels = channels or dt.CHANNELS
    files = sorted(Path(directory).glob("*.coffea"))
    if not files:
        raise FileNotFoundError(f"no .coffea files under {directory}")
    edges = np.array([0.0, MLJLJ_SPLIT_GEV, 2 * MLJLJ_SPLIT_GEV])

    out = {}
    for i, path in enumerate(files):
        if progress is not None:
            progress(i, len(files), path.name)
        for sample, sample_out in dt.read_coffea(path).items():
            per_channel = {}
            for ch_name, channel in channels.items():
                pair = MLJLJ_CHANNELS.get(ch_name)
                hist_name = OBSERVABLES["mulj_pt"].hists.get(ch_name)
                if pair is None or hist_name is None:
                    continue
                h = sample_out["hists"].get(hist_name)
                if h is None:
                    continue
                available = list(h.axes["channel"])
                if not all(sel in available for sel in pair.values()):
                    continue
                blind_sr = not allow_data_sr and not dt.is_simulation(sample, sample_out)
                regions = {}
                for region in ABCD_REGIONS:
                    if region == SR_ABCD_REGION and blind_sr:
                        continue
                    sumw, sumw2 = [], []
                    for key in ("low", "high"):
                        total = h[{"channel": pair[key],
                                   "abcd_region": hist.loc(region)}].sum(flow=flow)
                        sumw.append(float(total.value))
                        sumw2.append(float(total.variance))
                    regions[region] = Shape(edges, sumw, sumw2)
                if regions:
                    per_channel[ch_name] = regions
            if per_channel:
                out[sample] = per_channel
    return out


# --------------------------------------------------------------------------- #
# Grouping and totals
# --------------------------------------------------------------------------- #
def group_shapes(shapes, channels=None):
    """Sum per-sample shapes into the datacard background groups.

    Returns ``{group: {channel: {region: Shape}}}``, dropping groups that are
    empty everywhere so no datacard process has a zero rate in every bin.
    """
    channels = channels or dt.CHANNELS
    grouped = {}
    for sample, per_channel in shapes.items():
        group = dt.bkg_group(sample)
        for ch_name, regions in per_channel.items():
            for region, shape in regions.items():
                slot = grouped.setdefault(group, {}).setdefault(ch_name, {})
                slot[region] = slot[region] + shape if region in slot else shape
    return {g: v for g, v in grouped.items()
            if any(s.sumw.sum() > 0 for ch in v.values() for s in ch.values())}


def total_background_shape(grouped, channel_name, region):
    """Sum every background group in one region of one channel."""
    total = None
    for per_channel in grouped.values():
        shape = per_channel.get(channel_name, {}).get(region)
        if shape is None:
            continue
        total = shape if total is None else total + shape
    if total is None:
        raise ShapeError(f"no background shape for {channel_name} region {region}")
    return total


def rebin_all(shapes, edges):
    """Apply :meth:`Shape.rebinned` across a nested shape dictionary."""
    out = {}
    for key, value in shapes.items():
        if isinstance(value, Shape):
            out[key] = value.rebinned(edges)
        else:
            out[key] = rebin_all(value, edges)
    return out


def closure_table(grouped, channel_name):
    """Per-bin ``A`` versus ``B*C/D`` for one channel.

    Returns a list of dicts, one per bin, with the region contents, the
    prediction, and the effective simulated events behind each.  A bin where
    region D is empty has an undefined prediction and is reported as ``None``
    rather than infinity.
    """
    regions = {r: total_background_shape(grouped, channel_name, r)
               for r in ABCD_REGIONS}
    rows = []
    for i in range(regions[0].nbins):
        a, b, c, d = (regions[r].sumw[i] for r in range(4))
        n_eff = {dt.ABCD_REGIONS[r]: regions[r].n_eff[i] for r in range(4)}
        pred = (b * c / d) if d > 0 else None
        rows.append({
            "bin": i,
            "lo": regions[0].edges[i], "hi": regions[0].edges[i + 1],
            "A": a, "B": b, "C": c, "D": d,
            "prediction": pred,
            "ratio": (a / pred) if pred else None,
            "n_eff": n_eff,
        })
    return rows


# --------------------------------------------------------------------------- #
# Export for the Combine stage
# --------------------------------------------------------------------------- #
def export_shapes(path, bkg_grouped, signal_shapes, observable, channels=None,
                  meta=None):
    """Write everything the workspace builder needs as JSON.

    ROOT and ``RooParametricHist`` live in the Combine release while coffea
    lives in this one, so the two halves cannot run in the same process.  This
    is the handover: plain JSON, no pickles, readable by eye.
    """
    channels = channels or dt.CHANNELS
    obs = OBSERVABLES[observable] if isinstance(observable, str) else observable
    payload = {
        "observable": {"name": obs.name, "axis": obs.axis, "label": obs.label,
                       "note": obs.note},
        "campaign": dt.CAMPAIGN.name,
        "channels": sorted(channels),
        "regions": {str(k): v for k, v in dt.ABCD_REGIONS.items()},
        "lumi_pb": dt.LUMI_PB,
        "mc_only": True,
        "background": {
            group: {ch: {str(r): s.to_dict() for r, s in regions.items()}
                    for ch, regions in per_channel.items()}
            for group, per_channel in bkg_grouped.items()
        },
        "signal": {
            name: {ch: {str(r): s.to_dict() for r, s in regions.items()}
                   for ch, regions in per_channel.items()}
            for name, per_channel in signal_shapes.items()
        },
        "meta": meta or {},
    }
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1))
    return path


def load_shapes(path):
    """Read back an :func:`export_shapes` payload with ``Shape`` objects."""
    payload = json.loads(Path(path).read_text())
    for block in ("background", "signal"):
        payload[block] = {
            name: {ch: {int(r): Shape.from_dict(s) for r, s in regions.items()}
                   for ch, regions in per_channel.items()}
            for name, per_channel in payload[block].items()
        }
    return payload


# =========================================================================== #
# The 6D inclusive layout
# =========================================================================== #
# Later productions drop the pre-binned ``abcd_region`` axis and ship one
# histogram carrying every ABCD variable as an axis of its own:
#
#     abcd6d_<channel>: (channel, iso0, iso1, disp0, disp1, abs_dphi, ljlj_mass)
#
# Nothing about the split is baked in, so the regions are defined by slicing.
# That is what makes ``m_ljlj`` available as a fit observable, lets both ABCD
# planes be built from the same events with the same region A, and lets the
# region boundaries be moved to buy background statistics.
#
# Everything above this line works on the older pre-binned layout; the two sets
# of extractors return the same ``{region: Shape}`` structure, so everything
# downstream -- grouping, closure, datacards, plots -- is shared.

SIXD_HISTS = {"SR_2mu2e": "abcd6d_2mu2e", "SR_4mu": "abcd6d_4mu"}
SIXD_SELECTIONS = {
    "SR_2mu2e": "abcd6d_inclusive_2mu2e_dxySpread_cosAlpha_mu_veto_eHEM",
    "SR_4mu": "abcd6d_inclusive_4mu_dxySpread_cosAlpha_mu_veto_eHEM",
}

ISO_CUT = 0.25
DPHI_CUT = 2.0
MASS_CUT = 150.0
PIXEL_HIT_MAX = 2.5      # disp <= 2, cut placed on the edge between 2 and 3
LOST_HIT_MIN = 0.5       # disp >= 1


@dataclass(frozen=True)
class Selection:
    """Where the region boundaries sit, so they can be deliberately loosened.

    The nominal values reproduce the analysis selection.  Raising ``iso_cut``,
    lowering ``dphi_cut`` or ``mass_floor``, or switching ``apply_disp`` off
    all enlarge region A -- which is the point: with 1-11 effective simulated
    events in the nominal region A there is nothing to test closure against, so
    a looser working point trades signal purity for enough background to
    measure whether ``A = B*C/D`` holds bin by bin.

    A loose point is for *learning*, not for limits.  Nothing here is an
    optimisation.
    """

    name: str = "nominal"
    iso_cut: float = ISO_CUT
    dphi_cut: float = DPHI_CUT
    mass_floor: float = MASS_CUT
    apply_disp: bool = True
    pixel_hit_max: float = PIXEL_HIT_MAX
    lost_hit_min: float = LOST_HIT_MIN
    # Which lepton jet's displacement requirement to apply.  "both" is the
    # analysis selection; "leg0"/"leg1" apply it to one jet only, which is the
    # halfway house that keeps far more statistics than requiring both and so
    # lets the m_ljlj shape be compared with and without the cut.
    disp_legs: str = "both"

    def describe(self):
        bits = [f"iso<{self.iso_cut:g}", f"|dphi|>={self.dphi_cut:g}",
                f"m>={self.mass_floor:g}"]
        if not self.apply_disp:
            bits.append("disp OFF")
        elif self.disp_legs == "both":
            bits.append("disp both LJs")
        else:
            bits.append(f"disp {self.disp_legs} only")
        return ", ".join(bits)


NOMINAL = Selection()


@dataclass(frozen=True)
class Plane:
    """One choice of ABCD plane: which two variables split A/B/C/D.

    ``axis_a``/``axis_b`` are the plane axes.  ``low_is_signal`` says, per axis,
    whether the signal-like side is the low one -- it is for isolation and not
    for ``|dphi|``, which is the asymmetry that makes a single generic slicer
    easy to get backwards.

    Region convention, with the closure relation ``A = B*C/D`` in both cases::

        A: both signal-like        B: axis_a control, axis_b signal-like
        C: axis_a signal-like, axis_b control        D: both control
    """

    name: str
    axis_a: str
    axis_b: str
    cut_a: float
    cut_b: float
    low_is_signal_a: bool
    low_is_signal_b: bool
    label: str


PLANES = {
    "iso_iso": Plane("iso_iso", "iso0", "iso1", ISO_CUT, ISO_CUT, True, True,
                     r"iso$\times$iso"),
    "iso_dphi": Plane("iso_dphi", "iso0", "abs_dphi", ISO_CUT, DPHI_CUT, True, False,
                      r"iso$\times|\Delta\phi|$"),
}


def _edge_index(axis, value):
    """Index of the bin edge at ``value``; raises if it is not one.

    Deliberately strict.  In this production every cut we want lands exactly on
    an edge, so a near-miss means the binning changed and the caller should know
    rather than silently get a slightly different region.
    """
    edges = np.asarray(axis.edges)
    i = int(np.argmin(np.abs(edges - value)))
    if abs(edges[i] - value) > 1e-9:
        raise ShapeError(
            f"{value} is not a bin edge of {axis.name}; edges are {edges}")
    return i


def _side(axis, cut, low, keep_under, keep_over):
    """Flow-index slice for one side of a cut on one axis.

    In-range bin ``i`` sits at flow index ``i+1``, index 0 being underflow and
    ``size+1`` overflow.  Underflow follows the low side and overflow the high
    side, so the two sides partition every event exactly once.
    """
    n = len(axis)
    i = _edge_index(axis, cut)
    if low:
        return slice(0 if keep_under else 1, i + 1)
    return slice(i + 1, n + 2 if keep_over else n + 1)


def _sr_slices(hist, channel_name, selection=NOMINAL):
    """The cuts that define the signal region, common to every plane.

    These are the requirements that are *not* plane axes: the displacement
    variables.  ``abs_dphi`` and ``iso1`` are left out because whether each is
    an SR cut or a plane axis depends on the plane, and they are applied by
    :func:`region_shapes_6d`.

    With ``selection.apply_disp`` False nothing is cut here at all, which is the
    most effective single way to buy statistics: the displacement requirements
    alone take the 4mu background from ~15700 events to ~56.
    """
    if not selection.apply_disp:
        return {}
    ax = {a.name: a for a in hist.axes}
    out = {}
    want0 = selection.disp_legs in ("both", "leg0")
    want1 = selection.disp_legs in ("both", "leg1")
    # pixel hits: displaced muons cross fewer layers, so signal is the low side
    if want0:
        out["disp0"] = _side(ax["disp0"], selection.pixel_hit_max, low=True,
                             keep_under=True, keep_over=False)
    if want1:
        if channel_name == "SR_2mu2e":
            # electron lost hits: signal-like is high, and photon-only EGM
            # lepton jets carry a sentinel value that lands in the overflow --
            # keep it, it is 65% of the signal in this channel
            out["disp1"] = _side(ax["disp1"], selection.lost_hit_min, low=False,
                                 keep_under=False, keep_over=True)
        else:
            out["disp1"] = _side(ax["disp1"], selection.pixel_hit_max, low=True,
                                 keep_under=True, keep_over=False)
    return out


def region_shapes_6d(sample_out, channel_name, plane, sample_name="",
                     selection=NOMINAL, allow_data_sr=False):
    """``m_ljlj`` distributions in the four ABCD regions of one plane.

    Returns ``{region_index: Shape}`` over the ``ljlj_mass`` bins at or above
    ``selection.mass_floor`` (pass ``None`` there to keep the whole axis), or
    ``None`` if this sample has no 6D histogram for the channel.

    Region A is omitted unless the sample is positively identified as
    simulation, exactly as :func:`datacard_tools.region_yields` does.
    """
    plane = PLANES[plane] if isinstance(plane, str) else plane
    hist_name = SIXD_HISTS.get(channel_name)
    if hist_name is None:
        return None
    h = sample_out["hists"].get(hist_name)
    if h is None:
        return None
    # named channel_value, not `selection`: the Selection argument owns that
    # name, and shadowing it here silently disabled every cut
    channel_value = SIXD_SELECTIONS[channel_name]
    if channel_value not in list(h.axes["channel"]):
        return None

    ax = {a.name: a for a in h.axes}
    order = [a.name for a in h.axes]
    view = h.view(flow=True)
    values, variances = view["value"], view["variance"]

    base = _sr_slices(h, channel_name, selection)
    # the mass axis is the observable, so it is sliced but not summed over
    # The mass axis is the observable, so it is sliced but never summed over.
    # The slice runs to the overflow, which is folded back into the top bin
    # afterwards rather than given an edge of its own -- there is no upper edge
    # to give it, and a duplicated one makes a zero-width bin.
    mass_axis = ax["ljlj_mass"]
    mass_floor = selection.mass_floor
    first = 0 if mass_floor is None else _edge_index(mass_axis, mass_floor)
    mass_slice = slice(first + 1, len(mass_axis) + 2)
    edges = np.asarray(mass_axis.edges)[first:]

    # Any plane axis that is not in this plane stays an SR cut on its
    # signal-like side; that is what keeps region A identical between planes.
    fixed = dict(base)
    if plane.axis_b != "abs_dphi":
        fixed["abs_dphi"] = _side(ax["abs_dphi"], selection.dphi_cut, low=False,
                                  keep_under=False, keep_over=True)
    if plane.axis_b != "iso1" and plane.axis_a != "iso1":
        fixed["iso1"] = _side(ax["iso1"], selection.iso_cut, low=True,
                              keep_under=True, keep_over=False)

    blind_sr = not allow_data_sr and not dt.is_simulation(sample_name, sample_out)

    def side(axis_name, cut, low_is_signal, signal_side):
        low = low_is_signal if signal_side else not low_is_signal
        return _side(ax[axis_name], cut, low=low,
                     keep_under=low, keep_over=not low)

    cut_for = {"iso0": selection.iso_cut, "iso1": selection.iso_cut,
               "abs_dphi": selection.dphi_cut}

    quadrants = {
        0: (True, True),     # A
        1: (False, True),    # B: axis_a on its control side
        2: (True, False),    # C: axis_b on its control side
        3: (False, False),   # D
    }

    out = {}
    for region, (a_sig, b_sig) in quadrants.items():
        if region == SR_ABCD_REGION and blind_sr:
            continue
        sel = dict(fixed)
        sel[plane.axis_a] = side(plane.axis_a, cut_for.get(plane.axis_a, plane.cut_a),
                                 plane.low_is_signal_a, a_sig)
        sel[plane.axis_b] = side(plane.axis_b, cut_for.get(plane.axis_b, plane.cut_b),
                                 plane.low_is_signal_b, b_sig)
        sel["ljlj_mass"] = mass_slice
        index = tuple(sel.get(name, slice(None)) for name in order)
        block_v, block_w2 = values[index], variances[index]
        # every axis but the mass one is summed away
        mass_pos = order.index("ljlj_mass")
        axes_to_sum = tuple(i for i in range(len(order)) if i != mass_pos)
        sumw = block_v.sum(axis=axes_to_sum)
        sumw2 = block_w2.sum(axis=axes_to_sum)
        sumw[-2] += sumw[-1]        # fold the overflow into the top mass bin
        sumw2[-2] += sumw2[-1]
        out[region] = Shape(edges, sumw[:-1], sumw2[:-1])
    return out


def collect_6d(directory, plane, channels=None, selection=NOMINAL,
               progress=None, allow_data_sr=False):
    """Per-region ``m_ljlj`` shapes for every ``.coffea`` file in ``directory``.

    Returns ``{sample: {channel: {region: Shape}}}``.
    """
    resolved = str(Path(directory).resolve())
    for campaign in dt.CAMPAIGNS.values():
        if resolved == str(Path(campaign.data_dir).resolve()):
            raise BlindingError(
                f"shape_tools is MC-only and refuses to read the data merge "
                f"({campaign.name}).")
    channels = channels or list(SIXD_HISTS)
    files = sorted(Path(directory).glob("*.coffea"))
    if not files:
        raise FileNotFoundError(f"no .coffea files under {directory}")

    out = {}
    for i, path in enumerate(files):
        if progress is not None:
            progress(i, len(files), path.name)
        for sample, sample_out in dt.read_coffea(path).items():
            per_channel = {}
            for ch_name in channels:
                regions = region_shapes_6d(
                    sample_out, ch_name, plane, sample_name=sample,
                    selection=selection, allow_data_sr=allow_data_sr)
                if regions:
                    per_channel[ch_name] = regions
            if per_channel:
                out[sample] = per_channel
    return out


def grouped_yields(grouped):
    """Grouped shapes -> ``{channel: {region: {process: Yield}}}``.

    The structure the datacard writers expect.  Built directly rather than by
    feeding group names back through ``datacard_tools.group_backgrounds``: that
    re-applies ``bkg_group`` to names which are already groups, and "Diboson"
    does not match its own rule, so it would be silently re-labelled "other".
    """
    out = {}
    for group, per_channel in grouped.items():
        for channel, regions in per_channel.items():
            for region, shape in regions.items():
                out.setdefault(channel, {}).setdefault(region, {})[group] = shape.total
    return out


def integrated(shapes):
    """Collapse ``{sample: {channel: {region: Shape}}}`` to single yields.

    Gives the counting and plain-ABCD inputs from the same slicing the shape
    fit uses, so the methods differ only in how the datacard treats the bins.
    """
    return {sample: {ch: {r: s.total for r, s in regions.items()}
                     for ch, regions in per_channel.items()}
            for sample, per_channel in shapes.items()}
