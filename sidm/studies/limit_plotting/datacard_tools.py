"""Signal-region yields -> Combine counting datacards for the SIDM ABCD analysis.

The merged ``.coffea`` outputs store one ``Hist`` per observable with axes
``(channel, <observable>, abcd_region)``.  ``abcd_region`` is an ``IntCategory``
with ``0 = A``, the signal region of the ABCD plane, so the SR count for a
sample is simply that histogram summed over the observable axis at
``abcd_region == 0`` for the SR channel.

Two details matter for getting the count right:

* The histograms are already scaled to ``lumi * xs`` by
  ``sidm_processor.postprocess``, so the sums are yields, not raw entries.
  Signal has no entry in ``configs/cross_sections.yaml``, so
  ``utilities.get_xs`` falls back to 1 fb -- signal yields (and therefore the
  Combine signal strength ``r``) are relative to a 1 fb reference cross section.
* The observable axis is ``Regular`` and does overflow (a few percent for the
  pt axes), so sums must use ``flow=True``.  With flow included the SR sum
  reproduces the final row of the corresponding cutflow exactly.

The ``Weight`` storage carries the sum of squared weights, so each yield comes
with its MC statistical uncertainty, which is what the per-process MC
statistical nuisance in the datacards is built from.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field, fields
from pathlib import Path

import getpass
import platform
import socket
from datetime import datetime, timezone

import hist
import yaml
from coffea.util import load

from sidm import BASE_DIR

# --------------------------------------------------------------------------- #
# Sample locations and channel definitions
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Channel:
    """One counting bin: an SR selection plus the histogram to count it with.

    ``hist_name`` is a histogram filled once per selected event in this
    channel, so summing it over its observable axis gives the event yield.
    ``signal_prefix`` is the signal-sample name prefix whose final state this
    channel targets.
    """

    name: str            # datacard bin name
    selection: str       # `channel` axis value in the coffea histograms
    hist_name: str       # histogram filled once per event in this channel
    signal_prefix: str   # "2Mu2E" or "4Mu"


def _channels(sr_2mu2e, sr_4mu):
    """Build the two SR channels from one campaign's ``channel``-axis strings.

    Only the selection string changes from campaign to campaign: the histogram
    a channel is counted with and the final state it targets are properties of
    the analysis rather than of the production, so a campaign supplies just the
    two strings and the rest is filled in here.  Keeping this a function rather
    than spelling both ``Channel``s out per campaign means a new production is
    one line in ``CAMPAIGNS``, and that the datacard bin names -- which end up
    in every card, CSV and figure -- cannot drift between campaigns.
    """
    return {
        "SR_2mu2e": Channel(
            name="SR_2mu2e",
            selection=sr_2mu2e,
            hist_name="abcd_2mu2e_mulj_pt",
            signal_prefix="2Mu2E",
        ),
        "SR_4mu": Channel(
            name="SR_4mu",
            selection=sr_4mu,
            hist_name="abcd_4mu_mulj0_pt",
            signal_prefix="4Mu",
        ),
    }


@dataclass(frozen=True)
class Campaign:
    """One production campaign: a set of merged coffea outputs to analyse.

    ``channels`` carries the ``channel``-axis strings that campaign's merge was
    produced with.  They are *not* stable across productions -- the selections
    get renamed as cuts are added -- so they belong to the campaign rather than
    to the module, and the default is the pair the first two campaigns used.
    """

    name: str
    bkg_dir: str
    signal_dir: str
    data_dir: str
    note: str = ""
    # dict() so a campaign cannot mutate another's channels through a shared
    # default; the Channels themselves are frozen.
    channels: dict = field(
        default_factory=lambda: dict(_LEGACY_CHANNELS))


# The `channel`-axis strings both of the first two campaigns were produced
# with, and therefore the default for a Campaign that does not say otherwise.
_LEGACY_CHANNELS = _channels(
    sr_2mu2e="test_SR_2mu2e_spread_cosAlpha_mu_veto",
    sr_4mu="test_SR_4mu_spread_cosAlpha_mu_veto",
)


_EOS = "/eos/uscms/store/user/dlee3/sidm_condor"

CAMPAIGNS = {
    "cosmic_veto_v1": Campaign(
        name="cosmic_veto_v1",
        bkg_dir=f"{_EOS}/ABCD_cosmic_veto/ABCD_landing_10ch_cosmic_veto_v1_bkg_full_merged_samples_v1",
        signal_dir=f"{_EOS}/ABCD_cosmic_veto/ABCD_landing_10ch_cosmic_veto_v1_signal_full_merged_samples_v1",
        data_dir=f"{_EOS}/ABCD_cosmic_veto/ABCD_landing_10ch_cosmic_veto_v1_data_full_merged_samples_v1",
        note="first campaign analysed; its merge step left metadata['is_data'] empty",
    ),
    "golden_hotspot_iso025_v1": Campaign(
        name="golden_hotspot_iso025_v1",
        bkg_dir=f"{_EOS}/ABCD_golden_hotspot_iso025_v1/ABCD_golden_hotspot_iso025_v1_bkg_merged_samples_v1",
        signal_dir=f"{_EOS}/ABCD_golden_hotspot_iso025_v1/ABCD_golden_hotspot_iso025_v1_signal_merged_samples_v1",
        data_dir=f"{_EOS}/ABCD_golden_hotspot_iso025_v1/ABCD_golden_hotspot_iso025_v1_data_merged_samples_v1",
        note="adds the eta-phi hotspot veto and 0.25 isolation; metadata['is_data'] is populated",
    ),
    "muIso_dPhi_eHEM_skimv1_v1": Campaign(
        name="muIso_dPhi_eHEM_skimv1_v1",
        bkg_dir=f"{_EOS}/ABCD_muIso_dPhi_eHEM_skimv1_v1/ABCD_muIso_dPhi_eHEM_skimv1_v1_bkg_merged_samples_v1",
        signal_dir=f"{_EOS}/ABCD_muIso_dPhi_eHEM_skimv1_v1/ABCD_muIso_dPhi_eHEM_skimv1_v1_signal_merged_samples_v1",
        data_dir=f"{_EOS}/ABCD_muIso_dPhi_eHEM_skimv1_v1/ABCD_muIso_dPhi_eHEM_skimv1_v1_data_merged_samples_v1",
        # The ABCD plane itself changed: it is now mu-LJ isolation vs
        # |dPhi(LJ0,LJ1)| rather than the two lepton-jet isolations.  That is
        # why the muLJ-egmLJ / muLJ-muLJ dPhi >= 2 event cut is *gone* -- dPhi
        # is an axis of the plane now, not a cut -- and why an isolation cut on
        # the other lepton jet had to be added in its place, since that
        # isolation is no longer an axis.  Region A is still index 0 and the
        # closure relation is still A = B*C/D, so nothing downstream changes.
        note="ABCD plane is mu-LJ isolation vs |dPhi(LJ0,LJ1)|, not iso vs iso: "
             "the muLJ-egmLJ / muLJ-muLJ dPhi >= 2 event cut is dropped (dPhi is "
             "now a plane axis) and an isolation cut on the other lepton jet "
             "replaces it (leading egmLJ iso < 0.25 in 2mu2e, subleading muLJ "
             "iso < 0.25 in 4mu); adds a full-era HEM electron veto and a "
             "'pass flags' event cut; vxySpread renamed dxySpread; run over the "
             "v1 skim.  Six channels on the axis, not ten -- the invMass VR the "
             "shape study pairs with the SR is not in this production.  Signal "
             "grid extended to 180 points (m_bound 100 and 150 GeV added), of "
             "which only the original 120 have a theory cross section.",
        channels=_channels(
            sr_2mu2e="abcd_muIso_dPhi_SR_2mu2e_dxySpread_cosAlpha_mu_veto_eHEM",
            sr_4mu="abcd_muIso_dPhi_SR_4mu_dxySpread_cosAlpha_mu_veto_eHEM",
        ),
    ),
    "sixd_inclusive_v1": Campaign(
        name="sixd_inclusive_v1",
        bkg_dir=f"{_EOS}/ABCD_6D_inclusive_eHEM_skimv1/"
                f"ABCD_6D_inclusive_eHEM_skimv1_v1_bkg_merged_samples_v1",
        signal_dir=f"{_EOS}/ABCD_6D_inclusive_eHEM_skimv1/"
                   f"ABCD_6D_inclusive_eHEM_skimv1_v1_signal_merged_samples_v1",
        data_dir=f"{_EOS}/ABCD_6D_inclusive_eHEM_skimv1/"
                 f"ABCD_6D_inclusive_eHEM_skimv1_v1_data_merged_samples_v1",
        channels=_channels("abcd6d_inclusive_2mu2e_dxySpread_cosAlpha_mu_veto_eHEM",
                           "abcd6d_inclusive_4mu_dxySpread_cosAlpha_mu_veto_eHEM"),
        note="6D inclusive: every ABCD variable is an axis (iso0, iso1, disp0, "
             "disp1, abs_dphi, ljlj_mass) rather than a pre-binned abcd_region, "
             "so regions are sliced by shape_tools and m_ljlj can be the fit "
             "observable. No abcd_region axis at all -- the counting and ABCD "
             "helpers keyed on it do not apply; use shape_tools.collect_6d instead.",
    ),
}

DEFAULT_CAMPAIGN = "golden_hotspot_iso025_v1"

# Module-level pointers to the active campaign.  They stay module-level so every
# existing call site keeps working; use_campaign() repoints them.
CAMPAIGN = CAMPAIGNS[DEFAULT_CAMPAIGN]
BKG_DIR = CAMPAIGN.bkg_dir
SIGNAL_DIR = CAMPAIGN.signal_dir
DATA_DIR = CAMPAIGN.data_dir


def use_campaign(name):
    """Point the module at a different production campaign.

    Returns the ``Campaign``. Output directories are per-campaign
    (``campaigns/<name>/...``) so two campaigns can be compared side by side
    without one overwriting the other.
    """
    global CAMPAIGN, BKG_DIR, SIGNAL_DIR, DATA_DIR, CHANNELS
    if name not in CAMPAIGNS:
        raise KeyError(f"unknown campaign {name!r}; known: {sorted(CAMPAIGNS)}")
    CAMPAIGN = CAMPAIGNS[name]
    BKG_DIR, SIGNAL_DIR, DATA_DIR = CAMPAIGN.bkg_dir, CAMPAIGN.signal_dir, CAMPAIGN.data_dir
    # The `channel`-axis strings are renamed from production to production, so
    # they move with the campaign.  A stale CHANNELS would not raise: the
    # selection would simply be absent from the axis and region_yields() would
    # return None for every sample, i.e. an empty run rather than a wrong one.
    CHANNELS = CAMPAIGN.channels
    return CAMPAIGN


def campaign_outdir(study_dir=None, campaign=None):
    """Output root for a campaign: ``<study>/campaigns/<campaign>``."""
    study_dir = Path(study_dir or Path(__file__).resolve().parent)
    return study_dir / "campaigns" / (campaign or CAMPAIGN.name)


# Signal cross section assumed by utilities.get_xs for 2Mu2E/4Mu samples, in pb.
# Combine's `r` is a multiplier on this reference.
SIGNAL_REF_XS_PB = 0.001
LUMI_PB = 59830.0  # 2018, from configs/run_periods.yaml

# Which two variables span the ABCD plane is a property of the production, not
# of this module: the first two campaigns used the two lepton-jet isolations
# (the `abcd_iso_base` hist collection), muIso_dPhi_eHEM_skimv1_v1 uses the
# mu-LJ isolation against |dPhi(LJ0,LJ1)|.  What is common to all of them, and
# all this module needs, is the numbering: region A is index 0 and is the
# signal region, D is the corner diagonally opposite it, and the closure
# relation is A = B*C/D.
SR_ABCD_REGION = 0
ABCD_REGIONS = {0: "A", 1: "B", 2: "C", 3: "D"}


# The active campaign's channels.  Module-level, like BKG_DIR and friends, so
# every existing call site keeps working; use_campaign() repoints it.
CHANNELS = CAMPAIGN.channels

# Backgrounds are merged into these groups so no datacard process is empty and
# the card stays readable.  Anything unmatched falls through to "other".
BKG_GROUPS = {
    "QCD": lambda s: s.startswith("QCD"),
    "DY": lambda s: s.startswith("DY"),
    "TT": lambda s: s.startswith("TT"),
    "Diboson": lambda s: s in {"WW", "WZ", "ZZ"},
}


def bkg_group(sample):
    """Map a background sample name onto its datacard process name."""
    for group, matches in BKG_GROUPS.items():
        if matches(sample):
            return group
    return "other"


# `<prefix>_<mMediator>GeV_<mDarkPhoton>GeV_<ctau>mm`, e.g. 2Mu2E_1000GeV_1p2GeV_0p96mm
SIGNAL_NAME = re.compile(
    r"^(?P<prefix>2Mu2E|4Mu|Comb)_(?P<mzd>[\dp]+)GeV_(?P<mdp>[\dp]+)GeV_(?P<ctau>[\dp]+)mm$"
)


def parse_signal_name(name):
    """Split a signal sample name into its physics parameters.

    Returns a dict with the final state and the three grid coordinates as
    floats (``p`` is the decimal point in these names), or ``None`` if the
    name does not look like a signal point.
    """
    m = SIGNAL_NAME.match(name)
    if not m:
        return None
    num = lambda s: float(s.replace("p", "."))
    return {
        "final_state": m["prefix"],
        "m_mediator": num(m["mzd"]),
        "m_darkphoton": num(m["mdp"]),
        "ctau": num(m["ctau"]),
    }


# --------------------------------------------------------------------------- #
# Yield extraction
# --------------------------------------------------------------------------- #
@dataclass
class Yield:
    """A yield and its MC statistical uncertainty."""

    value: float
    variance: float

    @property
    def error(self):
        return math.sqrt(max(self.variance, 0.0))

    @property
    def rel_error(self):
        """Fractional MC stat uncertainty; 0 for an empty (or negative) yield."""
        return self.error / self.value if self.value > 0 else 0.0

    def __add__(self, other):
        return Yield(self.value + other.value, self.variance + other.variance)

    def scaled(self, factor):
        """Yield rescaled by ``factor`` (e.g. to a different luminosity).

        The variance scales as ``factor**2``, so the *relative* MC statistical
        error is unchanged -- which is correct: extrapolating to more
        luminosity does not create more simulated events.
        """
        return Yield(self.value * factor, self.variance * factor ** 2)


class BlindingError(RuntimeError):
    """Raised when an operation would expose the signal region of real data."""


# Primary-dataset prefixes, used only to make the error message specific.  The
# blinding decision does NOT rely on this list -- see is_simulation.
DATA_NAME_HINTS = ("DoubleMuon", "SingleMuon", "EGamma", "SingleElectron",
                   "DoubleEG", "MuonEG", "MET", "JetHT", "Charmonium", "Muon")


def is_data_flag(sample_out):
    """Read ``metadata["is_data"]`` if it carries a usable value.

    Returns ``True``/``False`` when the accumulator holds exactly one value, and
    ``None`` when it is missing, empty, or self-contradictory. An empty
    accumulator is the signature of the merge bug in the first campaign, where
    the flag was dropped for data and simulation alike -- so "empty" has to mean
    "no information", never "not data".
    """
    if not sample_out:
        return None
    values = {bool(v) for v in (sample_out.get("metadata", {}).get("is_data") or ())}
    return values.pop() if len(values) == 1 else None


def is_simulation(sample_name, sample_out=None, cfg="cross_sections.yaml"):
    """True only if this sample can be established to be simulation.

    Two tiers, in order:

    1. ``metadata["is_data"]``, when it carries a usable value. Campaigns from
       ``ABCD_golden_hotspot_iso025_v1`` onwards populate it correctly, so this
       is the authoritative answer where it exists.
    2. Otherwise, an **allow-list** fallback: the sample counts as simulation
       only if it is a known signal point or carries a cross section in
       ``cross_sections.yaml``. Anything unrecognised is treated as data.

    The fallback exists because the first campaign's merge step emptied the
    accumulator, making it an empty set for data and simulation alike -- a naive
    ``if metadata["is_data"]`` check read real data as MC. Keeping the
    allow-list underneath means a future regression of the same kind degrades to
    "withhold too much" rather than silently leaking the signal region.
    """
    flag = is_data_flag(sample_out)
    if flag is not None:
        return not flag
    if sample_name.startswith(("2Mu2E", "4Mu")):
        return True
    try:
        from sidm.tools import utilities
        return sample_name in utilities.load_yaml(f"{BASE_DIR}/configs/{cfg}")
    except Exception:
        # If the config cannot be read we cannot prove this is simulation.
        return False


def blinding_basis(sample_out):
    """Which tier decided: the metadata flag, or the allow-list fallback."""
    return "metadata_is_data" if is_data_flag(sample_out) is not None else "allow_list_fallback"


def blinding_reason(sample_name, sample_out=None):
    """Human-readable explanation of why a sample's SR is withheld."""
    if is_data_flag(sample_out) is True:
        return f"{sample_name!r} is flagged as data by metadata['is_data']"
    if any(sample_name.startswith(h) for h in DATA_NAME_HINTS):
        return (f"{sample_name!r} looks like a data primary dataset and metadata['is_data'] "
                f"carries no usable value")
    return (f"{sample_name!r} is not a known signal point and has no entry in "
            f"cross_sections.yaml, and metadata['is_data'] carries no usable value, "
            f"so it cannot be confirmed to be simulation")


def region_yields(sample_out, channel, sample_name="", flow=True,
                  allow_data_sr=False):
    """Yields in all four ABCD regions for one sample in one channel.

    Returns ``{region_index: Yield}``, or ``None`` if the file does not contain
    this channel's histogram.

    ``flow=True`` is the default and should stay that way: the observable axes
    are ``Regular`` and overflow at the few-percent level, and only with flow
    included does the sum over all four regions reproduce the final cutflow row.

    Blinding: unless the sample is positively identified as simulation (see
    ``is_simulation``) region A -- the signal region -- is omitted from the
    result, so a blinded run simply has no SR key to plot or tabulate.  Control
    regions B, C and D are always returned: they are what the background
    estimate is built from and are not blinded.  Pass ``allow_data_sr=True``
    only when the analysis has genuinely been unblinded.
    """
    h = sample_out["hists"].get(channel.hist_name)
    if h is None:
        return None
    if channel.selection not in list(h.axes["channel"]):
        return None

    blind_sr = not allow_data_sr and not is_simulation(sample_name, sample_out)
    out = {}
    for region in ABCD_REGIONS:
        if region == SR_ABCD_REGION and blind_sr:
            continue
        sliced = h[{"channel": channel.selection, "abcd_region": hist.loc(region)}]
        total = sliced.sum(flow=flow)
        out[region] = Yield(float(total.value), float(total.variance))
    return out


def sr_yield(sample_out, channel, sample_name="", flow=True, allow_data_sr=False):
    """Signal-region (ABCD region A) yield for one sample in one channel.

    Raises ``BlindingError`` unless the sample is known simulation or
    ``allow_data_sr`` is set.
    """
    if not allow_data_sr and not is_simulation(sample_name, sample_out):
        raise BlindingError(
            f"refusing to read the signal region of {channel.name}: "
            f"{blinding_reason(sample_name, sample_out)}. Pass allow_data_sr=True only "
            f"when the analysis has actually been unblinded."
        )
    regions = region_yields(sample_out, channel, sample_name=sample_name,
                            flow=flow, allow_data_sr=True)
    return None if regions is None else regions[SR_ABCD_REGION]


def read_coffea(path):
    """Load a merged coffea file and return ``{sample: sample_output}``."""
    return load(str(path))["out"]


def collect_yields(directory, channels=None, flow=True, progress=None,
                   allow_data_sr=False):
    """Extract per-region yields for every ``.coffea`` file in ``directory``.

    Returns ``{sample: {channel_name: {region_index: Yield}}}``.  Each merged
    file holds a single sample, but the loop tolerates several.

    Signal-region entries are omitted for any sample not positively identified
    as simulation unless ``allow_data_sr`` is set, so a blinded run simply has
    no region-0 key to plot or tabulate.
    """
    channels = channels or CHANNELS
    files = sorted(Path(directory).glob("*.coffea"))
    if not files:
        raise FileNotFoundError(f"no .coffea files under {directory}")

    out = {}
    for i, path in enumerate(files):
        if progress is not None:
            progress(i, len(files), path.name)
        for sample, sample_out in read_coffea(path).items():
            per_channel = {}
            for ch_name, channel in channels.items():
                y = region_yields(sample_out, channel, sample_name=sample,
                                  flow=flow, allow_data_sr=allow_data_sr)
                if y is not None:
                    per_channel[ch_name] = y
            out[sample] = per_channel
    return out


def sr_only(yields):
    """Reduce ``collect_yields`` output to ``{sample: {channel: Yield}}``.

    Samples or channels with no region-A entry (a blinded data sample) are
    dropped rather than reported as zero.
    """
    out = {}
    for sample, per_channel in yields.items():
        reduced = {ch: regions[SR_ABCD_REGION]
                   for ch, regions in per_channel.items()
                   if SR_ABCD_REGION in regions}
        if reduced:
            out[sample] = reduced
    return out


def group_backgrounds(bkg_yields, channels=None):
    """Sum per-sample background yields into datacard process groups.

    Returns ``{channel_name: {region_index: {process: Yield}}}``.
    """
    channels = channels or CHANNELS
    grouped = {ch: {} for ch in channels}
    for sample, per_channel in bkg_yields.items():
        process = bkg_group(sample)
        for ch_name, regions in per_channel.items():
            for region, y in regions.items():
                per_region = grouped[ch_name].setdefault(region, {})
                per_region[process] = per_region.get(process, Yield(0.0, 0.0)) + y
    return grouped


def total_background(grouped, channel_name, region):
    """Total background Yield in one channel and region."""
    processes = grouped.get(channel_name, {}).get(region, {})
    return sum(processes.values(), Yield(0.0, 0.0))


# --------------------------------------------------------------------------- #
# Datacard writing
# --------------------------------------------------------------------------- #
@dataclass
class DatacardConfig:
    """Knobs for the counting-experiment datacards.

    ``lumi_unc`` is the 2018 integrated-luminosity uncertainty (2.5%).
    ``bkg_norm_unc`` is an optional flat normalisation uncertainty applied to
    every background process -- a placeholder for the ABCD closure/transfer
    uncertainty, which is not derivable from the SR count alone.  Set it to
    ``None`` to write a card with only luminosity and MC statistical nuisances.
    ``min_bkg`` floors the total background so Combine has a non-zero
    expectation to build the Asimov dataset from.

    ``use_theory_xs`` switches the signal normalisation from the 1 fb reference
    to the point's own theory cross section, which turns Combine's signal
    strength from "the limit on sigma in fb" into "the limit relative to
    theory".  The two are related by a constant factor per point, so they carry
    the same information -- but with the theory normalisation the exclusion
    boundary is at r = 1, which is what a reader expects.

    ``mc_stat`` selects how the MC statistical uncertainty on each background
    is modelled.  The SR backgrounds here come from one or two raw simulated
    events, where a log-normal is a bad description of the uncertainty, so the
    default is ``"gmN"``: Combine is told the effective number of raw entries
    ``N = (rate / error)^2`` and the per-entry weight ``alpha = rate / N``, and
    profiles the true rate with a Gamma distribution.  ``"lnN"`` falls back to
    a log-normal built from the same relative error, and ``None`` writes no MC
    statistical nuisance at all.  Signal is always treated with a lnN, since
    its relative MC error is at the sub-percent level.
    """

    lumi_unc: float = 0.025
    # Luminosity the yields should be extrapolated to, in /pb.  The histograms
    # come out of sidm_processor scaled to LUMI_PB (2018), so setting this
    # rescales both signal and background by target/LUMI_PB.  None means "leave
    # the yields at the luminosity they were produced with".
    target_lumi_pb: float | None = None
    # How the signal is normalised.  False keeps the 1 fb reference the
    # histograms were produced with, so Combine's r is the limit on sigma in fb.
    # True rescales each signal point to its own cross section from
    # cross_sections.yaml, so r is directly sigma_limit/sigma_theory and r < 1
    # means the point is excluded.
    use_theory_xs: bool = False
    # What goes on the ABCD card's `observation` line for region A.
    # "mc"         -- the MC region-A count, so the card's data is the same as the
    #                 counting card's and any B*C/D non-closure shows up in the fit
    #                 as tension between the signal region and the control regions.
    # "prediction" -- the B*C/D value itself, which is self-consistent by
    #                 construction and therefore cannot show closure tension.
    abcd_observation: str = "mc"
    signal_unc: float | None = None
    bkg_norm_unc: float | None = None
    mc_stat: str | None = "gmN"
    mc_stat_cap: float = 2.0       # clip runaway lnN values from 1-2 MC events
    min_bkg: float = 1e-4
    min_process_rate: float = 1e-6  # processes at or below this are dropped
    observation: str = "bkg"        # "bkg" -> expected/blinded; or a float


def signal_xs_scale(signal_name, config):
    """Factor rescaling a signal yield from the 1 fb reference to the target.

    Returns 1.0 unless ``config.use_theory_xs`` is set, in which case it is
    ``sigma_theory / 1 fb`` for this point.  Raises if the point has no cross
    section, rather than silently falling back to the reference -- a card whose
    normalisation is not what its header claims is worse than no card.
    """
    if config is None or not config.use_theory_xs:
        return 1.0
    from sidm.tools import utilities
    xs_pb = utilities.get_xs(signal_name, use_signal_xs=True)
    if xs_pb is None or xs_pb <= 0:
        raise ValueError(f"no usable theory cross section for {signal_name}")
    # get_xs falls back to the 1 fb reference for unknown signal names, which
    # here would be a silent mis-normalisation, so check the name is really in
    # the config rather than trusting the returned value.
    if signal_name not in utilities.load_yaml(f"{BASE_DIR}/configs/cross_sections.yaml"):
        raise ValueError(
            f"{signal_name} has no entry in cross_sections.yaml; refusing to "
            f"normalise it to a theory cross section it does not have"
        )
    return xs_pb / SIGNAL_REF_XS_PB


def lumi_scale(config):
    """Factor to rescale yields from the production luminosity to the target.

    Returns 1.0 when no target is configured.  This scales the *yields* only:
    MC statistical uncertainties keep their relative size, because the number
    of simulated events does not change with the target luminosity.
    """
    if config is None or config.target_lumi_pb is None:
        return 1.0
    return config.target_lumi_pb / LUMI_PB


def _norm_comment(signal_name, config):
    """The datacard header line describing what the signal is normalised to."""
    lumi_fb = (config.target_lumi_pb or LUMI_PB) / 1000
    if config.use_theory_xs:
        xs_fb = signal_xs_scale(signal_name, config) * SIGNAL_REF_XS_PB * 1000
        return (f"# signal normalised to its theory cross section "
                f"{xs_fb:g} fb at {lumi_fb:g} /fb, so r = sigma_limit/sigma_theory "
                f"and r < 1 means excluded")
    return (f"# signal normalised to {SIGNAL_REF_XS_PB * 1000:g} fb at {lumi_fb:g} /fb, "
            f"so r = sigma / {SIGNAL_REF_XS_PB * 1000:g} fb")


def _fmt(x):
    return f"{x:.6g}"


def build_datacard(signal_name, signal_yield, bkg_processes, channel_name,
                   config=None, floored=None):
    """Render one single-bin counting datacard as text.

    ``bkg_processes`` is ``{process: Yield}``.  Processes at or below
    ``config.min_process_rate`` are dropped -- Combine cannot handle a
    zero-rate process carrying a nuisance.
    """
    config = config or DatacardConfig()

    kept = {
        name: y for name, y in bkg_processes.items()
        if y.value > config.min_process_rate
    }

    # Combine needs something non-zero to normalise the Asimov dataset to.  If
    # every background group is empty, keep one floored process so the card is
    # still usable, and let the caller know via `floored`.
    if not kept:
        kept = {"bkg": Yield(config.min_bkg, config.min_bkg ** 2)}
        if floored is not None:
            floored.append((channel_name, signal_name))

    processes = [("signal", signal_yield)] + list(kept.items())
    names = [name for name, _ in processes]

    # Resolve the MC statistical nuisances first: a gmN row constrains the rate
    # to be exactly N * alpha, so the rate column has to be built from the same
    # rounded numbers that go into the nuisance row.
    rates = [y.value for _, y in processes]
    stat_rows = []
    for i, (name, y) in enumerate(processes):
        if not config.mc_stat or y.rel_error <= 0:
            continue
        cells = ["-"] * len(processes)
        label = f"mcstat_{channel_name}_{name}"
        if config.mc_stat == "gmN" and name != "signal":
            n_raw = max(1, round(1.0 / y.rel_error ** 2))
            alpha = float(_fmt(y.value / n_raw))  # exactly as it is written out
            cells[i] = alpha
            rates[i] = n_raw * alpha
            stat_rows.append((label, f"gmN {n_raw}", cells))
        else:
            cells[i] = 1 + min(y.rel_error, config.mc_stat_cap)
            stat_rows.append((label, "lnN", cells))

    total_bkg = sum(rates[1:])
    observation = total_bkg if config.observation == "bkg" else float(config.observation)

    width = max(14, max(len(n) for n in names) + 2)
    col = lambda vals: "".join(f"{v:<{width}}" for v in vals)
    pad = 40  # keeps the header columns lined up with the nuisance rows

    lines = [
        f"# SIDM ABCD counting datacard -- {signal_name}, {channel_name}",
        f"# signal region = ABCD region A of selection "
        f"{CHANNELS[channel_name].selection}",
        _norm_comment(signal_name, config),
        "imax 1  number of bins",
        f"jmax {len(processes) - 1}  number of background processes",
        "kmax *  number of nuisance parameters",
        "-" * 80,
        f"bin          {channel_name}",
        f"observation  {_fmt(observation)}",
        "-" * 80,
        "bin".ljust(pad) + col([channel_name] * len(processes)),
        "process".ljust(pad) + col(names),
        # Combine's convention: signal <= 0, backgrounds >= 1.
        "process".ljust(pad) + col([str(i) for i in range(len(processes))]),
        "rate".ljust(pad) + col([_fmt(r) for r in rates]),
        "-" * 80,
    ]

    # Widths sized from the content so a long nuisance name can never run into
    # the type column (Combine would then read the first value as the type).
    def nuisance(name, kind, values):
        cells = [v if isinstance(v, str) else _fmt(v) for v in values]
        return f"{name:<{_name_w}}{kind:<{_kind_w}}" + col(cells)

    _all_names = ["lumi_13TeV", "signal_norm", "bkg_norm"] + [r[0] for r in stat_rows]
    _all_kinds = ["lnN"] + [r[1] for r in stat_rows]
    _name_w = max(len(n) for n in _all_names) + 2
    _kind_w = max(len(k) for k in _all_kinds) + 2

    lines.append(nuisance("lumi_13TeV", "lnN", [1 + config.lumi_unc] * len(processes)))

    if config.signal_unc:
        lines.append(nuisance(
            "signal_norm", "lnN",
            [1 + config.signal_unc] + ["-"] * (len(processes) - 1),
        ))

    if config.bkg_norm_unc:
        lines.append(nuisance(
            "bkg_norm", "lnN",
            ["-"] + [1 + config.bkg_norm_unc] * (len(processes) - 1),
        ))

    lines.extend(nuisance(*row) for row in stat_rows)

    return "\n".join(lines) + "\n"


def write_datacards(signal_yields, bkg_grouped, outdir, channels=None,
                    config=None, only_matching_channel=True):
    """Write one single-bin counting datacard per (signal point, channel).

    Consumes the per-region structures from ``collect_yields`` and
    ``group_backgrounds`` and uses region A only; for the version that predicts
    A from the control regions see ``write_abcd_datacards``.

    With ``only_matching_channel`` each signal point is written only for the
    channel targeting its final state (2Mu2E -> SR_2mu2e, 4Mu -> SR_4mu); the
    cross-final-state yields are ~1e-4 events and carry no sensitivity.

    Returns ``(written_paths, floored)`` where ``floored`` lists the
    ``(channel, signal)`` pairs whose background had to be floored.
    """
    channels = channels or CHANNELS
    config = config or DatacardConfig()
    scale = lumi_scale(config)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    written, floored = [], []
    for signal_name in sorted(signal_yields):
        info = parse_signal_name(signal_name)
        for ch_name, channel in channels.items():
            if only_matching_channel and info and info["final_state"] != channel.signal_prefix:
                continue
            sig_regions = signal_yields[signal_name].get(ch_name) or {}
            sig = sig_regions.get(SR_ABCD_REGION)
            if sig is None:
                continue
            sig = sig.scaled(scale * signal_xs_scale(signal_name, config))
            if sig.value <= config.min_process_rate:
                continue
            bkg = {name: y.scaled(scale) for name, y
                   in bkg_grouped.get(ch_name, {}).get(SR_ABCD_REGION, {}).items()}
            card = build_datacard(
                signal_name, sig, bkg, ch_name, config=config, floored=floored,
            )
            path = outdir / f"datacard_{ch_name}_{signal_name}.txt"
            path.write_text(card)
            written.append(path)

    write_datacard_metadata(outdir, datacard_metadata(
        config, channels, method="counting", n_cards=len(written),
        warnings=[f"floored background: {c}/{s}" for c, s in floored] or None,
        bkg_grouped=bkg_grouped))
    return written, floored


# --------------------------------------------------------------------------- #
# ABCD datacards: the SR background is predicted from the control regions
# --------------------------------------------------------------------------- #
def build_abcd_datacard(signal_name, signal_regions, bkg_regions, channel_name,
                        config=None, warnings=None):
    """Render a four-bin datacard that does the ABCD arithmetic inside Combine.

    ``signal_regions`` and ``bkg_regions`` are ``{region_index: Yield}`` and
    ``{region_index: {process: Yield}}`` respectively.

    Rather than handing Combine a pre-computed ``B*C/D``, each control region
    gets its own freely floating normalisation and region A's background is
    *defined* as the product::

        bNorm, cNorm, dNorm   free rateParams, one per control region
        aNorm = bNorm*cNorm/dNorm

    so the control-region counts constrain the signal-region background inside
    the fit.  The three CR normalisations are unconstrained, which is the point:
    the data (here, MC standing in for it) determines them.

    Signal is entered in **all four** bins at its own per-region yield, not just
    in A.  That matters here -- signal leaks into region C at the tens-of-percent
    level for the high-mass points -- because it lets the fit account for the
    signal that would otherwise inflate the background prediction and quietly
    weaken the limit.
    """
    config = config or DatacardConfig()
    scale = lumi_scale(config)
    sig_scale = scale * signal_xs_scale(signal_name, config)

    bins = [f"{channel_name}_{ABCD_REGIONS[r]}" for r in sorted(ABCD_REGIONS)]
    regions = sorted(ABCD_REGIONS)

    signal_rate, bkg_rate = {}, {}
    for r in regions:
        sig = signal_regions.get(r, Yield(0.0, 0.0)).scaled(sig_scale)
        bkg = sum(bkg_regions.get(r, {}).values(), Yield(0.0, 0.0)).scaled(scale)
        signal_rate[r] = sig
        bkg_rate[r] = bkg

    # A degenerate control region makes B*C/D meaningless, so say so rather than
    # emitting a card that silently divides by ~0.
    for r in (1, 2, 3):
        if bkg_rate[r].value <= config.min_process_rate:
            if warnings is not None:
                warnings.append(
                    f"{channel_name}/{signal_name}: control region "
                    f"{ABCD_REGIONS[r]} is empty ({bkg_rate[r].value:.3g}); "
                    f"the ABCD prediction is undefined"
                )
            return None

    predicted_a = bkg_rate[1].value * bkg_rate[2].value / bkg_rate[3].value

    # Combine multiplies the `rate` column by the rateParam, so the background
    # rate is written as 1 and the entire normalisation lives in the parameter.
    processes = ["signal", "bkg"]
    width = 14
    col = lambda vals: "".join(f"{v:<{width}}" for v in vals)
    pad = 34

    obs = {r: bkg_rate[r].value for r in regions}
    if config.abcd_observation == "prediction":
        # self-consistent by construction: observation equals what the model predicts
        obs[SR_ABCD_REGION] = predicted_a

    lines = [
        f"# SIDM ABCD datacard -- {signal_name}, {channel_name}",
        f"# selection {CHANNELS[channel_name].selection}",
        f"# 4 bins = ABCD regions; background in A is bNorm*cNorm/dNorm",
        _norm_comment(signal_name, config),
        f"# MC region A = {bkg_rate[0].value:.4g}, B*C/D = {predicted_a:.4g}"
        f"  (observation in A: {config.abcd_observation})",
        f"imax {len(bins)}  number of bins",
        "jmax 1  number of background processes",
        "kmax *  number of nuisance parameters",
        "-" * 96,
        "bin".ljust(pad) + col(bins),
        "observation".ljust(pad) + col([_fmt(obs[r]) for r in regions]),
        "-" * 96,
        "bin".ljust(pad) + col([b for b in bins for _ in processes]),
        "process".ljust(pad) + col([p for _ in bins for p in processes]),
        "process".ljust(pad) + col(["0", "1"] * len(bins)),
        # background rate is 1 everywhere; the rateParams below carry the yields
        "rate".ljust(pad) + col([v for r in regions
                                 for v in (_fmt(signal_rate[r].value), "1")]),
        "-" * 96,
    ]

    # Column widths are computed from the content: a name exactly as long as a
    # fixed width would run into the next column and Combine would read the
    # first value as the nuisance type.
    nuisance_rows = []

    def nuisance(name, kind, values):
        nuisance_rows.append((name, kind,
                              [v if isinstance(v, str) else _fmt(v) for v in values]))

    # Luminosity applies to signal only: the backgrounds are free-floating, so
    # a normalisation nuisance on them would be redundant with the rateParams.
    nuisance(
        "lumi_13TeV", "lnN",
        [v for _ in bins for v in (1 + config.lumi_unc, "-")],
    )
    if config.signal_unc:
        nuisance(
            "signal_norm", "lnN",
            [v for _ in bins for v in (1 + config.signal_unc, "-")],
        )
    if config.mc_stat:
        for r in regions:
            y = signal_rate[r]
            if y.rel_error <= 0:
                continue
            cells = []
            for rr in regions:
                cells += [1 + min(y.rel_error, config.mc_stat_cap)
                          if rr == r else "-", "-"]
            nuisance(
                f"mcstat_{channel_name}_{ABCD_REGIONS[r]}_signal", "lnN", cells)

    if nuisance_rows:
        name_w = max(len(n) for n, _, _ in nuisance_rows) + 2
        kind_w = max(len(k) for _, k, _ in nuisance_rows) + 2
        lines += [f"{n:<{name_w}}{k:<{kind_w}}" + col(c) for n, k, c in nuisance_rows]

    lines.append("-" * 96)
    # One free normalisation per control region, then A defined as their product.
    for r, param in ((1, "bNorm"), (2, "cNorm"), (3, "dNorm")):
        start = bkg_rate[r].value
        lines.append(
            f"{param:<24}rateParam  {bins[r]:<16} bkg  "
            f"{_fmt(start)}  [0,{_fmt(max(start * 50, 1.0))}]"
        )
    lines.append(
        f"{'aNorm':<24}rateParam  {bins[0]:<16} bkg  "
        f"(@0*@1/@2)  bNorm,cNorm,dNorm"
    )

    return "\n".join(lines) + "\n"


def write_abcd_datacards(signal_yields, bkg_grouped, outdir, channels=None,
                         config=None, only_matching_channel=True):
    """Write one four-bin ABCD datacard per (signal point, channel).

    Returns ``(written_paths, warnings)``.
    """
    channels = channels or CHANNELS
    config = config or DatacardConfig()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    written, warnings = [], []
    for signal_name in sorted(signal_yields):
        info = parse_signal_name(signal_name)
        for ch_name, channel in channels.items():
            if only_matching_channel and info and info["final_state"] != channel.signal_prefix:
                continue
            sig_regions = signal_yields[signal_name].get(ch_name)
            if not sig_regions:
                continue
            card = build_abcd_datacard(
                signal_name, sig_regions, bkg_grouped.get(ch_name, {}), ch_name,
                config=config, warnings=warnings,
            )
            if card is None:
                continue
            path = outdir / f"datacard_abcd_{ch_name}_{signal_name}.txt"
            path.write_text(card)
            written.append(path)

    write_datacard_metadata(outdir, datacard_metadata(
        config, channels, method="abcd", n_cards=len(written),
        warnings=warnings or None, bkg_grouped=bkg_grouped))
    return written, warnings


# --------------------------------------------------------------------------- #
# Provenance sidecars
# --------------------------------------------------------------------------- #
# The merged .coffea inputs each carry a .meta.yaml recording what produced
# them. The datacards and limits derived from them get the same treatment, so a
# set of limits can be traced back through the configuration that built its
# cards to the exact coffea files -- and so runs can be sorted by the conditions
# they were made under without opening the CSVs.
SIDECAR_SUFFIX = ".meta.yaml"


def _git_dirty(repo_root):
    """True if the working tree has uncommitted changes.

    Recorded alongside the commit: a commit hash alone is misleading provenance
    when the tree it was run from had local edits.
    """
    import subprocess
    try:
        r = subprocess.run(["git", "status", "--porcelain"], cwd=str(repo_root),
                           capture_output=True, text=True, timeout=5)
        return bool(r.stdout.strip()) if r.returncode == 0 else None
    except Exception:
        return None


def _utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _repo_root():
    return Path(BASE_DIR).parent


def _read_sidecar(path):
    try:
        with open(path) as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}


def selection_definitions(directory, wanted=None):
    """Pull the full cut definitions out of an input directory's sidecars.

    The coffea sidecars record every cut that produced them -- object cuts,
    post-lepton-jet object cuts and event cuts -- so those definitions are
    copied verbatim into the datacard sidecar rather than being described
    second-hand. ``wanted`` restricts the result to the selections actually
    used; pass ``None`` for all of them.
    """
    directory = Path(directory)
    sidecars = sorted(directory.glob("*" + SIDECAR_SUFFIX))
    if not sidecars:
        return {}, []
    meta = _read_sidecar(sidecars[0])
    out, available = {}, []
    for entry in meta.get("selections") or []:
        if not isinstance(entry, dict) or "name" not in entry:
            continue
        available.append(entry["name"])
        if wanted is None or entry["name"] in wanted:
            out[entry["name"]] = entry.get("definition")
    return out, available


def _selections_agree(wanted):
    """Check the background and signal inputs were made with the same cuts.

    If they were not, every yield in the datacards is comparing apples to
    oranges, so this is recorded explicitly rather than assumed.
    """
    bkg, _ = selection_definitions(BKG_DIR, wanted)
    sig, _ = selection_definitions(SIGNAL_DIR, wanted)
    if not bkg or not sig:
        return None, ["could not read selection definitions from both inputs"]
    mismatched = [name for name in wanted
                  if bkg.get(name) != sig.get(name)]
    return (not mismatched), mismatched


def _input_provenance(directory, max_files=3):
    """Summarise the .meta.yaml sidecars of the coffea files we read.

    Records the upstream ``sidm_commit`` and creation timestamps so a change in
    the inputs is visible here rather than having to be remembered.  Only a
    sample of sidecars is read; the commits are collected as a set so a mixed
    input set shows up as more than one entry.
    """
    directory = Path(directory)
    sidecars = sorted(directory.glob("*" + SIDECAR_SUFFIX))
    info = {
        "path": str(directory),
        "n_coffea_files": len(sorted(directory.glob("*.coffea"))),
        "n_sidecars": len(sidecars),
    }
    commits, created, schemas = set(), set(), set()
    for path in sidecars[:max_files]:
        meta = _read_sidecar(path)
        if not meta:
            continue
        for key, bucket in (("sidm_commit", commits), ("created_utc", created),
                            ("schema", schemas)):
            if meta.get(key) is not None:
                bucket.add(str(meta[key]))
    info["upstream_sidm_commit"] = sorted(commits)
    info["upstream_created_utc"] = sorted(created)
    info["upstream_schema"] = sorted(schemas)
    info["sidecars_sampled"] = min(len(sidecars), max_files)
    first = _read_sidecar(sidecars[0]) if sidecars else {}
    info["chunksize"] = first.get("chunksize")
    info["unweighted_hist"] = first.get("unweighted_hist")
    info["hist_collections"] = [h.get("name") for h in (first.get("hist_collections") or [])
                                if isinstance(h, dict)]
    return info


def datacard_metadata(config, channels=None, method="counting", n_cards=None,
                      warnings=None, bkg_grouped=None):
    """Assemble the provenance record for a set of datacards."""
    from sidm.tools import metadata as sidm_metadata
    channels = channels or CHANNELS
    config = config or DatacardConfig()

    record = {
        "created_utc": _utc_now(),
        "produced_by": "sidm/studies/limit_plotting/datacard_tools.py",
        "producer": getpass.getuser(),
        "host": socket.gethostname(),
        "sidm_commit": sidm_metadata._git_rev(_repo_root()),
        "sidm_working_tree_dirty": _git_dirty(_repo_root()),
        "python": platform.python_version(),
        "campaign": {
            "name": CAMPAIGN.name,
            "note": CAMPAIGN.note,
            "bkg_dir": CAMPAIGN.bkg_dir,
            "signal_dir": CAMPAIGN.signal_dir,
            "data_dir": CAMPAIGN.data_dir,
        },
        "method": method,
        "n_datacards": n_cards,
        "signal_normalisation": "theory_xs" if config.use_theory_xs else "1fb_reference",
        "inputs": {
            "background": _input_provenance(BKG_DIR),
            "signal": _input_provenance(SIGNAL_DIR),
        },
        "analysis": {
            "abcd_regions": dict(ABCD_REGIONS),
            "signal_region": ABCD_REGIONS[SR_ABCD_REGION],
            "abcd_closure_relation": "A = B*C/D",
            "flow_included_in_sums": True,
            "signal_ref_xs_pb": SIGNAL_REF_XS_PB,
            "production_lumi_pb": LUMI_PB,
            "effective_lumi_pb": config.target_lumi_pb or LUMI_PB,
            "channels": {
                name: {"selection": ch.selection, "hist": ch.hist_name,
                       "signal_prefix": ch.signal_prefix}
                for name, ch in channels.items()
            },
        },
        # The cuts that produced the inputs, copied verbatim out of the coffea
        # sidecars so the limits record what selection they were made under.
        "selection_cuts": _used_selection_cuts(channels),
        "datacard_config": {
            f.name: getattr(config, f.name) for f in fields(config)
        },
        "blinding": {
            "policy": "metadata['is_data'] where it carries a usable value; otherwise an "
                      "allow-list fallback releasing the signal region only for a known "
                      "signal point or a sample with a cross section in cross_sections.yaml",
            "basis_for_this_campaign": _campaign_blinding_basis(),
            "fails_closed": True,
            "note": "the fallback exists because the cosmic_veto_v1 merge emptied "
                    "metadata['is_data'], making it an empty set for data and simulation "
                    "alike; a naive check read real data as MC",
        },
    }
    if warnings:
        record["warnings"] = list(warnings)
    if bkg_grouped is not None:
        record["background_yields"] = {
            ch: {ABCD_REGIONS[r]: round(total_background(bkg_grouped, ch, r).value, 6)
                 for r in sorted(ABCD_REGIONS)}
            for ch in channels
        }
    return record


def _campaign_blinding_basis():
    """Whether this campaign's files carry a usable ``is_data`` flag."""
    try:
        files = sorted(Path(BKG_DIR).glob("*.coffea"))
        if not files:
            return "unknown (no input files found)"
        out = read_coffea(files[0])
        so = out[list(out)[0]]
        return blinding_basis(so)
    except Exception:
        return "unknown"


def _used_selection_cuts(channels):
    """Full cut definitions for the selections these datacards actually use.

    Copied verbatim from the input coffea sidecars, plus a check that the
    background and signal inputs were produced with identical cuts -- if they
    were not, the yields are not comparable and the sidecar says so.
    """
    wanted = sorted({ch.selection for ch in channels.values()})
    definitions, available = selection_definitions(BKG_DIR, wanted)
    agree, mismatched = _selections_agree(wanted)
    record = {
        "used": wanted,
        "background_signal_definitions_agree": agree,
        "available_in_inputs": available,
        "source": "copied from the .meta.yaml sidecars of the input coffea files",
        "definitions": definitions,
    }
    if mismatched:
        record["MISMATCHED_between_background_and_signal"] = mismatched
    return record


def write_datacard_metadata(outdir, record):
    """Write ``datacards.meta.yaml`` into a datacard directory."""
    path = Path(outdir) / ("datacards" + SIDECAR_SUFFIX)
    with open(path, "w") as f:
        yaml.safe_dump(record, f, sort_keys=False, default_flow_style=False)
    return path


def load_limit_metadata(path):
    """Read a ``limits.meta.yaml``, given it or the directory containing it."""
    path = Path(path)
    if path.is_dir():
        path = path / ("limits" + SIDECAR_SUFFIX)
    with open(path) as f:
        return yaml.safe_load(f)


def index_limit_runs(study_dir=None, pattern="limits*"):
    """Tabulate every set of limits under ``study_dir`` by how it was produced.

    Walks the ``limits*`` directories, reads each ``limits.meta.yaml`` and
    flattens the fields most useful for telling runs apart -- method,
    normalisation, blinding, the datacard configuration, the upstream coffea
    commit -- into one row per run.  Returns a ``pandas.DataFrame``, so runs can
    be sorted or filtered by the conditions they were made under.
    """
    import pandas as pd
    study_dir = Path(study_dir or Path(__file__).resolve().parent)
    rows = []
    for directory in sorted(study_dir.glob(pattern)):
        sidecar = directory / ("limits" + SIDECAR_SUFFIX)
        if sidecar.exists():
            try:
                with open(sidecar) as f:
                    meta = yaml.safe_load(f) or {}
            except Exception as exc:
                rows.append({"run": directory.name, "error": str(exc)})
                continue
            cards = meta.get("datacards") or {}
            cfg = cards.get("datacard_config") or {}
            combine = meta.get("combine") or {}
            run = meta.get("run") or {}
            summary = meta.get("results_summary") or {}
            inputs = (cards.get("inputs") or {}).get("background") or {}
            rows.append({
                "run": directory.name,
                "created_utc": meta.get("created_utc"),
                "method": cards.get("method"),
                "signal_norm": cards.get("signal_normalisation"),
                "blinded": combine.get("blinded"),
                "abcd_observation": cfg.get("abcd_observation"),
                "use_theory_xs": cfg.get("use_theory_xs"),
                "mc_stat": cfg.get("mc_stat"),
                "lumi_unc": cfg.get("lumi_unc"),
                "bkg_norm_unc": cfg.get("bkg_norm_unc"),
                "target_lumi_pb": cfg.get("target_lumi_pb"),
                "n_limits": run.get("n_limits"),
                "n_failed": run.get("n_failed"),
                "median_exp_r": summary.get("median_expected_r"),
                "median_obs_r": summary.get("median_observed_r"),
                "combine_version": (combine.get("version") or [None])[0],
                "sidm_commit": (cards.get("sidm_commit") or "")[:8] or None,
                "tree_dirty": cards.get("sidm_working_tree_dirty"),
                "upstream_commit": ((inputs.get("upstream_sidm_commit") or [""])[0])[:8] or None,
            })
    return pd.DataFrame(rows)


def load_campaign_limits(campaign, which="limits", study_dir=None):
    """Read one limit set for one campaign, with the derived plotting columns.

    ``which`` is the directory name under ``campaigns/<campaign>/`` --
    ``limits``, ``limits_abcd``, ``limits_obs``, and so on.
    """
    import pandas as pd
    path = campaign_outdir(study_dir, campaign) / which / "limits.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path).rename(columns={
        "final_state": "topology", "m_mediator": "m_bound", "m_darkphoton": "mzd",
        "exp_m2": "expected_2p5", "exp_m1": "expected_16", "exp": "expected_50",
        "exp_p1": "expected_84", "exp_p2": "expected_97p5"})
    df["campaign"] = campaign
    df["limit_set"] = which
    from sidm.tools import utilities
    df["theory_xs_fb"] = [utilities.get_xs(x, use_signal_xs=True) * 1000.0
                          for x in df.signal]
    df["r_theory"] = df.expected_50 / df.theory_xs_fb
    if "obs" in df and df.obs.notna().any():
        df["r_theory_obs"] = df.obs / df.theory_xs_fb
    return df


def compare_campaigns(which="limits", campaigns=None, study_dir=None):
    """Join the same limit set across campaigns, one row per signal point.

    Returns a DataFrame with ``expected_50`` from each campaign side by side
    plus their ratio, so the effect of a reprocessing can be read off directly.
    """
    import pandas as pd
    campaigns = campaigns or list(CAMPAIGNS)
    frames = {c: load_campaign_limits(c, which, study_dir) for c in campaigns}
    frames = {c: f for c, f in frames.items() if f is not None}
    if len(frames) < 2:
        return pd.DataFrame()
    keys = ["channel", "signal", "topology", "m_bound", "mzd", "ctau", "theory_xs_fb"]
    names = list(frames)
    merged = frames[names[0]][keys + ["expected_50", "r_theory"]].rename(
        columns={"expected_50": f"exp_{names[0]}", "r_theory": f"r_{names[0]}"})
    for name in names[1:]:
        right = frames[name][["channel", "signal", "expected_50", "r_theory"]].rename(
            columns={"expected_50": f"exp_{name}", "r_theory": f"r_{name}"})
        merged = merged.merge(right, on=["channel", "signal"])
    if len(names) == 2:
        merged["ratio"] = merged[f"exp_{names[1]}"] / merged[f"exp_{names[0]}"]
    return merged


# --------------------------------------------------------------------------- #
# Combining the two final states
# --------------------------------------------------------------------------- #
# The 2Mu2E and 4Mu samples are two decay modes of the same model, simulated at
# the same 60 (m_bound, m_ZD, ctau) points. A search would use both, so the
# combination pairs them into one datacard with a shared signal strength: the
# 2mu2e and 4mu signal regions become separate bins, and each bin receives the
# signal from BOTH samples (the 4Mu sample leaks ~1% into SR_2mu2e, which is
# small but no longer negligible in the newer production).
#
# ASSUMPTION: both samples are normalised to the same reference, so the
# combination assumes the model yields each final state at that cross section --
# i.e. it does not model the branching split between them. Change
# `final_state_weights` if a split is known.
COMBINED_PREFIX = "Comb"


def combined_grid_key(signal_name):
    """The grid part of a signal name, shared by both final states.

    Reuses the sample name's own suffix rather than reformatting the numbers, so
    ``5p0GeV`` does not silently become ``5GeV`` and the result round-trips
    through ``parse_signal_name``.
    """
    m = SIGNAL_NAME.match(signal_name)
    if not m:
        return None
    return signal_name[len(m["prefix"]) + 1:]


def combined_point_name(key):
    """Canonical name for a combined grid point, e.g. ``Comb_800GeV_5p0GeV_5p0mm``."""
    return f"{COMBINED_PREFIX}_{key}"


def combine_final_states(signal_yields, channels=None, final_state_weights=None):
    """Sum the 2Mu2E and 4Mu samples per grid point, per channel, per region.

    Returns ``{combined_name: {channel: {region: Yield}}}`` plus a dict mapping
    each combined name back to the sample names it was built from.
    """
    channels = channels or CHANNELS
    weights = final_state_weights or {}
    combined, provenance = {}, {}
    for sample, per_channel in signal_yields.items():
        key = combined_grid_key(sample)
        if key is None:
            continue
        info = parse_signal_name(sample)
        weight = weights.get(info["final_state"], 1.0)
        name = combined_point_name(key)
        target = combined.setdefault(name, {})
        provenance.setdefault(name, []).append(sample)
        for ch_name, regions in per_channel.items():
            if ch_name not in channels:
                continue
            bucket = target.setdefault(ch_name, {})
            for region, y in regions.items():
                scaled = y.scaled(weight)
                bucket[region] = bucket[region] + scaled if region in bucket else scaled
    # only keep points where both final states were present
    complete = {n: v for n, v in combined.items() if len(provenance[n]) == 2}
    return complete, {n: sorted(provenance[n]) for n in complete}


def combined_theory_xs_pb(combined_name, provenance):
    """Theory cross section for a combined point.

    Both final states at a given ``m_bound`` share one value, so this is
    unambiguous; it raises if that ever stops being true rather than silently
    picking one.
    """
    from sidm.tools import utilities
    values = {utilities.get_xs(s, use_signal_xs=True) for s in provenance[combined_name]}
    if len(values) != 1:
        raise ValueError(
            f"{combined_name}: its final states have different theory cross sections "
            f"({sorted(values)}), so the combined normalisation is ambiguous"
        )
    return values.pop()


def build_combined_datacard(name, signal_regions, bkg_grouped, channels=None,
                            config=None, theory_xs_pb=None, floored=None):
    """Two-bin counting card: both signal regions, one shared signal strength.

    ``signal_regions`` is ``{channel: {region: Yield}}`` for the combined point.
    Only region A of each channel is used; the control regions are what the ABCD
    variant needs.
    """
    channels = channels or CHANNELS
    config = config or DatacardConfig()
    scale = lumi_scale(config)
    if config.use_theory_xs:
        if not theory_xs_pb:
            raise ValueError(f"{name}: no theory cross section for a theory-normalised card")
        scale *= theory_xs_pb / SIGNAL_REF_XS_PB

    bins, sig_rate, bkg_procs = [], {}, {}
    for ch_name in channels:
        sig = signal_regions.get(ch_name, {}).get(SR_ABCD_REGION)
        if sig is None:
            continue
        kept = {p: y.scaled(scale if False else lumi_scale(config))
                for p, y in bkg_grouped.get(ch_name, {}).get(SR_ABCD_REGION, {}).items()
                if y.value > config.min_process_rate}
        if not kept:
            kept = {"bkg": Yield(config.min_bkg, config.min_bkg ** 2)}
            if floored is not None:
                floored.append((ch_name, name))
        bins.append(ch_name)
        sig_rate[ch_name] = sig.scaled(scale)
        bkg_procs[ch_name] = kept
    if not bins:
        return None

    # one column per (bin, process); signal first in each bin
    columns = []
    for b in bins:
        columns.append((b, "signal", sig_rate[b]))
        for pname, y in bkg_procs[b].items():
            columns.append((b, pname, y))

    width = max(14, max(len(p) for _, p, _ in columns) + 2)
    col = lambda vals: "".join(f"{v:<{width}}" for v in vals)
    pad = 40

    rates = [y.value for _, _, y in columns]
    stat_rows = []
    for i, (b, pname, y) in enumerate(columns):
        if not config.mc_stat or y.rel_error <= 0:
            continue
        cells = ["-"] * len(columns)
        label = f"mcstat_{b}_{pname}"
        if config.mc_stat == "gmN" and pname != "signal":
            n_raw = max(1, round(1.0 / y.rel_error ** 2))
            alpha = float(_fmt(y.value / n_raw))
            cells[i] = alpha
            rates[i] = n_raw * alpha
            stat_rows.append((label, f"gmN {n_raw}", cells))
        else:
            cells[i] = 1 + min(y.rel_error, config.mc_stat_cap)
            stat_rows.append((label, "lnN", cells))

    observation = {b: sum(r for (bb, pn, _), r in zip(columns, rates)
                          if bb == b and pn != "signal") for b in bins}

    lines = [
        f"# SIDM combined counting datacard -- {name}",
        f"# both signal regions as separate bins, one shared signal strength",
        f"# built from: {', '.join(sorted(signal_regions.get('__provenance__', [])) )}"
        if "__provenance__" in signal_regions else
        "# 2Mu2E and 4Mu final states combined; each bin receives signal from both",
        _norm_comment_combined(config, theory_xs_pb),
        f"imax {len(bins)}  number of bins",
        f"jmax *  number of background processes",
        "kmax *  number of nuisance parameters",
        "-" * 96,
        "bin".ljust(22) + "".join(f"{b:<{width}}" for b in bins),
        "observation".ljust(22) + "".join(f"{_fmt(observation[b]):<{width}}" for b in bins),
        "-" * 96,
        "bin".ljust(pad) + col([b for b, _, _ in columns]),
        "process".ljust(pad) + col([p for _, p, _ in columns]),
        "process".ljust(pad) + col([_process_index(b, p, bins, bkg_procs)
                                    for b, p, _ in columns]),
        "rate".ljust(pad) + col([_fmt(r) for r in rates]),
        "-" * 96,
    ]

    name_w = max([len("lumi_13TeV")] + [len(n) for n, _, _ in stat_rows]) + 2
    kind_w = max([len("lnN")] + [len(k) for _, k, _ in stat_rows]) + 2
    nuis = lambda n, k, v: f"{n:<{name_w}}{k:<{kind_w}}" + col(
        [x if isinstance(x, str) else _fmt(x) for x in v])

    lines.append(nuis("lumi_13TeV", "lnN", [1 + config.lumi_unc] * len(columns)))
    if config.signal_unc:
        lines.append(nuis("signal_norm", "lnN",
                          [1 + config.signal_unc if p == "signal" else "-"
                           for _, p, _ in columns]))
    if config.bkg_norm_unc:
        lines.append(nuis("bkg_norm", "lnN",
                          ["-" if p == "signal" else 1 + config.bkg_norm_unc
                           for _, p, _ in columns]))
    lines += [nuis(*row) for row in stat_rows]
    return "\n".join(lines) + "\n"


def _process_index(bin_name, process, bins, bkg_procs):
    """Combine's convention: signal <= 0, backgrounds >= 1, numbered per card."""
    if process == "signal":
        return "0"
    order = list(bkg_procs[bin_name])
    return str(order.index(process) + 1)


def _norm_comment_combined(config, theory_xs_pb):
    lumi_fb = (config.target_lumi_pb or LUMI_PB) / 1000
    if config.use_theory_xs and theory_xs_pb:
        return (f"# signal at its theory cross section {theory_xs_pb * 1000:g} fb "
                f"per final state at {lumi_fb:g} /fb, so r < 1 means excluded")
    return (f"# signal at {SIGNAL_REF_XS_PB * 1000:g} fb PER FINAL STATE at {lumi_fb:g} /fb; "
            f"the combination assumes the model yields both at that rate")


def write_combined_datacards(signal_yields, bkg_grouped, outdir, channels=None,
                             config=None, final_state_weights=None):
    """Write one two-bin combined counting card per grid point."""
    channels = channels or CHANNELS
    config = config or DatacardConfig()
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    combined, provenance = combine_final_states(signal_yields, channels, final_state_weights)

    written, floored = [], []
    for name in sorted(combined):
        xs = combined_theory_xs_pb(name, provenance) if config.use_theory_xs else None
        card = build_combined_datacard(name, combined[name], bkg_grouped, channels,
                                       config=config, theory_xs_pb=xs, floored=floored)
        if card is None:
            continue
        path = outdir / f"datacard_comb_{name}.txt"
        path.write_text(card)
        written.append(path)

    write_datacard_metadata(outdir, datacard_metadata(
        config, channels, method="combined", n_cards=len(written),
        warnings=[f"floored background: {c}/{s}" for c, s in floored] or None,
        bkg_grouped=bkg_grouped))
    return written, floored


# --------------------------------------------------------------------------- #
# Re-deriving the ABCD regions from the 2D isolation plane
# --------------------------------------------------------------------------- #
# The abcd_region axis is pre-binned at iso <= 0.25 on both axes, but the
# underlying 2D isolation distributions are also stored, so the region
# boundaries can be moved after the fact. That allows the control regions to be
# pushed away from the signal-rich corner to cut signal contamination, while
# region A itself stays exactly where it is.
#
# Verified: with cr_cut == a_cut == 0.25 this reproduces the abcd_region yields
# exactly, in every region, for signal and background.
#
# These two histograms exist only in the campaigns whose plane is iso vs iso.
# muIso_dPhi_eHEM_skimv1_v1 spans its plane with mu-LJ isolation against
# |dPhi(LJ0,LJ1)| and stores `abcd_corr_*_iso_vs_abs_dphi` instead, so this
# re-derivation does not apply there: the lookup below misses and
# region_yields_from_plane returns None rather than guessing at a plane it was
# not written for.  The pre-binned abcd_region axis -- which is what the
# datacards actually use -- is unaffected either way.
ISO_PLANE = {
    "SR_2mu2e": {"hist": "mulj_egmlj_iso", "axes": ("mu_lj_iso", "egm_lj_iso")},
    "SR_4mu": {"hist": "mulj_mulj_iso", "axes": ("mu_lj0_iso", "mu_lj1_iso")},
}
DEFAULT_A_CUT = 0.25


def _bin_index(axis, value):
    """Index of the bin edge at ``value``; raises if it is not on an edge."""
    edges = axis.edges
    idx = int(round((value - edges[0]) / (edges[1] - edges[0])))
    if not (0 <= idx <= axis.size) or abs(edges[idx] - value) > 1e-9:
        raise ValueError(
            f"{value} is not on a bin edge of {axis.name} "
            f"(edges {edges[0]:g}..{edges[-1]:g} step {edges[1]-edges[0]:g})"
        )
    return idx


def region_yields_from_plane(sample_out, channel, a_cut=DEFAULT_A_CUT, cr_cut=None,
                             sample_name="", allow_data_sr=False):
    """ABCD yields re-derived from the 2D isolation plane.

    ``a_cut`` bounds region A on both axes and is meant to stay fixed.
    ``cr_cut`` is where the control regions start; ``None`` means "same as
    ``a_cut``", reproducing the standard definition. Setting ``cr_cut > a_cut``
    leaves an unused buffer band between A and the control regions.

    The closure relation is unaffected by the buffer: if the two isolations
    factorise then B*C/D still equals A for any ``cr_cut``, because the
    ``cr_cut`` factors cancel between numerator and denominator.

    Returns ``{region_index: Yield}``; region A is withheld unless the sample is
    established simulation, exactly as ``region_yields`` does.
    """
    cr_cut = a_cut if cr_cut is None else cr_cut
    if cr_cut < a_cut:
        raise ValueError(f"cr_cut ({cr_cut}) must be >= a_cut ({a_cut})")
    spec = ISO_PLANE.get(channel.name if hasattr(channel, "name") else channel)
    if spec is None:
        return None
    h = sample_out["hists"].get(spec["hist"])
    if h is None or channel.selection not in list(h.axes["channel"]):
        return None

    sliced = h[{"channel": channel.selection}]
    view = sliced.view(flow=True)
    values, variances = view["value"], view["variance"]
    ax0, ax1 = sliced.axes[0], sliced.axes[1]
    a0, a1 = _bin_index(ax0, a_cut), _bin_index(ax1, a_cut)
    c0, c1 = _bin_index(ax0, cr_cut), _bin_index(ax1, cr_cut)

    # Index 0 is underflow and index size+1 is overflow, so in-range bin i sits
    # at index i+1. The underflow is kept on the LOW side: a few signal events
    # carry a sentinel isolation below the axis range, and the abcd_region axis
    # counts them as isolated. Dropping them loses up to 2.6% of the signal in
    # the worst point.
    low0, low1 = slice(0, a0 + 1), slice(0, a1 + 1)
    high0, high1 = slice(c0 + 1, ax0.size + 2), slice(c1 + 1, ax1.size + 2)
    quadrants = {
        0: (low0, low1),    # A: both isolated
        1: (high0, low1),   # B
        2: (low0, high1),   # C
        3: (high0, high1),  # D
    }

    blind_sr = not allow_data_sr and not is_simulation(sample_name, sample_out)
    out = {}
    for region, (s0, s1) in quadrants.items():
        if region == SR_ABCD_REGION and blind_sr:
            continue
        out[region] = Yield(float(values[s0, s1].sum()),
                            float(variances[s0, s1].sum()))
    return out


def collect_plane_yields(directory, channels=None, a_cut=DEFAULT_A_CUT, cr_cut=None,
                         progress=None, allow_data_sr=False):
    """``collect_yields`` equivalent, but re-deriving regions from the iso plane."""
    channels = channels or CHANNELS
    files = sorted(Path(directory).glob("*.coffea"))
    if not files:
        raise FileNotFoundError(f"no .coffea files under {directory}")
    out = {}
    for i, path in enumerate(files):
        if progress is not None:
            progress(i, len(files), path.name)
        for sample, sample_out in read_coffea(path).items():
            per_channel = {}
            for ch_name, channel in channels.items():
                y = region_yields_from_plane(sample_out, channel, a_cut, cr_cut,
                                             sample_name=sample,
                                             allow_data_sr=allow_data_sr)
                if y is not None:
                    per_channel[ch_name] = y
            out[sample] = per_channel
    return out


# --------------------------------------------------------------------------- #
# The isolation-vs-dPhi ABCD plane
# --------------------------------------------------------------------------- #
# The `muIso_dPhi_eHEM_skimv1_v1` production moved the second isolation into the
# event selection, which collapses the pre-binned `abcd_region` axis onto the
# isolated half-plane -- regions C and D come out exactly zero, so no ABCD card
# can be built from it (see NOTES.md, 2026-09-16).  That production instead
# ships a 2D isolation-vs-|dPhi| histogram, and dropping the dPhi cut from the
# SR is what makes |dPhi| usable as the second ABCD axis.  This is the plane the
# `abcd_region` axis should be repointed at upstream; until then the regions are
# re-derived here.
#
# The sense of the two axes differs, which is the one thing to be careful about:
# signal is *low* isolation but *high* |dPhi| (back-to-back lepton jets), so the
# signal region is the low-iso/high-dPhi corner and D is diagonally opposite it.
DPHI_PLANE = {
    "SR_2mu2e": {"hist": "abcd_corr_2mu2e_mu_iso_vs_abs_dphi",
                 "axes": ("2mu2e_mu_iso", "2mu2e_abs_dphi")},
    "SR_4mu": {"hist": "abcd_corr_4mu_mu0_iso_vs_abs_dphi",
               "axes": ("4mu_mu0_iso", "4mu_abs_dphi")},
}
DEFAULT_DPHI_CUT = 2.0


def _nearest_bin_index(axis, value, tol_bins=0.5):
    """Index of the bin edge nearest ``value``, with the snap reported.

    Unlike ``_bin_index`` this does not require the cut to land exactly on an
    edge.  The |dPhi| axis is 32 bins over 0..pi, so the natural cut at 2.0 sits
    between edges (1.9635 and 2.0617) and has to be snapped to one of them.
    Returns ``(index, effective_value)`` so callers can report the boundary that
    was actually used rather than the one that was asked for.
    """
    edges = axis.edges
    width = edges[1] - edges[0]
    idx = int(round((value - edges[0]) / width))
    idx = max(0, min(idx, axis.size))
    if abs(edges[idx] - value) > tol_bins * width:
        raise ValueError(
            f"{value} is more than {tol_bins} bins from any edge of "
            f"{axis.name} (nearest {edges[idx]:g}, width {width:g})")
    return idx, float(edges[idx])


def dphi_plane_boundaries(sample_out, channel, iso_cut=DEFAULT_A_CUT,
                          dphi_cut=DEFAULT_DPHI_CUT):
    """The boundaries actually used, after snapping to bin edges.

    Returns ``{"iso": value, "dphi": value}`` or ``None`` if this sample has no
    plane for the channel.  Worth recording in a sidecar: the requested dPhi cut
    is generally not the one applied.
    """
    spec = DPHI_PLANE.get(getattr(channel, "name", channel))
    if spec is None:
        return None
    h = sample_out["hists"].get(spec["hist"])
    if h is None or channel.selection not in list(h.axes["channel"]):
        return None
    sliced = h[{"channel": channel.selection}]
    _, iso_eff = _nearest_bin_index(sliced.axes[0], iso_cut)
    _, dphi_eff = _nearest_bin_index(sliced.axes[1], dphi_cut)
    return {"iso": iso_eff, "dphi": dphi_eff}


def region_yields_from_dphi_plane(sample_out, channel, iso_cut=DEFAULT_A_CUT,
                                  dphi_cut=DEFAULT_DPHI_CUT, sample_name="",
                                  allow_data_sr=False):
    """ABCD yields re-derived from the isolation-vs-|dPhi| plane.

    Region convention, with the closure relation still ``A = B*C/D``::

        A: iso <  iso_cut  and  |dPhi| >= dphi_cut     (signal region)
        B: iso >= iso_cut  and  |dPhi| >= dphi_cut
        C: iso <  iso_cut  and  |dPhi| <  dphi_cut
        D: iso >= iso_cut  and  |dPhi| <  dphi_cut

    Note the asymmetry against ``region_yields_from_plane``: there both axes are
    signal-like when low, here only the isolation is.

    Returns ``{region_index: Yield}``; region A is withheld unless the sample is
    established simulation, exactly as ``region_yields`` does.
    """
    spec = DPHI_PLANE.get(getattr(channel, "name", channel))
    if spec is None:
        return None
    h = sample_out["hists"].get(spec["hist"])
    if h is None or channel.selection not in list(h.axes["channel"]):
        return None

    sliced = h[{"channel": channel.selection}]
    view = sliced.view(flow=True)
    values, variances = view["value"], view["variance"]
    ax_iso, ax_dphi = sliced.axes[0], sliced.axes[1]
    i_iso, _ = _nearest_bin_index(ax_iso, iso_cut)
    i_dphi, _ = _nearest_bin_index(ax_dphi, dphi_cut)

    # Index 0 is underflow and size+1 overflow, so in-range bin i is at i+1.
    # Underflow goes with the low side and overflow with the high side on both
    # axes, so every event is counted exactly once and the four regions sum to
    # the channel total.
    iso_low, iso_high = slice(0, i_iso + 1), slice(i_iso + 1, ax_iso.size + 2)
    dphi_low, dphi_high = slice(0, i_dphi + 1), slice(i_dphi + 1, ax_dphi.size + 2)
    quadrants = {
        0: (iso_low, dphi_high),    # A: isolated and back-to-back
        1: (iso_high, dphi_high),   # B
        2: (iso_low, dphi_low),     # C
        3: (iso_high, dphi_low),    # D
    }

    blind_sr = not allow_data_sr and not is_simulation(sample_name, sample_out)
    out = {}
    for region, (s_iso, s_dphi) in quadrants.items():
        if region == SR_ABCD_REGION and blind_sr:
            continue
        out[region] = Yield(float(values[s_iso, s_dphi].sum()),
                            float(variances[s_iso, s_dphi].sum()))
    return out


def collect_dphi_plane_yields(directory, channels=None, iso_cut=DEFAULT_A_CUT,
                              dphi_cut=DEFAULT_DPHI_CUT, progress=None,
                              allow_data_sr=False):
    """``collect_yields`` equivalent for the isolation-vs-|dPhi| plane."""
    channels = channels or CHANNELS
    files = sorted(Path(directory).glob("*.coffea"))
    if not files:
        raise FileNotFoundError(f"no .coffea files under {directory}")
    out = {}
    for i, path in enumerate(files):
        if progress is not None:
            progress(i, len(files), path.name)
        for sample, sample_out in read_coffea(path).items():
            per_channel = {}
            for ch_name, channel in channels.items():
                y = region_yields_from_dphi_plane(
                    sample_out, channel, iso_cut, dphi_cut,
                    sample_name=sample, allow_data_sr=allow_data_sr)
                if y is not None:
                    per_channel[ch_name] = y
            if per_channel:
                out[sample] = per_channel
    return out
