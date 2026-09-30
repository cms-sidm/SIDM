"""Shared loading/parsing helpers for the final-state-anatomy and
trigger-context notebooks. The merged v3 output (5 channels x 9 hist
collections, all 180 samples, every file, unweighted) lives at V3_EOS_DIR;
load_v3() mirrors _lifetime_refit.load_truthkin() with a separate local cache.

Two follow-on productions add histograms the v3 run did not book. The
displacement-conditioned run (V3_COND_EOS_DIR, all 90 2Mu2E samples, channels
genOnly and genOnly_trigger) holds the joint l_xy-vs-opening-angle histogram;
the reconstruction self-consistency run (V3_SELFCON_EOS_DIR, the 12
mid-lifetime m_Zd = 1.2 GeV samples, channel base) holds the gen-matched
lepton-jet pair masses. load_v3_cond() and load_v3_selfcon() load them the
same way, each with its own local cache.
"""
import os
import re
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.abspath(os.path.join(_HERE, "..", "..", ".."))):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import _lifetime_refit as lr

V3_EOS_DIR = ("/store/group/lpcmetx/SIDM/coffea_outputs/murtazas/"
              "truth_kinematics_forAN/anatomy_v3")
V3_CACHE = os.environ.get("ANATOMY_CACHE",
                          os.path.expanduser("~/nobackup/truthkin_v3_cache"))

V3_COND_EOS_DIR = ("/store/group/lpcmetx/SIDM/coffea_outputs/murtazas/"
                   "truth_kinematics_forAN/anatomy_v3_cond")
V3_COND_CACHE = os.environ.get(
    "ANATOMY_CACHE_COND",
    os.path.expanduser("~/nobackup/truthkin_v3_cond_cache"))

V3_SELFCON_EOS_DIR = ("/store/group/lpcmetx/SIDM/coffea_outputs/murtazas/"
                      "truth_kinematics_forAN/anatomy_v3_selfcon")
V3_SELFCON_CACHE = os.environ.get(
    "ANATOMY_CACHE_SELFCON",
    os.path.expanduser("~/nobackup/truthkin_v3_selfcon_cache"))

CHANNELS = ["genOnly", "genOnly_trigger", "baseNoLj_noTrigger", "baseNoLj", "base"]

_SAMPLE_RE = re.compile(
    r"(?P<mode>4Mu|2Mu2E)_(?P<mbs>[0-9p]+)GeV_(?P<mzd>[0-9p]+)GeV_(?P<ctau>[0-9p]+)mm")

def _f(s):
    return float(s.replace("p", "."))

def parse_sample(name):
    """'4Mu_200GeV_1p2GeV_0p48mm' -> dict(mode, mbs, mzd, ctau_mm)."""
    m = _SAMPLE_RE.match(name)
    if m is None:
        raise ValueError(f"unparseable sample name: {name}")
    return {"mode": m["mode"], "mbs": _f(m["mbs"]), "mzd": _f(m["mzd"]),
            "ctau_mm": _f(m["ctau"])}

def load_v3(cache=V3_CACHE, eos_dir=V3_EOS_DIR):
    return lr.load_truthkin(cache=cache, eos_dir=eos_dir)

def load_v3_cond(cache=V3_COND_CACHE, eos_dir=V3_COND_EOS_DIR):
    """The displacement-conditioned run: 90 2Mu2E samples, channels genOnly
    and genOnly_trigger, holding genA_toMu_lxy_vs_daughters_dR."""
    return lr.load_truthkin(cache=cache, eos_dir=eos_dir)

def load_v3_selfcon(cache=V3_SELFCON_CACHE, eos_dir=V3_SELFCON_EOS_DIR):
    """The reconstruction self-consistency run: 12 mid-lifetime
    m_Zd = 1.2 GeV samples, channel base, holding the gen-matched
    lepton-jet pair masses."""
    return lr.load_truthkin(cache=cache, eos_dir=eos_dir)

def format_sample(name):
    """Human-readable figure label for a raw sample name."""
    p = parse_sample(name)
    return (rf"{p['mode']}, $m_{{B_s}}$ = {p['mbs']:g} GeV, "
            rf"$m_{{Z_d}}$ = {p['mzd']:g} GeV, $c\tau$ = {p['ctau_mm']:g} mm")

def format_sample_2line(name):
    """Two-line sample label, for panels too narrow for the one-line form."""
    p = parse_sample(name)
    return (rf"{p['mode']}, $m_{{B_s}}$ = {p['mbs']:g} GeV" "\n"
            rf"$m_{{Z_d}}$ = {p['mzd']:g} GeV, $c\tau$ = {p['ctau_mm']:g} mm")

def get_h(out, sample, hist_name, channel="genOnly"):
    """One sample's histogram, channel axis already selected."""
    h = out[sample]["hists"][hist_name]
    h = getattr(h, "hist", h)
    return h[{"channel": channel}]

def sum_h(out, samples, hist_name, channel="genOnly"):
    """Sum a histogram over several samples (e.g. the 5 lifetimes of a mass point)."""
    hs = [get_h(out, s, hist_name, channel) for s in samples]
    total = hs[0].copy()
    for h in hs[1:]:
        total += h
    return total

def select_samples(out, mode=None, mbs=None, mzd=None, ctau_mm=None):
    """All sample names matching the given grid coordinates (None = any)."""
    keep = []
    for name in sorted(out):
        p = parse_sample(name)
        if mode is not None and p["mode"] != mode:
            continue
        if mbs is not None and p["mbs"] != mbs:
            continue
        if mzd is not None and p["mzd"] != mzd:
            continue
        if ctau_mm is not None and p["ctau_mm"] != ctau_mm:
            continue
        keep.append(name)
    return keep

def mid_ctau(out, mode, mbs, mzd):
    """The middle-lifetime sample of a mass point (5 lifetimes per point)."""
    names = select_samples(out, mode=mode, mbs=mbs, mzd=mzd)
    names.sort(key=lambda n: parse_sample(n)["ctau_mm"])
    return names[len(names) // 2]

def hist_median(h):
    """Flow-aware median, linearly interpolated inside the containing bin.

    The cumulative distribution includes the underflow and the total includes
    both flow bins, so a distribution with substantial overflow is not
    silently truncated: NaN is returned when the median lies beyond the last
    edge."""
    v = h.values(flow=True)
    edges = h.axes[-1].edges
    tot = v.sum()
    if tot == 0:
        return float("nan")
    below = np.cumsum(v)[:len(edges)]  # mass below each edge, underflow included
    if below[-1] < 0.5 * tot:
        return float("nan")  # median beyond the last edge
    i = int(np.searchsorted(below, 0.5 * tot))
    if i == 0:
        return float(edges[0])  # median inside the underflow
    width = edges[i] - edges[i - 1]
    return float(edges[i - 1] + (0.5 * tot - below[i - 1]) / v[i] * width)

def hist_median_tiered(h_low, h_full):
    """Median from the coarse full-range histogram, refined with the fine
    low-range histogram when it falls inside the full-range first bin. The
    refinement falls back automatically when the low-range histogram is
    itself overflow-dominated (its median comes back NaN)."""
    m_full = hist_median(h_full)
    if np.isnan(m_full) or m_full > h_full.axes[-1].edges[1]:
        return m_full
    m_low = hist_median(h_low)
    return m_full if np.isnan(m_low) else m_low

def hist_mean(h):
    vals = h.values()
    if vals.sum() == 0:
        return np.nan
    return float((h.axes[-1].centers * vals).sum() / vals.sum())

def counts_and_edges(h):
    return h.values(), h.axes[-1].edges

def efficiency(num_h, den_h):
    """Per-bin efficiency and Clopper-Pearson-style uncertainty from two
    same-binned count histograms (unweighted fills)."""
    import hist.intervals
    num, den = num_h.values(), den_h.values()
    ok = den > 0
    eff = np.full_like(den, np.nan, dtype=float)
    eff[ok] = num[ok] / den[ok]
    lo = np.full_like(eff, np.nan)
    hi = np.full_like(eff, np.nan)
    band = hist.intervals.clopper_pearson_interval(num[ok], den[ok])
    lo[ok], hi[ok] = band[0], band[1]
    return eff, lo, hi, num_h.axes[-1].centers

def event_count(out, sample, channel):
    """Event count of a channel from a one-entry-per-event histogram."""
    h = get_h(out, sample, "genMu_n", channel)
    return float(h.sum(flow=True).value if hasattr(h.sum(flow=True), "value")
                 else h.sum(flow=True))

def trigger_eff(out, sample):
    """Per-sample trigger efficiency: genOnly_trigger events / genOnly events."""
    den = event_count(out, sample, "genOnly")
    num = event_count(out, sample, "genOnly_trigger")
    if den == 0:
        return np.nan, np.nan
    err = np.sqrt(num * (1 - num / den)) / den if den > 0 else np.nan
    return num / den, err

MBS_VALUES = [100.0, 150.0, 200.0, 500.0, 800.0, 1000.0]
MZD_VALUES = [0.25, 1.2, 5.0]
