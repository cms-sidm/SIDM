"""Every figure this study makes, as reusable functions.

Replaces what had grown into six separate plotting scripts with overlapping
helpers -- three different ``draw_closure``s and two copies of the Lxy lookup,
which is how those two copies drifted apart once already.

The organising idea is that **nothing here is specific to a campaign, a method
or a plane.**  Functions take either

* ``rows`` -- the flat per-signal-point records from :func:`as_rows`, or
* ``sets`` -- an ordered ``{label: rows}`` mapping, for anything comparative, or
* ``grouped`` -- ``{process: {channel: {region: Shape}}}`` from
  :func:`shape_tools.group_shapes`, for the region and closure diagnostics,

so the same call draws one campaign, five methods, or two planes depending only
on what is handed to it.  ``sidm/scripts/abcd_plots.py`` is a thin CLI over this
module; the notebooks call it directly.

Conventions kept throughout:

* Brazil bands are the standard convention -- green 1 sigma inside yellow
  2 sigma, black dashed median.
* Multi-panel figures build inside :data:`GRID_RC`.  ``hep.style.CMS`` sets
  ``font.size = 26``, right for one full-size figure and illegible across a grid.
* Every figure is stamped with the campaign it came from, so a plot lifted into
  a talk still says what produced it.
* Ratio panels follow the data/MC convention: points carry the numerator's
  statistical error, the band carries the denominator's.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import fsspec
import hist
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import yaml

STUDY_DIR = Path(__file__).resolve().parent
REPO_ROOT = STUDY_DIR.parents[2]

# --------------------------------------------------------------------------- #
# Style
# --------------------------------------------------------------------------- #
GRID_RC = {"font.size": 13, "axes.labelsize": 13, "axes.titlesize": 13,
           "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10}
BAND_STYLE = {"1sigma": "#5BC236", "2sigma": "#F5DF4D"}
THEORY_COLOUR = "#C0392B"
OVERLAY_COLOUR = "#2B5FAD"      # a second set of limits drawn over a Brazil band
GROUP_COLOURS = {"QCD": "#4c72b0", "DY": "#dd8452", "TT": "#55a868",
                 "Diboson": "#c44e52", "other": "#8172b3"}
GROUP_LABELS = {"DY": "DY+jets", "TT": r"$t\bar{t}$+jets", "QCD": "QCD",
                "Diboson": "diboson"}
REGION_NAMES = {0: "A", 1: "B", 2: "C", 3: "D"}
REGION_COLOURS = {0: "#333333", 1: "#4c72b0", 2: "#dd8452", 3: "#55a868"}
# enough distinct styles for the five-method comparison without repeats
SERIES_STYLES = [("#444444", "o", "-"), ("#4c72b0", "s", "-"),
                 ("#C0392B", "^", "-"), ("#55a868", "D", "--"),
                 ("#8172b3", "v", "--"), ("#dd8452", "P", ":"),
                 ("#9db8d8", "X", ":"), ("#e39b95", "*", ":")]


def style_for(i):
    """``(colour, marker, linestyle)`` for the i-th series in a comparison."""
    return SERIES_STYLES[i % len(SERIES_STYLES)]


def _series_kwargs(i, label, styles=None, **defaults):
    """Plot kwargs for one series: an explicit ``styles[label]`` wins, else the cycle.

    The cycle is fine for a handful of sets.  Past that a comparison needs
    composite encoding -- colour for one dimension, marker and dash for another
    -- and the caller states it per label rather than inheriting cycled hues.
    """
    if styles and label in styles:
        return {**defaults, **styles[label]}
    colour, marker, ls = style_for(i)
    return {**defaults, "color": colour, "marker": marker, "ls": ls}


def _factored_legend(fig, legend_keys, y=0.0, fontsize=9):
    """One legend block per encoding channel, side by side under the figure.

    ``legend_keys`` is ``{title: [(label, line_kwargs), ...]}``.  With colour
    meaning one thing and dash/marker another, a key per channel identifies
    every line with far fewer entries than one row per series.
    """
    from matplotlib.lines import Line2D
    legs = []
    for title, keys in legend_keys.items():
        handles = [Line2D([], [], **kw) for _, kw in keys]
        legs.append(fig.legend(handles, [lab for lab, _ in keys], title=title or None,
                               frameon=False, fontsize=fontsize,
                               title_fontsize=fontsize + 1, loc="upper center",
                               bbox_to_anchor=(0.5, y), ncol=min(len(keys), 3),
                               handlelength=3.2))
    # The blocks are independent legends, so lay them out by hand: render once,
    # measure, then re-anchor them side by side about the centre.  Fixed slots
    # leave a wide gap on a wide figure.
    fig.canvas.draw()
    gap = 0.025
    widths = [leg.get_window_extent().width / fig.bbox.width for leg in legs]
    x = 0.5 - (sum(widths) + gap * (len(legs) - 1)) / 2
    for leg, w in zip(legs, widths):
        leg.set_bbox_to_anchor((x + w / 2, y), transform=fig.transFigure)
        x += w + gap


# Run 2 (2018) conditions, for the mplhep CMS label on every figure.
LUMI_FB = 59.83
COM_TEV = 13


def cms_label(ax, lumi=LUMI_FB, com=COM_TEV, title_pad=22, data=False, **kwargs):
    """The standard CMS Simulation label with the Run 2 luminosity.

    Almost everything here is simulation, hence ``data=False``, which makes
    mplhep write "Simulation" under "CMS"; a panel showing (approved) data
    passes ``data=True`` and a ``text`` such as "Work in progress".

    Pass the figure's axes: a single Axes, a list, or the 2-D array that
    ``subplots(squeeze=False)`` returns.  For a row of panels the experiment
    label goes over the top-left panel and the luminosity over the top-right
    one, which is the CMS convention for multi-panel figures and is also what
    keeps the two apart when the first panel is narrow.
    """
    arr = np.asarray(ax, dtype=object)
    if arr.ndim == 0:
        top = [arr.item()]
    elif arr.ndim == 1:
        top = list(arr)
    else:
        top = list(arr[0])
    first, last = top[0], top[-1]
    # mplhep writes the label just above the axes box, which is exactly where a
    # panel title already sits.  Push the top row's titles up so the two do not
    # overprint, all of them together so a row of headings stays on one line.
    for a in top:
        if a.get_title():
            a.set_title(a.get_title(), pad=title_pad,
                        fontsize=a.title.get_fontsize())
    if first is last:
        return hep.cms.label(ax=first, data=data, lumi=lumi, com=com, **kwargs)
    hep.cms.label(ax=first, data=data, rlabel="", **kwargs)
    size = {k: v for k, v in kwargs.items() if k == "fontsize"}
    return hep.label.exp_label(exp="", llabel="", lumi=lumi, com=com, ax=last, **size)


def to_hist(edges, sumw, sumw2=None):
    """A weighted ``hist.Hist`` from the arrays a ``Shape`` carries."""
    h = hist.Hist(hist.axis.Variable(np.asarray(edges, dtype=float), name="x"),
                  storage=hist.storage.Weight())
    view = h.view()
    view["value"] = np.asarray(sumw, dtype=float)
    view["variance"] = (np.asarray(sumw2, dtype=float) if sumw2 is not None
                        else np.asarray(sumw, dtype=float))
    return h


def stamp(fig, text):
    """Campaign label above the top-right corner.

    Outside the figure box rather than inside it: at y=0.998 it collides with
    the suptitle on compact figures, and ``bbox_inches="tight"`` pulls it back
    into the saved image anyway.
    """
    fig.text(0.998, 1.015, text, ha="right", va="bottom", fontsize=8, color="0.45")


def save(fig, outdir, name, campaign=""):
    """Write a figure as both PDF and PNG, stamped with its campaign."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if campaign:
        stamp(fig, campaign)
    for suffix in (".pdf", ".png"):
        fig.savefig(outdir / f"{name}{suffix}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return outdir / f"{name}.pdf"


def _log_ylim(ax, values, lo=0.3, hi=30):
    """Set a log y-range from the positive values only.

    One empty bin on a log axis otherwise drags the range down through fifteen
    decades of nothing.
    """
    positive = np.concatenate([np.asarray(v)[np.asarray(v) > 0] for v in values
                               if np.size(v)])
    if positive.size:
        ax.set_ylim(positive.min() * lo, positive.max() * hi)


def _step(ax, edges, values, **kwargs):
    """Step plot on bin edges, repeating the last value to close the final bin."""
    return ax.step(edges, np.append(values, values[-1]), where="post", **kwargs)


def _band(ax, edges, lo, hi, **kwargs):
    """Filled band on bin EDGES.

    Always ``step="post"`` on the edges, matching :func:`_step`.  Passing bin
    centres with ``step="mid"`` puts the transitions halfway between centres,
    which are not the edges once the bins have unequal widths -- that is what
    made bands sit visibly offset from their own lines.
    """
    kwargs.setdefault("linewidth", 0)
    return ax.fill_between(edges, np.append(lo, lo[-1]), np.append(hi, hi[-1]),
                           step="post", **kwargs)


# --------------------------------------------------------------------------- #
# The signal grid
# --------------------------------------------------------------------------- #
SIGNAL_NAME = re.compile(
    r"^(?P<fs>2Mu2E|4Mu|Comb)_(?P<mb>[\dp]+)GeV_(?P<mzd>[\dp]+)GeV_(?P<ctau>[\dp]+)mm$")

# (m_bound, m_ZD) -> the five generated lifetimes, mapping onto the five mean
# lab-frame Lxy values below.  Supplied with the analysis; not derived here.
_LXY_ROWS = [
    (100, 0.25, 0.02, 0.2, 2, 10, 20), (100, 1.2, 0.096, 0.96, 9.6, 48, 96),
    (100, 5, 0.4, 4, 40, 200, 400), (150, 0.25, 0.013, 0.13, 1.3, 6.7, 13),
    (150, 1.2, 0.064, 0.64, 6.4, 32, 64), (150, 5, 0.27, 2.7, 27, 130, 270),
    (200, 0.25, 0.01, 0.1, 1, 5, 10), (200, 1.2, 0.048, 0.48, 4.8, 24, 48),
    (200, 5, 0.2, 2, 20, 100, 200), (500, 0.25, 0.004, 0.04, 0.4, 2, 4),
    (500, 1.2, 0.019, 0.19, 1.9, 9.6, 19), (500, 5, 0.08, 0.8, 8, 40, 80),
    (800, 0.25, 0.0025, 0.025, 0.25, 1.2, 2.5), (800, 1.2, 0.012, 0.12, 1.2, 6, 12),
    (800, 5, 0.05, 0.5, 5, 25, 50), (1000, 0.25, 0.002, 0.02, 0.2, 1, 2),
    (1000, 1.2, 0.0096, 0.096, 0.96, 4.8, 9.6), (1000, 5, 0.04, 0.4, 4, 20, 40),
]
LXY_CM = [0.3, 3.0, 30.0, 150.0, 300.0]
LXY_LOOKUP = {(float(r[0]), float(r[1]), round(float(ct), 6)): lxy
              for r in _LXY_ROWS for ct, lxy in zip(r[2:], LXY_CM)}


def parse_point(name):
    """``{fs, m_bound, mzd, ctau, lxy, eps, eps2}`` for a signal point name."""
    m = SIGNAL_NAME.match(name)
    if not m:
        return None
    num = lambda x: float(x.replace("p", "."))
    mb, mzd, ctau = num(m["mb"]), num(m["mzd"]), num(m["ctau"])
    # epsilon = sqrt(80 / m_ZD / ctau) * 1e-6, with m_ZD in GeV and ctau in mm
    eps = np.sqrt(80.0 / mzd / ctau) * 1e-6 if mzd and ctau else np.nan
    return {"fs": m["fs"], "m_bound": mb, "mzd": mzd, "ctau": ctau,
            "lxy": LXY_LOOKUP.get((mb, mzd, round(ctau, 6)), np.nan),
            "eps": eps, "eps2": eps ** 2}


def theory_xs_fb(config="cross_sections.yaml"):
    """``{signal: sigma in fb}``.

    The values depend only on the bound state mass, which is why the model
    prediction is a horizontal line on an Lxy axis and a curve on a mass axis.
    """
    raw = yaml.safe_load((REPO_ROOT / "sidm" / "configs" / config).read_text())
    return {k: v * 1000.0 for k, v in raw.items()}          # pb -> fb


def theory_by_mass(xs=None):
    """``{m_bound: sigma in fb}`` collapsed over final state, m_ZD and lifetime."""
    xs = xs if xs is not None else theory_xs_fb()
    out = {}
    for name, value in xs.items():
        p = parse_point(name)
        if p:
            out[p["m_bound"]] = value
    return dict(sorted(out.items()))


# --------------------------------------------------------------------------- #
# Limits
# --------------------------------------------------------------------------- #
QUANTILES = ("exp_m2", "exp_m1", "exp", "exp_p1", "exp_p2")


def load_limits(path, with_bands=True):
    """Read a ``limits.csv`` into ``{signal: {quantile: value}}``.

    Rows with no expected median are dropped -- for this study that means a
    point with no signal at all in the region, where a limit is undefined
    rather than failed.

    Every value returned is a cross section limit in fb.  Rows labelled ``1fb``
    predate the normalisation in ``datacard_tools.read_coffea``: their cards held
    each signal at the cross section the production applied -- its model cross
    section where it has one, 1 fb otherwise -- so their r is relative to that
    and is converted here.  ``1fb_ref`` rows are at 1 fb already.

    ``path`` may also be a URL such as ``root://cmseos.fnal.gov//store/...``,
    read through fsspec as ``datacard_tools.read_coffea`` reads outputs, so
    limits on EOS can be plotted where EOS is not mounted.
    """
    production_fb = None
    out = {}
    with fsspec.open(str(path), "rt") as fh:
        for row in csv.DictReader(fh):
            if not row.get("exp"):
                continue
            k = 1.0
            if row.get("normalisation") == "1fb":
                production_fb = production_fb or theory_xs_fb()
                k = production_fb.get(row["signal"], 1.0)      # production xs / 1 fb
            entry = {"exp": float(row["exp"]) * k}
            if with_bands:
                entry.update({q: float(row[q]) * k for q in QUANTILES if row.get(q)})
            if row.get("obs"):
                entry["obs"] = float(row["obs"]) * k
            out[row["signal"]] = entry
    return out


def as_rows(limits, xs=None):
    """Flatten ``load_limits`` output into records with grid coordinates.

    ``limits`` may be a path or an already-loaded dict, so callers can hand
    either a CSV or something they built themselves.
    """
    if isinstance(limits, (str, Path)):
        limits = load_limits(limits)
    xs = xs if xs is not None else theory_xs_fb()
    rows = []
    for name, entry in limits.items():
        p = parse_point(name)
        if not p:
            continue
        row = dict(p, signal=name, **entry)
        row["theory_fb"] = xs.get(name, np.nan)
        row["r_theory"] = (entry["exp"] / row["theory_fb"]
                           if row["theory_fb"] == row["theory_fb"] else np.nan)
        rows.append(row)
    return rows


def load_sets(paths, xs=None):
    """``{label: rows}`` from ``{label: limits.csv path}``, skipping what is absent.

    This is what makes every comparison function generic: give it two campaigns,
    five methods, or one of each, and the drawing code does not change.
    """
    xs = xs if xs is not None else theory_xs_fb()
    out = {}
    for label, path in paths.items():
        path = Path(path)
        if path.is_dir():
            path = path / "limits.csv"
        if not path.exists():
            print(f"  (missing) {path}")
            continue
        out[label] = as_rows(load_limits(path), xs)
    return out


def select(rows, **cuts):
    """Rows matching every ``field=value`` cut, sorted by Lxy then bound mass."""
    out = [r for r in rows if all(r.get(k) == v for k, v in cuts.items())]
    return sorted(out, key=lambda r: (r["lxy"], r["m_bound"]))


def best_by(rows, key, fs=None):
    """Best (smallest) expected limit at each distinct value of ``key``."""
    best = {}
    for r in rows:
        if fs and r["fs"] != fs:
            continue
        v = r.get(key)
        if v != v:                                   # NaN
            continue
        if v not in best or r["exp"] < best[v]:
            best[v] = r["exp"]
    return best


def summary_table(sets):
    """Printed one-line-per-set summary: best limit, its point, exclusion count."""
    width = max((len(k) for k in sets), default=10) + 2
    print(f"{'set':<{width}} {'best [fb]':>10}  {'best point':<28} {'excluded':>9}")
    print("-" * (width + 52))
    for label, rows in sets.items():
        if not rows:
            continue
        best = min(rows, key=lambda r: r["exp"])
        n_ex = sum(1 for r in rows if r["r_theory"] == r["r_theory"] and r["r_theory"] < 1)
        n_th = sum(1 for r in rows if r["r_theory"] == r["r_theory"])
        print(f"{label:<{width}} {best['exp']:>10.4g}  {best['signal']:<28} "
              f"{n_ex:>4d}/{n_th:<4d}")


# --------------------------------------------------------------------------- #
# Single result set
# --------------------------------------------------------------------------- #
def draw_brazil(ax, rows, theory=None, label_median="median expected", observed=False,
                x="lxy", overlay=None, overlay_label=""):
    """Brazil bands against Lxy, or against the field ``x`` of the rows, for rows
    that share a panel.

    ``observed`` adds the observed limit as the usual solid black line, for
    results where one exists (the validation-region fits); it is drawn only if
    every row in the panel has one.  ``theory`` is the model cross section: a
    number, drawn as a horizontal line, or ``{x: sigma}``, drawn as a curve,
    which is what it is against bound state mass.  ``overlay`` is a second set
    of rows for the same panel, drawn over the bands (:func:`_draw_overlay`).
    """
    rows = sorted([r for r in rows if r[x] == r[x]], key=lambda r: r[x])
    if not rows:
        return False
    xs = np.array([r[x] for r in rows])
    if all("exp_m2" in r and "exp_p2" in r for r in rows):
        ax.fill_between(xs, [r["exp_m2"] for r in rows], [r["exp_p2"] for r in rows],
                        color=BAND_STYLE["2sigma"], label=r"expected $\pm2\sigma$")
        ax.fill_between(xs, [r["exp_m1"] for r in rows], [r["exp_p1"] for r in rows],
                        color=BAND_STYLE["1sigma"], label=r"expected $\pm1\sigma$")
    ax.plot(xs, [r["exp"] for r in rows], color="black", ls="--", lw=1.8,
            marker="o", ms=4, mfc="white" if observed else None, label=label_median)
    if observed and all("obs" in r for r in rows):
        ax.plot(xs, [r["obs"] for r in rows], color="black", ls="-", lw=1.8,
                marker="o", ms=5, label="observed")
    if overlay:
        _draw_overlay(ax, overlay, x, overlay_label)
    if isinstance(theory, dict):
        if theory:
            masses = sorted(theory)
            ax.plot(masses, [theory[m] for m in masses], color=THEORY_COLOUR, lw=1.8,
                    marker="s", ms=3.5, label=r"$\sigma_\mathrm{theory}$")
    elif theory is not None and theory == theory:
        ax.axhline(theory, color=THEORY_COLOUR, lw=1.8,
                   label=r"$\sigma_\mathrm{theory}$")
    ax.set_xscale("log" if x == "lxy" else "linear")
    ax.set_yscale("log")
    return True


def _draw_overlay(ax, rows, x, label=""):
    """A second set of limits over a Brazil band, in one colour.

    The median is a solid line under the band's dashed one, so where the two
    sets agree the dashes show on top of it, and the edges of the +-1 and +-2
    sigma bands are dotted lines, since a second filled band would hide the first.
    """
    rows = sorted([r for r in rows if r[x] == r[x]], key=lambda r: r[x])
    if not rows:
        return
    xs = [r[x] for r in rows]
    suffix = f", {label}" if label else ""
    # zorder 1.9: above the bands (1), below the lines of the main set (2)
    ax.plot(xs, [r["exp"] for r in rows], color=OVERLAY_COLOUR, lw=2.6, marker="s",
            ms=5, zorder=1.9, label="median expected" + suffix)
    edges = [q for q in ("exp_m2", "exp_m1", "exp_p1", "exp_p2")
             if all(q in r for r in rows)]
    for i, q in enumerate(edges):
        ax.plot(xs, [r[q] for r in rows], color=OVERLAY_COLOUR, lw=1.1, ls=":",
                zorder=1.9, label=(r"expected $\pm1\sigma$ and $\pm2\sigma$" + suffix
                                   if i == 0 else None))


def _data_label(data, text="Work in progress"):
    """``cms_label`` keywords: "Simulation" by default, ``text`` after "CMS" for data."""
    return {"data": True, "text": text} if data else {}


def _legend_from_fullest(fig, axes, ncol=4, y=-0.04, loc="lower center"):
    """One figure legend, taken from the panel that has the most entries.

    A fixed corner panel may hold no points (or no observed line), which would
    silently drop entries from the legend.  ``loc="upper center"`` hangs the
    legend below ``y`` instead of standing it on ``y``, for a legend of several rows.
    """
    best = max((a.get_legend_handles_labels() for a in np.ravel(axes)),
               key=lambda hl: len(hl[1]))
    fig.legend(*best, frameon=False, ncol=ncol, loc=loc,
               bbox_to_anchor=(0.5, y))


def brazil_grid(rows, fs, outdir, campaign, theory=None, name=None, title="",
                observed=False, data=False, sharey=True, label_median="median expected",
                data_text="Work in progress", against="lxy", overlay=None,
                overlay_label=""):
    """Brazil bands across the grid, one panel per (m_bound, m_ZD).

    ``observed`` adds the observed limits; ``data`` labels the figure as data
    ("Work in progress") rather than simulation.  ``sharey=False`` gives every
    panel its own range, which is what shows a factor of two between observed
    and expected when the limits span five decades across the grid.  Pass
    ``theory={}`` to leave out the model cross section.

    ``against="m_bound"`` draws the limits against bound state mass instead, one
    panel per (Lxy, m_ZD), with the model cross section as a curve.  ``overlay``
    is a second set of rows, labelled ``overlay_label``, drawn over the bands of
    every panel (:func:`draw_brazil`), to compare two sets of limits.
    """
    theory = theory if theory is not None else theory_by_mass()
    sub = [r for r in rows if r["fs"] == fs]
    over = [r for r in overlay or [] if r["fs"] == fs]
    column = "m_bound" if against == "lxy" else "lxy"     # what changes along a row
    columns = sorted({r[column] for r in sub if r[column] == r[column]})
    masses = sorted({r["m_bound"] for r in sub})
    mzds = sorted({r["mzd"] for r in sub})
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(len(mzds), len(columns), sharex=True, sharey=sharey,
                                 figsize=(3.1 * len(columns), 2.5 * len(mzds)),
                                 squeeze=False, constrained_layout=True)
        for i, mzd in enumerate(mzds):
            for j, value in enumerate(columns):
                ax = axes[i][j]
                panel = [r for r in sub if r["mzd"] == mzd and r[column] == value]
                panel_overlay = [r for r in over if r["mzd"] == mzd and r[column] == value]
                if against == "lxy":
                    panel_theory = theory.get(value, np.nan)
                else:
                    panel_theory = {m: theory[m] for m in masses if m in theory}
                if not draw_brazil(ax, panel, panel_theory, label_median=label_median,
                                   observed=observed, x=against, overlay=panel_overlay,
                                   overlay_label=overlay_label):
                    ax.text(0.5, 0.5, "no points", transform=ax.transAxes,
                            ha="center", va="center", color="0.6")
                if i == 0:
                    ax.set_title(rf"$m_B={value:g}$ GeV" if against == "lxy"
                                 else rf"$L_{{xy}}={value:g}$ cm")
                if j == 0:
                    ax.set_ylabel(rf"$m_{{Z_d}}={mzd:g}$ GeV" + "\n"
                                  + r"$\sigma$ limit [fb]")
                if i == len(mzds) - 1:
                    ax.set_xlabel(r"$L_{xy}$ [cm]" if against == "lxy" else r"$m_B$ [GeV]")
                if against != "lxy":
                    ax.set_xticks(masses)
                ax.grid(alpha=0.25, which="both")
        if over:
            # one column for the set with the bands and one for the overlay, below the panels
            _legend_from_fullest(fig, axes, ncol=2, y=0.0, loc="upper center")
        else:
            _legend_from_fullest(fig, axes)
        fig.suptitle(title or f"{fs}: expected 95% CL limit (MC only, blinded)")
        cms_label(axes, fontsize=11, **_data_label(data, data_text))
        default = f"brazil_grid_{fs}" if against == "lxy" else f"brazil_grid_vs_mass_{fs}"
        return save(fig, outdir, name or default, campaign)


def compare_limits_grid(sets, fs, outdir, campaign, name, title="", colours=None,
                        observed=True, band=True, data=False, sharey=False, theory=None,
                        set_title="fit to", data_text="Work in progress"):
    """Several limit sets on one Brazil-style grid, one panel per (m_bound, m_ZD).

    Each set gets one colour: its median expected limit dashed (open markers),
    its +-1 sigma expected band shaded, and its observed limit solid (filled
    markers).  For comparing the same fit run on two inputs, such as data and
    simulation in a validation region.
    """
    sub = {k: [r for r in rows if r["fs"] == fs and r["lxy"] == r["lxy"]]
           for k, rows in sets.items()}
    masses = sorted({r["m_bound"] for rows in sub.values() for r in rows})
    mzds = sorted({r["mzd"] for rows in sub.values() for r in rows})
    colours = colours or {k: style_for(i)[0] for i, k in enumerate(sets)}
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(len(mzds), len(masses), sharex=True, sharey=sharey,
                                 figsize=(3.1 * len(masses), 2.5 * len(mzds)),
                                 squeeze=False, constrained_layout=True)
        for i, mzd in enumerate(mzds):
            for j, mb in enumerate(masses):
                ax = axes[i][j]
                drawn = False
                for label, rows in sub.items():
                    pts = sorted([r for r in rows if r["mzd"] == mzd and r["m_bound"] == mb],
                                 key=lambda r: r["lxy"])
                    if not pts:
                        continue
                    drawn = True
                    x, c = [r["lxy"] for r in pts], colours[label]
                    if band and all("exp_m1" in r and "exp_p1" in r for r in pts):
                        ax.fill_between(x, [r["exp_m1"] for r in pts],
                                        [r["exp_p1"] for r in pts], color=c, alpha=0.16, lw=0)
                    ax.plot(x, [r["exp"] for r in pts], color=c, ls="--", lw=1.5,
                            marker="o", ms=4, mfc="white")
                    if observed and all("obs" in r for r in pts):
                        ax.plot(x, [r["obs"] for r in pts], color=c, ls="-", lw=1.7,
                                marker="o", ms=4.5)
                if not drawn:
                    ax.text(0.5, 0.5, "no points", transform=ax.transAxes,
                            ha="center", va="center", color="0.6")
                if theory and theory.get(mb):
                    ax.axhline(theory[mb], color=THEORY_COLOUR, lw=1.6)
                ax.set_xscale("log")
                ax.set_yscale("log")
                if i == 0:
                    ax.set_title(rf"$m_B={mb:g}$ GeV")
                if j == 0:
                    ax.set_ylabel(rf"$m_{{Z_d}}={mzd:g}$ GeV" + "\n"
                                  + r"$\sigma$ limit [fb]")
                if i == len(mzds) - 1:
                    ax.set_xlabel(r"$L_{xy}$ [cm]")
                ax.grid(alpha=0.25, which="both")
        kinds = [("median expected", dict(color="0.3", ls="--", lw=1.5, marker="o", ms=4,
                                          mfc="white"))]
        if observed:
            kinds.append(("observed", dict(color="0.3", ls="-", lw=1.7, marker="o", ms=4.5)))
        if theory:
            kinds.append((r"$\sigma_\mathrm{theory}$", dict(color=THEORY_COLOUR, lw=1.6)))
        _factored_legend(fig, {set_title: [(k, dict(color=c, lw=6, alpha=0.5 if band else 1))
                                           for k, c in colours.items()],
                               "limit" + (r" (shaded: expected $\pm1\sigma$)" if band else ""):
                                   kinds}, y=-0.005)
        fig.suptitle(title)
        cms_label(axes, fontsize=11, **_data_label(data, data_text))
        return save(fig, outdir, name, campaign)


def limit_vs_mass(rows, fs, outdir, campaign, theory=None, name=None, title="",
                  data=False, data_text="Work in progress"):
    """Median limit against bound state mass, with the model cross section."""
    theory = theory if theory is not None else theory_by_mass()
    sub = [r for r in rows if r["fs"] == fs and r["lxy"] == r["lxy"]]
    lxys = sorted({r["lxy"] for r in sub})
    mzds = sorted({r["mzd"] for r in sub})
    colours = dict(zip(lxys, plt.cm.viridis(np.linspace(0, 0.9, len(lxys)))))
    styles = dict(zip(mzds, ["-", "--", ":", "-."]))
    with plt.rc_context(GRID_RC):
        fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
        # One curve per (m_ZD, Lxy): grouping on Lxy alone merges the m_ZD
        # values into a single polyline and zigzags between them.
        for mzd in mzds:
            for lxy in lxys:
                pts = sorted([r for r in sub if r["mzd"] == mzd and r["lxy"] == lxy],
                             key=lambda r: r["m_bound"])
                if len(pts) < 2:
                    continue
                ax.plot([p["m_bound"] for p in pts], [p["exp"] for p in pts],
                        color=colours[lxy], ls=styles[mzd], marker="o", ms=4, lw=1.4)
        if theory:
            ax.plot(list(theory), list(theory.values()), color=THEORY_COLOUR,
                    lw=2.4, marker="s", ms=6, label=r"$\sigma_\mathrm{theory}$")
        for lxy, colour in colours.items():
            ax.plot([], [], color=colour, marker="o", ms=4,
                    label=rf"$L_{{xy}}={lxy:g}$ cm")
        for mzd, style in styles.items():
            ax.plot([], [], color="0.35", ls=style, label=rf"$m_{{Z_d}}={mzd:g}$ GeV")
        ax.set_yscale("log")
        ax.set_xlabel(r"bound state mass $m_B$ [GeV]")
        ax.set_ylabel(r"median expected 95% CL limit on $\sigma$ [fb]")
        ax.set_title(title or fs)
        ax.grid(alpha=0.25, which="both")
        ax.legend(frameon=False, ncol=2, fontsize=9)
        cms_label(ax, **_data_label(data, data_text))
        return save(fig, outdir, name or f"limit_vs_mass_{fs}", campaign)


def exclusion_map(rows, fs, outdir, campaign, name=None, title="", field="r_theory",
                  label=r"$\sigma_\mathrm{limit}/\sigma_\mathrm{theory}$", outline_below=1.0):
    """sigma_limit / sigma_theory across the grid; cells below 1 are outlined.

    ``field`` maps any other value of the rows instead, such as the size of a
    systematic uncertainty, with ``label`` on the colour bar; ``outline_below=None``
    outlines nothing.
    """
    sub = [r for r in rows if r["fs"] == fs and r[field] == r[field]]
    if not sub:
        return None
    masses = sorted({r["m_bound"] for r in sub})
    mzds = sorted({r["mzd"] for r in sub})
    lxys = sorted({r["lxy"] for r in sub})
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, len(mzds), figsize=(4.4 * len(mzds), 3.6),
                                 squeeze=False, constrained_layout=True)
        grid = np.full((len(mzds), len(lxys), len(masses)), np.nan)
        for k, mzd in enumerate(mzds):
            for i, lxy in enumerate(lxys):
                for j, mb in enumerate(masses):
                    hit = [r[field] for r in sub if r["mzd"] == mzd
                           and r["lxy"] == lxy and r["m_bound"] == mb]
                    if hit:
                        grid[k, i, j] = min(hit)
        vmin, vmax = np.nanmin(grid), np.nanmax(grid)
        norm = matplotlib.colors.LogNorm(vmin=max(vmin, 1e-3), vmax=vmax)
        for k, mzd in enumerate(mzds):
            ax = axes[0][k]
            im = ax.imshow(grid[k], origin="lower", aspect="auto",
                           cmap="viridis_r", norm=norm)
            for i in range(len(lxys)):
                for j in range(len(masses)):
                    v = grid[k, i, j]
                    if v != v:
                        continue
                    ax.text(j, i, f"{v:.2g}", ha="center", va="center", fontsize=8,
                            color="white" if v > np.sqrt(vmin * vmax) else "black")
                    if outline_below is not None and v < outline_below:
                        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                                   fill=False, edgecolor="red", lw=2))
            ax.set_xticks(range(len(masses)))
            ax.set_xticklabels([f"{m:g}" for m in masses])
            ax.set_yticks(range(len(lxys)))
            ax.set_yticklabels([f"{l:g}" for l in lxys])
            ax.set_xlabel(r"$m_B$ [GeV]")
            if k == 0:
                ax.set_ylabel(r"$L_{xy}$ [cm]")
            ax.set_title(rf"$m_{{Z_d}}={mzd:g}$ GeV")
        fig.colorbar(im, ax=axes[0].tolist(), label=label)
        if not title and field == "r_theory":
            n_excl = int(np.nansum(grid < 1))
            title = f"{fs}: expected exclusion (outlined: below 1, {n_excl} cells)"
        fig.suptitle(title or f"{fs}: {label}")
        cms_label(axes, fontsize=11)
        return save(fig, outdir, name or f"exclusion_map_{fs}", campaign)


def kinetic_mixing(rows, outdir, campaign, name="eps2_vs_mdp", title="", ncols=None,
                   data=False, data_text="Work in progress"):
    """Expected limit in the eps^2 vs dark photon mass plane, one panel per mass.

    ``ncols`` wraps the panels onto several rows; the default is one row, which
    gets cramped past four bound state masses.
    """
    sub = [r for r in rows if r["eps2"] == r["eps2"]]
    masses = sorted({r["m_bound"] for r in sub})
    n = len(masses)
    ncols = min(ncols or n, n)
    nrows = -(-n // ncols)
    with plt.rc_context(GRID_RC):
        size = ((2.7 * n, 3.6) if nrows == 1
                else (3.0 * ncols + 1.0, 2.9 * nrows + 0.7))
        fig, axes = plt.subplots(nrows, ncols, figsize=size, sharex=True, sharey=True,
                                 squeeze=False, constrained_layout=True)
        vals = [r["exp"] for r in sub]
        norm = matplotlib.colors.LogNorm(vmin=min(vals), vmax=max(vals))
        for k in range(n, nrows * ncols):
            axes[k // ncols][k % ncols].set_visible(False)
        for j, mb in enumerate(masses):
            ax = axes[j // ncols][j % ncols]
            panel = [r for r in sub if r["m_bound"] == mb]
            sc = ax.scatter([r["mzd"] for r in panel], [r["eps2"] for r in panel],
                            c=[r["exp"] for r in panel], norm=norm, cmap="viridis_r",
                            s=42, edgecolor="0.3", linewidth=0.4)
            ax.set_xscale("log")
            ax.set_yscale("log")
            if j + ncols >= n:              # nothing below this panel
                ax.set_xlabel(r"$m_{Z_d}$ [GeV]")
                ax.tick_params(labelbottom=True)
            if j % ncols == 0:
                ax.set_ylabel(r"$\epsilon^2$")
            ax.set_title(rf"$m_B={mb:g}$ GeV")
            ax.grid(alpha=0.25, which="both")
        fig.colorbar(sc, ax=[a for a in axes.ravel() if a.get_visible()],
                     label=r"expected limit on $\sigma$ [fb]")
        fig.suptitle(title or r"Kinetic mixing plane "
                              r"($\epsilon=\sqrt{80/m_{Z_d}c\tau}\times10^{-6}$)")
        cms_label(axes, fontsize=11, **_data_label(data, data_text))
        return save(fig, outdir, name, campaign)


def result_set(rows, outdir, campaign, theory=None, prefix="", ncols=None):
    """The whole single-campaign result suite, for both final states."""
    theory = theory if theory is not None else theory_by_mass()
    made = []
    for fs in sorted({r["fs"] for r in rows}):
        made.append(brazil_grid(rows, fs, outdir, campaign, theory,
                                name=f"{prefix}brazil_grid_{fs}"))
        made.append(limit_vs_mass(rows, fs, outdir, campaign, theory,
                                  name=f"{prefix}limit_vs_mass_{fs}"))
        got = exclusion_map(rows, fs, outdir, campaign, name=f"{prefix}exclusion_map_{fs}")
        if got:
            made.append(got)
    made.append(kinetic_mixing(rows, outdir, campaign, name=f"{prefix}eps2_vs_mdp",
                               ncols=ncols))
    return [m for m in made if m]


# --------------------------------------------------------------------------- #
# Comparisons -- any number of sets, whatever they are
# --------------------------------------------------------------------------- #
def compare_vs_lxy(sets, fs, outdir, campaign, theory=None, name=None, title="",
                   bands=False, styles=None, legend_keys=None, ncols=None):
    """Best limit at each lifetime for every set, one panel per bound state mass.

    ``styles``/``legend_keys`` switch to composite encoding with a factored
    legend; see :func:`_series_kwargs` and :func:`_factored_legend`.
    """
    theory = theory if theory is not None else theory_by_mass()
    masses = sorted({r["m_bound"] for rows in sets.values() for r in rows
                     if r["fs"] == fs})
    with plt.rc_context(GRID_RC):
        n = len(masses)
        ncols = min(ncols or n, n)
        nrows = -(-n // ncols)
        size = (2.9 * n, 4.0) if nrows == 1 else (3.7 * ncols, 3.0 * nrows + 0.5)
        fig, axes = plt.subplots(nrows, ncols, figsize=size, sharex=True, sharey=True,
                                 squeeze=False, constrained_layout=True)
        for k in range(n, nrows * ncols):
            axes[k // ncols][k % ncols].set_visible(False)
        for j, mb in enumerate(masses):
            ax = axes[j // ncols][j % ncols]
            for i, (label, rows) in enumerate(sets.items()):
                kw = _series_kwargs(i, label, styles, ms=5, lw=1.7)
                sub = [r for r in rows if r["m_bound"] == mb and r["fs"] == fs]
                best = {}
                for r in sub:
                    if r["lxy"] != r["lxy"]:
                        continue
                    if r["lxy"] not in best or r["exp"] < best[r["lxy"]]["exp"]:
                        best[r["lxy"]] = r
                if not best:
                    continue
                xs_ = sorted(best)
                ax.plot(xs_, [best[x]["exp"] for x in xs_], label=label, **kw)
                if bands and all("exp_m1" in best[x] for x in xs_):
                    ax.fill_between(xs_, [best[x]["exp_m1"] for x in xs_],
                                    [best[x]["exp_p1"] for x in xs_],
                                    color=kw["color"], alpha=0.18)
            if theory.get(mb):
                ax.axhline(theory[mb], color="0.3", ls=(0, (1, 1)), lw=1.6)
            ax.set_xscale("log")
            ax.set_yscale("log")
            if j + ncols >= n:              # nothing below this panel
                ax.set_xlabel(r"$L_{xy}$ [cm]")
                ax.tick_params(labelbottom=True)
            ax.set_title(rf"$m_B={mb:g}$ GeV")
            ax.grid(alpha=0.25, which="both")
            if j % ncols == 0:
                ax.set_ylabel(r"best expected $\sigma$ limit [fb]")
        if legend_keys:
            _factored_legend(fig, {**legend_keys, "": [
                ("model cross section", dict(color="0.3", ls=(0, (1, 1)), lw=1.6))]},
                y=-0.02)
        else:
            handles, labels = axes[0][0].get_legend_handles_labels()
            fig.legend(handles, labels, frameon=False, ncol=3, loc="lower center",
                       bbox_to_anchor=(0.5, -0.13), fontsize=9)
        fig.suptitle(title or f"{fs}: dotted is the model cross section")
        cms_label(axes, fontsize=11)
        return save(fig, outdir, name or f"compare_vs_lxy_{fs}", campaign)


def compare_vs_mass(sets, outdir, campaign, theory=None, name="compare_vs_mass",
                    title="", styles=None, legend_keys=None):
    """Best limit per bound state mass for every set, one panel per final state."""
    theory = theory if theory is not None else theory_by_mass()
    states = sorted({r["fs"] for rows in sets.values() for r in rows})
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, len(states), figsize=(5.5 * len(states), 4.4),
                                 squeeze=False, constrained_layout=True)
        for ax, fs in zip(axes[0], states):
            for i, (label, rows) in enumerate(sets.items()):
                kw = _series_kwargs(i, label, styles, ms=6, lw=1.8)
                best = best_by(rows, "m_bound", fs)
                if not best:
                    continue
                ms = sorted(best)
                ax.plot(ms, [best[m] for m in ms], label=label, **kw)
            if theory:
                # under composite encoding the reds and oranges belong to the
                # methods, so the model line goes neutral, as on the Lxy plots
                tkw = (dict(color="0.25", lw=1.6, ls=(0, (1, 1))) if legend_keys else
                       dict(color=THEORY_COLOUR, lw=2.0, ls="--", marker="^", ms=5))
                ax.plot(list(theory), list(theory.values()),
                        label=r"$\sigma_\mathrm{theory}$", **tkw)
            ax.set_yscale("log")
            ax.set_xlabel(r"$m_B$ [GeV]")
            ax.set_ylabel(r"best expected $\sigma$ limit [fb]")
            ax.set_title(fs)
            ax.grid(alpha=0.25, which="both")
        if legend_keys:
            keys = dict(legend_keys)
            if theory:
                keys = {**keys, "": [("model cross section",
                                      dict(color="0.25", lw=1.6, ls=(0, (1, 1))))]}
            _factored_legend(fig, keys, y=-0.02)
        else:
            axes[0][0].legend(frameon=False, fontsize=9)
        fig.suptitle(title or "Best expected limit per bound state mass")
        cms_label(axes, fontsize=12)
        return save(fig, outdir, name, campaign)


def compare_ratio(sets, reference, outdir, campaign, against="lxy",
                  name=None, title="", styles=None):
    """Every set divided by a reference set, against Lxy or bound state mass.

    ``against="m_bound"`` is how to ask whether a method gains more at high
    mass; ``against="lxy"`` shows whether the gain depends on lifetime (it
    should not, since the background is the same at every lifetime and only the
    signal acceptance moves, so it cancels).
    """
    if reference not in sets:
        raise KeyError(f"reference {reference!r} not among {sorted(sets)}")
    ref = {r["signal"]: r["exp"] for r in sets[reference]}
    states = sorted({r["fs"] for rows in sets.values() for r in rows})
    label_x = {"lxy": r"$L_{xy}$ [cm]", "m_bound": r"$m_B$ [GeV]"}[against]
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, len(states), figsize=(5.8 * len(states), 4.4),
                                 sharey=True, squeeze=False, constrained_layout=True)
        for ax, fs in zip(axes[0], states):
            for i, (label, rows) in enumerate(sets.items()):
                if label == reference:
                    continue
                kw = _series_kwargs(i, label, styles, ms=5, lw=1.7)
                by_x = {}
                for r in rows:
                    if r["fs"] != fs or r["signal"] not in ref or not ref[r["signal"]]:
                        continue
                    x = r.get(against)
                    if x != x:
                        continue
                    by_x.setdefault(x, []).append(r["exp"] / ref[r["signal"]])
                if not by_x:
                    continue
                xs_ = sorted(by_x)
                allv = [v for vs in by_x.values() for v in vs]
                ax.plot(xs_, [float(np.median(by_x[x])) for x in xs_],
                        label=f"{label}  (median {np.median(allv):.2f})", **kw)
            ax.axhline(1.0, color="0.3", lw=1.2, ls="--")
            if against == "lxy":
                ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(label_x)
            ax.set_title(fs)
            ax.grid(alpha=0.25, which="both")
        axes[0][0].set_ylabel(f"limit / {reference} limit")
        axes[0][0].legend(frameon=False, fontsize=8)
        fig.suptitle(title or f"Relative to {reference} (below 1 = stronger)")
        cms_label(axes, fontsize=12)
        return save(fig, outdir, name or f"compare_ratio_vs_{against}", campaign)


def compare_scatter(sets, x_label, y_label, outdir, campaign,
                    name="compare_scatter", title=""):
    """Point-by-point scatter of two sets against each other."""
    xr = {r["signal"]: r["exp"] for r in sets[x_label]}
    yr = {r["signal"]: r["exp"] for r in sets[y_label]}
    shared = sorted(set(xr) & set(yr))
    with plt.rc_context(GRID_RC):
        fig, ax = plt.subplots(figsize=(5.6, 5.2), constrained_layout=True)
        for i, fs in enumerate(sorted({(parse_point(s) or {}).get("fs")
                                       for s in shared} - {None})):
            colour, marker, _ = style_for(i + 1)
            pts = [(xr[s], yr[s]) for s in shared
                   if (parse_point(s) or {}).get("fs") == fs]
            if pts:
                ax.scatter(*zip(*pts), s=28, color=colour, marker=marker, alpha=0.75,
                           label=f"{fs} ({len(pts)} points)")
        lims = [min(min(xr[s], yr[s]) for s in shared) * 0.6,
                max(max(xr[s], yr[s]) for s in shared) * 1.6]
        ax.plot(lims, lims, color="0.4", ls="--", lw=1.2, label="equal")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel(f"{x_label} expected limit [fb]")
        ax.set_ylabel(f"{y_label} expected limit [fb]")
        ax.set_title(title or f"Above the line: {y_label} is weaker")
        ax.grid(alpha=0.25, which="both")
        ax.legend(frameon=False, fontsize=9)
        cms_label(ax)
        return save(fig, outdir, name, campaign)


def compare_exclusion(sets, outdir, campaign, name="compare_exclusion",
                      per_mass=False, title="", styles=None, legend_keys=None,
                      groups=None, emphasise=()):
    """Points reaching expected exclusion, per set (optionally split by mass).

    With ``styles`` the chart is horizontal, one bar per set coloured by its
    style, and carries no number on each bar: with many sets the labels crowd,
    and the exact counts belong in a table beside the figure.
    """
    with plt.rc_context(GRID_RC):
        if styles:
            labels = list(sets)
            totals = [sum(1 for r in sets[k] if r["r_theory"] == r["r_theory"]
                          and r["r_theory"] < 1) for k in labels]
            fig, ax = plt.subplots(figsize=(7.6, 0.34 * len(labels) + 1.6),
                                   constrained_layout=True)
            ypos = np.arange(len(labels))[::-1]
            for y, k, n in zip(ypos, labels, totals):
                st_ = styles.get(k, {})
                ax.barh(y, n, height=0.62, color=st_.get("color", "0.5"),
                        edgecolor=st_.get("color", "0.5"),
                        hatch=st_.get("hatch"), linewidth=0.8, alpha=0.9)
            ax.set_yticks(ypos)
            ax.set_yticklabels(labels, fontsize=9)
            for tick, k in zip(ax.get_yticklabels(), labels):
                if k in emphasise:
                    tick.set_fontweight("bold")
            if groups:
                # a hairline wherever the campaign changes, so the bars read
                # in blocks rather than as one long list
                for y, (a, b) in zip(ypos[1:], zip(labels, labels[1:])):
                    if groups.get(a) != groups.get(b):
                        ax.axhline(y + 0.5, color="0.75", lw=0.8)
            ax.set_xlabel("points expected to be excluded")
            ax.grid(alpha=0.25, axis="x")
            ax.set_axisbelow(True)
            if legend_keys:
                _factored_legend(fig, legend_keys, y=-0.02)
            ax.set_title(title or "Expected exclusion")
            cms_label(ax)
            return save(fig, outdir, name, campaign)
        if per_mass:
            masses = sorted({r["m_bound"] for rows in sets.values() for r in rows
                             if r["r_theory"] == r["r_theory"]})
            fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
            width = 0.8 / max(len(sets), 1)
            for i, (label, rows) in enumerate(sets.items()):
                counts = [sum(1 for r in rows if r["m_bound"] == m
                              and r["r_theory"] == r["r_theory"] and r["r_theory"] < 1)
                          for m in masses]
                xs_ = np.arange(len(masses)) + (i - (len(sets) - 1) / 2) * width
                ax.bar(xs_, counts, width, color=style_for(i)[0], alpha=0.9,
                       label=f"{label}  ({sum(counts)})")
                for x, c in zip(xs_, counts):
                    if c:
                        ax.text(x, c + 0.3, str(c), ha="center", fontsize=8)
            ax.set_xticks(range(len(masses)))
            ax.set_xticklabels([f"{m:g}" for m in masses])
            ax.set_xlabel(r"$m_B$ [GeV]")
            ax.legend(frameon=False, fontsize=9)
        else:
            fig, ax = plt.subplots(figsize=(8.4, 4.2), constrained_layout=True)
            labels = list(sets)
            totals = [sum(1 for r in sets[k] if r["r_theory"] == r["r_theory"]
                          and r["r_theory"] < 1) for k in labels]
            bars = ax.bar(range(len(labels)), totals,
                          color=[style_for(i)[0] for i in range(len(labels))], alpha=0.9)
            for b, n in zip(bars, totals):
                ax.text(b.get_x() + b.get_width() / 2, n + 0.4, str(n),
                        ha="center", fontsize=10)
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=18, ha="right", fontsize=9)
        ax.set_ylabel("points expected to be excluded")
        ax.set_title(title or "Expected exclusion")
        cms_label(ax)
        return save(fig, outdir, name, campaign)


# --------------------------------------------------------------------------- #
# Region and closure diagnostics
# --------------------------------------------------------------------------- #
def region_arrays(grouped, channel, n_regions=4):
    """``(edges, {region: (sumw, err)}, prediction, prediction_err)``.

    The prediction error propagates the three control-region statistical errors
    through ``B*C/D`` as independent, which they are: disjoint sets of events.
    """
    import shape_tools as st
    shapes = {r: st.total_background_shape(grouped, channel, r)
              for r in range(n_regions)}
    edges = shapes[0].edges
    w = {r: shapes[r].sumw for r in shapes}
    e = {r: np.sqrt(shapes[r].sumw2) for r in shapes}
    with np.errstate(divide="ignore", invalid="ignore"):
        pred = np.where(w[3] > 0, w[1] * w[2] / w[3], np.nan)
        rel = np.sqrt(sum(np.where(w[k] > 0, (e[k] / w[k]) ** 2, 0) for k in (1, 2, 3)))
    return edges, {r: (w[r], e[r]) for r in shapes}, pred, pred * rel


def per_width(edges, sumw, sumw2=None, per=None):
    """Bin contents per ``per`` units of the axis, variances to match.

    For bins of unequal width, so a smoothly falling distribution looks smooth.
    ``per=None`` returns the contents unchanged.
    """
    if per is None:
        return sumw, sumw2
    f = per / np.diff(np.asarray(edges, dtype=float))
    return sumw * f, (None if sumw2 is None else sumw2 * f * f)


def region_distributions(grouped, channel, outdir, campaign, signal=None,
                         label=r"$m_{LJLJ}$ [GeV]", name=None, title="",
                         region_a_label="region A  (signal region)", legend_outside=False,
                         per=None, data=None):
    """The observable in all four ABCD regions, stacked by process group.

    ``per`` (e.g. 100) draws events per that many GeV instead of per bin.
    ``data`` is an optional ``{region: Shape}`` drawn as points; pass only what
    the analysis has approved looking at.

    Built from ``hist.Hist`` objects and drawn with ``mplhep.histplot``, so the
    stacking, the error band and the CMS label all follow the standard
    conventions.  ``signal`` is an optional ``{region: Shape}`` to overlay.

    Log scale, so a bin with no simulated background simply draws nothing --
    which is the thing worth seeing, and is left to the eye rather than
    annotated on the figure.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    # The legend sits in region A, but has to name every group drawn in any of
    # the four panels: a process present only in a control region would
    # otherwise be a coloured block with no entry.
    present = [g for g in sorted(grouped) if any(
        (grouped[g].get(channel, {}).get(r) is not None
         and grouped[g][channel][r].sumw.sum() > 0) for r in range(4))]
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, 4, figsize=(16, 4.0), sharey=True,
                                 constrained_layout=True)
        drew_signal = False
        for ax, region in zip(axes, range(4)):
            stack, labels, colours, total = [], [], [], None
            for group in present:
                shape = grouped[group].get(channel, {}).get(region)
                if shape is None or shape.sumw.sum() <= 0:
                    continue
                h = to_hist(shape.edges, *per_width(shape.edges, shape.sumw,
                                                    shape.sumw2, per))
                stack.append(h)
                labels.append(GROUP_LABELS.get(group, group))
                colours.append(GROUP_COLOURS.get(group, "0.6"))
                total = h if total is None else total + h
            if total is None:
                ax.set_visible(False)
                continue

            hep.histplot(stack, ax=ax, stack=True, histtype="fill",
                         label=labels, color=colours,
                         edgecolor="white", linewidth=0.4)
            hep.histplot(total, ax=ax, histtype="errorbar", yerr=True,
                         color="0.25", elinewidth=1.2, capsize=2,
                         markersize=0, label=None)

            sig = (signal or {}).get(region)
            top = float(total.view()["value"].max())
            obs = (data or {}).get(region)
            if obs is not None:
                ow, ow2 = per_width(obs.edges, obs.sumw, obs.sumw2, per)
                hep.histplot(to_hist(obs.edges, ow, ow2), ax=ax, histtype="errorbar",
                             yerr=np.sqrt(ow2), color="black", marker="o", markersize=3.5,
                             elinewidth=1.1, label="data")
                top = max(top, float(ow.max()))
            if sig is not None and sig.sumw.sum() > 0:
                sig_w = per_width(sig.edges, sig.sumw, per=per)[0]
                hep.histplot(to_hist(sig.edges, sig_w), ax=ax, yerr=False,
                             histtype="step", color="crimson", linewidth=1.8,
                             label="signal")
                top = max(top, float(sig_w.max()))
                drew_signal = True

            ax.set_yscale("log")
            ax.set_ylim(1e-3, max(top, 1e-2) * 10 ** 2.0)
            ax.set_xlabel(label)
            ax.set_ylabel((f"events / {per:g} GeV" if per else "events")
                          if ax is axes[0] else "")
            # a validation region has its own A, which is not the signal region
            name_r = region_a_label if region == 0 else f"region {REGION_NAMES[region]}"
            ax.legend([], [], title=name_r, frameon=False, title_fontsize=11,
                      loc="upper right", alignment="right")
        handles = [Patch(facecolor=GROUP_COLOURS.get(g, "0.6"), edgecolor="white",
                         label=GROUP_LABELS.get(g, g)) for g in present]
        if drew_signal:
            handles.append(Line2D([], [], color="crimson", lw=1.8, label="signal"))
        if data:
            handles.append(Line2D([], [], color="black", marker="o", ls="", ms=4, label="data"))
        if legend_outside:
            # under the figure, so a tall signal overlay is never hidden by it
            axes[0].legend([], [], title=region_a_label, frameon=False,
                           title_fontsize=11, loc="upper right", alignment="right")
            fig.legend(handles=handles, frameon=False, fontsize=10,
                       ncol=len(handles), loc="upper center", bbox_to_anchor=(0.5, 0.0))
        else:
            axes[0].legend(handles=handles, title=region_a_label,
                           frameon=False, fontsize=9, title_fontsize=11,
                           loc="upper right", ncol=2, alignment="right")
        # with data drawn, not "Simulation"
        cms_label(axes, fontsize=13, **({"data": True, "text": "Work in progress"} if data else {}))
        if title:
            fig.suptitle(title, y=1.06)
        return save(fig, outdir, name or f"regions_{channel}", campaign)


def closure_panel(grouped, channel, outdir, campaign, subtitle="",
                  label=r"$m_{LJLJ}$ [GeV]", name=None, title="", channel_label=None,
                  per=None):
    """Region A against the per-bin B*C/D prediction, with a ratio panel.

    ``channel_label`` names the selection in the panel titles when the channel
    name alone would mislead (a validation region filled from the SR channel).

    Ratio conventions, which are the thing to get right: the **points** carry
    the region-A statistical error only, and the **band around unity** is the
    prediction's relative error.  Merging both into the point errors hides which
    side is imprecise.  The chi2 does combine them, which is what a
    compatibility test needs.
    """
    edges, regions, pred, pred_err = region_arrays(grouped, channel)
    centres = 0.5 * (edges[:-1] + edges[1:])
    a, a_err = regions[0]
    # drawn per ``per`` GeV when asked; the ratio and chi2 use the counts
    dens = per / np.diff(edges) if per else np.ones(len(edges) - 1)
    ylab = f"events / {per:g} GeV" if per else "events"

    with plt.rc_context(GRID_RC):
        # Left is one tall axis; right is a plot above its ratio, sharing x.
        # A plain 2x2 with sharex would hide the left panel's tick labels.
        fig = plt.figure(figsize=(12.4, 6.0), constrained_layout=True)
        gs = fig.add_gridspec(2, 2, height_ratios=[2.4, 1])
        ax, ax_top = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1])
        rax = fig.add_subplot(gs[1, 1], sharex=ax_top)

        for r in range(4):
            w, err = regions[r]
            _step(ax, edges, w * dens, color=REGION_COLOURS[r], lw=1.8,
                  label=f"region {REGION_NAMES[r]}")
            ax.errorbar(centres, w * dens, yerr=err * dens, fmt="none",
                        ecolor=REGION_COLOURS[r], elinewidth=1.0, alpha=0.7)
        ax.set_yscale("log")
        ax.set_ylabel(ylab)
        ax.set_xlabel(label)
        ax.set_title(f"{channel_label or channel}: the four regions")
        ax.legend(frameon=False, fontsize=9, ncol=2)
        ax.grid(alpha=0.25, which="both")
        _log_ylim(ax, [regions[r][0] * dens for r in range(4)])

        _step(ax_top, edges, a * dens, color="black", lw=2, label="region A (simulation)")
        ax_top.errorbar(centres, a * dens, yerr=a_err * dens, fmt="none", ecolor="black",
                        elinewidth=1.2, capsize=2)
        _step(ax_top, edges, pred * dens, color=THEORY_COLOUR, lw=2,
              label=r"$B\,C/D$ prediction")
        _band(ax_top, edges, (pred - pred_err) * dens, (pred + pred_err) * dens,
              color=THEORY_COLOUR, alpha=0.2)
        ax_top.set_yscale("log")
        ax_top.set_ylabel(ylab)
        ax_top.set_title("closure per bin")
        ax_top.legend(frameon=False, fontsize=9)
        ax_top.grid(alpha=0.25, which="both")
        _log_ylim(ax_top, [a * dens, (pred * dens)[np.isfinite(pred)]], lo=0.2)

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(pred > 0, a / pred, np.nan)
            ratio_err = np.where((pred > 0) & (a > 0), a_err / pred, np.nan)
            pred_rel = np.where(pred > 0, pred_err / pred, np.nan)
        _band(rax, edges, 1 - pred_rel, 1 + pred_rel, color=THEORY_COLOUR,
              alpha=0.2, label="prediction uncertainty")
        rax.errorbar(centres, ratio, yerr=ratio_err, fmt="o", color="black",
                     ms=5, capsize=3, label="region A stat.")
        rax.axhline(1.0, color=THEORY_COLOUR, lw=1.6)
        rax.set_yscale("log")
        rax.set_ylim(0.05, 20)
        rax.set_ylabel(r"A / $B\,C/D$")
        rax.set_xlabel(label)
        rax.grid(alpha=0.25, which="both")
        rax.legend(frameon=False, fontsize=7.5, loc="lower left", ncol=2)

        good = np.isfinite(ratio)
        if good.sum():
            chi2 = np.nansum(((a[good] - pred[good])
                              / np.sqrt(a_err[good] ** 2 + pred_err[good] ** 2)) ** 2)
            rax.text(0.99, 0.93, rf"$\chi^2/\mathrm{{ndf}} = {chi2:.1f}/{good.sum()}$",
                     transform=rax.transAxes, ha="right", va="top", fontsize=10)
        if subtitle:
            fig.text(0.5, -0.02, subtitle, ha="center", va="top", fontsize=9,
                     family="monospace")
        fig.suptitle(title or f"{channel_label or channel} -- per-bin closure, background MC",
                     y=1.03)
        cms_label([ax, ax_top], fontsize=11)
        return save(fig, outdir, name or f"closure_{channel}", campaign)


def shape_factorisation(grouped, channel, outdir, campaign,
                        label=r"$m_{LJLJ}$ [GeV]", name=None, title=""):
    """The four regions area-normalised, and the A/B against C/D shape ratios.

    Factorisation would make those two ratios the same function of the
    observable -- that is the assumption the shape fit rests on, drawn directly.
    """
    edges, regions, _, _ = region_arrays(grouped, channel)
    centres = 0.5 * (edges[:-1] + edges[1:])
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)
        ax = axes[0]
        for r in range(4):
            w, err = regions[r]
            total = w.sum()
            if total <= 0:
                continue
            _step(ax, edges, w / total, color=REGION_COLOURS[r], lw=1.8,
                  label=f"region {REGION_NAMES[r]}")
            ax.errorbar(centres, w / total, yerr=err / total, fmt="none",
                        ecolor=REGION_COLOURS[r], elinewidth=1.0, alpha=0.7)
        ax.set_yscale("log")
        ax.set_xlabel(label)
        ax.set_ylabel("fraction of region")
        ax.set_title("area-normalised")
        ax.legend(frameon=False, fontsize=9, ncol=2)
        ax.grid(alpha=0.25, which="both")

        ax = axes[1]
        wA, wB, wC, wD = (regions[r][0] for r in range(4))
        with np.errstate(divide="ignore", invalid="ignore"):
            ab = np.where(wB > 0, (wA / wA.sum()) / (wB / wB.sum()), np.nan)
            cd = np.where(wD > 0, (wC / wC.sum()) / (wD / wD.sum()), np.nan)
        _step(ax, edges, ab, color="#4c72b0", lw=2, label="A/B (normalised)")
        _step(ax, edges, cd, color="#dd8452", lw=2, label="C/D (normalised)")
        ax.axhline(1.0, color="0.4", lw=1.0, ls="--")
        ax.set_yscale("log")
        ax.set_xlabel(label)
        ax.set_ylabel("shape ratio")
        ax.set_title("factorisation would make these two agree")
        ax.legend(frameon=False, fontsize=9)
        ax.grid(alpha=0.25, which="both")
        fig.suptitle(title or f"{channel} -- does the shape factorise?", y=1.04)
        cms_label(axes, fontsize=12)
        return save(fig, outdir, name or f"factorisation_{channel}", campaign)


def fit_closure_ratio(edges, ratio, err):
    """Fit log(ratio) against log(x): is the bias a constant, or does it slope?

    Returns the best-fit multiplicative ``factor`` (the intercept, which is what
    a correction would use), the ``slope`` with its error, and the chi2 of both
    hypotheses.  A constant offset is ``slope`` consistent with zero.
    """
    centres = 0.5 * (edges[:-1] + edges[1:])
    good = (np.isfinite(ratio) & (ratio > 0) & np.isfinite(err) & (err > 0))
    if good.sum() < 3:
        return None
    x, y = np.log(centres[good]), np.log(ratio[good])
    w = 1.0 / (err[good] / ratio[good]) ** 2

    const = np.sum(w * y) / np.sum(w)
    chi2_const = np.sum(w * (y - const) ** 2)
    design = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(design.T @ (design * w[:, None]))
    coef = cov @ (design.T @ (w * y))
    return {"factor": float(np.exp(const)),
            "factor_err": float(np.exp(const) / np.sqrt(np.sum(w))),
            "chi2_const": float(chi2_const), "ndf_const": int(good.sum() - 1),
            "slope": float(coef[1]), "slope_err": float(np.sqrt(cov[1, 1])),
            "chi2_line": float(np.sum(w * (y - design @ coef) ** 2)),
            "ndf_line": int(good.sum() - 2), "n": int(good.sum())}


def closure_ratio_overlay(variants, channel, outdir, campaign,
                          label=r"$m_{LJLJ}$ [GeV]", name=None, title=""):
    """Closure ratios from several ``{label: grouped}`` variants on one axis.

    Used to ask whether a bias is a stable constant across working points -- the
    dashed line per variant is its best-fit constant.
    """
    rows = []
    with plt.rc_context(GRID_RC):
        fig, ax = plt.subplots(figsize=(7.0, 4.4), constrained_layout=True)
        for i, (vlabel, grouped) in enumerate(variants.items()):
            edges, regions, pred, pred_err = region_arrays(grouped, channel)
            a, a_err = regions[0]
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(pred > 0, a / pred, np.nan)
                err = np.abs(ratio) * np.sqrt(
                    np.where(a > 0, (a_err / a) ** 2, 0)
                    + np.where(pred > 0, (pred_err / pred) ** 2, 0))
            fit = fit_closure_ratio(edges, ratio, err)
            if fit is None:
                continue
            rows.append(dict(variant=vlabel, channel=channel, **fit))
            colour = style_for(i + 1)[0]
            centres = 0.5 * (edges[:-1] + edges[1:])
            ax.errorbar(centres, ratio, yerr=err, fmt="o", ms=5, color=colour,
                        capsize=3, alpha=0.85,
                        label=f"{vlabel}  (factor {fit['factor']:.2f})")
            ax.axhline(fit["factor"], color=colour, ls="--", lw=1.2, alpha=0.8)
        ax.axhline(1.0, color="0.3", lw=1.4)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(label)
        ax.set_ylabel(r"A / $(B\,C/D)$")
        ax.set_title(title or f"{channel}: is the bias a constant factor?")
        ax.grid(alpha=0.25, which="both")
        ax.legend(frameon=False, fontsize=8)
        cms_label(ax)
        save(fig, outdir, name or f"closure_ratio_{channel}", campaign)
    return rows


def compare_region_shapes(variants, channel, region, outdir, campaign,
                          reference=None, label=r"$m_{LJLJ}$ [GeV]",
                          name=None, title=""):
    """One region's shape under several selections, area-normalised with a ratio.

    Answers "does this cut distort the distribution?" -- if it does not, the
    higher-statistics variant can measure a shape the tighter one cannot.
    """
    reference = reference or next(iter(variants))
    stats = {}
    with plt.rc_context(GRID_RC):
        fig, (ax, rax) = plt.subplots(2, 1, figsize=(7.0, 5.6), sharex=True,
                                      constrained_layout=True,
                                      gridspec_kw={"height_ratios": [2.2, 1]})
        ref = None
        for i, (vlabel, grouped) in enumerate(variants.items()):
            import shape_tools as st
            shape = st.total_background_shape(grouped, channel, region)
            edges, w, err = shape.edges, shape.sumw, np.sqrt(shape.sumw2)
            total = w.sum()
            if total <= 0:
                continue
            neff = total ** 2 / shape.sumw2.sum()
            frac, frac_err = w / total, err / total
            stats[vlabel] = (total, neff)
            colour = style_for(i)[0]
            _step(ax, edges, frac, color=colour, lw=1.9,
                  label=f"{vlabel}  ({total:.0f} ev, {neff:.0f} eff)")
            centres = 0.5 * (edges[:-1] + edges[1:])
            ax.errorbar(centres, frac, yerr=frac_err, fmt="none", ecolor=colour,
                        elinewidth=1.0, alpha=0.7)
            if vlabel == reference:
                ref = (frac, frac_err)
            elif ref is not None:
                with np.errstate(divide="ignore", invalid="ignore"):
                    r = np.where(ref[0] > 0, frac / ref[0], np.nan)
                    re_ = np.abs(r) * np.sqrt(
                        np.where(frac > 0, (frac_err / frac) ** 2, 0)
                        + np.where(ref[0] > 0, (ref[1] / ref[0]) ** 2, 0))
                rax.errorbar(centres, r, yerr=re_, fmt="o", ms=4, color=colour,
                             capsize=2, alpha=0.85)
        ax.set_yscale("log")
        ax.set_ylabel(f"fraction of region {REGION_NAMES[region]}")
        ax.grid(alpha=0.25, which="both")
        ax.legend(frameon=False, fontsize=8)
        rax.axhline(1.0, color="0.3", lw=1.4)
        rax.set_yscale("log")
        rax.set_ylim(0.05, 20)
        rax.set_xlabel(label)
        rax.set_ylabel(f"/ {reference}")
        rax.grid(alpha=0.25, which="both")
        fig.suptitle(title or f"{channel}, region {REGION_NAMES[region]}: "
                              f"shape vs selection", y=1.03)
        cms_label(ax)
        save(fig, outdir, name or f"shape_variants_{channel}", campaign)
    return stats


# --------------------------------------------------------------------------- #
# Systematic uncertainties and fit diagnostics
# --------------------------------------------------------------------------- #
VARIATION_COLOURS = {"up": THEORY_COLOUR, "down": OVERLAY_COLOUR}


def variation_plot(nominal, varied, outdir, campaign, name, title=""):
    """A template and its systematic variations in every region, with their ratio to it.

    ``nominal`` is ``{region: Shape}`` and ``varied`` is ``{label: {region: Shape}}``,
    such as the up and down variations of one nuisance.  The ratio panels shade the
    statistical uncertainty of the nominal template, for scale: a variation that
    follows the bin-to-bin fluctuations inside it is more likely noise than a shape.
    """
    regions = sorted(nominal)
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(2, len(regions), figsize=(3.6 * len(regions), 4.8),
                                 sharex=True, squeeze=False, constrained_layout=True,
                                 gridspec_kw={"height_ratios": [3, 1.4]})
        for j, region in enumerate(regions):
            top, bottom = axes[0][j], axes[1][j]
            nom = nominal[region]
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = np.where(nom.sumw > 0, np.sqrt(nom.sumw2) / nom.sumw, 0.0)
                _band(bottom, nom.edges, 1 - rel, 1 + rel, color="0.85",
                      label="statistical uncertainty of the nominal")
                _step(top, nom.edges, nom.sumw, color="black", lw=1.6, label="nominal")
                for i, (label, shapes) in enumerate(varied.items()):
                    colour = VARIATION_COLOURS.get(label, style_for(i + 1)[0])
                    values = shapes[region].sumw
                    _step(top, nom.edges, values, color=colour, lw=1.2, ls="--", label=label)
                    _step(bottom, nom.edges, np.where(nom.sumw > 0, values / nom.sumw, np.nan),
                          color=colour, lw=1.2)
            bottom.axhline(1.0, color="0.3", lw=0.8)
            top.set_title(f"region {REGION_NAMES.get(region, region)}")
            bottom.set_xlabel(r"$m_{LJLJ}$ [GeV]")
            if j == 0:
                top.set_ylabel("events")
                bottom.set_ylabel("ratio to nominal")
        handles, labels = axes[0][0].get_legend_handles_labels()
        band_handles, band_labels = axes[1][0].get_legend_handles_labels()
        fig.legend(handles + band_handles, labels + band_labels, frameon=False,
                   ncol=len(labels) + len(band_labels), loc="upper center",
                   bbox_to_anchor=(0.5, 0.0))
        fig.suptitle(title)
        cms_label(axes, fontsize=11)
        return save(fig, outdir, name, campaign)


def impact_plot(impacts, outdir, campaign, name, title=""):
    """Pulls and impacts of the nuisances of one fit, from combineTool.py -M Impacts.

    Left: each nuisance's best fit relative to its pre-fit value, with its post-fit
    uncertainty, both in units of the pre-fit uncertainty.  Right: the shift of r when
    the nuisance is fixed one post-fit standard deviation above (red) or below (blue)
    its best fit, relative to the best-fit r.  The largest impact is at the top.
    """
    r_lo, r_best, r_hi = impacts["POIs"][0]["fit"]

    def shifts(par):
        return [100 * (par["r"][k] - par["r"][1]) / r_best for k in (2, 0)]

    params = sorted(impacts["params"], key=lambda par: max(map(abs, shifts(par))))
    with plt.rc_context(GRID_RC):
        fig, (left, right) = plt.subplots(1, 2, sharey=True,
                                          figsize=(11, 1.8 + 0.42 * len(params)),
                                          gridspec_kw={"width_ratios": [1, 1.3]},
                                          constrained_layout=True)
        for width, shade in ((2, "0.93"), (1, "0.82")):
            left.axvspan(-width, width, color=shade, lw=0)
        for i, par in enumerate(params):
            lo, best, hi = par["fit"]
            pre_lo, pre, pre_hi = par["prefit"]
            sigma = (pre_hi - pre_lo) / 2 or 1.0
            left.errorbar([(best - pre) / sigma], [i],
                          xerr=[[(best - lo) / sigma], [(hi - best) / sigma]],
                          fmt="o", color="black", ms=4)
            up, down = shifts(par)
            right.barh(i, up, height=0.6, color=VARIATION_COLOURS["up"])
            right.barh(i, down, height=0.6, color=VARIATION_COLOURS["down"])
        left.set_xlim(-2.5, 2.5)
        left.set_yticks(range(len(params)))
        left.set_yticklabels([par["name"] for par in params])
        left.set_xlabel(r"$(\hat\theta-\theta_0)/\Delta\theta$")
        right.axvline(0, color="0.3", lw=0.8)
        largest = max([abs(v) for par in params for v in shifts(par)] + [1e-3])
        right.set_xlim(-1.15 * largest, 1.15 * largest)
        right.set_xlabel(r"$\Delta r/\hat r$ [%]")
        right.legend(handles=[matplotlib.patches.Patch(color=c, label=label)
                              for label, c in ((r"$+1\sigma$", VARIATION_COLOURS["up"]),
                                               (r"$-1\sigma$", VARIATION_COLOURS["down"]))],
                     frameon=False, loc="lower right")
        fig.suptitle(f"{title}\n" + rf"$\hat r = {r_best:.3g}"
                     rf"^{{+{r_hi - r_best:.2g}}}_{{-{r_best - r_lo:.2g}}}$")
        cms_label([left, right], fontsize=11)
        return save(fig, outdir, name, campaign)


def scan_plot(panels, outdir, campaign, name, title="", styles=None, ncols=3):
    """Likelihood scans, one panel per scanned parameter.

    ``panels`` is a list of ``(panel title, {label: scan}, prior)``, a scan being the
    output of run_combine_scan.py.  ``prior`` adds theta^2, what the unit Gaussian
    constraint of a nuisance contributes by itself: a scan on top of it means the
    data add nothing to the constraint.  The dotted lines are at 1 and 3.84, where a
    scan of one parameter crosses its 68% and 95% CL interval.
    """
    n = len(panels)
    ncols = min(ncols, n)
    nrows = -(-n // ncols)
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.7 * nrows),
                                 squeeze=False, constrained_layout=True)
        for k in range(n, nrows * ncols):
            axes[k // ncols][k % ncols].set_visible(False)
        for k, (panel_title, scans, prior) in enumerate(panels):
            ax = axes[k // ncols][k % ncols]
            xs = []
            for i, (label, scan) in enumerate(scans.items()):
                kw = {**_series_kwargs(i, label, styles, lw=1.8), "marker": None}
                ax.plot(scan["x"], scan["dnll2"], label=label, **kw)
                xs += scan["x"]
            if prior and xs:
                x = np.linspace(min(xs), max(xs), 200)
                ax.plot(x, x ** 2, color="0.45", ls=":", lw=1.6,
                        label=r"constraint alone, $\theta^2$")
            for level in (1.0, 3.84):
                ax.axhline(level, color="0.6", lw=0.8, ls=(0, (1, 2)))
            ax.set_ylim(0, 10)
            poi = next(iter(scans.values()))["poi"] if scans else ""
            ax.set_xlabel(r"$r$ [fb]" if poi == "r" else r"$\theta$")
            ax.set_ylabel(r"$-2\Delta\ln L$")
            ax.set_title(panel_title)
            ax.legend(frameon=False, fontsize=8)
        fig.suptitle(title)
        cms_label(axes, fontsize=11)
        return save(fig, outdir, name, campaign)


def correlation_plot(names, matrix, outdir, campaign, name, title="", annotate=0.1):
    """Correlation matrix of the parameters of a fit, in the order given; coefficients
    at least ``annotate`` in size are printed."""
    m = np.asarray(matrix, dtype=float)
    n = len(names)
    size = max(5.0, 0.3 * n + 2.0)
    with plt.rc_context(GRID_RC):
        fig, ax = plt.subplots(figsize=(size + 1.5, size), constrained_layout=True)
        im = ax.imshow(m, cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(n))
        ax.set_xticklabels(names, rotation=90, fontsize=7)
        ax.set_yticks(range(n))
        ax.set_yticklabels(names, fontsize=7)
        for i in range(n):
            for j in range(n):
                if i != j and abs(m[i, j]) >= annotate:
                    ax.text(j, i, f"{m[i, j]:.2f}", ha="center", va="center", fontsize=5,
                            color="white" if abs(m[i, j]) > 0.6 else "black")
        fig.colorbar(im, ax=ax, label="correlation", shrink=0.8)
        ax.set_title(title)
        cms_label(ax, fontsize=11)
        return save(fig, outdir, name, campaign)


def bias_plot(toys, outdir, campaign, name, title=""):
    """Signal-injection test of one card, from run_combine_toys.py.

    Left: the fitted r of the toys, against the injected one, with its mean.  Middle:
    its pull, (r - r_injected) over the uncertainty on the side of the injected value,
    against a unit Gaussian (pulls beyond 4 in size are drawn in the outer bins), with
    its median and the half-width of its central 68%, which are the mean and width for
    a Gaussian but are not thrown off by a few toys with a broken uncertainty; toys
    whose best fit would lie below the r = 0 boundary stop there and gather in the left
    tail.  Right: the mean and spread over the toys of every fitted nuisance, against
    the 0 +- 1 that a nuisance the data do not constrain follows.  Only converged fits
    (status 0) are drawn.
    """
    injected = toys["injected"]
    r, lo, hi, status = (np.asarray(toys[k], dtype=float) for k in ("r", "r_lo", "r_hi", "status"))
    ok = status == 0
    with np.errstate(divide="ignore", invalid="ignore"):
        pull = np.where(r < injected, (r - injected) / hi, (r - injected) / lo)
    pull = pull[ok & np.isfinite(pull)]
    at_zero = int(np.sum(r[ok] <= 1e-6 * max(injected, 1.0)))
    names = sorted(toys["nuisances"])
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, 3, figsize=(15.5, max(4.2, 1.5 + 0.35 * len(names))),
                                 constrained_layout=True, gridspec_kw={"width_ratios": [1, 1, 1.15]})
        fitted, pulls, nuisances = axes
        fitted.hist(r[ok], bins=30, color="0.65",
                    label=f"mean {r[ok].mean():.3g} $\\pm$ {r[ok].std() / np.sqrt(max(ok.sum(), 1)):.2g}")
        fitted.axvline(injected, color=THEORY_COLOUR, lw=1.8, label=f"injected r = {injected:.3g}")
        fitted.set_xlabel(r"fitted $r$ [fb]")
        fitted.set_ylabel("toys")
        fitted.legend(frameon=False, fontsize=9)
        bins = np.linspace(-4, 4, 33)
        if pull.size:
            q16, q50, q84 = np.percentile(pull, [16, 50, 84])
            pulls.hist(np.clip(pull, -3.99, 3.99), bins=bins, color="0.65",
                       label=f"median {q50:.2f}\n68% half-width {(q84 - q16) / 2:.2f}\n"
                             f"beyond $\\pm$4: {int(np.sum(np.abs(pull) > 4))}")
        x = np.linspace(-4, 4, 200)
        pulls.plot(x, pull.size * (bins[1] - bins[0]) * np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi),
                   color="black", lw=1.4, label="unit Gaussian")
        pulls.set_xlabel(r"$(\hat r - r_\mathrm{inj})/\sigma_r$")
        pulls.legend(frameon=False, fontsize=9)
        nuisances.axvspan(-1, 1, color="0.88", lw=0)
        for i, nuisance in enumerate(names):
            values = np.asarray(toys["nuisances"][nuisance], dtype=float)[ok]
            nuisances.errorbar([values.mean()], [i], xerr=[values.std()], fmt="o",
                               color="black", ms=4)
        nuisances.set_yticks(range(len(names)))
        nuisances.set_yticklabels(names, fontsize=9)
        nuisances.set_xlim(-2.5, 2.5)
        nuisances.set_xlabel(r"fitted $\theta$: mean $\pm$ spread over the toys")
        fig.suptitle(f"{title}: {int(ok.sum())} of {len(r)} fits converged, {at_zero} at r = 0")
        cms_label(axes, fontsize=11)
        return save(fig, outdir, name, campaign)


def breakdown_plot(changes, outdir, campaign, name, title="", ylabel="change of the limit [%]"):
    """Per-point changes of the expected limit, one column per label and one panel per
    final state: every signal point is a dot, spread sideways, and the median a bar.

    ``changes`` is ``{final state: {label: [change at each signal point]}}``.
    """
    states = list(changes)
    columns = max(len(per_label) for per_label in changes.values())
    jitter = np.random.default_rng(0)
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, len(states), figsize=(1.1 * columns * len(states) + 2, 4.8),
                                 sharey=True, squeeze=False, constrained_layout=True)
        for ax, fs in zip(axes[0], states):
            for i, (label, values) in enumerate(changes[fs].items()):
                v = np.asarray(values, dtype=float)
                ax.scatter(i + jitter.uniform(-0.25, 0.25, v.size), v, s=10,
                           color=style_for(i)[0], alpha=0.6, lw=0, label=label)
                if v.size:
                    ax.plot([i - 0.32, i + 0.32], [np.median(v)] * 2, color="black", lw=2)
            ax.set_xticks(range(len(changes[fs])))
            ax.set_xticklabels(list(changes[fs]), rotation=30, ha="right")
            ax.axhline(0, color="0.5", lw=0.8)
            ax.set_title(fs)
            ax.grid(alpha=0.25, axis="y")
        axes[0][0].set_ylabel(ylabel)
        fig.suptitle(title)
        cms_label(axes, fontsize=11)
        return save(fig, outdir, name, campaign)
