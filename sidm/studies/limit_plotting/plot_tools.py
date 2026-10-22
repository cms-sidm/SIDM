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
GROUP_COLOURS = {"QCD": "#4c72b0", "DY": "#dd8452", "TT": "#55a868",
                 "Diboson": "#c44e52", "other": "#8172b3"}
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


# Run 2 (2018) conditions, for the mplhep CMS label on every figure.
LUMI_FB = 59.83
COM_TEV = 13


def cms_label(ax, lumi=LUMI_FB, com=COM_TEV, title_pad=22, **kwargs):
    """The standard CMS Simulation label with the Run 2 luminosity.

    Everything here is simulation, so ``data=False`` -- which is what makes
    mplhep write "Simulation" under "CMS" rather than "Preliminary".  On a
    multi-panel figure the label belongs on the top-left axes only, as CMS
    style requires; an array of axes is accepted and reduced to its first
    element, since ``subplots(squeeze=False)`` hands back a 2-D array and
    mplhep needs a single Axes.
    """
    group = np.atleast_1d(np.asarray(ax, dtype=object)).ravel()
    # mplhep writes the label just above the axes box, which is exactly where a
    # panel title already sits.  Push the titles up so the two do not overprint
    # -- and push *every* panel's title, not only the labelled one, so a row of
    # column headings stays on one line.
    for a in group:
        if a.get_title():
            a.set_title(a.get_title(), pad=title_pad,
                        fontsize=a.title.get_fontsize())
    return hep.cms.label(ax=group[0], data=False, lumi=lumi, com=com, **kwargs)


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
    """
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if not row.get("exp"):
                continue
            entry = {"exp": float(row["exp"])}
            if with_bands:
                entry.update({k: float(row[k]) for k in QUANTILES if row.get(k)})
            if row.get("obs"):
                entry["obs"] = float(row["obs"])
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
def draw_brazil(ax, rows, theory=None, label_median="median expected"):
    """Brazil bands against Lxy for rows that share a panel."""
    rows = sorted([r for r in rows if r["lxy"] == r["lxy"]], key=lambda r: r["lxy"])
    if not rows:
        return False
    x = np.array([r["lxy"] for r in rows])
    if all("exp_m2" in r and "exp_p2" in r for r in rows):
        ax.fill_between(x, [r["exp_m2"] for r in rows], [r["exp_p2"] for r in rows],
                        color=BAND_STYLE["2sigma"], label=r"expected $\pm2\sigma$")
        ax.fill_between(x, [r["exp_m1"] for r in rows], [r["exp_p1"] for r in rows],
                        color=BAND_STYLE["1sigma"], label=r"expected $\pm1\sigma$")
    ax.plot(x, [r["exp"] for r in rows], color="black", ls="--", lw=1.8,
            marker="o", ms=4, label=label_median)
    if theory is not None and theory == theory:
        ax.axhline(theory, color=THEORY_COLOUR, lw=1.8,
                   label=r"$\sigma_\mathrm{theory}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    return True


def brazil_grid(rows, fs, outdir, campaign, theory=None, name=None, title=""):
    """Brazil bands across the grid, one panel per (m_bound, m_ZD)."""
    theory = theory if theory is not None else theory_by_mass()
    sub = [r for r in rows if r["fs"] == fs]
    masses = sorted({r["m_bound"] for r in sub})
    mzds = sorted({r["mzd"] for r in sub})
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(len(mzds), len(masses), sharex=True, sharey=True,
                                 figsize=(3.1 * len(masses), 2.5 * len(mzds)),
                                 squeeze=False, constrained_layout=True)
        for i, mzd in enumerate(mzds):
            for j, mb in enumerate(masses):
                ax = axes[i][j]
                panel = [r for r in sub if r["mzd"] == mzd and r["m_bound"] == mb]
                if not draw_brazil(ax, panel, theory.get(mb, np.nan)):
                    ax.text(0.5, 0.5, "no points", transform=ax.transAxes,
                            ha="center", va="center", color="0.6")
                if i == 0:
                    ax.set_title(rf"$m_B={mb:g}$ GeV")
                if j == 0:
                    ax.set_ylabel(rf"$m_{{Z_D}}={mzd:g}$ GeV" + "\n"
                                  + r"$\sigma$ limit [fb]")
                if i == len(mzds) - 1:
                    ax.set_xlabel(r"$L_{xy}$ [cm]")
                ax.grid(alpha=0.25, which="both")
        handles, labels = axes[0][-1].get_legend_handles_labels()
        fig.legend(handles, labels, frameon=False, ncol=4, loc="lower center",
                   bbox_to_anchor=(0.5, -0.04))
        fig.suptitle(title or f"{fs}: expected 95% CL limit (MC only, blinded)")
        cms_label(axes, fontsize=11)
        return save(fig, outdir, name or f"brazil_grid_{fs}", campaign)


def limit_vs_mass(rows, fs, outdir, campaign, theory=None, name=None, title=""):
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
            ax.plot([], [], color="0.35", ls=style, label=rf"$m_{{Z_D}}={mzd:g}$ GeV")
        ax.set_yscale("log")
        ax.set_xlabel(r"bound state mass $m_B$ [GeV]")
        ax.set_ylabel(r"median expected 95% CL limit on $\sigma$ [fb]")
        ax.set_title(title or fs)
        ax.grid(alpha=0.25, which="both")
        ax.legend(frameon=False, ncol=2, fontsize=9)
        cms_label(ax)
        return save(fig, outdir, name or f"limit_vs_mass_{fs}", campaign)


def exclusion_map(rows, fs, outdir, campaign, name=None, title=""):
    """sigma_limit / sigma_theory across the grid; cells below 1 are outlined."""
    sub = [r for r in rows if r["fs"] == fs and r["r_theory"] == r["r_theory"]]
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
                    hit = [r["r_theory"] for r in sub if r["mzd"] == mzd
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
                    if v < 1:
                        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                                   fill=False, edgecolor="red", lw=2))
            ax.set_xticks(range(len(masses)))
            ax.set_xticklabels([f"{m:g}" for m in masses])
            ax.set_yticks(range(len(lxys)))
            ax.set_yticklabels([f"{l:g}" for l in lxys])
            ax.set_xlabel(r"$m_B$ [GeV]")
            if k == 0:
                ax.set_ylabel(r"$L_{xy}$ [cm]")
            ax.set_title(rf"$m_{{Z_D}}={mzd:g}$ GeV")
        fig.colorbar(im, ax=axes[0].tolist(),
                     label=r"$\sigma_\mathrm{limit}/\sigma_\mathrm{theory}$")
        n_excl = int(np.nansum(grid < 1))
        fig.suptitle(title or f"{fs}: expected exclusion "
                              f"(outlined: below 1, {n_excl} cells)")
        cms_label(axes, fontsize=11)
        return save(fig, outdir, name or f"exclusion_map_{fs}", campaign)


def kinetic_mixing(rows, outdir, campaign, name="eps2_vs_mdp", title=""):
    """Expected limit in the eps^2 vs dark photon mass plane, one panel per mass."""
    sub = [r for r in rows if r["eps2"] == r["eps2"]]
    masses = sorted({r["m_bound"] for r in sub})
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, len(masses), figsize=(2.7 * len(masses), 3.6),
                                 sharey=True, squeeze=False, constrained_layout=True)
        vals = [r["exp"] for r in sub]
        norm = matplotlib.colors.LogNorm(vmin=min(vals), vmax=max(vals))
        for j, mb in enumerate(masses):
            ax = axes[0][j]
            panel = [r for r in sub if r["m_bound"] == mb]
            sc = ax.scatter([r["mzd"] for r in panel], [r["eps2"] for r in panel],
                            c=[r["exp"] for r in panel], norm=norm, cmap="viridis_r",
                            s=42, edgecolor="0.3", linewidth=0.4)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(r"$m_{Z_D}$ [GeV]")
            if j == 0:
                ax.set_ylabel(r"$\epsilon^2$")
            ax.set_title(rf"$m_B={mb:g}$ GeV")
            ax.grid(alpha=0.25, which="both")
        fig.colorbar(sc, ax=axes[0].tolist(), label=r"expected limit on $\sigma$ [fb]")
        fig.suptitle(title or r"Kinetic mixing plane "
                              r"($\epsilon=\sqrt{80/m_{Z_D}c\tau}\times10^{-6}$)")
        cms_label(axes, fontsize=11)
        return save(fig, outdir, name, campaign)


def result_set(rows, outdir, campaign, theory=None, prefix=""):
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
    made.append(kinetic_mixing(rows, outdir, campaign, name=f"{prefix}eps2_vs_mdp"))
    return [m for m in made if m]


# --------------------------------------------------------------------------- #
# Comparisons -- any number of sets, whatever they are
# --------------------------------------------------------------------------- #
def compare_vs_lxy(sets, fs, outdir, campaign, theory=None, name=None, title="",
                   bands=False):
    """Best limit at each lifetime for every set, one panel per bound state mass."""
    theory = theory if theory is not None else theory_by_mass()
    masses = sorted({r["m_bound"] for rows in sets.values() for r in rows
                     if r["fs"] == fs})
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, len(masses), figsize=(2.9 * len(masses), 4.0),
                                 sharey=True, squeeze=False, constrained_layout=True)
        for j, mb in enumerate(masses):
            ax = axes[0][j]
            for i, (label, rows) in enumerate(sets.items()):
                colour, marker, ls = style_for(i)
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
                ax.plot(xs_, [best[x]["exp"] for x in xs_], color=colour,
                        marker=marker, ls=ls, ms=5, lw=1.7, label=label)
                if bands and all("exp_m1" in best[x] for x in xs_):
                    ax.fill_between(xs_, [best[x]["exp_m1"] for x in xs_],
                                    [best[x]["exp_p1"] for x in xs_],
                                    color=colour, alpha=0.18)
            if theory.get(mb):
                ax.axhline(theory[mb], color="0.3", ls=(0, (1, 1)), lw=1.6)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(r"$L_{xy}$ [cm]")
            ax.set_title(rf"$m_B={mb:g}$ GeV")
            ax.grid(alpha=0.25, which="both")
            if j == 0:
                ax.set_ylabel(r"best expected $\sigma$ limit [fb]")
        handles, labels = axes[0][0].get_legend_handles_labels()
        fig.legend(handles, labels, frameon=False, ncol=3, loc="lower center",
                   bbox_to_anchor=(0.5, -0.13), fontsize=9)
        fig.suptitle(title or f"{fs}: dotted is the model cross section")
        cms_label(axes, fontsize=11)
        return save(fig, outdir, name or f"compare_vs_lxy_{fs}", campaign)


def compare_vs_mass(sets, outdir, campaign, theory=None, name="compare_vs_mass",
                    title=""):
    """Best limit per bound state mass for every set, one panel per final state."""
    theory = theory if theory is not None else theory_by_mass()
    states = sorted({r["fs"] for rows in sets.values() for r in rows})
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, len(states), figsize=(5.5 * len(states), 4.4),
                                 squeeze=False, constrained_layout=True)
        for ax, fs in zip(axes[0], states):
            for i, (label, rows) in enumerate(sets.items()):
                colour, marker, ls = style_for(i)
                best = best_by(rows, "m_bound", fs)
                if not best:
                    continue
                ms = sorted(best)
                ax.plot(ms, [best[m] for m in ms], color=colour, marker=marker,
                        ls=ls, ms=6, lw=1.8, label=label)
            if theory:
                ax.plot(list(theory), list(theory.values()), color=THEORY_COLOUR,
                        lw=2.0, ls="--", marker="^", ms=5,
                        label=r"$\sigma_\mathrm{theory}$")
            ax.set_yscale("log")
            ax.set_xlabel(r"$m_B$ [GeV]")
            ax.set_ylabel(r"best expected $\sigma$ limit [fb]")
            ax.set_title(fs)
            ax.grid(alpha=0.25, which="both")
        axes[0][0].legend(frameon=False, fontsize=9)
        fig.suptitle(title or "Best expected limit per bound state mass")
        cms_label(axes, fontsize=12)
        return save(fig, outdir, name, campaign)


def compare_ratio(sets, reference, outdir, campaign, against="lxy",
                  name=None, title=""):
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
                colour, marker, ls = style_for(i)
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
                ax.plot(xs_, [float(np.median(by_x[x])) for x in xs_], color=colour,
                        marker=marker, ls=ls, ms=5, lw=1.7,
                        label=f"{label}  (median {np.median(allv):.2f})")
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
                      per_mass=False, title=""):
    """Points reaching expected exclusion, per set (optionally split by mass)."""
    with plt.rc_context(GRID_RC):
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


def region_distributions(grouped, channel, outdir, campaign, signal=None,
                         label=r"$m_{LJLJ}$ [GeV]", name=None, title=""):
    """The observable in all four ABCD regions, stacked by process group.

    Built from ``hist.Hist`` objects and drawn with ``mplhep.histplot``, so the
    stacking, the error band and the CMS label all follow the standard
    conventions.  ``signal`` is an optional ``{region: Shape}`` to overlay.

    Log scale, so a bin with no simulated background simply draws nothing --
    which is the thing worth seeing, and is left to the eye rather than
    annotated on the figure.
    """
    with plt.rc_context(GRID_RC):
        fig, axes = plt.subplots(1, 4, figsize=(16, 4.0), sharey=True,
                                 constrained_layout=True)
        for ax, region in zip(axes, range(4)):
            stack, labels, colours, total = [], [], [], None
            for group in sorted(grouped):
                shape = grouped[group].get(channel, {}).get(region)
                if shape is None or shape.sumw.sum() <= 0:
                    continue
                h = to_hist(shape.edges, shape.sumw, shape.sumw2)
                stack.append(h)
                labels.append(group)
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
            if sig is not None and sig.sumw.sum() > 0:
                hep.histplot(to_hist(sig.edges, sig.sumw), ax=ax, yerr=False,
                             histtype="step", color="crimson", linewidth=1.8,
                             label="signal")
                top = max(top, float(sig.sumw.max()))

            ax.set_yscale("log")
            ax.set_ylim(1e-3, max(top, 1e-2) * 10 ** 2.0)
            ax.set_xlabel(label)
            ax.set_ylabel("events" if ax is axes[0] else "")
            name_r = (f"region {REGION_NAMES[region]}"
                      + ("  (signal region)" if region == 0 else ""))
            handles, labels_ = (ax.get_legend_handles_labels() if region == 0
                                else ([], []))
            ax.legend(handles, labels_, title=name_r, frameon=False,
                      fontsize=9, title_fontsize=11, loc="upper right",
                      ncol=2, alignment="right")
        cms_label(axes, fontsize=13)
        if title:
            fig.suptitle(title, y=1.06)
        return save(fig, outdir, name or f"regions_{channel}", campaign)


def closure_panel(grouped, channel, outdir, campaign, subtitle="",
                  label=r"$m_{LJLJ}$ [GeV]", name=None, title=""):
    """Region A against the per-bin B*C/D prediction, with a ratio panel.

    Ratio conventions, which are the thing to get right: the **points** carry
    the region-A statistical error only, and the **band around unity** is the
    prediction's relative error.  Merging both into the point errors hides which
    side is imprecise.  The chi2 does combine them, which is what a
    compatibility test needs.
    """
    edges, regions, pred, pred_err = region_arrays(grouped, channel)
    centres = 0.5 * (edges[:-1] + edges[1:])
    a, a_err = regions[0]

    with plt.rc_context(GRID_RC):
        # Left is one tall axis; right is a plot above its ratio, sharing x.
        # A plain 2x2 with sharex would hide the left panel's tick labels.
        fig = plt.figure(figsize=(12.4, 6.0), constrained_layout=True)
        gs = fig.add_gridspec(2, 2, height_ratios=[2.4, 1])
        ax, ax_top = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1])
        rax = fig.add_subplot(gs[1, 1], sharex=ax_top)

        for r in range(4):
            w, err = regions[r]
            _step(ax, edges, w, color=REGION_COLOURS[r], lw=1.8,
                  label=f"region {REGION_NAMES[r]}")
            ax.errorbar(centres, w, yerr=err, fmt="none", ecolor=REGION_COLOURS[r],
                        elinewidth=1.0, alpha=0.7)
        ax.set_yscale("log")
        ax.set_ylabel("events")
        ax.set_xlabel(label)
        ax.set_title(f"{channel}: the four regions")
        ax.legend(frameon=False, fontsize=9, ncol=2)
        ax.grid(alpha=0.25, which="both")
        _log_ylim(ax, [regions[r][0] for r in range(4)])

        _step(ax_top, edges, a, color="black", lw=2, label="region A (simulation)")
        ax_top.errorbar(centres, a, yerr=a_err, fmt="none", ecolor="black",
                        elinewidth=1.2, capsize=2)
        _step(ax_top, edges, pred, color=THEORY_COLOUR, lw=2,
              label=r"$B\,C/D$ prediction")
        _band(ax_top, edges, pred - pred_err, pred + pred_err,
              color=THEORY_COLOUR, alpha=0.2)
        ax_top.set_yscale("log")
        ax_top.set_ylabel("events")
        ax_top.set_title("closure per bin")
        ax_top.legend(frameon=False, fontsize=9)
        ax_top.grid(alpha=0.25, which="both")
        _log_ylim(ax_top, [a, pred[np.isfinite(pred)]], lo=0.2)

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
        fig.suptitle(title or f"{channel} -- per-bin closure, background MC", y=1.03)
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
