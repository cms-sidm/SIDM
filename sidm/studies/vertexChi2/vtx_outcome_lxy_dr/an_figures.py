"""AN figures: vertex outcome of gen-matched muon-type lepton jets vs Lxy and vs dR.

Reads outcome_vs_lxy.csv, outcome_vs_dR.csv and outcome_vs_gen_dR.csv next to this script,
reference channel base_ljDisplacementIso_cosAlpha only, and only the unweighted mean over the
signal samples (<outcome>_mean) and its standard error (<outcome>_sem). Writes
vtx_outcome_vs_lxy.pdf, vtx_outcome_vs_dR.pdf (reconstructed dR of the two leading muons) and
vtx_outcome_vs_gen_dR.pdf (generated dR of the matched dark photon's two daughter muons), vector,
plus .png previews at 100 dpi into OUTDIR.

Left panel: the three outcome fractions stacked per bin, all samples. Right panel: the fraction
failing the vertex requirement, 1 - pass, per bound-state mass and for all samples.

    python sidm/studies/vertexChi2/vtx_outcome_lxy_dr/an_figures.py OUTDIR

Needs matplotlib, mplhep and numpy.
"""
import csv
import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3]))   # repository root
from sidm.studies.truth_kinematics_forAN import an_style  # pylint: disable=wrong-import-position

SELECTION = "base_ljDisplacementIso_cosAlpha"
OUTCOMES = ("no_vertex", "chi2_ge5", "pass")
CHI2_GE5 = r"$\chi^2/\mathrm{ndf} \geq 5$"
CHI2_LT5 = r"$\chi^2/\mathrm{ndf} < 5$"

# (table, bin column, x-axis title, output name, right-panel legend location)
FIGURES = [
    ("outcome_vs_lxy.csv", "genA_lxy_cm", r"Generated $Z_d$ $L_{xy}$ [cm]",
     "vtx_outcome_vs_lxy", "upper left"),
    ("outcome_vs_dR.csv", "mumu_dR",
     r"Reconstructed $\Delta R$ of the two leading muons",
     "vtx_outcome_vs_dR", "upper left"),
    ("outcome_vs_gen_dR.csv", "gen_mumu_dR",
     r"Generated $\Delta R$ of the $Z_d$ daughter muons",
     "vtx_outcome_vs_gen_dR", "upper right"),
]

STACK = [   # bottom to top
    ("no_vertex", "No stored vertex", "#964a8b"),
    ("chi2_ge5", f"Vertex, {CHI2_GE5}", "#f89c20"),
    ("pass", f"Vertex, {CHI2_LT5}", "#5790fc"),
]
MARKERS = ["s", "^", "v", "D", "P", "X"]
LEGEND_PT = 20


def read_table(path, bin_col):
    """Reference-channel rows as {group: {bin: {column: value}}}, and the bins in file order.

    Rows whose n_samples_in_mean is empty or zero hold no mean and are dropped. A bin with a
    single sample in the mean has an empty sem, read as NaN (no error bar)."""
    groups, bins = {}, []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            n = row["n_samples_in_mean"].strip()
            if row["selection"] != SELECTION or not n or float(n) == 0:
                continue
            b = row[bin_col]
            if b not in bins:
                bins.append(b)
            groups.setdefault(row["group"], {})[b] = {
                f"{o}_{s}": float(row[f"{o}_{s}"].strip() or "nan")
                for o in OUTCOMES for s in ("mean", "sem")}
    return groups, bins


def column(rows, bins, key):
    return np.array([rows[b][key] if b in rows else np.nan for b in bins])


def stack_panel(ax, rows, bins):
    x = np.arange(len(bins))
    bottom = np.zeros(len(bins))
    for outcome, label, color in STACK:
        h = column(rows, bins, f"{outcome}_mean")
        ax.bar(x, h, bottom=bottom, width=0.8, color=color, edgecolor="white", linewidth=0.8,
               label=label)
        bottom = bottom + h
    # Internal boundaries: no_vertex, and no_vertex + chi2_ge5 = 1 - pass, whose error is pass_sem.
    nv = column(rows, bins, "no_vertex_mean")
    ge5 = column(rows, bins, "chi2_ge5_mean")
    for y, err in ((nv, column(rows, bins, "no_vertex_sem")),
                   (nv + ge5, column(rows, bins, "pass_sem"))):
        ax.errorbar(x, y, yerr=err, fmt="none", ecolor="black", elinewidth=1.8, capsize=5,
                    capthick=1.8)
    ax.set_ylim(0, 1.25)
    ax.set_ylabel("Fraction of muon-type LJs")
    handles, labels = ax.get_legend_handles_labels()
    # Stack order top to bottom, filled column-wise: the two "Vertex" entries share a column.
    return ax.legend(handles[::-1], labels[::-1], ncol=2, loc="upper center",
                     fontsize=LEGEND_PT, handlelength=1.2, columnspacing=1.0,
                     borderaxespad=0.3)


def fail_panel(ax, groups, bins, legend_loc):
    masses = sorted((g for g in groups if g != "all"), key=lambda g: float(g.split()[1]))
    assert len(masses) <= len(MARKERS), f"{len(masses)} mass groups, {len(MARKERS)} markers"
    colors = plt.cm.plasma(np.linspace(0.08, 0.82, len(masses)))
    series = [(g, g.split(maxsplit=1)[1], c, m, 1.8, 9)
              for g, c, m in zip(masses, colors, MARKERS)]
    series.append(("all", "all samples", "black", "o", 3.2, 11))
    x = np.arange(len(bins))
    offsets = np.linspace(-0.21, 0.21, len(series))   # small dodge so error bars stay apart
    for (group, label, color, marker, lw, ms), dx in zip(series, offsets):
        rows = groups[group]
        ax.errorbar(x + dx, 1 - column(rows, bins, "pass_mean"),
                    yerr=column(rows, bins, "pass_sem"), color=color, marker=marker,
                    markersize=ms, linewidth=lw, elinewidth=lw * 0.8, capsize=0, label=label,
                    zorder=3 if group == "all" else 2)
    ax.set_ylim(0, 1)
    # Two lines: on one line the title is longer than the axis.
    ax.set_ylabel(f"Fraction without a vertex\nof {CHI2_LT5}")
    return ax.legend(ncol=2, loc=legend_loc, title="Bound-state mass", fontsize=LEGEND_PT,
                     title_fontsize=LEGEND_PT, columnspacing=1.0, handlelength=1.6,
                     borderaxespad=0.4)


def make_figure(table, bin_col, xtitle, legend_loc):
    groups, bins = read_table(HERE / table, bin_col)
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), layout="constrained")
    legends = [stack_panel(axes[0], groups["all"], bins),
               fail_panel(axes[1], groups, bins, legend_loc)]
    for ax in axes:
        ax.set_xticks(np.arange(len(bins)), bins, rotation=35, ha="right",
                      rotation_mode="anchor")
        ax.tick_params(axis="x", which="minor", bottom=False, top=False)
        ax.set_xlim(-0.6, len(bins) - 0.4)
        ax.set_xlabel(xtitle)
    an_style.cms_sim_labels(axes)
    return fig, axes, legends, groups, bins


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    an_style.set_style()
    an_style.FIG_DIR = outdir   # an_style.save writes <FIG_DIR>/<name>.pdf
    for table, bin_col, xtitle, name, legend_loc in FIGURES:
        fig = make_figure(table, bin_col, xtitle, legend_loc)[0]
        an_style.save(fig, name)
        fig.savefig(os.path.join(outdir, f"{name}.png"), dpi=100)
        plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"usage: python {sys.argv[0]} OUTDIR")
    main(sys.argv[1])
