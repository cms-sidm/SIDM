"""Shared figure conventions for the truth-kinematics AN figures.

Every figure destined for the analysis note goes through these helpers so
sizing, labels, and output format stay consistent. Truth-level figures carry
the CMS Simulation label; save as vector PDF into figures/.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

FIG_DIR = os.path.join(os.path.dirname(__file__), "figures")

# single-column and full-width sizes for the cms-tdr note layout
SINGLE = (10, 8)
WIDE = (16, 8)


def set_style():
    plt.style.use(hep.style.CMS)
    plt.rcParams["savefig.bbox"] = "tight"


def cms_sim_labels(axes, com=13):
    """CMS Simulation label + sqrt(s) on every panel of a finished figure.

    Call this once, after every panel, colorbar and artist exists. The canvas
    is drawn so constrained layout settles with all of them in place, the
    geometry is then frozen, and only then are the labels placed: placing a
    label on a later panel must not reflow the layout and shift labels already
    placed on earlier ones. Labelling panel by panel as a figure is built
    freezes the layout too early, and a colorbar added afterwards gets no
    space allocated, so its label lands on the neighbouring axis."""
    axes = list(np.ravel(np.asarray(axes, dtype=object)))
    fig = axes[0].figure
    fig.canvas.draw()
    if fig.get_layout_engine() is not None:
        fig.set_layout_engine("none")
    for ax in axes:
        hep.cms.label(ax=ax, data=False, com=com)


def cms_sim_label(ax, com=13):
    """Single-panel form of cms_sim_labels."""
    cms_sim_labels([ax], com=com)


def save(fig, name):
    """Save an AN figure as vector PDF with a deterministic filename."""
    os.makedirs(FIG_DIR, exist_ok=True)
    path = os.path.join(FIG_DIR, f"{name}.pdf")
    fig.savefig(path)
    print(f"saved {path}")
    return path
