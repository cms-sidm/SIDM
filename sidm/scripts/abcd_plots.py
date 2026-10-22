#!/usr/bin/env python3
"""Every figure in the study, from one CLI.

A thin driver over ``limit_plotting/plot_tools.py``, which holds the drawing
code.  Nothing here is specific to a campaign, method or plane: the comparison
subcommands take ``label=path`` pairs, so the same call draws two campaigns,
five methods, or one of each.

Subcommands
-----------
``results``
    The single-set suite for one ``limits.csv``: Brazil bands across the grid,
    limit against bound state mass with the model cross section, the exclusion
    map, and the kinetic mixing plane.

``compare``
    Any number of limit sets on shared axes -- against lifetime, against bound
    state mass, ratio to a chosen reference, scatter, and exclusion counts.

``closure``
    Region distributions, per-bin closure with its ratio panel, and the shape
    factorisation check, for a shapes JSON.

MC only; every limit is expected and blinded.

Examples
--------
::

    python sidm/scripts/abcd_plots.py results \\
        --limits .../limits_6d_shape_iso_dphi --campaign sixd_inclusive_v1

    python sidm/scripts/abcd_plots.py compare --reference counting \\
        counting=.../limits_6d_counting_iso_iso \\
        "ABCD iso x iso"=.../limits_6d_abcd_iso_iso \\
        "shape m_ljlj"=.../limits_6d_shape_iso_dphi

    python sidm/scripts/abcd_plots.py closure --shapes .../shapes_iso_dphi.json
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

REPO = Path(__file__).resolve().parents[2]
STUDY = REPO / "sidm" / "studies" / "limit_plotting"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(STUDY))

import plot_tools as pt     # noqa: E402


def _pairs(items):
    """``["label=path", ...]`` -> ordered ``{label: path}``."""
    out = {}
    for item in items:
        if "=" not in item:
            raise SystemExit(f"expected label=path, got {item!r}")
        label, _, path = item.partition("=")
        out[label] = path
    return out


def cmd_results(args):
    rows = pt.as_rows(pt.load_limits(Path(args.limits) / "limits.csv"
                                     if Path(args.limits).is_dir() else args.limits))
    outdir = Path(args.outdir or STUDY / "plots" / args.campaign / "results")
    for path in pt.result_set(rows, outdir, args.campaign, prefix=args.prefix):
        print(" ", path)
    pt.summary_table({args.campaign: rows})
    return 0


def cmd_compare(args):
    sets = pt.load_sets(_pairs(args.sets))
    if not sets:
        raise SystemExit("no limit sets found")
    outdir = Path(args.outdir or STUDY / "plots" / "comparison")
    theory = pt.theory_by_mass()
    pt.summary_table(sets)

    states = sorted({r["fs"] for rows in sets.values() for r in rows})
    for fs in states:
        print(" ", pt.compare_vs_lxy(sets, fs, outdir, args.campaign, theory,
                                     name=f"{args.prefix}compare_vs_lxy_{fs}",
                                     bands=args.bands))
    print(" ", pt.compare_vs_mass(sets, outdir, args.campaign, theory,
                                  name=f"{args.prefix}compare_vs_mass"))
    if args.reference and args.reference in sets:
        for against in ("lxy", "m_bound"):
            print(" ", pt.compare_ratio(sets, args.reference, outdir, args.campaign,
                                        against=against,
                                        name=f"{args.prefix}ratio_vs_{against}"))
    if len(sets) == 2:
        a, b = list(sets)
        print(" ", pt.compare_scatter(sets, a, b, outdir, args.campaign,
                                      name=f"{args.prefix}scatter"))
    print(" ", pt.compare_exclusion(sets, outdir, args.campaign,
                                    name=f"{args.prefix}exclusion",
                                    per_mass=args.per_mass))
    return 0


def cmd_closure(args):
    import shape_tools as st
    payload = st.load_shapes(args.shapes)
    grouped = payload["background"]
    campaign = payload["campaign"]
    label = payload["observable"].get("label", "observable")
    subtitle = payload.get("meta", {}).get("selection", "")
    outdir = Path(args.outdir or STUDY / "plots" / campaign / "closure")

    signal_point = args.signal or next(iter(payload["signal"]), None)
    for channel in sorted({ch for per in grouped.values() for ch in per}):
        sig = payload["signal"].get(signal_point, {}).get(channel) if signal_point else None
        print(" ", pt.region_distributions(grouped, channel, outdir, campaign,
                                           signal=sig, label=label,
                                           name=f"{args.prefix}regions_{channel}"))
        print(" ", pt.closure_panel(grouped, channel, outdir, campaign,
                                    subtitle=subtitle, label=label,
                                    name=f"{args.prefix}closure_{channel}"))
        print(" ", pt.shape_factorisation(grouped, channel, outdir, campaign,
                                          label=label,
                                          name=f"{args.prefix}factorisation_{channel}"))
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="")
    p.add_argument("--campaign", default="sixd_inclusive_v1")
    p.add_argument("--prefix", default="", help="filename prefix for the outputs")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("results", help="the single-set suite for one limits.csv")
    r.add_argument("--limits", required=True)
    r.set_defaults(func=cmd_results)

    c = sub.add_parser("compare", help="several limit sets on shared axes")
    c.add_argument("sets", nargs="+", metavar="label=path")
    c.add_argument("--reference", default="",
                   help="set to divide the others by, e.g. counting")
    c.add_argument("--bands", action="store_true", help="draw +/-1 sigma bands")
    c.add_argument("--per-mass", action="store_true",
                   help="split the exclusion counts by bound state mass")
    c.set_defaults(func=cmd_compare)

    cl = sub.add_parser("closure", help="region, closure and factorisation figures")
    cl.add_argument("--shapes", required=True)
    cl.add_argument("--signal", default="", help="signal point to overlay")
    cl.set_defaults(func=cmd_closure)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
