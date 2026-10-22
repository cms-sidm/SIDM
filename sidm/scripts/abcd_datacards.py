#!/usr/bin/env python3
"""Extract ABCD yields and shapes, and write Combine datacards.

This builds Combine's *inputs*; it never runs Combine.  The pipeline is

    coffea  --abcd_datacards-->  shapes JSON + datacards
            --abcd_workspace-->  RooParametricHist workspace + shape cards
            --run_combine_limits-->  limits.csv
            --abcd_plots-->  figures

Subcommands:

``extract``
    Read a campaign's coffea files and write a shapes JSON -- per-ABCD-region
    distributions of one observable.  Works with both histogram layouts: the
    pre-binned ``abcd_region`` axis of the earlier productions, and the 6D
    inclusive layout where the regions are sliced from the plane.

``datacards``
    Build counting and four-bin ``rateParam`` ABCD cards from those shapes (or
    straight from the coffea files), at the 1 fb reference and, where a cross
    section exists, normalised to it.

The shape-fit workspace is a separate script, ``abcd_workspace.py``, because it
needs ROOT and so must run in the Combine release while these run in the release
where coffea lives.

MC only throughout: the extraction refuses the data merge outright and the
fail-closed blinding guard in ``datacard_tools`` applies on top.

Examples
--------
::

    # 6D campaign, m_ljlj shapes for both planes, plus counting and ABCD cards
    python sidm/scripts/abcd_datacards.py extract  --campaign sixd_inclusive_v1 \\
        --plane iso_iso,iso_dphi --datacards

    # older pre-binned campaign, two-bin m_ljlj from SR + VR_invMass
    python sidm/scripts/abcd_datacards.py extract --campaign golden_hotspot_iso025_v1 \\
        --observable mljlj_two_bin

    # then, in the Combine release:
    python sidm/scripts/abcd_workspace.py \\
        --shapes .../shapes_6d_iso_dphi.json --outdir .../datacards_shape_iso_dphi
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
STUDY = REPO / "sidm" / "studies" / "limit_plotting"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(STUDY))


# --------------------------------------------------------------------------- #
# extract
# --------------------------------------------------------------------------- #
def cmd_extract(args):
    import datacard_tools as dt
    import shape_tools as st

    dt.use_campaign(args.campaign)
    outdir = dt.campaign_outdir(STUDY, args.campaign)
    outdir.mkdir(parents=True, exist_ok=True)

    def progress(i, n, name):
        if not args.quiet:
            print(f"    [{i + 1}/{n}] {name}", flush=True)

    planes = args.plane.split(",") if args.plane else [None]
    for plane in planes:
        tag = plane or args.observable
        print(f"\n=== {args.campaign}: {tag}")

        if plane:                                    # 6D inclusive layout
            sel = st.Selection(
                args.selection_name, iso_cut=args.iso_cut, dphi_cut=args.dphi_cut,
                mass_floor=args.mass_floor, apply_disp=not args.no_disp,
                disp_legs=args.disp_legs)
            print(f"  selection: {sel.describe()}")
            get = lambda d: st.collect_6d(d, plane, selection=sel, progress=progress)
            observable = st.Observable(
                name=f"ljlj_mass_{plane}", hists={}, axis="ljlj_mass",
                label=r"$m_{LJLJ}$ [GeV]",
                note=f"6D slice, plane {plane}, {sel.describe()}")
            meta_extra = {"plane": plane, "selection": sel.describe()}
        elif args.observable == "mljlj_two_bin":     # pre-binned, SR + VR pair
            get = lambda d: st.mljlj_two_bin(d, progress=progress)
            observable = st.Observable(
                name="mljlj_two_bin", hists={}, axis="lj_lj_invmass",
                label=r"$m_{LJLJ}$ [GeV]",
                note=f"two bins split at {st.MLJLJ_SPLIT_GEV:g} GeV, "
                     f"from the SR and VR_invMass selections")
            meta_extra = {"observable": "mljlj_two_bin"}
        else:                                        # pre-binned, one observable
            if args.observable not in st.OBSERVABLES:
                raise SystemExit(f"unknown observable {args.observable!r}; known: "
                                 f"mljlj_two_bin, {', '.join(st.OBSERVABLES)}")
            if not st.OBSERVABLES[args.observable].available:
                raise SystemExit(st.missing_observable_help(args.observable))
            get = lambda d: st.collect_shapes(d, args.observable, progress=progress)
            observable = st.OBSERVABLES[args.observable]
            meta_extra = {"observable": args.observable}

        print(f"  background: {dt.BKG_DIR}")
        bkg = st.group_shapes(get(dt.BKG_DIR))
        print(f"  signal: {dt.SIGNAL_DIR}")
        signal = get(dt.SIGNAL_DIR)
        print(f"  groups: {', '.join(sorted(bkg))}; {len(signal)} signal points")

        path = outdir / f"shapes_{tag}.json"
        st.export_shapes(path, bkg, signal, observable,
                         meta={"normalisation": "signal normalised to 1 fb at "
                                                "59.83 /fb, so r = sigma / 1 fb",
                               "background_groups": sorted(bkg), **meta_extra})
        print(f"  -> {path.name}")

        for ch in sorted(dt.CHANNELS):
            try:
                y = {r: st.total_background_shape(bkg, ch, r).total for r in range(4)}
            except st.ShapeError:
                continue
            a, b, c, d = (y[r].value for r in range(4))
            pred = f"{b * c / d:.4f}" if d else "undefined"
            print(f"  {ch}: A {a:.4f}  B {b:.4f}  C {c:.4f}  D {d:.4f} "
                  f"| B*C/D {pred}")

        if args.datacards:
            _write_datacards(dt, st, bkg, signal, outdir, tag, args)
    return 0


def _write_datacards(dt, st, bkg_shapes, signal_shapes, outdir, tag, args):
    """Counting and rateParam ABCD cards from the integral of the same shapes.

    Building both from one extraction is what keeps the methods comparable:
    they differ only in what the datacard does with the bins, not in the events
    that went into them.
    """
    from sidm.tools import utilities
    xs = utilities.load_yaml(f"{dt.BASE_DIR}/configs/cross_sections.yaml")

    bkg = st.grouped_yields(bkg_shapes)
    signal = st.integrated(signal_shapes)
    # A point with no cross section is refused rather than silently normalised
    # to the reference, so the theory-normalised set covers only what it can.
    signal_theory = {k: v for k, v in signal.items() if k in xs}
    print(f"  theory-normalisable: {len(signal_theory)} of {len(signal)}")

    for suffix, writer, points, cfg in (
            ("counting", dt.write_datacards, signal, dt.DatacardConfig()),
            ("counting_theory", dt.write_datacards, signal_theory,
             dt.DatacardConfig(use_theory_xs=True)),
            ("abcd", dt.write_abcd_datacards, signal, dt.DatacardConfig()),
            ("abcd_theory", dt.write_abcd_datacards, signal_theory,
             dt.DatacardConfig(use_theory_xs=True))):
        name = f"datacards_{suffix}_{tag}"
        written, warnings = writer(points, bkg, outdir / name, config=cfg)
        note = f", {len(warnings)} warnings" if warnings else ""
        print(f"  {name}: {len(written)} cards{note}")
        for w in warnings[:2]:
            print(f"      {w}")


# --------------------------------------------------------------------------- #
# datacards (from an existing shapes JSON)
# --------------------------------------------------------------------------- #
def cmd_datacards(args):
    import datacard_tools as dt
    import shape_tools as st

    payload = st.load_shapes(args.shapes)
    dt.use_campaign(payload["campaign"])
    outdir = Path(args.outdir or dt.campaign_outdir(STUDY, payload["campaign"]))
    tag = payload["observable"]["name"]
    _write_datacards(dt, st, payload["background"], payload["signal"], outdir, tag, args)
    return 0


# --------------------------------------------------------------------------- #
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("extract", help="coffea -> shapes JSON (+ optional datacards)")
    e.add_argument("--campaign", default="sixd_inclusive_v1")
    e.add_argument("--plane", default="",
                   help="6D layout: comma-separated planes, e.g. iso_iso,iso_dphi")
    e.add_argument("--observable", default="mljlj_two_bin",
                   help="pre-binned layout: an entry in shape_tools.OBSERVABLES, "
                        "or mljlj_two_bin")
    e.add_argument("--datacards", action="store_true",
                   help="also write the counting and ABCD cards")
    e.add_argument("--iso-cut", type=float, default=0.25)
    e.add_argument("--dphi-cut", type=float, default=2.0)
    e.add_argument("--mass-floor", type=float, default=150.0)
    e.add_argument("--no-disp", action="store_true",
                   help="drop the displacement requirements (buys the most stats)")
    e.add_argument("--disp-legs", default="both", choices=["both", "leg0", "leg1"])
    e.add_argument("--selection-name", default="custom")
    e.add_argument("--quiet", action="store_true")
    e.set_defaults(func=cmd_extract)

    d = sub.add_parser("datacards", help="shapes JSON -> counting and ABCD cards")
    d.add_argument("--shapes", required=True)
    d.add_argument("--outdir", default="")
    d.set_defaults(func=cmd_datacards)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
