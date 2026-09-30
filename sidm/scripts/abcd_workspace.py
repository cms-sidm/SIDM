#!/usr/bin/env python3
"""RooParametricHist ABCD workspaces and shape datacards.

Second half of the shape-based chain, and the only part that must run in the
**Combine release** -- ROOT and ``RooParametricHist`` live there, coffea does
not.  Run it directly; :func:`build` is the entry point for importing it.

Structure follows the Combine ABCD tutorial: each control region gets one free
``RooRealVar`` per bin inside a ``RooParametricHist`` with a ``RooAddition`` for
its normalisation, and the signal region's background is not a parameter at all
but a per-bin ``RooFormulaVar`` in those.

Transfer-factor models
----------------------
``global_tf`` (default)
    ``A_i = B_i * (sum C / sum D)``.  Shape from region B, normalisation
    transfer from the integrated ratio, so the factor is measured from every
    event in C and D rather than bin by bin.  This is the one that survives our
    statistics.

``per_bin``
    The literal tutorial, ``A_i = B_i * (C_i / D_i)``.  Needs events in every
    bin of C and D; where they are empty the fit drives ``D_i`` to zero and
    Combine reports ``function value is NAN``.  Kept so the comparison can be
    made rather than asserted.

MC only: ``data_obs`` is filled from simulation, so an unblinded run here is a
closure test and not a measurement.  Expected limits use ``--run blind`` and
never read it.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

REGION_NAMES = {0: "A", 1: "B", 2: "C", 3: "D"}
# Each signal point is tested in the channel targeting its final state, the
# same pairing datacard_tools uses.  A 4Mu point does leak into SR_2mu2e, but
# building a card for that pairing too would double-count it against the
# counting and rateParam results, which pair one point with one channel.
FINAL_STATE_CHANNEL = {"2Mu2E": "SR_2mu2e", "4Mu": "SR_4mu"}
CONTROL_REGIONS = (1, 2, 3)          # B, C, D
SIGNAL_REGION = 0                    # A
LUMI_UNC = 1.025                     # 2.5%, 2018


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def load_payload(path):
    """Read the JSON written by ``shape_tools.export_shapes``."""
    payload = json.loads(Path(path).read_text())
    if not payload.get("mc_only", False):
        raise SystemExit(f"{path}: not marked mc_only; refusing to proceed")
    return payload


def summed_background(background, channel, region):
    """Total background shape in one region, summed over process groups."""
    edges, sumw, sumw2 = None, None, None
    for per_channel in background.values():
        regions = per_channel.get(channel, {})
        # JSON object keys are strings; the region index is an int everywhere else
        block = regions.get(region, regions.get(str(region)))
        if block is None:
            continue
        if edges is None:
            edges = list(block["edges"])
            sumw = [0.0] * len(block["sumw"])
            sumw2 = [0.0] * len(block["sumw2"])
        sumw = [a + b for a, b in zip(sumw, block["sumw"])]
        sumw2 = [a + b for a, b in zip(sumw2, block["sumw2"])]
    if edges is None:
        raise SystemExit(f"no background for {channel} region {region}")
    return edges, sumw, sumw2


def rebin(edges, sumw, sumw2, new_edges):
    """Merge onto coarser edges, folding the tails into the end bins."""
    idx = [min(range(len(edges)), key=lambda j: abs(edges[j] - e)) for e in new_edges]
    for j, e in zip(idx, new_edges):
        if abs(edges[j] - e) > 1e-9:
            raise SystemExit(f"rebin edge {e} is not an existing bin edge")
    out_w, out_w2 = [], []
    for i in range(len(new_edges) - 1):
        lo = idx[i] if i else 0
        hi = idx[i + 1] if i + 1 < len(new_edges) - 1 else len(sumw)
        out_w.append(sum(sumw[lo:hi]))
        out_w2.append(sum(sumw2[lo:hi]))
    return list(new_edges), out_w, out_w2


def to_th1(name, edges, sumw, sumw2=None):
    """A TH1D carrying the given contents, and errors from ``sumw2``."""
    import array
    h = ROOT.TH1D(name, name, len(sumw), array.array("d", edges))
    for i, w in enumerate(sumw, start=1):
        h.SetBinContent(i, w)
        if sumw2 is not None:
            h.SetBinError(i, math.sqrt(max(sumw2[i - 1], 0.0)))
    h.SetDirectory(0)
    return h


# --------------------------------------------------------------------------- #
# Workspace
# --------------------------------------------------------------------------- #
def build_workspace(out_root, channel, edges, bkg, signals, observable,
                    tf_model="global_tf", empty_floor=None, range_factor=2.0):
    """Create the workspace holding the parametric ABCD model for one channel.

    ``bkg`` maps region index onto ``(sumw, sumw2)``; ``signals`` maps a signal
    point name onto the same.  Returns ``(workspace, kept_objects)`` -- the
    second is only there to keep Python references alive while RooFit builds,
    since RooFit does not own the objects passed into its constructors.
    """
    nbins = len(edges) - 1
    ws = ROOT.RooWorkspace("wspace", "wspace")
    x = ROOT.RooRealVar(observable, observable, edges[0], edges[-1])
    keep = [x]

    # Every control-region bin is a free parameter.  A bin with no simulated
    # events would otherwise be pinned at exactly zero, which makes a per-bin
    # C/D transfer factor 0/0; the floor gives such a bin room to float up to
    # a level the simulation cannot exclude rather than freezing the fit.
    if empty_floor is None:
        nonzero = [w for region in CONTROL_REGIONS for w in bkg[region][0] if w > 0]
        empty_floor = (sum(nonzero) / len(nonzero)) if nonzero else 1.0

    bins = {}
    for region in CONTROL_REGIONS:
        name = REGION_NAMES[region]
        sumw, sumw2 = bkg[region]
        args = ROOT.RooArgList()
        per_bin = []
        for i in range(nbins):
            content = float(sumw[i])
            hi = max(range_factor * content, empty_floor)
            var = ROOT.RooRealVar(
                f"bkg_{channel}_{name}_bin{i}",
                f"background in region {name} bin {i}",
                content, 0.0, hi)
            per_bin.append(var)
            args.add(var)
        bins[region] = per_bin
        keep += per_bin + [args]

        shape = ROOT.RooParametricHist(
            f"bkg_{name}", f"background PDF in region {name}", x, args,
            to_th1(f"tmpl_{channel}_{name}", edges, sumw))
        norm = ROOT.RooAddition(
            f"bkg_{name}_norm", f"total background in region {name}", args)
        keep += [shape, norm]
        getattr(ws, "import")(shape, ROOT.RooFit.RecycleConflictNodes())
        getattr(ws, "import")(norm, ROOT.RooFit.RecycleConflictNodes())

    # Region A: not a parameter, a formula in the control-region parameters.
    a_args = ROOT.RooArgList()
    a_bins = []
    if tf_model == "global_tf":
        c_tot = ROOT.RooAddition(f"bkg_C_tot_{channel}", "sum over region C",
                                 _arglist(bins[2]))
        d_tot = ROOT.RooAddition(f"bkg_D_tot_{channel}", "sum over region D",
                                 _arglist(bins[3]))
        tf = ROOT.RooFormulaVar(
            "TF_global", "transfer factor, integrated C/D", "(@0/@1)",
            ROOT.RooArgList(c_tot, d_tot))
        keep += [c_tot, d_tot, tf]
        for i in range(nbins):
            expr = ROOT.RooFormulaVar(
                f"bkg_{channel}_A_bin{i}", f"region A bin {i} = B_i * TF",
                "@0*@1", ROOT.RooArgList(tf, bins[1][i]))
            a_bins.append(expr)
            a_args.add(expr)
    elif tf_model == "per_bin":
        for i in range(nbins):
            tf = ROOT.RooFormulaVar(
                f"TF_{channel}_bin{i}", f"transfer factor C/D bin {i}",
                "(@0/@1)", ROOT.RooArgList(bins[2][i], bins[3][i]))
            expr = ROOT.RooFormulaVar(
                f"bkg_{channel}_A_bin{i}", f"region A bin {i} = B_i * C_i / D_i",
                "@0*@1", ROOT.RooArgList(tf, bins[1][i]))
            keep.append(tf)
            a_bins.append(expr)
            a_args.add(expr)
    else:
        raise SystemExit(f"unknown transfer-factor model {tf_model!r}")
    keep += a_bins + [a_args]

    shape_a = ROOT.RooParametricHist(
        "bkg_A", "background PDF in region A", x, a_args,
        to_th1(f"tmpl_{channel}_A", edges, bkg[SIGNAL_REGION][0]))
    norm_a = ROOT.RooAddition("bkg_A_norm", "total background in region A", a_args)
    keep += [shape_a, norm_a]
    getattr(ws, "import")(shape_a, ROOT.RooFit.RecycleConflictNodes())
    getattr(ws, "import")(norm_a, ROOT.RooFit.RecycleConflictNodes())

    # Observations, from simulation: this is a closure test, not data.
    for region in (SIGNAL_REGION,) + CONTROL_REGIONS:
        name = REGION_NAMES[region]
        th1 = to_th1(f"obs_{channel}_{name}", edges, bkg[region][0])
        data = ROOT.RooDataHist(f"data_obs_{name}", f"MC as pseudo-data, region {name}",
                                ROOT.RooArgList(x), th1, 1.0)
        keep.append(data)
        getattr(ws, "import")(data, ROOT.RooFit.Rename(f"data_obs_{name}"))

    # Signal, entered in all four regions -- it leaks into the control regions
    # at the tens-of-percent level for the high-mass points, and a fit that did
    # not know that would mistake signal for background.
    for point, per_region in signals.items():
        for region in (SIGNAL_REGION,) + CONTROL_REGIONS:
            name = REGION_NAMES[region]
            sumw, sumw2 = per_region[region]
            th1 = to_th1(f"sig_{channel}_{point}_{name}", edges, sumw, sumw2)
            dh = ROOT.RooDataHist(f"{point}_{name}", f"{point}, region {name}",
                                  ROOT.RooArgList(x), th1, 1.0)
            keep.append(dh)
            getattr(ws, "import")(dh, ROOT.RooFit.Rename(f"{point}_{name}"))

    out_root.cd()
    ws.Write()
    return ws, keep


def _arglist(variables):
    args = ROOT.RooArgList()
    for v in variables:
        args.add(v)
    return args


# --------------------------------------------------------------------------- #
# Datacards
# --------------------------------------------------------------------------- #
def write_datacard(path, ws_name, channel, point, bkg, signal, tf_model,
                   normalisation, observable, campaign):
    """One datacard spanning all four ABCD regions of one channel.

    The rate column is ``1`` for the background, whose normalisation comes from
    the ``_norm`` object in the workspace, and ``-1`` for the signal, taken
    from its shape histogram.  Because ``run_combine_limits.py`` sizes its
    signal-strength scan from the rate column, and that column carries no
    yields here, the signal and background levels in region A are recorded in
    the header for it to read instead.
    """
    order = [(1, "B"), (2, "C"), (3, "D"), (0, "A")]
    bin_names = [f"{channel}_{name}" for _, name in order]

    sig_a = sum(signal[SIGNAL_REGION][0])
    bkg_a = sum(bkg[SIGNAL_REGION][0])

    lines = [
        f"# SIDM shape ABCD datacard -- {point}, {channel}",
        f"# campaign {campaign}; observable {observable}; "
        f"transfer-factor model {tf_model}",
        "# background in region A is a per-bin formula in the control-region "
        "parameters, not an input",
        "# MC ONLY: data_obs is simulation, so an unblinded run is a closure "
        "test and not a measurement",
        f"# {normalisation}",
        f"# signal_rate_A = {sig_a:.6g}",
        f"# background_rate_A = {bkg_a:.6g}",
        f"imax {len(order)}  number of bins",
        "jmax 1  number of background processes",
        "kmax *  number of nuisance parameters",
        "-" * 100,
    ]
    # A signal point can have no simulated events at all in a control region.
    # Combine rejects a process whose shape integrates to zero ("Null norm"),
    # and flooring it would invent a yield, so the process is simply left out
    # of that bin -- which is what a zero rate means anyway.
    has_signal = {region: sum(signal[region][0]) > 0 for region, _ in order}

    for region, name in order:
        b = f"{channel}_{name}"
        lines.append(f"shapes  bkg       {b}  {ws_name}  wspace:bkg_{name}")
        if has_signal[region]:
            lines.append(f"shapes  {point}  {b}  {ws_name}  wspace:{point}_{name}")
        lines.append(f"shapes  data_obs  {b}  {ws_name}  wspace:data_obs_{name}")
    lines.append("-" * 100)

    # (bin, process) columns, signal first in each bin that has any
    columns = []
    for region, name in order:
        if has_signal[region]:
            columns.append((f"{channel}_{name}", point, 0))
        columns.append((f"{channel}_{name}", "bkg", 1))

    width = max(max(len(b) for b in bin_names), len(point)) + 2
    col = lambda items: "".join(str(i).ljust(width) for i in items)

    lines += [
        "bin          " + col(bin_names),
        "observation  " + col(["-1"] * len(order)),
        "-" * 100,
        "bin      " + col([b for b, _, _ in columns]),
        "process  " + col([p for _, p, _ in columns]),
        "process  " + col([i for _, _, i in columns]),
        "rate     " + col(["-1" if i == 0 else "1" for _, _, i in columns]),
        "-" * 100,
    ]

    # The background has no rate parameter to attach a normalisation
    # uncertainty to -- it is defined by the free control-region bins -- so the
    # only nuisances are on the signal.
    nuis = [("lumi_13TeV", "lnN",
             [f"{LUMI_UNC}" if i == 0 else "-" for _, _, i in columns])]
    for region, name in order:
        sumw, sumw2 = signal[region]
        total, var = sum(sumw), sum(sumw2)
        rel = math.sqrt(var) / total if total > 0 else 0.0
        if rel <= 0:
            continue
        this_bin = f"{channel}_{name}"
        values = [f"{1 + min(rel, 1.0):.6g}" if (i == 0 and b == this_bin) else "-"
                  for b, _, i in columns]
        nuis.append((f"mcstat_{channel}_{name}_signal", "lnN", values))

    label = max(len(n) for n, _, _ in nuis) + 2
    for name, kind, values in nuis:
        lines.append(name.ljust(label) + kind.ljust(7) + col(values))

    # Combine needs to be told the control-region bins are free parameters of
    # the model rather than constants it should ignore.
    lines.append("-" * 100)
    for region in CONTROL_REGIONS:
        for i in range(len(bkg[region][0])):
            lines.append(f"bkg_{channel}_{REGION_NAMES[region]}_bin{i}  flatParam")

    Path(path).write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def build(shapes, outdir, tf_model="global_tf", rebin="", empty_floor=None,
          channels="", signals="", all_pairings=False):
    """Write one workspace per channel plus a datacard per signal point."""
    args = argparse.Namespace(
        shapes=shapes, outdir=outdir, tf_model=tf_model, rebin=rebin,
        empty_floor=empty_floor, channels=channels, signals=signals,
        all_pairings=all_pairings)

    payload = load_payload(args.shapes)
    background = payload["background"]
    signal_block = payload["signal"]
    observable = payload["observable"]["name"]
    campaign = payload["campaign"]
    normalisation = payload.get("meta", {}).get("normalisation", "signal normalised to 1 fb")

    channels = [c for c in payload["channels"]
                if not args.channels or c in args.channels.split(",")]
    new_edges = [float(e) for e in args.rebin.split(",")] if args.rebin else None

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []

    for channel in channels:
        bkg = {}
        edges = None
        for region in (SIGNAL_REGION,) + CONTROL_REGIONS:
            e, w, w2 = summed_background(background, channel, region)
            if new_edges:
                e, w, w2 = rebin(e, w, w2, new_edges)
            edges = e
            bkg[region] = (w, w2)

        points = {}
        for name, per_channel in signal_block.items():
            if args.signals and name not in args.signals.split(","):
                continue
            prefix = name.split("_")[0]
            if (not args.all_pairings
                    and FINAL_STATE_CHANNEL.get(prefix, channel) != channel):
                continue
            regions = per_channel.get(channel)
            if not regions:
                continue
            block = {}
            for region in (SIGNAL_REGION,) + CONTROL_REGIONS:
                r = regions.get(region) or regions.get(str(region))
                if r is None:
                    block = None
                    break
                e, w, w2 = list(r["edges"]), list(r["sumw"]), list(r["sumw2"])
                if new_edges:
                    e, w, w2 = rebin(e, w, w2, new_edges)
                block[region] = (w, w2)
            if block and sum(block[SIGNAL_REGION][0]) > 0:
                points[name] = block

        if not points:
            print(f"  {channel}: no signal points with a non-zero region A; skipped")
            continue

        ws_path = outdir / f"param_ws_{channel}.root"
        out_root = ROOT.TFile(str(ws_path), "RECREATE")
        build_workspace(out_root, channel, edges, bkg, points, observable,
                        tf_model=args.tf_model, empty_floor=args.empty_floor)
        out_root.Close()

        for point in points:
            card = outdir / f"datacard_shape_{channel}_{point}.txt"
            write_datacard(card, ws_path.name, channel, point, bkg, points[point],
                           args.tf_model, normalisation, observable, campaign)
            written.append(card)
        print(f"  {channel}: {len(edges) - 1} bins, {len(points)} signal points "
              f"-> {ws_path.name}")

    meta = {
        "tf_model": args.tf_model,
        "observable": observable,
        "campaign": campaign,
        "rebin": new_edges,
        "empty_floor": args.empty_floor,
        "n_cards": len(written),
        "mc_only": True,
        "source_shapes": str(Path(args.shapes).resolve()),
    }
    (outdir / "datacards.meta.json").write_text(json.dumps(meta, indent=1))
    print(f"wrote {len(written)} datacards to {outdir}")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--shapes", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--tf-model", default="global_tf", choices=["global_tf", "per_bin"])
    p.add_argument("--rebin", default="")
    p.add_argument("--empty-floor", type=float, default=None)
    p.add_argument("--channels", default="")
    p.add_argument("--signals", default="")
    p.add_argument("--all-pairings", action="store_true")
    a = p.parse_args(argv)
    return build(a.shapes, a.outdir, a.tf_model, a.rebin, a.empty_floor,
                 a.channels, a.signals, a.all_pairings)


if __name__ == "__main__":
    sys.exit(main())
