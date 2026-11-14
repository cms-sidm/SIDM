#!/usr/bin/env python3
"""Write the three committed transfer-factor records from the fold's full JSON outputs.

The records keep only the values section 10 of ../vertexChi2CosmicVetoStudy.ipynb reads.
The fold's channel keys are renamed (Before -> no_veto, CR -> spread_veto, VR -> vertex_veto)
and so are its table keys (T4_transfer -> transfer, T7_closure -> closure, inclusive ->
inclusive_closure, T7_sb_transport_failures -> sideband_transport_failures). Key order follows
the full output. Relative paths are taken from the repository root.
"""
import argparse
import json
import os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
TF = "sidm/studies/vertexChi2/transfer_factor"
CHANNEL = {"Before": "no_veto", "CR": "spread_veto", "VR": "vertex_veto"}

CLOSURE_SB = ["obs", "pred", "stat", "sig_circ", "sig_bin", "sig_chan", "R", "sR_full",
              "R_mc", "sR_mc", "neff_mc"]
CLOSURE_TAIL = ["obs", "pred", "stat", "sig_circ", "sig_bin", "sig_chan", "R", "sR_full",
                "obs_mc", "pred_mc", "R_mc", "sR_mc", "neff_mc"]
CLOSURE_NO_VETO = ["obs", "pred", "R", "sR_full"]      # only no_veto/nHits/0/TAIL is kept
PROVENANCE = ["header", "chunks_expected", "L_data_fb", "mc_present", "mc_missing", "data_samples"]
HALO_DATA = ["f_absdxyBS_gt2", "f_absdxyBS_gt10", "e_gt2", "f_cscOnly", "e_cscOnly"]
HALO_MC = ["f_absdxyBS_gt2", "e_gt2", "f_cscOnly", "e_cscOnly"]
HALO_SLICES = ["tail_lowNH", "tail_highNH", "low_lowNH"]


def rename(key):
    ch, rest = key.split("/", 1)
    return f"{CHANNEL[ch]}/{rest}"


def pick(d, keys):
    return {k: d[k] for k in keys}


def seg_record(F):
    closure = {}
    for k, v in F["T7_closure"].items():
        if k.startswith("Before/"):
            if k == "Before/nHits/0/TAIL":
                closure[rename(k)] = pick(v, CLOSURE_NO_VETO)
        else:
            closure[rename(k)] = pick(v, CLOSURE_TAIL if k.endswith("/TAIL") else CLOSURE_SB)
    return {
        "provenance": pick(F["provenance"], PROVENANCE),
        "transfer": {rename(k): [{"mc_neff_num": r["mc_neff_num"]} for r in rows]
                     for k, rows in F["T4_transfer"].items() if not k.startswith("Before/")},
        "inclusive_closure": {CHANNEL[ch]: v for ch, v in F["inclusive"].items() if ch != "Before"},
        "closure": closure,
        "sideband_transport_failures": [dict(f, slice=rename(f["slice"]))
                                        for f in F["T7_sb_transport_failures"]],
    }


def dr_record(F):
    return {
        "closure": {"spread_veto/nHits/0/TAIL": pick(F["T7_closure"]["CR/nHits/0/TAIL"], ["R_mc", "sR_mc"])},
        "sideband_transport_failures": [{"z_stat": f["z_stat"]} for f in F["T7_sb_transport_failures"]],
    }


def aux_record(A):
    cs = A["control_sizes"]["CR"]
    tb = A["tiebreak"]["CR"]
    halo = {}
    for src, keys in (("data", HALO_DATA), ("mc", HALO_MC)):
        halo[f"{src}/spread_veto"] = {s: pick(A["halo"][f"{src}/CR"][s], keys) for s in HALO_SLICES}
    return {
        "control_sizes": {"spread_veto": pick(cs, ["dR", "seg", "ratio", "dR_feedup_num", "seg_feedup_num",
                                                   "dR_src", "seg_src"])},
        "tiebreak": {"spread_veto": pick(tb, ["f_multiPartner", "rel_shift"])},
        "halo": halo,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seg", default=f"{TF}/full/tf_transfer_factor_seg68.json")
    ap.add_argument("--dr", default=f"{TF}/full/tf_transfer_factor_full68.json")
    ap.add_argument("--aux", default=f"{TF}/full/tf_seg_aux.json")
    ap.add_argument("--outdir", default=TF)
    args = ap.parse_args()
    jobs = [(seg_record, args.seg, "tf_transfer_factor_seg68.json"),
            (dr_record, args.dr, "tf_transfer_factor_full68.json"),
            (aux_record, args.aux, "tf_seg_aux.json")]
    for make, src, name in jobs:
        with open(os.path.join(REPO, src)) as fh:
            rec = make(json.load(fh))
        out = os.path.join(REPO, args.outdir, name)
        with open(out, "w") as fh:
            json.dump(rec, fh, indent=1)
            fh.write("\n")
        print(f"{src} -> {out}")


if __name__ == "__main__":
    main()
