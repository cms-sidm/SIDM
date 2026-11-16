#!/usr/bin/env python3
"""Auxiliary numbers for the segment-matched transfer-factor rerun.

(1) control-sample sizes: dR-matched vs segment-matched, same files, per channel
(2) segment-match partner multiplicity and the pT_PF tie-break spread on the feed-up rate
(3) beam-halo diagnostics for the mu-LJ DSA tail: |dxyBS|, CSC-only fraction, phi x |eta|
No plotting, read-only.  Writes JSON to --json.
"""
import argparse, glob, json, math, os
import numpy as np
from coffea.util import load

# Repository root; relative paths given on the command line are taken from here.
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
DATA_DIR = "RunOutputFiles/vtxResp2_2018Full"     # --data-dir
MC_DIR = "RunOutputFiles/vtxRespMC2_2018"     # --mc-dir
CH = {"Before": "cr_1egmOr1mu_lj",
      "CR": "cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu",
      "VR": "cr_1egmOr1mu_lj_vtxChi2_CosMu"}
WANT = ["dsaPfMatch_ptPf_ptDsa", "dsaSegMatch_ptPf_ptDsa", "dsaSegMatch_ptPfFewest_ptDsa",
        "dsaSegMatch_ptPf_ptDsa_nPartners",
        "mu_lj_dsaMuon_pt_nHits_absdxyBS", "mu_lj_dsaMuon_pt_nHits_nDTHits",
        "mu_lj_dsaMuon_pt_phi_abseta",
        "dsaSegMatch_ptDsa_absdxyBS_nHits", "dsaSegMatch_ptDsa_nDTHits_nHits"]

# flow-view index conventions (verified against the axis edges at run time)
PT_LOW, PT_SB, PT_TAIL = [1], [2], [3, 4]          # [0,100) ; [100,200) ; >=200 incl OF
NH_LOW, NH_HIGH = [1], [2, 3, 4]                   # [0,20) ; >=20 incl OF
PHI_H = [1, 8, 9, 16]                              # horizontal: phi near 0 and +-pi
PHI_V = [4, 5, 12, 13]                             # vertical:   phi near +-pi/2


class S:
    def __init__(self):
        self.v, self.ax, self.chan, self.n = {}, {}, None, 0
        self.names = []

    def add(self, path):
        obj = load(path)
        out = obj["out"]
        base = os.path.basename(path)[:-len(".coffea")]
        rec = out[base]
        hists = rec["hists"]
        self.names.append(base)
        self.nh = len(hists)
        chans = list(hists["mu_lj_dsaMuon_pt"].axes[0])
        if self.chan is None:
            self.chan = {c: i for i, c in enumerate(chans)}
        for hn in WANT:
            h = hists[hn]
            v = np.asarray(h.values(flow=True), dtype=float)
            a = [np.asarray(x.edges, dtype=float) for x in h.axes[1:]]
            if hn not in self.v:
                self.v[hn], self.ax[hn] = v.copy(), a
            else:
                self.v[hn] += v
        self.n += 1
        del obj

    def get(self, hn, ch):
        return self.v[hn][self.chan[ch]]


def frac(num, den):
    return float(num) / float(den) if den > 0 else float("nan")


def binom_err(p, n):
    if n <= 0 or not np.isfinite(p):
        return float("nan")
    return math.sqrt(max(p * (1 - p), 0.0) / n)


def halo_block(st, ch, tag):
    """Three regions on the (pT, nHits) plane, three diagnostics each."""
    dxy = st.get("mu_lj_dsaMuon_pt_nHits_absdxyBS", ch)      # pt x nHits x |dxyBS|
    ndt = st.get("mu_lj_dsaMuon_pt_nHits_nDTHits", ch)       # pt x nHits x nDTHits
    phe = st.get("mu_lj_dsaMuon_pt_phi_abseta", ch)          # pt x phi x |eta|
    nhx = st.get("mu_lj_dsaMuon_pt_nHits_nDTHits", ch)
    regions = {
        "tail_lowNH":  (PT_TAIL, NH_LOW),
        "tail_highNH": (PT_TAIL, NH_HIGH),
        "low_lowNH":   (PT_LOW,  NH_LOW),
        "low_highNH":  (PT_LOW,  NH_HIGH),
        "SB_lowNH":    (PT_SB,   NH_LOW),
    }
    o = {}
    for rn, (pi, ni) in regions.items():
        d = dxy[np.ix_(pi, ni)].sum(axis=(0, 1))             # |dxyBS| flow vector, len 8
        t = ndt[np.ix_(pi, ni)].sum(axis=(0, 1))             # nDTHits flow vector, len 7
        N = float(d.sum())
        Nt = float(t.sum())
        # |dxyBS| bins: UF,[0,.1),[.1,.5),[.5,2),[2,10),[10,50),[50,200),OF
        f_gt2 = frac(d[4:].sum(), N)
        f_gt10 = frac(d[5:].sum(), N)
        f_lt01 = frac(d[1], N)
        # nDTHits bins: UF,[0,1)=CSC-only,[1,6),[6,12),[12,24),[24,60),OF
        f_csc = frac(t[1], Nt)
        o[rn] = {"N": N, "f_absdxyBS_lt0p1": f_lt01, "f_absdxyBS_gt2": f_gt2,
                 "f_absdxyBS_gt10": f_gt10, "e_gt2": binom_err(f_gt2, N),
                 "N_ndt": Nt, "f_cscOnly": f_csc, "e_cscOnly": binom_err(f_csc, Nt),
                 "dxy_vec": [float(x) for x in d], "ndt_vec": [float(x) for x in t]}
    # phi x |eta| has no nHits axis: quote it in pT slices only
    for rn, pi in (("tail", PT_TAIL), ("SB", PT_SB), ("low", PT_LOW)):
        m = phe[pi].sum(axis=0)                              # phi x |eta| flow
        inner = m[1:17, 1:4]
        N = float(inner.sum())
        H = float(inner[[i - 1 for i in PHI_H], :].sum())
        V = float(inner[[i - 1 for i in PHI_V], :].sum())
        eta = inner.sum(axis=0)
        o[f"phi_{rn}"] = {"N": N, "H": H, "V": V, "HoverV": frac(H, V),
                          "f_H": frac(H, N), "e_f_H": binom_err(frac(H, N), N),
                          "eta_frac": [frac(x, N) for x in eta],
                          "eta_counts": [float(x) for x in eta]}
        # forward endcap (|eta|>1.2) x horizontal phi -- the halo corner
        cor = float(inner[[i - 1 for i in PHI_H], 2].sum())
        o[f"phi_{rn}"]["f_haloCorner"] = frac(cor, N)
    o["tag"] = tag
    return o


def main():
    global DATA_DIR, MC_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="sidm/studies/vertexChi2/transfer_factor/full/tf_seg_aux.json")
    ap.add_argument("--data-dir", default=DATA_DIR)
    ap.add_argument("--mc-dir", default=MC_DIR)
    args = ap.parse_args()
    DATA_DIR, MC_DIR = (os.path.join(REPO, p) for p in (args.data_dir, args.mc_dir))
    args.json = os.path.join(REPO, args.json)
    os.makedirs(os.path.dirname(args.json), exist_ok=True)
    R = {}
    D, M = S(), S()
    dps = sorted(glob.glob(f"{DATA_DIR}/*.coffea"))
    mps = sorted(glob.glob(f"{MC_DIR}/*.coffea"))
    print(f"data files: {len(dps)}   MC files: {len(mps)}")
    for p in dps:
        D.add(p)
    for p in mps:
        M.add(p)
    print(f"loaded data {D.n}  MC {M.n}   hists per file: {D.nh}")
    R["load"] = {"n_data": D.n, "n_mc": M.n, "hists_per_file": D.nh,
                 "data_samples": D.names, "mc_samples": M.names}

    # axis sanity
    e = D.ax["dsaSegMatch_ptPf_ptDsa"]
    print("matrix ptPf edges:", e[0])
    print("matrix ptDsa edges:", e[1])
    print("LJ cube pt edges:", D.ax["mu_lj_dsaMuon_pt_nHits_absdxyBS"][0],
          " nHits:", D.ax["mu_lj_dsaMuon_pt_nHits_absdxyBS"][1],
          " dxyBS:", D.ax["mu_lj_dsaMuon_pt_nHits_absdxyBS"][2])
    print("nDTHits edges:", D.ax["mu_lj_dsaMuon_pt_nHits_nDTHits"][2])
    print("phi edges:", D.ax["mu_lj_dsaMuon_pt_phi_abseta"][1])
    print("abseta edges:", D.ax["mu_lj_dsaMuon_pt_phi_abseta"][2])
    assert np.allclose(D.ax["dsaSegMatch_ptPf_ptDsa"][0],
                       D.ax["dsaSegMatch_ptPfFewest_ptDsa"][0])
    assert np.allclose(D.ax["dsaSegMatch_ptPf_ptDsa"][0],
                       D.ax["dsaPfMatch_ptPf_ptDsa"][0])

    # ---------------------------------------------------------------- (1) sample sizes
    print("\n=== CONTROL-SAMPLE SIZES (same 68 files, flow-inclusive) ===")
    print(f"{'channel':9s} {'dR-matched':>12s} {'seg-matched':>12s} {'seg/dR':>8s} "
          f"{'dR tail':>9s} {'seg tail':>9s}")
    sizes = {}
    for tag, ch in CH.items():
        a = D.get("dsaPfMatch_ptPf_ptDsa", ch)
        b = D.get("dsaSegMatch_ptPf_ptDsa", ch)
        ta, tb = a[:9, 11:].sum(), b[:9, 11:].sum()
        sizes[tag] = {"dR": float(a.sum()), "seg": float(b.sum()),
                      "ratio": frac(b.sum(), a.sum()),
                      "dR_feedup_num": float(ta), "seg_feedup_num": float(tb),
                      "dR_src": float(a[:9, :].sum()), "seg_src": float(b[:9, :].sum())}
        print(f"{tag:9s} {a.sum():12.0f} {b.sum():12.0f} {frac(b.sum(), a.sum()):8.4f} "
              f"{ta:9.0f} {tb:9.0f}")
    for tag, ch in CH.items():
        a = M.get("dsaPfMatch_ptPf_ptDsa", ch)
        b = M.get("dsaSegMatch_ptPf_ptDsa", ch)
        sizes[tag]["mc_dR"] = float(a.sum())
        sizes[tag]["mc_seg"] = float(b.sum())
        sizes[tag]["mc_ratio"] = frac(b.sum(), a.sum())
    R["control_sizes"] = sizes

    # ---------------------------------------------------------------- (2) tie-break
    print("\n=== SEGMENT-MATCH PARTNER MULTIPLICITY AND pT_PF TIE-BREAK ===")
    tb = {}
    for tag, ch in CH.items():
        np_ = D.get("dsaSegMatch_ptPf_ptDsa_nPartners", ch).sum(axis=(0, 1))
        inner = np_[1:6]
        tot = float(inner.sum())
        f_multi = frac(inner[1:].sum(), tot)
        a = D.get("dsaSegMatch_ptPf_ptDsa", ch)
        f = D.get("dsaSegMatch_ptPfFewest_ptDsa", ch)
        fa = frac(a[:9, 11:].sum(), a[:9, :].sum())
        ff = frac(f[:9, 11:].sum(), f[:9, :].sum())
        d = (ff - fa) / fa if fa > 0 else float("nan")
        tb[tag] = {"nPartners_inner": [float(x) for x in inner], "f_multiPartner": f_multi,
                   "feedup_most": fa, "feedup_fewest": ff, "rel_shift": d}
        print(f"{tag:9s} nPartners 1..5 = {[int(x) for x in inner]}  multi-partner frac "
              f"{f_multi:.4f}")
        print(f"          feed-up T(>=200 | pT_PF<100): most={fa:.5e}  fewest={ff:.5e}  "
              f"rel shift {d:+.4f}")
    R["tiebreak"] = tb

    # ---------------------------------------------------------------- (3) halo
    print("\n=== HALO DIAGNOSTICS (mu-LJ DSA muons) ===")
    halo = {}
    for tag, ch in CH.items():
        halo[f"data/{tag}"] = halo_block(D, ch, f"data/{tag}")
        halo[f"mc/{tag}"] = halo_block(M, ch, f"mc/{tag}")
    R["halo"] = halo
    hdr = (f"{'sample/chan':14s} {'region':12s} {'N':>9s} {'f|dxy|>2cm':>11s} "
           f"{'f|dxy|>10cm':>12s} {'f CSC-only':>11s}")
    print(hdr)
    for k in halo:
        for rn in ("tail_lowNH", "tail_highNH", "SB_lowNH", "low_lowNH", "low_highNH"):
            b = halo[k][rn]
            print(f"{k:14s} {rn:12s} {b['N']:9.1f} {b['f_absdxyBS_gt2']:11.4f} "
                  f"{b['f_absdxyBS_gt10']:12.4f} {b['f_cscOnly']:11.4f}")
    print(f"\n{'sample/chan':14s} {'pT region':10s} {'N':>9s} {'H':>8s} {'V':>8s} "
          f"{'H/V':>7s} {'f_H':>7s} {'eta frac (0-0.9/0.9-1.2/1.2-2.4)':>36s} "
          f"{'f halo corner':>14s}")
    for k in halo:
        for rn in ("phi_tail", "phi_SB", "phi_low"):
            b = halo[k][rn]
            ef = "/".join(f"{x:.3f}" for x in b["eta_frac"])
            print(f"{k:14s} {rn:10s} {b['N']:9.1f} {b['H']:8.1f} {b['V']:8.1f} "
                  f"{b['HoverV']:7.3f} {b['f_H']:7.4f} {ef:>36s} {b['f_haloCorner']:14.4f}")

    # segment-matched control, same corner
    print("\n=== SEGMENT-MATCHED CONTROL, same (pT_DSA, nHits) corner ===")
    ctl = {}
    for nm, ST in (("data", D), ("mc", M)):
        for tag, ch in CH.items():
            dd = ST.get("dsaSegMatch_ptDsa_absdxyBS_nHits", ch)   # ptDsa x dxy x nHits
            tt = ST.get("dsaSegMatch_ptDsa_nDTHits_nHits", ch)    # ptDsa x nDT x nHits
            for rn, (pi, ni) in (("tail_lowNH", (PT_TAIL, NH_LOW)),
                                 ("low_lowNH", (PT_LOW, NH_LOW)),
                                 ("tail_highNH", (PT_TAIL, NH_HIGH))):
                d = dd[np.ix_(pi, range(dd.shape[1]), ni)].sum(axis=(0, 2))
                t = tt[np.ix_(pi, range(tt.shape[1]), ni)].sum(axis=(0, 2))
                N, Nt = float(d.sum()), float(t.sum())
                ctl[f"{nm}/{tag}/{rn}"] = {"N": N, "f_absdxyBS_gt2": frac(d[4:].sum(), N),
                                           "f_cscOnly": frac(t[1], Nt), "N_ndt": Nt}
                print(f"{nm}/{tag:9s} {rn:12s} N={N:9.1f}  f|dxy|>2cm="
                      f"{frac(d[4:].sum(), N):7.4f}  f CSC-only={frac(t[1], Nt):7.4f}")
    R["control_halo"] = ctl

    with open(args.json, "w") as fh:
        json.dump(R, fh, indent=1)
    print("\nJSON written:", args.json)


if __name__ == "__main__":
    main()
