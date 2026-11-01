"""Part B.3: halo-sensitive fractions of the mu-LJ DSA muons in signal MC (vtxLxyDr_signal_2018),
both selections, robust subset (N_eff >= 200 as in the write-up), per mediator-mass group.
pooled = all samples of the group summed, each weighted to its yield at 1 fb (lumi*xs/sumw);
mean/sem = unweighted mean over samples with N_eff >= 20 in that cell, standard error on it."""
import glob
import math
import os
import pathlib
import sys

import coffea.util
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from halo_common import HALO_HISTS, HP, halo_rows  # pylint: disable=wrong-import-position

REPO = str(pathlib.Path(__file__).resolve().parents[4])   # repository root
NEW = sys.argv[1] if len(sys.argv) > 1 else f"{REPO}/RunOutputFiles/vtxLxyDr_signal_2018"
OUT = sys.argv[2] if len(sys.argv) > 2 else "."   # also where outcome_by_sample.csv is read
BYS = pd.read_csv(f"{OUT}/outcome_by_sample.csv")
ROBUST = set(BYS[BYS.robust & (BYS.selection == "base_ljDisplacementIso_cosAlpha")]["sample"])
SELS = ["base_ljDisplacementIso_cosAlpha", "base_ljDisplacementIso_cosAlpha_ljSpread"]
KEY = ["quantity", "pt", "nhits", "eta"]

per, pooled, edges = [], {}, None
for f in sorted(glob.glob(f"{NEW}/*.coffea")):
    n = os.path.basename(f)[:-7]
    if n not in ROBUST:
        continue
    m = float(n.split("_")[1].replace("GeV", ""))
    o = coffea.util.load(f)["out"][n]
    H = o["hists"]
    sc = float(BYS[(BYS["sample"] == n)].yield_scale.iloc[0])   # lumi*xs/sumw for unweighted runs
    edges = H[HP].axes["mu_lj_dsaMuon_phi"].edges
    for ch in SELS:
        v = {hn: H[hn][{"channel": ch}].values(flow=True) * sc for hn in HALO_HISTS}
        w = {hn: H[hn][{"channel": ch}].variances(flow=True) * sc ** 2 for hn in HALO_HISTS}
        for r in halo_rows(v, w, edges):
            per.append(dict(sample=n, m_med=m, selection=ch, **r))
        for g in ("all", m):
            key = (ch, g)
            pv, pw = pooled.setdefault(key, ({}, {}))
            for hn in HALO_HISTS:
                pv[hn] = pv.get(hn, 0) + v[hn]
                pw[hn] = pw.get(hn, 0) + w[hn]
P = pd.DataFrame(per)
rows = []
for (ch, g), (pv, pw) in pooled.items():
    sub = P[(P.selection == ch) & ((P.m_med == g) if g != "all" else True)]
    for r in halo_rows(pv, pw, edges):
        s = sub
        for k in KEY:
            s = s[s[k] == r[k]]
        s = s[(s.neff >= 20) & s.value.notna()]
        rows.append(dict(selection=ch, group=("all" if g == "all" else f"m_med {g:g} GeV"),
                         **{k: r[k] for k in KEY}, pooled=r["value"], pooled_err=r["err"],
                         pooled_neff=r["neff"], mean=s.value.mean() if len(s) else np.nan,
                         sem=s.value.std(ddof=1) / math.sqrt(len(s)) if len(s) > 1 else np.nan,
                         n_samples=len(s)))
T = pd.DataFrame(rows).sort_values(["selection", "group", "quantity", "pt", "nhits", "eta"])
T.to_csv(f"{OUT}/halo_fractions_signal.csv", index=False, float_format="%.6g")
print(f"robust samples: {P['sample'].nunique()}; wrote {OUT}/halo_fractions_signal.csv ({len(T)} rows)")
