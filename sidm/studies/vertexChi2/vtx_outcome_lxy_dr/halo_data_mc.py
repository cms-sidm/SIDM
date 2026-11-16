"""Part B.2: halo-sensitive fractions of all mu-LJ DSA muons in the 1-LJ control region,
data (vtxResp2_2018Full, 68 DoubleMuon 2018 skims) and MC (vtxRespMC2_2018, 8 samples, summed
with their lumi*xs weights; fractions within the MC set only, no data/MC normalisation).
Writes OUT_DIR/halo_fractions_data_mc.csv and one PNG per quantity (usage: halo_data_mc.py [OUT_DIR])."""
import glob
import os
import pathlib
import sys

import coffea.util
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from halo_common import HALO_HISTS, HP, halo_rows, PT, NH, ETA  # pylint: disable=wrong-import-position

matplotlib.use("Agg")

RUN = str(pathlib.Path(__file__).resolve().parents[4] / "RunOutputFiles")
W = sys.argv[1] if len(sys.argv) > 1 else "."   # output directory
CHS = {"cr_1egmOr1mu_lj": "1-LJ CR, no veto",
       "cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu": "1-LJ CR, Spread veto"}


def summed(tag):
    """Halo histograms of the two channels summed over every output of one run tag."""
    files = sorted(glob.glob(f"{RUN}/{tag}/*.coffea"))
    v, var, edges = {}, {}, None
    for f in files:
        n = os.path.basename(f)[:-7]
        H = coffea.util.load(f)["out"][n]["hists"]
        for hn in HALO_HISTS:
            for ch in CHS:
                h = H[hn][{"channel": ch}]
                v[(hn, ch)] = v.get((hn, ch), 0) + h.values(flow=True)
                var[(hn, ch)] = var.get((hn, ch), 0) + h.variances(flow=True)
        edges = H[HP].axes["mu_lj_dsaMuon_phi"].edges
    return files, v, var, edges


rows = []
for src, tag in (("data", "vtxResp2_2018Full"), ("MC", "vtxRespMC2_2018")):
    files, v, var, edges = summed(tag)
    print(f"{src}: {len(files)} files from {tag}")
    for ch, lab in CHS.items():
        for r in halo_rows({hn: v[(hn, ch)] for hn in HALO_HISTS},
                           {hn: var[(hn, ch)] for hn in HALO_HISTS}, edges):
            rows.append(dict(source=src, channel=ch, **r))
df = pd.DataFrame(rows)
wide = df.pivot_table(index=["channel", "quantity", "pt", "nhits", "eta"], columns="source",
                      values=["value", "err", "n", "neff"], aggfunc="first")
wide.columns = [f"{b}_{a}" for a, b in wide.columns]
wide = wide.reset_index()[["channel", "quantity", "pt", "nhits", "eta", "data_value", "data_err",
                           "data_n", "MC_value", "MC_err", "MC_n", "MC_neff"]]
out = f"{W}/halo_fractions_data_mc.csv"
wide.to_csv(out, index=False, float_format="%.6g")
print("wrote", out, len(wide), "rows")

# figures: one per quantity; panels = channels; x = pT bins, markers per nHits (or |eta|) bin
plt.style.use(hep.style.CMS)
QS = {"frac_absdxyBS_ge10": ("fraction with |dxy_BS| >= 10 cm", "nhits", list(NH)),
      "frac_zeroDT_cscOnly": ("fraction with zero DT hits (CSC-only)", "nhits", list(NH)),
      "ratio_H_over_V_phi": ("horizontal / vertical (phi)", "eta", list(ETA))}
ptl = [p for p in PT if p != "pT all"]
cols = ["tab:blue", "tab:orange", "tab:green", "k"]
for q, (ylab, key, cats) in QS.items():
    fig, axs = plt.subplots(1, 2, figsize=(15, 7.8), layout="constrained", sharey=True)
    for ax, (ch, lab) in zip(axs, CHS.items()):
        sub = wide[(wide.channel == ch) & (wide.quantity == q)]
        for j, c in enumerate(cats):
            s = sub[sub[key] == c].set_index("pt").loc[ptl]
            x = np.arange(len(ptl)) + (j - 1.5) * 0.09
            ax.errorbar(x - 0.03, s.data_value, s.data_err, fmt="o", color=cols[j], ms=8,
                        label=f"data, {c}")
            ax.errorbar(x + 0.03, s.MC_value, s.MC_err, fmt="s", mfc="none", color=cols[j],
                        ms=8, label=f"MC, {c}")
        ax.set_xticks(range(len(ptl)))
        ax.set_xticklabels([p.replace("pT ", "") + " GeV" for p in ptl], fontsize=16)
        ax.set_xlabel(r"DSA $\mu$ $p_T$", fontsize=18)
        ax.set_title(lab, fontsize=17)
        if q == "ratio_H_over_V_phi":
            ax.set_yscale("log")
    axs[0].set_ylabel(ylab, fontsize=18)
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="outside lower center", ncol=4, fontsize=15,
               title="filled: data (binomial errors); open: MC (sum w^2 errors)", title_fontsize=15)
    hep.cms.label(data=True, rlabel="2018 (13 TeV)", text="Work in Progress", ax=axs[0], fontsize=16)
    fig.savefig(f"{W}/halo_data_mc_{q}.png", dpi=110)
    print("wrote", f"{W}/halo_data_mc_{q}.png")
