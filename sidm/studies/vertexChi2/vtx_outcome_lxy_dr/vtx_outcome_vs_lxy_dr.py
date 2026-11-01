"""Vertex outcome of signal mu-type LJs vs gen Lxy and muon-pair dR (vtxLxyDr_signal_2018).

usage: python vtx_outcome_vs_lxy_dr.py [NEW_RUN_DIR] [OUT_DIR]

Reads NEW_RUN_DIR/*.coffea (default <repo>/RunOutputFiles/vtxLxyDr_signal_2018) and the existing
vtxChi2_signal_2018 outputs (oracle O2). Writes to OUT_DIR (default the working directory):
outcome_by_sample.csv, outcome_vs_lxy.csv, outcome_vs_dR.csv, outcome_vs_lxy_dR.csv,
oracles.txt and PNG figures.

Outcome of an LJ, from the vtx_chi2 axis edges [-1.5, -0.5, 0, 1, 2, 5, 10, 100, 1e4, 1e9]
(flow view, index 0 = underflow): no vertex = index 1 (the -1 sentinel); pass = [0, 5) =
indices 3-5; chi2 >= 5 = indices 6-10 (overflow included). Indices 0 and 2 hold no physical
value and are reported as "other" (expected 0).
"""
import glob
import math
import os
import pathlib
import sys

import coffea.util
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pandas as pd
import yaml

from sidm.tools import utilities

matplotlib.use("Agg")

REPO = str(pathlib.Path(__file__).resolve().parents[4])   # repository root
NEW = sys.argv[1] if len(sys.argv) > 1 else f"{REPO}/RunOutputFiles/vtxLxyDr_signal_2018"
OUT = sys.argv[2] if len(sys.argv) > 2 else "."
OLD = f"{REPO}/RunOutputFiles/vtxChi2_signal_2018"
REF = "base_ljDisplacementIso_cosAlpha"
SPREAD = "base_ljDisplacementIso_cosAlpha_ljSpread"
VTX = "base_ljDisplacementIso_cosAlpha_vtxChi2"
SELS = {REF: "reference (no cosmic veto)", SPREAD: "Spread veto"}
NEFF_MIN = 200          # robust subset, as in the write-up (section 4)
NEFF_BIN_MIN = 20       # a sample enters a bin's unweighted mean only with N_eff >= this there
H_LXY, H_DR, H_3D = ("mu_lj_vtxChi2_muonN_genA_lxy", "mu_lj_vtxChi2_muonN_mumu_dR",
                     "mu_lj_vtxChi2_genA_lxy_mumu_dR")
V_NOVTX, V_PASS, V_FAIL, V_OTHER = [1], [3, 4, 5], [6, 7, 8, 9, 10], [0, 2]
OUTCOMES = {"no_vertex": V_NOVTX, "chi2_ge5": V_FAIL, "pass": V_PASS}
MU_1, MU_GE2 = [1], [2, 3, 4]   # muonN axis [1, 2, 3, 11], flow view


def _num(s):
    return float(s.replace("GeV", "").replace("mm", "").replace("p", "."))


def parse(name):
    """4Mu_100GeV_1p2GeV_9p6mm -> ("4Mu", 100.0, 1.2, 9.6); "p" is the decimal point."""
    fs, m, ma, ct = name.split("_")
    return fs, _num(m), _num(ma), _num(ct)


def n_files(meta):
    """Number of input ROOT files recorded in a .meta.yaml sidecar."""
    with open(meta, encoding="utf-8") as fh:
        d = yaml.safe_load(fh)
    return sum(s.get("n_files", 0) for s in d.get("samples", {}).values()) if isinstance(
        d.get("samples"), dict) else d["samples"][0].get("n_files")


def hists(path, with_meta=False):
    """Sample name and hists (and metadata) of one per-sample .coffea output."""
    n = os.path.basename(path)[:-len(".coffea")]
    o = coffea.util.load(path)["out"][n]
    return (n, o["hists"], o["metadata"]) if with_meta else (n, o["hists"])


def yield_scale(n, md):
    """lumi*xs/sumw for an unweighted (condor) output, 1 for a weighted one: the factor
    postprocess would have applied, so pooled sums weight each sample by its yield."""
    if not bool(next(iter(md.get("unweighted_hist", {False})))):
        return 1.0
    return utilities.get_lumixs_weight(n, str(next(iter(md["year"]))), md["scaled_sum_weights"])


def vv(h, ch):
    """Values and variances, flow included, of one channel of a histogram."""
    s = h[{"channel": ch}]
    return s.values(flow=True), s.variances(flow=True)


def bin_labels(edges, over_fmt=">= {:g}"):
    """Labels of the in-range bins plus the overflow."""
    labs = [f"[{a:g}, {b:g})" for a, b in zip(edges[:-1], edges[1:])]
    return labs + [over_fmt.format(edges[-1])]


# ---------------------------------------------------------------- per-sample pass
new_files = sorted(glob.glob(f"{NEW}/*.coffea"))
rows, cubes, o1 = [], {}, []
edges = {}
for f in new_files:
    n, H, md = hists(f, with_meta=True)
    fs, m, ma, ct = parse(n)
    scale = yield_scale(n, md)
    _, Ho = hists(f"{OLD}/{n}.coffea")
    # N_eff as in the write-up (DSA-muon relPtErr histogram, reference channel), from each
    # histogram's own sum w and sum w^2: the new run's decides the robust flag
    v_rp, var_rp = vv(Ho["mu_lj_dsaMuon_pt_relPtErr"], REF)
    neff_old = float(v_rp.sum() ** 2 / var_rp.sum()) if var_rp.sum() > 0 else 0.0
    v_rp2, var_rp2 = vv(H["mu_lj_dsaMuon_pt_relPtErr"], REF)
    neff_new = float(v_rp2.sum() ** 2 / var_rp2.sum()) if var_rp2.sum() > 0 else 0.0
    eff_old = (Ho["mu_lj_pt"][{"channel": VTX}].values(flow=True).sum()
               / Ho["mu_lj_pt"][{"channel": REF}].values(flow=True).sum())
    nf_new = n_files(f[:-len(".coffea")] + ".meta.yaml")
    nf_old = n_files(f"{OLD}/{n}.meta.yaml")
    for ch in SELS:
        a, av = vv(H[H_LXY], ch)            # (vtx, muonN, lxy)
        d, dv = vv(H[H_DR], ch)             # (vtx, muonN, dR)
        c, cv = vv(H[H_3D], ch)             # (vtx, lxy coarse, dR coarse)
        cubes[(n, ch)] = (a * scale, av * scale ** 2, d * scale, dv * scale ** 2,
                          c * scale, cv * scale ** 2)
        if not edges:
            edges = {"lxy": H[H_LXY].axes["genA_lxy"].edges, "dR": H[H_DR].axes["mu_lj_mumu_dR"].edges,
                     "lxy3": H[H_3D].axes["genA_lxy"].edges, "dR3": H[H_3D].axes["mu_lj_mumu_dR"].edges}
        # O1: each new histogram's vtx_chi2 marginal against mu_lj_vtx_chi2 (Regular 52, -2, 50)
        b = H["mu_lj_vtx_chi2"][{"channel": ch}].values(flow=True)
        ref_grp = np.array([b[2], b[3:8].sum(), b[8:].sum(), b[0:2].sum()])
        for hn, arr in ((H_LXY, a), (H_DR, d), (H_3D, c)):
            marg = arr.sum(axis=tuple(range(1, arr.ndim)))
            grp = np.array([marg[V_NOVTX].sum(), marg[V_PASS].sum(), marg[V_FAIL].sum(),
                            marg[V_OTHER].sum()])
            o1.append(dict(sample=n, channel=ch, hist=hn, maxabs=float(np.abs(grp - ref_grp).max()),
                           rel=float(np.abs(grp - ref_grp).max() / max(ref_grp.sum(), 1e-300))))
        tot = a.sum()
        ge2 = a[:, MU_GE2, :].sum()
        out = {k: a[idx][:, MU_GE2, :].sum() / ge2 if ge2 > 0 else np.nan
               for k, idx in OUTCOMES.items()}
        rows.append(dict(
            sample=n, final_state=fs, m_med_GeV=m, m_Zd_GeV=ma, ctau_mm=ct, selection=ch,
            n_eff=neff_new, n_eff_oldrun=neff_old, robust=neff_new >= NEFF_MIN,
            robust_oldrun=neff_old >= NEFF_MIN, raw_or_weighted_mu_lj=tot,
            yield_mu_lj=tot * scale, yield_scale=scale,
            frac_1mu=a[:, MU_1, :].sum() / tot if tot > 0 else np.nan,
            frac_ge2mu=ge2 / tot if tot > 0 else np.nan,
            ge2_frac_no_vertex=out["no_vertex"], ge2_frac_chi2_ge5=out["chi2_ge5"],
            ge2_frac_pass=out["pass"],
            frac_other_vtx_bins=a[V_OTHER].sum() / tot if tot > 0 else np.nan,
            frac_gen_matched=a[:, :, 1:].sum() / tot if tot > 0 else np.nan,
            eff_pass_over_all=a[V_PASS].sum() / tot if tot > 0 else np.nan,
            eff_old_run_vtx_over_ref=eff_old if ch == REF else np.nan,
            n_files_new=nf_new, n_files_old=nf_old,
        ))
S = pd.DataFrame(rows)
S.to_csv(f"{OUT}/outcome_by_sample.csv", index=False, float_format="%.6g")
O1 = pd.DataFrame(o1)

# ---------------------------------------------------------------- oracles
ref = S[S.selection == REF].copy()
ref["o2_diff"] = ref.eff_pass_over_all - ref.eff_old_run_vtx_over_ref
lines = [f"samples read: {len(new_files)} from {NEW}"]
lines.append(f"O1: largest |difference| over {len(O1)} (sample, channel, hist) comparisons "
             f"of the grouped vtx_chi2 marginal: {O1.maxabs.max():.3g} (relative {O1.rel.max():.3g})")
lines.append(f"O1: 'other' vtx bins (underflow, [-0.5, 0)) largest fraction: "
             f"{S.frac_other_vtx_bins.max():.3g}")
lines.append(f"O2: per-sample |pass/all (new) - vtxChi2/ref yield (old)|: max {ref.o2_diff.abs().max():.3g}, "
             f"median {ref.o2_diff.abs().median():.3g}; > 1e-3: {(ref.o2_diff.abs() > 1e-3).sum()}; "
             f"> 1e-2: {(ref.o2_diff.abs() > 1e-2).sum()}")
for _, r in ref[ref.o2_diff.abs() > 1e-2].iterrows():
    lines.append(f"   off by > 1e-2: {r['sample']} new {r.eff_pass_over_all:.4f} old "
                 f"{r.eff_old_run_vtx_over_ref:.4f} files new {r.n_files_new} old {r.n_files_old}")
for _, r in ref[ref.n_files_new != ref.n_files_old].iterrows():
    lines.append(f"   file count differs: {r['sample']} new {r.n_files_new} old {r.n_files_old}")
rb = ref[ref.robust]
lines.append(f"O2: robust subset ({len(rb)} samples) unweighted mean pass/all: new "
             f"{rb.eff_pass_over_all.mean():.4f}, old {rb.eff_old_run_vtx_over_ref.mean():.4f} "
             f"(write-up 0.6897)")
lines.append(f"robust flag from the new run's own sums: {int(ref.robust.sum())} robust; "
             f"differs from the old weighted run's flag for "
             f"{int((ref.robust != ref.robust_oldrun).sum())} samples: "
             + ", ".join(ref[ref.robust != ref.robust_oldrun]['sample']))
lines.append(f"N_eff new (own sums) / old (weighted): median {(ref.n_eff / ref.n_eff_oldrun).median():.4f}, "
             f"min {(ref.n_eff / ref.n_eff_oldrun).min():.4f}, max {(ref.n_eff / ref.n_eff_oldrun).max():.4f}")
rbo = ref[ref.robust_oldrun]
lines.append(f"O2 with the old run's robust flag ({len(rbo)} samples): new "
             f"{rbo.eff_pass_over_all.mean():.4f}, old {rbo.eff_old_run_vtx_over_ref.mean():.4f}")
lines.append(f"frac_1mu: max over samples and selections {S.frac_1mu.max():.3g} "
             "(the channels cut mu_ljs on 'Mu >= 2')")
for grp, col in (("m_med_GeV", "mediator"), ("m_Zd_GeV", "m_Zd")):
    g = rb.groupby(grp).eff_pass_over_all.mean()
    lines.append(f"robust mean pass/all by {col}: " + ", ".join(f"{k:g}: {v:.3f}" for k, v in g.items()))
ctb = pd.cut(rb.ctau_mm, [0, 0.1, 1, 10, 100, np.inf], right=False)
lines.append("robust mean pass/all by ctau [0,0.1,1,10,100,inf) mm: " + ", ".join(
    f"{v:.3f}" for v in rb.groupby(ctb, observed=False).eff_pass_over_all.mean()))
for ex in ("2Mu2E_100GeV_0p25GeV_20p0mm", "2Mu2E_800GeV_0p25GeV_2p5mm", "4Mu_150GeV_5p0GeV_270p0mm",
           "4Mu_1000GeV_5p0GeV_40p0mm"):
    r = ref[ref["sample"] == ex]
    if len(r):
        r = r.iloc[0]
        lines.append(f"example {ex}: no-vertex fraction of all LJs "
                     f"{r.ge2_frac_no_vertex * r.frac_ge2mu:.3f}")
O1.to_csv(f"{OUT}/oracle_O1.csv", index=False, float_format="%.6g")
with open(f"{OUT}/oracles.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))

# ---------------------------------------------------------------- binned tables
GROUPS = {"all": None, **{f"m_med {m:g} GeV": m for m in sorted(S.m_med_GeV.unique())}}


def binned(sel, labels, getter):
    """getter(cube) -> (v, var) of shape (vtx, nbins) restricted to >= 2 muons (and matched)."""
    out = []
    sub = S[(S.selection == sel) & S.robust]
    for glab, gm in GROUPS.items():
        names = sub["sample"] if gm is None else sub[sub.m_med_GeV == gm]["sample"]
        per = {k: [] for k in OUTCOMES}
        pooled_v = 0.0
        pooled_var = 0.0
        for n in names:
            v, w = getter(cubes[(n, sel)])
            pooled_v = pooled_v + v
            pooled_var = pooled_var + w
            tot, totw = v.sum(axis=0), w.sum(axis=0)
            neff = np.divide(tot ** 2, totw, out=np.zeros_like(tot), where=totw > 0)
            for k, idx in OUTCOMES.items():
                fr = np.divide(v[idx].sum(axis=0), tot, out=np.full_like(tot, np.nan), where=tot > 0)
                fr[neff < NEFF_BIN_MIN] = np.nan
                per[k].append(fr)
        for j, lab in enumerate(labels):
            row = dict(selection=sel, group=glab, bin=lab)
            tot = pooled_v[:, j].sum()
            for k, idx in OUTCOMES.items():
                arr = np.array([p[j] for p in per[k]]) if per[k] else np.array([])
                arr = arr[~np.isnan(arr)]
                row[f"{k}_mean"] = arr.mean() if len(arr) else np.nan
                row[f"{k}_sem"] = arr.std(ddof=1) / math.sqrt(len(arr)) if len(arr) > 1 else np.nan
                a = pooled_v[idx, j].sum()
                aw = pooled_var[idx, j].sum()
                f = a / tot if tot > 0 else np.nan
                e = (math.sqrt(max(((1 - f) ** 2 * aw + f ** 2 * (pooled_var[:, j].sum() - aw)), 0))
                     / tot) if tot > 0 else np.nan
                row[f"{k}_pooled"] = f
                row[f"{k}_pooled_err"] = e
            row["n_samples_in_mean"] = int(len(arr))
            row["pooled_yield"] = tot
            tw = pooled_var[:, j].sum()
            row["pooled_neff"] = tot ** 2 / tw if tw > 0 else 0.0
            out.append(row)
    return out


lxy_lab = bin_labels(edges["lxy"])
dr_lab = bin_labels(edges["dR"])
# Lxy: >= 2 muons, matched (lxy flow 1..10); dR: >= 2 muons, dR flow 1..10
def get_lxy(c):
    """(vtx, Lxy) values and variances: >= 2 muons, gen-matched (Lxy flow 1..10)."""
    return c[0][:, MU_GE2, 1:].sum(axis=1), c[1][:, MU_GE2, 1:].sum(axis=1)


def get_dr(c):
    """(vtx, dR) values and variances: >= 2 muons (dR flow 1..10)."""
    return c[2][:, MU_GE2, 1:].sum(axis=1), c[3][:, MU_GE2, 1:].sum(axis=1)


def get_3d(i):
    """(vtx, dR) of coarse Lxy bin i (flow index 1 + i), dR flow 1..5."""
    return lambda c: (c[4][:, 1 + i, 1:], c[5][:, 1 + i, 1:])


T_lxy = pd.DataFrame(sum((binned(s, lxy_lab, get_lxy) for s in SELS), []))
T_dr = pd.DataFrame(sum((binned(s, dr_lab, get_dr) for s in SELS), []))
T_lxy.rename(columns={"bin": "genA_lxy_cm"}).to_csv(f"{OUT}/outcome_vs_lxy.csv", index=False, float_format="%.6g")
T_dr.rename(columns={"bin": "mumu_dR"}).to_csv(f"{OUT}/outcome_vs_dR.csv", index=False, float_format="%.6g")

# coarse 3D: lxy (flow 1..5, matched) x dR (flow 1..5, >= 2 muons), flattened into one bin axis
l3, d3 = bin_labels(edges["lxy3"]), bin_labels(edges["dR3"])
rows3 = []
for s in SELS:
    for i in range(5):
        for r in binned(s, d3, get_3d(i)):
            r["genA_lxy_cm"] = l3[i]
            r["mumu_dR"] = r.pop("bin")
            rows3.append(r)
T3 = pd.DataFrame(rows3)
T3 = T3[["selection", "group", "genA_lxy_cm", "mumu_dR"] + [c for c in T3.columns if c not in
        ("selection", "group", "genA_lxy_cm", "mumu_dR")]]
T3.to_csv(f"{OUT}/outcome_vs_lxy_dR.csv", index=False, float_format="%.6g")

# ---------------------------------------------------------------- figures
plt.style.use(hep.style.CMS)
COL = {"no_vertex": "#c44e52", "chi2_ge5": "#dd8452", "pass": "#4c72b0"}
LAB = {"no_vertex": "no stored vertex", "chi2_ge5": r"vertex, $\chi^2$/ndf $\geq$ 5",
       "pass": r"vertex, $\chi^2$/ndf < 5"}


def stack(ax, T, labels, col_suffix, fs):  # pylint: disable=invalid-name
    """Stacked outcome fractions of one table slice, with the pooled marker."""
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    for k in OUTCOMES:
        y = T[f"{k}_{col_suffix}"].to_numpy(float)
        ax.bar(x, np.nan_to_num(y), bottom=bottom, color=COL[k], width=0.85, label=LAB[k])
        if col_suffix == "mean":
            ax.errorbar(x, bottom + np.nan_to_num(y), T[f"{k}_sem"].to_numpy(float), fmt="none",
                        ecolor="k", capsize=2, lw=1)
        bottom = bottom + np.nan_to_num(y)
    p = T["pass_pooled"].to_numpy(float).copy()
    p[T["pooled_neff"].to_numpy(float) < NEFF_BIN_MIN] = np.nan
    ax.plot(x, 1 - p, "kD", ms=7, label="yield-weighted: no vertex + $\\chi^2$/ndf $\\geq$ 5")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=fs - 2)
    ax.set_ylim(0, 1.17)
    ax.set_yticks(np.linspace(0, 1, 6))


def panel_title(ax, text, fs):
    """Title drawn inside the top of the axes, clear of the CMS label."""
    ax.text(0.5, 0.975, text, transform=ax.transAxes, ha="center", va="top", fontsize=fs)


for var, T, labels, xl in (("lxy", T_lxy, lxy_lab, r"gen $Z_d$ $L_{xy}$ (cm), matched LJs"),
                           ("dR", T_dr, dr_lab, r"$\Delta R$ of the two leading LJ muons")):
    for s, slab in SELS.items():
        fig, axs = plt.subplots(2, 3, figsize=(15, 11.5), layout="constrained", sharey=True, sharex=True)
        for ax, g in zip(axs.flat, [g for g in GROUPS if g != "all"]):
            t = T[(T.selection == s) & (T.group == g)]
            stack(ax, t, labels, "mean", 18)
            panel_title(ax, g.replace("m_med", "mediator"), 17)
            empty = t.n_samples_in_mean.to_numpy() == 0
            for xi in np.where(empty)[0]:
                ax.text(xi, 0.5, "n/a", ha="center", fontsize=14, rotation=90)
        for ax in axs[:, 0]:
            ax.set_ylabel("fraction of LJs (>= 2 muons)", fontsize=17)
        for ax in axs[1]:
            ax.set_xlabel(xl, fontsize=17)
        h_, l_ = axs[0, 0].get_legend_handles_labels()
        fig.legend(h_, l_, loc="outside lower center", ncol=4, fontsize=16,
                   title=f"{slab}; robust samples; bars: unweighted mean over samples "
                         f"(N_eff >= {NEFF_BIN_MIN} in the bin) +- s.e.m.", title_fontsize=16)
        hep.cms.label(text="Simulation Work in Progress", ax=axs[0, 0], fontsize=16, data=True,
                      rlabel="")
        tag = "ref" if s == REF else "spread"
        fig.savefig(f"{OUT}/outcome_vs_{var}_by_mass_{tag}.png", dpi=100)
        plt.close(fig)
    # overall: unweighted mean and pooled, both selections
    fig, axs = plt.subplots(1, 2, figsize=(15, 8.3), layout="constrained", sharey=True)
    for ax, (s, slab) in zip(axs, SELS.items()):
        stack(ax, T[(T.selection == s) & (T.group == "all")], labels, "mean", 18)
        panel_title(ax, f"all robust samples, {slab}", 17)
        ax.set_xlabel(xl, fontsize=17)
    axs[0].set_ylabel("fraction of LJs (>= 2 muons)", fontsize=17)
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="outside lower center", ncol=2, fontsize=16,
               title=f"bars: unweighted mean over samples (N_eff >= {NEFF_BIN_MIN} in the bin) "
                     "+- s.e.m.", title_fontsize=16)
    hep.cms.label(text="Simulation Work in Progress", ax=axs[0], fontsize=16, data=True, rlabel="")
    fig.savefig(f"{OUT}/outcome_vs_{var}_overall.png", dpi=100)
    plt.close(fig)
print("wrote tables and figures to", OUT)
