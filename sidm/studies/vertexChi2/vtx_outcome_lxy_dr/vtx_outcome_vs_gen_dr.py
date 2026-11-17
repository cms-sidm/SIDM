"""Vertex outcome of signal mu-type LJs vs the generated dR of the matched dark photon's muons.

usage: python vtx_outcome_vs_gen_dr.py GENDR_RUN_DIR LXYDR_RUN_DIR OUT_DIR

GENDR_RUN_DIR holds the merged per-sample outputs of vtxGenDr_signal_2018 (hist collection
mu_lj_vtx_outcome_gen_dr), LXYDR_RUN_DIR those of vtxLxyDr_signal_2018. Both productions book
unweighted histograms, so every count written here is a raw LJ count. The robust flag of each
sample (N_eff >= 200) is read from outcome_by_sample.csv next to this script, which
vtx_outcome_vs_lxy_dr.py writes from vtxLxyDr_signal_2018; the new collection has no
DSA-muon histogram to compute it from.

Writes to OUT_DIR:
  closure_gen_dr_marginal.csv  per sample and channel, largest |difference| of two closures:
      "all_ljs": the grouped vtx_chi2 outcome of mu_lj_vtxChi2_genA_lxy_genMuMu_dR summed over
      both gen axes (flow included) against mu_lj_vtx_chi2; "matched_by_lxy": that histogram
      summed over the generated dR against mu_lj_vtxChi2_genA_lxy_mumu_dR summed over the
      reconstructed dR, bin for bin in (vtx_chi2, coarse Lxy) with flow; the gen-matched LJs
      are those with genA_lxy >= 0, the underflow holds the unmatched ones.
  counts_gen_dr.csv.gz  sample x channel x coarse Lxy bin x generated-dR bin, outcome counts.
      Underflow and overflow are their own bins: Lxy "< 0" = no matched dark photon; generated
      dR "< 0" = no matched dark photon, or a matched one with fewer than two daughters.
  counts_lxy.csv  from vtxLxyDr_signal_2018 (mu_lj_vtxChi2_muonN_genA_lxy): sample x channel x
      fine Lxy bin (the ten bins of outcome_vs_lxy.csv plus the "< 0" unmatched bin), LJs with
      >= 2 muons.
  counts_lxy_muonN.csv  the same histogram for Lxy [100, 300) and [300, 1000) only, per bin of
      the LJ's muon multiplicity (axis edges 1, 2, 3, 11).
  outcome_vs_gen_dR.csv  in the format of outcome_vs_dR.csv (groups, robust samples, unweighted
      mean / sem over the samples with N_eff >= 20 in the bin, yield-weighted pooled value),
      against the generated dR, gen-matched LJs of every Lxy.

Outcome of an LJ (vtx_chi2 flow index; edges -1.5, -0.5, 0, 1, 2, 5, 10, 100, 1e4, 1e9), as in
vtx_outcome_vs_lxy_dr.py: no vertex = 1-2, pass = 3-5, chi2 >= 5 = 6-10, other = 0 (below -1.5).
"""
import glob
import math
import os
import pathlib
import sys

import coffea.util
import numpy as np
import pandas as pd

from sidm.tools import utilities

HERE = pathlib.Path(__file__).resolve().parent
if len(sys.argv) != 4:
    sys.exit(f"usage: python {sys.argv[0]} GENDR_RUN_DIR LXYDR_RUN_DIR OUT_DIR")
NEW, OLD, OUT = sys.argv[1:4]
REF = "base_ljDisplacementIso_cosAlpha"
SPREAD = "base_ljDisplacementIso_cosAlpha_ljSpread"
SELS = (REF, SPREAD)
H_GEN, H_RECO3, H_VTX = ("mu_lj_vtxChi2_genA_lxy_genMuMu_dR", "mu_lj_vtxChi2_genA_lxy_mumu_dR",
                         "mu_lj_vtx_chi2")
H_LXY = "mu_lj_vtxChi2_muonN_genA_lxy"
OUTCOMES = {"no_vertex": [1, 2], "chi2_ge5": [6, 7, 8, 9, 10], "pass": [3, 4, 5]}
V_OTHER = [0]
COUNT_COLS = {"n_no_vertex": OUTCOMES["no_vertex"], "n_chi2_ge5": OUTCOMES["chi2_ge5"],
              "n_pass": OUTCOMES["pass"], "n_other": V_OTHER}
MU_GE2 = [2, 3, 4]      # muonN axis [1, 2, 3, 11], flow view
NEFF_BIN_MIN = 20       # as in vtx_outcome_vs_lxy_dr.py
NEFF_RTOL = 1e-9
LXY_SPLIT = ("[100, 300)", "[300, 1000)")


def _num(s):
    return float(s.replace("GeV", "").replace("mm", "").replace("p", "."))


def parse(name):
    """4Mu_100GeV_1p2GeV_9p6mm -> ("4Mu", 100.0, 1.2, 9.6); "p" is the decimal point."""
    fs, m, ma, ct = name.split("_")
    return fs, _num(m), _num(ma), _num(ct)


def load(path):
    """Sample name, hists and metadata of one per-sample .coffea output."""
    n = os.path.basename(path)[:-len(".coffea")]
    o = coffea.util.load(path)["out"][n]
    return n, o["hists"], o["metadata"]


def vv(h, ch):
    """Values and variances, flow included, of one channel of a histogram."""
    s = h[{"channel": ch}]
    v = s.values(flow=True)
    w = s.variances(flow=True)
    return v, (v if w is None else w)


def flow_labels(edges):
    """Underflow, in-range bins and overflow."""
    return ([f"< {edges[0]:g}"] + [f"[{a:g}, {b:g})" for a, b in zip(edges[:-1], edges[1:])]
            + [f">= {edges[-1]:g}"])


def as_int(a, what):
    """Unweighted histograms hold integer counts; refuse anything else."""
    r = np.rint(a)
    if not np.allclose(a, r, rtol=0, atol=1e-6):
        raise ValueError(f"{what}: non-integer counts, are these weighted histograms?")
    return r.astype(np.int64)


def production(md):
    """Signal production recorded in the metadata (v10 for both runs read here)."""
    p = md.get("production", set())
    return next(iter(p)) if len(p) == 1 and next(iter(p)) else "v10"


def yield_scale(n, md):
    """lumi*xs*gen_filter_eff/sumw of an unweighted output, as in vtx_outcome_vs_lxy_dr.py."""
    if not bool(next(iter(md.get("unweighted_hist", {False})))):
        raise ValueError(f"{n}: expected an unweighted output")
    return utilities.get_lumixs_weight(n, str(next(iter(md["year"]))), md["scaled_sum_weights"],
                                       production(md))


def outcome_counts(arr):
    """dict of count columns from an array whose first axis is vtx_chi2 (flow view)."""
    return {c: arr[idx].sum(axis=0) for c, idx in COUNT_COLS.items()}


# ---------------------------------------------------------------- new run: closure + counts
BS = pd.read_csv(HERE / "outcome_by_sample.csv")
ROBUST = dict(zip(BS[BS.selection == REF]["sample"], BS[BS.selection == REF].robust))
OLD_SCALE = dict(zip(BS[BS.selection == REF]["sample"], BS[BS.selection == REF].yield_scale))

new_files = sorted(glob.glob(f"{NEW}/*.coffea"))
closure, rows, cubes, meta_rows = [], [], {}, []
edges = {}
for f in new_files:
    n, H, md = load(f)
    fs, m, ma, ct = parse(n)
    scale = yield_scale(n, md)
    meta_rows.append(dict(sample=n, n_evts=int(md["n_evts"]), yield_scale=scale,
                          yield_scale_lxydr=OLD_SCALE.get(n, np.nan)))
    if not edges:
        ax = H[H_GEN].axes
        edges = {"lxy": ax["genA_lxy"].edges, "gdr": ax["genA_mumu_dR"].edges}
    l_lab, g_lab = flow_labels(edges["lxy"]), flow_labels(edges["gdr"])
    for ch in SELS:
        c, cv = vv(H[H_GEN], ch)                      # (vtx, lxy coarse, gen dR)
        c = as_int(c, f"{n} {ch} {H_GEN}")
        r3, _ = vv(H[H_RECO3], ch)                    # (vtx, lxy coarse, reco dR)
        b = H[H_VTX][{"channel": ch}].values(flow=True)  # Regular(52, -2, 50)
        ref_grp = np.array([b[2], b[3:8].sum(), b[8:].sum(), b[0:2].sum()])
        marg = c.sum(axis=(1, 2))
        grp = np.array([marg[OUTCOMES["no_vertex"]].sum(), marg[OUTCOMES["pass"]].sum(),
                        marg[OUTCOMES["chi2_ge5"]].sum(), marg[V_OTHER].sum()])
        closure.append(dict(sample=n, channel=ch, check="all_ljs", hist=f"{H_GEN} vs {H_VTX}",
                            maxabs=float(np.abs(grp - ref_grp).max()),
                            rel=float(np.abs(grp - ref_grp).max() / max(ref_grp.sum(), 1e-300)),
                            n_lj=int(ref_grp.sum())))
        d = np.abs(c.sum(axis=2) - r3.sum(axis=2))
        closure.append(dict(sample=n, channel=ch, check="matched_by_lxy",
                            hist=f"{H_GEN} vs {H_RECO3}",
                            maxabs=float(d.max()), rel=float(d.max() / max(r3.sum(), 1e-300)),
                            n_lj=int(r3.sum())))
        cubes[(n, ch)] = (c, cv, scale)
        oc = outcome_counts(c)
        for i, ll in enumerate(l_lab):
            for j, gl in enumerate(g_lab):
                rows.append(dict(sample=n, final_state=fs, m_med_GeV=m, m_Zd_GeV=ma, ctau_mm=ct,
                                 selection=ch, lxy_bin=ll, gen_dr_bin=gl,
                                 **{k: int(v[i, j]) for k, v in oc.items()}))
CL = pd.DataFrame(closure)
CL.to_csv(f"{OUT}/closure_gen_dr_marginal.csv", index=False, float_format="%.6g")
pd.DataFrame(rows).to_csv(f"{OUT}/counts_gen_dr.csv.gz", index=False)
META = pd.DataFrame(meta_rows)
META["scale_rel_diff"] = (META.yield_scale / META.yield_scale_lxydr - 1).abs()

# ---------------------------------------------------------------- earlier run: fine Lxy counts
rows_l, rows_mu, n_1mu = [], [], 0
old_files = sorted(glob.glob(f"{OLD}/*.coffea"))
for f in old_files:
    n, H, md = load(f)
    fs, m, ma, ct = parse(n)
    h = H[H_LXY]
    lab_l = flow_labels(h.axes["genA_lxy"].edges)
    lab_mu = flow_labels(h.axes["mu_lj_muonN"].edges)
    for ch in SELS:
        a = as_int(vv(h, ch)[0], f"{n} {ch} {H_LXY}")   # (vtx, muonN, lxy)
        n_1mu += int(a[:, [0, 1], :].sum())
        oc = outcome_counts(a[:, MU_GE2, :].sum(axis=1))
        for i, ll in enumerate(lab_l):
            rows_l.append(dict(sample=n, final_state=fs, m_med_GeV=m, m_Zd_GeV=ma, ctau_mm=ct,
                               selection=ch, lxy_bin=ll, **{k: int(v[i]) for k, v in oc.items()}))
        ocm = outcome_counts(a)
        for i, ll in enumerate(lab_l):
            if ll not in LXY_SPLIT:
                continue
            for k, ml in enumerate(lab_mu):
                rows_mu.append(dict(sample=n, final_state=fs, m_med_GeV=m, m_Zd_GeV=ma,
                                    ctau_mm=ct, selection=ch, lxy_bin=ll, muonN_bin=ml,
                                    **{c_: int(v[k, i]) for c_, v in ocm.items()}))
pd.DataFrame(rows_l).to_csv(f"{OUT}/counts_lxy.csv", index=False)
pd.DataFrame(rows_mu).to_csv(f"{OUT}/counts_lxy_muonN.csv", index=False)

# ---------------------------------------------------------------- group table vs generated dR
S = pd.DataFrame([dict(sample=n, m_med_GeV=parse(n)[1]) for n, _ in cubes if _ == REF])
GROUPS = {"all": None, **{f"m_med {m:g} GeV": m for m in sorted(S.m_med_GeV.unique())}}


def get_gen(cube):
    """(vtx, gen dR) values and variances of the gen-matched LJs (Lxy flow 1..), scaled to
    yield, generated-dR flow 1.. (in-range bins and overflow)."""
    c, cv, s = cube
    return c[:, 1:, 1:].sum(axis=1) * s, cv[:, 1:, 1:].sum(axis=1) * s ** 2


def binned(sel, labels, getter):
    """Rows of outcome_vs_dR.csv's format; the arithmetic of vtx_outcome_vs_lxy_dr.binned."""
    out = []
    names_all = [n for n in S["sample"] if ROBUST.get(n, False)]
    for glab, gm in GROUPS.items():
        names = [n for n in names_all if gm is None or parse(n)[1] == gm]
        per = {k: [] for k in OUTCOMES}
        pooled_v, pooled_var = 0.0, 0.0
        for n in names:
            v, w = getter(cubes[(n, sel)])
            pooled_v = pooled_v + v
            pooled_var = pooled_var + w
            tot, totw = v.sum(axis=0), w.sum(axis=0)
            neff = np.divide(tot ** 2, totw, out=np.zeros_like(tot), where=totw > 0)
            for k, idx in OUTCOMES.items():
                fr = np.divide(v[idx].sum(axis=0), tot, out=np.full_like(tot, np.nan),
                               where=tot > 0)
                fr[neff < NEFF_BIN_MIN * (1 - NEFF_RTOL)] = np.nan
                per[k].append(fr)
        for j, lab in enumerate(labels):
            row = dict(selection=sel, group=glab, gen_mumu_dR=lab)
            tot = pooled_v[:, j].sum()
            for k, idx in OUTCOMES.items():
                arr = np.array([p[j] for p in per[k]]) if per[k] else np.array([])
                arr = arr[~np.isnan(arr)]
                row[f"{k}_mean"] = arr.mean() if len(arr) else np.nan
                row[f"{k}_sem"] = arr.std(ddof=1) / math.sqrt(len(arr)) if len(arr) > 1 else np.nan
                a = pooled_v[idx, j].sum()
                aw = pooled_var[idx, j].sum()
                fv = a / tot if tot > 0 else np.nan
                e = (math.sqrt(max(((1 - fv) ** 2 * aw + fv ** 2 * (pooled_var[:, j].sum() - aw)),
                                   0)) / tot) if tot > 0 else np.nan
                row[f"{k}_pooled"] = fv
                row[f"{k}_pooled_err"] = e
            row["n_samples_in_mean"] = int(len(arr))
            row["pooled_yield"] = tot
            tw = pooled_var[:, j].sum()
            row["pooled_neff"] = tot ** 2 / tw if tw > 0 else 0.0
            out.append(row)
    return out


g_lab_in = flow_labels(edges["gdr"])[1:]
T = pd.DataFrame(sum((binned(s, g_lab_in, get_gen) for s in SELS), []))
T.to_csv(f"{OUT}/outcome_vs_gen_dR.csv", index=False, float_format="%.6g")

# ---------------------------------------------------------------- summary
CG = pd.DataFrame(rows)
lines = [f"vtxGenDr samples read: {len(new_files)}; vtxLxyDr samples read: {len(old_files)}"]
for k, g in CL.groupby("check"):
    lines.append(f"closure {k}: "
                 f"{int((g.maxabs == 0).sum())} of {len(g)} (sample, channel) pairs exact; "
                 f"largest |difference| {g.maxabs.max():.6g} (relative {g.rel.max():.3g})")
lines.append(f"'other' vtx_chi2 bin (below -1.5), total over samples and channels: "
             f"{int(CG.n_other.sum())}")
m = CG[(CG.lxy_bin != "< 0") & (CG.gen_dr_bin == "< 0")]
lines.append(f"gen-matched LJs with generated dR in the underflow (< 2 daughters): "
             f"{int(m[['n_no_vertex', 'n_chi2_ge5', 'n_pass', 'n_other']].sum().sum())}")
um = CG[(CG.lxy_bin == "< 0") & (CG.gen_dr_bin != "< 0")]
lines.append(f"unmatched LJs outside the generated-dR underflow (expected 0): "
             f"{int(um[['n_no_vertex', 'n_chi2_ge5', 'n_pass', 'n_other']].sum().sum())}")
lines.append(f"LJs with < 2 muons in vtxLxyDr (both channels, all samples): {n_1mu}")
lines.append(f"yield scale vs outcome_by_sample.csv: largest relative difference "
             f"{META.scale_rel_diff.max():.3g} over {META.scale_rel_diff.notna().sum()} samples")
with open(f"{OUT}/gen_dr_checks.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
