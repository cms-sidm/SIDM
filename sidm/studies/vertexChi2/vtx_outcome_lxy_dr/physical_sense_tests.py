"""Physical-sense tests of the vertex outcome of signal mu-type LJs, from the raw count tables.

usage: python physical_sense_tests.py TABLE_DIR OUT_DIR

Reads counts_gen_dr.csv.gz and counts_lxy.csv (written by vtx_outcome_vs_gen_dr.py) from
TABLE_DIR, channel base_ljDisplacementIso_cosAlpha only, raw unweighted LJ counts; every
fraction is a binomial ratio of counts and is written with its counts. Writes to OUT_DIR:

  TA_points.csv, TA_cells.csv  collapse test. Per cell (coarse generated Lxy bin, generated-dR
      bin) of the gen-matched LJs, and per (M, m_Zd) with lifetimes and both final states
      summed, the no-vertex and pass fractions where the point has >= 200 LJs; per cell the
      range across the points and the chi2 of the hypothesis that all points share one fraction
      (k x 2 Pearson chi2 about the pooled fraction, ndf = points - 1), also with the points
      split by final state. TA_order.csv: in each cell, the Spearman rank correlation of the
      fractions with M at fixed m_Zd (>= 3 points) and with m_Zd at fixed M (>= 2 points), and
      the fractions at the lowest and highest value. TA_final_state.csv: per (cell, M, m_Zd)
      with >= 200 LJs in each final state, 4Mu minus 2Mu2E and the 2 x 2 chi2.
  TB_samples.csv, TB_pooled.csv  Lxy [100, 300) against [300, 1000) cm (counts_lxy.csv, fine Lxy
      axis): per sample with >= 50 LJs in both bins, the three outcome fractions in each bin and
      the difference ([300, 1000) minus [100, 300)) with its binomial error; then the same with
      the counts of those samples summed per (M, m_Zd). TB_final_state.csv: the LJs and the
      no-vertex fraction of the two bins per final state and m_Zd, all samples summed.
  TC_floor.csv  generated Lxy < 10 cm (coarse bins [0, 1) and [1, 10)): the no-vertex fraction
      against the generated dR per bound-state mass M (m_Zd, lifetimes and final states
      summed), for the two Lxy bins together and for each. TC_dr_mix.csv: per M, for
      Lxy < 10 cm, the fraction over all generated dR, the share of the LJs with a generated dR
      below 0.002, the fraction that the dR distribution of that M gives with the fraction
      of each dR bin taken from all masses summed, and the measured fraction minus that, with
      the parts of the difference that come from the dR bin below 0.002 and from the bins
      between 0.002 and 0.01.
  TD_widest.csv  generated dR [0.2, 0.5) against [0.1, 0.2) at a fixed mass point: per coarse
      Lxy bin and (M, m_Zd), lifetimes summed, for both final states together and for each
      (rows with >= 100 LJs in both dR bins), the no-vertex fraction in each dR bin and the
      difference with its binomial error.
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

if len(sys.argv) != 3:
    sys.exit(f"usage: python {sys.argv[0]} TABLE_DIR OUT_DIR")
TAB, OUT = sys.argv[1:3]
REF = "base_ljDisplacementIso_cosAlpha"
CNT = ["n_no_vertex", "n_chi2_ge5", "n_pass", "n_other"]
FR = {"no_vertex": "n_no_vertex", "chi2_ge5": "n_chi2_ge5", "pass": "n_pass"}
N_CELL_MIN, N_TB_MIN, N_TD_MIN = 200, 50, 100
LXY_LO, LXY_HI = "[100, 300)", "[300, 1000)"
DR_NARROWEST, DR_WIDE, DR_WIDEST = "[0, 0.002)", "[0.1, 0.2)", "[0.2, 0.5)"
DR_NEXT = ("[0.002, 0.005)", "[0.005, 0.01)")


def order_of(labels):
    """Sort key of bin labels: underflow first, then by lower edge, overflow last."""
    def key(s):
        if s.startswith("<"):
            return -np.inf
        if s.startswith(">="):
            return np.inf
        return float(s[1:].split(",")[0])
    return {s: key(s) for s in labels}


def binom(k, n):
    """Fraction and binomial error."""
    f = k / n
    return f, np.sqrt(f * (1 - f) / n)


def homogeneity(k, n):
    """Pearson chi2 of a k x 2 table about the pooled fraction, ndf and p-value."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    fb = k.sum() / n.sum()
    if fb in (0.0, 1.0) or len(k) < 2:
        return 0.0, len(k) - 1, 1.0
    chi2 = float(((k - n * fb) ** 2 / (n * fb * (1 - fb))).sum())
    return chi2, len(k) - 1, float(stats.chi2.sf(chi2, len(k) - 1))


G = pd.read_csv(f"{TAB}/counts_gen_dr.csv.gz")
G = G[G.selection == REF].copy()
G["n"] = G[CNT].sum(axis=1)
L = pd.read_csv(f"{TAB}/counts_lxy.csv")
L = L[L.selection == REF].copy()
L["n"] = L[CNT].sum(axis=1)

# ---------------------------------------------------------------- T-A collapse
M = G[(G.lxy_bin != "< 0") & (G.gen_dr_bin != "< 0")]
lk, gk = order_of(M.lxy_bin.unique()), order_of(M.gen_dr_bin.unique())
keys = ["lxy_bin", "gen_dr_bin", "m_med_GeV", "m_Zd_GeV"]
P = M.groupby(keys)[CNT + ["n"]].sum().reset_index()
P = P[P.n >= N_CELL_MIN].copy()
for f, c in FR.items():
    P[f"frac_{f}"], P[f"err_{f}"] = binom(P[c], P.n)
P["_l"], P["_g"] = P.lxy_bin.map(lk), P.gen_dr_bin.map(gk)
P = P.sort_values(["_l", "_g", "m_med_GeV", "m_Zd_GeV"])
P.drop(columns=["_l", "_g"]).to_csv(f"{OUT}/TA_points.csv", index=False, float_format="%.5g")
PF = M.groupby(keys + ["final_state"])[CNT + ["n"]].sum().reset_index()
PF = PF[PF.n >= N_CELL_MIN]
cells, orders = [], []
for (lb, gb), g in P.groupby(["lxy_bin", "gen_dr_bin"]):
    row = dict(lxy_bin=lb, gen_dr_bin=gb, n_points=len(g), n_lj=int(g.n.sum()),
               n_min=int(g.n.min()))
    gf = PF[(PF.lxy_bin == lb) & (PF.gen_dr_bin == gb)]
    for f, c in FR.items():
        if f == "chi2_ge5":
            continue
        fr = g[c] / g.n
        i_lo, i_hi = fr.idxmin(), fr.idxmax()
        chi2, ndf, p = homogeneity(g[c], g.n)
        chi2f, ndff, pf = homogeneity(gf[c], gf.n)
        row.update({f"{f}_pooled": g[c].sum() / g.n.sum(), f"{f}_min": fr.min(),
                    f"{f}_min_at": f"{g.m_med_GeV[i_lo]:g}/{g.m_Zd_GeV[i_lo]:g}",
                    f"{f}_max": fr.max(), f"{f}_max_at": f"{g.m_med_GeV[i_hi]:g}/{g.m_Zd_GeV[i_hi]:g}",
                    f"{f}_chi2": chi2, f"{f}_ndf": ndf, f"{f}_p": p,
                    f"{f}_chi2_fs_split": chi2f, f"{f}_ndf_fs_split": ndff, f"{f}_p_fs_split": pf})
        for fixed, var, nmin in (("m_Zd_GeV", "m_med_GeV", 3), ("m_med_GeV", "m_Zd_GeV", 2)):
            for v, h in g.groupby(fixed):
                if len(h) < nmin:
                    continue
                rho = stats.spearmanr(h[var], h[c] / h.n).statistic
                orders.append(dict(lxy_bin=lb, gen_dr_bin=gb, fraction=f, fixed=fixed,
                                   fixed_value=v, against=var, n_points=len(h), spearman_rho=rho,
                                   frac_at_lowest=(h[c] / h.n).iloc[h[var].argmin()],
                                   frac_at_highest=(h[c] / h.n).iloc[h[var].argmax()]))
    cells.append(row)
TA = pd.DataFrame(cells)
TA["_l"], TA["_g"] = TA.lxy_bin.map(lk), TA.gen_dr_bin.map(gk)
TA = TA.sort_values(["_l", "_g"]).drop(columns=["_l", "_g"])
TA.to_csv(f"{OUT}/TA_cells.csv", index=False, float_format="%.5g")
TO = pd.DataFrame(orders)
TO.to_csv(f"{OUT}/TA_order.csv", index=False, float_format="%.4g")
# final state at fixed (cell, M, m_Zd): 2Mu2E against 4Mu where both have >= 200 LJs
fs_rows = []
for k_, h in PF.groupby(keys):
    if set(h.final_state) == {"2Mu2E", "4Mu"}:
        for f, c in (("no_vertex", "n_no_vertex"), ("pass", "n_pass")):
            chi2, _, _ = homogeneity(h[c], h.n)
            fs_rows.append(dict(zip(keys, k_), fraction=f, chi2=chi2,
                                d=(h[c] / h.n).iloc[list(h.final_state).index("4Mu")]
                                - (h[c] / h.n).iloc[list(h.final_state).index("2Mu2E")]))
FSD = pd.DataFrame(fs_rows)
FSD.to_csv(f"{OUT}/TA_final_state.csv", index=False, float_format="%.5g")

# ---------------------------------------------------------------- T-B downturn
W = L[L.lxy_bin.isin([LXY_LO, LXY_HI])].pivot_table(
    index=["sample", "final_state", "m_med_GeV", "m_Zd_GeV", "ctau_mm"], columns="lxy_bin",
    values=CNT + ["n"], aggfunc="sum")
W = W[(W[("n", LXY_LO)] >= N_TB_MIN) & (W[("n", LXY_HI)] >= N_TB_MIN)]


def downturn(W_):
    """Fractions in the two bins and their difference, from a frame of the pivoted counts."""
    o = pd.DataFrame(index=W_.index)
    o["n_100_300"], o["n_300_1000"] = W_[("n", LXY_LO)], W_[("n", LXY_HI)]
    for f, c in FR.items():
        f1, e1 = binom(W_[(c, LXY_LO)], W_[("n", LXY_LO)])
        f2, e2 = binom(W_[(c, LXY_HI)], W_[("n", LXY_HI)])
        o[f"k_{f}_100_300"], o[f"k_{f}_300_1000"] = W_[(c, LXY_LO)], W_[(c, LXY_HI)]
        o[f"{f}_100_300"], o[f"{f}_300_1000"] = f1, f2
        o[f"d_{f}"], o[f"d_{f}_err"] = f2 - f1, np.sqrt(e1 ** 2 + e2 ** 2)
    o["d_no_pass"], o["d_no_pass_err"] = -o.d_pass, o.d_pass_err
    return o


TBs = downturn(W).reset_index().sort_values(["m_med_GeV", "m_Zd_GeV", "final_state", "ctau_mm"])
TBs.to_csv(f"{OUT}/TB_samples.csv", index=False, float_format="%.5g")
Wp = W.groupby(level=["m_med_GeV", "m_Zd_GeV"]).sum()
TBp = downturn(Wp).reset_index()
TBp.insert(2, "n_samples", W.groupby(level=["m_med_GeV", "m_Zd_GeV"]).size().to_numpy())
TBp.to_csv(f"{OUT}/TB_pooled.csv", index=False, float_format="%.5g")
Wm = W.groupby(level=["m_med_GeV"]).sum()
TBm = downturn(Wm).reset_index()
TBm.insert(1, "n_samples", W.groupby(level=["m_med_GeV"]).size().to_numpy())
TBm.to_csv(f"{OUT}/TB_pooled_by_M.csv", index=False, float_format="%.5g")
TBf = (L[L.lxy_bin.isin([LXY_LO, LXY_HI])]
       .groupby(["final_state", "m_Zd_GeV", "lxy_bin"])[["n", "n_no_vertex"]].sum().reset_index())
TBf["frac_no_vertex"], TBf["err"] = binom(TBf.n_no_vertex, TBf.n)
TBf.to_csv(f"{OUT}/TB_final_state.csv", index=False, float_format="%.5g")

# ---------------------------------------------------------------- T-C floor inside 10 cm
C = G[G.lxy_bin.isin(["[0, 1)", "[1, 10)"]) & (G.gen_dr_bin != "< 0")]
tc = []
for lab, sub in (("[0, 10)", C), ("[0, 1)", C[C.lxy_bin == "[0, 1)"]),
                 ("[1, 10)", C[C.lxy_bin == "[1, 10)"])):
    t = sub.groupby(["m_med_GeV", "gen_dr_bin"])[["n_no_vertex", "n"]].sum().reset_index()
    t["lxy"] = lab
    tc.append(t)
TC = pd.concat(tc)
TC["frac_no_vertex"], TC["err"] = binom(TC.n_no_vertex, TC.n)
TC["_g"] = TC.gen_dr_bin.map(gk)
TC = TC.sort_values(["lxy", "m_med_GeV", "_g"]).drop(columns="_g")
TC = TC[["lxy", "m_med_GeV", "gen_dr_bin", "n", "n_no_vertex", "frac_no_vertex", "err"]]
TC.to_csv(f"{OUT}/TC_floor.csv", index=False, float_format="%.5g")
Ct = C.groupby("m_med_GeV")[["n_no_vertex", "n"]].sum()
print("T-C all generated dR, Lxy < 10 cm:", {m: f"{r.n_no_vertex}/{r.n} = {r.n_no_vertex / r.n:.4f}"
                                             for m, r in Ct.iterrows()})
Cb = C.groupby(["m_med_GeV", "gen_dr_bin"])[["n_no_vertex", "n"]].sum().reset_index()
Cp = Cb.groupby("gen_dr_bin")[["n_no_vertex", "n"]].sum()
Cb["k_from_dr_mix"] = Cb.n * Cb.gen_dr_bin.map(Cp.n_no_vertex / Cp.n)
Cm = Cb.groupby("m_med_GeV")[["n", "n_no_vertex", "k_from_dr_mix"]].sum()
Cm["frac_no_vertex"], Cm["err"] = binom(Cm.n_no_vertex, Cm.n)
Cm["n_dr_lt_0p002"] = Cb[Cb.gen_dr_bin == DR_NARROWEST].set_index("m_med_GeV").n
Cm["share_dr_lt_0p002"] = Cm.n_dr_lt_0p002 / Cm.n
Cm["frac_from_dr_mix"] = Cm.k_from_dr_mix / Cm.n
Cb["k_diff"] = Cb.n_no_vertex - Cb.k_from_dr_mix
Cm["diff_from_dr_mix"] = Cm.frac_no_vertex - Cm.frac_from_dr_mix
Cm["diff_in_dr_lt_0p002"] = Cb[Cb.gen_dr_bin == DR_NARROWEST].set_index("m_med_GeV").k_diff / Cm.n
Cm["diff_in_dr_0p002_0p01"] = (Cb[Cb.gen_dr_bin.isin(DR_NEXT)].groupby("m_med_GeV").k_diff.sum()
                               / Cm.n)
Cm = Cm.reset_index()[["m_med_GeV", "n", "n_no_vertex", "frac_no_vertex", "err", "n_dr_lt_0p002",
                       "share_dr_lt_0p002", "frac_from_dr_mix", "diff_from_dr_mix",
                       "diff_in_dr_lt_0p002", "diff_in_dr_0p002_0p01"]]
Cm.to_csv(f"{OUT}/TC_dr_mix.csv", index=False, float_format="%.5g")

# ---------------------------------------------------------------- T-D the widest generated dR
td = []
for fs in ("both", "2Mu2E", "4Mu"):
    sel = M.gen_dr_bin.isin([DR_WIDE, DR_WIDEST])
    if fs != "both":
        sel = sel & (M.final_state == fs)
    D2 = M[sel].pivot_table(index=["lxy_bin", "m_med_GeV", "m_Zd_GeV"], columns="gen_dr_bin",
                            values=["n", "n_no_vertex"], aggfunc="sum").fillna(0)
    D2 = D2[(D2[("n", DR_WIDE)] >= N_TD_MIN) & (D2[("n", DR_WIDEST)] >= N_TD_MIN)]
    d = pd.DataFrame(index=D2.index)
    for tag, b in (("0p1_0p2", DR_WIDE), ("0p2_0p5", DR_WIDEST)):
        d[f"n_{tag}"] = D2[("n", b)].astype(int)
        d[f"k_no_vertex_{tag}"] = D2[("n_no_vertex", b)].astype(int)
        d[f"no_vertex_{tag}"], d[f"err_{tag}"] = binom(D2[("n_no_vertex", b)], D2[("n", b)])
    d["d_no_vertex"] = d.no_vertex_0p2_0p5 - d.no_vertex_0p1_0p2
    d["d_no_vertex_err"] = np.sqrt(d.err_0p1_0p2 ** 2 + d.err_0p2_0p5 ** 2)
    d = d.reset_index()
    d.insert(0, "final_state", fs)
    td.append(d)
TD = pd.concat(td)
TD["_l"] = TD.lxy_bin.map(lk)
TD["_f"] = TD.final_state.map({"both": 0, "2Mu2E": 1, "4Mu": 2})
TD = TD.sort_values(["_f", "_l", "m_med_GeV", "m_Zd_GeV"]).drop(columns=["_l", "_f"])
TD.to_csv(f"{OUT}/TD_widest.csv", index=False, float_format="%.5g")
print("written to", OUT)
