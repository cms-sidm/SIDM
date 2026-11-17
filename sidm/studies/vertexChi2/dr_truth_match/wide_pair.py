"""Widest generated pairs (gen dR [0.2, 0.5) against [0.1, 0.2)): tables T0 to T5, from the
per-LJ truth-match records.

usage: python wide_pair.py ANALYZE_OUT COUNTS_CSV_GZ OUT_DIR

ANALYZE_OUT holds lj_truth.parquet and lj_classified.parquet written by
dr_truth_match/analyze.py --propagate. COUNTS_CSV_GZ is vtx_outcome_lxy_dr/counts_gen_dr.csv.gz,
the vertex-outcome counts of the condor histograms in bins of sample x Lxy x gen dR.
Writes to OUT_DIR (must exist):
T0_counts_gen_dr_by_final_state.csv (no-vertex fractions of COUNTS_CSV_GZ by final state);
T1_records_vs_counts_gen_dr.csv, T1_records_vs_counts_gen_dr_by_sample.csv (the records binned
as COUNTS_CSV_GZ, cell by cell);
T2_composition.csv; T2_dr_to_other.csv (quantiles of the smallest dR between m0, m1 and a daughter
of the other dark photon);
T3_cells.csv and T3_summary.csv (Lxy [100, 1000)), T3_cells_lxy10_100.csv and
T3_summary_lxy10_100.csv (Lxy [10, 100)); T4_cells.csv / T4_summary.csv (the two variables of T3
closing most of the gap); T4_alt_cells.csv / T4_alt_summary.csv (class x charge); T5_cells.csv /
T5_summary.csv (the [0.2, 0.5) bin standardised by soft_dau_unreco alone, both Lxy ranges);
T5_acceptance_cells.csv (soft_dau_unreco x acceptance); tables.txt.

The gen dR of the records is the dR of the final-state (status-1) daughter muons; the gen dR of
COUNTS_CSV_GZ (hists.mu_lj_gen_a_mumu_dr) is that of the dark photon's own children, before
final-state radiation. T1 compares the two binnings and is exact only where they agree.

Population of T2 to T4: gen-matched mu-type LJs (every LJ has >= 2 muons), raw counts.
Outcome: no vertex = vtx_chi2 < 0. Match radii: the pooled 99 % quantiles analyze.py wrote to
ANALYZE_OUT/table3.csv (rows PF and DSA), i.e. the radii of its classes.
T2 columns soft_dau_unreco (no LJ muon, PF or DSA, within its radius of the softer daughter of
the matched dark photon, propagated direction), acc_eta2p4_pt10 (both daughters |eta| < 2.4 and
the softer pT > 10 GeV) and lead_on_other (m0 or m1 within its radius of a daughter of another
dark photon, propagated from that dark photon's Lxy and phi, from the o1_*, o2_*, o_a_* fields).
Summaries of T3 to T5 carry binomial errors (err_*, the WIDE cell shares held fixed) and
bootstrap errors (*_bs_err: standard deviation over NBOOT replicas that resample the LJs of the
[0.1, 0.2) and of the [0.2, 0.5) bin separately, fixed seed; gap_closed_bs_q16 / _q84: its
16 % and 84 % quantiles).
"""
import sys

import awkward as ak
import numpy as np
import pandas as pd

IN, COUNTS, OUT = sys.argv[1:4]
REF = "base_ljDisplacementIso_cosAlpha"
LXY_EDGES = [0, 1, 10, 100, 1000]
LXY_LAB = ["[0, 1)", "[1, 10)", "[10, 100)", "[100, 1000)"]
DR_EDGES = [0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
DR_LAB = ["[0, 0.002)", "[0.002, 0.005)", "[0.005, 0.01)", "[0.01, 0.02)", "[0.02, 0.05)",
          "[0.05, 0.1)", "[0.1, 0.2)", "[0.2, 0.5)"]
lines = []


def emit(s=""):
    lines.append(str(s))
    print(s)


def flow_bin(x, edges, labs):
    """Histogram-style bin label with underflow '< 0' (value < first edge) and overflow."""
    x = np.asarray(x, dtype=float)
    out = np.empty(len(x), dtype=object)
    out[:] = ""
    out[x < edges[0]] = "< 0"
    out[x >= edges[-1]] = f">= {edges[-1]:g}"
    for lo, hi, lab in zip(edges[:-1], edges[1:], labs):
        out[(x >= lo) & (x < hi)] = lab
    return out


def dphi(a, b):
    return (a - b + np.pi) % (2 * np.pi) - np.pi


def phi_at_pca(phi0, pt, q, lxy_cm, phi_v):
    """As analyze.py: azimuth at the beamline PCA of a helix (B = 3.8 T along +z) through the
    transverse vertex lxy*(cos phi_v, sin phi_v) with azimuth phi0 there."""
    rad = pt / (0.2998 * 3.8)
    vx, vy = lxy_cm / 100 * np.cos(phi_v), lxy_cm / 100 * np.sin(phi_v)
    cx, cy = vx + q * rad * np.sin(phi0), vy - q * rad * np.cos(phi0)
    cn = np.hypot(cx, cy)
    px, py = cx - rad * cx / cn, cy - rad * cy / cn
    rx, ry = px - cx, py - cy
    return np.arctan2(-q * rx, q * ry)


# ------------------------------------------------------------------ T1: records vs COUNTS
R = ak.from_parquet(f"{IN}/lj_truth.parquet",
                    columns=["sample", "vtx_chi2", "gen_matched", "gen_a_lxy", "d1_eta",
                             "d1_phi", "d2_eta", "d2_phi"])
T = pd.DataFrame({k: ak.to_numpy(R[k]) for k in R.fields if k != "sample"})
T["sample"] = np.array(ak.to_list(R["sample"]))
gdr = np.sqrt((T.d1_eta - T.d2_eta) ** 2 + dphi(T.d1_phi, T.d2_phi) ** 2)
T["lxy_v"] = np.where(T.gen_matched, T.gen_a_lxy, -1.0)
T["gdr_v"] = np.where(T.gen_matched, gdr, -1.0)
T["lxy_bin"] = flow_bin(T.lxy_v, LXY_EDGES, LXY_LAB)
T["gen_dr_bin"] = flow_bin(T.gdr_v, DR_EDGES, DR_LAB)
T["nv"] = (T.vtx_chi2 < 0).astype(int)
T["oth"] = 1 - T.nv
mine = T.groupby(["sample", "lxy_bin", "gen_dr_bin"])[["nv", "oth"]].sum().reset_index()

O = pd.read_csv(COUNTS)
O = O[(O.selection == REF) & O["sample"].isin(T["sample"].unique())].copy()
O["nv_ref"] = O.n_no_vertex
O["oth_ref"] = O.n_chi2_ge5 + O.n_pass + O.n_other
C = O[["sample", "lxy_bin", "gen_dr_bin", "nv_ref", "oth_ref"]].merge(
    mine, on=["sample", "lxy_bin", "gen_dr_bin"], how="outer").fillna(0)
C["d_nv"] = C.nv - C.nv_ref
C["d_oth"] = C.oth - C.oth_ref
C.to_csv(f"{OUT}/T1_records_vs_counts_gen_dr.csv", index=False)
d = np.concatenate([C.d_nv.values, C.d_oth.values])
emit("# T1: records against counts_gen_dr.csv.gz (reference channel), cells = sample x Lxy bin "
     "(incl. '< 0' unmatched and '>= 1000') x gen dR bin (incl. flow) x {no vertex, other}")
emit(f"cells compared {len(d)}, exact {int((d == 0).sum())}, largest |diff| {int(np.abs(d).max())}")
core = C[C.lxy_bin.isin(LXY_LAB) & C.gen_dr_bin.isin(DR_LAB)]
dc = np.concatenate([core.d_nv.values, core.d_oth.values])
emit(f"coarse Lxy [0,1000) x gen dR [0,0.5) only: cells {len(dc)}, exact {int((dc == 0).sum())}, "
     f"largest |diff| {int(np.abs(dc).max())}")
bys = C.groupby("sample").agg(n=("nv", "sum"), n_ref=("nv_ref", "sum"),
                              m=("oth", "sum"), m_ref=("oth_ref", "sum"),
                              cells_diff=("d_nv", lambda x: 0)).reset_index()
bys["n_lj"] = bys.n + bys.m
bys["n_lj_ref"] = bys.n_ref + bys.m_ref
bys["cells_nonzero_diff"] = [int(((g.d_nv != 0) | (g.d_oth != 0)).sum())
                             for _, g in C.groupby("sample")]
bys["sum_abs_diff"] = [int(np.abs(g.d_nv).sum() + np.abs(g.d_oth).sum())
                       for _, g in C.groupby("sample")]
bys["max_abs_diff"] = [int(max(np.abs(g.d_nv).max(), np.abs(g.d_oth).max()))
                       for _, g in C.groupby("sample")]
bys = bys.rename(columns={"n": "n_novtx", "n_ref": "n_novtx_ref"}).drop(
    columns=["m", "m_ref", "cells_diff"])
bys.to_csv(f"{OUT}/T1_records_vs_counts_gen_dr_by_sample.csv", index=False)
emit(bys.to_string(index=False))
emit()

# ------------------------------------------------------------------ T0: COUNTS by final state
# straight from counts_gen_dr.csv.gz: m_Zd = 5 GeV, M = 100/150/200, reference channel
A = pd.read_csv(COUNTS)
A = A[(A.selection == REF) & (A.m_Zd_GeV == 5) & A.m_med_GeV.isin([100, 150, 200])
      & A.lxy_bin.isin(["[10, 100)", "[100, 1000)"]) & A.gen_dr_bin.isin(["[0.1, 0.2)",
                                                                        "[0.2, 0.5)"])].copy()
A["n"] = A.n_no_vertex + A.n_chi2_ge5 + A.n_pass + A.n_other
mine9 = set(T["sample"].unique())
rows = []
for pop, sel in (("both final states, all lifetimes", A.index == A.index),
                 ("2Mu2E, all lifetimes", A.final_state == "2Mu2E"),
                 ("4Mu, all lifetimes", A.final_state == "4Mu"),
                 ("the samples of this run", A["sample"].isin(mine9))):
    for m, g in A[sel].groupby("m_med_GeV"):
        r = dict(population=pop, M=int(m))
        for (lx, gb), c in g.groupby(["lxy_bin", "gen_dr_bin"]):
            r[f"{lx} {gb}"] = f"{c.n_no_vertex.sum() / c.n.sum():.3f} ({c.n_no_vertex.sum()}/{c.n.sum()})"
        rows.append(r)
T0 = pd.DataFrame(rows)
T0.to_csv(f"{OUT}/T0_counts_gen_dr_by_final_state.csv", index=False)
emit("# T0: no-vertex fraction (n_novtx/n) from counts_gen_dr.csv.gz, m_Zd = 5 GeV, reference channel")
emit(T0.to_string(index=False))
emit()

# ------------------------------------------------------------------ population
M = pd.read_parquet(f"{IN}/lj_classified.parquet")
M["nv"] = (M.vtx_chi2 < 0).astype(int)
M["lxy_c"] = flow_bin(M.lxy, LXY_EDGES, LXY_LAB)
M["gdr_c"] = flow_bin(M.gen_dr, DR_EDGES, DR_LAB)
M["content"] = np.select([(M.n_mu == 2) & (M.n_pf == 2), (M.n_mu == 2) & (M.n_pf == 1),
                          (M.n_mu == 2) & (M.n_dsa == 2)], ["PF-PF", "PF-DSA", "DSA-DSA"],
                         ">=3mu")
M["nmu"] = np.where(M.n_mu >= 3, ">=3", "2")
M["charge"] = np.where(M.same_charge, "same", "opposite")
PT_EDGES, PT_LAB = [0, 5, 10, 20, np.inf], ["[0,5)", "[5,10)", "[10,20)", ">=20"]
M["soft_reco_pt"] = np.minimum(M.m0_pt, M.m1_pt)
M["soft_gen_pt"] = np.minimum(M.d1_pt, M.d2_pt)
M["soft_reco_pt_bin"] = pd.cut(M.soft_reco_pt, PT_EDGES, right=False, labels=PT_LAB).astype(str)
M["soft_gen_pt_bin"] = pd.cut(M.soft_gen_pt, PT_EDGES, right=False, labels=PT_LAB).astype(str)
eta_pair = np.abs((M.d1_eta + M.d2_eta) / 2)
M["region"] = np.where(eta_pair < 0.9, "barrel", np.where(eta_pair > 1.2, "endcap", "overlap"))
acc = (np.abs(M.d1_eta) < 2.4) & (np.abs(M.d2_eta) < 2.4)
M["acc_pf"] = acc & (M.soft_gen_pt > 5)
M["acc_dsa"] = acc & (M.soft_gen_pt > 10)

# match radii of analyze.py's classes, and the two coverage flags of T2
T3A = pd.read_csv(f"{IN}/table3.csv").set_index("muons")
R_PF, R_DSA = float(T3A.loc["PF", "q99"]), float(T3A.loc["DSA", "q99"])
RT = ak.from_parquet(f"{IN}/lj_truth.parquet")
assert M.index.max() < len(RT) and (ak.to_numpy(RT.vtx_chi2)[M.index.values]
                                    == M.vtx_chi2.values).all(), "row alignment"
RM = RT[M.index.values]
mu = RM.muons
mu_r = ak.where(mu["type"] == 0, R_PF, R_DSA)
soft_is_d1 = (M.d1_pt < M.d2_pt).values
s_eta = np.where(soft_is_d1, M.d1_eta, M.d2_eta)
s_phim = np.where(soft_is_d1, M.d1_phim, M.d2_phim)
d_s = np.sqrt((mu.eta - ak.Array(s_eta)) ** 2 + dphi(mu.phi, ak.Array(s_phim)) ** 2)
M["soft_dau_unreco"] = ~ak.to_numpy(ak.any(d_s < mu_r, axis=1))
on_other = np.zeros(len(M), dtype=bool)
dr_other = np.full(len(M), np.inf)
for k in ("o1", "o2"):
    opt = ak.to_numpy(RM[f"{k}_pt"]).astype(float)
    oeta = ak.to_numpy(RM[f"{k}_eta"]).astype(float)
    ophi = ak.to_numpy(RM[f"{k}_phi"]).astype(float)
    oq = ak.to_numpy(RM[f"{k}_charge"]).astype(float)
    olxy = np.clip(ak.to_numpy(RM.o_a_lxy).astype(float), 0, None)
    ophim = phi_at_pca(ophi, opt, oq, olxy, ak.to_numpy(RM.o_a_phi).astype(float))
    exists = opt > 0
    for i in (0, 1):
        ri = np.where(M[f"m{i}_type"] == 0, R_PF, R_DSA)
        d = np.sqrt((M[f"m{i}_eta"].values - oeta) ** 2
                    + dphi(M[f"m{i}_phi"].values, ophim) ** 2)
        on_other |= exists & (d < ri)
        dr_other = np.where(exists, np.minimum(dr_other, d), dr_other)
M["lead_on_other"] = on_other
# smallest dR between m0, m1 and a daughter of the other dark photon; NaN without one
M["dr_lead_other"] = np.where(np.isfinite(dr_other), dr_other, np.nan)
emit(f"# match radii (analyze.py, table3.csv): PF {R_PF:.4f}, DSA {R_DSA:.4f}")
emit(f"# population: {len(M)} gen-matched mu-type LJs (n_mu >= 2: {(M.n_mu >= 2).all()}), "
     f"samples {M['sample'].nunique()}; per mass: "
     + ", ".join(f"{m}: {n}" for m, n in M.groupby("mass").size().items()))
emit()

ROWS3 = ["[0.05, 0.1)", "[0.1, 0.2)", "[0.2, 0.5)"]
MASSES = sorted(M["mass"].unique())


def mass_sets():
    yield "all", M
    for m in MASSES:
        yield str(m), M[M["mass"] == m]


# ------------------------------------------------------------------ T2: composition
rows = []
for lxy in ("[100, 1000)", "[10, 100)"):
    for mlab, X in mass_sets():
        for b in ROWS3:
            g = X[(X.lxy_c == lxy) & (X.gdr_c == b)]
            n = len(g)
            if n == 0:
                rows.append(dict(lxy=lxy, M=mlab, gen_dR=b, N=0))
                continue
            r = dict(lxy=lxy, M=mlab, gen_dR=b, N=n, n_novtx=int(g.nv.sum()),
                     f_novtx=g.nv.mean())
            for c in "ABCDE":
                r[f"cls_{c}"] = (g.cls == c).mean()
            for c in ("PF-PF", "PF-DSA", "DSA-DSA", ">=3mu"):
                r[f"cont_{c}"] = (g.content == c).mean()
            r["same_charge"] = g.same_charge.mean()
            r["med_soft_reco_pt"] = g.soft_reco_pt.median()
            r["med_soft_gen_pt"] = g.soft_gen_pt.median()
            r["acc_eta2p4_pt5"] = g.acc_pf.mean()
            r["acc_eta2p4_pt10"] = g.acc_dsa.mean()
            r["soft_dau_unreco"] = g.soft_dau_unreco.mean()
            r["lead_on_other"] = g.lead_on_other.mean()
            rows.append(r)
T2 = pd.DataFrame(rows)
T2.to_csv(f"{OUT}/T2_composition.csv", index=False)
emit("# T2: composition (>=3mu share = content 'other' = n_mu >= 3)")
emit(T2.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
emit()

# T2_dr_to_other: how far the other dark photon's daughters are from m0, m1
rows = []
for lxy in ("all", "[100, 1000)", "[10, 100)"):
    for mlab, X in mass_sets():
        if lxy != "all":
            X = X[X.lxy_c == lxy]
        for b in ["all"] + ROWS3:
            g = X if b == "all" else X[X.gdr_c == b]
            v = g.dr_lead_other.dropna().values
            r = dict(lxy=lxy, M=mlab, gen_dR=b, N=len(g), n_with_other=len(v))
            if len(v):
                r.update(min=v.min(), q01=np.quantile(v, 0.01), q05=np.quantile(v, 0.05),
                         q16=np.quantile(v, 0.16), q50=np.median(v))
            rows.append(r)
T2o = pd.DataFrame(rows)
T2o.to_csv(f"{OUT}/T2_dr_to_other.csv", index=False)
emit("# T2_dr_to_other: smallest dR between m0, m1 and a daughter of the other dark photon "
     "(propagated), over the LJs with one (n_with_other)")
emit(T2o[T2o.M == "all"].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
emit()

# ------------------------------------------------------------------ T3, T4: standardisation
VARS = {"class": "cls", "content": "content", "n_mu": "nmu", "charge": "charge",
        "soft_reco_pt": "soft_reco_pt_bin", "soft_gen_pt": "soft_gen_pt_bin",
        "region": "region"}
NARROW, WIDE = "[0.1, 0.2)", "[0.2, 0.5)"
NBOOT, SEED = 1000, 20261114


def bootstrap(X, keys):
    """Bootstrap of f_narrow, f_wide_obs, f_wide_pred and gap_closed: NBOOT replicas resample the
    LJs of NARROW and of WIDE separately (multinomial over cell x outcome), fixed seed."""
    both = X[X.gdr_c.isin([NARROW, WIDE])]
    isw = (both.gdr_c == WIDE).values
    nn, nw = int((~isw).sum()), int(isw.sum())
    if nn == 0 or nw == 0:
        return {}
    cid = both.groupby(keys).ngroup().values
    nc = int(cid.max()) + 1
    code = cid * 2 + both.nv.values
    cn = np.bincount(code[~isw], minlength=2 * nc)
    cw = np.bincount(code[isw], minlength=2 * nc)
    rng = np.random.default_rng(SEED)
    bn = rng.multinomial(nn, cn / nn, size=NBOOT).reshape(NBOOT, nc, 2)
    bw = rng.multinomial(nw, cw / nw, size=NBOOT).reshape(NBOOT, nc, 2)
    n_n, k_n = bn.sum(axis=2), bn[:, :, 1]
    m_w = bw.sum(axis=2)
    fn = k_n.sum(axis=1) / nn
    fo = bw[:, :, 1].sum(axis=1) / nw
    covered = n_n > 0
    f_c = np.where(covered, k_n / np.maximum(n_n, 1), 0.0)
    wcov = np.where(covered, m_w, 0).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        pred = (m_w * f_c).sum(axis=1) / wcov
        closed = (pred - fn) / (fo - fn)
    closed = closed[np.isfinite(closed)]
    return dict(f_narrow_bs_err=fn.std(), f_wide_obs_bs_err=fo.std(),
                f_wide_pred_bs_err=np.nanstd(pred), gap_closed_bs_err=closed.std(),
                gap_closed_bs_q16=np.quantile(closed, 0.16),
                gap_closed_bs_q84=np.quantile(closed, 0.84))


def standardise(X, keys):
    """Observed no-vertex fraction of WIDE, and the one predicted from NARROW's per-cell
    fractions with WIDE's cell shares. Binomial errors; the WIDE shares are held fixed."""
    xn, xw = X[X.gdr_c == NARROW], X[X.gdr_c == WIDE]
    gn = xn.groupby(keys).nv.agg(["sum", "count"])
    gw = xw.groupby(keys).nv.agg(["sum", "count"])
    nw, nn = len(xw), len(xn)
    fo = xw.nv.mean() if nw else np.nan
    fn = xn.nv.mean() if nn else np.nan
    eo = np.sqrt(fo * (1 - fo) / nw) if nw else np.nan
    en = np.sqrt(fn * (1 - fn) / nn) if nn else np.nan
    pred, var, uncovered, min_n, small = 0.0, 0.0, 0, np.inf, 0
    for k, row in gw.iterrows():
        w = row["count"] / nw
        if k not in gn.index or gn.loc[k, "count"] == 0:
            uncovered += int(row["count"])
            continue
        c = gn.loc[k]
        f = c["sum"] / c["count"]
        pred += w * f
        # binomial variance with f regularised to (k + 0.5) / (n + 1), so that a cell with 0 or
        # n no-vertex LJs out of few does not carry zero error
        fr = (c["sum"] + 0.5) / (c["count"] + 1)
        var += w ** 2 * fr * (1 - fr) / c["count"]
        min_n = min(min_n, int(c["count"]))
        if c["count"] < 20:
            small += int(row["count"])
    cov = 1 - uncovered / nw if nw else np.nan
    # with uncovered cells, renormalise the prediction to the covered weight
    pred = pred / cov if cov > 0 else np.nan
    ep = np.sqrt(var) / cov if cov > 0 else np.nan
    closed = (pred - fn) / (fo - fn) if nw and nn and fo != fn else np.nan
    return dict(N_narrow=nn, N_wide=nw, f_narrow=fn, err_narrow=en, f_wide_obs=fo,
                err_obs=eo, f_wide_pred=pred, err_pred=ep, gap=fo - fn, gap_closed=closed,
                wide_share_without_narrow_cell=1 - cov,
                wide_share_narrow_n_lt20=small / nw if nw else np.nan,
                min_narrow_cell_n=(min_n if np.isfinite(min_n) else 0),
                n_cells_wide=len(gw),
                n_cells_narrow_lt50=int(sum(1 for k in gw.index
                                            if k not in gn.index or gn.loc[k, "count"] < 50)),
                n_cells_wide_lt50=int((gw["count"] < 50).sum()), **bootstrap(X, keys))


L100 = "[100, 1000)"
T3S = {}
for lxy, suffix in ((L100, ""), ("[10, 100)", "_lxy10_100")):
    cells, summ = [], []
    for mlab, X in mass_sets():
        X = X[X.lxy_c == lxy]
        for vname, col in VARS.items():
            for b in ROWS3:
                g = X[X.gdr_c == b].groupby(col).nv.agg(["sum", "count"])
                for k, row in g.iterrows():
                    cells.append(dict(M=mlab, variable=vname, category=k, gen_dR=b,
                                      n=int(row["count"]), n_novtx=int(row["sum"]),
                                      f_novtx=row["sum"] / row["count"]))
            summ.append(dict(M=mlab, variable=vname, **standardise(X, [col])))
    T3c, T3s = pd.DataFrame(cells), pd.DataFrame(summ)
    T3c.to_csv(f"{OUT}/T3_cells{suffix}.csv", index=False)
    T3s.to_csv(f"{OUT}/T3_summary{suffix}.csv", index=False)
    emit(f"# T3 cells (Lxy {lxy}, M summed): no-vertex fraction (n_novtx / n)")
    piv = T3c[T3c.M == "all"].copy()
    piv["v"] = piv.apply(lambda r: f"{r.f_novtx:.3f} ({r.n_novtx}/{r.n})", axis=1)
    emit(piv.pivot_table(index=["variable", "category"], columns="gen_dR", values="v",
                         aggfunc="first").to_string())
    emit()
    emit(f"# T3 summary (Lxy {lxy}): direct standardisation of the [0.2, 0.5) bin")
    emit(T3s.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    emit()
    T3S[lxy] = T3s

# T4 (Lxy [100, 1000)): the two variables closing most of the gap (M summed); T4_alt: class x
# charge, a fixed choice
allv = T3S[L100][T3S[L100].M == "all"]
CHOICES = {"T4": allv.sort_values("gap_closed", ascending=False).variable.tolist()[:2],
           "T4_alt": ["class", "charge"]}
for tag, best in CHOICES.items():
    emit(f"# {tag}: joint standardisation with {best} (M summed and per M; Lxy [100, 1000))")
    cells, summ = [], []
    for mlab, X in mass_sets():
        X = X[X.lxy_c == L100]
        cols = [VARS[v] for v in best]
        for b in ROWS3:
            g = X[X.gdr_c == b].groupby(cols).nv.agg(["sum", "count"])
            for k, row in g.iterrows():
                cells.append(dict(M=mlab, cell=" | ".join(map(str, k)), gen_dR=b,
                                  n=int(row["count"]), n_novtx=int(row["sum"]),
                                  f_novtx=row["sum"] / row["count"]))
        summ.append(dict(M=mlab, variables=" x ".join(best), **standardise(X, cols)))
    T4c, T4s = pd.DataFrame(cells), pd.DataFrame(summ)
    T4c.to_csv(f"{OUT}/{tag}_cells.csv", index=False)
    T4s.to_csv(f"{OUT}/{tag}_summary.csv", index=False)
    piv = T4c[T4c.M == "all"].copy()
    piv["v"] = piv.apply(lambda r: f"{r.f_novtx:.3f} ({r.n_novtx}/{r.n})", axis=1)
    emit(piv.pivot_table(index="cell", columns="gen_dR", values="v",
                         aggfunc="first").to_string())
    emit(T4s.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    emit()

# T5: the [0.2, 0.5) bin standardised by soft_dau_unreco alone (softer daughter of the matched
# dark photon without an LJ muon within its radius), both Lxy ranges
cells, summ = [], []
for lxy in (L100, "[10, 100)"):
    for mlab, X in mass_sets():
        X = X[X.lxy_c == lxy]
        for b in ROWS3:
            g = X[X.gdr_c == b].groupby("soft_dau_unreco").nv.agg(["sum", "count"])
            for k, row in g.iterrows():
                cells.append(dict(Lxy=lxy, M=mlab, soft_dau_unreco=bool(k), gen_dR=b,
                                  n=int(row["count"]), n_novtx=int(row["sum"]),
                                  f_novtx=row["sum"] / row["count"]))
        summ.append(dict(Lxy=lxy, M=mlab, variable="soft_dau_unreco",
                         **standardise(X, ["soft_dau_unreco"])))
T5c, T5s = pd.DataFrame(cells), pd.DataFrame(summ)
T5c.to_csv(f"{OUT}/T5_cells.csv", index=False)
T5s.to_csv(f"{OUT}/T5_summary.csv", index=False)
emit("# T5: [0.2, 0.5) standardised by soft_dau_unreco (M summed and per M)")
emit(T5s.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
emit()

# T5_acceptance_cells: soft_dau_unreco x acceptance (both daughters |eta| < 2.4, softer pT > 10)
cells = []
for lxy in (L100, "[10, 100)"):
    for mlab, X in mass_sets():
        X = X[X.lxy_c == lxy]
        for b in ROWS3:
            xb = X[X.gdr_c == b]
            for (u, a), g in xb.groupby(["soft_dau_unreco", "acc_dsa"]):
                cells.append(dict(Lxy=lxy, M=mlab, gen_dR=b, soft_dau_unreco=bool(u),
                                  in_acceptance=bool(a), n=len(g), share_of_bin=len(g) / len(xb),
                                  n_novtx=int(g.nv.sum()), f_novtx=g.nv.mean()))
T5a = pd.DataFrame(cells)
T5a.to_csv(f"{OUT}/T5_acceptance_cells.csv", index=False)
emit("# T5_acceptance_cells: soft_dau_unreco x acceptance (both daughters |eta| < 2.4, softer "
     "gen pT > 10 GeV), M summed")
emit(T5a[T5a.M == "all"].to_string(index=False, float_format=lambda x: f"{x:.4f}"))
emit()

with open(f"{OUT}/tables.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
