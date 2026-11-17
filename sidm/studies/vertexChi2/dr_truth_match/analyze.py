"""Truth match of the two leading muons of signal mu-type LJs: the tables.

usage: python analyze.py PARQUET_DIR OUT_DIR [--propagate]

Reads PARQUET_DIR/*.parquet, the per-LJ records written by extract.py, and the by-sample vertex
outcomes of ../vtx_outcome_lxy_dr/outcome_by_sample.csv, which the records must reproduce.
Writes to OUT_DIR: lj_truth.parquet (the merged records), lj_classified.parquet (the gen-matched
LJs with every derived column), check_vs_outcome_by_sample.csv, table1.csv ... table14.csv, and
tables.md (all tables with their captions). The committed tables were made with --propagate.

--propagate: a reco muon is compared with the direction a generated daughter muon has at the
point of closest approach of its helix to the beamline (uniform 3.8 T field along z, starting
at the Zd decay vertex), which is where a reco::Track quotes its momentum. Without it, with the
daughter's direction at the decay vertex. The generated dR is the opening angle at the decay in
both cases.

Terms used in the tables:
- m0, m1, m2: the LJ's muons (PF and DSA pooled) in decreasing pT; m0 and m1 are the leading
  pair. reco dR = dR(m0, m1), the quantity of mu_lj_mumu_dr.
- d1, d2: the two generated Zd daughter muons of the event (the samples are 2Mu2E, one Zd ->
  mumu per event); gen dR = dR(d1, d2) at the decay vertex. Lxy: transverse decay length (cm)
  of the generated dark photon matched to the LJ (gen-matched: a genAs_toMu within dR < 0.4).
- outcome: no_vertex (vtx_chi2 < 0, no stored vertex), chi2_ge5 (stored vertex with
  normChi2 >= 5), pass (0 <= vtx_chi2 < 5).
- match radius: the 99 % quantile of dR(reco muon, nearest daughter) for PF and for DSA muons
  in gen-matched LJs with exactly two opposite-charge muons and a passing vertex, all masses
  pooled (table 3). A muon is "within its radius" of a daughter when their dR is below it.
- classes of a gen-matched LJ, from its leading pair:
  D: m0 or m1 is outside its radius of both daughters.
  A: m0 and m1 are within their radii of different daughters, have opposite charge, and the
     reco dR is consistent with the gen dR (ratio within [1/2, 2], or the two differ by at
     most the sum of the two radii): the decay pair.
  B: as A, but the reco dR is wider than the gen dR beyond that tolerance (empty by
     construction with daughter directions at the decay vertex; with --propagate a few LJs
     fall here, see the caption of table 4).
  C: not A or B, both nearest to the same daughter, and same charge or a PF-DSA pair that
     shares muon-system segments: two reconstructions of one daughter.
  E: the rest (sub-classes E1 to E4 in table 4).
"""
import argparse
import glob
import pathlib

import awkward as ak
import numpy as np
import pandas as pd

AP = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
AP.add_argument("parquet_dir", help="directory holding the per-LJ parquet files of extract.py")
AP.add_argument("out_dir", help="output directory (must exist)")
AP.add_argument("--propagate", action="store_true",
                help="match reco muons to the daughters' directions at the beamline PCA")
ARGS = AP.parse_args()
IN, OUT, PROPAGATE = ARGS.parquet_dir, ARGS.out_dir, ARGS.propagate
EDGES = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, np.inf]
BLAB = ["[0,0.005)", "[0.005,0.01)", "[0.01,0.02)", "[0.02,0.05)", "[0.05,0.1)", "[0.1,0.2)",
        "[0.2,0.3)", "[0.3,0.5)", ">=0.5"]
OUTCOME_BY_SAMPLE = (pathlib.Path(__file__).resolve().parents[1] / "vtx_outcome_lxy_dr"
                     / "outcome_by_sample.csv")

files = sorted(glob.glob(f"{IN}/*.parquet"))
R = ak.concatenate([ak.from_parquet(f) for f in files])   # one record per mu-type LJ
ak.to_parquet(R, f"{OUT}/lj_truth.parquet")
md = []


def emit(*a):
    """Append one line to tables.md."""
    md.append(" ".join(str(x) for x in a))


def dphi(a, b):
    return (a - b + np.pi) % (2 * np.pi) - np.pi


def dr(e1, p1, e2, p2):
    return np.sqrt((e1 - e2) ** 2 + dphi(p1, p2) ** 2)


def table(df, floatfmt="{:.3f}", fmts=None):
    """Markdown table; floats use fmts[column] if given, else floatfmt."""
    fmts = fmts or {}
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(fmts.get(c, floatfmt).format(v)
                                        if isinstance(v, (float, np.floating)) else str(v)
                                        for c, v in zip(cols, r.values)) + " |")
    return "\n".join(lines)


# ------------------------------------------------------------------ per-LJ quantities
mu = R.muons
n_mu = ak.to_numpy(ak.num(mu))
pad = ak.pad_none(mu, 3, axis=1)


def col(f, i, fill):
    return ak.to_numpy(ak.fill_none(pad[f][:, i], fill)).astype(float)


D = pd.DataFrame({   # all mu-type LJs
    "sample": np.array(ak.to_list(R.sample)),
    "run": ak.to_numpy(R.run), "lumi": ak.to_numpy(R.lumi), "event": ak.to_numpy(R.event),
    "vtx_chi2": ak.to_numpy(R.vtx_chi2), "n_mu": n_mu,
    "n_pf": ak.to_numpy(R.n_pf), "n_dsa": ak.to_numpy(R.n_dsa),
    "matched": ak.to_numpy(R.gen_matched), "lxy": ak.to_numpy(R.gen_a_lxy),
    "dr_ref": ak.to_numpy(R.mumu_dr_ref),
})
D["mass"] = D["sample"].str.split("_").str[1].str.replace("GeV", "").astype(int)
for i in range(3):
    for f in ("type", "pt", "eta", "phi", "charge"):
        D[f"m{i}_{f}"] = col(f, i, np.nan)
for k in ("d1", "d2"):
    for f in ("pt", "eta", "phi", "charge"):
        D[f"{k}_{f}"] = ak.to_numpy(R[f"{k}_{f}"]).astype(float)
D["gen_dr"] = dr(D.d1_eta, D.d1_phi, D.d2_eta, D.d2_phi)


def phi_at_pca(phi0, pt, q, lxy_cm, phi_v):
    """Azimuth of the momentum at the beamline PCA of a helix (B = 3.8 T along +z) through the
    transverse vertex lxy*(cos phi_v, sin phi_v) with azimuth phi0 there."""
    rad = pt / (0.2998 * 3.8)                     # m
    vx, vy = lxy_cm / 100 * np.cos(phi_v), lxy_cm / 100 * np.sin(phi_v)
    cx, cy = vx + q * rad * np.sin(phi0), vy - q * rad * np.cos(phi0)
    cn = np.hypot(cx, cy)
    px, py = cx - rad * cx / cn, cy - rad * cy / cn   # circle point nearest the origin
    rx, ry = px - cx, py - cy
    return np.arctan2(-q * rx, q * ry)                # velocity = q * (ry, -rx)


gen_a_phi = ak.to_numpy(R.gen_a_phi).astype(float)
for k in ("d1", "d2"):
    D[f"{k}_phim"] = (phi_at_pca(D[f"{k}_phi"].values, D[f"{k}_pt"].values,
                                 D[f"{k}_charge"].values, np.clip(D.lxy.values, 0, None),
                                 gen_a_phi)
                      if PROPAGATE else D[f"{k}_phi"])
D["reco_dr"] = dr(D.m0_eta, D.m0_phi, D.m1_eta, D.m1_phi)
for i in range(3):
    for k in ("d1", "d2"):
        D[f"dr_m{i}_{k}"] = dr(D[f"m{i}_eta"], D[f"m{i}_phi"], D[f"{k}_eta"], D[f"{k}_phim"])
    D[f"dr_m{i}_near"] = np.minimum(D[f"dr_m{i}_d1"], D[f"dr_m{i}_d2"])
    D[f"m{i}_nearest"] = np.where(D[f"dr_m{i}_d1"] <= D[f"dr_m{i}_d2"], 1, 2)
D["outcome"] = np.where(D.vtx_chi2 < 0, "no_vertex",
                        np.where(D.vtx_chi2 >= 5, "chi2_ge5", "pass"))
D["bin"] = pd.cut(D.reco_dr, EDGES, right=False, labels=BLAB)
D["gen_bin"] = pd.cut(D.gen_dr, EDGES, right=False, labels=BLAB)
lead_type = D.m0_type + D.m1_type
D["pair_type"] = np.select([lead_type == 0, lead_type == 1, lead_type == 2],
                           ["PF+PF", "PF+DSA", "DSA+DSA"], "?")
D["same_charge"] = D.m0_charge == D.m1_charge

# PF+DSA leading pair flagged as overlapping: the DSA muon lists the PF muon among its
# muonMatch1..5idx with muonMatchN >= 1 shared segments
mt = ak.to_numpy(ak.fill_none(pad["type"][:, :2], -1))
midx = ak.to_numpy(ak.fill_none(pad["idx"][:, :2], -9))
shared = np.zeros(len(D), dtype=int)
for j in range(1, 6):
    m_n = ak.to_numpy(ak.fill_none(pad[f"muonMatch{j}"][:, :2], -1))
    m_i = ak.to_numpy(ak.fill_none(pad[f"muonMatch{j}idx"][:, :2], -9))
    # DSA in slot 0 matched to the PF in slot 1, or the reverse
    hit01 = (mt[:, 0] == 1) & (mt[:, 1] == 0) & (m_i[:, 0] == midx[:, 1])
    hit10 = (mt[:, 1] == 1) & (mt[:, 0] == 0) & (m_i[:, 1] == midx[:, 0])
    shared = np.where(hit01, np.maximum(shared, m_n[:, 0]), shared)
    shared = np.where(hit10, np.maximum(shared, m_n[:, 1]), shared)
D["pfdsa_shared_seg"] = shared
D["pfdsa_overlap"] = (D.pair_type == "PF+DSA") & (shared >= 1)

emit("# dR-truth tables" + (" (matching to daughters propagated to the beamline PCA)"
                            if PROPAGATE else " (matching to daughter directions at the decay "
                            "vertex)") + "\n")
emit(f"Input: {len(files)} parquet files, {len(D)} mu-type LJ rows (all with >= 2 muons: "
     f"{(D.n_mu >= 2).all()}); gen-matched {D.matched.sum()}.\n")

# ------------------------------------------------------------------ check against the histogram run
emit("## Check against the histogram run (outcome_by_sample.csv)\n")
emit(f"reco dR from the two leading muons vs mumu_dr_ref (mu_lj_mumu_dr), all rows: "
     f"max |diff| = {np.max(np.abs(D.reco_dr - D.dr_ref)):.3g}\n")
orc = pd.read_csv(OUTCOME_BY_SAMPLE)
orc = orc[orc.selection == "base_ljDisplacementIso_cosAlpha"].set_index("sample")
rows = []
for s, g in D.groupby("sample"):
    f = g.outcome.value_counts(normalize=True)
    o = orc.loc[s]
    r = dict(sample=s, n=len(g), n_ref=int(o.raw_or_weighted_mu_lj),
             no_vtx=f.get("no_vertex", 0), no_vtx_ref=o.ge2_frac_no_vertex,
             chi2_ge5=f.get("chi2_ge5", 0), chi2_ge5_ref=o.ge2_frac_chi2_ge5,
             pass_=f.get("pass", 0), pass_ref=o.ge2_frac_pass)
    r["max_abs_diff"] = max(abs(r["no_vtx"] - r["no_vtx_ref"]),
                            abs(r["chi2_ge5"] - r["chi2_ge5_ref"]),
                            abs(r["pass_"] - r["pass_ref"]))
    rows.append(r)
O = pd.DataFrame(rows)
O.to_csv(f"{OUT}/check_vs_outcome_by_sample.csv", index=False)
emit(table(O))
emit(f"\nlargest |diff| over samples: {O.max_abs_diff.max():.4f} "
     f"({'PASS' if O.max_abs_diff.max() <= 0.005 else 'fails'} at 0.005)\n")

# population of the study: gen-matched
M = D[D.matched].copy()


def by_bin(df, fn):
    out = []
    for m in (100, 500, 1000):
        for b in BLAB:
            g = df[(df.mass == m) & (df.bin == b)]
            if len(g) == 0:
                continue
            out.append(dict(mass=m, reco_dR=b, n=len(g), **fn(g)))
    return pd.DataFrame(out)


# ------------------------------------------------------------------ table 1
T1 = by_bin(M, lambda g: dict(
    gen_dR_q10=g.gen_dr.quantile(.1), gen_dR_med=g.gen_dr.median(),
    gen_dR_q90=g.gen_dr.quantile(.9),
    frac_gen_lt_0p05=(g.gen_dr < 0.05).mean(), frac_gen_same_bin=(g.gen_bin == g.bin).mean()))
T1.to_csv(f"{OUT}/table1.csv", index=False)
emit("## Table 1: generated daughter dR per reco-dR bin\n")
emit(table(T1, "{:.4f}"))

# ------------------------------------------------------------------ table 2
T2 = by_bin(M, lambda g: dict(
    PFPF=(g.pair_type == "PF+PF").mean(), PFDSA=(g.pair_type == "PF+DSA").mean(),
    DSADSA=(g.pair_type == "DSA+DSA").mean(), n_mu_2=(g.n_mu == 2).mean(),
    n_mu_ge3=(g.n_mu >= 3).mean(), same_charge=g.same_charge.mean()))
T2.to_csv(f"{OUT}/table2.csv", index=False)
emit("\n## Table 2: leading-pair type, muon count, same-charge fraction\n")
emit(table(T2))

# ------------------------------------------------------------------ table 3: resolution
res = M[(M.n_mu == 2) & (~M.same_charge) & (M.outcome == "pass")]
pf = pd.concat([res.loc[res[f"m{i}_type"] == 0, f"dr_m{i}_near"] for i in (0, 1)])
dsa_rows = [(res.loc[res[f"m{i}_type"] == 1, f"dr_m{i}_near"],
             res.loc[res[f"m{i}_type"] == 1, "lxy"]) for i in (0, 1)]
dsa = pd.concat([a for a, _ in dsa_rows])
dsa_lxy = pd.concat([b for _, b in dsa_rows])
qs = [.5, .9, .99]
rows = [dict(muons="PF", n=len(pf), **{f"q{int(q*100)}": pf.quantile(q) for q in qs}),
        dict(muons="DSA", n=len(dsa), **{f"q{int(q*100)}": dsa.quantile(q) for q in qs})]
for lab, sel in (("DSA, Lxy < 10 cm", dsa_lxy < 10),
                 ("DSA, 10 <= Lxy < 100 cm", (dsa_lxy >= 10) & (dsa_lxy < 100)),
                 ("DSA, Lxy >= 100 cm", dsa_lxy >= 100)):
    x = dsa[sel.values]
    rows.append(dict(muons=lab, n=len(x), **{f"q{int(q*100)}": x.quantile(q) for q in qs}))
# per mass, for information
for m in (100, 500, 1000):
    r_m = res[res.mass == m]
    for t, name in ((0, "PF"), (1, "DSA")):
        x = pd.concat([r_m.loc[r_m[f"m{i}_type"] == t, f"dr_m{i}_near"] for i in (0, 1)])
        rows.append(dict(muons=f"{name}, {m} GeV only", n=len(x),
                         **{f"q{int(q*100)}": x.quantile(q) for q in qs}))
T3 = pd.DataFrame(rows)
T3.to_csv(f"{OUT}/table3.csv", index=False)
R_PF, R_DSA = float(pf.quantile(.99)), float(dsa.quantile(.99))
emit("\n## Table 3: dR(reco muon, nearest daughter), LJs with exactly 2 opposite-charge muons "
     "and a passing vertex, all three masses pooled\n")
emit(table(T3, "{:.4f}"))
emit(f"\nMatch radii (99 % quantile, pooled): PF {R_PF:.4f}, DSA {R_DSA:.4f}\n")


# ------------------------------------------------------------------ table 4: classes
def radius(i):
    return np.where(M[f"m{i}_type"] == 0, R_PF, R_DSA)


r0, r1 = radius(0), radius(1)
in0 = {k: M[f"dr_m0_{k}"] < r0 for k in ("d1", "d2")}
in1 = {k: M[f"dr_m1_{k}"] < r1 for k in ("d1", "d2")}
any0 = in0["d1"] | in0["d2"]
any1 = in1["d1"] | in1["d2"]
one2one = (in0["d1"] & in1["d2"]) | (in0["d2"] & in1["d1"])
opp = ~M.same_charge
tol_abs = r0 + r1
ratio_ok = (M.gen_dr >= M.reco_dr / 2) & (M.gen_dr <= 2 * M.reco_dr)
consistent = ratio_ok | (np.abs(M.gen_dr - M.reco_dr) <= tol_abs)
too_wide = (~consistent) & (M.reco_dr > M.gen_dr)
same_near = M.m0_nearest == M.m1_nearest
isD = ~(any0 & any1)
isAB = ~isD & one2one & opp
isA = isAB & consistent
isB = isAB & too_wide
isC = ~isD & ~isAB & same_near & (M.same_charge | M.pfdsa_overlap)
cls = np.select([isD, isA, isB, isC], ["D", "A", "B", "C"], "E")
# sub-classes of E
e_why = np.select(
    [isAB & ~consistent & ~too_wide,
     one2one & ~opp,
     ~one2one & same_near & opp & ~M.pfdsa_overlap],
    ["E1: decay pair, reco dR narrower than gen beyond tolerance",
     "E2: one-to-one but same charge (nearest different daughters)",
     "E3: both nearest same daughter, opposite charge, no PF-DSA segment overlap"],
    "E4: other")
M["cls"] = cls
M["e_why"] = np.where(cls == "E", e_why, "")
# both leading muons within radius of both daughters (the assignment cannot tell them apart)
M["ambiguous"] = in0["d1"] & in0["d2"] & in1["d1"] & in1["d2"]
M["r_lead0"], M["r_lead1"] = r0, r1

rows = []
for m in (100, 500, 1000):
    for b in BLAB:
        g = M[(M.mass == m) & (M.bin == b)]
        n = len(g)
        if n == 0:
            continue
        r = dict(mass=m, reco_dR=b, n=n)
        for c in "ABCDE":
            p = (g.cls == c).mean()
            r[c] = f"{p:.3f}±{np.sqrt(p*(1-p)/n):.3f}"
        r["ambig"] = f"{g.ambiguous.mean():.3f}"
        # supplementary, not a class: decay pair (A or B) whose reco dR exceeds 2x the gen dR
        r["AB_reco_gt_2gen"] = (
            f"{(g.cls.isin(['A', 'B']) & (g.reco_dr > 2 * g.gen_dr)).mean():.3f}")
        rows.append(r)
T4 = pd.DataFrame(rows)
T4.to_csv(f"{OUT}/table4.csv", index=False)
emit("\n## Table 4: class fractions (binomial errors). ambig = fraction of LJs whose two leading "
     "muons are each within their radius of BOTH daughters (the assignment cannot tell them "
     "apart). AB_reco_gt_2gen (supplementary, not a class) = fraction of all LJs in the bin that "
     "are A or B and have reco dR > 2 x gen dR. With daughter directions at the decay vertex B "
     "is empty by construction: both muons within radius of different daughters implies "
     "|reco dR - gen dR| < r0 + r1. With --propagate the match uses the directions at the "
     "beamline while the gen dR is the opening angle at the decay, and a few LJs fall in B "
     "(table 5)\n")
emit(table(T4))
emit("\nClass E composition, where E > 5 %:\n")
rows = []
for m in (100, 500, 1000):
    for b in BLAB:
        g = M[(M.mass == m) & (M.bin == b)]
        if len(g) and (g.cls == "E").mean() > 0.05:
            vc = g[g.cls == "E"].e_why.value_counts()
            rows.append(dict(mass=m, reco_dR=b, n=len(g), E=(g.cls == "E").sum(),
                             composition="; ".join(f"{k.split(':')[0]} {v}"
                                                   for k, v in vc.items())))
emit(table(pd.DataFrame(rows)) if rows else "(none)")
emit("\nE sub-classes: " + "; ".join(sorted(set(e_why))) + "\n")

# ------------------------------------------------------------------ table 5
rows = []
for m in (100, 500, 1000):
    for b in BLAB:
        for c in "ABCDE":
            g = M[(M.mass == m) & (M.bin == b) & (M.cls == c)]
            if len(g) == 0:
                continue
            vc = g.outcome.value_counts()
            rows.append(dict(
                mass=m, reco_dR=b, cls=c, n=len(g),
                no_vertex=f"{vc.get('no_vertex', 0)/len(g):.3f} ({vc.get('no_vertex', 0)})",
                chi2_ge5=f"{vc.get('chi2_ge5', 0)/len(g):.3f} ({vc.get('chi2_ge5', 0)})",
                pass_=f"{vc.get('pass', 0)/len(g):.3f} ({vc.get('pass', 0)})"))
T5 = pd.DataFrame(rows)
T5.to_csv(f"{OUT}/table5.csv", index=False)
emit("\n## Table 5: vertex outcome per class, fraction (count)\n")
emit(table(T5))


# ------------------------------------------------------------------ table 6
def rad_i(df, i):
    return np.where(df[f"m{i}_type"] == 0, R_PF, R_DSA)


# Distinct-muon coverage of the daughters, over all LJ constituents (not only three): 2 if two
# different LJ muons lie within radius of d1 and d2 respectively, 1 if some muon lies within
# radius of some daughter, else 0 (with collinear daughters one muon can be within radius of
# both, so a plain count would say 2 for a single muon).
mu_all = R.muons
r_all = ak.where(mu_all["type"] == 0, R_PF, R_DSA)


def _cov(k):
    d_ = np.sqrt((mu_all.eta - R[f"{k}_eta"]) ** 2
                 + dphi(mu_all.phi, ak.Array(D[f"{k}_phim"].values)) ** 2)
    return d_ < r_all


c1, c2 = _cov("d1"), _cov("d2")
n1, n2 = ak.to_numpy(ak.sum(c1, axis=1)), ak.to_numpy(ak.sum(c2, axis=1))
n12 = ak.to_numpy(ak.sum(c1 | c2, axis=1))
cov2 = (n1 >= 1) & (n2 >= 1) & (n12 >= 2)
D["n_dau_cov"] = np.where(cov2, 2, np.where(n12 >= 1, 1, 0))
M["n_dau_cov"] = D.loc[M.index, "n_dau_cov"]
wide = M[(M.reco_dr >= 0.1) & (M.mass.isin([500, 1000]))].copy()

rows = []
for m in (500, 1000):
    for c in ("C", "D"):
        g = wide[(wide.mass == m) & (wide.cls == c)]
        if len(g) == 0:
            rows.append(dict(mass=m, cls=c, n=0))
            continue
        if c == "C":
            # duplicate = the lower-pT leading muon (slot 1); the other daughter = not its nearest
            j = np.ones(len(g), dtype=int)
            other = np.where(g.m1_nearest == 1, 2, 1)
            extra = dict(frac_same_charge=float(g.same_charge.mean()),
                         frac_pfdsa_overlap=float(g.pfdsa_overlap.mean()))
        else:
            # extra = a leading muon outside its radius of both daughters (slot 1 if both are)
            out0 = ~((g.dr_m0_d1 < rad_i(g, 0)) | (g.dr_m0_d2 < rad_i(g, 0)))
            out1 = ~((g.dr_m1_d1 < rad_i(g, 1)) | (g.dr_m1_d2 < rad_i(g, 1)))
            j = np.where(out1, 1, 0)
            # the other leading muon's nearest daughter is "taken"; the third muon may be the other
            keep = np.where(j == 1, g.m0_nearest, g.m1_nearest)
            other = np.where(keep == 1, 2, 1)
            extra = dict(frac_both_lead_out=float((out0 & out1).mean()))
        typ = np.where(j == 0, g.m0_type, g.m1_type)
        pt = np.where(j == 0, g.m0_pt, g.m1_pt)
        eta = np.abs(np.where(j == 0, g.m0_eta, g.m1_eta))
        has3 = g.n_mu >= 3
        dr3_other = np.where(other == 1, g.dr_m2_d1, g.dr_m2_d2)
        third_is_other = has3 & (dr3_other < rad_i(g, 2))
        drn = np.where(j == 0, g.dr_m0_near, g.dr_m1_near)
        r = dict(mass=m, cls=c, n=len(g),
                 dR_near_q10=np.quantile(drn, .1), dR_near_q50=np.quantile(drn, .5),
                 dR_near_q90=np.quantile(drn, .9),
                 frac_DSA=float((typ == 1).mean()),
                 pt_q10=np.quantile(pt, .1), pt_q50=np.quantile(pt, .5),
                 pt_q90=np.quantile(pt, .9),
                 abseta_q50=np.quantile(eta, .5), frac_abseta_gt1p2=float((eta > 1.2).mean()),
                 frac_n_mu_ge3=float(has3.mean()),
                 frac_third_is_other_dau=float(third_is_other.mean()),
                 n_dau_cov_0=float((g.n_dau_cov == 0).mean()),
                 n_dau_cov_1=float((g.n_dau_cov == 1).mean()),
                 n_dau_cov_2=float((g.n_dau_cov == 2).mean()))
        r.update(extra)
        rows.append(r)
T6 = pd.DataFrame(rows)
T6.to_csv(f"{OUT}/table6.csv", index=False)
emit("\n## Table 6: classes C and D at reco dR >= 0.1, 500 and 1000 GeV. The muon described is "
     "the duplicate (C: the lower-pT leading muon) or the extra one (D: the leading muon outside "
     "its radius of both daughters, the subleading one if both are). 'third is other dau' = the "
     "LJ's third muon lies within its radius of the daughter not claimed by the other leading "
     "muon. n_dau_cov = number of daughters that distinct LJ muons (all constituents) lie within "
     "radius of; dR_near = the described muon's dR to its nearest daughter.\n")
emit(table(T6.T.reset_index().rename(columns={"index": "quantity"}), "{:.3f}"))

M.to_parquet(f"{OUT}/lj_classified.parquet")

# ------------------------------------------------------------------ tables 7 to 14
# Tables 7 to 13 use the gen-matched LJs of the 500 and 1000 GeV samples.
HM = M[M.mass.isin([500, 1000])].copy()
n_dau_matched = ak.to_numpy(R.n_gen_dau)[D.matched.values]
emit(f"\nGen-matched LJs whose event has exactly two generated Zd daughter muons (d1, d2): "
     f"{int((n_dau_matched == 2).sum())} of {len(n_dau_matched)}.")
LXY_FINE = [-1, 1, 10, 20, 30, 40, 50, 60, 70, 100, 300, 1e9]
LXY_COARSE = [-1, 1, 10, 100, 300, 1e9]
GEN_EDGES = [0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 10]
RECO_EDGES = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 10]
PHI_WINDOW = 0.03   # |dphi| window that ties a reco muon to a daughter in tables 8 and 10
TYPE_NAME = {0: "PF", 1: "DSA"}


def lxy_labels(edges):
    """Labels of Lxy slices closed on the right: <=1, (1,10], ..., >300 (cm)."""
    inner = [f"({lo:g},{hi:g}]" for lo, hi in zip(edges[1:-2], edges[2:-1])]
    return [f"<={edges[1]:g}"] + inner + [f">{edges[-2]:g}"]


def dr_labels(edges):
    """Labels of dR bins closed on the left: [0,0.002), ..."""
    return [f"[{lo:g},{hi:g})" for lo, hi in zip(edges[:-1], edges[1:])]


def outcome_fractions(df, keys, outcomes=("no_vertex", "chi2_ge5", "pass"), not_pass=False):
    """n and the fraction of each vertex outcome per group of keys (non-empty groups only)."""
    df = df.assign(**{f"is_{o}": df.outcome == o for o in outcomes},
                   is_not_pass=df.outcome != "pass")
    agg = {"n": ("outcome", "size")}
    agg.update({o: (f"is_{o}", "mean") for o in outcomes})
    if not_pass:
        agg["not_pass"] = ("is_not_pass", "mean")
    return df.groupby(keys, observed=True).agg(**agg).reset_index()


def q3(x, prefix, levels=(10, 50, 90)):
    """The 10/50/90 % quantiles of x as a dict."""
    return {f"{prefix}_q{q}": x.quantile(q / 100) for q in levels}


# ---------------------------------------------- table 7: size of the wide reco-dR population
n_novtx_hm = (HM.outcome == "no_vertex").sum()
rows = []
for lo in (0.1, 0.2):
    g = HM[HM.reco_dr >= lo]
    cf = g.cls.value_counts(normalize=True)
    gd = g[g.cls == "D"]
    out0, out1 = gd.dr_m0_near >= gd.r_lead0, gd.dr_m1_near >= gd.r_lead1
    rows.append(dict(reco_dR=f">={lo:g}", n=len(g), n_all=len(HM),
                     frac_of_all=len(g) / len(HM),
                     **{c: float(cf.get(c, 0.0)) for c in "ABCDE"},
                     no_vertex=(g.outcome == "no_vertex").mean(),
                     share_of_all_no_vertex=(g.outcome == "no_vertex").sum() / n_novtx_hm,
                     n_D=len(gd), D_m0_out_only=int((out0 & ~out1).sum()),
                     D_m1_out_only=int((~out0 & out1).sum()), D_both_out=int((out0 & out1).sum())))
T7 = pd.DataFrame(rows)
T7.to_csv(f"{OUT}/table7.csv", index=False)
emit("\n## Table 7: LJs at wide reco dR, 500 and 1000 GeV pooled. n_all = all gen-matched LJs "
     "of those masses; frac_of_all = n / n_all; A to E = class fractions; no_vertex = no-vertex "
     "fraction of the n LJs; share_of_all_no_vertex = their no-vertex LJs over all no-vertex "
     "LJs of those masses; n_D = class D LJs, split by which leading muon is outside its radius "
     "of both daughters (m0, m1, or both)\n")
emit(table(T7, "{:.4f}"))

# ---------------------------------------------- table 8: the muon outside its radius in class D
# Class D LJs at reco dR >= 0.1 (500 and 1000 GeV) with exactly one leading muon outside its
# radius of both daughters. "claimed" daughter = the nearest daughter of the leading muon that
# is inside; "other" daughter = the remaining one. Daughter azimuths as matched (propagated
# with --propagate).
gd = HM[(HM.cls == "D") & (HM.reco_dr >= 0.1)]
S = gd[(gd.dr_m0_near >= gd.r_lead0) ^ (gd.dr_m1_near >= gd.r_lead1)].copy()
i_out = np.where((S.dr_m0_near >= S.r_lead0).values, 0, 1)
i_in = 1 - i_out
for f in ("type", "pt", "eta", "phi", "charge"):
    S["out_" + f] = np.where(i_out == 0, S["m0_" + f].values, S["m1_" + f].values)
k_in = np.where(i_in == 0, S.m0_nearest.values, S.m1_nearest.values)
for f in ("pt", "eta", "phim", "charge"):
    S["claimed_" + f] = np.where(k_in == 1, S["d1_" + f].values, S["d2_" + f].values)
    S["other_" + f] = np.where(k_in == 1, S["d2_" + f].values, S["d1_" + f].values)
S["phi_other"] = dphi(S.out_phi, S.other_phim).abs() < PHI_WINDOW
S["phi_claimed"] = dphi(S.out_phi, S.claimed_phim).abs() < PHI_WINDOW
S["q_other"] = S.out_charge == S.other_charge
S["pt_other"] = (S.out_pt / S.other_pt - 1).abs() < 0.3
S["deta_other"] = (S.out_eta - S.other_eta).abs()
S["out_type_name"] = np.where(S.out_type == 0, "PF", "DSA")
rows = []
for t, g in S.groupby("out_type_name"):
    a = g.phi_other & g.q_other                                       # the other daughter
    b = a & g.pt_other
    c = g.phi_claimed & ~g.q_other & ~(g.phi_other & g.q_other)      # a second copy of claimed
    ga = g[a]
    rows.append(dict(muon_outside=t, n=len(g),
                     phi_near_other=g.phi_other.mean(), phi_near_claimed=g.phi_claimed.mean(),
                     phi_near_both=(g.phi_other & g.phi_claimed).mean(),
                     phi_near_neither=(~g.phi_other & ~g.phi_claimed).mean(),
                     other_phi_charge=a.mean(), other_phi_charge_pt=b.mean(),
                     claimed_phi_charge=c.mean(), none_of_these=(~a & ~c).mean(),
                     **q3(ga.deta_other, "other_deta"),
                     other_no_vertex=(ga.outcome == "no_vertex").mean(),
                     other_chi2_ge5=(ga.outcome == "chi2_ge5").mean(),
                     other_pass=(ga.outcome == "pass").mean(),
                     **q3(ga.lxy, "other_lxy")))
T8 = pd.DataFrame(rows)
T8.to_csv(f"{OUT}/table8.csv", index=False)
emit("\n## Table 8: the leading muon outside its radius in class D, reco dR >= 0.1, 500 and "
     "1000 GeV, LJs with exactly one such muon, by its type. 'claimed' = the daughter nearest to "
     "the other leading muon; 'other' = the remaining daughter. phi_near_X = |dphi| < 0.03 to "
     "daughter X; other_phi_charge = phi near the other daughter and the same charge as it; "
     "other_phi_charge_pt = also |pT / pT(other) - 1| < 0.3; claimed_phi_charge = phi near the "
     "claimed daughter, a charge different from the other daughter's, and not "
     "other_phi_charge; "
     "none_of_these = neither other_phi_charge nor claimed_phi_charge. Columns other_deta_*, "
     "other_<outcome> and other_lxy_* (cm) describe the other_phi_charge LJs\n")
emit(table(T8.T.reset_index().rename(columns={"index": "quantity"}), "{:.3f}"))


# ---------------------------------------------- tables 9, 10: eta resolution against Lxy
def fine_lxy(df):
    return pd.cut(df.lxy, LXY_FINE, labels=lxy_labels(LXY_FINE))


# 9: class A LJs with exactly two muons; each leading muon against its nearest daughter
ref = HM[(HM.cls == "A") & (HM.n_mu == 2)]
parts = []
for j in (0, 1):
    geta = np.where(ref[f"m{j}_nearest"] == 1, ref.d1_eta, ref.d2_eta)
    parts.append(pd.DataFrame(dict(type=ref[f"m{j}_type"].values,
                                   deta=np.abs(ref[f"m{j}_eta"].values - geta),
                                   lxy=ref.lxy.values)))
X = pd.concat(parts)
X["lxy_cm"] = fine_lxy(X)
T9 = X.groupby(["type", "lxy_cm"], observed=True).deta.agg(
    n="size", deta_q50=lambda x: x.quantile(.5), deta_q90=lambda x: x.quantile(.9),
    deta_q99=lambda x: x.quantile(.99)).reset_index()
T9["type"] = T9["type"].astype(int).map(TYPE_NAME)
T9.to_csv(f"{OUT}/table9.csv", index=False)
emit("\n## Table 9: |eta(reco) - eta(nearest daughter)| of the two muons of class A LJs with "
     "exactly two muons, 500 and 1000 GeV, by muon type and generated Lxy (cm, slices closed on "
     "the right): n and the 50/90/99 % quantiles\n")
emit(table(T9, "{:.4f}"))

# 10: every leading muon (m0, m1) of every class within 0.03 in phi of a daughter; eta compared
# with the daughter nearest in phi
parts = []
for j in (0, 1):
    p1 = dphi(HM[f"m{j}_phi"], HM.d1_phim).abs()
    p2 = dphi(HM[f"m{j}_phi"], HM.d2_phim).abs()
    geta = np.where(p1 <= p2, HM.d1_eta, HM.d2_eta)
    parts.append(pd.DataFrame(dict(type=HM[f"m{j}_type"].values,
                                   deta=np.abs(HM[f"m{j}_eta"].values - geta),
                                   dph=np.minimum(p1, p2).values, lxy=HM.lxy.values)))
Y = pd.concat(parts)
Y = Y[Y.dph < PHI_WINDOW].copy()
Y["lxy_cm"] = fine_lxy(Y)
T10 = Y.groupby(["type", "lxy_cm"], observed=True).deta.agg(
    n="size", frac_deta_gt_0p05=lambda x: (x > 0.05).mean()).reset_index()
T10["type"] = T10["type"].astype(int).map(TYPE_NAME)
T10.to_csv(f"{OUT}/table10.csv", index=False)
emit("\n## Table 10: leading muons (m0, m1) of all gen-matched LJs, 500 and 1000 GeV, that lie "
     "within |dphi| < 0.03 of a daughter: fraction with |eta(reco) - eta(that daughter)| > 0.05, "
     "by muon type and generated Lxy (cm)\n")
emit(table(T10, "{:.4f}"))

# ---------------------------------------------- tables 11 to 13: outcome against gen and reco dR
HM["gen_dR_bin"] = pd.cut(HM.gen_dr, GEN_EDGES, right=False, labels=dr_labels(GEN_EDGES))
HM["reco_dR_bin"] = pd.cut(HM.reco_dr, RECO_EDGES, right=False,
                          labels=dr_labels(RECO_EDGES))
HM["lxy_cm"] = pd.cut(HM.lxy, LXY_COARSE, labels=lxy_labels(LXY_COARSE))
BIN_NAMES = {"gen_dR_bin": "gen_dR", "reco_dR_bin": "reco_dR"}

T11 = pd.concat([outcome_fractions(HM[HM.mass == m], "gen_dR_bin").assign(mass=m)
                 for m in (500, 1000)]).rename(columns=BIN_NAMES)
T11 = T11[["mass", "gen_dR", "n", "no_vertex", "chi2_ge5", "pass"]]
T11.to_csv(f"{OUT}/table11.csv", index=False)
emit("\n## Table 11: vertex outcome against the generated dR of the daughters, all gen-matched "
     f"LJs, per mass. LJs outside the bins (gen dR >= 10): {int(HM.gen_dR_bin.isna().sum())}\n")
emit(table(T11))

T12 = outcome_fractions(HM, ["lxy_cm", "gen_dR_bin"], outcomes=("no_vertex",),
                        not_pass=True).rename(columns=BIN_NAMES)
T12.to_csv(f"{OUT}/table12.csv", index=False)
emit("\n## Table 12: no-vertex fraction and not-pass fraction (no vertex or chi2 >= 5) against "
     "the generated dR, per generated Lxy slice (cm), 500 and 1000 GeV pooled\n")
emit(table(T12))

T13 = outcome_fractions(HM, ["lxy_cm", "reco_dR_bin"], outcomes=("no_vertex",),
                        not_pass=True).rename(columns=BIN_NAMES)
T13.to_csv(f"{OUT}/table13.csv", index=False)
emit("\n## Table 13: the same against the reco dR of the leading pair, for comparison with "
     "table 12\n")
emit(table(T13))

# ---------------------------------------------- table 14: per-mass match radii
# Radii: the per-mass 99 % quantiles of table 3 ("PF, <m> GeV only", "DSA, <m> GeV only"),
# rounded to 4 decimals as table 3 prints them.
q99 = T3.set_index("muons").q99
rows = []
for m in (100, 500, 1000):
    rpf = round(float(q99[f"PF, {m} GeV only"]), 4)
    rdsa = round(float(q99[f"DSA, {m} GeV only"]), 4)
    g = M[M.mass == m]
    ra = np.where(g.m0_type == 0, rpf, rdsa)
    rb = np.where(g.m1_type == 0, rpf, rdsa)
    one2one_m = (((g.dr_m0_d1 < ra) & (g.dr_m1_d2 < rb))
                 | ((g.dr_m0_d2 < ra) & (g.dr_m1_d1 < rb)))
    anyout = ~((g.dr_m0_near < ra) & (g.dr_m1_near < rb))
    for lo, hi in ((0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5)):
        s = (g.reco_dr >= lo) & (g.reco_dr < hi)
        n = int(s.sum())
        if n:
            rows.append(dict(mass=m, reco_dR=f"[{lo:g},{hi:g})", n=n, r_PF=rpf, r_DSA=rdsa,
                             decay_pair=float((one2one_m & ~g.same_charge)[s].mean()),
                             one_muon_outside=float(anyout[s].mean()),
                             pooled_radii_A=float((g.cls == "A")[s].mean()),
                             pooled_radii_D=float((g.cls == "D")[s].mean()),
                             gen_dR_median=float(g.gen_dr[s].median())))
T14 = pd.DataFrame(rows)
T14.to_csv(f"{OUT}/table14.csv", index=False)
emit("\n## Table 14: control with per-mass match radii (r_PF, r_DSA: the per-mass 99 % "
     "quantiles of table 3). decay_pair = the leading muons within the per-mass radii of "
     "different daughters and of opposite charge; one_muon_outside = at least one leading muon "
     "outside its per-mass radius of both daughters; pooled_radii_A, pooled_radii_D = the class "
     "A and D fractions of table 4 (pooled radii) in the same bins; gen_dR_median = median "
     "generated dR\n")
emit(table(T14, "{:.3f}", fmts={"r_PF": "{:.4f}", "r_DSA": "{:.4f}", "gen_dR_median": "{:.4f}"}))

with open(f"{OUT}/tables.md", "w", encoding="utf-8") as fh:
    fh.write("\n".join(md) + "\n")
print("\n".join(md))
