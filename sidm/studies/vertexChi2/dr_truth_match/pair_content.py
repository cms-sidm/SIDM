"""Vertex outcome of signal mu-type LJs by the muon content of the LJ: tables 15 to 22.

usage: python pair_content.py CLASSIFIED_PARQUET OUT_DIR

Reads CLASSIFIED_PARQUET, the lj_classified.parquet written by analyze.py (one row per
gen-matched mu-type LJ of the 16 2Mu2E samples). Writes to OUT_DIR: table15.csv ... table22.csv
and pair_content.md (all tables with their captions and the input counts), and prints the same
text. LJs are counted without weights; the three mediator masses are pooled unless a table is
per mass. Every fraction is given with its numerator and denominator.

Terms used in the tables (the others as in analyze.py):
- content: the muon types of an LJ with exactly two muons, from its PF and DSA muon counts:
  PF-PF (two PF muons), PF-DSA (one PF and one DSA muon), DSA-DSA (two DSA muons). "other" = the
  LJs with three or more muons; they enter only the shares of table 17.
- charge: "opposite" or "same", the charges of the two leading muons m0 and m1 (for an LJ with
  exactly two muons, its two muons).
- outcome: no_vertex (vtx_chi2 < 0, no stored vertex), chi2_ge5 (stored vertex with
  normChi2 >= 5), pass (0 <= vtx_chi2 < 5). n_<outcome> = the LJs with that outcome, n = all LJs
  of the row; <outcome> = n_<outcome> / n.
- Lxy (lxy_cm): transverse decay length (cm) of the generated dark photon matched to the LJ, in
  slices closed on the left: [lo,hi) and >=lo for the last one.
- gen dR (gen_dR): dR(d1, d2) of the two generated Zd daughter muons at the decay vertex, in bins
  closed on the left.
- region (tables 20, 22): from eta = (eta(d1) + eta(d2)) / 2, barrel |eta| < 0.9, endcap |eta| > 1.2,
  overlap otherwise.
- soft-muon pT: the lower generated pT of d1 and d2 (GeV), in bins [0,5), [5,10), [10,20),
  [20,40), [40,80), >=80 (tables 21, 22).
- mass point: one (mediator mass, m_Zd) pair, its two lifetime samples pooled.
- consistency test (tables 21, 22): the LJs are split into cells of the listed variables. In a
  cell, each mass point holding at least min_lj_per_point LJs of the cell (100 or 50) enters
  with its no-vertex count k_i out of n_i LJs; p = sum k_i / N, N = sum n_i, and the cell
  contributes chi2 = sum_i (k_i - n_i p)^2 / (n_i p (1 - p)) with ndf = (number of such mass
  points) - 1. Cells with fewer than two such mass points, or with p = 0 or p = 1, are left out.
  chi2, ndf, the number of cells (n_cells) and the number of LJs in the test (n_lj, the sum of
  N) are summed over the cells: the test of one no-vertex fraction shared by all mass points
  inside every cell. extra = sqrt(max(0, S1) / S2) with S1 = sum_i [(k_i - n_i p)^2
  - n_i p (1 - p) (1 - n_i / N)] and S2 = sum_i [n_i^2 (1 - n_i / N)], both summed over all mass
  points of all cells: the spread of the no-vertex fraction between mass points beyond the
  binomial expectation.
- cells (tables 21, 22): Lxy slice x gen dR bin. lxy_split = k: each Lxy slice of table 16 cut
  in k parts, geometric for the slices with a lower edge above 0 and linear for [0,1); the open
  slice >= 300 is not cut. gen_dR_edges: "table 18" = the gen dR bins of table 18; "table 18 +
  0.0005, 0.001" = the same with [0,0.002) cut at 0.0005 and 0.001.
"""
import argparse
import math

import numpy as np
import pandas as pd

AP = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
AP.add_argument("classified_parquet", help="lj_classified.parquet written by analyze.py")
AP.add_argument("out_dir", help="output directory (must exist)")
ARGS = AP.parse_args()
OUT = ARGS.out_dir

CONTENTS = ["PF-PF", "PF-DSA", "DSA-DSA"]
LXY_CONTENT = [0, 1, 10, 30, 60, 110, 300, 1e9]          # tables 16, 17, 21, 22 (cm)
LXY_DSA = [0, 30, 60, 110, 150, 200, 250, 300, 350, 400, 500, 1e9]   # tables 19, 20 (cm)
GEN_EDGES = [0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 10]
GEN_EDGE_SETS = {"table 18": GEN_EDGES,                   # tables 21, 22
                 "table 18 + 0.0005, 0.001": [0, 0.0005, 0.001] + GEN_EDGES[1:]}
SOFT_PT_EDGES = [0, 5, 10, 20, 40, 80, 1e9]               # GeV
REGIONS = ["barrel", "overlap", "endcap"]
LXY_SPLITS = (1, 2, 4, 8)          # table 21; table 22 uses 4 and 8
MIN_POINTS = (100, 50)             # LJs a mass point needs inside a cell (table 21)
md = []


def emit(*a):
    """Append one line to pair_content.md."""
    md.append(" ".join(str(x) for x in a))


def table(df, floatfmt="{:.3f}"):
    """Markdown table; floats use floatfmt."""
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(floatfmt.format(v) if isinstance(v, (float, np.floating))
                                        else str(v) for v in r.values) + " |")
    return "\n".join(lines)


def left_labels(edges):
    """Labels of slices closed on the left: [lo,hi), ..., >=lo for the last (open) slice."""
    inner = [f"[{lo:g},{hi:g})" for lo, hi in zip(edges[:-2], edges[1:-1])]
    return inner + [f">={edges[-2]:g}"]


def dr_labels(edges):
    """Labels of dR bins closed on the left: [0,0.002), ..."""
    return [f"[{lo:g},{hi:g})" for lo, hi in zip(edges[:-1], edges[1:])]


def frac_cell(k, n):
    """'fraction (k/n)' for the markdown tables."""
    return f"{k / n:.3f} ({int(k)}/{int(n)})"


def no_vertex_table(df, keys):
    """n, n_no_vertex and no_vertex per group of keys (non-empty groups only)."""
    out = df.groupby(keys, observed=True).agg(n=("no_vtx", "size"),
                                              n_no_vertex=("no_vtx", "sum")).reset_index()
    out["no_vertex"] = out.n_no_vertex / out.n
    return out


def wide(df, rows, col, k="n_no_vertex", n="n"):
    """Markdown view of a long table: one cell 'fraction (k/n)' per (rows, col), '-' if empty."""
    w = df.assign(cell=[frac_cell(a, b) for a, b in zip(df[k], df[n])])
    w = w.set_index(rows + [col]).cell.unstack(col).fillna("-")
    return w.reset_index().rename_axis(columns=None)


def consistency(df, keys, min_point):
    """chi2, ndf, number of cells, number of LJs and extra of the consistency test over cells
    of keys, mass points with at least min_point LJs in the cell (docstring)."""
    tot, ndf, ncell, nlj, s1, s2 = 0.0, 0, 0, 0, 0.0, 0.0
    agg = df.groupby(keys + ["mass", "mzd"], observed=True).no_vtx.agg(["sum", "size"])
    agg = agg[agg["size"] >= min_point].reset_index()
    for _, c in agg.groupby(keys, observed=True):
        pts = list(zip(c["sum"].tolist(), c["size"].tolist()))
        if len(pts) < 2:
            continue
        k = sum(p[0] for p in pts)
        n = sum(p[1] for p in pts)
        p0 = k / n
        if p0 <= 0 or p0 >= 1:
            continue
        tot += sum((a - b * p0) ** 2 / (b * p0 * (1 - p0)) for a, b in pts)
        ndf += len(pts) - 1
        ncell += 1
        nlj += n
        s1 += sum((a - b * p0) ** 2 - b * p0 * (1 - p0) * (1 - b / n) for a, b in pts)
        s2 += sum(b ** 2 * (1 - b / n) for a, b in pts)
    extra = math.sqrt(max(s1, 0.0) / s2) if s2 > 0 else float("nan")
    return tot, ndf, ncell, nlj, extra


def region_of(df):
    """Region of each LJ from the mean eta of d1 and d2 (docstring)."""
    aeta = ((df.d1_eta + df.d2_eta) / 2).abs()
    return pd.Categorical(np.where(aeta < 0.9, "barrel", np.where(aeta > 1.2, "endcap",
                                                                  "overlap")),
                          categories=REGIONS)


def split(edges, k):
    """Each finite slice of edges cut in k parts: geometric for a lower edge > 0, linear for 0;
    the last (open) slice is kept."""
    out = [edges[0]]
    for lo, hi in zip(edges[:-2], edges[1:-1]):
        cut = np.geomspace(lo, hi, k + 1) if lo > 0 else np.linspace(lo, hi, k + 1)
        out.extend(cut[1:-1].tolist() + [hi])
    return out + [edges[-1]]


# ------------------------------------------------------------------ per-LJ quantities
d = pd.read_parquet(ARGS.classified_parquet)
d["mzd"] = (d["sample"].str.extract(r"GeV_(\d+p\d+)GeV")[0].str.replace("p", ".")
            .astype(float))
d["no_vtx"] = (d.vtx_chi2 < 0).astype(int)
d["outcome"] = np.where(d.vtx_chi2 < 0, "no_vertex",
                        np.where(d.vtx_chi2 >= 5, "chi2_ge5", "pass"))
content = np.select([(d.n_pf == 2) & (d.n_dsa == 0), (d.n_pf == 1) & (d.n_dsa == 1),
                     (d.n_pf == 0) & (d.n_dsa == 2)], CONTENTS, "other")
d["content"] = pd.Categorical(content, categories=CONTENTS + ["other"])
d["charge"] = pd.Categorical(np.where(d.same_charge, "same", "opposite"),
                             categories=["opposite", "same"])
d["lxy_cm"] = pd.cut(d.lxy, LXY_CONTENT, right=False, labels=left_labels(LXY_CONTENT))
d["gen_dR"] = pd.cut(d.gen_dr, GEN_EDGES, right=False, labels=dr_labels(GEN_EDGES))
d["soft_pt"] = pd.cut(np.minimum(d.d1_pt, d.d2_pt), SOFT_PT_EDGES, right=False,
                      labels=left_labels(SOFT_PT_EDGES))
P = d[d.content != "other"]   # the LJs with exactly two muons

emit("# Vertex outcome by the muon content of the LJ\n")
emit(f"Input: {len(d)} gen-matched mu-type LJs; content "
     + ", ".join(f"{c} {int((d.content == c).sum())}" for c in CONTENTS + ["other"])
     + f". 'other' is exactly the LJs with three or more muons: "
     f"{bool(((d.content == 'other') == (d.n_mu >= 3)).all())}.\n")

# ------------------------------------------------------------------ table 15: content x charge
T15 = no_vertex_table(P, ["content", "charge"])
T15.to_csv(f"{OUT}/table15.csv", index=False)
emit("## Table 15: no-vertex fraction by content and charge of the two muons, all masses and "
     "Lxy pooled\n")
emit(table(wide(T15, ["content"], "charge")))

# ------------------------------------------------------------------ table 16: content x Lxy
T16 = no_vertex_table(P, ["lxy_cm", "content"])
T16.to_csv(f"{OUT}/table16.csv", index=False)
emit("\n## Table 16: no-vertex fraction by Lxy slice (cm) and content, all masses pooled\n")
emit(table(wide(T16, ["lxy_cm"], "content")))

# ------------------------------------------------------------------ table 17: content shares
T17 = d.groupby(["mass", "lxy_cm", "content"], observed=True).size().rename("n").reset_index()
T17["n_slice"] = T17.groupby(["mass", "lxy_cm"], observed=True).n.transform("sum")
T17["share"] = T17.n / T17.n_slice
T17.to_csv(f"{OUT}/table17.csv", index=False)
emit("\n## Table 17: share of each content among all gen-matched LJs (other included) per "
     "mediator mass (GeV) and Lxy slice (cm); n_slice = all LJs of the mass and slice\n")
emit(table(wide(T17, ["mass", "lxy_cm"], "content", k="n", n="n_slice")))
near = T17[(T17.lxy_cm == "[0,1)") & (T17.content == "PF-DSA")]
emit("\nPF-DSA share at Lxy < 1 cm: "
     + "; ".join(f"{m} GeV {frac_cell(k, n)}" for m, k, n in zip(near.mass, near.n, near.n_slice)))

# ------------------------------------------------------------------ table 18: content x gen dR
T18 = no_vertex_table(P, ["gen_dR", "content"])
T18.to_csv(f"{OUT}/table18.csv", index=False)
emit("\n## Table 18: no-vertex fraction by gen dR and content, all masses and Lxy pooled\n")
emit(table(wide(T18, ["gen_dR"], "content")))

# ------------------------------------------------------------------ tables 19, 20: DSA-DSA
G = d[(d.content == "DSA-DSA") & (d.charge == "opposite")].copy()
G["lxy_cm"] = pd.cut(G.lxy, LXY_DSA, right=False, labels=left_labels(LXY_DSA))
G["region"] = region_of(G)
T19 = G.groupby("lxy_cm", observed=True).agg(
    n=("outcome", "size"), n_no_vertex=("outcome", lambda o: int((o == "no_vertex").sum())),
    n_chi2_ge5=("outcome", lambda o: int((o == "chi2_ge5").sum())),
    n_pass=("outcome", lambda o: int((o == "pass").sum()))).reset_index()
for o in ("no_vertex", "chi2_ge5", "pass"):
    T19[o] = T19[f"n_{o}"] / T19.n
T19 = T19[["lxy_cm", "n", "n_no_vertex", "no_vertex", "n_chi2_ge5", "chi2_ge5", "n_pass",
           "pass"]]
T19.to_csv(f"{OUT}/table19.csv", index=False)
emit("\n## Table 19: DSA-DSA LJs with opposite charge: the three vertex outcomes by Lxy slice "
     "(cm), all masses pooled\n")
v19 = T19[["lxy_cm", "n"]].copy()
for o in ("no_vertex", "chi2_ge5", "pass"):
    v19[o] = [frac_cell(k, n) for k, n in zip(T19[f"n_{o}"], T19.n)]
emit(table(v19))

T20 = no_vertex_table(G, ["lxy_cm", "region"])
T20.to_csv(f"{OUT}/table20.csv", index=False)
emit("\n## Table 20: DSA-DSA LJs with opposite charge: no-vertex fraction by Lxy slice (cm) and "
     "region (mean eta of the two daughters: barrel |eta| < 0.9, overlap, endcap |eta| > 1.2), "
     "all masses pooled\n")
emit(table(wide(T20, ["lxy_cm"], "region")))

# ------------------------------------------------------------------ table 21: consistency test
P = P.copy()
P["region"] = region_of(P)
for k in LXY_SPLITS:
    P[f"lxy_split{k}"] = pd.cut(P.lxy, split(LXY_CONTENT, k), right=False, labels=False)
for j, edges in enumerate(GEN_EDGE_SETS.values()):
    P[f"gen_dr_set{j}"] = pd.cut(P.gen_dr, edges, right=False, labels=False)
KEY_SETS = (("cells", []), ("cells x content", ["content"]),
            ("cells x soft-muon pT", ["soft_pt"]),
            ("cells x content x soft-muon pT", ["content", "soft_pt"]))


def ladder_view(t):
    """Markdown view of tables 21 and 22: chi2 to one decimal, extra to three."""
    v = t.copy()
    v["chi2"] = v.chi2.map("{:.1f}".format)
    v["extra"] = v.extra.map("{:.3f}".format)
    return table(v)


rows = []
for k in LXY_SPLITS:
    for j, set_name in enumerate(GEN_EDGE_SETS):
        for mp in MIN_POINTS:
            for name, extra_keys in KEY_SETS:
                chi2, ndf, ncell, nlj, ex = consistency(
                    P, [f"lxy_split{k}", f"gen_dr_set{j}"] + extra_keys, mp)
                rows.append(dict(lxy_split=k, gen_dR_edges=set_name, min_lj_per_point=mp,
                                 cells_of=name, n_cells=ncell, n_lj=nlj, chi2=chi2, ndf=ndf,
                                 extra=ex))
T21 = pd.DataFrame(rows)
T21.to_csv(f"{OUT}/table21.csv", index=False)
emit("\n## Table 21: consistency test of the no-vertex fraction across mass points (M, m_Zd) "
     "inside cells, LJs of the three contents (definition in the docstring). Cells: Lxy slices "
     "of table 16 each cut in lxy_split parts x gen dR bins (gen_dR_edges); "
     "min_lj_per_point = LJs a mass point needs in the cell; n_lj = LJs in the test; extra = "
     "spread of the no-vertex fraction between mass points beyond binomial\n")
emit(ladder_view(T21))

# ------------------------------------------------------------------ table 22: per content
KEY_SETS_CONTENT = (("cells", []), ("cells x charge", ["charge"]), ("cells x region", ["region"]),
                    ("cells x soft-muon pT", ["soft_pt"]), ("cells x class", ["cls"]))
rows = []
for ct in CONTENTS:
    for k in (4, 8):
        for j, set_name in enumerate(GEN_EDGE_SETS):
            for name, extra_keys in KEY_SETS_CONTENT:
                chi2, ndf, ncell, nlj, ex = consistency(
                    P[P.content == ct], [f"lxy_split{k}", f"gen_dr_set{j}"] + extra_keys, 50)
                rows.append(dict(content=ct, lxy_split=k, gen_dR_edges=set_name,
                                 min_lj_per_point=50, cells_of=name, n_cells=ncell, n_lj=nlj,
                                 chi2=chi2, ndf=ndf, extra=ex))
T22 = pd.DataFrame(rows)
T22.to_csv(f"{OUT}/table22.csv", index=False)
emit("\n## Table 22: the test of table 21 for each content separately, cells Lxy slice (table "
     "16 slices cut in lxy_split parts) x gen dR bin, alone and with the charge, the region, the "
     "soft-muon pT or the class (A to E of analyze.py) as a further cell variable; mass points "
     "with >= 50 LJs in the cell\n")
emit(ladder_view(T22))

with open(f"{OUT}/pair_content.md", "w", encoding="utf-8") as fh:
    fh.write("\n".join(md) + "\n")
print("\n".join(md))
