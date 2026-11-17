"""What the ten MET filters, the beam-halo filter above all, remove in the 1-LJ control region,
and what they cost the signal.

usage: python halo_filter_counts.py DATA_DIR MC_DIR SIGNAL_DIR OUT_DIR

Reads every per-sample .coffea in the three directories (the halo-filter study run: channels of
configs/selections_vertexChi2_halo.yaml with the histogram met_filter_fail_mask). The mask
histogram of a channel without "pass flags" holds, per event, the bit mask of failed filters
(bit i = MET_FILTERS[i], bit 1 the halo filter; 0 = all ten passed) against the highest mu-LJ DSA
muon pT. A channel family is a control region with its two copies (data, background MC):
<family> (with the filters), <family>_noFlags (the mask is read here) and <family>_failHalo.
Signal: 4mu / 2mu2e (mask) and 4mu_flags / 2mu2e_flags.

Writes to OUT_DIR:
  table1_counts.csv         per sample (and pooled group) and family: events without the filters,
                            with all ten, failing any, failing the halo filter in any combination,
                            failing only the halo filter; fractions with Clopper-Pearson 68%
  table2_counts_by_pt.csv   the same per highest mu-LJ DSA muon pT bin
  table3_per_filter.csv     events failing each filter, by name, at all pT and per pT bin
  table4_signal.csv         signal: fraction of selected events the filters remove, halo share
  checks.txt                the mask histogram against the cutflows, one pass/fail line per
                            sample and channel family; the uniformity of the halo-only rate
                            over the data samples
  tables.md                 the pooled rows of the tables above, as markdown

Units. Data: raw counts. Background MC (weighted histograms): weighted yields from the mask
histogram; raw counts only where a cutflow gives them exactly (without filters, with all ten,
failing only the halo filter); fractions from the weighted yields, Clopper-Pearson on N_eff.
Signal (unweighted histograms): raw counts; the weighted yield of a subset is its raw count times
the mean event weight of the sample's cutflow in that channel (exact for the channel totals),
summed over the samples of a pooled group."""
import glob
import os
import pathlib
import sys

import coffea.util
import numpy as np
import pandas as pd
from scipy.special import xlogy  # pylint: disable=no-name-in-module
from scipy.stats import beta

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[4]))   # repository root
from sidm.definitions.cuts import MET_FILTERS, HALO_FILTER  # pylint: disable=wrong-import-position

DATA_DIR, MC_DIR, SIG_DIR, OUT = sys.argv[1:5]
MASK = "met_filter_fail_mask"
FAMILIES = ["cr_1egmOr1mu_lj", "cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu"]
HALO_BIT = MET_FILTERS.index(HALO_FILTER)
MASKS = np.arange(1 << len(MET_FILTERS))
# mu_lj_dsaMuon_max_pt axis [0, 100, 200, 500], flow view
PT_BINS = {"none": [0], "[0, 100)": [1], "[100, 200)": [2], "[200, 500)": [3], ">= 500": [4]}
RTOL = 1e-6


def cp68(k, n):
    """Clopper-Pearson 68.27% interval of k/n."""
    if n <= 0:
        return np.nan, np.nan
    a = (1 - 0.6827) / 2
    lo = beta.ppf(a, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a, k + 1, n - k) if k < n else 1.0
    return lo, hi


def last(cutflow, ch):
    row = list(cutflow[ch].rows.values())[-1]
    return float(row["raw"]), float(row["weighted"])


def load_dir(d, kind):
    out = []
    for f in sorted(glob.glob(f"{d}/*.coffea")):
        n = os.path.basename(f)[:-len(".coffea")]
        o = coffea.util.load(f)
        s = o["out"][n]
        md = s["metadata"]
        out.append(dict(name=n, kind=kind, hists=s["hists"], cutflow=s["cutflow"],
                        unweighted=bool(next(iter(md.get("unweighted_hist", {False})))),
                        n_files=len({w.filename for w in o.get("processed", ())}),   # processed holds chunks
                        exception=o.get("exception")))
    return out


def mean_weight(s, ch):
    """Mean event weight of a channel's selected events, from its cutflow (0 if none selected)."""
    raw, weighted = last(s["cutflow"], ch)
    return weighted / raw if raw else 0.0


def families(s):
    """(family label, mask channel, with-filters channel, failHalo channel or None)."""
    if s["kind"] == "signal":
        ch = "4mu" if s["name"].startswith("4Mu") else "2mu2e"
        return [(ch, ch, ch + "_flags", None)]
    return [(f, f + "_noFlags", f, f + "_failHalo") for f in FAMILIES]


def mask_arrays(s, mch):
    """Values and variances (mask 0..1023, pT flow bins) of one channel's mask histogram."""
    h = s["hists"][MASK][{"channel": mch}]
    return h.values(flow=True)[1:-1], h.variances(flow=True)[1:-1]


def stats(v, w):
    """Counts of one (summed) mask spectrum: v, w per mask value."""
    halo_any = (MASKS >> HALO_BIT) & 1 == 1
    return dict(total=v.sum(), total_w2=w.sum(), mask0=v[0], fail_any=v[1:].sum(),
                halo_any=v[halo_any].sum(), halo_only=v[1 << HALO_BIT],
                **{f"fail_{fl}": v[(MASKS >> i) & 1 == 1].sum() for i, fl in enumerate(MET_FILTERS)})


def rows_for(label, kind, members):
    """Table 1/2/3 rows of one sample or pooled group. members: list of samples."""
    t1, t2, t3 = [], [], []
    for fam, mch, fch, hch in families(members[0]):
        v = sum(mask_arrays(s, mch)[0] for s in members)
        w = sum(mask_arrays(s, mch)[1] for s in members)
        raw_no = sum(last(s["cutflow"], mch)[0] for s in members)
        raw_with = sum(last(s["cutflow"], fch)[0] for s in members)
        # signal: each sample's raw spectrum times its own mean event weight
        vw = sum(mask_arrays(s, mch)[0] * mean_weight(s, mch) for s in members)
        raw_halo_only = sum(last(s["cutflow"], hch)[0] for s in members) if hch else np.nan
        weighted_hist = not members[0]["unweighted"] and kind != "data"
        for pl, idx in [("all", list(range(5)))] + list(PT_BINS.items()):
            st = stats(v[:, idx].sum(axis=1), w[:, idx].sum(axis=1))
            n = st["total"]
            neff = n * n / st["total_w2"] if st["total_w2"] > 0 else 0.0
            r = dict(sample=label, kind=kind, family=fam, pt_bin=pl, n_samples=len(members),
                     n_files=sum(s["n_files"] for s in members))
            if weighted_hist:       # background MC: weighted yields, raw from cutflows
                r.update(raw_no_filters=raw_no if pl == "all" else np.nan,
                         raw_all_ten=raw_with if pl == "all" else np.nan,
                         raw_halo_only=raw_halo_only if pl == "all" else np.nan, neff=neff)
                r.update({f"w_{k}": st[k] for k in ("total", "mask0", "fail_any", "halo_any",
                                                    "halo_only")})
            else:                   # raw-count histograms (data, signal)
                r.update({f"raw_{k}": st[k] for k in ("total", "mask0", "fail_any", "halo_any",
                                                      "halo_only")})
                if kind == "signal":
                    stw = stats(vw[:, idx].sum(axis=1), w[:, idx].sum(axis=1))
                    r.update({f"w_{k}": stw[k] for k in ("total", "mask0", "fail_any",
                                                         "halo_any", "halo_only")})
                neff = n
            for k in ("mask0", "fail_any", "halo_any", "halo_only"):
                f = st[k] / n if n > 0 else np.nan
                lo, hi = cp68(f * neff, neff) if n > 0 else (np.nan, np.nan)
                r[f"frac_{k}"], r[f"frac_{k}_lo68"], r[f"frac_{k}_hi68"] = f, lo, hi
            (t1 if pl == "all" else t2).append(r)
            t3.append(dict(sample=label, kind=kind, family=fam, pt_bin=pl, total=n,
                           **{fl: st[f"fail_{fl}"] for fl in MET_FILTERS}))
    return t1, t2, t3


def check(s):
    """Mask histogram against the cutflows for one sample: one line per channel family."""
    lines = []
    tol = 0.5 if (s["kind"] == "data" or s["unweighted"]) else None
    for fam, mch, fch, hch in families(s):
        v = mask_arrays(s, mch)[0].sum(axis=1)   # summed over the pT axis, flow included
        unit = 0 if (s["kind"] == "data" or s["unweighted"]) else 1   # cutflow raw / weighted
        tests = [("mask==0 vs with-filters cutflow", v[0], last(s["cutflow"], fch)[unit]),
                 ("total vs noFlags cutflow", v.sum(), last(s["cutflow"], mch)[unit])]
        if hch:
            tests.append(("mask==2 vs failHalo cutflow", v[1 << HALO_BIT],
                          last(s["cutflow"], hch)[unit]))
        ok = all(abs(a - b) <= (tol if tol is not None else RTOL * max(abs(b), 1.0))
                 for _, a, b in tests)
        lines.append(f"{'PASS' if ok else 'FAIL'} {s['kind']:6s} {s['name']} {fam}: "
                     + "; ".join(f"{t} {a:.6g} = {b:.6g}" for t, a, b in tests))
    return lines


samples = load_dir(DATA_DIR, "data") + load_dir(MC_DIR, "MC") + load_dir(SIG_DIR, "signal")
os.makedirs(OUT, exist_ok=True)
log = ["samples read (name, files processed):"]
log += [f"  {s['kind']:6s} {s['name']}: {s['n_files']} files"
        + (f"  EXCEPTION {s['exception']!r}" if s["exception"] else "") for s in samples]
for s in samples:
    log += check(s)
n_fail = sum(1 for line in log if line.startswith("FAIL"))
n_chk = sum(1 for line in log if line.startswith(("PASS", "FAIL")))
log.append(f"checks: {n_chk - n_fail} of {n_chk} (sample, channel family) pass")

T1, T2, T3 = [], [], []
groups = {}
for s in samples:
    a, b, c = rows_for(s["name"], s["kind"], [s])
    T1 += a; T2 += b; T3 += c
    if s["kind"] == "signal":
        fs, m = s["name"].split("_")[:2]
        groups.setdefault(("signal", f"pooled {fs}"), []).append(s)
        groups.setdefault(("signal", f"pooled {fs} {m}"), []).append(s)
    else:
        groups.setdefault((s["kind"], f"pooled {s['kind']}"), []).append(s)
for (kind, label), members in groups.items():
    a, b, c = rows_for(label, kind, members)
    T1 += a; T2 += b; T3 += c
    log.append(f"{label}: {len(members)} samples, {sum(s['n_files'] for s in members)} files: "
               + ", ".join(f"{s['name']} ({s['n_files']})" for s in members))
T1, T2, T3 = pd.DataFrame(T1), pd.DataFrame(T2), pd.DataFrame(T3)

# table 4: signal, 4mu / 2mu2e against their _flags variants
t4 = []
for _, r in T1[T1.kind == "signal"].iterrows():
    n, k = r.raw_total, r.raw_mask0
    lo, hi = cp68(n - k, n)
    t4.append(dict(sample=r["sample"], channel=r.family, raw_no_filters=n, raw_all_ten=k,
                   w_no_filters=r.w_total, w_all_ten=r.w_mask0,
                   frac_removed=(n - k) / n if n else np.nan, frac_removed_lo68=lo,
                   frac_removed_hi68=hi, frac_removed_by_halo=r.frac_halo_any,
                   halo_share_of_removed=r.raw_halo_any / (n - k) if n > k else np.nan))
T4 = pd.DataFrame(t4)



def deviance(k, n):
    """Poisson deviance of counts k in samples of size n against their common rate (last axis)."""
    mu = n * (k.sum(axis=-1, keepdims=True) / n.sum())
    return 2 * (xlogy(k, k / mu) - (k - mu)).sum(axis=-1)


# Is the halo-only rate the same in every data sample? Toys: Poisson counts at the pooled rate.
d1 = T1[(T1.kind == "data") & ~T1["sample"].str.startswith("pooled") & (T1.family == FAMILIES[0])
        & (T1.raw_total > 0)]
if len(d1) > 1 and d1.raw_halo_only.sum() > 0:
    k_obs, n_obs = d1.raw_halo_only.to_numpy(float), d1.raw_total.to_numpy(float)
    toys = np.random.default_rng(1).poisson(n_obs * k_obs.sum() / n_obs.sum(), size=(20000, len(n_obs)))
    toys = toys[toys.sum(axis=1) > 0].astype(float)
    log.append(f"halo-only rate over the {len(d1)} data samples with events ({FAMILIES[0]}): {k_obs.sum():.0f} of "
               f"{n_obs.sum():.0f}; Poisson deviance against one common rate {deviance(k_obs, n_obs):.1f} "
               f"for {len(d1) - 1} dof, p-value from {len(toys)} toys "
               f"{(deviance(toys, n_obs) >= deviance(k_obs, n_obs)).mean():.3f}")

for df, name in ((T1, "table1_counts"), (T2, "table2_counts_by_pt"), (T3, "table3_per_filter"),
                 (T4, "table4_signal")):
    df.to_csv(f"{OUT}/{name}.csv", index=False, float_format="%.10g")   # counts above 1e6 stay exact
with open(f"{OUT}/checks.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(log) + "\n")


def md(df, cols):
    def fmt(x):
        if isinstance(x, (float, np.floating)):
            return "" if np.isnan(x) else (f"{x:.0f}" if x == round(x) and abs(x) >= 1 else f"{x:.4g}")
        return str(x)
    d = df[[c for c in cols if c in df.columns]]
    return "\n".join(["| " + " | ".join(d.columns) + " |", "|" + "---|" * len(d.columns)]
                     + ["| " + " | ".join(fmt(x) for x in row) + " |" for row in d.itertuples(index=False)])


pooled = lambda df: df[df["sample"].str.startswith("pooled") | (df.kind == "MC")]
c1 = ["sample", "family", "n_files", "raw_total", "raw_no_filters", "w_total", "neff", "raw_mask0",
      "raw_all_ten", "w_mask0", "raw_fail_any", "w_fail_any", "raw_halo_any", "w_halo_any",
      "raw_halo_only", "w_halo_only", "frac_fail_any", "frac_fail_any_lo68", "frac_fail_any_hi68",
      "frac_halo_only", "frac_halo_only_lo68", "frac_halo_only_hi68"]
with open(f"{OUT}/tables.md", "w", encoding="utf-8") as fh:
    fh.write("# MET-filter counts (halo_filter_counts.py)\n\n" + "\n".join(
        line for line in log if not line.startswith(("PASS", "  "))) + "\n\n")
    fh.write("## Table 1, pooled data, MC samples, pooled signal\n\n" + md(pooled(T1), c1) + "\n\n")
    fh.write("## Table 2, by highest mu-LJ DSA muon pT\n\n"
             + md(pooled(T2), ["sample", "family", "pt_bin"] + c1[3:]) + "\n\n")
    fh.write("## Table 3, events failing each filter\n\n" + md(pooled(T3), list(T3.columns)) + "\n\n")
    fh.write("## Table 4, signal\n\n" + md(T4, list(T4.columns)) + "\n")
print("\n".join(log))
print("wrote", OUT)
