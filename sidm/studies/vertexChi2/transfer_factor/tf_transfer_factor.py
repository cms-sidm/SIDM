#!/usr/bin/env python3
"""Per-slice data-driven feed-up transfer fold for mu-type LJ DSA muons (2018).

METHOD
------
The mu-LJ DSA-muon pT spectrum shows a data/MC excess above 200 GeV.  The leading
hypothesis is FEED-UP: low-true-pT displaced-standalone (DSA) muons mismeasured to
high pT, at a higher rate in data than in MC.  This script tests that hypothesis
quantitatively with an independent, in-event control sample:

  * MATCHED CONTROL.  Raw DSA muons re-selected with the LJ-source ID cuts but WITHOUT
    the PF cross-cleaning, matched to a selected PF muon within dR < 0.2 (helper
    `matched_dsa`).  pT(PF) stands in for truth.  The booked migration matrices
    `dsaPfMatch_ptPf_ptDsa{,_abseta,_nHits}` are P(pT_DSA | pT_PF, slice).

  * TRANSFER.   T_{c,s}(j|i) = M_{c,s}[i,j] / sum_k M_{c,s}[i,k]  -- the probability that
    a muon of true (PF) pT bin i is MEASURED in target bin j, in channel c and slice s.
    eps_{c,s}(i) = sum_{k in S} T(k|i) is the probability it stays inside the source window.

  * FOLD.       Pred_{c,s}(j) = sum_{i in S} T_{c,s}(j|i) * N_true_{c,s}(i), where
    N_true is the LJ DSA muons' OWN low-pT population in the SAME channel and slice,
    corrected for the fact that the LJ spectrum is itself DSA-measured (see CIRCULARITY).

  * CLOSURE.    C = Obs/Pred, R = C - 1, computed separately for data (data T, data L)
    and for MC (MC T, MC L).  Every central number is a ratio formed INSIDE one sample,
    so the absolute data/MC normalisation and the luminosity never enter a central value.
    They enter only the conservative decision statistic's MC-genuine term and the
    secondary raw-yield table.

  * TARGETS.    TAIL = pT >= 200 GeV (including the >500 overflow) and the
    SIDEBAND SB = [100,200) GeV.  Both are tested in three channels:
    Before (no cosmic veto, a negative control), the production spread-veto CR, and the
    vertex-fit-veto VR (an independent veto algorithm).

  * NULL.       Obs = Pred(feed-up of the LJ's own low-pT muons) + genuine high-pT content.
    The genuine term is sized from MC and enters only the conservative statistic R_cons.

The matched population is a SIZING TOOL, not a substitute for the LJ population:
matched muons have a PF partner and are on average better measured, so T_matched <= T_LJ
and Pred is a LOWER BOUND on feed-up; and the two populations are NOT guaranteed disjoint
in code (the LJ cross-cleaning is a segment match, the control is a dR < 0.2 match), so a
dR-matched but segment-unmatched DSA muon sits in BOTH the matrix and the LJ spectrum.

FLOW-HANDLING RULES (every one of these is asserted or printed at run time)
--------------------------------------------------------------------------
Nothing is ever read with `values()`; every read is `values(flow=True)` /
`variances(flow=True)`, and every mask is derived from `axes[...].edges` at run time --
there is not one literal view index in the physics code.

  dsaPfMatch_ptPf     underflow ASSERTED zero (pT >= 0 by construction).
                      Overflow (pT_PF > 500) is PRINTED and EXCLUDED from every source
                      row: it is the ">= 500 GeV true" row and is never a source.
  dsaPfMatch_ptDsa    underflow ASSERTED zero.  Overflow (pT_DSA > 500) is INCLUDED in
                      every ">= 200" target -- it is the most feed-up-like population
                      (471 matched muons in the 2018C CR alone; dropping it under-counts
                      the tail by ~4 %).
  mu_lj_dsaMuon_pt    underflow ASSERTED zero; overflow INCLUDED in every ">= 200" target
                      (both on the fine 5 GeV / 20 GeV axes and on the coarse cubes).
  |eta|, phi          underflow and overflow ASSERTED zero (object cuts guarantee it).
  nHits               underflow ASSERTED zero; the >= 60 overflow is FOLDED INTO the top
                      nHits slice so the definition stays correct if a later chunk
                      populates it (it is empty today).
  dsaPfMatch_response underflow (r < 0.1) enters denominators only; overflow (r > 100) is
                      INCLUDED in every "r > threshold" numerator.  This axis is the only
                      one where a threshold falls inside a bin, so the inherited
                      log-linear within-bin split is used there and ONLY there -- the pT
                      thresholds 100/200/300 are exact edges on every pT axis and are
                      located with an exact-edge lookup, never interpolated.
  channel             a StrCategory with overflow=True: selection is by NAME
                      (`h[{"channel": name}]`), never by index.

ERROR CONVENTIONS (inherited from the 2018C harvest round)
----------------------------------------------------------
  * Data are unit-weighted (asserted element-wise: values == variances on every data
    file), so data errors are Poisson on raw counts.
  * MC variances are the stored sigma^2 (values scaled by w, variances by w^2 in the
    processor's postprocess); never sqrt(N) on MC.  N_eff = (sum w)^2 / sum w^2 is
    printed beside every MC number anything load-bearing rests on.
  * Subset fractions use the weighted-binomial form, never naive quadrature.
  * Uncertainties on Pred come from Gamma(N_eff + 1/2, wbar) toys that redraw every
    matrix cell and every LJ source bin, so the transfer and the input spectrum are
    varied coherently.
  * |1 - C| is never quoted as a precision.

USAGE
-----
  python tf_transfer_factor.py [--subset 2018C] [--only NAME] [--toys N] [--seed N]
                               [--json PATH]

No plotting.  Nothing under sidm/ is read for anything but the sample YAML.
"""

import argparse
import glob
import json
import math
import os
import re
import sys
from collections import OrderedDict

import numpy as np
import yaml
from coffea.util import load

try:
    from scipy.stats import norm as _norm
    def z_from_p(p):
        p = min(max(p, 1e-300), 1.0 - 1e-16)
        return float(_norm.isf(p))
except Exception:                                          # pragma: no cover
    def z_from_p(p):
        p = min(max(p, 1e-300), 1.0 - 1e-16)
        # Acklam's inverse-normal approximation on the upper tail
        q = 1.0 - p
        if q < 0.02425:
            r = math.sqrt(-2 * math.log(q))
            return -(((((-7.784894002430293e-03 * r - 3.223964580411365e-01) * r
                        - 2.400758277161838e+00) * r - 2.549732539343734e+00) * r
                      + 4.374664141464968e+00) * r + 2.938163982698783e+00) / \
                   ((((7.784695709041462e-03 * r + 3.224671290700398e-01) * r
                      + 2.445134137142996e+00) * r + 3.754408661907416e+00) * r + 1.0)
        r = q - 0.5
        if abs(r) <= 0.425:
            s = r * r
            return (((((-3.969683028665376e+01 * s + 2.209460984245205e+02) * s
                       - 2.759285104469687e+02) * s + 1.383577518672690e+02) * s
                     - 3.066479806614716e+01) * s + 2.506628277459239e+00) * r / \
                   (((((-5.447609879822406e+01 * s + 1.615858368580409e+02) * s
                       - 1.556989798598866e+02) * s + 6.680131188771972e+01) * s
                     - 1.328068155288572e+01) * s + 1.0)
        rr = math.sqrt(-2 * math.log(1.0 - q if r > 0 else q))
        v = (((((-7.784894002430293e-03 * rr - 3.223964580411365e-01) * rr
                - 2.400758277161838e+00) * rr - 2.549732539343734e+00) * rr
              + 4.374664141464968e+00) * rr + 2.938163982698783e+00) / \
            ((((7.784695709041462e-03 * rr + 3.224671290700398e-01) * rr
               + 2.445134137142996e+00) * rr + 3.754408661907416e+00) * rr + 1.0)
        return v if r > 0 else -v

def chi2_sf(x, k):
    """Upper-tail probability of a chi2 with k degrees of freedom (scipy if present,
    Wilson-Hilferty otherwise).  Used for the phi-fit goodness of fit."""
    if k <= 0:
        return 1.0
    if x <= 0:
        return 1.0
    try:
        from scipy.stats import chi2 as _chi2
        return float(_chi2.sf(x, k))
    except Exception:
        t = ((x / k) ** (1.0 / 3.0) - (1 - 2.0 / (9 * k))) / math.sqrt(2.0 / (9 * k))
        return float(0.5 * math.erfc(t / math.sqrt(2.0)))
# --------------------------------------------------------------------------- config
# Repository root; relative paths given on the command line are taken from here.
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
DATA_DIR = "RunOutputFiles/vtxResp_2018Full"     # --data-dir
MC_DIR = "RunOutputFiles/vtxRespMC_2018"     # --mc-dir
YAML_DATA = "sidm/configs/ntuples/data_skimmed.yaml"     # --yaml-data

# External CMS 2018 per-era integrated luminosity, /fb.  NOT in the repo: run_periods.yaml
# carries only the year total (59830 pb^-1).  Used ONLY for the provenance header and the
# secondary raw-yield table; no central number depends on it.
ERA_LUMI_FB = {"A": 14.03, "B": 7.07, "C": 6.90, "D": 31.83}
LUMI_MC_FB = 59.830          # what the processor normalises MC to (from run_periods.yaml)
N_DATA_CHUNKS_EXPECTED = 68
ERA_CHUNKS_EXPECTED = {"A": 3, "B": 14, "C": 1, "D": 50}
MC_SAMPLES_EXPECTED = [
    "DYJetsToLL_M50", "QCD_Pt80To120", "QCD_Pt120To170", "QCD_Pt170To300",
    "QCD_Pt300To470", "QCD_Pt470To600", "QCD_Pt600To800", "QCD_Pt800To1000",
]

CH_BEFORE = "cr_1egmOr1mu_lj"
CH_CR = "cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu"
CH_VR = "cr_1egmOr1mu_lj_vtxChi2_CosMu"
CH_COS_SPREAD = "cosmic_muons"
CH_COS_VTX = "cosmic_muons_vtxChi2"
FOLD_CHANNELS = [("Before", CH_BEFORE), ("CR", CH_CR), ("VR", CH_VR)]
COSMIC_FOR = {CH_BEFORE: CH_COS_SPREAD, CH_CR: CH_COS_SPREAD, CH_VR: CH_COS_VTX}

HISTS = [
    "dsaPfMatch_ptPf_ptDsa", "dsaPfMatch_ptPf_ptDsa_abseta", "dsaPfMatch_ptPf_ptDsa_nHits",
    "mu_lj_dsaMuon_pt", "mu_lj_dsaMuon_pt_abseta", "mu_lj_dsaMuon_pt_nHits",
    "mu_lj_dsaMuon_pt_abseta_nHits", "mu_lj_dsaMuon_pt_phi_nHits", "mu_lj_dsaMuon_pt_phi",
    "dsaPfMatch_response_abseta", "dsaPfMatch_response_nHits", "dsaPfMatch_response_relPtErr",
    "mu_lj_dsaMuon_pt_relPtErr",
]

G8 = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 200.0, 300.0, 500.0]   # common merge grid
NATIVE = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0, 200.0, 300.0, 500.0]
ABSETA_BOUNDS = [0.0, 0.9, 1.2, 2.4]
NHITS_BOUNDS = [0.0, 20.0, 30.0, 60.0]

# EM fixed point: cap, tolerance and the number of CONSECUTIVE iterations the tolerance
# must hold for.  See the stopping-rule comment in fold_core.
CONV_ITERS = 300000
CONV_TOL = 1e-9
CONV_STABLE = 5

# --------------------------------------------------------------------------- checks
CHECKS = []


# A check is one of three KINDS, and the distinction matters when reading a FAIL:
#   code        - implementation / bookkeeping.  A FAIL is a BUG and must be fixed.
#   method      - a numerical claim the method specification makes about itself.  A FAIL
#                 means the specification's own claim does not hold on this data.
#   expectation - a physics expectation stated in advance.  A FAIL is a RESULT: the data
#                 refused the expectation.  It is not a defect in the code.
KIND_BY_SUBSTR = [
    ("naive fold (iteration 0) is biased HIGH", "expectation"),
    ("MC self-closure", "expectation"),
    ("control: Before-", "expectation"),
    ("veto-independent", "expectation"),
    ("R_cons sign is STABLE", "expectation"),
    ("within-slice flatness", "expectation"),
    ("kappa_pop", "expectation"),
    ("phi positive control", "expectation"),
    ("era stability", "expectation"),
    ("circularity iteration converged", "method"),
    ("unit test 2", "method"),
    ("inclusive predictions", "method"),
    ("row guards", "method"),
]


def kind_of(name):
    for sub_, k in KIND_BY_SUBSTR:
        if sub_ in name:
            return k
    return "code"


def check(name, ok, detail="", kind=None):
    kind = kind or kind_of(name)
    CHECKS.append({"name": name, "pass": bool(ok), "detail": str(detail), "kind": kind})
    tag = {"code": "", "method": " method", "expectation": " expectation"}[kind]
    print(f"  [{'PASS' if ok else 'FAIL'}{tag}] {name}" + (f"  ::  {detail}" if detail else ""))
    return bool(ok)


def hr(t):
    print("\n" + "=" * 100)
    print(t)
    print("=" * 100)


def sub(t):
    print("\n--- " + t + " " + "-" * max(0, 95 - len(t)))


# --------------------------------------------------------------------------- axis maps
def edge_index(edges, x, tol=1e-9):
    """Index of the edge equal to x, with a tolerance (bin edges carry float noise:
    the |eta| axis stores 0.9 as 0.8999999999999999)."""
    j = int(np.argmin(np.abs(np.asarray(edges) - x)))
    if abs(edges[j] - x) > tol:
        raise AssertionError(f"edge {x} not present in axis (nearest {edges[j]})")
    return j


def merge_map(edges, targets):
    """Groups of FLOW-VIEW indices merging `edges` onto `targets`.

    Flow-view indexing: 0 = underflow, 1..n = inner bins, n+1 = overflow.
    The returned list has len(targets)-1 inner groups plus a final one-element group
    holding the overflow slot.  The underflow slot is never included (it is asserted
    zero wherever this is used).
    """
    edges = np.asarray(edges, dtype=float)
    idx = [edge_index(edges, t) for t in targets]
    assert idx == sorted(idx) and len(set(idx)) == len(idx), "targets not strictly ascending"
    groups = [list(range(a + 1, b + 1)) for a, b in zip(idx[:-1], idx[1:])]
    groups.append([len(edges)])                      # overflow slot
    return groups


def slice_map(edges, bounds, fold_overflow):
    g = merge_map(edges, bounds)
    of = g.pop()
    if fold_overflow:
        g[-1] = g[-1] + of
    return g, of


def merge_axis(a, axis, groups):
    return np.stack([np.take(a, g, axis=axis).sum(axis=axis) for g in groups], axis=axis)


def pick_axis(a, axis, idxs):
    return np.take(a, idxs, axis=axis).sum(axis=axis)


# --------------------------------------------------------------------------- loading
class Store:
    """Accumulated (values, variances) arrays, keyed by histogram name.

    Arrays keep the channel axis WITH its flow slot; channel selection is done by name
    through `self.chan[...]`.  One file is loaded at a time, its needed arrays added in,
    and the loaded object dropped -- nothing but these small arrays is held.
    """

    def __init__(self):
        self.v = {}
        self.w2 = {}
        self.axes = {}          # hn -> list of edge arrays for the non-channel axes
        self.chan = None        # channel name -> index
        self.n_files = 0
        self.samples = []
        self.meta = []

    def add(self, path, is_data, want_hists):
        obj = load(path)
        out = obj["out"]
        base = os.path.basename(path)[:-len(".coffea")]
        keys = list(out.keys())
        assert keys == [base], f"{base}: top-level keys {keys} != [{base}]"
        rec = out[base]
        hists = rec["hists"]
        assert len(hists) == 47, f"{base}: {len(hists)} hists, expected 47"
        md = rec.get("metadata", {})

        def scal(x):
            if isinstance(x, (set, frozenset)):
                return sorted(x)[0] if len(x) == 1 else sorted(x)
            try:
                return x if isinstance(x, (int, float, str)) else list(x)
            except Exception:
                return x
        meta = {"name": base,
                "n_evts": int(md.get("n_evts", 0) or 0),
                "n_removed_golden_json": int(md.get("n_removed_golden_json", 0) or 0),
                "scaled_sum_weights": float(md.get("scaled_sum_weights", 0.0) or 0.0),
                "year": scal(md.get("year", "")),
                "is_data": scal(md.get("is_data", is_data))}
        self.meta.append(meta)

        chans = list(hists["mu_lj_dsaMuon_pt"].axes[0])
        if self.chan is None:
            self.chan = {c: i for i, c in enumerate(chans)}
        else:
            assert list(self.chan.keys()) == chans, f"{base}: channel ordering differs"

        eq = True
        for hn in want_hists:
            h = hists[hn]
            v = np.asarray(h.values(flow=True), dtype=float)
            w2 = h.variances(flow=True)
            w2 = v.copy() if w2 is None else np.asarray(w2, dtype=float)
            ax = [np.asarray(a.edges, dtype=float) for a in h.axes[1:]]
            if hn not in self.v:
                self.v[hn] = v.copy()
                self.w2[hn] = w2.copy()
                self.axes[hn] = ax
            else:
                for a, b in zip(self.axes[hn], ax):
                    assert np.allclose(a, b), f"{base}: axis edges differ for {hn}"
                assert self.v[hn].shape == v.shape, f"{base}: shape differs for {hn}"
                self.v[hn] += v
                self.w2[hn] += w2
            if is_data and not np.array_equal(v, w2):
                eq = False
        self.n_files += 1
        self.samples.append(base)
        del obj, out, hists
        return meta, eq

    def merge(self, other):
        """Add another Store's accumulated arrays into this one (used to build the
        summed MC store from the per-sample stores without loading each file twice)."""
        if self.chan is None:
            self.chan = dict(other.chan)
        else:
            assert list(self.chan.keys()) == list(other.chan.keys()), "channel ordering differs"
        for hn in other.v:
            if hn not in self.v:
                self.v[hn] = other.v[hn].copy()
                self.w2[hn] = other.w2[hn].copy()
                self.axes[hn] = other.axes[hn]
            else:
                for a, b in zip(self.axes[hn], other.axes[hn]):
                    assert np.allclose(a, b), f"axis edges differ for {hn}"
                self.v[hn] += other.v[hn]
                self.w2[hn] += other.w2[hn]
        self.n_files += other.n_files
        self.samples += list(other.samples)
        self.meta += list(other.meta)
        return self
    def get(self, hn, ch):
        ci = self.chan[ch]
        return self.v[hn][ci], self.w2[hn][ci]


def read_sidecar(path):
    try:
        return yaml.safe_load(open(path))
    except Exception:
        return {}


def era_of(name):
    m = re.match(r"DoubleMuon_2018([A-D])", name)
    return m.group(1) if m else "?"


# --------------------------------------------------------------------------- fold core
def fold_core(M, L, S_idx, SB_idx, TAIL_idx, niter=3, conv_iters=0,
              conv_tol=CONV_TOL, conv_stable=CONV_STABLE):
    """Transfer fold with the bounded Bayesian circularity correction.

    M : (..., nb, nb)  merged migration matrix, [true/PF bin i, measured/DSA bin j]
    L : (..., nb)      LJ DSA-muon spectrum on the same merged grid (measured pT)
    Leading dimensions are broadcast (used for the toy ensemble).

    Iteration 0 is the NAIVE fold (N_true = L); the prior of step n is
    pi^(n) = normalise(N_true^(n)) and the transfer used in Bayes' rule is T (not the
    window-renormalised T~), which is the pairing whose fixed point pi ~ N_true
    reproduces the truth exactly.  Only the source block (i, k both < the source edge)
    is corrected; the tail is never unfolded.
    """
    n_row = M.sum(axis=-1)                                   # row totals, ptDsa OF included
    T = M / np.maximum(n_row[..., None], 1e-300)
    T_SS = T[..., S_idx, :][..., :, S_idx]                   # (...,nS,nS)  [i,k]
    eps = T_SS.sum(axis=-1)                                  # P(measured inside window | i)
    LS = L[..., S_idx]
    T_SB = T[..., S_idx, :][..., :, SB_idx].sum(axis=-1)     # (...,nS)
    T_TL = T[..., S_idx, :][..., :, TAIL_idx].sum(axis=-1)

    preds = {}
    Ntrue = LS.astype(float)
    preds[0] = ((T_SB * Ntrue).sum(axis=-1), (T_TL * Ntrue).sum(axis=-1))
    for it in range(1, niter + 1):
        pi = Ntrue / np.maximum(Ntrue.sum(axis=-1, keepdims=True), 1e-300)
        num = T_SS * pi[..., :, None]                        # [i,k]
        den = num.sum(axis=-2, keepdims=True)
        P = num / np.maximum(den, 1e-300)                    # P(i|k)
        Nmeas = (P * LS[..., None, :]).sum(axis=-1)
        Ntrue = Nmeas / np.maximum(eps, 1e-300)
        preds[it] = ((T_SB * Ntrue).sum(axis=-1), (T_TL * Ntrue).sum(axis=-1))
    extra = {}
    if conv_iters:
        # The pre-registered nominal is iteration `niter`; this continues the SAME EM
        # recursion to its fixed point so the distance still to run is measurable.  It is
        # never a central value -- it enters only as an extra term in sigma_circ.
        # STOPPING RULE.  The Pred sequence is NOT monotone: it falls, turns, and rises
        # back to its fixed point.  A test on the DERIVED SCALAR, |Pred_n - Pred_(n-1)|
        # <= tol * |Pred_n|, therefore fires at the TURNING POINT, where the derivative
        # is zero but the recursion is still moving.  At tol = 1e-4 that stopped the
        # deciding slice (CR, TAIL, nHits [0,20)) at iteration 61 on the MINIMUM of the
        # sequence, reporting 647.63 for a fixed point of 677.21.  The rule is now:
        #   (i)  judged on the STATE vector, max_i |N_i^(n) - N_i^(n-1)| / sum_i N_i^(n),
        #        which does not pass through zero at a turning point of a functional;
        #   (ii) the scalar test must hold too, and BOTH must hold for `conv_stable`
        #        CONSECUTIVE iterations, so a single small step cannot end it;
        #   (iii) the iteration cap is large enough to reach the fixed point (the
        #        deciding slice needs ~1.1e5 iterations at tol = 1e-9).
        # If the cap is hit the value returned is the iterate AT the cap and conv_ok is
        # False; it is reported as such rather than as a fixed point.
        Nc = Ntrue.copy()
        last = (T_TL * Nc).sum(axis=-1)
        nit = niter
        stable = 0
        conv_ok = False
        dstate = dpred = float("inf")
        for it in range(niter + 1, niter + conv_iters + 1):
            prev = Nc
            pi = Nc / np.maximum(Nc.sum(axis=-1, keepdims=True), 1e-300)
            num = T_SS * pi[..., :, None]
            den = num.sum(axis=-2, keepdims=True)
            P = num / np.maximum(den, 1e-300)
            Nmeas = (P * LS[..., None, :]).sum(axis=-1)
            Nc = Nmeas / np.maximum(eps, 1e-300)
            cur = (T_TL * Nc).sum(axis=-1)
            nit = it
            scale = np.maximum(np.abs(Nc).sum(axis=-1, keepdims=True), 1e-300)
            dstate = float(np.max(np.abs(Nc - prev) / scale))
            dpred = float(np.max(np.abs(cur - last)
                                 / np.maximum(np.abs(cur), 1e-300)))
            last = cur
            if dstate <= conv_tol and dpred <= conv_tol:
                stable += 1
                if stable >= conv_stable:
                    conv_ok = True
                    break
            else:
                stable = 0
        preds["conv"] = ((T_SB * Nc).sum(axis=-1), (T_TL * Nc).sum(axis=-1))
        # A source bin whose N_true is driven to zero is a DEGENERATE (boundary) fixed
        # point: the limit is well defined and is reported, but it is not interior.
        tot = np.maximum(np.abs(Nc).sum(axis=-1, keepdims=True), 1e-300)
        extra = {"Ntrue_conv": Nc, "n_conv": nit, "conv_ok": conv_ok,
                 "conv_dstate": dstate, "conv_dpred": dpred,
                 "conv_min_frac": float(np.min(Nc / tot))}
    return {"T": T, "eps": eps, "n_row": n_row, "T_SB": T_SB, "T_TL": T_TL,
            "Ntrue": Ntrue, "preds": preds, **extra}


def method_B(M, L, S_idx, SB_idx, TAIL_idx):
    """One-step correction with the MATCHED sample's own pT prior (a different
    population): Q(i|k) = M[i,k] / sum_i' M[i',k] over ALL rows."""
    n_row = M.sum(axis=-1)
    T = M / np.maximum(n_row[..., None], 1e-300)
    col = M.sum(axis=-2)                                     # (..., nb) over all rows
    Q = M / np.maximum(col[..., None, :], 1e-300)            # [i,k]
    eps = T[..., S_idx, :][..., :, S_idx].sum(axis=-1)
    nhat = (Q[..., S_idx, :][..., :, S_idx] * L[..., None, S_idx]).sum(axis=-1)
    nhat = nhat / np.maximum(eps, 1e-300)
    T_SB = T[..., S_idx, :][..., :, SB_idx].sum(axis=-1)
    T_TL = T[..., S_idx, :][..., :, TAIL_idx].sum(axis=-1)
    return (T_SB * nhat).sum(axis=-1), (T_TL * nhat).sum(axis=-1), nhat


def gamma_draw(rng, v, w2, nt, wbar_fallback):
    """Gamma(N_eff + 1/2, wbar) redraw.  Data (w2 == v) gives Gamma(N + 1/2, 1);
    empty cells get the Jeffreys spread instead of a zero-variance delta."""
    v = np.asarray(v, dtype=float)
    w2 = np.asarray(w2, dtype=float)
    pos = v > 0
    neff = np.where(pos, v * v / np.where(w2 > 0, w2, 1.0), 0.0)
    wbar = np.where(pos & (w2 > 0), w2 / np.where(v > 0, v, 1.0), wbar_fallback)
    shape = neff + 0.5
    return rng.gamma(shape, size=(nt,) + v.shape) * wbar


def wbar_of(v, w2):
    s = v.sum()
    return float(w2.sum() / s) if s > 0 else 1.0


def neff_of(v, w2):
    v = np.asarray(v, dtype=float).sum()
    w2 = np.asarray(w2, dtype=float).sum()
    return float(v * v / w2) if w2 > 0 else 0.0


def wbinom_err(num, num_var, den, den_var):
    """Weighted-binomial error on a subset fraction p = num/den (num inside den)."""
    if den <= 0:
        return 0.0
    p = num / den
    val = num_var * (1 - p) ** 2 + max(den_var - num_var, 0.0) * p ** 2
    return math.sqrt(max(val, 0.0)) / den


# --------------------------------------------------------------------------- response
def resp_weights(edges, thr):
    """Flow-view weight vector selecting r > thr on the log-spaced response axis,
    splitting the straddling bin LOG-LINEARLY (the axis has only 4 bins/decade, so
    snapping r > 2 to an edge moves it to 1.884 or 2.239 and shifts the answer by
    several per cent).  Underflow weight 0, overflow weight 1."""
    edges = np.asarray(edges, dtype=float)
    n = len(edges) - 1
    logr = np.log10(edges)
    lt = math.log10(thr)
    w = np.zeros(n + 2)
    w[-1] = 1.0                                       # overflow (r > 100) always above
    j = int(np.searchsorted(logr, lt, side="right")) - 1
    if j < 0:
        w[1:n + 1] = 1.0
    elif j >= n:
        pass
    else:
        w[j + 2:n + 1] = 1.0
        w[j + 1] = (logr[j + 1] - lt) / (logr[j + 1] - logr[j])
    return w


# --------------------------------------------------------------------------- phi test
def phi_partitions(edges):
    c = 0.5 * (np.asarray(edges[:-1]) + np.asarray(edges[1:]))
    d_v = np.minimum(np.abs(c - math.pi / 2), np.abs(c + math.pi / 2))
    d_h = np.minimum(np.abs(c), np.abs(np.abs(c) - math.pi))
    prim_V = np.where(d_v <= 0.5)[0]
    prim_H = np.where(d_h <= 0.5)[0]
    sec_V = np.where((np.abs(c) > math.pi / 4) & (np.abs(c) < 3 * math.pi / 4))[0]
    sec_H = np.where(~((np.abs(c) > math.pi / 4) & (np.abs(c) < 3 * math.pi / 4)))[0]
    return {"primary": (prim_V, prim_H), "secondary": (sec_V, sec_H)}


def bb_lite_nll(nobs, nu, sig2):
    """Poisson NLL profiled over one Barlow-Beeston-lite nuisance per bin."""
    nll = 0.0
    for n_, m_, s2 in zip(nobs, nu, sig2):
        if m_ <= 0:
            nll += 0.0 if n_ <= 0 else 50.0 * n_
            continue
        tau2 = s2 / (m_ * m_) if m_ > 0 else 0.0
        if tau2 <= 0:
            b = 1.0
        else:
            q = m_ * tau2 - 1.0
            b = (-q + math.sqrt(q * q + 4 * n_ * tau2)) / 2.0
            b = max(b, 1e-9)
        mu = b * m_
        nll += mu - (n_ * math.log(mu) if n_ > 0 else 0.0)
        if tau2 > 0:
            nll += (b - 1.0) ** 2 / (2 * tau2)
    return nll


def fit_fcos(nobs, cos_v, cos_w2, nc_v, nc_w2, ngrid=2001):
    """1-parameter binned-Poisson fit of the cosmic fraction with BB-lite template
    uncertainties.  Total normalisation fixed to the observed total."""
    N = float(np.sum(nobs))
    sc, snc = float(cos_v.sum()), float(nc_v.sum())
    if N <= 0 or sc <= 0 or snc <= 0:
        return None
    c = cos_v / sc
    c2 = cos_w2 / sc ** 2
    d = nc_v / snc
    d2 = nc_w2 / snc ** 2
    fs = np.linspace(0.0, 1.0, ngrid)
    nll = np.empty(ngrid)
    for k, f in enumerate(fs):
        nu = N * (f * c + (1 - f) * d)
        sig2 = (N ** 2) * (f * f * c2 + (1 - f) ** 2 * d2)
        nll[k] = bb_lite_nll(nobs, nu, sig2)
    k0 = int(np.argmin(nll))
    dn = nll - nll[k0]
    lo68 = fs[np.where(dn[:k0 + 1] <= 0.5)[0][0]] if np.any(dn[:k0 + 1] <= 0.5) else fs[0]
    hi_rel = np.where(dn[k0:] <= 0.5)[0]
    hi68 = fs[k0 + hi_rel[-1]] if len(hi_rel) else fs[-1]
    lo95 = fs[np.where(dn[:k0 + 1] <= 1.353)[0][0]] if np.any(dn[:k0 + 1] <= 1.353) else fs[0]
    hi_rel = np.where(dn[k0:] <= 1.353)[0]
    hi95 = fs[k0 + hi_rel[-1]] if len(hi_rel) else fs[-1]
    # ---- GOODNESS OF FIT (fix round, findings 2 and 9) -----------------------------
    # The two-template family nu = N*(f*c + (1-f)*d) can only produce shapes INSIDE the
    # span of the two templates.  If the data lie outside it the minimiser rails at a
    # boundary and the profile interval is the width of a REFUTED model, not a bound.
    # Nothing downstream may quote an interval from a fit whose GOF fails, so the GOF is
    # computed here and returned with the fit.  ndf = nbins - 1 (one fitted parameter;
    # the total normalisation is fixed to the observed total, so it costs no dof).
    fbest = float(fs[k0])
    nobs_a = np.asarray(nobs, dtype=float)
    nu_a = np.maximum(N * (fbest * c + (1 - fbest) * d), 1e-12)
    gof_chi2 = float(np.sum((nobs_a - nu_a) ** 2 / nu_a))
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(nobs_a > 0, nobs_a * np.log(np.maximum(nobs_a, 1e-300) / nu_a), 0.0)
    gof_dev = float(2.0 * np.sum(term - (nobs_a - nu_a)))
    gof_ndf = max(len(nobs_a) - 1, 1)
    gof_p = chi2_sf(gof_dev, gof_ndf)
    return {"fhat": fbest, "lo68": float(lo68), "hi68": float(hi68),
            "LB95": float(lo95), "UB95": float(hi95),
            "gof_chi2": gof_chi2, "gof_dev": gof_dev, "gof_ndf": gof_ndf, "gof_p": gof_p,
            "railed_at_boundary": bool(k0 == 0 or k0 == ngrid - 1),
            "rejected": bool(gof_p < 0.01)}


def jnum(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    return str(o)


# =========================================================================== MAIN
def main():
    global DATA_DIR, MC_DIR, YAML_DATA
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", default=None, help="e.g. 2018C: restrict data to DoubleMuon_2018C")
    ap.add_argument("--only", default=None, help="exact data sample name to keep")
    ap.add_argument("--toys", type=int, default=2000)
    ap.add_argument("--scan-toys", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20260901)
    ap.add_argument("--json", default="sidm/studies/vertexChi2/transfer_factor/full/tf_transfer_factor_full68.json")
    ap.add_argument("--data-dir", default=DATA_DIR)
    ap.add_argument("--mc-dir", default=MC_DIR)
    ap.add_argument("--yaml-data", default=YAML_DATA)
    args = ap.parse_args()
    DATA_DIR, MC_DIR, YAML_DATA = (os.path.join(REPO, p) for p in (args.data_dir, args.mc_dir, args.yaml_data))
    args.json = os.path.join(REPO, args.json)
    os.makedirs(os.path.dirname(args.json), exist_ok=True)
    rng = np.random.default_rng(args.seed)
    RES = OrderedDict()
    RES["args"] = {"subset": args.subset, "only": args.only, "toys": args.toys,
                   "scan_toys": args.scan_toys, "seed": args.seed}

    # ------------------------------------------------------------------ T1 provenance
    hr("T1  PROVENANCE")
    data_paths = sorted(glob.glob(f"{DATA_DIR}/*.coffea"))
    mc_paths = sorted(glob.glob(f"{MC_DIR}/*.coffea"))
    if args.only:
        data_paths = [p for p in data_paths if os.path.basename(p)[:-7] == args.only]
    elif args.subset:
        data_paths = [p for p in data_paths if args.subset in os.path.basename(p)]
    print(f"data .coffea globbed: {len(data_paths)}   MC .coffea globbed: {len(mc_paths)}")

    ycfg = yaml.safe_load(open(YAML_DATA))
    ykey = list(ycfg.keys())[0]
    ysamples = ycfg[ykey]["samples"]
    ytext = re.findall(r"^\s+(DoubleMuon_2018[A-D](?:_\d+)?):\s*$", open(YAML_DATA).read(), re.M)
    nfiles_all = {s: len(v["files"]) for s, v in ysamples.items()}
    present = [os.path.basename(p)[:-7] for p in data_paths]
    era_files_tot, era_files_have, era_chunks_have = {}, {}, {}
    for s, nf in nfiles_all.items():
        e = era_of(s)
        era_files_tot[e] = era_files_tot.get(e, 0) + nf
    for s in present:
        e = era_of(s)
        era_files_have[e] = era_files_have.get(e, 0) + nfiles_all.get(s, 0)
        era_chunks_have[e] = era_chunks_have.get(e, 0) + 1
    f_files = sum(era_files_have.values()) / max(sum(era_files_tot.values()), 1)
    L_data = sum(ERA_LUMI_FB[e] * era_files_have.get(e, 0) / max(era_files_tot.get(e, 1), 1)
                 for e in ERA_LUMI_FB)
    mc_present = [os.path.basename(p)[:-7] for p in mc_paths]
    mc_missing = [s for s in MC_SAMPLES_EXPECTED if s not in mc_present]
    print(f"yaml: {len(ytext)} textual sample keys, {len(ysamples)} unique after safe_load "
          f"(duplicate DoubleMuon_2018C -- textual parsers double-count it)")
    print(f"{'era':4s} {'chunks':>12s} {'input files':>18s} {'file frac':>10s} {'lumi /fb':>10s}")
    for e in "ABCD":
        print(f"{e:4s} {era_chunks_have.get(e,0):5d}/{ERA_CHUNKS_EXPECTED[e]:<6d} "
              f"{era_files_have.get(e,0):8d}/{era_files_tot.get(e,0):<9d} "
              f"{era_files_have.get(e,0)/max(era_files_tot.get(e,1),1):10.4f} "
              f"{ERA_LUMI_FB[e]*era_files_have.get(e,0)/max(era_files_tot.get(e,1),1):10.3f}")
    HDR = (f"{len(present)}/{N_DATA_CHUNKS_EXPECTED} chunks "
           f"(A {era_chunks_have.get('A',0)}/3, B {era_chunks_have.get('B',0)}/14, "
           f"C {era_chunks_have.get('C',0)}/1, D {era_chunks_have.get('D',0)}/50); "
           f"f_files={f_files:.4f}; L_data~{L_data:.2f}/fb of {sum(ERA_LUMI_FB.values()):.2f} "
           f"({L_data/sum(ERA_LUMI_FB.values()):.3f}) [external era table, not in repo]; "
           f"MC {len(mc_present)}/8, missing: {mc_missing if mc_missing else 'none'}")
    print("\nPROVENANCE HEADER:\n  " + HDR)
    RES["provenance"] = {"header": HDR, "chunks": len(present), "chunks_expected": 68,
                         "era_chunks": era_chunks_have, "era_files_have": era_files_have,
                         "era_files_tot": era_files_tot, "f_files": f_files,
                         "L_data_fb": L_data, "lumi_frac": L_data / sum(ERA_LUMI_FB.values()),
                         "mc_present": mc_present, "mc_missing": mc_missing,
                         "data_samples": present,
                         "yaml_keys_textual": len(ytext), "yaml_keys_unique": len(ysamples)}

    # ------------------------------------------------------------------ load
    hr("LOADING (one file at a time, arrays accumulated, object discarded)")
    D = Store()
    ERA = {e: Store() for e in "ABCD"}
    empty_chunks, eq_fail = [], []
    for p in data_paths:
        m, eq = D.add(p, True, HISTS)
        ERA[era_of(m["name"])].add(p, True, HISTS)
        if not eq:
            eq_fail.append(m["name"])
        if m["n_evts"] > 0 and m["n_removed_golden_json"] == m["n_evts"]:
            empty_chunks.append((m["name"], m["n_evts"]))
        print(f"  data {m['name']:24s} n_evts={m['n_evts']:>9d} "
              f"golden_removed={m['n_removed_golden_json']:>9d} unit_weights={eq}")
    M = Store()
    MCS = OrderedDict()          # per-sample MC stores: needed for the numerator-cell
                                 # dominant-sample share (fix round, finding 12)
    for p in mc_paths:
        _st = Store()
        m, _ = _st.add(p, False, HISTS)
        MCS[m["name"]] = _st
        M.merge(_st)
        sc = read_sidecar(p.replace(".coffea", ".meta.yaml"))
        xs = None
        for s in (sc.get("samples") or []):
            if isinstance(s, dict) and s.get("name") == m["name"]:
                xs = s.get("xsec_pb")
        v, w2 = M.get("mu_lj_dsaMuon_pt", CH_CR)   # per-file check is on the running sum
        print(f"  MC   {m['name']:24s} n_evts={m['n_evts']:>9d} "
              f"sum_w={m['scaled_sum_weights']:.4e} xs={xs}")
    print(f"\ndata files loaded: {D.n_files}   MC files loaded: {M.n_files}")
    print(f"data samples: {D.samples}")
    print(f"MC samples:   {M.samples}")
    print(f"golden-JSON-empty chunks (n_removed == n_evts), counted as PRESENT: "
          f"{empty_chunks if empty_chunks else 'none'}")
    RES["load"] = {"n_data_files": D.n_files, "n_mc_files": M.n_files,
                   "data_samples": D.samples, "mc_samples": M.samples,
                   "empty_golden_chunks": empty_chunks}

    sc_commits = set()
    for d, ps in ((DATA_DIR, data_paths), (MC_DIR, mc_paths)):
        for p in ps:
            s = read_sidecar(p.replace(".coffea", ".meta.yaml"))
            if s.get("sidm_commit"):
                sc_commits.add(s["sidm_commit"])
    print(f"sidm_commit across sidecars: {sorted(sc_commits)}")

    sub("structural checks")
    check("data unit weights (values == variances element-wise, every data file)",
          not eq_fail, f"failures: {eq_fail}" if eq_fail else f"{D.n_files} files")
    check("one sidm_commit across all sidecars", len(sc_commits) == 1, sorted(sc_commits))
    check("channel ordering identical across all files",
          list(D.chan.keys()) == list(M.chan.keys()), list(D.chan.keys()))
    for hn in HISTS:
        for a, b in zip(D.axes[hn], M.axes[hn]):
            assert np.allclose(a, b), f"data/MC axis mismatch on {hn}"
    check("data and MC axis edges identical on every used histogram", True, f"{len(HISTS)} hists")
    mcw = {}
    for hn in ["mu_lj_dsaMuon_pt"]:
        v, w2 = M.get(hn, CH_CR)
        mcw[hn] = w2.sum() / max(v.sum(), 1e-300)
    print(f"MC sum(var)/sum(val) in CR on mu_lj_dsaMuon_pt = {mcw['mu_lj_dsaMuon_pt']:.6g} "
          f"(a weighted mean of the per-sample lumi*xs weights; per-sample equality to the "
          f"stored weight was verified in the inventory probe)")

    # ------------------------------------------------------------------ grids
    hr("T2  AXIS MAPS AND FLOW AUDIT")
    ax = D.axes
    mm_ptPf = merge_map(ax["dsaPfMatch_ptPf_ptDsa"][0], G8)
    mm_ptDsa = merge_map(ax["dsaPfMatch_ptPf_ptDsa"][1], G8)
    mm_LJ20 = merge_map(ax["mu_lj_dsaMuon_pt_abseta"][0], G8)
    mm_LJ5_nat = merge_map(ax["mu_lj_dsaMuon_pt"][0], NATIVE)
    mm_nat_ptPf = merge_map(ax["dsaPfMatch_ptPf_ptDsa"][0], NATIVE)
    mm_nat_ptDsa = merge_map(ax["dsaPfMatch_ptPf_ptDsa"][1], NATIVE)
    mm_coarse = merge_map(ax["mu_lj_dsaMuon_pt_phi_nHits"][0], [0.0, 100.0, 200.0, 500.0])
    sl_mat_eta, of_mat_eta = slice_map(ax["dsaPfMatch_ptPf_ptDsa_abseta"][2], ABSETA_BOUNDS, False)
    sl_mat_nh, of_mat_nh = slice_map(ax["dsaPfMatch_ptPf_ptDsa_nHits"][2], NHITS_BOUNDS, True)
    sl_lj_eta, of_lj_eta = slice_map(ax["mu_lj_dsaMuon_pt_abseta"][1], ABSETA_BOUNDS, False)
    sl_lj_nh, of_lj_nh = slice_map(ax["mu_lj_dsaMuon_pt_nHits"][1], NHITS_BOUNDS, True)
    sl_cube_nh, _ = slice_map(ax["mu_lj_dsaMuon_pt_phi_nHits"][2], NHITS_BOUNDS, True)
    sl_cube_eta, _ = slice_map(ax["mu_lj_dsaMuon_pt_abseta_nHits"][1], ABSETA_BOUNDS, False)

    print(f"matrix ptPf -> G8 : {mm_ptPf}")
    print(f"matrix ptDsa-> G8 : {mm_ptDsa}")
    print(f"LJ 20 GeV  -> G8 : {mm_LJ20}")
    print(f"LJ 5 GeV   -> native 12-bin : {[len(g) for g in mm_LJ5_nat]} (group sizes)")
    print(f"coarse cube pT   : {mm_coarse}")
    print(f"matrix |eta| slices: {sl_mat_eta}   nHits slices: {sl_mat_nh} (nHits OF folded in)")
    print(f"LJ fine |eta| slices: {sl_lj_eta}   nHits slices: {sl_lj_nh}")

    sub("expected view-index maps (spec section 2)")
    check("matrix pT merge map", mm_ptPf == [[1, 2], [3, 4], [5, 6], [7], [8], [9, 10], [11], [12], [13]], str(mm_ptPf))
    check("LJ 20 GeV merge map", mm_LJ20 == [[1], [2], [3], [4], [5], [6, 7, 8, 9, 10],
                                             [11, 12, 13, 14, 15],
                                             [16, 17, 18, 19, 20, 21, 22, 23, 24, 25], [26]], "")
    check("LJ 5 GeV -> native group sizes", [len(g) for g in mm_LJ5_nat] ==
          [2, 2, 2, 2, 2, 2, 4, 4, 10, 10, 20, 40, 1], str([len(g) for g in mm_LJ5_nat]))
    check("coarse cube pT map (tail = last inner + OF)",
          mm_coarse == [[1], [2], [3], [4]], str(mm_coarse))
    check("|eta| fine slice map", sl_lj_eta == [list(range(1, 10)), [10, 11, 12], list(range(13, 25))], "")
    check("nHits fine slice map (OF folded into top slice)",
          sl_lj_nh == [[1, 2, 3, 4], [5, 6], [7, 8, 9, 10, 11, 12, 13]], str(sl_lj_nh))

    S_IDX, SB_IDX, TAIL_IDX = [0, 1, 2, 3, 4], [5], [6, 7, 8]
    TAIL300_IDX = [7, 8]
    print(f"\nG8 index sets: S(pT_PF<100)={S_IDX}  SB([100,200))={SB_IDX}  "
          f"TAIL(>=200 incl OF)={TAIL_IDX}  TAIL300={TAIL300_IDX}; "
          f"G8 row 8 = ptPf overflow (>500) -- never a source")

    sub("flow audit (values in the flow slots, summed over the loaded files)")
    fa = {}
    for tag, ST in (("data", D), ("MC", M)):
        for ch in (CH_BEFORE, CH_CR, CH_VR):
            v, _ = ST.get("dsaPfMatch_ptPf_ptDsa", ch)
            v3e, _ = ST.get("dsaPfMatch_ptPf_ptDsa_abseta", ch)
            v3n, _ = ST.get("dsaPfMatch_ptPf_ptDsa_nHits", ch)
            lj, _ = ST.get("mu_lj_dsaMuon_pt", ch)
            l2e, _ = ST.get("mu_lj_dsaMuon_pt_abseta", ch)
            l2n, _ = ST.get("mu_lj_dsaMuon_pt_nHits", ch)
            cph, _ = ST.get("mu_lj_dsaMuon_pt_phi_nHits", ch)
            row = {
                "ptPf_UF": v[0, :].sum(), "ptPf_OF": v[-1, :].sum(),
                "ptDsa_UF": v[:, 0].sum(), "ptDsa_OF": v[:, -1].sum(),
                "LJpt1D_UF": lj[0], "LJpt1D_OF": lj[-1],
                "LJpt20_UF": l2e[0, :].sum(), "LJpt20_OF": l2e[-1, :].sum(),
                "abseta_UF": v3e[:, :, 0].sum(), "abseta_OF": v3e[:, :, -1].sum(),
                "nHits_UF": v3n[:, :, 0].sum(), "nHits_OF": v3n[:, :, -1].sum(),
                "LJabseta_UF": l2e[:, 0].sum(), "LJabseta_OF": l2e[:, -1].sum(),
                "LJnHits_UF": l2n[:, 0].sum(), "LJnHits_OF": l2n[:, -1].sum(),
                "LJphi_UF": cph[:, 0, :].sum(), "LJphi_OF": cph[:, -1, :].sum(),
                "cube_ptOF": cph[-1].sum(),
            }
            fa[f"{tag}/{ch}"] = {k: float(x) for k, x in row.items()}
            print(f"  {tag:4s} {ch:46s} " + " ".join(f"{k}={x:.4g}" for k, x in row.items()))
    RES["flow_audit"] = fa
    zero_keys = ["ptPf_UF", "ptDsa_UF", "LJpt1D_UF", "LJpt20_UF", "abseta_UF", "abseta_OF",
                 "nHits_UF", "LJabseta_UF", "LJabseta_OF", "LJnHits_UF", "LJphi_UF", "LJphi_OF"]
    bad = {k: {kk: vv for kk, vv in r.items() if kk in zero_keys and abs(vv) > 1e-9}
           for k, r in fa.items()}
    bad = {k: v for k, v in bad.items() if v}
    check("all asserted-zero flow slots are zero (pT/|eta|/phi/nHits underflows, "
          "|eta| and phi overflows)", not bad, str(bad) if bad else "12 slots x 6 samples")
    nz = {k: r["nHits_OF"] for k, r in fa.items() if abs(r["nHits_OF"]) > 1e-9}
    check("nHits overflow (>=60) folded into the top nHits slice", True,
          f"non-zero in {nz}" if nz else "empty today; fold is in place regardless")

    # ------------------------------------------------------------------ additivity
    hr("T3  ADDITIVITY AND CONSISTENCY IDENTITIES")
    for tag, ST in (("data", D), ("MC", M)):
        for ch in (CH_BEFORE, CH_CR, CH_VR):
            p2, _ = ST.get("dsaPfMatch_ptPf_ptDsa", ch)
            for nm, hn, sl in (("abseta", "dsaPfMatch_ptPf_ptDsa_abseta", None),
                               ("nHits", "dsaPfMatch_ptPf_ptDsa_nHits", None)):
                v3, _ = ST.get(hn, ch)
                check(f"{tag}/{ch[:22]}: {hn} summed over slice axis == parent 2D",
                      np.allclose(v3.sum(axis=2), p2), f"tot {v3.sum():.6g} vs {p2.sum():.6g}")
            lj1, _ = ST.get("mu_lj_dsaMuon_pt", ch)
            l2e, _ = ST.get("mu_lj_dsaMuon_pt_abseta", ch)
            l2n, _ = ST.get("mu_lj_dsaMuon_pt_nHits", ch)
            cen, _ = ST.get("mu_lj_dsaMuon_pt_abseta_nHits", ch)
            cph, _ = ST.get("mu_lj_dsaMuon_pt_phi_nHits", ch)
            lph, _ = ST.get("mu_lj_dsaMuon_pt_phi", ch)
            c1 = merge_axis(lj1, 0, mm_coarse[:-1] + [mm_coarse[-1]])
            # coarse cube pT view has 3 inner + OF; merge the 1D onto [0,100,200,500]+OF
            mm1D_coarse = merge_map(ax["mu_lj_dsaMuon_pt"][0], [0.0, 100.0, 200.0, 500.0])
            c1 = merge_axis(lj1, 0, mm1D_coarse)
            check(f"{tag}/{ch[:22]}: LJ 1D == 2D(|eta|) == 2D(nHits) == coarse cubes (totals)",
                  np.allclose([lj1.sum(), l2e.sum(), l2n.sum(), cen.sum(), cph.sum(), lph.sum()],
                              lj1.sum()),
                  f"{lj1.sum():.6g}")
            mmL = merge_map(ax["mu_lj_dsaMuon_pt_abseta"][0], [0.0, 100.0, 200.0, 500.0])
            c2e = merge_axis(l2e.sum(axis=1), 0, mmL)
            c2n = merge_axis(l2n.sum(axis=1), 0, mmL)
            ccen = cen.sum(axis=(1, 2))[1:]
            ccph = cph.sum(axis=(1, 2))[1:]
            check(f"{tag}/{ch[:22]}: coarse-pT projections agree (1D, 2D x2, cubes x2)",
                  np.allclose(c1, c2e) and np.allclose(c1, c2n)
                  and np.allclose(c1, ccen) and np.allclose(c1, ccph),
                  f"{np.array2string(np.asarray(c1), precision=1)}")
            check(f"{tag}/{ch[:22]}: pt_phi_nHits summed over nHits == pt_phi rebinned",
                  np.allclose(cph[1:].sum(axis=2), merge_axis(lph, 0, mmL)),
                  "")
    # sum-of-per-file == file-of-sums, done on the era stores
    tot_era = sum(ERA[e].v["mu_lj_dsaMuon_pt"].sum() for e in "ABCD" if ERA[e].n_files)
    check("sum of per-era accumulators == global accumulator (LJ pT total)",
          abs(tot_era - D.v["mu_lj_dsaMuon_pt"].sum()) < 1e-6,
          f"{tot_era:.6g} vs {D.v['mu_lj_dsaMuon_pt'].sum():.6g}")

    # ------------------------------------------------------------------ 2018C regression
    hr("T3b  2018C REGRESSION (run with --subset 2018C to make this the only data loaded)")
    C_ONLY = Store()
    cpath = f"{DATA_DIR}/DoubleMuon_2018C.coffea"
    reg = {}
    if os.path.exists(cpath):
        C_ONLY.add(cpath, True, HISTS)
        anchors = [
            ("Before total/tail", CH_BEFORE, 63040, 2414),
            ("CR     total/tail", CH_CR, 40881, 309),
            ("VR     total/tail", CH_VR, 28303, 187),
            ("cosmic_muons_vtxChi2", CH_COS_VTX, 13781, 1336),
            ("cosmic_muons", CH_COS_SPREAD, 12514, 1242),
        ]
        mm1 = merge_map(C_ONLY.axes["mu_lj_dsaMuon_pt"][0], G8)
        print(f"{'quantity':26s} {'measured':>12s} {'anchor':>10s}")
        for nm, ch, tot_a, tail_a in anchors:
            lj, _ = C_ONLY.get("mu_lj_dsaMuon_pt", ch)
            m8 = merge_axis(lj, 0, mm1)
            tot, tail = lj.sum(), m8[TAIL_IDX].sum()
            print(f"{nm:26s} {tot:9.0f}/{tail:<5.0f} {tot_a:6d}/{tail_a:<5d}")
            check(f"2018C anchor {nm}", tot == tot_a and tail == tail_a,
                  f"{tot:.0f}/{tail:.0f} vs {tot_a}/{tail_a}")
            reg[nm] = {"total": float(tot), "tail": float(tail),
                       "anchor_total": tot_a, "anchor_tail": tail_a}
        mtots = []
        for hn in ("dsaPfMatch_ptPf_ptDsa", "dsaPfMatch_ptPf_ptDsa_abseta",
                   "dsaPfMatch_ptPf_ptDsa_nHits"):
            v, _ = C_ONLY.get(hn, CH_BEFORE)
            mtots.append(float(v.sum()))
        check("2018C matched total 402820 in all three dsaPfMatch_ptPf_ptDsa* hists",
              all(abs(x - 402820) < 1e-6 for x in mtots), str(mtots))
        reg["matched_total"] = mtots
        v, _ = C_ONLY.get("dsaPfMatch_ptPf_ptDsa", CH_BEFORE)
        m8 = merge_axis(merge_axis(v, 0, mm_ptPf), 1, mm_ptDsa)
        num = m8[np.ix_(S_IDX, TAIL_IDX)].sum()
        den = m8[S_IDX, :].sum()
        check("2018C feed-up bookkeeping on G8 == 3137/394257 (the merge-map test)",
              abs(num - 3137) < 1e-6 and abs(den - 394257) < 1e-6, f"{num:.0f}/{den:.0f}")
        reg["feedup_num_den"] = [float(num), float(den)]
        # The spec quotes this anchor on the RAW NATIVE FLOW VIEW (0 = underflow,
        # 1..12 = inner bins, 13 = overflow), which is exactly what `v` already is: the
        # matrix axes ARE the native grid.  It must NOT be applied to a merge_axis-ed
        # array, because merging drops the underflow slot, so the literal indices [:9] and
        # [11:] would then mean pT_PF < 150 and pT_DSA >= 300 (giving 1460/400785).
        check("2018C native-grid feed-up, raw flow view v[0:9,11:] == 3137/394257",
              abs(v[:9, 11:].sum() - 3137) < 1e-6 and abs(v[:9, :].sum() - 394257) < 1e-6,
              f"{v[:9,11:].sum():.0f}/{v[:9,:].sum():.0f}")
        native = merge_axis(merge_axis(v, 0, mm_nat_ptPf), 1, mm_nat_ptDsa)
        check("native-grid remerge is a no-op on the matrix (its axes are already native)",
              abs(native.sum() - v[1:, 1:].sum()) < 1e-6, f"{native.sum():.0f}")
        cen, _ = C_ONLY.get("mu_lj_dsaMuon_pt_abseta_nHits", CH_BEFORE)
        check("2018C SB (Before, coarse cube index 2) == 4706", abs(cen[2].sum() - 4706) < 1e-6,
              f"{cen[2].sum():.0f}")
        check("2018C LJ cube pT overflow (Before) == 375", abs(cen[-1].sum() - 375) < 1e-6,
              f"{cen[-1].sum():.0f}")
        reg["SB_before"] = float(cen[2].sum())
        reg["cube_ptOF_before"] = float(cen[-1].sum())
        vcr, _ = C_ONLY.get("dsaPfMatch_ptPf_ptDsa", CH_CR)
        check("2018C matched ptPf OF / ptDsa OF (CR) == 18 / 471",
              abs(vcr[-1, :].sum() - 18) < 1e-6 and abs(vcr[:, -1].sum() - 471) < 1e-6,
              f"{vcr[-1,:].sum():.0f}/{vcr[:,-1].sum():.0f}")
        l2n, _ = C_ONLY.get("mu_lj_dsaMuon_pt_nHits", CH_CR)
        mmL = merge_map(C_ONLY.axes["mu_lj_dsaMuon_pt_nHits"][0], G8)
        t_nh = merge_axis(merge_axis(l2n, 0, mmL), 1, sl_lj_nh + [of_lj_nh])[TAIL_IDX].sum(axis=0)
        check("2018C CR tail nHits split [0,20)=126, >=20=183",
              abs(t_nh[0] - 126) < 1e-6 and abs(t_nh[1] + t_nh[2] - 183) < 1e-6,
              f"{t_nh[0]:.0f} / {t_nh[1]+t_nh[2]:.0f}")
        reg["CR_tail_nHits_split"] = [float(t_nh[0]), float(t_nh[1] + t_nh[2])]
        cph, _ = C_ONLY.get("mu_lj_dsaMuon_pt_phi_nHits", CH_CR)
        parts = phi_partitions(C_ONLY.axes["mu_lj_dsaMuon_pt_phi_nHits"][1])
        tail_phi = cph[3:].sum(axis=(0, 2))[1:-1]
        for pn, (Vi, Hi) in parts.items():
            Vn, Hn = tail_phi[Vi].sum(), tail_phi[Hi].sum()
            print(f"  2018C CR tail phi partition {pn:9s}: V={Vn:.0f} H={Hn:.0f} "
                  f"(nV={len(Vi)} nH={len(Hi)})")
            if pn == "primary":
                check("2018C CR tail phi V/H == 77/82 with the 4/4 partition",
                      abs(Vn - 77) < 1e-6 and abs(Hn - 82) < 1e-6, f"{Vn:.0f}/{Hn:.0f}")
                reg["phi_VH_primary"] = [float(Vn), float(Hn)]
            else:
                check("the 8/8 partition does NOT give 77/82 (it is a different statistic)",
                      not (abs(Vn - 77) < 1e-6 and abs(Hn - 82) < 1e-6), f"{Vn:.0f}/{Hn:.0f}")
                reg["phi_VH_secondary"] = [float(Vn), float(Hn)]
        check("primary phi partition is 4 V bins and 4 H bins",
              len(parts["primary"][0]) == 4 and len(parts["primary"][1]) == 4, "")
    RES["regression_2018C"] = reg

    # ------------------------------------------------------------------ slicings
    SLICINGS = OrderedDict()
    SLICINGS["abseta"] = {
        "labels": ["|eta| 0.0-0.9", "|eta| 0.9-1.2", "|eta| 1.2-2.4"],
        "mat_hist": "dsaPfMatch_ptPf_ptDsa_abseta", "mat_sl": sl_mat_eta,
        "lj_hist": "mu_lj_dsaMuon_pt_abseta", "lj_sl": sl_lj_eta, "lj_mm": mm_LJ20,
        "resp_hist": "dsaPfMatch_response_abseta", "resp_sl": None}
    SLICINGS["nHits"] = {
        "labels": ["nHits 0-20", "nHits 20-30", "nHits 30-60(+OF)"],
        "mat_hist": "dsaPfMatch_ptPf_ptDsa_nHits", "mat_sl": sl_mat_nh,
        "lj_hist": "mu_lj_dsaMuon_pt_nHits", "lj_sl": sl_lj_nh, "lj_mm": mm_LJ20,
        "resp_hist": "dsaPfMatch_response_nHits", "resp_sl": None}
    SLICINGS["abseta"]["resp_sl"] = slice_map(ax["dsaPfMatch_response_abseta"][1],
                                              ABSETA_BOUNDS, False)[0]
    SLICINGS["nHits"]["resp_sl"] = slice_map(ax["dsaPfMatch_response_nHits"][1],
                                             NHITS_BOUNDS, True)[0]

    def get_M8(ST, ch, slicing, s):
        """Merged migration matrix + variances for one channel and slice.
        s is either an int (slice index) or a list of slice indices to merge."""
        cfg = SLICINGS[slicing]
        v, w2 = ST.get(cfg["mat_hist"], ch)
        idxs = cfg["mat_sl"][s] if isinstance(s, int) else \
            [i for k in s for i in cfg["mat_sl"][k]]
        v = pick_axis(v, 2, idxs)
        w2 = pick_axis(w2, 2, idxs)
        v = merge_axis(merge_axis(v, 0, mm_ptPf), 1, mm_ptDsa)
        w2 = merge_axis(merge_axis(w2, 0, mm_ptPf), 1, mm_ptDsa)
        return v, w2

    def get_L(ST, ch, slicing, s):
        cfg = SLICINGS[slicing]
        v, w2 = ST.get(cfg["lj_hist"], ch)
        idxs = cfg["lj_sl"][s] if isinstance(s, int) else \
            [i for k in s for i in cfg["lj_sl"][k]]
        v = pick_axis(v, 1, idxs)
        w2 = pick_axis(w2, 1, idxs)
        return merge_axis(v, 0, cfg["lj_mm"]), merge_axis(w2, 0, cfg["lj_mm"])

    def get_M8_incl(ST, ch):
        v, w2 = ST.get("dsaPfMatch_ptPf_ptDsa", ch)
        return (merge_axis(merge_axis(v, 0, mm_nat_ptPf), 1, mm_nat_ptDsa),
                merge_axis(merge_axis(w2, 0, mm_nat_ptPf), 1, mm_nat_ptDsa))

    def get_L_incl(ST, ch):
        v, w2 = ST.get("mu_lj_dsaMuon_pt", ch)
        return merge_axis(v, 0, mm_LJ5_nat), merge_axis(w2, 0, mm_LJ5_nat)

    # ------------------------------------------------------------------ T4 transfer
    hr("T4  TRANSFER PROBABILITIES  T(j|i) = P(measured in j | true/PF bin i)")
    print("rows: G8 source (PF) bins; targets: SB=[100,200), TAIL=>=200 incl OF, "
          "TAIL300=>=300 incl OF.  eps = P(measured < 100 | i).")
    print("Row guards: data n(i) >= 200 and MC N_eff(n(i)) >= 30.  headroom = "
          "84% Garwood upper limit on an empty/near-empty T(TAIL|i) times n_hat(i).")
    G8LAB = ["[0,20)", "[20,40)", "[40,60)", "[60,80)", "[80,100)",
             "[100,200)", "[200,300)", "[300,500)", ">500(OF)"]
    T4 = {}
    MONO = {}
    for slicing, cfg in SLICINGS.items():
        for si, lab in enumerate(cfg["labels"]):
            for cn, ch in FOLD_CHANNELS:
                dM, dMv = get_M8(D, ch, slicing, si)
                mM, mMv = get_M8(M, ch, slicing, si)
                # per-SAMPLE MC matrices for this (channel, slicing, slice): the row guard
                # protects n(i), but Pred is set by the NUMERATOR cell m(i->TAIL), whose
                # N_eff is far smaller and can be carried by a single high-weight sample.
                perS = {sn: get_M8(ST, ch, slicing, si)[0] for sn, ST in MCS.items()}
                dr = fold_core(dM, np.zeros(9), S_IDX, SB_IDX, TAIL_IDX, 0)
                mr = fold_core(mM, np.zeros(9), S_IDX, SB_IDX, TAIL_IDX, 0)
                key = f"{cn}/{slicing}/{lab}"
                T4[key] = []
                if cn == "CR":
                    sub(f"T4  {cn}  {slicing}  {lab}")
                    print(f"{'src':10s} {'data n':>9s} {'d T(SB)':>10s} {'d T(TAIL)':>11s} "
                          f"{'d T(300)':>10s} {'d eps':>8s} | {'MC n':>10s} {'Neff(n)':>8s} "
                          f"{'m T(TAIL)':>11s} {'ratio d/m':>10s} {'Neff(num)':>9s} "
                          f"{'dominant MC sample':>22s} {'share':>6s} {'headroom':>9s} {'guard':>6s}")
                for i in range(9):
                    dn = float(dM[i].sum())
                    mn = float(mM[i].sum())
                    mnv = float(mMv[i].sum())
                    neff = mn * mn / mnv if mnv > 0 else 0.0
                    dtl = float(dM[i, TAIL_IDX].sum())
                    dtlv = float(dMv[i, TAIL_IDX].sum())
                    mtl = float(mM[i, TAIL_IDX].sum())
                    mtlv = float(mMv[i, TAIL_IDX].sum())
                    dT = dtl / dn if dn > 0 else 0.0
                    mT = mtl / mn if mn > 0 else 0.0
                    edT = wbinom_err(dtl, dtlv, dn, float(dMv[i].sum()))
                    emT = wbinom_err(mtl, mtlv, mn, mnv)
                    ratio = dT / mT if mT > 0 else float("nan")
                    head = (1.84 / dn) if (dn > 0 and dtl == 0) else 0.0
                    # numerator-cell diagnostics
                    neff_num = mtl * mtl / mtlv if mtlv > 0 else 0.0
                    contrib = {sn: float(sv[i, TAIL_IDX].sum()) for sn, sv in perS.items()}
                    tot_c = sum(contrib.values())
                    if tot_c > 0:
                        dom_s = max(contrib, key=lambda k_: contrib[k_])
                        dom_f = contrib[dom_s] / tot_c
                    else:
                        dom_s, dom_f = "-", 0.0
                    rec = {"src": G8LAB[i], "data_n": dn, "data_T_SB": float(dr["T"][i, SB_IDX].sum()),
                           "data_T_TAIL": dT, "data_T_TAIL_err": edT,
                           "data_T_300": float(dr["T"][i, TAIL300_IDX].sum()),
                           "data_eps": float(dr["T"][i, S_IDX].sum()),
                           "mc_n": mn, "mc_neff": neff, "mc_T_TAIL": mT,
                           "mc_T_TAIL_err": emT, "ratio_d_over_m": ratio,
                           "mc_num": mtl, "mc_neff_num": neff_num,
                           "mc_dom_sample": dom_s, "mc_dom_share": dom_f,
                           "mc_num_by_sample": contrib,
                           "headroom_T": head,
                           "guard_ok": bool(dn >= 200 and neff >= 30),
                           "guard_num_ok": bool(neff_num >= 30)}
                    T4[key].append(rec)
                    if cn == "CR":
                        print(f"{G8LAB[i]:10s} {dn:9.0f} {dr['T'][i,SB_IDX].sum():10.5f} "
                              f"{dT:8.5f}+-{edT:5.5f} {dr['T'][i,TAIL300_IDX].sum():10.5f} "
                              f"{dr['T'][i,S_IDX].sum():8.5f} | {mn:10.1f} {neff:8.1f} "
                              f"{mT:8.5f}+-{emT:5.5f} {ratio:10.3f} {neff_num:9.1f} "
                              f"{dom_s:>22s} {dom_f:6.2f} {head:9.5f} "
                              f"{'ok' if (dn>=200 and neff>=30) else 'MERGE':>6s}")
                tt = [T4[key][i]["data_T_TAIL"] for i in S_IDX]
                tm = [T4[key][i]["mc_T_TAIL"] for i in S_IDX]
                mono_d = all(tt[k] <= tt[k + 1] * 1.0000001 for k in range(len(tt) - 1))
                mono_m = all(tm[k] <= tm[k + 1] * 1.0000001 for k in range(len(tm) - 1))
                MONO[key] = {"data": bool(mono_d), "mc": bool(mono_m),
                             "T_data": tt, "T_mc": tm}
                if cn == "CR":
                    print(f"  T(TAIL|i) monotone rising in i:  DATA {mono_d} :: "
                          f"{['%.5f' % x for x in tt]}")
                    print(f"  T(TAIL|i) monotone rising in i:  MC   {mono_m} :: "
                          f"{['%.5f' % x for x in tm]}   "
                          f"(a mismeasurement probability that FALLS with pT_PF is "
                          f"unphysical and flags a low-N_eff cell)")
    RES["T4_transfer"] = T4
    guard_bad = [(k, r["src"]) for k, rs in T4.items() for r in rs
                 if r["src"] in G8LAB[:5] and not r["guard_ok"]]
    check("row guards (data n >= 200 and MC N_eff >= 30) on every source row of every slice",
          not guard_bad, f"{len(guard_bad)} violations: {guard_bad[:8]}" if guard_bad else "45 rows")
    nonmono = [k for k, v in MONO.items() if not v["mc"]]
    check("MC T(TAIL|i) monotone rising in the source rows (same flag the data get)",
          not nonmono, f"non-monotone in MC: {nonmono}" if nonmono else "18 (channel,slicing,slice) sets",
          kind="expectation")
    nonmono_d = [k for k, v in MONO.items() if not v["data"]]
    check("data T(TAIL|i) monotone rising in the source rows",
          not nonmono_d, f"non-monotone in data: {nonmono_d}" if nonmono_d else "18 sets",
          kind="expectation")
    numguard = [(k, r["src"], r["mc_neff_num"], r["mc_dom_sample"], r["mc_dom_share"])
                for k, rs in T4.items() for r in rs
                if r["src"] in G8LAB[:5] and not r["guard_num_ok"]]
    check("MC NUMERATOR-cell guard N_eff(m(i->TAIL)) >= 30 on every source row "
          "(the quantity that actually sets Pred_MC, not the row total n(i))",
          not numguard,
          f"{len(numguard)} violations, e.g. " + "; ".join(
              f"{k}/{s}: Neff={n:.1f} dom={ds}({dsh:.2f})"
              for k, s, n, ds, dsh in numguard[:6]) if numguard else "45 rows",
          kind="method")
    domhi = [(k, r["src"], r["mc_dom_sample"], r["mc_dom_share"])
             for k, rs in T4.items() for r in rs
             if r["src"] in G8LAB[:5] and r["mc_dom_share"] > 0.8]
    check("no MC numerator cell is more than 80 % carried by a single sample",
          not domhi, f"{len(domhi)} cells: " + "; ".join(
              f"{k}/{s}:{ds}={f_:.2f}" for k, s, ds, f_ in domhi[:8]) if domhi else "45 rows",
          kind="method")
    RES["T4_mono"] = MONO

    # ---------------------------------------------------- T4c smoothed / anomalous cell
    sub("T4c  SMOOTHED TRANSFER  T(TAIL|i) = a / pT_PF  AND THE ANOMALOUS-CELL TEST "
        "(spec section 3; fix round, finding 12)")
    print("The nominal fold uses the cells as measured (pre-registered, no smoothing).  "
          "These are VARIANTS, printed to expose how much of Pred rests on one cell.")
    G8CEN = np.array([10.0, 30.0, 50.0, 70.0, 90.0])          # source-row pT_PF centres
    T4C = {}
    for slicing, cfg in SLICINGS.items():
        for si, lab in enumerate(cfg["labels"]):
            for tag, ST, wf in (("data", D, True), ("mc", M, False)):
                Mv, Mw = get_M8(ST, CH_CR, slicing, si)
                Lv, Lw = get_L(ST, CH_CR, slicing, si)
                r = fold_core(Mv, Lv, S_IDX, SB_IDX, TAIL_IDX, 3)
                Ntrue = r["Ntrue"]
                nrow = Mv.sum(axis=1)[S_IDX]
                num = Mv[np.ix_(S_IDX, TAIL_IDX)].sum(axis=1)
                numv = Mw[np.ix_(S_IDX, TAIL_IDX)].sum(axis=1)
                rowv = Mw[S_IDX].sum(axis=1)
                T = np.where(nrow > 0, num / np.maximum(nrow, 1e-300), 0.0)
                sT = np.array([wbinom_err(num[k], numv[k], nrow[k], rowv[k])
                               for k in range(len(T))])
                pred_nom = float((T * Ntrue).sum())
                # weighted least squares for T = a / pT
                w = np.where(sT > 0, 1.0 / sT ** 2, 0.0)
                den = float(np.sum(w / G8CEN ** 2))
                a = float(np.sum(w * T / G8CEN) / den) if den > 0 else 0.0
                Tsm = a / G8CEN
                pred_sm = float((Tsm * Ntrue).sum())
                # anomalous-cell test: any source row whose T exceeds its upper neighbour
                Tfix = T.copy()
                fixed = []
                for k in range(len(T) - 1):
                    if T[k] > T[k + 1] * 1.0000001:
                        if k == 0:
                            Tfix[k] = a / G8CEN[k]          # no lower neighbour: use the fit
                        else:
                            Tfix[k] = math.sqrt(max(T[k - 1], 1e-30) * max(T[k + 1], 1e-30))
                        fixed.append(G8LAB[S_IDX[k]])
                pred_fix = float((Tfix * Ntrue).sum())
                contrib = T * Ntrue
                top = int(np.argmax(contrib))
                T4C[f"{slicing}/{lab}/{tag}"] = {
                    "pred_nominal": pred_nom, "pred_smoothed": pred_sm,
                    "pred_anom_fixed": pred_fix, "a_fit": a,
                    "cells_fixed": fixed,
                    "top_row": G8LAB[S_IDX[top]],
                    "top_row_share": float(contrib[top] / max(pred_nom, 1e-9)),
                    "T": T.tolist(), "T_smoothed": Tsm.tolist()}
                print(f"  {slicing:7s} {lab:17s} {tag:4s} Pred(nominal)={pred_nom:9.1f}  "
                      f"Pred(T ~ a/pT)={pred_sm:9.1f}  Pred(anom cell repaired)={pred_fix:9.1f}  "
                      f"repaired={fixed if fixed else 'none'}  "
                      f"top row {G8LAB[S_IDX[top]]} carries {contrib[top]/max(pred_nom,1e-9):.2f} of Pred")
    RES["T4c_smoothed"] = T4C

    # ------------------------------------------------------------------ unit tests
    hr("T5  CIRCULARITY: UNIT TESTS AND THE ITERATION")
    dM, dMv = get_M8(D, CH_CR, "nHits", 0)
    colsS = dM[:, :].sum(axis=0)                       # all-row column sums
    colsSS = dM[np.ix_(S_IDX, range(9))].sum(axis=0)   # row-restricted column sums
    pB_sb, pB_tl, nhB = method_B(dM, colsS, S_IDX, SB_IDX, TAIL_IDX)
    target = dM[np.ix_(S_IDX, TAIL_IDX)].sum()
    rel = abs(pB_tl - target) / max(target, 1e-30)
    check("unit test 1: method B on the matched sample's own column sums returns "
          "sum_{i in S} M[i,TAIL] to machine precision", rel < 1e-12,
          f"{pB_tl:.10f} vs {target:.10f} rel={rel:.2e}")
    fr = fold_core(dM, colsSS, S_IDX, SB_IDX, TAIL_IDX, 3, conv_iters=20000, conv_tol=1e-10)
    ntrue3 = fr["Ntrue"]
    nrow = dM[S_IDX].sum(axis=1)
    dev = np.max(np.abs(ntrue3 - nrow) / np.maximum(nrow, 1e-30))
    devc = np.max(np.abs(fr["Ntrue_conv"] - nrow) / np.maximum(nrow, 1e-30))
    # n_row IS an exact fixed point of the recursion (pi ~ n_row makes P(i|k) = M[i,k]/L(k),
    # hence N_meas(i) = sum_k M[i,k] and N_true(i) = n_row(i) identically).  What this test
    # actually measures is whether 3 iterations REACH it; both distances are reported.
    check("unit test 2: the nominal iteration on the ROW-RESTRICTED matched sub-matrix "
          "recovers n(i) within 1 % by iteration 3", dev < 0.01,
          f"max rel dev at iteration 3 = {dev:.4f}; after {fr['n_conv']} iterations = "
          f"{devc:.4f} (the fixed point is exact; 3 iterations do not reach it)",
          kind="method")
    fr2 = fold_core(dM, colsS, S_IDX, SB_IDX, TAIL_IDX, 3)
    defc = fr2["Ntrue"] - nrow
    print("  unit test 3 (unrestricted matched columns): per-source-bin deficit "
          "= feed-DOWN size (genuine >=100 GeV muons measured <100), "
          + ", ".join(f"{G8LAB[i]}:{defc[k]:+.1f}" for k, i in enumerate(S_IDX)))
    check("unit test 3: the unrestricted fold does NOT close (feed-down exists and is exposed)",
          np.max(np.abs(defc / np.maximum(nrow, 1e-30))) > 0.01,
          f"max rel deficit {np.max(np.abs(defc/np.maximum(nrow,1e-30))):.4f}")
    dL, dLv = get_L(D, CH_CR, "nHits", 0)
    fl = fold_core(dM, dL, S_IDX, SB_IDX, TAIL_IDX, 3)
    check("unit test 4: the LJ fold differs from the matched self-fold",
          abs(fl["preds"][3][1] - fr["preds"][3][1]) / max(fr["preds"][3][1], 1e-9) > 0.01,
          f"LJ {fl['preds'][3][1]:.2f} vs matched-self {fr['preds'][3][1]:.2f}")
    print("\nCIRCULARITY AUDIT (R2): Pred is built from the MATCHED counts (T) and the LJ "
          "LOW-pT counts (L, pT<100); Obs is the LJ >=200 count.  No algebraic identity "
          "links the two sides.  CAVEAT: a dR<0.2-matched but segment-unmatched DSA muon "
          "sits in BOTH M and L/Obs; the overlap is not measurable from the booked hists.")

    # ------------------------------------------------------------------ the fold
    hr("T6/T7  FOLD PREDICTIONS AND CLOSURE")

    def run_fold(ST, ch, slicing, s, nt, wf_data, S_i=S_IDX, T_i=TAIL_IDX, nit=3,
                 mat_ch=None, rng_=None):
        Mv, Mw = get_M8(ST, mat_ch or ch, slicing, s)
        Lv, Lw = get_L(ST, ch, slicing, s)
        r = fold_core(Mv, Lv, S_i, SB_IDX, T_i, 3, conv_iters=CONV_ITERS)
        rB = method_B(Mv, Lv, S_i, SB_IDX, T_i)
        out = {"pred_sb": {k: float(v[0]) for k, v in r["preds"].items()},
               "pred_tl": {k: float(v[1]) for k, v in r["preds"].items()},
               "predB_sb": float(rB[0]), "predB_tl": float(rB[1]),
               "obs_sb": float(Lv[SB_IDX].sum()), "obs_tl": float(Lv[T_i].sum()),
               "obs_sb_var": float(Lw[SB_IDX].sum()), "obs_tl_var": float(Lw[T_i].sum()),
               "T_TAIL_from_SB": float((Mv[SB_IDX[0], T_i].sum())
                                       / max(Mv[SB_IDX[0]].sum(), 1e-300)),
               "n_conv": int(r.get("n_conv", 3)),
               "conv_ok": bool(r.get("conv_ok", False)),
               "conv_min_frac": float(r.get("conv_min_frac", float("nan"))),
               "conv_dstate": float(r.get("conv_dstate", float("nan"))),
               "eps": r["eps"].tolist(), "Ntrue": r["Ntrue"].tolist(),
               "L_src": Lv[S_i].tolist(), "L_bulk": float(Lv[:5].sum()),
               "row_contrib": (r["T_TL"] * r["Ntrue"]).tolist()}
        if nt and rng_ is not None:
            wf_m = wbar_of(Mv, Mw)
            wf_l = wbar_of(Lv, Lw)
            Mt = gamma_draw(rng_, Mv, Mw, nt, wf_m if not wf_data else 1.0)
            Lt = gamma_draw(rng_, Lv, Lw, nt, wf_l if not wf_data else 1.0)
            rt = fold_core(Mt, Lt, S_i, SB_IDX, T_i, 3)
            out["toy_pred_sb"] = rt["preds"][3][0]
            out["toy_pred_tl"] = rt["preds"][3][1]
            out["toy_pred_tl_i1"] = rt["preds"][1][1]
        return out

    def analytic_pieces(ST, ch, slicing, s):
        Mv, Mw = get_M8(ST, ch, slicing, s)
        Lv, Lw = get_L(ST, ch, slicing, s)
        n = Mv.sum(axis=1)
        T = Mv[:, TAIL_IDX].sum(axis=1) / np.maximum(n, 1e-300)
        mv = Mw[:, TAIL_IDX].sum(axis=1)
        varT = np.array([wbinom_err(Mv[i, TAIL_IDX].sum(), mv[i], n[i], Mw[i].sum()) ** 2
                         for i in range(9)])
        a = float((T[S_IDX] ** 2 * Lw[S_IDX]).sum())
        b = float((Lv[S_IDX] ** 2 * varT[S_IDX]).sum())
        return a, b

    FOLD = {}
    for cn, ch in FOLD_CHANNELS:
        for slicing, cfg in SLICINGS.items():
            for si, lab in enumerate(cfg["labels"]):
                key = (cn, slicing, si)
                d = run_fold(D, ch, slicing, si, args.toys, True, rng_=rng)
                m = run_fold(M, ch, slicing, si, args.toys, False, rng_=rng)
                dB = run_fold(D, ch, slicing, si, 0, True, mat_ch=CH_BEFORE)
                FOLD[key] = {"data": d, "mc": m, "data_beforeM": dB, "label": lab}

    # inclusive folds (three variants of one number)
    INCL = {}
    for cn, ch in FOLD_CHANNELS:
        preds_sb, preds_tl = [], []
        for slicing, cfg in SLICINGS.items():
            s_sb = sum(FOLD[(cn, slicing, i)]["data"]["pred_sb"][3] for i in range(3))
            s_tl = sum(FOLD[(cn, slicing, i)]["data"]["pred_tl"][3] for i in range(3))
            preds_sb.append(s_sb)
            preds_tl.append(s_tl)
        Mv, Mw = get_M8_incl(D, ch)
        Lv, Lw = get_L_incl(D, ch)
        nb = Mv.shape[0]
        Snat = list(range(0, 8))
        SBnat = [8, 9]
        TLnat = [10, 11, 12]
        rn = fold_core(Mv, Lv, Snat, SBnat, TLnat, 3)
        preds_sb.append(float(rn["preds"][3][0]))
        preds_tl.append(float(rn["preds"][3][1]))
        obs_sb = float(Lv[SBnat].sum())
        obs_tl = float(Lv[TLnat].sum())
        INCL[cn] = {"pred_sb": preds_sb, "pred_tl": preds_tl,
                    "obs_sb": obs_sb, "obs_tl": obs_tl,
                    "sigma_bin_sb": float(np.max(preds_sb) - np.min(preds_sb)),
                    "sigma_bin_tl": float(np.max(preds_tl) - np.min(preds_tl))}
        print(f"\nINCLUSIVE ({cn}) three predictions of one number "
              f"[sum|eta|, sum nHits, native unsliced]:")
        print(f"   TAIL: {['%.1f' % x for x in preds_tl]}  obs={obs_tl:.0f}  "
              f"sigma_bin={INCL[cn]['sigma_bin_tl']:.1f}")
        print(f"   SB  : {['%.1f' % x for x in preds_sb]}  obs={obs_sb:.0f}  "
              f"sigma_bin={INCL[cn]['sigma_bin_sb']:.1f}")
    RES["inclusive"] = INCL
    check("the three inclusive predictions (sum|eta|, sum nHits, native unsliced) agree "
          "within 25 % of their mean",
          all(INCL[c]["sigma_bin_tl"] / max(np.mean(INCL[c]["pred_tl"]), 1e-9) < 0.25
              for c in INCL),
          ", ".join(f"{c}:{INCL[c]['sigma_bin_tl']/max(np.mean(INCL[c]['pred_tl']),1e-9):.3f}"
                    for c in INCL))

    # closure table
    T7 = {}
    sub("T7  CLOSURE  (Obs, Pred, R = Obs/Pred - 1; MC self-closure with N_eff)")
    print(f"toys per entry = {args.toys}.  Pred_cv = the SAME EM recursion continued to its "
          f"fixed point (the pre-registered stop is iteration 3, which it does NOT reach).")
    print(f"R uses Pred(3); R_cv uses Pred(conv); the pair is the honest bracket.  "
          f"+-full = stat (+) circ (+) bin (+) chan.  R_bc = D - 1 = the closure after "
          f"transporting the method's own MC bias multiplicatively.")
    print(f"z_g / zc_g are GAUSSIAN equivalents of the toy ensemble (the tail-count p "
          f"saturates at 1/(N+1), which censors z at +-{z_from_p(1.0/(args.toys+1)):.2f}); "
          f"the censored values are printed as z / z_cons and marked with '*'.")
    print(f"{'chan':6s} {'slicing':7s} {'slice':17s} {'tgt':5s} {'Obs':>8s} {'Pred':>9s} "
          f"{'Pred_cv':>9s} {'+-stat':>7s} {'+-syst':>7s} {'R':>8s} {'+-full':>7s} "
          f"{'R_cv':>8s} {'z_g':>7s} {'R_cons':>8s} {'zc_g':>7s} {'R_MC':>8s} {'+-':>6s} "
          f"{'NeffMC':>7s} {'D':>7s} {'R_bc':>8s}")
    for cn, ch in FOLD_CHANNELS:
        for slicing, cfg in SLICINGS.items():
            for si, lab in enumerate(cfg["labels"]):
                F = FOLD[(cn, slicing, si)]
                d, m, dB = F["data"], F["mc"], F["data_beforeM"]
                sigbin = INCL[cn]["sigma_bin_tl"] * (d["pred_tl"][3]
                                                     / max(sum(FOLD[(cn, slicing, i)]["data"]["pred_tl"][3]
                                                               for i in range(3)), 1e-9))
                sigbin_sb = INCL[cn]["sigma_bin_sb"] * (d["pred_sb"][3]
                                                        / max(sum(FOLD[(cn, slicing, i)]["data"]["pred_sb"][3]
                                                                  for i in range(3)), 1e-9))
                for tgt in ("SB", "TAIL"):
                    pk = "pred_sb" if tgt == "SB" else "pred_tl"
                    ok = "obs_sb" if tgt == "SB" else "obs_tl"
                    tk = "toy_pred_sb" if tgt == "SB" else "toy_pred_tl"
                    Obs, Pred = d[ok], d[pk][3]
                    toy = d[tk]
                    stat = float(np.std(toy))
                    # ---- CIRCULARITY: TWO-SIDED BRACKET (fix round, finding 7) --------
                    # The pre-registered nominal is iteration 3, and it is KEPT as the
                    # nominal (moving a pre-registered choice after seeing the data would
                    # be worse).  But the EM recursion has not converged there and
                    # Pred(conv) < Pred(3) almost everywhere, so a ONE-SIDED UPWARD
                    # sigma_circ would be applied against the only signed, measured bias
                    # the study has.  The bracket now spans the whole iteration family
                    # {0 = naive, 1, 3, conv}: Pred_up (which can only lower R_cons, i.e.
                    # is conservative AGAINST triggering) takes its upper end, and the
                    # LOWER end is reported as an equal member of the pair.
                    pfam = [d[pk][0], d[pk][1], d[pk][3], d[pk].get("conv", d[pk][3])]
                    pred_hi = float(max(pfam))
                    pred_lo = float(max(min(pfam), 1e-9))
                    pred_conv = float(d[pk].get("conv", d[pk][3]))
                    sig_circ_up = max(pred_hi - Pred, 0.0)
                    sig_circ_dn = max(Pred - pred_lo, 0.0)
                    sig_circ_prior = abs(d[pk][3] - d[pk][1])
                    sig_circ_conv = abs(pred_conv - Pred)
                    sig_circ = max(sig_circ_up, sig_circ_dn)       # symmetrised, for sigma_R
                    sig_chan = abs(d[pk][3] - dB[pk][3])
                    sbn = sigbin if tgt == "TAIL" else sigbin_sb
                    R = Obs / Pred - 1 if Pred > 0 else float("nan")
                    R_conv = Obs / pred_conv - 1 if pred_conv > 0 else float("nan")
                    R_bracket = [Obs / pred_hi - 1 if pred_hi > 0 else float("nan"),
                                 Obs / pred_lo - 1 if pred_lo > 0 else float("nan")]
                    sR_pred = float(np.std(Obs / np.maximum(toy, 1e-9) - 1))
                    sR_obs = math.sqrt(max(Obs, 0)) / Pred if Pred > 0 else 0.0
                    sR = math.hypot(sR_pred, sR_obs)               # STAT ONLY
                    # ---- sigma_R now carries the systematics too (fix round, finding 11)
                    sR_syst = math.sqrt(sig_circ ** 2 + sbn ** 2 + sig_chan ** 2) \
                        / max(Pred, 1e-9)
                    sR_full = math.hypot(sR, sR_syst)
                    obs_toy = rng.poisson(np.maximum(toy, 1e-9))
                    p_raw = float((obs_toy >= Obs).mean())
                    clamped = bool(p_raw < 1.0 / (len(toy) + 1) or p_raw > 1 - 1.0 / (len(toy) + 1))
                    p = min(max(p_raw, 1.0 / (len(toy) + 1)), 1 - 1.0 / (len(toy) + 1))
                    z = z_from_p(p)
                    # ---- UNCENSORED significance (fix round, finding 3) ---------------
                    # A tail-count p from N toys saturates at 1/(N+1), so z is clipped to
                    # +-3.29 at N = 2000 by construction.  z_gauss is the Gaussian
                    # equivalent of the SAME posterior-predictive ensemble; it is what the
                    # decision gates use, and both are printed with the toy count.
                    mu_t, sd_t = float(np.mean(obs_toy)), float(np.std(obs_toy))
                    z_gauss = (Obs - mu_t) / sd_t if sd_t > 0 else float("nan")
                    # MC
                    ObsM, PredM = m[ok], m[pk][3]
                    RM = ObsM / PredM - 1 if PredM > 0 else float("nan")
                    toyM = m[tk]
                    sRM = float(np.std(ObsM / np.maximum(toyM, 1e-9) - 1))
                    ovar = m["obs_sb_var"] if tgt == "SB" else m["obs_tl_var"]
                    neffM = ObsM * ObsM / ovar if ovar > 0 else 0.0
                    sRM = math.hypot(sRM, math.sqrt(ovar) / PredM if PredM > 0 else 0.0)
                    Dd = (1 + R) / (1 + RM) if RM > -1 else float("nan")
                    # ---- METHOD-BIAS-CORRECTED closure (fix round, findings 1 and 8) ---
                    # R_MC != 0 means the fold does not reproduce its OWN MC.  R_MC < 0 is
                    # OVER-prediction (Pred_MC > Obs_MC); R_MC > 0 is UNDER-prediction.
                    # Transporting the bias multiplicatively, Pred_corr = Pred*(1 + R_MC),
                    # hence R_bc = D - 1 identically.  Reported BESIDE R, never instead of
                    # it, and the direction is written out in words in the readout.
                    R_bc = Dd - 1 if np.isfinite(Dd) else float("nan")
                    sR_bc = (abs(1 + R_bc) * math.hypot(sR_full / max(abs(1 + R), 1e-9),
                                                        sRM / max(abs(1 + RM), 1e-9))
                             if np.isfinite(Dd) else float("nan"))
                    sig_method = abs(R - R_bc)
                    mbias_flag = bool(abs(RM) > 0.3 and abs(RM) > 2 * max(sRM, 1e-9))
                    # conservative
                    Gmc = ObsM - PredM
                    sG = float(np.std(ObsM - toyM))
                    kap_bulk = d["L_bulk"] / max(m["L_bulk"], 1e-9)
                    Gup = kap_bulk * max(Gmc + sG, 0.0)
                    ff = d["T_TAIL_from_SB"] * max(d["obs_sb"] - d["pred_sb"][3], 0.0) \
                        if tgt == "TAIL" else 0.0
                    Pred_up = Pred + math.sqrt(sig_circ_up ** 2 + sbn ** 2 + sig_chan ** 2) + ff
                    Rc = (Obs - Pred_up - Gup) / Pred if Pred > 0 else float("nan")
                    obs_toy_c = rng.poisson(np.maximum(toy + (Pred_up - Pred) + Gup, 1e-9))
                    pc_raw = float((obs_toy_c >= Obs).mean())
                    clamped_c = bool(pc_raw < 1.0 / (len(toy) + 1)
                                     or pc_raw > 1 - 1.0 / (len(toy) + 1))
                    pc = min(max(pc_raw, 1.0 / (len(toy) + 1)), 1 - 1.0 / (len(toy) + 1))
                    zc = z_from_p(pc)
                    mu_c, sd_c = float(np.mean(obs_toy_c)), float(np.std(obs_toy_c))
                    zc_gauss = (Obs - mu_c) / sd_c if sd_c > 0 else float("nan")
                    # the excess (as a fraction of Pred) that condition (a) would need
                    trig = (Pred_up - Pred + Gup) / max(Pred, 1e-9) + 1.0 / 3.0
                    rec = {"obs": Obs, "pred": Pred, "pred_i0": d[pk][0], "pred_i1": d[pk][1],
                           "pred_i2": d[pk][2], "predB": d["predB_sb" if tgt == "SB" else "predB_tl"],
                           "pred_conv": pred_conv, "pred_hi": pred_hi, "pred_lo": pred_lo,
                           "n_conv": d.get("n_conv", 3),
                           "conv_ok": d.get("conv_ok", False),
                           "conv_min_frac": d.get("conv_min_frac", float("nan")),
                           "stat": stat, "sig_circ": sig_circ,
                           "sig_circ_up": sig_circ_up, "sig_circ_dn": sig_circ_dn,
                           "sig_circ_prior": sig_circ_prior, "sig_circ_conv": sig_circ_conv,
                           "sig_bin": sbn,
                           "sig_chan": sig_chan, "R": R, "R_conv": R_conv,
                           "R_bracket": R_bracket, "sR": sR, "sR_syst": sR_syst,
                           "sR_full": sR_full, "sR_pred": sR_pred,
                           "sR_obs": sR_obs, "p": p, "p_raw": p_raw, "z": z,
                           "z_clamped": clamped, "z_gauss": z_gauss, "n_toys": int(len(toy)),
                           "R_cons": Rc, "z_cons": zc, "z_cons_clamped": clamped_c,
                           "z_cons_gauss": zc_gauss, "p_cons_raw": pc_raw,
                           "trigger_excess_frac": trig,
                           "Pred_up": Pred_up, "Gup": Gup, "Gup_over_pred": Gup / max(Pred, 1e-9),
                           "kappa_bulk": kap_bulk, "feedfwd": ff,
                           "obs_mc": ObsM, "pred_mc": PredM, "R_mc": RM, "sR_mc": sRM,
                           "neff_mc": neffM, "D": Dd, "R_bc": R_bc, "sR_bc": sR_bc,
                           "sig_method": sig_method, "mc_bias_flag": mbias_flag,
                           "mc_direction": ("OVER-predicts" if RM < 0 else "under-predicts"),
                           "mc_factor": (PredM / max(ObsM, 1e-9)) if RM < 0
                                        else (ObsM / max(PredM, 1e-9)),
                           "toyR": (Obs / np.maximum(toy, 1e-9) - 1)}
                    T7[(cn, slicing, si, tgt)] = rec
                    print(f"{cn:6s} {slicing:7s} {lab:17s} {tgt:5s} {Obs:8.0f} {Pred:9.1f} "
                          f"{pred_conv:9.1f} {stat:7.1f} {sR_syst*Pred:7.1f} {R:8.3f} "
                          f"{sR_full:7.3f} {R_conv:8.3f} {z_gauss:7.2f} {Rc:8.3f} "
                          f"{zc_gauss:7.2f} {RM:8.3f} {sRM:6.3f} {neffM:7.1f} {Dd:7.3f} "
                          f"{R_bc:8.3f}"
                          + ("  [MC-BIAS]" if mbias_flag else "")
                          + ("  [z clamp]" if (clamped or clamped_c) else ""))
    RES["T7_closure"] = {f"{k[0]}/{k[1]}/{k[2]}/{k[3]}":
                         {kk: (vv if not isinstance(vv, np.ndarray) else None)
                          for kk, vv in v.items() if kk != "toyR"}
                         for k, v in T7.items()}
    # ---- spec section 8 HARD TRANSPORT FAILURE (fix round, finding 10) ---------------
    # "Pred_SB > Obs_SB in any slice is a hard transport failure there."  It was never
    # coded.  It is load-bearing twice: it is the floor on any tail-closure claim in the
    # same slice (the sideband has ~3x the statistics and is where the method should work
    # best), and it directly contradicts stated limit 1, which asserts Pred is a LOWER
    # BOUND on feed-up -- a lower bound cannot exceed the observation.
    sub("spec section 8: hard transport failure test  (Pred_SB > Obs_SB)")
    sb_fail = []
    for k, v in sorted(T7.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])):
        if k[3] != "SB":
            continue
        over = v["pred"] > v["obs"]
        zz = (v["obs"] - v["pred"]) / max(v["sR"] * v["pred"], 1e-9)
        lab_ = SLICINGS[k[1]]["labels"][k[2]]
        print(f"  {k[0]:6s} {k[1]:7s} {lab_:17s} Obs_SB={v['obs']:8.0f} "
              f"Pred_SB={v['pred']:9.1f}  R={v['R']:+7.3f} ({zz:+6.1f} sigma_stat)  "
              + ("<< HARD TRANSPORT FAILURE" if over else "ok"))
        if over:
            sb_fail.append((f"{k[0]}/{k[1]}/{lab_}", v["R"], zz))
    check("spec section 8: no hard transport failure (Pred_SB <= Obs_SB in every slice)",
          not sb_fail,
          f"{len(sb_fail)}/18 SB slices FAIL: " + "; ".join(
              f"{n}: R={r:+.3f} ({z_:+.1f}sd)" for n, r, z_ in sb_fail) if sb_fail
          else "18 SB slices", kind="method")
    RES["T7_sb_transport_failures"] = [{"slice": n, "R": r, "z_stat": z_}
                                       for n, r, z_ in sb_fail]

    sub("T6  fold budget: per-source-row contribution to Pred(TAIL), data, CR")
    for slicing in SLICINGS:
        for si, lab in enumerate(SLICINGS[slicing]["labels"]):
            d = FOLD[("CR", slicing, si)]["data"]
            print(f"  {slicing:7s} {lab:17s} L_src={['%.0f'%x for x in d['L_src']]}  "
                  f"Ntrue3={['%.0f'%x for x in d['Ntrue']]}  "
                  f"contrib={['%.2f'%x for x in d['row_contrib']]}  "
                  f"pred i0/i1/i3={d['pred_tl'][0]:.2f}/{d['pred_tl'][1]:.2f}/{d['pred_tl'][3]:.2f}  "
                  f"predB={d['predB_tl']:.2f}")
    conv_bad = []
    for k, F in FOLD.items():
        for tag in ("data", "mc"):
            for pk in ("pred_sb", "pred_tl"):
                p3, p2 = F[tag][pk][3], F[tag][pk][2]
                if p3 > 0 and abs(p3 - p2) / p3 > 0.01:
                    conv_bad.append((k, tag, pk, abs(p3 - p2) / p3))
    cdist = [abs(v["sig_circ_conv"]) / max(v["pred"], 1e-9) for v in T7.values()] or [0.0]
    check("circularity iteration converged: |Pred(3)-Pred(2)|/Pred < 1 % everywhere",
          not conv_bad,
          (f"{len(conv_bad)} violations, worst {max([x[3] for x in conv_bad]):.4f}; "
           f"|Pred(conv)-Pred(3)|/Pred spans [{min(cdist):.4f}, {max(cdist):.4f}], "
           f"folded into sigma_circ (conservative)")
          if conv_bad else "72 folds", kind="method")
    _tags = ("data", "mc", "data_beforeM")
    _nits = [F[t]["n_conv"] for F in FOLD.values() for t in _tags]
    _cap = [(k, t) for k, F in FOLD.items() for t in _tags if not F[t]["conv_ok"]]
    _mf = [F[t]["conv_min_frac"] for F in FOLD.values() for t in _tags]
    print(f"  EM fixed point (state tol {CONV_TOL:.0e}, {CONV_STABLE} consecutive "
          f"iterations, cap {CONV_ITERS}): iterations to convergence span "
          f"[{min(_nits)}, {max(_nits)}] over {len(_nits)} folds; {len(_cap)} hit the cap "
          f"without converging.  Smallest source-bin share of N_true(conv) spans "
          f"[{min(_mf):.2e}, {max(_mf):.2e}] -- where it is ~0 the limit is a DEGENERATE "
          f"(boundary) solution with a source bin driven to zero, and is reported as the "
          f"limit reached, not as an interior fixed point.")
    naive_hi = sum(1 for k, F in FOLD.items() if F["data"]["pred_tl"][0] > F["data"]["pred_tl"][3])
    ndev = [(F["data"]["pred_tl"][0] / max(F["data"]["pred_tl"][3], 1e-9) - 1)
            for F in FOLD.values()]
    check("naive fold (iteration 0) is biased HIGH relative to the corrected fold, as stated",
          naive_hi >= 0.8 * len(FOLD),
          f"{naive_hi}/{len(FOLD)} folds high; Pred(0)/Pred(3)-1 spans "
          f"[{min(ndev):+.3f}, {max(ndev):+.3f}] -- the sign is NOT universal",
          kind="expectation")
    a, b = analytic_pieces(D, CH_CR, "nHits", 0)
    print(f"\nanalytic first-order variance pieces, CR nHits[0,20) TAIL: "
          f"sum T^2 Var(L) = {a:.3f}; sum L^2 Var(T) = {b:.3f}  "
          f"(-> the {'transfer' if b > a else 'input-spectrum'} term dominates)")

    sub("self-check: MC self-closure and the Before controls")
    mcbad = [(k, v["R_mc"], v["sR_mc"], v["neff_mc"]) for k, v in T7.items()
             if abs(v["R_mc"]) > 0.3 and abs(v["R_mc"]) > 2 * max(v["sR_mc"], 1e-9)]
    check("MC self-closure |R_MC| <= 0.3 or not significant at 2 sigma (method bias)",
          not mcbad, f"{len(mcbad)} flagged: "
          + "; ".join(f"{k[0]}/{k[1]}/{k[2]}/{k[3]} R={r:+.2f}+-{s:.2f} Neff={n:.0f}"
                      for k, r, s, n in mcbad[:6]) if mcbad else "36 entries")
    print("  POWER: per-slice MC closure rests on N_eff ~ 10-20 for Obs_MC, i.e. a "
          "+-30-50 % test; |1-C| is NOT a precision.")
    bd = [T7[("Before", "nHits", i, "TAIL")]["R"] for i in range(3)]
    check("control: Before-DATA does NOT close (the un-vetoed cosmic excess is present)",
          all(x > 1.0 for x in bd), f"R(Before,nHits slices) = {['%.2f'%x for x in bd]}")
    bm = [T7[("Before", "nHits", i, "TAIL")]["R_mc"] for i in range(3)]
    bms = [T7[("Before", "nHits", i, "TAIL")]["sR_mc"] for i in range(3)]
    check("control: Before-MC closes (no cosmics in MC)",
          all(abs(x) < max(3 * s, 0.5) for x, s in zip(bm, bms)),
          f"R_MC = {['%.2f+-%.2f'%(x,s) for x,s in zip(bm,bms)]}")

    # ------------------------------------------------------------------ veto independence
    hr("T4b  VETO-INDEPENDENCE OF THE TRANSFER MATRIX (T_CR/T_Before, T_VR/T_Before)")
    vi = {}
    for slicing in SLICINGS:
        for si, lab in enumerate(SLICINGS[slicing]["labels"]):
            row = []
            Mb, Mbv = get_M8(D, CH_BEFORE, slicing, si)
            for cn, ch in (("CR", CH_CR), ("VR", CH_VR)):
                Mc, Mcv = get_M8(D, ch, slicing, si)
                # DEVIATION FROM SPEC (documented): the spec asks for T_c / T_Before with
                # "binomial-subset errors".  The vetoed channel is a SUBSET of Before, so
                # the subset-difference variance Var(T_c) - Var(T_Before) goes NEGATIVE
                # numerically whenever the two rates nearly coincide; clamping it produced
                # chi2 ~ 1e6 out of rates that agree to a per cent.  The clean and
                # equivalent test is on the DISJOINT complement X = Before \ channel:
                # T_c and T_X are statistically independent, so the pull (T_c - T_X)/sd is
                # well defined.  The DISPLAYED ratio stays T_c/T_Before, as the spec asks.
                chi2 = 0.0
                ndf = 0
                pulls = []
                rr = []
                for i in S_IDX:
                    nc, dc = Mc[i, TAIL_IDX].sum(), Mc[i].sum()
                    nb, db = Mb[i, TAIL_IDX].sum(), Mb[i].sum()
                    nx, dx = nb - nc, db - dc
                    vnc, vdc = Mcv[i, TAIL_IDX].sum(), Mcv[i].sum()
                    vnx = max(Mbv[i, TAIL_IDX].sum() - vnc, 0.0)
                    vdx = max(Mbv[i].sum() - vdc, 0.0)
                    tb = nb / max(db, 1e-9)
                    tc = nc / max(dc, 1e-9)
                    tx = nx / max(dx, 1e-9)
                    ec = wbinom_err(nc, vnc, dc, vdc)
                    ex = wbinom_err(max(nx, 0.0), vnx, max(dx, 1e-9), vdx)
                    r = tc / tb if tb > 0 else float("nan")
                    sd = math.hypot(ec, ex)
                    rr.append((r, ec / max(tb, 1e-9)))
                    if sd > 0 and dx > 0 and np.isfinite(tc) and np.isfinite(tx):
                        pull = (tc - tx) / sd
                        pulls.append(pull)
                        chi2 += pull ** 2
                        ndf += 1
                mx = max((abs(p) for p in pulls), default=0.0)
                row.append((cn, chi2 / max(ndf, 1), rr, mx, ndf))
                print(f"  {slicing:7s} {lab:17s} T_{cn}/T_Before per source row: "
                      + " ".join(f"{r:.3f}+-{s:.3f}" for r, s in rr)
                      + f"   complement chi2/{ndf} = {chi2/max(ndf,1):.2f}"
                      f"  max|pull| = {mx:.2f}")
            vi[f"{slicing}/{lab}"] = [{"chan": c, "chi2_per_ndf": c2, "max_pull": mx,
                                       "ndf": nd,
                                       "ratios": [[float(a), float(b)] for a, b in rr]}
                                      for c, c2, rr, mx, nd in row]
    RES["veto_independence"] = vi
    worst = max((r["max_pull"] for v in vi.values() for r in v), default=0.0)
    worstc = max((r["chi2_per_ndf"] for v in vi.values() for r in v), default=0.0)
    check("transfer matrix is veto-independent (no source row moves > 3 sigma on the "
          "disjoint-complement test)", worst < 3.0,
          f"max |pull| = {worst:.2f}, worst chi2/ndf = {worstc:.2f}", kind="expectation")

    # ------------------------------------------------------------------ scan
    hr("T8  SENSITIVITY SCAN (R1: a closure that closes in only one window is an artefact)")
    print("CR, nHits [0,20) slice.  source edge x target x iteration x matrix channel.")
    print(f"{'srcEdge':8s} {'target':8s} {'iter':5s} {'matrix':8s} {'Obs':>7s} {'Pred':>9s} "
          f"{'R':>8s} {'R_conv':>8s} {'R_cons':>8s} {'R_MC':>8s} {'D':>7s} {'A_low':>7s}")
    SCAN = []
    src_sets = {60: [0, 1, 2], 80: [0, 1, 2, 3], 100: [0, 1, 2, 3, 4]}
    tgt_sets = {"ge200": TAIL_IDX, "ge300": TAIL300_IDX}
    # ---- the scan must use the SAME systematic budget as T7 (fix round, finding 4) ----
    # The previous version used sigma_circ = |Pred(3)-Pred(1)| only, dropping the widened
    # convergence term, sigma_bin and sigma_chan, so its R_cons was ~+0.11 higher than
    # T7's for the identical configuration -- and condition (e) was evaluated on THAT,
    # which made it easier to pass.  sigma_bin has no per-configuration definition (it is
    # the spread of the three inclusive predictions at the nominal window), so it is
    # carried across as a RELATIVE quantity, sigma_bin/Pred from T7's central entry.
    _t7c = T7[("CR", "nHits", 0, "TAIL")]
    SIGBIN_REL = _t7c["sig_bin"] / max(_t7c["pred"], 1e-9)
    print(f"sigma_bin carried across the scan as a relative term, "
          f"sigma_bin/Pred = {SIGBIN_REL:.4f} (from the T7 central entry); "
          f"sigma_circ and sigma_chan are recomputed at every scan point.")

    def scan_one(ST, s, wf, mchan, Si, Ti):
        Mv, Mw = get_M8(ST, mchan or CH_CR, "nHits", s)
        Lv, Lw = get_L(ST, CH_CR, "nHits", s)
        r = fold_core(Mv, Lv, Si, SB_IDX, Ti, 3, conv_iters=CONV_ITERS)
        Mt = gamma_draw(rng, Mv, Mw, args.scan_toys, 1.0 if wf else wbar_of(Mv, Mw))
        Lt = gamma_draw(rng, Lv, Lw, args.scan_toys, 1.0 if wf else wbar_of(Lv, Lw))
        rt = fold_core(Mt, Lt, Si, SB_IDX, Ti, 3)
        return {"ptl": {k: float(v[1]) for k, v in r["preds"].items()},
                "psb": {k: float(v[0]) for k, v in r["preds"].items()},
                "toy": {k: v[1] for k, v in rt["preds"].items()},
                "obs": float(Lv[Ti].sum()), "obs_sb": float(Lv[SB_IDX].sum()),
                "T_from_SB": float(Mv[SB_IDX[0], Ti].sum()
                                   / max(Mv[SB_IDX[0]].sum(), 1e-300)),
                "Lbulk": float(Lv[Si].sum())}

    for se, Si in src_sets.items():
        for tn, Ti in tgt_sets.items():
            bymat = {}
            for mc_ch_name, mch in (("same", None), ("Before", CH_BEFORE)):
                bymat[mc_ch_name] = {"dd": scan_one(D, 0, True, mch, Si, Ti),
                                     "mm": scan_one(M, 0, False, mch, Si, Ti),
                                     "dh": scan_one(D, [1, 2], True, mch, Si, Ti),
                                     "mh": scan_one(M, [1, 2], False, mch, Si, Ti)}
            for it in (0, 1, 3):
                for mc_ch_name in ("same", "Before"):
                    B = bymat[mc_ch_name]
                    other = bymat["Before" if mc_ch_name == "same" else "same"]
                    dd, mm, dh, mh = B["dd"], B["mm"], B["dh"], B["mh"]
                    pred = dd["ptl"][it]
                    predsb = dd["psb"][it]
                    obs = dd["obs"]
                    R = obs / pred - 1 if pred > 0 else float("nan")
                    RM = mm["obs"] / mm["ptl"][it] - 1 if mm["ptl"][it] > 0 else float("nan")
                    Gmc = mm["obs"] - mm["ptl"][it]
                    sG = float(np.std(mm["obs"] - mm["toy"][it]))
                    kap = dd["Lbulk"] / max(mm["Lbulk"], 1e-9)
                    Gup = kap * max(Gmc + sG, 0.0)
                    pfam = [dd["ptl"][0], dd["ptl"][1], dd["ptl"][3], dd["ptl"]["conv"]]
                    sc_up = max(max(pfam) - pred, 0.0)
                    sc_dn = max(pred - min(pfam), 0.0)
                    s_chan = abs(pred - other["dd"]["ptl"][it])
                    s_bin = SIGBIN_REL * pred
                    ff = dd["T_from_SB"] * max(dd["obs_sb"] - predsb, 0.0)
                    Pred_up = pred + math.sqrt(sc_up ** 2 + s_bin ** 2 + s_chan ** 2) + ff
                    Rc = (obs - Pred_up - Gup) / pred if pred > 0 else float("nan")
                    R_cv = obs / dd["ptl"]["conv"] - 1 if dd["ptl"]["conv"] > 0 else float("nan")
                    Dl = (obs / max(pred, 1e-9)) / max(mm["obs"] / max(mm["ptl"][it], 1e-9), 1e-9)
                    Dh = (dh["obs"] / max(dh["ptl"][it], 1e-9)) / \
                         max(mh["obs"] / max(mh["ptl"][it], 1e-9), 1e-9)
                    Al = Dl / Dh if Dh > 0 else float("nan")
                    SCAN.append({"src_edge": se, "target": tn, "iter": it,
                                 "matrix": mc_ch_name, "obs": obs, "pred": pred,
                                 "pred_conv": dd["ptl"]["conv"],
                                 "R": R, "R_conv": R_cv, "R_cons": Rc, "R_mc": RM,
                                 "sig_circ_up": sc_up, "sig_circ_dn": sc_dn,
                                 "sig_bin": s_bin, "sig_chan": s_chan,
                                 "Pred_up": Pred_up, "Gup": Gup,
                                 "D": Dl, "A_low": Al})
                    print(f"{se:8d} {tn:8s} {it:5d} {mc_ch_name:8s} {obs:7.0f} "
                          f"{pred:9.1f} {R:8.3f} {R_cv:8.3f} {Rc:8.3f} {RM:8.3f} "
                          f"{Dl:7.3f} {Al:7.3f}")
    RES["T8_scan"] = SCAN
    sub_scan = [s for s in SCAN if s["iter"] in (1, 3)]
    signs = set(np.sign(s["R_cons"]) for s in sub_scan)
    check("R_cons sign is STABLE across the scan (iteration in {1,3}, srcEdge in {60,80,100}, "
          "target in {>=200,>=300}, matrix in {same,Before})",
          len(signs) == 1, f"signs seen: {sorted(signs)}; "
          f"R_cons range [{min(s['R_cons'] for s in sub_scan):.3f}, "
          f"{max(s['R_cons'] for s in sub_scan):.3f}]")
    sc_at_central = [s for s in SCAN if s["src_edge"] == 100 and s["target"] == "ge200"
                     and s["iter"] == 3 and s["matrix"] == "same"]
    if sc_at_central:
        s0 = sc_at_central[0]
        # R must match to machine precision (it is deterministic).  R_cons carries G_up,
        # whose sigma_G comes from the toy ensemble, and the scan uses fewer toys than
        # T7, so a tolerance of 5e-3 on R_cons is the toy noise, not a discrepancy.
        agree = (abs(s0["R"] - _t7c["R"]) < 1e-9
                 and abs(s0["R_cons"] - _t7c["R_cons"]) < 5e-3)
        check("the scan reproduces T7 at the central configuration "
              "(src<100, >=200, iteration 3, same-channel matrix; R exactly, R_cons to "
              "the toy noise of sigma_G)", agree,
              f"scan R={s0['R']:+.4f} R_cons={s0['R_cons']:+.4f} vs "
              f"T7 R={_t7c['R']:+.4f} R_cons={_t7c['R_cons']:+.4f}")
    print(f"  R range over the whole scan: [{min(s['R'] for s in SCAN):.3f}, "
          f"{max(s['R'] for s in SCAN):.3f}]; R_conv range "
          f"[{min(s['R_conv'] for s in SCAN):.3f}, {max(s['R_conv'] for s in SCAN):.3f}]"
          f"  (R1: a single number here would be an artefact)")

    # ------------------------------------------------------------------ flatness / kappa
    hr("T9  WITHIN-SLICE FLATNESS (R6) AND THE POPULATION-MISMATCH FACTOR kappa_pop")
    T9 = {}
    for slicing, cfg in SLICINGS.items():
        rh = cfg["resp_hist"]
        redges = D.axes[rh][0]
        w2r = resp_weights(redges, 2.0)
        w5r = resp_weights(redges, 5.0)
        dv, dvv = D.get(rh, CH_BEFORE)
        mv, mvv = M.get(rh, CH_BEFORE)
        nfine = dv.shape[1] - 2
        p2d = np.array([float((w2r * dv[:, b + 1]).sum()) for b in range(nfine)])
        n_d = np.array([float(dv[:, b + 1].sum()) for b in range(nfine)])
        p2m = np.array([float((w2r * mv[:, b + 1]).sum()) for b in range(nfine)])
        p2mv = np.array([float((w2r ** 2 * mvv[:, b + 1]).sum()) for b in range(nfine)])
        n_m = np.array([float(mv[:, b + 1].sum()) for b in range(nfine)])
        n_mv = np.array([float(mvv[:, b + 1].sum()) for b in range(nfine)])
        p5d = np.array([float((w5r * dv[:, b + 1]).sum()) for b in range(nfine)])
        p5m = np.array([float((w5r * mv[:, b + 1]).sum()) for b in range(nfine)])
        p5mv = np.array([float((w5r ** 2 * mvv[:, b + 1]).sum()) for b in range(nfine)])
        fd = np.where(n_d > 0, p2d / np.maximum(n_d, 1e-9), 0.0)
        fm = np.where(n_m > 0, p2m / np.maximum(n_m, 1e-9), 0.0)
        efd = np.array([wbinom_err(p2d[b], p2d[b], n_d[b], n_d[b]) for b in range(nfine)])
        efm = np.array([wbinom_err(p2m[b], p2mv[b], n_m[b], n_mv[b]) for b in range(nfine)])
        rat = np.where(fm > 0, fd / np.maximum(fm, 1e-12), np.nan)
        erat = rat * np.sqrt(np.where(fd > 0, (efd / np.maximum(fd, 1e-12)) ** 2, 0)
                             + np.where(fm > 0, (efm / np.maximum(fm, 1e-12)) ** 2, 0))
        fd5 = np.where(n_d > 0, p5d / np.maximum(n_d, 1e-9), 0.0)
        fm5 = np.where(n_m > 0, p5m / np.maximum(n_m, 1e-9), 0.0)
        rat5 = np.where(fm5 > 0, fd5 / np.maximum(fm5, 1e-12), np.nan)
        # LJ vs matched fine-bin weights, for kappa_pop
        ljv, _ = D.get(cfg["lj_hist"], CH_CR)
        mmL = merge_map(D.axes[cfg["lj_hist"]][0], G8)
        lj8 = merge_axis(ljv, 0, mmL)
        lj_bulk_fine = lj8[:5].sum(axis=0)[1:-1]          # inner fine slice bins
        matv, _ = D.get(cfg["mat_hist"] if False else rh, CH_CR)
        mat_fine = matv.sum(axis=0)[1:-1]
        for si, lab in enumerate(cfg["labels"]):
            bins = [b - 1 for b in cfg["resp_sl"][si] if 1 <= b <= nfine]
            r = rat[bins]
            e = erat[bins]
            good = np.isfinite(r) & (e > 0)
            wmean = float(np.sum(r[good] / e[good] ** 2) / max(np.sum(1 / e[good] ** 2), 1e-30))
            chi2 = float(np.sum(((r[good] - wmean) / e[good]) ** 2))
            # two-way binning
            half = len(bins) // 2
            def grp(idx):
                a1 = p2d[idx].sum(); a2 = n_d[idx].sum()
                b1 = p2m[idx].sum(); b2 = n_m[idx].sum()
                b1v = p2mv[idx].sum(); b2v = n_mv[idx].sum()
                fdx = a1 / max(a2, 1e-9); fmx = b1 / max(b2, 1e-9)
                edx = wbinom_err(a1, a1, a2, a2); emx = wbinom_err(b1, b1v, b2, b2v)
                rr = fdx / max(fmx, 1e-12)
                return rr, rr * math.hypot(edx / max(fdx, 1e-12), emx / max(fmx, 1e-12))
            g1 = grp(np.array(bins[:half])) if half else (np.nan, np.nan)
            g2 = grp(np.array(bins[half:])) if half else (np.nan, np.nan)
            chi2_2 = ((g1[0] - g2[0]) ** 2 / max(g1[1] ** 2 + g2[1] ** 2, 1e-30)) \
                if np.isfinite(g1[0]) and np.isfinite(g2[0]) else float("nan")
            wl = lj_bulk_fine[bins]
            wm_ = mat_fine[bins]
            wl = wl / max(wl.sum(), 1e-9)
            wm_ = wm_ / max(wm_.sum(), 1e-9)
            pb = fd[bins]
            kap = float(np.sum(pb * wl) / max(np.sum(pb * wm_), 1e-12))
            centers = 0.5 * (D.axes[rh][1][:-1] + D.axes[rh][1][1:])[bins]
            mean_lj = float(np.sum(wl * centers))
            mean_mt = float(np.sum(wm_ * centers))
            T9[f"{slicing}/{lab}"] = {
                "n_fine_bins": len(bins), "ratio_min": float(np.nanmin(r)),
                "ratio_max": float(np.nanmax(r)), "wmean": wmean,
                "chi2_flat_fine": chi2, "ndf_fine": int(good.sum() - 1),
                "chi2_flat_2way": float(chi2_2),
                "ratio_r5_min": float(np.nanmin(rat5[bins])),
                "ratio_r5_max": float(np.nanmax(rat5[bins])),
                "kappa_pop": kap, "mean_fine_LJ": mean_lj, "mean_fine_matched": mean_mt,
                "per_bin_ratio": [float(x) for x in r]}
            print(f"  {slicing:7s} {lab:17s} P(r>2) d/m over {len(bins):2d} fine bins: "
                  f"[{np.nanmin(r):.2f},{np.nanmax(r):.2f}] wmean={wmean:.2f} "
                  f"chi2/{int(good.sum()-1)}={chi2:.1f} (2-way chi2/1={chi2_2:.1f}) | "
                  f"P(r>5) [{np.nanmin(rat5[bins]):.2f},{np.nanmax(rat5[bins]):.2f}] | "
                  f"kappa_pop={kap:.3f} <fine>_LJ={mean_lj:.2f} <fine>_matched={mean_mt:.2f}")
    RES["T9_flatness"] = T9
    nonflat = [k for k, v in T9.items() if v["chi2_flat_fine"] > 3 * max(v["ndf_fine"], 1)]
    check("within-slice flatness of the data/MC transfer ratio (R6): reported, and the "
          "endcap slice is expected to FAIL", True,
          f"slices failing flatness at chi2 > 3*ndf: {nonflat}")
    kbig = {k: v["kappa_pop"] for k, v in T9.items() if abs(v["kappa_pop"] - 1) > 0.2}
    check("kappa_pop within 20 % of 1 (else the range must be quoted, and max(kappa,1) "
          "enters Pred_up)", not kbig, f"|kappa-1|>0.2 in: "
          f"{ {k: round(v,3) for k,v in kbig.items()} }" if kbig else "6 slices")

    # tail-vs-bulk shape (R7)
    sub("R7  tail-vs-bulk data/MC shape ratio of the slice variable (unit-normalised)")
    for slicing, cfg in SLICINGS.items():
        for cn, ch in FOLD_CHANNELS:
            dv, _ = D.get(cfg["lj_hist"], ch)
            mv, mvv = M.get(cfg["lj_hist"], ch)
            mmL = merge_map(D.axes[cfg["lj_hist"]][0], G8)
            d8 = merge_axis(dv, 0, mmL)
            m8 = merge_axis(mv, 0, mmL)
            m8v = merge_axis(mvv, 0, mmL)
            for reg_name, ridx in (("bulk<100", [0, 1, 2, 3, 4]), ("tail>=200", TAIL_IDX)):
                dd = np.array([d8[ridx][:, g].sum() for g in cfg["lj_sl"]])
                mm = np.array([m8[ridx][:, g].sum() for g in cfg["lj_sl"]])
                mmv = np.array([m8v[ridx][:, g].sum() for g in cfg["lj_sl"]])
                if dd.sum() <= 0 or mm.sum() <= 0:
                    continue
                fdv = dd / dd.sum()
                fmv = mm / mm.sum()
                vfd = fdv * (1 - fdv) / dd.sum()
                vfm = np.array([((1 - fmv[i]) ** 2 * mmv[i]
                                 + fmv[i] ** 2 * max(mmv.sum() - mmv[i], 0)) / mm.sum() ** 2
                                for i in range(len(mm))])
                chi2 = float(np.sum((fdv - fmv) ** 2 / np.maximum(vfd + vfm, 1e-30)))
                print(f"  {slicing:7s} {cn:6s} {reg_name:9s} shape d/m = "
                      + " ".join(f"{a/max(b,1e-9):.2f}" for a, b in zip(fdv, fmv))
                      + f"   chi2/{len(dd)-1} = {chi2:.1f}")

    # ------------------------------------------------------------------ stage 2 / floor
    hr("T10  STAGE-2 DECOMPOSITION AND THE BIAS-FREE 'NECESSARILY MIGRATED' FLOOR")
    print("stage 1 = sources below 100 GeV; stage 2 = T(TAIL|SB) x (Obs_SB - Pred_SB)_+, "
          "i.e. the genuine [100,200) content fed forward; rem = the remainder.")
    T10 = {}
    for slicing in SLICINGS:
        for si, lab in enumerate(SLICINGS[slicing]["labels"]):
            for cn in ("CR", "VR"):
                for tag, ST in (("data", D), ("mc", M)):
                    F = FOLD[(cn, slicing, si)][tag]
                    p1 = F["pred_tl"][3]
                    p2 = F["T_TAIL_from_SB"] * max(F["obs_sb"] - F["pred_sb"][3], 0.0)
                    obs = F["obs_tl"]
                    T10[f"{cn}/{slicing}/{lab}/{tag}"] = {
                        "pred1_over_obs": p1 / max(obs, 1e-9),
                        "pred2_over_obs": p2 / max(obs, 1e-9),
                        "rem_over_obs": (obs - p1 - p2) / max(obs, 1e-9)}
                    if cn == "CR":
                        print(f"  {slicing:7s} {lab:17s} {tag:4s} obs={obs:8.1f} "
                              f"pred1/obs={p1/max(obs,1e-9):6.3f} "
                              f"pred2/obs={p2/max(obs,1e-9):6.3f} "
                              f"rem/obs={(obs-p1-p2)/max(obs,1e-9):+6.3f}")
    floors = {}
    ptPf_e = D.axes["dsaPfMatch_ptPf_ptDsa"][0]
    ptDsa_e = D.axes["dsaPfMatch_ptPf_ptDsa"][1]
    for slicing, cfg in SLICINGS.items():
        for si, lab in enumerate(cfg["labels"]):
            v, _ = D.get(cfg["mat_hist"], CH_CR)
            v = pick_axis(v, 2, cfg["mat_sl"][si])
            nb = len(ptPf_e) - 1
            tot = v.sum()
            sel = 0.0
            for i in range(nb):
                hi = ptPf_e[i + 1]
                for j in range(nb):
                    if ptDsa_e[j] >= 2.0 * hi:
                        sel += v[i + 1, j + 1]
                if ptDsa_e[-1] >= 2.0 * hi:
                    sel += v[i + 1, -1]
            floors[f"{slicing}/{lab}"] = {"necessarily_r2_frac": float(sel / max(tot, 1e-9)),
                                          "counts": float(sel)}
            print(f"  one-sided bias-free floor (cells with ptDsa_lo >= 2*ptPf_hi), "
                  f"CR {slicing} {lab}: {sel:.0f} of {tot:.0f} = {sel/max(tot,1e-9):.4f}")
    RES["T10_stage2"] = T10
    RES["T10_floor"] = floors

    # ------------------------------------------------------------------ era stability
    hr("T12  ERA STABILITY  R_data(TAIL) per era")
    T12 = {}
    print(f"{'era':5s} {'chunks':>8s} " + " ".join(f"{l[:14]:>16s}"
                                                   for l in SLICINGS["nHits"]["labels"]))
    for e in "ABCD":
        if not ERA[e].n_files:
            continue
        row = []
        for si in range(3):
            Mv, Mw = get_M8(ERA[e], CH_CR, "nHits", si)
            Lv, Lw = get_L(ERA[e], CH_CR, "nHits", si)
            r = fold_core(Mv, Lv, S_IDX, SB_IDX, TAIL_IDX, 3)
            pred = float(r["preds"][3][1])
            obs = float(Lv[TAIL_IDX].sum())
            R = obs / pred - 1 if pred > 0 else float("nan")
            sR = math.sqrt(max(obs, 1)) / pred if pred > 0 else float("nan")
            row.append((R, sR, obs, pred))
        T12[e] = [{"R": float(a), "sR": float(b), "obs": c, "pred": d} for a, b, c, d in row]
        flag = (f" (PARTIAL {era_chunks_have.get(e,0)}/{ERA_CHUNKS_EXPECTED[e]})"
                if era_chunks_have.get(e, 0) < ERA_CHUNKS_EXPECTED[e] else "")
        print(f"{e:5s} {era_chunks_have.get(e,0):3d}/{ERA_CHUNKS_EXPECTED[e]:<4d} "
              + " ".join(f"{a:+8.3f}+-{b:5.3f}" for a, b, _, _ in row) + flag)
    RES["T12_era"] = T12
    # NB: the chi2 summary lives in its OWN dict.  Writing it back into T12 would put a
    # non-era key into the dict the next slice iterates over (it did: KeyError on slice 1).
    T12chi2 = {}
    for si in range(3):
        vals = [(T12[e][si]["R"], T12[e][si]["sR"]) for e in T12
                if np.isfinite(T12[e][si]["R"]) and T12[e][si]["sR"] > 0]
        if len(vals) > 1:
            w = np.array([1 / max(s, 1e-9) ** 2 for _, s in vals])
            mu = float(np.sum([v * ww for (v, _), ww in zip(vals, w)]) / w.sum())
            chi2 = float(np.sum([((v - mu) / s) ** 2 for v, s in vals]))
            print(f"  era chi2/{len(vals)-1} for nHits slice {si} = {chi2:.1f} "
                  f"(weighted mean R = {mu:+.3f})")
            T12chi2[f"slice{si}"] = {"chi2": chi2, "ndf": len(vals) - 1, "mean_R": mu}
    RES["T12_era_chi2"] = T12chi2
    check("era stability of R(TAIL) reported for every loaded era, with D flagged partial",
          len(T12) >= 1, f"eras: {[e for e in 'ABCD' if ERA[e].n_files]}")

    # ------------------------------------------------------------------ relPtErr
    hr("T9b  FAKE-MATCH / SELECTION-BIAS EXPOSURE  P(relPtErr < 0.3 | r > 2), data vs MC")
    redges = D.axes["dsaPfMatch_response_relPtErr"][0]
    peedges = D.axes["dsaPfMatch_response_relPtErr"][1]
    w2r = resp_weights(redges, 2.0)
    j03 = edge_index(peedges, 0.3)
    lowe = list(range(1, j03 + 1))
    out_rp = {}
    for cn, ch in FOLD_CHANNELS:
        dv, _ = D.get("dsaPfMatch_response_relPtErr", ch)
        mv, mvv = M.get("dsaPfMatch_response_relPtErr", ch)
        dn = float((w2r[:, None] * dv).sum())
        dl = float((w2r[:, None] * dv[:, lowe]).sum())
        mn = float((w2r[:, None] * mv).sum())
        ml = float((w2r[:, None] * mv[:, lowe]).sum())
        mnv = float((w2r[:, None] ** 2 * mvv).sum())
        mlv = float((w2r[:, None] ** 2 * mvv[:, lowe]).sum())
        pd = dl / max(dn, 1e-9)
        pm = ml / max(mn, 1e-9)
        epd = wbinom_err(dl, dl, dn, dn)
        epm = wbinom_err(ml, mlv, mn, mnv)
        rr = pd / max(pm, 1e-12)
        err = rr * math.hypot(epd / max(pd, 1e-12), epm / max(pm, 1e-12))
        out_rp[cn] = {"data": pd, "data_err": epd, "mc": pm, "mc_err": epm,
                      "ratio": rr, "ratio_err": err, "n_data_r2": dn, "n_mc_r2": mn}
        print(f"  {cn:6s} P(relPtErr<0.3 | r>2): data {pd:.4f}+-{epd:.4f} (n={dn:.0f})  "
              f"MC {pm:.4f}+-{epm:.4f}  ratio {rr:.3f}+-{err:.3f}")
    RES["relPtErr_exposure"] = out_rp
    ljr = "mu_lj_dsaMuon_pt_relPtErr"
    if ljr in D.v:
        mmL = merge_map(D.axes[ljr][0], G8)
        pe_e = D.axes[ljr][1]
        jj = edge_index(pe_e, 0.3)
        for cn, ch in (("CR", CH_CR),):
            dv, _ = D.get(ljr, ch)
            d8 = merge_axis(dv, 0, mmL)
            for nm, idx in (("bulk<100", [0, 1, 2, 3, 4]), ("tail>=200", TAIL_IDX)):
                tot = d8[idx].sum()
                lo = d8[idx][:, 1:jj + 1].sum()
                print(f"  LJ DSA muons (CR, {nm}): P(relPtErr<0.3) = {lo/max(tot,1e-9):.4f} "
                      f"(n={tot:.0f})  -- compare with the matched value above "
                      f"(LJ vs matched sigma(pT)/pT)")

    # ------------------------------------------------------------------ phi test
    hr("T11  PHI TEST on mu_lj_dsaMuon_pt_phi_nHits")
    phi_e = D.axes["mu_lj_dsaMuon_pt_phi_nHits"][1]
    parts = phi_partitions(phi_e)
    print(f"partition primary: V bins {parts['primary'][0].tolist()} "
          f"H bins {parts['primary'][1].tolist()}  (centres within 0.5 rad of +-pi/2 / of 0,+-pi; "
          f"8 diagonal bins excluded)")
    print(f"partition secondary (8/8): V {parts['secondary'][0].tolist()}")
    PT_TAIL_IDX = [3, 4]      # coarse cube pT view: [200,500) + OF
    PT_SB_IDX = [2]
    PT_BULK_IDX = [1]
    T11 = {}
    print("Every fit now carries a GOODNESS OF FIT (fix round, findings 2 and 9).  The "
          "two-template family can only produce phi shapes INSIDE the span of the cosmic "
          "and non-cosmic templates; when the data sit outside it the minimiser rails at "
          "f = 0 and the profile interval is the width of a REFUTED model, not a bound.  "
          "An interval from a fit with gof_p < 0.01 is NOT quoted and the slice is called "
          "INCONCLUSIVE.")
    print("The non-cosmic template (data, same channel, same nHits slice, pT < 100) is now "
          "cross-checked against the MC same-slice bulk AND the MC same-slice tail, as "
          "spec section 10 requires; a tail-vs-bulk phi shape difference is the direct "
          "explanation of a failed fit.")
    for cn, ch in FOLD_CHANNELS:
        cube, cubew = D.get("mu_lj_dsaMuon_pt_phi_nHits", ch)
        mcube, mcubew = M.get("mu_lj_dsaMuon_pt_phi_nHits", ch)
        cos_primary = COSMIC_FOR[ch]
        for si, lab in enumerate(SLICINGS["nHits"]["labels"]):
            nh = sl_cube_nh[si]
            for reg_name, pidx in (("TAIL", PT_TAIL_IDX), ("SB", PT_SB_IDX)):
                obs = cube[np.ix_(pidx, range(cube.shape[1]), nh)].sum(axis=(0, 2))[1:-1]
                nc = cube[np.ix_(PT_BULK_IDX, range(cube.shape[1]), nh)].sum(axis=(0, 2))[1:-1]
                # ---- MC cross-checks of the NON-COSMIC template (spec section 10) -----
                nc_mc = mcube[np.ix_(PT_BULK_IDX, range(mcube.shape[1]), nh)].sum(axis=(0, 2))[1:-1]
                nc_mcw = mcubew[np.ix_(PT_BULK_IDX, range(mcubew.shape[1]), nh)].sum(axis=(0, 2))[1:-1]
                tl_mc = mcube[np.ix_(pidx, range(mcube.shape[1]), nh)].sum(axis=(0, 2))[1:-1]
                tl_mcw = mcubew[np.ix_(pidx, range(mcubew.shape[1]), nh)].sum(axis=(0, 2))[1:-1]
                chi2_nc = float("nan")
                if nc.sum() > 0 and nc_mc.sum() > 0:
                    f_d = nc / nc.sum()
                    v_d = f_d * (1 - f_d) / nc.sum()
                    f_m = nc_mc / nc_mc.sum()
                    v_m = ((1 - f_m) ** 2 * nc_mcw
                           + f_m ** 2 * np.maximum(nc_mcw.sum() - nc_mcw, 0)) / nc_mc.sum() ** 2
                    chi2_nc = float(np.sum((f_d - f_m) ** 2 / np.maximum(v_d + v_m, 1e-30)))
                res = {}
                for tmpl_name, tch in (("cos_spread", CH_COS_SPREAD), ("cos_vtx", CH_COS_VTX)):
                    tc, tcw = D.get("mu_lj_dsaMuon_pt_phi_nHits", tch)
                    cosv = tc[np.ix_(pidx, range(tc.shape[1]), nh)].sum(axis=(0, 2))[1:-1]
                    cosw = tcw[np.ix_(pidx, range(tcw.shape[1]), nh)].sum(axis=(0, 2))[1:-1]
                    fallback = False
                    chi2_pt = float("nan")
                    if cosv.sum() < 50:
                        allpt = list(range(1, tc.shape[0]))
                        c2 = tc[np.ix_(allpt, range(tc.shape[1]), nh)].sum(axis=(0, 2))[1:-1]
                        if cosv.sum() > 0 and c2.sum() > 0:
                            f1 = cosv / cosv.sum()
                            f2 = c2 / c2.sum()
                            v1 = f1 * (1 - f1) / max(cosv.sum(), 1)
                            v2 = f2 * (1 - f2) / max(c2.sum(), 1)
                            chi2_pt = float(np.sum((f1 - f2) ** 2 / np.maximum(v1 + v2, 1e-30)))
                        cosv = c2
                        cosw = tcw[np.ix_(allpt, range(tcw.shape[1]), nh)].sum(axis=(0, 2))[1:-1]
                        fallback = True
                    fit = fit_fcos(obs, cosv, np.maximum(cosw, cosv), nc, nc) \
                        if (obs.sum() > 0 and cosv.sum() > 0 and nc.sum() > 0) else None
                    for pn, (Vi, Hi) in parts.items():
                        Vo, Ho = float(obs[Vi].sum()), float(obs[Hi].sum())
                        Vc, Hc = float(cosv[Vi].sum()), float(cosv[Hi].sum())
                        Vn, Hn = float(nc[Vi].sum()), float(nc[Hi].sum())
                        Vm, Hm = float(nc_mc[Vi].sum()), float(nc_mc[Hi].sum())
                        Vt, Ht = float(tl_mc[Vi].sum()), float(tl_mc[Hi].sum())
                        Vtv, Htv = float(tl_mcw[Vi].sum()), float(tl_mcw[Hi].sum())
                        po = Vo / max(Vo + Ho, 1e-9)
                        pc = Vc / max(Vc + Hc, 1e-9)
                        pn_ = Vn / max(Vn + Hn, 1e-9)
                        pm_bulk = Vm / max(Vm + Hm, 1e-9)
                        pm_tail = Vt / max(Vt + Ht, 1e-9)
                        neff_mtail = ((Vt + Ht) ** 2 / max(Vtv + Htv, 1e-30)) \
                            if (Vtv + Htv) > 0 else 0.0
                        epo = math.sqrt(max(po * (1 - po) / max(Vo + Ho, 1), 0))
                        epc = math.sqrt(max(pc * (1 - pc) / max(Vc + Hc, 1), 0))
                        epn = math.sqrt(max(pn_ * (1 - pn_) / max(Vn + Hn, 1), 0))
                        # how far the observed verticality sits from the non-cosmic
                        # template, in sigma.  NEGATIVE = the tail is LESS vertical than
                        # its own template, i.e. OUTSIDE the mixture family (f < 0).
                        anti = (po - pn_) / math.hypot(epo, epn) if math.hypot(epo, epn) > 0 \
                            else float("nan")
                        den = pc - pn_
                        f = (po - pn_) / den if abs(den) > 1e-6 else float("nan")
                        ef = (math.sqrt(epo ** 2 + epn ** 2 * (1 - f) ** 2
                                        + epc ** 2 * f ** 2) / abs(den)) \
                            if abs(den) > 1e-6 and np.isfinite(f) else float("nan")
                        res[f"{tmpl_name}/{pn}"] = {
                            "V": Vo, "H": Ho, "V_over_H": Vo / max(Ho, 1e-9),
                            "p_obs": po, "p_obs_err": epo, "p_cos": pc, "p_nc": pn_,
                            "p_nc_err": epn,
                            "p_nc_mc_bulk": pm_bulk, "p_nc_mc_tail": pm_tail,
                            "neff_mc_tail": neff_mtail,
                            "sigma_obs_vs_nc": anti,
                            "chi2_nc_data_vs_mc_bulk": chi2_nc,
                            "chi2_cos_template_pT_indep": chi2_pt,
                            "f_solve": f, "f_solve_err": ef,
                            "n_cos_template": float(cosv.sum()),
                            "template_fallback_all_pT": fallback,
                            "fit": fit}
                key = f"{cn}/{lab}/{reg_name}"
                T11[key] = res
                if reg_name == "TAIL":
                    pr = res[f"{'cos_spread' if cos_primary==CH_COS_SPREAD else 'cos_vtx'}/primary"]
                    ft = pr["fit"]
                    print(f"  {cn:6s} {lab:17s} TAIL  V={pr['V']:6.0f} H={pr['H']:6.0f} "
                          f"V/H={pr['V_over_H']:5.3f}  p_obs={pr['p_obs']:.3f}+-{pr['p_obs_err']:.3f} "
                          f"p_cos={pr['p_cos']:.3f} p_nc={pr['p_nc']:.3f}+-{pr['p_nc_err']:.3f}  "
                          f"(obs-nc)={pr['sigma_obs_vs_nc']:+.2f} sigma  "
                          f"f_solve={pr['f_solve']:+.3f}+-{pr['f_solve_err']:.3f}")
                    if ft:
                        print(f"         fit f={ft['fhat']:.3f} [{ft['lo68']:.3f},{ft['hi68']:.3f}] "
                              f"LB95={ft['LB95']:.3f} UB95={ft['UB95']:.3f}  "
                              f"GOF: dev={ft['gof_dev']:.1f}/{ft['gof_ndf']} p={ft['gof_p']:.2e} "
                              f"{'REJECTED -- interval NOT quotable' if ft['rejected'] else 'ok'}"
                              + ("  [railed at boundary]" if ft["railed_at_boundary"] else ""))
                    else:
                        print("         fit n/a")
                    print(f"         MC cross-check of the non-cosmic template: "
                          f"p_nc(data bulk)={pr['p_nc']:.3f}  p(MC bulk)={pr['p_nc_mc_bulk']:.3f}  "
                          f"p(MC tail)={pr['p_nc_mc_tail']:.3f} (Neff={pr['neff_mc_tail']:.0f})  "
                          f"phi-shape chi2/15 data-vs-MC bulk = {pr['chi2_nc_data_vs_mc_bulk']:.1f}  "
                          f"tmpl n={pr['n_cos_template']:.0f}"
                          + ("  [FALLBACK all-pT template]" if pr["template_fallback_all_pT"] else ""))
    rej = [f"{k}:{tk}" for k, v in T11.items() if k.endswith("/TAIL")
           for tk in ("cos_spread/primary", "cos_vtx/primary")
           if v.get(tk, {}).get("fit") and v[tk]["fit"]["rejected"]]
    check("phi mixture fit describes the data (goodness of fit p >= 0.01) in every "
          "TAIL slice", not rej,
          f"{len(rej)} REJECTED fits: {rej} -- no cosmic-fraction interval may be quoted "
          f"there" if rej else "9 TAIL slices", kind="method")
    RES["T11_phi"] = T11
    # nHits composition of the templates
    for tn, tch in (("cosmic_muons", CH_COS_SPREAD), ("cosmic_muons_vtxChi2", CH_COS_VTX)):
        tc, _ = D.get("mu_lj_dsaMuon_pt_phi_nHits", tch)
        comp = [float(tc[np.ix_(PT_TAIL_IDX, range(tc.shape[1]), g)].sum()) for g in sl_cube_nh]
        tot = sum(comp)
        print(f"  template {tn:22s} tail nHits composition = "
              + " ".join(f"{c/max(tot,1e-9):.3f}" for c in comp) + f"  (n={tot:.0f})")
    bpc = T11.get(f"Before/{SLICINGS['nHits']['labels'][0]}/TAIL", {}).get("cos_spread/primary")
    allb = []
    for si, lab in enumerate(SLICINGS["nHits"]["labels"]):
        r = T11.get(f"Before/{lab}/TAIL", {}).get("cos_spread/primary")
        if r and r["fit"]:
            allb.append(r["fit"]["fhat"])
    check("phi positive control: the Before channel is strongly cosmic-like "
          "(the veto removes ~2100 of 2414 tail muons on 2018C)",
          bool(allb) and max(allb) > 0.5, f"Before f_hat per nHits slice = "
          f"{['%.2f'%x for x in allb]}")

    # ------------------------------------------------------------------ T14
    hr("T14  (SECONDARY) RAW TAIL data/MC PER SLICE, MC SCALED TO L_data BY THE ERA TABLE")
    print(f"This is the ONLY table where the luminosity enters.  {len(mc_present)}/8 MC samples "
          + (f"(missing: {mc_missing}, incl. the highest-pT bin) " if mc_missing else "(all present) ")
          + "bulk coverage of the MC stack is "
          "only ~53 %, so this is continuity with the published 2.4, not a measurement.")
    lr = L_data / LUMI_MC_FB
    T14 = {}
    for slicing, cfg in SLICINGS.items():
        for si, lab in enumerate(cfg["labels"]):
            for cn, ch in FOLD_CHANNELS:
                d = FOLD[(cn, slicing, si)]["data"]
                m = FOLD[(cn, slicing, si)]["mc"]
                mo = m["obs_tl"] * lr
                mv = math.sqrt(m["obs_tl_var"]) * lr
                r = d["obs_tl"] / max(mo, 1e-9)
                er = r * math.hypot(math.sqrt(max(d["obs_tl"], 1)) / max(d["obs_tl"], 1e-9),
                                    mv / max(mo, 1e-9))
                T14[f"{cn}/{slicing}/{lab}"] = {"data": d["obs_tl"], "mc": mo, "mc_err": mv,
                                                "ratio": r, "ratio_err": er,
                                                "neff_mc": neff_of(np.array([m["obs_tl"]]),
                                                                   np.array([m["obs_tl_var"]]))}
                if cn == "CR":
                    print(f"  CR {slicing:7s} {lab:17s} data {d['obs_tl']:7.0f}  "
                          f"MC {mo:8.1f}+-{mv:6.1f} (Neff {m['obs_tl']**2/max(m['obs_tl_var'],1e-9):5.0f})  "
                          f"data/MC {r:6.3f}+-{er:.3f}")
    RES["T14_raw"] = T14

    # ------------------------------------------------------------------ decision
    hr("T13  DECISION READOUT  (pre-registered; PROVISIONAL until 68/68 chunks and 8/8 MC)")
    dec = {}
    ktl = ("CR", "nHits", 0, "TAIL")
    ksb = ("CR", "nHits", 0, "SB")
    tl, sb = T7[ktl], T7[ksb]
    bef = T7[("Before", "nHits", 0, "TAIL")]

    # ---- merged high-nHits slice, DATA *and* MC (fix round, finding 6) ---------------
    # The raw data difference R_low - R_high is confounded by the method's OWN
    # nHits-dependent bias, which the MC self-closure measures directly and which is
    # several times larger than the data effect and of the same sign.  Condition (c) is
    # therefore evaluated as the DOUBLE DIFFERENCE against the MC self-closure; the raw
    # difference and the MC swing are both printed so the confounder is visible.
    def merged_high(ST, wf):
        Mv, Mw = get_M8(ST, CH_CR, "nHits", [1, 2])
        Lv, Lw = get_L(ST, CH_CR, "nHits", [1, 2])
        r = fold_core(Mv, Lv, S_IDX, SB_IDX, TAIL_IDX, 3)
        Mt = gamma_draw(rng, Mv, Mw, args.toys, 1.0 if wf else wbar_of(Mv, Mw))
        Lt = gamma_draw(rng, Lv, Lw, args.toys, 1.0 if wf else wbar_of(Lv, Lw))
        rt = fold_core(Mt, Lt, S_IDX, SB_IDX, TAIL_IDX, 3)
        obs = float(Lv[TAIL_IDX].sum())
        ovar = float(Lw[TAIL_IDX].sum())
        pred = float(r["preds"][3][1])
        R = obs / pred - 1 if pred > 0 else float("nan")
        sR = math.hypot(float(np.std(obs / np.maximum(rt["preds"][3][1], 1e-9) - 1)),
                        math.sqrt(max(ovar, 1.0)) / max(pred, 1e-9))
        return R, sR, obs, pred
    R_h, sR_h, obs_h, pred_h = merged_high(D, True)
    R_hM, sR_hM, obs_hM, pred_hM = merged_high(M, False)

    # ---- PRE-REGISTERED STOP CONDITIONS, NOW GATES (fix round, finding 13) -----------
    # Spec 7.3: a flagged |R_MC| > 0.3 at > 2 sigma "must be explained by the section 4/6
    # variants BEFORE data is read".  Spec 7.4: "If Before-data closes, the method has no
    # power and the study stops there."  Both were informational check() lines; both fire.
    g73_ok = not tl["mc_bias_flag"]
    g74_ok = bef["R"] > 1.0
    validated = g73_ok and g74_ok

    # (a)
    a_ok = (tl["z_cons_gauss"] >= 3.0) and (tl["R_cons"] >= 1.0 / 3.0)
    # (b) tail-specific.  sigma_circ is FULLY CORRELATED between the SB and the TAIL of one
    # slice (spec 6f), so it is propagated by moving both to the same end of their
    # brackets rather than dropped (fix round, finding 11).
    dtoy = tl["toyR"] - sb["toyR"]
    b_stat = math.hypot(float(np.std(dtoy)), math.hypot(tl["sR_obs"], sb["sR_obs"]))
    def _dif(pT_, pS_):
        return (tl["obs"] / max(pT_, 1e-9) - 1) - (sb["obs"] / max(pS_, 1e-9) - 1)
    b_val = _dif(tl["pred"], sb["pred"])
    b_sys = max(abs(_dif(tl["pred_hi"], sb["pred_hi"]) - b_val),
                abs(_dif(tl["pred_lo"], sb["pred_lo"]) - b_val))
    sdiff = math.hypot(b_stat, b_sys)
    b_ok = b_val >= 2 * sdiff
    # (c) nHits-specific, as a double difference against the MC self-closure
    c_raw = tl["R"] - R_h
    c_raw_sd = math.hypot(tl["sR_full"], sR_h)
    c_mc = tl["R_mc"] - R_hM
    c_mc_sd = math.hypot(tl["sR_mc"], sR_hM)
    c_val = c_raw - c_mc
    c_sd = math.hypot(c_raw_sd, c_mc_sd)
    c_ok = c_val >= 2 * c_sd
    # (d) VR
    vr = T7[("VR", "nHits", 0, "TAIL")]
    d_ok = (vr["R_cons"] > 0) and (vr["z_cons_gauss"] >= 1.5)
    # (e) window stability -- now on the SAME R_cons definition as T7 (finding 4)
    e_ok = all(s["R_cons"] > 0 for s in SCAN if s["iter"] in (1, 3))
    refuses = a_ok and b_ok and c_ok and d_ok and e_ok
    method_fail = a_ok and not (b_ok and e_ok)

    sub("method validation gates (spec 7.3, 7.4) -- checked BEFORE the conditions")
    print(f"  7.3  MC self-closure in the slice under test: R_MC = {tl['R_mc']:+.3f} "
          f"+- {tl['sR_mc']:.3f} at N_eff = {tl['neff_mc']:.0f}.  The fold "
          f"{tl['mc_direction']} the MC low-nHits tail by a factor "
          f"{tl['mc_factor']:.2f} (Pred_MC = {tl['pred_mc']:.1f} vs Obs_MC = "
          f"{tl['obs_mc']:.1f}).  -> {'PASS' if g73_ok else 'FAIL'}")
    print(f"  7.4  Before-data must NOT close (the un-vetoed cosmic excess): "
          f"R(Before, nHits[0,20), TAIL) = {bef['R']:+.3f} vs the coded threshold "
          f"R > 1.0  -> {'PASS' if g74_ok else 'FAIL'}")
    print(f"  => method validated in the slice under test: {validated}")

    print(f"\n  slice under test: CR, TAIL, nHits [0,20); central config = source<100, "
          f"iteration 3, same-channel matrix, target >=200")
    print(f"  (a) z_cons >= 3 AND R_cons >= 1/3      : z_cons(gauss)={tl['z_cons_gauss']:.2f} "
          f"(toy-count z {'<=' if tl['z_cons_clamped'] else '='}{tl['z_cons']:+.2f}, "
          f"{tl['n_toys']} toys) R_cons={tl['R_cons']:+.3f}   -> {'PASS' if a_ok else 'FAIL'}")
    print(f"  (b) tail-specific (R_TAIL-R_SB) >= 2sd : {b_val:+.3f} vs 2sd={2*sdiff:.3f} "
          f"(stat {b_stat:.3f} (+) circ {b_sys:.3f})   -> {'PASS' if b_ok else 'FAIL'}")
    print(f"  (c) nHits-specific, DOUBLE DIFFERENCE  : "
          f"(R_lo-R_hi)_data {c_raw:+.3f}+-{c_raw_sd:.3f} MINUS "
          f"(R_lo-R_hi)_MC {c_mc:+.3f}+-{c_mc_sd:.3f} = {c_val:+.3f} vs 2sd={2*c_sd:.3f} "
          f"-> {'PASS' if c_ok else 'FAIL'}")
    print(f"        (raw data-only difference was {c_raw:+.3f}; the METHOD's own "
          f"nHits swing in MC is {c_mc:+.3f}, i.e. {abs(c_mc)/max(abs(c_raw),1e-9):.1f}x the "
          f"data effect and of the SAME sign -- the raw form is not interpretable)")
    print(f"        R_high(data)={R_h:+.3f}+-{sR_h:.3f}  R_high(MC)={R_hM:+.3f}+-{sR_hM:.3f}")
    print(f"  (d) VR corroborates R_cons>0, z>=1.5   : R_cons={vr['R_cons']:+.3f} "
          f"z_cons(gauss)={vr['z_cons_gauss']:.2f}   -> {'PASS' if d_ok else 'FAIL'}")
    print(f"  (e) window-stable R_cons>0 across scan : "
          f"{sum(1 for s in SCAN if s['iter'] in (1,3) and s['R_cons']>0)}"
          f"/{sum(1 for s in SCAN if s['iter'] in (1,3))} points positive "
          f"-> {'PASS' if e_ok else 'FAIL'}")
    print(f"  => 'REFUSES TO CLOSE' = {refuses}"
          + ("   [(a) holds but (b)/(e) fails: METHOD FAILURE in the slice, verdict = "
             "default, ESCALATE TO ANALYST]" if method_fail else ""))

    # ---- the closure, stated three ways (fix round, findings 1, 7, 8, 11) ------------
    sub("the deciding slice stated three ways -- none of them is 'R = X +- stat'")
    print(f"  1. AS MEASURED, iteration 3 (pre-registered): R = {tl['R']:+.3f} "
          f"+- {tl['sR']:.3f} (stat) +- {tl['sR_syst']:.3f} (syst) "
          f"= {tl['R']:+.3f} +- {tl['sR_full']:.3f} (full)")
    print(f"  2. ITERATION BRACKET: Pred(0 naive)={tl['pred_i0']:.1f} "
          f"Pred(1)={tl['pred_i1']:.1f} Pred(3)={tl['pred']:.1f} "
          f"Pred(conv, {tl['n_conv']} its)={tl['pred_conv']:.1f}  ->  "
          f"R spans [{min(tl['R_bracket']):+.3f}, {max(tl['R_bracket']):+.3f}]; "
          f"R(conv) = {tl['R_conv']:+.3f}.  "
          + ("THE SIGN OF R FLIPS INSIDE THE BRACKET: the headline closure is an artefact "
             "of where the recursion is stopped."
             if min(tl["R_bracket"]) * max(tl["R_bracket"]) < 0 else
             "The sign of R is stable across the bracket."))
    print(f"      Pred(conv) is the EM FIXED POINT, not a truncation: reached="
          f"{tl['conv_ok']} after {tl['n_conv']} iterations "
          f"(state tol {CONV_TOL:.0e}, "
          f"{CONV_STABLE} consecutive); smallest source-bin share of N_true(conv) = "
          f"{tl['conv_min_frac']:.2e}"
          + ("  -- a source bin is driven to ZERO, so the limit is a DEGENERATE "
             "(boundary) ML solution; it is quoted as the limit, not as an interior "
             "fixed point." if tl["conv_min_frac"] < 1e-3 else "."))
    print(f"  3. AFTER THE METHOD'S OWN MC BIAS: the fold {tl['mc_direction']} the MC tail "
          f"by {tl['mc_factor']:.2f}x, so the bias-transported closure is "
          f"R_bc = D - 1 = {tl['R_bc']:+.3f} +- {tl['sR_bc']:.3f}.  "
          f"|R_MC| = {abs(tl['R_mc']):.3f} is {abs(tl['R_mc'])/max(tl['sR_full'],1e-9):.0f}x "
          f"the full quoted precision, so THIS SLICE IS BIAS-LIMITED, NOT STAT-LIMITED.")
    print(f"  D (= C_data/C_MC) = {tl['D']:.3f}.  D was demoted from a trigger on the "
          f"grounds that its null departs from 1 whenever MC carries a genuine component "
          f"(R_MC > 0).  R_MC = {tl['R_mc']:+.3f} here, so that rationale does not apply in "
          f"this slice; D is printed, not suppressed.")

    # ---- phi verdict, now gated on the goodness of fit (findings 2 and 9) ------------
    pv = {}
    for si, lab in enumerate(SLICINGS["nHits"]["labels"]):
        cats = []
        gofs = []
        for tmpl in ("cos_spread", "cos_vtx"):
            for pn in ("primary", "secondary"):
                r = T11.get(f"CR/{lab}/TAIL", {}).get(f"{tmpl}/{pn}")
                if not r or not r["fit"]:
                    cats.append("inconclusive")
                    continue
                f = r["fit"]
                gofs.append(f["gof_p"])
                if f["rejected"]:
                    cats.append("inconclusive(model rejected)")
                elif f["fhat"] > 0.25 and f["LB95"] > 0.10:
                    cats.append("cosmic")
                elif f["UB95"] < 0.25:
                    cats.append("not_cosmic")
                else:
                    cats.append("inconclusive")
        vrr = T11.get(f"VR/{lab}/TAIL", {}).get("cos_vtx/primary")
        vr_pos = bool(vrr and vrr["fit"] and not vrr["fit"]["rejected"]
                      and vrr["fit"]["LB95"] > 0.0)
        cat = cats[0] if len(set(cats)) == 1 else "inconclusive"
        if cat == "cosmic" and not vr_pos:
            cat = "inconclusive"
        anti = T11.get(f"CR/{lab}/TAIL", {}).get("cos_spread/primary", {}).get(
            "sigma_obs_vs_nc", float("nan"))
        pv[lab] = {"categories": cats, "verdict": cat, "vr_positive": vr_pos,
                   "gof_p": gofs, "sigma_obs_vs_nc": anti}
        print(f"  phi verdict, CR TAIL {lab:17s}: {cat:28s}  "
              f"(obs-nc) = {anti:+.2f} sigma, GOF p = "
              f"{['%.1e' % g for g in gofs]}")
        if any(g < 0.01 for g in gofs):
            print(f"      the observed phi shape is inconsistent with its own same-slice "
                  f"pT<100 template; the two-template mixture has no solution with f >= 0, "
                  f"so THIS TEST CANNOT BOUND THE COSMIC FRACTION in this slice.  No UB95 "
                  f"is quoted.  The unexplained "
                  f"{'horizontal' if anti < 0 else 'vertical'} excess is itself a result.")
    # The phi test is THREE-STATE (cosmic / not_cosmic / inconclusive, the last with
    # "(model rejected)" when the mixture is refused by its own goodness of fit).  A
    # boolean cosmic_like collapses "the model was REJECTED and bounds nothing" onto the
    # same False as "measured, and not cosmic" -- in the one line a reader quotes.  The
    # status is therefore printed verbatim, and the boolean survives only as the ADOPT
    # trigger, which fires on the affirmative state alone and never on a rejected fit.
    phi_status_low = pv[SLICINGS["nHits"]["labels"][0]]["verdict"]
    phi_rejected_low = "model rejected" in phi_status_low
    cosmic_like = (phi_status_low == "cosmic") and not phi_rejected_low
    adopt = refuses or cosmic_like

    # ---- Before excess vs phi-fitted cosmic content (finding 13) ---------------------
    sub("cross-validation: the Before-channel fold excess vs the phi-fitted cosmic content")
    print("This is the strongest POSITIVE evidence in the study that the fold measures "
          "something real; it is also iteration-dependent, which is stated with it.")
    xval = {}
    for si, lab in enumerate(SLICINGS["nHits"]["labels"]):
        b = T7[("Before", "nHits", si, "TAIL")]
        r = T11.get(f"Before/{lab}/TAIL", {}).get("cos_spread/primary")
        fhat = r["fit"]["fhat"] if (r and r["fit"]) else float("nan")
        gp = r["fit"]["gof_p"] if (r and r["fit"]) else float("nan")
        exc3 = b["obs"] - b["pred"]
        excc = b["obs"] - b["pred_conv"]
        cosn = fhat * b["obs"] if np.isfinite(fhat) else float("nan")
        xval[lab] = {"obs": b["obs"], "excess_iter3": exc3, "excess_conv": excc,
                     "phi_fhat": fhat, "phi_cosmic_counts": cosn, "phi_gof_p": gp}
        print(f"  {lab:17s} Obs={b['obs']:7.0f}  excess(iter 3)={exc3:8.1f}  "
              f"excess(conv)={excc:8.1f}  phi f={fhat:.3f} -> {cosn:8.1f} cosmic muons "
              f"(GOF p={gp:.1e})  ratio excess/phi = "
              f"{exc3/max(cosn,1e-9):.2f} (iter 3) / {excc/max(cosn,1e-9):.2f} (conv)")

    print(f"\n  VERDICT (PROVISIONAL, {HDR}):")
    print(f"    refuses_to_close={refuses}   "
          f"phi_status(low-nHits)={phi_status_low}")
    if phi_rejected_low:
        print(f"        the phi status is NOT a boolean here: the two-template mixture is "
              f"REJECTED by its own goodness of fit in this slice, so it bounds NOTHING "
              f"and must not be read as 'not cosmic'.  No cosmic fraction is quoted.")
    if not validated:
        print(f"    *** VERDICT WITHHELD -- METHOD NOT VALIDATED IN THE SLICE UNDER TEST ***")
        print(f"        failing gates: "
              + ", ".join([g for g, okk in (("7.3 MC self-closure", g73_ok),
                                            ("7.4 Before-data power", g74_ok)) if not okk]))
        _same_side = ("yes" if (tl["R_cons"] < 1.0 / 3.0 and tl["R_bc"] < 1.0 / 3.0)
                      else "NO -- they disagree, escalate")
        print(f"        The pre-registered default (transfer factor only, no nHits cut) is "
              f"what the rule falls back to.  Both the uncorrected closure "
              f"(R = {tl['R']:+.3f}) and the MC-bias-corrected one (R_bc = {tl['R_bc']:+.3f}) "
              f"land on the same side of the MATERIALITY threshold of 1/3 used by "
              f"condition (a): {_same_side}. "
              f"The conclusion is NOT robust to the study's own validation and must go to "
              f"the analyst.")
        print(f"    -> ROUTE TO ANALYST.  Fallback default would be: "
              f"{'ADOPT nHits >= 20' if adopt else 'TRANSFER FACTOR ONLY, no nHits cut'}")
    else:
        print(f"    -> {'ADOPT nHits >= 20' if adopt else 'TRANSFER FACTOR ONLY, no nHits cut'}")
    if adopt:
        print("    signal cost of the adopted cut is 15 % of ALL signal DSA muons "
              "(the 8 % tail-only figure must never be quoted on its own).")
    print(f"    POWER (corrected, fix round findings 1 and 11): condition (a) requires "
          f"Obs >= Pred_up + G_up + Pred/3, i.e. an excess of "
          f"{tl['trigger_excess_frac']:.2f} of Pred -- not 3*sigma_R.  "
          f"sigma_R = {tl['sR']:.3f} (stat only) / {tl['sR_full']:.3f} (stat+syst); "
          f"the method's own demonstrated bias in this slice, |R_MC| = "
          f"{abs(tl['R_mc']):.3f}, is larger than either, so the low-nHits closure is "
          f"BIAS-LIMITED at the ~{abs(tl['R_mc'])*100:.0f} % level, not stat-limited at "
          f"~{tl['sR']*100:.0f} %.  G_up/Pred = {tl['Gup_over_pred']:.3f} "
          f"(a declared ~1.8x OVER-estimate: the {len(mc_present)}-sample MC covers ~96.6 % "
          f"of the tail but only ~53 % of the bulk, so it is conservative AGAINST "
          f"triggering).")
    dec = {"a": {"ok": a_ok, "z_cons": tl["z_cons"], "z_cons_gauss": tl["z_cons_gauss"],
                 "z_cons_clamped": tl["z_cons_clamped"], "R_cons": tl["R_cons"],
                 "n_toys": tl["n_toys"], "trigger_excess_frac": tl["trigger_excess_frac"]},
           "b": {"ok": b_ok, "diff": b_val, "sd": sdiff, "sd_stat": b_stat, "sd_circ": b_sys},
           "c": {"ok": c_ok, "double_diff": c_val, "sd": c_sd,
                 "raw_data_diff": c_raw, "raw_sd": c_raw_sd,
                 "mc_diff": c_mc, "mc_sd": c_mc_sd,
                 "R_high": R_h, "sR_high": sR_h, "R_high_mc": R_hM, "sR_high_mc": sR_hM},
           "d": {"ok": d_ok, "R_cons_VR": vr["R_cons"], "z_cons_VR": vr["z_cons"],
                 "z_cons_gauss_VR": vr["z_cons_gauss"]},
           "e": {"ok": e_ok},
           "gates": {"spec_7_3_mc_selfclosure": g73_ok, "spec_7_4_before_power": g74_ok,
                     "method_validated": validated},
           "closure_three_ways": {
               "as_measured_R": tl["R"], "sR_stat": tl["sR"], "sR_syst": tl["sR_syst"],
               "sR_full": tl["sR_full"], "R_bracket": tl["R_bracket"],
               "R_conv": tl["R_conv"], "R_bias_corrected": tl["R_bc"],
               "sR_bias_corrected": tl["sR_bc"], "D": tl["D"],
               "mc_direction": tl["mc_direction"], "mc_factor": tl["mc_factor"],
               "sign_flips_in_bracket": bool(min(tl["R_bracket"]) * max(tl["R_bracket"]) < 0)},
           "before_phi_cross_validation": xval,
           "refuses_to_close": refuses, "method_failure": method_fail,
           "verdict_withheld": (not validated),
           "phi_verdicts": pv, "phi_status_low_nHits": phi_status_low,
           "phi_model_rejected_low_nHits": phi_rejected_low,
           "cosmic_like_low_nHits": cosmic_like,
           "adopt_nHits_cut": adopt, "provisional": True, "header": HDR}
    RES["T13_decision"] = dec
    cons_ok = all(v["R_cons"] <= v["R"] + 1e-9 for v in T7.values())
    check("R_cons <= R in every slice (all conservative shifts are one-directional)",
          cons_ok, "36 entries")

    # ------------------------------------------------------------------ summary
    hr("SELF-CHECK SUMMARY")
    npass = sum(1 for c in CHECKS if c["pass"])
    for kd, hdr in (("code", "CODE checks -- a FAIL here is a BUG"),
                    ("method", "METHOD-CLAIM checks -- a FAIL means the specification's own "
                               "numerical claim does not hold on this data"),
                    ("expectation", "EXPECTATION checks -- a FAIL is a RESULT, not a defect")):
        sel = [c for c in CHECKS if c["kind"] == kd]
        bad = [c for c in sel if not c["pass"]]
        print(f"\n{hdr}: {len(sel)-len(bad)}/{len(sel)} pass")
        for c in bad:
            print(f"  FAIL: {c['name']} :: {c['detail']}")
    print(f"\n{npass}/{len(CHECKS)} checks PASS overall")
    RES["checks"] = CHECKS
    RES["n_checks_pass"] = npass
    RES["n_checks"] = len(CHECKS)

    hr("LIMITS THE WRITE-UP MUST STATE (spec section 13)")
    LIMITS = [
        "1. Matched muons have a PF partner; LJ muons by definition should not. The matching "
        "failure selects on track quality, so the HYPOTHESIS was T_matched <= T_LJ, making "
        "Pred a LOWER BOUND on feed-up. MEASURED HERE, THAT DIRECTION IS CONTRADICTED: the "
        "spec-section-8 test finds Pred_SB > Obs_SB in a number of slices (see the "
        "'hard transport failure' table), and a lower bound cannot exceed the observation. "
        "Limit 1 must be stated as a refuted hypothesis, not as a fact; the sideband "
        "non-closure in a slice is the floor on any tail-closure claim in that slice. "
        "Slicing reduces but cannot remove the population mismatch; a joint |eta| x nHits "
        "matrix is not booked.",
        "2. The two populations are NOT disjoint in code: the LJ cross-cleaning is a SEGMENT "
        "match, the control a dR < 0.2 match, so a dR-matched but segment-unmatched DSA muon "
        "sits in both M and L/Obs. The overlap is unmeasurable with the booked histograms; it "
        "biases closure toward 1 and under-states errors. A segment-match flag (or a dR axis) "
        "on the matched matrices is the recommended addition for the next rerun.",
        "3. The fold explains the LJ tail from the matched-sample rates; it does not explain "
        "why data feeds up ~3.5x more than MC -- that is a detector question.",
        "4. Stage 1 addresses only feed-up from below the source cut (~26 % of the matched "
        "sample's own tail); the residual BOUNDS, never measures, a genuine high-pT component.",
        "5. Circularity is corrected within the source block and bracketed by sigma_circ, not "
        "removed.  MEASURED HERE: the naive fold's bias does NOT have the universal high sign "
        "the specification assumed, and the EM recursion has not converged at the pre-registered "
        "iteration 3 (see the METHOD-CLAIM checks); the remaining distance to its fixed point is "
        "folded into sigma_circ, which can only lower R_cons.",
        "6. Per-slice MC closure is a +-30-50 % statement on ~10-20 effective events and "
        "confounds transport with genuine content. " + (
            f"MC samples missing: {mc_missing} (QCD_Pt800To1000 is the sample most populated at "
            "true high pT): G_MC is biased low, hence G_up and R_cons are biased toward "
            "triggering by an unknown amount until it lands."
            if mc_missing else
            f"All {len(mc_present)}/8 MC samples are present, so this bias is retired; the "
            "statement remains that the test rests on ~10-20 effective events."),
        f"7. No bulk-normalised rate may be compared across data and MC with the "
        f"{len(mc_present)}/8 tail-selected set (R4); G_up is a declared ~1.8x over-estimate "
        "for that reason.",
        "8. The [1.2,2.4) |eta| slice averages a factor-~4 internal spread in the data/MC "
        "transfer ratio (peak near |eta| 1.6-1.8, trough near 2.2-2.4); kappa_pop is the only handle.",
        "9. Wrong dR partners inflate T at the 1e-3 level where the low-row tail probabilities "
        "live; the relPtErr table exposes but cannot bound it (a 0.2 < dR < 0.4 control would).",
        "10. Cosmic templates are data-only and selection-biased (R8: the vertex-inverted "
        "template is ~94 % high-nHits and selects on cosAlpha, which a short track reconstructs "
        "worst); the phi test bounds only a VERTICAL component.",
        "11. CR and VR share most events; the |eta| and nHits slicings share the same muons; "
        "no quadrature combination of their pulls is legitimate. Two DSA muons of one LJ are "
        "correlated; the errors assume independence.",
        f"12. Data are {len(present)}/68 chunks. Closures need no luminosity, but T14 and the "
        "era table depend on an external era-lumi table until 68/68 land; the verdict is PROVISIONAL.",
        "13. The thresholds 3 sigma / 1/3 / 2 sigma / 0.25 / 0.10 are PRE-REGISTERED CHOICES, "
        "not derived numbers. 'Consistent with feed-up at factor-1.5 precision' remains the "
        "strongest honest formulation unless R_cons and A_low are flat across the scan.",
    ]
    for l in LIMITS:
        print("  " + l)
    RES["limits"] = LIMITS

    with open(args.json, "w") as f:
        json.dump(RES, f, indent=1, default=jnum)
    print(f"\nJSON written: {args.json}")


if __name__ == "__main__":
    main()
