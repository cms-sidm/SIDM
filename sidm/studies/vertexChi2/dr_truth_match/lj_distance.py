"""Which condition of the dimuon-vertex producer leaves a mu-type LJ without a stored vertex.

usage: python lj_distance.py REC_2MU2E REC_4MU REC_WP CLS_2MU2E CLS_4MU CLS_WP OUT_DIR

Reads the per-LJ records written by the extended extract.py (REC_2MU2E: the sixteen 2Mu2E
samples, REC_4MU: the nine 4Mu samples, REC_WP: the nine wide-pair 2Mu2E samples at the mass
points of the 4Mu samples; directories of parquet files, one record per mu-type LJ with its PF
and DSA constituents, the DSA track parameters and track-quality branches, the PF inner-track
parameters, and the stored vertex of the two leading muons, pv_*) and the lj_classified.parquet
files written by analyze.py for the same samples (CLS_*: Lxy, generated dR and class A to E of
each gen-matched LJ). Writes CSV tables and summary.txt (the checks and the main rows of every
table) to OUT_DIR (must exist), and prints the summary.

The producer. MuonVertexTableProducer (LLPNanoAOD, module muonVertexTable with dsaMuons =
displacedStandAloneMuons and patMuons = linkedObjects:muons) fits each muon pair with a
KalmanVertexFitter and stores no vertex when the fit is invalid, or when TwoTrackMinimumDistance
(TTMD) between the two tracks succeeds and gives a distance above 15 cm. The stored dca is that
distance and dcaStatus the TTMD status. TTMD is given a FreeTrajectoryState per track, built
from the reco::Track's reference point (vx, vy, vz), momentum and charge, in the field of
IdealMagneticFieldRecord. The track of a DSA muon is its displacedStandAloneMuons track; that of
a PF muon is the combined track of a global muon, else the standalone track of a standalone
muon, else tunePMuonBestTrack. Pairs are ordered PF before DSA, then by index.

The emulation (ttmd below) is a port of TTMD in FastMode (CMSSW_10_6_29): Newton steps on the
two helix phases from the two states (at most 5 accepted, a step out of more than 20 cm of
transverse arc fails); on failure, the transverse crossing (or closest approach) of the two
circles, ClosestApproachInRPhi, followed by the Newton steps from those points, the transverse
points being kept when that gives a larger distance. It returns the distance, the crossing point
(the mean of the two closest points, as TTMD::crossingPoint) and the path taken. The field is
the parabolic parametrisation Bz(r, z) of CMSSW, evaluated where TTMD evaluates it. A DSA leg
uses the stored DSA track, as the producer does. A PF leg uses the inner track's reference point
and momentum: the producer's track for a tracker-only muon, a stand-in for the combined or
standalone track otherwise; pairs with such a leg are not reproduced (validation_dca.csv), and
the distance enters the tables only for DSA-DSA LJs. validation_dca.csv compares the emulated
with the stored dca for the pairs with a stored vertex, in bins of the stored dca, per group,
with the parabolic field and with a uniform 3.8 T field. Checks: for every LJ with exactly two
muons, a stored vertex of the two (pv_n > 0) exists exactly when vtx_chi2 >= 0; the
opposite-charge LJ counts per content (table15_check.csv, table 15 of pair_content.py); the
emulated status against the stored one.

Per leading DSA-DSA pair, besides the emulated distance and crossing point: the 3D opening angle
of the two helices at the crossing point (the TTMD helices, each with Bz at the other track's
reference point; tangents at the projection of the crossing point on each circle); the global
minimum distance of the two helices (global_min: both closest points inside the fitter bounds
of the production, r < 740 cm and |z| < 960 cm, transverse arc |s| <= 1500 cm; candidates from a
20 cm grid, from exact one-parameter profiles at 5 cm steps and from 2 cm grids around the best
grid minima, each refined by a pattern search and Newton steps; gmin_stability compares it with
a finer variant and with an exhaustive 1 cm grid); and L = s / angle with s = 0.5 sqrt(dxyBSErr1^2
+ dxyBSErr2^2 + dzBSErr1^2 + dzBSErr2^2) (the records hold the beam-spot dxy and dz errors of the
DSA tracks, not dxyErr and dzErr; for these tracks the beam-spot width is negligible).

Populations. P(set): gen-matched LJs with exactly two muons, both DSA (DSA-DSA), of opposite
charge, of a sample set (2Mu2E: the sixteen 2Mu2E samples, the population of T1 and of the
tests unless stated; 2Mu2E wide-pair; 4Mu). T6 population: 4Mu, Lxy >= 100 cm, two-muon DSA-DSA
LJs of class A and of class D, generated dR [0.1, 0.2) (narrow) or [0.2, 0.5) (wide).
no-vertex = vtx_chi2 < 0. All LJs: every LJ whose two leading muons are DSA (no-vertex =
pv_n == 0), used for the status checks, the same-mass-point table, the 1 cm edge table and T7.
The two conditions: cut = emulated status 1 and distance > 15 cm; s15 = share of the LJs with
the cut; f_inv = no-vertex fraction of the other LJs (distance <= 15 cm: the invalid Kalman
fits); no-vertex fraction = s15 + (1 - s15) f_inv when every cut LJ has no vertex. Errors are
binomial, sqrt(f (1 - f) / n). Bins: distance [0, 0.1, 0.5, 1, 2, 5, 10, 15, inf) cm (coarse:
[0, 0.5, 2, 15, inf)); angle [0, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1, inf) rad; L
[0, 1, 3, 10, 30, 100, 300, 1e3, 3e3, 1e4, inf) cm.

Tables:
- T1_novtx_vs_distance.csv: no-vertex fraction against the distance, all and per mass point
  (M, m_Zd); T1_distance_shares_by_point.csv: the distance composition per mass point.
- T2_condition_share.csv: per sample set, share of the no-vertex LJs of P(set) with distance
  > 15 cm, by Lxy slice, region and mass point. T2_same_mass_points.csv: all LJs, per sample
  set and sample, the no-vertex LJs and those with distance > 15 cm.
- decomposition_cut_vs_invalid_fit.csv: s15 and f_inv of the two groups of each comparison,
  (a) (1000, 5) against (500, 5) GeV, (b) outer [150, 300) against inner [0, 60) cm, (c) endcap
  against barrel per Lxy slice below 200 cm and pooled, (d) 4Mu wide against narrow in class A
  and in class D; the difference of the no-vertex fractions split exactly into
  term_s15 = (1 - mean f_inv) (s15_b - s15_a) and term_finv = (1 - mean s15) (f_inv_b - f_inv_a)
  (residual: cut LJs with a stored vertex, emulation mismatches). diff = f_b - f_a, the
  no-vertex fractions; group a is the reference, the second-named group of the comparison label
  (columns group_a, group_b).
- invalid_fit_tests.csv: f_inv of the two groups of (a) to (d), LJs with distance <= 15 cm,
  compared inside cells (Mantel-Haenszel style: z = sum (k_B - n_B p) / sqrt(sum n_A n_B p
  (1 - p) / (n_A + n_B - 1)), weighted difference wdiff with weights n_A n_B / (n_A + n_B),
  cells with >= nmin LJs of each group, kept = share of the LJs of the two groups inside the
  cells used, base_wdiff and base_z the same comparison on the same LJs inside the base cells
  only). Base cells: (a) Lxy slice cut in 8 x generated dR (table 18 edges), nmin 30 and 5;
  (b) none; (c) Lxy slice; (d) none; nmin 5 for (b) to (d). Further cell variable: none,
  distance bin, angle bin, angle bin x coarse distance, and each of 13 track-quality variables
  (worst of the two muons, quintile bins over P(2Mu2E) or the T6 population).
- T4_standardised.csv: no-vertex fraction per region and Lxy slice, observed and standardised
  to the pooled distance composition of P(2Mu2E).
- invalid_fit_universality.csv: one f_inv shared by all groups in every cell (chi2, ndf, cells,
  LJs, kept, extra = the spread between groups beyond the binomial one; groups with >= 50 and
  >= 10 LJs in the cell), groups (i) the mass points of P(2Mu2E), (ii) barrel and endcap,
  (iii) the Lxy slices, (iv) the three sample sets; cells none, distance, angle, angle x
  distance.
- invalid_fit_candidates.csv: f_inv per sample set against the crossing-point radius and |z|,
  the fitter bounds, the TTMD path, the angle, dR of the two momenta at the reference points and
  the distance; invalid_fit_angle_x_distance.csv: against angle x distance (both distance bins).
- invalid_fit_deviance.csv: binomial deviance -2 sum [k ln f + (n - k) ln(1 - f)] of f_inv
  described per bin of distance, angle, angle x distance, L, L x distance (and one bin), on the
  same LJs (distance <= 15 cm, finite L), per sample set and pooled; aic = deviance + 2 bins.
- invalid_fit_angle_standardised.csv: the same LJs, per angle bin (and per angle_eff bin) and
  sample set: n, k, f_inv, binomial error, and z of the set against 2Mu2E in the bin (as in
  the tests, one cell); per set other than 2Mu2E: observed f_inv and f_inv predicted from the
  per-bin f_inv of 2Mu2E with the set's own bin shares (direct standardisation; bins empty in
  2Mu2E dropped, kept = share of the set's LJs in the other bins, f_obs_kept the observed f_inv
  on them), ratio = f_obs_kept / f_pred, with bootstrap errors as in wide_pair.py (NBOOT
  replicas resampling the LJs of 2Mu2E and of the set separately, multinomial over bin x
  outcome, fixed seed; *_bs_err the standard deviation, ratio_bs_q16 / _q84 the 16 % and 84 %
  quantiles).
- invalid_fit_sin_theta.csv: P(2Mu2E), distance <= 15 cm, inside the barrel and inside the
  endcap: f_inv in quintiles of sin(theta_pair), and the last quintile against the first inside
  (angle bin x coarse distance) cells. theta_pair: polar angle of the summed momentum of the two
  tracks at the crossing point; angle_eff = angle / sin(theta_pair), binned as the angle with one
  more bin, [1, 3) and >= 3 rad (specified from the fitter's transverse convergence criterion
  before it was computed).
- edge_1cm.csv: all LJs, no-vertex fraction in 1 cm bins of the distance from 10 to 20 cm.
- T7_global_minimum_cut_pairs.csv: all LJs without a vertex and distance > 15 cm (and the T6
  class A narrow and wide LJs among them): global minimum below 1, 5 and 15 cm.
  T7_stored_dca_vs_global_minimum.csv: all LJs with a stored vertex, per stored-dca bin, the
  share whose stored dca exceeds the global minimum by more than 1 cm (and falls below it by
  more than 0.1 cm). T7_minimiser_stability.csv: gmin_stability per sample set and the path-1
  check (the global minimum never above the TTMD distance of the same helices).
  T7_global_minimum_location.csv: per sample set, the cut pairs of T7 (all; global minimum
  below 15 cm; at or above 15 cm) and, as control, all LJs whose leading DSA-DSA pair has a
  stored vertex: where the point of the global minimum (the midpoint of the two closest points)
  lies against the generated decay. The generated quantity is gen_a_lxy, the transverse
  distance between the production and decay vertices of the dark photon matched to the LJ (the
  nearest genAs_toMu, as mu_lj_gen_a_lxy finds it; the sidm lxy helper); the LJ is matched to
  it as everywhere in this script, gen_matched = that dark photon within dR < 0.4 of the LJ
  (n_gen: gen-matched with gen_a_lxy >= 0; n_used: of those, a finite global minimum). The
  production vertex is the prompt mediator decay, so the radius r of a point is compared with
  gen_a_lxy directly: dr = |r - gen_a_lxy| (median, 16 % and 84 % quantiles, shares < 10 cm,
  < 20 cm, > 100 cm), the shares of points at r > 300 cm and at |z| > 600 cm, and the same dr
  for the emulated TTMD crossing point (n_cross: gen-matched pairs with a crossing point). The
  records hold the dark photon's Lxy, eta and phi but not its decay vertex in 3D, so no 3D
  distance is given; they hold no stored vertex position either (pv_*: n, dca, dcaStatus,
  isValid, normChi2), so the control row uses the same two points (global-minimum midpoint and
  crossing point) of the pairs with a stored vertex, not the stored vertex. All pairs are used.

The sample 2Mu2E_100GeV_5p0GeV_400mm is in both 2Mu2E record sets (REC_2MU2E and REC_WP): its
LJs count in both sets in invalid_fit_universality.csv (iv) and in the pooled rows of
invalid_fit_deviance.csv; they are not de-duplicated, and summary.txt gives the samples found in
both sets and the number of LJs of P(set) in them.
"""
import argparse
import glob
import math
import os

import awkward as ak
import numpy as np
import pandas as pd

C_GEV = 2.99792458e-3         # curvature in 1/cm per (field in T / momentum in GeV)
DCA_EDGES = [0, 0.1, 0.5, 1, 2, 5, 10, 15]
DIST_EDGES = [0, 0.1, 0.5, 1, 2, 5, 10, 15, np.inf]
DIST_COARSE = [0, 0.5, 2, 15, np.inf]
ANG_EDGES = [0, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1, np.inf]
AEFF_EDGES = ANG_EDGES[:-1] + [3, np.inf]   # angle_eff: the angle bins, one more at the top
L_EDGES = [0, 1, 3, 10, 30, 100, 300, 1e3, 3e3, 1e4, np.inf]
R_EDGES = [0, 50, 100, 200, 300, 400, 500, 600, 700, 740, 800, 1000, np.inf]
Z_EDGES = [0, 100, 300, 500, 700, 900, 960, 1100, np.inf]
R_MAX, Z_MAX, S_MAX = 740.0, 960.0, 1500.0   # fitter bounds of the production; arc range (cm)
NBOOT, SEED = 1000, 20261114                 # bootstrap of the standardisation, as wide_pair.py
LXY_DSA = [0, 30, 60, 110, 150, 200, 250, 300, 350, 400, 500, np.inf]
LXY_CONTENT = [0, 1, 10, 30, 60, 110, 300, 1e9]
GEN_EDGES = [0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 10]
GEN_EDGES_FINE = [0, 0.0005, 0.001] + GEN_EDGES[1:]
SETS = ("2Mu2E", "2Mu2E wide-pair", "4Mu")
SCALARS = ["sample", "run", "lumi", "event", "lj_pt", "lj_eta", "lj_phi", "n_pf", "n_dsa",
           "vtx_chi2", "gen_matched", "gen_a_lxy", "d1_pt", "d1_eta", "d1_phi", "d2_pt",
           "d2_eta", "d2_phi"]
PV = ["pv_n", "pv_dca", "pv_dcaStatus", "pv_isValid", "pv_normChi2"]
MU = ["type", "idx", "pt", "eta", "phi", "charge", "vx", "vy", "vz", "normChi2", "chi2", "ndof",
      "ptErr", "etaErr", "phiErr", "nSegments", "nDTSegments", "nCSCSegments", "trkNumHits",
      "trkNumDTHits", "trkNumCSCHits", "trkNumPlanes", "dxyBSErr", "dzBSErr", "displacedID",
      "innerVx", "innerVy", "innerVz", "innerPt", "innerEta", "innerPhi", "isGlobal",
      "isStandalone", "isTracker"]
# track-quality variables, worst of the two muons: (name, field or derived, "max" or "min")
QUALITY = [("normChi2", "normChi2", "max"), ("ptErr/pt", "ptrel", "max"),
           ("etaErr", "etaErr", "max"), ("phiErr", "phiErr", "max"),
           ("dxyBSErr", "dxyBSErr", "max"), ("dzBSErr", "dzBSErr", "max"),
           ("ndof", "ndof", "min"), ("trkNumHits", "trkNumHits", "min"),
           ("trkNumPlanes", "trkNumPlanes", "min"), ("nSegments", "nSegments", "min"),
           ("trkNumDTHits", "trkNumDTHits", "min"), ("trkNumCSCHits", "trkNumCSCHits", "min"),
           ("displacedID", "displacedID", "min")]
SUMMARY = []


def say(*a):
    """Print and keep a summary line."""
    line = " ".join(str(x) for x in a)
    print(line)
    SUMMARY.append(line)


def load(directory):
    """One row per LJ from every parquet file of directory: scalars, pv_* and the fields MU of
    the two leading muons (m0_*, m1_*, NaN when absent)."""
    files = sorted(glob.glob(os.path.join(directory, "*.parquet")))
    r = ak.concatenate([ak.from_parquet(f) for f in files])
    d = {c: (np.asarray(ak.to_list(r[c])) if c == "sample" else ak.to_numpy(r[c]))
         for c in SCALARS}
    d["n_mu"] = ak.to_numpy(ak.num(r.muons, axis=1))
    lead = ak.pad_none(r.muons, 2, axis=1)
    for i in (0, 1):
        for f in MU:
            d[f"m{i}_{f}"] = ak.to_numpy(ak.fill_none(lead[:, i][f], np.nan)).astype(float)
    for c in PV:
        d[c] = ak.to_numpy(r[c])
    df = pd.DataFrame(d)
    mm = df["sample"].str.extract(r"_(\d+)GeV_(\d+p\d+)GeV")
    df["mass"] = mm[0].astype(int)
    df["mzd"] = mm[1].str.replace("p", ".").astype(float)
    return df, len(files)


# ---------------------------------------------------------------- TwoTrackMinimumDistance port
def bz_parabolic(x, y, z):
    """Bz in T of ParabolicParametrizedMagneticField (CMSSW_10_6_29), without its validity cut
    (r < 115 cm, |z| < 280 cm)."""
    r2 = x * x + y * y
    return (-3.94991e-06 * z * z + 7.53701e-06 * z + 3.8114) * (2.43878e-11 * r2 + 1.0)


def bz_uniform(x, y, z):
    """Uniform 3.8 T."""
    return 3.8 + 0.0 * (x + y + z)


def _hh(G, H, qual, bz_tesla):
    """TwoTrackMinimumDistanceHelixHelix::calculate(G, H, qual). G, H: dicts of x, y, z, px,
    py, pz, q arrays. Returns (fail, point on G (3, n), point on H (3, n)). The field of H is
    taken at G's position and that of G at H's position, as in CMSSW."""
    bH = bz_tesla(G["x"], G["y"], G["z"]) * C_GEV
    bG = bz_tesla(H["x"], H["y"], H["z"]) * C_GEV
    Ht = np.hypot(H["px"], H["py"])
    Gt = np.hypot(G["px"], G["py"])
    fail = (Ht == 0) | (Gt == 0) | (H["q"] == 0) | (G["q"] == 0) | (bG == 0) | (bH == 0)
    with np.errstate(all="ignore"):
        h = Ht / (H["q"] * bH)
        sH0, cH0, tH = -H["py"] / Ht, -H["px"] / Ht, -H["pz"] / Ht
        pH0 = np.arctan2(sH0, cH0)
        g = Gt / (G["q"] * bG)
        sG0, cG0, tG = -G["py"] / Gt, -G["px"] / Gt, -G["pz"] / Gt
        pG0 = np.arctan2(sG0, cG0)
        a = H["x"] - G["x"] + g * sG0 - h * sH0
        b = H["y"] - G["y"] - g * cG0 + h * cH0
        c1, c2 = h * tH * tH, -g * tG * tG
        d1, d2 = -g * tG * tH, h * tG * tH
        e0 = H["z"] - G["z"] - h * pH0 * tH + g * pG0 * tG
        e1, e2 = tH * e0, tG * e0
        pH, pG = pH0.copy(), pG0.copy()
        active = ~fail
        for counter in range(6):
            if not active.any():
                break
            sH, cH, sG, cG = np.sin(pH), np.cos(pH), np.sin(pG), np.cos(pG)
            A11 = h * (-sH * (a - g * sG) + cH * (b + g * cG) + c1)
            A22 = -g * (-sG * (a + h * sH) + cG * (b - h * cH) + c2)
            A12 = h * (-g * cG * cH - g * sH * sG + d1)
            A21 = -g * (h * cG * cH + h * sH * sG + d2)
            det_i = 1.0 / (A11 * A22 - A12 * A21)
            z1 = h * (cH * (a - g * sG) + sH * (b + g * cG) + c1 * pH + d1 * pG + e1)
            z2 = -g * (cG * (a + h * sH) + sG * (b - h * cH) + c2 * pG + d2 * pH + e2)
            dH = (z1 * A22 - z2 * A12) * det_i
            dG = (z2 * A11 - z1 * A21) * det_i
            bad = (A11 < 0) | (A22 < 0)
            sing = ~bad & (np.abs(det_i) > 10.0)
            step = active & ~bad & ~sing
            pH = np.where(step, pH - dH, pH)
            pG = np.where(step, pG - dG, pG)
            ret = bad | sing | ~np.isfinite(dH) | ~np.isfinite(dG) | (counter > 4)
            fail |= active & ret
            conv = (np.abs(dH) <= qual) & (np.abs(dG) <= qual)
            active = active & ~ret & ~conv
        fail |= active
        fail |= np.abs(g * (pG - pG0)) > 20
        fail |= np.abs(h * (pH - pH0)) > 20
        PG = np.stack([G["x"] + g * (np.sin(pG) - sG0), G["y"] + g * (-np.cos(pG) + cG0),
                       G["z"] + g * tG * (pG - pG0)])
        PH = np.stack([H["x"] + h * (np.sin(pH) - sH0), H["y"] + h * (-np.cos(pH) + cH0),
                       H["z"] + h * tH * (pH - pH0)])
    fail |= ~np.all(np.isfinite(PG), axis=0) | ~np.all(np.isfinite(PH), axis=0)
    return fail, PG, PH


def _circle(t, bz):
    """Centre and radius of the transverse circle of state t (bz in curvature units)."""
    qob = t["q"] / bz
    pt = np.hypot(t["px"], t["py"])
    return t["x"] + qob * t["py"], t["y"] - qob * t["px"], np.abs(qob * pt)


def _zcoord(t, r, xc, yc, xg, yg):
    """ClosestApproachInRPhi::zCoord: z of state t at the transverse point (xg, yg)."""
    x, y, z, px, py, pz = t["x"], t["y"], t["z"], t["px"], t["py"], t["pz"]
    with np.errstate(all="ignore"):
        shp = np.sqrt((x - xg) ** 2 + (y - yg) ** 2) / (2 * r)
        cphi = np.clip(((x - xc) * (xg - xc) + (y - yc) * (yg - yc)) / (r * r), -1, 1)
        phi = np.where(shp < 0.383, 2 * np.arcsin(np.minimum(shp, 1)), np.abs(np.arccos(cphi)))
    sphi = np.where((x - xc) * (yg - yc) - (xg - xc) * (y - yc) > 0, 1.0, -1.0)
    som = np.where((x - xc) * py - (y - yc) * px > 0, 1.0, -1.0)
    return z + sphi * som * (pz / np.hypot(px, py)) * phi * r


def _cair(A, B, bz_tesla):
    """ClosestApproachInRPhi::calculate(A, B), one field (at A's position) for both circles:
    the transverse crossing with the smaller |dz|, or the closest points of two circles that do
    not cross. Returns (ok, point on A (3, n), point on B (3, n), bz)."""
    bz = bz_tesla(A["x"], A["y"], A["z"]) * C_GEV
    xa, ya, ra = _circle(A, bz)
    xb, yb, rb = _circle(B, bz)
    with np.errstate(all="ignore"):
        d = np.hypot(xb - xa, yb - ya)
        u, v = (xb - xa) / d, (yb - ya) / d
        f1 = (d <= ra + rb) & (d >= np.abs(rb - ra)) & (d > 0)
        f2o = (d > ra + rb) & (d > 0)
        f2i = (d < np.abs(rb - ra)) & (d > 0)
        cphi = (ra * ra - rb * rb + d * d) / (2 * ra * d)
        s2 = 1 - cphi * cphi
        cphi = np.where(s2 < 0, 1.0, cphi)
        sphi = np.sqrt(np.maximum(s2, 0))
        sign = np.where(ra <= rb, -1.0, 1.0)
        x1 = np.select([f1, f2o, f2i], [ra * cphi, ra, sign * ra], np.nan)
        y1 = np.select([f1, f2o, f2i], [ra * sphi, 0 * ra, 0 * ra], np.nan)
        x2 = np.select([f1, f2o, f2i], [ra * cphi, d - rb, d + sign * rb], np.nan)
        y2 = np.select([f1, f2o, f2i], [-ra * sphi, 0 * ra, 0 * ra], np.nan)
        xg1, yg1 = u * x1 - v * y1 + xa, v * x1 + u * y1 + ya
        xg2, yg2 = u * x2 - v * y2 + xa, v * x2 + u * y2 + ya
    ok = f1 | f2o | f2i
    za1 = _zcoord(A, ra, xa, ya, xg1, yg1)
    zb1 = _zcoord(B, rb, xb, yb, xg1, yg1)
    za2 = _zcoord(A, ra, xa, ya, xg2, yg2)
    zb2 = _zcoord(B, rb, xb, yb, xg2, yg2)
    first = np.abs(zb1 - za1) < np.abs(zb2 - za2)
    xA = np.where(f1, np.where(first, xg1, xg2), xg1)
    yA = np.where(f1, np.where(first, yg1, yg2), yg1)
    zA = np.where(f1, np.where(first, za1, za2), za1)
    xB = np.where(f1, xA, xg2)
    yB = np.where(f1, yA, yg2)
    zB = np.where(f1, np.where(first, zb1, zb2), _zcoord(B, rb, xb, yb, xg2, yg2))
    return ok, np.stack([xA, yA, zA]), np.stack([xB, yB, zB]), bz


def _new_traj(p, t, bz):
    """ClosestApproachInRPhi::newTrajectory: the state t moved along its circle to point p."""
    qob = t["q"] / bz
    xc, yc = t["x"] + qob * t["py"], t["y"] - qob * t["px"]
    return dict(x=p[0], y=p[1], z=p[2], px=(p[1] - yc) * (bz / t["q"]),
                py=(xc - p[0]) * (bz / t["q"]), pz=t["pz"], q=t["q"])


def _sub(t, i):
    """The entries i of every array of state t."""
    return {k: v[i] for k, v in t.items()}


def ttmd(A, B, bz_tesla):
    """TwoTrackMinimumDistance (FastMode) calculate(A, B) for two charged states. Returns
    (status, distance, crossing point (3, n), path): the crossing point is the mean of the two
    closest points (TTMD::crossingPoint); path 1 = Newton steps from the states, 2 = transverse
    step failed and Newton steps with qual 0.1, 3 = Newton steps from the transverse points,
    4 = transverse points kept."""
    n = len(A["x"])
    status = np.zeros(n, bool)
    dist = np.full(n, np.nan)
    path = np.zeros(n, int)
    cross = np.full((3, n), np.nan)
    f1, PG, PH = _hh(A, B, 1e-4, bz_tesla)
    ok = ~f1
    status[ok], path[ok] = True, 1
    dist[ok] = np.linalg.norm(PG[:, ok] - PH[:, ok], axis=0)
    cross[:, ok] = 0.5 * (PG[:, ok] + PH[:, ok])
    rest = np.where(f1)[0]
    if len(rest) == 0:
        return status, dist, cross, path
    A2, B2 = _sub(A, rest), _sub(B, rest)
    cok, pa, pb, bz = _cair(A2, B2, bz_tesla)
    i0 = rest[~cok]
    if len(i0):
        f2, PG2, PH2 = _hh(_sub(A, i0), _sub(B, i0), 0.1, bz_tesla)
        s = i0[~f2]
        status[s], path[s] = True, 2
        dist[s] = np.linalg.norm(PG2[:, ~f2] - PH2[:, ~f2], axis=0)
        cross[:, s] = 0.5 * (PG2[:, ~f2] + PH2[:, ~f2])
    j = np.where(cok)[0]
    i3 = rest[j]
    if len(i3):
        pa3, pb3, bz3 = pa[:, j], pb[:, j], bz[j]
        Aj, Bj = _sub(A2, j), _sub(B2, j)
        f3, PG3, PH3 = _hh(_new_traj(pa3, Aj, bz3), _new_traj(pb3, Bj, bz3), 1e-4, bz_tesla)
        dc = np.linalg.norm(pa3 - pb3, axis=0)
        d3 = np.linalg.norm(PG3 - PH3, axis=0)
        use3 = ~f3 & (d3 <= dc)
        status[i3] = True
        path[i3] = np.where(use3, 3, 4)
        dist[i3] = np.where(use3, d3, dc)
        cross[:, i3] = np.where(use3, 0.5 * (PG3 + PH3), 0.5 * (pa3 + pb3))
    return status, dist, cross, path


def leg(df, i):
    """State of leading muon i as the producer builds it, and the kind of leg: dsa (the
    producer's track), pf_inner (tracker-only PF muon: the producer's track), pf_global (inner
    track in place of the combined track), pf_sa (inner track in place of the standalone
    track)."""
    t = df[f"m{i}_type"].to_numpy()
    dsa = t == 1
    glb = (t == 0) & (df[f"m{i}_isGlobal"].to_numpy() > 0)
    sa = (t == 0) & ~glb & (df[f"m{i}_isStandalone"].to_numpy() > 0)
    inner = (t == 0) & ~glb & ~sa
    kind = np.select([dsa, inner, glb, sa], ["dsa", "pf_inner", "pf_global", "pf_sa"], "none")

    def pick(dsa_col, pf_col):
        return np.where(dsa, df[f"m{i}_{dsa_col}"], df[f"m{i}_{pf_col}"]).astype(float)

    pt, eta, phi = pick("pt", "innerPt"), pick("eta", "innerEta"), pick("phi", "innerPhi")
    pt = np.where(pt > 0, pt, np.nan)
    st = dict(x=pick("vx", "innerVx"), y=pick("vy", "innerVy"), z=pick("vz", "innerVz"),
              px=pt * np.cos(phi), py=pt * np.sin(phi), pz=pt * np.sinh(eta),
              q=df[f"m{i}_charge"].to_numpy().astype(float))
    return st, kind


def metrics(x, col):
    """Emulated (col) minus stored dca: median, 16 and 84% quantiles, shares within 0.1 cm and
    within 10% of the stored value, shares low and high by more than 1 cm (a failed emulation
    counts as neither within nor low or high), share with a failed emulation."""
    if len(x) == 0:
        return dict(n=0)
    diff = x[col] - x.pv_dca
    return dict(n=len(x), median=diff.median(), q16=diff.quantile(0.16),
                q84=diff.quantile(0.84), within_0p1=(diff.abs() < 0.1).mean(),
                within_10pc=(diff.abs() < 0.1 * x.pv_dca).mean(), low_gt1=(diff < -1).mean(),
                high_gt1=(diff > 1).mean(), failed=diff.isna().mean())


def validation(groups):
    """Emulated against stored dca for the pairs of each (name, frame) group with a stored
    vertex, dcaStatus 1 and both legs available, per bin of the stored dca, for both fields."""
    rows = []
    for name, x in groups:
        x = x[(x.pv_n > 0) & (x.pv_dcaStatus == 1) & (x.pv_dca >= 0) & (x.emu_status >= 0)]
        for field, col in (("parabolic", "emu"), ("uniform 3.8 T", "emu_u")):
            rows.append(dict(group=name, field=field, dca_bin="all", **metrics(x, col)))
            for b, g in x.groupby(pd.cut(x.pv_dca, DCA_EDGES, right=False), observed=False):
                rows.append(dict(group=name, field=field, dca_bin=f"[{b.left:g},{b.right:g})",
                                 **metrics(g, col)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- tables
def split(edges, k):
    """Each finite slice of edges cut in k parts (geometric above 0, linear from 0); the open
    last slice kept."""
    out = [edges[0]]
    for lo, hi in zip(edges[:-2], edges[1:-1]):
        cut = np.geomspace(lo, hi, k + 1) if lo > 0 else np.linspace(lo, hi, k + 1)
        out.extend(cut[1:-1].tolist() + [hi])
    return out + [edges[-1]]


def frac(k, n):
    """Fraction and binomial error."""
    f = k / n if n else float("nan")
    return f, (math.sqrt(f * (1 - f) / n) if n else float("nan"))


def rate_table(df, keys):
    """n, k (no vertex), f and err per group of keys."""
    t = df.groupby(keys, observed=True).no_vtx.agg(n="size", k="sum").reset_index()
    t["f"] = t.k / t.n
    t["err"] = np.sqrt(t.f * (1 - t.f) / t.n)
    return t


def mh(df, col, a, b, keys, nmin=30, baseline=True):
    """Group b against group a of column col inside cells of keys (docstring, T3). kept: share
    of the LJs of the two groups inside the cells used. With baseline, also the same comparison
    on the same LJs (those of the cells kept) inside the base cells only: the first two keys
    when they are (lc8, gc), else the first key (base_wdiff, base_z), so that a cell variable is
    judged on a fixed set of LJs."""
    empty = dict(cells=0, n_a=0, n_b=0, kept=0.0, f_a=np.nan, f_b=np.nan, wdiff=np.nan,
                 z=np.nan, base_wdiff=np.nan, base_z=np.nan)
    if not ((df[col] == a).any() and (df[col] == b).any()):
        return empty
    x = df[df[col].isin([a, b])]
    agg = x.groupby(keys + [col], observed=True).no_vtx.agg(["sum", "size"]).unstack(col)
    agg = agg.dropna()
    if agg.empty or ("size", a) not in agg.columns or ("size", b) not in agg.columns:
        return empty
    agg = agg[(agg[("size", a)] >= nmin) & (agg[("size", b)] >= nmin)]
    ka, na, kb, nb = agg[("sum", a)], agg[("size", a)], agg[("sum", b)], agg[("size", b)]
    p0 = (ka + kb) / (na + nb)
    ok = (p0 > 0) & (p0 < 1)
    ka, na, kb, nb, p0 = ka[ok], na[ok], kb[ok], nb[ok], p0[ok]
    if len(p0) == 0:
        return empty
    var = (nb * na * p0 * (1 - p0) / (na + nb - 1)).sum()
    w = na * nb / (na + nb)
    out = dict(cells=len(p0), n_a=int(na.sum()), n_b=int(nb.sum()),
               kept=float(na.sum() + nb.sum()) / len(x), f_a=ka.sum() / na.sum(),
               f_b=kb.sum() / nb.sum(), wdiff=(w * (kb / nb - ka / na)).sum() / w.sum(),
               z=(kb - nb * p0).sum() / math.sqrt(var), base_wdiff=np.nan, base_z=np.nan)
    base = keys[:2] if keys[:2] == ["lc8", "gc"] else keys[:1]
    if baseline and len(base) < len(keys):
        sub = x[x.set_index(keys).index.isin(p0.index)]
        r = mh(sub, col, a, b, base, 1, baseline=False)
        out["base_wdiff"], out["base_z"] = r["wdiff"], r["z"]
    return out


def consistency(df, keys, min_point=50, groups=("mass", "mzd")):
    """chi2, ndf, cells, LJs, kept and extra of one no-vertex fraction shared by all groups
    (default: the mass points) in every cell of keys, groups with >= min_point LJs in the cell
    (table 22 of pair_content.py)."""
    tot, ndf, ncell, nlj, s1, s2 = 0.0, 0, 0, 0, 0.0, 0.0
    agg = df.groupby(keys + list(groups), observed=True).no_vtx.agg(["sum", "size"])
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
        tot += sum((x - y * p0) ** 2 / (y * p0 * (1 - p0)) for x, y in pts)
        ndf += len(pts) - 1
        ncell += 1
        nlj += n
        s1 += sum((x - y * p0) ** 2 - y * p0 * (1 - p0) * (1 - y / n) for x, y in pts)
        s2 += sum(y ** 2 * (1 - y / n) for x, y in pts)
    extra = math.sqrt(max(s1, 0.0) / s2) if s2 > 0 else float("nan")
    return dict(cells=ncell, n_lj=nlj, kept=nlj / len(df) if len(df) else float("nan"),
                chi2=tot, ndf=ndf, extra=extra)


def standardised(df, keys, ref):
    """Per group of keys: observed no-vertex fraction and the fraction standardised to the
    distance composition ref (share per distance bin): sum_b ref_b f(group, b), with the
    share of ref on bins empty in the group (uncovered)."""
    rows = []
    for g, x in df.groupby(keys, observed=True):
        fb = x.groupby("dbin", observed=True).no_vtx.mean()
        cov = ref.index.isin(fb.index)
        f_std = (ref[cov] * fb.reindex(ref.index[cov])).sum() / ref[cov].sum()
        rows.append(dict(zip(keys, g if isinstance(g, tuple) else (g,)), n=len(x),
                         k=int(x.no_vtx.sum()), f=x.no_vtx.mean(), f_std=f_std,
                         uncovered=float(ref[~cov].sum())))
    return pd.DataFrame(rows)


def derive(df, cls):
    """Join the analyze.py columns (lxy, gen_dr, cls) and derive the cell variables."""
    # event numbers repeat across the files of a sample: the leading muons complete the key
    k = ["sample", "run", "lumi", "event", "m0_pt", "m0_eta", "m1_pt"]
    c = cls[k + ["lxy", "gen_dr", "cls"]].copy()
    m = df[df.gen_matched].merge(c, on=k, how="left", validate="one_to_one")
    m["no_vtx"] = (m.vtx_chi2 < 0).astype(int)
    m["content"] = np.select([(m.n_pf == 2) & (m.n_dsa == 0), (m.n_pf == 1) & (m.n_dsa == 1),
                              (m.n_pf == 0) & (m.n_dsa == 2)], ["PF-PF", "PF-DSA", "DSA-DSA"],
                             "other")
    m["opp"] = m.m0_charge != m.m1_charge
    aeta = ((m.d1_eta + m.d2_eta) / 2).abs()
    m["region"] = np.where(aeta < 0.9, "barrel", np.where(aeta > 1.2, "endcap", "overlap"))
    m["dbin"] = pd.cut(m.emu, DIST_EDGES, right=False, labels=False)
    m["dc"] = pd.cut(m.emu, DIST_COARSE, right=False, labels=False)
    m["lslice"] = pd.cut(m.lxy, LXY_DSA, right=False, labels=False)
    m["lc8"] = pd.cut(m.lxy, split(LXY_CONTENT, 8), right=False, labels=False)
    m["lc4"] = pd.cut(m.lxy, split(LXY_CONTENT, 4), right=False, labels=False)
    m["gc"] = pd.cut(m.gen_dr, GEN_EDGES, right=False, labels=False)
    m["gcf"] = pd.cut(m.gen_dr, GEN_EDGES_FINE, right=False, labels=False)
    m["gw"] = np.where((m.gen_dr >= 0.1) & (m.gen_dr < 0.2), "narrow",
                       np.where((m.gen_dr >= 0.2) & (m.gen_dr < 0.5), "wide", "other"))
    m["m0_ptrel"] = m.m0_ptErr / m.m0_pt
    m["m1_ptrel"] = m.m1_ptErr / m.m1_pt
    for name, f, how in QUALITY:
        a, b = m[f"m0_{f}"], m[f"m1_{f}"]
        m["q_" + name] = np.maximum(a, b) if how == "max" else np.minimum(a, b)
    return m


def qbins(m, pop):
    """Quintile bins of each quality variable over the population pop (index of m)."""
    for name, _, _ in QUALITY:
        v = m.loc[pop, "q_" + name]
        edges = np.unique(np.quantile(v.dropna(), [0, 0.2, 0.4, 0.6, 0.8, 1.0]))
        edges[0], edges[-1] = -np.inf, np.inf
        m["b_" + name] = pd.cut(m["q_" + name], edges, right=False, labels=False)


def status_checks(nm, r):
    """Emulated against stored TTMD status for leading DSA-DSA pairs, and the pairs that the
    emulation puts above 15 cm."""
    dd = r[(r.m0_type == 1) & (r.m1_type == 1) & (r.emu_status >= 0)]
    w = dd[dd.pv_n > 0]
    nv = dd[dd.pv_n == 0]
    say(f"status {nm} DSA-DSA pairs with a stored vertex {len(w)}: stored dcaStatus 1 "
        f"{int((w.pv_dcaStatus == 1).sum())}, emulated status 1 {int((w.emu_status == 1).sum())}, "
        f"emulated > 15 cm {int((w.emu > 15).sum())}; without a stored vertex {len(nv)}: "
        f"emulated status 0 {int((nv.emu_status == 0).sum())}, status 1 and > 15 cm "
        f"{int((nv.emu > 15).sum())}, status 1 and <= 15 cm {int((nv.emu <= 15).sum())}")
    say(f"path {nm} DSA-DSA (1 Newton from the states, 2 Newton qual 0.1, 3 Newton from the "
        "transverse points, 4 transverse points kept): "
        + str(dd.emu_path.value_counts().sort_index().to_dict()))


# ---------------------------------------------------------------- helices: angle, global minimum
def helices(A, B, bz_tesla):
    """Helix parameters of A and B as TTMD's helix-helix step sets them up: each helix with Bz
    at the other track's reference point. s is the forward transverse arc (cm)."""
    out = []
    for T, O in ((A, B), (B, A)):
        b = bz_tesla(O["x"], O["y"], O["z"]) * C_GEV
        pt = np.hypot(T["px"], T["py"])
        g = pt / (T["q"] * b)
        s0, c0 = -T["py"] / pt, -T["px"] / pt
        out.append(dict(x=T["x"], y=T["y"], z=T["z"], g=g, s0=s0, c0=c0,
                        p0=np.arctan2(s0, c0), tz=T["pz"] / pt))
    return out[0], out[1]


def _bc(P, nd):
    """Helix parameters broadcast over nd trailing axes."""
    return {k: v.reshape(v.shape + (1,) * nd) for k, v in P.items()}


def hpos(P, s):
    """Position (x, y, z) on helix P at arc s (s of shape (n, ...))."""
    Q = _bc(P, s.ndim - 1)
    ph = Q["p0"] - s / Q["g"]
    return (Q["x"] + Q["g"] * (np.sin(ph) - Q["s0"]),
            Q["y"] + Q["g"] * (-np.cos(ph) + Q["c0"]), Q["z"] + s * Q["tz"])


def _tangent(P, X, Y):
    """Unit tangent of helix P at the projection of the transverse point (X, Y) on its circle."""
    xc, yc = P["x"] - P["g"] * P["s0"], P["y"] + P["g"] * P["c0"]
    ph = np.arctan2((X - xc) / P["g"], -(Y - yc) / P["g"])
    t = np.stack([-np.cos(ph), -np.sin(ph), P["tz"]])
    return t / np.linalg.norm(t, axis=0)


def angle_at(A, B, cross, bz_tesla):
    """3D opening angle (rad) of the two helices at the point cross (3, n)."""
    PA, PB = helices(A, B, bz_tesla)
    ta, tb = _tangent(PA, cross[0], cross[1]), _tangent(PB, cross[0], cross[1])
    return np.arctan2(np.linalg.norm(np.cross(ta.T, tb.T), axis=1), np.sum(ta * tb, axis=0))


def _inside(x, y, z, s):
    return (x * x + y * y < R_MAX ** 2) & (np.abs(z) < Z_MAX) & (np.abs(s) <= S_MAX)


def _f2(PA, PB, sa, sb):
    """Squared distance of the helix points at sa and sb; inf when a point is out of range."""
    xa, ya, za = hpos(PA, sa)
    xb, yb, zb = hpos(PB, sb)
    f = (xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2
    return np.where(_inside(xa, ya, za, sa) & _inside(xb, yb, zb, sb), f, np.inf)


def _newton(PA, PB, sa, sb, iters=10, move="ab"):
    """Damped Newton steps on the squared distance, a step kept only when it lowers it inside
    the range; move "ab" both arcs, "a" or "b" one (profiles). sa, sb of shape (n, K)."""
    QA, QB = _bc(PA, 1), _bc(PB, 1)
    f = _f2(PA, PB, sa, sb)
    for _ in range(iters):
        xa, ya, za = hpos(PA, sa)
        xb, yb, zb = hpos(PB, sb)
        D = np.stack([xa - xb, ya - yb, za - zb])
        pha, phb = QA["p0"] - sa / QA["g"], QB["p0"] - sb / QB["g"]
        TA = np.stack([-np.cos(pha), -np.sin(pha), np.broadcast_to(QA["tz"], sa.shape)])
        TB = np.stack([-np.cos(phb), -np.sin(phb), np.broadcast_to(QB["tz"], sb.shape)])
        KA = np.stack([-np.sin(pha) / QA["g"], np.cos(pha) / QA["g"], 0 * sa])
        KB = np.stack([-np.sin(phb) / QB["g"], np.cos(phb) / QB["g"], 0 * sb])
        ga, gb = 2 * np.sum(D * TA, 0), -2 * np.sum(D * TB, 0)
        haa = 2 * (np.sum(TA * TA, 0) + np.sum(D * KA, 0))
        hbb = 2 * (np.sum(TB * TB, 0) - np.sum(D * KB, 0))
        hab = -2 * np.sum(TA * TB, 0)
        det = haa * hbb - hab * hab
        with np.errstate(all="ignore"):
            if move == "a":
                da, db, good = -ga / haa, 0 * gb, haa > 0
            elif move == "b":
                da, db, good = 0 * ga, -gb / hbb, hbb > 0
            else:
                da, db = -(hbb * ga - hab * gb) / det, -(haa * gb - hab * ga) / det
                good = (det > 0) & (haa > 0)
            good &= np.isfinite(da) & np.isfinite(db) & np.isfinite(f)
            scale = np.minimum(1.0, 50.0 / np.maximum(np.hypot(da, db), 1e-300))
        da, db = np.where(good, da * scale, 0), np.where(good, db * scale, 0)
        moved = np.zeros(sa.shape, bool)
        for lam in (1.0, 0.5, 0.25, 0.125):
            na, nb = sa + lam * da, sb + lam * db
            fn = _f2(PA, PB, na, nb)
            acc = ~moved & good & (fn < f)
            sa, sb, f = np.where(acc, na, sa), np.where(acc, nb, sb), np.where(acc, fn, f)
            moved |= acc
    return sa, sb, f


def _grid(pa, pb, ca, cb, half, step, K):
    """Grid of step `step` over [c - half, c + half] in both arcs around per-pair centres; the
    K smallest grid local minima (8 neighbours): sa, sb, f of shape (n, K), and the grid."""
    off = np.arange(-half, half + 0.5 * step, step)
    m, k = len(off), len(ca)
    ga, gb = ca[:, None] + off, cb[:, None] + off
    xa, ya, za = hpos(pa, ga)
    xb, yb, zb = hpos(pb, gb)
    f = ((xa[:, :, None] - xb[:, None, :]) ** 2 + (ya[:, :, None] - yb[:, None, :]) ** 2
         + (za[:, :, None] - zb[:, None, :]) ** 2)
    f = np.where(_inside(xa, ya, za, ga)[:, :, None] & _inside(xb, yb, zb, gb)[:, None, :],
                 f, np.inf)
    fp = np.pad(f, ((0, 0), (1, 1), (1, 1)), constant_values=np.inf)
    loc = np.isfinite(f)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                loc &= f <= fp[:, 1 + di:1 + di + m, 1 + dj:1 + dj + m]
    cand = np.where(loc, f, np.inf).reshape(k, -1)
    top = np.argpartition(cand, min(K, cand.shape[1] - 1), axis=1)[:, :K]
    return (np.take_along_axis(ga, top // m, 1), np.take_along_axis(gb, top % m, 1),
            np.take_along_axis(cand, top, 1), f, off)


def _profile(pa, pb, f, off, which, step=5.0, K=6):
    """Candidates from the profile: one arc on a grid of `step`, the other minimised by Newton
    steps from the coarse-grid row minimum; the K smallest local minima of the profile."""
    k, m = f.shape[0], len(off)
    fr = f if which == "a" else np.swapaxes(f, 1, 2)
    jbest = np.argmin(fr, axis=2)
    fine = np.arange(off[0], off[-1] + 0.5 * step, step)
    r = np.clip(np.rint((fine - off[0]) / (off[1] - off[0])).astype(int), 0, m - 1)
    other = off[jbest[:, r]] + (fine - off[r])
    s1 = np.broadcast_to(fine, (k, len(fine))).copy()
    sa, sb = (s1, other) if which == "a" else (other, s1)
    ok0 = np.isfinite(np.min(fr, axis=2))[:, r]
    sa, sb, F = _newton(pa, pb, sa, sb, iters=8, move="b" if which == "a" else "a")
    F = np.where(ok0, F, np.inf)
    Fp = np.pad(F, ((0, 0), (1, 1)), constant_values=np.inf)
    cand = np.where(np.isfinite(F) & (F <= Fp[:, :-2]) & (F <= Fp[:, 2:]), F, np.inf)
    top = np.argpartition(cand, min(K, cand.shape[1] - 1), axis=1)[:, :K]
    return (np.take_along_axis(sa, top, 1), np.take_along_axis(sb, top, 1),
            np.take_along_axis(cand, top, 1))


def _refine(pa, pb, sa, sb, fc, h):
    """Pattern search (12 levels of step h / 2^l, 5 x 5 offsets, two passes each), then Newton."""
    offs = np.array([(i, j) for i in range(-2, 3) for j in range(-2, 3)], float)
    for lev in range(12):
        for _ in range(2):
            na = sa[:, :, None] + offs[:, 0] * h / 2 ** lev
            nb = sb[:, :, None] + offs[:, 1] * h / 2 ** lev
            fn = _f2(pa, pb, na, nb)
            b = np.argmin(fn, axis=2)[:, :, None]
            fb = np.take_along_axis(fn, b, 2)[:, :, 0]
            acc = np.isfinite(fc) & (fb < fc)
            sa = np.where(acc, np.take_along_axis(na, b, 2)[:, :, 0], sa)
            sb = np.where(acc, np.take_along_axis(nb, b, 2)[:, :, 0], sb)
            fc = np.where(acc, fb, fc)
    return _newton(pa, pb, sa, sb)


def global_min(A, B, bz_tesla, h=20.0, K=10, fine=2.0, Kf=4, nfine=3, Kp=6, chunk=100,
               exhaustive=False):
    """Global minimum distance of the two helices (helices()), both closest points inside
    r < R_MAX, |z| < Z_MAX, |s| <= S_MAX. Candidates: the K smallest local minima of a grid of
    step h over the whole range, the Kf smallest of a grid of step `fine` over +-60 cm around
    each of the nfine best of them, the Kp smallest local minima of the two profiles; each
    refined by _refine. exhaustive: a single grid of step `fine` over the whole range (slow,
    for gmin_stability). Returns (distance, sA, sB), NaN when no point is in range."""
    PA, PB = helices(A, B, bz_tesla)
    n = len(A["x"])
    out = np.full((3, n), np.nan)
    for c0 in range(0, n, chunk):
        idx = np.arange(c0, min(n, c0 + chunk))
        pa, pb = _sub(PA, idx), _sub(PB, idx)
        z0 = np.zeros(len(idx))
        if exhaustive:
            sa, sb, fc, _, _ = _grid(pa, pb, z0, z0, S_MAX, fine, K)
        else:
            sa, sb, fc, fg, off = _grid(pa, pb, z0, z0, S_MAX, h, K)
            parts = [(sa, sb, fc), _profile(pa, pb, fg, off, "a", K=Kp),
                     _profile(pa, pb, fg, off, "b", K=Kp)]
            order = np.argsort(fc, axis=1)
            for r in range(min(nfine, K)):
                o = order[:, r:r + 1]
                ok = np.isfinite(np.take_along_axis(fc, o, 1)[:, 0])
                ca = np.where(ok, np.take_along_axis(sa, o, 1)[:, 0], 0.0)
                cb = np.where(ok, np.take_along_axis(sb, o, 1)[:, 0], 0.0)
                fa, fb, ff, _, _ = _grid(pa, pb, ca, cb, 60.0, fine, Kf)
                parts.append((fa, fb, np.where(ok[:, None], ff, np.inf)))
            sa, sb, fc = (np.concatenate([p[i] for p in parts], 1) for i in range(3))
        sa, sb, fc = _refine(pa, pb, sa, sb, fc, fine)
        best = np.argmin(fc, axis=1)[:, None]
        fbest = np.take_along_axis(fc, best, 1)[:, 0]
        out[0, idx] = np.where(np.isfinite(fbest), np.sqrt(fbest), np.nan)
        out[1, idx] = np.take_along_axis(sa, best, 1)[:, 0]
        out[2, idx] = np.take_along_axis(sb, best, 1)[:, 0]
    return out


def gmin_stability(r, rows_var, rows_exh):
    """Global minimum of the LJs rows_var with a finer variant (10 cm grid, 16 grid minima, four
    1 cm grids), and of rows_exh with an exhaustive 1 cm grid, against the default (column
    gmin of r)."""
    A, B, _, _ = ordered_states(r)
    res = {}
    for name, rows, kw in (("variant", rows_var, dict(h=10.0, K=16, fine=1.0, Kf=6, nfine=4)),
                           ("exhaustive", rows_exh, dict(fine=1.0, K=6, chunk=1,
                                                         exhaustive=True))):
        d = global_min(_sub(A, rows), _sub(B, rows), bz_parabolic, **kw)[0]
        res[name] = d - r.gmin.to_numpy()[rows]
    return res


# ---------------------------------------------------------------- per-pair columns
def ordered_states(df):
    """States of the two leading muons in the producer's order (PF before DSA, then by index),
    the rows where both are available, and the leg kinds."""
    s0, k0 = leg(df, 0)
    s1, k1 = leg(df, 1)
    t0, t1 = df.m0_type.to_numpy(), df.m1_type.to_numpy()
    swap = (t0 > t1) | ((t0 == t1) & (df.m0_idx.to_numpy() > df.m1_idx.to_numpy()))
    A = {k: np.where(swap, s1[k], s0[k]) for k in s0}
    B = {k: np.where(swap, s0[k], s1[k]) for k in s0}
    ok = np.all(np.isfinite(np.vstack(list(A.values()) + list(B.values()))), axis=0)
    ok &= df.n_mu.to_numpy() >= 2
    kinds = np.char.add(np.char.add(np.where(swap, k1, k0).astype(str), "+"),
                        np.where(swap, k0, k1).astype(str))
    return A, B, ok, kinds


def emulate(df, bz_tesla):
    """Emulated TTMD of the two leading muons in the producer's order: (status: 1, 0, or -1
    when a leg lacks a track; distance, NaN unless status 1; path; leg kinds; crossing point
    (3, n))."""
    A, B, ok, kinds = ordered_states(df)
    status = np.full(len(df), -1)
    dist = np.full(len(df), np.nan)
    path = np.zeros(len(df), int)
    cross = np.full((3, len(df)), np.nan)
    idx = np.where(ok)[0]
    s, d, c, p = ttmd(_sub(A, idx), _sub(B, idx), bz_tesla)
    status[idx], dist[idx], path[idx] = s.astype(int), np.where(s, d, np.nan), p
    cross[:, idx] = np.where(s, c, np.nan)
    return status, dist, path, kinds, cross


def pair_columns(r):
    """Emulated TTMD (parabolic and uniform field), crossing point, and for leading DSA-DSA
    pairs the opening angle at the crossing, L and the global minimum distance."""
    r["emu_status"], r["emu"], r["emu_path"], r["legs"], cross = emulate(r, bz_parabolic)
    r["emu_u"] = emulate(r, bz_uniform)[1]
    r["cx"], r["cy"], r["cz"] = cross
    r["cr"] = np.hypot(r.cx, r.cy)
    dd = ((r.m0_type == 1) & (r.m1_type == 1) & (r.emu_status >= 0)).to_numpy()
    A, B, _, _ = ordered_states(r)
    idx = np.where(dd)[0]
    As, Bs = _sub(A, idx), _sub(B, idx)
    ang, gm = np.full(len(r), np.nan), np.full(len(r), np.nan)
    with np.errstate(invalid="ignore"):
        ang[idx] = angle_at(As, Bs, cross[:, idx], bz_parabolic)
    g = global_min(As, Bs, bz_parabolic)
    gm[idx] = g[0]
    r["ang"], r["gmin"] = ang, gm
    PA, PB = helices(As, Bs, bz_parabolic)
    pa, pb = hpos(PA, g[1][:, None]), hpos(PB, g[2][:, None])
    for i, c in enumerate(("gx", "gy", "gz")):
        v = np.full(len(r), np.nan)
        v[idx] = 0.5 * (pa[i][:, 0] + pb[i][:, 0])
        r[c] = v
    dphi = (r.m0_phi - r.m1_phi + np.pi) % (2 * np.pi) - np.pi
    r["dr_ref"] = np.hypot(r.m0_eta - r.m1_eta, dphi)
    s = 0.5 * np.sqrt(r.m0_dxyBSErr ** 2 + r.m1_dxyBSErr ** 2 + r.m0_dzBSErr ** 2
                      + r.m1_dzBSErr ** 2)
    r["L"] = s / r.ang
    sin_theta_pair(r)


def sin_theta_pair(r):
    """sth: sin of the polar angle of the summed momentum of the two leading tracks at the crossing
    point (where the angle exists); aeff = angle / sth."""
    A, B, _, _ = ordered_states(r)
    idx = np.where(np.isfinite(r.ang.to_numpy()))[0]
    As, Bs = _sub(A, idx), _sub(B, idx)
    PA, PB = helices(As, Bs, bz_parabolic)
    cx, cy = r.cx.to_numpy()[idx], r.cy.to_numpy()[idx]
    p = (_tangent(PA, cx, cy) * np.sqrt(As["px"] ** 2 + As["py"] ** 2 + As["pz"] ** 2)
         + _tangent(PB, cx, cy) * np.sqrt(Bs["px"] ** 2 + Bs["py"] ** 2 + Bs["pz"] ** 2))
    sth = np.full(len(r), np.nan)
    sth[idx] = np.hypot(p[0], p[1]) / np.linalg.norm(p, axis=0)
    r["sth"] = sth
    r["aeff"] = r.ang / r.sth


# ---------------------------------------------------------------- the two conditions
def all_dd(r):
    """All LJs whose two leading muons are DSA with both tracks available, with no_vtx
    (pv_n == 0) and the bins of distance and angle."""
    x = r[(r.m0_type == 1) & (r.m1_type == 1) & (r.emu_status >= 0)]
    return x.assign(no_vtx=(x.pv_n == 0).astype(int),
                    ab=pd.cut(x.ang, ANG_EDGES, right=False, labels=False),
                    dbin=pd.cut(x.emu, DIST_EDGES, right=False, labels=False),
                    dc=pd.cut(x.emu, DIST_COARSE, right=False, labels=False))


def le15(x):
    """The LJs without the distance cut (distance <= 15 cm or TTMD failed)."""
    return x[~((x.emu_status == 1) & (x.emu > 15))]


def decompose(x):
    """n, k, f, the cut LJs (n_cut, k_cut with no vertex), s15, and f_inv of the others."""
    cut = ((x.emu_status == 1) & (x.emu > 15)).to_numpy()
    nv = x.no_vtx.to_numpy()
    n_le = int((~cut).sum())
    return dict(n=len(x), k=int(nv.sum()), f=nv.mean(), n_cut=int(cut.sum()),
                k_cut=int(nv[cut].sum()), s15=cut.mean(), n_le15=n_le,
                k_le15=int(nv[~cut].sum()), f_inv=nv[~cut].mean() if n_le else np.nan)


def decomposition(label, detail, xa, xb, la, lb):
    """Row of decomposition_cut_vs_invalid_fit.csv: group b against group a."""
    a, b = decompose(xa), decompose(xb)
    d = b["f"] - a["f"]
    ts = (1 - 0.5 * (a["f_inv"] + b["f_inv"])) * (b["s15"] - a["s15"])
    ti = (1 - 0.5 * (a["s15"] + b["s15"])) * (b["f_inv"] - a["f_inv"])
    row = dict(comparison=label, detail=detail, group_a=la, group_b=lb)
    row.update({f"{k}_a": v for k, v in a.items()})
    row.update({f"{k}_b": v for k, v in b.items()})
    row.update(diff=d, term_s15=ts, term_finv=ti, residual=d - ts - ti,
               share_s15=ts / d if d else np.nan, share_finv=ti / d if d else np.nan)
    return row


EXTRAS = [("none", []), ("distance", ["dbin"]), ("angle", ["ab"]),
          ("angle x coarse distance", ["ab", "dc"]), ("angle_eff", ["eb"]),
          ("angle_eff x coarse distance", ["eb", "dc"])]
EXTRAS_Q = EXTRAS + [("quality " + q[0], ["b_" + q[0]]) for q in QUALITY]


def deviance(x, keys):
    """Binomial deviance of the per-cell fractions of keys, and the number of cells."""
    g = x.groupby(keys, observed=True).no_vtx.agg(["sum", "size"])
    k, n = g["sum"].to_numpy(float), g["size"].to_numpy(float)
    f = k / n
    with np.errstate(divide="ignore", invalid="ignore"):
        ll = np.where(k > 0, k * np.log(f), 0.0) + np.where(n > k, (n - k) * np.log(1 - f), 0.0)
    return -2.0 * float(ll.sum()), len(g)


def by_bins(x, var, edges, label):
    """n, k, f of x per bin of var (or per value when edges is None)."""
    v = pd.cut(x[var], edges, right=False) if edges is not None else x[var]
    g = x.groupby(v, observed=True).no_vtx.agg(n="size", k="sum").reset_index()
    g["f"] = g.k / g.n
    g["bin"] = g[var].astype(str)
    return g.drop(columns=var).assign(variable=label)[["variable", "bin", "n", "k", "f"]]


def ctau(sample):
    """Proper decay length in mm from the sample name."""
    return sample.str.extract(r"GeV_([0-9p]+)mm")[0].str.replace("p", ".").astype(float)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    for a in ("rec_2mu2e", "rec_4mu", "rec_wp", "cls_2mu2e", "cls_4mu", "cls_wp", "out_dir"):
        ap.add_argument(a)
    a = ap.parse_args()
    out = a.out_dir
    rec = dict(zip(SETS, (a.rec_2mu2e, a.rec_wp, a.rec_4mu)))
    cls = dict(zip(SETS, (a.cls_2mu2e, a.cls_wp, a.cls_4mu)))

    # ---------------------------------------------------------------- records and checks
    R, M, P = {}, {}, {}
    for nm in SETS:
        r, nf = load(rec[nm])
        say(f"files {nm}: {nf}; LJs {len(r)}")
        two = r.n_mu == 2
        say(f"pv check {nm}: two-muon LJs {int(two.sum())}; (pv_n > 0) == (vtx_chi2 >= 0) in "
            f"{int(((r.pv_n > 0) == (r.vtx_chi2 >= 0))[two].sum())}; pv_n > 1 in "
            f"{int((r.pv_n > 1).sum())}; stored isValid != 1 in "
            f"{int(((r.pv_n > 0) & (r.pv_isValid != 1)).sum())}")
        pair_columns(r)
        status_checks(nm, r)
        R[nm] = r
        M[nm] = derive(r, pd.read_parquet(cls[nm]))
        say(f"join {nm}: gen-matched {len(M[nm])}, without analyze.py columns "
            f"{int(M[nm].cls.isna().sum())}")
        M[nm]["ab"] = pd.cut(M[nm].ang, ANG_EDGES, right=False, labels=False)
        M[nm]["lb"] = pd.cut(M[nm].L, L_EDGES, right=False, labels=False)
        M[nm]["eb"] = pd.cut(M[nm].aeff, AEFF_EDGES, right=False, labels=False)
        M[nm]["point"] = M[nm].mass.astype(str) + "," + M[nm].mzd.astype(str)
        M[nm]["one"] = 0
        P[nm] = M[nm][(M[nm].content == "DSA-DSA") & M[nm].opp].copy()
        say(f"population P({nm}): {len(P[nm])}, no vertex {int(P[nm].no_vtx.sum())}, "
            f"emulated distance missing {int(P[nm].emu.isna().sum())}")
    both = sorted(set(R["2Mu2E"]["sample"]) & set(R["2Mu2E wide-pair"]["sample"]))
    say(f"samples in both 2Mu2E record sets: {both}; LJs of P in them (counted in both sets in "
        f"universality (iv) and the pooled deviance, not de-duplicated): " + ", ".join(
            f"{nm} {int(P[nm]['sample'].isin(both).sum())} (distance <= 15 cm "
            f"{int(le15(P[nm])['sample'].isin(both).sum())})" for nm in SETS[:2]))
    m2, m4, p = M["2Mu2E"], M["4Mu"], P["2Mu2E"]
    t15 = m2[m2.content != "other"].groupby(["content", "opp"]).no_vtx.agg(["size", "sum"])
    t15.to_csv(f"{out}/table15_check.csv")
    say("table 15 opposite-charge n: " + ", ".join(
        f"{c} {int(t15.loc[(c, True), 'size'])}" for c in ("PF-PF", "PF-DSA", "DSA-DSA")))
    dd_all = m2[m2.content == "DSA-DSA"].copy()
    qbins(m2, p.index)
    for name, _, _ in QUALITY:
        p["b_" + name] = m2.loc[p.index, "b_" + name]
        dd_all["b_" + name] = m2.loc[dd_all.index, "b_" + name]
    base = m4[(m4.lxy >= 100) & (m4.gw != "other") & m4.cls.isin(["A", "D"])]
    q = base[base.content == "DSA-DSA"].copy()
    say("T6 population: class A/D LJs at Lxy >= 100, narrow+wide " + str(len(base))
        + "; two-muon DSA-DSA " + str(len(q)) + "; by class and bin "
        + str(q.groupby(["cls", "gw"]).size().to_dict()))
    qbins(m4, q.index)
    for name, _, _ in QUALITY:
        q["b_" + name] = m4.loc[q.index, "b_" + name]

    validation_table(R, p, dd_all, m4, q, out)
    t1_t2(R, P, out)
    pr = p.assign(rad=np.where(p.lxy < 60, "inner", np.where((p.lxy >= 150) & (p.lxy < 300),
                                                              "outer", "mid")))
    pb = p[p.lxy < 200]
    y_opp = p[(p.mzd == 5.0) & p.mass.isin([500, 1000])]
    y_all = dd_all[(dd_all.mzd == 5.0) & dd_all.mass.isin([500, 1000])]
    decomposition_table(y_opp, y_all, pr, pb, q, out)
    tests_table(y_opp, y_all, pr, pb, q, out)
    ref = p.dbin.value_counts(normalize=True).sort_index()
    t4s = pd.concat([standardised(p.assign(region="all"), ["region", "lslice"], ref),
                     standardised(p, ["region", "lslice"], ref)])
    t4s["lxy_cm"] = t4s.lslice.map(lambda i: f"[{LXY_DSA[int(i)]:g},{LXY_DSA[int(i) + 1]:g})")
    t4s.to_csv(f"{out}/T4_standardised.csv", index=False)
    say("T4 standardised (barrel, endcap, all)\n" + t4s[t4s.region.isin(["all", "barrel", "endcap"])]
        [["region", "lxy_cm", "n", "k", "f", "f_std", "uncovered"]].round(3).to_string(index=False))
    universality_table(P, out)
    candidates_tables(R, P, out)
    deviance_table(P, out)
    angle_standardised_table(P, out)
    sin_theta_table(P, out)
    edge_table(R, out)
    t7_tables(R, q, out)
    t7_location_table(R, out)
    with open(f"{out}/summary.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(SUMMARY) + "\n")


def validation_table(R, p, dd_all, m4, q, out):
    """validation_dca.csv: per group (set, content, charge, region or mass point), emulated
    against stored dca."""
    groups = []
    for fs in ("2Mu2E", "4Mu"):
        r = R[fs]
        pair = np.select([(r.m0_type == 1) & (r.m1_type == 1), (r.m0_type + r.m1_type) == 1],
                         ["DSA-DSA", "PF-DSA"], "PF-PF")
        exact = ~r.legs.str.contains("pf_global") & ~r.legs.str.contains("pf_sa")
        for ct in ("DSA-DSA", "PF-DSA", "PF-PF"):
            groups.append((f"{fs} all LJs, leading pair {ct}", r[pair == ct]))
        for ct in ("PF-DSA", "PF-PF"):
            groups.append((f"{fs} all LJs, leading pair {ct}, PF legs tracker-only",
                           r[(pair == ct) & exact]))
    m4o = m4[(m4.content == "DSA-DSA") & m4.opp]
    groups += [("2Mu2E gen-matched DSA-DSA opposite", p),
               ("2Mu2E gen-matched DSA-DSA both charges", dd_all),
               ("2Mu2E gen-matched DSA-DSA opposite, barrel", p[p.region == "barrel"]),
               ("2Mu2E gen-matched DSA-DSA opposite, endcap", p[p.region == "endcap"]),
               ("2Mu2E gen-matched DSA-DSA opposite, (500, 5)", p[(p.mass == 500) & (p.mzd == 5.0)]),
               ("2Mu2E gen-matched DSA-DSA opposite, (1000, 5)", p[(p.mass == 1000) & (p.mzd == 5.0)]),
               ("4Mu gen-matched DSA-DSA opposite, barrel", m4o[m4o.region == "barrel"]),
               ("4Mu gen-matched DSA-DSA opposite, endcap", m4o[m4o.region == "endcap"]),
               ("4Mu class A narrow, all contents", m4[(m4.cls == "A") & (m4.gw == "narrow")]),
               ("4Mu class A wide, all contents", m4[(m4.cls == "A") & (m4.gw == "wide")])]
    for c in ("A", "D"):
        for gw in ("narrow", "wide"):
            groups.append((f"T6 population, class {c} {gw}", q[(q.cls == c) & (q.gw == gw)]))
    r = R["2Mu2E wide-pair"]
    groups.append(("2Mu2E wide-pair all LJs, leading pair DSA-DSA",
                   r[(r.m0_type == 1) & (r.m1_type == 1)]))
    val = validation(groups)
    val.to_csv(f"{out}/validation_dca.csv", index=False)
    vb = val[(val.dca_bin != "all") & (val.n >= 20)]
    worst = vb.groupby(["group", "field"], sort=False).agg(
        worst_low=("low_gt1", "max"), worst_high=("high_gt1", "max"))
    va = val[val.dca_bin == "all"].set_index(["group", "field"])[
        ["n", "within_0p1", "within_10pc", "low_gt1", "high_gt1", "failed"]]
    say("emulated against stored dca, pairs with a stored vertex (all bins; worst stored-dca "
        "bin among those with >= 20 pairs)\n" + va.join(worst).round(3).to_string())


def t1_t2(R, P, out):
    """T1 (P(2Mu2E)), T2 per sample set, and the same-mass-point table."""
    p = P["2Mu2E"]
    t1 = pd.concat([rate_table(p.assign(point="all"), ["point", "dbin"]),
                    rate_table(p, ["point", "dbin"])])
    t1["dist_cm"] = t1.dbin.map(lambda i: f"[{DIST_EDGES[int(i)]:g},{DIST_EDGES[int(i) + 1]:g})")
    t1.to_csv(f"{out}/T1_novtx_vs_distance.csv", index=False)
    say("T1 (all points)\n" + t1[t1.point == "all"][["dist_cm", "n", "k", "f", "err"]]
        .round(4).to_string(index=False))
    sh = p.groupby("point").dbin.value_counts(normalize=True).unstack().round(3)
    sh.to_csv(f"{out}/T1_distance_shares_by_point.csv")
    rows = []
    for nm in SETS:
        nv = P[nm][P[nm].no_vtx == 1].assign(gt15=lambda x: (x.emu > 15).astype(int))
        for var in ("lslice", "region", "point"):
            g = nv.groupby(var, observed=True).gt15.agg(["size", "sum"]).reset_index()
            for _, rr in g.iterrows():
                f, e = frac(rr["sum"], rr["size"])
                rows.append(dict(set=nm, by=var, value=rr[var], n_novtx=int(rr["size"]),
                                 n_gt15=int(rr["sum"]), share_gt15=f, err=e))
        f, e = frac(nv.gt15.sum(), len(nv))
        rows.append(dict(set=nm, by="all", value="all", n_novtx=len(nv),
                         n_gt15=int(nv.gt15.sum()), share_gt15=f, err=e))
    t2 = pd.DataFrame(rows)
    t2.to_csv(f"{out}/T2_condition_share.csv", index=False)
    say("T2 (all)\n" + t2[t2.by == "all"].round(4).to_string(index=False))
    rows = []
    for nm in SETS:
        r = R[nm]
        dd = r[(r.m0_type == 1) & (r.m1_type == 1) & (r.emu_status >= 0)]
        dd = dd.assign(nv=(dd.pv_n == 0).astype(int), ctau=ctau(dd["sample"]),
                       gt15=((dd.pv_n == 0) & (dd.emu > 15)).astype(int))
        g = dd.groupby(["mass", "mzd", "ctau"]).agg(n_lj=("nv", "size"), n_novtx=("nv", "sum"),
                                                    n_gt15=("gt15", "sum")).reset_index()
        g = pd.concat([g, pd.DataFrame([dict(mass="all", mzd="all", ctau="all",
                                             n_lj=len(dd), n_novtx=int(dd.nv.sum()),
                                             n_gt15=int(dd.gt15.sum()))])])
        rows.append(g.assign(set=nm))
    sm = pd.concat(rows)[["set", "mass", "mzd", "ctau", "n_lj", "n_novtx", "n_gt15"]]
    sm["share_gt15"] = sm.n_gt15 / sm.n_novtx
    sm["err"] = np.sqrt(sm.share_gt15 * (1 - sm.share_gt15) / sm.n_novtx)
    sm.to_csv(f"{out}/T2_same_mass_points.csv", index=False)
    say("T2 same mass points (all LJs)\n" + sm.round(3).to_string(index=False))


def decomposition_table(y_opp, y_all, pr, pb, q, out):
    """decomposition_cut_vs_invalid_fit.csv."""
    rows = []
    for lab, y in (("opposite", y_opp), ("both charges", y_all)):
        rows.append(decomposition("(a) (1000, 5) vs (500, 5)", lab, y[y.mass == 500],
                                  y[y.mass == 1000], "(500, 5)", "(1000, 5)"))
    rows.append(decomposition("(b) outer vs inner", "opposite", pr[pr.rad == "inner"],
                              pr[pr.rad == "outer"], "[0, 60) cm", "[150, 300) cm"))
    rows.append(decomposition("(c) endcap vs barrel", "Lxy < 200 cm pooled",
                              pb[pb.region == "barrel"], pb[pb.region == "endcap"], "barrel",
                              "endcap"))
    for i, x in pb.groupby("lslice"):
        rows.append(decomposition("(c) endcap vs barrel",
                                  f"Lxy [{LXY_DSA[int(i)]:g},{LXY_DSA[int(i) + 1]:g})",
                                  x[x.region == "barrel"], x[x.region == "endcap"], "barrel",
                                  "endcap"))
    for c in ("A", "D"):
        for lab, x in (("both charges", q[q.cls == c]), ("opposite", q[(q.cls == c) & q.opp])):
            rows.append(decomposition(f"(d) 4Mu class {c} wide vs narrow", lab,
                                      x[x.gw == "narrow"], x[x.gw == "wide"], "narrow", "wide"))
    t = pd.DataFrame(rows)
    t.to_csv(f"{out}/decomposition_cut_vs_invalid_fit.csv", index=False)
    say("decomposition (no-vertex fraction = s15 + (1 - s15) f_inv)\n" + t[
        ["comparison", "detail", "n_a", "s15_a", "f_inv_a", "f_a", "n_b", "s15_b", "f_inv_b",
         "f_b", "diff", "term_s15", "term_finv", "residual"]].round(4).to_string(index=False))


def tests_table(y_opp, y_all, pr, pb, q, out):
    """invalid_fit_tests.csv: the comparisons on f_inv (LJs with distance <= 15 cm)."""
    rows = []

    def add(test, charge, nmin, x, col, ab, keys, extras):
        for cl, ex in extras:
            rows.append(dict(test=test, charge=charge, nmin=nmin, cells_of=cl,
                             **mh(le15(x), col, ab[0], ab[1], keys + ex, nmin)))

    for lab, y in (("opposite", y_opp), ("both charges", y_all)):
        for nmin in (30, 5):
            add("(a) (1000, 5) vs (500, 5)", lab, nmin, y, "mass", (500, 1000), ["lc8", "gc"],
                EXTRAS_Q if nmin == 5 else EXTRAS)
    add("(b) outer vs inner", "opposite", 5, pr, "rad", ("inner", "outer"), ["one"], EXTRAS_Q)
    add("(c) endcap vs barrel", "opposite", 5, pb, "region", ("barrel", "endcap"), ["lslice"],
        EXTRAS_Q)
    for c in ("A", "D"):
        for lab, x in (("both charges", q[q.cls == c]), ("opposite", q[(q.cls == c) & q.opp])):
            add(f"(d) 4Mu class {c} wide vs narrow", lab, 5, x, "gw", ("narrow", "wide"),
                ["one"], EXTRAS_Q)
    t = pd.DataFrame(rows)
    t.to_csv(f"{out}/invalid_fit_tests.csv", index=False)
    say("invalid-fit tests (f_inv, distance <= 15 cm)\n" + t[t.cells_of.isin([e[0] for e in EXTRAS])]
        [["test", "charge", "nmin", "cells_of", "cells", "n_a", "n_b", "kept", "f_a", "f_b",
          "wdiff", "z", "base_wdiff", "base_z"]].round(4).to_string(index=False))


def universality_table(P, out):
    """invalid_fit_universality.csv."""
    p2 = le15(P["2Mu2E"])
    sets = pd.concat([le15(P[nm]).assign(set=nm) for nm in SETS])
    groupings = [("(i) mass points of the 2Mu2E set", p2, ["mass", "mzd"]),
                 ("(ii) barrel and endcap, 2Mu2E set", p2[p2.region.isin(["barrel", "endcap"])],
                  ["region"]),
                 ("(iii) Lxy slices, 2Mu2E set", p2, ["lslice"]),
                 ("(iv) the three sample sets", sets, ["set"])]
    rows = []
    for lab, x, g in groupings:
        for cl, keys in (("none", ["one"]), ("distance", ["dbin"]), ("angle", ["ab"]),
                         ("angle x distance", ["ab", "dbin"]), ("angle_eff", ["eb"]),
                         ("angle_eff x distance", ["eb", "dbin"])):
            for mp in (50, 10):
                rows.append(dict(groups=lab, cells_of=cl, min_lj_per_group=mp,
                                 **consistency(x, keys, mp, g)))
    t = pd.DataFrame(rows)
    t.to_csv(f"{out}/invalid_fit_universality.csv", index=False)
    say("invalid-fit universality\n" + t.round(3).to_string(index=False))


def candidates_tables(R, P, out):
    """invalid_fit_candidates.csv and invalid_fit_angle_x_distance.csv, for P(set) and for all
    LJs with a leading DSA-DSA pair."""
    rows, rows2 = [], []
    for nm, pop, x in [(nm, pop, x) for nm in SETS
                       for pop, x in (("P", le15(P[nm])), ("all LJs", le15(all_dd(R[nm]))))]:
        x = x.assign(az=x.cz.abs(), bounds=np.where((x.cr < R_MAX) & (x.cz.abs() < Z_MAX),
                                                    "inside", "outside"))
        for var, edges, label in (("cr", R_EDGES, "crossing r (cm)"),
                                  ("az", Z_EDGES, "crossing |z| (cm)"),
                                  ("bounds", None, "r < 740 cm and |z| < 960 cm"),
                                  ("emu_path", None, "TTMD path"),
                                  ("ang", ANG_EDGES, "angle at the crossing (rad)"),
                                  ("dr_ref", ANG_EDGES, "dR of the momenta at the reference points"),
                                  ("emu", DIST_EDGES, "distance (cm)")):
            rows.append(by_bins(x, var, edges, label).assign(set=nm, population=pop))
        for dl, dcol, edges in (("distance", "dbin", DIST_EDGES), ("coarse distance", "dc",
                                                                   DIST_COARSE)):
            g = x.groupby(["ab", dcol], observed=True).no_vtx.agg(n="size", k="sum").reset_index()
            g["angle_rad"] = g.ab.map(lambda i: f"[{ANG_EDGES[int(i)]:g},{ANG_EDGES[int(i) + 1]:g})")
            g["distance_cm"] = g[dcol].map(lambda i, e=edges: f"[{e[int(i)]:g},{e[int(i) + 1]:g})")
            rows2.append(g.assign(set=nm, population=pop, distance_bins=dl, f=g.k / g.n)
                         [["set", "population", "distance_bins", "angle_rad", "distance_cm", "n",
                           "k", "f"]])
    t = pd.concat(rows)[["set", "population", "variable", "bin", "n", "k", "f"]]
    t.to_csv(f"{out}/invalid_fit_candidates.csv", index=False)
    say("invalid-fit candidates (LJs with distance <= 15 cm)\n" + t[
        t.variable.isin(["angle at the crossing (rad)", "r < 740 cm and |z| < 960 cm"])]
        .round(3).to_string(index=False))
    t2 = pd.concat(rows2)
    t2.to_csv(f"{out}/invalid_fit_angle_x_distance.csv", index=False)
    say("invalid fits, P(set), angle x coarse distance\n" + t2[
        (t2.distance_bins == "coarse distance") & (t2.population == "P")]
        .round(3).to_string(index=False))


def deviance_table(P, out):
    """invalid_fit_deviance.csv."""
    rows = []
    allx = pd.concat([le15(P[nm]).assign(set=nm) for nm in SETS])
    for nm, x in [(nm, allx[allx.set == nm]) for nm in SETS] + [("pooled", allx)]:
        ok = np.isfinite(x.L) & x.ab.notna() & x.lb.notna() & x.dbin.notna() & x.eb.notna()
        y = x[ok]
        for lab, keys in (("one bin", ["one"]), ("distance", ["dbin"]), ("angle", ["ab"]),
                          ("angle x distance", ["ab", "dbin"]), ("L", ["lb"]),
                          ("L x distance", ["lb", "dbin"]), ("angle_eff", ["eb"]),
                          ("angle_eff x distance", ["eb", "dbin"])):
            dev, nb = deviance(y, keys)
            rows.append(dict(set=nm, description=lab, n=len(y), k=int(y.no_vtx.sum()),
                             excluded=int((~ok).sum()), bins=nb, deviance=dev,
                             aic=dev + 2 * nb))
    t = pd.DataFrame(rows)
    t["rank"] = t.groupby("set").deviance.rank(method="min").astype(int)
    t.to_csv(f"{out}/invalid_fit_deviance.csv", index=False)
    say("invalid-fit deviance\n" + t.round(1).to_string(index=False))


def _z1(ka, na, kb, nb):
    """z of group b against group a in one cell, as mh."""
    p0 = (ka + kb) / (na + nb) if na and nb else float("nan")
    if not 0 < p0 < 1:
        return float("nan")
    return (kb - nb * p0) / math.sqrt(na * nb * p0 * (1 - p0) / (na + nb - 1))


def boot_ratio(nr, kr, nt, kt):
    """Bootstrap of the standardisation (as wide_pair.py): NBOOT replicas resample the LJs of
    the reference and of the target separately (multinomial over bin x outcome), fixed seed."""
    cr, ct = np.stack([nr - kr, kr], 1).ravel(), np.stack([nt - kt, kt], 1).ravel()
    rng = np.random.default_rng(SEED)
    br = rng.multinomial(nr.sum(), cr / nr.sum(), size=NBOOT).reshape(NBOOT, -1, 2)
    bt = rng.multinomial(nt.sum(), ct / nt.sum(), size=NBOOT).reshape(NBOOT, -1, 2)
    n_r, m_t = br.sum(axis=2), bt.sum(axis=2)
    cov = n_r > 0
    f_r = np.where(cov, br[:, :, 1] / np.maximum(n_r, 1), 0.0)
    wcov = np.where(cov, m_t, 0).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        pred = (m_t * f_r).sum(axis=1) / wcov
        obs = np.where(cov, bt[:, :, 1], 0).sum(axis=1) / wcov
        ratio = obs / pred
    ratio = ratio[np.isfinite(ratio)]
    return dict(f_obs_kept_bs_err=np.nanstd(obs), f_pred_bs_err=np.nanstd(pred),
                ratio_bs_err=ratio.std(), ratio_bs_q16=np.quantile(ratio, 0.16),
                ratio_bs_q84=np.quantile(ratio, 0.84))


def angle_standardised_table(P, out):
    """invalid_fit_angle_standardised.csv: f_inv per angle (angle_eff) bin and set, and the
    sets other than 2Mu2E standardised to the per-bin f_inv of 2Mu2E."""
    allx = pd.concat([le15(P[nm]).assign(set=nm) for nm in SETS])
    ok = np.isfinite(allx.L) & allx.ab.notna() & allx.lb.notna() & allx.dbin.notna() & allx.eb.notna()
    y = allx[ok]
    rows = []
    for var, col, edges in (("angle", "ab", ANG_EDGES), ("angle_eff", "eb", AEFF_EDGES)):
        yb = y.assign(b=y[col].astype(int))
        bins = sorted(yb.b.unique())
        g = yb.groupby(["set", "b"]).no_vtx.agg(n="size", k="sum")

        def counts(nm, g=g, bins=bins):
            t = g.loc[nm].reindex(bins, fill_value=0)
            return t.n.to_numpy(), t.k.to_numpy()
        nr, kr = counts(SETS[0])
        for nm in SETS:
            n, k = counts(nm)
            for i, b in enumerate(bins):
                if n[i]:
                    f, err = frac(k[i], n[i])
                    rows.append(dict(table="per bin", variable=var, set=nm, bin=int(b),
                                     lo=edges[b], hi=edges[b + 1], n=int(n[i]), k=int(k[i]),
                                     f=f, err=err, z_vs_2Mu2E=(_z1(kr[i], nr[i], k[i], n[i])
                                                               if nm != SETS[0] else np.nan)))
        for nm in SETS[1:]:
            n, k = counts(nm)
            cov = nr > 0
            fr = np.where(cov, kr / np.maximum(nr, 1), 0.0)
            f, err = frac(k.sum(), n.sum())
            f_pred = (n * fr)[cov].sum() / n[cov].sum()
            f_obs = k[cov].sum() / n[cov].sum()
            rows.append(dict(table="standardised", variable=var, set=nm, n=int(n.sum()),
                             k=int(k.sum()), f=f, err=err, kept=n[cov].sum() / n.sum(),
                             f_obs_kept=f_obs, f_pred=f_pred, ratio=f_obs / f_pred,
                             **boot_ratio(nr, kr, n, k)))
    t = pd.DataFrame(rows)
    t.to_csv(f"{out}/invalid_fit_angle_standardised.csv", index=False)
    say("f_inv per angle bin and set, standardised to 2Mu2E\n" + t.round(4).to_string(index=False))


def sin_theta_table(P, out):
    """invalid_fit_sin_theta.csv: inside each region of P(2Mu2E), distance <= 15 cm, f_inv per
    quintile of sin(theta_pair), and the last quintile against the first in cells."""
    x = le15(P["2Mu2E"])
    rows = []
    for reg in ("barrel", "endcap"):
        y = x[x.region == reg].copy()
        y["sb"], edges = pd.qcut(y.sth, 5, labels=False, retbins=True)
        for b, g in y.groupby("sb"):
            rows.append(dict(region=reg, row="quintile", sin_bin=int(b), sin_lo=edges[int(b)],
                             sin_hi=edges[int(b) + 1], n=len(g), k=int(g.no_vtx.sum()),
                             f=g.no_vtx.mean()))
        for cl, keys in (("none", ["one"]), ("angle", ["ab"]),
                         ("angle x coarse distance", ["ab", "dc"])):
            rows.append(dict(region=reg, row="last against first quintile", cells_of=cl,
                             **mh(y, "sb", 0, 4, keys, 5)))
    t = pd.DataFrame(rows)
    t.to_csv(f"{out}/invalid_fit_sin_theta.csv", index=False)
    say("f_inv against sin(theta_pair) inside each region, P(2Mu2E)\n" + t.round(4).to_string(index=False))


def edge_table(R, out):
    """edge_1cm.csv: all LJs, 1 cm bins of the distance from 10 to 20 cm."""
    rows = []
    for nm in SETS:
        r = R[nm]
        dd = r[(r.m0_type == 1) & (r.m1_type == 1) & (r.emu_status == 1)]
        dd = dd[(dd.emu >= 10) & (dd.emu < 20)].assign(no_vtx=lambda x: (x.pv_n == 0).astype(int))
        rows.append(by_bins(dd, "emu", list(range(10, 21)), "distance (cm)").assign(set=nm))
        w = r[(r.m0_type == 1) & (r.m1_type == 1) & (r.pv_n > 0) & (r.emu > 15)]
        say(f"edge {nm}: LJs with a stored vertex and emulated distance > 15 cm {len(w)}; "
            f"paths {w.emu_path.value_counts().to_dict()}; stored dca max "
            f"{w.pv_dca.max() if len(w) else float('nan'):.3f}")
    t = pd.concat(rows)[["set", "bin", "n", "k", "f"]]
    t.to_csv(f"{out}/edge_1cm.csv", index=False)
    say("edge, 1 cm bins\n" + t.round(3).to_string(index=False))


def t7_tables(R, q, out):
    """T7: global minimum of the cut pairs, stored dca against it, minimiser stability."""
    rows, rows2, rows3 = [], [], []
    rng = np.random.default_rng(20261114)
    groups = []
    for nm in SETS:
        r = R[nm]
        dd = r[(r.m0_type == 1) & (r.m1_type == 1) & (r.emu_status >= 0)]
        groups.append((nm, dd[(dd.pv_n == 0) & (dd.emu > 15)]))
    for c in ("A narrow", "A wide"):
        x = q[(q.cls == c[0]) & (q.gw == c[2:]) & (q.no_vtx == 1) & (q.emu > 15)]
        groups.append((f"4Mu T6 class {c}", x))
    for nm, x in groups:
        d = x.gmin
        row = dict(group=nm, n=len(x), gmin_nan=int(d.isna().sum()))
        for t in (1, 5, 15):
            row[f"k_lt{t}"] = int((d < t).sum())
            row[f"f_lt{t}"] = (d < t).mean()
        row.update(median_gmin=d.median(), median_distance=x.emu.median())
        rows.append(row)
    for nm in SETS:
        r = R[nm]
        w = r[(r.m0_type == 1) & (r.m1_type == 1) & (r.emu_status >= 0) & (r.pv_n > 0)]
        w = w.assign(dbin=pd.cut(w.pv_dca, DCA_EDGES, right=False, labels=False))
        for b, y in [("all", w)] + [(f"[{DCA_EDGES[int(k)]:g},{DCA_EDGES[int(k) + 1]:g})", y)
                                    for k, y in w.groupby("dbin")]:
            ex = y.pv_dca - y.gmin
            rows2.append(dict(set=nm, stored_dca_cm=b, n=len(y), k_gt1=int((ex > 1).sum()),
                              f_gt1=(ex > 1).mean(), k_below=int((ex < -0.1).sum()),
                              gmin_nan=int(y.gmin.isna().sum())))
        dd = np.where(((r.m0_type == 1) & (r.m1_type == 1) & (r.emu_status >= 0)).to_numpy())[0]
        nv = dd[((r.pv_n.to_numpy()[dd] == 0) & (r.emu.to_numpy()[dd] > 15))]
        wv = dd[r.pv_n.to_numpy()[dd] > 0]
        rv = np.sort(np.concatenate([nv, rng.choice(wv, min(len(wv), 1000), replace=False)]))
        re = np.sort(np.concatenate([rng.choice(nv, min(len(nv), 30), replace=False),
                                     rng.choice(wv, min(len(wv), 30), replace=False)]))
        res = gmin_stability(r, rv, re)
        dv, de = res["variant"], res["exhaustive"]
        g0 = r.gmin.to_numpy()
        flip = sum(int(((g0[nv] < t) != (g0[nv] + dv[np.isin(rv, nv)] < t)).sum())
                   for t in (1, 5, 15))
        p1 = r.iloc[dd][r.emu_path.to_numpy()[dd] == 1]
        rows3.append(dict(set=nm, n_variant=len(rv), variant_max_abs_diff=np.nanmax(np.abs(dv)),
                          variant_n_gt_0p01=int((np.abs(dv) > 0.01).sum()),
                          variant_n_gt_0p1=int((np.abs(dv) > 0.1).sum()),
                          variant_class_changes_1_5_15=flip, n_exhaustive=len(re),
                          default_minus_exhaustive_max=-np.nanmin(de),
                          default_minus_exhaustive_min=-np.nanmax(de), n_path1=len(p1),
                          path1_gmin_above_ttmd=int((p1.gmin > p1.emu + 1e-3).sum())))
    for name, t in (("T7_global_minimum_cut_pairs", rows), ("T7_stored_dca_vs_global_minimum",
                                                            rows2),
                    ("T7_minimiser_stability", rows3)):
        t = pd.DataFrame(t)
        t.to_csv(f"{out}/{name}.csv", index=False)
        say(name + "\n" + t.round(4).to_string(index=False))


def t7_location_table(R, out):
    """T7_global_minimum_location.csv (docstring): the global-minimum midpoint and the crossing
    point of the cut pairs, and of the pairs with a stored vertex, against gen_a_lxy."""
    rows = []
    for nm in SETS:
        r = R[nm]
        dd = r[(r.m0_type == 1) & (r.m1_type == 1) & (r.emu_status >= 0)]
        cut = dd[(dd.pv_n == 0) & (dd.emu > 15)]
        for lab, x in (("cut pairs, all", cut),
                       ("cut pairs, global minimum < 15 cm", cut[cut.gmin < 15]),
                       ("cut pairs, global minimum >= 15 cm", cut[cut.gmin >= 15]),
                       ("control: pairs with a stored vertex", dd[dd.pv_n > 0])):
            gen = x.gen_matched.astype(bool) & (x.gen_a_lxy >= 0)
            u = x[gen & np.isfinite(x.gx)]
            gl = u.gen_a_lxy.to_numpy()
            rg = np.hypot(u.gx, u.gy).to_numpy()
            d = np.abs(rg - gl)
            row = dict(set=nm, group=lab, n=len(x), n_gen=int(gen.sum()), n_used=len(u))
            if len(u):
                q16, q50, q84 = np.quantile(d, [0.16, 0.5, 0.84])
                row.update(gen_lxy_median=np.median(gl), dr_median=q50, dr_q16=q16, dr_q84=q84,
                           f_dr_lt10=np.mean(d < 10), f_dr_lt20=np.mean(d < 20),
                           f_dr_gt100=np.mean(d > 100), f_r_gt300=np.mean(rg > 300),
                           f_absz_gt600=np.mean(np.abs(u.gz.to_numpy()) > 600))
            c = x[gen & np.isfinite(x.cr)]
            dc = np.abs(c.cr - c.gen_a_lxy).to_numpy()
            row.update(n_cross=len(c), cross_dr_median=np.median(dc) if len(c) else np.nan,
                       cross_f_dr_lt20=np.mean(dc < 20) if len(c) else np.nan)
            rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(f"{out}/T7_global_minimum_location.csv", index=False)
    say("T7_global_minimum_location\n" + t.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
