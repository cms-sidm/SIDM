"""Halo-sensitive fractions from the dsa_tail_diag histograms (flow-aware index groups).

Histograms (after selecting a channel), all Weight storage, read with flow=True:
  HD = mu_lj_dsaMuon_pt_nHits_absdxyBS  (pT [0,100,200,500], nHits [0,20,30,60], |dxyBS| [0,.1,.5,2,10,50,200])
  HT = mu_lj_dsaMuon_pt_nHits_nDTHits   (pT, nHits, nDTHits [0,1,6,12,24,60])
  HP = mu_lj_dsaMuon_pt_phi_abseta      (pT, phi Regular(16,-pi,pi), |eta| [0,.9,1.2,2.4])
"""
import math
import numpy as np

HD, HT, HP = ("mu_lj_dsaMuon_pt_nHits_absdxyBS", "mu_lj_dsaMuon_pt_nHits_nDTHits",
              "mu_lj_dsaMuon_pt_phi_abseta")
HALO_HISTS = (HD, HT, HP)
# flow-view index groups
PT = {"pT 0-100": [1], "pT 100-200": [2], "pT >=200": [3, 4], "pT all": [0, 1, 2, 3, 4]}
NH = {"nHits 0-20": [0, 1], "nHits 20-30": [2], "nHits >=30": [3, 4], "nHits all": [0, 1, 2, 3, 4]}
ETA = {"barrel |eta|<0.9": [1], "overlap 0.9-1.2": [2], "endcap 1.2-2.4": [3, 4],
       "|eta| all": [0, 1, 2, 3, 4]}
DXY_GE10 = [5, 6, 7]          # |dxyBS| >= 10 cm incl. overflow
DT_ZERO = [1]                 # nDTHits in [0, 1): CSC-only


def phi_masks(edges):
    """Vertical: bin centre within 0.5 rad of +-pi/2; horizontal: within 0.5 rad of 0 or +-pi
    (vertexChi2CosmicVetoStudy.ipynb sections 7 and 10.4, judged on bin centres)."""
    c = 0.5 * (edges[1:] + edges[:-1])
    vert = np.minimum(abs(c - math.pi / 2), abs(c + math.pi / 2)) < 0.5
    horz = np.minimum(abs(c), math.pi - abs(c)) < 0.5
    return vert, horz


def frac(a, va, n, vn):
    """f = a/n with a a subset of n; variance from the independent pieces a and n-a
    (equals f(1-f)/n for unit weights)."""
    if n <= 0:
        return float("nan"), float("nan")
    f = a / n
    var = ((1 - f) ** 2 * va + f ** 2 * (vn - va)) / n ** 2
    return f, math.sqrt(max(var, 0.0))


def ratio(num, vnum, den, vden):
    """num/den with independent variances."""
    if num <= 0 or den <= 0:
        return float("nan"), float("nan")
    r = num / den
    return r, r * math.sqrt(vnum / num ** 2 + vden / den ** 2)


def _neff(s, s2):
    return s * s / s2 if s2 > 0 else 0.0


def _take(a, idx, axis):
    return np.take(a, idx, axis=axis).sum(axis=axis)


def halo_rows(v, var, phi_edges):
    """v, var: dict hist name -> flow array for ONE channel (pt, y, z). Returns rows of dicts."""
    vert, horz = phi_masks(phi_edges)
    vmask = np.concatenate(([False], vert, [False]))
    hmask = np.concatenate(([False], horz, [False]))
    rows = []
    for pl, pi in PT.items():
        for nl, ni in NH.items():
            d, dv = _take(_take(v[HD], pi, 0), ni, 0), _take(_take(var[HD], pi, 0), ni, 0)
            t, tv = _take(_take(v[HT], pi, 0), ni, 0), _take(_take(var[HT], pi, 0), ni, 0)
            fd = frac(d[DXY_GE10].sum(), dv[DXY_GE10].sum(), d.sum(), dv.sum())
            ft = frac(t[DT_ZERO].sum(), tv[DT_ZERO].sum(), t.sum(), tv.sum())
            rows.append(dict(quantity="frac_absdxyBS_ge10", pt=pl, nhits=nl, eta="|eta| all",
                             value=fd[0], err=fd[1], n=d.sum(), neff=_neff(d.sum(), dv.sum())))
            rows.append(dict(quantity="frac_zeroDT_cscOnly", pt=pl, nhits=nl, eta="|eta| all",
                             value=ft[0], err=ft[1], n=t.sum(), neff=_neff(t.sum(), tv.sum())))
        p, pv = _take(v[HP], pi, 0), _take(var[HP], pi, 0)      # (phi, eta)
        for el, ei in ETA.items():
            q, qv = _take(p, ei, 1), _take(pv, ei, 1)
            r = ratio(q[hmask].sum(), qv[hmask].sum(), q[vmask].sum(), qv[vmask].sum())
            rows.append(dict(quantity="ratio_H_over_V_phi", pt=pl, nhits="nHits all", eta=el,
                             value=r[0], err=r[1], n=q[1:-1].sum(),
                             neff=_neff(q[1:-1].sum(), qv[1:-1].sum())))
    return rows
