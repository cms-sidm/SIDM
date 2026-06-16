"""Pure-NumPy port of the official MSTW 2008 PDF grid reader and interpolator.

This is a line-by-line translation of the MSTW collaboration's ``mstwpdf`` code
(G. Watt, mstwpdf.cc / mstwpdf.f / mstwpdf.m) so that PDF values reproduce the
original Fortran/Mathematica implementation to machine precision.  Only the
central set (eigenvector ih=0) is needed here, so the Hessian-error machinery is
dropped; everything else -- the bicubic spline over a (log10 x, log10 Q^2) grid,
the heavy-quark threshold handling, and the low-x / high-Q^2 extrapolation -- is
kept faithfully.

Flavour code ``f`` follows the PDG convention except the gluon is 0 (not 21):
    f = -5..-1, 0, 1..5  ->  bbar..dbar, g, d..b   (and 6/-6 = top = 0).
"""

import re
import numpy as np

# --- fixed grid layout (identical to mstwpdf) -------------------------------
XMIN, XMAX = 1.0e-6, 1.0
QSQMIN, QSQMAX = 1.0, 1.0e9
EPS = 1.0e-6
NP, NX, NQ = 12, 64, 48

XX = np.array([
    1e-6, 2e-6, 4e-6, 6e-6, 8e-6, 1e-5, 2e-5, 4e-5, 6e-5, 8e-5,
    1e-4, 2e-4, 4e-4, 6e-4, 8e-4, 1e-3, 2e-3, 4e-3, 6e-3, 8e-3,
    1e-2, 1.4e-2, 2e-2, 3e-2, 4e-2, 6e-2, 8e-2, 0.1, 0.125, 0.15,
    0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4,
    0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65,
    0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
    0.925, 0.95, 0.975, 1.0])

# QQ holds Q^2 nodes (GeV^2); zeros are placeholders filled at the charm/bottom
# thresholds when the grid is read.
QQ_ORIG = np.array([
    1e0, 1.25e0, 1.5e0, 0., 0., 2.5e0, 3.2e0, 4e0, 5e0, 6.4e0,
    8e0, 1e1, 1.2e1, 0., 0., 2.6e1, 4e1, 6.4e1, 1e2, 1.6e2,
    2.4e2, 4e2, 6.4e2, 1e3, 1.8e3, 3.2e3, 5.6e3, 1e4, 1.8e4, 3.2e4,
    5.6e4, 1e5, 1.8e5, 3.2e5, 5.6e5, 1e6, 1.8e6, 3.2e6, 5.6e6, 1e7,
    1.8e7, 3.2e7, 5.6e7, 1e8, 1.8e8, 3.2e8, 5.6e8, 1e9])

# Numerical-Recipes bcucof weight matrix (16x16).
_IWT = np.array([
    [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
    [-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0],
    [2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0],
    [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
    [0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1],
    [0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1],
    [-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0],
    [9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2],
    [-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2],
    [2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0],
    [-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1],
    [4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1],
], dtype=float)


# --- three-point derivative estimators (mstwpdf polderiv1/2/3) ---------------
def _polderiv1(x1, x2, x3, y1, y2, y3):
    return (x3*x3*(y1-y2) + 2.*x1*(x3*(-y1+y2)+x2*(y1-y3))
            + x2*x2*(-y1+y3) + x1*x1*(-y2+y3)) / ((x1-x2)*(x1-x3)*(x2-x3))


def _polderiv2(x1, x2, x3, y1, y2, y3):
    return (x3*x3*(y1-y2) - 2.*x2*(x3*(y1-y2)+x1*(y2-y3))
            + x2*x2*(y1-y3) + x1*x1*(y2-y3)) / ((x1-x2)*(x1-x3)*(x2-x3))


def _polderiv3(x1, x2, x3, y1, y2, y3):
    return (x3*x3*(-y1+y2) + 2.*x2*x3*(y1-y3) + x1*x1*(y2-y3)
            + x2*x2*(-y1+y3) + 2.*x1*x3*(-y2+y3)) / ((x1-x2)*(x1-x3)*(x2-x3))


class MSTWPDF:
    """Reads one MSTW2008 grid file and interpolates xf(x, Q)."""

    def __init__(self, grid_path):
        self._read_grid(grid_path)
        self._init_all()

    # ---- grid file parsing (mstwpdf ReadPDFGrid) ---------------------------
    def _read_grid(self, path):
        text = open(path).read()
        def _hdr(name):
            m = re.search(name + r"\s*=\s*([-+0-9.eE]+)", text)
            return float(m.group(1))
        self.mCharm = _hdr("mCharm")
        self.mBottom = _hdr("mBottom")
        self.alphaSorder = int(_hdr(r"alphaSorder,alphaSnfmax"))
        self.alphaSMZ = _hdr(r"alphaS\(MZ\)")
        mex = re.search(r"nExtraFlavours\s*=\s*([0-9]+)", text)
        self.nExtraFlavours = int(mex.group(1))

        lines = text.splitlines()
        hi = next(i for i, l in enumerate(lines) if "xg" in l and "xd" in l)
        nums = []
        for l in lines[hi + 1:]:
            p = l.split()
            if p:
                nums.extend(p)
        nums = np.array(nums, dtype=float)

        # Build the Q^2 node array with charm/bottom thresholds (mstwpdf logic).
        qq = QQ_ORIG.copy()
        mc2, mb2 = self.mCharm**2, self.mBottom**2
        # charm: indices are 1-based in mstwpdf; convert to 0-based here.
        if mc2 <= qq[0] or mc2 + EPS >= qq[7]:
            raise ValueError("invalid mCharm")
        elif mc2 < qq[1]:
            self.nqc0 = 2; qq[3] = qq[1]; qq[4] = qq[2]
        elif mc2 < qq[2]:
            self.nqc0 = 3; qq[4] = qq[2]
        elif mc2 < qq[5]:
            self.nqc0 = 4
        elif mc2 < qq[6]:
            self.nqc0 = 5; qq[3] = qq[5]
        else:
            self.nqc0 = 6; qq[3] = qq[5]; qq[4] = qq[6]
        if mb2 <= qq[11] or mb2 + EPS >= qq[16]:
            raise ValueError("invalid mBottom")
        elif mb2 < qq[12]:
            self.nqb0 = 13; qq[14] = qq[12]
        elif mb2 < qq[15]:
            self.nqb0 = 14
        else:
            self.nqb0 = 15; qq[13] = qq[15]
        # nqc0/nqb0 are 1-based node numbers; place thresholds (0-based slots).
        qq[self.nqc0 - 1] = mc2
        qq[self.nqc0] = mc2 + EPS
        qq[self.nqb0 - 1] = mb2
        qq[self.nqb0] = mb2 + EPS
        self.qq = qq
        self.xxlog = np.log10(XX)
        self.qqlog = np.log10(qq)

        # Fill ff[ip, ix, iq] with the read order ix (outer), iq (mid), ip (inner).
        ff = np.zeros((NP, NX, NQ))
        k = 0
        for ix in range(NX):
            for iq in range(NQ):
                for ip in range(NP):
                    forced = (ix == NX - 1
                              or (self.alphaSorder < 2 and ip in (9, 10))
                              or (self.nExtraFlavours < 1 and ip == 11))
                    if forced:
                        ff[ip, ix, iq] = 0.0
                    else:
                        ff[ip, ix, iq] = nums[k]; k += 1
        if k != len(nums):
            raise ValueError(f"grid read mismatch: used {k} of {len(nums)} numbers")
        self.ff = ff

    # ---- bicubic coefficient setup (mstwpdf InitialisePDF) -----------------
    def _init_all(self):
        # cc[ip, n, m, k, j] : bicubic coeffs per cell (n,m), 0-based.
        self.cc = np.zeros((NP, NX - 1, NQ - 1, 4, 4))
        for ip in range(NP):
            self._init_pdf(ip)

    def _dq_index(self, m):
        """Return derivative scheme for q-node m (0-based): 'f','b','c'."""
        nqc0, nqb0 = self.nqc0, self.nqb0  # 1-based
        m1 = m + 1  # 1-based
        if m1 == 1 or m1 == nqc0 + 1 or m1 == nqb0 + 1:
            return 'f'
        if m1 == NQ or m1 == nqc0 or m1 == nqb0:
            return 'b'
        return 'c'

    def _q_deriv(self, arr):
        """d(arr)/d(qqlog) over a (NX,NQ) array, with threshold schemes."""
        ql = self.qqlog
        out = np.zeros_like(arr)
        nqc0 = self.nqc0
        for m in range(NQ):
            sch = self._dq_index(m)
            if nqc0 == 2 and m == 0:
                out[:, m] = (arr[:, m + 1] - arr[:, m]) / (ql[m + 1] - ql[m])
            elif nqc0 == 2 and m == 1:
                out[:, m] = (arr[:, m] - arr[:, m - 1]) / (ql[m] - ql[m - 1])
            elif sch == 'f':
                out[:, m] = _polderiv1(ql[m], ql[m+1], ql[m+2],
                                       arr[:, m], arr[:, m+1], arr[:, m+2])
            elif sch == 'b':
                out[:, m] = _polderiv3(ql[m-2], ql[m-1], ql[m],
                                       arr[:, m-2], arr[:, m-1], arr[:, m])
            else:
                out[:, m] = _polderiv2(ql[m-1], ql[m], ql[m+1],
                                       arr[:, m-1], arr[:, m], arr[:, m+1])
        return out

    def _x_deriv(self, arr):
        """d(arr)/d(xxlog) over a (NX,NQ) array (3-point, endpoints one-sided)."""
        xl = self.xxlog
        out = np.zeros_like(arr)
        out[0, :] = _polderiv1(xl[0], xl[1], xl[2], arr[0, :], arr[1, :], arr[2, :])
        out[NX-1, :] = _polderiv3(xl[NX-3], xl[NX-2], xl[NX-1],
                                  arr[NX-3, :], arr[NX-2, :], arr[NX-1, :])
        for n in range(1, NX - 1):
            out[n, :] = _polderiv2(xl[n-1], xl[n], xl[n+1],
                                   arr[n-1, :], arr[n, :], arr[n+1, :])
        return out

    def _init_pdf(self, ip):
        f = self.ff[ip]
        ff1 = self._x_deriv(f)
        ff2 = self._q_deriv(f)
        ff12 = self._x_deriv(ff2)
        ff21 = self._q_deriv(ff1)
        ff12 = 0.5 * (ff12 + ff21)

        xl, ql = self.xxlog, self.qqlog
        for n in range(NX - 1):
            d1 = xl[n+1] - xl[n]
            for m in range(NQ - 1):
                d2 = ql[m+1] - ql[m]
                d1d2 = d1 * d2
                yy0 = np.array([f[n, m], f[n+1, m], f[n+1, m+1], f[n, m+1]])
                yy1 = np.array([ff1[n, m], ff1[n+1, m], ff1[n+1, m+1], ff1[n, m+1]])
                yy2 = np.array([ff2[n, m], ff2[n+1, m], ff2[n+1, m+1], ff2[n, m+1]])
                yy12 = np.array([ff12[n, m], ff12[n+1, m], ff12[n+1, m+1], ff12[n, m+1]])
                z = np.empty(16)
                z[0:4] = yy0
                z[4:8] = yy1 * d1
                z[8:12] = yy2 * d2
                z[12:16] = yy12 * d1d2
                cl = _IWT.dot(z)
                self.cc[ip, n, m] = cl.reshape(4, 4)

    # ---- node location (mstwpdf locx / locq) -------------------------------
    @staticmethod
    def _loc(grid, val):
        n = len(grid)
        if val == grid[0]:
            return 0
        if val == grid[-1]:
            return n - 2
        ju, jl = n, -1
        while ju - jl > 1:
            jm = (ju + jl) // 2
            if val >= grid[jm]:
                jl = jm
            else:
                ju = jm
        return jl

    def _interpolate(self, ip, xlog, qlog):
        n = self._loc(self.xxlog, xlog)
        m = self._loc(self.qqlog, qlog)
        t = (xlog - self.xxlog[n]) / (self.xxlog[n+1] - self.xxlog[n])
        u = (qlog - self.qqlog[m]) / (self.qqlog[m+1] - self.qqlog[m])
        c = self.cc[ip, n, m]
        z = 0.0
        for l in range(3, -1, -1):
            z = t * z + ((c[l, 3]*u + c[l, 2])*u + c[l, 1])*u + c[l, 0]
        return z

    def _extrapolate(self, ip, xlog, qlog):
        xl, ql = self.xxlog, self.qqlog
        n = self._loc(xl, xlog)
        m = self._loc(ql, qlog)

        def loglin(f0, f1, a, a0, a1):
            if f0 >= 1e-3 and f1 >= 1e-3:
                return np.exp(np.log(f0) + (np.log(f1) - np.log(f0)) / (a1 - a0) * (a - a0))
            return f0 + (f1 - f0) / (a1 - a0) * (a - a0)

        if n == -1 and 0 <= m < NQ - 1:   # low-x edge (mstwpdf n==0, 1-based)
            f0 = self._interpolate(ip, xl[0], qlog)
            f1 = self._interpolate(ip, xl[1], qlog)
            return loglin(f0, f1, xlog, xl[0], xl[1])
        if n >= 0 and m == NQ - 1:        # high-Q^2 edge
            f0 = self._interpolate(ip, xlog, ql[NQ-1])
            f1 = self._interpolate(ip, xlog, ql[NQ-2])
            if f0 >= 1e-3 and f1 >= 1e-3:
                return np.exp(np.log(f0) + (np.log(f0) - np.log(f1)) /
                              (ql[NQ-1] - ql[NQ-2]) * (qlog - ql[NQ-1]))
            return f0 + (f0 - f1) / (ql[NQ-1] - ql[NQ-2]) * (qlog - ql[NQ-1])
        if n == -1 and m == NQ - 1:       # both edges
            def hi_q(xn):
                f0 = self._interpolate(ip, xn, ql[NQ-1])
                f1 = self._interpolate(ip, xn, ql[NQ-2])
                if f0 >= 1e-3 and f1 >= 1e-3:
                    return np.exp(np.log(f0) + (np.log(f0) - np.log(f1)) /
                                  (ql[NQ-1] - ql[NQ-2]) * (qlog - ql[NQ-1]))
                return f0 + (f0 - f1) / (ql[NQ-1] - ql[NQ-2]) * (qlog - ql[NQ-1])
            z0, z1 = hi_q(xl[0]), hi_q(xl[1])
            return loglin(z0, z1, xlog, xl[0], xl[1])
        raise RuntimeError("ExtrapolatePDF: unexpected region")

    def _interp_vec(self, ip, xlogs, qlog):
        """Vectorised bicubic interpolation over an array of xlog at fixed qlog.

        Interpolation only (assumes in-grid); used for the parton-luminosity
        integrals where x and Q always lie inside the grid.
        """
        n = np.searchsorted(self.xxlog, xlogs, side='right') - 1
        n = np.clip(n, 0, NX - 2)
        m = self._loc(self.qqlog, qlog)
        t = (xlogs - self.xxlog[n]) / (self.xxlog[n + 1] - self.xxlog[n])
        u = (qlog - self.qqlog[m]) / (self.qqlog[m + 1] - self.qqlog[m])
        c = self.cc[ip, n, m]            # (len, 4, 4)
        z = np.zeros_like(xlogs)
        for l in range(3, -1, -1):
            z = t * z + ((c[:, l, 3] * u + c[:, l, 2]) * u + c[:, l, 1]) * u + c[:, l, 0]
        return z

    def _ip0(self, f):
        if f == 0:
            return 0
        if 1 <= f <= 5:
            return f
        if -5 <= f <= -1:
            return -f
        if 7 <= f <= 11:
            return f - 1
        if f == 13:
            return 11
        return None

    def xf_grid(self, xs, q, f):
        """Vectorised x*f(x,Q): xs an array of in-grid momentum fractions.

        Matches xf() to machine precision in the interpolation region (validated).
        """
        xs = np.asarray(xs, dtype=float)
        qsq = q * q
        nqc0, nqb0 = self.nqc0, self.nqb0
        if self.qq[nqc0 - 1] < qsq < self.qq[nqc0]:
            qsq = self.qq[nqc0]
        if self.qq[nqb0 - 1] < qsq < self.qq[nqb0]:
            qsq = self.qq[nqb0]
        if abs(f) == 6 or f == 12:
            return np.zeros_like(xs)
        ip0 = self._ip0(f)
        xlogs = np.log10(xs)
        qlog = np.log10(qsq)
        res = self._interp_vec(ip0, xlogs, qlog)
        if -5 <= f <= -1:
            res = res - self._interp_vec(ip0 + 5, xlogs, qlog)
        return res

    # ---- public PDF call (mstwpdf xf) --------------------------------------
    def xf(self, x, q, f):
        """x*f(x,Q) for flavour f at momentum fraction x and scale q (GeV)."""
        x = float(x); q = float(q); qsq = q * q
        nqc0, nqb0 = self.nqc0, self.nqb0
        if self.qq[nqc0-1] < qsq < self.qq[nqc0]:
            qsq = self.qq[nqc0]
        if self.qq[nqb0-1] < qsq < self.qq[nqb0]:
            qsq = self.qq[nqb0]
        xlog = np.log10(x)
        qlog = np.log10(qsq)

        if f == 0:
            ip = 1
        elif 1 <= f <= 5:
            ip = f + 1
        elif -5 <= f <= -1:
            ip = -f + 1
        elif 7 <= f <= 11:
            ip = f
        elif f == 13:
            ip = 12
        elif abs(f) != 6 and f != 12:
            raise ValueError(f"invalid f={f}")
        else:
            ip = None
        ip0 = (ip - 1) if ip is not None else None  # 0-based

        if x <= 0. or x > XMAX or q <= 0.:
            raise ValueError(f"invalid (x,q)=({x},{q})")
        if abs(f) == 6 or f == 12:
            return 0.0

        if qsq < QSQMIN:
            qmin_log = np.log10(QSQMIN)
            qmin_log2 = np.log10(1.01 * QSQMIN)
            interp = self._extrapolate if x < XMIN else self._interpolate
            res = interp(ip0, xlog, qmin_log)
            res1 = interp(ip0, xlog, qmin_log2)
            if -5 <= f <= -1:
                res -= interp(ip0 + 5, xlog, qmin_log)
                res1 -= interp(ip0 + 5, xlog, qmin_log2)
            anom = max(-2.5, (res1 - res) / res / 0.01) if abs(res) >= 1e-5 else 1.0
            res *= (qsq / QSQMIN) ** (anom * qsq / QSQMIN + 1. - qsq / QSQMIN)
            return res
        if x < XMIN or qsq > QSQMAX:
            res = self._extrapolate(ip0, xlog, qlog)
            if -5 <= f <= -1:
                res -= self._extrapolate(ip0 + 5, xlog, qlog)
            return res
        res = self._interpolate(ip0, xlog, qlog)
        if -5 <= f <= -1:
            res -= self._interpolate(ip0 + 5, xlog, qlog)
        return res
