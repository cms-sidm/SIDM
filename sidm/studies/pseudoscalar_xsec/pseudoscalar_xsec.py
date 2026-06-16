"""Pseudoscalar dark-matter bound-state production cross section at the LHC.

Direct port of the Mathematica notebook ``PscalarXsec.nb``.  Computes the pp ->
(darkonium) production cross section for a Coulombic dark-matter bound state that
mixes with a CP-odd pseudoscalar mediator ``A`` in a dark-matter simplified
model, via gluon fusion (top/bottom loops) and qqbar fusion, convolved with
MSTW2008 NNLO PDFs and multiplied by QCD K-factors.

All physics constants, the A-B mass mixing, the pseudoscalar loop form factor,
and the parton-luminosity integrals reproduce the original notebook.

Normalization and uncertainties
-------------------------------
Absolute cross sections are returned in fb, with two caveats on their overall
scale. (1) Following the original notebook, the GeV^-2 -> fb conversion uses
hbar*c ~ 2e-14 cm.GeV (``GEV2_TO_FB`` below), which rounds the exact (hbar*c)^2
up by ~2.7%; every absolute rate is high by that factor (divide by 1.027 for the
exact normalization). (2) This is a single-PDF-member (MSTW2008 NNLO central),
fixed-K-factor calculation: it carries no PDF or QCD-scale uncertainty -- the
missing gluon-fusion scale uncertainty alone is O(10-20%), far larger than the
2.7% offset. The numbers are a central-value estimate, not a precision prediction.
"""

from functools import lru_cache

import numpy as np
from mstwpdf import MSTWPDF


@lru_cache(maxsize=None)
def _gauss_legendre(n):
    """Cached Gauss-Legendre nodes/weights (the n=400 rule is reused heavily)."""
    return np.polynomial.legendre.leggauss(n)

# --- fixed inputs (notebook "constants" cell) -------------------------------
SQRT_S = 13000.0     # pp centre-of-mass energy, GeV
ALPHA_S = 0.12       # strong coupling used in the ggA effective vertex
VEV = 174.0          # v = 246/sqrt(2), GeV  (so y_q = m_q / VEV)
MT = 174.0           # top pole mass, GeV (loops + A width)
MB = 4.7             # bottom pole mass, GeV (loops + A width)
YB = MB / VEV
YT = MT / VEV
# running quark masses used for the tree-level qqbar->A couplings:
MC, MS, MU, MD = 1.28, 0.095, 0.0023, 0.0048

ZETA3 = 1.2020569031595942854   # Apery's constant, zeta(3)
# GeV^-2 -> fb conversion as used in the original notebook: (hbar*c ~ 2e-14 cm.GeV).
# This rounds the exact (hbar*c)^2 = 3.894e11 up by ~2.7%; absolute rates carry that
# offset (kept for faithfulness to the original -- see the module docstring).
GEV2_TO_FB = (2.0 / 1e14) ** 2 / 1e-39   # = 4.0e11 fb per GeV^-2 (exact: 3.894e11)

S = SQRT_S ** 2


def floop(tau):
    """Pseudoscalar one-loop function f(tau) (complex above threshold)."""
    if tau <= 1.0:
        return np.arcsin(np.sqrt(tau)) ** 2
    r = np.sqrt(1.0 - 1.0 / tau)
    return -0.25 * (np.log((1.0 + r) / (1.0 - r)) - 1j * np.pi) ** 2


def A12A(tau):
    """Pseudoscalar fermion-loop form factor A_{1/2}^A(tau) = (2/tau) f(tau)."""
    return (2.0 / tau) * floop(tau)


def _gammaA(mA, tanb, yq):
    bb = (3.0 / (8 * np.pi)) * ((YB * yq) / np.sqrt(2)) ** 2 * tanb ** 2 * mA \
        * np.sqrt(1.0 - 4.0 * MB ** 2 / mA ** 2)
    tt = 0.0
    if mA > 2 * MT:
        tt = (3.0 / (8 * np.pi)) * ((YT * yq) / np.sqrt(2)) ** 2 * (1.0 / tanb ** 2) * mA \
            * np.sqrt(1.0 - 4.0 * MT ** 2 / mA ** 2)
    return abs(bb + tt)


def _eigen(mDM, mA, tanb, alphaD, yq):
    GammaA = _gammaA(mA, tanb, yq)
    GammaBS = (alphaD ** 5 / 4.0) * (2 * mDM)
    f2 = alphaD ** 3 * mDM ** 2 / (8 * np.pi)
    M = np.array([[mA ** 2 - 1j * mA * GammaA, -f2],
                  [-f2, (2 * mDM) ** 2 - 1j * (2 * mDM) * GammaBS]], dtype=complex)
    ev = np.linalg.eigvals(M)
    order = np.argsort(-np.abs(ev))     # descending |eigenvalue| (Mathematica order)
    ev = ev[order]
    mM2, mB2 = ev[0].real, ev[1].real
    GammaMmM, GammaBmB = ev[0].imag, ev[1].imag
    return mM2, mB2, GammaMmM, GammaBmB


def _lumi_integral(pdf, mB2, fa, fb, n=400):
    """Parton luminosity Integral_{tau0}^{1} xf_a(x) xf_b(tau0/x)/(x*tau0) dx,
    where xf returns x*f(x).  Evaluated in y=ln(x) (dx = x dy) with an
    n-point Gauss-Legendre rule; matches the notebook NIntegrate to ~1e-7.
    """
    Q = np.sqrt(mB2)
    tau0 = mB2 / S
    ylo, yhi = np.log(tau0), 0.0
    nodes, weights = _gauss_legendre(n)
    ys = 0.5 * (yhi - ylo) * nodes + 0.5 * (yhi + ylo)
    jac = 0.5 * (yhi - ylo)
    xs = np.exp(ys)
    x2s = tau0 / xs
    fa_vals = pdf.xf_grid(xs, Q, fa)
    fb_vals = pdf.xf_grid(x2s, Q, fb)
    integrand = (fa_vals * fb_vals) / (xs * tau0) * xs   # * x from dx = x dy
    return jac * np.sum(weights * integrand)


def cross_gg(pdf, mDM, mA, tanb, alphaD, yq, ychi, n=400):
    """Gluon-fusion bound-state production cross section [fb], pre-K-factor.

    Absolute scale carries the +2.7% normalization offset and has no PDF/scale
    uncertainty; see the module docstring.
    """
    mM2, mB2, GMmM, GBmB = _eigen(mDM, mA, tanb, alphaD, yq)
    lam_inv2 = (ALPHA_S ** 2 * yq ** 2 / (16 ** 2 * (np.sqrt(2) * VEV) ** 2 * np.pi ** 2)) * \
        abs(A12A(mB2 / (4 * MT ** 2)) * (1.0 / tanb)
            + A12A(mB2 / (4 * MB ** 2)) * tanb) ** 2
    bw = (mB2 - mM2) ** 2 + (GMmM - GBmB) ** 2
    lumi = _lumi_integral(pdf, mB2, 0, 0, n=n)
    sigma = (ZETA3 * np.pi / 256.0) * (alphaD ** 3 * ychi ** 2 * mB2 ** 3 / bw) \
        * lam_inv2 * (1.0 / S) * lumi * GEV2_TO_FB
    return sigma


def _ykwq(tanb):
    """Tree-level qqbar->A couplings g_q for d,u,s,c,b,t (notebook ykwq list)."""
    return np.array([
        (48 * tanb) / (1e4 * VEV),       # d  (down-type ~ tanb)
        23 / (1e4 * tanb * VEV),         # u  (up-type   ~ 1/tanb)
        (95 * tanb) / (1e3 * VEV),       # s
        128 / (100 * tanb * VEV),        # c
        (418 * tanb) / (100 * VEV),      # b
        174 / (tanb * VEV),              # t  (PDF=0, no contribution)
    ])


def cross_qq(pdf, mDM, mA, tanb, alphaD, yq, ychi, n=400):
    """qqbar-fusion bound-state production cross section [fb], pre-K-factor.

    Absolute scale carries the +2.7% normalization offset and has no PDF/scale
    uncertainty; see the module docstring.
    """
    mM2, mB2, GMmM, GBmB = _eigen(mDM, mA, tanb, alphaD, yq)
    bw = (mB2 - mM2) ** 2 + (GMmM - GBmB) ** 2
    yk = _ykwq(tanb)
    ssum = 0.0
    for i in range(1, 7):                 # flavours d,u,s,c,b,t -> f = i
        lumi = _lumi_integral(pdf, mB2, -i, i, n=n)
        ssum += yk[i - 1] ** 2 * lumi
    sigma = (ZETA3 / 48.0) * alphaD ** 3 * ychi ** 2 * (yq ** 2 / 2.0) \
        * (mB2 ** 2 / bw) * (2.0 / S) * ssum * GEV2_TO_FB
    return sigma


# --- K-factor tables (notebook ggAKfactorTab / qqHLOTab / qqHNNLOTab) --------
_ggA_K = np.array([[99,3.241497],[103.,3.241497],[107.5,3.207483],[116.5,3.159864],[124.,3.119048],[131.5,3.071429],[139.,3.017007],[145.,2.97619],[154.459459,2.943433],[169.559033,2.886787],[186.525097,2.841108],[197.5,2.79932],[213.438669,2.759988],[230.5,2.710884],[251.556757,2.672495],[270.745174,2.64551],[291.513514,2.629987],[311.135135,2.603879],[330.756757,2.581449],[350.378378,2.561224],[370.,2.542471],[389.621622,2.520776],[409.243243,2.503493],[428.864865,2.487314],[448.486486,2.472973],[468.108108,2.459735],[487.72973,2.452749],[507.351351,2.443188],[526.972973,2.431789],[546.594595,2.422596],[566.216216,2.411197],[585.837838,2.400165],[605.459459,2.395017],[625.081081,2.391708],[644.702703,2.386192],[664.324324,2.379206],[683.945946,2.367807],[703.567568,2.361555],[723.189189,2.354569],[742.810811,2.351995],[762.432432,2.346479],[782.054054,2.340963],[801.675676,2.336919],[821.297297,2.334345],[841.,2.323129],[860.540541,2.323313],[879.864865,2.312098],[899.783784,2.314856],[919.405405,2.313017],[939.027027,2.306766],[958.648649,2.306031],[978.27027,2.302721],[993.432432,2.297941],[1000.432432,2.297941]])
_qqH_LO = np.array([[99,0.214313],[102.299138,0.214313],[110.271782,0.192763],[117.611676,0.171119],[125.605412,0.149659],[132.291465,0.135182],[139.631359,0.119706],[146.971254,0.106265],[153.931498,0.095985],[161.271393,0.085631],[168.421462,0.076584],[176.225373,0.067648],[182.911426,0.061408],[189.871671,0.055468],[197.211565,0.049731],[204.17181,0.045144],[211.511704,0.040676],[218.471949,0.0372],[225.432193,0.033769],[233.278287,0.030201],[241.124381,0.02762],[248.084626,0.025072],[255.044871,0.023101],[262.005115,0.021232],[268.96536,0.019514],[275.925605,0.018025],[282.885849,0.016566],[290.383931,0.015122],[296.975072,0.014228],[400.,0.014228],[1001,0.014228]])
_qqH_NNLO = np.array([[99,0.917276],[100.183486,0.917276],[103.691187,0.831307],[111.227937,0.691283],[118.188182,0.587503],[125.078121,0.503276],[132.389893,0.427868],[138.919516,0.372513],[146.591604,0.317204],[153.551849,0.276425],[160.512093,0.241433],[166.909894,0.212995],[174.362277,0.186079],[181.471921,0.163745],[188.282767,0.145189],[195.250042,0.128957],[202.03628,0.117135],[208.440367,0.105925],[216.013265,0.09176],[223.224601,0.082049],[230.536373,0.073279],[237.707534,0.06599],[244.588685,0.059755],[251.390742,0.054109],[258.350987,0.049994],[262.752294,0.047315],[266.055046,0.045316],[271.192661,0.04217],[275.59633,0.039811],[279.469002,0.039057],[286.429247,0.035728],[293.389491,0.032817],[298.451487,0.030837],[400.,0.030837],[1001,0.030837]])


def Kgg(mass):
    """gluon-fusion QCD K-factor at scale ``mass`` (=2 m_DM), linear interp."""
    return float(np.interp(mass, _ggA_K[:, 0], _ggA_K[:, 1]))


def Kqq(mass):
    """qqbar QCD K-factor = sigma_NNLO/sigma_LO at scale ``mass`` (=2 m_DM)."""
    lo = np.interp(mass, _qqH_LO[:, 0], _qqH_LO[:, 1])
    nnlo = np.interp(mass, _qqH_NNLO[:, 0], _qqH_NNLO[:, 1])
    return float(nnlo / lo)


def total_xsec(pdf, mDM, mA, tanb, alphaD, yq, ychi, n=400):
    """Total pp -> bound-state cross section [fb], K-factors included.

    Absolute scale carries the +2.7% normalization offset and has no PDF/scale
    uncertainty; see the module docstring.
    """
    sgg = cross_gg(pdf, mDM, mA, tanb, alphaD, yq, ychi, n=n)
    sqq = cross_qq(pdf, mDM, mA, tanb, alphaD, yq, ychi, n=n)
    return sgg * Kgg(2 * mDM) + sqq * Kqq(2 * mDM)


# --- convenience wrappers for the notebook ----------------------------------
def mdm_from_mb(mB, alphaD):
    """Dark-matter mass from the Coulombic ground-state (bound-state) mass mB.

    The 1S binding energy lowers the bound-state mass to mB = mDM (2 - alphaD^2/4),
    so mDM = mB / (2 - alphaD^2/4).
    """
    return mB / (2.0 - alphaD ** 2 / 4.0)


def gamma_A(mA, tanb, yq):
    """Total width of the pseudoscalar mediator A (A->bbar + A->ttbar) [GeV]."""
    return _gammaA(mA, tanb, yq)


def mixing(mDM, mA, tanb, alphaD, yq):
    """Diagonalise the A-B mass-squared matrix.

    Returns dict with the two physical masses (mass_M = mediator-like,
    mass_B = bound-state-like) [GeV] and their widths [GeV].
    """
    mM2, mB2, GMmM, GBmB = _eigen(mDM, mA, tanb, alphaD, yq)
    return {
        "mM2": mM2, "mB2": mB2,
        "mass_M": np.sqrt(mM2), "mass_B": np.sqrt(mB2),
        "Gamma_M": abs(GMmM) / np.sqrt(mM2),
        "Gamma_B": abs(GBmB) / np.sqrt(mB2),
    }


def breakdown(pdf, mDM, mA, tanb, alphaD, yq, ychi, n=400):
    """Full intermediate + final breakdown for one parameter point.

    Absolute scale carries the +2.7% normalization offset and has no PDF/scale
    uncertainty; see the module docstring.
    """
    mix = mixing(mDM, mA, tanb, alphaD, yq)
    sgg = cross_gg(pdf, mDM, mA, tanb, alphaD, yq, ychi, n=n)
    sqq = cross_qq(pdf, mDM, mA, tanb, alphaD, yq, ychi, n=n)
    kgg, kqq = Kgg(2 * mDM), Kqq(2 * mDM)
    return {
        "mDM": mDM, "mA": mA, "tanb": tanb, "alphaD": alphaD,
        "yq": yq, "ychi": ychi,
        "Gamma_A": _gammaA(mA, tanb, yq),
        "mass_boundstate": mix["mass_B"], "mass_mediatorlike": mix["mass_M"],
        "Gamma_boundstate_like": mix["Gamma_B"],
        "sigma_gg_LO": sgg, "sigma_qq_LO": sqq, "Kgg": kgg, "Kqq": kqq,
        "sigma_gg": sgg * kgg, "sigma_qq": sqq * kqq,
        "sigma_total": sgg * kgg + sqq * kqq,
    }


def scan_mB(pdf, mB_vals, mA, tanb, alphaD, yq, ychi, n=400):
    """Cross sections vs bound-state mass mB_vals (mediator mass mA fixed).

    Returns (gg, qq, total) arrays [fb], K-factors included. Absolute scale carries
    the +2.7% normalization offset and has no PDF/scale uncertainty (module docstring).
    """
    mB_vals = np.asarray(mB_vals, dtype=float)
    gg = np.empty_like(mB_vals); qq = np.empty_like(mB_vals); tot = np.empty_like(mB_vals)
    for i, mB in enumerate(mB_vals):
        mDM = mdm_from_mb(mB, alphaD)
        k_gg, k_qq = Kgg(2 * mDM), Kqq(2 * mDM)
        gg[i] = cross_gg(pdf, mDM, mA, tanb, alphaD, yq, ychi, n=n) * k_gg
        qq[i] = cross_qq(pdf, mDM, mA, tanb, alphaD, yq, ychi, n=n) * k_qq
        tot[i] = gg[i] + qq[i]
    return gg, qq, tot


def scan_grid(pdf, mB_vals, mA_vals, tanb, alphaD, yq, ychi, n=400):
    """Total cross section over a (mB, mA) grid. Returns 2D array [len(mA), len(mB)]."""
    mB_vals = np.asarray(mB_vals, dtype=float)
    mA_vals = np.asarray(mA_vals, dtype=float)
    Z = np.empty((len(mA_vals), len(mB_vals)))
    for j, mA in enumerate(mA_vals):
        _, _, Z[j] = scan_mB(pdf, mB_vals, mA, tanb, alphaD, yq, ychi, n=n)
    return Z


# --- exposed pieces of the physics (for the explanatory plots) --------------
def binding_energy(mDM, alphaD):
    """1S Coulomb binding energy of the bound state [GeV]: (alphaD^2/4) mDM."""
    return 0.25 * alphaD ** 2 * mDM


def gamma_BS(mDM, alphaD):
    """Bound-state annihilation width [GeV]: (alphaD^5/4)(2 mDM)."""
    return (alphaD ** 5 / 4.0) * (2 * mDM)


def gamma_A_parts(mA, tanb, yq):
    """Heavy-quark partial widths of A: (Gamma_bb, Gamma_tt) [GeV]."""
    bb = (3.0 / (8 * np.pi)) * ((YB * yq) / np.sqrt(2)) ** 2 * tanb ** 2 * mA \
        * np.sqrt(1.0 - 4.0 * MB ** 2 / mA ** 2)
    tt = 0.0
    if mA > 2 * MT:
        tt = (3.0 / (8 * np.pi)) * ((YT * yq) / np.sqrt(2)) ** 2 * (1.0 / tanb ** 2) * mA \
            * np.sqrt(1.0 - 4.0 * MT ** 2 / mA ** 2)
    return bb, tt


def A12A_loop(tau):
    """Pseudoscalar triangle form factor A_{1/2}^A(tau) (complex). Public alias."""
    return A12A(tau)


def lambda_inv2(mB2, tanb, yq):
    """Effective ggA coupling-squared Lambda^-2 [GeV^-2] (top + bottom loops)."""
    return (ALPHA_S ** 2 * yq ** 2 / (16 ** 2 * (np.sqrt(2) * VEV) ** 2 * np.pi ** 2)) * \
        abs(A12A(mB2 / (4 * MT ** 2)) * (1.0 / tanb)
            + A12A(mB2 / (4 * MB ** 2)) * tanb) ** 2


# parton flavour groups for the luminosity plot. One orientation per channel so
# the cross-channel comparison is on a consistent basis (the luminosity integral
# is symmetric under x <-> tau0/x). "qqbar" is the light quarks d,u,s,c; bottom
# is shown separately as its own "bb" curve.
_LUMI_CHANNELS = {
    "gg": [(0, 0)],
    "bb": [(-5, 5)],
    "qqbar": [(-i, i) for i in range(1, 5)],
}


def parton_lumi(pdf, mhat, channel, n=400):
    """Parton luminosity entering sigma at partonic mass mhat [GeV].

    Sum over the channel's (a,b) flavour pairs of
    Integral xf_a(x) xf_b(tau0/x)/(x*tau0) dx, tau0 = mhat^2/s -- i.e. exactly
    the luminosity factor that multiplies the partonic cross section.
    channel in {'gg','qqbar','bb'}.
    """
    mB2 = mhat ** 2
    return sum(_lumi_integral(pdf, mB2, fa, fb, n=n) for fa, fb in _LUMI_CHANNELS[channel])


def scan_param(pdf, name, values, base, n=400):
    """Total cross section while sweeping one parameter.

    base is a dict with keys mDM, mA, tanb, alphaD, yq, ychi; `name` is the key
    to sweep over `values`. Returns a numpy array of sigma_total [fb].
    """
    out = np.empty(len(values))
    for i, v in enumerate(values):
        p = dict(base); p[name] = v
        out[i] = total_xsec(pdf, p["mDM"], p["mA"], p["tanb"], p["alphaD"],
                            p["yq"], p["ychi"], n=n)
    return out
