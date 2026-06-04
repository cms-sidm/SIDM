"""Shared analysis helpers for the dark-photon lifetime / proper-cτ study.

Loads the full-grid coffea and provides the cτ estimators used by the companion
notebooks:
  - the histogram mean of the proper-lxyz distribution (faithful regime),
  - the core-slope fit (exponential slope on the un-truncated low-x core),
  - the acceptance-corrected fit (folds in the lab-decay-length cap).

Each fit also returns its covariance, so the notebooks can draw a ±1σ error band
on the fitted curve and quote cτ with its uncertainty.

The grid was produced by _fullgrid_run.py and lives at
root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/coffea_outputs/<user>/lifetime_study/.
"""
import os
import math
import subprocess
from collections import defaultdict

import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import coffea.util

CHANNEL = "baseNoLj_noTrigger"   # gen-level channel, no HLT-trigger sculpting
EOS_COFFEA = (
    "root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/"
    "coffea_outputs/murtazas/lifetime_study/lifetime_fullgrid.coffea"
)


def load_grid(local_path="lifetime_fullgrid.coffea", eos_url=EOS_COFFEA):
    """Load the full-grid coffea; xrdcp it from the shared EOS area if not cached locally."""
    if not os.path.exists(local_path):
        subprocess.run(["xrdcp", "-f", eos_url, local_path], check=True)
    return coffea.util.load(local_path)


def ctau_cm(name):
    """Nominal proper cτ in cm, parsed from the sample-name token (in mm)."""
    return float(name.split("_")[-1].replace("mm", "").replace("p", ".")) / 10.0


def mass_gev(token):
    """Numeric mass in GeV from a name token, e.g. '500GeV' -> 500.0, '5p0GeV' -> 5.0."""
    return float(token.replace("GeV", "").replace("p", "."))


def mass_point(name):
    """(channel, BS-mass, Dp-mass) key, e.g. ('4Mu', '500GeV', '5p0GeV')."""
    p = name.split("_")
    return (p[0], p[1], p[2])


def mass_label(key):
    """Reviewer-facing label for a mass point, e.g. '4Mu  m(Phi)=500, m(A)=5.0 GeV'."""
    ch, bs, dp = key
    return f"{ch}   m($\\Phi$)={mass_gev(bs):g}, m(A)={mass_gev(dp):g} GeV"


def proper_density(output, sample, channel=CHANNEL):
    """Proper-lxyz dN/dx, effective counts (value**2/variance), and the histogram-mean cτ."""
    h = output[sample]["hists"]["genAs_lxyz_proper"][{"channel": channel}]
    edges = h.axes[-1].edges
    centers = h.axes[-1].centers
    vals = h.values()
    var = h.variances()
    with np.errstate(divide="ignore", invalid="ignore"):
        effN = np.where(var > 0, vals**2 / var, 0.0)
        dens = np.where(vals > 0, vals / np.diff(edges), np.nan)
    mean = (vals * centers).sum() / max(vals.sum(), 1)
    return centers, dens, effN, mean


def core_slope_fit(output, sample, channel=CHANNEL, k=3.0):
    """Fit log(dN/dx) = a - x/cτ on the un-truncated low-x core (cτ = -1/slope),
    iteratively trimming the systematic high-x deficit (the lab-cap truncation).

    Returns a fit dict with cτ and its uncertainty (propagated from the slope
    covariance), plus the linear log-model, its (slope, intercept) and 2x2
    covariance, the fit x-range, and the underlying distribution arrays."""
    centers, dens, effN, mean = proper_density(output, sample, channel)
    ly = np.where(dens > 0, np.log(dens), np.nan)
    mask = effN >= 5
    slope, ic, cov = np.nan, np.nan, np.full((2, 2), np.nan)
    if mask.sum() >= 4:
        hi_x = np.median(centers[mask])   # freeze the upper-half trim gate from the initial core
        for _ in range(30):
            # cov='unscaled': the weights are exact 1/sigma (Poisson), so the covariance is
            # the pure statistical one, not rescaled by reduced chi-square (which the
            # one-sided trim would make ill-defined)
            (slope, ic), cov = np.polyfit(
                centers[mask], ly[mask], 1, w=np.sqrt(effN[mask]), cov="unscaled")
            with np.errstate(divide="ignore", invalid="ignore"):
                sig = np.where(effN > 0, 1 / np.sqrt(effN), np.inf)
            trim = mask & (ly < ic + slope * centers - k * sig) & (centers > hi_x)
            if not trim.any():
                break
            mask = mask & ~trim
    ctau = -1 / slope if (np.isfinite(slope) and slope < 0) else np.nan
    # cτ = -1/slope  =>  d(cτ)/d(slope) = 1/slope**2 = cτ**2
    ctau_err = ctau**2 * math.sqrt(cov[0, 0]) if (np.isfinite(ctau) and np.isfinite(cov[0, 0])) else np.nan
    xfit = (centers[mask].min(), centers[mask].max()) if mask.any() else (np.nan, np.nan)
    return dict(kind="core-slope fit", ctau=ctau, ctau_err=ctau_err, mean=mean,
                logmodel=lambda x, s, b: b + s * x, popt=np.array([slope, ic]), pcov=cov,
                xfit=xfit, centers=centers, dens=dens, effN=effN)


def betagamma_cdf(output, anchor, channel=CHANNEL):
    """Empirical F_βγ and ⟨βγ⟩ from an (un-truncated) anchor sample's βγ histogram.

    The cumulative weight through bin i is the CDF at that bin's RIGHT EDGE, so F_βγ is
    interpolated on the bin edges (with a leading 0), not the centers — otherwise the CDF
    is shifted by half a bin and biases the fitted R_max."""
    h = output[anchor]["hists"]["genAs_betagamma"][{"channel": channel}]
    edges = h.axes[-1].edges
    w = h.values()
    tot = max(w.sum(), 1)
    cdf_edges = np.concatenate([[0.0], np.cumsum(w)]) / tot
    F_bg = interp1d(edges, cdf_edges, bounds_error=False, fill_value=(0.0, 1.0))
    return F_bg, (h.axes[-1].centers * w).sum() / tot


def acceptance_fit(output, sample, F_bg, mean, channel=CHANNEL, p0_ctau=None):
    """Fit log(dN/dx) = logA - x/cτ + log F_βγ(R_max/x) for (cτ, R_max).

    ε(x) = F_βγ(R_max/x) is the probability a decay with proper length x survives the
    lab-decay-length cap R_max, given the (lifetime-independent) boost distribution F_βγ.

    p0_ctau optionally seeds cτ — the core-slope value is less truncation-biased than the
    histogram mean (the bias being corrected here), so it is a safer start; it falls back
    to `mean`. absolute_sigma=True keeps the covariance purely statistical (the bin errors
    are exact Poisson).

    Returns a fit dict with cτ, R_max and their uncertainties (from the curve_fit
    covariance), plus the log-model, its parameters and full covariance, the fit
    x-range, and the underlying distribution arrays."""
    centers, dens, effN, _ = proper_density(output, sample, channel)
    good = effN >= 5
    nan_fit = dict(kind="acceptance-corrected fit", ctau=np.nan, ctau_err=np.nan,
                   Rmax=np.nan, Rmax_err=np.nan, logmodel=None, popt=None, pcov=None,
                   xfit=(np.nan, np.nan), centers=centers, dens=dens, effN=effN)
    if good.sum() < 4:
        return nan_fit
    x = centers[good]
    y = np.log(dens[good])
    yerr = 1 / np.sqrt(effN[good])

    def logmodel(x, logA, ctau, logR):
        return logA - x / ctau + np.log(np.clip(F_bg(np.exp(logR) / x), 1e-12, 1))

    seed = p0_ctau if (p0_ctau is not None and np.isfinite(p0_ctau) and p0_ctau > 0) else mean
    try:
        popt, pcov = curve_fit(
            logmodel, x, y, p0=[y.max(), seed, np.log(800.0)], sigma=yerr,
            absolute_sigma=True,
            bounds=([-50, 1e-6, np.log(10)], [50, 1e5, np.log(1e6)]), maxfev=40000,
        )
    except Exception:
        return nan_fit
    perr = np.sqrt(np.diag(pcov))
    Rmax = math.exp(popt[2])
    return dict(kind="acceptance-corrected fit", ctau=popt[1], ctau_err=perr[1],
                Rmax=Rmax, Rmax_err=Rmax * perr[2], logmodel=logmodel, popt=popt, pcov=pcov,
                xfit=(x.min(), x.max()), centers=centers, dens=dens, effN=effN)


def model_band(fit, xv):
    """Best-fit dN/dx and its ±1σ band on the grid xv, propagating the fit covariance
    through a numerical Jacobian of the log-model (works for both fit kinds)."""
    lm, popt, pcov = fit["logmodel"], fit["popt"], fit["pcov"]
    if lm is None or popt is None or pcov is None:
        nan = np.full_like(xv, np.nan, dtype=float)
        return nan, nan, nan
    f0 = lm(xv, *popt)
    J = np.zeros((len(xv), len(popt)))
    for j in range(len(popt)):
        step = 1e-4 * max(abs(popt[j]), 1e-3)
        dp, dm = np.array(popt, dtype=float), np.array(popt, dtype=float)
        dp[j] += step
        dm[j] -= step
        # central difference: captures R_max sensitivity at the turnover, which a forward
        # step into the saturated region (log F_βγ = 0) would miss
        J[:, j] = (lm(xv, *dp) - lm(xv, *dm)) / (2 * step)
    sd = np.sqrt(np.clip(np.einsum("ik,kl,il->i", J, pcov, J), 0, None))
    return np.exp(f0), np.exp(f0 - sd), np.exp(f0 + sd)


def group_by_mass_point(output):
    """{(channel, BS, Dp): [samples sorted by cτ]}."""
    g = defaultdict(list)
    for s in output:
        g[mass_point(s)].append(s)
    for k in g:
        g[k].sort(key=ctau_cm)
    return dict(g)


def compute_grid(output, channel=CHANNEL):
    """Per-sample nominal / mean / core-slope / acceptance-corrected cτ (each with its
    fit error) + R_max, and ⟨βγ⟩ per mass point.

    F_βγ for each mass point is taken from its shortest-cτ (faithful) anchor, since βγ
    depends on the production kinematics (masses), not the lifetime."""
    groups = group_by_mass_point(output)
    rows = {}
    bg_mean = {}
    for key, samples in groups.items():
        F_bg, bgm = betagamma_cdf(output, samples[0], channel)
        bg_mean[key] = bgm
        for s in samples:
            fB = core_slope_fit(output, s, channel)
            fC = acceptance_fit(output, s, F_bg, fB["mean"], channel, p0_ctau=fB["ctau"])
            rows[s] = dict(mass_point=key, nominal=ctau_cm(s), mean=fB["mean"],
                           core_slope=fB["ctau"], core_slope_err=fB["ctau_err"],
                           acceptance=fC["ctau"], acceptance_err=fC["ctau_err"],
                           Rmax=fC["Rmax"], Rmax_err=fC["Rmax_err"])
    return rows, groups, bg_mean


def betagamma_pct(output, anchor, p, channel=CHANNEL):
    """β-γ at percentile p (0..1) from the anchor sample's β-γ histogram — a robust
    high-boost scale. (βγ_max is a noisy single-event order statistic and, on this grid,
    is also pinned to the histogram's upper-edge axis ceiling at high boost.)"""
    h = output[anchor]["hists"]["genAs_betagamma"][{"channel": channel}]
    edges = h.axes[-1].edges
    w = h.values()
    cdf = np.cumsum(w) / max(w.sum(), 1)
    return edges[min(int(np.searchsorted(cdf, p)) + 1, len(edges) - 1)]


def insulated_core_fit(output, sample, x_insul, channel=CHANNEL):
    """Strictly-insulated core cτ: fit log(dN/dx) = a - x/cτ over x < x_insul only — the
    region where every event survives (ε=1, no truncation at all).

    DIAGNOSTIC, not a production estimator: once x_insul << cτ the window sits on the flat
    top of the exponential, so it has no lever arm and a huge error. Returns cτ, its error,
    the effective event count in the window (n), and the window's fitted log-fall in nats
    (lever_nats = (x_max - x_min)/cτ over the used window) as a lever-arm measure."""
    centers, dens, effN, _ = proper_density(output, sample, channel)
    sel = (centers < x_insul) & (effN >= 5)
    n = float(effN[sel].sum())
    if sel.sum() < 4:
        return dict(ctau=np.nan, ctau_err=np.nan, n=n, lever_nats=np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        (slope, ic), cov = np.polyfit(
            centers[sel], np.log(dens[sel]), 1, w=np.sqrt(effN[sel]), cov="unscaled")
    ctau = -1 / slope if (np.isfinite(slope) and slope < 0) else np.nan
    ctau_err = (ctau**2 * math.sqrt(cov[0, 0])
                if (np.isfinite(ctau) and np.isfinite(cov[0, 0])) else np.nan)
    lever_nats = (centers[sel].max() - centers[sel].min()) / ctau if np.isfinite(ctau) else np.nan
    return dict(ctau=ctau, ctau_err=ctau_err, n=n, lever_nats=lever_nats)


def plot_fit_grid(output, groups, phys_channel, kind, channel=CHANNEL, ncols=3):
    """One panel per mass point (cτ scan overlaid) showing every signal sample's
    proper-lxyz dN/dx with its fit + ±1σ band; cτ_fit ± error annotated per curve.

    kind: "core_slope" or "acceptance". Returns the Figure."""
    import matplotlib.pyplot as plt
    import mplhep as hep

    keys = sorted([k for k in groups if k[0] == phys_channel],
                  key=lambda k: (mass_gev(k[1]), mass_gev(k[2])))
    nrows = math.ceil(len(keys) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.6 * ncols, 5.8 * nrows), squeeze=False)
    fig.set_dpi(78)

    for idx, key in enumerate(keys):
        ax = axes[idx // ncols][idx % ncols]
        samples = groups[key]
        F_bg, _ = betagamma_cdf(output, samples[0], channel)
        colors = plt.cm.viridis(np.linspace(0, 0.92, len(samples)))
        for si, s in enumerate(samples):
            centers, dens, effN, mean = proper_density(output, s, channel)
            if kind == "core_slope":
                fit = core_slope_fit(output, s, channel)
            else:
                cs = core_slope_fit(output, s, channel)
                fit = acceptance_fit(output, s, F_bg, mean, channel, p0_ctau=cs["ctau"])
            good = effN > 0
            yerr = dens[good] / np.sqrt(effN[good])
            ax.errorbar(centers[good], dens[good], yerr=yerr, fmt="o", ms=3,
                        color=colors[si], alpha=0.45, lw=1, zorder=2)
            lab = f"{ctau_cm(s):g}"
            if np.isfinite(fit["ctau"]):
                xv = np.geomspace(*fit["xfit"], 200)
                m, lo, hi = model_band(fit, xv)
                ax.plot(xv, m, "-", color=colors[si], lw=2.2, zorder=3)
                ax.fill_between(xv, lo, hi, color=colors[si], alpha=0.2, lw=0, zorder=1)
                err = fit["ctau_err"]
                lab += f"  →  {fit['ctau']:.2g} ± {err:.1g}" if np.isfinite(err) \
                    else f"  →  {fit['ctau']:.2g}"
            ax.plot([], [], "-o", color=colors[si], ms=4, label=lab)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.minorticks_on()
        ax.set_xlabel("Proper decay length  $x = \\ell_{xyz}/\\beta\\gamma$  [cm]")
        ax.set_ylabel("dN/dx  [cm$^{-1}$]")
        ax.legend(title="cτ:  nom → fit ± err  [cm]", fontsize=9.5,
                  title_fontsize=10.5, loc="lower left", framealpha=0.9)
        hep.cms.label(ax=ax, data=False, fontsize=13)
        ax.text(0.97, 0.95, mass_label(key), transform=ax.transAxes,
                ha="right", va="top", fontsize=12)

    for idx in range(len(keys), nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)
    fig.tight_layout()
    return fig
