"""Module to define miscellaneous helper methods"""

import math
import os
import json
import warnings
import yaml
import numpy as np
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import scipy.optimize as opt
import hist
import hist.intervals
from sidm import BASE_DIR
import coffea.util
from coffea.processor import accumulate

def print_list(l):
    """Print one list element per line"""
    print('\n'.join(l))

def print_debug(name, val, print_mode=True):
    """Print variable name and value"""
    if print_mode:
        print(f"{name}: {val}")

def partition_list(l, condition):
    """Given a single list, return separate lists of elements that pass or fail a condition"""
    passes = []
    fails = []
    for x in l:
        if condition(x):
            passes.append(x)
        else:
            fails.append(x)
    return passes, fails

def flatten(x):
    """Flatten arbitrarily nested list or dict"""
    # https://stackoverflow.com/questions/2158395/
    flattened_list = []
    def loop(sublist):
        if isinstance(sublist, dict):
            sublist = sublist.values()
        for item in sublist:
            if isinstance(item, (dict, list)):
                loop(item)
            else:
                flattened_list.append(item)
    loop(x)
    return flattened_list

def add_unique_and_flatten(flattened_list, x):
    """Flatten arbitrarily nested list or dict, keeping only unique items"""
    # https://stackoverflow.com/questions/2158395/
    def loop(sublist):
        if isinstance(sublist, dict):
            sublist = sublist.values()
        for item in sublist:
            if isinstance(item, (dict, list)):
                loop(item)
            elif item not in flattened_list:
                flattened_list.append(item)
    loop(x)
    return flattened_list

def as_int(array):
    """Return array with values converted to ints"""
    return ak.values_astype(array, "int64")

def dR(obj1, obj2):
    """Return dR between obj1 and the nearest obj2; returns inf if no obj2 is found"""
    dr = obj1.nearest(obj2, return_metric=True)[1]
    return ak.fill_none(dr, np.inf)

def dR_general(obj1, obj2):
    """Return ΔR between obj1 and obj2, filling None with inf"""
    return ak.fill_none(obj1.delta_r(obj2), np.inf)

def dR_outer(obj1, obj2):
    """Return dR between outer tracks of obj1 and obj2"""
    return ak.fill_none(np.sqrt((obj1.outerEta - obj2.outerEta)**2 + (obj1.outerPhi - obj2.outerPhi)**2), np.inf)

def drop_none(obj):
    """Remove None entries from an array (not available in Awkward 1)"""
    return obj[~ak.is_none(obj, axis=1)] # fixme: not clear why axis=1 works and axis=-1 doesn't

def matched(obj1, obj2, r):
    """Return set of obj1 that have >=1 obj2 within r; remove None entries before returning"""
    return drop_none(obj1[dR(obj1, obj2) < r])

def add_matched_dsamuon_mass(obj):
    obj["mass"] = ak.full_like(obj.pt, 0.105712890625)
    return obj

def lj_combination_dR(obj):
    pair = ak.combinations(obj, 2, axis=1, fields=["lj1", "lj2"])
    dR = dR_general(pair["lj1"], pair["lj2"])
    min_dR = ak.min(dR_general(pair["lj1"], pair["lj2"]), axis=1)
    max_dR = ak.max(dR_general(pair["lj1"], pair["lj2"]), axis=1)
    return dR, min_dR, max_dR

def rho(obj, ref=None, use_v=False):
    """Return transverse distance between object and reference (default reference is 0,0)"""
    if use_v:
        obj_x = obj.vx
        obj_y = obj.vy
        ref_x = ref.vx if ref is not None else 0.0
        ref_y = ref.vy if ref is not None else 0.0
    else:
        obj_x = obj.x
        obj_y = obj.y
        ref_x = ref.x if ref is not None else 0.0
        ref_y = ref.y if ref is not None else 0.0
    return np.sqrt((obj_x - ref_x)**2 + (obj_y - ref_y)**2)

def dxy(obj, ref=None):
    """Return transverse distance between obj and ref at their point of closest approach"""
    # caveats discussed here apply: https://github.com/cms-sw/cmssw/blob/1bd97a649226ce2c2585f8b61f210aab6d0d4c44/DataFormats/TrackReco/interface/TrackBase.h#L678-L683
    shape = ak.ones_like(obj.vx)
    x_val = ak.flatten(ref.x) if ref is not None else 0.0
    y_val = ak.flatten(ref.y) if ref is not None else 0.0
    ref_x = x_val*shape
    ref_y = y_val*shape
    return (-(obj.vx - ref_x)*obj.py + (obj.vy - ref_y)*obj.px)/obj.pt

def lxy(obj):
    """Return transverse distance between production and decay vertices"""
    return rho(obj, ak.firsts(obj.children, axis=2), use_v=True)

def lxyz(obj):
    """Return 3D distance between an object's production and decay vertices"""
    ref = ak.firsts(obj.children, axis=2)
    return np.sqrt((obj.vx - ref.vx)**2 + (obj.vy - ref.vy)**2 + (obj.vz - ref.vz)**2)

def betagamma(obj):
    """Return beta*gamma (= p/m = pt*cosh(eta)/mass), guarding unphysical mass==0 gen copies"""
    mass = ak.where(obj.mass > 0, obj.mass, np.nan)
    return (obj.pt*np.cosh(obj.eta))/mass

def lxyz_proper(obj):
    """Return proper decay length: 3D flight distance deboosted by beta*gamma (= p/m)"""
    return lxyz(obj)/betagamma(obj)

def set_plot_style(style='cms', dpi=50):
    """Set plotting style using mplhep"""
    if style == 'cms':
        plt.style.use(hep.style.CMS)
    else:
        raise NotImplementedError
    plt.rcParams['figure.dpi'] = dpi

def plot(hists, skip_label=False, yerr=None, **kwargs):
    """Plot using hep.hist(2d)plot and add cms labels.

    yerr: error bars for a 1D histogram. **Pass this explicitly for any derived
        quantity** -- an efficiency, ratio, scale factor, or any hist filled via
        ``.view()[:] = ...`` rather than from event counts -- otherwise mplhep draws
        sqrt(bin_content) bars, which are meaningless for a non-count value.
        ``get_eff_hist`` returns the matching ``errors`` to pass here. Leave as None
        for an ordinary count histogram, where mplhep's default is correct:
        asymmetric Poisson (Garwood) bars for unweighted/unit-weight counts, and
        symmetric sqrt(sum w^2) for weighted hists.
    """

    # set default arguments
    default_kwargs = {
        'flow': "sum",
    }
    kwargs = {**default_kwargs, **kwargs}

    dim = len(hists[0].axes) if isinstance(hists, list) else len(hists.axes)
    if dim == 1:
        h = hep.histplot(hists, yerr=yerr, **kwargs)
    elif dim == 2:
        h = hep.hist2dplot(hists, **kwargs)
    else:
        raise NotImplementedError(f"Cannot plot {dim}-dimensional hist")
    if not skip_label:
        hep.cms.label()
    return h

def get_eff_hist(num_hist, denom_hist):
    """Returns the histogram of num_hist/denom_hist and a 2D numpy array of the up/down errors on the efficiency. Plot the errors using yerr=errors when plotting. """
    # make efficiency hist
    denom_vals = denom_hist.values()
    num_vals = num_hist.values()
    eff_values = num_vals/denom_vals
    eff_hist = hist.Hist(*num_hist.axes, storage=hist.storage.Weight())
    eff_hist.view().value[:] = eff_values

    # approximate weighted-hist errors based on avg (per-bin) weights
    num_counts = num_vals**2 / num_hist.variances()
    denom_counts = denom_vals**2 / denom_hist.variances()
    errors = hist.intervals.ratio_uncertainty(num_counts, denom_counts, 'efficiency')

    # Store the symmetrized Clopper-Pearson error^2 as the bin variance, so that
    # plotting this hist WITHOUT yerr=errors yields a right-magnitude (symmetric)
    # bar instead of mplhep's default sqrt(efficiency). Callers should still pass
    # yerr=errors for the correct asymmetric interval (e.g. plot_ratio does).
    sym_err = 0.5 * (errors[0] + errors[1])
    eff_hist.view().variance[:] = sym_err ** 2

    return eff_hist, errors

def load_yaml(cfg):
    """Load yaml files and return corresponding dict"""
    with open(cfg, encoding="utf8") as yaml_cfg:
        return yaml.safe_load(yaml_cfg)

def load_census_skip(census_skip):
    """Load a census skip-list into {sample: set(basenames-to-skip)}. census_skip is a path to a
    skip JSON (the `_skip.json` written by file_census.cleaned_filelists, or a promoted copy under
    sidm/configs/census/), or a bare name resolved there. The skip JSON is {"skip": {sample: [...]}}."""
    path = census_skip
    if not os.path.isfile(path):
        path = f"{BASE_DIR}/configs/census/{census_skip}"
    with open(path, encoding="utf8") as f:
        data = json.load(f)
    return {sample: set(files) for sample, files in data.get("skip", {}).items()}


def apply_census_veto(file_list, skip_basenames):
    """Additively drop files whose basename is in skip_basenames (ON TOP of any the location YAML
    already comments out -- the YAML commenting is honored upstream; this only removes more).
    Returns (kept_files, n_dropped)."""
    if not skip_basenames:
        return file_list, 0
    kept = [f for f in file_list if os.path.basename(f) not in skip_basenames]
    return kept, len(file_list) - len(kept)


def make_fileset(samples, ntuple_version, max_files=-1, location_cfg="signal_v8.yaml", fileset=None,
                 replace_xcache=False, census_skip=None):
    """Make fileset to pass to processor.runner.

    replace_xcache: when True, rewrite 'root://xcache//' to 'root://cmseos.fnal.gov//' so the
    fileset is usable from LPC (xcache is the coffea-casa cache and does not resolve elsewhere).
    census_skip: optional census skip-list (a path, or a name under sidm/configs/census/). Files
    the census flagged (bad / unreachable / empty -- see file_census) are dropped IN ADDITION to
    those the location YAML comments out; the YAML's manual commenting is always honored.
    """
    # assume location_cfg is stored in sidm/configs/ntuples/
    location_cfg = f"{BASE_DIR}/configs/ntuples/" + location_cfg
    locations = load_yaml(location_cfg)[ntuple_version]
    skip = load_census_skip(census_skip) if census_skip else {}
    if not fileset:
        fileset = {}
    for sample in samples:
        sample_yaml = locations["samples"][sample]
        base_path = locations["path"] + sample_yaml["path"]
        file_list = [base_path + f for f in sample_yaml["files"]]
        if replace_xcache:
            file_list = [f.replace("root://xcache//", "root://cmseos.fnal.gov//") for f in file_list]
        file_list, n_vetoed = apply_census_veto(file_list, skip.get(sample))
        if n_vetoed:
            print(f"  census veto: dropped {n_vetoed} flagged file(s) from {sample}")
        if max_files != -1:
            file_list = file_list[:max_files]
        fileset[sample] = {
            "files": file_list,
            "metadata": {
                "skim_factor": sample_yaml.get("skim_factor", 1.0),
                "is_data": sample_yaml.get("is_data", False),
                "year": sample_yaml.get("year", "2018"),
            },
        }
    return fileset

def check_bit(array, bit_num):
    """Return boolean stored in the bit_numth bit of array"""
    return (array & pow(2, bit_num)) > 0

def check_bits(array, bit_nums):
    result = (array & pow(2, bit_nums[0]))>0
    for x in bit_nums[1:]:
        result = (result & ((array & pow(2, x))>0))>0
    return result

def get_hist_mean(h):
    """Return mean of 1D histogram"""
    return np.atleast_1d(h.profile(axis=0).view())[0].value

def _gaussian_ratio(num_hist, den_hist):
    """Generic bin-by-bin ratio num/den with Gaussian-propagated errors, for ratios
    that are NOT efficiencies (num not a subset of den), where the binomial
    Clopper-Pearson interval used by get_eff_hist is invalid. Returns
    (ratio_hist, errors) with errors a symmetric (2, N) array for plot()."""
    n = num_hist.values()
    d = den_hist.values()
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(d != 0, n / d, np.nan)
        rel2 = (np.where(n != 0, num_hist.variances() / n**2, 0.0)
                + np.where(d != 0, den_hist.variances() / d**2, 0.0))
        rerr = np.nan_to_num(np.abs(ratio) * np.sqrt(rel2))
    ratio_hist = hist.Hist(*num_hist.axes, storage=hist.storage.Weight())
    ratio_hist.view().value[:] = np.nan_to_num(ratio)
    ratio_hist.view().variance[:] = rerr**2
    return ratio_hist, np.array([rerr, rerr])

def plot_ratio(num, den, **kwargs):
    # Lower-panel label/limits default to the efficiency convention but can be
    # overridden for a generic ratio; popped so they do not reach hep.histplot.
    ratio_ylabel = kwargs.pop("ratio_ylabel", "Efficiency")
    ratio_ylim = kwargs.pop("ratio_ylim", (0, 1.2))
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12, 12), sharex=True,
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}
    )
    plt.sca(ax1)
    label = kwargs["legend"][0] if "legend" in kwargs else None
    plot(den, flow='none', color="k", skip_label=True, label=label)

    if not isinstance(num, list):
        num = [num]

    for i, x in enumerate(num):
        plot(x, flow='none', label=label)

    if "legend" in kwargs:
        ax1.legend(title = kwargs["text"], alignment="left")

    if "ylim" in kwargs:
        plt.ylim(kwargs["ylim"])
    if "ylabel" in kwargs:
        plt.ylabel(kwargs["ylabel"])
    plt.tight_layout()
    plt.sca(ax2)
    for x in num:
        try:
            eff, errors = get_eff_hist(x, den)
        except ValueError:
            # num is not a subset of den (a shape ratio, not an efficiency): the
            # binomial Clopper-Pearson interval is invalid and ratio_uncertainty
            # raises -- fall back to a Gaussian-propagated ratio.
            eff, errors = _gaussian_ratio(x, den)
        plot(eff, histtype='errorbar', yerr=errors, skip_label=True)

    ax2.set_ylabel(ratio_ylabel)
    if ratio_ylim is not None:
        ax2.set_ylim(*ratio_ylim)

def round_sigfig(val, digits=1):
    """Return a number rounded to a given number of significant figures. Uses magic copied from
    https://stackoverflow.com/questions/3410976/how-to-round-a-number-to-significant-figures-in-python"""
    return float('{:g}'.format(float('{:.{p}g}'.format(val, p=digits))))

def proper_ctau(bs, zd, lab_ct, grid_cfg=f"{BASE_DIR}/configs/signal_grid.yaml"):
    """Convert average lab-frame transverse decay length in cm to proper decay
    length in mm for SIDM signals"""
    grid = load_yaml(grid_cfg)
    # handle goofy edge cases that I suspect stems from Weinan rounding errors
    if (float(bs), float(zd), float(lab_ct)) == (150, 0.25, 150):
        proper_ct = 6.7
    elif (float(bs), float(zd), float(lab_ct)) == (150, 5, 150):
        proper_ct = 130.0
    elif (float(bs), float(zd), float(lab_ct)) == (800, 0.25, 150):
        proper_ct = 1.2
    else:
        proper_ct = lab_ct/grid[bs][zd]["labframe_factor"]
    return round_sigfig(proper_ct, digits=2)

def lab_ctau(bs, zd, proper_ct, grid_cfg=f"{BASE_DIR}/configs/signal_grid.yaml"):
    """Convert proper decay length in mm to average lab-frame transverse decay
    length in cm for SIDM signals"""
    grid = load_yaml(grid_cfg)
    # handle goofy edge case that I suspect stems from Weinan rounding errors
    if (float(bs), float(zd), float(proper_ct)) == (150, 0.25, 6.7):
        lab_ct = 150.0
    elif (float(bs), float(zd), float(proper_ct)) == (150, 5, 130):
        lab_ct = 150.0
    elif (float(bs), float(zd), float(proper_ct)) == (800, 0.25, 1.2):
        lab_ct = 150.0
    else:
        lab_ct = proper_ct*grid[bs][zd]["labframe_factor"]
    return round_sigfig(lab_ct, digits=2)

def get_xs(dataset, cfg="cross_sections.yaml"):
    """Fetch dataset xs from cfg"""
    # assume location_cfg is stored in sidm/configs/
    xs_menu = load_yaml(f"{BASE_DIR}/configs/" + cfg)
    try:
        return xs_menu[dataset]
    except KeyError:
        if dataset.startswith(("2Mu2E", "4Mu")):
            print("Signal not in xs cfg, assuming 1fb")
            return 0.001
        else:
            raise

def get_lumi(year, cfg="run_periods.yaml"):
    """Fetch run period lumi from cfg"""
    # assume location_cfg is stored in sidm/configs/
    lumi_menu = load_yaml(f"{BASE_DIR}/configs/" + cfg)
    return lumi_menu[year]["lumi"]

def get_lumixs_weight(dataset, year, sum_weights):
    """Get weights to scale n_evts to lumi*xs"""
    # n_evts: sum of weights from processed events
    lumi = get_lumi(year)
    xs = get_xs(dataset)
    return lumi*xs/sum_weights

def check_variablePhoton(value, min_val=0b01):
    """
    Function to check if the variable is at least `min_val`
    """
    return value >= min_val

def select_numbersPhoton(number, var1, var2):
    """
    Function to select the numbers where each variable is at least 0b010 except for variable of choice
    """
    selected = True

    # number will have 14 bits (2 bits per each cut)
    # starting with MinPtCut at the LSB
    # and ending with PhoIsoWithEALinScalingCut at the MSB
    variables = [
        ('MinPtCut', 0),  # 2 bits for MinPtCut (LSB)
        ('PhoSCEtaMultiRangeCut', 2),
        ('PhoSingleTowerHadOverEmCut', 4),
        ('PhoFull5x5SigmaIEtaIEtaCut', 6),
        ('ChHadIsoWithEALinScalingCut', 8),
        ('NeuHadIsoWithEAQuadScalingCut', 10),
        ('PhoIsoWithEALinScalingCut', 12),
    ]

    # Check each variable except variable of choice
    for var, start_bit in variables:
        # Get the 2 bits corresponding to this variable
        value = (number >> start_bit) & 0b11  # Extract 2 bits
        if (var != f'{var1}') and (var != f'{var2}'):
            if not check_variablePhoton(value):
                selected = False
                break

    return selected

def returnBitMapTArrayPhoton(bitMap, var1, var2):
    tList = []
    for i in range(len(bitMap)):
        temp = []
        if len(bitMap[i]) == 0:
            tList.append(temp)
            continue
        for j in range(len(bitMap[i])):
            if select_numbersPhoton(bitMap[i][j], var1, var2):
                temp.append(True)
            else:
                temp.append(False)
        tList.append(temp)
    return ak.Array(tList)

def lepton_dxy_resolution(leptons, pvs, rank="all", diff=False):
    matched = leptons.matched_gen

    if rank == "all":
        valid = (~ak.is_none(matched)) & (matched.status == 1)
        leptons = leptons[valid]
        matched = matched[valid]

        dxy_gen = dxy(matched, ref=pvs)
        dxy_reco = leptons.dxy
        nonzero = dxy_gen != 0

        result = dxy_reco - dxy_gen
        return result if diff else result[nonzero] / dxy_gen[nonzero]

    # rank is an integer: one lepton per event
    enough = ak.num(leptons) > rank
    leptons = leptons[enough]
    pvs = pvs[enough]
    matched = leptons.matched_gen

    lep = leptons[:, rank]
    gen = matched[:, rank]
    dxy_gen = dxy(gen, ref=pvs)
    dxy_reco = lep.dxy
    result = dxy_reco - dxy_gen

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = ak.where((~ak.is_none(gen)) & (gen.status == 1) & (dxy_gen != 0), result / dxy_gen, np.nan)

    return result if diff else ratio

def spin1_model(x, A, alpha):
    """Physics model: dN/dCosTheta ~ A * (1 + alpha * cos^2(theta))"""
    return A * (1 + alpha * x**2)

def plot_and_fit_polarization(hist, ax=None, color='black', label_prefix="Data", fit_range=(0, 0.8), density=False):
    """
    Extracts data from a CosTheta histogram, fits the Spin-1 model,
    and plots Data + Fit + 1-Sigma Band.

    Args:
        hist: The Coffea histogram object
        ax: The matplotlib axis to plot on (creates new if None)
        color: Color for the markers and fit line
        label_prefix: String for the legend (e.g., "Gen Muons")
        fit_range: Tuple (min, max) to restrict the fit (avoiding acceptance effects)
        density (bool): If True, normalizes the histogram to unit area (probability density).
                        Errors are scaled correctly to preserve statistical significance.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    raw_counts = hist.values().flatten()
    edges = hist.axes[-1].edges
    centers = (edges[:-1] + edges[1:]) / 2

    min_len = min(len(raw_counts), len(centers))
    raw_counts = raw_counts[:min_len]
    centers = centers[:min_len]
    edges = edges[:min_len + 1]

    if density:
        widths = edges[1:] - edges[:-1]
        integral = np.sum(raw_counts * widths)
        scale_factor = 1.0 / (integral if integral > 0 else 1.0)
    else:
        scale_factor = 1.0

    y_values = raw_counts * scale_factor
    y_err = np.sqrt(raw_counts) * scale_factor
    y_err[raw_counts == 0] = scale_factor

    # Draw the SAME symmetric sqrt(N) errors the fit uses below, so the plotted
    # bars equal the curve_fit sigmas. (Passing the raw hist with yerr=True would
    # let mplhep draw its asymmetric Garwood interval, which the fit does not use;
    # density is applied here too so the markers match the fit input exactly.)
    hep.histplot(y_values, bins=edges, ax=ax, yerr=y_err, color=color,
                 histtype='errorbar', marker='o', markersize=4, capsize=2,
                 label=label_prefix)

    mask = (centers >= fit_range[0]) & (centers <= fit_range[1])
    x_fit = centers[mask]
    y_fit = y_values[mask]
    y_err_fit = y_err[mask]

    p0 = [np.max(y_fit), 0.5]

    try:
        popt, pcov = opt.curve_fit(
            spin1_model, x_fit, y_fit, sigma=y_err_fit, absolute_sigma=True, p0=p0
        )
        A_opt, alpha_opt = popt
        perr = np.sqrt(np.diag(pcov))

        x_model = np.linspace(0, 1, 100)
        y_model = spin1_model(x_model, *popt)

        label_fit = f"Fit ($\\alpha={alpha_opt:.2f} \\pm {perr[1]:.2f}$)"
        ax.plot(x_model, y_model, '-', color=color, linewidth=2, label=label_fit)

        # Confidence Band
        jac = np.vstack([1 + alpha_opt * x_model**2, A_opt * x_model**2]).T
        y_sigma = np.sqrt(np.sum((jac @ pcov) * jac, axis=1))

        ax.fill_between(x_model, y_model - y_sigma, y_model + y_sigma,
                        color=color, alpha=0.2)

        ax.axvline(fit_range[1], color=color, linestyle=':', alpha=0.3)

    except Exception as e:
        print(f"Fit failed for {label_prefix}: {e}")

    hep.cms.label()
    return ax

def gaussian_model(x, A, mu, sigma):
    """Standard Gaussian with a norm, mean, and sigma param"""
    return A * np.exp(-0.5 * ((x - mu) / sigma)**2)

def plot_and_fit_gaussian(hist, ax=None, color='black', label_prefix="Data", fit_range=(-3, 3), density=False):
    """
    Extracts data from a histogram, fits the standard Gaussian model,
    and plots Data + Fit.

    Args:
        hist: The Coffea histogram object
        ax: The matplotlib axis to plot on (creates new if None)
        color: Color for the markers and fit line
        label_prefix: String for the legend (e.g., "Gen Muons")
        fit_range: Tuple (min, max) to restrict the fit (avoiding acceptance effects)
        density (bool): If True, normalizes the histogram to unit area (probability density).
                        Errors are scaled correctly to preserve statistical significance.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    counts = hist.values().flatten()
    edges = hist.axes[-1].edges
    centers = (edges[:-1] + edges[1:]) / 2

    if density:
        widths = edges[1:] - edges[:-1]
        integral = np.sum(counts * widths)
        scale_factor = 1.0 / (integral if integral > 0 else 1.0)
    else:
        scale_factor = 1.0

    y_values = counts * scale_factor
    y_err = np.sqrt(counts) * scale_factor
    y_err[counts == 0] = scale_factor

    # Draw the same symmetric errors the fit uses below (see plot_and_fit_polarization).
    hep.histplot(y_values, bins=edges, ax=ax, yerr=y_err, color=color,
                 histtype='errorbar', marker='o', markersize=4, capsize=2,
                 label=label_prefix)

    mask = (centers >= fit_range[0]) & (centers <= fit_range[1])
    x_fit = centers[mask]
    y_fit = y_values[mask]
    y_err_fit = y_err[mask]

    if len(x_fit) > 0 and np.sum(y_fit) > 0:
        mean_guess = np.average(x_fit, weights=y_fit)
        sigma_guess = np.sqrt(np.average((x_fit - mean_guess)**2, weights=y_fit))
        amp_guess = np.max(y_fit)
    else:
        mean_guess, sigma_guess, amp_guess = 0, 1, 1

    p0 = [amp_guess, mean_guess, sigma_guess]
    try:
        popt, pcov = opt.curve_fit(
            gaussian_model, x_fit, y_fit, sigma=y_err_fit, absolute_sigma=True, p0=p0
        )
        A_opt, mu_opt, sigma_opt = popt

        x_model = np.linspace(edges[0], edges[-1], 200)
        y_model = gaussian_model(x_model, *popt)

        label_fit = rf"Fit: $\mu={mu_opt:.2f}, \sigma={abs(sigma_opt):.2f}$"
        ax.plot(x_model, y_model, '-', color=color, linewidth=2, label=label_fit)

    except Exception as e:
        print(f"Fit failed: {e}")

    hep.cms.label()
    return ax

def sum_hist(samples_list, folder_name):
    summed_out = None
    for x in samples_list:
        output = coffea.util.load(f"{folder_name}/{x}.coffea")
        hists = output["out"][x]["hists"]
        if summed_out is None:
            summed_out = hists.copy()
        else:
            summed_out = accumulate ([hists, summed_out])
    return summed_out

def cosAlpha(objs):
    pairs = ak.combinations(objs, 2, axis=1)
    v1, v2 = ak.unzip(pairs)
    return np.cos(v1.deltaangle(v2))

def get_cut_yield(out_dict, sample, cut_name, channel, weighted=True, cumulative=True):
    """Return a weighted or unweighted yield for a cutflow cut."""
    if sample not in out_dict:
        raise KeyError(f"Sample '{sample}' not found in output dictionary")
    if "cutflow" not in out_dict[sample]:
        raise KeyError(f"Sample '{sample}' has no 'cutflow' entry")
    if channel not in out_dict[sample]["cutflow"]:
        raise KeyError(f"Channel '{channel}' not found in sample '{sample}' cutflow")

    cf = out_dict[sample]["cutflow"][channel]

    if hasattr(cf, "rows"):
        if not cumulative:
            raise AttributeError(
                "SimpleCutflow stores cumulative sequential yields only; "
                "individual cut yields are not available."
            )
        if cut_name not in cf.rows:
            raise KeyError(f"Cut '{cut_name}' not found in sample '{sample}' for channel '{channel}'")
        key = "weighted" if weighted else "raw"
        return cf.rows[cut_name][key]

    flow = cf.flow if weighted else cf.unweighted_flow

    for elem in flow:
        if elem.cut == cut_name:
            if cumulative:
                return elem.n_all
            if not hasattr(elem, "n_ind"):
                raise AttributeError(
                    f"Cut '{cut_name}' in sample '{sample}' channel '{channel}' "
                    "does not store individual cut yield 'n_ind'"
                )
            return elem.n_ind

    raise KeyError(f"Cut '{cut_name}' not found in sample '{sample}' for channel '{channel}'")

def print_CR_yield(out_dict, signals, backgrounds, cut_name, channel, weighted=True, cumulative=True):
    """Print signal/background yields and signal contamination fraction S/(S+B)."""
    S = 0.0
    B = 0.0

    print("=== Signal yields ===")
    for sig in signals:
        y = get_cut_yield(
            out_dict,
            sig,
            cut_name=cut_name,
            channel=channel,
            weighted=weighted,
            cumulative=cumulative,
        )
        S += y
        print(f"{sig:35s} : {y:.6f}")

    print("\n=== Background yields ===")
    for bkg in backgrounds:
        y = get_cut_yield(
            out_dict,
            bkg,
            cut_name=cut_name,
            channel=channel,
            weighted=weighted,
            cumulative=cumulative,
        )
        B += y
        print(f"{bkg:35s} : {y:.6f}")

    print("\n=== Totals ===")
    print(f"S = {S:.6f}")
    print(f"B = {B:.6f}")

    if (S + B) > 0:
        frac = S / (S + B)
        print(f"S/(S+B) = {frac:.6e}")
    else:
        frac = None
        print("S/(S+B) undefined (S+B = 0)")

    return {
        "S": S,
        "B": B,
        "S_over_SplusB": frac,
    }

def extract_cutflow_yields(cutflow_obj, unweighted=False, cumulative=True):
    """Return cut names and yields from a Cutflow or SimpleCutflow object."""
    if hasattr(cutflow_obj, "rows"):
        if not cumulative:
            raise AttributeError(
                "SimpleCutflow stores cumulative sequential yields only; "
                "individual cut yields are not available."
            )
        key = "raw" if unweighted else "weighted"
        cut_names = list(cutflow_obj.rows.keys())
        yields = [cutflow_obj.rows[cut][key] for cut in cut_names]
        return cut_names, yields

    flow = cutflow_obj.unweighted_flow if unweighted else cutflow_obj.flow

    cut_names = []
    yields = []
    for elem in flow:
        cut_names.append(elem.cut)
        if cumulative:
            yields.append(elem.n_all)
        else:
            if not hasattr(elem, "n_ind"):
                raise AttributeError(
                    f"Cut '{elem.cut}' does not store individual cut yield 'n_ind'"
                )
            yields.append(elem.n_ind)

    return cut_names, yields

def plot_cutflow_panels(
    out,
    sample_groups,
    channel,
    titles=None,
    unweighted=False,
    cumulative=True,
    normalize=False,
    logy=False,
    figsize=None,
    ncols=None,
    ncol_legend=1,
    linewidth=1.8,
    markersize=4,
    sharey=True,
    show=True,
):
    """Plot cutflow yields for an arbitrary number of sample groups."""
    n_panels = len(sample_groups)
    if n_panels == 0:
        raise ValueError("sample_groups must contain at least one panel")

    if titles is None:
        titles = [f"Panel {i+1}" for i in range(n_panels)]
    if len(titles) != n_panels:
        raise ValueError("titles must have the same length as sample_groups")

    if ncols is None:
        ncols = min(3, n_panels)
    nrows = int(np.ceil(n_panels / ncols))

    if figsize is None:
        figsize = (12 * ncols, 8 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharey=sharey)
    axes = np.atleast_1d(axes).flatten()

    global_cut_names = None

    for ax, samples, panel_title in zip(axes, sample_groups, titles):
        if len(samples) == 0:
            ax.set_title(panel_title)
            ax.axis("off")
            continue

        panel_cut_names = None

        for sample in samples:
            cutflow_obj = out[sample]["cutflow"][channel]
            cut_names, yields = extract_cutflow_yields(
                cutflow_obj,
                unweighted=unweighted,
                cumulative=cumulative,
            )

            y = np.array(yields, dtype=float)

            if normalize:
                if y[0] != 0:
                    y = y / y[0]
                else:
                    y = np.zeros_like(y)

            if panel_cut_names is None:
                panel_cut_names = cut_names
            elif cut_names != panel_cut_names:
                raise ValueError(
                    f"Cut order mismatch in panel '{panel_title}' for sample '{sample}'"
                )

            if global_cut_names is None:
                global_cut_names = cut_names
            elif cut_names != global_cut_names:
                raise ValueError(
                    f"Global cut order mismatch in panel '{panel_title}' for sample '{sample}'"
                )

            ax.plot(
                range(len(cut_names)),
                y,
                marker="o",
                linestyle="--",
                linewidth=linewidth,
                markersize=markersize,
                label=sample,
            )

        ax.set_title(panel_title)
        ax.set_xticks(range(len(panel_cut_names)))
        ax.set_xticklabels(panel_cut_names, rotation=90)
        ax.grid(True, axis="y", alpha=0.3)

        if logy:
            ax.set_yscale("log")

        ax.legend(fontsize=20, ncol=ncol_legend)

    for ax in axes[n_panels:]:
        ax.axis("off")

    ylabel = "Efficiency" if normalize else "all cut N" if cumulative else "individual cut N"
    for row_start in range(0, n_panels, ncols):
        axes[row_start].set_ylabel(ylabel)

    plt.tight_layout()
    if show:
        plt.show()

    return fig, axes[:n_panels]

def get_hist_1d_from_variable(out_sample, variable, channel):
    """Return a 1D hist, projecting eta/phi from a corresponding 2D eta_phi hist if needed."""
    hists = out_sample["hists"]

    if variable in hists:
        return hists[variable][channel, :]

    if variable.endswith("_eta"):
        base_variable = variable.removesuffix("_eta") + "_eta_phi"
        if base_variable not in hists:
            raise KeyError(
                f"Cannot find {variable!r} or corresponding 2D hist {base_variable!r}."
            )
        return hists[base_variable][channel, :, :].project(0)

    if variable.endswith("_phi"):
        base_variable = variable.removesuffix("_phi") + "_eta_phi"
        if base_variable not in hists:
            raise KeyError(
                f"Cannot find {variable!r} or corresponding 2D hist {base_variable!r}."
            )
        return hists[base_variable][channel, :, :].project(1)

    raise KeyError(
        f"Cannot find variable {variable!r}, and it does not look like an eta/phi projection."
    )

def sum_hists_projectable(out, samples, variable, channel):
    """Sum sample histograms for normal 1D variables and projected eta/phi variables."""
    h_sum = None

    for sample in samples:
        if sample not in out:
            raise KeyError(f"Sample {sample!r} not found in output dictionary")

        h = get_hist_1d_from_variable(out[sample], variable, channel)

        if h_sum is None:
            h_sum = h.copy()
        else:
            h_sum += h

    return h_sum

def plot_MC_sig_vs_bkg_panels(
    out,
    channels,
    backgrounds,
    signal_maps,
    variable=None,
    channel_idx=0,
    xlabel=None,
    masses=None,
    background_colors=None,
    signal_colors=None,
    signal_linestyles=None,
    ncols=3,
    figsize=None,
    yscale="log",
    ylim=(1e-3, 1e9),
    show=True,
):
    """Plot stacked MC backgrounds and signal overlays in panels grouped by mass.

    Warns if a background's per-bin variances are integers, which means it looks
    unweighted -- its bars are then asymmetric Poisson rather than sqrt(sum w^2) and its
    normalization is probably not lumi*xsec-weighted. See the README, "Error bars on plots".
    """
    if variable is None:
        raise ValueError("variable must be specified")
    if len(backgrounds) == 0:
        raise ValueError("backgrounds must contain at least one background group")
    if len(signal_maps) == 0:
        raise ValueError("signal_maps must contain at least one signal group")

    channel = channels[channel_idx]

    if masses is None:
        mass_set = set()
        for sample_map in signal_maps.values():
            mass_set.update(sample_map.keys())
        masses = sorted(mass_set)
    if len(masses) == 0:
        raise ValueError("masses must contain at least one mass point")

    n_panels = len(masses)
    ncols = min(ncols, n_panels)
    nrows = int(np.ceil(n_panels / ncols))

    if figsize is None:
        figsize = (12 * ncols, 10 * nrows)

    background_hists = []
    background_labels = []
    for label, samples in backgrounds.items():
        h = sum_hists_projectable(out, samples, variable, channel)
        if h is None:
            continue
        background_hists.append(h)
        background_labels.append(label)

    if len(background_hists) == 0:
        raise ValueError("No background histograms were found to plot")

    # Warn if a background MC hist looks unweighted (integer variances): the yerr=True
    # below then draws asymmetric Poisson (Garwood) bars rather than sqrt(sum w^2), and
    # the normalization is probably not lumi*xsec-weighted. Catches a real MC mistake.
    for _lbl, _h in zip(background_labels, background_hists):
        _v = _h.variances()
        if _v is not None and np.allclose(_v, np.round(_v)):
            warnings.warn(
                f"plot_MC_sig_vs_bkg_panels: background '{_lbl}' has integer variances "
                f"(looks unweighted) -- its error bars will be asymmetric Poisson, not "
                f"sqrt(sum w^2); check that it is lumi*xsec-weighted.",
                stacklevel=2,
            )

    background_plot_kwargs = {}
    if background_colors is not None:
        background_plot_kwargs["color"] = [background_colors[label] for label in background_labels]

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_1d(axes).flatten()

    for i, mass in enumerate(masses):
        ax = axes[i]
        plt.sca(ax)

        plot(
            background_hists,
            density=False,
            histtype="fill",
            label=background_labels,
            alpha=0.6,
            stack=True,
            yerr=True,
            linewidth=1,
            **background_plot_kwargs,
        )

        for signal_label, signal_map in signal_maps.items():
            if mass not in signal_map:
                continue

            linestyle = None
            if signal_linestyles is not None:
                linestyle = signal_linestyles.get(signal_label)

            for j, sample in enumerate(signal_map[mass]):
                signal_plot_kwargs = {}
                if signal_colors is not None:
                    signal_plot_kwargs["color"] = [signal_colors[j % len(signal_colors)]]
                if linestyle is not None:
                    signal_plot_kwargs["linestyle"] = linestyle

                signal_name = sample.replace(f"{signal_label}_", "", 1)
                mass_prefix = f"{mass}GeV_"
                if signal_name.startswith(mass_prefix):
                    signal_name = signal_name[len(mass_prefix):]

                plot(
                    get_hist_1d_from_variable(out[sample], variable, channel),
                    density=False,
                    yerr=True,
                    linewidth=3,
                    label=f"{signal_name} ({signal_label})",
                    **signal_plot_kwargs,
                )

        if yscale is not None:
            plt.yscale(yscale)
        if ylim is not None:
            plt.ylim(*ylim)

        plt.text(0.02, 0.92, f"BS:{mass}GeV", fontsize=30, transform=ax.transAxes)

        if xlabel is not None:
            plt.xlabel(xlabel)
        elif variable.endswith("_eta"):
            plt.xlabel(r"$\eta$")
        elif variable.endswith("_phi"):
            plt.xlabel(r"$\phi$")
        else:
            plt.xlabel(variable)

        plt.ylabel("Events")
        plt.legend(fontsize=20)

    for ax in axes[n_panels:]:
        ax.axis("off")

    plt.tight_layout()
    if show:
        plt.show()

    return fig, axes[:n_panels]

def fold_hist_flow(h):
    """Fold underflow and overflow bins into the first and last regular bins."""
    h_fold = h.copy()

    vals = h.values(flow=True)
    vars_ = h.variances(flow=True)

    new_vals = vals[1:-1].copy()
    new_vals[0] += vals[0]
    new_vals[-1] += vals[-1]

    view = h_fold.view()
    view.value[...] = new_vals

    if vars_ is not None:
        new_vars = vars_[1:-1].copy()
        new_vars[0] += vars_[0]
        new_vars[-1] += vars_[-1]
        view.variance[...] = new_vars

    return h_fold

def auto_rebin_hist(h, target_bins=50, rebin_threshold=50):
    """Automatically rebin a histogram when the number of bins is larger than rebin_threshold."""
    nbins = h.axes[0].size

    if nbins <= rebin_threshold:
        return h

    min_factor = math.ceil(nbins / target_bins)

    rebin_factor = None
    for factor in range(min_factor, nbins + 1):
        if nbins % factor == 0:
            rebin_factor = factor
            break

    if rebin_factor is None:
        print(f"[auto_rebin_hist] Could not find valid rebin factor for {nbins} bins. Skipping rebin.")
        return h

    return h[::hist.rebin(rebin_factor)]

def plot_data_mc(
    out,
    out_data,
    variable,
    channel,
    data_samples,
    backgrounds,
    signal_maps=None,
    year=2018,
    lumi=21.89,
    logy=True,
    xlabel=None,
    background_colors=None,
    signal_colors=["r", "orange", "gold", "lime", "cyan", "magenta"],
    signal_linestyles=None,
    signal_scale=1.0,
    fold_flow=True,
    auto_rebin=True,
    target_bins=50,
    rebin_threshold=50,
    mc_lumi_scale=1.0,
    ratio_ylim=(0.5, 1.5),
    yscale_min=1e-1,
    yscale_max=1e10,
    figsize=(18, 18),
    legend_ncol=3,
    legend_fontsize=19,
    show=True,
):
    """Plot data vs stacked MC backgrounds with optional signal overlays."""
    if not hasattr(hep, "comp") or not hasattr(hep.comp, "data_model"):
        raise ImportError(
            "plot_data_mc requires mplhep with hep.comp.data_model. "
            "Please install mplhep>=1.2.0 and restart the kernel."
        )

    if len(data_samples) == 0:
        raise ValueError("data_samples must contain at least one data sample")
    if len(backgrounds) == 0:
        raise ValueError("backgrounds must contain at least one background group")

    def prepare_hist(h):
        if fold_flow:
            h = fold_hist_flow(h)
        if auto_rebin:
            h = auto_rebin_hist(
                h,
                target_bins=target_bins,
                rebin_threshold=rebin_threshold,
            )
        return h

    background_hists = []
    background_labels = []
    for label, samples in backgrounds.items():
        h = sum_hists_projectable(out, samples, variable, channel)
        h = prepare_hist(h)
        h *= mc_lumi_scale
        background_hists.append(h)
        background_labels.append(label)

    h_data = None
    for sample in data_samples:
        if sample not in out_data:
            raise KeyError(f"Data sample {sample!r} not found in data output dictionary")

        h = get_hist_1d_from_variable(out_data[sample], variable, channel)
        if h_data is None:
            h_data = h.copy()
        else:
            h_data += h

    h_data = prepare_hist(h_data)

    background_plot_kwargs = {}
    if background_colors is not None:
        background_plot_kwargs["stacked_colors"] = [
            background_colors[label] for label in background_labels
        ]

    fig, ax_main, ax_ratio = hep.comp.data_model(
        data_hist=h_data,
        stacked_components=background_hists,
        stacked_labels=background_labels,
        xlabel=xlabel if xlabel is not None else variable,
        ylabel="Events",
        model_uncertainty=True,
        stacked_kwargs={"alpha": 0.6},
        **background_plot_kwargs,
    )

    if signal_maps is not None:
        default_colors = ["r", "orange", "gold", "lime", "cyan", "magenta"]
        colors = signal_colors if signal_colors is not None else default_colors
        color_idx = 0

        for signal_label, samples in signal_maps.items():
            linestyle = None
            if signal_linestyles is not None:
                linestyle = signal_linestyles.get(signal_label)

            for sample in samples:
                h_sig = get_hist_1d_from_variable(out[sample], variable, channel).copy()
                h_sig = prepare_hist(h_sig)
                h_sig *= signal_scale
                h_sig *= mc_lumi_scale

                signal_name = sample.replace(f"{signal_label}_", "", 1)
                signal_name = signal_name.replace("_", ", ")

                signal_plot_kwargs = {}
                if linestyle is not None:
                    signal_plot_kwargs["linestyle"] = linestyle

                hep.histplot(
                    h_sig,
                    ax=ax_main,
                    histtype="step",
                    linewidth=3,
                    color=colors[color_idx % len(colors)],
                    label=f"{signal_label} ({signal_name})",
                    **signal_plot_kwargs,
                )
                color_idx += 1

    hep.cms.label(ax=ax_main, data=True, year=year, lumi=lumi)

    if logy:
        ax_main.set_yscale("log")
        ax_main.set_ylim(yscale_min, yscale_max)

    ax_ratio.set_ylabel("Data/MC")
    ax_ratio.set_ylim(*ratio_ylim)

    ax_main.legend(
        loc="upper left",
        ncol=legend_ncol,
        fontsize=legend_fontsize,
    )

    fig.set_size_inches(*figsize)
    fig.subplots_adjust(hspace=0.05)

    if show:
        plt.show()

    return fig, ax_main, ax_ratio
def get_pairs(obj):
    pairs = ak.combinations(obj, 2, axis=1)
    return (pairs)
