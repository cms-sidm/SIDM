"""Module to define miscellaneous helper methods"""

import yaml
import numpy as np
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import scipy.optimize as opt
import hist.intervals
from sidm import BASE_DIR
import coffea.util

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

def set_plot_style(style='cms', dpi=50):
    """Set plotting style using mplhep"""
    if style == 'cms':
        plt.style.use(hep.style.CMS)
    else:
        raise NotImplementedError
    plt.rcParams['figure.dpi'] = dpi

def plot(hists, skip_label=False, **kwargs):
    """Plot using hep.hist(2d)plot and add cms labels"""

    # set default arguments
    default_kwargs = {
        'flow': "sum",
    }
    kwargs = {**default_kwargs, **kwargs}

    dim = len(hists[0].axes) if isinstance(hists, list) else len(hists.axes)
    if dim == 1:
        h = hep.histplot(hists, **kwargs)
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
    eff_hist = hist.Hist(*num_hist.axes)
    eff_hist.values()[:] = eff_values

    # approximate weighted-hist errors based on avg (per-bin) weights
    num_counts = num_vals**2 / num_hist.variances()
    denom_counts = denom_vals**2 / denom_hist.variances()
    errors = hist.intervals.ratio_uncertainty(num_counts, denom_counts, 'efficiency')

    return eff_hist, errors

def load_yaml(cfg):
    """Load yaml files and return corresponding dict"""
    with open(cfg, encoding="utf8") as yaml_cfg:
        return yaml.safe_load(yaml_cfg)

def make_fileset(samples, ntuple_version, max_files=-1, location_cfg="signal_v8.yaml", fileset=None):
    """Make fileset to pass to processor.runner"""
    # assume location_cfg is stored in sidm/configs/ntuples/
    location_cfg = f"{BASE_DIR}/configs/ntuples/" + location_cfg
    locations = load_yaml(location_cfg)[ntuple_version]
    if not fileset:
        fileset = {}
    for sample in samples:
        sample_yaml = locations["samples"][sample]
        base_path = locations["path"] + sample_yaml["path"]
        file_list = [base_path + f for f in sample_yaml["files"]]
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

def plot_ratio(num, den, **kwargs):
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
        eff, errors = get_eff_hist(x, den)
        plot(eff, histtype='errorbar', yerr=errors, skip_label=True)

    ax2.set_ylabel("Efficiency")
    ax2.set_ylim(0, 1.2)

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

    hep.histplot(hist, ax=ax, yerr=True, density=density, color=color, 
                 histtype='errorbar', marker='o', markersize=4, capsize=2, 
                 label=label_prefix)

    raw_counts = hist.values().flatten()
    edges = hist.axes[-1].edges
    centers = (edges[:-1] + edges[1:]) / 2

    min_len = min(len(raw_counts), len(centers))
    raw_counts = raw_counts[:min_len]
    centers = centers[:min_len]

    if density:
        widths = edges[1:] - edges[:-1]
        widths = widths[:min_len]
        integral = np.sum(raw_counts * widths)
        scale_factor = 1.0 / (integral if integral > 0 else 1.0)
    else:
        scale_factor = 1.0
        
    y_values = raw_counts * scale_factor
    y_err = np.sqrt(raw_counts) * scale_factor 
    y_err[y_err == 0] = scale_factor 

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

    hep.histplot(hist, ax=ax, yerr=True, density=density, color=color, 
                 histtype='errorbar', marker='o', markersize=4, capsize=2, 
                 label=label_prefix)

    counts = hist.values().flatten()
    edges = hist.axes[-1].edges
    centers = (edges[:-1] + edges[1:]) / 2

    if density:
        widths = edges[1:] - edges[:-1]
        integral = np.sum(counts * widths)
        scale_factor = 1.0 / (integral if integral > 0 else 1.0)
    else:
        scale_factor = 1.0

    y_fit = counts * scale_factor
    y_err_fit = np.sqrt(counts) * scale_factor
    y_err_fit[y_err_fit == 0] = scale_factor

    mask = (centers >= fit_range[0]) & (centers <= fit_range[1])
    x_fit = centers[mask]
    y_fit = y_fit[mask]
    y_err_fit = y_err_fit[mask]

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

def sum_hist(bgs_list, histogram_name, channel_name, folder_name=None, ratio=1, binning=1j):
    """
    Sum histograms for a list of backgrounds/data
    For data give ratio =1
    When we run all the data or compare background with all data, give ratio=1,
    """
    summed_hist = None
    for x in bgs_list:
        output_bg = coffea.util.load(f"{folder_name}/{x}.coffea")
        try:
            hist = output_bg["out"][x]["hists"][histogram_name][channel_name, ::binning]
        except:#for 2D histograms
            hist = output_bg["out"][x]["hists"][histogram_name][channel_name, ::binning, ::binning]
        if summed_hist is None:
            summed_hist = hist.copy()
        else:
            summed_hist += hist
    return summed_hist
