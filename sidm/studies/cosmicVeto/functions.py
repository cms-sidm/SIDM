import awkward as ak
import coffea.util
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import yaml

## get list of data file names
def get_data_files():
    yaml_file_path = '../../configs/ntuples/data_skimmed.yaml'
    with open(yaml_file_path, 'r') as file:
        yaml_data = yaml.safe_load(file)
    data_files = list(yaml_data["llpNanoAOD_v2"]["samples"].keys())
    return data_files

## get the list of bg file name
def get_bg_file_names(name):
    yaml_file_path = '../../configs/composite_samples.yaml'
    with open(yaml_file_path, 'r') as file:
        yaml_data = yaml.safe_load(file)
    file_names = list(yaml_data[name])
    return file_names

# get the list of 2mu2e signal
def get_signal_2mu2e_list():
    yaml_file_path = '../../configs/ntuples/signal_2mu2e_v10.yaml'
    with open(yaml_file_path, 'r') as file:
        data = yaml.safe_load(file)
    signals = list(data["llpNanoAOD_v2"]["samples"].keys())
    return signals

## get the list of 4mu signals
def get_signal_4mu_list():
    yaml_file_path = '../../configs/ntuples/signal_4mu_v10.yaml'
    with open(yaml_file_path, 'r') as file:
        data = yaml.safe_load(file)
    signals = list(data["llpNanoAOD_v2"]["samples"].keys())
    return signals


## load output files from Output folder for each study
def load_output(file_name, channel_name):
    output = coffea.util.load(f"OutputFiles/{channel_name}/{file_name}.coffea")
    return output

## get the cuts, cuflow list for easy plotting
def get_cutflow_list(cutflow, raw=False):
    rows = cutflow.rows
    cuts = rows.keys()
    N = []
    for x in cuts:
        if raw:
            N.append(rows[x]["raw"])
        else:
            N.append(rows[x]["weighted"])
    return cuts, N

## summing cutflow for backgrounds
def sum_bg_cutflow(bg_list, channel_name, raw=False):
    summed_cutflow = None
    cut = None
    for s in bg_list:
        output = load_output(s, channel_name)
        cutflow = output["out"][s]["cutflow"][channel_name]
        if summed_cutflow is None:
            summed_cutflow = cutflow
        else:
            summed_cutflow += cutflow
    return summed_cutflow

## get a plotting label for each signal
def get_signal_label(signal):
    allowed_lxy = [0.3, 3, 30, 150, 300]

    parts = signal.split("_")
    if parts[0] == "2Mu2E":
        label =r"2$\mu2e$"
    else:
        label =r"4$\mu$"
    mass = parts[1]
    zd_mass = parts[2].replace("p", ".")
    ctau = float(parts[3].replace("p", ".").replace("mm", ""))
    yaml_file_path = '../../configs/signal_grid.yaml'
    with open(yaml_file_path, 'r') as file:
        data = yaml.safe_load(file)
    factor = data[int(mass.split("G")[0])][float(zd_mass.split("G")[0])]["labframe_factor"]
    lxy = ctau *factor
    closest_lxy = min(allowed_lxy, key=lambda x: abs(x - lxy))  
    label  = f"{label}, {mass}, {zd_mass}, {closest_lxy}cm"
    return label

def isM(histogram):
    """Return a copy of histogram with negative bin contents made positive.

    NOTE: negative bins are physical -- they come from negative MC event
    weights -- so flipping their sign biases the background upwards. Prefer
    leaving them alone, or clipping to zero, unless the plotting backend
    genuinely cannot take them.
    """
    values = histogram.values().copy()
    isMinus = values < 0
    if not np.any(isMinus):
        return histogram
    print("found negative values, setting them positive")
    values[isMinus] = abs(values[isMinus])
    hist_corrected = histogram.copy()
    hist_corrected.values()[...] = values
    return hist_corrected

## plots data/MC for all the histograms provided in a list
def plot_DataMC(histogram_list, channel_name,
                summed_QCD, 
                summed_TT, summed_DY, 
                summed_DB, summed_data, 
                ratio, binning=1j, 
                ranges=None, file_name = None,
               title = None):
    columns=  len(histogram_list)
    fig, axs = plt.subplots(2, columns,figsize=(15*columns, 15),
    gridspec_kw={"height_ratios": [3, 1],"hspace": 0.05,})
    for i, histogram_name in enumerate(histogram_list):
        print(histogram_name)
        if columns == 1:
            ax_main = axs[0]
            ax_comp = axs[1]
        else:
            ax_main = axs[0, i]
            ax_comp = axs[1, i]

        if ranges:
            sum_bg_qcd = ratio * isM(summed_QCD[histogram_name][channel_name, ranges[0]:ranges[1]:binning])
            sum_bg_tt = ratio * isM(summed_TT[histogram_name][channel_name,ranges[0]:ranges[1]:binning])
            sum_bg_dy = ratio * isM(summed_DY[histogram_name][channel_name,ranges[0]:ranges[1]:binning])
            sum_bg_db = ratio * isM(summed_DB[histogram_name][channel_name,ranges[0]:ranges[1]:binning])
            sum_data = summed_data[histogram_name][channel_name,ranges[0]:ranges[1]:binning]
        else:
            sum_bg_qcd = ratio * isM(summed_QCD[histogram_name][channel_name, ::binning])
            sum_bg_tt = ratio * isM(summed_TT[histogram_name][channel_name, ::binning])
            sum_bg_dy = ratio * isM(summed_DY[histogram_name][channel_name, ::binning])
            sum_bg_db = ratio * isM(summed_DB[histogram_name][channel_name, ::binning])
            sum_data = summed_data[histogram_name][channel_name, ::binning]
        hep.comp.data_model(
        data_hist=sum_data,
        stacked_components=[sum_bg_qcd, sum_bg_tt, sum_bg_dy, sum_bg_db],
        stacked_labels=["QCD", "TT","DY", "Diboson"],
        xlabel=sum_data.axes[0].label,
        ylabel="Events",
        data_w2method="poisson",
        fig=fig,
        ax_main=ax_main,
        ax_comparison=ax_comp,
         )
        hep.cms.label(data=True, lumi=59.83, ax=ax_main)
        ax_main.set_yscale("log")
        if title:
            ax_main.legend(title=title)
        else:
            ax_main.legend()
    if file_name:
        plt.savefig(f"Plots/{file_name}.png")
    plt.show()
    plt.close()

## get 1D histograms from 2D
def plot_1d (eta_phi_histogram, axs_name):
    histogram_1d = eta_phi_histogram.project(axs_name)
    return(histogram_1d)

def plot_fraction_less (signals, histogram_name, channel_name,
                        output_signal_2mu, output_signal_4mu,
                        thresholds, colors, ylabel=None):
    """plot fraction error when we use less than or equal"""
    plt.figure(figsize =(15, 12))
    for i, s in enumerate(signals):
        label = get_signal_label(s)
        if label[0] == "2":
            histogram = output_signal_2mu["out"][s]["hists"][histogram_name][channel_name, :]
        elif label[0] == "4":
            histogram = output_signal_4mu["out"][s]["hists"][histogram_name][channel_name, :]
        else:
            raise ValueError(f"cannot tell 2mu2e and 4mu apart for sample {s}")
        values = histogram.values()
        values_flow = histogram.values(flow=True)
        overflow = values_flow[-1]
        underflow = values_flow[0]
        bin_edges = histogram.axes[0].edges
        n_total = values_flow.sum()
        n_events = []
        for threshold in thresholds:
            idx = np.searchsorted(bin_edges, threshold, side="right") - 1
            idx = max(idx, 0)
            n_events_thres = values[:idx + 1].sum() + underflow
            n_events.append(n_events_thres)

            n_fraction = np.array(n_events) / n_total
       
        plt.plot(thresholds, n_fraction ,color=colors[i],marker="o", label=label )
    if ylabel:
        plt.ylabel(ylabel)
    else:
        plt.ylabel("Fraction of Events")
    plt.xlabel(histogram_name)
    plt.ylim(0, 1.2)
    plt.legend(ncols =2)
    plt.savefig(f"Plots/fraction_{histogram_name}.png")

def plot_fraction_great (signals, histogram_name, channel_name,
                        output_signal_2mu, output_signal_4mu,
                        thresholds, colors, ylabel=None):
    """plot fraction error when we use greater than or equal"""
    plt.figure(figsize =(15, 12))
    for i, s in enumerate(signals):
        label = get_signal_label(s)
        if label[0] == "2":
            histogram = output_signal_2mu["out"][s]["hists"][histogram_name][channel_name, :]
        elif label[0] == "4":
            histogram = output_signal_4mu["out"][s]["hists"][histogram_name][channel_name, :]
        else:
            raise ValueError(f"cannot tell 2mu2e and 4mu apart for sample {s}")
        values = histogram.values()
        values_flow = histogram.values(flow=True)
        overflow = values_flow[-1]
        underflow = values_flow[0]
        bin_edges = histogram.axes[0].edges
        n_total = values_flow.sum()
        n_events = []
        for threshold in thresholds:
            idx = np.searchsorted(bin_edges, threshold, side="right") - 1
            idx = max(idx, 0)
            n_events_thres = values[idx:].sum() + overflow
            n_events.append(n_events_thres)

            n_fraction = np.array(n_events) / n_total
        
        plt.plot(thresholds, n_fraction, label=label, color=colors[i],marker="o" )
        if ylabel:
            plt.ylabel(ylabel)
        else:
            plt.ylabel("Fraction of Events")
        plt.xlabel(histogram_name)
        plt.ylim(0, 1.2)
        plt.legend(ncols =2)
        plt.savefig(f"Plots/fraction_{histogram_name}.png")