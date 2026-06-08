""" Tool to add ntuple locations to sidm/configs/ntuple_locations.yaml.

Developed and tested on FNAL LPC. Will be updated after issue
https://github.com/CoffeaTeam/coffea-casa/issues/374 is resolved. Note that cmsenv or equivalent
is needed to import XRootD.

Usage: python add_ntuples.py -o OUTPUT_CONFIG -n NTUPLE_NAME -c NTUPLE_COMMENT -d NTUPLE_ROOT_DIR
"""

from __future__ import print_function
import argparse
import yaml
import ROOT
from XRootD import client


parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output-cfg", dest="cfg", required=True,
                    help="Path to output config, e.g. '../configs/ntuple_locations.yaml'")
parser.add_argument("-n", "--name", dest="name", required=True,
                    help="Name of group of ntuples, e.g. 'ffntuple_v4'")
parser.add_argument("-c", "--comment", dest="comment", required=True,
                    help=("Comment to describe group of ntuples, e.g. "
                          "'Most recent ntuples from Weinan -- only includes 2mu2e'"))
parser.add_argument("-d", "--directory", dest="directory", required=True,
                    help=("Path to ntuple root directory, e.g. "
                    "'root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/ffNtupleV4/2018/'"))
parser.add_argument("-f", "--first-dir", dest="first_dir", action='store_true',
                    help="Choose first option when encountering unexpected directory structure")
parser.add_argument("-s", "--skim", action='store_true',
                    help="Identify skimmed ntuples")
# fixme: add option to associate multiple subdirectories with one process name
args = parser.parse_args()


def parse_name(name):
    """Parse sample directory name to produce simplified name

    Assumes structure like "SIDM_XXTo2ATo2Mu2E_mXX-100_mA-1p2_ctau-0p096_TuneCP..."
    """

    name = name.removeprefix("CutDecayFalse_")
    name = name.removeprefix("LLPnanoAODv1_")

    process_names = {
        "SIDM_XXTo2ATo2Mu2E_mXX": "2Mu2E_",
        "SIDM_XXTo2ATo4Mu_mXX" : "4Mu_",
        "SIDM_BsTo2DpTo2Mu2e_MBs" : "2Mu2E_",
        "SIDM_BsTo2DpTo4Mu_MBs" : "4Mu_",
        "DYJetsToLL_M" : "DYJetsToLL_M",
        "DYJetsToMuMu_M" : "DYJetsToMuMu_M",
        "QCD_Pt" : "QCD_Pt",
        "TTJets_TuneCP5_13TeV" : "TTJets",
        "TTJets_TuneCP5" : "TTJets",
        "WW_TuneCP5_13TeV" : "WW",
        "WZ_TuneCP5_13TeV" : "WZ",
        "ZZ_TuneCP5_13TeV" : "ZZ",
        "Run2018C" : "DoubleMuon_2018C",
        "Run2018C_15Feb2022_UL2018": "DoubleMuon_2018C",

        "Run2018B_15Feb2022_UL2018_Block1": "DoubleMuon_2018B_0",
        "Run2018B_15Feb2022_UL2018_Block2": "DoubleMuon_2018B_1",
        "Run2018B_15Feb2022_UL2018_Block3": "DoubleMuon_2018B_2",
        "Run2018B_15Feb2022_UL2018_Block4": "DoubleMuon_2018B_3",
        "Run2018B_15Feb2022_UL2018_Block5": "DoubleMuon_2018B_4",
        "Run2018B_15Feb2022_UL2018_Block6": "DoubleMuon_2018B_5",
        "Run2018B_15Feb2022_UL2018_Block7": "DoubleMuon_2018B_6",
        "Run2018B_15Feb2022_UL2018_Block8": "DoubleMuon_2018B_7",
        "Run2018B_15Feb2022_UL2018_Block9": "DoubleMuon_2018B_8",
        "Run2018B_15Feb2022_UL2018_Block10": "DoubleMuon_2018B_9",
        "Run2018B_15Feb2022_UL2018_Block11": "DoubleMuon_2018B_10",
        "Run2018B_15Feb2022_UL2018_Block12": "DoubleMuon_2018B_11",
        "Run2018B_15Feb2022_UL2018_Block13": "DoubleMuon_2018B_12",
        "Run2018B_15Feb2022_UL2018_Block14": "DoubleMuon_2018B_13",

        "Run2018D_15Feb2022_UL2018_Block_0aface09": "DoubleMuon_2018D_0",
        "Run2018D_15Feb2022_UL2018_Block_18f48385": "DoubleMuon_2018D_1",
        "Run2018D_15Feb2022_UL2018_Block_1cf7400c": "DoubleMuon_2018D_2",
        "Run2018D_15Feb2022_UL2018_Block_227e0956": "DoubleMuon_2018D_3",
        "Run2018D_15Feb2022_UL2018_Block_22e16122": "DoubleMuon_2018D_4",
        "Run2018D_15Feb2022_UL2018_Block_2592308d": "DoubleMuon_2018D_5",
        "Run2018D_15Feb2022_UL2018_Block_27d5e0fb": "DoubleMuon_2018D_6",
        "Run2018D_15Feb2022_UL2018_Block_2f365cb7": "DoubleMuon_2018D_7",
        "Run2018D_15Feb2022_UL2018_Block_31b39d2d": "DoubleMuon_2018D_8",
        "Run2018D_15Feb2022_UL2018_Block_327a81b6": "DoubleMuon_2018D_9",
        "Run2018D_15Feb2022_UL2018_Block_3d030fc3": "DoubleMuon_2018D_10",
        "Run2018D_15Feb2022_UL2018_Block_48a90814": "DoubleMuon_2018D_11",
        "Run2018D_15Feb2022_UL2018_Block_51c6e6eb": "DoubleMuon_2018D_12",
        "Run2018D_15Feb2022_UL2018_Block_51c8c04b": "DoubleMuon_2018D_13",
        "Run2018D_15Feb2022_UL2018_Block_5314eda4": "DoubleMuon_2018D_14",
        "Run2018D_15Feb2022_UL2018_Block_57cf89e5": "DoubleMuon_2018D_15",
        "Run2018D_15Feb2022_UL2018_Block_5f23e071": "DoubleMuon_2018D_16",
        "Run2018D_15Feb2022_UL2018_Block_65f115dd": "DoubleMuon_2018D_17",
        "Run2018D_15Feb2022_UL2018_Block_69710812": "DoubleMuon_2018D_18",
        "Run2018D_15Feb2022_UL2018_Block_697a4e54": "DoubleMuon_2018D_19",
        "Run2018D_15Feb2022_UL2018_Block_6b2ace48": "DoubleMuon_2018D_20",
        "Run2018D_15Feb2022_UL2018_Block_7d8d8539": "DoubleMuon_2018D_21",
        "Run2018D_15Feb2022_UL2018_Block_7dae55a6": "DoubleMuon_2018D_22",
        "Run2018D_15Feb2022_UL2018_Block_81881293": "DoubleMuon_2018D_23",
        "Run2018D_15Feb2022_UL2018_Block_82c1a5f9": "DoubleMuon_2018D_24",
        "Run2018D_15Feb2022_UL2018_Block_895cfb04": "DoubleMuon_2018D_25",
        "Run2018D_15Feb2022_UL2018_Block_8a364826": "DoubleMuon_2018D_26",
        "Run2018D_15Feb2022_UL2018_Block_8bcd32b9": "DoubleMuon_2018D_27",
        "Run2018D_15Feb2022_UL2018_Block_8efc6d7f": "DoubleMuon_2018D_28",
        "Run2018D_15Feb2022_UL2018_Block_972e3c94": "DoubleMuon_2018D_29",
        "Run2018D_15Feb2022_UL2018_Block_98f6a4bc": "DoubleMuon_2018D_30",
        "Run2018D_15Feb2022_UL2018_Block_a258a556": "DoubleMuon_2018D_31",
        "Run2018D_15Feb2022_UL2018_Block_a4e18704": "DoubleMuon_2018D_32",
        "Run2018D_15Feb2022_UL2018_Block_b37fdf5c": "DoubleMuon_2018D_33",
        "Run2018D_15Feb2022_UL2018_Block_b3f1aa9e": "DoubleMuon_2018D_34",
        "Run2018D_15Feb2022_UL2018_Block_b8ae72fe": "DoubleMuon_2018D_35",
        "Run2018D_15Feb2022_UL2018_Block_ba86209f": "DoubleMuon_2018D_36",
        "Run2018D_15Feb2022_UL2018_Block_c701125b": "DoubleMuon_2018D_37",
        "Run2018D_15Feb2022_UL2018_Block_c94cb5b6": "DoubleMuon_2018D_38",
        "Run2018D_15Feb2022_UL2018_Block_cd71ffd0": "DoubleMuon_2018D_39",
        "Run2018D_15Feb2022_UL2018_Block_d12b2ef6": "DoubleMuon_2018D_40",
        "Run2018D_15Feb2022_UL2018_Block_d1a63b87": "DoubleMuon_2018D_41",
        "Run2018D_15Feb2022_UL2018_Block_d904c1b4": "DoubleMuon_2018D_42",
        "Run2018D_15Feb2022_UL2018_Block_dae69e31": "DoubleMuon_2018D_43",
        "Run2018D_15Feb2022_UL2018_Block_e71160f8": "DoubleMuon_2018D_44",
        "Run2018D_15Feb2022_UL2018_Block_e766d64d": "DoubleMuon_2018D_45",
        "Run2018D_15Feb2022_UL2018_Block_ee09f7a6": "DoubleMuon_2018D_46",
        "Run2018D_15Feb2022_UL2018_Block_f037f5e1": "DoubleMuon_2018D_47",
        "Run2018D_15Feb2022_UL2018_Block_f18716d8": "DoubleMuon_2018D_48",
        "Run2018D_15Feb2022_UL2018_Block_fec58313": "DoubleMuon_2018D_49",

        "0000": "DoubleMuon_2018A_0",
        "0001": "DoubleMuon_2018A_1",
        "0002": "DoubleMuon_2018A_2",

        # fixme: add backgrounds and data as necessary
    }
    chunks = name.split("-")
    try:
        simplified_name = process_names[chunks[0]] # process name
    except KeyError:
        print("Unrecognized process name. Skipping {}".format(name))
        return None

    # further simplify names as necessary
    if name.startswith("SIDM"):
        #simplified_name += chunks[1].replace("_mA", "GeV_") # bound state mass (weinan notation)
        simplified_name += chunks[1].replace("_MDp", "GeV_") # bound state mass
        simplified_name += chunks[2].replace("_ctau", "GeV_") # dark photon mass
        simplified_name += chunks[3].split("_")[0] + "mm" # dark photon ctau
    elif name.startswith("DYJetsToLL_M"):
        simplified_name += chunks[1].split("_")[0] # mass range
    elif name.startswith("DYJetsToMuMu_M"):
        simplified_name += chunks[1].split("_")[0] # mass range
    elif name.startswith("QCD_Pt"):
        simplified_name += chunks[1].split("_")[0] # pT range

    return simplified_name


def descend(ntuple_path, sample_path, choose_first_dir=False):
    path = ntuple_path + "/" + sample_path
    dir_contents = xrd_client.dirlist(path)[1]
    num_found = dir_contents.size

    if [r for r in dir_contents if r.name.endswith("root")]:
        print("Root files found at this layer. Assuming that these are the ntuples")
        return sample_path

    # Handle emtpy directories
    if num_found == 0:
        print("Found zero objects in {}. Skipping.".format(path))
        return None

    # Allow user to choose directory if more than one is found
    if num_found > 1 and not choose_first_dir:
        print("Unexpected directory structure. Found {} objects in {}".format(num_found, path))
        print("Please type the number of the directory you would like to use. Options are:")
        print("S", "SKIP DIRECTORY")
        for i, x in enumerate(dir_contents):
            print(i, x.name)
        dir_ix = input() # fixme: check input
    else:
        dir_ix = 0

    if dir_ix == "S":
        return None

    return sample_path + "/" + dir_contents.dirlist[int(dir_ix)].name


# Set up xrd client
redirector = args.directory.split("//store")[0]
xrd_client = client.FileSystem(redirector)

coffea_casa_dir = args.directory.replace("cmseos.fnal.gov", "xcache")
output = {
    args.name: {
        "path": coffea_casa_dir,
        "samples": {},
    }
}
ntuple_path = args.directory.split(redirector)[1]
samples = xrd_client.dirlist(ntuple_path)[1]

# Traverse ntuple directory and construct output dictionary
# Assumes same structure as root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/ffNtupleV4/2018/
for sample in samples:
    simple_name = parse_name(sample.name)
    print(f"{sample.name} --> {simple_name}")
    if simple_name is None:
        continue
    output[args.name]["samples"][simple_name] = {}
    sample_path = sample.name

    # Descend one layer, expecting to find a single directory
    try:
        for _ in range(2):
            sample_path = descend(ntuple_path, sample_path, args.first_dir)
            if sample_path is None:
                raise StopIteration()
    except StopIteration:
        continue

    # If traversal was successful, add path and files to output dictionary
    try:
        files = [f.name for f in xrd_client.dirlist(ntuple_path + sample_path)[1]]
        # Handle cases with additional directory layer
        if len(files) == 1 and "0000" in files:
            sample_path += "/0000"
            files = [f.name for f in xrd_client.dirlist(ntuple_path + sample_path)[1]]
    except TypeError:
        print("Unexpected directory structure. Skipping {}".format(sample_path))
    # remove non-root files
    files = [f for f in files if f.endswith("root")]
    output[args.name]["samples"][simple_name]["path"] = sample_path + "/"
    output[args.name]["samples"][simple_name]["files"] = files

    # get skim factor
    skim_factor = 1.0
    if args.skim:
        skimmed_evts = 0
        original_evts = 0
        # check first 10 files for valid OriginalEventIndex and sufficient stats
        for file in files[:10000]:
            file_path = f"{redirector}//{ntuple_path}{sample_path}/{file}"
            try:
                in_file = ROOT.TFile.Open(file_path)
            except Exception as e:
                print(e)
                print(f"failed reading {file_path}")
                print("skipping")
                continue
            if not in_file or in_file.IsZombie():
                print(f"failed reading {file_path}")
                print("skipping")
                continue

            # tree = in_file.Events
            tree = in_file.Get("Events")
            if tree is None:
                print(f"No Events tree in {file_path}")
                print("skipping")
                in_file.Close()
                continue

            for entry in tree:
                file_skimmed_evts = tree.GetEntries()
                file_original_evts = entry.OriginalEventIndex
                break
            in_file.Close()

            # Go to next file if OriginalEventIndex is invalid
            if file_original_evts == 0:
                print("invalid OriginalEventIndex, going to next file")
                #continue
            else:
                skimmed_evts += file_skimmed_evts
                original_evts += file_original_evts
                print(skimmed_evts, original_evts, file_skimmed_evts/file_original_evts, skimmed_evts/original_evts)
            if skimmed_evts > 100000:
                break
        try:
            skim_factor = skimmed_evts / original_evts
        except ZeroDivisionError:
            print("No valid OriginalEventIndex found, setting skim_factor to 1.0")
    print(f"Setting skim factor to {skim_factor}")
    output[args.name]["samples"][simple_name]["skim_factor"] = skim_factor

    if "DoubleMuon_2018" in simple_name:
        output[args.name]["samples"][simple_name]["is_data"] = True
        output[args.name]["samples"][simple_name]["year"] = "2018"


# Avoid yaml references, a la stackoverflow.com/questions/13518819
yaml.Dumper.ignore_aliases = lambda *args: True

with open(args.cfg, 'a') as out_file:
    out_file.write("\n\n# " + args.comment + "\n")
    yaml.dump(output, out_file, default_flow_style=False)
    out_file.write("\n")

