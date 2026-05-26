# SIDM
SIDM analysis at coffea-casa.  
Inspired by github.com/phylsix/Firefighter and/or github.com/phylsix/FireROOT

## Getting started
- Fork this repository ([here's a nice guide to follow](https://gist.github.com/Chaser324/ce0505fbed06b947d962))
- Log in to coffea.casa as described [here](https://coffea-casa.readthedocs.io/en/latest/cc_user.html#cms-authz-authentication-instance)
- Clone your fork of this repository using the [coffea-casa git interface](https://coffea-casa.readthedocs.io/en/latest/cc_user.html#using-git). Note the following:
  - Use the https link, not the ssh one (e.g. `https://github.com/btcardwell/SIDM.git`, not `git@github.com:btcardwell/SIDM.git`)
  - You will likely need to generate a [github personal access token](https://docs.github.com/en/enterprise-server@3.4/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) to log in to your github account on coffea-casa
- Navigate to the newly created SIDM directory using the coffea-casa file browser
- From the SIDM directory, run setup.sh to pip install the sidm package
- You should be good to go! If you want to test that your environment is set up correctly, try running any of the existing notebooks in SIDM/studies and comparing your output with the output shown [here](https://github.com/btcardwell/SIDM/tree/main/sidm/studies)

## Getting started on LPC

This is **one** way to get a working interactive environment on Fermilab LPC (`cmslpc-el9.fnal.gov`). Other paths exist — e.g. the [Elastic Analysis Facility (EAF)](https://analytics-hub.fnal.gov), JupyterHub on LPC, or plain shell sessions — and the venv-setup steps below (1–4) transfer directly to those. Only step 5 (how to actually open a notebook) is specific to running VS Code on a laptop over SSH; on EAF or JupyterHub you skip it and pick the registered kernel from the web UI instead.

This walkthrough covers the **interactive / notebook** environment. For running large batch jobs on LPC Condor, see [condor/README.md](condor/README.md) once steps 1–4 are done.

### 1. Clone the repo

SSH to LPC and clone into your work area:

    ssh cmslpc-el9.fnal.gov
    cd /uscms_data/d3/$USER
    git clone <your-fork-of-SIDM>
    cd SIDM

### 2. Build the venv on LCG_107 Python 3.11

The repo's `requirements.txt` pins versions (`dask`, `awkward`, …) that need Python ≥ 3.10. LPC's system Python is 3.9, so we use Python 3.11 from cvmfs:

    source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh
    unset PYTHONPATH
    python -m venv sidm_venv
    source sidm_venv/bin/activate
    python -m pip install --upgrade pip setuptools wheel

After creation the venv's `bin/python` is a symlink to the cvmfs binary, so you don't need to re-source LCG_107 to run it later.

### 3. Install requirements (skipping the xrootd PyPI build)

The `xrootd` PyPI wheel tries to compile the xrootd C library from source and fails on LPC. We get the same functionality by symlinking LCG's prebuilt `XRootD` and `pyxrootd` packages into the venv:

    grep -v "^xrootd" requirements.txt | python -m pip install -r /dev/stdin
    python -m pip install "distributed==2025.3.0"   # required by sidm/tools/scaleout.py
    python -m pip install -e .
    python -m pip install jupyter ipykernel pyarrow

    cd sidm_venv/lib/python3.11/site-packages
    ln -sfn /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/lib/python3.11/site-packages/XRootD   XRootD
    ln -sfn /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/lib/python3.11/site-packages/pyxrootd pyxrootd
    cd -

Smoke-test that uproot can read a remote ROOT file through xrootd:

    sidm_venv/bin/python -c "import uproot; print(uproot.open('root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/ULSignalSamples/2018_v10/BsTo2DpTo2Mu2e/CutDecayFalse_SIDM_BsTo2DpTo2Mu2e_MBs-500_MDp-1p2_ctau-1p9_v3/LLPnanoAODv2/CutDecayFalse_SIDM_BsTo2DpTo2Mu2e_MBs-500_MDp-1p2_ctau-1p9_v3_part-0.root')['Events'].num_entries)"

### 4. Register the Jupyter kernel

    /uscms_data/d3/$USER/SIDM/sidm_venv/bin/python -m ipykernel install --user \
        --name sidm_venv \
        --display-name "SIDM (LCG_107 Py3.11)"

This drops a kernelspec under `~/.local/share/jupyter/kernels/sidm_venv/` on LPC. Any Jupyter Server running on LPC — including EAF or JupyterHub — will offer this kernel.

### 5. Open notebooks (one approach: VS Code over SSH)

VS Code on your laptop can't directly see LPC-side kernels through an SSHFS mount, since the kernel binary is x86_64 Linux. One pattern that works: start a Jupyter Server on LPC, port-forward 8888 to your laptop, point VS Code at `http://localhost:8888/` as a remote server.

In a terminal on your laptop:

    ssh -L 8888:localhost:8888 cmslpc-el9.fnal.gov \
        "cd /uscms_data/d3/$USER/SIDM && source sidm_venv/bin/activate && \
         jupyter server --no-browser --port=8888 --ip=127.0.0.1"

Copy the URL it prints (`http://localhost:8888/?token=…`). In VS Code: open any notebook → kernel selector (top-right) → **Select Another Kernel…** → **Existing Jupyter Server…** → paste URL → pick **SIDM (LCG_107 Py3.11)**.

If you're on EAF, JupyterHub, or running Jupyter directly on LPC with your own tunneling, skip this step entirely and just pick the kernel from your usual UI.

### 6. Reading remote ROOT files

The location YAMLs (`sidm/configs/ntuples/*.yaml`) use `root://xcache//…` URLs. `xcache` is the coffea-casa cache and does not resolve on LPC. Pass `replace_xcache=True` to `utilities.make_fileset` to substitute the FNAL EOS redirector at fileset-build time:

    fileset = utilities.make_fileset(
        samples, "llpNanoAOD_v2", location_cfg="signal_2mu2e_v10.yaml",
        replace_xcache=True,
    )

Default is `False` so existing coffea-casa notebooks remain unchanged.

For batch processing on LPC Condor, see [condor/README.md](condor/README.md).

## Code structure
All the interesting code in this repository is in `SIDM/sidm/`, which is organized into the following subdirectories:
- `tools` contains the classes and methods that form the backbone of SIDM. In particular, `sidm_processor.py` defines how events are analyzed. We use the NanoAODSchema for the data in our files, as defined within Coffea [here](https://coffea-hep.readthedocs.io/en/latest/api/coffea.nanoevents.NanoAODSchema.html)
- `definitions` contains files that define the current set of histograms one can make, the cuts one can apply, and some of the objects one can use when making histograms and applying cuts.
- `configs` contains yaml configuration files that define how histograms are grouped together into collections and cuts are grouped together into selections. When running `sidm_processor`, one provides names of selections and names of histogram collections to choose which cuts to apply and which histograms to make.
- `test_notebooks` contains notebooks to test new classes or functionalities as they are added. These notebooks can also serve as a form  of unit test: if you edit some code in a way that you think shouldn't affect the behavior, you can run these notebooks to confirm the output is unchanged.
- `studies` is where the physics happens. The notebooks in this directory are meant to serve effectively as pages in a lab notebook. The intention is to create a new notebook for each unique physics study and to include markdown comments to describe the intentions and observations of the person performing the study. These notebooks can also serve as unit tests in the same way as those in `test_notebooks`.

## Analysis how-tos

### General workflow
I suggest the following workflow for performing a physics study:
1. Create a new branch with a descriptive name (e.g. `leptonIsolation`)
2. Create a new notebook in `studies` with a descriptive name (e.g. `study_lepton_isolation.ipynb`)
3. Following the examples of the existing notebooks in `studies`, use `sidm_processor` to apply cuts and make histograms starting from an Firefighter ntuple of your choosing. Make sure to describe your reasoning and observations in text as you go.
4. If you find you need to define new selections, new cuts, or new histograms, follow the guides below.
5. Commit your changes as you go and submit a Pull Request once you have a reasonable standalone study or have added new selections, cuts, histograms, objects, classes, or features.

### How to define a new histogram and add it to a collection
1. Add an entry to the `hist_defs` dictionary inside [hists.py](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/definitions/hists.py). One can potentially do this by mimicking the structure of the existing histograms, but here are some details for those who are interested:
    - Each histogram in `hist_defs` requires a name and Histogram object, which is created using the [Histogram constructor](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/tools/histogram.py#L14-L18).
    - The only required argument when creating a Hist is a list of Axis objects, which are created using the [Axis constructor](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/tools/histogram.py#L47-L50).
    - When creating an Axis, one must provide a [Hist axis](https://hist.readthedocs.io/en/latest/user-guide/axes.html) and a fill function, which is a lambda expression that defines the quantity used to fill the histogram axis. Note that the fill function will usually take the object from which the quantity is derived as an argument (this is the magic bit that allows us to define how histograms will be filled before running the analysis).
2. After defining your new histogram in `hists.py`, add the name of the histogram to one of the histogram collections in [hist_collections.yaml](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/configs/hist_collections.yaml). Alternatively, you can also create an entirely new histogram collection that includes the new cut. Note that the new histogram will automatically be included in any collections that include the collection to which it was added (e.g. any histograms added to [electron_base](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/configs/hist_collections.yaml#L12-L15) will automatically be included in [base](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/configs/hist_collections.yaml#L77-L88))
3. That's it! Running `sidm_processor` while specifying the proper histogram collection will now produce the new histogram.

### How to define a new cut and add it to a selection
1. Add a new entry to either the [obj_cut_defs](https://github.com/btcardwell/SIDM/blob/4e6685669067429e8492d4dcfc87f463c86b96d7/analysis/definitions/cuts.py#L10-L25) or [evt_cut_defs](https://github.com/btcardwell/SIDM/blob/4e6685669067429e8492d4dcfc87f463c86b96d7/analysis/definitions/cuts.py#L27-L36) dictionary in [cuts.py](https://github.com/btcardwell/SIDM/blob/4e6685669067429e8492d4dcfc87f463c86b96d7/analysis/definitions/cuts.py). Each cut must have a name and a lambda expression that defines the cut. The difference between object-level and event-level cuts is as follows:
    - Object-level cuts slim object collections, i.e. remove objects from a given collection without accepting or rejecting the whole event. For example, one might apply an object-level cut to remove low-momentum electrons from events.
    - Event-level cuts accept or reject whole events and are applied after object-level cuts. For example, one could apply an event-level cut that rejects events without at least two electrons that pass the above-mentioned object-level cut.
2. After defining your new cut in `cuts.py`, add it's name to an existing selection in [selections.yaml](https://github.com/btcardwell/SIDM/blob/4e6685669067429e8492d4dcfc87f463c86b96d7/analysis/configs/selections.yaml). Alternatively, you can also create an entirely new selection. Note that object-level and event-level cuts must be listed separately within a selection, and object-level cuts are further organized by object type. As with histogram collections, any selections that include the selection to which you added your cut will also include your cut.
3. That's it! Running `sidm_processor` while specifying the proper selection will now apply your new cut.

## Miscellaneous how-tos

### How to update requirements.txt
```
cd SIDM/
pip install pipreqs
pipreqs . --force # overwrites current requirements.txt
```

### How to use DASK
1. In your jupyter notebook make sure you import `scaleout` from `sidm.tools`
2. Modify the dependencies in  `scaleout.py` from `sidm/tools` so that DASK points to your branch from your fork. i.e.
```
dependencies = [
        "git+https://github.com/YOURUSERNAME/SIDM.git@YOURBRANCH",
    ]
```
3. In a cell of your Jupyter Notebook, instantiate the DASK Client:
```
client = scaleout.make_dask_client("tls://localhost:8786")
```
4. In your `processor.Runner()` function make sure `executor=processor.DaskExecutor(client=client)`