# SIDM Condor Workflow on CMSLPC

This directory contains the scripts needed to run the SIDM Coffea processor on CMSLPC Condor.

The workflow is:

    samples + YAML locations
            ↓
    make_job_args.py
            ↓
    many Condor jobs, each processing a ROOT-file chunk
            ↓
    chunk .coffea outputs on EOS
            ↓
    merge_coffea_chunks_eos.py
            ↓
    one merged .coffea file per sample on EOS

This setup replaces the Dask workflow with Condor-scale parallelism.

---

## 1. Repository assumptions

From the repository parent:

    SIDM/
      sidm/
        tools/
          sidm_processor.py
          utilities.py
          llpnanoaodschema.py
      condor/
        make_job_args.py
        run_sidm_chunk.py
        run_job.sh
        submit.sub
        requirements.txt
        signal_samples.txt
        background_samples.txt

Most commands assume you are starting from:

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM

---

## 2. Setup CMSSW environment

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM

    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=el9_amd64_gcc12

This workflow was tested with:

    CMSSW_14_1_0_pre4
    SCRAM_ARCH=el9_amd64_gcc12
    Python 3.9.14

---

## 3. Condor directory setup

    mkdir -p condor/filelists
    mkdir -p condor/logs
    mkdir -p condor/output

---

## 4. Condor requirements file

The Condor jobs create a temporary Python virtual environment inside the Condor scratch directory and install the packages listed in:

    condor/requirements.txt

Recommended contents:

    awkward==2.8.2
    coffea==2025.5.0rc2
    fastjet==3.4.3.1
    hist==2.8.0
    numpy==1.26.4
    PyYAML==6.0.2
    setuptools==80.1.0
    tabulate==0.9.0
    vector==1.6.2
    uproot==5.6.2

---

## 5. Sample lists

Sample-list files should contain one sample name per line, with no quotes and no commas.

Example background file:

    condor/background_samples.txt

Example contents:

    TTJets
    QCD_Pt170To300
    QCD_Pt300To470
    QCD_Pt470To600
    QCD_Pt600To800
    DYJetsToMuMu_M10to50
    DYJetsToMuMu_M50
    QCD_Pt1000
    QCD_Pt120To170
    QCD_Pt15To20
    QCD_Pt20To30
    QCD_Pt30To50
    QCD_Pt50To80
    QCD_Pt80To120
    QCD_Pt800To1000

Example signal file:

    condor/signal_samples.txt

Example contents:

    2Mu2E_200GeV_0p25GeV_0p1mm
    2Mu2E_200GeV_0p25GeV_1p0mm
    2Mu2E_200GeV_0p25GeV_5p0mm

---

## 6. Generate Condor job arguments

The script:

    condor/make_job_args.py

uses:

    utilities.make_fileset(samples, dataset, location_cfg=...)

to build the list of ROOT files, split them into chunks, and write:

    condor/job_args.txt
    condor/filelists/*.txt

### Backgrounds

For backgrounds, use chunks of around 10 ROOT files per job:

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM

    rm -f condor/filelists/*.txt condor/job_args.txt

    python3 condor/make_job_args.py \
      --samples-file condor/background_samples.txt \
      --dataset llpNanoAOD_v2 \
      --location-cfg backgrounds.yaml \
      --files-per-job 10 \
      --replace-xcache

This creates many jobs per background sample. For example, if a sample has 3000 ROOT files and `--files-per-job 10`, that sample becomes about 300 Condor jobs.

### Signals

For Dask-like parallelism over signal samples, use chunks of around 5 ROOT files per job:

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM

    rm -f condor/filelists/*.txt condor/job_args.txt

    python3 condor/make_job_args.py \
      --samples-file condor/signal_samples.txt \
      --dataset llpNanoAOD_v2 \
      --location-cfg signal_2mu2e_v10.yaml \
      --files-per-job 5 \
      --replace-xcache

This creates multiple Condor jobs per signal point. The chunks are merged afterward.

---

## 7. Important path fixes

### Fix the `condor/filelists/` double-prefix issue

Because we submit from inside the `condor/` directory, the filelist paths inside `job_args.txt` must look like:

    filelists/QCD_Pt800To1000_407.txt

not:

    condor/filelists/QCD_Pt800To1000_407.txt

If needed, fix with:

    sed -i 's#condor/filelists/#filelists/#g' condor/job_args.txt

Check:

    grep -n "condor/filelists" condor/job_args.txt

This should print nothing.

### Fix `root://xcache//` paths

Some YAML configs may produce paths like:

    root://xcache//store/...

These may fail on Condor workers. Replace them with:

    root://cmseos.fnal.gov//store/...

Run:

    sed -i 's#root://xcache//#root://cmseos.fnal.gov//#g' condor/filelists/*.txt

Check:

    grep -R "root://xcache" condor/filelists | head

This should print nothing.

---

## 8. Package the SIDM code

Whenever files under `sidm/` or `condor/run_sidm_chunk.py` change, remake the tarball:

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM

    tar \
      --exclude="*.root" \
      --exclude="__pycache__" \
      --exclude=".git" \
      --exclude="condor/logs" \
      --exclude="condor/output" \
      --exclude="condor/filelists" \
      --exclude="condor/filelists_retry" \
      --exclude="signal_chunks" \
      --exclude="signal_merged" \
      --exclude="background_chunks" \
      --exclude="background_merged" \
      --exclude="sidm_check_venv" \
      --exclude="py39_packages" \
      -czf condor/sidm_code.tar.gz \
      sidm condor/run_sidm_chunk.py

Check:

    ls -lh condor/sidm_code.tar.gz

You do not need to re-tar if you only changed:

    condor/requirements.txt
    condor/submit.sub
    condor/run_job.sh
    condor/job_args.txt
    condor/filelists/*.txt

Those are transferred separately.

---

## 9. Create EOS output directories

Use separate directories for signal chunks, background chunks, and merged outputs.

    xrdfs root://cmseos.fnal.gov mkdir -p /store/user/$USER/sidm_condor/SignalChunks_v1
    xrdfs root://cmseos.fnal.gov mkdir -p /store/user/$USER/sidm_condor/BackgroundChunks_v1

    xrdfs root://cmseos.fnal.gov mkdir -p /store/user/$USER/sidm_condor/SignalMerged_v1
    xrdfs root://cmseos.fnal.gov mkdir -p /store/user/$USER/sidm_condor/BackgroundMerged_v1

---

## 10. Condor submit file

File:

    condor/submit.sub

Example:

    universe = vanilla

    executable = run_job.sh

    arguments = $(sample) $(chunk) $(filelist) root://cmseos.fnal.gov//store/user/YOUR_USERNAME/sidm_condor/BackgroundChunks_v1

    should_transfer_files = YES
    when_to_transfer_output = ON_EXIT

    transfer_input_files = sidm_code.tar.gz, requirements.txt, $(filelist)
    transfer_output_files = ""

    output = logs/$(sample)_$(chunk)_$(Cluster)_$(Process).out
    error  = logs/$(sample)_$(chunk)_$(Cluster)_$(Process).err
    log    = logs/$(sample)_$(chunk)_$(Cluster)_$(Process).log

    request_cpus = 1
    request_memory = 4GB
    request_disk = 4GB

    x509userproxy = $ENV(X509_USER_PROXY)

    queue sample, chunk, filelist from job_args.txt

For signals, change the output path to:

    root://cmseos.fnal.gov//store/user/YOUR_USERNAME/sidm_condor/SignalChunks_v1

For backgrounds, use:

    root://cmseos.fnal.gov//store/user/YOUR_USERNAME/sidm_condor/BackgroundChunks_v1

Replace `YOUR_USERNAME` with your LPC username.

---

## 11. Proxy setup

Before submitting jobs:

    voms-proxy-init --valid 192:00 -voms cms
    export X509_USER_PROXY=$(voms-proxy-info -path)

    voms-proxy-info -timeleft

---

## 12. Submit jobs

Submit from inside the `condor/` directory:

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM/condor

    export X509_USER_PROXY=$(voms-proxy-info -path)

    condor_submit submit.sub

Example output:

    Attempting to submit jobs to lpcschedd5.fnal.gov
    5787 job(s) submitted to cluster 59092403.

Record:

    schedd:  lpcschedd5.fnal.gov
    cluster: 59092403

---

## 13. Check job status

Use the schedd and cluster ID printed by `condor_submit`.

    condor_q -name lpcschedd5.fnal.gov 59092403

Summary:

    condor_q -name lpcschedd5.fnal.gov 59092403 -totals

Held jobs:

    condor_q -name lpcschedd5.fnal.gov 59092403 -hold

Watch live:

    watch -n 10 'condor_q -name lpcschedd5.fnal.gov 59092403 -totals'

If jobs disappear from `condor_q`, they are no longer active. Check logs and EOS outputs.

---

## 14. Check logs

From the `condor/` directory:

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM/condor

    grep -R "return value 0" logs/*59092403*.log | wc -l
    grep -R "return value 1" logs/*59092403*.log | wc -l
    grep -R "return value 2" logs/*59092403*.log | wc -l

    grep -R "Traceback\|ERROR\|Exception\|No such file" logs/*59092403*.err | tail -100

A successful job should have:

    Normal termination (return value 0)

---

## 15. Compare expected outputs to actual EOS outputs

From the repo parent:

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM

    awk '{print $1"_"$2".coffea"}' condor/job_args.txt | sort > expected_outputs.txt

    xrdfs root://cmseos.fnal.gov ls /store/user/$USER/sidm_condor/BackgroundChunks_v1 \
      | xargs -n1 basename \
      | sort > actual_outputs.txt

    comm -23 expected_outputs.txt actual_outputs.txt > missing_outputs.txt

    echo "Expected:"
    wc -l expected_outputs.txt

    echo "Actual:"
    wc -l actual_outputs.txt

    echo "Missing:"
    wc -l missing_outputs.txt

    head missing_outputs.txt

For signals, replace `BackgroundChunks_v1` with `SignalChunks_v1`.

---

## 16. Retry failed jobs from missing ROOT files

A common failure is:

    OSError: Failed to open file:
    No such file or directory

This means one ROOT file in a chunk is missing on EOS. One missing ROOT file causes the entire chunk to fail.

### Get failed ProcIds

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM

    condor_history -name lpcschedd5.fnal.gov \
      -constraint 'ClusterId == 59092403 && ExitCode != 0' \
      -af ProcId ExitCode > failed_procs_59092403.txt

    wc -l failed_procs_59092403.txt
    head failed_procs_59092403.txt

### Build filtered retry filelists

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM

    mkdir -p condor/filelists_retry
    rm -f condor/failed_retry_job_args_59092403.txt

    while read PROC EXITCODE; do
        LINE=$(awk -v p="$PROC" 'NR==p+1 {print}' condor/job_args.txt)

        SAMPLE=$(echo "$LINE" | awk '{print $1}')
        CHUNK=$(echo "$LINE" | awk '{print $2}')
        FILELIST=$(echo "$LINE" | awk '{print $3}')

        OLD_FILELIST="condor/${FILELIST}"
        NEW_FILELIST="condor/filelists_retry/${SAMPLE}_${CHUNK}.txt"

        echo "Checking $SAMPLE chunk $CHUNK"

        > "$NEW_FILELIST"

        while read ROOTFILE; do
            EOSPATH=$(echo "$ROOTFILE" | sed 's#root://cmseos.fnal.gov//#/#' | sed 's#^//store#/store#')

            if xrdfs root://cmseos.fnal.gov stat "$EOSPATH" >/dev/null 2>&1; then
                echo "$ROOTFILE" >> "$NEW_FILELIST"
            else
                echo "  MISSING: $ROOTFILE"
            fi
        done < "$OLD_FILELIST"

        NGOOD=$(wc -l < "$NEW_FILELIST")

        if [ "$NGOOD" -gt 0 ]; then
            echo "$SAMPLE $CHUNK filelists_retry/${SAMPLE}_${CHUNK}.txt" >> condor/failed_retry_job_args_59092403.txt
        else
            echo "  WARNING: no good files left for $SAMPLE chunk $CHUNK"
        fi

    done < failed_procs_59092403.txt

### Submit retry jobs

Create a retry submit file:

    cp condor/submit.sub condor/submit_retry_59092403.sub

Edit the final queue line:

    queue sample, chunk, filelist from failed_retry_job_args_59092403.txt

Submit:

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM/condor

    export X509_USER_PROXY=$(voms-proxy-info -path)

    condor_submit submit_retry_59092403.sub

---

## 17. Merge `.coffea` chunks from EOS

Use the EOS-aware merge script:

    merge_coffea_chunks_eos.py

This script:

    input EOS chunk directory
            ↓
    copies chunks to temporary local workdir
            ↓
    groups chunks by sample
            ↓
    uses processor.accumulate(outputs)
            ↓
    writes merged .coffea locally
            ↓
    copies merged .coffea to output EOS directory
            ↓
    optionally deletes temporary local files

### Merge backgrounds

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM
    source sidm_check_venv/bin/activate

    python merge_coffea_chunks_eos.py \
      --input-eos-dir /store/user/$USER/sidm_condor/BackgroundChunks_v1 \
      --output-eos-dir /store/user/$USER/sidm_condor/BackgroundMerged_v1 \
      --workdir merge_tmp_background

### Merge signals

    cd /uscms_data/d3/$USER/CMSSW_14_1_0_pre4/src/SIDM
    source sidm_check_venv/bin/activate

    python merge_coffea_chunks_eos.py \
      --input-eos-dir /store/user/$USER/sidm_condor/SignalChunks_v1 \
      --output-eos-dir /store/user/$USER/sidm_condor/SignalMerged_v1 \
      --workdir merge_tmp_signal

### Check merged outputs

    xrdfs root://cmseos.fnal.gov ls /store/user/$USER/sidm_condor/BackgroundMerged_v1
    xrdfs root://cmseos.fnal.gov ls /store/user/$USER/sidm_condor/SignalMerged_v1

    xrdfs root://cmseos.fnal.gov ls /store/user/$USER/sidm_condor/BackgroundMerged_v1 | wc -l
    xrdfs root://cmseos.fnal.gov ls /store/user/$USER/sidm_condor/SignalMerged_v1 | wc -l

Expected:

    BackgroundMerged_v1: one .coffea per background sample
    SignalMerged_v1: one .coffea per signal sample

---

## 18. Access merged files in Coffea-Casa or CMSLPC notebooks

The merged EOS paths are:

    root://cmseos.fnal.gov//store/user/$USER/sidm_condor/SignalMerged_v1/
    root://cmseos.fnal.gov//store/user/$USER/sidm_condor/BackgroundMerged_v1/

For `.coffea` files, it is usually safer to copy locally before loading with `coffea.util.load()`.

Example from CMSLPC:

    mkdir -p signal_merged background_merged

    for f in $(xrdfs root://cmseos.fnal.gov ls /store/user/$USER/sidm_condor/SignalMerged_v1 | grep ".coffea"); do
        base=$(basename "$f")
        xrdcp -f "root://cmseos.fnal.gov/${f}" "signal_merged/${base}"
    done

    for f in $(xrdfs root://cmseos.fnal.gov ls /store/user/$USER/sidm_condor/BackgroundMerged_v1 | grep ".coffea"); do
        base=$(basename "$f")
        xrdcp -f "root://cmseos.fnal.gov/${f}" "background_merged/${base}"
    done

Then in Python:

    import glob
    from coffea.util import load

    signal_files = sorted(glob.glob("signal_merged/*.coffea"))
    background_files = sorted(glob.glob("background_merged/*.coffea"))

    sig = load(signal_files[0])
    bkg = load(background_files[0])

    print(sig.keys())
    print(bkg.keys())

---

## 19. Files that should be committed

Commit reusable scripts and templates:

    git add condor/README.md
    git add condor/make_job_args.py
    git add condor/run_sidm_chunk.py
    git add condor/run_job.sh
    git add condor/submit.sub
    git add condor/requirements.txt
    git add condor/signal_samples.txt
    git add condor/background_samples.txt
    git add merge_coffea_chunks_eos.py

Also commit relevant processor changes:

    git add sidm/tools/sidm_processor.py

Do not commit generated files or credentials.

---

## 20. Files that should not be committed

Do not commit:

    *.coffea
    *.root
    *.parquet
    *.json
    *.pem
    x509up_*
    sidm_check_venv/
    py39_packages/
    condor/logs/
    condor/filelists/
    condor/filelists_retry/
    condor/job_args.txt
    condor/sidm_code.tar.gz
    actual_outputs.txt
    expected_outputs.txt
    missing_outputs.txt
    failed_procs_*.txt
    signal_chunks/
    signal_merged/
    background_chunks/
    background_merged/

Add these to `.gitignore`.

---

## 21. Common errors

### `condor/condor/filelists/... No such file or directory`

Cause: `job_args.txt` contains `condor/filelists/...` while submitting from inside `condor/`.

Fix:

    sed -i 's#condor/filelists/#filelists/#g' condor/job_args.txt

### `root://xcache//... Invalid address`

Cause: worker cannot resolve `root://xcache//`.

Fix:

    sed -i 's#root://xcache//#root://cmseos.fnal.gov//#g' condor/filelists/*.txt

### `No module named coffea`

Cause: Condor worker Python does not have Coffea.

Fix: ensure `run_job.sh` creates the venv and installs:

    python -m pip install --no-cache-dir -r requirements.txt

### `KeyError: skim_factor`

Cause: processor expects `events.metadata["skim_factor"]`.

Fix: either add metadata in `run_sidm_chunk.py`:

    "metadata": {
        "dataset": args.sample,
        "skim_factor": 1.0,
    }

or make the processor use a default:

    skim_factor = events.metadata.get("skim_factor", 1.0)

### Missing ROOT file

Cause: file listed in YAML/filelist is not present on EOS.

Fix: use the retry workflow to remove missing ROOT files from failed chunks and resubmit.

