# SIDM Condor Workflow on CMSLPC

This directory contains the scripts needed to run the SIDM Coffea processor on CMSLPC Condor.

The workflow is:

```
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
```

This setup replaces the Dask workflow with Condor-scale parallelism.

A runnable walkthrough of the full pipeline (with every shell command, file description, and a sanity plot at the end) lives at [`sidm/test_notebooks/lpc_condor_example.ipynb`](../sidm/test_notebooks/lpc_condor_example.ipynb). Read it first if you are new — this README is the detailed reference.

---

## 1. Repository assumptions

From the repository parent:

```
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
```

Most commands assume you are starting from:

```bash
cd /uscms_data/d3/$USER/SIDM
```

This is the same standalone checkout built in the [top-level README](../README.md) — a directory in your data area with `sidm_venv/` inside it. The venv records this absolute path, so keep the checkout where its venv was built; if you relocate it, rebuild the venv (top-README steps 2–4) before submitting.

---

## 2. Activate your sidm_venv

The submitter side (where you run `make_job_args.py` and `condor_submit`) just needs the LCG_107 `sidm_venv` from the top-level README — no CMSSW required, because `make_job_args.py` only imports `coffea` / `awkward` / `PyYAML` which the venv already has.

```bash
cd /uscms_data/d3/$USER/SIDM
source sidm_venv/bin/activate
```

The Condor *workers* (where the actual processor runs) do still need CMSSW; that is set up automatically inside `run_job.sh` and you do not need to do anything for it on the submitter side.

This workflow was tested with:

```
Submitter: LCG_107 Python 3.11.9, coffea==2025.5.0rc2
Worker:    CMSSW_14_1_0_pre4, SCRAM_ARCH=el9_amd64_gcc12
```

---

## 3. Condor directory setup

```bash
mkdir -p condor/filelists
mkdir -p condor/logs
mkdir -p condor/output
```

---

## 4. Condor requirements file

The Condor jobs create a temporary Python virtual environment inside the Condor scratch directory and install the packages listed in:

```
condor/requirements.txt
```

Recommended contents:

```
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
```

---

## 5. Sample lists

Sample-list files should contain one sample name per line, with no quotes and no commas.

Example background file:

```
condor/background_samples.txt
```

Example contents:

```
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
```

Example signal file:

```
condor/signal_samples.txt
```

Example contents:

```
2Mu2E_200GeV_0p25GeV_0p1mm
2Mu2E_200GeV_0p25GeV_1p0mm
2Mu2E_200GeV_0p25GeV_5p0mm
```

---

## 6. Generate Condor job arguments

The script:

```
condor/make_job_args.py
```

uses:

```python
utilities.make_fileset(samples, dataset, location_cfg=...)
```

to build the list of ROOT files, split them into chunks, and write:

```
condor/job_args.txt
condor/filelists/*.txt
```

### Backgrounds

For backgrounds, use chunks of around 10 ROOT files per job:

```bash
cd /uscms_data/d3/$USER/SIDM

rm -f condor/filelists/*.txt condor/job_args.txt

python3 condor/make_job_args.py \
  --samples-file condor/background_samples.txt \
  --dataset llpNanoAOD_v2 \
  --location-cfg backgrounds.yaml \
  --files-per-job 10 \
  --replace-xcache
```

This creates many jobs per background sample. For example, if a sample has 3000 ROOT files and `--files-per-job 10`, that sample becomes about 300 Condor jobs.

### Signals

For Dask-like parallelism over signal samples, use chunks of around 5 ROOT files per job:

```bash
cd /uscms_data/d3/$USER/SIDM

rm -f condor/filelists/*.txt condor/job_args.txt

python3 condor/make_job_args.py \
  --samples-file condor/signal_samples.txt \
  --dataset llpNanoAOD_v2 \
  --location-cfg signal_2mu2e_v10.yaml \
  --files-per-job 5 \
  --replace-xcache
```

This creates multiple Condor jobs per signal point. The chunks are merged afterward.

---

## 7. Path fixes (no longer needed)

Earlier revisions of this guide had two manual `sed` workarounds in this section: one to strip a `condor/filelists/` double-prefix from `job_args.txt`, and one to rewrite `root://xcache//` to `root://cmseos.fnal.gov//` in every filelist. Both are now fixed at source:

- `make_job_args.py` writes filelist paths in `job_args.txt` relative to `condor/`, so the double-prefix never appears.
- Passing `--replace-xcache` to `make_job_args.py` (see §6) rewrites xcache URLs at fileset-build time.

Skip this section if `make_job_args.py` was run with `--replace-xcache` and `condor/job_args.txt` shows paths starting with `filelists/`.

---

## 8. Package the SIDM code

Whenever files under `sidm/` or `condor/run_sidm_chunk.py` change, remake the tarball.
The `sidm/studies/` exclude matters: it holds committed notebooks with embedded plot
outputs (hundreds of MB) the workers never use -- without it the tarball balloons ~100x
and per-job input transfer slows accordingly:

```bash
cd /uscms_data/d3/$USER/SIDM

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
  --exclude="sidm_venv" \
  --exclude="py39_packages" \
  --exclude="sidm/studies" \
  -czf condor/sidm_code.tar.gz \
  sidm condor/run_sidm_chunk.py
```

Check:

```bash
ls -lh condor/sidm_code.tar.gz
```

You do not need to re-tar if you only changed:

```
condor/requirements.txt
condor/submit.sub
condor/run_job.sh
condor/job_args.txt
condor/filelists/*.txt
```

Those are transferred separately.

---

## 9. Create EOS output directories

Use separate directories for signal chunks, background chunks, and merged outputs.

```bash
xrdfs root://cmseos.fnal.gov mkdir -p /store/user/$USER/sidm_condor/SignalChunks_v1
xrdfs root://cmseos.fnal.gov mkdir -p /store/user/$USER/sidm_condor/BackgroundChunks_v1

xrdfs root://cmseos.fnal.gov mkdir -p /store/user/$USER/sidm_condor/SignalMerged_v1
xrdfs root://cmseos.fnal.gov mkdir -p /store/user/$USER/sidm_condor/BackgroundMerged_v1
```

---

## 10. Condor submit file

File:

```
condor/submit.sub
```

Example:

```
universe = vanilla

executable = run_job.sh

arguments = $(sample) $(chunk) $(filelist) root://cmseos.fnal.gov//store/user/$ENV(USER)/sidm_condor/BackgroundChunks_v1

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
```

`$ENV(USER)` is expanded at submit time, so the same file works for any LPC user. The only thing to edit between runs is the trailing "run tag" — for signals use:

```
root://cmseos.fnal.gov//store/user/$ENV(USER)/sidm_condor/SignalChunks_v1
```

For backgrounds:

```
root://cmseos.fnal.gov//store/user/$ENV(USER)/sidm_condor/BackgroundChunks_v1
```

---

## 11. Proxy setup

Before submitting jobs:

```bash
voms-proxy-init --valid 192:00 -voms cms
```

`voms-proxy-init` writes to `/tmp/x509up_u<UID>` by default, but `/tmp` is **per-node** — the schedd that runs your jobs lives on a different host (`lpcscheddN.fnal.gov`) and cannot read your interactive node's `/tmp`. Jobs submitted with `x509userproxy=/tmp/...` will land in HOLD with:

```
Transfer input files failure ... reading from file /tmp/x509up_u<UID>: (errno 2) No such file or directory
```

Copy the proxy to a shared NFS location and point `X509_USER_PROXY` at the copy instead — that's what `submit.sub` references via `$ENV(X509_USER_PROXY)`:

```bash
cp /tmp/x509up_u$(id -u) /uscms_data/d3/$USER/x509_proxy.pem
chmod 600 /uscms_data/d3/$USER/x509_proxy.pem
export X509_USER_PROXY=/uscms_data/d3/$USER/x509_proxy.pem

voms-proxy-info -file $X509_USER_PROXY -timeleft
```

(The dask-from-notebook path via `lpcjobqueue` does **not** need this copy — its submit description uses `use_x509userproxy = true`, which lets condor handle the spooling automatically.)

---

## 12. Submit jobs

Submit from inside the `condor/` directory:

```bash
cd /uscms_data/d3/$USER/SIDM/condor

export X509_USER_PROXY=/uscms_data/d3/$USER/x509_proxy.pem   # the shared-NFS copy from §11

condor_submit submit.sub
```

Example output:

```
Attempting to submit jobs to lpcschedd5.fnal.gov
5787 job(s) submitted to cluster 59092403.
```

Record:

```
schedd:  lpcschedd5.fnal.gov
cluster: 59092403
```

---

## 13. Check job status

Use the schedd and cluster ID printed by `condor_submit`.

```bash
condor_q -name lpcschedd5.fnal.gov 59092403
```

Summary:

```bash
condor_q -name lpcschedd5.fnal.gov 59092403 -totals
```

Held jobs:

```bash
condor_q -name lpcschedd5.fnal.gov 59092403 -hold
```

Watch live:

```bash
watch -n 10 'condor_q -name lpcschedd5.fnal.gov 59092403 -totals'
```

If jobs disappear from `condor_q`, they are no longer active. Check logs and EOS outputs.

---

## 14. Check logs

From the `condor/` directory:

```bash
cd /uscms_data/d3/$USER/SIDM/condor

grep -R "return value 0" logs/*59092403*.log | wc -l
grep -R "return value 1" logs/*59092403*.log | wc -l
grep -R "return value 2" logs/*59092403*.log | wc -l

grep -R "Traceback\|ERROR\|Exception\|No such file" logs/*59092403*.err | tail -100
```

A successful job should have:

```
Normal termination (return value 0)
```

---

## 15. Compare expected outputs to actual EOS outputs

From the repo parent:

```bash
cd /uscms_data/d3/$USER/SIDM

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
```

For signals, replace `BackgroundChunks_v1` with `SignalChunks_v1`.

---

## 16. Retry failed jobs from missing ROOT files

A common failure is:

```
OSError: Failed to open file:
No such file or directory
```

This means one ROOT file in a chunk is missing on EOS. One missing ROOT file causes the entire chunk to fail.

### Get failed ProcIds

```bash
cd /uscms_data/d3/$USER/SIDM

condor_history -name lpcschedd5.fnal.gov \
  -constraint 'ClusterId == 59092403 && ExitCode != 0' \
  -af ProcId ExitCode > failed_procs_59092403.txt

wc -l failed_procs_59092403.txt
head failed_procs_59092403.txt
```

### Build filtered retry filelists

```bash
cd /uscms_data/d3/$USER/SIDM

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
```

### Submit retry jobs

Create a retry submit file:

```bash
cp condor/submit.sub condor/submit_retry_59092403.sub
```

Edit the final queue line:

```
queue sample, chunk, filelist from failed_retry_job_args_59092403.txt
```

Submit:

```bash
cd /uscms_data/d3/$USER/SIDM/condor

export X509_USER_PROXY=$(voms-proxy-info -path)

condor_submit submit_retry_59092403.sub
```

---

## 17. Merge `.coffea` chunks from EOS into the shared lpcmetx area

Chunks stay in `/store/user/$USER/sidm_condor/...`. Merged outputs go to `/store/group/lpcmetx/SIDM/coffea_outputs/$USER/<study>/`, which is group-writable and world-readable under `lpcmetx`. Pass `--filelists-dir` + `--selections` + `--hist-collections` to also emit a `.meta.yaml` sidecar describing each merged file:

```bash
cd /uscms_data/d3/$USER/SIDM
source sidm_venv/bin/activate
```

### Merge backgrounds

```bash
python sidm/scripts/merge_coffea_chunks_eos.py \
  --input-eos-dir  /store/user/$USER/sidm_condor/BackgroundChunks_v1 \
  --output-eos-dir /store/group/lpcmetx/SIDM/coffea_outputs/$USER/backgrounds_v1 \
  --workdir merge_tmp_background \
  --filelists-dir condor/filelists \
  --selections base \
  --hist-collections muon_base \
  --schema LLPNanoAODSchema --chunksize 50000 --unweighted-hist
```

### Merge signals

```bash
python sidm/scripts/merge_coffea_chunks_eos.py \
  --input-eos-dir  /store/user/$USER/sidm_condor/SignalChunks_v1 \
  --output-eos-dir /store/group/lpcmetx/SIDM/coffea_outputs/$USER/signals_v1 \
  --workdir merge_tmp_signal \
  --filelists-dir condor/filelists \
  --selections base \
  --hist-collections muon_base \
  --schema LLPNanoAODSchema --chunksize 50000 --unweighted-hist
```

### What the sidecar contains

Each `.meta.yaml` records the input ROOT file list, the full selection + hist-collection definitions, the schema, chunksize, `unweighted_hist` flag, the SIDM git commit, the coffea version, and a UTC timestamp. Loaded via `sidm.tools.metadata.load_run_metadata()`.

### Check merged outputs

```bash
xrdfs root://cmseos.fnal.gov ls /store/group/lpcmetx/SIDM/coffea_outputs/$USER/backgrounds_v1
xrdfs root://cmseos.fnal.gov ls /store/group/lpcmetx/SIDM/coffea_outputs/$USER/signals_v1
```

Expected: one `.coffea` and one `.meta.yaml` per merged sample.

### Chunk retention and cleanup

The per-chunk `.coffea` files under `/store/user/$USER/sidm_condor/<v>Chunks_<n>/` are not deleted by the merge script by default. Re-merging is cheap if the chunks remain, but they accumulate EOS storage over time.

Pass `--delete-chunks-after-merge` to the merge invocation to remove the input chunks once the merged file and sidecar have been xrdcp'd to EOS. The flag is opt-in because the deletion is irreversible; run a merge first and verify the output before passing it.

---

## 18. Access merged files in CMSLPC notebooks (and across collaborators)

Merged outputs live under `/store/group/lpcmetx/SIDM/coffea_outputs/$USER/<study>/` on EOS, world-readable under the `lpcmetx` group area. `coffea.util.load` / `load_run_metadata` open **local files** (not `root://` URLs), and LPC recommends accessing EOS over **xrootd** rather than the `/eos/uscms` FUSE mount — so `xrdcp` the files down first, then load the local copies. Set the producer to your own `$USER`, or a collaborator's name to read theirs:

```python
import os, glob, subprocess
from coffea.util import load
from sidm.tools.metadata import load_run_metadata

producer, study = os.environ["USER"], "signals_v1"   # collaborator: set producer to their name
redir   = "root://cmseos.fnal.gov"
eos_dir = f"/store/group/lpcmetx/SIDM/coffea_outputs/{producer}/{study}"
local   = f"/uscms_data/d3/{os.environ['USER']}/SIDM/imported/{producer}_{study}"
os.makedirs(local, exist_ok=True)

for r in subprocess.run(["xrdfs", redir, "ls", eos_dir], check=True,
                        capture_output=True, text=True).stdout.split():
    if r.endswith((".coffea", ".meta.yaml")):          # copy the .coffea AND its sidecar
        subprocess.run(["xrdcp", "-f", f"{redir}/{r}", f"{local}/{os.path.basename(r)}"], check=True)

for path in sorted(glob.glob(f"{local}/*.coffea")):
    out  = load(path)                 # local file; dict {out, processed, exception}
    meta = load_run_metadata(path)    # reads the sibling local .meta.yaml
    print(path.split("/")[-1], "->", list(out.get("out", out).keys()), "| sidm_commit:", meta["sidm_commit"])
    for s in meta["samples"]:         # n_files / xsec_pb are per-sample, under "samples"
        print("   ", s["name"], s["n_files"], "files, xsec_pb =", s["xsec_pb"])
```

Reading these world-readable files works with your Kerberos ticket; EOS writes need a valid VOMS proxy. The two example notebooks (`sidm/test_notebooks/lpc_dask_example.ipynb`, `lpc_condor_example.ipynb`) use this exact pattern, with a guard that raises a clear error if `producer`/`study` is wrong or the directory has no `.coffea`.

---

## 19. Files that should be committed

Commit reusable scripts and templates:

```bash
git add condor/README.md
git add condor/make_job_args.py
git add condor/run_sidm_chunk.py
git add condor/run_job.sh
git add condor/submit.sub
git add condor/requirements.txt
git add condor/signal_samples.txt
git add condor/background_samples.txt
git add merge_coffea_chunks_eos.py
```

Also commit relevant processor changes:

```bash
git add sidm/tools/sidm_processor.py
```

Do not commit generated files or credentials.

---

## 20. Files that should not be committed

Do not commit:

```
*.coffea
*.root
*.parquet
*.json
*.pem
x509up_*
sidm_venv/
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
```

Add these to `.gitignore`.

---

## 21. Common errors

### `No module named coffea`

Cause: Condor worker Python does not have Coffea.

Fix: ensure `run_job.sh` creates the venv and installs:

```bash
python -m pip install --no-cache-dir -r requirements.txt
```

### `KeyError: is_data` (or `skim_factor`, `year`)

Cause: the fileset handed to `processor.Runner.run` lacks a `metadata` block. `SidmProcessor.process` reads `events.metadata["is_data"]` / `["skim_factor"]` / `["year"]` without `.get()` defaults.

Fix: in `run_sidm_chunk.py`, build the fileset as

```python
fileset = {args.sample: {"files": files, "metadata": {"is_data": False, "skim_factor": 1.0, "year": "2018"}}}
```

Note that `"dataset"` is a reserved key in `metadata` for newer coffea (≥ 2025.x) and must not be supplied — coffea sets it from the fileset key.

### Missing ROOT file

Cause: file listed in YAML/filelist is not present on EOS.

Fix: use the retry workflow to remove missing ROOT files from failed chunks and resubmit.

### `XRootD: [ERROR] Operation expired` / `ValueError: Empty list provided to reduction`

Cause: EOS is timing out file opens (degraded or congested), so the reads fail. This is not a code or cluster fault, and it hits every read path — Condor jobs, dask workers, and local runs all read the same EOS. `skipbadfiles=True` only rescues a run when *some* files are bad; when every open fails there are no chunks left to process.

Fix: confirm EOS is the problem with `xrdfs cmseos.fnal.gov stat /store/group/lpcmetx/SIDM` (it hangs when EOS is degraded), wait for it to recover, then resubmit (or use the retry workflow in §16).


---

## 22. Orchestrated reconcile and report (`condor_campaign.py`)

Sections 13–16 walk through checking a campaign by hand: `condor_q`, grepping the
logs for `return value`, an `awk`/`comm` of expected vs actual EOS outputs, and a
copy-paste loop to rebuild failed chunks. `condor/condor_campaign.py` automates the
*read-only* half of that — it reconciles what you submitted against the logs and
EOS, classifies every failure by root cause, and writes one report — in a single
command. It submits and resubmits nothing (that is a later phase), so it is safe to
run at any time, including while jobs are still going.

It is standard-library only; run it on the submitter in the `sidm_venv` (a bare
`python3` works too):

```bash
cd /uscms_data/d3/$USER/SIDM
python condor/condor_campaign.py reconcile \
    --job-args      condor/job_args.txt \
    --logs-dir      condor/logs \
    --eos-chunk-dir /store/user/$USER/sidm_condor/SignalChunks_v1 \
    --run-id        signals_v1
```

This writes `condor/campaigns/signals_v1/report.json` and `report.md` (and prints
the report). Re-render or summarize an existing report without re-scanning:

```bash
python condor/condor_campaign.py report --report condor/campaigns/signals_v1/report.json
python condor/condor_campaign.py status --report condor/campaigns/signals_v1/report.json
```

### What the report tells you

- Campaign totals: expected / DONE (present on EOS, non-empty) / failed / terminal.
- A breakdown of the non-DONE chunks by failure class. Each class is the root cause,
  read from the durable per-job `.log` and `.err` (so it still works after
  `condor_history` has aged out):
  - `MISSING_ROOT_FILE` — an input ROOT file is **confirmed absent** on EOS now (an
    `xrdfs stat` returned "No such file"); a timeout or any other error is treated as
    `EOS_TIMEOUT`, never as a missing file.
  - `EOS_TIMEOUT` — a transient EOS read failure; inputs are still present (or could
    not be probed). Resubmit the whole chunk unchanged once EOS is healthy.
  - `OOM` — killed for memory; resubmit at higher `request_memory`.
  - `HELD_PROXY` — held because the x509 proxy was unreadable by the schedd (the
    `/tmp` proxy problem; see §11). Renew + copy the proxy, then resubmit unchanged.
  - `CODE_BUG` — a deterministic python exception (e.g. a `KeyError`). A plain retry
    re-fails; fix the processor and rebuild `sidm_code.tar.gz` first.
  - `STALLED` / `HELD_TRANSFER` / `XRDCP_FAILED` / `MISSING_OUTPUT` /
    `ZERO_BYTE_OUTPUT` / `UNKNOWN` / `LOST` — see each chunk reason in the report.
- Per-sample completeness, listing the missing chunk indices, so you never merge a
  sample that is not yet complete.
- A "what to do next" block derived from the class histogram.
- A per-chunk **NON-DONE CHUNK DETAIL** section listing each failed or stalled chunk, its class and root-cause reason (the "needs a human" ones first), with a tail of the job stdout/stderr so you can see what it was doing before it failed -- e.g. a stalled chunk shows how long it made no progress before removal.

`--no-stat-missing` skips the per-file `xrdfs stat`, which is faster but cannot tell
`MISSING_ROOT_FILE` from `EOS_TIMEOUT` (everything that open-failed is reported as
`EOS_TIMEOUT`, i.e. nothing is ever blamed on a missing file). A unit suite covering
the parsers and every failure class is in `condor/test_campaign_lib.py`
(`python condor/test_campaign_lib.py`).

Generated reports live under `condor/campaigns/<run-id>/` and are git-ignored.

### Exit codes and mid-campaign safety

`reconcile` / `status` / `report` exit `0` only when **every** chunk is DONE, `1`
when something is not DONE (confirmed failures, in-flight, or presence-unknown), `2`
when the EOS dir could not be listed, and `3` for both — so
`condor_campaign.py status ... && merge ...` will not merge an incomplete campaign.
Chunks still queued or running are reported `IN_PROGRESS` and are **not** counted as
failures, so the reconcile is safe to run mid-campaign. `--eos-chunk-dir` takes a bare
`/store/...` path (a full `root://…` URL is accepted and stripped); a typo is reported
as a "not-found" directory, not an EOS outage.
