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

This is **one** way to get a working interactive environment on Fermilab LPC (`cmslpc-el9.fnal.gov`). Other paths exist — e.g. the [Elastic Analysis Facility (EAF)](https://analytics-hub.fnal.gov), JupyterHub on LPC, or plain shell sessions — and the venv-setup steps below (1–4) transfer directly to those. Only step 5 (how to actually open a notebook) is laptop-specific — it offers two options, **JupyterLab in a browser** (5a) or **VS Code over SSH** (5b); on EAF or JupyterHub you skip it and pick the registered kernel from the web UI instead.

This walkthrough covers the **interactive / notebook** environment. For running large batch jobs on LPC Condor, see [condor/README.md](condor/README.md) once steps 1–4 are done.

**Which path to use, by scale:**

| Scale | Path | Section |
|---|---|---|
| 1–few files, quick iteration | Interactive notebook (one worker on cmslpc-el9) | Steps 1–6 |
| 10s–100s of files, want to keep iterating | Dask-from-notebook (elastic Condor workers, single notebook session) | Step 7 |
| 1000s of files, fire-and-forget | Condor scripts (unattended overnight runs) | [condor/README.md](condor/README.md) |

You can climb this ladder as a study grows — the same notebook can graduate from interactive → dask without restructuring (see step 7.1, "Adapting an existing studies/ notebook").

**Tutorial notebooks** cover the full setup (terminal commands, processor run, sanity plot) for both scaled paths:

- Dask-from-notebook: [`sidm/test_notebooks/lpc_dask_example.ipynb`](sidm/test_notebooks/lpc_dask_example.ipynb)
- Condor batch: [`sidm/test_notebooks/lpc_condor_example.ipynb`](sidm/test_notebooks/lpc_condor_example.ipynb)

Both are committed with outputs.

**Shared output area + sidecar metadata.** Both tutorials push their merged `.coffea` to `/store/group/lpcmetx/SIDM/coffea_outputs/$USER/<study>/` on EOS (CMS-LPC group-writable, world-readable) and write a `.meta.yaml` sidecar next to each `.coffea` recording the input file list, the full selection + hist-collection definitions, the SIDM git commit, the coffea version, and a UTC timestamp. The helpers are in [`sidm/tools/metadata.py`](sidm/tools/metadata.py): `write_run_metadata` and `load_run_metadata`.

### Prerequisites (one-time, from your laptop)

Three pieces of CMS/FNAL access are usually arranged through your institution's CMS computing rep before you can SSH to LPC:

1. **FNAL services account** with `cmslpc` computing access — gives you a username (e.g. `jsmith@FNAL.GOV`) and a home directory under `/uscms_data/d3/jsmith/`.
2. **CMS VO membership** — required for VOMS proxies and `/store/group/lpcmetx/...` access.
3. **Grid certificate** loaded into your laptop's certificate store (Keychain on macOS, Certificate Store on Windows; typically obtained from [CILogon](https://cilogon.org)) — used on the LPC side when you run `voms-proxy-init`.

You also need a Kerberos client on your laptop with an FNAL principal:

```bash
# macOS has Kerberos built in; Linux: install krb5-user (Debian/Ubuntu) or krb5-workstation (RHEL/Alma).
kinit YOUR_FNAL_USERNAME@FNAL.GOV
ssh cmslpc-el9.fnal.gov   # should now succeed without a password prompt
```

> **On Windows, run all laptop-side steps inside [WSL2](https://learn.microsoft.com/windows/wsl/install) (a real Ubuntu).** There `sudo apt install krb5-user` gives a working `kinit`, SSH's Kerberos/GSSAPI auth behaves, and every bash block below — including the `\` line-continuations — runs verbatim; `localhost` port-forwarding still reaches your Windows browser. Native PowerShell/cmd trip on both: there is no out-of-the-box `kinit`, and `\` is **not** a line-continuation (PowerShell uses `` ` ``), so a pasted multi-line `ssh -L … \` command breaks apart and the `-L` forward silently never takes effect — the server starts on LPC but the browser can't reach it. If you must use native PowerShell, put each command on a single line and keep the server in the same session that holds the `-L` forward (or let VS Code **Remote-SSH** manage the tunnel for you).

The canonical starting page for new CMS-LPC users is [https://uscms.org/uscms_at_work/computing/getstarted/](https://uscms.org/uscms_at_work/computing/getstarted/). Once these are in place, every step below happens on `cmslpc-el9.fnal.gov`.

> **Refresh your VOMS proxy at the start of each work session** (lasts 192 h / 8 days):
> ```bash
> voms-proxy-init --valid 192:00 -voms cms
> ```
> Needed for reading remote ROOT files over XRootD, submitting Condor jobs, and spawning dask workers via `lpcjobqueue`. The notebook helpers in `sidm/tools/scaleout.py` call `check_voms_proxy()` first and fail with this exact command if your proxy is missing or about to expire.

### 1. Clone the repo

SSH to LPC and clone into your work area:

```bash
ssh cmslpc-el9.fnal.gov
cd /uscms_data/d3/$USER
git clone <your-fork-of-SIDM>
cd SIDM
```

`/uscms_data/d3/$USER/SIDM` is the permanent home for this checkout — a standalone directory in your data area with `sidm_venv/` built inside it. The venv records this absolute path in every console script (`pip`, `jupyter`, …) and in the Jupyter kernelspec, so **build the venv here and do not move or rename the folder afterward**; a relocated venv fails with `bad interpreter: …/sidm_venv/bin/python: No such file or directory`. If you must relocate, move the checkout first, then delete `sidm_venv/` and redo steps 2–4 (re-registering the kernel in step 4).

### 2. Build the venv on LCG_107 Python 3.11

The repo's `requirements.txt` pins versions (`dask`, `awkward`, …) that need Python ≥ 3.10. LPC's system Python is 3.9, so we use Python 3.11 from cvmfs:

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh
unset PYTHONPATH
python -m venv sidm_venv
source sidm_venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
```

After creation the venv's `bin/python` is a symlink to the cvmfs binary, so you don't need to re-source LCG_107 to run it later.

### 3. Install requirements

```bash
python -m pip install -c constraints.txt -r requirements.txt
python -m pip install -c constraints.txt "distributed==2025.3.0"   # required by sidm/tools/scaleout.py
python -m pip install "htcondor<25" "git+https://github.com/CoffeaTeam/lpcjobqueue.git@4c6dae0b85838be0c42f7fa7d7dad0527920d56b"  # needed by step 7
python -m pip install -e .
python -m pip install jupyter ipykernel pyarrow
```

`xrootd` installs as an ordinary wheel and needs no special handling. It used to:
version 5.7.0 was published as a source distribution only, so pip tried to compile
the XRootD C library and failed on LPC, and the recipe worked around that by
stripping the line out of `requirements.txt` and symlinking LCG's prebuilt
`XRootD` and `pyxrootd` into the venv. Upstream has shipped manylinux wheels since
6.x, so both workarounds are gone. `fsspec-xrootd` pulls `xrootd` in anyway, which
is why the pinned version is the one that was already being installed.

The `-c constraints.txt` flag pins every transitive dependency as well as the
direct ones, so you get the same environment the CI tests against rather than
whatever PyPI happens to serve today. `requirements.txt` remains the file to
edit when the analysis needs a new package; `constraints.txt` is generated, and
its header documents how to regenerate it after a deliberate bump. Condor
workers use their own `condor/constraints.txt`, because they run CMSSW's Python
3.9 and resolve to a different set.

This matters more than it looks. Before these files existed, roughly 68 of the
84 packages in this venv were unpinned, and on 2026-08-11 a numba release began
miscompiling one of coffea's kernels and broke every CI run until it was pinned
by hand.

Verify the scale-out dependencies used by step 7 are all installed (no `MISSING:` list = success):

```bash
python -c "import importlib.util as u; miss=[x for x in ['bokeh','lpcjobqueue','dask_jobqueue','distributed','htcondor'] if u.find_spec(x) is None]; print('MISSING: '+', '.join(miss) if miss else 'scale-out deps OK')"
```

A name in the `MISSING:` list means that install line above did not land (commonly the `lpcjobqueue` git install) — rerun it. This checks the packages are present without importing them, so it stays fast and avoids the harmless "Condor configuration not found!" notice that `import lpcjobqueue` prints when `CONDOR_CONFIG` is unset.

Confirm `sidm` installed **editable** (an in-place link to your checkout, not a copy):

```bash
sidm_venv/bin/pip show sidm | grep Editable    # -> "Editable project location: .../SIDM"
```

If that line is missing, a non-editable copy of `sidm` is shadowing your checkout in `site-packages`. `import sidm` then resolves to stale code whenever you run a standalone script — `python script.py` puts the script's directory on `sys.path`, not the repo root, so the working tree isn't found (notebooks and `python -c` from the repo root are unaffected). Reinstall editable:

```bash
sidm_venv/bin/pip uninstall sidm -y && sidm_venv/bin/pip install -e .
```

Smoke-test that uproot can read a remote ROOT file through xrootd:

```bash
sidm_venv/bin/python -c "import uproot; print(uproot.open('root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/ULSignalSamples/2018_v10/BsTo2DpTo2Mu2e/CutDecayFalse_SIDM_BsTo2DpTo2Mu2e_MBs-500_MDp-1p2_ctau-1p9_v3/LLPnanoAODv2/CutDecayFalse_SIDM_BsTo2DpTo2Mu2e_MBs-500_MDp-1p2_ctau-1p9_v3_part-0.root')['Events'].num_entries)"
```

### 4. Register the Jupyter kernel

```bash
/uscms_data/d3/$USER/SIDM/sidm_venv/bin/python -m ipykernel install --user \
    --name sidm_venv \
    --display-name "SIDM (LCG_107 Py3.11)"
```

This drops a kernelspec under `~/.local/share/jupyter/kernels/sidm_venv/` on LPC. Any Jupyter Server running on LPC — including EAF or JupyterHub — will offer this kernel.

### 5. Open notebooks

The kernel runs on LPC (its binary is x86_64 Linux), so you open notebooks through a frontend that connects to a Jupyter server on LPC. Two patterns work, and both reach the **SIDM (LCG_107 Py3.11)** kernel from step 4: **JupyterLab in a browser** (5a — simplest) and **VS Code over SSH** (5b).

#### 5a. JupyterLab in a browser (simplest)

From a terminal on your laptop, start a JupyterLab server on LPC and forward its ports:

```bash
ssh -L 8888:localhost:8888 -L 8787:localhost:8787 cmslpc-el9.fnal.gov \
    "export X509_USER_PROXY=/uscms_data/d3/$USER/x509_proxy.pem && \
     cd /uscms_data/d3/$USER/SIDM && source sidm_venv/bin/activate && \
     jupyter lab --no-browser --port=8888 --ip=127.0.0.1"
```

Open the URL it prints — `http://127.0.0.1:8888/lab?token=…` (`localhost` also works over the tunnel) — in your laptop browser and pick the **SIDM (LCG_107 Py3.11)** kernel. `--no-browser` is required (LPC is headless) and the token in the URL logs you in. The second `-L 8787` forwards the Dask dashboard (step 7) to `http://localhost:8787/status` — drop it if you don't need it.

> **VOMS proxy — mint it from a JupyterLab terminal.** Your dask/condor cells need a valid CMS proxy that the **kernel** can find. Two things bite: the `ssh "…"` launch runs a non-interactive shell that does **not** source your `.bashrc` (so the kernel won't inherit a proxy path set there), and the default `/tmp/x509up_u<UID>` is per-node. So mint the proxy to the fixed NFS path that the launch's `X509_USER_PROXY` already points at — from a terminal **inside JupyterLab** (File ▸ New ▸ Terminal), which runs on your kernel's node:
> ```bash
> voms-proxy-init --valid 192:00 -voms cms -out /uscms_data/d3/$USER/x509_proxy.pem
> chmod 600 /uscms_data/d3/$USER/x509_proxy.pem
> ```
> `scaleout.check_voms_proxy()` and the Condor workers then resolve the same file. Forgot the launch export? Set it from the kernel in the notebook before `check_voms_proxy()` (it runs on LPC, so `os.environ['USER']` resolves): `import os; os.environ["X509_USER_PROXY"] = f"/uscms_data/d3/{os.environ['USER']}/x509_proxy.pem"`.

**List and stop servers.** `jupyter server list` (equivalently `jupyter lab list`) prints each running server with its port and token. To shut one down, press **Ctrl+C** in the launch terminal, or from a shell on the same node run `jupyter lab stop 8888` (equivalently `jupyter server stop 8888`).

> **Port and same-node gotchas.** Port 8888 is per-node. If a previous launch still holds it, Jupyter quietly starts on the next free port (8889, 8890, …) — which your `-L 8888` tunnel isn't forwarding; run `jupyter server list`, stop the stray, and relaunch (change the port in *both* the `-L` and `--port` if you change it). And because `cmslpc-el9.fnal.gov` load-balances across many nodes (`cmslpc101`, `cmslpc336`, …), `jupyter … stop` and the in-Lab terminal only reach a server that's on their node — keeping the tunnel and server in one `ssh -L … "…"` command keeps everything together.

#### 5b. VS Code over SSH

VS Code attaches to the same kind of server — launch it like 5a but with `jupyter server` in place of `jupyter lab` (the `-L 8787` is optional unless you want the dashboard):

```bash
ssh -L 8888:localhost:8888 cmslpc-el9.fnal.gov \
    "export X509_USER_PROXY=/uscms_data/d3/$USER/x509_proxy.pem && \
     cd /uscms_data/d3/$USER/SIDM && source sidm_venv/bin/activate && \
     jupyter server --no-browser --port=8888 --ip=127.0.0.1"
```

Copy the URL it prints (`http://localhost:8888/?token=…`). In VS Code: open any notebook → kernel selector (top-right) → **Select Another Kernel…** → **Existing Jupyter Server…** → paste URL → pick **SIDM (LCG_107 Py3.11)**. The VOMS-proxy note, the list/stop note, and the port/same-node gotchas from 5a apply identically — except there is no in-Lab terminal in the bare `jupyter server` flow, so mint the proxy from a separate `ssh` shell on the **same node** as this server (or use the in-notebook `os.environ["X509_USER_PROXY"]` override).

**Alternative — let VS Code manage the connection (Remote-SSH).** Instead of the manual `ssh -L` tunnel above, install the **Remote - SSH** extension and run **Remote-SSH: Connect to Host…** → `cmslpc-el9.fnal.gov` from the Command Palette. VS Code runs its server component on LPC, so it sees the x86_64 kernel natively: open the `SIDM` folder, open a notebook, pick **SIDM (LCG_107 Py3.11)** — no port-forward, no token to paste, and none of the PowerShell pitfalls above. The integrated terminal also becomes a Linux shell on LPC, so the bash blocks in steps 1–4 run there directly. (Auth is still Kerberos — on Windows the WSL2 note above still applies — and since `cmslpc-el9` is round-robin you'll land on whichever node it picks.) Because VS Code — not an `ssh "…"` command — launches the kernel here, the `X509_USER_PROXY` export does not reach it; set the proxy from inside the kernel instead, via the `os.environ["X509_USER_PROXY"]` override in the dask notebook's proxy-check cell before `check_voms_proxy()`.

If you're on EAF, JupyterHub, or running Jupyter directly on LPC with your own tunneling, skip this step entirely and just pick the kernel from your usual UI.

### 6. Reading remote ROOT files

The location YAMLs (`sidm/configs/ntuples/*.yaml`) use `root://xcache//…` URLs. `xcache` is the coffea-casa cache and does not resolve on LPC. Pass `replace_xcache=True` to `utilities.make_fileset` to substitute the FNAL EOS redirector at fileset-build time:

```python
fileset = utilities.make_fileset(
    samples, "llpNanoAOD_v2", location_cfg="signal_2mu2e_v10.yaml",
    replace_xcache=True,
)
```

Default is `False` so existing coffea-casa notebooks remain unchanged.

### 7. Scaling jobs from a notebook with dask + HTCondor

For studies that exceed what a single interactive node can handle (≳ 4 GB memory, a handful of cores), HTCondor workers can be fanned out from the notebook. The pattern matches coffea-casa: a `Client` connected to a Dask cluster, passed to `processor.DaskExecutor`:

```python
from sidm.tools import scaleout
from coffea import processor

scaleout.check_voms_proxy()                       # fails loudly if proxy is missing/expiring

cluster, client = scaleout.make_lpc_client(
    min_workers=2, max_workers=10,
    memory="4GB", disk="4GB",
)

runner = processor.Runner(
    executor=processor.DaskExecutor(client=client),
    schema=...,
    chunksize=50_000,
)
output = runner.run(fileset, treename="Events", processor_instance=p)

cluster.close()
```

Under the hood: workers are HTCondor jobs that run inside the `coffea-dask-almalinux9:2025.5.0.rc2-py3.11` apptainer image on the LPC pool; the notebook itself stays in `sidm_venv` outside the apptainer. The key convenience over coffea-casa: your local `sidm/` code (including uncommitted edits) is shipped to each worker via `UploadDirectory`, so there is no commit/push-to-GitHub roundtrip (and no branch to pin) during iteration — re-run `make_lpc_client()` after editing anything under `sidm/` to refresh the workers. Only the importable tree travels: `sidm/tools`, `sidm/definitions`, `sidm/scripts`, `sidm/configs` and `sidm/data`. Notebooks and their committed outputs (`sidm/studies`, `sidm/test_notebooks`) are excluded, because the scheduler runs inside the notebook process and re-sends the payload to every worker: shipping `sidm/studies` has OOM-killed kernels on shared interactive nodes. **Code that runs inside a dask task must live in one of the uploaded directories, not next to a notebook in `sidm/studies/`.**

The required Python packages (`htcondor<25`, `lpcjobqueue`) were installed as part of step 3, so no extra setup is needed. A runnable end-to-end example lives at [sidm/test_notebooks/lpc_dask_example.ipynb](sidm/test_notebooks/lpc_dask_example.ipynb).

`scaleout.check_voms_proxy()` (called at the top of the example, and again inside `make_lpc_client`) needs `X509_USER_PROXY` visible to the **kernel** — see the proxy note in step 5. The same exported file is spooled to the Condor workers via `use_x509userproxy`, so one export covers both.

To watch the Dask dashboard from your laptop browser, forward its port too: the 5a launch already includes `-L 8787:localhost:8787` (if you started from 5b or dropped it, add it to the `ssh` launch), then open `http://localhost:8787/status`. `make_lpc_client()` pins the dashboard link to `localhost` so the URL shown by the cluster (and printed by the example notebook's `print('dashboard:', cluster.dashboard_link)` cell) matches this tunnel — otherwise `lpcjobqueue` advertises a relative `/proxy/8787/status` link that only resolves under a JupyterHub, not over a plain SSH forward. The repr shows the actual port if 8787 was already taken. To move the dashboard off 8787 deliberately — e.g. when a coffea-casa Dask dashboard is already on that port — pass `make_lpc_client(dashboard_port=8790)` (any free port) and forward that **same** port in your launch (`-L 8790:localhost:8790`). The two sides of the tunnel must use the same number: the printed link uses the scheduler's port, so remapping only the `ssh -L` local port leaves the link pointing at the old one. A pinned port is still only a request — if it too is busy on the node, distributed warns and rebinds to a random free port, so forward whatever port `cluster.dashboard_link` prints. The dashboard needs the `bokeh` package (pinned in `requirements.txt`); without it `localhost:8787` still loads but shows a *"Dask needs bokeh>=3.1.0 for the dashboard"* notice instead of the live charts. One caveat: the **Groups** tab (Task Group Graph) can log a harmless `KeyError: 'automatic_retries'` — a bug in distributed's dashboard layout, triggered by coffea's retry wrapper rather than by your run; the main **Status**, task-stream, and worker panels are unaffected, as is the computation itself. The dashboard is optional — the run completes without it.

**Saving the output.** Write the `.coffea` (and a `.meta.yaml` metadata sidecar via `sidm.tools.metadata`) to the group-writable lpcmetx EOS area so it persists and other CMS users can re-read it:

```bash
xrdfs cmseos.fnal.gov mkdir -p /store/group/lpcmetx/SIDM/coffea_outputs/$USER/<study>
xrdcp -f out.coffea root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/coffea_outputs/$USER/<study>/out.coffea
```

**Reading a saved output back.** `coffea.util.load` does a local `open()` and rejects `root://` URLs, so `xrdcp` the `.coffea` to local scratch first, then load the local copy. The redirector to use depends on where you run:

```python
import os, subprocess, tempfile
from coffea.util import load

lfn   = "/store/group/lpcmetx/SIDM/coffea_outputs/<user>/<study>/out.coffea"
redir = "root://cmseos.fnal.gov/"     # on LPC / EAF (needs a VOMS proxy or a Kerberos ticket)
# redir = "root://xcache/"            # on coffea.casa instead (auto-authenticated via SciToken)
local = os.path.join(tempfile.mkdtemp(), os.path.basename(lfn))
subprocess.run(["xrdcp", "-f", redir + lfn, local], check=True)
output = load(local)
```

The working redirector **flips between facilities**: on LPC use `root://cmseos.fnal.gov//…`; on **coffea.casa** use `root://xcache//…` (the same `cmseos`↔`xcache` swap as the inputs in §6). On coffea.casa the direct `cmseos` door fails (`Auth failed` — no proxy in the session) and the global `cms-xrd-global` redirector times out (group space is not advertised cross-site), so `xcache` is the one that works there — confirmed by running the notebook below. A runnable walkthrough that probes all three redirectors and loads the output plus its `.meta.yaml` sidecar is [`sidm/test_notebooks/coffea_casa_read_output.ipynb`](sidm/test_notebooks/coffea_casa_read_output.ipynb); to *find* which outputs exist in the first place, see [`sidm/test_notebooks/output_browser.ipynb`](sidm/test_notebooks/output_browser.ipynb).

**Troubleshooting — XRootD/EOS read failures.** A run that dies with `ValueError: Empty list provided to reduction` (or workers logging `[ERROR] Operation expired` while opening files) means EOS is timing out the reads, not a fault in the code or the cluster. It affects every executor — the local `IterativeExecutor`/`FuturesExecutor` and the Dask-over-Condor workers all read the same EOS. `skipbadfiles=True` only helps when *some* files are bad; if every open fails there are no chunks left to reduce. Confirm with `xrdfs cmseos.fnal.gov stat /store/group/lpcmetx/SIDM` (it hangs when EOS is degraded), then rerun once EOS recovers.

### 7.1. Adapting an existing studies/ notebook from coffea-casa to LPC

Most notebooks in `sidm/studies/` were written for coffea-casa. Three swaps move them to LPC:

1. **Drop the `importlib.reload(scaleout)` line.** That was a workaround for editing `sidm/tools/scaleout.py` between runs. With `make_lpc_client()` your local `sidm/` code is uploaded to workers automatically, so there is no source to edit. Note that helper modules sitting beside a notebook in `sidm/studies/` are **not** uploaded; if a dask task needs one, move it under `sidm/tools/` or `sidm/definitions/` first.

2. **Replace the client constructor.** Where the notebook has:
   ```python
   client = scaleout.make_dask_client("tls://localhost:8786")
   ```
   write:
   ```python
   scaleout.check_voms_proxy()
   cluster, client = scaleout.make_lpc_client(min_workers=2, max_workers=10)
   ```
   Add `cluster.close()` near the end of the notebook (after you've saved/loaded your output) so workers are released.

3. **Add `replace_xcache=True` to every `make_fileset` call.** Without it the URLs in the location YAMLs still point at `root://xcache//`, which does not resolve on LPC.
   ```python
   # Before:
   fileset = utilities.make_fileset(samples, "llpNanoAOD_v2", max_files=2)
   # After:
   fileset = utilities.make_fileset(samples, "llpNanoAOD_v2", max_files=2, replace_xcache=True)
   ```

Everything else (the `processor.Runner` setup, `DaskExecutor(client=client)`, the `runner.run(...)` call, all the cuts/hists logic) stays exactly the same.

### 8. Long-running batch jobs (no notebook attached)

For unattended runs without an open notebook, use the Condor scripts path: [condor/README.md](condor/README.md) is the reference, and [`sidm/test_notebooks/lpc_condor_example.ipynb`](sidm/test_notebooks/lpc_condor_example.ipynb) is the runnable walkthrough (sample list → tarball → submit → merge → load and plot).

The merge step is also where per-sample **normalization** is made correct: each Condor job normalizes only its own files (`lumi*xs/sumw_chunk`), and `merge_coffea_chunks_eos.py` re-stitches them into the full-sample normalization rather than inflating yields by the number of jobs — see [condor/README.md](condor/README.md) §17. (The dask path in §7 runs a single pass and normalizes once, so it needs no separate merge.)

## Code structure
All the interesting code in this repository is in `SIDM/sidm/`, which is organized into the following subdirectories:
- `tools` contains the classes and methods that form the backbone of SIDM. In particular, `sidm_processor.py` defines how events are analyzed. We use the NanoAODSchema for the data in our files, as defined within Coffea [here](https://coffea-hep.readthedocs.io/en/latest/api/coffea.nanoevents.NanoAODSchema.html)
- `definitions` contains files that define the current set of histograms one can make, the cuts one can apply, and some of the objects one can use when making histograms and applying cuts.
- `configs` contains yaml configuration files that define how histograms are grouped together into collections and cuts are grouped together into selections. When running `sidm_processor`, one provides names of selections and names of histogram collections to choose which cuts to apply and which histograms to make.
- `test_notebooks` contains notebooks to test new classes or functionalities as they are added. These notebooks can also serve as a form  of unit test: if you edit some code in a way that you think shouldn't affect the behavior, you can run these notebooks to confirm the output is unchanged.
- `studies` is where the physics happens. The notebooks in this directory are meant to serve effectively as pages in a lab notebook. The intention is to create a new notebook for each unique physics study and to include markdown comments to describe the intentions and observations of the person performing the study. These notebooks can also serve as unit tests in the same way as those in `test_notebooks`.

Outside `sidm/`, the repository-level `tests/` directory holds the per-PR chain-report CI harness and its committed fixture — see [`tests/README.md`](tests/README.md).

## Analysis how-tos

### General workflow
I suggest the following workflow for performing a physics study:
1. Create a new branch with a descriptive name (e.g. `leptonIsolation`)
2. Create a new notebook in `studies` with a descriptive name (e.g. `study_lepton_isolation.ipynb`)
3. Following the examples of the existing notebooks in `studies`, use `sidm_processor` to apply cuts and make histograms starting from an Firefighter ntuple of your choosing. Make sure to describe your reasoning and observations in text as you go.
4. If you find you need to define new selections, new cuts, or new histograms, follow the guides below.
5. Commit your changes as you go and submit a Pull Request once you have a reasonable standalone study or have added new selections, cuts, histograms, objects, classes, or features. Every PR automatically gets an advisory chain report comparing your branch to `main` (cutflows, hist collections, warnings); see [`tests/README.md`](tests/README.md).

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

### Error bars on plots

`utilities.plot` follows one convention for error bars:

- **Count histograms** (filled from event counts): leave `yerr` alone -- mplhep's default is correct. That is asymmetric Poisson (Garwood) bars for unweighted/unit-weight counts and symmetric `sqrt(sum w^2)` for lumi*xsec-weighted hists, so a raw-count plot showing *asymmetric* bars is expected, not a bug.
- **Derived quantities** (efficiency, ratio, scale factor, or any hist filled via `.view()`): **pass explicit `yerr`** -- otherwise mplhep draws `sqrt(bin_content)`, which is meaningless for a non-count value. `utilities.get_eff_hist(num, den)` returns the efficiency hist plus the matching `errors` to pass as `yerr`; `utilities.plot_ratio` does this for you (and falls back to a Gaussian interval when the ratio exceeds 1, where the binomial interval is undefined).

One runtime guard to know: `plot_MC_sig_vs_bkg_panels` warns *"background '<name>' has integer variances (looks unweighted) ... check that it is lumi*xsec-weighted"* when a background's per-bin variances are integers -- the hist looks like raw counts rather than lumi*xsec-weighted MC. If you meant to show normalized MC, confirm the weights were applied. It is a guard, not an error; the plot still renders.

## Continuous integration

Three checks run on pull requests and one runs on a schedule. They are split so
that each answers exactly one question.

* **The pull-request checks ask "did this change break something?"** They
  install from the frozen `constraints.txt`, so a dependency released overnight
  cannot colour the result. A red check is about the diff.
* **The weekly check asks "has anything upstream moved?"** It deliberately
  installs *without* the constraints. It is about whether a future upgrade is
  safe, not about anything being broken now.

The two used to be a single signal, which is why the separation exists: in
August 2026 a numba release broke the chain report on every open pull request,
and it took a full reproduction to show that the failing PR was not at fault.

### On a pull request

| check | runs on | red means |
|---|---|---|
| `report` | every PR | the diff introduces a new error in the processing chain |
| `pylint` | every PR | score below 8, new errors in changed files, or the coffea versions disagree |
| `resolve` | PRs touching `condor/requirements.txt` | the Condor pins no longer install on CMSSW's Python 3.9 |

**`report`** runs the processor over a small committed fixture twice, once with
`main` and once with the PR, and writes a current-vs-PR diff to the run's **job
summary** (Actions tab, then the run). Open that summary first: it names what
moved — a hist collection that lost a definition, a cut that stopped evaluating,
a channel whose cutflow changed. It is advisory and does not block a merge.

If the summary reports that the chain crashed on **both** sides, the PR is not
at fault: either the environment broke or `main` is broken. The constraints
files make that rare.

**`pylint`** also runs `tests/check_worker_image.py`, which fails when
`requirements.txt`, `constraints.txt` and the apptainer image tag in
`sidm/tools/scaleout.py` disagree about the coffea version. See the upgrade
how-to below.

**`resolve`** appears only on pull requests that touch `condor/requirements.txt`;
its absence elsewhere is normal rather than a skipped check. Condor workers
build their venv from CMSSW, whose Python is 3.9, while everything else here is
3.11 — so a pin that needs Python 3.10 or newer installs cleanly everywhere it
gets tested and then fails on every worker.

### The weekly check

`scheduled-env-check` runs on Tuesdays and can also be triggered by hand from
the Actions tab. It installs the stack **unfrozen** and runs the same chain on
`main`.

**Green** means a fresh, unpinned resolve still works, so regenerating
`constraints.txt` would be safe right now. There is nothing to do, but this is
the signal to consult *before* upgrading anything.

**Red** means something released upstream that breaks the stack. Nothing is
broken for you: the CI, your venv and Condor jobs all install from the frozen
files and are unaffected. It is a warning about the future, not an incident.

1. Do not regenerate `constraints.txt` for now.
2. To find the culprit, open the failing run's `pip freeze` step, which records
   exactly what PyPI served that day, and compare it against the last green run.
3. If the upgrade is needed anyway, pin around the offending package in
   `requirements.txt`, as was done for `numba` and `llvmlite`.

A permanently red weekly check is worse than no check, because it stops being
read. If it stays red, either pin past the problem or decide consciously to live
with it.

Failures show up as a red run in the Actions tab and an email from GitHub to
whoever last committed to the workflow file. No issue is opened, because Issues
are disabled on this repository.

### What the checks do not cover

* **Physics correctness.** `report` verifies that the chain executes and that
  nothing regressed against `main`. A wrong cut or a mis-scaled histogram runs
  cleanly and passes.
* **Dask workers.** No job starts an apptainer worker. `check_worker_image.py`
  verifies that the versions agree, not that workers run.
* **coffea-casa**, where `make_dask_client` installs SIDM from an unpinned git
  branch.

## Miscellaneous how-tos

### How to change or upgrade a dependency

`requirements.txt` says what the analysis needs. `constraints.txt` pins every
package, transitive ones included, and is generated. Condor workers have their
own pair under `condor/`, because they run CMSSW's Python 3.9 and resolve to a
different set. Edit the requirements file; never hand-edit a constraints file.

1. Check that the last weekly `scheduled-env-check` run was green. Red means a
   fresh resolve is currently broken, so this is not the moment to upgrade.
2. Edit `requirements.txt`.
3. Regenerate `constraints.txt` using the recipe in its own header.
4. If coffea moved, update `_DEFAULT_LPC_IMAGE` in `sidm/tools/scaleout.py` to
   the matching tag, and confirm that tag exists under
   `/cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/`. Dask workers
   take coffea from that image rather than from these files, so otherwise the
   client and the workers end up on different versions and fail at runtime. The
   `pylint` check fails if you forget this step.
5. Rebuild the venv (step 3 above) and run the chain locally, or let the pull
   request's `report` check do it for you.
6. Condor needs nothing extra: every job builds its venv fresh from the files it
   is sent.

`pipreqs` is useful for spotting a package you forgot to declare —
`pipreqs . --force` derives requirements from `import` statements — but do not
use it to regenerate `requirements.txt` wholesale. It drops runtime
dependencies that are never imported directly (`bokeh`, loaded by `distributed`
for the Dask dashboard) along with the scale-out packages installed by hand in
step 3 (`lpcjobqueue`, `htcondor`, `distributed`) and the pins added to work
around upstream breakage (`numba`, `llvmlite`), and it leaves `constraints.txt`
stale.

### How to use DASK

**On coffea-casa** (connecting to a pre-existing scheduler):
1. `from sidm.tools import scaleout` in your notebook.
2. Edit `dependencies` in `sidm/tools/scaleout.py` to point at your fork/branch:
   ```python
   dependencies = ["git+https://github.com/YOURUSERNAME/SIDM.git@YOURBRANCH"]
   ```
3. Instantiate the client:
   ```python
   client = scaleout.make_dask_client("tls://localhost:8786")
   ```
4. Pass it to the Runner: `executor=processor.DaskExecutor(client=client)`.

**On LPC** (spinning up an HTCondor-backed cluster from the notebook itself), see [section 7 above](#7-scaling-jobs-from-a-notebook-with-dask--htcondor) — `scaleout.make_lpc_client(...)` returns a ready `(cluster, client)` pair and uploads your local `sidm/` code to workers automatically, so no edit-the-source step is needed. What travels is the importable tree only (`tools`, `definitions`, `scripts`, `configs`, `data`); see [section 7](#7-scaling-jobs-from-a-notebook-with-dask--htcondor).