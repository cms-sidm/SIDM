"""Generate the full-grid lifetime coffea via Dask-over-Condor (lpcjobqueue).

One-time generation helper (not part of the committed notebooks): the notebooks
load the resulting lifetime_fullgrid.coffea. Run from the repo root:
    python sidm/studies/lifetime_study/_fullgrid_run.py
"""
import sys
import os
# Force the working-tree sidm (with the new lifetime infra) ahead of the stale
# non-editable copy in site-packages: `python script.py` puts the script dir on
# sys.path[0], not the repo root, so insert the repo root explicitly.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
import time
import subprocess
import yaml
import coffea.util
from coffea import processor
from sidm.tools import utilities, sidm_processor, llpnanoaodschema, scaleout
from sidm.tools.metadata import write_run_metadata

t0 = time.time()
def log(msg):
    print(f"[{time.time()-t0:6.0f}s] {msg}", flush=True)

# Full v10 grid: every 4Mu and 2Mu2E sample (BS mass x Dp mass x ctau)
s4 = list(yaml.safe_load(open("sidm/configs/ntuples/signal_4mu_v10.yaml"))["llpNanoAOD_v2"]["samples"].keys())
s2 = list(yaml.safe_load(open("sidm/configs/ntuples/signal_2mu2e_v10.yaml"))["llpNanoAOD_v2"]["samples"].keys())
fs = utilities.make_fileset(s4, "llpNanoAOD_v2", max_files=1,
                            location_cfg="signal_4mu_v10.yaml", replace_xcache=True)
fs = utilities.make_fileset(s2, "llpNanoAOD_v2", max_files=1,
                            location_cfg="signal_2mu2e_v10.yaml", fileset=fs, replace_xcache=True)
log(f"fileset built: {len(fs)} samples ({len(s4)} 4Mu + {len(s2)} 2Mu2E)")

cluster, client = scaleout.make_lpc_client(min_workers=10, max_workers=40, memory="4GB", disk="4GB")
log(f"dashboard: {cluster.dashboard_link}")
client.wait_for_workers(1, timeout=600)
log(f">=1 worker connected; launching full-grid run")

runner = processor.Runner(
    executor=processor.DaskExecutor(client=client, status=False),
    schema=llpnanoaodschema.LLPNanoAODSchema,
    skipbadfiles=True, chunksize=50_000,
)
p = sidm_processor.SidmProcessor(["baseNoLj_noTrigger"], ["genA_lifetime"], unweighted_hist=True)
out = runner.run(fs, treename="Events", processor_instance=p)["out"]
local_path = "sidm/studies/lifetime_study/lifetime_fullgrid.coffea"
coffea.util.save(out, local_path)
meta_path = write_run_metadata(
    local_path, fileset=fs,
    selections=["baseNoLj_noTrigger"], hist_collections=["genA_lifetime"],
    schema="LLPNanoAODSchema", chunksize=50_000, unweighted_hist=True,
)
log(f"SAVED local {local_path} (+ {os.path.basename(meta_path)}): {len(out)}/{len(fs)} samples")
cluster.close(); client.close()

# durable copy of the coffea + sidecar to the shared lpcmetx EOS area (re-readable later)
eos_dir = "/store/group/lpcmetx/SIDM/coffea_outputs/murtazas/lifetime_study"
subprocess.run(["xrdfs", "cmseos.fnal.gov", "mkdir", "-p", eos_dir], check=False)
for lp in (local_path, meta_path):
    eos_path = f"root://cmseos.fnal.gov/{eos_dir}/{os.path.basename(lp)}"
    r = subprocess.run(["xrdcp", "-f", lp, eos_path], capture_output=True, text=True)
    log(f"xrdcp -> EOS rc={r.returncode}  {eos_path}  {r.stderr.strip()[:160]}")
log("done")
