"""Generate the reco-migration coffea output over the full v10 signal grid.

Documents the reconstruction-migration configuration and serves for small Dask
reruns; the canonical production runs through the condor pipeline (see
README.md). Measures how dark-photon decays populate the four
reconstructed-object categories (GED electron, photon, PF muon, DSA muon) as a
function of Lxy, via the lepton_genA_base nearGenA histograms. Outputs go to a
dask_reruns subdirectory so they never mix with the canonical merged files. Channels keep
the LJ-source object definitions and PV filter but deliberately drop the HLT,
so the migration story is trigger-independent; the NoLjsource channel gives
the same maps before object quality cuts. Run from the repo root:
    python sidm/studies/truth_kinematics_forAN/_recomigration_run.py
Smoke test:
    python sidm/studies/truth_kinematics_forAN/_recomigration_run.py \
        --local --max-files 1 --limit 1 --tag smoke --no-eos
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
import time
import argparse
import subprocess
import yaml
import coffea.util
from coffea import processor
from sidm.tools import utilities, sidm_processor, llpnanoaodschema, scaleout
from sidm.tools.metadata import write_run_metadata

CHANNELS = ["baseNoLj_noTrigger", "baseNoLjNoLjsource_noTrigger"]
COLLECTIONS = ["genA_base", "lepton_genA_base"]
EOS_DIR = "/store/group/lpcmetx/SIDM/coffea_outputs/murtazas/truth_kinematics_forAN/dask_reruns"

parser = argparse.ArgumentParser()
parser.add_argument("--max-files", type=int, default=-1)
parser.add_argument("--limit", type=int, default=0, help="take only the first N samples per channel yaml")
parser.add_argument("--local", action="store_true", help="iterative executor, no Dask")
parser.add_argument("--no-eos", action="store_true", help="skip the EOS upload")
parser.add_argument("--tag", default="fullgrid", help="output filename tag")
args = parser.parse_args()

t0 = time.time()
def log(msg):
    print(f"[{time.time()-t0:6.0f}s] {msg}", flush=True)

s4 = list(yaml.safe_load(open("sidm/configs/ntuples/signal_4mu_v10.yaml"))["llpNanoAOD_v2"]["samples"].keys())
s2 = list(yaml.safe_load(open("sidm/configs/ntuples/signal_2mu2e_v10.yaml"))["llpNanoAOD_v2"]["samples"].keys())
if args.limit:
    s4, s2 = s4[:args.limit], s2[:args.limit]
fs = utilities.make_fileset(s4, "llpNanoAOD_v2", max_files=args.max_files,
                            location_cfg="signal_4mu_v10.yaml", replace_xcache=True)
fs = utilities.make_fileset(s2, "llpNanoAOD_v2", max_files=args.max_files,
                            location_cfg="signal_2mu2e_v10.yaml", fileset=fs, replace_xcache=True)
log(f"fileset built: {len(fs)} samples ({len(s4)} 4Mu + {len(s2)} 2Mu2E)")

if args.local:
    executor = processor.IterativeExecutor()
else:
    cluster, client = scaleout.make_lpc_client(min_workers=10, max_workers=100, memory="4GB", disk="4GB")
    log(f"dashboard: {cluster.dashboard_link}")
    client.wait_for_workers(1, timeout=600)
    log(">=1 worker connected; launching run")
    executor = processor.DaskExecutor(client=client, status=False)

runner = processor.Runner(
    executor=executor,
    schema=llpnanoaodschema.LLPNanoAODSchema,
    skipbadfiles=True, chunksize=50_000,
)
p = sidm_processor.SidmProcessor(CHANNELS, COLLECTIONS, unweighted_hist=True)
out = runner.run(fs, treename="Events", processor_instance=p)["out"]
wrapped = {"out": out}

local_path = f"sidm/studies/truth_kinematics_forAN/truthkin_recomigration_{args.tag}.coffea"
coffea.util.save(wrapped, local_path)
meta_path = write_run_metadata(
    local_path, fileset=fs,
    selections=CHANNELS, hist_collections=COLLECTIONS,
    schema="LLPNanoAODSchema", chunksize=50_000, unweighted_hist=True,
)
log(f"SAVED local {local_path} (+ {os.path.basename(meta_path)}): {len(out)}/{len(fs)} samples")
if not args.local:
    cluster.close(); client.close()

if not args.no_eos:
    subprocess.run(["xrdfs", "cmseos.fnal.gov", "mkdir", "-p", EOS_DIR], check=False)
    for lp in (local_path, meta_path):
        eos_path = f"root://cmseos.fnal.gov/{EOS_DIR}/{os.path.basename(lp)}"
        r = subprocess.run(["xrdcp", "-f", lp, eos_path], capture_output=True, text=True)
        log(f"xrdcp -> EOS rc={r.returncode}  {eos_path}  {r.stderr.strip()[:160]}")
log("done")
