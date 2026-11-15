"""Halo-filter study run: the 1-LJ control region with the MET filters recorded one by one, on
data and background MC, and the signal channels with and without "pass flags".

    D=sidm/studies/vertexChi2/vtx_outcome_lxy_dr
    python $D/run_halo_filter.py OUT_DIR {data,mc,signal} [...]

Runs with Dask on the LPC condor pool (sidm.tools.scaleout.make_lpc_client), one sample at a
time, and writes OUT_DIR/<run tag>/<sample>.coffea with its .meta.yaml sidecar; a sample whose
output exists is skipped, so an interrupted run can be restarted. Channels are those of
configs/selections_vertexChi2_halo.yaml.

  data    68 DoubleMuon 2018 skims -> halo_data_2018: cr_1egmOr1mu_lj and
          cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu, each as is (with "pass flags"), _noFlags
          and _failHalo; histogram collection halo_filter
  mc      the eight background samples -> halo_mc_2018, same channels and histograms
  signal  24 samples, the first five files of each -> halo_signal_2018: 4mu / 2mu2e and their
          _flags variants; histogram collection met_filter_mask; raw-count histograms

halo_filter_counts.py reads the three directories. Needs X509_USER_PROXY at a proxy the condor
workers can read."""
import os
import pathlib
import re
import sys
import time

import coffea.util
import yaml
from coffea import processor

SIDM = str(pathlib.Path(__file__).resolve().parents[4])   # repository root
sys.path.insert(0, SIDM)
# pylint: disable=wrong-import-position
from sidm.tools import utilities, scaleout, sidm_processor, llpnanoaodschema
from sidm.tools.metadata import write_run_metadata
# pylint: enable=wrong-import-position

SELECTIONS_CFG = "configs/selections_vertexChi2_halo.yaml"
CR = ["cr_1egmOr1mu_lj", "cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu"]
CR_CHANNELS = [c + s for c in CR for s in ("", "_noFlags", "_failHalo")]
SKIP_JSON = f"{SIDM}/sidm/studies/vertexChi2/signal_skip.json"

DATA = (["DoubleMuon_2018A_%d" % i for i in range(3)]
        + ["DoubleMuon_2018B_%d" % i for i in range(14)]
        + ["DoubleMuon_2018C"]
        + ["DoubleMuon_2018D_%d" % i for i in range(50)])
BGS = ["QCD_Pt170To300", "QCD_Pt300To470", "DYJetsToLL_M50", "QCD_Pt80To120", "QCD_Pt470To600",
       "QCD_Pt120To170", "QCD_Pt600To800", "QCD_Pt800To1000"]


def signal_samples():
    """2Mu2E and 4Mu at 100/500/1000 GeV, Zd mass 0.25 and 5 GeV, shortest and longest ctau."""
    out = []
    for fs, cfg in (("2Mu2E", "signal_2mu2e_v10.yaml"), ("4Mu", "signal_4mu_v10.yaml")):
        with open(f"{SIDM}/sidm/configs/ntuples/{cfg}", encoding="utf-8") as fh:
            names = list(yaml.safe_load(fh)["llpNanoAOD_v2"]["samples"])
        for m in ("100GeV", "500GeV", "1000GeV"):
            for ma in ("0p25GeV", "5p0GeV"):
                pts = [n for n in names if re.fullmatch(f"{fs}_{m}_{ma}_[0-9p]+mm", n)]
                ct = lambda n: float(n.split("_")[-1][:-2].replace("p", "."))
                pts.sort(key=ct)
                out += [(pts[0], cfg), (pts[-1], cfg)]
    return out


def save(output, out_dir, x, run_tag, fileset, channels, collections, unweighted=None):
    """Write OUT_DIR/<run_tag>/<x>.coffea and its sidecar."""
    d = f"{out_dir}/{run_tag}"
    os.makedirs(d, exist_ok=True)
    coffea.util.save(output, f"{d}/{x}.coffea")
    write_run_metadata(f"{d}/{x}.coffea", fileset=fileset, selections=channels,
                       selections_cfg=SELECTIONS_CFG, hist_collections=collections,
                       schema="LLPNanoAODSchema", chunksize=10000,
                       unweighted_hist=unweighted, extra={"run_tag": run_tag})


def main(out_dir, parts):
    cluster, client = scaleout.make_lpc_client(min_workers=10, max_workers=40, memory="4GB",
                                               disk="4GB")
    client.wait_for_workers(1, timeout=6 * 3600)
    print("first worker connected:", cluster, flush=True)
    runner = processor.Runner(executor=processor.DaskExecutor(client=client, status=False),
                              schema=llpnanoaodschema.LLPNanoAODSchema, skipbadfiles=True,
                              chunksize=10000)
    cr_cols = ["halo_filter"]
    p_cr = sidm_processor.SidmProcessor(CR_CHANNELS, cr_cols, selections_cfg=SELECTIONS_CFG)
    t0 = time.time()
    try:
        if "mc" in parts:
            for x in BGS:
                if os.path.exists(f"{out_dir}/halo_mc_2018/{x}.coffea"):
                    continue
                fs = utilities.make_fileset([x], "skimmed_llpNanoAOD_v2", max_files=-1,
                                            location_cfg="backgrounds.yaml", replace_xcache=True)
                out = runner.run(fs, treename="Events", processor_instance=p_cr)
                save(out, out_dir, x, "halo_mc_2018", fs, CR_CHANNELS, cr_cols)
                print(f"[{time.time() - t0:7.0f}s] mc {x} done", flush=True)
        if "data" in parts:
            for x in DATA:
                if os.path.exists(f"{out_dir}/halo_data_2018/{x}.coffea"):
                    continue
                fs = utilities.make_fileset([x], "llpNanoAOD_v2", max_files=-1,
                                            location_cfg="data_skimmed.yaml", replace_xcache=True)
                out = runner.run(fs, treename="Events", processor_instance=p_cr)
                save(out, out_dir, x, "halo_data_2018", fs, CR_CHANNELS, cr_cols)
                print(f"[{time.time() - t0:7.0f}s] data {x} done", flush=True)
        if "signal" in parts:
            for x, cfg in signal_samples():
                if os.path.exists(f"{out_dir}/halo_signal_2018/{x}.coffea"):
                    continue
                chans = ["4mu", "4mu_flags"] if x.startswith("4Mu") else ["2mu2e", "2mu2e_flags"]
                p_sig = sidm_processor.SidmProcessor(chans, ["met_filter_mask"],
                                                     selections_cfg=SELECTIONS_CFG,
                                                     unweighted_hist=True)
                fs = utilities.make_fileset([x], "llpNanoAOD_v2", max_files=5, location_cfg=cfg,
                                            replace_xcache=True, census_skip=SKIP_JSON)
                out = runner.run(fs, treename="Events", processor_instance=p_sig)
                save(out, out_dir, x, "halo_signal_2018", fs, chans, ["met_filter_mask"],
                     unweighted=True)
                print(f"[{time.time() - t0:7.0f}s] signal {x} done", flush=True)
    finally:
        client.close()
        cluster.close()
    print("ALL DONE", time.time() - t0, flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 3 or not set(sys.argv[2:]) <= {"data", "mc", "signal"}:
        sys.exit(f"usage: python {sys.argv[0]} OUT_DIR {{data,mc,signal}} [...]")
    main(sys.argv[1], sys.argv[2:])
