"""Trigger efficiency of signal events against the generated decay radius of the dark photons
that decay to muons.

usage: python trigger_vs_lxy.py OUT_CSV FILELIST [FILELIST ...]

Each FILELIST holds one ntuple URL per line (the lists condor/make_job_args.py writes); its
basename is <sample>_<chunk>.txt. Only the GenPart and HLT branches are read, with no event
selection, so the denominator is every stored event (all passed the generator filter).

- trigger: the OR of the four DoubleL2Mu NoVtx paths of the "pass triggers" cut.
- Lxy: transverse distance between a dark photon's production vertex and the production vertex
  of its muon daughters. A 2Mu2E event has one muon-decaying dark photon, a 4Mu event two.
- OUT_CSV: one row per sample and (lxy_far_bin, lxy_near_bin), the bins of the larger and the
  smaller Lxy of the event's muon-decaying dark photons (lxy_near_bin is "-" in 2Mu2E), with the
  number of events and the number that pass the trigger. Bin edges 0, 10, 100, 300 cm.
"""
import collections
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor

import awkward as ak
import numpy as np
import pandas as pd
import uproot

HLT = ["HLT_DoubleL2Mu23NoVtx_2Cha", "HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed",
       "HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4", "HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4"]
BRANCHES = ["GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_vx", "GenPart_vy"] + HLT
EDGES = [0, 10, 100, 300]
LABELS = ["[0, 10)", "[10, 100)", "[100, 300)", ">= 300"]
N_THREADS = 8


def read(job):
    """(sample, larger Lxy, smaller Lxy, trigger decision) of one file."""
    sample, url = job
    for attempt in (0, 1):
        try:
            a = uproot.open(url)["Events"].arrays(BRANCHES)
            break
        except OSError:
            if attempt:
                raise
    pid, mom = a.GenPart_pdgId, a.GenPart_genPartIdxMother
    mother = ak.where(mom >= 0, mom, 0)
    mu = (abs(pid) == 13) & (mom >= 0) & (pid[mother] == 32)
    lxy = np.hypot(a.GenPart_vx - a.GenPart_vx[mother], a.GenPart_vy - a.GenPart_vy[mother])[mu]
    trig = a[HLT[0]]
    for path in HLT[1:]:
        trig = trig | a[path]
    return (sample, ak.to_numpy(ak.max(lxy, axis=1)), ak.to_numpy(ak.min(lxy, axis=1)),
            ak.to_numpy(trig))


def main(out_csv, filelists):
    jobs = []
    for fl in filelists:
        sample = re.sub(r"_\d+\.txt$", "", os.path.basename(fl))
        with open(fl, encoding="utf-8") as fh:
            jobs += [(sample, line.strip()) for line in fh if line.strip()]
    with ThreadPoolExecutor(N_THREADS) as pool:
        results = list(pool.map(read, jobs))
    n_all, n_trig = collections.Counter(), collections.Counter()
    for sample, far, near, trig in results:
        two = sample.startswith("4Mu")
        i_far = np.digitize(far, EDGES) - 1
        i_near = np.digitize(near, EDGES) - 1
        for i, far_label in enumerate(LABELS):
            for j, near_label in (enumerate(LABELS) if two else [(0, "-")]):
                sel = i_far == i
                if two:
                    sel = sel & (i_near == j)
                key = (sample, far_label, near_label)
                n_all[key] += int(sel.sum())
                n_trig[key] += int((sel & trig).sum())
    rows = []
    for (sample, far, near), n in sorted(n_all.items()):
        if n == 0:
            continue
        m = re.match(r"(2Mu2E|4Mu)_(\d+)GeV_([\dp]+)GeV_([\dp]+)mm", sample)
        rows.append(dict(sample=sample, final_state=m.group(1), m_med_GeV=float(m.group(2)),
                         m_Zd_GeV=float(m.group(3).replace("p", ".")), lxy_far_bin=far,
                         lxy_near_bin=near, n_events=n, n_trigger=n_trig[(sample, far, near)],
                         efficiency=n_trig[(sample, far, near)] / n))
    pd.DataFrame(rows).to_csv(out_csv, index=False, float_format="%.5g")
    print(f"{len(jobs)} files, {sum(n_all.values())} events, written to {out_csv}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(f"usage: python {sys.argv[0]} OUT_CSV FILELIST [FILELIST ...]")
    main(sys.argv[1], sys.argv[2:])
