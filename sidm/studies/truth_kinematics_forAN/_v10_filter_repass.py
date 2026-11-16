"""Tabulate the gen-filter re-pass counts of the 180 v10 signal samples.

Reads the merged per-sample truth_kinematics_forAN outputs, copied from EOS
    /store/group/lpcmetx/SIDM/coffea_outputs/murtazas/truth_kinematics_forAN/<sample>.coffea
into a local ./coffea/ directory, and writes ./v10_filter_repass.json: for each
sample and each emulation channel (genFilterEmulation,
genFilterEmulation_isHardProcess) the first two cutflow rows, the raw event
counts before (n_all) and after (n_pass) the "pass gen filter" cut, and their
ratio (repass). It then prints the re-pass range and the median per mXX.

The emulation channels apply the central-production cut string (pT > 5 GeV);
the v10 samples were produced with the same string at pT > 4 GeV, so repass is
the fraction of stored events without a 4-5 GeV lepton deciding the count.

Run from a directory holding ./coffea/:
    xrdcp -r root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/coffea_outputs/murtazas/truth_kinematics_forAN/ coffea
    python _v10_filter_repass.py
"""
import glob
import json

import numpy as np
from coffea.util import load

CHANNELS = ("genFilterEmulation", "genFilterEmulation_isHardProcess")
MASSES = ["100GeV", "150GeV", "200GeV", "500GeV", "800GeV", "1000GeV"]


def main():
    """Write v10_filter_repass.json from ./coffea/*.coffea and print its summary."""
    counts = {}
    for path in sorted(glob.glob("coffea/*.coffea")):
        output = load(path)
        if isinstance(output, dict) and "out" in output:
            output = output["out"]
        for sample, sample_out in output.items():
            if not isinstance(sample_out, dict) or "cutflow" not in sample_out:
                continue
            entry = {}
            for channel in CHANNELS:
                rows = sample_out["cutflow"][channel].rows
                names = list(rows)
                n_all, n_pass = rows[names[0]]["raw"], rows[names[1]]["raw"]
                entry[channel] = {"cuts": names[:2], "n_all": int(n_all), "n_pass": int(n_pass),
                                  "repass": float(n_pass / n_all)}
            counts[sample] = entry
    with open("v10_filter_repass.json", "w", encoding="utf8") as f:
        json.dump(counts, f, indent=1, sort_keys=True)
    repass = np.array([v["genFilterEmulation"]["repass"] for v in counts.values()])
    print(len(counts), "samples; range", round(repass.min(), 3), round(repass.max(), 3),
          "median", round(float(np.median(repass)), 3))
    for mass in MASSES:
        sel = [v["genFilterEmulation"]["repass"] for s, v in counts.items() if f"_{mass}_" in s]
        print(mass, round(float(np.median(sel)), 3), len(sel))


if __name__ == "__main__":
    main()
