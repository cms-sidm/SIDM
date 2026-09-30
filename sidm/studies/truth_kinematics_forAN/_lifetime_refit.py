"""Refit the dark-photon proper lifetime from the full-statistics truthkin output.

Loads the 180 per-sample merged coffea files from EOS, adapts them to the
lifetime_analysis API (channel="genOnly": generator truth, no PV filter, unlike
the original baseNoLj_noTrigger run), reruns the grid fits, and prints the
headline comparison numbers. Run from the repo root.
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lifetime_study")))
import glob
import subprocess
import numpy as np
import coffea.util

EOS_DIR = "/store/group/lpcmetx/SIDM/coffea_outputs/murtazas/truth_kinematics_forAN"
CACHE = os.environ.get(
    "TRUTHKIN_CACHE",
    os.path.expanduser("~/nobackup/truthkin_cache"),
)

def load_truthkin(cache=CACHE, eos_dir=EOS_DIR):
    """Download (once) and load all per-sample merged files into one dict."""
    os.makedirs(cache, exist_ok=True)
    names = subprocess.run(["xrdfs", "cmseos.fnal.gov", "ls", eos_dir],
                           capture_output=True, text=True).stdout.split()
    coffeas = [n for n in names if n.endswith(".coffea")]
    out = {}
    for n in coffeas:
        base = os.path.basename(n)
        local = os.path.join(cache, base)
        if not os.path.exists(local):
            subprocess.run(["xrdcp", "-sf", f"root://cmseos.fnal.gov/{n}", local], check=True)
        top = coffea.util.load(local)
        for s, d in top["out"].items():
            out[s] = d
    return out

if __name__ == "__main__":
    import lifetime_analysis as la
    out = load_truthkin()
    print(f"loaded {len(out)} samples", flush=True)
    rows, groups, bg_mean = la.compute_grid(out, channel="genOnly")
    import numpy as np
    mean_r, acc_r, rmaxes = [], [], []
    faithful = []
    for smp, r in rows.items():
        nom = r["nominal"]
        if not nom:
            continue
        m, a = r["mean"] / nom, r["acceptance"] / nom if r["acceptance"] else np.nan
        mean_r.append(m); acc_r.append(a)
        if nom <= 1.0:
            faithful.append(m)
        if np.isfinite(r.get("Rmax", np.nan)):
            rmaxes.append(r["Rmax"])
    acc_r = [x for x in acc_r if np.isfinite(x)]
    print(f"faithful regime (ctau<=1cm): N={len(faithful)}, median mean/nominal = {np.median(faithful):.4f}")
    print(f"acceptance-corrected: N={len(acc_r)}, median = {np.median(acc_r):.4f}, "
          f"within +-5%: {sum(1 for x in acc_r if abs(x-1)<0.05)}/{len(acc_r)}")
    print(f"median fitted R_max = {np.median(rmaxes):.0f} cm  (N={len(rmaxes)})")
    out_npy = os.path.join(CACHE, "refit_grid.npy")
    np.save(out_npy, {"rows": rows, "bg_mean": bg_mean}, allow_pickle=True)
    print(f"saved {out_npy}")
    print("done")
