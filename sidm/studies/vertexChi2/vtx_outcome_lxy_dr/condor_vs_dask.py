"""Condor-merged 4Mu_100GeV_0p25GeV_0p02mm against the Dask-produced one in dask_check/, flow
included, for the three new histograms and mu_lj_vtx_chi2, both channels.

usage: python condor_vs_dask.py [REPORT_FILE]   (default condor_vs_dask.txt in the working directory)

Condor chunks are unweighted (run_job.sh passes --unweighted-hist); the Dask output carries
genWeight * lumi*xs/sumw. So the test is: Dask = c * condor in every bin, with c the
lumi*xs/sumw of the merged metadata, and Dask variances = c^2 * condor counts (true only if
every event has the same genWeight). The report is written to REPORT_FILE (and echoed), since
the cross-section lookup prints a notice of its own to stdout."""
import pathlib
import sys
import coffea.util
import numpy as np
from sidm.tools import utilities
R = str(pathlib.Path(__file__).resolve().parents[4] / "RunOutputFiles" / "vtxLxyDr_signal_2018")
N = "4Mu_100GeV_0p25GeV_0p02mm"
REPORT = sys.argv[1] if len(sys.argv) > 1 else "condor_vs_dask.txt"
d = coffea.util.load(f"{R}/dask_check/{N}.coffea")["out"][N]
c = coffea.util.load(f"{R}/{N}.coffea")["out"][N]
md = c["metadata"]
scale = utilities.get_lumixs_weight(N, str(next(iter(md["year"]))), md["scaled_sum_weights"])
lines = [f"condor unweighted_hist={md.get('unweighted_hist')}; n_evts dask {d['metadata']['n_evts']} "
         f"condor {md['n_evts']}; sumw dask {d['metadata']['scaled_sum_weights']:.6g} condor "
         f"{md['scaled_sum_weights']:.6g}; c = lumi*xs/sumw = {scale:.8g}"]
worst = 0.0
fits = []
for hn in ("mu_lj_vtxChi2_muonN_genA_lxy", "mu_lj_vtxChi2_muonN_mumu_dR",
           "mu_lj_vtxChi2_genA_lxy_mumu_dR", "mu_lj_vtx_chi2"):
    for ch in ("base_ljDisplacementIso_cosAlpha", "base_ljDisplacementIso_cosAlpha_ljSpread"):
        hd, hc = d["hists"][hn][{"channel": ch}], c["hists"][hn][{"channel": ch}]
        vd, wd = hd.values(flow=True), hd.variances(flow=True)
        vc, wc = hc.values(flow=True), hc.variances(flow=True)
        assert vd.shape == vc.shape
        nz = vc > 0
        c_fit = float((vd[nz] * vc[nz]).sum() / (vc[nz] ** 2).sum())      # least squares, vd = c*vc
        fits.append(c_fit)
        rel_v = np.abs(vd - scale * vc).max() / vd.max()
        rel_w = np.abs(wd - scale ** 2 * vc).max() / wd.max()
        same_support = np.array_equal(vd > 0, vc > 0)
        counts_int = np.allclose(vc, np.round(vc)) and np.allclose(wc, vc)
        worst = max(worst, rel_v, rel_w)
        lines.append(f"{hn:34s} {ch[-12:]:12s} bins {vd.size:5d} nonzero {int((vc > 0).sum()):4d} "
                     f"raw total {vc.sum():9.0f}  max|dask - c*condor|/max {rel_v:.2e}  "
                     f"max|var_dask - c^2*condor|/max {rel_w:.2e}  same nonzero bins {same_support}  "
                     f"condor = integer counts {counts_int}")
lines.append(f"fitted c over the 8 (hist, channel) pairs: min {min(fits):.10g} max {max(fits):.10g}; "
             f"lumi*xs/sumw {scale:.10g}; max |c_fit/c - 1| {max(abs(f / scale - 1) for f in fits):.3g}")
lines.append(f"worst relative difference: {worst:.3g}")
with open(REPORT, "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
