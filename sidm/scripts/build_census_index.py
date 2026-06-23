#!/usr/bin/env python3
"""Generate sidm/configs/census/README.md -- the index of published file censuses -- from the
committed *.census.summary.json files. Re-run after publishing a new census. Run from repo root."""
import glob
import json
import os

REDIR_DIR = "/store/group/lpcmetx/SIDM/census"
rows = []
for f in sorted(glob.glob("sidm/configs/census/*.census.summary.json")):
    s = json.load(open(f))
    m = s["meta"]
    # the per-sample rollup always carries the true version block; meta.version can be a
    # blanket "ALL" for older/merge-built manifests, so prefer the sample value.
    samp_ver = s["samples"][0].get("version") if s.get("samples") else None
    rows.append({
        "run_id": s["run_id"],
        "summary": os.path.basename(f),
        "yaml": os.path.basename(m.get("source_yaml", "?")),
        "version": samp_ver or m.get("version", "?"),
        "n_files": s["n_files"],
        "n_samples": len(s["samples"]),
        "n_dropped": len(s["dropped"]),
        "n_anomaly": sum(x.get("n_runs_anomaly", 0) for x in s["samples"]),
        "complete": m.get("complete"),
        "backend": m.get("backend", "threaded"),
        "date": (m.get("ended_utc") or "")[:10],
    })

lines = []
lines.append("# Published file censuses\n")
lines.append("Each row is a completed **shallow (\"light\") census** of one location YAML. The")
lines.append("committed `*.census.summary.json` in this directory holds the per-sample rollup, the")
lines.append("dropped-file list (with reasons), and the provenance `meta` (source-YAML sha256,")
lines.append("`complete` + `n_enumerated`/`n_probed`, backend). The **full per-file manifest** and the")
lines.append("**cleaned per-sample filelists** that jobs run over live on EOS under")
lines.append(f"`{REDIR_DIR}/<run-id>/` as `manifest.json` + `cleaned_filelists.tar.gz`.\n")
lines.append("- How to run a census: [`condor/CENSUS_README.md`](../../../condor/CENSUS_README.md)")
lines.append("- How to read one: [`sidm/test_notebooks/file_census.ipynb`]"
             "(../../test_notebooks/file_census.ipynb)\n")
lines.append("| run-id | source YAML | version | files | samples | dropped | runs-anomaly | complete | backend | date |")
lines.append("|---|---|---|---:|---:|---:|---:|:---:|---|---|")
for r in rows:
    chk = "yes" if r["complete"] else "**NO**"
    lines.append(f"| `{r['run_id']}` | {r['yaml']} | {r['version']} | {r['n_files']:,} | "
                 f"{r['n_samples']} | {r['n_dropped']} | {r['n_anomaly']} | {chk} | {r['backend']} | {r['date']} |")
tot_f = sum(r["n_files"] for r in rows)
tot_d = sum(r["n_dropped"] for r in rows)
tot_a = sum(r["n_anomaly"] for r in rows)
lines.append(f"\n**Total: {tot_f:,} files** across {len(rows)} censuses, {tot_d} dropped "
             f"(bad / unreachable / empty -- excluded from the cleaned filelists), {tot_a} "
             "runs-count anomalies (kept).\n")
lines.append("**runs-anomaly** = files where the Runs-tree counters (genEventCount + genEventSumw) "
             "undercount the Events tree (Events > genEventCount) -- a producer counting convention in "
             "the signal ntuples, NOT corruption. These files are KEPT (they process cleanly); "
             "normalize from the Events genWeight sum, not the Runs Sw.\n")
lines.append("Each run is `complete`: the dask driver and the Condor merge both reconcile probed-")
lines.append("vs-enumerated counts and refuse to publish a partial census as complete. To load any")
lines.append("run, `xrdcp` its `manifest.json` from the EOS path above (see the notebook).")

out = "sidm/configs/census/README.md"
with open(out, "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"wrote {out} indexing {len(rows)} censuses ({tot_f:,} files total)")
