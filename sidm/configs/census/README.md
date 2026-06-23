# Published file censuses

Each row is a completed **shallow ("light") census** of one location YAML. The
committed `*.census.summary.json` in this directory holds the per-sample rollup, the
dropped-file list (with reasons), and the provenance `meta` (source-YAML sha256,
`complete` + `n_enumerated`/`n_probed`, backend). The **full per-file manifest** and the
**cleaned per-sample filelists** that jobs run over live on EOS under
`/store/group/lpcmetx/SIDM/census/<run-id>/` as `manifest.json` + `cleaned_filelists.tar.gz`.

- How to run a census: [`condor/CENSUS_README.md`](../../../condor/CENSUS_README.md)
- How to read one: [`sidm/test_notebooks/file_census.ipynb`](../../test_notebooks/file_census.ipynb)

| run-id | source YAML | version | files | samples | dropped | runs-anomaly | complete | backend | date |
|---|---|---|---:|---:|---:|---:|:---:|---|---|
| `backgrounds_skimmed_light` | backgrounds.yaml | skimmed_llpNanoAOD_v2 | 56,296 | 18 | 256 | 0 | yes | condor (113 shards) | 2026-06-23 |
| `data_light` | data.yaml | llpNanoAOD_v2 | 636 | 1 | 0 | 0 | yes | dask (max_workers=20) | 2026-06-23 |
| `data_skimmed_light` | data_skimmed.yaml | llpNanoAOD_v2 | 127,415 | 68 | 324 | 0 | yes | condor (255 shards) | 2026-06-23 |
| `signal_2mu2e_v10_light` | signal_2mu2e_v10.yaml | llpNanoAOD_v2 | 4,082 | 90 | 5 | 57 | yes | condor (9 shards) | 2026-06-23 |
| `signal_4mu_v10_light` | signal_4mu_v10.yaml | llpNanoAOD_v2 | 4,960 | 90 | 1 | 24 | yes | condor (10 shards) | 2026-06-23 |

**Total: 193,389 files** across 5 censuses, 586 dropped (bad / unreachable / empty -- excluded from the cleaned filelists), 81 runs-count anomalies (kept).

**runs-anomaly** = files where the Runs-tree counters (genEventCount + genEventSumw) undercount the Events tree (Events > genEventCount) -- a producer counting convention in the signal ntuples, NOT corruption. These files are KEPT (they process cleanly); normalize from the Events genWeight sum, not the Runs Sw.

Each run is `complete`: the dask driver and the Condor merge both reconcile probed-
vs-enumerated counts and refuse to publish a partial census as complete. To load any
run, `xrdcp` its `manifest.json` from the EOS path above (see the notebook).
