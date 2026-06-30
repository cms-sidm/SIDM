# SIDM file census

Open **every** input ROOT file of every sample — including the ones hand-commented-out in the
location YAMLs (which `make_fileset` never sees, because `yaml.safe_load` drops the `# - …root`
lines) — verify each one, and roll the verdicts into a self-describing **manifest** and **cleaned
per-sample filelists**. It tells you which files are good / bad / unreachable / empty, restores
files that were wrongly retired, and gives a denominator-honest event/Σw census.

The tool is `sidm/tools/file_census.py`. The verdict discipline is shared with the campaign
reconciler (`condor/campaign_lib.py`): **a file is condemned only on a positive, deterministic
failure** (a reproducible open/decompress error, or a server "no such file"); a timeout or
redirector hiccup is `inconclusive` and retried, never blacklisted.

## Verdicts

| verdict | meaning |
|---|---|
| `reachable` | opens and reads (the analysis can use it). `events_entries == 0` ⇒ also flagged **empty** |
| `bad` | opens but is corrupt/truncated, or hits a deterministic decompress error — a plain open-check misses these |
| `unreachable` | a positive server "no such file" |
| `inconclusive` | could not be determined (timeout / probe error) — re-probed next run |

**Empty files are valid, not corrupt** (0 events survived the skim). They matter because a Condor
chunk built *entirely* of empties makes coffea raise `TypeError: cannot unpack non-iterable NoneType`
(fixed by the `run_sidm_chunk.py` guard); the census flags them so they can be dropped or spread.

**Reachable files can also carry non-fatal flags** and stay usable — they are *kept*, never
condemned: a **runs-count anomaly** (`Events > genEventCount`, a producer counting convention, not
corruption — normalize from the Events `genWeight` sum, not the `Runs` `Sw`) and **GenPart
self-reference** (see below). Both are surfaced in the report + manifest but never drop a file by
themselves.

## Two depths

- **shallow** (default): open + read the per-file `Runs` tree (`genEventSumw`, `genEventCount`) +
  peek `Events`, plus the GenPart self-reference scan on MC (below) — no event loop. Fast; this is
  the publishable light census. Uses the same xrootd transport dask/Condor workers use, so it
  genuinely checks file access.
- **deep** (`--deep`): additionally run the **full `SidmProcessor`** on each reachable non-empty
  file to certify it processes (per-file `process_status`). Expensive (≈ running the analysis per
  file); for validation runs, **not** the light census. (Self-reference files are processed too,
  not skipped — the `gen.py` guard keeps the processor from hanging on them.)

## GenPart self-reference (the distinctParent hang)

Some skimmed QCD outputs contain corrupt, **self-referential** `GenPart` entries — a particle whose
`genPartIdxMother` points to its own index. coffea computes `distinctParent` by climbing the
same-pdgId mother chain in an unguarded loop, so that self-loop never terminates and the **job
hangs** on a single such event. The census reads `GenPart_genPartIdxMother` on every **MC** file by
default — a cheap extra branch read (≈ +5 %; the file open dominates), `--no-check-genpart` to turn
it off, data has no GenPart and is skipped — and reports the offenders under **GENPART INTEGRITY**
with per-file event counts.

These files are **kept** (`reachable`, not condemned): the corrupt entries are a handful of inert
filler particles, and the real gens and weights are intact. The hang itself is fixed independently
by the cycle-safe `distinctParent` kernel in [`sidm/tools/gen.py`](../sidm/tools/gen.py), installed
process-wide whenever the LLPNanoAOD schema is used, so the processor (and any gen mother-tracking
study) runs over them without hanging. If you must exclude them anyway, `filelists
--drop-genpart-corrupt` adds them to the skip-list.

## Running it

```bash
cd /uscms_data/d3/$USER/SIDM
source sidm_venv/bin/activate

# serial / threaded from the submitter node (shallow probe is I/O-bound -> threads scale well):
python sidm/tools/file_census.py census --location-cfg signal_2mu2e_v10.yaml \
    --version llpNanoAOD_v2 --workers 48 --run-id signal_2mu2e_v10 \
    --out census_signal_2mu2e_v10.json

# dask, on LPC Condor workers (cheap shallow census at larger scale):
python sidm/scripts/census_dask.py --location-cfg backgrounds.yaml \
    --version skimmed_llpNanoAOD_v2 --batch-size 50 --max-workers 40 \
    --out census_backgrounds.json --run-id backgrounds

# Condor (the authoritative, detached run; supports deep) -- FOUR steps:
# 1. chunk: writes condor/census_job_args.txt, condor/filelists_census/*.txt, census_expected.json
python condor/make_census_args.py --location-cfg backgrounds.yaml \
    --version skimmed_llpNanoAOD_v2 --files-per-job 200
# 2. rebuild the code tarball INCLUDING the census worker + campaign_lib, and mkdir the EOS dir:
tar -czf condor/sidm_code.tar.gz sidm \
    condor/run_sidm_chunk.py condor/run_census_chunk.py condor/campaign_lib.py
xrdfs root://cmseos.fnal.gov mkdir -p /store/user/$USER/sidm_census/run_v1
# 3. submit (edit the EOS shard dir + shallow|deep in submit_census.sub to match step 2/4):
ssh lpc 'cd /uscms_data/d3/$USER/SIDM/condor && condor_submit submit_census.sub'
# 4. merge WITH the expected sidecar so a missing shard / partial run is caught, not published:
python condor/merge_census_shards.py --shard-eos-dir /store/user/$USER/sidm_census/run_v1 \
    --expected condor/census_expected.json --out census_backgrounds.json --run-id backgrounds
```

The tarball **must** contain `run_census_chunk.py` + `campaign_lib.py` (step 2) or every job dies on
its first command; the EOS shard dir **must** exist (step 2) or every `xrdcp` fails; and `--expected`
(step 4) is what lets `merge` refuse to publish a partial census (it stamps `meta.complete`).

## Cleaned filelists (what people run over)

```bash
python sidm/tools/file_census.py filelists --manifest census_signal_2mu2e_v10.json \
    --out-dir cleaned_signal_2mu2e_v10
```

Writes `<sample>.txt` of the **good** files (reachable, non-empty — live *or* commented-but-good, so
a wrongly retired file is restored), plus `_dropped.tsv` (with reasons), **`_skip.json`** (the
machine-readable veto for `make_fileset`, see below), and `_census_ref.json` (provenance pointer).
Publish the manifest + filelists to the shared `/store/group/lpcmetx/SIDM/census/<run-id>/` area on
EOS; jobs point at the cleaned lists. Pass `--drop-genpart-corrupt` to also drop the GenPart
self-reference files (default: keep them — the `gen.py` guard handles them).

## Skip flagged files in your jobs (the runner veto)

`filelists` writes `_skip.json` — the per-sample list of flagged files (bad / unreachable / empty,
plus the self-reference files if you passed `--drop-genpart-corrupt`). Point any run at it and those
files are dropped **in addition to** the ones the location YAML comments out (the YAML's manual
commenting is always honored — the veto only removes more):

```python
# notebook / dask:
fileset = utilities.make_fileset(samples, version, location_cfg="backgrounds.yaml",
                                 census_skip="backgrounds_skimmed.skip.json")
```
```bash
# Condor:
python condor/make_job_args.py --location-cfg backgrounds.yaml --version skimmed_llpNanoAOD_v2 \
    --census-skip backgrounds_skimmed.skip.json
```

`census_skip` is a path or a bare name resolved under `sidm/configs/census/`. Matching is by
basename per sample, so it is robust to xcache-vs-cmseos URLs.

## Manifest

Self-describing. `meta` records when/what/where/how/by-whom: `mode`, start/end + duration, user +
host, sidm commit + branch, uproot version, the source YAML **+ its sha256**, version, samples,
`max_live`, redirector, probe config. `files[]` is one row per file (verdict, events, `Runs` info,
`n_genpart_selfref_events`, `process_status`); `samples[]` is the per-sample rollup (counts,
`n_empty`, `n_genpart_anomaly`, `skim_basis`, `genEventSumw_reachable_norm`). The light shallow run
is the artifact published with the PR; the deep run is operational validation.

## Published censuses

The censuses already run for this analysis are indexed in
[`sidm/configs/census/README.md`](../sidm/configs/census/README.md) — one row per run with its
file / sample / dropped counts, backend, completeness, and date. For each run the committed
`sidm/configs/census/<run-id>.census.summary.json` holds the per-sample rollup + the dropped-file
list, and the **full manifest + cleaned filelists** live on EOS at
`/store/group/lpcmetx/SIDM/census/<run-id>/`. The example notebook
[`sidm/test_notebooks/file_census.ipynb`](../sidm/test_notebooks/file_census.ipynb) walks through
loading one and reading its verdicts. Re-generate the index with
`python sidm/scripts/build_census_index.py` after publishing a new run.

## Refreshing the committed census

A census run writes its outputs to wherever you point `--out` / `--out-dir`; it does **not** touch
the committed defaults under `sidm/configs/census/`. To refresh the blessed default after a new run:
copy your `<run-id>.census.summary.json` (and, to update the veto, the `_skip.json` as
`<run-id>.skip.json`) into `sidm/configs/census/`, `xrdcp` the full manifest + cleaned filelists to
`/store/group/lpcmetx/SIDM/census/<run-id>/`, re-run `python sidm/scripts/build_census_index.py`,
and commit. That keeps the repo default reproducible and reviewable in the diff.

## Adding a new check

Each verdict/flag is a small, self-contained detector: a probe (`probe_file`, which already opens
every file) sets a field on `ProbeResult`, `_rollup` counts it, `render_census` prints a section,
and — if it should drop a file — `cleaned_filelists` adds it to the skip-list. The runs-count and
GenPart-self-reference flags follow exactly that shape, so a new failure or warning mode discovered
on, e.g., data is a few lines along the same path, not a re-plumb.
