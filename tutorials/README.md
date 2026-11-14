# Tutorials

Guides for students and analysts new to the `sidm` package. Each file covers one topic,
so you can go straight to the one you need. The first time through, read them in order.

## Topics

- [`01_first_run.ipynb`](01_first_run.ipynb): run `SidmProcessor` on a 200-event signal
  file committed to the repository, pick a channel and a histogram collection, print the
  cutflow, draw a histogram, and read everything the processor returns. No grid proxy or
  EOS access needed.
- [`02_normalization.md`](02_normalization.md): what the weighted numbers mean.
  Luminosity × cross section × generator-filter efficiency / sum of weights, why signal is
  at 1 fb, how the production tag reaches the processor, and what the Condor merge checks.
- [`03_data_and_blinding.md`](03_data_and_blinding.md): where the code stops you from
  processing collision data in a signal region and where it does not, the golden JSON, and
  the working rule. Read it before you run on data.

## Before you start

You need a working environment on LPC. The top [README](../README.md) covers it; this
folder does not repeat it:

- accounts, grid certificate and Kerberos: [Prerequisites](../README.md#prerequisites-one-time-from-your-laptop);
- the venv and the `SIDM (LCG_107 Py3.11)` kernel: steps 1 to 4 of
  [Getting started on LPC](../README.md#getting-started-on-lpc);
- opening notebooks on LPC: step 5 of the same section;
- reading the real samples (every topic after 01): a VOMS proxy (same section) and
  [6. Reading remote ROOT files](../README.md#6-reading-remote-root-files).

## Which way to run, by scale

The processor and its output are the same on every path; what changes is who runs the
chunks and how the output is put together.

| scale | way to run | normalization | data guard | read |
|---|---|---|---|---|
| one file or a few, quick checks | locally in a notebook, coffea `IterativeExecutor` | one `postprocess` over the whole run | none | `01_first_run.ipynb`; README [section 6](../README.md#6-reading-remote-root-files) for real files |
| tens to hundreds of files, while you iterate | dask from a notebook, workers on LPC Condor (`scaleout.make_lpc_client`, coffea `DaskExecutor`) | one `postprocess` over the whole run | none | README [section 7](../README.md#7-scaling-jobs-from-a-notebook-with-dask--htcondor); `sidm/test_notebooks/lpc_dask_example.ipynb` |
| thousands of files, unattended | Condor jobs (`condor/run_job.sh` runs `condor/run_sidm_chunk.py` per chunk), then `sidm/scripts/merge_coffea_chunks_eos.py` | each job normalizes its own chunk; the merge corrects this | `run_sidm_chunk.py` refuses data without `--allow-data` | `condor/README.md`; `sidm/test_notebooks/lpc_condor_example.ipynb` |

The normalization and data-guard columns are explained in `02_normalization.md` and
`03_data_and_blinding.md`.

## Older notebooks these replace

The older notebooks stay in `sidm/test_notebooks/` for now; where a tutorial here covers the
same ground, start with the tutorial.

| topic | replaces in `sidm/test_notebooks/` |
|---|---|
| 01 | `test_singleFile.ipynb` and `test_SidmProcessor.ipynb` (one file through the processor, cutflow, plots), `save_output.ipynb` (saving and loading an output) |
| 02 | none; it collects what was spread over the top README (section 8), `condor/README.md` (section 17) and `sidm/studies/limit_plotting/README.md` |
| 03 | none |

`lpc_dask_example.ipynb` and `lpc_condor_example.ipynb` stay as the walkthroughs of the
two scaled paths (dask from a notebook, and Condor jobs). Other notebooks that remain the reference
for their subject: `sidm/studies/corrections/00_correctionlib_tutorial.ipynb` (the standard
corrections), `sidm/test_notebooks/debug_branch_demo.ipynb` (per-event debug output),
`sidm/test_notebooks/file_census.ipynb` (which input files are good) and
`sidm/test_notebooks/output_browser.ipynb` (finding shared outputs on EOS).

## To come

These topics wait for open work that changes the code they would describe:

- **Condor campaigns from start to finish**: waits for PR #85, which changes how
  campaigns are configured and submitted. Until then: `condor/README.md`.
- **Loading outputs and the plotting helpers**: waits for issue #76, which brings the
  loaders and plotting style into one place. Until then: `output_browser.ipynb` and the
  top README section "Error bars on plots".
- **Defining selections and histograms**: waits for PR #71, which changes
  `sidm/definitions/cuts.py`, `sidm/definitions/hists.py` and the processor. Until then:
  the top README section "Analysis how-tos".
