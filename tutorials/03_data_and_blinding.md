# 03 Collision data and blinding

The analysis is blinded: nobody looks at collision data in a signal region until the
group decides to unblind. This page says which parts of the code stop you from doing it
by accident and which do not, how the golden JSON is applied, and the working rule that
follows. Read it before you put a data sample into any fileset.

## How the processor treats data

Data samples are listed in `sidm/configs/ntuples/data.yaml` and
`sidm/configs/ntuples/data_skimmed.yaml`, and every entry there carries `is_data: true`.
That flag reaches `SidmProcessor` through the fileset metadata (from
`utilities.make_fileset`, or from `condor/run_sidm_chunk.py`, `sample_metadata`) and
switches the processor to its data path:

- the golden JSON is applied (next section);
- every event has weight 1 and no corrections are applied;
- `postprocess` does not scale the sample.

## The golden JSON

For data, `SidmProcessor.process` first drops every event outside the certified
luminosity sections, before any object is built or any cut is evaluated:

- the JSON file is named by `golden_json` for the sample's `year` in
  `sidm/configs/run_periods.yaml`, and read from `sidm/data/`. For 2018 it is
  `Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt`; `sidm/data/README.txt` says
  where it was copied from;
- the mask is coffea's `lumi_tools.LumiMask`, evaluated on `run` and `luminosityBlock`;
- each chunk prints `Applying 2018 Golden JSON: <N> events removed`, and the total is
  stored in `out[sample]["metadata"]["n_removed_golden_json"]`.

So the first cutflow row (`None`) of a data sample counts certified events only.

Two consequences:

- The JSON is applied only when `is_data` is true. A data file in a hand-built fileset
  with `is_data: False` gets no golden JSON and is sent down the simulation path, which
  reads `genWeight`, a branch collision data does not have.
- Only 2018 has an entry in `run_periods.yaml`. Data with any other `year` stops with a
  `KeyError` in `process`.

A chunk can lose every event to the mask. The committed data file
`tests/data/events_DoubleMuon2018C_200ev.root` is built so that its first 100-event chunk
is entirely outside the JSON, and CI runs it to check that this case does not crash
(`tests/README.md`).

## Guards that exist

**Condor jobs.** `condor/run_sidm_chunk.py` refuses collision data unless it is given
`--allow-data`. A job's sample counts as data if its location YAML entry has
`is_data: true`, if that YAML's file name starts with `data`, or if any input path
contains `/SIDM/Data/`. The job then stops before reading any file:

```
ERROR: <sample> is collision data (<yaml>:<version>). The analysis is blinded; pass --allow-data to process it on purpose.
```

`condor/run_job.sh`, the executable of the Condor submit files, never passes
`--allow-data`.
The data submit files in `condor/` (`submit_*_data.sub`) say so in their headers: they
record how earlier data inputs were made and do not run as they are. This guard works per
sample, not per channel: it stops data for control regions as well as signal regions.

**Limit and ABCD readers.** The background estimate is an ABCD method: two variables split
the selected events into four regions, and the background in the signal-like region, A, is
predicted from the counts in the other three. In `sidm/studies/limit_plotting/datacard_tools.py`,
`region_yields` and `collect_yields` leave out region A of the ABCD plane
(`SR_ABCD_REGION`, the signal region) for every sample that `is_simulation` cannot
establish as simulation, and `sr_yield` raises `BlindingError` for such a sample.
`is_simulation` trusts `metadata["is_data"]` when it holds one value; otherwise it accepts
only signal names and samples listed in `cross_sections.yaml`, and treats everything else
as data. In `sidm/studies/limit_plotting/shape_tools.py`, `region_shapes`,
`collect_shapes`, `mljlj_two_bin`, `region_shapes_6d` and `collect_6d` apply the same rule,
and the directory readers refuse to open a campaign's data merge directory at all. All of
these functions take `allow_data_sr=False`; the `datacard_tools` docstrings say to set it
to `True` only when the analysis has been unblinded.

## Where there is no guard

- `SidmProcessor` runs any channel on data, including the signal-region selections, and
  fills every histogram you ask for. Nothing in it compares the channels with `is_data`.
- `utilities.make_fileset` builds filesets from `data.yaml` and `data_skimmed.yaml` like
  any other location file.
- The dask path (`sidm/tools/scaleout.py`, `make_lpc_client`, with coffea's
  `DaskExecutor`) adds no check.
- `sidm/scripts/merge_coffea_chunks_eos.py` merges data chunks like any others.
- Loading and plotting (`coffea.util.load`, `sidm/tools/output_browser.py`,
  `utilities.plot` and the other plotting helpers in `sidm/tools/utilities.py`) show
  whatever the output contains. `output_browser.query` can filter on `is_data`, but does
  not hide anything.

In short: in a notebook or on dask, the only protection is you.

## Working rule

1. Do not run a signal-region channel on data, in any notebook, on dask or on Condor.
   There is no list of signal-region channels in the code. The selections named `*_sr` in
   `sidm/configs/selections.yaml` are signal regions
   (`python3 sidm/scripts/print_cuts.py --list -s '*_sr'`), and so is region A of the ABCD
   histograms the limit code reads. A study's own selections file can define more; check
   each channel with `print_cuts.py` before running it on data.
2. Before adding a data sample to a fileset, read the channel list and the histogram
   collections once more. A histogram that spans the ABCD plane contains region A even
   when the channel is not called `_sr`.
3. Read data outputs for ABCD and limit work through the `datacard_tools` and
   `shape_tools` functions above, with `allow_data_sr` left at `False`.
4. `--allow-data` is an option of `run_sidm_chunk.py` only. No submit file and no version
   of `run_job.sh` passes it, so a data run on Condor needs a deliberate change to the job
   script: agree on it in the group first, and use channels that exclude every signal
   region. How Condor campaigns are configured is changing (PR #85).
5. Data are 2018 only, and the golden JSON is applied for you as long as the fileset says
   `is_data: true`. Build data filesets with `make_fileset` from the data location files so
   that it does.
