# Working in this repository

How changes to the SIDM analysis code are made. These rules are for everyone: human
contributors and coding agents alike (Codex reads this file; Claude Code reads
`CLAUDE.md`, which imports it). When a rule here conflicts with an older habit in the
code, follow the rule and say so in the pull request.

## The analysis is blinded

- Never run on, open, or print anything derived from collision data unless the analysis
  group has approved that specific use. Collision data means the DoubleMuon samples in
  `sidm/configs/ntuples/data*.yaml`, anything under `/store/.../SIDM/Data/`, and the CI
  fixture `tests/data/events_DoubleMuon2018C_200ev.root`.
- Studies use simulation by default. Checks that need data (data/MC comparisons, closure
  in control regions) are written up as proposals for the group, not run.
- `condor/run_sidm_chunk.py` refuses data unless given `--allow-data`, and Combine runs
  blinded (`--run blind`) by default. Do not change either without approval.

## Extend the framework; do not build a parallel one

Before writing anything new, look for it: `sidm/definitions/hists.py` holds thousands of
histogram definitions, and `sidm/tools/utilities.py` and the `*_study.py` modules of the
studies hold shared helpers. Reuse and extend them.

- **Histograms.** Every plot or number comes from a histogram defined in
  `sidm/definitions/hists.py`, grouped in a collection in
  `sidm/configs/hist_collections.yaml`, and filled by `SidmProcessor` (see "How to define
  a new histogram and add it to a collection" in `README.md`). Reuse the helpers and fill
  functions already there. Debug branches are only for one-off checks.
- **Objects, cuts, selections.** A new object goes in `sidm/definitions/objects.py`, a new
  cut in `sidm/definitions/cuts.py`, and a selection is composed in
  `sidm/configs/selections.yaml` by inheriting from the existing anchors (`base`, `2mu2e`,
  `4mu`, ...). `sidm/scripts/print_cuts.py` shows what a selection resolves to.
- **Corrections and weight systematics.** Use `sidm/configs/corrections.yaml` and
  `sidm/tools/corrections.py`. Corrections are event weights, and their variations come
  out on the `syst` axis with `syst_hists=True`. A correction of our own is a
  correctionlib file registered there (`sidm/studies/corrections/README.md`).
- **Systematics in the limits.** Add one entry to `sidm/studies/systematics/systematics.yaml`,
  then rerun, in this order, `03_limits.ipynb` (limits with and without), and
  `04_systematic_effects.ipynb` and `05_fit_checks.ipynb` (what each systematic does, and
  checks of the fit), all in `sidm/studies/systematics/`. The last two read the outputs of
  the first and stop if they are not from the current registry.
- **Producing samples.** Always a Condor campaign; see the next section.
- **Where outputs live.** Inputs and outputs of studies (merged coffea files, datacards,
  limits, large tables) live on EOS, not in the repository, and notebooks reach them through
  xrootd URLs (`root://cmseos.fnal.gov//store/...`, read and written with fsspec), so the
  same notebook runs at the LPC and on CoffeaCasa, where EOS is not mounted. Scratch work
  goes in a temporary directory outside the repository.
- **Limits.** Use the tools in `sidm/studies/limit_plotting` (`datacard_tools.py`,
  `shape_tools.py`, `plot_tools.py`) and `sidm/scripts` (`abcd_datacards.py`,
  `abcd_workspace.py`, and the `run_combine_*.py` scripts for limits, pulls, impacts,
  likelihood scans and toys), and extend them where needed.

## Producing samples: always a Condor campaign

Every run of `SidmProcessor` over full samples is a campaign of `condor/campaign.py`
(`condor/README.md`, section 23). Do not write job lists, submit files, retry lists or
merge commands by hand, and do not run the processor over full samples locally or with
Dask; local runs are for tests on a few files. Sections 8 (packaging the code) and 11
(the proxy) of `condor/README.md` still apply; `campaign.py` replaces the manual steps of
the other sections and `condor_campaign.py`.

1. Copy `condor/campaign_configs/template.yaml` to `condor/campaign_configs/<study>_v1.yaml`
   and fill it in: channels from `sidm/configs/selections.yaml`, histogram collections
   from `sidm/configs/hist_collections.yaml`, corrections and `hist_mode`, the sample
   lists, and the chunk and merged directories on EOS. To test a config, give a copy its
   own name and directories and set `max_files`.
2. `python condor/campaign.py <name> jobs` makes the job lists.
3. Remake the code tarball if anything under `sidm/` or `condor/run_sidm_chunk.py` changed
   since it was last made (section 8), and set up the proxy with its copy on the shared
   disk (section 11).
4. `python condor/campaign.py <name> submit`
5. `python condor/campaign.py <name> status`, then `retry --submit` and `skipped --submit`,
   until every chunk is DONE and no skipped file is left to retry.
6. `python condor/campaign.py <name> merge`. Merged outputs go to
   `/store/group/lpcmetx/SIDM/coffea_outputs/$USER/<study>/`, or to the user area on EOS
   until `lpcmetx` write access is granted.

Notebooks and limit scripts read the merged outputs over xrootd. A production for the
limits is registered as a `Campaign` in `sidm/studies/limit_plotting/datacard_tools.py`.
Changing what a campaign produces (selection, histograms, corrections) means a new config
with the next version (`_v2`), so outputs already in use do not change. If the campaign
tool lacks something a study needs, extend `condor/campaign.py` and `template.yaml`, with
a test in `condor/test_campaign.py`, rather than working around it. Commit the config;
the job lists, filelists, logs and `condor/campaigns/<name>/` stay out of git. A campaign
over collision data needs the group's approval first (see "The analysis is blinded").

## Code

- Prefer the simplest solution that works. Readable, maintainable code beats clever
  code. Keep changes minimal and focused on the task.
- Note possible speedups and refactors in the pull request instead of making them
  unasked.
- Match the style of the surrounding code. Give every function a short docstring saying
  what it does and, where it is not obvious, why.
- Before committing, run pylint on the changed files:
  `pylint --rcfile .pylintrc <files>`. CI fails on new errors and on a whole-repository
  score below 8.
- A study gets its own folder under `sidm/studies/`, with a `README.md` saying what it
  does, how to run it, and what it found, and with its shared notebook code in a
  `*_study.py` module.
- A study in progress also keeps a `PROGRESS.md` in its folder: the current state at the
  top (what exists, what is running, what waits on whom, the next steps, how to resume) and
  a dated log below, newest first. Read it before working on the study, and update it as
  part of every piece of work on it, so that anyone, or any new agent session, can pick the
  work up from there.
- In documentation, comments and notebooks, write plain, concise English and use no em
  dashes.

## Environment

- The Python stack is pinned in `requirements.txt` and `constraints.txt`. The top-level
  `README.md` explains the setup. Correction files are read from CVMFS; CI uses the copies
  in `tests/data/corrections` (`SIDM_CORRECTIONS_ROOT`), so a correction file newly read
  in CI needs a copy there.
- Combine runs in a separate CMSSW release, which the scripts in `sidm/scripts` dispatch
  to themselves.

## Commits and pull requests

- Coding agents commit or push only when asked.
- Stage files by explicit path. Never use `git add -A` or `git add .`: working trees hold
  large untracked outputs (coffea files, condor logs, notes, slides) that must not be
  committed.
- Do not commit `.coffea` or `.root` files, condor logs, job outputs, or anything derived
  from data.
- CI runs the chain report (the processor over the simulated fixture, compared with
  `main`; `tests/README.md`) and pylint. Read both on every pull request.
