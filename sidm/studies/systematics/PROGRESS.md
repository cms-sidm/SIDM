# Systematics study: progress log

Read this first when picking the study up again, then [`README.md`](README.md) for the
methods and results. Update it as part of every piece of work on the study: the current
state below, and a dated entry in the log (newest first). [`AGENTS.md`](../../../AGENTS.md)
asks for this.

## Current state (2026-11-17, 11:10)

**Step 2 of the lepton-jet efficiency done for the prompt regime (`lj_reco_efficiency.md`).** The
user set the 4mu working-point deck's safety criteria as the guide for reading data (now a general
rule in AGENTS.md) and allowed committing results from approved data uses.
- **Simulation:** campaign `lj_reco_cr_v2`, notebook 07: every region of both control samples is
  safe (with the sparse-bin rule, switch `STRICT_SPARSE_BINS`, used until the advisor rules).
- **Data:** campaign `lj_reco_cr_data_v1` complete (2,590 chunks: 2,580 plus 10 extra chunks for 65
  files that failed to read at first; 463 empty files skipped, as in `systematics_data_v1`; 68
  samples merged). Notebook 08, read through the safety table: muon-type (J/ψ, above 50 GeV)
  data/simulation 0.988 ± 0.006 (endcaps 0.976 ± 0.008); egamma-type (Z → μμγ, unconverted
  photons above 30 GeV) 1.01 ± 0.04, 0.80 ± 0.08 above 60 GeV (112 probes).
- **AN Sec. 7:** the lepton-jet reconstruction subsection now holds the measurement (two figures,
  `efficiency_jpsi_above50_data_simulation.pdf` and `efficiency_zmmg_data_simulation.pdf` in
  `an_systematics/`), replacing the authors' placeholder, in the advisor's style; compiles.
- **Deck** `slides/systematics.tex`: 40 pages with the data results (local).
- **Upstream `main` merged** (2026-11-17): it holds its own version of the campaign tool, with the
  same `allow_data` opt-in, so the merge takes `main`'s campaign tool and re-applies only
  `file_stride`. `main` also selects the leading lepton-jet pair that survives the subtype cuts
  (commits 0404af3, f0c4836 and later, S. Abbott). The control-sample histograms are identical
  before and after the merge (checked on Drell-Yan, QCD and a signal file), but the signal-region
  yields change (`abcd6d_2mu2e` -1.3% for 2Mu2E 500/1.2/30 cm, `abcd6d_4mu` -3.3% for 4Mu at the
  same point): the campaigns `systematics_v2` and `systematics_data_v1`, and every number of
  notebooks 03 to 05 and AN Sec. 7, predate it. A rerun (`systematics_v3`, `systematics_data_v2`)
  would align them with `main`; the new pair selection may also settle the spread-cut part of open
  question 4 (to check).
- **Commits** on `systematicsStudy`, not pushed: `25fe827` (campaign tool), `0a5fd77` (AGENTS.md),
  `3dd4114` (the study up to notebook 07), and the data step (notebook 08 and the documents).

**The user's workflow (2026-11-16):** after running anything, update the AN (in the advisor's
style, `an_systematics/WRITING_GUIDE.md`), then the summary deck, then commit and paste an updated
PR description in chat. Keep `slides/` and `an_systematics/` local.

**Next:** step 3 (per-lepton-jet scale factors for the prompt regime, `Corrections` extension,
two registry entries, rerun of 03 to 05), after the group's answers to questions 2 to 6 of
`lj_reco_efficiency.md`; a purer Z → μμγ probe for converted photons (a new version of the sample).

## State before step 2 (2026-11-16, 17:10)

**Done: the advisor's comments on the PR, with data-driven results, and AN Sec. 7.** The
advisor asked (via the user) for a data vs simulation plot of the number of primary vertices in
the AN, for the systematics section to be fixed up in general, and for the systematic results
and Combine checks to be based on data-driven expected limits; he approved looking at data in
the control regions we suggested. The user approved: data through `Z_to_mumu_control_region`,
good-PV count only; and data through `SR_2mu2e`/`SR_4mu`, regions B, C and D only (m_LJLJ
included), region A never read, Combine blinded. Every data campaign runs on Condor and saves on
EOS (the user's rule).
- **Campaigns complete:** `npv_data_v1` (172 of 172 chunks, every 100th file of the 2018
  SkimNano V2 DoubleMuon data; 5 skipped files, all empty) and `systematics_data_v1` (2,580 of
  2,580 chunks, all 127,415 files; 463 skipped files, all empty); 68 data samples merged in each.
  30 slow jobs moved (nodes cmswn2176 and cmswn2236 excluded for every queued job).
- **Notebooks rerun:** 01 section 6 (NPV in data against the simulation), 03 (data-driven
  cards, integrated and 8-bin estimates), 04 and 05 (on the integrated cards), each with
  rewritten interpretations; the README is updated to match.
- **AN Sec. 7:** `an_systematics/` (git-ignored as a whole) holds the finished `systematics.tex`,
  its three figures (`npv_data_vs_simulation.pdf` and the two impacts plots), `references.bib`
  (three new references) and `WRITING_GUIDE.md`. Checked against the current `results.tex`
  (labels, macros, conventions) and compiled with stub macros: no errors, the only undefined
  references are those defined elsewhere in the AN. It states that the selection here is not
  the one of Sec. 8 (open question 4) and gives the predicted backgrounds (53.6 events in 2mu2e
  against 38.7 in Sec. 8, and 108 in 4mu at 3 pixel hits).
- **Code, uncommitted:** `campaign.py` `allow_data` (passes `--allow-data` through
  `submit_campaign.sub` and `run_job.sh`; default false) and `file_stride` (every Nth file),
  with tests (61 pass); the configs `npv_data_v1.yaml` and `systematics_data_v1.yaml` and the
  sample list `condor/samples_data_doublemuon_2018.txt`; collections `abcd_limits_data` and
  `pv_npvs_good`; `abcd_workspace.py --sr-data-blinded` (B, C, D from the payload's data block,
  region A's observation the prediction, each workspace checked bin by bin; free bins span 50
  times their content); `systematics_study.data_control_regions`, `write_payloads(data=...)`,
  `an_hypotheses`, `LIMITS_SUFFIX` (`_limits_data`), `LIMITS_METHOD` (`integrated`) and the
  NPV helpers; `systematics.yaml` with the three lepton scale-factor uncertainties; the
  `systematics_v2` Campaign's `data_dir`; notebooks 01, 03, 04 and 05; the README, AGENTS.md and
  `condor/README.md` notes on data campaigns.

**Next steps**
1. The user: review `an_systematics/` and say what to commit (the code, configs and notebooks
   above; `an_systematics/` and `slides/` stay untracked).
2. The lepton-jet efficiency with data (`lj_reco_efficiency.md`, step 2): scope each control
   sample before its first read; the group questions below still apply.

**Earlier, still open:** the slides (`slides/`, untracked, built on the simulation-driven v2
results) await the user's decision; a campaign `lj_reco_eff_v2` splitting the muon-type lepton
jets by constituent was offered.

**Notebooks** (inputs and outputs on EOS, read and written through xrootd)

| notebook | state | input | outputs |
|---|---|---|---|
| `01_pileup` | executed 2026-11-16 (section 6 with data) | campaigns `corrections_v1`, `npv_data_v1` | in the notebook, NPV figure on EOS |
| `02_jet_energy` | executed 2026-11-12 | campaign `corrections_v1` | in the notebook |
| `03_limits` | executed 2026-11-16, data-driven, with interpretations | campaigns `systematics_v2`, `systematics_data_v1` | `.../systematics_v2_limits_data/` |
| `04_systematic_effects` | executed 2026-11-16, with interpretations | outputs of 03 | `.../systematics_v2_limits_data/effects/` |
| `05_fit_checks` | executed 2026-11-16, with interpretations | outputs of 03 | `.../systematics_v2_limits_data/fit_checks/` |
| `06_lj_reco_efficiency` | executed 2026-11-14, with interpretations | campaign `lj_reco_eff_v1` | `.../lj_reco_eff_v1_results/` |

`...` is `root://cmseos.fnal.gov//store/user/scampbel/sidm_coffea_outputs`. After a change to
`systematics.yaml`, rerun 03, then 04 and 05; 04 and 05 stop if 03 is out of date. The
simulation-driven outputs of 2026-11-14 stay in `.../systematics_v2_limits/`.

**Campaigns** (configs in `condor/campaign_configs/`, run with `condor/campaign.py`, working
state in `condor/campaigns/<name>/`)
- `corrections_v1`: input of 01 and 02 (`.../corrections_v1`).
- `systematics_v1`: complete on 2026-11-13; superseded by v2, since it lacks the
  generator-filter efficiency.
- `systematics_v2`: complete on 2026-11-14 (3,158 jobs, 200 samples merged in
  `.../systematics_v2`); signal and simulated background of 03 to 05.
- `npv_data_v1`: complete on 2026-11-16 (`.../npv_data_v1`); data of 01 section 6.
- `systematics_data_v1`: complete on 2026-11-16 (`.../systematics_data_v1`); data of 03 to 05,
  read for regions B, C and D only.
- `lj_reco_eff_v1`: complete on 2026-11-14 (422 jobs, 180 samples merged in
  `.../lj_reco_eff_v1`); input of 06.

**Registry** (`systematics.yaml`, seven entries): `lumi_13TeV` flat 0.83%, `CMS_pileup` from
the pileup weight variations, `CMS_l1_prefiring` flat 1%, `CMS_jet_energy` from the isolation
scaled by 1.10, and `CMS_eff_m_id`, `CMS_eff_e_reco`, `CMS_eff_g_pixelveto` from the scale-factor
variations.

**Results to keep in mind** (details in the README and in the notebooks' interpretations)
- Data control regions: 2mu2e B, C, D = 79, 19, 28 (prediction 53.6); 4mu at 3 pixel hits 45,
  12, 5 (prediction 108). AN Sec. 8 has 66, 17, 29 (38.7) in 2mu2e, with its own selection.
- NPV in the Z region: mean 24.8 in data and 24.4 in the simulation, inside the 23.4 to 25.4 of
  the pileup variations; the data are wider. Reweighting the signal by the data/simulation
  ratio changes it by a median 0.2%, at most 1.0%.
- With corrections and systematics the data-driven expected limits are higher by a median of
  5.2% (2Mu2E) and 1.4% (4Mu) over the AN's 60 hypotheses per final state, nearly all from the
  corrections (signal in region A 4.6% and 1.4% lower); the seven systematics alone change no
  limit of the grid by more than 0.27%.
- Impacts about each nuisance's own size, at most 1%, against 60 to 100% statistical
  uncertainty. Signal injection: every fit converges in both channels, median pulls −0.15 to
  0.01, half-widths 0.96 to 1.02; in 4Mu a quarter of the toys at the expected limit fit r near
  zero (region D holds 5 events), so the mean fitted r is 11 to 12% low.
- At the shortest lifetime, five 2Mu2E points have 22 to 38% signal statistics in region A, and
  their pileup values reach 22% (one limit moves by 4.7%).
- Lepton-jet efficiency, step 1: on the limit grid only the shortest lifetime decays mostly
  inside the pixel layers; the 4Mu muon-type efficiency is 0.91 below 10 cm, 0.82 from 10 to
  100 cm and 0.33 beyond 300 cm; egamma-type 0.89, 0.80, zero beyond 140 cm.

**Open questions for the user or the group**
1. The 4mu method and working point: AN Sec. 8 uses at most 2 pixel hits (C and D empty in
   data); here 3 pixel hits, where the estimate rests on 5 events in region D: group.
2. Lepton-jet reconstruction efficiency (`lj_reco_efficiency.md`, "What is needed"): scale
   factors per lepton jet or per constituent, the provisional value until step 2, and how to set
   the uncertainty of the most collimated muon pairs (m_Zd = 0.25 GeV at 800 and 1000 GeV),
   which no control sample reaches: group. Data control samples may now be read (the advisor's
   approval), one at a time and scoped first.
3. The shortest-lifetime points: cap their per-point values or take pileup from the
   neighbouring lifetimes, if they stay in the limits; and whether to cover the wider NPV
   distribution in data explicitly (README open item 4): group.
4. The selection here is not the one of AN Sec. 8: the HEM electron veto and the `mu_veto` are
   not defined in the repository, and the dz and dxy spread cuts of `SR_2mu2e`/`SR_4mu` filter
   only the `mu_ljs` view, so they never reach the `abcd6d_*` histograms. A fix needs the
   vetoes' definitions, a decision on the spread cuts (drop the lepton jet or veto the event)
   and new campaigns: told the user, not changed.
5. Whether 03 should carry the agreed 4mu bin merging: offered to the user, not done (the
   20% per-bin uncertainty for 2mu2e is in the binned data-driven cards).
6. Write access to the `lpcmetx` group area: until then merged outputs stay in the user area.
7. Lepton-jet control samples that our inputs cannot give (`lj_reco_efficiency.md`, question 5):
   cosmics (no cosmic simulation; both legs always reconstructed in our triggered, skimmed data)
   and J/psi -> ee (needs EGamma data), so the displaced muon-type regime has no data
   measurement: produce those inputs or cover the regime otherwise: group.
8. The sparse-bin rule of the safety check (`lj_reco_efficiency.md`, question 6): used with the
   user's agreement (2026-11-16) until the advisor rules; `STRICT_SPARSE_BINS` in
   `systematics_study.py` switches it off: advisor.

**Committed** (2026-11-14, branch `systematicsStudy`, for a pull request to `main`, not
pushed): the study folder, the framework and limit-tool changes it needs, the campaign configs,
and AGENTS.md and CLAUDE.md; upstream `main` merged in (`39a2883`), the fix for its
normalization (`8e337fc`), the rerun (`2a396f0`) and its log (`298a412`). The branch also
carries the campaign-configs commit `430542a`, not yet in `main`. Left out on purpose
(AGENTS.md): job lists, filelists, logs, `condor/campaigns/`, outputs, and the local notes in
`sidm/studies/limit_plotting/`. Everything of 2026-11-16 is uncommitted; stage by explicit path
only.

**Environment**
- `cmsenv` in `CMSSW_17_0_0_pre4/src`, with `PYTHONUSERBASE` set to the repository's
  `.python_userbase` (`/usr/bin/python3` lacks numba).
- Combine: `/uscms_data/d3/scampbel/CMSSW_14_1_0_pre4` (v10.2.1, with `combineTool.py`); the
  scripts in `sidm/scripts` set it up themselves.
- Grid proxy for Condor and EOS: `X509_USER_PROXY` set to its copy on the shared disk
  (`condor/README.md`, section 11); the current one is valid to about 2026-11-21.
- Executing a notebook, from this folder: `python3 -m jupyter nbconvert --to notebook --execute
  --inplace --ExecutePreprocessor.timeout=7200 --ExecutePreprocessor.kernel_name=python3 <notebook>`.

## Log (newest first)

### 2026-11-17

- (11:20) Upstream `main` had moved (15 commits, with its own campaign tool, reviewed in #85, and a
  new lepton-jet pair selection). Merged it: `main`'s campaign tool kept, `file_stride` re-applied
  with its test; `sidm_processor.py` keeps both `needs_probe_inputs` and `select_lepton_jets`.
  Tests pass (151). Regression check before and after the merge: control samples identical, signal
  region changed (see the current state).
- (11:10) Data campaign merged after the 10 extra chunks; notebook 08 written and executed (data
  read only through `lj_cr_data`, gated by the safety table): results above. The J/ψ figure for
  the AN shows the pT dependence from 30 GeV and the opening angle and |η| above 50 GeV, where the
  threshold no longer matters; the comparison figures now share their x axis between the panels
  (they did not at first). AN subsection, deck, README and proposal updated.
- (09:40) The data campaign finished overnight; 65 skipped files in 10 chunks resubmitted. The user
  allowed committing results from approved data uses: AGENTS.md now says so, and states the
  group's safety criteria for reading a data region as a general rule. Committed the work of
  2026-11-16 (campaign tool, AGENTS.md, the study).

### 2026-11-16

- (23:10) `lj_reco_cr_v2` complete at 22:35 (all 1,926 chunks, no failures; 72 moves). Its watcher
  had judged every job against the campaign-wide median (18 min, set by the fast signal chunks),
  so at 20:41 to 20:47 it restarted 72 DYJetsToLL_M10to50 chunks that were running normally for
  this campaign (52 min; the probe histograms make chunks 1.2 to 1.8 times slower than in
  `systematics_v2`) and excluded 30 nodes; fixed to per-sample limits and the exclusions reset
  (about 40 min lost). Notebook 07 executed: every region safe; the first Z → μμγ efficiencies
  were above one because the non-Z photons rise steeply under the peak, so the Z windows are now
  85 to 95 GeV with adjacent sidebands (`LJ_CR`), and section 4 shows the Z efficiency above 30 GeV
  with a Drell-Yan curve. Data campaign `lj_reco_cr_data_v1` submitted at 22:57.
- (20:05) The advisor's writing rule (from the review of `results.tex`) generalized in
  `WRITING_GUIDE.md` ("The AN describes the analysis as it is", with a search-word checklist) and
  applied to the parts of Sec. 7 written here: the selection caveat, the comparison with the
  Sec. 8 background (53.6 against 38.7), the "loosened" working point, the forward reference of
  the scale-factor paragraph to the unfinished lepton-jet section and the shortest-lifetime limit
  shift (4.7%) moved out of the AN into the deck; the shortest-lifetime values now justify quoting
  the sizes for L_xy >= 3 cm. Compiles. The trigger and lepton-jet placeholders are the authors'
  text and stay. Summary deck drafted (33 pages, no overfull frames).
- (19:50) The user agreed to the sparse-bin rule until the advisor rules (switch
  `STRICT_SPARSE_BINS`, question 6) and asked that the control samples we cannot use be a group
  question (question 5). Added the `_alone` versions of both control samples, so that no event of
  the signal region enters what data will read (a J/ψ event with a second lepton jet could sit in
  region A); the refactored v1 histograms reproduce the earlier local output bin for bin. Since
  the safety check must judge what data will read, `lj_reco_cr_v1` was stopped at 24% and
  `lj_reco_cr_v2` (both versions in one pass) submitted at 19:40. Notebook 07 written. The
  watcher crashed once on a job without a start date yet (fixed, restarted).
- (19:20) Lepton-jet efficiency, step 2 started. Scoped the four control samples against our
  inputs (two feasible, two not; `lj_reco_efficiency.md`). Simulation scan of a few skimmed files:
  in Drell-Yan about 0.1% of the Z-region events have an FSR photon probe; QCD and tt̄ events
  passing the skim often hold a J/ψ candidate (b-hadron decays), Drell-Yan none with a generator
  J/ψ. Wrote the two channels, probes, histograms and collection; local test on one file each of
  Drell-Yan, QCD and a 4Mu point (the DSA cross-cleaning warning is pre-existing: it also appears
  for `SR_2mu2e`). The user's guide for data reads is the 4mu working-point deck's criteria; the
  safety check implements them per probe-pT bin with the stress test and the union check, and adds
  a sparse-bin verdict because the J/ψ test put about 5 signal events in bins without simulated
  background while the whole histogram passed the decks' verdicts. Campaign `lj_reco_cr_v1`
  submitted (1,926 jobs); pylint: no E or F message in the changed files.
- (17:10) AN Sec. 7 finished and proofread against the current `results.tex`: the summary
  table lists median and largest values for L_xy >= 3 cm, re-derived from the EOS outputs
  (the pixel-veto maximum is 0.35%, so 0.3%; the electron and photon rows are below 0.1% in
  4mu); the signal statistics enter as in Sec. 8 (one parameter per region), not as one common
  parameter; the evaluation states that its selection is not that of Sec. 8 and gives the
  predicted backgrounds. Compiled with stub macros: no errors. The same imprecise shortest-
  lifetime wording ("22 to 38%" holds for five of the twelve points) fixed in the interpretation
  of notebook 04. README results of 04 and 05 rewritten for the data-driven run; open items 2
  and 3 updated.
- (16:50) Notebooks 03, 04 and 05 executed on the data-driven inputs, with interpretations: the
  limits with corrections and systematics are a median 5.2% (2Mu2E) and 1.4% (4Mu) higher, nearly
  all from the corrections; impacts at most 1%; every toy fit converges; in 4Mu a quarter of the
  toys at the expected limit fit r near zero, which a counting model of the four regions
  reproduces from the fluctuations of region D (5 events). Notebook 04 reports the limit
  breakdown over the AN's hypotheses (`card_name()` maps them to the existing cards).
- (15:20) Both data campaigns complete and merged (`npv_data_v1`: 172 chunks, 5 empty files
  skipped; `systematics_data_v1`: 2,580 chunks, 463 empty files skipped); 30 slow jobs moved in
  all. Notebook 01 section 6 executed: NPV in data against the simulation, and the signal
  reweighted by their ratio (median 0.2%, at most 1.0%).
- The advisor's comments (NPV data vs simulation, fix up Sec. 7, data-driven expected limits
  for the systematic results and Combine checks) and the user's approvals (above). Two data
  campaigns submitted after one-file local tests (deleted afterwards): the Z region gives about
  5000 events per data file, so every 100th file suffices for the NPV shape; one data file has
  no event in B, C or D of either plane, as expected for the tight selection.
  `abcd_workspace.py --sr-data-blinded` tested end to end with synthetic control counts (B, C, D
  = 21, 7, 10: transfer factor 0.7, prediction 14.7; check passed). The skim of the data
  (SkimNano V2) is the basic skim of `selections.yaml` (`trigger_muPtCut`: triggers and two
  muons with pT > 26 GeV), which both the signal-region channels and the Z region require.

### 2026-11-14

- Slides rebuilt with the v2 figures (21 replaced) and numbers, 40 pages, no overfull frames.
  The rerun committed as `2a396f0`.
- Notebooks 03, 04 and 05 rerun on `systematics_v2` with the 4mu working point at 3 pixel hits
  (20:03 to 20:38), interpretations and README rewritten. 03: background unchanged from v1 (4mu
  B·C/D 0.38; region A 20.9 events, 19.3 Drell-Yan, which no control region sees); every limit
  is the v1 one divided by its filter efficiency (checked to three decimals at all 192 grid
  limits); exclusion 12 of 48 (2Mu2E) and 23 of 48 (4Mu), with and without the corrections.
  04: impacts and templates as before; the limit breakdown shows four 2Mu2E points at 200 GeV
  that move only without the jet energy (0.1 to 0.25%), where it is largest, besides jumps of
  the limit calculation; the 2Mu2E breakdown is identical in v1 and v2, as the scaling
  requires. 05 failed in its correlations: FitDiagnostics evaluates the likelihood at r = 0
  before any fit, which is undefined where the injected signal is a bin's only content (4mu C
  300 to 400 GeV, A 1000 to 1500 GeV); `run_combine_pulls.py` now starts from the injected r,
  and both fits converge with status 0 (2Mu2E had status 1 before). 4Mu injection: a third of
  the fits fail, the converged ones 2.5 to 4% low (13 to 15% at 2 pixel hits).
- `systematics_v2` complete at 20:01: all 3,158 chunks done, nothing to retry, all 200 samples
  merged (194 at 19:00, 3 at 19:15, the last 3 at 20:10, each round after the skipped-file
  check of its chunks). 10 jobs moved off slow nodes in all. 03 started on it.
- `systematics_v2` at 19:00: 3,147 of 3,158 chunks done, no failures. 8 slow jobs moved off
  cmswn2236, cmswn2176 and cmswn2254 (for example 4Mu_150GeV_1p2GeV_0p64mm chunk 1: 47 min on
  cmswn2176 and not finished, then 11 min on the new node, 10 min in v1). Two of the 18:51 moves
  swapped onto the other slow node through claim reuse, so the watcher now excludes every slow
  node seen so far. Skipped-file check: the same 4 files without events as v1, nothing to retry.
  Merge of the 194 fully checked samples started.
- Campaign `systematics_v2` submitted at 17:01 (3,158 jobs; same chunks as v1, checked) after a
  local test on one signal and one background file: the signal output records
  `gen_filter_eff` 0.7795, the background 1.0, and the fixed `read_coffea` reads both. Code
  tarball remade. The watcher moves slow jobs (Current state).
- Merged upstream main into the branch (`39a2883`; 60 upstream commits, the user's request so
  the PR can merge). Conflicts in `hists.py` and `hist_collections.yaml` (both sides added;
  both kept) and in the `datacard_tools.py` docstring (rewritten to cover both conventions).
  Upstream now applies the generator-filter efficiency of the v10 signal samples in the
  normalization, which broke `production_xs_pb` for new outputs (the first cutflow row is no
  longer lumi * xs); fixed in `8e337fc` (the recorded efficiency is divided out before the
  match). Checked: no duplicate keys or names, every collection's hists defined, `tests/` and
  the campaign tests pass (120).
- The 4mu working point loosened to at most 3 pixel hits on both lepton jets (the user's
  correction). 03 rerun on v1 with it (4mu simulated background A 20.9, B 26.3, C 0.039, D
  2.75, B*C/D 0.38 instead of 0.003; 4Mu limits a median 11% lower; exclusion 24 to 27 of 48),
  superseded by v2.
- Summary deck in [`slides/`](slides/): `systematics.tex` (beamer, the style of the corrections
  deck), 40 pages: the four systematics, the limits with and without, what each systematic does,
  the fit checks, step 1 of the lepton-jet efficiency, the questions for the group, open items
  and backup. Build with `pdflatex -interaction=nonstopmode systematics.tex` (twice) in
  `slides/`. Its figures in `slides/figures/` are copies of the PDFs that 03 to 06 wrote to EOS
  and of PNG outputs of 01 and 02; like the deck PDF they are gitignored. The deck is not
  committed (made after the commits).
- Committed for a pull request, at the user's request, in four commits on `systematicsStudy`:
  the user's limit-tool changes of the local `limit_updates` branch (`6f04fbc`, until now
  applied uncommitted), the `campaign.py skipped` fix, AGENTS.md and CLAUDE.md, and the study.
  pylint: whole repository 8.76/10; no E or F message in the changed files apart from the
  ROOT `no-member` ones, which appear only where ROOT is importable (not in CI). The campaign
  tests pass.
- Lepton-jet reconstruction efficiency, step 1 done. Campaign `lj_reco_eff_v1` complete: 422
  jobs, no failures, nothing to retry (4 input files without events skipped), 180 samples
  merged. Its last job ran 82 minutes on a slow worker (cmswn2236: 22 minutes to preprocess
  25 files that read in under a second each from cmslpc, and 11 minutes for the pip install of
  its restart); it was excluded from that node with `condor_qedit` of `Requirements`, held and
  released, and finished in 11 minutes elsewhere. `06_lj_reco_efficiency.ipynb` executed (31 s)
  with interpretation cells for sections 2 to 4 and section 5 (what step 1 tells the plan, the
  four open questions for the group, a possible `lj_reco_eff_v2`); results added to
  `lj_reco_efficiency.md` (step 1 marked done, questions 3 and 4 added) and to the README.
  Fixed `systematics_study.summed`: it summed weighted histograms over signal points of
  different cross sections, whose effective counts can put more events in a numerator bin than
  in its denominator (Clopper-Pearson then fails); it now sums counts of simulated dark photons
  (new helper `counts`), exact since every event of a signal sample has the same weight
  (checked on all samples). The map titles were moved above the CMS label.
- Started this progress log; AGENTS.md now asks every study in progress to keep one.
- Lepton-jet reconstruction efficiency, step 1: the campaign `lj_reco_eff_v1` was tested
  locally on one file per final state and submitted (422 jobs). `06_lj_reco_efficiency.ipynb`
  and its helpers in `systematics_study.py` (`lj_efficiency_inputs`, `regime_fractions`,
  `plot_lj_efficiency`, `plot_lj_efficiency_map`, `plot_regimes`) were written; the notebook
  takes the opening angle from the per-dark-photon map, since the `genMu_genMu_dR_*`
  histograms are per event. In `lj_reco_efficiency.md`, the mention of the BoostedJPsi
  productions in the user's EOS area and the question about them were removed (they belong to
  a different analysis), and the step 1 go-ahead question was dropped (given by the user).
- Interpretation cells added to 03, 04 and 05 at the user's request (01 and 02 already had
  them); 05 rerun with a robust pull summary (median and central 68%), since a few 4Mu toys
  with broken uncertainties made the mean and width meaningless.
- 4Mu signal-injection failures traced with saved toys: 34 of the 35 toys with both C and D
  empty fail, against 19 of the 41 with only C empty (README open item 2).
- New notebooks `04_systematic_effects` (templates, size maps, impacts, limit breakdown) and
  `05_fit_checks` (signal injection, likelihood scans, correlations), rerun after 03 for every
  new systematic. New scripts `run_combine_impacts.py` (needs minimizer tolerances 100 times
  tighter than Combine's defaults), `run_combine_scan.py` and `run_combine_toys.py`;
  `run_combine_pulls.py` gained `--expect-signal` and correlation matrices;
  `run_combine_limits.add_release_options` is shared by the scripts; new `plot_tools` figure
  functions; `systematics_study.limits_outputs` stops 04 and 05 when 03 is out of date.
- 03 redesigned at the user's request: Brazil grids against the bound state mass and against
  L_xy with the model cross section, comparing "without corrections or systematics"
  (uncorrected weights) with "with corrections and systematics" (`plot_tools.brazil_grid`
  gained `against` and `overlay`; `load_limits` reads xrootd URLs; `upload` clears the
  fsspec-xrootd listing cache, which otherwise breaks reading back).
- AGENTS.md: section "Producing samples: always a Condor campaign".
- Answered which working point, plane and event cuts the limits use (README, "Limits with and
  without the corrections and systematics"); found the `mu_ljs` spread cuts that never reach
  the limits (open question 4).

### 2026-11-13

- Campaign `systematics_v1` (configs, job lists, tarball) run to completion: 3,158 jobs, no
  failures, 200 samples merged, 4 input files without events skipped. Fixed
  `condor/campaign.py skipped`, which could not unpickle chunk outputs (`sidm` not on the
  path).
- `03_limits` built and run: ABCD shape fit in m_LJLJ, with and without the four systematics
  (at most 0.3% apart), one pull plot; its inputs and outputs moved to EOS through xrootd so it
  runs on CoffeaCasa.
- Limit selection and inputs added: `SR_2mu2e`/`SR_4mu` in `selections.yaml`, the
  `abcd6d_*` histograms (with `_iso_up` copies for the jet energy), the `syst` and `variant`
  slicing in `shape_tools`, `abcd_workspace.py --signal-systematics`, `run_combine_pulls.py`,
  the `systematics_v1` Campaign in `datacard_tools`; `systematics.yaml` registry.
- Step 1 of the lepton-jet reconstruction efficiency prepared: `baseNoLj_A_ee`,
  `genA_toE_lxy_vs_daughters_dR`, the collection `lj_reco_effi_systematics`, the campaign
  config `lj_reco_eff_v1`, and the proposal `lj_reco_efficiency.md`.
- AGENTS.md and CLAUDE.md created: the rules for any coding agent in this repository.

### 2026-11-12

- Study started: `01_pileup` and `02_jet_energy` on the campaign `corrections_v1`, and the
  values for luminosity (0.83%), pileup (per point), L1 prefiring (1%) and jet energy (per
  point), with their reasoning in the README.
