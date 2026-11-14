# Systematics study: progress log

Read this first when picking the study up again, then [`README.md`](README.md) for the
methods and results. Update it as part of every piece of work on the study: the current
state below, and a dated entry in the log (newest first). [`AGENTS.md`](../../../AGENTS.md)
asks for this.

## Current state (2026-11-14, 13:00)

**Running:** nothing. Step 1 of the lepton-jet reconstruction efficiency
([`lj_reco_efficiency.md`](lj_reco_efficiency.md)) is done: the campaign `lj_reco_eff_v1` is
complete (422 of 422 chunks, 180 samples merged) and
[`06_lj_reco_efficiency.ipynb`](06_lj_reco_efficiency.ipynb) is executed with its
interpretations and its section 5 (what step 1 tells the plan, and the questions left open).

**Next steps**
1. Wait on the group for the lepton-jet efficiency (open question 2 below, four parts); steps 2
   and 3 of the plan start from their answers.
2. Offered to the user, not done: a campaign `lj_reco_eff_v2` that splits the muon-type lepton
   jets by constituent (PF, DSA, both; the channels exist), to see which control sample covers
   the decays at 10 to 100 cm (section 5 of 06).
3. The rest waits on the user or the group (open questions below).

**Notebooks** (inputs and outputs on EOS, read and written through xrootd)

| notebook | state | input | outputs |
|---|---|---|---|
| `01_pileup` | executed 2026-11-12 | campaign `corrections_v1` | in the notebook |
| `02_jet_energy` | executed 2026-11-12 | campaign `corrections_v1` | in the notebook |
| `03_limits` | executed 2026-11-14, with interpretations | campaign `systematics_v1` | `.../systematics_v1_limits/` |
| `04_systematic_effects` | executed 2026-11-14, with interpretations | outputs of 03 | `.../systematics_v1_limits/effects/` |
| `05_fit_checks` | executed 2026-11-14, with interpretations | outputs of 03 | `.../systematics_v1_limits/fit_checks/` |
| `06_lj_reco_efficiency` | executed 2026-11-14, with interpretations | campaign `lj_reco_eff_v1` | `.../lj_reco_eff_v1_results/` |

`...` is `root://cmseos.fnal.gov//store/user/scampbel/sidm_coffea_outputs`. After a change to
`systematics.yaml`, rerun 03, then 04 and 05; 04 and 05 stop if 03 is out of date.

**Campaigns** (configs in `condor/campaign_configs/`, run with `condor/campaign.py`, working
state in `condor/campaigns/<name>/`)
- `corrections_v1`: input of 01 and 02 (`.../corrections_v1`).
- `systematics_v1`: complete on 2026-11-13 (3,158 jobs, 200 samples merged in
  `.../systematics_v1`); input of 03.
- `lj_reco_eff_v1`: complete on 2026-11-14 (422 jobs, 180 samples merged in
  `.../lj_reco_eff_v1`); input of 06.

**Registry** (`systematics.yaml`): `lumi_13TeV` flat 0.83%, `CMS_pileup` from the pileup
weight variations, `CMS_l1_prefiring` flat 1%, `CMS_jet_energy` from the isolation scaled by
1.10.

**Results to keep in mind** (details in the README and in the notebooks' interpretations)
- With corrections and systematics the expected limits are 2.0% (2Mu2E) and 1.3% (4Mu) higher
  than without; the exclusion is the same (13 of 48 and 24 of 48 points).
- The four systematics do not change the expected limits; each moves the fitted r by about its
  own size, 1% or less.
- The signal-injection test validates the 2Mu2E fit but not the 4Mu one: regions C and D are
  often empty in 4mu (README open item 2).
- The 2mu2e simulated background in B and C rests on about 1.5 effective QCD events each, so
  the simulation-only 2Mu2E limits test the machinery, not the sensitivity.
- Lepton-jet efficiency, step 1: on the limit grid only the shortest lifetime decays mostly
  inside the pixel layers; the other three take 41 to 73% of their lepton jets from 10 to
  100 cm. The 2Mu2E muon-type efficiency (0.97 to 0.99) is conditioned by the base selection,
  which already requires both muons; the 4Mu one (0.91 below 10 cm, 0.82 from 10 to 100 cm,
  0.33 beyond 300 cm) shows the size. Egamma-type: 0.89, 0.80, zero beyond 140 cm.

**Open questions for the user or the group**
1. The 4mu working point or method, given the empty C and D (the 3-pixel-hit working point
   and the transfer factor from the validation regions are under discussion): group.
2. Lepton-jet reconstruction efficiency (`lj_reco_efficiency.md`, "What is needed"): which data
   control samples may be read (cosmics, Z → μμγ, J/ψ → μμ, J/ψ → ee), scale factors per
   lepton jet or per constituent, the provisional value until step 2, and (new from step 1)
   how to set the uncertainty of the most collimated muon pairs (m_Zd = 0.25 GeV at 800 and
   1000 GeV), which no control sample reaches: group.
3. Whether 03 should carry the agreed 20% per-bin shape uncertainty of the region-A prediction
   for 2mu2e (`abcd_workspace.py --shape-sys 0.2`) and the agreed 4mu bin merging: offered to
   the user, not done.
4. The `mu_ljs` cuts of `SR_2mu2e`/`SR_4mu` (the dz and dxy spread cuts) never reach the
   limits: they filter only the `mu_ljs` view, while the category cut and the `abcd6d_*`
   histograms read `ljs`. A fix needs a decision (drop the lepton jet or veto the event) and a
   new campaign: told the user, not changed.
5. The HEM electron veto and the `mu_veto` of the earlier limit selection are not defined in
   the repository.
6. Write access to the `lpcmetx` group area: until then merged outputs stay in the user area.

**Committed** (2026-11-14, branch `systematicsStudy`, for a pull request to `main`): the study
folder, the framework and limit-tool changes it needs, the campaign configs, and AGENTS.md and
CLAUDE.md. The branch sits on the campaign-configs commit `430542a`, not yet in `main`. Left
out on purpose (AGENTS.md): job lists, filelists, logs, `condor/campaigns/`, outputs, and the
local notes in `sidm/studies/limit_plotting/`. Later changes: stage by explicit path only.

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

### 2026-11-14

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
