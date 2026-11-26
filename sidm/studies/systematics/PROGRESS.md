# Systematics study: progress log

Read this first when picking the study up again, then [`README.md`](README.md) for the
methods and results. Update it as part of every piece of work on the study: the current
state below, and a dated entry in the log (newest first). [`AGENTS.md`](../../../AGENTS.md)
asks for this.

## Current state (2026-11-26, 12:30): step 3 done and pushed (bab7cf0); the advisor's answers implemented

**The advisor's answers (2026-11-26).** (1) The FSR closure: MadGraph+Pythia is the most representative of
the signal simulation, so everything stays as it is (the pixel-veto factor per photon, 4.4% for
`CMS_ljreco_egm`). (2) The size rule makes sense: the whole range sets the size unless there is
meaningful evidence of a dependence on the kinematics. (3) Low priority: look at the stability of the
scale factors against the run period where the statistics allow (the DSA question). (4) An uncertainty
that covers the possible trigger bias is reasonable, as it hardly affects the results.
- **Done (2026-11-26).** Notebook 12, section 3.1: the pixel-seed check and the second generator from
  `fsr_pixel_seed_v1`, run (numbers as below), with the decision. Helpers in `systematics_study.py`
  (`pixel_veto_sf`, `fsr_seed_inputs`, `seed_outcomes`, `seed_passing`, `pixel_veto_factors`) and two
  synthetic tests in `tests/test_systematics_helpers.py`. README, `lj_reco_efficiency.md` (questions 2,
  5 and 6 answered), a comment at `CMS_ljreco_egm` in `systematics.yaml` (the registry is unchanged),
  and the AN (local): the generator choice in one sentence of the FSR paragraph, and the
  two-standard-deviation criterion of the size rule stated in the method.
- **Done (2026-11-26): the DSA factor per run period** (notebook 10, section 7; campaign
  `dsa_prompt_data_v2`, collision data approved 2026-11-26 for this use: 2,580 jobs over the skimmed
  DoubleMuon 2018 A to D, analysis triggers, collection `dsa_prompt_probes_z`, 68 blocks merged, the
  463 skipped files are empty). The skim keeps 99% of the 2018C probes with the same data/DY (0.8840).
  Data/DY with the analysis triggers: 0.891 (A), 0.886 (B), 0.884 (C), 0.881 (D), 0.884 for 2018. In
  common bins (those of 2018B after merging) every period agrees with the year within the factor's
  total uncertainty in all 32 bins (largest 1.3% below |eta| 0.9, 4.7% above 2.1), and 2018C sits at
  the average: the factor of 2018C applies to all of 2018. Helpers: `lj_cr_data(samples=)`,
  `in_period`, `period_spread`, `dsa_sf_groups`, `dsa_scale_factors(trigger=, groups=)`.
- **Deck and message (local):** `slides/dsa_periods_2026-11-26.tex` (9 slides) and
  `slides/dsa_periods_2026-11-26_mattermost.md`. The AN (local) has one sentence on the periods in the
  DSA paragraph.
- **Next:** the commit of all of this (asked for 2026-11-26), then the user pushes; the PR
  description is updated in chat.

## State of 2026-11-25, 17:30

**The advisor's answer on step 3 (2026-11-24):** proceed with the corrections and systematics,
taking the value from the nearest populated bin, with scale factors **per constituent** and their
uncertainties, then evaluate the lepton-jet reconstruction and selection with the corrections
applied, and assign an additional systematic if needed. He also asked whether the points of one bin
are offset in x in the deck's figures: they were, on purpose (10 to 12% of the bin width); the bins
are the same, and every point now sits at its bin centre (`CURVE_STYLES`).
- **Done.** The DSA-muon scale factor per DSA muon (notebook 10, section 6,
  `sidm/data/corrections/dsa_muon_sf_2018.json`, `dsa_muon_id` in the default corrections);
  `Histogram(entry_factors=...)` and the `_sf` probe histograms (collection `lj_reco_closure`);
  campaign `lj_reco_cr_v4` (1,699 jobs, no file skipped, 18 samples merged);
  `12_lj_corrected_closure.ipynb` run and interpreted; `CMS_ljreco_mu` (1.7% per muon-type lepton
  jet) and `CMS_ljreco_egm` (4.4% per egamma-type) in `systematics.yaml` (kind `per_lj`);
  `systematics_v4` (signal, 422 jobs); notebooks 03, 04, 05 rerun on it with ten systematics.
- **Results.** With every scale factor applied, data/simulation is 1.009 ± 0.006 for J/ψ lepton
  jets above 50 GeV (1.017 ± 0.003 at displaced vertices), 1.008 ± 0.017 for J/ψ pairs with two
  DSA muons, 1.044 ± 0.041 for FSR photons above 30 GeV. The corrections lower the signal in region
  A by a median 18% (2Mu2E) and 26% (4Mu), most of it the DSA factor; the expected limits rise by
  23% and 37%; the two lepton-jet terms add about 0.5%. Exclusion: 3 and 1 of 48 integrated, 14 and
  17 binned.
- **Checks on a few simulated files (scratch, not committed):** a quarter of the passing J/ψ
  lepton jets above 50 GeV are one PF plus one DSA muon (so the J/ψ closure tests the DSA factor
  inside lepton jets); in the signal a DSA muon of a selected muon-type lepton jet is nearly always
  needed for its two muons (one besides two PF muons in 0.4 to 1.1%). `ljs.muons` holds PF and DSA
  muons together; PF muons alone are `ljs.pfMuons`.
- **Waits on the group:** the size rule (nearest regions regardless would give up to 8.3% and 12%);
  2018A/D for the DSA factor; the trigger part of its uncertainty.
- **Done (2026-11-25): `fsr_pixel_seed_v1`** (see the log; analysed with scratch scripts, not
  yet in a notebook). The pixel seed: 11.5% of the FSR photons above 30 GeV have one; they still
  pass 42% (MiNNLO+PHOTOS) or 32% (MadGraph+Pythia) of the time, always through an electron-only
  lepton jet, so the veto factor applied as an efficiency reweighting lowers Drell-Yan by 2.6 to
  3.0% instead of the 4.9 to 5.0% of the per-photon factor (data/DY x0.976 to 0.979). The
  generator matters more: without scale factors the M50 samples pass 0.894 ± 0.010 (MiNNLO) and
  0.829 ± 0.024 (MadGraph); data/DY above 30 GeV, no SF / SF per photon / veto reweighted, is
  0.991 / 1.044 / 1.022 (MadGraph) and 0.924 / 0.973 / 0.949 (MiNNLO M50 with MadGraph M10to50).
  The 4.4% of `CMS_ljreco_egm` covers both generators as the factors are applied now; with the
  veto reweighted the envelope would be 5.1%.
- **Waits on the user and the advisor:** the follow-up on question 2 (drafted in chat); whether to
  keep the per-photon pixel-veto factor and 4.4% (proposed) or move to the reweighting.
- **Next:** on the answer, put the analysis into notebook 12 and `systematics_study.py` (and the
  README), then commit the new histogram, its test and the campaign config.

## State of 2026-11-24, 14:30

**The review of PR #90 and the feedback on step 2 are done and committed, not yet pushed.**
Part 1 (the limits on the analysis pair) and part 2 (version 3 of the control samples, in
simulation and in data) are finished, and the follow-ups to the reply on step 2 (2026-11-24) are
done: the DSA muons of prompt muons in data and the coverage map. The user pushes, answers the two
review comments, opens the MET-filter PR, and shares the summary deck.
- **Part 1, done.** `systematics_v3` (3,158 chunks, 4 empty files skipped) and
  `systematics_data_v2` (2,580 chunks, 463 empty files skipped) merged; notebooks 03, 04, 05 rerun
  on them with rewritten interpretations; README and AN Sec. 7 (local) updated. Data B, C, D:
  2mu2e 80, 17, 30 (prediction 45.3), 4mu 42, 5, 3 (prediction 70); 79, 19, 28 and 45, 12, 5 on the
  leading pair. The limits are a median 11% (2Mu2E) and 15% (4Mu) lower than on the leading pair,
  mostly from the background; the signal in region A is a median 1.7% and 4.9% lower (up to 26% for
  4Mu at 5 GeV and 300 cm), from the pair selection of PR #71: the axes alone change region A by at
  most 0.15% at eight points (three at 5 GeV, 150 to 300 cm; scratch script
  `abcd_pair_compare_5gev.py`). The 4Mu injection test has more toys near zero (31% at the limit),
  which a counting model of the four regions reproduces (C and D hold 5 and 3 events).
  `sidm/studies/sr_mljlj` (untracked) rerun too.
- **MET filters:** branch `passFlagsInBase` (`28a34f9`, off upstream `main`) for a separate PR,
  not pushed yet; its description is ready for the user. Whichever of it and PR #90 merges second
  drops the explicit `"pass flags"` of `SR_2mu2e`/`SR_4mu`.
- **Part 2, simulation done.** `lj_reco_cr_v3` complete (1,927 chunks, 114 samples merged at
  08:30; 2 empty files skipped). Ten signal files of the v10 production fail every attempt: each
  has a damaged basket (lzma data error) in at least one muon-vertex branch (`PatMuonVertex_*`,
  `PatDSAMuonVertex_*`, `DSAMuonVertex_*`), four of them in other branches too (`LHEPdfWeight`,
  `TrigObj_pt`, ...), found by reading every branch of each file. Only the v3 histograms read
  those branches, so `systematics_v3` read the same files without error. The campaign tool
  dropped them for good (`condor/campaigns/lj_reco_cr_v3/bad_files.txt`); the merge normalizes
  each sample to the events it processed, so this costs statistics only. Notebooks 07 (section 6)
  and 09 run with interpretations; README and `lj_reco_efficiency.md` (step 2, version 3;
  question 5 on DSA muons in data) updated.
- **Safety of version 3:** every region of the FSR Z sample is safe; every J/ψ vertex bin below
  15 cm is safe; the bin above 15 cm is a caution for shape uses (4Mu, 800 GeV, 5 GeV, 5 mm: 4.1
  events against 9,300, binned Z 1.28) and is never read (`lj_cr_data` zeroes it).
- **Part 2, data done.** `lj_reco_cr_data_v2` complete (2,580 chunks, the 463 empty files
  skipped, 68 samples merged at 11:45; the watcher moved 12 slow jobs, which then ran 17 to 27
  minutes). Notebook 08, section 6, run and interpreted. J/ψ above 50 GeV: the passing fraction in
  data is 0.979 to 0.985 at every vertex L_xy from 0 to 15 cm, data/simulation 0.995 ± 0.003
  above 1 mm; below 0.5 mm the simulated passing fraction exceeds one (few simulated prompt
  J/ψ), which the integrated ratio of v2 (0.988 ± 0.006) contains. Converted photons (FSR
  selection, above 30 GeV): data/Drell-Yan 0.94 ± 0.10, limited by the simulated statistics;
  unconverted 1.02 ± 0.04 (0.86 ± 0.09 above 60 GeV). `efficiency_ratio` now gives a positive
  uncertainty where the subtraction leaves a negative passing fraction (below the 30 GeV
  threshold), which had stopped the plot. Sections 1 to 5 print the same as before.
- **Documents:** README, `lj_reco_efficiency.md`, AN Sec. 7 (local: the J/ψ displacement with a
  new figure, the egamma comparison rewritten with the FSR selection for converted and
  unconverted photons) and the summary deck `slides/summary_2026-11-24.tex` (24 pages, local)
  updated.

- **Follow-ups to the reply on step 2 (2026-11-24).** Points 1 to 3 (J/ψ production radius and
  vertex, converted photons, closure) were already done in version 3. Point 4, DSA muons in data:
  approved for 2018C; the safety check of both DSA regions (notebook 09, section 5) found every
  region safe; campaign `dsa_prompt_data_v1` (636 jobs over the unskimmed DoubleMuon 2018C,
  `request_memory: 6GB`, since whole unskimmed files pass 2 GB; jobs that started at 2 GB were
  released with 6 GB) merged at 13:40; notebook 10: data/simulation 0.862 per prompt muon with
  Mu37_TkMu27 (0.80 above 120 GeV, 0.77 above |η| = 2.1; 0.87 to 0.88 with the other triggers),
  0.77 ± 0.01 for J/ψ pairs above 40 GeV with the analysis triggers, the square of the single-muon
  ratio. Point 5: the coverage map (notebook 11, from `lj_reco_eff_v1`): over the limit grid a
  median 97 to 99.7% of the muon-type lepton jets lie below ΔR = 0.03, 65 to 68% beyond 15 cm, 28
  to 29% above 300 GeV; 97% of the egamma-type ones above 120 GeV. README, `lj_reco_efficiency.md`,
  AN Sec. 7 (local: a DSA paragraph and figure) and the deck (25 pages, local) updated.
- Unskimmed DoubleMuon ntuples are on EOS for 2018A (2,762 files) and, in blocks, 2018D (not in
  `data.yaml`); 2018C has 636 of its 637 files.

**Next:** the user pushes `systematicsStudy` and `passFlagsInBase` (replies, PR texts and a
Mattermost summary were given in chat); step 3 (scale factors, with the DSA muons per muon) waits
on the group's answers in `lj_reco_efficiency.md` ("What is needed").

## State of 2026-11-17 (step 2 for the prompt regime)

- **Simulation:** campaign `lj_reco_cr_v2`, notebook 07: every region of both control samples is
  safe (with the sparse-bin rule, switch `STRICT_SPARSE_BINS`).
- **Data:** campaign `lj_reco_cr_data_v1` complete (2,590 chunks; 463 empty files skipped; 68
  samples merged). Notebook 08, read through the safety table: muon-type (J/ψ, above 50 GeV)
  data/simulation 0.988 ± 0.006 (endcaps 0.976 ± 0.008); egamma-type (Z → μμγ, unconverted
  photons above 30 GeV) 1.01 ± 0.04, 0.80 ± 0.08 above 60 GeV (112 probes).
- **AN Sec. 7** holds the measurement (two figures in `an_systematics/`); the deck
  `slides/systematics.tex` (38 pages, local) presents the study only: no AN updates and no
  references to the advisor (the user's rules).
- **Upstream `main` merged** (2026-11-17), with its lepton-jet pair selection (S. Abbott): the
  control-sample histograms are unchanged by it, the signal-region yields are not (−1.3% for 2Mu2E
  500/1.2/30 cm, −3.3% for 4Mu), which `systematics_v3` and `systematics_data_v2` take in.

**The user's workflow (2026-11-16):** after running anything, update the AN (in the advisor's
style, `an_systematics/WRITING_GUIDE.md`), then the summary deck (the study's results, without
the AN updates), then commit and paste an updated PR description in chat. Keep `slides/` and `an_systematics/` local.

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

### 2026-11-26

- (12:30) `dsa_prompt_data_v2` submitted 10:00 (cluster 86619898, lpcschedd6), the nine known-slow nodes
  excluded at submission and the 30-minute mover running; done 11:44. The last 21 jobs sat on a new slow
  node, cmswn2246 (58 min, period medians 20 to 31 min), and were moved by hand. Merged (68 blocks).
  Notebook 10, section 7: a first comparison against the factor's systematic alone, with each period's
  own merging of sparse bins, flagged 3 of 32 bins; two were sparse bins within the factor's total
  uncertainty, and one compared differently merged bins (80-200 GeV of 2018B against 120-200 GeV of
  2018). With common bins and the total uncertainty, which is what the correction applies, no bin
  departs. The rerun rewrites `dsa_muon_sf_2018.json` byte for byte.
- Advisor's answers (above). Notebook 12, section 3.1 added and run: it reproduces the scratch
  analysis of 2026-11-25 (seeded 11.2% / 11.5%, recovered through electrons 32% / 42%, veto factor
  per photon x0.9505 / x0.9509 against x0.9705 / x0.9744 reweighted, M10to50 -77.0 +- 24.2 probes
  from 22 effective events of mean weight 9.4, data/DY 0.991 / 1.044 / 1.022 and 0.924 / 0.973 /
  0.949). The pixel-veto factors are read through `Corrections` (0.9623, 0.8833, the same from CVMFS
  and the CI copy). The scratchpad was cleared again overnight; nothing in use depended on it.
- `dsa_prompt_data_v2` prepared (see above); its submission was blocked by the permission check and
  waits for the user. The 2018C safety table already covers it: notebook 09 judged the DSA probes at
  the luminosity of all of 2018 (64 million background events), with the trigger axes kept.

### 2026-11-25

- (18:30) Deck `slides/fsr_pixel_seed_2026-11-25.tex` (10 slides, local like every deck): seeded
  FSR photons, the veto per photon against reweighted, the two generators against the data, the
  M10to50 sample in the subtraction, data/Drell-Yan for every prediction, the egamma-type size and
  open points. Figures in `slides/figures_fsr_seed/`, drawn at slide size by the scratch script
  `fsr_deck_figures.py`; it shows control-sample data, so it stays off the MC-only EOS area.
- (17:30) `fsr_pixel_seed_v1` analysed. The last chunk (`DYJetsToLL_M50` 125) sat 3 h 10 min
  on cmswn2245, a known slow node that was not excluded at submission (sample median 9 min);
  moved at 16:34, done 16:46; merged 16:51, no file skipped. Checks: the campaign's
  `DYJetsToLL_M50` reproduces `lj_reco_cr_v4` exactly (0.8288 ± 0.0238 without factors), and the
  seed histogram gives the closure histogram's efficiencies (the m_mmg bins are the peak and
  sidebands). Gen-matched FSR probes above 30 GeV, MiNNLO+PHOTOS / MadGraph+Pythia: seeded 11.5% /
  11.2%; seeded pass 0.422 ± 0.024 / 0.321 ± 0.055, all through electron-only lepton jets;
  unseeded probes that pass the other photon cuts always pass; the other photon cuts pass 91.0% /
  88.4%. The veto reweighted (pass x SF, fail x (1 - SF e)/(1 - e), e the simulated fraction
  without a seed, 0.94 barrel and 0.72 endcap) against per photon, all probes after the
  subtraction: 2.6% / 3.0% against 4.9% / 5.0%. `DYJetsToLL_M10to50` is about 22 events of weight
  about 10, mostly in the upper sideband: -77 ± 24 probes at 0.19 ± 0.10 against 1,482 for
  `DYJetsToLL_M50`, which lowers data/DY by about 4%. To add it to `DYJetsToMuMu_M50`, the latter
  goes through `corrections_study.drop_genweight_outliers` (the 1.45e14 event, x632). Data/DY
  above 30 GeV, no SF / SF / reweighted: MadGraph M50 only 1.033 / 1.088 / 1.066; MadGraph
  M10to50 + M50 0.991 ± 0.039 / 1.044 ± 0.041 / 1.022; MiNNLO M50 only 0.958 ± 0.022 / 1.008 ±
  0.023 / 0.984; MiNNLO M50 + MadGraph M10to50 0.924 ± 0.025 / 0.973 ± 0.026 / 0.949. Scratch
  scripts (scratchpad, may be cleared): `fsr_seed_analysis.py`, `fsr_seed_closure_reweighted.py`,
  `fsr_closure_by_generator.py`.
- (13:15) Campaign `fsr_pixel_seed_v1` submitted (1,412 jobs, cluster 5088141, lpcschedd4;
  simulation only): `DYJetsToLL_M50` and `DYJetsToMuMu_M50` (independent, about three times the
  skimmed events, PHOTOS), channel `lj_cr_Z_to_mumugamma_alone`, collection `lj_reco_pixel_seed`:
  the new `lj_cr_zmmg_fsr_probes_alone_seed` (`zmmg_probe_seed_hist`: the probe's pixel seed,
  whether it passes the other lepton-jet photon cuts, a selected electron within dR < 0.025, and
  the photons and electrons of the lepton jet that makes it pass, `lj_composition`) and the FSR
  closure histogram with the scale factors (a second Drell-Yan sample for the closure).
  `zmmg_probe_gen_fsr` now also accepts a photon whose mother is the Z: PHOTOS attaches the
  radiation to the decaying Z, Pythia to the muon (checked on one file of each). Tested locally
  on 46 files (fills, no warnings); unit test for `lj_composition`; code tarball remade.
- First look from `lj_reco_cr_v3` (Drell-Yan, gen-matched FSR probes above 30 GeV): 11.6% have
  a pixel seed; with a seed they still pass 32% of the time (n_eff 72), without 88% (n_eff 577).
  About a third of the photons that fail the veto are recovered, so the per-photon pixel-veto
  factor (about 0.95 on these photons) would over-correct by about 2% of its 5%.
- Advisor's follow-up on step 3 (reply drafted in chat, not yet sent): (2) why the EGM scale
  factors do not improve the FSR closure (0.991 ± 0.039 without, 1.044 ± 0.041 with; mostly the
  pixel-seed veto, EGM 2018 Loose 0.962 barrel, 0.883 endcap, inclusive in R9 as applied): within
  the 4% precision of the test, plus the recovery above; (3) the derivation of the additional
  systematic, explained step by step. The scratchpad was cleared overnight (old drafts, the
  watcher, the unused AN figures).

### 2026-11-24

- (23:00) Step 3 finished: `lj_reco_cr_v4` merged (no file skipped), notebook 12 run and
  interpreted, `CMS_ljreco_mu`/`CMS_ljreco_egm` registered, notebooks 03, 04 and 05 rerun on
  `systematics_v4` with ten systematics and their interpretations rewritten; README,
  `lj_reco_efficiency.md`, AN Sec. 7 (local), the deck `slides/lj_corrections_2026-11-24.tex`
  (results, questions) updated; a variable in `nearest_factors` renamed for pylint;
  `systematics_study.acts_in`, so that notebook 05 skips `CMS_ljreco_egm` in 4mu (a row of ones,
  which Combine leaves out of the workspace; the scan failed on it); committed.
- (17:30) Step 3 code and campaigns: the DSA scale factor, the `sf` axis, `lj_reco_cr_v4` and
  `systematics_v4` submitted; the plots without x offsets.
- (14:30) Follow-ups to the reply on step 2: DSA safety check (notebook 09, section 5),
  `dsa_prompt_data_v1` on the unskimmed 2018C (approved), notebook 10 (DSA data/simulation 0.862
  per muon, 0.77 per J/ψ pair), notebook 11 (the coverage map). `integrated_efficiency` takes a
  `passing` argument, `plot_mass_data_mc` a `scale`; `LJ_COVERAGE`, `fraction_inside` and
  `coverage_fractions` are new. README, `lj_reco_efficiency.md`, AN Sec. 7 and the deck updated;
  committed.
- (12:20) `lj_reco_cr_data_v2` complete and merged; notebook 08 section 6 run with interpretations
  (J/ψ flat in the vertex displacement out to 15 cm, data/simulation 0.995 ± 0.003 above 1 mm;
  converted photons data/Drell-Yan 0.94 ± 0.10; unconverted 1.02 ± 0.04). `efficiency_ratio`
  fixed for negative passing fractions. README, `lj_reco_efficiency.md`, AN Sec. 7 and the deck
  (two data slides) updated; committed.
- (09:00) `lj_reco_cr_v3` complete and merged (the last QCD chunk ran 175 minutes on a normal
  node after a move off a slow one). Notebooks 07 and 09 run and interpreted: the sideband
  subtraction closes for the J/ψ probes above ΔR = 0.03 and, with the FSR selection, for the
  photons (0.06 high on average); the FSR selection removes QCD from the Z peak and raises the
  purity of converted photons to 0.65 to 0.82; the J/ψ passing fraction is flat in the production
  radius and the vertex displacement out to 5 cm (1.00 to 0.98); 81% of prompt muons have a
  selected DSA muon, both muons of a boosted J/ψ 44 to 52%, and the analysis triggers bias that
  up. Safety: all safe except the J/ψ vertex bin above 15 cm. `lj_reco_cr_data_v2` submitted.
  Deck: injection slide redrawn as one legible panel, `\pt` defined.
- (06:20) Part 1 of the review: `systematics_v3` and `systematics_data_v2` merged, notebooks 03,
  04 and 05 rerun and their interpretations rewritten, README and AN Sec. 7 (local) updated; the
  v3/v2 comparison and the axes-only test at three 5 GeV points (no change in any region) done.
  `lj_reco_cr_v3`: 9 failed signal chunks retried without their damaged files, slow QCD chunks
  moved, skipped files resubmitted.

### 2026-11-23

- (22:40) Review of PR #90 (two comments) and the feedback on step 2: fixed the ABCD axes
  (`43e9ae1`), moved the MET filters into the base selection on a separate branch (`28a34f9`),
  wrote the v3 control samples and the DSA study (`1aa8f3a`), and submitted `systematics_v3`,
  `systematics_data_v2` and `lj_reco_cr_v3`. Notebooks 07 (section 6) and 09 written, not run.

### 2026-11-17

- (15:10) The limit breakdown figure was unreadable on its slide: a 4:1 strip of columns with
  the nuisance names, its scale set by the 132% point at the shortest lifetime.
  `plot_tools.breakdown_plot` now draws one row per left-out nuisance with readable names, 2Mu2E
  and 4Mu in one panel (a colour-blind-safe blue and orange) and the medians as black ticks;
  notebook 04 plots the limit grid only (the shortest lifetime stays in its table). Notebook 04
  rerun: every printed number unchanged. The deck no longer refers to the advisor (the
  sparse-bin rule left its questions for the group). AGENTS.md: agents do not list themselves
  as authors or contributors (the user's rule).
- (13:25) At the user's request the deck no longer covers the AN: the two AN Sec. 7 slides and
  the AN lines of the title, summary, what-changed and step-2 slides are gone (38 pages). The
  caveats of the removed "kept out of the AN" slide were already on the study slides; the
  AN Sec. 8 control-region counts moved to the background slide. The backup slides were brought
  up to date with the merge and notebook 08. Rule added to `an_systematics/WRITING_GUIDE.md`.
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
