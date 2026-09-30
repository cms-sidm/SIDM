# Truth kinematics for the analysis note

Generator-level characterization of the 2018 v10 signal samples for the
`signal_kinematics` section of AN-23-107: proper-lifetime validation, lepton
impact parameter, boost, dark-photon polarization, reconstruction migration,
and the production gen-filter story. Three executed notebooks with their
figures in `figures/` (vector PDF, CMS Simulation style; the repo-wide
`*.pdf` ignore keeps them untracked; they regenerate on notebook execution,
render inline in the committed notebooks, and are copied into the analysis
note repository).

## The canonical output

All notebooks read one production:
`/store/group/lpcmetx/SIDM/coffea_outputs/murtazas/truth_kinematics_forAN/`
(180 per-sample `.coffea` files + metadata sidecars; full statistics, every
file of every sample; unweighted histograms). Load it with
`_lifetime_refit.load_truthkin()`, which caches locally.

Channels, defined in `sidm/configs/selections.yaml`:

- `genOnly`: status-1 generator leptons as the only object definitions;
  **no event cuts at all** (no trigger, no PV filter, no vertex requirements).
- `genOnly_born`: the status-23 variant, for rest-frame comparisons.
- `genFilterEmulation` / `genFilterEmulation_isHardProcess`: re-apply the
  central production gen-filter cut string (two statusFlags variants). Their
  cutflows are re-pass rates on already-filtered events (0.84-1.00 across the
  grid, median ~0.97) that validate the filter form, never a filter
  efficiency (see below).
- `baseNoLj_noTrigger` / `baseNoLjNoLjsource_noTrigger`: analysis LJ-source
  object definitions with/without quality cuts, PV filter, **no HLT**; used
  for the reconstruction-migration maps.

Hist collections: `genA_lifetime`, `gen_truth` (new: the pure-generator subset
of `genE_base`+`genMu_base`, plus `fromA` pT/dxy variants), 
`genBS_genA_kinematics`, `genA_base`, `lepton_genA_base`.

Produced with the standard condor pipeline: `condor/submit_truthkin.sub`
plus job arguments regenerated from the committed sample lists with
`condor/make_job_args.py` (`--files-per-job 5 --replace-xcache`, once per
channel YAML; 1833 jobs), merged with
`sidm/scripts/merge_coffea_chunks_eos.py`. The Dask scripts
`_truthkin_gen_run.py` / `_recomigration_run.py` document the same
configuration and serve for small reruns; a full-grid Dask driver does not
survive the login node.

## The anatomy extension (v3)

A second production for the final-state-anatomy and trigger figures:
`/store/group/lpcmetx/SIDM/coffea_outputs/murtazas/truth_kinematics_forAN/anatomy_v3/`
(same format and statistics as v2: 180 merged per-sample `.coffea` + sidecars,
every file, unweighted). Load with
`_anatomy_lib.load_v3()` (local cache, `ANATOMY_CACHE` env to relocate).

Channels: `genOnly`, **`genOnly_trigger`** (new: identical generator-level
object definitions plus one event cut, the analysis dimuon-trigger OR, so the
ratio of a per-event histogram to its `genOnly` counterpart is the
absolute signal trigger efficiency vs that variable (per-object
histograms fold in the partner pair in 4Mu, so the trigger notebook
uses the 2Mu2E samples for per-dark-photon efficiencies)), `baseNoLj_noTrigger`, `baseNoLj`, and
`base` (the reconstruction-level figures state their selections explicitly).
This production predates the `>=2 muons (pf or dsa) pT > 26 GeV` cut that
`base` has since acquired, so its `base` channel corresponds to the current
`base_noMuPtCut`; regenerate the reconstruction-level figures against that
channel to reproduce them.

Hist collections: `gen_truth`, `genA_lifetime`, `genBS_genA_kinematics`,
`genA_base`, `genA_ratio`, `lj_base`, `lj_lj_base`, `muon_base`, and
`anatomy_extra` (new: per-dark-photon daughter-pair `dR` on linear and log
axes and the log-log `dR` vs `pT` map, built from the recorded decay
daughters so 4Mu pairing is unambiguous; status-1 four-lepton invariant
masses `gen4Mu_invmass`/`gen2Mu2E_invmass`; leading/sub-leading PF and DSA
muon `pT`; the `mulj_egmlj_invmass` reco pair mass). Known-empty in this
production: the two `genA_*Lj_lxyRatio` histograms of `genA_ratio`
fail to fill in every job (`AttributeError: no field named 'kinvtx'`
-- they reference an LJ kinematic-vertex field the current LJ
reconstruction does not produce; pre-existing upstream, flagged for
an upstream report).

Produced with `condor/submit_truthkin_v3.sub` and
`condor/job_args_truthkin_v3.txt` regenerated exactly as for v2 (same sample
lists, `--files-per-job 5 --replace-xcache`, once per channel YAML, then
concatenated; 1833 jobs; nine chunks failed on transient xrootd
errors and were resubmitted with a copy of the submit file pointed at
the failing (sample, chunk) lines), merged with the same
`sidm/scripts/merge_coffea_chunks_eos.py` invocation into `anatomy_v3/`.

Two follow-on productions add histograms v3 did not book, written to
`/store/group/lpcmetx/SIDM/coffea_outputs/murtazas/truth_kinematics_forAN/anatomy_v3_cond/`
and `.../anatomy_v3_selfcon/`. The displacement-conditioned run covers all 90
2Mu2E samples in channels `genOnly` and `genOnly_trigger` with the
`trigger_conditioned` collection, whose one histogram
`genA_toMu_lxy_vs_daughters_dR` holds the joint distribution of dark-photon
decay radius and daughter opening angle (both axes filled from the same
decayed-daughter selection, so an entry is one dark photon and the pair of
axes is a genuine joint distribution) alongside the one-dimensional `dR` and
`lxy` histograms it is checked against. Produced with
`condor/submit_conditioned.sub` and `condor/job_args_conditioned.txt`
(`--files-per-job 10`; 435 jobs), merged with the same
`sidm/scripts/merge_coffea_chunks_eos.py` invocation into `anatomy_v3_cond/`.

The reconstruction self-consistency run covers the 12 mid-lifetime
`m_Zd = 1.2 GeV` samples, six bound-state masses in each channel, in channel
`base` with the `selfconsistency_matched` collection:
`genA_matched_lj_lj_invmass` and `genA_matched_mulj_egmlj_invmass` pair each
dark photon with its own nearest lepton jet within `dR = 0.4`, against which
the blind leading-pT pair masses `lj_lj_invmass` and `mulj_egmlj_invmass`
already in v3 are compared. Produced with `condor/submit_selfcon.sub` and
`condor/job_args_selfcon.txt` (the per-channel `job_args_selfcon_4mu.txt` and
`job_args_selfcon_2mu2e.txt` concatenated; `--files-per-job 20`; 37 jobs),
merged with the same script into `anatomy_v3_selfcon/`.

## Notebooks

- `lifetime_forAN.ipynb`: proper-lifetime faithfulness and the
  acceptance-corrected recovery, at ~50x the statistics of the original
  `lifetime_study` pass. Where the lab cap leaves the distribution intact
  (108/180 samples) the corrected fit closes within 0.5% of nominal; over the
  full grid the median is 1.002 with 119/180 within 5%, the 36 points beyond
  +10% being the longest-ctau sample of each mass point (the single-radius
  approximation to the filter's rho/z cylinder).
  Fitted lab cap R_max = 814 cm over the heavily-truncated samples; the
  truncation-onset test covers all 36 mass points.
- `truth_kinematics_forAN.ipynb`: gen lepton |d0| across the lifetime grid
  (the variable the NoVtx triggers and displaced reconstruction respond to),
  the boost map across the mass grid, and the gen-filter efficiency figure.
- `polarization_migration_forAN.ipynb`: the polarization-fit summary
  over all five m_XX points (transverse alpha ~ 1 at M_Zd >= 1.2 GeV; muon
  velocity suppression at M_Zd = 0.25 GeV, directional only above m_XX = 500
  GeV, where the rest-frame angle is numerically ill-conditioned at the boost)
  and the reconstruction-migration maps
  (electron-to-photon and PF-to-DSA handoff vs Lxy).
- `final_state_anatomy_forAN.ipynb`: the first-principles anatomy, from v3:
  production system, back-to-back topology and the Zd momentum scale,
  gen and reco self-consistency masses (m(4l) and m(LJ,LJ) at m_Bs, pair
  masses at m_Zd, LJ/Zd pT response at unity), the gen-matched pair-mass
  companion that separates the pairing from the matching requirement (the
  blind leading-pT width reaches 0.373 of the median at m_Bs = 1000 GeV in
  4Mu against 0.134 gen-matched, with combinatorics capped at 1.9%), the
  collimation scan with the
  2 m/pT law and the grid median-dR map, lepton pT vs the mass scan, the
  pT-asymmetry vs cos(theta*) relation with the muon-velocity floor, and the
  displacement ladder of the display corners.
- `trigger_context_forAN.ipynb`: the signal muons against the four L2 NoVtx
  dimuon paths, from v3: sub-leading-muon spectra vs the 23/25 GeV
  thresholds, efficiency turn-ons from the stored HLT bits
  (`genOnly_trigger`/`genOnly` ratios), efficiency vs displacement and pair
  opening angle, both conditioned on the other from the joint
  lxy-vs-dR histogram (efficiency halves by lxy = 95 cm for the tightest
  pairs against 205 cm for the widest, and the best opening angle moves out
  to dR = 0.13 beyond lxy = 150 cm), median-efficiency grid maps per channel,
  efficiency vs
  lifetime, and the truth-level retention maps for the 26 GeV plateau cut
  (Allie Hall's trigger-efficiency study; now part of the `base` selection).
- `event_displays_forAN.ipynb`: generator-level eta-phi and R-z displays of
  one deterministically chosen typical event from each of four grid corners,
  a decade of displacement apart (0.4 cm to 2 m), reading the ntuples
  directly (no v2/v3 dependence).


Each notebook is generated deterministically by a builder script, kept on the
development branch `truth-kinematics-forAN-dev` so that this directory holds
only what a reader needs. The committed notebooks are self-contained: given
EOS access they re-execute with
`jupyter nbconvert --to notebook --execute --inplace <name>.ipynb`, and the
figures regenerate from the stored histograms.

## The production gen filter

The v10 samples were produced **with the central gen filter applied**: at
least 4 e/mu with pT > 5 GeV, |eta| < 2.4, and production vertex within
rho < 740 cm, |z| < 960 cm. Evidence: the samples' own `GenFilterInfo`
bookkeeping (e.g. tried 4500 / passed 2325, eff = 0.517, in the retained
`4Mu MBs-500 MDp-1.2 ctau-19mm` AODSIM), and the kinematic walls in the
stored events sitting exactly on the filter thresholds. Re-applying the cut
string at the ntuple level passes 0.84-1.00 of stored events depending on the
mass point (lowest at the lightest masses, where leptons sit near the pT
threshold), a form validation on filtered events, not an efficiency. The `CutDecayFalse`
sample-name prefix refers to the MadGraph run-card flag (which removes
LHE-level cuts on the decay leptons and is indeed False); it does not mean
the samples are unfiltered. Note `edmProvDump` does not reveal this filter;
`GenFilterInfo` does.

Two consequences:

1. The lifetime truncation modeled in `lifetime_forAN.ipynb` (R_max ~ 8 m) is
   this filter's vertex cylinder, approximated by a single radius.
2. **Normalization**: the stored sum of generator weights is post-filter. Any
   absolute signal normalization must use sigma x eff_filter with the
   per-point efficiency; otherwise yields are overstated by 1.1-5.2x.

Per-point central-production efficiencies (with uncertainties) are committed
here as `central_genFilterEfficiencies.yml`
(source: `phylsix/Firefighter`, `ffConfig/python/production/Autumn18/sigmc/central/`).
Two points were re-measured from scratch with the exact production
configuration (gridpack from cvmfs, CMSSW_10_2_16_UL, GenFilterInfo counting)
and agree: 0.8607 +- 0.0035 vs 0.8593 recorded, and 0.6060 +- 0.0049 vs
0.5994, for the 2Mu2E mXX-1000 mA-0.25 prompt and lxy-300 points.

## Framework changes in this branch

- `sidm/definitions/cuts.py`: `_gen_filter_count` helper + the two
  gen-filter event cuts.
- `sidm/configs/selections.yaml`: the four new channels above.
- `sidm/configs/hist_collections.yaml`: the `gen_truth` collection.
- `sidm/definitions/hists.py`: fix. Eight leading/subleading gen dxy
  histograms masked the leptons but not the PV reference, so their fills
  silently skipped whenever any event failed the event mask; also new
  `fromA` pT/dxy histogram definitions.
- `sidm/tools/selection.py`: fix. The cut-application excepts now re-raise
  I/O errors instead of swallowing them, so a transient xrootd failure during
  lazy-array materialization fails the chunk loudly rather than silently
  corrupting the selection state. Downstream handling is the pipeline's:
  with `skipbadfiles=True` coffea skips chunks whose error text matches its
  bad-file patterns and otherwise kills the job for manual resubmission;
  it does not retry. The event-cut warning now prints the exception.
- `sidm/definitions/hists.py`: fix. The four 2D `*_nearGenA_n_genA_lxy`
  histograms mixed a per-dark-photon axis with a per-event axis, silently
  dropping every processing chunk containing an event with other than one
  matching dark photon (all 4Mu content and ~10% of 2Mu2E events). Both axes
  now count per dark photon; the canonical production carries the fixed
  definitions.
- `sidm/configs/selections.yaml`: the `genOnly_trigger` channel (anatomy
  extension).
- `sidm/definitions/hists.py`: the `anatomy_extra` histogram definitions and
  the `decayed_daughter_pairs`/`daughters_dR` helpers (per-dark-photon
  daughter pairing via the generator children links; `ak.num` needs the
  explicit `axis=2` there).
- `sidm/configs/hist_collections.yaml`: the `anatomy_extra` collection.
- `sidm/tools/histogram.py`: the fill-failure warning now prints the caught
  exception, so a broken fill function is diagnosable from job logs instead
  of failing silently per chunk.
