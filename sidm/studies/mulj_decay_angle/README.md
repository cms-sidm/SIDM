# Decay angle of the leading muon pair in a muon lepton jet

For every muon-type lepton jet that holds at least two muons, this study takes the two
leading-transverse-momentum constituent muons, boosts the leading one into the rest frame
of the pair, and records the angle between that boosted momentum and the direction of the
pair in the laboratory. The absolute cosine of that angle, together with the laboratory
ratio of the two muon transverse momenta, is then compared between simulated signal,
simulated background and collision data under a common set of control and validation
regions.

The question the study is built to answer is whether that angle separates the muon pairs
of a real dark-photon decay from the pairs a fake lepton jet produces: well enough to be
worth using in the selection or in the background estimate, and consistently enough across
the regions, the mediator masses and the dark-photon masses that the separation measured in
one of them can be carried to another.

The inputs were produced with the channel definitions kept in
`sidm/configs/selections_mulj_decay_angle_as_run.yaml`, which are not main's; "Selections used
in this study" below gives the differences and how the file is tied to the inputs.

## Contents

| file | what it is |
| --- | --- |
| `00_produce_inputs.ipynb` | **step 0**: makes the `.coffea` inputs; committed with its outputs |
| `_build_produce_inputs.py` | writes `00_produce_inputs.ipynb` |
| `mulj_decay_angle.ipynb` | the study itself; committed with its outputs |
| `_build_mulj_notebook.py` | writes `mulj_decay_angle.ipynb` |
| `_mulj_lib.py` | all loading, normalisation, statistics and drawing; one function per figure |
| `check_gen_boost_precision.py` | an offline check of one unselected 4mu file per bound-state mass at m_Zd = 0.25 GeV, read over xrootd: the generated \|cos(theta*)\| with the processor's single-precision boost, in double precision and without a boost; its numbers are quoted in section b and caveat 9 |
| `../../configs/selections_mulj_decay_angle_as_run.yaml` | the channel definitions the three campaigns ran, passed to the processor by the submit files and by `00_produce_inputs.ipynb` |
| `figures/` | vector PDFs of every figure, written by the study notebook; untracked |

Read them in that order: `00_produce_inputs.ipynb` produces one `.coffea` per sample, and
`mulj_decay_angle.ipynb` reads a directory of them.

The study notebook runs a to j, and three campaigns are behind it.

| sections | what they measure | inputs |
| --- | --- | --- |
| a - h | one muon lepton jet: the decay angle of its two hardest muons, what it correlates with, what a cut on it would buy, the caveats and the verdict | `INPUT_DIR`, the first campaign |
| i | the two muon lepton jets of a four-muon event together: whether their angles are correlated, and whether the angle still separates inside a mass window | `lib.configure_pair(...)`, a second campaign |
| j | the lepton-jet pair as a decay-plane system: the signed decay angle, the electron pair of an electromagnetic lepton jet, the azimuthal angle between the two decay planes, the transverse balance and the missing momentum | `lib.configure_planes(...)`, `LjAnglePlanes_Run2018_v1` |

Sections i and j each set their own input directory in their own configuration cell, so
changing `INPUT_DIR` does not move them.

Both notebooks are generated. Edit the builder or the library, regenerate, re-execute; hand
edits to a `.ipynb` put it out of step with its source.

## Quantities

All of them are per muon lepton jet, from its two hardest constituent muons. The
constituent list mixes particle-flow and displaced-standalone muons and is sorted by
transverse momentum inside each lepton jet; lepton jets with fewer than two muons are
excluded explicitly rather than through a channel's muon-multiplicity cut.

* Each muon is given a four-vector carrying the muon mass, 0.105658 GeV. The pair is the
  sum of the two, and its invariant mass is the pair mass used throughout.
* The decay angle is the angle between the leading muon after it is boosted into the pair
  rest frame and the laboratory direction of the pair. Its absolute cosine is stored,
  because exchanging the two muons flips the sign.
* The momentum ratio is the sub-leading over the leading transverse momentum in the
  laboratory.
* Recorded alongside, for each pair: lepton-jet transverse momentum and isolation, pair
  mass, angular separation of the two muons, the smaller and the larger of their
  transverse impact parameters, how many of the two are displaced-standalone muons, each
  muon's transverse momentum, and the lepton-jet impact-parameter spreads that the
  cleaning cuts use.

Every axis of a given histogram is filled from one shared jagged array, so each pair
contributes exactly one entry per axis. The computation is a single helper in
`sidm/definitions/hists.py`; the per-axis functions only pick fields off it.

A muon lepton jet may hold more than two muons. This study always takes the two hardest
and does not attempt to choose the pair that best reconstructs a dark photon, so in a
lepton jet with three or more muons the pair need not be the true decay pair. The number
of displaced-standalone muons in the pair is recorded so that population can be separated.

## Regions

| channel | named in the notebook | what it selects |
| --- | --- | --- |
| `data_control_region_1muLj` | uncut single-lepton-jet control region | exactly one lepton jet, a muon lepton jet, no cleaning cuts; the largest and least biased fake sample, the only region that spans the full range of displacement, and the one with the most pairs containing displaced-standalone muons in number |
| `data_control_region_1muLj_cosmic_veto` | cosmic-veto single-lepton-jet control region | the same with the cosmic-muon veto applied |
| `data_control_region_1muLj_spread_cosAlpha_mu_veto` | cleaned single-lepton-jet control region | the same with the lepton-jet impact-parameter spread cuts and the muon-pair opening-angle veto applied |
| `test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto` | two-muon-two-electron validation region | the two-lepton-jet selection with one muon and one electromagnetic lepton jet, the displacement requirement inverted on every lepton jet, plus the spread and opening-angle cleaning |
| `test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto` | four-muon validation region | the two-muon-lepton-jet selection with the displacement requirement inverted on every muon lepton jet, plus the same cleaning |
| `4mu` | four-muon channel | the inclusive four-muon channel: two muon-type lepton jets and two muons above 26 GeV, with no isolation, displacement or separation requirement (`4mu_sr` adds those). **Simulation only, never data** |
| `2mu2e` | two-muon-two-electron channel | the inclusive two-muon-two-electron channel: one muon and one electromagnetic lepton jet, likewise. **Simulation only, never data** |

The last two rows are where the signal efficiencies of the cut scan are measured, and they
are the only place the signal shape is seen under its own selection.  Data is never read in
them; nothing in this study looks at the signal region in data.

The two validation regions are orthogonal to the signal region by construction: the signal
region requires every lepton jet to be displaced and these require the exact negation. For a
muon lepton jet the negation is that at least one of its muons has more than two pixel hits.
A displaced-standalone muon is given zero pixel hits, so that muon is always a particle-flow
muon, but the requirement does not exclude displaced-standalone muons from the lepton jet,
nor from the two leading muons the angle is computed from: 29 percent of the data pairs in
the four-muon validation region and 10 percent in the two-muon-two-electron one contain one
(caveat 2). The dependence of the angle on displacement is measured in the uncut
single-lepton-jet control region, which has no requirement on displacement.

Simulation is also processed in the signal-region variants of those two channels and in
variants without the muon momentum cut, which are kept for cross-checks and are not drawn
in any figure here.

## Selections used in this study

The inputs come from three Condor campaigns on 2018 data and simulation, run from the branch
`mulj-decay-angle-dev`:

| campaign | ran | merged outputs | read by |
| --- | --- | --- | --- |
| `MuljDecayAngle_Run2018_v1` | 2026-08-28 to 2026-08-31 | `mulj_merged` | sections a to h |
| `MuljPairAngle_Run2018_v1` | 2026-08-31 to 2026-09-01 | `mulj_pair_merged` | section i |
| `LjAnglePlanes_Run2018_v1` | 2026-09-02 to 2026-09-05 | `lj_angle_planes_merged` | section j |

All three ran the channel definitions of `sidm/configs/selections.yaml` at commit `61224901`
of that branch. Those definitions are kept in
`sidm/configs/selections_mulj_decay_angle_as_run.yaml`, with every YAML anchor resolved, every
cut list in the order it is applied, and every channel under the name stored in the
histograms. Load it with

```python
SidmProcessor(channels, collections,
              selections_cfg="configs/selections_mulj_decay_angle_as_run.yaml")
```

or pass `--selections-cfg configs/selections_mulj_decay_angle_as_run.yaml` to
`condor/run_sidm_chunk.py`, which `condor/run_job.sh` forwards as its ninth argument. The
study's submit files and `00_produce_inputs.ipynb` already do. Without it the processor reads
`configs/selections.yaml`, whose channels of the same names select different events.

Three pieces of evidence tie the file to the inputs. The jobs ran from a working copy with
uncommitted changes, and commit `61224901` was made on 2026-09-07, after all three campaigns
had run. The metadata file written next to each merged output (152, 101 and 154 such files)
records a definition for every channel, and every recorded definition agrees with the as-run
file; those definitions are read from the checkout that did the merging, not from the jobs,
and the commit recorded there is likewise the merging checkout's. For the third campaign the
code shipped to the jobs was kept and is identical, file by file, to commit `61224901`; the
code shipped to the first two campaigns was not kept. The strongest evidence is a rerun:
stored Condor outputs of all three campaigns (signal chunks of each, data chunks of the first
and third) rerun from this branch with the as-run file reproduce the stored histograms and
cutflows bin for bin.

The merged outputs are on EOS under `/store/group/lpcmetx/SIDM/coffea_outputs/murtazas/`, in
`mulj_decay_angle`, `mulj_pair_angle` and `lj_angle_planes` for the three campaigns in that
order, and the unmerged chunks under `/store/user/murtazas/sidm_condor/<campaign>`. The local
directories in the table are copies of the first.

**How these differ from the analysis-wide definitions.** `sidm/configs/selections.yaml` on
main (commit `5bb94cd`) changed the base cuts every channel inherits (PR #56) after these
campaigns ran:

* main applies an eta-phi veto to electrons, photons, muons and displaced-standalone muons;
  none of the channels here has it;
* main has `>=2 muons (pf or dsa) pT > 26 GeV` and `all cosAlpha(mu, mu) > -0.99` among its
  base event cuts. Here the 26 GeV requirement is in `4mu`, `4mu_sr`, `2mu2e`, `2mu2e_sr` and
  the four spread-veto regions, and absent from the two `_noMuPtCut` variants, the uncut and
  cosmic-veto single-lepton-jet control regions and the electromagnetic-lepton-jet
  control region; the opening-angle veto is in the four spread-veto regions only;
* main cuts on `dzSpread` and `dxySpread` in the muon-lepton-jet cuts of its signal channels.
  None of the signal-like channels here has a spread cut; the four spread-veto regions cut on
  `dzSpread` and `vxySpread`, with main's thresholds (50 for all muons, 10 for particle-flow
  muons).

The four spread-veto regions, `data_control_region_1muLj_spread_cosAlpha_mu_veto`,
`cr_2lj_invDisplaced_spread_cosAlpha_mu_veto` and the two `test_VR_*` regions, are not defined
on main at all. The notebook renders, when it runs, the table below and a channel-by-channel
comparison with the checkout's `configs/selections.yaml`. The table lists what each channel
adds to the cuts every channel shares; generator-level particle selections are left out.

| channel | object cuts | lepton-jet cuts | event cuts |
| --- | --- | --- | --- |
| **baseline, every channel** | electrons: pT > 10 GeV, \|eta\| < 2.4, MVANonIsoWPL<br>muons: looseID, pT > 5 GeV, \|eta\| < 2.4<br>photons: pT > 20 GeV, \|eta\| < 2.5, Custom Cutbased, pixelSeed, Photon DR Veto 0p025<br>dsaMuons: pT > 10 GeV, \|eta\| < 2.4, displaced ID, all | ljs: pT > 30 GeV, \|eta\| < 2.4 | pass triggers, PV filter |
| ***signal-like (simulation only)*** | | | |
| `4mu` | none beyond the baseline | mu_ljs: Mu >= 2 | >=2 LJs, >=2 muons (pf or dsa) pT > 26 GeV, 4mu |
| `4mu_sr` | none beyond the baseline | mu_ljs: Mu >= 2, mu_lj_iso < 0.1, displaced<br>egm_ljs: egm_lj_iso < 0.2, displaced | >=2 LJs, >=2 muons (pf or dsa) pT > 26 GeV, LJ-LJ dPhi > 2, 4mu |
| `4mu_noMuPtCut` | none beyond the baseline | mu_ljs: Mu >= 2 | >=2 LJs, 4mu |
| `2mu2e` | none beyond the baseline | mu_ljs: Mu >= 2 | >=2 LJs, >=2 muons (pf or dsa) pT > 26 GeV, 2mu2e |
| `2mu2e_sr` | none beyond the baseline | mu_ljs: Mu >= 2, mu_lj_iso < 0.1, displaced<br>egm_ljs: egm_lj_iso < 0.2, displaced | >=2 LJs, >=2 muons (pf or dsa) pT > 26 GeV, LJ-LJ dPhi > 2, 2mu2e |
| `2mu2e_noMuPtCut` | none beyond the baseline | mu_ljs: Mu >= 2 | >=2 LJs, 2mu2e |
| ***control*** | | | |
| `data_control_region_1muLj` | none beyond the baseline | mu_ljs: Mu >= 2 | pass flags, = 1 LJs, = 1 muLJs |
| `data_control_region_1muLj_cosmic_veto` | dsaMuons: \|dxy\| <= 40, \|dz\| <= 60 | mu_ljs: Mu >= 2 | pass flags, = 1 LJs, = 1 muLJs, all cos_alpha(dsa, dsa) > -0.9 |
| `data_control_region_1muLj_spread_cosAlpha_mu_veto` | none beyond the baseline | mu_ljs: Mu >= 2, dzSpread_mu <= 50, vxySpread_mu <= 50, dzSpread_pf <= 10, vxySpread_pf <= 10 | pass flags, >=2 muons (pf or dsa) pT > 26 GeV, = 1 LJs, = 1 muLJs, all cosAlpha(mu, mu) > -0.99 |
| `data_control_region_1egmLj_cosmic_veto` | dsaMuons: \|dxy\| <= 40, \|dz\| <= 60 | none beyond the baseline | pass flags, = 1 LJs, = 1 egmLJs, all cos_alpha(dsa, dsa) > -0.9 |
| `cr_2lj_invDisplaced_spread_cosAlpha_mu_veto` | none beyond the baseline | mu_ljs: Mu >= 2, inverse displaced, dzSpread_mu <= 50, vxySpread_mu <= 50, dzSpread_pf <= 10, vxySpread_pf <= 10<br>egm_ljs: inverse displaced | >=2 LJs, >=2 muons (pf or dsa) pT > 26 GeV, all cosAlpha(mu, mu) > -0.99 |
| ***validation*** | | | |
| `test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto` | none beyond the baseline | mu_ljs: Mu >= 2, inverse displaced, dzSpread_mu <= 50, vxySpread_mu <= 50, dzSpread_pf <= 10, vxySpread_pf <= 10 | >=2 LJs, >=2 muons (pf or dsa) pT > 26 GeV, >=2 mu_ljs, muLJ-muLJ dPhi >= 2, muLJ-muLJ invariant mass >= 150 GeV, all cosAlpha(mu, mu) > -0.99 |
| `test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto` | none beyond the baseline | mu_ljs: Mu >= 2, inverse displaced, dzSpread_mu <= 50, vxySpread_mu <= 50, dzSpread_pf <= 10, vxySpread_pf <= 10<br>egm_ljs: inverse displaced | >=2 LJs, >=2 muons (pf or dsa) pT > 26 GeV, >=1 mu_ljs, >=1 egm_ljs, muLJ-egmLJ dPhi >= 2, muLJ-egmLJ invariant mass >= 150 GeV, all cosAlpha(mu, mu) > -0.99 |

## Reproducing the inputs

`00_produce_inputs.ipynb` is the executable version of this section. It runs two files of
one signal sample and two of one data sample end to end, with the same processor, channels
and histogram collections the campaign used, and its committed outputs are the evidence
that the recipe works. The full campaign is written out there as commands, the same ones
repeated below.

Everything runs on the LPC. Each `ssh` is a fresh shell, so chain the environment setup
into one invocation.

### 1. Build the job arguments, then submit

The file lists and job-argument files are generated, not stored. Build them with the
committed generator, one call per category; each writes one file list per chunk under
`--outdir` and one line per chunk to `--job-args`.

```bash
cd <repository root>
source /uscms_data/d3/murtazas/SIDM/sidm_venv/bin/activate

python condor/make_job_args.py --samples-file condor/mulj_signal_4mu_samples.txt \
  --location-cfg signal_4mu_v10.yaml --files-per-job 5 --replace-xcache \
  --outdir condor/filelists_mulj_signal --job-args condor/job_args_mulj_signal_4mu.txt

python condor/make_job_args.py --samples-file condor/mulj_signal_2mu2e_samples.txt \
  --location-cfg signal_2mu2e_v10.yaml --files-per-job 5 --replace-xcache \
  --outdir condor/filelists_mulj_signal --job-args condor/job_args_mulj_signal_2mu2e.txt

python condor/make_job_args.py --samples-file condor/mulj_background_samples.txt \
  --location-cfg backgrounds.yaml --dataset skimmed_llpNanoAOD_v2 \
  --files-per-job 30 --replace-xcache \
  --outdir condor/filelists_mulj_bkg --job-args condor/job_args_mulj_bkg.txt

python condor/make_job_args.py --samples-file condor/mulj_data_samples.txt \
  --location-cfg data_skimmed.yaml --files-per-job 30 --replace-xcache \
  --outdir condor/filelists_mulj_data --job-args condor/job_args_mulj_data.txt
```

The chunking above is what the campaign ran: 5 files per job for signal, 30 for background
and data, giving 701, 2,109 and 4,278 jobs over 152 samples and 193,770 input files.

Submit one cluster per category from inside `condor/`, with a proxy on shared NFS rather
than in `/tmp`, since the schedd cannot read `/tmp`:

```bash
cd condor
export X509_USER_PROXY=/uscms_data/d3/murtazas/x509_proxy.pem
condor_submit submit_mulj_signal.sub    #   701 jobs
condor_submit submit_mulj_bkg.sub       # 2,109 jobs
condor_submit submit_mulj_data.sub      # 4,278 jobs
```

Each submit file is `condor/submit.sub` with its `job_args` and `filelists` paths pointed
at one category, and carries the channel list, the histogram collections and the as-run
selections file as arguments.
The sample lists, the three submit files, the job-argument generator and the merge script
are all in the repository; only the generated file lists and job-argument files are not,
which is why the step above regenerates them. A drained queue is not success: read the
exit-code histogram and the chunk count on EOS. When a cluster drains, reconcile it with
the committed tool:

```bash
python condor/condor_campaign.py reconcile \
  --job-args condor/job_args_mulj_data.txt \
  --logs-dir condor/logs_mulj_data \
  --eos-chunk-dir /store/user/murtazas/sidm_condor/MuljDecayAngle_Run2018_v1 \
  --run-id mulj_data_v1
```

(swap `data` for `bkg` and `signal`).

### 2. Merge

All three categories write into the same EOS directory and sample names are unique across
them, so one merge pass covers everything. The chunks were produced with unweighted
histograms, so the merge has to be told:

```bash
python sidm/scripts/merge_coffea_chunks_eos.py \
  --input-eos-dir  /store/user/murtazas/sidm_condor/MuljDecayAngle_Run2018_v1 \
  --output-eos-dir /store/group/lpcmetx/SIDM/coffea_outputs/murtazas/mulj_decay_angle \
  --unweighted-hist \
  --filelists-dir condor/filelists_mulj_data \
  --hist-collections mulj_decay_angle,mulj_decay_angle_gen \
  --selections-cfg configs/selections_mulj_decay_angle_as_run.yaml \
  --selections data_control_region_1muLj,data_control_region_1muLj_cosmic_veto,\
data_control_region_1muLj_spread_cosAlpha_mu_veto,test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto,\
test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto,4mu,2mu2e,4mu_sr,2mu2e_sr,4mu_noMuPtCut,2mu2e_noMuPtCut
```

One merged `<sample>.coffea` comes out per sample, which is exactly what the notebook
expects. `--filelists-dir` only feeds the metadata sidecar, so run the merge once per
category if a complete sidecar is wanted for each. The chunk metadata already records that
the run was unweighted; `--unweighted-hist` is kept for outputs whose metadata predates
that key, and passing it is harmless. The channel names and `--selections-cfg` only feed
the sidecar, which then records the definitions the chunks were run with rather than those
of `configs/selections.yaml`.

### 3. Run the notebook

Point the notebook at the merged outputs by editing one cell. If they are read from EOS,
either give the `root://` URL directly, in which case files are cached locally on first
read, or copy them down first with `xrdcp`; `coffea.util.load` cannot open a URL itself.

```python
INPUT_DIR = "root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/coffea_outputs/murtazas/mulj_decay_angle"
LUMI_FB = 59.8
```

Then rebuild and execute:

```bash
cd <repository root>
source /uscms_data/d3/murtazas/SIDM/sidm_venv/bin/activate
export PYTHONPATH=$PWD
cd sidm/studies/mulj_decay_angle
python _build_mulj_notebook.py
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=sidm_venv \
    --ExecutePreprocessor.timeout=1800 mulj_decay_angle.ipynb
```

The sample lists are read from whatever is in `INPUT_DIR`, so no list has to be edited.
Signal samples are recognised from names of the form `4Mu_500GeV_1p2GeV_1p9mm`, data from
`DoubleMuon_2018C_0`, and a simulated background is anything else that has a cross section
configured in `sidm/configs/cross_sections.yaml`. Panels with no input degrade to a
labelled empty panel and a note in the run record at the end of the notebook, rather than
raising.

Figures are written to `figures/` as vector PDFs. That directory is untracked, and the
repository ignores PDFs; the rendered figures a reviewer sees are the ones committed
inside the notebook.

### 4. The lepton-jet pair inputs (section i only)

Section i reads a second, separate set of merged outputs, produced by a later campaign that
filled two event-level histograms: `mulj_pair_absCosTheta0_vs_absCosTheta1`, one entry per
event holding two qualifying muon lepton jets, and `mulj_pair_absCosTheta_vs_mjj`, one entry
per qualifying lepton jet against the muon-lepton-jet-pair mass of its event. Its channels
are `4mu`, `4mu_sr`, `4mu_noMuPtCut` and the four-muon validation region, and collision
data was written in that validation region only. The directory is set
in the section's own configuration cell, independently of `INPUT_DIR`:

```python
lib.configure_pair("/uscms_data/d3/murtazas/mulj_pair_merged")
# or root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/coffea_outputs/murtazas/mulj_pair_angle
```

### 5. The decay-plane inputs (section j only)

Section j reads a third set of merged outputs, from the campaign `LjAnglePlanes_Run2018_v1`,
submitted on 2026-09-02 as 7,096 jobs over 154 samples and merged on 2026-09-05. It adds three
groups of histograms to the ones above:

* the **signed** decay angle of a muon lepton jet, `mu_lj_mumu_cosTheta` and
  `mu_lj_mumu_cosTheta_vs_mass`, where the sign is set by the negatively charged muon;
* the same quantities for the two leading **electrons** of an electromagnetic lepton jet,
  `egm_lj_ee_absCosTheta`, `egm_lj_ee_cosTheta`, `egm_lj_ee_ptRatio`, `egm_lj_ee_mass`,
  `egm_lj_ee_dR` and the four planes among them;
* the **event-level decay-plane observables** of the two leading qualifying lepton jets of
  either type, `lj_pair_cosThetaA0_vs_cosThetaA1`, `lj_pair_absCosThetaA_vs_mjj`,
  `lj_pair_Phi`, `lj_pair_planeAngle`, `lj_pair_planeAngle_vs_mjj`, `lj_pair_mjj`,
  `lj_pair_type`, `lj_pair_hasCharge`, `lj_pair_ptAsym`, `lj_pair_ptVecSum`, `lj_pair_met`,
  `lj_pair_metCorr` and `lj_pair_metCorr_vs_ptVecSum`, together with their generator-level
  counterparts `gen_pair_Phi`, `gen_pair_planeAngle` and `gen_pair_cosThetaA0_vs_cosThetaA1`
  on the signal samples.

Its collections are `lj_decay_angle` and `lj_pair_decay_planes` everywhere, plus
`lj_pair_decay_planes_gen` on the two signal categories. Its channels are the four signal
channels `4mu`, `4mu_sr`, `2mu2e`, `2mu2e_sr` and five control and validation regions:

| region | what it selects |
| --- | --- |
| `data_control_region_1muLj_spread_cosAlpha_mu_veto` | exactly one lepton jet, of muon type, with the spread cleaning and the opening-angle veto |
| `data_control_region_1egmLj_cosmic_veto` | exactly one lepton jet, of electromagnetic type, with the cosmic veto |
| `cr_2lj_invDisplaced_spread_cosAlpha_mu_veto` | two or more lepton jets of either type, the displacement requirement inverted on all of them, with the spread cleaning and the opening-angle veto; no requirement on the mass of the lepton-jet pair or on the azimuthal separation of the two jets, so it is the widest two-lepton-jet fake sample |
| `test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto` | the same for two muon lepton jets, with a muon-lepton-jet-pair mass of at least 150 GeV and an azimuthal separation of at least 2, which the four-muon channel does not require |
| `test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto` | the same for one muon and one electromagnetic lepton jet, with the same pair-mass and separation requirements |

**Collision data was written in those five regions only**, and `lib.planes_counts` raises rather
than returning anything if a data sample is passed with any other channel, so a blinding mistake
cannot reach a figure.

The signal grid is the 4Mu and 2Mu2E points at m_XX = 100, 500 and 1000 GeV, every lifetime;
the search uses bound-state masses of 200 GeV and above, so section j leads every number with
500 and 1000 GeV and says where 200 to 250 GeV is not covered.

Point the section at its inputs in its own configuration cell:

```python
lib.configure_planes("/uscms_data/d3/murtazas/lj_angle_planes_merged")
# or root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/coffea_outputs/murtazas/lj_angle_planes
```

### 6. Rerunning the second and third campaigns

They are built and submitted as in step 1, one `make_job_args.py` call per row with the same
location configs (`signal_4mu_v10.yaml`, `signal_2mu2e_v10.yaml`, `backgrounds.yaml` with
`--dataset skimmed_llpNanoAOD_v2`, and `data_skimmed.yaml`), `--outdir condor/filelists_<name>`
and `--job-args condor/job_args_<name>.txt`, where `<name>` is the submit file's name without
`submit_`. Every submit file passes the as-run selections file.

| submit file | sample list | samples | files per job | jobs | channels | collections |
| --- | --- | --- | --- | --- | --- | --- |
| `submit_muljpair_signal.sub` | `mulj_signal_4mu_samples.txt` | 33 | 10 | 199 | `4mu`, `4mu_sr`, `4mu_noMuPtCut`, four-muon validation region | `mulj_decay_angle_pair` |
| `submit_muljpair_data.sub` | `mulj_data_samples.txt` | 68 | 60 | 2,156 | four-muon validation region | `mulj_decay_angle_pair` |
| `submit_ljangle_signal4mu.sub` | `mulj_signal_4mu_samples.txt` | 33 | 10 | 199 | `4mu`, `4mu_sr`, `2mu2e`, `2mu2e_sr` and the five regions of step 5 | `lj_decay_angle`, `lj_pair_decay_planes`, `lj_pair_decay_planes_gen` |
| `submit_ljangle_signal2mu2e.sub` | `lj_signal_2mu2e_samples.txt` | 33 | 10 | 163 | as above | as above |
| `submit_ljangle_bkg.sub` | `lj_bkg_samples.txt` | 20 | 20 | 3,520 | as above | `lj_decay_angle`, `lj_pair_decay_planes` |
| `submit_ljangle_data.sub` | `mulj_data_samples.txt` | 68 | 40 | 3,214 | the five regions of step 5 only | `lj_decay_angle`, `lj_pair_decay_planes` |

The first two rows write to `/store/user/murtazas/sidm_condor/MuljPairAngle_Run2018_v1`, the
other four to `LjAnglePlanes_Run2018_v1/{signal4mu,signal2mu2e,bkg,data}` next to it. Merge as
in step 2, with that campaign's chunk directories, collections and channels, into the merged
directories named under "Selections used in this study".

## Caveats

These apply to every reading of the figures and are repeated in the notebook.

1. **The reconstructed angle is shaped by the selection.** The high-to-low population
   ratio printed in the legends is not what an unsculpted 1 + cos^2 decay would give, and
   the two-muon-two-electron final state loses the region near an absolute cosine of one,
   because both muons that fire the dimuon trigger have to come from the single muon
   lepton jet and the sub-leading muon of a wide-angle pair falls below threshold. The
   generator-level distribution in the same channel shows the same shape at m_Zd = 1.2 and
   5.0 GeV, which is how the sculpting is identified as acceptance rather than physics; at
   m_Zd = 0.25 GeV the generated histogram is not a sound reference (caveat 9). The study
   compares signal and fake shapes under an identical selection; it is not a measurement of
   the angular distribution of the decay, which would need the efficiency corrected for.
2. **The inverted-displacement requirement does not exclude displaced-standalone muons.**
   In the two validation regions and the two-lepton-jet control region every muon lepton
   jet must contain at least one muon with more than two pixel hits. The vector builder
   sets the pixel-hit count of a displaced-standalone track to zero, so that muon is a
   particle-flow muon, but the two leading muons the angle is computed from can include
   displaced-standalone ones. In data the pairs split 71, 29 and 0.2 percent into
   particle-flow, mixed and displaced-standalone pairs in the four-muon validation region,
   90, 10 and 0.2 percent in the two-muon-two-electron one, 59, 41 and below 0.1 percent in
   the two-lepton-jet control region, and 93, 6.5 and 0.8 percent in the uncut
   single-lepton-jet control region; the notebook prints the table in section a, with the
   fraction of each class above |cos(theta*)| = 0.4. The region dependence found in
   sections e and h is not created by composition. Within particle-flow pairs alone the
   fraction above 0.4 runs from 0.12 in the cleaned single-lepton-jet control region to
   0.59 and 0.62 in the four-muon validation region and the two-lepton-jet control region,
   while mixed pairs stay between 0.13 and 0.27 in every region. Any single number quoted
   for an inverted-displacement region is therefore an average over two populations with
   unlike shapes. Section e separates pairs by the number of displaced-standalone muons in
   the single-lepton-jet control regions; the table in section a gives the fraction above
   0.4 for each class in every region. Every figure that shows a validation region states
   its requirement in the region label.
3. **Impact parameters do not all share a reference.** The angle helper uses
   primary-vertex referenced transverse and longitudinal impact parameters for both muon
   types. The lepton-jet spread quantities, kept exactly as the cleaning cuts define them,
   mix primary-vertex referenced particle-flow values with beamspot referenced
   displaced-standalone values; for a mixed pair the transverse spread is close to the
   displaced-standalone muon's own impact parameter rather than a separation between two
   vertices. Do not call it a vertex separation, and slice by the number of
   displaced-standalone muons where the distinction matters.
4. **These are the as-run definitions, not main's.** Every channel here is defined by
   `selections_mulj_decay_angle_as_run.yaml`. Main's channels of the same names add an
   eta-phi veto and the base-cut and spread-cut differences listed under "Selections used in
   this study", and under them the shapes would change as well as the yields.
5. **Histograms are raw counts; cutflows are not.** The campaign ran with unweighted
   histograms, so nothing in them carries a cross-section or generator weight, while the
   cutflows in the same files are luminosity times cross-section scaled. Never combine the
   two. A normalised simulation prediction needs luminosity times cross section divided by
   the number of *generated* events; the sum of weights in the metadata counts events that
   survived the skim, and the signal skim keeps roughly a tenth of what was generated, so
   the generated count from the file census has to be substituted before a yield is
   quoted. `_mulj_lib.lumi_xs_weight()` implements the lookup and says so in its
   docstring. Shapes are unaffected, because the factor cancels in a unit-area
   normalisation, and every comparison in this study is unit-area.
6. **Signal-region channels evaluate their event cuts on unfiltered lepton jets**, so an
   event whose lepton jets are all prompt can pass a signal-region channel with no muon
   lepton jets left in it. Orthogonality between the validation and signal regions holds
   at the level of histogram entries, not events.
7. **Two axes merge things that are not the same.** The lepton-jet count in the collection
   counts unfiltered lepton jets, that is, the event multiplicity. Isolation fills a
   missing value with zero, so "no matched jet" and "isolation exactly zero" land in the
   same bin and in the first isolation slice.
8. **One warning to watch for in the campaign logs.** "Unable to apply event cuts to
   evt_weights. Skipping." means an event cut was not applied to the event weights. Where it
   appears, the cutflows of that channel should not be trusted until it is understood.
9. **Numerical precision.** The generated histograms build the generator-level momenta from
   the GenPart px, py and pz branches, which are full single precision (the pt, eta, phi and
   mass branches keep 8 mantissa bits, and of those only the dark photon's mass is read), and
   every boost is computed in single precision. At m_Zd = 0.25 GeV and m_XX = 500 or
   1000 GeV the dark photon's Lorentz factor is 1,000 to 2,000 and the boost into its rest
   frame pushes the generated |cos(theta*)| towards one, so that histogram is not a reference
   there; the reconstructed angle at that mass carries the same arithmetic. The plane angle,
   and every generated reference at m_Zd = 1.2 and 5.0 GeV, are unaffected. The notebook's
   section b gives a toy that sizes this point by point. The overflow bin of each angle
   histogram holds the entries this arithmetic places outside the axis: a few in a hundred
   thousand in data, at most 0.2 percent of the simulated entries at m_Zd = 1.2 and 5.0 GeV,
   and 2 to 18 percent at m_Zd = 0.25 GeV with m_XX = 500 and 1000 GeV, where the generated
   ones are cosines the single-precision boost returns as not a number (an offline check of one
   unselected file per point, `check_gen_boost_precision.py`). Section b prints these overflow
   fractions from the stored histograms. Every figure is built from the in-range bins.

## Figure conventions

Figures are designed for placement at 7.5 in of text width. `_mulj_lib.Canvas` derives the
print scale from the figure width and sizes every piece of in-figure text so that it
prints at 9 pt or more; it also refuses to let a grid grow taller than 1.15 times its
width, draws the CMS label exactly once after all panels are complete, and freezes the
layout at that point. Legends always take an explicit location and get an opaque backing
when they sit over data. Efficiency points are drawn only where the denominator holds at
least 20 entries, since a Clopper-Pearson interval on one or two entries spans almost the
whole range. Log axes are given explicit minor ticks.
