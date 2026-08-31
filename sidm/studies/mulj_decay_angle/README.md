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

## Contents

| file | what it is |
| --- | --- |
| `00_produce_inputs.ipynb` | **step 0**: makes the `.coffea` inputs; committed with its outputs |
| `_build_produce_inputs.py` | writes `00_produce_inputs.ipynb` |
| `mulj_decay_angle.ipynb` | the study itself; committed with its outputs |
| `_build_mulj_notebook.py` | writes `mulj_decay_angle.ipynb` |
| `_mulj_lib.py` | all loading, normalisation, statistics and drawing; one function per figure |
| `figures/` | vector PDFs of every figure, written by the study notebook; untracked |

Read them in that order: `00_produce_inputs.ipynb` produces one `.coffea` per sample, and
`mulj_decay_angle.ipynb` reads a directory of them.

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

| region | what it selects |
| --- | --- |
| `data_control_region_1muLj` | exactly one muon-type lepton jet, no cleaning cuts; the largest and least biased fake sample, and the only region that spans the full range of displacement and both muon types |
| `data_control_region_1muLj_cosmic_veto` | the same with the cosmic-muon veto applied |
| `data_control_region_1muLj_spread_cosAlpha_mu_veto` | the same with the lepton-jet impact-parameter spread cuts and the muon-pair opening-angle veto applied |
| `test_VR_2mu2e_invDisplaced_spread_cosAlpha_mu_veto` | the two-lepton-jet selection with one muon and one electromagnetic lepton jet, the displacement requirement inverted on every lepton jet, plus the spread and opening-angle cleaning |
| `test_VR_4mu_invDisplaced_spread_cosAlpha_mu_veto` | the two-muon-lepton-jet selection with the displacement requirement inverted on every lepton jet, plus the same cleaning |
| `4mu` | the four-muon signal channel: two muon-type lepton jets, full selection. **Simulation only, never data** |
| `2mu2e` | the two-muon-two-electron signal channel: one muon and one electromagnetic lepton jet, full selection. **Simulation only, never data** |

The last two rows are where the signal efficiencies of the cut scan are measured, and they
are the only place the signal shape is seen under its own selection.  Data is never read in
them; nothing in this study looks at the signal region in data.

The two validation regions are orthogonal to the signal region by construction: the signal
region requires every lepton jet to be displaced and these require the exact negation. In
both, the muon lepton jet under study is prompt by construction, which is why the
dependence of the angle on displacement is measured in the uncut single-lepton-jet control
region instead.

Simulation is also processed in the signal-region variants of those two channels and in
variants without the muon momentum cut, which are kept for cross-checks and are not drawn
in any figure here.

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
cd /uscms_data/d3/murtazas/SIDM-wt-muljangle
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
at one category, and carries the channel list and the histogram collections as arguments.
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
  --hist-collections mulj_decay_angle,mulj_decay_angle_gen
```

One merged `<sample>.coffea` comes out per sample, which is exactly what the notebook
expects. `--filelists-dir` only feeds the metadata sidecar, so run the merge once per
category if a complete sidecar is wanted for each. The chunk metadata already records that
the run was unweighted; `--unweighted-hist` is kept for outputs whose metadata predates
that key, and passing it is harmless.

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
cd /uscms_data/d3/murtazas/SIDM-wt-muljangle
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

## Caveats

These apply to every reading of the figures and are repeated in the notebook.

1. **The reconstructed angle is shaped by the selection.** The high-to-low population
   ratio printed in the legends is not what an unsculpted 1 + cos^2 decay would give, and
   the two-muon-two-electron final state loses the region near an absolute cosine of one,
   because both muons that fire the dimuon trigger have to come from the single muon
   lepton jet and the sub-leading muon of a wide-angle pair falls below threshold. The
   generator-level distribution in the same channel shows the same shape, which is how the
   sculpting is identified as acceptance rather than physics. The study compares signal
   and fake shapes under an identical selection; it is not a measurement of the angular
   distribution of the decay, which would need the efficiency corrected for.
2. **The inverted-displacement validation regions contain no displaced-standalone muons.**
   The vector builder sets the pixel-hit count of a displaced-standalone track to zero, so
   the inverted displacement requirement can only be met by a particle-flow muon with more
   than two pixel hits. Pairs containing displaced-standalone muons are therefore studied
   in the uncut single-lepton-jet control region only, and every figure that shows one of
   those regions says so on the figure.
3. **Impact parameters do not all share a reference.** The angle helper uses
   primary-vertex referenced transverse and longitudinal impact parameters for both muon
   types. The lepton-jet spread quantities, kept exactly as the cleaning cuts define them,
   mix primary-vertex referenced particle-flow values with beamspot referenced
   displaced-standalone values; for a mixed pair the transverse spread is close to the
   displaced-standalone muon's own impact parameter rather than a separation between two
   vertices. Do not call it a vertex separation, and slice by the number of
   displaced-standalone muons where the distinction matters.
4. **Yields will not match other checkouts exactly.** Other working copies of these regions
   apply a tighter displaced-standalone cross-cleaning and an eta-phi veto that is not
   applied here.
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
   evt_weights. Skipping." was seen once per validation-region channel on data before a
   duplicate muon-momentum line was removed from the selection. If it reappears, the
   simulation cutflows should not be trusted until it is understood.

## Figure conventions

Figures are designed for placement at 7.5 in of text width. `_mulj_lib.Canvas` derives the
print scale from the figure width and sizes every piece of in-figure text so that it
prints at 9 pt or more; it also refuses to let a grid grow taller than 1.15 times its
width, draws the CMS label exactly once after all panels are complete, and freezes the
layout at that point. Legends always take an explicit location and get an opaque backing
when they sit over data. Efficiency points are drawn only where the denominator holds at
least 20 entries, since a Clopper-Pearson interval on one or two entries spans almost the
whole range. Log axes are given explicit minor ticks.
