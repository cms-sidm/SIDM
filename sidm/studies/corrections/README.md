# Standard corrections with correctionlib

This study is the first step towards real systematics: apply the standard CMS object
corrections through [correctionlib](https://cms-nanoaod.github.io/correctionlib/),
make sure they are applied correctly, and find the places where the standard
corrections do not describe our objects. It also documents every change to the
framework that came with it.

**Only simulation is used here.** The analysis is blinded, so no collision data is
read anywhere in this study, and the Condor jobs refuse data unless explicitly told
otherwise. The data/MC checks that normally go with these corrections are listed under
[open items](#open-items), for the group to decide where and when to look.

## In short

* `SidmProcessor(..., corrections=True)` applies the standard corrections to
  simulation. With the default, `corrections=None`, the output is identical to
  `main`, cutflows and histograms bin by bin.
* The default set is **pileup reweighting, the L1 prefiring weight, the PF-muon
  Loose ID SF and the electron reconstruction SF**, read from dated CVMFS snapshots
  pinned in [`sidm/configs/corrections.yaml`](../../configs/corrections.yaml).
* Every cutflow row carries the weighted count of every **weight category**
  (uncorrected, one correction at a time, every up/down variation), `syst_hists=True`
  adds a `syst` axis with the same categories to every histogram, and the output
  metadata records what was applied.
* A **Condor campaign** ran the `corrections_study` histogram collection over every
  v10 signal point (180) and every skimmed background (20 samples), 79k files and
  about 280M events. The merged outputs are on EOS for anyone in the group, and the
  notebooks read them from there.
* In the SR the default corrections change the signal yields by $-2.2$ to $+0.7$%
  (2Mu2E, median $-0.7$%) and $-1.8$ to $+0.2$% (4Mu, median $-1.4$%), mostly through
  the prefiring weight, and bring uncertainties of 0.15–4.3% and 0.2–3.0% (median 0.6%
  for both). The per-point lnN values are in
  [`tables/sr_correction_nuisances.csv`](tables/sr_correction_nuisances.csv).
* The standard corrections miss most of what makes our objects special: 79–100% of
  the lepton-jet PF muons are outside the phase space the muon SFs were measured in
  (7% of DY muons), and DSA muons, the displaced-muon triggers and our electron and
  photon IDs have no SF at all. That is where the real systematics will come from.
* Found along the way, and worth a look by everyone: the Condor merge script inflated
  every merged weighted yield by the number of chunks (fixed); the Condor jobs
  normalised skimmed backgrounds without their skim factor (fixed); and one event of
  the `DYJetsToMuMu_M50` skim has a generator weight of 1.45e14, which sets that
  sample's normalisation ([found along the way](#found-along-the-way)).

## Using it

```python
from sidm.tools import sidm_processor

p = sidm_processor.SidmProcessor(["2mu2e_sr", "4mu_sr"], ["lj_lj_base"],
                                 corrections=True,     # the default list
                                 syst_hists=True,      # optional: syst axis on every hist
                                 monitor_corrections=["photon_pixel_veto"])  # optional
out = runner.run(fileset, treename="Events", processor_instance=p)["out"]

cf = out[sample]["cutflow"]["2mu2e_sr"]
cf.variation_yields()          # {"nominal": ..., "uncorrected": ..., "pileupUp": ..., ...}
h = out[sample]["hists"]["lj_lj_invmass"][{"channel": "2mu2e_sr", "syst": "muon_idUp"}]
out[sample]["metadata"]["object_corrections"]   # which object SFs each channel used
```

The weight categories, in the cutflows and on the `syst` axis:

| category | weight |
|---|---|
| `nominal` | genWeight × every applied correction |
| `uncorrected` | genWeight only |
| `only_<name>` | genWeight × one correction: what it does alone |
| `<name>Up`, `<name>Down` | nominal with one correction shifted |
| `with_<name>`, `with_<name>Up`, `with_<name>Down` | nominal × a monitored correction (evaluated, not applied) and its shifts |

* `corrections=["pileup", "l1_prefiring"]` applies a subset; any entry of the config
  can be named, including ones that are not in the default list
  (`photon_pixel_veto`, `electron_id`, `photon_id`).
* Data is never corrected. With `syst_hists=True` data fills only the `nominal`
  category.
* The luminosity × cross section normalisation still uses the plain genWeight sum,
  so a correction that changes the total (pileup reweighting does, slightly) shows up
  in the yields. The metadata key `scaled_sum_weights_corrected` holds the sums with
  each event-level correction, before any selection.
* `syst_hists=True` stores one copy of every histogram per category, 14 with the
  default list and 17 with the pixel veto monitored, so use it with small
  collections.
* **Histograms of the corrections themselves**: with corrections on, fill functions
  can read `objs["weight_factors"]` (one column per factor and shift, e.g. `pileup`,
  `pileupUp`, for the selected events), `objs["object_corrections"]` (the working point
  of each object SF in the channel) and `objs["year"]`. The `corrections_study`
  collection in `hist_collections.yaml` is built on them; its histograms are defined
  at the end of `sidm/definitions/hists.py`.
* **Condor**: `condor/run_job.sh` takes the corrections as its 7th argument, `syst` as
  the 8th (weighted histograms with the `syst` axis; the default is still unweighted)
  and the monitored corrections as the 9th; see [`condor/README.md`](../../../condor/README.md).
* **Metadata sidecars:** `metadata.write_run_metadata(..., corrections=True)` records
  each correction's file and its md5 checksum; the Condor merge script forwards what
  the processor recorded, monitored corrections included.
* **Your own correction** (e.g. a trigger or DSA-muon SF measured in the analysis):
  write it as a correctionlib file, register it in `corrections.yaml` with the
  selection cut it corrects, and name it in `corrections=[...]`. The tutorial,
  section 9, does this end to end.

## The Condor campaign

Notebooks 01–03 read the merged outputs of one campaign, so that anyone can rerun them
without processing anything.

**What ran.** `condor/submit_corrections.sub`: `SidmProcessor` with the channels
`base`, `2mu2e`, `4mu`, `2mu2e_sr`, `4mu_sr` and `Z_to_mumu_control_region`, the
`corrections_study` histogram collection, `corrections=True`,
`monitor_corrections=["photon_pixel_veto"]` and `syst_hists=True`, over

* every v10 signal point (`signal_2mu2e_v10.yaml`, `signal_4mu_v10.yaml`,
  `llpNanoAOD_v2`): 180 points, 8,799 files, 25.7M events;
* every skimmed background (`backgrounds.yaml`, `skimmed_llpNanoAOD_v2`): DY→μμ and DY→ℓℓ
  (both mass ranges), the 12 QCD μ-enriched bins, $t\bar t$, WW, WZ, ZZ; 69,990 files
  (the 240 files the background census flagged are left out), about 250M events.

The jobs also had two speed-ups of `SidmProcessor` that are not part of this code (they
come in a separate change): channels with the same object cuts sharing one lepton-jet
clustering, and a vectorised mask in `make_constituent`. They change no output: with the
campaign's options, this code and the code of the jobs give identical histograms and
cutflows on a 4Mu signal file, so the outputs below are what this code produces. Without
them the processing takes about 3 times longer, so the job times quoted below are lower
than a rerun with this code would see.

**Where the outputs are.** One `<sample>.coffea` and its `<sample>.meta.yaml` sidecar
per sample in

```
/store/user/scampbel/sidm_coffea_outputs/corrections_v1/
```

which is readable by everyone. (The `lpcmetx` group area would be the natural home,
but its EOS ACL does not give this account write access; once it does,
`xrdcp -r` the directory there and change `corrections_study.CAMPAIGN_DIR`.) In a
notebook:

```python
from sidm.studies.corrections import corrections_study as study
out = study.load(["4Mu_500GeV_1p2GeV_1p9mm", "DYJetsToMuMu_M50"])   # copies once into outputs/
out["4Mu_500GeV_1p2GeV_1p9mm"]["cutflow"]["4mu_sr"].variation_yields()
```

A Kerberos ticket is enough to read them on LPC. `study.load` also takes the
generator-weight outlier of `DYJetsToMuMu_M50` out of the normalisation (see
[found along the way](#found-along-the-way)); a plain `coffea.util.load` does not.

**How it was run** (to rerun it, or to run another collection the same way):

```bash
# job lists: files per job chosen for ~100k events per job
python3 condor/make_job_args.py --samples-file condor/samples_corrections_2mu2e.txt \
    --dataset llpNanoAOD_v2 --location-cfg signal_2mu2e_v10.yaml --files-per-job 25 \
    --replace-xcache --outdir condor/filelists_corrections --job-args condor/job_args_corrections_2mu2e.txt
#   (the same for 4mu; backgrounds with --dataset skimmed_llpNanoAOD_v2 --location-cfg backgrounds.yaml
#    --files-per-job-map condor/files_per_job_corrections.yaml --census-skip backgrounds_skimmed.skip.json)
# package the code (condor/README.md section 8), then from condor/:
condor_submit JOB_ARGS=job_args_corrections_2mu2e.txt submit_corrections.sub   # and 4mu, background
# then, from the repository root, until nothing is left:
python sidm/scripts/corrections_campaign.py status
python sidm/scripts/corrections_campaign.py retry --submit
python sidm/scripts/corrections_campaign.py skipped --submit
python sidm/scripts/corrections_campaign.py merge
```

`corrections_campaign.py` reconciles the job lists against the Condor logs and the
chunk outputs (`condor/condor_campaign.py`), resubmits failed chunks, finds input files
that coffea skipped, and merges every complete sample.

**What it took.** 2,955 jobs and six resubmission rounds. 99% of the chunks were done in
six hours; the rest took another six, a tail of slow chunks on slow worker nodes and an
hour-long EOS outage.

* **13 corrupt signal files** failed in every attempt, all 4Mu: one to four branches
  each do not decompress (`event`, `GenJet_eta`, `LHEPdfWeight` and a few others). They
  are dropped; the list, with the chunk each was in, is in
  [`campaign/bad_files.txt`](campaign/bad_files.txt). These branches are new to the
  analysis (`event` seeds the JER smearing, `GenJet` is the generator match of the
  jets), which is why earlier campaigns and the file census did not trip on them; the
  census should be extended to read them.
* **4 signal files** have no events and are skipped by coffea
  ([`campaign/empty_files.txt`](campaign/empty_files.txt)).
* No other input file was skipped.
* A handful of jobs were slow on a few worker nodes (up to six times slower than
  elsewhere) and were restarted on other nodes. Four chunks that stayed stuck were
  processed on an interactive node with the code of the job tarball and the same options
  (`campaign/local_chunks.txt`). The submit file now holds and restarts a job that runs
  for more than twelve hours (six during the campaign, raised for the slower processing
  without the speed-ups).
* The tail of the campaign was the chunks of 300 tiny files (low-pT QCD, DY→ℓℓ
  M10–50, a few events per file): each file is a separate work item with a fixed cost
  of about 30 seconds, so these took 2 to 3.5 hours while 100k-event chunks of large
  files took 10 to 40 minutes. Next time, cap the files per job at about 100
  (`condor/files_per_job_corrections.yaml`).

## What changed in the repository

| file | change |
|---|---|
| `sidm/configs/corrections.yaml` | **new.** Per year: which POG file and correction to use for each correction, the pinned snapshot, the selection cut an object SF corrects, the names of its variations, and a validity note |
| `sidm/tools/corrections.py` | **new.** `Corrections` (loads the config and the JSON files, evaluates event and object corrections, records provenance), `EventWeights` (the event weight as a product of named factors with up/down variations, and the weight categories), jet and b-tagging helpers for studies, and `import_correctionlib()` (see [environment](#environment)) |
| `sidm/tools/test_corrections.py` | **new.** Standalone test in the style of `test_gen.py`: 49 checks, with synthetic corrections and, when CVMFS is mounted, the pinned POG files. `python sidm/tools/test_corrections.py` |
| `sidm/tools/sidm_processor.py` | `corrections=`, `syst_hists=` and `monitor_corrections=` options: builds the event weights for simulation, adds each channel's object SFs, carries every weight category through the event selection, adds the `syst` axis, passes the factors to fill functions, and records the corrections in the metadata. |
| `sidm/tools/selection.py` | the cutflow rows also record the weighted count of every category; a cut listed twice in a selection is still reported and skipped as before |
| `sidm/tools/cutflow.py` | `SimpleCutflow` rows can carry `variations`, which add, scale and survive `identity()`; new `variation_yields()` |
| `sidm/tools/histogram.py` | optional `syst` axis: `make_hist(..., systematics=[...])`, and `fill()` takes `{category: weights}` |
| `sidm/definitions/hists.py`, `sidm/configs/hist_collections.yaml` | **new histograms and collections** `corrections_weights`, `corrections_objects`, `corrections_applicability`, `corrections_study`: pileup and vertex distributions, the weight factors, the SF of every lepton-jet constituent, the inputs of the applicability study, and the lepton-jet isolation recomputed with raw, JES-shifted and JER-smeared jets |
| `sidm/definitions/objects.py` | `pileup` (simulation only) and `event_number` as standard objects; `genPhotons` as a derived object |
| `sidm/tools/metadata.py` | `write_run_metadata(..., corrections=..., monitor_corrections=...)` records the corrections, with file checksums |
| `sidm/scripts/merge_coffea_chunks_eos.py` | **bug fix**: the per-chunk normalisation correction was never applied to real chunk files ([found along the way](#found-along-the-way)). Also: accumulates one chunk at a time (a sample of 1,000 chunks no longer needs 16 GB of memory), frees each sample's local chunks once merged, `--samples` to merge a subset, and forwards the corrections to the sidecar |
| `sidm/scripts/corrections_campaign.py` | **new.** Status, retry, skipped-file and merge steps of the campaign above |
| `condor/run_sidm_chunk.py` | reads `skim_factor`, `is_data` and `year` from the location YAML of the files (it used `skim_factor = 1` for every sample); refuses collision data unless `--allow-data`; `--corrections`, `--syst-hists`, `--monitor-corrections` |
| `condor/run_job.sh` | optional 7th (corrections), 8th (`unweighted`, `weighted` or `syst`) and 9th (monitored corrections) arguments; the defaults run as before |
| `condor/make_job_args.py` | `--files-per-job-map`: files per job per sample |
| `condor/submit_corrections.sub`, `condor/samples_corrections_*.txt`, `condor/job_args_corrections_*.txt`, `condor/files_per_job_corrections.yaml` | **new.** The campaign above |
| `requirements.txt`, `condor/requirements.txt` | `correctionlib==2.9.0` pinned as a direct dependency (it was already in both constraints files) |
| `sidm/scripts/install_correctionlib_cmssw.sh` | **new.** Builds correctionlib into `.python_userbase` for the CMSSW_17 kernel (see [environment](#environment)) |
| `README.md`, `condor/README.md` | how to turn the corrections on; the new job arguments; the merge fix |
| `.gitignore` | the campaign's filelists and logs |
| `sidm/studies/corrections/` | **new.** This README, four notebooks, `corrections_study.py` (helpers shared by the notebooks, including the campaign loader), `tables/`, `campaign/` and `slides/` |

Nothing in the default path changed: with `corrections=None`, a 4Mu signal file and a
DY file give identical cutflows, histograms (values, variances and flow bins), metadata
and warnings on this branch and on `main`, in weighted and unweighted mode.

## The corrections

All files come from `/cvmfs/cms-griddata.cern.ch/cat/metadata/<POG>/Run2-2018-UL-NanoAODv9/<date>/`,
the CMS analysis-corrections area that replaced jsonpog-integration.

| name | applied by default | source (snapshot) | what it corrects | validity for our objects |
|---|---|---|---|---|
| `pileup` | yes | LUM `puWeights` (2021-09-10) | pileup profile of the simulation | standard |
| `l1_prefiring` | yes | NanoAOD `L1PreFiringWeight_Nom/Up/Dn` | L1 prefiring; in 2018 muons only | standard; DSA-only lepton jets do not enter it |
| `muon_id` | yes | MUO `muon_JPsi` below 20 GeV, `muon_Z` above (2024-07-02) | PF muon `looseID` (or `tightID`, picked from the channel's cuts) | approximate: measured on prompt, isolated muons |
| `electron_reco` | yes | EGM `electron`, RecoBelow20/RecoAbove20 (2024-07-02) | GSF electron reconstruction | approximate: measured on prompt electrons |
| `photon_pixel_veto` | **no** (monitored in the campaign) | EGM `photon`, PixVeto, Loose (2024-07-02) | photon pixel-seed veto | not applicable: our "photons" are mostly electron pairs |
| `electron_id` | no | EGM `electron`, Loose (2024-07-02) | cut-based Loose ID | not used by the lepton-jet selections |
| `photon_id` | no | EGM `photon`, Loose (2024-07-02) | cut-based Loose ID | not used by the lepton-jet selections |
| `jet_energy` | helper | JME `jet_jerc` (2026-06-10), `jer_smear` (2026-09-18) | JEC check, JES uncertainty, JER smearing | see below |
| `jet_veto` | helper | JME `jetvetomaps` (2026-06-10) | jet veto map, including HEM | standard |
| `btag` | helper | BTV `btagging` (2025-08-19) | DeepJet fixed-WP SFs | standard, needs MC efficiencies |

**How object SFs are applied.** Each object SF is tied to the selection cut whose
efficiency it corrects (`working_points` in the config). In a channel whose `muons`
cuts include `looseID`, every PF muon inside the channel's lepton jets (after the
lepton-jet cuts) gets the Loose ID SF, and the event weight is the product over those
muons. A channel asking for `tightID` gets the Tight SF instead. Objects outside the
measured range get an SF of one. Up and down move all the SFs of an event together.
`out[sample]["metadata"]["object_corrections"]` lists what each channel used.

**Already in NanoAODv9, not to be applied again:** the jet energy corrections
(`Summer19UL18_V5`, the current tag; notebook 00 recomputes them and finds agreement
within NanoAOD's float precision) and the electron and photon energy scale and
smearing (`eCorr` differs from one in simulation).

## Where the standard corrections stop describing our objects

From `02_applicability.ipynb`, over all 180 signal points at the `base` selection,
compared with DY muons in the Z $\to\mu\mu$ control region:

* **PF muons.** 79–100% of the lepton-jet PF muons lie outside the phase space of the
  Muon POG measurement (prompt, full pixel track, no muon nearby), against 7% of DY
  muons; 98% or more at every lifetime of the benchmark ($m_\Phi$ = 500 GeV,
  $m_{Z_d}$ = 1.2 GeV). The main reason is the partner muon: 87–100% of the benchmark
  muons have another muon within $\Delta R$ = 0.3. At $c\tau$ = 1.9 mm three quarters
  also have at most two pixel hits.
* **DSA muons** have no standard correction at all. They are in 1–23% of the
  muon-type lepton jets for a median dark-photon $L_{xy}$ below 5 cm, 24–66% between 5
  and 50 cm, and 48–88% beyond.
* **Electrons.** The analysis ID (`mvaFall17V2noIso_WPL`) has no SF; the closest
  working point with one, wp90noiso, keeps 15–87% of our electrons (57% at the
  shortest benchmark lifetime, 24–31% at the longer ones).
* **Photons.** For displaced decays, 97–100% of the lepton-jet "photons" are within
  $\Delta R$ < 0.1 of a generator electron, and they outnumber the electrons (1.7 to
  4.9 per electron at the displaced benchmark points). The photon SFs, measured with
  real photons, do not apply, which is why the pixel-veto SF (0.96 in the barrel, 0.88
  in the endcap) is not in the default list.
* **Lepton-jet isolation** is computed from the JEC-corrected energy of the matched
  AK4 jet. For the 4Mu benchmark at $c\tau$ = 1.9 mm, nine in ten of those jets are
  more than 90% leptons (46–98% over the grid), and the JEC scales them up by about 7%.
  Using the raw energy instead changes the isolation decision for 0.1–2.6% of lepton
  jets; JER smearing for up to 1.3%; the JES variation for up to 1.2%.
* **Triggers.** There is no POG SF for `DoubleL2Mu23NoVtx_2Cha` and its siblings.

The last section of notebook 02 has the full table: every requirement of the
selection, the correction that covers it, and its status.

## Effect on the results

From `03_impact_on_results.ipynb`: all 180 signal points in `2mu2e_sr` and `4mu_sr`
(full samples). Changes of the SR signal yield (min to max over the grid):

| | 2Mu2E | 4Mu |
|---|---|---|
| **all default corrections** | **−2.2 to +0.7%** (median −0.7%) | **−1.8 to +0.2%** (median −1.4%) |
| L1 prefiring | −0.7 to −0.1% | −1.2 to −0.4% |
| muon ID SF | −0.2 to −0.05% | −1.0 to −0.2% |
| electron reconstruction SF | −1.3 to +1.8% | 0 |
| pileup | −0.3 to +1.0% | −0.2 to +1.0% |
| photon pixel veto (monitored, not applied) | −5.6 to −0.05% | −0.2 to 0% |

Uncertainties (half-width of the up/down variation of the SR yield):

| | 2Mu2E | 4Mu |
|---|---|---|
| pileup | 0.0 to 2.1% | 0.0 to 2.9% |
| L1 prefiring | 0.04 to 0.16% | 0.10 to 0.29% |
| muon ID SF | 0.01 to 0.19% | 0.11 to 1.7% |
| electron reconstruction SF | 0.05 to 4.2% | 0 |
| **default, in quadrature** | **0.15 to 4.3%** (median 0.6%) | **0.2 to 3.0%** (median 0.6%) |

* With the background fixed, a cross-section limit scales as the inverse of the
  signal yield, so the expected limits move by the same fraction, mostly by about 1%.
* The extreme values are all at low mass. At $m_\Phi$ = 100 GeV the electrons are
  soft, and the reconstruction SF below 20 GeV can be above one and has a 6%
  uncertainty per electron (0.6% above), hence the positive 2Mu2E shifts and the 4.2%;
  the 4Mu muons are soft too, hence the 1.7% from the J/$\psi$-based muon SF.
* Pileup is above 1% for 12 points, ten of them at $m_\Phi$ = 100 GeV. With the full
  samples the 4Mu values are 2.5–8 times their statistical precision; the 2Mu2E ones,
  with 133–1610 SR events, only 1.4–3 times.
* The effects of prefiring and the muon SF shrink at long lifetimes, where more of the
  muons are DSA muons, which enter neither.
* The simulated backgrounds go down by 0.8–1.0% at the `base` selection (DY, QCD,
  $t\bar t$ and diboson alike) and by 1.7% (`2mu2e_sr`) and 0.5% (`4mu_sr`) in the
  SRs, whichever DY sample is used; QCD is nine tenths of the SR background. The pileup
  variations of 4–11% in the `4mu` channels come from a few events with large
  weights. The two DY samples overlap and are never added; which one the analysis uses
  is open, so notebook 03 gives the total with each.

Simulated backgrounds (change of the yield ± the default corrections' uncertainty,
from notebook 03, which also shows them per channel in `impact_backgrounds`):

| | `base` | `2mu2e_sr` | `4mu_sr` | SR yields, 2Mu2E / 4Mu (59.8 fb⁻¹) |
|---|---|---|---|---|
| DY→μμ | −0.9 ± 0.2% | −1.8 ± 0.7% | +1.2 ± 8.6% | 824 / 390 |
| DY→ℓℓ | −0.9 ± 0.2% | +0.2 ± 5.3% | +1.6 ± 10% | 586 / 438 |
| QCD | −0.8 ± 0.5% | −1.7 ± 0.6% | −0.6 ± 4.1% | 17400 / 14100 |
| $t\bar t$ | −1.0 ± 0.2% | −1.9 ± 0.8% | −0.5 ± 4.7% | 1100 / 173 |
| diboson | −1.0 ± 0.2% | −0.8 ± 2.4% | −2.9 ± 3.0% | 17 / 3.5 |
| **total, with DY→μμ** | **−0.9 ± 0.2%** | **−1.7 ± 0.6%** | **−0.5 ± 4.3%** | 19300 / 14600 |
| **total, with DY→ℓℓ** | **−0.9 ± 0.2%** | **−1.6 ± 0.5%** | **−0.5 ± 4.3%** | 19100 / 14700 |

At `base` the prefiring weight takes 0.6–0.7% from every group and the muon ID SF
0.2–0.3%; the electron SF adds −0.6 to −0.8% in the 2Mu2E channels. The large 4Mu
variations come from a few simulated events with large weights. The SR yields are
large because the SR selections do not apply their isolation and displacement cuts to
the lepton jets the categories count (found along the way).

For the datacards, each correction becomes one lnN nuisance correlated across
channels and signal points; `tables/sr_correction_nuisances.csv` has, for every point,
the SR yields with and without the corrections, and each correction's shift and
$\kappa_{\rm up}$/$\kappa_{\rm down}$. The limit notebooks use the ABCD selections of the
`abcd-plane-study` work rather than `2mu2e_sr`/`4mu_sr`; the numbers should be
recomputed with those selections by running the processor with `corrections=True`,
but they are driven by the lepton content and should be similar.

## Validation

* `sidm/tools/test_corrections.py`: 49 checks, including a zero factor that must not
  produce NaN variations, the pT stitching of the muon SFs, the working-point
  selection from cuts, cutflow variations through `add`/`scale`/`identity`, and the
  `syst` axis.
* With `corrections=None` the branch reproduces `main` exactly (see above).
* With the campaign's options, this code and the code the campaign jobs ran (which had
  the processor speed-ups, see [the Condor campaign](#the-condor-campaign)) give identical
  histograms and cutflows on the MC fixture and a 4Mu signal file.
* Closure (notebook 00): the event weights the processor applies are rebuilt from
  their inputs (the NanoAOD branches and the correction files) with the `Corrections`
  API, to machine precision.
* Consistency (notebook 03): in every sample, channel and weight category of the
  campaign, the integral of a one-entry-per-event histogram equals the cutflow yield,
  to $10^{-15}$.
* The merge was checked by hand on two-chunk samples (identical to the sumw-weighted
  combination) and on a single chunk merged four times (1.000×; before the fix, 4×).
* The JEC recomputed from the raw jet $p_T$ agrees with NanoAOD (notebook 00).
* The correctionlib import was tested under CMSSW_17_0_0_pre4 (the kernel) and under
  CMSSW_14_1_0_pre4 with the PyPI wheel (the Condor jobs).

## Environment

The analysis runs under `cmsenv` in CMSSW_17_0_0_pre4 (Python 3.12), with packages in
`.python_userbase`. Under `cmsenv`, the PyPI correctionlib wheel cannot be imported:

```
ImportError: .../correctionlib/_core.cpython-312-x86_64-linux-gnu.so: undefined symbol: __cxa_call_terminate
```

The wheel's compiled core finds its bundled `libcorrectionlib.so` through a RUNPATH,
but cmsenv puts CMSSW's own, older copy on `LD_LIBRARY_PATH`, which the dynamic
loader searches first. Two fixes are in place:

* `bash sidm/scripts/install_correctionlib_cmssw.sh` (run once, after `cmsenv`)
  rebuilds correctionlib 2.9.0 from source into `.python_userbase` with CMSSW's
  compiler and zlib and an RPATH, so a plain `import correctionlib` works in the
  kernel. This has been done for the current `.python_userbase`.
* `sidm/tools/corrections.py` loads correctionlib through `import_correctionlib()`,
  which binds the core to its own library first. That is what makes the Condor jobs,
  which install the wheel under CMSSW_14, work without the rebuild.

The LCG_107 venv described in the top-level README is not affected. `pylint` 4.0.6
was also installed into `.python_userbase` to run the CI lint locally.

## Open items

In rough order of how much they matter for the limits:

1. **Trigger efficiency.** No standard SF exists for the displaced-dimuon triggers;
   it has to be measured (orthogonal triggers, or the cosmic-seeded paths) and
   entered as a correction of our own.
2. **DSA muons**: identification and reconstruction efficiency, e.g. from cosmics or
   standalone-track tag-and-probe.
3. **Displaced and close-by PF muons, displaced electrons**: the applied SFs are
   provisional and their quoted uncertainties (per mille level) describe prompt,
   isolated objects. Until a dedicated measurement exists, they need an uncertainty
   of their own.
4. **Electron ID (WPL) and the photon ID and pixel veto** for lepton-jet objects.
5. **Jets and lepton-jet isolation**: decide whether the isolation should use the raw
   jet energy; then propagate JES and JER, which needs a rerun of the selection per
   variation (the processor varies only weights today).
6. **Data/MC validation** of the corrections that are applied (for pileup, the
   vertex multiplicity in a control region), once the group decides which regions can
   be looked at.
7. **Electron and photon energy scale**: in these ntuples the scale variations
   (`dEscaleUp/Down`) are zero for every object in simulation; the smearing
   variations are filled. The scale uncertainty has to come from elsewhere.
8. **b-tagging** in the $t\bar t$ control regions: the SFs are available
   (`Corrections.btag_sf`), the MC tagging efficiencies needed for an event weight are
   not measured yet.
9. **Muon momentum scale** (Rochester corrections): not in correctionlib for Run 2,
   small at our momenta.
10. **Datacards**: add the correction uncertainties as lnN nuisances in
    `limit_plotting/datacard_tools.py` from `tables/sr_correction_nuisances.csv`.
11. **Run 3**: `corrections.yaml` has only 2018; the other years need their own block
    with the matching POG files.
12. **Generator-weight outliers**: decide how the analysis treats the
    `DYJetsToMuMu_M50` event below (drop it in the skim, or cap the weights), and add a
    weight check to the file census.

## Found along the way

* **The Condor merge inflated weighted yields by the number of chunks.**
  `merge_coffea_chunks_eos.py` re-weights each chunk by its sum of weights before
  adding them, but it looked for the metadata at the top of the chunk file, while
  `run_sidm_chunk.py` saves what `Runner.run` returns,
  `{"out": {sample: ...}, "processed": ..., "exception": ...}`. The re-weighting was
  therefore skipped, and every merged weighted cutflow and weighted histogram came out
  multiplied by the number of chunks of its sample (exactly 4.000× for four copies of
  one chunk). Raw counts and unweighted histograms were right. Fixed; **existing merged
  outputs should be re-merged** (the chunks are enough). This includes, if they were
  made with this script, the `*Merged*` directories under
  `/store/user/scampbel/sidm_condor/`.
* **Skimmed backgrounds were normalised without their skim factor on Condor.**
  `run_sidm_chunk.py` set `skim_factor = 1` for every sample, which overestimates a
  skimmed sample's yields by 1/skim_factor (e.g. 150× for `DYJetsToLL_M10to50`). It now
  reads the skim factor from the location YAML. Condor outputs of skimmed backgrounds
  made before this branch need to be reprocessed, or rescaled by the skim factor.
* **One `DYJetsToMuMu_M50` event has a generator weight of 1.45e14** (the LHE weight
  too; the other events have about 2400): event 147047969, luminosity block 390048, in
  `/store/group/lpcmetx/SIDM/Backgrounds/2018_v2/Skims/DYJetsToMuMu_M-50_massWgtFix/skimmed_output_6129.root`.
  It alone is 99.8% of the sample's sum of weights, so any normalisation of this sample
  by its genWeight sum comes out about 630 times too small. It fails the trigger, so it
  is in no selected yield, and taking it out of the normalisation is exact;
  `corrections_study.load` does that (`GENWEIGHT_OUTLIERS`), and a scan of every chunk
  of the campaign found no other such event.
* The SR selections `2mu2e_sr` and `4mu_sr` apply their isolation and displacement
  cuts to `mu_ljs`/`egm_ljs`, but the event cuts that define the categories count
  `ljs`, which those cuts do not slim, so on this branch the SR differs from
  `2mu2e`/`4mu` only by the $\Delta\phi$ cut (notebook 03, section 6). The limit
  notebooks use different SR selections, from the `abcd-plane-study` work.
* **Integer quantities on bin edges.** A `hist.axis.Regular` whose edges fall on the
  integers computes the bin of an integer value by a division that can round down:
  with 100 bins on 0–100, 29, 57 and 58 land one bin low, and `pv_n` (50 bins on
  0–100) puts 58 vertices in the bin of 56–57. `pileup_nTrueInt` had this in the
  campaign (notebook 01 draws it in bins of two, which keeps 28 and 29 together) and is
  now centred on the integers (−0.5 to 99.5); other integer histograms are worth
  centring the same way.
* `lj_lj_absdphi` and the ABCD-plane histograms built on it
  (`lj_lj_absdphi_invmass` and its siblings in `hists.py`) fill
  `abs(phi1 - phi0)` without wrapping, so they run from 0 to 2π and a value $x$ above
  π is the same configuration as $2\pi - x$. A cut on that axis, such as
  $|\Delta\phi| > 2$, would also keep pairs with a true $|\Delta\phi|$ below 2. The
  event cut `LJ-LJ dPhi > 2` uses the wrapped `delta_phi` and is not affected.

## Notebooks

| notebook | what it shows | input | runtime |
|---|---|---|---|
| [`00_correctionlib_tutorial.ipynb`](00_correctionlib_tutorial.ipynb) | correctionlib from scratch: files on CVMFS, inspecting, evaluating, flow and variations, compound corrections, reproducible smearing, writing and registering a correction of our own, the processor options and weight categories, and the closure of the weights | the MC fixture in `tests/data` | 1 min |
| [`01_standard_corrections.ipynb`](01_standard_corrections.ipynb) | each default correction on signal and DY: pileup weights, prefiring, object SFs, and the corrections already in NanoAOD | campaign outputs | 2 min the first time |
| [`02_applicability.ipynb`](02_applicability.ipynb) | where the standard corrections stop describing lepton-jet objects, over the whole grid, with the full coverage table | campaign outputs | 1 min |
| [`03_impact_on_results.ipynb`](03_impact_on_results.ipynb) | effect of the corrections and their uncertainties on the SR yields of every signal point and on the simulated backgrounds; writes `tables/sr_correction_nuisances.csv` | campaign outputs | 2 min |

Notebooks 01–03 make every plot from the `corrections_study` histograms and the
cutflows: no event loop, no processor run. The campaign outputs are copied into
`outputs/` (gitignored) on first use; delete it to fetch them again, or set
`study.CAMPAIGN_DIR` to read another production.

Every histogram in the notebooks is a `hist.Hist`, drawn with `mplhep` (`histplot`,
`hist2dplot`) in the CMS style, and every figure carries the CMS Simulation header
through `corrections_study.cms_label()` (with the 2018 luminosity on yield plots).
Distributions on log-binned axes are shown as fractions per bin
(`corrections_study.fractions()`), since a density there would divide by the linear
bin width. `corrections_study.save()` writes each figure as a PDF into
`slides/figures/`, where the deck picks it up.

## Slides

`slides/corrections.tex` is a summary of this study for the group: the framework, the
campaign, the status of each correction, the applicability study, the effect on the
SR yields, the findings along the way and the open items, with the notebooks' figures.
Like the other decks in the repository, the PDF and the figures are not committed
(`*.pdf` is gitignored), so rebuild it with:

```bash
# run the four notebooks first: they write slides/figures/*.pdf
cd sidm/studies/corrections/slides
pdflatex -interaction=nonstopmode corrections.tex
pdflatex -interaction=nonstopmode corrections.tex   # second pass for the frame references
```
