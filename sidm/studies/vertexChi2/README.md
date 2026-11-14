# Dimuon-vertex-fit cosmic veto, 2018C

This study compares two cosmic-muon vetoes for the lepton jets (LJs) of the 1-LJ control
region, on 2018C `DoubleMuon` data, the 2018 background stack and the v10 signal grid:

- the **vertex veto**, new here. It keeps a mu-type LJ (an LJ that contains muons) only if
  `0 <= vtx_chi2 < 5`, where `mu_ljs.vtx_chi2` is the smallest normalized chi2 of any stored
  dimuon vertex built from two of the LJ's own muons;
- the **Spread veto**, the `*Spread_*` veto already used in the 1-LJ control region: four cuts
  on how far apart the LJ's muons are in impact parameter and vertex position, together with the
  **cosAlpha** event cut `all cosAlpha(mu, mu) > -0.99`, which rejects back-to-back muon pairs
  (`cosAlpha` is the cosine of the opening angle between two muons).

It then characterizes the **tail**: the DSA (displaced standalone) muons in mu-type LJs with
pT >= 200 GeV that survive either veto.

**Outcome.** The vertex veto tags the same cosmics as the Spread veto: it removes 98% of the
mu-type LJs in the spread-inverted cosmic region (defined below). It keeps 0.690 of the signal
LJs against 0.961 for the Spread veto, and it leaves the same data/MC above 200 GeV in DSA-muon
pT (2.38 against 2.31 on the full background stack). Spread + cosAlpha, `main`'s baseline since
PR #56, therefore stays, and the vertex veto is not adopted. The tail that survives the Spread
veto stays as sections 5-8 of the write-up leave it: not cosmic (below 6.2% at 95% CL), and of
the size DSA momentum mismeasurement would produce at roughly factor-1.5 precision. Section 10
tests a transfer factor at full-2018 statistics: it takes the DSA mismeasurement rate of a
control sample of DSA muons matched to PF muons and folds it (spreads each low-pT bin over the measured-pT bins with the
control's migration probabilities) onto the
lepton jets' own low-pT DSA muons to predict the tail. The same prediction made entirely in MC
does not reproduce the MC tail (it fails its MC self-closure), and nothing from it is applied.

**Start here: [`vertexChi2CosmicVetoStudy.ipynb`](vertexChi2CosmicVetoStudy.ipynb).** It is the
write-up and reads only the `.coffea` outputs of the run notebooks below and the three JSON
records in `transfer_factor/`; it runs no processor of its own.

## Selections used in this study

Every channel of this study is defined in
[`sidm/configs/selections_vertexChi2_as_run.yaml`](../../configs/selections_vertexChi2_as_run.yaml),
with each cut list written out in full; the write-up renders it as a table in its section
"Selections used in this study".

**Where the outputs come from.** The outputs were produced between 2026-08-25 and 2026-09-05 on
the fork `murtaza-safdari/SIDM`, branch `vertex-cosmic-veto-data2018` (kept for provenance).
That was before the cosmic-veto baseline of PR #56 was merged into `main`. The commits named
here are on that branch. The 2018C, background and signal runs did not record their commit:
their sidecars were written afterwards and carry `sidm_commit: null`. The four full-2018
response runs record commit `f6e23f1`, from a working tree that already contained the histogram
collections committed in `5aa6988`.

**The channel definitions were the same for all runs.** `selections.yaml` did not change on that
branch between `defc8a0` and its tip `26c6d1c`, and the run-time sidecars of the four full-2018
directories record the same definitions. The six back-filled sidecars (see "Provenance
sidecars") take their channel list from each output's cutflow and the definitions from that
branch's `selections.yaml` at the time they were written. They therefore agree with it by
construction and confirm nothing on their own. The stored-output reruns described under
"Notebooks" confirm the definitions independently for one sample of each of two directories, and
for a third through its ZZ sample.

**Why the channels are kept as run.** The study compares cosmic-veto options. It therefore needs
a baseline with no cosmic veto at all, and each veto applied on its own to that same baseline.
On `main` the baseline itself has changed, so `main`'s channels cannot provide the no-veto
reference, and the channels are kept here as they were run. How the as-run baseline differs
from `main` today:

- Unchanged: the eta-phi veto on all four LJ sources (electrons, photons, PF muons, DSA muons)
  and the event cut `>=2 muons (pf or dsa) pT > 26 GeV`.
- As run, the PF muons at the LJ source have no explicit DT + CSC hit-count requirement. `main` now also
  requires `DT + CSC hits > 3` on them (`trkNumDTHits + trkNumCSCHits > 3`), added in commit
  `65b46a86` to the `muons_ljsource_cuts` anchor of `selections.yaml`, which now reads
  looseID, pT > 5 GeV, |eta| < 2.4, eta-phi_veto, DT + CSC hits > 3. The as-run baseline lacks
  this cut, and its effect on the veto comparison was not measured.
- As run, the base event cuts do not contain `all cosAlpha(mu, mu) > -0.99`; only the channels
  that use it as part of a veto apply it. On `main` it is in the base event cuts. `main`'s
  `cr_1egmOr1mu_lj` and `base_ljDisplacementIso` therefore include it, and `main`'s
  `cosmic_muons` requires both it and its inverse `any cosAlpha(mu, mu) <= -0.99`.
- As run, the muon LJs of `base_ljDisplacementIso` carry no spread cuts. On `main` they require
  `dzSpread_mu <= 50`, `dxySpread_mu <= 50`, `dzSpread_pf <= 10` and `dxySpread_pf <= 10`.
  `cr_1egmOr1mu_lj` has no spread cuts in either version.
- The Spread veto as run (`cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu`,
  `base_ljDisplacementIso_cosAlpha_ljSpread`) cuts on `vxySpread`, the transverse spread of the
  muon vertex positions; `main`'s spread cuts use `dxySpread`, the transverse impact-parameter
  spread. Both use `dzSpread` for the longitudinal cuts. The `dxySpread` version was also run,
  as `cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu_dxy`, on data and the background stack only.
  The signal grid was run with the `vxySpread` form alone, so every signal efficiency of the
  study is for that form; `main`'s `dxySpread` form was not run on signal.
- `inverse_spread` is the mu-LJ cut that defines the spread-inverted cosmic region
  `cosmic_muons`. At production time it was the OR of `vzSpread_pf > 10`, `vxySpread_pf > 10`
  and `vxySpread_mu > 50`, the last term written twice. `main` has since redefined that name as
  the OR of `dzSpread_pf > 10`, `dxySpread_pf > 10`, `dxySpread_mu > 50` and
  `dzSpread_mu > 50`. The production-time definition is kept in `sidm/definitions/cuts.py`
  under the name `inverse_spread (vz/vxy, as run)`, which the as-run file uses; the run-time
  sidecars list it as `inverse_spread`. It is an object cut, so it is not a cutflow row. As
  run, `cosmic_muons` is therefore only approximately the inverse of the Spread veto: it uses
  `vzSpread_pf` where the veto uses `dzSpread_pf`, and it has no longitudinal muon term, so a
  lepton jet that fails only `dzSpread_mu <= 50` is in neither region.
- Four of the channel names also exist in `main`'s `selections.yaml`. `cr_1egmOr1mu_lj`,
  `cosmic_muons` and `base_ljDisplacementIso` select different events there.
  `cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu_dxy` applies the same cuts there in a different
  order, and `main`'s version lists `all cosAlpha(mu, mu) > -0.99` twice, once through its base
  event cuts and once explicitly. The other channels exist only in the as-run file.

## Two cosmic-enriched regions, and why they are never interchangeable

The study uses two different cosmic-enriched selections. They are named apart everywhere (in
the write-up prose, in the figure legends and here), because they are not the same region and
the difference decides what one of the results is conditional on.

| name used throughout | channel | mu-LJ cuts | uses `vtx_chi2`? |
|---|---|---|---|
| **spread-inverted cosmic region** | `cosmic_muons` | `lj_iso < 0.2`, `inverse_spread (vz/vxy, as run)` (see above), cosAlpha event cut inverted | no |
| **vertex-inverted cosmic region** | `cosmic_muons_vtxChi2` | `lj_iso < 0.2`, `no vtx or vtx_chi2 >= 5`, cosAlpha event cut inverted | yes |

The spread-inverted region makes the section-2b comparison a test of `vtx_chi2` rather than a
tautology: the variable is asked to flag objects selected with entirely different quantities.
The vertex-inverted region is where the azimuthal cosmic template of section 7 is measured, so
the 95% CL residual-cosmic bound is conditional on it. The write-up shows that the bound moves
only from 6.2% to 7.8% even if the template's vertical fraction is degraded from 0.99 to 0.90.

## Notebooks

Every run notebook opens its own Dask cluster on LPC Condor, writes one `.coffea` per sample
to `<repo>/RunOutputFiles/<run_tag>/` together with a `.meta.yaml` provenance sidecar, and
copies both to `root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/coffea_outputs/$USER/<run_tag>/`.
`sidm_path` is derived from the working directory, so run them from this directory. Per-sample
outputs already on disk are skipped, so an interrupted run can be resumed.

Every stored histogram output except `vtxGenDr_signal_2018` predates the correction of `dR_outer`
in `sidm/tools/utilities.py` (PR #65,
merged 2026-11-02: the outer-track phi difference is wrapped across +-pi, and the (-5, -5)
coordinates of a missing outer track count as missing). `dR_outer` enters the DSA cross-cleaning
cut `all` of every channel here and the segment-matched control of the write-up's section 10.
Rerunning `vtxResp2_2018Full/DoubleMuon_2018D_22` and
`vtxChi2_signal_2018/2Mu2E_200GeV_0p25GeV_10p0mm` from their sidecars with the as-run selections on
the corrected code reproduced every histogram and cutflow of the stored outputs exactly, as did the
code with `dR_outer` as it was before the fix: the 1 (data) and 4 (signal) DSA muons whose
cross-cleaning decision the wrap flips are below the DSA pT cut, matched to a PF muon that fails the
loose ID, or removed through another match, and no pair has both outer tracks missing. A third
rerun does change: in chunk 0 of `vtxLxyDr_signal_2018/2Mu2E_500GeV_1p2GeV_1p9mm` (ten files) the
corrected `dR_outer` removes 2 of the 12,905 DSA muons in muon LJs and moves 2 of the 21,857 muon
LJs to another bin, with every cutflow count unchanged, and the code with `dR_outer` as it was
before the fix reproduces the stored histograms exactly. So the stored outputs can differ from a
rerun at the 10^-4 level in the DSA-muon content. Nothing else has been rerun.
`vtxGenDr_signal_2018` (below) was produced with the correction, over the samples, files,
channels and selections of `vtxLxyDr_signal_2018`. In the reference channel
(`base_ljDisplacementIso_cosAlpha`, the selection without a cosmic veto) it holds 68 fewer
mu-type LJs with >= 2 muons than that run, of 5,764,154 (57 of the 180 samples differ; summed
over the samples, 82 fewer LJs pass, 7 more have no vertex and 7 more have chi2 >= 5), and 63
fewer of 5,581,391 in the Spread-veto channel (`base_ljDisplacementIso_cosAlpha_ljSpread`, 53
samples). These differences were not traced LJ by LJ. The per-LJ records of `dr_truth_match/`
(below) were also produced with the correction: their LJ counts equal those of
`vtxGenDr_signal_2018` in each of their 16 samples.

Every stored output but `vtxGenDr_signal_2018` also predates the standard corrections (PR #58,
merged 2026-11-12), which
`SidmProcessor` applies to simulation unless it is given `corrections=None`. The run notebooks do
not set that argument, so a simulation notebook rerun on main weights its cutflows and weighted
histograms with the corrections; `corrections=None` gives the weights of the stored outputs. Data
are not affected. The condor submission of `vtxGenDr_signal_2018` passes `none` for the
corrections and books unweighted histograms, as that of `vtxLxyDr_signal_2018` does.

Every run notebook passes `selections_cfg="configs/selections_vertexChi2_as_run.yaml"` to
`SidmProcessor` and to `write_run_metadata`. That argument was added after the notebooks were
executed, when these channels were still defined in `sidm/configs/selections.yaml`; the stored
outputs are those of the original runs. With the argument, the notebooks rerun on main with the
selections they were run with. Three stored outputs were rerun this way: one 2018D data sample
(`DoubleMuon_2018D_22`, one of the 68, of the full-2018 segment-matched run `vtxResp2_2018Full`), one signal
sample (`2Mu2E_200GeV_0p25GeV_10p0mm` of `vtxChi2_signal_2018`) and the ZZ background (`ZZ` of
`vtxChi2_vs_spread_2018C`). The data and signal histograms and cutflows were reproduced bit for
bit. For the simulated ZZ background the raw counts are identical and the weighted sums agree to
single-precision rounding (largest relative difference 2.5e-7), consistent with a different
order of summation.

`cosmic_veto` in `sidm/configs/hist_collections.yaml` is as on `main`. The nine histograms this
study added to it (the `vtx_chi2` distributions and the 2D DSA-muon pT histograms of the tail)
are in a separate collection, `mu_lj_vtx_chi2_dsa_pt`. Every run notebook requests both
collections, so a rerun books every histogram its stored outputs contain. In the table below,
`cosmic_veto` means the collection as it was when each output was produced, with these
histograms as far as they existed then.

| notebook | what it runs | `run_tag` / output directory |
|---|---|---|
| `runDataBg2018.ipynb` | 2018C `DoubleMuon` plus the background stack (TTJets, DYJetsToLL M10to50 and M50, DYJetsToMuMu M10to50, eleven QCD `Pt` bins, WW/WZ/ZZ) through six channels: Before, the `vxySpread` and `dxySpread` variants of the Spread veto with cosAlpha, the vertex veto with and without cosAlpha, and the vertex-inverted cosmic region. Hist collection `cosmic_veto`. | `vtxChi2_vs_spread_2018C` (19 samples) |
| `runDataCosmic2018.ipynb` | 2018C `DoubleMuon` in the `cosmic_muons` channel only: the **spread-inverted** cosmic region, defined by inverted spread cuts and inverted cosAlpha, which uses no vertex information. Hist collection `cosmic_veto`. | `vtxCosmic_2018C` (1 sample) |
| `runDataTail2018.ipynb` | 2018C `DoubleMuon` through Before, Spread + cosAlpha, vertex + cosAlpha and `cosmic_muons_vtxChi2` (the **vertex-inverted** cosmic region), with hist collections `cosmic_veto` and `dsa_pf_response`. Adds the 2D DSA-muon pT versus nHits, \|eta\|, phi, \|dxy\| and sigma(pT)/pT histograms, and the matched DSA-to-PF response. | `vtxTail2_2018C` (1 sample) |
| `runMCTail2018.ipynb` | The MC counterpart of the tail check, same two hist collections, Before and Spread + cosAlpha only: the eight samples that carry the >= 200 GeV DSA tail (DYJetsToLL_M50 and QCD_Pt80To120 through QCD_Pt800To1000). | `vtxTailMC2_2018` (8 samples) |
| `runSignalVtx.ipynb` | The full v10 signal grid, 180 samples (90 4Mu, 90 2Mu2E), through `base_ljDisplacementIso`, `+ cosAlpha`, `+ Spread veto` and `+ vertex veto`. Hist collection `cosmic_veto` as it stood before the 2D nHits and \|eta\| histograms were added, so those are absent from these outputs. | `vtxChi2_signal_2018` (180 samples) |
| `runSignalVtxNhits.ipynb` | A 26-sample subset of the same grid (13 4Mu, 13 2Mu2E, mirrored mass and lifetime points), re-run after `mu_lj_dsaMuon_pt_nHits` and `mu_lj_dsaMuon_pt_abseta` joined the `cosmic_veto` collection. Used only for the signal cost of an nHits cut. | `vtxChi2_signal_2018_nhits` (26 samples) |
| `runDataFull2018.ipynb` | All 68 2018A-D `DoubleMuon` skims through the five cosmic-veto channels, with hist collections `cosmic_veto`, `dsa_pf_response` and `dsa_response_matrix` (the dR-matched DSA-to-PF response). Feeds the dR-matched comparison numbers of section 10. | `vtxResp_2018Full` (68 samples) |
| `runMCRespFull.ipynb` | The MC counterpart of `runDataFull2018.ipynb`: the same eight background samples used from section 5 on, same hist collections. | `vtxRespMC_2018` (8 samples) |
| `runDataSegMatch2018.ipynb` | All 68 2018A-D `DoubleMuon` skims again, this time with `dsa_response_segmatch` (control muons matched by the segment-based cross-cleaning the analysis actually uses, not `dR`) and `dsa_tail_diag` (beamspot displacement, DT/CSC hit split, phi x \|eta\|). The primary input to section 10. | `vtxResp2_2018Full` (68 samples) |
| `runMCRespSegMatch2018.ipynb` | The MC counterpart of `runDataSegMatch2018.ipynb`, same eight background samples, segment-matched response and halo diagnostics. | `vtxRespMC2_2018` (8 samples) |
| `vertexChi2CosmicVetoStudy.ipynb` | The write-up. Reads the ten output directories above and the three JSON records in `transfer_factor/`, and recomputes every figure and printed number from them. Some values quoted in the section-10 text are not in the committed records: the dR-matched values of 10.1 and 10.2 and the no-veto column of the 10.1 table, the ratios, fractions and significance quoted in 10.3, and in 10.4 the no-veto, MC `\|eta\|` and control-muon fractions, the azimuth significance, the fit probability and the halo-corner fractions. They are in the fold logs, `transfer_factor/logs/`, and, for the MC `\|eta\|` fractions, come from the `vtxRespMC2_2018` outputs. The fold scripts that wrote the JSON records are in `transfer_factor/`; its README says how to rerun them. | reads only; writes nothing |
| `runSignalVtxLxyDr.ipynb` | The Dask form of the `vtxLxyDr_signal_2018` run: the v10 signal grid through `+ cosAlpha` and `+ Spread veto`, hist collections `mu_lj_vtx_outcome_lxy_dr` (vertex outcome of every mu-type LJ, no vertex / normChi2 >= 5 / pass, against its muon count, the Lxy of the matched gen dark photon and the dR of its two leading muons), `mu_lj_vtx_chi2_dsa_pt` and `dsa_tail_diag`. It writes to `RunOutputFiles/vtxLxyDr_signal_2018/`, the directory the condor merge below fills. At about 18 h for the grid it was run for one sample only, `4Mu_100GeV_0p25GeV_0p02mm`, whose output was then moved by hand to `dask_check/` as the cross-check of the condor production (`vtx_outcome_lxy_dr/condor_vs_dask.txt`); the 180 samples themselves come from condor, not from this notebook. Not read by the write-up. | `vtxLxyDr_signal_2018/dask_check` (1 sample, moved by hand) |

### Condor production of `vtxLxyDr_signal_2018`

The 180 samples of `vtxLxyDr_signal_2018` were produced with the condor chunk pipeline
(`condor/README.md`), not by the notebook: `condor/submit_vtxlxydr_signal.sub` (same channels,
hist collections and selections file as the notebook; it passes `none` for the standard
corrections, which the production predates) over `condor/job_args_vtxlxydr.txt`,
964 jobs of 10 files. The per-chunk filelists are not committed; from the repository root,
with the `sidm_venv` active, they and the job list are regenerated by

```bash
for fs in 4mu 2mu2e; do
  python condor/make_job_args.py --samples-file condor/samples_vtxlxydr_${fs}.txt \
    --dataset llpNanoAOD_v2 --location-cfg signal_${fs}_v10.yaml --files-per-job 10 \
    --replace-xcache --census-skip sidm/studies/vertexChi2/signal_skip.json \
    --outdir condor/filelists_vtxlxydr --job-args /tmp/job_args_vtxlxydr_${fs}.txt
done
cat /tmp/job_args_vtxlxydr_4mu.txt /tmp/job_args_vtxlxydr_2mu2e.txt > condor/job_args_vtxlxydr.txt
```

The census skip list `signal_skip.json` (next to this README; the Dask notebook reads the same
file) drops 17 known-corrupt files from 15 samples. One of them, `part-1.root`
of `4Mu_200GeV_5p0GeV_200p0mm` (1718 events), had been read by the earlier `vtxChi2_signal_2018`
run, so that sample holds 92544 events here against 94262 there; `vtx_outcome_lxy_dr/closure_checks.txt`
lists it as the only sample whose event count differs between the two runs. `run_job.sh` books
**unweighted** histograms and processes chunks of 50000 events, so the merged histograms hold raw
counts, not the lumi x xs weighted sums of the Dask outputs; a per-LJ fraction is the same
either way when every event of a sample carries the same generator weight, and the yield of a
sample is its count times `lumi*xs/sumw` from the merged metadata. The chunks are merged with

```bash
python sidm/scripts/merge_coffea_chunks_eos.py \
  --input-eos-dir /store/user/$USER/sidm_condor/vtxLxyDr_signal_2018 \
  --output-eos-dir /store/group/lpcmetx/SIDM/coffea_outputs/$USER/vtxLxyDr_signal_2018 \
  --workdir <scratch>/merge_vtxlxydr --keep-local --filelists-dir condor/filelists_vtxlxydr \
  --selections base_ljDisplacementIso_cosAlpha,base_ljDisplacementIso_cosAlpha_ljSpread \
  --hist-collections mu_lj_vtx_outcome_lxy_dr,mu_lj_vtx_chi2_dsa_pt,dsa_tail_diag \
  --selections-cfg configs/selections_vertexChi2_as_run.yaml --chunksize 50000 --unweighted-hist
```

where `<scratch>` is any local scratch directory (`--keep-local` leaves the chunk copies there),
and the merged files are copied to `RunOutputFiles/vtxLxyDr_signal_2018/`.

### Condor production of `vtxGenDr_signal_2018`

`vtxGenDr_signal_2018` repeats that production with one more histogram: the vertex outcome
against the generated dR of the matched dark photon's two daughter muons. The hist collection
`mu_lj_vtx_outcome_gen_dr` holds it (`mu_lj_vtxChi2_genA_lxy_genMuMu_dR`: `vtx_chi2` x coarse
generated Lxy x generated dR) together with `mu_lj_vtxChi2_genA_lxy_mumu_dR` and
`mu_lj_vtx_chi2` for the closures. Samples, channels, selections file, census skip list and job
layout are those of `vtxLxyDr_signal_2018`: `condor/submit_vtxgendr_signal.sub` over
`condor/job_args_vtxgendr.txt`, 964 jobs of 10 files, unweighted, no standard corrections. The
job list and the filelists are regenerated with the recipe above, with `vtxgendr` for `vtxlxydr`
in `--outdir` and in the job-list names (the two sample lists `condor/samples_vtxlxydr_*.txt`
serve both productions). The code was `main` at 9d9b387c plus the two changes that book the
histogram and add the submission. The provenance sidecars record commit 78b7404b, which is not
in this repository: it is 9d9b387c with those same two changes. The chunks are merged with

```bash
python sidm/scripts/merge_coffea_chunks_eos.py \
  --input-eos-dir /store/user/$USER/sidm_condor/vtxGenDr_signal_2018 \
  --output-eos-dir /store/group/lpcmetx/SIDM/coffea_outputs/$USER/vtxGenDr_signal_2018 \
  --workdir <scratch>/merge_vtxgendr --keep-local --filelists-dir condor/filelists_vtxgendr \
  --selections base_ljDisplacementIso_cosAlpha,base_ljDisplacementIso_cosAlpha_ljSpread \
  --hist-collections mu_lj_vtx_outcome_gen_dr \
  --selections-cfg configs/selections_vertexChi2_as_run.yaml --chunksize 50000 --unweighted-hist
```

### Vertex outcome vs Lxy and dR, and the halo fractions (`vtx_outcome_lxy_dr/`)

The tables in `vtx_outcome_lxy_dr/` are read from `vtxLxyDr_signal_2018` (signal), the existing
`vtxChi2_signal_2018` (the reference of the efficiency closure), and `vtxResp2_2018Full` / `vtxRespMC2_2018` (data and
MC halo fractions). The scripts sit next to them; from the repository root, with the
`sidm_venv` active and `PYTHONPATH` at the repository root, `OUT` an output directory:

```bash
D=sidm/studies/vertexChi2/vtx_outcome_lxy_dr
python $D/vtx_outcome_vs_lxy_dr.py RunOutputFiles/vtxLxyDr_signal_2018 $OUT   # first: halo_signal reads its by-sample table
python $D/halo_signal.py RunOutputFiles/vtxLxyDr_signal_2018 $OUT
python $D/halo_data_mc.py $OUT
python $D/condor_vs_dask.py $OUT/condor_vs_dask.txt
```

They also write PNG figures, which are not committed. `an_figures.py OUTDIR` draws the
figures for the analysis note from `outcome_vs_lxy.csv`, `outcome_vs_dR.csv` and
`outcome_vs_gen_dR.csv` (next section), reference channel, mean over the samples: the outcome
fractions against the generated Lxy, against the reconstructed dR and against the generated dR,
in the style of `truth_kinematics_forAN/an_style.py`, as vector PDF with
PNG previews (the "bound-state mass" of their legends is the mediator mass of this README). Outcome of a mu-type LJ: **no vertex**
(`vtx_chi2 < 0`: the -1 sentinel of an LJ without a stored within-LJ vertex, and the three LJs
in the grid whose stored normChi2 lies in [-0.5, 0), which the vertex cut and `mu_lj_vtx_chi2`
count as no vertex too), **chi2 >= 5** (a stored vertex with normChi2 >= 5), **pass**
(`0 <= vtx_chi2 < 5`); the underflow bin below -1.5 is reported as `frac_other_vtx_bins` and is
empty. Every LJ of the two channels has >= 2 muons (`Mu >= 2`). The
**robust** samples (162 of 180) have `N_eff = (sum w)^2 / sum w^2 >= 200`, computed from the
reference channel's `mu_lj_dsaMuon_pt_relPtErr`, flow included. Two ways of combining samples:
**mean** / **sem** is the unweighted mean of the per-sample fractions, and its standard error,
over the samples with `N_eff >= 20` in that bin; **pooled** sums the LJs of all samples, each
sample weighted to its yield at 1 fb (`lumi*xs*eff/sumw`, since the condor histograms are raw
counts; `eff` is the efficiency of the generator-level filter the samples were made with, 0.23
to 0.90 across the v10 samples, PR #81). The mean and sem do not depend on that weight.

| file | what it holds | written by |
|---|---|---|
| `outcome_by_sample.csv` | Per sample and channel: mass point, `N_eff` from this run and from `vtxChi2_signal_2018`, robust flag, LJ count and yield, 1-muon and >= 2-muon fractions, the three outcome fractions among >= 2-muon LJs, the gen-matched fraction, pass/all, the old run's per-LJ vertex-veto efficiency (reference channel), and the input-file and event counts of both runs. | `vtx_outcome_vs_lxy_dr.py` |
| `outcome_vs_lxy.csv` | Robust samples, gen-matched LJs (nearest `genAs_toMu` within dR < 0.4): outcome fractions in bins of the dark photon's Lxy (cm), mean / sem and pooled with its error, for all samples and per mediator mass, in both channels. | `vtx_outcome_vs_lxy_dr.py` |
| `outcome_vs_dR.csv` | The same in bins of dR between the LJ's two highest-pT muons (PF and DSA pooled). At mediator masses of 500 and 1000 GeV a reconstructed dR above 0.1 is mostly not the opening angle of the decay pair: one of the two muons lies away from both generated daughters, while the generated pair is narrow (`dr_truth_match/` below). | `vtx_outcome_vs_lxy_dr.py` |
| `outcome_vs_lxy_dR.csv` | The same in coarse Lxy x dR bins. | `vtx_outcome_vs_lxy_dr.py` |
| `closure_checks.txt`, `closure_vtx_chi2_marginal.csv` | Closure checks. The marginal closure (the CSV, summarised in the text file): in every sample, channel and new histogram, the `vtx_chi2` marginal, grouped as no vertex / pass / chi2 >= 5, against `mu_lj_vtx_chi2`. The efficiency closure: per-sample pass/all against the old run's per-LJ efficiency, the samples whose input-file or event count differs between the two runs, and the robust means by mediator mass, m_Zd and ctau. | `vtx_outcome_vs_lxy_dr.py` |
| `condor_vs_dask.txt` | `4Mu_100GeV_0p25GeV_0p02mm` from condor against the Dask run in `dask_check/`, every bin of the three new histograms and `mu_lj_vtx_chi2`, both channels. | `condor_vs_dask.py` |
| `halo_fractions_data_mc.csv` | 2018 data (68 skims) and the eight-sample MC, `cr_1egmOr1mu_lj` and its Spread-veto version, all DSA muons of mu-type LJs. Per DSA pT bin and nHits bin: the fraction with `|dxy_BS| >= 10` cm and the fraction with zero DT hits. Per pT bin and \|eta\| bin: the horizontal-over-vertical phi ratio. Data errors binomial, MC errors from sum w^2. | `halo_data_mc.py` |
| `halo_fractions_signal.csv` | The same three quantities for the robust signal samples in both signal channels, all samples and per mediator mass, pooled and mean / sem. | `halo_signal.py` |

In the halo tables, horizontal and vertical are as in sections 7 and 10.4 of the write-up: a phi
bin is vertical if its centre is within 0.5 rad of +-pi/2, horizontal if within 0.5 rad of 0
or +-pi. `halo_common.py` holds these index groups and the fraction and ratio helpers shared by
`halo_signal.py` and `halo_data_mc.py`. A horizontal-over-vertical ratio with no horizontal
entries is 0 with error 0; one with no vertical entries is left blank.

### Vertex outcome against the generated dR, and tests of the trends (`vtx_outcome_lxy_dr/`)

At mediator masses of 500 and 1000 GeV a reconstructed dR above 0.1 between the two leading muons
is mostly not the opening angle of the decay pair (`dr_truth_match/` below), so the outcome is
also tabulated against the generated dR of the matched dark photon's two daughter muons, from
`vtxGenDr_signal_2018`. With `GENDR` and `LXYDR` local copies of the merged outputs of
`vtxGenDr_signal_2018` and `vtxLxyDr_signal_2018`:

```bash
D=sidm/studies/vertexChi2/vtx_outcome_lxy_dr
python $D/vtx_outcome_vs_gen_dr.py $GENDR $LXYDR $OUT     # reads outcome_by_sample.csv next to it
python $D/physical_sense_tests.py $OUT $OUT/physical_sense
python $D/trigger_vs_lxy.py $OUT/trigger_vs_lxy.csv condor/filelists_vtxgendr/<sample>_0.txt ...
```

Both condor productions are unweighted, so every number of these tables is a raw LJ count or a
ratio of raw counts. A mass point is one pair (mediator mass, m_Zd), its lifetimes and both final
states summed unless a table splits them. The generated dR is that of the dark photon the LJ is
matched to (nearest `genAs_toMu` within dR < 0.4), between its own two daughters, from their
momenta as generated.

| file | what it holds | written by |
|---|---|---|
| `counts_gen_dr.csv.gz` | Per sample and channel: LJs by outcome in bins of the coarse generated Lxy ([0, 1), [1, 10), [10, 100), [100, 1000) cm) and of the generated dR (edges 0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5). Underflow and overflow are bins of their own; "< 0" is an LJ without a matched dark photon. | `vtx_outcome_vs_gen_dr.py` |
| `counts_lxy.csv`, `counts_lxy_muonN.csv` | From `vtxLxyDr_signal_2018`: LJs by outcome per sample, channel and bin of the fine generated Lxy (the ten bins of `outcome_vs_lxy.csv`); the second file per muon multiplicity, for [100, 300) and [300, 1000) cm only. | `vtx_outcome_vs_gen_dr.py` |
| `outcome_vs_gen_dR.csv` | `outcome_vs_dR.csv` with the generated dR in place of the reconstructed one: robust samples, gen-matched LJs of every Lxy, mean / sem and pooled, for all samples and per mediator mass, in both channels. | `vtx_outcome_vs_gen_dr.py` |
| `closure_gen_dr_marginal.csv`, `gen_dr_checks.txt` | Per sample and channel, the new histogram summed over both generated axes against `mu_lj_vtx_chi2`, and summed over the generated dR against the reconstructed-dR histogram of the same run summed over the reconstructed dR. Both are exact in all 360 (sample, channel) pairs. No gen-matched LJ has fewer than two daughters. | `vtx_outcome_vs_gen_dr.py` |
| `physical_sense/TC_floor.csv`, `TC_dr_mix.csv` | Decays within 10 cm: the no-vertex fraction against the generated dR, per mediator mass; and per mediator mass the fraction over all dR, the share of the LJs below 0.002, the fraction that its dR distribution gives with the fraction of each dR bin taken from all masses summed, and the measured fraction minus that, with the parts of the difference from the dR bin below 0.002 and from the bins between 0.002 and 0.01. | `physical_sense_tests.py` |
| `physical_sense/TA_points.csv`, `TA_cells.csv`, `TA_order.csv`, `TA_final_state.csv` | Whether the mass points share one fraction inside a (coarse Lxy, generated dR) cell. Per cell and mass point (points with >= 200 LJs): the no-vertex and pass fractions; per cell their range over the mass points and the chi2 of one common fraction; the rank correlation of the fractions with the mediator mass at fixed m_Zd and with m_Zd at fixed mediator mass; 4Mu against 2Mu2E. The cells are coarse (one spans 10 to 100 cm), and a different Lxy or dR distribution inside a cell also makes the points disagree; tables 21 and 22 of `dr_truth_match/` repeat the test in finer cells. | `physical_sense_tests.py` |
| `physical_sense/TB_samples.csv`, `TB_pooled.csv`, `TB_pooled_by_M.csv` | Lxy [100, 300) against [300, 1000) cm: the three outcome fractions in each bin and their difference, per sample with >= 50 LJs in both bins, and with the counts of those samples summed per mass point and per mediator mass. | `physical_sense_tests.py` |
| `physical_sense/TB_final_state.csv` | The LJs and the no-vertex fraction of those two Lxy bins per final state and m_Zd, all samples summed. | `physical_sense_tests.py` |
| `physical_sense/TD_widest.csv` | Generated dR [0.2, 0.5) against [0.1, 0.2) at a fixed mass point and coarse Lxy bin: the no-vertex fraction in each and their difference, for both final states together and for each (rows with >= 100 LJs in both dR bins). | `physical_sense_tests.py` |
| `trigger_vs_lxy.csv` | Twelve samples (both final states at 150 and 500 GeV, m_Zd of 0.25, 1.2 and 5 GeV, the longest lifetime of each; the first two filelists, 20 files, 537,470 events): events and the fraction that pass the trigger (the OR of the "pass triggers" cut, no other selection) by the decay radius of the dark photons that decay to muons, in 4Mu the larger and the smaller of the two. | `trigger_vs_lxy.py` |

Numbers, reference channel. The first, second and fourth items are from `vtxGenDr_signal_2018`
(5,761,625 gen-matched LJs in the 180 samples), the third from `vtxLxyDr_signal_2018`
(5,761,693).

- **Collimated pairs.** Among the 1,288,870 LJs whose dark photon decays within 10 cm, the
  no-vertex fraction is 0.18 to 0.27 for a generated dR below 0.002 at every mediator mass (239
  of 1,059 at 100 GeV, 23,301 of 85,269 at 1000 GeV), 0.083 to 0.153 in [0.002, 0.005), 0.055
  to 0.140 in [0.005, 0.01), 0.025 to 0.089 in [0.01, 0.02) and 0.004 to 0.011 in [0.1, 0.2)
  (0.009 to 0.014 in [0.2, 0.5) at 100 to 200 GeV, 1 of 96 at 500 GeV). Over all dR inside
  10 cm it is 0.050, 0.068, 0.075, 0.122, 0.157 and 0.139 at 100, 150, 200, 500, 800 and
  1000 GeV, and the share of the LJs below 0.002 is 0.03, 0.06, 0.10, 0.39, 0.43 and 0.30
  (`TC_dr_mix.csv`). With the fraction of each dR bin taken from all masses summed, the dR
  distribution of each mass gives 0.074, 0.080, 0.081, 0.130, 0.147 and 0.125: the step from
  100 to 200 GeV up to 500 to 1000 GeV, the maximum at 800 GeV, and 0.073 of the range of
  0.108. The measured fraction is below that at 100 to 500 GeV (by 0.025, 0.012, 0.006 and
  0.008) and above it at 800 and 1000 GeV (by 0.010 and 0.013): the fraction inside a dR bin
  changes with the mediator mass as well. At 100 and 150 GeV the difference comes from the
  bins between 0.002 and 0.01 (0.023 and 0.009 of it), where the fraction rises with the
  mediator mass (0.083 at 100 GeV and 0.145 at 1000 GeV in [0.002, 0.005), 0.055 and 0.140 in
  [0.005, 0.01)). From 200 GeV up it comes mostly from the bin below 0.002 itself (-0.006,
  -0.010, +0.009 and +0.009 at 200, 500, 800 and 1000 GeV), where the fraction is 0.183 at
  200 GeV and 0.273 at 1000 GeV. A bin holds a different m_Zd at each mediator mass (next
  item).
- **Mass points inside a cell.** Of the 32 (coarse Lxy, generated dR) cells, 26 reject one
  common no-vertex fraction for all mass points at p < 0.001 (27 of 32 for the pass
  fraction). In Lxy [10, 100) cm and dR [0.02, 0.05), the fraction runs from 0.071
  (500 GeV, m_Zd = 5 GeV) to 0.150 (200 GeV, 1.2 GeV) over nine mass points and 891,217 LJs.
  Within a cell the fraction is lower at the higher m_Zd in 26 of the 30 pairs of mass points
  of equal mediator mass, and its rank correlation with the mediator mass at fixed m_Zd is
  positive in 32 of 44 groups of three to six masses. These cells are coarse. In the sixteen
  2Mu2E samples of `dr_truth_match/` the disagreement does not go away when each Lxy slice is
  cut in up to eight: the spread of the fraction between mass points beyond the binomial
  expectation stays at 0.03. It falls to between 0.013 and 0.021 with the muon content of the
  pair as a further cell variable, or with the dR bin below 0.002 cut in three, and what is
  left of it is in the DSA-DSA pairs (tables 21 and 22 there).
- **Largest radii.** The no-vertex fraction falls from [100, 300) to [300, 1000) cm in 37 of the
  43 samples with at least 50 LJs in both bins, in 20 of them by more than twice its error (it
  rises by more than that in one), so the fall of `outcome_vs_lxy.csv` is not made by the change
  of the samples that enter the bin (the fraction without a passing vertex falls in 30 of the
  43). With the counts summed per mass point the change is -0.18 to -0.28 at m_Zd = 0.25 GeV
  from 150 GeV up (-0.05 +- 0.03 at 100 GeV), -0.14 to -0.23 at 1.2 GeV from 500 GeV up
  (+0.06 +- 0.03, -0.03 +- 0.02 and -0.04 +- 0.02 at 100, 150 and 200 GeV), and at 5 GeV 0.00
  at 100 to 200 GeV and -0.04 to -0.06 from 500 GeV up. 35 of the 43 samples are 4Mu: beyond
  300 cm, 3,785 of the 3,833 LJs at 0.25 and 1.2 GeV are 4Mu (`TB_final_state.csv`). At 5 GeV
  both final states populate the bin and differ there, 0.059 (96 of 1,621) in 2Mu2E against
  0.108 (1,403 of 12,944) in 4Mu, after 0.150 and 0.153 in [100, 300) cm. The trigger selects
  differently in the two (`trigger_vs_lxy.csv`). A 2Mu2E event whose muon pair decays beyond
  300 cm passes it in 0.001 to 0.004 of the cases at 0.25 and 1.2 GeV (10 of 11,333 to 39 of
  8,938) and in 0.014 to 0.044 at 5 GeV, against 0.13 to 0.23 in [100, 300) cm. A 4Mu event
  passes in 0.20 to 0.34 with both pairs beyond 300 cm and in 0.52 to 0.75 with the other pair
  within 100 cm. The 2Mu2E LJs beyond 300 cm are therefore a small, trigger-selected part of
  the decays there. Whether that selection makes their no-vertex fraction lower, and what
  makes the fraction fall beyond 300 cm in 4Mu, is not established.
- **The widest generated dR.** In `outcome_vs_gen_dR.csv` the no-vertex fraction (mean over the
  samples) falls from 0.315 +- 0.010 below 0.002 to 0.046 +- 0.004 in [0.1, 0.2) and is
  0.124 +- 0.017 in [0.2, 0.5). The 82,405 LJs of the last bin are from m_Zd = 5 GeV at 100 to
  200 GeV but for 851. At a fixed mass point (`TD_widest.csv`) the last bin is above
  [0.1, 0.2) beyond 100 cm: 0.164 (1168 of 7130) against 0.076 at 100 GeV, 0.192 against 0.078
  at 150 GeV, 0.208 (383 of 1841) against 0.083 at 200 GeV. Within 100 cm the two bins agree to
  0.018 with both final states summed. The final states differ. In 2Mu2E the last bin is high
  in [10, 100) cm as well: 0.194 (75 of 386) against 0.055 at 100 GeV, 0.230 (38 of 165)
  against 0.037 at 150 GeV (0.191 and 0.305 beyond 100 cm). In 4Mu, which
  holds 13,109 of the 13,752 LJs of the last bin beyond 100 cm, it is high only there. These
  tables do not establish the cause; a truth match of these pairs is in
  `dr_truth_match/wide_pair/` (below).

### Truth match of the two leading muons of a mu-type LJ (`dr_truth_match/`)

Against the dR of the LJ's two leading muons, the no-vertex fraction at mediator masses 500 and
1000 GeV is highest in the widest bins, [0.2, 0.3) and [0.3, 0.5) (`outcome_vs_dR.csv` above;
table 13 here for the samples of this study). This study asks whether that dR measures the
opening angle of the Zd decay pair. It writes one record per mu-type LJ, with its PF and DSA
constituents and the generated Zd daughter muons of the event, and matches the two leading
muons to the daughters in direction. Signal MC only: 16 2Mu2E samples of the 2018 v10 grid,
listed in `condor/samples.txt` (mediator mass 100, 500 and 1000 GeV; m_Zd of 0.25, 1.2 and
5 GeV, without the 0.25 GeV point at 100 GeV; two lifetimes each), reference channel
`base_ljDisplacementIso_cosAlpha` of `configs/selections_vertexChi2_as_run.yaml` (no cosmic
veto), no corrections.

Terms used below and in `tables.md`:

- **m0, m1**: the LJ's two highest-pT muons, PF and DSA pooled (m2 is the third). **reco dR** =
  dR(m0, m1), the quantity of `mu_lj_mumu_dr` and of `outcome_vs_dR.csv`.
- **d1, d2**: the two generated muons from the decay of the dark photon matched to the LJ (the
  generated muons whose first ancestor with a different pdgId is that Zd). In 2Mu2E they are the
  event's only Zd daughter muons (`tables.md`, line above table 7); a 4Mu event has four, two per
  dark photon (`wide_pair/` below). **gen dR** = dR(d1, d2) at the decay vertex. **Lxy**:
  transverse decay length (cm) of the generated dark photon matched to the LJ. **Gen-matched**:
  a generated dark photon decaying to muons (the `genAs_toMu` collection) within dR < 0.4 of the
  LJ, as in `outcome_vs_lxy.csv`.
- **Propagated match**: a reco muon is compared with the direction a daughter has at the point of
  closest approach of its helix to the beamline, where a track quotes its momentum (helix in a
  uniform 3.8 T field along z, from the generated decay vertex; no material, no field map). The
  committed tables use it (`--propagate`). The gen dR stays the opening angle at the decay.
- **Outcome**: **no vertex** (`vtx_chi2 < 0`), **chi2 >= 5**, **pass** (`0 <= vtx_chi2 < 5`), as in
  `vtx_outcome_lxy_dr/`; **not pass** = no vertex or chi2 >= 5.
- **Match radius**: the 99% quantile of dR(reco muon, nearest daughter) in gen-matched LJs with
  exactly two opposite-charge muons and a passing vertex, all three masses pooled: 0.0556 for PF
  and 0.0761 for DSA muons (table 3). A muon is **within its radius** of a daughter when their dR
  is below it.
- **Classes** of a gen-matched LJ, from m0 and m1 (table 4): **D**, m0 or m1 is outside its radius
  of both daughters; **A**, m0 and m1 are within their radii of different daughters, have
  opposite charge, and the reco dR agrees with the gen dR (ratio within [1/2, 2], or a difference
  of at most the sum of the two radii), i.e. the decay pair; **B**, as A but with the reco dR wider
  than that tolerance; **C**, both nearest to the same daughter, with the same charge or a PF-DSA
  pair sharing muon-system segments, i.e. two reconstructions of one daughter; **E**, the rest.

From the repository root, with the `sidm_venv` active, the records are produced by 32 condor
jobs of up to 25 files (`condor/job_args.txt`). Build `condor/sidm_code.tar.gz` as in
`condor/README.md` (section 8), then

```bash
D=sidm/studies/vertexChi2/dr_truth_match
python condor/make_job_args.py --samples-file $D/condor/samples.txt --dataset llpNanoAOD_v2 \
  --location-cfg signal_2mu2e_v10.yaml --files-per-job 25 --replace-xcache \
  --census-skip sidm/studies/vertexChi2/signal_skip.json \
  --outdir $D/condor/filelists --job-args /tmp/job_args_dr_truth.txt
cd $D/condor && mkdir -p logs out && condor_submit submit.sub && cd -
```

The filelists are not committed; this command writes the same lists the committed tables were
made from. The job list it writes is discarded: its paths are relative to `condor/`, and
`condor/job_args.txt` here is the same list relative to `dr_truth_match/condor/`. `run.sh`
builds the worker environment as `condor/run_job.sh` does, from the repository's
`condor/requirements.txt` and `condor/constraints.txt`, and runs `extract.py` on one filelist;
each job writes one parquet file to `out/`. `extract.py` also runs by hand on a single filelist
(`python $D/extract.py SAMPLE FILELIST OUT.parquet`, from the repository root). The tables are
then made, in under a minute, by

```bash
python $D/analyze.py $D/condor/out $OUT --propagate
```

which also writes `lj_truth.parquet` (the merged records) and `lj_classified.parquet` (the
gen-matched LJs with every derived column) to `$OUT`; neither is committed.

| file | what it holds |
|---|---|
| `tables.md` | Every table below with its caption, the input counts and the check against the histograms of `vtxLxyDr_signal_2018`. |
| `check_vs_outcome_by_sample.csv` | Per sample: the LJ count and the three outcome fractions from the records, against `vtx_outcome_lxy_dr/outcome_by_sample.csv` (reference channel), and the largest difference. |
| `table1.csv` | Per mass and reco-dR bin: 10/50/90% quantiles of the gen dR, the fraction with gen dR < 0.05, and the fraction whose gen dR falls in the same bin. |
| `table2.csv` | Per mass and reco-dR bin: PF+PF / PF+DSA / DSA+DSA leading pair, 2 and >= 3 muons, same-charge leading pair. |
| `table3.csv` | dR(reco muon, nearest daughter) quantiles for PF and DSA muons (pooled, by DSA Lxy range, per mass): the match radii. |
| `table4.csv` | Per mass and reco-dR bin: class fractions A to E with binomial errors, the ambiguous fraction, and the A-or-B fraction with reco dR > 2 x gen dR. |
| `table5.csv` | Per mass, reco-dR bin and class: outcome fractions and counts. |
| `table6.csv` | Classes C and D at reco dR >= 0.1, 500 and 1000 GeV: the duplicate or extra muon (type, pT, \|eta\|, dR to its nearest daughter), the third muon, and how many daughters the LJ's muons cover. |
| `table7.csv` | 500 and 1000 GeV pooled, reco dR >= 0.1 and >= 0.2: LJ count and share, class fractions, no-vertex fraction and share of all no-vertex LJs, and the class D split by which leading muon is outside. |
| `table8.csv` | The class D LJs of table 7 at reco dR >= 0.1 with exactly one leading muon outside: what that muon is, by type (azimuth within 0.03 of a daughter, charge and pT match to it, \|eta\| difference, and the outcome and Lxy of those that match the other daughter). |
| `table9.csv` | \|eta(reco) - eta(nearest daughter)\| quantiles of the two muons of two-muon class A LJs, 500 and 1000 GeV, by type and Lxy. |
| `table10.csv` | Leading muons of all gen-matched LJs, 500 and 1000 GeV, within 0.03 in azimuth of a daughter: fraction with \|eta\| difference > 0.05, by type and Lxy. |
| `table11.csv` | Outcome fractions against the gen dR, per mass (500, 1000). |
| `table12.csv` | No-vertex and not-pass fractions against the gen dR per Lxy slice, 500 and 1000 GeV pooled. |
| `table13.csv` | The same against the reco dR. |
| `table14.csv` | Per mass and reco-dR bin, with per-mass match radii: decay-pair fraction and fraction with a leading muon outside, next to the class A and D fractions of table 4, and the median gen dR. |

The records reproduce the histograms of `vtxLxyDr_signal_2018`: the LJ count equals that of
`outcome_by_sample.csv` in 15 samples and is one LJ short in `2Mu2E_500GeV_0p25GeV_4p0mm`
(17911 against 17912; that run predates the `dR_outer` correction, and against
`vtxGenDr_signal_2018`, made with it as the records are, the gen-matched LJ count is equal in
all 16 samples), the three outcome fractions agree to 3e-5 in that sample and
to 1e-6 in the others (`check_vs_outcome_by_sample.csv`), and the reco dR recomputed from the
records equals `mu_lj_mumu_dr` to 4.17e-07 (`tables.md`). The records hold 85438 LJs, 85355 of
them gen-matched (`tables.md`).

**Result.** At 500 and 1000 GeV a reconstructed dR above 0.1 is mostly not the opening angle
of the decay pair: one of the two leading muons lies outside its radius of both daughters, while
the generated pair is narrow. Against the generated dR the no-vertex fraction falls.

- Such LJs are few: of the 84383 gen-matched LJs at 500 and 1000 GeV, 3016 (3.6%) have reco
  dR >= 0.1 and 555 (0.7%) have reco dR >= 0.2. They hold 0.0535 of the no-vertex LJs of
  those masses (table 7).
- Most of them are class D: 0.647 at reco dR >= 0.1 and 0.941 at >= 0.2 (table 7); per bin
  0.519, 0.897 and 1.000 in [0.1, 0.2), [0.2, 0.3) and [0.3, 0.5) at 500 GeV, and 0.806, 0.982
  and 1.000 at 1000 GeV (table 4). The generated pairs in these bins are narrow: the median gen
  dR of all LJs of the bin is 0.0447, 0.0110 and 0.0105 at 500 GeV, 0.0168, 0.0113 and 0.0081 at
  1000 GeV (table 1). The rest at reco dR >= 0.1 is the decay pair (class A, 0.337), nearly all
  of it in [0.1, 0.2).
- What the outside muon is (table 8; the 1901 class D LJs at reco dR >= 0.1 with exactly one
  leading muon outside, 1129 PF and 772 DSA). A PF outside muon sits at the azimuth of the
  pair with a shifted eta: 0.959 are within 0.03 in azimuth of a daughter (0.041 of neither),
  while by the class definition they are more than 0.0556 in dR from both. Which daughter it
  is rests on the charge alone, since 0.836 are within 0.03 of both daughters: 0.663 have the
  charge of the daughter not claimed by the other leading muon (0.541 also within 30% of its
  pT; the pT window has no control sample), with an \|eta\| difference to it of 0.141
  (median; 0.089 and 0.257 at 10% and 90%) and Lxy of 19.6 to 58.9 cm (10% to 90%, median
  43.7), and 0.283 have the charge of the claimed daughter. Of the 772 DSA outside muons,
  0.437 match the other daughter in azimuth and charge (\|eta\| difference 0.175, median),
  0.088 the claimed one, and 0.475 neither; 0.935 of the first group have no vertex. Taken
  together, the outside muon has the azimuth and charge of the other daughter in 1085 of the
  1901 LJs (57%), which is 36% of the 3016 LJs at reco dR >= 0.1.
- The PF eta measurement degrades with Lxy in the well-matched LJs too: the 90% quantile of
  \|eta(reco) - eta(daughter)\| in class A grows from 0.0007 at Lxy <= 1 cm to 0.0054 at 20 to
  30 cm and 0.018 to 0.025 from 30 to 100 cm (table 9). Class A requires dR < 0.0556, so table
  9 cannot hold offsets of the size seen in table 8; it shows the core widening, not the tail.
  The fraction of azimuth-matched PF muons with an \|eta\| difference above 0.05 grows from
  0.0028 to 0.0136 and then 0.060 to 0.109 over the same ranges (table 10; this population
  contains the class D muons themselves, so it restates table 8 against Lxy and is not an
  independent check). The DSA eta width does not depend on Lxy: the 90% quantile stays at
  0.0135 to 0.0155 below 300 cm (table 9).
- Against the gen dR the no-vertex fraction falls (not monotonically in every bin): at 500 GeV
  from 0.359 at gen dR < 0.002 to 0.073 in [0.05, 0.1) and 0.057 in [0.1, 0.2) (402 LJs); at
  1000 GeV from 0.304 to 0.071 and 0.048 (21 LJs) (table 11). The highest fractions are those of
  the most collimated pairs. In every Lxy slice below 300 cm it is lower in [0.05, 0.1) than
  below 0.002: 0.062 against 0.424 at 10 to 100 cm, 0.117 against 0.480 at 100 to 300 cm
  (table 12). Against the reco dR, in the same slices, it is 0.462 and 0.574 in [0.2, 0.3) and
  [0.3, 0.5) at 10 to 100 cm, and 0.646 and 0.905 (21 LJs) at 100 to 300 cm (table 13).
- With per-mass match radii (500 GeV: 0.0614 PF, 0.0667 DSA; 1000 GeV: 0.0409, 0.0603) the
  fraction with a leading muon outside both daughters is 0.532 and 0.848 in [0.1, 0.2) at 500
  and 1000 GeV, against class D fractions of 0.519 and 0.806 with the pooled radii, and 0.897 and
  0.991 in [0.2, 0.3) against 0.897 and 0.982 (table 14).
- At 100 GeV the generated pairs are themselves wide: the median gen dR is 0.162 in [0.1, 0.2)
  and 0.216 in [0.2, 0.3), and 0.49 and 0.58 of the LJs have their gen dR in the same bin as
  their reco dR (table 1; 321 and 139 LJs). There the reco dR does follow the opening angle for
  about half of the LJs. The class fractions at this mass depend on the match radius: with the
  pooled radii class A is 0.586 and 0.317 and class D 0.396 and 0.647 (table 4); with the
  100 GeV radii, 0.120 for PF and 0.183 for DSA, which come from 269 and 857 muons and are as
  large as the gen dR itself, the decay-pair fraction is 0.766 and 0.597 (table 14).

Caveats:

- Sixteen 2Mu2E samples, no 4Mu; two dark-photon masses per mediator mass at 100 GeV, three at
  500 and 1000 GeV (`condor/samples.txt`), reference channel only. In these samples one Zd
  decays to muons, so d1 and d2 are the event's only Zd daughter muons; for 4Mu `extract.py`
  assigns the daughters per LJ, which the tables 1 to 22 do not use.
- The propagation places the decay vertex at Lxy along the dark photon's momentum from the
  origin, so the beam-spot offset is ignored.
- LJs are pooled by raw count, without cross-section weights, so each sample counts by its number
  of LJs, from 18 to 17911 (`check_vs_outcome_by_sample.csv`).
- The match uses direction only. A muon with the right azimuth and a wrong eta is outside its
  radius; table 8 is what ties such muons to the decay pair, with a fixed 0.03 azimuth window
  and a 30% pT window. It does not show whether the outside muon is a daughter reconstructed
  with a wrong eta or another track at the same azimuth; the 0.475 of DSA and 0.055 of PF
  outside muons that match neither daughter are not identified. The radii come from LJs with a
  passing vertex, so they describe well-measured muons.
- Tables 7 to 14 give fractions with their counts and no errors; several bins hold a handful of
  LJs (gen dR >= 0.2 holds 1 LJ per mass, table 11). Tables 12 and 13 pool 500 and 1000 GeV, and
  the gen-dR dependence in them is not separated from the dark-photon mass and pT that set the
  gen dR.
- Class B holds 24 LJs over the three masses (table 5). With daughter directions at the decay
  vertex it is empty by construction; the propagated match compares with directions at the
  beamline, while the gen dR is the opening angle at the decay.

#### Muon content of the pair (`pair_content.py`, tables 15 to 22)

`pair_content.py` reads the `lj_classified.parquet` that `analyze.py` writes and tabulates the
vertex outcome by the muon content of the LJ:

```bash
python $D/pair_content.py $OUT/lj_classified.parquet $OUT     # tables 15 to 22, pair_content.md
```

**Content** is defined for an LJ with exactly two muons, from its PF and DSA muon counts: PF-PF
(35,602 LJs), PF-DSA (18,319), DSA-DSA (26,768); the other 4,666 of the 85,355 gen-matched LJs
have three or more muons and enter only the shares of table 17. **Charge** is that of the two
muons, opposite or same: table 15 splits by it, tables 19 and 20 take the opposite-charge pairs,
and tables 16, 18, 21 and 22 take both together. The **region** of tables 20 and 22 is from the
mean eta of d1 and d2: barrel below 0.9, endcap above 1.2 in absolute value. LJs are raw counts,
the three mediator masses (eight mass points) pooled unless a table is per mass.

| table | what it holds |
|---|---|
| 15 | No-vertex fraction by content and charge. |
| 16 | No-vertex fraction by Lxy slice and content. |
| 17 | Share of each content, and of the LJs with three or more muons, per mediator mass and Lxy slice. |
| 18 | No-vertex fraction by gen dR and content. |
| 19 | DSA-DSA LJs of opposite charge: the three outcomes per Lxy slice. |
| 20 | DSA-DSA LJs of opposite charge: the no-vertex fraction per Lxy slice and region. |
| 21 | Whether the mass points share one no-vertex fraction inside cells of (Lxy, gen dR), and with the content, the lower generated pT of d1 and d2, or both as further cell variables: chi2 and ndf summed over the cells, the LJs in the test, and the spread of the fraction between mass points beyond the binomial expectation (`extra`; definitions in the script's docstring). One row per cell size (each Lxy slice of table 16 cut in 1, 2, 4 or 8), set of gen dR edges (those of table 18, or the same with the bin below 0.002 cut at 0.0005 and 0.001) and minimum number of LJs of a mass point in a cell (100 or 50). |
| 22 | The same test for each content separately (Lxy slices cut in 4 and 8, at least 50 LJs), with the cells alone and with the charge, the region, the lower generated pT or the class of `analyze.py` as a further cell variable. |

- Opposite-charge pairs have no stored vertex in 0.067 of the PF-PF LJs (2360 of 35,444), 0.446
  of the PF-DSA LJs (7246 of 16,256) and 0.357 of the DSA-DSA LJs (9341 of 26,167) (table 15).
- Against the gen dR (table 18) the PF-PF fraction is 0.165 below 0.002, 0.108 in
  [0.002, 0.005) and 0.011 to 0.022 from 0.005 to 0.2. With a DSA muon in the pair it is 0.39 to
  0.52 below 0.02 and 0.07 to 0.23 from 0.02 to 0.2 (67 PF-DSA and 352 DSA-DSA LJs in
  [0.1, 0.2)).
- The content follows the decay radius (table 17): at 500 GeV the PF-PF share is 0.863 within
  1 cm and 0.011 in [60, 110) cm, where DSA-DSA holds 0.818. It also changes with the mass at a
  given radius. Within 1 cm the PF-DSA share is 0.115 at 500 GeV (1179 of 10,266) and 0.179 at
  1000 GeV (1284 of 7173) (the 100 GeV samples have 39 LJs there), and there 0.460 of the PF-DSA
  LJs have no vertex (1134 of 2464) against 0.079 of the PF-PF LJs (1153 of 14,628; table 16).
- For DSA-DSA pairs of opposite charge the no-vertex fraction falls with the decay radius: 0.502
  within 30 cm (129 of 257), 0.414 in [60, 110) cm, 0.277 in [150, 200), 0.048 in [250, 300) and
  0.018 in [300, 350) (3 of 168), while the pass fraction rises from 0.471 to 0.911 (table 19).
  In the barrel it is 0.588 (97 of 165), 0.447, 0.292, 0.045 and 0.020 in those slices, and in
  the endcap 0.341 (28 of 82), 0.242, 0.158 and 0.042 (2 of 48) in the first four; the barrel is
  above the endcap out to 200 cm (table 20). Neither the fall with the radius nor the
  difference between barrel and endcap is explained here.
- In cells of (Lxy, gen dR) the mass points do not share one no-vertex fraction (table 21).
  With the Lxy slices of table 16, the dR bins of table 18 and mass points of at least 100 LJs,
  chi2 = 167.3 for 33 degrees of freedom (97.8 for 41 with the content as a third cell
  variable, 132.7 for 43 with the lower generated pT of the two daughters instead, 84.9 for 41
  with both), and the spread of the fraction between mass points beyond the binomial
  expectation is 0.028. Cutting each Lxy slice in 2, 4 and 8 leaves that spread at 0.027 to
  0.029, so it is not made by different Lxy distributions inside a slice. Two changes reduce
  it: the content as a further cell variable (0.019, 0.017, 0.017 and 0.013 for the four cell
  sizes), and the dR bin below 0.002 cut at 0.0005 and 0.001 (0.021, 0.020, 0.020 and 0.014).
  With those finer dR bins the content adds nothing (0.016 to 0.021).
- Per content (table 22) the PF-PF mass points agree: chi2 = 54.2 for
  42 degrees of freedom with the Lxy slices cut in 4 and the dR bins of table 18, 36.6 for 42
  with the finer dR bins, 38.8 for 38 and 39.7 for 39 with the slices cut in 8. The PF-DSA
  mass points agree in every case (15.9 for 20, 10.3 for 17). The DSA-DSA mass points do not:
  74.6 for 30, 72.5 for 26, 47.3 for 21 and 44.3 for 19, a spread of 0.029 to 0.035, and it
  stays at 0.027 to 0.038 with the charge, the region, the lower generated pT or the class as
  a further cell variable. What separates the mass points for DSA-DSA pairs is not identified
  by these tables.

The producer of the vertex tables (`MuonVertexTableProducer` of `cms-sidm/LLPNanoAOD`, read at
the head of its main branch, commit 14c5287454) writes no vertex for a pair in two cases: the
Kalman vertex fit of the two tracks is not valid, or the distance of closest approach of the
two tracks, where its calculation succeeds, is above 15 cm. Which case applies to an LJ is not
stored. That the 2018 v10 samples were made with this version of the producer was not checked.

These tables have the limits of the caveats above: sixteen 2Mu2E samples, three mediator
masses, two lifetimes per mass point, raw counts. They give fractions with their counts and no
errors.

#### Widest generated pairs (`wide_pair.py`, `wide_pair/`)

In `vtx_outcome_lxy_dr/counts_gen_dr.csv.gz` (m_Zd = 5 GeV, reference channel), at Lxy >= 100 cm
the LJs whose generated pair has gen dR in [0.2, 0.5) lose the vertex more often than those in
[0.1, 0.2): 0.164 against 0.076 at 100 GeV, 0.192 against 0.078 at 150 GeV, 0.208 against 0.083
at 200 GeV (both final states summed, `T0_counts_gen_dr_by_final_state.csv`). This study asks
which property of the LJ goes with it, from the per-LJ truth records.

**Samples.** m_Zd = 5 GeV, mediator 100, 150 and 200 GeV, the three longest lifetimes of each
(100 GeV: 40, 200, 400 mm; 150 GeV: 27, 130, 270 mm; 200 GeV: 20, 100, 200 mm), in both final
states: 2Mu2E (`condor/samples_wide_pair_2mu2e.txt`, 19 jobs in
`condor/job_args_wide_pair_2mu2e.txt`) and 4Mu (`condor/samples_4mu.txt`, 23 jobs in
`condor/job_args_4mu.txt`). Signal MC only, reference channel `base_ljDisplacementIso_cosAlpha`.

**Records.** In `extract.py`, d1 and d2 are the final-state (status-1) daughter muons of the dark
photon matched to the LJ, and o1, o2 (with `n_other_dau`, and `o_a_lxy`, `o_a_phi` of that dark
photon) the Zd daughter muons of the event that are not d1, d2. In 2Mu2E d1 and d2 are the
event's only two Zd daughter muons, so the committed tables do not change (`analyze.py
--propagate` and `pair_content.py` reproduce all 25 committed outputs byte for byte from the new
records). In 4Mu, where the event has four Zd daughter muons, d1 and d2 are those of the dark
photon matched to the LJ and o1, o2 those of the other one. Absent generated muons have pt -1,
eta and phi -99 and charge 0 (d1, d2, o1, o2); without another dark photon `o_a_lxy` is -1 and
`o_a_phi` -99.

`extract.py` now also writes, per muon of the LJ, the DSA track parameters at the reference point
and the DSA track-quality branches (the list `DSA_TRACK`), the PF muon's inner-track reference
point and momentum and its muon-type flags (`PF_TRACK`), each -999 for a muon of the other type,
and for the two leading muons the stored vertex of that pair (`pv_n`, `pv_dca`, `pv_dcaStatus`,
`pv_isValid`, `pv_normChi2`; `pv_n` is 0 and the others -999 without one). No table in this
directory uses these fields.

From the repository root:

```bash
D=sidm/studies/vertexChi2/dr_truth_match
python condor/make_job_args.py --samples-file $D/condor/samples_4mu.txt --dataset llpNanoAOD_v2 \
  --location-cfg signal_4mu_v10.yaml --files-per-job 25 --replace-xcache \
  --census-skip sidm/studies/vertexChi2/signal_skip.json \
  --outdir $D/condor/filelists --job-args /tmp/job_args_4mu.txt
while read s c f; do
  python $D/extract.py $s $D/condor/$f $REC/${s}_${c}.parquet --chunksize 50000 --workers 1
done < $D/condor/job_args_4mu.txt
python $D/analyze.py $REC $AN --propagate
python $D/wide_pair.py $AN sidm/studies/vertexChi2/vtx_outcome_lxy_dr/counts_gen_dr.csv.gz \
  $D/wide_pair/4mu
```

and the same for 2Mu2E with `samples_wide_pair_2mu2e.txt`, `signal_2mu2e_v10.yaml`,
`job_args_wide_pair_2mu2e.txt` and `$D/wide_pair/2mu2e`. Rerun on 2026-11-14, the make_job_args
command writes the same filelists (23 and 19 files, identical byte for byte) as those the records
were made from, and the job lists here are its output with the filelist paths relative to
`dr_truth_match/condor/`. The records were made by hand on cmslpc, not on condor. `analyze.py`
stops at table 7, which needs the 500 and 1000 GeV samples, after writing `lj_truth.parquet`,
`lj_classified.parquet` and `table3.csv`, which `wide_pair.py` reads. Its match radii, computed
from these samples: 2Mu2E PF 0.0510, DSA 0.1106; 4Mu PF 0.0607, DSA 0.1355.

Terms:

- **narrow** = gen dR [0.1, 0.2), **wide** = gen dR [0.2, 0.5); no vertex = `vtx_chi2 < 0`.
- **gen dR** here is that of the final-state daughters d1, d2. `counts_gen_dr.csv.gz` bins the
  dark photon's own children (`hists.mu_lj_gen_a_mumu_dr`), before final-state radiation; see
  "Records against counts_gen_dr.csv.gz" below.
- **soft_dau_unreco**: no muon of the LJ, PF or DSA, within its radius of the softer daughter
  (lower generated pT of d1, d2), compared with the daughter's direction at the point of closest
  approach to the beamline (the propagated match of `analyze.py`).
- **acceptance**: both daughters in abs(eta) < 2.4 and the softer above 10 GeV of generated pT.
- **lead_on_other**: m0 or m1 within its radius of a daughter of another dark photon (o1, o2,
  propagated from that dark photon's Lxy and phi).
- **Standardisation**: the wide bin's no-vertex fraction predicted from the narrow bin's fraction
  in each cell of a variable, weighted by the wide bin's cell shares. **gap closed** =
  (predicted - narrow) / (observed wide - narrow).
- **Errors**: `err_*` are binomial with the wide bin's cell shares held fixed; `*_bs_err` are the
  standard deviations of a bootstrap (1000 replicas, fixed seed, the LJs of the narrow and of the
  wide bin resampled separately), and `gap_closed_bs_q16`, `_q84` its 16 % and 84 % quantiles.
  Values below quote gap closed +- its bootstrap standard deviation, or its 16 to 84 % range where
  that is asymmetric.

| file (in `wide_pair/2mu2e/` and `wide_pair/4mu/`) | what it holds |
|---|---|
| `T0_counts_gen_dr_by_final_state.csv` | From `counts_gen_dr.csv.gz`: no-vertex fraction of narrow and wide, Lxy [10, 100) and [100, 1000), per mediator mass, for both final states, each final state, and the samples of this run. |
| `T1_records_vs_counts_gen_dr.csv`, `T1_records_vs_counts_gen_dr_by_sample.csv` | The records binned as `counts_gen_dr.csv.gz`, per sample, Lxy, gen dR and outcome, with the differences. |
| `T2_composition.csv` | Per Lxy range, mass and gen dR bin ([0.05, 0.1), narrow, wide): LJs, no-vertex fraction, class shares A to E, content shares, same-charge share, median soft reco and gen pT, generated acceptance (both daughters in abs(eta) < 2.4 and softer pT > 5 or > 10 GeV), soft_dau_unreco, lead_on_other. |
| `T2_dr_to_other.csv` | Quantiles (minimum, 1, 5, 16, 50 %) of the smallest dR between m0, m1 and a daughter of the other dark photon. |
| `T3_cells.csv`, `T3_summary.csv`; `T3_cells_lxy10_100.csv`, `T3_summary_lxy10_100.csv` | Standardisation by one variable at a time (class, content, number of muons, charge, soft reco pT, soft gen pT, region), Lxy [100, 1000) and [10, 100). |
| `T4_cells.csv`, `T4_summary.csv` | The two variables of T3 closing most of the gap together (4Mu: class x soft gen pT; 2Mu2E: soft gen pT x class). |
| `T4_alt_cells.csv`, `T4_alt_summary.csv` | class x charge. |
| `T5_cells.csv`, `T5_summary.csv` | soft_dau_unreco alone, both Lxy ranges, per mass. |
| `T5_acceptance_cells.csv` | soft_dau_unreco x acceptance per bin: LJs, share of the bin, no-vertex fraction. |
| `tables.txt` | All tables as printed. |

**4Mu, Lxy >= 100 cm, masses summed**

- Narrow 0.0834 (4994 of 59,863), wide 0.178 (2317 of 13,003): a gap of 0.095 (`T5_summary.csv`,
  `T2_composition.csv`).
- The softer daughter has no LJ muon within its radius in 0.054 of the narrow and 0.308 of the
  wide LJs; per mass 0.054 to 0.251 (100 GeV), 0.054 to 0.339 (150), 0.054 to 0.440 (200)
  (`T2_composition.csv`).
- No-vertex fraction, narrow to wide: with the softer daughter found 0.071 (4008 of 56,613) to
  0.097 (873 of 9000); without it 0.303 (986 of 3250) to 0.361 (1444 of 4003) (`T5_cells.csv`).
- Standardised by soft_dau_unreco the wide bin is predicted at 0.142 +- 0.003, closing
  0.62 +- 0.03 of the gap; per mass 0.54 +- 0.12, 0.59 +- 0.05, 0.74 +- 0.07 (`T5_summary.csv`).
- Missing daughter and acceptance (`T5_acceptance_cells.csv`): of the 4003 wide LJs with the
  softer daughter missing, 2820 (0.70) are outside the acceptance (no-vertex fraction 0.385) and
  1183 (0.30; 0.091 of the bin) inside it (0.303). In the narrow bin 544 of 3250 (0.17) are
  outside. The wide LJs outside the acceptance with the daughter found, 1323 (0.102 of the bin),
  lose the vertex at 0.106, those inside it with the daughter found at 0.096. So most of the
  missing-daughter share of the wide bin (0.217 of its 0.308) is a daughter outside the acceptance,
  and the no-vertex fraction goes with the missing daughter, inside the acceptance or not.
- Other variables (`T3_summary.csv`): class closes 0.44 +- 0.02 (0.26 +- 0.07, 0.37 +- 0.03,
  0.51 +- 0.05 per mass); soft gen pT 0.36 +- 0.04; content 0.18 +- 0.01; charge 0.17 +- 0.01;
  soft reco pT 0.11 +- 0.01; region 0.03 +- 0.004; number of muons -0.05 +- 0.004. class x charge
  0.47 +- 0.02 (`T4_alt_summary.csv`); class x soft gen pT 0.60 +- 0.05, with 3 of 15 cells on
  narrow cells under 50 LJs, 0.023 of the wide bin on narrow cells under 20 and 0.009 without a
  narrow cell (`T4_summary.csv`).
- Class A share 0.841 to 0.548, class D 0.103 to 0.343 (`T2_composition.csv`); no-vertex
  fraction in class A 0.063 (3156 of 50,319) to 0.097 (692 of 7125), class D 0.220 (1353 of
  6139) to 0.307 (1367 of 4460) (`T3_cells.csv`).
- The softer generated pT has median 23.9 GeV (narrow) and 12.1 GeV (wide); in the acceptance:
  0.979 and 0.681 (`T2_composition.csv`).
- lead_on_other is 0 in every row (`T2_composition.csv`). That is fixed by the event topology:
  the smallest dR between m0, m1 and a daughter of the other dark photon is at least 0.39 in
  every 4Mu LJ, 1.84 at the 1 % quantile and 2.98 at the median (`T2_dr_to_other.csv`), far
  beyond the match radii (at most 0.14). It shows that the other dark photon is far away, nothing
  more.

**4Mu, Lxy [10, 100) cm.** Narrow 0.0488 (8145 of 166,759), wide 0.0553 (2772 of 50,164), a gap
of 0.0064. The soft_dau_unreco share rises from 0.031 to 0.094 (`T2_composition.csv`), and the
standardised prediction, 0.0619 +- 0.0008, is above the observed 0.0553 +- 0.0010
(`T5_summary.csv`; gap closed 2.04, 16 to 84 % range 1.76 to 2.51): with the softer daughter
found the fraction falls from 0.042 to 0.033, without it 0.251 to 0.268 (`T5_cells.csv`). Of the
4703 wide LJs with the softer daughter missing, 3623 (0.77) are outside the acceptance
(`T5_acceptance_cells.csv`).

**2Mu2E.** At Lxy >= 100 cm narrow 0.0615 (516 of 8385), wide 0.218 (139 of 638), a gap of 0.156;
at [10, 100) cm 0.0376 (719 of 19,127) and 0.217 (127 of 586) (`T3_summary.csv`,
`T3_summary_lxy10_100.csv`). soft_dau_unreco 0.049 to 0.492 at >= 100 cm and 0.033 to 0.369 at
[10, 100) cm (`T2_composition.csv`).

- Standardised by soft_dau_unreco: 0.140 +- 0.010 at >= 100 cm, closing 0.50 +- 0.08 of the gap
  (100 GeV 0.62, range 0.36 to 0.90; 150 GeV 0.33 +- 0.07; the 200 GeV value rests on 47 wide
  LJs, range 0.62 to 2.40, and is not quoted), and 0.097 +- 0.007 at [10, 100) cm against the
  observed 0.217 +- 0.017, closing 0.33 +- 0.05 (`T5_summary.csv`).
- No-vertex fraction, narrow to wide, at >= 100 cm: with the softer daughter found 0.053 (422 of
  7975) to 0.049 (16 of 324), without it 0.229 (94 of 410) to 0.392 (123 of 314). At
  [10, 100) cm: with the softer daughter found 0.032 (589 of 18,504) to 0.124 (46 of 370), without
  it 0.209 (130 of 623) to 0.375 (81 of 216) (`T5_cells.csv`).
- Missing daughter and acceptance: of the wide LJs with the softer daughter missing, 169 of 314
  (0.54) at >= 100 cm and 116 of 216 (0.54) at [10, 100) cm are outside the acceptance, with
  no-vertex fractions 0.49 (outside) and 0.28 (inside) at >= 100 cm (`T5_acceptance_cells.csv`).
- By class 0.39 +- 0.07 (>= 100 cm) and 0.22 +- 0.03 ([10, 100) cm) (`T3_summary.csv`,
  `T3_summary_lxy10_100.csv`). By soft gen pT 0.96 (range 0.71 to 1.23) at >= 100 cm is not
  quoted as a result: 0.26 of the wide bin sits on narrow cells under 20 LJs and 0.024 has no
  narrow cell; the same holds for soft gen pT x class (1.01, 0.25 on cells under 20, 0.033 without
  a cell; `T4_summary.csv`). At [10, 100) cm soft gen pT closes 0.70 +- 0.17, with 0.022 of the
  wide bin without a narrow cell.
- In `counts_gen_dr.csv.gz` the 2Mu2E wide pairs lose the vertex more often already at
  [10, 100) cm (0.194, 0.230, 0.312 against 0.055, 0.037, 0.033), 4Mu wide pairs do not (0.063,
  0.050, 0.042 against 0.048, 0.046, 0.046) (`T0_counts_gen_dr_by_final_state.csv`).

**Records against counts_gen_dr.csv.gz** (`T1_records_vs_counts_gen_dr.csv`)

Every LJ is in both (330,322 in 4Mu, 34,325 in 2Mu2E), but the gen dR bin differs for a few: of
the 1080 cells (sample x Lxy x gen dR x outcome) 956 agree exactly in 4Mu (largest difference
106 LJs) and 1036 in 2Mu2E (largest 26). The records use the final-state daughters, the histogram
the dark photon's own children before final-state radiation. In 4Mu at Lxy [10, 100) cm the
records hold 266 more LJs in gen dR [0.05, 0.1) and 201 fewer in [0.2, 0.5); at >= 100 cm the wide
bin holds 13,003 LJs against 13,050 (2Mu2E: 638 against 640). A test extraction of
4Mu_200GeV_5p0GeV_20p0mm that also stored `hists.mu_lj_gen_a_mumu_dr` (not a staged file)
reproduced all 120 cells of that sample with it, against 99 of 120 with the final-state
daughters; the two dR differ by more than 10 % in 0.75 % of its LJs, and in 70 % of those the
final-state dR is the smaller.

**What is and is not explained**

- In 4Mu at Lxy >= 100 cm the softer daughter is missing from the LJ far more often in the wide
  bin (0.308 against 0.054), and those LJs lose the vertex at 0.30 to 0.36. Standardising for it
  alone closes 0.62 +- 0.03 of the gap (0.54 to 0.74 by mass); the class of `analyze.py` closes
  0.44 +- 0.02. About 70 % of the missing daughters are outside the acceptance (0.217 of the wide
  bin, against 0.319 of the wide bin outside the acceptance in all).
- Inside each group the wide bin is still higher (0.071 to 0.097 with the softer daughter found,
  0.303 to 0.361 without); that part is not explained.
- In 2Mu2E at Lxy >= 100 cm the same variable closes 0.50 +- 0.08 of a larger gap (0.156); there
  the LJs with the softer daughter found do not rise (0.053 to 0.049), and the rest of the gap is
  inside the group without it (0.229 to 0.392).
- At Lxy [10, 100) cm the two final states differ. In 4Mu the missing-daughter share rises in the
  wide bin without a matching rise of the no-vertex fraction (prediction 0.062 above observation
  0.055). In 2Mu2E the prediction is below the observation (0.097 against 0.217, closing 0.33),
  and the no-vertex fraction rises inside the group with the softer daughter found (0.032 to
  0.124, 46 of 370 wide LJs).
- The standardisation is arithmetic over two groups and does not show that the missing daughter
  causes the vertex loss.
- lead_on_other = 0 says only that the other dark photon's daughters are far from the LJ's
  leading muons (dR >= 0.39 in every 4Mu LJ).

**Limits**

- The match radii are computed from these samples by `analyze.py` and differ from those of the
  committed tables (PF 0.0556, DSA 0.0761); soft_dau_unreco and the classes depend on them.
- Values per mass rest on few wide LJs at 200 GeV in 2Mu2E (47 at >= 100 cm, 45 at [10, 100) cm)
  and on 532 (2Mu2E) and 2146 (4Mu) narrow LJs at 100 GeV, >= 100 cm; read them with the
  bootstrap ranges.
- A daughter reconstructed as a muon outside the LJ, or matched beyond its radius, counts as
  missing.
- The gen dR of the records is the final-state one, so a few LJs sit in a neighbouring gen dR bin
  of `counts_gen_dr.csv.gz` (above).

### Slides (`slides/`)

`slides/dsa_tail.tex` is a group-meeting deck on sections 8 and 10 of the write-up, with context
from sections 3 to 9. Every number on it is quoted from the write-up's markdown, the records in
`transfer_factor/` or the fold logs in `transfer_factor/logs/`, and each slide names its source.
From `slides/`, `python extract_figures.py` writes the write-up's own PNG outputs to `figures/`
(not tracked), named by section and cell; then `pdflatex dsa_tail.tex`, run twice, builds
`dsa_tail.pdf` (not tracked).

## Provenance sidecars

Each `.coffea` is written with a `<sample>.meta.yaml` next to it via `sidm.tools.metadata`,
recording the full selection definitions, the hist-collection contents, the input ROOT file
list, the per-sample cross section, the SIDM git commit, the coffea version, the schema and
chunksize, and a UTC timestamp. Read one back with

```python
from sidm.tools.metadata import load_run_metadata
meta = load_run_metadata("<path to the .coffea>")   # it derives the .meta.yaml path itself
```

The outputs that already existed when the sidecar write was added carry **back-filled**
sidecars, written without reprocessing by a one-off back-fill script that is not in the
repository. They are labelled as such: they contain a `backfill:` block naming where each field
came from, and their `sidm_commit` is `null`, because the commit that produced them was never
recorded and is not recoverable. Eight output directories carry back-filled sidecars: those of
the first six run notebooks in the table above, and `vtxTail_2018C` and `vtxTailMC_2018`,
superseded first passes that no current run notebook produces. The sidecars of those two say so
and derive the hist-collection list from the outputs themselves.

## Run order

1. `runDataBg2018.ipynb`: the control-region yields and the `vtx_chi2` threshold scan.
2. `runDataCosmic2018.ipynb`: the spread-inverted cosmic region, used as a `vtx_chi2` cross-check.
3. `runDataTail2018.ipynb`: the tail histograms and the DSA-PF response, data.
4. `runMCTail2018.ipynb`: the same histograms in MC. Must use the same channel list as 3.
5. `runSignalVtx.ipynb`: the 180-sample signal efficiency of both vetoes.
6. `runSignalVtxNhits.ipynb`: the 26-sample nHits subset. Needs the 2D nHits and \|eta\|
   histograms, which are in `mu_lj_vtx_chi2_dsa_pt`.
7. `vertexChi2CosmicVetoStudy.ipynb`: the write-up, once 1 to 6 are on disk.

Steps 1 to 6 are independent of each other and can be run in any order or in parallel; step 7
needs all of them. Signal samples read `../../configs/ntuples/signal_4mu_v10.yaml` and
`signal_2mu2e_v10.yaml`. `runSignalVtx.ipynb` and `runSignalVtxNhits.ipynb` also read a JSON
veto list of corrupt input files and append to it any file that fails to read (17 files across
15 samples after both runs).

Steps 1, 3, 4 and 5 take hours each on a hundred Condor workers. All steps require a valid VOMS
proxy readable by the Condor schedd, and `replace_xcache=True` in `make_fileset`, since the
location YAMLs carry `root://xcache//` URLs that only resolve on coffea-casa. The proxy path
(`X509_USER_PROXY`) and the veto-list path (`SKIP_JSON`) in the run notebooks are the author's;
whoever reruns them sets their own.

Section 10 of the write-up (the transfer-factor fold) needs four more runs, independent of 1-7
above and of each other: `runDataFull2018.ipynb` + `runMCRespFull.ipynb` (the dR-matched
comparison) and `runDataSegMatch2018.ipynb` + `runMCRespSegMatch2018.ipynb` (the
segment-matched control, the primary input). All four read the full 2018 sample (all 68
skimmed 2018A-D `DoubleMuon` samples on the data side), not 2018C alone. The fold scripts in
`transfer_factor/` run on their `.coffea` outputs and write the three JSON records there, which
the write-up re-reads (`transfer_factor/README.md` has the command lines); the write-up also
re-opens the segment-matched data `.coffea` files directly for its `phi` x `nHits` and
`|eta|` x `nHits` figure. Rerunning the notebooks and then the fold scripts therefore
reproduces every input of the write-up.

## A note on the two MC sets

Two different MC sets appear in the study. The distinction matters when reading it, because
several headline numbers change completely between them:

- The **full background stack** (output directory `vtxChi2_vs_spread_2018C`) is what the
  control-region yields and the threshold scan use.
- The **eight-sample tail set** (output directory `vtxTailMC2_2018`) carries only the eight
  samples that populate the >= 200 GeV DSA tail, but is the only MC set with the 2D
  histograms.

Numbers that must always carry the name of the set they came from:

| quantity | full stack | eight-sample tail set |
|---|---|---|
| pre-veto MC tail (`pT >= 200` GeV) | 293.6 | 133.7 (46% of the stack's) |
| post-veto MC tail | 133.7 | 129.5 (97% of the stack's) |
| fraction of the MC tail the veto removes | **54%** | **3%** |
| pre-veto tail data/MC | 8.22 | 18.06 |
| post-veto tail data/MC | 2.31 | 2.385 +/- 0.426 |
| MC yield over all pT, Spread channel | 108137.7 | 57091.1 (52.8% of the stack's) |

The two sets agree after the veto because the tail set covers 97% of the MC tail there. They
disagree before it because the full stack's pre-veto tail contains one high-weight
`QCD_Pt30To50` event worth about 156 weighted entries, which the veto removes. (That the
eight-sample set's pre-veto tail, 133.7, equals the full stack's post-veto tail, also 133.7, is
a coincidence of two unrelated quantities.)

The tail set does **not** cover the bulk (the low-pT reference): over all pT it carries 52.8% of
the full stack, the missing half being the soft-QCD bins. Absolute data/MC in the bulk is
therefore meaningless with it, and so is anything built from its total yield; only shapes across
bins are used there. Section 8 of the write-up marks the one place where that restriction bites,
and gives the argument in a form (a ratio of rates) that does not depend on it.

## One more normalization caveat

The control region itself is not normalized, and no scale factor is applied anywhere in this
study. In the Before channel (no cosmic veto) it holds 0.770 data mu-type lepton jets per MC
one, but only 0.556 data DSA muons per MC one. The ratio of those two is a DSA-muon
multiplicity mismodelling: MC carries 0.351 DSA muons per mu-type lepton jet against 0.254 in
data, 38% more. Veto efficiencies and the azimuthal cosmic test are ratios within one sample and
are untouched by it. Absolute tail data/MC ratios are not, which is why the write-up's
load-bearing tail statement is a ratio of rates rather than the raw 2.4 +/- 0.4 of the
eight-sample tail set. Section 3.1 of the write-up sets this out in full.

## The dimuon-vertex-fit variable and its cuts

`mu_ljs.vtx_chi2` is the best (smallest) `normChi2` of any stored dimuon vertex built from two
of the LJ's own muons, taken from the LLPNanoAOD `DSAMuonVertex`, `PatMuonVertex` and
`PatDSAMuonVertex` tables. The matching lives in
[lj_vertex_chi2.py](../../tools/lj_vertex_chi2.py). `SidmProcessor` fills the field only when a
requested channel cuts on `vtx_chi2` or a requested histogram reads it (those histograms are
listed in `vtx_chi2_or_raw_dsa_hists` in `sidm/definitions/hists.py`). Otherwise it does not
load the vertex tables, and the LJs carry no `vtx_chi2` field.

A value of -1 means no such vertex was stored, and it *fails* the keep-cuts by design. The
producer only writes pairs whose fit converged, so an LJ whose muons never form a vertex is the
signature of a fake or cosmic pairing. The variable is only defined for LJs with at least two
muons, so pair any cut or plot with `Mu >= 2`.

One asymmetry matters when comparing with the `*Spread_*` cuts. The Spread variables are a
maximum over the LJ's muon pairs (the worst pair must be compact), while `vtx_chi2` is a minimum
(one good dimuon suffices). The two coincide only for LJs with exactly 2 muons. For LJs with
>= 3 muons the vertex cut is the looser criterion by construction: it asks "does this LJ contain
a real dimuon", not "is every pairing consistent". The `mu_lj_muonN` histogram measures how much
of a given region has >= 3 muons.

The cuts in `cuts.py` are `0 <= vtx_chi2 < 5` and its complement at the mu-LJ-cut level,
`no vtx or vtx_chi2 >= 5`. Four channels use them, mirroring the `vxySpread`-based Spread-veto
channels (a `dxySpread` variant of the Spread-veto channel was also run,
`cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu_dxy`). These channels, and every other channel of
this study, are defined in
[`sidm/configs/selections_vertexChi2_as_run.yaml`](../../configs/selections_vertexChi2_as_run.yaml),
not in `selections.yaml`:

- `cr_1egmOr1mu_lj_vtxChi2_CosMu`: `cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu` with the four
  Spread cuts replaced by `0 <= vtx_chi2 < 5`
- `cr_1egmOr1mu_lj_vtxChi2`: the vertex cut without the cosAlpha event cut
- `cosmic_muons_vtxChi2`: the **vertex-inverted** cosmic region, `no vtx or vtx_chi2 >= 5` with
  `lj_iso < 0.2` and the cosAlpha cut inverted. It mirrors `cosmic_muons`, the
  **spread-inverted** cosmic region, which inverts the Spread cuts instead and so uses no vertex
  information anywhere in its definition. The two are different selections, they are named apart
  wherever either is used, and a template measured in one does not transfer to the other.
- `base_ljDisplacementIso_cosAlpha_vtxChi2`: the signal channel,
  `base_ljDisplacementIso_cosAlpha` with `0 <= vtx_chi2 < 5` on the mu-type LJs

The distributions are `mu_lj_vtx_chi2` (linear axis) and `mu_lj_vtx_chi2_log` (log axis, for
threshold scans), both in the `mu_lj_vtx_chi2_dsa_pt` hist collection. On the linear axis the -1
entries sit in their own bin at the left edge, and because `PatMuonVertex` normChi2 spans ~8
decades, a large fraction of PF-PF entries lands in the overflow. The threshold 5 was not
optimized; section 3.2 scans it.

When it is computed, `vtx_chi2` is all-or-nothing. Unless all three vertex tables are present in
the input, it is NaN and every vtx_chi2 cut, keep and inverse alike, selects nothing; the
processor then prints `vtx_chi2 disabled: missing vertex table(s) ...`. Without this rule a
partial table set would put whole LJ categories into the inverse-cut region as fake "no vertex"
entries. An exactly empty vtxChi2 channel therefore means missing tables; check the worker
logs. The samples in `sidm/configs/ntuples/` (`signal_2mu2e_v10.yaml`, `signal_4mu_v10.yaml`,
`data_skimmed.yaml`, `backgrounds.yaml`) all carry them. The brute-force closure test lives at
[tests/test_lj_vertex_chi2.py](../../../tests/test_lj_vertex_chi2.py). It also reproduces the
stored `dR` of every PF-PF vertex in the two NanoAOD files in `tests/data/` (496 rows) from the
PF muons its indices point to.

The DSA-leg index convention is supported, but not closed by a `dR` test. In the same two
files, for 2096 of 2379 DSA legs (88%, against 31% by chance) the DSA muon nearest in direction
to the leg's stored refitted track is the one its index points to. The stored `dR` of the 1884
vertices with a DSA leg, however, matches neither the NanoAOD DSA-muon directions nor the
refitted tracks (none within 4e-4), so it is computed from a track state NanoAOD does not store.

## Result of the 2018C comparison (August 2026)

Inputs are 2018C `DoubleMuon`, the 2018 background stack (MC scaled by 6.97/59.83, the 2018C
share of the 2018 integrated luminosity) and the full 180-sample v10 signal grid.

The vertex veto tags the same cosmics as the Spread veto and costs more signal, so the Spread
veto with the cosAlpha event cut stays as the selection. `0 <= vtx_chi2 < 5` flags 98% of the
mu-type lepton jets in the spread-inverted cosmic region, which is defined without any vertex
information. But it has no data-versus-MC separating power: the data/MC ratio of the lepton jets
it keeps is flat at 0.75 across four decades of threshold, against 0.77 with no cut. It leaves
the same residual data/MC above 200 GeV in DSA-muon pT (2.38 against 2.31 for the Spread veto,
both on the full background stack). And over the signal grid it keeps 28% fewer signal lepton
jets (0.690 vs 0.961). The signal efficiencies are for the Spread veto in its as-run `vxySpread`
form; `main`'s `dxySpread` form was not run on signal. The one corner where the ordering
reverses is m_Zd = 5 GeV, mediator mass 100-150 GeV, ctau >= 130 mm: 8 of the 180 signal
samples when efficiency is counted per DSA muon (counted per lepton jet it is 6 samples, one
of them at ctau = 40 mm). There the
max-over-pairs Spread variables start failing on wide displaced pairs and the min-over-pairs
vertex requirement does not.

The second half of the study characterizes what survives either veto. The Spread veto removes
87% of the >= 200 GeV DSA-muon tail in data (2414 -> 309). In MC it removes 54% of the full
background stack's tail (293.6 -> 133.7) and 3% of the eight-sample tail set's (133.7 -> 129.5).
These are two different MC sets, and neither number travels without its label. The removed
component is unambiguously cosmic: the pre-veto phi vertical/horizontal ratio is 12.2 +/- 1.3
against an MC expectation of 1.3 +/- 0.2, the tail is 93% barrel against 71%, and two
independent routes put the pre-veto tail's cosmic content at 87% and 84%. What survives is not
cosmic: the post-veto ratio is 0.94 +/- 0.15, bounding a residual cosmic contribution below
6.2% of the tail at 95% CL, a bound that moves only to 7.8% if the template's vertical fraction
is degraded from 0.99 to 0.90.

The surviving tail exceeds MC by 2.4 +/- 0.4 (eight-sample tail set), but that ratio is not the
load-bearing number: the control region is not normalized, MC carrying 38% more DSA muons per
mu-type lepton jet than data, and no scale factor is applied. The statement that carries is a
ratio of rates, which is immune to that mismatch. The **feed-up** rate is the rate at which a
muon with PF pT below 100 GeV is reconstructed as a DSA muon at or above 200 GeV. An independent
DSA-to-PF matched control gives a data/MC feed-up ratio of 3.47 +/- 0.25, and the post-veto
tail rate divided by that feed-up rate is 0.95 +/- 0.06 in data and 0.99 +/- 0.18 in MC, both
consistent with one, against 4.81 +/- 0.13 for data before the veto. That is a consistency
statement of about factor-1.5 precision, not a demonstration that the two agree: the control
rate is measured before the veto on a different population, is not flat in |eta| or nHits, and
the closure moves between 0.70 and 1.41 in data (0.57 and 1.35 in MC) as the pT_PF window moves.
The surviving tail therefore stays as sections 5-8 leave it: not cosmic (below 6.2% at 95% CL),
and of the size DSA momentum mismeasurement would produce at roughly factor-1.5 precision.

The recommendation is to keep Spread + cosAlpha (`main`'s baseline since PR #56), not to adopt
the vertex veto, and to add no further object cut. `nHits >= 20` (nHits: DT plus CSC hits on the
DSA track) removes 41% of the post-veto data tail (309 -> 183) and takes its data/MC
(eight-sample tail set) from 2.39 to 1.61. But it costs 15% of the signal DSA muons at all pT
(8% of those already in the tail, a number not to quote on its own), and it is itself
mismodelled: it removes 24.5% of the data bulk against 14.1% of the MC bulk.
`sigma(pT)/pT < 0.3` moves the ratio the wrong way (it keeps 0.084 of the data tail against
0.145 of the MC tail) while costing 44% of the signal at a 1 TeV mediator. The section-10
transfer factor fails its MC self-closure, so no correction from the DSA-PF control is
available; the residual DSA high-pT rate has to be set from data in the background estimate.
The vertex-veto code, its channels and its histograms stay in the repository as a documented
negative result.