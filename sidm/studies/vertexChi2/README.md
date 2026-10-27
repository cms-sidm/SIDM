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
control sample of DSA muons matched to PF muons and folds it (applies it bin by bin) onto the
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
"Notebooks" confirm the definitions independently for two directories, and for a third through
its ZZ sample.

**Why the channels are kept as run.** The study compares cosmic-veto options. It therefore needs
a baseline with no cosmic veto at all, and each veto applied on its own to that same baseline.
On `main` the baseline itself has changed, so `main`'s channels cannot provide the no-veto
reference, and the channels are kept here as they were run. How the as-run baseline differs
from `main` today:

- Unchanged: the eta-phi veto on all four LJ sources (electrons, photons, PF muons, DSA muons)
  and the event cut `>=2 muons (pf or dsa) pT > 26 GeV`.
- As run, the PF muons at the LJ source have no muon-chamber hit requirement. `main` now also
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

Every run notebook passes `selections_cfg="configs/selections_vertexChi2_as_run.yaml"` to
`SidmProcessor` and to `write_run_metadata`. That argument was added after the notebooks were
executed, when these channels were still defined in `sidm/configs/selections.yaml`; the stored
outputs are those of the original runs. With the argument, the notebooks rerun on main with the
selections they were run with. Three stored outputs were rerun this way: one 2018D data file of
the full-2018 segment-matched run (`DoubleMuon_2018D_22` of `vtxResp2_2018Full`), one signal
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
| `vertexChi2CosmicVetoStudy.ipynb` | The write-up. Reads the ten output directories above and the three JSON records in `transfer_factor/`, and recomputes every figure and printed number from them. Some values quoted in the section-10 text are not in the committed records: the dR-matched values of 10.1 and 10.2 and the no-veto column of the 10.1 table, the ratios, fractions and significance quoted in 10.3, and in 10.4 the no-veto, MC `\|eta\|` and control-muon fractions, the azimuth significance, the fit probability and the halo-corner fractions. They come from the fold's full output and, for the MC `\|eta\|` fractions, from the `vtxRespMC2_2018` outputs. The fold script that wrote the JSON records is not in the repository. | reads only; writes nothing |

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
skimmed 2018A-D `DoubleMuon` samples on the data side), not 2018C alone. The fold itself is not
in the repository. It ran on their `.coffea` outputs and wrote the three JSON records in
`transfer_factor/`, which the write-up re-reads; the write-up also re-opens the segment-matched
data `.coffea` files directly for its `phi` x `nHits` and `|eta|` x `nHits` figure. Rerunning
the notebooks therefore reproduces every input of the write-up except those three JSON
records.

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
samples when efficiency is counted per DSA muon (6 counted per lepton jet). There the
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