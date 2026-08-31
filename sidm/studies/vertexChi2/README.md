# Dimuon-vertex-fit cosmic veto, 2018C

Compares a per-lepton-jet dimuon-vertex-fit cosmic veto (`mu_ljs.vtx_chi2`) against the
`*Spread_*` veto used in the 1-LJ control region, on 2018C `DoubleMuon` data, the 2018
background stack and the v10 signal grid, and characterizes the >= 200 GeV DSA-muon tail
that survives either veto.

**Start here: [`vertexChi2CosmicVetoStudy.ipynb`](vertexChi2CosmicVetoStudy.ipynb).** It is the
write-up and reads only the `.coffea` outputs of the run notebooks below; it runs no
processor of its own.

## Two cosmic-enriched regions, and why they are never interchangeable

The study uses two different cosmic-enriched selections. They are named apart everywhere — in
the write-up prose, in the figure legends and here — because they are not the same region and
the difference decides what one of the results is conditional on.

| name used throughout | channel | mu-LJ cuts | uses `vtx_chi2`? |
|---|---|---|---|
| **spread-inverted cosmic region** | `cosmic_muons` | `lj_iso < 0.2`, `inverse_spread`, cosAlpha event cut inverted | no |
| **vertex-inverted cosmic region** | `cosmic_muons_vtxChi2` | `lj_iso < 0.2`, `no vtx or vtx_chi2 >= 5`, cosAlpha event cut inverted | yes |

The spread-inverted region is what makes the section-2b comparison a test of `vtx_chi2` rather
than a tautology: the variable is asked to flag objects selected with entirely different
quantities. The vertex-inverted region is where the azimuthal cosmic template of section 7 is
measured, so the 95% CL residual-cosmic bound is conditional on it — the write-up shows the
bound moves only from 6.2% to 7.8% even if the template's vertical fraction is degraded from
0.99 to 0.90.

## Notebooks

Every run notebook opens its own Dask cluster on LPC Condor, writes one `.coffea` per sample
to `<repo>/RunOutputFiles/<run_tag>/` together with a `.meta.yaml` provenance sidecar, and
copies both to `root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/coffea_outputs/$USER/<run_tag>/`.
`sidm_path` is derived from the working directory, so run them from this directory.
Per-sample outputs already on disk are skipped, so a run can be resumed.

| notebook | what it runs | `run_tag` / output directory |
|---|---|---|
| `runDataBg2018.ipynb` | 2018C `DoubleMuon` plus the background stack (TTJets, DYJetsToLL M10to50 and M50, DYJetsToMuMu M10to50, eleven QCD `Pt` bins, WW/WZ/ZZ) through six channels: Before, the `vxySpread` and `dxySpread` variants of the Spread veto with cosAlpha, the vertex veto with and without cosAlpha, and the vertex-inverted cosmic region. Hist collection `cosmic_veto`. | `vtxChi2_vs_spread_2018C` (19 samples) |
| `runDataCosmic2018.ipynb` | 2018C `DoubleMuon` in the `cosmic_muons` channel only — the **spread-inverted** cosmic region, defined by inverted spread cuts and inverted cosAlpha, which uses no vertex information. Hist collection `cosmic_veto`. | `vtxCosmic_2018C` (1 sample) |
| `runDataTail2018.ipynb` | 2018C `DoubleMuon` through Before, Spread + cosAlpha, vertex + cosAlpha and `cosmic_muons_vtxChi2` (the **vertex-inverted** cosmic region), with hist collections `cosmic_veto` and `dsa_pf_response`. Adds the 2D DSA-muon pT versus nHits, \|eta\|, phi, \|dxy\| and sigma(pT)/pT histograms, and the matched DSA-to-PF response. | `vtxTail2_2018C` (1 sample) |
| `runMCTail2018.ipynb` | The MC counterpart of the tail check, same two hist collections, Before and Spread + cosAlpha only: the eight samples that carry the >= 200 GeV DSA tail (DYJetsToLL_M50 and QCD_Pt80To120 through QCD_Pt800To1000). | `vtxTailMC2_2018` (8 samples) |
| `runSignalVtx.ipynb` | The full v10 signal grid, 180 samples (90 4Mu, 90 2Mu2E), through `base_ljDisplacementIso`, `+ cosAlpha`, `+ Spread veto` and `+ vertex veto`. Hist collection `cosmic_veto` as it stood before the 2D nHits and \|eta\| histograms were added, so those are absent from these outputs. | `vtxChi2_signal_2018` (180 samples) |
| `runSignalVtxNhits.ipynb` | A 26-sample subset of the same grid (13 4Mu, 13 2Mu2E, mirrored mass and lifetime points), re-run after `mu_lj_dsaMuon_pt_nHits` and `mu_lj_dsaMuon_pt_abseta` joined the `cosmic_veto` collection. Used only for the signal cost of an nHits cut. | `vtxChi2_signal_2018_nhits` (26 samples) |
| `vertexChi2CosmicVetoStudy.ipynb` | The write-up. Loads all six output directories above and produces every number and figure of the study. | reads only; writes nothing |

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
sidecars, written by `runlogs/backfill_meta.py` with no reprocessing. Those are labelled as
such: they contain a `backfill:` block naming where each field came from, and their
`sidm_commit` is `null`, because the commit that produced them was never recorded and is not
recoverable. Two of the eight output directories, `vtxTail_2018C` and `vtxTailMC_2018`, are
superseded first passes that no current run notebook produces; their sidecars say so and
derive the hist-collection list from the outputs themselves.

## Run order for a full reproduction from the samples

1. `runDataBg2018.ipynb` — the control-region yields and the `vtx_chi2` threshold scan.
2. `runDataCosmic2018.ipynb` — the spread-inverted cosmic region used as a `vtx_chi2` cross-check.
3. `runDataTail2018.ipynb` — the tail histograms and the DSA-PF response, data.
4. `runMCTail2018.ipynb` — the same histograms in MC. Must use the same channel list as 3.
5. `runSignalVtx.ipynb` — the 180-sample signal efficiency of both vetoes.
6. `runSignalVtxNhits.ipynb` — the 26-sample nHits subset. Needs the extended `cosmic_veto`
   collection, so it must be run from a checkout that has it.
7. `vertexChi2CosmicVetoStudy.ipynb` — the write-up, once 1 to 6 are on disk.

Steps 1 to 6 are independent of each other and can be run in any order or in parallel;
step 7 needs all of them. Signal samples read `../../configs/ntuples/signal_4mu_v10.yaml`
and `signal_2mu2e_v10.yaml`; `runSignalVtx.ipynb` and `runSignalVtxNhits.ipynb` also read a
JSON veto list of corrupt input files and append to it any file that fails to read (17 files
across 15 samples after both runs).

Steps 1, 3, 4 and 5 take hours each on a hundred Condor workers. All of them require a valid
VOMS proxy readable by the Condor schedd and `replace_xcache=True` in `make_fileset`, since
the location YAMLs carry `root://xcache//` URLs that only resolve on coffea-casa.
The proxy path (`X509_USER_PROXY`) and the veto-list path (`SKIP_JSON`) in the run notebooks
are the author's; whoever reruns them sets their own.

## A note on the two MC sets

Two different MC selections appear in the study, and the distinction matters when reading it,
because several headline numbers change completely between them:

- `vtxChi2_vs_spread_2018C` carries the full background stack and is what the control-region
  yields and the threshold scan use.
- `vtxTailMC2_2018` carries only the eight samples that populate the >= 200 GeV DSA tail, but
  is the only MC set with the 2D histograms.

Numbers that must always carry the name of the set they came from:

| quantity | full stack | eight-sample tail set |
|---|---|---|
| pre-veto MC tail (`pT >= 200` GeV) | 293.6 | 133.7 (46% of the stack's) |
| post-veto MC tail | 133.7 | 129.5 (97% of the stack's) |
| fraction of the MC tail the veto removes | **55%** | **3%** |
| pre-veto tail data/MC | 8.22 | 18.06 |
| post-veto tail data/MC | 2.31 | 2.385 +/- 0.426 |
| MC yield over all pT, Spread channel | 108137.7 | 57091.1 (52.8% of the stack's) |

The two sets agree after the veto because the tail set covers 97% of the MC tail there; they
disagree before it because the full stack's pre-veto tail contains one high-weight
`QCD_Pt30To50` event worth about 156 weighted entries, which the veto removes. (That the
eight-sample set's pre-veto tail, 133.7, equals the full stack's post-veto tail, also 133.7,
is a coincidence of two unrelated quantities.)

The tail set does **not** cover the bulk: over all pT it carries 52.8% of the full stack, the
missing half being the soft-QCD bins. Absolute data/MC in the bulk is therefore meaningless
with it, and so is anything built from its total yield — only shapes across bins are used
there. Section 8 of the write-up marks the one place where that restriction bites, and gives
the argument in a form (a ratio of rates) that does not depend on it.

## One more normalization caveat

The control region itself is not normalized, and no scale factor is applied anywhere in this
study. In the Before channel it holds 0.770 data mu-type lepton jets per MC one but only 0.556
data DSA muons per MC one; the ratio of those two is a DSA-muon multiplicity mismodelling, MC
carrying 0.351 DSA muons per mu-type lepton jet against 0.254 in data, 38% more. Veto
efficiencies and the azimuthal cosmic test are ratios within one sample and are untouched by
it; absolute tail data/MC ratios are not, which is why the write-up's load-bearing tail
statement is a ratio of rates rather than the raw 2.4 +/- 0.4 of the eight-sample tail set.
Section 3.1 of the write-up sets this out in full.

## The dimuon-vertex-fit variable and its cuts

`mu_ljs.vtx_chi2` is the best (smallest) `normChi2` of any stored dimuon vertex built from two of the LJ's own muons, taken from the LLPNanoAOD `DSAMuonVertex`, `PatMuonVertex` and `PatDSAMuonVertex` tables. The matching lives in [lj_vertex_chi2.py](../../tools/lj_vertex_chi2.py) and `sidm_processor` fills the field for every LJ, so there is nothing to switch on.

A value of -1 means no such vertex was stored, and it *fails* the keep-cuts by design: the producer only writes pairs whose fit converged, so an LJ whose muons never form a vertex is the signature of a fake or cosmic pairing. The variable is only defined for LJs with at least two muons, so pair any cut or plot with `Mu >= 2`.

One asymmetry to keep in mind when comparing against the `*Spread_*` cuts: the Spread variables are a max over the LJ's muon pairs (the worst pair must be compact), while `vtx_chi2` is a min (one good dimuon suffices). The two only coincide for exactly-2-muon LJs; for >= 3-muon LJs the vertex cut is the looser criterion by construction -- it asks "does this LJ contain a real dimuon", not "is every pairing consistent". The `mu_lj_muonN` histogram quantifies how much of a given region is >= 3-muon.

The cuts in `cuts.py` are `0 <= vtx_chi2 < 5` and its complement at the mu-LJ-cut level, `no vtx or vtx_chi2 >= 5`. Three channels use them, mirroring the vxySpread-based cosmic-veto control regions (the ones on the August 2026 slides; a dxySpread variant of the baseline was also run, `cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu_dxy`). These channels, and every other channel of this study, are defined in [`sidm/configs/selections_vertexChi2_as_run.yaml`](../../configs/selections_vertexChi2_as_run.yaml), not in `selections.yaml`:

- `cr_1egmOr1mu_lj_vtxChi2_CosMu`: `cr_1egmOr1mu_lj_spreadMu_spreadPfMu_CosMu` with the four Spread cuts replaced by `0 <= vtx_chi2 < 5`
- `cr_1egmOr1mu_lj_vtxChi2`: the vertex cut without the cosAlpha event cut
- `cosmic_muons_vtxChi2`: the **vertex-inverted** cosmic region -- `no vtx or vtx_chi2 >= 5` with `lj_iso < 0.2` and the cosAlpha cut inverted. It mirrors `cosmic_muons`, the **spread-inverted** cosmic region, which inverts the Spread cuts instead and so uses no vertex information anywhere in its definition. The two are different selections, they are named apart wherever either is used, and a template measured in one does not transfer to the other.

The distributions are `mu_lj_vtx_chi2` (linear axis; the -1 entries sit in their own bin at the left edge, and note that `PatMuonVertex` normChi2 spans ~8 decades so a large fraction of PF-PF entries lands in the overflow) and `mu_lj_vtx_chi2_log` (log axis, for threshold scans), both in the `cosmic_veto` hist collection. The threshold is not yet optimized: choose it from the log-axis shape together with the signal efficiency of the same channels on the v10 signal samples.

`vtx_chi2` is all-or-nothing: unless ALL THREE vertex tables are present in the input, it is NaN and every vtx_chi2 cut (keep and inverse alike) selects nothing -- a partial table set would otherwise dump whole LJ categories into the inverse-cut region as fake "no vertex" entries, and the processor prints `vtx_chi2 disabled: missing vertex table(s) ...` in that case. An exactly-empty vtxChi2 channel therefore means missing tables; check the worker logs. The samples in `sidm/configs/ntuples/` (`signal_2mu2e_v10.yaml`, `signal_4mu_v10.yaml`, `data_skimmed.yaml`, `backgrounds.yaml`) all carry them. The brute-force closure test lives at [tests/test_lj_vertex_chi2.py](../../../tests/test_lj_vertex_chi2.py).

## Result of the 2018C comparison (August 2026)

Inputs are 2018C `DoubleMuon`, the 2018 background stack (MC scaled by 6.97/59.83) and the
full 180-sample v10 signal grid.

The vertex veto is an equivalent cosmic tagger and a worse cut, so the Spread veto with the
cosAlpha event cut stays as the selection. `0 <= vtx_chi2 < 5` flags 98% of the mu-type lepton
jets in the spread-inverted cosmic region, which is defined without any vertex information, but
it has no data-versus-MC separating power (the data/MC ratio of the lepton jets it keeps is
flat at 0.75 across four decades of threshold, against 0.77 with no cut), it leaves the same
residual data/MC above 200 GeV in DSA-muon pT (2.38 against 2.31 for the Spread veto, both on
the full background stack), and over the signal grid it keeps only 0.69 of the signal lepton
jets against 0.96 - about 30% more signal cost. The one corner where the ordering reverses is
m_Zd = 5 GeV, mediator mass 100-150 GeV, ctau >= 130 mm, exactly 8 of the 180 samples per DSA
muon, where the max-over-pairs Spread variables start failing on genuinely wide displaced pairs
and the min-over-pairs vertex requirement does not.

The second half of the study characterizes what survives either veto. The Spread veto removes
87% of the >= 200 GeV DSA-muon tail in data (2414 -> 309). In MC it removes 55% of the full
background stack's tail (293.6 -> 133.7) and 3% of the eight-sample tail set's
(133.7 -> 129.5) -- two different MC sets, and neither number travels without its label. The
removed component is unambiguously cosmic (pre-veto phi vertical/horizontal ratio
12.2 +/- 1.3 against an MC expectation of 1.3 +/- 0.2, and 93% barrel against 71%; two
independent routes put the pre-veto tail's cosmic content at 87% and 84%), and what survives is
not (post-veto ratio 0.94 +/- 0.15, bounding a residual cosmic contribution below 6% of the
tail at 95% CL, a bound that moves only to 7.8% if the template's vertical fraction is degraded
from 0.99 to 0.90).

The surviving tail exceeds MC by 2.4 +/- 0.4 (eight-sample tail set), but that ratio is not the
load-bearing number: the control region is not normalized, MC carrying 38% more DSA muons per
mu-type lepton jet than data, and no scale factor is applied. The statement that carries is a
ratio of rates, which is immune to that mismatch: an independent DSA-to-PF matched control
gives a data/MC feed-up ratio of 3.47 +/- 0.25, and the post-veto tail rate divided by that
feed-up rate is 0.95 +/- 0.06 in data and 0.99 +/- 0.18 in MC, both consistent with one,
against 4.81 +/- 0.13 for data before the veto. That is a consistency statement of about
factor-1.5 precision and not a demonstration that the two agree, since the control rate is
measured before the veto on a different population, is not flat in |eta| or nHits, and the
closure moves between 0.70 and 1.41 in data (0.57 and 1.35 in MC) as the pT_PF window moves.
The recommendation is a data-driven correction sized on that control rather than a further
object cut: `nHits >= 20` halves the residual but costs about 15% of the signal DSA muons at
all pT (8% of those already in the tail, a number not to quote on its own) and is itself a
variable data and MC disagree about at the 15-sigma level in the bulk; `sigma(pT)/pT < 0.3`
moves the ratio the wrong way while costing 45% of the signal at 1 TeV. The vertex-veto code,
its channels and its histograms stay in the repository as a documented negative result.
