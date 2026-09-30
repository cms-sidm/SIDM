# Limit setting from the ABCD signal-region counts

Turns the ABCD-plane yields in the merged coffea outputs into Combine datacards, runs Combine
over them, and plots the expected limits.

Two flavours of datacard are produced: a **counting** card (one bin, region-A background taken
straight from MC) and an **ABCD** card (four bins, region-A background defined inside Combine
as `bNorm*cNorm/dNorm` so the control regions constrain it in the fit). Each is written twice,
with the signal normalised either to the 1 fb reference or to its own theory cross section —
four sets of limits in total:

| directory | method | signal normalised to | `r` means |
|---|---|---|---|
| `limits/` | counting | 1 fb reference | limit on sigma in fb |
| `limits_abcd/` | ABCD | 1 fb reference | limit on sigma in fb |
| `limits_theory/` | counting | its theory sigma | sigma_limit / sigma_theory |
| `limits_abcd_theory/` | ABCD | its theory sigma | sigma_limit / sigma_theory |

Set the normalisation with `DatacardConfig(use_theory_xs=True)`. The two routes agree to a
median of 0.00% (signal strength scales linearly); the few-percent tail is Combine's default 5%
scan tolerance. `run_combine_limits.py` records which normalisation each card used, read from
the card's own header rather than the directory name.

## Contents

Three modules, three notebooks, and three command-line entry points. Everything is
campaign-, method- and plane-agnostic: what changes between studies is the arguments, not
the code.

| path | what it is |
|---|---|
| `datacard_tools.py` | campaigns, yields, the blinding guard, and the counting / ABCD / combined datacard writers |
| `shape_tools.py` | per-region *shape* extraction, for **both** coffea layouts — the pre-binned `abcd_region` axis and the 6D inclusive one where regions are sliced here |
| `plot_tools.py` | every figure in the study, as reusable functions |
| `01_yields_and_datacards.ipynb` | coffea → shapes JSON → counting and ABCD cards |
| `02_limits_and_comparison.ipynb` | `limits.csv` → the result suite and every comparison |
| `03_closure_studies.ipynb` | per-bin closure, constant-factor tests, selection variants |
| `NOTES.md` | dated running log of what was done and what was found |
| `campaigns/<name>/` | per-campaign `shapes_*.json`, `datacards*/`, `limits*/` |
| `plots/<campaign>/` | figures, each stamped with its campaign name; gitignored |
| `slides/` | the decks (gitignored) |
| `reference_limits/` | drop digitised external dark photon contours here |
| `datacards*/datacards.meta.yaml`, `limits*/limits.meta.yaml` | provenance sidecars |

### Command line

| script | release | what it does |
|---|---|---|
| `sidm/scripts/abcd_datacards.py` | SIDM (coffea) | `extract` coffea → shapes JSON (+ `--datacards`); `datacards` from an existing JSON |
| `sidm/scripts/abcd_workspace.py` | **Combine** (ROOT) | shapes JSON → `RooParametricHist` workspace + shape datacards |
| `sidm/scripts/run_combine_limits.py` | either | runs `combine -M AsymptoticLimits` over a card directory, collects `limits.csv` |
| `sidm/scripts/abcd_plots.py` | SIDM | `results`, `compare`, `closure` |

The split exists because ROOT and `RooParametricHist` live in the Combine release while coffea
does not. Note that **nothing here runs Combine except `run_combine_limits.py`** — the others
build its inputs or plot its outputs.

```bash
# 1. extract + counting/ABCD cards          (SIDM release)
python sidm/scripts/abcd_datacards.py extract --campaign sixd_inclusive_v1 \
    --plane iso_iso,iso_dphi --datacards

# 2. shape workspace + cards                (Combine release)
python sidm/scripts/abcd_workspace.py \
    --shapes  .../shapes_iso_dphi.json \
    --outdir  .../datacards_shape_iso_dphi

# 3. limits
python sidm/scripts/run_combine_limits.py -j 8 \
    --datacards .../datacards_shape_iso_dphi \
    --pattern  'datacard_shape_*.txt' --outdir .../limits_shape_iso_dphi

# 4. figures
python sidm/scripts/abcd_plots.py compare --reference counting \
    counting=.../limits_counting "ABCD"=.../limits_abcd "shape"=.../limits_shape
```

## Campaigns

Two productions are analysed, with outputs side by side under `campaigns/<name>/`:

| campaign | notes |
|---|---|
| `cosmic_veto_v1` | original; its merge left `metadata["is_data"]` empty |
| `golden_hotspot_iso025_v1` | **default**; adds the eta-phi hotspot veto and 0.25 isolation; `is_data` populated |

```python
datacard_tools.use_campaign("golden_hotspot_iso025_v1")
datacard_tools.compare_campaigns("limits_abcd_obs")     # cross-campaign join
```

Each campaign directory holds `sr_yields.pkl`, `datacards*/` and `limits*/`. Both notebooks
take a `CAMPAIGN` variable at the top.

Figures are written to `plots/<campaign>/` so two campaigns cannot overwrite each other, and
**every figure carries its campaign name in the top-right corner** — a plot lifted into a talk
still says which production it came from. Cross-campaign comparisons go to `plots/campaigns/`
and are stamped `old vs new`.

## Two ABCD planes

`muIso_dPhi_eHEM_skimv1_v1` ships a degenerate `abcd_region` axis — regions C and D are
*exactly* zero, because the second isolation moved into the event selection and that cut was
the old region boundary. No ABCD card can be built from it. The regions are instead re-derived
from the 2D `abcd_corr_*_iso_vs_abs_dphi` histogram:

```python
dt.use_campaign("muIso_dPhi_eHEM_skimv1_v1")
bkg = dt.group_backgrounds(dt.collect_dphi_plane_yields(dt.BKG_DIR))
dt.dphi_plane_boundaries(sample_out, channel)   # the cuts actually applied
```

* **The axes have opposite sense** to the iso-vs-iso plane: signal is low isolation but *high*
  |dPhi|, so A is the low-iso/high-dPhi corner. Hence a separate function rather than a flag.
* **The |dPhi| cut snaps to a bin edge.** 32 bins over 0..pi means 2.0 is not representable;
  the effective cut is **1.9635**, and `dphi_plane_boundaries()` reports it.
* **The iso-vs-iso plane does not close** (`A/pred` = 15.8 and 57.2); the iso-vs-dPhi plane does
  to within 15%, with far better balanced control regions. Its limits are 4.4x *weaker* — which
  confirms that the old ABCD "improvement" was an artefact of an unrealistically low prediction.
* Comparison figures: `python sidm/scripts/plot_plane_comparison.py`; deck
  `slides/abcd_plane_comparison.tex`.

## Shape-based ABCD (RooParametricHist)

A second, shape-based estimate lives alongside the counting and `rateParam` ABCD cards: the
control regions supply one free parameter *per bin* and region A's background is a per-bin
formula in them, following the Combine `RooParametricHist` ABCD tutorial. **MC only** —
`shape_tools._reject_data_dir` refuses the data merge outright, on top of the usual blinding
guard.

ROOT and `RooParametricHist` live in the Combine release while coffea lives in this one, so the
chain is in two stages with a JSON handover:

```bash
# 1. extract per-ABCD-region distributions (SIDM release)
python sidm/scripts/export_abcd_shapes.py --observable mljlj_two_bin

# 2. build the workspace and datacards (Combine release)
python sidm/scripts/build_parametric_ws.py \
    --shapes  .../campaigns/<c>/shapes_mljlj_two_bin.json \
    --outdir  .../campaigns/<c>/datacards_shape_mljlj \
    --tf-model global_tf

# 3. limits, through the usual driver
python sidm/scripts/run_combine_limits.py -j 8 \
    --datacards .../datacards_shape_mljlj --pattern 'datacard_shape_*.txt' \
    --outdir    .../limits_shape_mljlj

# 4. figures
python sidm/scripts/plot_abcd_shapes.py --shapes .../shapes_mljlj_two_bin.json
```

* **There is no `m_ljlj` histogram per ABCD region yet.** The `abcd_*` collection carries per-LJ
  kinematics only. `--observable mljlj_two_bin` recovers a genuine *two-bin* `m_ljlj`
  distribution by pairing each SR selection with the `VR_invMass` selection that inverts only
  its mass cut. Asking for the real thing prints what the production would have to add.
* **`--tf-model global_tf` is the default and `per_bin` is the tutorial.** `per_bin`
  (`A_i = B_i*C_i/D_i`) needs events in every bin of C and D; `SR_4mu` region C is empty below
  150 GeV, and at four bins Combine dies with `function value is NAN`. `global_tf`
  (`A_i = B_i * sum C / sum D`) takes the shape from B and the normalisation transfer from the
  integrated ratio, and converges on every card.
* **Validated against the existing cards.** Collapsed to one bin it reproduces the `rateParam`
  four-bin ABCD limit exactly — 0.0200 (`SR_4mu`) and 0.0801 (`SR_2mu2e`).
* **Deck**: `slides/shape_abcd.tex` (16 pages) documents the machinery, the validation and the
  naive look. Rebuild with `pdflatex shape_abcd.tex` twice, after regenerating the figures.
* **The two-bin `m_ljlj` fit buys nothing**, a few percent in either direction, because the
  signal puts 0.01–0.6% of its region-A yield below 150 GeV against 19% of the background. The
  150 GeV cut is already well placed; any gain has to come from shape *above* it.

## Provenance

Every set of datacards and limits carries a `.meta.yaml` sidecar, the same convention the
merged coffea inputs use, so the chain **coffea → datacards → limits** is recoverable from one
file and runs can be sorted by the conditions they were made under.

`datacards*/datacards.meta.yaml` records the full `DatacardConfig`, the ABCD convention, the
per-region background yields, the blinding policy, this repo's commit **and whether the working
tree was dirty**, and the upstream `sidm_commit`/timestamps of the coffea files read.

It also carries a **`selection_cuts`** block with the complete cut definitions — object cuts,
post-lepton-jet object cuts and event cuts — copied verbatim out of the input coffea sidecars
for the two SR selections actually used, so the limits record exactly which selection produced
them. `background_signal_definitions_agree` records whether the background and signal inputs
were made with identical cuts; if they were not, the yields are not comparable and the field
says so (currently `true`).

`limits*/limits.meta.yaml` records the Combine version and options, blinded vs not, how `rMax`
was chosen, how many cards ran and failed, a results summary, and embeds the datacard sidecar
whole.

To compare runs:

```python
import datacard_tools
datacard_tools.index_limit_runs()          # one row per limits*/ directory
datacard_tools.load_limit_metadata("limits_abcd_obs")   # the full record
```

## Campaigns

Two productions are analysed, with outputs side by side under `campaigns/<name>/`:

| campaign | notes |
|---|---|
| `cosmic_veto_v1` | original; its merge left `metadata["is_data"]` empty |
| `golden_hotspot_iso025_v1` | **default**; adds the eta-phi hotspot veto and 0.25 isolation; `is_data` populated |

```python
datacard_tools.use_campaign("golden_hotspot_iso025_v1")
datacard_tools.compare_campaigns("limits_abcd_obs")     # cross-campaign join
```

Each campaign directory holds `sr_yields.pkl`, `datacards*/` and `limits*/`. Both notebooks
take a `CAMPAIGN` variable at the top.

Figures are written to `plots/<campaign>/` so two campaigns cannot overwrite each other, and
**every figure carries its campaign name in the top-right corner** — a plot lifted into a talk
still says which production it came from. Cross-campaign comparisons go to `plots/campaigns/`
and are stamped `old vs new`.

## Provenance

Every set of datacards and limits carries a `.meta.yaml` sidecar, the same convention the
merged coffea inputs use, so the chain **coffea → datacards → limits** is recoverable from one
file and runs can be sorted by the conditions they were made under.

`datacards*/datacards.meta.yaml` records the full `DatacardConfig`, the ABCD convention, the
per-region background yields, the blinding policy, this repo's commit **and whether the working
tree was dirty**, and the upstream `sidm_commit`/timestamps of the coffea files read.

It also carries a **`selection_cuts`** block with the complete cut definitions — object cuts,
post-lepton-jet object cuts and event cuts — copied verbatim out of the input coffea sidecars
for the two SR selections actually used, so the limits record exactly which selection produced
them. `background_signal_definitions_agree` records whether the background and signal inputs
were made with identical cuts; if they were not, the yields are not comparable and the field
says so (currently `true`).

`limits*/limits.meta.yaml` records the Combine version and options, blinded vs not, how `rMax`
was chosen, how many cards ran and failed, a results summary, and embeds the datacard sidecar
whole.

To compare runs:

```python
import datacard_tools
datacard_tools.index_limit_runs()          # one row per limits*/ directory
datacard_tools.load_limit_metadata("limits_abcd_obs")   # the full record
```

## Things worth knowing

* **The SR count needs `flow=True`.** The observable axes are `Regular(100, 0, 700)` and
  overflow at the few-percent level. With flow included, the SR sum reproduces the final row
  of the corresponding cutflow *exactly* (verified in both channels, signal and background);
  without it the yield is low by a few percent. The `abcd_region` axis carries no flow
  content, so nothing is lost there.

* **`r` is a cross section in fb.** `utilities.get_xs` returns a 1 fb reference for signal
  unless `use_signal_xs=True`, and `sidm_processor.postprocess` scales the histograms by
  `lumi * xs`. So the datacards are built at 1 fb and Combine's `r` *is* the limit on sigma in
  fb. The theory cross sections are applied afterwards in `limit_plots.ipynb` — no refit, and
  nothing about the coffea processing changes.

* **One fixed bug in the theory cross sections, and one confirmed oddity.** A duplicate YAML
  key meant `2Mu2E_500GeV_5p0GeV_0p8mm` had no cross section at all (fixed — YAML silently
  keeps the last of a repeated key). The values are non-monotonic in mass — 0.0499, 1.264,
  8.812, 0.3416 fb for 200/500/800/1000 GeV — which was flagged for checking and has since been
  **confirmed correct** (2026-09-22). Every expected-excluded point sitting at 800 GeV is
  therefore a property of the model, not a transcription error.

* **Comparing to the published dark photon limits.** `limit_plots.ipynb` will overlay
  digitised external contours from `reference_limits/`, but **none are shipped** — they are
  published results and must come from a source the group agrees on, with provenance recorded.
  Note also that only 40/120 of our points fall on that plot's canvas (it stops at
  `m_A' = 1 GeV`, so our 1.2 and 5 GeV columns are off the edge), and that those limits
  constrain *direct* dark photon production whereas ours constrain a bound state decaying to
  dark photons. Same plane, different quantity.

* **Every figure names its background estimate.** `plots/by_method/` holds the `Lxy` and
  `m_bound` views for `counting`, for `abcd`, and for the two compared on the same axes, with
  the method in the filename. The older figures under `plots/summary/` are counting-only.

* **How the theory line is drawn.** The cross sections in `cross_sections.yaml` depend only on
  `m_bound` — every `m_ZD` and every lifetime at a given mass shares one value (0.0499, 1.264,
  8.812, 0.3416 fb for 200/500/800/1000 GeV). So on an `Lxy` axis, where `m_bound` is fixed
  within a panel, the theory line is a single **horizontal line**; on an `m_bound` axis the same
  four numbers become a **curve**. No interpolation or fitting — it is those four tabulated
  values joined.

* **`epsilon` on the plots.** Derived as `sqrt(80 / m_ZD / ctau) * 1e-6` with `m_ZD` in GeV and
  `ctau` in mm. Since `ctau` depends on both `m_ZD` and `m_bound`, one `m_ZD` column contains
  points from all four bound state masses at different `epsilon^2`; the panels are kept
  separate because those are different models. The figures are scatters of the simulated
  points, not interpolated contours — the grid has only three `m_ZD` values.

* **The background is still MC.** Every ABCD region rests on 1–7 raw simulated events, which
  is why the counting cards carry a `gmN` nuisance rather than a log-normal. When the
  data-driven prediction exists, substitute it for the MC rate.

* **The ABCD card observes the MC region-A count.** `DatacardConfig.abcd_observation`
  defaults to `"mc"`, so the four-bin card's `observation` in region A is the same number the
  counting card observes while the model predicts `B*C/D`. Setting it to `"prediction"` makes
  the card self-consistent by construction, which cannot show closure tension. Note the blinded
  expected limit is insensitive to this — `--run blind` never reads `observation` — so the
  effect only appears in the `*_obs/` directories, produced with `--unblind`.

* **Unblinded-on-MC runs are closure tests, not limits.** `limits_obs/`, `limits_abcd_obs/` and
  their theory-normalised twins run the fit against the MC signal-region count standing in for
  data. The counting card's observation equals its background by construction, so only the ABCD
  card can show tension. Blinded, ABCD looks 3–5x stronger than counting; observed, it is
  1.1–1.3x *weaker*, and both exclude the same 32/120 points.

* **The MC cannot test ABCD closure.** `B*C/D` and the MC region-A count differ by under
  1 sigma in both channels (0.89 and 0.99), with 45–111% uncertainties. The test has no
  statistical power. In particular, the ABCD cards give limits 3–5x stronger than the counting
  cards only because `B*C/D` happens to land below the MC region-A count — a fluctuation, not
  a gain in sensitivity. Do not quote the improvement. Confronted with the MC count as data the
  ABCD exclusion drops from 42/120 to 32/120, close to the counting result of 28/120.

* **Signal contaminates the control regions.** 19/120 points exceed 10% contamination in some
  control region at the 1 fb reference (32/120 at theory cross sections), concentrated in
  region C and reaching 89%. Not because leakage is large — `S_C/S_A` is only 1–3% for 4Mu —
  but because region C holds just 1.67 background events. This is why the ABCD cards enter
  signal in all four bins rather than only in A.

* **The signal region is blinded by construction, in two tiers.** `datacard_tools` first
  trusts `metadata["is_data"]` where it carries a usable value — the new production populates
  it correctly. Where it does not (an *empty* accumulator counts as no information, never as
  "not data"), it falls back to an allow-list: a sample counts as simulation only if it is a
  known signal point or has a cross section in `cross_sections.yaml`. Anything else has its SR
  omitted and `sr_yield` raises `BlindingError`. The guard fails toward withholding too much,
  never toward leaking, and tier 1 additionally catches a data file *named* like signal.

* **The ABCD cards do the arithmetic inside Combine.** `bNorm`, `cNorm` and `dNorm` are
  unconstrained `rateParam`s on the three control regions and region A's background is defined
  as `bNorm*cNorm/dNorm`. Signal is entered in all four bins, not just A, because it leaks into
  region C at the tens-of-percent level for the high-mass points.

* **`target_lumi_pb` extrapolates to a different dataset.** It scales signal and background
  yields by `target/59830` while leaving relative MC statistical errors alone — more
  luminosity does not create more simulated events. Run 3 placeholders are in
  `run_periods.yaml`, unfilled.

* **The cards are blinded.** `observation` is set to the total background, and
  `run_combine_limits.py` passes `--run blind` by default, so only expected limits are
  meaningful. Use `--unblind` once there is real data to unblind to.
