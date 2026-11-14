# Systematic uncertainties in the signal yield

This study sets the systematic uncertainties of the signal yield that enter the limits
(AN Sec. 7, `systematics.tex`). The background is estimated from data, so none of them
apply to it. It covers the basic sources first: the integrated luminosity, pileup, L1
prefiring and the jet energy. The efficiencies without a standard correction (trigger, DSA
muons, the lepton-jet constituents) come later; they are listed under open items in
[`../corrections/README.md`](../corrections/README.md), which this study builds on.

**Only simulation is used.** The analysis is blinded; the checks that need collision data
are listed under [open items](#open-items) for the group to decide.

**Where things stand:** [`PROGRESS.md`](PROGRESS.md) has the current state, what is running,
the next steps and a dated log of the work. Read it first, and update it with every change.

## In short

| source | nuisance | value | why | notebook |
|---|---|---|---|---|
| integrated luminosity | `lumi_13TeV` | 0.83%, flat | 2018 legacy luminosity, CMS-LUM-20-001 | none, [below](#integrated-luminosity) |
| pileup | `CMS_pileup` | per point, asymmetric: 0.2 to 0.4% for a typical point, at most 1.5% on the limit grid | minimum-bias cross section varied by its 4.6% uncertainty | [`01_pileup.ipynb`](01_pileup.ipynb) |
| L1 prefiring | `CMS_l1_prefiring` | 1%, flat | NanoAOD variation (at most 0.3%) plus the DSA muons the weight leaves out (at most 0.7%) | none, [below](#l1-prefiring) |
| jet energy | `CMS_jet_energy` | per point, symmetric: 0.1 to 0.5% for a typical point, at most 2.2% on the limit grid | isolation scaled by 10%, the energy scale of the few-GeV hadronic part of the matched jets | [`02_jet_energy.ipynb`](02_jet_energy.ipynb) |

* Each is one log-normal (`lnN`) nuisance on the signal only, with the same name in every
  card: fully correlated between 2mu2e and 4mu, across signal points, and across the ABCD
  regions and shape bins.
* "Limit grid" means the 48 points per final state the limits use: $m_\Phi$ = 200 to 1000 GeV
  (the masses with a model cross section), without the shortest lifetime (median
  $L_{xy}$ about 0.2 cm), which the group no longer considers
  (`systematics_study.limit_grid`).
* In quadrature they come to about 1.5% at a typical point and at most about 3%, small
  next to the efficiencies that are still to be measured.
* The pileup and jet energy values are those of region A at the working points, from
  [`03_limits.ipynb`](03_limits.ipynb) on the `systematics_v1` campaign. The two notebooks
  that set the method gave slightly smaller extremes, 1.2% and 1.9%, from the corrections
  campaign's looser categories. Where flat values are wanted, 1.5% for pileup and 2.3% for
  the jet energy cover every point of the limit grid.

## From these numbers to the datacards

They enter the datacards as signal-only log-normal rows. The path is now in place and
runs end to end in [`03_limits.ipynb`](03_limits.ipynb):

1. **The nominal signal yields include the corrections.** Sec. 7 says the pileup,
   prefiring and lepton SFs are applied, and the uncertainties here are uncertainties of
   those corrections. `SidmProcessor` applies them by default (`corrections=True`). On the
   SR signal yields they move 2Mu2E by $-6.0$ to $+0.3$% and 4Mu by $-1.8$ to $+0.2$%
   (corrections README, "Effect on the results"). The earlier limit inputs
   (`sixd_inclusive_v1`) were made outside this repository with the plain genWeight.
2. **The limit selection is in this repository.**
   - `SR_2mu2e` and `SR_4mu` in `selections.yaml`: the `2mu2e` and `4mu` categories with
     the MET filters.
   - The `abcd6d_*` histograms in `hists.py`: the 6D layout that
     `limit_plotting/shape_tools.py` slices into regions A to D. They carry the `syst` axis,
     plus `_iso_up` copies with the isolation scaled by 1.10 for the jet energy.
   - The campaign [`condor/campaign_configs/systematics_v1.yaml`](../../../condor/campaign_configs/systematics_v1.yaml)
     produces them for every signal point and background.
   - Not reproduced: the HEM electron veto and the `mu_veto` of the `sixd_inclusive_v1`
     selection, because neither is defined in this repository.
3. **The shape cards take per-point values.** `abcd_workspace.py --signal-systematics`
   reads, for every point, one value per nuisance and ABCD region: a number, or
   `[down, up]`. It writes one log-normal row per nuisance on the signal. Values may move
   in opposite directions in different regions, as the jet energy moves signal from A to B.
   [`systematics.yaml`](systematics.yaml) lists the nuisances, and
   `systematics_study.abcd_inputs` computes the values from the campaign outputs.

The counting and plain ABCD writers in `datacard_tools.py` still have only the flat
`lumi_unc`, `signal_unc` and `bkg_norm_unc`. They would take the same table if those
methods are kept.

## Integrated luminosity

**Assign 0.83%, flat** (`lumi_13TeV`, $\kappa = 1.0083$), on the signal of every card.

* The CMS Run 2 legacy luminosity (CMS-LUM-20-001) gives
  $59.56 \pm 0.49\ \mathrm{fb^{-1}}$ for 2018 (0.83%), replacing the 2.5% of
  CMS-PAS-LUM-18-002. Sec. 7 already quotes it.
* No notebook: the luminosity scales every simulated yield by the same factor, so there is
  nothing analysis-specific to measure. It is one nuisance for the whole analysis (one year).
* The simulation is normalized to $59.83\ \mathrm{fb^{-1}}$
  (`sidm/configs/run_periods.yaml`), the pre-legacy value, while the 0.83% belongs to the
  legacy value. Moving to $59.56\ \mathrm{fb^{-1}}$ lowers every simulated yield by 0.45%; that
  choice is the group's, and Sec. 8 quotes $59.83\ \mathrm{fb^{-1}}$.

## Pileup

**Assign the up and down change of each point's yield** (`CMS_pileup`, asymmetric
`kappa_down/kappa_up`): 0.2 to 0.4% for a typical point and at most 1.5% on the limit grid
(region A at the working points, `03_limits.ipynb`).
Details in [`01_pileup.ipynb`](01_pileup.ipynb).

* The pileup weights (LUM POG, `Collisions18_UltraLegacy_goldenJSON`) reweight the simulated
  number of true interactions to the 2018 profile at a minimum-bias cross section of
  69.2 mb; the up and down variations use 72.4 and 66.0 mb, its 4.6% uncertainty. They move
  the mean number of interactions of the selected events by about $\pm 1.5$ (from about 32).
* The signal samples were generated with the profile the weights assume: their mean pileup
  weight before any selection is one to within 0.1% at every point, so the weights apply to
  them as they are and do not change their normalization.
* The selection hardly depends on pileup. The yields change by a few tenths of a percent,
  and for most points that is within twice its statistical precision (computed exactly from
  `pileup_nTrueInt`, since the weight depends on nothing else). The largest values, 2 to 3%,
  are at $m_\Phi$ = 100 GeV, which is not on the limit grid; on it only two 4Mu points at
  1000 GeV reach 1.0 and 1.2%.
* The isolation, the cut most exposed to pileup (neutral pileup energy in the matched jet),
  changes the efficiency of region A by a few hundredths of a percent for a typical point,
  since the signal lepton jets are far below the cuts, and by several tenths at a few points
  with few events. Including it leaves the largest values on the limit grid unchanged.
* The values must come from the limit selection (item 2 above). They are in every cutflow
  row of a run with the corrections on, so no extra histogram is needed.

## L1 prefiring

**Assign 1%, flat** (`CMS_l1_prefiring`, $\kappa = 1.01$).

* In 2018 only the muon part of the prefiring is active (the ECAL part is one). The weight
  comes with NanoAOD (`L1PreFiringWeight_Nom`, `_Up`, `_Dn`), the product over the
  reconstructed muons of one minus their prefiring probability, and the variations carry the
  statistical and systematic uncertainty of the muon prefiring maps.
* On its own it lowers the SR signal yield by 0.14 to 0.71% (2Mu2E) and 0.37 to 1.22% (4Mu),
  and its up and down variations change the yield by at most 0.16% (2Mu2E) and 0.29% (4Mu)
  (`tables/sr_correction_nuisances.csv` of the corrections study).
* The weight is computed from loose PF muons only, so DSA muons without a PF match do not
  enter it, and at long lifetimes many lepton-jet muons are such DSA muons. The effect of
  the weight shrinks with the lifetime accordingly, from $-0.66$ to $-0.32$% (2Mu2E) and from
  $-1.15$ to $-0.52$% (4Mu) between the shortest and the longest lifetime (medians over the
  masses). If DSA muons prefire as often as PF muons, which is reasonable since the L1 muon
  is built from the same muon chambers, the difference from the shortest lifetime of the
  same mass point is what the weight misses: up to 0.40% (2Mu2E) and 0.72% (4Mu) on the limit
  grid.
* The variation and the missing DSA muons, added in quadrature, give at most 0.41% (2Mu2E)
  and 0.73% (4Mu) on the limit grid. A flat 1% covers both at every point; per-point values
  would claim more precision than the DSA estimate has.
* No notebook: the weight and its variation come with NanoAOD, the effect is well below 1%,
  and the one analysis-specific question, the DSA muons, is answered from the existing
  per-point table:

```python
import numpy as np, pandas as pd
t = pd.read_csv("sidm/studies/corrections/tables/sr_correction_nuisances.csv")
point = ["final_state", "m_phi_GeV", "m_zd_GeV"]
t["rank"] = t.groupby(point).median_lxy_cm.rank().astype(int)         # 1 = shortest lifetime
t["effect"] = 100 * t.l1_prefiring_shift                               # what the weight does (%)
t["half"] = 100 * (t.l1_prefiring_kappa_up - t.l1_prefiring_kappa_down).abs() / 2
shortest = t[t["rank"] == 1].set_index(point).effect
t["dsa"] = t.effect - shortest.reindex(pd.MultiIndex.from_frame(t[point])).to_numpy()
t["total"] = np.hypot(t.half, t.dsa)
grid = t[(t.m_phi_GeV >= 200) & (t["rank"] >= 2)]
print(grid.groupby("final_state")[["half", "dsa", "total"]].max())
```

## Jet energy

**Assign the loss of each point's yield when the isolation is scaled by 1.10, symmetric**
(`CMS_jet_energy`): 0.1 to 0.5% for a typical point and at most 2.2% on the limit grid
(region A at the working points, `03_limits.ipynb`).
Details in [`02_jet_energy.ipynb`](02_jet_energy.ipynb).

* Jets enter the signal selection only through the lepton-jet isolation,
  $I = (E_{\rm jet}/E_{\rm LJ})(1 - f_{\rm lep})$, with the AK4 CHS jet nearest the lepton
  jet; the `jets` and `bjets` cuts are used only in the $t\bar t$ control regions. The energy
  fractions are fractions of the uncorrected jet energy, so $I$ scales with the jet energy.
* The jet energy corrections (Summer19UL18_V5) are applied in NanoAODv9 as corrections of
  the jet momentum, not as weights. The JER smearing of simulated jets is not applied.
* The standard variations, the total JES uncertainty of the matched jet and the JER smearing,
  move the signal yield by about 0.1% for a typical point and by at most about 1%.
* But the matched jets are mostly the lepton jets themselves: about 92% of their energy is
  leptons. The energy the isolation measures is a hadronic remainder of a few GeV, below the
  lowest calibrated $p_T$ of the JES (10 GeV), where its uncertainty is 5% in the barrel and
  11% at $|\eta| = 2.4$. A scale of 1.10 on the isolation is a conservative envelope for it.
  It is larger than the standard JES and JER variations at every point of the limit grid (and
  at all but one of the other points), so one nuisance covers both.
* With the cuts at 0.3 and 0.5 and most signal lepton jets far below them, the envelope
  stays small: 0.51% (2Mu2E) and 0.26% (4Mu) for a typical point of the limit grid, at most
  1.9% and 0.85%.
* The values must come from the limit selection (item 2 above). The processor varies only
  weights, but the isolation can be scaled in a histogram's fill function, as the
  `lj_isolation_*` histograms of the corrections campaign do. A per-event histogram of the
  two region-A isolations, filled with the nominal and the scaled isolation, would give
  them exactly; the per-lepton-jet sum used in the notebook is exact to first order.
* Whether the isolation should use the uncorrected jet energy (the `raw` variant) is a
  question of the isolation's definition, not a systematic; it is open item 5 of the
  corrections study.

## Notes for Sec. 7 of the AN

Sec. 7 belongs to its owner; these are points it may want to take up.

* **Pileup**: the variation applied is 4.6% of the minimum-bias cross section (69.2 mb), the
  LUM POG recommendation, rather than 5%.
* **L1 prefiring**: the loss reaches 1.2% for 4mu at short lifetimes, so "below 1%" holds for
  2mu2e only. The weight leaves out DSA muons without a PF match; the uncertainty above covers it.
* **Jet energy**: the jet energy corrections are applied to the jet momenta in NanoAODv9,
  not as per-object weights, and the JER smearing is not applied. Jets enter only through the
  lepton-jet isolation.
* **Lepton-jet constituents**: the photon pixel-veto SF is applied by default too, alongside
  the muon ID and electron reconstruction SFs.

## Limits with and without the corrections and systematics

[`03_limits.ipynb`](03_limits.ipynb) sets the expected limits with the ABCD method and a
shape fit in $m_{LJLJ}$, from simulation only. It compares two sets of limits:

- **without corrections or systematics:** signal and background weighted by the genWeight
  alone (the `uncorrected` category of the `syst` axis), with no systematic uncertainty on
  the signal;
- **with corrections and systematics:** the standard corrections applied (the `nominal`
  category) and one log-normal row on the signal for every entry of `systematics.yaml`.

Both keep the statistical uncertainty of the simulated signal and the free background
parameters of the ABCD model. The notebook works in five steps:

1. One pass over the `systematics_v1` outputs gives the region shapes at the working
   points with both event weights, and each signal point's per-region values for every
   entry of `systematics.yaml`.
2. It builds the shape cards of both sets (`abcd_workspace.py` in the Combine release).
3. It runs blinded expected limits on both sets (`run_combine_limits.py`, at a relative
   precision of $10^{-4}$).
4. It draws both sets as Brazil plots with the model cross section, in two grids per final
   state: against the bound state mass, one panel per dark photon mass and $L_{xy}$, and
   against $L_{xy}$, one panel per dark photon mass and bound state mass
   (`plot_tools.brazil_grid`). The set with corrections and systematics is the band (green
   and yellow, black dashed median); the set without them is drawn over it in blue (solid
   median, dotted band edges).
5. It shows the nuisance pulls of one point with corrections and systematics
   (`sidm/scripts/run_combine_pulls.py`: `FitDiagnostics` on the Asimov dataset,
   background-only and signal plus background).

All of its inputs and outputs are on EOS, read and written through xrootd, so it runs at the
LPC and on CoffeaCasa alike (see [Notebooks](#notebooks)).

The simulated background uses `DYJetsToLL` for Drell-Yan, as the earlier limits did. The
campaign also processed `DYJetsToMuMu`, which overlaps it and is therefore left out
(`systematics_study.EXCLUDED_BACKGROUNDS`).

**The notebook follows the registry.** To add a systematic, add an entry to
[`systematics.yaml`](systematics.yaml) and rerun the notebook, then notebooks 04 and 05
([below](#what-each-systematic-does-and-checks-of-the-fit)). An entry is one of:

- `flat`: a fixed relative value;
- `syst`: a correction whose `Up`/`Down` weight categories the campaign stored;
- `variant`: a histogram filled with a varied input, such as `abcd6d_*_iso_up`.

New weight variations or histograms need a campaign rerun first.

**Results (`systematics_v1`, 2026-11-14).** The campaign ran in full: 3,158 jobs, no
failures, all 200 samples merged. Four of its input files hold no events and were skipped.

Each set has 175 expected limits for the 180 signal points (five have no signal in region
A), the 96 points of the limit grid among them. The results:

- **Corrections and systematics together:** the expected limits with them are higher than
  without them by a median 2.0% for 2Mu2E (from 2.6% lower to 4.6% higher) and 1.3% for
  4Mu (0.4 to 2.0% higher). The expected exclusion is the same in both sets: 13 of 48
  points for 2Mu2E and 24 of 48 for 4Mu.
- **What moves them:** the corrections lower the signal in region A by a median 5.2% for
  2Mu2E (1.7 to 6.4%) and 1.2% for 4Mu (0.5 to 1.8%). For 2Mu2E they also lower the
  predicted background (B·C/D from 476 to 460 events), which offsets part of that.
- **Systematics alone:** a run on 2026-11-13 compared the corrected limits with and
  without the systematics. They changed the limits by at most 0.3% (2Mu2E) and not at all
  (4Mu), since signal uncertainties of 1 to 2% enter an expected limit only at second
  order. The change above therefore comes from the corrections.
- **Pulls:** in the Asimov fit every systematic stays at zero with its full pre-fit
  uncertainty, as it must for signal-only nuisances against a background-only dataset.
- **4mu background:** the simulated background agrees with the earlier production
  (region A 3.6 events, from one DY event; B·C/D = 0.003).
- **2mu2e background (caveat):** the simulated background at the working point is QCD,
  with about 1.5 effective events in each of B and C. B·C/D is therefore 460 events,
  against 25 in region A of the simulation; the `sixd_inclusive_v1` simulation gave 102.
  The simulation-only 2Mu2E limits are therefore a test of the machinery, not an estimate
  of the sensitivity, which needs B, C and D from data.

## What each systematic does, and checks of the fit

Two notebooks study the systematics on the cards of `03_limits.ipynb`, for one signal point
per final state (2Mu2E and 4Mu at $m_B$ = 500 GeV, $m_{Z_d}$ = 1.2 GeV, $L_{xy}$ = 30 cm)
unless a plot covers the grid. They read the outputs of 03_limits from EOS and stop if those
are not from the current `systematics.yaml`, so the order after adding a systematic is
03, 04, 05. Their figures and fit outputs go to EOS, under `effects/` and `fit_checks/` in
the output directory of 03_limits.

[`04_systematic_effects.ipynb`](04_systematic_effects.ipynb), what each systematic does:

1. **Up and down templates.** The signal $m_{LJLJ}$ in regions A to D with each
   systematic's variations and their ratio to the nominal, against the statistical
   uncertainty of the nominal template. A flat ratio means one log-normal number per region
   (what the cards use) is enough; a sloped one is a shape effect.
2. **Size across the signal grid.** Maps of each per-point systematic, and of the
   statistical uncertainty of the simulated signal, in region A over $m_B$, $m_{Z_d}$ and
   $L_{xy}$.
3. **Impacts.** The shift of $r$ when each nuisance moves by one post-fit standard deviation
   (`sidm/scripts/run_combine_impacts.py`, `combineTool.py -M Impacts`), next to its pull.
   The fit is to the Asimov dataset with the signal injected at the point's expected limit:
   without a signal, a nuisance on the signal alone has no impact. The fits use minimizer
   tolerances 100 times tighter than Combine's defaults, which leave $r$ uncertain by about
   1% of its own uncertainty, as much as the impact of a 1% signal uncertainty.
4. **Limit breakdown.** The expected limits of the limit grid with each systematic fixed
   (`--freezeParameters`, the same as leaving it out of the cards), all of them fixed, and
   every constrained nuisance fixed ("statistics only", which also fixes the statistical
   uncertainty of the simulated signal), as the change each makes to the limit.

[`05_fit_checks.ipynb`](05_fit_checks.ipynb), checks that the fit behaves with the
systematics in:

5. **Signal injection.** Toys with the signal at the expected limit and at three times it,
   fitted back (`sidm/scripts/run_combine_toys.py`, frequentist toys from the pre-fit model,
   so the cards' data_obs is never used). The pull of $r$ should follow a unit Gaussian and
   each unconstrained nuisance should come back as $0 \pm 1$. The limits keep $r \geq 0$, so
   at the expected limit some toys stop at $r = 0$; the test at three times the expected
   limit is free of that.
6. **Likelihood scans.** $-2\Delta\ln L$ against $r$ with every nuisance and with statistics
   only, on the Asimov dataset without signal and with the signal at the expected limit,
   and against each systematic with $r$ floating, compared with $\theta^2$, what its
   constraint gives by itself (`sidm/scripts/run_combine_scan.py`).
7. **Correlations.** The correlation matrix of every parameter of the signal-plus-background
   fit, the free background parameters included (`sidm/scripts/run_combine_pulls.py
   --expect-signal`).

Each section of both notebooks ends with an interpretation of its results and what they mean for
the analysis.

**Results (`systematics_v1`, four systematics, 2026-11-14).**

- **Templates:** all four systematics act on the normalisation of each region; none changes
  the shape of $m_{LJLJ}$ within a region, so log-normal values per region are an adequate
  model. The jet energy moves a little signal from region A to region B.
- **Sizes:** pileup shows no trend over the grid and is largely statistical (notebook 01); the
  jet energy grows towards heavy dark photons and light bound states; the statistical
  uncertainty of the simulated signal (0.4 to 5.3% in region A, largest at 3 and 300 cm) is
  the largest signal uncertainty at many points.
- **Impacts:** each systematic moves the fitted $r$ by about its own size (prefiring 1%,
  luminosity 0.8%, the others 0.1 to 0.6%), against a statistical uncertainty of $r$ of 60 to
  90% at the expected limit.
- **Limit breakdown:** the systematics do not change the expected limits. They enter at second
  order, about 0.01%, below the four decimals Combine prints; the few points that move by 0.1
  to 0.3% move by the same amount whichever systematic is removed, so they are jumps of the
  limit calculation, not effects of the systematics.
- **Signal injection, 2Mu2E:** unbiased to within a few percent with the right uncertainty
  (mean fitted $r$ 3% low at three times the expected limit, median pull −0.08, 68% half-width
  1.00); at the expected limit 10% of the toys stop at the $r = 0$ boundary.
- **Signal injection, 4Mu:** not validated. 35 to 54% of the toy fits fail, nearly all of them in
  toys with an empty region C (0.003 expected background events) and often D (0.75); the
  converged fits return $r$ 13 to 15% low. The 4mu ABCD prediction is undefined whenever C and D
  are empty (open item 2).
- **Scans and correlations:** the systematics neither widen the likelihood nor correlate with
  anything; $r$ correlates with the background of the control regions where the signal leaks
  into them, fully so in 4mu region C, which holds no background.

## Lepton-jet reconstruction efficiency

The next and probably largest signal uncertainty. A plan for measuring it, separately for
muon-type and egamma-type lepton jets, is in
[`lj_reco_efficiency.md`](lj_reco_efficiency.md). Its first step, in simulation only, is done:
[`06_lj_reco_efficiency.ipynb`](06_lj_reco_efficiency.ipynb) maps the efficiency against
$L_{xy}$, $p_T$ and the opening angle over all 180 signal points (campaign `lj_reco_eff_v1`),
and gives per point the fraction of its lepton jets from each decay-radius regime. Steps 2 and
3 wait on the group (data approval, scale-factor granularity, provisional value, and a new
question about the most collimated muon pairs).

### Results (2026-11-14)

- **Regimes.** On the limit grid only the shortest lifetime of each mass pair decays mostly
  inside the pixel layers; the other three take 41 to 73% of their lepton jets from the strip
  tracker (10 to 100 cm), and the muon-type ones of the two longest lifetimes up to 53% from
  beyond 100 cm. The displaced control samples matter most.
- **Muon-type.** 0.91 below 10 cm, 0.82 from 10 to 100 cm and 0.33 from 300 to 500 cm in the 4Mu
  samples. In the 2Mu2E samples it is 0.97 to 0.99 out to 200 cm, because the base selection
  already requires both muons of the muon dark photon (a 2mu2e signal event has no other
  muons): there the dependence on the displacement sits in the event acceptance, so the scale
  factor has to be measured without that condition and applied per lepton jet in both channels.
- **Egamma-type.** 0.89 below 10 cm, 0.80 from 10 to 100 cm, zero beyond 140 cm, where the ECAL
  starts.
- **Collimated muon pairs.** $m_{Z_d}$ = 0.25 GeV at 800 and 1000 GeV (muons about 0.0006 apart
  in $\Delta R$) has the lowest prompt muon-type efficiency, 0.87 and 0.77 in 4Mu, and no
  dimuon resonance reaches that opening angle.

## Open items

1. The 2mu2e simulated background in regions B and C (about 1.5 effective events each,
   from heavily weighted QCD events). With B, C and D from data, as in the approved
   expected limits, this simulation is not needed for the background.
2. The 4mu ABCD prediction at this working point: region C expects almost no background and D
   less than one event, so the transfer factor is undefined in most toys and the fits of the
   signal-injection test fail or come out biased (`05_fit_checks.ipynb`). The working point or
   method for 4mu is a group decision; rerun notebook 05 on whatever is chosen.
3. The HEM electron veto and the `mu_veto` of the earlier limit selection: define them in
   `cuts.py` and add them to `SR_2mu2e`/`SR_4mu`, if the group keeps them.
4. Pileup in data: the number of reconstructed vertices in a control region, data against
   the reweighted simulation, as Sec. 7 plans. Needs collision data, so it waits for the
   group.
5. Prefiring of DSA muons: computing their prefiring probability from the muon prefiring
   maps would replace the estimate above.
6. The JER smearing: apply it to the matched jets in the nominal simulation or keep it as an
   uncertainty, together with the isolation definition (corrections open item 5).
7. The lepton-jet reconstruction efficiency: step 1 is done (`06_lj_reco_efficiency.ipynb`);
   steps 2 and 3 wait on the four questions for the group in
   [`lj_reco_efficiency.md`](lj_reco_efficiency.md).

## Notebooks

| notebook | what it shows | input | runtime |
|---|---|---|---|
| [`01_pileup.ipynb`](01_pileup.ipynb) | the profile check, the yield change under the pileup variations with its statistical precision, and the effect of the isolation | `corrections_v1` | 2 min |
| [`02_jet_energy.ipynb`](02_jet_energy.ipynb) | the matched jets, which JES uncertainty applies, the standard JES and JER variations, and the envelope | `corrections_v1` | 2 min |
| [`03_limits.ipynb`](03_limits.ipynb) | expected limits (ABCD, $m_{LJLJ}$ shape fit) with and without the corrections and every systematic of `systematics.yaml`, as Brazil grids against bound state mass and $L_{xy}$, and the nuisance pulls of one point | `systematics_v1` | about 10 min |
| [`04_systematic_effects.ipynb`](04_systematic_effects.ipynb) | what each systematic does: up and down templates, sizes across the grid, impacts, limit breakdown | outputs of 03 | about 12 min |
| [`05_fit_checks.ipynb`](05_fit_checks.ipynb) | checks of the fit: signal injection with toys, likelihood scans, correlations | outputs of 03 | about 5 min |
| [`06_lj_reco_efficiency.ipynb`](06_lj_reco_efficiency.ipynb) | step 1 of the lepton-jet reconstruction efficiency: the efficiency against $L_{xy}$, $p_T$ and the opening angle, its maps, and the decay-radius regimes of every signal point | `lj_reco_eff_v1` | about 1 min |

Notebooks 01 and 02 read the merged outputs of the corrections campaign
(`corrections_study.CAMPAIGN_DIR`) through `corrections_study.load`, which copies them once
into `sidm/studies/corrections/outputs/` (gitignored).

Notebook 03 keeps everything on EOS:
- It reads the `systematics_v1` merged directory and writes its shapes, systematics, cards,
  workspaces, limits, pulls and logs next to it, in `.../systematics_v1_limits/`.
- Both go through xrootd (`systematics_study.eos_url`, `list_coffea` and `upload`), so it
  runs the same at the LPC and on CoffeaCasa, where EOS is not mounted.
- Combine works in a temporary directory that is removed at the end. Combine itself is
  taken from the `PATH` if it is there, otherwise from the CMSSW release in
  `COMBINE_CMSSW_BASE` or the LPC default of `run_combine_limits.py`. Where neither exists,
  the notebook stops after writing the shapes and systematics to EOS.

Notebooks 04 and 05 do the same: they read the cards, shapes, systematics and limits of
notebook 03 from EOS (`systematics_study.limits_outputs` checks that these come from the
current `systematics.yaml`), copy the cards into a temporary directory for Combine
(`download`), and write their figures and fit outputs to `effects/` and `fit_checks/` in
`.../systematics_v1_limits/`. Every Combine fit goes through a script in `sidm/scripts`,
which sets up the Combine release itself (`systematics_study.run_script`).

Use the `CMSSW_17_0_0_pre4 (SIDM)` kernel at the LPC. The shared helpers are in
[`systematics_study.py`](systematics_study.py).
