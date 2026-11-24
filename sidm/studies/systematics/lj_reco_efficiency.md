# Proposal: the lepton-jet reconstruction efficiency uncertainty

This is the uncertainty AN Sec. 7 announces as "Lepton Jet reconstruction and selection":
how well the simulation describes the probability that a dark photon gives a lepton jet
(LJ) that is reconstructed and passes the LJ selection. It is expected to be the largest
signal uncertainty, larger than the four in [`README.md`](README.md).

## What is covered today, and what is not

The default corrections include three standard scale factors for LJ constituents: the
PF-muon Loose ID, the electron reconstruction, and the photon pixel veto. All three were
measured with prompt, isolated objects. Our objects are different:

- 79 to 100% of the lepton-jet PF muons are outside the phase space of the muon SF
  (displaced, or with another muon closer than ΔR = 0.3; corrections README).
- DSA muons, the electron ID of the analysis (MVA WPL with missing inner hits) and the
  custom photon ID have no SF at all.
- Nothing describes the LJ clustering itself or the LJ-level cuts (pT, η, the dxy and dz
  spread of the muons).

So the provisional SFs stay as the nominal correction, and this study measures what they
miss.

## Why 2mu2e and 4mu need different values

- **Different objects.** A muon-type LJ is made of PF and DSA muons, reconstructed with
  the tracker and the muon system. An egamma-type LJ is made of GSF electrons with missing
  inner hits and photons (conversions and merged electron pairs), reconstructed in the
  ECAL with the tracker. They fail in different ways and need different control samples.
- **Different multiplicity.** A 4mu event has two muon-type LJs, so the muon-type
  uncertainty enters twice, fully correlated (about twice the per-LJ value). A 2mu2e event
  has one of each, and the two uncertainties come from independent measurements, so they
  add in quadrature.

This gives two nuisances, `CMS_ljreco_mu` and `CMS_ljreco_egm`. Each is a per-LJ weight
with an up and a down shift, so the processor turns it into the right per-event value in
each channel, as it does for the existing object SFs.

## Plan

### Step 1: where the efficiency comes from, in simulation (done, 2026-11-14)

For every dark photon decaying to muons (electrons), the efficiency is the fraction with a
selected muon-type (egamma-type) LJ within ΔR < 0.4. It is mapped against L_xy, pT and the
opening angle ΔR(ℓℓ), over the whole signal grid:

- **Reused:** the efficiency chain of [`../ljRecoEffi_studies`](../ljRecoEffi_studies). It
  consists of:
  - the channels `baseNoLj_A_mumu` (denominator) and `baseNoLj_A_mumu_matched_lj`
    (numerator), with the constituent splits `..._matched_pflj_highPT`,
    `..._matched_dsalj_highPT` and `..._matched_pfDsalj_highPT`;
  - the egamma-type `baseNoLj_A_ee_matched_lj` with its electron and photon splits;
  - the `lj_reco_effi_base` histograms;
  - `utilities.get_eff_hist`.
- **Added for this study:**
  - `baseNoLj_A_ee`, the egamma-type denominator, which did not exist;
  - `genA_toE_lxy_vs_daughters_dR`, the electron counterpart of the existing
    `genA_toMu_lxy_vs_daughters_dR`;
  - the collection `lj_reco_effi_systematics`;
  - the campaign config [`condor/campaign_configs/lj_reco_eff_v1.yaml`](../../../condor/campaign_configs/lj_reco_eff_v1.yaml)
    (all 180 signal points, 422 jobs);
  - the notebook [`06_lj_reco_efficiency.ipynb`](06_lj_reco_efficiency.ipynb) and its helpers
    in `systematics_study.py`.
- **Outcome:** efficiency maps for each LJ type. They show which regime each signal point
  depends on (prompt and collimated, displaced in the tracker, beyond the tracker), and so
  which control samples are needed and how to weight them per point.

**Results** (run of 2026-11-14; figures, interpretations and the table per point in the
notebook; outputs in `/store/user/scampbel/sidm_coffea_outputs/lj_reco_eff_v1_results/`):

- **Regimes.** On the limit grid only the shortest lifetime of each mass pair decays mostly
  inside the pixel layers (93 to 97% of its matched dark photons below 10 cm). The other three
  take 41 to 73% of their LJs from the strip tracker (10 to 100 cm), and the muon-type LJs of
  the two longest lifetimes up to 53% from beyond 100 cm. The displaced rows of step 2 matter
  most. `efficiency_regimes.csv` gives the fractions per point, the weights for a provisional
  value and later for the scale factors.
- **Muon-type.** In the 4Mu samples the efficiency is 0.91 below 10 cm, 0.82 from 10 to 100 cm,
  0.62 from 100 to 200 cm and 0.33 from 300 to 500 cm. In the 2Mu2E samples it is 0.97 to 0.99
  out to 200 cm, because the base selection requires two muons above 26 GeV and a 2mu2e signal
  event has no others than the two of its muon dark photon: the dependence on the displacement
  sits in the event acceptance there. The muon-type SF must therefore be measured without that
  condition and applied per LJ in both channels, as planned. It also depends on pT or the
  opening angle, not on L_xy alone.
- **Egamma-type.** 0.89 below 10 cm, 0.80 from 10 to 100 cm, and zero beyond 140 cm, where the
  ECAL starts; L_xy (and pT near the 30 GeV threshold) describes it.
- **Collimated muon pairs.** The heaviest bound states with the lightest dark photon
  (m_Zd = 0.25 GeV at 800 and 1000 GeV, muons about 0.0006 apart in ΔR) have the lowest prompt
  muon-type efficiency (0.87 and 0.77 in 4Mu), and no dimuon resonance reaches that opening
  angle (step 3, "Bins without a measurement").

### Step 2: data and simulation in control samples (done for the prompt regime, 2026-11-17)

The advisor approved reading data in the control samples proposed here, asking that each be
thought through before its first look, especially distributions whose tails may be
signal-dominated. Scoping them against our inputs (the 2018 DoubleMuon SkimNano V2 data and the
simulation already processed) gave:

| LJ type | regime | control sample | with our inputs |
|---|---|---|---|
| egamma-type | prompt photons, converted or not (the displaced egamma-type lepton jets are mostly photons without a pixel seed) | Z → μμγ, the photon from final-state radiation | **feasible**: DoubleMuon data against Drell-Yan simulation; channel `lj_cr_Z_to_mumugamma` |
| muon-type | prompt, collimated | boosted J/ψ → μμ | **feasible**: DoubleMuon data; the simulation has only J/ψ from b-hadron decays (QCD, tt̄); channel `lj_cr_JPsi_to_mumu` |
| muon-type | displaced, DSA | cosmic muons | **not as a scale factor**: there is no cosmic simulation, and in our data both legs are always reconstructed (the triggers need two L2 muons, the skim two muons above 26 GeV), so at large dxy the second leg is not a free probe. A data-only look at the DSA identification against dxy would remain, as a cross-check |
| egamma-type | prompt, collimated | boosted J/ψ → ee | **not possible**: needs electron-triggered (EGamma) data, which are not in our ntuples |

So the displaced muon-type regime, the one that matters most for the long lifetimes, has no
data/simulation measurement with the present inputs (step 3, "Bins without a measurement").

**The probes.** In each channel the lepton jets are built and selected as in step 1; the probes
come from the uncut collections (`objects.zmmg_probe_photons`, `objects.jpsi_probe_dimuons`), so
the probe definition does not depend on the lepton jets, and each histogram entry is a probe
flagged when a selected lepton jet of the tested type lies within ΔR < 0.2:

- Z → μμγ: exactly two selected PF muons of opposite charge with 35 < m(μμ) < 85 GeV, and any
  reconstructed photon with pT > 20 GeV, |η| < 2.4, ΔR > 0.5 from both muons (so the clustering
  keeps it apart from them) and 60 < m(μμγ) < 120 GeV. Histogram `lj_cr_zmmg_probes`: photon pT,
  |η|, R9 (a conversion proxy), pixel seed, m(μμγ), passing.
- J/ψ → μμ: two tracker muons of opposite charge with |η| < 2.4, 2.6 < m < 3.6 GeV, ΔR < 0.4 and
  pT(μμ) > 20 GeV. Only this window is kept, and no other part of the dimuon mass spectrum is
  filled. The signal's dark photon masses (0.25, 1.2 and 5 GeV) are outside it, but mismeasured
  pairs of displaced 5 GeV dark photons reach it (below). Histogram `lj_cr_jpsi_probes`: pair pT,
  ΔR(μμ), |η|, the softer muon's pT, m(μμ), passing.

**What data will read.** The control samples do not exclude events of the signal region by
themselves: a J/ψ event with a second lepton jet could sit in region A. Data therefore read only
the `_alone` versions (channels `lj_cr_Z_to_mumugamma_alone`, `lj_cr_JPsi_to_mumu_alone`,
collection `lj_reco_control_samples_alone`), whose probes have no other lepton jet of the signal
selection's kinds (in `ljs`, two or more muons or none) farther than ΔR = 0.4 away
(`objects._alone`). An event of the signal region holds two such lepton jets, so none of its
probes enters; the condition is the same for passing and failing probes. The versions without
it stay in the simulation, to check that it does not change the passing fraction.

**The safety check before any data.** Each control sample is judged in simulation with the
group's criteria (as in the working-point decks): worst case over the 96 points of the limit
grid at their model cross sections, safe if S < 1 event, or S/(S+B) < 5% and Z < 1; not safe if
S/(S+B) > 10% or Z > 2. It is applied to the whole probe histogram and to each probe-pT bin
(the regions a data step reads), as a count, as a distribution (binned Z), and, stricter than the
decks, to the bins with under one expected background event taken together, because a fine
histogram can hold a signal-dominated tail there (the decks report that signal but leave it out
of the verdict). The stricter rule is used with the user's agreement (2026-11-16) until the
advisor rules on it (question 4); `STRICT_SPARSE_BINS = False` in `systematics_study.py` returns
to the decks' verdicts. Every verdict is repeated with the heaviest
simulated background event removed, and the regions that pass only through S < 1 are judged
again as one. A data step reads only regions that pass everything:
`systematics_study.require_safe` stops it otherwise.

**Results in simulation** (campaign `lj_reco_cr_v2`, notebook
[`07_lj_reco_control_samples.ipynb`](07_lj_reco_control_samples.ipynb); outputs and the safety
table in `/store/user/scampbel/sidm_coffea_outputs/lj_reco_cr_v2_results/`):

- **Safety: every region of both samples is safe**, with the sparse-bin rule and in the stress
  test, so both are read in data in full. Z → μμγ: the worst signal point puts 0.010 events into
  the whole sample. J/ψ → μμ: displaced 5 GeV dark photons with mismeasured pairs reach the window
  at high pT, but the worst point puts 3.0 events against 899 background events above 300 GeV
  (S/(S+B) = 0.33%, binned Z = 0.57) and at most 0.22 events in any bin below 150 GeV.
- **J/ψ → μμ, muon-type:** a peak on a flat background, 98% QCD (J/ψ from b-hadron decays), about
  380,000 probes after the subtraction. The passing fraction is 0.75 ± 0.11 from 30 to 40 GeV (the
  lepton-jet threshold) and 0.94 to 1.00 above; flat near one for ΔR(μμ) from 0.03 to 0.12, 0.87 ±
  0.06 from 0.02 to 0.03 and 0.32 ± 0.25 from 0.01 to 0.02. The pairs cover the 5 GeV dark photons,
  the 1.2 GeV ones in part, not the 0.25 GeV ones.
- **Z → μμγ, egamma-type:** photons that do not come from the Z make a background that rises
  steeply with m(μμγ), so the peak is 85 to 95 GeV with the adjacent sidebands 75 to 85 and 95 to
  105 GeV (wider windows gave efficiencies above one). Above 30 GeV, Drell-Yan photons pass in 0.86 ±
  0.04 of the cases from 30 to 35 GeV and 0.95 to 0.99 from 35 to 50 GeV. The estimate with all
  backgrounds agrees with Drell-Yan alone only for unconverted photons (R9 ≥ 0.94: 0.83 ± 0.06
  against 0.86 ± 0.03); converted photons, the closest analogue of displaced egamma-type lepton
  jets, need a purer probe (for example m(μμ) + m(μμγ) < 180 GeV and the photon nearer a muon),
  that is a new version of the sample.
- **The alone condition** keeps 98% of the J/ψ probes and does not change the passing fraction.

**Results in data** (campaign `lj_reco_cr_data_v1`, all 2018 DoubleMuon files, the alone versions
only, read through the safety table; notebook [`08_lj_reco_data.ipynb`](08_lj_reco_data.ipynb);
outputs in `/store/user/scampbel/sidm_coffea_outputs/lj_reco_cr_data_v1_results/`):

- **Muon-type (J/ψ → μμ):** above 50 GeV the passing fraction is 0.9819 ± 0.0004 in data and
  0.9936 ± 0.0062 in simulation, a ratio of **0.988 ± 0.006** (barrel 0.991 ± 0.009, overlap 0.994
  ± 0.007, endcaps 0.976 ± 0.008; ΔR from 0.05 to 0.12: 0.993 to 0.994 ± 0.003). Near the 30 GeV
  threshold the data turn on more slowly (0.87 ± 0.04 from 40 to 50 GeV), where the spectrum inside
  a bin differs between the prompt and non-prompt J/ψ of the data and the non-prompt ones of the
  simulation. Below ΔR = 0.02 both samples hold few probes.
- **Egamma-type (Z → μμγ, unconverted photons above 30 GeV):** 0.854 ± 0.016 in data, 0.842 ± 0.025
  in Drell-Yan, a ratio of **1.01 ± 0.04** (1.04 ± 0.04 from 30 to 60 GeV, which holds 90% of the
  probes; 0.80 ± 0.08 above 60 GeV, from 112 probes, 2.2 standard deviations).
- **Not measured:** displaced DSA muons, converted photons (a purer probe), prompt electron pairs,
  and the most collimated muon pairs. Their scale factors come from the nearest measured bin, each
  case recorded (step 3, "Bins without a measurement").

### Step 3: from ratios to a correction

- **Correction.** The data/simulation ratios become per-LJ scale factors with an up and a
  down shift. They are written as a correctionlib file under `sidm/` and registered in
  `sidm/configs/corrections.yaml`, as described in "your own correction" in the
  corrections README.
- **Bins without a measurement.** The control samples cannot reach every corner of the signal
  phase space (the most displaced and most collimated lepton jets, the highest pT). There the
  scale factor is taken from the nearest measured bin (agreed with the advisor, 2026-11-23), and
  every such case is recorded in a coverage map: for each LJ type and each bin of pT,
  displacement and collimation of the signal, the measured bin its scale factor comes from, with
  the fraction of each signal point's lepton jets it covers. The group then assesses the cases one
  by one.
- **Code needed.** `Corrections` applies object SFs per constituent today, so a per-LJ SF
  needs a small extension there.
- **Into the limits.** Two entries in [`systematics.yaml`](systematics.yaml),
  `CMS_ljreco_mu: {syst: lj_reco_mu}` and `CMS_ljreco_egm: {syst: lj_reco_egm}`, then a
  rerun of [`03_limits.ipynb`](03_limits.ipynb).

### Until step 2 exists

A provisional value per LJ type could come from the data/simulation differences CMS has
published for the nearest objects (displaced muons, photon conversions). It would be
applied to the fraction of each point's LJs in the regime where it holds, which step 1
provides (`efficiency_regimes.csv`). For muon-type LJs it would be judged on the 4Mu
efficiency, since the 2Mu2E one is conditioned by the base selection. This is a placeholder for
the group to agree on, not a measurement.

## What is needed

These are the questions for the group. Step 1 is done and did not depend on them; steps 2 and 3
do.

1. **Data approval.** Given by the advisor (2026-11-16) for the control samples proposed here,
   each to be thought through before its first look; the group's safety criteria decide which
   regions are read. As with the earlier approvals, the signal region is never read.
2. **SF granularity.** Per LJ or per constituent. This decides the change to `Corrections`.
3. **Provisional value.** Which published data/simulation differences to take per regime until
   step 2 exists.
4. **The sparse-bin rule** (raised by step 2, for the advisor). Whether a region's bins that expect
   under one background event should be judged together as one more region, as done here, or
   left out of the verdict, as in the working-point decks.
