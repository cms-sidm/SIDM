# Proposal: the lepton-jet reconstruction efficiency uncertainty

This is the uncertainty AN Sec. 7 announces as "Lepton Jet reconstruction and selection":
how well the simulation describes the probability that a dark photon gives a lepton jet
(LJ) that is reconstructed and passes the LJ selection. It is expected to be the largest
signal uncertainty, larger than the four in [`README.md`](README.md).

## What is covered today, and what is not

The default corrections include four scale factors for LJ constituents: the PF-muon Loose ID,
the electron reconstruction and the photon pixel veto (POG), and since step 3 the DSA muon,
measured here. All were measured with prompt objects. Our objects are different:

- 79 to 100% of the lepton-jet PF muons are outside the phase space of the muon SF
  (displaced, or with another muon closer than ΔR = 0.3; corrections README).
- The electron ID of the analysis (MVA WPL with missing inner hits) and the custom photon ID
  have no SF at all, and the DSA SF is measured on prompt muons only.
- Nothing describes the LJ clustering itself or the LJ-level cuts (pT, η, the dxy and dz
  spread of the muons).

So the SFs per constituent are the nominal correction, and this study measures what they miss
(step 3).

## Why 2mu2e and 4mu need different values

- **Different objects.** A muon-type LJ is made of PF and DSA muons, reconstructed with
  the tracker and the muon system. An egamma-type LJ is made of GSF electrons with missing
  inner hits and photons (conversions and merged electron pairs), reconstructed in the
  ECAL with the tracker. They fail in different ways and need different control samples.
- **Different multiplicity.** A 4mu event has two muon-type LJs, so the muon-type
  uncertainty enters twice, fully correlated (about twice the per-LJ value). A 2mu2e event
  has one of each, and the two uncertainties come from independent measurements, so they
  add in quadrature.

This gives two nuisances, `CMS_ljreco_mu` and `CMS_ljreco_egm`, each a size per LJ (kind
`per_lj` in `systematics.yaml`) that a 2mu2e event takes once per type and a 4mu event twice for
the muon type.

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

**Version 3 in simulation** (campaign `lj_reco_cr_v3`, 2026-11-24; notebook 07, section 6, and
[`09_dsa_prompt_muons.ipynb`](09_dsa_prompt_muons.ipynb); outputs in
`/store/user/scampbel/sidm_coffea_outputs/lj_reco_cr_v3_results/`): the generator truth of the
probes, a purer Z → μμγ sample, the J/ψ probes against the displacement of their vertex, and the
DSA muons of prompt muons.

- **Closure of the sideband subtraction.** Muon-type: it agrees with the generator truth to 0.02
  where the statistics are large and within the uncertainties from 30 GeV and ΔR = 0.03 up, but not
  below ΔR = 0.03 (0.29 ± 0.25 and 0.83 ± 0.07 against 0.96 and 0.94). Egamma-type: with the
  selection of v2 the sidebands lie 0.11 above the truth on average (2.8 standard deviations
  together), with the FSR selection 0.06 (1.6), every bin within 2 standard deviations;
  unconverted photons agree in both.
- **Production radius.** Above 50 GeV, 97% of the simulated J/ψ come from b-hadron decays, most made
  0.2 to 2 cm from the beam line. The passing fraction is 1.00 within 0.5 mm and 0.98 to 0.99 out
  to 5 cm, so the mix of prompt and non-prompt J/ψ in data matters at 2% at most.
- **Converted photons.** The FSR selection (both muons isolated, m(μμ) + m(μμγ) < 180 GeV) removes
  the QCD probes from the Z peak and keeps 84% of the Drell-Yan ones. With every background the
  fraction of FSR photons for R9 from 0.5 to 0.94 rises from 0.22 to 0.64 to 0.65 to 0.82; about
  680 such photons are expected in data.
- **J/ψ vertex displacement.** The passing fraction is flat, 0.98 to 1.00 out to 5 cm (0.96 ± 0.11
  from 5 to 15 cm).
- **Safety.** Every region of the FSR sample is safe (0.002 signal events in the worst case), and
  every J/ψ vertex bin below 15 cm (at most 0.55 events, binned Z at most 0.40). Above 15 cm, 4Mu at
  800 GeV, 5 GeV, 5 mm puts 4.1 events against 9,300 (binned Z 1.28, a caution for shape uses), so
  that bin is never read. The data campaign is `lj_reco_cr_data_v2`.
- **DSA muons of prompt muons.** The displaced ID has no impact-parameter cut, so prompt muons test
  it. 81% of the prompt muons of Z → μμ have a selected DSA muon (0.85 in the barrel, 0.56 above
  |η| = 2.1); both muons of a boosted J/ψ have one only 44 to 52% of the time above 40 GeV. The
  analysis triggers raise the fraction by 3.3 points (0.14 to 0.23 for pairs), Mu37_TkMu27 by 0.6
  points, IsoMu24 not at all. Our skimmed ntuples keep only events of the analysis triggers, so the
  measurement in data (below) uses Mu37_TkMu27 on the unskimmed DoubleMuon ntuples of 2018C.

**Version 3 in data** (campaign `lj_reco_cr_data_v2`, all 2018 DoubleMuon files, the alone
versions only, read through the safety table of version 3; notebook
[`08_lj_reco_data.ipynb`](08_lj_reco_data.ipynb), section 6; outputs in
`/store/user/scampbel/sidm_coffea_outputs/lj_reco_cr_data_v2_results/`). Every region is read
except the J/ψ vertex bin above 15 cm, which is set to zero in data and in the simulation.

- **Muon-type against the vertex displacement.** Above 50 GeV the passing fraction in data is 0.979
  to 0.985 in every bin of L_xy from 0 to 15 cm. For L_xy above 1 mm, where the simulation is
  precise, data/simulation is **0.995 ± 0.003** (0.989 to 0.998 per bin out to 5 cm, 1.02 ± 0.12
  from 5 to 15 cm). Below 0.5 mm the simulation holds few of its J/ψ (1.6% within 0.1 mm, against
  7.9% of the data) and its subtraction gives passing fractions at or above one, so the ratio of
  0.98 there, which the integrated ratio of version 2 (0.988 ± 0.006) contains, reflects the
  simulation rather than the data.
- **Converted photons** (the FSR selection, R9 < 0.94, above 30 GeV). 0.908 ± 0.048 in data (512
  probes) and 0.966 ± 0.094 in Drell-Yan, a ratio of **0.94 ± 0.10** (0.95 ± 0.11 from 30 to 60
  GeV; 30 data probes above), with the precision set by the simulated statistics. Against R9 the
  ratio is 0.76 ± 0.19, 0.91 ± 0.19 and 0.91 ± 0.17 from 0.5 to 0.94. The sidebands remove about
  half of the peak in data (0.54) as in the simulation (0.58), where the closure test agrees within
  its uncertainties.
- **Unconverted photons** (the FSR selection, R9 ≥ 0.94, above 30 GeV). 0.834 ± 0.013 in data and
  0.816 ± 0.026 in Drell-Yan, a ratio of **1.02 ± 0.04** (1.04 ± 0.04 from 30 to 60 GeV, 0.86 ± 0.09
  above 60 GeV), in line with the selection of version 2.
- **Still not measured:** displaced DSA muons, prompt electron pairs, muon pairs below ΔR = 0.03,
  converted photons above 60 GeV. Their scale factors come from the nearest measured bin (step 3).

**DSA muons of prompt muons in data** (campaign `dsa_prompt_data_v1`, the unskimmed DoubleMuon
ntuples of 2018C, the alone versions only, read through the safety table of notebook 09, which
found both regions safe; notebook [`10_dsa_prompt_data.ipynb`](10_dsa_prompt_data.ipynb); outputs
in `/store/user/scampbel/sidm_coffea_outputs/dsa_prompt_data_v1_results/`).

- **Single muons** (Z → μμ tag and probe, Mu37_TkMu27). 70.4% of the prompt muons have a selected
  DSA muon in data and 81.7% in Drell-Yan, a ratio of **0.862**: 0.87 from 20 to 60 GeV falling to
  0.80 above 120 GeV, 0.87 in the barrel and in the endcaps up to |η| = 2.1, 0.83 in the overlap,
  0.77 above |η| = 2.1. With the analysis triggers it is 0.88, with IsoMu24 0.87, for tight probes
  0.87. The simulation overestimates the DSA reconstruction and selection by about 14%.
- **Two collimated muons** (J/ψ pairs above 40 GeV, analysis triggers). Both have a DSA muon in
  50.3% of the pairs in data and 65.6% in the simulation, a ratio of **0.77 ± 0.01**, the square of
  the single-muon ratio with the same triggers (0.88² = 0.78), so a scale factor per DSA muon also
  describes the pairs from ΔR = 0.03 to 0.4.

### Step 3: corrections per constituent, and the lepton-jet modelling with them applied (2026-11-24)

The advisor's answer (2026-11-24): proceed with the corrections and systematics, taking the value
from the nearest populated bin, with scale factors **per constituent** and their uncertainties;
then evaluate the lepton-jet reconstruction and selection with the corrections applied, and assign
an additional systematic where needed.

- **Scale factors per constituent.** Every PF muon, DSA muon, electron and photon of a selected
  lepton jet carries a scale factor (`sidm/configs/corrections.yaml`, all in the default set): the
  PF-muon Loose ID, the electron reconstruction and the photon pixel veto from the POGs, and the
  DSA muon from this study (`dsa_muon_id`). Their uncertainties are the nuisances `CMS_eff_m_id`,
  `CMS_eff_m_dsa`, `CMS_eff_e_reco` and `CMS_eff_g_pixelveto`.
- **The DSA-muon scale factor** ([`10_dsa_prompt_data.ipynb`](10_dsa_prompt_data.ipynb), section
  6; `sidm/data/corrections/dsa_muon_sf_2018.json`): the Z → μμ tag and probe with Mu37_TkMu27 in
  the bins of muon pT (10 to 200 GeV) and |η| (0, 0.9, 1.2, 2.1, 2.4), evaluated at the DSA muon's
  own pT and |η|. It is 0.80 to 0.89 up to |η| = 2.1 and 0.76 to 0.89 above, falling with pT in
  every region. Bins with fewer than 2,000 probes in data take the value of their merged neighbours
  (below 20 GeV, and above 80 GeV for |η| > 2.1), and a DSA muon above 200 GeV takes the last bin.
  The uncertainty is the statistical one (mostly below 0.01) and the largest change with the
  analysis triggers, IsoMu24 or tight probes: 0.01 to 0.03 per muon up to |η| = 2.1 and 0.05 to
  0.08 above, mostly from the analysis triggers, which select probes more likely to have a DSA muon.
- **The lepton-jet modelling with the corrections applied** (campaign `lj_reco_cr_v4`, every
  background; [`12_lj_corrected_closure.ipynb`](12_lj_corrected_closure.ipynb)). The probes of the
  control samples in simulation, each passing probe weighted by the scale factors of the
  constituents of the lepton jet that makes it pass (or of its DSA muons for the DSA probes), as the
  analysis weights its selected lepton jets; the `_sf` histograms carry this on an `sf` axis with
  the shifted scale factors. The passing fraction with the corrections is the weighted passing
  probes over all the probes unweighted, since a scale factor changes the probability to pass, not
  the number of probes. Data are those already read (no new data).
- **The additional systematic** (`CMS_ljreco_mu`, `CMS_ljreco_egm`, kind `per_lj` in
  [`systematics.yaml`](systematics.yaml)): per lepton-jet type, max(|1 − r|, σ_r) of the
  data/simulation ratio r with the corrections applied, over the measured range above the turn-on,
  or in the measured bin nearest to the signal where that differs by more than two standard
  deviations; the muon-type value is the larger of its PF-muon and DSA-muon regimes. One log-normal
  term per type, which a 2mu2e event takes once per type and a 4mu event twice (muon-type, fully
  correlated).
- **Into the limits.** Campaign `systematics_v4` (the signal with `dsa_muon_id` in the default
  corrections; the simulated background of `systematics_v3`, which the data-driven limits do not
  use), then notebooks 03, 04 and 05.
- **Bins without a measurement.** The control samples cannot reach every corner of the signal
  phase space (the most displaced and most collimated lepton jets, the highest pT). There the
  scale factor is taken from the nearest measured bin (agreed with the advisor, 2026-11-23), and
  every such case is recorded in a coverage map: for each LJ type and each bin of pT,
  displacement and collimation of the signal, the measured bin its scale factor comes from, with
  the fraction of each signal point's lepton jets it covers. The group then assesses the cases one
  by one. The map exists in the generator quantities of step 1
  ([`11_lj_coverage_map.ipynb`](11_lj_coverage_map.ipynb), 2026-11-24): over the limit grid, a
  median 97 to 99.7% of the muon-type lepton jets lie below ΔR = 0.03 (all for 0.25 GeV dark
  photons, 10% for 5 GeV ones), 65 to 68% beyond L_xy = 15 cm and 28 to 29% above 300 GeV, and 97%
  of the egamma-type ones above 120 GeV. The nearest-bin rule therefore applies to most of the
  signal, and the per-point table on EOS lists every case. The DSA scale factor takes its last pT
  bin above 200 GeV; the additional systematic compares the measured regions nearest to the
  signal with the whole range, by the rule above.

**Results of step 3** ([`12_lj_corrected_closure.ipynb`](12_lj_corrected_closure.ipynb), campaign
`lj_reco_cr_v4`, 2026-11-24). Data/simulation with every scale factor applied:

- **Muon-type, J/ψ above 50 GeV:** 1.009 ± 0.006 (0.988 ± 0.006 without), 1.00 to 1.025 in every
  pT bin, 1.017 ± 0.003 at vertex displacements from 1 mm to 15 cm. A quarter of these lepton jets
  hold one PF and one DSA muon (a check on a few simulated QCD files), so the sample also tests the
  DSA factor inside lepton jets, and that factor makes most of the change.
- **Two collimated DSA muons, J/ψ pairs above 40 GeV:** 1.008 ± 0.017 (0.767 ± 0.013 without): the
  factor per DSA muon describes the pairs. On its own measurement (Z → μμ) it gives 0.999.
- **Egamma-type, FSR photons above 30 GeV:** 1.044 ± 0.041 (0.991 ± 0.039 without; the pixel-seed
  veto factor, 0.964 per barrel photon, lowers Drell-Yan by about 5%).
- **The additional systematic:** no region nearest the signal differs from its whole range by more
  than 1.3 standard deviations, so `CMS_ljreco_mu` is 1.7% per muon-type lepton jet (the DSA
  regime, which also covers the displaced J/ψ) and `CMS_ljreco_egm` 4.4% per egamma-type one. They
  raise the expected limits by about 0.5%; the corrections, most of all the DSA factor, by a median
  23% (2Mu2E) and 37% (4Mu).

## What is needed

These were the questions for the group; the answers so far are marked.

1. **Data approval.** Given by the advisor (2026-11-16) for the control samples proposed here,
   each to be thought through before its first look; the group's safety criteria decide which
   regions are read. As with the earlier approvals, the signal region is never read.
2. **SF granularity.** Answered (2026-11-24): per constituent, with the lepton-jet modelling
   checked with the corrections applied and an additional systematic where needed (step 3); the
   nearest populated bin where the control samples do not reach. Done (step 3). Open: the rule
   takes the whole measured range unless a region nearest the signal disagrees with it by more than
   two standard deviations; taking the nearest regions regardless would give up to 8.3% (muon-type)
   and 12% (egamma-type), set by their statistics, instead of 1.7 and 4.4%.
3. **Provisional value.** No longer needed: step 2 measured the prompt regime and step 3 applies
   it.
4. **The sparse-bin rule** (raised by step 2, for the advisor). Whether a region's bins that expect
   under one background event should be judged together as one more region, as done here, or
   left out of the verdict, as in the working-point decks.
5. **DSA muons in data.** Measured on 2018C (notebook 10) and applied per DSA muon in bins of pT
   and |η| (step 3, answered 2026-11-24). Open: whether to add the other periods of 2018
   (unskimmed DoubleMuon ntuples are on EOS for period A and, in blocks, for period D, not yet in
   `data.yaml`), and whether the analysis-trigger difference belongs in its uncertainty (it sets
   most of it).
