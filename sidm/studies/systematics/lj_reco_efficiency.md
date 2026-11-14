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
  angle (question 4 below).

### Step 2: data and simulation in control samples (needs the group's approval for each data sample)

| LJ type | regime | control sample | measured, data against simulation |
|---|---|---|---|
| muon-type | prompt, collimated | boosted J/ψ → μμ (and Υ), whose muons are within ΔR ≈ 2m/pT | probability that a dimuon found with tracker tracks gives a selected muon-type LJ, against ΔR(μμ) and pT |
| muon-type | displaced, DSA | cosmic muons, one leg tagging the other | DSA reconstruction and ID efficiency against dxy |
| egamma-type | displaced, collimated | photon conversions in Z → μμγ, the photon from final-state radiation | probability that a converted photon gives a selected egamma-type LJ, against conversion radius and pT |
| egamma-type | prompt, collimated | boosted J/ψ → ee | as for J/ψ → μμ |

The simulated side uses the backgrounds already processed (DY, QCD, tt̄); the J/ψ rows would
also need simulated boosted J/ψ events.

### Step 3: from ratios to a correction

- **Correction.** The data/simulation ratios become per-LJ scale factors with an up and a
  down shift. They are written as a correctionlib file under `sidm/` and registered in
  `sidm/configs/corrections.yaml`, as described in "your own correction" in the
  corrections README.
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

1. **Data approval.** Which control samples may be read: cosmics, Z → μμγ, J/ψ → μμ and
   J/ψ → ee. As with the earlier approvals, the signal region is never read.
2. **SF granularity.** Per LJ or per constituent. This decides the change to `Corrections`.
3. **Provisional value.** Which published data/simulation differences to take per regime until
   step 2 exists.
4. **Collimated muon pairs** (raised by step 1). How to set the uncertainty of the muon-type LJs
   from the most collimated pairs (m_Zd = 0.25 GeV at 800 and 1000 GeV, 8 points of the limit
   grid per final state), which no control sample of step 2 reaches.
