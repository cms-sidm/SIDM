# Transfer-factor records read by section 10 of the study notebook

These three files hold the values section 10 of `../vertexChi2CosmicVetoStudy.ipynb` reads, and nothing else. The fold script that wrote them is not in this repository. Channel keys: `spread_veto` (Spread veto + cosAlpha), `vertex_veto` (vertex veto + cosAlpha) and `no_veto`. Targets: `TAIL` is DSA-muon pT >= 200 GeV with overflow, and `SB` is 100-200 GeV.

- `tf_transfer_factor_seg68.json`: the fold with the segment-matched DSA-PF control. It was computed from the `vtxResp2_2018Full` outputs of `runDataSegMatch2018.ipynb` (68 DoubleMuon 2018 chunks, 59.83/fb) and the `vtxRespMC2_2018` outputs of `runMCRespSegMatch2018.ipynb` (8 MC samples). Its keys:
  - `provenance`: luminosity, data chunks and MC samples used.
  - `transfer`: per channel and slice, the transfer-matrix rows in pT_PF source-bin order, from [0,20) up to > 500 GeV. Only the MC effective numerator count `mc_neff_num` is kept.
  - `closure`: per channel, slicing (`abseta` or `nHits`), slice index and target. It holds the observed and predicted counts in data and MC, the statistical and three systematic terms (`sig_circ`, `sig_bin`, `sig_chan`), and the data and MC closures `R`, `R_mc` (Obs/Pred - 1) with their uncertainties.
  - `inclusive_closure`: the same closure summed over slices.
  - `sideband_transport_failures`: the sideband slices where the prediction exceeds the observation.
- `tf_transfer_factor_full68.json`: the same fold with the dR-matched control. It was computed from the same 68 chunks and 8 MC samples, in the `vtxResp_2018Full` (`runDataFull2018.ipynb`) and `vtxRespMC_2018` (`runMCRespFull.ipynb`) outputs. It keeps only the MC closure of the Spread-veto `nHits < 20` tail (`closure`) and the sideband failures.
- `tf_seg_aux.json`: auxiliary numbers from the same outputs. Its keys:
  - `control_sizes`: control-muon counts and feed-up numerators and denominators for both control definitions.
  - `tiebreak`: the fraction of control muons with several PF partners, and the relative shift of the feed-up rate between partner choices.
  - `halo`: the |dxy_BS| and CSC-only fractions of the mu-type lepton-jet DSA muons, in the bulk and the two tail slices, for data and MC.
