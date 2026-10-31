# Transfer-factor fold of section 10, and the records the study notebook reads

Section 10 of `../vertexChi2CosmicVetoStudy.ipynb` reads three JSON records in this directory. They are extracts of the full outputs of the fold scripts here, and keep only the values the notebook reads. Channel keys: `spread_veto` (Spread veto + cosAlpha), `vertex_veto` (vertex veto + cosAlpha) and `no_veto`. Targets: `TAIL` is DSA-muon pT >= 200 GeV with overflow, and `SB` is 100-200 GeV.

## Scripts

| script | what it does | inputs (run tag under `RunOutputFiles/`) | full output | log |
|---|---|---|---|---|
| `tf_transfer_factor.py` | the fold with the dR-matched DSA-PF control | `vtxResp_2018Full` (`runDataFull2018.ipynb`), `vtxRespMC_2018` (`runMCRespFull.ipynb`) | `full/tf_transfer_factor_full68.json` | `logs/tf_transfer_factor.log` |
| `tf_transfer_factor_seg.py` | the same fold with the segment-matched control | `vtxResp2_2018Full` (`runDataSegMatch2018.ipynb`), `vtxRespMC2_2018` (`runMCRespSegMatch2018.ipynb`) | `full/tf_transfer_factor_seg68.json` | `logs/tf_transfer_factor_seg.log` |
| `tf_seg_aux.py` | control-sample sizes for both control definitions, the segment-match tie-break, and the beam-halo diagnostics of section 10.4 | `vtxResp2_2018Full`, `vtxRespMC2_2018` | `full/tf_seg_aux.json` | `logs/tf_seg_aux.log` |
| `extract_records.py` | writes the three committed records from the three full outputs | `full/*.json` | the three records below | none |

The two fold scripts also read the sample list `sidm/configs/ntuples/data_skimmed.yaml` for the provenance header. Both use 2000 toys per closure entry, 500 per scan point and seed 20260901 (their defaults). They need `numpy`, `scipy`, `yaml` and `coffea`; the notebook venv has them.

Run from the repository root, with `RunOutputFiles/` holding the four run tags (a symlink is fine) and the worktree first on `PYTHONPATH`:

```bash
T=sidm/studies/vertexChi2/transfer_factor
PYTHONPATH=$PWD python $T/tf_transfer_factor.py     > $T/logs/tf_transfer_factor.log
PYTHONPATH=$PWD python $T/tf_transfer_factor_seg.py > $T/logs/tf_transfer_factor_seg.log
PYTHONPATH=$PWD python $T/tf_seg_aux.py             > $T/logs/tf_seg_aux.log
python $T/extract_records.py
```

Each fold takes about 2-3 minutes on an LPC login node and 0.4 GB of memory; `tf_seg_aux.py` takes seconds. Every path is a command-line argument (`--data-dir`, `--mc-dir`, `--yaml-data`, `--json`; `--seg`, `--dr`, `--aux`, `--outdir` for the extractor), and a relative path is taken from the repository root wherever the script is run from. The full outputs (about 390 kB each) go to `full/`, which is not tracked; the logs are.

The last line of each log names the absolute path the full JSON was written to; it is the only line that depends on where the scripts run.

## Provenance

The scripts were written and run outside the repository and are committed here unchanged except for their paths, which were hard-coded (`/uscms_data/d3/murtazas/SIDM-wt-maria/RunOutputFiles/<run tag>` for the inputs, the same tree's `sidm/configs/ntuples/data_skimmed.yaml`, and an output JSON outside the repository). Original dates: `tf_transfer_factor.py` 2026-09-01, records of 2026-09-02; `tf_transfer_factor_seg.py` and `tf_seg_aux.py` 2026-09-07, records of the same day. The committed records were cut from those outputs; `extract_records.py` was written afterwards to make the cut reproducible.

Rerun on 2026-10-31 against the same inputs (all `.coffea` files predate the original records), the three scripts reproduce the original full JSON outputs byte for byte, and `extract_records.py` reproduces the three committed records byte for byte. The logs here are those of the rerun.

## The records

- `tf_transfer_factor_seg68.json`: the fold with the segment-matched control, 68 DoubleMuon 2018 samples (59.83/fb) and 8 MC samples. Its keys:
  - `provenance`: luminosity, data chunks and MC samples used.
  - `transfer`: per channel and slice, the transfer-matrix rows in pT_PF source-bin order, from [0,20) up to > 500 GeV. Only the MC effective numerator count `mc_neff_num` is kept.
  - `closure`: per channel, slicing (`abseta` or `nHits`), slice index and target. It holds the observed and predicted counts in data and MC, the statistical and three systematic terms (`sig_circ`, `sig_bin`, `sig_chan`), and the data and MC closures `R`, `R_mc` (Obs/Pred - 1) with their uncertainties. For `no_veto` only the `nHits < 20` tail is kept.
  - `inclusive_closure`: the same closure summed over slices.
  - `sideband_transport_failures`: the sideband slices where the prediction exceeds the observation.
- `tf_transfer_factor_full68.json`: the same fold with the dR-matched control, same samples. It keeps only the MC closure of the Spread-veto `nHits < 20` tail (`closure`) and the sideband failures.
- `tf_seg_aux.json`: from `tf_seg_aux.py`. Its keys:
  - `control_sizes`: control-muon counts and feed-up numerators and denominators for both control definitions.
  - `tiebreak`: the fraction of control muons with several PF partners, and the relative shift of the feed-up rate between partner choices.
  - `halo`: the |dxy_BS| and CSC-only fractions of the mu-type lepton-jet DSA muons, in the bulk and the two tail slices, for data and MC.

## Section-10 numbers that are in the logs but not in the records

The notebook's section-10 text quotes some values the records do not carry. They are in the logs:

- 10.1, no-veto column and dR-matched values of the sideband table: the `T7 CLOSURE` table of `logs/tf_transfer_factor_seg.log` (rows `Before`) and of `logs/tf_transfer_factor.log`.
- 10.2, the stat/syst split `R = -0.033 +- 0.033 (stat) +- 0.096 (syst)` and its dR-matched counterpart: the `AS MEASURED` line of each fold log. The factor 2.45 (dR-matched 2.56): the `MC self-closure in the slice under test` line of the segment-matched log and the `pred1/obs` column of each log.
- 10.3, the fit-quality ratios 0.974 +/- 0.025 etc. (dR-matched 0.873 +/- 0.014 etc.) and the 9.5% of lepton-jet tail muons with `sigma(pT)/pT < 0.3`: the `P(relPtErr<0.3 | r>2)` and `LJ DSA muons` lines of each fold log.
- 10.4, the no-veto and control-muon fractions and the halo-corner fraction (10.8%): `HALO DIAGNOSTICS` and `SEGMENT-MATCHED CONTROL` in `logs/tf_seg_aux.log`; the azimuth significance (-3.58 sigma) and the fit probability (3.2e-06): the `CR nHits 0-20 TAIL` phi lines of `logs/tf_transfer_factor_seg.log`. The 10.6% expected under factorization is the product of two logged numbers (`f_H` 0.2712 and the `|eta| 1.2-2.4` fraction 0.389 of the `data/CR phi_tail` row). The MC `|eta|` fractions (69% against 54%) are not in any log; they come from the `vtxRespMC2_2018` outputs.
