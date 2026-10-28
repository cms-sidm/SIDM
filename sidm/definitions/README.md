# Lepton-jet collections

Use `candidate_ljs` for quantities that describe the analysis pair, and the
baseline or typed collections for inclusive object studies. The processor
prepares these collections before evaluating event cuts and filling histograms.

Object views are defined in [objects.py](objects.py), cuts in [cuts.py](cuts.py),
and histograms in [hists.py](hists.py). Select the cuts and histograms to run in
[selections.yaml](../configs/selections.yaml) and
[hist_collections.yaml](../configs/hist_collections.yaml).

## Choosing a collection

| Collection | Contents |
|---|---|
| `ljs` | All LJs passing generic LJ cuts, in reconstruction pT order. Type-specific cuts do not remove jets from this collection. |
| `mu_ljs`, `egm_ljs` | Muon-type (`muon_n > 0`) and egamma-type (`muon_n == 0`) LJs passing their configured cuts. Each collection can contain more than two jets. |
| `pfmu_ljs`, `dsamu_ljs` | Pure PF-muon and pure DSA-muon diagnostic views of `mu_ljs`. |
| `electron_ljs`, `photon_ljs` | Pure electron and pure photon diagnostic views of `egm_ljs`. |
| `candidate_ljs` | The two highest-pT eligible surviving LJs, or an empty list if fewer than two survive. Muon candidates must contain at least two muons. |
| `candidate_mu_ljs`, `candidate_egm_ljs` | Muon and egamma members of that same complete pair. |

Mixed PF/DSA jets belong to `mu_ljs` but neither pure-muon view; mixed
electron/photon jets likewise belong to `egm_ljs` but neither pure-egamma view.
Use the parent collection when a selection should include mixed jets.

## Selecting the analysis pair

For each selection, `SidmProcessor.select_lepton_jets` applies generic LJ cuts,
builds the typed views, applies their cuts, and selects the leading eligible
surviving pair. Extra jets are allowed. For example, three muon LJs with pT
100, 80 and 60 GeV select the 80 and 60 GeV jets if the 100 GeV jet fails its
type-specific cuts, provided each survivor contains at least two muons.

The `4mu` event cut requires two muon members of this pair; `2mu2e` requires one
muon and one egamma member. The two channels are mutually exclusive within a
selection. Fewer than two eligible survivors fails both channel cuts and fills
no pair histograms. Passing a baseline LJ-count cut such as `>=2 LJs` alone
does not guarantee a valid analysis pair.

The LJ-LJ dPhi and mass cuts, `lj_lj_*` histograms, mixed muon-egamma mass
histograms and `mu_lj_sum_pt_score` use this pair. When adding pair observables,
use candidate collections for both their masks and values.

Inclusive and leading-object diagnostics use their named collections:
`lj0_*` and `lj1_*` refer to the leading baseline jets, which may differ from
the analysis pair. Gen-matched histograms use the jets matched to generated
particles. Single-muon LJs remain available in `mu_ljs` diagnostics but cannot
fill candidate mass histograms, including the `Muon_selection` mass plots.
A pair-histogram integral therefore need not equal the selection's event yield.

Tighter cuts can select a different pair, changing an event's pair mass or
channel. Pair distributions from different selections need not be nested.

## Adding cuts or using the collections directly

Place cuts that should affect the analysis pair on `ljs`, `mu_ljs` or
`egm_ljs`. A parent cut also removes rejected jets from its diagnostic children.
Cuts on a pure diagnostic child affect only that view; they do not select its
parent or the analysis pair. An LJ cut retains only entries whose mask is true;
false and `None` reject the jet. Define the intended behavior of reductions over
empty constituent collections when adding a custom cut.

Every reconstructed LJ has a stable per-event `idx`. Preserve it through cuts;
it identifies membership even after filtering leaves gaps or jets have equal
pT. Custom diagnostic descendants can be registered in `lj_parents` with an
acyclic parent hierarchy and the same indices. Synchronization preserves cuts
already applied to children and does not rerun their construction callbacks.

Outside the processor, provide pT-ordered `ljs` and the selected `mu_ljs` and
`egm_ljs` views with matching indices, then finalize the pair:

```python
from sidm.definitions.objects import add_lj_candidates

objs = add_lj_candidates(objs)
pair = objs["candidate_ljs"]
```

Call this before evaluating pair cuts or histograms. `JaggedSelection` refreshes
an existing pair after further object cuts; recompute it explicitly after
manually changing its parent collections.
