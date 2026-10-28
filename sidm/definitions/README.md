# Lepton-jet collections and the analysis pair

The processor keeps diagnostic collections separate from the pair used to define
an analysis channel. `SidmProcessor.select_lepton_jets` applies generic LJ cuts,
builds typed views, applies their cuts, and calls `add_lj_candidates`.

| Collection | Meaning |
|---|---|
| `ljs` | All baseline LJs after generic LJ cuts, in reconstruction pT order. Subtype cuts do not change this diagnostic collection. |
| `mu_ljs`, `egm_ljs` | Baseline mu-type (`muon_n > 0`) and egamma-type (`muon_n == 0`) views after their configured cuts. They can contain more than two LJs. |
| `pfmu_ljs`, `dsamu_ljs`, `electron_ljs`, `photon_ljs` | Pure-subtype diagnostic children, restricted to surviving parents and any cuts applied directly to the child. Mixed constituent jets are not covered by the union of the pure children. Child-only cuts do not select the parent or the analysis pair. |
| `candidate_ljs` | The leading two eligible LJs after all `mu_ljs`/`egm_ljs` cuts, retaining baseline pT order. Each muon candidate must have at least two muons. Fewer than two eligible LJs gives an empty list. |
| `candidate_mu_ljs`, `candidate_egm_ljs` | Muon and egamma members of that same complete pair. |

**Explicit analysis choice for review:** apply all type-specific cuts before
choosing the leading surviving pair. A third/fourth LJ can replace a rejected
leading jet. Extra jets are still allowed; this is not an exactly-two-LJ
requirement. For example, baseline `[mu(100), mu(80), mu(60)]` selects the
80 and 60 GeV jets when the 100 GeV jet fails isolation. Generic cuts and
type-specific cuts both precede pair selection.

`4mu` requires two muon members of the candidate pair; `2mu2e` requires one muon
and one egamma member. Thus they remain exclusive even with extra LJs. A
single-muon leading jet is ineligible and can be replaced by a later eligible jet.
The LJ-LJ dPhi and mass event cuts, all `lj_lj_*` histograms, mixed muon-egamma
mass variants, and `mu_lj_sum_pt_score` use this same pair. Their masks and values
both refer to candidate collections. Incomplete pairs fail the event cuts and
fill no pair histograms, including when the whole input batch is empty.

The generic `>=2 LJs`/exact baseline-count cuts, `lj_n`/inclusive LJ histograms,
leading-object diagnostics and nearest/gen-matched diagnostics retain their
explicit baseline or inclusive typed meaning. In particular, `>=2 LJs` alone
does not assert that an analysis pair survived. Gen-matched masses intentionally
use the jets matched to the generated particles, rather than the analysis pair.
Even in a channel with pair cuts, `lj0_*`/`lj1_*` can describe baseline jets
outside the analysis pair, including jets rejected by subtype cuts.

Tighter subtype cuts can promote a different pair and move an event between
`4mu` and `2mu2e`, or into a channel it did not pass with looser cuts. Channel
yields and pair histograms are therefore not necessarily nested across
selections. Ratios between them must not be interpreted as selection efficiencies.

Candidate eligibility applies to pair histograms in **every** channel. For
example, `baseNoLj_Only_One_Muon` selects `mu_ljs` with `Mu == 1`; it has no
pair event cut, so adding candidates does not impose an event veto there.
Its selected-muon distributions remain available, but no single-muon LJ
can fill `mulj_egmlj_invmass` in `Muon_selection_base` or
`Muon_selection_bkg_base`. Those menus therefore no longer provide the previous
single-muon-plus-egamma mass diagnostic. That diagnostic is intentionally
dropped, with no separately paired replacement. The `Muon_selection` summary
excludes the one-muon channel from its mass comparisons; its other muon
distributions remain available. A pair-histogram integral need not equal the
channel's event yield. Current two-LJ data control selections already
require `Mu >= 2`, as do all current selections using the LJ-LJ dPhi cut.

An LJ criterion must evaluate to true to retain the jet: false or `None` rejects
it. In particular, an undefined reduction over empty constituents is not treated
as a vacuous pass. The current SR `displaced` cuts reduce over **all** muons or
egamma constituents; DSA muons have zero pixel hits and photons have `lostHits`
set to 999 by reconstruction, so DSA-only and photon-only candidates pass these
displacement criteria. The separate `pf_pixelhits <= 2` and `lostHits >= 1` cut
definitions reduce over PF muons and electrons, respectively, and may be
undefined for empty constituents. Neither is used in the current default
`selections.yaml`; historical/custom configurations must specify any different
empty-constituent rule explicitly. No default selection cuts the pure child
views; those remain diagnostic views even when given custom cuts.

Every reconstructed LJ has a stable per-event `idx`, assigned after pT sorting.
Cuts preserve it, including gaps and pT ties. Descendants are intersected by this
identity after a parent cut; they are not rebuilt in a way that resurrects a
previously removed child. Custom diagnostic descendants can be declared in
`lj_parents`, with the same stable-index requirement and an acyclic hierarchy.
No custom post-LJ callback is rerun during synchronization. If a dictionary
already has candidates, further `JaggedSelection` calls refresh them afterward.
Additional generic or subtype cuts can change which surviving jets define the
pair. The processor always applies all generic cuts before building typed views,
then all configured subtype cuts before selecting candidates.
Missing identity during synchronization and missing candidate setup in pair
cuts are errors, not silently skipped selections.

For direct use outside the processor, construct the typed views with stable
indices, apply their cuts, and call `add_lj_candidates(objs)` before evaluating
pair cuts or histograms. Do not reuse a cached candidate dictionary after manually
mutating its parent collections.

Run the dataset-free regressions from the repository root in the analysis
Python environment:

```bash
PYTHONPATH=. python tests/test_lj_collections.py
```

These tests check software membership, identity and pair agreement using literal
records. They do not measure physics acceptance or validate real inputs.
