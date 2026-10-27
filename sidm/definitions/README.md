# Lepton-jet collections and the analysis pair

The processor keeps diagnostic collections separate from the pair used to define
an analysis channel. `SidmProcessor.select_lepton_jets` applies generic LJ cuts,
builds typed views, applies their cuts, and calls `add_lj_candidates`.

| Collection | Meaning |
|---|---|
| `ljs` | All baseline LJs after generic LJ cuts, in reconstruction pT order. Subtype cuts do not change this diagnostic collection. |
| `mu_ljs`, `egm_ljs` | Baseline mu-type (`muon_n > 0`) and egamma-type (`muon_n == 0`) views after their configured cuts. They can contain more than two LJs. |
| `pfmu_ljs`, `dsamu_ljs`, `electron_ljs`, `photon_ljs` | Pure-subtype diagnostic children, restricted to surviving parents and any cuts applied directly to the child. Mixed constituent jets are not covered by the union of the pure children. Child-only cuts do not select the parent or the analysis pair. |
| `candidate_ljs` | The original leading two baseline LJs, if **both** survive the `mu_ljs`/`egm_ljs` cuts and every muon candidate has at least two muons; otherwise an empty list. |
| `candidate_mu_ljs`, `candidate_egm_ljs` | Muon and egamma members of that same complete pair. |

**Explicit analysis choice for review:** require the original leading pair to
survive. Do not promote a third/fourth LJ after a leading LJ fails a subtype cut.
Extra jets are still allowed; this is not an exactly-two-LJ requirement. For
example, baseline `[mu(100), mu(80), mu(60)]` fails pair selection if the 100 GeV
jet fails isolation, even when the other two survive. The baseline pair is chosen
after generic LJ cuts, so a jet rejected by those cuts was never part of it.

`4mu` requires two muon members of the candidate pair; `2mu2e` requires one muon
and one egamma member. Thus they remain exclusive even with extra LJs. A
single-muon leading jet cannot be a candidate or be bypassed by a later jet.
The LJ-LJ dPhi and mass event cuts, all `lj_lj_*` histograms, mixed muon-egamma
mass variants, and `mu_lj_sum_pt_score` use this same pair. Their masks and values
both refer to candidate collections. Incomplete pairs fail the event cuts and
fill no pair histograms, including when the whole input batch is empty.

The generic `>=2 LJs`/exact baseline-count cuts, `lj_n`/inclusive LJ histograms,
leading-object diagnostics and nearest/gen-matched diagnostics retain their
explicit baseline or inclusive typed meaning. In particular, `>=2 LJs` alone
does not assert that an analysis pair survived. Gen-matched masses intentionally
use the jets matched to the generated particles, rather than the analysis pair.

Candidate eligibility applies to pair histograms in **every** channel. For
example, `baseNoLj_Only_One_Muon` selects `mu_ljs` with `Mu == 1`; it has no
pair event cut, so adding candidates does not impose an event veto there.
Its selected-muon distributions remain available, but a leading single-muon LJ
cannot fill `mulj_egmlj_invmass` in `Muon_selection_base` or
`Muon_selection_bkg_base`. Those menus therefore no longer provide the previous
single-muon-plus-egamma mass diagnostic. This affects comparisons in the
`Muon_selection` study: a pair-histogram integral need not equal the channel's
event yield. The pair name always means the analysis candidates; a diagnostic
with different eligibility or independently chosen typed leaders needs a
separate explicit definition. Current two-LJ data control selections already
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
An additional generic `ljs` cut changes the baseline from which the leading two
are defined; a subtype cut only tests survival of that baseline pair. The
processor always applies all generic cuts before building typed views.
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
