# Chain report (CI)

A reviewer-facing report that runs automatically on every pull request. It runs the
SidmProcessor over a small committed fixture for **both `main` and your PR** and shows,
side by side, how the processing + plotting chain behaves — so a reviewer can see at a
glance what your change does to the analysis.

## What it shows
- **Cutflows** — per channel, the cumulative count of events passing each **event-level**
  cut (object-level cuts slim collections but add no cutflow rows), **raw and weighted**,
  shown `main` → your PR. A *regression* check: it flags when a surviving count moves, not
  whether the count is physically correct. The weighted column is what moves when a
  weight/normalization change leaves the raw counts untouched.
- **Selections** — how many of the ~135 channels executed, plus a ❌ banner if the chain crashes.
- **Hist collections / hists** — how many collections are valid vs. broken (referencing
  undefined hists), **which hists each broken collection is missing**, and which hists
  could not be filled at runtime.
- **Warnings** — the processor's "Unable to apply … skipping" messages, and which ones your
  PR **adds** or **removes**.
- **Channel health** — which channels this fixture cannot regression-test (see Scope below).

## What you do as a contributor
**Nothing special.** Open or update a PR and the report runs on its own (~10–15 min). To read it:

1. On the PR, open the **Checks** tab → **chain-report**.
2. Read the **job summary** (rendered at the top): a ✅/❌ banner, a current-vs-PR summary
   table, and a **"What this PR changes"** section.

The report **does not block merging** — it informs the reviewer, who decides. The one
exception is a **soft red ❌**: the check goes red only if your PR introduces a **new
error**, which means any of:

- a **newly broken or worsened hist collection** (a collection referencing an undefined
  hist it didn't reference before — including adding a missing hist to a collection that
  was already broken);
- a **chain crash** on your PR (or a failure to build the hist-collection inventory) —
  a crash on **both** sides is also red: either an environment/harness change in your PR
  broke the chain (only `sidm/` is swapped for the base run), or `main` itself is broken
  (in which case every open PR's report goes red until it's fixed);
- a **hist that newly fails to evaluate/fill at runtime** (skipped when the two sides'
  hist passes covered different channels/collections — then the delta is shown but not
  used for the red check, since it isn't attributable to your PR);
- a **cut that newly fails to evaluate or apply** — matched against the processor's
  cut-failure warnings by exact class, so removing a cut *on purpose* (which emits no
  warning) does not go red; it just shows in the cutflow diff.

**New warnings outside those failure classes do not fail it** (they're shown for
awareness, the check stays green).

If your PR goes red: the ❌ section at the top of the job summary names the exact
collection, cut, hist, or crashed side. Reproduce locally with the two commands below
and iterate.

If you **intentionally** changed selections, cuts, or hist collections, that's fine — the
report just shows the resulting cutflow/warning deltas for the reviewer to confirm. There
is **no baseline to update**: the "current" side is recomputed from `main` every run (the
first parent of the PR's merge commit, i.e. the exact `main` your PR is being merged into).

## Pre-existing issues
The report currently surfaces some pre-existing rot. These are **not** introduced by your
PR and will **not** fail it; they're shown so they can be fixed over time:

- **6 hist collections reference removed `*_lj_isolation` hists**, including the main
  `base` collection — which is why `test_SidmProcessor.ipynb` (which runs `base`) is
  currently broken.
- The `good_matched_muons` skip (`sidm_processor.py:178`) is **not harmless**: that
  collection feeds the dsaMuon cross-cleaning veto (`cuts.py:208-213`), so in the
  `matched_ljsource_leptons` and `baseNoLj_matchedLeptons` channels the veto runs on an
  unmatched-muon input (791 vs 749 muons on this fixture). A fix exists in the form of
  the nested `dR(mu, A) < 0.5` cut variant at `cuts.py:163`. Only the
  `good_matched_dsa_muons` sibling (line 179) is diagnostic-only.
- An undefined `pT < 650 GeV` ljs cut (`selections.yaml:1475`) — the one channel using it
  (`baseNoLj_A_mumu_matched_lj_noPfMatch_highLjPt`) drops it with a catch-and-skip
  warning (visible in this report's baseline warnings).
- Benign control-region `event cuts to …` skips.

## Run it locally (optional)
In the analysis venv (it runs the processor):
```bash
python tests/chain_report.py compute my_state.json            # this checkout's chain state
python tests/chain_report.py render base.json my_state.json   # markdown diff of two states
```

## Scope
This verifies **execution + regression against `main`**, not physics correctness — a
silently wrong value or a mis-set cut can still execute cleanly.

All ~135 channels *execute* in the cutflow pass, but on this fixture not all of them are
*regression-sensitive*: ~10 channels (the 4mu family and a few high-purity selections)
end at 0 surviving events (the fixture is a 2Mu2E signal sample), and ~15 (the
object-cut-only families, e.g. `barrelE`, `ljsource_cuts`) apply no event-level cuts, so
their cutflow holds only the initial event total. The report's channel-health lines list
exactly which channels are in those two states on each side. Separately, channels whose
*distinguishing* cuts are all object-level (e.g. the electron-ID `sieie`/`hoe` family)
carry only generic trigger/PV cutflow rows — their ID content never appears in a cutflow
and is exercised only by the hist pass. A small 4Mu fixture to cover the muon families is
a planned follow-up. Warnings are collected from **both** passes (the all-channel cutflow
pass — which also fills its one collection, so it can emit hist warnings for it across
all channels — and the 14-representative-channel × valid-collections pass). The committed
fixture (`tests/data/events_2mu2e_500GeV_200ev.root`, a 200-event 2Mu2E skim — a
full-width snapshot, all branches retained) lets the chain run on a stock CI runner with
no EOS/XRootD/cvmfs access. The harness runs with a deliberately non-unit `skim_factor`
so the skim-scaling path of the weighted cutflow column is genuinely tested. The
processor applies the standard corrections to simulation by default, so
`tests/data/corrections/` holds copies of the correction files that the default list and
the `corrections_*` histograms read (about 800 kB, the CVMFS layout under
`SIDM_CORRECTIONS_ROOT`); a correction file that the processor or a histogram starts
reading needs copying there too.

The data path (`is_data: True`) is exercised by a second fixture,
`tests/data/events_DoubleMuon2018C_200ev.root`: a 200-event slice of a real DoubleMuon
2018C skim (`skimmed_output_583.root`, entries 11649-11849) chosen to span a golden-JSON
dead-run boundary. Its first 100 events are run 320007, which is absent from the 2018
golden JSON entirely; the last 100 are run 319993, which is fully golden. Both fixtures
are processed at `chunksize=100`, so the data fixture's first chunk is 100% non-golden:
any data-only code that filters on the golden JSON produces a zero-event chunk here, the
exact configuration that crashed real 2018C processing and that an MC fixture can never
probe (the mask only runs when `is_data` is true). The report carries a separate
data-fixture cutflow diff and goes red if the data pass crashes on the PR side. MC counts
are unaffected by the chunking: accumulation is a sum over chunks.
