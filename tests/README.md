# Chain report (CI)

A reviewer-facing report that runs automatically on every pull request. It runs the
SidmProcessor over a small committed fixture for **both `main` and your PR** and shows,
side by side, how the processing + plotting chain behaves — so a reviewer can see at a
glance what your change does to the analysis.

## What it shows
- **Cutflows** — raw per-channel event counts, current (`main`) → your PR, with any changes called out.
- **Selections** — how many of the ~135 channels executed; anything that errors out.
- **Hist collections / hists** — how many collections are valid vs. broken (referencing undefined hists), and which hists could not be filled.
- **Warnings** — the processor's "Unable to apply … skipping" messages, and which ones your PR **adds** or **removes**.

## What you do as a contributor
**Nothing special.** Open or update a PR and the report runs on its own (~10–15 min). To read it:

1. On the PR, open the **Checks** tab → **chain-report**.
2. Read the **job summary** (rendered at the top): a ✅/❌ banner, a current-vs-PR summary table, and a **"What this PR changes"** section.

The report **does not block merging** — it informs the reviewer, who decides. The one exception is a **soft red ❌**: the check goes red only if your PR introduces a **new error** — a newly *broken* hist collection, or a chain *crash*. **New warnings do not fail it** (they're shown for awareness, the check stays green).

If you **intentionally** changed selections, cuts, or hist collections, that's fine — the report just shows the resulting cutflow/warning deltas for the reviewer to confirm. There is **no baseline to update**: the "current" side is recomputed from `main` every run.

## Pre-existing issues
The report currently surfaces some pre-existing rot (a handful of hist collections that
reference removed `*_lj_isolation` hists — including `base` — and a set of harmless
catch-and-skip warnings). These are **not** introduced by your PR and will **not** fail it;
they're shown so they can be fixed over time.

## Run it locally (optional)
In the analysis venv (it runs the processor):
```bash
python tests/chain_report.py compute my_state.json            # this checkout's chain state
python tests/chain_report.py render base.json my_state.json   # markdown diff of two states
```

## Scope
This verifies **execution + regression against `main`**, not physics correctness — a
silently wrong value or a mis-set cut can still execute cleanly. The committed fixture
(`tests/data/events_2mu2e_500GeV_200ev.root`, a 200-event 2Mu2E skim) lets the chain run
on a stock CI runner with no EOS/XRootD/cvmfs access.
