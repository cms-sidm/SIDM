"""Per-PR chain report for SIDM.

Runs the SidmProcessor over a small committed fixture and reports the state of the
processing + plotting chain, so a reviewer can see -- side by side, current (main)
vs. this PR -- what executes, what errors, what warns, and how the cutflows move.
It does NOT gate the merge; it informs the reviewer.

Two modes:
    python tests/chain_report.py compute  > state.json          # run the chain, dump state
    python tests/chain_report.py render base.json pr.json        # diff two states -> markdown

CI runs `compute` once on the base commit and once on the PR, then `render`s the
diff into the GitHub Actions job summary.

NOTE: this reports EXECUTION + regression, not physics correctness -- a silently
wrong value or cut can still execute cleanly.
"""
import json
import os
import re
import sys
import traceback

import matplotlib
matplotlib.use("Agg")  # headless: utilities imports pyplot at import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _e2e_common as e2e

_SKIPPED_HIST_RE = re.compile(r"histogram with the name (\w+) could not be evaluated")


def compute_state():
    """Run the chain over the fixture and return a JSON-able state dict."""
    state = {
        "n_channels_run": 0, "cutflows": {}, "broken_collections": {},
        "valid_collection_count": 0, "warnings": [], "skipped_hists": [], "error": None,
    }
    try:
        out_sel, warns_sel = e2e.run_chain(e2e.all_channels(), ["pv_base"])
        state["cutflows"] = e2e.cutflow_counts(out_sel)
        state["n_channels_run"] = len(state["cutflows"])
        _, warns_hist = e2e.run_chain(e2e.REPRESENTATIVE_CHANNELS, e2e.valid_collections())
        warnings = sorted(warns_sel | warns_hist)
        state["warnings"] = warnings
        state["skipped_hists"] = sorted({m.group(1) for w in warnings
                                         for m in [_SKIPPED_HIST_RE.search(w)] if m})
        state["broken_collections"] = e2e.broken_collections()
        state["valid_collection_count"] = len(e2e.valid_collections())
    except Exception:
        state["error"] = traceback.format_exc()
    return state


def _fmt_list(items, empty="_none_"):
    return "\n".join(f"- `{x}`" for x in items) if items else empty


def new_errors(base, pr):
    """Reasons this PR introduces NEW errors -> drives the soft red check.
    New *warnings* do NOT count (they stay informational/green)."""
    reasons = []
    if pr.get("error") and not base.get("error"):
        reasons.append("the chain crashed on this PR (traceback in the report)")
    newly_broken = sorted(set(pr.get("broken_collections", {})) - set(base.get("broken_collections", {})))
    if newly_broken:
        reasons.append(f"newly broken hist collection(s): {newly_broken}")
    return reasons


def render(base, pr):
    """Render a markdown diff report (current=base vs after=pr)."""
    out = ["## SIDM chain report", "", "_Current (`main`) → this PR. Reviewer-facing and advisory "
           "(does not hard-block merge). Verifies execution + regression, not physics correctness._", ""]
    _errs = new_errors(base, pr)
    out += [("### ❌ This PR introduces new errors\n\n" + "\n".join(f"- {r}" for r in _errs))
            if _errs else "### ✅ No new errors introduced by this PR", ""]

    if pr.get("error"):
        out += ["### ❌ The chain ERRORED on this PR", "",
                "```", pr["error"].strip()[-2500:], "```", ""]
    if base.get("error"):
        out += ["> Note: the base (`main`) chain also errored; comparison is partial.", ""]

    b_brk, p_brk = base.get("broken_collections", {}), pr.get("broken_collections", {})
    b_warn, p_warn = set(base.get("warnings", [])), set(pr.get("warnings", []))
    b_skip, p_skip = set(base.get("skipped_hists", [])), set(pr.get("skipped_hists", []))
    cf_b, cf_p = base.get("cutflows", {}), pr.get("cutflows", {})
    changed = sorted(c for c in set(cf_b) & set(cf_p) if cf_b[c] != cf_p[c])
    added_ch = sorted(set(cf_p) - set(cf_b))
    removed_ch = sorted(set(cf_b) - set(cf_p))

    def arrow(a, b_):
        return f"{a}" if a == b_ else f"**{a} → {b_}**"

    out += [
        "### Summary", "",
        f"| metric | current | this PR |",
        f"|---|--:|--:|",
        f"| channels executed | {base.get('n_channels_run',0)} | {pr.get('n_channels_run',0)} |",
        f"| valid hist collections | {base.get('valid_collection_count',0)} | {pr.get('valid_collection_count',0)} |",
        f"| broken collections | {len(b_brk)} | {len(p_brk)} |",
        f"| catch-and-skip warnings | {len(b_warn)} | {len(p_warn)} |",
        f"| hists skipped (could not evaluate) | {len(b_skip)} | {len(p_skip)} |",
        f"| channels with changed cutflow | — | {len(changed)} |",
        "",
    ]

    # Highlight what THIS PR changes
    new_warn, gone_warn = sorted(p_warn - b_warn), sorted(b_warn - p_warn)
    new_brk = sorted(set(p_brk) - set(b_brk))
    fixed_brk = sorted(set(b_brk) - set(p_brk))
    out += ["### What this PR changes", ""]
    if not (new_warn or gone_warn or new_brk or fixed_brk or changed or added_ch or removed_ch or pr.get("error")):
        out += ["No change to chain execution, collections, warnings, or cutflows. ✅", ""]
    else:
        if new_brk:
            out += [f"**Newly broken collections ({len(new_brk)}):**", _fmt_list(new_brk), ""]
        if fixed_brk:
            out += [f"**Collections repaired ({len(fixed_brk)}):**", _fmt_list(fixed_brk), ""]
        if new_warn:
            out += [f"**New warnings ({len(new_warn)}):**", _fmt_list(new_warn), ""]
        if gone_warn:
            out += [f"**Warnings no longer emitted ({len(gone_warn)}):**", _fmt_list(gone_warn), ""]
        if added_ch or removed_ch:
            out += [f"**Channels added:** {added_ch or '_none_'}  **removed:** {removed_ch or '_none_'}", ""]
        if changed:
            out += [f"**Cutflow changed in {len(changed)} channel(s):**", "",
                    "| channel | cut | current | this PR |", "|---|---|--:|--:|"]
            for c in changed[:50]:
                cuts = sorted(set(cf_b[c]) | set(cf_p[c]))
                for cut in cuts:
                    a, b_ = cf_b[c].get(cut), cf_p[c].get(cut)
                    if a != b_:
                        out.append(f"| {c} | {cut} | {a} | {b_} |")
            if len(changed) > 50:
                out.append(f"\n_…and {len(changed) - 50} more changed channels._")
            out.append("")

    # Always-present current-state detail (collapsed)
    out += ["<details><summary>Current broken collections (baseline state)</summary>", ""]
    out += [f"- `{c}` → missing {sorted(v)}" for c, v in sorted(b_brk.items())] or ["_none_"]
    out += ["", "</details>", ""]
    out += ["<details><summary>Current catch-and-skip warnings (baseline state)</summary>", ""]
    out += [_fmt_list(sorted(b_warn))]
    out += ["", "</details>", ""]
    return "\n".join(out)


def main(argv):
    # compute writes JSON to a FILE (not stdout): the chain emits banner/progress
    # text to stdout at the C level (fastjet), which would corrupt a stdout dump.
    if len(argv) == 2 and argv[0] == "compute":
        with open(argv[1], "w") as fh:
            json.dump(compute_state(), fh, indent=2, sort_keys=True)
    elif len(argv) == 3 and argv[0] == "render":
        base = json.load(open(argv[1]))
        pr = json.load(open(argv[2]))
        print(render(base, pr))
        if new_errors(base, pr):   # soft red check -- report is already printed
            sys.exit(1)
    else:
        sys.exit("usage: chain_report.py compute OUT.json | render base.json pr.json")


if __name__ == "__main__":
    main(sys.argv[1:])
