#!/usr/bin/env python3
"""Print SIDM YAML cuts after resolving anchors, merges, and nested lists.

Requires only PyYAML; does not import Coffea or run any cut functions.
Run with no arguments to print every selection in configs/selections.yaml.
"""

import argparse
from collections import Counter
from difflib import SequenceMatcher
import fnmatch
import json
import os
from pathlib import Path
import sys

try:
    import yaml
except ImportError:
    raise SystemExit("PyYAML is required. Install it with: python3 -m pip install PyYAML")


DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "selections.yaml"
SECTIONS = {
    "obj_cuts": "Object cuts (before LJ reconstruction)",
    "preLj_obj_cuts": "Additional pre-LJ cuts (configured; not applied by current process())",
    "postLj_obj_cuts": "LJ and post-LJ object cuts",
    "evt_cuts": "Event cuts",
}


def use_color(mode):
    """Explicit modes override auto detection (including NO_COLOR)."""
    if mode != "auto":
        return mode == "always"
    return sys.stdout.isatty() and "NO_COLOR" not in os.environ and os.environ.get("TERM") != "dumb"


def paint(text, style, color):
    codes = {"title": "1;35", "section": "1;36", "object": "1;33",
             "same": "2", "add": "32", "remove": "31", "warning": "33"}
    return f"\033[{codes[style]}m{text}\033[0m" if color else text


def flatten_cuts(value, context="cuts", active=None):
    """Flatten list/dict values like utilities.flatten, retaining repeats/order.

    Dict values matter: an event list can contain an entire aliased selection,
    e.g. {'evt_cuts': ['pass triggers']}. Reject cycles and non-string leaves.
    """
    if isinstance(value, str):
        return [value]
    if not isinstance(value, (list, dict)):
        raise ValueError(f"{context}: expected a cut string, list, or mapping; got {value!r}")
    active = set() if active is None else active
    if id(value) in active:
        raise ValueError(f"{context}: recursive YAML alias")
    active.add(id(value))
    try:
        items = value.values() if isinstance(value, dict) else value
        return [cut for item in items for cut in flatten_cuts(item, context, active)]
    finally:
        active.remove(id(value))


def read_selections(path):
    """Use the same safe YAML loader as SIDM, including merge/override behavior."""
    with Path(path).open(encoding="utf-8") as stream:
        menu = yaml.safe_load(stream)
    if not isinstance(menu, dict):
        raise ValueError("Expected a mapping of selection names to cut definitions")
    selections = {
        name: value for name, value in menu.items()
        if isinstance(value, dict) and any(key in value for key in SECTIONS)
    }
    if not selections:
        raise ValueError("No selections found (expected obj_cuts or evt_cuts sections)")
    if any(not isinstance(name, str) for name in selections):
        raise ValueError("Selection names must be strings")
    return selections


def expand_selection(selection, name):
    expanded = {}
    for section in SECTIONS:
        if section not in selection:
            continue
        value = selection[section]
        context = f"{name}.{section}"
        if section == "evt_cuts":
            expanded[section] = flatten_cuts(value, context)
        else:
            if not isinstance(value, dict):
                raise ValueError(f"{context}: expected a mapping of object names to cuts")
            expanded[section] = {
                obj: flatten_cuts(cuts, f"{context}.{obj}")
                for obj, cuts in value.items()
            }
    return expanded


def format_selection(name, selection, color=False):
    lines = [paint(f"=== {name} ===", "title", color)]

    def numbered(cuts, indent):
        if not cuts:
            lines.append(f"{indent}(none)")
        for index, cut in enumerate(cuts, 1):
            lines.append(f"{indent}{index}. {cut}")

    for section, value in selection.items():
        lines.append(paint(f"  {SECTIONS[section]} [{section}]", "section", color))
        if section == "evt_cuts":
            numbered(value, "    ")
        elif not value:
            lines.append("    (none)")
        else:
            for obj, cuts in value.items():
                lines.append(paint(f"    {obj}:", "object", color))
                numbered(cuts, "      ")
    return "\n".join(lines)


def compare_selections(left_name, left, right_name, right):
    """Compare labels per section/collection, including multiplicity and order.

    Missing sections/collections and empty ones both impose zero cuts. A moved
    cut appears as a removal/addition in the ordered diff, not a new condition.
    """
    groups = []
    for section in SECTIONS:
        if section not in left and section not in right:
            continue
        if section == "evt_cuts":
            objects = [None]
        else:
            objects = list(dict.fromkeys([*left.get(section, {}), *right.get(section, {})]))
        for obj in objects:
            a = left.get(section, []) if obj is None else left.get(section, {}).get(obj, [])
            b = right.get(section, []) if obj is None else right.get(section, {}).get(obj, [])
            rows = []
            for tag, i, end_i, j, end_j in SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
                if tag == "equal":
                    rows.extend({"kind": "same", "left_index": x + 1,
                                 "right_index": y + 1, "cut": a[x]}
                                for x, y in zip(range(i, end_i), range(j, end_j)))
                else:
                    rows.extend({"kind": "remove", "left_index": x + 1,
                                 "right_index": None, "cut": a[x]} for x in range(i, end_i))
                    rows.extend({"kind": "add", "left_index": None,
                                 "right_index": y + 1, "cut": b[y]} for y in range(j, end_j))
            groups.append({"section": section, "object": obj,
                           "same_cuts": Counter(a) == Counter(b), "same_order": a == b,
                           "rows": rows})
    return {"left": left_name, "right": right_name,
            "same_cuts": all(group["same_cuts"] for group in groups),
            "same_order": all(group["same_order"] for group in groups), "groups": groups}


def format_comparison(comparison, color=False):
    left, right = comparison["left"], comparison["right"]
    lines = [paint(f"=== {left} vs {right} ===", "title", color)]
    for key, label in [("same_cuts", "Same cuts per group (including repeats)"),
                       ("same_order", "Same cut order within each group")]:
        answer = comparison[key]
        lines.append(paint(f"{label}: {'YES' if answer else 'NO'}",
                           "add" if answer else "warning", color))
    lines.extend(["", "Ordered diff: = shared, - left removal, + right addition; [left:right] positions.",
                  "Moved cuts appear as -/+; absent and empty groups both mean no cuts.",
                  "Compares configured labels, not function behavior or event-level acceptance."])
    current_section = None
    for group in comparison["groups"]:
        section = group["section"]
        if section != current_section:
            lines.append(paint(f"\n  {SECTIONS[section]} [{section}]", "section", color))
            current_section = section
        if group["object"] is not None:
            lines.append(paint(f"    {group['object']}:", "object", color))
        if group["same_cuts"] and not group["same_order"]:
            lines.append(paint("      Same cuts, different order", "warning", color))
        if not group["rows"]:
            lines.append("      (none in either selection)")
        for row in group["rows"]:
            marker = {"same": "=", "remove": "-", "add": "+"}[row["kind"]]
            position = f"{row['left_index'] or '-'}:{row['right_index'] or '-'}"
            lines.append(paint(f"      {marker} [{position}] {row['cut']}", row["kind"], color))
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", nargs="?", type=Path, default=DEFAULT_CONFIG,
                        help="selection YAML file (default: bundled selections.yaml)")
    parser.add_argument("-s", "--selection", action="append", default=[], metavar="NAME",
                        help="selection name or quoted glob; repeat to select several")
    parser.add_argument("--list", action="store_true", help="list matching selection names only")
    parser.add_argument("--json", action="store_true", help="output expanded cuts as JSON")
    parser.add_argument("--compare", nargs=2, metavar=("LEFT", "RIGHT"),
                        help="compare two exact selection names (e.g. denominator numerator)")
    parser.add_argument("--color", choices=("auto", "always", "never"), default="auto",
                        help="terminal colors (default: auto; honors NO_COLOR; JSON stays plain)")
    args = parser.parse_args(argv)
    if args.compare and (args.selection or args.list):
        parser.error("--compare cannot be combined with --selection or --list")
    try:
        menu = read_selections(args.config)
        color = use_color(args.color)
        if args.compare:
            for name in args.compare:
                if name not in menu:
                    raise ValueError(f"Unknown selection {name!r}; use --list to see names")
            left, right = args.compare
            comparison = compare_selections(left, expand_selection(menu[left], left),
                                            right, expand_selection(menu[right], right))
            if args.json:
                print(json.dumps(comparison, indent=2))
            else:
                print(f"Config: {args.config.resolve()}\n")
                print(format_comparison(comparison, color))
            return
        for pattern in args.selection:
            if not any(fnmatch.fnmatchcase(name, pattern) for name in menu):
                raise ValueError(f"No selection matches {pattern!r}; use --list to see names")
        names = [name for name in menu if not args.selection or any(
            fnmatch.fnmatchcase(name, pattern) for pattern in args.selection)]
        if args.list:
            print("\n".join(paint(name, "title", color) for name in names))
            return
        expanded = {name: expand_selection(menu[name], name) for name in names}
        if args.json:
            print(json.dumps(expanded, indent=2))
        else:
            print(f"Config: {args.config.resolve()}")
            print("Configured cut labels; inheritance expanded, order and repeats preserved.\n")
            print("\n\n".join(format_selection(name, cuts, color) for name, cuts in expanded.items()))
    except (OSError, ValueError, yaml.YAMLError, RecursionError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # Allow piping a long report into head/less without a traceback.
        raise SystemExit(0)
