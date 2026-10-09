#!/usr/bin/env python3
"""Print SIDM YAML cuts after resolving anchors, merges, and nested lists.

Requires only PyYAML; does not import Coffea or run any cut functions.
Run with no arguments to print every selection in configs/selections.yaml.
"""

import argparse
import fnmatch
import json
from pathlib import Path

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


def format_selection(name, selection):
    lines = [f"=== {name} ==="]

    def numbered(cuts, indent):
        if not cuts:
            lines.append(f"{indent}(none)")
        for index, cut in enumerate(cuts, 1):
            lines.append(f"{indent}{index}. {cut}")

    for section, value in selection.items():
        lines.append(f"  {SECTIONS[section]} [{section}]")
        if section == "evt_cuts":
            numbered(value, "    ")
        elif not value:
            lines.append("    (none)")
        else:
            for obj, cuts in value.items():
                lines.append(f"    {obj}:")
                numbered(cuts, "      ")
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", nargs="?", type=Path, default=DEFAULT_CONFIG,
                        help="selection YAML file (default: bundled selections.yaml)")
    parser.add_argument("-s", "--selection", action="append", default=[], metavar="NAME",
                        help="selection name or quoted glob; repeat to select several")
    parser.add_argument("--list", action="store_true", help="list matching selection names only")
    parser.add_argument("--json", action="store_true", help="output expanded cuts as JSON")
    args = parser.parse_args(argv)
    try:
        menu = read_selections(args.config)
        for pattern in args.selection:
            if not any(fnmatch.fnmatchcase(name, pattern) for name in menu):
                raise ValueError(f"No selection matches {pattern!r}; use --list to see names")
        names = [name for name in menu if not args.selection or any(
            fnmatch.fnmatchcase(name, pattern) for pattern in args.selection)]
        if args.list:
            print("\n".join(names))
            return
        expanded = {name: expand_selection(menu[name], name) for name in names}
        if args.json:
            print(json.dumps(expanded, indent=2))
        else:
            print(f"Config: {args.config.resolve()}")
            print("Configured cut labels; inheritance expanded, order and repeats preserved.\n")
            print("\n\n".join(format_selection(name, cuts) for name, cuts in expanded.items()))
    except (OSError, ValueError, yaml.YAMLError, RecursionError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # Allow piping a long report into head/less without a traceback.
        raise SystemExit(0)
