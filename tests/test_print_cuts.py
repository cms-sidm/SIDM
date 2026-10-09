"""Tests requiring only the standard library and PyYAML.

Run: python3 -m unittest discover -s tests -p test_print_cuts.py
"""

import ast
import importlib.util
import io
import json
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
import unittest
from unittest.mock import patch

import yaml


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("print_cuts", ROOT / "sidm/scripts/print_cuts.py")
printer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(printer)


class PrintCutsTests(unittest.TestCase):
    def test_merge_override_nested_dict_and_repeats(self):
        menu = yaml.safe_load('''
trigger: &trigger
  evt_cuts: [trigger]
base: &base
  obj_cuts: &objects
    muons: [loose, pt]
    electrons: [id]
  evt_cuts: &events
    - *trigger
    - [pv, pv]
child:
  <<: *base
  obj_cuts:
    <<: *objects
    muons: [tight]
  evt_cuts: [*events, channel]
''')
        cuts = printer.expand_selection(menu["child"], "child")
        self.assertEqual(cuts["obj_cuts"], {"muons": ["tight"], "electrons": ["id"]})
        self.assertEqual(cuts["evt_cuts"], ["trigger", "pv", "pv", "channel"])

    def test_all_real_selections_match_processor_flatten(self):
        # Extract the actual pure helper without importing its physics dependencies.
        source = (ROOT / "sidm/tools/utilities.py").read_text()
        tree = ast.parse(source)
        node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "flatten")
        namespace = {}
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(ROOT / "sidm/tools/utilities.py"), "exec"), namespace)
        flatten = namespace["flatten"]
        menu = printer.read_selections(printer.DEFAULT_CONFIG)
        self.assertGreater(len(menu), 100)
        for name, selection in menu.items():
            with self.subTest(selection=name):
                expanded = printer.expand_selection(selection, name)
                for section, value in selection.items():
                    if section == "evt_cuts":
                        self.assertEqual(expanded[section], flatten(value))
                    elif section in printer.SECTIONS:
                        for obj, cuts in value.items():
                            self.assertEqual(expanded[section][obj], flatten(cuts))

    def test_shared_alias_is_not_a_cycle(self):
        shared = ["pt"]
        self.assertEqual(printer.flatten_cuts([shared, shared]), ["pt", "pt"])

    def test_cycle_is_rejected(self):
        value = yaml.safe_load("&cycle [*cycle]")
        with self.assertRaisesRegex(ValueError, "recursive YAML alias"):
            printer.flatten_cuts(value)

    def test_invalid_leaf_is_rejected(self):
        with self.assertRaises(ValueError):
            printer.flatten_cuts([True])

    def test_comparison_distinguishes_reordering_from_changed_multiplicity(self):
        left = {"evt_cuts": ["trigger", "pv", "pv"]}
        reordered = {"evt_cuts": ["pv", "trigger", "pv"]}
        comparison = printer.compare_selections("a", left, "b", reordered)
        self.assertTrue(comparison["same_cuts"])
        self.assertFalse(comparison["same_order"])
        rows = comparison["groups"][0]["rows"]
        self.assertEqual([r["cut"] for r in rows if r["kind"] != "add"], left["evt_cuts"])
        self.assertEqual([r["cut"] for r in rows if r["kind"] != "remove"], reordered["evt_cuts"])
        fewer = printer.compare_selections("a", left, "b", {"evt_cuts": ["trigger", "pv"]})
        self.assertFalse(fewer["same_cuts"])

    def test_comparison_preserves_object_and_section_scope(self):
        left = {"obj_cuts": {"muons": ["pt"]}, "evt_cuts": ["trigger"]}
        right = {"obj_cuts": {"electrons": ["pt"]}, "postLj_obj_cuts": {"ljs": ["trigger"]}}
        result = printer.compare_selections("a", left, "b", right)
        self.assertFalse(result["same_cuts"])
        self.assertEqual([(g["section"], g["object"]) for g in result["groups"]],
                         [("obj_cuts", "muons"), ("obj_cuts", "electrons"),
                          ("postLj_obj_cuts", "ljs"), ("evt_cuts", None)])
        empty = printer.compare_selections("a", {}, "b", {"obj_cuts": {"muons": []}, "evt_cuts": []})
        self.assertTrue(empty["same_order"])

    def test_real_comparison_finds_added_channel_cut(self):
        menu = printer.read_selections(printer.DEFAULT_CONFIG)
        base = printer.expand_selection(menu["base_sr"], "base_sr")
        channel = printer.expand_selection(menu["4mu_sr"], "4mu_sr")
        result = printer.compare_selections("base_sr", base, "4mu_sr", channel)
        differences = [r for g in result["groups"] for r in g["rows"] if r["kind"] != "same"]
        self.assertEqual([(r["kind"], r["cut"]) for r in differences], [("add", "4mu")])

    def test_color_detection_and_plain_json(self):
        with patch.dict(printer.os.environ, {"TERM": "xterm"}, clear=True):
            with patch.object(printer.sys.stdout, "isatty", return_value=True):
                self.assertTrue(printer.use_color("auto"))
                with patch.dict(printer.os.environ, {"NO_COLOR": ""}):
                    self.assertFalse(printer.use_color("auto"))
                    self.assertTrue(printer.use_color("always"))
                self.assertFalse(printer.use_color("never"))
            with patch.object(printer.sys.stdout, "isatty", return_value=False):
                self.assertFalse(printer.use_color("auto"))
        for flags in [["-s", "4mu_sr"], ["--compare", "4mu_sr", "4mu_sr"]]:
            output = io.StringIO()
            with redirect_stdout(output):
                printer.main(flags + ["--json", "--color", "always"])
            self.assertNotIn("\033[", output.getvalue())
            result = json.loads(output.getvalue())
            if "--compare" in flags:
                self.assertTrue(result["same_order"])

    def test_comparison_cli_text_colors_and_errors(self):
        for mode in ["auto", "never", "always"]:
            output = io.StringIO()
            with redirect_stdout(output):
                printer.main(["--compare", "base_sr", "4mu_sr", "--color", mode])
            self.assertEqual("\033[" in output.getvalue(), mode == "always")
            self.assertIn("+ [-:6] 4mu", output.getvalue())
        for flags in [["--compare", "missing", "4mu_sr"],
                      ["--compare", "base", "4mu", "--list"],
                      ["--compare", "base", "4mu", "-s", "4mu"]]:
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                printer.main(flags)
            self.assertEqual(error.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
