"""Tests requiring only the standard library and PyYAML.

Run: python3 -m unittest discover -s tests -p test_print_cuts.py
"""

import ast
import importlib.util
from pathlib import Path
import unittest

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


if __name__ == "__main__":
    unittest.main()
