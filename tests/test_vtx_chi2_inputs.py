#!/usr/bin/env python3
"""The histograms that read mu_ljs.vtx_chi2 or objs["dsaMuonsRaw"] are listed explicitly.

    python tests/test_vtx_chi2_inputs.py

SidmProcessor loads dsaMuonsRaw and the three dimuon-vertex tables, and computes
mu_ljs.vtx_chi2, only when a requested channel cuts on vtx_chi2 or a requested histogram is
in sidm.definitions.hists.vtx_chi2_or_raw_dsa_hists. A histogram that reads either input but
is missing from that set would be filled without it. This test finds every such histogram
from the source of hists.py: a definition that names the attribute or key directly, or calls
a module-level helper that does (helpers calling helpers included).
"""
import ast
import os
import sys

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)
HISTS_PY = os.path.join(_REPO, "sidm", "definitions", "hists.py")
TOKENS = {"vtx_chi2", "dsaMuonsRaw"}


def _reads_token(node, readers=()):
    """Does this AST node name a token (attribute, or string constant equal to it), or call one
    of the helper functions in readers?"""
    for sub in ast.walk(node):
        if isinstance(sub, ast.Attribute) and sub.attr in TOKENS:
            return True
        if isinstance(sub, ast.Constant) and sub.value in TOKENS:
            return True
        if isinstance(sub, ast.Name) and sub.id in readers:
            return True
    return False


def hists_reading_tokens(path=HISTS_PY):
    """Names of the hist_defs entries that read vtx_chi2 or dsaMuonsRaw."""
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read())
    funcs = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    readers = {name for name, fn in funcs.items() if _reads_token(fn)}
    grew = True
    while grew:
        new = {name for name, fn in funcs.items()
               if name not in readers and _reads_token(fn, readers)}
        readers |= new
        grew = bool(new)
    hist_dict = next(n.value for n in tree.body if isinstance(n, ast.Assign)
                     and any(isinstance(t, ast.Name) and t.id == "hist_defs" for t in n.targets))
    return {k.value for k, v in zip(hist_dict.keys, hist_dict.values) if _reads_token(v, readers)}


def test_every_reader_is_listed():
    from sidm.definitions.hists import vtx_chi2_or_raw_dsa_hists
    missing = hists_reading_tokens() - set(vtx_chi2_or_raw_dsa_hists)
    assert not missing, f"read vtx_chi2 or dsaMuonsRaw but not in vtx_chi2_or_raw_dsa_hists: {sorted(missing)}"


def test_listed_names_are_histograms():
    from sidm.definitions.hists import hist_defs, vtx_chi2_or_raw_dsa_hists
    unknown = set(vtx_chi2_or_raw_dsa_hists) - set(hist_defs)
    assert not unknown, f"not in hist_defs: {sorted(unknown)}"


if __name__ == "__main__":
    test_every_reader_is_listed()
    test_listed_names_are_histograms()
    print("2 passed")
