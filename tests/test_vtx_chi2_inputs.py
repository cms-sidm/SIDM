#!/usr/bin/env python3
"""The histograms that read on-demand inputs are listed explicitly.

    python tests/test_vtx_chi2_inputs.py

SidmProcessor loads dsaMuonsRaw and the three dimuon-vertex tables, and computes
mu_ljs.vtx_chi2, only when a requested channel cuts on vtx_chi2 or a requested histogram is
in sidm.definitions.hists.vtx_chi2_or_raw_dsa_hists; likewise it loads the uncut photon and
muon copies of the lepton-jet control samples only for a channel with a probe cut or a
histogram in probe_hists. A histogram that reads such an input but is missing from its set
would be filled without it. This test finds every such histogram from the source of
hists.py: a definition that names the attribute or key directly, or calls a module-level
helper that does (helpers calling helpers included).
"""
import ast
import os
import sys

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)
HISTS_PY = os.path.join(_REPO, "sidm", "definitions", "hists.py")
TOKENS = frozenset({"vtx_chi2", "dsaMuonsRaw"})
# the probe objects of the lepton-jet control samples read the uncut copies photonsRaw and muonsRaw
PROBE_TOKENS = frozenset({"photonsRaw", "muonsRaw", "zmmg_probe_photons", "jpsi_probe_dimuons",
                          "dsa_probe_muons"})


def _reads_token(node, readers=(), tokens=TOKENS):
    """Does this AST node name a token (attribute, or string constant equal to it), or call one
    of the helper functions in readers?"""
    for sub in ast.walk(node):
        if isinstance(sub, ast.Attribute) and sub.attr in tokens:
            return True
        if isinstance(sub, ast.Constant) and sub.value in tokens:
            return True
        if isinstance(sub, ast.Name) and sub.id in readers:
            return True
    return False


def hists_reading_tokens(path=HISTS_PY, tokens=TOKENS):
    """Names of the hist_defs entries that read one of tokens (by default vtx_chi2 or dsaMuonsRaw)."""
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read())
    funcs = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    readers = {name for name, fn in funcs.items() if _reads_token(fn, tokens=tokens)}
    grew = True
    while grew:
        new = {name for name, fn in funcs.items()
               if name not in readers and _reads_token(fn, readers, tokens)}
        readers |= new
        grew = bool(new)
    hist_dict = next(n.value for n in tree.body if isinstance(n, ast.Assign)
                     and any(isinstance(t, ast.Name) and t.id == "hist_defs" for t in n.targets))
    return {k.value for k, v in zip(hist_dict.keys, hist_dict.values)
            if _reads_token(v, readers, tokens)}


def test_every_reader_is_listed():
    from sidm.definitions.hists import vtx_chi2_or_raw_dsa_hists
    missing = hists_reading_tokens() - set(vtx_chi2_or_raw_dsa_hists)
    assert not missing, f"read vtx_chi2 or dsaMuonsRaw but not in vtx_chi2_or_raw_dsa_hists: {sorted(missing)}"


def test_listed_names_are_histograms():
    from sidm.definitions.hists import hist_defs, vtx_chi2_or_raw_dsa_hists
    unknown = set(vtx_chi2_or_raw_dsa_hists) - set(hist_defs)
    assert not unknown, f"not in hist_defs: {sorted(unknown)}"


def test_every_probe_reader_is_listed():
    """Every histogram that reads the probe objects is in probe_hists."""
    from sidm.definitions.hists import probe_hists
    missing = hists_reading_tokens(tokens=PROBE_TOKENS) - set(probe_hists)
    assert not missing, f"read the probe objects but not in probe_hists: {sorted(missing)}"


def test_listed_probe_names_are_histograms():
    """Every name in probe_hists is a histogram."""
    from sidm.definitions.hists import hist_defs, probe_hists
    unknown = set(probe_hists) - set(hist_defs)
    assert not unknown, f"not in hist_defs: {sorted(unknown)}"


if __name__ == "__main__":
    test_every_reader_is_listed()
    test_listed_names_are_histograms()
    test_every_probe_reader_is_listed()
    test_listed_probe_names_are_histograms()
    print("4 passed")
