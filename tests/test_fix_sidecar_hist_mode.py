"""Sidecar repair: only the hist mode of a sidecar that disagrees with its output changes.

Run from the repository root:
    python -m pytest tests/test_fix_sidecar_hist_mode.py -v
"""

import shutil
import sys

import pytest
from coffea import processor
from coffea.util import save

from sidm.scripts import fix_sidecar_hist_mode

SIDECAR = """created_utc: '2026-01-01T00:00:00Z'
unweighted_hist: {mode}
samples:
- name: s
  unweighted_hist: kept
"""


def output(unweighted):
    """A merged output that records its hist mode, wrapped as the Runner leaves it"""
    metadata = {"unweighted_hist": processor.set_accumulator([unweighted])}
    return {"out": {"s": {"metadata": metadata}}, "processed": set(), "exception": 0}


@pytest.fixture(name="eos")
def fixture_eos(tmp_path, monkeypatch):
    """A local directory standing in for EOS"""
    monkeypatch.setattr(fix_sidecar_hist_mode, "eos_find",
                        lambda eos_dir: sorted(str(path) for path in tmp_path.rglob("*.*")))
    monkeypatch.setattr(fix_sidecar_hist_mode, "xrdcp_from_eos", shutil.copy)
    monkeypatch.setattr(fix_sidecar_hist_mode, "xrdcp_to_eos", shutil.copy)
    # wrong: the output is unweighted, the sidecar says weighted
    save(output(True), str(tmp_path / "wrong.coffea"))
    (tmp_path / "wrong.meta.yaml").write_text(SIDECAR.format(mode="false"))
    save(output(False), str(tmp_path / "right.coffea"))
    (tmp_path / "right.meta.yaml").write_text(SIDECAR.format(mode="false"))
    (tmp_path / "alone.meta.yaml").write_text(SIDECAR.format(mode="false"))
    # a sidecar that states no mode is not a wrong one
    save(output(True), str(tmp_path / "unstated.coffea"))
    (tmp_path / "unstated.meta.yaml").write_text(SIDECAR.format(mode="null"))
    return tmp_path


def run(eos, monkeypatch, *flags):
    """Run the script on the local directory and return the text of each sidecar"""
    monkeypatch.setattr(sys, "argv", ["fix_sidecar_hist_mode.py", "--eos-dir", str(eos), *flags])
    fix_sidecar_hist_mode.main()
    return {name: (eos / f"{name}.meta.yaml").read_text()
            for name in ("wrong", "right", "alone", "unstated")}


def test_report_only_by_default(eos, monkeypatch, capsys):
    """Without --apply the wrong sidecar is named and nothing is written"""
    before = {path: path.read_bytes() for path in eos.iterdir()}
    run(eos, monkeypatch)
    assert {path: path.read_bytes() for path in eos.iterdir()} == before
    printed = capsys.readouterr().out
    assert "WRONG: " + str(eos / "wrong.meta.yaml") in printed
    assert "SIDECAR STATES NO MODE: " + str(eos / "unstated.meta.yaml") in printed
    assert "WRONG: " + str(eos / "right.meta.yaml") not in printed


def test_apply_rewrites_only_the_field_of_the_wrong_sidecar(eos, monkeypatch):
    """With --apply the wrong sidecar gets the output's mode, every other byte stays"""
    texts = run(eos, monkeypatch, "--apply")
    assert texts["wrong"] == SIDECAR.format(mode="true")
    assert texts["right"] == SIDECAR.format(mode="false")
    assert texts["alone"] == SIDECAR.format(mode="false")
    assert texts["unstated"] == SIDECAR.format(mode="null")


def test_fill_missing_writes_the_mode_into_the_sidecar_that_states_none(eos, monkeypatch):
    """With --fill-missing the sidecar with a null field gets the output's mode, every
    other sidecar stays byte for byte"""
    before = {name: (eos / f"{name}.meta.yaml").read_bytes()
              for name in ("wrong", "right", "alone")}
    texts = run(eos, monkeypatch, "--fill-missing")
    assert texts["unstated"] == SIDECAR.format(mode="true")
    assert {name: (eos / f"{name}.meta.yaml").read_bytes() for name in before} == before
