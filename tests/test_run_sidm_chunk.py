"""Tests of the sample-metadata lookup in condor/run_sidm_chunk.py.

Run: python -m pytest tests/test_run_sidm_chunk.py
"""

import importlib.util
import sys
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("run_sidm_chunk", ROOT / "condor/run_sidm_chunk.py")
chunk = importlib.util.module_from_spec(spec)
spec.loader.exec_module(chunk)

SAMPLE = "DYJetsToLL_M10to50"
SKIM_DIR = "/store/group/test/Backgrounds/Skims/DY_skimmed/"
FULL_DIR = "/store/group/test/Backgrounds/DY_full/260630_045316/"


def write_yaml(path, cfg):
    path.write_text(yaml.safe_dump(cfg))
    return path


@pytest.fixture(autouse=True)
def ntuple_dir(tmp_path, monkeypatch):
    """Location YAMLs searched in place of sidm/configs/ntuples: a background listed in a
    skimmed and an unskimmed section, and a data config."""
    write_yaml(tmp_path / "backgrounds.yaml", {
        "skimmed": {
            "path": "root://xcache//store/group/test/Backgrounds/Skims/",
            "samples": {SAMPLE: {"path": "DY_skimmed/", "skim_factor": 0.0067,
                                 "files": ["skimmed_output_1.root"]}},
        },
        "unskimmed": {
            "path": "root://xcache//store/group/test/Backgrounds/",
            "samples": {SAMPLE: {"path": "DY_full/260630_045316/", "skim_factor": 1.0,
                                 "files": ["0000/output_1.root"]}},
        },
    })
    write_yaml(tmp_path / "data.yaml", {
        "data_2018": {
            "path": "root://xcache//store/group/test/Data/",
            "samples": {"DoubleMuon_2018A": {"path": "A/", "year": "2018"}},
        },
    })
    monkeypatch.setattr(chunk, "NTUPLE_DIR", tmp_path)


def test_skimmed_files_get_the_skim_factor_of_the_skimmed_section():
    files = [f"root://cmseos.fnal.gov/{SKIM_DIR}skimmed_output_{i}.root" for i in (1, 2)]
    metadata, source = chunk.sample_metadata(SAMPLE, files)
    assert source == "backgrounds.yaml:skimmed"
    assert metadata["skim_factor"] == pytest.approx(0.0067)
    assert metadata["is_data"] is False


def test_unskimmed_files_get_the_skim_factor_of_the_unskimmed_section():
    files = [f"root://cmseos.fnal.gov/{FULL_DIR}0000/output_{i}.root" for i in (1, 2)]
    metadata, source = chunk.sample_metadata(SAMPLE, files)
    assert source == "backgrounds.yaml:unskimmed"
    assert metadata["skim_factor"] == 1.0


def test_sample_of_a_data_config_is_data():
    files = ["root://cmseos.fnal.gov//store/group/test/Data/A/out_1.root"]
    metadata, source = chunk.sample_metadata("DoubleMuon_2018A", files)
    assert source == "data.yaml:data_2018"
    assert metadata["is_data"] is True
    assert metadata["year"] == "2018"


@pytest.mark.parametrize("prefix", [
    "root://cmsxrootd.fnal.gov/",
    "root://xcache/",
    "root://cmseos.fnal.gov//eos/uscms",
    "",
])
def test_any_xrootd_prefix_matches(prefix):
    files = [f"{prefix}{SKIM_DIR}skimmed_output_7.root"]
    _, source = chunk.sample_metadata(SAMPLE, files)
    assert source == "backgrounds.yaml:skimmed"


def test_store_path_strips_the_server_part():
    assert chunk.store_path("root://xcache//store/a/b.root") == "/store/a/b.root"
    assert chunk.store_path("root://cmseos.fnal.gov//store/a/b.root") == "/store/a/b.root"
    assert chunk.store_path("local/b.root") == "local/b.root"


def test_files_outside_every_listed_directory_have_no_metadata():
    files = ["root://cmseos.fnal.gov//store/group/elsewhere/DY/out_1.root"]
    assert chunk.sample_metadata(SAMPLE, files) == (None, None)


def test_unlisted_sample_has_no_metadata():
    files = [f"root://cmseos.fnal.gov/{SKIM_DIR}x.root"]
    assert chunk.sample_metadata("NotASample", files) == (None, None)
