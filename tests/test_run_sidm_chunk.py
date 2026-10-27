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
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(cfg))
    return path


@pytest.fixture(name="configs")
def fixture_configs(tmp_path):
    """A primary YAML listing the sample in a skimmed and an unskimmed section, and an
    ntuple directory holding a data config that does not list it."""
    primary = write_yaml(tmp_path / "ntuples" / "backgrounds.yaml", {
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
    write_yaml(tmp_path / "ntuples" / "data.yaml", {
        "data_2018": {
            "path": "root://xcache//store/group/test/Data/",
            "samples": {"DoubleMuon_2018A": {"path": "A/", "is_data": True, "year": "2018"}},
        },
    })
    return primary, tmp_path / "ntuples"


def resolve(configs, files, sample=SAMPLE, section="skimmed"):
    primary, ntuple_dir = configs
    return chunk.find_sample_cfg(sample, files, primary, section, ntuple_dir=ntuple_dir)


def test_skimmed_files_resolve_to_the_skimmed_section(configs):
    files = [f"root://cmseos.fnal.gov/{SKIM_DIR}skimmed_output_{i}.root" for i in (1, 2)]
    entry, _, section = resolve(configs, files)
    assert section == "skimmed"
    assert entry["skim_factor"] == pytest.approx(0.0067)


def test_unskimmed_files_resolve_to_the_unskimmed_section(configs):
    files = [f"root://cmseos.fnal.gov/{FULL_DIR}0000/output_{i}.root" for i in (1, 2)]
    entry, _, section = resolve(configs, files)
    assert section == "unskimmed"
    assert entry["skim_factor"] == 1.0


def test_primary_section_does_not_override_the_files_directory(configs):
    files = [f"root://cmseos.fnal.gov/{FULL_DIR}0000/output_1.root"]
    _, _, section = resolve(configs, files, section="skimmed")
    assert section == "unskimmed"


def test_data_sample_resolves_from_the_ntuple_directory(configs):
    files = ["root://cmseos.fnal.gov//store/group/test/Data/A/out_1.root"]
    entry, yaml_path, section = resolve(configs, files, sample="DoubleMuon_2018A")
    assert (Path(yaml_path).name, section) == ("data.yaml", "data_2018")
    assert entry["is_data"] is True
    assert entry["year"] == "2018"


@pytest.mark.parametrize("prefix", [
    "root://cmsxrootd.fnal.gov/",
    "root://xcache/",
    "root://cmseos.fnal.gov//eos/uscms",
    "",
])
def test_any_xrootd_prefix_matches(configs, prefix):
    files = [f"{prefix}{SKIM_DIR}skimmed_output_7.root"]
    _, _, section = resolve(configs, files)
    assert section == "skimmed"


def test_store_path_strips_the_server_part():
    assert chunk._store_path("root://xcache//store/a/b.root") == "/store/a/b.root"  # pylint: disable=protected-access
    assert chunk._store_path("root://cmseos.fnal.gov//store/a/b.root") == "/store/a/b.root"  # pylint: disable=protected-access
    assert chunk._store_path("local/b.root") == "local/b.root"  # pylint: disable=protected-access


def test_no_matching_section_raises_naming_the_candidates(configs):
    files = ["root://cmseos.fnal.gov//store/group/elsewhere/DY/out_1.root"]
    with pytest.raises(KeyError) as error:
        resolve(configs, files)
    message = str(error.value)
    assert "no section" in message
    assert "backgrounds.yaml:skimmed -> " + SKIM_DIR in message
    assert "backgrounds.yaml:unskimmed -> " + FULL_DIR in message


def test_unlisted_sample_raises(configs):
    with pytest.raises(KeyError, match="Candidates: none"):
        resolve(configs, [f"root://cmseos.fnal.gov/{SKIM_DIR}x.root"], sample="NotASample")


def test_files_spanning_two_sections_raise(configs):
    files = [f"root://cmseos.fnal.gov/{SKIM_DIR}skimmed_output_1.root",
             f"root://cmseos.fnal.gov/{FULL_DIR}0000/output_1.root"]
    with pytest.raises(KeyError, match="no section"):
        resolve(configs, files)


def test_two_matching_sections_raise(configs):
    _, ntuple_dir = configs
    write_yaml(ntuple_dir / "backgrounds_copy.yaml", {
        "skimmed_again": {
            "path": "root://xcache//store/group/test/Backgrounds/Skims/",
            "samples": {SAMPLE: {"path": "DY_skimmed/", "skim_factor": 0.5}},
        },
    })
    files = [f"root://cmseos.fnal.gov/{SKIM_DIR}skimmed_output_1.root"]
    with pytest.raises(KeyError) as error:
        resolve(configs, files)
    message = str(error.value)
    assert "2 sections" in message
    assert "backgrounds_copy.yaml:skimmed_again" in message
