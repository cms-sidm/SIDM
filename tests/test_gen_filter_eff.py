"""Checks of the v10 signal gen-filter efficiency table and its use in the
normalization; reads only repo files.

Run from the repository root:
    python -m pytest tests/test_gen_filter_eff.py -v
"""

import os

import numpy as np
import pytest
import yaml

from sidm import BASE_DIR
from sidm.scripts import make_signal_gen_filter_eff as make_table
from sidm.tools import utilities

V10_NTUPLES = ["signal_2mu2e_v10.yaml", "signal_4mu_v10.yaml"]
# Median re-pass of the 5 GeV cut string by mXX, from the merged
# truth_kinematics_forAN outputs (genFilterEmulation channel).
REPASS_MEDIANS = {100: 0.847, 150: 0.898, 200: 0.923, 500: 0.969, 800: 0.980, 1000: 0.985}


def _load(path):
    with open(path, encoding="utf8") as f:
        return yaml.safe_load(f)


def _samples(ntuple_cfgs):
    names = set()
    for cfg in ntuple_cfgs:
        for version in _load(f"{BASE_DIR}/configs/ntuples/{cfg}").values():
            names |= set(version["samples"])
    return names


TABLE = _load(f"{BASE_DIR}/configs/signal_gen_filter_eff.yaml")


def test_table_covers_exactly_the_v10_samples():
    """Every v10 sample is in the table and nothing else is."""
    v10 = _samples(V10_NTUPLES)
    assert len(v10) == 180
    assert set(TABLE) == v10


def test_efficiencies_ordered():
    """0 < central <= v10 <= 1 for every sample."""
    for name, entry in TABLE.items():
        assert 0 < entry["central"]["eff"] <= entry["v10"]["eff"] <= 1, name


def test_repass_medians_by_mxx():
    """Re-pass medians by mXX match the merged outputs to 3 decimals."""
    for mxx, median in REPASS_MEDIANS.items():
        vals = [e["repass_5gev"]["n_pass"] / e["repass_5gev"]["n_all"]
                for name, e in TABLE.items() if name.split("_")[1] == f"{mxx}GeV"]
        assert round(float(np.median(vals)), 3) == median, mxx


def test_table_regenerates_from_inputs():
    """The committed table equals a fresh build from the committed inputs."""
    table, missing = make_table.build_from_inputs()
    assert not missing
    assert table == TABLE


def test_name_parsing_by_value():
    """'p' is the decimal point; central and sample names agree on values, not strings."""
    assert make_table.sample_point("2Mu2E_1000GeV_0p25GeV_0p002mm") == ("2Mu2E", 1000.0, 0.25, 0.002)
    assert make_table.sample_point("4Mu_500GeV_1p2GeV_19p0mm") == ("4Mu", 500.0, 1.2, 19.0)
    assert make_table.sample_point("4Mu_1000GeV_1p2GeV_0p0096mm") == ("4Mu", 1000.0, 1.2, 0.0096)
    central = make_table.central_by_point(os.path.join(make_table.REPO, make_table.CENTRAL))
    # mA-5 in the central key, 5p0GeV in ours; ctau-19 vs 19p0mm
    point = make_table.sample_point("4Mu_500GeV_5p0GeV_80p0mm")
    assert point == ("4Mu", 500.0, 5.0, 80.0) and point in central
    assert make_table.sample_point("4Mu_500GeV_1p2GeV_19p0mm") in central


def test_production_tag_selects_efficiency(capsys):
    """v10 -> v10.eff, central -> central.eff, None or unknown -> 1.0 with one warning."""
    name = "4Mu_500GeV_1p2GeV_19p0mm"
    assert utilities.get_production("signal_4mu_v10.yaml") == "v10"
    assert utilities.get_production("ffntuples.yaml") == "central"
    assert utilities.get_production("signal_v8.yaml") is None
    assert utilities.get_gen_filter_eff(name, "v10") == TABLE[name]["v10"]["eff"]
    assert utilities.get_gen_filter_eff(name, "central") == TABLE[name]["central"]["eff"]
    capsys.readouterr()
    assert utilities.get_gen_filter_eff(name) == 1.0
    assert name in capsys.readouterr().out
    assert utilities.get_gen_filter_eff(name, "v8") == 1.0
    assert capsys.readouterr().out == ""
    assert utilities.get_gen_filter_eff("DYJetsToLL_M50") == 1.0


def test_make_fileset_sets_production():
    """make_fileset tags each sample with its location config's production."""
    for cfg, tag in [("signal_2mu2e_v10.yaml", "v10"), ("signal_v8.yaml", None)]:
        version = next(iter(_load(f"{BASE_DIR}/configs/ntuples/{cfg}")))
        sample = next(iter(_load(f"{BASE_DIR}/configs/ntuples/{cfg}")[version]["samples"]))
        fileset = utilities.make_fileset([sample], version, max_files=1, location_cfg=cfg)
        assert fileset[sample]["metadata"]["production"] == tag


def test_lumixs_weight_includes_gen_filter_eff():
    """Signal weights carry the table efficiency; backgrounds carry none."""
    sum_weights = 12345.6
    lumi = utilities.get_lumi("2018")

    signal = "2Mu2E_500GeV_1p2GeV_1p9mm"
    eff = TABLE[signal]["v10"]["eff"]
    assert eff < 1
    expected = lumi * utilities.get_xs(signal) * eff / sum_weights
    assert utilities.get_lumixs_weight(signal, "2018", sum_weights, "v10") == pytest.approx(expected,
                                                                                         rel=1e-12)

    background = "DYJetsToLL_M50"
    expected = lumi * utilities.get_xs(background) / sum_weights
    assert utilities.get_lumixs_weight(background, "2018", sum_weights) == pytest.approx(expected, rel=1e-12)


def test_unlisted_signal_falls_back_and_warns_once(capsys):
    """An unlisted signal dataset gets 1.0 and one warning."""
    name = "2Mu2E_500GeV_1p2GeV_123p0mm"
    assert name not in TABLE
    assert utilities.get_gen_filter_eff(name, "v10") == 1.0
    assert name in capsys.readouterr().out
    assert utilities.get_gen_filter_eff(name, "v10") == 1.0
    assert capsys.readouterr().out == ""
