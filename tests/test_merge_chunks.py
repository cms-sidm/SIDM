"""Chunk merge: Runner-wrapped chunks are re-normalized, and inconsistent chunks refuse to merge.

Run from the repository root:
    python -m pytest tests/test_merge_chunks.py -v
"""

import hist
import pytest
from coffea import processor

from sidm.scripts.merge_coffea_chunks_eos import combine_normalized_chunks, sample_outputs
from sidm.tools.cutflow import SimpleCutflow

SAMPLE = "2Mu2E_500GeV_1p2GeV_1p9mm"


def chunk(sumw, eff=0.7594, production="v10", wrapped=True, unweighted=False):
    """One chunk as postprocess leaves it: a weighted hist already scaled to lumi*xs*eff/sumw.
    With unit raw content per chunk, the normalized value is 1/sumw (taking lumi*xs*eff = 1)."""
    h = hist.Hist.new.Reg(1, 0, 1).Weight()
    h.fill([0.5], weight=[1.0 if unweighted else 1.0 / sumw])
    cutflow = SimpleCutflow()
    cutflow.add_row("None", raw=1, weighted=1.0 / sumw)
    out = {SAMPLE: {
        "metadata": {
            "is_data": processor.set_accumulator([False]),
            "year": processor.set_accumulator(["2018"]),
            "unweighted_hist": processor.set_accumulator([unweighted]),
            "production": processor.set_accumulator([production]),
            "gen_filter_eff": processor.set_accumulator([eff]),
            "scaled_sum_weights": sumw,
        },
        "cutflow": {"base": cutflow},
        "hists": {"x": h},
    }}
    if not wrapped:
        return out
    # the Runner wrapper of a run that succeeded
    return {"out": out, "processed": set(), "exception": 0}


@pytest.mark.parametrize("wrapped", [True, False])
@pytest.mark.parametrize("unweighted", [False, True])
def test_merge_normalizes_to_full_sample(wrapped, unweighted):
    """Two chunks of sumw 10 and 30, one raw entry each: lumi*xs*eff * 2 / 40 = 0.05."""
    merged = sample_outputs(combine_normalized_chunks(
        [chunk(10.0, wrapped=wrapped, unweighted=unweighted),
         chunk(30.0, wrapped=wrapped, unweighted=unweighted)]))
    h = merged[SAMPLE]["hists"]["x"]
    assert h.values()[0] == pytest.approx(2.0 if unweighted else 0.05)
    assert h.variances()[0] == pytest.approx(2.0 if unweighted else 2.0 / 40.0**2)
    row = merged[SAMPLE]["cutflow"]["base"].rows["None"]
    assert row["raw"] == 2
    assert row["weighted"] == pytest.approx(0.05)
    assert merged[SAMPLE]["metadata"]["scaled_sum_weights"] == pytest.approx(40.0)


@pytest.mark.parametrize("kwargs", [{"eff": 0.5}, {"production": None}])
def test_merge_refuses_inconsistent_chunks(kwargs):
    """Chunks of one sample normalized with different settings raise instead of summing."""
    with pytest.raises(ValueError, match="chunks disagree"):
        combine_normalized_chunks([chunk(10.0), chunk(30.0, **kwargs)])


@pytest.mark.parametrize("missing", [("gen_filter_eff",), ("production",),
                                    ("gen_filter_eff", "production")])
@pytest.mark.parametrize("legacy_first", [False, True])
def test_merge_refuses_mixed_legacy_metadata_before_scaling(missing, legacy_first):
    """A missing key must not disappear into the other chunk's recorded value."""
    current = chunk(10.0)
    legacy = chunk(30.0)
    for key in missing:
        del legacy["out"][SAMPLE]["metadata"][key]
    outputs = [legacy, current] if legacy_first else [current, legacy]
    with pytest.raises(ValueError, match="chunks disagree"):
        combine_normalized_chunks(outputs)
    for output, sumw in [(current, 10.0), (legacy, 30.0)]:
        sample = output["out"][SAMPLE]
        assert sample["cutflow"]["base"].rows["None"]["weighted"] == pytest.approx(1.0 / sumw)
        assert sample["hists"]["x"].values()[0] == pytest.approx(1.0 / sumw)


def test_merge_accepts_all_legacy_chunks():
    """Keep old campaigns usable when every chunk lacks the normalization keys."""
    outputs = [chunk(10.0), chunk(30.0)]
    for output in outputs:
        for key in ("gen_filter_eff", "production"):
            del output["out"][SAMPLE]["metadata"][key]
    merged = sample_outputs(combine_normalized_chunks(outputs))
    assert merged[SAMPLE]["cutflow"]["base"].rows["None"]["weighted"] == pytest.approx(0.05)
    assert "gen_filter_eff" not in merged[SAMPLE]["metadata"]
    assert "production" not in merged[SAMPLE]["metadata"]


def test_normalization_settings_are_checked_per_sample():
    """Different samples may use different production tags and efficiencies."""
    first = chunk(10.0)
    second = chunk(30.0, eff=1.0, production=None)
    second["out"]["other"] = second["out"].pop(SAMPLE)
    merged = sample_outputs(combine_normalized_chunks([first, second]))
    assert merged[SAMPLE]["metadata"]["gen_filter_eff"] == {0.7594}
    assert merged["other"]["metadata"]["gen_filter_eff"] == {1.0}


@pytest.mark.parametrize("values", [None, set(), {0.5, 1.0}])
def test_merge_refuses_ambiguous_normalization_in_one_chunk(values):
    output = chunk(10.0)
    output["out"][SAMPLE]["metadata"]["gen_filter_eff"] = values
    with pytest.raises(ValueError, match="exactly one value"):
        combine_normalized_chunks([output])
