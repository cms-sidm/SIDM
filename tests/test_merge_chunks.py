"""Chunk merge: Runner-wrapped chunks are re-normalized, and inconsistent chunks refuse to merge.

Run from the repository root:
    python -m pytest tests/test_merge_chunks.py -v
"""

import hist
import pytest
from coffea import processor

from sidm.scripts.merge_coffea_chunks_eos import combine_normalized_chunks, sample_outputs

SAMPLE = "2Mu2E_500GeV_1p2GeV_1p9mm"


def chunk(sumw, eff=0.7594, production="v10", wrapped=True):
    """One chunk as postprocess leaves it: a weighted hist already scaled to lumi*xs*eff/sumw.
    With unit raw content per chunk, the normalized value is 1/sumw (taking lumi*xs*eff = 1)."""
    h = hist.Hist.new.Reg(1, 0, 1).Weight()
    h.fill([0.5], weight=[1.0 / sumw])
    out = {SAMPLE: {
        "metadata": {
            "is_data": processor.set_accumulator([False]),
            "year": processor.set_accumulator(["2018"]),
            "unweighted_hist": processor.set_accumulator([False]),
            "production": processor.set_accumulator([production]),
            "gen_filter_eff": processor.set_accumulator([eff]),
            "scaled_sum_weights": sumw,
        },
        "cutflow": {},
        "hists": {"x": h},
    }}
    if not wrapped:
        return out
    # the Runner wrapper of a run that succeeded
    return {"out": out, "processed": set(), "exception": 0}


@pytest.mark.parametrize("wrapped", [True, False])
def test_merge_normalizes_to_full_sample(wrapped):
    """Two chunks of sumw 10 and 30, one raw entry each: lumi*xs*eff * 2 / 40 = 0.05."""
    merged = sample_outputs(combine_normalized_chunks(
        [chunk(10.0, wrapped=wrapped),
         chunk(30.0, wrapped=wrapped)]))
    assert merged[SAMPLE]["hists"]["x"].values()[0] == pytest.approx(0.05)
    assert merged[SAMPLE]["metadata"]["scaled_sum_weights"] == pytest.approx(40.0)


@pytest.mark.parametrize("kwargs", [{"eff": 0.5}, {"production": None}])
def test_merge_refuses_inconsistent_chunks(kwargs):
    """Chunks of one sample normalized with different settings raise instead of summing."""
    with pytest.raises(ValueError, match="chunks disagree"):
        combine_normalized_chunks([chunk(10.0), chunk(30.0, **kwargs)])
