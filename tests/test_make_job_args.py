"""job_args.txt line format of condor/make_job_args.py; no EOS access.

Run from the repository root:
    python -m pytest tests/test_make_job_args.py -v
"""

import importlib.util
import os

from sidm.tools import utilities

_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "condor", "make_job_args.py")
_SPEC = importlib.util.spec_from_file_location("make_job_args", _PATH)
make_job_args = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(make_job_args)


def test_line_carries_production():
    """The fourth column is the production tag, as submit.sub queues it."""
    line = make_job_args.job_args_line("2Mu2E_500GeV_1p2GeV_1p9mm", 3,
                                       "filelists/2Mu2E_500GeV_1p2GeV_1p9mm_3.txt", "v10")
    assert line == "2Mu2E_500GeV_1p2GeV_1p9mm 3 filelists/2Mu2E_500GeV_1p2GeV_1p9mm_3.txt v10\n"


def test_line_writes_none_for_untagged():
    """A location config without a production tag gives the literal 'none'."""
    production = utilities.get_production("backgrounds.yaml")
    assert production is None
    line = make_job_args.job_args_line("DYJetsToLL_M50", 0, "filelists/DYJetsToLL_M50_0.txt", production)
    assert line.split() == ["DYJetsToLL_M50", "0", "filelists/DYJetsToLL_M50_0.txt", "none"]
