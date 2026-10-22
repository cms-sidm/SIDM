#!/bin/bash
# Build correctionlib from source into this repo's .python_userbase, for use
# under cmsenv (the CMSSW_17_0_0_pre4 Jupyter kernel).
#
#     cmsenv                                   # inside CMSSW_17_0_0_pre4/src
#     bash sidm/scripts/install_correctionlib_cmssw.sh [version]
#
# Why this exists: the PyPI wheel cannot be imported under cmsenv.
#
#   ImportError: .../correctionlib/_core.cpython-312-x86_64-linux-gnu.so:
#       undefined symbol: __cxa_call_terminate
#
# The wheel's _core extension finds its bundled lib/libcorrectionlib.so through
# a RUNPATH. cmsenv puts CMSSW's own libcorrectionlib.so (2.7.0 in
# CMSSW_17_0_0_pre4) on LD_LIBRARY_PATH, and the dynamic loader searches
# LD_LIBRARY_PATH before RUNPATH, so _core binds the CMSSW copy and misses a
# symbol that only the bundled copy provides. Building here, with CMSSW's
# compiler and zlib, and linking with --disable-new-dtags writes a DT_RPATH
# instead, which the loader searches before LD_LIBRARY_PATH.
#
# sidm/tools/corrections.py works around the same problem at runtime (it
# preloads the bundled library before importing correctionlib), so the analysis
# code also runs on a stock wheel, e.g. on Condor. This script is what makes a
# plain `import correctionlib` work in a notebook.
#
# The default version follows constraints.txt, so a later
# `pip install -c constraints.txt ...` sees the requirement already satisfied
# and leaves this build alone.
set -euo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
DEFAULT_VERSION=$(sed -n 's/^correctionlib==//p' "$REPO/constraints.txt")
VERSION=${1:-$DEFAULT_VERSION}

if [ -z "${CMSSW_BASE:-}" ]; then
    echo "ERROR: run cmsenv first (this build must use CMSSW's compiler and zlib)" >&2
    exit 1
fi
PYVER=$(python3 -c 'import sys; print("%d.%d" % sys.version_info[:2])')
if [ "$PYVER" != "3.12" ]; then
    echo "ERROR: python3 is $PYVER; the SIDM kernel uses CMSSW_17_0_0_pre4's Python 3.12" >&2
    exit 1
fi

ZLIB_BASE=$(cd "$CMSSW_BASE/src" && scram tool info zlib | sed -n 's/^ZLIB_BASE=//p')
if [ -z "$ZLIB_BASE" ]; then
    echo "ERROR: could not find CMSSW's zlib with 'scram tool info zlib'" >&2
    exit 1
fi

export PYTHONUSERBASE="$REPO/.python_userbase"
RPATH_FLAG="-Wl,--disable-new-dtags"
export CMAKE_ARGS="-DZLIB_ROOT=$ZLIB_BASE -DCMAKE_MODULE_LINKER_FLAGS=$RPATH_FLAG -DCMAKE_SHARED_LINKER_FLAGS=$RPATH_FLAG"

echo "Building correctionlib $VERSION into $PYTHONUSERBASE"
echo "  compiler: $(command -v c++) ($(c++ -dumpfullversion))"
echo "  zlib:     $ZLIB_BASE"

# --no-cache-dir: a wheel cached from an earlier build with other flags would
# otherwise be reused as-is.
python3 -m pip install --user --no-deps --force-reinstall --no-cache-dir \
    --no-binary correctionlib "correctionlib==$VERSION"

SITE="$PYTHONUSERBASE/lib/python3.12/site-packages"
CORE=$(ls "$SITE"/correctionlib/_core*.so)
if ! readelf -d "$CORE" | grep -q '(RPATH)'; then
    echo "ERROR: $CORE was linked with a RUNPATH, not an RPATH; the import will still fail under cmsenv" >&2
    exit 1
fi

# Import from outside the repo so the check exercises the installed copy.
(cd / && PYTHONPATH="$SITE${PYTHONPATH:+:$PYTHONPATH}" python3 - "$VERSION" <<'PY'
import sys
import correctionlib
import correctionlib._core
assert correctionlib.__version__ == sys.argv[1], (correctionlib.__version__, sys.argv[1])
print(f"OK: correctionlib {correctionlib.__version__} from {correctionlib.__file__}")
PY
)
