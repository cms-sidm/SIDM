#!/bin/bash
# Condor worker bootstrap for a file-census chunk. Mirrors run_job.sh: set up CMSSW +
# a throwaway venv from requirements.txt, unpack the SIDM code, probe the chunk's files
# with run_census_chunk.py, and xrdcp the JSON shard to EOS. The full stack is installed
# so the same script supports deep mode (full SidmProcessor); shallow mode only needs
# uproot but big chunks amortize the setup cost.
set -euo pipefail
set -x

CHUNK=$1
FILELIST=$2
EOS_OUTDIR=$3
MODE=${4:-shallow}

hostname; date
cd "${_CONDOR_SCRATCH_DIR}"
ls -lah

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
scramv1 project CMSSW CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
cd "${_CONDOR_SCRATCH_DIR}"

python3 -m venv sidm_venv --system-site-packages
source sidm_venv/bin/activate
export PIP_CACHE_DIR="${_CONDOR_SCRATCH_DIR}/pip_cache"
python -m pip install --upgrade pip setuptools wheel
python -m pip install --no-cache-dir -r requirements.txt

tar -xzf sidm_code.tar.gz
export PYTHONPATH="${_CONDOR_SCRATCH_DIR}:${PYTHONPATH:-}"

FILELIST_BASENAME=$(basename "${FILELIST}")
OUTFILE="census_${CHUNK}.json"
DEEPFLAG=""
if [[ "${MODE}" == "deep" ]]; then DEEPFLAG="--deep"; fi

python condor/run_census_chunk.py \
    --filelist "${FILELIST_BASENAME}" \
    --output "${OUTFILE}" \
    ${DEEPFLAG}

echo "Copying shard to EOS: ${EOS_OUTDIR}/${OUTFILE}"
xrdcp -f "${OUTFILE}" "${EOS_OUTDIR}/${OUTFILE}"
rm -f "${OUTFILE}"
date
echo "Done"
