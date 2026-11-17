#!/bin/bash
# Condor job of the dR truth match: builds the worker venv (CMSSW_14_1_0_pre4 python plus
# condor/requirements.txt and condor/constraints.txt of the repository, as condor/run_job.sh
# does), unpacks the SIDM code and writes the per-LJ records of one filelist with extract.py.
# usage (by submit.sub): run.sh SAMPLE CHUNK FILELIST
set -euo pipefail
set -x
SAMPLE=$1; CHUNK=$2; FILELIST=$3
cd "${_CONDOR_SCRATCH_DIR}"
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
scramv1 project CMSSW CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src; eval `scramv1 runtime -sh`; cd "${_CONDOR_SCRATCH_DIR}"
python3 -m venv sidm_venv --system-site-packages
source sidm_venv/bin/activate
export PIP_CACHE_DIR="${_CONDOR_SCRATCH_DIR}/pip_cache"
python -m pip install --upgrade pip setuptools wheel
python -m pip install --no-cache-dir -r requirements.txt -c constraints.txt
tar -xzf sidm_code.tar.gz
export PYTHONPATH="${_CONDOR_SCRATCH_DIR}:${PYTHONPATH:-}"
python extract.py "${SAMPLE}" "$(basename ${FILELIST})" "${SAMPLE}_${CHUNK}.parquet" \
  --chunksize 50000 --workers 1
ls -l "${SAMPLE}_${CHUNK}.parquet"
