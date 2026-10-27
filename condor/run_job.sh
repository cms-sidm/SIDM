#!/bin/bash
set -euo pipefail
set -x

SAMPLE=$1
CHUNK=$2
FILELIST=$3
EOS_OUTDIR=$4
# Optional positional args 5/6: comma-separated channels + hist collections.
# Default to the historical "base" / "muon_base" so older submit.sub files
# that don't pass them keep working.
CHANNELS=${5:-base}
HIST_COLLECTIONS=${6:-muon_base}
# Optional positional args 7/8: an ntuple location YAML and a section in it, passed on to
# run_sidm_chunk.py. They do not choose the sample metadata (is_data, skim_factor, year):
# that comes from the one section, of this YAML or of the configs in sidm/configs/ntuples/,
# whose directory holds the chunk's files. The YAML is searched first and the section only
# orders the candidates named in the error raised when no section, or more than one, matches.
NTUPLE_YAML=${7:-sidm/configs/ntuples/backgrounds.yaml}
NTUPLE_SECTION=${8:-skimmed_llpNanoAOD_v2}
# Optional positional arg 9: the selections file the channels are read from, relative to
# sidm/ (SidmProcessor's selections_cfg). Studies frozen to the definitions they ran with
# pass their own file, e.g. configs/selections_mulj_decay_angle_as_run.yaml.
SELECTIONS_CFG=${9:-configs/selections.yaml}

echo "Host:"
hostname
date

echo "CONDOR scratch:"
echo "${_CONDOR_SCRATCH_DIR:-undefined}"

cd "${_CONDOR_SCRATCH_DIR}"

echo "Initial files:"
ls -lah

# ------------------------------------------------------------
# 1. Setup CMSSW base environment
# ------------------------------------------------------------
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12

scramv1 project CMSSW CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
cd "${_CONDOR_SCRATCH_DIR}"

echo "Base Python:"
which python3
python3 --version

# ------------------------------------------------------------
# 2. Create venv inside Condor scratch
# ------------------------------------------------------------
python3 -m venv sidm_venv --system-site-packages
source sidm_venv/bin/activate

echo "Venv Python:"
which python
python --version

# Keep pip cache inside scratch
export PIP_CACHE_DIR="${_CONDOR_SCRATCH_DIR}/pip_cache"

python -m pip install --upgrade pip setuptools wheel
python -m pip install --no-cache-dir -r requirements.txt -c constraints.txt

# ------------------------------------------------------------
# 3. Verify dependencies before running analysis
# ------------------------------------------------------------
python - <<'PY'
import coffea
import awkward
import uproot
import fastjet
import vector
import yaml
print("coffea:", coffea.__version__)
print("All required imports worked")
PY

# ------------------------------------------------------------
# 4. Unpack SIDM code
# ------------------------------------------------------------
tar -xzf sidm_code.tar.gz

export PYTHONPATH="${_CONDOR_SCRATCH_DIR}:${PYTHONPATH:-}"

echo "After unpacking:"
ls -lah
ls -lah sidm
ls -lah condor

FILELIST_BASENAME=$(basename "${FILELIST}")
OUTFILE="${SAMPLE}_${CHUNK}.coffea"

echo "Using filelist:"
cat "${FILELIST_BASENAME}"

# ------------------------------------------------------------
# 5. Run processor using venv python
# ------------------------------------------------------------
python condor/run_sidm_chunk.py \
    --sample "${SAMPLE}" \
    --filelist "${FILELIST_BASENAME}" \
    --output "${OUTFILE}" \
    --chunksize 50000 \
    --workers 1 \
    --channels "${CHANNELS}" \
    --hist-collections "${HIST_COLLECTIONS}" \
    --backgrounds-yaml "${NTUPLE_YAML}" \
    --backgrounds-section "${NTUPLE_SECTION}" \
    --selections-cfg "${SELECTIONS_CFG}" \
    --unweighted-hist

echo "Output:"
ls -lh "${OUTFILE}"

# ------------------------------------------------------------
# 6. Copy output to EOS
# ------------------------------------------------------------
echo "Copying to EOS:"
echo "${EOS_OUTDIR}/${OUTFILE}"

xrdcp -f "${OUTFILE}" "${EOS_OUTDIR}/${OUTFILE}"

XRDEXIT=$?
if [[ ${XRDEXIT} -ne 0 ]]; then
    echo "xrdcp failed with exit code ${XRDEXIT}"
    exit ${XRDEXIT}
fi

rm -f "${OUTFILE}"

date
echo "Done"