#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd ${SCRIPT_DIR}
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval $(scramv1 runtime -sh)

export G4REALSURFACEDATA=$CMSSW_BASE/src/HGCalTileSim/Tile/data/SurfaceData
${CMSSW_BASE}/bin/${SCRAM_ARCH}/LYSim_SquareTrigger $@
