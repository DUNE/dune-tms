#!/bin/bash
# @Liam on Slack when this breaks

echo "Make sure you source this in the top level directory of dune-tms/"

# Source Liam's spack environment that includes dune_spack repository and package binaries
source /exp/dune/data/users/losulliv/september-spack/share/spack/setup-env.sh
# Activate our environment
spack env activate tms-env

# Load up a load of schtuff with spack
#  - Load root, geant4, clhep, expatn and edep-sim
spack load cmake    && echo "Setup cmake..."
spack load root     && echo "Setup root..."
spack load geant4   && echo "Setup geant4..."
spack load clhep    && echo "Setup clhep..."
spack load expat    && echo "Setup expat..."
spack load edep-sim && echo "Setup edep-sim..."

export TMS_DIR=${PWD}
export PATH=${PATH}:${TMS_DIR}/bin

if [ "${LD_LIBRARY_PATH}x" == "x" ]; then
    LD_LIBRARY_PATH=${TMS_DIR}/lib
else
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TMS_DIR}/lib
fi

# Add edep-sim library path if available
if [ "${EDEP_ROOT}x" != "x" ]; then
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${EDEP_ROOT}/lib
fi

# Add edep-sim library path if available
which clhep-config > /dev/null 2>&1
if [ $? ]; then
    CLHEP_ROOT=`clhep-config --prefix | sed 's/\"//g'` # (:
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CLHEP_ROOT}/lib
fi
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}

echo "TMS_DIR = ${TMS_DIR}"
echo "LD_LIBRARY_PATH = ${LD_LIBRARY_PATH}"
echo "Setup TMS environment :)"
