#!/bin/bash
# @Liam on Slack when this breaks

echo "Make sure you source this in the top level directory of dune-tms/"

# Setup FNAL spack
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
# Source Liam's spack environment that includes dune_spack repository and package binaries
source /exp/dune/data/users/losulliv/spack/my_spack/setup-env.sh

# Load up a load of schtuff with spack
#  - All the good stuff is built with gcc 12.2.0 toolchain, load newer cmake
#  - Load root, geant4, clhep, expatn and edep-sim built with above toolchain
spack load gcc@12.2.0%gcc@11.4.1    && echo "Setup gcc..."
spack load cmake@3.27.7%gcc@12.2.0  && echo "Setup cmake..."
spack load root@6.28.06%gcc@12.2.0  && echo "Setup root..."
spack load geant4@10.6.1%gcc@12.2.0 && echo "Setup geant4..."
spack load clhep@2.4.6.4%gcc@12.2.0 && echo "Setup clhep..."
spack load expat@2.5.0%gcc@12.2.0   && echo "Setup expat..."
spack load edep-sim@3.2%gcc@12.2.0  && echo "Setup edep-sim..."

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

echo "TMS_DIR = ${TMS_DIR}"
echo "LD_LIBRARY_PATH = ${LD_LIBRARY_PATH}"
echo "Setup TMS environment :)"
