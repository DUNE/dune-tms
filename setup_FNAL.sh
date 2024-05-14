#!/bin/bash
# @Liam on Slack when this breaks

echo "Make sure you source this in the top level directory of dune-tms/"

# Setup FNAL spack
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
# Source Liam's spack environment that includes dune_spack repository and package binaries
source /exp/dune/data/users/losulliv/spack/my_spack/setup-env.sh

# Load up a load of schtuff with spack, check if packages are loaded first
#  - All the good stuff is built with gcc 12.2.0 toolchain, load newer cmake
#  - Load root, geant4, clhep, expatn and edep-sim built with above toolchain
#  - We first check the list of loaded packages using grep to see if we have them loaded
#    (e.g. already sourced this file) and load it with spack if not

spack load --list | grep "gcc@12.2.0" >/dev/null 2>&1 # # 0(true) if it's loaded, 1(false) otherwise
if [ $? -ne 0 ]; then spack load gcc@12.2.0%gcc@11.4.1 && echo "Setup gcc..."; else echo "gcc already loaded."; fi
spack load --list | grep "cmake@3.27.7" >/dev/null 2>&1 # # 0(true) if it's loaded, 1(false) otherwise
if [ $? -ne 0 ]; then spack load cmake@3.27.7%gcc@12.2.0 && echo "Setup cmake..."; else echo "cmake already loaded."; fi
spack load --list | grep "root@6.28.06" >/dev/null 2>&1 # # 0(true) if it's loaded, 1(false) otherwise
if [ $? -ne 0 ]; then spack load root@6.28.06%gcc@12.2.0 && echo "Setup root..."; else echo "root already loaded."; fi
spack load --list | grep "geant4@10.6.1" >/dev/null 2>&1 # # 0(true) if it's loaded, 1(false) otherwise
if [ $? -ne 0 ]; then spack load geant4@10.6.1%gcc@12.2.0 && echo "Setup geant4..."; else echo "geant4 already loaded."; fi
spack load --list | grep "clhep@2.4.6.4" >/dev/null 2>&1 # # 0(true) if it's loaded, 1(false) otherwise
if [ $? -ne 0 ]; then spack load clhep@2.4.6.4%gcc@12.2.0 && echo "Setup clhep..."; else echo "clhep already loaded."; fi
spack load --list | grep "expat@2.5.0" >/dev/null 2>&1 # # 0(true) if it's loaded, 1(false) otherwise
if [ $? -ne 0 ]; then spack load expat@2.5.0%gcc@12.2.0 && echo "Setup expat..."; else echo "expat already loaded."; fi
spack load --list | grep "edep-sim@3.2" >/dev/null 2>&1 # # 0(true) if it's loaded, 1(false) otherwise
if [ $? -ne 0 ]; then spack load edep-sim@3.2%gcc@12.2.0 && echo "Setup edep-sim..."; else echo "edep-sim already loaded."; fi

export TMS_DIR=${PWD}
# check if TMS is already in our path, assume PATH is never empty
echo ${PATH} | grep ${TMS_DIR} > /dev/null 2>&1
if [ $? ]; then
    export PATH=${PATH}:${TMS_DIR}/bin
fi

# Check if LD_LIBRARY_PATH is set, create if not, append if so
if [ "${LD_LIBRARY_PATH}x" == "x" ]; then
    LD_LIBRARY_PATH=${TMS_DIR}/lib
else
    echo ${PATH} | grep ${TMS_DIR} > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TMS_DIR}/lib
    fi
fi

# Add edep-sim library path if available
if [ "${EDEP_ROOT}x" != "x" ]; then
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${EDEP_ROOT}/lib
fi

# Do some housekeeping here to keep users happy, nominally these _should_ already be set but *insert permissible spack joke*

# Add edep-sim library path if available
which clhep-config > /dev/null 2>&1 # check if clhep-config executable is available
if [  $? ]; then
    CLHEP_ROOT=`clhep-config --prefix | sed 's/\"//g'` # (:
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CLHEP_ROOT}/lib
fi

# Check if we have Geant4 libraries in the LD_LIBRARY_PATH
echo ${LD_LIBRARY_PATH} | grep geant >/dev/null 2>&1 # 0(true) if 'geant' found in path, 1(false) otherwise
if [ $? -ne 0 ]; then
    echo "Geant4 library path not found in LD_LIBRARY_PATH, adding it..."
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:`geant4-config --prefix`/lib64
fi

# Check if we have root libraries in the LD_LIBRARY_PATH
echo ${LD_LIBRARY_PATH} | grep "root-" >/dev/null 2>&1 # 0(true) if 'geant' found in path, 1(false) otherwise
if [ $? -ne 0 ]; then
    echo "ROOT library path not found in LD_LIBRARY_PATH, adding it..."
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:`root-config --prefix`/lib/root
fi

echo "TMS_DIR = ${TMS_DIR}"
#echo "PATH = ${PATH}" # Commented as this long
#echo "LD_LIBRARY_PATH = ${LD_LIBRARY_PATH}" # Commented as this long
echo ""
echo "Setup TMS environment :)"
