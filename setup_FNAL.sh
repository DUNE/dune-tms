#!/bin/bash
# @Liam on Slack when this breaks

echo "Make sure you source this in the top level directory of dune-tms/"

# Setup FNAL spack
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
# Source Liam's spack environment that includes dune_spack repository and package binaries
source /exp/dune/data/users/losulliv/spack/my_spack/setup-env.sh

# Load up a load of schtuff with spack, check if packages are loaded first
#  - All the good stuff is built with gcc 12.2.0 toolchain, load newer cmake
#  - Load root, geant4, clhep, expat, and edep-sim built with above toolchain
#  - We first check the list of loaded packages using grep to see if we have them loaded
#    (e.g. already sourced this file) and load it with spack if not

LOADED_SPACKAGES=`spack load --list` # Save to a variable so we only call this once, it's not fast
echo ${LOADED_SPACKAGES} | grep "gcc@12.2.0" >/dev/null 2>&1 # # 0(true) if it's loaded, 1(false) otherwise
if [ $? -ne 0 ]; then spack load gcc@12.2.0%gcc@11.4.1 && echo "Setup gcc..."; else echo "gcc already loaded."; fi
echo ${LOADED_SPACKAGES} | grep "cmake@3.27.7" >/dev/null 2>&1
if [ $? -ne 0 ]; then spack load cmake@3.27.7%gcc@12.2.0 && echo "Setup cmake..."; else echo "cmake already loaded."; fi
echo ${LOADED_SPACKAGES} | grep "root@6.28.06" >/dev/null 2>&1
if [ $? -ne 0 ]; then spack load root@6.28.06%gcc@12.2.0 && echo "Setup root..."; else echo "root already loaded."; fi
echo ${LOADED_SPACKAGES} | grep "geant4@10.6.1" >/dev/null 2>&1
if [ $? -ne 0 ]; then spack load geant4@10.6.1%gcc@12.2.0 && echo "Setup geant4..."; else echo "geant4 already loaded."; fi
echo ${LOADED_SPACKAGES} | grep "clhep@2.4.6.4" >/dev/null 2>&1
if [ $? -ne 0 ]; then spack load clhep@2.4.6.4%gcc@12.2.0 && echo "Setup clhep..."; else echo "clhep already loaded."; fi
echo ${LOADED_SPACKAGES} | grep "expat@2.5.0" >/dev/null 2>&1
if [ $? -ne 0 ]; then spack load expat@2.5.0%gcc@12.2.0 && echo "Setup expat..."; else echo "expat already loaded."; fi
echo ${LOADED_SPACKAGES} | grep "edep-sim@3.2" >/dev/null 2>&1
if [ $? -ne 0 ]; then spack load edep-sim@3.2%gcc@12.2.0 && echo "Setup edep-sim..."; else echo "edep-sim already loaded."; fi
unset LOADED_SPACKAGES

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

# check if expat in LD_LIBRARY_PATH
echo ${LD_LIBRARY_PATH} | grep expat > /dev/null 2>&1
if [ $? -ne 0 ]; then
    # Add expat library path 
    EXPAT_PATH=`spack find --paths --loaded expat | sed "2q;d" | cut -c14-` # Just trust in spack
    echo ${EXPAT_PATH} | grep my_spack > /dev/null 2>&1
    if [ $? -ne 0 ]; then
    echo "Adding expat to LD_LIBRARY_PATH..."
        LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${EXPAT_PATH}/lib
    fi
fi

# check if edep-sim in LD_LIBRARY_PATH
echo ${LD_LIBRARY_PATH} | grep edep-sim > /dev/null 2>&1
if [ $? -ne 0 ]; then
    # Add edep-sim library path if available
    if [ "${EDEP_ROOT}x" != "x" ]; then
        echo "Adding edep-sim to LD_LIBRARY_PATH..."
        LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${EDEP_ROOT}/lib
    fi
fi

echo ${LD_LIBRARY_PATH} | grep edep-sim > /dev/null 2>&1
if [ $? -ne 0 ]; then
    # Add edep-sim library path if available
    which clhep-config > /dev/null 2>&1 # check if clhep-config executable is available
    if [  $? -ne 0 ]; then
        CLHEP_ROOT=`clhep-config --prefix | sed 's/\"//g'` # (:
        echo "Adding clhep to LD_LIBRARY_PATH..."
        LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CLHEP_ROOT}/lib
    fi
fi

# Check if we have Geant4 libraries in the LD_LIBRARY_PATH
echo ${LD_LIBRARY_PATH} | grep geant >/dev/null 2>&1 # 0(true) if 'geant' found in path, 1(false) otherwise
if [ $? -ne 0 ]; then
    echo "Adding geant4 to LD_LIBRARY_PATH..."
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:`geant4-config --prefix`/lib64
fi

# Check if we have root libraries in the LD_LIBRARY_PATH
echo ${LD_LIBRARY_PATH} | grep "root-" >/dev/null 2>&1 # 0(true) if 'geant' found in path, 1(false) otherwise
if [ $? -ne 0 ]; then
    echo "Adding root to LD_LIBRARY_PATH..."
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:`root-config --prefix`/lib/root
fi

echo "TMS_DIR = ${TMS_DIR}"
#echo "PATH = ${PATH}" # Commented as this long
#echo "LD_LIBRARY_PATH = ${LD_LIBRARY_PATH}" # Commented as this long
echo ""
echo "Setup TMS environment :)"
