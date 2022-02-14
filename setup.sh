#!/bin/bash

# Example of setup needed for SSRI dependencies
# Here we use our own version of edep-sim, but you should be able to use one set up by ups

MY_SSRI_DIR=/dune/app/users/cwret/SSRI

# Copied from Gavin's SSRI setup
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

#setup dk2nu        v01_05_01b   -q e15:prof
#setup genie        v2_12_10c    -q e15:prof
#setup genie_xsec   v2_12_10     -q DefaultPlusValenciaMEC
#setup genie_phyopt v2_12_10     -q dkcharmtau
#setup geant4       v4_10_3_p01b -q e15:prof
#setup ifdhc

# Specific requirements for duneanaobj
setup root v6_12_06a -f Linux64bit+3.10-2.17 -q e15:prof
#setup root v6_22_08b -f Linux64bit+3.10-2.17 -q e20:p383b:prof
setup cmake v3_17_3 -f Linux64bit+3.10-2.17
#setup duneanaobj v01_00_00 -f Linux64bit+3.10-2.17 -q debug:e17:gv2

# Setup 
source /dune/app/users/cwret/ups/setup.sh
# Flag for enabling linking against DUNEANAOBJ
if [ ! -z "${DUNEANAOBJ_INC}" ]; then
  export DUNEANAOBJ_ENABLED=1
  echo "Enabling DUNEANAOBJ..."
fi


# edep-sim needs to find geant4.sh to build
#export PATH=$PATH:${GEANT4_FQ_DIR}/bin

# edep-sim also needs cmake for setup
#export GXMLPATH=$MY_SSRI_DIR/gevgen:${GXMLPATH}
#export GNUMIXML="GNuMIFlux.xml"

# Add custom edep-sim (can also use UPS product)
export EDEP_SIM=${MY_SSRI_DIR}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu
export LD_LIBRARY_PATH=${EDEP_SIM}/lib:${LD_LIBRARY_PATH}
export PATH=${EDEP_SIM}/bin:${PATH}

# Add TMS execs and library directory to env
export TMS_DIR=${PWD}
export PATH=${PATH}:${TMS_DIR}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TMS_DIR}/lib


echo "Setup TMS environment"
