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

#setup root v6_22_08b -f Linux64bit+3.10-2.17 -q e20:p383b:prof
setup cmake v3_17_3 -f Linux64bit+3.10-2.17

# Official duneanaobj UPS products doesn't include TMS info yet...
#setup duneanaobj v01_00_00 -f Linux64bit+3.10-2.17 -q debug:e17:gv2

# Setup custom UPS products, for CAFana format development
#source /dune/app/users/cwret/ups/setup.sh
#setup duneanaobj v02_01_01 -q e20:gv3:prof
# Flag for enabling linking against DUNEANAOBJ
#if [ ! -z "${DUNEANAOBJ_INC}" ]; then
  #export DUNEANAOBJ_ENABLED=1
  #echo "Enabling DUNEANAOBJ..."
#fi

# Need CLHEP for unit conversion in edep-sim info
setup clhep v2_4_4_1 -f Linux64bit+3.10-2.17 -q debug:e20
# Specific requirements for duneanaobj
setup root v6_12_06a -f Linux64bit+3.10-2.17 -q e15:prof

# Set up edep-sim
setup edepsim v3_2_0 -f Linux64bit+3.10-2.17 -q e20:prof

# Add TMS execs and library directory to env
export TMS_DIR=${PWD}
export PATH=${PATH}:${TMS_DIR}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TMS_DIR}/lib

# Debugging options passed to makefile
export DEBUG=1
#export VERBOSE=1


echo "Setup TMS environment"
