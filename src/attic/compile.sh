#!/bin/bash

# --evelibs needed for TGeoNode and TGeoManager
# Link directly to edep-sim

g++ -g -O0 -o ConvertToTMS.exe ConvertToTMSProducts.cpp \
  -I${EDEP_SIM}/include -L ${EDEP_SIM}/lib -ledepsim_io \
  $(root-config --cflags --evelibs) \
