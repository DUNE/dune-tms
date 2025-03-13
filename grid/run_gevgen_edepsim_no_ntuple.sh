#! /usr/bin/env bash

##################################################

HORN=$1
FIRST=$2
NPER=$3
GEOMETRY=$4
TEST=$5
if [ "${HORN}" != "FHC" ] && [ "${HORN}" != "RHC" ]; then
echo "Invalid beam mode ${HORN}"
echo "Must be FHC or RHC"
kill -INT $$
fi

MODE="neutrino"
RHC=""
if [ "${HORN}" = "RHC" ]; then
MODE="antineutrino"
RHC=" --rhc"
fi

if [ "${NPER}" = "" ]; then
echo "Number of events per job not specified, using 1000"
NPER=1000
fi

if [ "${FIRST}" = "" ]; then
echo "First run number not specified, using 0"
FIRST=0
fi

CP="ifdh cp"
if [ "${TEST}" = "test" ]; then
echo "Test mode"
PROCESS=0
mkdir -p test
cd test
fi

echo "Running edepsim for ${HORN} mode, ${NPER} events"

RNDSEED=$((${PROCESS}+${FIRST}))
# 5E16 is about 15000 events which runs in ~6 hours
NEVENTS="-n ${NPER}"      # -n XXXX number of events, -e XE16 for POT

TOPVOL="volArgonCubeActive"

TOPDIR="/pnfs/dune/scratch/users/gsdavies/ssri"

##################################################

## Setup UPS and required products

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup dk2nu        v01_05_01b   -q e15:prof
setup genie        v2_12_10c    -q e15:prof
setup genie_xsec   v2_12_10     -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10     -q dkcharmtau
setup geant4       v4_10_3_p01b -q e15:prof
setup ifdhc

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

## Copy GENIE stuff to the local node
echo "Fetching GENIE stuff"
${CP} ${TOPDIR}/gevgen.tar.gz gevgen.tar.gz
tar -xzf gevgen.tar.gz

# tarball has contents in a folder to avoid tarbombing while testing
mv gevgen/* .

# Get flux files to local node
echo "Retrieving flux files"
chmod +x copy_dune_ndtf_flux
./copy_dune_ndtf_flux --top /pnfs/dune/persistent/users/dbrailsf/flux/nd/gsimple/v2_8_6d --output local_flux_files --flavor ${MODE} --base OptimizedEngineeredNov2017 --maxmb=60

####################

## Add the location of the GNuMIFlux.xml to the GENIE xml path

export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

## Run GENIE
echo "Running gevgen"
gevgen_fnal \
    -f local_flux_files/gsimple*.root,DUNEND \
    -g ${GEOMETRY}.gdml \
    -t ${TOPVOL} \
    -m ${GEOMETRY}.${TOPVOL}.maxpl.xml \
    -L cm -D g_cm3 \
    ${NEVENTS} \
    --seed ${RNDSEED} \
    -r ${RNDSEED} \
    -o ${MODE} \
    --message-thresholds Messenger_production.xml \
    --cross-sections ${GENIEXSECPATH}/gxspl-FNALsmall.xml \
    --event-record-print-level 0 \
    --event-generator-list Default+CCMEC

##################################################

## Convert the genie output to rootracker

echo "Fetching edep-sim stuff"
${CP} ${MODE}.${RNDSEED}.ghep.root input_file.ghep.root

## Copy a tarball of all the code we need
${CP} ${TOPDIR}/edep-sim.tar.gz ${PWD}/edep-sim.tar.gz
tar -xzf edep-sim.tar.gz

export LD_LIBRARY_PATH=${PWD}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/lib:${LD_LIBRARY_PATH}
export PATH=${PWD}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/bin:${PATH}
export EDEPSIM_ROOT=$(dirname $(which edep-sim))/..

echo "Converting ghep to rootracker"
gntpc -i input_file.ghep.root -f rootracker \
      --event-record-print-level 0 \
      --message-thresholds Messenger_production.xml

# run all the events in a file
NPER=$(echo "std::cout << gtree->GetEntries() << std::endl;" | genie -l -b input_file.ghep.root 2>/dev/null  | tail -1)
echo "Running edep-sim with all ${NPER} events"

##################################################

## Run edep-sim
echo "Running edep-sim"
edep-sim \
    -C \
    -g ${GEOMETRY}.gdml \
    -o ${PWD}/edep.${RNDSEED}.root \
    -u \
    -e ${NPER} \
    dune-nd.mac

## Copy the output files

#${CP} ${MODE}.${RNDSEED}.ghep.root ${TOPDIR}/outputs/ghep/${GEOMETRY}/${MODE}.${RNDSEED}.ghep.root
${CP} edep.${RNDSEED}.root ${TOPDIR}/outputs/edep/orig/${MODE}.${RNDSEED}.edepsim.root
#${CP} dump.root ${TOPDIR}/outputs/ntuple/${GEOMETRY}/${MODE}.${RNDSEED}.dump.root
