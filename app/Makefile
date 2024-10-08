CXX=g++
AR=ar
LD_SHARED=g++

# Need evelibs for TGeoNode and TGeoManager for geometry
ROOT_LIBS = $(shell root-config --evelibs)
ROOT_INC = $(shell root-config --cflags)

ifndef CLHEP_INC # If not defined, normal build. 
CLHEP_INC = $(shell clhep-config --prefix)/include
CLHEP_LIBS = $(shell clhep-config --libs)
else # github CI passes CLHEP_INC as an environment var, use this to set lib dir
  ifndef CLHEP_LIBS # If this is defined here it was set in src/Makefile
  CLHEP_LIBS = -L$(CLHEP_INC)/../lib64
  endif
endif

# Check if geant4-config is available
GEANT4_CONFIG := $(shell command -v geant4-config 2>/dev/null)

# Set GEANT4 and EDEP variables based on the availability of geant4-config
ifneq ($(GEANT4_CONFIG),)
GEANT4_LIBS = $(shell geant4-config --libs)
GEANT4_INC = $(shell geant4-config --prefix)/include/
EDEP_LIBS = -L $(EDEPSIM_LIB) -ledepsim_io  $(CLHEP_LIBS)
EDEP_INC = -I$(EDEPSIM_INC) -I$(GEANT4_INC) -I$(CLHEP_INC)
else
EDEP_LIBS = -L $(EDEPSIM_LIB) -ledepsim_io
# Need CLHEP for unit conversion, can probably just move this into a standalone header to remove dependency?
EDEP_INC = -I$(EDEPSIM_INC) -I$(CLHEP_INC)
endif

# The TMS includes
TMS_INC = -I../src

TOML_INC = -I../toml11

CXXFLAGS = -Wall -Wextra -g -fPIC
ifdef DEBUG
  CXXFLAGS+=-Og
else
  CXXFLAGS+=-O3
endif
ifdef VERBOSE
  CXXFLAGS+=-DDEBUG
endif
CXXFLAGS += $(ROOT_INC) $(EDEP_INC) $(TMS_INC) $(TOML_INC)

# Our lovely collection of dependencies
LIB_OBJ = $(EDEP_LIBS) $(ROOT_LIBS)

# Library directory
LIB_DIR = ../lib
LIB_TMS = -lTMS_Prod

# Exectuable directory
EXE_DIR = ../bin
MKDIR_P := mkdir -p

all: directories $(EXE_DIR)/ConvertToTMSTree.exe $(EXE_DIR)/BetheBloch_Example.exe $(EXE_DIR)/DBSCAN_test.exe $(EXE_DIR)/DrawEvents.exe $(EXE_DIR)/TrackLengthTester.exe $(EXE_DIR)/ShootRay.exe $(EXE_DIR)/BField_tester.exe $(EXE_DIR)/CherryPickEvents.exe

# Sanitize target
sanitize: CXXFLAGS += -fsanitize=address
sanitize: all

directories: $(EXE_DIR)

$(EXE_DIR):
	${MKDIR_P} $(EXE_DIR)

$(EXE_DIR)/ConvertToTMSTree.exe: ConvertToTMSTree.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ -L $(LIB_DIR) $(LIB_TMS) $(LIB_OBJ) $(BUILD_FLAGS)

$(EXE_DIR)/CherryPickEvents.exe: CherryPickEvents.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ -L $(LIB_DIR) $(LIB_TMS) $(LIB_OBJ) $(BUILD_FLAGS)

$(EXE_DIR)/BetheBloch_Example.exe: BetheBloch_Example.cpp
	$(CXX) $(CXXFLAGS) $(ROOT_INC) -o $@ $^ $(ROOT_LIBS)

$(EXE_DIR)/DBSCAN_test.exe: DBSCAN_test.cpp
	$(CXX) $(CXXFLAGS) $(ROOT_INC) -o $@ $^ $(ROOT_LIBS)

$(EXE_DIR)/DrawEvents.exe: DrawEvents.cpp
	$(CXX) $(CXXFLAGS) $(ROOT_INC) -o $@ $^ -L $(LIB_DIR) $(LIB_TMS) $(LIB_OBJ) $(BUILD_FLAGS)

$(EXE_DIR)/TrackLengthTester.exe: TrackLengthTester.cpp
	$(CXX) $(CXXFLAGS) $(ROOT_INC) -o $@ $^  -L $(LIB_DIR) $(LIB_TMS) $(LIB_OBJ) $(BUILD_FLAGS)

$(EXE_DIR)/ShootRay.exe: ShootRay.cpp
	$(CXX) $(CXXFLAGS) $(ROOT_INC) -o $@ $^  -L $(LIB_DIR) $(LIB_TMS) $(LIB_OBJ) $(BUILD_FLAGS)

$(EXE_DIR)/BField_tester.exe: BField_tester.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^  -L $(LIB_DIR) $(LIB_TMS) $(LIB_OBJ) $(BUILD_FLAGS)

clean:
	rm -f $(EXE_DIR)/*.exe
