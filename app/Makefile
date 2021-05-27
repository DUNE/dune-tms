CXX=g++
AR=ar
LD_SHARED=g++

# Need evelibs for TGeoNode and TGeoManager
ROOT_LIBS = $(shell root-config --evelibs)
ROOT_INC = $(shell root-config --cflags)

EDEP_LIBS = -L $(EDEP_SIM)/lib -ledepsim_io
# Need CLHEP for unit conversion, can probably just move this into a standalone header to remove dependency?
EDEP_INC = -I$(EDEP_SIM)/include -I$(CLHEP_INC)

# The TMS includes
TMS_INC = -I../src

CXXFLAGS = -Wall -Wextra -g -fPIC -Werror
ifdef DEBUG
  CXXFLAGS+=-O0
else
  CXXFLAGS+=-O3
endif
CXXFLAGS += $(ROOT_INC) $(EDEP_INC) $(TMS_INC)

# Our lovely collection of dependencies
LIB_OBJ = $(EDEP_LIBS) $(ROOT_LIBS)

# Library directory
LIB_DIR = ../lib
LIB_TMS = -lTMS_Prod

# Exectuable directory
EXE_DIR = ../bin
MKDIR_P := mkdir -p

all: directories $(EXE_DIR)/ConvertToTMSTree.exe $(EXE_DIR)/BetheBloch_Example.exe $(EXE_DIR)/DBSCAN_test.exe

directories: $(EXE_DIR)

$(EXE_DIR):
	${MKDIR_P} $(EXE_DIR)

$(EXE_DIR)/ConvertToTMSTree.exe: ConvertToTMSTree.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ -L $(LIB_DIR) $(LIB_TMS) $(LIB_OBJ) $(BUILD_FLAGS)

$(EXE_DIR)/BetheBloch_Example.exe: BetheBloch_Example.cpp
	$(CXX) $(CXXFLAGS) $(ROOT_INC) -o $@ $^ $(ROOT_LIBS)

$(EXE_DIR)/DBSCAN_test.exe: DBSCAN_test.cpp
	$(CXX) $(CXXFLAGS) $(ROOT_INC) -o $@ $^ $(ROOT_LIBS)

clean:
	rm -f $(EXE_DIR)/*.exe
