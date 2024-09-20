CXX=g++
AR=ar
LD_SHARED=g++

# Need evelibs for TGeoNode and TGeoManager
ROOT_LIBS = $(shell root-config --evelibs)
ROOT_INC = $(shell root-config --cflags)

ifndef CLHEP_INC # If not defined, normal build. 
CLHEP_INC = $(shell clhep-config --prefix)/include
CLHEP_LIBS = $(shell clhep-config --libs)
else # github CI passes CLHEP_INC as an environment var, use this to set lib dir
CLHEP_LIBS = -L$(CLHEP_INC)/../../lib64
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

TOML_INC = -I../toml11/

CXXFLAGS = -Wall -Wextra -g -fPIC
ifdef DEBUG
  CXXFLAGS+=-Og
else
  CXXFLAGS+=-O3
endif
ifdef VERBOSE
  CXXFLAGS+=-DDEBUG
endif
CXXFLAGS += $(ROOT_INC) $(EDEP_INC) $(TOML_INC)

# Replace and create
ARFLAGS = rc

# Get the ROOT CXX flags in the shared lib
SO_FLAGS = -shared $(ROOT_INC)

# Include this directory when looking for libs
BUILD_FLAGS = -Wl,-R -Wl,.

MKDIR_P := mkdir -p

# Library directory where we build to
LIB_DIR=../lib

TMS_OBJ = TMS_Bar.o TMS_Hit.o TMS_TrueHit.o TMS_Event.o TMS_Track.o TMS_TrueParticle.o TMS_EventViewer.o TMS_Reco.o TMS_ChargeID.o TMS_Kalman.o TMS_TreeWriter.o TMS_ReadoutTreeWriter.o TMS_Manager.o TMS_Readout_Manager.o BField_Handler.o TMS_Utils.o

# Our lovely collection of libs
LIB_OBJ = $(EDEP_LIBS) $(ROOT_LIBS)

# If linking against DUNE anaobj
ifdef DUNEANAOBJ_ENABLED
  CXXFLAGS += -I$(DUNEANAOBJ_INC)
  CXXFLAGS += -DDUNEANAOBJ_ENABLED
  #TMS_OBJ += TMS_CAF_converter.o
  LIB_OBJ += -L $(DUNEANAOBJ_LIB) -lduneanaobj_StandardRecord
endif

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $<

all: directories libTMS_Prod

# Sanitize target
sanitize: CXXFLAGS := $(CXXFLAGS) -fsanitize=address
sanitize: all

directories: $(LIB_DIR)

$(LIB_DIR):
	${MKDIR_P} $(LIB_DIR)

libTMS_Prod: $(LIB_DIR)/libTMS_Prod.so

# Make the archive or static lib
libTMS_Prod.a: $(TMS_OBJ)
	$(AR) $(ARFLAGS) $@ $^

# Make the shared object
$(LIB_DIR)/libTMS_Prod.so: libTMS_Prod.a
	$(LD_SHARED) $(SO_FLAGS) $(CXXFLAGS) -o $@ -Wl,--whole-archive $^ -Wl,--no-whole-archive $(LIB_OBJ)

clean:
	rm -f $(LIB_DIR)/lib* *.o *.a *.exe
