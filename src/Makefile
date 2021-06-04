CXX=g++
AR=ar
LD_SHARED=g++

# Need evelibs for TGeoNode and TGeoManager
ROOT_LIBS = $(shell root-config --evelibs)
ROOT_INC = $(shell root-config --cflags)

EDEP_LIBS = -L $(EDEP_SIM)/lib -ledepsim_io
# Need CLHEP for unit conversion, can probably just move this into a standalone header to remove dependency?
EDEP_INC = -I$(EDEP_SIM)/include -I$(CLHEP_INC)

CXXFLAGS = -Wall -Wextra -g -fPIC -Werror
ifdef DEBUG
  CXXFLAGS+=-O0
else
  CXXFLAGS+=-O3
endif
ifdef VERBOSE
  CXXFLAGS+=-DDEBUG
endif
CXXFLAGS += $(ROOT_INC) $(EDEP_INC)

# Replace and create
ARFLAGS = rc

# Get the ROOT CXX flags in the shared lib
SO_FLAGS = -shared $(ROOT_INC)

# Include this directory when looking for libs
BUILD_FLAGS = -Wl,-R -Wl,.

MKDIR_P := mkdir -p

# Library directory where we build to
LIB_DIR=../lib

# Our lovely collection of libs
LIB_OBJ = $(EDEP_LIBS) $(ROOT_LIBS)

TMS_OBJ = TMS_Bar.o TMS_Hit.o TMS_TrueHit.o TMS_Event.o TMS_TrueParticle.o TMS_EventViewer.o TMS_Reco.o TMS_Kalman.o TMS_TreeWriter.o

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $<

all: directories libTMS_Prod

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
