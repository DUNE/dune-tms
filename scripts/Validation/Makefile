CXX = g++
ROOTCONFIG := $(shell which root-config)
ROOTCFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS := $(shell $(ROOTCONFIG) --libs)
CXXFLAGS := -g -Wall -Wextra -DNDEBUG -O0

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
EXEC = $(SRCS:.cpp=)

BRANCH_NAME := $(shell git rev-parse --abbrev-ref HEAD)
GIT_STATUS := $(shell git diff-index --quiet HEAD -- src/ app/ || echo "-modified")

CXXFLAGS += -DGIT_BRANCH_NAME=\"$(BRANCH_NAME)$(GIT_STATUS)\"

.PHONY: all clean

all: $(EXEC)

#$(EXEC): $(OBJS)
#$(EXEC): %: %.o 
$(EXEC): %: %.o
	$(CXX) $*.o -o $@ $(ROOTLIBS)

 %.o: %.cpp
	$(CXX) $(ROOTCFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
