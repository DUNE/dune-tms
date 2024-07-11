CXX = g++
ROOTCONFIG := $(shell which root-config)
ROOTCFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS := $(shell $(ROOTCONFIG) --libs)
CXXFLAGS := -g -Wall -Wextra -DNDEBUG -O0

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
EXEC = $(SRCS:.cpp=)

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
