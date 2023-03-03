CC=gcc
CXX=mpicxx
CFLAGS=-O3 -std=c99
CXXFLAGS=-O3
EXECS=heated-plate heated-plate-parallel-synchronous heated-plate-parallel-asynchronous

all: $(EXECS)

heated-plate: HeatedPlate_Sequential.cpp
	$(CXX) -o $@ $<  $(CXXFLAGS)

heated-plate-parallel-synchronous: HeatedPlate_Synchronous.cpp
	$(CXX) -o $@ $<  $(CXXFLAGS)

heated-plate-parallel-asynchronous: HeatedPlate_Asynchronous.cpp
	$(CXX) -o $@ $<  $(CXXFLAGS)


clean:
	rm -f $(EXECS) *.o