CC=gcc
CXX=g++
CFLAGS=-O3 -fopenmp -std=c99
CXXFLAGS=-O3 -fopenmp
EXECS=daxpy_static

all: $(EXECS)

matmul_serial-sol: daxpy_static.cpp
	$(CXX) -o $@ $<  $(CXXFLAGS)

clean:
	rm -f $(EXECS) *.o
