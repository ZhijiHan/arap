CXX=g++
BIN_ROOT=bin
RESEARCH_ROOT=/home/taodu/research

# IGL include.
IGL_INC=-I$(RESEARCH_ROOT)/libigl/include
INC+=$(IGL_INC)

# Eigen include.
EIGEN_INC=-I$(RESEARCH_ROOT)/eigen -I$(RESEARCH_ROOT)/eigen/unsupported
INC+=$(EIGEN_INC)

CXXFLAGS=-Wall -c $(INC)

all: arap

arap: main.o
	$(CXX) main.o -o $(BIN_ROOT)/arap

main.o: main.cc
	$(CXX) $(CXXFLAGS) main.cc

clean:
	rm -rf *.o arap
