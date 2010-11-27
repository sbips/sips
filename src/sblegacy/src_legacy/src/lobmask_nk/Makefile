#Compile with gcc-2.95

#CXX=/home116/egibson/openSuSE9.0/gcc-2.95.2.1-bin/bin/g++
CXXFLAGS=-O2 

all: lobmask lobmask_thalamus

lobmask:
	$(CXX) lobmask.cpp $(CXXFLAGS) -o $@

lobmask_thalamus:
	$(CXX) lobmask_thalamus.cpp $(CXXFLAGS) -o $@

clean:
	rm -rf lobmask lobmask_thalamus

