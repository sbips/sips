#Compile with gcc-2.95

#CXX=/home116/egibson/openSuSE9.0/gcc-3.3-bin/bin/g++
CXXFLAGS=-O2 
INCLUDES=../../../include_from_cortex/MyVolume.cpp ../../../include_from_cortex/connexe.c -I/cygdrive/c/erindocs/work/legacy/include_from_cortex/ ../../../include_from_cortex/MyImage.cpp

all: Lsegment fillHoles cleanEdgeSeg 

Lsegment: 
	$(CXX) Lsegment.cpp $(CXXFLAGS)  $(INCLUDES) -o $@

fillHoles:
	$(CXX) -o fillHoles fillHoles.cpp $(CXXFLAGS) 

cleanEdgeSeg: 
	$(CXX) cleanEdgeSeg.cpp $(CXXFLAGS) $(INCLUDES) -o $@


clean:
	rm -rf Lsegment fillHoles cleanEdgeSeg

