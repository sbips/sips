#Compile with gcc-2.95

#CXX=/home116/egibson/openSuSE9.0/gcc-3.3-bin/bin/g++
CXXFLAGS=-O2 
SOURCES=3T_T1seg_v4.cpp
EXECUTABLE=T1seg_3T

$(EXECUTABLE): 
	$(CXX) $(SOURCES) $(CXXFLAGS) -o $@

clean:
	rm -rf $(EXECUTABLE)

