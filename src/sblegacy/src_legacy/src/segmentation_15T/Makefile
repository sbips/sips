#Compile with gcc-2.95

#CXX=/home116/egibson/openSuSE9.0/gcc-2.95.2.1-bin/bin/g++
CXXFLAGS=-O2 
SOURCES=T1seg.cpp
EXECUTABLE=T1seg

$(EXECUTABLE): 
	$(CXX) $(SOURCES) $(CXXFLAGS) -o $@

clean:
	rm -rf $(EXECUTABLE)

