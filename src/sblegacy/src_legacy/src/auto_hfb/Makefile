#Compile with gcc-2.95

#CXX=/home116/egibson/openSuSE9.0/gcc-2.95.2.1-bin/bin/g++
CXXFLAGS=-O2 
SOURCES=auto_hfb_revised.cpp
EXECUTABLE=auto_hfb

$(EXECUTABLE): 
	$(CXX) $(SOURCES) $(CXXFLAGS) -o $@

clean:
	rm -rf $(EXECUTABLE)

