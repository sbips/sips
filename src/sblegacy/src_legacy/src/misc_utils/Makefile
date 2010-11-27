#Compile with gcc-2.95

#CXX=/home116/egibson/openSuSE9.0/gcc-2.95.2.1-bin/bin/g++
CXXFLAGS=-O2 

all: count8 count16 mask_img

count8: 
	$(CXX) count8.cpp $(CXXFLAGS) -o $@

count16: 
	$(CXX) count16.cpp $(CXXFLAGS) -o $@
	
mask_img: 
	$(CXX) mask_img.cpp $(CXXFLAGS) -o $@	

clean:
	rm -rf count8 count16 mask_img
