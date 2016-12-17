CXX = g++
CXXFLAGS = -s -O3 -DNDEBUG -fopenmp

all: qad

qad: qad.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	$(RM) qad
