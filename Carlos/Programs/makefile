CXX = g++
CXXFLAGS = -std=c++17 -O2
SPICE = -I/home/wind/naif/cspice/include /home/wind/naif/cspice/lib/cspice.a -lm





all: main.x


main.x: main.cpp tools.h OrbitPropagator.h 
	$(CXX) $(CXXFLAGS) $< -o $@ $(SPICE)
	./main.x 


clean:
	rm -rf *.x a.out *.x.*
