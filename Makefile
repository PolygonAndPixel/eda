CXX=g++-4.8
CXXFLAGS= -O3 -I /opt/OpenBLAS/include/ -L/usr/lib/lapack -L/opt/OpenBLAS/lib -lpthread -std=c++11 -g
FILES= main.cpp Minimizer/*.cpp likelihood/*.cpp helper/*

all: main

main: main.cpp
	$(CXX) $(CXXFLAGS) $(FILES)  -llapack -llapacke -o main

clean:
	rm -f main
