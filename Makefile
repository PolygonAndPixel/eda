CXX=g++-4.8
CXXFLAGS= -O3 -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -llapack -std=c++11 -g
FILES= main.cpp Minimizer/*.cpp likelihood/*.cpp helper/*

all: main

main: main.cpp
	$(CXX) $(CXXFLAGS) $(FILES)  -o main 

clean:
	rm -f main
