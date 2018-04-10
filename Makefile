CXX=g++-4.8
CXXFLAGS= -O3 -I /opt/OpenBLAS/include/ -L/usr/lib/lapack -L/opt/OpenBLAS/lib -lpthread -std=c++11 -g
FILES= main.cpp Minimizer/*.cpp likelihood/*.cpp helper/*

all: fortran main

fortran:
	cd Minimizer/polychord && $(MAKE)

main: main.cpp
	$(CXX) $(CXXFLAGS) $(FILES)  -llapack -llapacke -o main.exe

clean:
	rm -f main.exe
