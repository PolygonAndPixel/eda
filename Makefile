GCC=g++
GCCFLAGS= -O3 -I /home/mhierony/lib/lapack-3.8.0/LAPACKE/include -std=c++11
FILES= Minimizer/SampleSpace.cpp Minimizer/SampleSpace.h likelihood/*.cpp helper/*

all: test

test: main.cpp
	$(GCC) $(GCCFLAGS) main.cpp $(FILES)  -o main 

clean:
	rm -f main
