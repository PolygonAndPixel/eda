GCC=g++
GCCFLAGS= -O3 -I /home/mhierony/lib/lapack-3.8.0/LAPACKE/include

all: test

test: main.cpp
	$(GCC) $(GCCFLAGS) main.cpp Minimizer/*.cpp likelihood/*.cpp helper/* -o main 

clean:
	rm -f main
