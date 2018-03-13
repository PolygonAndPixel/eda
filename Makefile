GCC=g++
GCCFLAGS= -O2 

all: test

test: main.cpp
	$(GCC) $(GCCFLAGS) main.cpp Minimizer/*.cpp likelihood/*.cpp -o main

clean:
	rm -f main
