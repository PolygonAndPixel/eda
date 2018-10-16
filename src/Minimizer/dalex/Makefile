CXX = scorep g++
CXXFLAGS = -std=c++14 -fopenmp -w -larpack
LDFLAGS = -L/usr/lib -L/opt/OpenBLAS/lib -L/home/mhierony/ARPACK -lpthread -llapack -llapacke -fopenmp -lgfortran -larpack

BUILD = ../../../build
OBJ_DIR = $(BUILD)/objects/src/Minimizer/dalex
TARGET_DIR = $(OBJ_DIR)/dalex
APP_DIR = $(BUILD)
TARGET = dalex
INCLUDE = -I /opt/OpenBLAS/include/ -Iinclude \
	-I ../../../include/Minimizer/dalex/include

SRC = \
	$(wildcard src/analysis/*.cpp) \
	$(wildcard src/controls/*.cpp) \
	$(wildcard src/dalex/*.cpp) \
	$(wildcard src/utils/*.cpp) \

OBJECTS = $(SRC:%.cpp=$(OBJ_DIR)/%.o)

all: build $(OBJECTS)

$(OBJ_DIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

.PHONY: all build clean debug release

build:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OBJ_DIR)/src
	@mkdir -p $(OBJ_DIR)/src/analysis
	@mkdir -p $(OBJ_DIR)/src/controls
	@mkdir -p $(OBJ_DIR)/src/dalex
	@mkdir -p $(OBJ_DIR)/src/utils

debug: CXXFLAGS += -g
debug: all

release: CXXFLAGS += -O2 
release: all

clean:
	-@rm -rvf *.o
	-@rm -rvf $(OBJ_DIR)/*