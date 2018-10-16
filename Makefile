CXX = scorep g++
CXXFLAGS = -std=c++14 -fopenmp  -w
CXX_THIRD_PARTY =
LDFLAGS = -L/usr/lib -L/opt/OpenBLAS/lib -L/home/mhierony/ARPACK -lpthread -llapack -llapacke -fopenmp -lgfortran -larpack
FFLAG =

BUILD = build
OBJ_DIR = $(BUILD)/objects
APP_DIR = $(BUILD)/apps
TARGET = eda
INCLUDE = -I /opt/OpenBLAS/include/ -Iinclude
SRC = \
	$(wildcard src/*.cpp) \
	$(wildcard src/IceCubeToy/*.cpp) \
	$(wildcard src/likelihood/*.cpp) \
	$(wildcard src/Minimizer/*.cpp) \

OBJECTS = $(SRC:%.cpp=$(OBJ_DIR)/%.o)
FOBJECTS = $(OBJ_DIR)/src/Minimizer/polychord/*.o $(OBJ_DIR)/src/Minimizer/multinest/*.o
THIRD_OBJECTS = $(OBJ_DIR)/src/Minimizer/dalex/src/analysis/*.o \
	$(OBJ_DIR)/src/Minimizer/dalex/src/controls/*.o \
	$(OBJ_DIR)/src/Minimizer/dalex/src/dalex/*.o \
	$(OBJ_DIR)/src/Minimizer/dalex/src/utils/*.o \
	# $(OBJ_DIR)/src/Minimizer/dalex/src/examples/*.o \
	# $(OBJ_DIR)/src/Minimizer/dalex/src/tests/*.o \
	

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $(APP_DIR)/$(TARGET) $(OBJECTS) $(FOBJECTS) $(THIRD_OBJECTS) $(LDFLAGS)


.PHONY: all build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OBJ_DIR)/src
	@mkdir -p $(OBJ_DIR)/src/Minimizer
	cd src/Minimizer/polychord && $(MAKE) $(FFLAG)
	cd src/Minimizer/multinest && $(MAKE) $(FFLAG)
	cd src/Minimizer/dalex && $(MAKE) $(CXX_THIRD_PARTY)

debug: CXXFLAGS += -g
debug: FFLAG += debug
debug: CXX_THIRD_PARTY += debug
debug: all

release: CXXFLAGS += -O2 -Wno-conversion-null 
release: FFLAG += release
release: CXX_THIRD_PARTY += release
release: all

clean:
	cd src/Minimizer/polychord && $(MAKE) clean
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/$(TARGET)
	-@rm -rvf $(APP_DIR)/tmp/*
