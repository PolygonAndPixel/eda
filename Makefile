CXX = g++
CXXFLAGS = -std=c++11
LDFLAGS = -L/usr/lib -L/opt/OpenBLAS/lib -lpthread -llapack -llapacke -lgfortran

BUILD = build
OBJ_DIR = $(BUILD)/objects
APP_DIR = $(BUILD)/apps
TARGET = eda
INCLUDE = -I /opt/OpenBLAS/include/ -Iinclude
SRC = \
	$(wildcard src/*.cpp) \
	$(wildcard src/likelihood/*.cpp) \
	$(wildcard src/Minimizer/*.cpp) \

OBJECTS = $(SRC:%.cpp=$(OBJ_DIR)/%.o)
FOBJECTS = $(OBJ_DIR)/src/Minimizer/polychord/*.o

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $(APP_DIR)/$(TARGET) $(OBJECTS) $(FOBJECTS) $(LDFLAGS)


.PHONY: all build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OBJ_DIR)/src
	@mkdir -p $(OBJ_DIR)/src/Minimizer
	cd src/Minimizer/polychord && $(MAKE) release

debug: CXXFLAGS += -g
debug: all

release: CXXFLAGS += -O3
release: all

clean:
	cd src/Minimizer/polychord && $(MAKE) clean
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/$(TARGET)
