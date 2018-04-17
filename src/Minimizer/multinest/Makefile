FC = gfortran
BUILD = ../../../build
OBJ_DIR = $(BUILD)/objects
TARGET_DIR = $(OBJ_DIR)/src/Minimizer/multinest
FFLAGS = -ffree-line-length-none -cpp -fPIC -c -I $(TARGET_DIR)

all: build

.PHONY: all build clean debug release

build:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OBJ_DIR)/src
	@mkdir -p $(OBJ_DIR)/src/Minimizer
	@mkdir -p $(TARGET_DIR)

	$(FC) $(FFLAGS) utils.f90 -o $(TARGET_DIR)/utils.o
	$(FC) $(FFLAGS) utils1.f90 -o $(TARGET_DIR)/utils1.o
	$(FC) $(FFLAGS) priors.f90 -o $(TARGET_DIR)/priors.o
	$(FC) $(FFLAGS) kmeans_clstr.f90 -o $(TARGET_DIR)/kmeans_clstr.o
	$(FC) $(FFLAGS) xmeans_clstr.f90 -o $(TARGET_DIR)/xmeans_clstr.o
	$(FC) $(FFLAGS) posterior.F90 -o $(TARGET_DIR)/posterior.o
	$(FC) $(FFLAGS) nested.F90 -o $(TARGET_DIR)/nested.o
	$(FC) $(FFLAGS) cwrapper.f90 -o $(TARGET_DIR)/cwrapper.o
	rm *.mod

debug: FFLAGS += -DDEBUG -g
debug: all

release: FFLAGS += -O3
release: all

clean:
	-@rm -rvf *.o *.mod
	-@rm -rvf $(OBJ_DIR)/*
