FC = gfortran
BUILD = ../../../build
OBJ_DIR = $(BUILD)/objects
TARGET_DIR = $(OBJ_DIR)/src/Minimizer/polychord
FFLAGS = -ffree-line-length-none -cpp -fPIC -c -I $(TARGET_DIR)

all: build

.PHONY: all build clean debug release

build:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OBJ_DIR)/src
	@mkdir -p $(OBJ_DIR)/src/Minimizer
	@mkdir -p $(TARGET_DIR)

	$(FC) $(FFLAGS) utils.f90 -o $(TARGET_DIR)/utils.o
	$(FC) $(FFLAGS) run_time_info_type.f90 -o $(TARGET_DIR)/run_time_info_type.o
	$(FC) $(FFLAGS) abort.F90 -o $(TARGET_DIR)/abort.o
	$(FC) $(FFLAGS) array_utils.f90 -o $(TARGET_DIR)/array_utils.o
	$(FC) $(FFLAGS) settings.f90 -o $(TARGET_DIR)/settings.o
	$(FC) $(FFLAGS) calculate.f90 -o $(TARGET_DIR)/calculate.o
	$(FC) $(FFLAGS) mpi_utils.F90 -o $(TARGET_DIR)/mpi_utils.o
	$(FC) $(FFLAGS) params.f90 -o $(TARGET_DIR)/params.o
	$(FC) $(FFLAGS) priors.f90 -o $(TARGET_DIR)/priors.o
	$(FC) $(FFLAGS) random_utils.F90 -o $(TARGET_DIR)/random_utils.o
	$(FC) $(FFLAGS) run_time_info.f90 -o $(TARGET_DIR)/run_time_info.o
	$(FC) $(FFLAGS) read_write.f90 -o $(TARGET_DIR)/read_write.o
	$(FC) $(FFLAGS) ini.f90 -o $(TARGET_DIR)/ini.o
	$(FC) $(FFLAGS) feedback.f90 -o $(TARGET_DIR)/feedback.o
	$(FC) $(FFLAGS) clustering.f90 -o $(TARGET_DIR)/clustering.o
	$(FC) $(FFLAGS) chordal_sampling.f90 -o $(TARGET_DIR)/chordal_sampling.o
	$(FC) $(FFLAGS) generate.F90 -o $(TARGET_DIR)/generate.o
	$(FC) $(FFLAGS) nested_sampling.F90 -o $(TARGET_DIR)/nested_sampling.o
	$(FC) $(FFLAGS)	loglikelihood_c.F90 -o $(TARGET_DIR)/loglikelihood_c.o
	$(FC) $(FFLAGS) interfaces.F90 -o $(TARGET_DIR)/interfaces.o
	rm *.mod

debug: FFLAGS += -DDEBUG -g
debug: all

release: FFLAGS += -O3
release: all

clean:
	-@rm -rvf *.o *.mod
	-@rm -rvf $(OBJ_DIR)/*
