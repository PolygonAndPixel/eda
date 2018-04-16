FC = gfortran
FFLAGS = -ffree-line-length-none -cpp -fPIC -c
BUILD = ../../../build
OBJ_DIR = $(BUILD)/objects

all: build

.PHONY: all build clean debug release

build:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OBJ_DIR)/src
	@mkdir -p $(OBJ_DIR)/src/Minimizer
	@mkdir -p $(OBJ_DIR)/src/Minimizer/polychord
	$(FC) $(FFLAGS) utils.f90
	$(FC) $(FFLAGS) abort.F90
	$(FC) $(FFLAGS) array_utils.f90
	$(FC) $(FFLAGS) settings.f90
	$(FC) $(FFLAGS) calculate.f90
	$(FC) $(FFLAGS) mpi_utils.F90
	$(FC) $(FFLAGS) params.f90
	$(FC) $(FFLAGS) priors.f90
	$(FC) $(FFLAGS) random_utils.F90
	$(FC) $(FFLAGS) run_time_info.f90
	$(FC) $(FFLAGS) read_write.f90
	$(FC) $(FFLAGS) ini.f90
	$(FC) $(FFLAGS) feedback.f90
	$(FC) $(FFLAGS) clustering.f90
	$(FC) $(FFLAGS) chordal_sampling.f90
	$(FC) $(FFLAGS) generate.F90
	$(FC) $(FFLAGS) nested_sampling.F90
	$(FC) $(FFLAGS)	loglikelihood_c.F90
	$(FC) $(FFLAGS) interfaces.F90
	mv *.o $(OBJ_DIR)/src/Minimizer/polychord
	mv *.mod $(OBJ_DIR)/src/Minimizer/polychord

debug: FFLAGS += -DDEBUG -g
debug: all

release: FFLAGS += -O3
release: all

clean:
	-@rm -rvf *.o *.mod
	-@rm -rvf $(OBJ_DIR)/*
