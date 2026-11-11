# make/create.mk -- KCL CREATE HPC
# Load the modules first:  source hpc/create/modules.sh
#   (gcc/12.5.0, openblas/0.3.33, fftw/3.3.11 as of Nov 2026; check: module spider fftw)
# The spack modules put the include paths on CPATH or C_INCLUDE_PATH (fftw uses the
# latter) and the libraries on LIBRARY_PATH, so no -I/-L is needed here. They do NOT
# set LD_LIBRARY_PATH: modules.sh adds the fftw and openblas lib folders for run time.
# The spack libraries are built per CPU type (login node: Intel cascadelake; cpu
# partition: AMD zen2/zen3/zen5). A binary built on the login node dies on the compute
# nodes with "Illegal instruction": hpc/create/job.slurm builds on the node it runs on.
# If FFTW is not found, pass its prefix: make PLATFORM=create FFTW_DIR=/path/to/fftw
FC         = gfortran
FFLAGS     = -O2 -fopenmp
DEBUGFLAGS = -O0 -g -fopenmp -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall
MODFLAG    = -J
ifdef FFTW_DIR
  INCLUDES = -I$(FFTW_DIR)/include
  LIBS     = -L$(FFTW_DIR)/lib -Wl,-rpath,$(FFTW_DIR)/lib -lfftw3 -lopenblas
else
  INCLUDES = $(addprefix -I,$(subst :, ,$(CPATH) $(C_INCLUDE_PATH)))
  LIBS     = -lfftw3 -lopenblas
endif
