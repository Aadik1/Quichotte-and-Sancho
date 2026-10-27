# make/create.mk -- KCL CREATE HPC
# Load the modules first (or: source hpc/create/modules.sh):
#   module load gcc/13.2.0
#   module load openblas/0.3.24-gcc-13.2.0-multithreading-openmp
#   module load fftw          <- check the exact name with: module avail fftw
# The spack modules put the libraries on LIBRARY_PATH/CPATH, so no -L/-I needed.
# If FFTW is not found, pass its prefix: make PLATFORM=create FFTW_DIR=/path/to/fftw
FC         = gfortran
FFLAGS     = -O2 -fopenmp
DEBUGFLAGS = -O0 -g -fopenmp -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall
MODFLAG    = -J
ifdef FFTW_DIR
  INCLUDES = -I$(FFTW_DIR)/include
  LIBS     = -L$(FFTW_DIR)/lib -Wl,-rpath,$(FFTW_DIR)/lib -lfftw3 -lopenblas
else
  INCLUDES = $(addprefix -I,$(subst :, ,$(CPATH)))
  LIBS     = -lfftw3 -lopenblas
endif
