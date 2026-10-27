# make/create.mk -- KCL CREATE HPC
# Load the modules first (or: source hpc/create/modules.sh):
#   module load gcc/13.2.0
#   module load openblas/0.3.24-gcc-13.2.0-multithreading-openmp
# The spack module puts OpenBLAS on LIBRARY_PATH, so no -L needed.
FC         = gfortran
FFLAGS     = -O2 -fopenmp
DEBUGFLAGS = -O0 -g -fopenmp -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall
MODFLAG    = -J
INCLUDES   =
LIBS       = -lopenblas
