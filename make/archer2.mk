# make/archer2.mk -- ARCHER2 (HPE Cray)
# Load the modules first (or: source hpc/archer2/modules.sh):
#   module load PrgEnv-gnu
#   module load cray-fftw
# The ftn wrapper links Cray LibSci (BLAS/LAPACK) and cray-fftw automatically.
FC         = ftn
FFLAGS     = -O2 -fopenmp
DEBUGFLAGS = -O0 -g -fopenmp -fcheck=all -fbacktrace -Wall
MODFLAG    = -J
INCLUDES   =
LIBS       =
