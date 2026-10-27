# make/archer2.mk -- ARCHER2 (HPE Cray)
# Load the modules first (or: source hpc/archer2/modules.sh):
#   module load PrgEnv-gnu
# The ftn wrapper links Cray LibSci (BLAS/LAPACK) automatically.
FC         = ftn
FFLAGS     = -O2 -fopenmp
DEBUGFLAGS = -O0 -g -fopenmp -fcheck=all -fbacktrace -Wall
MODFLAG    = -J
INCLUDES   =
LIBS       =
