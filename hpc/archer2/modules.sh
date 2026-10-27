#!/bin/bash
# source hpc/archer2/modules.sh   -- modules to BUILD and RUN QUICHOTTE on ARCHER2
module restore
module load PrgEnv-gnu        # ftn then links Cray LibSci (BLAS/LAPACK) itself
module load cray-fftw
module list
