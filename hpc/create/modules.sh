#!/bin/bash
# source hpc/create/modules.sh   -- modules to BUILD and RUN on KCL CREATE
# (module names as of Nov 2026; check with: module spider gcc / openblas / fftw)
module purge
module load gcc/12.5.0
module load openblas/0.3.33
module load fftw/3.3.11
module list

# make/create.mk takes the FFTW include path from CPATH: warn if it is missing
if [[ ":$CPATH:" != *fftw* ]]; then
   echo "WARNING: no fftw directory on CPATH -- build with: make PLATFORM=create FFTW_DIR=<fftw prefix>"
   echo "         (the prefix is shown by: module show fftw/3.3.11)"
fi
