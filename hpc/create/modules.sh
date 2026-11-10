#!/bin/bash
# source hpc/create/modules.sh   -- modules to BUILD and RUN on KCL CREATE
# (module names as of Nov 2026; check with: module spider gcc / openblas / fftw)
module purge
module load gcc/12.5.0
module load openblas/0.3.33
module load fftw/3.3.11
module list

# The spack modules set LIBRARY_PATH (link time) but not LD_LIBRARY_PATH (run time):
# without this the binary links, then fails with "cannot open libfftw3.so.3"
export LD_LIBRARY_PATH=${LIBRARY_PATH}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

# make/create.mk takes the FFTW include path from CPATH or C_INCLUDE_PATH
if [[ ":$CPATH:$C_INCLUDE_PATH:" != *fftw* ]]; then
   echo "WARNING: no fftw directory on CPATH or C_INCLUDE_PATH -- build with:"
   echo "         make PLATFORM=create FFTW_DIR=<fftw prefix>   (prefix: module show fftw/3.3.11)"
fi
