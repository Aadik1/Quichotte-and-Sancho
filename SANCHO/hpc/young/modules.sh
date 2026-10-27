#!/bin/bash
# source hpc/young/modules.sh   -- modules to BUILD and RUN on Young (MMM Hub)
module purge
module load gcc-libs/10.2.0
module load openblas/0.3.13-openmp/gnu-10.2.0
module list
