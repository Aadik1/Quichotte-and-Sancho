#!/bin/bash
# source hpc/create/modules.sh   -- modules to BUILD and RUN on KCL CREATE
module purge
module load gcc/13.2.0
module load openblas/0.3.24-gcc-13.2.0-multithreading-openmp
module list
