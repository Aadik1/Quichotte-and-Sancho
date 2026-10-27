# make/linux.mk -- generic Linux workstation with gfortran
#   Ubuntu/Debian:  sudo apt install gfortran libopenblas-dev
FC         = gfortran
FFLAGS     = -O2 -fopenmp
DEBUGFLAGS = -O0 -g -fopenmp -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall
MODFLAG    = -J
INCLUDES   =
LIBS       = -lopenblas
