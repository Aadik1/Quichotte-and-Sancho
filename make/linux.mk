# make/linux.mk -- generic Linux workstation with gfortran
#   Ubuntu/Debian:  sudo apt install gfortran libopenblas-dev libfftw3-dev
FC         = gfortran
FFLAGS     = -O2 -fopenmp
DEBUGFLAGS = -O0 -g -fopenmp -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall
MODFLAG    = -J
FFTW_DIR  ?= /usr
INCLUDES   = -I$(FFTW_DIR)/include
LIBS       = -L$(FFTW_DIR)/lib -lfftw3 -lopenblas
