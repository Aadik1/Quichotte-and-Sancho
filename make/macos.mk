# make/macos.mk -- macOS with Homebrew
#   brew install gcc openblas fftw
FC         = gfortran
FFLAGS     = -O2 -fopenmp
DEBUGFLAGS = -O0 -g -fopenmp -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall
MODFLAG    = -J
FFTW_DIR     ?= $(shell brew --prefix fftw 2>/dev/null || echo /usr/local/opt/fftw)
OPENBLAS_DIR ?= $(shell brew --prefix openblas 2>/dev/null || echo /usr/local/opt/openblas)
INCLUDES   = -I$(FFTW_DIR)/include
LIBS       = -L$(FFTW_DIR)/lib -lfftw3 -L$(OPENBLAS_DIR)/lib -lopenblas
