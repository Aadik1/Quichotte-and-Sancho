# make/macos.mk -- macOS with Homebrew
#   brew install gcc openblas
FC         = gfortran
FFLAGS     = -O2 -fopenmp
DEBUGFLAGS = -O0 -g -fopenmp -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall
MODFLAG    = -J
OPENBLAS_DIR ?= $(shell brew --prefix openblas 2>/dev/null || echo /usr/local/opt/openblas)
INCLUDES   =
LIBS       = -L$(OPENBLAS_DIR)/lib -lopenblas
