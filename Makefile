# =============================================================================
# Makefile  --  QUICHOTTE (QUantum Interacting CHains: Out-of-equilibrium TransporT)
#
#   make                      build bin/quichotte.x     (platform auto-detected)
#   make PLATFORM=create      ... on KCL CREATE         (see make/create.mk)
#   make PLATFORM=archer2     ... on ARCHER2            (see make/archer2.mk)
#   make tools                helper programs in tools/ -> bin/
#   make test                 build, then run the regression tests in tests/
#   make update-references    overwrite the test references (deliberate physics changes ONLY)
#   make debug                build with bounds/NaN checking (slow; for hunting bugs)
#   make check-common         check common/ is identical to SANCHO/common/
#   make sancho               build SANCHO too (same as: make -C SANCHO)
#   make test-all             test QUICHOTTE and SANCHO, and run check-common
#   make clean                remove build/ and bin/ (never touches run output)
#
# Platform settings (compiler, flags, BLAS/LAPACK, FFTW) live in make/<name>.mk.
# Default: macos on a Mac, linux everywhere else.
# =============================================================================

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  PLATFORM ?= macos
else
  PLATFORM ?= linux
endif
include make/$(PLATFORM).mk

EXE    = bin/quichotte.x
BUILD  = build
export OMP_NUM_THREADS ?= 2

# Modules in dependency order (common/ first, then src/)
COMMON  = constants linear_algebra input_parser broyden_mixing scf_checks
SOURCES = parameters hamiltonian green_functions self_energies scf observables input_output quichotte
OBJ     = $(addprefix $(BUILD)/,$(addsuffix .o,$(COMMON) $(SOURCES)))
vpath %.f90 common src

.PHONY: all tools test update-references debug check-common sancho test-all clean help FORCE
all: $(EXE)

$(EXE): $(OBJ) | bin
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIBS)
	@echo "Built $@  (PLATFORM=$(PLATFORM), version $$(sed "s/.*= '\(.*\)'/\1/" $(BUILD)/version.inc))"

$(BUILD)/%.o: %.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(INCLUDES) -I$(BUILD) $(MODFLAG)$(BUILD) -c $< -o $@

# Code version stamped into the binary (git describe); rewritten only when it changes
$(BUILD)/version.inc: FORCE | $(BUILD)
	@v=$$(git describe --always --dirty --tags 2>/dev/null || echo unknown); \
	 echo "   character(len=*), parameter :: code_version = '$$v'" > $@.tmp; \
	 cmp -s $@.tmp $@ && rm -f $@.tmp || mv $@.tmp $@

# Module dependencies
$(BUILD)/linear_algebra.o:  $(BUILD)/constants.o
$(BUILD)/input_parser.o:    $(BUILD)/constants.o
$(BUILD)/broyden_mixing.o:  $(BUILD)/constants.o
$(BUILD)/scf_checks.o:      $(BUILD)/constants.o
$(BUILD)/parameters.o:      $(BUILD)/constants.o
$(BUILD)/hamiltonian.o:     $(BUILD)/linear_algebra.o $(BUILD)/parameters.o
$(BUILD)/green_functions.o: $(BUILD)/linear_algebra.o $(BUILD)/parameters.o $(BUILD)/hamiltonian.o
$(BUILD)/self_energies.o:   $(BUILD)/green_functions.o
$(BUILD)/scf.o:             $(BUILD)/broyden_mixing.o $(BUILD)/scf_checks.o $(BUILD)/self_energies.o
$(BUILD)/observables.o:     $(BUILD)/green_functions.o
$(BUILD)/input_output.o:    $(BUILD)/input_parser.o $(BUILD)/scf.o
$(BUILD)/quichotte.o:       $(BUILD)/input_output.o $(BUILD)/observables.o $(BUILD)/scf.o $(BUILD)/version.inc

# Helper programs
tools: bin/helix_xyz.x bin/dIdV.x
bin/%.x: tools/%.f90 | bin
	$(FC) $(FFLAGS) -o $@ $<

# Regression tests (the runner is itself a Fortran program)
bin/run_tests.x: tests/run_tests.f90 | bin
	$(FC) -O2 -o $@ $<
test: $(EXE) bin/run_tests.x
	@bin/run_tests.x "$(CURDIR)/$(EXE)" "$(CURDIR)/tests"
update-references: $(EXE) bin/run_tests.x
	@bin/run_tests.x "$(CURDIR)/$(EXE)" "$(CURDIR)/tests" -u

debug:
	$(MAKE) clean
	$(MAKE) FFLAGS="$(DEBUGFLAGS)"

# The shared files must stay identical in both codes
check-common:
	@ok=1; for f in common/*.f90 tests/run_tests.f90; do \
	   cmp -s $$f SANCHO/$$f || { echo "DIFFERS: $$f vs SANCHO/$$f"; ok=0; }; done; \
	 for f in SANCHO/common/*.f90; do [ -f $${f#SANCHO/} ] || { echo "only in SANCHO: $$f"; ok=0; }; done; \
	 [ $$ok = 1 ] && echo "common/ and SANCHO/common/ are identical" || exit 1

sancho:
	$(MAKE) -C SANCHO PLATFORM=$(PLATFORM)
test-all: test check-common
	$(MAKE) -C SANCHO PLATFORM=$(PLATFORM) test

$(BUILD) bin:
	@mkdir -p $@

clean:
	rm -rf $(BUILD) bin

help:
	@sed -n '2,18p' Makefile
