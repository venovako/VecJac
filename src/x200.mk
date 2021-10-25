SHELL=/bin/bash
ARCH=$(shell uname)
ifndef ABI
ABI=lp64
endif # !ABI
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
ifndef FPU
FPU=precise
endif # !FPU
ifndef WP
WP=q
endif # !WP
RM=rm -rfv
AR=xiar
ARFLAGS=-qnoipo -lib rsv
CC=icc
FC=ifort
CXX=icpx
CPUFLAGS=-DKNL -fPIC -fexceptions -fno-omit-frame-pointer -rdynamic
ifdef NDEBUG
CPUFLAGS += -qopenmp
SUFX=-$(ABI)_$(NDEBUG)$(WP)
else # DEBUG
SUFX=-$(ABI)_$(DEBUG)$(WP)
endif # ?NDEBUG
DBGFLAGS=-traceback
FPUFLAGS=-fp-model $(FPU) -fprotect-parens -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fimf-use-svml=true
ifeq ($(WP),l)
FPUFLAGS += -DUSE_EXTENDED
endif # ?WP
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -xHost -qopt-multi-version-aggressive -qopt-zmm-usage=high
DBGFLAGS += -DNDEBUG -qopt-report=5
else # DEBUG
OPTFLAGS=-O0 -xHost -qopt-multi-version-aggressive -qopt-zmm-usage=high
DBGFLAGS += -$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug parallel -debug pubnames -DPRINTOUT=stderr
FPUFLAGS += -fp-stack-check
endif # ?NDEBUG
LIBFLAGS=-static-libgcc -D_GNU_SOURCE -D_LARGEFILE64_SOURCE -I. -I../../JACSD/jstrat
ifndef LAPACK
LIBFLAGS += -DUSE_MKL
ifeq ($(ABI),ilp64)
LIBFLAGS += -DMKL_ILP64
endif # ilp64
LIBFLAGS += -I${MKLROOT}/include/intel64/$(ABI) -I${MKLROOT}/include
endif # MKL
LDFLAGS=-L. -lvecjac$(SUFX) -L../../JACSD -ljstrat$(DEBUG)
ifdef LAPACK
LDFLAGS += -L$(LAPACK) -ltmglib -llapack -lrefblas -lifcoremt
else # MKL
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_$(ABI) -lmkl_sequential -lmkl_core
endif # ?LAPACK
LDFLAGS += -lpthread -lm -ldl -lmemkind
CFLAGS=-std=c18 $(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS)
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS) -standard-semantics -threads -assume ieee_fpe_flags
CXXFLAGS=-std=gnu++20 -qtbb $(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS)
ifdef NDEBUG
CFLAGS += -w3
CXXFLAGS += -w3
else # DEBUG
CFLAGS += -check=stack,uninit -w3
FFLAGS += -check all
CXXFLAGS += -check=stack,uninit -w3
endif # ?NDEBUG
ifeq ($(ABI),ilp64)
FFLAGS += -i8
endif # ilp64
