SHELL=/bin/bash
ARCH=$(shell uname)
ifndef ABI
ABI=ilp64
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
ifdef SLEEF
CXX=icpx
endif # SLEEF
CPUFLAGS=-DKNL -fPIC -fexceptions -fno-omit-frame-pointer -rdynamic
ifdef NDEBUG
ifdef MKL
ifeq ($(MKL),intel_thread)
CPUFLAGS += -qopenmp
endif # intel_thread
else # !MKL
CPUFLAGS += -qopenmp
endif # ?MKL
SUFX=-$(ABI)_$(NDEBUG)$(WP)
else # DEBUG
SUFX=-$(ABI)_$(DEBUG)$(WP)
endif # ?NDEBUG
ifndef MKL
ifndef LAPACK
MKL=sequential
endif # !LAPACK
endif # !MKL
DBGFLAGS=-traceback -diag-disable=10397 -DJTRACE
FPUFLAGS=-fp-model $(FPU) -fprotect-parens -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -qsimd-honor-fp-model -qsimd-serialize-fp-reduction #-fimf-use-svml=true
ifeq ($(WP),l)
FPUFLAGS += -DUSE_EXTENDED
endif # ?WP
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -xHost -qopt-multi-version-aggressive -qopt-zmm-usage=high -inline-level=2 -vec-threshold0
DBGFLAGS += -DNDEBUG -qopt-report=5
else # DEBUG
OPTFLAGS=-O0 -xHost -qopt-multi-version-aggressive -qopt-zmm-usage=high
DBGFLAGS += -$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug parallel -debug pubnames -DPRINTOUT=stderr
FPUFLAGS += -fp-stack-check
endif # ?NDEBUG
LIBFLAGS=-static-libgcc -D_GNU_SOURCE -D_LARGEFILE64_SOURCE -I. -I../../JACSD/jstrat -DUSE_INL -DUSE_2SUM #-DUSE_SECANTS
ifdef SLEEF
LIBFLAGS += -DDZNRME_SEQRED -DSCNRME_SEQRED -DUSE_SLEEF -I$(SLEEF)/include
endif # SLEEF
ifdef CR_MATH
LIBFLAGS += -DUSE_CR_MATH
endif # CR_MATH
ifndef LAPACK
LIBFLAGS += -DUSE_MKL -I${MKLROOT}/include/intel64/$(ABI) -I${MKLROOT}/include
endif # MKL
ifeq ($(ABI),ilp64)
LIBFLAGS += -DMKL_ILP64
endif # ilp64
ifdef CR_MATH
LDFLAGS=$(CR_MATH)/src/binary32/hypot/hypotf.o $(CR_MATH)/src/binary64/hypot/hypot.o -L. -lvecjac$(SUFX)
else # !CR_MATH
LDFLAGS=
endif # ?CR_MATH
LDFLAGS += -L. -lvecjac$(SUFX) -L../../JACSD -ljstrat$(DEBUG)
ifdef SLEEF
LDFLAGS += -L$(SLEEF)/lib64 -Wl,-rpath=$(SLEEF)/lib64 -lsleefquad
endif # SLEEF
ifdef LAPACK
LDFLAGS += -L$(LAPACK) -ltmglib -llapack -lrefblas -lifcoremt
else # MKL
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_$(ABI) -lmkl_$(MKL) -lmkl_core
endif # ?LAPACK
ifdef MKL
ifeq ($(MKL),intel_thread)
ifeq ($(findstring qopenmp,$(CPUFLAGS)),)
LDFLAGS += -liomp5
endif # ?qopenmp
endif # intel_thread
endif # MKL
LDFLAGS += -lpthread -lm -ldl
CFLAGS=-std=gnu18 $(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS)
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS) -standard-semantics -recursive -threads
ifndef NDEBUG
CFLAGS += -check=stack,uninit
FFLAGS += -check all -assume ieee_fpe_flags
endif # !NDEBUG
ifdef SLEEF
CXXFLAGS=-std=gnu++20 -qtbb $(OPTFLAGS) $(subst -qopt-report=5,,$(subst -debug pubnames,,$(DBGFLAGS))) $(LIBFLAGS) $(CPUFLAGS) $(subst -qsimd-serialize-fp-reduction,,$(subst -qsimd-honor-fp-model,,$(subst -no-ftz,,$(FPUFLAGS))))
ifndef NDEBUG
CXXFLAGS += -fcheck=stack,uninit
endif # !NDEBUG
endif # SLEEF
ifeq ($(ABI),ilp64)
FFLAGS += -i8
endif # ilp64
FLFLAGS=-lifcoremt
