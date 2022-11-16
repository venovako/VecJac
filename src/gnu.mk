ifdef SLEEF
$(error SLEEF is not yet supported with GCC)
endif # SLEEF
SHELL=/bin/bash
ARCH=$(shell uname)
ifneq ($(ARCH),Linux)
$(error GCC build is only supported on Linux)
endif # !Linux
ifndef CPU
CPU=native
endif # !CPU
ifndef ABI
ABI=ilp64
endif # !ABI
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
ifndef WP
WP=l
endif # !WP
RM=rm -rfv
AR=ar
ARFLAGS=rsv
CC=gcc$(GNU)
FC=gfortran$(GNU)
CPUFLAGS=-fPIC -fexceptions -fno-omit-frame-pointer -rdynamic
ifdef NDEBUG
ifdef MKL
ifeq ($(MKL),gnu_thread)
CPUFLAGS += -fopenmp
endif # gnu_thread
else # !MKL
CPUFLAGS += -fopenmp
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
OPTFLAGS=-march=$(CPU)
DBGFLAGS=-DJTRACE
FPUFLAGS=-DUSE_EXTENDED -ffp-contract=fast
ifdef NDEBUG
OPTFLAGS += -O$(NDEBUG) -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller
DBGFLAGS += -DNDEBUG -fopt-info-optimized-vec
FPUFLAGS += -fno-math-errno
else # DEBUG
OPTFLAGS += -O$(DEBUG)
DBGFLAGS += -$(DEBUG) -DPRINTOUT=stderr
endif # ?NDEBUG
LIBFLAGS=-I. -I../../JACSD/jstrat -DUSE_INL -DUSE_2SUM #-DUSE_SECANTS
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
ifdef LAPACK
LDFLAGS += -L$(LAPACK) -ltmglib -llapack -lrefblas
else # MKL
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_gf_$(ABI) -lmkl_$(MKL) -lmkl_core
endif # ?LAPACK
LIBFLAGS += -static-libgcc -D_GNU_SOURCE -D_LARGEFILE64_SOURCE
ifdef MKL
ifeq ($(MKL),gnu_thread)
ifeq ($(findstring fopenmp,$(CPUFLAGS)),)
LDFLAGS += -lgomp
endif # ?qopenmp
endif # intel_thread
endif # MKL
LDFLAGS += -lpthread -lm -ldl
CFLAGS=-std=gnu18 $(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS) -Wno-incompatible-pointer-types
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS) -fprotect-parens
ifndef NDEBUG
FFLAGS += -fcheck=all -finit-local-zero -finit-real=snan -finit-derived -pedantic
endif # DEBUG
ifeq ($(ABI),ilp64)
FFLAGS += -fdefault-integer-8
endif # ilp64
FLFLAGS=-lgfortran
