SHELL=/bin/bash
ARCH=$(shell uname)
ifndef CPU
CPU=Host
endif # !CPU
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
CC=icx
FC=ifx
ifdef SLEEF
CXX=icpx
endif # SLEEF
CPUFLAGS=-fPIC -fexceptions -fno-omit-frame-pointer
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
DBGFLAGS=-traceback -DJTRACE
FPUFLAGS=-fp-model=$(FPU) -fp-speculation=safe -fprotect-parens -fma -no-ftz -fimf-precision=high #-qsimd-honor-fp-model -qsimd-serialize-fp-reduction
ifeq ($(WP),l)
FPUFLAGS += -DUSE_EXTENDED
endif # ?WP
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -x$(CPU) -inline-level=2 -mprefer-vector-width=512 -vec-threshold0
DBGFLAGS += -DNDEBUG -qopt-report=3
else # DEBUG
OPTFLAGS=-O0 -x$(CPU) -mprefer-vector-width=512 -vec-threshold0
DBGFLAGS += -$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -debug parallel -DPRINTOUT=stderr
endif # ?NDEBUG
LIBFLAGS=-I. -I../../JACSD/jstrat -DUSE_INL -DUSE_2SUM #-DUSE_SECANTS
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
LDFLAGS=-rdynamic -static-libgcc
ifdef CR_MATH
LDFLAGS += $(CR_MATH)/src/binary32/hypot/hypotf.o $(CR_MATH)/src/binary64/hypot/hypot.o
endif # CR_MATH
LDFLAGS += -L. -lvecjac$(SUFX) -L../../JACSD -ljstrat$(DEBUG)
ifdef SLEEF
ifeq ($(wildcard $(SLEEF)/lib64),)
LDFLAGS += -L$(SLEEF)/lib -Wl,-rpath=$(SLEEF)/lib
else # lib64
LDFLAGS += -L$(SLEEF)/lib64 -Wl,-rpath=$(SLEEF)/lib64
endif # ?lib64
LDFLAGS += -lsleefquad
endif # SLEEF
ifdef LAPACK
LDFLAGS += -L$(LAPACK) -ltmglib -llapack -lrefblas -lifcoremt
else # MKL
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_$(ABI) -lmkl_$(MKL) -lmkl_core
endif # ?LAPACK
LIBFLAGS += -D_GNU_SOURCE -D_LARGEFILE64_SOURCE
ifdef MKL
ifeq ($(MKL),intel_thread)
ifeq ($(findstring qopenmp,$(CPUFLAGS)),)
LDFLAGS += -liomp5
endif # ?qopenmp
endif # intel_thread
endif # MKL
LDFLAGS += -lpthread -lm -ldl
CFLAGS=-std=gnu18 $(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS) -Wno-parentheses -Wno-unused-command-line-argument -Wno-pointer-sign -Wno-incompatible-pointer-types
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS) -standard-semantics -recursive -threads
ifdef SLEEF
CXXFLAGS=-std=gnu++20 -qtbb $(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS) -Wno-unused-command-line-argument
endif # SLEEF
ifndef NDEBUG
CFLAGS += -fcheck=stack,uninit
FFLAGS += -check all -assume ieee_fpe_flags
endif # !NDEBUG
ifdef SLEEF
ifndef NDEBUG
CXXFLAGS += -check=stack,uninit
endif # !NDEBUG
endif # SLEEF
ifeq ($(ABI),ilp64)
FFLAGS += -i8
endif # ilp64
FLFLAGS=-lifcoremt
