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
CC=icc
FC=ifort
ifdef SLEEF
CXX=icpc
endif # SLEEF
CPUFLAGS=-fPIC -fexceptions -fno-omit-frame-pointer -rdynamic
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
DBGFLAGS=-traceback #-DJTRACE
FPUFLAGS=-fp-model $(FPU) -fprotect-parens -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt #-qsimd-honor-fp-model -qsimd-serialize-fp-reduction -fimf-use-svml=true
ifeq ($(WP),l)
FPUFLAGS += -DUSE_EXTENDED
endif # ?WP
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -x$(CPU) -qopt-multi-version-aggressive -inline-level=2 -vec-threshold0
DBGFLAGS += -DNDEBUG -qopt-report=5
else # DEBUG
OPTFLAGS=-O0 -x$(CPU) -qopt-multi-version-aggressive
DBGFLAGS += -$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -DPRINTOUT=stderr
ifneq ($(ARCH),Darwin) # Linux
DBGFLAGS += -debug parallel
endif # Linux
FPUFLAGS += -fp-stack-check
endif # ?NDEBUG
LIBFLAGS=-I. -I../../JACSD/jstrat -DUSE_INL
ifdef SLEEF
LIBFLAGS += -DDZNRME_SEQRED -DSCNRME_SEQRED -DUSE_SLEEF -I$(SLEEF)/include
endif # SLEEF
ifndef LAPACK
LIBFLAGS += -DUSE_MKL -I${MKLROOT}/include/intel64/$(ABI) -I${MKLROOT}/include
endif # MKL
ifeq ($(ABI),ilp64)
LIBFLAGS += -DMKL_ILP64
endif # ilp64
LDFLAGS=-L. -lvecjac$(SUFX) -L../../JACSD -ljstrat$(DEBUG)
ifdef SLEEF
ifdef ($(ARCH),Darwin)
LDFLAGS += -L$(SLEEF)/lib -Wl,-rpath,$(SLEEF)/lib
else # Linux
ifeq ($(wildcard $(SLEEF)/lib64),)
LDFLAGS += -L$(SLEEF)/lib -Wl,-rpath=$(SLEEF)/lib
else # lib64
LDFLAGS += -L$(SLEEF)/lib64 -Wl,-rpath=$(SLEEF)/lib64
endif # ?lib64
endif # ?Darwin
LDFLAGS += -lsleefquad
endif # SLEEF
ifdef LAPACK
LDFLAGS += -L$(LAPACK) -ltmglib -llapack -lrefblas -lifcoremt
else # MKL
ifeq ($(ARCH),Darwin)
LDFLAGS += -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_$(ABI) -lmkl_sequential -lmkl_core
else # Linux
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_$(ABI) -lmkl_sequential -lmkl_core
endif # ?Darwin
endif # ?LAPACK
ifneq ($(ARCH),Darwin) # Linux
LIBFLAGS += -static-libgcc -D_GNU_SOURCE -D_LARGEFILE64_SOURCE
endif # Linux
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
ifdef SLEEF
CXXFLAGS=-std=gnu++20 -qtbb $(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS)
endif # SLEEF
ifndef NDEBUG
CFLAGS += -check=stack,uninit
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
