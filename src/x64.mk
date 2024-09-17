SHELL=/bin/bash
ARCH=$(shell uname)
ifndef CPU
CPU=Host
# COMMON-AVX512 for KNLs
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
CPUFLAGS=-fPIC -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -qopt-multi-version-aggressive -qopt-zmm-usage=high -vec-threshold0
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
DBGFLAGS=-traceback -diag-disable=10397,10441 -DJTRACE
FPUFLAGS=-fp-model $(FPU) -fprotect-parens -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -qsimd-honor-fp-model -qsimd-serialize-fp-reduction
ifeq ($(WP),l)
FPUFLAGS += -DUSE_EXTENDED
endif # ?WP
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -x$(CPU) -fno-math-errno -inline-level=2
DBGFLAGS += -DNDEBUG -qopt-report=5
else # DEBUG
OPTFLAGS=-O0 -x$(CPU)
DBGFLAGS += -$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -DPRINTOUT=stderr -D__INTEL_COMPILER_USE_INTRINSIC_PROTOTYPES
ifneq ($(ARCH),Darwin) # Linux
DBGFLAGS += -debug parallel
endif # Linux
FPUFLAGS += -fp-stack-check
endif # ?NDEBUG
LIBFLAGS=-I. -I../../JACSD/jstrat -DUSE_INL -DUSE_2SUM #-DUSE_SECANTS
LDFLAGS=-rdynamic -static-libgcc
ifdef SLEEF
LIBFLAGS += -DDZNRME_SEQRED -DSCNRME_SEQRED -DUSE_SLEEF -I$(SLEEF)/include
endif # SLEEF
ifndef LAPACK
LIBFLAGS += -DUSE_MKL -I${MKLROOT}/include/intel64/$(ABI) -I${MKLROOT}/include
endif # MKL
ifeq ($(ABI),ilp64)
LIBFLAGS += -DMKL_ILP64
endif # ilp64
LDFLAGS += -L. -lvecjac$(SUFX) -L../../JACSD -ljstrat$(DEBUG)
ifdef CR_MATH
LIBFLAGS += -DUSE_CR_MATH
LDFLAGS += $(CR_MATH)/src/binary32/hypot/hypotf.o $(CR_MATH)/src/binary64/hypot/hypot.o
endif # CR_MATH
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
LDFLAGS += -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_$(ABI) -lmkl_$(MKL) -lmkl_core
else # Linux
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_$(ABI) -lmkl_$(MKL) -lmkl_core
endif # ?Darwin
endif # ?LAPACK
ifneq ($(ARCH),Darwin) # Linux
LIBFLAGS += -D_GNU_SOURCE -D_LARGEFILE64_SOURCE
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
FLFLAGS=-lifcoremt
