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
RM=rm -rfv
AR=xiar
ARFLAGS=-qnoipo -lib rsv
CC=icc
CPUFLAGS=-fPIC -fexceptions -fno-omit-frame-pointer -qopenmp -rdynamic -DTSC_FREQ_HZ=0ull
ifdef NDEBUG
SUFX=-$(ABI)_$(NDEBUG)
else # DEBUG
SUFX=-$(ABI)_$(DEBUG)
endif # ?NDEBUG
DBGFLAGS=-traceback -w3
FPUFLAGS=-fp-model $(FPU) -fprotect-parens -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fimf-use-svml=true
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -xHost -qopt-multi-version-aggressive
DBGFLAGS += -DNDEBUG -qopt-report=5
else # DEBUG
OPTFLAGS=-O0 -xHost -qopt-multi-version-aggressive
DBGFLAGS += -$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -check=stack,uninit
ifneq ($(ARCH),Darwin) # Linux
DBGFLAGS += -debug parallel
endif # ?Linux
FPUFLAGS += -fp-stack-check
endif # ?NDEBUG
LIBFLAGS=-I. -DUSE_MKL
ifeq ($(ABI),ilp64)
LIBFLAGS += -DMKL_ILP64
endif # ilp64
LIBFLAGS += -I${MKLROOT}/include/intel64/$(ABI) -I${MKLROOT}/include
LDFLAGS=-L. -lvecjac$(SUFX)
ifeq ($(ARCH),Darwin)
LDFLAGS += -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_$(ABI) -lmkl_sequential -lmkl_core
else # Linux
LIBFLAGS += -static-libgcc -D_GNU_SOURCE
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_$(ABI) -lmkl_sequential -lmkl_core
endif # ?Darwin
LDFLAGS += -lpthread -lm -ldl
CFLAGS=-std=c18 $(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(CPUFLAGS) $(FPUFLAGS)
