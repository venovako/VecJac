SHELL=/bin/bash
ARCH=$(shell uname)
ifndef ABI
ABI=ilp64
endif # !ABI
ifndef WP
WP=l
endif # !WP
SUFX=-$(ABI)_$(WP)
RM=rm -rfv
AR=ar
ARFLAGS=rsv
include ../../libpvn/src/pvn.mk
CC=$(PVN_CC)
CXX=$(PVN_CXX)
FC=$(PVN_FC)
CPPFLAGS=$(PVN_CPPFLAGS)
CFLAGS=$(PVN_CFLAGS)
CXXFLAGS=$(PVN_CXXFLAGS)
FFLAGS=$(PVN_FCFLAGS)
LDFLAGS=$(PVN_LDFLAGS)
ifndef MKL
ifndef LAPACK
MKL=sequential
endif # !LAPACK
endif # !MKL
CPPFLAGS += -DJTRACE -I. -DUSE_INL -DUSE_2SUM #-DUSE_SECANTS
ifeq ($(WP),l)
CPPFLAGS += -DUSE_EXTENDED
endif # ?WP
ifeq ($(findstring PRINTOUT,$(PVN_CPPFLAGS)),PRINTOUT)
CPPFLAGS += -DPRINTOUT=stderr
endif # PRINTOUT
ifeq ($(findstring SLEEF,$(PVN_CPPFLAGS)),SLEEF)
CPPFLAGS += -DDZNRME_SEQRED -DSCNRME_SEQRED -DUSE_SLEEF
ifdef PVN_CXX
SLEEF=cpp
else # !PVN_CXX
SLEEF=c
endif # ?PVN_CXX
endif # SLEEF
ifndef LAPACK
CPPFLAGS += -DUSE_MKL -I${MKLROOT}/include/intel64/$(ABI) -I${MKLROOT}/include
endif # MKL
ifeq ($(ABI),ilp64)
CPPFLAGS += -DMKL_ILP64
endif # ilp64
ifeq ($(findstring CR_MATH,$(PVN_CPPFLAGS)),CR_MATH)
CPPFLAGS += -DUSE_CR_MATH
endif # CR_MATH
LDFLAGS += -L. -lvecjac$(SUFX)
ifdef LAPACK
LDFLAGS += -L$(LAPACK) -ltmglib -llapack -lrefblas
else # MKL
LDFLAGS += -L${MKLROOT}/lib -Wl,-rpath=${MKLROOT}/lib -lmkl_gf_$(ABI) -lmkl_$(MKL) -lmkl_core
endif # ?LAPACK
LDFLAGS += $(PVN_LIBS)
CFLAGS += $(CPPFLAGS)
CXXFLAGS += $(CPPFLAGS)
ifeq ($(ABI),ilp64)
FFLAGS += -fdefault-integer-8
endif # ilp64
FFLAGS += $(CPPFLAGS)
FLFLAGS=-lgfortran
