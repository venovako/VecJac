ifndef NVCC
NVCC=nvcc
endif # !NVCC
ifndef SM
SM=52
endif # !SM
NVCCFLAGS=-std=c++17 -expt-relaxed-constexpr -extended-lambda -restrict -arch=sm_$(SM) -extra-device-vectorization -ewp -res-usage -D_LARGEFILE64_SOURCE
ifdef NDEBUG
NVCCFLAGS += -O$(NDEBUG) -opt-info=inline -dlto
else # !NDEBUG
NVCCFLAGS += -g -G -v -keep
endif # ?NDEBUG
ifndef MARCH
MARCH=native
endif # !MARCH
NVCCFLAGS += -Xcompiler=-fopenmp,-march=$(MARCH) -Xptxas=-v,-warn-lmem-usage,-warn-spills

all: tcubjac2.exe

common.o: common.c common.h GNUmakefile
	$(NVCC) $(NVCCFLAGS) -c common.c

wnrme.o: wnrme.c wnrme.h common.h GNUmakefile
	$(NVCC) $(NVCCFLAGS) -c wnrme.c

tcubjac2.exe: tcubjac2.cu device_code.h timer.h wnrme.o common.o GNUmakefile
	$(NVCC) $(NVCCFLAGS) tcubjac2.cu wnrme.o common.o -o $@

clean:
	-rm -rfv tcubjac2.exe wnrme.o common.o
