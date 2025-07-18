!IFNDEF NVCC
NVCC=nvcc.exe
!ENDIF # !NVCC
!IFNDEF SM
SM=native
!ENDIF # !SM
NVCCFLAGS=-std=c++20 -expt-relaxed-constexpr -extended-lambda -restrict -arch=$(SM) -extra-device-vectorization -ewp -res-usage -Wno-deprecated-gpu-targets
!IFDEF NDEBUG
NVCCFLAGS=$(NVCCFLAGS) -O$(NDEBUG) -opt-info=inline -dlto
!ELSE # !NDEBUG
NVCCFLAGS=$(NVCCFLAGS) -g -G -v -keep
!ENDIF # ?NDEBUG
NVCCFLAGS=$(NVCCFLAGS) -Xptxas=-v,-warn-lmem-usage,-warn-spills

all: tcubjac2.exe

common.obj: common.c common.h Makefile
	$(NVCC) $(NVCCFLAGS) -c common.c

wnrme.obj: wnrme.c wnrme.h common.h Makefile
	$(NVCC) $(NVCCFLAGS) -c wnrme.c

tcubjac2.exe: tcubjac2.cu device_code.h timer.h wnrme.obj common.obj Makefile
	$(NVCC) $(NVCCFLAGS) tcubjac2.cu wnrme.obj common.obj -o $@

clean:
	-del /f /q tcubjac2.exe wnrme.obj common.obj
