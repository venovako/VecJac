CC=icx
ifndef CPU
CPU=Host
endif # !CPU
ifndef NDEBUG
NDEBUG=3
endif # !NDEBUG
CFLAGS=-D_GNU_SOURCE -D_LARGEFILE64_SOURCE -DNDEBUG -std=gnu18 -O$(NDEBUG) -x$(CPU) -fno-math-errno -inline-level=2 -qopt-report=3 -qopenmp -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -fp-model=precise -fp-speculation=safe -fprotect-parens -fma -no-ftz -fimf-precision=high -mprefer-vector-width=512 -traceback -vec-threshold0
LDFLAGS=-rdynamic -static-libgcc -lpthread -lm -ldl

all: tsjac2.exe

common.o: common.c common.h GNUmakefile
	$(CC) $(CFLAGS) -c common.c

wnrme.o: wnrme.c wnrme.h common.h GNUmakefile
	$(CC) $(CFLAGS) -c wnrme.c

serial.o: serial.c serial.h GNUmakefile
	$(CC) $(CFLAGS) -c serial.c

tsjac2.exe: tsjac2.c timer.h serial.o wnrme.o common.o GNUmakefile
	$(CC) $(CFLAGS) tsjac2.c serial.o wnrme.o common.o -o $@ $(LDFLAGS)

clean:
	-rm -rfv tsjac2.exe tsjac2.opt* serial.o serial.opt* wnrme.o wnrme.opt* common.o common.opt*
