CC=icc
ifndef CPU
CPU=Host
endif # !CPU
ifndef NDEBUG
NDEBUG=3
endif # !NDEBUG
CFLAGS=-D_LARGEFILE64_SOURCE -DNDEBUG -std=c18 -O$(NDEBUG) -x$(CPU) -qopenmp -qopt-multi-version-aggressive -qopt-report=5 -fp-model precise -fprotect-parens -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fPIC -fexceptions -fno-omit-frame-pointer -rdynamic -traceback
LDFLAGS=-lpthread -lm -ldl
ifneq ($(shell uname),Darwin) # Linux
CFLAGS += -D_GNU_SOURCE
LDFLAGS += -static-libgcc
endif # Linux

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
	-rm -rfv tsjac2.exe tsjac2.optrpt serial.o serial.optrpt wnrme.o wnrme.optrpt common.o common.optrpt
