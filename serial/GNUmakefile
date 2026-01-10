# MARCH=common-avx512 for KNLs
ifndef LIBPVN
LIBPVN=../../libpvn
endif # !LIBPVN
include $(LIBPVN)/src/pvn.mk
CC=$(PVN_CC)
CFLAGS=$(PVN_CFLAGS) $(PVN_CPPFLAGS)
LDFLAGS=$(PVN_LDFLAGS) $(PVN_LIBS)

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
	-rm -rfv *.exe *.o *.opt*
