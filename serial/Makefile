!INCLUDE ..\..\libpvn\src\pvn.mk
CC=$(PVN_CC)
CFLAGS=$(PVN_CFLAGS) $(PVN_CPPFLAGS)
LDFLAGS=$(PVN_LDFLAGS) $(PVN_LIBS)

all: tsjac2.exe

common.obj: common.c common.h Makefile
	$(CC) $(CFLAGS) /c common.c

wnrme.obj: wnrme.c wnrme.h common.h Makefile
	$(CC) $(CFLAGS) /c wnrme.c

serial.obj: serial.c serial.h Makefile
	$(CC) $(CFLAGS) /c serial.c

tsjac2.exe: tsjac2.c timer.h serial.obj wnrme.obj common.obj Makefile
	$(CC) $(CFLAGS) tsjac2.c serial.obj wnrme.obj common.obj /Fe$@ $(LDFLAGS)

clean:
	-del /f /q tsjac2.exe tsjac2.optrpt serial.obj serial.optrpt wnrme.obj wnrme.optrpt common.obj common.optrpt
