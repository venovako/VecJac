CC=icx.exe
!IFNDEF MARCH
MARCH=Host
!ENDIF # !MARCH
!IFNDEF NDEBUG
NDEBUG=3
!ENDIF # !NDEBUG
CFLAGS=/nologo /DNDEBUG /Qstd=c18 /O$(NDEBUG) /Qx$(MARCH) /Qlong-double /Qopenmp /Qopt-report:3 /fp:precise /Qprotect-parens /Qfma /Qftz- /Ob2 /Qvec-threshold:0 /Wno-deprecated-declarations /traceback

all: tsjac2.exe

common.obj: common.c common.h Makefile
	$(CC) $(CFLAGS) /c common.c

wnrme.obj: wnrme.c wnrme.h common.h Makefile
	$(CC) $(CFLAGS) /c wnrme.c

serial.obj: serial.c serial.h Makefile
	$(CC) $(CFLAGS) /c serial.c

tsjac2.exe: tsjac2.c timer.h serial.obj wnrme.obj common.obj Makefile
	$(CC) $(CFLAGS) tsjac2.c serial.obj wnrme.obj common.obj /Fe$@ /link /RELEASE

clean:
	-del /f /q tsjac2.exe tsjac2.optrpt serial.obj serial.optrpt wnrme.obj wnrme.optrpt common.obj common.optrpt
