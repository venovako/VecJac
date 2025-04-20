CC=icx.exe
!IFNDEF MARCH
MARCH=Host
!ENDIF # !MARCH
!IFNDEF NDEBUG
NDEBUG=3
!ENDIF # !NDEBUG
CFLAGS=/nologo /Qstd=c18 /Qlong-double /Qoption,cpp,--extended_float_types /DNDEBUG=$(NDEBUG) /I. /EHsc /O$(NDEBUG) /Qx$(MARCH) /Ob2 /Qopenmp /fp:precise /Qfma /Qftz- /Qopt-report:3 /Qvec-threshold:0 /traceback /Wno-deprecated-declarations /Wno-overriding-option

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
