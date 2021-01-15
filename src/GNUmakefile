ifndef COMPILER
COMPILER=x64
endif # !COMPILER
include $(COMPILER).mk
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean

all: libvecjac$(SUFX).a

help:
	@echo "make [COMPILER=x64|x200] [NDEBUG=0|1|2|3|4|5] [ABI=lp64|ilp64] [all|clean|help]"

libvecjac$(SUFX).a: common$(SUFX).o rnd$(SUFX).o vec$(SUFX).o $(MKFS)
	$(AR) $(ARFLAGS) $@ common$(SUFX).o rnd$(SUFX).o vec$(SUFX).o

common$(SUFX).o: common.c common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

rnd$(SUFX).o: rnd.c rnd.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

vec$(SUFX).o: vec.c vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-$(RM) common$(SUFX).o
	-$(RM) rnd$(SUFX).o
	-$(RM) vec$(SUFX).o
	-$(RM) libvecjac$(SUFX).a
