ifndef COMPILER
COMPILER=x64
endif # !COMPILER
include $(COMPILER).mk
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean

all: libvecjac$(SUFX).a gen2rand$(SUFX).exe proba$(SUFX).exe

help:
	@echo "make [COMPILER=x64|x200] [NDEBUG=0|1|2|3|4|5] [ABI=lp64|ilp64] [FPU=precise|strict] [all|clean|help]"

libvecjac$(SUFX).a: common$(SUFX).o ddpscl$(SUFX).o dnorme$(SUFX).o dnormx$(SUFX).o rnd$(SUFX).o timer$(SUFX).o vec$(SUFX).o $(MKFS)
	$(AR) $(ARFLAGS) $@ common$(SUFX).o ddpscl$(SUFX).o dnorme$(SUFX).o dnormx$(SUFX).o rnd$(SUFX).o timer$(SUFX).o vec$(SUFX).o

common$(SUFX).o: common.c common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

ddpscl$(SUFX).o: ddpscl.c ddpscl.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnorme$(SUFX).o: dnorme.c dnorme.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnormx$(SUFX).o: dnormx.c dnormx.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

rnd$(SUFX).o: rnd.c rnd.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

timer$(SUFX).o: timer.c timer.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

vec$(SUFX).o: vec.c vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

gen2rand$(SUFX).exe: gen2rand.c gen2rand.h libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

proba$(SUFX).exe: proba.c proba.h libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	-$(RM) common$(SUFX).o common$(SUFX).optrpt
	-$(RM) ddpscl$(SUFX).o ddpscl$(SUFX).optrpt
	-$(RM) dnorme$(SUFX).o dnorme$(SUFX).optrpt
	-$(RM) dnormx$(SUFX).o dnormx$(SUFX).optrpt
	-$(RM) rnd$(SUFX).o rnd$(SUFX).optrpt
	-$(RM) timer$(SUFX).o timer$(SUFX).optrpt
	-$(RM) vec$(SUFX).o vec$(SUFX).optrpt
	-$(RM) libvecjac$(SUFX).a
	-$(RM) gen2rand$(SUFX).exe gen2rand.optrpt
	-$(RM) proba$(SUFX).exe proba.optrpt
