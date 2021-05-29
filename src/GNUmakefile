ifndef COMPILER
COMPILER=x64
endif # !COMPILER
include $(COMPILER).mk
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean

all: libvecjac$(SUFX).a gen2rand$(SUFX).exe tzjac2$(SUFX).exe

help:
	@echo "make [COMPILER=x64|x200] [NDEBUG=0|1|2|3|4|5] [ABI=lp64|ilp64] [FPU=precise|strict] [all|clean|help]"

libvecjac$(SUFX).a: aalloc$(SUFX).o common$(SUFX).o ddpscl$(SUFX).o djac2$(SUFX).o dnorme$(SUFX).o dnormx$(SUFX).o dscale$(SUFX).o mtxio$(SUFX).o rnd$(SUFX).o timer$(SUFX).o vec$(SUFX).o wnrme$(SUFX).o zdpscl$(SUFX).o zjac2$(SUFX).o znorme$(SUFX).o znormx$(SUFX).o zmerge$(SUFX).o zscale$(SUFX).o zsplit$(SUFX).o $(MKFS)
	$(AR) $(ARFLAGS) $@ aalloc$(SUFX).o common$(SUFX).o ddpscl$(SUFX).o djac2$(SUFX).o dnorme$(SUFX).o dnormx$(SUFX).o dscale$(SUFX).o mtxio$(SUFX).o rnd$(SUFX).o timer$(SUFX).o vec$(SUFX).o wnrme$(SUFX).o zdpscl$(SUFX).o zjac2$(SUFX).o znorme$(SUFX).o znormx$(SUFX).o zmerge$(SUFX).o zscale$(SUFX).o zsplit$(SUFX).o

aalloc$(SUFX).o: aalloc.c aalloc.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

common$(SUFX).o: common.c common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

ddpscl$(SUFX).o: ddpscl.c ddpscl.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djac2$(SUFX).o: djac2.c djac2.h dzjac2.h vecdef.h vec.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnorme$(SUFX).o: dnorme.c dnorme.h dkvsrt.h defops.h dznrme.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnormx$(SUFX).o: dnormx.c dnormx.h dznrmx.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dscale$(SUFX).o: dscale.c dscale.h dzscal.h dnormx.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

mtxio$(SUFX).o: mtxio.c mtxio.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

rnd$(SUFX).o: rnd.c rnd.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

timer$(SUFX).o: timer.c timer.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

vec$(SUFX).o: vec.c vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

wnrme$(SUFX).o: wnrme.c wnrme.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zdpscl$(SUFX).o: zdpscl.c zdpscl.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zjac2$(SUFX).o: zjac2.c zjac2.h dzjac2.h vecdef.h vec.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znorme$(SUFX).o: znorme.c znorme.h dkvsrt.h defops.h dznrme.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znormx$(SUFX).o: znormx.c znormx.h dznrmx.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zmerge$(SUFX).o: zmerge.c zmerge.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zscale$(SUFX).o: zscale.c zscale.h dzscal.h znormx.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zsplit$(SUFX).o: zsplit.c zsplit.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

gen2rand$(SUFX).exe: gen2rand.c gen2rand.h libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tzjac2$(SUFX).exe: tzjac2.c tzjac2.h libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	-$(RM) aalloc$(SUFX).o aalloc$(SUFX).optrpt
	-$(RM) common$(SUFX).o common$(SUFX).optrpt
	-$(RM) ddpscl$(SUFX).o ddpscl$(SUFX).optrpt
	-$(RM) djac2$(SUFX).o djac2$(SUFX).optrpt
	-$(RM) dnorme$(SUFX).o dnorme$(SUFX).optrpt
	-$(RM) dnormx$(SUFX).o dnormx$(SUFX).optrpt
	-$(RM) dscale$(SUFX).o dscale$(SUFX).optrpt
	-$(RM) mtxio$(SUFX).o mtxio$(SUFX).optrpt
	-$(RM) rnd$(SUFX).o rnd$(SUFX).optrpt
	-$(RM) timer$(SUFX).o timer$(SUFX).optrpt
	-$(RM) vec$(SUFX).o vec$(SUFX).optrpt
	-$(RM) zdpscl$(SUFX).o zdpscl$(SUFX).optrpt
	-$(RM) zjac2$(SUFX).o zjac2$(SUFX).optrpt
	-$(RM) znorme$(SUFX).o znorme$(SUFX).optrpt
	-$(RM) znormx$(SUFX).o znormx$(SUFX).optrpt
	-$(RM) zmerge$(SUFX).o zmerge$(SUFX).optrpt
	-$(RM) zscale$(SUFX).o zscale$(SUFX).optrpt
	-$(RM) zsplit$(SUFX).o zsplit$(SUFX).optrpt
	-$(RM) libvecjac$(SUFX).a
	-$(RM) gen2rand$(SUFX).exe gen2rand.optrpt
	-$(RM) tzjac2$(SUFX).exe tzjac2.optrpt
