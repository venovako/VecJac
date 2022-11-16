ifndef COMPILER
COMPILER=x64
endif # !COMPILER
include $(COMPILER).mk
MKFS=GNUmakefile $(COMPILER).mk
ifneq ($(wildcard local.mk),)
include local.mk
MKFS += local.mk
endif # local.mk

.PHONY: all help clean

LIBS=libvecjac$(SUFX).a

SEXES=bslevd2$(SUFX).exe tsbjac2$(SUFX).exe tsgesvj$(SUFX).exe tsvjsvd$(SUFX).exe sflops$(SUFX).exe
DEXES=bdlevd2$(SUFX).exe tdbjac2$(SUFX).exe tdgesvj$(SUFX).exe tdvjsvd$(SUFX).exe dflops$(SUFX).exe
CEXES=bclevd2$(SUFX).exe tcbjac2$(SUFX).exe tcgesvj$(SUFX).exe tcvjsvd$(SUFX).exe
ZEXES=bzlevd2$(SUFX).exe tzbjac2$(SUFX).exe tzgesvj$(SUFX).exe tzvjsvd$(SUFX).exe
XEXES=gen2rand$(SUFX).exe mkl_ver$(SUFX).exe tpjs$(SUFX).exe
ifdef SLEEF
DEXES += tdnorme$(SUFX).exe
endif # SLEEF
EXES=$(SEXES) $(DEXES) $(CEXES) $(ZEXES) $(XEXES)

SOBJS=sbjac2$(SUFX).o sdpscl$(SUFX).o sgsscl$(SUFX).o sjrot$(SUFX).o sjrotf$(SUFX).o snorm2$(SUFX).o snormx$(SUFX).o snssvj$(SUFX).o sscale$(SUFX).o sswp$(SUFX).o svjsvd$(SUFX).o
DOBJS=dbjac2$(SUFX).o ddpscl$(SUFX).o dgsscl$(SUFX).o djrot$(SUFX).o djrotf$(SUFX).o dnorm2$(SUFX).o dnormx$(SUFX).o dnssvj$(SUFX).o dscale$(SUFX).o dswp$(SUFX).o dvjsvd$(SUFX).o
COBJS=cbjac2$(SUFX).o cdpscl$(SUFX).o cgsscl$(SUFX).o cjrot$(SUFX).o cjrotf$(SUFX).o cmerge$(SUFX).o cnorm2$(SUFX).o cnormx$(SUFX).o cnssvj$(SUFX).o cscale$(SUFX).o csplit$(SUFX).o cvjsvd$(SUFX).o
ZOBJS=zbjac2$(SUFX).o zdpscl$(SUFX).o zgsscl$(SUFX).o zjrot$(SUFX).o zjrotf$(SUFX).o zmerge$(SUFX).o znorm2$(SUFX).o znormx$(SUFX).o znssvj$(SUFX).o zscale$(SUFX).o zsplit$(SUFX).o zvjsvd$(SUFX).o
XOBJS=aalloc$(SUFX).o common$(SUFX).o laev2$(SUFX).o mtxio$(SUFX).o pjs$(SUFX).o rnd$(SUFX).o timer$(SUFX).o vec$(SUFX).o vecjac$(SUFX).o wnrme$(SUFX).o
ifdef SLEEF
DOBJS += dnormb$(SUFX).o dnorme$(SUFX).o dnormf$(SUFX).o dnorms$(SUFX).o
endif # SLEEF
OBJS=$(SOBJS) $(DOBJS) $(COBJS) $(ZOBJS) $(XOBJS)

all: $(LIBS) $(EXES)

help:
	@echo "make [COMPILER=x64|x64x|x200|gnu] [x64[x]|gnu ? CPU=...] [NDEBUG=0|1|2|3|4|5] [ABI=ilp64|lp64] [FPU=precise|strict] [WP=q|l] [MKL=sequential|intel_thread] [SLEEF=/path/to/sleef] [CR_MATH=/path/to/core-math] [all|clean|help]"
	@echo "Setting explicitly MKL to sequential (default) turns off OpenMP (otherwise it is turned on, unless debugging)."

libvecjac$(SUFX).a: $(OBJS) $(MKFS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

aalloc$(SUFX).o: aalloc.c aalloc.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cbjac2$(SUFX).o: cbjac2.c cbjac2.h scjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cdpscl$(SUFX).o: cdpscl.c cdpscl.h s2sum.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cgsscl$(SUFX).o: cgsscl.c cgsscl.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cjrot$(SUFX).o: cjrot.c cjrot.h sjrot.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cjrotf$(SUFX).o: cjrotf.c cjrotf.h sjrotf.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cmerge$(SUFX).o: cmerge.c cmerge.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cnorm2$(SUFX).o: cnorm2.c cnorm2.h scnrm2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cnormx$(SUFX).o: cnormx.c cnormx.h scnrmx.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cnssvj$(SUFX).o: cnssvj.f $(MKFS)
	$(FC) $(FFLAGS) -c $< -o $@

common$(SUFX).o: common.c common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cscale$(SUFX).o: cscale.c cscale.h scscal.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

csplit$(SUFX).o: csplit.c csplit.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cvjsvd$(SUFX).o: cvjsvd.c cvjsvd.h cnormx.h cscale.h cnorm2.h scnrm2.h cdpscl.h cgsscl.h cbjac2.h cjrot.h cjrotf.h sswp.h vecdef.h sefops.h vec.h timer.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dbjac2$(SUFX).o: dbjac2.c dbjac2.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

ddpscl$(SUFX).o: ddpscl.c ddpscl.h d2sum.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dgsscl$(SUFX).o: dgsscl.c dgsscl.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djrot$(SUFX).o: djrot.c djrot.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djrotf$(SUFX).o: djrotf.c djrotf.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnorm2$(SUFX).o: dnorm2.c dnorm2.h dznrm2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

ifdef SLEEF
dnormb$(SUFX).o: dnormb.f90 $(MKFS)
	$(FC) $(FFLAGS) -c $< -o $@

dnorme$(SUFX).o: dnorme.c dnorme.h dznrme.h dkvsrt.h defops.h d8sort.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnormf$(SUFX).o: dnormf.c dnormf.h dznrmf.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnorms$(SUFX).o: dnorms.c dnorms.h dznrms.h d8sort.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@
endif # SLEEF

dnormx$(SUFX).o: dnormx.c dnormx.h dznrmx.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnssvj$(SUFX).o: dnssvj.f $(MKFS)
	$(FC) $(FFLAGS) -c $< -o $@

dscale$(SUFX).o: dscale.c dscale.h dzscal.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dswp$(SUFX).o: dswp.c dswp.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dvjsvd$(SUFX).o: dvjsvd.c dvjsvd.h dnormx.h dscale.h dnorm2.h dznrm2.h ddpscl.h dgsscl.h dbjac2.h djrot.h djrotf.h dswp.h vecdef.h defops.h vec.h timer.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

laev2$(SUFX).o: laev2.c laev2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

mtxio$(SUFX).o: mtxio.c mtxio.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

pjs$(SUFX).o: pjs.c pjs.h common.h ../../JACSD/jstrat/jstrat.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

rnd$(SUFX).o: rnd.c rnd.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sbjac2$(SUFX).o: sbjac2.c sbjac2.h scjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sdpscl$(SUFX).o: sdpscl.c sdpscl.h s2sum.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sgsscl$(SUFX).o: sgsscl.c sgsscl.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sjrot$(SUFX).o: sjrot.c sjrot.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sjrotf$(SUFX).o: sjrotf.c sjrotf.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

snorm2$(SUFX).o: snorm2.c snorm2.h scnrm2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

snormx$(SUFX).o: snormx.c snormx.h scnrmx.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

snssvj$(SUFX).o: snssvj.f $(MKFS)
	$(FC) $(FFLAGS) -c $< -o $@

sscale$(SUFX).o: sscale.c sscale.h scscal.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sswp$(SUFX).o: sswp.c sswp.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

svjsvd$(SUFX).o: svjsvd.c svjsvd.h snormx.h sscale.h snorm2.h scnrm2.h sdpscl.h sgsscl.h sbjac2.h sjrot.h sjrotf.h sswp.h vecdef.h sefops.h vec.h timer.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

timer$(SUFX).o: timer.c timer.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

vec$(SUFX).o: vec.c vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

vecjac$(SUFX).o vecjac.mod: vecjac.F90 $(MKFS)
	$(FC) $(FFLAGS) -c $< -o $@

ifdef SLEEF
wdp$(SUFX).o: wdp.c wdp.h dnormb.h dnorme.h dnormf.h dnorms.h d8sort.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@
endif # SLEEF

wnrme$(SUFX).o: wnrme.c wnrme.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zbjac2$(SUFX).o: zbjac2.c zbjac2.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zdpscl$(SUFX).o: zdpscl.c zdpscl.h d2sum.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zgsscl$(SUFX).o: zgsscl.c zgsscl.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zjrot$(SUFX).o: zjrot.c zjrot.h djrot.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zjrotf$(SUFX).o: zjrotf.c zjrotf.h djrotf.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zmerge$(SUFX).o: zmerge.c zmerge.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znorm2$(SUFX).o: znorm2.c znorm2.h dznrm2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znormx$(SUFX).o: znormx.c znormx.h dznrmx.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znssvj$(SUFX).o: znssvj.f $(MKFS)
	$(FC) $(FFLAGS) -c $< -o $@

zscale$(SUFX).o: zscale.c zscale.h dzscal.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zsplit$(SUFX).o: zsplit.c zsplit.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zvjsvd$(SUFX).o: zvjsvd.c zvjsvd.h znormx.h zscale.h znorm2.h dznrm2.h zdpscl.h zgsscl.h zbjac2.h zjrot.h zjrotf.h dswp.h vecdef.h defops.h vec.h timer.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

bclevd2$(SUFX).exe: bclevd2.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bdlevd2$(SUFX).exe: bdlevd2.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bslevd2$(SUFX).exe: bslevd2.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bzlevd2$(SUFX).exe: bzlevd2.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

gen2rand$(SUFX).exe: gen2rand.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

mkl_ver$(SUFX).exe: mkl_ver.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tcbjac2$(SUFX).exe: tcbjac2.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tdbjac2$(SUFX).exe: tdbjac2.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tsbjac2$(SUFX).exe: tsbjac2.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tzbjac2$(SUFX).exe: tzbjac2.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

ifdef SLEEF
psort$(SUFX).o: psort.cpp psort.h $(MKFS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

tdnorme$(SUFX).o: tdnorme.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

tdnorme$(SUFX).exe: tdnorme$(SUFX).o psort$(SUFX).o wdp$(SUFX).o $(LIBS) $(MKFS)
	$(CXX) $(CXXFLAGS) tdnorme$(SUFX).o psort$(SUFX).o wdp$(SUFX).o -o $@ $(LDFLAGS) $(FLFLAGS)
endif # SLEEF

tpjs$(SUFX).exe: tpjs.c $(LIBS) ../../JACSD/libjstrat$(DEBUG).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tcgesvj$(SUFX).exe: tcgesvj.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(FLFLAGS)

tdgesvj$(SUFX).exe: tdgesvj.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(FLFLAGS)

tsgesvj$(SUFX).exe: tsgesvj.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(FLFLAGS)

tzgesvj$(SUFX).exe: tzgesvj.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(FLFLAGS)

tcvjsvd$(SUFX).exe: tcvjsvd.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tdvjsvd$(SUFX).exe: tdvjsvd.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tsvjsvd$(SUFX).exe: tsvjsvd.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tzvjsvd$(SUFX).exe: tzvjsvd.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

dflops$(SUFX).exe: dflops.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

sflops$(SUFX).exe: sflops.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	-$(RM) $(EXES)
	-$(RM) $(LIBS)
	-$(RM) $(OBJS)
	-$(RM) $(EXES:.exe=.exe.dSYM)
	-$(RM) $(EXES:$(SUFX).exe=.optrpt)
	-$(RM) $(OBJS:.o=.optrpt)
	-$(RM) vecjac.mod
ifdef SLEEF
	-$(RM) wdp$(SUFX).o wdp$(SUFX).optrpt
	-$(RM) psort$(SUFX).o psort$(SUFX).optrpt
	-$(RM) tdnorme$(SUFX).o tdnorme$(SUFX).optrpt
endif # SLEEF
