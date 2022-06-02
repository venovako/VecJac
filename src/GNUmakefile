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

EXES=bclevd2$(SUFX).exe bdlevd2$(SUFX).exe bslevd2$(SUFX).exe bzlevd2$(SUFX).exe gen2rand$(SUFX).exe mkl_ver$(SUFX).exe tcbjac2$(SUFX).exe tdbjac2$(SUFX).exe tsbjac2$(SUFX).exe tzbjac2$(SUFX).exe tpjs$(SUFX).exe tdgesvj$(SUFX).exe tzgesvj$(SUFX).exe tdvjsvd$(SUFX).exe tzvjsvd$(SUFX).exe
ifdef SLEEF
EXES += tdnorme$(SUFX).exe
endif # SLEEF

OBJS=aalloc$(SUFX).o cbjac2$(SUFX).o common$(SUFX).o dbjac2$(SUFX).o ddpscl$(SUFX).o dgsscl$(SUFX).o djrot$(SUFX).o djrotf$(SUFX).o dnorm2$(SUFX).o dnormx$(SUFX).o dnssvj$(SUFX).o dscale$(SUFX).o dswp$(SUFX).o dvjsvd$(SUFX).o laev2$(SUFX).o mtxio$(SUFX).o pjs$(SUFX).o rnd$(SUFX).o sbjac2$(SUFX).o timer$(SUFX).o vec$(SUFX).o vecjac$(SUFX).o wnrme$(SUFX).o zbjac2$(SUFX).o zdpscl$(SUFX).o zjrot$(SUFX).o zjrotf$(SUFX).o zmerge$(SUFX).o znorm2$(SUFX).o znormx$(SUFX).o znssvj$(SUFX).o zscale$(SUFX).o zsplit$(SUFX).o zvjsvd$(SUFX).o
ifdef SLEEF
OBJS += dnormb$(SUFX).o dnorme$(SUFX).o dnormf$(SUFX).o dnorms$(SUFX).o
endif # SLEEF

all: $(LIBS) $(EXES)

help:
	@echo "make [COMPILER=x64|x200] [x64 ? CPU=...] [NDEBUG=0|1|2|3|4|5] [ABI=ilp64|lp64] [FPU=precise|strict] [WP=q|l] [MKL=sequential|intel_thread] [SLEEF=/path/to/sleef] [all|clean|help]"
	@echo "Setting explicitly MKL to sequential (default) turns off OpenMP (otherwise it is turned on, unless debugging)."

libvecjac$(SUFX).a: $(OBJS) $(MKFS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

aalloc$(SUFX).o: aalloc.c aalloc.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cbjac2$(SUFX).o: cbjac2.c cbjac2.h scjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

common$(SUFX).o: common.c common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dbjac2$(SUFX).o: dbjac2.c dbjac2.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

ddpscl$(SUFX).o: ddpscl.c ddpscl.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dgsscl$(SUFX).o: dgsscl.c dgsscl.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djrot$(SUFX).o: djrot.c djrot.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djrotf$(SUFX).o: djrotf.c djrotf.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnorm2$(SUFX).o: dnorm2.c dnorm2.h dznrm2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

ifdef SLEEF
dnormb$(SUFX).o: dnormb.f90 $(MKFS)
	$(FC) $(FFLAGS) -c $< -o $@

dnorme$(SUFX).o: dnorme.c dnorme.h dznrme.h dkvsrt.h d8sort.h defops.h vecdef.h vec.h common.h $(MKFS)
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

dvjsvd$(SUFX).o: dvjsvd.c dvjsvd.h dnormx.h dscale.h dnorm2.h dznrm2.h ddpscl.h dgsscl.h dbjac2.h djrotf.h djrot.h dswp.h vecdef.h defops.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

laev2$(SUFX).o: laev2.c laev2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

mtxio$(SUFX).o: mtxio.c mtxio.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

pjs$(SUFX).o: pjs.c pjs.h common.h ../../JACSD/jstrat/jstrat.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

rnd$(SUFX).o: rnd.c rnd.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sbjac2$(SUFX).o: sbjac2.c sbjac2.h scjac2.h vecdef.h vec.h common.h $(MKFS)
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

zdpscl$(SUFX).o: zdpscl.c zdpscl.h vec.h common.h $(MKFS)
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

zvjsvd$(SUFX).o: zvjsvd.c zvjsvd.h znormx.h zscale.h znorm2.h dznrm2.h zdpscl.h zbjac2.h zjrotf.h zjrot.h dswp.h vecdef.h defops.h vec.h common.h $(MKFS)
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
	$(CXX) $(CXXFLAGS) tdnorme$(SUFX).o psort$(SUFX).o wdp$(SUFX).o -o $@ $(LDFLAGS) -lifcoremt
endif # SLEEF

tpjs$(SUFX).exe: tpjs.c $(LIBS) ../../JACSD/libjstrat$(DEBUG).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tdgesvj$(SUFX).exe: tdgesvj.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) -lifcoremt

tzgesvj$(SUFX).exe: tzgesvj.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) -lifcoremt

tdvjsvd$(SUFX).exe: tdvjsvd.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tzvjsvd$(SUFX).exe: tzvjsvd.c $(LIBS) $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	-$(RM) aalloc$(SUFX).o aalloc$(SUFX).optrpt
	-$(RM) cbjac2$(SUFX).o cbjac2$(SUFX).optrpt
	-$(RM) common$(SUFX).o common$(SUFX).optrpt
	-$(RM) dbjac2$(SUFX).o dbjac2$(SUFX).optrpt
	-$(RM) ddpscl$(SUFX).o ddpscl$(SUFX).optrpt
	-$(RM) dgsscl$(SUFX).o dgsscl$(SUFX).optrpt
	-$(RM) djrot$(SUFX).o djrot$(SUFX).optrpt
	-$(RM) djrotf$(SUFX).o djrotf$(SUFX).optrpt
	-$(RM) dnorm2$(SUFX).o dnorm2$(SUFX).optrpt
	-$(RM) dnormx$(SUFX).o dnormx$(SUFX).optrpt
	-$(RM) dnssvj$(SUFX).o dnssvj$(SUFX).optrpt
	-$(RM) dscale$(SUFX).o dscale$(SUFX).optrpt
	-$(RM) dswp$(SUFX).o dswp$(SUFX).optrpt
	-$(RM) dvjsvd$(SUFX).o dvjsvd$(SUFX).optrpt
	-$(RM) laev2$(SUFX).o laev2$(SUFX).optrpt
	-$(RM) mtxio$(SUFX).o mtxio$(SUFX).optrpt
	-$(RM) pjs$(SUFX).o pjs$(SUFX).optrpt
	-$(RM) rnd$(SUFX).o rnd$(SUFX).optrpt
	-$(RM) sbjac2$(SUFX).o sbjac2$(SUFX).optrpt
	-$(RM) timer$(SUFX).o timer$(SUFX).optrpt
	-$(RM) vec$(SUFX).o vec$(SUFX).optrpt
	-$(RM) vecjac$(SUFX).o vecjac$(SUFX).optrpt vecjac.mod
	-$(RM) wnrme$(SUFX).o wnrme$(SUFX).optrpt
	-$(RM) zbjac2$(SUFX).o zbjac2$(SUFX).optrpt
	-$(RM) zdpscl$(SUFX).o zdpscl$(SUFX).optrpt
	-$(RM) zjrot$(SUFX).o zjrot$(SUFX).optrpt
	-$(RM) zjrotf$(SUFX).o zjrotf$(SUFX).optrpt
	-$(RM) zmerge$(SUFX).o zmerge$(SUFX).optrpt
	-$(RM) znorm2$(SUFX).o znorm2$(SUFX).optrpt
	-$(RM) znormx$(SUFX).o znormx$(SUFX).optrpt
	-$(RM) znssvj$(SUFX).o znssvj$(SUFX).optrpt
	-$(RM) zscale$(SUFX).o zscale$(SUFX).optrpt
	-$(RM) zsplit$(SUFX).o zsplit$(SUFX).optrpt
	-$(RM) zvjsvd$(SUFX).o zvjsvd$(SUFX).optrpt
	-$(RM) $(LIBS)
	-$(RM) bclevd2$(SUFX).exe bclevd2.optrpt
	-$(RM) bdlevd2$(SUFX).exe bdlevd2.optrpt
	-$(RM) bslevd2$(SUFX).exe bslevd2.optrpt
	-$(RM) bzlevd2$(SUFX).exe bzlevd2.optrpt
	-$(RM) gen2rand$(SUFX).exe gen2rand.optrpt
	-$(RM) mkl_ver$(SUFX).exe mkl_ver.optrpt
	-$(RM) tcbjac2$(SUFX).exe tcbjac2.optrpt
	-$(RM) tdbjac2$(SUFX).exe tdbjac2.optrpt
	-$(RM) tsbjac2$(SUFX).exe tsbjac2.optrpt
	-$(RM) tzbjac2$(SUFX).exe tzbjac2.optrpt
	-$(RM) tpjs$(SUFX).exe tpjs.optrpt
	-$(RM) tdgesvj$(SUFX).exe tdgesvj.optrpt
	-$(RM) tzgesvj$(SUFX).exe tzgesvj.optrpt
	-$(RM) tdvjsvd$(SUFX).exe tdvjsvd.optrpt
	-$(RM) tzvjsvd$(SUFX).exe tzvjsvd.optrpt
ifdef SLEEF
	-$(RM) dnormb$(SUFX).o dnormb$(SUFX).optrpt
	-$(RM) dnorme$(SUFX).o dnorme$(SUFX).optrpt
	-$(RM) dnormf$(SUFX).o dnormf$(SUFX).optrpt
	-$(RM) dnorms$(SUFX).o dnorms$(SUFX).optrpt
	-$(RM) wdp$(SUFX).o wdp$(SUFX).optrpt
	-$(RM) psort$(SUFX).o psort$(SUFX).optrpt psort$(SUFX).opt.yaml
	-$(RM) tdnorme$(SUFX).exe tdnorme$(SUFX).o tdnorme$(SUFX).optrpt
endif # SLEEF
