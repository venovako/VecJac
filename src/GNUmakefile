ifndef COMPILER
COMPILER=x64
endif # !COMPILER
include $(COMPILER).mk
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean

all: libvecjac$(SUFX).a bdjac2f$(SUFX).exe bdjac2s$(SUFX).exe bdlevd2$(SUFX).exe bzjac2f$(SUFX).exe bzjac2s$(SUFX).exe bzlevd2$(SUFX).exe gen2rand$(SUFX).exe tdnorme$(SUFX).exe tpjs$(SUFX).exe mkl_ver$(SUFX).exe

help:
	@echo "make [COMPILER=x64|x200] [NDEBUG=0|1|2|3|4|5] [ABI=lp64|ilp64] [FPU=precise|strict] [WP=q|l] [MKL=sequential|intel_thread] [all|clean|help]"
	@echo "Setting explicitly MKL to sequential (default) turns off OpenMP (otherwise it is turned on, unless debugging)."

libvecjac$(SUFX).a: aalloc$(SUFX).o cdpscl$(SUFX).o cjrot$(SUFX).o cnorm2$(SUFX).o cnormx$(SUFX).o cscale$(SUFX).o common$(SUFX).o ddpscl$(SUFX).o djac2f$(SUFX).o djac2s$(SUFX).o djrot$(SUFX).o dnorm2$(SUFX).o dnormb$(SUFX).o dnorme$(SUFX).o dnormf$(SUFX).o dnorms$(SUFX).o dnormx$(SUFX).o dscale$(SUFX).o dvjsvd$(SUFX).o laev2$(SUFX).o mtxio$(SUFX).o pjs$(SUFX).o rnd$(SUFX).o sdpscl$(SUFX).o sjac2f$(SUFX).o sjac2s$(SUFX).o sjrot$(SUFX).o snorm2$(SUFX).o snormx$(SUFX).o sscale$(SUFX).o timer$(SUFX).o vec$(SUFX).o vecjac$(SUFX).o wdp$(SUFX).o wnrme$(SUFX).o zdpscl$(SUFX).o zjac2f$(SUFX).o zjac2s$(SUFX).o zjrot$(SUFX).o zmerge$(SUFX).o znorme$(SUFX).o znormf$(SUFX).o znorms$(SUFX).o znormx$(SUFX).o zscale$(SUFX).o zsplit$(SUFX).o zvjsvd$(SUFX).o $(MKFS)
	$(AR) $(ARFLAGS) $@ aalloc$(SUFX).o cdpscl$(SUFX).o cjrot$(SUFX).o cnorm2$(SUFX).o cnormx$(SUFX).o cscale$(SUFX).o common$(SUFX).o ddpscl$(SUFX).o djac2f$(SUFX).o djac2s$(SUFX).o djrot$(SUFX).o dnorm2$(SUFX).o dnormb$(SUFX).o dnorme$(SUFX).o dnormf$(SUFX).o dnorms$(SUFX).o dnormx$(SUFX).o dscale$(SUFX).o dvjsvd$(SUFX).o laev2$(SUFX).o mtxio$(SUFX).o pjs$(SUFX).o rnd$(SUFX).o sdpscl$(SUFX).o sjac2f$(SUFX).o sjac2s$(SUFX).o sjrot$(SUFX).o snorm2$(SUFX).o snormx$(SUFX).o sscale$(SUFX).o timer$(SUFX).o vec$(SUFX).o vecjac$(SUFX).o wdp$(SUFX).o wnrme$(SUFX).o zdpscl$(SUFX).o zjac2f$(SUFX).o zjac2s$(SUFX).o zjrot$(SUFX).o zmerge$(SUFX).o znorme$(SUFX).o znormf$(SUFX).o znorms$(SUFX).o znormx$(SUFX).o zscale$(SUFX).o zsplit$(SUFX).o zvjsvd$(SUFX).o

aalloc$(SUFX).o: aalloc.c aalloc.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cdpscl$(SUFX).o: cdpscl.c cdpscl.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cjrot$(SUFX).o: cjrot.c cjrot.h sjrot.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cnorm2$(SUFX).o: cnorm2.c cnorm2.h scnrm2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cnormx$(SUFX).o: cnormx.c cnormx.h scnrmx.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

cscale$(SUFX).o: cscale.c cscale.h scscal.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

common$(SUFX).o: common.c common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

ddpscl$(SUFX).o: ddpscl.c ddpscl.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djac2f$(SUFX).o: djac2f.c djac2f.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djac2s$(SUFX).o: djac2s.c djac2s.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djrot$(SUFX).o: djrot.c djrot.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnorm2$(SUFX).o: dnorm2.c dnorm2.h dznrm2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnormb$(SUFX).o: dnormb.f90 $(MKFS)
	$(FC) $(FFLAGS) -c $< -o $@

dnorme$(SUFX).o: dnorme.c dnorme.h dznrme.h dkvsrt.h d8sort.h defops.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnormf$(SUFX).o: dnormf.c dnormf.h dznrmf.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnorms$(SUFX).o: dnorms.c dnorms.h dznrms.h d8sort.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnormx$(SUFX).o: dnormx.c dnormx.h dznrmx.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dscale$(SUFX).o: dscale.c dscale.h dzscal.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dvjsvd$(SUFX).o: dvjsvd.c dvjsvd.h dnormx.h dscale.h dnorm2.h dznrm2.h ddpscl.h djac2f.h djrot.h vecdef.h defops.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

laev2$(SUFX).o: laev2.c laev2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

mtxio$(SUFX).o: mtxio.c mtxio.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

pjs$(SUFX).o: pjs.c pjs.h common.h ../../JACSD/jstrat/jstrat.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

rnd$(SUFX).o: rnd.c rnd.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sdpscl$(SUFX).o: sdpscl.c sdpscl.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sjac2f$(SUFX).o: sjac2f.c sjac2f.h scjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sjac2s$(SUFX).o: sjac2s.c sjac2s.h scjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sjrot$(SUFX).o: sjrot.c sjrot.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

snorm2$(SUFX).o: snorm2.c snorm2.h scnrm2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

snormx$(SUFX).o: snormx.c snormx.h scnrmx.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

sscale$(SUFX).o: sscale.c sscale.h scscal.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

timer$(SUFX).o: timer.c timer.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

vec$(SUFX).o: vec.c vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

vecjac$(SUFX).o vecjac.mod: vecjac.F90 $(MKFS)
	$(FC) $(FFLAGS) -c $< -o $@

wdp$(SUFX).o: wdp.c wdp.h dnormb.h dnorme.h dnormf.h dnorms.h d8sort.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

wnrme$(SUFX).o: wnrme.c wnrme.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zdpscl$(SUFX).o: zdpscl.c zdpscl.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zjac2f$(SUFX).o: zjac2f.c zjac2f.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zjac2s$(SUFX).o: zjac2s.c zjac2s.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zjrot$(SUFX).o: zjrot.c zjrot.h djrot.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zmerge$(SUFX).o: zmerge.c zmerge.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znorm2$(SUFX).o: znorm2.c znorm2.h dznrm2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znorme$(SUFX).o: znorme.c znorme.h dznrme.h dkvsrt.h d8sort.h defops.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znormf$(SUFX).o: znormf.c znormf.h dznrmf.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znorms$(SUFX).o: znorms.c znorms.h dznrms.h d8sort.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znormx$(SUFX).o: znormx.c znormx.h dznrmx.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zscale$(SUFX).o: zscale.c zscale.h dzscal.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zsplit$(SUFX).o: zsplit.c zsplit.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zvjsvd$(SUFX).o: zvjsvd.c zvjsvd.h znormx.h zscale.h znorm2.h dznrm2.h zdpscl.h zjac2f.h zjrot.h vecdef.h defops.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

bdjac2f$(SUFX).exe: bdjac2f.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bdjac2s$(SUFX).exe: bdjac2s.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bdlevd2$(SUFX).exe: bdlevd2.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bzjac2f$(SUFX).exe: bzjac2f.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bzjac2s$(SUFX).exe: bzjac2s.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bzlevd2$(SUFX).exe: bzlevd2.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

gen2rand$(SUFX).exe: gen2rand.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

psort$(SUFX).o: psort.cpp psort.h $(MKFS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

tdnorme$(SUFX).o: tdnorme.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

tdnorme$(SUFX).exe: tdnorme$(SUFX).o psort$(SUFX).o libvecjac$(SUFX).a $(MKFS)
	$(CXX) $(CXXFLAGS) tdnorme$(SUFX).o psort$(SUFX).o -o $@ $(LDFLAGS) -lifcoremt

tpjs$(SUFX).exe: tpjs.c libvecjac$(SUFX).a ../../JACSD/libjstrat$(DEBUG).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

mkl_ver$(SUFX).exe: mkl_ver.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	-$(RM) aalloc$(SUFX).o aalloc$(SUFX).optrpt
	-$(RM) cdpscl$(SUFX).o cdpscl$(SUFX).optrpt
	-$(RM) cjrot$(SUFX).o cjrot$(SUFX).optrpt
	-$(RM) cnorm2$(SUFX).o cnorm2$(SUFX).optrpt
	-$(RM) cnormx$(SUFX).o cnormx$(SUFX).optrpt
	-$(RM) cscale$(SUFX).o cscale$(SUFX).optrpt
	-$(RM) common$(SUFX).o common$(SUFX).optrpt
	-$(RM) ddpscl$(SUFX).o ddpscl$(SUFX).optrpt
	-$(RM) djac2f$(SUFX).o djac2f$(SUFX).optrpt
	-$(RM) djac2s$(SUFX).o djac2s$(SUFX).optrpt
	-$(RM) djrot$(SUFX).o djrot$(SUFX).optrpt
	-$(RM) dnorm2$(SUFX).o dnorm2$(SUFX).optrpt
	-$(RM) dnormb$(SUFX).o dnormb$(SUFX).optrpt
	-$(RM) dnorme$(SUFX).o dnorme$(SUFX).optrpt
	-$(RM) dnormf$(SUFX).o dnormf$(SUFX).optrpt
	-$(RM) dnorms$(SUFX).o dnorms$(SUFX).optrpt
	-$(RM) dnormx$(SUFX).o dnormx$(SUFX).optrpt
	-$(RM) dscale$(SUFX).o dscale$(SUFX).optrpt
	-$(RM) dvjsvd$(SUFX).o dvjsvd$(SUFX).optrpt
	-$(RM) laev2$(SUFX).o laev2$(SUFX).optrpt
	-$(RM) mtxio$(SUFX).o mtxio$(SUFX).optrpt
	-$(RM) pjs$(SUFX).o pjs$(SUFX).optrpt
	-$(RM) rnd$(SUFX).o rnd$(SUFX).optrpt
	-$(RM) sdpscl$(SUFX).o sdpscl$(SUFX).optrpt
	-$(RM) sjac2f$(SUFX).o sjac2f$(SUFX).optrpt
	-$(RM) sjac2s$(SUFX).o sjac2s$(SUFX).optrpt
	-$(RM) sjrot$(SUFX).o sjrot$(SUFX).optrpt
	-$(RM) snorm2$(SUFX).o snorm2$(SUFX).optrpt
	-$(RM) snormx$(SUFX).o snormx$(SUFX).optrpt
	-$(RM) sscale$(SUFX).o sscale$(SUFX).optrpt
	-$(RM) timer$(SUFX).o timer$(SUFX).optrpt
	-$(RM) vec$(SUFX).o vec$(SUFX).optrpt
	-$(RM) vecjac$(SUFX).o vecjac$(SUFX).optrpt vecjac.mod
	-$(RM) wdp$(SUFX).o wdp$(SUFX).optrpt
	-$(RM) wnrme$(SUFX).o wnrme$(SUFX).optrpt
	-$(RM) zdpscl$(SUFX).o zdpscl$(SUFX).optrpt
	-$(RM) zjac2f$(SUFX).o zjac2f$(SUFX).optrpt
	-$(RM) zjac2s$(SUFX).o zjac2s$(SUFX).optrpt
	-$(RM) zjrot$(SUFX).o zjrot$(SUFX).optrpt
	-$(RM) znorme$(SUFX).o znorme$(SUFX).optrpt
	-$(RM) znormf$(SUFX).o znormf$(SUFX).optrpt
	-$(RM) znorms$(SUFX).o znorms$(SUFX).optrpt
	-$(RM) znormx$(SUFX).o znormx$(SUFX).optrpt
	-$(RM) zmerge$(SUFX).o zmerge$(SUFX).optrpt
	-$(RM) zscale$(SUFX).o zscale$(SUFX).optrpt
	-$(RM) zsplit$(SUFX).o zsplit$(SUFX).optrpt
	-$(RM) zvjsvd$(SUFX).o zvjsvd$(SUFX).optrpt
	-$(RM) libvecjac$(SUFX).a
	-$(RM) bdjac2f$(SUFX).exe bdjac2f.optrpt
	-$(RM) bdjac2s$(SUFX).exe bdjac2s.optrpt
	-$(RM) bdlevd2$(SUFX).exe bdlevd2.optrpt
	-$(RM) bzjac2f$(SUFX).exe bzjac2f.optrpt
	-$(RM) bzjac2s$(SUFX).exe bzjac2s.optrpt
	-$(RM) bzlevd2$(SUFX).exe bzlevd2.optrpt
	-$(RM) gen2rand$(SUFX).exe gen2rand.optrpt
	-$(RM) psort$(SUFX).o psort$(SUFX).optrpt psort$(SUFX).opt.yaml
	-$(RM) tdnorme$(SUFX).exe tdnorme$(SUFX).o tdnorme$(SUFX).optrpt
	-$(RM) tpjs$(SUFX).exe tpjs.optrpt
	-$(RM) mkl_ver$(SUFX).exe mkl_ver.optrpt
