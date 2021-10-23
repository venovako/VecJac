ifndef COMPILER
COMPILER=x64
endif # !COMPILER
include $(COMPILER).mk
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean

all: libvecjac$(SUFX).a bdjac2$(SUFX).exe bdjac2f$(SUFX).exe bdlevd2$(SUFX).exe bzjac2$(SUFX).exe bzjac2f$(SUFX).exe bzlevd2$(SUFX).exe gen2rand$(SUFX).exe tdjac2$(SUFX).exe tdnorme$(SUFX).exe tpjs$(SUFX).exe tzjac2$(SUFX).exe

help:
	@echo "make [COMPILER=x64|x200] [NDEBUG=0|1|2|3|4|5] [ABI=lp64|ilp64] [FPU=precise|strict] [WP=q|l] [all|clean|help]"

libvecjac$(SUFX).a: aalloc$(SUFX).o common$(SUFX).o ddpscl$(SUFX).o djac2$(SUFX).o djac2f$(SUFX).o djrot$(SUFX).o dnorme$(SUFX).o dnormx$(SUFX).o dscale$(SUFX).o dvjsvd$(SUFX).o dzjac2$(SUFX).o laev2$(SUFX).o mtxio$(SUFX).o pjs$(SUFX).o rnd$(SUFX).o timer$(SUFX).o vec$(SUFX).o wdp$(SUFX).o wnrme$(SUFX).o zdpscl$(SUFX).o zjac2$(SUFX).o zjac2f$(SUFX).o zjrot$(SUFX).o znorme$(SUFX).o znormx$(SUFX).o zmerge$(SUFX).o zscale$(SUFX).o zsplit$(SUFX).o zvjsvd$(SUFX).o $(MKFS)
	$(AR) $(ARFLAGS) $@ aalloc$(SUFX).o common$(SUFX).o ddpscl$(SUFX).o djac2$(SUFX).o djac2f$(SUFX).o djrot$(SUFX).o dnorme$(SUFX).o dnormx$(SUFX).o dscale$(SUFX).o dvjsvd$(SUFX).o dzjac2$(SUFX).o laev2$(SUFX).o mtxio$(SUFX).o pjs$(SUFX).o rnd$(SUFX).o timer$(SUFX).o vec$(SUFX).o wdp$(SUFX).o wnrme$(SUFX).o zdpscl$(SUFX).o zjac2$(SUFX).o zjac2f$(SUFX).o zjrot$(SUFX).o znorme$(SUFX).o znormx$(SUFX).o zmerge$(SUFX).o zscale$(SUFX).o zsplit$(SUFX).o zvjsvd$(SUFX).o

aalloc$(SUFX).o: aalloc.c aalloc.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

common$(SUFX).o: common.c common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

ddpscl$(SUFX).o: ddpscl.c ddpscl.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djac2$(SUFX).o: djac2.c djac2.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djac2f$(SUFX).o: djac2f.c djac2f.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

djrot$(SUFX).o: djrot.c djrot.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnorme$(SUFX).o: dnorme.c dnorme.h dznrme.h dkvsrt.h defops.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dnormx$(SUFX).o: dnormx.c dnormx.h dznrmx.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dscale$(SUFX).o: dscale.c dscale.h dzscal.h dnormx.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dvjsvd$(SUFX).o: dvjsvd.c dvjsvd.h dscale.h dnorme.h ddpscl.h djac2.h djrot.h vecdef.h defops.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

dzjac2$(SUFX).o: dzjac2.c dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

laev2$(SUFX).o: laev2.c laev2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

mtxio$(SUFX).o: mtxio.c mtxio.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

pjs$(SUFX).o: pjs.c pjs.h common.h ../../JACSD/jstrat/jstrat.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

rnd$(SUFX).o: rnd.c rnd.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

timer$(SUFX).o: timer.c timer.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

vec$(SUFX).o: vec.c vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

wdp$(SUFX).o: wdp.c wdp.h dnorme.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

wnrme$(SUFX).o: wnrme.c wnrme.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zdpscl$(SUFX).o: zdpscl.c zdpscl.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zjac2$(SUFX).o: zjac2.c zjac2.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zjac2f$(SUFX).o: zjac2f.c zjac2f.h dzjac2.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zjrot$(SUFX).o: zjrot.c zjrot.h djrot.h vecdef.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znorme$(SUFX).o: znorme.c znorme.h dznrme.h dkvsrt.h defops.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

znormx$(SUFX).o: znormx.c znormx.h dznrmx.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zmerge$(SUFX).o: zmerge.c zmerge.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zscale$(SUFX).o: zscale.c zscale.h dzscal.h znormx.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zsplit$(SUFX).o: zsplit.c zsplit.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

zvjsvd$(SUFX).o: zvjsvd.c zvjsvd.h zscale.h znorme.h zdpscl.h zjac2.h zjrot.h vecdef.h defops.h vec.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

bdjac2$(SUFX).exe: bdjac2.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bdjac2f$(SUFX).exe: bdjac2f.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bdlevd2$(SUFX).exe: bdlevd2.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bzjac2$(SUFX).exe: bzjac2.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bzjac2f$(SUFX).exe: bzjac2f.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

bzlevd2$(SUFX).exe: bzlevd2.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

gen2rand$(SUFX).exe: gen2rand.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tdjac2$(SUFX).exe: tdjac2.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

psort$(SUFX).o: psort.cpp psort.h $(MKFS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

tdnorme$(SUFX).o: tdnorme.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

tdnorme$(SUFX).exe: tdnorme$(SUFX).o psort$(SUFX).o libvecjac$(SUFX).a $(MKFS)
	$(CXX) $(CXXFLAGS) tdnorme$(SUFX).o psort$(SUFX).o -o $@ $(LDFLAGS)

tpjs$(SUFX).exe: tpjs.c libvecjac$(SUFX).a ../../JACSD/libjstrat$(DEBUG).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

tzjac2$(SUFX).exe: tzjac2.c libvecjac$(SUFX).a $(MKFS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	-$(RM) aalloc$(SUFX).o aalloc$(SUFX).optrpt
	-$(RM) common$(SUFX).o common$(SUFX).optrpt
	-$(RM) ddpscl$(SUFX).o ddpscl$(SUFX).optrpt
	-$(RM) djac2$(SUFX).o djac2$(SUFX).optrpt
	-$(RM) djac2f$(SUFX).o djac2f$(SUFX).optrpt
	-$(RM) djrot$(SUFX).o djrot$(SUFX).optrpt
	-$(RM) dnorme$(SUFX).o dnorme$(SUFX).optrpt
	-$(RM) dnormx$(SUFX).o dnormx$(SUFX).optrpt
	-$(RM) dscale$(SUFX).o dscale$(SUFX).optrpt
	-$(RM) dvjsvd$(SUFX).o dvjsvd$(SUFX).optrpt
	-$(RM) dzjac2$(SUFX).o dzjac2$(SUFX).optrpt
	-$(RM) laev2$(SUFX).o laev2$(SUFX).optrpt
	-$(RM) mtxio$(SUFX).o mtxio$(SUFX).optrpt
	-$(RM) pjs$(SUFX).o pjs$(SUFX).optrpt
	-$(RM) rnd$(SUFX).o rnd$(SUFX).optrpt
	-$(RM) timer$(SUFX).o timer$(SUFX).optrpt
	-$(RM) vec$(SUFX).o vec$(SUFX).optrpt
	-$(RM) wdp$(SUFX).o wdp$(SUFX).optrpt
	-$(RM) wnrme$(SUFX).o wnrme$(SUFX).optrpt
	-$(RM) zdpscl$(SUFX).o zdpscl$(SUFX).optrpt
	-$(RM) zjac2$(SUFX).o zjac2$(SUFX).optrpt
	-$(RM) zjac2f$(SUFX).o zjac2f$(SUFX).optrpt
	-$(RM) zjrot$(SUFX).o zjrot$(SUFX).optrpt
	-$(RM) znorme$(SUFX).o znorme$(SUFX).optrpt
	-$(RM) znormx$(SUFX).o znormx$(SUFX).optrpt
	-$(RM) zmerge$(SUFX).o zmerge$(SUFX).optrpt
	-$(RM) zscale$(SUFX).o zscale$(SUFX).optrpt
	-$(RM) zsplit$(SUFX).o zsplit$(SUFX).optrpt
	-$(RM) zvjsvd$(SUFX).o zvjsvd$(SUFX).optrpt
	-$(RM) libvecjac$(SUFX).a
	-$(RM) bdjac2$(SUFX).exe bdjac2.optrpt
	-$(RM) bdjac2f$(SUFX).exe bdjac2f.optrpt
	-$(RM) bdlevd2$(SUFX).exe bdlevd2.optrpt
	-$(RM) bzjac2$(SUFX).exe bzjac2.optrpt
	-$(RM) bzjac2f$(SUFX).exe bzjac2f.optrpt
	-$(RM) bzlevd2$(SUFX).exe bzlevd2.optrpt
	-$(RM) gen2rand$(SUFX).exe gen2rand.optrpt
	-$(RM) tdjac2$(SUFX).exe tdjac2.optrpt
	-$(RM) psort$(SUFX).o psort$(SUFX).optrpt psort$(SUFX).opt.yaml
	-$(RM) tdnorme$(SUFX).exe tdnorme$(SUFX).o tdnorme$(SUFX).optrpt
	-$(RM) tpjs$(SUFX).exe tpjs.optrpt
	-$(RM) tzjac2$(SUFX).exe tzjac2.optrpt
