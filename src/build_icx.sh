#!/bin/bash
# !!! has to be run from this directory !!!
if [ -z "${TGT}" ]
then
	cd ../../libpvn/src
	make COMPILER=icx NDEBUG=3 SAFE=SV2,NRM OPENMP=0 GMP=/opt/gmp MPFR=/opt/mpfr SLEEF=/opt/sleef sleef=0 clean
	make COMPILER=icx NDEBUG=3 SAFE=SV2,NRM OPENMP=0 GMP=/opt/gmp MPFR=/opt/mpfr SLEEF=/opt/sleef sleef=0 -j all
	cd ../../VecJac/src
	make clean all
else
	cd ../../libpvn/src
	make COMPILER=icx MARCH=${TGT} NDEBUG=3 SAFE=SV2,NRM OPENMP=0 GMP=/opt/gmp MPFR=/opt/mpfr SLEEF=/opt/sleef sleef=0 clean
	make COMPILER=icx MARCH=${TGT} NDEBUG=3 SAFE=SV2,NRM OPENMP=0 GMP=/opt/gmp MPFR=/opt/mpfr SLEEF=/opt/sleef sleef=0 -j all
	cd ../../VecJac/src
	make clean all
fi
