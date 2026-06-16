#!/bin/bash
# !!! has to be run from this directory !!!
if [ -z "${GNU}" ]
then
	cd ../../libpvn/src
	make COMPILER=gcc NDEBUG=3 SAFE=DET,SV2,NRM OPENMP=0 GMP=/opt/gmp MPFR=/opt/mpfr SLEEF=/opt/sleef sleef=0 clean
	make COMPILER=gcc NDEBUG=3 SAFE=DET,SV2,NRM OPENMP=0 GMP=/opt/gmp MPFR=/opt/mpfr SLEEF=/opt/sleef sleef=0 -j all
	cd ../../VecJac/src
	make COMPILER=gnu clean all
else
	cd ../../libpvn/src
	make COMPILER=gcc COMPILER_SUFFIX=${GNU} NDEBUG=3 SAFE=DET,SV2,NRM OPENMP=0 GMP=/opt/gmp MPFR=/opt/mpfr SLEEF=/opt/sleef sleef=0 clean
	make COMPILER=gcc COMPILER_SUFFIX=${GNU} NDEBUG=3 SAFE=DET,SV2,NRM OPENMP=0 GMP=/opt/gmp MPFR=/opt/mpfr SLEEF=/opt/sleef sleef=0 -j all
	cd ../../VecJac/src
	make COMPILER=gnu clean all
fi
