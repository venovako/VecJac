#!/bin/bash
# !!! has to be run from this directory !!!
if [ -z "${GNU}" ]
then
	cd ../../libpvn/src
	make COMPILER=gcc NDEBUG=3 SAFE=SV2,NRM OPENMP=0 GMP=$HOME/gnu/lib MPFR=$HOME/gnu/lib SLEEF=$HOME/gnu/lib sleef=0 clean
	make COMPILER=gcc NDEBUG=3 SAFE=SV2,NRM OPENMP=0 GMP=$HOME/gnu/lib MPFR=$HOME/gnu/lib SLEEF=$HOME/gnu/lib sleef=0 -j
	cd -
	make COMPILER=gnu clean all
else
	cd ../../libpvn/src
	make COMPILER=gcc COMPILER_SUFFIX=${GNU} NDEBUG=3 SAFE=SV2,NRM OPENMP=0 GMP=$HOME/gnu/lib MPFR=$HOME/gnu/lib SLEEF=$HOME/gnu/lib sleef=0 clean
	make COMPILER=gcc COMPILER_SUFFIX=${GNU} NDEBUG=3 SAFE=SV2,NRM OPENMP=0 GMP=$HOME/gnu/lib MPFR=$HOME/gnu/lib SLEEF=$HOME/gnu/lib sleef=0 -j
	cd -
	make COMPILER=gnu clean all
fi
