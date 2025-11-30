# VecJac
Vectorization of the Jacobi-type methods for the SVD and the 2x2 EVD

This software is a supplementary material for the paper
doi:[10.1137/22M1478847](https://doi.org/10.1137/22M1478847 "Vectorization of a thread-parallel Jacobi singular value decomposition method")
(arXiv:[2202.08361](https://arxiv.org/abs/2202.08361 "Vectorization of a thread-parallel Jacobi singular value decomposition method") \[math.NA\]).

## Building

### Prerequisites

A recent Intel CPU (with at least the AVX512F instruction subset) and oneAPI Base and HPC Toolkits (e.g., 2021.2 or newer) on a 64-bit Linux (e.g., CentOS 7.9) are required.

First, clone and build [libpvn](https://github.com/venovako/libpvn) repository, with the same parent directory as this one has (e.g., ``venovako/libpvn`` and ``venovako/VecJac``), with the same compilers and the (no-)debug mode as it is meant to be used here.
Also, set the make variable ``OPENMP`` to at least ``0``.

### Make options

Run ``make`` in the ``src`` subdirectory as follows:
```bash
make [COMPILER=x64x|x200|gnu] [MARCH=...] [NDEBUG=0|1|2|3|...|g] [ABI=ilp64|lp64] [FPU=precise|strict] [WP=q|l] [MKL=sequential|intel_thread] [SLEEF=/path/to/sleef] [CR_MATH=/path/to/core-math] [all|clean|help]
```
where ``COMPILER`` should be set to ``x64x`` for Xeon processors or to ``x200`` for Xeon Phi KNLs, and ``NDEBUG`` to the desired optimization level (``3`` is a sensible choice).

The option ``gnu`` selects the GCC family of compilers, for which the [SLEEF](https://sleef.org)-related parts cannot yet be built, and ``SLEEF`` has to remain unset.
This option might cause the code to miscompile and should *not* be used at the moment!
Neither ``x64x`` is guaranteed to be safe, since the testing was performed only with the obsolete non-LLVM Intel compilers (see ``src/x64.mk``).

As a hack, setting ``MKL`` *explicitly* to ``sequential`` turns off OpenMP (otherwise it is turned on, unless debugging).
The ``MARCH`` option might be set to a particular CPU generation (e.g., ``ICELAKE-SERVER``, or its lowercase variant for ``gnu``) or left undefined to take the default value of ``Host`` (or ``common-avx512`` for ``x200``).
Other options should not be used unless their consequences are fully understood.
For example, ``make COMPILER=x64x NDEBUG=3 clean all`` will trigger a full, release-mode rebuild for Xeons.

## Running

The ``etc`` subdirectory contains ``env.sh``, the script that should be sourced *before* launching the executables.
The script should be modified for any Linux distribution different from RedHat/CentOS 7.
On newer Intel Xeon CPUs it might be sufficient to query the CPUID leaf ``0x15`` for the TSC frequency in run-time, but since it is not guaranteed to work, a safer timing procedure, using ``CLOCK_MONOTONIC_RAW``, is default there.
The OpenMP-enabled executables require setting ``OMP_NUM_THREADS`` environment variable for the upper bound on the number of threads to be used.

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
