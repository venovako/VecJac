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
make [ABI=ilp64|lp64] [WP=q|l] [MKL=sequential|intel_thread] [all|clean|help]
```

The definitions from ``src/x64x.mk`` will be used, even though they are not guaranteed to be safe, since the testing was performed only with the obsolete non-LLVM Intel compilers (see ``src/x64.mk``).

## Running

The ``etc`` subdirectory contains ``env.sh``, the script that should be sourced *before* launching the executables.
The script should be modified for any Linux distribution different from RedHat/CentOS 7.
On newer Intel Xeon CPUs it might be sufficient to query the CPUID leaf ``0x15`` for the TSC frequency in run-time, but since it is not guaranteed to work, a safer timing procedure, using ``CLOCK_MONOTONIC_RAW``, is default there.
The OpenMP-enabled executables require setting ``OMP_NUM_THREADS`` environment variable for the upper bound on the number of threads to be used.

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
