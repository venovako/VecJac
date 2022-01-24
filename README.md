# VecJac
Vectorization of the Jacobi-type methods for the SVD and the EVD

(... work in progress ...)

## Building

### Prerequisites

A recent Intel CPU (with at least the AVX512F instruction subset) and oneAPI Base and HPC Toolkits (e.g., 2021.2) on a 64-bit Linux (e.g., CentOS 7.9) are required; macOS (e.g., Big Sur) should also be supported (not tested).

First, clone and build [JACSD](https://github.com/venovako/JACSD) repository, with the same parent directory as this one.  In fact, only the ``jstrat`` library is required to be built there.

### Make options

Run ``make`` in the ``src`` subdirectory as follows:
```bash
make [COMPILER=x64|x200] [NDEBUG=0|1|2|3|4|5] [ABI=ilp64|lp64] [FPU=precise|strict] [WP=q|l] [MKL=sequential|intel_thread] [SLEEF=/path/to/sleef] [all|clean|help]
```
where ``COMPILER`` should be set to ``x64`` for Xeons (not tested), or to ``x200`` for Xeon Phi KNLs, respectively.
Here, ``NDEBUG`` should be set to the desired optimization level (``3`` is a sensible choice).
As a hack, setting ``MKL`` *explicitly* to ``sequential`` turns off OpenMP (otherwise it is turned on, unless debugging).
Other options should not be used unless their consequences are fully understood.
For example, ``make COMPILER=x200 NDEBUG=3 clean all`` will trigger a full, release-mode rebuild for the KNLs.

## Running

The ``etc`` subdirectory contains a script that should be sourced *before* launching the executables.
The script should be modified for any Linux distribution different from RedHat/CentOS 7.
On newer Intel Xeon CPUs it should probably be sufficient to query the CPUID leaf ``0x15`` for the TSC frequency in run-time.
The OpenMP-enabled executables require setting ``OMP_NUM_THREADS`` environment variable for the upper bound on the number of threads to be used.

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
