#!/bin/bash
if [ -n "${TSC_FREQ_kHZ}" ]
then
    unset TSC_FREQ_kHZ
fi
# A dirty hack to get the TSC frequency (works on CentOS 7)...
if [ `if [ -r /etc/redhat-release ]; then grep -c 'release 7' /etc/redhat-release; else echo 0; fi` = 1 ]
then
    TSC_FREQ_kHZ=`dmesg | grep 'TSC clocksource calibration' | cut -d':' -f3 | cut -d' ' -f2 | sed 's/\.//g'`
else
    TSC_FREQ_kHZ=0
fi
if [ "${TSC_FREQ_kHZ}" = "0" ]
then
    export TSC_FREQ_HZ=0
else
    export TSC_FREQ_HZ=${TSC_FREQ_kHZ}000
fi
unset TSC_FREQ_kHZ
if [ -n "${KMP_DETERMINISTIC_REDUCTION}" ]
then
    unset KMP_DETERMINISTIC_REDUCTION
fi
export KMP_DETERMINISTIC_REDUCTION=TRUE
if [ -n "${OMP_PLACES}" ]
then
    unset OMP_PLACES
fi
export OMP_PLACES=CORES
if [ -n "${OMP_PROC_BIND}" ]
then
    unset OMP_PROC_BIND
fi
export OMP_PROC_BIND=SPREAD
if [ -n "${OMP_DYNAMIC}" ]
then
    unset OMP_DYNAMIC
fi
export OMP_DYNAMIC=FALSE
if [ -n "${MKL_DYNAMIC}" ]
then
    unset MKL_DYNAMIC
fi
export MKL_DYNAMIC=FALSE
