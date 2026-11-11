#!/bin/sh
# Reruns the 2026-10-04 Ch. 7 step calibration (exact_over_proxy_step_2033 = 1.23). Six exact runs of ~40 min each
# (two threads each, run in parallel), the ring proxy and the 5^3 proxy, then analyze.py and ratio.py.
set -e
cd "$(dirname "$0")"
for enc in fixed chi; do for lg in "1 1" "3 0.5" "1 2"; do set -- $lg
  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 exact_step.py $1 $2 $enc 0.05,0.025,0.0125 > log_${enc}_$1_$2.txt 2>&1 &
done; done
OMP_NUM_THREADS=1 python3 proxy_red.py > log_proxy.txt 2>&1 &
OMP_NUM_THREADS=1 python3 proxy_full.py > log_proxy_full.txt 2>&1 &
wait
python3 analyze.py > analyze.txt
python3 ratio.py
