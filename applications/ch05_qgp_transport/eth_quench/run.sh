#!/usr/bin/env bash
# Rerun the Ch. 5 quench/ETH exact-evolution study. Outputs (out_*.json) and logs land in the current
# directory; the archived copies are in output/ and logs/. The largest runs take hours (2T 2x2: a 28,576-state
# dense diagonalization; Z2 5x3; the 16-site Z2 + fermion chain), so run them on a many-core machine.
# This is the job list of the original runs (2026-10-05); Python >= 3.12 (g2t_ramp2.py uses PEP 701 f-strings).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
P="${PYTHON:-python3}"
export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}"
run() { local log="$1"; shift; ( cd "$HERE" && "$P" "$@" ) > "$log" 2>&1; }
( export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
  run log_z2_33.txt z2.py 3 3 0.5 24; run log_z2_43.txt z2.py 4 3 0.5 24
  run log_2t_21.txt g2t.py 2 1 1.2;   run log_2t_31.txt g2t.py 3 1 1.2
  run log_z2_53.txt z2.py 5 3 0.5 16 ) &
( export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
  run log_f12.txt z2f.py 12 0.5 1.0 24; run log_f14.txt z2f.py 14 0.5 1.0 24; run log_f16.txt z2f.py 16 0.5 1.0 16 ) &
( export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8; run log_2t_22.txt g2t.py 2 2 1.2 ) &
( export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
  run log_trot44.txt z2trot.py 4 4 0.5 20
  run log_ramp22.txt g2t_ramp.py 2 2 1.2
  TAG=R48 RTYP=48 run log_ramp22b.txt g2t_ramp2.py 2 2 1.2 ) &
wait
