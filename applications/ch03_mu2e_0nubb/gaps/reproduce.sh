#!/bin/sh
# Regenerates every Ch. 3 input the model marks COMPUTED from these drivers (recovered 2026-11-11, see README.md).
# Needs the pinned interactions: run ../kb3g/fetch_snt.sh first. Outputs land in this folder.
# Cost on one 32-core node: minutes for everything except KB3G 48Ti in full pf (M-dim 634744), about an hour.
set -e
cd "$(dirname "$0")"
for f in kb3g gxpf1a usdb; do ln -sf ../kb3g/snt/$f.snt $f.snt; done
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
# 0+ gaps of the truncated Hamiltonians (gap_0plus_* fields; GXPF1A is the record)
for snt in kb3g gxpf1a; do for sp in f7p3 f7p3p1; do for nuc in 48Ca 48Ti; do python3 run_gaps.py $snt $sp $nuc 0 0 3; done; done; done
python3 run_gaps.py gxpf1a pf 48Ca 0 0 3
python3 run_gaps.py kb3g pf 48Ca 0 0 3
python3 run_gaps.py kb3g pf 48Ti 0 0 3
# 27Al in sd with USDB (gap_al27_usdb_mev, gap_al27_usdb_same_j_mev), then the readout quantities
python3 run_gaps.py usdb sd 27Al 5 2.5 3
python3 run_gaps.py usdb sd 27Al 1 -1 8
python3 al27_sd.py
# 4He data size D (he4_d_* fields): Minnesota NN, Lawson beta = 2 and 10
python3 he4D.py 1.0 1.0 2.0
python3 he4D.py 1.0 1.0 10.0
python3 compare.py
