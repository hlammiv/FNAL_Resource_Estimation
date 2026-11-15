# Ch. 5 quench thermalization: exact-evolution study

These are the exploratory numerics behind the chapter's quench paragraph (method (c), "Thermalization after a quench"): exact real-time evolution of small lattice gauge theories from the strong-coupling (electric) vacuum, compared with the canonical ensemble at the energy of the quench. They were run on 2026-10-05 and are not part of the resource model; no number in `ch05_qgp_transport.py` or in a benchmark box is read from them. The model's `quench_c` and `quench_ramp_a_s` inputs are the assumptions these runs were meant to support.

## Models

| script | model | lattices run |
|---|---|---|
| `z2.py` | 2+1D Z2 gauge theory on an Lx x Ly torus, zero-winding sector, through its exact Ising dual | 3x3, 4x3, 5x3 at h = 0.5 |
| `g2t.py` | 2+1D pure gauge on the binary tetrahedral group 2T (the Ch. 5 2028 group), open plaquette patch, maximal-tree gauge fixing, singlet sector | 2x1, 3x1, 2x2 at g^2 = 1.2 |
| `z2f.py` | 1+1D Z2 gauge theory with one staggered fermion, open chain, Gauss's law solved, neutral sector | N = 12, 14, 16 at m = 0.5, h_E = 1.0 |
| `g2t_ramp.py`, `g2t_ramp2.py` | 2T 2x2: a finite-speed coupling ramp from the E = 0 product state, to reach temperatures below the sudden quench's | 2x2 at g^2 = 1.2 |
| `z2trot.py` | Z2 dual on 4x4: second-order Trotterized against exact quench evolution at the thermalization step | 4x4 at h = 0.5 |
| `z2k.py` | Krylov / typicality version of `z2.py` for sectors too large for dense diagonalization (uses `qk.py`) | 4x4 started, not finished |

`qa.py` holds the shared analysis: dense diagonalization of the physical sector, canonical E(T) and local observables, quench ensembles (the electric vacuum and Fock or Gibbs-weighted product states at higher energy), their local relaxation and two-time correlators at several waiting times, the diagonal ensemble, and 1/e decay times. `qk.py` is the Krylov variant.

## Outputs

`output/` holds the JSON files the runs wrote (`out_<model>_<lattice>_<couplings>.json`) and `logs/` their console logs. Each `qa.analyze` file has the spectrum bounds, a canonical table (`canonical`, `curves`), and per quench ensemble (`quench`) the matched temperature `T_matched`, the local observables against time with their canonical values at that temperature (`local`, `local_canon_at_T`), and the correlators with their decay times (`corr`). `out_z2_4x3_h0.5.json` and `out_z2f_N14_m0.5_h1.0.json` lack the `curves` table (the canonical values at the matched temperature are present). `out_z2k_*` and `out_z2_4x4_*` were not written (the runs were stopped).

The chapter's summary figures (local observables within 1-4% of canonical, 6% on 3x3, relaxation within t of about 0.6-2.5/T, the k_min decay rate of a non-conserved operator to 5%, and the matched temperatures Ta = 0.75-1.4) were read from these outputs. The script that extracted them was not kept. The 1-4%, 6% and 5% compare the diagonal ensemble of the `vac` quench with the canonical ensemble at `T_matched` (`local_diag` against `local_canon_at_T`, `t_e_diag` against `t_e_canon`). At fixed times (`local`, `t_e_tw`) the deviations are 15-50% and up to about 40%.

## Rerunning

```
bash applications/ch05_qgp_transport/eth_quench/run.sh      # from anywhere; writes into this folder
```

`run.sh` is the original job list. It needs numpy and scipy and Python 3.12 or later. The largest jobs (2T 2x2, Z2 5x3, the 16-site chain) take hours and several GB of memory; the small lattices (Z2 3x3, 2T 2x1 and 3x1, chain N = 12) run in minutes.
