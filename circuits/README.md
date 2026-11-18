# Explicit circuits for selected 2028 benchmarks

This directory holds explicit logical-level circuits for four of the report's 2028 benchmarks (Chs. 2, 4, 5 and 7). Each circuit is written out gate by gate by a generator, exported as OpenQASM 3, checked against a numpy reference that does not use the circuit code, and counted by walking its gate list. Every arbitrary rotation is kept as an `rz` with its angle, and every multi-controlled gate is decomposed to `ccx`, so a circuit can be resynthesized, recompiled or rescheduled without rebuilding it. The counts are set next to the line items of the chapter models in `estimates/`. Each circuit was then reviewed in a second pass that recounted the QASM text with its own parser, wrote its own reference, and ran single-fault mutations of the generator. The per-circuit `README.md` has the full detail. This file is the index.

Each subdirectory has the same layout: `build.py` (generator), `count.py` and `counts.json` (walked counts, model values, ratios), `test_circuit.py` (pytest suite), the QASM exports, and a `README.md`. Chs. 2, 5 and 7 also have a second checker from the review (`independent_check.py` or `verify_independent.py`).

## Conventions

- Counts are per shot of one circuit on one machine, serial.
- T total = 7 T per Toffoli (`common.toffoli_t`, textbook) + c_T per rotation + the T gates already in the circuit. c_T is `common.t_per_rotation(common.eps_rot_for(N_rot, eps_syn))` with N_rot taken from the circuit. This is the repeat-until-success fit, an expected count. No rotation was synthesized.
- eps_syn = 1e-2 per shot (`common.EPS_SYN`) for Chs. 2, 5 and 7. Ch. 4 uses 1.784e-4, the chapter's own departure (model intermediate `eps_syn`).
- "Fixed T" counts `t` and `tdg` gates in the circuit.
- Logical qubits are the allocated register with all workspace. No synthesis ancilla is included in any circuit.
- T-depth is a longest-path count under a convention that differs by directory, given in the second table. None is a lattice-surgery schedule.

## The four circuits

| directory | chapter and instance | logical qubits, circuit / model | Toffolis | rotations | fixed T | T total, circuit / model | ratio | verified, and at what size | verdict of the review |
|---|---|---|---|---|---|---|---|---|---|
| `ch02_four_nucleon` | Ch. 2, neutrino-nucleus response. Four-nucleon charge correlator on a periodic 3^3 lattice. Complete Hadamard test: 8 preparation layers, current insertion, 6 controlled second-order steps, the (p0, p0) real part. | 36 / 36 | 2,240 | 1,192 | 744 | 38,948 / 38,948 | 1.000 | The full 36-qubit file, simulated exactly from the QASM text by splitting on the ancilla (2^24 amplitudes per branch), error 3.1e-15. It reproduces the chapter's biases at 2, 4 and 6 steps. Complete circuits as statevectors at 12, 18, 20 and 24 qubits. Every block at full-instance size. 39 tests. | sound with caveats |
| `ch04_logdet_staggered` | Ch. 4, hybrid lattice QCD. log det W of free staggered fermions at V = 4^4, m0 = 0.4, K = 1. Hadamard test on a degree-84 QET of the block encoding of arXiv:2407.13080 Sec. II.2. One shot, `literal` variant. | 18 / 18 | 4,032 | 336 | 0 | 35,342 / 35,442 | 0.997 | The full 18-qubit, d = 84 circuit as a statevector, all three variants, four sites, deviation about 1e-13. log det from the noiseless circuit is 150.1774 against the exact 150.5528. Oracles and block encoding up to L = 8, D = 4 (22 qubits). 86 tests. | sound with caveats |
| `ch05_bt_link_step` | Ch. 5, QGP transport. 2+1D pure-gauge SU(2) on the binary tetrahedral group, 4x4 periodic lattice, H_KS. One Trotter step (the first result): one electric layer on 32 links, one magnetic layer on 16 plaquettes. | 176 / 177 to 178 | 1,408 | 416 | 640 | 17,994 / 94,729 | 0.190 | Both 176-qubit files (one step and six steps) block by block: 48 and 288 contiguous blocks that exhaust the gate list, each electric block as a dense unitary on 7 qubits, each plaquette block on all 24^4 configurations, errors at or below 6e-15. One open plaquette step as a 22-qubit statevector. Not run as one 176-qubit statevector. 29 tests. | sound with caveats |
| `ch07_scalar_ds_step` | Ch. 7, quantum fields in the early universe. Free scalar on de Sitter, periodic 4^3 lattice, K = 8 (3 qubits per site). One second-order step at h = 0.5/H_inf with both end shears. Evolution step only. | 192 / 192 system, 222 to 242 with N_anc | 0 | 3,584 | 1,664 | 72,661 / 64,753 | 1.122 | The 192-qubit file as an exact factorization into site factors and link ZZ rotations, each against a dense reference (7e-15, and 3e-16 on the 1,728 link angles). The full step as a unitary on 9 qubits and as a statevector at 12, 18 and 24 qubits. Not run as one 192-qubit statevector. 36 tests. | sound with caveats |

Other variants that are built and counted, for recompilation studies:

| directory | variant | Toffolis | rotations | fixed T | T total | T / model | T-depth | T-depth convention |
|---|---|---|---|---|---|---|---|---|
| `ch02_four_nucleon` | baseline, 6 steps | 2,240 | 1,192 | 744 | 38,948 | 1.000 | 15,757 | rz = c_T sequential T, Toffoli = 7-T circuit of T-depth 4. With the Toffoli as one layer: 10,078. |
| | 4 steps | 1,896 | 1,010 | 648 | 32,866 | 1.000 | | |
| | 2 steps | 1,552 | 828 | 552 | 26,812 | 1.000 | | |
| | ancilla phases merged, 6 steps | 2,240 | 1,037 | 744 | 35,900 | 0.922 | 15,635 | same |
| `ch04_logdet_staggered` | `literal` (model's upper end) | 4,032 | 336 | 0 | 35,342 | 0.997 | 7,455 | Toffoli one layer, rz = c_T layers, commuting gates overlap. In strict gate order: 10,899. |
| | `merged` (model's lower end) | 1,680 | 336 | 0 | 18,878 | 0.995 | 7,455 | same |
| | `compact` (not in the model) | 1,168 | 335 | 0 | 15,272 | 0.80 of the model's lower end | 7,407 | same |
| `ch05_bt_link_step` | one step, default | 1,408 | 416 | 640 | 17,994 | 0.190 | 936 | Toffoli one layer, T one layer, rz = c_T layers, as soon as possible on the circuit's qubits |
| | six steps (full benchmark) | 8,448 | 2,496 | 3,840 | 111,671 | 0.188 | 5,963 | same |
| | one step, U_mul run backwards in the uncompute | 1,280 | 416 | 640 | 17,098 | | 932 | same |
| | one step, second order (two magnetic half layers) | 2,560 | 544 | 640 | 28,486 | | 1,301 | same |
| | one step, 2 or 64 workspace qubits (162 or 224 qubits) | 1,408 | 416 | 640 | 17,994 | 0.190 | 5,028 or 498 | same |
| `ch07_scalar_ds_step` | `full` (both end shears) | 0 | 3,584 | 1,664 | 72,661 | 1.122 | 847 | rz = c_T, T = 1, any number of rotations in flight |
| | `in_shear_only` | 0 | 3,392 | 1,664 | 68,703 | 1.061 | 786 | same |
| | `no_outer_shear` (the model's accounting) | 0 | 3,200 | 1,664 | 64,753 | 1.000 | 741 | same |
| | `full_merged_rz` | 0 | 3,456 | 1,344 | 69,701 | 1.076 | 803 | same |
| | `no_outer_shear_merged_rz` | 0 | 3,072 | 1,344 | 61,806 | 0.954 | 700 | same |

Model values, all from running the model code:

- Ch. 2: `counts('2028', n)`, `circuit(Assumptions(), '2028', n)` and `model(a, '2028')` in `estimates/ch02_nu_nucleus.py`. `lq` = (36, 36), 164,424 shots, wall time 6,420 s. The model defines no T-depth.
- Ch. 4: `model(Assumptions(), '2028')` in `estimates/ch04_hybrid_lqcd.py`. `lq` = (18, 18), `hard_ops` = (18,978.4, 35,442.4), `projector_phases` = 85, `rotations_per_shot` = 338, `t_depth_per_shot` = 7,417.4.
- Ch. 5: `evolution_2028(a, 1)` and `evolution_2028(a, 6)` in `estimates/ch05_qgp_transport.py`, 94,729 and 593,918 T. `model(a, '2028').lq` = (177, 178). `step_depth_2028` gives 1,626 to 8,721 for one step and 10,407 to 55,714 for six.
- Ch. 7: `shot_counts(Assumptions(), '2028', 1)` in `estimates/ch07_curved_space.py`, 3,200 rotations and 64,753 T for the step. `register(a, '2028')` gives 192 system qubits and `lq` = (222, 242). `model(a, '2028').hard_ops` = 94,753 with the stated 30,000 T preparation. `Assumptions.t_depth_2028` = (346, 1308).

## Discrepancies with the model that matter

1. **Ch. 5, T count.** The circuit is 0.19 of the model at one step and at six. The six-step full benchmark is 1.12e5 T, which is 1.12 times the 1e5 reference count where the model has 5.9 times. The circuit uses constructions that differ from the published ones the model quotes: an 8-Toffoli U_mul with no ancilla (model 22), an 8-rotation U_Tr (model 11), and an irrep-sorting transform in place of the full Fourier transform (2 Toffoli, 10 T and 2 rotations against 14 Toffoli and 42 rotations). They are exact for this Hamiltonian. The irrep-sorting transform does not carry over to terms that need matrix elements inside an irrep.
2. **Ch. 5, U_mul multiplicity in `groups.PRIMCOST`.** The table value 6(d-1) per link-step is twice what the plaquette figures of arXiv:2203.02823 and arXiv:1903.08807 show (six multiplications per plaquette, 3(d-1) per link). The review confirmed this reading from the papers' sources. In the Ch. 5 model it alone moves the magnetic term from 1,124 to 662 T per link-step. It applies to every model that reads `PRIMCOST`. Needs an author check.
3. **Ch. 5, missing U_phase.** `groups.GROUPS["2T"]` has no electric phase gate. The circuit needs one per link-step, 4 Toffoli and 5 rotations, about 126 T. Here the model is low, by 4 percent of its own link-step.
4. **Ch. 5, Trotter order.** The model calls the steps second order but counts one electric and one magnetic layer per step, a first-order product. The circuit follows the model's count.
5. **Ch. 2, counts.** All line items equal the model at 2, 4 and 6 steps. This is by construction: the circuit is the study's `larger_circuit.py` written out gate by gate, and the model literals were tallied from that script. The independent evidence that the circuit is right is the simulation, not the count match.
6. **Ch. 2, qubits.** With the report's repeat-until-success synthesis (`common.RUS_ANCILLA = 1`) the register is 37, not 36. 80 rotations on the comparison flag sit inside triple contacts with all 11 workspace qubits busy and do not commute out. The alternative is ancilla-free synthesis for those 80.
7. **Ch. 2, wall time.** (Adopted in the paper 2026-11-18, card CIRC1118-4: the 2028 box quotes 4.6-7.2 h with the circuit's T-depth.) The walked T-depth is 15,757. `common.depth_exports` then gives F* = 2.47, below the 10 factories the wall-time convention assumes, and a campaign floor of 25,909 s (7.2 h) against the model's 6,420 s (1.8 h). With the Toffoli as one layer the floor is 16,570 s (4.6 h). The paper does not label the 2028 wall time a lower bound. This is for one layout, with contacts in sequence on shared workspace.
8. **Ch. 4, number of phases.** (Fixed in the model 2026-11-18, card CIRC1118-5: 84 phases, 18,878-35,342 T.) The model counts d + 1 = 85 projector phases. A QSVT of even degree d applies 84. The model is high by 8 Toffolis, 2 rotations and 100 T per shot at each end.
9. **Ch. 4, T-depth.** The model's depth assumes a measurement-based uncompute of the AND that its T count does not take. The unitary circuit gives 7,455 against 7,417. The model's statement that O_c runs under the O_A rotations holds only with commutation-aware scheduling. In strict gate order the `literal` variant has depth 10,899.
10. **Ch. 4, synthesis error.** (Stated in the paper 2026-11-18, card CIRC1118-5.) The chapter's 0.1% synthesis error holds for randomized synthesis redrawn from shot to shot (mean shift of <Z> -3.9e-5). One fixed pattern of rotation errors moves <Z> by 4.6e-3 rms, which is 2.6% of log det W. This is from an angle-error model, not from synthesized sequences.
11. **Ch. 4, qubits.** The register is 18 as the model says, but only 14 qubits are acted on, because the low bit of each coordinate is a spectator for free W. The repeat-until-success ancilla is in neither the circuit nor the model's 18. With it the register is 19 while a rotation runs.
12. **Ch. 7, end shears.** (Adopted 2026-11-18, card CIRC1118-6: the model counts the first shear, 68,703 T per step and 98,703 T per shot.) The headline circuit has 384 rotations the model does not count, the frame shears at the two ends of the shot. With both, the step is 72,661 T (+12%), and with the stated 30,000 T preparation the shot is 102,661 T, above the 1e5 target. Without them the circuit equals the model line for line. The last shear is a diagonal phase before a field-basis readout and can be dropped. The first can be absorbed only if the preparation carries that site phase, which is open because the preparation is not built. With the first shear kept the step is 68,703 T and the shot 98,703 T.
13. **Ch. 7, ancilla and depth.** (Adopted 2026-11-18, card CIRC1118-6: `t_depth_2028` = 1,364-2,253 from the `in_shear_only` width-capped schedule.) The circuit uses 192 qubits and no ancilla. The model's N_anc = 30 to 50 is an allowance, not derived from a circuit. The widest layer has 128 rotations. With at most 50 or 30 rotations in flight a greedy schedule gives a T-depth of about 1.5e3 to 2.4e3, above the model's upper input of 1,308. That input comes from `factory.json ch 7` for the retired K = 4 four-step circuit, a file that is not in the repository.

## Not yet verified

- **Rotation synthesis, all four circuits.** The `rz` gates are not compiled to Clifford+T. Every T total is the `common.py` fit applied to the walked rotation count. gridsynth is not installed, so the fit was not checked against compiled sequences. In Ch. 4, 12 of the 84 QET phases lie within eps_rot of a multiple of pi/4 and are counted as arbitrary rotations.
- **The synthesis error budget.** Not simulated in Chs. 2, 5 and 7. In Ch. 4 only through the angle-error model of `synthesis_sensitivity.py`.
- **Layout, routing, lattice-surgery schedule, factories and logical noise.** The circuits assume all-to-all connectivity. The T-depths are longest-path counts under the stated conventions.
- **Sampling.** No shots were drawn. Shot counts, estimator statistics and wall times of the models are not tested here.
- **Parts of the benchmarks that are not built.** Ch. 5: the Hadamard test, the initial state and the observable (the circuit has 176 qubits against the chapter's 178). Ch. 7: the in-vacuum preparation (30,000 T, stated in the model, a third of the shot) and the three readout settings. Thus the Ch. 7 shot count of 94,753 T is confirmed only for its 64,753 T step.
- **Whole-register simulation of the two large circuits.** Ch. 5 (176 qubits) is verified as a product of blocks, each exact on its physical domain. Ch. 7 (192 qubits) is verified as an exact factorization. Neither was run as one statevector. No periodic Ch. 5 lattice was run as a statevector (2x2 needs 40 qubits), and no Ch. 7 lattice with six distinct neighbours per site fits in 24 qubits.
- **Ch. 2.** The pytest suite stops at 24 qubits. The exact 36-qubit check is in `independent_check.py` only. Lattices with L other than 3 per axis are not supported by the generator.
- **Ch. 4.** The QET at L > 4 only at d = 2 and 8. Nothing here covers the gauged circuits of the 2033 and later rows.
- **Ch. 5.** g^2 = 1 is a stand-in in the exported angles. The electric term (Cayley-graph Laplacian over the 8 elements with Re Tr = 1) is the builder's choice. Trotter error was not examined. The variants other than the default were not checked block by block on the lattice. The published circuits of arXiv:2208.12309 and arXiv:2408.00075 were not rebuilt.
- **Ch. 7.** The accuracy of the K = 8 grid on the 4^3 lattice. The 0.14% figure is the four-site ring check, which the circuit reproduces. The model's T-depth input could not be reproduced.
- **Independent counts by another tool.** Qualtran and pyLIQTR are not installed. pyzx and stim are installed and were not used, so no automated T-count optimization was run.

## How to run the tests

Run from the repository root, one pytest session per directory. The four directories use the same module names (`build`, `count`, `test_circuit`), and pytest stops at collection when two of them are given in one session.

```
cd FNAL_Resource_Estimation
for d in ch02_four_nucleon ch04_logdet_staggered ch05_bt_link_step ch07_scalar_ds_step
do
    PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 python3 -m pytest circuits/$d/test_circuit.py -q -p no:cacheprovider
done
```

Result when this index was written, on one machine with 4 threads: 39 passed in 75 s (Ch. 2), 86 passed in 35 s (Ch. 4), 29 passed in 37 s (Ch. 5), 36 passed in 38 s (Ch. 7). The largest statevector is 24 qubits. The Ch. 7 suite peaks at 1.6 GB.

Counts and exports, and the second checkers from the review:

```
# Ch. 2 (run inside the directory)
cd circuits/ch02_four_nucleon
python3 build.py                        # QASM 3 exports
python3 count.py                        # counts.json
python3 independent_check.py --count    # recount of the QASM text, seconds
python3 independent_check.py --pairs    # exact 36-qubit check of all four current pairs, 6 to 7 min, 1.5 GB

# Chs. 4, 5 and 7 (run from the repository root)
python3 circuits/ch04_logdet_staggered/count.py --qasm        # counts.json and the three QASM files
python3 circuits/ch05_bt_link_step/build.py                   # QASM files
python3 circuits/ch05_bt_link_step/count.py                   # counts.json
python3 circuits/ch05_bt_link_step/verify_independent.py      # exported files block by block, and a recount
python3 circuits/ch07_scalar_ds_step/count.py                 # counts.json and the QASM file
python3 circuits/ch07_scalar_ds_step/independent_check.py     # recount and factorization check, about 3 s
```

Needs Python 3, numpy, scipy, qiskit and qiskit-aer. The model modules under `estimates/` are imported read-only.
